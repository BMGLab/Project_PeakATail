#!/usr/bin/env python3
"""TASK E item 1 — the clip-anchored cleavage-offset census (offline, read-only).

One pass over a BAM with the CALLER'S OWN primitives (`read_check`, `clip_site`,
`cluster_clip_sites`) so the clusters this script sees are the clusters the
caller called.  For every tier-1 cluster it records the member clip positions
relative to the cluster's reported cleavage (the read-weighted mode), in
TRANSCRIPT orientation (positive = downstream of the call), and reports the
candidate per-library offset statistics.

Also counts `check_clip_rate`'s own numerator/denominator over the WHOLE file,
which is the truth the strided QC estimator is validated against (item 2).

Writes a TSV of per-cluster statistics and a text summary.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam")
    ap.add_argument("--seq-len", type=int, required=True)
    ap.add_argument("--min-clip", type=int, default=6)
    ap.add_argument("--min-purity", type=float, default=0.8)
    ap.add_argument("--seed-window", type=int, default=25)
    ap.add_argument("--barcode-tag", default="CB")
    ap.add_argument("--cb-len", type=int, default=16)
    ap.add_argument("--ignore-chro", default="MT")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = whole file")
    args = ap.parse_args()

    import pysam
    from ema.countmatrix.read import read_check
    from ema.countmatrix.polya import (
        clip_site, read_umi, clip_read_ok, ClipAccumulator, cluster_clip_sites,
    )
    import ema
    sys.stderr.write("[code] ema from %s\n" % ema.__file__)

    ignore = set(x for x in args.ignore_chro.split(",") if x)
    # per (chrom, strand) accumulator
    accs: dict[tuple[str, bool], ClipAccumulator] = {}
    n_cb = n_clip = 0          # check_clip_rate's own two counters, whole file
    n_seen = 0
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam:
            n_seen += 1
            if args.max_reads and n_seen > args.max_reads:
                break
            # --- check_clip_rate's denominator, over the whole file ----------
            if not (read.is_unmapped or read.is_secondary or read.is_supplementary) \
                    and read.has_tag(args.barcode_tag):
                n_cb += 1
                if clip_site(read, args.min_clip, args.min_purity) is not None:
                    n_clip += 1
            # --- the caller's own acceptance + clip channel -------------------
            direction = bool(read.is_reverse)
            chro1, start1, end1, _s, cb = read_check(
                read=read, direction=direction, barcode=args.barcode_tag,
                barcode_len=args.cb_len, seq_len=args.seq_len,
                ignore_chro=ignore, sample_id="default",
            )
            if chro1 == 0:
                continue
            site = clip_site(read, args.min_clip, args.min_purity)
            if site is None:
                continue
            key = (read.reference_name, direction)
            acc = accs.get(key)
            if acc is None:
                acc = accs[key] = ClipAccumulator()
            acc.add(site, cb, read_umi(read), None, clip_read_ok(read))

    sys.stderr.write("[scan] records=%d cb_reads=%d clip_reads=%d rate=%.6f\n"
                     % (n_seen, n_cb, n_clip, (n_clip / n_cb) if n_cb else 0.0))

    rows = []
    pooled = Counter()          # read-weighted d histogram over all clusters
    pooled_hi = Counter()       # ... over clusters with >=10 molecules
    for (chrom, rev), acc in sorted(accs.items()):
        sites = acc.sites if hasattr(acc, "sites") else acc._sites
        for cl in cluster_clip_sites(sites, args.seed_window):
            mode = cl["mode"]
            members = cl["sites"]
            nmol = cl["numis"]
            # transcript orientation: positive d = DOWNSTREAM of the call
            ds = []
            for p in members:
                w = sites[p][0]
                d = (p - mode) if not rev else (mode - p)
                ds.append((d, w))
            tot = sum(w for _, w in ds)
            mean_d = sum(d * w for d, w in ds) / tot if tot else 0.0
            flat = sorted(d for d, w in ds for _ in range(w))
            med_d = flat[len(flat) // 2] if flat else 0
            p05 = flat[max(0, int(0.05 * len(flat)) - 1)] if flat else 0
            mind = min(d for d, _ in ds) if ds else 0
            rows.append((chrom, "-" if rev else "+", mode, len(members), nmol,
                         cl["nreads"], mean_d, med_d, p05, mind))
            for d, w in ds:
                pooled[d] += w
                if nmol >= 10:
                    pooled_hi[d] += w

    with open(args.out_prefix + ".clusters.tsv", "w") as fh:
        fh.write("chrom\tstrand\tmode\tn_positions\tn_molecules\tn_reads"
                 "\tmean_d\tmedian_d\tp05_d\tmin_d\n")
        for r in rows:
            fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\n" % r)

    def summarise(name, hist, fh):
        tot = sum(hist.values())
        if not tot:
            return
        flat = []
        for d in sorted(hist):
            flat.extend([d] * hist[d])
        mean = sum(d * c for d, c in hist.items()) / tot
        med = flat[len(flat) // 2]
        fh.write("\n[%s] clip reads=%d  mean_d=%+.3f  median_d=%+d  mode_d=%+d\n"
                 % (name, tot, mean, med, max(hist, key=lambda k: (hist[k], -k))))
        fh.write("  d :  " + " ".join("%+d" % d for d in range(-10, 11)) + "\n")
        fh.write("  %% :  " + " ".join("%.2f" % (100.0 * hist.get(d, 0) / tot)
                                       for d in range(-10, 11)) + "\n")
        fh.write("  frac |d|<=2 = %.4f ; d==0 = %.4f ; d>0 = %.4f ; d<0 = %.4f\n"
                 % (sum(c for d, c in hist.items() if abs(d) <= 2) / tot,
                    hist.get(0, 0) / tot,
                    sum(c for d, c in hist.items() if d > 0) / tot,
                    sum(c for d, c in hist.items() if d < 0) / tot))

    def bin_of(n):
        return "1" if n == 1 else "2-4" if n <= 4 else "5-10" if n <= 10 \
            else "11-50" if n <= 50 else ">50"

    with open(args.out_prefix + ".summary.txt", "w") as fh:
        fh.write("# %s\n" % os.path.abspath(args.bam))
        fh.write("# records=%d  check_clip_rate CB reads=%d  qualifying clips=%d"
                 "  WHOLE-FILE rate=%.6f (%.4f%%)\n"
                 % (n_seen, n_cb, n_clip, (n_clip / n_cb) if n_cb else 0.0,
                    100.0 * (n_clip / n_cb) if n_cb else 0.0))
        fh.write("# tier-1 clusters=%d  (>=2 mol %d, >=10 mol %d)\n"
                 % (len(rows), sum(1 for r in rows if r[4] >= 2),
                    sum(1 for r in rows if r[4] >= 10)))
        summarise("pooled, all clusters", pooled, fh)
        summarise("pooled, clusters with >=10 molecules", pooled_hi, fh)
        fh.write("\n## per-cluster offset statistics, by molecule bin "
                 "(transcript orientation, + = downstream of the call)\n")
        fh.write("%-8s %8s %10s %10s %10s %10s %10s\n" % (
            "mol", "n", "mean(mean_d)", "med(med_d)", "med(p05_d)",
            "med(min_d)", "med(npos)"))
        for b in ("1", "2-4", "5-10", "11-50", ">50"):
            sub = [r for r in rows if bin_of(r[4]) == b]
            if not sub:
                continue
            def med(vals):
                v = sorted(vals)
                return v[len(v) // 2]
            fh.write("%-8s %8d %10.3f %10d %10d %10d %10d\n" % (
                b, len(sub), sum(r[6] for r in sub) / len(sub),
                med([r[7] for r in sub]), med([r[8] for r in sub]),
                med([r[9] for r in sub]), med([r[3] for r in sub])))
        fh.write("\n## by strand (clusters with >=2 molecules)\n")
        for st in ("+", "-"):
            sub = [r for r in rows if r[1] == st and r[4] >= 2]
            if not sub:
                continue
            fh.write("%s  n=%6d  mean(mean_d)=%+.3f  median(median_d)=%+d\n" % (
                st, len(sub), sum(r[6] for r in sub) / len(sub),
                sorted(r[7] for r in sub)[len(sub) // 2]))
    print(open(args.out_prefix + ".summary.txt").read())


if __name__ == "__main__":
    main()
