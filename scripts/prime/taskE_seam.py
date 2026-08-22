#!/usr/bin/env python3
"""TASK E item 4 — decompose the post-calling seam that drops tier-1 clip clusters.

Reads a finished `ema run` output tree and answers, for every PAS the peak
caller emitted, WHY it is or is not in `pasbed.bed`.  Nothing here changes the
caller; it re-derives the two gates the caller itself applies:

  gate A  the count matrix: a PAS whose every count came from a cell the
          `--min-read` CB filter dropped has an empty row and is dropped by
          `make_dataframe` (rows of `filterdmatrix.mtx` with no entry).
  gate B  gene assignment: `find_close` keeps TIER_1 + TIER_2 only, so a PAS
          whose nearest same-strand GENE BODY is further than the gene's
          annotated 3'UTR length (x --utr-multiplier), or whose gene has no
          annotated UTR at all, is dropped -- however much clip evidence it
          carries.

Gate B is recomputed here with the same `bedtools closest -s -D b -t first`
call and the same `assign_tier` rule the caller uses, from the caller's own
`gene_end.bed` / `utr_lengths.tsv` cache, so the classification is the
caller's, not a re-implementation of it.

Writes a fate table and (with --emit-beds) the scorer-ready point BEDs for the
rescue arms so they can be scored with scripts/prime/score_pas.sh.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile

TIER_KEEP = {"TIER_1", "TIER_2"}


def assign_tier(distance, utr_length, utr_multiplier=2.0, max_distance=5000):
    """Verbatim ema/annotate/find_close.py::assign_tier."""
    abs_dist = abs(distance)
    if utr_length > 0:
        if abs_dist <= utr_length:
            return "TIER_1"
        if abs_dist <= utr_length * utr_multiplier:
            return "TIER_2"
    if abs_dist <= max_distance:
        return "TIER_3"
    return "INTERGENIC"


def read_beds(run):
    """{pas_id: (chrom, start, end, score, strand)} over the raw caller BEDs."""
    pas = {}
    pk = os.path.join(run, "peakcalling")
    files = sorted(f for f in os.listdir(pk) if f.endswith((".pos.bed", ".neg.bed")))
    for f in files:
        with open(os.path.join(pk, f)) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                pas[int(p[3])] = (p[0], int(p[1]), int(p[2]), float(p[4]), p[5])
    return pas, files


def read_support(run):
    """{pas_id: dict} from the run-root pas_support.tsv (header-driven)."""
    path = os.path.join(run, "pas_support.tsv")
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d = dict(zip(hdr, p))
            out[int(d["pas_id"])] = d
    return out, hdr


def nonzero_rows(mtx):
    """PAS ids with >=1 entry in a MatrixMarket file (rows are pas ids)."""
    nz = set()
    with open(mtx) as fh:
        seen_header = False
        for line in fh:
            if line.startswith("%"):
                continue
            if not seen_header:
                seen_header = True
                continue
            nz.add(int(line.split(" ", 1)[0]))
    return nz


def ids_in_bed(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) > 3:
                out.add(int(p[3]))
    return out


def utr_lengths(run):
    path = os.path.join(run, "utr_lengths.tsv")
    if not os.path.exists(path):
        path = os.path.join(run, "gtf_cache", "utr_lengths.tsv")
    d = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                d[p[0]] = int(p[1])
    return d


def closest_tiers(run, pas, utr, utr_multiplier, max_distance, workdir):
    """{pas_id: (gene_id, distance, tier)} via the caller's own bedtools call."""
    gene_end = os.path.join(run, "gene_end.bed")
    if not os.path.exists(gene_end):
        gene_end = os.path.join(run, "gtf_cache", "gene_end.bed")
    q = os.path.join(workdir, "q.bed")
    with open(q, "w") as out:
        for pid, (c, s, e, sc, st) in sorted(
                pas.items(), key=lambda kv: (kv[1][0], kv[1][1], kv[1][2])):
            out.write("%s\t%d\t%d\t%d\t%g\t%s\n" % (c, s, e, pid, sc, st))
    q_sorted = os.path.join(workdir, "q.sorted.bed")
    g_sorted = os.path.join(workdir, "g.sorted.bed")
    env = dict(os.environ, LC_ALL="C")
    subprocess.run("sort -k1,1 -k2,2n %s > %s" % (q, q_sorted),
                   shell=True, check=True, env=env)
    subprocess.run("sort -k1,1 -k2,2n %s > %s" % (gene_end, g_sorted),
                   shell=True, check=True, env=env)
    res = subprocess.run(
        ["bedtools", "closest", "-s", "-D", "b", "-t", "first",
         "-a", q_sorted, "-b", g_sorted],
        capture_output=True, text=True, check=True, env=env)
    out = {}
    for line in res.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 13:
            continue
        pid = int(p[3])
        gid = p[9]
        dist = int(p[12])
        if gid == ".":
            out[pid] = (".", dist, "NOGENE")
            continue
        tier = assign_tier(dist, utr.get(gid, 0), utr_multiplier, max_distance)
        out[pid] = (gid, dist, tier)
    return out


def bins(n):
    if n <= 0:
        return "0"
    if n == 1:
        return "1"
    if n <= 4:
        return "2-4"
    if n <= 10:
        return "5-10"
    if n <= 50:
        return "11-50"
    return ">50"


BIN_ORDER = ["0", "1", "2-4", "5-10", "11-50", ">50"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", help="the run/ directory of a finished ema run")
    ap.add_argument("--utr-multiplier", type=float, default=2.0)
    ap.add_argument("--max-distance", type=int, default=5000)
    ap.add_argument("--out", default=None, help="TSV fate table")
    ap.add_argument("--emit-beds", default=None,
                    help="directory for scorer-ready point BEDs of each arm")
    ap.add_argument("--rescue-min-mol", type=int, default=10,
                    help="molecule floor for the rescue arm BEDs")
    args = ap.parse_args()

    run = args.run
    pas, bedfiles = read_beds(run)
    sup, hdr = read_support(run)
    nz = nonzero_rows(os.path.join(run, "filterdmatrix.mtx"))
    surv = ids_in_bed(os.path.join(run, "pasbed.bed"))
    # post-filter (IP) per-strand BEDs: what actually entered gene assignment
    entered = set()
    for name in ("posbed.bed", "negbed.bed"):
        p = os.path.join(run, name)
        if os.path.exists(p):
            entered |= ids_in_bed(p)
    utr = utr_lengths(run)

    with tempfile.TemporaryDirectory(prefix="seam.") as wd:
        tiers = closest_tiers(run, pas, utr, args.utr_multiplier,
                              args.max_distance, wd)

    rows = []
    fate = {}
    for pid, (c, s, e, sc, st) in pas.items():
        tier_call = int(sup[pid]["tier"]) if pid in sup else 0
        mol = int(sup[pid]["clip_umis"]) if pid in sup else 0
        gid, dist, gtier = tiers.get(pid, (".", 10 ** 9, "NOGENE"))
        if pid in surv:
            f = "kept"
        elif pid not in entered:
            f = "dropped_ip_filter"
        elif pid not in nz:
            f = "dropped_no_counts"
        elif gtier not in TIER_KEEP:
            f = "dropped_gene_tier"
        else:
            f = "dropped_other"
        fate[pid] = f
        rows.append((pid, c, s, e, sc, st, tier_call, mol, gid, dist, gtier, f))

    def table(pred, label, fh):
        sel = [r for r in rows if pred(r)]
        n = len(sel)
        counts = {}
        for r in sel:
            counts[r[11]] = counts.get(r[11], 0) + 1
        fh.write("\n## %s (n=%d)\n" % (label, n))
        for k in sorted(counts):
            fh.write("%-22s %7d  %6.2f%%\n" % (k, counts[k], 100.0 * counts[k] / n if n else 0))
        # by molecule bin x fate
        fh.write("%-8s %8s " % ("clip_mol", "n"))
        keys = sorted(counts)
        for k in keys:
            fh.write("%20s " % k)
        fh.write("\n")
        for b in BIN_ORDER:
            sub = [r for r in sel if bins(r[7]) == b]
            if not sub:
                continue
            fh.write("%-8s %8d " % (b, len(sub)))
            for k in keys:
                m = sum(1 for r in sub if r[11] == k)
                fh.write("%12d %6.2f%% " % (m, 100.0 * m / len(sub)))
            fh.write("\n")

    out = open(args.out, "w") if args.out else sys.stdout
    out.write("# seam decomposition for %s\n" % os.path.abspath(run))
    out.write("# caller BEDs: %s\n" % ", ".join(bedfiles))
    out.write("# PAS emitted by the caller: %d (tier1 %d, tier2 %d)\n" % (
        len(pas),
        sum(1 for r in rows if r[6] == 1),
        sum(1 for r in rows if r[6] == 2)))
    out.write("# entered gene assignment: %d ; nonzero matrix rows: %d ; in pasbed.bed: %d\n"
              % (len(entered), len(nz), len(surv)))
    table(lambda r: r[6] == 1, "TIER 1 (clip clusters)", out)
    table(lambda r: r[6] == 2, "TIER 2 (coverage only)", out)
    table(lambda r: r[6] == 1 and r[7] >= 2, "TIER 1 with >=2 clip molecules", out)
    # gene-tier detail for the dropped tier-1
    out.write("\n## why gene assignment dropped a tier-1 cluster\n")
    d = [r for r in rows if r[6] == 1 and r[11] == "dropped_gene_tier"]
    sub = {}
    for r in d:
        gid = r[8]
        key = (r[10], "utr=0" if utr.get(gid, 0) == 0 else "utr>0")
        sub[key] = sub.get(key, 0) + 1
    for k in sorted(sub):
        out.write("%-12s %-6s %7d\n" % (k[0], k[1], sub[k]))
    if d:
        ds = sorted(abs(r[9]) for r in d)
        out.write("distance to nearest same-strand gene body: median %d bp, "
                  "p10 %d, p90 %d, frac dist==0 %.4f\n" % (
                      ds[len(ds) // 2], ds[len(ds) // 10], ds[9 * len(ds) // 10],
                      sum(1 for x in ds if x == 0) / len(ds)))
    if args.out:
        out.close()

    # per-PAS detail, for downstream joins
    if args.out:
        det = args.out.replace(".txt", "") + ".detail.tsv"
        with open(det, "w") as fh:
            fh.write("pas_id\tchrom\tstart\tend\tscore\tstrand\ttier\tclip_umis"
                     "\tgene_id\tdistance\tgene_tier\tfate\n")
            for r in sorted(rows):
                fh.write("\t".join(str(x) for x in r) + "\n")

    if args.emit_beds:
        os.makedirs(args.emit_beds, exist_ok=True)

        def point(r):
            c, s, e, st = r[1], r[2], r[3], r[5]
            return (c, e - 1, e, r[0], r[4], st) if st == "+" else (c, s, s + 1, r[0], r[4], st)

        def write(name, pred):
            path = os.path.join(args.emit_beds, name)
            n = 0
            with open(path + ".raw", "w") as fh:
                for r in rows:
                    if pred(r):
                        fh.write("%s\t%d\t%d\t%d\t%g\t%s\n" % point(r))
                        n += 1
            env = dict(os.environ, LC_ALL="C")
            subprocess.run("sort -k1,1 -k2,2n %s.raw > %s && rm -f %s.raw"
                           % (path, path, path), shell=True, check=True, env=env)
            print("%-40s %7d" % (name, n))

        m = args.rescue_min_mol
        write("baseline_t1ge2.bed",
              lambda r: r[6] == 1 and r[7] >= 2 and r[11] == "kept")
        write("dropped_t1ge2.bed",
              lambda r: r[6] == 1 and r[7] >= 2 and r[11].startswith("dropped"))
        write("dropped_gene_t1ge2.bed",
              lambda r: r[6] == 1 and r[7] >= 2 and r[11] == "dropped_gene_tier")
        write("dropped_gene_t1ge%d.bed" % m,
              lambda r: r[6] == 1 and r[7] >= m and r[11] == "dropped_gene_tier")
        write("rescue%d_t1ge2.bed" % m,
              lambda r: r[6] == 1 and r[7] >= 2 and (
                  r[11] == "kept" or (r[11] == "dropped_gene_tier" and r[7] >= m)))


if __name__ == "__main__":
    main()
