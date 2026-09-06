#!/usr/bin/env python3
"""Extract per-stage, per-base read coverage at switching genes (cached step).

For each selected gene, for each mouse, for each spermatogenic stage, walk the
reads over the locus and build a per-base depth vector.  Read acceptance matches
the caller: `-F 3844` (unmapped, secondary, QC-fail, duplicate, supplementary all
excluded, the same CLIP_EXCLUDE_FLAGS the tool uses) plus a CB tag belonging to a
cell labelled with that stage.  Coverage is kept per strand so the plotting step
can verify which strand actually carries the locus rather than assuming one.

Output: results/figures/manuscript/switch_tracks_coverage.tsv.gz  (long form)
        results/figures/manuscript/switch_tracks_meta.tsv         (per-track totals)
"""
import gzip, os, re, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUT = WD / "results/figures/manuscript"
BAM = {m: WD / f"data/benchmark/gse104556/starsolo/Mouse{m[-1]}_scRNAseq/Aligned.sortedByCoord.out.bam"
       for m in ("mouse1", "mouse2")}
LAB = {m: WD / f"results/stage3_spermatogenesis_v2/{m}/input/stage3_labels.tsv" for m in BAM}
FLANK = 1500
STAGES = ["SPC", "RS", "ES"]

# 3 shortening + 1 lengthening, matching the 66/34 split of all 3,035 replicated
# opposite-direction switches.  All ES_vs_SPC, both strands represented.
GENES = [
    dict(symbol="Rragc", gene="ENSMUSG00000028646", chrom="4", strand="+",
         prox={"mouse1": 123935899, "mouse2": 123935899},
         dist={"mouse1": 123936992, "mouse2": 123936992},
         direction="shortening", pair="ES_vs_SPC", qmax=2.782e-56,
         dprox=(0.8755, 0.8738), ddist=(-0.9318, -0.9329)),
    dict(symbol="Pttg1ip", gene="ENSMUSG00000009291", chrom="10", strand="+",
         prox={"mouse1": 77597270, "mouse2": 77597270},
         dist={"mouse1": 77598731, "mouse2": 77598731},
         direction="shortening", pair="ES_vs_SPC", qmax=2.921e-61,
         dprox=(0.6648, 0.6646), ddist=(-0.8225, -0.8729)),
    dict(symbol="Pom121", gene="ENSMUSG00000053293", chrom="5", strand="-",
         prox={"mouse1": 135377832, "mouse2": 135377832},
         dist={"mouse1": 135376139, "mouse2": 135376144},
         direction="shortening", pair="ES_vs_SPC", qmax=3.105e-20,
         dprox=(0.5726, 0.5637), ddist=(-0.6218, -0.8456)),
    dict(symbol="Map3k11", gene="ENSMUSG00000004054", chrom="19", strand="+",
         prox={"mouse1": 5702861, "mouse2": 5702861},
         dist={"mouse1": 5707100, "mouse2": 5707100},
         direction="lengthening", pair="ES_vs_SPC", qmax=5.148e-33,
         dprox=(-0.9412, -0.9221), ddist=(0.8687, 0.8951)),
]

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def barcodes_by_stage(path):
    d = pd.read_csv(path, sep="\t")
    assert {"cb", "stage"} <= set(d.columns), d.columns.tolist()
    return {s: set(g.cb) for s, g in d.groupby("stage")}, d


def coverage(bam, chrom, lo, hi, cb_of_stage):
    """Per-base depth per (stage, strand) over [lo, hi), 0-based half-open."""
    width = hi - lo
    acc = {(s, t): np.zeros(width, dtype=np.int64) for s in cb_of_stage for t in "+-"}
    used = {s: 0 for s in cb_of_stage}
    lookup = {}
    for s, bcs in cb_of_stage.items():
        for b in bcs:
            lookup[b] = s
    cmd = ["samtools", "view", "-F", "3844", str(bam), f"{chrom}:{lo + 1}-{hi}"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 11)
        cb = None
        for tag in f[11].split("\t") if len(f) > 11 else ():
            if tag.startswith("CB:Z:"):
                cb = tag[5:].strip()
                break
        if cb is None:
            continue
        st = lookup.get(cb)
        if st is None:
            continue
        strand = "-" if int(f[1]) & 16 else "+"
        pos = int(f[3]) - 1
        a = acc[(st, strand)]
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "M=X":
                x0, x1 = max(pos - lo, 0), min(pos + n - lo, width)
                if x1 > x0:
                    a[x0:x1] += 1
                pos += n
            elif op in "DN":
                pos += n
        used[st] += 1
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(f"samtools failed on {chrom}:{lo}-{hi}")
    return acc, used


def main():
    rows, meta = [], []
    for m, bam in BAM.items():
        by_stage, lab = barcodes_by_stage(LAB[m])
        ncell = {s: len(v) for s, v in by_stage.items()}
        print(f"{m}: cells per stage {ncell}", flush=True)
        for g in GENES:
            sites = [g["prox"][m], g["dist"][m]]
            lo = min(sites) - FLANK
            hi = max(sites) + FLANK
            acc, used = coverage(bam, g["chrom"], lo, hi, by_stage)
            for (s, strand), vec in acc.items():
                if not vec.any():
                    continue
                nz = np.nonzero(vec)[0]
                for i in range(nz[0], nz[-1] + 1):     # contiguous span, keeps zeros inside
                    rows.append((g["symbol"], m, s, strand, lo + i, int(vec[i])))
            for s in STAGES:
                meta.append(dict(symbol=g["symbol"], gene=g["gene"], mouse=m, stage=s,
                                 n_cells=ncell.get(s, 0), reads_used=used.get(s, 0),
                                 chrom=g["chrom"], strand=g["strand"], lo=lo, hi=hi,
                                 prox=g["prox"][m], dist=g["dist"][m], direction=g["direction"],
                                 pair=g["pair"], qmax=g["qmax"],
                                 dprox=g["dprox"][0 if m == "mouse1" else 1],
                                 ddist=g["ddist"][0 if m == "mouse1" else 1]))
            print(f"  {g['symbol']:8s} {m} {g['chrom']}:{lo}-{hi}  reads/stage "
                  + " ".join(f"{s}={used.get(s,0)}" for s in STAGES), flush=True)
    cov = pd.DataFrame(rows, columns=["symbol", "mouse", "stage", "strand", "pos", "depth"])
    cov.to_csv(OUT / "switch_tracks_coverage.tsv.gz", sep="\t", index=False,
               compression="gzip")
    pd.DataFrame(meta).to_csv(OUT / "switch_tracks_meta.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT/'switch_tracks_coverage.tsv.gz'}  ({len(cov):,} rows)")
    print(f"wrote {OUT/'switch_tracks_meta.tsv'}")


if __name__ == "__main__":
    main()
