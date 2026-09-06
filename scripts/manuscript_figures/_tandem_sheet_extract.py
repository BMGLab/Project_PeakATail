#!/usr/bin/env python3
"""Coverage + an objective switch score for every tandem 3'UTR event.

Draws nothing.  For each of the 130 genes whose two replicated sites sit inside one
3'UTR of one transcript, this pulls mouse-1 per-stage coverage over the locus and
scores two things that decide whether the call looks real in the reads:

  step_prox / step_dist  read depth 200 bp INSIDE the transcript over 200 bp
                         OUTSIDE, at each site.  A genuine 3' end drops off a
                         cliff; a call sitting mid-transcript does not.
  peak_shift             does the coverage maximum actually move between the two
                         sites across stages?  Compares which site's +/-150 bp
                         window holds more depth in the first stage against the
                         last, so a switch that is real in the reads scores True.

Output: results/figures/manuscript/tandem_sheet_{coverage.tsv.gz,meta.tsv}
"""
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUT = WD / "results/figures/manuscript"
BAM = WD / "data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam"
LAB = WD / "results/stage3_spermatogenesis_v2/mouse1/input/stage3_labels.tsv"
TAB = OUT / "tandem_utr_switches.tsv"
STAGES = ["SPC", "RS", "ES"]
FLANK, STEPW, PEAKW = 1200, 200, 150
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    T = pd.read_csv(TAB, sep="\t")
    T = T[T.utr_class == "same 3'UTR interval"].reset_index(drop=True)
    lab = pd.read_csv(LAB, sep="\t")
    by = {s: set(g.cb) for s, g in lab.groupby("stage")}
    ncell = {s: len(v) for s, v in by.items()}
    lookup = {b: s for s, bs in by.items() for b in bs}
    print(f"{len(T)} tandem 3'UTR genes; mouse1 cells per stage {ncell}", flush=True)

    rows, meta = [], []
    for i, g in T.iterrows():
        lo = min(g.prox, g.dist) - FLANK
        hi = max(g.prox, g.dist) + FLANK
        width = hi - lo
        acc = {s: np.zeros(width, np.int64) for s in STAGES}
        used = 0
        p = subprocess.Popen(["samtools", "view", "-F", "3844", str(BAM),
                              f"{g.chrom}:{lo + 1}-{hi}"],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t", 11)
            if len(f) < 12:
                continue
            cb = None
            for t in f[11].split("\t"):
                if t.startswith("CB:Z:"):
                    cb = t[5:].strip()
                    break
            st = lookup.get(cb) if cb else None
            if st is None:
                continue
            if ("-" if int(f[1]) & 16 else "+") != g.strand:
                continue
            rp = int(f[3]) - 1
            a = acc[st]
            for n, op in CIG.findall(f[5]):
                n = int(n)
                if op in "M=X":
                    x0, x1 = max(rp - lo, 0), min(rp + n - lo, width)
                    if x1 > x0:
                        a[x0:x1] += 1
                    rp += n
                elif op in "DN":
                    rp += n
            used += 1
        p.stdout.close()
        p.wait()

        pooled = sum(acc.values())
        def step(pos):
            j = pos - lo
            ins, out = ((pooled[max(j - STEPW, 0):j], pooled[j:j + STEPW]) if g.strand == "+"
                        else (pooled[j:j + STEPW], pooled[max(j - STEPW, 0):j]))
            a_, b_ = (ins.mean() if len(ins) else 0), (out.mean() if len(out) else 0)
            return (a_ + .01) / (b_ + .01), a_

        sp, dep_p = step(int(g.prox))
        sd, dep_d = step(int(g.dist))

        def frac_at_prox(st):
            v = acc[st]
            def win(pos):
                j = pos - lo
                return v[max(j - PEAKW, 0):j + PEAKW].sum()
            a_, b_ = win(int(g.prox)), win(int(g.dist))
            return a_ / (a_ + b_) if (a_ + b_) else np.nan

        f_first, f_last = frac_at_prox(STAGES[0]), frac_at_prox(STAGES[-1])
        shift = (f_last - f_first) if np.isfinite(f_last) and np.isfinite(f_first) else np.nan
        expect_up = (g.direction == "shortening")
        peak_shift = bool(np.isfinite(shift) and ((shift > 0) == expect_up) and abs(shift) >= 0.15)

        for s_, vec in acc.items():
            if vec.any():
                nz = np.nonzero(vec)[0]
                for k in range(nz[0], nz[-1] + 1):
                    if vec[k]:
                        rows.append((g["name"], s_, lo + k, int(vec[k])))
        meta.append(dict(name=g["name"], gene=g.gene, chrom=g.chrom, strand=g.strand,
                         lo=lo, hi=hi, prox=int(g.prox), dist=int(g.dist), span=int(g.span),
                         direction=g.direction, max_abs_dprop=g.max_abs_dprop, qmax=g.qmax,
                         tx=g.tx, reads=used, step_prox=round(sp, 2), step_dist=round(sd, 2),
                         depth_prox=round(dep_p, 2), depth_dist=round(dep_d, 2),
                         frac_prox_first=f_first, frac_prox_last=f_last,
                         peak_shift=peak_shift,
                         **{f"n_cells_{s_}": ncell.get(s_, 0) for s_ in STAGES}))
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(T)}", flush=True)

    pd.DataFrame(rows, columns=["name", "stage", "pos", "depth"]).to_csv(
        OUT / "tandem_sheet_coverage.tsv.gz", sep="\t", index=False, compression="gzip")
    M = pd.DataFrame(meta)
    M.to_csv(OUT / "tandem_sheet_meta.tsv", sep="\t", index=False)
    both = (M.step_prox >= 3) & (M.step_dist >= 3)
    print(f"\nboth sites at a clean 3' end (step >= 3) : {both.sum()}/{len(M)} ({both.mean():.0%})")
    print(f"coverage maximum actually shifts          : {M.peak_shift.sum()}/{len(M)} ({M.peak_shift.mean():.0%})")
    print(f"both of the above                         : {(both & M.peak_shift).sum()}/{len(M)}")
    print("wrote tandem_sheet_{coverage.tsv.gz,meta.tsv}")


if __name__ == "__main__":
    main()
