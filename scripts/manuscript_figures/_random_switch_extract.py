#!/usr/bin/env python3
"""Coverage at a RANDOM, UNFILTERED sample of replicated spermatogenesis switches.

The curated figure (_switch_track_extract.py) selects for legibility: exactly two
sites, annotation-clean, large effect.  This one deliberately does not.  It draws a
uniform random sample of genes from the replicated set with a fixed seed and keeps
whatever comes out -- three or more sites, sites hundreds of kilobases apart, small
effects, genes overlapping other genes.  It exists to answer "do PeakATail's switch
calls look real in general", which a curated panel cannot answer.

Output: results/figures/manuscript/random_switch_{coverage.tsv.gz,meta.tsv,sites.tsv}
"""
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUT = WD / "results/figures/manuscript"
REPL = WD / "results/stage3_spermatogenesis_v2/summary/switch_replicated_pas_true.tsv"
GTF = WD / "data/references/mouse/Mus_musculus.GRCm38.102.gtf"
MOUSE = "mouse1"
BAM = WD / "data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam"
LAB = WD / f"results/stage3_spermatogenesis_v2/{MOUSE}/input/stage3_labels.tsv"
STAGES = ["SPC", "RS", "ES"]
SEED, N_GENES, FLANK = 20260906, 12, 1500
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    d = pd.read_csv(REPL, sep="\t")
    rng = np.random.default_rng(SEED)
    genes = np.sort(d.gene_id.unique())
    pick = rng.choice(genes, size=N_GENES, replace=False)
    print(f"replicated set: {len(d):,} rows over {len(genes):,} genes; "
          f"uniform sample of {N_GENES} with seed {SEED}")

    sym = {}
    with open(GTF) as fh:
        for line in fh:
            if line[0] == "#" or "\tgene\t" not in line:
                continue
            m = re.search(r'gene_id "([^"]+)"', line)
            if m and m.group(1) in set(pick):
                n = re.search(r'gene_name "([^"]+)"', line)
                sym[m.group(1)] = n.group(1) if n else m.group(1)
                if len(sym) == len(pick):
                    break

    lab = pd.read_csv(LAB, sep="\t")
    by_stage = {s: set(g.cb) for s, g in lab.groupby("stage")}
    ncell = {s: len(v) for s, v in by_stage.items()}
    lookup = {b: s for s, bs in by_stage.items() for b in bs}
    print(f"{MOUSE} cells per stage: {ncell}")

    rows, meta, sites = [], [], []
    for gid in pick:
        s = d[d.gene_id == gid]
        chrom, strand = str(s.chrom.iloc[0]), s.strand.iloc[0]
        pos = sorted(s.start_m1.unique())
        lo, hi = min(pos) - FLANK, max(pos) + FLANK
        width = hi - lo
        acc = {st: np.zeros(width, np.int64) for st in STAGES}
        used = 0
        p = subprocess.Popen(["samtools", "view", "-F", "3844", str(BAM),
                              f"{chrom}:{lo + 1}-{hi}"],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t", 11)
            if len(f) < 12:
                continue
            cb = None
            for tag in f[11].split("\t"):
                if tag.startswith("CB:Z:"):
                    cb = tag[5:].strip()
                    break
            st = lookup.get(cb) if cb else None
            if st is None:
                continue
            if ("-" if int(f[1]) & 16 else "+") != strand:
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

        name = sym.get(gid, gid)
        for st, vec in acc.items():
            if vec.any():
                nz = np.nonzero(vec)[0]
                for i in range(nz[0], nz[-1] + 1):
                    if vec[i]:
                        rows.append((name, st, lo + i, int(vec[i])))
        for pp in pos:
            r = s[s.start_m1 == pp]
            sites.append(dict(symbol=name, pos=int(pp),
                              max_abs_dprop=float(r.dprop_m1.abs().max()),
                              sign=int(np.sign(r.dprop_m1.iloc[0])),
                              min_q=float(r.q_m1.min()), n_pairs=int(r.pair.nunique())))
        meta.append(dict(symbol=name, gene=gid, chrom=chrom, strand=strand,
                         lo=lo, hi=hi, span=int(max(pos) - min(pos)), n_sites=len(pos),
                         n_pairs=int(s.pair.nunique()), reads=used,
                         max_abs_dprop=float(s.dprop_m1.abs().max()),
                         min_q=float(s.q_m1.min()),
                         **{f"n_cells_{st}": ncell.get(st, 0) for st in STAGES}))
        print(f"  {name:16s} {chrom}:{lo}-{hi}  sites={len(pos)}  span={max(pos)-min(pos):>7,}  "
              f"reads={used:,}", flush=True)

    pd.DataFrame(rows, columns=["symbol", "stage", "pos", "depth"]).to_csv(
        OUT / "random_switch_coverage.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(meta).to_csv(OUT / "random_switch_meta.tsv", sep="\t", index=False)
    pd.DataFrame(sites).to_csv(OUT / "random_switch_sites.tsv", sep="\t", index=False)
    print("wrote random_switch_{coverage.tsv.gz,meta.tsv,sites.tsv}")


if __name__ == "__main__":
    main()
