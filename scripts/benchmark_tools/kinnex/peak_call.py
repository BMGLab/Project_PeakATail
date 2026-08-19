#!/usr/bin/env python3
"""peak_call.py <pos_counts.tsv> <mode:truth|decoy> <peaks.bed> <assign.tsv> [--win 25] [--min 5]

GREEDY 3'-terminus peak calling, replacing naive single-linkage clustering.

WHY (documented deviation from kinnex_truth_plan.md Step 4): at 104 M dedup molecules,
single-linkage joining of termini <=25 nt apart CHAINS through highly expressed loci --
MALAT1 collapsed into one 5,955 bp "PAS" and the MT rRNA locus into one 3,762 bp "PAS".
Greedy peak calling with a fixed +/-25 nt exclusion window keeps the same 25 nt merge
radius but cannot chain: repeatedly take the highest-count remaining terminus position,
emit it as a PAS, absorb every terminus within +/-25 nt, and remove them from competition.

Input : chrom pos strand n_umi ip     (one row per unique terminus position, sorted
                                       by chrom,strand,pos -- the uniq_key.tsv order)
Output: peaks.bed  chrom start end id n_umi strand rep_site n_termini_positions modal_umi width
        assign.tsv chrom pos strand peak_id|.   (SAME ORDER as input, for lockstep joining)
"""
import sys
import numpy as np

src, mode, peaks_out, assign_out = sys.argv[1:5]
WIN = 25
MINUMI = 5
for i, a in enumerate(sys.argv):
    if a == "--win":
        WIN = int(sys.argv[i + 1])
    if a == "--min":
        MINUMI = int(sys.argv[i + 1])
want_ip = 1 if mode == "decoy" else 0

chrom, pos, strand, cnt, ip = [], [], [], [], []
with open(src) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom.append(f[0]); pos.append(int(f[1])); strand.append(f[2])
        cnt.append(int(f[3])); ip.append(int(f[4]))
n = len(pos)
pos = np.asarray(pos, dtype=np.int64)
cnt = np.asarray(cnt, dtype=np.int64)
ip = np.asarray(ip, dtype=np.int8)
sys.stderr.write(f"[peak_call:{mode}] input positions={n:,} umis={cnt.sum():,}\n")

# group boundaries: input is sorted by (chrom, strand, pos)
keys = [c + "\x00" + s for c, s in zip(chrom, strand)]
bounds = [0]
for i in range(1, n):
    if keys[i] != keys[i - 1]:
        bounds.append(i)
bounds.append(n)

assign = np.full(n, -1, dtype=np.int64)
rows = []
pid = 0
n_used = 0
for gi in range(len(bounds) - 1):
    lo, hi = bounds[gi], bounds[gi + 1]
    sel = np.nonzero(ip[lo:hi] == want_ip)[0]
    if sel.size == 0:
        continue
    idx = sel + lo                      # global indices, ascending position
    p = pos[idx]; c = cnt[idx]
    n_used += idx.size
    taken = np.zeros(idx.size, dtype=bool)
    order = np.argsort(-c, kind="stable")
    for j in order:
        if taken[j]:
            continue
        a = np.searchsorted(p, p[j] - WIN, side="left")
        b = np.searchsorted(p, p[j] + WIN, side="right")
        newly = np.nonzero(~taken[a:b])[0] + a
        taken[a:b] = True
        tot = int(c[newly].sum())
        if tot < MINUMI:
            continue                    # absorbed but too weak to be a truth PAS
        pid += 1
        gid = idx[newly]
        assign[gid] = pid
        rows.append((chrom[idx[j]], int(p[newly].min()), int(p[newly].max()) + 1,
                     pid, tot, strand[idx[j]], int(p[j]), int(newly.size), int(c[j])))
sys.stderr.write(f"[peak_call:{mode}] positions_in_class={n_used:,} peaks(>={MINUMI} UMI)={len(rows):,} "
                 f"umis_in_peaks={sum(r[4] for r in rows):,}\n")

with open(peaks_out, "w") as fo:
    for r in rows:
        fo.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{mode}_{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[8]}\t{r[2]-r[1]}\n")
with open(assign_out, "w") as fo:
    for i in range(n):
        fo.write(f"{chrom[i]}\t{pos[i]}\t{strand[i]}\t"
                 f"{(mode + '_' + str(assign[i])) if assign[i] > 0 else '.'}\n")
