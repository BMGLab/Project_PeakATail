#!/usr/bin/env python3
"""ip_flag.py <genome.fa> <uniq_key.tsv> <out.tsv>
Internal-priming call for each unique 3' terminus, in input order (lockstep-safe).
Criterion (validation plan 3.1): the 18 genomic nt immediately DOWNSTREAM of the cleavage
site, read in transcript direction, contains >=12 A or >=6 consecutive A.
For '-' strand the downstream window is the forward slice [t-18, t) and the transcript-strand
A's are forward-strand T's -- so no reverse complement is needed.
"""
import sys

gpath, kpath, opath = sys.argv[1:4]
seqs, name, buf = {}, None, []
with open(gpath, "rb") as fh:
    for line in fh:
        if line[:1] == b">":
            if name is not None:
                seqs[name] = b"".join(buf).upper()
            name = line[1:].split()[0].decode()
            buf = []
        else:
            buf.append(line.rstrip())
    if name is not None:
        seqs[name] = b"".join(buf).upper()
sys.stderr.write(f"[ip_flag] loaded {len(seqs)} contigs\n")

n = nip = nshort = 0
with open(kpath) as fi, open(opath, "w") as fo:
    w = fo.write
    for line in fi:
        c, p, s = line.rstrip("\n").split("\t")
        t = int(p)
        g = seqs.get(c)
        if g is None:
            w(f"{c}\t{p}\t{s}\t0\n"); n += 1; continue
        if s == "+":
            win = g[t + 1: t + 19]
            ip = 1 if (win.count(b"A") >= 12 or b"AAAAAA" in win) else 0
        else:
            win = g[max(0, t - 18): t]
            ip = 1 if (win.count(b"T") >= 12 or b"TTTTTT" in win) else 0
        if len(win) < 18:
            nshort += 1
        w(f"{c}\t{p}\t{s}\t{ip}\n")
        n += 1; nip += ip
sys.stderr.write(f"[ip_flag] unique_positions={n} internal_priming={nip} "
                 f"({nip/max(n,1):.4f}) truncated_windows={nshort}\n")
