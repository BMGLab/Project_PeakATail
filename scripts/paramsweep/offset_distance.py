#!/usr/bin/env python3
"""scripts/paramsweep/offset_distance.py -- LEAD B: the signed-distance evidence.

For a PAS point BED, report the SIGNED distance from each call to its nearest
same-strand truth site, in TRANSCRIPT orientation:

    d = position(truth) - position(call)      (+ = the truth lies DOWNSTREAM)

so a coverage peak that stops short of the real cleavage site produces a
POSITIVE mode at the offset, and a correct call produces a mode at 0.  The
raw number comes from `bedtools closest -s -D b -t first`, whose sign
convention (verified on a two-line fixture: a truth 100 bp downstream of a
'+' call returns -100) is the negative of this one.

Also sweeps a post-hoc shift so the positional effect of an offset is
separable from the call-set change a real --cleavage-offset run also causes
(the shift happens before the internal-priming veto and gene assignment).
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
from collections import Counter

WD = "/mnt/ssd1/Projects/PeakATail_wd"
R = WD + "/data/references"
KIN = WD + "/results/benchmark_tools/kinnex_truth/x3p"
ENV = dict(os.environ, LC_ALL="C")


def keep(src, dst, contigs):
    ks = set(contigs.split(","))
    with open(src) as fh, open(dst + ".raw", "w") as out:
        for line in fh:
            if line.split("\t", 1)[0] in ks:
                out.write(line)
    subprocess.run("sort -k1,1 -k2,2n %s.raw > %s && rm -f %s.raw" % (dst, dst, dst),
                   shell=True, check=True, env=ENV)
    return sum(1 for _ in open(dst))


def shifted(src, dst, shift):
    with open(src) as fh, open(dst + ".raw", "w") as out:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            b = int(p[1]) + shift if p[5] == "+" else int(p[1]) - shift
            if b < 0:
                b = 0
            out.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (p[0], b, b + 1, p[3], p[4], p[5]))
    subprocess.run("sort -k1,1 -k2,2n %s.raw > %s && rm -f %s.raw" % (dst, dst, dst),
                   shell=True, check=True, env=ENV)


def signed_dists(q, t):
    res = subprocess.run(["bedtools", "closest", "-s", "-D", "b", "-t", "first",
                          "-a", q, "-b", t], capture_output=True, text=True,
                         check=True, env=ENV)
    out = []
    for line in res.stdout.splitlines():
        p = line.split("\t")
        if p[7] == "-1":            # no B feature on this contig
            continue
        out.append(-int(p[-1]))     # transcript orientation: + = truth downstream
    return out


def describe(ds, name, fh, binw=10):
    if not ds:
        fh.write("[%s] no distances\n" % name)
        return {}
    s = sorted(ds)
    n = len(s)
    hist = Counter((d // binw) * binw for d in s)
    mode_bin = max(hist, key=lambda k: hist[k])
    stats = {
        "n": n, "median": s[n // 2], "mean": sum(s) / n,
        "p10": s[n // 10], "p90": s[9 * n // 10],
        "mode_bin": mode_bin,
        "frac_within_10": sum(1 for d in s if abs(d) <= 10) / n,
        "frac_within_25": sum(1 for d in s if abs(d) <= 25) / n,
        "frac_within_100": sum(1 for d in s if abs(d) <= 100) / n,
        "frac_pos": sum(1 for d in s if d > 0) / n,
    }
    fh.write("\n[%s] n=%d  median=%+d  mean=%+.1f  p10=%+d  p90=%+d  "
             "modal %d-bp bin=[%+d,%+d)  |d|<=10 %.4f  |d|<=25 %.4f  d>0 %.4f\n"
             % (name, n, stats["median"], stats["mean"], stats["p10"], stats["p90"],
                binw, mode_bin, mode_bin + binw, stats["frac_within_10"],
                stats["frac_within_25"], stats["frac_pos"]))
    lo, hi = -200, 300
    fh.write("  bin :  " + " ".join("%+5d" % b for b in range(lo, hi, binw * 2)) + "\n")
    fh.write("  %%   :  " + " ".join("%5.2f" % (100.0 * (hist.get(b, 0) + hist.get(b + binw, 0)) / n)
                                     for b in range(lo, hi, binw * 2)) + "\n")
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed")
    ap.add_argument("--name", required=True)
    ap.add_argument("--species", default="human")
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--shifts", default="0")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--work", default="/mnt/ssd0/emaout/peakatail_benchmark/paramsweep/_scratch")
    a = ap.parse_args()

    wd = tempfile.mkdtemp(prefix="offdist.", dir=a.work)
    if a.species == "mouse":
        atl = [R + "/atlases/polyasite2.GRCm38.96.rep_sites.bed6",
               R + "/atlases/tes.protein_coding.GRCm38.102.bed6"]
    else:
        atl = [R + "/atlases/polyasite2.GRCh38.96.rep_sites.bed6",
               R + "/atlases/tes.protein_coding.GRCh38.99.bed6"]
    parts = []
    for i, p in enumerate(atl):
        keep(p, "%s/a%d.bed" % (wd, i), a.contigs)
        parts.append("%s/a%d.bed" % (wd, i))
    subprocess.run("sort -k1,1 -k2,2n %s > %s/atlas.bed" % (" ".join(parts), wd),
                   shell=True, check=True, env=ENV)
    truths = {"atlas": wd + "/atlas.bed"}
    if a.species == "human":
        keep(KIN + "/x3p_truth_t5.point.bed", wd + "/kin5.bed", a.contigs)
        truths["kinnex_t5"] = wd + "/kin5.bed"

    keep(a.bed, wd + "/q0.bed", a.contigs)
    rows = []
    with open(a.out_prefix + ".txt", "w") as fh:
        fh.write("# %s  (%s, contigs %s)\n" % (a.name, a.bed, a.contigs))
        fh.write("# d = position(truth) - position(call), transcript orientation;"
                 " + means the truth lies DOWNSTREAM of the call\n")
        for sh in [int(x) for x in a.shifts.split(",")]:
            q = wd + "/q_%d.bed" % sh
            shifted(wd + "/q0.bed", q, sh)
            for tname, tpath in truths.items():
                st = describe(signed_dists(q, tpath),
                              "%s | %s | post-hoc shift %+d" % (a.name, tname, sh), fh)
                if st:
                    st.update(name=a.name, truth=tname, shift=sh)
                    rows.append(st)
    import pandas as pd
    cols = ["name", "truth", "shift", "n", "median", "mean", "p10", "p90",
            "mode_bin", "frac_within_10", "frac_within_25", "frac_within_100", "frac_pos"]
    pd.DataFrame(rows)[cols].to_csv(a.out_prefix + ".tsv", sep="\t", index=False,
                                    float_format="%.6f")
    subprocess.run(["rm", "-rf", wd])
    print(open(a.out_prefix + ".txt").read())


if __name__ == "__main__":
    main()
