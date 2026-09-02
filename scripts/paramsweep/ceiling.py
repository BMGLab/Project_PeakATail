#!/usr/bin/env python3
"""scripts/paramsweep/ceiling.py -- LEAD A: does the candidate ceiling move?

For one run, against the slice-restricted DETECTED-GENE atlas:
  ceiling_R@100      share of detected-atlas sites with >=1 RAW CANDIDATE of any
                     tier within 100 bp (same strand) -- the most any downstream
                     threshold, filter or tier policy could ever recover
  default_R@100      share reached by the default arm (tier-1, IP-pass, >=2 mol)
  peakless_of_miss   share of the sites the DEFAULT arm misses that have no
                     candidate of any kind within 100 bp -- manuscript 25 §4's
                     class (e), recomputed per arm
  reachable_of_miss  1 - peakless_of_miss
"""
from __future__ import annotations
import argparse, os, subprocess, tempfile

WD = "/mnt/ssd1/Projects/PeakATail_wd"
ENV = dict(os.environ, LC_ALL="C")


def keep_sorted(src, dst, contigs):
    ks = set(contigs.split(","))
    with open(src) as fh, open(dst + ".raw", "w") as out:
        for line in fh:
            if line.split("\t", 1)[0] in ks:
                out.write(line)
    subprocess.run("sort -k1,1 -k2,2n %s.raw > %s && rm -f %s.raw" % (dst, dst, dst),
                   shell=True, check=True, env=ENV)


def hit_set(det, q):
    """indices of det sites with a q feature within 100 bp, same strand"""
    res = subprocess.run(["bedtools", "closest", "-s", "-d", "-t", "first",
                          "-a", det, "-b", q], capture_output=True, text=True,
                         check=True, env=ENV)
    hits = set()
    for i, line in enumerate(res.stdout.splitlines()):
        p = line.split("\t")
        d = int(p[-1])
        if 0 <= d <= 100:
            hits.add((p[0], p[1], p[5]))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", required=True); ap.add_argument("--default-arm", required=True)
    ap.add_argument("--name", required=True); ap.add_argument("--species", default="human")
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--work", default="/mnt/ssd0/emaout/peakatail_benchmark/paramsweep/_scratch")
    a = ap.parse_args()
    det = (WD + "/results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed"
           if a.species == "mouse" else
           WD + "/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed")
    wd = tempfile.mkdtemp(prefix="ceil.", dir=a.work)
    keep_sorted(det, wd + "/det.bed", a.contigs)
    keep_sorted(a.cand, wd + "/cand.bed", a.contigs)
    keep_sorted(a.default_arm, wd + "/def.bed", a.contigs)
    ndet = sum(1 for _ in open(wd + "/det.bed"))
    ncand = sum(1 for _ in open(wd + "/cand.bed"))
    hc = hit_set(wd + "/det.bed", wd + "/cand.bed")
    hd = hit_set(wd + "/det.bed", wd + "/def.bed")
    miss = ndet - len(hd)
    peakless = miss - len(hc - hd)
    print("\t".join(str(x) for x in [
        a.name, a.species, ndet, ncand, len(hc), "%.6f" % (len(hc) / ndet),
        len(hd), "%.6f" % (len(hd) / ndet), miss, peakless,
        "%.6f" % (peakless / miss) if miss else "nan",
        "%.6f" % (1 - peakless / miss) if miss else "nan"]))
    subprocess.run(["rm", "-rf", wd])


if __name__ == "__main__":
    main()
