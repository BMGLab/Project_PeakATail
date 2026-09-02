#!/usr/bin/env python3
"""scripts/paramsweep/resolution.py -- LEAD C: resolution and per-PAS quantification.

Recall as the scorer measures it credits BOTH members of a merged doublet
(nearest-neighbour matching is many-to-one), so an arm that splits merged calls
can gain resolution without gaining recall.  This measures the thing recall
cannot see:

  n_calls                     calls in the arm
  matched@25                  calls with an atlas site within 25 bp (same strand)
  distinct_atlas@25           DISTINCT atlas sites those calls matched
  calls_per_atlas_site        matched@25 / distinct_atlas@25  (1.0 = one call per
                              site; >1 = several calls collapsing onto one site)
  recall_sites@100            distinct DETECTED-GENE atlas sites with >=1 call
                              within 100 bp  (the scorer's recall numerator)
  many_to_one@100             share of the arm's matched calls that share their
                              nearest detected-atlas site with another call
  median_gap                  median distance to the nearest same-strand call in
                              the same arm -- the direct resolution read-out
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
from collections import Counter

WD = "/mnt/ssd1/Projects/PeakATail_wd"
R = WD + "/data/references"
ENV = dict(os.environ, LC_ALL="C")


def keep_sorted(src, dst, contigs):
    ks = set(contigs.split(","))
    with open(src) as fh, open(dst + ".raw", "w") as out:
        for line in fh:
            if line.split("\t", 1)[0] in ks:
                out.write(line)
    subprocess.run("sort -k1,1 -k2,2n %s.raw > %s && rm -f %s.raw" % (dst, dst, dst),
                   shell=True, check=True, env=ENV)


def closest(q, b, extra=()):
    res = subprocess.run(["bedtools", "closest", "-s", "-d", "-t", "first",
                          *extra, "-a", q, "-b", b],
                         capture_output=True, text=True, check=True, env=ENV)
    return res.stdout.splitlines()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed"); ap.add_argument("--name", required=True)
    ap.add_argument("--species", default="human"); ap.add_argument("--contigs", required=True)
    ap.add_argument("--work", default="/mnt/ssd0/emaout/peakatail_benchmark/paramsweep/_scratch")
    a = ap.parse_args()
    wd = tempfile.mkdtemp(prefix="res.", dir=a.work)
    if a.species == "mouse":
        atl = [R + "/atlases/polyasite2.GRCm38.96.rep_sites.bed6",
               R + "/atlases/tes.protein_coding.GRCm38.102.bed6"]
        det = WD + "/results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed"
    else:
        atl = [R + "/atlases/polyasite2.GRCh38.96.rep_sites.bed6",
               R + "/atlases/tes.protein_coding.GRCh38.99.bed6"]
        det = WD + "/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed"
    parts = []
    for i, p in enumerate(atl):
        keep_sorted(p, "%s/a%d.bed" % (wd, i), a.contigs); parts.append("%s/a%d.bed" % (wd, i))
    subprocess.run("sort -k1,1 -k2,2n %s > %s/atlas.bed" % (" ".join(parts), wd),
                   shell=True, check=True, env=ENV)
    keep_sorted(det, wd + "/det.bed", a.contigs)
    keep_sorted(a.bed, wd + "/q.bed", a.contigs)

    n = sum(1 for _ in open(wd + "/q.bed"))
    matched25 = 0
    sites25 = set()
    for line in closest(wd + "/q.bed", wd + "/atlas.bed"):
        p = line.split("\t")
        d = int(p[-1])
        if 0 <= d <= 25:
            matched25 += 1
            sites25.add((p[6], p[7], p[11]))
    # nearest DETECTED-atlas site per call, at 100 bp -> many-to-one structure
    site_hits = Counter()
    matched100 = 0
    for line in closest(wd + "/q.bed", wd + "/det.bed"):
        p = line.split("\t")
        d = int(p[-1])
        if 0 <= d <= 100:
            matched100 += 1
            site_hits[(p[6], p[7], p[11])] += 1
    shared = sum(c for c in site_hits.values() if c > 1)
    # self-spacing
    gaps = []
    for line in closest(wd + "/q.bed", wd + "/q.bed", extra=("-io", "-N")) if False else []:
        pass
    res = subprocess.run(["bedtools", "closest", "-s", "-d", "-io", "-t", "first",
                          "-a", wd + "/q.bed", "-b", wd + "/q.bed"],
                         capture_output=True, text=True, check=True, env=ENV)
    for line in res.stdout.splitlines():
        p = line.split("\t")
        d = int(p[-1])
        if d >= 0:
            gaps.append(d)
    gaps.sort()
    row = {
        "name": a.name, "n_calls": n, "matched@25": matched25,
        "distinct_atlas@25": len(sites25),
        "calls_per_atlas_site": (matched25 / len(sites25)) if sites25 else float("nan"),
        "matched@100": matched100, "recall_sites@100": len(site_hits),
        "many_to_one@100": (shared / matched100) if matched100 else float("nan"),
        "median_gap": gaps[len(gaps) // 2] if gaps else -1,
        "p10_gap": gaps[len(gaps) // 10] if gaps else -1,
        "frac_gap_lt100": (sum(1 for g in gaps if g < 100) / len(gaps)) if gaps else float("nan"),
    }
    subprocess.run(["rm", "-rf", wd])
    print("\t".join(str(v) for v in row.values()))
    sys.stderr.write("\t".join(row.keys()) + "\n")


if __name__ == "__main__":
    main()
