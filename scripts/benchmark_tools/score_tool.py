#!/usr/bin/env python3
"""
score_tool.py <pas.bed> <label> -- score ANY tool's PAS point BED exactly like
scripts/manuscript_figures/benchmark_curated.py scores the PeakATail arms.

The machinery below is adapted verbatim from benchmark_curated.py (same awk,
same bedtools invocations, same row schema) so that numbers are directly
comparable row-for-row with results/figures/manuscript/benchmark_curated.tsv.

WHAT IT COMPUTES (point mode, strand-matched `bedtools closest -s -d -t first`)
  precision  vs full curated PolyASite 2.0 rep sites + protein-coding TES,
             at 10/25/50/100/200 bp
  null       3 seeds of width-preserving `bedtools shuffle -chrom
             -noOverlapping -incl <merged gene bodies>`, scored identically
             (NOTE: for a 1-bp point input the shuffled widths are 1 bp; the
             PeakATail arms shuffled their full peak intervals then reduced to
             the 3' base -- equivalent nulls of 'random genic points')
  recall     (a) full atlas (569,005 rep sites)
             (b) SAME 14,851-detected-gene restricted atlas file
                 benchmark_curated.py used (285,220 sites; byte-checked) so
                 denominators are comparable across tools
             TES recall as sanity
  f1/roadmap same bars (P>=0.70, R>=0.60, F1>=0.65)

INPUT: BED6 whose intervals are the inferred PAS. '#' comments allowed.
  Intervals wider than 1 bp are reduced to the 3'-most base by strand
  ('+': end-1; '-': start), identical to benchmark_curated.py.

OUTPUT: <outdir>/score_<label>.tsv with the benchmark_curated.tsv columns
  (panel arm series reference cutoff_bp replicate n_query n_matched value bar
  passes); arm == <label>.

REFERENCES: default to the human GRCh38 set (unchanged behavior, byte-identical
  output).  For other genomes override with --atlas/--tes/--genome/--genebodies
  and pass --detected-atlas none (or a path) -- 'none' skips the
  detected-gene-restricted recall flavor (and its f1/roadmap rows) entirely.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
REF_DIR = WD / "data/references"
PAS2 = REF_DIR / "atlases/polyasite2.GRCh38.96.rep_sites.bed6"
TES = REF_DIR / "atlases/tes.protein_coding.GRCh38.99.bed6"
GENOME = REF_DIR / "chrom.sizes.nochr.filt"
GENE_END = REF_DIR / "gene_end.bed"
ANN = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/"
           "runs/grid/lg_annotate/annotatedpas.bed")
SHARED = WD / "results/benchmark_tools/shared_refs"
REF_DET = SHARED / "pas2.in_detected_genes.bed"
GENIC = SHARED / "genebodies.merged.bed"
GENES_TXT = SHARED / "detected_genes.txt"

CUTOFFS = [10, 25, 50, 100, 200]
N_SEEDS = 3
BAR_P, BAR_R, BAR_F1 = 0.70, 0.60, 0.65

ENV = dict(os.environ, LC_ALL="C")  # tr_TR locale silently corrupts BED sorting

AWK_HIST = (
    "awk -F'\\t' -v CUT=\"" + ",".join(map(str, CUTOFFS)) + "\" '"
    "BEGIN{nc=split(CUT,cc,\",\")}"
    "{d=$NF; n++; if(d>=0){for(i=1;i<=nc;i++) if(d<=cc[i]) h[i]++}}"
    "END{printf \"%d\", n; for(i=1;i<=nc;i++) printf \"\\t%d\", h[i]+0; printf \"\\n\"}'"
)


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


def make_point(iv, pt):
    """Collapse each interval to its 3'-most base by strand ('+': end-1; '-': start)."""
    sh("awk -F'\\t' 'BEGIN{OFS=\"\\t\"}"
       "{if($6==\"+\"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,$4,$5,$6}' "
       f"{iv} | sort -k1,1 -k2,2n > {pt}")


def closest_hist(a, b):
    out = sh_out(f"bedtools closest -s -d -t first -a {a} -b {b} 2>/dev/null "
                 f"| {AWK_HIST}").strip().split("\t")
    return int(out[0]), dict(zip(CUTOFFS, (int(x) for x in out[1:])))


def ensure_shared_refs():
    """Regenerate shared_refs exactly as benchmark_curated.py did, if missing."""
    SHARED.mkdir(parents=True, exist_ok=True)
    if not GENIC.exists():
        sh(f"sort -k1,1 -k2,2n {GENE_END} | bedtools merge -i - > {GENIC}")
    if not REF_DET.exists():
        det_bed = SHARED / "detected_gene_bodies.bed"
        sh(f"cut -f5 {ANN} | sort -u > {GENES_TXT}")
        sh(f"awk -F'\\t' 'NR==FNR{{g[$1]=1;next}} g[$4]' {GENES_TXT} {GENE_END} "
           f"| sort -k1,1 -k2,2n > {det_bed}")
        sh(f"bedtools intersect -a {PAS2} -b {det_bed} -s -u > {REF_DET}")
    n_genes = int(sh_out(f"wc -l < {GENES_TXT}").strip())
    n_det = int(sh_out(f"wc -l < {REF_DET}").strip())
    assert n_genes == 14851, f"expected 14,851 detected genes, got {n_genes}"
    assert n_det == 285220, f"expected 285,220 restricted atlas sites, got {n_det}"
    return n_det


def f1(p, r):
    return 0.0 if p + r == 0 else 2 * p * r / (p + r)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("pas_bed", help="tool's standard-form PAS BED6 (points)")
    ap.add_argument("label", help="tool label, becomes the 'arm' column")
    ap.add_argument("--outdir", default=None,
                    help="default: directory of pas_bed")
    ap.add_argument("--workdir", default=None,
                    help="scratch dir for intermediates (default: <outdir>/.work_<label>)")
    ap.add_argument("--atlas", default=None,
                    help=f"rep-site atlas BED6 (default: {PAS2})")
    ap.add_argument("--tes", default=None,
                    help=f"protein-coding TES BED6 (default: {TES})")
    ap.add_argument("--genome", default=None,
                    help=f"chrom sizes for filtering + shuffle (default: {GENOME})")
    ap.add_argument("--genebodies", default=None,
                    help=f"merged gene-body BED for the null shuffle "
                         f"(default: {GENIC}, auto-built)")
    ap.add_argument("--detected-atlas", default=None,
                    help=f"detected-gene-restricted atlas BED, or 'none' to "
                         f"skip that recall flavor (default: {REF_DET}, "
                         f"auto-built + count-asserted, human only)")
    a = ap.parse_args()

    # Safety guard: a non-default --atlas with --detected-atlas/--genebodies
    # left at their HUMAN defaults would silently score another species
    # against human references. Refuse unless overrides are explicit.
    if a.atlas is not None and (a.detected_atlas is None or a.genebodies is None):
        ap.error("--atlas overridden but --detected-atlas/--genebodies left at human "
                 "defaults; pass matching references (or --detected-atlas none) explicitly.")

    src = Path(a.pas_bed).resolve()
    label = a.label
    outdir = Path(a.outdir) if a.outdir else src.parent
    work = Path(a.workdir) if a.workdir else outdir / f".work_{label}"
    for d in (outdir, work):
        d.mkdir(parents=True, exist_ok=True)

    atlas = Path(a.atlas) if a.atlas else PAS2
    tes = Path(a.tes) if a.tes else TES
    genome = Path(a.genome) if a.genome else GENOME
    genic = Path(a.genebodies) if a.genebodies else GENIC
    if a.detected_atlas is None:
        det = REF_DET
    elif a.detected_atlas.lower() == "none":
        det = None
    else:
        det = Path(a.detected_atlas)

    # human shared_refs machinery (build + count assertions) only when the
    # human defaults are actually in play
    if genic == GENIC or det == REF_DET:
        ensure_shared_refs()
    for p in [atlas, tes, genome, genic] + ([det] if det else []):
        if not p.exists():
            sys.exit(f"missing reference file: {p}")
    refs = {"atlas_full": atlas, "tes": tes}
    if det is not None:
        refs["atlas_detected"] = det
    refs = {k: refs[k] for k in ("atlas_full", "atlas_detected", "tes")
            if k in refs}  # preserve original human iteration order

    # 1. genome-filter (drop chroms absent from GENOME), sort, reduce to point
    raw_n = int(sh_out(f"grep -vc '^#' {src} || true").strip())
    iv, pt = work / f"{label}.iv.bed", work / f"{label}.pt.bed"
    sh(f"grep -v '^#' {src} | "
       f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {genome} - "
       f"| sort -k1,1 -k2,2n > {iv}")
    make_point(iv, pt)
    n = int(sh_out(f"wc -l < {pt}").strip())
    print(f"[real] {label}: n={n:,} points ({raw_n - n} off-genome peaks dropped)")

    # 2. nulls: width-preserving shuffle inside merged annotated gene bodies
    null_pt = {}
    for seed in range(1, N_SEEDS + 1):
        siv = work / f"{label}.null.s{seed}.iv.bed"
        spt = work / f"{label}.null.s{seed}.pt.bed"
        sh(f"bedtools shuffle -i {iv} -g {genome} -chrom "
           f"-noOverlapping -maxTries 5000 -seed {seed} -incl {genic} "
           f"2>/dev/null | sort -k1,1 -k2,2n > {siv}")
        make_point(siv, spt)
        null_pt[seed] = spt
    print(f"[null] {label}: {N_SEEDS} seeds shuffled")

    rows = []

    def record(**kw):
        rows.append(kw)

    # 3. precision: real vs atlas_full + tes; null vs atlas_full
    P = {}
    for refkey in ("atlas_full", "tes"):
        nq, hits = closest_hist(pt, refs[refkey])
        P[refkey] = {c: hits[c] / nq for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="precision", arm=label, series="real", reference=refkey,
                   cutoff_bp=c, replicate=0, n_query=nq, n_matched=hits[c],
                   value=hits[c] / nq)
    nullp = {}
    for seed in range(1, N_SEEDS + 1):
        nq, hits = closest_hist(null_pt[seed], atlas)
        nullp[seed] = {c: hits[c] / nq for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="precision", arm=label, series="null_genic",
                   reference="atlas_full", cutoff_bp=c, replicate=seed,
                   n_query=nq, n_matched=hits[c], value=hits[c] / nq)
    print(f"[prec] {label}: @100bp real {P['atlas_full'][100]:.4f} "
          f"(TES {P['tes'][100]:.4f}), null "
          + "/".join(f"{nullp[s][100]:.4f}" for s in range(1, N_SEEDS + 1)))

    # 4. recall: both atlas flavors + TES sanity
    R = {}
    for refkey, refp in refs.items():
        nq, hits = closest_hist(refp, pt)
        R[refkey] = {c: hits[c] / nq for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="recall", arm=label, series="real", reference=refkey,
                   cutoff_bp=c, replicate=0, n_query=nq, n_matched=hits[c],
                   value=hits[c] / nq)
    det_msg = (f"detected {R['atlas_detected'][100]:.4f}"
               if "atlas_detected" in R else "detected SKIPPED")
    print(f"[recall] {label}: @100bp full {R['atlas_full'][100]:.4f} | "
          f"{det_msg} | TES {R['tes'][100]:.4f}")

    # 5. F1 + roadmap pass/fail (F1 pairs atlas_full precision w/ each recall flavor)
    for flav in [f for f in ("atlas_full", "atlas_detected") if f in R]:
        for c in CUTOFFS:
            v = f1(P["atlas_full"][c], R[flav][c])
            record(panel="f1", arm=label, series="real", reference=flav,
                   cutoff_bp=c, replicate=0, n_query=np.nan, n_matched=np.nan,
                   value=v)
            for metric, val, bar in (("precision", P["atlas_full"][c], BAR_P),
                                     ("recall", R[flav][c], BAR_R),
                                     ("f1", v, BAR_F1)):
                record(panel="roadmap", arm=label, series=metric, reference=flav,
                       cutoff_bp=c, replicate=0, n_query=np.nan,
                       n_matched=np.nan, value=val, bar=bar,
                       passes=int(val >= bar))

    record(panel="meta", arm=label, series="n_points", reference="",
           cutoff_bp=np.nan, replicate=0, n_query=raw_n, n_matched=n,
           value=raw_n - n)

    df = pd.DataFrame(rows)
    out = outdir / f"score_{label}.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6f")
    print(f"[out] {out}  ({len(df)} rows)")


if __name__ == "__main__":
    sys.exit(main())
