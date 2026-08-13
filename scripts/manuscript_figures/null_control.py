#!/usr/bin/env python3
"""
null_control.py -- NULL CONTROL for the PolyASite-3.0 atlas benchmark.

The pipeline reports precision ~0.9986 at a 100 bp cutoff for called PAS against
/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed
(18,432,135 entries).  With a reference that dense, is high precision evidence
of accuracy, or an artifact of reference density?  This script quantifies it.

Design
------
REAL sets   runs/grid/lg_annotate/pasbed.bed (primary, n=22,633)
            runs/B1_cohort_full/pasbed.bed   (secondary, n=16,500)

NULL sets   `bedtools shuffle -chrom -noOverlapping` -> identical n, identical
            per-chromosome counts, identical interval widths, strand carried
            over.  Two null models:
              N_genome  anywhere on the same chromosome
              N_genic   anywhere inside an annotated gene body on the same
                        chromosome (-incl merged gene_end.bed).  The harsher,
                        biologically fair null: a PAS caller only ever sees
                        transcribed DNA.
            5 replicates (primary) / 3 replicates (secondary).

MATCHING    bedtools closest -s -d -t first  (strand-matched nearest atlas
            site; distance 0 = overlap).
              mode "interval"  whole called peak -- the geometry the pipeline
                               scores
              mode "point"     the peak's 3'-most base by strand, i.e. the
                               actual inferred cleavage site
            NOTE, and this matters for every sentence written about these
            numbers: this is NOT the pipeline's own matcher.
            tools/PeakATail/ema/benchmark/metrics.py scores with
              BedTool(predicted).window(reference, w=cutoff)
            -- a symmetric +/-cutoff window around the whole called interval
            with NO strand requirement.  Adding -s is a deliberate tightening
            here (it makes the null harder to beat and is the right control),
            but it means the numbers in this figure are systematically a shade
            below the pipeline's published ones and must never be described as
            "the pipeline's" or "as benchmarked".  Section 6b reconciles the
            two explicitly and both appear in the tsv.
ATLAS TIERS full (18.43 M) | col10 class in {TE,AL,EX} | col5 average TPM >= 1
RECALL      recomputed on a reference restricted to atlas sites inside the
            strand-matched gene body of a gene detected in this cohort.
CLASS TEST  composition (atlas col10) of the class of the nearest matched atlas
            site -- real vs null vs the atlas's own background composition.

Everything streams; the 1.2 GB atlas is never loaded into pandas.
Distance histograms are cached under results/figures/manuscript/.cache_null_control/
so re-plotting is cheap; delete that directory to force a full recompute.
"""

import os
import subprocess
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- manuscript deliverables: PNG + vector PDF into manuscript/figures/ ---
import os as _os
import matplotlib as _mpl
_mpl.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (journal requirement)
_mpl.rcParams['ps.fonttype'] = 42
FIGDIR = '/mnt/ssd1/Projects/PeakATail_wd/manuscript/figures'
_os.makedirs(FIGDIR, exist_ok=True)
def save_manuscript(fig, name, **kw):
    for ext in ('png', 'pdf'):
        p = _os.path.join(FIGDIR, f'{name}.{ext}')
        fig.savefig(p, **({'dpi': 300} if ext == 'png' else {}), **kw)
        print('wrote', p)


# ----------------------------------------------------------------------------
# paths / constants
# ----------------------------------------------------------------------------
RUNS = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs")
ATLAS = Path("/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed")
GENOME = Path("/mnt/ssd1/Projects/PeakATail_wd/data/references/chrom.sizes.nochr.filt")
GENE_END = Path("/mnt/ssd1/Projects/PeakATail_wd/data/references/gene_end.bed")

OUTDIR = Path("/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript")
NAME = "null_control"
CACHE = OUTDIR / f".cache_{NAME}"
WORK = Path(
    "/tmp/claude-1000/-mnt-ssd1-Projects-PeakATail-wd/"
    "93f9b511-1339-4cbc-9fdb-148fd4b2f08f/scratchpad/nullctl_run"
)

CUTOFFS = [50, 100, 200, 500, 1000]
CLASS_CUTOFF = 100
CLASSES = ["TE", "AL", "EX", "IN", "DI", "UI"]
CLASS_DESC = {"TE": "terminal exon", "AL": "any last exon", "EX": "exonic",
              "IN": "intronic", "DI": "downstream", "UI": "upstream/intergenic"}
N_REP_PRIMARY = 5
N_REP_SECONDARY = 3

ENV = dict(os.environ, LC_ALL="C")  # a Turkish locale would corrupt sort order

# validated palette -- one colour per ENTITY, fixed across every panel
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"
COL = {
    "real": C[0],           # real called PAS
    "N_genome": C[1],       # null, shuffled on chromosome
    "N_genic": C[2],        # null, shuffled inside gene bodies
    "atlas": C[3],          # the atlas reference itself
    "ref_full": C[3],
    "ref_detected": C[4],   # atlas within detected gene bodies
    "ref_te": C[5],         # + terminal-exon site classes
    "ref_tpm1": INK,        # + average TPM >= 1 (strictest)
}
LABEL = {"real": "Real called PAS",
         "N_genome": "Null: shuffled on chromosome",
         "N_genic": "Null: shuffled in gene bodies",
         "atlas": "Atlas background"}
SHORT = {"real": "Real called PAS",
         "N_genome": "shuffled on chromosome",
         "N_genic": "shuffled in gene bodies"}


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


def log(*a):
    print(*a, flush=True)


def cached(key: str, cmd: str) -> str:
    """Run a streaming bedtools|awk pipeline once; memoise its small stdout."""
    f = CACHE / f"{key}.txt"
    if not f.exists():
        f.write_text(sh_out(cmd))
    return f.read_text()


# ----------------------------------------------------------------------------
# 0. workspace + sanity checks
# ----------------------------------------------------------------------------
for d in (WORK, OUTDIR, CACHE):
    d.mkdir(parents=True, exist_ok=True)

genome_rows = [l.split("\t") for l in GENOME.read_text().strip().split("\n")]
genome_chroms = {r[0] for r in genome_rows}
GENOME_BP = sum(int(r[1]) for r in genome_rows)
atlas_chroms = set(cached("atlas_chroms", f"cut -f1 {ATLAS} | uniq -c"
                          ).split()[1::2])
assert atlas_chroms <= genome_chroms, (
    f"genome file does not cover atlas chrom naming: "
    f"{sorted(atlas_chroms - genome_chroms)}")
N_ATLAS = int(cached("atlas_n", f"wc -l < {ATLAS}").strip())
log(f"[check] genome file matches atlas naming (ensembl); {len(genome_chroms)} "
    f"chroms, {GENOME_BP/1e9:.3f} Gb; atlas entries {N_ATLAS:,}")


# ----------------------------------------------------------------------------
# 1. real PAS sets
# ----------------------------------------------------------------------------
def make_point(iv: Path, pt: Path):
    """Collapse each peak to its 3'-most base by strand = inferred cleavage site."""
    sh("awk -F'\\t' 'BEGIN{OFS=\"\\t\"}"
       "{if($6==\"+\"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,$4,$5,$6}' "
       f"{iv} | sort -k1,1 -k2,2n > {pt}")


def prep_real(src: Path, tag: str):
    raw_n = int(sh_out(f"wc -l < {src}").strip())
    iv, pt = WORK / f"{tag}.iv.bed", WORK / f"{tag}.pt.bed"
    if not iv.exists():
        sh(f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {GENOME} {src} "
           f"| sort -k1,1 -k2,2n > {iv}")
        make_point(iv, pt)
    n = int(sh_out(f"wc -l < {iv}").strip())
    w = sh_out("awk '{s+=$3-$2; if($3-$2>mx)mx=$3-$2}END"
               "{printf \"%.1f %d\", s/NR, mx}' " + str(iv)).split()
    return dict(tag=tag, raw_n=raw_n, n=n, dropped=raw_n - n, iv=iv, pt=pt,
                mean_width=float(w[0]), max_width=int(w[1]))


REAL = {"lg_annotate": prep_real(RUNS / "grid/lg_annotate/pasbed.bed", "real_lg"),
        "B1_cohort_full": prep_real(RUNS / "B1_cohort_full/pasbed.bed", "real_b1")}
for k, v in REAL.items():
    log(f"[real] {k}: n={v['n']:,} ({v['dropped']} off-genome scaffold peaks "
        f"dropped), mean width {v['mean_width']:.0f} bp, max {v['max_width']:,} bp")


# ----------------------------------------------------------------------------
# 2. null sets
# ----------------------------------------------------------------------------
GENIC = WORK / "genebodies.merged.bed"
if not GENIC.exists():
    sh(f"sort -k1,1 -k2,2n {GENE_END} | bedtools merge -i - > {GENIC}")
GENIC_BP = int(sh_out("awk '{s+=$3-$2}END{print s}' " + str(GENIC)).strip())
log(f"[null] merged gene bodies span {GENIC_BP/1e9:.3f} Gb "
    f"({100*GENIC_BP/GENOME_BP:.1f}% of the genome)")


def shuffle_set(real_iv: Path, tag: str, model: str, seed: int):
    out_iv = WORK / f"{tag}.{model}.s{seed}.iv.bed"
    out_pt = WORK / f"{tag}.{model}.s{seed}.pt.bed"
    incl = f"-incl {GENIC}" if model == "N_genic" else ""
    if not out_pt.exists():
        sh(f"bedtools shuffle -i {real_iv} -g {GENOME} -chrom -noOverlapping "
           f"-maxTries 5000 -seed {seed} {incl} 2>/dev/null "
           f"| sort -k1,1 -k2,2n > {out_iv}")
        make_point(out_iv, out_pt)
    return out_iv, out_pt


# ----------------------------------------------------------------------------
# 3. atlas tiers
# ----------------------------------------------------------------------------
te_bed, tpm_bed = WORK / "atlas.te_al_ex.bed", WORK / "atlas.tpm1.bed"
if not te_bed.exists():
    sh(f"awk -F'\\t' '$10==\"TE\"||$10==\"AL\"||$10==\"EX\"' {ATLAS} > {te_bed}")
if not tpm_bed.exists():
    sh(f"awk -F'\\t' '($5+0)>=1' {ATLAS} > {tpm_bed}")
TIERS = {
    "full": (ATLAS, N_ATLAS, "every atlas entry"),
    "te_al_ex": (te_bed, int(sh_out(f"wc -l < {te_bed}").strip()),
                 "col10 class in {TE, AL, EX}"),
    "tpm1": (tpm_bed, int(sh_out(f"wc -l < {tpm_bed}").strip()),
             "col5 average TPM >= 1"),
}
for k, (_, n, d) in TIERS.items():
    log(f"[atlas tier] {k}: {n:,} sites -- {d}")


# ----------------------------------------------------------------------------
# 4. precision engine
# ----------------------------------------------------------------------------
AWK_HIST = (
    "awk -F'\\t' -v CUT=\"" + ",".join(map(str, CUTOFFS)) + "\" '"
    "BEGIN{nc=split(CUT,cc,\",\")}"
    "{d=$NF; n++; if(d>=0){for(i=1;i<=nc;i++) if(d<=cc[i]) h[i]++}}"
    "END{printf \"%d\", n; for(i=1;i<=nc;i++) printf \"\\t%d\", h[i]+0; printf \"\\n\"}'"
)
AWK_CLASS = (
    "awk -F'\\t' -v CUT=" + str(CLASS_CUTOFF) + " '"
    "{d=$NF; if(d>=0 && d<=CUT){c[$16]++; n++}}"
    "END{for(k in c) printf \"%s\\t%d\\t%d\\n\", k, c[k], n}'"
)

rows = []


def record(**kw):
    rows.append(kw)


def precision(query: Path, tier: str):
    ref = TIERS[tier][0]
    out = cached(f"prec__{query.stem}__{tier}",
                 f"bedtools closest -s -d -t first -a {query} -b {ref} "
                 f"2>/dev/null | {AWK_HIST}").strip().split("\t")
    n, hits = int(out[0]), [int(x) for x in out[1:]]
    return n, dict(zip(CUTOFFS, hits))


def match_classes(query: Path):
    """Composition of the col10 site class of the nearest atlas match <=100 bp."""
    out = cached(f"class__{query.stem}",
                 f"bedtools closest -s -d -t first -a {query} -b {ATLAS} "
                 f"2>/dev/null | {AWK_CLASS}")
    counts, tot = {}, 0
    for line in out.strip().split("\n"):
        k, c, n = line.split("\t")
        counts[k] = int(c)
        tot = int(n)
    return tot, counts


for dataset, meta in REAL.items():
    nrep = N_REP_PRIMARY if dataset == "lg_annotate" else N_REP_SECONDARY
    tiers_here = list(TIERS) if dataset == "lg_annotate" else ["full"]

    for mode, q in (("interval", meta["iv"]), ("point", meta["pt"])):
        for tier in tiers_here:
            n, hits = precision(q, tier)
            for c in CUTOFFS:
                record(panel="precision", dataset=dataset, series="real",
                       mode=mode, atlas_tier=tier, cutoff_bp=c, replicate=0,
                       n_query=n, n_matched=hits[c], value=hits[c] / n)
        log(f"[prec] {dataset} real {mode:8s} full  "
            + " ".join(f"{c}:{precision(q,'full')[1][c]/meta['n']:.4f}"
                       for c in CUTOFFS))

    for model in ("N_genome", "N_genic"):
        for rep in range(1, nrep + 1):
            siv, spt = shuffle_set(meta["iv"], meta["tag"], model, rep)
            for mode, q in (("interval", siv), ("point", spt)):
                for tier in tiers_here:
                    n, hits = precision(q, tier)
                    for c in CUTOFFS:
                        record(panel="precision", dataset=dataset, series=model,
                               mode=mode, atlas_tier=tier, cutoff_bp=c,
                               replicate=rep, n_query=n, n_matched=hits[c],
                               value=hits[c] / n)
            log(f"[prec] {dataset} {model} rep{rep} done")

# ---- class composition of the nearest match (primary dataset, point mode) ---
tot, cnt = match_classes(REAL["lg_annotate"]["pt"])
for cl in CLASSES:
    record(panel="match_class", dataset="lg_annotate", series="real",
           mode="point", atlas_tier=f"class_{cl}", cutoff_bp=CLASS_CUTOFF,
           replicate=0, n_query=tot, n_matched=cnt.get(cl, 0),
           value=cnt.get(cl, 0) / tot)
for model in ("N_genome", "N_genic"):
    for rep in range(1, N_REP_PRIMARY + 1):
        _, spt = shuffle_set(REAL["lg_annotate"]["iv"], "real_lg", model, rep)
        tot, cnt = match_classes(spt)
        for cl in CLASSES:
            record(panel="match_class", dataset="lg_annotate", series=model,
                   mode="point", atlas_tier=f"class_{cl}",
                   cutoff_bp=CLASS_CUTOFF, replicate=rep, n_query=tot,
                   n_matched=cnt.get(cl, 0), value=cnt.get(cl, 0) / tot)
bg = cached("atlas_class_bg",
            f"awk -F'\\t' '{{c[$10]++;n++}}END{{for(k in c) "
            f"printf \"%s\\t%d\\t%d\\n\", k, c[k], n}}' {ATLAS}")
bgc, bgn = {}, 0
for line in bg.strip().split("\n"):
    k, c, n = line.split("\t")
    bgc[k], bgn = int(c), int(n)
for cl in CLASSES:
    record(panel="match_class", dataset="atlas", series="atlas", mode="background",
           atlas_tier=f"class_{cl}", cutoff_bp=CLASS_CUTOFF, replicate=0,
           n_query=bgn, n_matched=bgc.get(cl, 0), value=bgc.get(cl, 0) / bgn)
log("[class] match-class composition done")


# ----------------------------------------------------------------------------
# 5. reference-density control
#    bp of the genome within d of ANY strand-matched atlas site, per strand
# ----------------------------------------------------------------------------
DENS_AWK = (
    "awk -F'\\t' -v CUT=\"" + ",".join(map(str, CUTOFFS)) + "\" '"
    "BEGIN{nc=split(CUT,cc,\",\")}"
    "{c=$1; for(i=1;i<=nc;i++){s=$2-cc[i]; if(s<0)s=0; e=$3+cc[i];"
    "  if(c!=pc[i]){ if(pc[i]!=\"\") tot[i]+=pe[i]-ps[i]; pc[i]=c; ps[i]=s; pe[i]=e }"
    "  else if(s>pe[i]){ tot[i]+=pe[i]-ps[i]; ps[i]=s; pe[i]=e }"
    "  else if(e>pe[i]) pe[i]=e }}"
    "END{for(i=1;i<=nc;i++){tot[i]+=pe[i]-ps[i]; printf \"%s\\t%d\\n\", cc[i], tot[i]}}'"
)
for strand, sym in (("plus_strand", "+"), ("minus_strand", "-")):
    out = cached(f"density__{strand}",
                 f"awk -F'\\t' '$6==\"{sym}\"' {ATLAS} | {DENS_AWK}")
    for line in out.strip().split("\n"):
        c, bp = line.split("\t")
        record(panel="reference_density", dataset="atlas", series=strand,
               mode="genome_coverage", atlas_tier="full", cutoff_bp=int(c),
               replicate=0, n_query=GENOME_BP, n_matched=int(bp),
               value=int(bp) / GENOME_BP)
    log(f"[density] {sym} strand done")


# ----------------------------------------------------------------------------
# 6. corrected recall on a cohort-restricted reference
#
# RESTRICTION, stated exactly:
#   detected genes := unique Ensembl gene IDs in column 5 of
#     runs/grid/lg_annotate/annotatedpas.bed -- genes to which this cohort's run
#     assigned at least one PAS.
#   gene bodies    := rows of data/references/gene_end.bed whose gene ID is in
#     that set (full gene span, strand kept).
#   restricted ref := atlas sites overlapping such a gene body ON THE SAME
#     STRAND (bedtools intersect -s -u).  Two further nestings add the col10
#     TE/AL/EX class filter, then the col5 average-TPM >= 1 filter.
# Recall := fraction of reference sites whose strand-matched nearest called peak
#   is within the cutoff.  Reported in interval mode (the geometry the pipeline
#   uses) AND in point mode, because interval mode inflates recall for exactly
#   the reason it inflates precision: an 821 bp mean peak swallows atlas sites
#   whole.  Both go in the tsv; the point-mode number is the honest one.
# ----------------------------------------------------------------------------
ann = RUNS / "grid/lg_annotate/annotatedpas.bed"
genes_txt, det_bed = WORK / "detected_genes.txt", WORK / "detected_gene_bodies.bed"
if not det_bed.exists():
    sh(f"cut -f5 {ann} | sort -u > {genes_txt}")
    sh(f"awk -F'\\t' 'NR==FNR{{g[$1]=1;next}} g[$4]' {genes_txt} {GENE_END} "
       f"| sort -k1,1 -k2,2n > {det_bed}")
n_detected = int(sh_out(f"wc -l < {genes_txt}").strip())
n_gtf_genes = int(sh_out(f"cut -f4 {GENE_END} | sort -u | wc -l").strip())
det_bp = int(sh_out("sort -k1,1 -k2,2n " + str(det_bed) +
                    " | bedtools merge -i - | awk '{s+=$3-$2}END{print s}'").strip())
log(f"[recall] detected genes {n_detected:,} / {n_gtf_genes:,} annotated "
    f"({100*n_detected/n_gtf_genes:.0f}%), gene bodies span {det_bp/1e9:.3f} Gb")

ref_det = WORK / "atlas.in_detected_genes.bed"
ref_det_te = WORK / "atlas.in_detected_genes.te.bed"
ref_det_tpm = WORK / "atlas.in_detected_genes.tpm1.bed"
if not ref_det.exists():
    sh(f"bedtools intersect -a {ATLAS} -b {det_bed} -s -u > {ref_det}")
if not ref_det_te.exists():
    sh(f"awk -F'\\t' '$10==\"TE\"||$10==\"AL\"||$10==\"EX\"' {ref_det} > {ref_det_te}")
if not ref_det_tpm.exists():
    sh(f"awk -F'\\t' '($5+0)>=1' {ref_det} > {ref_det_tpm}")

REFS = {
    "ref_full": (ATLAS, "Full atlas"),
    "ref_detected": (ref_det, "In detected gene bodies"),
    "ref_te": (ref_det_te, "+ TE/AL/EX classes"),
    "ref_tpm1": (ref_det_tpm, "+ average TPM \u2265 1"),
}


def nlab(n):
    """Legend n, derived from the data rather than hard-coded."""
    return f"{n/1e6:.3g} M" if n >= 1e6 else (f"{n/1e3:.0f} k" if n >= 1e3
                                              else str(n))


REF_N = {}
real_iv, real_pt = REAL["lg_annotate"]["iv"], REAL["lg_annotate"]["pt"]
for key, (refp, desc) in REFS.items():
    for mode, target in (("interval", real_iv), ("point", real_pt)):
        pfx = "recall__" if mode == "interval" else "recall_point__"
        out = cached(f"{pfx}{key}",
                     f"bedtools closest -s -d -t first -a {refp} -b {target} "
                     f"2>/dev/null | {AWK_HIST}").strip().split("\t")
        n, hits = int(out[0]), [int(x) for x in out[1:]]
        REF_N[key] = n
        for c, h in zip(CUTOFFS, hits):
            record(panel="recall", dataset="lg_annotate", series=key, mode=mode,
                   atlas_tier=key, cutoff_bp=c, replicate=0, n_query=n,
                   n_matched=h, value=h / n)
        log(f"[recall] {key:13s} {mode:8s} n_ref={n:>10,}  "
            + " ".join(f"{c}:{h/n:.4f}" for c, h in zip(CUTOFFS, hits)))
REF_N_LAB = {k: nlab(v) for k, v in REF_N.items()}


# ----------------------------------------------------------------------------
# 6b. AUDIT: what the pipeline actually measured
#
# tools/PeakATail/ema/benchmark/metrics.py:
#     hits     = predicted.window(reference, w=cutoff)   # NO -sm / -sw
#     hits_rev = reference.window(predicted, w=cutoff)
# -> symmetric window around the whole interval, strand-agnostic, denominators
# 22,633 predicted and 18,432,135 reference (scaffold peaks included).
# Re-running `closest` WITHOUT -s on the pipeline's own denominators reproduces
# runs/grid/lg_annotate/benchmark_vs_polyasite_v3.json exactly, which is what
# licenses the sentence "the strand requirement, not the scaffold peaks, is why
# this figure reads 0.9985 where the pipeline reads 0.9986".
# ----------------------------------------------------------------------------
import json  # noqa: E402

pj = json.loads((RUNS / "grid/lg_annotate/benchmark_vs_polyasite_v3.json").read_text())
for c in CUTOFFS:
    d = pj["cutoffs"][str(c)]
    record(panel="pipeline_audit", dataset="lg_annotate", series="pipeline_json",
           mode="window_nostrand_precision", atlas_tier="full", cutoff_bp=c,
           replicate=0, n_query=pj["n_predicted"],
           n_matched=d["matched_predicted"],
           value=d["matched_predicted"] / pj["n_predicted"])
    record(panel="pipeline_audit", dataset="lg_annotate", series="pipeline_json",
           mode="window_nostrand_recall", atlas_tier="full", cutoff_bp=c,
           replicate=0, n_query=pj["n_reference"],
           n_matched=d["matched_reference"],
           value=d["matched_reference"] / pj["n_reference"])

RAW_N = REAL["lg_annotate"]["raw_n"]           # 22,633, scaffold peaks included
out = cached("audit_prec_nostrand",
             f"bedtools closest -d -t first -a {real_iv} -b {ATLAS} "
             f"2>/dev/null | {AWK_HIST}").strip().split("\t")
ns_hits = dict(zip(CUTOFFS, [int(x) for x in out[1:]]))
out = cached("audit_recall_nostrand",
             f"bedtools closest -d -t first -a {ATLAS} -b {real_iv} "
             f"2>/dev/null | {AWK_HIST}").strip().split("\t")
ns_rec_n, ns_rec = int(out[0]), dict(zip(CUTOFFS, [int(x) for x in out[1:]]))
for c in CUTOFFS:
    record(panel="pipeline_audit", dataset="lg_annotate", series="reproduced",
           mode="window_nostrand_precision", atlas_tier="full", cutoff_bp=c,
           replicate=0, n_query=RAW_N, n_matched=ns_hits[c],
           value=ns_hits[c] / RAW_N)
    record(panel="pipeline_audit", dataset="lg_annotate", series="reproduced",
           mode="window_nostrand_recall", atlas_tier="full", cutoff_bp=c,
           replicate=0, n_query=ns_rec_n, n_matched=ns_rec[c],
           value=ns_rec[c] / ns_rec_n)
PIPE_P100 = pj["cutoffs"]["100"]["matched_predicted"] / pj["n_predicted"]
PIPE_R100 = pj["cutoffs"]["100"]["matched_reference"] / pj["n_reference"]
for c in CUTOFFS:
    ok_p = ns_hits[c] == pj["cutoffs"][str(c)]["matched_predicted"]
    ok_r = ns_rec[c] == pj["cutoffs"][str(c)]["matched_reference"]
    log(f"[audit] @{c:>4}bp  precision matched {ns_hits[c]} vs json "
        f"{pj['cutoffs'][str(c)]['matched_predicted']} {'OK' if ok_p else 'DIFF'} | "
        f"recall matched {ns_rec[c]} vs json "
        f"{pj['cutoffs'][str(c)]['matched_reference']} {'OK' if ok_r else 'DIFF'}")


# ----------------------------------------------------------------------------
# 7. tidy table
# ----------------------------------------------------------------------------
df = pd.DataFrame(rows)
df.to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False, float_format="%.6f")
log(f"[out] {OUTDIR / (NAME + '.tsv')}  ({len(df)} rows)")


def get(panel, series, mode, tier, dataset="lg_annotate"):
    s = df[(df.panel == panel) & (df.dataset == dataset) & (df.series == series)
           & (df["mode"] == mode) & (df.atlas_tier == tier)]
    g = s.groupby("cutoff_bp")["value"]
    return (g.mean().reindex(CUTOFFS).values, g.min().reindex(CUTOFFS).values,
            g.max().reindex(CUTOFFS).values)


i100 = CUTOFFS.index(100)
p_real_iv = get("precision", "real", "interval", "full")[0]
p_ng_iv = get("precision", "N_genome", "interval", "full")[0]
p_gc_iv = get("precision", "N_genic", "interval", "full")[0]
p_real_pt = get("precision", "real", "point", "full")[0]
p_ng_pt = get("precision", "N_genome", "point", "full")[0]
p_gc_pt = get("precision", "N_genic", "point", "full")[0]


# ----------------------------------------------------------------------------
# 8. figure
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 7.8, "axes.titlesize": 8.6,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 6.8,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.linewidth": 0.8, "font.family": "DejaVu Sans",
})
fig, axes = plt.subplots(2, 3, figsize=(11.0, 7.45))
ax_a, ax_b, ax_c = axes[0]
ax_d, ax_e, ax_f = axes[1]
X = np.arange(len(CUTOFFS))


def style(ax, ylab, xlab):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)


def frac_axis(ax):
    ax.set_ylim(0, 1.06)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_yticklabels([f"{v:.1f}" for v in np.arange(0, 1.01, 0.2)])


def cutoff_axis(ax, xlab):
    ax.set_xticks(X)
    ax.set_xticklabels([str(c) for c in CUTOFFS])
    ax.set_xlabel(xlab)


# --- A. reference density ---------------------------------------------------
dens = {s: df[(df.panel == "reference_density") & (df.series == s)]
        .set_index("cutoff_bp")["value"].reindex(CUTOFFS).values
        for s in ("plus_strand", "minus_strand")}
dm = (dens["plus_strand"] + dens["minus_strand"]) / 2
ax_a.fill_between(X, np.minimum(*dens.values()), np.maximum(*dens.values()),
                  color=COL["atlas"], alpha=0.18, lw=0)
ax_a.plot(X, dm, color=COL["atlas"], lw=1.8, marker="o", ms=6, mec=SURFACE,
          mew=0.9, zorder=3)
ax_a.text(0.04, 0.80, "Genome within reach of a\nstrand-matched atlas site",
          transform=ax_a.transAxes, color=COL["atlas"], fontsize=6.8,
          fontweight="bold", va="top", ha="left")
for xi, v in zip(X, dm):
    ax_a.text(xi, v + 0.032, f"{100*v:.0f}%", ha="center", va="bottom",
              fontsize=6.6, color=COL["atlas"], fontweight="bold")
frac_axis(ax_a)
cutoff_axis(ax_a, "Distance tolerance (bp)")
style(ax_a, "Fraction of the 3.09 Gb genome", "Distance tolerance (bp)")
ax_a.set_title(f"{100*dm[i100]:.0f}% of the genome lies within 100 bp\n"
               f"of a strand-matched atlas site")
ax_a.text(0.03, 0.04,
          f"The PolyASite 3.0 file used here has {N_ATLAS/1e6:.2f} M entries,\n"
          f"one every {GENOME_BP/N_ATLAS:.0f} bp of genome. Line is the mean of\n"
          f"the two strands; the shaded band spans + vs − strand.",
          transform=ax_a.transAxes, fontsize=6.2, color=MUTED, va="bottom",
          linespacing=1.35)


# --- B/C. precision curves --------------------------------------------------
YLAB = {"interval": "Precision over the whole peak interval",
        "point": "Precision scoring only the peak 3′ base"}


def curves(ax, mode, offsets):
    for series in ("real", "N_genome", "N_genic"):
        m, lo, hi = get("precision", series, mode, "full")
        ax.fill_between(X, lo, hi, color=COL[series], alpha=0.18, lw=0)
        ax.plot(X, m, color=COL[series], lw=1.8, marker="o", ms=6, mec=SURFACE,
                mew=0.9, zorder=3, label=LABEL[series])
        ax.annotate(SHORT[series], (X[-1], m[-1] + offsets[series]), ha="right",
                    va="bottom", color=COL[series], fontsize=6.8,
                    fontweight="bold")
    frac_axis(ax)
    cutoff_axis(ax, "Match cutoff (bp to nearest atlas site)")
    style(ax, YLAB[mode], "Match cutoff (bp to nearest atlas site)")


curves(ax_b, "interval", {"real": 0.035, "N_genome": -0.105, "N_genic": 0.035})
ax_b.set_title(f"Shuffled peaks already score {p_ng_iv[i100]:.2f}\u2013"
               f"{p_gc_iv[i100]:.2f} when\nthe whole peak interval is scored")
ax_b.text(0.03, 0.04,
          f"A called peak is a wide interval (mean "
          f"{REAL['lg_annotate']['mean_width']:.0f} bp, max "
          f"{REAL['lg_annotate']['max_width']/1000:.1f} kb);\n"
          f"distance 0 only means the interval contains an atlas site.",
          transform=ax_b.transAxes, fontsize=6.2, color=MUTED, va="bottom")
ax_b.legend(loc="lower left", bbox_to_anchor=(0.02, 0.18), frameon=False,
            handlelength=1.6, borderpad=0.2, labelspacing=0.3)

curves(ax_c, "point", {"real": 0.035, "N_genome": -0.105, "N_genic": 0.035})
ax_c.set_title(f"Scoring the peak 3\u2032 base alone drops the null\n"
               f"to {p_ng_pt[i100]:.2f}\u2013{p_gc_pt[i100]:.2f} while real "
               f"PAS hold {p_real_pt[i100]:.3f}")
ax_c.text(0.99, 0.28,
          "The peak's 3\u2032-most base by strand is the\n"
          "inferred cleavage site itself.",
          transform=ax_c.transAxes, fontsize=6.2, color=MUTED, va="bottom",
          ha="right")
ax_c.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02), frameon=False,
            handlelength=1.6, borderpad=0.2, labelspacing=0.3)


# --- D. precision at 100 bp by atlas stringency (point mode) ----------------
tier_order = ["full", "te_al_ex", "tpm1"]
tier_lab = [f"Full atlas\n{TIERS['full'][1]/1e6:.2f} M",
            f"TE/AL/EX class\n{TIERS['te_al_ex'][1]/1e6:.2f} M",
            f"Avg TPM \u2265 1\n{TIERS['tpm1'][1]/1e3:.0f} k"]
xs, w = np.arange(len(tier_order)), 0.26
for j, series in enumerate(("real", "N_genome", "N_genic")):
    vals = np.array([get("precision", series, "point", t)[0][i100]
                     for t in tier_order])
    los = np.array([get("precision", series, "point", t)[1][i100] for t in tier_order])
    his = np.array([get("precision", series, "point", t)[2][i100] for t in tier_order])
    pos = xs + (j - 1) * (w + 0.03)
    ax_d.bar(pos, vals, width=w, color=COL[series], label=LABEL[series],
             zorder=3, edgecolor=SURFACE, linewidth=0.6)
    if series != "real":
        ax_d.errorbar(pos, vals, yerr=[vals - los, his - vals], fmt="none",
                      ecolor=INK, elinewidth=0.8, capsize=2, zorder=4)
    for p, v in zip(pos, vals):
        # never print a bar as "0.00" when it is not zero; stagger the two
        # 4-decimal labels so they cannot collide above near-zero bars
        tiny = v < 0.01
        ax_d.text(p, v + 0.02 + (0.06 if tiny and series == "N_genome" else 0),
                  f"{v:.2f}" if not tiny else f"{v:.4f}",
                  ha="center", va="bottom", fontsize=6.3 if not tiny else 5.8,
                  color=COL[series], fontweight="bold")
ax_d.set_xticks(xs)
ax_d.set_xticklabels(tier_lab)
frac_axis(ax_d)
ax_d.set_ylim(0, 1.12)
style(ax_d, "Precision at 100 bp (peak 3\u2032 base)", "Atlas subset used as truth")
r_tpm = get("precision", "real", "point", "tpm1")[0][i100]
n_tpm = get("precision", "N_genic", "point", "tpm1")[0][i100]
ax_d.set_title(f"On a 40 k high-confidence atlas the null collapses\n"
               f"to {n_tpm:.3f} while real PAS still reach {r_tpm:.2f}")
ax_d.legend(loc="upper right", frameon=False, handlelength=1.2, borderpad=0.2,
            labelspacing=0.28)
ax_d.text(0.99, 0.55, "Sparser truth set = lower precision for everyone,\n"
                      "but only the real calls survive it.",
          transform=ax_d.transAxes, fontsize=6.2, color=MUTED, va="bottom",
          ha="right")


# --- E. class of the nearest atlas match ------------------------------------
def cls_vals(series, dataset="lg_annotate", mode="point"):
    s = df[(df.panel == "match_class") & (df.series == series) &
           (df.dataset == dataset) & (df["mode"] == mode)]
    g = s.groupby("atlas_tier")["value"].mean()
    return np.array([g.get(f"class_{c}", 0.0) for c in CLASSES])


xs, w = np.arange(len(CLASSES)), 0.21
for j, series in enumerate(("real", "N_genome", "N_genic", "atlas")):
    ds = "atlas" if series == "atlas" else "lg_annotate"
    md = "background" if series == "atlas" else "point"
    v = cls_vals(series, ds, md)
    ax_e.bar(xs + (j - 1.5) * (w + 0.015), v, width=w, color=COL[series],
             label=LABEL[series], zorder=3, edgecolor=SURFACE, linewidth=0.5)
ax_e.set_xticks(xs)
ax_e.set_xticklabels(CLASSES, fontsize=7.6, fontweight="bold")
ax_e.set_ylim(0, 0.88)
ax_e.set_yticks(np.arange(0, 0.81, 0.2))
style(ax_e, "Fraction of matches at 100 bp",
      "PolyASite site class (col 10) of the nearest match")
te_real = cls_vals("real")[CLASSES.index("TE")]
te_bg = cls_vals("atlas", "atlas", "background")[CLASSES.index("TE")]
noise = sum(cls_vals("atlas", "atlas", "background")[CLASSES.index(c)]
            for c in ("IN", "DI", "UI"))
ng_ratio = cls_vals("N_genic") / cls_vals("atlas", "atlas", "background")
ax_e.set_title(f"Null matches track the atlas background ({ng_ratio.min():.1f}\u2013"
               f"{ng_ratio.max():.1f}\u00d7);\nreal matches favour terminal exons "
               f"{te_real/te_bg:.0f}\u00d7")
ax_e.legend(loc="upper left", bbox_to_anchor=(0.0, 0.82), frameon=False,
            handlelength=1.2, borderpad=0.2, labelspacing=0.28, ncol=1)
ax_e.text(0.99, 0.99,
          "TE terminal exon \u00b7 AL any last exon \u00b7 EX exonic\n"
          "IN intronic \u00b7 DI downstream \u00b7 UI upstream/intergenic\n"
          f"IN+DI+UI are {100*noise:.0f}% of all atlas entries.",
          transform=ax_e.transAxes, fontsize=6.2, color=MUTED, va="top",
          ha="right", linespacing=1.35)


# --- F. corrected recall ----------------------------------------------------
# label anchors chosen so the two near-overlapping bottom curves stay legible
lab_at = {"ref_full": (0, -0.035, "top", "left"),
          "ref_detected": (len(X) - 1, 0.02, "bottom", "right"),
          "ref_te": (len(X) - 1, 0.02, "bottom", "right"),
          "ref_tpm1": (1, -0.03, "top", "center")}
for key, (_, desc) in REFS.items():
    s = df[(df.panel == "recall") & (df.series == key) & (df["mode"] == "interval")]
    m = s.set_index("cutoff_bp")["value"].reindex(CUTOFFS).values
    assert int(s["n_query"].iloc[0]) == REF_N[key]   # legend n comes from the data
    ax_f.plot(X, m, color=COL[key], lw=1.8, marker="o", ms=6, mec=SURFACE,
              mew=0.9, zorder=3, label=f"{desc}, n={REF_N_LAB[key]}")
    xi, dy, va, ha = lab_at[key]
    ax_f.annotate(desc, (X[xi], m[xi] + dy), ha=ha, va=va,
                  color=COL[key], fontsize=6.6, fontweight="bold")
frac_axis(ax_f)
ax_f.set_ylim(-0.06, 1.06)
cutoff_axis(ax_f, "Match cutoff (bp to nearest called peak)")
style(ax_f, "Recall (fraction of reference sites found)",
      "Match cutoff (bp to nearest called peak)")


def rec_at(key, mode="interval", cut=100):
    return df[(df.panel == "recall") & (df.series == key) & (df["mode"] == mode)
              & (df.cutoff_bp == cut)]["value"].iloc[0]


rec = {k: rec_at(k) for k in REFS}
rec_pt = {k: rec_at(k, "point") for k in REFS}
ax_f.set_title(f"Recall at 100 bp rises from {rec['ref_full']:.3f} to "
               f"{rec['ref_tpm1']:.2f}\nonce the reference is made cohort-relevant")
ax_f.legend(loc="upper left", frameon=False, handlelength=1.6, borderpad=0.2,
            labelspacing=0.3)
ax_f.text(0.99, 0.74,
          f"Restriction: atlas sites inside the strand-matched gene\n"
          f"body of one of the {n_detected:,} genes this cohort assigned a\n"
          f"PAS to; then site class; then average TPM.",
          transform=ax_f.transAxes, fontsize=6.2, color=MUTED, va="top",
          ha="right")

fig.suptitle("A precision of 0.9986 against PolyASite 3.0 is mostly reference "
             "density \u2014 but real PAS still beat a matched null on every test",
             x=0.006, ha="left", fontsize=11, fontweight="bold", color=INK, y=0.995)
FOOTNOTE = (
    "Every panel above scores with strand-matched `bedtools closest -s -d`. The pipeline's own benchmark "
    "(ema/benchmark/metrics.py) does NOT require a strand match: it uses a symmetric \u00b1cutoff window, and reports precision "
    f"{PIPE_P100:.4f} and recall {PIPE_R100:.4f} at 100 bp on denominators {pj['n_predicted']:,} and {pj['n_reference']:,}. "
    "Its exact match counts are reproduced in the tsv (panel = pipeline_audit), so the strand requirement \u2014 not the "
    f"{REAL['lg_annotate']['dropped']} off-genome scaffold peaks \u2014 is why panel b reads {p_real_iv[i100]:.4f} and panel f "
    f"reads {rec['ref_full']:.3f} here. Panel f is interval-mode and so carries the same width inflation as panel b: scored on "
    "the peak 3\u2032 base its four recalls at 100 bp are "
    + " / ".join(f"{rec_pt[k]:.3f}" for k in REFS) + ".")
fig.text(0.006, 0.028, textwrap.fill(FOOTNOTE, 197),
         fontsize=6.4, color=MUTED, va="bottom", ha="left", linespacing=1.5)
fig.tight_layout(rect=[0, 0.107, 1, 0.957])
fig.savefig(OUTDIR / f"{NAME}.png", dpi=300, facecolor=SURFACE)
save_manuscript(fig, NAME, facecolor=SURFACE)
log(f"[out] {OUTDIR / (NAME + '.png')}")


# ----------------------------------------------------------------------------
# 9. console summary
# ----------------------------------------------------------------------------
log("\n================ HEADLINE ================")
log(f"reference density: {100*dm[0]:.0f}% / {100*dm[i100]:.0f}% / {100*dm[-1]:.0f}% "
    f"of genome within 50 / 100 / 1000 bp of a strand-matched atlas site")
for lbl, r, a, b in (("interval (whole peak)", p_real_iv, p_ng_iv, p_gc_iv),
                     ("point (peak 3' base)", p_real_pt, p_ng_pt, p_gc_pt)):
    for ci, c in enumerate(CUTOFFS):
        log(f"  {lbl:22s} @{c:>4}bp  real {r[ci]:.4f} | null_genome {a[ci]:.4f} "
            f"| null_genic {b[ci]:.4f} | gap {r[ci]-b[ci]:+.4f}")
log(f"tier tpm1 @100bp point: real {r_tpm:.4f} vs null_genic {n_tpm:.4f} "
    f"({r_tpm/n_tpm:.1f}x)")
log("match-class @100bp (point): " + "  ".join(
    f"{c}: real {cls_vals('real')[i]:.3f} / null {cls_vals('N_genic')[i]:.3f} "
    f"/ bg {cls_vals('atlas','atlas','background')[i]:.3f}"
    for i, c in enumerate(CLASSES)))
log("recall @100bp interval: " + "  ".join(f"{k} {rec[k]:.4f}" for k in REFS))
log("recall @100bp point   : " + "  ".join(f"{k} {rec_pt[k]:.4f}" for k in REFS))
log(f"PIPELINE (window, no strand, n=22633/18.43M) @100bp: precision "
    f"{PIPE_P100:.4f}, recall {PIPE_R100:.4f}  <- the published numbers; this "
    f"figure's strand-matched equivalents are {p_real_iv[i100]:.4f} and "
    f"{rec['ref_full']:.4f}")
b_iv = get("precision", "real", "interval", "full", "B1_cohort_full")[0]
b_pt = get("precision", "real", "point", "full", "B1_cohort_full")[0]
b_ng = get("precision", "N_genic", "point", "full", "B1_cohort_full")[0]
b_ng_iv = get("precision", "N_genic", "interval", "full", "B1_cohort_full")[0]
log(f"B1_cohort_full @100bp: real interval {b_iv[i100]:.4f} (null_genic "
    f"{b_ng_iv[i100]:.4f}) | real point {b_pt[i100]:.4f} (null_genic {b_ng[i100]:.4f})")
log("==========================================")
