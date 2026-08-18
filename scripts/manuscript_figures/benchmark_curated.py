#!/usr/bin/env python3
"""
benchmark_curated.py -- the defensible accuracy benchmark against the CURATED
PolyASite 2.0 atlas (569,005 clusters), replacing the 18.4M-entry dump that
null_control.py showed to be mostly reference density.

Design
------
REAL sets   5 grid arms + the cohort set (pasbed.bed each).  lg_ip_off is
            byte-identical to lg_annotate (asserted below) -> scored/plotted
            once, labelled 'lg_annotate = lg_ip_off'.

POINT MODE ONLY: each called peak is reduced to its 3'-most base by strand
            ('+': end-1; '-': start) = the inferred cleavage site.  Peaks on
            scaffold contigs absent from the reference genome file are dropped.

MATCHING    strand-matched `bedtools closest -s -d -t first`, point vs point.
            Cutoffs 10/25/50/100/200 bp.  NOTE: this is NOT the pipeline's own
            matcher (ema/benchmark/metrics.py uses a strand-agnostic window
            around the whole interval); numbers here are deliberately stricter.

NULLS       3 seeds of width-preserving `bedtools shuffle -chrom -noOverlapping
            -incl <merged gene bodies>` per arm (same approach as
            null_control.py), reduced to the 3' base and scored identically.

RECALL      two flavors, stated exactly:
              (a) vs the full curated atlas (569,005 representative sites);
              (b) vs the atlas restricted to strand-matched gene bodies of the
                  14,851 genes detected in this cohort (unique IDs in column 5
                  of runs/grid/lg_annotate/annotatedpas.bed, mapped to
                  data/references/gene_end.bed spans, bedtools intersect -s -u).
            Recall = fraction of reference sites with a called PAS 3'-base
            within the cutoff.

SANITY      the same precision/recall engine vs the protein-coding TES
            reference (Ensembl 99, 73,439 unique ends).

ROADMAP     publication bars: precision >= 0.70, recall >= 0.60, F1 >= 0.65 at
            the 100 bp (malleable) cutoff.  Pass/fail is evaluated per arm, per
            cutoff, per recall flavor, and written to the tsv + console.

Everything streams through bedtools/awk; no reference is loaded into pandas.
Cached under results/figures/manuscript/.cache_benchmark_curated/ -- delete to
force recompute.
"""

import hashlib
import os
import subprocess
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
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
REF_DIR = Path("/mnt/ssd1/Projects/PeakATail_wd/data/references")
PAS2 = REF_DIR / "atlases/polyasite2.GRCh38.96.rep_sites.bed6"   # PRIMARY ref
TES = REF_DIR / "atlases/tes.protein_coding.GRCh38.99.bed6"      # sanity ref
GENOME = REF_DIR / "chrom.sizes.nochr.filt"
GENE_END = REF_DIR / "gene_end.bed"
ANN = RUNS / "grid/lg_annotate/annotatedpas.bed"

OUTDIR = Path("/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript")
NAME = "benchmark_curated"
CACHE = OUTDIR / f".cache_{NAME}"
WORK = Path("/mnt/ssd2/claude-tmp/claude-1000/-mnt-ssd1-Projects-PeakATail-wd/"
            "93f9b511-1339-4cbc-9fdb-148fd4b2f08f/scratchpad/bench_curated")

CUTOFFS = [10, 25, 50, 100, 200]
N_SEEDS = 3
BAR_P, BAR_R, BAR_F1 = 0.70, 0.60, 0.65
BAR_CUT = 100

ENV = dict(os.environ, LC_ALL="C")  # tr_TR locale silently corrupts BED sorting

# validated palette -- one colour per ENTITY (arm), fixed order across panels
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

ARMS = ["lg_annotate", "lp_annotate", "si_annotate", "lg_ip_filter",
        "B1_cohort_full"]
SRC = {a: RUNS / f"grid/{a}/pasbed.bed" for a in ARMS[:4]}
SRC["B1_cohort_full"] = RUNS / "B1_cohort_full/pasbed.bed"
COL = dict(zip(ARMS, C))
ARM_LAB = {"lg_annotate": "lg_annotate = lg_ip_off",
           "lp_annotate": "lp_annotate",
           "si_annotate": "si_annotate",
           "lg_ip_filter": "lg_ip_filter",
           "B1_cohort_full": "B1_cohort_full"}


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


def log(*a):
    print(*a, flush=True)


def cached(key, cmd):
    f = CACHE / f"{key}.txt"
    if not f.exists():
        f.write_text(sh_out(cmd))
    return f.read_text()


for d in (WORK, OUTDIR, CACHE):
    d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------
# 0. sanity: lg_ip_off byte-identical to lg_annotate; chrom naming compatible
# ----------------------------------------------------------------------------
md5 = {t: hashlib.md5((RUNS / f"grid/{t}/pasbed.bed").read_bytes()).hexdigest()
       for t in ("lg_annotate", "lg_ip_off")}
assert md5["lg_annotate"] == md5["lg_ip_off"], "lg_ip_off differs from lg_annotate!"
log(f"[check] lg_ip_off byte-identical to lg_annotate (md5 {md5['lg_annotate']})")

genome_chroms = {l.split("\t")[0] for l in GENOME.read_text().strip().split("\n")}
N_PAS2 = int(sh_out(f"wc -l < {PAS2}").strip())
N_TES = int(sh_out(f"wc -l < {TES}").strip())
n_pas2_offgenome = int(sh_out(
    f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} !ok[$1]' {GENOME} {PAS2} | wc -l").strip())
log(f"[check] PolyASite 2.0 rep sites: {N_PAS2:,} ({n_pas2_offgenome:,} on "
    f"scaffolds absent from the genome file); TES ref: {N_TES:,}")

# ----------------------------------------------------------------------------
# 1. real sets: genome-filter, sort, reduce to 3'-most base
# ----------------------------------------------------------------------------
def make_point(iv, pt):
    """Collapse each peak to its 3'-most base by strand ('+': end-1; '-': start)."""
    sh("awk -F'\\t' 'BEGIN{OFS=\"\\t\"}"
       "{if($6==\"+\"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,$4,$5,$6}' "
       f"{iv} | sort -k1,1 -k2,2n > {pt}")


REAL = {}
for arm in ARMS:
    src = SRC[arm]
    raw_n = int(sh_out(f"wc -l < {src}").strip())
    iv, pt = WORK / f"{arm}.iv.bed", WORK / f"{arm}.pt.bed"
    if not pt.exists():
        sh(f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {GENOME} {src} "
           f"| sort -k1,1 -k2,2n > {iv}")
        make_point(iv, pt)
    n = int(sh_out(f"wc -l < {pt}").strip())
    REAL[arm] = dict(raw_n=raw_n, n=n, dropped=raw_n - n, iv=iv, pt=pt)
    log(f"[real] {arm}: n={n:,} points ({raw_n - n} scaffold peaks dropped)")

# ----------------------------------------------------------------------------
# 2. nulls: width-preserving shuffle inside merged annotated gene bodies
# ----------------------------------------------------------------------------
GENIC = WORK / "genebodies.merged.bed"
if not GENIC.exists():
    sh(f"sort -k1,1 -k2,2n {GENE_END} | bedtools merge -i - > {GENIC}")
GENIC_BP = int(sh_out("awk '{s+=$3-$2}END{print s}' " + str(GENIC)).strip())
log(f"[null] merged gene bodies span {GENIC_BP/1e9:.3f} Gb")

NULL_PT = {}   # (arm, seed) -> point bed
for arm in ARMS:
    for seed in range(1, N_SEEDS + 1):
        siv = WORK / f"{arm}.null.s{seed}.iv.bed"
        spt = WORK / f"{arm}.null.s{seed}.pt.bed"
        if not spt.exists():
            sh(f"bedtools shuffle -i {REAL[arm]['iv']} -g {GENOME} -chrom "
               f"-noOverlapping -maxTries 5000 -seed {seed} -incl {GENIC} "
               f"2>/dev/null | sort -k1,1 -k2,2n > {siv}")
            make_point(siv, spt)
        NULL_PT[(arm, seed)] = spt
    log(f"[null] {arm}: {N_SEEDS} seeds shuffled")

# ----------------------------------------------------------------------------
# 3. references: full atlas, detected-gene-restricted atlas, TES
# ----------------------------------------------------------------------------
genes_txt = WORK / "detected_genes.txt"
det_bed = WORK / "detected_gene_bodies.bed"
ref_det = WORK / "pas2.in_detected_genes.bed"
if not ref_det.exists():
    sh(f"cut -f5 {ANN} | sort -u > {genes_txt}")
    sh(f"awk -F'\\t' 'NR==FNR{{g[$1]=1;next}} g[$4]' {genes_txt} {GENE_END} "
       f"| sort -k1,1 -k2,2n > {det_bed}")
    sh(f"bedtools intersect -a {PAS2} -b {det_bed} -s -u > {ref_det}")
N_DET_GENES = int(sh_out(f"wc -l < {genes_txt}").strip())
N_REF_DET = int(sh_out(f"wc -l < {ref_det}").strip())
assert N_DET_GENES == 14851, f"expected 14,851 detected genes, got {N_DET_GENES}"
log(f"[ref] detected genes {N_DET_GENES:,}; atlas sites in their strand-matched "
    f"gene bodies: {N_REF_DET:,} ({100*N_REF_DET/N_PAS2:.1f}% of the atlas)")

REFS = {"atlas_full": (PAS2, "full curated atlas"),
        "atlas_detected": (ref_det, "atlas in detected-gene bodies"),
        "tes": (TES, "protein-coding TES (sanity)")}

# ----------------------------------------------------------------------------
# 4. engine: strand-matched nearest-site histogram at the cutoffs
# ----------------------------------------------------------------------------
AWK_HIST = (
    "awk -F'\\t' -v CUT=\"" + ",".join(map(str, CUTOFFS)) + "\" '"
    "BEGIN{nc=split(CUT,cc,\",\")}"
    "{d=$NF; n++; if(d>=0){for(i=1;i<=nc;i++) if(d<=cc[i]) h[i]++}}"
    "END{printf \"%d\", n; for(i=1;i<=nc;i++) printf \"\\t%d\", h[i]+0; printf \"\\n\"}'"
)


def closest_hist(a, b, key):
    out = cached(key, f"bedtools closest -s -d -t first -a {a} -b {b} "
                      f"2>/dev/null | {AWK_HIST}").strip().split("\t")
    return int(out[0]), dict(zip(CUTOFFS, (int(x) for x in out[1:])))


rows = []


def record(**kw):
    rows.append(kw)


# --- precision: real + null, vs PolyASite 2.0 and TES ------------------------
P = {}       # (arm, ref) -> {cutoff: value}   real precision
NULLP = {}   # (arm, seed) -> {cutoff: value}  null precision vs atlas_full
for arm in ARMS:
    pt = REAL[arm]["pt"]
    for refkey in ("atlas_full", "tes"):
        n, hits = closest_hist(pt, REFS[refkey][0], f"prec__{arm}__{refkey}")
        P[(arm, refkey)] = {c: hits[c] / n for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="precision", arm=arm, series="real", reference=refkey,
                   cutoff_bp=c, replicate=0, n_query=n, n_matched=hits[c],
                   value=hits[c] / n)
    for seed in range(1, N_SEEDS + 1):
        n, hits = closest_hist(NULL_PT[(arm, seed)], PAS2,
                               f"prec__{arm}.null.s{seed}__atlas_full")
        NULLP[(arm, seed)] = {c: hits[c] / n for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="precision", arm=arm, series="null_genic",
                   reference="atlas_full", cutoff_bp=c, replicate=seed,
                   n_query=n, n_matched=hits[c], value=hits[c] / n)
    log(f"[prec] {arm}: @100bp real {P[(arm,'atlas_full')][100]:.4f} "
        f"(TES {P[(arm,'tes')][100]:.4f}), null "
        + "/".join(f"{NULLP[(arm,s)][100]:.4f}" for s in range(1, N_SEEDS + 1)))

# --- recall: both flavors + TES sanity ---------------------------------------
R = {}   # (arm, ref) -> {cutoff: value}
for arm in ARMS:
    pt = REAL[arm]["pt"]
    for refkey, (refp, _) in REFS.items():
        n, hits = closest_hist(refp, pt, f"rec__{refkey}__{arm}")
        R[(arm, refkey)] = {c: hits[c] / n for c in CUTOFFS}
        for c in CUTOFFS:
            record(panel="recall", arm=arm, series="real", reference=refkey,
                   cutoff_bp=c, replicate=0, n_query=n, n_matched=hits[c],
                   value=hits[c] / n)
    log(f"[recall] {arm}: @100bp full {R[(arm,'atlas_full')][100]:.4f} | "
        f"detected {R[(arm,'atlas_detected')][100]:.4f} | "
        f"TES {R[(arm,'tes')][100]:.4f}")

# --- F1 + roadmap pass/fail --------------------------------------------------
def f1(p, r):
    return 0.0 if p + r == 0 else 2 * p * r / (p + r)


F1 = {}
for arm in ARMS:
    for flav in ("atlas_full", "atlas_detected"):
        for c in CUTOFFS:
            v = f1(P[(arm, "atlas_full")][c], R[(arm, flav)][c])
            F1[(arm, flav, c)] = v
            record(panel="f1", arm=arm, series="real", reference=flav,
                   cutoff_bp=c, replicate=0, n_query=np.nan, n_matched=np.nan,
                   value=v)
            for metric, val, bar in (("precision", P[(arm, "atlas_full")][c], BAR_P),
                                     ("recall", R[(arm, flav)][c], BAR_R),
                                     ("f1", v, BAR_F1)):
                record(panel="roadmap", arm=arm, series=metric, reference=flav,
                       cutoff_bp=c, replicate=0, n_query=np.nan,
                       n_matched=np.nan, value=val, bar=bar,
                       passes=int(val >= bar))

for arm in ARMS:
    record(panel="meta", arm=arm, series="n_points", reference="",
           cutoff_bp=np.nan, replicate=0, n_query=REAL[arm]["raw_n"],
           n_matched=REAL[arm]["n"], value=REAL[arm]["dropped"])

df = pd.DataFrame(rows)
df.to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False, float_format="%.6f")
log(f"[out] {OUTDIR / (NAME + '.tsv')}  ({len(df)} rows)")

# ----------------------------------------------------------------------------
# 5. figure
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
fig = plt.figure(figsize=(11.4, 8.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.86], hspace=0.44, wspace=0.24,
                      left=0.06, right=0.985, top=0.9, bottom=0.16)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, :])
X = np.arange(len(CUTOFFS))
i100, i50 = CUTOFFS.index(100), CUTOFFS.index(50)


def style(ax, ylab, xlab):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_xticks(X)
    ax.set_xticklabels([str(c) for c in CUTOFFS])


def spread(vals, minsep):
    """Nudge label y-positions apart (ascending pass) so they never collide."""
    vals = np.asarray(vals, float)
    order = np.argsort(vals)
    out = vals.copy()
    for i in range(1, len(order)):
        lo, hi = order[i - 1], order[i]
        if out[hi] - out[lo] < minsep:
            out[hi] = out[lo] + minsep
    return out


# --- A. precision vs cutoff, per arm, with pooled null band ------------------
null_mat = np.array([[NULLP[(a, s)][c] for c in CUTOFFS]
                     for a in ARMS for s in range(1, N_SEEDS + 1)])
ax_a.fill_between(X, null_mat.min(0), null_mat.max(0), color=MUTED, alpha=0.22,
                  lw=0, zorder=1)
ax_a.plot(X, null_mat.mean(0), color=MUTED, lw=1.4, ls="--", zorder=2)
ax_a.annotate("null: width-preserving shuffle\nin gene bodies "
              f"({len(ARMS)} arms \u00d7 {N_SEEDS} seeds)",
              (X[1], null_mat[:, 1].max() + 0.03), color=MUTED, fontsize=6.4,
              ha="left", va="bottom", fontweight="bold")
ax_a.axhline(BAR_P, color=INK, lw=0.8, ls=":", zorder=2)
ax_a.text(0.02, BAR_P + 0.012, f"ROADMAP precision bar {BAR_P:.2f}",
          color=INK, fontsize=6.2, ha="left", va="bottom",
          transform=ax_a.get_yaxis_transform())
pv = {a: [P[(a, "atlas_full")][c] for c in CUTOFFS] for a in ARMS}
ends = spread([pv[a][-1] for a in ARMS], 0.052)
for a, ye in zip(ARMS, ends):
    ax_a.plot(X, pv[a], color=COL[a], lw=1.8, marker="o", ms=5.5, mec=SURFACE,
              mew=0.9, zorder=3, label=ARM_LAB[a])
    ax_a.annotate(ARM_LAB[a], (X[-1] + 0.13, ye), color=COL[a], fontsize=6.5,
                  ha="left", va="center", fontweight="bold",
                  annotation_clip=False)
ax_a.set_xlim(-0.35, len(X) - 1 + 1.85)
ax_a.set_ylim(0, 1.06)
ax_a.set_yticks(np.arange(0, 1.01, 0.2))
style(ax_a, "Precision (fraction of called 3\u2032 bases matched)",
      "Match cutoff (bp to nearest curated atlas site)")
p100 = [P[(a, "atlas_full")][100] for a in ARMS]
null100 = null_mat[:, i100]
ax_a.set_title(f"Precision on the curated atlas is {min(p100):.2f}\u2013"
               f"{max(p100):.2f} at 100 bp;\nthe matched null sits at "
               f"{null100.mean():.2f}")
ax_a.legend(loc="lower right", frameon=False, handlelength=1.5, borderpad=0.2,
            labelspacing=0.28, bbox_to_anchor=(1.0, 0.02))

# --- B. recall (both flavors) vs cutoff --------------------------------------
rb = {a: [R[(a, "atlas_detected")][c] for c in CUTOFFS] for a in ARMS}
ra = {a: [R[(a, "atlas_full")][c] for c in CUTOFFS] for a in ARMS}
ends = spread([rb[a][-1] for a in ARMS], 0.038)
for a, ye in zip(ARMS, ends):
    ax_b.plot(X, ra[a], color=COL[a], lw=1.2, ls="--", alpha=0.75, zorder=2)
    ax_b.plot(X, rb[a], color=COL[a], lw=1.8, marker="o", ms=5.5, mec=SURFACE,
              mew=0.9, zorder=3)
    ax_b.annotate(ARM_LAB[a], (X[-1] + 0.13, ye), color=COL[a], fontsize=6.5,
                  ha="left", va="center", fontweight="bold",
                  annotation_clip=False)
ax_b.axhline(BAR_R, color=INK, lw=0.8, ls=":", zorder=2)
ax_b.text(0.02, BAR_R + 0.012, f"ROADMAP recall bar {BAR_R:.2f}", color=INK,
          fontsize=6.2, ha="left", va="bottom",
          transform=ax_b.get_yaxis_transform())
rb100 = [R[(a, "atlas_detected")][100] for a in ARMS]
ra100 = [R[(a, "atlas_full")][100] for a in ARMS]
ymax = max(0.72, max(max(rb[a][-1] for a in ARMS), ends.max()) * 1.12 + 0.04)
ax_b.set_xlim(-0.35, len(X) - 1 + 1.85)
ax_b.set_ylim(0, ymax)
style(ax_b, "Recall (fraction of reference sites found)",
      "Match cutoff (bp to nearest called 3\u2032 base)")
if max(rb100) < BAR_R:
    ax_b.set_title(f"Recall at 100 bp reaches only {max(rb100):.2f} even on the\n"
                   f"detected-gene reference \u2014 every arm is far below the bar")
else:
    ax_b.set_title(f"Recall at 100 bp: {min(rb100):.2f}\u2013{max(rb100):.2f} on "
                   f"the detected-gene reference")
ax_b.legend(handles=[
    Line2D([], [], color=MUTED, lw=1.8, marker="o", ms=5, mec=SURFACE,
           label=f"flavor b: atlas in detected-gene bodies (n={N_REF_DET:,})"),
    Line2D([], [], color=MUTED, lw=1.2, ls="--",
           label=f"flavor a: full curated atlas (n={N_PAS2:,})")],
    loc="upper left", frameon=False, handlelength=1.9, borderpad=0.2,
    labelspacing=0.32)

# --- C. table: P / R / F1 at 50 and 100 bp vs the roadmap bars ---------------
ax_c.axis("off")
ax_c.set_title("Scorecard at 50 and 100 bp against the ROADMAP publication bars "
               f"(P \u2265 {BAR_P:.2f}, R \u2265 {BAR_R:.2f}, F1 \u2265 {BAR_F1:.2f} "
               "at the 100 bp malleable cutoff)")
GOOD, BAD = "#009E73", "#D55E00"
col_hdr = ["P", "R (a)", "R (b)", "F1 (b)"] * 2 + ["TES P"]
col_x = np.concatenate([np.linspace(0.20, 0.44, 4),
                        np.linspace(0.52, 0.80, 4), [0.865]])
name_x, verdict_x = 0.002, 1.0
y_group, y_hdr = 0.985, 0.895
row_y = np.linspace(0.76, 0.18, len(ARMS))
row_h = (row_y[0] - row_y[1])

ax_c.text((col_x[0] + col_x[3]) / 2, y_group, "at 50 bp", ha="center",
          va="center", fontsize=7.4, color=INK, fontweight="bold")
ax_c.text((col_x[4] + col_x[7]) / 2, y_group, "at 100 bp (roadmap cutoff)",
          ha="center", va="center", fontsize=7.4, color=INK, fontweight="bold")
ax_c.text(verdict_x, y_group, "bars @100", ha="right", va="center",
          fontsize=7.4, color=INK, fontweight="bold")
for xh, h in zip(col_x, col_hdr):
    ax_c.text(xh, y_hdr, h, ha="center", va="center", fontsize=6.8,
              color=MUTED, fontweight="bold")
ax_c.text(verdict_x, y_hdr, "(flavor b)", ha="right", va="center",
          fontsize=6.4, color=MUTED)
ax_c.plot([0.48, 0.48], [0.10, 1.0], color=GRID, lw=0.8, clip_on=False)
ax_c.plot([0.835, 0.835], [0.10, 1.0], color=GRID, lw=0.8, clip_on=False)
ax_c.plot([0.92, 0.92], [0.10, 1.0], color=GRID, lw=0.8, clip_on=False)
# text/plot above use data coords; pin them to the unit square so autoscale
# from the separator lines cannot clip the table
ax_c.set_xlim(0, 1)
ax_c.set_ylim(0, 1)

npass = {"P": 0, "R": 0, "F1": 0}
for i, (a, y) in enumerate(zip(ARMS, row_y)):
    if i % 2 == 0:
        ax_c.add_patch(Rectangle((0, y - row_h * 0.42), 1.0, row_h * 0.84,
                                 transform=ax_c.transAxes, fc="#F0F4F5", ec="none",
                                 zorder=0, clip_on=False))
    ax_c.text(name_x, y + 0.028, ARM_LAB[a], ha="left", va="center",
              fontsize=6.9, color=COL[a], fontweight="bold")
    ax_c.text(name_x, y - 0.058, f"n = {REAL[a]['n']:,} called 3\u2032 bases",
              ha="left", va="center", fontsize=6.0, color=MUTED)
    vals = [(P[(a, "atlas_full")][50], None),
            (R[(a, "atlas_full")][50], None),
            (R[(a, "atlas_detected")][50], None),
            (F1[(a, "atlas_detected", 50)], None),
            (P[(a, "atlas_full")][100], BAR_P),
            (R[(a, "atlas_full")][100], BAR_R),
            (R[(a, "atlas_detected")][100], BAR_R),
            (F1[(a, "atlas_detected", 100)], BAR_F1),
            (P[(a, "tes")][100], None)]
    for xh, (v, bar) in zip(col_x, vals):
        c = INK if bar is None else (GOOD if v >= bar else BAD)
        ax_c.text(xh, y, f"{v:.3f}", ha="center", va="center", fontsize=6.9,
                  color=c, fontweight="bold" if bar is not None else "normal")
    okP = P[(a, "atlas_full")][100] >= BAR_P
    okR = R[(a, "atlas_detected")][100] >= BAR_R
    okF = F1[(a, "atlas_detected", 100)] >= BAR_F1
    npass["P"] += okP; npass["R"] += okR; npass["F1"] += okF
    verdict = (f"P {'\u2713' if okP else '\u2717'}  "
               f"R {'\u2713' if okR else '\u2717'}  "
               f"F1 {'\u2713' if okF else '\u2717'}")
    ax_c.text(verdict_x, y, verdict, ha="right", va="center", fontsize=6.9,
              color=GOOD if (okP and okR and okF) else BAD, fontweight="bold")
verdict_line = (
    f"Verdict at the 100 bp roadmap cutoff: {npass['P']}/{len(ARMS)} arms "
    f"clear the precision bar; {npass['R']}/{len(ARMS)} clear recall and "
    f"{npass['F1']}/{len(ARMS)} clear F1 even on the favourable detected-gene "
    f"reference (flavor b). Recall flavor (a), the full curated atlas, is "
    f"lower still \u2014 an arm calling "
    f"~{int(round(np.mean([REAL[a]['n'] for a in ARMS]) / 1000))}k sites "
    f"cannot recover a {N_PAS2 / 1000:.0f}k-site atlas.")
ax_c.text(0, 0.065, textwrap.fill(verdict_line, 185), ha="left", va="top",
          fontsize=6.8, color=INK, linespacing=1.45)

fold = np.mean(p100) / null100.mean()
if npass["P"] + npass["R"] + npass["F1"] == 0:
    headline = (f"No arm clears any ROADMAP bar on the curated PolyASite 2.0 "
                f"atlas: precision {min(p100):.2f}\u2013{max(p100):.2f} at 100 bp "
                f"(\u2248{fold:.0f}\u00d7 null), recall \u2264 {max(rb100):.2f}")
else:
    headline = (f"Curated PolyASite 2.0 benchmark: precision "
                f"{min(p100):.2f}\u2013{max(p100):.2f} at 100 bp "
                f"(\u2248{fold:.0f}\u00d7 null), recall {min(rb100):.2f}\u2013"
                f"{max(rb100):.2f} on the detected-gene reference \u2014 "
                f"{npass['P']}/P {npass['R']}/R {npass['F1']}/F1 arms pass")
fig.suptitle(headline, x=0.006, ha="left", fontsize=10.5, fontweight="bold",
             color=INK, y=0.99)

FOOT = (
    "Matching: strand-matched `bedtools closest -s -d` between each called peak's 3\u2032-most base ('+': end\u22121; '\u2212': start) and the "
    f"PolyASite 2.0 representative-site point BED ({N_PAS2:,} curated clusters, GRCh38.96) \u2014 point mode only; this is stricter than the "
    "pipeline's own strand-agnostic whole-interval window matcher, and none of these numbers is comparable to the 0.9986 previously quoted "
    "against the 18.4M-entry dump (see null_control). "
    f"Null: {N_SEEDS} seeds of width-preserving `bedtools shuffle -chrom -noOverlapping` inside merged annotated gene bodies "
    f"({GENIC_BP/1e9:.2f} Gb), reduced to the 3\u2032 base and scored identically; band = min\u2013max over all arms \u00d7 seeds. "
    f"Recall flavor (a) denominator is all {N_PAS2:,} curated sites and includes {n_pas2_offgenome:,} sites on scaffolds no call can reach; "
    f"flavor (b) restricts to sites inside strand-matched gene bodies of the {N_DET_GENES:,} genes detected in this cohort "
    f"(n={N_REF_DET:,}). F1 pairs precision on the full curated atlas with the stated recall flavor. TES P = precision vs the "
    f"{N_TES:,} protein-coding transcript ends (Ensembl 99), a coarser sanity reference. lg_ip_off/pasbed.bed is byte-identical to "
    f"lg_annotate (md5 {md5['lg_annotate'][:8]}\u2026), scored once. Scaffold peaks dropped per arm: "
    + ", ".join(f"{a} {REAL[a]['dropped']}" for a in ARMS) + ".")
fig.text(0.006, 0.012, textwrap.fill(FOOT, 208), fontsize=6.15, color=MUTED,
         va="bottom", ha="left", linespacing=1.45)

fig.savefig(OUTDIR / f"{NAME}.png", dpi=300, facecolor=SURFACE)
save_manuscript(fig, NAME, facecolor=SURFACE)
log(f"[out] {OUTDIR / (NAME + '.png')}")

# ----------------------------------------------------------------------------
# 6. console summary (plain pass/fail statements)
# ----------------------------------------------------------------------------
log("\n================ HEADLINE ================")
for a in ARMS:
    log(f"{ARM_LAB[a]} (n={REAL[a]['n']:,}):")
    for c in CUTOFFS:
        log(f"  @{c:>3}bp  P={P[(a,'atlas_full')][c]:.4f}  "
            f"Ra={R[(a,'atlas_full')][c]:.4f}  Rb={R[(a,'atlas_detected')][c]:.4f}  "
            f"F1a={F1[(a,'atlas_full',c)]:.4f}  F1b={F1[(a,'atlas_detected',c)]:.4f}  "
            f"TES_P={P[(a,'tes')][c]:.4f}  null_P={np.mean([NULLP[(a,s)][c] for s in range(1,N_SEEDS+1)]):.4f}")
    okP = P[(a, "atlas_full")][BAR_CUT] >= BAR_P
    okRb = R[(a, "atlas_detected")][BAR_CUT] >= BAR_R
    okFb = F1[(a, "atlas_detected", BAR_CUT)] >= BAR_F1
    okRa = R[(a, "atlas_full")][BAR_CUT] >= BAR_R
    okFa = F1[(a, "atlas_full", BAR_CUT)] >= BAR_F1
    p_cuts = [c for c in CUTOFFS if P[(a, "atlas_full")][c] >= BAR_P]
    r_cuts_b = [c for c in CUTOFFS if R[(a, "atlas_detected")][c] >= BAR_R]
    log(f"  ROADMAP @100bp on PolyASite2 curated: precision "
        f"{'PASS' if okP else 'FAIL'}"
        + (f" (passes at cutoffs {p_cuts})" if p_cuts else " (no cutoff passes)")
        + f"; recall flavor b {'PASS' if okRb else 'FAIL'}"
        + (f" (passes at {r_cuts_b})" if r_cuts_b else " (no cutoff passes)")
        + f"; recall flavor a {'PASS' if okRa else 'FAIL'}; "
        f"F1 b {'PASS' if okFb else 'FAIL'}; F1 a {'PASS' if okFa else 'FAIL'}")
log("==========================================")
