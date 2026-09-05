#!/usr/bin/env python3
"""
spermatogenesis_control.py -- KNOWN-BIOLOGY POSITIVE CONTROL for PeakATail's
APA quantification, on the GSE104556 adult-mouse-testis replicates
(Mouse1_scRNAseq = SRR6129050, Mouse2_scRNAseq = SRR6129051).

The claim under test
--------------------
3'UTRs shorten progressively along spermatogenesis
(spermatogonia > spermatocyte > round spermatid > elongating/condensing).
This is the strongest documented APA gradient in any tissue (Li 2016;
Berkovits & Mayr 2015; MacDonald 2019).  If PeakATail cannot recover it, its
APA quantification is not trustworthy.

What had to be fixed before the control could be run
----------------------------------------------------
BLOCKER 1 (fatal, in the tool).  ``ema/matrixfilter.py::make_dataframe``
returns the FULL-height count matrix (one row per pas_id, 1..max_id) but
returns ``pas_ids = np.sort(np.unique(coo.row + 1))``, i.e. only the rows that
carry at least one non-zero.  ``ema/annotate/annotate.py`` then builds
``pas_id_to_row = {pid: i for i, pid in enumerate(pas_ids)}`` and slices
``sparse_matrix[keep_rows]``.  Whenever any PAS row is all-zero the two index
spaces diverge, so PAS *p* is written out carrying the count vector of the PAS
``p - (# all-zero ids below p)``.  On mouse1 1,540 of 98,574 rows are empty, so
the drift runs 0 -> 1,539 across the genome; on mouse2 it runs 0 -> 1,208.
Evidence: per-PAS matrix totals vs. samtools read depth over the SAME interval
(200 stratified PAS, mouse1) -> Spearman 0.058 as shipped, 0.730 after the
repair.  Concretely: the peak over Prm1's 3' end (16:10795807-10796877,
4.30 M BAM reads) carries 642 counts in the shipped matrix, while
16:45157472-45157558 (55 BAM reads) carries 726,873.
REPAIR used here: read ``filterdmatrix.mtx`` directly -- its MatrixMarket row
index IS the pas_id -- and re-key rows by pas_id.

BLOCKER 2.  ``07_clustering/*/clusters.h5ad`` stores the TF-IDF-transformed
matrix in ``.X`` (values are non-integer, min ~0.2), not counts.  The
documented ``peakatail switch length -i clusters.h5ad`` therefore computes PDUI on
TF-IDF weights; proximal and distal PAS of the same gene carry different IDF
factors, so the ratio is biased.  Raw counts survive in
``06_preprocessing/*/preprocessed.h5ad`` (identical obs/var order, verified).

BLOCKER 3.  ``peakatail switch length --strategy classic --isoform-agg per_gene``
is strand-blind.  The synthesised per-gene map gives every PAS rank=1, so the
strategy's ``sort_values("rank")`` degenerates to pas_id order (= ascending
genomic coordinate) and proximal/distal are swapped for every MINUS-strand
gene.  Verified by unit test: for a 2-PAS minus-strand gene the package's own
reference ``ema.quantification.pdui.calculate_pdui`` returns PDUI 0.9 where
``ClassicPDUIStrategy.compute(..., aggregation="per_gene")`` returns 0.1.

BLOCKER 4 (not a bug, but it decides the analysis).  PeakATail's CB filter
(min_read=1500 on PAS reads) keeps 9,564 / 9,074 barcodes where STARsolo calls
1,294 / 1,364 cells -- a strict superset, so ~86% of the "cells" are ambient
droplets (median 963-993 UMI vs 12.9-15.7 k in real cells).  All analyses here
are restricted to the STARsolo-called cells.

Design
------
STAGES are assigned twice and independently:
  (a) PAS-derived  -- the repaired PeakATail PAS matrix summed to gene level,
      TF-IDF + LSI + leiden, marker panels scored with scanpy.tl.score_genes,
      each cluster given the arg-max panel (germ panels z-scored across
      clusters; a somatic panel wins only if its ABSOLUTE mean marker
      expression exceeds 0.35, because a pure z-argmax forces every cluster
      into some panel).  This is the "PAS-only" labelling the control asks for.
  (b) GEX-derived  -- the orthogonal STARsolo gene matrix from the same BAM,
      standard scanpy pipeline, same panels.  These labels are NOT derived
      from any PeakATail output.
Concordance (a) vs (b): 93.9% (mouse1), 93.8% (mouse2).

3'UTR LENGTH INDEX.  For every gene with >=2 PAS in the analysed subset, PAS
are ordered in transcription order (strand-aware; a PAS's position is the
3'-most base of its interval).  With rank r_i in 0..n-1,
    relative position   p_i = r_i / (n-1)                    (0 = proximal)
    WDI_g(cell)  = sum_i x_ig * p_i / sum_i x_ig
    WUL_g(cell)  = sum_i x_ig * d_i / sum_i x_ig    (d_i = |pos_i - pos_prox|, bp)
A cell's index is the unweighted mean of WDI_g over genes with >=5 UMI in that
cell (>=20 such genes required).  WDI falls when usage moves proximal, i.e.
when 3'UTRs SHORTEN.  Cluster/stage pseudobulk uses the same formulae on
summed counts, on a gene set fixed across stages.

PAS SUBSETS (the sensitivity analysis the control demands):
    all       every gene-assigned PAS
    atlas     within 100 bp (strand-matched) of polyasite2.GRCm38.96 rep sites
    noatlas   the complement of `atlas`
    utr3      inside the PAS's OWN gene's annotated 3'UTR (Ensembl 102) + 1 kb

Outputs
-------
  results/figures/manuscript/spermatogenesis_control.png        (300 dpi)
  manuscript/figures/spermatogenesis_control.{png,pdf}          (fonttype 42)
  results/figures/manuscript/spermatogenesis_control*.tsv       (every plotted value)
Analysis steps that produce the cached intermediates live in
  scripts/manuscript_figures/spermatogenesis_control_steps/
"""

import os
os.environ.setdefault("LC_ALL", "C")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "8"

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# --- manuscript deliverables: PNG + vector PDF into manuscript/figures/ ---
# (helper copied verbatim from scripts/manuscript_figures/null_control.py)
import matplotlib as _mpl
_mpl.rcParams["pdf.fonttype"] = 42   # TrueType, not Type 3 (journal requirement)
_mpl.rcParams["ps.fonttype"] = 42
FIGDIR = "/mnt/ssd1/Projects/PeakATail_wd/manuscript/figures"
os.makedirs(FIGDIR, exist_ok=True)
def save_manuscript(fig, name, **kw):
    for ext in ("png", "pdf"):
        p = os.path.join(FIGDIR, f"{name}.{ext}")
        fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), **kw)
        print("wrote", p)

SCRATCH = Path("/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl")
OUTDIR  = Path("/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript")
OUTDIR.mkdir(parents=True, exist_ok=True)
NAME = "spermatogenesis_control"

# validated palette -- one colour per ENTITY, fixed across every panel
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
STAGES = ["Spermatogonia", "Spermatocyte", "RoundSpermatid", "Elongating"]
SHORT  = {"Spermatogonia": "SPG", "Spermatocyte": "SPC",
          "RoundSpermatid": "RS", "Elongating": "ES"}
SCOL   = {"Spermatogonia": C[0], "Spermatocyte": C[2],
          "RoundSpermatid": C[3], "Elongating": C[1]}
SUBCOL = {"all": MUTED, "atlas": C[0], "noatlas": C[4], "utr3": C[2]}
SUBLAB = {"all": "all called PAS", "atlas": "atlas-supported (<=100 bp)",
          "noatlas": "atlas-UNsupported", "utr3": "inside own 3'UTR (+1 kb)"}
PANELS = ["Spermatogonia", "Spermatocyte", "RoundSpermatid", "Elongating",
          "Sertoli", "Leydig", "Immune"]
MICE = ["mouse1", "mouse2"]


def style(ax, ylab=None, xlab=None, title=None):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    if ylab:  ax.set_ylabel(ylab, color=INK, fontsize=9)
    if xlab:  ax.set_xlabel(xlab, color=INK, fontsize=9)
    if title: ax.set_title(title, color=INK, fontsize=10, loc="left", pad=6)


# ---------------------------------------------------------------- load
labels = {m: pd.read_csv(SCRATCH / f"{m}_cell_labels.tsv", sep="\t") for m in MICE}
trend  = pd.read_csv(SCRATCH / "trend_stats_v2.tsv", sep="\t")
percell = pd.read_csv(SCRATCH / "percell_v2.tsv.gz", sep="\t")
pbfix  = pd.read_csv(SCRATCH / "pseudobulk_fixedgene.tsv", sep="\t")
shortg = {s: pd.read_csv(SCRATCH / f"gene_shortening_{s}.tsv", sep="\t", index_col=0)
          for s in ("utr3", "atlas")}
mstage = {m: pd.read_csv(SCRATCH / f"{m}_PAS_cluster_stage.tsv", sep="\t", index_col=0)
          for m in MICE}
repval = pd.read_csv(SCRATCH / "matrix_repair_validation.tsv", sep="\t")
ipdiag = pd.read_csv(SCRATCH / "internal_priming_diag.tsv", sep="\t")

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(16.5, 14.4), facecolor="white")
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 1.0, 1.15],
              hspace=0.42, wspace=0.28,
              left=0.055, right=0.985, top=0.925, bottom=0.055)

fig.suptitle(
    "Spermatogenesis positive control: PeakATail recovers progressive 3'UTR shortening "
    "only for PAS it places in real 3'UTRs",
    x=0.055, y=0.975, ha="left", color=INK, fontsize=15, fontweight="bold")
fig.text(0.055, 0.949,
         "GSE104556 adult mouse testis, two mice analysed independently  |  "
         "PeakATail 0.2.0 count matrix re-keyed by pas_id (shipped matrix rows are mis-indexed; see caveats)  |  "
         "cells restricted to the 1,294 / 1,364 STARsolo-called cells",
         ha="left", color=MUTED, fontsize=9.5)

# ---- (a) UMAP by assigned stage -------------------------------------------
umap_rows = []
for k, m in enumerate(MICE):
    ax = fig.add_subplot(gs[0, k])
    d = labels[m]
    for st in STAGES:
        s = d[d["stage_pas"] == st]
        if len(s) == 0:
            continue
        ax.scatter(s["umap1_pas"], s["umap2_pas"], s=5, lw=0, alpha=0.85,
                   color=SCOL[st], label=f"{SHORT[st]}  n={len(s)}")
    style(ax, "UMAP-2", "UMAP-1",
          f"(a{k+1}) {m}: stage from the PAS matrix alone")
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.legend(frameon=False, fontsize=7.5, loc="best", labelcolor=INK,
              handletextpad=0.2, borderpad=0.2)
    agree = (d["stage_pas"] == d["stage_gex"]).mean()
    ax.text(0.0, -0.125,
            f"agreement with orthogonal STARsolo-GEX labels: {agree:.1%}\n"
            "PAS-only clustering merges the rare spermatogonia into SPC,\n"
            "so panels b use the 4-stage STARsolo-GEX labels",
            transform=ax.transAxes, fontsize=7.0, color=MUTED, va="top")
    u = d[["cb", "stage_pas", "stage_gex", "umap1_pas", "umap2_pas", "umi_pas"]].copy()
    u["mouse"] = m
    umap_rows.append(u)
pd.concat(umap_rows).to_csv(OUTDIR / f"{NAME}_panelA_umap.tsv", sep="\t", index=False)

# ---- (d) marker-score matrix (audit) --------------------------------------
ax = fig.add_subplot(gs[0, 2])
blocks, ylabs, rowsrc = [], [], []
for m in MICE:
    t = mstage[m].copy()
    t.index = [f"{m[-1]}·{i}" for i in t.index]
    z = t[[f"z_s_{p}" for p in PANELS]]
    z.columns = PANELS
    order = t["stage"].map({s: i for i, s in enumerate(STAGES)}).fillna(9).values
    idx = np.argsort(order, kind="stable")
    z = z.iloc[idx]; t = t.iloc[idx]
    blocks.append(z)
    ylabs += [f"{i}  {SHORT.get(s, s[:3])}  n={n}"
              for i, s, n in zip(z.index, t["stage"], t["n_cells"])]
    r = z.copy(); r["mouse"] = m; r["assigned_stage"] = t["stage"].values
    r["n_cells"] = t["n_cells"].values
    rowsrc.append(r)
Z = pd.concat(blocks)
im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
ax.set_xticks(range(len(PANELS)))
ax.set_xticklabels([SHORT.get(p, p) for p in PANELS], rotation=40, ha="right", fontsize=7.6)
ax.set_yticks(range(len(ylabs)))
ax.set_yticklabels(ylabs, fontsize=6.2)
ax.axhline(len(blocks[0]) - 0.5, color=INK, lw=1.2)
style(ax, None, None, "(d) marker-panel score per PAS-derived cluster (audit)")
ax.grid(False)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cb.set_label("panel score, z across clusters", fontsize=7.5, color=INK)
cb.ax.tick_params(labelsize=7, colors=MUTED)
ax.text(0.0, -0.185, "rows: mouse·cluster, ordered by assigned stage",
        transform=ax.transAxes, ha="left", fontsize=7, color=MUTED)
pd.concat(rowsrc).to_csv(OUTDIR / f"{NAME}_panelD_marker_scores.tsv", sep="\t")

# ---- (b) 3'UTR length index by stage --------------------------------------
box_rows = []
for k, m in enumerate(MICE):
    ax = fig.add_subplot(gs[1, k])
    d = percell[(percell["mouse"] == m) & (percell["pas_subset"] == "utr3") &
                (percell["depth"] == "raw") & (percell["labels"] == "GEXlabels")]
    st_here = [s for s in STAGES if (d["stage"] == s).sum() >= 30]
    data = [d[d["stage"] == s]["wdi"].values for s in st_here]
    bp = ax.boxplot(data, positions=range(len(st_here)), widths=0.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color=INK, lw=1.4),
                    whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED),
                    boxprops=dict(lw=0.8, edgecolor=MUTED))
    for patch, s in zip(bp["boxes"], st_here):
        patch.set_facecolor(SCOL[s]); patch.set_alpha(0.55)
    means = [np.mean(x) for x in data]
    ax.plot(range(len(st_here)), means, "-o", color=INK, lw=1.6, ms=5, zorder=5,
            label="stage mean")
    ax.set_xticks(range(len(st_here)))
    ax.set_xticklabels([f"{SHORT[s]}\nn={len(x)}" for s, x in zip(st_here, data)], fontsize=8)
    r = trend[(trend.mouse == m) & (trend.pas_subset == "utr3") &
              (trend.depth == "raw") & (trend.labels == "GEXlabels")].iloc[0]
    style(ax, "distal-usage index (WDI)\nlow = shorter 3'UTR", None,
          f"(b{k+1}) {m}: 3'UTR length index across ordered stages (GEX labels)")
    ax.text(0.02, 0.05,
            f"Spearman $\\rho$ = {r.wdi_rho:.3f}   "
            f"p {'< 1e-300' if r.wdi_p == 0 else '= %.1e' % r.wdi_p}\n"
            f"depth-adjusted $\\rho$ = {r.wdi_partial_rho_depthadj:.3f}\n"
            f"Cliff's $\\delta$ (ES vs SPG) = {r.wdi_cliffs:.3f}\n"
            f"label-shuffle null |$\\rho$|max = {r.wdi_null_rho_absmax:.3f}\n"
            f"monotone decreasing: {'YES' if r.wdi_monotone_down else 'no'}",
            transform=ax.transAxes, fontsize=7.6, color=INK, va="bottom",
            bbox=dict(fc="white", ec=GRID, lw=0.6, pad=3))
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", labelcolor=INK)
    for s, x in zip(st_here, data):
        box_rows.append(dict(mouse=m, pas_subset="utr3", labels="GEXlabels", stage=s,
                             n_cells=len(x), mean=float(np.mean(x)), sd=float(np.std(x, ddof=1)),
                             q1=float(np.percentile(x, 25)), median=float(np.median(x)),
                             q3=float(np.percentile(x, 75))))
pd.DataFrame(box_rows).to_csv(OUTDIR / f"{NAME}_panelB_wdi_by_stage.tsv", sep="\t", index=False)

# ---- (c) sensitivity by PAS subset ----------------------------------------
ax = fig.add_subplot(gs[1, 2])
subs = ["all", "atlas", "noatlas", "utr3"]
w = 0.36
crows = []
for j, m in enumerate(MICE):
    for i, sub in enumerate(subs):
        r = trend[(trend.mouse == m) & (trend.pas_subset == sub) &
                  (trend.depth == "raw") & (trend.labels == "GEXlabels")].iloc[0]
        x = i + (j - 0.5) * w
        ax.bar(x, r.wdi_rho, width=w * 0.9, color=SUBCOL[sub],
               alpha=1.0 if j == 0 else 0.55,
               edgecolor=INK if j == 1 else "none", lw=0.7)
        ax.text(x, r.wdi_rho + (0.035 if r.wdi_rho >= 0 else -0.055),
                f"{r.wdi_rho:.2f}", ha="center", fontsize=7,
                color=INK, va="bottom" if r.wdi_rho >= 0 else "top")
        crows.append(dict(mouse=m, pas_subset=sub, n_pas=int(r.n_pas), n_cells=int(r.n_cells),
                          rho=float(r.wdi_rho), p=float(r.wdi_p),
                          rho_depth_adjusted=float(r.wdi_partial_rho_depthadj),
                          rho_depth_matched=float(
                              trend[(trend.mouse == m) & (trend.pas_subset == sub) &
                                    (trend.depth == "matched") &
                                    (trend.labels == "GEXlabels")]["wdi_rho"].iloc[0]),
                          null_abs_rho_max=float(r.wdi_null_rho_absmax),
                          monotone_down=bool(r.wdi_monotone_down)))
nullmax = trend[trend.depth == "raw"]["wdi_null_rho_absmax"].max()
ax.axhspan(-nullmax, nullmax, color=GRID, alpha=0.7, zorder=0)
ax.axhline(0, color=INK, lw=1.0)
ax.set_xticks(range(len(subs)))
ax.set_xticklabels([SUBLAB[s].replace(" (", "\n(") for s in subs], fontsize=7.5)
style(ax, "Spearman $\\rho$ (WDI vs stage rank)", None,
      "(c) the trend lives ONLY in atlas-/3'UTR-supported PAS")
ax.set_ylim(-1.0, 1.0)
ax.text(0.015, 0.965, "$\\rho$ < 0  =  3'UTRs SHORTEN  (the literature-expected direction)",
        transform=ax.transAxes, fontsize=7.6, color=INK, va="top")
ax.text(0.015, 0.915, "grey band = |$\\rho$| from 200 stage-label shuffles",
        transform=ax.transAxes, fontsize=7, color=MUTED, va="top")
ax.legend(handles=[Line2D([], [], marker="s", ls="", color=MUTED, label="mouse1 (solid)"),
                   Line2D([], [], marker="s", ls="", mfc="none", mec=INK, color=INK,
                          label="mouse2 (outlined)")],
          frameon=False, fontsize=7.5, loc="lower left", labelcolor=INK)
pd.DataFrame(crows).to_csv(OUTDIR / f"{NAME}_panelC_subset_sensitivity.tsv", sep="\t", index=False)

# ---- (e) top shortening genes, dumbbell -----------------------------------
ax = fig.add_subplot(gs[2, :2])
g = shortg["utr3"].copy()
top = g.sort_values("mean_delta").head(22).iloc[::-1]
y = np.arange(len(top))
ax.hlines(y, top["m1_delta"], top["m2_delta"], color=GRID, lw=2.4, zorder=1)
ax.scatter(top["m1_delta"], y, s=42, color=C[0], zorder=3, label="mouse1")
ax.scatter(top["m2_delta"], y, s=42, color=C[1], zorder=3, label="mouse2")
ax.axvline(0, color=INK, lw=1.0)
ax.set_yticks(y); ax.set_yticklabels(top.index, fontsize=8)
style(ax, None, "$\\Delta$WDI  (elongating - spermatocyte);  negative = 3'UTR shortening",
      "(e) top shortening genes, replicated in both mice  [3'UTR-restricted PAS]")
ax.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK)
n_com = len(g)
from scipy.stats import spearmanr as _sp
rep_rho = _sp(g["m1_delta"], g["m2_delta"]).statistic
both_s = int(((g.m1_delta < 0) & (g.m2_delta < 0)).sum())
both_l = int(((g.m1_delta > 0) & (g.m2_delta > 0)).sum())
ax.text(0.985, 0.045,
        f"{n_com} genes with >=2 3'UTR PAS measurable in both mice\n"
        f"per-gene $\\Delta$WDI replication: Spearman $\\rho$ = {rep_rho:.3f}\n"
        f"shorten in BOTH mice: {both_s} ({both_s/n_com:.0%})   "
        f"lengthen in both: {both_l} ({both_l/n_com:.0%})",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=INK,
        bbox=dict(fc="white", ec=GRID, lw=0.6, pad=3.5))
g.to_csv(OUTDIR / f"{NAME}_panelE_gene_shortening.tsv", sep="\t")

# ---- caveats block ---------------------------------------------------------
axc = fig.add_subplot(gs[2, 2]); axc.axis("off")
r1 = repval[repval.matrix == "as_shipped"]; r2 = repval[repval.matrix == "repaired"]
cav = (
 "CAVEATS -- read before quoting any number here\n"
 "\n"
 "1. The shipped matrix could not be used.  ema/matrixfilter.py\n"
 "   make_dataframe() returns a full-height matrix but a compacted\n"
 "   pas_id list; annotate() then indexes rows by list POSITION, so\n"
 "   every PAS carries another PAS's counts (drift 0->1,539 mouse1,\n"
 "   0->1,208 mouse2).  Per-PAS counts vs samtools depth over the\n"
 "   same interval: Spearman 0.058 as shipped -> 0.730 re-keyed.\n"
 "   Everything above uses the re-keyed matrix.\n"
 "\n"
 "2. clusters.h5ad .X is TF-IDF, not counts, so `peakatail switch length`\n"
 "   as documented computes PDUI on IDF-weighted values.  And its\n"
 "   classic/per_gene path is strand-blind (proximal/distal swapped\n"
 "   for every minus-strand gene; unit-tested).  The index used here\n"
 "   is computed directly from raw counts, strand-aware.\n"
 "\n"
 "3. Stages are MARKER-INFERRED, not sorted.  The PAS-derived labels\n"
 "   agree with orthogonal STARsolo-GEX labels at 93.9% / 93.8%, and\n"
 "   every trend is reported under both labellings.\n"
 "   Spermatogonia are rare (64 / 68 cells) and the PAS-only\n"
 "   clustering merges them into spermatocytes.\n"
 "\n"
 "4. THE RED FLAG.  With all called PAS the trend is absent or\n"
 "   REVERSED; only 37% of PeakATail PAS are atlas-supported and\n"
 "   36% sit in an annotated 3'UTR.  The atlas-UNsupported peaks\n"
 "   carry a strong opposite-direction signal in both mice.  They are\n"
 "   modestly enriched for an internal-priming signature (A-run>=6 in\n"
 "   +10..+40 nt: 5.4% vs 3.6%) and depleted of a canonical hexamer\n"
 "   (4.2% vs 9.1%) -- suggestive, not conclusive, because PeakATail\n"
 "   peaks are wide intervals rather than point cleavage sites.\n"
 "\n"
 "5. Depth is a confounder (spermatocytes carry ~2x the UMIs of\n"
 "   elongating spermatids).  Trends survive depth-adjusted partial\n"
 "   correlation and binomial thinning to 6,000 UMI/cell.\n"
 "\n"
 "6. PeakATail's CB filter keeps 9,564 / 9,074 barcodes vs STARsolo's\n"
 "   1,294 / 1,364 cells; the extra ~86% are ambient droplets and are\n"
 "   excluded here.\n"
)
axc.text(0.0, 1.0, cav, transform=axc.transAxes, va="top", ha="left",
         fontsize=6.6, color=INK, family="monospace", linespacing=1.30)

fig.savefig(OUTDIR / f"{NAME}.png", dpi=300, facecolor="white")
print("wrote", OUTDIR / f"{NAME}.png")
save_manuscript(fig, NAME, facecolor="white")

# ---- consolidated stats tsv ------------------------------------------------
trend.to_csv(OUTDIR / f"{NAME}_all_trend_stats.tsv", sep="\t", index=False)
pbfix.to_csv(OUTDIR / f"{NAME}_pseudobulk_fixedgeneset.tsv", sep="\t", index=False)
repval.to_csv(OUTDIR / f"{NAME}_matrix_repair_validation.tsv", sep="\t", index=False)
ipdiag.to_csv(OUTDIR / f"{NAME}_internal_priming_diag.tsv", sep="\t", index=False)
shortg["atlas"].to_csv(OUTDIR / f"{NAME}_gene_shortening_atlas.tsv", sep="\t")
print("wrote tsvs to", OUTDIR)
