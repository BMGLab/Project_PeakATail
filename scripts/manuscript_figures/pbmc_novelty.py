#!/usr/bin/env python3
"""
Manuscript figure: pbmc_novelty -- does poly(A)-site space carry cell-type
structure of its own, on independent public data?

The claim under test (the manuscript's central remaining positive claim):
clustering cells on poly(A)-site (PAS) profiles ALONE recovers -- and possibly
extends -- the gene-expression cell-type partition. It was established on the
Laughney cohort (median AMI 0.662 / ARI 0.463 across 17 samples). Here it is
re-tested on 10x pbmc_10k_v3, a public dataset that no part of the method was
tuned on, and then pushed: is any PAS-only population REAL, or is every one of
them a technical or expression artefact?

Panels
  a  UMAPs, GEX-derived and PAS-derived embeddings, both coloured by the same
     marker-derived GEX cell-type label.
  b  contingency PAS cluster x marker-derived GEX type (row-normalised), with
     the four detected splits marked.
  c  concordance of the PAS partition with the GEX partition, PBMC vs the
     Laughney cohort median and range.
  d  what the PAS feature space actually encodes: the identical clustering
     pipeline re-run on feature spaces that strip one kind of information each.
  e  the novelty test: every candidate PAS-only split and the pre-registered
     criterion that killed it.

Inputs are all pre-computed by
  scripts/manuscript_figures/pbmc_novelty_gex.py       (stage 1)
  scripts/manuscript_figures/pbmc_novelty_analysis.py  (stage 2)
  scripts/manuscript_figures/pbmc_novelty_spaces.py    (stage 3)
  scripts/manuscript_figures/pbmc_novelty_atlas_test.py(stage 4)
This script only draws. Every plotted value is written to the companion .tsv.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

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


ROOT = "/mnt/ssd1/Projects/PeakATail_wd"
SRC = os.path.join(ROOT, "results/pbmc_novelty")
OUTDIR = os.path.join(ROOT, "results/figures/manuscript")
NAME = "pbmc_novelty"
os.makedirs(OUTDIR, exist_ok=True)
PNG = os.path.join(OUTDIR, NAME + ".png")
TSV = os.path.join(OUTDIR, NAME + ".tsv")

# ----------------------------------------------------------------------------
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK = "#1B2429"
MUTED = "#5A6B73"
GRID = "#D8E0E3"
SURFACE = "#FFFFFF"

CT_COLOR = {
    "CD4-T": C[0], "CD8-T": C[5], "NK": C[2], "B": C[1],
    "CD14-Mono": C[3], "FCGR3A-Mono": C[4], "DC": MUTED, "Platelet": INK,
}
CT_ORDER = ["CD4-T", "CD8-T", "NK", "B", "CD14-Mono", "FCGR3A-Mono",
            "DC", "Platelet"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
})

# ----------------------------------------------------------------------------
# load
# ----------------------------------------------------------------------------
J = json.load(open(os.path.join(SRC, "novelty_results.json")))
CONC = J["concordance"]
LAU = J["laughney"]
CRIT = J["criteria"]
M = pd.read_csv(os.path.join(SRC, "matched_cells.tsv"), sep="\t")
CT = pd.read_csv(os.path.join(SRC, "contingency_pas_x_celltype.tsv"),
                 sep="\t", index_col=0)
CT.index = CT.index.astype(str)
AT = pd.read_csv(os.path.join(SRC, "atlas_support_tests.tsv"), sep="\t")
FS = pd.read_csv(os.path.join(SRC, "feature_spaces.tsv"), sep="\t")
KM = pd.read_csv(os.path.join(SRC, "feature_spaces_kmatched.tsv"), sep="\t")
VM = pd.read_csv(os.path.join(SRC, "variance_matched.tsv"), sep="\t")
UF = json.load(open(os.path.join(SRC, "usage_features_summary.json")))
PM = pd.read_csv(os.path.join(SRC, "power_matched.tsv"), sep="\t")
PWR = json.load(open(os.path.join(SRC, "power_matched.json")))
LAUGH = pd.read_csv(
    os.path.join(OUTDIR, "clustering_concordance.tsv"), sep="\t")
LAUGH = LAUGH[LAUGH.panel == "a"].pivot_table(
    index="gsm", columns="metric", values="value")

rows_tsv = []


def rec(panel, record, key, value, note=""):
    rows_tsv.append({"panel": panel, "record": record, "key": key,
                     "value": value, "note": note})


# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(12.2, 12.9), facecolor=SURFACE)
gs = fig.add_gridspec(
    3, 3, height_ratios=[0.98, 1.05, 1.00], width_ratios=[1.0, 1.0, 1.16],
    left=0.075, right=0.978, top=0.925, bottom=0.175, hspace=0.46, wspace=0.62)

# ============================================================ panel a: UMAPs
def umap_panel(ax, xk, yk, title, sub):
    for ct in CT_ORDER:
        s = M[M.celltype == ct]
        if not len(s):
            continue
        ax.scatter(s[xk], s[yk], s=0.55, lw=0, alpha=0.55,
                   color=CT_COLOR[ct], rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_title(title, fontsize=8.5, color=INK, pad=3, loc="left")
    ax.text(0.0, 1.005, "", transform=ax.transAxes)
    ax.set_xlabel(sub, fontsize=6.8, color=MUTED, labelpad=2)


axA1 = fig.add_subplot(gs[0, 0])
umap_panel(axA1, "gex_umap1", "gex_umap2",
           "a   GEX-derived embedding",
           "10x gene counts → HVG/PCA/Leiden;  colour = marker-derived label")
axA2 = fig.add_subplot(gs[0, 1])
umap_panel(axA2, "pas_umap1", "pas_umap2",
           "     PAS-derived embedding",
           "275,370 PAS counts → TF-IDF/LSI/Leiden;  no gene expression used")

handles = [Line2D([], [], marker="o", ls="", ms=4.2, color=CT_COLOR[c],
                  label=f"{c} ({int((M.celltype == c).sum())})")
           for c in CT_ORDER if (M.celltype == c).any()]
axA2.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.075),
            ncol=4, frameon=False, handletextpad=0.3, columnspacing=1.0,
            fontsize=6.3, borderpad=0.0)
for ct in CT_ORDER:
    rec("a", ct, "n_cells", int((M.celltype == ct).sum()),
        "marker-derived label, matched cells")

# ==================================================== panel c: concordance
axC = fig.add_subplot(gs[0, 2])
pairs = [("AMI\nvs cell type", "AMI_pas_vs_celltype", "AMI_celltype_vs_pas"),
         ("ARI\nvs cell type", "ARI_pas_vs_celltype", "ARI_celltype_vs_pas"),
         ("AMI vs\nGEX Leiden", "AMI_pas_vs_gexleiden", None),
         ("ARI vs\nGEX Leiden", "ARI_pas_vs_gexleiden", None)]
LAU_MED = {"AMI_celltype_vs_pas": LAU["median_AMI_celltype_vs_pas"],
           "ARI_celltype_vs_pas": LAU["median_ARI_celltype_vs_pas"],
           None: None}
lau_gex = {"AMI_pas_vs_gexleiden": LAU["median_AMI_gexleiden_vs_pas"],
           "ARI_pas_vs_gexleiden": LAU["median_ARI_gexleiden_vs_pas"]}
x = np.arange(len(pairs))
w = 0.36
for i, (lab, kk, lk) in enumerate(pairs):
    v = CONC[kk]
    axC.bar(i - w / 2, v, w, color=C[0], lw=0, zorder=3)
    axC.text(i - w / 2, v + 0.012, f"{v:.3f}", ha="center", va="bottom",
             fontsize=6.6, color=C[0], fontweight="bold")
    rec("c", "pbmc_10k_v3", kk, round(float(v), 4), "this dataset")
    if lk is not None:
        med = LAU[f"median_{lk}"]
        lo, hi = LAU[f"min_{lk}"], LAU[f"max_{lk}"]
        vals = LAUGH[lk].values
        axC.bar(i + w / 2, med, w, color=GRID, edgecolor=MUTED, lw=0.6,
                zorder=3)
        axC.vlines(i + w / 2, lo, hi, color=MUTED, lw=1.0, zorder=4)
        axC.scatter(np.full(len(vals), i + w / 2) +
                    np.random.RandomState(0).uniform(-0.08, 0.08, len(vals)),
                    vals, s=5, color=MUTED, alpha=0.7, lw=0, zorder=5)
        axC.text(i + w / 2, hi + 0.012, f"{med:.3f}", ha="center",
                 va="bottom", fontsize=6.6, color=MUTED)
        rec("c", "laughney_n17", f"median_{lk}", round(float(med), 4),
            f"range {lo:.3f}-{hi:.3f}")
    else:
        med = lau_gex[kk]
        axC.bar(i + w / 2, med, w, color=GRID, edgecolor=MUTED, lw=0.6,
                zorder=3)
        axC.text(i + w / 2, med + 0.012, f"{med:.3f}", ha="center",
                 va="bottom", fontsize=6.6, color=MUTED)
        rec("c", "laughney_n17", kk + "_median", round(float(med), 4),
            "cohort median")
axC.set_xticks(x)
axC.set_xticklabels([p[0] for p in pairs], fontsize=6.8)
axC.set_ylim(0, 1.06)
axC.set_ylabel("agreement of PAS clusters with the GEX partition")
axC.grid(axis="y", color=GRID, lw=0.6, zorder=0)
axC.set_axisbelow(True)
for sp in ("top", "right"):
    axC.spines[sp].set_visible(False)
axC.set_title("c   the Laughney claim replicates", fontsize=8.5, loc="left",
              pad=3)
axC.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=5, color=C[0],
           label="pbmc_10k_v3 (this work)"),
    Line2D([], [], marker="s", ls="", ms=5, color=GRID,
           markeredgecolor=MUTED, label="Laughney median (n=17, dots=samples)")],
    loc="upper right", frameon=False, fontsize=6.3, handletextpad=0.4)

# =========================================== panel b: contingency heatmap
axB = fig.add_subplot(gs[1, :2])
ctm = CT.reindex(columns=[c for c in CT_ORDER if c in CT.columns])
frac = ctm.div(ctm.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
im = axB.imshow(frac.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
axB.set_xticks(range(frac.shape[1]))
axB.set_xticklabels(frac.columns, rotation=38, ha="right", fontsize=6.8)
axB.set_yticks(range(frac.shape[0]))
axB.set_yticklabels([f"{i}  (n={int(ctm.loc[i].sum())})" for i in frac.index],
                    fontsize=6.2)
axB.set_ylabel("PAS cluster (Leiden, PAS counts only)")
for i in range(frac.shape[0]):
    for j in range(frac.shape[1]):
        n = int(ctm.values[i, j])
        if n == 0:
            continue
        axB.text(j, i, f"{n}", ha="center", va="center", fontsize=5.6,
                 color=("white" if frac.values[i, j] > 0.55 else INK))
        rec("b", f"pas{frac.index[i]}", frac.columns[j], n, "cells")
splits = {s["celltype"]: s["pas_clusters"][:2] for s in J["splits"]}
for ct, cls in splits.items():
    j = list(frac.columns).index(ct)
    for cl in cls:
        i = list(frac.index).index(str(cl))
        axB.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                edgecolor=C[1], lw=1.6, zorder=6))
axB.set_title("b   PAS cluster × marker-derived GEX type   "
              "(orange = the four detected splits)", fontsize=8.5, loc="left",
              pad=4)
cax = axB.inset_axes([0.66, -0.235, 0.32, 0.022])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("row fraction (composition of the PAS cluster)", fontsize=6.0,
             labelpad=2)
cb.ax.tick_params(labelsize=5.8, length=2, pad=1)
cb.outline.set_edgecolor(GRID)

# =================================== panel d: what the feature space encodes
axD = fig.add_subplot(gs[1, 2])
fs = FS.set_index("space")
km = KM.set_index("space")
vm = VM.set_index("space")
pm = PM.set_index("space")


def best_of(tbl, a, b, col):
    """Generous reading: report the better of the raw and column-scaled
    variant of a continuous-feature space, so no result here can be blamed on
    having crippled the representation."""
    va = tbl.loc[a, col] if a in tbl.index else -np.inf
    vb = tbl.loc[b, col] if b in tbl.index else -np.inf
    return float(max(va, vb))


BARS = [
    # label, AMI at the pipeline's own resolution, AMI at k matched to 8 types,
    # colour key, n features, note
    ("all PAS, counts  (shipped)", float(CONC["AMI_pas_vs_celltype"]),
     float(km.loc["pas_counts", "AMI_vs_celltype_kmatched"]), "counts",
     275370, "the published result"),
    ("summed per gene = EXPRESSION",
     float(fs.loc["gene_level", "AMI_vs_celltype"]),
     float(km.loc["gene_level", "AMI_vs_celltype_kmatched"]), "expr",
     14891, "site resolution discarded"),
    ("real PAS only, counts",
     float(fs.loc["atlas_only", "AMI_vs_celltype"]),
     float(km.loc["atlas_only", "AMI_vs_celltype_kmatched"]), "real",
     32361, ""),
    ("NON-PAS peaks only, counts",
     float(fs.loc["nonatlas_only", "AMI_vs_celltype"]),
     float(km.loc["nonatlas_only", "AMI_vs_celltype_kmatched"]), "fake",
     243009, ""),
    ("gene-DETECTION mask only",
     float(fs.loc["detection_mask", "AMI_vs_celltype"]),
     float(km.loc["detection_mask", "AMI_vs_celltype_kmatched"]), "mask",
     10000, "no usage value at all"),
    ("USAGE, top 10k sites",
     best_of(fs, "usage_top10k", "usage_top10k_scaled", "AMI_vs_celltype"),
     best_of(km, "usage_top10k", "usage_top10k_scaled",
             "AMI_vs_celltype_kmatched"), "fake", 10000,
     f"{1 - UF['atlas_frac_selected']:.0%} of them are NOT poly(A) sites"),
    ("USAGE, REAL poly(A) sites only",
     best_of(vm, "real_pas_usage_top", "real_pas_usage_top_scaled", "AMI_res1"),
     best_of(vm, "real_pas_usage_top", "real_pas_usage_top_scaled",
             "AMI_kmatched"), "real", 10000, "true APA signal"),
    ("USAGE, NON-PAS (var-matched)",
     best_of(vm, "decoy_usage_var_matched", "decoy_usage_var_matched_scaled",
             "AMI_res1"),
     best_of(vm, "decoy_usage_var_matched", "decoy_usage_var_matched_scaled",
             "AMI_kmatched"), "fake", 10000, "same n, same usage variance"),
    ("USAGE, NON-PAS (var + depth)",
     best_of(pm, "decoy_usage_depth_and_var_matched",
             "decoy_usage_depth_and_var_matched_scaled", "AMI_res1"),
     best_of(pm, "decoy_usage_depth_and_var_matched",
             "decoy_usage_depth_and_var_matched_scaled", "AMI_kmatched"),
     "fake", 10000,
     f"non-PAS counts thinned to p={PWR['thinning_p']:.2f} so per-gene depth "
     f"matches the real-PAS set"),
    ("USAGE shuffled (null)",
     best_of(fs, "usage_shuffled", "usage_shuffled_scaled", "AMI_vs_celltype"),
     best_of(km, "usage_shuffled", "usage_shuffled_scaled",
             "AMI_vs_celltype_kmatched"), "null", 10000, ""),
]
KEYCOL = {"counts": C[0], "expr": C[3], "real": C[2], "fake": C[4],
          "mask": MUTED, "null": "#B9C4C9"}
ypos = np.arange(len(BARS))[::-1]
for y, (lab, v1, vk, key, nf, note) in zip(ypos, BARS):
    axD.barh(y, v1, 0.62, color=KEYCOL[key], lw=0, zorder=3)
    axD.scatter([vk], [y], s=17, facecolor="white", edgecolor=INK, lw=0.9,
                zorder=5)
    axD.text(max(v1, vk) + 0.012, y, f"{v1:.2f}", va="center", fontsize=6.3,
             color=INK, fontweight="bold")
    rec("d", lab, "AMI_vs_celltype_res1.0", round(v1, 4),
        f"n_features={nf}; {note}")
    rec("d", lab, "AMI_vs_celltype_k_matched_to_8_types", round(vk, 4), "")
axD.set_yticks(ypos)
axD.set_yticklabels([b[0] for b in BARS], fontsize=6.4)
axD.set_xlabel("AMI vs marker-derived cell type")
axD.set_xlim(0, 1.0)
axD.grid(axis="x", color=GRID, lw=0.6, zorder=0)
axD.set_axisbelow(True)
for sp in ("top", "right"):
    axD.spines[sp].set_visible(False)
axD.axvline(float(CONC["AMI_pas_vs_celltype"]), color=C[0], lw=0.8, ls=":",
            zorder=2)
axD.set_title("d   what the PAS space actually encodes\n"
              "identical pipeline, one kind of information\n"
              "removed each time; Leiden resolution 1.0",
              fontsize=8.5, loc="left", pad=4)
axD.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=5, color=C[2],
           label="features are real PAS (atlas ≤100 bp)"),
    Line2D([], [], marker="s", ls="", ms=5, color=C[4],
           label="features are NOT PAS"),
    Line2D([], [], marker="o", ls="", ms=4.2, markerfacecolor="white",
           markeredgecolor=INK, label="same space, k matched to 8 types")],
    loc="upper left", bbox_to_anchor=(0.0, -0.135), ncol=1, frameon=False,
    fontsize=6.0, handletextpad=0.4, borderpad=0.0, labelspacing=0.35)
axD.text(0.0, -0.315,
         "Dotted line = the shipped result.",
         transform=axD.transAxes, fontsize=6.0, color=MUTED, va="top")

# ================================================= panel e: the novelty test
axE = fig.add_subplot(gs[2, :])
axE.axis("off")
axE.set_title("e   novelty test — does any PAS-only split survive? "
              "Pre-registered criteria, applied in order",
              fontsize=8.5, loc="left", pad=2)

CRITS = [
    ("depth\n(PAS UMI)", "depth"),
    ("depth\n(GEX UMI)", "gexdepth"),
    ("mito %", "mito"),
    ("doublet\nscore", "doublet"),
    ("ambient /\noff-lineage", "ambient"),
    ("plain\nexpression", "expression"),
    ("GEX Leiden\nalready splits", "gexres"),
    ("discriminating\nsites are real PAS", "atlas"),
]
at = AT.set_index("split")
tags = [t for t in J["apa"] if t != "POOLED"]
n_r, n_c = len(tags), len(CRITS)
cw, ch = 1.0 / (n_c + 2.45), 1.0 / (n_r + 3.4)
x0, y0 = 0.245, 0.80

for j, (lab, _) in enumerate(CRITS):
    axE.text(x0 + (j + 0.5) * cw, y0 + 0.055, lab, ha="center", va="bottom",
             fontsize=6.1, color=MUTED, linespacing=1.15)

for i, tag in enumerate(tags):
    r = J["apa"][tag]
    a = r["audit"]
    y = y0 - (i + 0.5) * ch
    axE.text(x0 - 0.012, y, f"{r['celltype']}   PAS {r['pas_A']} (n={r['nA']}) "
             f"vs PAS {r['pas_B']} (n={r['nB']})",
             ha="right", va="center", fontsize=6.6, color=INK)
    ar = at.loc[tag]
    cells = [
        (a["depth_confounded"],
         f"{2**abs(a['log2ratio_pas_counts']):.1f}×\nAUC {ar['depth_AUC_pas']:.2f}"),
        (a["gexdepth_confounded"], f"{2**abs(a['log2ratio_gex_counts']):.1f}×"),
        (a["mito_confounded"], f"{2**abs(a['log2ratio_pct_mt']):.1f}×"),
        (a["doublet_confounded"],
         f"{2**abs(a['log2ratio_doublet_score']):.1f}×"),
        (a["ambient_confounded"],
         f"{2**abs(a['ambient_max_abs_log2']):.1f}×\n{a['ambient_worst_panel'][4:]}"),
        (a["expression_confounded"], f"{a['n_gex_de_q05_lfc1']} DE\ngenes"),
        (a["gex_already_resolves"], f"ARI\n{a['ARI_arm_vs_gexleiden']:.2f}"),
        (float(ar["atlas_frac_apa_driven"]) <
         float(ar["atlas_frac_candidate_universe"]),
         f"{ar['n_apa_driven_atlas']:.0f}/{ar['n_apa_driven']:.0f}\n"
         f"={ar['atlas_frac_apa_driven']:.0%}"),
    ]
    for j, (bad, txt) in enumerate(cells):
        xc = x0 + j * cw
        axE.add_patch(Rectangle((xc + 0.004, y - ch * 0.44), cw - 0.008,
                                ch * 0.88,
                                facecolor=("#F6DCCE" if bad else "#DFF0EA"),
                                edgecolor=(C[1] if bad else C[2]), lw=0.7,
                                transform=axE.transAxes, zorder=2))
        axE.text(xc + cw / 2, y, txt, ha="center", va="center", fontsize=5.5,
                 color=(C[1] if bad else "#0B5E4A"), linespacing=1.1, zorder=3)
        rec("e", tag, CRITS[j][1], ("FAIL" if bad else "pass"),
            txt.replace("\n", " "))
    rec("e", tag, "n_candidates", r["n_candidates"], "")
    rec("e", tag, "n_significant_q05", r["n_significant"], "per-cell MWU, BH")
    rec("e", tag, "n_apa_driven", r["n_apa_driven"], "sig & expr comparable")
    rec("e", tag, "n_apa_driven_atlas_supported", r["n_apa_driven_atlas"], "")
    rec("e", tag, "null_sig_mean", round(r["null_sig_mean"], 2),
        f"{CRIT['N_PERM']} within-split relabellings")

pool = AT[AT.split == "POOLED"].iloc[0]
verdict = (
    "VERDICT: no split survives.  Every one of the four PAS-only splits fails at "
    "least one pre-registered criterion.\n"
    f"Pooled across splits, {int(pool['n_apa_driven'])} PAS pass the APA-driven "
    f"filter, but only {int(pool['n_apa_driven_atlas'])} "
    f"({pool['atlas_frac_apa_driven']:.1%}) are within {CRIT['ATLAS_WIN']} bp of a "
    f"PolyASite 2.0 site — significantly FEWER than the "
    f"{pool['atlas_frac_candidate_universe']:.1%} of the sites they were drawn "
    f"from\n(Fisher OR {pool['OR_vs_candidates']:.2f}, p = "
    f"{pool['p_fisher_vs_candidates']:.3f}), and no better than the "
    f"{pool['atlas_frac_all_pas']:.1%} background of all PeakATail PBMC calls.  "
    "The 'APA-driven' signal is carried by peaks that are not poly(A) sites.")
axE.add_patch(Rectangle((0.0, 0.0), 1.0, 0.205, transform=axE.transAxes,
                        facecolor="#FBF0EA", edgecolor=C[1], lw=0.9, zorder=1))
axE.text(0.012, 0.1025, verdict, transform=axE.transAxes, fontsize=6.7,
         color=INK, va="center", ha="left", linespacing=1.45, zorder=2)
rec("e", "POOLED", "n_apa_driven", int(pool["n_apa_driven"]), "")
rec("e", "POOLED", "n_apa_driven_atlas_supported",
    int(pool["n_apa_driven_atlas"]), "")
rec("e", "POOLED", "atlas_frac_apa_driven",
    round(float(pool["atlas_frac_apa_driven"]), 4), "")
rec("e", "POOLED", "atlas_frac_candidate_universe",
    round(float(pool["atlas_frac_candidate_universe"]), 4), "")
rec("e", "POOLED", "fisher_OR_vs_candidates",
    round(float(pool["OR_vs_candidates"]), 4), "")
rec("e", "POOLED", "fisher_p_vs_candidates",
    float(pool["p_fisher_vs_candidates"]), "")

crit_txt = (
    f"split: ≥2 PAS clusters each ≥{CRIT['SPLIT_MIN_PURITY']:.0%} one type, "
    f"each ≥{CRIT['SPLIT_MIN_FRAC']:.0%} of it and ≥{CRIT['SPLIT_MIN_CELLS']} "
    f"cells  |  candidate site: |Δusage| ≥ {CRIT['CAND_MIN_DPROP']}, gene "
    f"≥{CRIT['CAND_MIN_GENE_COUNTS']} counts/arm  |  significant: per-CELL "
    f"usage, Mann–Whitney, BH q<{CRIT['FDR_Q']}\n"
    f"APA-driven: significant AND |log2FC| ≤ {CRIT['EXPR_MAX_ABS_LOG2FC']} for "
    f"the gene's TOTAL expression in the independent 10x GEX matrix  |  "
    f"confound flagged at >1.5× median ratio, expression at ≥10 DE genes "
    f"(q<0.05, |log2FC|>1), GEX-resolved at ARI ≥ {CRIT['GEX_RESOLVED_ARI']}")
axE.text(0.0, -0.075, crit_txt, transform=axE.transAxes, fontsize=5.9,
         color=MUTED, va="top", ha="left", linespacing=1.45)

# ============================================================== title/caveats
fig.text(0.075, 0.988,
         "Poly(A)-site clustering reproduces gene-expression cell types on "
         "public PBMC data.",
         fontsize=12.0, color=INK, fontweight="bold", ha="left", va="top")
fig.text(0.075, 0.968,
         "The signal is 3′-end expression, not poly(A)-site choice, and no "
         "PAS-only population survives.",
         fontsize=12.0, color=C[1], fontweight="bold", ha="left", va="top")
fig.text(0.075, 0.949,
         "10x pbmc_10k_v3 (11,836 PAS cells / 11,043 GEX cells; "
         f"{CONC['n_shared_cells']:,} matched, "
         f"{CONC['coverage_pas_by_gex']:.1%} of PAS cells). "
         "PeakATail PAS calls, TF-IDF/LSI/Leiden on PAS counts only. "
         "Independent replication of the Laughney-cohort result.",
         fontsize=7.6, color=MUTED, ha="left", va="top")

import textwrap as _tw
caveats = _tw.fill(
    "CAVEATS ON EVERY NUMBER ABOVE.  (1) The reference partition is "
    "MARKER-DERIVED, not curated ground truth: canonical PBMC marker panels "
    "scored per GEX-Leiden cluster, modal label; every error of that panel is "
    "inherited here.  (2) The shipped PAS features are per-site COUNTS, which "
    "carry gene abundance; panel d shows that collapsing them to gene totals "
    "costs nothing, so the concordance is a 3′-end EXPRESSION result and is "
    "not evidence of isoform choice.  (3) On this exact dataset PeakATail is "
    "LAST among de novo callers for PAS accuracy (F1 0.134 vs polyApipe "
    f"0.261) and only {CONC['atlas_frac_all_pas']:.1%} of its 275,370 PBMC "
    "sites lie within 100 bp of a PolyASite 2.0 representative site; that is "
    "why the atlas column of panel e and the real-PAS bars of panel d are "
    "decisive.  (4) Absence of surviving splits is not proof that none exist: "
    "splits are scored on the two largest type-dominated PAS clusters of each "
    "type, and the atlas is itself incomplete.  (5) One dataset, one donor, "
    "one parameter set.  (6) In panel d the two continuous usage spaces are "
    "each shown at the better of their raw and column-scaled representation, "
    "so the collapse of the real-poly(A)-site bar is not an artefact of "
    "representation; the open circles repeat every space with its Leiden "
    "resolution retuned so the cluster count matches the 8 cell types, which "
    "removes granularity as an explanation.  (7) Real poly(A) sites are "
    "sparser here (3 sites/gene, gene detected in "
    f"{PWR['real_det_rate']:.1%} of cells) than PeakATail's non-PAS peaks "
    f"(10 sites/gene, {PWR['nonpas_det_rate']:.1%}), so the last bar repeats "
    "the non-PAS decoy after binomially thinning its counts to the same "
    f"per-gene depth ({PWR['nonpas_thinned_det_rate']:.1%} detection, median "
    f"gene total {PWR['nonpas_thinned_med_gene_total']:.0f} vs "
    f"{PWR['real_med_gene_total']:.0f}): it still doubles the real-PAS score, "
    "so measurement power is not the explanation.", width=178)
fig.text(0.075, 0.108, caveats, fontsize=6.4, color=MUTED, ha="left",
         va="top", linespacing=1.55)

fig.savefig(PNG, dpi=300, facecolor=SURFACE)
print("wrote", PNG)
save_manuscript(fig, NAME, facecolor=SURFACE)

pd.DataFrame(rows_tsv).to_csv(TSV, sep="\t", index=False)
print("wrote", TSV)
