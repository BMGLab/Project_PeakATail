#!/usr/bin/env python3
"""
Manuscript figure: GEX-vs-PAS clustering concordance (PeakATail / ema).

Core claim carried by this figure: cells clustered on poly(A)-site (PAS)
profiles ALONE -- no gene-level expression matrix -- recover the cell-type
structure that gene-expression marker scoring assigns.

NOTE ON SCOPE: the PAS feature space is per-site COUNTS (TF-IDF + LSI, see
07_clustering/*/clusters.h5ad), not within-gene usage fractions. Per-site counts
carry gene abundance, so this figure shows that PAS-resolved quantification
recovers the GEX partition -- it does NOT isolate isoform-choice signal from
expression-level signal. Both limitations are stated on the figure itself.

Panels
  (a) per-sample AMI and ARI of PAS clusters vs GEX cell types, sorted by AMI,
      cohort medians annotated.
  (b) n_pas_clusters vs n_celltypes with the y = x line: does PAS space find
      structure of comparable granularity?
  (c) the honest caveat: coverage = matched_cells / pas_cells per sample. PAS
      clusters contain many cells the GEX side never labelled, so concordance is
      scored on a subset only.

Everything plotted is also written to the companion .tsv (tidy/long) so every
value on the page is auditable.

Inputs are READ-ONLY (/mnt/ssd2). Outputs go to
/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript/
"""

import json
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

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
# paths
# ----------------------------------------------------------------------------
SRC = (
    "/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/"
    "runs/B1_cohort_full/B2_gex_celltyping"
)
CSV = os.path.join(SRC, "concordance.csv")
SUMMARY = os.path.join(SRC, "summary.json")

OUTDIR = "/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript"
NAME = "clustering_concordance"
os.makedirs(OUTDIR, exist_ok=True)
PNG = os.path.join(OUTDIR, NAME + ".png")
TSV = os.path.join(OUTDIR, NAME + ".tsv")

# ----------------------------------------------------------------------------
# validated palette (fixed hexes, fixed entity assignment)
# ----------------------------------------------------------------------------
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK = "#1B2429"
MUTED = "#5A6B73"
GRID = "#D8E0E3"
SURFACE = "#FFFFFF"

STAGE_ORDER = ["Normal", "Stage I–II", "Stage IV primary", "Metastasis"]
STAGE_COLOR = {
    "Normal": C[0],
    "Stage I–II": C[1],
    "Stage IV primary": C[2],
    "Metastasis": C[3],
}

LABEL_TO_GROUP = {
    "Normal": "Normal",
    "StageIA": "Stage I–II",
    "StageIB": "Stage I–II",
    "StageIIA": "Stage I–II",
    "StageIVprimary": "Stage IV primary",
    "MetBone": "Metastasis",
    "MetBrain": "Metastasis",
    "MetAdrenal": "Metastasis",
}

plt.rcParams.update(
    {
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
    }
)


def tidy_axes(ax, value_axis="y"):
    """Hide top/right spines; gridlines on the VALUE axis only, below the data."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(
        True,
        axis=value_axis,
        color=GRID,
        linewidth=0.6,
        zorder=0,
    )
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.8)


# ----------------------------------------------------------------------------
# load + derive
# ----------------------------------------------------------------------------
df = pd.read_csv(CSV)
assert len(df) == 17, f"expected 17 samples, got {len(df)}"

df["stage_label"] = df["gsm"].str.split("-").str[1]
unknown = set(df["stage_label"]) - set(LABEL_TO_GROUP)
assert not unknown, f"unmapped stage labels: {unknown}"
df["stage_group"] = df["stage_label"].map(LABEL_TO_GROUP)

# integrity checks against the raw numbers
df["coverage_recomputed"] = df["matched_cells"] / df["pas_cells"]
cov_dev = float((df["coverage_recomputed"] - df["coverage"]).abs().max())
matched_eq_gex = bool((df["matched_cells"] == df["gex_cells"]).all())

with open(SUMMARY) as fh:
    summ = json.load(fh)
sj = pd.DataFrame(
    [
        {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}
        for rec in summ.values()
    ]
)
merged = df.merge(sj, on="gsm", suffixes=("", "_json"))
json_agrees = bool(
    np.allclose(
        merged["AMI_celltype_vs_pas"], merged["AMI_celltype_vs_pas_json"], atol=1e-9
    )
    and np.allclose(
        merged["ARI_celltype_vs_pas"], merged["ARI_celltype_vs_pas_json"], atol=1e-9
    )
)

# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------
med_ami = float(df["AMI_celltype_vs_pas"].median())
med_ari = float(df["ARI_celltype_vs_pas"].median())
med_cov = float(df["coverage"].median())
med_ami_leiden = float(df["AMI_gexleiden_vs_pas"].median())
med_ari_leiden = float(df["ARI_gexleiden_vs_pas"].median())

groups_all = [
    df.loc[df["stage_group"] == g, "AMI_celltype_vs_pas"].to_numpy()
    for g in STAGE_ORDER
]
groups_n3 = [g for g, name in zip(groups_all, STAGE_ORDER) if len(g) >= 2]
kw_ami = stats.kruskal(*groups_n3)
kw_ari = stats.kruskal(
    *[
        df.loc[df["stage_group"] == g, "ARI_celltype_vs_pas"].to_numpy()
        for g in STAGE_ORDER
        if (df["stage_group"] == g).sum() >= 2
    ]
)
mw_norm_met = stats.mannwhitneyu(
    df.loc[df["stage_group"] == "Normal", "AMI_celltype_vs_pas"],
    df.loc[df["stage_group"] == "Metastasis", "AMI_celltype_vs_pas"],
    alternative="two-sided",
)

# n_pas_clusters in concordance.csv counts PAS clusters SURVIVING the restriction
# to GEX-labelled cells. Recover the count in the full PAS clustering so panel (b)
# can state that the plotted value is the conservative one.
CLUSTDIR = os.path.join(os.path.dirname(SRC), "07_clustering")
mean_npas_full = None
try:
    import h5py
    from anndata.io import read_elem

    full_counts = []
    for g in df["gsm"]:
        with h5py.File(os.path.join(CLUSTDIR, g, "clusters.h5ad"), "r") as h:
            full_counts.append(read_elem(h["obs"])["leiden"].astype(str).nunique())
    mean_npas_full = float(np.mean(full_counts))
except Exception as exc:  # source volume unavailable -> omit the clause
    print(f"note: full PAS cluster counts unavailable ({exc})")

rho_cov_ami = stats.spearmanr(df["coverage"], df["AMI_celltype_vs_pas"])
rho_cov_ari = stats.spearmanr(df["coverage"], df["ARI_celltype_vs_pas"])
rho_clusters = stats.spearmanr(df["n_celltypes"], df["n_pas_clusters"])
n_pas_above = int((df["n_pas_clusters"] > df["n_celltypes"]).sum())
n_pas_equal = int((df["n_pas_clusters"] == df["n_celltypes"]).sum())
n_pas_below = int((df["n_pas_clusters"] < df["n_celltypes"]).sum())

# ----------------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 7.8), facecolor=SURFACE)
gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.40, 1.00],
    left=0.190,
    right=0.965,
    top=0.820,
    bottom=0.070,
    hspace=0.32,
    wspace=0.24,
)
ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])

# ---------------------------------------------------------------- panel (a) --
a = df.sort_values("AMI_celltype_vs_pas", ascending=True).reset_index(drop=True)
ypos = np.arange(len(a))

for i, r in a.iterrows():
    col = STAGE_COLOR[r["stage_group"]]
    ax_a.plot(
        [r["ARI_celltype_vs_pas"], r["AMI_celltype_vs_pas"]],
        [i, i],
        color=col,
        lw=1.8,
        solid_capstyle="round",
        zorder=2,
        alpha=0.85,
    )
    # ARI = open marker, AMI = filled marker (shape carries the metric)
    ax_a.plot(
        r["ARI_celltype_vs_pas"],
        i,
        "o",
        ms=6.0,
        mfc=SURFACE,
        mec=col,
        mew=1.6,
        zorder=3,
    )
    ax_a.plot(
        r["AMI_celltype_vs_pas"],
        i,
        "o",
        ms=6.5,
        mfc=col,
        mec=SURFACE,
        mew=0.9,
        zorder=4,
    )

ax_a.axvline(med_ami, color=INK, lw=1.0, ls=(0, (4, 2)), zorder=1, alpha=0.75)
ax_a.axvline(med_ari, color=MUTED, lw=1.0, ls=(0, (2, 2)), zorder=1, alpha=0.85)

ax_a.set_yticks(ypos)
ax_a.set_yticklabels(a["gsm"], fontsize=6)
for tick, grp in zip(ax_a.get_yticklabels(), a["stage_group"]):
    tick.set_color(STAGE_COLOR[grp])
ax_a.set_ylim(-0.85, len(a) - 0.15)
# x starts at 0 because the axis label declares 0 = chance and the claim is
# "no sample near chance" -- the reader must be able to see the real distance.
ax_a.set_xlim(0.0, 0.86)
ax_a.set_xticks(np.arange(0.0, 0.81, 0.1))
ax_a.set_xlabel(
    "agreement with GEX cell-type labels (0 = chance, 1 = identical partitions)"
)
tidy_axes(ax_a, value_axis="x")

# median call-outs, sitting above the top row
ax_a.text(
    med_ami,
    len(a) - 0.42,
    f" median AMI {med_ami:.2f}",
    fontsize=6.8,
    color=INK,
    ha="left",
    va="center",
)
ax_a.text(
    med_ari,
    len(a) - 0.42,
    f"median ARI {med_ari:.2f} ",
    fontsize=6.8,
    color=MUTED,
    ha="right",
    va="center",
)

# direct labels on the top (best) sample, so the two marks are named in place
top = a.iloc[-1]
ax_a.annotate(
    "AMI",
    xy=(top["AMI_celltype_vs_pas"], len(a) - 1),
    xytext=(top["AMI_celltype_vs_pas"] + 0.035, len(a) - 1.95),
    fontsize=6.8,
    color=MUTED,
    ha="left",
    va="center",
    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=0, shrinkB=3),
)
ax_a.annotate(
    "ARI",
    xy=(top["ARI_celltype_vs_pas"], len(a) - 1),
    xytext=(top["ARI_celltype_vs_pas"] - 0.035, len(a) - 1.95),
    fontsize=6.8,
    color=MUTED,
    ha="right",
    va="center",
    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=0, shrinkB=3),
)

ax_a.legend(
    handles=[
        Line2D(
            [], [], marker="o", ls="", ms=6.5, mfc=MUTED, mec=SURFACE, mew=0.9,
            label="AMI (adjusted mutual information)",
        ),
        Line2D(
            [], [], marker="o", ls="", ms=6.0, mfc=SURFACE, mec=MUTED, mew=1.6,
            label="ARI (adjusted Rand index)",
        ),
    ],
    loc="lower right",
    frameon=True,
    facecolor=SURFACE,
    edgecolor=GRID,
    framealpha=1.0,
    fontsize=6.5,
    handletextpad=0.4,
    borderpad=0.5,
)

ax_a.set_title(
    "a   PAS-only clusters recover GEX cell types in all 17 samples "
    f"(median AMI {med_ami:.2f}, ARI {med_ari:.2f})",
    loc="left",
    fontsize=8.5,
    color=INK,
    pad=8,
)
ax_a.text(
    0.012,
    0.845,
    f"Kruskal–Wallis across stage groups: AMI p = {kw_ami.pvalue:.2f}\n"
    "no stage effect detected (n = 4/7/5; low power)",
    transform=ax_a.transAxes,
    fontsize=6.3,
    color=MUTED,
    ha="left",
    va="top",
    linespacing=1.45,
    bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.6),
)

# ---------------------------------------------------------------- panel (b) --
lo, hi = 5, 18
ax_b.plot(
    [lo, hi], [lo, hi], color=MUTED, lw=1.0, ls=(0, (4, 2)), zorder=1
)
ax_b.text(
    16.6, 16.9, "y = x", fontsize=6.5, color=MUTED, ha="center", va="bottom",
    rotation=45, rotation_mode="anchor",
)

rng = np.random.default_rng(0)
jx = rng.uniform(-0.13, 0.13, len(df))
jy = rng.uniform(-0.13, 0.13, len(df))
for grp in STAGE_ORDER:
    m = (df["stage_group"] == grp).to_numpy()
    ax_b.scatter(
        df["n_celltypes"].to_numpy()[m] + jx[m],
        df["n_pas_clusters"].to_numpy()[m] + jy[m],
        s=64,
        facecolor=STAGE_COLOR[grp],
        edgecolor=SURFACE,
        linewidth=0.9,
        zorder=3,
        label=grp,
    )

# direct labels for the two most identifiable points
ib = int(df.index[df["gsm"] == "GSM3516664-MetBone"][0])
ax_b.annotate(
    "MetBone",
    xy=(df.at[ib, "n_celltypes"] + jx[ib], df.at[ib, "n_pas_clusters"] + jy[ib]),
    xytext=(6.5, 8.5),
    fontsize=6.2,
    color=STAGE_COLOR["Metastasis"],
    ha="left",
    va="center",
    arrowprops=dict(arrowstyle="-", color=STAGE_COLOR["Metastasis"], lw=0.7,
                    shrinkA=1, shrinkB=4),
)
iv = int(df.index[df["gsm"] == "GSM3516665-StageIVprimary"][0])
ax_b.annotate(
    "StageIVprimary",
    xy=(df.at[iv, "n_celltypes"] + jx[iv], df.at[iv, "n_pas_clusters"] + jy[iv]),
    xytext=(17.7, 10.4),
    fontsize=6.2,
    color=STAGE_COLOR["Stage IV primary"],
    ha="right",
    va="center",
    arrowprops=dict(arrowstyle="-", color=STAGE_COLOR["Stage IV primary"], lw=0.7,
                    shrinkA=1, shrinkB=4),
)

ax_b.set_xlim(lo, hi)
ax_b.set_ylim(lo, hi)
ax_b.set_xticks(np.arange(6, 19, 2))
ax_b.set_yticks(np.arange(6, 19, 2))
ax_b.set_aspect("equal", adjustable="box")
ax_b.set_xlabel("GEX cell types (count)")
ax_b.set_ylabel("PAS Leiden clusters (count)")
tidy_axes(ax_b, value_axis="y")
ax_b.grid(True, axis="x", color=GRID, linewidth=0.6)
ax_b.set_title(
    "b   PAS space finds structure of the\n"
    "      same granularity, slightly finer",
    loc="left",
    fontsize=8.5,
    color=INK,
    pad=6,
)
ax_b.text(
    0.035,
    0.965,
    f"above y = x: {n_pas_above}/17  on: {n_pas_equal}  below: {n_pas_below}\n"
    f"Spearman ρ = {rho_clusters.statistic:.2f} (p = {rho_clusters.pvalue:.2f})",
    transform=ax_b.transAxes,
    fontsize=6.0,
    color=INK,
    ha="left",
    va="top",
    linespacing=1.45,
)
ax_b.text(
    0.995,
    0.028,
    "counts jittered ±0.13 to separate ties;\n"
    "both counts on GEX-labelled cells only"
    + (
        ""
        if mean_npas_full is None
        else f"\n(conservative: full PAS clustering\n"
        f"averages {mean_npas_full:.1f}, not {df['n_pas_clusters'].mean():.1f} plotted)"
    ),
    transform=ax_b.transAxes,
    fontsize=5.4,
    color=MUTED,
    ha="right",
    va="bottom",
    linespacing=1.4,
)

# ---------------------------------------------------------------- panel (c) --
c = df.sort_values("coverage", ascending=True).reset_index(drop=True)
ypc = np.arange(len(c))
ax_c.barh(
    ypc,
    c["coverage"],
    height=0.62,
    color=[STAGE_COLOR[g] for g in c["stage_group"]],
    edgecolor=SURFACE,
    linewidth=0.7,
    zorder=2,
)
ax_c.axvline(med_cov, color=INK, lw=1.0, ls=(0, (4, 2)), zorder=3, alpha=0.75)

short = [g.replace("GSM351", "") for g in c["gsm"]]
ax_c.set_yticks(ypc)
ax_c.set_yticklabels(short, fontsize=5.0)
for tick, grp in zip(ax_c.get_yticklabels(), c["stage_group"]):
    tick.set_color(STAGE_COLOR[grp])
ax_c.set_ylim(-0.7, len(c) + 0.5)
ax_c.set_xlim(0, 1.0)
ax_c.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_c.set_xlabel("coverage = matched / PAS cells (fraction)")
tidy_axes(ax_c, value_axis="x")

for i, v in enumerate(c["coverage"]):
    # 3 dp, not 2: at 2 dp the 0.498 sample prints as "0.50", which contradicts
    # this panel's own "9 of 17 below 0.5" title (a reader counts only 8).
    ax_c.text(
        v + 0.015,
        i,
        f"{v:.3f}",
        fontsize=5.0,
        color=MUTED,
        ha="left",
        va="center",
        family="DejaVu Sans Mono",
    )
ax_c.text(
    med_cov + 0.014,
    len(c) + 0.05,
    f"median {med_cov:.3f}",
    fontsize=6.3,
    color=INK,
    ha="left",
    va="center",
)
ax_c.set_title(
    "c   Caveat: concordance is scored on a\n"
    f"      minority of PAS cells in {int((df['coverage'] < 0.5).sum())} of 17 samples",
    loc="left",
    fontsize=8.5,
    color=INK,
    pad=6,
)
ax_c.text(
    0.995,
    0.02,
    "matched = GEX-labelled cells;\n"
    "matched = gex_cells in all 17\n"
    "samples, so GEX limits coverage.\n"
    "Coverage does not drive\n"
    f"agreement (Spearman ρ = {rho_cov_ami.statistic:.2f},\n"
    f"p = {rho_cov_ami.pvalue:.2f}).",
    transform=ax_c.transAxes,
    fontsize=5.4,
    color=MUTED,
    ha="right",
    va="bottom",
    linespacing=1.55,
    bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.6),
)

# ------------------------------------------------------------ figure chrome --
fig.text(
    0.024,
    0.986,
    "Clustering on poly(A)-site profiles alone recovers gene-expression cell-type structure",
    fontsize=10,
    color=INK,
    ha="left",
    va="top",
    weight="bold",
)
fig.text(
    0.024,
    0.963,
    "17 Laughney lung-cancer samples · PAS clusters use no gene-level expression matrix · "
    "colour = clinical stage group",
    fontsize=7,
    color=MUTED,
    ha="left",
    va="top",
)
# provenance of the reference partition -- it is NOT curated ground truth, and
# the reader must know that before reading any agreement number on this page.
fig.text(
    0.024,
    0.9445,
    "Reference partition = marker-signature labels assigned per GEX Leiden cluster "
    "— an automated partition, not curated ground truth.",
    fontsize=6.3,
    color=MUTED,
    ha="left",
    va="top",
    style="italic",
)
# The feature space is per-site ABUNDANCE (TF-IDF + LSI), not within-gene usage
# fractions, so agreement cannot be attributed to isoform choice alone.
fig.text(
    0.024,
    0.9275,
    "Features are per-site counts (TF-IDF + LSI), which carry gene abundance "
    "— this figure does not isolate isoform-usage signal.",
    fontsize=6.3,
    color=MUTED,
    ha="left",
    va="top",
    style="italic",
)
fig.legend(
    handles=[
        Line2D([], [], marker="o", ls="", ms=6.5, mfc=STAGE_COLOR[g],
               mec=SURFACE, mew=0.9, label=f"{g} (n={int((df['stage_group']==g).sum())})")
        for g in STAGE_ORDER
    ],
    loc="upper left",
    bbox_to_anchor=(0.020, 0.9095),
    ncol=4,
    frameon=False,
    fontsize=7,
    handletextpad=0.35,
    columnspacing=1.3,
)

fig.savefig(PNG, dpi=300, facecolor=SURFACE)
save_manuscript(fig, NAME, facecolor=SURFACE)
plt.close(fig)

# ----------------------------------------------------------------------------
# tidy TSV: every number that appears on the page
# ----------------------------------------------------------------------------
rows = []


def add(panel, record, gsm, stage_label, stage_group, ypos_, metric, value):
    rows.append(
        dict(
            panel=panel,
            record=record,
            gsm=gsm,
            stage_label=stage_label,
            stage_group=stage_group,
            plot_position=ypos_,
            metric=metric,
            value=value,
        )
    )


for i, r in a.iterrows():
    for metric in ("AMI_celltype_vs_pas", "ARI_celltype_vs_pas"):
        add("a", "plotted", r["gsm"], r["stage_label"], r["stage_group"],
            int(i), metric, float(r[metric]))
for i, r in df.iterrows():
    for metric in ("n_celltypes", "n_pas_clusters"):
        add("b", "plotted", r["gsm"], r["stage_label"], r["stage_group"],
            "", metric, int(r[metric]))
for i, r in c.iterrows():
    add("c", "plotted", r["gsm"], r["stage_label"], r["stage_group"],
        int(i), "coverage", float(r["coverage"]))

# context columns (not drawn as marks but needed to audit the drawn values)
for _, r in df.iterrows():
    for metric in ("pas_cells", "gex_cells", "matched_cells", "n_gex_clusters",
                   "AMI_gexleiden_vs_pas", "ARI_gexleiden_vs_pas"):
        add("-", "source", r["gsm"], r["stage_label"], r["stage_group"],
            "", metric, float(r[metric]))

for panel, metric, value in [
    ("a", "median_AMI_celltype_vs_pas", med_ami),
    ("a", "median_ARI_celltype_vs_pas", med_ari),
    ("a", "min_AMI_celltype_vs_pas", float(df["AMI_celltype_vs_pas"].min())),
    ("a", "max_AMI_celltype_vs_pas", float(df["AMI_celltype_vs_pas"].max())),
    ("a", "min_ARI_celltype_vs_pas", float(df["ARI_celltype_vs_pas"].min())),
    ("a", "max_ARI_celltype_vs_pas", float(df["ARI_celltype_vs_pas"].max())),
    ("a", "kruskal_H_AMI_by_stage_3groups", float(kw_ami.statistic)),
    ("a", "kruskal_p_AMI_by_stage_3groups", float(kw_ami.pvalue)),
    ("-", "kruskal_H_ARI_by_stage_3groups", float(kw_ari.statistic)),
    ("-", "kruskal_p_ARI_by_stage_3groups", float(kw_ari.pvalue)),
    ("-", "mannwhitney_U_AMI_Normal_vs_Met", float(mw_norm_met.statistic)),
    ("-", "mannwhitney_p_AMI_Normal_vs_Met", float(mw_norm_met.pvalue)),
    ("-", "median_AMI_gexleiden_vs_pas", med_ami_leiden),
    ("-", "median_ARI_gexleiden_vs_pas", med_ari_leiden),
    ("b", "n_samples_pas_gt_celltypes", n_pas_above),
    ("b", "n_samples_pas_eq_celltypes", n_pas_equal),
    ("b", "n_samples_pas_lt_celltypes", n_pas_below),
    ("b", "mean_n_pas_clusters", float(df["n_pas_clusters"].mean())),
    ("b", "mean_n_celltypes", float(df["n_celltypes"].mean())),
    ("b", "median_n_pas_clusters", float(df["n_pas_clusters"].median())),
    ("b", "median_n_celltypes", float(df["n_celltypes"].median())),
    ("b", "mean_n_pas_clusters_full_clustering",
     float("nan") if mean_npas_full is None else mean_npas_full),
    ("b", "spearman_rho_ncelltypes_vs_npas", float(rho_clusters.statistic)),
    ("b", "spearman_p_ncelltypes_vs_npas", float(rho_clusters.pvalue)),
    ("c", "median_coverage", med_cov),
    ("c", "min_coverage", float(df["coverage"].min())),
    ("c", "max_coverage", float(df["coverage"].max())),
    ("c", "n_samples_coverage_below_0.5", int((df["coverage"] < 0.5).sum())),
    ("c", "spearman_rho_coverage_vs_AMI", float(rho_cov_ami.statistic)),
    ("c", "spearman_p_coverage_vs_AMI", float(rho_cov_ami.pvalue)),
    ("c", "spearman_rho_coverage_vs_ARI", float(rho_cov_ari.statistic)),
    ("c", "spearman_p_coverage_vs_ARI", float(rho_cov_ari.pvalue)),
    ("-", "max_abs_coverage_recompute_error", cov_dev),
    ("-", "matched_equals_gex_cells_all_samples", int(matched_eq_gex)),
    ("-", "summary_json_agrees_with_csv", int(json_agrees)),
]:
    add(panel, "summary", "", "", "", "", metric, value)

for grp in STAGE_ORDER:
    sub = df[df["stage_group"] == grp]
    add("a", "summary", "", "", grp, "", "n_samples", int(len(sub)))
    add("a", "summary", "", "", grp, "", "median_AMI_celltype_vs_pas",
        float(sub["AMI_celltype_vs_pas"].median()))
    add("a", "summary", "", "", grp, "", "median_ARI_celltype_vs_pas",
        float(sub["ARI_celltype_vs_pas"].median()))
    add("c", "summary", "", "", grp, "", "median_coverage",
        float(sub["coverage"].median()))

pd.DataFrame(rows).to_csv(TSV, sep="\t", index=False)

# ----------------------------------------------------------------------------
# console report
# ----------------------------------------------------------------------------
print(f"PNG  {PNG}")
print(f"TSV  {TSV}  ({len(rows)} rows)")
print(f"median AMI (celltype vs PAS) = {med_ami:.4f}  "
      f"range {df['AMI_celltype_vs_pas'].min():.4f}-{df['AMI_celltype_vs_pas'].max():.4f}")
print(f"median ARI (celltype vs PAS) = {med_ari:.4f}  "
      f"range {df['ARI_celltype_vs_pas'].min():.4f}-{df['ARI_celltype_vs_pas'].max():.4f}")
print(f"median AMI (gex leiden vs PAS) = {med_ami_leiden:.4f}; "
      f"median ARI = {med_ari_leiden:.4f}")
print(f"median coverage = {med_cov:.4f}  range "
      f"{df['coverage'].min():.4f}-{df['coverage'].max():.4f}; "
      f"{int((df['coverage']<0.5).sum())}/17 below 0.5")
print(f"n_pas>n_celltypes {n_pas_above}, == {n_pas_equal}, < {n_pas_below}; "
      f"mean n_pas {df['n_pas_clusters'].mean():.2f} vs mean n_ct "
      f"{df['n_celltypes'].mean():.2f}")
print(f"Kruskal-Wallis AMI by stage (3 groups, StageIVprimary n=1 excluded): "
      f"H={kw_ami.statistic:.3f} p={kw_ami.pvalue:.4f}")
print(f"Kruskal-Wallis ARI by stage (3 groups): H={kw_ari.statistic:.3f} "
      f"p={kw_ari.pvalue:.4f}")
print(f"Mann-Whitney AMI Normal vs Metastasis: U={mw_norm_met.statistic:.1f} "
      f"p={mw_norm_met.pvalue:.4f}")
print(f"Spearman coverage vs AMI: rho={rho_cov_ami.statistic:.3f} "
      f"p={rho_cov_ami.pvalue:.4f}")
print(f"Spearman coverage vs ARI: rho={rho_cov_ari.statistic:.3f} "
      f"p={rho_cov_ari.pvalue:.4f}")
print(f"Spearman n_celltypes vs n_pas: rho={rho_clusters.statistic:.3f} "
      f"p={rho_clusters.pvalue:.4f}")
print("per-stage medians (AMI / ARI / coverage):")
for grp in STAGE_ORDER:
    sub = df[df["stage_group"] == grp]
    print(f"  {grp:<18} n={len(sub)}  AMI={sub['AMI_celltype_vs_pas'].median():.4f}  "
          f"ARI={sub['ARI_celltype_vs_pas'].median():.4f}  "
          f"cov={sub['coverage'].median():.4f}")
print(f"integrity: max|coverage recompute error| = {cov_dev:.6f}; "
      f"matched==gex in all samples: {matched_eq_gex}; "
      f"summary.json agrees with csv: {json_agrees}")
