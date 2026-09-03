#!/usr/bin/env python3
"""
figS10_clustering.py -- manuscript Fig S10 ("figS10_clustering"): the DROPPED
clustering claim stated honestly -- rebuild of the retired
`clustering_concordance` WITH the gene-sum ablation as a mandatory co-equal
panel (21 D5 / 21 section 4.3-S7 / 05 section S10; executed 2026-09-02).

BINDING RULE (21 section 4.3-S7, 05 R2/R6): "A figure of the positive result
without the ablation on the same axes must never ship."  Panel a therefore
carries the headline positive number (PBMC AMI 0.708) ONLY beside the ablation
bar (per-gene totals, AMI 0.698) on one shared axis; panel b (the Laughney
cohort replication) exists downstream of that verdict and repeats it in its
title.  The claim ships as an OPTIONAL UTILITY for label-free exploration,
nothing more (01 section S4).

SOURCES OF TRUTH (nothing numeric typed by hand except EXPECT constants,
each cited to a verified document)
  panel a  results/figures/manuscript/pbmc_novelty.tsv, panel-d ablation rows
           (the rows 01 section S4 quotes as verified): "all PAS, counts
           (shipped)" AMI 0.708 @ n_features=275,370 vs "summed per gene =
           EXPRESSION" AMI 0.6978 @ n_features=14,891 (plus the k-matched-to-8
           secondary clustering of the same record); PBMC ARI 0.5023 from the
           panel-c rows.  The USAGE-fraction (isoform-only) rows of the same
           TSV are NOT plotted: computed on mis-keyed matrices, re-run pending
           (01 section S4; 10 section "Novelty / AMI ablations") -- the claim
           does not depend on their outcome.
  panel b  the original clustering_concordance input, verified FIXED in
           05 R2 (AMI/ARI recomputed from the raw per-cell labels of all 34
           .h5ad files, all nine headline numbers reproduced):
           /mnt/ssd2/.../RERUN_2026-08_fixed/runs/B1_cohort_full/
           B2_gex_celltyping/concordance.csv  (READ-ONLY)
  panel c  results/figures/manuscript/parameter_sweep.tsv (verified FIXED,
           05 R3) -- the retired sweep's ONE durable finding: Leiden
           resolution dominates cluster count (+112.5% median, 0.5 -> 2.0,
           17/17 datasets).  The caption carries the BINDING VERDICTS section 1
           label on the retired sweep (results/paramsweep/VERDICTS.md).

OUTPUTS
  manuscript/figures/figS10_clustering.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS10_clustering.caption.md  (sidecar, script-written --
      single source of the journal legend; the figure image carries no prose)
  results/figures/manuscript/figS10_clustering.tsv        (panels a+b values)
  results/figures/manuscript/figS10_clustering_sweep.tsv  (panel c values)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS10_clustering.py
"""
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TYPE, apply_rc, sentence_case   # shared publication style

apply_rc()   # DESIGN_DIRECTIVES.md item 1: one type scale across every figure
ANN = TYPE["annotation_min"]   # 6 pt floor for on-figure annotation


def sc_title(s):
    """Sentence-case a panel title (directive 6) while keeping the lower-case
    panel letter that prefixes it: 'a   the ...' -> 'a   The ...'."""
    m = re.match(r"^([a-z])(\s+)(.*)$", s, flags=re.S)
    return m.group(1) + m.group(2) + sentence_case(m.group(3)) if m else sentence_case(s)

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS10_clustering"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito colourblind-safe palette (figures/README.md convention); fixed
# entity assignment, never cycled.
C_SITES = "#0072B2"     # site-level PAS feature space / Laughney AMI
C_GENE = "#D55E00"      # per-gene-total (ablation) feature space / PBMC row
C_RES = "#009E73"       # Leiden-resolution knob (the one that matters)
C_OTHER = "#999999"     # every other swept knob

NOVELTY_TSV = OUTDIR / "pbmc_novelty.tsv"
SWEEP_TSV = OUTDIR / "parameter_sweep.tsv"
CONCORDANCE_CSV = Path(
    "/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/"
    "runs/B1_cohort_full/B2_gex_celltyping/concordance.csv")   # READ-ONLY

# ---------------------------------------------------------------------------
# load panel a: the verified ablation rows of pbmc_novelty.tsv (01 section S4)
# ---------------------------------------------------------------------------
nov = pd.read_csv(NOVELTY_TSV, sep="\t")


def nov_val(panel, record, key):
    q = nov[(nov.panel == panel) & (nov.record == record) & (nov.key == key)]
    assert len(q) == 1, (panel, record, key, len(q))
    return q.iloc[0]


ABL_SITES = "all PAS, counts  (shipped)"
ABL_GENE = "summed per gene = EXPRESSION"
ami_sites = float(nov_val("d", ABL_SITES, "AMI_vs_celltype_res1.0").value)
ami_gene = float(nov_val("d", ABL_GENE, "AMI_vs_celltype_res1.0").value)
ami_sites_k = float(nov_val("d", ABL_SITES, "AMI_vs_celltype_k_matched_to_8_types").value)
ami_gene_k = float(nov_val("d", ABL_GENE, "AMI_vs_celltype_k_matched_to_8_types").value)
n_sites = int(nov_val("d", ABL_SITES, "AMI_vs_celltype_res1.0").note.split("n_features=")[1].split(";")[0])
n_gene = int(nov_val("d", ABL_GENE, "AMI_vs_celltype_res1.0").note.split("n_features=")[1].split(";")[0])
pbmc_ami_c = float(nov_val("c", "pbmc_10k_v3", "AMI_pas_vs_celltype").value)
pbmc_ari_c = float(nov_val("c", "pbmc_10k_v3", "ARI_pas_vs_celltype").value)
lau_med_ami_c = float(nov_val("c", "laughney_n17", "median_AMI_celltype_vs_pas").value)
lau_med_ari_c = float(nov_val("c", "laughney_n17", "median_ARI_celltype_vs_pas").value)

# EXPECT -- verified headline values, 01 section S4 / 21 section 4.3-S7 / manifest S10 row
assert round(ami_sites, 3) == 0.708 and n_sites == 275370, (ami_sites, n_sites)
assert round(ami_gene, 4) == 0.6978 and n_gene == 14891, (ami_gene, n_gene)
assert round(pbmc_ami_c, 3) == 0.708 and round(pbmc_ari_c, 4) == 0.5023
assert abs(ami_sites - pbmc_ami_c) < 1e-9, "panel-c and panel-d site-space AMI must agree"
assert round(ami_sites_k, 4) == 0.7878 and round(ami_gene_k, 4) == 0.7602
d_abl = ami_sites - ami_gene           # what site resolution buys: 0.0102

# ---------------------------------------------------------------------------
# load panel b: Laughney concordance (verified FIXED, 05 R2)
# ---------------------------------------------------------------------------
lau = pd.read_csv(CONCORDANCE_CSV)
assert len(lau) == 17, f"expected 17 samples, got {len(lau)}"
med_ami = float(lau.AMI_celltype_vs_pas.median())
med_ari = float(lau.ARI_celltype_vs_pas.median())
min_ami, max_ami = float(lau.AMI_celltype_vs_pas.min()), float(lau.AMI_celltype_vs_pas.max())
min_ari, max_ari = float(lau.ARI_celltype_vs_pas.min()), float(lau.ARI_celltype_vs_pas.max())
med_cov = float(lau.coverage.median())
n_cov_lo = int((lau.coverage < 0.5).sum())
from scipy import stats
rho_cov = stats.spearmanr(lau.coverage, lau.AMI_celltype_vs_pas)

# EXPECT -- 05 R2 verified headlines, and agreement with pbmc_novelty.tsv panel c
assert round(med_ami, 4) == 0.6616 and round(med_ari, 4) == 0.4627, (med_ami, med_ari)
assert abs(med_ami - lau_med_ami_c) < 1e-4 and abs(med_ari - lau_med_ari_c) < 1e-4
assert abs(min_ami - 0.367) < 5e-4 and abs(max_ami - 0.780) < 6e-4, (min_ami, max_ami)
assert abs(min_ari - 0.200) < 5e-4 and abs(max_ari - 0.671) < 6e-4, (min_ari, max_ari)
assert (lau.AMI_celltype_vs_pas > lau.ARI_celltype_vs_pas).all(), "AMI > ARI must hold 17/17"
assert round(med_cov, 3) == 0.498 and n_cov_lo == 9, (med_cov, n_cov_lo)

# ---------------------------------------------------------------------------
# load panel c: the retired sweep's one durable finding (verified FIXED, 05 R3)
# ---------------------------------------------------------------------------
sw = pd.read_csv(SWEEP_TSV, sep="\t")
res_counts = sw[(sw.panel == "c") & (sw.branch.isin(["A3_res0.5", "A2_trim_default", "A3_res2.0"]))]
res_counts = res_counts.set_index("branch")
n05 = res_counts.loc["A3_res0.5"]
n10 = res_counts.loc["A2_trim_default"]     # res 1.0 (default)
n20 = res_counts.loc["A3_res2.0"]
assert (n05["median"], n10["median"], n20["median"]) == (11.0, 15.0, 23.0)

dn = sw[(sw.panel == "d") & (sw.metric == "n_clusters")].copy()
KNOBS = [  # (branch pair in TSV, display label)
    ("A3_res0.5 -> A3_res2.0", "Leiden res 0.5→2.0"),
    ("A3_nn15 -> A3_nn50", "n_neighbors 15→50"),
    ("A2_trim_default -> A3_libsize", "TF-IDF→libsize"),
    ("A2_trim_default -> A2_trim_d5000_ext", "incl. extended on"),
    ("A2_trim_d1000_ext -> A2_trim_d10000_ext", "gene dist. 1k→10k"),
    ("A2_trim_mult1.5 -> A2_trim_mult3.0", "utr mult. 1.5→3.0"),
]
dn = dn.set_index("branch")
res_row = dn.loc["A3_res0.5 -> A3_res2.0"]
# EXPECT -- the durable finding: +112.5% median, 17/17 datasets changed (05 R3)
assert res_row["median"] == 112.5 and res_row["n_datasets_changed"] == 17.0, dict(res_row)
libsize_row = dn.loc["A2_trim_default -> A3_libsize"]

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# canvas: the figure-level title/description chrome and the footer caption moved
# to the legend sidecar; height shrinks by the freed bands.  Width stays 8.3 in
# (211 mm): a 7.1 in (180 mm) attempt clipped panel c's y tick labels and panel
# b's title -- legibility outranks the double-column width target here.
fig = plt.figure(figsize=(8.3, 6.65))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.12], height_ratios=[1.0, 1.18],
                      left=0.152, right=0.965, top=0.935, bottom=0.083,
                      hspace=0.55, wspace=0.44)
axA = fig.add_subplot(gs[0, 0])
axC = fig.add_subplot(gs[1, 0])
axB = fig.add_subplot(gs[:, 1])


def style(ax, axis="x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    ax.grid(True, axis=axis, color=GRID, lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)


# ---- panel a: THE ABLATION -- positive result and its refutation, one axis --
axA.set_label("a")
rows_a = [
    (1.0, ami_sites, ami_sites_k, C_SITES,
     f"275,370 PAS sites\n(the published positive result)"),
    (0.0, ami_gene, ami_gene_k, C_GENE,
     f"14,891 per-gene totals\n(site resolution DISCARDED)"),
]
for y, v, vk, c, lab in rows_a:
    axA.barh(y, v, height=0.52, color=c, edgecolor="white", lw=0.8, zorder=3)
    axA.plot(vk, y, marker="D", ms=5, mfc="white", mec=c, mew=1.2, ls="none", zorder=5)
    axA.text(v - 0.012, y, f"{v:.3f}", ha="right", va="center", fontsize=7.2,
             color="white", fontweight="bold", zorder=4)
    axA.text(0.012, y + 0.38, lab, ha="left", va="bottom", fontsize=ANN + 0.3, color=c, zorder=5)
# Delta bracket
xb = max(ami_sites, ami_gene) + 0.022
axA.plot([xb, xb], [0, 1], color=INK, lw=0.9, zorder=4)
axA.plot([xb - 0.008, xb], [1, 1], color=INK, lw=0.9, zorder=4)
axA.plot([xb - 0.008, xb], [0, 0], color=INK, lw=0.9, zorder=4)
axA.text(xb + 0.015, 0.5, f"Δ AMI = {d_abl:+.3f}\nsite resolution\nadds nothing",
         ha="left", va="center", fontsize=6.4, color=INK, fontweight="bold")
axA.plot(pbmc_ari_c, 1.0, marker="o", ms=5, mfc="white", mec=C_SITES, mew=1.2,
         ls="none", zorder=6)
axA.annotate(f"ARI {pbmc_ari_c:.3f} (site space)", (pbmc_ari_c, 1.0),
             xytext=(0.535, 1.62), fontsize=ANN, color=MUTED,
             ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=3))
axA.set_xlim(0, 1.0)
axA.set_ylim(-0.75, 1.95)
axA.set_yticks([])
axA.set_xlabel(sentence_case("AMI vs GEX marker cell types (0 = chance)"))
axA.set_title(sc_title("a   the ablation that dropped the claim (PBMC)"), loc="left", fontweight="bold")
axA.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=7, mfc=MUTED, mec=MUTED, label="Leiden res 1.0 (bars)"),
    Line2D([], [], marker="D", ls="", ms=5, mfc="white", mec=MUTED, mew=1.2,
           label="k matched to 8 types"),
    Line2D([], [], marker="o", ls="", ms=5, mfc="white", mec=MUTED, mew=1.2, label="ARI"),
], loc="lower left", bbox_to_anchor=(0.0, -0.02), frameon=False, fontsize=ANN,
    handletextpad=0.35, labelspacing=0.3)
style(axA, axis="x")

# ---- panel b: the cohort replication of the (dropped) positive result -------
axB.set_label("b")
b = lau.sort_values("AMI_celltype_vs_pas", ascending=True).reset_index(drop=True)
for i, r in b.iterrows():
    axB.plot([r.ARI_celltype_vs_pas, r.AMI_celltype_vs_pas], [i, i], color=C_SITES,
             lw=1.5, solid_capstyle="round", alpha=0.75, zorder=2)
    axB.plot(r.ARI_celltype_vs_pas, i, "o", ms=4.6, mfc="white", mec=C_SITES, mew=1.3, zorder=3)
    axB.plot(r.AMI_celltype_vs_pas, i, "o", ms=5.0, mfc=C_SITES, mec="white", mew=0.8, zorder=4)
# PBMC as its own top row, distinct colour (different dataset, same recipe)
y_p = len(b) + 1.0
axB.plot([pbmc_ari_c, pbmc_ami_c], [y_p, y_p], color=C_GENE, lw=1.5,
         solid_capstyle="round", alpha=0.85, zorder=2)
axB.plot(pbmc_ari_c, y_p, "o", ms=4.6, mfc="white", mec=C_GENE, mew=1.3, zorder=3)
axB.plot(pbmc_ami_c, y_p, "o", ms=5.0, mfc=C_GENE, mec="white", mew=0.8, zorder=4)
axB.axvline(med_ami, color=INK, lw=0.9, ls=(0, (4, 2)), alpha=0.75, zorder=1)
axB.axvline(med_ari, color=MUTED, lw=0.9, ls=(0, (1, 2)), alpha=0.85, zorder=1)
axB.text(med_ami + 0.006, len(b) - 0.25, f"median AMI {med_ami:.3f}", fontsize=ANN,
         color=INK, ha="left", va="center")
axB.text(med_ari - 0.006, len(b) - 0.25, f"median ARI {med_ari:.3f}", fontsize=ANN,
         color=MUTED, ha="right", va="center")
axB.set_yticks(list(range(len(b))) + [y_p])
axB.set_yticklabels([g.replace("GSM351", "…") for g in b.gsm] + ["PBMC 10k v3"], fontsize=ANN)
for tick in axB.get_yticklabels():
    tick.set_color(MUTED)
axB.get_yticklabels()[-1].set_color(C_GENE)
axB.set_ylim(-0.8, y_p + 0.9)
axB.set_xlim(0.0, 0.9)
axB.set_xlabel(sentence_case("agreement with GEX marker cell types\n(0 = chance, 1 = identical)"))
# the '— but gene totals suffice (a)' corrective moved to the Legend sidecar (directive 2,
# no-loss: the legend opens with 'collapsing the sites to per-gene totals recovers them
# just as well' and (a) states the ablation numerically)
axB.set_title(sc_title("b   the (dropped) recovery replicates\n      in 17/17 Laughney samples"),
              loc="left", fontweight="bold")
axB.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=5.0, mfc=C_SITES, mec="white", mew=0.8, label="AMI (Laughney)"),
    Line2D([], [], marker="o", ls="", ms=4.6, mfc="white", mec=C_SITES, mew=1.3, label="ARI (Laughney)"),
    Line2D([], [], marker="o", ls="", ms=5.0, mfc=C_GENE, mec="white", mew=0.8, label="PBMC 10k v3"),
], loc="lower right", frameon=False, fontsize=ANN, handletextpad=0.35, labelspacing=0.3)
# the coverage/scoring caveat block moved to the legend (sidecar); a short data
# annotation keeps the 17/17 ordering fact visible
axB.text(0.02, 0.845, f"AMI range {min_ami:.3f}–{max_ami:.3f};\nAMI > ARI in 17/17",
         transform=axB.transAxes, fontsize=ANN, color=MUTED, ha="left", va="top",
         linespacing=1.5)
style(axB, axis="x")

# ---- panel c: the retired sweep's one durable finding -----------------------
axC.set_label("c")
yk = np.arange(len(KNOBS))[::-1]
for (pair, lab), y in zip(KNOBS, yk):
    r = dn.loc[pair]
    is_res = pair.startswith("A3_res")
    c = C_RES if is_res else C_OTHER
    axC.barh(y, r["median"], height=0.55, color=c, edgecolor="white", lw=0.8, zorder=3)
    axC.errorbar(r["median"], y, xerr=[[r["median"] - r["p25"]], [r["p75"] - r["median"]]],
                 color=INK, lw=0.9, capsize=2.0, zorder=4, ls="none")
    va = "center"
    if r["median"] >= 0:
        axC.text(max(r["p75"], r["median"]) + 4, y, f"{r['median']:+.1f}%",
                 ha="left", va=va, fontsize=ANN + 0.2,
                 color=C_RES if is_res else MUTED,
                 fontweight="bold" if is_res else "normal")
    else:
        axC.text(min(r["p25"], r["median"]) - 4, y, f"{r['median']:+.1f}%",
                 ha="right", va=va, fontsize=ANN + 0.2, color=MUTED)
axC.axvline(0, color=MUTED, lw=0.8, zorder=2)
axC.set_yticks(yk)
axC.set_yticklabels([sentence_case(lab) for _, lab in KNOBS], fontsize=TYPE["tick"])
# floor lowered so the cluster-count data note sits clear of the bottom bar's
# value label at the shared 6 pt floor (design pass 2026-09-03)
axC.set_ylim(-1.85, yk[0] + 0.62)
axC.set_xlim(-45, 165)
axC.set_xlabel(sentence_case("paired % change in Leiden cluster count per dataset\n"
                             "(median across 17 samples, whiskers = IQR)"))
# the "(retired sweep's sole durable finding)" aside moved to the Legend sidecar
# (directive 2, no-loss: it is the '(c) The retired parameter sweep's one durable
# finding: ...' sentence there)
axC.set_title(sc_title("c   one knob governs granularity"), loc="left", fontweight="bold")
# the library-size caveat moved to the legend (sidecar); the cluster counts stay
# as a data annotation
axC.text(0.97, 0.06,
         f"resolution 0.5 → 2.0: clusters {n05['median']:.0f} → {n20['median']:.0f}\n"
         f"(res 1.0 default: {n10['median']:.0f}); 17/17 datasets move,\n"
         f"per-dataset {res_row['q1_or_min']:+.0f}% to {res_row['q3_or_max']:+.0f}%",
         transform=axC.transAxes, fontsize=ANN, color=MUTED, ha="right", va="bottom",
         linespacing=1.45)
style(axC, axis="x")

# ---- legend (moved OFF the image; the .caption.md sidecar is its single source)
# The former on-figure title and description chrome open the legend verbatim.
legend = (
    "Figure S10 | PAS-profile clustering: the dropped claim, stated honestly (optional utility only). "
    "Cells cluster on poly(A)-site count profiles (TF-IDF + LSI + Leiden) and recover GEX cell types — "
    "but collapsing the sites to per-gene totals recovers them just as well: the signal is 3'-end gene "
    "expression re-encoded, not isoform choice. Offered only as an optional utility for label-free "
    "exploration when no GEX labels exist (01 S4). The reference partition is marker-signature labels "
    "assigned per GEX Leiden cluster — an automated partition, not curated ground truth; features are "
    "per-site counts (TF-IDF + LSI), which carry gene abundance. "
    f"(a) PBMC 10k v3: clustering on all {n_sites:,} PAS-site counts recovers GEX marker cell types "
    f"(AMI {ami_sites:.3f}, ARI {pbmc_ari_c:.3f}) — but summing the same counts to {n_gene:,} per-gene totals, "
    f"discarding all site-level resolution, recovers them just as well (AMI {ami_gene:.4f}; Δ {d_abl:+.4f}; "
    f"k-matched-to-8-types clusterings {ami_sites_k:.4f} vs {ami_gene_k:.4f}). The site resolution adds "
    "nothing, so the recovery is 3'-end expression signal re-encoded, not isoform choice (01 S4). "
    "The isoform-only USAGE-fraction spaces of the same experiment are NOT shown: they were computed on "
    "mis-keyed matrices and their re-run is pending; the ablation verdict does not depend on them (01 S4). "
    f"(b) The (dropped) positive result replicates: 17/17 Laughney samples, median AMI {med_ami:.3f} "
    f"(range {min_ami:.3f}-{max_ami:.3f}), median ARI {med_ari:.3f} (range {min_ari:.3f}-{max_ari:.3f}), "
    "AMI > ARI in 17/17; the reference partition is marker-derived, not curated truth, and concordance is "
    f"scored on GEX-labelled cells only (median coverage {med_cov:.3f}, {n_cov_lo}/17 below 0.5; coverage "
    f"does not drive agreement, Spearman rho {rho_cov.statistic:.2f}, p {rho_cov.pvalue:.2f}). "
    f"(c) The retired parameter sweep's one durable finding: Leiden resolution dominates cluster count "
    f"({res_row['median']:+.1f}% median, 0.5 to 2.0, 17/17 datasets; clusters {n05['median']:.0f} to "
    f"{n20['median']:.0f}) while every annotation-trim knob is at least 8x smaller. "
    f"The library-size bar's median must not be quoted alone (per-dataset range "
    f"{libsize_row['q1_or_min']:+.1f}% to {libsize_row['q3_or_max']:+.1f}%). "
    "BINDING LABEL (results/paramsweep/VERDICTS.md 1): --min-pas-prominence is PROVEN byte-identically "
    "inert on the lambda strategies, so no prominence arm of the retired sweep could ever have been "
    "informative; the sweep ran on the pre-clip-seeded lambda_gradient cohort run and only this "
    "clustering-knob finding survives."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=600)
print("wrote", p)

# ---------------------------------------------------------------------------
# TSVs: every plotted value, with a source column
# ---------------------------------------------------------------------------
SRC_NOV = "results/figures/manuscript/pbmc_novelty.tsv (verified ablation rows; 01 S4)"
SRC_LAU = str(CONCORDANCE_CSV) + " (verified FIXED, 05 R2)"
SRC_SW = "results/figures/manuscript/parameter_sweep.tsv (verified FIXED, 05 R3)"
rows = []
for space, n_feat, v, vk in [("all_PAS_site_counts_275370", n_sites, ami_sites, ami_sites_k),
                             ("per_gene_totals_14891", n_gene, ami_gene, ami_gene_k)]:
    rows.append(dict(panel="a", record=space, metric="AMI_vs_celltype_res1.0",
                     value=v, n_features=n_feat, source=SRC_NOV))
    rows.append(dict(panel="a", record=space, metric="AMI_vs_celltype_k_matched_to_8_types",
                     value=vk, n_features=n_feat, source=SRC_NOV))
rows.append(dict(panel="a", record="all_PAS_site_counts_275370", metric="ARI_vs_celltype",
                 value=pbmc_ari_c, n_features=n_sites, source=SRC_NOV))
rows.append(dict(panel="a", record="delta_sites_minus_gene_totals", metric="delta_AMI_res1.0",
                 value=d_abl, n_features="", source=SRC_NOV))
for _, r in lau.iterrows():
    rows.append(dict(panel="b", record=r.gsm, metric="AMI_celltype_vs_pas",
                     value=float(r.AMI_celltype_vs_pas), n_features="", source=SRC_LAU))
    rows.append(dict(panel="b", record=r.gsm, metric="ARI_celltype_vs_pas",
                     value=float(r.ARI_celltype_vs_pas), n_features="", source=SRC_LAU))
    rows.append(dict(panel="b", record=r.gsm, metric="coverage",
                     value=float(r.coverage), n_features="", source=SRC_LAU))
rows.append(dict(panel="b", record="pbmc_10k_v3", metric="AMI_pas_vs_celltype",
                 value=pbmc_ami_c, n_features="", source=SRC_NOV))
rows.append(dict(panel="b", record="pbmc_10k_v3", metric="ARI_pas_vs_celltype",
                 value=pbmc_ari_c, n_features="", source=SRC_NOV))
for metric, v in [("median_AMI", med_ami), ("median_ARI", med_ari),
                  ("min_AMI", min_ami), ("max_AMI", max_ami),
                  ("min_ARI", min_ari), ("max_ARI", max_ari),
                  ("median_coverage", med_cov), ("n_coverage_below_0.5", n_cov_lo),
                  ("spearman_rho_coverage_vs_AMI", float(rho_cov.statistic)),
                  ("spearman_p_coverage_vs_AMI", float(rho_cov.pvalue))]:
    rows.append(dict(panel="b", record="laughney_n17_summary", metric=metric,
                     value=v, n_features="", source=SRC_LAU))
p = OUTDIR / f"{NAME}.tsv"
pd.DataFrame(rows).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

srows = []
for branch, label in [("A3_res0.5", "res 0.5"), ("A2_trim_default", "res 1.0 (default)"),
                      ("A3_res2.0", "res 2.0")]:
    r = res_counts.loc[branch]
    srows.append(dict(panel="c", record=branch, x_label=label, metric="n_clusters_median",
                      value=float(r["median"]), q1=float(r.q1_or_min), q3=float(r.q3_or_max),
                      n_datasets=17, source=SRC_SW))
for pair, lab in KNOBS:
    r = dn.loc[pair]
    srows.append(dict(panel="c", record=pair, x_label=lab,
                      metric="paired_pct_change_n_clusters_median", value=float(r["median"]),
                      q1=float(r["p25"]), q3=float(r["p75"]),
                      n_datasets=int(r["n_datasets_changed"]), source=SRC_SW))
p = OUTDIR / f"{NAME}_sweep.tsv"
pd.DataFrame(srows).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S10 — `figS10_clustering` sidecar (generated by `scripts/manuscript_figures/figS10_clustering.py`; rebuild of retired `clustering_concordance` + ablation, executed 2026-09-02; single source of the legend)

## Legend

{legend}

**Binding rebuild rule honoured (21 §4.3-S7 / 05 §S10):** the positive result never appears without the
ablation on the same axes — panel a carries the headline PBMC number only beside the per-gene-total bar on
one shared AMI axis, and panels b/c are downstream of that verdict.

**Caveats that travel with this figure (05 R2/R3, verbatim obligations):**
- The reference partition is marker-signature labels assigned per GEX Leiden cluster — an automated
  partition, not curated ground truth.
- Features are per-site counts (TF-IDF + LSI); per-site counts carry gene abundance, and this figure's own
  panel a shows the abundance signal is sufficient — the figure does not demonstrate isoform-choice signal.
- Concordance is scored on GEX-labelled cells only (median coverage {med_cov:.3f}; {n_cov_lo}/17 samples
  below 0.5); coverage does not drive agreement (Spearman ρ = {rho_cov.statistic:.2f}, p = {rho_cov.pvalue:.2f}).
- No stage effect is claimed for panel b (the retired figure's Kruskal–Wallis was null; n = 4/7/5/1, low power).
- The isoform-only USAGE-fraction spaces of the pbmc_novelty experiment are excluded: computed on mis-keyed
  matrices, re-run pending (01 §S4); the ablation verdict does not depend on their outcome.
- Panel c binding label (`results/paramsweep/VERDICTS.md` §1): `--min-pas-prominence` is PROVEN
  byte-identically inert on the lambda strategies (`compute_lambda(heights)` overrides it), so any retired
  `parameter_sweep` arm sweeping it was a no-op; the retired sweep ran on the pre-clip-seeded
  `lambda_gradient` cohort run — the wrong caller for every current claim — and only the clustering-knob
  finding shown here survives. The library-size bar's median must not be quoted alone (per-dataset range
  {libsize_row['q1_or_min']:+.1f}% to {libsize_row['q3_or_max']:+.1f}%). Annotation-knob neutrality is
  median-neutrality: 7–9 of 17 datasets shift by 1–2 clusters.
- Cluster count is a granularity proxy only; no cross-branch ARI was computed, so nothing here says cell
  *assignments* are stable across knobs.

## Provenance

Sources: `{SRC_NOV}`; `{SRC_LAU}`; `{SRC_SW}`; claim framing `manuscript/01_outline_and_journals.md` §S4.
Every plotted value: `results/figures/manuscript/figS10_clustering.tsv` and `figS10_clustering_sweep.tsv`.
Script: `scripts/manuscript_figures/figS10_clustering.py`; rendered at PNG 600 dpi / PDF fonttype 42,
8.3 in (211 mm) wide (a 7.1 in narrowing clipped panel b/c labels and was rejected).
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure annotation,
key entry and tick label now sits at or above the 6 pt floor. Axis labels, panel titles and prose tick labels
are sentence-cased through `_pubstyle.sentence_case()`, canonical identifiers preserved and the lower-case panel
letters kept. Two title clauses left the image for the Legend above: panel b's *— but gene totals suffice (a)*
(no-loss: the legend's opening sentence and panel a state it) and panel c's *(retired sweep's sole durable
finding)* (no-loss: it is the legend's \"(c) The retired parameter sweep's one durable finding\" sentence). Two
layout fixes: panel c's floor was lowered so its cluster-count note clears the bottom bar's value label, and the
left margin was widened so panel c's knob labels clear the canvas edge. No panel, number or audit TSV changed;
both TSVs regenerate byte-identical.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS10_clustering: sites AMI {ami_sites:.4f} (n={n_sites:,}) vs gene totals "
      f"{ami_gene:.4f} (n={n_gene:,}), delta {d_abl:+.4f}; Laughney median AMI {med_ami:.4f} / "
      f"ARI {med_ari:.4f} (17/17); Leiden res +{res_row['median']:.1f}%")
