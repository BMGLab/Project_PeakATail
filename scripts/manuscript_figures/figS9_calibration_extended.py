#!/usr/bin/env python3
"""
figS9_calibration_extended.py -- manuscript Fig S9 ("figS9_calibration_extended"):
the calibration deep-dive behind Fig 4, entirely over the existing verified
outputs of the fdr_calibration_v2 run.  NOTHING is re-run: no permutation, no
switch test -- every panel is read (or deterministically aggregated) from
files already behind manuscript/14 (verified SOUND 2026-08-21).

Per FIGURES_MANIFEST.md row S9 and manuscript/21 sec.4.3 row S8 (plan):
per-stage-pair and per-expression-stratum null rates (2.6-4.0% for the
calibrated arm), the permutation-calibrated q comparison, reads-vs-cells
count mode, the NB dispersion-floor diagnostic, and (21) the runtime of the
permutation harness.  Cited with Fig 4.

SOURCES OF TRUTH (verified)
  * results/figures/manuscript/fig4_calibration_stats.tsv     headline stats/arm
  * results/figures/manuscript/fig4_calibration_per_pair.tsv  per-stage-pair
  * results/fdr_calibration_v2/fig4_calibration_null_pvalues_all_arms.tsv.gz
        every null p-value of every arm (the file the verifier reproduced 14
        from) -- used ONLY for the B0 expression-stratum aggregation
  * results/fdr_calibration_v2/input/stage_labelled.h5ad      the run's input
        matrix (1,230 cells x 76,711 PAS) -- expressing-cells per PAS
  * manuscript/14_switch_calibration_v2.md                    EXPECT values

STRATUM CAVEAT (stated in the legend of the .caption.md sidecar)
  14 records "6 strata, <50 to >=600 expressing cells gives 2.6-4.0%
  everywhere" but the per-stratum values and exact bin edges were never
  persisted.  Panel b is a deterministic recomputation from the two archived
  files above under the declared method (expressing cells = cells with
  count > 0 among all 1,230 labelled cells, label-independent; edges
  50/100/200/400/600).  It gives 2.5-4.0%; 14's conclusion -- the calibrated
  arm is never anti-conservative at any expression level -- reproduces
  exactly, and the doc's quoted range is drawn as a reference band.  No other
  panel involves any recomputation.

PANELS
  a  Per-stage-pair null p<0.05 of all six arms (no pair drives any verdict;
     B0 2.94-3.13% everywhere, B 12.6-13.7% everywhere).
  b  B0 null p<0.05 by expression stratum (recomputed; see caveat).
  c  reads-vs-cells count mode at both marker settings (the pseudoreplication
     mechanism, 14 sec. mechanism 3): 20.3->13.0% with pre-selection on,
     9.6->3.0% off; mean false q<0.05 hits per null run annotated.
  d  Permutation-calibrated q on the TRUE runs: nominal BH hits vs q_perm<0.05
     vs q_perm + effect floor, per arm (log axis).  For B0 the permutation
     null is conservative, so q_perm ADDS hits (66,630 vs 60,332).
  e  NB dispersion-floor diagnostic: tests at the 1e-4 plug-in floor are 8.5%
     of C0's null tests but carry 67% of its null q<0.05 hits (13% in C).

OUTPUTS
  manuscript/figures/figS9_calibration_extended.{png,pdf,caption.md}
      (PNG 600 dpi; PDF fonttype 42; legend lives in the sidecar, not on the image)
  results/figures/manuscript/figS9_calibration_extended_per_pair.tsv   (a)
  results/figures/manuscript/figS9_calibration_extended_strata.tsv    (b)
  results/figures/manuscript/figS9_calibration_extended_countmode.tsv (c)
  results/figures/manuscript/figS9_calibration_extended_permq.tsv     (d)
  results/figures/manuscript/figS9_calibration_extended_dispfloor.tsv (e)
  results/figures/manuscript/figS9_calibration_extended_runtime.tsv   (caption)

Env: FIGS9_REUSE_STRATA=1 reuses an existing strata TSV instead of re-reading
the 4.3M-row gz + h5ad (layout iteration only; asserts still run either way).

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS9_calibration_extended.py
"""
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TYPE, apply_rc, plain_log, sentence_case   # shared publication style

apply_rc()   # DESIGN_DIRECTIVES.md item 1: one type scale across every figure
ANN = TYPE["annotation_min"]   # 6 pt floor for on-figure annotation


def sc_title(s):
    """Sentence-case a panel title (directive 6) while keeping the lower-case
    panel letter that prefixes it: 'a   the ...' -> 'a   The ...'."""
    m = re.match(r"^([a-z])(\s+)(.*)$", s, flags=re.S)
    return m.group(1) + m.group(2) + sentence_case(m.group(3)) if m else sentence_case(s)

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
TSVDIR = WD / "results/figures/manuscript"
V2 = WD / "results/fdr_calibration_v2"
FIGDIR = WD / "manuscript/figures"
NAME = "figS9_calibration_extended"

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito; colour follows the TEST (fig4 convention): fisher-reads BLUE,
# fisher-cells GREEN, nb_pairwise VERMILION; marker-off = hollow.
BLUE, GREEN, VERM = "#0072B2", "#009E73", "#D55E00"
ARMS = {
    "A_fisher_reads": dict(label="A  fisher, count-mode reads, top-200 markers", short="fisher\nreads", color=BLUE, markers=True),
    "B_fisher_cells": dict(label="B  fisher, count-mode cells, top-200 markers", short="fisher\ncells", color=GREEN, markers=True),
    "C_nb_pairwise": dict(label="C  nb_pairwise, top-200 markers", short="nb\npairwise", color=VERM, markers=True),
    "A0_fisher_reads_nomarker": dict(label="A0  fisher reads, no marker pre-sel.", short="fisher\nreads\nno marker", color=BLUE, markers=False),
    "B0_fisher_cells_nomarker": dict(label="B0  fisher cells, no marker pre-sel. (SHIPS)", short="fisher\ncells\nno marker", color=GREEN, markers=False),
    "C0_nb_pairwise_nomarker": dict(label="C0  nb_pairwise, no marker pre-sel.", short="nb\npairwise\nno marker", color=VERM, markers=False),
}
# 'short' is used only for panel a's x tick labels; the no-marker arms were re-wrapped
# from two lines to three in the 2026-09-03 design pass so that adjacent labels clear
# each other at the shared 6 pt floor (same words, no content dropped).
ORDER = list(ARMS)
PAIRS = ["ES_vs_RS", "ES_vs_SPC", "RS_vs_SPC"]
PAIR_MARK = {"ES_vs_RS": "o", "ES_vs_SPC": "s", "RS_vs_SPC": "^"}

# ---------------------------------------------------------------------------
# verified tables
# ---------------------------------------------------------------------------
stats = pd.read_csv(TSVDIR / "fig4_calibration_stats.tsv", sep="\t").set_index("arm")
per_pair = pd.read_csv(TSVDIR / "fig4_calibration_per_pair.tsv", sep="\t")
SRC_STATS = "results/figures/manuscript/fig4_calibration_stats.tsv"
SRC_PP = "results/figures/manuscript/fig4_calibration_per_pair.tsv"

def pct1(x):
    return round(float(x) * 1000) / 10

# EXPECT: the 14 sec."Arms and verdicts" table (verified SOUND)
assert pct1(stats.loc["A_fisher_reads", "frac_null_p_lt_05"]) == 20.3
assert pct1(stats.loc["B_fisher_cells", "frac_null_p_lt_05"]) == 13.0
assert pct1(stats.loc["C_nb_pairwise", "frac_null_p_lt_05"]) == 24.7
assert pct1(stats.loc["A0_fisher_reads_nomarker", "frac_null_p_lt_05"]) == 9.6
assert pct1(stats.loc["B0_fisher_cells_nomarker", "frac_null_p_lt_05"]) == 3.0
assert pct1(stats.loc["C0_nb_pairwise_nomarker", "frac_null_p_lt_05"]) == 5.1
assert stats.loc["B0_fisher_cells_nomarker", "frac_null_q_lt_fdr"] == 0.0
assert int(stats.loc["B0_fisher_cells_nomarker", "n_runs_with_q_hit"]) == 0
assert bool(stats.loc["B0_fisher_cells_nomarker", "calibrated_by_prereg_rule"]) is True
assert bool(stats.loc["B0_fisher_cells_nomarker", "conservative"]) is True
assert round(float(stats.loc["A0_fisher_reads_nomarker", "mean_hits_per_null_run"])) == 1487
assert round(float(stats.loc["C0_nb_pairwise_nomarker", "mean_hits_per_null_run"])) == 139

# EXPECT: 14 "Per-pair breakdown": B0 2.94-3.13%, B 12.6-13.7%, no pair excepted
b0p = per_pair[per_pair.arm == "B0_fisher_cells_nomarker"]
bp = per_pair[per_pair.arm == "B_fisher_cells"]
assert len(b0p) == 3 and (b0p.null_frac_q_lt_fdr == 0).all()
assert (round(b0p.null_frac_p_lt_05.min() * 10000) / 100,
        round(b0p.null_frac_p_lt_05.max() * 10000) / 100) == (2.94, 3.13)
assert (pct1(bp.null_frac_p_lt_05.min()), pct1(bp.null_frac_p_lt_05.max())) == (12.6, 13.7)

# EXPECT: 14 sec."Mechanism" 4 -- NB dispersion floor (C0 8.5% of null tests,
# 67% of null q hits; 13% in C) and min null p 4e-182
assert pct1(stats.loc["C0_nb_pairwise_nomarker", "null_share_dispersion_floor"]) == 8.5
assert round(float(stats.loc["C0_nb_pairwise_nomarker", "null_hits_share_dispersion_floor"]) * 100) == 67
assert round(float(stats.loc["C_nb_pairwise", "null_hits_share_dispersion_floor"]) * 100) == 13
assert f"{float(stats.loc['C0_nb_pairwise_nomarker', 'min_null_p']):.0e}" == "4e-182"

# EXPECT: 14 sec."Permutation-calibrated q": A 186/307 (85 with floor),
# B 640/668 (572), C 1,733/1,754 (1,588), B0 66,630 vs 60,332 (46,230);
# resolution ~0.00046 for the marker arms
for arm, nom, permq, floor in [("A_fisher_reads", 307, 186, 85), ("B_fisher_cells", 668, 640, 572),
                               ("C_nb_pairwise", 1754, 1733, 1588),
                               ("B0_fisher_cells_nomarker", 60332, 66630, 46230)]:
    s = stats.loc[arm]
    assert (int(s.true_q_hits), int(s.true_perm_q_hits), int(s.true_perm_q_hits_with_effect_floor)) \
        == (nom, permq, floor), arm
assert round(float(stats.loc["A_fisher_reads", "min_p_emp"]), 5) == 0.00046
assert f"{float(stats.loc['B0_fisher_cells_nomarker', 'min_p_emp']):.0e}" == "7e-07"

# ---------------------------------------------------------------------------
# panel b: B0 null p<0.05 by expression stratum (deterministic recomputation)
# ---------------------------------------------------------------------------
STRATA_TSV = TSVDIR / f"{NAME}_strata.tsv"
EDGES = [50, 100, 200, 400, 600]
SLAB = ["<50", "50-99", "100-199", "200-399", "400-599", ">=600"]
STRAT_METHOD = ("expressing cells per PAS = cells with count>0 in "
                "results/fdr_calibration_v2/input/stage_labelled.h5ad layers['counts'] "
                "(all 1,230 labelled cells, label-independent); tests = every "
                "B0_fisher_cells_nomarker row of fig4_calibration_null_pvalues_all_arms.tsv.gz "
                "(20 permutations x 3 pairs); bin edges 50/100/200/400/600")
if os.environ.get("FIGS9_REUSE_STRATA") == "1" and STRATA_TSV.exists():
    strata = pd.read_csv(STRATA_TSV, sep="\t")
else:
    import warnings
    warnings.filterwarnings("ignore")
    import anndata as ad
    A = ad.read_h5ad(V2 / "input/stage_labelled.h5ad")
    assert A.shape == (1230, 76711)
    assert int(A.layers["counts"].sum()) == 32702229   # report.json counts_layer_sum
    nz = np.asarray((A.layers["counts"] > 0).sum(axis=0)).ravel()
    m = dict(zip(A.var.pas_id.astype(int), nz))
    it = pd.read_csv(V2 / "fig4_calibration_null_pvalues_all_arms.tsv.gz", sep="\t",
                     usecols=["arm", "pas_id", "pvalue"], chunksize=2_000_000)
    b0 = pd.concat([c[c.arm == "B0_fisher_cells_nomarker"][["pas_id", "pvalue"]] for c in it],
                   ignore_index=True)
    assert len(b0) == int(stats.loc["B0_fisher_cells_nomarker", "n_null_tests"])
    e = b0.pas_id.map(m).to_numpy()
    assert not np.isnan(e.astype(float)).any()
    bins = np.digitize(e, EDGES)
    rows = []
    for i, lab in enumerate(SLAB):
        sel = bins == i
        rows.append(dict(panel="b", arm="B0_fisher_cells_nomarker",
                         stratum=lab, n_null_tests=int(sel.sum()),
                         n_null_p_lt_05=int((b0.pvalue[sel] < 0.05).sum()),
                         frac_null_p_lt_05=float((b0.pvalue[sel] < 0.05).mean()),
                         method=STRAT_METHOD,
                         doc_range_manuscript14="2.6-4.0% (6 strata, <50 to >=600 expressing cells; "
                                                "verifier's exact bin edges not persisted)",
                         source=("results/fdr_calibration_v2/fig4_calibration_null_pvalues_all_arms.tsv.gz"
                                 " + results/fdr_calibration_v2/input/stage_labelled.h5ad")))
    strata = pd.DataFrame(rows)
# EXPECT gates on the recomputation:
#  * pooled rate must equal the verified per-arm stat (same rows, so exactly)
pooled = strata.n_null_p_lt_05.sum() / strata.n_null_tests.sum()
assert abs(pooled - float(stats.loc["B0_fisher_cells_nomarker", "frac_null_p_lt_05"])) < 1e-9
#  * every stratum below nominal 5% and far below the 7% pre-registered gate
assert (strata.frac_null_p_lt_05 < 0.05).all()
#  * the doc's upper end reproduces exactly (4.0%); the min lands at 2.5% under
#    the declared binning vs the doc's 2.6% -- consistent to 0.1 pp (caveated)
assert pct1(strata.frac_null_p_lt_05.max()) == 4.0
assert 2.4 <= pct1(strata.frac_null_p_lt_05.min()) <= 2.6
assert strata.n_null_tests.sum() == 4281158

# ---------------------------------------------------------------------------
# remaining audit tables
# ---------------------------------------------------------------------------
pp_rows = []
for _, r in per_pair.iterrows():
    pp_rows.append(dict(panel="a", arm=r.arm, arm_label=ARMS[r.arm]["label"], pair=r.pair,
                        n_null_tests=int(r.n_null_tests),
                        null_frac_p_lt_05=float(r.null_frac_p_lt_05),
                        null_frac_q_lt_fdr=float(r.null_frac_q_lt_fdr),
                        null_mean_hits_per_run=float(r.null_mean_hits_per_run),
                        source=SRC_PP))
cm_rows = []
for arm, presel in [("A_fisher_reads", "top-200 wilcoxon"), ("B_fisher_cells", "top-200 wilcoxon"),
                    ("A0_fisher_reads_nomarker", "off"), ("B0_fisher_cells_nomarker", "off")]:
    s = stats.loc[arm]
    cm_rows.append(dict(panel="c", arm=arm, count_mode=("reads" if "reads" in arm else "cells"),
                        marker_preselection=presel,
                        frac_null_p_lt_05=float(s.frac_null_p_lt_05),
                        mean_false_q_hits_per_null_run=float(s.mean_hits_per_null_run),
                        n_runs_with_q_hit=int(s.n_runs_with_q_hit), source=SRC_STATS))
pq_rows = []
for arm in ORDER:
    s = stats.loc[arm]
    pq_rows.append(dict(panel="d", arm=arm, arm_label=ARMS[arm]["label"],
                        true_nominal_q_hits=int(s.true_q_hits),
                        true_perm_q_hits=int(s.true_perm_q_hits),
                        true_perm_q_hits_with_effect_floor=int(s.true_perm_q_hits_with_effect_floor),
                        effect_kind=s.effect_kind, effect_floor=float(s.effect_floor),
                        perm_q_resolution_min_p_emp=float(s.min_p_emp), source=SRC_STATS))
df_rows = []
for arm in ["C_nb_pairwise", "C0_nb_pairwise_nomarker"]:
    s = stats.loc[arm]
    df_rows.append(dict(panel="e", arm=arm, arm_label=ARMS[arm]["label"],
                        null_share_tests_at_dispersion_floor=float(s.null_share_dispersion_floor),
                        null_share_q_hits_from_floored_tests=float(s.null_hits_share_dispersion_floor),
                        dispersion_floor=1e-4, min_null_p=float(s.min_null_p), source=SRC_STATS))
rt_rows = [dict(panel="caption", arm=a, true_run_elapsed_s=float(stats.loc[a, "true_elapsed_s"]),
                perm_run_elapsed_mean_s=float(stats.loc[a, "perm_elapsed_mean_s"]),
                n_perms=int(stats.loc[a, "n_perms"]), source=SRC_STATS) for a in ORDER]
rt = pd.DataFrame(rt_rows)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# canvas: the old footer band (bottom=0.235 of a 10.6 in canvas) is gone -- the
# legend lives in the sidecar; height shrinks by the freed space, axes keep size
# panel d grew (extra row height + wider bar spacing below): at the shared 6 pt
# annotation floor its three per-arm count labels collided (design pass 2026-09-03)
fig = plt.figure(figsize=(8.3, 9.0))
gs = fig.add_gridspec(3, 2, left=0.138, right=0.965, top=0.958, bottom=0.064,
                      height_ratios=[1.0, 1.0, 1.2], hspace=0.56, wspace=0.30)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axE = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[2, :])


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    ax.grid(True, color=GRID, lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def rate_lines(ax, x0=0.0):
    ax.axhline(5.0, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=2)
    ax.axhline(7.0, color=VERM, lw=0.8, ls=(0, (1, 1.5)), zorder=2)


# ---- panel a: per-stage-pair null rates ----------------------------------
for xi, arm in enumerate(ORDER):
    d = ARMS[arm]
    for k, pair in enumerate(PAIRS):
        r = per_pair[(per_pair.arm == arm) & (per_pair.pair == pair)]
        assert len(r) == 1
        v = float(r.null_frac_p_lt_05.iloc[0]) * 100
        kw = dict(color=d["color"]) if d["markers"] else dict(mfc="none", mec=d["color"], mew=0.9)
        axA.plot([xi + (k - 1) * 0.20], [v], PAIR_MARK[pair], ms=4.2, zorder=4, **kw)
rate_lines(axA)
axA.text(-0.45, 5.35, "5% nominal", fontsize=ANN, color=INK, ha="left")
axA.text(-0.45, 7.35, "7% pre-registered gate", fontsize=ANN, color=VERM, ha="left")
axA.annotate("B0: 2.94-3.13%\nin every pair,\n0/20 runs with a\nq<0.05 hit", xy=(4.0, 3.4),
             xytext=(3.62, 12.0), fontsize=ANN, color=INK, ha="center", linespacing=1.3,
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
axA.set_xticks(range(6))
axA.set_xticklabels([ARMS[a]["short"] for a in ORDER], fontsize=ANN)
axA.set_ylabel(sentence_case("null p<0.05 per stage pair (%)"))
axA.set_ylim(0, 27.5)
axA.set_title(sc_title("a   no stage pair drives any verdict"), loc="left", fontweight="bold")
axA.legend([Line2D([], [], marker=PAIR_MARK[p], ls="none", ms=4.2, color=MUTED) for p in PAIRS]
           + [Line2D([], [], marker="o", ls="none", ms=4.2, color=MUTED),
              Line2D([], [], marker="o", ls="none", ms=4.2, mfc="none", mec=MUTED, mew=0.9)],
           [p.replace("_vs_", " vs ") for p in PAIRS] + ["marker pre-sel. on", "marker pre-sel. off"],
           loc="upper right", frameon=False, fontsize=ANN, handlelength=1.0, labelspacing=0.25)
style(axA)

# ---- panel b: B0 expression strata ---------------------------------------
xs = np.arange(len(SLAB))
axB.axhspan(2.6, 4.0, color=GRID, alpha=0.55, zorder=1)
axB.bar(xs, strata.frac_null_p_lt_05 * 100, width=0.62, facecolor="none",
        edgecolor=GREEN, lw=1.2, zorder=3)
for x, v in zip(xs, strata.frac_null_p_lt_05 * 100):
    axB.text(x, v + 0.12, f"{v:.1f}", ha="center", fontsize=ANN, color=INK)
rate_lines(axB)
axB.text(5.35, 5.15, "5% nominal", fontsize=ANN, color=INK, ha="right")
axB.text(5.35, 7.15, "7% pre-registered gate", fontsize=ANN, color=VERM, ha="right")
# the recomputation caveat behind the band moved to the legend (sidecar); the
# band keeps a short name on the image
axB.text(0.02, 0.99, "band: verified range 2.6-4.0%", transform=axB.transAxes,
         fontsize=ANN, color=MUTED, va="top")
axB.set_xticks(xs)
axB.set_xticklabels([f"{l}\n{n/1e6:.2f}M" if n >= 1e6 else f"{l}\n{n/1000:.0f}k"
                     for l, n in zip(SLAB, strata.n_null_tests)], fontsize=ANN)
axB.set_xlabel(sentence_case("expressing cells per PAS (bin; n null tests)"))
axB.set_ylabel(sentence_case("B0 null p<0.05 (%)"))
axB.set_ylim(0, 8.7)
axB.set_title(sc_title("b   calibrated arm, by expression stratum"), loc="left", fontweight="bold")
style(axB)

# ---- panel c: reads vs cells count mode ----------------------------------
groups = [("top-200 markers\n(shipped default)", "A_fisher_reads", "B_fisher_cells"),
          ("no marker pre-selection", "A0_fisher_reads_nomarker", "B0_fisher_cells_nomarker")]
for gi, (glab, a_reads, a_cells) in enumerate(groups):
    for k, (arm, col) in enumerate([(a_reads, BLUE), (a_cells, GREEN)]):
        s = stats.loc[arm]
        v = float(s.frac_null_p_lt_05) * 100
        filled = ARMS[arm]["markers"]
        kw = dict(color=col) if filled else dict(facecolor="none", edgecolor=col, lw=1.2)
        x = gi * 1.15 + (k - 0.5) * 0.42
        axC.bar(x, v, width=0.38, zorder=3, **kw)
        axC.text(x, v + 0.4, f"{v:.1f}%", ha="center", fontsize=ANN, color=INK)
        mh = float(s.mean_hits_per_null_run)
        mh_lab = "0" if mh == 0 else (f"{mh:,.0f}" if mh >= 100 else f"{mh:.1f}")
        axC.text(x, 0.55, mh_lab, ha="center", fontsize=ANN,
                 color=("white" if filled else col), zorder=4, fontweight="bold")
rate_lines(axC)
axC.text(1.72, 5.28, "5% nominal", fontsize=ANN, color=INK, ha="right")
axC.text(1.72, 7.28, "7% gate", fontsize=ANN, color=VERM, ha="right")
axC.set_xticks([0, 1.15])
axC.set_xticklabels([sentence_case(g[0]) for g in groups], fontsize=TYPE["tick"])
axC.set_ylabel(sentence_case("null p<0.05 (%)"))
axC.set_ylim(0, 23.5)
axC.set_title(sc_title("c   count-mode reads vs cells"), loc="left", fontweight="bold")
axC.legend([Patch(color=BLUE), Patch(color=GREEN)],
           ["count-mode reads (pseudoreplicates UMIs within cells)", "count-mode cells"],
           loc="upper right", frameon=False, fontsize=ANN, handlelength=1.1, labelspacing=0.3)
axC.text(0.99, 0.74, "numbers at bar bases:\nmean false q<0.05 hits\nper null run",
         transform=axC.transAxes, fontsize=ANN, color=MUTED, ha="right", va="top", linespacing=1.3)
style(axC)
axC.grid(False, axis="x")

# ---- panel e: NB dispersion floor ----------------------------------------
labels_e = ["share of null tests\nat the 1e-4 floor", "share of null q<0.05 hits\nfrom floored tests"]
for gi, arm in enumerate(["C_nb_pairwise", "C0_nb_pairwise_nomarker"]):
    s = stats.loc[arm]
    filled = ARMS[arm]["markers"]
    for k, v in enumerate([float(s.null_share_dispersion_floor) * 100,
                           float(s.null_hits_share_dispersion_floor) * 100]):
        kw = dict(color=VERM) if filled else dict(facecolor="none", edgecolor=VERM, lw=1.2)
        x = k * 1.1 + (gi - 0.5) * 0.42
        axE.bar(x, v, width=0.38, zorder=3, **kw)
        axE.text(x, v + 1.3, f"{v:.1f}%", ha="center", fontsize=ANN, color=INK)
axE.set_xticks([0, 1.1])
axE.set_xticklabels([sentence_case(s) for s in labels_e], fontsize=TYPE["tick"])
axE.set_ylabel(sentence_case("% of C-arm null tests / hits"))
axE.set_ylim(0, 84)
axE.set_title(sc_title("e   NB plug-in dispersion floor inflates the tail"), loc="left", fontweight="bold")
axE.legend([Patch(color=VERM), Patch(facecolor="none", edgecolor=VERM, lw=1.2)],
           ["C  nb_pairwise, top-200 markers", "C0  nb_pairwise, no marker pre-sel."],
           loc="upper left", frameon=False, fontsize=ANN, handlelength=1.1, labelspacing=0.3)
# the mechanism/decision prose moved to the legend (sidecar); the two shares are
# already printed on the bars
style(axE)
axE.grid(False, axis="x")

# ---- panel d: permutation-calibrated q -----------------------------------
H = 0.30   # was 0.24: the three per-arm value labels touched at the 6 pt floor
for yi, arm in enumerate(ORDER):
    s = stats.loc[arm]
    d = ARMS[arm]
    y = len(ORDER) - 1 - yi
    vals = [int(s.true_q_hits), int(s.true_perm_q_hits), int(s.true_perm_q_hits_with_effect_floor)]
    for k, (v, hatch) in enumerate(zip(vals, [None, "////", None])):
        if d["markers"]:
            kw = dict(color=d["color"]) if k != 2 else dict(color=d["color"], alpha=0.45)
            if hatch:
                kw.update(hatch=hatch, edgecolor="white", lw=0)
        else:
            kw = dict(facecolor="none", edgecolor=d["color"], lw=1.1)
            if hatch:
                kw.update(hatch=hatch)
            if k == 2:
                kw.update(ls=(0, (2, 1.2)))
        axD.barh(y + (1 - k) * H, v, height=H * 0.92, zorder=3, **kw)
        axD.text(v * 1.06, y + (1 - k) * H, f"{v:,}", va="center", fontsize=ANN, color=INK)
DSHORT = {"A_fisher_reads": "A  fisher reads\n+ top-200 markers",
          "B_fisher_cells": "B  fisher cells\n+ top-200 markers",
          "C_nb_pairwise": "C  nb_pairwise\n+ top-200 markers",
          "A0_fisher_reads_nomarker": "A0  fisher reads\nno marker pre-sel.",
          "B0_fisher_cells_nomarker": "B0  fisher cells\nno marker (SHIPS)",
          "C0_nb_pairwise_nomarker": "C0  nb_pairwise\nno marker pre-sel."}
axD.set_yticks(range(len(ORDER)))
axD.set_yticklabels([DSHORT[a] for a in reversed(ORDER)], fontsize=ANN, linespacing=1.05)
axD.set_xscale("log")
axD.set_xlim(55, 4.2e5)
plain_log(axD, "x", compact=True)   # 100 / 1k / 10k / 100k, not 10^n: a mathtext
                                    # exponent lands at 4.55 pt (judge finding 3)
axD.set_xlabel(sentence_case("significant tests in the TRUE run (log scale)"))
axD.set_title(sc_title("d   permutation-calibrated q vs nominal BH q (TRUE runs)"),
              loc="left", fontweight="bold")
axD.legend([Patch(color=MUTED), Patch(color=MUTED, hatch="////", edgecolor="white", lw=0),
            Patch(color=MUTED, alpha=0.45)],
           ["top bar of each arm: nominal BH q<0.05",
            "middle bar: permutation-calibrated q_perm<0.05",
            "bottom bar: q_perm<0.05 AND effect floor (|Δprop|≥0.1 / |log2FC|≥1)"],
           loc="upper right", bbox_to_anchor=(0.985, 0.70), frameon=False, fontsize=ANN,
           handlelength=1.2, labelspacing=0.3)
# the conservative-null / resolution / comparability prose moved to the legend
# (sidecar); a short label keeps the B0 anomaly visible on the image
axD.text(0.985, 0.97, "B0: q_perm ADDS hits", transform=axD.transAxes, fontsize=ANN,
         color=INK, ha="right", va="top", fontweight="bold")
style(axD)
axD.grid(False, axis="y")

# ---- legend (moved OFF the image; the .caption.md sidecar is its single source)
rt_line = (f"Permutation-harness cost (fig4_calibration_stats.tsv): TRUE runs "
           f"{rt.true_run_elapsed_s.min():.0f}-{rt.true_run_elapsed_s.max():.0f} s per arm, "
           f"mean permutation run {rt.perm_run_elapsed_mean_s.min():.0f}-"
           f"{rt.perm_run_elapsed_mean_s.max():.0f} s (20 permutations per arm, shared seeds).")
legend = (
    "Figure S9 | Calibration deep-dive behind Fig 4 (fdr_calibration_v2 run; testis mouse 1, 1,230 cells, "
    "76,711 PAS, 3 stage pairs, 20 shared label permutations; manuscript/14, verified SOUND). "
    "(a) Null p<0.05 per stage pair for all six switch-test configurations: no pair drives any verdict; "
    "the shipping configuration B0 (Fisher, count-mode cells, no marker pre-selection) is 2.94-3.13% with "
    "0/20 null runs carrying any q<0.05 hit. (b) B0 stratified by a label-independent expression measure: "
    "every stratum is below the 5% nominal rate -- the test is never anti-conservative at any expression "
    "level. 14's verified range (2.6-4.0%) is the shaded band; the bars recompute it from the archived "
    "per-test null p-values and input matrix under the declared binning (2.5-4.0%; the verifier's exact "
    "bin edges were not persisted -- stated caveat, method in the audit TSV). "
    "(c) count-mode reads pseudoreplicates UMIs within cells: 20.3% vs 13.0% null p<0.05 with marker "
    "pre-selection on, 9.6% vs 3.0% with it off, and ~1,500 false q<0.05 hits per null run in reads mode "
    "even without pre-selection. (d) Permutation-calibrated q (empirical p = rank in the pooled same-pair "
    "null, BH per pair) vs nominal BH q on the TRUE runs; for B0 the permutation null is conservative and "
    "q_perm adds hits. (e) nb_pairwise's plug-in dispersion floor (1e-4): 8.5% of C0's null tests carry "
    "67% of its false q<0.05 hits (13% in C; min null p 4e-182) -- whenever marker pre-selection or "
    "nb_pairwise is used, permutation-calibrated q is required (14 decision). "
    "In panel d, B0's permutation null is conservative, so q_perm ADDS hits (66,630 vs 60,332 nominal); "
    "the q_perm resolution is 1/(N_null+1): "
    f"{float(stats.loc['A_fisher_reads','min_p_emp']):.5f} for the marker arms, "
    f"{float(stats.loc['B0_fisher_cells_nomarker','min_p_emp']):.1e} without pre-selection. TRUE hit "
    "counts are NOT comparable across marker modes: pre-selection changes the tested universe AND the "
    "Fisher gene denominator (14 mechanism 2) -- it is a different test. " + rt_line + " "
    "Caveats travel from 14: single mouse, single tissue, very large true stage effects; the "
    "label-permutation null tests exchangeability only."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = TSVDIR / f"{NAME}.png"
fig.savefig(p, dpi=600)
print("wrote", p)

# ---------------------------------------------------------------------------
# audit TSVs
# ---------------------------------------------------------------------------
outs = [("per_pair", pd.DataFrame(pp_rows)), ("strata", strata),
        ("countmode", pd.DataFrame(cm_rows)), ("permq", pd.DataFrame(pq_rows)),
        ("dispfloor", pd.DataFrame(df_rows)), ("runtime", rt)]
for suffix, df in outs:
    p = TSVDIR / f"{NAME}_{suffix}.tsv"
    df.to_csv(p, sep="\t", index=False, float_format="%.8g")
    print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
smin, smax = pct1(strata.frac_null_p_lt_05.min()), pct1(strata.frac_null_p_lt_05.max())
cap_md = f"""# Fig S9 — `{NAME}` sidecar (generated by `scripts/manuscript_figures/{NAME}.py`; single source of the legend)

## Legend

{legend}

**Panel a** — per-stage-pair null p<0.05 from `fig4_calibration_per_pair.tsv`: B0 spans
{b0p.null_frac_p_lt_05.min()*100:.2f}-{b0p.null_frac_p_lt_05.max()*100:.2f}% with `null_frac_q_lt_fdr` 0 in
every pair; B spans {bp.null_frac_p_lt_05.min()*100:.1f}-{bp.null_frac_p_lt_05.max()*100:.1f}% — matching 14
§Per-pair breakdown ("no pair drives the result").
**Panel b** — recomputed stratification of B0's {strata.n_null_tests.sum():,} null tests. Method (declared in
full in `{NAME}_strata.tsv`): {STRAT_METHOD}. Result {smin}-{smax}% across the six strata (pooled rate equals
the verified 3.03% exactly, asserted). 14 records 2.6-4.0% for "6 strata, <50 to ≥600 expressing cells"
without persisting the per-stratum values or exact edges; the recomputation reproduces the upper end exactly
and the lower end to 0.1 pp, and reproduces the verified conclusion — the calibrated arm is never
anti-conservative at any expression level (every stratum < 5% nominal, far under the 7% gate).
**Panel c** — the reads-vs-cells mechanism (14 §Mechanism 3, pseudoreplication): null p<0.05 20.3% (A) vs
13.0% (B) with top-200 marker pre-selection, 9.6% (A0) vs 3.0% (B0) without; A0 still produces
~{stats.loc['A0_fisher_reads_nomarker','mean_hits_per_null_run']:,.0f} false q<0.05 hits per null run while
B0 produces 0 in 20 runs.
**Panel d** — permutation-calibrated q (14 §Permutation-calibrated q): TRUE-run hits nominal → q_perm →
q_perm+effect-floor: A 307 → 186 → 85; B 668 → 640 → 572; C 1,754 → 1,733 → 1,588; A0 49,923 → 33,325 →
16,983; B0 60,332 → 66,630 → 46,230 (q_perm *adds* hits — the permutation null is conservative for B0);
C0 43,127 → 42,974 → 33,067. Resolution 1/(N_null+1): {float(stats.loc['A_fisher_reads','min_p_emp']):.5f}
(marker arms) / {float(stats.loc['B0_fisher_cells_nomarker','min_p_emp']):.1e} (no-marker arms). TRUE hit
counts are not comparable across marker modes (different tested universe AND a different Fisher gene
denominator — 14 §Mechanism 2 / verifier addition).
**Panel e** — NB dispersion-floor diagnostic (14 §Mechanism 4): tests whose plug-in dispersion hits the 1e-4
floor are {stats.loc['C0_nb_pairwise_nomarker','null_share_dispersion_floor']*100:.1f}% of C0's null tests but
carry {stats.loc['C0_nb_pairwise_nomarker','null_hits_share_dispersion_floor']*100:.0f}% of its null q<0.05
hits ({stats.loc['C_nb_pairwise','null_hits_share_dispersion_floor']*100:.0f}% in C); min null p
{float(stats.loc['C0_nb_pairwise_nomarker','min_null_p']):.0e}.

What ships (14 §decision): Fisher + `--count-mode cells` + `--marker-top-n 0`, nominal BH q (conservative
FDR control in every pair and stratum tested), then the ≥2-sample replication filter and the effect floor.
Whenever marker pre-selection or `nb_pairwise` is used, permutation-calibrated q is required.

Caveats: single mouse, single tissue, very large true stage effects (TRUE counts are detectability upper
bounds, not precision estimates); the label-permutation null tests exchangeability under the global null
only; KS-vs-uniform rejects for every arm (discrete Fisher mass at p=1) — reported, not gated; "calibrated"
is the pre-registered operational rule of 14, not a standard. Panel b is the one recomputed panel — nothing
else on this figure required computation beyond reading the verified TSVs, and no permutation was re-run.

## Provenance

Sources: `manuscript/14_switch_calibration_v2.md` (SOUND); `results/figures/manuscript/fig4_calibration_stats.tsv`,
`fig4_calibration_per_pair.tsv`; `results/fdr_calibration_v2/{{report.json,fig4_calibration_null_pvalues_all_arms.tsv.gz,input/stage_labelled.h5ad}}`.
Every plotted value: `results/figures/manuscript/{NAME}_*.tsv` (source column per row).
Script: `scripts/manuscript_figures/{NAME}.py` (env `FIGS9_REUSE_STRATA=1` reuses an existing strata TSV;
EXPECT asserts run either way). Rendered at PNG 600 dpi / PDF fonttype 42.
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure
annotation, key entry and tick label now sits at or above the 6 pt floor. Axis labels, panel titles and prose
tick labels are sentence-cased through `_pubstyle.sentence_case()`, canonical identifiers preserved and the
lower-case panel letters kept. The figure carried no on-image subtitle to move. Three collisions that the
larger type exposed were fixed by geometry, never by dropping content: panel **d** gained row height and its
per-arm bar spacing went 0.24 -> 0.30 so its three count labels clear each other; panel **a**'s no-marker arm
labels are wrapped on three lines instead of two (same words); and the left margin was widened for panel d's
arm labels. No panel, number or audit TSV changed; all six TSVs regenerate byte-identical.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

# ---------------------------------------------------------------------------
# edge check
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS9: B0 3.0% null p<0.05 (per-pair 2.94-3.13), strata {smin}-{smax}% "
      f"(doc 2.6-4.0), perm-q B0 66,630 vs 60,332 nominal, C0 floor 8.5%->67% of hits")
