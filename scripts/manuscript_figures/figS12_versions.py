#!/usr/bin/env python3
"""
figS12_versions.py -- the four-way version-comparison figure:
PeakATail  shipped -> v1 -> v2 -> prime, on four libraries, one scorer.

WHAT THE FOUR VERSIONS ARE (manuscript/24_prime_preregistration.md section 1,
results/prime_bench/README.md section 1)
  shipped  pre-Stage-1 develop, coverage-only `lambda_gradient` caller.  No clip
           evidence, no tiers, no internal-priming veto.  Its BED score column is
           0 for every call, so it has no ranked default output at all.
  v1       4efeb125, Stage 1b/1c clip-seeded, pre-IP-strand-fix.
           manuscript/15_final_gate.md (verified SOUND).
  v2       9dfdefb3 = + IP-strand fix (#96) + performance (#97).
           **THE MANUSCRIPT'S CURRENT RECORD.**  manuscript/19_final_gate_v2.md
           (verifier verdict FIXED); pbmc4k added by manuscript/26.
  prime    6954082d, branch `peakAtail-prime`, **a development branch, NOT merged
           into develop and NOT the manuscript's default.**  Every prime number is
           EXPLORATORY with respect to the manuscript's gates (24 section 3.4).

SOURCE OF TRUTH -- nothing numeric is typed by hand.
  results/prime_bench/fourway.tsv          all four versions, one scorer
                                           (score_tool.py, stage2_final_launch.sh
                                           reference args), 2,950 tidy rows, each
                                           carrying its own source file.
  results/figures/manuscript/fig2_accuracy.tsv   the competitor panel, i.e. the
                                           audit TSV of the manuscript's Fig 2
                                           (same scorer, same denominators).
Consistency asserts below re-check the headline v2 numbers against 19/24, re-derive
the pre-registered criterion of 24 section 3.1 from the table, and refuse to run if
prime and v2 are not exactly equal on the flagged arm.

PANELS
  a  progression on the FULL call set -- the output every version emits.  x = version,
     solid = atlas-agreement precision @100 bp, dashed = detected-gene recall @100 bp,
     colour = dataset.  shipped -> v1 is the Stage-1 clip evidence; v1 -> v2 is the
     minus-strand IP fix (#96); v2 -> prime is exactly zero.
  b  progression at the PRE-REGISTERED PRECISION DEFAULT (tier-1, IP-pass, >=2 clip
     molecules).  shipped cannot produce this output.  The hollow markers are v2 run
     with NO behaviour flag -- what a user gets by typing nothing -- and the arrow to
     prime is the only thing prime's default actually changes.
  c  precision/recall plane, PBMC 10k v3, with the competitor panel as context, the
     F1_det isolines and the pre-registered gate P@100 >= 0.50 (13 section 1).
  d  the same for GSE104556 testis (two mice, mean +- range).
  e  MATCHED CALL COUNT -- the honest way to compare call sets of different sizes.
     Every version's full set ranked under its own key and truncated to a common N.
     Two keys, because the four versions do not share one: UMI depth (all four) and
     clip molecules (v1/v2/prime only -- shipped has no such column).
  f  compute: wall time (upper strip, lollipops, log) over peak RSS (bars, log).
     Two stacked axes, not a twin axis -- on a twin axis the wall-time markers land
     inside the bars and cover the value labels.
  g  the reliability axis: does the curated atlas agree with the long reads?  Atlas
     precision and Kinnex x3p concordance / internal-priming decoy rate, per version,
     on the precision default.

OUTPUTS
  manuscript/figures/figS12_versions.{png,pdf}      (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS12_versions.caption.md     (sidecar: single source of the
      legend, + provenance and the index para; the image carries no caption prose)
  results/figures/manuscript/figS12_versions.tsv               panels a-e, g
  results/figures/manuscript/figS12_versions_compute.tsv       panel f
  results/figures/manuscript/figS12_versions_criterion.tsv     24 section 3.1
  results/figures/manuscript/figS12_versions_reference_lines.tsv

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/figS12_versions.py
"""
import os
import json
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
NAME = "figS12_versions"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"

# Okabe-Ito colourblind-safe palette (manuscript/figures/README.md convention).
VCOL = {                       # colour = code version (panels c-g)
    "shipped": "#999999",
    "v1":      "#E69F00",
    "v2":      "#0072B2",
    "prime":   "#009E73",
    "v2_noIP": "#CC79A7",      # v2 with no behaviour flag = what a user gets
}
DCOL = {                       # colour = dataset (panels a-b)
    "pbmc10k": "#0072B2",
    "mouse1":  "#D55E00",
    "mouse2":  "#E69F00",
    "pbmc4k":  "#56B4E9",
}
TOOLCOL = {                    # competitor panel, same as Fig 2
    "polyApipe":  "#D55E00",
    "SCAPTURE":   "#009E73",
    "Sierra":     "#CC79A7",
    "scAPAtrap":  "#E69F00",
    "scUTRquant": "#56B4E9",
}
DLABEL = {"pbmc10k": "PBMC 10k v3", "mouse1": "testis mouse 1",
          "mouse2": "testis mouse 2", "pbmc4k": "pbmc4k (donor 2)"}
VERSIONS = ["shipped", "v1", "v2", "prime"]
VTICK = {"shipped": "shipped\n(pre-Stage-1)", "v1": "v1\n4efeb125",
         "v2": "v2\n9dfdefb3", "prime": "prime\n6954082d"}
GATE_P = 0.50                       # 13 section 1, pre-registered precision-first gate
F1_ISO = [0.1, 0.2, 0.3, 0.4]

# ---------------------------------------------------------------------------
# sources
# ---------------------------------------------------------------------------
FOURWAY = WD / "results/prime_bench/fourway.tsv"
FIG2 = OUTDIR / "fig2_accuracy.tsv"
# the adversarial verifier's matched-budget control for the flagless pair: the ONLY
# pair on this figure whose two arms differ, compared at matched atlas precision,
# matched atlas recall and matched call count (results/prime_bench/README.md section 10
# caveat 4 records that the benchmark itself measured matched CALL COUNT only, and
# deferred matched precision to a dev-slice run).  Panel g needs it because its bars
# are NOT at equal call count.
VETO = WD / "results/prime_bench/adversarial/veto_matched_pbmc10k.tsv"
assert FOURWAY.exists() and FIG2.exists() and VETO.exists()
FW = pd.read_csv(FOURWAY, sep="\t")
F2 = pd.read_csv(FIG2, sep="\t")
VT = pd.read_csv(VETO, sep="\t")

plotted = []          # every value that ends up on the page


def val(dataset, version, arm, metric, kind="recomputed"):
    """One value out of fourway.tsv, with its own provenance."""
    q = FW[(FW.dataset == dataset) & (FW.version == version) & (FW.arm == arm)
           & (FW.metric == metric) & (FW.source_kind == kind)]
    assert len(q) == 1, (dataset, version, arm, metric, kind, len(q))
    return float(q.value.iloc[0]), str(q.source_file.iloc[0])


def have(dataset, version, arm, metric, kind="recomputed"):  # noqa: D401
    q = FW[(FW.dataset == dataset) & (FW.version == version) & (FW.arm == arm)
           & (FW.metric == metric) & (FW.source_kind == kind)]
    return len(q) == 1


def rec(panel, dataset, version, arm, metric, value, source, ranking="", note=""):
    plotted.append(dict(panel=panel, dataset=dataset, version=version, arm=arm,
                        ranking=ranking, metric=metric, value=value, note=note,
                        source_file=source))
    return value


def get(panel, dataset, version, arm, metric, ranking="", note="", kind="recomputed"):
    v, s = val(dataset, version, arm, metric, kind=kind)
    return rec(panel, dataset, version, arm, metric, v, s, ranking, note)


def veto(comparison, arm, metric):
    """One value out of the matched-budget control TSV, with its own provenance."""
    q = VT[(VT.comparison == comparison) & (VT.arm == arm)]
    assert len(q) == 1, (comparison, arm, len(q))
    return float(q[metric].iloc[0]), str(VETO)


def getv(panel, comparison, arm, metric, note=""):
    v, s = veto(comparison, arm, metric)
    return rec(panel, "pbmc10k", arm, comparison, metric, v, s, note=note)


# ---------------------------------------------------------------------------
# consistency with the verified record before anything is drawn
# ---------------------------------------------------------------------------
# 1. the manuscript's current headline (19 section 1 / 24 section 2), v2 default
EXPECT_V2 = {"pbmc10k": (46524, 0.706195, 0.175443),
             "mouse1":  (26255, 0.745001, 0.204790),
             "mouse2":  (26526, 0.757219, 0.208002)}
for ds, (n, p, r) in EXPECT_V2.items():
    assert int(val(ds, "v2", "precision_default", "n_scored")[0]) == n
    assert round(val(ds, "v2", "precision_default", "P@100")[0], 6) == p
    assert round(val(ds, "v2", "precision_default", "R_det@100")[0], 6) == r
# 2. the 44 values fourway.tsv carries as `quoted` from 15/19/24/26 must equal the
#    re-scored twins to document rounding -- the benchmark's own cross-check, re-run here
#    (the quoted values carry the document's own rounding, so each is compared at
#     half a unit in its own last printed decimal place -- nothing is loosened further)
FW_TXT = pd.read_csv(FOURWAY, sep="\t", dtype=str)
_nq = 0
for i, r in FW_TXT[FW_TXT.source_kind == "quoted"].iterrows():
    m = FW[(FW.dataset == r.dataset) & (FW.version == r.version) & (FW.arm == r.arm)
           & (FW.metric == r.metric) & (FW.source_kind == "recomputed")]
    if len(m) != 1:
        continue
    txt = r.value.rstrip("0").rstrip(".") if "." in r.value else r.value
    dec = len(txt.split(".")[1]) if "." in txt else 0
    tol = max(1e-9, 0.5 * 10 ** (-dec) + 1e-9)
    assert abs(float(m.value.iloc[0]) - float(r.value)) <= tol, (r.to_dict(), float(m.value.iloc[0]), tol)
    _nq += 1
assert _nq >= 40, _nq
# 3. prime must be EXACTLY v2 on every flagged arm -- if this ever stops being true
#    the figure's central claim has changed and the script must not silently redraw.
ZERO_ARMS = ["full", "tier1", "precision_default"]
for ds in ("pbmc10k", "mouse1", "mouse2", "pbmc4k"):
    for arm in ZERO_ARMS:
        for m in ("n_scored", "P@100", "R_det@100", "F1_det@100"):
            if have(ds, "prime", arm, m) and have(ds, "v2", arm, m):
                assert val(ds, "prime", arm, m)[0] == val(ds, "v2", arm, m)[0], (ds, arm, m)
# 4. v1_noIP and v2_noIP are the same call set (19 section 2)
for m in ("n_scored", "P@100", "R_det@100"):
    assert val("pbmc10k", "v1_noIP", "full", m)[0] == val("pbmc10k", "v2_noIP", "full", m)[0]

# ---- the pre-registered PRIMARY criterion, 24 section 3.1, re-derived ------
CRIT = []
for ds in ("pbmc10k", "mouse1", "mouse2"):
    dP = val(ds, "prime", "precision_default", "P@100")[0] - val(ds, "v2", "precision_default", "P@100")[0]
    dR = val(ds, "prime", "precision_default", "R_det@100")[0] - val(ds, "v2", "precision_default", "R_det@100")[0]
    dF = val(ds, "prime", "precision_default", "F1_det@100")[0] - val(ds, "v2", "precision_default", "F1_det@100")[0]
    CRIT.append(dict(dataset=ds, dP=dP, dR=dR, dF1=dF,
                     i_dP_ge_minus0p005=bool(dP >= -0.005),
                     ii_dR_ge_plus0p010=bool(dR >= 0.010),
                     iii_dF1_gt_0=bool(dF > 0.0),
                     verdict="PASS" if (dP >= -0.005 and dR >= 0.010 and dF > 0.0) else "FAIL"))
CRIT_OVERALL = "PASS" if all(c["verdict"] == "PASS" for c in CRIT) else "FAIL"
assert CRIT_OVERALL == "FAIL", "24 section 3.1 verdict changed -- rewrite the figure text"
assert all(c["dP"] == 0.0 and c["dR"] == 0.0 and c["dF1"] == 0.0 for c in CRIT)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# canvas: the footer caption band and the long panel subtitles moved to the
# legend sidecar; height shrinks by the freed footer/subtitle space (hspace was
# sized for 4-5-line subtitles, now at most one line), axes keep their size
fig = plt.figure(figsize=(9.6, 11.45))
gs = fig.add_gridspec(4, 2, height_ratios=[1.00, 1.14, 0.92, 0.92],
                      left=0.072, right=0.984, top=0.958, bottom=0.037,
                      hspace=0.55, wspace=0.22)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, :])
# panel f is two stacked sub-axes: the wall-time strip must not share a y-axis with the
# RSS bars -- on a twin log axis the wall-time markers land INSIDE the bars (measured:
# 8 of 14 markers) and cover the in-bar value labels.
gsF = gs[3, 0].subgridspec(2, 1, height_ratios=[0.40, 0.60], hspace=0.14)
axFw = fig.add_subplot(gsF[0])
axF = fig.add_subplot(gsF[1])
axG = fig.add_subplot(gs[3, 1])
ref_lines = []
TITLES = []
# Panel subtitles are no longer drawn (DESIGN_DIRECTIVES.md item 2).  Every one is
# collected here verbatim and written into the caption sidecar's '## Legend' under
# 'Panel notes' -- the no-loss rule: nothing is dropped, only relocated.
SUBTITLES = []


def style(ax, xgrid=True, ygrid=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    if ygrid:
        ax.grid(True, axis="y", color=GRID, lw=0.5, alpha=0.8)
    if xgrid:
        ax.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.8)
    ax.set_axisbelow(True)


def titled(ax, title, subtitle):
    """Bold, sentence-cased panel title -- and nothing else on the image.

    Descriptive subtitles left the canvas in the 2026-09-03 design pass
    (DESIGN_DIRECTIVES.md item 2); each one is recorded in SUBTITLES and written
    verbatim into the caption sidecar's '## Legend' -> 'Panel notes'."""
    ax.set_title(sc_title(title), loc="left", fontweight="bold", pad=7.0)
    TITLES.append((ax, title))
    if subtitle:
        SUBTITLES.append((title.strip().split()[0], title, subtitle))


# ===========================================================================
# panels a / b -- the progression
# ===========================================================================
DS_A = ["pbmc10k", "mouse1", "mouse2", "pbmc4k"]
LAB_DY = {"pbmc10k": -7, "mouse1": -7, "mouse2": 7, "pbmc4k": 7}


def progression(ax, panel, arm, datasets, versions):
    xs = {v: i for i, v in enumerate(versions)}
    for ds in datasets:
        c = DCOL[ds]
        for metric, ls, mk, mfc in (("P@100", "-", "o", c),
                                    ("R_det@100", (0, (3, 1.6)), "s", "white")):
            X, Y = [], []
            for v in versions:
                if not have(ds, v, arm, metric):
                    continue
                X.append(xs[v])
                Y.append(get(panel, ds, v, arm, metric))
            if not X:
                continue
            ax.plot(X, Y, ls=ls, color=c, lw=1.3, zorder=4)
            ax.plot(X, Y, ls="none", marker=mk, ms=5.0, mfc=mfc, mec=c, mew=1.2, zorder=5)
    ax.set_xticks(range(len(versions)))
    ax.set_xticklabels([VTICK[v] for v in versions], fontsize=TYPE["tick"])
    ax.set_xlim(-0.48, len(versions) - 1 + 0.62)
    ax.set_ylabel(sentence_case("value @100 bp"))
    style(ax, xgrid=False)


# ---- panel a: the full call set -------------------------------------------
progression(axA, "a", "full", DS_A, VERSIONS)
axA.set_ylim(0, 0.80)
axA.axvspan(2 - 0.13, 3 + 0.13, color="#EFF2F3", zorder=0)
axA.annotate("prime = v2 — Δ = 0.000000\non every metric, every dataset",
             xy=(2.42, 0.565), ha="center", va="bottom", fontsize=6.0, color=INK, zorder=6)
for ds in DS_A:
    for v, dx, ha in (("shipped", -6, "right"), ("prime", 6, "left")):
        if have(ds, v, "full", "P@100"):
            y = val(ds, v, "full", "P@100")[0]
            axA.annotate(f"{y:.3f}", (VERSIONS.index(v), y), xytext=(dx, LAB_DY[ds]),
                         textcoords="offset points", fontsize=ANN, color=DCOL[ds],
                         ha=ha, va="center")
# the flag-matching explanation moved to the legend (sidecar); the metric
# encoding is already carried by the in-panel legend
titled(axA, "a   Progression on the full call set", "")
axA.legend([Line2D([], [], color=DCOL[d], lw=1.3, marker="o", ms=4.5, mfc=DCOL[d], mec=DCOL[d])
            for d in DS_A]
           + [Line2D([], [], color=MUTED, lw=1.3, marker="o", ms=4.5, mfc=MUTED, mec=MUTED),
              Line2D([], [], color=MUTED, lw=1.3, ls=(0, (3, 1.6)), marker="s", ms=4.5,
                     mfc="white", mec=MUTED)],
           [DLABEL[d] for d in DS_A]
           + ["atlas-agreement precision", "detected-gene recall"],
           loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, ncol=2,
           handlelength=1.9, labelspacing=0.28, columnspacing=0.9, fontsize=ANN)

# ---- panel b: the pre-registered precision default -------------------------
progression(axB, "b", "precision_default", DS_A, VERSIONS)
axB.set_ylim(0, 1.02)
axB.axvspan(-0.48, 0.5, color="#F5F2EF", zorder=0)
axB.axvspan(2 - 0.13, 3 + 0.13, color="#EFF2F3", zorder=0)
# full explanation (and the SCAPTURE parallel) moved to the legend (sidecar)
# same words on four narrower lines: at the 6 pt floor the third line reached the
# boxed v2-no-flag note on its right (design pass 2026-09-03)
axB.text(-0.42, 0.30, "shipped: no ranked\ndefault output\n(BED score\ncolumn = 0)",
         fontsize=ANN, color=MUTED, ha="left", va="center")
# what a user gets by typing nothing: v2's flagless >=2-molecule output (PBMC only)
p_nf = get("b", "pbmc10k", "v2_noIP", "ge2mol_noIP_POSTHOC", "P@100", note="v2 run with no behaviour flag")
r_nf = get("b", "pbmc10k", "v2_noIP", "ge2mol_noIP_POSTHOC", "R_det@100", note="v2 run with no behaviour flag")
p_pr = val("pbmc10k", "prime", "precision_default", "P@100")[0]
r_pr = val("pbmc10k", "prime", "precision_default", "R_det@100")[0]
axB.plot([2], [p_nf], marker="o", ms=6.5, mfc="white", mec=VCOL["v2_noIP"], mew=1.5, zorder=7)
axB.plot([2], [r_nf], marker="s", ms=6.0, mfc="white", mec=VCOL["v2_noIP"], mew=1.5, zorder=7)
for y0, y1 in ((p_nf, p_pr), (r_nf, r_pr)):
    axB.annotate("", xy=(2.94, y1), xytext=(2.06, y0), zorder=7,
                 arrowprops=dict(arrowstyle="-|>", color=VCOL["v2_noIP"], lw=1.1,
                                 shrinkA=4, shrinkB=4))
axB.text(0.62, 0.355, f"v2 with NO behaviour flag → prime\n"
         f"(PBMC 10k v3, ≥2 molecules):\n"
         f"P  {p_nf:.4f} → {p_pr:.4f}  ({(p_pr/p_nf-1)*100:+.1f} %)\n"
         f"R  {r_nf:.4f} → {r_pr:.4f}  ({(r_pr/r_nf-1)*100:+.1f} %)",
         fontsize=ANN, color=VCOL["v2_noIP"], ha="left", va="center", zorder=8,
         bbox=dict(facecolor="white", edgecolor=VCOL["v2_noIP"], lw=0.5, pad=2.4, alpha=0.95))
axB.text(2.5, 0.995, "prime = v2\nΔ = 0.000000", ha="center", va="top", fontsize=6.0,
         color=INK, zorder=6)
_dP = {ds: val(ds, "v2", "precision_default", "P@100")[0] - val(ds, "v1", "precision_default", "P@100")[0]
       for ds in ("pbmc10k", "mouse1", "mouse2")}
axB.axhline(GATE_P, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
axB.text(3.58, GATE_P + 0.012, "pre-registered gate\nP@100 ≥ 0.50", fontsize=ANN,
         color=INK, ha="right", va="bottom")
ref_lines.append(dict(panel="b", kind="gate", label="pre-registered gate (13 section 1)", value=GATE_P))
# the hollow-marker explanation and the v1 -> v2 deltas moved to the legend
# (sidecar); the default's definition stays as the panel's one-line subtitle
titled(axB, "b   Progression at the pre-registered default",
       "tier-1 ∩ IP-pass ∩ ≥2 clip molecules")
axB.legend([Line2D([], [], color=VCOL["v2_noIP"], lw=0, marker="o", ms=5, mfc="white",
                   mec=VCOL["v2_noIP"], mew=1.5)],
           ["v2 as a user runs it (no behaviour flag)"], loc="lower right",
           bbox_to_anchor=(1.0, 0.0), frameon=False, handlelength=1.4, fontsize=ANN)


# ===========================================================================
# panels c / d -- precision/recall plane against the field
# ===========================================================================
def f1_iso(ax, xmax):
    R = np.linspace(0.005, xmax, 400)
    for f in F1_ISO:
        with np.errstate(divide="ignore", invalid="ignore"):
            P = f * R / (2 * R - f)
        ok = (2 * R - f > 0) & (P <= 1.0) & (P >= 0)
        ax.plot(R[ok], P[ok], color="#C9D1D5", lw=0.7, zorder=1)
        # offset off the contour's own end point: placed on it, the curve ran through
        # the label (design pass 2026-09-03)
        ax.annotate(f"F1 {f:.1f}", (R[ok][-1], P[ok][-1]), xytext=(-2.5, -3.0),
                    textcoords="offset points", color="#9AA7AD", fontsize=ANN,
                    ha="right", va="top", zorder=1)
        ref_lines.append(dict(panel=ax.get_label(), kind="F1_det isoline", label=f"F1 = {f}", value=f))


def plane(ax, panel, f2_dataset, ds_keys, xmax):
    ax.set_label(panel)
    f1_iso(ax, xmax)
    ax.axhline(GATE_P, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.text(xmax * 0.995, GATE_P + 0.009, "pre-registered gate P@100 ≥ 0.50",
            ha="right", va="bottom", fontsize=ANN, color=INK)
    ref_lines.append(dict(panel=panel, kind="gate", label="pre-registered gate (13 section 1)", value=GATE_P))
    sub = F2[(F2.dataset == f2_dataset) & (F2.tool != "PeakATail")]
    hands = []
    for tool in ["scUTRquant", "SCAPTURE", "polyApipe", "Sierra", "scAPAtrap"]:
        q = sub[sub.tool == tool]
        if q.empty:
            continue
        c = TOOLCOL[tool]
        hollow = tool == "scUTRquant"
        x = float(q.recall_detected_genes_100bp.mean())
        y = float(q.atlas_agreement_precision_100bp.mean())
        if len(q) == 2:
            ax.errorbar(x, y, xerr=[[x - q.recall_detected_genes_100bp.min()],
                                    [q.recall_detected_genes_100bp.max() - x]],
                        yerr=[[y - q.atlas_agreement_precision_100bp.min()],
                              [q.atlas_agreement_precision_100bp.max() - y]],
                        color=c, lw=0.8, capsize=1.5, zorder=3, ls="none")
        h, = ax.plot(x, y, marker="^", ms=5.8, mfc="white" if hollow else c, mec=c, mew=1.2,
                     ls="none", zorder=4)
        nm = "scUTRquant* (catalog)" if hollow else tool
        if f2_dataset == "gse104556" and len(q) == 1:
            nm += " (m1 only)"
        hands.append((h, nm))
        for r in q.itertuples():
            for met, v in (("P@100", r.atlas_agreement_precision_100bp),
                           ("R_det@100", r.recall_detected_genes_100bp)):
                plotted.append(dict(panel=panel, dataset=f2_dataset, version=f"competitor:{tool}",
                                    arm="tool call set", ranking="", metric=met, value=v,
                                    note=r.replicate, source_file=r.source))
        ax.plot(x, float(q.null_P_mean.mean()), marker="|", ms=5, mec=c, mew=1.0, ls="none", zorder=3)

    def pt(version, arm):
        P = float(np.mean([get(panel, d, version, arm, "P@100") for d in ds_keys]))
        R = float(np.mean([get(panel, d, version, arm, "R_det@100") for d in ds_keys]))
        N = float(np.mean([get(panel, d, version, arm, "null_P@100_mean3") for d in ds_keys]))
        return R, P, N

    for arm, mk, ms in (("full", "o", 6.2), ("precision_default", "D", 7.4)):
        path = []
        for v in VERSIONS:
            if not all(have(d, v, arm, "P@100") for d in ds_keys):
                continue
            R, P, Nl = pt(v, arm)
            path.append((v, R, P))
            ax.plot(R, Nl, marker="|", ms=5, mec=VCOL[v], mew=1.0, ls="none", zorder=3)
        for (_, x0, y0), (_, x1, y1) in zip(path[:-1], path[1:]):
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=5,
                        arrowprops=dict(arrowstyle="-|>", color="#7E8C93", lw=0.9,
                                        shrinkA=6, shrinkB=8))
        for v, R, P in path:
            if v == "prime":          # prime lands exactly on v2: draw it as a ring
                ax.plot(R, P, marker=mk, ms=ms + 4.6, mfc="none", mec=VCOL["prime"],
                        mew=1.5, zorder=7)
            else:
                ax.plot(R, P, marker=mk, ms=ms, mfc=VCOL[v], mec="white", mew=0.9, zorder=6)
    if f2_dataset == "pbmc_10k_v3":
        # NO arrow here: the flagless point and the prime point are ~21 pt apart while the
        # two markers are ~12 pt across, so any arrow between them degenerates into a
        # stranded arrowhead.  Panel b draws that move, with the numbers.
        ax.plot(r_nf, p_nf, marker="D", ms=7.4, mfc="white", mec=VCOL["v2_noIP"], mew=1.6, zorder=6)
        ax.annotate("v2, no flags", (r_nf, p_nf), xytext=(-7, 0), textcoords="offset points",
                    fontsize=ANN, color=VCOL["v2_noIP"], ha="right", va="center", zorder=8)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 0.95)
    ax.set_xlabel(sentence_case("detected-gene recall R_det @100 bp"))
    ax.set_ylabel(sentence_case("atlas-agreement precision @100 bp"))
    style(ax)
    return hands


# the ring/arrow/scTail explanations moved to the legend (sidecar); a one-line
# key naming the marker shapes stays
hC = plane(axC, "c", "pbmc_10k_v3", ["pbmc10k"], 0.42)
titled(axC, "c   Against the field — PBMC 10k v3",
       "circles = full call set · diamonds = precision default · | = genic-shuffle null (mean)")
axC.annotate("shipped", (val("pbmc10k", "shipped", "full", "R_det@100")[0],
                         val("pbmc10k", "shipped", "full", "P@100")[0]),
             xytext=(0, -9), textcoords="offset points", fontsize=6.0,
             color=VCOL["shipped"], ha="center", va="top")
axC.annotate("v2 = prime\n(full set)", (val("pbmc10k", "v2", "full", "R_det@100")[0],
                                        val("pbmc10k", "v2", "full", "P@100")[0]),
             xytext=(9, -3), textcoords="offset points", fontsize=6.0, color=VCOL["v2"],
             ha="left", va="top")
axC.annotate("v2 = prime\nprecision default", (r_pr, p_pr), xytext=(10, 4),
             textcoords="offset points", fontsize=6.0, color=VCOL["v2"], ha="left",
             va="bottom", fontweight="bold")

# the no-flagless-mouse-arm and SCAPTURE-mouse-2 notes moved to the legend
hD = plane(axD, "d", "gse104556", ["mouse1", "mouse2"], 0.42)
titled(axD, "d   Against the field — GSE104556 testis",
       "two mice, mean ± range")
axD.annotate("shipped", (np.mean([val(d, "shipped", "full", "R_det@100")[0] for d in ("mouse1", "mouse2")]),
                         np.mean([val(d, "shipped", "full", "P@100")[0] for d in ("mouse1", "mouse2")])),
             xytext=(0, -9), textcoords="offset points", fontsize=6.0,
             color=VCOL["shipped"], ha="center", va="top")
axD.annotate("v2 = prime\nprecision default",
             (np.mean([val(d, "v2", "precision_default", "R_det@100")[0] for d in ("mouse1", "mouse2")]),
              np.mean([val(d, "v2", "precision_default", "P@100")[0] for d in ("mouse1", "mouse2")])),
             xytext=(-11, 0), textcoords="offset points", fontsize=6.0, color=VCOL["v2"],
             ha="right", va="center", fontweight="bold")

_leg = hC + [(Line2D([], [], marker="o", ms=5.5, mfc=VCOL[v], mec="white", ls="none"), v)
             for v in VERSIONS]
_leg += [(Line2D([], [], marker="D", ms=6.2, mfc="none", mec=VCOL["prime"], mew=1.5, ls="none"),
          "prime = a ring on v2"),
         (Line2D([], [], marker="D", ms=6.2, mfc="white", mec=VCOL["v2_noIP"], mew=1.6, ls="none"),
          "v2, no behaviour flag")]
fig.legend([h for h, _ in _leg], [n for _, n in _leg], loc="upper center",
           bbox_to_anchor=(0.53, axC.get_position().y0 - 0.033), ncol=6, frameon=False,
           handletextpad=0.4, columnspacing=1.2, labelspacing=0.35, fontsize=6.1)


# ===========================================================================
# panel e -- matched call count
# ===========================================================================
axE.set_label("e")
BLOCKS = [
    ("UMI-depth key, N = shipped's own call count", "umi",
     {"pbmc10k": "matched_umi_N277164", "mouse1": "matched_umi_N45921",
      "mouse2": "matched_umi_N46672"}, VERSIONS),
    ("UMI-depth key, N = v2's precision-default budget", "umiref",
     {"pbmc10k": "matched_umiref_N46544", "mouse1": "matched_umiref_N26263",
      "mouse2": "matched_umiref_N26533"}, VERSIONS),
    ("clip-molecule key, N = v2's precision-default budget", "molref",
     {"pbmc10k": "matched_molref_N46544", "mouse1": "matched_molref_N26263",
      "mouse2": "matched_molref_N26533"}, ["v1", "v2", "prime"]),
]
DS_E = ["pbmc10k", "mouse1", "mouse2"]
bw = 0.20
xpos, xlab, block_span = [], [], []
cursor = 0.0
n_tie_flagged = 0
for btitle, ranking, arms, vers in BLOCKS:
    start = cursor
    for ds in DS_E:
        arm = arms[ds]
        for k, v in enumerate(vers):
            if not have(ds, v, arm, "P@100"):
                continue
            off = (k - (len(vers) - 1) / 2) * bw
            P = get("e", ds, v, arm, "P@100", ranking=ranking)
            R = get("e", ds, v, arm, "R_det@100", ranking=ranking)
            get("e", ds, v, arm, "n_scored", ranking=ranking)
            # tie exposure = sites dropped purely by the tie-break, per output site.
            # >0.5 means more than half an output's worth of tied sites was cut by the
            # tie-break rather than by evidence (25 section 2.3) -- not like-for-like.
            tie = False
            if have(ds, v, arm, "rank_n_tied_at_cut_value", kind="measured"):
                tied = get("e", ds, v, arm, "rank_n_tied_at_cut_value", ranking=ranking, kind="measured")
                kept = get("e", ds, v, arm, "rank_n_kept_at_cut_value", ranking=ranking, kind="measured")
                nw = get("e", ds, v, arm, "rank_n_written", ranking=ranking, kind="measured")
                tie = (tied - kept) / nw > 0.5
            n_tie_flagged += int(tie)
            axE.bar(cursor + off, P, width=bw * 0.92, color=VCOL[v], edgecolor="white",
                    lw=0.6, zorder=3, hatch="////" if tie else None, alpha=0.55 if tie else 1.0)
            axE.plot(cursor + off, R, marker="D", ms=3.6, mfc="white", mec=INK, mew=0.9, zorder=5)
            axE.text(cursor + off, max(P, R) + 0.055, f"{P:.3f}", ha="center", va="bottom",
                     fontsize=ANN, color=INK, rotation=90, zorder=5)
            if tie:
                axE.text(cursor + off, 0.016, "tie-\nbreak", ha="center", va="bottom",
                         fontsize=ANN, color="#7a4b00", zorder=6)
        xpos.append(cursor)
        xlab.append(f"{DLABEL[ds]}\nN = {int(arms[ds].split('_N')[1]):,}")
        cursor += 1.05
    block_span.append((start, cursor - 1.05, btitle))
    cursor += 0.80
assert n_tie_flagged == 3, n_tie_flagged   # exactly the three v1 clip-molecule rows
axE.set_xticks(xpos)
axE.set_xticklabels([sentence_case(t) for t in xlab], fontsize=TYPE["tick"])
axE.set_ylim(0, 1.22)
axE.set_xlim(-0.66, cursor - 1.19)
axE.set_ylabel(sentence_case("atlas-agreement precision\n@100 bp (bars)"))
for s0, s1, btitle in block_span:
    axE.plot([s0 - 0.50, s1 + 0.50], [1.13, 1.13], color=MUTED, lw=0.8)
    axE.text((s0 + s1) / 2, 1.145, sentence_case(btitle), ha="center", va="bottom",
             fontsize=TYPE["tick"], color=INK)
_pu = val("pbmc10k", "v2", "matched_umiref_N46544", "P@100")[0]
_pm = val("pbmc10k", "v2", "matched_molref_N46544", "P@100")[0]
# the two-key rationale, tie-break numbers and the clip-vs-UMI comparison moved
# to the legend (sidecar); the hatch is named in the in-panel legend
titled(axE, "e   Matched call count — each version's own ranking, truncated to a common N", "")
axE.legend([Patch(facecolor=VCOL[v], edgecolor="white") for v in VERSIONS]
           + [Line2D([], [], marker="D", ms=4, mfc="white", mec=INK, ls="none"),
              Patch(facecolor="#bbbbbb", edgecolor="white", hatch="////", alpha=0.55)],
           VERSIONS + ["detected-gene recall R_det @100 bp", "decided by the tie-break"],
           loc="upper left", bbox_to_anchor=(0.003, 0.895), frameon=False, ncol=3,
           handlelength=1.5, labelspacing=0.3, columnspacing=1.3, fontsize=6.0)
style(axE, xgrid=False)


# ===========================================================================
# panel f -- compute
# ===========================================================================
axF.set_label("f")
axFw.set_label("f_wall")
DS_F = ["pbmc10k", "mouse1", "mouse2", "pbmc4k"]
compute_rows = []
bw = 0.20
WALL_FLOOR = 6.0            # bottom of the wall-time strip, minutes (min value is 9.08)
for i, ds in enumerate(DS_F):
    for k, v in enumerate(VERSIONS):
        if not have(ds, v, "RUN", "peak_rss_gb", kind="measured"):
            continue
        rss, src = val(ds, v, "RUN", "peak_rss_gb", kind="measured")
        wall, _ = val(ds, v, "RUN", "wall_seconds", kind="measured")
        cpu, _ = val(ds, v, "RUN", "cpu_percent", kind="measured")
        off = (k - 1.5) * bw
        # --- peak RSS, the quotable figure: bars, value ABOVE the bar (the strip above
        #     carries the wall time now, so nothing can land on top of a bar any more)
        axF.bar(i + off, rss, width=bw * 0.9, color=VCOL[v], edgecolor="white", lw=0.6, zorder=3)
        axF.text(i + off, rss * 1.10, f"{rss:,.1f}", ha="center", va="bottom", fontsize=ANN,
                 color=INK, rotation=90, zorder=5)
        # --- wall time, a run record only: lollipop in its own strip, drawn lighter
        wmin = wall / 60.0
        axFw.plot([i + off, i + off], [WALL_FLOOR, wmin], color=VCOL[v], lw=0.9, zorder=3,
                  solid_capstyle="butt")
        axFw.plot(i + off, wmin, marker="D", ms=3.8, mfc="white", mec=VCOL[v], mew=1.0, zorder=5)
        # offset in POINTS, not in data units: a multiplicative offset on a log axis this
        # short is ~1 pt and the label lands on its own marker
        axFw.annotate(f"{wmin:,.0f}" if wmin >= 20 else f"{wmin:,.1f}", (i + off, wmin),
                      xytext=(0, 4.5), textcoords="offset points", ha="center", va="bottom",
                      fontsize=ANN, color=MUTED, rotation=90, zorder=6)
        compute_rows.append(dict(dataset=ds, version=v, peak_rss_gb=rss, wall_seconds=wall,
                                 wall_minutes=wall / 60.0, cpu_percent=cpu, source_file=src))
axF.set_yscale("log")
axF.set_ylim(1.0, 2500)
axF.set_xticks(range(len(DS_F)))
axF.set_xticklabels([sentence_case(DLABEL[d]) for d in DS_F], fontsize=TYPE["tick"])
axF.set_xlim(-0.5, len(DS_F) - 0.5)
# two lines: rotated, the one-line form is taller than this short sub-axes (the
# script's own y-label guard below catches it) -- same words (design pass 2026-09-03)
axF.set_ylabel(sentence_case("peak RSS\n(GB, log)"), fontsize=TYPE["axis_label"])
axFw.set_yscale("log")
axFw.set_ylim(WALL_FLOOR, 6000)
axFw.set_xlim(-0.5, len(DS_F) - 0.5)
axFw.set_xticks(range(len(DS_F)))
axFw.set_xticklabels([])
axFw.tick_params(labelbottom=False)
# the strip is short: a longer y-label overruns it and collides with the RSS label
# below (asserted at the bottom of this script)
axFw.set_ylabel(sentence_case("wall (min)"), color=MUTED, fontsize=TYPE["axis_label"])
_r1 = val("pbmc10k", "v1", "RUN", "peak_rss_gb", kind="measured")[0]
_r2 = val("pbmc10k", "v2", "RUN", "peak_rss_gb", kind="measured")[0]
_r1n = val("pbmc10k", "v1_noIP", "RUN", "peak_rss_gb", kind="measured")[0]
_r2n = val("pbmc10k", "v2_noIP", "RUN", "peak_rss_gb", kind="measured")[0]
_m1v2 = val("mouse1", "v2", "RUN", "peak_rss_gb", kind="measured")[0]
_m1pr = val("mouse1", "prime", "RUN", "peak_rss_gb", kind="measured")[0]
# the #97 RSS story and the concurrency detail moved to the legend (sidecar);
# the binding wall-time qualifier stays as a one-line label of the upper strip
titled(axFw, "f   What each step cost — compute",
       "wall time (upper strip) is a run record, not a benchmark; both strips log")
axF.text(0.985, 0.985, f"mouse 1 peak RSS moves {_m1v2:.2f} → {_m1pr:.2f} GB ({_m1pr/_m1v2:.2f}×)\n"
         "v2 → prime, inside the pre-registered 1.5× guard rail (24 §3.2)",
         transform=axF.transAxes, fontsize=ANN, color=INK, ha="right", va="top")
style(axF, xgrid=False)
style(axFw, xgrid=False)
axFw.tick_params(axis="x", length=0)
axFw.legend([Line2D([], [], color=VCOL[v], lw=0.9, marker="D", ms=3.6, mfc="white",
                    mec=VCOL[v], mew=1.0) for v in VERSIONS], VERSIONS,
            loc="upper right", bbox_to_anchor=(1.0, 1.04), frameon=False, ncol=4,
            handlelength=1.3, columnspacing=1.0, fontsize=ANN)


# ===========================================================================
# panel g -- atlas vs long reads
# ===========================================================================
axG.set_label("g")
GMET = [("P@100", "atlas agreement\n@100 bp", False),
        ("kinnex_t5_frac25", "Kinnex ≥5 UMI\n@25 bp", False),
        ("kinnex_t20_frac25", "Kinnex ≥20 UMI\n@25 bp", False),
        ("kinnex_decoy_frac25", "Kinnex IP decoy\n@25 bp", True)]
GVERS = [("v1", "precision_default"), ("v2", "precision_default"),
         ("prime", "precision_default"), ("v2_noIP", "ge2mol_noIP_POSTHOC")]
bw = 0.19
for i, (metric, mlab, lower_better) in enumerate(GMET):
    for k, (v, arm) in enumerate(GVERS):
        y = get("g", "pbmc10k", v, arm, metric)
        off = (k - 1.5) * bw
        axG.bar(i + off, y, width=bw * 0.9, color=VCOL[v], edgecolor="white", lw=0.6,
                zorder=3, hatch="\\\\\\" if v == "prime" else None)
        axG.text(i + off, y + 0.014, f"{y:.3f}", ha="center", va="bottom", fontsize=ANN,
                 color=INK, rotation=90, zorder=5)
    if lower_better:
        axG.text(i, 0.52, "lower is better ↓", ha="center", va="bottom", fontsize=ANN, color=INK)
axG.set_xticks(range(len(GMET)))
axG.set_xticklabels([sentence_case(m[1]) for m in GMET], fontsize=TYPE["tick"])
axG.set_ylim(0, 1.35)
axG.set_ylabel(sentence_case("fraction"))
# the four bars are NOT at equal call count -- say so in the legend, because a
# fraction-of-calls metric falls with n and the pink bar's set is 33 % larger
_gn = {v: int(get("g", "pbmc10k", v, arm, "n_scored")) for v, arm in GVERS}
# matched-budget control: truncate v2-no-flags to prime's own N and re-read all four
_mc = {m: getv("g", "matched_calls", "v2_flagless_topN46524", m,
               note="matched-budget control: v2 no-flags truncated by clip-molecule rank "
                    "to prime's own N = 46,524")
       for m in ("P_atlas@100", "kinnex_t5@25", "kinnex_t20@25", "kinnex_decoy@25")}
_mp_r = getv("g", "matched_precision", "v2_flagless_topN40687", "R_det@100",
             note="matched-precision control: v2 no-flags truncated to P@100 = prime's 0.706195")
# the interpretation and matched-budget-control numbers moved to the legend
# (sidecar); the dataset scope stays as a one-line subtitle
titled(axG, "g   Does the atlas agree with the long reads?",
       "PBMC 10k v3, ≥2-molecule output — human-only (the mice have no long-read truth)")
axG.legend([Patch(facecolor=VCOL["v1"], edgecolor="white"),
            Patch(facecolor=VCOL["v2"], edgecolor="white"),
            Patch(facecolor=VCOL["prime"], edgecolor="white", hatch="\\\\\\"),
            Patch(facecolor=VCOL["v2_noIP"], edgecolor="white")],
           ["v1 default", "v2 default = the manuscript", "prime default (= v2, hatched)",
            "v2 typed with no flags"],
           loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, ncol=1,
           handlelength=1.5, labelspacing=0.32, fontsize=ANN)
_sh_t5 = get("g", "pbmc10k", "shipped", "full", "kinnex_t5_frac25",
             note="shipped has no >=2-molecule arm; its FULL call set is shown for scale only")
# the scale numbers and the human-only note moved to the legend / subtitle; a
# short absence label stays
axG.text(0.995, 1.0, "shipped: no ≥2-molecule output",
         transform=axG.transAxes, fontsize=ANN, color=MUTED, ha="right", va="top")
style(axG, xgrid=False)


# ===========================================================================
# legend (moved OFF the image; the .caption.md sidecar is its single source)
# ===========================================================================
_v2p = [val(d, "v2", "precision_default", "P@100")[0] for d in ("pbmc10k", "mouse1", "mouse2")]
_v1p = [val(d, "v1", "precision_default", "P@100")[0] for d in ("pbmc10k", "mouse1", "mouse2")]
_k5_nf = val("pbmc10k", "v2_noIP", "ge2mol_noIP_POSTHOC", "kinnex_t5_frac25")[0]
_k5_pr = val("pbmc10k", "prime", "precision_default", "kinnex_t5_frac25")[0]
_dc_nf = val("pbmc10k", "v2_noIP", "ge2mol_noIP_POSTHOC", "kinnex_decoy_frac25")[0]
_dc_pr = val("pbmc10k", "prime", "precision_default", "kinnex_decoy_frac25")[0]

legend_txt = (
    "Figure S12 | Four versions of the PeakATail caller on four libraries, scored once. shipped (pre-Stage-1, "
    "coverage-only) → v1 (4efeb125, clip-seeded; 15) → v2 (9dfdefb3, + minus-strand IP fix #96 + performance "
    "#97; 19, verifier FIXED, pbmc4k added by 26) → prime (6954082d, branch peakAtail-prime). "
    "**v2 is the manuscript's current record. prime is a development branch: not merged, not the default, and "
    "every prime value here is EXPLORATORY with respect to the manuscript's gates (24 §3.4).** "
    "Flags (PBMC arm): shipped had no clip evidence, no tiers and no internal-priming veto; v1 and v2 ran "
    "`--peak-strategy clip_seeded --ip-filter --genome-fasta FA --ip-filter-mode filter`; prime ran "
    "`--peak-strategy clip_seeded --genome-fasta FA` — two flags fewer, none added, because "
    "`--ip-filter-default auto` and `--ip-filter-mode auto→filter` now supply them. "
    "HEADLINE: prime's default reproduces v2's manuscript arm byte for byte (md5-identical pas.bed, pas_tier1.bed, "
    "pas_tier2.bed and precision default on all four libraries), so the pre-registered PRIMARY criterion of 24 §3.1 "
    "(ΔP ≥ −0.005 AND ΔR ≥ +0.010 AND ΔF1 > 0, on all three gate datasets at once) **FAILS with every Δ exactly "
    f"0.000000** — no new accuracy on these libraries, and no regression either. What prime's default does change is "
    f"what a user gets by typing nothing: PBMC P@100 {p_nf:.4f} → {p_pr:.4f} ({(p_pr/p_nf-1)*100:+.1f} %) for R_det "
    f"{r_nf:.4f} → {r_pr:.4f} ({(r_pr/r_nf-1)*100:+.1f} %), Kinnex ≥5-UMI concordance {_k5_nf:.4f} → {_k5_pr:.4f} and "
    f"internal-priming decoy rate {_dc_nf:.4f} → {_dc_pr:.4f}; prime's set is a strict SUBSET of v2's flagless set "
    f"(0 added, 68,855 removed). Precision default P@100 v1 {_v1p[0]:.3f}/{_v1p[1]:.3f}/{_v1p[2]:.3f} → v2 = prime "
    f"{_v2p[0]:.3f}/{_v2p[1]:.3f}/{_v2p[2]:.3f} (gate P ≥ 0.50 PASS on all three); v1 → v2 is a BUG FIX, not an "
    "accuracy gain, and PBMC precision goes slightly DOWN with it. PROVENANCE: shipped / v1 / v2 are the runs "
    "recorded and verified in 15, 19 and 26, re-scored here with the same scorer (240 recomputed-vs-archived "
    "comparisons, 0 mismatches); prime's four arms are NEW runs from an unmerged branch, checked by this benchmark "
    "only and not yet through an independent manuscript verification pass. Atlas precision is agreement with "
    "PolyASite 2.0, not ground truth; Kinnex truth is a different donor and the mice have none; shipped and v1 were "
    "never run on pbmc4k and no such run was fabricated; wall time is not comparable across version sets (f)."
)

# a panel title must fit inside its own panel -- measured on the real renderer,
# not estimated, so a longer title in a future edit fails here instead of clipping
fig.canvas.draw()
_rend = fig.canvas.get_renderer()
for _ax, _t in TITLES:
    _w = _ax.title.get_window_extent(renderer=_rend).width
    _aw = _ax.get_window_extent(renderer=_rend).width
    assert _w <= _aw, f"panel title {_t!r} is {_w:.0f} px wide in a {_aw:.0f} px panel"
# a rotated y-label must fit inside its own axes -- panel f's two sub-axes are short and
# a label that overruns one of them lands on the other's label
for _ax in (axA, axB, axC, axD, axE, axF, axFw, axG):
    _lab = _ax.yaxis.label
    if not _lab.get_text():
        continue
    _h = _lab.get_window_extent(renderer=_rend).height
    _ah = _ax.get_window_extent(renderer=_rend).height
    assert _h <= _ah, (f"y-label {_lab.get_text()!r} is {_h:.0f} px tall in a {_ah:.0f} px axes")

# a panel's title + subtitle block grows UPWARD from its axes: it must not climb into the
# panel above it.  Panel g sits under panel e and panel f under panel e too, so a longer
# subtitle in a future edit fails here instead of printing over panel e's tick labels.
for _ax, _above in ((axG, axE), (axFw, axE)):
    _blk = _ax.title.get_window_extent(renderer=_rend)
    _abv = _above.get_window_extent(renderer=_rend)
    assert _blk.y1 <= _abv.y0, (
        f"panel {_ax.get_label() or _ax.title.get_text()[:1]!r} title/subtitle block reaches "
        f"y={_blk.y1:.0f} px, above panel e's axes bottom y={_abv.y0:.0f} px -- shorten the subtitle")

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
fig.savefig(OUTDIR / f"{NAME}.png", dpi=600)
print("wrote", OUTDIR / f"{NAME}.png")

# ---------------------------------------------------------------------------
# TSVs -- every plotted value
# ---------------------------------------------------------------------------
pl = pd.DataFrame(plotted).drop_duplicates(
    subset=["panel", "dataset", "version", "arm", "ranking", "metric", "note"])
pl = pl[["panel", "dataset", "version", "arm", "ranking", "metric", "value", "note", "source_file"]]
pl = pl.sort_values(["panel", "dataset", "version", "arm", "metric"])
p = OUTDIR / f"{NAME}.tsv"; pl.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p, len(pl), "rows")

p = OUTDIR / f"{NAME}_compute.tsv"
pd.DataFrame(compute_rows).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

crit = pd.DataFrame(CRIT)
crit["criterion"] = ("24 section 3.1 PRIMARY: dP>=-0.005 AND dR>=+0.010 AND dF1>0, "
                     "on all three gate datasets simultaneously")
crit["overall"] = CRIT_OVERALL
crit["baseline"] = "v2 precision_default"
crit["test"] = "prime precision_default"
crit["source_file"] = str(FOURWAY)
p = OUTDIR / f"{NAME}_criterion.tsv"; crit.to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

p = OUTDIR / f"{NAME}_reference_lines.tsv"
pd.DataFrame(ref_lines).drop_duplicates().to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption + the paragraph for manuscript/05_figure_index.md
# ---------------------------------------------------------------------------
index_para = f"""## figS12_versions — Fig S12, the four-way version comparison: shipped → v1 → v2 → prime (2026-08-22; prime EXPLORATORY, unmerged)

Seven panels, four libraries, one scorer: (a) progression on the full call set — the arm every version emits —
and (b) progression at the pre-registered precision default, which `shipped` cannot produce at all (its BED score
column is 0 for every call); (c, d) the four versions in the precision/recall plane against the competitor panel,
with F1 isolines, the pre-registered gate P@100 ≥ 0.50 and the genic-shuffle nulls, PBMC and testis; (e) matched
call count under each version's own ranking, two keys (UMI depth for all four, clip molecules for v1/v2/prime only)
with the tie-break exposure drawn; (f) compute — PBMC peak RSS {_r1:,.1f} → {_r2:,.1f} GB across v1 → v2
({_r1/_r2:.1f}×; {_r1n/_r2n:.1f}× on the flagless pair), wall time as a run record only; (g) atlas precision beside
Kinnex long-read concordance and the internal-priming decoy rate. **Headline: prime's default reproduces v2's
manuscript arm byte for byte**, so the pre-registered PRIMARY criterion of
[24](24_prime_preregistration.md) §3.1 **FAILS with every Δ exactly 0.000000** — no new accuracy on these four
libraries, and no regression. What prime changes is the *flagless* default: typed with no behaviour flag, PBMC
P@100 {p_nf:.4f} → {p_pr:.4f} for R_det {r_nf:.4f} → {r_pr:.4f}, Kinnex ≥5-UMI concordance {_k5_nf:.4f} →
{_k5_pr:.4f}, decoy rate {_dc_nf:.4f} → {_dc_pr:.4f} — the one axis in the whole benchmark where the curated atlas
and the long reads move together. **Caveats that travel with it:** v2 is the manuscript's record and prime is an
unmerged development branch whose every number is EXPLORATORY (24 §3.4 — adoption needs its own pre-registration,
a fresh independent verification pass and a statement of which figures move); shipped/v1/v2 arms are the runs
verified in 15, 19 and 26 and were re-scored here (240 comparisons, 0 mismatches) while prime's are new; v1 → v2 is
a bug fix, not an accuracy gain, and PBMC precision falls {abs(_dP['pbmc10k']):.4f} across it; atlas precision is
atlas *agreement*; Kinnex truth is a different donor and the mice have no long-read truth; `shipped` cannot be
ranked by clip molecules; shipped and v1 were never run on pbmc4k; wall times are not comparable across version
sets. Full caption: `figures/figS12_versions.caption.md`; benchmark record `results/prime_bench/README.md`;
pre-registration [24_prime_preregistration.md](24_prime_preregistration.md)."""

cap_md = f"""# Fig S12 — `figS12_versions` sidecar (generated by `scripts/manuscript_figures/figS12_versions.py`; renamed 2026-09-02, see FIGURE_MAP.tsv; single source of the legend)

## Legend

{legend_txt}

**Panel a** — progression on the **full call set**, the output every version emits: atlas-agreement precision
@100 bp (solid, filled circles) and detected-gene recall @100 bp (dashed, open squares); colour = dataset. PBMC
precision {val('pbmc10k','shipped','full','P@100')[0]:.4f} (shipped) → {val('pbmc10k','v1','full','P@100')[0]:.4f}
(v1) → {val('pbmc10k','v2','full','P@100')[0]:.4f} (v2) → identical (prime); recall
{val('pbmc10k','shipped','full','R_det@100')[0]:.4f} → {val('pbmc10k','v1','full','R_det@100')[0]:.4f} →
{val('pbmc10k','v2','full','R_det@100')[0]:.4f} → identical. The shaded v2 → prime band is exactly zero on every
metric and every dataset. These arms are **not flag-matched by choice of this figure**: shipped's code has no
internal-priming veto, v1 and v2 were given one with `--ip-filter`, prime turns it on by itself.
**Panel b** — the same at the **pre-registered precision default** (tier-1 ∩ IP-pass ∩ ≥2 clip molecules;
12 CORRECTION, 13 §1). `shipped` is greyed out because it cannot emit this output at all. The hollow markers are
v2 run with **no behaviour flag** — what a user gets by typing nothing — and the arrows to prime are the only thing
prime's default actually changes. The v1 → v2 precision deltas are printed at true size
({_dP['pbmc10k']:+.4f} PBMC, {_dP['mouse1']:+.4f} / {_dP['mouse2']:+.4f} mice): they are small, and they are drawn
small.
**Panels c, d** — the same four versions in the precision/recall plane with the competitor panel as context
(competitor points read from the Fig 2 audit TSV — same scorer, same denominators), F1_det isolines 0.1–0.4, the
pre-registered gate P@100 ≥ 0.50 (13 §1) and the 3-seed genic-shuffle nulls as ticks. Circles = full call set,
diamonds = precision default. **prime is drawn as a ring because it lands exactly on v2** — a ring is the honest
way to draw a zero difference. The hollow pink diamond in c is v2 typed with **no behaviour flag**; no arrow is
drawn from it to prime because at this scale the two points are ~21 pt apart while the markers are ~12 pt across,
so any arrow between them degenerates into a stranded arrowhead — panel b draws that move, with the numbers.
scUTRquant is catalog-based (hollow triangle, not ranked); scTail is absent (not
runnable on this BAM, R1 = 28 bp); SCAPTURE: mouse 1 plotted (mouse-2 sites scored 0.672; 15 §3, 21 §11 fix).
**Panel e** — **matched call count**: each version's full set ranked under its own shipped ranking and truncated to
a common N with the tie-break of `d3_rank_sweep.py` (descending key, then contig / start / strand / name). Two
keys, because the four versions do not share one: `shipped` has **no clip-molecule column** (BED field 5 is 0 for
all 277,164 / 45,921 / 46,672 of its calls) — exactly the limitation 25 §2.1 records for SCAPTURE — so the
four-way rows use UMI depth. The three hatched bars are the rows where more than one output's worth of tied sites
was cut by the tie-break rather than by evidence: on PBMC at N = 46,544 v1 cuts at 1 clip molecule with 117,243
sites tied and only 2,131 kept (1.8 %), while v2 and prime cut at 2 molecules and keep 16,718 of 16,718. **Do not
quote a hatched bar as like-for-like.** Same library, same budget (N = 46,544), same code, only the ranking key
differs: clip molecules {_pm:.4f} vs UMI depth {_pu:.4f} atlas precision — the Stage-1 clip evidence is what makes
the ranking work.
**Panel f** — compute, two stacked axes: the upper strip is wall time (lollipops, log), the lower is peak
RSS (bars, log). The #97 performance work is
the story: PBMC peak RSS {_r1:,.1f} → {_r2:,.1f} GB ({_r1/_r2:.1f}×) on the arms shown and {_r1n:,.1f} →
{_r2n:,.1f} GB ({_r1n/_r2n:.1f}×) on the flagless pair. **Wall time is a run record, not a benchmark**: the v2 arms
ran concurrently while prime's ran one at a time with OMP/MKL/OPENBLAS/NUMBA pinned to 1 thread. Peak RSS is the
quotable figure; the one movement outside noise from v2 to prime is mouse 1, {_m1v2:.2f} → {_m1pr:.2f} GB
({_m1pr/_m1v2:.2f}×), inside the pre-registered 1.5× guard rail.
**Panel g** — the reliability axis, PBMC 10k v3 ≥2-molecule output: atlas-agreement precision beside Kinnex x3p
long-read concordance at ≥5 and ≥20 UMI and the fraction of calls sitting at a Kinnex internal-priming decoy
(lower is better). v1 → v2 moves both truths slightly **down** together (atlas {_dP['pbmc10k']:+.4f}, Kinnex ≥5
{val('pbmc10k','v2','precision_default','kinnex_t5_frac25')[0]-val('pbmc10k','v1','precision_default','kinnex_t5_frac25')[0]:+.4f}) —
the fix was for correctness, not for score. The v2-no-flags → prime step moves all four read-outs the right way at
once, and it is the only change in this benchmark that does. **The four bars are not at equal call count**
(n = {_gn['v1']:,} / {_gn['v2']:,} / {_gn['prime']:,} / {_gn['v2_noIP']:,}) and a fraction-of-calls read-out falls
with n, so the flagless bar is read against three matched-budget controls, all on the same scorer
(`results/prime_bench/adversarial/veto_matched_pbmc10k.tsv`): **at matched call count** (v2-no-flags truncated by
clip-molecule rank to prime's own N = {_gn['prime']:,}) it reads atlas {_mc['P_atlas@100']:.4f} / Kinnex ≥5
{_mc['kinnex_t5@25']:.4f} / ≥20 {_mc['kinnex_t20@25']:.4f} / decoy {_mc['kinnex_decoy@25']:.4f} against prime's
{val('pbmc10k','prime','precision_default','P@100')[0]:.4f} / {_k5_pr:.4f} /
{val('pbmc10k','prime','precision_default','kinnex_t20_frac25')[0]:.4f} / {_dc_pr:.4f} — prime is ahead on **all
four**, and on recall too ({veto('matched_calls','v2_flagless_topN46524','R_det@100')[0]:.4f} vs
{val('pbmc10k','prime','precision_default','R_det@100')[0]:.4f}); **at matched atlas precision** (v2-no-flags
truncated to P@100 {veto('matched_precision','v2_flagless_topN40687','P_atlas@100')[0]:.6f}, prime
{val('pbmc10k','prime','precision_default','P@100')[0]:.6f}) prime holds
{val('pbmc10k','prime','precision_default','R_det@100')[0]:.4f} recall against
{_mp_r:.4f} (+{(val('pbmc10k','prime','precision_default','R_det@100')[0]/_mp_r-1)*100:.1f} %) and wins Kinnex ≥5
and the decoy rate, while Kinnex ≥20 goes the other way by
{veto('matched_precision','v2_flagless_topN40687','kinnex_t20@25')[0]-val('pbmc10k','prime','precision_default','kinnex_t20_frac25')[0]:+.4f}
— the single sign disagreement anywhere in this comparison, and it disappears at matched n; **at matched atlas
recall** (v2-no-flags at N = {int(veto('matched_recall','v2_flagless_topN51900','n')[0]):,}) prime leads on precision
{val('pbmc10k','prime','precision_default','P@100')[0]:.4f} vs
{veto('matched_recall','v2_flagless_topN51900','P_atlas@100')[0]:.4f} and on every Kinnex read-out. One place the
veto does **not** win: forced to v2-no-flags' own budget of {_gn['v2_noIP']:,} calls prime must reach into its
1-molecule tier and lands slightly behind
({veto('crossing_check','prime_topN62110','P_atlas@100')[0]:.4f} /
{veto('crossing_check','prime_topN62110','R_det@100')[0]:.4f} vs
{veto('crossing_check','v2_flagless_ge2mol','P_atlas@100')[0]:.4f} /
{veto('crossing_check','v2_flagless_ge2mol','R_det@100')[0]:.4f}) — that row is
{100*(1-veto('crossing_check','prime_topN62110','n_kept_at_cut')[0]/veto('crossing_check','prime_topN62110','n_tied_at_cut')[0]):.0f} %
tie-decided at the cut and is not like-for-like. shipped has no ≥2-molecule output at all; its full call set
reaches Kinnex ≥5-UMI {_sh_t5:.4f} against {val('pbmc10k','prime','full','kinnex_t5_frac25')[0]:.4f} for prime's
full set — quoted for scale only, not drawn as a bar.

**Must travel with it.**
- **v2 is the manuscript's record; prime is not.** `peakAtail-prime` (6954082d) is unmerged. 24 §3.4 is explicit
  that prime cannot retroactively become "the pre-registered default": adoption needs its own pre-registration
  written before an adoption run, a fresh independent verification pass, and a statement of which figures move.
  Until those exist the reported default is v2 and **nothing on this figure changes a manuscript number**.
- **The pre-registered criterion FAILED, and it failed at exactly zero.** 24 §3.1 asked for ΔR ≥ +0.010 and
  ΔF1 > 0 against **v2's default output**, which the pre-registration defines as the flagged arm. Every Δ is
  0.000000 because prime's default output *is* that file. A FAIL here means **no new accuracy on these four
  libraries**, not a regression. Do not restate the criterion to make it pass.
- **The flagless comparison (panels b, c, g) is NOT pre-registered** and exists on PBMC 10k v3 only:
  `stage2_final_launch.sh` always passed `--ip-filter` on the mouse arms, so no v2 no-IP mouse run exists and none
  was fabricated. Quote it as "what changes for a user who types nothing", never as a gain over the manuscript arm.
- **prime = v2 is a byte-for-byte identity, not an approximation.** md5-identical `pas.bed`, `pas_tier1.bed`,
  `pas_tier2.bed` and `pas_PRESPEC_precision_default.bed` on all four libraries (PBMC md5
  2a7f9f76d6a70ba2daacc8aab0c2817d); the whole 57-file PBMC run tree differs only in sidecar enrichment
  (`pas_support.tsv` gains 26 columns; `cut -f1-7` of prime's file is byte-identical to the whole of v2's).
- **v1 → v2 is a bug fix, not an accuracy gain.** The corrected minus-strand IP window (#96) moves the PBMC
  precision default **down** {val('pbmc10k','v1','precision_default','P@100')[0]:.4f} →
  {val('pbmc10k','v2','precision_default','P@100')[0]:.4f} and Kinnex ≥5-UMI concordance down
  {val('pbmc10k','v1','precision_default','kinnex_t5_frac25')[0]:.4f} →
  {val('pbmc10k','v2','precision_default','kinnex_t5_frac25')[0]:.4f}, while both mice move up
  {_dP['mouse1']:+.4f} / {_dP['mouse2']:+.4f}. Draw and quote it as a correctness change.
- **`shipped` cannot be ranked by clip molecules** (BED field 5 = 0 everywhere), so it appears in the matched
  comparison only under the UMI key and has no ≥2-molecule arm anywhere on this figure. That is a property of the
  pre-Stage-1 caller, not a gap in the benchmark.
- **Atlas precision is atlas *agreement*, not truth** — an atlas-novel true site counts as a false positive. The
  atlas-independent read-out is Kinnex, whose own limits stand (16): a **different donor** from both PBMC
  libraries, site-level metric. The two mouse libraries have no long-read truth, so panel g is human-only.
- **shipped and v1 were never run on pbmc4k** and no such run was fabricated; that library compares v2 and prime
  only (and they are identical). pbmc4k's protocol, denominator and gate were fixed by 26 before it was scored.
- **Wall time is not comparable across version sets** (v2 concurrent, prime sequential with BLAS pinned to one
  thread); peak RSS is. The two things prime ships **off** — `--read-geometry true` and `--pas-score calibrated` —
  were not exercised here and are not on this figure.

## Provenance

Sources: `results/prime_bench/fourway.tsv` + `results/prime_bench/README.md` (the four-way benchmark record; every
row carries its own source file) for panels a, b, e, f, g and the PeakATail points of c and d;
`results/figures/manuscript/fig2_accuracy.tsv` (the Fig 2 audit TSV, same scorer and denominators) for the
competitor points of c and d. Verified records behind the non-prime arms: `manuscript/15_final_gate.md` (v1,
SOUND), `manuscript/19_final_gate_v2.md` + `results/benchmark_tools/final_v2_verify/VERIFIED_v2.md` (v2, FIXED),
`manuscript/26` (pbmc4k). Pre-registration and the status of every prime number:
`manuscript/24_prime_preregistration.md` §3.1, §3.4. Every plotted value:
`results/figures/manuscript/figS12_versions.tsv`; panel f `figS12_versions_compute.tsv`; the criterion
evaluation `figS12_versions_criterion.tsv`; gates and isolines `figS12_versions_reference_lines.tsv`.
Palette: Okabe-Ito (`manuscript/figures/README.md`). Script: `scripts/manuscript_figures/figS12_versions.py`;
rendered at PNG 600 dpi / PDF fonttype 42.

### Paragraph for `manuscript/05_figure_index.md` (this script does NOT write it there)

{index_para}
"""
# DESIGN_DIRECTIVES.md item 2, no-loss: the panel subtitles that used to sit under the
# panel titles are written into the LEGEND section verbatim -- its only home is this
# sidecar, and the no-loss rule puts a removed on-figure sentence in the Legend.
_pn = ("**Panel notes** (each was a subtitle line on the image until the 2026-09-03 design "
       "pass; reproduced verbatim, nothing dropped):\n\n")
for _letter, _title, _sub in SUBTITLES:
    _pn += f"- **{_title.strip()}** — {_sub}\n"
assert "\n## Provenance\n" in cap_md, "caption template lost its Provenance heading"
cap_md = cap_md.replace("\n## Provenance\n", "\n" + _pn + "\n## Provenance\n", 1)
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and all 19 on-figure
annotation call sites that sat below the 6 pt floor — the smallest were 4.7 pt — were raised to it. Axis labels,
panel titles and prose tick labels are sentence-cased through `_pubstyle.sentence_case()`, canonical identifiers
preserved and the lower-case panel letters kept. Every panel subtitle left the image for the **Panel notes**
list above (that is the whole of the no-loss move). The iso-F1 contour labels in **c** and **d** are offset off
their own curves, which used to run through them. No panel, number or audit TSV changed; all four TSVs
regenerate byte-identical.
"""
p = FIGDIR / f"{NAME}.caption.md"; p.write_text(cap_md); print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: nothing may be clipped; the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"\n24 section 3.1 verdict: {CRIT_OVERALL} (dP/dR/dF1 = 0.000000 on all three gate datasets)")
print(f"v2 = prime precision default P@100 {_v2p[0]:.4f} / {_v2p[1]:.4f} / {_v2p[2]:.4f}")
print(f"flagless PBMC P@100 {p_nf:.4f} -> {p_pr:.4f}, R_det {r_nf:.4f} -> {r_pr:.4f}")
