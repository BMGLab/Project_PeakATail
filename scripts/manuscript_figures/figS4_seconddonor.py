#!/usr/bin/env python3
"""
figS4_seconddonor.py -- manuscript Fig S4 ("figS4_seconddonor"): the second
human donor validation of manuscript/26_second_donor_preregistration.md
(pre-registered 2026-08-21 23:47-23:53, run 2026-08-22 00:07-00:23, verifier
verdict FIXED).

SOURCES OF TRUTH (FIGURES_MANIFEST.md row S4: "26 (same scorer/nulls)")
  * pbmc4k arms: the run's own score_tool.py TSVs under
    results/benchmark_tools/pbmc4k_donor2/run_ipfilt/ (the files 26 R2 quotes;
    verifier re-scored all four arms to 4 dp with an independent pipeline).
  * the other three libraries: results/figures/manuscript/fig2_accuracy.tsv --
    the verified Fig 2 audit TSV (19 FIXED), read programmatically.
  * cross-donor concordance: results/benchmark_tools/pbmc4k_donor2/
    concordance.txt (the 26 section 5 protocol output; verifier-exact).
  * matched-call-count control: the verifier's table in 26 ("Matched-call-count
    control, added by the verifier, 2026-08-22") -- cited constants with EXPECT
    asserts; there is no standalone TSV for it in the results tree.

WHAT THE FIGURE MUST SAY (the verifier's corrections travel with the numbers)
  a  The unchanged pre-registered gate P@100 >= 0.50 now holds on FOUR
     libraries; pbmc4k = 0.8279. The >=1-molecule arm is shown beside each
     default because 26 R4.3 re-establishes that >=2 molecules is a
     reliability choice, not the F1 optimum.
  b  Cross-donor concordance, BOTH directions as pre-registered (26 section 5):
     84.2% of donor-2 default calls reproduce within 100 bp in donor 1
     (81.3% within 25 bp), ~361x the genic-shuffle null; the reverse direction
     is 39.9% against its own null and its own arithmetic ceiling of 44.4%
     (20,672 / 46,524) -- the asymmetry is call-count arithmetic, not
     disagreement (39.9% is 90% of the ceiling).
  c  The verifier's matched-N caveat, drawn so it cannot be missed: at matched
     call count donor 1 DOMINATES donor 2 on both axes (P 0.886-0.907 /
     R_det 0.112-0.124 by every reduction route), so donor 2's higher headline
     precision is an operating-point effect of the fixed rule on a shallower
     v2 library. "Precision generalises; it does not improve."

NUMBER POLICY -- panels a/b are read from result TSVs with EXPECT asserts
against 26; panel c's five donor-1 reductions exist only in 26's verified
table and are cited constants (source column says so). Nothing is approximated.

OUTPUTS
  manuscript/figures/figS4_seconddonor.{png,pdf}  (PNG 300 dpi, PDF fonttype 42)
  manuscript/figures/figS4_seconddonor.caption.md
  results/figures/manuscript/figS4_seconddonor.tsv             (panel a)
  results/figures/manuscript/figS4_seconddonor_concordance.tsv (panel b)
  results/figures/manuscript/figS4_seconddonor_matchedN.tsv    (panel c)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS4_seconddonor.py
"""
import os
import textwrap
from pathlib import Path

os.environ.setdefault("LC_ALL", "C")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 7.5

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "figS4_seconddonor"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito, one colour per library (same assignment as figS1: donor 2 green)
C_D2, C_D1, C_M1, C_M2 = "#009E73", "#0072B2", "#E69F00", "#CC79A7"
GATE_P = 0.50   # 13 section 1 / 26 section 6 -- the unchanged pre-registered gate

# ---------------------------------------------------------------------------
# pbmc4k arms from the run's own score TSVs (the files 26 R2 quotes)
# ---------------------------------------------------------------------------
D2 = WD / "results/benchmark_tools/pbmc4k_donor2/run_ipfilt"


def read_score(tsv: Path) -> dict:
    df = pd.read_csv(tsv, sep="\t")
    at = df[df.cutoff_bp == 100.0]

    def one(panel, series, reference, col="value", replicate=None):
        q = at[(at.panel == panel) & (at.series == series) & (at.reference == reference)]
        if replicate is not None:
            q = q[q.replicate == replicate]
        assert len(q) == 1, (tsv, panel, series, reference, replicate, len(q))
        return float(q[col].iloc[0])

    p25 = df[(df.cutoff_bp == 25.0) & (df.panel == "precision") & (df.series == "real")
             & (df.reference == "atlas_full")].value.iloc[0]
    p10 = df[(df.cutoff_bp == 10.0) & (df.panel == "precision") & (df.series == "real")
             & (df.reference == "atlas_full")].value.iloc[0]
    nulls = [one("precision", "null_genic", "atlas_full", replicate=s) for s in (1, 2, 3)]
    return dict(
        n=int(one("precision", "real", "atlas_full", col="n_query")),
        P100=one("precision", "real", "atlas_full"), P25=float(p25), P10=float(p10),
        R_det=one("recall", "real", "atlas_detected"),
        R_det_denominator=int(one("recall", "real", "atlas_detected", col="n_query")),
        F1_det=one("f1", "real", "atlas_detected"),
        null_P_seed1=nulls[0], null_P_seed2=nulls[1], null_P_seed3=nulls[2],
        null_P_mean=float(np.mean(nulls)),
    )


d2_def = read_score(D2 / "score_pbmc4k_final_v2_ipfilt__PRESPEC_precision_default.tsv")
d2_t1 = read_score(D2 / "score_pbmc4k_final_v2_ipfilt__tier1.tsv")

# EXPECT vs manuscript/26 R2 (verifier re-scored to 4 dp)
assert d2_def["n"] == 20672 and round(d2_def["P100"], 4) == 0.8279, d2_def
assert round(d2_def["P25"], 4) == 0.7614 and round(d2_def["P10"], 4) == 0.6208, d2_def
assert round(d2_def["R_det"], 4) == 0.1104 and round(d2_def["F1_det"], 4) == 0.1948, d2_def
assert d2_def["R_det_denominator"] == 268097, d2_def   # 26 R1: 13,763 genes -> 268,097 sites
assert d2_t1["n"] == 62911 and round(d2_t1["P100"], 4) == 0.5042, d2_t1
assert round(d2_t1["R_det"], 4) == 0.1784 and round(d2_t1["F1_det"], 4) == 0.2635, d2_t1

# ---------------------------------------------------------------------------
# the other three libraries from the verified Fig 2 audit TSV (v2, 19 FIXED)
# ---------------------------------------------------------------------------
F2 = OUTDIR / "fig2_accuracy.tsv"
f2 = pd.read_csv(F2, sep="\t")


def f2_row(dataset, role, replicate="single"):
    q = f2[(f2.dataset == dataset) & (f2.role == role) & (f2.replicate == replicate)]
    assert len(q) == 1, (dataset, role, replicate, len(q))
    r = q.iloc[0]
    return dict(n=int(r.n), P100=float(r.atlas_agreement_precision_100bp),
                R_det=float(r.recall_detected_genes_100bp),
                F1_det=float(r.F1_detected_genes_100bp),
                null_P_seed1=float(r.null_P_seed1), null_P_seed2=float(r.null_P_seed2),
                null_P_seed3=float(r.null_P_seed3), null_P_mean=float(r.null_P_mean))


d1_def = f2_row("pbmc_10k_v3", "path4")
d1_t1 = f2_row("pbmc_10k_v3", "path3")
m1_def = f2_row("gse104556", "path4", "mouse1")
m2_def = f2_row("gse104556", "path4", "mouse2")
m1_t1 = f2_row("gse104556", "path3", "mouse1")
m2_t1 = f2_row("gse104556", "path3", "mouse2")

# EXPECT vs 19 section 1 / 26 R2 / 19 section 6
assert (d1_def["n"], round(d1_def["P100"], 4), round(d1_def["R_det"], 4)) == (46524, 0.7062, 0.1754)
assert (d1_t1["n"], round(d1_t1["P100"], 4)) == (167565, 0.3520)
assert (m1_def["n"], round(m1_def["P100"], 4)) == (26255, 0.7450)
assert (m2_def["n"], round(m2_def["P100"], 4)) == (26526, 0.7572)
assert round(m1_t1["P100"], 4) == 0.5686 and round(m2_t1["P100"], 4) == 0.5880

# gate holds on every default arm (and 26's tier-1 note: pbmc4k tier-1 also
# clears 0.50 numerically but is NOT the gated arm)
for d in (d2_def, d1_def, m1_def, m2_def):
    assert d["P100"] >= GATE_P, d

# ---------------------------------------------------------------------------
# cross-donor concordance (26 section 5 protocol output, verifier-exact)
# ---------------------------------------------------------------------------
CONC = WD / "results/benchmark_tools/pbmc4k_donor2/concordance.txt"
conc_rows = []
arm = None
for line in CONC.read_text().splitlines():
    if line.startswith("############ ARM="):
        arm = line.split("ARM=")[1].strip()
    elif line and not line.startswith(("#", "direction")):
        p = line.split("\t")
        conc_rows.append(dict(arm=arm, direction=p[0], n_query=int(p[1]),
                              n_le25=int(p[2]), n_le100=int(p[3]),
                              frac_le25=float(p[4]), frac_le100=float(p[5]),
                              source=str(CONC.relative_to(WD))))
conc = pd.DataFrame(conc_rows)


def cpick(arm, direction):
    q = conc[(conc.arm == arm) & (conc.direction == direction)]
    assert len(q) == 1, (arm, direction)
    return q.iloc[0]


fwd = cpick("PRESPEC_precision_default", "pbmc4k->pbmc10k_REAL")
rev = cpick("PRESPEC_precision_default", "pbmc10k->pbmc4k_REAL")
fwd_null = conc[(conc.arm == "PRESPEC_precision_default")
                & conc.direction.str.startswith("pbmc4k->pbmc10k_null")]
rev_null = conc[(conc.arm == "PRESPEC_precision_default")
                & conc.direction.str.startswith("pbmc10k->pbmc4k_null")]
assert len(fwd_null) == 3 and len(rev_null) == 3
fwd_null100 = float(fwd_null.frac_le100.mean())
rev_null100 = float(rev_null.frac_le100.mean())

# EXPECT vs 26 R3 (verifier recomputed both directions exactly)
assert (fwd.n_query, round(fwd.frac_le25, 4), round(fwd.frac_le100, 4)) == (20672, 0.8132, 0.8423)
assert (rev.n_query, round(rev.frac_le25, 4), round(rev.frac_le100, 4)) == (46524, 0.3613, 0.3993)
FOLD = fwd.frac_le100 / fwd_null100
assert round(FOLD) == 361, FOLD                       # "~361x the genic-shuffle null"
CEIL = d2_def["n"] / d1_def["n"]                      # 20,672 / 46,524
assert round(CEIL * 100, 1) == 44.4, CEIL             # 26 V: reverse-direction ceiling
assert 0.89 < rev.frac_le100 / CEIL < 0.91            # "39.9% is 90% of that ceiling"

# ---------------------------------------------------------------------------
# matched-call-count control -- 26's verified table (cited constants; no TSV
# exists in the results tree for this control)
# ---------------------------------------------------------------------------
MN_SRC = "manuscript/26_second_donor_preregistration.md, 'Matched-call-count control (added by the verifier, 2026-08-22)'"
matched = pd.DataFrame([
    dict(reduction="top 20,672 by clip-molecule rank (coordinate tie-break)", n=20672, P100=0.9049, R_det=0.1135),
    dict(reduction="top 20,672 by clip-molecule rank (reverse tie-break)",    n=20672, P100=0.9049, R_det=0.1134),
    dict(reduction=">=4 molecules", n=23679, P100=0.8856, R_det=0.1241),
    dict(reduction=">=5 molecules", n=20320, P100=0.9074, R_det=0.1122),
])
matched["arm"] = "donor 1 default, reduced to ~donor 2's call count"
matched["source"] = MN_SRC
# the caveat the panel exists to show: donor 1 dominates on BOTH axes at matched N
for r in matched.itertuples():
    assert r.P100 > d2_def["P100"] and r.R_det > d2_def["R_det"], r

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8.3, 9.2))
gs = fig.add_gridspec(2, 2, height_ratios=[0.80, 1.0],
                      left=0.105, right=0.975, top=0.955, bottom=0.230,
                      hspace=0.52, wspace=0.34)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=6.8)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


# ---- panel a: the gate on four libraries, default + >=1-molecule arm -------
LIBS = [
    ("pbmc4k\n(donor 2, v2 chem.)", C_D2, d2_def, d2_t1,
     "results/benchmark_tools/pbmc4k_donor2/run_ipfilt/score_*.tsv"),
    ("PBMC 10k v3\n(donor 1, v3 chem.)", C_D1, d1_def, d1_t1,
     "results/figures/manuscript/fig2_accuracy.tsv (19 FIXED)"),
    ("testis mouse 1\n(GSE104556)", C_M1, m1_def, m1_t1,
     "results/figures/manuscript/fig2_accuracy.tsv (19 FIXED)"),
    ("testis mouse 2\n(GSE104556)", C_M2, m2_def, m2_t1,
     "results/figures/manuscript/fig2_accuracy.tsv (19 FIXED)"),
]
tsvA = []
xs = np.arange(len(LIBS)) * 1.0
bw = 0.32
for x, (lab, c, dflt, t1, src) in zip(xs, LIBS):
    axA.bar(x - bw / 2 - 0.02, dflt["P100"], width=bw, color=c, edgecolor=c, zorder=3)
    axA.bar(x + bw / 2 + 0.02, t1["P100"], width=bw, color="white", edgecolor=c,
            hatch="////", linewidth=1.0, zorder=3)
    axA.text(x - bw / 2 - 0.02, dflt["P100"] + 0.012, f"{dflt['P100']:.4f}\nn={dflt['n']:,}",
             ha="center", va="bottom", fontsize=5.9, color=INK, linespacing=1.15)
    axA.text(x + bw / 2 + 0.02, t1["P100"] + 0.012, f"{t1['P100']:.4f}\nn={t1['n']:,}",
             ha="center", va="bottom", fontsize=5.9, color=MUTED, linespacing=1.15)
    for dx, arm_name, d in ((-bw / 2 - 0.02, "pre-registered default (tier-1, IP, >=2 mol)", dflt),
                            (+bw / 2 + 0.02, ">=1 molecule arm (tier-1, IP; not gated)", t1)):
        axA.plot([x + dx - bw * 0.4, x + dx + bw * 0.4], [d["null_P_mean"]] * 2,
                 color=INK, lw=1.0, zorder=4)
        tsvA.append(dict(panel="a", library=lab.replace("\n", " "), arm=arm_name,
                         n=d["n"], P100=d["P100"], R_det=d["R_det"], F1_det=d["F1_det"],
                         null_P_seed1=d["null_P_seed1"], null_P_seed2=d["null_P_seed2"],
                         null_P_seed3=d["null_P_seed3"], null_P_mean=d["null_P_mean"],
                         gate="P@100 >= 0.50" if "default" in arm_name else "not gated",
                         verdict=("PASS" if d["P100"] >= GATE_P else "FAIL") if "default" in arm_name else "",
                         source=src))
axA.axhline(GATE_P, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
# gate label sits in the free gap between the donor-1 and mouse-1 groups
axA.text(1.5, GATE_P + 0.025, "pre-registered gate\nP@100 ≥ 0.50\n(13 §1; unchanged, 26 §6)",
         ha="center", va="bottom", fontsize=5.9, color=INK, linespacing=1.2,
         bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
axA.annotate("PASS on a 4th library --\nbut a more conservative\noperating point (see c)",
             xy=(xs[0] - bw / 2 - 0.02, d2_def["P100"]), xytext=(xs[0] + 0.32, 0.985),
             fontsize=6.0, color=C_D2, ha="left", va="top",
             arrowprops=dict(arrowstyle="-", color=C_D2, lw=0.7))
axA.text(0.995, 0.985, "solid = pre-registered default (gated)   hatched = ≥1-molecule arm (not gated)\n"
         "— = 3-seed genic-shuffle null mean", transform=axA.transAxes,
         ha="right", va="top", fontsize=5.8, color=MUTED, linespacing=1.35)
axA.set_xticks(xs)
axA.set_xticklabels([l[0].replace("\n", "\n") for l in LIBS], fontsize=6.6)
axA.set_xlim(-0.62, xs[-1] + 0.62)
axA.set_ylim(0, 1.0)
axA.set_ylabel("atlas-agreement precision @100 bp")
axA.set_title("a   The unchanged gate P@100 ≥ 0.50, now on four libraries",
              loc="left", fontweight="bold", fontsize=7.6)
style(axA)

# ---- panel b: cross-donor concordance, both directions ---------------------
bars = [
    ("d2→d1 ≤100 bp", fwd.frac_le100, C_D2, fwd_null100, fwd.n_query),
    ("d2→d1 ≤25 bp", fwd.frac_le25, C_D2, float(fwd_null.frac_le25.mean()), fwd.n_query),
    ("d1→d2 ≤100 bp", rev.frac_le100, C_D1, rev_null100, rev.n_query),
    ("d1→d2 ≤25 bp", rev.frac_le25, C_D1, float(rev_null.frac_le25.mean()), rev.n_query),
]
ys = np.arange(len(bars))[::-1] * 0.9
for (lab, v, c, nl, nq), y in zip(bars, ys):
    axB.barh(y, v, height=0.52, color=c, edgecolor=c, zorder=3,
             alpha=1.0 if "100" in lab else 0.55)
    axB.text(v + 0.015, y, f"{v * 100:.1f}%", va="center", fontsize=6.4, color=INK, zorder=5)
    axB.plot([nl] * 2, [y - 0.30, y + 0.30], color=INK, lw=1.2, zorder=4)
axB.set_yticks(ys)
axB.set_yticklabels([b[0] for b in bars], fontsize=6.3)
# the reverse direction's arithmetic ceiling, drawn only across its two bars
axB.plot([CEIL] * 2, [ys[2] - 0.45, ys[3] + 0.45], color=C_D1, lw=1.0, ls=(0, (4, 2)), zorder=4)
axB.text(CEIL + 0.015, (ys[2] + ys[3]) / 2, "arithmetic ceiling 44.4%\n(20,672 / 46,524):\n"
         "39.9% = 90% of it.\nAsymmetry is call-count\narithmetic, not disagreement",
         fontsize=5.8, color=C_D1, va="center", ha="left", linespacing=1.25)
axB.annotate(f"{FOLD:.0f}× the genic-\nshuffle null (|)",
             xy=(fwd_null100, ys[0] + 0.26), xytext=(0.22, ys[0] + 0.33),
             fontsize=5.9, color=INK, va="bottom",
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
axB.set_xlim(0, 1.0)
axB.set_ylim(ys[-1] - 0.62, ys[0] + 0.75)
axB.set_xlabel("fraction of default calls reproduced in the other donor\n(strand-matched, 26 §5 protocol)",
               fontsize=6.8)
axB.text(0.985, 0.03, "d1 = donor 1 (PBMC 10k v3, n 46,524)\nd2 = donor 2 (pbmc4k, n 20,672)",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.25)
axB.set_title("b   Cross-donor concordance,\nboth directions (pre-registered)",
              loc="left", fontweight="bold", fontsize=7.6)
style(axB)
axB.grid(True, axis="x", color=GRID, lw=0.5, alpha=0.7)
axB.grid(False, axis="y")

# ---- panel c: the matched-N caveat, drawn -----------------------------------
# shade the region that dominates donor 2 on both axes
axC.add_patch(matplotlib.patches.Rectangle((d2_def["R_det"], d2_def["P100"]),
              0.20 - d2_def["R_det"], 1.0 - d2_def["P100"],
              facecolor="#DCE9F2", edgecolor="none", zorder=0))
axC.text(0.1975, 0.845, "region that dominates\ndonor 2 on both axes", fontsize=5.9, color=C_D1,
         ha="right", va="bottom", style="italic")
# donor 1 full default, for context; one arrow to its matched-N reductions
axC.plot(d1_def["R_det"], d1_def["P100"], marker="D", ms=7, mfc="white", mec=C_D1,
         mew=1.4, ls="none", zorder=6)
axC.annotate("donor 1 default\n(n = 46,524)", (d1_def["R_det"], d1_def["P100"]),
             xytext=(2, 10), textcoords="offset points", ha="left", va="bottom",
             fontsize=6.1, color=C_D1)
axC.annotate("", xy=(0.1205, 0.893), xytext=(d1_def["R_det"], d1_def["P100"]),
             arrowprops=dict(arrowstyle="-|>", color=C_D1, lw=1.0, alpha=0.7,
                             shrinkA=8, shrinkB=10), zorder=2)
axC.text(0.153, 0.790, "same rule, cut to\ndonor 2's call count\n(4 routes, 26)", fontsize=5.9,
         color=C_D1, ha="left", va="top", linespacing=1.25)
mk = dict(marker="o", ms=6, mfc=C_D1, mec=C_D1, ls="none", zorder=6)
lab_off = {">=4 molecules": (7, -1, "left", "center"),
           ">=5 molecules": (-4, 8, "right", "bottom"),
           "top 20,672 by clip-molecule rank (coordinate tie-break)": (7, 6, "left", "center"),
           "top 20,672 by clip-molecule rank (reverse tie-break)": (7, -5, "left", "center")}
short_lab = {">=4 molecules": "≥4 mol (n 23,679)",
             ">=5 molecules": "≥5 mol (n 20,320)",
             "top 20,672 by clip-molecule rank (coordinate tie-break)": "top 20,672 by mol. rank",
             "top 20,672 by clip-molecule rank (reverse tie-break)": "(reverse tie-break)"}
for r in matched.itertuples():
    axC.plot(r.R_det, r.P100, **mk)
    dx, dy, ha, va = lab_off[r.reduction]
    axC.annotate(short_lab[r.reduction], (r.R_det, r.P100), xytext=(dx, dy),
                 textcoords="offset points", fontsize=5.8, color=C_D1, ha=ha, va=va)
axC.plot(d2_def["R_det"], d2_def["P100"], marker="D", ms=8, mfc=C_D2, mec=C_D2, ls="none", zorder=7)
axC.annotate("donor 2 default\n(n = 20,672)\nheadline 0.8279", (d2_def["R_det"], d2_def["P100"]),
             xytext=(-8, -10), textcoords="offset points", ha="left", va="top",
             fontsize=6.1, color=C_D2, fontweight="bold", linespacing=1.25)
axC.set_xlim(0.10, 0.20)
axC.set_ylim(0.68, 0.95)
axC.set_xlabel("detected-gene recall R_det @100 bp", fontsize=6.8)
axC.set_ylabel("atlas-agreement precision @100 bp")
axC.set_title("c   The matched-N caveat: donor 2's\nheadline is an operating-point effect",
              loc="left", fontweight="bold", fontsize=7.6)
axC.text(0.1025, 0.684, "matched-call-count control: 26 (verifier,\n2026-08-22); every route lands up-and-right",
         ha="left", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.3)
style(axC)

# ---------------------------------------------------------------------------
# caption (footer)
# ---------------------------------------------------------------------------
caption = (
    "Fig S4 | Second-donor validation (pre-registered in 26 on 2026-08-21, BEFORE any pbmc4k number "
    "existed; run 2026-08-22; verifier verdict FIXED). (a) The 10x public pbmc4k library (v2 chemistry, "
    "CellRanger 2.1.0, 4,340 cells) scored by the UNCHANGED pre-registered precision default under the "
    "identical frozen code (9dfdefb3), scorer, atlas and nulls: 107 of 108 recorded parameters identical, "
    "the one change being --seq-len 91→98, a library descriptor. Gate P@100 ≥ 0.50: PASS at 0.8279 "
    "(38× its null) — the gate now holds on four libraries. The ≥1-molecule arm (hatched) again has the "
    "higher F1_det (0.2635 vs 0.1948): ≥2 molecules is a reliability choice, not the F1 optimum. "
    "(b) Cross-donor concordance per the pre-registered §5 protocol, both directions, never one alone: "
    "84.2% of donor-2 default calls have a strand-matched donor-1 call within 100 bp (81.3% within 25 bp), "
    "~361× the genic-shuffle null; the reverse direction is 39.9% (36.1% at 25 bp) against its own null "
    "and an arithmetic ceiling of 44.4% — donor 1 emits 46,524 default calls to donor 2's 20,672, so the "
    "asymmetry is call-count arithmetic, not disagreement. (c) The verifier's matched-N control, drawn so it "
    "cannot be buried: cut donor 1's default to donor 2's call count by ANY route (its own clip-molecule rank "
    "under both tie-breaks, ≥4 or ≥5 molecules) and donor 1 lands at P@100 0.886–0.907 with R_det "
    "0.112–0.124 — dominating donor 2 on both axes. Donor 2's higher headline precision is therefore an "
    "operating-point effect: the fixed rule lands at a more conservative point on the shallower v2 library. "
    "Precision GENERALISES; it does not improve. Scope (26): a second library, chemistry and CellRanger "
    "version — only presumptively a second individual (10x publishes no donor id). Recall is genuinely lower "
    "on donor 2 (R_det 0.1104 vs 0.1754; survives every denominator swap). "
    "Values: results/figures/manuscript/figS4_seconddonor*.tsv."
)
_cap = textwrap.fill(caption, 158)
_n_lines = _cap.count("\n") + 1
assert _n_lines <= 13, f"caption is {_n_lines} wrapped lines; collides with panel x labels"
fig.text(0.03, 0.012, _cap, fontsize=5.8, color=INK, va="bottom", ha="left", linespacing=1.35)

for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=300); print("wrote", p)

# ---------------------------------------------------------------------------
# audit TSVs
# ---------------------------------------------------------------------------
p = OUTDIR / f"{NAME}.tsv"
pd.DataFrame(tsvA).to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

conc_out = conc.copy()
conc_out["panel"] = "b"
conc_out["note"] = np.where(conc_out.direction.str.contains("null"),
                            "genic-shuffle null (same seeds/inclusion file as score_tool.py)", "")
p = OUTDIR / f"{NAME}_concordance.tsv"
conc_out.to_csv(p, sep="\t", index=False, float_format="%.6f"); print("wrote", p)

matched_out = matched.copy()
matched_out["panel"] = "c"
ctx = pd.DataFrame([
    dict(reduction="donor 2 default (comparison target)", n=d2_def["n"], P100=d2_def["P100"],
         R_det=d2_def["R_det"], arm="pbmc4k pre-registered default", panel="c",
         source="run_ipfilt/score_pbmc4k_final_v2_ipfilt__PRESPEC_precision_default.tsv"),
    dict(reduction="donor 1 default, full (context)", n=d1_def["n"], P100=d1_def["P100"],
         R_det=d1_def["R_det"], arm="pbmc_10k_v3 pre-registered default", panel="c",
         source="results/figures/manuscript/fig2_accuracy.tsv"),
])
p = OUTDIR / f"{NAME}_matchedN.tsv"
pd.concat([matched_out, ctx], ignore_index=True).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S4 — `figS4_seconddonor` caption (generated by `scripts/manuscript_figures/figS4_seconddonor.py`)

{caption}

**Pre-registration provenance (26, verifier verdict FIXED).** Written 2026-08-21 23:47–23:53, before any
accuracy number from pbmc4k existed; the driver `run_donor2.sh` was frozen at 23:48:30, 19 minutes before the
caller started (2026-08-22 00:07:45). The gate is quoted verbatim from the committed, already-verified 19, so it
demonstrably predates donor 2. Nothing was retuned — measured, not promised: a field-by-field diff of the two
runs' `run_config.json` finds 107/108 recorded parameters identical, the single difference being
`--seq-len 91→98` (the R2 read length, one correct value per library; `--threads` 16→8 is output-irrelevant
per the byte-identity test of `chrom_parallel.py`).

**What may and may not be claimed (26 §R4 + V).**
- The reliability claim now holds on two human 10x 3' libraries spanning two chemistries (v2/v3) and two
  CellRanger versions (2.1.0/3.0.0), plus two mice — with the caveat that 10x publishes no donor identifier
  for pbmc4k, so it is demonstrably a second *library and chemistry* and only presumptively a second
  *individual*. "Second donor" must always carry that qualification.
- **Never** "precision improves on donor 2": at matched call count donor 1 dominates on both axes (panel c).
  The defensible claim is that the precision-first behaviour is not an artefact of one donor.
- Recall is lower on donor 2 (R_det 0.1104 vs 0.1754 on the default) and this is real, not a denominator
  artefact: donor 2's detected-gene denominator is *smaller* (268,097 vs 285,136 sites, which inflates R_det)
  and its R_det is lower anyway; the gap survives scoring both arms against both denominators (26 §"Denominator
  swap"). The clip-evidence channel was 1.64× *richer* on donor 2 (3.6965% vs 2.2565% head-sample estimator),
  so channel poverty cannot be invoked as an excuse.
- No claim about v2-vs-v3 chemistry as such — that needs its own design with >1 library per chemistry.

**Panel a** — solid bars: the pre-registered default (tier-1 ∩ IP-pass ∩ ≥2 clip molecules), the arm the
gate is defined over; hatched: the ≥1-molecule sensitivity arm, reported (not gated) because it is again the
F1 optimum (26 §R4.3). Black dashes: 3-seed genic-shuffle null means (~0.022 everywhere, as they must be —
same inclusion file). pbmc4k tier-1 clears 0.50 numerically (0.5042) but is not the gated arm.
**Panel b** — both pre-registered directions with their own nulls; the ceiling explains the asymmetry
(with 20,672 donor-2 calls, at most 44.4% of donor 1's 46,524 can have a counterpart; 39.9% is 90% of that).
**Panel c** — the verifier's matched-call-count control (26): the numbers exist only in the verified document
table, cited as constants with EXPECT asserts (`figS4_seconddonor_matchedN.tsv`, source column).

Sources: `manuscript/26_second_donor_preregistration.md` (FIXED) — panel a pbmc4k bars from the run's own
`score_tool.py` TSVs (`results/benchmark_tools/pbmc4k_donor2/run_ipfilt/`), the other three libraries from the
verified `fig2_accuracy.tsv` (19 FIXED); panel b from `pbmc4k_donor2/concordance.txt`; every plotted value in
`results/figures/manuscript/figS4_seconddonor*.tsv` with a source column.
"""
p = FIGDIR / f"{NAME}.caption.md"; p.write_text(cap_md); print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: outer 8 px must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS4_seconddonor: gate PASS x4 ({d2_def['P100']:.4f} / {d1_def['P100']:.4f} / "
      f"{m1_def['P100']:.4f} / {m2_def['P100']:.4f}); concordance {fwd.frac_le100:.4f} fwd "
      f"({FOLD:.0f}x null) / {rev.frac_le100:.4f} rev (ceiling {CEIL:.4f}); "
      f"matched-N: donor 1 dominates on both axes")
