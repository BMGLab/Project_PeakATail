#!/usr/bin/env python3
"""
PeakATail manuscript figure: PARAMETER ROBUSTNESS across the 13 reannotate branches.

Source (READ-ONLY):
  /mnt/ssd2/.../RERUN_2026-08_fixed/runs/reannotate/*/branch_manifest.json   (params + final_cells/final_pas)
  /mnt/ssd2/.../reannotate/*/07_clustering/<ds>/stage_stats.json             (input_pas denominator)
  /mnt/ssd2/.../reannotate/*/07_clustering/<ds>/clusters.h5ad                (obs/leiden codes ONLY, via h5py)

Outputs:
  results/figures/manuscript/parameter_sweep.png            (300 dpi)
  results/figures/manuscript/parameter_sweep.tsv            (exact plotted values)
  results/figures/manuscript/parameter_sweep_per_dataset.tsv (full per-branch x per-dataset audit table)

Panels
  a  final_pas vs max_gene_distance (annotation trim sweep), median + IQR over 17 datasets
  b  final_cells vs max_gene_distance, same
  c  n leiden clusters across the clustering arm (A3_* branches + default)
  d  paired per-dataset % change across each knob's full span -> which knob dominates
"""
import json
import os
import glob
import itertools

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = "/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs/reannotate"
OUTDIR = "/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript"
NAME = "parameter_sweep"

# ── validated palette (do not invent colors) ──────────────────────────────────
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"
# Color follows the ENTITY = the measured outcome, fixed across every panel:
CPAS, CCELL, CCLUST = C[0], C[1], C[2]
# Slots 1-3 all-pairs CVD dE (Machado sev-1.0, OKLab x100): 21.9 / 18.0 / 11.0 -> all PASS;
# contrast vs white 5.19 / 3.87 / 3.42 -> all >= 3.0. No floor-level pair is used.

BASELINE = "A2_trim_default"   # d5000, mult2.0, ext=off, leiden_tfidf, res1.0, nn30


# ── 1. extract ────────────────────────────────────────────────────────────────
def load() -> pd.DataFrame:
    rows = []
    for mf in sorted(glob.glob(os.path.join(ROOT, "*/branch_manifest.json"))):
        bdir = os.path.dirname(mf)
        branch = os.path.basename(bdir)
        m = json.load(open(mf))
        trim, clus, filt = m["trim"], m["clustering"], m["filters"]
        for ds in m["datasets"]:
            dsid = ds["dataset_id"]
            rec = dict(
                branch=branch, dataset_id=dsid,
                max_gene_distance=trim["max_gene_distance"],
                utr_multiplier=trim["utr_multiplier"],
                include_extended=bool(trim["include_extended"]),
                method=clus["method"], resolution=clus["resolution"],
                n_neighbors=clus["n_neighbors"],
                min_read=filt["min_read"], min_cells=filt["min_cells"],
                min_pas_per_cell=filt["min_pas_per_cell"],
                final_cells=ds["final_cells"], final_pas=ds["final_pas"],
                input_pas=np.nan, annotated_pas=np.nan, n_clusters=np.nan,
            )
            ss = os.path.join(bdir, "07_clustering", dsid, "stage_stats.json")
            if os.path.exists(ss):
                s = json.load(open(ss))
                rec["input_pas"] = s.get("input_matrix", {}).get("input_pas", np.nan)
                rec["annotated_pas"] = s.get("annotated", {}).get("annotated_pas", np.nan)
            h5 = os.path.join(bdir, "07_clustering", dsid, "clusters.h5ad")
            if os.path.exists(h5):
                # read ONLY the leiden code vector -- never touches X / the matrix
                with h5py.File(h5, "r") as f:
                    codes = f["obs/leiden/codes"][:]
                rec["n_clusters"] = int(len(np.unique(codes[codes >= 0])))
            rows.append(rec)
    df = pd.DataFrame(rows)
    df["pas_retention"] = df["final_pas"] / df["input_pas"]
    return df


def mq(df, branch, col):
    """median, q1, q3 of `col` across the 17 datasets of `branch`."""
    s = df.loc[df.branch == branch, col].astype(float)
    return s.median(), s.quantile(0.25), s.quantile(0.75)


def paired_pct(df, lo, hi, col):
    """per-dataset paired % change from branch `lo` to branch `hi`.

    Returns median, min, max, q1, q3.  min/max are quoted in the captions;
    q1/q3 are the whiskers actually drawn in panel d.
    """
    a = df[df.branch == lo].set_index("dataset_id")[col].astype(float)
    b = df[df.branch == hi].set_index("dataset_id")[col].astype(float)
    pc = 100.0 * (b - a) / a
    return pc.median(), pc.min(), pc.max(), pc.quantile(0.25), pc.quantile(0.75)


def paired_n_changed(df, lo, hi, col):
    """how many of the 17 datasets change AT ALL between the two branches."""
    a = df[df.branch == lo].set_index("dataset_id")[col].astype(float)
    b = df[df.branch == hi].set_index("dataset_id")[col].astype(float)
    return int((b - a).ne(0).sum())


df = load()
N_DS = df.dataset_id.nunique()
assert df.branch.nunique() == 13 and N_DS == 17, (df.branch.nunique(), N_DS)

# trim sweep: max_gene_distance at include_extended=True
TRIM = [("A2_trim_d1000_ext", 1000), ("A2_trim_d2000_ext", 2000),
        ("A2_trim_d3000_ext", 3000), ("A2_trim_d5000_ext", 5000),
        ("A2_trim_d10000_ext", 10000)]
# clustering arm: label -> branch  (default sits inside the resolution ladder)
CLUS = [("res 0.5", "A3_res0.5"), ("res 1.0\n(default)", BASELINE), ("res 2.0", "A3_res2.0"),
        ("nn 15", "A3_nn15"), ("nn 50", "A3_nn50"), ("libsize\nnorm.", "A3_libsize")]
# knob spans for the dominance panel: (short label, low branch, high branch, arm)
SPANS = [("max_gene\ndistance\n1k→2k…10k", "A2_trim_d1000_ext", "A2_trim_d10000_ext", "annotation"),
         ("include\nextended\noff→on", BASELINE, "A2_trim_d5000_ext", "annotation"),
         ("utr\nmultiplier\n1.5→3.0", "A2_trim_mult1.5", "A2_trim_mult3.0", "annotation"),
         ("leiden\nresolution\n0.5→2.0", "A3_res0.5", "A3_res2.0", "clustering"),
         ("n_neighbors\n\n15→50", "A3_nn15", "A3_nn50", "clustering"),
         ("cluster\nmethod\ntfidf→libsize", BASELINE, "A3_libsize", "clustering")]

# ── 2. assemble every number that will be drawn ───────────────────────────────
plotted = []
pa = {"x": [], "med": [], "q1": [], "q3": []}
pb = {"x": [], "med": [], "q1": [], "q3": []}
for br, dist in TRIM:
    m, q1, q3 = mq(df, br, "final_pas")
    pa["x"].append(dist); pa["med"].append(m); pa["q1"].append(q1); pa["q3"].append(q3)
    plotted.append(dict(panel="a", series="include_extended = on", branch=br, x_label=dist,
                        metric="final_pas", median=m, q1=q1, q3=q3, n_datasets=N_DS))
    m, q1, q3 = mq(df, br, "final_cells")
    pb["x"].append(dist); pb["med"].append(m); pb["q1"].append(q1); pb["q3"].append(q3)
    plotted.append(dict(panel="b", series="include_extended = on", branch=br, x_label=dist,
                        metric="final_cells", median=m, q1=q1, q3=q3, n_datasets=N_DS))

ref_pas = mq(df, BASELINE, "final_pas")
ref_cells = mq(df, BASELINE, "final_cells")
plotted.append(dict(panel="a", series="include_extended = off (default)", branch=BASELINE, x_label=5000,
                    metric="final_pas", median=ref_pas[0], q1=ref_pas[1], q3=ref_pas[2], n_datasets=N_DS))
plotted.append(dict(panel="b", series="include_extended = off (default)", branch=BASELINE, x_label=5000,
                    metric="final_cells", median=ref_cells[0], q1=ref_cells[1], q3=ref_cells[2], n_datasets=N_DS))

pc = []
for lab, br in CLUS:
    m, q1, q3 = mq(df, br, "n_clusters")
    pc.append((lab.replace("\n", " "), br, m, q1, q3))
    plotted.append(dict(panel="c", series="n leiden clusters", branch=br,
                        x_label=lab.replace("\n", " "), metric="n_clusters",
                        median=m, q1=q1, q3=q3, n_datasets=N_DS))

pd_rows = []
for lab, lo, hi, arm in SPANS:
    r = {"label": lab.replace("\n", " ").strip(), "arm": arm, "lo": lo, "hi": hi}
    for metric, key in (("final_pas", "pas"), ("n_clusters", "clu"), ("final_cells", "cell")):
        med, lo_, hi_, q1_, q3_ = paired_pct(df, lo, hi, metric)
        r[key] = med; r[key + "_min"] = lo_; r[key + "_max"] = hi_
        r[key + "_q1"] = q1_; r[key + "_q3"] = q3_
        r[key + "_nchg"] = paired_n_changed(df, lo, hi, metric)
        plotted.append(dict(panel="d", series=f"paired % change in {metric}",
                            branch=f"{lo} -> {hi}", x_label=r["label"], metric=metric,
                            median=med, q1=lo_, q3=hi_, p25=q1_, p75=q3_,
                            n_datasets_changed=r[key + "_nchg"], n_datasets=N_DS))
    pd_rows.append(r)

# ── numbers that appear only in panel captions -- recorded so the .tsv really is
#    "every value printed on the figure", not just every value drawn as a mark ──
ret = df.assign(r=100.0 * df.final_pas / df.input_pas)
ret_def = ret.loc[ret.branch == BASELINE, "r"]
plotted.append(dict(panel="a", series="caption number", branch=BASELINE, x_label="pas_retention (%)",
                    metric="pas_retention_pct", median=ret_def.median(),
                    q1=ret_def.min(), q3=ret_def.max(), n_datasets=N_DS))
_m1000 = mq(df, "A2_trim_d1000_ext", "final_pas")
_iqr_w = float(_m1000[2] - _m1000[1])
_dist_eff = float(mq(df, "A2_trim_d10000_ext", "final_pas")[0] - mq(df, "A2_trim_d1000_ext", "final_pas")[0])
plotted.append(dict(panel="a", series="caption number", branch="A2_trim_d1000_ext",
                    x_label="between-dataset IQR / distance effect", metric="ratio",
                    median=_iqr_w / _dist_eff, q1=_iqr_w, q3=_dist_eff, n_datasets=N_DS))
# worst-case cell movement over ALL ordered branch pairs (the "<= 3.8%" in panel b)
_cell_worst, _cell_medmax = 0.0, 0.0
for _a, _b in itertools.permutations(sorted(df.branch.unique()), 2):
    _m, _lo, _hi, _, _ = paired_pct(df, _a, _b, "final_cells")
    _cell_worst = max(_cell_worst, abs(_lo), abs(_hi))
    _cell_medmax = max(_cell_medmax, abs(_m))
plotted.append(dict(panel="b", series="caption number", branch="all ordered branch pairs",
                    x_label="final_cells paired change (%)", metric="final_cells_worst_case",
                    median=_cell_medmax, q1=0.0, q3=_cell_worst, n_datasets=N_DS))

os.makedirs(OUTDIR, exist_ok=True)
out = pd.DataFrame(plotted)
out.columns = [c if c not in ("q1", "q3") else c for c in out.columns]
out = out.rename(columns={"q1": "q1_or_min", "q3": "q3_or_max"})
out.to_csv(os.path.join(OUTDIR, f"{NAME}.tsv"), sep="\t", index=False, float_format="%.4f")
df.to_csv(os.path.join(OUTDIR, f"{NAME}_per_dataset.tsv"), sep="\t", index=False, float_format="%.6f")

# ── 3. draw ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8, "axes.titlesize": 8.8, "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5, "legend.fontsize": 7.2, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": GRID,
})


def style(ax):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.8)
        ax.spines[s].set_color(GRID)
    ax.tick_params(length=3, width=0.8)


TAGS = []


def panel(ax, letter, title):
    ax.set_title(title, loc="left", color=INK, pad=8)
    TAGS.append((ax, letter))


def place_tags(fig):
    """Letters at a fixed offset left of each axes box (immune to y-label width)."""
    fig.canvas.draw()
    for ax, letter in TAGS:
        bb = ax.get_position()
        fig.text(bb.x0 - 0.082, bb.y1 + 0.012, letter, fontsize=11,
                 fontweight="bold", color=INK, va="bottom", ha="left")


fig = plt.figure(figsize=(7.4, 8.0), facecolor=SURFACE)
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.02, 0.92, 1.0],
              hspace=1.05, wspace=0.34, left=0.125, right=0.985, top=0.885, bottom=0.085)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, :])
ax_d = fig.add_subplot(gs[2, :])

X = np.array(pa["x"], float)
XT, XL = [1000, 2000, 3000, 5000, 10000], ["1k", "2k", "3k", "5k", "10k"]


def sweep_panel(ax, d, ref, color, ylab, scale=1.0):
    style(ax)
    q1 = np.array(d["q1"], float) / scale
    q3 = np.array(d["q3"], float) / scale
    med = np.array(d["med"], float) / scale
    ax.fill_between(X, q1, q3, color=color, alpha=0.15, linewidth=0, zorder=2)
    ax.plot(X, med, color=color, lw=1.8, marker="o", ms=6, mfc=color,
            mec=SURFACE, mew=1.2, zorder=5)
    ax.plot([5000], [ref[0] / scale], marker="D", ms=7.5, color=color, mfc=SURFACE,
            mec=color, mew=1.8, ls="none", zorder=6)
    ax.set_xscale("log")
    ax.set_xticks(XT)
    ax.set_xticklabels(XL)
    ax.minorticks_off()
    ax.set_xlabel("max_gene_distance (bp)")
    ax.set_ylabel(ylab)
    lo = min(q1.min(), ref[1] / scale)
    hi = max(q3.max(), ref[2] / scale)
    pad = 0.10 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    return med


# ---- panel a : PAS vs max_gene_distance --------------------------------------
ax = ax_a
med_a = sweep_panel(ax, pa, ref_pas, CPAS, "PAS retained (thousands)", 1e3)
ax.annotate("extended on", xy=(10000, med_a[-1]), xytext=(-3, 9),
            textcoords="offset points", ha="right", color=CPAS, fontsize=7)
ax.annotate("extended off\n(default)", xy=(5000, ref_pas[0] / 1e3), xytext=(9, -3),
            textcoords="offset points", ha="left", va="top", color=CPAS, fontsize=7)
ax.annotate("", xy=(5000, ref_pas[0] / 1e3), xytext=(5000, med_a[3]),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9, shrinkA=4, shrinkB=4))
ax.text(4750, (ref_pas[0] / 1e3 + med_a[3]) / 2, "+13%", color=INK, fontsize=7,
        ha="right", va="center")
panel(ax, "a", "A 10× wider gene-distance window\nadds only 6% more PAS")
ax.text(0.0, -0.42, "band = between-dataset IQR, ~6× wider than the entire\n"
                    f"distance effect · n = {N_DS} datasets · at the pipeline default,\n"
                    f"{ret_def.median():.1f}% of called peaks are retained ({ret_def.min():.1f}–{ret_def.max():.1f}%);\n"
                    "the filter thresholds that set this were NOT swept",
        transform=ax.transAxes, fontsize=6.3, color=MUTED, va="top")

# ---- panel b : cells vs max_gene_distance ------------------------------------
ax = ax_b
med_b = sweep_panel(ax, pb, ref_cells, CCELL, "cells retained (n)")
ax.annotate("extended on", xy=(10000, med_b[-1]), xytext=(-3, 9),
            textcoords="offset points", ha="right", color=CCELL, fontsize=7)
ax.annotate("extended off (default)\nlies under the same line",
            xy=(5000, ref_cells[0]), xytext=(0, -10), textcoords="offset points",
            ha="center", va="top", color=CCELL, fontsize=7)
panel(ax, "b", "Cell recovery is untouched by the\nannotation window")
ax.text(0.0, -0.42, "band = between-dataset IQR · median is flat to the pixel:\n"
                    f"worst single dataset moves {_cell_worst:.1f}% over all 13 branches,\n"
                    f"and no branch has a median shift above {_cell_medmax:.2f}%\n"
                    f"· n = {N_DS} datasets",
        transform=ax.transAxes, fontsize=6.3, color=MUTED, va="top")

# shared legend for a+b: shape encodes the condition, colour encodes the metric
fig.legend(handles=[
    Line2D([], [], color=INK, lw=1.8, marker="o", ms=6, mec=SURFACE, mew=1.2,
           label="include_extended = on (distance sweep)"),
    Line2D([], [], color=INK, lw=0, marker="D", ms=7.5, mfc=SURFACE, mec=INK,
           mew=1.8, label="include_extended = off (pipeline default)")],
    loc="upper center", bbox_to_anchor=(0.55, 0.995), ncol=2, frameon=False,
    handlelength=1.6, columnspacing=1.8)

# ---- panel c : clusters across the clustering arm ----------------------------
ax = ax_c
style(ax)
labels = [l for l, _ in CLUS]
meds = np.array([r[2] for r in pc], float)
q1 = np.array([r[3] for r in pc], float)
q3 = np.array([r[4] for r in pc], float)
xs = np.arange(len(labels))
ax.bar(xs, meds, width=0.60, color=CCLUST, edgecolor=SURFACE, linewidth=1.2, zorder=3)
ax.errorbar(xs, meds, yerr=[meds - q1, q3 - meds], fmt="none", ecolor=INK,
            elinewidth=1.0, capsize=3.5, capthick=1.0, zorder=4)
base = meds[1]
ax.axhline(base, color=MUTED, lw=0.9, zorder=2)
ax.text(-0.55, base + 0.6, f"default = {base:.0f} clusters", color=MUTED,
        fontsize=6.8, ha="left", va="bottom")
for xi, m, hi in zip(xs, meds, q3):
    ax.text(xi, hi + 1.0, f"{m:.0f}", ha="center", va="bottom", color=INK, fontsize=8)
ax.set_xticks(xs)
ax.set_xticklabels(labels)
ax.set_xlim(-0.62, len(labels) - 0.38)
ax.set_ylim(0, 29)
ax.set_ylabel("leiden clusters (n)")
ax.set_xlabel("clustering branch — all six share one identical preprocessed matrix")
ax.text(5.44, 28.7, "libsize: median near default,\nbut per-dataset range −42% to +93%",
        ha="right", va="top", fontsize=6.4, color=MUTED, linespacing=1.4)
panel(ax, "c", "Resolution alone doubles the cluster count; neighbours and normalisation shift it far less")
ax.text(0.0, -0.44, f"bars = median, whiskers = IQR across {N_DS} datasets · counted as unique obs['leiden'] "
                    "labels in 07_clustering/<dataset>/clusters.h5ad\n"
                    "cluster COUNT is a granularity proxy only — equal counts do not imply equal cell "
                    "assignments; partition agreement is not tested here (see clustering_concordance)",
        transform=ax.transAxes, fontsize=6.3, color=MUTED, va="top")

# ---- panel d : which knob dominates ------------------------------------------
ax = ax_d
style(ax)
short = [l for l, _, _, _ in SPANS]
vals_pas = np.array([r["pas"] for r in pd_rows], float)
vals_clu = np.array([r["clu"] for r in pd_rows], float)
q1_pas = np.array([r["pas_q1"] for r in pd_rows], float)
q3_pas = np.array([r["pas_q3"] for r in pd_rows], float)
q1_clu = np.array([r["clu_q1"] for r in pd_rows], float)
q3_clu = np.array([r["clu_q3"] for r in pd_rows], float)
n_pas = [r["pas_nchg"] for r in pd_rows]
n_clu = [r["clu_nchg"] for r in pd_rows]
xs = np.arange(len(pd_rows))
w = 0.33
xp, xc = xs - w / 2 - 0.015, xs + w / 2 + 0.015
ax.bar(xp, vals_pas, width=w, color=CPAS, edgecolor=SURFACE,
       linewidth=1.2, zorder=3, label="PAS retained")
ax.bar(xc, vals_clu, width=w, color=CCLUST, edgecolor=SURFACE,
       linewidth=1.2, zorder=3, label="leiden clusters")
# IQR whiskers: a bar whose median is 0 but whose whiskers are NOT is a knob that
# moves individual datasets while cancelling out on the median. Without these the
# median-zero bars are indistinguishable from the exactly-zero ones.
for xv, v, lo_, hi_ in list(zip(xp, vals_pas, q1_pas, q3_pas)) + \
                       list(zip(xc, vals_clu, q1_clu, q3_clu)):
    ax.errorbar([xv], [v], yerr=[[v - lo_], [hi_ - v]], fmt="none", ecolor=INK,
                elinewidth=1.0, capsize=3.0, capthick=1.0, zorder=5)
ax.axhline(0, color=MUTED, lw=0.9, zorder=4)
# value label sits above the whisker, not the bar, so the two never collide
for xv, v, hi_, lo_, col, nch in \
        [(x, v, h, l, CPAS, n) for x, v, h, l, n in zip(xp, vals_pas, q3_pas, q1_pas, n_pas)] + \
        [(x, v, h, l, CCLUST, n) for x, v, h, l, n in zip(xc, vals_clu, q3_clu, q1_clu, n_clu)]:
    up = v >= 0
    anchor = max(v, hi_) if up else min(v, lo_)
    ax.annotate(f"{v:+.0f}%", xy=(xv, anchor), xytext=(0, 3 if up else -3),
                textcoords="offset points", ha="center",
                va="bottom" if up else "top",
                color=col if abs(v) > 0.05 else MUTED, fontsize=6.9)
    if abs(v) <= 0.05:
        # a flat bar is ambiguous: either provably identical on every dataset, or
        # a median that cancels real per-dataset movement. Say which.
        ax.annotate(f"{nch}/{N_DS}\nchanged",
                    xy=(xv, anchor), xytext=(0, 13 if up else -13),
                    textcoords="offset points", ha="center",
                    va="bottom" if up else "top",
                    color=MUTED if nch == 0 else INK, fontsize=5.6, linespacing=1.15)
ax.set_xticks(xs)
ax.set_xticklabels(short, fontsize=6.9, linespacing=1.4)
ax.set_ylabel("median paired change (%)")
ax.set_xlabel("parameter knob, swept end to end")
ax.set_ylim(-40, 165)
ax.set_xlim(-0.62, len(pd_rows) - 0.38)
ax.axvline(2.5, color=GRID, lw=1.2, zorder=1)
ax.text(1.0, 155, "annotation knobs", color=MUTED, fontsize=7, ha="center")
ax.text(4.0, 155, "clustering knobs", color=MUTED, fontsize=7, ha="center")
ax.legend(loc="upper left", frameon=False, handlelength=1.3, borderaxespad=0.3)
panel(ax, "d", "Leiden resolution is the one swept knob that dominates: +112% clusters vs ≤13% for any annotation knob")
_ann_lo = min(r["clu_min"] for r in pd_rows[:3])
_ann_hi = max(r["clu_max"] for r in pd_rows[:3])
ax.text(0.0, -0.50,
        f"bars = median per-dataset paired change between the two ends of that knob, whiskers = IQR over {N_DS} datasets; "
        "the n/17 tag counts datasets that changed at all.\n"
        "The two arms are NOT symmetric: clustering knobs leave PAS and cell counts bit-identical (0/17 datasets move, "
        "exactly 0.00%), whereas annotation knobs leave the\nMEDIAN cluster count unchanged but do shift 7–9 of 17 individual "
        f"datasets by 1–2 clusters ({_ann_lo:+.0f}% to {_ann_hi:+.0f}%) — a median-zero bar is not a null result.\n"
        "Also: tfidf→libsize has a small median (+6%) but the widest per-dataset spread of any knob (−42% to +93%), so it is "
        "stable only on average; and the filter\nthresholds (min_read/min_cells/min_pas_per_cell), plausibly the larger lever "
        "on PAS counts, are held fixed in all 13 branches and are not tested here.",
        transform=ax.transAxes, fontsize=6.3, color=MUTED, va="top")

place_tags(fig)
fig.savefig(os.path.join(OUTDIR, f"{NAME}.png"), dpi=300, facecolor=SURFACE,
            bbox_inches="tight", pad_inches=0.16)
print("wrote", os.path.join(OUTDIR, f"{NAME}.png"))
for r in pd_rows:
    print(f"{r['label']:34s} pas {r['pas']:+7.2f}%  clusters {r['clu']:+7.2f}% "
          f"[{r['clu_min']:+.1f},{r['clu_max']:+.1f}]  cells {r['cell']:+.2f}%")
