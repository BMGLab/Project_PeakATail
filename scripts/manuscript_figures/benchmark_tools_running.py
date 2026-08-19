#!/usr/bin/env python3
"""
benchmark_tools_running.py -- the RUNNING cross-tool comparison table+figure.
External tools so far: polyApipe 0.1.0 and Sierra 0.99.27 on pbmc_10k_v3
(each with two reduction variants), scored by
scripts/benchmark_tools/score_tool.py with byte-identical machinery/references
to benchmark_curated.py.  Sierra variants: 3' of the fitted peak interval
(Fit.start..Fit.end -> strand-aware 3' base; primary, matches how all tools
are reduced) and MaxPosition raw-coverage summit (sensitivity).  PeakATail rows (lg_annotate, B1_cohort_full) are
pulled from benchmark_curated.tsv and are PROVISIONAL: they come from the
Laughney cohort -- the matched CellRanger pbmc_10k_v3 PeakATail run is still
in progress.  Re-run this script as each new tool lands.

Every plotted number is in benchmark_tools_running.tsv.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "benchmark_tools_running"
for d in (OUTDIR, FIGDIR):
    d.mkdir(parents=True, exist_ok=True)


def save_manuscript(fig, name, **kw):
    fig.savefig(OUTDIR / f"{name}.png", dpi=300, **kw)
    print("wrote", OUTDIR / f"{name}.png")
    for ext in ("png", "pdf"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), **kw)
        print("wrote", p)


INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

# entities (fixed order + fixed hue; shade-pairs group tool families)
TOOLS = ["polyapipe", "polyapipe_all", "sierra", "sierra_summit",
         "lg_annotate", "B1_cohort_full"]
LAB = {"polyapipe": "polyApipe\n(misprime excl.)",
       "polyapipe_all": "polyApipe\n(all peaks)",
       "sierra": "Sierra\n(3' of fit)",
       "sierra_summit": "Sierra\n(summit)",
       "lg_annotate": "PeakATail\nlg_annotate*",
       "B1_cohort_full": "PeakATail\nB1_cohort_full*"}
COLOR = {"polyapipe": "#0072B2", "polyapipe_all": "#56B4E9",
         "sierra": "#009E73", "sierra_summit": "#66C6A3",
         "lg_annotate": "#D55E00", "B1_cohort_full": "#E69F00"}
DATASET = {"polyapipe": "pbmc_10k_v3", "polyapipe_all": "pbmc_10k_v3",
           "sierra": "pbmc_10k_v3", "sierra_summit": "pbmc_10k_v3",
           "lg_annotate": "laughney_cohort", "B1_cohort_full": "laughney_cohort"}
PROVISIONAL = {"lg_annotate", "B1_cohort_full"}

# ---------------------------------------------------------------- load scores
PA_DIR = WD / "results/benchmark_tools/pbmc_10k_v3/polyapipe"
SI_DIR = WD / "results/benchmark_tools/pbmc_10k_v3/sierra"
parts = []
for t in ("polyapipe", "polyapipe_all"):
    d = pd.read_csv(PA_DIR / f"score_{t}.tsv", sep="\t")
    parts.append(d)
for t in ("sierra", "sierra_summit"):
    d = pd.read_csv(SI_DIR / f"score_{t}.tsv", sep="\t")
    parts.append(d)
bc = pd.read_csv(OUTDIR / "benchmark_curated.tsv", sep="\t")
parts.append(bc[bc["arm"].isin(["lg_annotate", "B1_cohort_full"])].copy())
df = pd.concat(parts, ignore_index=True)
df["tool"] = df["arm"]
df["dataset"] = df["arm"].map(DATASET)
df["provisional"] = df["arm"].isin(PROVISIONAL).astype(int)
df["note"] = np.where(df["provisional"] == 1,
                      "Laughney cohort; pbmc_10k_v3 PeakATail rerun pending", "")

# ------------------------------------------------------------- resource rows
# polyApipe: /usr/bin/time -v in results/benchmark_tools/pbmc_10k_v3/polyapipe/run.log
#   Elapsed 3:26:37 = 12397 s; Max RSS 13,723,976 kbytes (KiB) = 13.088 GiB.
#   One run covers both variants (peaks GFF + counting).
# Sierra: /usr/bin/time -v in results/benchmark_tools/pbmc_10k_v3/sierra/run.log
#   Elapsed 2:25:04 = 8704 s; Max RSS 10,050,780 kbytes (KiB) = 9.585 GiB.
#   One run (FindPeaks 16 cores + CountPeaks) covers both reduction variants.
# PeakATail: resources.jsonl (5 s sampling of the pipeline process, GiB), full
#   pipeline incl. clustering/differential -- much wider scope than polyApipe.
RES = {  # tool -> (runtime_s, rss_gib, how)
    "polyapipe": (12397.0, 13723976 * 1024 / 1024**3,
                  "/usr/bin/time -v (whole run; shared by both variants)"),
    "sierra": (8704.0, 10050780 * 1024 / 1024**3,
               "/usr/bin/time -v (whole run; shared by both variants)"),
    "lg_annotate": (26323.989, 3.6043,
                    "resources.jsonl 5s samples (full pipeline incl. downstream)"),
    "B1_cohort_full": (62405.595, 9.2226,
                       "resources.jsonl 5s samples (full pipeline incl. downstream)"),
}
res_rows = []
for t, (rt, rss, how) in RES.items():
    for series, val in (("runtime_s", rt), ("runtime_h", rt / 3600.0),
                        ("max_rss_gib", rss)):
        res_rows.append(dict(panel="resource", arm=t, series=series,
                             reference=how, cutoff_bp=np.nan, replicate=0,
                             n_query=np.nan, n_matched=np.nan, value=val,
                             tool=t, dataset=DATASET[t],
                             provisional=int(t in PROVISIONAL),
                             note={"polyapipe": "polyApipe run also covers "
                                                "the _all variant",
                                   "sierra": "Sierra run also covers the "
                                             "_summit variant"}.get(
                                   t, "Laughney cohort; wider pipeline scope; "
                                      "pbmc_10k_v3 rerun pending")))
for t, kb in (("polyapipe", 13723976), ("sierra", 10050780)):
    res_rows.append(dict(panel="resource", arm=t, series="max_rss_kbytes",
                         reference="/usr/bin/time -v raw", cutoff_bp=np.nan,
                         replicate=0, n_query=np.nan, n_matched=np.nan,
                         value=kb, tool=t, dataset="pbmc_10k_v3",
                         provisional=0, note="raw kbytes (KiB) from run.log"))
df = pd.concat([df, pd.DataFrame(res_rows)], ignore_index=True)
df.to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False, float_format="%.6f")
print("wrote", OUTDIR / f"{NAME}.tsv", f"({len(df)} rows)")


def val(tool, panel, ref, cutoff, series="real"):
    m = df[(df.tool == tool) & (df.panel == panel) & (df.series == series)
           & (df.reference == ref) & (df.cutoff_bp == cutoff)]
    return float(m["value"].iloc[0])


def nulls(tool):
    m = df[(df.tool == tool) & (df.panel == "precision")
           & (df.series == "null_genic") & (df.cutoff_bp == 100)]
    return m["value"].to_list()


def ncalled(tool):
    m = df[(df.tool == tool) & (df.panel == "meta") & (df.series == "n_points")]
    return int(m["n_matched"].iloc[0])   # points actually scored


# ---------------------------------------------------------------------- figure
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 7.8, "axes.titlesize": 8.6,
    "xtick.labelsize": 6.8, "ytick.labelsize": 7.2, "legend.fontsize": 6.8,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.linewidth": 0.8, "font.family": "DejaVu Sans",
})

fig = plt.figure(figsize=(11.4, 7.4))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.78], hspace=0.62,
                      wspace=1.05, left=0.065, right=0.985, top=0.90,
                      bottom=0.175)
axP = fig.add_subplot(gs[0, 0:3])
axR = fig.add_subplot(gs[0, 3:6])
axN = fig.add_subplot(gs[1, 0:2])
axT = fig.add_subplot(gs[1, 2:4])
axM = fig.add_subplot(gs[1, 4:6])


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def bars(ax, tools, values, fmt, unit_lab=None):
    xs = np.arange(len(tools))
    for x, t, v in zip(xs, tools, values):
        ax.bar(x, v, 0.62, color=COLOR[t], zorder=3,
               hatch="//" if t in PROVISIONAL else None,
               edgecolor=SURFACE, linewidth=0.8)
        ax.annotate(fmt(v), (x, v), xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=6.8, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([LAB[t] for t in tools])
    if unit_lab:
        ax.set_ylabel(unit_lab)
    style(ax)


# A precision@100 vs full atlas + genic-null seeds + roadmap bar
pv = [val(t, "precision", "atlas_full", 100) for t in TOOLS]
bars(axP, TOOLS, pv, lambda v: f"{v:.3f}", "precision @100 bp")
for x, t in enumerate(TOOLS):
    ns = nulls(t)
    axP.scatter([x] * len(ns), ns, marker="D", s=9, color=INK, zorder=4)
axP.axhline(0.70, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=2)
axP.text(len(TOOLS) - 0.42, 0.705, "roadmap 0.70", fontsize=6.2, color=MUTED,
         ha="right", va="bottom")
axP.set_ylim(0, 0.80)
axP.set_title("Precision vs PolyASite 2.0 rep sites (point mode, strand-matched)")

# B recall@100 vs detected-gene-restricted atlas + roadmap bar
rv = [val(t, "recall", "atlas_detected", 100) for t in TOOLS]
bars(axR, TOOLS, rv, lambda v: f"{v:.3f}", "recall @100 bp")
axR.axhline(0.60, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=2)
axR.text(len(TOOLS) - 0.42, 0.605, "roadmap 0.60", fontsize=6.2, color=MUTED,
         ha="right", va="bottom")
axR.set_ylim(0, 0.68)
axR.set_title("Recall vs atlas in 14,851 detected-gene bodies (shared denominator)")

# C n_called
nv = [ncalled(t) for t in TOOLS]
bars(axN, TOOLS, nv, lambda v: f"{v/1000:.0f}k", "PAS called (n scored)")
axN.set_xticklabels([LAB[t].replace("\n", " ") for t in TOOLS],
                    rotation=30, ha="right", fontsize=6.0)
axN.yaxis.set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:.0f}k"))
axN.set_title("Sites called")

# D runtime  /  E peak RSS  (one bar per RUN; polyApipe run covers both variants)
RTOOLS = ["polyapipe", "sierra", "lg_annotate", "B1_cohort_full"]
RLAB = dict(LAB, polyapipe="polyApipe\n(one run,\nboth variants)",
            sierra="Sierra\n(one run,\nboth variants)")
tv = [RES[t][0] / 3600.0 for t in RTOOLS]
bars(axT, RTOOLS, tv, lambda v: f"{v:.1f}h", "wall time (h)")
axT.set_xticklabels([RLAB[t] for t in RTOOLS])
axT.set_title("Runtime (not comparable: scope differs)")
mv = [RES[t][1] for t in RTOOLS]
bars(axM, RTOOLS, mv, lambda v: f"{v:.1f}", "peak RSS (GiB)")
axM.set_xticklabels([RLAB[t] for t in RTOOLS])
axM.set_title("Peak memory")

handles = [Patch(facecolor=COLOR[t], edgecolor=SURFACE,
                 hatch="//" if t in PROVISIONAL else None,
                 label=LAB[t].replace("\n", " ").replace("*", "")) for t in TOOLS]
handles += [Line2D([], [], marker="D", ls="", color=INK, markersize=3,
                   label="genic-shuffle null (3 seeds)"),
            Patch(facecolor="#FFFFFF", edgecolor=MUTED, hatch="//",
                  label="* provisional (Laughney cohort)")]
fig.legend(handles=handles, ncol=3, loc="lower center",
           bbox_to_anchor=(0.5, 0.055), frameon=False)

fig.suptitle("Tool benchmark -- RUNNING comparison (2 external tools in)",
             x=0.065, ha="left", fontsize=10, fontweight="bold", color=INK)
fig.text(0.065, 0.008,
         "CAVEATS: PeakATail bars (*) are PROVISIONAL -- scored on the Laughney cohort, not pbmc_10k_v3; the matched CellRanger pbmc_10k_v3 PeakATail run is still in progress. "
         "polyApipe and Sierra scored on pbmc_10k_v3 (CellRanger BAM); Sierra reduced to the strand-aware 3' base of Fit.start..Fit.end (primary) or MaxPosition summit (sensitivity). Cross-dataset precision/recall differences partly reflect dataset depth/composition, not tool quality alone.\n"
         "Runtime/RSS scope differs: polyApipe/Sierra = PAS calling + per-cell counting (/usr/bin/time, one run for both variants); PeakATail = full pipeline incl. clustering/differential (5 s process sampling, may miss subprocess peaks). "
         "Recall denominator identical for all tools: 285,220 PolyASite 2.0 rep sites in the 14,851 Laughney-detected genes (recall flavor b of benchmark_curated.py).",
         fontsize=5.9, color=MUTED, va="bottom", wrap=True)

save_manuscript(fig, NAME)
plt.close(fig)
