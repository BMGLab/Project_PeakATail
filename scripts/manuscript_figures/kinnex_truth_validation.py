#!/usr/bin/env python3
"""
kinnex_truth_validation.py -- TIER-1, ATLAS-INDEPENDENT validation of the pbmc_10k_v3
head-to-head using PacBio Kinnex single-cell Iso-Seq long reads as ground truth.

Truth substrate: `isoseq groupdedup` FLNC molecules (1 record = 1 CB/UMI molecule that
passed `isoseq refine --require-polya`), aligned to GRCh38 (Ensembl no-chr), 3' terminus
taken per molecule, internal-priming-filtered, single-linkage clustered at 25 nt, >=5 UMIs.

CAVEAT carried on every panel: the Kinnex donors are NOT the pbmc_10k_v3 donor and the
chemistries differ, so this is SITE-LEVEL (not cell-matched) truth.

Inputs : results/benchmark_tools/kinnex_truth/<truth_tag>/score_<tool>_vs_<truth_tag>.tsv
         (written by scripts/benchmark_tools/score_tool.py with --atlas <truth point bed>)
         plus the filter-cascade reports from build_truth.sh.
Outputs: PNG 300dpi -> results/figures/manuscript/
         PNG+PDF    -> manuscript/figures/
         every plotted value -> results/figures/manuscript/kinnex_truth_validation.tsv
"""
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as _mpl
_mpl.rcParams["pdf.fonttype"] = 42
_mpl.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
SCORE_ROOT = WD / "results/benchmark_tools/kinnex_truth"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "kinnex_truth_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"

# de novo tools ranked by the ATLAS benchmark (manuscript/09_headtohead_results.md, F1@100
# vs PolyASite 2.0 restricted to detected genes). scUTRquant is annotation-based (catalog).
ATLAS_F1 = {"polyapipe": 0.261, "scapture": 0.199, "sierra": 0.177,
            "scapatrap": 0.150, "peakatail": 0.134, "scutrquant": 0.290}
DENOVO = ["polyapipe", "scapture", "sierra", "scapatrap", "peakatail"]
LABEL = {"polyapipe": "polyApipe", "scapture": "SCAPTURE", "sierra": "Sierra",
         "scapatrap": "scAPAtrap", "peakatail": "PeakATail", "scutrquant": "scUTRquant*"}
COLOR = {t: PALETTE[i] for i, t in enumerate(
    ["peakatail", "polyapipe", "scapture", "sierra", "scapatrap", "scutrquant"])}
CUT = 100


def save_manuscript(fig, name, **kw):
    for ext in ("png", "pdf"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), **kw)
        print("wrote", p)
    p = OUTDIR / f"{name}.png"
    fig.savefig(p, dpi=300, **kw)
    print("wrote", p)


def load(truth_tag):
    """Return tidy frame of every score row for one truth set."""
    d = SCORE_ROOT / truth_tag
    rows = []
    for f in sorted(d.glob(f"score_*_vs_{truth_tag}.tsv")):
        tool = re.sub(rf"^score_(.*)_vs_{truth_tag}\.tsv$", r"\1", f.name)
        df = pd.read_csv(f, sep="\t")
        df["tool"] = tool
        df["truth"] = truth_tag
        rows.append(df)
    if not rows:
        sys.exit(f"no score TSVs under {d}")
    return pd.concat(rows, ignore_index=True)


def cascade(truth_tag):
    """Parse the filter-cascade report emitted by build_truth.sh."""
    rep = SCORE_ROOT / truth_tag / f"{truth_tag}_filter_report.txt"
    txt = rep.read_text() if rep.exists() else ""
    g = lambda p: (int(re.search(p, txt).group(1)) if re.search(p, txt) else np.nan)
    return {
        "molecules_kept_cb_clip": g(r"termini_total.*?=\s*(\d+)"),
        "unique_positions": g(r"unique_terminus_positions\s*=\s*(\d+)"),
        "ip_positions": g(r"IP-flagged unique positions\s*=\s*(\d+)"),
        "termini_in": g(r"\[cluster:truth\] termini_in=(\d+)"),
        "termini_ip": g(r"\[cluster:truth\].*internal_priming_flagged=(\d+)"),
        "termini_used": g(r"\[cluster:truth\].*used=(\d+)"),
        "truth_pas": g(r"\[cluster:truth\].*clusters_kept\(>=5 UMI\)=(\d+)"),
        "clusters_dropped": g(r"\[cluster:truth\].*clusters_dropped_lt5UMI=(\d+)"),
        "decoy_pas": g(r"\[cluster:decoy\].*clusters_kept\(>=5 UMI\)=(\d+)"),
        "hex_canonical": (float(re.search(r"AATAAA/ATTAAA = ([0-9.]+)", txt).group(1))
                          if "AATAAA/ATTAAA" in txt else np.nan),
        "hex_any": (float(re.search(r"15 PAS hexamers = ([0-9.]+)", txt).group(1))
                    if "15 PAS hexamers" in txt else np.nan),
    }


def main():
    truths = [t for t in sys.argv[1:]] or ["gemx"]
    frames = {t: load(t) for t in truths}
    casc = {t: cascade(t) for t in truths}
    prim = truths[0]                      # panel B/C/D use the first (primary) truth
    df = frames[prim]

    def val(tool, panel, series, ref, cut):
        s = df[(df.tool == f"{tool}_vs_{prim}") & (df.panel == panel)
               & (df.series == series) & (df.reference == ref) & (df.cutoff_bp == cut)]
        return float(s.value.mean()) if len(s) else np.nan

    tools = [t for t in ["polyapipe", "scapture", "sierra", "scapatrap",
                         "peakatail", "scutrquant"]
             if f"{t}_vs_{prim}" in set(df.tool)]

    tidy = []
    for t in tools:
        for c in [10, 25, 50, 100, 200]:
            p = val(t, "precision", "real", "atlas_full", c)
            r = val(t, "recall", "real", "atlas_full", c)
            n = val(t, "precision", "null_genic", "atlas_full", c)
            f1 = 0.0 if (p + r) == 0 else 2 * p * r / (p + r)
            tidy.append(dict(truth=prim, tool=t, cutoff_bp=c, precision_vs_truth=p,
                             recall_of_truth=r, f1_vs_truth=f1, null_precision=n,
                             precision_over_null=(p / n if n else np.nan),
                             atlas_f1_at100=ATLAS_F1.get(t, np.nan)))
    tidy = pd.DataFrame(tidy)

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(13.6, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26,
                          left=0.075, right=0.985, top=0.895, bottom=0.115)

    def style(ax):
        ax.set_facecolor("white")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)

    # A: truth-set construction cascade
    axA = fig.add_subplot(gs[0, 0]); style(axA)
    c = casc[prim]
    bars = [c["termini_in"], c["termini_used"]]
    names = ["molecule 3' termini\n(real-cell CB, 3' clip <=30 nt)",
             "after internal-priming\nfilter (>=12 A/18 nt or >=6 A run)"]
    axA.bar(range(len(bars)), [b / 1e6 for b in bars], color=[PALETTE[0], PALETTE[2]],
            width=0.6)
    for i, b in enumerate(bars):
        axA.text(i, b / 1e6, f"{b/1e6:,.1f} M", ha="center", va="bottom",
                 fontsize=9, color=INK)
    axA.set_xticks(range(len(bars))); axA.set_xticklabels(names, fontsize=8.5, color=INK)
    axA.set_ylabel("million molecules (= UMIs)", fontsize=9.5, color=INK)
    axA.set_title(f"A  Long-read truth construction  ({int(c['truth_pas']):,} truth PAS, "
                  f"{int(c['decoy_pas']):,} internal-priming decoys)",
                  fontsize=10.5, color=INK, loc="left")
    axA.text(0.02, 0.86, f"canonical AATAAA/ATTAAA in -40..-10: "
                         f"{c['hex_canonical']*100:.1f}%\nany of 15 PAS hexamers: "
                         f"{c['hex_any']*100:.1f}%",
             transform=axA.transAxes, fontsize=8.5, color=MUTED, va="top")

    # B: precision vs truth
    axB = fig.add_subplot(gs[0, 1]); style(axB)
    order = tidy[tidy.cutoff_bp == CUT].sort_values("precision_vs_truth", ascending=False)
    xs = np.arange(len(order))
    axB.bar(xs, order.precision_vs_truth, color=[COLOR[t] for t in order.tool], width=0.62)
    axB.plot(xs, order.null_precision, ls="none", marker="_", ms=18, mew=2.2,
             color=INK, label="gene-body null (3 seeds)")
    for x, (p, n) in enumerate(zip(order.precision_vs_truth, order.null_precision)):
        axB.text(x, p, f"{p:.3f}", ha="center", va="bottom", fontsize=8.5, color=INK)
    axB.set_xticks(xs); axB.set_xticklabels([LABEL[t] for t in order.tool],
                                            fontsize=8.5, color=INK, rotation=20, ha="right")
    axB.set_ylabel(f"precision @ {CUT} bp vs long-read truth", fontsize=9.5, color=INK)
    axB.set_title("B  Precision against orthogonal long-read PAS", fontsize=10.5,
                  color=INK, loc="left")
    axB.legend(frameon=False, fontsize=8.5, labelcolor=MUTED, loc="upper right")

    # C: recall of truth
    axC = fig.add_subplot(gs[1, 0]); style(axC)
    orderC = tidy[tidy.cutoff_bp == CUT].sort_values("recall_of_truth", ascending=False)
    xs = np.arange(len(orderC))
    axC.bar(xs, orderC.recall_of_truth, color=[COLOR[t] for t in orderC.tool], width=0.62)
    for x, r in enumerate(orderC.recall_of_truth):
        axC.text(x, r, f"{r:.3f}", ha="center", va="bottom", fontsize=8.5, color=INK)
    axC.set_xticks(xs); axC.set_xticklabels([LABEL[t] for t in orderC.tool],
                                            fontsize=8.5, color=INK, rotation=20, ha="right")
    axC.set_ylabel(f"recall of truth PAS @ {CUT} bp", fontsize=9.5, color=INK)
    axC.set_title("C  Recall of long-read truth PAS", fontsize=10.5, color=INK, loc="left")

    # D: does the atlas ranking hold?
    axD = fig.add_subplot(gs[1, 1]); style(axD)
    dn = tidy[(tidy.cutoff_bp == CUT) & (tidy.tool.isin(DENOVO))]
    ra = dn.sort_values("atlas_f1_at100", ascending=False).tool.tolist()
    rk = dn.sort_values("f1_vs_truth", ascending=False).tool.tolist()
    for t in dn.tool:
        y0, y1 = len(ra) - ra.index(t), len(rk) - rk.index(t)
        axD.plot([0, 1], [y0, y1], color=COLOR[t], lw=2.4, marker="o", ms=8)
        axD.text(-0.06, y0, LABEL[t], ha="right", va="center", fontsize=9, color=COLOR[t])
        axD.text(1.06, y1, LABEL[t], ha="left", va="center", fontsize=9, color=COLOR[t])
    axD.set_xlim(-0.55, 1.55); axD.set_xticks([0, 1])
    axD.set_xticklabels(["PolyASite 2.0 atlas\n(manuscript/09)",
                         "Kinnex long-read truth\n(this figure)"], fontsize=9, color=INK)
    axD.set_yticks(range(1, len(ra) + 1))
    axD.set_yticklabels([f"#{len(ra) - i + 1}" for i in range(1, len(ra) + 1)],
                        fontsize=9, color=MUTED)
    axD.yaxis.grid(False)
    axD.set_ylabel("F1 rank among de novo tools", fontsize=9.5, color=INK)
    axD.set_title("D  Does the atlas ranking survive orthogonal truth?", fontsize=10.5,
                  color=INK, loc="left")

    fig.suptitle("Atlas-independent validation of the pbmc_10k_v3 head-to-head against "
                 "PacBio Kinnex single-cell long reads", fontsize=13, color=INK, x=0.075,
                 ha="left", y=0.965)
    fig.text(0.075, 0.028,
             "CAVEAT  The Kinnex donor is NOT the pbmc_10k_v3 donor and the 10x chemistries differ, so this is SITE-LEVEL truth, not cell-matched: it validates PAS POSITIONS, not per-cell usage.\n"
             "Truth = dedup FLNC molecules (isoseq refine --require-polya), minimap2 -ax splice:hq -uf, real-cell CB, 3' softclip <=30 nt, internal-priming filtered, 25 nt single-linkage, >=5 UMIs.\n"
             "Truth is limited to PBMC-expressed genes at >=5 UMIs, so recall is bounded by long-read depth, not only by tool sensitivity. scUTRquant* is annotation-based (fixed catalog) and is shown but not ranked with de novo tools.",
             fontsize=7.6, color=MUTED, ha="left", va="bottom")

    save_manuscript(fig, NAME, facecolor="white")
    plt.close(fig)

    allc = pd.DataFrame([dict(truth=t, **casc[t]) for t in truths])
    out = OUTDIR / f"{NAME}.tsv"
    with open(out, "w") as fh:
        fh.write("# panel B/C/D values (one row per tool x cutoff), truth = "
                 + ",".join(truths) + "\n")
        tidy.to_csv(fh, sep="\t", index=False, float_format="%.6f")
        fh.write("\n# panel A: truth-set construction cascade\n")
        allc.to_csv(fh, sep="\t", index=False, float_format="%.6f")
    print("wrote", out)
    for t in truths[1:]:
        print(f"[replicate {t}] cascade:", casc[t])


if __name__ == "__main__":
    main()
