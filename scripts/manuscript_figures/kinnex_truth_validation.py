#!/usr/bin/env python3
"""
kinnex_truth_validation.py [tag ...]   (default: gemx)

TIER-1, ATLAS-INDEPENDENT validation of the pbmc_10k_v3 head-to-head (manuscript/09) using
PacBio Kinnex single-cell Iso-Seq long reads as ground truth.

Truth substrate: `isoseq groupdedup` FLNC molecules (1 record = 1 CB/UMI molecule that passed
`isoseq refine --require-polya`), genome-aligned, per-molecule 3' terminus, real-cell CB,
3' softclip <= 30 nt, internal-priming split, greedy +/-25 nt peak calling, molecular-support
sweep (>=5 / 20 / 100 / 500 UMIs).

CAVEAT carried on the figure: the Kinnex donor is NOT the pbmc_10k_v3 donor and the 10x
chemistries differ -- SITE-LEVEL truth, not cell-matched.

Inputs  results/benchmark_tools/kinnex_truth/<tag>t{5,20,100,500}/score_<tool>_vs_<tag>t<T>.tsv
        results/benchmark_tools/kinnex_truth/<tag>/<tag>_{call_support_decomposition,
                                              truth_calibration,truth_threshold_qc}.tsv
Outputs PNG 300dpi -> results/figures/manuscript/ ; PNG+PDF(fonttype 42) -> manuscript/figures/
        every plotted value -> results/figures/manuscript/kinnex_truth_validation.tsv
"""
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
from scipy.stats import spearmanr

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
ROOT = WD / "results/benchmark_tools/kinnex_truth"
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
NAME = "kinnex_truth_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
THRESH = [5, 20, 100, 500]
CUT = 100

# manuscript/09_headtohead_results.md, PolyASite 2.0 / detected-gene denominator, @100 bp
ATLAS = {"scutrquant": (0.793, 0.178, 0.290), "polyapipe": (0.380, 0.199, 0.261),
         "scapture": (0.652, 0.118, 0.199), "sierra": (0.256, 0.135, 0.177),
         "scapatrap": (0.100, 0.300, 0.150), "peakatail": (0.118, 0.156, 0.134)}
DENOVO = ["polyapipe", "scapture", "sierra", "scapatrap", "peakatail"]
TOOLS = ["peakatail", "polyapipe", "scapture", "sierra", "scapatrap", "scutrquant"]
LABEL = {"polyapipe": "polyApipe", "scapture": "SCAPTURE", "sierra": "Sierra",
         "scapatrap": "scAPAtrap", "peakatail": "PeakATail", "scutrquant": "scUTRquant*"}
COLOR = dict(zip(TOOLS, PALETTE))
# PeakATail first in TOOLS so it takes the primary blue; keep it visually dominant.


def save_manuscript(fig, name, **kw):
    for ext in ("png", "pdf"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), **kw)
        print("wrote", p)
    p = OUTDIR / f"{name}.png"
    fig.savefig(p, dpi=300, **kw)
    print("wrote", p)


def f1(p, r):
    return 0.0 if (p + r) == 0 else 2 * p * r / (p + r)


def load_scores(tag):
    rows = []
    for T in THRESH:
        d = ROOT / f"{tag}t{T}"
        for tool in TOOLS:
            f = d / f"score_{tool}_vs_{tag}t{T}.tsv"
            if not f.exists():
                continue
            df = pd.read_csv(f, sep="\t")

            def v(panel, series, ref):
                s = df[(df.panel == panel) & (df.series == series)
                       & (df.reference == ref) & (df.cutoff_bp == CUT)]
                return float(s.value.mean()) if len(s) else np.nan

            p = v("precision", "real", "atlas_full")
            r = v("recall", "real", "atlas_full")
            nl = v("precision", "null_genic", "atlas_full")
            nq = df[(df.panel == "precision") & (df.series == "real")
                    & (df.reference == "atlas_full")].n_query.iloc[0]
            rows.append(dict(truth=tag, umi_threshold=T, tool=tool, n_calls=int(nq),
                             precision_vs_truth=p, recall_of_truth=r, f1_vs_truth=f1(p, r),
                             null_precision=nl, precision_over_null=p / nl if nl else np.nan,
                             atlas_precision=ATLAS[tool][0], atlas_recall=ATLAS[tool][1],
                             atlas_f1=ATLAS[tool][2]))
    return pd.DataFrame(rows)


def style(ax):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    tags = sys.argv[1:] or ["gemx"]
    tag = tags[0]
    sc = pd.concat([load_scores(t) for t in tags], ignore_index=True)
    if sc.empty:
        sys.exit("no score TSVs found")
    dec = pd.read_csv(ROOT / tag / f"{tag}_call_support_decomposition.tsv",
                      sep="\t", comment="#")
    cal = pd.read_csv(ROOT / tag / f"{tag}_truth_calibration.tsv", sep="\t", comment="#")
    thq = pd.read_csv(ROOT / tag / f"{tag}_truth_threshold_qc.tsv", sep="\t", comment="#")

    fig = plt.figure(figsize=(14.2, 11.6))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.24,
                          left=0.072, right=0.986, top=0.895, bottom=0.215)

    # ---- A: sequence-only calibration of the truth set -----------------------
    axA = fig.add_subplot(gs[0, 0]); style(axA)
    order = ["5-9", "10-19", "20-49", "50-99", "100-199", "200-499",
             "500-999", "1000-4999", "5000-inf"]
    x = np.arange(len(order))
    for lab, key, col in (("long-read PAS (not internally primed)",
                           "truth_non_internal_priming", PALETTE[2]),
                          ("internal-priming termini (negative control)",
                           "decoy_internal_priming", PALETTE[1])):
        sub = cal[cal.set == key].set_index("umi_bin").reindex(order)
        axA.plot(x, sub.hexamer_canonical, marker="o", ms=5.5, lw=2.2, color=col, label=lab)
    axA.axhline(0.4111, ls="--", lw=1.3, color=MUTED)
    axA.text(len(order) - 1.1, 0.4111, "PolyASite 2.0 rep sites (0.411)", fontsize=8,
             color=MUTED, va="bottom", ha="right")
    axA.axhline(0.3799, ls=":", lw=1.3, color=MUTED)
    axA.text(len(order) - 1.1, 0.3799, "protein-coding TES (0.380)", fontsize=8,
             color=MUTED, va="top", ha="right")
    axA.axvline(0.5, color=PALETTE[1], lw=1.1, alpha=0.45)
    axA.annotate("the plan's >=5 UMI cut sits at the noise floor at 104 M\n"
                 "molecules -- so the truth is reported as a support sweep",
                 xy=(0.5, 0.14), xycoords="data",
                 xytext=(0.035, 0.70), textcoords="axes fraction",
                 fontsize=8.2, color=MUTED, va="top", ha="left",
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0,
                                 connectionstyle="angle3,angleA=0,angleB=90"))
    axA.set_xticks(x); axA.set_xticklabels(order, rotation=35, ha="right", fontsize=8, color=INK)
    axA.set_xlabel("molecular support of the long-read peak (UMIs)", fontsize=9.5, color=INK)
    axA.set_ylabel("fraction with AATAAA / ATTAAA at -40..-10", fontsize=9.5, color=INK)
    axA.set_title("A  The truth set is real, but its support threshold is depth-dependent",
                  fontsize=10.5, color=INK, loc="left")
    axA.legend(frameon=False, fontsize=8.2, labelcolor=MUTED, loc="upper left")

    # ---- B: precision -- atlas vs long-read truth ----------------------------
    axB = fig.add_subplot(gs[0, 1]); style(axB)
    s5 = sc[(sc.truth == tag) & (sc.umi_threshold == 5)].set_index("tool")
    present = [t for t in TOOLS if t in s5.index]

    def spread(vals, gap):
        """nudge label y-positions apart, preserving order"""
        idx = sorted(range(len(vals)), key=lambda i: vals[i])
        out = list(vals)
        for k in range(1, len(idx)):
            a, b = idx[k - 1], idx[k]
            if out[b] - out[a] < gap:
                out[b] = out[a] + gap
        return out

    lft = spread([ATLAS[t][0] for t in present], 0.030)
    rgt = spread([s5.loc[t, "precision_vs_truth"] for t in present], 0.030)
    for t, ly, ry in zip(present, lft, rgt):
        axB.plot([0, 1], [ATLAS[t][0], s5.loc[t, "precision_vs_truth"]],
                 color=COLOR[t], lw=2.4, marker="o", ms=8)
        axB.text(-0.045, ly, LABEL[t], ha="right", va="center",
                 fontsize=8.8, color=COLOR[t])
        axB.text(1.045, ry, LABEL[t], ha="left", va="center",
                 fontsize=8.8, color=COLOR[t])
    dn = [t for t in DENOVO if t in s5.index]
    rho = spearmanr([ATLAS[t][0] for t in dn],
                    [s5.loc[t, "precision_vs_truth"] for t in dn]).statistic
    axB.set_xlim(-0.62, 1.62); axB.set_xticks([0, 1])
    axB.set_xticklabels(["PolyASite 2.0 atlas\n(manuscript/09)",
                         "Kinnex long-read truth\n(>=5 UMI, this figure)"],
                        fontsize=9, color=INK)
    axB.set_ylabel(f"precision @ {CUT} bp", fontsize=9.5, color=INK)
    axB.set_title(f"B  Precision ordering replicates (de novo Spearman rho = {rho:.2f})",
                  fontsize=10.5, color=INK, loc="left")

    # ---- C: what the calls actually sit on ----------------------------------
    axC = fig.add_subplot(gs[1, 0]); style(axC)
    d = dec.set_index("tool").reindex(
        [t for t in ["scutrquant", "scapture", "polyapipe", "sierra", "scapatrap", "peakatail"]
         if t in set(dec.tool)])
    y = np.arange(len(d))[::-1]
    left = np.zeros(len(d))
    for col, c, lab in ((("genuine_pas"), PALETTE[2], "genuine PAS  (long-read 3' end, clean downstream sequence)"),
                        (("ip_artifact"), PALETTE[1], "internal-priming 3' end  (long-read supported, but genomic A-rich)"),
                        (("no_lr_support"), "#B9C4C9", "no reproducible long-read 3' end  (>=5 UMI) nearby")):
        axC.barh(y, d[col], left=left, color=c, height=0.62, label=lab)
        for yy, (v, l) in zip(y, zip(d[col], left)):
            if v > 0.055:
                axC.text(l + v / 2, yy, f"{v*100:.0f}%", ha="center", va="center",
                         fontsize=8.4, color="white" if c != "#B9C4C9" else INK)
        left = left + d[col].values
    axC.set_yticks(y); axC.set_yticklabels([LABEL[t] for t in d.index], fontsize=9, color=INK)
    axC.set_xlim(0, 1); axC.set_xlabel("fraction of the tool's PAS calls", fontsize=9.5, color=INK)
    axC.xaxis.grid(True, color=GRID, lw=0.8); axC.yaxis.grid(False)
    axC.set_title("C  What each tool's calls actually sit on (100 bp, strand-matched)",
                  fontsize=10.5, color=INK, loc="left")
    axC.legend(frameon=False, fontsize=8.4, labelcolor=MUTED, loc="upper left",
               bbox_to_anchor=(0.0, -0.155), ncol=1, handlelength=1.4,
               handleheight=0.9, borderpad=0.0, labelspacing=0.42)

    # ---- D: F1 rank is not stable -------------------------------------------
    axD = fig.add_subplot(gs[1, 1]); style(axD)
    cols = ["atlas"] + [str(t) for t in THRESH]
    ranks = {}
    atl = sorted(DENOVO, key=lambda t: -ATLAS[t][2])
    ranks["atlas"] = {t: atl.index(t) + 1 for t in DENOVO}
    for T in THRESH:
        s = sc[(sc.truth == tag) & (sc.umi_threshold == T)].set_index("tool")
        o = sorted([t for t in DENOVO if t in s.index], key=lambda t: -s.loc[t, "f1_vs_truth"])
        ranks[str(T)] = {t: o.index(t) + 1 for t in o}
    xs = np.arange(len(cols))
    for t in DENOVO:
        ys = [len(DENOVO) + 1 - ranks[c].get(t, np.nan) for c in cols]
        axD.plot(xs, ys, color=COLOR[t], lw=2.4, marker="o", ms=7.5)
        axD.text(xs[-1] + 0.12, ys[-1], LABEL[t], fontsize=8.8, color=COLOR[t], va="center")
    axD.axvline(0.5, color=MUTED, lw=1.0, ls="--")
    axD.set_xlim(-0.35, len(cols) + 0.9)
    axD.set_xticks(xs)
    axD.set_xticklabels(["PolyASite\natlas"] + [f">={t} UMI" for t in THRESH],
                        fontsize=8.6, color=INK)
    axD.set_xlabel("long-read truth stringency", fontsize=9.5, color=INK)
    axD.set_yticks(range(1, len(DENOVO) + 1))
    axD.set_yticklabels([f"#{len(DENOVO) + 1 - i}" for i in range(1, len(DENOVO) + 1)],
                        fontsize=9, color=MUTED)
    axD.set_ylim(0.45, 5.5)
    axD.yaxis.grid(True, color=GRID, lw=0.8)
    axD.set_ylabel("F1 rank among de novo tools", fontsize=9.5, color=INK)
    npas = dict(zip(thq.umi_threshold, thq.n_pas))
    axD.set_title("D  The F1 ranking does NOT survive; only PeakATail's position does",
                  fontsize=10.5, color=INK, loc="left")
    axD.text(0.0, -0.235, "truth PAS:   "
             + "    ".join(f">={t}: {npas[t]:,}" for t in THRESH),
             transform=axD.transAxes, fontsize=8.0, color=MUTED)

    fig.suptitle("Atlas-independent validation of the pbmc_10k_v3 head-to-head against "
                 "PacBio Kinnex single-cell long reads", fontsize=13.5, color=INK,
                 x=0.072, ha="left", y=0.962)
    fig.text(0.072, 0.016,
             "CAVEAT  The Kinnex donor is NOT the pbmc_10k_v3 donor and the 10x chemistries differ (Kinnex GEM-X 3' v4 vs pbmc_10k_v3 3' v3): this is SITE-LEVEL truth and validates PAS POSITIONS, not per-cell usage.\n"
             "Truth = 104.0 M dedup FLNC molecules (isoseq refine --require-polya), real-cell CB, 3' softclip <=30 nt, internal priming = >=12 A or >=6 consecutive A in the 18 nt downstream, greedy +/-25 nt peak calling.\n"
             "F1 compares a fixed call set against truth sets that differ 53-fold in size, so it is not a stable ranking statistic here; precision and recall (panel B) are. scUTRquant* is annotation-based (fixed catalog) and is shown but not ranked with de novo tools.",
             fontsize=7.6, color=MUTED, ha="left", va="bottom")

    save_manuscript(fig, NAME, facecolor="white")
    plt.close(fig)

    # ---- every plotted value ------------------------------------------------
    out = OUTDIR / f"{NAME}.tsv"
    rankdf = pd.DataFrame([{"tool": t, **{f"f1_rank_{c}": ranks[c].get(t) for c in cols}}
                           for t in DENOVO])
    corr = []
    for T in THRESH:
        s = sc[(sc.truth == tag) & (sc.umi_threshold == T)].set_index("tool")
        dn = [t for t in DENOVO if t in s.index]
        corr.append(dict(umi_threshold=T,
                         spearman_precision=spearmanr([ATLAS[t][0] for t in dn],
                                                      [s.loc[t, "precision_vs_truth"] for t in dn]).statistic,
                         spearman_recall=spearmanr([ATLAS[t][1] for t in dn],
                                                   [s.loc[t, "recall_of_truth"] for t in dn]).statistic,
                         spearman_f1=spearmanr([ATLAS[t][2] for t in dn],
                                               [s.loc[t, "f1_vs_truth"] for t in dn]).statistic))
    with open(out, "w") as fh:
        fh.write(f"# {NAME}: Kinnex long-read truth vs pbmc_10k_v3 tool calls. truth tag(s): "
                 + ",".join(tags) + "\n")
        fh.write("# CAVEAT: different donor and 10x chemistry from pbmc_10k_v3 -- site-level truth.\n")
        fh.write("# panel B/D source: per-tool precision, recall, F1 at each truth stringency\n")
        sc.to_csv(fh, sep="\t", index=False, float_format="%.6f")
        fh.write("\n# panel A: sequence-only truth calibration (hexamer vs support)\n")
        cal.to_csv(fh, sep="\t", index=False, float_format="%.6f")
        fh.write("\n# panel C: long-read support decomposition of each tool's calls\n")
        dec.to_csv(fh, sep="\t", index=False, float_format="%.6f")
        fh.write("\n# truth-set size / independent atlas agreement per stringency\n")
        thq.to_csv(fh, sep="\t", index=False, float_format="%.6f")
        fh.write("\n# panel D: F1 ranks among de novo tools\n")
        rankdf.to_csv(fh, sep="\t", index=False)
        fh.write("\n# rank correlation with the PolyASite-2.0 head-to-head (de novo tools only)\n")
        pd.DataFrame(corr).to_csv(fh, sep="\t", index=False, float_format="%.4f")
    print("wrote", out)
    print(pd.DataFrame(corr).to_string(index=False))
    print(rankdf.to_string(index=False))


if __name__ == "__main__":
    main()
