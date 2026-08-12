#!/usr/bin/env python3
"""Manuscript figure: peak-calling strategy benchmark vs PolyASite 3.0.

Panels
  (a) precision vs distance cutoff, one line per arm (log-x), with the
      random-placement chance floor drawn in (width-matched genome shuffle).
  (b) number of PAS called per arm (bar).
  (c) recall vs distance cutoff (log-x).

HONESTY: the reference is 18,432,135 PolyASite 3.0 3'-end clusters (NOT
single nucleotides: mean 9.0 bp wide, 165 Mb in total) spanning the whole
genome. "Recall" is therefore coverage of the entire atlas (including genes
not expressed in lung), not sensitivity, and it tracks the total sequence an
arm searches -- n_peaks x (peak width + 2 x cutoff) -- rather than the number
of peaks called: lambda_poisson calls MORE peaks than the IP-filtered arm yet
covers LESS atlas at every cutoff <= 1 kb. Precision is likewise inflated by
reference density: atlas +/-50 bp windows already tile 38.6% of the genome,
and called peaks are 0.4-0.6 kb wide, so width-matched RANDOM peaks score
~0.67-0.74 precision at +/-50 bp. All of this is annotated on the panels.

Inputs are read-only from /mnt/ssd2. Outputs:
  results/figures/manuscript/benchmark_strategies.png (300 dpi)
  results/figures/manuscript/benchmark_strategies.tsv (every plotted value)
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths ----
GRID = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/"
            "RERUN_2026-08_fixed/runs/grid")
ATLAS = Path("/mnt/ssd2/Laugney_Aligned/refs/"
             "polyasite_3.0_GRCh38_ensembl_sorted.bed")
BAM_FOR_GENOME = GRID / "lg_annotate" / "GSM3516663-StageIA_merged_sorted.bam"
OUTDIR = Path("/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript")
CACHE = OUTDIR / ".cache_benchmark_strategies"
NAME = "benchmark_strategies"

CUTOFFS = [50, 100, 200, 500, 1000, 5000]
SHUFFLE_SEEDS = [1, 2, 3]

# ------------------------------------------------------- design tokens ----
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRIDC, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

# One colour per ENTITY, fixed across every panel.
ARMS = [
    # key             colour  legend label                    bar label
    ("lg_annotate",   C[0], "lambda_gradient (default)",  "lambda_gradient\n(default)"),
    ("lp_annotate",   C[1], "lambda_poisson",             "lambda_poisson"),
    ("si_annotate",   C[2], "sierra_iterative",           "sierra_iterative"),
    ("lg_ip_filter",  C[3], "lambda_gradient + IP filter", "lambda_gradient\n(IP filter ON)"),
    ("lg_ip_off",     C[4], "lambda_gradient, IP off",    "lambda_gradient\n(IP filter OFF)"),
]
ORDER = [a[0] for a in ARMS]
COLOR = {a[0]: a[1] for a in ARMS}
LEGLAB = {a[0]: a[2] for a in ARMS}
BARLAB = {a[0]: a[3] for a in ARMS}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.linewidth": 0.8,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
})


def sh(cmd: str) -> str:
    r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {cmd}\n{r.stderr[:2000]}")
    return r.stdout


# ------------------------------------------------------------- load in ----
def load_benchmarks() -> dict:
    out = {}
    for arm in ORDER:
        p = GRID / arm / "benchmark_vs_polyasite_v3.json"
        with open(p) as fh:
            out[arm] = json.load(fh)
    return out


def bed_facts() -> pd.DataFrame:
    """Line count (confirms n_predicted), md5 (detects duplicate arms), widths."""
    rows = []
    for arm in ORDER:
        bed = GRID / arm / "pasbed.bed"
        n = int(sh(f"wc -l < {bed}").strip())
        md5 = sh(f"md5sum {bed}").split()[0]
        w = pd.read_csv(bed, sep="\t", header=None, usecols=[1, 2],
                        names=["start", "end"])
        width = (w["end"] - w["start"]).to_numpy()
        rows.append(dict(arm=arm, bed_lines=n, md5=md5,
                         median_width_bp=int(np.median(width)),
                         mean_width_bp=round(float(width.mean()), 1),
                         total_peak_bp=int(width.sum())))
    return pd.DataFrame(rows).set_index("arm")


def atlas_site_stats() -> dict:
    """Width of the atlas features themselves.

    The atlas is NOT a set of single nucleotides: PolyASite 3.0 entries are
    3'-end clusters. Read the real widths so the figure text cannot drift.
    """
    cache = CACHE / "atlas_site_stats.tsv"
    if not cache.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        sh("awk -F'\\t' '{w=$3-$2; s+=w; n++; if(w==1)o++} "
           "END{printf \"%d\\t%.4f\\t%.4f\\t%d\\n\", n, s/n, o/n, s}' "
           f"{ATLAS} > {cache}")
    n, mean_w, frac1, total = cache.read_text().split()
    return dict(n=int(n), mean_width_bp=float(mean_w),
                frac_single_nt=float(frac1), total_bp=int(total))


def atlas_union_coverage() -> pd.DataFrame:
    """Fraction of the 24 atlas chromosomes covered by union(atlas site +/- c).

    This is the probability that a single RANDOM nucleotide 'matches' the
    reference at each cutoff -- the density-inflation term behind precision.
    """
    cache = CACHE / "atlas_union_cov.tsv"
    if not cache.exists():
        cus = ",".join(map(str, CUTOFFS))
        awk = r'''BEGIN{nc=split(CUS,CU,",")}
{ for(i=1;i<=nc;i++){ c=CU[i]; s=$2-c; if(s<0)s=0; e=$3+c;
    if($1!=pc[i]){ if(pc[i]!=""){tot[i]+=pe[i]-ps[i]} pc[i]=$1; ps[i]=s; pe[i]=e; continue }
    if(s>pe[i]){ tot[i]+=pe[i]-ps[i]; ps[i]=s; pe[i]=e } else { if(e>pe[i]) pe[i]=e } } }
END{ for(i=1;i<=nc;i++){ tot[i]+=pe[i]-ps[i]; print CU[i]"\t"tot[i] } }'''
        sh(f"awk -v CUS='{cus}' '{awk}' {ATLAS} > {cache}")
    cov = pd.read_csv(cache, sep="\t", header=None,
                      names=["cutoff_bp", "covered_bp"])
    gbp = int(sh(f"awk 'NR==FNR{{ok[$1]=1;next}} ok[$1]{{s+=$2}} "
                 f"END{{print s}}' {CACHE/'atlas_chroms.txt'} "
                 f"{CACHE/'genome.txt'}").strip())
    cov["genome_bp"] = gbp
    cov["frac_genome_within_cutoff_of_atlas"] = cov["covered_bp"] / gbp
    return cov


def prep_genome_files() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    g, ac, ga = CACHE / "genome.txt", CACHE / "atlas_chroms.txt", CACHE / "genome_atlas.txt"
    if not g.exists():
        sh(f"samtools view -H {BAM_FOR_GENOME} | awk '$1==\"@SQ\"{{c=\"\";l=\"\";"
           f"for(i=2;i<=NF;i++){{if($i~/^SN:/){{c=substr($i,4)}};"
           f"if($i~/^LN:/){{l=substr($i,4)}}}}; print c\"\\t\"l}}' > {g}")
    if not ac.exists():
        sh(f"cut -f1 {ATLAS} | uniq > {ac}")
    if not ga.exists():
        sh(f"awk 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {ac} {g} | sort -k1,1 > {ga}")
    for arm in ORDER:
        s = CACHE / f"{arm}.sorted.bed"
        if not s.exists():
            sh(f"sort -k1,1 -k2,2n {GRID/arm/'pasbed.bed'} > {s}")


def null_precision() -> pd.DataFrame:
    """Chance floor: shuffle the SAME peaks (widths preserved) genome-wide,
    then score them against the atlas exactly as the real benchmark does."""
    cache = CACHE / "null_precision.tsv"
    if not cache.exists():
        prep_genome_files()
        lines = []
        for arm in ORDER:
            for seed in SHUFFLE_SEEDS:
                d = CACHE / f"nulldist_{arm}_{seed}.txt"
                if not d.exists():
                    sh(f"bedtools shuffle -i {CACHE/f'{arm}.sorted.bed'} "
                       f"-g {CACHE/'genome_atlas.txt'} -seed {seed} "
                       f"| sort -k1,1 -k2,2n "
                       f"| bedtools closest -a - -b {ATLAS} -d -t first "
                       f"-g {CACHE/'genome_atlas.txt'} "
                       f"| awk -F'\\t' '{{print $NF}}' > {d}")
                dist = pd.read_csv(d, header=None).iloc[:, 0].to_numpy()
                for c in CUTOFFS:
                    m = int(((dist >= 0) & (dist <= c)).sum())
                    lines.append(f"{arm}\t{seed}\t{c}\t{m}\t{len(dist)}")
        cache.write_text("\n".join(lines) + "\n")
    df = pd.read_csv(cache, sep="\t", header=None,
                     names=["arm", "seed", "cutoff_bp", "matched", "n"])
    df["null_precision"] = df["matched"] / df["n"]
    return df


def verify_precision_reproduces(bm: dict) -> str:
    """Independently recompute precision with bedtools closest and compare."""
    prep_genome_files()
    msgs = []
    for arm in ORDER:
        d = CACHE / f"realdist_{arm}.txt"
        if not d.exists():
            sh(f"bedtools closest -a {CACHE/f'{arm}.sorted.bed'} -b {ATLAS} "
               f"-d -t first -g <(sort -k1,1 {CACHE/'genome.txt'}) 2>/dev/null "
               f"| awk -F'\\t' '{{print $NF}}' > {d}")
        dist = pd.read_csv(d, header=None).iloc[:, 0].to_numpy()
        for c in CUTOFFS:
            mine = int(((dist >= 0) & (dist <= c)).sum())
            theirs = bm[arm]["cutoffs"][str(c)]["matched_predicted"]
            if mine != theirs:
                msgs.append(f"{arm}@{c}: recomputed {mine} != json {theirs}")
    return "OK" if not msgs else "; ".join(msgs)


# ---------------------------------------------------------------- plot ----
def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    bm = load_benchmarks()
    facts = bed_facts()
    nullp = null_precision()
    cov = atlas_union_coverage()
    astats = atlas_site_stats()
    check = verify_precision_reproduces(bm)

    # Atlas coverage ("recall") tracks the total sequence an arm searches,
    # n_peaks * (peak width + 2*cutoff), NOT the number of peaks called.
    # Spearman rank correlation of that footprint with panel-c recall, per
    # cutoff, over the 4 distinct arms (lg_ip_off duplicates lg_annotate).
    uniq_arms = [a for a in ORDER if a != "lg_ip_off"]
    rho_fp, rho_n = {}, {}
    for c in CUTOFFS:
        fp = np.array([facts.loc[a, "total_peak_bp"]
                       + facts.loc[a, "bed_lines"] * 2 * c for a in uniq_arms],
                      dtype=float)
        nn = np.array([facts.loc[a, "bed_lines"] for a in uniq_arms], dtype=float)
        rc = np.array([bm[a]["cutoffs"][str(c)]["recall"] for a in uniq_arms])
        rank = lambda v: pd.Series(v).rank().to_numpy()
        rho_fp[c] = float(np.corrcoef(rank(fp), rank(rc))[0, 1])
        rho_n[c] = float(np.corrcoef(rank(nn), rank(rc))[0, 1])

    n_ref = bm[ORDER[0]]["n_reference"]
    assert all(bm[a]["n_reference"] == n_ref for a in ORDER)

    # duplicate-arm detection (md5 of pasbed.bed)
    dup_pairs = [list(idx) for idx in facts.groupby("md5").groups.values()
                 if len(idx) > 1]

    # mean null across seeds, per arm
    nmean = (nullp.groupby(["arm", "cutoff_bp"])["null_precision"]
             .mean().unstack())
    null_lo = nmean.min(axis=0).reindex(CUTOFFS).to_numpy()
    null_hi = nmean.max(axis=0).reindex(CUTOFFS).to_numpy()

    x = np.array(CUTOFFS, dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 5.9))
    fig.patch.set_facecolor(SURFACE)
    ax_a, ax_b, ax_c = axes

    def style(ax, ylab):
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRIDC, linewidth=0.6)
        ax.xaxis.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel(ylab, fontsize=8.8)

    # ---- (a) precision -----------------------------------------------
    band = ax_a.fill_between(x, 0, null_hi, color=GRIDC, alpha=0.9, zorder=1)
    ax_a.plot(x, null_hi, color=MUTED, lw=1.2, ls=(0, (4, 2)), zorder=2)
    handles, labels = [], []
    for arm in ORDER:
        y = [bm[arm]["cutoffs"][str(c)]["precision"] for c in CUTOFFS]
        ls = (0, (1.2, 1.8)) if arm == "lg_ip_off" else "-"
        ln, = ax_a.plot(x, y, color=COLOR[arm], lw=1.8, ls=ls, marker="o",
                        ms=6, mec=SURFACE, mew=0.8, zorder=4)
        handles.append(ln)
        labels.append(LEGLAB[arm])
    handles.append(band)
    labels.append("chance floor: same peaks shuffled genome-wide (3 seeds)")
    ax_a.set_xscale("log")
    ax_a.set_xticks(CUTOFFS)
    ax_a.set_xticklabels([str(c) for c in CUTOFFS], fontsize=8)
    ax_a.set_xlim(40, 6400)
    ax_a.set_ylim(0.55, 1.025)
    ax_a.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    style(ax_a, "Precision (called PAS within cutoff\nof any atlas site)")
    ax_a.set_xlabel("Distance cutoff (bp, log scale)", fontsize=8.8)
    ax_a.set_title("a  Every strategy saturates precision \u2014 but so does\n"
                   "    randomly shuffled noise", fontsize=9.5, loc="left",
                   color=INK, pad=8)
    ax_a.annotate(f"peaks shuffled to random genomic\npositions still score "
                  f"{null_hi[0]:.2f} at $\\pm$50 bp",
                  xy=(50, null_hi[0]), xytext=(52, 0.885),
                  fontsize=7.2, color=MUTED, ha="left", va="bottom",
                  linespacing=1.4,
                  arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7,
                                  shrinkA=2, shrinkB=3))
    ax_a.text(0.98, 0.455, "zoom \u2014 the arms differ only in the 3rd decimal",
              transform=ax_a.transAxes, ha="right", va="bottom",
              fontsize=7.0, color=MUTED)

    # zoom inset: the arms only separate in the 3rd decimal
    axi = ax_a.inset_axes([0.44, 0.15, 0.54, 0.28])
    axi.set_facecolor(SURFACE)
    for arm in ORDER:
        y = [bm[arm]["cutoffs"][str(c)]["precision"] for c in CUTOFFS]
        ls = (0, (1.2, 1.8)) if arm == "lg_ip_off" else "-"
        axi.plot(x, y, color=COLOR[arm], lw=1.4, ls=ls, marker="o", ms=3.4,
                 mec=SURFACE, mew=0.5)
    axi.set_xscale("log")
    axi.set_ylim(0.9952, 1.0004)
    axi.set_yticks([0.996, 0.998, 1.000])
    axi.set_yticklabels(["0.996", "0.998", "1.000"], fontsize=6)
    axi.set_xticks([50, 500, 5000])
    axi.set_xticklabels(["50", "500", "5000"], fontsize=6)
    axi.tick_params(length=2, pad=1.5)
    axi.set_axisbelow(True)
    axi.yaxis.grid(True, color=GRIDC, linewidth=0.5)
    for s in ("top", "right"):
        axi.spines[s].set_visible(False)

    # ---- (b) PAS called ----------------------------------------------
    bx = np.arange(len(ORDER))
    counts = [bm[a]["n_predicted"] for a in ORDER]
    ax_b.bar(bx, counts, width=0.66, color=[COLOR[a] for a in ORDER], zorder=3)
    for i, v in enumerate(counts):
        ax_b.text(i, v + 900, f"{v:,}", ha="center", va="bottom",
                  fontsize=7.8, color=INK)
    ax_b.set_xticks(bx)
    ax_b.set_xticklabels([LEGLAB[a] for a in ORDER], fontsize=7.2, color=INK,
                         rotation=26, ha="right", rotation_mode="anchor")
    ax_b.set_ylim(0, 57000)
    ax_b.set_yticks([0, 10000, 20000, 30000, 40000, 50000])
    ax_b.set_yticklabels(["0", "10,000", "20,000", "30,000", "40,000",
                          "50,000"], fontsize=8)
    style(ax_b, "PAS called (n peaks, 17-sample cohort)")
    n_si = bm["si_annotate"]["n_predicted"]
    r_lo = min(n_si / bm[a]["n_predicted"] for a in ("lg_annotate", "lp_annotate"))
    r_hi = max(n_si / bm[a]["n_predicted"] for a in ("lg_annotate", "lp_annotate"))
    ip_drop = 100 * (1 - bm["lg_ip_filter"]["n_predicted"]
                     / bm["lg_annotate"]["n_predicted"])
    ax_b.set_title(f"b  sierra_iterative calls {r_lo:.2f}–{r_hi:.2f}× as many PAS\n"
                   f"    as the lambda methods; the IP filter removes {ip_drop:.1f}%",
                   fontsize=9.5, loc="left", color=INK, pad=8)
    if dup_pairs:
        pair = " and ".join(dup_pairs[0])
        ax_b.text(0.99, 0.985,
                  f"anomaly: {pair}\n"
                  "are the same file (identical md5), so the\n"
                  "'IP off' arm merely reproduces the default.\n"
                  "The only real IP contrast is ON vs OFF:\n"
                  "20,609 vs 22,633 peaks (-2,024, -8.9%)",
                  transform=ax_b.transAxes, ha="right", va="top",
                  fontsize=7.0, color=MUTED, linespacing=1.5)

    # ---- (c) recall ---------------------------------------------------
    for arm in ORDER:
        y = [bm[arm]["cutoffs"][str(c)]["recall"] for c in CUTOFFS]
        ls = (0, (1.2, 1.8)) if arm == "lg_ip_off" else "-"
        ax_c.plot(x, y, color=COLOR[arm], lw=1.8, ls=ls, marker="o", ms=6,
                  mec=SURFACE, mew=0.8, zorder=4)
    # direct labels (lg_annotate and lg_ip_off coincide -> one shared label)
    lab_at = {"si_annotate": ("sierra_iterative", 0.000),
              "lg_annotate": ("lambda_gradient (default),\nidentical to 'IP off' arm", 0.012),
              "lg_ip_filter": ("lambda_gradient + IP filter", -0.011),
              "lp_annotate": ("lambda_poisson", -0.028)}
    for arm, (txt, dy) in lab_at.items():
        yv = bm[arm]["cutoffs"]["5000"]["recall"]
        ax_c.annotate(txt, xy=(5000, yv), xytext=(6300, yv + dy),
                      color=COLOR[arm], fontsize=7.2, va="center", ha="left",
                      linespacing=1.35)
    ax_c.set_xscale("log")
    ax_c.set_xticks(CUTOFFS)
    ax_c.set_xticklabels([str(c) for c in CUTOFFS], fontsize=8)
    ax_c.set_xlim(40, 62000)
    ax_c.set_ylim(0, 0.30)
    style(ax_c, "Fraction of the 18.4 M-cluster atlas covered\n(this is NOT sensitivity)")
    ax_c.set_xlabel("Distance cutoff (bp, log scale)", fontsize=8.8)
    ax_c.set_title("c  'Recall' here tracks how much sequence each arm\n"
                   "    searches \u2014 it is not a sensitivity estimate",
                   fontsize=9.5, loc="left", color=INK, pad=8)
    ax_c.text(0.02, 0.985,
              f"denominator = {n_ref:,} PolyASite 3.0\n"
              f"3'-end clusters (mean {astats['mean_width_bp']:.0f} bp wide, "
              f"{astats['total_bp'] / 1e6:.0f} Mb\n"
              "in total) across the whole genome,\n"
              "including genes not expressed in lung.\n"
              "A perfect lung PAS caller could not\n"
              "approach 1.0 on this axis.",
              transform=ax_c.transAxes, ha="left", va="top",
              fontsize=7.0, color=MUTED, linespacing=1.5)
    n_lp, n_ipf = bm["lp_annotate"]["n_predicted"], bm["lg_ip_filter"]["n_predicted"]
    ax_c.text(0.025, 0.62,
              "...and not even peak COUNT:\n"
              "lambda_poisson calls MORE peaks\n"
              f"({n_lp:,}) than the IP-filtered arm\n"
              f"({n_ipf:,}) yet covers LESS atlas.\n"
              "It tracks searched BASES\n"
              "(foot-note 1).",
              transform=ax_c.transAxes, ha="left", va="top",
              fontsize=6.9, color=MUTED, linespacing=1.5)

    # ---- shared legend + caveat strip --------------------------------
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, 0.128), ncol=6, frameon=False,
               fontsize=7.2, handlelength=2.0, columnspacing=1.5,
               handletextpad=0.6)

    wmin = int(facts["median_width_bp"].min())
    wmax = int(facts["median_width_bp"].max())
    c50 = cov.loc[cov.cutoff_bp == 50, "frac_genome_within_cutoff_of_atlas"].iloc[0]
    c5k = cov.loc[cov.cutoff_bp == 5000, "frac_genome_within_cutoff_of_atlas"].iloc[0]
    fig.text(
        0.008, 0.105,
        "What this benchmark cannot show  |  "
        f"(1) The reference is {n_ref:,} PolyASite 3.0 3'-end clusters (mean {astats['mean_width_bp']:.0f} bp wide, {astats['total_bp'] / 1e6:.0f} Mb), genome-wide, so panel c is atlas coverage, not sensitivity: across arms it tracks the bases searched — rank corr. {min(rho_fp.values()):.2f} at every cutoff vs {min(rho_n.values()):.2f}–{max(rho_n.values()):.2f} for peak count.\n"
        f"(2) Precision is inflated by reference density \u2014 atlas $\\pm$50 bp windows already tile {c50 * 100:.1f}% of the genome ({c5k * 100:.1f}% at $\\pm$5 kb) \u2014 and called peaks are {wmin}\u2013{wmax} bp wide (median), so the effective match window is peak width + 2\u00d7cutoff, not 2\u00d7cutoff.\n"
        f"(3) Precision of ~0.999 must therefore be read against the {null_hi[0]:.2f} chance floor in panel a; it is NOT a clean accuracy claim. (4) Matching is strand-agnostic (bedtools window, no -s). (5) One cohort-level run per arm: no replicates, no CI on the arm curves.",
        fontsize=6.5, color=MUTED, ha="left", va="top", linespacing=1.7)

    fig.tight_layout(rect=(0, 0.175, 0.995, 1.0))
    png = OUTDIR / f"{NAME}.png"
    fig.savefig(png, dpi=300, facecolor=SURFACE)
    plt.close(fig)

    # ----------------------------------------------------------- tsv ----
    rows = []
    for arm in ORDER:
        b = bm[arm]
        rows.append(dict(panel="b", arm=arm, arm_label=LEGLAB[arm],
                         series="n_pas_called", cutoff_bp="", value=b["n_predicted"],
                         numerator=b["n_predicted"], denominator="",
                         bed_lines=facts.loc[arm, "bed_lines"],
                         pasbed_md5=facts.loc[arm, "md5"],
                         median_peak_width_bp=facts.loc[arm, "median_width_bp"],
                         mean_peak_width_bp=facts.loc[arm, "mean_width_bp"],
                             total_peak_bp=facts.loc[arm, "total_peak_bp"]))
        for c in CUTOFFS:
            m = b["cutoffs"][str(c)]
            rows.append(dict(panel="a", arm=arm, arm_label=LEGLAB[arm],
                             series="precision", cutoff_bp=c,
                             value=m["precision"],
                             numerator=m["matched_predicted"],
                             denominator=b["n_predicted"],
                             bed_lines=facts.loc[arm, "bed_lines"],
                             pasbed_md5=facts.loc[arm, "md5"],
                             median_peak_width_bp=facts.loc[arm, "median_width_bp"],
                             mean_peak_width_bp=facts.loc[arm, "mean_width_bp"],
                             total_peak_bp=facts.loc[arm, "total_peak_bp"]))
            rows.append(dict(panel="c", arm=arm, arm_label=LEGLAB[arm],
                             series="recall_atlas_coverage", cutoff_bp=c,
                             value=m["recall"],
                             numerator=m["matched_reference"],
                             denominator=b["n_reference"],
                             bed_lines=facts.loc[arm, "bed_lines"],
                             pasbed_md5=facts.loc[arm, "md5"],
                             median_peak_width_bp=facts.loc[arm, "median_width_bp"],
                             mean_peak_width_bp=facts.loc[arm, "mean_width_bp"],
                             total_peak_bp=facts.loc[arm, "total_peak_bp"]))
            sub = nullp[(nullp.arm == arm) & (nullp.cutoff_bp == c)]
            rows.append(dict(panel="a", arm=arm, arm_label=LEGLAB[arm],
                             series="null_precision_shuffled_mean3seeds",
                             cutoff_bp=c,
                             value=round(float(sub["null_precision"].mean()), 4),
                             numerator=int(sub["matched"].mean().round()),
                             denominator=b["n_predicted"],
                             bed_lines=facts.loc[arm, "bed_lines"],
                             pasbed_md5=facts.loc[arm, "md5"],
                             median_peak_width_bp=facts.loc[arm, "median_width_bp"],
                             mean_peak_width_bp=facts.loc[arm, "mean_width_bp"],
                             total_peak_bp=facts.loc[arm, "total_peak_bp"]))
    for _, r in cov.iterrows():
        rows.append(dict(panel="a", arm="(reference)", arm_label="PolyASite 3.0 atlas",
                         series="frac_genome_within_cutoff_of_any_atlas_site",
                         cutoff_bp=int(r["cutoff_bp"]),
                         value=round(float(r["frac_genome_within_cutoff_of_atlas"]), 4),
                         numerator=int(r["covered_bp"]),
                         denominator=int(r["genome_bp"]),
                         bed_lines="", pasbed_md5="",
                         median_peak_width_bp="", mean_peak_width_bp="",
                         total_peak_bp=""))
    tsv = pd.DataFrame(rows, columns=[
        "panel", "arm", "arm_label", "series", "cutoff_bp", "value",
        "numerator", "denominator", "bed_lines", "pasbed_md5",
        "median_peak_width_bp", "mean_peak_width_bp", "total_peak_bp"])
    tsv_path = OUTDIR / f"{NAME}.tsv"
    tsv.to_csv(tsv_path, sep="\t", index=False)

    print(f"n_reference = {n_ref:,}")
    print(f"pasbed line-count vs n_predicted: "
          f"{'ALL MATCH' if all(facts.loc[a,'bed_lines']==bm[a]['n_predicted'] for a in ORDER) else 'MISMATCH'}")
    print(f"precision recomputation check: {check}")
    print(f"duplicate pasbed arms: {dup_pairs}")
    print(facts.to_string())
    print(cov.to_string(index=False))
    print(nmean.round(4).to_string())
    print(f"wrote {png}\nwrote {tsv_path}")


if __name__ == "__main__":
    main()
