#!/usr/bin/env python3
"""
Figure: cohort QC / PAS landscape for the 17-sample Laughney lung-cancer cohort.

Panels
  (a) cells retained per sample, coloured by clinical stage category
  (b) PAS-per-gene distribution (only multi-PAS genes can yield a classic PDUI)
  (c) PAS annotation TIER distribution (cohort PAS set vs all annotated PAS)
  (d) PAS width distribution (end - start)

All sources under /mnt/ssd2 are READ-ONLY.
Outputs: <OUT>/cohort_qc.png (300 dpi) and <OUT>/cohort_qc.tsv (every plotted value).
"""

import json
import os
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import blended_transform_factory

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


# ----------------------------------------------------------------------------- paths
RUN = ("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/"
       "runs/B1_cohort_full")
PASBED = os.path.join(RUN, "pasbed.bed")
ANNBED = os.path.join(RUN, "annotatedpas.bed")
UTRLEN = os.path.join(RUN, "utr_lengths.tsv")
CLUST = os.path.join(RUN, "07_clustering")
AMAT = os.path.join(RUN, "05_annotated_matrix")
# The un-annotated called-peak universe: every peak that entered gene assignment.
POSBED = os.path.join(RUN, "posbed.bed")
NEGBED = os.path.join(RUN, "negbed.bed")

OUT = "/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript"
os.makedirs(OUT, exist_ok=True)
PNG = os.path.join(OUT, "cohort_qc.png")
TSV = os.path.join(OUT, "cohort_qc.tsv")

# ----------------------------------------------------------------------------- style
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK = "#1B2429"
MUTED = "#5A6B73"
GRID = "#D8E0E3"
SURFACE = "#FFFFFF"

# Entity -> colour, fixed across every panel of this figure.
CAT_ORDER = ["Adjacent normal", "Primary, stage I-II", "Primary, stage IV",
             "Metastasis"]
CAT_SHORT = {CAT_ORDER[0]: "Normal", CAT_ORDER[1]: "Stage I-II",
             CAT_ORDER[2]: "Stage\nIV", CAT_ORDER[3]: "Metastasis"}
CAT_COLOR = {CAT_ORDER[0]: C[0], CAT_ORDER[1]: C[1],
             CAT_ORDER[2]: C[2], CAT_ORDER[3]: C[3]}
COHORT_PAS_COLOR = C[0]        # the 16.5k cohort PAS set, panels b/c/d
ALL_PAS_COLOR = C[5]           # the full annotated PAS universe, panel c
UNUSABLE_COLOR = MUTED         # single-PAS genes: cannot yield a PDUI

plt.rcParams.update({
    "font.size": 7.4,
    "axes.titlesize": 8.2,
    "axes.labelsize": 7.6,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.fontsize": 6.4,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.linewidth": 0.8,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.tick_params(length=2.5, width=0.8)


def stage_category(sample_dir):
    """Map a GSM directory name (GSM3516662-StageIA) to a clinical category."""
    suffix = sample_dir.split("-", 1)[1] if "-" in sample_dir else sample_dir
    if suffix.lower().startswith("normal"):
        return CAT_ORDER[0]
    if suffix.lower().startswith("met"):
        return CAT_ORDER[3]
    m = re.match(r"Stage(IV|III|II|I)", suffix)
    if m:
        return CAT_ORDER[2] if m.group(1) == "IV" else CAT_ORDER[1]
    raise ValueError("unparsed stage in %s" % sample_dir)


# ============================================================ (a) cells per sample
rows = []
for d in sorted(os.listdir(CLUST)):
    sd = os.path.join(CLUST, d)
    if not os.path.isdir(sd):
        continue
    with open(os.path.join(sd, "clustering_stats.json")) as fh:
        stats = json.load(fh)
    with open(os.path.join(sd, "filtered_cb.tsv")) as fh:
        n_cb = sum(1 for line in fh if line.strip())
    rows.append({
        "sample": d,
        "gsm": d.split("-")[0],
        "stage_label": d.split("-", 1)[1],
        "category": stage_category(d),
        "cells_after_cb_filter": n_cb,
        "cells_retained": int(stats["final_cells"]),
        "pas_retained": int(stats["final_pas"]),
    })
samples = pd.DataFrame(rows)
samples["retention_pct"] = (100.0 * samples["cells_retained"]
                            / samples["cells_after_cb_filter"])
samples["cat_rank"] = samples["category"].map({c: i for i, c in enumerate(CAT_ORDER)})
samples = samples.sort_values(["cat_rank", "cells_retained"],
                              ascending=[True, False]).reset_index(drop=True)

n_samples = len(samples)
cells_min = int(samples["cells_retained"].min())
cells_max = int(samples["cells_retained"].max())
cells_total = int(samples["cells_retained"].sum())
median_retention = float(samples["retention_pct"].median())
min_retention = float(samples["retention_pct"].min())

# ================================================== cohort PAS set + annotation
pas = pd.read_csv(PASBED, sep="\t", header=None,
                  names=["chrom", "start", "end", "pas_id", "score", "strand"],
                  dtype={"chrom": str})
pas["width"] = pas["end"] - pas["start"]
n_pas = len(pas)

ann = pd.read_csv(ANNBED, sep="\t", header=None,
                  names=["chrom", "start", "end", "pas_id", "gene_id",
                         "gene_symbol", "strand", "score", "tier"],
                  usecols=["pas_id", "gene_id", "gene_symbol", "tier"],
                  dtype={"pas_id": np.int64})
n_ann = len(ann)
ann_tier_all = ann["tier"].value_counts()

cohort = pas.merge(ann, on="pas_id", how="left", validate="one_to_one")
n_unannotated = int(cohort["gene_id"].isna().sum())
cohort_tier = cohort["tier"].value_counts()


def _count_lines(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


# The honest denominator for panel c. find_close.py keeps TIER_1 + TIER_2 only
# (include_extended=False by default), so annotatedpas.bed cannot contain a
# TIER_3 or INTERGENIC peak: quoting a tier purity against it is circular.
# posbed + negbed are the peaks as called, before any gene/tier filtering.
n_called = _count_lines(POSBED) + _count_lines(NEGBED)
n_dropped = n_called - n_ann          # TIER_3 + INTERGENIC + no gene within range

# What the 16,500-PAS "cohort set" actually is: the strict 17/17 intersection of
# the per-sample PAS sets, not a union and not the pipeline's PAS output.
detect = {}
for d in sorted(os.listdir(AMAT)):
    f = os.path.join(AMAT, d, "annotated_pas_ids.tsv")
    if not os.path.isfile(f):
        continue
    for pid in pd.read_csv(f, sep="\t")["pas_id"]:
        detect[pid] = detect.get(pid, 0) + 1
n_ids_any = len(detect)
n_samples_detect = len({d for d in os.listdir(AMAT)
                        if os.path.isfile(os.path.join(AMAT, d,
                                                       "annotated_pas_ids.tsv"))})
ids_all_samples = {p for p, c in detect.items() if c == n_samples_detect}
# Hard check: the plotted "cohort PAS set" IS that intersection, exactly.
assert ids_all_samples == set(pas["pas_id"]), (
    "pasbed.bed is not the %d/%d intersection (%d vs %d)"
    % (n_samples_detect, n_samples_detect, len(ids_all_samples), len(pas)))
cohort_frac_of_any = n_pas / n_ids_any

# ============================================================ (b) PAS per gene
per_gene = (cohort.dropna(subset=["gene_id"])
            .groupby("gene_id")["pas_id"].size().rename("n_pas"))
n_genes = int(per_gene.size)
b_labels = ["1", "2", "3", "4+"]
b_counts = [int((per_gene == 1).sum()), int((per_gene == 2).sum()),
            int((per_gene == 3).sum()), int((per_gene >= 4).sum())]
n_multi = n_genes - b_counts[0]
multi_frac = n_multi / n_genes
max_pas_per_gene = int(per_gene.max())

# 3'UTR length context for the panel-b caveat (single- vs multi-PAS genes)
utr = pd.read_csv(UTRLEN, sep="\t").set_index("gene_id")["utr_length"]
utr_single_med = float(utr.reindex(per_gene[per_gene == 1].index).dropna().median())
utr_multi_med = float(utr.reindex(per_gene[per_gene >= 2].index).dropna().median())

# ============================================================ (d) PAS widths
w = cohort["width"].to_numpy()
w_med = float(np.median(w))
bins = np.logspace(np.log10(max(w.min(), 1)), np.log10(w.max() * 1.02), 31)
w_hist, w_edges = np.histogram(w, bins=bins)

# ================================================================ figure
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
ax_a, ax_b, ax_c, ax_d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
caveats = {}

# ---------------------------------------------------------------------- (a)
x = np.arange(n_samples)
yv = samples["cells_retained"].to_numpy()
ax_a.bar(x, yv, width=0.76, color=[CAT_COLOR[c] for c in samples["category"]],
         edgecolor="none")
style_axis(ax_a)
ax_a.set_xticks(x)
ax_a.set_xticklabels([g.replace("GSM35166", "") for g in samples["gsm"]],
                     rotation=90)
ax_a.set_xlim(-0.7, n_samples - 0.3)
ax_a.set_ylim(0, yv.max() * 1.30)
ax_a.set_xlabel("Sample (GSM35166xx), grouped by clinical category", labelpad=25)
ax_a.set_ylabel("Cells retained after QC (n)")
ax_a.set_title("Cell yield spans 1,036-4,923 across 17 samples", loc="left",
               pad=6)
# direct labels: coloured bracket under the axis for each category
bt = blended_transform_factory(ax_a.transData, ax_a.transAxes)
for cat in CAT_ORDER:
    idx = np.flatnonzero((samples["category"] == cat).to_numpy())
    lo, hi = idx.min() - 0.38, idx.max() + 0.38
    ax_a.plot([lo, hi], [-0.135, -0.135], transform=bt, color=CAT_COLOR[cat],
              linewidth=1.8, clip_on=False, solid_capstyle="butt")
    ax_a.text((lo + hi) / 2.0, -0.175, CAT_SHORT[cat], transform=bt,
              ha="center", va="top", fontsize=6.2, color=CAT_COLOR[cat],
              fontweight="bold", linespacing=1.05)
ax_a.legend([plt.Rectangle((0, 0), 1, 1, color=CAT_COLOR[c]) for c in CAT_ORDER],
            ["%s (n=%d)" % (c, int((samples["category"] == c).sum()))
             for c in CAT_ORDER],
            frameon=False, loc="upper right", handlelength=0.85,
            handleheight=0.85, borderpad=0.1, labelspacing=0.22,
            borderaxespad=0.1)
caveats[ax_a] = ("Post-QC cells; median %.0f%% (min %.0f%%) of barcode-filtered "
                 "cells survive QC" % (median_retention, min_retention))

# ---------------------------------------------------------------------- (b)
xb = np.arange(4)
ax_b.bar(xb, b_counts, width=0.72,
         color=[UNUSABLE_COLOR] + [COHORT_PAS_COLOR] * 3, edgecolor="none")
style_axis(ax_b)
ax_b.set_xticks(xb)
ax_b.set_xticklabels(b_labels)
ax_b.set_xlabel("PAS per gene (n)")
ax_b.set_ylabel("Genes (n)")
ax_b.set_title("Only %.0f%% of 17/17-shared genes have >1 PAS"
               % (100 * multi_frac), loc="left", pad=6)
ax_b.set_ylim(0, max(b_counts) * 1.42)
for xi, ci in zip(xb, b_counts):
    ax_b.text(xi, ci + max(b_counts) * 0.02, "{:,}".format(ci), ha="center",
              va="bottom", fontsize=6.6, color=INK)
ax_b.text(0, max(b_counts) * 1.14, "Single-PAS\n%s genes (%.0f%%)"
          % ("{:,}".format(b_counts[0]), 100 * (1 - multi_frac)),
          ha="center", va="bottom", fontsize=6.2, color=UNUSABLE_COLOR,
          fontweight="bold", linespacing=1.1)
ax_b.text(2.4, max(b_counts) * 0.52, "Multi-PAS (APA-testable)\n%s genes (%.0f%%)"
          % ("{:,}".format(n_multi), 100 * multi_frac), ha="center", va="bottom",
          fontsize=6.2, color=COHORT_PAS_COLOR, fontweight="bold",
          linespacing=1.1)
ax_b.legend([plt.Rectangle((0, 0), 1, 1, color=UNUSABLE_COLOR),
             plt.Rectangle((0, 0), 1, 1, color=COHORT_PAS_COLOR)],
            ["Single-PAS gene", "Multi-PAS gene"], frameon=False,
            loc="upper right", handlelength=0.85, handleheight=0.85,
            borderpad=0.1, labelspacing=0.22, borderaxespad=0.1)
caveats[ax_b] = ("Set = PAS found in all %d/%d samples (%.1f%% of %s)\n"
                 "Median 3'UTR %s bp single-PAS vs %s bp multi-PAS"
                 % (n_samples_detect, n_samples_detect, 100 * cohort_frac_of_any,
                    "{:,}".format(n_ids_any),
                    "{:,.0f}".format(utr_single_med),
                    "{:,.0f}".format(utr_multi_med)))

# ---------------------------------------------------------------------- (c)
tiers = ["TIER_1", "TIER_2"]
tier_names = ["TIER 1\nin annotated 3'UTR", "TIER 2\n3'UTR extension",
              "Not retained\nTIER 3/intergenic"]
# Each series is a share of ITS OWN set. The third category is what makes the
# comparison honest: the cohort set has 0% there by construction, the called-peak
# universe has 38.9%, which is exactly the filter that produces the 99.4% purity.
cohort_pct = [100.0 * cohort_tier.get(t, 0) / n_pas for t in tiers] + [0.0]
all_pct = ([100.0 * ann_tier_all.get(t, 0) / n_called for t in tiers]
           + [100.0 * n_dropped / n_called])
xc = np.arange(3)
bw = 0.36
ax_c.bar(xc - bw / 2 - 0.015, cohort_pct, width=bw, color=COHORT_PAS_COLOR,
         edgecolor="none")
ax_c.bar(xc + bw / 2 + 0.015, all_pct, width=bw, color=ALL_PAS_COLOR,
         edgecolor="none")
style_axis(ax_c)
ax_c.set_xticks(xc)
ax_c.set_xticklabels(tier_names, fontsize=6.2)
ax_c.set_xlim(-0.62, 2.62)
ax_c.set_ylabel("Share of own PAS set (%)")
ax_c.set_xlabel("Annotation confidence tier")
ax_c.set_title("TIER 1 purity is produced by an upstream filter", loc="left",
               pad=6)
ax_c.set_ylim(0, 132)
ax_c.set_yticks([0, 25, 50, 75, 100])
for xi, v in zip(xc - bw / 2 - 0.015, cohort_pct):
    ax_c.text(xi, v + 2.5, "%.1f%%" % v, ha="center", va="bottom", fontsize=6.6,
              color=COHORT_PAS_COLOR, fontweight="bold")
for xi, v in zip(xc + bw / 2 + 0.015, all_pct):
    ax_c.text(xi, v + 2.5, "%.1f%%" % v, ha="center", va="bottom", fontsize=6.6,
              color=ALL_PAS_COLOR, fontweight="bold")
# Direct labels, each sitting over its own bar in the TIER 1 group.
ax_c.text(-bw / 2 - 0.015, cohort_pct[0] + 8.5,
          "PAS in all\n%d/%d samples\nn=%s" % (n_samples_detect, n_samples_detect,
                                               "{:,}".format(n_pas)),
          fontsize=6.2, color=COHORT_PAS_COLOR, fontweight="bold", ha="center",
          va="bottom", linespacing=1.15)
ax_c.text(bw / 2 + 0.015, all_pct[0] + 8.5,
          "All called\npeaks\nn=%s" % "{:,}".format(n_called), fontsize=6.2,
          color=ALL_PAS_COLOR, fontweight="bold", ha="center", va="bottom",
          linespacing=1.15)
caveats[ax_c] = ("find_close keeps TIER 1+2 only, so %s of %s called peaks "
                 "(%.1f%%) were dropped\nThe cohort set's 0%% in the third "
                 "category is by construction, not a result"
                 % ("{:,}".format(n_dropped), "{:,}".format(n_called),
                    100.0 * n_dropped / n_called))

# ---------------------------------------------------------------------- (d)
ax_d.bar(w_edges[:-1], w_hist, width=np.diff(w_edges), align="edge",
         color=COHORT_PAS_COLOR, edgecolor=SURFACE, linewidth=0.35)
ax_d.set_xscale("log")
style_axis(ax_d)
ax_d.set_ylim(0, w_hist.max() * 1.28)
ax_d.axvline(w_med, color=INK, linewidth=1.0, linestyle=(0, (4, 2)))
ax_d.annotate("median %.0f bp\n(IQR %.0f-%.0f bp)"
              % (w_med, np.percentile(w, 25), np.percentile(w, 75)),
              xy=(w_med, w_hist.max() * 1.03),
              xytext=(w_med * 1.55, w_hist.max() * 1.24),
              fontsize=6.4, color=INK, fontweight="bold", ha="left", va="top",
              arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.7))
ax_d.set_xlabel("PAS peak width, end - start (bp, log scale)")
ax_d.set_ylabel("PAS (n)")
ax_d.set_title("PAS are broad peak intervals, not single bases", loc="left",
               pad=6)
ax_d.set_xticks([100, 300, 1000, 3000, 10000])
ax_d.set_xticklabels(["100", "300", "1,000", "3,000", "10,000"])
ax_d.minorticks_off()
caveats[ax_d] = "Width is the called peak interval; the cleavage site lies inside it"

fig.tight_layout(pad=1.0, h_pad=4.2, w_pad=2.4, rect=[0.0, 0.075, 1.0, 0.985])

# panel tags + caveat strips, placed under each panel's rendered extent
fig.canvas.draw()
rend = fig.canvas.get_renderer()
fh = fig.get_figheight() * fig.dpi
for ax, letter in zip([ax_a, ax_b, ax_c, ax_d], ["a", "b", "c", "d"]):
    p = ax.get_position()
    fig.text(max(p.x0 - 0.058, 0.004), p.y1 + 0.012, letter, fontsize=10,
             fontweight="bold", va="bottom", ha="left", color=INK)
    bottom = ax.get_tightbbox(rend).y0 / fh
    fig.text(p.x0, min(bottom, p.y0) - 0.012, caveats[ax], fontsize=5.9,
             color=MUTED, ha="left", va="top")

fig.savefig(PNG, dpi=300, facecolor=SURFACE)
save_manuscript(fig, "cohort_qc", facecolor=SURFACE)
plt.close(fig)

# ================================================================ audit table
recs = []
for _, r in samples.iterrows():
    recs.append(dict(panel="a", metric="cells_retained_after_qc",
                     category=r["category"], label=r["sample"],
                     value=r["cells_retained"],
                     denominator=r["cells_after_cb_filter"],
                     value_pct=round(r["retention_pct"], 2)))
for lab, cnt in zip(b_labels, b_counts):
    recs.append(dict(panel="b", metric="genes_with_n_pas",
                     category=("single_PAS" if lab == "1" else "multi_PAS"),
                     label="%s PAS" % lab, value=cnt, denominator=n_genes,
                     value_pct=round(100.0 * cnt / n_genes, 2)))
tier_keys = tiers + ["NOT_RETAINED_TIER_3_OR_INTERGENIC"]
cohort_counts_c = [int(cohort_tier.get(t, 0)) for t in tiers] + [0]
all_counts_c = [int(ann_tier_all.get(t, 0)) for t in tiers] + [int(n_dropped)]
for t, v, p in zip(tier_keys, cohort_counts_c, cohort_pct):
    recs.append(dict(panel="c", metric="pas_by_tier",
                     category="pas_detected_in_all_17_samples",
                     label=t, value=v, denominator=n_pas, value_pct=round(p, 3)))
for t, v, p in zip(tier_keys, all_counts_c, all_pct):
    recs.append(dict(panel="c", metric="pas_by_tier", category="all_called_peaks",
                     label=t, value=v, denominator=n_called,
                     value_pct=round(p, 3)))
for lo, hi, cnt in zip(w_edges[:-1], w_edges[1:], w_hist):
    recs.append(dict(panel="d", metric="pas_width_histogram",
                     category="cohort_pas_set",
                     label="%.1f-%.1f bp" % (lo, hi), value=int(cnt),
                     denominator=n_pas, value_pct=round(100.0 * cnt / n_pas, 3)))

summary = [
    ("total_cohort_pas", "pas_detected_in_all_17_samples", "pasbed.bed", n_pas,
     n_pas, 100.0),
    ("called_peaks_total", "all", "posbed.bed + negbed.bed", n_called, n_called,
     100.0),
    ("gene_assigned_tier1_2_pas", "all", "annotatedpas.bed", n_ann, n_called,
     round(100.0 * n_ann / n_called, 2)),
    ("called_peaks_dropped_tier3_intergenic", "all", "not in annotatedpas.bed",
     n_dropped, n_called, round(100.0 * n_dropped / n_called, 2)),
    ("pas_detected_in_ge1_sample", "all", "union over 17 samples", n_ids_any,
     n_ids_any, 100.0),
    ("cohort_set_share_of_union", "pas_detected_in_all_17_samples",
     "16,500 / union", n_pas, n_ids_any, round(100.0 * cohort_frac_of_any, 2)),
    ("cohort_pas_unannotated", "all", "no gene_id", n_unannotated, n_pas,
     round(100.0 * n_unannotated / n_pas, 3)),
    ("genes_with_ge1_pas", "all", "genes", n_genes, n_genes, 100.0),
    ("multi_pas_genes", "multi_PAS", "genes with >=2 PAS", n_multi, n_genes,
     round(100.0 * multi_frac, 2)),
    ("max_pas_per_gene", "all", "genes", max_pas_per_gene, n_genes, np.nan),
    ("samples", "all", "GSM dirs", n_samples, n_samples, 100.0),
    ("cells_retained_total", "all", "cohort", cells_total,
     int(samples["cells_after_cb_filter"].sum()),
     round(100.0 * cells_total / samples["cells_after_cb_filter"].sum(), 2)),
    ("cells_retained_min", "all", "per sample", cells_min, np.nan, np.nan),
    ("cells_retained_max", "all", "per sample", cells_max, np.nan, np.nan),
    ("cells_retained_median", "all", "per sample",
     float(samples["cells_retained"].median()), np.nan, np.nan),
    ("pas_width_median_bp", "cohort_pas_set", "bp", w_med, np.nan, np.nan),
    ("pas_width_p25_bp", "cohort_pas_set", "bp", float(np.percentile(w, 25)),
     np.nan, np.nan),
    ("pas_width_p75_bp", "cohort_pas_set", "bp", float(np.percentile(w, 75)),
     np.nan, np.nan),
    ("pas_width_min_bp", "cohort_pas_set", "bp", float(w.min()), np.nan, np.nan),
    ("pas_width_max_bp", "cohort_pas_set", "bp", float(w.max()), np.nan, np.nan),
    ("utr_median_bp_single_pas_genes", "single_PAS", "bp", utr_single_med,
     np.nan, np.nan),
    ("utr_median_bp_multi_pas_genes", "multi_PAS", "bp", utr_multi_med,
     np.nan, np.nan),
]
for m, c, l, v, d, pc in summary:
    recs.append(dict(panel="summary", metric=m, category=c, label=l, value=v,
                     denominator=d, value_pct=pc))

pd.DataFrame(recs)[["panel", "metric", "category", "label", "value",
                    "denominator", "value_pct"]].to_csv(TSV, sep="\t",
                                                        index=False)

print("wrote", PNG)
print("wrote", TSV)
print("samples=%d cells %d-%d (total %d), PAS=%d, genes=%d, multi-PAS=%d (%.1f%%)"
      % (n_samples, cells_min, cells_max, cells_total, n_pas, n_genes, n_multi,
         100 * multi_frac))
print("tier cohort:", dict(cohort_tier), " tier all:", dict(ann_tier_all))
print("called peaks %d -> gene-assigned TIER1/2 %d (%.2f%%), dropped %d (%.2f%%)"
      % (n_called, n_ann, 100.0 * n_ann / n_called, n_dropped,
         100.0 * n_dropped / n_called))
print("cohort set is the %d/%d intersection: %d of %d PAS seen in >=1 sample "
      "(%.2f%%)" % (n_samples_detect, n_samples_detect, n_pas, n_ids_any,
                    100 * cohort_frac_of_any))
print("width median %.0f IQR %.0f-%.0f max %.0f"
      % (w_med, np.percentile(w, 25), np.percentile(w, 75), w.max()))
print("unannotated cohort PAS:", n_unannotated, "max PAS/gene:", max_pas_per_gene)
print("retention median %.1f%% min %.1f%%" % (median_retention, min_retention))
