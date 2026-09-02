#!/usr/bin/env python3
"""
figS7_novelfunnel.py -- the pre-registered "trusted de novo PAS" definition (13 §3) on the
final PBMC precision default (IP arm) against Kinnex long-read truth.

VERSION SWITCH -- env var TRUSTED_NOVEL_VERSION picks the source run:
  v2  (default; THE PAPER NUMBERS)  code 9dfdefb (#96 IP-filter strand fix + #97),
      results/reliability/trusted_novel_final_v2_pbmc/ -- REPORT_PROVISIONAL.md
      STATUS: VERIFIED (adversarial verifier, verdict FIXED); manuscript/19 §4.
  v1  (TRUSTED_NOVEL_VERSION=v1; kept selectable for the record)  code 4efeb125,
      results/reliability/trusted_novel_final_pbmc/, manuscript/16 (verified FIXED).

Verdict: the pre-registered definition reaches 0.485 (v2; v1 0.521) within 25 bp of a Kinnex
x3p poly(A)-verified 3' end at >=5 UMI; the pre-registered target was >=0.70. NEGATIVE RESULT.
Every number on the figure is read from files under the selected run directory (nothing typed in):

  call_primary/funnel.tsv                         panel a  counts per stage
  validate_incl_genebodies/validation.tsv         panels a,b  concordance + Wilson CI + null per stage/set/threshold
  validate_incl_genebodies/validation_null_seeds.tsv  panel b  null band (10 gene-body shuffles)
  feature_breakdown.tsv                           panel c  3'UTR / other exon / intron / intergenic
  distance_crosscheck.tsv                         panel c  fraction within 25 bp of a Kinnex IP-decoy terminus
  stratified_concordance_trusted_novel.tsv        panel d  POST-HOC strata (feature class, clip-molecule support)
  verifier_crosscheck/RESULTS.txt                 footer   robustness (GEM-X truth, pooled truth, 50/100 bp)

Panel d is exploratory: the strata were chosen after seeing the data and are NOT a definition.
Any revised definition must be pre-registered and validated on data not used here.

Outputs
-------
manuscript/figures/figS7_novelfunnel.{png,pdf}        600 dpi / fonttype 42 (no Type 3);
    all on-figure prose (suptitle, subtitle, caveat block) moved to the caption
    sidecar's '## Legend' in the 2026-09-02 submission pass
results/figures/manuscript/figS7_novelfunnel.png
results/figures/manuscript/figS7_novelfunnel.tsv                  panel a (funnel + concordance per stage)
results/figures/manuscript/figS7_novelfunnel_concordance.tsv      panel b (every point, CI, null)
results/figures/manuscript/figS7_novelfunnel_composition.tsv      panel c (feature classes + decoy proximity)
results/figures/manuscript/figS7_novelfunnel_strata.tsv           panel d (post-hoc strata)
results/figures/manuscript/figS7_novelfunnel_reference_lines.tsv  target + robustness values + sources
"""
import os, re

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch
import numpy as np
import pandas as pd

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

WD = "/mnt/ssd1/Projects/PeakATail_wd"
FIGDIR = f"{WD}/manuscript/figures"
TSVDIR = f"{WD}/results/figures/manuscript"
os.makedirs(FIGDIR, exist_ok=True); os.makedirs(TSVDIR, exist_ok=True)
NAME = "figS7_novelfunnel"

# ---------------------------------------------------------------------------
# Which caller run feeds the funnel: v2 (default, the paper numbers) or v1 (record).
# Only the source directory, the code commit and the expected caller arm differ;
# the definition, target, window, truth sets and null are identical in both.
# ---------------------------------------------------------------------------
VERSION = os.environ.get("TRUSTED_NOVEL_VERSION", "v2").lower()
assert VERSION in ("v1", "v2"), f"TRUSTED_NOVEL_VERSION must be v1 or v2, got {VERSION!r}"
VERSIONS = {
    "v2": dict(src="results/reliability/trusted_novel_final_v2_pbmc", code="9dfdefb",
               arm="peakatail_clipseeded_final_v2_ipfilt",
               writeup="manuscript/19_final_gate_v2.md §4 (+ 16 for the v1 record)",
               status="REPORT_PROVISIONAL.md STATUS: VERIFIED (adversarial verifier, verdict FIXED)",
               prev=dict(label="v1 run (pre-fix caller, code 4efeb125)", frac=0.521, n_tn=6628, decoy=0.245)),
    "v1": dict(src="results/reliability/trusted_novel_final_pbmc", code="4efeb125",
               arm="peakatail_clipseeded_final_ipfilt",
               writeup="manuscript/16_trusted_novel_kinnex.md (verified FIXED)",
               status="REPORT_PROVISIONAL.md §1-9 (verifier-reproduced)",
               prev=None),
}[VERSION]
SRC = f"{WD}/{VERSIONS['src']}"
CODE = VERSIONS["code"]      # frozen snapshot of the input run; checked against run_manifest below
TARGET = 0.70                # pre-registered target (13 §3), committed before any number existed
WINDOW_BP = 25
Z95 = 1.959964

# Okabe-Ito (task-mandated; six-checks validated 2026-08-21 in this order: worst adjacent CVD dE 11.4 protan /
# 8.5 tritan, normal-vision floor 16.4 -- PASS; contrast WARN for sky/orange/pink -> every series carries a
# marker + line style + legend entry as secondary encoding).  Colour follows the SET, fixed across panels.
BLUE, SKY, GREEN, ORANGE, PINK, VERM, GREY = "#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#999999"
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"
# feature classes: an ordinal single-hue ramp (validated --ordinal: light end 2.20:1, monotone L, hue spread 7 deg)
FEAT_COLS = {"3utr": "#1B2429", "other_exon": "#43525A", "intronic": "#6E8089", "intergenic": "#9FAFB6"}
FEAT_LABEL = {"3utr": "3'UTR", "other_exon": "other exon", "intronic": "intron", "intergenic": "intergenic"}

SETS = {  # validation.tsv query id -> plotting meta (label n filled in from the file)
    "trusted_novel":            dict(label="trusted-novel (pre-registered)", color=BLUE,   ls="-",  marker="o", lw=1.9, mfc=BLUE,   z=8),
    "trusted_novel_strong":     dict(label="trusted-novel, strong hexamer",  color=SKY,    ls=(0, (4, 2)), marker="s", lw=1.2, mfc="none", z=7),
    "CAL_atlas_known_hexpass":  dict(label="atlas-known, hexamer-pass (calibration)", color=GREEN, ls="-", marker="D", lw=1.2, mfc=GREEN, z=6),
    "CAL_default_output_all":   dict(label="full precision default", color=ORANGE, ls="-", marker="^", lw=1.2, mfc=ORANGE, z=5),
    "CAL_atlas_novel_hexfail":  dict(label="atlas-novel, hexamer-FAIL (diagnostic)", color=PINK, ls=(0, (1.2, 1.2)), marker="v", lw=1.2, mfc="none", z=6),
    "SENS_min1_trusted_novel":  dict(label=">=1-molecule trusted-novel (sensitivity)", color=VERM, ls=(0, (4, 1.5, 1, 1.5)), marker="x", lw=1.2, mfc=VERM, z=5),
}
TRUTHS = {"t5": 5, "t20": 20, "t100": 100, "t500": 500}


def wilson(k, n, z=Z95):
    k, n = float(k), float(n)
    if n == 0:
        return np.nan, np.nan
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def fmt_ci(p, lo, hi):
    return f"{p:.3f} [{lo:.3f}-{hi:.3f}]"

# ----------------------------------------------------------------------------
# 1. read sources
# ----------------------------------------------------------------------------
funnel = pd.read_csv(f"{SRC}/call_primary/funnel.tsv", sep="\t").set_index("stage")
val = pd.read_csv(f"{SRC}/validate_incl_genebodies/validation.tsv", sep="\t")
null_seeds = pd.read_csv(f"{SRC}/validate_incl_genebodies/validation_null_seeds.tsv", sep="\t")
dist = pd.read_csv(f"{SRC}/distance_crosscheck.tsv", sep="\t")
feat = pd.read_csv(f"{SRC}/feature_breakdown.tsv", sep="\t").set_index("set")
strat = pd.read_csv(f"{SRC}/stratified_concordance_trusted_novel.tsv", sep="\t", header=None,
                    names=["stratum", "level", "n", "n_hit_t5", "frac_t5", "n_hit_t20", "frac_t20"], dtype=str)
strat = strat[strat["stratum"] != "stratum"].copy()      # the header line sits mid-file (sorted output)
for c in ("n", "n_hit_t5", "n_hit_t20"):
    strat[c] = strat[c].astype(int)
for c in ("frac_t5", "frac_t20"):
    strat[c] = strat[c].astype(float)
verifier_txt = open(f"{SRC}/verifier_crosscheck/RESULTS.txt").read()
report_txt = open(f"{SRC}/REPORT_PROVISIONAL.md").read()
manifest_txt = open(f"{SRC}/call_primary/run_manifest.json").read()
assert f"{VERSIONS['arm']}/pas_PRESPEC_precision_default.bed" in manifest_txt, \
    f"input is not the {VERSION} pre-registered precision default ({VERSIONS['arm']})"
assert re.search(rf"commit={CODE}[0-9a-f]*", report_txt), f"REPORT_PROVISIONAL.md does not record commit {CODE} for the input run"
assert (val["window_bp"] == WINDOW_BP).all() and (val["strand_aware"] == 1).all()
assert (val["null_seeds"] == 10).all()

V = val.set_index(["query", "truth"])
def vrow(q, t):
    return V.loc[(q, t)]

# sanity: Wilson implementation reproduces the tool's CI columns
r = vrow("trusted_novel", "t5")
lo, hi = wilson(r["n_hit"], r["n_query"])
assert abs(lo - r["wilson95_lo"]) < 2e-4 and abs(hi - r["wilson95_hi"]) < 2e-4, (lo, hi, r["wilson95_lo"], r["wilson95_hi"])

m = (re.search(r"\(([\d,]+) distinct genes;", report_txt)                        # v1 prose
     or re.search(r"distinct genes \|\s*trusted-novel ([\d,]+)", report_txt))     # v2 verifier table
N_GENES = m.group(1) if m else None
assert N_GENES, "could not read the distinct-gene count from REPORT_PROVISIONAL.md"

def verifier_frac(pattern):
    """parse 'label: k/n=frac' lines from the verifier's RESULTS.txt"""
    m = re.search(pattern + r":\s*(\d+)/(\d+)=([\d.]+)", verifier_txt)
    if not m:
        return None
    return dict(k=int(m.group(1)), n=int(m.group(2)), frac=float(m.group(3)))
rob = {
    "x3p_t5_25bp": verifier_frac(r"x3p t5 trusted_novel"),
    "gemx_t5_25bp": verifier_frac(r"GEM-X t5 trusted_novel"),
    "pooled_t5_25bp": verifier_frac(r"x3p\+GEM-X union t5 trusted_novel window 25"),
    "pooled_t5_50bp": verifier_frac(r"x3p\+GEM-X union t5 trusted_novel window 50"),
    "pooled_t5_100bp": verifier_frac(r"x3p\+GEM-X union t5 trusted_novel window 100"),
}
assert rob["x3p_t5_25bp"]["k"] == int(r["n_hit"]) and rob["x3p_t5_25bp"]["n"] == int(r["n_query"])

# ----------------------------------------------------------------------------
# 2. panel a data: funnel + concordance per stage
# ----------------------------------------------------------------------------
STAGES = [  # (funnel.tsv stage, validation.tsv query, label, colour)
    ("on_listed_contigs", "01_on_listed_contigs", "input: precision default, IP arm\n(sites on primary contigs)", GREY),
    ("not_internal_priming", "05_not_internal_priming", "clip-supported, >=2 molecules,\nnot IP (no-op: caller enforces)", GREY),
    ("hexamer_any12_-40..-5", "06_hexamer_any12_-40..-5", "canonical hexamer (any of 12)\nat -40..-5", GREY),
    ("trusted_novel", "08_trusted_novel", ">=100 bp from any PolyASite 2.0\nsite = TRUSTED-NOVEL", BLUE),
    ("trusted_novel_strong", "09_trusted_novel_strong", "AATAAA / ATTAAA only\n= strong subset", SKY),
]
a_rows = []
prev = None
for st, q, label, col in STAGES:
    n = int(funnel.loc[st, "n_pass"]); v = vrow(q, "t5")
    assert int(v["n_query"]) == n, (st, n, v["n_query"])
    d = dict(stage=st, validation_query=q, label=label.replace("\n", " "), n=n,
             frac_of_input=float(funnel.loc[st, "frac_of_input"]),
             n_hit_t5_25bp=int(v["n_hit"]), frac_t5=float(v["frac_hit"]),
             wilson95_lo=float(v["wilson95_lo"]), wilson95_hi=float(v["wilson95_hi"]),
             delta_frac_vs_previous=(float(v["frac_hit"]) - prev) if prev is not None else np.nan,
             null_mean_t5=float(v["null_mean"]), enrichment_t5=float(v["enrichment"]), color=col)
    a_rows.append(d); prev = float(v["frac_hit"])
A = pd.DataFrame(a_rows)
A.drop(columns="color").to_csv(f"{TSVDIR}/{NAME}.tsv", sep="\t", index=False)
N_INPUT = int(A.loc[0, "n"]); N_TN = int(A.loc[3, "n"]); N_STRONG = int(A.loc[4, "n"])
TN = vrow("trusted_novel", "t5"); FRAC_TN = float(TN["frac_hit"])

# ----------------------------------------------------------------------------
# 3. panel b data: concordance vs truth stringency, all sets
# ----------------------------------------------------------------------------
b_rows = []
for q, meta in SETS.items():
    for t, umi in TRUTHS.items():
        v = vrow(q, t)
        ns = null_seeds[(null_seeds["query"] == q) & (null_seeds["truth"] == t)]
        b_rows.append(dict(set=q, label=meta["label"], truth=t, umi_threshold=umi, n_query=int(v["n_query"]),
                           n_truth=int(v["n_truth"]), n_hit=int(v["n_hit"]), frac_hit=float(v["frac_hit"]),
                           wilson95_lo=float(v["wilson95_lo"]), wilson95_hi=float(v["wilson95_hi"]),
                           null_mean=float(v["null_mean"]), null_sd=float(v["null_sd"]),
                           null_min=float(ns["frac_hit"].min()), null_max=float(ns["frac_hit"].max()),
                           n_null_seeds=int(len(ns)), enrichment=float(v["enrichment"]),
                           empirical_p=float(v["empirical_p"]), meets_target_0_70=int(v["meets_target_0.70"])))
B = pd.DataFrame(b_rows)
B.to_csv(f"{TSVDIR}/{NAME}_concordance.tsv", sep="\t", index=False)
for q, meta in SETS.items():
    meta["n"] = int(B[B.set == q]["n_query"].iloc[0])

# ----------------------------------------------------------------------------
# 4. panel c data: composition + decoy proximity
# ----------------------------------------------------------------------------
COMP_SETS = [("trusted_novel", "trusted-novel", BLUE), ("CAL_atlas_known_hexpass", "atlas-known, hexamer-pass", GREEN)]
DECOY_SETS = [("trusted_novel", "trusted-novel", BLUE), ("trusted_novel_strong", "trusted-novel, strong", SKY),
              ("default_output_all", "full precision default", ORANGE), ("CAL_atlas_known_hexpass", "atlas-known, hexamer-pass", GREEN)]
c_rows = []
for s, lab, col in COMP_SETS:
    f = feat.loc[s]
    tot = 0
    for cls, ncol in (("3utr", "n_3utr"), ("other_exon", "n_other_exon"), ("intronic", "n_intronic"), ("intergenic", "n_intergenic_sense")):
        tot += int(f[ncol])
        c_rows.append(dict(panel="c_feature_class", set=s, label=lab, n=int(f["n"]), item=cls, count=int(f[ncol]),
                           frac=int(f[ncol]) / int(f["n"]), wilson95_lo=np.nan, wilson95_hi=np.nan))
    assert tot == int(f["n"]), (s, tot, f["n"])
D = dist.set_index(["query", "ref"])
for s, lab, col in DECOY_SETS:
    d = D.loc[(s, "decoy")]
    lo, hi = wilson(d["le25"], d["n"])
    c_rows.append(dict(panel="c_decoy_within25bp", set=s, label=lab, n=int(d["n"]), item="within_25bp_of_kinnex_ip_decoy",
                       count=int(d["le25"]), frac=float(d["frac_le25"]), wilson95_lo=lo, wilson95_hi=hi))
C = pd.DataFrame(c_rows)
C.to_csv(f"{TSVDIR}/{NAME}_composition.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 5. panel d data: POST-HOC strata of the trusted-novel set (exploratory only)
# ----------------------------------------------------------------------------
S = strat.set_index(["stratum", "level"])
def stratum(kind, level, file_levels):
    """sum the file's counts over one or more file levels (>=5 molecules = 5-9 + 10+)."""
    n = sum(int(S.loc[(kind, l), "n"]) for l in file_levels)
    k5 = sum(int(S.loc[(kind, l), "n_hit_t5"]) for l in file_levels)
    k20 = sum(int(S.loc[(kind, l), "n_hit_t20"]) for l in file_levels)
    lo5, hi5 = wilson(k5, n); lo20, hi20 = wilson(k20, n)
    return dict(group=kind, level=level, file_levels="+".join(file_levels), n=n, frac_of_trusted_novel=n / N_TN,
                n_hit_t5=k5, frac_t5=k5 / n, t5_wilson95_lo=lo5, t5_wilson95_hi=hi5,
                n_hit_t20=k20, frac_t20=k20 / n, t20_wilson95_lo=lo20, t20_wilson95_hi=hi20)
d_rows = [
    stratum("feature", "3'UTR", ["3utr"]), stratum("feature", "other exon", ["other_exon"]),
    stratum("feature", "intron", ["intronic"]), stratum("feature", "intergenic", ["intergenic"]),
    stratum("support", "2", ["2"]), stratum("support", "3-4", ["3-4"]), stratum("support", ">=5", ["5-9", "10+"]),
]
for lvl in ("5-9", "10+"):   # components of the pooled >=5 bin, recorded but not plotted
    d_rows.append({**stratum("support", lvl, [lvl]), "group": "support_component_not_plotted"})
Dd = pd.DataFrame(d_rows)
assert Dd[Dd.group == "feature"]["n"].sum() == N_TN and Dd[Dd.group == "support"]["n"].sum() == N_TN
Dd.insert(0, "note", "POST HOC / exploratory: strata chosen after seeing the data; not a definition; not validated on held-out data")
Dd.to_csv(f"{TSVDIR}/{NAME}_strata.tsv", sep="\t", index=False)

# reference lines / robustness values and where each comes from
ref_rows = [
    dict(item="pre_registered_target", value=TARGET, source="manuscript/13_reliability_positioning.md §3 (committed before any number)"),
    dict(item="trusted_novel_x3p_t5_25bp", value=rob["x3p_t5_25bp"]["frac"], source="verifier_crosscheck/RESULTS.txt; = validation.tsv trusted_novel/t5"),
    dict(item="trusted_novel_gemx_t5_25bp", value=rob["gemx_t5_25bp"]["frac"], source="verifier_crosscheck/RESULTS.txt (GEM-X truth, other donor)"),
    dict(item="trusted_novel_pooled_x3p_gemx_t5_25bp", value=rob["pooled_t5_25bp"]["frac"], source="verifier_crosscheck/RESULTS.txt (union of both donors' truth)"),
    dict(item="trusted_novel_pooled_x3p_gemx_t5_50bp", value=rob["pooled_t5_50bp"]["frac"], source="verifier_crosscheck/RESULTS.txt (NOT the pre-registered window)"),
    dict(item="trusted_novel_pooled_x3p_gemx_t5_100bp", value=rob["pooled_t5_100bp"]["frac"], source="verifier_crosscheck/RESULTS.txt (NOT the pre-registered window)"),
    dict(item="null_mean_trusted_novel_t5", value=float(TN["null_mean"]), source="validation.tsv (gene-body-shuffled null, 10 seeds)"),
    dict(item="enrichment_trusted_novel_t5", value=float(TN["enrichment"]), source="validation.tsv"),
    dict(item="empirical_p_floor", value=float(TN["empirical_p"]), source="validation.tsv (1/(10 seeds + 1))"),
    dict(item="n_genes_trusted_novel", value=(int(N_GENES.replace(",", "")) if N_GENES else np.nan), source="REPORT_PROVISIONAL.md §5 (gene_id in run/annotatedpas.bed; subtitle only)"),
    dict(item="n_sites_default_incl_offcontig", value=int(funnel.loc["input", "n_pass"]), source=f"call_primary/funnel.tsv 'input' (denominator of frac_of_input; {int(funnel.loc['input', 'n_pass']) - int(funnel.loc['on_listed_contigs', 'n_pass'])} off-contig sites excluded before the funnel)"),
]
pd.DataFrame(ref_rows).to_csv(f"{TSVDIR}/{NAME}_reference_lines.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 6. figure
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 7.4, "axes.labelsize": 7.2, "axes.titlesize": 8.0,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.0,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 5, "axes.linewidth": 0.7,
    "font.family": "DejaVu Sans", "hatch.linewidth": 0.6,
})
# canvas: the pre-submission 7.5 x 8.7 in canvas carried a suptitle+subtitle
# band on top and a 12-line caveat block below; both now live in the caption
# sidecar ('## Legend') and the canvas height drops by the freed space.
# figsize narrowed 7.5 -> 7.1 in; because this script saves with
# bbox_inches='tight', the funnel's long y labels (left) and panel b's legend
# (right) overhang the declared canvas, so the DELIVERED image is ~8.9 in
# (225 mm) wide, down from 9.3 in pre-surgery.  Narrowing further crowds the
# 5.2-6.5 pt annotations (legibility outranks the 180 mm width target).
fig = plt.figure(figsize=(7.1, 7.3))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.05, 1.0], hspace=0.46, wspace=0.30,
                      left=0.058, right=0.985, top=0.965, bottom=0.105)
gs_a = gs[0, 0].subgridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.05)
ax_a1 = fig.add_subplot(gs_a[0, 0]); ax_a2 = fig.add_subplot(gs_a[0, 1], sharey=ax_a1)
ax_b = fig.add_subplot(gs[0, 1])
gs_c = gs[1, 0].subgridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.75)
ax_c1 = fig.add_subplot(gs_c[0, 0]); ax_c2 = fig.add_subplot(gs_c[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

def style(ax, axis="y"):
    ax.set_axisbelow(True); ax.grid(True, axis=axis, color=GRID, lw=0.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

# --- a: funnel (left) + concordance per stage (right) ------------------------
y = np.arange(len(A))
ax_a1.barh(y, A["n"], color=A["color"], height=0.62, zorder=3)
for i, rw in A.iterrows():
    ax_a1.text(rw["n"] + N_INPUT * 0.015, i, f"{rw['n']:,} ({rw['frac_of_input']:.1%})".replace("100.0%", "100%"), va="center", ha="left", fontsize=6.0, color=INK)
ax_a1.set_yticks(y); ax_a1.set_yticklabels([s[2] for s in STAGES], fontsize=5.9, linespacing=1.15)
ax_a1.invert_yaxis(); ax_a1.set_xlim(0, N_INPUT * 1.95)
ax_a1.set_xticks([0, 20000, 40000]); ax_a1.set_xticklabels(["0", "20k", "40k"])
N_DEFAULT = int(funnel.loc["input", "n_pass"])   # 44,413 incl. 19 off-contig sites; frac_of_input is relative to this
style(ax_a1, "x"); ax_a1.set_xlabel(f"sites retained (% of\nthe {N_DEFAULT:,}-site default)", fontsize=6.3)
ax_a1.set_title("a  Pre-registered funnel (13 §3) and per-stage concordance")

ax_a2.axvline(TARGET, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax_a2.text(TARGET + 0.01, len(A) - 0.45, "target 0.70", fontsize=5.8, color=INK, ha="left", va="center")
ax_a2.errorbar(A["frac_t5"], y, xerr=[A["frac_t5"] - A["wilson95_lo"], A["wilson95_hi"] - A["frac_t5"]],
               fmt="none", ecolor=INK, elinewidth=0.8, capsize=1.5, zorder=4)
ax_a2.scatter(A["frac_t5"], y, s=26, c=A["color"], edgecolor=SURFACE, linewidth=0.6, zorder=5)
for i, rw in A.iterrows():
    ax_a2.text(rw["frac_t5"] + 0.022, i, f"{rw['frac_t5']:.3f}", ha="left", va="center", fontsize=6.0, color=INK, zorder=6)
    if i > 0:
        d = rw["delta_frac_vs_previous"]
        big = d < -0.1
        ax_a2.text(0.415, i - 0.5, (f"{d:+.3f}" if abs(d) > 0 else "no change"), ha="left", va="center", fontsize=5.8,
                   color=(VERM if big else MUTED), fontweight=("bold" if big else "normal"), zorder=6)
ax_a2.text(0.415, -0.5, "change vs\nprevious stage", ha="left", va="center", fontsize=5.2, color=MUTED)
ax_a2.set_xlim(0.40, 1.0); ax_a2.set_xticks([0.5, 0.7, 0.9]); ax_a2.set_ylim(len(A) - 0.35, -0.85)
plt.setp(ax_a2.get_yticklabels(), visible=False); ax_a2.tick_params(axis="y", length=0)
style(ax_a2, "x"); ax_a2.set_xlabel("fraction within 25 bp of a\nKinnex 3' end (>=5 UMI),\nWilson 95% CI", fontsize=6.3)

# --- b: concordance vs truth stringency --------------------------------------
xs = np.array(list(TRUTHS.values()), dtype=float)
nb = B[B.set == "trusted_novel"].sort_values("umi_threshold")
ax_b.fill_between(xs, nb["null_min"].values, nb["null_max"].values, color=GREY, alpha=0.6, lw=0, zorder=1)   # true 10-seed range, not widened
ax_b.plot(xs, nb["null_mean"].values, color=GREY, lw=1.3, ls="-", zorder=2)
ax_b.axhline(TARGET, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax_b.text(560, TARGET + 0.012, "pre-registered target 0.70", fontsize=5.9, color=INK, ha="right", va="bottom")
handles = []
for q, meta in SETS.items():
    g = B[B.set == q].sort_values("umi_threshold")
    ax_b.errorbar(g["umi_threshold"], g["frac_hit"], yerr=[g["frac_hit"] - g["wilson95_lo"], g["wilson95_hi"] - g["frac_hit"]],
                  fmt="none", ecolor=meta["color"], elinewidth=0.7, capsize=1.4, zorder=meta["z"])
    ax_b.plot(g["umi_threshold"], g["frac_hit"], ls=meta["ls"], lw=meta["lw"], color=meta["color"], zorder=meta["z"])
    ax_b.plot(g["umi_threshold"], g["frac_hit"], ls="none", marker=meta["marker"], ms=4.6 if q == "trusted_novel" else 3.8,
              mfc=meta["mfc"], mec=meta["color"], mew=0.9, zorder=meta["z"] + 1)
    lab = f"{meta['label']}, n = {meta['n']:,}"
    if q == "trusted_novel":
        lab += f": {fmt_ci(FRAC_TN, TN['wilson95_lo'], TN['wilson95_hi'])} at >=5 UMI"
    handles.append(Line2D([], [], ls=meta["ls"], lw=meta["lw"], color=meta["color"], marker=meta["marker"], ms=3.8,
                          mfc=meta["mfc"], mec=meta["color"], mew=0.9, label=lab))
ax_b.set_xscale("log"); ax_b.set_xticks(xs); ax_b.set_xticklabels([f">={int(v)}" for v in xs]); ax_b.minorticks_off()
ax_b.set_xlim(4.2, 620); ax_b.set_ylim(0, 1.36); ax_b.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax_b.spines["left"].set_bounds(0, 1.0)
ax_b.set_xlabel("Kinnex truth stringency: UMI per 3' end (log)"); ax_b.set_ylabel("fraction within 25 bp, strand-matched (Wilson 95% CI)")
handles.append(Line2D([], [], color=GREY, lw=1.3, label=f"gene-body-shuffled null, 10 seeds (mean {float(nb['null_mean'].iloc[0]):.3f}, "
                      f"range {float(nb['null_min'].iloc[0]):.3f}-{float(nb['null_max'].iloc[0]):.3f} at >=5 UMI)"))
style(ax_b); ax_b.legend(handles=handles, loc="upper left", frameon=False, handlelength=2.4, labelspacing=0.3, borderaxespad=0.1, fontsize=5.6)
ax_b.set_title("b  Concordance vs truth stringency: target missed throughout")

# --- c1: feature-class composition --------------------------------------------
rows_c = C[C.panel == "c_feature_class"]
yc = np.arange(len(COMP_SETS))
for j, (s, lab, col) in enumerate(COMP_SETS):
    left = 0.0
    g = rows_c[rows_c.set == s].set_index("item")
    for cls in ("3utr", "other_exon", "intronic", "intergenic"):
        fr = float(g.loc[cls, "frac"])
        ax_c1.barh(j, fr, left=left, color=FEAT_COLS[cls], height=0.58, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        if fr >= 0.06:
            ax_c1.text(left + fr / 2, j, f"{fr:.1%}", ha="center", va="center", fontsize=5.9,
                       color=(SURFACE if cls in ("3utr", "other_exon", "intronic") else INK), zorder=4)
        left += fr
ax_c1.set_yticks(yc); ax_c1.set_yticklabels([f"{lab}\nn = {int(feat.loc[s, 'n']):,}" for s, lab, _ in COMP_SETS], fontsize=6.2, linespacing=1.15)
for t, (_, _, col) in zip(ax_c1.get_yticklabels(), COMP_SETS):
    t.set_color(col)
ax_c1.invert_yaxis(); ax_c1.set_xlim(0, 1); ax_c1.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax_c1.set_xticklabels(["0", "25", "50", "75", "100%"])
style(ax_c1, "x"); ax_c1.set_xlabel("share of sites (Ensembl 99 GTF, strand-matched, hierarchical)", labelpad=2)
ax_c1.legend(handles=[Patch(facecolor=FEAT_COLS[c], label=FEAT_LABEL[c]) for c in FEAT_COLS], loc="upper center",
             bbox_to_anchor=(0.5, -0.46), ncol=4, frameon=False, handlelength=1.2, columnspacing=1.2, handletextpad=0.5)
ax_c1.set_title("c  Feature class and internal-priming decoy proximity")

# --- c2: decoy proximity ---------------------------------------------------------
rows_d = C[C.panel == "c_decoy_within25bp"]
yd = np.arange(len(DECOY_SETS))
cols = [c for _, _, c in DECOY_SETS]
ax_c2.barh(yd, rows_d["frac"], color=cols, height=0.58, zorder=3)
ax_c2.errorbar(rows_d["frac"], yd, xerr=[rows_d["frac"] - rows_d["wilson95_lo"], rows_d["wilson95_hi"] - rows_d["frac"]],
               fmt="none", ecolor=INK, elinewidth=0.7, capsize=1.4, zorder=4)
for j, (_, rw) in enumerate(rows_d.iterrows()):
    ax_c2.text(rw["wilson95_hi"] + 0.006, j, f"{rw['frac']:.1%} ({rw['count']:,} / {rw['n']:,})", va="center", ha="left", fontsize=5.9, color=INK)
ax_c2.set_yticks(yd); ax_c2.set_yticklabels([lab for _, lab, _ in DECOY_SETS], fontsize=6.2)
ax_c2.invert_yaxis(); ax_c2.set_xlim(0, 0.40); ax_c2.set_xticks([0, 0.1, 0.2, 0.3, 0.4]); ax_c2.set_xticklabels(["0", "10", "20", "30", "40%"])
style(ax_c2, "x")
ax_c2.set_xlabel("within 25 bp (same strand) of a Kinnex x3p internal-priming decoy terminus\n(Kinnex: >=12 A in +1..+18; caller's IP rule: -10..+30, >=6 A or >=70% A)", labelpad=2, fontsize=6.0)

# --- d: POST-HOC strata -------------------------------------------------------------
plot_d = Dd[Dd.group.isin(["feature", "support"])].reset_index(drop=True)
xpos = []; x = 0.0
for i, rw in plot_d.iterrows():
    if i > 0 and rw["group"] != plot_d.loc[i - 1, "group"]:
        x += 0.9
    xpos.append(x); x += 1.0
xpos = np.array(xpos); w = 0.36
T20 = float(vrow("trusted_novel", "t20")["frac_hit"])
XR = xpos[-1] + 0.62
ax_d.axhline(TARGET, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax_d.axhline(FRAC_TN, color=BLUE, lw=0.7, ls="-", alpha=0.6, zorder=2)
ax_d.axhline(T20, color=BLUE, lw=0.7, ls=(0, (1, 1.5)), alpha=0.6, zorder=2)
LBOX = dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none")
ax_d.text(XR, TARGET, "target\n0.70", fontsize=5.3, color=INK, ha="left", va="center", linespacing=1.1, bbox=LBOX, zorder=6)
ax_d.text(XR, FRAC_TN, f"all, >=5 UMI\n{FRAC_TN:.3f}", fontsize=5.3, color=BLUE, ha="left", va="center", linespacing=1.1, bbox=LBOX, zorder=6)
ax_d.text(XR, T20, f"all, >=20 UMI\n{T20:.3f}", fontsize=5.3, color=BLUE, ha="left", va="center", linespacing=1.1, bbox=LBOX, zorder=6)
ax_d.bar(xpos - w / 2, plot_d["frac_t5"], width=w * 0.94, color=BLUE, zorder=3)
ax_d.bar(xpos + w / 2, plot_d["frac_t20"], width=w * 0.94, facecolor="none", edgecolor=BLUE, linewidth=1.0, zorder=3)
ax_d.errorbar(xpos - w / 2, plot_d["frac_t5"], yerr=[plot_d["frac_t5"] - plot_d["t5_wilson95_lo"], plot_d["t5_wilson95_hi"] - plot_d["frac_t5"]],
              fmt="none", ecolor=INK, elinewidth=0.7, capsize=1.3, zorder=4)
ax_d.errorbar(xpos + w / 2, plot_d["frac_t20"], yerr=[plot_d["frac_t20"] - plot_d["t20_wilson95_lo"], plot_d["t20_wilson95_hi"] - plot_d["frac_t20"]],
              fmt="none", ecolor=INK, elinewidth=0.7, capsize=1.3, zorder=4)
for xi, (_, rw) in zip(xpos, plot_d.iterrows()):
    ax_d.text(xi - w / 2, rw["t5_wilson95_hi"] + 0.012, f"{rw['frac_t5']:.2f}", ha="center", va="bottom", fontsize=5.6, color=INK)
    ax_d.text(xi + w / 2, rw["t20_wilson95_hi"] + 0.012, f"{rw['frac_t20']:.2f}", ha="center", va="bottom", fontsize=5.6, color=MUTED)
ax_d.set_xticks(xpos)
TICK_NAME = {"other exon": "other\nexon", "intergenic": "inter-\ngenic"}
ax_d.set_xticklabels([f"{TICK_NAME.get(rw['level'], rw['level'])}\n{rw['n']:,}\n({rw['frac_of_trusted_novel']:.0%})" for _, rw in plot_d.iterrows()],
                     fontsize=5.5, linespacing=1.12)
for grp, lab in (("feature", "by feature class"), ("support", "by clip-molecule support")):
    sel = xpos[plot_d.group.values == grp]
    ax_d.text(sel.mean(), -0.255, lab, transform=ax_d.get_xaxis_transform(), ha="center", va="top", fontsize=6.2, color=INK, fontweight="bold")
ax_d.set_xlim(xpos[0] - 0.6, xpos[-1] + 1.75); ax_d.set_ylim(0, 1.0)
ax_d.set_ylabel("fraction within 25 bp (Wilson 95% CI)")
style(ax_d)
ax_d.legend(handles=[Patch(facecolor=BLUE, label="Kinnex >=5 UMI"), Patch(facecolor="none", edgecolor=BLUE, lw=1.0, label="Kinnex >=20 UMI")],
            loc="upper left", frameon=False, ncol=1, handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)
ax_d.set_title("d  POST HOC, not a definition: trusted-novel sites by stratum", color=VERM)
ax_d.text(0.99, 0.985, "exploratory: strata chosen AFTER seeing the data;\n>=5 molecules = 5-9 and 10+ bins pooled", transform=ax_d.transAxes,
          ha="right", va="top", fontsize=5.2, color=VERM, linespacing=1.2)

# dashed frame around panel d so it cannot be read as part of the pre-registered result
fig.canvas.draw()
bb = ax_d.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
pad_x, pad_y = 0.008, 0.010
frame = FancyBboxPatch((bb.x0 - pad_x, bb.y0 - pad_y), bb.width + 2 * pad_x, bb.height + 2 * pad_y,
                       boxstyle="round,pad=0.002,rounding_size=0.006", transform=fig.transFigure, fill=False,
                       edgecolor=VERM, linewidth=0.9, linestyle=(0, (4, 2)), zorder=0, clip_on=False)
fig.add_artist(frame)

# --- title + caveats: moved off the image (2026-09-02 submission pass) --------
# The figure-level suptitle, the definition/result subtitle and the 12-line
# truth-caveat block now live in the caption sidecar's '## Legend' (substance
# verbatim; the sidecar describes the v2 default render).  A v1 record render
# (TRUSTED_NOVEL_VERSION=v1) therefore carries no on-figure prose either; its
# numbers remain in the run directory and manuscript/16.

fig.savefig(f"{TSVDIR}/{NAME}.png", dpi=600, bbox_inches="tight", pad_inches=0.06)
for ext in ("png", "pdf"):
    p = f"{FIGDIR}/{NAME}.{ext}"
    fig.savefig(p, bbox_inches="tight", pad_inches=0.06, **({"dpi": 600} if ext == "png" else {}))
    print("wrote", p)
# edge check: nothing may be clipped; the outer 8 px of the PNG must be blank
try:
    from PIL import Image
    im = np.asarray(Image.open(f"{FIGDIR}/{NAME}.png").convert("L"))
    edge = 8
    border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                             im[:, :edge].ravel(), im[:, -edge:].ravel()])
    ink = int((border < 250).sum())
    print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
    assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"
except ImportError:
    print("edge check skipped (no PIL)")

print(A[["stage", "n", "frac_of_input", "frac_t5", "delta_frac_vs_previous"]].to_string(index=False))
print(f"VERSION={VERSION} code={CODE} src={VERSIONS['src']}")
print(f"VERDICT: trusted-novel {FRAC_TN:.3f} [{TN['wilson95_lo']:.3f}-{TN['wilson95_hi']:.3f}] at >=5 UMI vs target {TARGET:.2f}: "
      f"{'MET' if FRAC_TN >= TARGET else 'NOT MET'}")

# ---------------------------------------------------------------------------
# sidecar caption -- script-written as of the 2026-09-02 rename pass (the
# figure convention requires the script, never a hand edit, to own the
# caption; the text below is the verified 2026-08-21 caption, stems
# renumbered, content unchanged).  It describes the v2 (default) render;
# it is only written when VERSION == 'v2' so a v1 record render cannot
# overwrite the manuscript caption with mismatched prose.
# ---------------------------------------------------------------------------
CAPTION_MD = r"""# Fig S7 — `figS7_novelfunnel` caption (generated by `scripts/manuscript_figures/figS7_novelfunnel.py`; renamed 2026-09-02, see FIGURE_MAP.tsv — demoted from 21's main Fig 6 slot, prominence kept in R6/abstract/Author Summary; submission pass 2026-09-02: the on-figure suptitle, subtitle and caveat block moved into the Legend below, image at 600 dpi)

## Legend

Figure S7 | Pre-registered trusted-novel definition vs Kinnex long-read truth: 49% concordance, target 70% not met (negative result). Input = the final PBMC 10k v3 precision default (IP arm, PeakATail code 9dfdefb; 46,524 sites on primary contigs). Definition (13 §3, committed before any number existed): clip-supported ∧ ≥2 molecules ∧ not internal priming ∧ canonical hexamer (any of 12) at −40..−5 ∧ ≥100 bp from any PolyASite 2.0 site → 7,259 trusted-novel sites (15.6%, 3,578 genes; 4,329 with AATAAA/ATTAAA). Truth = Kinnex x3p poly(A)-verified 3′ ends (different donor); hit = within 25 bp, same strand; Wilson 95% CIs; null = sites shuffled within gene bodies (10 seeds). Result: 0.485 [0.474–0.497] at ≥5 UMI (3,523 / 7,259); 0.239 / 0.066 / 0.012 at ≥20 / ≥100 / ≥500 UMI — the pre-registered ≥0.70 target is not met at any threshold. Robustness (verifier, not part of the pre-registered metric): GEM-X truth 0.515; x3p + GEM-X pooled 0.582 at 25 bp (0.669 at 50 bp, 0.714 at 100 bp). The sites are far from random (enrichment ~65× over the null at ≥5 UMI, empirical p at the 1/11 floor), but the target was an absolute fraction. Every plotted value: `results/figures/manuscript/figS7_novelfunnel*.tsv`.

**Panel a** — the pre-registered funnel as horizontal bars (sites retained, % of input) with, on the right, each
stage's concordance at ≥5 UMI: input 0.765 → clip / ≥2 molecules / not-IP 0.765 (no-op: the caller already
enforces them) → hexamer 0.811 (+0.046) → atlas-novel = trusted-novel 0.485 (−0.326) → strong hexamer 0.502
(+0.017). The atlas-novelty stage is the one that moves concordance; the dashed line is the 0.70 target.
**Panel b** — concordance vs truth stringency (x = Kinnex UMI threshold 5 / 20 / 100 / 500, log) for
trusted-novel (n = 7,259), its strong-hexamer subset (4,329), the atlas-known hexamer-pass complement as
calibration (28,453: 0.894 / 0.707 / 0.389 / 0.142), the full precision default (46,524: 0.765 / 0.558 /
0.281 / 0.098), the atlas-novel hexamer-FAIL set as diagnostic (6,432: 0.487 — the hexamer adds nothing among
atlas-novel sites) and the ≥1-molecule trusted-novel sensitivity arm (40,077: 0.289); grey line = null mean
(0.007 at ≥5 UMI) with the true 10-seed range (0.006–0.010) as a band — the band is not widened for visibility. **Panel c** — feature class (Ensembl 99, strand-matched, hierarchical 3′UTR > other
exon > intron > intergenic): trusted-novel 14.9% 3′UTR / 7.5% / 69.4% intron / 8.2%, atlas-known hexamer-pass
69.9% / 5.0% / 18.0% / 7.1%; below, the fraction of sites within 25 bp of a Kinnex x3p terminus flagged as
internal priming (decoy): trusted-novel 27.0% (1,959 / 7,259), strong 16.8%, full default 13.0%, atlas-known
hexamer-pass 7.6% (2,171 / 28,453). **Panel d (boxed, POST HOC)** — trusted-novel concordance at ≥5 and ≥20 UMI split by
feature class (3′UTR 0.77 / 0.55, n = 1,081; other exon 0.61 / 0.38, n = 545; intron 0.41 / 0.16, n = 5,037;
intergenic 0.47 / 0.20, n = 596) and by clip-molecule support (2: 0.40 / 0.14, n = 4,235; 3–4: 0.52 / 0.24,
n = 1,723; ≥5: 0.71 / 0.56, n = 1,301 = the file's 5–9 and 10+ bins pooled). Exploratory only.

**Must travel with it.**
- **This is a negative result and must be written as one.** The paper makes no "trusted novel" claim: the
  pre-registered definition reaches 0.485, not ≥0.70, and it is missed under every truth choice (x3p 0.485,
  GEM-X 0.515, pooled 0.582 at the pre-registered 25 bp). Do not quote the 50-bp (0.669) or 100-bp (0.714)
  pooled numbers as meeting the target — neither is the pre-registered window.
- **v1 → v2.** The pre-fix run (code 4efeb125, 16) gave 0.521 on 6,628 trusted-novel sites with 24.5% decoy
  proximity; the corrected minus-strand IP filter (#96) flags fewer raw peaks, so more A-rich-downstream sites
  survive into the default and the negative result stands and slightly strengthens (−0.036 concordance,
  +0.025 decoy proximity; 19 §4). The same upward decoy shift appears in every reference set (atlas-known
  hexamer-pass 6.6% → 7.6%, full default 11.2% → 13.0%), so it is a property of the new input, not of the
  trusted-novel subset. The two site sets overlap, so the v1-vs-v2 Δ is descriptive, not a test.
- **Panel d is post hoc.** The 3′UTR (0.77) and ≥5-molecule (0.71) strata were chosen after seeing the data;
  they are diagnostic input for a separately pre-registered v2 definition (3′UTR-restricted and/or ≥5
  molecules and/or a Kinnex-style IP rule, `--ip-rule kinnex`, Stage-1d) that must be validated on data not
  used here (e.g. a held-out truth set). Never call a de novo site "trusted" on the strength of
  panel d. The ≥5 bin pools the file's 5–9 (0.692 / 0.491, n = 672) and 10+ (0.728 / 0.625, n = 629) bins.
- Truth is site-level, not sample-level: both Kinnex sets come from donors other than the PBMC 10k v3 donor,
  so a site used here but absent (or <5 UMI) in the long-read donor counts as a miss. This caps the attainable
  fraction (atlas-known hexamer-pass: 0.894 x3p at ≥5 UMI; ≈0.39 at ≥100 UMI, where the
  target is unattainable for any set). Donor mismatch cannot explain the 0.41 gap because the atlas-known
  complement scores 0.894 under the same truth. The target is meaningfully testable only against the
  ≥5-UMI (or pooled) truth.
- Enrichment over the shuffled null (65–179×, empirical p at the 1/11 floor with 10 seeds) says the sites
  are not random; it does not rescue the target, which was an absolute fraction.
- The hexamer criterion adds +0.046 on the whole default output but nothing among atlas-novel sites
  (hexamer-FAIL 0.487 vs trusted-novel 0.485); the atlas-novelty stage costs −0.326. Say "the novelty
  criterion selects the sites the long reads see least", not "the hexamer filter validates the novel sites".
- Funnel stages clip-supported / ≥2 molecules / not-IP remove nothing because the pre-registered default
  already enforces them (the IP flag lives in `run/03_gtf_annotation/default/annotatedpas.bed:10`, not in the
  9-column top-level file; re-deriving IP from the genome flags 36 sites under the tool-default +10..+30
  window and 661 under the caller's own −10..+29 window — the raw-peak vs refined-cleavage-point offset —
  changing the trusted-novel count by ≤3%). PolyA_DB could not be
  used (hg19 only); "≥100 bp from any atlas site" means PolyASite 2.0 only.
- 27.0% of trusted-novel sites coincide with a Kinnex internal-priming decoy terminus (7.6% for atlas-known
  hexamer-pass): the caller's IP rule (−10..+30, ≥6 A or ≥70% A) is looser than Kinnex's (≥12 of 18 A in
  +1..+18). The residual false-positive class is internal priming the filter does not catch (issue 10
  addendum, `--ip-rule kinnex`, NOT superseded by #96). Decoy-proximal sites are less concordant
  (0.39 vs 0.52), so the decoy overlap is a plausible mechanism, not a proven one.
- Single PBMC donor, single CellRanger BAM; Kinnex x3p t5 truth is itself permissive (only ~15% of t5 peaks
  are within 25 bp of PolyASite); signed distances are asymmetric (Kinnex termini sit a few bp upstream of the
  PeakATail cleavage coordinate) — no effect within the 25-bp window.
- Panel a's percentages are relative to the 46,544-site default (x-axis label; 20 sites on unlisted contigs are
  excluded before the funnel, so the 46,524-site "input" row is 99.96%, shown as 100%); all concordance
  denominators are the on-contig counts. The 3,578-gene count (subtitle only) and the 46,544 denominator are
  recorded in `figS7_novelfunnel_reference_lines.tsv`.

## Provenance

Sources: `manuscript/19_final_gate_v2.md` §4 (v2 gate, verifier verdict FIXED) and
`manuscript/16_trusted_novel_kinnex.md` **§v2** (the v2 record for every number on this figure,
verified FIXED; the pre-§v2 body of 16 is the v1 record);
`results/reliability/trusted_novel_final_v2_pbmc/` (`REPORT_PROVISIONAL.md` STATUS: VERIFIED, adversarial
verifier 2026-08-21) — `call_primary/funnel.tsv`,
`validate_incl_genebodies/validation.tsv` + `validation_null_seeds.tsv`, `feature_breakdown.tsv`,
`distance_crosscheck.tsv`, `stratified_concordance_trusted_novel.tsv`, `verifier_crosscheck/RESULTS.txt`.
Tool: `scripts/reliability/trusted_novel_pas.py` (md5-identical to the v1 run).
The v1 figure remains reproducible: `TRUSTED_NOVEL_VERSION=v1 python3 scripts/manuscript_figures/figS7_novelfunnel.py`.
Palette: Okabe-Ito, six-checks validated 2026-08-21 (adjacent CVD ΔE ≥ 8.5; feature-class ramp validated
`--ordinal`).
"""
if VERSION == "v2":
    with open(f"{FIGDIR}/{NAME}.caption.md", "w") as fh:
        fh.write(CAPTION_MD)
    print("wrote", f"{FIGDIR}/{NAME}.caption.md")
else:
    print("caption sidecar NOT rewritten (v1 record render; the sidecar describes v2)")
