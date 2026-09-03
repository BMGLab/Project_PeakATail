#!/usr/bin/env python3
"""
figS3_nulldesign.py -- manuscript Fig S3 ("figS3_nulldesign"): why the
benchmark is built this way -- REBUILD of the retired `null_control` against
the curated PolyASite 2.0 point reference (21 D5 / 21 section 4.3-S4 /
05 section S3; executed 2026-09-02).

THE TWO MANDATORY PRECONDITIONS (21 section 4.3-S4, manifest S3 row), and how
each is honoured:

(i) NULL RECONCILIATION (gap 7).  The two retired figures shipped two
    conflicting nulls for the same arm (lg_annotate @100 bp, interval mode,
    vs the 18,432,135-entry PolyASite dump): `benchmark_strategies.py`
    0.753 (strand-agnostic window, genome-wide width-preserving shuffle,
    no -chrom) vs `null_control.py` 0.615 (strand-matched closest,
    same-chromosome shuffle; its gene-body variant 0.749).  05 R4/R5 record
    the definitional differences; 07 section 1 (verified SOUND) records that
    randomly placed width-matched peaks score 0.62-0.75 under that regime and
    that round 2 replaced both with ONE reconciled null: width-preserving
    `bedtools shuffle -chrom -noOverlapping -incl <merged gene bodies>`
    (1.80 Gb), 3 seeds, scored strand-matched in point mode against the
    curated reference.  That is byte-for-byte the null implemented in
    scripts/benchmark_tools/score_tool.py and used by EVERY verified score
    TSV (benchmark_curated.tsv, 5 arms x 3 seeds, 0.020-0.024; the
    fig2_accuracy.tsv null_P columns for every tool on all three datasets).
    This figure shows the retired pair ONLY as the superseded discrepancy
    being reconciled (panel a, marked RETIRED); the single null it ships is
    score_tool.py's, read from the verified TSVs.

(ii) DENSITY RECOMPUTATION.  Every density figure here is recomputed fresh
    on the curated PolyASite 2.0 point reference
    (data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6, 569,005
    representative sites, provenance md5-logged in that directory's README).
    The retired 23/33/45/61/73% and "one entry every 168 bp" figures are
    properties of the quarantined 18.4M dump and are NOT carried over
    (05 R5 binding caveat); the dump appears in this figure only as an entry
    count and a retired-regime label, never as a density.

DROPPED PANELS (house rule: a value with no verified source is dropped and
the reason recorded, never approximated):
  * match-class composition -- the only verified class-composition record
    (05 R5: TE 0.310 real vs 0.022 null etc.) is measured against the dump's
    col-10 classes and is forbidden by (ii); the curated rep-site BED6
    carries no class column, and no verified class analysis on the curated
    reference exists.  Unblocked by: one bedtools pass of the five arms
    against atlas.clusters.2.0.GRCh38.96.bed col 10 + verification.
  * quantitative interval-vs-point matching on the curated reference -- the
    old 0.9985 (interval) vs 0.9940 (point) pair is dump-derived; no
    interval-mode score against the curated reference has been run or
    verified.  The point is carried qualitatively by panel a's regime
    contrast (07 section 1's component table).

SOURCES OF TRUTH (nothing numeric typed by hand except EXPECT constants,
each cited to a verified document):
  panel a  results/figures/manuscript/benchmark_curated.tsv (405 rows;
           verified SOUND, 07 section 2 primary table) -- five arms real +
           3-seed null @100 bp; retired-regime constants from 07 section 1
           (0.62-0.75 floor; 18,432,135 entries; 0.996-1.000 / spread
           <0.004) and 05 R4/R5 (0.753 / 0.615 / 0.749 conflict values).
  panel b  computed HERE from the curated reference + chrom.sizes +
           results/benchmark_tools/shared_refs/genebodies.merged.bed
           (precondition ii recompute; every value in the density TSV).
  panel c  results/figures/manuscript/fig2_accuracy.tsv (verified v2 chain,
           19 FIXED) -- real P@100 and the score_tool.py 3-seed null per
           tool on PBMC 10k v3 + testis mouse 1/2.

OUTPUTS
  manuscript/figures/figS3_nulldesign.{png,pdf}   (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/figS3_nulldesign.caption.md  (sidecar, script-written;
      '## Legend' = the journal legend incl. every caveat moved off the image
      in the 2026-09-02 submission pass; '## Provenance' = record-keeping)
  results/figures/manuscript/figS3_nulldesign.tsv          (panels a+c values)
  results/figures/manuscript/figS3_nulldesign_density.tsv  (panel b values)

Run:  export LC_ALL=C OMP_NUM_THREADS=1; python3 scripts/manuscript_figures/figS3_nulldesign.py
Cache: results/figures/manuscript/.cache_figS3_nulldesign/ (delete to force
       a full density recompute; <=1 thread throughout, bedtools+awk stream).
"""
import os
import re
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pubstyle import PAL, TOOL_STYLE, TYPE, apply_rc, plain_log, sentence_case   # shared publication style

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
NAME = "figS3_nulldesign"
CACHE = OUTDIR / f".cache_{NAME}"
for d in (OUTDIR, FIGDIR, CACHE):
    d.mkdir(parents=True, exist_ok=True)

ATLAS = WD / "data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6"
GENOME = WD / "data/references/chrom.sizes.nochr.filt"
GENIC = WD / "results/benchmark_tools/shared_refs/genebodies.merged.bed"
BC_TSV = OUTDIR / "benchmark_curated.tsv"
F2_TSV = OUTDIR / "fig2_accuracy.tsv"

ENV = dict(os.environ, LC_ALL="C")  # tr_TR locale silently corrupts BED sorting

INK, MUTED, GRID = "#1B2429", "#5A6B73", "#D8E0E3"
# Okabe-Ito colourblind-safe palette (figures/README.md convention); fixed
# entity assignment, never cycled.  Panel c tools reuse fig2_accuracy.py's map.
C_REAL = "#0072B2"      # real calls under the curated regime
C_NULL = "#009E73"      # the one reconciled (score_tool.py) null + its density
C_RETIRED = "#D55E00"   # the retired dump regime
C_GENOME = "#999999"    # genome-wide context curve
# Panel c tools take the ONE tool identity of `_pubstyle.TOOL_STYLE` (the map Figs 2
# and 3 draw from) instead of a local restatement of it: the restated copy had
# scUTRquant in light blue where the shared map makes it neutral grey (judge finding 1).
TOOL_COLOR = {tool: st["color"] for tool, st in TOOL_STYLE.items()}

# ---------------------------------------------------------------------------
# EXPECT constants -- every one cited to a verified document; used ONLY in
# asserts and retired-regime annotation, never silently as plotted data
# ---------------------------------------------------------------------------
EXP_N_ATLAS = 569005          # 07 section 1 (curated cluster count; file recounted below)
EXP_N_DUMP = 18432135         # 07 section 1 (the quarantined dump's entry count)
EXP_GENIC_GB = 1.80           # 07 section 1 (merged gene bodies, the null's space)
EXP_GENOME_GB = 3.088         # 05 R4 verified record ("3.088 Gb genome")
DUMP_REAL_LO, DUMP_REAL_HI = 0.996, 1.000   # 07 section 5 item 5 / 05 R4 (FIXED)
DUMP_SPREAD_MAX = 0.004                     # 05 R4 ("full spread < 0.004") / manifest S3
DUMP_NULL_LO, DUMP_NULL_HI = 0.62, 0.75     # 07 section 1 (verified SOUND)
CONFLICT = [  # the gap-7 pair (+ null_control's own second model), all RETIRED
    (0.615, "chrom-shuffle, strand-matched closest (null_control)",
     "manuscript/05_figure_index.md R5 (FIXED record)"),
    (0.749, "gene-body shuffle, strand-matched closest (null_control)",
     "manuscript/05_figure_index.md R5 (FIXED record)"),
    (0.753, "genome shuffle, strand-agnostic window (benchmark_strategies)",
     "manuscript/05_figure_index.md R4/R5 (FIXED records; gap 7)"),
]
EXP_ARM_P100 = {  # 07 section 2 primary table (verified SOUND), 3 dp
    "lg_annotate": 0.447, "lp_annotate": 0.492, "si_annotate": 0.355,
    "lg_ip_filter": 0.457, "B1_cohort_full": 0.475}
EXP_NULL_LO, EXP_NULL_HI = 0.020, 0.024     # 07 section 2 (5 arms x 3 seeds @100 bp)
EXP_SPREAD_2DP = 0.14                       # 07 section 2 ("spread 0.14")
EXP_DEFAULT_P = {"pbmc": 0.7062, "mouse1": 0.7450, "mouse2": 0.7572}  # 19 section 1
EXP_PBMC_NULL_4DP = 0.0217                  # 19/21 L1 (default arm 3-seed mean)
EXP_MOUSE_SEED_RANGE = (0.0127, 0.0150)     # 21 L1 (testis default per-seed range)
CUTOFFS = [10, 25, 50, 100, 200, 500, 1000]
I100 = CUTOFFS.index(100)


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


# ---------------------------------------------------------------------------
# panel b: DENSITY, recomputed on the curated point reference (precondition ii)
# fraction of (genome | merged gene-body) bp within d of a same-strand curated
# site; per strand, then averaged (point-mode matching is strand-matched, so a
# base can only be claimed by a same-strand site)
# ---------------------------------------------------------------------------
GENOME_BP = sum(int(l.split("\t")[1]) for l in GENOME.read_text().strip().split("\n"))
GENIC_BP = int(sh_out(f"awk '{{s+=$3-$2}}END{{print s}}' {GENIC}").strip())
N_ATLAS_FILE = int(sh_out(f"wc -l < {ATLAS}").strip())
N_ON_GENOME = int(sh_out(
    f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {GENOME} {ATLAS} | wc -l").strip())

assert N_ATLAS_FILE == EXP_N_ATLAS, (N_ATLAS_FILE, EXP_N_ATLAS)
assert abs(GENOME_BP / 1e9 - EXP_GENOME_GB) < 0.005, GENOME_BP
assert round(GENIC_BP / 1e9, 2) == EXP_GENIC_GB, GENIC_BP
N_SCAFFOLD = N_ATLAS_FILE - N_ON_GENOME   # 07 section 2: 397 scaffold-only sites
assert N_SCAFFOLD == 397, N_SCAFFOLD

DENS_CACHE = CACHE / "density.tsv"
if not DENS_CACHE.exists():
    rows = []
    for sym, strand in (("+", "plus"), ("-", "minus")):
        for cut in CUTOFFS:
            merged = CACHE / "tmp_merged.bed"
            # pad each same-strand point by +-cut, clamp to [0, chrom_len], merge
            sh_out(
                "awk -F'\\t' 'NR==FNR{len[$1]=$2;next} "
                f"($1 in len) && $6==\"{sym}\" "
                f"{{s=$2-{cut}; if(s<0)s=0; e=$3+{cut}; if(e>len[$1])e=len[$1]; "
                "print $1\"\\t\"s\"\\t\"e}' "
                f"{GENOME} {ATLAS} | bedtools merge -i - > {merged}")
            gw = int(sh_out(f"awk '{{s+=$3-$2}}END{{print s}}' {merged}").strip())
            gb = int(sh_out(f"bedtools intersect -a {merged} -b {GENIC} "
                            f"| awk '{{s+=$3-$2}}END{{print s+0}}'").strip())
            rows.append((strand, cut, gw, gb))
            merged.unlink()
    with open(DENS_CACHE, "w") as fh:
        fh.write("strand\tcutoff_bp\tcovered_bp_genome\tcovered_bp_genic\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
dens = pd.read_csv(DENS_CACHE, sep="\t")
assert len(dens) == 2 * len(CUTOFFS)
dmean = dens.groupby("cutoff_bp").mean(numeric_only=True).reindex(CUTOFFS)
frac_gw = (dmean.covered_bp_genome / GENOME_BP).values   # genome-wide density
frac_gb = (dmean.covered_bp_genic / GENIC_BP).values     # gene-body density
# structural sanity: monotone in cutoff; genic denser than genome-wide
assert (np.diff(frac_gw) > 0).all() and (np.diff(frac_gb) > 0).all()
assert (frac_gb > frac_gw).all()
BP_PER_SITE = GENOME_BP / N_ON_GENOME   # replaces the forbidden "every 168 bp"

# ---------------------------------------------------------------------------
# panel a: the five arms under the curated regime (benchmark_curated.tsv,
# verified SOUND by 07 section 2) + the retired dump regime as constants
# ---------------------------------------------------------------------------
bc = pd.read_csv(BC_TSV, sep="\t")
bc100 = bc[(bc.panel == "precision") & (bc.cutoff_bp == 100.0)
           & (bc.reference == "atlas_full")]
arm_real = (bc100[bc100.series == "real"]
            .set_index("arm")["value"].to_dict())
null_seeds_bc = bc100[bc100.series == "null_genic"]["value"].values
assert set(arm_real) == set(EXP_ARM_P100)
for k, v in EXP_ARM_P100.items():
    assert round(arm_real[k], 3) == v, (k, arm_real[k], v)
assert len(null_seeds_bc) == 15  # 5 arms x 3 seeds
assert null_seeds_bc.min() >= EXP_NULL_LO and null_seeds_bc.max() <= EXP_NULL_HI, (
    null_seeds_bc.min(), null_seeds_bc.max())
CUR_SPREAD = max(arm_real.values()) - min(arm_real.values())
assert round(CUR_SPREAD, 2) == EXP_SPREAD_2DP, CUR_SPREAD

# consistency observation (reported, and softly asserted): the reconciled null
# is random genic points, so it should score ~ the gene-body density at 100 bp
null_mean_bc = float(null_seeds_bc.mean())
dens_ratio = null_mean_bc / frac_gb[I100]
assert 0.5 < dens_ratio < 2.0, (null_mean_bc, frac_gb[I100])

# ---------------------------------------------------------------------------
# panel c: the reconciled null across the verified chain (fig2_accuracy.tsv)
# ---------------------------------------------------------------------------
f2 = pd.read_csv(F2_TSV, sep="\t")
f2 = f2[f2.cutoff_bp == 100]
BLOCKS = [("pbmc_10k_v3", "single", "PBMC 10k v3 (human)"),
          ("gse104556", "mouse1", "testis mouse 1"),
          ("gse104556", "mouse2", "testis mouse 2")]
TOOL_ORDER = [  # fixed order, default first, catalog tool last (hollow)
    ("PeakATail", "path4", "PeakATail default"),
    ("SCAPTURE", "competitor", "SCAPTURE"),
    ("polyApipe", "competitor", "polyApipe"),
    ("Sierra", "competitor", "Sierra"),
    ("scAPAtrap", "competitor", "scAPAtrap"),
    ("scUTRquant", "catalog", "scUTRquant*"),
]
chain = []
for ds, rep, block_lab in BLOCKS:
    for tool, role, lab in TOOL_ORDER:
        q = f2[(f2.dataset == ds) & (f2.replicate == rep)
               & (f2.tool == tool) & (f2.role == role)]
        if rep == "mouse2" and tool == "SCAPTURE":
            # fig2_accuracy convention (21 map amendment, RESOLVED 2026-09-02):
            # SCAPTURE mouse-2 completed at site level (24,076 points, P@100
            # 0.672, 15 section 3) but only mouse 1 carries a verified
            # null-bearing row; noted on-figure and in the caption, not plotted.
            assert len(q) == 0, "SCAPTURE mouse2 row appeared; drop this skip"
            continue
        assert len(q) == 1, (ds, rep, tool, role, len(q))
        r = q.iloc[0]
        chain.append(dict(
            dataset=ds, block=block_lab, tool=tool, label=lab, role=role,
            real=float(r.atlas_agreement_precision_100bp),
            null_mean=float(r.null_P_mean),
            null_seeds=[float(r.null_P_seed1), float(r.null_P_seed2),
                        float(r.null_P_seed3)],
            src=str(r.source)))
ch = pd.DataFrame(chain)
defaults = ch[ch.role == "path4"].set_index("block")
assert round(defaults.loc["PBMC 10k v3 (human)", "real"], 4) == EXP_DEFAULT_P["pbmc"]
assert round(defaults.loc["testis mouse 1", "real"], 4) == EXP_DEFAULT_P["mouse1"]
assert round(defaults.loc["testis mouse 2", "real"], 4) == EXP_DEFAULT_P["mouse2"]
assert round(defaults.loc["PBMC 10k v3 (human)", "null_mean"], 4) == EXP_PBMC_NULL_4DP
m_seeds = np.array([s for b in ("testis mouse 1", "testis mouse 2")
                    for s in defaults.loc[b, "null_seeds"]])
assert (round(m_seeds.min(), 4), round(m_seeds.max(), 4)) == EXP_MOUSE_SEED_RANGE, (
    m_seeds.min(), m_seeds.max())
folds = (defaults.real / defaults.null_mean)
assert 32 < folds.min() < 34 and 53 < folds.max() < 56, folds.to_dict()  # 19: "33-55x"

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# canvas: the pre-submission 8.3 x 8.9 in canvas carried a title/intro band on
# top and a 15-line footer caption below; both now live in the caption sidecar
# ('## Legend'), and the canvas height drops by the freed space.  Width stays
# 8.3 in (211 mm): panels a/b carry right-aligned annotation blocks outside
# their axes and panel c a 6-row per-block dot plot -- narrowing to the 180 mm
# double-column norm crowds them (legibility outranks the width target).
fig = plt.figure(figsize=(8.3, 6.1))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.02], height_ratios=[1.05, 1.0],
                      left=0.085, right=0.975, top=0.925, bottom=0.090,
                      hspace=0.52, wspace=0.34)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[:, 1])


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=TYPE["tick"])
    ax.set_axisbelow(True)


# ---- panel a: the dump could not rank; the curated point reference can -----
axA.yaxis.grid(True, color=GRID, lw=0.5, alpha=0.8)
X_DUMP, X_CUR = 0.0, 1.0
BW = 0.62
# retired regime: real band + null floor band + the three conflicting nulls
axA.add_patch(Rectangle((X_DUMP - BW / 2, DUMP_REAL_LO), BW,
                        DUMP_REAL_HI - DUMP_REAL_LO, facecolor=C_RETIRED,
                        alpha=0.85, edgecolor="none", zorder=3))
axA.add_patch(Rectangle((X_DUMP - BW / 2, DUMP_NULL_LO), BW,
                        DUMP_NULL_HI - DUMP_NULL_LO, facecolor=C_RETIRED,
                        alpha=0.22, edgecolor="none", zorder=2))
for v, lab, _src in CONFLICT:
    axA.plot([X_DUMP - BW / 2, X_DUMP + BW / 2], [v, v], color=C_RETIRED,
             lw=1.1, ls=(0, (4, 2)), zorder=4)
axA.text(X_DUMP - BW / 2 - 0.06, 0.998, "five strategies\n0.996–1.000\nspread < 0.004",
         ha="right", va="center", fontsize=6.0, color=C_RETIRED, fontweight="bold",
         linespacing=1.3)
axA.text(X_DUMP - BW / 2 - 0.06, 0.687,
         "\"null\" 0.62–0.75:\nmoves 0.14 with\narbitrary definition\nchoices (gap 7)",
         ha="right", va="center", fontsize=6.0, color=C_RETIRED, linespacing=1.3)
axA.annotate("0.615", (X_DUMP + BW / 2, 0.615), xytext=(X_DUMP + BW / 2 + 0.05, 0.585),
             fontsize=ANN, color=C_RETIRED, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=C_RETIRED, lw=0.5))
axA.annotate("0.749 / 0.753", (X_DUMP + BW / 2, 0.751),
             xytext=(X_DUMP + BW / 2 + 0.05, 0.79), fontsize=ANN, color=C_RETIRED,
             ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=C_RETIRED, lw=0.5))
# curated regime: five arms + 15 null seeds
rng = np.random.default_rng(7)
arm_order = sorted(arm_real, key=arm_real.get)
jit = np.linspace(-0.16, 0.16, len(arm_order))
ARM_LAB = {"lg_annotate": "lg", "lp_annotate": "lp", "si_annotate": "si",
           "lg_ip_filter": "lg+IP", "B1_cohort_full": "B1"}
for i, (x0, arm) in enumerate(zip(jit, arm_order)):
    v = arm_real[arm]
    axA.plot(X_CUR + x0, v, "o", ms=5.5, mfc=C_REAL, mec="white", mew=0.8, zorder=5)
    dy = 4.5 if i % 2 == 0 else 12.5   # stagger labels: the arms sit 0.02 apart
    axA.annotate(ARM_LAB[arm], (X_CUR + x0, v), xytext=(0, dy),
                 textcoords="offset points", fontsize=ANN, color=C_REAL,
                 ha="center", va="bottom")
axA.plot(X_CUR + rng.uniform(-0.16, 0.16, len(null_seeds_bc)), null_seeds_bc,
         "o", ms=3.6, mfc=C_NULL, mec="white", mew=0.6, alpha=0.9, zorder=5)
# spread bracket for the curated arms
xb = X_CUR + 0.30
axA.plot([xb, xb], [min(arm_real.values()), max(arm_real.values())],
         color=INK, lw=0.9, zorder=4)
for y in (min(arm_real.values()), max(arm_real.values())):
    axA.plot([xb - 0.025, xb], [y, y], color=INK, lw=0.9, zorder=4)
axA.text(xb + 0.05, 0.42, f"spread {CUR_SPREAD:.2f}\nranks the\nstrategies",
         ha="left", va="center", fontsize=6.0, color=INK, fontweight="bold",
         linespacing=1.3)
axA.text(X_CUR - 0.42, 0.062, "one reconciled null\n(3 seeds × 5 arms)\n0.020–0.024",
         ha="right", va="center", fontsize=6.0, color=C_NULL, linespacing=1.3)
axA.set_xlim(-1.05, 1.75)
axA.set_ylim(0, 1.06)
axA.set_xticks([X_DUMP, X_CUR])
axA.set_xticklabels([sentence_case("RETIRED regime:\nunfiltered PolyASite file\n(18.4M entries),\nwhole-interval match"),
                     sentence_case("curated PolyASite 2.0\npoints, strand-matched\n3′ base, one null")],
                    fontsize=TYPE["tick"])
axA.get_xticklabels()[0].set_color(C_RETIRED)
axA.get_xticklabels()[1].set_color(INK)
axA.set_ylabel(sentence_case("atlas-agreement precision @100 bp"))
# Short identifier, not a sentence (judge finding 4): the claim it used to state
# on the image -- "the dump could not rank the strategies; the curated point
# reference can" -- is a sentence of the Legend below, verbatim, and the internal
# word "dump" is replaced everywhere on the image by what that reference set is:
# the unfiltered 18,432,135-entry PolyASite download (07 section 1).
axA.set_title(sc_title("a   Two reference regimes"), loc="left", fontweight="bold")
style(axA)

# ---- panel b: density recomputed on the curated reference ------------------
axB.yaxis.grid(True, color=GRID, lw=0.5, alpha=0.8)
Xd = np.arange(len(CUTOFFS))
axB.plot(Xd, frac_gb, color=C_NULL, lw=1.8, marker="o", ms=5, mec="white",
         mew=0.8, zorder=4, label="within merged gene bodies (the null's space)")
axB.plot(Xd, frac_gw, color=C_GENOME, lw=1.6, marker="s", ms=4.2, mec="white",
         mew=0.8, zorder=3, label="genome-wide")
axB.plot(np.full(len(null_seeds_bc), Xd[I100])
         + rng.uniform(-0.1, 0.1, len(null_seeds_bc)), null_seeds_bc, "o",
         ms=3.4, mfc="white", mec=C_NULL, mew=1.0, zorder=5,
         label="measured null (15 seeds, human arms)")
axB.set_yscale("log")
axB.set_ylim(8e-4, 1.0)
plain_log(axB, "y")     # 1 / 0.1 / 0.01 / 0.001, not 10^n: a mathtext exponent lands at
                        # 4.55 pt under the 6.5 pt tick size (judge finding 3)
axB.set_xticks(Xd)
axB.set_xticklabels([str(c) for c in CUTOFFS])
axB.set_xlabel(sentence_case("distance to nearest same-strand curated site (bp)"))
axB.set_ylabel(sentence_case("fraction of bp within reach (log)"))
# the parenthetical subtitle line '(recomputed here; dump densities retired)' moved
# to the Legend sidecar (directive 2, no-loss: it is the '(b) Reference density
# recomputed on the curated point reference (dump densities retired)' sentence)
axB.set_title(sc_title("b   chance = local density of the curated reference"),
              loc="left", fontweight="bold")
axB.annotate(f"@100 bp: {frac_gb[I100]:.4f} of gene-body bp\n"
             f"measured null mean {null_mean_bc:.4f}\n"
             f"(ratio {dens_ratio:.2f})",
             (Xd[I100], frac_gb[I100]), xytext=(Xd[I100] - 0.35, 0.125),
             fontsize=ANN, color=C_NULL, ha="right", va="center", linespacing=1.35,
             arrowprops=dict(arrowstyle="-", color=C_NULL, lw=0.6))
# short data label only; the scaffold-exclusion and retired-dump-count
# sentences moved to the legend (2026-09-02 submission pass)
axB.text(0.985, 0.04,
         f"{N_ON_GENOME:,} sites on the {GENOME_BP/1e9:.2f} Gb genome\n"
         f"— one per {BP_PER_SITE/1e3:.1f} kb",
         transform=axB.transAxes, fontsize=ANN, color=MUTED, ha="right",
         va="bottom", linespacing=1.4)
axB.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), frameon=False,
           fontsize=ANN, handletextpad=0.4, labelspacing=0.35)
style(axB)

# ---- panel c: the reconciled null across the verified chain ----------------
axC.xaxis.grid(True, color=GRID, lw=0.5, alpha=0.8)
y = 0
yticks, ylabs, ycols = [], [], []
block_bounds = []
for ds, rep, block_lab in BLOCKS:
    y_start = y
    sub = ch[ch.block == block_lab]
    for _, r in sub.iterrows():
        col = TOOL_COLOR[r.tool]
        hollow = r.role == "catalog"
        axC.plot([r.null_mean, r.real], [y, y], color=col, lw=1.4, alpha=0.7,
                 solid_capstyle="round", zorder=2)
        axC.plot(r.null_mean, y, "o", ms=4.2, mfc="white", mec=col, mew=1.2, zorder=3)
        axC.plot(r.real, y, "o", ms=5.2, mfc="white" if hollow else col,
                 mec=col if hollow else "white", mew=1.4 if hollow else 0.8, zorder=4)
        if r.role == "path4":
            axC.annotate(f"{r.real/r.null_mean:.0f}× its null",
                         (r.real, y), xytext=(2, 4), textcoords="offset points",
                         fontsize=6.0, color=col, fontweight="bold",
                         ha="right", va="bottom")
        yticks.append(y)
        ylabs.append(r.label)
        ycols.append(col)
        y += 1
    block_bounds.append((block_lab, y_start, y - 1))
    y += 1.4  # gap between dataset blocks
for block_lab, y0, y1 in block_bounds:
    axC.text(0.0115, y1 + 0.85, block_lab, fontsize=6.6, color=INK,
             fontweight="bold", ha="left", va="center")
axC.set_yticks(yticks)
axC.set_yticklabels(ylabs, fontsize=TYPE["tick"])
for t, c in zip(axC.get_yticklabels(), ycols):
    t.set_color(c)
axC.set_xscale("log")
axC.set_xlim(0.009, 1.05)
axC.set_xticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
axC.set_xticklabels(["0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1"])
axC.set_ylim(y + 0.5, -0.9)   # inverted: first block on top; bottom pad for notes
axC.set_xlabel(sentence_case("atlas-agreement precision @100 bp (log)"))
axC.set_title(sc_title("c   one null everywhere: every verified score TSV\n"
                       "      carries score_tool.py's 3-seed genic shuffle"),
              loc="left", fontweight="bold")
axC.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=4.2, mfc="white", mec=MUTED, mew=1.2,
           label="null (3-seed mean)"),
    Line2D([], [], marker="o", ls="", ms=5.2, mfc=MUTED, mec="white", mew=0.8,
           label="real call set"),
    Line2D([], [], marker="o", ls="", ms=5.2, mfc="white", mec=MUTED, mew=1.4,
           label="catalog-based (not ranked)"),
], loc="lower right", bbox_to_anchor=(1.0, -0.005), frameon=False, fontsize=ANN,
    handletextpad=0.35, labelspacing=0.3)
# (the multi-sentence null-definition note that sat here moved to the legend
# in the 2026-09-02 submission pass; its content is in the caption below)
style(axC)

# ---------------------------------------------------------------------------
# the journal legend (single source: the caption sidecar's '## Legend').
# The former on-figure title, intro paragraph, italic retired-record note and
# 15-line footer all live here now, substance verbatim.
# ---------------------------------------------------------------------------
caption = (
    "Figure S3 | Why the benchmark is built this way: curated points, one reconciled null "
    "(rebuild of the retired null_control against the curated reference only). "
    "(a) Two reference regimes, and only one of them can rank the strategies. "
    f"Against the unfiltered {EXP_N_DUMP:,}-entry PolyASite download (the retired 'dump' reference) "
    "with whole-peak-interval matching, the five "
    f"peak-calling strategies scored {DUMP_REAL_LO:.3f}-{DUMP_REAL_HI:.3f} (spread < {DUMP_SPREAD_MAX}), "
    "while width-matched shuffled peaks already scored 0.62-0.75 (07 1) - and the two retired figures "
    "computed different nulls for the same arm (0.615 strand-matched chrom-shuffle vs 0.753 strand-agnostic "
    "genome-shuffle - a 0.14 disagreement produced by definition choices alone; gene-body variant 0.749; "
    "gap 7): under that regime the benchmark can rank nothing. "
    "Scored at the strand-matched 3'-most base against the curated PolyASite 2.0 representative sites "
    f"({EXP_N_ATLAS:,} clusters, GRCh38.96), the same five arms separate (P@100 "
    f"{min(arm_real.values()):.3f}-{max(arm_real.values()):.3f}, spread {CUR_SPREAD:.2f}) over one "
    f"reconciled null ({null_seeds_bc.min():.3f}-{null_seeds_bc.max():.3f}, 3 seeds x 5 arms). "
    f"(b) Reference density recomputed on the curated point reference (dump densities retired): "
    f"{frac_gb[I100]*100:.2f}% of gene-body bp lies within 100 bp of a same-strand curated site, and the "
    f"measured null mean ({null_mean_bc:.4f}) equals that local density within {abs(dens_ratio-1)*100:.0f}% "
    f"- chance in this benchmark IS reference density, which is why the null is indispensable. One site per "
    f"{BP_PER_SITE/1e3:.1f} kb of genome ({N_ON_GENOME:,} sites; {N_SCAFFOLD} scaffold-only excluded; the "
    f"retired dump held {EXP_N_DUMP / N_ON_GENOME:.0f}x as many entries). "
    "(c) The same score_tool.py null (width-preserving 3-seed shuffle inside 1.80 Gb of merged gene bodies, "
    "scored identically to the real calls) accompanies every verified score TSV: the pre-registered "
    f"PeakATail default scores {EXP_DEFAULT_P['pbmc']:.4f} / {EXP_DEFAULT_P['mouse1']:.4f} / "
    f"{EXP_DEFAULT_P['mouse2']:.4f} (PBMC / mouse 1 / mouse 2) against nulls of "
    f"{defaults.loc['PBMC 10k v3 (human)', 'null_mean']:.4f} / "
    f"{defaults.loc['testis mouse 1', 'null_mean']:.4f} / "
    f"{defaults.loc['testis mouse 2', 'null_mean']:.4f} - "
    f"{folds.min():.0f}-{folds.max():.0f}x above chance - and every competitor is scored with its own "
    "identically-built null. scUTRquant (*) is catalog-based and not ranked; SCAPTURE mouse 2 is not "
    "plotted (site-level run scored 0.672 (15 3) but only mouse 1 carries a verified null-bearing row - "
    "fig2_accuracy convention). The rebuilt scoring is at the strand-matched 3'-most base; retired values "
    "(a, left) appear only as the superseded record (05 R4/R5, 07 1), no dump density carried over, and "
    "the match-class panel is omitted (the only verified record is dump-derived)."
)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
p = OUTDIR / f"{NAME}.png"
fig.savefig(p, dpi=600)
print("wrote", p)

# ---------------------------------------------------------------------------
# TSVs: every plotted value, with a source column
# ---------------------------------------------------------------------------
SRC_BC = "results/figures/manuscript/benchmark_curated.tsv (verified SOUND, 07 §2)"
SRC_07 = "manuscript/07_curated_benchmark_report.md §1 (verified SOUND; retired-regime constant)"
SRC_05 = "manuscript/05_figure_index.md R4/R5 (FIXED records; retired-regime constant, gap 7)"
rows = []
rows.append(dict(panel="a", record="dump_regime_real_band", metric="P100_interval_dump",
                 value=f"{DUMP_REAL_LO}-{DUMP_REAL_HI}", detail="five strategies, spread<0.004; RETIRED",
                 source=SRC_07))
rows.append(dict(panel="a", record="dump_regime_null_band", metric="null_P100_interval_dump",
                 value=f"{DUMP_NULL_LO}-{DUMP_NULL_HI}",
                 detail="width-matched shuffled peaks under the dump regime; RETIRED", source=SRC_07))
for v, lab, src in CONFLICT:
    rows.append(dict(panel="a", record="retired_null_conflict", metric="null_P100_interval_dump",
                     value=v, detail=lab + "; RETIRED, reconciled by score_tool.py's null",
                     source=src))
for arm, v in arm_real.items():
    rows.append(dict(panel="a", record=arm, metric="P100_point_curated", value=v,
                     detail=f"label {ARM_LAB[arm]}", source=SRC_BC))
for i, v in enumerate(null_seeds_bc):
    rows.append(dict(panel="a,b", record="curated_null_seed", metric="null_P100_point_curated",
                     value=v, detail=f"seed point {i+1}/15 (5 arms x 3 seeds)", source=SRC_BC))
rows.append(dict(panel="a", record="curated_spread", metric="P100_spread_5_arms",
                 value=CUR_SPREAD, detail="max-min across the five arms", source=SRC_BC))
for _, r in ch.iterrows():
    rows.append(dict(panel="c", record=f"{r.block}|{r.label}", metric="P100_point_curated",
                     value=r.real, detail=f"role {r.role}", source=r.src + " via fig2_accuracy.tsv (19 FIXED)"))
    rows.append(dict(panel="c", record=f"{r.block}|{r.label}", metric="null_P100_mean_3seeds",
                     value=r.null_mean, detail="seeds " + "/".join(f"{s:.6f}" for s in r.null_seeds),
                     source=r.src + " via fig2_accuracy.tsv (19 FIXED)"))
for b, f in folds.items():
    rows.append(dict(panel="c", record=f"{b}|PeakATail default", metric="fold_over_null",
                     value=f, detail="real / null_mean (19: 33-55x below)",
                     source="derived from the two rows above"))
p = OUTDIR / f"{NAME}.tsv"
pd.DataFrame(rows).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

SRC_DENS = ("computed 2026-09-02 from data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6 "
            "(569,005 sites, provenance md5-logged in atlases/README.md) + chrom.sizes.nochr.filt + "
            "results/benchmark_tools/shared_refs/genebodies.merged.bed -- precondition (ii) recompute; "
            "dump-derived densities forbidden (05 R5)")
drows = []
for _, r in dens.iterrows():
    for scope, cov, denom in (("genome", r.covered_bp_genome, GENOME_BP),
                              ("gene_bodies", r.covered_bp_genic, GENIC_BP)):
        drows.append(dict(panel="b", strand=r.strand, scope=scope,
                          cutoff_bp=int(r.cutoff_bp), covered_bp=int(cov),
                          denominator_bp=denom, fraction=cov / denom, source=SRC_DENS))
for cut, gw, gb in zip(CUTOFFS, frac_gw, frac_gb):
    drows.append(dict(panel="b", strand="mean_of_strands", scope="genome", cutoff_bp=cut,
                      covered_bp="", denominator_bp=GENOME_BP, fraction=gw, source=SRC_DENS))
    drows.append(dict(panel="b", strand="mean_of_strands", scope="gene_bodies", cutoff_bp=cut,
                      covered_bp="", denominator_bp=GENIC_BP, fraction=gb, source=SRC_DENS))
drows.append(dict(panel="b", strand="both", scope="genome", cutoff_bp="",
                  covered_bp=N_ON_GENOME, denominator_bp=GENOME_BP,
                  fraction=BP_PER_SITE, source=SRC_DENS +
                  " | fraction column here = bp of genome per curated site (replaces the retired 'every 168 bp')"))
p = OUTDIR / f"{NAME}_density.tsv"
pd.DataFrame(drows).to_csv(p, sep="\t", index=False, float_format="%.6f")
print("wrote", p)

# ---------------------------------------------------------------------------
# sidecar caption
# ---------------------------------------------------------------------------
cap_md = f"""# Fig S3 — `figS3_nulldesign` caption (generated by `scripts/manuscript_figures/figS3_nulldesign.py`; rebuild of retired `null_control` against the curated reference, executed 2026-09-02; submission pass 2026-09-02: all on-figure prose — title, intro, retired-record note, footer — moved into the Legend below, image at 600 dpi)

## Legend

{caption}

**Caveats that travel with this figure (legend content, binding):**
- Panel a's left group is the RETIRED record (05 R4/R5, verdicts FIXED), shown solely as the discrepancy
  being reconciled; none of its values may be quoted as benchmark results.
- The null-equals-density reading of panel b (ratio {dens_ratio:.2f} at 100 bp) is a consistency
  observation of this recompute, not an independently verified claim; the small departure from 1 reflects
  the width-preserving, non-overlapping shuffle geometry and gene-body edge effects (07 §5 item 6).
- Panel c mouse rows are scored against the mouse curated atlas via the same `score_tool.py` machinery
  (same shuffle, same matcher); folds are computed against each call set's own 3-seed null mean.
- scUTRquant is catalog-based (hollow marker) and not ranked (fig2_accuracy conventions).
- The gene-body density denominator is the merged annotated gene-body space (1.80 Gb) actually used by
  the null shuffle — not the whole genome; the genome-wide curve is context only.

## Provenance

**Binding preconditions honoured (21 §4.3-S4 / manifest S3 row):**
1. *Null reconciliation (gap 7).* The conflicting 0.753 (benchmark_strategies: strand-agnostic window,
   genome-wide shuffle) and 0.615 (null_control: strand-matched closest, same-chromosome shuffle; its
   gene-body variant 0.749) were both interval-mode nulls against the quarantined 18.4M dump — 07 §1
   (SOUND) records the 0.62–0.75 floor of that regime and replaces both with ONE reconciled null
   (width-preserving `bedtools shuffle -chrom -noOverlapping -incl` merged gene bodies, 3 seeds,
   point-mode strand-matched scoring), implemented in `scripts/benchmark_tools/score_tool.py` and read
   here from the verified score TSVs. The retired pair appears only as the superseded record (panel a,
   RETIRED marks); the figure ships exactly one null.
2. *Density recompute.* Every density plotted (panel b) is recomputed on the curated PolyASite 2.0 point
   reference; the retired 23/33/45/61/73% and "one entry every 168 bp" dump figures are NOT carried
   over. The dump appears only as an entry count ({EXP_N_DUMP:,}) and a retired-regime label.

**Panels dropped, with reasons (no verified source; nothing approximated):**
- *Match-class composition* — the only verified record (05 R5) is measured against the dump's class
  column and is forbidden by precondition (ii); the curated rep-site BED6 has no class column. Unblock:
  one bedtools pass of the arms against `atlas.clusters.2.0.GRCh38.96.bed` col 10, plus verification.
- *Quantitative interval-vs-point matching on the curated reference* — the old 0.9985/0.9940 pair is
  dump-derived; no interval-mode score against the curated reference has been run. The point is carried
  qualitatively by panel a's regime contrast (07 §1 component table).

Sources: `{SRC_BC}`; `fig2_accuracy.tsv` (19 FIXED, v2 chain); density `{SRC_DENS}`;
retired-regime constants `{SRC_07}` and `{SRC_05}`.
Every plotted value: `results/figures/manuscript/figS3_nulldesign.tsv` and `figS3_nulldesign_density.tsv`.
"""
cap_md += """
**Design pass 2026-09-03** (`manuscript/figures/DESIGN_DIRECTIVES.md`, supplement light pass). Type comes from
the shared style module `scripts/manuscript_figures/_pubstyle.py` (`apply_rc()`), and every on-figure annotation
sits at or above the 6 pt floor. Axis labels, panel titles and prose tick labels are sentence-cased through
`_pubstyle.sentence_case()` with canonical identifiers preserved and the lower-case panel letters kept. Panel b's
parenthetical subtitle line *(recomputed here; dump densities retired)* left the image for the Legend above
(no-loss: it is the \"(b) Reference density recomputed on the curated point reference (dump densities retired)\"
sentence). No panel, number or audit TSV changed; both TSVs regenerate byte-identical.
"""
p = FIGDIR / f"{NAME}.caption.md"
p.write_text(cap_md)
print("wrote", p)

# ---------------------------------------------------------------------------
# edge check: the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
from PIL import Image
im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
edge = 8
border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                         im[:, :edge].ravel(), im[:, -edge:].ravel()])
ink = int((border < 250).sum())
print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"

print(f"figS3_nulldesign: curated arms {min(arm_real.values()):.3f}-{max(arm_real.values()):.3f} "
      f"(spread {CUR_SPREAD:.3f}) vs null {null_seeds_bc.min():.4f}-{null_seeds_bc.max():.4f}; "
      f"gene-body density@100 {frac_gb[I100]:.4f} vs null mean {null_mean_bc:.4f} "
      f"(ratio {dens_ratio:.2f}); defaults {folds.min():.1f}-{folds.max():.1f}x over null")
