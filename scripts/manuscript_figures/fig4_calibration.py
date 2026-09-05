#!/usr/bin/env python3
"""
fig4_calibration.py -- manuscript Fig 4: re-calibration of `peakatail switch diff` on
a CORRECT input.

Why v2: the first calibration (fdr_calibration.py, retired) ran on a mis-keyed matrix
(bug 0a), TF-IDF values instead of counts (bug 0c), and before --count-mode
cells (#86) existed.  v2 uses the Stage-1b acceptance run (clip-seeded v2,
mouse1, verified keying, layers['counts']), germ-cell STAGE labels derived from
the orthogonal STARsolo GEX matrix (step4_stage.py marker-argmax), restricted to
the 1,294 STARsolo cells and the robust SPC/RS/ES ordering (SPG dropped, n=64;
manuscript/11) -> 1,230 cells x 76,711 PAS, 3 stage pairs per run.

MAIN ARMS (TRUE + 20 label permutations, seeds 1..20, identical permutations
across arms, full `switch diff` pipeline per run incl. the default top-200
wilcoxon marker pre-selection, so label double-dipping is inside the null):
  A  fisher --count-mode reads   (legacy, pseudoreplicated)
  B  fisher --count-mode cells   (D4, de-pseudoreplicated)
  C  nb_pairwise                 (per-cell NB GLM, Wald)
DIAGNOSTIC ARMS (same, but --marker-top-n 0: marker pre-selection OFF) isolate
the test's own calibration from the marker double-dip:
  A0 / B0 / C0 (C0 only if its runs exist -- NB on all 76,711 PAS is slow)

Metrics per arm: % null p<0.05, % null q<0.05, null runs with >=1 q<0.05 hit,
mean hits/null run vs TRUE, QQ + KS vs uniform; PLUS a permutation-calibrated
alternative: for each TRUE test the empirical p is its rank in the pooled null
p-values of the same stage pair (all completed perms), BH per pair -> q_perm.

The script reads whatever runs are complete (DONE.ok) and says so on the figure.

Outputs
-------
manuscript/figures/fig4_calibration.{png,pdf}           600 dpi / fonttype 42
manuscript/figures/fig4_calibration.caption.md          sidecar ('## Legend' = the journal
                                                        legend, single source of the caption;
                                                        '## Provenance'; script-written)
results/figures/manuscript/fig4_calibration.png
results/figures/manuscript/fig4_calibration.tsv         per-run summary
results/figures/manuscript/fig4_calibration_stats.tsv   headline stats per arm
results/figures/manuscript/fig4_calibration_per_pair.tsv
results/figures/manuscript/fig4_calibration_null_pvalues.tsv
results/figures/manuscript/fig4_calibration_hist.tsv
results/figures/manuscript/fig4_calibration_true_permcal.tsv  TRUE tests + q_perm
results/fdr_calibration_v2/fig4_calibration_{null_pvalues,true_permcal}_all_arms.tsv.gz
results/fdr_calibration_v2/report.json                  machine-readable

The INPUT run directory stays results/fdr_calibration_v2/ (it names the
calibration run, not the figure stem, and is not renamed).

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/fig4_calibration.py
"""
import glob, json, os, re, sys, textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import offset_copy
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, kstest, false_discovery_control

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pubstyle import PAL, TYPE, apply_rc, sentence_case   # shared publication style

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

WD = "/mnt/ssd1/Projects/PeakATail_wd"
V2 = f"{WD}/results/fdr_calibration_v2"
FIGDIR = f"{WD}/manuscript/figures"
TSVDIR = f"{WD}/results/figures/manuscript"
os.makedirs(FIGDIR, exist_ok=True); os.makedirs(TSVDIR, exist_ok=True)
NAME = "fig4_calibration"
FDR = 0.05
N_PERMS_PLANNED = 20
EFFECT_FLOOR = 0.10     # |delta_proportion| floor (manuscript/13 sec.2); nb: |log2FC| >= 1

# 2026-09-02 design pass (DESIGN_DIRECTIVES.md item 4, recorded exception).  The
# figure is a property of TEST CONFIGURATIONS a user chooses -- test x count unit x
# marker pre-selection -- not a development story.  Colour therefore follows the
# STATUS the panel is about (PAL['good'] = controls the FDR under the label-
# permutation null, PAL['bad'] = anti-conservative), which is the one legitimate
# use of the reserved status tokens; every configuration is also named in full on
# its own row, so nothing depends on hue alone.  Six-checks (validate_palette.js,
# 2026-09-02): good/bad adjacent CVD deltaE 11.0 deutan / 31.4 tritan,
# normal-vision 25.8, contrast on the surface >= 3:1 -- ALL PASS.
BLUE, GREEN, VERM = PAL["peakatail"], PAL["good"], PAL["bad"]
ARMS = {  # arm dir -> meta ; main arms first, diagnostics after
    "A_fisher_reads": dict(label="fisher, count-mode reads", short="fisher\nreads", prefix="fisher", color=BLUE, markers=True,
                           test="Fisher", unit="reads"),
    "B_fisher_cells": dict(label="fisher, count-mode cells", short="fisher\ncells", prefix="fisher", color=GREEN, markers=True,
                           test="Fisher", unit="cells"),
    "C_nb_pairwise": dict(label="nb_pairwise", short="nb\npairwise", prefix="nb_pairwise", color=VERM, markers=True,
                          test="NB pairwise", unit="cells"),
    "A0_fisher_reads_nomarker": dict(label="fisher reads, no marker pre-sel.", short="fisher\nreads\nno-mk", prefix="fisher", color=BLUE, markers=False,
                                     test="Fisher", unit="reads"),
    "B0_fisher_cells_nomarker": dict(label="fisher cells, no marker pre-sel.", short="fisher\ncells\nno-mk", prefix="fisher", color=GREEN, markers=False,
                                     test="Fisher", unit="cells"),
    "C0_nb_pairwise_nomarker": dict(label="nb_pairwise, no marker pre-sel.", short="nb\npairwise\nno-mk", prefix="nb_pairwise", color=VERM, markers=False,
                                    test="NB pairwise", unit="cells"),
}
INK, MUTED, GRID, SURFACE = PAL["ink"], PAL["muted"], PAL["grid"], "#FFFFFF"

# ----------------------------------------------------------------------------
# 1. harvest
# ----------------------------------------------------------------------------
def read_run(diff_dir, prefix):
    frames = []
    for f in sorted(glob.glob(os.path.join(diff_dir, f"{prefix}_*_vs_*.tsv"))):
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        pair = re.sub(rf"^{prefix}_|\.tsv$", "", os.path.basename(f))
        frames.append(pd.DataFrame({
            "pair": pair, "pas_id": df["pas_id"].astype(str),
            "gene_id": df["gene_id"].astype(str) if "gene_id" in df else "",
            "pvalue": df["pvalue"].astype(float), "qvalue": df["qvalue"].astype(float),
            "effect": (df["delta_proportion"].astype(float) if "delta_proportion" in df
                       else df["log2fc"].astype(float)),
            "effect_kind": "delta_proportion" if "delta_proportion" in df else "log2fc",
            "dispersion": df["dispersion"].astype(float) if "dispersion" in df else np.nan,
        }))
    if not frames:
        return pd.DataFrame(columns=["pair", "pas_id", "gene_id", "pvalue", "qvalue", "effect", "effect_kind", "dispersion"])
    return pd.concat(frames, ignore_index=True)

def elapsed(d):
    try:
        return json.loads(open(os.path.join(d, "DONE.ok")).read())["elapsed_s"]
    except Exception:
        return np.nan

data = {}
for arm, meta in ARMS.items():
    tdir = f"{V2}/{arm}/true"
    if not os.path.exists(f"{tdir}/DONE.ok"):
        if meta["markers"]:
            raise SystemExit(f"{arm}: TRUE run not complete (main arm)")
        print(f"{arm}: not run -- skipped (diagnostic arm)"); continue
    true_df = read_run(f"{tdir}/differential", meta["prefix"])
    assert len(true_df), f"{arm}: TRUE differential TSVs empty"
    perm_frames, perm_ids, perm_el = [], [], []
    for pdir in sorted(glob.glob(f"{V2}/{arm}/null/perm_*")):
        if not os.path.exists(f"{pdir}/DONE.ok"):
            continue
        k = int(os.path.basename(pdir).split("_")[1])
        df = read_run(f"{pdir}/differential", meta["prefix"])
        df.insert(0, "perm", k)
        perm_frames.append(df); perm_ids.append(k); perm_el.append(elapsed(pdir))
    if not perm_ids:
        if meta["markers"]:
            raise SystemExit(f"{arm}: no completed permutations (main arm)")
        print(f"{arm}: TRUE only, no permutations yet -- skipped"); continue
    data[arm] = dict(true_df=true_df, null_df=pd.concat(perm_frames, ignore_index=True),
                     perm_ids=perm_ids, true_elapsed=elapsed(tdir),
                     perm_elapsed_mean=float(np.nanmean(perm_el)), **meta)
    print(f"{arm}: TRUE {len(true_df)} tests; {len(perm_ids)} perms, {len(data[arm]['null_df'])} null tests")
main_arms = [a for a in data if data[a]["markers"]]
diag_arms = [a for a in data if not data[a]["markers"]]
info = json.load(open(f"{V2}/input/input_info.json"))

# ----------------------------------------------------------------------------
# 2. per-run summary, per-pair, pooled null
# ----------------------------------------------------------------------------
rows = []
def summarize(arm, tag, seed, df):
    by_pair = df.groupby("pair")["qvalue"].apply(lambda q: int((q < FDR).sum()))
    rows.append(dict(
        arm=arm, run=tag, seed=seed, n_tests=len(df), n_pairs_tested=df["pair"].nunique(),
        n_q_lt_fdr=int((df["qvalue"] < FDR).sum()),
        n_pair_families_with_hit=int((by_pair > 0).sum()),
        frac_p_lt_05=float((df["pvalue"] < 0.05).mean()) if len(df) else np.nan,
        frac_q_lt_fdr=float((df["qvalue"] < FDR).mean()) if len(df) else np.nan,
        min_pvalue=float(df["pvalue"].min()) if len(df) else np.nan,
        median_pvalue=float(df["pvalue"].median()) if len(df) else np.nan,
        **{f"hits_{p}": int(v) for p, v in by_pair.items()},
    ))
for arm, d in data.items():
    summarize(arm, "true", "", d["true_df"])
    for k in d["perm_ids"]:
        summarize(arm, f"perm_{k:02d}", k, d["null_df"][d["null_df"]["perm"] == k])
run_tab = pd.DataFrame(rows)
run_tab.to_csv(f"{TSVDIR}/{NAME}.tsv", sep="\t", index=False)

# pooled null p-values: main arms (small) next to the figure; ALL arms (the
# no-marker arms carry ~2M null tests each) gzipped in the run directory.
pooled_cols = ["arm", "perm", "pair", "pas_id", "gene_id", "pvalue", "qvalue", "effect"]
pooled = pd.concat([d["null_df"].assign(arm=arm)[pooled_cols] for arm, d in data.items() if d["markers"]], ignore_index=True)
pooled.to_csv(f"{TSVDIR}/{NAME}_null_pvalues.tsv", sep="\t", index=False)
pd.concat([d["null_df"].assign(arm=arm)[pooled_cols] for arm, d in data.items()], ignore_index=True).to_csv(
    f"{V2}/{NAME}_null_pvalues_all_arms.tsv.gz", sep="\t", index=False, compression="gzip")

pair_rows = []
for arm, d in data.items():
    for pair, g in d["null_df"].groupby("pair"):
        t = d["true_df"][d["true_df"]["pair"] == pair]
        pair_rows.append(dict(arm=arm, pair=pair, n_null_tests=len(g),
                              null_frac_p_lt_05=float((g["pvalue"] < 0.05).mean()),
                              null_frac_q_lt_fdr=float((g["qvalue"] < FDR).mean()),
                              null_mean_hits_per_run=float(g.groupby("perm")["qvalue"].apply(lambda q: (q < FDR).sum()).mean()),
                              true_n_tests=len(t), true_q_hits=int((t["qvalue"] < FDR).sum())))
pd.DataFrame(pair_rows).to_csv(f"{TSVDIR}/{NAME}_per_pair.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 3. permutation-calibrated q for the TRUE run
# ----------------------------------------------------------------------------
# For each TRUE test (pair, PAS): p_emp = (1 + #{null p of the same pair <= p}) /
# (1 + N_null_pair), null pooled over the completed permutations (a global-null
# reference distribution, SAM-style: null tests are on whatever PAS the
# permuted-label pipeline tested -- the same pipeline the TRUE run went
# through).  q_perm = BH on p_emp within each pair (mirrors the tool's per-pair
# BH family).  Resolution floor: 1/(N_null_pair+1).
permcal_frames = []
for arm, d in data.items():
    t = d["true_df"].copy()
    t["p_emp"] = np.nan; t["n_null_ref"] = 0
    for pair, g in t.groupby("pair"):
        null_p = np.sort(d["null_df"].loc[d["null_df"]["pair"] == pair, "pvalue"].values)
        le = np.searchsorted(null_p, g["pvalue"].values, side="right")
        t.loc[g.index, "p_emp"] = (1.0 + le) / (1.0 + len(null_p))
        t.loc[g.index, "n_null_ref"] = len(null_p)
    t["q_perm"] = np.nan
    for pair, g in t.groupby("pair"):
        t.loc[g.index, "q_perm"] = false_discovery_control(g["p_emp"].values, method="bh")
    t["nominal_hit"] = t["qvalue"] < FDR
    t["perm_hit"] = t["q_perm"] < FDR
    is_dp = (t["effect_kind"] == "delta_proportion").all()
    t["effect_floor_pass"] = t["effect"].abs() >= (EFFECT_FLOOR if is_dp else 1.0)
    t.insert(0, "arm", arm)
    d["permcal"] = t
    permcal_frames.append(t)
permcal = pd.concat(permcal_frames, ignore_index=True)
permcal[permcal["arm"].isin(main_arms)].to_csv(f"{TSVDIR}/{NAME}_true_permcal.tsv", sep="\t", index=False)
permcal.to_csv(f"{V2}/{NAME}_true_permcal_all_arms.tsv.gz", sep="\t", index=False, compression="gzip")

# ----------------------------------------------------------------------------
# 4. headline stats per arm
# ----------------------------------------------------------------------------
def calibrated(st):
    """Operational rule for 'FDR-calibrated' = not anti-conservative.  One-sided
    on purpose: a discrete test (Fisher) is conservative (mass at p=1) and fails
    a two-sided KS-vs-uniform trivially without ever inflating the FDR; KS is
    reported as a diagnostic, not gated."""
    return (st["frac_null_p_lt_05"] <= 0.07 and st["frac_null_p_lt_01"] <= 0.015
            and st["frac_null_q_lt_fdr"] <= 0.05 and st["n_runs_with_q_hit"] <= 0.25 * st["n_perms"])
def conservative(st):
    """Verifier (2026-08-21): an exactly calibrated test with 3 BH families per run
    would show >=1 q<0.05 hit in roughly 1-3 of 20 null runs; 0/20 (0/60 families)
    together with a null p<0.05 rate at the lower edge (~3%) is conservative, not
    exactly calibrated.  Label it so."""
    return st["frac_null_p_lt_05"] < 0.035 or st["n_runs_with_q_hit"] == 0
PREREG_RULE = ("calibrated iff null p<0.05 rate <= 7% AND null p<0.01 rate <= 1.5% AND null q<0.05 rate <= 5% AND "
               "<= 25% of null runs carry any q<0.05 hit (anti-conservative checks only; KS reported, not gated)")

stats = {}
for arm, d in data.items():
    nt = run_tab[(run_tab.arm == arm) & (run_tab.run != "true")]
    tt = run_tab[(run_tab.arm == arm) & (run_tab.run == "true")].iloc[0]
    p = d["null_df"]["pvalue"].values; q = d["null_df"]["qvalue"].values
    ks = kstest(p, "uniform"); pc = d["permcal"]
    n_true_hits = int(tt["n_q_lt_fdr"]); mean_null_hits = float(nt["n_q_lt_fdr"].mean())
    st = dict(
        arm=arm, label=d["label"], marker_preselection=("top-200 wilcoxon" if d["markers"] else "off"),
        n_perms=len(d["perm_ids"]), n_null_tests=int(len(p)), null_tests_per_run=float(nt["n_tests"].mean()),
        frac_null_p_lt_05=float((p < 0.05).mean()), frac_null_p_lt_01=float((p < 0.01).mean()),
        frac_null_q_lt_fdr=float((q < FDR).mean()),
        n_runs_with_q_hit=int((nt["n_q_lt_fdr"] > 0).sum()),
        mean_hits_per_null_run=mean_null_hits, median_hits_per_null_run=float(nt["n_q_lt_fdr"].median()),
        max_hits_per_null_run=int(nt["n_q_lt_fdr"].max()),
        min_null_p=float(p.min()), ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
        n_fam=int(nt["n_pairs_tested"].sum()), n_fam_hit=int(nt["n_pair_families_with_hit"].sum()),
        true_n_tests=int(tt["n_tests"]), true_q_hits=n_true_hits,
        true_frac_p_lt_05=float(tt["frac_p_lt_05"]), true_frac_q_lt_fdr=float(tt["frac_q_lt_fdr"]),
        empirical_fdr_at_nominal=(mean_null_hits / n_true_hits) if n_true_hits else np.nan,
        expected_false_in_true_at_nominal=float((q < FDR).mean()) * int(tt["n_tests"]),
        true_perm_q_hits=int(pc["perm_hit"].sum()),
        true_nominal_hits_surviving_perm_q=int((pc["nominal_hit"] & pc["perm_hit"]).sum()),
        true_perm_q_hits_with_effect_floor=int((pc["perm_hit"] & pc["effect_floor_pass"]).sum()),
        true_nominal_hits_with_effect_floor=int((pc["nominal_hit"] & pc["effect_floor_pass"]).sum()),
        effect_kind=str(pc["effect_kind"].iloc[0]),
        effect_floor=(EFFECT_FLOOR if pc["effect_kind"].iloc[0] == "delta_proportion" else 1.0),
        min_p_emp=float(pc["p_emp"].min()), true_elapsed_s=d["true_elapsed"], perm_elapsed_mean_s=d["perm_elapsed_mean"],
    )
    # nb only: share of tests whose plug-in dispersion sits at the 1e-4 clip floor
    # (near-Poisson variance) -- among all null tests vs among null q<0.05 hits
    disp = d["null_df"]["dispersion"]
    if disp.notna().any():
        at_floor = disp <= 1e-4 + 1e-9
        st["null_share_dispersion_floor"] = float(at_floor.mean())
        hits = d["null_df"]["qvalue"] < FDR
        st["null_hits_share_dispersion_floor"] = float(at_floor[hits].mean()) if hits.any() else np.nan
    else:
        st["null_share_dispersion_floor"] = np.nan; st["null_hits_share_dispersion_floor"] = np.nan
    st["calibrated_by_prereg_rule"] = bool(calibrated(st))
    st["conservative"] = bool(conservative(st))
    stats[arm] = st
stat_tab = pd.DataFrame(stats.values())
stat_tab.to_csv(f"{TSVDIR}/{NAME}_stats.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 5. NOT PLOTTED SINCE 2026-09-02 -- the null p-value histogram and the QQ-vs-
#    uniform arrays.  The design pass (DESIGN_DIRECTIVES.md items 1 + 4) cut the
#    figure to the three panels that carry the message; the per-arm distribution
#    shape is figS9's job.  The TSV is still written for the record and every
#    value it used to plot is asserted in the sidecar-generation check below.
# ----------------------------------------------------------------------------
nbins = 20; hist_rows = []
for arm, d in data.items():
    p = d["null_df"]["pvalue"].values
    cnt, edges = np.histogram(p, bins=nbins, range=(0, 1))
    dens = cnt / len(p) * nbins
    for lo, hi, c, dv in zip(edges[:-1], edges[1:], cnt, dens):
        hist_rows.append(dict(arm=arm, bin_lo=lo, bin_hi=hi, count=int(c), density=dv))
HIST = pd.DataFrame(hist_rows)
HIST.to_csv(f"{TSVDIR}/{NAME}_hist.tsv", sep="\t", index=False)   # not plotted since 2026-09-02
QQ_FLOOR = 1e-16
QQ_CLIPPED = {arm: int((d["null_df"]["pvalue"].values < QQ_FLOOR).sum()) for arm, d in data.items()}

# ----------------------------------------------------------------------------
# 6. figure -- "which test configurations control the FDR"
#    Three panels on ONE shared row axis of the six test configurations
#    (test x count unit x marker pre-selection).  One message per panel:
#      a  the test-level false-positive rate, against the honest-null 5%
#      b  how many false calls a single null run produces
#      c  what that costs in the real comparison
# ----------------------------------------------------------------------------
apply_rc()
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 6,
})
ANN, ANN_MIN = TYPE["annotation"], TYPE["annotation_min"]
short1 = lambda d: d["short"].replace("\n", " ")

# rows: worst-controlled at the bottom, so the eye runs down the failure gradient
ROWS = sorted(data, key=lambda a: stats[a]["frac_null_p_lt_05"])
NROW = len(ROWS)
YS = {a: NROW - 1 - i for i, a in enumerate(ROWS)}          # top row = best controlled
CTRL = {a: bool(stats[a]["calibrated_by_prereg_rule"]) for a in ROWS}
COL = {a: (GREEN if CTRL[a] else VERM) for a in ROWS}
ROWLAB = [sentence_case(f"{data[a]['test']} · {data[a]['unit']}") + "\n"
          + sentence_case("markers: off" if not data[a]["markers"] else "markers: top-200")
          for a in ROWS]

fig = plt.figure(figsize=(7.09, 3.32))
gs = fig.add_gridspec(1, 3, width_ratios=[1.00, 1.12, 1.00], wspace=0.26,
                      left=0.160, right=0.985, top=0.855, bottom=0.150)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
ax_c = fig.add_subplot(gs[0, 2], sharey=ax_a)

def panel_tag(ax, letter, text, dy=1.035):
    """House panel tag (shared with Fig 2/3/5/6): a bold panel letter at
    TYPE['panel_letter'] and, offset a fixed 13 pt to its right so the gap does
    not scale with panel width, the sentence-cased title at TYPE['panel_title']."""
    ax.text(0.0, dy, letter, transform=ax.transAxes, fontsize=TYPE["panel_letter"],
            fontweight="bold", va="bottom", ha="left", color=INK)
    ax.text(0.0, dy, sentence_case(text), va="bottom", ha="left", color=INK,
            fontsize=TYPE["panel_title"],
            transform=offset_copy(ax.transAxes, fig=ax.figure, x=13.0, y=0.0, units="points"))


def style(ax, first=False):
    ax.set_axisbelow(True); ax.grid(True, axis="x", color=GRID, lw=0.5)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", length=2.5, pad=2)
    if not first:
        plt.setp(ax.get_yticklabels(), visible=False)

# --- a: test-level false-positive rate -------------------------------------
ax_a.axvline(100 * FDR, color=INK, lw=0.8, zorder=2)
for a in ROWS:
    y, v = YS[a], 100 * stats[a]["frac_null_p_lt_05"]
    ax_a.plot([100 * FDR, v], [y, y], color=COL[a], lw=1.4, solid_capstyle="butt", zorder=3)
    ax_a.plot([v], [y], "o", ms=5.0, color=COL[a], mec=SURFACE, mew=0.8, zorder=4)
    if v >= 100 * FDR:                       # label outside the stem end
        ax_a.text(v + 1.0, y, f"{v:.1f}%", ha="left", va="center", fontsize=ANN, color=COL[a], zorder=5)
    else:                                    # below the reference: label sits over the dot
        ax_a.text(v, y + 0.24, f"{v:.1f}%", ha="center", va="bottom", fontsize=ANN, color=COL[a], zorder=5)
ax_a.set_xlim(0, 29.5); ax_a.set_xticks([0, 5, 10, 15, 20, 25])
ax_a.set_ylim(-0.72, NROW + 0.60)
ax_a.set_yticks(range(NROW)); ax_a.set_yticklabels(ROWLAB[::-1], linespacing=1.35)
ax_a.set_xlabel(sentence_case(f"null tests with p < {FDR:g}  (%)"))
ax_a.set_ylabel(sentence_case("test configuration"), labelpad=6)
style(ax_a, first=True)
ax_a.text(100 * FDR + 0.7, NROW + 0.48, f"{100 * FDR:.0f}% expected under\nan honest null",
          ha="left", va="top", fontsize=ANN, color=INK, linespacing=1.3, zorder=5)
ax_a.legend(handles=[Line2D([], [], ls="none", marker="o", ms=5.0, color=GREEN, label=sentence_case("controls the FDR")),
                     Line2D([], [], ls="none", marker="o", ms=5.0, color=VERM, label=sentence_case("anti-conservative"))],
            loc="lower right", bbox_to_anchor=(1.02, -0.035), frameon=False,
            handletextpad=0.35, labelspacing=0.28, borderaxespad=0.0)
panel_tag(ax_a, "a", "False-positive rate")

# --- b: false calls per null run -------------------------------------------
rng = np.random.default_rng(0)
bmax = max(int(run_tab[(run_tab.arm == a) & (run_tab.run != "true")]["n_q_lt_fdr"].max()) for a in ROWS)
for a in ROWS:
    y = YS[a]
    nt = run_tab[(run_tab.arm == a) & (run_tab.run != "true")]["n_q_lt_fdr"].to_numpy()
    ax_b.plot(nt, y + rng.uniform(-0.17, 0.17, nt.size), "o", ms=3.4, ls="none",
              color=COL[a], mec=SURFACE, mew=0.5, alpha=0.85, zorder=3)
    ax_b.plot([np.median(nt)] * 2, [y - 0.30, y + 0.30], color=INK, lw=1.0, zorder=4)
ax_b.set_xscale("symlog", linthresh=1.0, linscale=0.45)
ax_b.set_xlim(-0.35, bmax * 4.0)
ax_b.set_xticks([0, 1, 10, 100, 1000]); ax_b.set_xticklabels(["0", "1", "10", "100", "1,000"])
ax_b.set_xlabel(sentence_case(f"PAS called q < {FDR:g} in one null run"))
style(ax_b)
_best = ROWS[0]
ax_b.annotate(sentence_case(f"no false call in any of the {stats[_best]['n_perms']} runs"),
              xy=(0, YS[_best]), xytext=(1.6, YS[_best] + 0.62), fontsize=ANN, color=GREEN,
              ha="left", va="bottom", zorder=6,
              arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.7, shrinkA=0, shrinkB=3))
panel_tag(ax_b, "b", "False calls per null run")

# --- c: consequence for the real comparison ---------------------------------
for a in ROWS:
    y, st = YS[a], stats[a]
    n0 = max(st["true_q_hits"], 1)
    v1 = 100.0 * st["true_nominal_hits_surviving_perm_q"] / n0
    v2 = 100.0 * st["true_perm_q_hits_with_effect_floor"] / n0
    ax_c.barh(y, v1, height=0.50, color=COL[a], edgecolor="none", zorder=3)
    ax_c.plot([v2], [y], marker="D", ms=4.2, mfc=SURFACE, mec=INK, mew=0.9, zorder=5)
    ax_c.text(v1 + 2.5, y, "100%" if v1 >= 99.95 else f"{v1:.1f}%",
              ha="left", va="center", fontsize=ANN, color=COL[a], zorder=6)
ax_c.set_xlim(0, 120); ax_c.set_xticks([0, 25, 50, 75, 100])
ax_c.set_xlabel(sentence_case(f"real-run q < {FDR:g} calls kept after\npermutation calibration  (%)"))
style(ax_c)
ax_c.text(0.0, NROW + 0.50, "◇ " + sentence_case("also clears the effect floor"),
          ha="left", va="top", fontsize=ANN, color=INK, zorder=6)
panel_tag(ax_c, "c", "Cost in the real run")

# --- title + caveats --------------------------------------------------------
calib = [a for a in stats if stats[a]["calibrated_by_prereg_rule"]]
notcal = [a for a in stats if a not in calib]
if not calib:
    verdict = "no arm is FDR-calibrated at face value"
else:
    verdict = ("valid FDR control: " + " / ".join(short1(data[a]) + (" (conservative)" if stats[a]["conservative"] else " (calibrated)") for a in calib)
               + "; anti-conservative: " + " / ".join(short1(data[a]) for a in notcal))
verdict_sentence = (
    # (unreachable while any configuration passes the rule; kept so the figure can
    #  still state the all-fail case without dev-history framing -- directive 4)
    ("No test configuration is FDR-calibrated under the label-permutation null; permutation-calibrated q-values are required for any switch call a user reports."
     if not calib else
     "Only " + " and ".join(data[a]["label"] for a in calib) + " gives valid FDR control under the label-permutation null"
     + (" (conservative: " + ", ".join(f"{stats[a]['frac_null_p_lt_05']:.1%} null p<0.05, {stats[a]['frac_null_q_lt_fdr']:.2%} null q<0.05" for a in calib) + ")")
     + "; every arm that keeps the default top-200 marker pre-selection is anti-conservative ("
     + ", ".join(f"{short1(data[a])} {stats[a]['frac_null_p_lt_05']:.0%} null p<0.05" for a in notcal if data[a]["markers"]) + ")"
     + (", as is " + ", ".join(short1(data[a]) for a in notcal if not data[a]["markers"]) + " without pre-selection" if any(not data[a]["markers"] for a in notcal) else "")
     + "; permutation-calibrated q is NOT required for the calibrated arm but IS required whenever marker pre-selection or an anti-conservative test is used."))
# (the suptitle -- the verdict line -- moved to the sidecar Legend opener,
#  journal style: 'Figure 4 | <title>.'; the verdict text is unchanged there)
prelim = "" if all(stats[a]["n_perms"] == N_PERMS_PLANNED for a in main_arms) else " (PRELIMINARY: main arms incomplete)"
diag_note = ("; ".join(f"{data[a]['label']} {stats[a]['n_perms']}/{N_PERMS_PLANNED} perms" for a in diag_arms)
             if diag_arms else "no diagnostic (no-marker) arms complete")
caveat = (
    f"Input: Stage-1b acceptance run clusters.h5ad (verified keying, layers['counts'], Stage 0c); labels = STARsolo-GEX marker-argmax stage "
    f"(step4_stage.py, orthogonal to the PAS matrix), restricted to the 1,294 STARsolo cells, SPG dropped (n={info['n_spg_dropped']}, fragile; manuscript/11) "
    f"-> {info['n_cells']:,} cells (SPC {info['stage_counts']['SPC']}, RS {info['stage_counts']['RS']}, ES {info['stage_counts']['ES']}), "
    f"{info['n_pas']:,} PAS; 3-stage agreement with PAS-derived labels {info['agreement_with_pas_labels_3stage']:.1%}. "
    f"Null = obs['stage'] permuted across cells (label counts kept; numpy seeds 1-{N_PERMS_PLANNED}, the SAME permutations for every arm), full `switch diff` "
    f"re-run per permutation; main arms include the default top-200 wilcoxon marker pre-selection on .X (label double-dip deliberately inside the null), "
    f"the configurations with marker pre-selection off switch it off (--marker-top-n 0; {diag_note}); BH families are per stage pair. "
    + "; ".join(f"{data[a]['label']}: {stats[a]['frac_null_p_lt_05']:.1%} null p<0.05, {stats[a]['frac_null_q_lt_fdr']:.2%} null q<{FDR}, "
                f"{stats[a]['n_runs_with_q_hit']}/{stats[a]['n_perms']} null runs with hits" for a in stats) + f"{prelim}. "
    f"TRUE vs null hit COUNTS are not directly comparable (permuted labels change the marker set / the >=10-cells filter); per-test null rates are the calibrated quantities. "
    f"Permutation-calibrated q: empirical p of each TRUE test = rank within the pooled null p-values of the same pair (global-null reference, resolution 1/(N_null+1)), BH per pair; "
    f"effect floor |delta proportion| >= {EFFECT_FLOOR} (fisher) or |log2FC| >= 1 (nb). Rule for 'calibrated': {PREREG_RULE}. "
    f"Caveats: one mouse, one tissue with very large true stage effects (TRUE hit counts are upper bounds, not precision); label permutation tests "
    f"exchangeability under the global null of no stage effect -- it does not model within-stage heterogeneity; nb_pairwise dispersion is a per-PAS plug-in "
    f"clipped to [1e-4, 10] without shrinkage -- "
    + "; ".join(f"{short1(data[a])}: {stats[a]['null_share_dispersion_floor']:.1%} of null tests but {stats[a]['null_hits_share_dispersion_floor']:.0%} of null q<{FDR} hits sit at the floor"
                for a in stats if not np.isnan(stats[a].get("null_share_dispersion_floor", np.nan)))
    + "; fisher p-values are discrete (mass at p=1). Palette six-checks re-validated 2026-09-02."
)
# (the caveat block that hung below the panels moved to the sidecar Legend,
#  verbatim -- the `caveat` string is written there; nothing is drawn here)

# ----------------------------------------------------------------------------
# 6b. bounding-box discipline -- no text may overlap another text artist, no
# annotation may fall below the minimum print size, and nothing may run off the
# canvas (DESIGN_DIRECTIVES.md item 1).  The audit runs BEFORE the save, and the
# save is deliberately NOT `bbox_inches="tight"`: the canvas is laid out at the
# final print width (7.09 in / 180 mm) and must stay there, so an artist that
# overflows has to be fixed rather than silently grow the page.
# ----------------------------------------------------------------------------
def bbox_audit(figure, min_pt=ANN_MIN, name=NAME, margin=1.5):
    figure.canvas.draw()
    rend = figure.canvas.get_renderer()
    items, small = [], []

    def in_view(ax_, axis):
        """Tick labels for ticks outside the current view are laid out but never
        drawn; keep only the ones the reader actually sees."""
        lo, hi = (ax_.get_xlim() if axis == "x" else ax_.get_ylim())
        lo, hi = min(lo, hi), max(lo, hi)
        locs = ax_.get_xticks() if axis == "x" else ax_.get_yticks()
        labs = ax_.get_xticklabels() if axis == "x" else ax_.get_yticklabels()
        return [t for loc, t in zip(locs, labs) if lo - 1e-9 <= loc <= hi + 1e-9]

    for ax in figure.axes:
        arts = list(ax.texts) + in_view(ax, "x") + in_view(ax, "y")
        arts += [ax.xaxis.label, ax.yaxis.label, ax.title]
        lg = ax.get_legend()
        if lg is not None:
            arts += list(lg.get_texts())
        for t in arts:
            if not t.get_visible() or not t.get_text().strip():
                continue
            if t.get_fontsize() < min_pt - 1e-9:
                small.append((t.get_text()[:34], t.get_fontsize()))
            items.append((t.get_text()[:34], t.get_window_extent(rend)))
    for t in figure.texts:
        if t.get_visible() and t.get_text().strip():
            if t.get_fontsize() < min_pt - 1e-9:
                small.append((t.get_text()[:34], t.get_fontsize()))
            items.append((t.get_text()[:34], t.get_window_extent(rend)))
    assert not small, f"{name}: annotation below {min_pt} pt -- {small}"
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            ov_w = min(a.x1, b.x1) - max(a.x0, b.x0)
            ov_h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ov_w > 1.0 and ov_h > 1.0:
                bad.append((items[i][0], items[j][0], round(ov_w, 1), round(ov_h, 1)))
    assert not bad, f"{name}: overlapping text -- {bad}"
    fb = figure.bbox
    out = [(t, [round(v, 1) for v in (bb.x0, bb.y0, bb.x1, bb.y1)]) for t, bb in items
           if bb.x0 < fb.x0 + margin or bb.y0 < fb.y0 + margin
           or bb.x1 > fb.x1 - margin or bb.y1 > fb.y1 - margin]
    assert not out, (f"{name}: text runs off the {fb.x1 / figure.dpi:.2f} x "
                     f"{fb.y1 / figure.dpi:.2f} in canvas -- {out}")
    print(f"bbox audit: {len(items)} text artists, no overlap, none below {min_pt} pt, "
          f"none off-canvas")

bbox_audit(fig)

for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = f"{FIGDIR}/{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
fig.savefig(f"{TSVDIR}/{NAME}.png", dpi=600)

# ----------------------------------------------------------------------------
# 6c. SIDECAR-GENERATION CHECKS -- the value-asserts for everything the design
# pass took OFF the image.  The numbers must stay pinned somewhere reproducible,
# so they are asserted here and then written into the Legend below.
# ----------------------------------------------------------------------------
assert len(HIST) == len(data) * nbins, (len(HIST), len(data), nbins)      # ex-panel a
for arm in data:
    _h = HIST[HIST.arm == arm]
    assert abs(_h["density"].mean() - 1.0) < 1e-9, arm                    # density normalisation
    assert int(_h["count"].sum()) == len(data[arm]["null_df"]), arm
    _first = float(_h.iloc[0]["density"])
    assert np.isfinite(_first) and _first > 0, arm
    assert np.isfinite(stats[arm]["ks_stat"]) and np.isfinite(stats[arm]["ks_p"]), arm   # ex-panel b
    assert 0 < stats[arm]["min_null_p"] <= 1, arm
    assert QQ_CLIPPED[arm] >= 0
HIST_LEAD = {a: float(HIST[(HIST.arm == a) & (HIST.bin_lo == 0.0)]["density"].iloc[0]) for a in data}
print("sidecar checks: hist density normalised, KS + min-null-p finite for all",
      len(data), "configurations")

# ----------------------------------------------------------------------------
# 6d. sidecar caption (added at the 2026-09-02 rename pass -- the figure
# convention requires the script, never a hand edit, to write the caption;
# this was the one asset without a sidecar).  Every number below is a live
# value computed above from the run TSVs, or a cited fact from the verified
# write-up manuscript/14_switch_calibration_v2.md (verdict SOUND).
# ----------------------------------------------------------------------------
_cfg_name = {a: f"{data[a]['test']} on {data[a]['unit']}, marker pre-selection "
                f"{'top-200 wilcoxon' if data[a]['markers'] else 'off'}" for a in data}
# compact form for the value lists below, where the full name would drown the numbers
_cfg_short = {a: f"{data[a]['test']}/{data[a]['unit']} markers "
                 f"{'top-200' if data[a]['markers'] else 'off'}" for a in data}
_arm_lines = "\n".join(
    f"- **{_cfg_name[a]}** (`{data[a]['label']}`): "
    f"{stats[a]['frac_null_p_lt_05']:.1%} null p<0.05, {stats[a]['frac_null_q_lt_fdr']:.2%} null q<{FDR}, "
    f"{stats[a]['n_runs_with_q_hit']}/{stats[a]['n_perms']} null runs with a q<{FDR} hit "
    f"(mean {stats[a]['mean_hits_per_null_run']:.0f} hits/null run); TRUE run {stats[a]['true_q_hits']:,} nominal hits, "
    f"{stats[a]['true_perm_q_hits']:,} after permutation-calibrated q, "
    f"{stats[a]['true_perm_q_hits_with_effect_floor']:,} after the effect floor -- "
    f"{'CALIBRATED' + (' (conservative)' if stats[a]['conservative'] else '') if stats[a]['calibrated_by_prereg_rule'] else 'ANTI-CONSERVATIVE'}"
    for a in ROWS)
cap_md = f"""# Fig 4 — `fig4_calibration` caption (generated by `scripts/manuscript_figures/fig4_calibration.py`; renamed 2026-09-02, see FIGURE_MAP.tsv)

## Legend

Figure 4 | **Which `switch diff` test configurations control the false-discovery rate.** Each row of the figure is one
configuration a user chooses on the command line — statistical test × count unit (`--count-mode reads|cells`) ×
marker pre-selection (`--marker-top-n 200`, the default, or `0`). Configurations are scored against a
label-permutation null run through the identical pipeline, so the score is a property of the configuration, not of
any one dataset. **{sentence_case(verdict_sentence.rstrip("."))}.**
Testis mouse1 clip-seeded v2, {info['n_cells']:,} STARsolo cells
(SPC {info['stage_counts']['SPC']} / RS {info['stage_counts']['RS']} / ES {info['stage_counts']['ES']}),
{info['n_pas']:,} PAS, 3 stage pairs per run, {N_PERMS_PLANNED} label permutations shared identically across all
configurations, the full `switch diff` pipeline re-run per permutation. Green = controls the FDR under the
pre-registered rule; vermilion = anti-conservative. Verdict in short form: {verdict}.

**Panel a** — share of null tests reaching p < {FDR:g} in each configuration, against the {100 * FDR:.0f}% an honest
null produces (vertical rule); the stem runs from that reference to the observed rate, so its length and side are the
size and sign of the miscalibration. **Panel b** — false calls made by a *single* null run: one dot per label
permutation ({N_PERMS_PLANNED} per configuration), bar = median, symmetric-log axis so exact zeros are drawn at 0.
**Panel c** — the cost in the real (unpermuted) comparison: the share of that configuration's nominal q < {FDR:g}
calls that survive permutation-calibrated q (bar; label = surviving count / all real-run calls), with the diamond
marking the share that also clears the pre-registered effect floor (|Δproportion| ≥ {EFFECT_FLOOR} for Fisher;
|log2FC| ≥ 1 for NB).

**Moved off the image at the 2026-09-02 design pass (no-loss rule), with the numbers kept here and in the audit
TSVs.** The per-configuration *null p-value distribution* is no longer drawn: its leading-bin density (0 ≤ p < 0.05,
uniform = 1) is {", ".join(f"{_cfg_short[a]} {HIST_LEAD[a]:.2f}" for a in ROWS)}, and the full 20-bin histogram is
`{NAME}_hist.tsv`. The *QQ-vs-uniform* panel is likewise gone: the smallest null p-value per configuration is
{", ".join(f"{_cfg_short[a]} {stats[a]['min_null_p']:.0e}" for a in ROWS)}; the number of null p-values at or below
the 1e-16 plotting floor was {", ".join(f"{_cfg_short[a]} {QQ_CLIPPED[a]:,}" for a in ROWS)}; and the
Kolmogorov–Smirnov statistic against uniform is
{", ".join(f"{_cfg_short[a]} D = {stats[a]['ks_stat']:.3f}" for a in ROWS)} — KS rejects for every configuration
because Fisher p-values are discrete, which is why it is reported and never gated. Per-configuration diagnostic
depth (per-stage-pair rates, expression strata, count-mode contrast, the NB dispersion floor, permutation-calibrated
vs nominal q) lives in **figS9**, and is not duplicated here.

**Per-configuration result** (rule for "controls the FDR": {PREREG_RULE}):
{_arm_lines}

{caveat}

**Caveats that travel with it (05_figure_index / 14, verified SOUND):** marker mode also changes the Fisher
gene denominator (restricted matrix), so marker-on vs marker-off are different tests, not subsets; single
mouse/tissue; TRUE hit counts are detectability upper bounds, not precision; the label-permutation null tests
the global null only; KS rejects for every arm (discrete Fisher) and is not gated.

## Provenance

Sources: `manuscript/14_switch_calibration_v2.md` (verified SOUND) + `results/fdr_calibration_v2/report.json`.
Every plotted value: `results/figures/manuscript/{NAME}.tsv` (per-run summary), `{NAME}_stats.tsv` (headline
stats per configuration), `{NAME}_per_pair.tsv`, `{NAME}_null_pvalues.tsv`, `{NAME}_true_permcal.tsv`.
Still written, **not plotted since 2026-09-02**: `{NAME}_hist.tsv` (the retired null p-value histogram); its
values, and the retired QQ panel's min-null-p / KS statistics, are asserted in the script's
sidecar-generation check and quoted in the Legend above.

Design: shared publication style `scripts/manuscript_figures/_pubstyle.py` (PAL / TYPE / `apply_rc()` /
`sentence_case()`), so Fig 1–6 read as one system. Colour follows the calibration STATUS
(`PAL['good']` = controls the FDR, `PAL['bad']` = anti-conservative) and every configuration is named in full on
its own row, so nothing depends on hue alone; six-checks re-validated 2026-09-02 (`validate_palette.js`, light
mode): lightness band PASS, chroma floor PASS, CVD separation ΔE 11.0 deutan / 31.4 tritan PASS, normal-vision
floor 25.8 PASS, contrast ≥ 3:1 PASS. Canvas 7.09 in (180 mm) at final print width; PNG 600 dpi, PDF vector with
subsetted TrueType (fonttype 42, no Type 3); the render passes an automated text-overlap, minimum-type-size
({TYPE['annotation_min']:.0f} pt floor) and off-canvas audit plus the 8-px edge check, and is saved at that fixed
canvas (no tight bounding box), so an overflowing artist fails the audit instead of quietly widening the figure.
"""
with open(f"{FIGDIR}/{NAME}.caption.md", "w") as fh:
    fh.write(cap_md)
print("wrote", f"{FIGDIR}/{NAME}.caption.md")

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
    assert im.shape[1] <= 600 * 7.10, (
        f"canvas {im.shape[1] / 600:.2f} in wide -- past the 180 mm (7.09 in) print target")
except ImportError:
    print("edge check skipped (no PIL)")

# ----------------------------------------------------------------------------
# 6. machine-readable report
# ----------------------------------------------------------------------------
report = dict(
    name=NAME, input=info, fdr=FDR, n_perms_planned=N_PERMS_PLANNED, main_arms=main_arms, diagnostic_arms=diag_arms,
    arms={a: {k: (v.item() if isinstance(v, (np.floating, np.integer, np.bool_)) else v) for k, v in stats[a].items()} for a in stats},
    verdict=verdict, verdict_sentence=verdict_sentence, prereg_rule=PREREG_RULE,
    outputs=dict(figure_png=f"{FIGDIR}/{NAME}.png", figure_pdf=f"{FIGDIR}/{NAME}.pdf",
                 caption_md=f"{FIGDIR}/{NAME}.caption.md",
                 run_tsv=f"{TSVDIR}/{NAME}.tsv", stats_tsv=f"{TSVDIR}/{NAME}_stats.tsv",
                 per_pair_tsv=f"{TSVDIR}/{NAME}_per_pair.tsv", null_pvalues_tsv=f"{TSVDIR}/{NAME}_null_pvalues.tsv",
                 hist_tsv=f"{TSVDIR}/{NAME}_hist.tsv", true_permcal_tsv=f"{TSVDIR}/{NAME}_true_permcal.tsv",
                 null_pvalues_all_arms_tsv_gz=f"{V2}/{NAME}_null_pvalues_all_arms.tsv.gz",
                 true_permcal_all_arms_tsv_gz=f"{V2}/{NAME}_true_permcal_all_arms.tsv.gz"),
)
json.dump(report, open(f"{V2}/report.json", "w"), indent=1, default=float)
cols = ["arm", "n_perms", "null_tests_per_run", "frac_null_p_lt_05", "frac_null_p_lt_01", "frac_null_q_lt_fdr", "n_runs_with_q_hit",
        "mean_hits_per_null_run", "ks_p", "true_n_tests", "true_q_hits", "true_perm_q_hits", "true_perm_q_hits_with_effect_floor", "calibrated_by_prereg_rule"]
print(stat_tab[cols].to_string(index=False))
print("VERDICT:", verdict)
print("VERDICT SENTENCE:", verdict_sentence)
