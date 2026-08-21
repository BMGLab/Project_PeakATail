#!/usr/bin/env python3
"""
fdr_calibration_v2.py -- re-calibration of `ema switch diff` on a CORRECT input.

Why v2: the first calibration (fdr_calibration.py) ran on a mis-keyed matrix
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
manuscript/figures/fdr_calibration_v2.{png,pdf}           300 dpi / fonttype 42
results/figures/manuscript/fdr_calibration_v2.png
results/figures/manuscript/fdr_calibration_v2.tsv         per-run summary
results/figures/manuscript/fdr_calibration_v2_stats.tsv   headline stats per arm
results/figures/manuscript/fdr_calibration_v2_per_pair.tsv
results/figures/manuscript/fdr_calibration_v2_null_pvalues.tsv
results/figures/manuscript/fdr_calibration_v2_hist.tsv
results/figures/manuscript/fdr_calibration_v2_true_permcal.tsv  TRUE tests + q_perm
results/fdr_calibration_v2/report.json                    machine-readable
"""
import glob, json, os, re, textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, kstest, false_discovery_control

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

WD = "/mnt/ssd1/Projects/PeakATail_wd"
V2 = f"{WD}/results/fdr_calibration_v2"
FIGDIR = f"{WD}/manuscript/figures"
TSVDIR = f"{WD}/results/figures/manuscript"
os.makedirs(FIGDIR, exist_ok=True); os.makedirs(TSVDIR, exist_ok=True)
NAME = "fdr_calibration_v2"
FDR = 0.05
N_PERMS_PLANNED = 20
EFFECT_FLOOR = 0.10     # |delta_proportion| floor (manuscript/13 sec.2); nb: |log2FC| >= 1

# validated palette (task-mandated; six-checks validated 2026-08-19: worst
# adjacent-pair CVD deltaE 11.0 deutan / 8.6 tritan, normal-vision floor 18.7,
# contrast on white >= 3:1 -- PASS).  Colour follows the TEST (entity), fixed
# order; "marker pre-selection off" is a secondary encoding (dashed / hollow).
BLUE, GREEN, VERM = "#0072B2", "#009E73", "#D55E00"
ARMS = {  # arm dir -> meta ; main arms first, diagnostics after
    "A_fisher_reads": dict(label="fisher, count-mode reads", short="fisher\nreads", prefix="fisher", color=BLUE, markers=True),
    "B_fisher_cells": dict(label="fisher, count-mode cells", short="fisher\ncells", prefix="fisher", color=GREEN, markers=True),
    "C_nb_pairwise": dict(label="nb_pairwise", short="nb\npairwise", prefix="nb_pairwise", color=VERM, markers=True),
    "A0_fisher_reads_nomarker": dict(label="fisher reads, no marker pre-sel.", short="fisher\nreads\nno-mk", prefix="fisher", color=BLUE, markers=False),
    "B0_fisher_cells_nomarker": dict(label="fisher cells, no marker pre-sel.", short="fisher\ncells\nno-mk", prefix="fisher", color=GREEN, markers=False),
    "C0_nb_pairwise_nomarker": dict(label="nb_pairwise, no marker pre-sel.", short="nb\npairwise\nno-mk", prefix="nb_pairwise", color=VERM, markers=False),
}
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

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
# 5. figure (2 x 2, double-column width)
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 7.6, "axes.labelsize": 7.4, "axes.titlesize": 8.0,
    "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6.1,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 5, "axes.linewidth": 0.7,
    "font.family": "DejaVu Sans", "hatch.linewidth": 0.6,
})
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.9))
ax_a, ax_b, ax_c, ax_d = axes.ravel()
def style(ax):
    ax.set_axisbelow(True); ax.grid(True, axis="y", color=GRID, lw=0.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
def ls_for(d):  # secondary encoding: marker pre-selection off = dashed/hollow
    return "-" if d["markers"] else (0, (2.2, 1.4))
short1 = lambda d: d["short"].replace("\n", " ")

# --- a: null p histograms ---------------------------------------------------
nbins = 20; hist_rows = []
for arm, d in data.items():
    p = d["null_df"]["pvalue"].values
    cnt, edges = np.histogram(p, bins=nbins, range=(0, 1))
    dens = cnt / len(p) * nbins
    ax_a.stairs(dens, edges, color=d["color"], lw=1.3 if d["markers"] else 1.0, ls=ls_for(d), zorder=4,
                label=f"{d['label']}: {stats[arm]['frac_null_p_lt_05']:.1%}")
    if d["markers"]:
        ax_a.stairs(dens, edges, fill=True, alpha=0.10, color=d["color"], zorder=3)
    for lo, hi, c, dv in zip(edges[:-1], edges[1:], cnt, dens):
        hist_rows.append(dict(arm=arm, bin_lo=lo, bin_hi=hi, count=int(c), density=dv))
pd.DataFrame(hist_rows).to_csv(f"{TSVDIR}/{NAME}_hist.tsv", sep="\t", index=False)
ax_a.axhline(1.0, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=5, label="uniform (honest null): 5.0%")
style(ax_a); ax_a.set_xlim(0, 1)
ymax_a = max(np.histogram(d["null_df"]["pvalue"].values, bins=nbins, range=(0, 1))[0].max() / len(d["null_df"]) * nbins for d in data.values())
ax_a.set_ylim(0, ymax_a * 1.55)
ax_a.set_xlabel("null p-value"); ax_a.set_ylabel("density (uniform = 1)")
ax_a.legend(loc="upper right", frameon=False, handlelength=1.8, title="share of null p < 0.05", title_fontsize=6.2, alignment="left")
ax_a.set_title("a  Null p-value distribution (label permutations)")

# --- b: QQ vs uniform -------------------------------------------------------
FLOOR = 1e-16
def qq_arrays(p):
    ps = np.sort(p); n = len(ps)
    obs = -np.log10(np.maximum(ps, FLOOR)); ii = np.arange(1, n + 1)
    exp = -np.log10((ii - 0.5) / n)
    lo = -np.log10(beta_dist.ppf(0.975, ii, n - ii + 1)); hi = -np.log10(beta_dist.ppf(0.025, ii, n - ii + 1))
    keep = np.unique(np.r_[np.arange(min(600, n)), np.linspace(0, n - 1, 1200).astype(int)])
    idx = np.sort((n - 1) - keep)
    return exp[idx], obs[idx], lo[idx], hi[idx], int((ps < FLOOR).sum())
lim = 0.0; clip_notes = []
for arm, d in data.items():
    exp, obs, lo, hi, n_clip = qq_arrays(d["null_df"]["pvalue"].values)
    ax_b.fill_between(exp, lo, hi, color=GRID, alpha=0.5, lw=0, zorder=2)
    if d["markers"]:
        ax_b.plot(exp, obs, ".", ms=2.6, color=d["color"], zorder=4, label=d["label"])
    else:
        ax_b.plot(exp, obs, "o", ms=2.6, mfc="none", mec=d["color"], mew=0.6, zorder=4, label=d["label"])
    lim = max(lim, exp.max())
    if n_clip:
        clip_notes.append(f"{n_clip} {short1(d)} p<1e-16 at top")
lim *= 1.05
ax_b.plot([0, lim], [0, lim], color=INK, lw=0.9, ls=(0, (3, 2)), zorder=3, label="y = x")
ax_b.legend(loc="upper left", frameon=False, handletextpad=0.5)
ax_b.text(0.98, 0.03, "gray: 95% pointwise uniform band" + ("\n" + "\n".join(clip_notes) if clip_notes else ""),
          transform=ax_b.transAxes, ha="right", va="bottom", fontsize=5.6, color=MUTED, linespacing=1.3)
style(ax_b); ax_b.set_xlim(0, lim)
ax_b.set_xlabel(r"expected $-\log_{10}p$ (uniform)"); ax_b.set_ylabel(r"observed null $-\log_{10}p$")
ax_b.set_title("b  QQ vs uniform (min null p: " + ", ".join(f"{short1(d)} {stats[a]['min_null_p']:.0e}" for a, d in data.items() if d["markers"]) + ")")

# --- c: q<0.05 hits per null run vs TRUE ------------------------------------
rng = np.random.default_rng(0)
arm_order = list(data.keys())
for xi, arm in enumerate(arm_order):
    d, st = data[arm], stats[arm]
    nt = run_tab[(run_tab.arm == arm) & (run_tab.run != "true")]
    xs = xi + rng.uniform(-0.16, 0.16, len(nt))
    if d["markers"]:
        ax_c.plot(xs, nt["n_q_lt_fdr"] + 0.5, "o", ms=3.6, mew=0.6, mec=SURFACE, color=d["color"], alpha=0.9, zorder=3, ls="none")
    else:
        ax_c.plot(xs, nt["n_q_lt_fdr"] + 0.5, "o", ms=3.6, mfc="none", mec=d["color"], mew=0.8, zorder=3, ls="none")
    ax_c.plot([xi - 0.3, xi + 0.3], [st["true_q_hits"] + 0.5] * 2, color=INK, lw=1.1, zorder=4)
    ax_c.plot(xi, st["true_q_hits"] + 0.5, "D", ms=4.2, color=INK, mec=SURFACE, mew=0.6, zorder=5)
    ax_c.text(xi, (st["true_q_hits"] + 0.5) * (1.9 if xi % 2 == 0 else 5.5), f"TRUE {st['true_q_hits']:,}", ha="center", va="bottom", fontsize=5.9, color=INK)
    tr = mtransforms.blended_transform_factory(ax_c.transData, ax_c.transAxes)
    ax_c.text(xi, 0.015, f"{st['n_runs_with_q_hit']}/{st['n_perms']} runs\nmean {st['mean_hits_per_null_run']:.0f}",
              transform=tr, ha="center", va="bottom", fontsize=5.5, color=INK)
ax_c.set_yscale("log"); ax_c.set_ylim(0.09, 8e6)
ax_c.set_yticks([0.5, 1.5, 10.5, 100.5, 1000.5, 10000.5]); ax_c.set_yticklabels(["0", "1", "10", "100", "1,000", "10,000"])
ax_c.set_xlim(-0.6, len(arm_order) - 0.4); ax_c.set_xticks(range(len(arm_order)))
ax_c.set_xticklabels([data[a]["short"] for a in arm_order])
style(ax_c); ax_c.set_ylabel(f"PAS with q < {FDR} per run (log; 0 at floor)")
ax_c.text(0.99, 0.015, "bottom text: null runs with >=1 hit / mean hits per null run", transform=ax_c.transAxes,
          ha="right", va="bottom", fontsize=5.2, color=MUTED) if False else None
h = [plt.Line2D([], [], marker="o", ls="none", ms=3.6, color=MUTED),
     plt.Line2D([], [], marker="o", ls="none", ms=3.6, mfc="none", mec=MUTED),
     plt.Line2D([], [], marker="D", ls="none", ms=4.2, color=INK)]
ax_c.legend(h, ["null run (arm colour)", "null run, no marker pre-sel.", "TRUE labels"], frameon=False,
            loc="upper left", ncol=2, columnspacing=0.8, handletextpad=0.4, borderaxespad=0.1)
ax_c.set_title(f"c  q<{FDR} hits per null run vs TRUE run")

# --- d: share of TRUE nominal hits surviving permutation calibration ---------
w = 0.36
for xi, arm in enumerate(arm_order):
    d, st = data[arm], stats[arm]
    n0 = max(st["true_q_hits"], 1)
    v1 = 100.0 * st["true_nominal_hits_surviving_perm_q"] / n0
    v2 = 100.0 * st["true_perm_q_hits_with_effect_floor"] / n0
    kw = dict(color=d["color"], lw=0) if d["markers"] else dict(facecolor="none", edgecolor=d["color"], lw=1.0)
    ax_d.bar(xi - w / 2, v1, width=w * 0.92, zorder=3, **kw)
    ax_d.bar(xi + w / 2, v2, width=w * 0.92, zorder=3, hatch="////", **({**kw, "edgecolor": SURFACE} if d["markers"] else kw))
    ax_d.text(xi - w / 2, v1 + 1.5, f"{st['true_nominal_hits_surviving_perm_q']:,}", ha="center", va="bottom", fontsize=5.7, color=INK, rotation=90)
    ax_d.text(xi + w / 2, v2 + 1.5, f"{st['true_perm_q_hits_with_effect_floor']:,}", ha="center", va="bottom", fontsize=5.7, color=INK, rotation=90)
ax_d.set_xticks(range(len(arm_order)))
ax_d.set_xticklabels([f"{data[a]['short']}\n({stats[a]['true_q_hits']:,})" for a in arm_order])
style(ax_d); ax_d.set_ylim(0, 148); ax_d.set_yticks([0, 25, 50, 75, 100])
ax_d.set_ylabel(f"% of TRUE nominal q<{FDR} hits retained")
hd = [Patch(facecolor=MUTED, lw=0), Patch(facecolor=MUTED, hatch="////", edgecolor=SURFACE, lw=0),
      Patch(facecolor="none", edgecolor=MUTED, lw=1.0)]
ax_d.legend(hd, [f"permutation-calibrated q<{FDR}", f"+ effect floor (|dprop|>={EFFECT_FLOOR}; nb |log2FC|>=1)", "no marker pre-selection"],
            frameon=False, loc="upper left", ncol=1, handletextpad=0.4, borderaxespad=0.1)
ax_d.set_title("d  TRUE hits surviving permutation calibration (counts on bars)")
ax_d.set_xlabel("test arm (nominal q<0.05 hits in TRUE run)")

# --- title + caveats --------------------------------------------------------
calib = [a for a in stats if stats[a]["calibrated_by_prereg_rule"]]
notcal = [a for a in stats if a not in calib]
if not calib:
    verdict = "no arm is FDR-calibrated at face value"
else:
    verdict = ("valid FDR control: " + " / ".join(short1(data[a]) + (" (conservative)" if stats[a]["conservative"] else " (calibrated)") for a in calib)
               + "; anti-conservative: " + " / ".join(short1(data[a]) for a in notcal))
verdict_sentence = (
    ("No test is FDR-calibrated under the label-permutation null; permutation-calibrated q-values are required for any shipped switch call."
     if not calib else
     "Only " + " and ".join(data[a]["label"] for a in calib) + " gives valid FDR control under the label-permutation null"
     + (" (conservative: " + ", ".join(f"{stats[a]['frac_null_p_lt_05']:.1%} null p<0.05, {stats[a]['frac_null_q_lt_fdr']:.2%} null q<0.05" for a in calib) + ")")
     + "; every arm that keeps the default top-200 marker pre-selection is anti-conservative ("
     + ", ".join(f"{short1(data[a])} {stats[a]['frac_null_p_lt_05']:.0%} null p<0.05" for a in notcal if data[a]["markers"]) + ")"
     + (", as is " + ", ".join(short1(data[a]) for a in notcal if not data[a]["markers"]) + " without pre-selection" if any(not data[a]["markers"] for a in notcal) else "")
     + "; permutation-calibrated q is NOT required for the calibrated arm but IS required whenever marker pre-selection or an anti-conservative test is used."))
fig.suptitle(textwrap.fill(f"`switch diff` calibration on a correctly keyed count matrix -- {verdict} "
             f"(testis mouse1 clip-seeded v2, {info['n_cells']:,} STARsolo cells, SPC/RS/ES, 3 stage pairs per run, 20 label permutations)", 118),
             x=0.01, ha="left", fontsize=8.0, fontweight="bold", y=1.012)
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
    f"diagnostic arms (dashed/hollow) switch it off (--marker-top-n 0; {diag_note}); BH families are per stage pair. "
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
    + "; fisher p-values are discrete (mass at p=1). Palette six-checks validated 2026-08-19."
)
fig.text(0.01, -0.015, textwrap.fill(caveat, 168), ha="left", va="top", fontsize=5.6, color=MUTED, linespacing=1.4)
fig.tight_layout(rect=(0, 0.0, 1, 0.985), h_pad=1.8, w_pad=1.4)
fig.savefig(f"{TSVDIR}/{NAME}.png", dpi=300, bbox_inches="tight")
for ext in ("png", "pdf"):
    p = f"{FIGDIR}/{NAME}.{ext}"
    fig.savefig(p, bbox_inches="tight", **({"dpi": 300} if ext == "png" else {}))
    print("wrote", p)

# ----------------------------------------------------------------------------
# 6. machine-readable report
# ----------------------------------------------------------------------------
report = dict(
    name=NAME, input=info, fdr=FDR, n_perms_planned=N_PERMS_PLANNED, main_arms=main_arms, diagnostic_arms=diag_arms,
    arms={a: {k: (v.item() if isinstance(v, (np.floating, np.integer, np.bool_)) else v) for k, v in stats[a].items()} for a in stats},
    verdict=verdict, verdict_sentence=verdict_sentence, prereg_rule=PREREG_RULE,
    outputs=dict(figure_png=f"{FIGDIR}/{NAME}.png", figure_pdf=f"{FIGDIR}/{NAME}.pdf",
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
