#!/usr/bin/env python3
"""
fdr_calibration.py -- is `ema switch diff` FDR honest under a label-permutation
null?  Side-by-side calibration of BOTH strategies: fisher and nb_pairwise.

Design
------
TRUE runs   results/fdr_calibration/true/            fisher, real obs['stage']
            results/fdr_calibration/true_nb*/        nb_pairwise, real labels
NULL runs   results/fdr_calibration/null/perm_XX/    fisher
            results/fdr_calibration/null_nb/perm_XX/ nb_pairwise
            obs['stage'] permuted across cells (label counts preserved, numpy
            seed = XX), then the FULL `switch diff` pipeline re-run end-to-end
            -- including the top-200-per-cluster marker pre-selection, so any
            double-dipping from selecting markers on the labels being tested is
            part of what the null measures.

The script reads WHATEVER permutations are complete (perm dirs containing
DONE.ok) and states the count on the figure, so it stays re-runnable.
Dataset: Laughney TREG per-celltype h5ad (3,790 cells x 115,450 PAS), 6 stage
pairs per run; BH is applied per pair (per-pair families).

Outputs
-------
results/figures/manuscript/fdr_calibration.png              (300 dpi)
manuscript/figures/fdr_calibration.{png,pdf}                (300 dpi / fonttype 42)
results/figures/manuscript/fdr_calibration.tsv              (per-run summary, both methods)
results/figures/manuscript/fdr_calibration_null_pvalues.tsv (pooled null p/q, both methods)
results/figures/manuscript/fdr_calibration_hist.tsv         (panel-a bin densities)
results/figures/manuscript/fdr_calibration_stats.tsv        (headline stats per method)
"""
import glob
import os
import re
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist

# --- manuscript deliverables: PNG + vector PDF into manuscript/figures/ ---
import matplotlib as _mpl
_mpl.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (journal requirement)
_mpl.rcParams['ps.fonttype'] = 42
WD = '/mnt/ssd1/Projects/PeakATail_wd'
FIGDIR = f'{WD}/manuscript/figures'
TSVDIR = f'{WD}/results/figures/manuscript'
RUNDIR = f'{WD}/results/fdr_calibration'
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TSVDIR, exist_ok=True)

def save_manuscript(fig, name, **kw):
    for ext in ('png', 'pdf'):
        p = os.path.join(FIGDIR, f'{name}.{ext}')
        fig.savefig(p, **({'dpi': 300} if ext == 'png' else {}), **kw)
        print('wrote', p)

FDR = 0.05
N_PERMS_PLANNED = 20
NAME = "fdr_calibration"

# validated palette (Okabe-Ito subset; re-validated 2026-08-19 with the dataviz
# six-checks script: worst adjacent pair CVD deltaE(OKLab*100) = 11.0 deutan /
# 8.6 tritan, normal-vision floor 18.7, contrast on white >= 3:1 -- all PASS)
BLUE, VERM, GREEN = "#0072B2", "#D55E00", "#009E73"     # fisher, TRUE/ref, nb
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

METHODS = {  # tag -> (tsv prefix, TRUE dir glob, null perm glob, color, label)
    "fisher": ("fisher", f"{RUNDIR}/true",
               f"{RUNDIR}/null/perm_*", BLUE, "fisher"),
    "nb": ("nb_pairwise", f"{RUNDIR}/true_nb*",
           f"{RUNDIR}/null_nb/perm_*", GREEN, "nb_pairwise"),
}

# ----------------------------------------------------------------------------
# 1. harvest
# ----------------------------------------------------------------------------
def read_run(diff_dir: str, prefix: str) -> pd.DataFrame:
    """All per-pair TSVs of one run -> long df (pair, pas_id, p, q)."""
    frames = []
    for f in sorted(glob.glob(os.path.join(diff_dir, f"{prefix}_*_vs_*.tsv"))):
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        pair = re.sub(rf"^{prefix}_|\.tsv$", "", os.path.basename(f))
        frames.append(pd.DataFrame({
            "pair": pair,
            "pas_id": df["pas_id"].astype(str),
            "gene_id": df.get("gene_id", ""),
            "pvalue": df["pvalue"].astype(float),
            "qvalue": df["qvalue"].astype(float),
        }))
    if not frames:
        return pd.DataFrame(columns=["pair", "pas_id", "gene_id", "pvalue", "qvalue"])
    return pd.concat(frames, ignore_index=True)

def true_dir_for(pattern: str, prefix: str) -> str:
    """Newest dir matching `pattern` that actually holds differential TSVs
    (`ema` appends a timestamp to -o unless PEAKATAIL_NO_TIMESTAMP=1, so the
    nb TRUE run landed in true_nb_<ts>/)."""
    cands = [d for d in sorted(glob.glob(pattern))
             if glob.glob(os.path.join(d, "differential", f"{prefix}_*_vs_*.tsv"))]
    assert cands, f"no TRUE run with {prefix} TSVs matches {pattern}"
    return cands[-1]

data = {}          # tag -> dict(true_df, null_df, perm_ids, ...)
for tag, (prefix, tpat, npat, color, label) in METHODS.items():
    tdir = true_dir_for(tpat, prefix)
    true_df = read_run(os.path.join(tdir, "differential"), prefix)
    assert len(true_df), f"TRUE differential TSVs empty in {tdir}"
    perm_frames, perm_ids = [], []
    for pdir in sorted(glob.glob(npat)):
        if not os.path.exists(os.path.join(pdir, "DONE.ok")):
            continue   # incomplete / failed permutation: skip, say so on figure
        k = int(os.path.basename(pdir).split("_")[1])
        df = read_run(os.path.join(pdir, "differential"), prefix)
        df.insert(0, "perm", k)
        perm_frames.append(df)
        perm_ids.append(k)
    assert perm_ids, f"no completed {tag} permutations under {npat}"
    data[tag] = dict(true_dir=tdir, true_df=true_df, perm_ids=perm_ids,
                     null_df=pd.concat(perm_frames, ignore_index=True),
                     color=color, label=label)
    print(f"{tag}: TRUE={tdir} ({len(true_df)} tests), "
          f"{len(perm_ids)} null perms ({len(data[tag]['null_df'])} tests)")

# ----------------------------------------------------------------------------
# 2. per-run summary + pooled tables
# ----------------------------------------------------------------------------
rows = []
def summarize(method, tag, seed, df):
    by_pair = df.groupby("pair")["qvalue"].apply(lambda q: (q < FDR).sum())
    rows.append(dict(
        method=method, run=tag, seed=seed, n_tests=len(df),
        n_pairs_tested=df["pair"].nunique(),
        n_q_lt_fdr=int((df["qvalue"] < FDR).sum()),
        n_pair_families_with_hit=int((by_pair > 0).sum()),
        frac_p_lt_05=float((df["pvalue"] < 0.05).mean()) if len(df) else np.nan,
        min_pvalue=float(df["pvalue"].min()) if len(df) else np.nan,
        median_pvalue=float(df["pvalue"].median()) if len(df) else np.nan,
    ))
for tag, d in data.items():
    summarize(tag, "true", "", d["true_df"])
    for k in d["perm_ids"]:
        summarize(tag, f"perm_{k:02d}", k, d["null_df"][d["null_df"]["perm"] == k])
run_tab = pd.DataFrame(rows)
run_tab.to_csv(f"{TSVDIR}/{NAME}.tsv", sep="\t", index=False)

pooled = pd.concat(
    [d["null_df"].assign(method=tag)[["method", "perm", "pair", "pas_id", "pvalue", "qvalue"]]
     for tag, d in data.items()], ignore_index=True)
pooled.to_csv(f"{TSVDIR}/{NAME}_null_pvalues.tsv", sep="\t", index=False)
print(f"wrote {TSVDIR}/{NAME}.tsv ({len(run_tab)} runs)")
print(f"wrote {TSVDIR}/{NAME}_null_pvalues.tsv ({len(pooled)} null tests)")

stats = {}
for tag, d in data.items():
    nt = run_tab[(run_tab.method == tag) & (run_tab.run != "true")]
    tt = run_tab[(run_tab.method == tag) & (run_tab.run == "true")].iloc[0]
    p = d["null_df"]["pvalue"].values
    q = d["null_df"]["qvalue"].values
    stats[tag] = dict(
        method=tag, n_perms=len(d["perm_ids"]), n_null_tests=len(p),
        frac_null_p_lt_05=float((p < 0.05).mean()),
        frac_null_q_lt_fdr=float((q < FDR).mean()),
        n_runs_with_q_hit=int((nt["n_q_lt_fdr"] > 0).sum()),
        mean_hits_per_null_run=float(nt["n_q_lt_fdr"].mean()),
        min_null_p=float(p.min()),
        n_fam=int(nt["n_pairs_tested"].sum()),
        n_fam_hit=int(nt["n_pair_families_with_hit"].sum()),
        true_n_tests=int(tt["n_tests"]), true_q_hits=int(tt["n_q_lt_fdr"]),
    )
stat_tab = pd.DataFrame(stats.values())
stat_tab.to_csv(f"{TSVDIR}/{NAME}_stats.tsv", sep="\t", index=False)
print(f"wrote {TSVDIR}/{NAME}_stats.tsv")

# ----------------------------------------------------------------------------
# 3. figure
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 7.8, "axes.titlesize": 8.6,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 6.8,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.linewidth": 0.8, "font.family": "DejaVu Sans",
})
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(11.0, 3.9))

def style(ax):
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# --- A: null p-value histograms, both methods ----------------------------
nbins = 20
hist_rows = []
for tag, d in data.items():
    p = d["null_df"]["pvalue"].values
    cnt, edges = np.histogram(p, bins=nbins, range=(0, 1))
    dens = cnt / len(p) * nbins        # density: uniform => 1.0 everywhere
    ax_a.stairs(dens, edges, fill=True, alpha=0.22, color=d["color"], zorder=3)
    ax_a.stairs(dens, edges, color=d["color"], lw=1.5, zorder=4,
                label=f"{d['label']} null p")
    for lo, hi, c, dv in zip(edges[:-1], edges[1:], cnt, dens):
        hist_rows.append(dict(method=tag, bin_lo=lo, bin_hi=hi,
                              count=int(c), density=dv))
pd.DataFrame(hist_rows).to_csv(f"{TSVDIR}/{NAME}_hist.tsv", sep="\t", index=False)
print(f"wrote {TSVDIR}/{NAME}_hist.tsv")
ax_a.axhline(1.0, color=VERM, lw=1.4, ls="--", zorder=5,
             label="uniform (honest null)")
f_st, n_st = stats["fisher"], stats["nb"]
ax_a.text(0.08, stats["fisher"]["frac_null_p_lt_05"] / 0.05 - 0.32,
          f"fisher: {f_st['frac_null_p_lt_05']:.1%} of null p < 0.05",
          fontsize=6.6, color=INK, ha="left", va="top")
ax_a.text(0.08, stats["nb"]["frac_null_p_lt_05"] / 0.05 + 0.12,
          f"nb: {n_st['frac_null_p_lt_05']:.1%} (expect 5%)",
          fontsize=6.6, color=INK, ha="left", va="bottom")
style(ax_a)
ax_a.set_xlim(0, 1)
ax_a.set_xlabel("null p-value")
ax_a.set_ylabel("density (uniform = 1)")
ax_a.legend(loc="upper right", frameon=False, handlelength=1.6)
ax_a.set_title("a  Null p-values far from uniform for both tests\n"
               f"    (fisher {f_st['n_perms']} perms/{f_st['n_null_tests']:,} tests; "
               f"nb {n_st['n_perms']} perms/{n_st['n_null_tests']:,} tests)")

# --- B: QQ vs uniform, both methods --------------------------------------
FLOOR = 1e-30                      # display floor (smallest p ~1e-36 here)
def qq_arrays(p):
    p_sorted = np.sort(p)
    n = len(p_sorted)
    obs = -np.log10(np.maximum(p_sorted, FLOOR))
    ii = np.arange(1, n + 1)
    exp = -np.log10((ii - 0.5) / n)
    lo = -np.log10(beta_dist.ppf(0.975, ii, n - ii + 1))
    hi = -np.log10(beta_dist.ppf(0.025, ii, n - ii + 1))
    # decimate for rendering only (full p-values live in the *_null_pvalues
    # TSV): keep the 1,000 most significant + ~1,500 evenly rank-spaced rest
    keep = np.unique(np.r_[np.arange(min(1000, n)),
                           np.linspace(0, n - 1, 1500).astype(int)])
    idx = (n - 1) - keep            # ranks from most significant end
    idx = np.sort(idx)
    return exp[idx], obs[idx], lo[idx], hi[idx], int((p_sorted < FLOOR).sum())

lim = 0.0
for tag, d in data.items():
    exp, obs, lo, hi, n_clip = qq_arrays(d["null_df"]["pvalue"].values)
    ax_b.fill_between(exp, lo, hi, color=GRID, alpha=0.6, lw=0, zorder=2)
    ax_b.plot(exp, obs, ".", ms=2.8, color=d["color"], zorder=4,
              label=f"{d['label']} null p", rasterized=False)
    lim = max(lim, exp.max())
    if n_clip:
        ax_b.text(0.02, 0.98, f"{n_clip} nb p < 1e-30 clipped to axis top",
                  transform=ax_b.transAxes, ha="left", va="top",
                  fontsize=6.4, color=MUTED)
lim *= 1.05
ax_b.plot([0, lim], [0, lim], color=VERM, lw=1.2, ls="--", zorder=3,
          label="y = x (honest null)")
h_, l_ = ax_b.get_legend_handles_labels()
ax_b.legend(h_, l_, loc="center right", frameon=False, handletextpad=0.6)
ax_b.text(0.98, 0.02, "gray: 95% pointwise uniform band",
          transform=ax_b.transAxes, ha="right", va="bottom",
          fontsize=6.2, color=MUTED)
style(ax_b)
ax_b.set_xlim(0, lim)
ax_b.set_xlabel(r"expected $-\log_{10}p$ (uniform)")
ax_b.set_ylabel(r"observed $-\log_{10}p$")
ax_b.set_title("b  QQ vs uniform: both inflated\n"
               f"    (min null p: fisher {f_st['min_null_p']:.0e}, nb {n_st['min_null_p']:.0e})")

# --- C: q<0.05 hits per null run, TRUE marked ----------------------------
rng = np.random.default_rng(0)
for xi, tag in enumerate(data):
    d, st = data[tag], stats[tag]
    nt = run_tab[(run_tab.method == tag) & (run_tab.run != "true")]
    xs = xi + rng.uniform(-0.13, 0.13, len(nt))
    ax_c.plot(xs, nt["n_q_lt_fdr"], "o", ms=4.2, mew=0.6, mec=SURFACE,
              color=d["color"], alpha=0.9, zorder=3, ls="none")
    # TRUE run marker + label
    ax_c.plot([xi - 0.24, xi + 0.24], [st["true_q_hits"]] * 2,
              color=VERM, lw=1.6, zorder=4)
    ax_c.plot(xi, st["true_q_hits"], "D", ms=5, color=VERM, mec=SURFACE,
              mew=0.6, zorder=5)
    ax_c.text(xi + 0.27, st["true_q_hits"], f"TRUE: {st['true_q_hits']}",
              ha="left", va="center", fontsize=6.6, color=INK)
    import matplotlib.transforms as mtransforms
    tr = mtransforms.blended_transform_factory(ax_c.transData, ax_c.transAxes)
    ax_c.text(xi, 0.025,
              f"{st['n_runs_with_q_hit']}/{st['n_perms']} runs \u2265 1 hit\n"
              f"mean {st['mean_hits_per_null_run']:.0f}/run",
              transform=tr, ha="center", va="bottom", fontsize=6.4, color=INK)
ax_c.set_yscale("log")
ax_c.set_ylim(2, 1300)
ax_c.set_xlim(-0.55, 1.75)
ax_c.set_xticks([0, 1])
ax_c.set_xticklabels([data[t]["label"] for t in data])
style(ax_c)
ax_c.set_xlabel("test strategy")
ax_c.set_ylabel(f"PAS with q < {FDR} per run (log scale)")
handles = [plt.Line2D([], [], marker="o", ls="none", ms=4.2, color=BLUE),
           plt.Line2D([], [], marker="o", ls="none", ms=4.2, color=GREEN),
           plt.Line2D([], [], marker="D", ls="none", ms=5, color=VERM)]
ax_c.legend(handles, ["fisher null run", "nb null run", "TRUE labels"],
            frameon=False, loc="upper left", borderaxespad=0.2)
ax_c.set_title(f"c  Every null run reports q<{FDR} hits under both tests\n"
               "    (nb null runs average more hits than the TRUE run)")

fig.suptitle(
    "Neither `switch diff` strategy is FDR-calibrated under a stage-label permutation null "
    "-- fisher badly, nb_pairwise ~4x nominal -- Laughney TREG, 3,790 cells, 6 stage pairs per run",
    x=0.005, ha="left", fontsize=9.2, fontweight="bold", y=1.045)

nb_null = data["nb"]["null_df"]
caveat = (
    f"Null = obs['stage'] shuffled across cells (label counts kept), full pipeline re-run per permutation incl. top-200 marker "
    f"pre-selection (double-dip deliberately inside the null), seeds 1-{N_PERMS_PLANNED}; "
    f"fisher {f_st['n_perms']}/{N_PERMS_PLANNED} and nb {n_st['n_perms']}/{N_PERMS_PLANNED} permutations available when drawn"
    + ("" if (f_st['n_perms'] == N_PERMS_PLANNED and n_st['n_perms'] == N_PERMS_PLANNED) else " (PRELIMINARY -- re-run when all complete)")
    + f". BH families are per stage pair: fisher {f_st['n_fam_hit']}/{f_st['n_fam']} and nb {n_st['n_fam_hit']}/{n_st['n_fam']} "
    "null families contain >=1 q<0.05 hit (expect ~5%). "
    "Fisher drivers: the test treats reads as independent (pseudoreplication; reads within a cell are correlated) and Fisher p-values are discrete. "
    "nb_pairwise is a per-cell NB GLM (Wald z on the cluster coefficient, offset log library size) so it fixes pseudoreplication, yet "
    f"{n_st['frac_null_p_lt_05']:.0%} of null p<0.05 and {n_st['frac_null_q_lt_fdr']:.0%} of null q<{FDR} remain: per-PAS plug-in dispersion "
    "(clipped to [1e-4,10]) hits the 1e-4 floor in 21% of null tests -- near-Poisson variance -- and those tests carry ~50% of null p<0.05; no "
    "dispersion shrinkage; marker double-dip included. Hit counts per run are not comparable across TRUE vs null: permuted labels homogenize "
    f"clusters, so more PAS pass the >=10-expressing-cells-per-group filter (null nb runs test ~{len(nb_null)//n_st['n_perms']:,} PAS/run vs "
    f"{n_st['true_n_tests']} in TRUE); the per-test rates above are the calibrated quantities. Verdict: nb_pairwise q-values are also "
    "anti-conservative and cannot be taken at face value -- manuscript significance needs D4 cells-mode or permutation-calibrated thresholds."
)
fig.text(0.005, -0.05, textwrap.fill(caveat, 210), ha="left", va="top",
         fontsize=6.3, color=MUTED, linespacing=1.45)
fig.tight_layout(rect=(0, 0.02, 1, 0.97))

fig.savefig(f"{TSVDIR}/{NAME}.png", dpi=300, bbox_inches="tight")
print("wrote", f"{TSVDIR}/{NAME}.png")
save_manuscript(fig, NAME, bbox_inches="tight")

# terminal summary
for tag, st in stats.items():
    print(f"\n[{tag}] TRUE: {st['true_n_tests']} tests, {st['true_q_hits']} q<{FDR}")
    print(f"[{tag}] NULL: {st['n_perms']} runs, {st['n_null_tests']} tests, "
          f"{st['frac_null_p_lt_05']:.1%} p<0.05, {st['frac_null_q_lt_fdr']:.1%} q<{FDR}, "
          f"{st['n_runs_with_q_hit']}/{st['n_perms']} runs with >=1 q<{FDR} hit "
          f"(mean {st['mean_hits_per_null_run']:.1f}/run), min p {st['min_null_p']:.1e}")
