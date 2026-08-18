#!/usr/bin/env python3
"""
fdr_calibration.py -- is `ema switch diff` (fisher strategy) FDR honest under
a label-permutation null?

Design
------
TRUE run    results/fdr_calibration/true/          : real obs['stage'] labels.
NULL runs   results/fdr_calibration/null/perm_XX/  : obs['stage'] permuted
            across cells (label counts preserved, numpy seed = XX), then the
            FULL `switch diff` pipeline re-run end-to-end -- including the
            top-200-per-cluster marker pre-selection, so any double-dipping
            from selecting markers on the same labels being tested is part of
            what the null measures.

The script reads WHATEVER permutations are complete (perm dirs containing
DONE.ok) and states the count on the figure, so it can be re-run once all 20
exist. Dataset: Laughney TREG per-celltype h5ad (3,790 cells x 115,450 PAS),
6 stage pairs per run; BH is applied per pair (per-pair families).

Outputs
-------
manuscript/figures/fdr_calibration.{png,pdf}     (300 dpi / fonttype 42)
results/figures/manuscript/fdr_calibration.tsv                (per-run summary)
results/figures/manuscript/fdr_calibration_null_pvalues.tsv   (pooled null p/q)
"""
import glob
import os
import re

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

FDR = 0.05
N_PERMS_PLANNED = 20

# validated palette (Okabe-Ito subset; CVD deltaE(OKLab*100) >= 8.8 all pairs,
# normal-vision >= 18.7, contrast on white >= 3.4)
BLUE, VERM, GREEN = "#0072B2", "#D55E00", "#009E73"
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"

# ----------------------------------------------------------------------------
# 1. harvest
# ----------------------------------------------------------------------------
def read_run(diff_dir: str) -> pd.DataFrame:
    """All per-pair fisher TSVs of one run -> long df (pair, pas_id, p, q)."""
    frames = []
    for f in sorted(glob.glob(os.path.join(diff_dir, "fisher_*_vs_*.tsv"))):
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        pair = re.sub(r"^fisher_|\.tsv$", "", os.path.basename(f))
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

true_df = read_run(os.path.join(RUNDIR, "true", "differential"))
assert len(true_df), f"TRUE run differential TSVs not found under {RUNDIR}/true"

perm_frames, perm_ids = [], []
for pdir in sorted(glob.glob(os.path.join(RUNDIR, "null", "perm_*"))):
    if not os.path.exists(os.path.join(pdir, "DONE.ok")):
        continue  # incomplete / failed permutation: skip, and say so on the figure
    k = int(os.path.basename(pdir).split("_")[1])
    df = read_run(os.path.join(pdir, "differential"))
    df.insert(0, "perm", k)
    perm_frames.append(df)
    perm_ids.append(k)
n_perms = len(perm_ids)
assert n_perms > 0, "no completed permutations found"
null_df = pd.concat(perm_frames, ignore_index=True)

# ----------------------------------------------------------------------------
# 2. per-run summary table
# ----------------------------------------------------------------------------
rows = []
def summarize(tag, seed, df):
    by_pair = df.groupby("pair")["qvalue"].apply(lambda q: (q < FDR).sum())
    rows.append(dict(
        run=tag, seed=seed, n_tests=len(df),
        n_pairs_tested=df["pair"].nunique(),
        n_q_lt_fdr=int((df["qvalue"] < FDR).sum()),
        n_pair_families_with_hit=int((by_pair > 0).sum()),
        min_pvalue=float(df["pvalue"].min()) if len(df) else np.nan,
        median_pvalue=float(df["pvalue"].median()) if len(df) else np.nan,
    ))
summarize("true", "", true_df)
for k in perm_ids:
    summarize(f"perm_{k:02d}", k, null_df[null_df["perm"] == k])
run_tab = pd.DataFrame(rows)
run_tab.to_csv(f"{TSVDIR}/fdr_calibration.tsv", sep="\t", index=False)
null_df.to_csv(f"{TSVDIR}/fdr_calibration_null_pvalues.tsv", sep="\t", index=False)
print(f"wrote {TSVDIR}/fdr_calibration.tsv ({len(run_tab)} runs)")
print(f"wrote {TSVDIR}/fdr_calibration_null_pvalues.tsv ({len(null_df)} null tests)")

null_runs = run_tab[run_tab["run"] != "true"]
frac_any = float((null_runs["n_q_lt_fdr"] > 0).mean())
n_fam = int(null_runs["n_pairs_tested"].sum())
n_fam_hit = int(null_runs["n_pair_families_with_hit"].sum())
p_null = null_df["pvalue"].values
n_null = len(p_null)
frac_p05 = float((p_null < 0.05).mean())

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

# --- A: null p-value histogram -------------------------------------------
nbins = 20
cnt, edges = np.histogram(p_null, bins=nbins, range=(0, 1))
dens = cnt / n_null * nbins            # density: uniform => 1.0 everywhere
ax_a.bar(edges[:-1], dens, width=1 / nbins * 0.92, align="edge",
         color=BLUE, edgecolor=SURFACE, linewidth=0.5, zorder=3)
ax_a.axhline(1.0, color=VERM, lw=1.4, ls="--", zorder=4)
ax_a.text(0.985, 1.06, "uniform (honest null)", color=VERM, fontsize=6.8,
          ha="right", va="bottom")
ax_a.text(edges[1] + 0.03, dens[0] * 0.99,
          f"{cnt[0]:,} of {n_null:,} in first bin\n(uniform expectation: {n_null // nbins})",
          ha="left", va="top", fontsize=6.8, color=INK)
style(ax_a)
ax_a.set_xlim(0, 1)
ax_a.set_xlabel("null p-value")
ax_a.set_ylabel("density (uniform = 1)")
ax_a.set_title(f"a  Null p-values are far from uniform\n    ({n_perms} permutation runs pooled, {n_null:,} tests)")

# --- B: QQ vs uniform ----------------------------------------------------
p_sorted = np.sort(p_null)
FLOOR = 1e-30                          # display floor; smallest p can be ~1e-70
obs = -np.log10(np.maximum(p_sorted, FLOOR))
exp = -np.log10((np.arange(1, n_null + 1) - 0.5) / n_null)
ii = np.arange(1, n_null + 1)
band_lo = -np.log10(beta_dist.ppf(0.975, ii, n_null - ii + 1))
band_hi = -np.log10(beta_dist.ppf(0.025, ii, n_null - ii + 1))
order = np.argsort(exp)
ax_b.fill_between(exp[order], band_lo[order], band_hi[order],
                  color=GRID, alpha=0.65, lw=0, zorder=2,
                  label="95% pointwise band (uniform)")
lim = max(exp.max(), 1.0) * 1.05
ax_b.plot([0, lim], [0, lim], color=VERM, lw=1.2, ls="--", zorder=3,
          label="y = x (honest null)")
ax_b.plot(exp, obs, ".", ms=3.5, color=BLUE, zorder=4, label="observed null p")
n_clip = int((p_sorted < FLOOR).sum())
if n_clip:
    ax_b.text(0.02, 0.98, f"{n_clip} p-values < 1e-30 clipped to axis top",
              transform=ax_b.transAxes, ha="left", va="top",
              fontsize=6.6, color=MUTED)
style(ax_b)
ax_b.set_xlim(0, lim)
ax_b.set_xlabel(r"expected $-\log_{10}p$ (uniform)")
ax_b.set_ylabel(r"observed $-\log_{10}p$")
ax_b.legend(loc="lower right", frameon=False, handletextpad=0.6)
ax_b.set_title("b  QQ vs uniform: extreme inflation")

# --- C: q<0.05 hits per run ----------------------------------------------
xs = np.arange(len(run_tab))
colors = [VERM if r == "true" else BLUE for r in run_tab["run"]]
ax_c.bar(xs, run_tab["n_q_lt_fdr"], color=colors, width=0.8,
         edgecolor=SURFACE, linewidth=0.5, zorder=3)
for x, v in zip(xs, run_tab["n_q_lt_fdr"]):
    ax_c.text(x, v + 0.6, str(v), ha="center", va="bottom",
              fontsize=6.2, color=INK)
ax_c.set_xticks(xs)
ax_c.set_xticklabels(
    ["TRUE"] + [f"{k:02d}" for k in perm_ids],
    rotation=90 if len(xs) > 12 else 0, fontsize=6.2)
style(ax_c)
ax_c.set_xlabel("run (TRUE labels vs permutation seed)")
ax_c.set_ylabel(f"PAS with q < {FDR}")
ax_c.set_title(f"c  Every null run reports q<{FDR} hits\n    ({frac_any:.0%} of runs; "
               f"{n_fam_hit}/{n_fam} BH families with ≥ 1 hit, expect ≈ 5%)")
ax_c.set_ylim(0, run_tab["n_q_lt_fdr"].max() * 1.34)
handles = [plt.Rectangle((0, 0), 1, 1, color=VERM),
           plt.Rectangle((0, 0), 1, 1, color=BLUE)]
ax_c.legend(handles, ["real stage labels", "permuted labels (null)"],
            frameon=False, loc="upper right", borderaxespad=0.2)

fig.suptitle(
    "ema switch diff (fisher) is NOT FDR-calibrated under a stage-label permutation null "
    "-- Laughney TREG, 3,790 cells, 6 stage pairs per run",
    x=0.005, ha="left", fontsize=9.2, fontweight="bold", y=1.045)
caveat = (
    f"Null = obs['stage'] shuffled across cells (label counts kept), full pipeline re-run per permutation, seeds 1-{N_PERMS_PLANNED}. "
    f"{n_perms}/{N_PERMS_PLANNED} permutations available when drawn"
    + (" (PRELIMINARY -- re-run when all complete). " if n_perms < N_PERMS_PLANNED else ". ")
    + f"{frac_p05:.0%} of null p-values < 0.05 (expect 5%). Drivers of inflation: Fisher on read counts pseudoreplicates "
    "(reads within a cell are correlated; documented anti-conservative in ema code), and the top-200-per-cluster marker "
    "pre-selection uses the same labels that are then tested (selection bias is deliberately included in this null). "
    "Fisher-exact p-values are discrete, so mild deviation from uniform is expected at low counts -- not 40+ orders of magnitude."
)
fig.text(0.005, -0.045, caveat, ha="left", va="top", fontsize=6.4,
         color=MUTED, wrap=True)
fig.tight_layout(rect=(0, 0.02, 1, 0.97))

for ext in ("png", "pdf"):
    p = os.path.join(FIGDIR, f"fdr_calibration.{ext}")
    fig.savefig(p, **({"dpi": 300} if ext == "png" else {}), bbox_inches="tight")
    print("wrote", p)

# terminal summary
print(f"\nTRUE run: {run_tab.iloc[0]['n_tests']} tests, {run_tab.iloc[0]['n_q_lt_fdr']} q<{FDR}")
print(f"NULL: {n_perms} runs, {n_null} tests, {frac_p05:.1%} p<0.05, "
      f"{frac_any:.0%} of runs with any q<{FDR} hit, "
      f"{n_fam_hit}/{n_fam} BH families with a hit")
