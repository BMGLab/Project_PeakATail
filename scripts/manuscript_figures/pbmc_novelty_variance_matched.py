#!/usr/bin/env python3
"""
pbmc_novelty_variance_matched.py -- STAGE 7, the fair version of stage 6.

Stage 6 compared ALL 28,381 atlas-supported sites against the TOP 28,381
non-atlas sites by usage variance, so the decoy enjoyed a selection advantage.
This removes that advantage two ways:

  real_top          the top-N atlas-supported sites by usage variance
  decoy_var_matched N non-atlas sites chosen to match the real set's usage
                    variance distribution one-for-one (nearest variance,
                    without replacement)

If the decoy still recovers cell types and the real poly(A) sites still do not,
the difference is not feature count and not feature variance -- it is what the
features are.
"""
import os
os.environ.setdefault("LC_ALL", "C")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

import json
import sys
from pathlib import Path

import anndata as ad_mod
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

sys.path.insert(0, "/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail")
from ema.clustering.strategies import get_strategy       # noqa: E402

sc.settings.n_jobs = 8
ROOT = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUT = ROOT / "results/pbmc_novelty"
PAS_CLUST = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                    "07_clustering/default/clusters.h5ad")
PAS_RAW = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                  "06_preprocessing/default/preprocessed.h5ad")
SEED, ATLAS_WIN, N_TOP = 42, 100, 10000
PARAMS = dict(resolution=1.0, n_components=50, n_dims=40, n_neighbors=30,
              scale_factor=1e4, depth_corr_threshold=0.75, random_seed=SEED)
RES_GRID = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]


def log(*a):
    print(*a, flush=True)


pas = sc.read_h5ad(PAS_CLUST)
with h5py.File(PAS_RAW, "r") as f:
    counts = sp.csr_matrix((f["X/data"][:].astype(np.float64),
                            f["X/indices"][:], f["X/indptr"][:]),
                           shape=tuple(f["X"].attrs["shape"]))
gene_codes = pas.var["gene_id"].cat.codes.values.astype(np.int32)
n_genes = len(pas.var["gene_id"].cat.categories)
pas_id_arr = pas.var["pas_id"].values
atl = pd.read_csv(OUT / "pas_atlas_distance.tsv", sep="\t")
d = pd.Series(atl["atlas_dist"].values,
              index=atl["pas_id"].values).reindex(pas_id_arr).values
atlas_ok = np.nan_to_num(np.where(d < 0, np.nan, d), nan=1e9) <= ATLAS_WIN
M = pd.read_csv(OUT / "matched_cells.tsv", sep="\t")
obs_names = pas.obs_names.values
ROWS = pd.Series(np.arange(len(obs_names)),
                 index=obs_names).loc[obs_names[M["pas_row"].values]].values
CT = M["celltype"].values
cell_tot = np.asarray(counts.sum(axis=1)).ravel()


def usage_stats(mask):
    """Within-gene usage variance for every site in `mask`, usage taken as the
    site's share of its gene's total OVER THE SAME SITE SET."""
    cols = np.where(mask)[0]
    gc = gene_codes[cols]
    n_per = np.bincount(gc, minlength=n_genes)
    ok = n_per[gc] >= 2
    cols, gc = cols[ok], gc[ok]
    Gs = sp.csr_matrix((np.ones(cols.size), (gc, np.arange(cols.size))),
                       shape=(n_genes, cols.size))
    sub = counts[:, cols]
    genec = np.asarray((sub @ Gs.T).todense(), np.float32)
    det = genec >= 3
    det_n = det.sum(axis=0)
    co = sub.tocoo()
    tv = genec[co.row, gc[co.col]]
    keep = tv >= 3
    f = np.zeros(co.nnz)
    f[keep] = co.data[keep] / tv[keep]
    n_obs = det_n[gc].astype(float)
    s1 = np.bincount(co.col, weights=f, minlength=cols.size)
    s2 = np.bincount(co.col, weights=f ** 2, minlength=cols.size)
    mean_f = np.divide(s1, n_obs, out=np.zeros_like(s1), where=n_obs > 0)
    var_f = np.divide(s2, n_obs, out=np.zeros_like(s2),
                      where=n_obs > 0) - mean_f ** 2
    var_f[n_obs <= 50] = 0.0
    return dict(cols=cols, gc=gc, det=det, co=co, keep=keep, f=f,
                mean_f=mean_f, var_f=var_f)


def build(S, take):
    colmap = np.full(S["cols"].size, -1, np.int64)
    colmap[take] = np.arange(take.size)
    Om = S["det"][:, S["gc"][take]]
    U = np.tile(S["mean_f"][take].astype(np.float32), (counts.shape[0], 1))
    U[Om] = 0.0
    cc = colmap[S["co"].col]
    m = S["keep"] & (cc >= 0)
    U[S["co"].row[m], cc[m]] = S["f"][m].astype(np.float32)
    return U


def run(name, U, scaled=False):
    X = U
    if scaled:
        X = U - U.mean(0, keepdims=True)
        sd = X.std(0, keepdims=True)
        sd[sd == 0] = 1
        X = X / sd
    A = ad_mod.AnnData(np.asarray(X, np.float64))
    A.obs_names = obs_names
    A.obs["total_counts"] = cell_tot
    strat = get_strategy("leiden_tfidf", **PARAMS)
    A = strat.reduce_dims(A)
    A = strat.cluster(A)
    rows = []
    for r in RES_GRID:
        sc.tl.leiden(A, resolution=r, key_added=f"_r{r}", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
        l = A.obs[f"_r{r}"].astype(str).values[ROWS]
        rows.append({"space": name, "resolution": r,
                     "n_clusters": int(len(set(l))),
                     "AMI_vs_celltype": float(adjusted_mutual_info_score(l, CT)),
                     "ARI_vs_celltype": float(adjusted_rand_score(l, CT))})
    l1 = A.obs["leiden"].astype(str).values[ROWS]
    log(f"[{name}] res1.0 k={len(set(l1))} "
        f"AMI={adjusted_mutual_info_score(l1, CT):.4f} "
        f"ARI={adjusted_rand_score(l1, CT):.4f}")
    return rows


SR = usage_stats(atlas_ok)
SD = usage_stats(~atlas_ok)
log(f"[sets] real {SR['cols'].size} sites | decoy {SD['cols'].size} sites")

take_r = np.argsort(-SR["var_f"])[:N_TOP]
vr = SR["var_f"][take_r]
# nearest-variance decoy, without replacement
ordd = np.argsort(SD["var_f"])
vd_sorted = SD["var_f"][ordd]
used = np.zeros(ordd.size, bool)
take_d = np.empty(N_TOP, np.int64)
for i, v in enumerate(vr):
    j = int(np.searchsorted(vd_sorted, v))
    lo, hi = j - 1, j
    pick = -1
    while pick < 0:
        cl = lo >= 0 and not used[lo]
        ch = hi < ordd.size and not used[hi]
        if cl and ch:
            pick = lo if abs(vd_sorted[lo] - v) <= abs(vd_sorted[hi] - v) else hi
        elif cl:
            pick = lo
        elif ch:
            pick = hi
        else:
            lo -= 1
            hi += 1
            if lo < 0 and hi >= ordd.size:
                raise RuntimeError("decoy pool exhausted")
            continue
    used[pick] = True
    take_d[i] = ordd[pick]
vd = SD["var_f"][take_d]
log(f"[match] real usage variance median {np.median(vr):.4f} "
    f"[{vr.min():.4f}-{vr.max():.4f}] | decoy median {np.median(vd):.4f} "
    f"[{vd.min():.4f}-{vd.max():.4f}] | max |diff| {np.abs(vr - vd).max():.2e}")

rows = []
Ur = build(SR, take_r)
rows += run("real_pas_usage_top", Ur)
rows += run("real_pas_usage_top_scaled", Ur, scaled=True)
del Ur
Ud = build(SD, take_d)
rows += run("decoy_usage_var_matched", Ud)
rows += run("decoy_usage_var_matched_scaled", Ud, scaled=True)
del Ud

SW = pd.DataFrame(rows)
SW.to_csv(OUT / "variance_matched_sweep.tsv", sep="\t", index=False)
n_ct = len(set(CT))
km = []
for name, g in SW.groupby("space"):
    g = g.copy()
    g["dk"] = (g["n_clusters"] - n_ct).abs()
    b = g.sort_values(["dk", "resolution"]).iloc[0]
    r1 = g[g.resolution == 1.0].iloc[0]
    km.append({"space": name, "n_features": N_TOP,
               "k_res1": int(r1["n_clusters"]),
               "AMI_res1": r1["AMI_vs_celltype"],
               "ARI_res1": r1["ARI_vs_celltype"],
               "res_kmatched": b["resolution"], "k_matched": int(b["n_clusters"]),
               "AMI_kmatched": b["AMI_vs_celltype"],
               "ARI_kmatched": b["ARI_vs_celltype"],
               "AMI_best": g["AMI_vs_celltype"].max(),
               "ARI_best": g["ARI_vs_celltype"].max()})
KM = pd.DataFrame(km)
log("\n" + KM.to_string(index=False))
KM.to_csv(OUT / "variance_matched.tsv", sep="\t", index=False)
with open(OUT / "variance_matched.json", "w") as fh:
    json.dump({"n_features": N_TOP,
               "n_real_pool": int(SR["cols"].size),
               "n_decoy_pool": int(SD["cols"].size),
               "var_real_median": float(np.median(vr)),
               "var_decoy_median": float(np.median(vd)),
               "var_max_abs_diff": float(np.abs(vr - vd).max()),
               "table": KM.to_dict("records")}, fh, indent=2)
log("[done]")
