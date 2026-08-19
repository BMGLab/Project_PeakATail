#!/usr/bin/env python3
"""
pbmc_novelty_power_matched.py -- STAGE 8: kill the power explanation.

Real poly(A) sites are sparser than PeakATail's non-PAS peaks: 3 sites per
gene vs 10, and a gene's atlas-supported total is >=3 counts in only 3.4% of
cells vs 13.2% for the non-PAS set (results/pbmc_novelty/usage_power_check.json).
So the collapse of the real-PAS usage space could be a MEASUREMENT-POWER
result rather than a biological one.

Control: binomially thin the NON-PAS counts so that their per-gene sequencing
depth matches the real-PAS set, then rebuild the decoy usage space exactly as
before (same variance-matched selection procedure, same feature count, same
pipeline). If the depth-matched decoy still recovers cell types, power is not
the explanation.

Reported alongside is the reverse framing that matters for the manuscript: the
reason real poly(A) sites are under-powered here is that PeakATail assigns most
of each gene's reads to peaks that are not poly(A) sites.
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
SEED, ATLAS_WIN, N_TOP = 42, 100, 10000
PARAMS = dict(resolution=1.0, n_components=50, n_dims=40, n_neighbors=30,
              scale_factor=1e4, depth_corr_threshold=0.75, random_seed=SEED)
RES_GRID = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]


def log(*a):
    print(*a, flush=True)


pas = sc.read_h5ad(ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                           "07_clustering/default/clusters.h5ad"))
with h5py.File(ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                       "06_preprocessing/default/preprocessed.h5ad"), "r") as f:
    counts = sp.csr_matrix((f["X/data"][:].astype(np.float64),
                            f["X/indices"][:], f["X/indptr"][:]),
                           shape=tuple(f["X"].attrs["shape"]))
gene_codes = pas.var["gene_id"].cat.codes.values.astype(np.int32)
n_genes = len(pas.var["gene_id"].cat.categories)
atl = pd.read_csv(OUT / "pas_atlas_distance.tsv", sep="\t")
d = pd.Series(atl["atlas_dist"].values,
              index=atl["pas_id"].values).reindex(pas.var["pas_id"].values).values
atlas_ok = np.nan_to_num(np.where(d < 0, np.nan, d), nan=1e9) <= ATLAS_WIN
M = pd.read_csv(OUT / "matched_cells.tsv", sep="\t")
obs_names = pas.obs_names.values
ROWS = pd.Series(np.arange(len(obs_names)),
                 index=obs_names).loc[obs_names[M["pas_row"].values]].values
CT = M["celltype"].values
cell_tot = np.asarray(counts.sum(axis=1)).ravel()
rng = np.random.default_rng(SEED)


def stats(mask, mat):
    cols = np.where(mask)[0]
    g = gene_codes[cols]
    n = np.bincount(g, minlength=n_genes)
    k = n[g] >= 2
    cols, g = cols[k], g[k]
    Gs = sp.csr_matrix((np.ones(cols.size), (g, np.arange(cols.size))),
                       shape=(n_genes, cols.size))
    sub = mat[:, cols]
    genec = np.asarray((sub @ Gs.T).todense(), np.float32)
    det = genec >= 3
    det_n = det.sum(axis=0)
    co = sub.tocoo()
    tv = genec[co.row, g[co.col]]
    keep = tv >= 3
    f = np.zeros(co.nnz)
    f[keep] = co.data[keep] / tv[keep]
    n_obs = det_n[g].astype(float)
    s1 = np.bincount(co.col, weights=f, minlength=cols.size)
    s2 = np.bincount(co.col, weights=f ** 2, minlength=cols.size)
    mean_f = np.divide(s1, n_obs, out=np.zeros_like(s1), where=n_obs > 0)
    var_f = np.divide(s2, n_obs, out=np.zeros_like(s2),
                      where=n_obs > 0) - mean_f ** 2
    var_f[n_obs <= 50] = 0.0
    genes = np.unique(g)
    return dict(cols=cols, g=g, det=det, co=co, keep=keep, f=f, mean_f=mean_f,
                var_f=var_f,
                det_rate=float(det[:, genes].mean()),
                med_gene_tot=float(np.median(genec[:, genes].sum(0))))


def build(S, take):
    cm = np.full(S["cols"].size, -1, np.int64)
    cm[take] = np.arange(take.size)
    Om = S["det"][:, S["g"][take]]
    U = np.tile(S["mean_f"][take].astype(np.float32), (counts.shape[0], 1))
    U[Om] = 0.0
    cc = cm[S["co"].col]
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
    A = get_strategy("leiden_tfidf", **PARAMS).reduce_dims(A)
    A = get_strategy("leiden_tfidf", **PARAMS).cluster(A)
    rows = []
    for r in RES_GRID:
        sc.tl.leiden(A, resolution=r, key_added=f"_r{r}", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
        l = A.obs[f"_r{r}"].astype(str).values[ROWS]
        rows.append({"space": name, "resolution": r,
                     "n_clusters": int(len(set(l))),
                     "AMI_vs_celltype": float(adjusted_mutual_info_score(l, CT)),
                     "ARI_vs_celltype": float(adjusted_rand_score(l, CT))})
    r1 = [x for x in rows if x["resolution"] == 1.0][0]
    log(f"[{name}] res1.0 k={r1['n_clusters']} AMI={r1['AMI_vs_celltype']:.4f} "
        f"ARI={r1['ARI_vs_celltype']:.4f}")
    return rows


SR = stats(atlas_ok, counts)
SD_full = stats(~atlas_ok, counts)
p = SR["med_gene_tot"] / SD_full["med_gene_tot"]
log(f"[power] real det_rate {SR['det_rate']:.4f} med gene total "
    f"{SR['med_gene_tot']:.0f} | non-PAS det_rate {SD_full['det_rate']:.4f} "
    f"med gene total {SD_full['med_gene_tot']:.0f} -> thinning p = {p:.4f}")

# binomially thin the non-PAS counts
nz = counts[:, ~atlas_ok].tocoo()
thin = rng.binomial(nz.data.astype(np.int64), p)
Cd = sp.csr_matrix(
    (thin.astype(np.float64), (nz.row, nz.col)),
    shape=(counts.shape[0], int((~atlas_ok).sum())))
Cd.eliminate_zeros()
full = sp.csr_matrix((counts.shape[0], counts.shape[1]))
idx_non = np.where(~atlas_ok)[0]
full = sp.csr_matrix(
    (Cd.data, idx_non[Cd.indices], Cd.indptr), shape=counts.shape)
SDt = stats(~atlas_ok, full)
log(f"[power] thinned non-PAS det_rate {SDt['det_rate']:.4f} "
    f"med gene total {SDt['med_gene_tot']:.0f}")

take_r = np.argsort(-SR["var_f"])[:N_TOP]
vr = SR["var_f"][take_r]
ordd = np.argsort(SDt["var_f"])
vd_sorted = SDt["var_f"][ordd]
used = np.zeros(ordd.size, bool)
take_d = np.empty(N_TOP, np.int64)
for i, v in enumerate(vr):
    j = int(np.searchsorted(vd_sorted, v))
    lo, hi, pick = j - 1, j, -1
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
                raise RuntimeError("pool exhausted")
    used[pick] = True
    take_d[i] = ordd[pick]
log(f"[match] variance median real {np.median(vr):.4f} decoy "
    f"{np.median(SDt['var_f'][take_d]):.4f} max|diff| "
    f"{np.abs(vr - SDt['var_f'][take_d]).max():.2e}")

rows = []
Ud = build(SDt, take_d)
rows += run("decoy_usage_depth_and_var_matched", Ud)
rows += run("decoy_usage_depth_and_var_matched_scaled", Ud, scaled=True)
SW = pd.DataFrame(rows)
SW.to_csv(OUT / "power_matched_sweep.tsv", sep="\t", index=False)
n_ct = len(set(CT))
km = []
for name, g in SW.groupby("space"):
    g = g.copy()
    g["dk"] = (g["n_clusters"] - n_ct).abs()
    b = g.sort_values(["dk", "resolution"]).iloc[0]
    r1 = g[g.resolution == 1.0].iloc[0]
    km.append({"space": name, "k_res1": int(r1["n_clusters"]),
               "AMI_res1": r1["AMI_vs_celltype"],
               "ARI_res1": r1["ARI_vs_celltype"],
               "k_matched": int(b["n_clusters"]),
               "AMI_kmatched": b["AMI_vs_celltype"],
               "ARI_kmatched": b["ARI_vs_celltype"],
               "AMI_best": g["AMI_vs_celltype"].max()})
KM = pd.DataFrame(km)
log("\n" + KM.to_string(index=False))
KM.to_csv(OUT / "power_matched.tsv", sep="\t", index=False)
with open(OUT / "power_matched.json", "w") as fh:
    json.dump({"thinning_p": float(p),
               "real_det_rate": SR["det_rate"],
               "nonpas_det_rate": SD_full["det_rate"],
               "nonpas_thinned_det_rate": SDt["det_rate"],
               "real_med_gene_total": SR["med_gene_tot"],
               "nonpas_med_gene_total": SD_full["med_gene_tot"],
               "nonpas_thinned_med_gene_total": SDt["med_gene_tot"],
               "table": KM.to_dict("records")}, fh, indent=2)
log("[done]")
