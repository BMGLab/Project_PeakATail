#!/usr/bin/env python3
"""
pbmc_novelty_real_pas_usage.py -- STAGE 6, the decisive experiment.

Within-gene USAGE space recovers cell types (stage 3). But 89.4% of the sites
carrying that signal are not within 100 bp of a PolyASite 2.0 representative
site, and PeakATail calls up to 71 "PAS" per gene here -- so "within-gene
usage" over that feature set is a within-gene COVERAGE-SHAPE signature, not
poly(A)-site choice.

This restricts the usage space to sites that are actually poly(A) sites:
genes with >=2 atlas-supported PAS, usage computed as each atlas-supported
site's share of that gene's atlas-supported total. If cell-type structure
survives, alternative polyadenylation genuinely carries cell identity here.
If it collapses while the unrestricted (mostly non-PAS) version does not, the
signal is coverage shape and the APA interpretation is unsupported.

Matched controls, all through the identical PeakATail pipeline:
  real_pas_usage        usage over atlas-supported sites  (n features = N_R)
  decoy_usage_matched   usage over NON-atlas sites, same n features, matched
                        on usage variance rank -- the like-for-like decoy
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
PAS_CLUST = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                    "07_clustering/default/clusters.h5ad")
PAS_RAW = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                  "06_preprocessing/default/preprocessed.h5ad")
OUT = ROOT / "results/pbmc_novelty"
SEED, ATLAS_WIN = 42, 100
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
pos = pd.Series(np.arange(len(obs_names)), index=obs_names)
ROWS = pos.loc[obs_names[M["pas_row"].values]].values
CT = M["celltype"].values
cell_tot = np.asarray(counts.sum(axis=1)).ravel()


def usage_space(mask, label, n_want=None, rank_by=None):
    """Build a within-gene usage matrix over the PAS selected by `mask`.
    Usage is each site's share of its gene's total OVER THE SAME SITE SET, so
    the atlas-restricted space is a self-consistent APA measurement."""
    cols = np.where(mask)[0]
    gc = gene_codes[cols]
    n_per = np.bincount(gc, minlength=n_genes)
    ok = n_per[gc] >= 2                     # gene needs >=2 sites in this set
    cols = cols[ok]
    gc = gc[ok]
    log(f"[{label}] {cols.size} sites in "
        f"{int((np.bincount(gc, minlength=n_genes) >= 2).sum())} genes with >=2")
    Gs = sp.csr_matrix((np.ones(cols.size), (gc, np.arange(cols.size))),
                       shape=(n_genes, cols.size))
    sub = counts[:, cols]
    genec = np.asarray((sub @ Gs.T).todense(), np.float32)
    det = genec >= 3
    det_n = det.sum(axis=0)
    subc = sub.tocoo()
    tv = genec[subc.row, gc[subc.col]]
    keep = tv >= 3
    f_nz = np.zeros(subc.nnz)
    f_nz[keep] = subc.data[keep] / tv[keep]
    n_obs = det_n[gc].astype(float)
    s1 = np.bincount(subc.col, weights=f_nz, minlength=cols.size)
    s2 = np.bincount(subc.col, weights=f_nz ** 2, minlength=cols.size)
    mean_f = np.divide(s1, n_obs, out=np.zeros_like(s1), where=n_obs > 0)
    var_f = np.divide(s2, n_obs, out=np.zeros_like(s2),
                      where=n_obs > 0) - mean_f ** 2
    var_f[n_obs <= 50] = 0.0
    order = np.argsort(-var_f)
    if n_want is not None:
        order = order[:n_want]
    top = order
    log(f"[{label}] kept {top.size} features, usage variance "
        f"{var_f[top].min():.2e}-{var_f[top].max():.2e}")
    colmap = np.full(cols.size, -1, np.int64)
    colmap[top] = np.arange(top.size)
    Om = det[:, gc[top]]
    U = np.tile(mean_f[top].astype(np.float32), (counts.shape[0], 1))
    U[Om] = 0.0
    cc = colmap[subc.col]
    m = keep & (cc >= 0)
    U[subc.row[m], cc[m]] = f_nz[m].astype(np.float32)
    return U, cols[top], var_f[top]


def run(name, X, prenormalised=True):
    A = ad_mod.AnnData(np.asarray(X, np.float64))
    A.obs_names = obs_names
    strat = get_strategy("leiden_tfidf", **PARAMS)
    A.obs["total_counts"] = cell_tot
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
    base = {"space": name, "resolution": 1.0,
            "n_clusters": int(len(set(l1))),
            "AMI_vs_celltype": float(adjusted_mutual_info_score(l1, CT)),
            "ARI_vs_celltype": float(adjusted_rand_score(l1, CT))}
    log(f"[{name}] res=1.0 k={base['n_clusters']} "
        f"AMI={base['AMI_vs_celltype']:.4f} ARI={base['ARI_vs_celltype']:.4f}")
    return base, rows


# --- real PAS usage ---------------------------------------------------------
Ur, cols_r, var_r = usage_space(atlas_ok, "real_pas_usage")
N_R = cols_r.size
base_r, sw_r = run("real_pas_usage", Ur)
zr = (Ur - Ur.mean(0, keepdims=True))
sdr = zr.std(0, keepdims=True)
sdr[sdr == 0] = 1
zr /= sdr
base_rs, sw_rs = run("real_pas_usage_scaled", zr)
del zr

# --- matched decoy: NON-atlas sites, same number of features ----------------
Ud, cols_d, var_d = usage_space(~atlas_ok, "decoy_usage_matched", n_want=N_R)
base_d, sw_d = run("decoy_usage_matched", Ud)
zd = (Ud - Ud.mean(0, keepdims=True))
sdd = zd.std(0, keepdims=True)
sdd[sdd == 0] = 1
zd /= sdd
base_ds, sw_ds = run("decoy_usage_matched_scaled", zd)
del zd

SW = pd.DataFrame(sw_r + sw_rs + sw_d + sw_ds)
SW.to_csv(OUT / "real_pas_usage_sweep.tsv", sep="\t", index=False)
B = pd.DataFrame([base_r, base_rs, base_d, base_ds])
n_ct = len(set(CT))
km = []
for name, g in SW.groupby("space"):
    g = g.copy()
    g["dk"] = (g["n_clusters"] - n_ct).abs()
    b = g.sort_values(["dk", "resolution"]).iloc[0]
    km.append({"space": name, "res_kmatched": b["resolution"],
               "k_matched": int(b["n_clusters"]),
               "AMI_kmatched": b["AMI_vs_celltype"],
               "ARI_kmatched": b["ARI_vs_celltype"],
               "AMI_best": g["AMI_vs_celltype"].max(),
               "ARI_best": g["ARI_vs_celltype"].max()})
KM = pd.DataFrame(km)
OUTT = B.merge(KM, on="space")
OUTT["n_features"] = [N_R, N_R, cols_d.size, cols_d.size]
log("\n" + OUTT.to_string(index=False))
OUTT.to_csv(OUT / "real_pas_usage.tsv", sep="\t", index=False)
with open(OUT / "real_pas_usage.json", "w") as fh:
    json.dump({"n_real_pas_features": int(N_R),
               "n_decoy_features": int(cols_d.size),
               "median_var_real": float(np.median(var_r)),
               "median_var_decoy": float(np.median(var_d)),
               "table": OUTT.to_dict("records")}, fh, indent=2)
log("[done]")
