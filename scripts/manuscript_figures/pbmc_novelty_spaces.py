#!/usr/bin/env python3
"""
pbmc_novelty_spaces.py -- STAGE 3 of the pbmc_novelty analysis.

THE SKEPTICAL CONTROL the concordance claim needs.

"Cells clustered on poly(A)-site profiles alone recover gene-expression cell
types" is only interesting if the recovery comes from poly(A)-site BIOLOGY.
The shipped PAS feature space is per-site UMI COUNTS (TF-IDF -> LSI -> Leiden).
Per-site counts carry gene abundance, so a 3'-end count matrix is close to a
gene-expression matrix with extra rows. This script decomposes the signal by
re-running the IDENTICAL clustering pipeline (PeakATail's LeidenTfidfStrategy, the same
parameters the run used: TF-IDF sf 1e4, 50 SVD comps, depth-corr filter 0.75,
40 dims, 30 neighbours, resolution 1.0, seed 42) on feature spaces that differ
only in what information they retain:

  published        the shipped labels (reference)
  pas_counts       all 275,370 PAS, counts        -- rerun sanity check
  gene_level       PAS collapsed to 14,891 genes  -- pure 3'-end EXPRESSION
  atlas_only       the 32,361 PAS within 100 bp of a PolyASite 2.0 site
  nonatlas_only    the 243,009 PAS that are NOT   -- are the calls even PAS?
  counts_top10k    10,000 selected multi-PAS-gene sites, as counts
  usage_top10k     the SAME 10,000 sites as within-gene USAGE FRACTIONS,
                   i.e. gene abundance divided out -- APA information only

Read the table this writes as: if gene_level ~= pas_counts, the concordance is
an expression result. If usage_top10k collapses relative to counts_top10k on
identical features, the concordance is NOT isoform choice. If nonatlas_only
~= pas_counts, the result does not depend on the peaks being real poly(A) sites.

Output: results/pbmc_novelty/feature_spaces.tsv / .json
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
SEED = 42
ATLAS_WIN = 100
N_TOP = 10000
PARAMS = dict(resolution=1.0, n_components=50, n_dims=40, n_neighbors=30,
              scale_factor=1e4, depth_corr_threshold=0.75, random_seed=SEED)


def log(*a):
    print(*a, flush=True)


def tfidf_sparse(X, scale_factor=1e4):
    """Signac Method 1 TF-IDF, sparse. Exactly equivalent to PeakATail's dense
    implementation: the transform maps 0 -> log1p(0) = 0, so the sparsity
    pattern is preserved and only stored non-zeros need transforming."""
    X = sp.csr_matrix(X, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    cw = np.asarray((X > 0).sum(axis=0)).ravel().astype(np.float64)
    cw[cw == 0] = 1.0
    idf = X.shape[0] / cw
    Y = X.copy()
    Y.data = Y.data / np.repeat(tot, np.diff(Y.indptr))
    Y = Y.multiply(sp.csr_matrix(idf.reshape(1, -1)))
    Y = sp.csr_matrix(Y)
    Y.data = np.log1p(Y.data * scale_factor)
    return Y


def run_space(name, X, obs_names, counts_for_depth=None, prenormalised=False):
    """Identical PeakATail pipeline on an arbitrary feature matrix."""
    A = ad_mod.AnnData(sp.csr_matrix(X, dtype=np.float64) if sp.issparse(X)
                       else np.asarray(X, np.float64))
    A.obs_names = obs_names
    strat = get_strategy("leiden_tfidf", **PARAMS)
    if prenormalised:
        # usage fractions are already normalised; only the depth vector and
        # the LSI/Leiden half of the pipeline apply
        A.obs["total_counts"] = counts_for_depth
    else:
        A.obs["total_counts"] = (np.asarray(A.X.sum(axis=1)).ravel()
                                 if sp.issparse(A.X) else A.X.sum(axis=1))
        A.X = tfidf_sparse(A.X, PARAMS["scale_factor"])
    A = strat.reduce_dims(A)
    A = strat.cluster(A)
    lab = A.obs["leiden"].astype(str).values
    log(f"[{name}] n_features={X.shape[1]} n_clusters={len(set(lab))} "
        f"lsi_kept={len(A.uns['lsi_kept_components'])}")
    EMB[name] = A.obsm["X_lsi"]
    NEIGH[name] = A
    return lab


EMB = {}
NEIGH = {}
RES_GRID = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]


def sweep(name, ct):
    """Leiden resolution sweep on the already-built kNN graph. Cluster count is
    a free parameter of every one of these spaces; comparing AMI/ARI at a
    single resolution compares partitions of different granularity, and ARI in
    particular punishes over-splitting. This reports the whole curve and the
    k-matched point."""
    A = NEIGH[name]
    out = []
    for r in RES_GRID:
        sc.tl.leiden(A, resolution=r, key_added=f"_r{r}", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
        l = A.obs[f"_r{r}"].astype(str).values[ROWS]
        out.append({"space": name, "resolution": r,
                    "n_clusters": int(len(set(l))),
                    "AMI_vs_celltype": float(adjusted_mutual_info_score(l, ct)),
                    "ARI_vs_celltype": float(adjusted_rand_score(l, ct))})
    return out


# ---------------------------------------------------------------------------
log("[load]")
pas = sc.read_h5ad(PAS_CLUST)
with h5py.File(PAS_RAW, "r") as f:
    counts = sp.csr_matrix(
        (f["X/data"][:].astype(np.float64), f["X/indices"][:],
         f["X/indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
log(f"[load] {counts.shape} nnz={counts.nnz:,}")

M = pd.read_csv(OUT / "matched_cells.tsv", sep="\t")
atl = pd.read_csv(OUT / "pas_atlas_distance.tsv", sep="\t")
pas_id_arr = pas.var["pas_id"].values
d = pd.Series(atl["atlas_dist"].values,
              index=atl["pas_id"].values).reindex(pas_id_arr).values
d = np.where(d < 0, np.nan, d)
atlas_ok = np.nan_to_num(d, nan=1e9) <= ATLAS_WIN
log(f"[atlas] supported {atlas_ok.sum():,} / not {(~atlas_ok).sum():,}")

gene_codes = pas.var["gene_id"].cat.codes.values.astype(np.int32)
n_genes = len(pas.var["gene_id"].cat.categories)
n_pas_per_gene = np.bincount(gene_codes, minlength=n_genes)
multi = n_pas_per_gene[gene_codes] >= 2
G = sp.csr_matrix((np.ones(pas.n_vars), (gene_codes, np.arange(pas.n_vars))),
                  shape=(n_genes, pas.n_vars))

obs_names = pas.obs_names.values
cell_tot = np.asarray(counts.sum(axis=1)).ravel()

# ---- feature selection for the counts-vs-usage head-to-head ---------------
# The decisive comparison: the SAME features, once as counts and once as
# within-gene usage fractions. Everything here is vectorised; the only dense
# object ever built is the final n_cells x N_TOP usage matrix.
log("[select] top-variance within-gene usage features (multi-PAS genes)")
genec = np.asarray((counts @ G.T).todense(), np.float32)     # cells x genes
det = (genec >= 3)
gene_det_n = det.sum(axis=0)
gene_det = gene_det_n / counts.shape[0]
ok_gene = (gene_det >= 0.25) & (n_pas_per_gene >= 2)
cand = np.where(multi & ok_gene[gene_codes])[0]
log(f"[select] {int(ok_gene.sum())} genes detected (>=3 counts) in >=25% of "
    f"cells -> {cand.size} candidate PAS")

sub = counts[:, cand].tocoo()
cand_gene = gene_codes[cand]
tvals = genec[sub.row, cand_gene[sub.col]]
keep_nz = tvals >= 3                       # only cells where the gene is called
f_nz = np.zeros(sub.nnz, np.float64)
f_nz[keep_nz] = sub.data[keep_nz] / tvals[keep_nz]

n_obs = gene_det_n[cand_gene].astype(np.float64)             # per feature
sum_f = np.bincount(sub.col, weights=f_nz, minlength=cand.size)
sum_f2 = np.bincount(sub.col, weights=f_nz ** 2, minlength=cand.size)
mean_f = np.divide(sum_f, n_obs, out=np.zeros_like(sum_f), where=n_obs > 0)
var_f = np.divide(sum_f2, n_obs, out=np.zeros_like(sum_f2),
                  where=n_obs > 0) - mean_f ** 2
var_f[n_obs <= 50] = 0.0
top = np.argsort(-var_f)[:N_TOP]
sel = cand[top]
log(f"[select] kept {sel.size} features; usage variance range "
    f"{var_f[top].min():.2e}-{var_f[top].max():.2e}")

# dense usage matrix for the selected features only:
#   observed cell (gene called)  -> the real within-gene fraction
#   unobserved cell              -> the feature's mean over observed cells
#                                   (i.e. "no information", not "usage 0")
colmap = np.full(cand.size, -1, np.int64)
colmap[top] = np.arange(top.size)
Om = det[:, cand_gene[top]]                                   # n_cells x N_TOP
U = np.tile(mean_f[top].astype(np.float32), (counts.shape[0], 1))
U[Om] = 0.0
cc_all = colmap[sub.col]
m = keep_nz & (cc_all >= 0)
U[sub.row[m], cc_all[m]] = f_nz[m].astype(np.float32)
frac_observed = float(Om.mean())
log(f"[select] mean fraction of cells with the gene detected: "
    f"{frac_observed:.3f}")
OM = Om.copy()          # kept: the detection mask is itself a confounder
del genec, det, sub, tvals, f_nz, Om

# ---------------------------------------------------------------------------
M_ = pd.read_csv(OUT / "matched_cells.tsv", sep="\t")
_pos = pd.Series(np.arange(len(obs_names)), index=obs_names)
ROWS = _pos.loc[pas.obs_names.values[M_["pas_row"].values]].values
CT_ = M_["celltype"].values

labels = {"published": pas.obs["leiden"].astype(str).values}
labels["pas_counts"] = run_space("pas_counts", counts, obs_names)
labels["gene_level"] = run_space("gene_level", (counts @ G.T).tocsr(), obs_names)
labels["atlas_only"] = run_space("atlas_only", counts[:, atlas_ok], obs_names)
labels["nonatlas_only"] = run_space("nonatlas_only", counts[:, ~atlas_ok],
                                    obs_names)
labels["counts_top10k"] = run_space("counts_top10k", counts[:, sel], obs_names)
labels["usage_top10k"] = run_space("usage_top10k", U, obs_names,
                                   counts_for_depth=cell_tot,
                                   prenormalised=True)
# generous variant: column-centred + unit-variance, i.e. TruncatedSVD == PCA.
# Fractions are continuous, so this is the representation usage-space deserves;
# running it removes the "you crippled the usage space" objection.
Uz = U - U.mean(axis=0, keepdims=True)
sd = Uz.std(axis=0, keepdims=True)
sd[sd == 0] = 1.0
Uz /= sd
labels["usage_top10k_scaled"] = run_space("usage_top10k_scaled", Uz, obs_names,
                                          counts_for_depth=cell_tot,
                                          prenormalised=True)
del Uz

# ---- the control the usage space needs -----------------------------------
# U imputes "gene not detected" cells with the feature mean, so the PATTERN OF
# IMPUTATION is a gene-detection matrix -- i.e. expression -- and a clustering
# could ride on it without using any usage value at all. Two controls:
#   detection_mask  the binary mask alone, no usage values whatsoever
#   usage_shuffled  the same mask, but observed usage values permuted within
#                   each feature across cells (marginals and mask preserved,
#                   the cell<->usage link destroyed)
labels["detection_mask"] = run_space("detection_mask",
                                     OM.astype(np.float64), obs_names)
rng = np.random.default_rng(SEED)
Ush = U.copy()
for k in range(Ush.shape[1]):
    idx = np.where(OM[:, k])[0]
    if idx.size > 1:
        Ush[idx, k] = Ush[rng.permutation(idx), k]
labels["usage_shuffled"] = run_space("usage_shuffled", Ush, obs_names,
                                     counts_for_depth=cell_tot,
                                     prenormalised=True)
Ushz = Ush - Ush.mean(axis=0, keepdims=True)
sdh = Ushz.std(axis=0, keepdims=True)
sdh[sdh == 0] = 1.0
Ushz /= sdh
labels["usage_shuffled_scaled"] = run_space("usage_shuffled_scaled", Ushz,
                                            obs_names,
                                            counts_for_depth=cell_tot,
                                            prenormalised=True)
del Ush, Ushz

# ---------------------------------------------------------------------------
log("\n[score] vs marker-derived GEX cell types and GEX Leiden "
    "(matched cells only)")
M = M_
rows = ROWS
ct = CT_
gl = M["gex_leiden"].astype(str).values
pub = labels["published"][rows]

rec = []
for name, lab in labels.items():
    l = lab[rows]
    rec.append({
        "space": name,
        "n_features": {"published": pas.n_vars, "pas_counts": pas.n_vars,
                       "gene_level": n_genes,
                       "atlas_only": int(atlas_ok.sum()),
                       "nonatlas_only": int((~atlas_ok).sum()),
                       "counts_top10k": int(sel.size),
                       "usage_top10k": int(sel.size),
                       "usage_top10k_scaled": int(sel.size),
                       "detection_mask": int(sel.size),
                       "usage_shuffled": int(sel.size),
                       "usage_shuffled_scaled": int(sel.size)}[name],
        "n_clusters": int(len(set(lab))),
        "AMI_vs_celltype": float(adjusted_mutual_info_score(l, ct)),
        "ARI_vs_celltype": float(adjusted_rand_score(l, ct)),
        "AMI_vs_gexleiden": float(adjusted_mutual_info_score(l, gl)),
        "ARI_vs_gexleiden": float(adjusted_rand_score(l, gl)),
        "AMI_vs_published": float(adjusted_mutual_info_score(l, pub)),
        "ARI_vs_published": float(adjusted_rand_score(l, pub)),
    })
df = pd.DataFrame(rec)
log(df.to_string(index=False))
df.to_csv(OUT / "feature_spaces.tsv", sep="\t", index=False)

log("\n[sweep] Leiden resolution sweep per space (k-matched comparison)")
sw = []
for name in labels:
    if name in NEIGH:
        sw += sweep(name, ct)
SW = pd.DataFrame(sw)
SW.to_csv(OUT / "feature_spaces_resolution_sweep.tsv", sep="\t", index=False)
km = []
n_ct = len(set(ct))
for name, g in SW.groupby("space"):
    g = g.copy()
    g["dk"] = (g["n_clusters"] - n_ct).abs()
    b = g.sort_values(["dk", "resolution"]).iloc[0]
    km.append({"space": name, "resolution": b["resolution"],
               "n_clusters": int(b["n_clusters"]),
               "AMI_vs_celltype_kmatched": b["AMI_vs_celltype"],
               "ARI_vs_celltype_kmatched": b["ARI_vs_celltype"],
               "AMI_best_over_sweep": g["AMI_vs_celltype"].max(),
               "ARI_best_over_sweep": g["ARI_vs_celltype"].max()})
KM = pd.DataFrame(km).set_index("space").reindex(
    [s for s in df["space"] if s in set(SW["space"])])
log(KM.to_string())
KM.to_csv(OUT / "feature_spaces_kmatched.tsv", sep="\t")
pd.DataFrame({k: v for k, v in labels.items()},
             index=obs_names).to_csv(OUT / "feature_space_labels.tsv", sep="\t")
with open(OUT / "feature_spaces.json", "w") as fh:
    json.dump({"table": rec, "params": PARAMS, "N_TOP": N_TOP,
               "ATLAS_WIN": ATLAS_WIN,
               "usage_frac_cells_observed": frac_observed,
               "n_atlas_supported": int(atlas_ok.sum()),
               "n_not_atlas_supported": int((~atlas_ok).sum())}, fh, indent=2)
log(f"[out] {OUT/'feature_spaces.tsv'}")
log("[done]")
