#!/usr/bin/env python3
"""
pbmc_novelty_usage_features.py -- STAGE 5.

If the within-gene USAGE space turns out to carry cell-type identity, the very
next question is whether the sites carrying it are real poly(A) sites. This
redoes the (deterministic) feature selection of pbmc_novelty_spaces.py and
reports, for the 10,000 selected usage features: atlas support, the genes they
belong to, and how atlas support compares with the matched background.
"""
import os
os.environ.setdefault("LC_ALL", "C")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import fisher_exact

ROOT = Path("/mnt/ssd1/Projects/PeakATail_wd")
PAS_CLUST = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                    "07_clustering/default/clusters.h5ad")
PAS_RAW = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                  "06_preprocessing/default/preprocessed.h5ad")
PAS_BED = ROOT / "results/benchmark_tools/pbmc_10k_v3/peakatail/run/annotatedpas.bed"
OUT = ROOT / "results/pbmc_novelty"
N_TOP, ATLAS_WIN = 10000, 100

pas = sc.read_h5ad(PAS_CLUST)
with h5py.File(PAS_RAW, "r") as f:
    counts = sp.csr_matrix((f["X/data"][:].astype(np.float64),
                            f["X/indices"][:], f["X/indptr"][:]),
                           shape=tuple(f["X"].attrs["shape"]))
gene_codes = pas.var["gene_id"].cat.codes.values.astype(np.int32)
gene_cats = np.array(pas.var["gene_id"].cat.categories)
n_genes = len(gene_cats)
n_pas_per_gene = np.bincount(gene_codes, minlength=n_genes)
multi = n_pas_per_gene[gene_codes] >= 2
G = sp.csr_matrix((np.ones(pas.n_vars), (gene_codes, np.arange(pas.n_vars))),
                  shape=(n_genes, pas.n_vars))
pas_id_arr = pas.var["pas_id"].values

atl = pd.read_csv(OUT / "pas_atlas_distance.tsv", sep="\t")
d = pd.Series(atl["atlas_dist"].values,
              index=atl["pas_id"].values).reindex(pas_id_arr).values
atlas_ok = np.nan_to_num(np.where(d < 0, np.nan, d), nan=1e9) <= ATLAS_WIN

bed = pd.read_csv(PAS_BED, sep="\t", header=None, low_memory=False,
                  names=["chrom", "start", "end", "pas_id", "gene_id",
                         "symbol", "strand", "score", "tier"])
sym = bed.drop_duplicates("gene_id").set_index("gene_id")["symbol"].to_dict()

genec = np.asarray((counts @ G.T).todense(), np.float32)
det = genec >= 3
gene_det_n = det.sum(axis=0)
ok_gene = (gene_det_n / counts.shape[0] >= 0.25) & (n_pas_per_gene >= 2)
cand = np.where(multi & ok_gene[gene_codes])[0]
sub = counts[:, cand].tocoo()
cand_gene = gene_codes[cand]
tvals = genec[sub.row, cand_gene[sub.col]]
keep = tvals >= 3
f_nz = np.zeros(sub.nnz)
f_nz[keep] = sub.data[keep] / tvals[keep]
n_obs = gene_det_n[cand_gene].astype(float)
s1 = np.bincount(sub.col, weights=f_nz, minlength=cand.size)
s2 = np.bincount(sub.col, weights=f_nz ** 2, minlength=cand.size)
mean_f = np.divide(s1, n_obs, out=np.zeros_like(s1), where=n_obs > 0)
var_f = np.divide(s2, n_obs, out=np.zeros_like(s2), where=n_obs > 0) - mean_f**2
var_f[n_obs <= 50] = 0.0
top = np.argsort(-var_f)[:N_TOP]
sel = cand[top]

a_sel = int(atlas_ok[sel].sum())
a_cand = int(atlas_ok[cand].sum())
a_all = int(atlas_ok.sum())
orv, p = fisher_exact([[a_sel, len(sel) - a_sel],
                       [a_cand - a_sel, (len(cand) - len(sel)) - (a_cand - a_sel)]])
orv2, p2 = fisher_exact([[a_sel, len(sel) - a_sel],
                         [a_all, pas.n_vars - a_all]])
res = {
    "n_selected": int(len(sel)),
    "n_candidate_universe": int(len(cand)),
    "n_all_pas": int(pas.n_vars),
    "atlas_frac_selected": a_sel / len(sel),
    "atlas_frac_candidate_universe": a_cand / len(cand),
    "atlas_frac_all_pas": a_all / pas.n_vars,
    "fisher_OR_selected_vs_candidates": float(orv),
    "fisher_p_selected_vs_candidates": float(p),
    "fisher_OR_selected_vs_all": float(orv2),
    "fisher_p_selected_vs_all": float(p2),
    "n_genes_selected": int(len(set(gene_codes[sel]))),
}
print(json.dumps(res, indent=2))
tab = pd.DataFrame({
    "pas_col": sel, "pas_id": pas_id_arr[sel],
    "gene_id": gene_cats[gene_codes[sel]],
    "gene_symbol": [sym.get(g, "") for g in gene_cats[gene_codes[sel]]],
    "usage_variance": var_f[top], "mean_usage": mean_f[top],
    "atlas_dist": d[sel], "atlas_supported": atlas_ok[sel],
    "n_pas_in_gene": n_pas_per_gene[gene_codes[sel]],
})
tab.sort_values("usage_variance", ascending=False).to_csv(
    OUT / "usage_features.tsv", sep="\t", index=False)
with open(OUT / "usage_features_summary.json", "w") as fh:
    json.dump(res, fh, indent=2)
print("top 25 usage features:")
print(tab.head(25)[["gene_symbol", "usage_variance", "atlas_dist",
                    "atlas_supported", "n_pas_in_gene"]].to_string(index=False))
