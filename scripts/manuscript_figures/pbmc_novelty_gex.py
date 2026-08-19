#!/usr/bin/env python3
"""
pbmc_novelty_gex.py -- STAGE 1 of the pbmc_novelty analysis.

Standard scanpy gene-expression pipeline on the public 10x pbmc_10k_v3
filtered_feature_bc_matrix, plus MARKER-DERIVED cell-type labels.

These labels are NOT curated ground truth. They are the modal assignment of a
canonical PBMC marker panel scored per Leiden cluster; they inherit every error
of the marker panel and of the GEX clustering itself. Everything downstream that
calls them "cell types" means "marker-derived label".

Output: results/pbmc_novelty/gex_labeled.h5ad
"""
import os
os.environ.setdefault("LC_ALL", "C")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.n_jobs = 8

GEX_DIR = "/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/filtered_feature_bc_matrix"
OUT = Path("/mnt/ssd1/Projects/PeakATail_wd/results/pbmc_novelty")
OUT.mkdir(parents=True, exist_ok=True)
H5AD = OUT / "gex_labeled.h5ad"
SEED = 42

# ---------------------------------------------------------------------------
# canonical PBMC marker panel (the sets named in the task, lightly extended)
# lineage level first, then two sub-splits.  Every gene is a SYMBOL; the 10x
# features.tsv carries both ENSG and symbol so we map on symbol.
# ---------------------------------------------------------------------------
LINEAGE = {
    "T":          ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"],
    "NK":         ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"],
    "B":          ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Mono":       ["LYZ", "CD14", "S100A8", "S100A9", "VCAN", "FCGR3A", "MS4A7"],
    # NOTE: CST3 is in the task's DC set but is pan-myeloid (mean log-norm 2.9-3.3
    # in every monocyte cluster here). Including it mislabels the FCGR3A/MS4A7/
    # CDKN1C non-classical monocyte cluster as DC. It is scored and reported
    # (sc_DC_withCST3) but excluded from the set that drives the call.
    "DC":         ["FCER1A", "CLEC10A", "CD1C", "CLEC9A", "LILRA4", "IL3RA"],
    "Platelet":   ["PPBP", "PF4", "ITGA2B", "TUBB1"],
    "Erythroid":  ["HBB", "HBA1", "HBA2", "ALAS2"],
}
EXTRA = {  # scored for the record, never used for the call
    "DC_withCST3": ["FCER1A", "CST3", "CLEC10A", "CD1C", "CLEC9A", "LILRA4"],
}
T_SUB = {
    "CD4-T": ["IL7R", "CD4", "CCR7", "LTB", "TCF7"],
    "CD8-T": ["CD8A", "CD8B", "GZMK", "CCL5", "CD8B2"],
}
MONO_SUB = {
    "CD14-Mono":  ["CD14", "S100A8", "S100A9", "VCAN", "LYZ"],
    "FCGR3A-Mono": ["FCGR3A", "MS4A7", "CDKN1C", "LST1", "AIF1"],
}

# CD3 gate. NKG7/GNLY/KLRD1 are shared by NK cells and cytotoxic CD8 T cells,
# so an NK-vs-T argmax on those genes alone calls CD3D/CD3E/TRAC-high effector
# CD8 T cells "NK". The classical discriminator is CD3 surface expression:
# NK cells are CD3-negative. Applied as a hard gate on the per-cluster mean
# log-normalised CD3 module.
CD3_MODULE = ["CD3D", "CD3E", "CD3G", "TRAC"]
CD3_GATE = 1.0        # mean log-norm CD3 module above this => T lineage
NK_MODULE = ["NKG7", "GNLY", "KLRD1", "KLRF1"]
NK_GATE = 1.0         # ... below the CD3 gate, this much NK module => NK


def log(*a):
    print(*a, flush=True)


def main():
    log("[gex] reading 10x matrix")
    ad = sc.read_10x_mtx(GEX_DIR, var_names="gene_symbols", cache=False)
    ad.var_names_make_unique()
    log(f"[gex] raw {ad.shape[0]} cells x {ad.shape[1]} genes")
    ad.layers["counts_raw"] = ad.X.copy()

    # ---- QC ---------------------------------------------------------------
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    ad.var["ribo"] = ad.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt", "ribo"], percent_top=None,
                               log1p=False, inplace=True)
    n0 = ad.n_obs
    keep = (ad.obs["n_genes_by_counts"] >= 200) & (ad.obs["pct_counts_mt"] < 20)
    ad.obs["qc_pass"] = keep.values
    qc_drop = int((~keep).sum())
    log(f"[gex] QC: {qc_drop}/{n0} cells fail (n_genes<200 or pct_mt>=20)")
    ad = ad[keep].copy()
    sc.pp.filter_genes(ad, min_cells=3)
    log(f"[gex] after QC {ad.shape[0]} cells x {ad.shape[1]} genes")

    # ---- doublet score (scores only; scanpy's scrublet needs an explicit
    #      threshold because skimage is absent -- we use the SCORE, not the call)
    try:
        sc.pp.scrublet(ad, threshold=0.25, random_state=SEED, verbose=False)
        log("[gex] scrublet doublet scores computed")
    except Exception as e:          # pragma: no cover
        log(f"[gex] scrublet failed ({e}); filling NaN")
        ad.obs["doublet_score"] = np.nan

    # ---- normalise / HVG / PCA / neighbours / leiden ----------------------
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw = ad                      # keep full log-normalised for scoring
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
    lognorm = ad.X.copy()            # all genes, log-normalised
    ad = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=50, svd_solver="arpack", random_state=SEED)
    # match the PAS side's graph settings (n_neighbors=30, 40 comps, res 1.0)
    sc.pp.neighbors(ad, n_neighbors=30, n_pcs=40, random_state=SEED)
    sc.tl.leiden(ad, resolution=1.0, key_added="leiden",
                 random_state=SEED, flavor="igraph", n_iterations=2,
                 directed=False)
    sc.tl.umap(ad, random_state=SEED)
    log(f"[gex] leiden: {ad.obs['leiden'].nunique()} clusters")

    # ---- marker scoring ---------------------------------------------------
    # score on the FULL log-normalised matrix (ad.raw), not the HVG subset
    full = ad.raw.to_adata()
    full.obs = ad.obs.copy()

    def score_sets(sets, prefix):
        present = {}
        cols = []
        for name, genes in sets.items():
            g = [x for x in genes if x in full.var_names]
            present[name] = g
            key = f"{prefix}{name}"
            if not g:
                full.obs[key] = np.nan
            else:
                sc.tl.score_genes(full, g, score_name=key, random_state=SEED)
            cols.append(key)
        return cols, present

    lin_cols, lin_present = score_sets(LINEAGE, "sc_")
    t_cols, t_present = score_sets(T_SUB, "sc_")
    m_cols, m_present = score_sets(MONO_SUB, "sc_")
    x_cols, _ = score_sets(EXTRA, "sc_")

    S = full.obs[lin_cols + t_cols + m_cols + x_cols].copy()
    # z-score each set across cells so sets of different size are comparable
    Z = (S - S.mean()) / S.std(ddof=0)
    Z["leiden"] = ad.obs["leiden"].values
    per_cluster = Z.groupby("leiden", observed=True).mean()

    # per-cluster mean log-norm expression of the gate modules
    import scipy.sparse as sp

    def module_mean(genes):
        g = [x for x in genes if x in full.var_names]
        M = full[:, g].X
        M = M.toarray() if sp.issparse(M) else np.asarray(M)
        s = pd.Series(M.mean(axis=1), index=full.obs_names)
        return s.groupby(ad.obs["leiden"].values, observed=True).mean()

    cd3 = module_mean(CD3_MODULE)
    nkm = module_mean(NK_MODULE)

    lineage_call = per_cluster[lin_cols].idxmax(axis=1).str.replace("sc_", "",
                                                                   regex=False)
    gated = {}
    label = {}
    for cl, lin in lineage_call.items():
        gate_note = ""
        if lin == "NK" and cd3[cl] >= CD3_GATE:
            lin, gate_note = "T", f"NK->T (CD3 module {cd3[cl]:.2f})"
        elif lin == "T" and cd3[cl] < CD3_GATE and nkm[cl] >= NK_GATE:
            lin, gate_note = "NK", f"T->NK (CD3 {cd3[cl]:.2f}, NK {nkm[cl]:.2f})"
        gated[cl] = gate_note
        if lin == "T":
            label[cl] = per_cluster.loc[cl, t_cols].idxmax().replace("sc_", "")
        elif lin == "Mono":
            label[cl] = per_cluster.loc[cl, m_cols].idxmax().replace("sc_", "")
        else:
            label[cl] = lin
    ad.obs["celltype"] = ad.obs["leiden"].map(label).astype("category")
    ad.obs["lineage"] = ad.obs["leiden"].map(
        {c: (label[c].split("-")[-1] if "-" in label[c] else label[c])
         for c in label}).astype("category")
    for c in lin_cols + t_cols + m_cols + x_cols:
        ad.obs[c] = full.obs[c].values

    # ---- QC flags per GEX cluster (doublet / low-depth), reported not dropped
    dbl = ad.obs.groupby("leiden", observed=True)["doublet_score"].median()
    umi = ad.obs.groupby("leiden", observed=True)["total_counts"].median()
    dbl_flag = dbl > 0.15
    lowdepth_flag = umi < 0.5 * float(umi.median())
    ad.obs["gex_cluster_doublet_suspect"] = ad.obs["leiden"].map(dbl_flag).values
    ad.obs["gex_cluster_low_depth"] = ad.obs["leiden"].map(lowdepth_flag).values

    log("[gex] per-cluster marker z-scores (mean over cells):")
    tab = per_cluster.round(3).copy()
    tab.insert(0, "n_cells", ad.obs["leiden"].value_counts().reindex(tab.index))
    tab.insert(1, "CD3mod", cd3.round(2))
    tab.insert(2, "NKmod", nkm.round(2))
    tab.insert(3, "med_dblscore", dbl.round(3))
    tab.insert(4, "med_nUMI", umi.astype(int))
    tab["LABEL"] = pd.Series(label)
    tab["gate"] = pd.Series(gated)
    tab["dbl_suspect"] = dbl_flag
    tab["low_depth"] = lowdepth_flag
    log(tab.to_string())
    log("[gex] label composition:")
    log(ad.obs["celltype"].value_counts().to_string())
    per_cluster = tab

    # ---- persist ----------------------------------------------------------
    ad.uns["marker_sets_used"] = json.dumps(
        {"lineage": lin_present, "T_sub": t_present, "Mono_sub": m_present})
    ad.uns["qc"] = json.dumps({"n_cells_input": int(n0), "n_cells_qc_fail":
                               qc_drop, "n_cells_kept": int(ad.n_obs)})
    # store the full log-normalised matrix as a layer-free companion for the
    # novelty stage (gene-level totals) -- write separately to keep h5ad lean
    ad.write(H5AD)
    log(f"[gex] wrote {H5AD}")

    per_cluster.to_csv(OUT / "gex_cluster_marker_scores.tsv", sep="\t")
    # log-normalised full matrix for the APA test
    import scipy.sparse as sp
    sp.save_npz(OUT / "gex_lognorm.npz", sp.csr_matrix(lognorm))
    pd.Series(full.var_names).to_csv(OUT / "gex_lognorm_genes.tsv", sep="\t",
                                     index=False, header=["gene_symbol"])
    log("[gex] done")


if __name__ == "__main__":
    main()
