#!/usr/bin/env python3
"""
pbmc_novelty_analysis.py -- STAGE 2 of the pbmc_novelty analysis.

Barcode matching, PAS-vs-GEX clustering concordance, and the NOVELTY TEST:
are there populations resolved in poly(A)-site space that gene-expression
clustering does not resolve, and if so is the split driven by APA (within-gene
site choice) rather than by depth / ambient / doublets / expression level?

All intermediate results are written to results/pbmc_novelty/ as tsv/json so
the figure script never recomputes anything.

Pre-registered criteria (fixed before looking at results, printed on the figure)
  SPLIT       a marker-derived GEX type T is "split" in PAS space if >=2 PAS
              clusters are each >=60% type-T by composition and each holds
              >=10% of T's cells and >=100 cells.
  CANDIDATE   within a split (A vs B, the two largest such PAS clusters), a PAS
              p of gene g is a candidate if pseudobulk within-gene usage differs
              by |dProp| >= 0.20, with gene pseudobulk >= 100 PAS counts in each
              arm and >= 2 PAS in the gene.
  SIGNIFICANT per-CELL usage of p (cells with >=5 counts on g, >=30 such cells
              per arm) differs between arms, Mann-Whitney, BH q < 0.05.
              Cells -- not reads -- are the replicate unit; a read-level Fisher
              test on the same data is anti-conservative (see fdr_calibration).
  APA-DRIVEN  significant AND the gene's TOTAL expression measured independently
              in the 10x GEX matrix is comparable across the arms,
              |log2FC(mean CP10K)| <= 0.5.
  NULL        the entire pipeline is re-run on 20 random relabellings of the
              same cells within the split -> empirical false-positive count.
"""
import os
os.environ.setdefault("LC_ALL", "C")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

import json
from pathlib import Path

import anndata as ad_mod
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

sc.settings.n_jobs = 8

ROOT = Path("/mnt/ssd1/Projects/PeakATail_wd")
PAS_CLUST = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                    "07_clustering/default/clusters.h5ad")
PAS_RAW = ROOT / ("results/benchmark_tools/pbmc_10k_v3/peakatail/run/"
                  "06_preprocessing/default/preprocessed.h5ad")
PAS_BED = ROOT / "results/benchmark_tools/pbmc_10k_v3/peakatail/run/annotatedpas.bed"
OUT = ROOT / "results/pbmc_novelty"
GEX_H5AD = OUT / "gex_labeled.h5ad"

SEED = 42
N_PERM = 20

# ---- pre-registered thresholds --------------------------------------------
SPLIT_MIN_PURITY = 0.60
SPLIT_MIN_FRAC = 0.10
SPLIT_MIN_CELLS = 100
SPLIT_MIN_TYPE_CELLS = 300
CAND_MIN_DPROP = 0.20
CAND_MIN_GENE_COUNTS = 100
CELL_MIN_GENE_COUNTS = 5
CELL_MIN_N = 30
FDR_Q = 0.05
EXPR_MAX_ABS_LOG2FC = 0.5

# Laughney reference (results/figures/manuscript/clustering_concordance.tsv)
LAUGHNEY = {
    "median_AMI_celltype_vs_pas": 0.6616,
    "median_ARI_celltype_vs_pas": 0.4627,
    "median_AMI_gexleiden_vs_pas": 0.6565,
    "median_ARI_gexleiden_vs_pas": 0.4564,
    "min_AMI_celltype_vs_pas": 0.3667,
    "max_AMI_celltype_vs_pas": 0.7795,
    "min_ARI_celltype_vs_pas": 0.2002,
    "max_ARI_celltype_vs_pas": 0.6705,
    "n_samples": 17,
}


def log(*a):
    print(*a, flush=True)


def bh(p):
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


# ===========================================================================
# 1. load
# ===========================================================================
log("[load] PAS clusters")
pas = sc.read_h5ad(PAS_CLUST)                 # X = TF-IDF, obsm X_lsi/X_umap
log(f"[load] PAS {pas.shape}")
log("[load] PAS raw counts")
with h5py.File(PAS_RAW, "r") as f:
    raw = sp.csr_matrix(
        (f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
        shape=tuple(f["X"].attrs["shape"]))
    raw_obs = np.array([x.decode() for x in f["obs/_index"][:]])
    raw_var = np.array([x.decode() for x in f["var/_index"][:]])
assert (raw_obs == pas.obs_names.values).all(), "raw/clusters obs order differs"
assert (raw_var == pas.var_names.values).all(), "raw/clusters var order differs"
pas.layers["counts"] = raw
log(f"[load] raw counts total = {raw.sum():,.0f}")

log("[load] GEX labelled")
gex = sc.read_h5ad(GEX_H5AD)
lognorm = sp.load_npz(OUT / "gex_lognorm.npz").tocsc()
lognorm_genes = pd.read_csv(OUT / "gex_lognorm_genes.tsv", sep="\t")[
    "gene_symbol"].values
assert lognorm.shape[0] == gex.n_obs
gene_pos = {g: i for i, g in enumerate(lognorm_genes)}

# PAS gene annotation
bed = pd.read_csv(PAS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "pas_id", "gene_id",
                         "symbol", "strand", "score", "tier"])
sym = bed.drop_duplicates("gene_id").set_index("gene_id")["symbol"].to_dict()
pas_info = bed.drop_duplicates("pas_id").set_index("pas_id")
pas.var["gene_symbol"] = pas.var["gene_id"].map(sym).astype(str)

# ===========================================================================
# 2. barcode matching
# ===========================================================================
log("\n[match] barcode reconciliation")
pas_bc_raw = pas.obs["barcode"].astype(str).values
prefixes = {b.split("_")[0] for b in pas_bc_raw}
log(f"[match] PAS barcode prefixes: {prefixes}")
pas_bc = np.array([b.split("_", 1)[1] if "_" in b else b for b in pas_bc_raw])
assert all(len(b) == 16 for b in pas_bc[:100]), "unexpected CB length"
pas_bc_cr = np.array([b + "-1" for b in pas_bc])   # restore CellRanger suffix
pas.obs["cb"] = pas_bc_cr

gex_bc = gex.obs_names.values
shared = np.intersect1d(pas_bc_cr, gex_bc)
n_pas, n_gex, n_sh = len(pas_bc_cr), len(gex_bc), len(shared)
log(f"[match] PAS cells {n_pas} | GEX cells (post-QC) {n_gex} | shared {n_sh}")
log(f"[match] coverage of PAS by GEX = {n_sh/n_pas:.4f}; "
    f"of GEX by PAS = {n_sh/n_gex:.4f}")
assert len(set(pas_bc_cr)) == n_pas, "duplicate PAS barcodes"

# raw (pre-QC) 10x barcode list, to separate 'absent' from 'QC-dropped'
raw_gex_bc = pd.read_csv(
    ROOT / "data/benchmark/pbmc_10k_v3/filtered_feature_bc_matrix/barcodes.tsv.gz",
    header=None)[0].values
n_in_raw = int(np.isin(pas_bc_cr, raw_gex_bc).sum())
log(f"[match] PAS cells present in the 10x filtered list (pre-QC): {n_in_raw} "
    f"({n_in_raw/n_pas:.4f}); lost to GEX QC: {n_in_raw - n_sh}; "
    f"never in the 10x call set: {n_pas - n_in_raw}")

pi = pd.Series(np.arange(n_pas), index=pas_bc_cr).loc[shared].values
gi = pd.Series(np.arange(n_gex), index=gex_bc).loc[shared].values

M = pd.DataFrame({
    "cb": shared,
    "pas_leiden": pas.obs["leiden"].values[pi].astype(str),
    "gex_leiden": gex.obs["leiden"].values[gi].astype(str),
    "celltype": gex.obs["celltype"].values[gi].astype(str),
    "pas_counts": pas.obs["total_counts"].values[pi],
    "pas_n_genes": pas.obs["n_genes"].values[pi],
    "gex_counts": gex.obs["total_counts"].values[gi],
    "gex_n_genes": gex.obs["n_genes_by_counts"].values[gi],
    "pct_mt": gex.obs["pct_counts_mt"].values[gi],
    "doublet_score": gex.obs["doublet_score"].values[gi],
    "dbl_cluster": gex.obs["gex_cluster_doublet_suspect"].values[gi],
    "lowdepth_cluster": gex.obs["gex_cluster_low_depth"].values[gi],
    "pas_umap1": pas.obsm["X_umap"][pi, 0],
    "pas_umap2": pas.obsm["X_umap"][pi, 1],
    "gex_umap1": gex.obsm["X_umap"][gi, 0],
    "gex_umap2": gex.obsm["X_umap"][gi, 1],
})
M["pas_row"] = pi
M["gex_row"] = gi
M.to_csv(OUT / "matched_cells.tsv", sep="\t", index=False)

# ---- does the PAS LSI space encode sequencing depth? ----------------------
lsi = pas.obsm["X_lsi"]
depth = np.log10(pas.obs["total_counts"].values.astype(float))
lsi_depth_r = np.array([np.corrcoef(lsi[:, k], depth)[0, 1]
                        for k in range(lsi.shape[1])])
log(f"[qc] |corr(LSI_k, log10 depth)| max = {np.abs(lsi_depth_r).max():.3f} "
    f"at component {int(np.abs(lsi_depth_r).argmax())}; "
    f"LSI1 r = {lsi_depth_r[0]:.3f}; "
    f"n components with |r|>0.5: {int((np.abs(lsi_depth_r) > 0.5).sum())}")

# ===========================================================================
# 3. concordance
# ===========================================================================
log("\n[concordance] on shared cells")
conc = {}
for name, a, b in [
    ("pas_vs_gexleiden", M["pas_leiden"], M["gex_leiden"]),
    ("pas_vs_celltype", M["pas_leiden"], M["celltype"]),
    ("gexleiden_vs_celltype", M["gex_leiden"], M["celltype"]),
]:
    conc[f"AMI_{name}"] = float(adjusted_mutual_info_score(a, b))
    conc[f"ARI_{name}"] = float(adjusted_rand_score(a, b))
conc["n_shared_cells"] = int(n_sh)
conc["n_pas_cells"] = int(n_pas)
conc["n_gex_cells_postqc"] = int(n_gex)
conc["n_pas_clusters"] = int(M["pas_leiden"].nunique())
conc["n_gex_clusters"] = int(M["gex_leiden"].nunique())
conc["n_celltypes"] = int(M["celltype"].nunique())
conc["coverage_pas_by_gex"] = float(n_sh / n_pas)
conc["max_abs_corr_lsi_depth"] = float(np.abs(lsi_depth_r).max())
conc["corr_lsi1_depth"] = float(lsi_depth_r[0])

# sensitivity: drop the doublet-suspect and low-depth GEX clusters
clean = M[~(M["dbl_cluster"] | M["lowdepth_cluster"])]
conc["n_shared_cells_clean"] = int(len(clean))
conc["AMI_pas_vs_celltype_clean"] = float(
    adjusted_mutual_info_score(clean["pas_leiden"], clean["celltype"]))
conc["ARI_pas_vs_celltype_clean"] = float(
    adjusted_rand_score(clean["pas_leiden"], clean["celltype"]))
conc["AMI_pas_vs_gexleiden_clean"] = float(
    adjusted_mutual_info_score(clean["pas_leiden"], clean["gex_leiden"]))
conc["ARI_pas_vs_gexleiden_clean"] = float(
    adjusted_rand_score(clean["pas_leiden"], clean["gex_leiden"]))
for k, v in conc.items():
    log(f"  {k:38s} {v}")
log("  Laughney medians: AMI(celltype) %.4f  ARI(celltype) %.4f  "
    "AMI(gexleiden) %.4f  ARI(gexleiden) %.4f" % (
        LAUGHNEY["median_AMI_celltype_vs_pas"],
        LAUGHNEY["median_ARI_celltype_vs_pas"],
        LAUGHNEY["median_AMI_gexleiden_vs_pas"],
        LAUGHNEY["median_ARI_gexleiden_vs_pas"]))

# ===========================================================================
# 4. contingency + split detection
# ===========================================================================
ct = pd.crosstab(M["pas_leiden"], M["celltype"])
ct = ct.loc[sorted(ct.index, key=int)]
ct.to_csv(OUT / "contingency_pas_x_celltype.tsv", sep="\t")
log("\n[contingency] PAS cluster x marker-derived GEX type")
log(ct.to_string())
purity = ct.div(ct.sum(axis=1), axis=0)          # composition of each PAS clust
share = ct.div(ct.sum(axis=0), axis=1)           # share of each type per clust

splits = []
for t in ct.columns:
    n_t = int(ct[t].sum())
    if n_t < SPLIT_MIN_TYPE_CELLS:
        continue
    ok = [c for c in ct.index
          if purity.loc[c, t] >= SPLIT_MIN_PURITY
          and share.loc[c, t] >= SPLIT_MIN_FRAC
          and ct.loc[c, t] >= SPLIT_MIN_CELLS]
    if len(ok) >= 2:
        ok = sorted(ok, key=lambda c: -ct.loc[c, t])
        splits.append({"celltype": t, "n_type_cells": n_t,
                       "pas_clusters": ok,
                       "sizes": [int(ct.loc[c, t]) for c in ok]})
log("\n[splits] GEX types split across >=2 type-dominated PAS clusters:")
for s in splits:
    log(f"  {s['celltype']:12s} n={s['n_type_cells']:5d} -> PAS clusters "
        f"{s['pas_clusters']} sizes {s['sizes']}")
if not splits:
    log("  (none)")

# ===========================================================================
# 5. the APA test
# ===========================================================================
counts = pas.layers["counts"].tocsc()
gene_codes = pas.var["gene_id"].cat.codes.values
gene_cats = np.array(pas.var["gene_id"].cat.categories)
n_pas_per_gene = np.bincount(gene_codes, minlength=len(gene_cats))
multi = n_pas_per_gene[gene_codes] >= 2
log(f"\n[apa] {int(multi.sum())}/{pas.n_vars} PAS lie in multi-PAS genes "
    f"({int((n_pas_per_gene >= 2).sum())}/{len(gene_cats)} genes)")

# gene x PAS incidence for fast pseudobulk gene totals
G = sp.csr_matrix((np.ones(pas.n_vars), (gene_codes, np.arange(pas.n_vars))),
                  shape=(len(gene_cats), pas.n_vars))

cp10k = None  # lazily built GEX CP10K means


def gex_mean_cp10k(rows, gene_symbol):
    j = gene_pos.get(gene_symbol)
    if j is None:
        return np.nan
    v = np.asarray(lognorm[rows, j].todense()).ravel()
    return float(np.expm1(v).mean())


def screen(rowsA, rowsB):
    """Pseudobulk within-gene usage screen. Returns candidate DataFrame."""
    sA = np.asarray(counts[rowsA].sum(axis=0)).ravel()
    sB = np.asarray(counts[rowsB].sum(axis=0)).ravel()
    gA = G @ sA
    gB = G @ sB
    keep_gene = (gA >= CAND_MIN_GENE_COUNTS) & (gB >= CAND_MIN_GENE_COUNTS) & \
                (n_pas_per_gene >= 2)
    m = multi & keep_gene[gene_codes]
    idx = np.where(m)[0]
    if idx.size == 0:
        return pd.DataFrame()
    pA = sA[idx] / gA[gene_codes[idx]]
    pB = sB[idx] / gB[gene_codes[idx]]
    d = pA - pB
    sel = np.abs(d) >= CAND_MIN_DPROP
    return pd.DataFrame({
        "pas_col": idx[sel], "gene_code": gene_codes[idx][sel],
        "propA": pA[sel], "propB": pB[sel], "dprop": d[sel],
        "countsA": sA[idx][sel], "countsB": sB[idx][sel],
        "gene_countsA": gA[gene_codes[idx]][sel],
        "gene_countsB": gB[gene_codes[idx]][sel],
    })


def percell_test(cand, rowsA, rowsB):
    """Per-cell usage test; cells are the replicate unit."""
    if cand.empty:
        return cand.assign(p=[], q=[], nA=[], nB=[])
    ps, nAs, nBs, mA, mB = [], [], [], [], []
    gsum_cache = {}
    for _, r in cand.iterrows():
        gc = int(r["gene_code"])
        if gc not in gsum_cache:
            cols = np.where(gene_codes == gc)[0]
            gsum_cache[gc] = (
                np.asarray(counts[rowsA][:, cols].sum(axis=1)).ravel(),
                np.asarray(counts[rowsB][:, cols].sum(axis=1)).ravel())
        tA, tB = gsum_cache[gc]
        j = int(r["pas_col"])
        cA = np.asarray(counts[rowsA][:, j].todense()).ravel()
        cB = np.asarray(counts[rowsB][:, j].todense()).ravel()
        okA = tA >= CELL_MIN_GENE_COUNTS
        okB = tB >= CELL_MIN_GENE_COUNTS
        nA, nB = int(okA.sum()), int(okB.sum())
        nAs.append(nA)
        nBs.append(nB)
        if nA < CELL_MIN_N or nB < CELL_MIN_N:
            ps.append(np.nan)
            mA.append(np.nan)
            mB.append(np.nan)
            continue
        fA = cA[okA] / tA[okA]
        fB = cB[okB] / tB[okB]
        mA.append(float(fA.mean()))
        mB.append(float(fB.mean()))
        try:
            ps.append(float(mannwhitneyu(fA, fB, alternative="two-sided")[1]))
        except ValueError:
            ps.append(np.nan)
    out = cand.copy()
    out["nA_cells"] = nAs
    out["nB_cells"] = nBs
    out["cell_propA"] = mA
    out["cell_propB"] = mB
    out["p"] = ps
    v = ~out["p"].isna()
    out["q"] = np.nan
    if v.any():
        out.loc[v, "q"] = bh(out.loc[v, "p"].values)
    return out


results = {}
rng = np.random.default_rng(SEED)
for s in splits:
    t = s["celltype"]
    cA_id, cB_id = s["pas_clusters"][0], s["pas_clusters"][1]
    selA = (M["celltype"] == t) & (M["pas_leiden"] == cA_id)
    selB = (M["celltype"] == t) & (M["pas_leiden"] == cB_id)
    rowsA = M.loc[selA, "pas_row"].values
    rowsB = M.loc[selB, "pas_row"].values
    gA = M.loc[selA, "gex_row"].values
    gB = M.loc[selB, "gex_row"].values
    tag = f"{t}_PAS{cA_id}_vs_PAS{cB_id}"
    log(f"\n[apa] === {tag}: nA={len(rowsA)} nB={len(rowsB)} ===")

    # ---- confounder audit -------------------------------------------------
    aud = {}
    for k in ["pas_counts", "pas_n_genes", "gex_counts", "gex_n_genes",
              "pct_mt", "doublet_score"]:
        a, b = M.loc[selA, k].astype(float), M.loc[selB, k].astype(float)
        aud[f"medA_{k}"] = float(a.median())
        aud[f"medB_{k}"] = float(b.median())
        aud[f"log2ratio_{k}"] = float(np.log2((a.median() + 1e-9) /
                                              (b.median() + 1e-9)))
        aud[f"mwu_p_{k}"] = float(mannwhitneyu(a, b, alternative="two-sided")[1])
    aud["fracA_dbl_cluster"] = float(M.loc[selA, "dbl_cluster"].mean())
    aud["fracB_dbl_cluster"] = float(M.loc[selB, "dbl_cluster"].mean())
    aud["depth_confounded"] = bool(abs(aud["log2ratio_pas_counts"]) >
                                   np.log2(1.5))
    log("  depth/QC audit: " + "  ".join(
        f"{k.replace('log2ratio_','L2:')}={aud[k]:+.2f}"
        for k in aud if k.startswith("log2ratio_")))
    log(f"  medians PAS counts {aud['medA_pas_counts']:.0f} vs "
        f"{aud['medB_pas_counts']:.0f} | GEX counts "
        f"{aud['medA_gex_counts']:.0f} vs {aud['medB_gex_counts']:.0f} | "
        f"doublet {aud['medA_doublet_score']:.3f} vs "
        f"{aud['medB_doublet_score']:.3f} | depth_confounded="
        f"{aud['depth_confounded']}")

    # ---- real screen + test ----------------------------------------------
    cand = screen(rowsA, rowsB)
    log(f"  candidates (|dprop|>={CAND_MIN_DPROP}): {len(cand)}")
    res = percell_test(cand, rowsA, rowsB)
    if not res.empty:
        res["gene_id"] = gene_cats[res["gene_code"].astype(int)]
        res["gene_symbol"] = res["gene_id"].map(sym).fillna("")
        res["pas_id"] = pas.var["pas_id"].values[res["pas_col"].astype(int)]
        info = pas_info.reindex(res["pas_id"].values)
        res["chrom"] = info["chrom"].values
        res["start"] = info["start"].values
        res["end"] = info["end"].values
        res["strand"] = info["strand"].values
        res["pas_class"] = info["tier"].values
        res["sig"] = res["q"] < FDR_Q
        # independent GEX total expression check
        eA, eB = [], []
        for gs in res["gene_symbol"].values:
            eA.append(gex_mean_cp10k(gA, gs))
            eB.append(gex_mean_cp10k(gB, gs))
        res["gex_cp10k_A"] = eA
        res["gex_cp10k_B"] = eB
        res["gex_log2fc"] = np.log2((np.asarray(eA) + 1) /
                                    (np.asarray(eB) + 1))
        res["expr_comparable"] = res["gex_log2fc"].abs() <= EXPR_MAX_ABS_LOG2FC
        res["apa_driven"] = res["sig"] & res["expr_comparable"]
        # PAS-side gene abundance (depth-normalised) for completeness
        res["pas_gene_log2fc"] = np.log2(
            (res["gene_countsA"] / max(M.loc[selA, "pas_counts"].sum(), 1)) /
            (res["gene_countsB"] / max(M.loc[selB, "pas_counts"].sum(), 1)))
        n_sig = int(res["sig"].sum())
        n_apa = int(res["apa_driven"].sum())
        genes_apa = sorted(set(res.loc[res["apa_driven"], "gene_symbol"]) - {""})
        log(f"  significant (q<{FDR_Q}): {n_sig} | APA-driven "
            f"(|log2FC gene GEX|<={EXPR_MAX_ABS_LOG2FC}): {n_apa}")
        log(f"  APA-driven genes ({len(genes_apa)}): "
            f"{', '.join(genes_apa[:40])}")
    else:
        n_sig = n_apa = 0
        genes_apa = []
    res.to_csv(OUT / f"apa_test_{tag}.tsv", sep="\t", index=False)

    # ---- permutation null -------------------------------------------------
    allrows = np.concatenate([rowsA, rowsB])
    nA = len(rowsA)
    null_cand, null_sig, null_apa = [], [], []
    for k in range(N_PERM):
        perm = rng.permutation(allrows)
        pA, pB = perm[:nA], perm[nA:]
        c = screen(pA, pB)
        null_cand.append(len(c))
        r = percell_test(c, pA, pB)
        ns = int((r["q"] < FDR_Q).sum()) if not r.empty else 0
        null_sig.append(ns)
        null_apa.append(ns)   # expression check is ~vacuous under permutation
    log(f"  NULL over {N_PERM} within-split relabellings: candidates "
        f"mean {np.mean(null_cand):.1f} (max {max(null_cand)}), "
        f"q<{FDR_Q} mean {np.mean(null_sig):.2f} (max {max(null_sig)})")

    results[tag] = {
        "celltype": t, "pas_A": cA_id, "pas_B": cB_id,
        "nA": int(len(rowsA)), "nB": int(len(rowsB)),
        "n_candidates": int(len(cand)), "n_significant": n_sig,
        "n_apa_driven": n_apa, "genes_apa_driven": genes_apa,
        "null_candidates_mean": float(np.mean(null_cand)),
        "null_candidates_max": int(max(null_cand)),
        "null_sig_mean": float(np.mean(null_sig)),
        "null_sig_max": int(max(null_sig)),
        "empirical_p_vs_null": float(
            (np.sum(np.asarray(null_sig) >= n_sig) + 1) / (N_PERM + 1)),
        "audit": aud,
    }

with open(OUT / "novelty_results.json", "w") as fh:
    json.dump({"concordance": conc, "laughney": LAUGHNEY,
               "splits": splits, "apa": results,
               "criteria": {
                   "SPLIT_MIN_PURITY": SPLIT_MIN_PURITY,
                   "SPLIT_MIN_FRAC": SPLIT_MIN_FRAC,
                   "SPLIT_MIN_CELLS": SPLIT_MIN_CELLS,
                   "SPLIT_MIN_TYPE_CELLS": SPLIT_MIN_TYPE_CELLS,
                   "CAND_MIN_DPROP": CAND_MIN_DPROP,
                   "CAND_MIN_GENE_COUNTS": CAND_MIN_GENE_COUNTS,
                   "CELL_MIN_GENE_COUNTS": CELL_MIN_GENE_COUNTS,
                   "CELL_MIN_N": CELL_MIN_N, "FDR_Q": FDR_Q,
                   "EXPR_MAX_ABS_LOG2FC": EXPR_MAX_ABS_LOG2FC,
                   "N_PERM": N_PERM}}, fh, indent=2)
log(f"\n[out] {OUT/'novelty_results.json'}")
log("[done]")
