#!/usr/bin/env python3
"""
pbmc_novelty_analysis.py -- STAGE 2 of the pbmc_novelty analysis.

Barcode matching, PAS-vs-GEX clustering concordance, and the NOVELTY TEST:
are there populations resolved in poly(A)-site space that gene-expression
clustering does not resolve, and if so is the split driven by APA (within-gene
site choice) rather than by depth / ambient / doublets / expression level --
and are the discriminating sites real poly(A) sites at all?

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
  REAL        ... and the discriminating site is atlas-supported: within 100 bp
              (strand-matched, 3'-most base) of a PolyASite 2.0 representative
              site. The head-to-head benchmark put PeakATail's PBMC precision
              at 0.118 @100bp, so this is the decisive filter, not a formality.
  CONFOUNDS   a split is only "novel" if it is not explained by depth, mito%,
              doublet score, ambient/off-lineage content, plain expression, or
              by GEX-Leiden already resolving it. Each is measured and reported.
  NULL        the entire pipeline is re-run on 20 random relabellings of the
              same cells within the split -> empirical false-positive count.
"""
import os
os.environ.setdefault("LC_ALL", "C")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

import json
import subprocess
import tempfile
from pathlib import Path

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
ATLAS = ROOT / "data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6"
OUT = ROOT / "results/pbmc_novelty"
GEX_H5AD = OUT / "gex_labeled.h5ad"
SCRATCH = Path("/mnt/ssd0/emaout/peakatail_benchmark/scratch")

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
ATLAS_WIN = 100
DEPTH_MAX_ABS_LOG2 = np.log2(1.5)
GEX_RESOLVED_ARI = 0.20

# ambient / soup panel: abundant, lineage-restricted transcripts. Their
# presence inside a cell of another lineage is soup, not biology.
AMBIENT_PANEL = {
    "Mono": ["LYZ", "S100A8", "S100A9", "VCAN"],
    "Platelet": ["PPBP", "PF4"],
    "Erythroid": ["HBB", "HBA1", "HBA2"],
    "B": ["IGKC", "IGHM", "CD79A"],
    "T": ["CD3D", "CD3E", "TRAC", "IL7R"],
    "NK": ["GNLY", "NKG7"],
}
TYPE_TO_LINEAGE = {"CD4-T": "T", "CD8-T": "T", "NK": "NK", "B": "B",
                   "CD14-Mono": "Mono", "FCGR3A-Mono": "Mono", "DC": "Mono",
                   "Platelet": "Platelet", "Erythroid": "Erythroid"}

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
        (f["X/data"][:].astype(np.float32), f["X/indices"][:], f["X/indptr"][:]),
        shape=tuple(f["X"].attrs["shape"]))
    raw_obs = np.array([x.decode() for x in f["obs/_index"][:]])
    raw_var = np.array([x.decode() for x in f["var/_index"][:]])
assert (raw_obs == pas.obs_names.values).all(), "raw/clusters obs order differs"
assert (raw_var == pas.var_names.values).all(), "raw/clusters var order differs"
counts_csr = raw
log(f"[load] raw counts total = {counts_csr.sum():,.0f}")

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
pas_id_arr = pas.var["pas_id"].values


# ===========================================================================
# 1b. atlas support for EVERY PAS in the matrix (point mode, strand-matched)
# ===========================================================================
ATLAS_TSV = OUT / "pas_atlas_distance.tsv"
if ATLAS_TSV.exists():
    log("[atlas] cached")
    atl = pd.read_csv(ATLAS_TSV, sep="\t")
else:
    log("[atlas] bedtools closest -s -d -t first, 3'-most base vs PolyASite 2.0")
    SCRATCH.mkdir(parents=True, exist_ok=True)
    info = pas_info.reindex(pas_id_arr)
    strand = info["strand"].values
    st = info["start"].values.astype(np.int64)
    en = info["end"].values.astype(np.int64)
    point = np.where(strand == "+", en - 1, st)
    q = pd.DataFrame({"chrom": info["chrom"].values.astype(str),
                      "start": point, "end": point + 1,
                      "name": pas_id_arr, "score": 0, "strand": strand})
    q = q.dropna(subset=["chrom"])
    with tempfile.TemporaryDirectory(dir=str(SCRATCH)) as td:
        qf = os.path.join(td, "q.bed")
        af = os.path.join(td, "a.bed")
        of = os.path.join(td, "o.bed")
        q.sort_values(["chrom", "start"]).to_csv(qf, sep="\t", header=False,
                                                 index=False)
        subprocess.run(f"sort -k1,1 -k2,2n {ATLAS} > {af}", shell=True,
                       check=True, env={**os.environ, "LC_ALL": "C"})
        subprocess.run(
            f"bedtools closest -a {qf} -b {af} -s -d -t first > {of}",
            shell=True, check=True, env={**os.environ, "LC_ALL": "C"})
        cl = pd.read_csv(of, sep="\t", header=None, usecols=[3, 12],
                         names=["pas_id", "atlas_dist"],
                         dtype={3: np.int64, 12: np.int64})
    cl = cl.drop_duplicates("pas_id")
    atl = cl
    atl.to_csv(ATLAS_TSV, sep="\t", index=False)
atlas_dist = pd.Series(atl["atlas_dist"].values,
                       index=atl["pas_id"].values).reindex(pas_id_arr).values
atlas_dist = np.where(atlas_dist < 0, np.nan, atlas_dist)   # -1 = no match
atlas_ok = atlas_dist <= ATLAS_WIN
bg_atlas_frac = float(np.nanmean(atlas_ok))
log(f"[atlas] {int(np.nansum(atlas_ok)):,}/{pas.n_vars:,} "
    f"({bg_atlas_frac:.4f}) of all PAS in the matrix are within {ATLAS_WIN} bp "
    f"of a PolyASite 2.0 representative site (strand-matched, point mode)")

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
    "pct_ribo": gex.obs["pct_counts_ribo"].values[gi],
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

# ambient / off-lineage score per matched cell (CP10K fraction)
def cp10k_sum(rows, genes):
    js = [gene_pos[g] for g in genes if g in gene_pos]
    if not js:
        return np.zeros(len(rows))
    v = np.asarray(lognorm[rows][:, js].todense())
    return np.expm1(v).sum(axis=1)


for lin, genes in AMBIENT_PANEL.items():
    M[f"amb_{lin}"] = cp10k_sum(M["gex_row"].values, genes) / 1e4
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
conc["atlas_frac_all_pas"] = bg_atlas_frac

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
purity = ct.div(ct.sum(axis=1), axis=0)
share = ct.div(ct.sum(axis=0), axis=1)

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
gene_codes = pas.var["gene_id"].cat.codes.values.astype(np.int32)
gene_cats = np.array(pas.var["gene_id"].cat.categories)
n_pas_per_gene = np.bincount(gene_codes, minlength=len(gene_cats))
multi = n_pas_per_gene[gene_codes] >= 2
log(f"\n[apa] {int(multi.sum())}/{pas.n_vars} PAS lie in multi-PAS genes "
    f"({int((n_pas_per_gene >= 2).sum())}/{len(gene_cats)} genes)")
log(f"[apa] atlas-supported fraction among multi-PAS-gene PAS: "
    f"{np.nanmean(atlas_ok[multi]):.4f}")

G = sp.csr_matrix((np.ones(pas.n_vars, np.float32),
                   (gene_codes, np.arange(pas.n_vars))),
                  shape=(len(gene_cats), pas.n_vars))


def gex_mean_cp10k(rows, gene_symbol):
    j = gene_pos.get(gene_symbol)
    if j is None:
        return np.nan
    v = np.asarray(lognorm[rows, j].todense()).ravel()
    return float(np.expm1(v).mean())


def screen(sub_csr, mA, mB):
    """Pseudobulk within-gene usage screen on a split submatrix."""
    sA = mA.astype(np.float32) @ sub_csr
    sB = mB.astype(np.float32) @ sub_csr
    sA = np.asarray(sA).ravel()
    sB = np.asarray(sB).ravel()
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


def percell_test(cand, sub_csc, cellgene_csc, mA, mB):
    """Per-cell usage test; cells are the replicate unit."""
    if cand.empty:
        return cand.assign(nA_cells=[], nB_cells=[], cell_propA=[],
                           cell_propB=[], p=[], q=[])
    ps, nAs, nBs, mAv, mBv = [], [], [], [], []
    tcache = {}
    for _, r in cand.iterrows():
        gc = int(r["gene_code"])
        if gc not in tcache:
            tcache[gc] = np.asarray(cellgene_csc[:, gc].todense()).ravel()
        tv = tcache[gc]
        j = int(r["pas_col"])
        cv = np.asarray(sub_csc[:, j].todense()).ravel()
        okA = mA & (tv >= CELL_MIN_GENE_COUNTS)
        okB = mB & (tv >= CELL_MIN_GENE_COUNTS)
        nA, nB = int(okA.sum()), int(okB.sum())
        nAs.append(nA)
        nBs.append(nB)
        if nA < CELL_MIN_N or nB < CELL_MIN_N:
            ps.append(np.nan)
            mAv.append(np.nan)
            mBv.append(np.nan)
            continue
        fA = cv[okA] / tv[okA]
        fB = cv[okB] / tv[okB]
        mAv.append(float(fA.mean()))
        mBv.append(float(fB.mean()))
        try:
            ps.append(float(mannwhitneyu(fA, fB, alternative="two-sided")[1]))
        except ValueError:
            ps.append(np.nan)
    out = cand.copy()
    out["nA_cells"] = nAs
    out["nB_cells"] = nBs
    out["cell_propA"] = mAv
    out["cell_propB"] = mBv
    out["p"] = ps
    v = ~out["p"].isna()
    out["q"] = np.nan
    if v.any():
        out.loc[v, "q"] = bh(out.loc[v, "p"].values)
    return out


def gex_de(gA_rows, gB_rows):
    """Plain-expression check: how different are the arms in the GEX matrix?"""
    A = lognorm[gA_rows]
    B = lognorm[gB_rows]
    dtA = np.asarray((A > 0).sum(axis=0)).ravel() / A.shape[0]
    dtB = np.asarray((B > 0).sum(axis=0)).ravel() / B.shape[0]
    keep = np.where((dtA >= 0.10) | (dtB >= 0.10))[0]
    if keep.size == 0:
        return 0, 0, []
    Ad = np.asarray(A[:, keep].todense())
    Bd = np.asarray(B[:, keep].todense())
    mA = np.expm1(Ad).mean(axis=0)
    mB = np.expm1(Bd).mean(axis=0)
    l2 = np.log2((mA + 1) / (mB + 1))
    with np.errstate(invalid="ignore"):
        st, p = mannwhitneyu(Ad, Bd, alternative="two-sided", axis=0)
    p = np.nan_to_num(p, nan=1.0)
    q = bh(p)
    sig = (q < 0.05)
    strong = sig & (np.abs(l2) >= 1.0)
    order = np.argsort(-np.abs(l2) * strong)
    top = [(str(lognorm_genes[keep[i]]), float(l2[i]))
           for i in order[:15] if strong[i]]
    return int(sig.sum()), int(strong.sum()), top


results = {}
rng = np.random.default_rng(SEED)
for s in splits:
    t = s["celltype"]
    cA_id, cB_id = s["pas_clusters"][0], s["pas_clusters"][1]
    selA = ((M["celltype"] == t) & (M["pas_leiden"] == cA_id)).values
    selB = ((M["celltype"] == t) & (M["pas_leiden"] == cB_id)).values
    rowsA = M.loc[selA, "pas_row"].values
    rowsB = M.loc[selB, "pas_row"].values
    gA = M.loc[selA, "gex_row"].values
    gB = M.loc[selB, "gex_row"].values
    tag = f"{t}_PAS{cA_id}_vs_PAS{cB_id}"
    log(f"\n[apa] === {tag}: nA={len(rowsA)} nB={len(rowsB)} ===")

    # ---- split submatrix (once) ------------------------------------------
    allrows = np.concatenate([rowsA, rowsB])
    nA = len(rowsA)
    sub_csr = counts_csr[allrows].tocsr()
    sub_csc = sub_csr.tocsc()
    cellgene_csc = (sub_csr @ G.T).tocsc()
    mA0 = np.zeros(len(allrows), bool)
    mA0[:nA] = True
    mB0 = ~mA0

    # ---- confounder audit -------------------------------------------------
    aud = {}
    audit_keys = ["pas_counts", "pas_n_genes", "gex_counts", "gex_n_genes",
                  "pct_mt", "pct_ribo", "doublet_score"] + \
                 [c for c in M.columns if c.startswith("amb_")]
    for k in audit_keys:
        a = M.loc[selA, k].astype(float)
        b = M.loc[selB, k].astype(float)
        aud[f"medA_{k}"] = float(a.median())
        aud[f"medB_{k}"] = float(b.median())
        aud[f"log2ratio_{k}"] = float(np.log2((a.median() + 1e-9) /
                                              (b.median() + 1e-9)))
        aud[f"mwu_p_{k}"] = float(mannwhitneyu(a, b, alternative="two-sided")[1])
    aud["fracA_dbl_cluster"] = float(M.loc[selA, "dbl_cluster"].mean())
    aud["fracB_dbl_cluster"] = float(M.loc[selB, "dbl_cluster"].mean())
    aud["depth_confounded"] = bool(abs(aud["log2ratio_pas_counts"]) >
                                   DEPTH_MAX_ABS_LOG2)
    aud["gexdepth_confounded"] = bool(abs(aud["log2ratio_gex_counts"]) >
                                      DEPTH_MAX_ABS_LOG2)
    aud["mito_confounded"] = bool(abs(aud["log2ratio_pct_mt"]) >
                                  DEPTH_MAX_ABS_LOG2)
    aud["doublet_confounded"] = bool(abs(aud["log2ratio_doublet_score"]) >
                                     DEPTH_MAX_ABS_LOG2)
    own = TYPE_TO_LINEAGE.get(t, t)
    amb_off = [f"amb_{l}" for l in AMBIENT_PANEL if l != own]
    aud["ambient_max_abs_log2"] = float(max(abs(aud[f"log2ratio_{k}"])
                                            for k in amb_off))
    aud["ambient_worst_panel"] = max(amb_off,
                                     key=lambda k: abs(aud[f"log2ratio_{k}"]))
    aud["ambient_confounded"] = bool(aud["ambient_max_abs_log2"] >
                                     DEPTH_MAX_ABS_LOG2)

    # GEX already resolves the split?
    arm = np.where(selA[selA | selB], "A", "B")
    sub_gexleiden = M.loc[selA | selB, "gex_leiden"].values
    ari_arm_gexleiden = float(adjusted_rand_score(arm, sub_gexleiden))
    ami_arm_gexleiden = float(adjusted_mutual_info_score(arm, sub_gexleiden))
    modeA = pd.Series(M.loc[selA, "gex_leiden"]).value_counts()
    modeB = pd.Series(M.loc[selB, "gex_leiden"]).value_counts()
    aud["ARI_arm_vs_gexleiden"] = ari_arm_gexleiden
    aud["AMI_arm_vs_gexleiden"] = ami_arm_gexleiden
    aud["modal_gexleiden_A"] = str(modeA.index[0])
    aud["modal_gexleiden_B"] = str(modeB.index[0])
    aud["modal_fracA"] = float(modeA.iloc[0] / modeA.sum())
    aud["modal_fracB"] = float(modeB.iloc[0] / modeB.sum())
    aud["gex_already_resolves"] = bool(ari_arm_gexleiden >= GEX_RESOLVED_ARI)
    n_de, n_de_strong, top_de = gex_de(gA, gB)
    aud["n_gex_de_q05"] = n_de
    aud["n_gex_de_q05_lfc1"] = n_de_strong
    aud["top_gex_de"] = top_de
    aud["expression_confounded"] = bool(n_de_strong >= 10)

    log("  audit  " + "  ".join(
        f"{k[len('log2ratio_'):]}={aud[k]:+.2f}"
        for k in aud if k.startswith("log2ratio_")))
    log(f"  depth_conf={aud['depth_confounded']} "
        f"mito_conf={aud['mito_confounded']} "
        f"dbl_conf={aud['doublet_confounded']} "
        f"ambient_conf={aud['ambient_confounded']} "
        f"({aud['ambient_worst_panel']} {aud['ambient_max_abs_log2']:+.2f}) "
        f"expr_conf={aud['expression_confounded']} "
        f"(n_DE q<0.05 |lfc|>1 = {n_de_strong}) "
        f"gex_resolves={aud['gex_already_resolves']} "
        f"(ARI arm~gexleiden {ari_arm_gexleiden:.3f})")
    if top_de:
        log("  top GEX DE: " + ", ".join(f"{g}({l:+.1f})" for g, l in top_de[:10]))

    # ---- real screen + test ----------------------------------------------
    cand = screen(sub_csr, mA0, mB0)
    log(f"  candidates (|dprop|>={CAND_MIN_DPROP}): {len(cand)}")
    res = percell_test(cand, sub_csc, cellgene_csc, mA0, mB0)
    if not res.empty:
        res["gene_id"] = gene_cats[res["gene_code"].astype(int)]
        res["gene_symbol"] = res["gene_id"].map(sym).fillna("")
        res["pas_id"] = pas_id_arr[res["pas_col"].astype(int)]
        info = pas_info.reindex(res["pas_id"].values)
        res["chrom"] = info["chrom"].values
        res["start"] = info["start"].values
        res["end"] = info["end"].values
        res["strand"] = info["strand"].values
        res["pas_class"] = info["tier"].values
        res["atlas_dist"] = atlas_dist[res["pas_col"].astype(int)]
        res["atlas_supported"] = res["atlas_dist"] <= ATLAS_WIN
        res["sig"] = res["q"] < FDR_Q
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
        res["apa_driven_atlas"] = res["apa_driven"] & res["atlas_supported"]
        res["pas_gene_log2fc"] = np.log2(
            (res["gene_countsA"] / max(M.loc[selA, "pas_counts"].sum(), 1)) /
            (res["gene_countsB"] / max(M.loc[selB, "pas_counts"].sum(), 1)))
        n_sig = int(res["sig"].sum())
        n_apa = int(res["apa_driven"].sum())
        n_apa_atlas = int(res["apa_driven_atlas"].sum())
        genes_apa = sorted(set(res.loc[res["apa_driven"], "gene_symbol"]) - {""})
        genes_apa_atlas = sorted(
            set(res.loc[res["apa_driven_atlas"], "gene_symbol"]) - {""})
        af_cand = float(res["atlas_supported"].mean()) if len(res) else np.nan
        af_sig = float(res.loc[res["sig"], "atlas_supported"].mean()) \
            if n_sig else np.nan
        af_apa = float(res.loc[res["apa_driven"], "atlas_supported"].mean()) \
            if n_apa else np.nan
        log(f"  significant (q<{FDR_Q}): {n_sig} | APA-driven "
            f"(|log2FC gene GEX|<={EXPR_MAX_ABS_LOG2FC}): {n_apa} | "
            f"AND atlas-supported (<={ATLAS_WIN}bp): {n_apa_atlas}")
        log(f"  atlas-supported fraction: candidates {af_cand:.3f} | "
            f"significant {af_sig:.3f} | APA-driven {af_apa:.3f} | "
            f"background(all PAS) {bg_atlas_frac:.3f}")
        log(f"  APA-driven genes ({len(genes_apa)}): "
            f"{', '.join(genes_apa[:40])}")
        log(f"  APA-driven AND atlas-supported genes ({len(genes_apa_atlas)}): "
            f"{', '.join(genes_apa_atlas[:40])}")
    else:
        n_sig = n_apa = n_apa_atlas = 0
        genes_apa = genes_apa_atlas = []
        af_cand = af_sig = af_apa = float("nan")
    res.to_csv(OUT / f"apa_test_{tag}.tsv", sep="\t", index=False)

    # ---- permutation null -------------------------------------------------
    null_cand, null_sig = [], []
    for k in range(N_PERM):
        perm = rng.permutation(len(allrows))
        mA = np.zeros(len(allrows), bool)
        mA[perm[:nA]] = True
        mB = ~mA
        c = screen(sub_csr, mA, mB)
        null_cand.append(len(c))
        r = percell_test(c, sub_csc, cellgene_csc, mA, mB)
        ns = int((r["q"] < FDR_Q).sum()) if not r.empty else 0
        null_sig.append(ns)
    log(f"  NULL over {N_PERM} within-split relabellings: candidates "
        f"mean {np.mean(null_cand):.1f} (max {max(null_cand)}), "
        f"q<{FDR_Q} mean {np.mean(null_sig):.2f} (max {max(null_sig)})")

    results[tag] = {
        "celltype": t, "pas_A": cA_id, "pas_B": cB_id,
        "nA": int(len(rowsA)), "nB": int(len(rowsB)),
        "n_candidates": int(len(cand)), "n_significant": n_sig,
        "n_apa_driven": n_apa, "n_apa_driven_atlas": n_apa_atlas,
        "genes_apa_driven": genes_apa,
        "genes_apa_driven_atlas": genes_apa_atlas,
        "atlas_frac_candidates": af_cand,
        "atlas_frac_significant": af_sig,
        "atlas_frac_apa_driven": af_apa,
        "atlas_frac_background": bg_atlas_frac,
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
                   "ATLAS_WIN": ATLAS_WIN,
                   "DEPTH_MAX_ABS_LOG2": float(DEPTH_MAX_ABS_LOG2),
                   "GEX_RESOLVED_ARI": GEX_RESOLVED_ARI,
                   "N_PERM": N_PERM}}, fh, indent=2, default=str)
log(f"\n[out] {OUT/'novelty_results.json'}")
log("[done]")
