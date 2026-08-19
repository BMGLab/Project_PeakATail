#!/usr/bin/env python3
"""Step 1: gene-level aggregation of the RAW PAS count matrix (06_preprocessing),
carrying the leiden labels + UMAP from 07_clustering (identical obs/var order, verified)."""
import os, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="8"
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
RUN="/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/gse104556/peakatail/{m}/run"
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
res={}
for m in ("mouse1","mouse2"):
    r=RUN.format(m=m)
    C=ad.read_h5ad(f"{r}/07_clustering/default/clusters.h5ad")
    P=ad.read_h5ad(f"{r}/06_preprocessing/default/preprocessed.h5ad")
    assert list(C.obs_names)==list(P.obs_names) and list(C.var_names)==list(P.var_names)
    P.obs["leiden"]=C.obs["leiden"].values
    P.obs["total_counts_tfidf_stage"]=C.obs["total_counts"].values
    P.obsm["X_umap"]=C.obsm["X_umap"]
    P.obs["umi_total"]=np.asarray(P.X.sum(1)).ravel().astype(int)
    ann=pd.read_csv(f"{r}/annotatedpas.bed",sep="\t",header=None,
        names=["chrom","start","end","pas_id","gene_id","symbol","strand","score","tier"],
        dtype={"chrom":str,"pas_id":str,"gene_id":str,"symbol":str,"strand":str,"tier":str})
    keep=ann.set_index("pas_id").reindex(P.var_names.astype(str))
    P.var["symbol"]=keep["symbol"].values; P.var["gene_ens"]=keep["gene_id"].values
    P.var["chrom"]=keep["chrom"].values; P.var["pstart"]=keep["start"].values
    P.var["pend"]=keep["end"].values; P.var["strand"]=keep["strand"].values
    P.var["tier"]=keep["tier"].values
    P.write_h5ad(f"{OUT}/{m}_pas_raw.h5ad")
    # gene-level sum
    sym=P.var["symbol"].astype(str).values
    ok=(sym!="nan")&(sym!="")
    idx=np.nonzero(ok)[0]
    cat=pd.Categorical(sym[idx]); genes=cat.categories
    S=sp.csr_matrix((np.ones(len(idx)),(np.arange(len(idx)),cat.codes)),shape=(len(idx),len(genes)))
    G=(P.X.tocsc()[:,idx].tocsr()@S).tocsr()
    ga=ad.AnnData(X=G,obs=P.obs.copy(),var=pd.DataFrame(index=pd.Index(genes,name="symbol")))
    ga.write_h5ad(f"{OUT}/{m}_gene.h5ad")
    print(m,"raw PAS",P.shape,"genes",G.shape,"median UMI/cell",int(np.median(P.obs['umi_total'])),flush=True)
    res[m]=dict(n_cells=int(P.n_obs),n_pas=int(P.n_vars),n_genes=int(G.shape[1]),
                n_clusters=int(P.obs["leiden"].nunique()),median_umi=float(np.median(P.obs['umi_total'])))
json.dump(res,open(f"{OUT}/step1_stats.json","w"),indent=1); print(json.dumps(res,indent=1))
