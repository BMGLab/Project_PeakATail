#!/usr/bin/env python3
"""Step 4: cluster + stage-annotate, independently (a) from the REPAIRED PeakATail
PAS matrix aggregated to gene level and (b) from the orthogonal STARsolo GEX matrix."""
import os, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="8"
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
from scipy.io import mmread
sc.settings.n_jobs=8; sc.settings.verbosity=0
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
SS="/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/gse104556/starsolo/{S}/Solo.out/Gene/filtered"

MARKERS={
 "Spermatogonia":  ["Dazl","Kit","Uchl1","Zbtb16","Gfra1","Stra8","Crabp1","Dmrt1","Egr4","Sall4","Nanos3","Utf1","Lin28a"],
 "Spermatocyte":   ["Sycp3","Sycp1","Sycp2","Piwil1","Spo11","Hormad1","Meiob","Tex101","Ccnb3","Rad51ap2","Insl6"],
 "RoundSpermatid": ["Acr","Acrv1","Txndc8","Tssk1","Spata19","Catsper3","Izumo1","Sun5","Lyzl1","Tex21","Prss21","Meig1","Tsga8","Spaca1"],
 "Elongating":     ["Prm1","Prm2","Prm3","Tnp1","Tnp2","Oaz3","Smcp","Hils1","Actl7b","Tssk6","Gapdhs","Odf1"],
 "Sertoli":        ["Sox9","Amh","Wt1","Rhox8","Cldn11","Shbg"],
 "Leydig":         ["Cyp17a1","Cyp11a1","Hsd3b1","Insl3","Star","Hsd17b3"],
 "Immune":         ["Ptprc","Cd52","Laptm5","C1qa","C1qb","Adgre1","Tyrobp"],
}
GERM=["Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"]
SOMATIC_ABS_MIN=0.35

SOMATIC=["Sertoli","Leydig","Immune"]
def score_clusters(A, key):
    scols=[]
    for k,gl in MARKERS.items():
        gg=[g for g in gl if g in A.var_names]
        if gg:
            sc.tl.score_genes(A,gg,score_name=f"s_{k}",ctrl_size=50,random_state=0)
        else:
            A.obs[f"s_{k}"]=np.nan
        scols.append(f"s_{k}")
    cl=A.obs[key].astype(str)
    per=A.obs.groupby(cl,observed=True)[scols].mean()
    per["n_cells"]=cl.value_counts().reindex(per.index).values
    # absolute mean marker expression per panel (log1p CP10K) -- used for the
    # somatic guard, because a pure z-argmax forces every cluster into a panel
    # even when that panel's genes are not expressed at all.
    absm={}
    for k,gl in MARKERS.items():
        gg=[g for g in gl if g in A.var_names]
        absm[f"abs_{k}"]=(np.asarray(A[:,gg].X.todense()).mean(1) if gg else np.zeros(A.n_obs))
    absdf=pd.DataFrame(absm,index=A.obs_names).groupby(cl.values).mean()
    per=per.join(absdf)
    z=per[scols].apply(lambda c:(c-c.mean())/(c.std(ddof=0)+1e-12))
    zg=z[[f"s_{k}" for k in GERM]]
    stage=zg.idxmax(axis=1).str.replace("s_","",regex=False)
    # somatic override: a somatic panel must be expressed in ABSOLUTE terms
    som_abs=per[[f"abs_{k}" for k in SOMATIC]]
    som_hit=som_abs.max(axis=1)>SOMATIC_ABS_MIN
    stage=stage.where(~som_hit, som_abs.idxmax(axis=1).str.replace("abs_","",regex=False))
    per["stage"]=stage
    per["somatic_abs_max"]=som_abs.max(axis=1)
    per["margin"]=zg.max(axis=1)-zg.apply(lambda r:r.nlargest(2).iloc[-1],axis=1)
    return per, z.add_prefix("z_")

def tfidf_lsi(A,n=50,seed=42):
    X=A.X.tocsr().astype(np.float64)
    cs=np.asarray(X.sum(1)).ravel(); cs[cs==0]=1
    tf=sp.diags(1e4/cs)@X
    nz=np.asarray((X>0).sum(0)).ravel(); idf=np.log(1+X.shape[0]/np.maximum(nz,1))
    T=(tf@sp.diags(idf)).log1p()
    from sklearn.utils.extmath import randomized_svd
    U,S,_=randomized_svd(T,n_components=n,random_state=seed)
    L=U*S
    return L[:,1:]   # drop comp 0 (depth)

res={}
for m,S in [("mouse1","Mouse1_scRNAseq"),("mouse2","Mouse2_scRNAseq")]:
    # ---------- (a) PAS-derived ----------
    P=ad.read_h5ad(f"{OUT}/{m}_pas_repaired.h5ad")
    P=P[P.obs["starsolo_cell"].values].copy()
    sym=P.var["symbol"].astype(str).values; ok=(sym!="nan")&(sym!="")
    idx=np.nonzero(ok)[0]; cat=pd.Categorical(sym[idx])
    Smat=sp.csr_matrix((np.ones(len(idx)),(np.arange(len(idx)),cat.codes)),shape=(len(idx),len(cat.categories)))
    Gm=(P.X.tocsc()[:,idx].tocsr()@Smat).tocsr()
    Gp=ad.AnnData(X=Gm,obs=P.obs.copy(),var=pd.DataFrame(index=pd.Index(cat.categories,name="symbol")))
    Gp.obsm["X_lsi"]=tfidf_lsi(Gp)
    sc.pp.neighbors(Gp,use_rep="X_lsi",n_neighbors=15,random_state=42)
    sc.tl.leiden(Gp,resolution=1.0,key_added="clu",random_state=42,flavor="igraph",n_iterations=2,directed=False)
    sc.tl.umap(Gp,random_state=42)
    Gp.layers["counts"]=Gp.X.copy(); sc.pp.normalize_total(Gp,target_sum=1e4); sc.pp.log1p(Gp)
    perP,zP=score_clusters(Gp,"clu")

    # ---------- (b) STARsolo GEX ----------
    feat=pd.read_csv(f"{SS.format(S=S)}/features.tsv",sep="\t",header=None,names=["ens","sym","t"])
    bc=pd.read_csv(f"{SS.format(S=S)}/barcodes.tsv",header=None)[0].astype(str).values
    M=mmread(f"{SS.format(S=S)}/matrix.mtx").tocsr().T.tocsr()   # cells x genes
    Gg=ad.AnnData(X=M,obs=pd.DataFrame(index=pd.Index(bc,name="cb")),
                  var=pd.DataFrame(index=pd.Index(feat["sym"].values,name="symbol")))
    Gg.var_names_make_unique()
    Gg=Gg[Gp.obs_names].copy()
    sc.pp.filter_genes(Gg,min_cells=3)
    Gg.layers["counts"]=Gg.X.copy()
    sc.pp.normalize_total(Gg,target_sum=1e4); sc.pp.log1p(Gg)
    sc.pp.highly_variable_genes(Gg,n_top_genes=2000)
    Gg.raw=Gg
    sc.pp.pca(Gg,n_comps=30,mask_var="highly_variable",random_state=42)
    sc.pp.neighbors(Gg,n_neighbors=15,random_state=42)
    sc.tl.leiden(Gg,resolution=1.0,key_added="clu",random_state=42,flavor="igraph",n_iterations=2,directed=False)
    sc.tl.umap(Gg,random_state=42)
    perG,zG=score_clusters(Gg,"clu")

    for tag,per,z,Aobj in [("PAS",perP,zP,Gp),("GEX",perG,zG,Gg)]:
        out=pd.concat([per,z],axis=1); out.index.name="cluster"; out["mouse"]=m; out["source"]=tag
        out.to_csv(f"{OUT}/{m}_{tag}_cluster_stage.tsv",sep="\t")
        Aobj.obs["stage"]=Aobj.obs["clu"].astype(str).map(per["stage"]).values
        print(f"--- {m} {tag}")
        print(per[["n_cells","stage","margin"]].sort_index(key=lambda i:i.astype(int)).to_string(),flush=True)
        print("   stage cell counts:",Aobj.obs["stage"].value_counts().to_dict(),flush=True)
    lab=pd.DataFrame({"cb":Gp.obs_names,"clu_pas":Gp.obs["clu"].astype(str).values,
                      "stage_pas":Gp.obs["stage"].values,
                      "clu_gex":Gg.obs["clu"].astype(str).values,"stage_gex":Gg.obs["stage"].values,
                      "umi_pas":Gp.obs["umi_total"].values})
    lab["umap1_pas"]=Gp.obsm["X_umap"][:,0]; lab["umap2_pas"]=Gp.obsm["X_umap"][:,1]
    lab["umap1_gex"]=Gg.obsm["X_umap"][:,0]; lab["umap2_gex"]=Gg.obsm["X_umap"][:,1]
    lab.to_csv(f"{OUT}/{m}_cell_labels.tsv",sep="\t",index=False)
    ct=pd.crosstab(lab["stage_pas"],lab["stage_gex"])
    print(f"=== {m} stage concordance PAS vs GEX\n{ct.to_string()}",flush=True)
    agree=(lab["stage_pas"]==lab["stage_gex"]).mean()
    res[m]=dict(agreement=float(agree),n=int(len(lab)))
    print(f"=== {m} raw agreement {agree:.3f}",flush=True)
    # marker expression matrices for the audit panel
    for tag,Aobj in [("PAS",Gp),("GEX",Gg)]:
        mk=[g for gl in MARKERS.values() for g in gl if g in Aobj.var_names]
        E=pd.DataFrame(np.asarray(Aobj[:,mk].X.todense()),columns=mk,index=Aobj.obs_names)
        E["cluster"]=Aobj.obs["clu"].astype(str).values
        E.groupby("cluster").mean().to_csv(f"{OUT}/{m}_{tag}_marker_expr.tsv",sep="\t")
json.dump(res,open(f"{OUT}/step4_stats.json","w"),indent=1)
