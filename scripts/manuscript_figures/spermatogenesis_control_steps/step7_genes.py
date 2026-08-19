#!/usr/bin/env python3
"""Step 7: per-gene 3'UTR shortening across the three well-populated stages."""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="8"
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from scipy.stats import spearmanr, binomtest
import importlib.util
spec=importlib.util.spec_from_file_location("s6","/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl/step6_length2.py")
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
ORDER3=["Spermatocyte","RoundSpermatid","Elongating"]
MIN_PAS_UMI=20; MIN_STAGE_UMI=50
def load(m):
    A=ad.read_h5ad(f"{OUT}/{m}_pas_repaired.h5ad"); A=A[A.obs["starsolo_cell"].values].copy()
    lab=pd.read_csv(f"{OUT}/{m}_cell_labels.tsv",sep="\t").set_index("cb"); A.obs=A.obs.join(lab)
    v=A.var.copy()
    d=pd.read_csv(f"{OUT}/{m}_pas_atlasdist.tsv",sep="\t",header=None,names=["pas_id","dd"]).drop_duplicates("pas_id").set_index("pas_id")
    v["atlas_ok"]=((d["dd"]>=0)&(d["dd"]<=100)).reindex(v["pas_id"].values).fillna(False).values
    u=pd.read_csv(f"{OUT}/{m}_pas_in_utr3.tsv",sep="\t",header=None,names=["pas_id","g"])
    v["utr3_ok"]=(v["pas_id"].astype(str)+"|"+v["symbol"].astype(str)).isin(
        set(u["pas_id"].astype(str)+"|"+u["g"].astype(str))).values
    v["site"]=np.where(v["strand"].values=="+",v["end"].values-1,v["start"].values)
    return A,v
res={}
for sub in ("utr3","atlas"):
    tabs={}
    for m in ("mouse1","mouse2"):
        A,v=load(m)
        tot=np.asarray(A.X.sum(0)).ravel(); v["tot"]=tot
        k=(v["tot"].values>=MIN_PAS_UMI)&(v["symbol"].astype(str).values!="nan")
        k&= v["utr3_ok"].values if sub=="utr3" else v["atlas_ok"].values
        v=v[k].copy(); A=A[:,v.index].copy()
        v["ordkey"]=np.where(v["strand"].values=="+",v["site"].values,-v["site"].values)
        v=v.sort_values(["symbol","ordkey"]); gg=v.groupby("symbol",sort=False)
        v["rank"]=gg.cumcount(); v["ngene"]=gg["rank"].transform("size")
        v["prox"]=gg["site"].transform("first")
        v["relpos"]=v["rank"]/np.maximum(v["ngene"]-1,1)
        v["dist_bp"]=np.abs(v["site"].values-v["prox"].values).astype(float)
        v=v[v["ngene"]>=2]; A=A[:,v.index].copy()
        symc=pd.Categorical(v["symbol"].astype(str).values)
        S=sp.csr_matrix((np.ones(len(symc)),(np.arange(len(symc)),symc.codes)),shape=(len(symc),len(symc.categories))).toarray()
        X=A.X.tocsr(); pb={}
        for s in ORDER3:
            x=np.asarray(X[(A.obs["stage_gex"].values==s)].sum(0)).ravel()
            t=x@S; nr=(x*v["relpos"].values)@S; nb=(x*v["dist_bp"].values)@S
            pb[s]=pd.DataFrame({"tot":t,"wdi":np.where(t>0,nr/np.maximum(t,1e-9),np.nan),
                                "wul":np.where(t>0,nb/np.maximum(t,1e-9),np.nan)},index=symc.categories)
        ok=np.all([pb[s]["tot"].values>=MIN_STAGE_UMI for s in ORDER3],axis=0)
        T=pd.DataFrame({f"wdi_{s}":pb[s]["wdi"].values for s in ORDER3},index=symc.categories)[ok]
        for s in ORDER3: T[f"umi_{s}"]=pb[s]["tot"].values[ok]
        T["n_pas"]=v.groupby("symbol")["rank"].size().reindex(T.index).values
        T["delta_wdi"]=T[f"wdi_{ORDER3[-1]}"]-T[f"wdi_{ORDER3[0]}"]
        T["monotone_down"]=(T[f"wdi_{ORDER3[0]}"]>T[f"wdi_{ORDER3[1]}"])&(T[f"wdi_{ORDER3[1]}"]>T[f"wdi_{ORDER3[2]}"])
        T.index.name="symbol"; tabs[m]=T
        print(f"{sub} {m}: genes with >=2 {sub} PAS and >=50 UMI in all 3 stages = {len(T)}; "
              f"frac shortening = {(T['delta_wdi']<0).mean():.3f}",flush=True)
    a,b=tabs["mouse1"],tabs["mouse2"]; com=a.index.intersection(b.index)
    both=pd.DataFrame({"m1_delta":a.loc[com,"delta_wdi"],"m2_delta":b.loc[com,"delta_wdi"],
                       "m1_mono":a.loc[com,"monotone_down"],"m2_mono":b.loc[com,"monotone_down"],
                       "m1_npas":a.loc[com,"n_pas"],"m2_npas":b.loc[com,"n_pas"]})
    for s in ORDER3:
        both[f"m1_wdi_{s}"]=a.loc[com,f"wdi_{s}"]; both[f"m2_wdi_{s}"]=b.loc[com,f"wdi_{s}"]
    both["mean_delta"]=both[["m1_delta","m2_delta"]].mean(1)
    r=spearmanr(both["m1_delta"],both["m2_delta"])
    ns=int(((both.m1_delta<0)&(both.m2_delta<0)).sum()); nl=int(((both.m1_delta>0)&(both.m2_delta>0)).sum())
    print(f"  {sub}: common={len(com)} replication rho={r.statistic:.3f} p={r.pvalue:.1e} "
          f"| both-shorten={ns} ({ns/len(com):.1%}) both-lengthen={nl} ({nl/len(com):.1%}) "
          f"binom p={binomtest(ns,ns+nl,0.5).pvalue:.2e}",flush=True)
    both.sort_values("mean_delta").to_csv(f"{OUT}/gene_shortening_{sub}.tsv",sep="\t")
    print("  top20 shortening:", ", ".join(f"{i}({v:.2f})" for i,v in both.sort_values("mean_delta")["mean_delta"].head(20).items()),flush=True)
    LIT=["Cstf2t","Ybx2","Prm1","Prm2","Prm3","Tnp1","Tnp2","Smcp","Odf1","Akap4","Ppp1cc","Crem","Camk4",
         "Hnrnpa2b1","Ldhc","Meig1","Acr","Acrv1","Gapdhs","Spata19","Oaz3","Tssk6","Pgk2","Hils1","Elavl1","Nsun7","Spem1"]
    hit=[x for x in LIT if x in both.index]
    print(f"  curated spermatogenesis-APA genes measurable: {len(hit)}/{len(LIT)}",hit,flush=True)
    if hit: print(both.loc[hit,["m1_delta","m2_delta","mean_delta","m1_mono","m2_mono"]].round(3).to_string(),flush=True)
    res[sub]=dict(n_common=len(com),rep_rho=float(r.statistic),both_short=ns,both_long=nl)
import json; json.dump(res,open(f"{OUT}/step7.json","w"),indent=1)
