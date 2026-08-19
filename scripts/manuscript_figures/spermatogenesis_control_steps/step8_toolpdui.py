#!/usr/bin/env python3
"""Step 8: repeat the trend with the TOOL'S OWN strand-aware PDUI implementation
(ema.quantification.pdui.calculate_pdui / calculate_pdui_per_cluster), on stage
pseudobulk built from the repaired count matrix."""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="8"
sys.path.insert(0,"/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail")
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import spearmanr, wilcoxon
from ema.quantification.pdui import calculate_pdui
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
ORD=["Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"]
rows=[]; per=[]
for m in ("mouse1","mouse2"):
    A=ad.read_h5ad(f"{OUT}/{m}_pas_repaired.h5ad"); A=A[A.obs["starsolo_cell"].values].copy()
    lab=pd.read_csv(f"{OUT}/{m}_cell_labels.tsv",sep="\t").set_index("cb"); A.obs=A.obs.join(lab)
    v=A.var.copy()
    d=pd.read_csv(f"{OUT}/{m}_pas_atlasdist.tsv",sep="\t",header=None,names=["pas_id","dd"]).drop_duplicates("pas_id").set_index("pas_id")
    v["atlas_ok"]=((d["dd"]>=0)&(d["dd"]<=100)).reindex(v["pas_id"].values).fillna(False).values
    u=pd.read_csv(f"{OUT}/{m}_pas_in_utr3.tsv",sep="\t",header=None,names=["pas_id","g"])
    v["utr3_ok"]=(v["pas_id"].astype(str)+"|"+v["symbol"].astype(str)).isin(
        set(u["pas_id"].astype(str)+"|"+u["g"].astype(str))).values
    for sub in ("all","atlas","utr3","noatlas"):
        k=(v["symbol"].astype(str).values!="nan")
        if sub=="atlas": k&=v["atlas_ok"].values
        elif sub=="utr3": k&=v["utr3_ok"].values
        elif sub=="noatlas": k&=~v["atlas_ok"].values
        vv=v[k]; AA=A[:,vv.index]
        # stage pseudobulk count matrix, PAS x stage  (the tool wants PAS rows)
        cnt={}
        for s in ORD:
            sel=(AA.obs["stage_gex"].values==s)
            if sel.sum()<30: continue
            cnt[s]=np.asarray(AA.X[sel].sum(0)).ravel()
        stages=list(cnt); cm=pd.DataFrame(cnt,index=vv["pas_id"].astype(int).values)
        pas_info=pd.DataFrame({"chrom":vv["chrom"].astype(str).values,"start":vv["start"].values,
                               "end":vv["end"].values,"pas_id":vv["pas_id"].astype(int).values,
                               "score":0,"strand":vv["strand"].astype(str).values,
                               "gene_id":vv["symbol"].astype(str).values})
        pd_df=calculate_pdui(cm,pas_info,pseudocount=0.0)
        # require >=50 UMI at the prox+distal pair in every stage
        keep=pd_df.notna().all(axis=1)
        pd_df=pd_df[keep][stages]
        # depth guard: gene must have >=50 UMI in each stage
        tot=cm.groupby(pas_info.set_index("pas_id").loc[cm.index,"gene_id"].values).sum()
        ok=(tot.reindex(pd_df.index)[stages]>=50).all(axis=1)
        P=pd_df[ok.values]
        means=[float(P[s].mean()) for s in stages]
        rho,p=spearmanr(np.repeat(np.arange(len(stages)),len(P)),
                        np.concatenate([P[s].values for s in stages]))
        w=wilcoxon(P[stages[-1]],P[stages[0]])
        rows.append(dict(mouse=m,pas_subset=sub,n_genes=int(len(P)),stages=",".join(stages),
                         **{f"PDUI_{s}":mm for s,mm in zip(stages,means)},
                         monotone_down=bool(all(means[i]>means[i+1] for i in range(len(means)-1))),
                         delta_last_first=means[-1]-means[0],
                         gene_rho=float(rho),gene_p=float(p),
                         wilcoxon_p=float(w.pvalue),
                         frac_genes_shortening=float((P[stages[-1]]<P[stages[0]]).mean())))
        pp=P.copy(); pp["mouse"]=m; pp["pas_subset"]=sub; pp.index.name="symbol"; per.append(pp.reset_index())
R=pd.DataFrame(rows); R.to_csv(f"{OUT}/tool_pdui_trend.tsv",sep="\t",index=False)
pd.concat(per).to_csv(f"{OUT}/tool_pdui_per_gene.tsv.gz",sep="\t",index=False)
pd.set_option("display.width",250)
print(R.round(4).to_string())
