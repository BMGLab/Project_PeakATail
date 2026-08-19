#!/usr/bin/env python3
"""Step 6: final 3'UTR length trend test.
Subsets: all / atlas-supported / atlas-unsupported / inside-own-3'UTR.
Depth-matched (binomial thinning) sensitivity + fixed-gene-set pseudobulk."""
import os, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="8"
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from scipy.stats import spearmanr, mannwhitneyu
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
ORDER=["Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"]
MIN_PAS_UMI=50; MIN_GENE_UMI_CELL=5; MIN_GENES_CELL=20; MIN_STAGE_UMI=100
TARGET_UMI=6000
rng=np.random.default_rng(42)

def partial_spearman(x,y,covs):
    """Spearman of x vs y after linearly removing rank-transformed covariates."""
    from scipy.stats import rankdata, t as tdist
    import numpy as np
    X=rankdata(x).astype(float); Y=rankdata(y).astype(float)
    C=np.column_stack([np.ones(len(X))]+[rankdata(c).astype(float) for c in covs])
    bx=np.linalg.lstsq(C,X,rcond=None)[0]; by=np.linalg.lstsq(C,Y,rcond=None)[0]
    rx=X-C@bx; ry=Y-C@by
    r=float(np.corrcoef(rx,ry)[0,1])
    n=len(X); k=C.shape[1]-1; dof=n-2-k
    tv=r*np.sqrt(dof/max(1e-12,1-r**2)); p=float(2*tdist.sf(abs(tv),dof))
    return r,p

def cliffs(a,b):
    n=len(a)*len(b)
    if n==0: return np.nan
    return 2*mannwhitneyu(a,b,alternative="two-sided").statistic/n-1

def thin(X,target,rs):
    X=X.tocsr().astype(np.int64).copy()
    tot=np.asarray(X.sum(1)).ravel()
    keep=tot>=target
    X=X[keep]; tot=tot[keep]
    p=np.minimum(target/np.maximum(tot,1),1.0)
    pr=np.repeat(p,np.diff(X.indptr))
    X.data=rs.binomial(X.data,pr)
    X.eliminate_zeros()
    return X,keep

def load(m):
    A=ad.read_h5ad(f"{OUT}/{m}_pas_repaired.h5ad")
    A=A[A.obs["starsolo_cell"].values].copy()
    lab=pd.read_csv(f"{OUT}/{m}_cell_labels.tsv",sep="\t").set_index("cb")
    A.obs=A.obs.join(lab)
    v=A.var.copy()
    dist=pd.read_csv(f"{OUT}/{m}_pas_atlasdist.tsv",sep="\t",header=None,
                     names=["pas_id","atlas_dist"]).drop_duplicates("pas_id").set_index("pas_id")
    v["atlas_dist"]=dist["atlas_dist"].reindex(v["pas_id"].values).values
    v["atlas_ok"]=(v["atlas_dist"]>=0)&(v["atlas_dist"]<=100)
    u=pd.read_csv(f"{OUT}/{m}_pas_in_utr3.tsv",sep="\t",header=None,names=["pas_id","utr_gene"])
    u["k"]=u["pas_id"].astype(str)+"|"+u["utr_gene"].astype(str)
    key=v["pas_id"].astype(str)+"|"+v["symbol"].astype(str)
    v["utr3_ok"]=key.isin(set(u["k"])).values
    v["site"]=np.where(v["strand"].values=="+", v["end"].values-1, v["start"].values)
    return A,v

def prep(A,v,subset):
    tot=np.asarray(A.X.sum(0)).ravel(); v=v.copy(); v["tot"]=tot
    keep=(v["tot"].values>=MIN_PAS_UMI)&(v["symbol"].astype(str).values!="nan")
    if subset=="atlas": keep&=v["atlas_ok"].values
    elif subset=="noatlas": keep&=~v["atlas_ok"].values
    elif subset=="utr3": keep&=v["utr3_ok"].values
    v=v[keep].copy()
    v["ordkey"]=np.where(v["strand"].values=="+", v["site"].values, -v["site"].values)
    v=v.sort_values(["symbol","ordkey"])
    g=v.groupby("symbol",sort=False)
    v["rank"]=g.cumcount(); v["ngene"]=g["rank"].transform("size")
    v["prox_site"]=g["site"].transform("first")
    v["relpos"]=v["rank"]/np.maximum(v["ngene"]-1,1)
    v["dist_bp"]=np.abs(v["site"].values-v["prox_site"].values).astype(float)
    v=v[v["ngene"]>=2]
    return A[:,v.index].copy(), v

def gene_ops(v):
    symc=pd.Categorical(v["symbol"].astype(str).values)
    S=sp.csr_matrix((np.ones(len(symc)),(np.arange(len(symc)),symc.codes)),
                    shape=(len(symc),len(symc.categories)))
    return S,symc.categories

def percell(X,v,S):
    Xc=X.tocsc().astype(float)
    tot=(Xc@S).toarray()
    nr=(Xc.multiply(sp.csr_matrix(v["relpos"].values.reshape(1,-1)))@S).toarray()
    nb=(Xc.multiply(sp.csr_matrix(v["dist_bp"].values.reshape(1,-1)))@S).toarray()
    use=tot>=MIN_GENE_UMI_CELL
    with np.errstate(invalid="ignore",divide="ignore"):
        wdi=np.where(use,nr/np.maximum(tot,1e-9),np.nan); wul=np.where(use,nb/np.maximum(tot,1e-9),np.nan)
    n=use.sum(1)
    return (np.where(n>=MIN_GENES_CELL,np.nanmean(wdi,axis=1),np.nan),
            np.where(n>=MIN_GENES_CELL,np.nanmean(wul,axis=1),np.nan), n)

rows=[]; cells=[]; genes=[]; pbrows=[]
for m in ("mouse1","mouse2"):
    A0,v0=load(m)
    for subset in ("all","atlas","noatlas","utr3"):
        A,v=prep(A0,v0,subset); S,gnames=gene_ops(v)
        for depth in ("raw","matched"):
            if depth=="raw":
                X=A.X; obs=A.obs
            else:
                X,kp=thin(A.X,TARGET_UMI,np.random.default_rng(42)); obs=A.obs[kp]
            wdi,wul,nu=percell(X,v,S)
            for lc,tag in [("stage_pas","PASlabels"),("stage_gex","GEXlabels")]:
                d=pd.DataFrame({"cb":obs.index,"stage":obs[lc].values,"wdi":wdi,"wul":wul,
                                "n_genes_used":nu,"umi":np.asarray(X.sum(1)).ravel()})
                stages=[s for s in ORDER if (d["stage"]==s).sum()>=30]
                d=d[d["stage"].isin(stages)].dropna(subset=["wdi"]).copy()
                d["rank"]=d["stage"].map({s:i for i,s in enumerate(stages)})
                d["mouse"]=m; d["pas_subset"]=subset; d["depth"]=depth; d["labels"]=tag
                cells.append(d)
                r={"mouse":m,"pas_subset":subset,"depth":depth,"labels":tag,
                   "n_cells":len(d),"n_pas":int(A.n_vars),"n_genes":len(gnames),
                   "stages":",".join(stages)}
                for metric in ("wdi","wul"):
                    rho,p=spearmanr(d["rank"],d[metric]); r[f"{metric}_rho"]=rho; r[f"{metric}_p"]=p
                    mn=d.groupby("stage")[metric].mean().reindex(stages)
                    sd=d.groupby("stage")[metric].std().reindex(stages)
                    for s in stages: r[f"{metric}_mean_{s}"]=mn[s]; r[f"{metric}_sd_{s}"]=sd[s]
                    a=d[d["stage"]==stages[0]][metric].values; b=d[d["stage"]==stages[-1]][metric].values
                    r[f"{metric}_delta_last_minus_first"]=float(np.mean(b)-np.mean(a))
                    r[f"{metric}_cliffs"]=float(cliffs(b,a))
                    ps=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
                    r[f"{metric}_cohens_d"]=float((b.mean()-a.mean())/ps)
                    r[f"{metric}_monotone_down"]=bool(all(mn.values[i]>mn.values[i+1] for i in range(len(stages)-1)))
                # label-shuffle null
                pr,pp=partial_spearman(d["rank"].values,d["wdi"].values,
                                       [np.log10(d["umi"].values+1),d["n_genes_used"].values])
                r["wdi_partial_rho_depthadj"]=pr; r["wdi_partial_p_depthadj"]=pp
                nul=[spearmanr(rng.permutation(d["rank"].values),d["wdi"].values).statistic for _ in range(200)]
                r["wdi_null_rho_absmax"]=float(np.max(np.abs(nul)))
                rows.append(r)
        # --- fixed-gene-set pseudobulk on RAW counts, GEX labels ---
        stages=[s for s in ORDER if (A.obs["stage_gex"]==s).sum()>=30]
        Xr=A.X.tocsr(); pb={}
        for s in stages:
            x=np.asarray(Xr[(A.obs["stage_gex"].values==s)].sum(0)).ravel()
            tot=x@S.toarray(); nr=(x*v["relpos"].values)@S.toarray(); nb=(x*v["dist_bp"].values)@S.toarray()
            pb[s]=pd.DataFrame({"tot":tot,"wdi":np.where(tot>0,nr/np.maximum(tot,1e-9),np.nan),
                                "wul":np.where(tot>0,nb/np.maximum(tot,1e-9),np.nan)},index=gnames)
        ok=np.all([pb[s]["tot"].values>=MIN_STAGE_UMI for s in stages],axis=0)
        gg=pd.DataFrame({f"wdi_{s}":pb[s]["wdi"].values for s in stages},index=gnames)[ok]
        for s in stages: gg[f"tot_{s}"]=pb[s]["tot"].values[ok]
        gg["delta_wdi"]=gg[f"wdi_{stages[-1]}"]-gg[f"wdi_{stages[0]}"]
        gg["mouse"]=m; gg["pas_subset"]=subset; gg.index.name="symbol"
        genes.append(gg.reset_index())
        # fixed-gene-set stage means (removes gene-composition drift)
        rr={"mouse":m,"pas_subset":subset,"n_genes_fixed":int(ok.sum()),"stages":",".join(stages)}
        for s in stages: rr[f"pb_wdi_{s}"]=float(np.nanmean(gg[f"wdi_{s}"]))
        vals=[rr[f"pb_wdi_{s}"] for s in stages]
        rr["pb_monotone_down"]=bool(all(vals[i]>vals[i+1] for i in range(len(vals)-1)))
        rr["pb_delta_last_minus_first"]=vals[-1]-vals[0]
        rr["pb_frac_genes_shortening"]=float((gg["delta_wdi"]<0).mean())
        rho,p=spearmanr(np.repeat(np.arange(len(stages)),ok.sum()),
                        np.concatenate([gg[f"wdi_{s}"].values for s in stages]))
        rr["pb_rho"]=float(rho); rr["pb_p"]=float(p)
        pbrows.append(rr)
        print(f"{m} {subset}: fixedgenes={int(ok.sum())} pb_means={[round(x,4) for x in vals]} mono={rr['pb_monotone_down']} fracShort={rr['pb_frac_genes_shortening']:.3f}",flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/trend_stats_v2.tsv",sep="\t",index=False)
pd.concat(cells).to_csv(f"{OUT}/percell_v2.tsv.gz",sep="\t",index=False)
pd.concat(genes).to_csv(f"{OUT}/gene_wdi_v2.tsv.gz",sep="\t",index=False)
pd.DataFrame(pbrows).to_csv(f"{OUT}/pseudobulk_fixedgene.tsv",sep="\t",index=False)
R=pd.DataFrame(rows)
print(R[["mouse","pas_subset","depth","labels","n_cells","n_pas","wdi_rho","wdi_p","wdi_monotone_down","wdi_cliffs","wdi_null_rho_absmax"]].to_string())
