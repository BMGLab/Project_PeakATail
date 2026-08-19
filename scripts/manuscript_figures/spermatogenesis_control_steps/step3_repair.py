#!/usr/bin/env python3
"""Step 3: rebuild the PAS x cell matrix with rows correctly keyed by pas_id.

Bug being worked around (verified at source):
  ema/matrixfilter.py::make_dataframe returns `sparse_csc` with ALL
  max(pas_id) rows but `pas_ids = np.sort(np.unique(coo.row+1))` --
  i.e. only the rows that carry >=1 non-zero.  ema/annotate/annotate.py
  then builds `pas_id_to_row = {pid: i for i, pid in enumerate(pas_ids)}`
  and slices `sparse_matrix[keep_rows]`.  Whenever any PAS row is all-zero
  the position index and the pas_id index diverge, so PAS `p` receives the
  count vector of the PAS `p - (#all-zero ids < p)`.
Repair: read filterdmatrix.mtx (row index IS the pas_id) directly.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[v]="8"
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, json
RUN="/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/gse104556/peakatail/{m}/run"
SS="/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/gse104556/starsolo/{S}/Solo.out/Gene/filtered"
OUT="/mnt/ssd0/emaout/peakatail_benchmark/scratch/sperm_ctl"
stats={}
for m,S in [("mouse1","Mouse1_scRNAseq"),("mouse2","Mouse2_scRNAseq")]:
    r=RUN.format(m=m)
    hdr=open(f"{r}/filterdmatrix.mtx").readlines()[1].split()
    nrow,ncol,nnz=int(hdr[0]),int(hdr[1]),int(hdr[2])
    d=pd.read_csv(f"{r}/filterdmatrix.mtx",sep=" ",skiprows=2,header=None,names=["pas","cb","cnt"],
                  dtype={"pas":np.int32,"cb":np.int32,"cnt":np.int32})
    cb=pd.read_csv(f"{r}/filtered_cb.tsv",header=None)[0].astype(str).values
    assert len(cb)==ncol, (len(cb),ncol)
    ann=pd.read_csv(f"{r}/annotatedpas.bed",sep="\t",header=None,
        names=["chrom","start","end","pas_id","gene_id","symbol","strand","score","tier"],
        dtype={"chrom":str,"symbol":str,"gene_id":str,"strand":str,"tier":str})
    ann=ann.drop_duplicates("pas_id").sort_values("pas_id")
    kept=ann["pas_id"].values
    row_of=np.full(nrow+1,-1,dtype=np.int64); row_of[kept]=np.arange(len(kept))
    ok=(d["pas"].values<=nrow)&(row_of[d["pas"].values]>=0)
    rr=row_of[d["pas"].values[ok]]; cc=d["cb"].values[ok]-1; vv=d["cnt"].values[ok]
    M=sp.csr_matrix((vv,(cc,rr)),shape=(ncol,len(kept)))    # cells x PAS
    A=ad.AnnData(X=M,obs=pd.DataFrame(index=pd.Index([c.split("_",1)[-1] for c in cb],name="cb")),
                 var=ann.set_index(ann["pas_id"].astype(str)).rename_axis("pas"))
    A.obs["cb_full"]=cb
    A.obs["umi_total"]=np.asarray(M.sum(1)).ravel().astype(int)
    ssbc=pd.read_csv(f"{SS.format(S=S)}/barcodes.tsv",header=None)[0].astype(str).values
    A.obs["starsolo_cell"]=A.obs_names.isin(set(ssbc))
    A.write_h5ad(f"{OUT}/{m}_pas_repaired.h5ad")
    stats[m]=dict(mtx_rows_declared=nrow,nonzero_rows=int(d['pas'].nunique()),
                  empty_rows=int(nrow-d['pas'].nunique()),kept_pas=int(len(kept)),
                  cells=int(ncol),starsolo_cells=int(A.obs['starsolo_cell'].sum()),
                  umi_in_kept_pas=int(vv.sum()))
    print(m,stats[m],flush=True)
json.dump(stats,open(f"{OUT}/step3_stats.json","w"),indent=1)
