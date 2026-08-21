#!/usr/bin/env python
"""Stage 3 (Laughney) -- label confirmation, implementing scripts/stage3/LABEL_POLICY.md
(pre-registered addendum item 5 + item 6 of manuscript/13_reliability_positioning.md).

Inputs (read-only):
  --curated      laughney_celltypes_curated.parquet   (cell = full '<GSM>_<16nt>' id, laughney_celltype)
  --celltypist   B2_celltypist/percell_labels.csv      (barcode = full id, new_label, conf, imm_label)
  --gex-dir      B2_gex_celltyping/<GSM>_gex_labeled.h5ad   (layers['counts'], obs_names = full id)
  --sample-table scripts/stage3/laughney_samples.tsv  (dataset_id, patient, group)
Outputs (--out-dir): confirmed_labels.tsv, per_gsm_retained.tsv, per_type_summary.tsv, purity_check.tsv,
  manifest.json.  Every join is on the FULL cell id.  Every number is PROVISIONAL until verified.
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import h5py  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

sc.settings.verbosity = 0

# ----------------------------------------------------------------------------- policy (LABEL_POLICY.md s2)
CONF_MIN = 0.5
MIN_CELLS = 20
P25_Q = 25.0

IMMUNE_ACCEPT = {  # curated -> accepted Human_Lung_Atlas majority-voting labels
    "B_cell": {"B cells"},
    "Plasma_IG": {"Plasma cells"},
    "Monocyte": {"Classical monocytes", "Non-classical monocytes"},
    "Macrophage": {"Alveolar macrophages", "Alveolar Mph CCL3+", "Alveolar Mph MT-positive",
                   "Alveolar Mph proliferating", "Monocyte-derived Mph", "Interstitial Mph perivascular"},
    "Dendritic": {"DC1", "DC2", "Migratory DCs", "Plasmacytoid DCs"},
    "NK": {"NK cells"},
    "T_cell": {"CD4 T cells", "CD8 T cells", "T cells proliferating"},
    "Treg": {"CD4 T cells", "T cells proliferating"},
    "Mast": {"Mast cells"},
}
IMM_TREG = {"Regulatory T cells", "Treg(diff)"}   # Immune_All_Low labels that define Treg

NONIMMUNE_MARKERS = {  # LABEL_POLICY.md s3.1 (pre-registered, verbatim)
    "Epithelial_Tumor": ["EPCAM", "KRT8", "KRT18"],
    "Fibroblast": ["COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Pericyte": ["RGS5", "ACTA2", "PDGFRB"],
}
NONIMMUNE = list(NONIMMUNE_MARKERS)

# diagnostic only (imm_agree column): curated -> Immune_All_Low labels of the same lineage
IMM_DIAG = {
    "B_cell": {"B cells", "Memory B cells", "Naive B cells", "Age-associated B cells", "Germinal center B cells",
               "Proliferative germinal center B cells", "Follicular B cells", "Transitional B cells", "Cycling B cells"},
    "Plasma_IG": {"Plasma cells", "Plasmablasts"},
    "Monocyte": {"Classical monocytes", "Non-classical monocytes", "Monocytes", "Cycling monocytes", "Monocyte precursor"},
    "Macrophage": {"Macrophages", "Alveolar macrophages", "Intermediate macrophages", "Erythrophagocytic macrophages",
                   "Mono-mac", "Kupffer cells", "Intestinal macrophages", "Kidney-resident macrophages", "Hofbauer cells"},
    "Dendritic": {"DC", "DC1", "DC2", "DC3", "pDC", "Migratory DCs", "DC precursor", "Cycling DCs", "Transitional DC", "pDC precursor"},
    "NK": {"NK cells", "CD16+ NK cells", "CD16- NK cells", "Cycling NK cells", "Transitional NK"},
    "T_cell": {"Tem/Effector helper T cells", "Tem/Effector helper T cells PD1+", "Tem/Trm cytotoxic T cells",
               "Tem/Temra cytotoxic T cells", "Trm cytotoxic T cells", "Tcm/Naive helper T cells",
               "Tcm/Naive cytotoxic T cells", "Follicular helper T cells", "CRTAM+ gamma-delta T cells",
               "gamma-delta T cells", "Cycling gamma-delta T cells", "Cycling T cells", "MAIT cells",
               "Memory CD4+ cytotoxic T cells", "Type 1 helper T cells", "Type 17 helper T cells", "NKT cells",
               "CD8a/a", "CD8a/b(entry)", "T(agonist)", "Double-negative thymocytes", "Double-positive thymocytes"},
    "Treg": set(IMM_TREG),
    "Mast": {"Mast cells"},
    "Epithelial_Tumor": {"Epithelial cells"},
    "Fibroblast": {"Fibroblasts"},
    "Endothelial": {"Endothelial cells"},
    "Pericyte": set(),
}

# purity check genes (task spec) + extras flagged as such
CHECK_GENES = {  # gene -> (types it is the 'own' marker for, in_prereg_list)
    "MS4A1": ({"B_cell"}, True), "CD79A": ({"B_cell"}, True),
    "CD3E": ({"T_cell", "Treg"}, True),
    "LYZ": ({"Monocyte", "Macrophage"}, True),
    "NKG7": ({"NK"}, True), "GNLY": ({"NK"}, True),
    "JCHAIN": ({"Plasma_IG"}, True),
    "TPSAB1": ({"Mast"}, True),
    "EPCAM": ({"Epithelial", "Epithelial_Tumor"}, True),
    "COL1A2": ({"Fibroblast"}, True),
    "PECAM1": ({"Endothelial"}, True),
    "RGS5": ({"Pericyte"}, False),
    "FOXP3": ({"Treg"}, False),
}


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def finfo(path: str) -> dict:
    st = os.stat(path)
    return {"path": path, "bytes": st.st_size, "mtime": time.strftime("%F %T", time.localtime(st.st_mtime))}


def read_gex(path: str):
    """Element-wise read (the files carry a newer anndata 'null' encoding in uns that 0.11 cannot read)."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        counts = read_elem(f["layers/counts"])
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    counts = counts.tocsr()
    d = counts.data
    if not np.allclose(d, np.round(d)):
        sys.exit(f"{path}: layers['counts'] is not integral")
    return obs.index.astype(str), var.index.astype(str), counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curated", default="/mnt/ssd2/Laugney_Aligned/laughney_ref/laughney_celltypes_curated.parquet")
    ap.add_argument("--celltypist", default="/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs/B1_cohort_full/B2_celltypist/percell_labels.csv")
    ap.add_argument("--gex-dir", default="/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs/B1_cohort_full/B2_gex_celltyping")
    ap.add_argument("--sample-table", default="/mnt/ssd1/Projects/PeakATail_wd/scripts/stage3/laughney_samples.tsv")
    ap.add_argument("--policy", default="/mnt/ssd1/Projects/PeakATail_wd/scripts/stage3/LABEL_POLICY.md")
    ap.add_argument("--out-dir", default="/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/labels")
    ap.add_argument("--min-cells", type=int, default=MIN_CELLS)
    args = ap.parse_args()
    if args.out_dir.startswith("/mnt/ssd2"):
        sys.exit("REFUSING to write under /mnt/ssd2")
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    man = {"status": "PROVISIONAL -- not verified", "policy": {**finfo(args.policy), "md5": md5(args.policy)},
           "inputs": {"curated": {**finfo(args.curated), "md5": md5(args.curated)},
                      "celltypist": {**finfo(args.celltypist), "md5": md5(args.celltypist)},
                      "sample_table": {**finfo(args.sample_table), "md5": md5(args.sample_table)}},
           "versions": {"scanpy": sc.__version__, "anndata": ad.__version__, "numpy": np.__version__,
                        "pandas": pd.__version__, "python": sys.version.split()[0]},
           "params": {"conf_min": CONF_MIN, "min_cells": args.min_cells, "p25_q": P25_Q,
                      "score_genes": {"ctrl_size": 50, "n_bins": 25, "random_state": 0,
                                      "input": "layers['counts'] -> normalize_total(1e4) -> log1p"},
                      "immune_accept": {k: sorted(v) for k, v in IMMUNE_ACCEPT.items()},
                      "imm_treg": sorted(IMM_TREG), "nonimmune_markers": NONIMMUNE_MARKERS},
           "join": "FULL cell id '<GSM>_<16nt>' (parquet cell == csv barcode == h5ad obs_names)"}

    # ---------------------------------------------------------------- tables
    tab = pd.read_csv(args.sample_table, sep="\t", dtype=str)
    gsm2patient = dict(zip(tab["dataset_id"], tab["patient"]))
    gsm2group = dict(zip(tab["dataset_id"], tab["group"]))
    cur = pd.read_parquet(args.curated)
    cur = cur.rename(columns={"laughney_celltype": "curated"})[["cell", "curated"]]
    cur["cell"] = cur["cell"].astype(str)
    if cur["cell"].duplicated().any():
        sys.exit("duplicate cell ids in the curated parquet")
    ct = pd.read_csv(args.celltypist, dtype={"gsm": str, "barcode": str, "new_label": str, "imm_label": str})
    ct = ct.rename(columns={"barcode": "cell", "new_label": "celltypist"})[["cell", "celltypist", "conf", "imm_label"]]
    if ct["cell"].duplicated().any():
        sys.exit("duplicate cell ids in the CellTypist csv")
    if set(ct["cell"]) != set(cur["cell"]):
        sys.exit(f"cell-id sets differ: curated {len(cur)} vs celltypist {len(ct)}; intersection {len(set(ct['cell']) & set(cur['cell']))}")
    df = cur.merge(ct, on="cell", how="inner", validate="one_to_one")
    df["gsm"] = df["cell"].str.split("_", n=1).str[0]
    bad = sorted(set(df["gsm"]) - set(gsm2patient))
    if bad:
        sys.exit(f"GSM prefixes not in sample table: {bad}")
    df["patient"] = df["gsm"].map(gsm2patient)
    df["group"] = df["gsm"].map(gsm2group)
    df["conf"] = df["conf"].astype(float)
    man["n_cells"] = int(len(df))
    man["n_gsm"] = int(df["gsm"].nunique())
    man["n_patients"] = int(df["patient"].nunique())

    # ---------------------------------------------------------------- immune rule (s2)
    is_imm = df["curated"].isin(IMMUNE_ACCEPT)
    acc = df.apply(lambda r: r["celltypist"] in IMMUNE_ACCEPT.get(r["curated"], set()), axis=1)
    conf_ok = df["conf"] >= CONF_MIN
    imm_is_treg = df["imm_label"].isin(IMM_TREG)
    extra_ok = pd.Series(np.where(df["curated"] == "T_cell", ~imm_is_treg,
                                  np.where(df["curated"] == "Treg", imm_is_treg, True)), index=df.index).astype(bool)
    df["confirmed"] = False
    df["reason"] = ""
    m = is_imm
    df.loc[m & ~conf_ok, "reason"] = "conf_lt_0.5"
    df.loc[m & conf_ok & ~acc, "reason"] = "celltypist_label_not_accepted"
    df.loc[m & conf_ok & acc & ~extra_ok & (df["curated"] == "T_cell"), "reason"] = "imm_label_is_treg"
    df.loc[m & conf_ok & acc & ~extra_ok & (df["curated"] == "Treg"), "reason"] = "imm_label_not_treg"
    ok_imm = m & conf_ok & acc & extra_ok
    df.loc[ok_imm, "confirmed"] = True
    df.loc[ok_imm, "reason"] = "ok"
    df["imm_agree"] = df.apply(lambda r: (r["imm_label"] in IMM_DIAG[r["curated"]]) if IMM_DIAG.get(r["curated"]) else np.nan, axis=1)

    # ---------------------------------------------------------------- marker scores + CP10k per GSM
    for t in NONIMMUNE:
        df[f"score_{t}"] = np.nan
    for g in CHECK_GENES:
        df[f"cp10k_{g}"] = np.nan
    df = df.set_index("cell")
    absent = {}
    gex_files = {}
    for gsm in sorted(df["gsm"].unique()):
        path = os.path.join(args.gex_dir, f"{gsm}_gex_labeled.h5ad")
        if not os.path.exists(path):
            sys.exit(f"missing GEX object {path}")
        gex_files[gsm] = finfo(path)
        obs_idx, var_idx, counts = read_gex(path)
        cells = df.index[df["gsm"] == gsm]
        pos = pd.Index(obs_idx).get_indexer(cells)
        if (pos < 0).any():
            sys.exit(f"{gsm}: {(pos < 0).sum()} label cells absent from the GEX object (full-id join)")
        A = ad.AnnData(X=counts.astype(np.float32))
        A.obs_names = obs_idx
        A.var_names = var_idx
        A.var_names_make_unique()
        tot = np.asarray(counts.sum(axis=1)).ravel().astype(float)
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        vset = set(A.var_names)
        absent[gsm] = {}
        for t, genes in NONIMMUNE_MARKERS.items():
            present = [g for g in genes if g in vset]
            absent[gsm][t] = [g for g in genes if g not in vset]
            if not present:
                continue  # stays NaN for every cell of this GSM (s3.4)
            sc.tl.score_genes(A, present, score_name=f"score_{t}", ctrl_size=50, n_bins=25, random_state=0)
            df.loc[cells, f"score_{t}"] = A.obs.loc[cells, f"score_{t}"].astype(float).values
        for g in CHECK_GENES:
            if g in vset:
                j = A.var_names.get_loc(g)
                col = np.asarray(counts[:, j].todense()).ravel().astype(float)
                cp = np.where(tot > 0, col / np.where(tot > 0, tot, 1) * 1e4, np.nan)
                df.loc[cells, f"cp10k_{g}"] = cp[pos]
        print(f"[{gsm}] cells {len(cells)} genes {A.n_vars} absent markers {absent[gsm]}", flush=True)
    man["gex_files"] = gex_files
    man["absent_markers_per_gsm"] = absent
    df = df.reset_index()

    # ---------------------------------------------------------------- non-immune gate (s3.3)
    S = df[[f"score_{t}" for t in NONIMMUNE]].to_numpy(dtype=float)
    own_col = {t: i for i, t in enumerate(NONIMMUNE)}
    passA = np.zeros(len(df), dtype=bool)
    own_score = np.full(len(df), np.nan)
    for t, i in own_col.items():
        rows = np.where(df["curated"].values == t)[0]
        if len(rows) == 0:
            continue
        own = S[rows, i]
        others = np.delete(S[rows], i, axis=1)
        with np.errstate(invalid="ignore"):
            omax = np.nanmax(np.where(np.isnan(others), -np.inf, others), axis=1)  # -inf when all NaN
        a = np.isfinite(own) & (own > omax)
        passA[rows] = a
        own_score[rows] = own
    p25 = {}
    for t in NONIMMUNE:
        rows = np.where((df["curated"].values == t) & passA)[0]
        p25[t] = float(np.percentile(own_score[rows], P25_Q)) if len(rows) else None
    man["p25_thresholds"] = p25
    man["n_pass_argmax_for_p25"] = {t: int(((df["curated"].values == t) & passA).sum()) for t in NONIMMUNE}
    df["p25_threshold"] = df["curated"].map(lambda t: p25.get(t, np.nan) if t in NONIMMUNE else np.nan)
    for t in NONIMMUNE:
        rows = df["curated"].values == t
        thr = p25[t]
        own_nan = rows & ~np.isfinite(own_score)
        df.loc[own_nan, "reason"] = "own_markers_absent"
        notA = rows & np.isfinite(own_score) & ~passA
        df.loc[notA, "reason"] = "not_argmax"
        if thr is None:
            df.loc[rows & passA, "reason"] = "no_p25_defined"
            continue
        below = rows & passA & (own_score < thr)
        df.loc[below, "reason"] = "below_p25"
        okk = rows & passA & (own_score >= thr)
        df.loc[okk, "confirmed"] = True
        df.loc[okk, "reason"] = "ok"
    unknown_types = sorted(set(df["curated"]) - set(IMMUNE_ACCEPT) - set(NONIMMUNE))
    if unknown_types:
        df.loc[df["curated"].isin(unknown_types), "reason"] = "curated_class_not_in_policy"
    man["curated_classes_not_in_policy"] = unknown_types

    # ---------------------------------------------------------------- item 6 relabel + final_label
    df["report_label"] = df["curated"]
    normal_epi = (df["curated"] == "Epithelial_Tumor") & (df["group"] == "Normal")
    df.loc[normal_epi, "report_label"] = "Epithelial"
    df["final_label"] = np.where(df["confirmed"], df["report_label"], "")
    df["confirmed"] = df["confirmed"].astype(bool)

    # ---------------------------------------------------------------- per-GSM retained table
    grp = df.groupby(["gsm", "report_label"], sort=True)
    ret = grp.agg(curated_n=("cell", "size"), confirmed_n=("confirmed", "sum")).reset_index()
    ret["confirmed_n"] = ret["confirmed_n"].astype(int)
    ret["fraction"] = (ret["confirmed_n"] / ret["curated_n"]).round(4)
    ret["retained"] = ret["confirmed_n"] >= args.min_cells
    ret["patient"] = ret["gsm"].map(gsm2patient)
    ret["group"] = ret["gsm"].map(gsm2group)
    ret = ret.rename(columns={"report_label": "celltype"})
    ret = ret[["gsm", "patient", "group", "celltype", "curated_n", "confirmed_n", "fraction", "retained"]]
    ret.to_csv(os.path.join(args.out_dir, "per_gsm_retained.tsv"), sep="\t", index=False)

    # ---------------------------------------------------------------- per-type summary
    rows = []
    for t, g in df.groupby("report_label", sort=True):
        r = ret[(ret["celltype"] == t) & ret["retained"]]
        rows.append({"celltype": t, "curated_n": int(len(g)), "confirmed_n": int(g["confirmed"].sum()),
                     "fraction": round(float(g["confirmed"].mean()), 4),
                     "dropped_n": int((~g["confirmed"]).sum()),
                     "n_gsm_with_curated": int(g["gsm"].nunique()),
                     "n_gsm_retained": int(len(r)), "n_patients_retained": int(r["patient"].nunique()),
                     "patients_retained": ",".join(sorted(r["patient"].unique())),
                     "imm_agree_frac_all": (round(float(g["imm_agree"].mean()), 4) if g["imm_agree"].notna().any() else None),
                     "imm_agree_frac_confirmed": (round(float(g.loc[g["confirmed"], "imm_agree"].mean()), 4) if g.loc[g["confirmed"], "imm_agree"].notna().any() else None),
                     "reasons": json.dumps({k: int(v) for k, v in g["reason"].value_counts().items()})})
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(args.out_dir, "per_type_summary.tsv"), sep="\t", index=False)

    # ---------------------------------------------------------------- purity check (>= min_cells confirmed)
    prow = []
    for (gsm, t), g in df.groupby(["gsm", "report_label"], sort=True):
        gc = g[g["confirmed"]]
        if len(gc) < args.min_cells:
            continue
        for gene, (own_types, prereg) in CHECK_GENES.items():
            b = g[f"cp10k_{gene}"]; a = gc[f"cp10k_{gene}"]
            if b.isna().all():
                prow.append({"gsm": gsm, "patient": gsm2patient[gsm], "celltype": t, "gene": gene,
                             "own_marker": t in own_types, "in_prereg_list": prereg, "n_curated": len(g),
                             "n_confirmed": len(gc), "gene_absent_in_gsm": True})
                continue
            prow.append({"gsm": gsm, "patient": gsm2patient[gsm], "celltype": t, "gene": gene,
                         "own_marker": t in own_types, "in_prereg_list": prereg,
                         "n_curated": int(len(g)), "n_confirmed": int(len(gc)), "gene_absent_in_gsm": False,
                         "before_median_cp10k": round(float(b.median()), 3), "after_median_cp10k": round(float(a.median()), 3),
                         "before_mean_cp10k": round(float(b.mean()), 3), "after_mean_cp10k": round(float(a.mean()), 3),
                         "before_frac_expr": round(float((b > 0).mean()), 4), "after_frac_expr": round(float((a > 0).mean()), 4)})
    pur = pd.DataFrame(prow)
    pur.to_csv(os.path.join(args.out_dir, "purity_check.tsv"), sep="\t", index=False)

    # ---------------------------------------------------------------- per-cell table
    cols = ["cell", "gsm", "patient", "group", "curated", "celltypist", "conf", "imm_label", "imm_agree",
            "confirmed", "reason", "final_label"] + [f"score_{t}" for t in NONIMMUNE] + ["p25_threshold"]
    out = df[cols].copy()
    out["conf"] = out["conf"].round(6)
    for t in NONIMMUNE:
        out[f"score_{t}"] = out[f"score_{t}"].round(5)
    out["p25_threshold"] = out["p25_threshold"].round(5)
    out.to_csv(os.path.join(args.out_dir, "confirmed_labels.tsv"), sep="\t", index=False)

    man["counts"] = {"confirmed_total": int(df["confirmed"].sum()), "dropped_total": int((~df["confirmed"]).sum()),
                     "per_type": summ.set_index("celltype")[["curated_n", "confirmed_n", "fraction", "n_gsm_retained", "n_patients_retained"]].to_dict(orient="index"),
                     "reasons_total": {k: int(v) for k, v in df["reason"].value_counts().items()},
                     "gsm_x_type_retained": int(ret["retained"].sum()), "gsm_x_type_total": int(len(ret))}
    man["elapsed_s"] = round(time.time() - t0, 1)
    man["outputs"] = {k: os.path.join(args.out_dir, k) for k in ("confirmed_labels.tsv", "per_gsm_retained.tsv", "per_type_summary.tsv", "purity_check.tsv")}
    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump(man, f, indent=1)

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
    print("\n=== PROVISIONAL per-type summary ===")
    print(summ.drop(columns=["reasons"]).to_string(index=False))
    print("\n=== PROVISIONAL reasons per type ===")
    for r in rows:
        print(f"{r['celltype']:18s} {r['reasons']}")
    print("\n=== PROVISIONAL P25 thresholds ===", p25)
    print(f"\nconfirmed {df['confirmed'].sum()} / {len(df)} cells; (GSM,type) retained {ret['retained'].sum()} / {len(ret)}; elapsed {man['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
