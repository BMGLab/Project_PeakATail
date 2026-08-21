#!/usr/bin/env python
"""Stage 3 (Laughney) -- inject GEX-derived cell-type labels into ONE sample's
clusters.h5ad and write the label-shuffle null inputs.

Join rule (manuscript/14 orthogonality + discovery caveat on inject_labels.py):
    obs_name  ==  label-table cell id   (exact string, '<dataset_id>_<16nt barcode>')
No bare-barcode join, no regex.  Cells without a label are DROPPED (not a class).

Label source (mutually exclusive):
    --labels PARQUET            the curated parquet (cell, <label-col>, score)   -- unconfirmed, dry-run only
    --confirmed-labels TSV      confirmed_labels.tsv from stage3_confirm_labels.py (LABEL_POLICY.md):
                                only rows with confirmed == True are used, label = final_label
                                (Epithelial in Normal GSMs, Epithelial_Tumor elsewhere; pairs are
                                defined on final_label and never pooled across the two).
Cell types with < --min-cells cells in this sample are DROPPED and recorded.
Exit code 3 (after writing label_report.json) when fewer than 2 cell types remain (no pair possible).

Cohort mode (design A, UNIVERSE_POLICY.md): --universe-tsv restricts var to the pre-registered universe
(ip_pass AND any_tier1 AND cohort clip_umis_sum >= 2) and asserts that every h5ad PAS is IP-pass;
--pasbed must then be the run's unified/multi_sample_merged.bed (100 % coverage asserted).

Null inputs: --n-perm permutations of obs['celltype'] with
numpy.random.default_rng(seed=k), k = 1..N, permutation of the cell order applied
to the TRUE label vector (identical rule to results/fdr_calibration_v2/build_input.py).

Outputs (--out-dir):
    <sample>.labelled.h5ad          TRUE labels  (obs: celltype, celltype_score, sample, patient, group, stage)
    <sample>.labels.tsv             obs_name -> celltype (TRUE)
    perms/perm_XX.h5ad              permuted labels, same cells / same matrix
    perms/perm_XX_labels.tsv
    label_report.json               counts for the run report (provisional until verified)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import anndata as ad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True, help="<run>/07_clustering/<sample>/clusters.h5ad")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--labels", help="cohort parquet with columns cell, <label-col>, score (UNCONFIRMED curated labels; dry runs only)")
    src.add_argument("--confirmed-labels", help="confirmed_labels.tsv from stage3_confirm_labels.py; uses confirmed==True rows, label = final_label")
    ap.add_argument("--label-col", default="laughney_celltype", help="label column of the parquet (ignored with --confirmed-labels: final_label)")
    ap.add_argument("--allow-unconfirmed", action="store_true",
                    help="required together with --labels: acknowledges that the UNCONFIRMED curated labels are used "
                         "(dry runs only; LABEL_POLICY.md forbids them for any switch statistic that is reported)")
    ap.add_argument("--sample", required=True, help="dataset id == obs_name prefix, e.g. GSM3516665-StageIVprimary")
    ap.add_argument("--sample-table", required=True, help="scripts/stage3/laughney_samples.tsv")
    ap.add_argument("--pasbed", default=None, help="<run>/pasbed.bed; checked for 100%% pas_id coverage")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=5)
    ap.add_argument("--universe-tsv", default=None,
                    help="COHORT MODE (UNIVERSE_POLICY.md): universe.tsv from stage3_build_universe.py; the tested PAS are "
                         "var pas_id with in_universe == True (ip_pass AND any_tier1 AND cohort clip_umis_sum >= 2). "
                         "Requires --tier-filter tier1_ge2. Also asserts every var pas_id is ip_pass.")
    ap.add_argument("--tier-filter", choices=["none", "tier1", "tier1_ge2"], default="none",
                    help="restrict the tested PAS to pasbed col5>0 (clip-supported tier-1) or col5>=2 "
                         "(the pre-registered precision-first default, manuscript/13 s1). Default none = "
                         "every PAS in clusters.h5ad (both tiers). Decide BEFORE looking at results.")
    args = ap.parse_args()

    if args.labels and not args.allow_unconfirmed:
        sys.exit("--labels injects UNCONFIRMED curated labels (LABEL_POLICY.md); pass --confirmed-labels "
                 "confirmed_labels.tsv, or --allow-unconfirmed for an explicitly unconfirmed dry run")
    S = args.sample
    os.makedirs(os.path.join(args.out_dir, "perms"), exist_ok=True)
    label_src = args.confirmed_labels or args.labels
    label_col = "final_label" if args.confirmed_labels else args.label_col
    rep: dict = {"sample": S, "h5ad": args.h5ad, "labels": label_src, "label_col": label_col,
                 "label_source": "confirmed_labels.tsv (LABEL_POLICY.md; confirmed==True only)" if args.confirmed_labels else "curated parquet (UNCONFIRMED)",
                 "min_cells": args.min_cells, "n_perm": args.n_perm,
                 "status": "PROVISIONAL -- not verified"}

    tab = pd.read_csv(args.sample_table, sep="\t", dtype=str)
    row = tab.loc[tab["dataset_id"] == S]
    if len(row) != 1:
        sys.exit(f"sample {S!r} not found exactly once in {args.sample_table}")
    row = row.iloc[0]

    A = ad.read_h5ad(args.h5ad)
    rep["n_cells_h5ad"] = int(A.n_obs)
    rep["n_pas_h5ad"] = int(A.n_vars)
    if "counts" not in A.layers:
        sys.exit("clusters.h5ad has no layers['counts'] -- refuse to test normalised .X")
    L = A.layers["counts"]
    data = L.data if hasattr(L, "data") else np.asarray(L)
    if not np.allclose(data, np.round(data)):
        sys.exit("layers['counts'] is not integral")
    rep["counts_layer_sum"] = int(data.sum())

    # --- prefix check: obs_names must be '<S>_<bc>' or the join is meaningless ---
    prefixes = pd.Series([n.split("_", 1)[0] for n in A.obs_names]).value_counts()
    rep["obs_prefixes"] = {str(k): int(v) for k, v in prefixes.items()}
    if list(prefixes.index) != [S]:
        sys.exit(f"obs_name prefix(es) {dict(prefixes)} != dataset id {S!r}; "
                 "the run YAML dataset id / @RG must equal the label prefix. Not rewriting ids.")

    # --- exact join on the full cell id ---
    if args.confirmed_labels:
        lab = pd.read_csv(args.confirmed_labels, sep="\t", dtype={"cell": str, "final_label": str, "confirmed": str,
                                                                "reason": str, "patient": str, "group": str})
        for c in ("cell", "confirmed", "final_label"):
            if c not in lab.columns:
                sys.exit(f"--confirmed-labels lacks column {c!r}")
        lab_all_s = lab.loc[lab["cell"].str.startswith(S + "_")]
        rep["n_cells_in_confirmed_table_for_sample"] = int(len(lab_all_s))
        conf_mask = lab_all_s["confirmed"].astype(str).str.lower().isin(("true", "1")) & lab_all_s["final_label"].fillna("").ne("")
        rep["n_unconfirmed_cells_for_sample"] = int((~conf_mask).sum())
        rep["unconfirmed_reasons_for_sample"] = {str(k): int(v) for k, v in lab_all_s.loc[~conf_mask, "reason"].value_counts().items()} if "reason" in lab_all_s else None
        lab = lab_all_s.loc[conf_mask].copy()
        if "patient" in lab.columns and (lab["patient"].astype(str) != str(row["patient"])).any():
            sys.exit("patient column in confirmed_labels.tsv disagrees with the sample table")
    else:
        lab = pd.read_parquet(args.labels)
    lab_s = lab.loc[lab["cell"].str.startswith(S + "_")].copy()
    if lab_s["cell"].duplicated().any():
        sys.exit("duplicate cell ids in the label table for this sample")
    lab_s = lab_s.set_index("cell")
    rep["n_labelled_cells_in_parquet_for_sample"] = int(len(lab_s))
    rep["n_labelled_cells_for_sample"] = int(len(lab_s))
    matched = A.obs_names.isin(lab_s.index)
    rep["n_cells_matched"] = int(matched.sum())
    rep["n_cells_unlabelled_dropped"] = int((~matched).sum())
    rep["n_parquet_cells_absent_from_h5ad"] = int(len(lab_s) - matched.sum())
    rep["frac_h5ad_cells_labelled"] = float(matched.mean()) if A.n_obs else 0.0
    if matched.sum() == 0:
        sys.exit("0 cells matched the label table")

    B = A[matched].copy()
    B.obs["celltype"] = lab_s.loc[B.obs_names, label_col].astype(str).values
    if "score" in lab_s.columns:
        B.obs["celltype_score"] = lab_s.loc[B.obs_names, "score"].astype(float).values
    if args.confirmed_labels:
        for c in ("curated", "celltypist", "conf", "reason"):
            if c in lab_s.columns:
                B.obs[f"label_{c}"] = lab_s.loc[B.obs_names, c].values
    B.obs["sample"] = S
    B.obs["patient"] = str(row["patient"])
    B.obs["group"] = str(row["group"])
    B.obs["stage"] = str(row["stage"])

    vc = B.obs["celltype"].value_counts()
    rep["celltype_counts_matched"] = {str(k): int(v) for k, v in vc.items()}
    keep = sorted(vc.index[vc >= args.min_cells].tolist())
    dropped = {str(k): int(v) for k, v in vc.items() if v < args.min_cells}
    rep["celltypes_dropped_lt_min_cells"] = dropped
    rep["celltypes_kept"] = keep
    if len(keep) < 2:
        rep["no_pairs"] = True
        with open(os.path.join(args.out_dir, "label_report.json"), "w") as f:
            json.dump(rep, f, indent=1)
        print(f"NOPAIRS: fewer than 2 cell types with >= {args.min_cells} cells: {dict(vc)}")
        return 3
    B = B[B.obs["celltype"].isin(keep)].copy()
    B.obs["celltype"] = pd.Categorical(B.obs["celltype"].astype(str), categories=keep)
    rep["n_cells_final"] = int(B.n_obs)
    rep["celltype_counts_final"] = {str(k): int(v) for k, v in B.obs["celltype"].value_counts().items()}
    rep["n_pairs_expected"] = len(keep) * (len(keep) - 1) // 2

    # --- pasbed coverage check (the rerun trap: bed must cover 100% of tested pas_ids) ---
    rep["tier_filter"] = args.tier_filter
    if args.tier_filter != "none" and not args.pasbed:
        sys.exit("--tier-filter needs --pasbed (col5 = clip molecules)")
    if args.universe_tsv:
        # ---- cohort mode: universe from UNIVERSE_POLICY.md s1, never from the unified bed's col 5 ----
        if args.tier_filter != "tier1_ge2":
            sys.exit("--universe-tsv encodes TIER_FILTER=tier1_ge2; pass --tier-filter tier1_ge2")
        own = f"umis__{S}"
        U = pd.read_csv(args.universe_tsv, sep="\t", dtype={"new_pas_id": str},
                        usecols=lambda c: c in ("new_pas_id", "ip_pass", "any_tier1", "in_universe", "clip_umis_sum", "clip_umis_max", "n_datasets", own))
        U = U.set_index("new_pas_id")
        for c in ("ip_pass", "any_tier1", "in_universe"):
            U[c] = U[c].astype(str).str.lower().eq("true")
        vid = B.var["pas_id"].astype(str)
        missing = ~vid.isin(U.index)
        if missing.any():
            sys.exit(f"{int(missing.sum())} h5ad pas_ids are absent from universe.tsv -- h5ad and universe come from different runs")
        ipp = U.loc[vid.values, "ip_pass"].values
        rep["n_pas_h5ad_ip_pass"] = int(ipp.sum())
        if not ipp.all():
            sys.exit(f"{int((~ipp).sum())} h5ad pas_ids are IP-flagged (not in posbed/negbed) -- the IP filter did not remove them")
        keep_mask = U.loc[vid.values, "in_universe"].values
        rep["universe_file"] = args.universe_tsv
        rep["universe_rule"] = "ip_pass AND any_tier1 AND cohort clip_umis_sum >= 2 (UNIVERSE_POLICY.md s1)"
        rep["n_universe_cohort"] = int(U.in_universe.sum())
        rep["n_pas_before_tier_filter"] = int(B.n_vars)
        B = B[:, keep_mask].copy()
        rep["n_pas_after_tier_filter"] = int(B.n_vars)
        rep["frac_h5ad_pas_in_universe"] = float(keep_mask.mean())
        rep["n_universe_cohort_absent_from_this_h5ad"] = int(rep["n_universe_cohort"] - B.n_vars)
        if own in U.columns:
            ov = U.loc[B.var["pas_id"].astype(str).values, own].values
            rep["diag_n_pas_after_tier_filter_with_own_umis_ge2"] = int((ov >= 2).sum())
            rep["diag_n_pas_after_tier_filter_with_own_umis_ge1"] = int((ov >= 1).sum())
            rep["diag_n_pas_after_tier_filter_with_own_umis_0"] = int((ov == 0).sum())
        rep["diag_n_pas_after_tier_filter_n_datasets_hist"] = {str(k): int(v) for k, v in pd.Series(U.loc[B.var["pas_id"].astype(str).values, "n_datasets"].values).value_counts().sort_index().items()}
        if B.n_vars == 0:
            sys.exit("universe filter left 0 PAS")
    if args.pasbed:
        bed = pd.read_csv(args.pasbed, sep="\t", header=None, dtype=str, usecols=[3, 4])
        ids = set(bed[3])
        if args.tier_filter != "none" and not args.universe_tsv:
            thr = 1 if args.tier_filter == "tier1" else 2
            keep_ids = set(bed.loc[pd.to_numeric(bed[4], errors="coerce").fillna(0) >= thr, 3])
            vid = B.var["pas_id"].astype(str) if "pas_id" in B.var.columns else pd.Series(B.var_names.astype(str), index=B.var_names)
            mask = vid.isin(keep_ids).values
            rep["n_pas_before_tier_filter"] = int(B.n_vars)
            B = B[:, mask].copy()
            rep["n_pas_after_tier_filter"] = int(B.n_vars)
            if B.n_vars == 0:
                sys.exit("tier filter left 0 PAS")
        var_ids = set(B.var["pas_id"].astype(str)) if "pas_id" in B.var.columns else set(B.var_names.astype(str))
        cov = len(var_ids & ids) / max(1, len(var_ids))
        rep["pasbed_coverage_of_var_pas_id"] = cov
        rep["pasbed_n_ids"] = len(ids)
        if cov < 1.0:
            sys.exit(f"pasbed covers {cov:.4f} of h5ad pas_ids -- wrong bed for this h5ad")

    # --- string hygiene for anndata writeback ---
    B.obs.index = B.obs.index.astype(object)
    B.var.index = B.var.index.astype(object)
    for df in (B.obs, B.var):
        for c in df.columns:
            if str(df[c].dtype).startswith("string"):
                df[c] = df[c].astype(object)

    true_path = os.path.join(args.out_dir, f"{S}.labelled.h5ad")
    B.write_h5ad(true_path)
    B.obs[["celltype", "sample", "patient", "group", "stage"]].to_csv(
        os.path.join(args.out_dir, f"{S}.labels.tsv"), sep="\t")
    rep["labelled_h5ad"] = true_path

    orig = np.asarray(B.obs["celltype"].astype(str).values)
    perms = []
    for k in range(1, args.n_perm + 1):
        rng = np.random.default_rng(seed=k)
        perm = rng.permutation(len(orig))
        B.obs["celltype"] = pd.Categorical(orig[perm], categories=keep)
        assert B.obs["celltype"].value_counts().to_dict() == pd.Series(orig).value_counts().to_dict()
        p = os.path.join(args.out_dir, "perms", f"perm_{k:02d}.h5ad")
        B.write_h5ad(p)
        B.obs[["celltype"]].to_csv(os.path.join(args.out_dir, "perms", f"perm_{k:02d}_labels.tsv"), sep="\t")
        perms.append({"k": k, "seed": k, "h5ad": p,
                      "frac_labels_unchanged": float((orig[perm] == orig).mean())})
    rep["perms"] = perms
    rep["perm_seed_rule"] = "numpy default_rng(seed=k), k=1..N_PERM, permutation of cell order applied to the TRUE celltype vector"

    with open(os.path.join(args.out_dir, "label_report.json"), "w") as f:
        json.dump(rep, f, indent=1)
    print(json.dumps({k: rep.get(k) for k in ("sample", "n_cells_h5ad", "n_labelled_cells_in_parquet_for_sample",
                                              "n_cells_matched", "n_cells_final", "n_pas_h5ad", "n_pas_before_tier_filter",
                                              "n_pas_after_tier_filter", "n_universe_cohort", "pasbed_coverage_of_var_pas_id",
                                              "celltypes_kept", "celltype_counts_final", "celltypes_dropped_lt_min_cells")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
