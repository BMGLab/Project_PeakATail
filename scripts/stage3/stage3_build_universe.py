#!/usr/bin/env python
"""Stage 3 (Laughney, cohort run / design A) -- build the pre-registered tested-PAS universe.

Implements scripts/stage3/UNIVERSE_POLICY.md s1 (read it first):
    in_universe(u) = ip_pass(u) AND any_tier1(u) AND clip_umis_sum(u) >= --min-molecules (2)
computed per UNIFIED PAS u from the per-dataset support sidecars (peakcalling/<ds>_0.{pos,neg}.support.tsv,
keyed by OLD pasnumber) joined through unified/multi_sample_pas_mapping.tsv, with ip_pass = presence in the
post-IP-filter unified set (run/posbed.bed + run/negbed.bed).

It also VERIFIES (and records) how the unified bed's col 5 is carried (one member's clip_umis, not a sum),
that tier-2 members carry 0 clip molecules, that every support row maps to a unified id and vice versa,
and that every dataset in the sample table contributed.

Outputs (--out-dir): universe.tsv (one row per unified PAS, all columns of s1 + s2 diagnostics, per-dataset
umis__<ds>), universe_ids.txt (in_universe ids), universe.bed (BED6 subset of the unified bed),
universe_summary.json. Every number is PROVISIONAL until the verifier passes.
"""
import argparse
import glob
import hashlib
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd


def md5(p: str) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="<cohort_run>/run (ema output dir of the cohort run)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-table", required=True, help="scripts/stage3/laughney_samples.tsv")
    ap.add_argument("--policy", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "UNIVERSE_POLICY.md"))
    ap.add_argument("--min-molecules", type=int, default=2)
    a = ap.parse_args()
    if a.out_dir.startswith("/mnt/ssd2"):
        sys.exit("REFUSING to write under /mnt/ssd2")
    os.makedirs(a.out_dir, exist_ok=True)
    run = a.run
    t0 = time.time()
    S: dict = {"status": "PROVISIONAL -- not verified", "run": run, "min_molecules": a.min_molecules,
               "policy": {"path": a.policy, "md5": md5(a.policy), "mtime": time.strftime("%F %T", time.localtime(os.stat(a.policy).st_mtime))},
               "built_at": time.strftime("%F %T")}

    uni_p = f"{run}/unified/multi_sample_merged.bed"
    map_p = f"{run}/unified/multi_sample_pas_mapping.tsv"
    for p in (uni_p, map_p, f"{run}/posbed.bed", f"{run}/negbed.bed"):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
    S["inputs"] = {os.path.relpath(p, run): {"md5": md5(p), "bytes": os.path.getsize(p)} for p in (uni_p, map_p, f"{run}/posbed.bed", f"{run}/negbed.bed")}

    uni = pd.read_csv(uni_p, sep="\t", header=None, dtype=str, names=["chrom", "start", "end", "new_pas_id", "score", "strand"])
    if uni.new_pas_id.duplicated().any():
        sys.exit("duplicate ids in unified bed")
    uni["score"] = pd.to_numeric(uni.score, errors="coerce")
    S["n_unified"] = int(len(uni))

    mp = pd.read_csv(map_p, sep="\t", dtype=str)
    mp["old_pasnumber"] = mp.old_pasnumber.astype(int)
    S["n_mapping_rows"] = int(len(mp))
    if not set(mp.new_pas_id) <= set(uni.new_pas_id):
        sys.exit("mapping refers to ids absent from the unified bed")

    tab = pd.read_csv(a.sample_table, sep="\t", dtype=str)
    expected_ds = sorted(tab.dataset_id)
    sup_list = []
    for f in sorted(glob.glob(f"{run}/peakcalling/*.support.tsv")):
        m = re.match(r"(.+)_(\d+)\.(pos|neg)\.support\.tsv$", os.path.basename(f))
        if not m:
            sys.exit(f"unexpected sidecar name {f}")
        d = pd.read_csv(f, sep="\t")
        d["dataset_id"] = m.group(1)
        d["strand"] = "+" if m.group(3) == "pos" else "-"
        sup_list.append(d)
    sup = pd.concat(sup_list, ignore_index=True)
    got_ds = sorted(sup.dataset_id.unique())
    S["datasets_in_sidecars"] = got_ds
    if got_ds != expected_ds:
        sys.exit(f"sidecar datasets {got_ds} != sample table {expected_ds}")
    if sorted(mp.dataset_id.unique()) != expected_ds:
        sys.exit("mapping datasets != sample table")
    # tier-2 == 0 clip molecules (policy s0); asserted
    t2 = sup.loc[sup.tier == 2, "clip_umis"]
    S["tier2_rows"] = int(len(t2)); S["tier2_rows_with_clip_umis_gt0"] = int((t2 > 0).sum())
    t1 = sup.loc[sup.tier == 1, "clip_umis"]
    S["tier1_rows"] = int(len(t1)); S["tier1_rows_with_clip_umis_0"] = int((t1 == 0).sum())
    if S["tier2_rows_with_clip_umis_gt0"]:
        print("WARNING: tier-2 rows with clip_umis > 0 exist; any_tier1 is then NOT implied by the molecule count (rule still applied as written)")

    j = sup.merge(mp, left_on=["dataset_id", "strand", "pas_id"], right_on=["dataset_id", "strand", "old_pasnumber"], how="outer", indicator=True)
    S["support_rows"] = int(len(sup))
    S["support_rows_without_mapping"] = int((j._merge == "left_only").sum())
    S["mapping_rows_without_support"] = int((j._merge == "right_only").sum())
    if S["support_rows_without_mapping"] or S["mapping_rows_without_support"]:
        sys.exit(f"support/mapping mismatch: {S['support_rows_without_mapping']} support rows unmapped, {S['mapping_rows_without_support']} mapping rows without support")
    j = j.drop(columns=["_merge"])

    g = j.groupby("new_pas_id")
    U = pd.DataFrame({
        "n_members": g.size(),
        "n_datasets": g.dataset_id.nunique(),
        "n_tier1_members": g.tier.apply(lambda t: int((t == 1).sum())),
        "clip_umis_sum": g.clip_umis.sum(),
        "clip_umis_max": g.clip_umis.max(),
        "clip_umis_f3844_sum": g.clip_umis_f3844.sum(),
        "window_reads_sum": g.window_reads.sum(),
    })
    U["any_tier1"] = U.n_tier1_members > 0
    U["all_tier1"] = U.n_tier1_members == U.n_members
    per_ds = j.pivot_table(index="new_pas_id", columns="dataset_id", values="clip_umis", aggfunc="sum", fill_value=0)
    per_ds.columns = [f"umis__{c}" for c in per_ds.columns]
    U = uni.set_index("new_pas_id").join(U, how="left").join(per_ds, how="left")
    if U.n_members.isna().any():
        sys.exit("unified ids without any member")
    for c in ("n_members", "n_datasets", "n_tier1_members", "clip_umis_sum", "clip_umis_max", "clip_umis_f3844_sum", "window_reads_sum"):
        U[c] = U[c].astype(int)
    for c in per_ds.columns:
        U[c] = U[c].fillna(0).astype(int)

    ip = set(pd.read_csv(f"{run}/posbed.bed", sep="\t", header=None, dtype=str, usecols=[3])[3]) | \
         set(pd.read_csv(f"{run}/negbed.bed", sep="\t", header=None, dtype=str, usecols=[3])[3])
    U["ip_pass"] = U.index.isin(ip)
    S["ip_pass_ids"] = int(len(ip)); S["ip_pass_in_unified"] = int(U.ip_pass.sum())
    if not ip <= set(U.index):
        sys.exit("posbed/negbed contain ids absent from the unified bed")
    pf = f"{run}/04_pas_gene_assignment/peak_filters_stats.json"
    if os.path.exists(pf):
        s = json.load(open(pf))
        S["peak_filters_stats"] = {k: s.get(k) for k in ("ip_filter", "ip_filter_mode", "n_ip_flagged", "ip_flag_rate")}
        S["ip_flagged_consistent"] = (S["n_unified"] - S["ip_pass_in_unified"]) == s.get("n_ip_flagged")
    ann_p = f"{run}/annotatedpas.bed"
    if os.path.exists(ann_p):
        ann_ids = set(pd.read_csv(ann_p, sep="\t", header=None, dtype=str, usecols=[3])[3])
        U["in_annotated"] = U.index.isin(ann_ids)
        S["annotated_ids"] = int(len(ann_ids)); S["annotated_not_ip_pass"] = int(len(ann_ids - ip))
    else:
        U["in_annotated"] = False

    U["in_universe"] = U.ip_pass & U.any_tier1 & (U.clip_umis_sum >= a.min_molecules)

    # how col5 is carried: one member's clip_umis (bedtools -o first), not the sum
    member_vals = j.groupby("new_pas_id").clip_umis.apply(set)
    score_is_member = np.array([(sc in member_vals.get(i, set())) if pd.notna(sc) else False for i, sc in zip(U.index, U.score)])
    multi = U.n_members > 1
    S["col5_check"] = {
        "score_equals_some_member_clip_umis_frac": float(score_is_member.mean()),
        "score_equals_sum_frac_all": float((U.score == U.clip_umis_sum).mean()),
        "multi_member_pas": int(multi.sum()),
        "multi_member_score_ne_sum_frac": float((U.loc[multi, "score"] != U.loc[multi, "clip_umis_sum"]).mean()) if multi.any() else None,
        "pas_with_score_ge2": int((U.score >= 2).sum()), "pas_with_sum_ge2": int((U.clip_umis_sum >= a.min_molecules).sum()),
    }
    S["counts"] = {
        "n_unified": int(len(U)), "ip_pass": int(U.ip_pass.sum()), "any_tier1": int(U.any_tier1.sum()),
        "all_tier1": int(U.all_tier1.sum()), "sum_ge_min": int((U.clip_umis_sum >= a.min_molecules).sum()),
        "ip_pass_and_any_tier1": int((U.ip_pass & U.any_tier1).sum()),
        "in_universe": int(U.in_universe.sum()),
        "in_universe_and_annotated": int((U.in_universe & U.in_annotated).sum()),
        "in_universe_all_tier1": int((U.in_universe & U.all_tier1).sum()),
        "in_universe_max_ge_min": int((U.in_universe & (U.clip_umis_max >= a.min_molecules)).sum()),
        "n_datasets_hist_in_universe": {int(k): int(v) for k, v in U.loc[U.in_universe, "n_datasets"].value_counts().sort_index().items()},
        "n_members_hist_in_universe": {int(k): int(v) for k, v in U.loc[U.in_universe, "n_members"].value_counts().sort_index().items()},
    }
    S["per_dataset"] = {}
    for ds in expected_ds:
        c = f"umis__{ds}"
        S["per_dataset"][ds] = {"caller_pas": int((j.dataset_id == ds).sum()),
                                "caller_tier1": int(((j.dataset_id == ds) & (j.tier == 1)).sum()),
                                "caller_clip_umis_ge_min": int(((j.dataset_id == ds) & (j.clip_umis >= a.min_molecules)).sum()),
                                "unified_with_own_member": int((U[c] > 0).sum() + int(((U[c] == 0) & U.index.isin(j.loc[j.dataset_id == ds, "new_pas_id"])).sum())),
                                "in_universe_with_own_umis_ge_min": int((U.in_universe & (U[c] >= a.min_molecules)).sum()),
                                "in_universe_with_own_umis_ge1": int((U.in_universe & (U[c] >= 1)).sum())}

    out_tsv = os.path.join(a.out_dir, "universe.tsv")
    U.reset_index().to_csv(out_tsv, sep="\t", index=False)
    ids = U.index[U.in_universe]
    with open(os.path.join(a.out_dir, "universe_ids.txt"), "w") as f:
        f.write("\n".join(ids) + "\n")
    U.loc[U.in_universe].reset_index()[["chrom", "start", "end", "new_pas_id", "clip_umis_sum", "strand"]].to_csv(
        os.path.join(a.out_dir, "universe.bed"), sep="\t", index=False, header=False)
    S["outputs"] = {"universe_tsv": out_tsv, "universe_tsv_md5": md5(out_tsv), "rows": int(len(U)),
                    "universe_ids_txt": os.path.join(a.out_dir, "universe_ids.txt"), "n_ids": int(len(ids))}
    S["elapsed_s"] = round(time.time() - t0, 1)
    json.dump(S, open(os.path.join(a.out_dir, "universe_summary.json"), "w"), indent=1, default=str)
    print(json.dumps({k: S[k] for k in ("n_unified", "counts", "col5_check", "tier2_rows", "tier2_rows_with_clip_umis_gt0", "support_rows", "ip_pass_in_unified", "elapsed_s")}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
