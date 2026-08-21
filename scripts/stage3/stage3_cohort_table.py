#!/usr/bin/env python
"""Aggregate per-GSM summary.json of the Stage-3 cohort switch stage into SUMMARY_ALL.tsv / summary_all.json
(one row per GSM: cells per type, universe size, pairs, TRUE q<0.05 hits, null p<0.05 rate, hits per perm, flags).
PROVISIONAL until the verifier passes."""
import argparse, glob, json, os, sys
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--sample-table", required=True)
    a = ap.parse_args()
    tab = pd.read_csv(a.sample_table, sep="\t", dtype=str)
    rows, per_pair = [], []
    for _, r in tab.iterrows():
        s = r.dataset_id; d = os.path.join(a.out_root, s)
        row = {"sample": s, "patient": r.patient, "group": r.group,
               "status": "NOPAIRS" if os.path.exists(f"{d}/NOPAIRS.ok") else "DONE" if os.path.exists(f"{d}/DONE.ok") else
                         "FAILED" if os.path.exists(f"{d}/FAILED.err") else "PENDING"}
        sj = f"{d}/summary.json"
        if os.path.exists(sj):
            S = json.load(open(sj)); ns = S.get("null_summary", {}); t = S.get("true", {})
            row.update({"cells_h5ad": S.get("n_cells_h5ad"), "cells_confirmed_matched": S.get("n_cells_matched"), "cells_final": S.get("n_cells_final"),
                        "types_kept": ",".join(S.get("celltypes_kept") or []), "n_types": len(S.get("celltypes_kept") or []),
                        "cells_per_type": json.dumps(S.get("celltype_counts_final")), "types_dropped_lt20": json.dumps(S.get("celltypes_dropped_lt_min_cells")),
                        "pas_h5ad": S.get("n_pas_h5ad"), "universe_size": S.get("n_pas_after_tier_filter"), "universe_cohort": S.get("n_universe_cohort"),
                        "diag_own_umis_ge2": S.get("diag_n_pas_after_tier_filter_with_own_umis_ge2"), "diag_own_umis_0": S.get("diag_n_pas_after_tier_filter_with_own_umis_0"),
                        "pasbed_coverage": S.get("pasbed_coverage_of_var_pas_id"),
                        "pairs": t.get("n_pairs"), "true_tests": t.get("n_tests"), "true_pas_tested": t.get("n_pas_tested_union"), "true_genes_tested": t.get("n_genes_tested_union"),
                        "true_q05_hits": t.get("n_q_lt_fdr"), "true_q05_hits_dprop_ge0.1": sum((v.get("n_q_lt_fdr_and_abs_dprop_ge_0.1") or 0) for v in (t.get("per_pair") or {}).values()) if t.get("per_pair") else None,
                        "true_frac_p05": t.get("frac_p_lt_0.05"),
                        "null_perms": ns.get("n_perms_done"), "null_pooled_p05": ns.get("pooled_frac_p_lt_0.05"), "null_pooled_p01": ns.get("pooled_frac_p_lt_0.01"),
                        "null_pooled_q05": ns.get("pooled_frac_q_lt_fdr"), "null_hits_per_perm": json.dumps(ns.get("hits_per_perm")),
                        "null_perms_with_qhit": ns.get("perms_with_any_hit"), "null_mean_hits": ns.get("mean_hits"),
                        "flag_p05_gt7": ns.get("flag_null_p05_gt_7pct"), "flag_qhit_gt25pct_perms": ns.get("flag_qhits_in_gt_25pct_perms"), "calibration_flag": ns.get("calibration_flag"),
                        "low_clip_warning": (S.get("low_clip_rate_warning") or {}).get("raised_for_this_bam"),
                        "switch_true_wall": (S.get("runtime_switch_true") or {}).get("elapsed_wall"), "switch_true_rss_gb": (S.get("runtime_switch_true") or {}).get("max_rss_gb")})
        rows.append(row)
        pp = f"{d}/summary_per_pair.tsv"
        if os.path.exists(pp):
            per_pair.append(pd.read_csv(pp, sep="\t"))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out_root, "SUMMARY_ALL.tsv"), sep="\t", index=False)
    if per_pair:
        pd.concat(per_pair, ignore_index=True).to_csv(os.path.join(a.out_root, "SUMMARY_PER_PAIR_ALL.tsv"), sep="\t", index=False)
    done = df[df.status == "DONE"]
    agg = {"status": "PROVISIONAL -- not verified", "n_done": int((df.status == "DONE").sum()), "n_nopairs": int((df.status == "NOPAIRS").sum()),
           "n_failed": int((df.status == "FAILED").sum()), "n_pending": int((df.status == "PENDING").sum()),
           "n_flagged": int(done["calibration_flag"].eq(True).sum()) if "calibration_flag" in done else None,
           "flagged": done.loc[done["calibration_flag"].eq(True), "sample"].tolist() if "calibration_flag" in done else None,
           "total_true_q05_hits": int(done["true_q05_hits"].fillna(0).sum()) if "true_q05_hits" in done else None,
           "total_null_tests": None, "rows": rows}
    json.dump(agg, open(os.path.join(a.out_root, "summary_all.json"), "w"), indent=1, default=str)
    cols = [c for c in ("sample", "status", "cells_final", "n_types", "universe_size", "pairs", "true_tests", "true_q05_hits", "true_q05_hits_dprop_ge0.1", "null_pooled_p05", "null_hits_per_perm", "null_perms_with_qhit", "calibration_flag", "low_clip_warning", "switch_true_wall") if c in df.columns]
    print("=== PROVISIONAL cohort switch table (not verified) ===")
    print(df[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
