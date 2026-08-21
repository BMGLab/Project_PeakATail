#!/usr/bin/env python
"""Stage 3 (Laughney) -- per-sample summary of TRUE vs label-shuffle null switch tests.

Reads  <sample-dir>/true/differential/<strategy>_<c1>_vs_<c2>.tsv  and
       <sample-dir>/null/perm_XX/differential/...,  labelled/label_report.json, runtime files.
Writes <sample-dir>/summary.json and summary_per_pair.tsv.  Every number is
PROVISIONAL until the verifier passes (the files say so).

Calibration check (task rule; anti-conservative direction only, manuscript/14 arm B0 as the reference):
    FLAG if pooled null p<0.05 rate > 0.07  OR  > 25 % of null perms carry any q<fdr hit.
"""
import argparse
import glob
import json
import os
import re
import sys

import pandas as pd

P05_MAX = 0.07
PERM_HIT_FRAC_MAX = 0.25


def read_dir(d: str, fdr: float, want_ids: bool = False) -> dict:
    files = sorted(glob.glob(os.path.join(d, "differential", "*_vs_*.tsv")))
    per_pair, n_tests, n_hits, n_p05, n_p01, n_q05 = {}, 0, 0, 0, 0, 0
    pas_ids, genes = set(), set()
    for f in files:
        base = os.path.basename(f)[:-4]
        cols = ("pas_id", "gene_id", "pvalue", "qvalue", "delta_proportion")
        df = pd.read_csv(f, sep="\t", usecols=lambda c: c in cols, dtype={"pas_id": str, "gene_id": str})
        n = len(df)
        q = (df["qvalue"] < fdr)
        hits = int(q.sum())
        eff = int((q & (df["delta_proportion"].abs() >= 0.1)).sum()) if "delta_proportion" in df else None
        per_pair[base] = {"n_tests": n, "n_q_lt_fdr": hits, "n_q_lt_fdr_and_abs_dprop_ge_0.1": eff,
                          "frac_p_lt_0.05": float((df["pvalue"] < 0.05).mean()) if n else None,
                          "n_pas": int(df["pas_id"].nunique()) if "pas_id" in df else None,
                          "n_genes": int(df["gene_id"].nunique()) if "gene_id" in df else None}
        n_tests += n; n_hits += hits; n_q05 += hits
        n_p05 += int((df["pvalue"] < 0.05).sum()); n_p01 += int((df["pvalue"] < 0.01).sum())
        if want_ids and "pas_id" in df:
            pas_ids |= set(df["pas_id"]); genes |= set(df["gene_id"].dropna()) if "gene_id" in df else set()
    out = {"n_pairs": len(files), "n_tests": n_tests, "n_q_lt_fdr": n_hits, "n_p_lt_0.05": n_p05, "n_p_lt_0.01": n_p01,
           "frac_p_lt_0.05": (n_p05 / n_tests) if n_tests else None,
           "frac_p_lt_0.01": (n_p01 / n_tests) if n_tests else None,
           "frac_q_lt_fdr": (n_q05 / n_tests) if n_tests else None,
           "per_pair": per_pair, "done": os.path.exists(os.path.join(d, "DONE.ok"))}
    if want_ids:
        out["n_pas_tested_union"] = len(pas_ids); out["n_genes_tested_union"] = len(genes)
    return out


def time_file(p: str) -> dict:
    out = {}
    if not os.path.exists(p):
        return out
    for line in open(p):
        if "Elapsed (wall clock)" in line:
            out["elapsed_wall"] = line.split("):")[-1].strip()
        if "Maximum resident set size" in line:
            out["max_rss_kb"] = int(line.split(":")[-1].strip())
            out["max_rss_gb"] = round(out["max_rss_kb"] / 1024 / 1024, 2)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--cohort-log", default=None, help="cohort ema_run.log; records whether this sample's BAM raised the LOW POLY(A) CLIP RATE warning")
    a = ap.parse_args()
    D = a.sample_dir
    summ = {"sample": a.sample, "status": "PROVISIONAL -- not verified; not a result", "fdr": a.fdr,
            "calibration_rule": {"null_p05_max": P05_MAX, "perms_with_qhit_frac_max": PERM_HIT_FRAC_MAX}}
    lr = os.path.join(D, "labelled", "label_report.json")
    if os.path.exists(lr):
        rep = json.load(open(lr))
        for k in ("n_cells_h5ad", "n_labelled_cells_in_parquet_for_sample", "n_cells_matched", "n_unconfirmed_cells_for_sample",
                  "n_cells_unlabelled_dropped", "n_cells_final", "n_pas_h5ad", "n_pas_before_tier_filter", "n_pas_after_tier_filter",
                  "n_universe_cohort", "frac_h5ad_pas_in_universe", "universe_rule", "tier_filter", "celltypes_kept",
                  "celltypes_dropped_lt_min_cells", "celltype_counts_matched", "celltype_counts_final", "pasbed_coverage_of_var_pas_id",
                  "n_pas_h5ad_ip_pass", "diag_n_pas_after_tier_filter_with_own_umis_ge2", "diag_n_pas_after_tier_filter_with_own_umis_ge1",
                  "diag_n_pas_after_tier_filter_with_own_umis_0", "n_pairs_expected", "no_pairs"):
            summ[k] = rep.get(k)
    rm = os.path.join(D, "run_manifest.json")
    if os.path.exists(rm):
        m = json.load(open(rm))
        summ["manifest"] = {k: m.get(k) for k in ("mode", "tool_commit", "patient", "group", "universe_md5", "universe_policy_md5", "labels_md5", "switch_extra_args", "pasbed_md5")}
    summ["nopairs"] = os.path.exists(os.path.join(D, "NOPAIRS.ok"))
    if a.cohort_log and os.path.exists(a.cohort_log):
        # The ema log is written by rich at 80 columns: long BAM paths are HARD-WRAPPED across lines (verifier 2026-08-21:
        # the previous whitespace-collapsing regex captured only '.../RER' and matched no GSM). Dropping ALL whitespace makes
        # every path contiguous again (paths carry no whitespace), so the flagged BAM can be identified by its dataset id.
        txt = re.sub(r"\x1b\[[0-9;]*m", "", open(a.cohort_log, errors="replace").read())
        compact = re.sub(r"\s+", "", txt)
        flagged = {m.group(2): float(m.group(1)) for m in
                   re.finditer(r"LOWPOLY\(A\)CLIPRATE:([0-9.]+)%ofthefirst\d+CBreadsin(/.+?\.bam)carryapoly\(A\)softclip", compact)}
        sampled = {m.group(1): float(m.group(2)) for m in
                   re.finditer(r"poly\(A\)clipratefor(/.+?\.bam):([0-9.]+)%of\d+sampledCBreads", compact)}
        bam = json.load(open(rm)).get("bam", "") if os.path.exists(rm) else ""
        def _is_mine(path):
            return a.sample in path or (bam and path == bam)   # SRR-named BAMs (no GSM in the path) match via the manifest
        mine = {k: v for k, v in flagged.items() if _is_mine(k)}
        mine_rate = {k: v for k, v in sampled.items() if _is_mine(k)}
        summ["low_clip_rate_warning"] = {"raised_for_this_bam": bool(mine), "clip_rate_pct": (list(mine.values())[0] if mine else None),
                                         "clip_rate_pct_sampled": (list(mine_rate.values())[0] if mine_rate else None),
                                         "n_bams_flagged_in_cohort_log": len(flagged),
                                         "flagged_bams": sorted(os.path.basename(k) for k in flagged),
                                         "n_bams_with_clip_rate_in_log": len(sampled)}
    summ["runtime_ema_run"] = time_file(os.path.join(D, "runtime_mem_run.txt"))
    summ["runtime_switch_true"] = time_file(os.path.join(D, "runtime_mem_switch_true.txt"))
    summ["runtime_length"] = time_file(os.path.join(D, "runtime_mem_length.txt"))
    summ["true"] = read_dir(os.path.join(D, "true"), a.fdr, want_ids=True)
    nulls = {}
    for d in sorted(glob.glob(os.path.join(D, "null", "perm_*"))):
        nulls[os.path.basename(d)] = read_dir(d, a.fdr)
    summ["null"] = nulls
    done = [v for v in nulls.values() if v["done"]]
    hits = [v["n_q_lt_fdr"] for v in done]
    tot_tests = sum(v["n_tests"] for v in done); tot_p05 = sum(v["n_p_lt_0.05"] for v in done)
    tot_p01 = sum(v["n_p_lt_0.01"] for v in done); tot_q = sum(v["n_q_lt_fdr"] for v in done)
    ns = {"n_perms_done": len(done), "hits_per_perm": hits,
          "mean_hits": (sum(hits) / len(hits)) if hits else None, "max_hits": max(hits) if hits else None,
          "perms_with_any_hit": sum(1 for h in hits if h > 0),
          "frac_perms_with_any_hit": (sum(1 for h in hits if h > 0) / len(hits)) if hits else None,
          "pooled_frac_p_lt_0.05": (tot_p05 / tot_tests) if tot_tests else None,
          "pooled_frac_p_lt_0.01": (tot_p01 / tot_tests) if tot_tests else None,
          "pooled_frac_q_lt_fdr": (tot_q / tot_tests) if tot_tests else None,
          "per_perm_frac_p_lt_0.05": [v["frac_p_lt_0.05"] for v in done],
          "mean_frac_p_lt_0.05": (sum(v["frac_p_lt_0.05"] for v in done if v["frac_p_lt_0.05"] is not None) / len(done)) if done else None,
          "n_null_tests_total": tot_tests}
    if done:
        ns["flag_null_p05_gt_7pct"] = bool(ns["pooled_frac_p_lt_0.05"] is not None and ns["pooled_frac_p_lt_0.05"] > P05_MAX)
        ns["flag_qhits_in_gt_25pct_perms"] = bool(ns["frac_perms_with_any_hit"] > PERM_HIT_FRAC_MAX)
        ns["calibration_flag"] = ns["flag_null_p05_gt_7pct"] or ns["flag_qhits_in_gt_25pct_perms"]
        ns["calibration_verdict"] = "FLAGGED (anti-conservative by the task rule)" if ns["calibration_flag"] else "passes the task rule (conservative/calibrated like manuscript/14 arm B0)"
    summ["null_summary"] = ns
    # per-pair null pooled p<0.05
    pair_null = {}
    for v in done:
        for pair, pp in v["per_pair"].items():
            r = pair_null.setdefault(pair, {"n_tests": 0, "n_p05": 0, "q_hits": 0, "perms": 0, "perms_with_hit": 0})
            r["n_tests"] += pp["n_tests"]; r["n_p05"] += int(round((pp["frac_p_lt_0.05"] or 0) * pp["n_tests"])); r["q_hits"] += pp["n_q_lt_fdr"]
            r["perms"] += 1; r["perms_with_hit"] += int(pp["n_q_lt_fdr"] > 0)
    for tag in ("LENGTH.ok", "LENGTH.skipped", "LENGTH.failed"):
        if os.path.exists(os.path.join(D, tag)):
            summ["length_status"] = tag + ": " + open(os.path.join(D, tag)).read().strip()
    json.dump(summ, open(os.path.join(D, "summary.json"), "w"), indent=1)
    rows = []
    for pair, v in summ["true"]["per_pair"].items():
        pn = pair_null.get(pair, {})
        r = {"sample": a.sample, "pair": pair, "true_n_tests": v["n_tests"], "true_n_pas": v["n_pas"], "true_n_genes": v["n_genes"],
             "true_q_hits": v["n_q_lt_fdr"], "true_q_hits_dprop_ge_0.1": v["n_q_lt_fdr_and_abs_dprop_ge_0.1"], "true_frac_p05": v["frac_p_lt_0.05"],
             "null_perms": pn.get("perms"), "null_pooled_frac_p05": (pn["n_p05"] / pn["n_tests"]) if pn.get("n_tests") else None,
             "null_q_hits_total": pn.get("q_hits"), "null_perms_with_qhit": pn.get("perms_with_hit")}
        for pk, pv in nulls.items():
            pp = pv["per_pair"].get(pair, {})
            r[f"{pk}_n_tests"] = pp.get("n_tests"); r[f"{pk}_q_hits"] = pp.get("n_q_lt_fdr")
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(D, "summary_per_pair.tsv"), sep="\t", index=False)
    print("=== PROVISIONAL SUMMARY", a.sample, "(not verified; not a result) ===")
    for k in ("n_cells_h5ad", "n_labelled_cells_in_parquet_for_sample", "n_cells_matched", "n_cells_final", "celltype_counts_final",
              "n_pas_h5ad", "n_pas_after_tier_filter", "n_universe_cohort", "pasbed_coverage_of_var_pas_id", "celltypes_kept",
              "celltypes_dropped_lt_min_cells", "runtime_switch_true", "nopairs", "low_clip_rate_warning"):
        print(f"{k}: {summ.get(k)}")
    t = summ["true"]
    print(f"TRUE: pairs={t['n_pairs']} tests={t['n_tests']} pas_union={t.get('n_pas_tested_union')} q<{a.fdr}={t['n_q_lt_fdr']} frac_p<0.05={t['frac_p_lt_0.05']}")
    for pk, pv in nulls.items():
        print(f"{pk}: pairs={pv['n_pairs']} tests={pv['n_tests']} q<{a.fdr}={pv['n_q_lt_fdr']} frac_p<0.05={pv['frac_p_lt_0.05']} done={pv['done']}")
    print("NULL SUMMARY:", json.dumps({k: ns.get(k) for k in ("n_perms_done", "hits_per_perm", "perms_with_any_hit", "pooled_frac_p_lt_0.05", "pooled_frac_q_lt_fdr", "calibration_verdict")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
