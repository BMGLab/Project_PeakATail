#!/usr/bin/env python
"""Stage 3 (Laughney) -- per-sample summary of TRUE vs label-shuffle null switch tests.

Reads  <sample-dir>/true/differential/<strategy>_<c1>_vs_<c2>.tsv  and
       <sample-dir>/null/perm_XX/differential/...,  label_report.json, runtime files.
Writes <sample-dir>/summary.json and summary_per_pair.tsv.  Every number is
PROVISIONAL until the verifier passes (the files say so).
"""
import argparse
import glob
import json
import os
import re
import sys

import pandas as pd


def read_dir(d: str, fdr: float) -> dict:
    files = sorted(glob.glob(os.path.join(d, "differential", "*_vs_*.tsv")))
    per_pair, n_tests, n_hits, n_p05, n_p01 = {}, 0, 0, 0, 0
    for f in files:
        base = os.path.basename(f)[:-4]
        df = pd.read_csv(f, sep="\t", usecols=lambda c: c in ("pvalue", "qvalue", "delta_proportion"))
        n = len(df)
        q = (df["qvalue"] < fdr)
        hits = int(q.sum())
        eff = int((q & (df["delta_proportion"].abs() >= 0.1)).sum()) if "delta_proportion" in df else None
        per_pair[base] = {"n_tests": n, "n_q_lt_fdr": hits, "n_q_lt_fdr_and_abs_dprop_ge_0.1": eff,
                          "frac_p_lt_0.05": float((df["pvalue"] < 0.05).mean()) if n else None}
        n_tests += n; n_hits += hits
        n_p05 += int((df["pvalue"] < 0.05).sum()); n_p01 += int((df["pvalue"] < 0.01).sum())
    return {"n_pairs": len(files), "n_tests": n_tests, "n_q_lt_fdr": n_hits,
            "frac_p_lt_0.05": (n_p05 / n_tests) if n_tests else None,
            "frac_p_lt_0.01": (n_p01 / n_tests) if n_tests else None,
            "per_pair": per_pair, "done": os.path.exists(os.path.join(d, "DONE.ok"))}


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
    a = ap.parse_args()
    D = a.sample_dir
    summ = {"sample": a.sample, "status": "PROVISIONAL DRY-RUN -- not verified; not a result",
            "fdr": a.fdr}
    lr = os.path.join(D, "labelled", "label_report.json")
    if os.path.exists(lr):
        rep = json.load(open(lr))
        for k in ("n_cells_h5ad", "n_labelled_cells_in_parquet_for_sample", "n_cells_matched",
                  "n_cells_unlabelled_dropped", "n_cells_final", "n_pas_h5ad", "celltypes_kept",
                  "celltypes_dropped_lt_min_cells", "celltype_counts_final", "pasbed_coverage_of_var_pas_id"):
            summ[k] = rep.get(k)
    summ["runtime_ema_run"] = time_file(os.path.join(D, "runtime_mem_run.txt"))
    summ["runtime_switch_true"] = time_file(os.path.join(D, "runtime_mem_switch_true.txt"))
    summ["runtime_length"] = time_file(os.path.join(D, "runtime_mem_length.txt"))
    summ["true"] = read_dir(os.path.join(D, "true"), a.fdr)
    nulls = {}
    for d in sorted(glob.glob(os.path.join(D, "null", "perm_*"))):
        nulls[os.path.basename(d)] = read_dir(d, a.fdr)
    summ["null"] = nulls
    hits = [v["n_q_lt_fdr"] for v in nulls.values() if v["done"]]
    summ["null_summary"] = {"n_perms_done": len(hits),
                            "hits_per_perm": hits,
                            "mean_hits": (sum(hits) / len(hits)) if hits else None,
                            "perms_with_any_hit": sum(1 for h in hits if h > 0),
                            "mean_frac_p_lt_0.05": (sum(v["frac_p_lt_0.05"] for v in nulls.values() if v["done"] and v["frac_p_lt_0.05"] is not None) / len(hits)) if hits else None}
    # length step status
    for tag in ("LENGTH.ok", "LENGTH.skipped", "LENGTH.failed"):
        if os.path.exists(os.path.join(D, tag)):
            summ["length_status"] = tag + ": " + open(os.path.join(D, tag)).read().strip()
    json.dump(summ, open(os.path.join(D, "summary.json"), "w"), indent=1)
    rows = []
    for pair, v in summ["true"]["per_pair"].items():
        r = {"pair": pair, "true_n_tests": v["n_tests"], "true_q_hits": v["n_q_lt_fdr"],
             "true_q_hits_dprop_ge_0.1": v["n_q_lt_fdr_and_abs_dprop_ge_0.1"], "true_frac_p05": v["frac_p_lt_0.05"]}
        for pk, pv in nulls.items():
            pp = pv["per_pair"].get(pair, {})
            r[f"{pk}_n_tests"] = pp.get("n_tests"); r[f"{pk}_q_hits"] = pp.get("n_q_lt_fdr")
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(D, "summary_per_pair.tsv"), sep="\t", index=False)
    print("=== PROVISIONAL DRY-RUN SUMMARY", a.sample, "===")
    for k in ("n_cells_h5ad", "n_labelled_cells_in_parquet_for_sample", "n_cells_matched", "n_cells_final",
              "n_pas_h5ad", "celltypes_kept", "celltypes_dropped_lt_min_cells", "runtime_ema_run",
              "runtime_switch_true", "length_status"):
        print(f"{k}: {summ.get(k)}")
    t = summ["true"]
    print(f"TRUE: pairs={t['n_pairs']} tests={t['n_tests']} q<{a.fdr}={t['n_q_lt_fdr']} frac_p<0.05={t['frac_p_lt_0.05']}")
    for pk, pv in nulls.items():
        print(f"{pk}: pairs={pv['n_pairs']} tests={pv['n_tests']} q<{a.fdr}={pv['n_q_lt_fdr']} frac_p<0.05={pv['frac_p_lt_0.05']} done={pv['done']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
