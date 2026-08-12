#!/usr/bin/env python3
"""Extract per-branch / per-dataset parameter-sweep numbers into a tidy TSV.

Reads ONLY from /mnt/ssd2 (read-only). Writes the tidy table to the scratchpad
so the plotting script can be re-run without re-touching /mnt/ssd2.
"""
import json
import os
import glob
import h5py
import numpy as np
import pandas as pd

ROOT = "/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs/reannotate"
OUT = "/tmp/claude-1000/-mnt-ssd1-Projects-PeakATail-wd/93f9b511-1339-4cbc-9fdb-148fd4b2f08f/scratchpad/sweep_tidy.tsv"

rows = []
for mf in sorted(glob.glob(os.path.join(ROOT, "*/branch_manifest.json"))):
    branch_dir = os.path.dirname(mf)
    branch = os.path.basename(branch_dir)
    m = json.load(open(mf))
    trim, clus, filt = m["trim"], m["clustering"], m["filters"]
    for ds in m["datasets"]:
        dsid = ds["dataset_id"]
        rec = dict(
            branch=branch,
            dataset_id=dsid,
            max_gene_distance=trim["max_gene_distance"],
            utr_multiplier=trim["utr_multiplier"],
            include_extended=bool(trim["include_extended"]),
            method=clus["method"],
            resolution=clus["resolution"],
            n_neighbors=clus["n_neighbors"],
            min_read=filt["min_read"],
            min_cells=filt["min_cells"],
            min_pas_per_cell=filt["min_pas_per_cell"],
            final_cells=ds["final_cells"],
            final_pas=ds["final_pas"],
        )
        # stage_stats gives the denominators (input_pas before annotation/filtering)
        ss_path = os.path.join(branch_dir, "07_clustering", dsid, "stage_stats.json")
        rec["input_pas"] = np.nan
        rec["annotated_pas"] = np.nan
        rec["cells_kept_cb"] = np.nan
        if os.path.exists(ss_path):
            ss = json.load(open(ss_path))
            rec["input_pas"] = ss.get("input_matrix", {}).get("input_pas", np.nan)
            rec["annotated_pas"] = ss.get("annotated", {}).get("annotated_pas", np.nan)
            rec["cells_kept_cb"] = ss.get("cb_filter", {}).get("cells_kept", np.nan)
        # n clusters: read ONLY obs/leiden codes from the h5ad (no X, no matrix load)
        h5 = os.path.join(branch_dir, "07_clustering", dsid, "clusters.h5ad")
        rec["n_clusters"] = np.nan
        rec["n_leiden_categories"] = np.nan
        if os.path.exists(h5):
            try:
                with h5py.File(h5, "r") as f:
                    g = f["obs/leiden"]
                    codes = g["codes"][:]
                    cats = g["categories"][:]
                    rec["n_leiden_categories"] = int(len(cats))
                    rec["n_clusters"] = int(len(np.unique(codes[codes >= 0])))
            except Exception as e:  # noqa: BLE001
                rec["n_clusters_error"] = repr(e)[:200]
        rows.append(rec)

df = pd.DataFrame(rows)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
df.to_csv(OUT, sep="\t", index=False)
print("branches:", df.branch.nunique(), "datasets:", df.dataset_id.nunique(), "rows:", len(df))
print("missing n_clusters rows:", int(df.n_clusters.isna().sum()))
print("categories != observed clusters:",
      int((df.n_leiden_categories != df.n_clusters).sum()))
print(df.groupby("branch")[["final_pas", "final_cells", "n_clusters"]].median().to_string())
print("\ninput_pas constant per dataset across branches?")
chk = df.groupby("dataset_id")["input_pas"].nunique()
print(chk.value_counts().to_string())
print("wrote", OUT)
