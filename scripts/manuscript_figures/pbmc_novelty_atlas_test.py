#!/usr/bin/env python3
"""
pbmc_novelty_atlas_test.py -- STAGE 4: are the discriminating sites real PAS?

For every split, the PAS that pass the APA-driven filter are tested for atlas
support (PolyASite 2.0 representative sites, strand-matched, 3'-most base,
<=100 bp) against TWO denominators:
  (i)  the split's own candidate universe (the sites that entered the test) --
       the fair, well-expressed-multi-PAS-gene background;
  (ii) all 275,370 PAS in the matrix.
Fisher exact, two-sided. Also: an AUC of log10 PAS depth separating the two
arms, which quantifies how far each split is a library-depth axis.

Output: results/pbmc_novelty/atlas_support_tests.tsv
"""
import os
os.environ.setdefault("LC_ALL", "C")
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.metrics import roc_auc_score

OUT = Path("/mnt/ssd1/Projects/PeakATail_wd/results/pbmc_novelty")
J = json.load(open(OUT / "novelty_results.json"))
M = pd.read_csv(OUT / "matched_cells.tsv", sep="\t")
BG = J["concordance"]["atlas_frac_all_pas"]
N_ALL = J["concordance"]["n_pas_cells"]

rows = []
tot = {"apa_hit": 0, "apa_n": 0, "cand_hit": 0, "cand_n": 0}
for tag, r in J["apa"].items():
    f = OUT / f"apa_test_{tag}.tsv"
    d = pd.read_csv(f, sep="\t")
    cand_n, cand_hit = len(d), int(d["atlas_supported"].sum())
    sig = d[d["sig"]]
    apa = d[d["apa_driven"]]
    apa_n, apa_hit = len(apa), int(apa["atlas_supported"].sum())
    tot["apa_hit"] += apa_hit; tot["apa_n"] += apa_n
    tot["cand_hit"] += cand_hit; tot["cand_n"] += cand_n
    # Fisher: APA-driven vs the rest of the split's candidate universe
    rest_n = cand_n - apa_n
    rest_hit = cand_hit - apa_hit
    orv, p_cand = fisher_exact([[apa_hit, apa_n - apa_hit],
                                [rest_hit, rest_n - rest_hit]])
    # Fisher: APA-driven vs all PAS in the matrix
    n_all_pas = 275370
    all_hit = int(round(BG * n_all_pas))
    orv2, p_all = fisher_exact([[apa_hit, apa_n - apa_hit],
                                [all_hit, n_all_pas - all_hit]])
    # depth AUC between arms
    ct, a, b = r["celltype"], str(r["pas_A"]), str(r["pas_B"])
    selA = (M["celltype"] == ct) & (M["pas_leiden"].astype(str) == a)
    selB = (M["celltype"] == ct) & (M["pas_leiden"].astype(str) == b)
    y = np.r_[np.ones(selA.sum()), np.zeros(selB.sum())]
    x = np.r_[np.log10(M.loc[selA, "pas_counts"]),
              np.log10(M.loc[selB, "pas_counts"])]
    auc = float(roc_auc_score(y, x))
    xg = np.r_[np.log10(M.loc[selA, "gex_counts"]),
               np.log10(M.loc[selB, "gex_counts"])]
    aucg = float(roc_auc_score(y, xg))
    rows.append({
        "split": tag, "celltype": ct, "nA": r["nA"], "nB": r["nB"],
        "n_candidates": cand_n, "n_significant": len(sig),
        "n_apa_driven": apa_n, "n_apa_driven_atlas": apa_hit,
        "atlas_frac_apa_driven": apa_hit / apa_n if apa_n else np.nan,
        "atlas_frac_candidate_universe": cand_hit / cand_n if cand_n else np.nan,
        "atlas_frac_all_pas": BG,
        "OR_vs_candidates": orv, "p_fisher_vs_candidates": p_cand,
        "OR_vs_all_pas": orv2, "p_fisher_vs_all_pas": p_all,
        "depth_AUC_pas": auc, "depth_AUC_gex": aucg,
        "median_dprop_apa": float(apa["dprop"].abs().median()) if apa_n else np.nan,
    })
df = pd.DataFrame(rows)
# pooled
apa_n, apa_hit = tot["apa_n"], tot["apa_hit"]
cand_n, cand_hit = tot["cand_n"], tot["cand_hit"]
orv, p = fisher_exact([[apa_hit, apa_n - apa_hit],
                       [cand_hit - apa_hit, (cand_n - apa_n) - (cand_hit - apa_hit)]])
n_all_pas = 275370
all_hit = int(round(BG * n_all_pas))
orv2, p2 = fisher_exact([[apa_hit, apa_n - apa_hit],
                         [all_hit, n_all_pas - all_hit]])
df.loc[len(df)] = {
    "split": "POOLED", "celltype": "all", "nA": np.nan, "nB": np.nan,
    "n_candidates": cand_n, "n_significant": np.nan,
    "n_apa_driven": apa_n, "n_apa_driven_atlas": apa_hit,
    "atlas_frac_apa_driven": apa_hit / apa_n,
    "atlas_frac_candidate_universe": cand_hit / cand_n,
    "atlas_frac_all_pas": BG, "OR_vs_candidates": orv,
    "p_fisher_vs_candidates": p, "OR_vs_all_pas": orv2,
    "p_fisher_vs_all_pas": p2, "depth_AUC_pas": np.nan,
    "depth_AUC_gex": np.nan, "median_dprop_apa": np.nan}
df.to_csv(OUT / "atlas_support_tests.tsv", sep="\t", index=False)
print(df.to_string(index=False))
