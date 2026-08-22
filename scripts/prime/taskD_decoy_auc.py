#!/usr/bin/env python3
"""TASK D -- the bar `PRIME_PLAN.md` Change 2 point 2 sets for the score.

The verifier found that on separating genuine Kinnex long-read termini from
internal-priming decoys, **the tool's own `ip_tool_afrac`, inverted, scores AUC
0.7790 -- beating the entire A3 model's 0.7655** (`VERIFY` §6(c)).  A score that
cannot beat one covariate the caller already computes is not a score, it is a
rename.  This computes that comparison for the shipped model.

Pool and labels are the verifier's, unchanged: candidates inside
`tier == 1 & ip_tool_flag == 0` that are within 25 bp of a Kinnex t5 truth
terminus XOR within 25 bp of an internal-priming decoy terminus.  AUC is the
rank statistic (Mann-Whitney), computed with numpy so nothing here needs
scikit-learn either.

Kinnex truth SELECTS NOTHING (manuscript/24 §3.3); this is a read-out.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, "/mnt/ssd1/Projects/PeakATail_wd/tools/pa-prime")
from taskD_eval import Evaluator                          # noqa: E402
from taskD_transfer import (FEATURESETS, TRANSFORM_OVERRIDE,  # noqa: E402
                            design, pool_of, prep)

from ema.countmatrix.pas_score import PasScoreModel        # noqa: E402


def auc(labels: np.ndarray, score: np.ndarray) -> float:
    """Mann-Whitney AUC with mid-ranks for ties."""
    order = np.argsort(score, kind="stable")
    s = score[order]
    r = np.empty(len(s), dtype=np.float64)
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        r[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    ranks = np.empty(len(s))
    ranks[order] = r
    n1 = int(labels.sum())
    n0 = len(labels) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((ranks[labels == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    spec = json.loads(Path(a.model).read_text())
    m = PasScoreModel(spec)
    ev = Evaluator(Path(a.table))
    c = prep(ev.c)
    if "d_kin_t5" not in c:
        sys.exit("this dataset has no Kinnex truth columns")
    s = m.predict_proba_matrix(
        design(Path(a.table), FEATURESETS[spec["featureset"]],
               TRANSFORM_OVERRIDE.get(spec["featureset"]))[0])

    pool = pool_of(c)
    yk = ((c.d_kin_t5 >= 0) & (c.d_kin_t5 <= 25)).values
    yd = ((c.d_kin_decoy >= 0) & (c.d_kin_decoy <= 25)).values
    sel = pool & (yk ^ yd)
    lab = yk[sel].astype(int)
    print(f"pool n={sel.sum():,} (truth {lab.sum():,} / decoy {(1 - lab).sum():,})")

    rows = []
    for name, v in [
        (f"pas_score ({spec['name']}, trained on {spec['trained_on']})", s),
        ("molecules (raw)", c.molecules.values.astype(float)),
        ("clip_reads (raw)", c.clip_reads.values.astype(float)),
        ("ip_tool_afrac (raw, INVERTED)", -c.ip_tool_afrac.values.astype(float)),
        ("a_count_d18 (raw, INVERTED)", -c.a_count_d18.values.astype(float)),
        ("hex_strong (raw)", c.hex_strong.values.astype(float)),
    ]:
        rows.append({"feature": name, "AUC_truth_vs_decoy": auc(lab, v[sel])})
    d = pd.DataFrame(rows).sort_values("AUC_truth_vs_decoy", ascending=False)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out, sep="\t", index=False, float_format="%.6f")
    print(d.to_string(index=False, float_format=lambda x: "%.4f" % x))
    return 0


if __name__ == "__main__":
    sys.exit(main())
