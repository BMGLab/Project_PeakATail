#!/usr/bin/env python3
"""TASK D -- does the TOOL's ``pas_score`` column equal the offline score?

Every accuracy number in this task is computed offline, from the caller's own
``pas_support.tsv``, because that is the only way to sweep 80 thresholds without
80 runs.  That is only legitimate while the offline computation and the tool's
own are the same number.  This checks it on a real run: it recomputes the score
from the run's sidecar with `ema.countmatrix.pas_score` and compares, row by
row, against the `pas_score` column the tool wrote during that run.

They should agree to the sidecar's own printing precision (``%.6f``).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/ssd1/Projects/PeakATail_wd/tools/pa-prime")
from ema.countmatrix.pas_score import load_model          # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--support", required=True, help="a run's pas_support.tsv WITH pas_score")
    ap.add_argument("--model", default="prime1")
    a = ap.parse_args()

    m = load_model(a.model)
    sup = pd.read_csv(a.support, sep="\t", dtype={"pas_id": str}, na_values=["NA"])
    if "pas_score" not in sup.columns:
        sys.exit("this sidecar has no pas_score column (run with --pas-score calibrated)")
    cols = {f: pd.to_numeric(sup[f], errors="coerce").values.astype(float)
            for f in m.features}
    mine = m.predict_proba(cols)
    theirs = sup.pas_score.values.astype(float)
    ok = ~np.isnan(theirs)
    d = np.abs(mine[ok] - theirs[ok])
    print(f"rows {len(sup):,}  scored {int(ok.sum()):,}  NA {int((~ok).sum()):,}")
    print(f"max |offline - tool| = {d.max():.3e}   mean {d.mean():.3e}")
    print(json.dumps({"rows": int(len(sup)), "scored": int(ok.sum()),
                      "max_abs_diff": float(d.max())}))
    assert d.max() <= 5e-7, "the tool's pas_score and the offline score disagree"
    print("AGREE to the sidecar's %.6f printing precision")
    return 0


if __name__ == "__main__":
    sys.exit(main())
