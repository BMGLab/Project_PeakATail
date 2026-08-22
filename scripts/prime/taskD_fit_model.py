#!/usr/bin/env python3
"""TASK D step 2 -- fit the SHIPPED model offline, export it as numpy constants,
and choose its default threshold.  Run only if `taskD_transfer.py` passed.

Non-negotiable 3 of `tools/pa-prime/PRIME_PLAN.md`: scikit-learn never enters the
tool's runtime import path.  It is used HERE, offline, and the artefact that
ships is a JSON of constants that `ema.countmatrix.pas_score` evaluates with
numpy alone.  This script asserts the exported constants reproduce scikit-learn's
own probabilities to 1e-9 on every training row before it writes anything.

THRESHOLD RULE -- DECLARED HERE, BEFORE ANY THRESHOLD'S EFFECT WAS LOOKED AT
---------------------------------------------------------------------------
`manuscript/24_prime_preregistration.md` §3.3 Rule T: a shipped default
threshold is fixed on **GSE104556 testis mouse 1** and evaluated unchanged on
PBMC and mouse 2; Kinnex long-read truth selects nothing, ever.  The task brief
adds: *choose the cut from a decoy/null set, never by tuning the reported
metric.*  Two candidate rules, in preference order, both computed on mouse 1
alone:

  **T1 (primary) -- the calibrated decision boundary, p >= 0.50.**  "Keep a site
  when the model says it is more likely than not to be a real polyadenylation
  site."  It reads no precision, no recall and no call count.  It is only
  meaningful because the model is calibrated, so its calibration on mouse 1 is
  reported next to it.

  **T2 (fallback, used only if T1's mouse-1 call count falls outside
  [0.5x, 2.0x] of the incumbent rule's) -- the decoy-referenced cut**: the 95th
  percentile of the score over mouse 1's tier-1 candidates that the
  internal-priming veto FLAGS.  Those are a decoy population defined by genomic
  sequence, not by any atlas, so the cut is chosen against a null rather than
  against a metric.

Whichever fires is recorded in the exported JSON together with the numbers that
selected it, and mouse 1's own post-hoc P/R are labelled
`THRESHOLD-FITTING SET -- not evidence`.

CALIBRATION (task brief item 4)
-------------------------------
The probability ships as a sidecar column whatever the threshold does, so its
calibration is a result in its own right.  Reliability tables are written
against BOTH truths -- the training label (the atlas) and, on PBMC, the Kinnex
long-read truth the model never saw -- with the expected calibration error of
each.  The gap between them is the publishable part.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from taskD_eval import Evaluator                        # noqa: E402
from taskD_transfer import (FEATURESETS, TRANSFORM_OVERRIDE,  # noqa: E402
                            build_matrix, design, gb, kind_of, lr, pool_of, prep)

T1_P = 0.50
T2_DECOY_PCT = 95.0
T1_COUNT_BAND = (0.5, 2.0)


# ------------------------------------------------------------------ export --
def transform_spec(cols, overrides=None):
    """The per-feature transform, as data rather than as code, so the tool and
    this script cannot drift apart."""
    return [kind_of(k, overrides) for k in cols]


def export_linear(pipe, cols, ov=None) -> dict:
    sc, clf = pipe.named_steps["sc"], pipe.named_steps["lr"]
    # standardise-then-linear collapses to one linear form:
    #   z = b + sum_j w_j (x_j - mu_j) / s_j
    w = clf.coef_[0] / sc.scale_
    b = float(clf.intercept_[0] - np.dot(clf.coef_[0], sc.mean_ / sc.scale_))
    return {"family": "linear", "features": list(cols),
            "transform": transform_spec(cols, ov),
            "coef": [float(x) for x in w], "intercept": b}


def export_hgb(model, cols, ov=None) -> dict:
    """Flatten a HistGradientBoostingClassifier into plain arrays.

    Only the raw-threshold (non-binned, non-categorical) predictor form is
    exported; the assertion at the end of `main` is what guarantees the export
    is faithful, so an unsupported node shape fails loudly rather than silently.
    """
    feat, thr, left, right, isleaf, value, offs = [], [], [], [], [], [], [0]
    for stage in model._predictors:
        assert len(stage) == 1, "multiclass HGB is not supported here"
        nodes = stage[0].nodes
        assert not nodes["is_categorical"].any(), "categorical splits not supported"
        for n in nodes:
            feat.append(int(n["feature_idx"]))
            thr.append(float(n["num_threshold"]))
            left.append(int(n["left"]))
            right.append(int(n["right"]))
            isleaf.append(int(n["is_leaf"]))
            value.append(float(n["value"]))
        offs.append(len(feat))
    return {"family": "hgb", "features": list(cols),
            "transform": transform_spec(cols, ov),
            "baseline": float(np.ravel(model._baseline_prediction)[0]),
            "node_feature": feat, "node_threshold": thr, "node_left": left,
            "node_right": right, "node_is_leaf": isleaf, "node_value": value,
            "tree_offset": offs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--train", default="mouse1")
    ap.add_argument("--family", required=True, choices=("lr", "gb"))
    ap.add_argument("--featureset", required=True, choices=tuple(FEATURESETS))
    ap.add_argument("--name", required=True, help="model name, e.g. prime1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    T, out = Path(a.tables), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cols = FEATURESETS[a.featureset]
    ov = TRANSFORM_OVERRIDE.get(a.featureset)

    ev = Evaluator(T / a.train)
    c = prep(ev.c)
    y = ((c.d_atlas >= 0) & (c.d_atlas <= 25)).values.astype(int)
    X = design(T / a.train, cols, ov)[0]
    model = (lr() if a.family == "lr" else gb())
    model.fit(X, y)
    p_sk = model.predict_proba(X)[:, 1]

    spec = (export_linear(model, cols, ov) if a.family == "lr"
            else export_hgb(model, cols, ov))
    spec.update({"name": a.name, "trained_on": a.train, "featureset": a.featureset,
                 "label": "PolyASite 2.0 within 25 bp, same strand",
                 "n_train": int(len(y)), "n_positive": int(y.sum())})

    # the export must reproduce scikit-learn exactly, or nothing ships
    sys.path.insert(0, "/mnt/ssd1/Projects/PeakATail_wd/tools/pa-prime")
    from ema.countmatrix.pas_score import PasScoreModel      # noqa: E402
    m = PasScoreModel(spec)
    p_np = m.predict_proba_matrix(X)
    err = float(np.max(np.abs(p_np - p_sk)))
    print(f"[export] {a.family}/{a.featureset}: max |numpy - sklearn| = {err:.3e}")
    assert err < 1e-9, f"export is not faithful: max abs error {err}"

    # ---- threshold, on the training dataset only ---------------------------
    pool = pool_of(c)
    mol = pd.to_numeric(c.molecules, errors="coerce").fillna(0).values
    n_rule = int((pool & (mol >= 2)).sum())
    n_t1 = int((pool & (p_sk >= T1_P)).sum())
    ipflag = pd.to_numeric(c.get("ip_tool_flag"), errors="coerce").fillna(0).values
    decoy_pool = (c.tier.values == 1) & (ipflag == 1)
    t2 = float(np.percentile(p_sk[decoy_pool], T2_DECOY_PCT)) if decoy_pool.sum() else float("nan")
    lo, hi = T1_COUNT_BAND[0] * n_rule, T1_COUNT_BAND[1] * n_rule
    use_t1 = lo <= n_t1 <= hi
    thr = T1_P if use_t1 else t2
    rule = "T1_p>=0.50" if use_t1 else f"T2_decoy_p{T2_DECOY_PCT:g}"
    print(f"[threshold] rule n(k>=2) = {n_rule:,}; T1 p>=0.50 -> n {n_t1:,} "
          f"(band {lo:,.0f}-{hi:,.0f}); T2 decoy p{T2_DECOY_PCT:g} = {t2:.4f} "
          f"-> n {int((pool & (p_sk >= t2)).sum()):,}; CHOSEN {rule} = {thr:.6f}")
    spec.update({"threshold": float(thr), "threshold_rule": rule,
                 "threshold_evidence": {
                     "train_n_rule_k2": n_rule, "train_n_T1": n_t1,
                     "T1_count_band": list(T1_COUNT_BAND),
                     "T2_decoy_percentile": T2_DECOY_PCT,
                     "T2_value": None if np.isnan(t2) else t2,
                     "n_decoy_pool": int(decoy_pool.sum())}})

    (out / f"pas_score_model_{a.name}.json").write_text(json.dumps(spec) + "\n")
    print(f"[out] {out}/pas_score_model_{a.name}.json")

    # ---- calibration against every available truth -------------------------
    rows = []
    for ds in sorted(p.name for p in T.iterdir() if (p / "candidates.parquet").exists()):
        e2 = Evaluator(T / ds)
        c2 = prep(e2.c)
        s = m.predict_proba_matrix(design(T / ds, cols, ov)[0])
        pool2 = pool_of(c2)
        truths = {"atlas25": ((c2.d_atlas >= 0) & (c2.d_atlas <= 25)).values}
        if "d_kin_t5" in c2:
            truths["kinnex_t5_25"] = ((c2.d_kin_t5 >= 0) & (c2.d_kin_t5 <= 25)).values
        edges = np.arange(0, 1.0001, 0.1)
        for tname, tv in truths.items():
            b = np.clip(np.digitize(s[pool2], edges) - 1, 0, len(edges) - 2)
            sp, tp = s[pool2], tv[pool2]
            for k in range(len(edges) - 1):
                sel = b == k
                if sel.sum() == 0:
                    continue
                rows.append({"dataset": ds, "truth": tname,
                             "bin_lo": edges[k], "bin_hi": edges[k + 1],
                             "n": int(sel.sum()), "pred": float(sp[sel].mean()),
                             "obs": float(tp[sel].mean())})
    cal = pd.DataFrame(rows)
    cal["absgap"] = (cal.pred - cal.obs).abs()
    ece = (cal.groupby(["dataset", "truth"])
           .apply(lambda g: float((g.n * g.absgap).sum() / g.n.sum()), include_groups=False)
           .rename("ECE").reset_index())
    cal.to_csv(out / "calibration.tsv", sep="\t", index=False, float_format="%.6f")
    ece.to_csv(out / "calibration_ece.tsv", sep="\t", index=False, float_format="%.6f")
    print("\n== expected calibration error (tier-1 & IP-pass pool) ==")
    print(ece.to_string(index=False, float_format=lambda x: "%.4f" % x))
    return 0


if __name__ == "__main__":
    sys.exit(main())
