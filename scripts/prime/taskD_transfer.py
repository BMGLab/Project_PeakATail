#!/usr/bin/env python3
"""TASK D step 1 -- THE TRANSFER TEST.  Run and judged BEFORE any tool code.

`manuscript/23_algorithm_roadmap.md` §3 Step 1 and `manuscript/24_prime_
preregistration.md` §3.3 put this first on purpose: the result that decides
whether a per-site score is worth shipping is not how well a model fits, it is
whether a model fitted on one library still ranks on another.

PROTOCOL (fixed here, before the first number was produced)
-----------------------------------------------------------
* Recipe: the verifier's, unchanged --
  ``HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06,
  max_leaf_nodes=31, min_samples_leaf=100, l2_regularization=1.0,
  early_stopping=True, validation_fraction=0.1, random_state=0)``
  (`results/algo_headroom/VERIFY/code/v3_a3_decontaminate.py`).
* Features: only what the caller itself can compute at run time -- v2's
  ``pas_support.tsv`` columns plus the ``--pas-features on`` columns, plus
  ``width`` from the BED.  The **decontaminated** variant additionally drops
  every downstream-A feature, which is the variant the verifier found transfers
  BETTER because it stops the model reciting the Kinnex internal-priming rule.
  No annotation feature is used at all (A3's F_noann performs the same, and a
  GTF-derived feature would not survive a species change cleanly).
* Train on **GSE104556 testis mouse 1**, label = mouse PolyASite 2.0 within
  25 bp.  Mouse has no long-read truth, so the mouse side is atlas-trained; the
  mouse->PBMC direction is therefore a GENERALISATION test, never a headline.
* Evaluate the SAME fitted model, unchanged, on **PBMC 10k v3** and **testis
  mouse 2**.  Neither is ever used to fit anything.
* Deployment pool on every dataset: ``tier == 1 AND ip_tool_flag == 0`` -- the
  score is a re-ranker INSIDE the tier-1 and internal-priming gates and replaces
  ONLY the ">= 2 clip molecules" threshold.  Both gates stay hard.

THE BAR, DECLARED BEFORE THE RUN
--------------------------------
The incumbent is not merely "k >= 2": inside the pool, the natural incumbent
ranking is **by molecule count** (the verifier showed that ranking beats the
boundary statistic it withdrew).  The score must beat THAT, on data it never
saw.  TRANSFER PASSES iff on BOTH held-out datasets (PBMC and mouse 2):

  (a) at matched call count n = n(rule tier1 & IP-pass & >=2 molecules):
      dP@100 > 0 AND dR_det@100 > 0 against the rule, and
  (b) at matched precision P@100 = P@100(rule): relative recall gain >= +5 %,
      and the same two comparisons against the molecule-ranking control are
      not negative.

FAIL => stop; report the negative result; ship no model.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from taskD_eval import Evaluator, interp_at            # noqa: E402

# ---------------------------------------------------------------- features --
F_CLIP = ["clip_reads", "clip_umis", "clip_reads_f3844", "clip_umis_f3844",
          "tier1", "clip_positions", "clip_span"]
F_COV = ["window_reads", "width"]
F_HEX = ["hex_strong", "hex_any12", "hex_n_types", "hex_best_off",
         "hex_strong_off", "seq_ok"]
F_CTX = ["d_prev_cand", "d_next_cand", "n_cand_100", "n_cand_500",
         "mol_500_sum", "is_local_mol_max", "mol_frac_local"]
#: every downstream-A covariate the caller emits; dropped by the decontaminated
#: variant (the verifier's ADDOWN set, mapped onto this branch's column names).
F_ADOWN = ["ip_tool_flag", "ip_tool_afrac", "ip_tool_arun", "a_count_d18",
           "a_frac_d18", "a_run_d18", "a_frac_d30", "a_run_d30", "kin_ip_flag"]

FEATURESETS = {"NOA": F_CLIP + F_COV + F_HEX + F_CTX,
               "ALL": F_CLIP + F_COV + F_HEX + F_CTX + F_ADOWN}

LOGCOUNT = {"clip_reads", "clip_umis", "clip_reads_f3844", "clip_umis_f3844",
            "window_reads", "width", "clip_positions", "clip_span",
            "mol_500_sum", "n_cand_100", "n_cand_500", "hex_n_types",
            "a_count_d18", "a_run_d18", "a_run_d30", "ip_tool_arun"}
SIGNDIST = {"hex_best_off", "hex_strong_off"}
POSDIST = {"d_prev_cand", "d_next_cand"}


def prep(c: pd.DataFrame) -> pd.DataFrame:
    c = c.copy()
    c["tier1"] = (c.tier == 1).astype(np.int8)
    return c


def build_matrix(c: pd.DataFrame, cols) -> np.ndarray:
    X = np.empty((len(c), len(cols)), dtype=np.float64)
    for j, k in enumerate(cols):
        v = pd.to_numeric(c[k], errors="coerce").values.astype(np.float64)
        v = np.nan_to_num(v, nan=0.0)
        if k in LOGCOUNT:
            v = np.log1p(np.maximum(v, 0))
        elif k in SIGNDIST:
            v = np.sign(v) * np.log1p(np.abs(v))
        elif k in POSDIST:
            v = np.log1p(np.minimum(np.maximum(v, 0), 1e6))
        X[:, j] = v
    return X


def gb(seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(
        max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
        min_samples_leaf=100, l2_regularization=1.0,
        early_stopping=True, validation_fraction=0.1, random_state=seed)


def lr():
    """The SHIPPABLE alternative: ~20 coefficients evaluated as one dot product.

    A gradient-boosted ensemble can be shipped as numpy constants (arrays of
    node thresholds and children), but it is a ~1 MB opaque blob; a linear model
    is a table a reviewer can read.  It is fitted here so the choice between the
    two is made on measured transfer, not on taste.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    return Pipeline([("sc", StandardScaler()),
                     ("lr", LogisticRegression(C=1.0, max_iter=5000, solver="lbfgs"))])


MODELS = {"gb": gb, "lr": lr}


def pool_of(c: pd.DataFrame) -> np.ndarray:
    """tier-1 AND internal-priming-pass: the hard gates in front of the score."""
    ip = pd.to_numeric(c.get("ip_tool_flag"), errors="coerce").fillna(0).values
    return (c.tier.values == 1) & (ip == 0)


def sweep(ev, c, pool, score, label, n_points=60):
    order = np.flatnonzero(pool)[np.argsort(-score[pool], kind="stable")]
    rows = []
    for n in np.unique(np.round(np.geomspace(2000, len(order), n_points)).astype(int)):
        m = np.zeros(len(c), bool)
        m[order[:n]] = True
        rows.append(ev.row(label, m, {"kind": "sweep", "model": label,
                                      "thresh": float(score[order[n - 1]])}))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, help="dir holding <name>/candidates.parquet")
    ap.add_argument("--train", default="mouse1")
    ap.add_argument("--eval", nargs="+", default=["pbmc", "mouse2"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    T = Path(a.tables)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    names = [a.train] + [e for e in a.eval if e != a.train]

    ev, cand, pools = {}, {}, {}
    for nm in names:
        ev[nm] = Evaluator(T / nm)
        cand[nm] = prep(ev[nm].c)
        pools[nm] = pool_of(cand[nm])
        y = ((cand[nm].d_atlas >= 0) & (cand[nm].d_atlas <= 25)).values
        print(f"[{nm}] candidates {len(cand[nm]):,}  pool(tier1&IPpass) {pools[nm].sum():,}  "
              f"atlas+@25 {y.mean():.4f}")

    rows, models = [], {}
    ytr = ((cand[a.train].d_atlas >= 0) & (cand[a.train].d_atlas <= 25)).values.astype(int)
    for fs, cols in FEATURESETS.items():
        Xtr = build_matrix(cand[a.train], cols)
        for mk, mfn in MODELS.items():
            key = f"{mk}_{fs}"
            m = mfn()
            m.fit(Xtr, ytr)
            models[key] = (m, cols)
            print(f"[fit] {key}: {len(cols)} features on {a.train}, {ytr.sum():,} positives")
            for nm in names:
                s = m.predict_proba(build_matrix(cand[nm], cols))[:, 1]
                np.save(out / f"score_{key}_{nm}.npy", s)
                for r in sweep(ev[nm], cand[nm], pools[nm], s, f"MODEL_{key}"):
                    r["dataset"] = nm
                    r["is_train"] = nm == a.train
                    rows.append(r)

    # ---- incumbent arms on every dataset -----------------------------------
    for nm in names:
        c, p = cand[nm], pools[nm]
        mol = pd.to_numeric(c.molecules, errors="coerce").fillna(0).values
        for k in (1, 2, 3, 5, 10):
            r = ev[nm].row(f"RULE k>={k}", p & (mol >= k), {"kind": "rule", "k": k,
                                                            "model": "RULE"})
            r["dataset"] = nm
            r["is_train"] = nm == a.train
            rows.append(r)
        # the honest control: rank by molecule count (ties broken by clip reads)
        cr = pd.to_numeric(c.clip_reads, errors="coerce").fillna(0).values
        ctrl = mol + cr / (cr.max() + 1.0)
        for r in sweep(ev[nm], c, p, ctrl, "CTRL_molecules"):
            r["dataset"] = nm
            r["is_train"] = nm == a.train
            rows.append(r)

    d = pd.DataFrame(rows)
    d.to_csv(out / "transfer_curves.tsv", sep="\t", index=False, float_format="%.6f")

    # ---- the verdict table -------------------------------------------------
    res = []
    for nm in names:
        g = d[d.dataset == nm]
        base = g[(g.kind == "rule") & (g.k == 2)].iloc[0]
        for model in [f"MODEL_{k}" for k in models] + ["CTRL_molecules"]:
            cu = g[g.model == model].sort_values("n")
            if cu.empty:
                continue
            n0 = float(base["n"])
            res.append({
                "dataset": nm, "is_train": nm == a.train, "model": model,
                "rule_n": int(base["n"]), "rule_P@100": base["P@100"],
                "rule_R_det@100": base["R_det@100"],
                "P@100_at_matched_n": interp_at(cu, "n", "P@100", n0),
                "R_det@100_at_matched_n": interp_at(cu, "n", "R_det@100", n0),
                "P@10_at_matched_n": interp_at(cu, "n", "P@10", n0),
                "KinP_t5@25_at_matched_n": interp_at(cu, "n", "KinP_t5@25", n0),
                "KinDecoy@25_at_matched_n": interp_at(cu, "n", "KinDecoy@25", n0),
                "n_at_matched_P": interp_at(cu.sort_values("P@100", ascending=False),
                                            "P@100", "n", base["P@100"]),
                "R_det@100_at_matched_P": interp_at(cu.sort_values("P@100", ascending=False),
                                                    "P@100", "R_det@100", base["P@100"]),
            })
    v = pd.DataFrame(res)
    v["dR_matched_n"] = v["R_det@100_at_matched_n"] - v["rule_R_det@100"]
    v["dP_matched_n"] = v["P@100_at_matched_n"] - v["rule_P@100"]
    v["rel_dR_matched_P"] = v["R_det@100_at_matched_P"] / v["rule_R_det@100"] - 1.0
    v.to_csv(out / "transfer_verdict.tsv", sep="\t", index=False, float_format="%.6f")
    pd.set_option("display.width", 220)
    print("\n== TRANSFER VERDICT (train = %s) ==" % a.train)
    print(v[["dataset", "is_train", "model", "rule_n", "rule_P@100", "rule_R_det@100",
             "P@100_at_matched_n", "R_det@100_at_matched_n", "dP_matched_n",
             "dR_matched_n", "R_det@100_at_matched_P", "rel_dR_matched_P"]]
          .to_string(index=False, float_format=lambda x: "%.4f" % x))
    json.dump({"train": a.train, "eval": a.eval,
               "featuresets": {k: v_ for k, v_ in FEATURESETS.items()}},
              open(out / "protocol.json", "w"), indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
