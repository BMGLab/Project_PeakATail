#!/usr/bin/env python3
"""TASK D step 3 -- measure the scored arm against v2's default on every dataset.

The arms are reconstructed from the SAME no-internal-priming-filter runs the
transfer test used, because the reconstruction is exact: on GSE104556 mouse 1,
``tier == 1 & ip_tool_flag == 0 & molecules >= 2`` gives n 26,255 with
P@10 0.579661 / P@25 0.682042 / P@100 0.745001 / R_det@100 0.204790 /
F1 0.321268 -- digit-for-digit the published v2 default arm
(`manuscript/24_prime_preregistration.md` §2).  Reconstructing rather than
re-running means the two arms differ ONLY in the selection rule, with the same
peak calling, the same candidate set and the same veto underneath.

Reported, per `manuscript/24` §3.1:
  * PRIMARY   default vs default: dP@100, dR_det@100, dF1 against v2, with the
    pass/fail verdict spelled out per dataset.
  * SUPPORTING matched precision (|dP@100| <= 0.002) and matched call count
    (n +/- 1 %), plus the whole trade curve.
  * SECONDARY  P@10 / P@25 and their retention ratios; Kinnex t5/t20 P@25 and
    the internal-priming decoy rate (PBMC only); cross-mouse concordance.

Every arm it selects is also written out as a BED so `scripts/prime/score_pas.sh`
can re-derive the numbers with score_tool.py itself.
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
from taskD_eval import Evaluator, interp_at              # noqa: E402
from taskD_transfer import (FEATURESETS, TRANSFORM_OVERRIDE,  # noqa: E402
                            design, pool_of, prep)

from ema.countmatrix.pas_score import PasScoreModel      # noqa: E402


def write_bed(c, mask, path):
    d = c.loc[mask, ["chrom", "start", "end", "pas_id", "molecules", "strand"]]
    d = d.sort_values(["chrom", "start"], kind="stable")
    d.to_csv(path, sep="\t", header=False, index=False)
    return len(d)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--datasets", nargs="+", default=["mouse1", "pbmc", "mouse2"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    T, out = Path(a.tables), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(a.model).read_text())
    m = PasScoreModel(spec)
    thr = float(spec["threshold"])
    cols = FEATURESETS[spec["featureset"]]
    ov = TRANSFORM_OVERRIDE.get(spec["featureset"])
    print(f"[model] {spec['name']} family={spec['family']} featureset={spec['featureset']} "
          f"trained_on={spec['trained_on']} threshold={thr:.6f} ({spec['threshold_rule']})")

    rows, curves, verdict, scores = [], [], [], {}
    for ds in a.datasets:
        ev = Evaluator(T / ds)
        c = prep(ev.c)
        pool = pool_of(c)
        s = m.predict_proba_matrix(design(T / ds, cols, ov)[0])
        scores[ds] = (c, s, pool, ev)
        mol = c.molecules.values

        v2 = pool & (mol >= 2)
        pr = pool & (s >= thr)
        r_v2 = ev.row("v2_default", v2, {"dataset": ds})
        r_pr = ev.row("prime_scored", pr, {"dataset": ds})
        rows += [r_v2, r_pr]
        write_bed(c, v2, out / f"{ds}_v2_default.bed")
        write_bed(c, pr, out / f"{ds}_prime_scored.bed")

        # trade curve over the score, inside the same pool
        order = np.flatnonzero(pool)[np.argsort(-s[pool], kind="stable")]
        for n in np.unique(np.round(np.geomspace(2000, len(order), 80)).astype(int)):
            mk = np.zeros(len(c), bool)
            mk[order[:n]] = True
            r = ev.row("curve", mk, {"dataset": ds, "thresh": float(s[order[n - 1]])})
            curves.append(r)
        cu = pd.DataFrame([r for r in curves if r["dataset"] == ds])

        d = {"dataset": ds,
             "n_v2": r_v2["n"], "n_prime": r_pr["n"],
             "dP@100": r_pr["P@100"] - r_v2["P@100"],
             "dR_det@100": r_pr["R_det@100"] - r_v2["R_det@100"],
             "dF1": r_pr["F1_det@100"] - r_v2["F1_det@100"],
             "dP@10": r_pr["P@10"] - r_v2["P@10"],
             "dP@25": r_pr["P@25"] - r_v2["P@25"],
             "R_at_matched_P": interp_at(cu, "P@100", "R_det@100", r_v2["P@100"]),
             "n_at_matched_P": interp_at(cu, "P@100", "n", r_v2["P@100"]),
             "P_at_matched_n": interp_at(cu, "n", "P@100", r_v2["n"]),
             "R_at_matched_n": interp_at(cu, "n", "R_det@100", r_v2["n"])}
        d["rel_dR_at_matched_P"] = d["R_at_matched_P"] / r_v2["R_det@100"] - 1.0
        d["prereg_i_dP>=-0.005"] = bool(d["dP@100"] >= -0.005)
        d["prereg_ii_dR>=+0.010"] = bool(d["dR_det@100"] >= 0.010)
        d["prereg_iii_dF1>0"] = bool(d["dF1"] > 0)
        d["PASS"] = bool(d["prereg_i_dP>=-0.005"] and d["prereg_ii_dR>=+0.010"]
                         and d["prereg_iii_dF1>0"])
        verdict.append(d)

    pd.DataFrame(rows).to_csv(out / "arms.tsv", sep="\t", index=False, float_format="%.6f")
    pd.DataFrame(curves).to_csv(out / "trade_curves.tsv", sep="\t", index=False,
                                float_format="%.6f")
    V = pd.DataFrame(verdict)
    V.to_csv(out / "prereg_verdict.tsv", sep="\t", index=False, float_format="%.6f")

    # ---- cross-mouse concordance (24 §3.2.3) --------------------------------
    if {"mouse1", "mouse2"} <= set(a.datasets):
        conc = []
        for arm, mk in (("v2_default", lambda c, s, p, mol: p & (mol >= 2)),
                        ("prime_scored", lambda c, s, p, mol: p & (s >= thr))):
            pts = {}
            for ds in ("mouse1", "mouse2"):
                c, s, p, _ = scores[ds]
                sel = mk(c, s, p, c.molecules.values)
                pts[ds] = c.loc[sel, ["chrom", "start", "strand"]]
            for aa, bb in (("mouse1", "mouse2"), ("mouse2", "mouse1")):
                A, B = pts[aa], pts[bb]
                hit = 0
                for (ch, st), ga in A.groupby(["chrom", "strand"], sort=False):
                    gb_ = B[(B.chrom == ch) & (B.strand == st)].start.values
                    if len(gb_) == 0:
                        continue
                    gb_ = np.sort(gb_)
                    i = np.searchsorted(gb_, ga.start.values)
                    dl = np.where(i > 0, ga.start.values - gb_[np.maximum(i - 1, 0)], 10 ** 9)
                    dr = np.where(i < len(gb_), gb_[np.minimum(i, len(gb_) - 1)] - ga.start.values,
                                  10 ** 9)
                    hit += int((np.minimum(dl, dr) <= 25).sum())
                conc.append({"arm": arm, "from": aa, "to": bb, "n": len(A),
                             "frac_within_25bp": hit / max(len(A), 1)})
        pd.DataFrame(conc).to_csv(out / "cross_mouse_concordance.tsv", sep="\t",
                                  index=False, float_format="%.6f")
        print("\n== cross-mouse concordance (25 bp, same strand) ==")
        print(pd.DataFrame(conc).to_string(index=False, float_format=lambda x: "%.4f" % x))

    pd.set_option("display.width", 240)
    print("\n== ARMS ==")
    print(pd.DataFrame(rows)[["dataset", "arm", "n", "P@10", "P@25", "P@100",
                              "R_det@100", "F1_det@100", "KinP_t5@25",
                              "KinP_t20@25", "KinDecoy@25"]]
          .to_string(index=False, float_format=lambda x: "%.4f" % x))
    print("\n== manuscript/24 §3.1 VERDICT (default vs default) ==")
    print(V.to_string(index=False, float_format=lambda x: "%.4f" % x))
    print("\nOVERALL: %s" % ("PASS on all datasets" if V.PASS.all()
                             else "FAIL -- the flag stays defaulted OFF"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
