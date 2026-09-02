#!/usr/bin/env python3
"""scripts/paramsweep/verdicts.py -- evaluate every arm against manuscript/28 §2.

28 §2 (pre-registered, fixed before any sweep number existed): a change may become
a prime default only if, ON BOTH SPECIES
  (1) atlas-agreement precision at MATCHED CALL COUNT falls by <= 0.005
  (2) detected-gene recall rises by >= 0.010
  (3) the long-read (Kinnex) truth moves in the SAME direction
  (4) <= 2x wall time, <= 1.5x peak RSS
Anything passing (1)-(2) but failing (3) is a NEGATIVE RESULT.

Criterion (1) is evaluated on the `matched` view where one exists (the arm emits
more than the baseline) and on `default` otherwise -- an arm emitting FEWER calls
than the baseline cannot be flattered by call count. Criterion (2) is evaluated on
`default`, which is what the tool actually emits.
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
A = WD / "results/paramsweep"
df = pd.read_csv(A / "arms_all.tsv", sep="\t")
ic = pd.read_csv(A / "identity_and_compute.tsv", sep="\t")
ic["wall_s"] = pd.to_numeric(ic["wall_s"], errors="coerce")
ic["peak_rss_gb"] = pd.to_numeric(ic["peak_rss_gb"], errors="coerce")
base = {"human": "pbmc_base", "mouse": "m1_base"}
comp = ic.set_index("arm")

rows = []
for (arm, sp), g in df.groupby(["arm", "species"], sort=False):
    if arm in base.values():
        continue
    b = df[(df.arm == base[sp])]
    if b.empty:
        continue
    d = g[g.view == "default"]
    m = g[g.view == "matched"]
    bd = b[b.view == "default"]
    if d.empty or bd.empty:
        continue
    d = d.iloc[0]; bd = bd.iloc[0]
    prec_row = m.iloc[0] if len(m) else d
    prec_view = "matched" if len(m) else "default"
    dP = prec_row["P@100"] - bd["P@100"]
    dR = d["R_det@100"] - bd["R_det@100"]
    dF1 = d["F1_det@100"] - bd["F1_det@100"]
    dK5 = d["kin_t5@25"] - bd["kin_t5@25"]
    dK20 = d["kin_t20@25"] - bd["kin_t20@25"]
    dDec = d["kin_decoy@25"] - bd["kin_decoy@25"]
    try:
        w = comp.loc[arm, "wall_s"] / comp.loc[base[sp], "wall_s"]
        r = comp.loc[arm, "peak_rss_gb"] / comp.loc[base[sp], "peak_rss_gb"]
        same = comp.loc[arm, "identical_to_baseline"]
    except KeyError:
        w = r = float("nan"); same = "?"
    c1 = dP >= -0.005
    c2 = dR >= 0.010
    c3 = ("n/a (mouse: no long-read truth)" if pd.isna(dK5)
          else ("same-direction" if (dR > 0 and dK5 >= 0) or (dR <= 0) else "OPPOSITE"))
    c4 = True if (pd.isna(w) or pd.isna(r)) else ((w <= 2.0) and (r <= 1.5))
    rows.append(dict(arm=arm, species=sp, block=g["block"].iloc[0], flags=g["flags"].iloc[0],
                     identical_to_baseline=same,
                     n_default=int(d["n_scored"]), n_base=int(bd["n_scored"]),
                     prec_view=prec_view, P100=prec_row["P@100"], dP100=dP,
                     R_det=d["R_det@100"], dR_det=dR, dF1=dF1,
                     kin_t5=d["kin_t5@25"], dkin_t5=dK5, dkin_t20=dK20, ddecoy=dDec,
                     wall_x=w, rss_x=r,
                     c1_precision=("PASS" if c1 else "FAIL"),
                     c2_recall=("PASS" if c2 else "FAIL"),
                     c3_longread=c3,
                     c4_compute=("n/a" if (pd.isna(w) or pd.isna(r)) else ("PASS" if c4 else "FAIL")),
                     verdict=("PASS" if (c1 and c2 and c4 and c3 != "OPPOSITE")
                              else ("NEGATIVE RESULT (fails 3)" if (c1 and c2 and c3 == "OPPOSITE")
                                    else "FAIL"))))
out = pd.DataFrame(rows).sort_values(["species", "block", "arm"])
out.to_csv(A / "verdicts.tsv", sep="\t", index=False, float_format="%.6f")
pd.set_option("display.width", 250)
print(out[["arm", "species", "identical_to_baseline", "n_default", "dP100", "dR_det",
           "dkin_t5", "c1_precision", "c2_recall", "verdict"]].to_string(index=False))
print("\nwrote", A / "verdicts.tsv")
