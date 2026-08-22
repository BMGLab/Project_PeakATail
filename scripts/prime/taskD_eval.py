#!/usr/bin/env python3
"""TASK D -- the evaluator, one dataset per instance.

Reproduces ``scripts/benchmark_tools/score_tool.py``'s metrics on arbitrary
SUBSETS of a candidate table without re-running bedtools per subset:

  P@W        fraction of SELECTED calls whose nearest same-strand PolyASite 2.0
             rep site is <= W bp  (score_tool: precision / real / atlas_full / W)
  R_det@W    fraction of the detected-gene-restricted atlas sites with at least
             one SELECTED call within W bp, same strand
             (score_tool: recall / real / atlas_detected / W)
  F1_det@W   harmonic mean of the two, on the detected flavour
  KinP_*@25  human only: fraction of SELECTED calls within 25 bp of a Kinnex
             long-read 3' end (t5 / t20) or of an internal-priming decoy

`score_pas.sh` remains the quoted authority; this module exists so a 60-point
threshold sweep costs seconds instead of an hour, and every arm it reports at a
fixed threshold is re-checked against `score_pas.sh` before it is believed.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


class Evaluator:
    def __init__(self, table_dir):
        d = Path(table_dir)
        self.meta = json.loads((d / "meta.json").read_text())
        self.c = pd.read_parquet(d / "candidates.parquet")
        self.idx = {p: i for i, p in enumerate(self.c.pas_id.values)}
        self.n_det = int(self.meta["n_det_atlas"])
        self.d_atlas = self.c.d_atlas.values
        self.pairs = {"det": self._pairs(d / "pairs_det_atlas.tsv")}
        p = d / "pairs_kin_t20.tsv"
        if p.exists():
            self.pairs["kin_t20"] = self._pairs(p)
        self.kin = {k[6:]: self.c[k].values for k in self.c.columns
                    if k.startswith("d_kin_")}

    def _pairs(self, path):
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["truth", "pas_id", "dist"], dtype={"pas_id": str})
        ci = df.pas_id.map(self.idx)
        keep = ci.notna().values
        return (df.truth.values[keep].astype(np.int64),
                ci.values[keep].astype(np.int64),
                df.dist.values[keep].astype(np.int64))

    # -- metrics -------------------------------------------------------------
    def precision(self, mask, W=100):
        d = self.d_atlas[mask]
        return float(np.sum((d >= 0) & (d <= W)) / len(d)) if len(d) else float("nan")

    def recall_det(self, mask, W=100):
        tr, ci, dd = self.pairs["det"]
        sel = mask[ci] & (dd <= W)
        return float(np.unique(tr[sel]).size / self.n_det)

    def kin_precision(self, mask, truth="t5", W=25):
        if truth not in self.kin:
            return float("nan")
        d = self.kin[truth][mask]
        return float(np.sum((d >= 0) & (d <= W)) / len(d)) if len(d) else float("nan")

    def kin_recall(self, mask, W=25):
        if "kin_t20" not in self.pairs:
            return float("nan")
        tr, ci, dd = self.pairs["kin_t20"]
        sel = mask[ci] & (dd <= W)
        return float(np.unique(tr[sel]).size / self.meta["n_kin_t20"])

    def row(self, label, mask, extra=None):
        p100 = self.precision(mask, 100)
        r100 = self.recall_det(mask, 100)
        r = {"arm": label, "n": int(mask.sum()),
             "P@10": self.precision(mask, 10),
             "P@25": self.precision(mask, 25),
             "P@100": p100,
             "R_det@100": r100,
             "F1_det@100": 0.0 if (p100 + r100) == 0 else 2 * p100 * r100 / (p100 + r100),
             "KinP_t5@25": self.kin_precision(mask, "t5", 25),
             "KinP_t20@25": self.kin_precision(mask, "t20", 25),
             "KinDecoy@25": self.kin_precision(mask, "decoy", 25),
             "KinR_t20@25": self.kin_recall(mask, 25)}
        if extra:
            r.update(extra)
        return r


def interp_at(curve: pd.DataFrame, xcol: str, ycol: str, x: float) -> float:
    """Linear interpolation of *ycol* at *xcol* == x along a monotone-ish sweep."""
    g = curve[[xcol, ycol]].dropna().sort_values(xcol)
    return float(np.interp(x, g[xcol].values, g[ycol].values))
