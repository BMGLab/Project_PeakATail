#!/usr/bin/env python3
"""
replication_filter.py -- cross-sample replication filter for APA switches.

Stage-3 reliability program, requirement "replication" in
manuscript/13_reliability_positioning.md: a cell-type switch is reported only
if it is called, with the SAME direction, in >= K independent samples
(Laughney: 17 patients; testis: 2 mice).  This is the PI's strongest
false-positive control for cell-type-specific switches.

INPUTS (one per sample, repeat --sample)
    * an `ema switch diff` output directory (contains differential/), or the
      differential/ directory itself.  Per-pair TSVs are
      differential/<strategy>_<c1>_vs_<c2>.tsv with columns
          pas_id gene_id chrom start end strand cluster1 cluster2
          pvalue qvalue ... odds_ratio delta_proportion log2fc
      (ema/switch_test/runner.py::_augment_diff_df + strategies/fisher.py).
      `nb_multi` writes an omnibus table without pairs/effects and is NOT a
      valid input for a directional replication filter.
    * an `ema switch length` per-cell PDUI table (pdui_classic.tsv, or the
      directory holding it).  The tool summarises it to per-(feature, pair)
      dPDUI = mean PDUI(c1) - mean PDUI(c2), Mann-Whitney U p, BH q per pair.
    * a pre-summarised PDUI table with columns
      gene_id cluster1 cluster2 (dpdui|delta_pdui) qvalue [pvalue ...].

SIGN CONVENTIONS (confirmed in the PeakATail code, 2026-08-21)
    delta_proportion = prop(c1) - prop(c2)        (+ => higher in cluster1)
    log2fc           = log2(prop(c2) / prop(c1))  (+ => higher in cluster2)
    nb_pairwise log2fc = cluster2-indicator coefficient (+ => higher in c2)
    dpdui            = mean PDUI(c1) - mean PDUI(c2) (+ => cluster1 longer)
  Direction consistency is evaluated on ONE effect column across all samples,
  so the convention only matters for the human-readable `utr_direction`
  column and the (c1,c2) orientation canonicalisation: when a sample wrote
  the pair as <c2>_vs_<c1>, the effect sign is flipped before comparison.

OUTPUTS (--out DIR)
    all_features.tsv      every (pair, feature) seen in >= 1 sample, with
                          per-sample q__<s>, effect__<s>, called__<s> columns,
                          replication counts and flags
    replicated.tsv        rows with passes_replication == True
    null_control.tsv      one row per null combination (if nulls given)
    summary.json          machine-readable report (params, counts, null)

Run the unit tests:
    .venv/bin/python -m pytest scripts/reliability/tests -v
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
import math
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

log = logging.getLogger("replication_filter")

__version__ = "0.1.0"

KNOWN_STRATEGIES = ("nb_pairwise", "nb_multi", "fisher")
# effect column -> (+1 if positive means "higher in cluster1", -1 if cluster2)
EFFECT_SIGN = {
    "delta_proportion": +1,
    "dpdui": +1,
    "delta_pdui": +1,
    "log2fc": -1,
}
EFFECT_PREFERENCE = ("delta_proportion", "dpdui", "delta_pdui", "log2fc")
META_COLS = ["feature_id", "pas_id", "gene_id", "chrom", "start", "end", "strand"]
SKIP_FILES = ("_omnibus.tsv", "switch_diff_long")


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def _finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe, monotone)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def parse_pair_filename(name: str, strategy: str | None) -> tuple[str, str, str] | None:
    """Return (strategy, c1, c2) from `<strategy>_<c1>_vs_<c2>.tsv`, or None.

    Cluster labels and strategy names may both contain underscores, so the
    strategy prefix is stripped first (explicit --strategy, else the longest
    known strategy that matches, else the first token) and the remainder is
    split on the single `_vs_`.
    """
    if not name.endswith(".tsv") or "_vs_" not in name:
        return None
    stem = name[:-4]
    strat = None
    if strategy and stem.startswith(strategy + "_"):
        strat = strategy
    elif strategy:
        return None
    else:
        for s in sorted(KNOWN_STRATEGIES, key=len, reverse=True):
            if stem.startswith(s + "_"):
                strat = s
                break
        if strat is None:
            strat = stem.split("_", 1)[0]
    rest = stem[len(strat) + 1:]
    parts = rest.split("_vs_")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        return None
    return strat, parts[0], parts[1]


def canonical_pair(c1: str, c2: str) -> tuple[str, str, bool]:
    """Lexicographic canonical orientation; flipped=True when (c1,c2) swapped."""
    a, b = str(c1), str(c2)
    if a <= b:
        return a, b, False
    return b, a, True


def pick_effect_col(columns: Iterable[str], requested: str = "auto") -> str:
    cols = list(columns)
    if requested != "auto":
        if requested not in cols:
            raise KeyError(f"effect column {requested!r} not in table columns {cols}")
        return requested
    for c in EFFECT_PREFERENCE:
        if c in cols:
            return c
    raise KeyError(
        "no signed effect column found (looked for "
        f"{EFFECT_PREFERENCE}); columns were {cols}"
    )


# --------------------------------------------------------------------------
# loading: switch diff output dirs
# --------------------------------------------------------------------------
def resolve_differential_dir(path: str | Path) -> Path:
    p = Path(path)
    if p.is_file():
        return p.parent
    if (p / "differential").is_dir():
        return p / "differential"
    if p.is_dir():
        return p
    raise FileNotFoundError(f"{path}: not a switch-diff output dir")


def list_pair_files(diff_dir: Path, strategy: str | None) -> list[tuple[Path, str, str, str]]:
    out = []
    for f in sorted(diff_dir.glob("*.tsv")):
        if any(s in f.name for s in SKIP_FILES):
            continue
        parsed = parse_pair_filename(f.name, strategy)
        if parsed is None:
            continue
        out.append((f, *parsed))
    return out


def load_diff_sample(
    sample: str,
    path: str | Path,
    strategy: str | None = None,
    effect_col: str = "auto",
    pairs: set[tuple[str, str]] | None = None,
) -> tuple[pd.DataFrame, str, str]:
    """Load one sample's per-pair switch-diff TSVs into the long schema.

    Returns (long_df, strategy_used, effect_col_used).  long_df columns:
    sample, c1, c2 (canonical), flipped, feature_id, pas_id, gene_id, chrom,
    start, end, strand, pvalue, qvalue, effect.
    """
    diff_dir = resolve_differential_dir(path)
    files = list_pair_files(diff_dir, strategy)
    if not files:
        raise FileNotFoundError(
            f"{sample}: no <strategy>_<c1>_vs_<c2>.tsv files in {diff_dir}"
            + (f" for strategy {strategy!r}" if strategy else "")
        )
    strategies = {s for _, s, _, _ in files}
    if len(strategies) > 1:
        raise ValueError(
            f"{sample}: several strategies in {diff_dir}: {sorted(strategies)}; "
            "pass --strategy"
        )
    strat = strategies.pop()
    frames = []
    eff_used = None
    for f, _, fc1, fc2 in files:
        df = pd.read_csv(f, sep="\t", dtype={"pas_id": str, "gene_id": str,
                                             "chrom": str, "cluster1": str,
                                             "cluster2": str})
        if df.empty:
            continue
        # self-describing cluster columns win over the filename parse
        if "cluster1" in df.columns and "cluster2" in df.columns:
            c1s = df["cluster1"].astype(str).unique()
            c2s = df["cluster2"].astype(str).unique()
            if len(c1s) != 1 or len(c2s) != 1:
                raise ValueError(f"{f}: mixed cluster1/cluster2 labels in one file")
            fc1, fc2 = str(c1s[0]), str(c2s[0])
        a, b, flipped = canonical_pair(fc1, fc2)
        if pairs is not None and (a, b) not in pairs:
            continue
        eff = pick_effect_col(df.columns, effect_col)
        if eff_used is None:
            eff_used = eff
        elif eff != eff_used:
            raise ValueError(f"{sample}: effect column changed between files ({eff_used} vs {eff})")
        if "qvalue" not in df.columns:
            raise ValueError(f"{f}: no qvalue column")
        out = pd.DataFrame({
            "sample": sample,
            "c1": a, "c2": b, "flipped": flipped,
            "pas_id": df["pas_id"].astype(str) if "pas_id" in df.columns else "",
            "gene_id": df["gene_id"].fillna("").astype(str) if "gene_id" in df.columns else "",
            "chrom": df["chrom"].fillna("").astype(str) if "chrom" in df.columns else "",
            "start": pd.to_numeric(df["start"], errors="coerce") if "start" in df.columns else np.nan,
            "end": pd.to_numeric(df["end"], errors="coerce") if "end" in df.columns else np.nan,
            "strand": df["strand"].fillna("").astype(str) if "strand" in df.columns else "",
            "pvalue": pd.to_numeric(df["pvalue"], errors="coerce") if "pvalue" in df.columns else np.nan,
            "qvalue": pd.to_numeric(df["qvalue"], errors="coerce"),
            "effect": pd.to_numeric(df[eff], errors="coerce") * (-1.0 if flipped else 1.0),
        })
        out["feature_id"] = out["pas_id"]
        frames.append(out)
    if not frames:
        return _empty_long(), strat, eff_used or "delta_proportion"
    long_df = pd.concat(frames, ignore_index=True)
    # duplicate feature within one (sample, pair): keep the smallest q
    long_df = (long_df.sort_values("qvalue", na_position="last")
               .drop_duplicates(["sample", "c1", "c2", "feature_id"], keep="first"))
    return long_df.reset_index(drop=True), strat, eff_used


def _empty_long() -> pd.DataFrame:
    return pd.DataFrame(columns=["sample", "c1", "c2", "flipped", "pas_id", "gene_id",
                                 "chrom", "start", "end", "strand", "pvalue", "qvalue",
                                 "effect", "feature_id"])


# --------------------------------------------------------------------------
# loading: PDUI tables
# --------------------------------------------------------------------------
def resolve_pdui_file(path: str | Path) -> Path:
    p = Path(path)
    if p.is_file():
        return p
    for cand in ("pdui_classic.tsv",):
        if (p / cand).is_file():
            return p / cand
    hits = sorted(p.glob("**/pdui_classic.tsv")) if p.is_dir() else []
    if len(hits) == 1:
        return hits[0]
    raise FileNotFoundError(f"{path}: no pdui_classic.tsv found (or ambiguous: {hits})")


def shuffle_cell_labels(df: pd.DataFrame, cluster_col: str, rng: np.random.Generator,
                        cell_col: str = "cell") -> pd.DataFrame:
    """Label-shuffle null at the CELL level for a long gene x cell table.

    One permutation of the cell -> cluster map is drawn (label multiset over
    cells preserved) and broadcast to every row of that cell, so a cell keeps a
    single label across all genes -- exactly what `ema switch diff`'s null
    driver does on ``obs[label]``.  Permuting the column row-wise instead would
    give each gene an independent, artificially balanced label split.
    """
    if cell_col not in df.columns:
        raise ValueError(f"label shuffle needs a {cell_col!r} column")
    cells = df.drop_duplicates(cell_col)[[cell_col, cluster_col]]
    perm = rng.permutation(cells[cluster_col].to_numpy())
    mapping = dict(zip(cells[cell_col].to_numpy(), perm))
    out = df.copy()
    out[cluster_col] = out[cell_col].map(mapping).to_numpy()
    return out


def summarise_pdui_cells(
    cells: pd.DataFrame,
    pairs: set[tuple[str, str]] | None = None,
    min_cells: int = 10,
    cluster_col: str = "cluster",
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Per-cell PDUI table -> per-(feature, canonical pair) dPDUI / p / q.

    feature_id = gene_id, or gene_id:transcript_id for per-isoform tables.
    dpdui = mean PDUI(c1) - mean PDUI(c2) (+ => cluster1 longer).  p from a
    two-sided Mann-Whitney U on the per-cell PDUI values; BH within each pair
    (mirrors switch diff, which BH-corrects per pair).  If `rng` is given the
    cluster labels are permuted across cells first (label-shuffle null).
    """
    from scipy.stats import mannwhitneyu

    need = {"gene_id", "pdui", cluster_col}
    missing = need - set(cells.columns)
    if missing:
        raise ValueError(f"PDUI table lacks columns {sorted(missing)}")
    df = cells.copy()
    df["pdui"] = pd.to_numeric(df["pdui"], errors="coerce")
    df = df[np.isfinite(df["pdui"])]
    df[cluster_col] = df[cluster_col].astype(str)
    if rng is not None:
        df = shuffle_cell_labels(df, cluster_col, rng)
    tid = df["transcript_id"].astype(str) if "transcript_id" in df.columns else pd.Series("_gene_", index=df.index)
    df["feature_id"] = np.where(tid.isin(["_gene_", "", "nan"]), df["gene_id"].astype(str),
                                df["gene_id"].astype(str) + ":" + tid)
    clusters = sorted(df[cluster_col].unique())
    all_pairs = [(a, b) for i, a in enumerate(clusters) for b in clusters[i + 1:]]
    if pairs is not None:
        all_pairs = [p for p in all_pairs if p in pairs]
    coord_cols = [c for c in ("distal_chrom", "distal_start", "distal_end", "distal_strand",
                              "proximal_pas_id", "distal_pas_id") if c in df.columns]
    meta = df.drop_duplicates("feature_id").set_index("feature_id")
    rows = []
    grouped = {k: g for k, g in df.groupby(["feature_id", cluster_col])["pdui"]}
    features = df["feature_id"].unique()
    for a, b in all_pairs:
        for fid in features:
            x = grouped.get((fid, a))
            y = grouped.get((fid, b))
            if x is None or y is None or len(x) < min_cells or len(y) < min_cells:
                continue
            xv, yv = x.to_numpy(), y.to_numpy()
            d = float(xv.mean() - yv.mean())
            if np.all(xv == xv[0]) and np.all(yv == yv[0]) and xv[0] == yv[0]:
                p = 1.0
            else:
                p = float(mannwhitneyu(xv, yv, alternative="two-sided").pvalue)
            m = meta.loc[fid]
            rows.append({
                "c1": a, "c2": b, "feature_id": fid,
                "gene_id": str(m["gene_id"]),
                "pas_id": str(m["distal_pas_id"]) if "distal_pas_id" in coord_cols else "",
                "chrom": str(m["distal_chrom"]) if "distal_chrom" in coord_cols else "",
                "start": pd.to_numeric(m.get("distal_start", np.nan), errors="coerce"),
                "end": pd.to_numeric(m.get("distal_end", np.nan), errors="coerce"),
                "strand": str(m["distal_strand"]) if "distal_strand" in coord_cols else "",
                "n_cells_c1": int(len(xv)), "n_cells_c2": int(len(yv)),
                "mean_pdui_c1": float(xv.mean()), "mean_pdui_c2": float(yv.mean()),
                "pvalue": p, "dpdui": d,
            })
    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=["c1", "c2", "feature_id", "gene_id", "pas_id", "chrom",
                                    "start", "end", "strand", "n_cells_c1", "n_cells_c2",
                                    "mean_pdui_c1", "mean_pdui_c2", "pvalue", "dpdui", "qvalue"])
        return out
    out["qvalue"] = np.nan
    for (a, b), idx in out.groupby(["c1", "c2"]).groups.items():
        out.loc[idx, "qvalue"] = bh_qvalues(out.loc[idx, "pvalue"].to_numpy())
    return out


def load_pdui_sample(
    sample: str,
    path: str | Path,
    pairs: set[tuple[str, str]] | None = None,
    min_cells: int = 10,
    effect_col: str = "auto",
    rng: np.random.Generator | None = None,
    _cache: dict | None = None,
) -> tuple[pd.DataFrame, str]:
    """Load a per-cell or pre-summarised PDUI table into the long schema."""
    f = resolve_pdui_file(path)
    if _cache is not None and f in _cache:
        raw = _cache[f]
    else:
        raw = pd.read_csv(f, sep="\t", dtype={"gene_id": str, "transcript_id": str,
                                              "cluster": str, "cluster1": str,
                                              "cluster2": str})
        if _cache is not None:
            _cache[f] = raw
    if {"cell", "pdui"} <= set(raw.columns) and "cluster1" not in raw.columns:
        summ = summarise_pdui_cells(raw, pairs=pairs, min_cells=min_cells, rng=rng)
        eff = "dpdui"
        summ["flipped"] = False
    else:
        if rng is not None:
            raise ValueError(f"{sample}: internal label shuffles need a per-cell PDUI table")
        eff = pick_effect_col(raw.columns, effect_col)
        need = {"gene_id", "cluster1", "cluster2", "qvalue"}
        if not need <= set(raw.columns):
            raise ValueError(f"{f}: pre-summarised PDUI table needs {sorted(need)}")
        summ = raw.copy()
        canon = summ.apply(lambda r: canonical_pair(r["cluster1"], r["cluster2"]), axis=1,
                           result_type="expand")
        summ["c1"], summ["c2"], summ["flipped"] = canon[0], canon[1], canon[2]
        if pairs is not None:
            summ = summ[[(a, b) in pairs for a, b in zip(summ["c1"], summ["c2"])]]
        summ["feature_id"] = summ["feature_id"] if "feature_id" in summ.columns else summ["gene_id"].astype(str)
        summ[eff] = pd.to_numeric(summ[eff], errors="coerce") * np.where(summ["flipped"], -1.0, 1.0)
    out = pd.DataFrame({
        "sample": sample,
        "c1": summ["c1"].astype(str), "c2": summ["c2"].astype(str), "flipped": summ["flipped"],
        "pas_id": summ["pas_id"].astype(str) if "pas_id" in summ.columns else "",
        "gene_id": summ["gene_id"].astype(str),
        "chrom": summ["chrom"].astype(str) if "chrom" in summ.columns else "",
        "start": pd.to_numeric(summ["start"], errors="coerce") if "start" in summ.columns else np.nan,
        "end": pd.to_numeric(summ["end"], errors="coerce") if "end" in summ.columns else np.nan,
        "strand": summ["strand"].astype(str) if "strand" in summ.columns else "",
        "pvalue": pd.to_numeric(summ["pvalue"], errors="coerce") if "pvalue" in summ.columns else np.nan,
        "qvalue": pd.to_numeric(summ["qvalue"], errors="coerce"),
        "effect": pd.to_numeric(summ[eff], errors="coerce"),
        "feature_id": summ["feature_id"].astype(str),
    })
    out = (out.sort_values("qvalue", na_position="last")
           .drop_duplicates(["sample", "c1", "c2", "feature_id"], keep="first"))
    return out.reset_index(drop=True), eff


def detect_input_kind(path: str | Path) -> str:
    p = Path(path)
    if p.is_file():
        head = pd.read_csv(p, sep="\t", nrows=0).columns
        if "pdui" in head or "dpdui" in head or "delta_pdui" in head:
            return "pdui"
        return "diff"
    if (p / "differential").is_dir():
        return "diff"
    if (p / "pdui_classic.tsv").is_file():
        return "pdui"
    if p.is_dir() and any(parse_pair_filename(f.name, None) for f in p.glob("*.tsv")):
        return "diff"
    raise FileNotFoundError(f"{path}: cannot tell whether this is a switch-diff dir or a PDUI table")


# --------------------------------------------------------------------------
# the replication computation
# --------------------------------------------------------------------------
@dataclass
class Params:
    min_samples: int = 2
    fdr: float = 0.05
    effect_floor: float = 0.1
    discordant_policy: str = "exclude"   # or "allow"
    level: str = "pas"                   # or "gene"
    effect_sign: int = +1                # +1: + means higher in c1


def replication_table(long_df: pd.DataFrame, prm: Params, samples: list[str]) -> pd.DataFrame:
    """Collapse the long per-sample table to one row per (pair, feature).

    Columns added:
      n_samples_tested, n_samples_called, n_pos, n_neg, replication_count,
      consensus_direction, direction_consistent, passes_replication,
      n_pos_floor, n_neg_floor, replication_count_floor,
      direction_consistent_floor, effect_floor_pass, passes_replication_floor,
      min_q_called, mean_effect_consensus, utr_direction,
      q__<s>, effect__<s>, called__<s> for every sample.
    """
    if long_df.empty:
        return pd.DataFrame()
    df = long_df.copy()
    df["called"] = (df["qvalue"] < prm.fdr) & np.isfinite(df["effect"])
    df["pos"] = df["called"] & (df["effect"] > 0)
    df["neg"] = df["called"] & (df["effect"] < 0)
    df["floor_ok"] = np.isfinite(df["effect"]) & (df["effect"].abs() >= prm.effect_floor)
    df["pos_f"] = df["pos"] & df["floor_ok"]
    df["neg_f"] = df["neg"] & df["floor_ok"]
    df["q_called"] = df["qvalue"].where(df["called"])

    keys = ["c1", "c2", "feature_id"]
    agg = df.groupby(keys, sort=True).agg(
        n_samples_tested=("sample", "nunique"),
        n_samples_called=("called", "sum"),
        n_pos=("pos", "sum"), n_neg=("neg", "sum"),
        n_pos_floor=("pos_f", "sum"), n_neg_floor=("neg_f", "sum"),
        min_q_called=("q_called", "min"),
    ).reset_index()
    for c in ("n_samples_called", "n_pos", "n_neg", "n_pos_floor", "n_neg_floor"):
        agg[c] = agg[c].astype(int)

    # first non-empty metadata per feature (coordinates may be absent in some samples)
    meta = (df.sort_values(["qvalue"], na_position="last")
            .groupby(keys, sort=True)[["pas_id", "gene_id", "chrom", "start", "end", "strand"]]
            .first().reset_index())
    agg = agg.merge(meta, on=keys, how="left")

    agg["replication_count"] = agg[["n_pos", "n_neg"]].max(axis=1)
    agg["consensus_direction"] = np.select(
        [agg["n_pos"] > agg["n_neg"], agg["n_neg"] > agg["n_pos"],
         (agg["n_pos"] == agg["n_neg"]) & (agg["n_pos"] > 0)],
        ["+", "-", "tie"], default="")
    agg["direction_consistent"] = (agg["replication_count"] > 0) & (agg[["n_pos", "n_neg"]].min(axis=1) == 0)
    agg["replication_count_floor"] = agg[["n_pos_floor", "n_neg_floor"]].max(axis=1)
    agg["direction_consistent_floor"] = (agg["replication_count_floor"] > 0) & (
        agg[["n_pos_floor", "n_neg_floor"]].min(axis=1) == 0)

    passes = agg["replication_count"] >= prm.min_samples
    passes_f = agg["replication_count_floor"] >= prm.min_samples
    if prm.discordant_policy == "exclude":
        passes &= agg["direction_consistent"]
        passes_f &= agg["direction_consistent_floor"]
    agg["passes_replication"] = passes
    agg["effect_floor_pass"] = passes_f      # the effect-size-floor flag
    agg["passes_replication_floor"] = passes & passes_f

    # mean effect over the called samples in the consensus direction
    cons = agg.set_index(keys)["consensus_direction"]
    d2 = df.merge(cons.rename("cons").reset_index(), on=keys, how="left")
    in_cons = d2["called"] & (((d2["cons"] == "+") & (d2["effect"] > 0)) |
                              ((d2["cons"] == "-") & (d2["effect"] < 0)))
    me = d2[in_cons].groupby(keys)["effect"].mean().rename("mean_effect_consensus").reset_index()
    agg = agg.merge(me, on=keys, how="left")

    # human-readable UTR direction for the consensus sign (effect_sign: +1
    # means '+' == higher in cluster1).  Only meaningful at the distal PAS /
    # for dPDUI; reported as a convenience, never used for filtering.
    up = np.where(prm.effect_sign >= 0, "higher_in_c1", "higher_in_c2")
    down = np.where(prm.effect_sign >= 0, "higher_in_c2", "higher_in_c1")
    agg["utr_direction"] = np.select(
        [agg["consensus_direction"] == "+", agg["consensus_direction"] == "-"],
        [str(up), str(down)], default="")

    # per-sample wide columns.  (sample, pair, feature) is unique after the
    # loaders' drop_duplicates, so a plain unstack is exact and keeps dtypes
    # (pivot_table silently drops / coerces boolean value columns).
    idx = pd.MultiIndex.from_frame(agg[keys])
    stacked = df.set_index(keys + ["sample"])
    wide_q = stacked["qvalue"].unstack("sample")
    wide_e = stacked["effect"].unstack("sample")
    wide_c = stacked["called"].unstack("sample")
    for s in samples:
        agg[f"q__{s}"] = wide_q[s].reindex(idx).to_numpy() if s in wide_q.columns else np.nan
        agg[f"effect__{s}"] = wide_e[s].reindex(idx).to_numpy() if s in wide_e.columns else np.nan
        if s in wide_c.columns:
            # nullable boolean: True/False where tested, <NA> where the sample lacks the row
            agg[f"called__{s}"] = pd.array(wide_c[s].reindex(idx).to_numpy(), dtype="boolean")
        else:
            agg[f"called__{s}"] = pd.array([pd.NA] * len(agg), dtype="boolean")
    front = ["c1", "c2", "feature_id", "pas_id", "gene_id", "chrom", "start", "end", "strand",
             "n_samples_tested", "n_samples_called", "n_pos", "n_neg", "replication_count",
             "consensus_direction", "direction_consistent", "passes_replication",
             "n_pos_floor", "n_neg_floor", "replication_count_floor", "direction_consistent_floor",
             "effect_floor_pass", "passes_replication_floor", "min_q_called",
             "mean_effect_consensus", "utr_direction"]
    rest = [c for c in agg.columns if c not in front]
    return agg[front + rest].sort_values(["c1", "c2", "replication_count", "min_q_called"],
                                         ascending=[True, True, False, True]).reset_index(drop=True)


def collapse_to_gene(pas_table: pd.DataFrame) -> pd.DataFrame:
    """Gene-level view: per (pair, gene) keep the best-replicating PAS.

    A gene replicates iff one of its PAS replicates (same PAS, same direction,
    >= K samples).  Adds representative_pas_id, n_pas_tested, n_pas_replicated
    and `distal_pas_replicated` / `distal_utr_direction` from the strand-aware
    3'-most PAS of the gene (needs coordinates; '' otherwise).
    """
    if pas_table.empty:
        return pas_table
    t = pas_table.copy()
    t = t[t["gene_id"].astype(str) != ""]
    if t.empty:
        return t
    t["_rank"] = t["replication_count"].astype(float)
    t = t.sort_values(["c1", "c2", "gene_id", "passes_replication", "_rank", "min_q_called"],
                      ascending=[True, True, True, False, False, True])
    keys = ["c1", "c2", "gene_id"]
    best = t.drop_duplicates(keys, keep="first").copy()
    best = best.rename(columns={"pas_id": "representative_pas_id"})
    stats = t.groupby(keys).agg(n_pas_tested=("feature_id", "nunique"),
                                n_pas_replicated=("passes_replication", "sum")).reset_index()
    best = best.merge(stats, on=keys, how="left")
    best["feature_id"] = best["gene_id"]
    # distal PAS per gene (3'-most by strand), over all pairs
    coords = t[(t["strand"].isin(["+", "-"])) & np.isfinite(pd.to_numeric(t["start"], errors="coerce"))]
    distal_dir, distal_rep = {}, {}
    if not coords.empty:
        for (a, b, g), grp in coords.groupby(keys):
            strand = grp["strand"].iloc[0]
            if strand == "+":
                row = grp.loc[pd.to_numeric(grp["end"]).idxmax()]
            else:
                row = grp.loc[pd.to_numeric(grp["start"]).idxmin()]
            distal_rep[(a, b, g)] = bool(row["passes_replication"])
            cd = row["consensus_direction"]
            if not row["passes_replication"] or cd not in ("+", "-"):
                distal_dir[(a, b, g)] = ""
            else:
                # '+' on the distal PAS == higher distal usage in the cluster
                # the effect convention calls 'higher' => that cluster is longer
                distal_dir[(a, b, g)] = {"higher_in_c1": "c1_longer",
                                         "higher_in_c2": "c2_longer"}[row["utr_direction"]]
    idx = list(zip(best["c1"], best["c2"], best["gene_id"]))
    best["distal_pas_replicated"] = [distal_rep.get(k, False) for k in idx]
    best["distal_utr_direction"] = [distal_dir.get(k, "") for k in idx]
    best = best.drop(columns=["_rank"])
    front = ["c1", "c2", "feature_id", "gene_id", "representative_pas_id", "n_pas_tested",
             "n_pas_replicated", "distal_pas_replicated", "distal_utr_direction"]
    rest = [c for c in best.columns if c not in front]
    return best[front + rest].reset_index(drop=True)


def count_replicated(table: pd.DataFrame) -> dict:
    if table.empty:
        return {"n_features": 0, "n_replicated": 0, "n_replicated_floor": 0, "per_pair": {}}
    per_pair = {}
    for (a, b), g in table.groupby(["c1", "c2"]):
        per_pair[f"{a}_vs_{b}"] = {
            "n_features": int(len(g)),
            "n_replicated": int(g["passes_replication"].sum()),
            "n_replicated_floor": int(g["passes_replication_floor"].sum()),
        }
    return {
        "n_features": int(len(table)),
        "n_replicated": int(table["passes_replication"].sum()),
        "n_replicated_floor": int(table["passes_replication_floor"].sum()),
        "per_pair": per_pair,
    }


# --------------------------------------------------------------------------
# null control
# --------------------------------------------------------------------------
def expand_null_specs(specs: list[str], samples: list[str]) -> dict[str, list[Path]]:
    """--null-dirs specs -> {sample: [perm_dir, ...]} (sorted, de-duplicated).

    Accepted forms (repeatable):
      NAME=GLOB           null dirs of sample NAME
      GLOB with {sample}  the token is substituted with each sample name
      GLOB (bare)         positional: the i-th bare glob belongs to the i-th sample
    """
    out: dict[str, list[Path]] = {s: [] for s in samples}
    positional = 0
    for spec in specs:
        if "=" in spec and not spec.startswith("="):
            name, pat = spec.split("=", 1)
            if name not in out:
                raise KeyError(f"--null-dirs {spec!r}: unknown sample {name!r}; samples are {samples}")
            out[name].extend(sorted(Path(p) for p in glob.glob(pat)))
        elif "{sample}" in spec:
            for s in samples:
                out[s].extend(sorted(Path(p) for p in glob.glob(spec.replace("{sample}", s))))
        else:
            if positional >= len(samples):
                raise ValueError(f"more bare --null-dirs globs than samples ({len(samples)})")
            out[samples[positional]].extend(sorted(Path(p) for p in glob.glob(spec)))
            positional += 1
    for s in samples:
        out[s] = sorted(dict.fromkeys(out[s]))
    return out


def null_combinations(n_per_sample: dict[str, int], n_combos: int | None,
                      rng: np.random.Generator) -> list[dict[str, int]]:
    """Which null perm index to use for each sample in each combination.

    Default (n_combos None): index-aligned, min(n_perms) combos (perm i of
    every sample).  Otherwise n_combos independent random draws per sample.
    Samples with 0 nulls are dropped from the null (reported in summary).
    """
    avail = {s: n for s, n in n_per_sample.items() if n > 0}
    if not avail:
        return []
    if n_combos is None:
        m = min(avail.values())
        return [{s: i for s in avail} for i in range(m)]
    return [{s: int(rng.integers(0, n)) for s, n in avail.items()} for _ in range(n_combos)]


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_sample_spec(spec: str) -> tuple[str, str]:
    if "=" in spec and not spec.startswith("=") and not os.path.exists(spec):
        name, path = spec.split("=", 1)
        return name, path
    p = Path(spec)
    name = p.stem if p.is_file() else p.name
    if name == "differential":
        name = p.parent.name
    return name, spec


def parse_pairs(spec: str | None) -> set[tuple[str, str]] | None:
    if not spec:
        return None
    out = set()
    for item in spec.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = [x.strip() for x in item.split(",")]
        if len(parts) != 2:
            raise ValueError(f"--pairs item {item!r} is not 'c1,c2'")
        a, b, _ = canonical_pair(*parts)
        out.add((a, b))
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="replication_filter.py",
        description=__doc__.split("\n\n")[0] + "\n\n" + "See module docstring / README for the full contract.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--sample", action="append", required=True, metavar="NAME=PATH|PATH",
                    help="per-sample switch-diff output dir (or differential/ dir), "
                         "or PDUI table/dir. Repeat once per sample.")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--input-kind", choices=["auto", "diff", "pdui"], default="auto")
    ap.add_argument("--strategy", default=None,
                    help="TSV prefix (fisher, nb_pairwise). Auto-detected when unique.")
    ap.add_argument("--effect-col", default="auto",
                    help="signed effect column (auto: delta_proportion > dpdui > log2fc)")
    ap.add_argument("--level", choices=["pas", "gene"], default="pas",
                    help="report PAS-level rows or collapse to genes (PDUI input is gene/"
                         "transcript-level already)")
    ap.add_argument("--min-samples", "-k", type=int, default=2,
                    help="K: minimum number of samples called in the same direction [2]")
    ap.add_argument("--fdr", type=float, default=0.05, help="per-sample q threshold [0.05]")
    ap.add_argument("--effect-floor", type=float, default=0.1,
                    help="|effect| floor for the effect_floor_pass flag [0.1]")
    ap.add_argument("--discordant-policy", choices=["exclude", "allow"], default="exclude",
                    help="exclude: a feature called in the OPPOSITE direction in any sample "
                         "never passes [exclude]")
    ap.add_argument("--pairs", default=None, help="restrict to 'c1,c2;c3,c4' (orientation-free)")
    ap.add_argument("--min-cells", type=int, default=10,
                    help="PDUI mode: min cells with finite PDUI per cluster [10]")
    ap.add_argument("--null-dirs", action="append", default=[], metavar="[NAME=]GLOB",
                    help="per-sample null (label-shuffle) switch outputs; NAME=GLOB, a GLOB "
                         "with {sample}, or one bare GLOB per sample in --sample order")
    ap.add_argument("--null-perms", type=int, default=0,
                    help="PDUI per-cell mode only: number of internal per-sample label "
                         "shuffles to run as the null [0]")
    ap.add_argument("--null-universe", choices=["real", "all"], default="real",
                    help="real: score each null combination only on the (pair, feature) set "
                         "tested in the real run, so the null and the real replication counts "
                         "are over the same universe (a label shuffle changes which features "
                         "pass per-sample filters); all: keep every null feature [real]")
    ap.add_argument("--null-combos", type=int, default=None,
                    help="number of random cross-sample null combinations (default: "
                         "index-aligned, min n_perms)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--note", default="", help="free text embedded in summary.json")
    ap.add_argument("-v", "--verbose", action="store_true")
    return ap


def load_all(samples: list[tuple[str, str]], kind: str, args, prm_pairs, rng=None,
             cache=None) -> tuple[pd.DataFrame, str, str | None]:
    frames, eff_used, strat_used = [], None, None
    for name, path in samples:
        if kind == "diff":
            ldf, strat, eff = load_diff_sample(name, path, strategy=args.strategy,
                                               effect_col=args.effect_col, pairs=prm_pairs)
            strat_used = strat
        else:
            ldf, eff = load_pdui_sample(name, path, pairs=prm_pairs, min_cells=args.min_cells,
                                        effect_col=args.effect_col, rng=rng, _cache=cache)
        if eff_used is None:
            eff_used = eff
        elif eff != eff_used:
            raise ValueError(f"effect column differs between samples ({eff_used} vs {eff})")
        frames.append(ldf)
    long_df = pd.concat(frames, ignore_index=True) if frames else _empty_long()
    return long_df, eff_used, strat_used


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    t0 = time.time()
    samples = [parse_sample_spec(s) for s in args.sample]
    names = [n for n, _ in samples]
    if len(set(names)) != len(names):
        raise SystemExit(f"duplicate sample names: {names}; use NAME=PATH")
    kinds = {detect_input_kind(p) for _, p in samples} if args.input_kind == "auto" else {args.input_kind}
    if len(kinds) != 1:
        raise SystemExit(f"mixed input kinds {kinds}; pass --input-kind")
    kind = kinds.pop()
    pairs = parse_pairs(args.pairs)
    rng = np.random.default_rng(args.seed)
    cache: dict = {}

    long_df, eff, strat = load_all(samples, kind, args, pairs, cache=cache)
    if kind == "diff" and strat == "nb_multi":
        raise SystemExit("nb_multi is an omnibus test without pairs or a signed effect; "
                         "replication needs fisher or nb_pairwise outputs")
    if eff == "log2fc" and args.effect_floor == 0.1:
        log.warning("effect column is log2fc; the default --effect-floor 0.1 was chosen for "
                    "proportion/PDUI deltas -- set it explicitly for log2 fold changes")
    prm = Params(min_samples=args.min_samples, fdr=args.fdr, effect_floor=args.effect_floor,
                 discordant_policy=args.discordant_policy, level=args.level,
                 effect_sign=EFFECT_SIGN.get(eff, +1))
    pas_table = replication_table(long_df, prm, names)
    table = collapse_to_gene(pas_table) if (args.level == "gene" and kind == "diff") else pas_table

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "all_features.tsv", sep="\t", index=False)
    rep = table[table["passes_replication"]] if not table.empty else table
    rep.to_csv(out / "replicated.tsv", sep="\t", index=False)
    real = count_replicated(table)
    log.info("%d samples, %d (pair, feature) rows, %d replicated (K>=%d, same direction), "
             "%d also pass |effect|>=%g", len(names), real["n_features"], real["n_replicated"],
             prm.min_samples, real["n_replicated_floor"], prm.effect_floor)

    # ---- null control -----------------------------------------------------
    null_report: dict = {"enabled": False}
    null_rows = []
    if args.null_dirs or args.null_perms > 0:
        if args.null_dirs:
            null_map = expand_null_specs(args.null_dirs, names)
            n_per = {s: len(v) for s, v in null_map.items()}
            combos = null_combinations(n_per, args.null_combos, rng)
            mode = "external"
        else:
            if kind != "pdui":
                raise SystemExit("--null-perms (internal label shuffles) is only valid for "
                                 "per-cell PDUI input; switch-diff input needs --null-dirs")
            n_per = {s: args.null_perms for s in names}
            combos = [{s: i for s in names} for i in range(args.null_perms)]
            mode = "internal_label_shuffle"
        missing = [s for s, n in n_per.items() if n == 0]
        if missing:
            log.warning("no null outputs for samples %s; they are dropped from the null "
                        "(null replication is then computed over fewer samples and is "
                        "NOT comparable unless K is met by the remaining ones)", missing)
        real_keys = (set(zip(long_df["c1"], long_df["c2"], long_df["feature_id"]))
                     if not long_df.empty else set())
        for ci, combo in enumerate(combos):
            if mode == "external":
                specs = [(s, str(null_map[s][i])) for s, i in combo.items()]
                ndf, _, _ = load_all(specs, kind, args, pairs, cache=cache)
            else:
                sub_rng = np.random.default_rng([args.seed, ci])
                specs = [(s, p) for s, p in samples]
                ndf, _, _ = load_all(specs, kind, args, pairs, rng=sub_rng, cache=cache)
            if args.null_universe == "real" and not ndf.empty:
                ndf = restrict_to_universe(ndf, real_keys)
            nt = replication_table(ndf, prm, names)
            if args.level == "gene" and kind == "diff":
                nt = collapse_to_gene(nt)
            c = count_replicated(nt)
            null_rows.append({"combo": ci, **{f"perm__{s}": i for s, i in combo.items()},
                              "n_features": c["n_features"], "n_replicated": c["n_replicated"],
                              "n_replicated_floor": c["n_replicated_floor"]})
            log.info("null combo %d/%d: %d replicated, %d with floor", ci + 1, len(combos),
                     c["n_replicated"], c["n_replicated_floor"])
        nulldf = pd.DataFrame(null_rows)
        nulldf.to_csv(out / "null_control.tsv", sep="\t", index=False)
        null_report = summarise_null(nulldf, real, mode, n_per, missing)
        null_report["universe"] = args.null_universe
        null_report["n_features_real"] = real["n_features"]
        null_report["n_features_null_mean"] = (float(nulldf["n_features"].mean())
                                               if not nulldf.empty else None)

    summary = {
        "tool": "replication_filter.py", "version": __version__,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "runtime_s": round(time.time() - t0, 2),
        "note": args.note,
        "params": {"min_samples": prm.min_samples, "fdr": prm.fdr, "effect_floor": prm.effect_floor,
                   "discordant_policy": prm.discordant_policy, "level": args.level,
                   "input_kind": kind, "strategy": strat, "effect_col": eff,
                   "effect_sign_convention": ("+ means higher in cluster1 (c1 - c2)"
                                              if prm.effect_sign > 0 else
                                              "+ means higher in cluster2 (c2 vs c1)"),
                   "pairs": sorted(f"{a}_vs_{b}" for a, b in pairs) if pairs else "all",
                   "min_cells": args.min_cells, "seed": args.seed},
        "samples": [{"name": n, "path": str(Path(p).resolve()),
                     "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(p)))}
                    for n, p in samples],
        "n_samples": len(names),
        "pairs_seen": sorted({f"{a}_vs_{b}" for a, b in zip(long_df["c1"], long_df["c2"])}) if not long_df.empty else [],
        "per_sample_called": {
            s: int(((long_df["sample"] == s) & (long_df["qvalue"] < prm.fdr)).sum()) for s in names
        } if not long_df.empty else {},
        "real": real,
        "null": null_report,
        "outputs": {k: str(out / v) for k, v in
                    {"all_features": "all_features.tsv", "replicated": "replicated.tsv",
                     "summary": "summary.json",
                     **({"null_control": "null_control.tsv"} if null_rows else {})}.items()},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=_json_default))
    log.info("wrote %s", out / "summary.json")
    return 0


def restrict_to_universe(long_df: pd.DataFrame, keys: set) -> pd.DataFrame:
    """Keep only rows whose (c1, c2, feature_id) is in `keys` (the real universe)."""
    if long_df.empty or not keys:
        return long_df.iloc[0:0]
    mask = [k in keys for k in zip(long_df["c1"], long_df["c2"], long_df["feature_id"])]
    return long_df[np.asarray(mask, dtype=bool)].reset_index(drop=True)


def summarise_null(nulldf: pd.DataFrame, real: dict, mode: str, n_per: dict, missing: list) -> dict:
    if nulldf.empty:
        return {"enabled": True, "mode": mode, "n_combos": 0, "n_null_per_sample": n_per,
                "samples_without_null": missing}
    n = len(nulldf)
    rep = nulldf["n_replicated"].to_numpy()
    repf = nulldf["n_replicated_floor"].to_numpy()

    def _block(vals, obs):
        obs = int(obs)
        ge = int((vals >= obs).sum())
        return {
            "observed": obs,
            "null_mean": float(vals.mean()), "null_sd": float(vals.std(ddof=1)) if n > 1 else 0.0,
            "null_median": float(np.median(vals)), "null_max": int(vals.max()),
            "null_per_combo": [int(v) for v in vals],
            "empirical_p_ge_observed": (ge + 1) / (n + 1),
            "expected_false_replicated_fraction": (float(vals.mean()) / obs) if obs > 0 else None,
        }
    return {
        "enabled": True, "mode": mode, "n_combos": n, "n_null_per_sample": n_per,
        "samples_without_null": missing,
        "replicated": _block(rep, real["n_replicated"]),
        "replicated_floor": _block(repf, real["n_replicated_floor"]),
        "definition": ("each null combination pairs one independent label-shuffle output per "
                       "sample and runs the identical filter; expected_false_replicated_fraction "
                       "= mean null replicated / observed replicated (an empirical FDR estimate "
                       "for the replicated set, valid to the extent the per-sample null mimics "
                       "the real pipeline)"),
    }


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


if __name__ == "__main__":
    sys.exit(main())
