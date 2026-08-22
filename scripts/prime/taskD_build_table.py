#!/usr/bin/env python3
"""TASK D step 1 -- build one candidate feature/truth table per dataset.

A row is one PeakATail candidate PAS, exactly as ``score_tool.py`` would see it
(the run's ``pas.bed``: 3'-most base by strand, Ensembl contig names, restricted
to the contigs of the scoring ``chrom.sizes``).  Columns are

  * every column of the run's ``pas_support.tsv`` -- v2's seven plus the 24
    ``--pas-features on`` columns, i.e. EXACTLY what the caller can compute at
    run time.  Nothing here is derivable only offline: that is the point.  A
    feature the tool cannot produce cannot be shipped.
  * ``width`` from ``pasbed.bed`` (the only per-call quantity in the BED that is
    not in the sidecar), and ``molecules`` = BED column 5.
  * truth distances, which are NEVER features: ``d_atlas`` (nearest same-strand
    PolyASite 2.0 rep site) and, human only, ``d_kin_t5`` / ``d_kin_t20`` /
    ``d_kin_decoy`` (nearest Kinnex long-read 3' end / internal-priming decoy).

It also writes the recall PAIR lists (`bedtools window -sm -w 100`, one row per
truth-site x candidate pair) so recall can be recomputed for any SUBSET of the
candidates without another bedtools call -- the same trick A3's ``evalkit`` used,
re-derived here so nothing depends on A3's scratch directory surviving.

Distances are computed the way ``score_tool.py`` computes them: strand-matched
``bedtools closest -s -d -t first`` on 1-bp points.  For 1-bp points bedtools'
distance is |start_a - start_b|, which is what the pair lists store.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
ENV = dict(os.environ, LC_ALL="C")

REFS = {
    "human": dict(
        atlas=WD / "data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6",
        sizes=WD / "data/references/chrom.sizes.nochr.filt",
        det=WD / "results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed",
    ),
    "mouse": dict(
        atlas=WD / "data/references/atlases/polyasite2.GRCm38.96.rep_sites.bed6",
        sizes=WD / "data/references/mouse/chrom.sizes.filt",
        det=WD / "results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed",
    ),
}
KIN = WD / "results/benchmark_tools/kinnex_truth/x3p"
KIN_FILES = {"t5": KIN / "x3p_truth_t5.point.bed",
             "t20": KIN / "x3p_truth_t20.point.bed",
             "decoy": KIN / "x3p_decoy.point.bed"}


def sh(cmd: str) -> str:
    r = subprocess.run(cmd, shell=True, env=ENV, executable="/bin/bash",
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {cmd}\n{r.stderr[-4000:]}")
    return r.stdout


def closest_dist(calls: Path, ref: Path, work: Path, tag: str) -> pd.Series:
    """pas_id -> distance to the nearest same-strand feature of *ref* (-1 = none)."""
    out = work / f"closest_{tag}.tsv"
    sh(f"bedtools closest -s -d -t first -a {calls} -b {ref} 2>/dev/null "
       f"| awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $4, $NF}}' > {out}")
    d = pd.read_csv(out, sep="\t", header=None, names=["pas_id", tag],
                    dtype={0: str, 1: np.int64})
    return d.groupby("pas_id")[tag].min()


def pair_list(truth: Path, calls: Path, work: Path, tag: str) -> Path:
    """(truth_row, pas_id, distance) for every truth x call pair within 100 bp,
    same strand.  Truth rows are numbered by input line so duplicated names in a
    reference file cannot collapse two distinct sites into one."""
    numbered = work / f"truth_{tag}.bed"
    sh(f"awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,NR,$5,$6}}' {truth} "
       f"| sort -k1,1 -k2,2n > {numbered}")
    out = work / f"pairs_{tag}.tsv"
    sh(f"bedtools window -sm -w 100 -a {numbered} -b {calls} 2>/dev/null "
       f"| awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}}{{d=$2-$8; if(d<0)d=-d; print $4,$10,d}}' > {out}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--run", required=True, help="run tree from run_slice.sh (has pas.bed and run/)")
    ap.add_argument("--species", required=True, choices=("human", "mouse"))
    ap.add_argument("--out", required=True, help="output directory for the table")
    ap.add_argument("--kinnex", action="store_true", help="also join the Kinnex truths (human PBMC only)")
    a = ap.parse_args()

    run = Path(a.run)
    out = Path(a.out)
    work = out / ".work"
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    ref = REFS[a.species]

    # ---- candidates, exactly as score_tool.py sees them --------------------
    calls = work / "calls.bed"
    sh(f"grep -v '^#' {run}/pas.bed | awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' "
       f"{ref['sizes']} - | sort -k1,1 -k2,2n > {calls}")
    bed = pd.read_csv(calls, sep="\t", header=None,
                      names=["chrom", "start", "end", "pas_id", "molecules", "strand"],
                      dtype={"chrom": str, "pas_id": str})
    assert bed.pas_id.is_unique, "pas_id is not unique in pas.bed"

    # width comes from the INTERVAL bed (pas.bed is already collapsed to a point)
    iv = pd.read_csv(run / "run/pasbed.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "pas_id", "molecules", "strand"],
                     dtype={"chrom": str, "pas_id": str})
    iv["width"] = iv.end - iv.start
    bed = bed.merge(iv[["pas_id", "width"]], on="pas_id", how="left")

    sup = pd.read_csv(run / "run/pas_support.tsv", sep="\t", dtype={"pas_id": str},
                      na_values=["NA"], keep_default_na=True)
    # The FULL sidecar is kept, not just the scored subset.  It is the exact
    # population the caller scores at run time (pas_support.tsv is written
    # BEFORE any filter, so it holds every candidate peak calling emitted --
    # 195,940 on mouse 1 against the 85,816 that survive gene assignment into
    # pasbed.bed).  Any within-run statistic (the `rankpct` transform) must be
    # computed over THAT population offline too, or the offline score and the
    # tool's score are two different numbers with one name.
    sup.to_parquet(out / "support_full.parquet", index=False)
    c = bed.merge(sup, on="pas_id", how="left", validate="one_to_one")
    print(f"[cand] {len(c):,} candidates ({len(bed):,} in pas.bed after contig filter); "
          f"sidecar {len(sup):,} rows x {len(sup.columns)} columns")

    # ---- truth (never a feature) -------------------------------------------
    c["d_atlas"] = c.pas_id.map(closest_dist(calls, ref["atlas"], work, "d_atlas")).fillna(-1).astype(np.int64)
    det_pairs = pair_list(ref["det"], calls, work, "det_atlas")
    n_det = int(sh(f"wc -l < {ref['det']}").split()[0])

    meta = {"species": a.species, "run": str(run), "n_candidates": int(len(c)),
            "n_support_rows": int(len(sup)),
            "n_det_atlas": n_det, "atlas": str(ref["atlas"]), "det": str(ref["det"])}

    if a.kinnex:
        for k, p in KIN_FILES.items():
            c["d_kin_" + k] = c.pas_id.map(closest_dist(calls, p, work, "d_kin_" + k)).fillna(-1).astype(np.int64)
        pair_list(KIN_FILES["t20"], calls, work, "kin_t20")
        meta["n_kin_t20"] = int(sh(f"wc -l < {KIN_FILES['t20']}").split()[0])

    c.to_parquet(out / "candidates.parquet", index=False)
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    for src, dst in [(det_pairs, out / "pairs_det_atlas.tsv")] + (
            [(work / "pairs_kin_t20.tsv", out / "pairs_kin_t20.tsv")] if a.kinnex else []):
        os.replace(src, dst)
    print(f"[out] {out}/candidates.parquet  rows {len(c):,}  "
          f"atlas-positive@25bp {int(((c.d_atlas >= 0) & (c.d_atlas <= 25)).sum()):,}")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
