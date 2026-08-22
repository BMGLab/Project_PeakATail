#!/usr/bin/env python3
"""Isolate the cost of ``--pas-features on`` at the internal-priming seam.

A whole-run wall-clock delta cannot resolve this: on the PBMC chr19+21 slice
the run-to-run spread is ~30 s, an order of magnitude larger than anything a
per-PAS feature computation can cost.  So this measures the seam directly --
the pass over the real caller BEDs with the real genome, with and without the
collector -- and reports wall time and peak RSS for each arm, one arm per
process invocation so ``ru_maxrss`` is attributable.

Usage: bench_features.py --bed A.bed [--bed B.bed] --fasta FA
                         --arm {plain,features,context,append} [--reps N]
"""
from __future__ import annotations

import argparse
import os
import resource
import shutil
import sys
import tempfile
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", action="append", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--arm", required=True,
                    choices=("plain", "features", "context", "append"))
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--support")
    a = ap.parse_args()

    from ema.countmatrix.pas_features import FeatureCollector, append_columns
    from ema.experimental.internal_priming import filter_internal_priming

    n_rows = sum(sum(1 for _ in open(b)) for b in a.bed)
    times = []
    tmpdir = tempfile.mkdtemp(prefix="benchfeat_")
    try:
        for _ in range(a.reps):
            coll = FeatureCollector() if a.arm != "plain" else None
            t0 = time.perf_counter()
            if a.arm in ("plain", "features"):
                for i, bed in enumerate(a.bed):
                    filter_internal_priming(
                        bed, a.fasta, os.path.join(tmpdir, f"out{i}.bed"),
                        mode="filter", features=coll)
            elif a.arm == "context":
                for i, bed in enumerate(a.bed):
                    filter_internal_priming(
                        bed, a.fasta, os.path.join(tmpdir, f"out{i}.bed"),
                        mode="filter", features=coll)
                coll.finish()
            else:  # append
                for i, bed in enumerate(a.bed):
                    filter_internal_priming(
                        bed, a.fasta, os.path.join(tmpdir, f"out{i}.bed"),
                        mode="filter", features=coll)
                feats = coll.finish()
                dst = os.path.join(tmpdir, "pas_support.tsv")
                shutil.copyfile(a.support, dst)
                append_columns(dst, feats)
            times.append(time.perf_counter() - t0)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # kB on Linux
    print(f"{a.arm}\t{n_rows}\t{min(times):.3f}\t"
          f"{sum(times) / len(times):.3f}\t{max(times):.3f}\t{rss}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
