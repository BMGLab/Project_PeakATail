#!/usr/bin/env python3
"""scripts/prime/identity_check.py REF_DIR NEW_DIR [-o report.tsv] [--strict-volatile]

Byte-identity checker for the peakAtail-prime three-way comparison.

The cardinal rule of this branch is that a compatibility mode must reproduce the
v2 (commit 9dfdefb) outputs BYTE-FOR-BYTE.  This script is how that is checked.
It walks two caller output trees, pairs files by relative path and reports, for
every path, one of:

  SAME            md5 identical
  DIFF            md5 differs                                <- FAILS the check
  ONLY_IN_REF     present in REF_DIR only                    <- FAILS the check
  ONLY_IN_NEW     present in NEW_DIR only                    <- FAILS the check
  SAME_NORM       byte-different but identical after normalising the run path and
                  ISO timestamps out of the text -- only ever applied to the
                  CONFIG set below, and reported, never silently accepted
  DIFF_NORM       config file that still differs after normalisation      <- FAILS
  SAME_LOG        run journal identical once timings and worker order are normalised
  LOG_DIFF        run journal still differs -- printed with a sample, NOT fatal

VOLATILE (path/clock-bearing by construction, never a statement about the data):
  peakatail_*.log, run_config.json, run_manifest.json, *.log
For these the raw bytes cannot match across two runs in different directories, so
they are compared after substituting the two run roots with a placeholder and
blanking ISO-8601 timestamps, elapsed seconds and run ids.  A DIFF_NORM there is
a real finding (e.g. a resolved config value drifted) and fails the check.

HDF5 (.h5ad) files are compared byte-wise first; if they differ the script drops
to a dataset-by-dataset numeric comparison with h5py (HDF5 containers are not
guaranteed byte-stable) and reports SAME_H5/DIFF_H5.  SAME_H5 does not fail.

Exit status: 0 iff every file is SAME / SAME_NORM / SAME_H5 / SAME_LOG / LOG_DIFF
and neither side has an unmatched file.  Only run journals may differ.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

# Three tiers, on purpose:
#   DATA    -- everything not listed below.  Must be byte-identical.  A difference FAILS.
#   CONFIG  -- carries semantic content (resolved config, artifact hashes, job list) wrapped in
#              paths/clock values that can never match.  Must be identical after normalisation;
#              a residual difference FAILS, because it means a real setting drifted.
#   LOG     -- the human-readable run journal.  Full of elapsed times, per-worker RSS and
#              nondeterministic worker completion order.  Reported, never fatal; the residual
#              normalised line difference is printed so drift is still visible.
VOLATILE_GLOBS = ("run_config.json", "run_manifest.json", "*.peak_jobs.json")
LOG_GLOBS = ("peakatail_*.log", "*.log")
ISO = re.compile(rb"\d{4}-\d{2}-\d{2}[T_ ]\d{2}[:-]\d{2}[:-]\d{2}(\.\d+)?")
RUNID = re.compile(rb"run_\d{4}-\d{2}-\d{2}_\d{6}")
FLOATSEC = re.compile(rb"\d+\.\d{3,}s")
# per-job telemetry the caller writes into the output tree (peakcalling/*.peak_jobs.json):
# wall-clock seconds and peak RSS can never be byte-stable between two runs.
TELEMETRY = re.compile(rb'"(wall_s|maxrss_mb|elapsed_s|rss_mb|peak_rss_mb)":\s*[0-9.eE+-]+')
# a timestamped log FILENAME (peakatail_2026-08-21_223727.log) -- two runs never share it,
# so pair such files by their normalised name instead of reporting both as unmatched.
TS_NAME = re.compile(r"_\d{4}-\d{2}-\d{2}_\d{6}")


def md5(p: Path) -> str:
    h = hashlib.md5()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _matches(rel: Path, globs) -> bool:
    return any(rel.match(g) or rel.name == g for g in globs)


def is_volatile(rel: Path) -> bool:
    return _matches(rel, VOLATILE_GLOBS)


def is_log(rel: Path) -> bool:
    return _matches(rel, LOG_GLOBS)


ELAPSED = re.compile(rb"elapsed: [0-9.]+s")
RSSLINE = re.compile(rb"\d+(\.\d+)? (s|GB RSS|MB RSS)")


def log_lines(raw: bytes, roots: list[bytes]) -> list[bytes]:
    """Normalise a run log to the set of lines that carry meaning, order-insensitively."""
    raw = normalise(raw, roots)
    raw = ELAPSED.sub(b"elapsed: <SEC>", raw)
    raw = RSSLINE.sub(b"<NUM>", raw)
    return sorted(ln.rstrip() for ln in raw.splitlines() if ln.strip())


def normalise(raw: bytes, roots: list[bytes]) -> bytes:
    for r in roots:
        raw = raw.replace(r, b"<RUNROOT>")
    raw = ISO.sub(b"<TS>", raw)
    raw = RUNID.sub(b"<RUNID>", raw)
    raw = FLOATSEC.sub(b"<SEC>", raw)
    raw = TELEMETRY.sub(lambda m: b'"' + m.group(1) + b'": <NUM>', raw)
    # sha256 content hashes inside run_manifest.json cover files we compare directly
    raw = re.sub(rb"sha256:[0-9a-f]{64}", b"sha256:<H>", raw)
    return raw


def h5_equal(a: Path, b: Path) -> tuple[bool, str]:
    try:
        import h5py
        import numpy as np
    except Exception as exc:                                   # pragma: no cover
        return False, f"h5py unavailable ({exc})"
    diffs: list[str] = []

    def walk(ga, gb, prefix=""):
        ka, kb = set(ga.keys()), set(gb.keys())
        if ka != kb:
            diffs.append(f"{prefix}: key sets differ {sorted(ka ^ kb)[:5]}")
            return
        for k in sorted(ka):
            va, vb = ga[k], gb[k]
            if hasattr(va, "keys"):
                walk(va, vb, f"{prefix}/{k}")
            else:
                aa, bb = va[()], vb[()]
                try:
                    same = bool(np.array_equal(aa, bb))
                except Exception:
                    same = bytes(aa) == bytes(bb)
                if not same:
                    diffs.append(f"{prefix}/{k}")

    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        walk(fa, fb)
    return (not diffs), ("; ".join(diffs[:5]) if diffs else "all datasets equal")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("-o", "--out", default=None, help="write a TSV report here")
    ap.add_argument("--strict-volatile", action="store_true",
                    help="require raw byte equality even for the volatile set")
    a = ap.parse_args()
    ref, new = Path(a.ref).resolve(), Path(a.new).resolve()
    for d in (ref, new):
        if not d.is_dir():
            sys.exit(f"not a directory: {d}")
    roots = [str(ref).encode(), str(new).encode()]

    def rels(root: Path) -> dict[str, Path]:
        # key = relative path with any embedded run timestamp normalised away, so
        # peakatail_<ts>.log from two runs is compared rather than reported twice.
        out: dict[str, Path] = {}
        for p in root.rglob("*"):
            if p.is_file() and "__pycache__" not in p.parts:
                rel = p.relative_to(root)
                out[TS_NAME.sub("_<TS>", str(rel))] = rel
        return out

    ma, mb = rels(ref), rels(new)
    ra, rb = set(ma), set(mb)
    rows: list[tuple[str, str, str]] = []
    for key in sorted(ra - rb):
        rows.append(("ONLY_IN_REF", str(ma[key]), ""))
    for key in sorted(rb - ra):
        rows.append(("ONLY_IN_NEW", str(mb[key]), ""))

    for key in sorted(ra & rb):
        rel = ma[key]
        pa, pb = ref / ma[key], new / mb[key]
        ha, hb = md5(pa), md5(pb)
        if ha == hb:
            rows.append(("SAME", str(rel), ha))
            continue
        relb = Path(mb[key])
        if is_log(relb) and not a.strict_volatile:
            la, lb = log_lines(pa.read_bytes(), roots), log_lines(pb.read_bytes(), roots)
            if la == lb:
                rows.append(("SAME_LOG", str(rel), "identical after normalisation"))
            else:
                extra = len(set(la) ^ set(lb))
                sample = sorted(set(la) ^ set(lb))[:3]
                rows.append(("LOG_DIFF", str(rel),
                             "%d normalised lines differ, e.g. %s" %
                             (extra, b" | ".join(sample)[:160].decode("utf8", "replace"))))
            continue
        if is_volatile(relb) and not a.strict_volatile:
            na = normalise(pa.read_bytes(), roots)
            nb = normalise(pb.read_bytes(), roots)
            rows.append(("SAME_NORM" if na == nb else "DIFF_NORM", str(rel),
                         "%s != %s" % (ha[:8], hb[:8])))
            continue
        if rel.suffix == ".h5ad":
            ok, why = h5_equal(pa, pb)
            rows.append(("SAME_H5" if ok else "DIFF_H5", str(rel), why))
            continue
        rows.append(("DIFF", str(rel), "%s != %s (%d vs %d bytes)"
                     % (ha[:8], hb[:8], pa.stat().st_size, pb.stat().st_size)))

    order = ["DIFF", "DIFF_NORM", "DIFF_H5", "ONLY_IN_REF", "ONLY_IN_NEW",
             "LOG_DIFF", "SAME_NORM", "SAME_H5", "SAME_LOG", "SAME"]
    counts = {k: 0 for k in order}
    for st, _, _ in rows:
        counts[st] = counts.get(st, 0) + 1

    print(f"REF {ref}\nNEW {new}\n")
    for st in order:
        if counts.get(st):
            print(f"  {st:<12} {counts[st]}")
    bad = [r for r in rows if r[0] in ("DIFF", "DIFF_NORM", "DIFF_H5",
                                       "ONLY_IN_REF", "ONLY_IN_NEW")]
    if bad:
        print("\nfiles that are NOT identical:")
        for st, rel, why in bad:
            print(f"  {st:<12} {rel}  {why}")
    info = [r for r in rows if r[0] == "LOG_DIFF"]
    if info:
        print("\nrun-journal differences (timings / worker order; NOT fatal):")
        for st, rel, why in info:
            print(f"  {st:<12} {rel}  {why}")
    soft = [r for r in rows if r[0] in ("SAME_NORM", "SAME_H5", "SAME_LOG")]
    if soft:
        print("\nidentical only after normalisation / structural comparison:")
        for st, rel, why in soft:
            print(f"  {st:<12} {rel}  {why}")

    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write("status\tpath\tdetail\n")
            for st, rel, why in rows:
                fh.write(f"{st}\t{rel}\t{why}\n")
        print(f"\nreport: {a.out}")

    verdict = "IDENTICAL" if not bad else "NOT IDENTICAL"
    print(f"\nVERDICT: {verdict}  ({len(rows)} files compared)")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
