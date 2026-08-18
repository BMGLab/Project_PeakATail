#!/usr/bin/env python3
"""normalize_chroms.py -- bidirectional UCSC('chr') <-> Ensembl chromosome-name mapper.

Part of the benchmark_tools harness (see README.md section 2). Locale-safe:
no locale-dependent operations are used, so it behaves identically under
tr_TR and C locales.

Why: the benchmark contract says every tool's pas.bed uses bare Ensembl names
(1..22, X, Y, MT, accessioned scaffolds like GL000225.1 / KI270728.1), because
the pbmc_10k_v3 BAM, the shared GTF, and every PolyASite reference on this
machine are Ensembl-named. Competitor tools that bundle UCSC annotation print
chr1 / chrM / chrUn_GL000195v1 instead. This script streams a BED/TSV and
rewrites one column.

Mapping rules
  --to ensembl   chr1 -> 1 ... chrX -> X, chrM -> MT,
                 chr14_GL000009v2_random / chrUn_GL000195v1 / chr1_KI270762v1_alt
                   -> GL000009.2 / GL000195.1 / KI270762.1  (accession derived
                   from the UCSC name itself: '<acc>v<n>' -> '<acc>.<n>')
                 already-Ensembl names pass through unchanged (idempotent).
  --to ucsc      1 -> chr1 ... X -> chrX, MT -> chrM,
                 already-UCSC names pass through unchanged (idempotent).
                 Ensembl scaffold -> UCSC is REFUSED (unmapped): the UCSC name
                 embeds placement info (chrUn_ vs chrN_*_random vs *_alt) that
                 the bare accession does not carry; guessing would corrupt
                 coordinates silently. Use an assembly report if you ever need
                 scaffolds in UCSC naming.

Unmapped names are handled per --unmapped: fail (default, exit 65), drop
(remove the line, count to stderr), keep (pass through untouched).

Usage
  ./normalize_chroms.py --to ensembl [--col 1] [--unmapped fail|drop|keep] [IN.bed] > out.bed
  cat in.bed | ./normalize_chroms.py --to ucsc > out.bed
  ./normalize_chroms.py --selftest        # assertion suite, prints OK
  ./normalize_chroms.py --list            # print the primary-chromosome table

Also importable:  from normalize_chroms import to_ensembl, to_ucsc
"""
from __future__ import annotations

import argparse
import re
import sys

PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y"]

# UCSC GRCh38 scaffold names: chr14_GL000009v2_random, chrUn_KI270752v1,
# chr1_KI270762v1_alt, chrX_ML143385v1_fix ...
_UCSC_SCAFFOLD = re.compile(
    r"^(?:\d{1,2}|X|Y|M|MT|Un)_([A-Za-z]{1,4}\d+)v(\d+)(?:_random|_alt|_fix|_decoy)?$"
)
# Ensembl/GenBank accessioned scaffold: GL000009.2, KI270728.1, ML143385.1 ...
_ENS_SCAFFOLD = re.compile(r"^[A-Za-z]{1,4}\d+\.\d+$")


def to_ensembl(name: str) -> str | None:
    """UCSC -> Ensembl. Idempotent on Ensembl names. None if unmappable."""
    if not name.startswith("chr"):
        if name in PRIMARY or name == "MT" or _ENS_SCAFFOLD.match(name):
            return name
        return None
    base = name[3:]
    if base in ("M", "MT"):
        return "MT"
    if base in PRIMARY:
        return base
    m = _UCSC_SCAFFOLD.match(base)
    if m:
        return f"{m.group(1)}.{m.group(2)}"
    return None


def to_ucsc(name: str) -> str | None:
    """Ensembl -> UCSC. Idempotent on UCSC names. None if unmappable
    (notably Ensembl scaffolds -- see module docstring)."""
    if name.startswith("chr"):
        return name
    if name == "MT":
        return "chrM"
    if name in PRIMARY:
        return "chr" + name
    return None


def _selftest() -> None:
    assert to_ensembl("chr1") == "1" and to_ensembl("chr22") == "22"
    assert to_ensembl("chrX") == "X" and to_ensembl("chrY") == "Y"
    assert to_ensembl("chrM") == "MT" and to_ensembl("chrMT") == "MT"
    assert to_ensembl("chrUn_GL000195v1") == "GL000195.1"
    assert to_ensembl("chr14_GL000009v2_random") == "GL000009.2"
    assert to_ensembl("chr1_KI270762v1_alt") == "KI270762.1"
    assert to_ensembl("1") == "1" and to_ensembl("MT") == "MT"     # idempotent
    assert to_ensembl("KI270728.1") == "KI270728.1"                # idempotent
    assert to_ensembl("chrEBV") is None and to_ensembl("weird") is None
    assert to_ucsc("1") == "chr1" and to_ucsc("X") == "chrX"
    assert to_ucsc("MT") == "chrM"
    assert to_ucsc("chr5") == "chr5" and to_ucsc("chrM") == "chrM"  # idempotent
    assert to_ucsc("GL000195.1") is None                            # refused
    # round trips on the primary assembly
    for c in PRIMARY + ["MT"]:
        assert to_ensembl(to_ucsc(c)) == c
    print("normalize_chroms.py selftest OK", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Full docs: open the file, the module docstring is the manual.")
    ap.add_argument("infile", nargs="?", default="-",
                    help="BED/TSV to translate ('-' or absent = stdin)")
    ap.add_argument("--to", choices=("ensembl", "ucsc"), dest="direction",
                    help="target naming scheme")
    ap.add_argument("--col", type=int, default=1,
                    help="1-based column holding the chromosome name (default 1)")
    ap.add_argument("--unmapped", choices=("fail", "drop", "keep"),
                    default="fail",
                    help="what to do with untranslatable names (default fail)")
    ap.add_argument("--selftest", action="store_true", help="run assertions and exit")
    ap.add_argument("--list", action="store_true",
                    help="print the primary-chromosome mapping table and exit")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return 0
    if args.list:
        for c in PRIMARY + ["MT"]:
            print(f"{to_ucsc(c)}\t{c}")
        return 0
    if not args.direction:
        ap.error("--to {ensembl,ucsc} is required (or use --selftest / --list)")

    conv = to_ensembl if args.direction == "ensembl" else to_ucsc
    ci = args.col - 1
    fh = sys.stdin if args.infile == "-" else open(args.infile)
    n_in = n_out = n_drop = 0
    try:
        for line in fh:
            if line.startswith(("#", "track ", "browser ")):
                sys.stdout.write(line)
                continue
            n_in += 1
            fields = line.rstrip("\n").split("\t")
            if ci >= len(fields):
                sys.exit(f"normalize_chroms.py: line {n_in} has < {args.col} columns")
            new = conv(fields[ci])
            if new is None:
                if args.unmapped == "fail":
                    sys.exit(f"normalize_chroms.py: unmappable chromosome "
                             f"{fields[ci]!r} on line {n_in} "
                             f"(--unmapped drop|keep to override); exit 65")
                if args.unmapped == "drop":
                    n_drop += 1
                    continue
                new = fields[ci]  # keep
            fields[ci] = new
            sys.stdout.write("\t".join(fields) + "\n")
            n_out += 1
    finally:
        if fh is not sys.stdin:
            fh.close()
    print(f"normalize_chroms.py: {n_in} rows in, {n_out} out, {n_drop} dropped "
          f"({args.direction})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
