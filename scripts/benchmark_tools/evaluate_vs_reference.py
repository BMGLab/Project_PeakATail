#!/usr/bin/env python3
"""evaluate_vs_reference.py -- uniform scorer for the benchmark_tools harness.

Scores ANY tool's pas.bed (BED6, see README.md section 2) against a reference
PAS point BED: POINT-mode, STRAND-MATCHED `bedtools closest -s -d -t first`,
at cutoffs 10/25/50/100/200 bp (configurable), with an optional shuffled-null
control (`bedtools shuffle -chrom -noOverlapping`, genome-wide and/or
restricted to gene bodies).

PROVENANCE / CITATION
    The matching, point-collapse, histogram, and shuffle logic is adapted from
    /mnt/ssd1/Projects/PeakATail_wd/scripts/manuscript_figures/null_control.py
    (functions make_point / precision / shuffle_set and the AWK_HIST pipeline).
    That script's header documents WHY this geometry is the honest one: with a
    dense reference (PolyASite 3.0 full = 18.4 M sites) most of the genome sits
    within 100 bp of some site, so interval-mode and/or unstranded scoring
    saturates and a shuffled null is mandatory context for any precision
    number. Keep the two scripts consistent if either changes.

STATUS: skeleton for the harness -- the core scoring path below is
implemented and runnable; the deliberately-omitted extras (atlas TPM/class
tiers, detected-gene-restricted recall, class composition of matches) live in
null_control.py sections 3/5/6 and can be ported when the manuscript needs
them per-tool. TODO markers below.

INPUTS
    pas_bed            tool output, BED6. Interval peaks are allowed; they are
                       collapsed to the 3'-most base by strand unless
                       --mode interval is requested.
    --reference        reference PAS POINT BED (BED6+, strand in col 6).
                       Ensembl chrom naming, like every reference on this box.
    --genome           chrom.sizes (two columns). Used to (a) drop off-genome
                       rows, exactly like null_control.prep_real, and (b) feed
                       bedtools shuffle. Default: the 24-chrom Ensembl file in
                       data/references.
    --null             comma list from {genome,genic}, or 'none' (default).
    --gene-bodies      BED of gene spans (merged internally); required for
                       --null genic. Default: data/references/gene_end.bed.

OUTPUT: tidy TSV (stdout or --out), columns
    query series replicate mode metric cutoff_bp n n_matched value
where series is real / N_genome / N_genic, metric is precision (fraction of
called PAS within cutoff of a strand-matched reference site) and, with
--recall, recall (fraction of reference sites within cutoff of a call).

Everything shells out to bedtools/awk/sort with LC_ALL=C (tr_TR-safe) and
streams; the reference is never loaded into python.

Example
  ./evaluate_vs_reference.py results/benchmark_tools/pbmc_10k_v3/peakatail/pas.bed \\
      --reference /mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_hiconf_score0.1.bed \\
      --null genome,genic --recall --out eval.tsv
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
DEFAULT_GENOME = WD / "data/references/chrom.sizes.nochr.filt"
DEFAULT_GENE_BODIES = WD / "data/references/gene_end.bed"
DEFAULT_CUTOFFS = [10, 25, 50, 100, 200]

ENV = dict(os.environ, LC_ALL="C")  # tr_TR locale would corrupt sort order


def sh(cmd: str) -> None:                       # adapted: null_control.sh
    subprocess.run(cmd, shell=True, check=True, env=ENV, executable="/bin/bash")


def sh_out(cmd: str) -> str:                    # adapted: null_control.sh_out
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


def log(*a) -> None:
    print("[eval]", *a, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# core geometry -- adapted from null_control.py
# --------------------------------------------------------------------------
def make_point(iv: Path, pt: Path) -> None:
    """Collapse each row to its 3'-most base by strand = inferred cleavage site.
    Verbatim adaptation of null_control.make_point."""
    sh("awk -F'\\t' 'BEGIN{OFS=\"\\t\"}"
       "{if($6==\"+\"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,$4,$5,$6}' "
       f"{iv} | sort -k1,1 -k2,2n > {pt}")


def hist_awk(cutoffs: list[int]) -> str:
    """awk that turns `bedtools closest -d` output into: n  hits@c1  hits@c2 ...
    Adaptation of null_control.AWK_HIST (cutoffs parameterised)."""
    return (
        "awk -F'\\t' -v CUT=\"" + ",".join(map(str, cutoffs)) + "\" '"
        "BEGIN{nc=split(CUT,cc,\",\")}"
        "{d=$NF; n++; if(d>=0){for(i=1;i<=nc;i++) if(d<=cc[i]) h[i]++}}"
        "END{printf \"%d\", n; for(i=1;i<=nc;i++) printf \"\\t%d\", h[i]+0;"
        " printf \"\\n\"}'")


def closest_hist(a: Path, b: Path, cutoffs: list[int]) -> tuple[int, dict[int, int]]:
    """n rows of `a`; how many have a strand-matched nearest site in `b`
    within each cutoff. Adaptation of null_control.precision()."""
    out = sh_out(f"bedtools closest -s -d -t first -a {a} -b {b} 2>/dev/null "
                 f"| {hist_awk(cutoffs)}").strip().split("\t")
    n, hits = int(out[0]), [int(x) for x in out[1:]]
    return n, dict(zip(cutoffs, hits))


def shuffle_set(iv: Path, out_iv: Path, genome: Path, seed: int,
                incl: Path | None) -> None:
    """Null set with identical n / per-chrom counts / widths / strand.
    Adaptation of null_control.shuffle_set."""
    inc = f"-incl {incl}" if incl else ""
    sh(f"bedtools shuffle -i {iv} -g {genome} -chrom -noOverlapping "
       f"-maxTries 5000 -seed {seed} {inc} 2>/dev/null "
       f"| sort -k1,1 -k2,2n > {out_iv}")


# --------------------------------------------------------------------------
# input preparation
# --------------------------------------------------------------------------
def prep_bed(src: Path, genome: Path, work: Path, tag: str,
             collapse: bool) -> tuple[Path, Path | None, int, int]:
    """Filter to on-genome chroms + sort (null_control.prep_real); optionally
    also make the point-collapsed version. Returns (iv, pt, n_raw, n_kept)."""
    iv = work / f"{tag}.iv.bed"
    n_raw = int(sh_out(f"grep -vc '^#' {src} || true").strip() or 0)
    sh(f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} !/^#/ && ok[$1]' {genome} {src} "
       f"| sort -k1,1 -k2,2n > {iv}")
    n_kept = int(sh_out(f"wc -l < {iv}").strip())
    pt = None
    if collapse:
        pt = work / f"{tag}.pt.bed"
        make_point(iv, pt)
    return iv, pt, n_raw, n_kept


def check_naming(n_raw: int, n_kept: int, src: Path, genome: Path) -> None:
    """Loud failure when the query looks like a chrom-naming mismatch."""
    if n_raw == 0:
        sys.exit(f"evaluate_vs_reference.py: {src} has no data rows")
    if n_kept == 0 or n_kept / n_raw < 0.5:
        sys.exit(
            f"evaluate_vs_reference.py: only {n_kept}/{n_raw} rows of {src} are on "
            f"chromosomes named in {genome}. That is almost certainly a chr-vs-"
            f"Ensembl naming mismatch -- pipe the bed through "
            f"`normalize_chroms.py --to ensembl` first (see README.md sec 2).")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Score a pas.bed against a reference PAS point BED "
                    "(point-mode, strand-matched, shuffled-null option).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="The module docstring is the manual; logic adapted from "
               "scripts/manuscript_figures/null_control.py.")
    ap.add_argument("pas_bed", type=Path, help="tool pas.bed (BED6)")
    ap.add_argument("--reference", type=Path, required=True,
                    help="reference PAS point BED (strand in col 6)")
    ap.add_argument("--genome", type=Path, default=DEFAULT_GENOME,
                    help=f"chrom.sizes for filtering + shuffling "
                         f"(default {DEFAULT_GENOME})")
    ap.add_argument("--cutoffs", default=",".join(map(str, DEFAULT_CUTOFFS)),
                    help="comma list of bp cutoffs (default %(default)s)")
    ap.add_argument("--mode", choices=("point", "interval"), default="point",
                    help="collapse query to 3'-most base (point, default; the "
                         "honest geometry) or score raw intervals")
    ap.add_argument("--null", default="none",
                    help="comma list from {genome,genic}, or 'none' (default)")
    ap.add_argument("--n-shuffles", type=int, default=3,
                    help="replicates per null model (default 3)")
    ap.add_argument("--seed", type=int, default=1,
                    help="first shuffle seed; replicate r uses seed+r-1")
    ap.add_argument("--gene-bodies", type=Path, default=DEFAULT_GENE_BODIES,
                    help="gene-span BED for --null genic (merged internally)")
    ap.add_argument("--collapse-reference", action="store_true",
                    help="also 3'-collapse the reference (use when the "
                         "reference has interval rows, e.g. raw cluster beds)")
    ap.add_argument("--recall", action="store_true",
                    help="additionally report recall (reference -> query)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output TSV (default stdout)")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="keep intermediates here (default: a temp dir, removed)")
    args = ap.parse_args()

    cutoffs = sorted({int(c) for c in args.cutoffs.split(",")})
    nulls = [] if args.null in ("", "none") else args.null.split(",")
    for m in nulls:
        if m not in ("genome", "genic"):
            ap.error(f"unknown null model {m!r} (allowed: genome, genic, none)")
    if not shutil.which("bedtools"):
        sys.exit("bedtools not on PATH")

    own_tmp = args.workdir is None
    work = Path(tempfile.mkdtemp(prefix="eval_vs_ref.")) if own_tmp else args.workdir
    work.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []

    def record(series: str, rep: int, metric: str, n: int,
               hits: dict[int, int]) -> None:
        for c in cutoffs:
            rows.append(dict(query=str(args.pas_bed), series=series,
                             replicate=rep, mode=args.mode, metric=metric,
                             cutoff_bp=c, n=n, n_matched=hits[c],
                             value=hits[c] / n if n else float("nan")))

    try:
        # ---- query ----------------------------------------------------------
        iv, pt, n_raw, n_kept = prep_bed(args.pas_bed, args.genome, work,
                                         "query", collapse=True)
        check_naming(n_raw, n_kept, args.pas_bed, args.genome)
        query = pt if args.mode == "point" else iv
        log(f"query {args.pas_bed}: {n_kept:,} rows "
            f"({n_raw - n_kept} off-genome dropped), mode={args.mode}")

        # ---- reference ------------------------------------------------------
        riv, rpt, rn_raw, rn_kept = prep_bed(args.reference, args.genome, work,
                                             "ref",
                                             collapse=args.collapse_reference)
        check_naming(rn_raw, rn_kept, args.reference, args.genome)
        ref = rpt if args.collapse_reference else riv
        log(f"reference {args.reference}: {rn_kept:,} sites")

        # ---- real precision (and recall) ------------------------------------
        n, hits = closest_hist(query, ref, cutoffs)
        record("real", 0, "precision", n, hits)
        log("precision  real      " +
            " ".join(f"{c}:{hits[c]/n:.4f}" for c in cutoffs))
        if args.recall:
            n, hits = closest_hist(ref, query, cutoffs)
            record("real", 0, "recall", n, hits)
            log("recall     real      " +
                " ".join(f"{c}:{hits[c]/n:.4f}" for c in cutoffs))

        # ---- shuffled nulls --------------------------------------------------
        genic = None
        if "genic" in nulls:
            genic = work / "genebodies.merged.bed"
            sh(f"sort -k1,1 -k2,2n {args.gene_bodies} "
               f"| bedtools merge -i - > {genic}")
        for model in nulls:
            series = "N_genome" if model == "genome" else "N_genic"
            for rep in range(1, args.n_shuffles + 1):
                siv = work / f"null.{model}.s{rep}.iv.bed"
                shuffle_set(iv, siv, args.genome, args.seed + rep - 1,
                            genic if model == "genic" else None)
                if args.mode == "point":
                    spt = work / f"null.{model}.s{rep}.pt.bed"
                    make_point(siv, spt)
                    sq = spt
                else:
                    sq = siv
                n, hits = closest_hist(sq, ref, cutoffs)
                record(series, rep, "precision", n, hits)
                log(f"precision  {series:9s} rep{rep} " +
                    " ".join(f"{c}:{hits[c]/n:.4f}" for c in cutoffs))

        # TODO(harness): port atlas tiers (TPM>=1 / TE-AL-EX class filters),
        # detected-gene-restricted recall, and match-class composition from
        # null_control.py sections 3/5/6 when per-tool figures need them.

        # ---- write ----------------------------------------------------------
        cols = ["query", "series", "replicate", "mode", "metric", "cutoff_bp",
                "n", "n_matched", "value"]
        out = sys.stdout if args.out is None else open(args.out, "w")
        try:
            out.write("\t".join(cols) + "\n")
            for r in rows:
                r["value"] = f"{r['value']:.6f}"
                out.write("\t".join(str(r[c]) for c in cols) + "\n")
        finally:
            if out is not sys.stdout:
                out.close()
                log(f"wrote {args.out} ({len(rows)} rows)")
    finally:
        if own_tmp:
            shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
