#!/usr/bin/env python3
"""
convert_scapture.py -- turn SCAPTURE `*.peaks.evaluated.bed` (BED12+3) into
1-bp PAS point BED6 files for scripts/benchmark_tools/score_tool.py.

GEOMETRY (proven from SCAPTURE source, not assumed -- see the header block
written into every output file, and the docstring below).

Outputs
  pas.bed      PRIMARY. DeepPASS-positive peaks from the exonic + intronic
               classes = the tool's own documented recommended PAS selection
               (README.md "Run scapture PASquant module", options 1/2).
  pas_all.bed  PERMISSIVE. Every evaluated peak, all three classes
               (exonic + intronic + 3primeExtended), no DeepPASS filter.

Both are BED6, 0-based half-open, 1-bp intervals, coordinate-sorted, and
deduplicated on (chrom, start, end, strand) -- matching the convention already
used by sierra/pas.bed, peakatail/pas.bed and polyapipe/pas.bed in this repo.
"""
import argparse
import gzip
import sys
from pathlib import Path

CLASSES = ("exonic", "intronic", "3primeExtended")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def load(path, cls):
    """Yield (chrom, pas_start, pas_end, name, strand, deeppass) per peak."""
    with opener(path) as fh:
        for ln, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                sys.exit(f"{path}:{ln}: expected >=15 columns, got {len(f)}")
            chrom, start, end, name, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]
            deeppass = f[13]
            if strand == "+":
                s, e = end - 1, end          # 3'-most base of the peak
            elif strand == "-":
                s, e = start, start + 1      # 3'-most base of the peak
            else:
                sys.exit(f"{path}:{ln}: unexpected strand {strand!r}")
            if s < 0:
                sys.exit(f"{path}:{ln}: negative PAS start")
            yield chrom, s, e, f"{name}|{cls}|{deeppass}", strand, deeppass


def dedup(recs):
    """Coordinate-sort; collapse identical (chrom,start,end,strand).

    Tie-break inside a point group: DeepPASS-positive before negative, then
    exonic < intronic < 3primeExtended, then peak name. So a surviving record
    keeps the positive call, and keeps the exonic/intronic label whenever any
    exonic/intronic peak shares the point. That makes pas.bed both a subset of
    pas_all.bed's points AND exactly reproducible from pas_all.bed by filtering
    on the name's |<class>|<call> suffix.
    """
    rank = {c: i for i, c in enumerate(CLASSES)}
    recs = sorted(recs, key=lambda r: (r[0], r[1], r[2], r[4],
                                       0 if r[5] == "positive" else 1,
                                       rank[r[3].split("|")[-2]], r[3]))
    out, seen = [], set()
    for r in recs:
        k = (r[0], r[1], r[2], r[4])
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out


def write(path, recs, header):
    with open(path, "w") as fh:
        fh.write(header)
        for chrom, s, e, name, strand, _ in recs:
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\n")


COMMON = """\
# GEOMETRY PROOF -- SCAPTURE v1.0 (2021/01/25), repo tools/SCAPTURE.
#
# 1) THE PEAK. Columns 1-12 of *.peaks.evaluated.bed are a BED12 spliced peak
#    region (README.md: "the 1-12 columns represent the spliced peak region of
#    PAS"). Columns 13/14/15 are, per README.md:
#      13 = "number of refernce poly(A) sites supporting the PAS"
#      14 = "DeepPASS prediction of the PAS"  (positive|negative)
#      15 = "+/-100 nt sequence around the CLEAVAGE SITE of PAS"
#
# 2) WHERE THE CLEAVAGE SITE IS. scapture_evaluate line 95 builds the col-15
#    window that README.md calls "around the cleavage site":
#      if($F[5] eq "+"){ $F[1]=$F[2]-100; $F[2]=$F[2]+100; }
#      if($F[5] eq "-"){ $F[2]=$F[1]+100; $F[1]=$F[1]-100; }
#    i.e. the window is anchored on BED chromEnd for '+' and BED chromStart for
#    '-' -- the 3' terminus of the peak in transcript orientation. The
#    known-PAS overlap window (scapture_evaluate line 108) is anchored the same
#    way and is deliberately asymmetric about it (+: [end-50,end+25);
#    -: [start-25,start+50) in transcript orientation), i.e. 50 nt upstream /
#    25 nt downstream of that same terminus -- consistent with it being the
#    cleavage position.
#    => PAS = 3' terminus of the peak, by strand. Confirmed empirically on 120
#    randomly sampled peaks (40 per class, seed 7, both strands, 11 of them
#    spliced with blockCount>1): the col-15 string is byte-identical to the
#    genome slice this rule predicts in 120/120 cases, while the off-by-one
#    alternative (window centred on the terminal base rather than on its 3'
#    edge) matched 0/120. The window is genomic, not spliced, even for
#    multi-block peaks -- line 95 forces $F[9]=1 before bedtools getfasta.
#
# 3) BED CONVERSION. The cleavage boundary sits at the outer edge of the
#    peak's 3'-most base, so that base as a 0-based half-open 1-bp interval is
#      strand '+': [chromEnd-1, chromEnd)
#      strand '-': [chromStart, chromStart+1)
#    identical to the reduction score_tool.py make_point() applies
#    ('+': end-1; '-': start), so these files are idempotent under it.
#    NOTE thickStart/thickEnd (cols 7-8) are NOT the PAS: scapture_callpeak
#    line 212 sets both to chromEnd for intronic/3primeExtended peaks on BOTH
#    strands, and exonic peaks inherit thick == chromStart/chromEnd from
#    genePredToBed. They are placeholders; do not use them.
#
# COLUMNS: chrom, start(0-based), end, name, score, strand.
#   name  = <SCAPTURE peak name>|<peak class>|<DeepPASS call>, where the
#           SCAPTURE peak name is gene|index|gene|gene_type|transcript|region
#           (region = 3UTR|5UTR|CDS|exon|intron|3primeExtended) and peak class
#           is the source file (exonic|intronic|3primeExtended).
#   score = 0 for every row. SCAPTURE's own BED score column is 0 for all
#           162,346 evaluated peaks, and no per-PAS expression is available
#           because this run's PASquant quantification produced an EMPTY UMI
#           matrix (see CAVEAT below). Score carries no information -- do not
#           rank on it.
# Deduplicated on (chrom,start,end,strand). Tie-break within a point:
# DeepPASS-positive before negative, then exonic < intronic < 3primeExtended,
# then peak name -- so pas.bed is exactly reproducible from pas_all.bed by
# filtering on the name's trailing |<class>|<call>.
# Chromosomes are Ensembl-style (no 'chr'); scaffolds and MT are retained,
# score_tool.py drops chroms absent from its genome file itself.
#
# CAVEAT (does not affect these files): PAScall succeeded, but the PASquant
# module silently failed. `samtools sort` never merged its output -- 364
# pbmc10k.PASquant.KeepCell.reassigned.bam.tmp.NNNN.bam shards (36 GB) remain
# and the merged pbmc10k.PASquant.KeepCell.reassigned.bam was never written
# (/mnt/ssd1 is at 100% capacity, the likely cause). samtools index and then
# umi_tools raised FileNotFoundError, but scapture_quant discards those errors
# (&> /dev/null) and still exited 0, so DONE.ok is misleading.
# pbmc10k.PASquant.KeepCell.UMIs.tsv.gz therefore has 36,204 PAS rows and ZERO
# cell columns. PAS CALLING (what these files encode) is unaffected;
# single-cell QUANTIFICATION for this run is unusable and must be rerun.
#
# RUN: scapture -m PAScall -l 91 -p 16 --species human, polyaDB = NULL.
#   BAM   data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
#   GTF   Ensembl GRCh38 (annotation module, --extend 2000)
#   FASTA Homo_sapiens.GRCh38.dna.primary_assembly.fa
#   PAScall wall 8:57:58, maxRSS 23.0 GiB (24,160,272 kB), 16 threads.
# Because --polyaDB was NULL, column 13 is "NA" for every peak, so SCAPTURE's
# recommended selection (README option 1, "$F[12] > 0 | $F[13] eq positive")
# collapses exactly onto option 2 ("$F[13] eq positive"). They are the same set
# here; no known-PAS filter was or could be applied.
# Evaluated peaks by class: exonic 54,258 (23,962 DeepPASS-positive),
#   intronic 100,148 (12,242 positive), 3primeExtended 7,940 (885 positive).
"""

PRIMARY = """\
# SCAPTURE v1.0 PAS points -- PRIMARY (tool's documented recommended output).
# SELECTION: DeepPASS-positive peaks from the exonic + intronic classes only.
#   perl -alne 'print @F[0..11] if $F[13] eq "positive";' \\
#     pbmc10k.exonic.peaks.evaluated.bed pbmc10k.intronic.peaks.evaluated.bed
# This is verbatim README.md "Run scapture PASquant module" option 2, which
# here equals recommended option 1 (see CAVEAT on polyaDB below), and is
# exactly the set this run fed to PASquant (36,204 peaks -> KeepPAS.bed).
# The 3primeExtended class is EXCLUDED because all three documented selection
# recipes in README.md read only the exonic and intronic files. Its 885
# DeepPASS-positive peaks land on 821 points that no exonic/intronic peak
# already covers; a scope-matched permissive/primary pair can be built by
# appending them:
#   awk -F'\\t' '$4 ~ /\\|3primeExtended\\|positive$/' pas_all.bed
# yields exactly those 821 rows (35,773 + 821 = 36,594 points).
"""

ALL = """\
# SCAPTURE v1.0 PAS points -- PERMISSIVE (all evaluated peaks).
# SELECTION: every peak in all three evaluated files (exonic + intronic +
# 3primeExtended), with NO DeepPASS filter. Superset of pas.bed.
# Reconstruct README.md option 3 ("select all raw PASs", exonic + intronic
# only) with:  awk -F'\\t' '$4 !~ /\\|3primeExtended\\|/' pas_all.bed
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="SCAPTURE output directory")
    ap.add_argument("--prefix", required=True, help="run prefix, e.g. pbmc10k")
    a = ap.parse_args()
    d = Path(a.dir)

    peaks = {}
    for cls in CLASSES:
        p = d / f"{a.prefix}.{cls}.peaks.evaluated.bed"
        if not p.exists():
            sys.exit(f"missing {p}")
        peaks[cls] = list(load(p, cls))
        print(f"  {cls}: {len(peaks[cls]):,} peaks, "
              f"{sum(1 for r in peaks[cls] if r[5] == 'positive'):,} positive")

    prim = dedup([r for cls in ("exonic", "intronic") for r in peaks[cls]
                  if r[5] == "positive"])
    allr = dedup([r for cls in CLASSES for r in peaks[cls]])

    write(d / "pas.bed", prim,
          PRIMARY + COMMON + f"# n_points={len(prim)} (deduplicated)\n")
    write(d / "pas_all.bed", allr,
          ALL + COMMON + f"# n_points={len(allr)} (deduplicated)\n")

    pts = {(r[0], r[1], r[2], r[4]) for r in allr}
    assert all((r[0], r[1], r[2], r[4]) in pts for r in prim), \
        "pas.bed is not a subset of pas_all.bed"
    print(f"  pas.bed     {len(prim):,} points")
    print(f"  pas_all.bed {len(allr):,} points")


if __name__ == "__main__":
    main()
