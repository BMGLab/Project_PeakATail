#!/usr/bin/env bash
# scripts/prime/make_slice.sh [out.bam] -- ensure the team's standard dev slice exists.
#
# The slice is PBMC 10k v3 restricted to contigs 19 and 21 (~7% of the BAM,
# 41,608,052 + 9,290,404 records).  The perf work already built one at
# /mnt/ssd0/emaout/peakatail_benchmark/perf/slice/pbmc_chr19_21.bam; this script
# REUSES it after verifying the contig content, and only re-extracts if it is
# missing.  Verification, not assumption: samtools idxstats must show exactly
# contigs 19 and 21 carrying reads.
set -uo pipefail
source "$(dirname "$0")/common.sh"
OUT=${1:-$SLICE_BAM}

verify() {
  local bam=$1
  [ -s "$bam" ] && [ -s "$bam.bai" ] || return 1
  local got
  got=$(samtools idxstats "$bam" | awk -F'\t' '$3>0{print $1}' | sort | tr '\n' ',')
  echo "[slice] contigs with reads: $got"
  [ "$got" = "19,21," ]
}

if verify "$OUT"; then
  echo "[slice] REUSING existing slice: $OUT"
  samtools idxstats "$OUT" | awk -F'\t' '$3>0{printf "  %s\t%d reads\n",$1,$3}'
  mkdir -p "$ARTROOT"; md5sum "$OUT" | tee "$ARTROOT/slice_bam.md5"
  exit 0
fi

echo "[slice] building $OUT from $PBMC_BAM (contigs 19,21)"
mkdir -p "$(dirname "$OUT")"
samtools view -@ 4 -b -o "$OUT" "$PBMC_BAM" 19 21 || exit 1
samtools index -@ 4 "$OUT" || exit 1
verify "$OUT" || { echo "[slice] VERIFY FAILED"; exit 1; }
mkdir -p "$ARTROOT"; md5sum "$OUT" | tee "$ARTROOT/slice_bam.md5"
