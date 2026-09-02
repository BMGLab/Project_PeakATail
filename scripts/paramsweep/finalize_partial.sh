#!/usr/bin/env bash
# scripts/paramsweep/finalize_partial.sh <label> -- derive the scorer-ready BEDs of an
# arm whose run was stopped AFTER run/pasbed.bed was final.
#
# Why this exists: `--min-read 0` keeps every barcode, and the *clustering* stage then
# runs Leiden on ~240 k cells for hours.  pasbed.bed is written BEFORE preprocessing /
# clustering (ema/main.py: make_dataframe -> annotate -> write_pas_gene_artifacts, then
# preprocessing), and no stage after it can add, drop or move a PAS.  So the call set is
# already final and the remaining compute is irrelevant to every number in this sweep.
# The derivation below is byte-for-byte run_slice.sh's.
set -uo pipefail
source "$(dirname "$0")/common.sh"
L=${1:?label}; OUT=$OUTROOT/$L
[ -s "$OUT/run/pasbed.bed" ] || { echo "no pasbed.bed for $L"; exit 1; }
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}' "$OUT/run/pasbed.bed" \
  | "$PY" "$NORM" --to ensembl --unmapped drop | sort -k1,1 -k2,2n > "$OUT/pas.bed"
awk -F'\t' '$5>0'  "$OUT/pas.bed" > "$OUT/pas_tier1.bed"
awk -F'\t' '$5==0' "$OUT/pas.bed" > "$OUT/pas_tier2.bed"
awk -F'\t' '$5>=2' "$OUT/pas.bed" > "$OUT/pas_tier1_ge2mol.bed"
printf 'pas %d | tier1 %d | tier2 %d | tier1>=2mol %d\n' \
  "$(wc -l < "$OUT/pas.bed")" "$(wc -l < "$OUT/pas_tier1.bed")" \
  "$(wc -l < "$OUT/pas_tier2.bed")" "$(wc -l < "$OUT/pas_tier1_ge2mol.bed")" | tee "$OUT/counts.txt"
echo "PARTIAL: run stopped after pasbed.bed was final (clustering only); call set is complete. $(date -Is)" > "$OUT/DONE.ok"
