#!/usr/bin/env bash
# scripts/paramsweep/drive.sh <grid.tsv> [species_filter] -- run every arm in the
# grid through scripts/prime/run_slice.sh (symlinked here so the output roots are
# the paramsweep sandbox), two arms at a time, 8 threads each == the 16-thread cap.
set -uo pipefail
source "$(dirname "$0")/common.sh"
GRID=${1:?grid.tsv}; FILT=${2:-}
MM_BAM=/mnt/ssd0/emaout/peakatail_benchmark/prime/_slice_mouse/mouse1_chr18_19.bam
MM_GTF=$WD/data/references/mouse/Mus_musculus.GRCm38.102.gtf
MM_FA=$WD/data/references/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa
RUN=$WD/scripts/paramsweep/run_slice.sh
CONC=${CONC:-2}
run_one() {
  local label=$1 species=$2 flags=$3
  [ -f "$OUTROOT/$label/DONE.ok" ] && { echo "[skip] $label already done"; return 0; }
  if [ "$species" = mouse ]; then
    bash "$RUN" --label "$label" --code "$V2_CODE" --threads 8 \
      --bam "$MM_BAM" --gtf "$MM_GTF" --seq-len 98 -- \
      --ip-filter --genome-fasta "$MM_FA" --ip-filter-mode filter $flags \
      > "$OUTROOT/$label.log" 2>&1
  else
    bash "$RUN" --label "$label" --code "$V2_CODE" --threads 8 -- \
      --ip-filter --genome-fasta "$HUM_FA" --ip-filter-mode filter $flags \
      > "$OUTROOT/$label.log" 2>&1
  fi
  echo "[done] $label rc=$? $(cat "$OUTROOT/$label/counts.txt" 2>/dev/null)"
}
while IFS=$'\t' read -r label species block flags; do
  case "$label" in \#*|"") continue;; esac
  [ -n "$FILT" ] && [ "$species" != "$FILT" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge "$CONC" ]; do wait -n; done
  echo "[launch] $label ($species, $block): $flags"
  run_one "$label" "$species" "$flags" &
done < "$GRID"
wait
echo "ALL ARMS COMPLETE $(date -Is)"
