#!/usr/bin/env bash
# score_vs_kinnex.sh <truth_point_bed> <truth_tag> <outroot>
# Atlas-INDEPENDENT scoring of every pbmc_10k_v3 tool PAS BED against Kinnex long-read truth.
set -euo pipefail
export LC_ALL=C
TRUTH=$1; TTAG=$2; OUTROOT=$3
WD=/mnt/ssd1/Projects/PeakATail_wd
GEN=$WD/data/references/chrom.sizes.nochr.filt
GB=$WD/results/benchmark_tools/shared_refs/genebodies.merged.bed
TES=$WD/data/references/atlases/tes.protein_coding.GRCh38.99.bed6
mkdir -p "$OUTROOT"
# restrict truth to the 24 scorable contigs (chrom.sizes.nochr.filt) so precision and
# recall denominators are on the same contig set the tool BEDs are filtered to.
TF="$OUTROOT/$(basename ${TRUTH%.bed}).scorable.bed"
awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' "$GEN" "$TRUTH" | sort -k1,1 -k2,2n > "$TF"
echo "truth PAS on scorable contigs: $(wc -l < "$TF") of $(grep -vc '^#' "$TRUTH")"
for t in peakatail polyapipe scapatrap scapture scutrquant sierra; do
  BED=$WD/results/benchmark_tools/pbmc_10k_v3/$t/pas.bed
  [ -s "$BED" ] || { echo "SKIP $t (no pas.bed)"; continue; }
  echo "=== $t vs $TTAG ==="
  python3 $WD/scripts/benchmark_tools/score_tool.py "$BED" "${t}_vs_${TTAG}" \
    --outdir "$OUTROOT" --workdir "$OUTROOT/.work_${t}" \
    --atlas "$TF" --detected-atlas none --tes "$TES" --genome "$GEN" --genebodies "$GB"
done
