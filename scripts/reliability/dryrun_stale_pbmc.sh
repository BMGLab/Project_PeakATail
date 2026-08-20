#!/usr/bin/env bash
# dryrun_stale_pbmc.sh [outdir] -- MECHANICS DRY-RUN of trusted_novel_pas.py on the STALE
# Stage-2 PBMC output (results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded).
#
# The stale BED carries READS in column 5 (pre-Stage-1c), so criterion 2 is not the
# pre-registered ">= 2 distinct molecules"; the run exists only to prove the pipeline runs
# end to end.  NOTHING it prints is a result.  The final Stage-3 caller run replaces it.
set -euo pipefail
export LC_ALL=C
WD=/mnt/ssd1/Projects/PeakATail_wd
PY=$WD/tools/PeakATail/.venv/bin/python
TOOL=$WD/scripts/reliability/trusted_novel_pas.py
IN=$WD/results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded
OUT=${1:-$WD/results/reliability/_dryrun_stale_stage2_pbmc}
GENOME=$WD/data/references/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ATLAS=$WD/data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6
SIZES=$WD/data/references/chrom.sizes.nochr.filt
TRUTH=$WD/results/benchmark_tools/kinnex_truth/x3p
SEEDS=${SEEDS:-10}
mkdir -p "$OUT"
cat > "$OUT/DRYRUN_STALE_INPUTS_NOT_RESULTS.txt" <<EOF
This directory is a MECHANICS DRY-RUN of scripts/reliability/trusted_novel_pas.py on the
STALE Stage-2 PBMC output ($IN).
- BED column 5 = READS (pre-Stage-1c), not distinct molecules -> --support-unit reads.
- No internal-priming annotation in that run -> flags recomputed from the genome.
Numbers here are NOT results and must not be quoted.  Generated $(date +'%F %T').
EOF
echo "[$(date +'%F %T')] call"
"$PY" "$TOOL" call \
  --pas-bed "$IN/pas.bed" --clip-bed "$IN/pas_tier1.bed" --support-unit reads \
  --genome "$GENOME" --atlas "$ATLAS" --restrict-chroms "$SIZES" \
  --outdir "$OUT/call" 2>&1 | tee "$OUT/call.log"
echo "[$(date +'%F %T')] validate (x3p truth t5/t20/t100/t500, $SEEDS null seeds)"
Q=()
for b in "$OUT"/call/stages/*.bed; do Q+=(--query "$(basename "${b%.bed}")=$b"); done
"$PY" "$TOOL" validate "${Q[@]}" \
  --truth t5="$TRUTH/x3p_truth_t5.point.bed" --truth t20="$TRUTH/x3p_truth_t20.point.bed" \
  --truth t100="$TRUTH/x3p_truth_t100.point.bed" --truth t500="$TRUTH/x3p_truth_t500.point.bed" \
  --chrom-sizes "$SIZES" --seeds "$SEEDS" --stale-dry-run \
  --outdir "$OUT/validate" 2>&1 | tee "$OUT/validate.log"
echo "[$(date +'%F %T')] DRY-RUN DONE (stale inputs; not results)"
