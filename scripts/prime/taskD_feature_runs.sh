#!/usr/bin/env bash
# scripts/prime/taskD_feature_runs.sh -- produce the three FULL-BAM feature tables
# the TASK D transfer test needs (manuscript/24 §3.3 primary protocol).
#
# Each arm runs the peakAtail-prime caller with its branch defaults
# (--read-geometry fixed == v2 geometry, --pas-features on) plus a --genome-fasta
# and WITHOUT --ip-filter, so that:
#   * the call set is v2's (features change no call -- TASK C),
#   * internally-primed candidates are still PRESENT, carrying `ip_tool_flag`,
#     exactly as A3's PBMC universe was built (the veto is then applied by the
#     offline analysis at SELECTION time, which is where it belongs),
#   * the sequence + context columns are real (they ride the single FASTA pass
#     the collector makes when --ip-filter is off).
# Runs are SEQUENTIAL at 8 threads: the box cap for this work is 8.
set -uo pipefail
source "$(dirname "$0")/common.sh"

MOU_GTF=$WD/data/references/mouse/Mus_musculus.GRCm38.102.gtf
MOU_FA=$WD/data/references/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa
M1=$WD/data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam
M2=$WD/data/benchmark/gse104556/starsolo/Mouse2_scRNAseq/Aligned.sortedByCoord.out.bam
RS=$(dirname "$0")/run_slice.sh

set -x
"$RS" --label taskD_feat_m1    --bam "$M1"       --gtf "$MOU_GTF" --seq-len 98 --threads 8 -- --genome-fasta "$MOU_FA"
"$RS" --label taskD_feat_m2    --bam "$M2"       --gtf "$MOU_GTF" --seq-len 98 --threads 8 -- --genome-fasta "$MOU_FA"
"$RS" --label taskD_feat_pbmc  --bam "$PBMC_BAM" --gtf "$HUM_GTF" --seq-len 91 --threads 8 -- --genome-fasta "$HUM_FA"
set +x
echo "TASK D feature runs finished $(date -Is)"
