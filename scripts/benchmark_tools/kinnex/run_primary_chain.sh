#!/usr/bin/env bash
# Full primary (Kinnex PBMC 10x 3' v3.1, 81,678,354 dedup FLNC molecules) Tier-1 chain.
# Fires as soon as the minimap2 alignment marker appears.
set -euo pipefail
export LC_ALL=C
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
WD=/mnt/ssd1/Projects/PeakATail_wd
R=$WD/results/benchmark_tools/kinnex_truth
TAG=x3p
BAM=$K/kinnex_primary.dedup.GRCh38.bam
rm -f $K/PRIMARY_CHAIN.DONE.ok $K/PRIMARY_CHAIN.FAILED.err
trap 'echo "FAILED rc=$? at line $LINENO" > $K/PRIMARY_CHAIN.FAILED.err' ERR
while [ ! -f $K/ALIGN_PRIMARY.DONE.ok ]; do
  [ -f $K/ALIGN_PRIMARY.FAILED.err ] && { echo "alignment failed" > $K/PRIMARY_CHAIN.FAILED.err; exit 1; }
  sleep 60
done
echo "$(date +'%F %T') alignment done, starting extraction"
# Step 2/3: per-molecule 3' termini. Ensembl contigs already (STRIP=0); real cells = the
# 12,852 genes_seurat barcodes, CB orientation AS-IS (verified against the whitelist).
samtools view -@ 3 -F 0x904 -q 1 "$BAM" \
 | gawk -v WL=$K/work/whitelist_asis.txt -v STRIP=0 -v CLIP=30 -f $K/termini.awk \
   > $K/work/primary_termini.tsv 2> $K/logs/primary_termini.stats
cat $K/logs/primary_termini.stats
echo "$(date +'%F %T') prep"
$K/build_prep.sh $K/work/primary_termini.tsv $K/scratch/${TAG}_bt
echo "$(date +'%F %T') truth"
$K/build_truth2.sh $TAG $K/scratch/${TAG}_bt $R/$TAG
cat $K/logs/primary_termini.stats >> $R/$TAG/${TAG}_filter_report.txt
echo "$(date +'%F %T') scoring"
for T in 5 20 100 500; do
  $K/score_vs_kinnex.sh $R/$TAG/${TAG}_truth_t${T}.point.bed ${TAG}t${T} $R/${TAG}t${T}
done
$K/score_vs_kinnex.sh $R/$TAG/${TAG}_truth_t5_plus_decoy.point.bed ${TAG}t5ip $R/${TAG}t5ip
echo "$(date +'%F %T') decomposition"
$WD/scripts/benchmark_tools/kinnex_decompose.sh $TAG $R/$TAG $K/scratch/${TAG}_bt/pos_counts.tsv
touch $K/PRIMARY_CHAIN.DONE.ok
echo "$(date +'%F %T') PRIMARY CHAIN DONE"
