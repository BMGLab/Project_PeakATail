#!/usr/bin/env bash
# Kinnex PBMC 10x3p (v3.1) primary dedup FLNC BAM -> GRCh38 (Ensembl no-chr) splice alignment.
# Tag pass-through: CB (cell barcode, whitelist orientation AS-IS), XM (UMI), rc (real-cell flag).
set -euo pipefail
export LC_ALL=C
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
IN=/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/kinnex_pbmc_10x3p/primary_10x3p/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam
REF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
OUT=$K/kinnex_primary.dedup.GRCh38.bam
rm -f $K/ALIGN_PRIMARY.DONE.ok $K/ALIGN_PRIMARY.FAILED.err
trap 'echo "FAILED rc=$?" > $K/ALIGN_PRIMARY.FAILED.err' ERR
date +'%F %T start' > $K/logs/align_primary.progress
samtools fastq -@ 3 -T CB,XM,rc "$IN" 2> $K/logs/fastq.err \
 | minimap2 -ax splice:hq -uf -y --secondary=no -t 8 "$REF" - 2> $K/logs/minimap2.err \
 | samtools sort -@ 4 -m 6G -T $K/scratch/sortp -o "$OUT" - 2> $K/logs/sort.err
samtools index -@ 4 "$OUT"
samtools idxstats "$OUT" > $K/logs/primary_idxstats.txt
samtools flagstat -@ 4 "$OUT" > $K/logs/primary_flagstat.txt
date +'%F %T done' >> $K/logs/align_primary.progress
touch $K/ALIGN_PRIMARY.DONE.ok
