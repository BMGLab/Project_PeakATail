#!/usr/bin/env bash
# Secondary (GEM-X 10x 3' v4) Kinnex BAM is ALREADY genome-aligned (pbmm2, hg38 GENCODE v39,
# chr-prefixed) -> strip to Ensembl no-chr on the fly (STRIP=1). Cross-chemistry replicate.
set -euo pipefail
export LC_ALL=C
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
BAM=/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/kinnex_pbmc_10x3p/secondary_gemx3p/scisoseq.mapped.bam
rm -f $K/EXTRACT_SECONDARY.DONE.ok $K/EXTRACT_SECONDARY.FAILED.err
trap 'echo "FAILED rc=$?" > $K/EXTRACT_SECONDARY.FAILED.err' ERR
samtools view -@ 3 -F 0x904 -q 1 "$BAM" \
 | gawk -v WL=$K/work/secondary_realcells.txt -v STRIP=1 -v CLIP=30 -f $K/termini.awk \
   > $K/work/secondary_termini.tsv 2> $K/logs/secondary_termini.stats
wc -l < $K/work/secondary_termini.tsv > $K/logs/secondary_termini.n
touch $K/EXTRACT_SECONDARY.DONE.ok
