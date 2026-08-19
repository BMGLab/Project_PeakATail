#!/usr/bin/env bash
# truth_qc.sh <truth_point_bed> <tag> <report>
# Acceptance checks from kinnex_truth_plan.md section 3, INDEPENDENT of how the truth was built.
set -euo pipefail
export LC_ALL=C
T=$1; TAG=$2; REP=$3
WD=/mnt/ssd1/Projects/PeakATail_wd
PAS2=$WD/data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6
TES=$WD/data/references/atlases/tes.protein_coding.GRCh38.99.bed6
GB=$WD/data/references/gene_end.bed
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex; S=$K/scratch/${TAG}_qc; mkdir -p $S
{
echo "== truth QC: $TAG =="
N=$(wc -l < "$T"); echo "truth PAS (point) n = $N"
for ref in "$PAS2:PolyASite2.0" "$TES:protein-coding TES"; do
  f=${ref%%:*}; nm=${ref##*:}
  sort -k1,1 -k2,2n "$f" > $S/ref.bed
  bedtools closest -s -d -t first -a "$T" -b $S/ref.bed 2>/dev/null \
  | gawk -F'\t' -v NM="$nm" -v N=$N '{d=$NF; if(d>=0){ if(d<=25)a++; if(d<=50)b++; if(d<=100)c++ }}
      END{printf("  within 25/50/100 bp of %s: %.4f / %.4f / %.4f\n", NM, a/N,b/N,c/N)}'
done
# strand sanity: truth PAS antisense to an overlapping annotated gene body
sort -k1,1 -k2,2n "$GB" > $S/gb.bed
S_SENSE=$(bedtools intersect -a "$T" -b $S/gb.bed -s -u | wc -l)
S_ANTI=$(bedtools intersect -a "$T" -b $S/gb.bed -S -u | wc -l)
gawk -v s=$S_SENSE -v a=$S_ANTI -v n=$N 'BEGIN{printf("  in a sense gene body: %.4f ; ONLY antisense-overlapping: n/a (sense=%d antisense=%d)\n", s/n, s, a)}'
gawk -v a=$S_ANTI -v s=$S_SENSE 'BEGIN{ if(a+s>0) printf("  antisense fraction among gene-body-overlapping truth PAS: %.4f\n", a/(a+s)) }'
} | tee -a "$REP"
