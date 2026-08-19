#!/usr/bin/env bash
# kinnex_peakatail_selfrank.sh <tag> <kinnex_truth_dir>
# Does PeakATail's OWN confidence signal (total-UMI rank, TIER label, peak width) predict
# whether a call sits on a real long-read 3' end?  This is the orthogonal-evidence replication
# of manuscript/09 point 5 ("the TIER label gives false confidence").
set -euo pipefail
export LC_ALL=C
TAG=$1; G=$2
WD=/mnt/ssd1/Projects/PeakATail_wd
P=$WD/results/benchmark_tools/pbmc_10k_v3/peakatail
GEN=$WD/data/references/chrom.sizes.nochr.filt
S=/mnt/ssd0/emaout/peakatail_benchmark/kinnex/scratch/${TAG}_selfrank; mkdir -p $S
awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN $G/${TAG}_truth_t5.point.bed | sort -k1,1 -k2,2n > $S/truth5.bed
awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN $G/${TAG}_decoy.point.bed    | sort -k1,1 -k2,2n > $S/decoy.bed
OUT=$G/${TAG}_peakatail_selfrank.tsv
{ echo "# Does PeakATail's own ranking predict long-read support? (100 bp, strand-matched)"
  echo "# truth = Kinnex $TAG >=5 UMI non-internally-primed peaks; decoy = >=5 UMI internal-priming peaks"
  echo -e "subset\tn_calls\tgenuine_pas\tip_artifact\tno_lr_support"; } > $OUT
sup(){ paste <(bedtools closest -s -d -t first -a "$1" -b $S/truth5.bed 2>/dev/null|awk -F'\t' '{print $NF}') \
             <(bedtools closest -s -d -t first -a "$1" -b $S/decoy.bed 2>/dev/null|awk -F'\t' '{print $NF}') \
   | awk -v L="$2" -v OFS='\t' '{n++;g=($1>=0&&$1<=100);i=($2>=0&&$2<=100); if(g)a++; else if(i)b++; else c++}
       END{printf "%s\t%d\t%.6f\t%.6f\t%.6f\n", L,n,a/n,b/n,c/n}' >> $OUT; }
mk(){ grep -v '^#' "$1" | awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN - \
  | awk -F'\t' -v OFS='\t' '{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,NR,0,$6}' | sort -k1,1 -k2,2n; }
mk $P/pas.bed > $S/q.bed; sup $S/q.bed "all_calls"
for N in 22000 36000 106000; do mk $P/topN/pas_top${N}.bed > $S/q.bed; sup $S/q.bed "top${N}_by_umi"; done
for T in TIER_1 TIER_2; do
  awk -F'\t' -v OFS='\t' -v t=$T '$9==t{if($7=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,NR,0,$7}' $P/run/annotatedpas.bed \
   | awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN - | sort -k1,1 -k2,2n > $S/q.bed; sup $S/q.bed "$T"; done
for W in "0 50" "50 150" "150 400" "400 100000"; do set -- $W
  awk -F'\t' -v OFS='\t' -v lo=$1 -v hi=$2 '{w=$3-$2; if(w>=lo&&w<hi){if($7=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,NR,0,$7}}' $P/run/annotatedpas.bed \
   | awk -F'\t' 'NR==FNR{ok[$1]=1;next} ok[$1]' $GEN - | sort -k1,1 -k2,2n > $S/q.bed; sup $S/q.bed "peak_width_${1}_${2}bp"; done
echo "wrote $OUT"; cat $OUT
