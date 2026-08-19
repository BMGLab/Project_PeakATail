#!/usr/bin/env bash
# build_truth2.sh <tag> <scratchdir_with_termini.sorted.tsv+uniq_ip.tsv> <outdir>
# Greedy-peak version of Step 4 (see peak_call.py for why single-linkage was replaced).
set -euo pipefail
export LC_ALL=C
TAG=$1; SCR=$2; OUT=$3
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
REF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
mkdir -p "$OUT"
REP=$OUT/${TAG}_filter_report.txt

echo "[A] per-position UMI counts (lockstep with internal-priming calls)" | tee -a "$REP"
gawk -F'\t' -v OFS='\t' -v IPF="$SCR/uniq_ip.tsv" '
  function flush(){ if(havep){ if((getline il < IPF)<=0){print "IPF underrun">"/dev/stderr"; exit 3}
      split(il,ia,"\t")
      if(ia[1]!=pc||ia[2]+0!=pp+0||ia[3]!=ps){print "IPF desync "pc","pp","ps" vs "il >"/dev/stderr"; exit 3}
      print pc,pp,ps,cnt,ia[4] } }
  { if($1!=pc||$3!=ps||$2!=pp){ flush(); pc=$1; pp=$2; ps=$3; cnt=0; havep=1 } cnt++ }
  END{ flush() }' "$SCR/termini.sorted.tsv" > "$SCR/pos_counts.tsv"
echo "  unique positions = $(wc -l < "$SCR/pos_counts.tsv")" | tee -a "$REP"

for MODE in truth decoy; do
  echo "[B:$MODE] greedy peak call (+/-25 nt exclusion, >=5 UMI)" | tee -a "$REP"
  python3 $K/peak_call.py "$SCR/pos_counts.tsv" $MODE \
      "$SCR/${MODE}.peaks.bed" "$SCR/${MODE}.assign.tsv" --win 25 --min 5 2>&1 | tee -a "$REP"
  sort -k1,1 -k2,2n "$SCR/${MODE}.peaks.bed" > "$OUT/${TAG}_${MODE}_pas.bed"
done

echo "[C] truth PAS x cell-barcode counts" | tee -a "$REP"
gawk -F'\t' -v OFS='\t' -v AF="$SCR/truth.assign.tsv" '
  { if($1!=pc||$3!=ps||$2!=pp){ if((getline al < AF)<=0){print "assign underrun">"/dev/stderr"; exit 3}
      split(al,aa,"\t")
      if(aa[1]!=$1||aa[2]+0!=$2+0||aa[3]!=$3){print "assign desync">"/dev/stderr"; exit 3}
      cur=aa[4]; pc=$1; pp=$2; ps=$3 }
    if(cur!=".") print cur, $4 }' "$SCR/termini.sorted.tsv" \
 | sort -S 24G --parallel=8 -T "$SCR" | uniq -c \
 | gawk -v OFS='\t' '{print $2,$3,$1}' > "$OUT/${TAG}_truth_pas_by_cb.triplets.tsv"
echo "  nonzero PAS x CB entries = $(wc -l < "$OUT/${TAG}_truth_pas_by_cb.triplets.tsv")" | tee -a "$REP"
cut -f1 "$OUT/${TAG}_truth_pas_by_cb.triplets.tsv" | uniq -c | gawk -v OFS='\t' '{print $2,$1}' \
  | sort -k1,1 > "$SCR/ncells.tsv"

echo "[D] point BEDs at a range of molecular-support thresholds" | tee -a "$REP"
sort -k4,4 "$OUT/${TAG}_truth_pas.bed" > "$SCR/truth.byid.bed"
join -1 4 -2 1 -t$'\t' -a1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2 -e 0 \
  "$SCR/truth.byid.bed" "$SCR/ncells.tsv" > "$SCR/truth.withcells.bed"
mv "$SCR/truth.withcells.bed" "$SCR/tmp" && sort -k1,1 -k2,2n "$SCR/tmp" > "$OUT/${TAG}_truth_pas.bed"
for T in 5 20 100 500; do
  gawk -F'\t' -v OFS='\t' -v T=$T '$5>=T{print $1,$7,$7+1,$4,$5,$6}' "$OUT/${TAG}_truth_pas.bed" \
   | sort -k1,1 -k2,2n > "$OUT/${TAG}_truth_t${T}.point.bed"
  echo "  truth >=${T} UMI : $(wc -l < "$OUT/${TAG}_truth_t${T}.point.bed") PAS" | tee -a "$REP"
done
gawk -F'\t' -v OFS='\t' '{print $1,$7,$7+1,$4,$5,$6}' "$OUT/${TAG}_decoy_pas.bed" \
 | sort -k1,1 -k2,2n > "$OUT/${TAG}_decoy.point.bed"
cat "$OUT/${TAG}_truth_t5.point.bed" "$OUT/${TAG}_decoy.point.bed" | sort -k1,1 -k2,2n \
 > "$OUT/${TAG}_truth_t5_plus_decoy.point.bed"
echo "  decoy (internal-priming) PAS : $(wc -l < "$OUT/${TAG}_decoy.point.bed")" | tee -a "$REP"
echo "DONE $TAG (greedy)" | tee -a "$REP"
