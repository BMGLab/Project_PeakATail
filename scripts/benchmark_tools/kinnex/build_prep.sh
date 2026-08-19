#!/usr/bin/env bash
# build_prep.sh <termini.tsv> <scratchdir>   -- Steps 1-3 of the plan (sort, unique positions,
# internal-priming call). Shared by the greedy peak-calling truth build (build_truth2.sh).
set -euo pipefail
export LC_ALL=C
IN=$1; SCR=$2
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
REF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
mkdir -p "$SCR"
sort -S 24G --parallel=8 -T "$SCR" -t$'\t' -k1,1 -k3,3 -k2,2n "$IN" > "$SCR/termini.sorted.tsv"
echo "  termini_total = $(wc -l < "$SCR/termini.sorted.tsv")"
gawk -F'\t' -v OFS='\t' '{ k=$1 SUBSEP $3 SUBSEP $2; if(k==p) next; p=k; print $1,$2,$3 }' \
  "$SCR/termini.sorted.tsv" > "$SCR/uniq_key.tsv"
python3 $K/ip_flag.py "$REF" "$SCR/uniq_key.tsv" "$SCR/uniq_ip.tsv"
