#!/usr/bin/env bash
# scripts/paramsweep/offset_precision_sweep.sh <pas.bed> <label> <species> <contigs> <outdir> [shifts]
#
# LEAD B: precision of a call set against BOTH truths as a function of a POST-HOC
# transcript-oriented shift, over the wide range the ~90-105 nt coverage-offset
# hypothesis lives in (the prime harness's taskE_offset_sweep.sh only covers
# -8..+8, which is the base-pair-resolution question, not this one).
# Post-hoc means POSITION ONLY: a real --cleavage-offset run also changes the
# call set, because the shift happens before the internal-priming veto and gene
# assignment.  Both are reported; this one isolates the geometry.
set -uo pipefail
source "$(dirname "$0")/common.sh"
BED=${1:?}; LABEL=${2:?}; SPECIES=${3:?}; CONTIGS=${4:?}; OUT=${5:?}
SHIFTS=${6:-"0 10 20 30 40 50 60 70 80 85 90 95 100 105 110 120 140 160 200"}
mkdir -p "$OUT"; W=$(mktemp -d -p "$OUT" .sweep.XXXXXX)
keep() { awk -F'\t' -v C="$CONTIGS" 'BEGIN{n=split(C,a,",");for(i=1;i<=n;i++)k[a[i]]=1}($1 in k)' "$1" | sort -k1,1 -k2,2n > "$2"; }
if [ "$SPECIES" = mouse ]; then
  keep "$R/atlases/polyasite2.GRCm38.96.rep_sites.bed6" "$W/a1.bed"
  keep "$R/atlases/tes.protein_coding.GRCm38.102.bed6"  "$W/a2.bed"
else
  keep "$R/atlases/polyasite2.GRCh38.96.rep_sites.bed6" "$W/a1.bed"
  keep "$R/atlases/tes.protein_coding.GRCh38.99.bed6"   "$W/a2.bed"
fi
sort -k1,1 -k2,2n "$W/a1.bed" "$W/a2.bed" > "$W/atlas.bed"
HAVE_KIN=0
if [ "$SPECIES" = human ] && [ -d "$KINNEX" ]; then
  keep "$KINNEX/x3p_truth_t5.point.bed" "$W/kin5.bed"; keep "$KINNEX/x3p_decoy.point.bed" "$W/dec.bed"; HAVE_KIN=1
fi
TSV=$OUT/offsetsweep_$LABEL.tsv
{ printf 'label\tshift\tn\tatlas_P10\tatlas_P25\tatlas_P50\tatlas_P100'
  [ "$HAVE_KIN" = 1 ] && printf '\tkin5_P10\tkin5_P25\tkin5_P100\tdecoy_P25'
  printf '\n'; } > "$TSV"
keep "$BED" "$W/src.bed"; N=$(wc -l < "$W/src.bed")
for s in $SHIFTS; do
  awk -F'\t' -v S="$s" 'BEGIN{OFS="\t"}{if($6=="+"){b=$2+S}else{b=$2-S}; if(b<0)b=0; print $1,b,b+1,$4,$5,$6}' "$W/src.bed" | sort -k1,1 -k2,2n > "$W/q.bed"
  A=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/atlas.bed" 2>/dev/null \
      | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=10)p10++;if(d>=0&&d<=25)p25++;if(d>=0&&d<=50)p50++;if(d>=0&&d<=100)p100++}
                    END{printf "%.4f\t%.4f\t%.4f\t%.4f",p10/n,p25/n,p50/n,p100/n}')
  ROW=$(printf '%s\t%s\t%s\t%s' "$LABEL" "$s" "$N" "$A")
  if [ "$HAVE_KIN" = 1 ]; then
    K=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/kin5.bed" 2>/dev/null \
        | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=10)p10++;if(d>=0&&d<=25)p25++;if(d>=0&&d<=100)p100++}END{printf "%.4f\t%.4f\t%.4f",p10/n,p25/n,p100/n}')
    D=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/dec.bed" 2>/dev/null \
        | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=25)p25++}END{printf "%.4f",p25/n}')
    ROW="$ROW	$K	$D"
  fi
  printf '%s\n' "$ROW" >> "$TSV"
done
rm -rf "$W"; column -t "$TSV"
