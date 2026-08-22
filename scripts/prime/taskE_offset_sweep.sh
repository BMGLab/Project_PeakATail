#!/usr/bin/env bash
# scripts/prime/taskE_offset_sweep.sh <pas.bed> <label> <species> <contigs> <outdir>
#
# TASK E item 1: how a signed cleavage-offset shift moves precision at every
# window, against BOTH truths.  Reports a TSV -- one row per shift -- so the
# argmax and the flatness at W >= 25 are readable rather than asserted.
#
# The shift is applied in TRANSCRIPT orientation (+ = downstream): on '+' the
# point moves to a higher coordinate, on '-' to a lower one -- the same
# convention ema/countmatrix/cleavage_offset.py::shift_cleavage_point uses.
#
# Precision here is computed with `bedtools closest` directly rather than with
# score_tool.py, because a 17-point sweep x 3 arms is 51 scorer runs; the
# headline arms are ALWAYS re-scored with score_tool.py, and the two agree
# (see results/prime/taskE/TASK_E.md).
set -uo pipefail
source "$(dirname "$0")/common.sh"

BED=${1:?usage: taskE_offset_sweep.sh <pas.bed> <label> <species> <contigs> <outdir>}
LABEL=${2:?}; SPECIES=${3:?}; CONTIGS=${4:?}; OUT=${5:?}
mkdir -p "$OUT"
W=$(mktemp -d -p "$OUT" .sweep.XXXXXX)

keep() { awk -F'\t' -v C="$CONTIGS" 'BEGIN{n=split(C,a,",");for(i=1;i<=n;i++)k[a[i]]=1}($1 in k)' "$1" \
         | sort -k1,1 -k2,2n > "$2"; }

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
  keep "$KINNEX/x3p_truth_t5.point.bed"  "$W/kin5.bed"
  keep "$KINNEX/x3p_truth_t20.point.bed" "$W/kin20.bed"
  HAVE_KIN=1
fi

TSV=$OUT/offset_sweep_$LABEL.tsv
{
  printf 'shift\tn\tatlas_P1\tatlas_P5\tatlas_P10\tatlas_P25\tatlas_P100'
  [ "$HAVE_KIN" = 1 ] && printf '\tkin5_P1\tkin5_P10\tkin5_P25\tkin20_P25'
  printf '\n'
} > "$TSV"

N=$(wc -l < "$BED")
for s in $(seq -8 1 8); do
  awk -F'\t' -v S="$s" 'BEGIN{OFS="\t"}
    { if($6=="+"){b=$2+S}else{b=$2-S}; if(b<0)b=0; print $1,b,b+1,$4,$5,$6 }' "$BED" \
    | sort -k1,1 -k2,2n > "$W/q.bed"
  A=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/atlas.bed" 2>/dev/null \
      | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=1)p1++;if(d>=0&&d<=5)p5++;if(d>=0&&d<=10)p10++;
                     if(d>=0&&d<=25)p25++;if(d>=0&&d<=100)p100++}
                    END{printf "%.4f\t%.4f\t%.4f\t%.4f\t%.4f",p1/n,p5/n,p10/n,p25/n,p100/n}')
  ROW=$(printf '%s\t%s\t%s' "$s" "$N" "$A")
  if [ "$HAVE_KIN" = 1 ]; then
    K5=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/kin5.bed" 2>/dev/null \
        | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=1)p1++;if(d>=0&&d<=10)p10++;if(d>=0&&d<=25)p25++}
                      END{printf "%.4f\t%.4f\t%.4f",p1/n,p10/n,p25/n}')
    K20=$(bedtools closest -s -d -t first -a "$W/q.bed" -b "$W/kin20.bed" 2>/dev/null \
        | awk -F'\t' '{d=$NF;n++;if(d>=0&&d<=25)p25++}END{printf "%.4f",p25/n}')
    ROW="$ROW	$K5	$K20"
  fi
  printf '%s\n' "$ROW" >> "$TSV"
done
rm -rf "$W"
column -t "$TSV"
echo "sweep -> $TSV"
