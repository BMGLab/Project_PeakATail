#!/usr/bin/env bash
# scripts/paramsweep/identity_compute.sh -- per-arm inertness + compute table.
# An arm whose run/pasbed.bed is byte-identical to its species baseline changed
# NOTHING; that is the decisive test for whether a flag is wired into the
# default clip_seeded path at all.
set -uo pipefail
source "$(dirname "$0")/common.sh"
OUT=$WD/results/paramsweep/identity_and_compute.tsv
printf 'arm\tspecies\tflags\tmd5_pasbed\tidentical_to_baseline\tn_candidates\tn_pas\tn_tier1\tn_tier2\tn_default_arm\twall_s\tpeak_rss_gb\n' > $OUT
B_H=$(md5sum "$OUTROOT/pbmc_base/run/pasbed.bed" | cut -d' ' -f1)
B_M=$(md5sum "$OUTROOT/m1_base/run/pasbed.bed" | cut -d' ' -f1)
emit() {
  local label=$1 species=$2 flags=$3
  local d=$OUTROOT/$label
  [ -f "$d/DONE.ok" ] || return 0
  local m; m=$(md5sum "$d/run/pasbed.bed" | cut -d' ' -f1)
  local base=$B_H; [ "$species" = mouse ] && base=$B_M
  local same=no; [ "$m" = "$base" ] && same=YES
  local nc; nc=$(cat "$d"/run/peakcalling/*.pos.bed "$d"/run/peakcalling/*.neg.bed | wc -l)
  local w r
  w=$(awk -F': ' '/Elapsed \(wall clock\)/{split($NF,a,":");
        if(length(a)==3) printf "%.0f", a[1]*3600+a[2]*60+a[3]; else printf "%.0f", a[1]*60+a[2]}' "$d/runtime_mem.txt")
  r=$(awk '/Maximum resident set size/{printf "%.3f", $NF/1048576}' "$d/runtime_mem.txt")
  printf '%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n' \
    "$label" "$species" "$flags" "$m" "$same" "$nc" \
    "$(wc -l < "$d/pas.bed")" "$(wc -l < "$d/pas_tier1.bed")" \
    "$(wc -l < "$d/pas_tier2.bed")" "$(wc -l < "$d/pas_tier1_ge2mol.bed")" "$w" "$r" >> $OUT
}
emit pbmc_base human "(v2 default)"
emit m1_base mouse "(v2 default)"
for g in "$WD/scripts/paramsweep/GRID.tsv" "$WD/scripts/paramsweep/GRID_extra.tsv"; do
  [ -f "$g" ] || continue
  while IFS=$'\t' read -r label species block flags; do
    case "$label" in \#*|"") continue;; esac
    emit "$label" "$species" "$flags"
  done < "$g"
done
column -t -s$'\t' $OUT
echo "-> $OUT"
