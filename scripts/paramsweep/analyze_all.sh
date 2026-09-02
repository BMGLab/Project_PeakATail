#!/usr/bin/env bash
# scripts/paramsweep/analyze_all.sh -- every per-lead analysis, after the arms have run.
set -uo pipefail
source "$(dirname "$0")/common.sh"
A=$WD/results/paramsweep
mkdir -p $A/leadA $A/leadB $A/leadC $A/seam
PY=$WD/tools/PeakATail/.venv/bin/python
arms() { { cat $WD/scripts/paramsweep/GRID.tsv $WD/scripts/paramsweep/GRID_extra.tsv 2>/dev/null | grep -v '^#'; \
           printf 'pbmc_base\thuman\tbase\t(v2 default)\nm1_base\tmouse\tbase\t(v2 default)\n'; } ; }

# ---- LEAD A: the candidate ceiling, per arm -------------------------------
printf 'name\tspecies\tn_det_atlas\tn_cand\thit_by_cand\tceiling_R@100\thit_by_default\tdefault_R@100\tn_miss\tn_peakless\tpeakless_of_miss\treachable_of_miss\n' > $A/leadA/ceiling.tsv
arms | while IFS=$'\t' read -r label species block flags; do
  d=$OUTROOT/$label; [ -f "$d/DONE.ok" ] || continue
  case "$block" in AD|base) ;; *) continue;; esac
  C="19,21"; [ "$species" = mouse ] && C="18,19"
  [ -s "$d/score/candidates.bed" ] || continue
  $PY $WD/scripts/paramsweep/ceiling.py --cand "$d/score/candidates.bed" \
      --default-arm "$d/pas_tier1_ge2mol.bed" --name "$label" --species "$species" --contigs "$C" >> $A/leadA/ceiling.tsv
done

# ---- LEAD C: resolution / per-PAS quantification --------------------------
printf 'name\tn_calls\tmatched@25\tdistinct_atlas@25\tcalls_per_atlas_site\tmatched@100\trecall_sites@100\tmany_to_one@100\tmedian_gap\tp10_gap\tfrac_gap_lt100\n' > $A/leadC/resolution.tsv
arms | while IFS=$'\t' read -r label species block flags; do
  d=$OUTROOT/$label; [ -f "$d/DONE.ok" ] || continue
  case "$block" in C|base|I) ;; *) continue;; esac
  C="19,21"; [ "$species" = mouse ] && C="18,19"
  $PY $WD/scripts/paramsweep/resolution.py "$d/pas_tier1_ge2mol.bed" --name "${label}__default" --species "$species" --contigs "$C" >> $A/leadC/resolution.tsv 2>/dev/null
done

# ---- LEAD B: post-hoc offset sweeps + signed-distance histograms ----------
for a in pbmc_base m1_base; do
  sp=human; C="19,21"; [ "$a" = m1_base ] && { sp=mouse; C="18,19"; }
  [ -f "$OUTROOT/$a/DONE.ok" ] || continue
  bash $WD/scripts/paramsweep/offset_precision_sweep.sh "$OUTROOT/$a/pas_tier2.bed" "${a}_tier2" $sp "$C" $A/leadB > /dev/null
  bash $WD/scripts/paramsweep/offset_precision_sweep.sh "$OUTROOT/$a/pas_tier1_ge2mol.bed" "${a}_default" $sp "$C" $A/leadB > /dev/null
done
echo "analyses done -> $A"
