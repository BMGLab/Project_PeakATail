#!/usr/bin/env bash
# scripts/prime/score_pas.sh <pas.bed> <label> [--species human|mouse] [--outdir DIR]
#                            [--slice-contigs 19,21] [--no-kinnex]
#
# Scores any PAS point BED with scripts/benchmark_tools/score_tool.py using the
# EXACT reference arguments of scripts/benchmark_tools/stage2_final_launch.sh, then
# prints the six headline numbers and (human only) the Kinnex long-read concordance
# at 25 bp.  Nothing here re-implements the metric: the P/R/F1 numbers are READ OUT
# of score_tool.py's own TSV, so they are row-for-row comparable with the v2
# manuscript arms.
#
#   P@k       = precision / real / atlas_full     / cutoff k   (curated atlas + TES)
#   R_det@100 = recall    / real / atlas_detected / cutoff 100 (detected-gene atlas)
#   F1        = f1        / real / atlas_detected / cutoff 100 (the Stage-2 gate metric)
#   Kinnex    = fraction of calls with an x3p long-read 3' end within 25 bp, same
#               strand (t5 and t20 truth), plus the internal-priming decoy rate
#
# --slice-contigs restricts EVERY reference (atlas, TES, chrom.sizes, gene bodies,
# detected atlas, Kinnex truth/decoy) to the named contigs, so a chr19+21 slice run
# is scored against slice-sized denominators instead of genome-wide ones.  Without
# it the recall denominator stays genome-wide (correct for full-BAM runs, and still
# valid for comparing two slice runs to each other).
set -uo pipefail
source "$(dirname "$0")/common.sh"

BED=${1:?usage: score_pas.sh <pas.bed> <label> [opts]}; LABEL=${2:?label required}; shift 2
SPECIES=human; OUTDIR=""; CONTIGS=""; DO_KIN=1
while [ $# -gt 0 ]; do
  case "$1" in
    --species) SPECIES=$2; shift 2;;
    --outdir)  OUTDIR=$2; shift 2;;
    --slice-contigs) CONTIGS=$2; shift 2;;
    --no-kinnex) DO_KIN=0; shift;;
    *) echo "unknown option $1"; exit 2;;
  esac
done
[ -s "$BED" ] || { echo "empty or missing: $BED"; exit 1; }
OUTDIR=${OUTDIR:-$(dirname "$BED")}
mkdir -p "$OUTDIR"

REFS=()
if [ "$SPECIES" = mouse ]; then REFS=("${MOUSE_REFS[@]}"); else REFS=("${HUMAN_DET[@]}"); fi

# ---- optional slice restriction of every reference file -------------------
if [ -n "$CONTIGS" ]; then
  W=$OUTDIR/.slicerefs_$(echo "$CONTIGS" | tr ',' '_'); mkdir -p "$W"
  keep() { awk -F'\t' -v C="$CONTIGS" 'BEGIN{n=split(C,a,",");for(i=1;i<=n;i++)k[a[i]]=1}($1 in k)' "$1" > "$2"; }
  if [ "$SPECIES" = mouse ]; then
    A=$R/atlases/polyasite2.GRCm38.96.rep_sites.bed6; T=$R/atlases/tes.protein_coding.GRCm38.102.bed6
    S=$R/mouse/chrom.sizes.filt; GB=$G/shared_refs/genebodies.merged.bed; DA=$G/shared_refs/pas2.in_detected_genes.bed
  else
    A=$R/atlases/polyasite2.GRCh38.96.rep_sites.bed6; T=$R/atlases/tes.protein_coding.GRCh38.99.bed6
    S=$SIZES; GB=$WD/results/benchmark_tools/shared_refs/genebodies.merged.bed
    DA=$WD/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed
  fi
  keep "$A" "$W/atlas.bed"; keep "$T" "$W/tes.bed"; keep "$S" "$W/chrom.sizes"
  keep "$GB" "$W/genebodies.bed"; keep "$DA" "$W/detected_atlas.bed"
  REFS=(--atlas "$W/atlas.bed" --tes "$W/tes.bed" --genome "$W/chrom.sizes"
        --genebodies "$W/genebodies.bed" --detected-atlas "$W/detected_atlas.bed")
  echo "[refs] slice-restricted to $CONTIGS: atlas $(wc -l < "$W/atlas.bed") tes $(wc -l < "$W/tes.bed") det_atlas $(wc -l < "$W/detected_atlas.bed") genebodies $(wc -l < "$W/genebodies.bed")"
fi

"$PY" "$SCORE" "$BED" "$LABEL" --outdir "$OUTDIR" "${REFS[@]}" > "$OUTDIR/score_$LABEL.log" 2>&1 \
  || { echo "score_tool.py FAILED"; tail -20 "$OUTDIR/score_$LABEL.log"; exit 1; }
TSV=$OUTDIR/score_$LABEL.tsv
[ -s "$TSV" ] || { echo "no score TSV produced"; exit 1; }

# ---- Kinnex x3p concordance at 25 bp (human only) ------------------------
if [ "$DO_KIN" = 1 ] && [ "$SPECIES" = human ] && [ -d "$KINNEX" ]; then
  KW=$(mktemp -d -p "$OUTDIR" .kin.XXXXXX)
  CFILT=${CONTIGS:-}
  filt() { if [ -n "$CFILT" ]; then awk -F'\t' -v C="$CFILT" 'BEGIN{n=split(C,a,",");for(i=1;i<=n;i++)k[a[i]]=1}($1 in k)'; \
           else awk -F'\t' 'NR==FNR{ok[$1]=1;next}($1 in ok)' "$SIZES" -; fi; }
  grep -v '^#' "$BED" | filt | sort -k1,1 -k2,2n > "$KW/q.bed"
  NQ=$(wc -l < "$KW/q.bed")
  for tr in t5 t20 decoy; do
    src=$KINNEX/x3p_truth_$tr.point.bed; [ "$tr" = decoy ] && src=$KINNEX/x3p_decoy.point.bed
    filt < "$src" | sort -k1,1 -k2,2n > "$KW/$tr.bed"
    H=$(bedtools closest -s -d -t first -a "$KW/q.bed" -b "$KW/$tr.bed" 2>/dev/null \
        | awk -F'\t' '{d=$NF; if(d>=0 && d<=25) h++} END{print h+0}')
    eval "HIT_$tr=$H"
  done
  rm -rf "$KW"
fi

"$PY" - "$TSV" "$LABEL" "${HIT_t5:--1}" "${HIT_t20:--1}" "${HIT_decoy:--1}" "${NQ:--1}" <<'PYEOF'
import sys
import pandas as pd
tsv, label, h5, h20, hd, nq = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
d = pd.read_csv(tsv, sep="\t")
def g(panel, ref, cut):
    r = d[(d.panel == panel) & (d.series == "real") & (d.reference == ref) & (d.cutoff_bp == cut)]
    return float(r.value.iloc[0]) if len(r) else float("nan")
meta = d[d.panel == "meta"]
n = int(meta.n_matched.iloc[0]) if len(meta) else -1
raw = int(meta.n_query.iloc[0]) if len(meta) else -1
null100 = d[(d.panel == "precision") & (d.series == "null_genic") & (d.cutoff_bp == 100)].value.mean()
rows = [("n (scored / raw)", "%d / %d" % (n, raw)),
        ("P@10   (atlas_full)", "%.4f" % g("precision", "atlas_full", 10)),
        ("P@25   (atlas_full)", "%.4f" % g("precision", "atlas_full", 25)),
        ("P@100  (atlas_full)", "%.4f" % g("precision", "atlas_full", 100)),
        ("R_det@100 (detected)", "%.4f" % g("recall", "atlas_detected", 100)),
        ("F1@100 (detected)", "%.4f" % g("f1", "atlas_detected", 100)),
        ("null P@100 (3 seeds)", "%.4f" % null100)]
if nq > 0:
    rows += [("Kinnex t5  P@25", "%.4f (%d/%d)" % (h5 / nq, h5, nq)),
             ("Kinnex t20 P@25", "%.4f (%d/%d)" % (h20 / nq, h20, nq)),
             ("Kinnex decoy@25", "%.4f (%d/%d)" % (hd / nq, hd, nq))]
w = max(len(k) for k, _ in rows)
print("---- %s ----" % label)
for k, v in rows:
    print("%-*s  %s" % (w, k, v))
PYEOF
echo "score TSV: $TSV"
