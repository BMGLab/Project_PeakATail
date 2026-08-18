#!/bin/bash
# Publishes the prepared PeakATail collaboration items to GitHub.
# RUN THIS YOURSELF after reviewing the bodies in this directory.
# NOTE: posts under whichever account `gh auth status` shows (currently @yasinkaymaz on this box).
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"

gh label create manuscript -R $R -c '0072B2' -d 'blocks or supports the manuscript' 2>/dev/null || true

footer=$'\n\n---\n*Filed by Claude (AI assistant) on behalf of the PeakATail manuscript effort (Ebru), from the BioLab workstation; posted via the gh account shown above.*'

declare -A titles=(
 [1]="Sweep pipeline: SWITCH_CELLTYPE interpolates null pasbed — all 24 per-celltype switch tasks failed"
 [2]="ema benchmark: precision vs PolyASite is inflated by reference density and interval matching — needs point-mode + strand-matched scoring and a shuffled-null baseline"
 [3]="Grid arm lg_ip_off produced output byte-identical to lg_annotate — internal-priming contrast incomplete"
 [4]="Release engineering for publication: CI running pytest, PyPI/bioconda, CITATION.cff, Zenodo DOI"
 [5]="Docs inconsistencies to resolve before Methods is written"
)
declare -A files=(
 [1]=issue1_switch_celltype_bug.md [2]=issue2_benchmark_defensibility.md
 [3]=issue3_lg_ip_off_duplicate.md [4]=issue4_release_engineering.md
 [5]=issue5_docs_inconsistencies.md
)
declare -A urls
for i in 1 2 3 4 5; do
  body="$(cat "${files[$i]}")$footer"
  urls[$i]=$(gh issue create -R $R --label manuscript -t "${titles[$i]}" -b "$body")
  echo "issue $i -> ${urls[$i]}"
done
BENCH_NUM=$(basename "${urls[2]}")

# --- PR: point-mode strand-matched benchmark (branch prepared locally by Claude) ---
TOOL=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail-prwork
if [ -d "$TOOL" ] && git -C "$TOOL" log --oneline -1 >/dev/null 2>&1; then
  git -C "$TOOL" push -u origin feat/benchmark-point-strand
  PRBODY_FILE=pr_benchmark_body.md
  sed "s/__BENCH_ISSUE__/${BENCH_NUM}/" "$PRBODY_FILE" > /tmp/prbody.$$
  gh pr create -R $R -B develop -H feat/benchmark-point-strand \
     -t "feat(benchmark): point-mode, strand-matched scoring with shuffled-null baseline" \
     --body-file /tmp/prbody.$$ 
  gh pr edit -R $R feat/benchmark-point-strand --add-reviewer TRextabat || true
  gh issue comment -R $R "$BENCH_NUM" -b "PR implementing this: see the linked feat/benchmark-point-strand pull request."
  rm -f /tmp/prbody.$$
else
  echo "PR branch not found at $TOOL — skipping PR step."
fi
echo "ALL DONE"
