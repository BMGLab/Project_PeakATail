#!/bin/bash
# Publishes the Stage 0 fix branches as four PRs. REVIEW THE BODIES FIRST, then run yourself.
# Posts under whichever account `gh auth status` shows.
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"
T=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
footer=$'\n\n---\n*Prepared by Claude (AI assistant) for the Stage-0 fix campaign (PI-approved); adversarially verified before staging. gh account per `gh auth status`.*\n\n🤖 Generated with [Claude Code](https://claude.com/claude-code)'

# fast-forward the original read-path branch onto its repaired head
if [ "$(git -C "$T" branch --show-current)" = "fix/cellranger-input-compat" ]; then git -C "$T" merge --ff-only fix/cellranger-input-compat-repair; else git -C "$T" branch -f fix/cellranger-input-compat fix/cellranger-input-compat-repair; fi

declare -A HEAD=( [matrix]=fix/matrix-pas-id-keying [switch]=fix/switch-length-strand-counts
                  [readpath]=fix/cellranger-input-compat [ci]=chore/ci-and-fixture )
declare -A TITLE=(
 [matrix]="fix(matrix): key count-matrix rows by PAS ID, not list position (99.87% of PAS carried another PAS's counts)"
 [switch]="fix(switch-length): strand-aware fallback rank + PDUI on raw counts, not TF-IDF (closes the #67 blocker)"
 [readpath]="fix(read): stock CellRanger BAM compatibility — CB suffix, unmapped reads, bijective RG sanitiser"
 [ci]="ci: pytest gate + CellRanger regression fixture + CITATION.cff (stacked on the read-path PR)"
)
for k in matrix switch readpath ci; do
  git -C "$T" push -u origin "${HEAD[$k]}"
  gh pr create -R $R -B develop -H "${HEAD[$k]}" -t "${TITLE[$k]}" \
     -b "$(cat pr_stage0_$k.md)$footer"
  gh pr edit -R $R "${HEAD[$k]}" --add-reviewer TRextabat || true
done
gh issue comment -R $R 67 -b "Stage-0 fixes are up as PRs: the switch-length PR closes this issue's blocker. Please HOLD the SWITCH_CELLTYPE rerun until it merges — rerunning today returns inverted PDUI computed on TF-IDF weights (1667/1667 fallback rows inverted, 10.6% of dPDUI signs flip). Half-day unblock path in the PR body.$footer"
echo "STAGE 0 PUBLISHED"
