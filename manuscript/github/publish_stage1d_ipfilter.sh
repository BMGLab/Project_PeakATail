#!/bin/bash
# Publish the Stage 1d-A PR (IP-filter minus-strand fix). Review, then run:  bash manuscript/github/publish_stage1d_ipfilter.sh
# Requires the Stage 1b+1c PR to exist first (publish_stage1bc.sh) — this branch is based on feat/polya-evidence.
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"; T=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
footer=$'\n\n---\n*Prepared by Claude (AI assistant) for the PeakATail manuscript effort (PI: Ebru Kocakaya); adversarially verified (FIXED → ready) before staging; posted via the gh account shown by `gh auth status`.*\n\n🤖 Generated with [Claude Code](https://claude.com/claude-code)'
git -C "$T" push -u origin fix/ip-filter-strand:fix/ip-filter-strand
PR=$(gh pr create -R $R -B develop -H fix/ip-filter-strand \
  -t "fix(ip-filter): test the downstream side of the cleavage site on the '-' strand (Stage 1d-A)" \
  -b "$(cat pr_stage1d_ipfilter_strand.md)$footer")
gh pr edit -R $R fix/ip-filter-strand --add-reviewer TRextabat || true
echo "- PR Stage 1d-A (IP-filter strand): $PR" >> PUBLISHED.md
echo "STAGE 1d-A PUBLISHED: $PR"
