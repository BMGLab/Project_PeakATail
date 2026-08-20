#!/bin/bash
# Publish Stage 1 after PRs #81-#84 are merged (or at least #83). Review, then run.
set -euo pipefail
cd "$(dirname "$0")"
T=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
footer=$'\n\n---\n*Prepared by Claude (AI assistant); adversarially verified (SOUND) before staging.*\n\n🤖 Generated with [Claude Code](https://claude.com/claude-code)'
git -C "$T" push -u origin feat/polya-evidence
gh pr create -R BMGLab/PeakATail -B develop -H feat/polya-evidence \
  -t "feat(caller): poly(A) read evidence + clip-seeded PAS calling — slice F1 0.177 -> 0.291" \
  -b "$(cat pr_stage1_polya.md)$footer"
gh pr edit -R BMGLab/PeakATail feat/polya-evidence --add-reviewer TRextabat || true
echo "STAGE 1 PUBLISHED"
