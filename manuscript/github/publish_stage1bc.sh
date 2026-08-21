#!/bin/bash
# Publish the Stage 1b+1c PR and issues 9/10. Review, then run:  bash manuscript/github/publish_stage1bc.sh
# (Claude's attempt to run this was blocked by the permission classifier; the commands are the same as publish_stage1.sh.)
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"; T=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
footer=$'\n\n---\n*Prepared by Claude (AI assistant) for the PeakATail manuscript effort (PI: Ebru Kocakaya); adversarially verified before staging; posted via the gh account shown by `gh auth status`.*\n\n🤖 Generated with [Claude Code](https://claude.com/claude-code)'
# 1. push the branch ref (no checkout needed; tools/pa-polya may be on another branch)
git -C "$T" push -u origin feat/polya-evidence:feat/polya-evidence
# 2. PR feat/polya-evidence -> develop
PR=$(gh pr create -R $R -B develop -H feat/polya-evidence \
  -t "Stage 1b+1c: count clip_seeded tier-1 PAS from the cleavage site; BED col5 = UMI-deduplicated molecules; pas_support.tsv sidecar" \
  -b "$(cat pr_stage1bc_tier1_quant_umi.md)$footer")
gh pr edit -R $R feat/polya-evidence --add-reviewer TRextabat || true
# 3. issues referenced by the PR's follow-ups
I9=$(gh issue create -R $R --label manuscript \
  -t "switch diff: default marker pre-selection is a label double-dip (anti-conservative in all tests) and changes the Fisher denominator; only fisher --count-mode cells --marker-top-n 0 controls FDR" \
  -b "$(cat issue9_switch_marker_preselection.md)$footer")
I10=$(gh issue create -R $R --label bug \
  -t "Stage 1d: internal-priming filter tests the wrong side on the minus strand; 280 GiB peak RSS / ~1.4 cores on PBMC 10k; output bookkeeping (IP flag, zero-count PAS, unified bed col5)" \
  -b "$(cat issue10_stage1d_ipfilter_memory.md)$footer")
{ echo "# Published $(date '+%Y-%m-%d %H:%M')"; echo "- PR Stage 1b+1c: $PR"; echo "- Issue 9: $I9"; echo "- Issue 10: $I10"; } >> PUBLISHED.md
echo "STAGE 1b+1c PUBLISHED: $PR | $I9 | $I10"
