#!/bin/bash
# Follow-up posts after the first publish round. Review, then run yourself.
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"
footer=$'\n\n---\n*Filed by Claude (AI assistant) on behalf of the PeakATail manuscript effort (Ebru); posted via the gh account shown by `gh auth status`.*'

# Issue 7: FDR calibration failure
gh issue create -R $R --label manuscript \
  -t "switch diff fisher is not FDR-calibrated: 100% of permutation-null runs report q<0.05 hits (mean 42/run); D4 cells-mode is not CLI-exposed" \
  -b "$(cat issue7_fdr_calibration.md)$footer"

# Comment on #67: rerun command should not rely on fisher defaults
gh issue comment -R $R 67 -b "Heads-up before the rerun: permutation-null calibration (new issue above) shows \`--strategy fisher\` q-values are anti-conservative to the point that null label shuffles yield as many q<0.05 hits as real labels (mean 42 vs 40). Suggest the fixed SWITCH_CELLTYPE either adds an \`nb_pairwise\` arm or waits for D4 \`count_mode=cells\` CLI exposure; fisher output should be treated as a ranking screen only.$footer"

# Comment on #69: ip-filter live verification verdict
gh issue comment -R $R 69 -b "Live verification done (small-BAM, three runs: no flags / annotate / filter — full write-up in PeakATail_wd/manuscript/08_ipfilter_live_check.md): the internal-priming filter IS wired and works — annotate mode writes True/False into the internal_priming column (3.8% flagged on the test BAM), filter mode drops exactly the flagged PAS (strict subset). And the default is genuinely ip-OFF: the no-flag run executes no ip code and is byte-identical to annotate's coordinates. So lg_ip_off duplicating lg_annotate is expected behavior — the sweep grid just needs relabeling, and the stale 'currently no-op' rows in docs/cli/run.md should be deleted.$footer"
echo "FOLLOW-UPS POSTED"
