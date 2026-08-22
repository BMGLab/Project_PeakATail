#!/usr/bin/env bash
# scripts/prime/taskD_slice_validation.sh -- TASK D's slice-level validation.
#
# Three arms on the PBMC chr19+21 dev slice:
#   1. taskD_slice_default   branch defaults, no flags.  Must reproduce the v2
#      reference run byte-for-byte (`--pas-score none` is the default), which is
#      this branch's cardinal rule re-validated on real data.
#   2. taskD_slice_cal       --pas-score calibrated (+ --genome-fasta): the
#      sidecar gains pas_score and the BED must not move.
#   3. taskD_slice_sel       --ip-filter filter + --pas-score select: the arm the
#      tool itself would emit as the scored default, scored against the same
#      references as the v2 ipfilt arm.
set -uo pipefail
source "$(dirname "$0")/common.sh"
RS=$(dirname "$0")/run_slice.sh

set -x
"$RS" --label taskD_slice_default --threads 8
"$RS" --label taskD_slice_cal     --threads 8 -- --genome-fasta "$HUM_FA" --pas-score calibrated
"$RS" --label taskD_slice_sel     --threads 8 -- --genome-fasta "$HUM_FA" \
      --ip-filter --ip-filter-mode filter --pas-score select
set +x

"$PY" "$(dirname "$0")/identity_check.py" "$OUTROOT/slice_v2_ref" "$OUTROOT/taskD_slice_default" \
  -o "$ARTROOT/identity_taskD_default_vs_v2.tsv" | tee "$ARTROOT/identity_taskD_default_vs_v2.txt"
echo "TASK D slice validation finished $(date -Is)"
