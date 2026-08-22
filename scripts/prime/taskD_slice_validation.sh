#!/usr/bin/env bash
# scripts/prime/taskD_slice_validation.sh -- TASK D's slice-level validation.
#
# Three arms on the PBMC chr19+21 dev slice:
#   1. taskD_slice_compat    the FULL v2 pin (--pas-features off, and
#      --pas-score none, which is the default).  Must reproduce the v2 reference
#      run byte-for-byte -- this branch's cardinal rule, re-validated on real
#      data.  NOTE the pin needs --pas-features off as well: TASK C made `on`
#      the branch default, so a bare no-flag run legitimately writes a wider
#      sidecar.  taskD_slice_default (no flags at all) is kept as the arm that
#      shows exactly that, and nothing else, is what differs.
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
"$RS" --label taskD_slice_compat  --threads 8 -- --pas-features off
"$RS" --label taskD_slice_cal     --threads 8 -- --genome-fasta "$HUM_FA" --pas-score calibrated
"$RS" --label taskD_slice_sel     --threads 8 -- --genome-fasta "$HUM_FA" \
      --ip-filter --ip-filter-mode filter --pas-score select
set +x

# compare the caller's OUTPUT TREE, not the harness's bookkeeping around it
"$PY" "$(dirname "$0")/identity_check.py" "$OUTROOT/slice_v2_ref/run" "$OUTROOT/taskD_slice_compat/run" \
  -o "$ARTROOT/identity_taskD_compat_vs_v2.tsv" | tee "$ARTROOT/identity_taskD_compat_vs_v2.txt"
"$PY" "$(dirname "$0")/identity_check.py" "$OUTROOT/slice_v2_ref/run" "$OUTROOT/taskD_slice_default/run" \
  -o "$ARTROOT/identity_taskD_default_vs_v2.tsv" | tee "$ARTROOT/identity_taskD_default_vs_v2.txt"
echo "TASK D slice validation finished $(date -Is)"
