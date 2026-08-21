#!/bin/bash
# Stage 3 (Laughney) -- run EVERY sample in laughney_samples.tsv through laughney_per_sample.sh.
# Concurrency 3 x 8 threads = 24 threads (the cap while the Stage-2 final run occupies 56).
# Skips samples whose DONE.ok exists. Logs: $OUT_ROOT/logs/<sample>.log + run_all.log.
# *** NOT launched automatically. Check `uptime`, `free -g`, and that results/benchmark_tools/*final* are
#     finished or that 24 more threads fit, before starting:  nohup scripts/stage3/run_all_samples.sh &
set -uo pipefail
export LC_ALL=C
WD=/mnt/ssd1/Projects/PeakATail_wd
TABLE=$WD/scripts/stage3/laughney_samples.tsv
DRIVER=$WD/scripts/stage3/laughney_per_sample.sh
export STAGE3_OUT_ROOT=${STAGE3_OUT_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney}
CONC=${CONC:-3}; THREADS=${THREADS:-8}; N_PERM=${N_PERM:-5}
# per-sample policy knobs, passed through to laughney_per_sample.sh (verifier addition 2026-08-21):
#   STAGE3_STOP_AFTER=run|label|all  -- 'run' = peak-call only (label-independent); switch tests deferred until the
#                                       label gate / TIER_FILTER / PAS-id design are decided
#   TIER_FILTER=none|tier1|tier1_ge2, MIN_CELLS, LENGTH_TIMEOUT  -- see the driver header
export STAGE3_STOP_AFTER=${STAGE3_STOP_AFTER:-all} TIER_FILTER=${TIER_FILTER:-none} MIN_CELLS=${MIN_CELLS:-20} LENGTH_TIMEOUT=${LENGTH_TIMEOUT:-600}
case "$STAGE3_OUT_ROOT" in /mnt/ssd2*) echo "REFUSING: output root on /mnt/ssd2"; exit 2;; esac
mkdir -p "$STAGE3_OUT_ROOT/logs"
MASTER=$STAGE3_OUT_ROOT/logs/run_all.log
{
  echo "=== run_all start $(date) | conc $CONC x threads $THREADS = $((CONC*THREADS)) | n_perm $N_PERM | stop_after $STAGE3_STOP_AFTER | tier_filter $TIER_FILTER | min_cells $MIN_CELLS | length_timeout $LENGTH_TIMEOUT"
  echo "load: $(uptime) | free_gb: $(free -g | awk '/Mem:/{print $4}') | nproc $(nproc)"
  echo "ssd0: $(df -h /mnt/ssd0 | tail -1)"
} | tee -a "$MASTER"
PENDING=$(awk -F'\t' 'NR>1{print $1}' "$TABLE" | while read -r s; do
  if [[ -f $STAGE3_OUT_ROOT/$s/DONE.ok ]]; then echo "skip $s (DONE.ok)" >&2; else echo "$s"; fi; done)
echo "pending: $(echo "$PENDING" | wc -w) samples" | tee -a "$MASTER"
if [[ -z $PENDING ]]; then echo "nothing pending (every sample has DONE.ok)" | tee -a "$MASTER"; exit 0; fi
export DRIVER THREADS N_PERM
printf '%s\n' $PENDING | xargs -P "$CONC" -I{} bash -c \
  'echo "start {} $(date)"; "$DRIVER" {} "$THREADS" "$N_PERM" > "$STAGE3_OUT_ROOT/logs/{}.log" 2>&1; rc=$?; echo "end {} exit $rc $(date)"' \
  2>&1 | tee -a "$MASTER"
echo "=== run_all finished $(date)" | tee -a "$MASTER"
for s in $(awk -F'\t' 'NR>1{print $1}' "$TABLE"); do
  if [[ -f $STAGE3_OUT_ROOT/$s/DONE.ok ]]; then st=DONE; elif [[ -f $STAGE3_OUT_ROOT/$s/FAILED.err ]]; then st="FAILED($(cat $STAGE3_OUT_ROOT/$s/FAILED.err))"; elif ls $STAGE3_OUT_ROOT/$s/STOPPED_AFTER.* >/dev/null 2>&1; then st="STOPPED($(ls $STAGE3_OUT_ROOT/$s/STOPPED_AFTER.* | sed 's/.*STOPPED_AFTER\.//'))"; else st=PENDING; fi
  echo "$s $st"; done | tee -a "$MASTER"
