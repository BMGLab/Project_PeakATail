#!/bin/bash
# Stage 3 (Laughney, design A) -- switch tables + nulls for EVERY GSM of the cohort run (addendum items 3-4).
#   1. waits for <cohort_run>/RUN.ok (aborts on FAILED.err)            4. laughney_per_sample.sh --from-cohort-run, CONC x THREADS
#   2. runs verify_cohort.py (non-fatal, for the verifier)               5. re-queues every FAILED GSM once
#   3. builds the universe (stage3_build_universe.py, UNIVERSE_POLICY.md)  6. aggregates switch/SUMMARY_ALL.{tsv,json}
# Whole tree pinned to CPUSET (24 logical CPUs; <= 24 concurrent threads on the shared box). BLAS = 1 thread.
# Usage: nohup scripts/stage3/run_cohort_switch_all.sh > <switch>/logs/run_all.out 2>&1 &
# Every number this produces is PROVISIONAL until the verifier passes.
set -uo pipefail
export LC_ALL=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
WD=/mnt/ssd1/Projects/PeakATail_wd
TABLE=$WD/scripts/stage3/laughney_samples.tsv
DRIVER=$WD/scripts/stage3/laughney_per_sample.sh
BUILD=$WD/scripts/stage3/stage3_build_universe.py
AGG=$WD/scripts/stage3/stage3_cohort_table.py
VENV_PY=$WD/tools/PeakATail/.venv/bin/python
SNAP=$WD/tools/pa-polya-run-4efeb125
export STAGE3_COHORT_RUN=${STAGE3_COHORT_RUN:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/cohort_run}
export STAGE3_OUT_ROOT=${STAGE3_OUT_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/switch}
export MIN_CELLS=${MIN_CELLS:-20} TIER_FILTER=tier1_ge2 STAGE3_STOP_AFTER=all
unset STAGE3_SWITCH_EXTRA
CONC=${CONC:-3}; THREADS=${THREADS:-8}; N_PERM=${N_PERM:-10}; CPUSET=${CPUSET:-28-51}
MIN_FREE_GB=${MIN_FREE_GB:-120}; WAIT_POLL=${WAIT_POLL:-300}; MAX_WAIT_H=${MAX_WAIT_H:-20}
export TMPDIR=/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/tmp
case "$STAGE3_OUT_ROOT" in /mnt/ssd2*) echo "REFUSING: output root on /mnt/ssd2"; exit 2;; esac
mkdir -p "$STAGE3_OUT_ROOT/logs" "$TMPDIR"
MASTER=$STAGE3_OUT_ROOT/logs/run_all.log
log() { echo "[$(date '+%F %T')] $*" | tee -a "$MASTER"; }
CR=$STAGE3_COHORT_RUN
log "=== run_cohort_switch_all start | conc $CONC x threads $THREADS = $((CONC*THREADS)) | n_perm $N_PERM | cpuset $CPUSET | min_cells $MIN_CELLS | tier_filter $TIER_FILTER | cohort_run $CR | out $STAGE3_OUT_ROOT"

# ---- 1. wait for the cohort run -------------------------------------------------------------------------------
t0=$(date +%s)
while true; do
  if [[ -s $CR/RUN.ok ]]; then log "cohort RUN.ok: $(cat $CR/RUN.ok)"; break; fi
  if [[ -s $CR/FAILED.err ]]; then log "cohort run FAILED: $(cat $CR/FAILED.err) -- nothing to do"; echo "cohort-failed $(date)" > "$STAGE3_OUT_ROOT/ABORTED.err"; exit 3; fi
  if (( $(date +%s) - t0 > MAX_WAIT_H*3600 )); then log "waited > ${MAX_WAIT_H} h for RUN.ok; giving up"; echo "wait-timeout $(date)" > "$STAGE3_OUT_ROOT/ABORTED.err"; exit 4; fi
  if ! pgrep -f "run -c $CR/cohort_17.yaml" > /dev/null; then
    sleep 120   # give the launcher time to write RUN.ok / FAILED.err after the process exits
    if [[ ! -s $CR/RUN.ok && ! -s $CR/FAILED.err ]]; then log "cohort ema process gone without RUN.ok/FAILED.err"; echo "cohort-process-lost $(date)" > "$STAGE3_OUT_ROOT/ABORTED.err"; exit 5; fi
  fi
  echo "[$(date '+%F %T')] waiting for $CR/RUN.ok (peakcalling done: $(ls $CR/run/peakcalling/*.neg.bed 2>/dev/null | wc -l)/17; h5ad: $(ls $CR/run/07_clustering/*/clusters.h5ad 2>/dev/null | wc -l))" >> "$MASTER"
  sleep "$WAIT_POLL"
done
rm -f "$STAGE3_OUT_ROOT/ABORTED.err"

# ---- 2. verify_cohort (non-fatal; the cohort task's own completion step) ------------------------------------------
export PYTHONPATH=$SNAP
if [[ -s $CR/verify_cohort.py && ! -s $CR/verify_cohort.json ]]; then
  log "verify_cohort.py -> $CR/verify_cohort.json"
  "$VENV_PY" "$CR/verify_cohort.py" "$CR/run" "$CR/verify_cohort.json" "$CR/runtime_mem_cohort.txt" "$TABLE" > "$STAGE3_OUT_ROOT/logs/verify_cohort.out" 2>&1 || log "verify_cohort.py exited non-zero (non-fatal; see logs/verify_cohort.out)"
fi

# ---- 3. universe --------------------------------------------------------------------------------------------------
if [[ ! -s $CR/universe/universe.tsv ]]; then
  log "building universe -> $CR/universe"
  "$VENV_PY" "$BUILD" --run "$CR/run" --out-dir "$CR/universe" --sample-table "$TABLE" > "$STAGE3_OUT_ROOT/logs/build_universe.out" 2>&1 \
    || { log "universe build FAILED (see logs/build_universe.out)"; echo "universe-failed $(date)" > "$STAGE3_OUT_ROOT/ABORTED.err"; exit 6; }
  tail -40 "$STAGE3_OUT_ROOT/logs/build_universe.out" | tee -a "$MASTER"
else
  log "universe present: $CR/universe/universe.tsv"
fi

# ---- 4. resources, then run -------------------------------------------------------------------------------------
FREE=$(free -g | awk '/Mem:/{print $7}')
log "load: $(uptime | sed 's/.*average: //') | free_gb: $FREE | nproc $(nproc) | ssd0 $(df -h /mnt/ssd0 | awk 'NR==2{print $4" free"}')"
while (( FREE < MIN_FREE_GB )); do log "free_gb $FREE < $MIN_FREE_GB; waiting 10 min"; sleep 600; FREE=$(free -g | awk '/Mem:/{print $7}'); done
PENDING=$(awk -F'\t' 'NR>1{print $1}' "$TABLE" | while read -r s; do
  if [[ -f $STAGE3_OUT_ROOT/$s/DONE.ok ]]; then echo "skip $s (DONE.ok)" >&2; else echo "$s"; fi; done)
log "pending: $(echo "$PENDING" | wc -w) samples"
export DRIVER THREADS N_PERM STAGE3_OUT_ROOT CPUSET
run_batch() {
  printf '%s\n' "$@" | xargs -P "$CONC" -I{} bash -c \
    'echo "start {} $(date)"; taskset -c "$CPUSET" "$DRIVER" --from-cohort-run {} "$THREADS" "$N_PERM" >> "$STAGE3_OUT_ROOT/logs/{}.log" 2>&1; rc=$?; echo "end {} exit $rc $(date)"' \
    2>&1 | tee -a "$MASTER"
}
[[ -n $PENDING ]] && run_batch $PENDING
# ---- 5. re-queue failures once ----------------------------------------------------------------------------------
FAILED=$(for s in $(awk -F'\t' 'NR>1{print $1}' "$TABLE"); do [[ -f $STAGE3_OUT_ROOT/$s/FAILED.err && ! -f $STAGE3_OUT_ROOT/$s/DONE.ok ]] && echo "$s"; done)
if [[ -n $FAILED ]]; then
  log "re-queue (once): $FAILED"
  for s in $FAILED; do mv "$STAGE3_OUT_ROOT/$s/FAILED.err" "$STAGE3_OUT_ROOT/$s/FAILED.first_attempt" 2>/dev/null; done
  run_batch $FAILED
fi
# ---- 6. aggregate ---------------------------------------------------------------------------------------------------
"$VENV_PY" "$AGG" --out-root "$STAGE3_OUT_ROOT" --sample-table "$TABLE" 2>&1 | tee -a "$MASTER"
log "=== run_cohort_switch_all finished"
for s in $(awk -F'\t' 'NR>1{print $1}' "$TABLE"); do
  if [[ -f $STAGE3_OUT_ROOT/$s/NOPAIRS.ok ]]; then st=NOPAIRS; elif [[ -f $STAGE3_OUT_ROOT/$s/DONE.ok ]]; then st=DONE; elif [[ -f $STAGE3_OUT_ROOT/$s/FAILED.err ]]; then st="FAILED($(cat $STAGE3_OUT_ROOT/$s/FAILED.err))"; else st=PENDING; fi
  echo "$s $st"; done | tee -a "$MASTER"
echo "finished $(date)" > "$STAGE3_OUT_ROOT/ALL.done"
