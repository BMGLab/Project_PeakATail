#!/bin/bash
# Stage 3 (Laughney, design A) -- cross-PATIENT replication filter over the cohort switch tables
# (manuscript/13 addendum item 2: unit of replication = patient; K>=2 patients same direction; discordant = exclude;
#  |dprop| >= 0.1 floor reported alongside; nulls = every GSM's label-shuffle perms, index-aligned combos).
# Runs scripts/reliability/replication_filter.py at PAS and gene level for K=2 (pre-registered) and K=3 (sensitivity),
# then stage3_replication_report.py.  Waits for <switch>/ALL.done unless STAGE3_REPL_NOWAIT=1.
# Usage: nohup scripts/stage3/run_replication.sh > <repl>/logs/run_replication.out 2>&1 &
# Every number this produces is PROVISIONAL until the verifier passes.
set -uo pipefail
export LC_ALL=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
WD=/mnt/ssd1/Projects/PeakATail_wd
VENV_PY=$WD/tools/PeakATail/.venv/bin/python
FILTER=$WD/scripts/reliability/replication_filter.py
REPORT=$WD/scripts/stage3/stage3_replication_report.py
TABLE=${STAGE3_SAMPLE_TABLE:-$WD/scripts/stage3/laughney_samples.tsv}
SWITCH=${STAGE3_SWITCH_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/switch}
REPL=${STAGE3_REPL_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/replication}
GENE_NAMES=${STAGE3_GENE_NAMES:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/replication/gene_names.GRCh38.99.tsv}
NOTE=${STAGE3_REPL_NOTE:-"Stage-3 Laughney cohort (design A, tier1_ge2 universe, confirmed labels); unit of replication = patient (addendum item 2); PROVISIONAL until verified"}
WAIT_POLL=${WAIT_POLL:-300}; MAX_WAIT_H=${MAX_WAIT_H:-30}
export TMPDIR=/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/tmp
case "$REPL" in /mnt/ssd2*) echo "REFUSING: output root on /mnt/ssd2"; exit 2;; esac
mkdir -p "$REPL/logs" "$TMPDIR"
LOG=$REPL/logs/run_replication.log
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
log "=== run_replication start | switch $SWITCH | out $REPL | table $TABLE"

# ---- 1. wait for the switch stage ------------------------------------------------------------------------------
t0=$(date +%s)
if [[ ${STAGE3_REPL_NOWAIT:-0} != 1 ]]; then
  while [[ ! -f $SWITCH/ALL.done ]]; do
    if [[ -f $SWITCH/ABORTED.err ]]; then log "switch stage ABORTED: $(cat $SWITCH/ABORTED.err)"; echo "switch-aborted $(date)" > "$REPL/ABORTED.err"; exit 3; fi
    if (( $(date +%s) - t0 > MAX_WAIT_H*3600 )); then log "waited > ${MAX_WAIT_H} h for ALL.done; giving up"; echo "wait-timeout $(date)" > "$REPL/ABORTED.err"; exit 4; fi
    echo "[$(date '+%F %T')] waiting for $SWITCH/ALL.done (DONE: $(ls $SWITCH/*/DONE.ok 2>/dev/null | wc -l) NOPAIRS: $(ls $SWITCH/*/NOPAIRS.ok 2>/dev/null | wc -l) FAILED: $(ls $SWITCH/*/FAILED.err 2>/dev/null | wc -l))" >> "$LOG"
    sleep "$WAIT_POLL"
  done
  log "switch ALL.done: $(cat $SWITCH/ALL.done)"
fi
rm -f "$REPL/ABORTED.err"

# ---- 2. sample list + patient map --------------------------------------------------------------------------------
# a GSM enters iff DONE.ok, not NOPAIRS, and true/differential has >= 1 fisher pair file; its nulls are null/perm_*/
SAMPLES=(); GROUPS_TSV=$REPL/sample_groups.tsv
printf 'sample\tgroup\tdataset_group\tstatus\tn_true_pairs\tn_null_perms\n' > "$GROUPS_TSV"
while IFS=$'\t' read -r ds patient grp rest; do
  [[ $ds == dataset_id ]] && continue
  d=$SWITCH/$ds
  ntrue=$(ls "$d"/true/differential/fisher_*_vs_*.tsv 2>/dev/null | wc -l)
  nnull=$(ls -d "$d"/null/perm_*/DONE.ok 2>/dev/null | wc -l)
  if [[ -f $d/NOPAIRS.ok ]]; then st=NOPAIRS
  elif [[ -f $d/DONE.ok && $ntrue -ge 1 ]]; then st=INCLUDED; SAMPLES+=("$ds=$d/true")
  elif [[ -f $d/FAILED.err ]]; then st=FAILED
  else st=MISSING; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ds" "$patient" "$grp" "$st" "$ntrue" "$nnull" >> "$GROUPS_TSV"
done < "$TABLE"
log "samples included: ${#SAMPLES[@]}  (map: $GROUPS_TSV)"; cat "$GROUPS_TSV" | tee -a "$LOG"
if (( ${#SAMPLES[@]} < 2 )); then log "fewer than 2 GSMs with switch tables; nothing to replicate"; echo "too-few-samples $(date)" > "$REPL/ABORTED.err"; exit 5; fi
SARGS=(); for s in "${SAMPLES[@]}"; do SARGS+=(--sample "$s"); done

# ---- 3. four configurations -------------------------------------------------------------------------------------
rc_all=0
for level in pas gene; do for K in 2 3; do
  cfg=${level}_K${K}; out=$REPL/$cfg
  if [[ -f $out/DONE.ok ]]; then log "$cfg: DONE.ok present, skipping"; continue; fi
  rm -rf "$out"; mkdir -p "$out"
  log "$cfg: start"
  /usr/bin/time -v -o "$out/runtime_mem.txt" "$VENV_PY" "$FILTER" "${SARGS[@]}" --out "$out" \
      --input-kind diff --strategy fisher --effect-col delta_proportion --level "$level" \
      --min-samples "$K" --fdr 0.05 --effect-floor 0.1 --discordant-policy exclude \
      --sample-group-tsv "$GROUPS_TSV" \
      --null-dirs "$SWITCH/{sample}/null/perm_*" --null-universe real \
      --note "$NOTE; config $cfg" > "$out/run.log" 2>&1
  rc=$?
  if [[ $rc -eq 0 && -s $out/summary.json ]]; then
    echo "exit=0 $(date) PROVISIONAL-until-verified" > "$out/DONE.ok"; log "$cfg: DONE ($(grep -m1 'replicated (K' "$out/run.log" | sed 's/.*INFO //'))"
  else
    echo "exit=$rc $(date)" > "$out/FAILED.err"; log "$cfg: FAILED rc=$rc (see $out/run.log)"; rc_all=1
  fi
done; done

# ---- 4. report ---------------------------------------------------------------------------------------------------
"$VENV_PY" "$REPORT" --rep-root "$REPL" --gene-names "$GENE_NAMES" > "$REPL/logs/report.out" 2>&1 || { log "report FAILED (see logs/report.out)"; rc_all=1; }
log "=== run_replication finished rc=$rc_all"
if [[ $rc_all -eq 0 ]]; then echo "finished $(date) PROVISIONAL-until-verified" > "$REPL/ALL.done"; else echo "finished-with-failures $(date)" > "$REPL/FAILED.err"; fi
exit $rc_all
