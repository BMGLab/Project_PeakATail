#!/bin/bash
# Stage 3 (Laughney) -- ONE sample end to end on the FROZEN snapshot (tools/pa-polya-run-4efeb125):
#   (a) ema run  clip_seeded + --ip-filter filter (pre-registered default, manuscript/13 s1)   [SKIPPED in cohort mode]
#   (b) inject CONFIRMED GEX cell-type labels (confirmed_labels.tsv, LABEL_POLICY.md; exact obs_name == full cell id),
#       drop unconfirmed/unlabelled cells, drop types < MIN_CELLS; restrict var to the tested universe
#   (c) ema switch diff  RELIABLE config (manuscript/14): fisher, --count-mode cells, --marker-top-n 0,
#       --counts-layer counts, --cluster-key celltype, --fdr 0.05, all pairwise cell-type pairs
#   (d) N_PERM label-shuffle nulls through the IDENTICAL switch diff  -> <sample>/null/perm_XX/
#   (e) ema switch length (PDUI classic) if it finishes within LENGTH_TIMEOUT s, else skipped + recorded  [SKIPPED in cohort mode]
# Usage:  laughney_per_sample.sh [--from-cohort-run] <dataset_id> [THREADS=8] [N_PERM=5]
# Env:    STAGE3_OUT_ROOT (default /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney), MIN_CELLS=20, LENGTH_TIMEOUT=600,
#         TIER_FILTER=none|tier1|tier1_ge2 (restrict tested PAS to clip-supported tier-1 / >=2 molecules; default none = both tiers)
#         STAGE3_STOP_AFTER=run|label|all (default all). 'run' = ema run only (label-independent, the expensive step);
#             'label' = run + label injection/report (no switch tests).
# COHORT MODE (design A, manuscript/13 addendum item 1; --from-cohort-run or STAGE3_COHORT_RUN=<cohort_run dir>):
#         no per-sample ema run; H5 = <cohort_run>/run/07_clustering/<id>/clusters.h5ad (unified PAS ids),
#         PASBED = <cohort_run>/run/unified/multi_sample_merged.bed (BED6 superset; 100% coverage asserted),
#         universe = <cohort_run>/universe/universe.tsv (stage3_build_universe.py, UNIVERSE_POLICY.md; TIER_FILTER forced to tier1_ge2),
#         OUT_ROOT default /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/switch, no switch length step.
#         STAGE3_SWITCH_EXTRA = extra `switch diff` args (SMOKE TESTS ONLY, recorded in the manifest; must be empty for real runs).
# Layout: $OUT_ROOT/<sample>/{run, labelled, true, null/perm_XX, length, logs, run_manifest.json, summary.json, DONE.ok|FAILED.err|NOPAIRS.ok}
# Resumable: RUN.ok / LABEL.ok / DONE.ok markers per step; a sample with DONE.ok is skipped entirely.
# Every number this produces is PROVISIONAL until the verifier passes.
set -uo pipefail
export LC_ALL=C PEAKATAIL_NO_TIMESTAMP=1

COHORT_MODE=0
if [[ ${1:-} == --from-cohort-run ]]; then COHORT_MODE=1; shift; fi
S=${1:?usage: laughney_per_sample.sh [--from-cohort-run] <dataset_id> [THREADS] [N_PERM]}
THREADS=${2:-8}
N_PERM=${3:-5}
MIN_CELLS=${MIN_CELLS:-20}
TIER_FILTER=${TIER_FILTER:-none}     # none | tier1 (pasbed col5>0) | tier1_ge2 (pre-registered precision default); decide BEFORE results
LENGTH_TIMEOUT=${LENGTH_TIMEOUT:-600}
STOP_AFTER=${STAGE3_STOP_AFTER:-all}
case "$STOP_AFTER" in run|label|all) ;; *) echo "STAGE3_STOP_AFTER must be run|label|all"; exit 2;; esac
SEQLEN=98                       # 10x 3' v2 R2 = 98 bp (verified on 20,000 reads, discovery)
COHORT_RUN=${STAGE3_COHORT_RUN:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/cohort_run}
[[ -n ${STAGE3_COHORT_RUN:-} ]] && COHORT_MODE=1
SWITCH_EXTRA=${STAGE3_SWITCH_EXTRA:-}

WD=/mnt/ssd1/Projects/PeakATail_wd
COMMIT=4efeb1252e6d7c7d51b252b443b5be5947ae000f
SNAP=$WD/tools/pa-polya-run-${COMMIT:0:8}
export PYTHONPATH=$SNAP
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
VENV_PY=$WD/tools/PeakATail/.venv/bin/python
GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
FA=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# CONFIRMED labels (LABEL_POLICY.md, stage3_confirm_labels.py); the unconfirmed curated parquet is no longer injected
# (verifier fix 2026-08-21: the dry run used it; addendum item 5 forbids it for any reported switch statistic)
LABELS=${CONFIRMED_LABELS:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/labels/confirmed_labels.tsv}
TABLE=$WD/scripts/stage3/laughney_samples.tsv
INJECT=$WD/scripts/stage3/stage3_label_inject.py
SUMMARISE=$WD/scripts/stage3/stage3_summarise.py
UNIVERSE_POLICY=$WD/scripts/stage3/UNIVERSE_POLICY.md
if [[ $COHORT_MODE -eq 1 ]]; then
  OUT_ROOT=${STAGE3_OUT_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/switch}
  LINK_ROOT=""                              # WD/results/stage3_laughney_v2 -> /mnt/ssd0/.../stage3_laughney_v2 already covers switch/
  TIER_FILTER=tier1_ge2
  UNIVERSE=$COHORT_RUN/universe/universe.tsv
else
  OUT_ROOT=${STAGE3_OUT_ROOT:-/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney}
  LINK_ROOT=$WD/results/stage3_laughney
  UNIVERSE=""
fi

# ---- hard guards -----------------------------------------------------------------
case "$OUT_ROOT" in /mnt/ssd2*) echo "REFUSING: output root on /mnt/ssd2 (99% full, read-only inputs)"; exit 2;; esac
[[ -d $SNAP ]] || { echo "frozen snapshot missing: $SNAP"; exit 2; }
[[ $(git -C "$SNAP" rev-parse HEAD) == "$COMMIT" ]] || { echo "snapshot HEAD != $COMMIT"; exit 2; }
BAM=$(awk -F'\t' -v s="$S" '$1==s{print $6}' "$TABLE")
[[ -n $BAM && -s $BAM && -s $BAM.bai ]] || { echo "sample $S not in $TABLE or BAM/index missing: $BAM"; exit 2; }
PATIENT=$(awk -F'\t' -v s="$S" '$1==s{print $2}' "$TABLE"); GROUP=$(awk -F'\t' -v s="$S" '$1==s{print $3}' "$TABLE")
for f in "$GTF" "$FA" "$FA.fai" "$LABELS" "$INJECT" "$SUMMARISE"; do [[ -s $f ]] || { echo "missing input: $f"; exit 2; }; done
if [[ $COHORT_MODE -eq 1 ]]; then
  [[ -s $COHORT_RUN/RUN.ok || ${STAGE3_ALLOW_NO_RUNOK:-0} == 1 ]] || { echo "cohort run not finished: $COHORT_RUN/RUN.ok missing"; exit 2; }
  [[ -s $COHORT_RUN/run/07_clustering/$S/clusters.h5ad ]] || { echo "missing cohort h5ad for $S"; exit 2; }
  [[ -s $COHORT_RUN/run/unified/multi_sample_merged.bed ]] || { echo "missing unified bed"; exit 2; }
  [[ -s $UNIVERSE ]] || { echo "missing universe: $UNIVERSE (run stage3_build_universe.py first)"; exit 2; }
  [[ -s $UNIVERSE_POLICY ]] || { echo "missing $UNIVERSE_POLICY"; exit 2; }
fi

OUT=$OUT_ROOT/$S
mkdir -p "$OUT/logs"
if [[ -n $LINK_ROOT ]]; then mkdir -p "$LINK_ROOT"; [[ -e $LINK_ROOT/$S ]] || ln -s "$OUT" "$LINK_ROOT/$S"; fi
if [[ -f $OUT/DONE.ok ]]; then echo "$S: DONE.ok present, skipping"; exit 0; fi
rm -f "$OUT/FAILED.err"
fail() { echo "$S: FAILED at $1 $(date)" | tee -a "$OUT/logs/driver.log"; echo "step=$1 $(date)" > "$OUT/FAILED.err"; exit 1; }
log() { echo "[$(date '+%F %T')] $S: $*" | tee -a "$OUT/logs/driver.log"; }

# the import-path assertion from stage2_final_launch.sh, run before EVERY ema job
# (verifier 2026-08-21: output tee'd into driver.log so the assertion is recorded next to the run, not only on stdout; pipefail is on)
assert_code() { "$VENV_PY" -c "import ema.strategies.clip_seeded as m, ema.switch_test.runner as r; assert '$SNAP' in m.__file__ and '$SNAP' in r.__file__, (m.__file__, r.__file__); print('[' + __import__('time').strftime('%F %T') + '] code:', m.__file__, r.__file__)" 2>&1 | tee -a "$OUT/logs/driver.log" || fail "import-path-assert"; }
EMA=("$VENV_PY" -c 'from ema.cli import main; main()')

log "start | commit $COMMIT | snapshot $SNAP | cohort_mode $COHORT_MODE | threads $THREADS | n_perm $N_PERM | min_cells $MIN_CELLS | tier_filter $TIER_FILTER | stop_after $STOP_AFTER | bam $BAM"
log "load: $(uptime | sed 's/.*load/load/') | free_gb: $(free -g | awk '/Mem:/{print $7}') | switch_extra '${SWITCH_EXTRA}'"

# ---- run manifest (ours; ema's own run_manifest.json lives inside run/) -----------
"$VENV_PY" - "$OUT/run_manifest.json" <<PY
import json, os, sys, socket, time, hashlib
def md5(p):
    h=hashlib.md5(); h.update(open(p,"rb").read()); return h.hexdigest()
bam="$BAM"; st=os.stat(bam)
cohort=int("$COHORT_MODE")==1
m=dict(sample="$S", patient="$PATIENT", group="$GROUP", bam=bam, bam_bytes=st.st_size, bam_mtime=time.strftime("%F %T", time.localtime(st.st_mtime)),
  mode="cohort (design A: one cohort-mode ema run, unified PAS ids)" if cohort else "per-sample (design B, superseded for replication)",
  tool="PeakATail (ema) frozen snapshot", tool_commit="$COMMIT", snapshot_path="$SNAP", python="$VENV_PY",
  gtf="$GTF", genome_fasta="$FA", labels_file="$LABELS", labels_kind="confirmed_labels.tsv (LABEL_POLICY.md)", labels_mtime=time.strftime("%F %T", time.localtime(os.stat("$LABELS").st_mtime)),
  labels_md5=md5("$LABELS"),
  label_join="exact obs_name == confirmed_labels.tsv cell ('<dataset_id>_<barcode>', confirmed==True only); unconfirmed/unlabelled cells dropped; cell types < $MIN_CELLS cells dropped",
  tier_filter="$TIER_FILTER",
  run_flags="--peak-strategy clip_seeded --ip-filter --ip-filter-mode filter --barcode-tag CB --cb-len 16 --seq-len $SEQLEN --ignore-chro MT --plot-engine none --no-progress --threads $THREADS (tool defaults otherwise: --polya-min-umis 1, both tiers emitted)",
  switch_diff_flags="--strategy fisher --count-mode cells --marker-top-n 0 --counts-layer counts --cluster-key celltype --fdr 0.05 --no-plots (all pairwise cell-type pairs; --min-cells-per-group tool default 10; --isoform-agg tool default per_gene)" + (" EXTRA(smoke only): $SWITCH_EXTRA" if "$SWITCH_EXTRA".strip() else ""),
  switch_extra_args="$SWITCH_EXTRA",
  null="label-shuffle: $N_PERM perms of obs['celltype'], numpy default_rng(seed=k) k=1..$N_PERM, identical switch diff",
  length=("not run in cohort mode" if cohort else "ema switch length --strategy classic, skipped if > $LENGTH_TIMEOUT s"),
  threads=$THREADS, env=dict(OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS"), PYTHONPATH=os.environ.get("PYTHONPATH")),
  host=socket.gethostname(), started=time.strftime("%F %T"), status="PROVISIONAL -- not verified")
if cohort:
    cr="$COHORT_RUN"
    m["cohort_run"]=cr
    m["cohort_run_ok"]=open(cr+"/RUN.ok").read().strip() if os.path.exists(cr+"/RUN.ok") else "MISSING (STAGE3_ALLOW_NO_RUNOK smoke)"
    m["h5ad"]=cr+"/run/07_clustering/$S/clusters.h5ad"; m["h5ad_mtime"]=time.strftime("%F %T", time.localtime(os.stat(m["h5ad"]).st_mtime))
    m["pasbed"]=cr+"/run/unified/multi_sample_merged.bed"; m["pasbed_md5"]=md5(m["pasbed"])
    m["universe_tsv"]="$UNIVERSE"; m["universe_md5"]=md5("$UNIVERSE")
    m["universe_policy"]="$UNIVERSE_POLICY"; m["universe_policy_md5"]=md5("$UNIVERSE_POLICY")
    us=os.path.join(os.path.dirname("$UNIVERSE"),"universe_summary.json")
    if os.path.exists(us):
        u=json.load(open(us)); m["universe_counts"]=u.get("counts"); m["universe_policy_md5_at_build"]=u.get("policy",{}).get("md5")
json.dump(m, open(sys.argv[1],"w"), indent=1)
PY

# ---- (a) ema run ---------------------------------------------------------------------
if [[ $COHORT_MODE -eq 1 ]]; then
  H5=$COHORT_RUN/run/07_clustering/$S/clusters.h5ad
  PASBED=$COHORT_RUN/run/unified/multi_sample_merged.bed
  log "(a) cohort mode: no per-sample ema run; h5ad $H5 ; pasbed $PASBED ; universe $UNIVERSE"
else
RUN=$OUT/run
if [[ ! -f $OUT/RUN.ok ]]; then
  rm -rf "$RUN"; mkdir -p "$RUN"
  cat > "$OUT/run_config.yaml" <<YML
# Stage 3 Laughney per-sample run -- written by laughney_per_sample.sh; CLI flags override keys
gtf: $GTF
seqlen: $SEQLEN
cb_len: 16
barcode_tag: CB
datasets:
- id: $S
  merge_strategy: none
  bams:
  - $BAM
YML
  assert_code
  log "(a) ema run -> $RUN"
  /usr/bin/time -v -o "$OUT/runtime_mem_run.txt" \
    "${EMA[@]}" run -c "$OUT/run_config.yaml" --output "$RUN" --threads "$THREADS" \
      --gtf "$GTF" --barcode-tag CB --cb-len 16 --seq-len "$SEQLEN" --ignore-chro MT \
      --plot-engine none --no-progress --peak-strategy clip_seeded \
      --ip-filter --genome-fasta "$FA" --ip-filter-mode filter > "$OUT/logs/ema_run.log" 2>&1
  rc=$?
  H5=$RUN/07_clustering/$S/clusters.h5ad
  [[ $rc -eq 0 && -s $RUN/pasbed.bed && -s $H5 ]] || { log "(a) ema exit $rc; pasbed/h5ad present? $(ls $RUN/pasbed.bed $H5 2>&1 | tr '\n' ' ')"; fail "ema-run"; }
  echo "exit=0 commit=$COMMIT $(date) $(grep -E 'Elapsed|Maximum resident' "$OUT/runtime_mem_run.txt" | tr -s ' ' | tr '\n' ';')" > "$OUT/RUN.ok"
  log "(a) done: $(cat "$OUT/RUN.ok")"
else
  log "(a) RUN.ok present, skipping ema run"
fi
H5=$RUN/07_clustering/$S/clusters.h5ad
PASBED=$RUN/pasbed.bed
if [[ $STOP_AFTER == run ]]; then echo "stopped after ema run $(date) (STAGE3_STOP_AFTER=run)" > "$OUT/STOPPED_AFTER.run"; log "STOP_AFTER=run: ema run complete, switch tests deferred"; exit 0; fi
fi

# ---- (b) label injection + universe restriction + null inputs ---------------------------
LAB=$OUT/labelled
if [[ ! -f $OUT/LABEL.ok ]]; then
  rm -rf "$LAB" "$OUT/NOPAIRS.ok"; mkdir -p "$LAB"
  log "(b) label injection"
  UARG=(); [[ -n $UNIVERSE ]] && UARG=(--universe-tsv "$UNIVERSE")
  "$VENV_PY" "$INJECT" --h5ad "$H5" --confirmed-labels "$LABELS" --sample "$S" --sample-table "$TABLE" \
      --pasbed "$PASBED" --out-dir "$LAB" --min-cells "$MIN_CELLS" --n-perm "$N_PERM" --tier-filter "$TIER_FILTER" "${UARG[@]}" > "$OUT/logs/label_inject.log" 2>&1
  rc=$?
  cat "$OUT/logs/label_inject.log" | tee -a "$OUT/logs/driver.log"
  if [[ $rc -eq 3 ]]; then
    echo "no pair: fewer than 2 cell types with >= $MIN_CELLS confirmed cells $(date)" > "$OUT/NOPAIRS.ok"
    "$VENV_PY" "$SUMMARISE" --sample-dir "$OUT" --sample "$S" --fdr 0.05 --cohort-log "$COHORT_RUN/logs/ema_run.log" 2>&1 | tee -a "$OUT/logs/driver.log"
    echo "exit=0 commit=$COMMIT $(date) NOPAIRS PROVISIONAL-until-verified" > "$OUT/DONE.ok"; log "DONE (NOPAIRS)"; exit 0
  fi
  [[ $rc -eq 0 ]] || fail "label-inject"
  echo "ok $(date)" > "$OUT/LABEL.ok"
else
  log "(b) LABEL.ok present, skipping"
fi
TRUE_H5=$LAB/$S.labelled.h5ad
if [[ $STOP_AFTER == label ]]; then echo "stopped after label injection $(date) (STAGE3_STOP_AFTER=label)" > "$OUT/STOPPED_AFTER.label"; log "STOP_AFTER=label: labels injected, switch tests deferred"; exit 0; fi

# ---- (c)+(d) switch diff, reliable config; identical call for TRUE and every null perm ---
run_diff() {  # h5ad outdir tag
  local h5=$1 od=$2 tag=$3
  if [[ -f $od/DONE.ok ]]; then log "$tag: DONE.ok, skip"; return 0; fi
  rm -rf "$od"; mkdir -p "$od"
  assert_code
  local t0=$(date +%s)
  /usr/bin/time -v -o "$OUT/runtime_mem_switch_${tag//\//_}.txt" \
  "${EMA[@]}" switch diff -i "$h5" --pasbed "$PASBED" --gtf "$GTF" \
      --cluster-key celltype --counts-layer counts \
      --strategy fisher --count-mode cells --marker-top-n 0 --fdr 0.05 \
      --no-plots --threads "$THREADS" $SWITCH_EXTRA -o "$od" > "$od/run.log" 2>&1
  local rc=$?
  local n=$(ls "$od"/differential/fisher_*_vs_*.tsv 2>/dev/null | wc -l)
  if [[ $rc -eq 0 && $n -gt 0 ]]; then
    echo "{\"sample\":\"$S\",\"run\":\"$tag\",\"rc\":0,\"elapsed_s\":$(( $(date +%s) - t0 )),\"n_pair_files\":$n,\"h5ad\":\"$h5\",\"pasbed\":\"$PASBED\",\"commit\":\"$COMMIT\"}" > "$od/DONE.ok"
    log "$tag: OK ($n pair files, $(( $(date +%s) - t0 )) s)"; return 0
  fi
  echo "rc=$rc n_pair_files=$n $(date)" > "$od/FAILED.err"; log "$tag: FAILED rc=$rc (see $od/run.log)"; return 1
}
run_diff "$TRUE_H5" "$OUT/true" "true" || fail "switch-diff-true"
nfail=0
for k in $(seq 1 "$N_PERM"); do
  kk=$(printf '%02d' "$k")
  run_diff "$LAB/perms/perm_$kk.h5ad" "$OUT/null/perm_$kk" "null/perm_$kk" || nfail=$((nfail+1))
done
[[ $nfail -eq 0 ]] || fail "switch-diff-null($nfail failed)"

# ---- (e) switch length (PDUI), bounded; not in cohort mode --------------------------------
if [[ $COHORT_MODE -eq 1 ]]; then
  [[ -f $OUT/LENGTH.skipped ]] || echo "skipped: not run in cohort mode (switch tables + nulls only) $(date)" > "$OUT/LENGTH.skipped"
elif [[ ! -f $OUT/LENGTH.ok && ! -f $OUT/LENGTH.skipped ]]; then
  rm -rf "$OUT/length"; mkdir -p "$OUT/length"; assert_code
  log "(e) switch length (timeout ${LENGTH_TIMEOUT}s)"
  /usr/bin/time -v -o "$OUT/runtime_mem_length.txt" \
  timeout "$LENGTH_TIMEOUT" "${EMA[@]}" switch length -i "$TRUE_H5" --pasbed "$PASBED" \
      --cluster-key celltype --counts-layer counts --strategy classic \
      --no-plots --threads "$THREADS" -o "$OUT/length" > "$OUT/length/run.log" 2>&1
  rc=$?
  if [[ $rc -eq 0 && -n $(ls "$OUT/length"/pdui_*.tsv 2>/dev/null) ]]; then
    echo "ok $(date) $(grep Elapsed "$OUT/runtime_mem_length.txt")" > "$OUT/LENGTH.ok"; log "(e) length OK"
  elif [[ $rc -eq 124 ]]; then
    echo "skipped: exceeded ${LENGTH_TIMEOUT}s $(date)" > "$OUT/LENGTH.skipped"; log "(e) length SKIPPED (timeout)"
  else
    echo "failed rc=$rc $(date)" > "$OUT/LENGTH.failed"; log "(e) length FAILED rc=$rc (non-fatal, see length/run.log)"
  fi
fi

# ---- summary -------------------------------------------------------------------------
"$VENV_PY" "$SUMMARISE" --sample-dir "$OUT" --sample "$S" --fdr 0.05 --cohort-log "$COHORT_RUN/logs/ema_run.log" 2>&1 | tee -a "$OUT/logs/driver.log" || fail "summarise"
rm -f "$OUT"/STOPPED_AFTER.*
echo "exit=0 commit=$COMMIT $(date) PROVISIONAL-until-verified" > "$OUT/DONE.ok"
log "DONE"
