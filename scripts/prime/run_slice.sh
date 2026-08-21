#!/usr/bin/env bash
# scripts/prime/run_slice.sh --label L [--code DIR] [--threads N] [--bam BAM]
#                            [--gtf GTF] [--seq-len N] [--] <extra ema flags...>
#
# Runs the PeakATail caller on the PBMC chr19+21 dev slice and writes everything
# under $OUTROOT/<label>/.  Fixed flags are copied verbatim from
# scripts/benchmark_tools/stage2_final_launch.sh so a slice run differs from the
# manuscript arms only in the input BAM (and whatever extra flags are passed).
#
#   --code    tree to import `ema` from.  Default $PRIME_CODE (branch
#             peakAtail-prime).  Use $V2_CODE for a v2 reference run.
#   --threads worker ceiling.  Default 8 (hard cap: the Stage-3 chain owns the box).
#
# Outputs:
#   $OUTROOT/<label>/run/          the caller's output tree (PEAKATAIL_NO_TIMESTAMP=1,
#                                  so the leaf is exactly `run`, not `run_<ts>`)
#   $OUTROOT/<label>/pas.bed       scorer-ready 1-bp points, produced by the SAME
#                                  awk + normalize_chroms.py the launcher uses
#   $OUTROOT/<label>/pas_tier1.bed pas_tier2.bed pas_tier1_ge2mol.bed
#   $OUTROOT/<label>/{cmd.txt,code_version.txt,modpath.txt,runtime_mem.txt,log,DONE.ok}
set -uo pipefail
source "$(dirname "$0")/common.sh"

LABEL=""; CODE=$PRIME_CODE; THREADS=8; BAM=$SLICE_BAM; GTF=$HUM_GTF; SEQLEN=91
while [ $# -gt 0 ]; do
  case "$1" in
    --label)   LABEL=$2; shift 2;;
    --code)    CODE=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --bam)     BAM=$2; shift 2;;
    --gtf)     GTF=$2; shift 2;;
    --seq-len) SEQLEN=$2; shift 2;;
    --)        shift; break;;
    *)         break;;
  esac
done
[ -n "$LABEL" ] || { echo "usage: run_slice.sh --label L [opts] [-- extra ema flags]"; exit 2; }
[ "$THREADS" -le 8 ] || { echo "REFUSING --threads $THREADS: harness cap is 8"; exit 2; }
EXTRA=("$@")

OUT=$OUTROOT/$LABEL
mkdir -p "$OUT"
rm -rf "$OUT/run" "$OUT/DONE.ok" "$OUT/FAILED.err"

guard_load || exit 3
assert_code "$CODE" | tee "$OUT/modpath.txt" || { echo "code assertion FAILED"; exit 4; }
{ echo "code_dir=$CODE"; git -C "$CODE" rev-parse HEAD 2>/dev/null; \
  git -C "$CODE" rev-parse --abbrev-ref HEAD 2>/dev/null; \
  git -C "$CODE" status --short 2>/dev/null; } > "$OUT/code_version.txt"

export PYTHONPATH=$CODE
export PEAKATAIL_NO_TIMESTAMP=1          # stable output leaf name; required for identity checks
export PYTHONDONTWRITEBYTECODE=1         # never write .pyc into a read-only frozen snapshot
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2

ARGS=(run --bam-dir "$BAM" --gtf "$GTF" --output "$OUT/run" --threads "$THREADS"
      --barcode-tag CB --cb-len 16 --seq-len "$SEQLEN" --ignore-chro 'MT'
      --plot-engine none --no-progress --peak-strategy clip_seeded "${EXTRA[@]}")
printf '%q ' "$PY" -c 'from ema.cli import main; main()' "${ARGS[@]}" > "$OUT/cmd.txt"; echo >> "$OUT/cmd.txt"

echo "=== $LABEL start $(date -Is) | code $CODE | extra: ${EXTRA[*]:-<none>} ==="
cd "$OUT" || exit 1                       # the CLI drops switch_combined/ logs in CWD
/usr/bin/time -v -o "$OUT/runtime_mem.txt" \
  "$PY" -c 'from ema.cli import main; main()' "${ARGS[@]}" > "$OUT/log" 2>&1
RC=$?
echo "=== $LABEL ema exit $RC $(date -Is) ==="
if [ $RC -ne 0 ] || [ ! -s "$OUT/run/pasbed.bed" ]; then
  echo "exit=$RC" > "$OUT/FAILED.err"; tail -20 "$OUT/log"; exit 1
fi

# scorer-ready points: 3'-most base by strand, Ensembl contig names, sorted.
# Byte-for-byte the launcher's derivation.
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}' "$OUT/run/pasbed.bed" \
  | "$PY" "$NORM" --to ensembl --unmapped drop | sort -k1,1 -k2,2n > "$OUT/pas.bed"
awk -F'\t' '$5>0'  "$OUT/pas.bed" > "$OUT/pas_tier1.bed"
awk -F'\t' '$5==0' "$OUT/pas.bed" > "$OUT/pas_tier2.bed"
awk -F'\t' '$5>=2' "$OUT/pas.bed" > "$OUT/pas_tier1_ge2mol.bed"

printf 'pas %d | tier1 %d | tier2 %d | tier1>=2mol %d\n' \
  "$(wc -l < "$OUT/pas.bed")" "$(wc -l < "$OUT/pas_tier1.bed")" \
  "$(wc -l < "$OUT/pas_tier2.bed")" "$(wc -l < "$OUT/pas_tier1_ge2mol.bed")" | tee "$OUT/counts.txt"
grep -E 'Elapsed|Maximum resident' "$OUT/runtime_mem.txt"
echo "exit=0 code=$CODE commit=$(git -C "$CODE" rev-parse --short HEAD 2>/dev/null) $(date -Is)" > "$OUT/DONE.ok"
echo "=== $LABEL DONE -> $OUT ==="
