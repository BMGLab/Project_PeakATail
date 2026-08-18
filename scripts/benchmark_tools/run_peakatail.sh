#!/usr/bin/env bash
# run_peakatail.sh -- REFERENCE WRAPPER for the benchmark_tools harness.
#
# Runs PeakATail (`ema run`) on the pbmc_10k_v3 CellRanger BAM and emits the
# standard result directory defined in README.md section 2. Every competitor
# wrapper should copy this file's shape.
#
# STATUS 2026-08-18: written and reviewed, NOT YET EXECUTED -- the 41 GiB BAM
#   was still downloading (~14%) when this was written. Preflight checks below
#   refuse to start on an incomplete BAM.
#
# TOOL   /mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
#        branch biolab-manuscript, commit 18678ef, version 0.2.0 (`ema --version`)
#
# VENV CAVEAT (verified 2026-08-18, details in status/peakatail.md):
#   tools/PeakATail/.venv is Python 3.10.12 with a STALE console script: the
#   old egg-info entry point `ema = ema.main:main` crashes on any invocation
#   (tries to open BAM "None"). `pip install -e . --no-deps` cannot refresh it
#   because pyproject now pins requires-python >=3.11. The current click CLI
#   (`ema.cli:main`) imports and runs fine on 3.10, so this wrapper bypasses
#   the console script and calls the click group directly through the venv
#   python (the `ema_cli` function below). Missing deps of the biolab-manuscript
#   branch (sortedcontainers, click, rich, questionary, tqdm, PyYAML, psutil,
#   statsmodels, pyarrow) were pip-installed into the venv on 2026-08-18.
#
# DATASET FACTS the flags below depend on (verified from the BAM itself):
#   - chrom naming: bare Ensembl (1..22, X, Y, MT + GL/KI scaffolds, 194 @SQ).
#     Matches the GTF and the PolyASite references directly -- so the chr-name
#     translation below is an idempotent no-op safety net, not a conversion.
#   - R2 read length 91 bp (CellRanger 3.0.0, 10x v3 chemistry) -> --seq-len 91
#   - cell barcode: CB tag, 16 bp -> --barcode-tag CB --cb-len 16
#
# FAIRNESS: no --atlas snapping. The benchmark scores de-novo calls; snapping
#   PeakATail's coordinates onto the reference atlas that evaluate_vs_reference
#   later scores against would be circular.

set -euo pipefail
export LC_ALL=C                    # tr_TR locale corrupts sort order
export PEAKATAIL_NO_TIMESTAMP=1    # deterministic output dir (ema/cli/common.py)

# ---------------------------------------------------------------------------
# knobs (env-overridable; THREADS must be identical across tools per dataset)
# ---------------------------------------------------------------------------
WD=/mnt/ssd1/Projects/PeakATail_wd
DATASET=${DATASET:-pbmc_10k_v3}
THREADS=${THREADS:-16}

BAM=$WD/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
BAM_BYTES_EXPECTED=44198457238     # server Content-Length; 10x publishes no md5
GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf

TOOLDIR=$WD/tools/PeakATail
VENV_PY=$TOOLDIR/.venv/bin/python
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

OUT=$WD/results/benchmark_tools/$DATASET/peakatail
RUN=$OUT/run                       # ema's native output tree, kept for audit

# click CLI via the venv python; see VENV CAVEAT above.
ema_cli() { "$VENV_PY" -c 'from ema.cli import main; main()' "$@"; }

# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------
[[ -x $VENV_PY ]]  || { echo "FATAL: $VENV_PY missing" >&2; exit 1; }
[[ -r $GTF ]]      || { echo "FATAL: GTF $GTF unreadable" >&2; exit 1; }
[[ -r $BAM ]]      || { echo "FATAL: BAM $BAM missing" >&2; exit 1; }
bytes=$(stat -c%s "$BAM")
[[ $bytes -eq $BAM_BYTES_EXPECTED ]] || {
  echo "FATAL: BAM is $bytes B, expected $BAM_BYTES_EXPECTED (download incomplete?)" >&2
  exit 1; }
samtools quickcheck -v "$BAM" || { echo "FATAL: samtools quickcheck failed" >&2; exit 1; }
command -v /usr/bin/time >/dev/null || { echo "FATAL: /usr/bin/time missing" >&2; exit 1; }

mkdir -p "$OUT" "$RUN"

# ---------------------------------------------------------------------------
# tool_version.txt
# ---------------------------------------------------------------------------
{ ema_cli --version
  echo "git: $(git -C "$TOOLDIR" rev-parse HEAD) ($(git -C "$TOOLDIR" rev-parse --abbrev-ref HEAD))"
  echo "python: $("$VENV_PY" --version 2>&1)"
} > "$OUT/tool_version.txt"

# ---------------------------------------------------------------------------
# the run  (/usr/bin/time -v -> runtime_mem.txt; stdout+stderr -> log)
# ---------------------------------------------------------------------------
/usr/bin/time -v -o "$OUT/runtime_mem.txt" \
  "$VENV_PY" -c 'from ema.cli import main; main()' run \
    --bam-dir "$BAM" \
    --gtf "$GTF" \
    --output "$RUN" \
    --threads "$THREADS" \
    --barcode-tag CB \
    --cb-len 16 \
    --seq-len 91 \
    --ignore-chro 'MT' \
    --plot-engine none \
    --no-progress \
  > "$OUT/log" 2>&1

# ---------------------------------------------------------------------------
# standard outputs
#   ema writes (verified on the Laughney runs): pasbed.bed = BED6 peak
#   INTERVALS (chrom start end id score strand, chrom naming follows the BAM),
#   annotated_matrix.mtx = cell x PAS counts, plus per-stage dirs and logs.
# ---------------------------------------------------------------------------
PASBED=$RUN/pasbed.bed
[[ -s $PASBED ]] || { echo "FATAL: $PASBED missing/empty -- see $OUT/log" >&2; exit 1; }

# interval -> 3'-most base by strand (README sec 2; same awk as null_control.py),
# then the harness chr-name helper (idempotent here: BAM is already Ensembl;
# --unmapped drop guards against oddball contigs), then LC_ALL=C sort.
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}' \
    "$PASBED" \
  | "$HERE/normalize_chroms.py" --to ensembl --unmapped drop \
  | sort -k1,1 -k2,2n > "$OUT/pas.bed"

if [[ -s $RUN/annotated_matrix.mtx ]]; then
  cp "$RUN/annotated_matrix.mtx" "$OUT/counts.mtx"
fi

echo "DONE: $(wc -l < "$OUT/pas.bed") PAS -> $OUT/pas.bed" >&2
