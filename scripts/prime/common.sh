#!/usr/bin/env bash
# scripts/prime/common.sh -- shared paths and guards for the peakAtail-prime dev harness.
# Source it; do not execute it.  Every prime script starts with `export LC_ALL=C`
# because the box's tr_TR locale silently corrupts GNU sort / bedtools ordering.
export LC_ALL=C

WD=/mnt/ssd1/Projects/PeakATail_wd
PY=$WD/tools/PeakATail/.venv/bin/python              # the ONLY interpreter used by this project
PRIME_CODE=$WD/tools/pa-prime                        # branch peakAtail-prime
V2_CODE=$WD/tools/pa-polya-run-9dfdefb3              # frozen 9dfdefb snapshot == "v2" (READ-ONLY)
V1_CODE=$WD/tools/pa-polya-run-4efeb125              # frozen 4efeb125 snapshot == "v1" (READ-ONLY)

OUTROOT=/mnt/ssd0/emaout/peakatail_benchmark/prime   # all prime run outputs live here
ARTROOT=$WD/results/prime                            # small artefacts (scores, identity reports)
SLICE_BAM=/mnt/ssd0/emaout/peakatail_benchmark/perf/slice/pbmc_chr19_21.bam

HUM_GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
HUM_FA=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
PBMC_BAM=$WD/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam

NORM=$WD/scripts/benchmark_tools/normalize_chroms.py
SCORE=$WD/scripts/benchmark_tools/score_tool.py
SIZES=$WD/data/references/chrom.sizes.nochr.filt
KINNEX=$WD/results/benchmark_tools/kinnex_truth/x3p

# score_tool.py reference arguments -- copied verbatim from
# scripts/benchmark_tools/stage2_final_launch.sh so prime numbers are row-for-row
# comparable with the v2 manuscript numbers.
HUMAN_DET=(--detected-atlas $WD/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed)
R=$WD/data/references; G=$WD/results/benchmark_tools/gse104556
MOUSE_REFS=(--atlas $R/atlases/polyasite2.GRCm38.96.rep_sites.bed6
            --tes $R/atlases/tes.protein_coding.GRCm38.102.bed6
            --genome $R/mouse/chrom.sizes.filt
            --genebodies $G/shared_refs/genebodies.merged.bed
            --detected-atlas $G/shared_refs/pas2.in_detected_genes.bed)

# assert_code <code_dir> -- refuse to run unless `ema` really resolves inside <code_dir>.
# This project has five worktrees of the same package; importing the wrong one has
# silently invalidated results before.  Every run calls this first.
assert_code() {
  local code=$1
  PYTHONPATH=$code PYTHONDONTWRITEBYTECODE=1 "$PY" - "$code" <<'PYEOF'
import sys, ema, ema.strategies.clip_seeded as cs, ema.countmatrix.polya as pol
code = sys.argv[1].rstrip('/')
for m in (ema, cs, pol):
    assert m.__file__.startswith(code + '/'), f'WRONG TREE: {m.__name__} -> {m.__file__} (wanted {code})'
print('code OK:', ema.__file__)
PYEOF
}

# guard_load -- refuse to start heavy work when the box is already busy.
# The Stage-3 replication chain owns 24 cores under stage3_laughney_v3/; this
# harness must never compete with it.
guard_load() {
  local maxload=${1:-60}
  local load; load=$(awk '{printf "%d", $1}' /proc/loadavg)
  local freeg; freeg=$(free -g | awk '/^Mem:/{print $7}')
  echo "[guard] load1=$load available_GB=$freeg"
  if [ "$load" -gt "$maxload" ]; then echo "[guard] REFUSING: load1 $load > $maxload"; return 1; fi
  if [ "$freeg" -lt 20 ]; then echo "[guard] REFUSING: only ${freeg}G available"; return 1; fi
  return 0
}
