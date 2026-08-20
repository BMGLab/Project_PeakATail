#!/bin/bash
# Stage 2: ONE full re-run of the clip-seeded caller (feat/polya-evidence, rebased) on both datasets.
# Four arms run concurrently; each post-processes to pas.bed (tier tag preserved in col 5),
# scores itself with score_tool.py on the per-dataset denominator, and writes DONE.ok/FAILED.err.
set -uo pipefail
export LC_ALL=C PEAKATAIL_NO_TIMESTAMP=1
WD=/mnt/ssd1/Projects/PeakATail_wd
export PYTHONPATH=$WD/tools/pa-polya            # the rebased Stage-1 worktree (f0370f7)
VENV_PY=$WD/tools/PeakATail/.venv/bin/python
NORM=$WD/scripts/benchmark_tools/normalize_chroms.py
SCORE=$WD/scripts/benchmark_tools/score_tool.py
HUM_GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
HUM_FA=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
MOU_GTF=$WD/data/references/mouse/Mus_musculus.GRCm38.102.gtf
PBMC=$WD/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
M1=$WD/data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam
M2=$WD/data/benchmark/gse104556/starsolo/Mouse2_scRNAseq/Aligned.sortedByCoord.out.bam
R=$WD/data/references; G=$WD/results/benchmark_tools/gse104556
MOUSE_REFS=(--atlas $R/atlases/polyasite2.GRCm38.96.rep_sites.bed6 --tes $R/atlases/tes.protein_coding.GRCm38.102.bed6
            --genome $R/mouse/chrom.sizes.filt --genebodies $G/shared_refs/genebodies.merged.bed
            --detected-atlas $G/shared_refs/pas2.in_detected_genes.bed)
HUMAN_DET=(--detected-atlas $WD/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed)

# run_arm OUTDIR LABEL BAM GTF SEQLEN THREADS SPECIES [extra ema flags...]
run_arm() {
  local out=$1 label=$2 bam=$3 gtf=$4 seqlen=$5 threads=$6 species=$7; shift 7
  mkdir -p "$out"; local run=$out/run; rm -rf "$run"
  {
    echo "=== $label start $(date) | strategy clip_seeded | extra: $* ==="
    echo "code commit: $(git -C $WD/tools/pa-polya rev-parse HEAD)  (tree must stay frozen for the run)"
    "$VENV_PY" -c 'import ema.strategies.clip_seeded as m; assert "pa-polya" in m.__file__, m.__file__; print("code:", m.__file__)'
    /usr/bin/time -v -o "$out/runtime_mem.txt" \
      "$VENV_PY" -c 'from ema.cli import main; main()' run \
        --bam-dir "$bam" --gtf "$gtf" --output "$run" --threads "$threads" \
        --barcode-tag CB --cb-len 16 --seq-len "$seqlen" --ignore-chro 'MT' \
        --plot-engine none --no-progress --peak-strategy clip_seeded "$@" > "$out/log" 2>&1
    local rc=$?
    echo "=== $label ema exit $rc $(date) ==="
    [[ $rc -eq 0 && -s $run/pasbed.bed ]] || { echo "exit=$rc" > "$out/FAILED.err"; return 1; }
    # point-form pas.bed: strand-aware 3' base; KEEP col5 = clip-read count (tier tag: 0 = coverage-only)
    awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}' "$run/pasbed.bed" \
      | "$VENV_PY" "$NORM" --to ensembl --unmapped drop | sort -k1,1 -k2,2n > "$out/pas.bed"
    awk -F'\t' '$5>0'  "$out/pas.bed" > "$out/pas_tier1.bed"
    awk -F'\t' '$5==0' "$out/pas.bed" > "$out/pas_tier2.bed"
    [[ -s $run/annotated_matrix.mtx ]] && cp "$run/annotated_matrix.mtx" "$out/counts.mtx"
    echo "pas.bed $(wc -l < $out/pas.bed) | tier1 $(wc -l < $out/pas_tier1.bed) | tier2 $(wc -l < $out/pas_tier2.bed)"
    # score all three (tier2 may be empty-ish; score_tool tolerates)
    for v in pas pas_tier1 pas_tier2; do
      local lab=${label}_${v#pas}; lab=${lab%_}
      if [[ $species == human ]]; then "$VENV_PY" "$SCORE" "$out/$v.bed" "$lab" --outdir "$out" "${HUMAN_DET[@]}"
      else "$VENV_PY" "$SCORE" "$out/$v.bed" "$lab" --outdir "$out" "${MOUSE_REFS[@]}"; fi
    done
    echo "=== $label DONE $(date) ==="; echo "exit=0 $(date)" > "$out/DONE.ok"
  } > "$out/run.log" 2>&1
}

B=$WD/results/benchmark_tools
run_arm $B/pbmc_10k_v3/peakatail_clipseeded        pbmc_clipseeded        "$PBMC" "$HUM_GTF" 91 16 human &
run_arm $B/pbmc_10k_v3/peakatail_clipseeded_ipfilt pbmc_clipseeded_ipfilt "$PBMC" "$HUM_GTF" 91 16 human \
        --ip-filter --genome-fasta "$HUM_FA" --ip-filter-mode filter &
run_arm $B/gse104556/peakatail_clipseeded/mouse1   testis_m1_clipseeded   "$M1"   "$MOU_GTF" 98 12 mouse &
run_arm $B/gse104556/peakatail_clipseeded/mouse2   testis_m2_clipseeded   "$M2"   "$MOU_GTF" 98 12 mouse &
wait
echo "STAGE 2 ALL ARMS FINISHED $(date)"
