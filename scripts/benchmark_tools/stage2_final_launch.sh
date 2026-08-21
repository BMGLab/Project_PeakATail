#!/bin/bash
# FINAL Stage-2 re-run (after Stage 1b+1c). Code tree is SNAPSHOTTED at launch so in-flight edits cannot leak.
# Defaults are PRE-SPECIFIED here, before any number is seen:
#   * default output = clip_seeded at --polya-min-umis 1 (tool defaults), both tiers emitted, tiers scored separately
#   * secondary, LABELLED POST-HOC SENSITIVITY arm = tier-1 subset with >=2 distinct molecules (col5 >= 2), scored alongside
#   * PRE-REGISTERED precision-first default (PI, 2026-08-21, BEFORE this run): ipfilt arm, tier-1, >=2 molecules;
#     gate for it: P@100 >= 0.50 on pbmc AND both testis mice (F1 reported, not gated)
# Gate (unchanged, two-sided): F1@100 > 0.261 AND P@100 >= 0.38 on pbmc_10k_v3 default output. Tiers reported separately.
set -uo pipefail
export LC_ALL=C PEAKATAIL_NO_TIMESTAMP=1
WD=/mnt/ssd1/Projects/PeakATail_wd
SRC=$WD/tools/pa-polya
COMMIT=$(git -C "$SRC" rev-parse HEAD)
SNAP=$WD/tools/pa-polya-run-${COMMIT:0:8}
[ -d "$SNAP" ] || git -C "$WD/tools/PeakATail" worktree add --detach "$SNAP" "$COMMIT" >/dev/null
export PYTHONPATH=$SNAP
VENV_PY=$WD/tools/PeakATail/.venv/bin/python
NORM=$WD/scripts/benchmark_tools/normalize_chroms.py
SCORE=$WD/scripts/benchmark_tools/score_tool.py
HUM_GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
HUM_FA=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa
MOU_GTF=$WD/data/references/mouse/Mus_musculus.GRCm38.102.gtf
MOU_FA=$WD/data/references/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa
PBMC=$WD/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
M1=$WD/data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam
M2=$WD/data/benchmark/gse104556/starsolo/Mouse2_scRNAseq/Aligned.sortedByCoord.out.bam
R=$WD/data/references; G=$WD/results/benchmark_tools/gse104556
MOUSE_REFS=(--atlas $R/atlases/polyasite2.GRCm38.96.rep_sites.bed6 --tes $R/atlases/tes.protein_coding.GRCm38.102.bed6
            --genome $R/mouse/chrom.sizes.filt --genebodies $G/shared_refs/genebodies.merged.bed
            --detected-atlas $G/shared_refs/pas2.in_detected_genes.bed)
HUMAN_DET=(--detected-atlas $WD/results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed)

run_arm() {
  local out=$1 label=$2 bam=$3 gtf=$4 seqlen=$5 threads=$6 species=$7; shift 7
  mkdir -p "$out"; local run=$out/run; rm -rf "$run"
  {
    echo "=== $label start $(date) | code commit $COMMIT (frozen snapshot $SNAP) | extra: $* ==="
    "$VENV_PY" -c "import ema.strategies.clip_seeded as m; assert '$SNAP' in m.__file__, m.__file__; print('code:', m.__file__)"
    /usr/bin/time -v -o "$out/runtime_mem.txt" \
      "$VENV_PY" -c 'from ema.cli import main; main()' run \
        --bam-dir "$bam" --gtf "$gtf" --output "$run" --threads "$threads" \
        --barcode-tag CB --cb-len 16 --seq-len "$seqlen" --ignore-chro 'MT' \
        --plot-engine none --no-progress --peak-strategy clip_seeded "$@" > "$out/log" 2>&1
    local rc=$?; echo "=== $label ema exit $rc $(date) ==="
    [[ $rc -eq 0 && -s $run/pasbed.bed ]] || { echo "exit=$rc" > "$out/FAILED.err"; return 1; }
    awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}' "$run/pasbed.bed" \
      | "$VENV_PY" "$NORM" --to ensembl --unmapped drop | sort -k1,1 -k2,2n > "$out/pas.bed"
    awk -F'\t' '$5>0'  "$out/pas.bed" > "$out/pas_tier1.bed"
    awk -F'\t' '$5==0' "$out/pas.bed" > "$out/pas_tier2.bed"
    awk -F'\t' '$5>=2' "$out/pas.bed" > "$out/pas_tier1_ge2umi_POSTHOC.bed"
    # PRE-REGISTERED precision-first default (PI decision 2026-08-21, manuscript/13_reliability_positioning.md):
    # tier-1 only, >=2 distinct molecules; meaningful on the --ip-filter arm (IP-flagged sites already dropped there)
    # Only an --ip-filter arm can produce the pre-registered default; non-IP arms get a POSTHOC-labelled file instead.
    if [[ " $* " == *" --ip-filter "* ]]; then awk -F'\t' '$5>=2' "$out/pas.bed" > "$out/pas_PRESPEC_precision_default.bed"
    else awk -F'\t' '$5>=2' "$out/pas.bed" > "$out/pas_tier1_ge2mol_noIP_POSTHOC.bed"; fi
    [[ -s $run/annotated_matrix.mtx ]] && cp "$run/annotated_matrix.mtx" "$out/counts.mtx"
    echo "pas.bed $(wc -l < $out/pas.bed) | tier1 $(wc -l < $out/pas_tier1.bed) | tier2 $(wc -l < $out/pas_tier2.bed) | tier1>=2umi $(wc -l < $out/pas_tier1_ge2umi_POSTHOC.bed)"
    for v in pas pas_tier1 pas_tier2 pas_tier1_ge2umi_POSTHOC pas_PRESPEC_precision_default pas_tier1_ge2mol_noIP_POSTHOC; do
      [[ -s $out/$v.bed ]] || continue
      local lab=${label}_${v#pas}; lab=${lab%_}
      if [[ $species == human ]]; then "$VENV_PY" "$SCORE" "$out/$v.bed" "$lab" --outdir "$out" "${HUMAN_DET[@]}"
      else "$VENV_PY" "$SCORE" "$out/$v.bed" "$lab" --outdir "$out" "${MOUSE_REFS[@]}"; fi
    done
    echo "=== $label DONE $(date) ==="; echo "exit=0 commit=$COMMIT $(date)" > "$out/DONE.ok"
  } > "$out/run.log" 2>&1
}
B=$WD/results/benchmark_tools
run_arm $B/pbmc_10k_v3/peakatail_clipseeded_final        pbmc_final        "$PBMC" "$HUM_GTF" 91 16 human &
run_arm $B/pbmc_10k_v3/peakatail_clipseeded_final_ipfilt pbmc_final_ipfilt "$PBMC" "$HUM_GTF" 91 16 human \
        --ip-filter --genome-fasta "$HUM_FA" --ip-filter-mode filter &
# mouse arms run WITH the IP filter: the pre-registered precision default (13 §1) is tier-1 ∩ IP-pass ∩ >=2 molecules
# and is gated on BOTH mice; the no-IP mouse baseline already exists (peakatail_clipseeded_v3/mouse1, Stage-1c acceptance).
run_arm $B/gse104556/peakatail_clipseeded_final/mouse1   testis_m1_final   "$M1"   "$MOU_GTF" 98 12 mouse \
        --ip-filter --genome-fasta "$MOU_FA" --ip-filter-mode filter &
run_arm $B/gse104556/peakatail_clipseeded_final/mouse2   testis_m2_final   "$M2"   "$MOU_GTF" 98 12 mouse \
        --ip-filter --genome-fasta "$MOU_FA" --ip-filter-mode filter &
wait
echo "STAGE 2 FINAL: ALL ARMS FINISHED $(date) (code $COMMIT)"
