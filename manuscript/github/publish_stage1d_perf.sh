#!/bin/bash
# Publish the Stage 1d-B PR (memory / CPU). Review, then run:  bash manuscript/github/publish_stage1d_perf.sh
set -euo pipefail
cd "$(dirname "$0")"
R="BMGLab/PeakATail"; T=/mnt/ssd1/Projects/PeakATail_wd/tools/PeakATail
footer=$'\n\n---\n*Prepared by Claude (AI assistant) for the PeakATail manuscript effort (PI: Ebru Kocakaya); adversarially verified (FIXED → ready; outputs byte-identical on full PBMC, mouse1 and slice) before staging; posted via the gh account shown by `gh auth status`.*\n\n🤖 Generated with [Claude Code](https://claude.com/claude-code)'
git -C "$T" push -u origin perf/clip-memory:perf/clip-memory
PR=$(gh pr create -R $R -B develop -H perf/clip-memory \
  -t "perf: sparse TF-IDF + per-chromosome parallel peak calling — PBMC 10k 293.7 GB → 12.5 GB, 3h46m → 28 min, outputs byte-identical (Stage 1d-B)" \
  -b "$(cat pr_stage1d_memory_cpu.md)$footer")
gh pr edit -R $R perf/clip-memory --add-reviewer TRextabat || true
echo "- PR Stage 1d-B (memory/CPU): $PR" >> PUBLISHED.md
gh issue comment -R $R 95 -b "PR for §2 (memory / CPU): $PR — the 280 GiB was the dense TF-IDF in clustering (4 float64 copies of cells × PAS), not the caller; sparse TF-IDF is bit-identical; peak calling now runs per (contig, strand) with a deterministic merge. PBMC 10k: 293.7 GB → 12.5 GB, 3h46m → 28 min; all 18 output files byte-identical to the final benchmark run.$footer" || true
echo "STAGE 1d-B PUBLISHED: $PR"
