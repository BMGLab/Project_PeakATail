# PeakATail manuscript — planning hub

Goal: **Q1 methods + application paper** on PeakATail (`ema`), the single-cell PAS/APA
toolkit in `tools/PeakATail` (branch `biolab-manuscript`, v0.2.0). Three pillars:

1. **Prove it's right** — PAS calls validated against ground truth; sensitivity/specificity
   quantified; FDR calibration of switch tests (→ `02_validation_plan.md`).
2. **Prove it's competitive** — head-to-head benchmarks vs Sierra, scAPAtrap, SCAPTURE,
   polyApipe, scUTRquant (panel follows the NAR 2026 benchmark, gkag490) (→ `03_benchmarking_plan.md`).
3. **Show it matters** — recover known cell-type APA biology, then novel insight in the
   Laughney LUAD cohort, tumor-associated macrophages, PBMC (→ `04_datasets.md`).

## Documents

| File | Contents |
|---|---|
| `01_outline_and_journals.md` | Ranked journals (Genome Biology → Genome Research → BiB), title candidates, abstract draft, figure-by-figure plan (6 figs), section outline with per-claim evidence, honest novelty positioning |
| `02_validation_plan.md` | 4-tier ground truth (matched long-read Kinnex PBMC > bulk 3'-seq > atlases > simulation), positive-control gene checklist, negative controls / internal-priming specificity, metric definitions & pass/fail bars |
| `03_benchmarking_plan.md` | Survey table of ~18 scAPA tools + maintenance status, published benchmarks, recommended head-to-head panel, concrete protocol (precision/recall vs PolyASite/PolyA_DB at distance windows, downsampling, runtime/memory, diff-APA concordance) |
| `04_datasets.md` | Verified accessions: validation (GSE130708 ScNaUmi-seq, FLAMES scmixology, Kinnex PBMC), benchmarking (pbmc_10k_v3, pbmc4k/8k, GSE104556 testis, scAPA's exact inputs), applications (GSE123904 in hand, LuCA components with FASTQ, GBM GSE163120, HCC GSE156625) — with chemistry/access flags |

## What already exists (don't redo)

- **RERUN_2026-08_fixed** (`results/laughney_rerun_2026-08_fixed` → ssd2): full 17-sample
  Laughney cohort run (lambda_gradient), 5-arm strategy grid **with
  `benchmark_vs_polyasite_v3.json` per arm** (e.g. lg_annotate precision@100bp ≈ 0.999),
  13 reannotate branches, GEX-vs-PAS concordance (ARI/AMI per sample + Sankey inputs).
- **VALIDATION_de26d16** (ssd2): 11/11 PASS audit of the barcode-doubling fix, strategy
  sanity runs, determinism check.
- July analysis assets on ssd2 (`bench_summit_vs_atlas.py`, `crosscheck_diff_vs_length.py`,
  `sankey_gex_vs_pas.py`, `recover_switch_direction.py`, …) — reusable, but hardcode the
  old July run paths.
- Technical-report figures in `tools/PeakATail/docs/assets/figures/` (peak precision,
  dual-DB validation, clustering ARI/AMI sweeps) can seed manuscript panels.
- `tools/PeakATail/ROADMAP.md` already carries the publication plan + validation bars
  (precision ≥70%, recall ≥60%, F1 ≥0.65; strong: 80/70).

## Immediate next steps

1. **Unblock the cell-type switches:** RERUN's 24 SWITCH_CELLTYPE tasks all failed on a
   `--pasbed null/pasbed.bed` interpolation bug in the pipeline's `main.nf` (ssd2, Amir's
   area). One-line fix + `nextflow -resume` reruns only those ~24 short tasks. This yields
   the per-celltype `switch diff/length/trend` tables — the core biological result.
2. Verify `--ip-filter` wiring live (docs contradict themselves; HANDOFF says D9 landed) —
   never claim it in Methods before confirming.
3. Download the PacBio Kinnex PBMC 3' set + pbmc_10k_v3 (BAM ready-made) to start Tier-1
   validation and the benchmark panel.
4. Positive controls first pass: spermatogenesis GSE104556 (strongest known APA gradient)
   + T-cell activation genes from the checklist in `02_validation_plan.md`.
