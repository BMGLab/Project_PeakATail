# Head-to-head benchmark results (2026-08-19)

Six tools, two datasets, one scoring path (`scripts/benchmark_tools/score_tool.py`), per-dataset
detected-gene denominators, 3-seed gene-body-constrained nulls. Point-mode, strand-matched,
vs curated PolyASite 2.0. **These results are unfavourable to PeakATail and are reported as measured.**

## Human — 10x pbmc_10k_v3 (CellRanger BAM, 11.8k cells; denominator 285,136 sites)

| tool | n called | P@100 | recall_det | **F1** | TES P@100 |
|---|---:|---:|---:|---:|---:|
| scUTRquant\* | 40,519 | 0.793 | 0.178 | **0.290** | 0.693 |
| polyApipe | 120,916 | 0.380 | 0.199 | **0.261** | 0.118 |
| SCAPTURE | 35,759 | 0.652 | 0.118 | **0.199** | 0.386 |
| Sierra | 106,170 | 0.256 | 0.135 | **0.177** | 0.114 |
| scAPAtrap | 787,138 | 0.100 | 0.300 | **0.150** | 0.024 |
| **PeakATail** | 277,111 | 0.118 | 0.156 | **0.134** | 0.041 |

\*annotation-based (fixed UTRome catalog filtered by detection) — near-tautological, not ranked with
de novo tools. Null precision ~0.022 for every arm.

**PeakATail ranks last among de novo tools on this dataset.**

## Mouse — GSE104556 testis (STARsolo BAMs, 2 mice; denominator 126,686 sites)

| tool | n called (m1/m2) | P@100 | recall_det | **F1** | replicate concordance@100 |
|---|---:|---:|---:|---:|---:|
| scUTRquant\* | 35.7k/35.9k | 0.782/0.784 | 0.259 | **0.389** | n/a (catalog) |
| polyApipe | 89.5k/88.0k | 0.401/0.412 | 0.250/0.251 | **0.308/0.312** | 0.49/0.50 |
| **PeakATail** | 45.9k/46.7k | 0.369/0.368 | 0.198/0.201 | **0.258/0.260** | 0.75/0.73 |
| SCAPTURE | 23.6k (m1 only) | 0.694 | 0.155 | **0.253** | m2 invalid (disk) |
| scAPAtrap | 135.5k/151.1k | 0.250/0.239 | 0.245/0.259 | **0.248/0.248** | 0.88/0.79 |
| Sierra | 22.6k/21.4k | 0.509/0.581 | 0.107/0.117 | **0.177/0.195** | 0.79/0.83 |

PeakATail is mid-pack here — second among de novo tools by F1 — but it is **not** the reproducibility
leader. Ranked by replicate concordance@100: scAPAtrap 0.88/0.79 > Sierra 0.79/0.83 > PeakATail
0.75/0.73 > polyApipe 0.49/0.50. scAPAtrap achieves this at *higher* recall than PeakATail
(0.245–0.259 vs 0.198–0.201), so the honest scope of PeakATail's reproducibility claim is **against
polyApipe only**. The qualifier that must travel with scAPAtrap's 0.88: it is measured on a set its
own `reducePeaks(min.cells=10, min.count=10)` has already depth-cleaned (literally zero depth-1
calls), and it buys that reproducibility at precision 0.24–0.25 vs PeakATail's 0.37.

## Why PeakATail underperforms on PBMC — what we tested

1. **Not a threshold/ranking problem.** Ranking its 277k calls by total UMI count and scoring the top
   22k / 36k / 106k gives precision 0.119 / 0.120 / 0.118 — *identical to the full set*. Its internal
   ranking does not enrich for atlas-supported sites at all.
   > **⚠ numbers pending re-verification (harness bug 0g).** The depth used for this ranking was
   > mis-keyed (coordinate-ordered `pas.bed` zipped onto pas_id-ordered matrix rows), so the top-N
   > sets were wrong — the shipped top-22k shares only 7.89% of its members with the correctly-ranked
   > one. Re-scored with correct keying the precisions are **0.120 / 0.118 / 0.118**, i.e. the
   > *conclusion is unchanged and strengthened*. The three digits above are left as published pending
   > sign-off; see the PENDING RE-VERIFICATION block in `manuscript/05_figure_index.md` §8 for the
   > full delta table. `conc_all` replicate concordance (0.73–0.75, below) is depth-free and
   > unaffected. **The re-keyed 0.120/0.118/0.118 are themselves provisional: they are read
   > correctly from a matrix that bug 0a (uncommitted) says is 99.87% mis-keyed at the tool level, so
   > this ranking is not yet attributable to PeakATail. Do not publish either triple until 0a lands.**
2. **Not the cleavage offset.** Relaxing the cutoff to 200 bp lifts it only 0.118 → 0.166, while
   polyApipe goes 0.380 → 0.410. The gap is not a coordinate shift.
3. **Not simply "far from gene ends" — that is partly the dataset.** On the *same* BAM, median distance
   to the nearest annotated TES is PeakATail 14.1 kb, Sierra 14.2 kb, polyApipe 13.1 kb, scAPAtrap
   29.2 kb; only annotation-guided SCAPTURE is close (1.8 kb). Deep 10x data yields abundant internal
   coverage signal for every de novo caller.
4. **The discriminator is evidence type.** Tools that use orthogonal evidence win: polyApipe requires
   non-templated poly(A) soft-clips; SCAPTURE applies a DeepPASS sequence classifier; scUTRquant matches
   a curated catalog. PeakATail (and scAPAtrap) call peaks from coverage shape alone, and on deep data
   that admits many non-PAS peaks.
5. **PeakATail's own TIER label does not flag this**: 98.4% of PBMC calls are TIER_1, essentially the
   same as testis (98.0%) — the tier reflects gene assignment, not 3'-end proximity, and so gives
   false confidence.

## Implications

- **Accuracy leadership is not claimable.** Any framing built on "PeakATail calls PAS more accurately"
  is contradicted by this data and must be dropped.
- **Concrete tool improvement** (the actionable output): add poly(A) soft-clip read evidence and/or a
  sequence-model filter to peak calling — precisely the axis separating the winners here. The existing
  internal-priming machinery is adjacent but tests genomic A-richness, not read-level tails.
- **Dataset sensitivity is itself a finding**: PeakATail's precision spans 0.118 (PBMC/CellRanger),
  0.369 (testis/STARsolo), 0.447 (Laughney/STARsolo). Reporting a single accuracy number for any tool
  in this field is misleading; the benchmark should report per-dataset.
- **Where PeakATail does stand out**: replicate reproducibility *relative to polyApipe* (0.73–0.75 vs
  0.49–0.50, whose F1 lead rides on depth-1 singletons — 62% of its mouse calls are depth-1 and those
  reproduce only 0.31 of the time vs 0.81 for its depth ≥2 calls). It is **not** the reproducibility
  leader overall: scAPAtrap (0.79–0.88) and Sierra (0.79–0.83) both sit above it. The remaining
  unambiguous differentiator is the peak-based clustering novelty, an orthogonal axis no competitor
  offers.

## Reproduce

```bash
export LC_ALL=C
python3 scripts/benchmark_tools/score_tool.py <pas.bed> <label> \
    --detected-atlas results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed   # human
# mouse: add --atlas/--tes/--genome/--genebodies for GRCm38 (see scripts/benchmark_tools/README_STATUS.md)
```
Per-arm TSVs: `results/benchmark_tools/<dataset>/<tool>/score_*.tsv`.

---

## Atlas-independent confirmation (PacBio Kinnex long reads, 2026-08-19)

104 M deduplicated poly(A)-verified FLNC molecules, donor-mismatched so this is **site-level** truth.
Figure: `manuscript/figures/kinnex_truth_validation.*`.

**The precision ordering replicates against orthogonal evidence** (de novo Spearman ρ = 0.90):

| tool | P@100 vs PolyASite | P@100 vs long-read truth | genuine PAS | internal priming | no long-read support |
|---|---:|---:|---:|---:|---:|
| scUTRquant\* | 0.793 | 0.64 | 64% | 5% | 31% |
| SCAPTURE | 0.652 | 0.58 | 58% | 22% | 20% |
| polyApipe | 0.380 | 0.46 | 46% | 7% | 47% |
| Sierra | 0.256 | 0.41 | 41% | **50%** | 9% |
| scAPAtrap | 0.100 | 0.21 | 21% | 26% | 53% |
| **PeakATail** | 0.118 | 0.20 | **20%** | **30%** | **51%** |

\*annotation-based, shown but not ranked with de novo tools.

**Three findings.**
1. **PeakATail's position is confirmed, not an atlas artefact** — last/tied-last among de novo tools under
   truth that owes nothing to PolyASite. Its precision does rise (0.118 → 0.20) because the long-read truth
   set is denser, but its rank does not move.
2. **30% of PeakATail's calls are internal-priming artefacts** — now measured directly (long-read-supported
   3' ends with genomic A-rich downstream), not inferred. `--ip-filter` is off by default. This is a second,
   independent, cheap lever alongside poly(A)-clip seeding.
3. **Sierra's atlas precision was half internal priming** (50% of its calls) — those artefacts are present in
   the atlas, so atlas-only benchmarking flattered it. A caution for the field, not just for us.

**Methodological note:** F1 rank is *not* stable here (panel D) — a fixed call set scored against truth sets
differing 53-fold in size reorders every tool except PeakATail. Precision and recall are the reportable
quantities; quote F1 only with its truth set fixed and named.
