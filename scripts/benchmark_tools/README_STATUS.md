# Benchmark tools — consolidated status

Consolidated: 2026-08-18 23:25 +03, from `status/*.md` (per-tool install reports) plus a live
check of the BAM download. Details, exact install commands, and full invocations live in the
per-tool files under `status/`.

## Tool status table

| Tool | State | Version | Smoke test | Blocking issue | Next step |
|---|---|---|---|---|---|
| PeakATail (reference) | Ready via workaround invocation | 0.2.0 (branch `biolab-manuscript` @ `18678ef`) | PASS — `.venv/bin/python -c 'from ema.cli import main; main()'` gives `--help` / `--version` / `run --help` | Bundled `.venv` is Python 3.10 but pyproject needs >=3.11, so `pip install -e .` fails and `.venv/bin/ema` is a stale entry point that crashes | Run via `run_peakatail.sh` (uses the workaround) when BAM lands; optionally build `.venv311` for the proper fix |
| polyApipe | READY | 0.1.0 (py 3.13.9, pysam 0.24.0; conda `bench_polyapipe`: featureCounts 2.1.1, umi_tools 1.1.6) | PASS — full end-to-end run on bundled demo BAM (peaks GFF + counts + polyA BAM, exit 0) | None (`pip install .` of repo broken; script installed directly — done. R polyApiper optional, skipped) | Launch on the BAM (`--umi_tag UB`) as soon as it verifies |
| Sierra | READY | 0.99.27 (R 4.5.3, regtools 1.0.0; conda `bench_sierra`) | PASS — library loads, `FindPeaks` exported; regtools verified on test BAM (329,972 junctions) | None; must use Ensembl-named GTF (`Homo_sapiens.GRCh38.99.gtf`) to match the chr-less BAM. Hours-scale runtime | Start FIRST when BAM lands (slowest step): `regtools junctions extract` then `FindPeaks`/`CountPeaks` |
| scAPAtrap | READY after 1-line param fix | 0.2.0 (R 4.5.3, samtools 1.24, featureCounts 2.1.1, umi_tools 1.1.6; conda `bench_scapatrap`) | PASS — package loads, `setTools(check=TRUE)` validates all 4 binaries | Example config in its status file sets `tp$chrs` to `chr1`-style — WRONG for this BAM (bare Ensembl names); `readlength` must be 91 | Set `tp$chrs <- c(1:22,"X","Y")`, `readlength=91`, then run when BAM lands (outputDir must not pre-exist) |
| SCAPTURE | Installed; annotation prebuild pending | v1.0 (conda `bench_scapture`: py 3.7.8, R 3.6.3, TF 2.0.0, featureCounts 1.6.4 added) | PASS — usage prints; DeepPASS end-to-end on CPU wrote `Predict_Result.txt` | GTF needs `gene_biotype`→`gene_type` sed (untested); CUDA 10.0 predates the RTX A4000, so DeepPASS must run with `CUDA_VISIBLE_DEVICES=""` | Build annotation NOW (no BAM needed): sed the GTF + `scapture -m annotation`; then `PAScall -l 91` when BAM lands |
| scUTRquant | PARTIAL — orchestrator only | pipeline v0.5.1 (snakemake 7.32.4, py 3.11.15; conda `bench_scutrquant`) | PASS — dry run on the real pbmc_10k_v3 BAM config builds the full 14-job DAG, zero errors | Per-rule conda envs (~GBs, incl. patched `kallisto=0.46.2sq`) and 358 MB UTRome index not yet fetched/built | Run `snakemake --conda-create-envs-only` NOW in <=440 s chunks (no BAM needed); then full run when BAM lands |
| scTail | Installed; NOT usable for de novo PAS on pbmc_10k_v3 | 0.1.8 (venv `bench_envs/sctail`, py 3.10, torch CPU-only) | PASS — all 3 subcommands print usage; import + CPU tensor OK | Needs long (>100 bp) read-1 alignments; pbmc_10k_v3 R1 is 28 bp (CB+UMI only), so PAS *detection* is impossible on this dataset | Decide: run `scTail-count` only (quantification vs an external PAS BED), and/or add a long-R1 public dataset for detection |

## Data-download state (10x pbmc_10k_v3)

Target dir: `/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3`

| File | Status |
|---|---|
| BAM (44,198,457,238 B) | IN PROGRESS — 57% (25.3 GB) at 23:22, wget PID 809842 alive, ~16 min ETA; resumable with the same `wget -c` |
| BAM index (.bai) | COMPLETE |
| Filtered matrix tar.gz | COMPLETE (`gzip -t` OK; contains barcodes.tsv.gz needed by several tools) |
| Web summary | COMPLETE |

- No publisher md5 exists; verify by exact size (`stat -c%s` == 44198457238) + `samtools quickcheck -v`.
- BAM facts (verified from the partial file): **bare Ensembl chromosome names** (`1..22, X, Y, MT` — no `chr` prefix), CellRanger 3.0.0 / GRCh38-3.0.0 (Ensembl 93), R2 = 91 bp, `CB` on ~98% of reads. Matches the on-disk Ensembl GTF and PolyASite directly; `normalize_chroms.py` is only for tools emitting UCSC names.
- 11,769 filtered cells; `export LC_ALL=C` before any sort (tr_TR locale breaks sorting).
- Dataset 2 (GSE104556 mouse testis): fetch script written (`fetch_gse104556.sh`), NOT run; needs CellRanger + mm10.

## Prioritized next actions for the head-to-head

**Now, while the BAM finishes (~16 min) — none of these need the BAM:**
1. scUTRquant: prebuild the rule conda envs (`--conda-create-envs-only`, chunked) — the only tool whose install is still materially incomplete.
2. SCAPTURE: sed the GTF (`gene_biotype`→`gene_type`) and run `scapture -m annotation` — validates the one untested step off the critical path.
3. Extract `barcodes.tsv` from the filtered-matrix tarball; make both variants (with `-1` for Sierra/SCAPTURE, stripped for scAPAtrap).
4. Fix the scAPAtrap run script: `tp$chrs` = bare Ensembl names, `tp$readlength = 91`.
5. (Optional) PeakATail proper fix: build `.venv311` — not required, the workaround invocation works.

**The moment the BAM lands** (after `stat -c%s` == 44198457238 and `samtools quickcheck -v`):
1. **Sierra** — start first, it is the slowest (hours): `regtools junctions extract`, then `FindPeaks` with `ncores`.
2. **PeakATail** — `run_peakatail.sh` (its preflight refuses an incomplete BAM; no `--atlas` to avoid circularity).
3. **polyApipe** — ready as-is; `--cell_barcode_tag CB --umi_tag UB`.
4. **scAPAtrap** — after action 3+4 above.
5. **SCAPTURE PAScall** (`-l 91`, `CUDA_VISIBLE_DEVICES=""`) — after its annotation build; then PASquant.
6. **scUTRquant** full snakemake run — after its envs are built (its own DAG fetches UTRome + whitelists).

**Ready-to-run summary**: polyApipe, Sierra, and PeakATail (workaround) can launch the moment the
BAM verifies; scAPAtrap needs only the 1-line chrs/readlength fix; SCAPTURE needs its annotation
prebuild; scUTRquant needs its env prebuild; scTail cannot do de novo detection on this dataset at all.

**Later / decisions:**
- scTail: either include as quantification-only (scTail-count vs a common PAS BED) with an explicit caveat, or add a long-R1 dataset (scTail paper's GSE deposits) for a fair detection comparison.
- Dataset 2: run `fetch_gse104556.sh` + CellRanger/mm10 once dataset 1 head-to-head is underway.
- Scoring: `evaluate_vs_reference.py` + `normalize_chroms.py` are written; remember PolyASite is bare-Ensembl like the BAM, so only UCSC-style tool outputs need translation.

## Evaluation template + first head-to-head (2026-08-19)

polyApipe FINISHED on pbmc_10k_v3 (3:26:37 wall, 13.09 GiB peak RSS, exit 0) and is the
first externally scored tool. The per-tool evaluation path is now standardized:

1. **Standard form**: `results/benchmark_tools/pbmc_10k_v3/<tool>/pas.bed` — BED6, one 1-bp
   POINT per inferred PAS, `#` header documents the geometry derivation. polyApipe converter:
   `scripts/benchmark_tools/convert_polyapipe.py` (PAS = GFF **end** for `+`/`f`, GFF **start**
   for `-`/`r`; verified against polyApipe source and asserted per-feature against the peak
   name). Variants: `pas.bed` = misprime-excluded DEFAULT (121,013), `pas_all.bed` (211,844).
2. **Scoring**: `scripts/benchmark_tools/score_tool.py <pas.bed> <label>` — machinery adapted
   verbatim from `scripts/manuscript_figures/benchmark_curated.py` (point mode, strand-matched
   closest, cutoffs 10/25/50/100/200, 3-seed genic-shuffle null, recall vs full atlas AND vs
   the shared 285,220-site detected-gene reference in `results/benchmark_tools/shared_refs/`).
   VALIDATED: rerunning it on `lg_annotate/pasbed.bed` reproduces every benchmark_curated.tsv
   precision/recall/F1 row (incl. null seeds) exactly.
3. **Running comparison**: `scripts/manuscript_figures/benchmark_tools_running.py` regenerates
   `results/figures/manuscript/benchmark_tools_running.{tsv,png}` (+ PNG/PDF in
   `manuscript/figures/`). Re-run after each new tool's `score_<label>.tsv` lands (add the
   label to `TOOLS`/`RES` in the script).

First numbers @100 bp (PeakATail rows PROVISIONAL — Laughney cohort, pbmc rerun pending):

| tool | n scored | precision | recall (detected-gene) | null prec (mean) |
|---|---|---|---|---|
| polyApipe (misprime excl.) | 120,916 | 0.380 | 0.197 | 0.022 |
| polyApipe (all peaks) | 211,706 | 0.277 | 0.242 | 0.022 |
| PeakATail lg_annotate* | 22,629 | 0.447 | 0.066 | 0.022 |
| PeakATail B1_cohort_full* | 16,497 | 0.475 | 0.050 | 0.022 |

## GSE104556 mouse testis — consolidated scoring round (2026-08-19)

Four tools finished and scored on both biological replicates (Mouse1/Mouse2 STARsolo BAMs,
GRCm38): **PeakATail** (`results/benchmark_tools/gse104556/peakatail/`), **Sierra** 0.99.27
(`.../sierra/`, geometry-proof conversion via `make_pas_bed.sh`), **scUTRquant** v0.5.1
(`.../scutrquant/`, target `utrome_mm10_v2`; STARsolo→CellRanger-tag BAM prep required because
patched kallisto segfaults on raw STARsolo BAMs), and **polyApipe** 0.1.0 (`.../polyapipe/`,
finished 09:38 +03 2026-08-19; converted with `convert_polyapipe.py` — same proven geometry as
pbmc, per-feature name/coordinate assertions all passed; default = misprime-excluded `pas.bed`,
`pas_all.bed` scored as the `_all` arms). scAPAtrap/SCAPTURE not finished — score with the SAME
command when they land.

**Detected-gene restricted reference (denominator warning RESOLVED).** Built per
`data/references/atlases/README.md` §6 into `results/benchmark_tools/gse104556/shared_refs/`
(README + md5s there): detected genes = UNION of gene IDs in the two PeakATail runs'
`annotatedpas.bed` = **14,014** ENSMUSG ids (mouse1 13,406, mouse2 13,332; all present exactly
once in `gene_end.GRCm38.102.bed`); restricted atlas = strand-matched PolyASite2 mouse rep
sites inside those gene bodies = **126,686** sites (42.1% of the 301,006-site atlas) =
`pas2.in_detected_genes.bed`. Every mouse arm below now has `recall`/`f1`/`roadmap` rows vs
`atlas_detected` on this ONE shared denominator — mouse recall flavors now mirror the pbmc
figure and are directly comparable across tools.

**Scoring command** (all eleven arms; the `--detected-atlas` guard in `score_tool.py` demands
the full explicit mouse flag set):

    python3 scripts/benchmark_tools/score_tool.py <pas.bed> <label> --outdir <tooldir> \
      --atlas data/references/atlases/polyasite2.GRCm38.96.rep_sites.bed6 \
      --tes data/references/atlases/tes.protein_coding.GRCm38.102.bed6 \
      --genome data/references/mouse/chrom.sizes.filt \
      --genebodies results/benchmark_tools/gse104556/shared_refs/genebodies.merged.bed \
      --detected-atlas results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed

Sierra re-scored with the new flag: all pre-existing rows (precision, null, full/TES recall,
F1, roadmap, meta) **byte-identical** to the previously committed `score_sierra_*.tsv`; each
arm only GAINED 25 `atlas_detected` rows.

### Mouse table (tool × mouse; @100 bp unless noted; null = 3-seed genic-shuffle mean)

| tool arm | n called (scored) | P@50 | P@100 | null P@100 | fold/null | R@100 full (301,006) | R@100 restricted (126,686) | R@100 TES | F1@100 restr. | wall / peak RSS | concordance ≤100 bp |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PeakATail mouse1 | 45,921 (45,907) | 0.226 | 0.369 | 0.0139 | 27× | 0.093 | **0.198** | 0.255 | **0.258** | 2:06:48 / 14.0 GiB | 0.746 (m1→m2) |
| PeakATail mouse2 | 46,672 (46,654) | 0.227 | 0.368 | 0.0143 | 26× | 0.095 | **0.201** | 0.257 | **0.260** | 1:26:46 / 13.6 GiB | 0.734 (m2→m1) |
| Sierra mouse1 | 22,583 (22,550) | 0.355 | 0.509 | 0.0144 | 35× | 0.057 | 0.107 | 0.170 | 0.177 | 1:12:48 / 4.39 GiB | 0.790 (m1→m2) |
| Sierra mouse2 | 21,365 (21,335) | 0.396 | 0.581 | 0.0142 | 41× | 0.062 | 0.117 | 0.181 | 0.195 | 57:50 / 4.36 GiB | 0.834 (m2→m1) |
| Sierra pooled | 37,863 (37,814) | 0.315 | 0.505 | 0.0139 | 36× | 0.068 | 0.130 | 0.193 | 0.207 | — | — |
| scUTRquant mouse1 † | 35,745 (35,732) | 0.757 | 0.782 | 0.0138 | 57× | 0.151 | 0.259 | 0.603 | 0.389 | 11:26 (+1:04:32 BAM prep) / 4.0 GiB | n/a † |
| scUTRquant mouse2 † | 35,858 (35,845) | 0.758 | 0.784 | 0.0142 | 55× | 0.152 | 0.259 | 0.604 | 0.389 | 10:16 (+53:34 BAM prep) / 4.0 GiB | n/a † |
| polyApipe mouse1 | 89,588 (89,527) | 0.385 | 0.401 | 0.0142 | 28× | 0.168 | **0.250** | 0.294 | **0.308** | 1:18:36 / 6.77 GiB | 0.494 (m1→m2) |
| polyApipe mouse2 | 88,037 (87,983) | 0.396 | 0.412 | 0.0145 | 28× | 0.170 | **0.251** | 0.295 | **0.312** | 1:02:59 / 6.76 GiB | 0.503 (m2→m1) |
| polyApipe mouse1 all ‡ | 100,049 (99,984) | 0.364 | 0.380 | 0.0141 | 27× | 0.179 | 0.267 | 0.311 | 0.313 | (same run) | — |
| polyApipe mouse2 all ‡ | 96,929 (96,871) | 0.378 | 0.394 | 0.0146 | 27× | 0.179 | 0.265 | 0.311 | 0.317 | (same run) | — |

† scUTRquant is ANNOTATION-BASED: its sites are the fixed mm10 UTRome catalog filtered by
detection (`mouse{1,2}/pas.bed` headers carry the full caveat +
`scripts/benchmark_tools/scutrquant/gse104556/extract_pas_bed_mouse.py` derivation, identical
to the pbmc one). Precision vs an annotation-derived atlas — and cross-replicate concordance
of a fixed catalog — are near-tautological; footnote in any figure. Its recall is also
bounded by the catalog. Runtime is dominated by the STARsolo→CR-tag BAM prep this benchmark
imposed, not the pipeline itself.

‡ polyApipe `_all` arms = `pas_all.bed` incl. misprime="True" peaks (sensitivity variant);
the misprime-excluded `pas.bed` is the tool's DEFAULT and the arm to quote. Concordance was
computed on the default variant only.

### Headline read (identical denominators, restricted recall)

- **PeakATail vs Sierra** (both de novo): Sierra is the precision arm (0.51–0.58 vs 0.37 @100)
  but calls half as many sites and reaches barely half PeakATail's restricted recall
  (0.107–0.117 vs 0.198–0.201) — **PeakATail wins F1 on the shared denominator**
  (0.258/0.260 vs 0.177/0.195; Sierra pooled 0.207 still below either PeakATail mouse).
- **scUTRquant** posts the best raw numbers (P 0.78, R_restr 0.259, F1 0.389) but they are
  catalog-tautological (see †) — report separated from the de novo tools.
- **polyApipe** (de novo) lands ABOVE both on F1@100 restricted: **0.308/0.312**
  (vs PeakATail 0.258/0.260, Sierra 0.177/0.195; `_all` variant 0.313/0.317). It pairs
  PeakATail-like precision (0.40–0.41 vs 0.37 @100) with the best de novo restricted recall
  (0.250/0.251 vs 0.198/0.201) on a call set ~2× PeakATail / ~4× Sierra (88–90k scored). The
  trade-off shows up in replicate concordance: only ~0.49–0.50 of calls reproduce within
  100 bp (vs 0.73–0.75 PeakATail, 0.79–0.83 Sierra) — 62–63% of its misprime-excluded peaks
  are depth-1 singletons.
- All real precisions are 26–57× the genic-shuffle null (~0.014 @100) — mouse panel is
  signal-dominated, same conclusion as pbmc.

### Replicate concordance (PeakATail vs Sierra vs polyApipe, side by side)

Same method all tools (strand-matched `bedtools closest -s -d -t first` between the
genome-filtered point sets; `<tool>/replicate_concordance.tsv`; polyApipe on its
misprime-excluded default):

| tool | direction | n | ≤10 | ≤25 | ≤50 | ≤100 | ≤200 | exact |
|---|---|---|---|---|---|---|---|---|
| PeakATail | m1→m2 | 45,907 | 0.322 | 0.489 | 0.631 | 0.746 | 0.769 | 0.078 |
| PeakATail | m2→m1 | 46,654 | 0.317 | 0.481 | 0.620 | 0.734 | 0.757 | — |
| Sierra | m1→m2 | 22,550 | 0.442 | 0.626 | 0.747 | 0.790 | 0.811 | 0.269 |
| Sierra | m2→m1 | 21,335 | 0.467 | 0.661 | 0.788 | 0.834 | 0.854 | — |
| polyApipe | m1→m2 | 89,527 | 0.436 | 0.461 | 0.475 | 0.494 | 0.517 | 0.255 |
| polyApipe | m2→m1 | 87,983 | 0.444 | 0.469 | 0.483 | 0.503 | 0.526 | — |

Read: ~73–75% of PeakATail calls reproduce within 100 bp across biological replicates vs
~79–83% for Sierra — but on a call set twice the size (45.9k vs 22.6k); Sierra's higher
per-call concordance mirrors its higher precision / lower recall trade-off. PeakATail's low
exact-point fraction (7.8% vs Sierra 26.9%) reflects read-density peak summits vs Sierra's
fitted-interval ends. NOTE the PeakATail matched counts are identical in both directions at
every cutoff — verified real, not a bug: every ≤200 bp cross-replicate match is a MUTUAL
nearest-neighbor pair (pair sets byte-identical both ways), so the cumulative histograms
coincide; only the denominators differ. polyApipe shows the same mutual-NN signature at
cutoffs ≤100 (matched counts identical both directions; they differ by 5 at ≤200) and a
distinctive flat curve: 0.44 already at ≤10 bp (exact-coordinate peak calls, 25.5% exact)
but only ~0.49–0.50 by ≤100 — the half of its calls that do reproduce match near-exactly,
while the depth-1 singleton half (62–63% of peaks) mostly has no counterpart in the other
replicate.

### Remaining for the mouse panel

- scAPAtrap (mouse1 mid-flight), SCAPTURE — score with the same command + add rows to the
  table above.
- `benchmark_tools_running.py` still has no dataset facet (hardcoded pbmc tool list); add a
  `dataset` facet and pull in the eleven `score_*.tsv` above in the consolidated figure pass.

**Next-analysis note (consolidation pass):** depth-stratified reproducibility panel — polyApipe's
F1 lead rides on depth-1 singletons (62-63% of its calls; concordance 0.31 vs 0.81 for depth>1).
Plot per-tool concordance-vs-depth and F1-vs-reproducibility; report tool DEFAULT arms as primary,
post-hoc depth filters only as a labeled sensitivity view.

## Benchmark matrix COMPLETE (2026-08-20)

SCAPTURE mouse2 rescued: PAScall succeeded on the ssd0 rerun (24,076 PAS points; the earlier
"failure" was PASquant-only, not needed for site scoring). Final mouse de novo F1@100:
polyApipe 0.308/0.312 > PeakATail 0.258/0.260 > SCAPTURE 0.253/0.241 ≈ scAPAtrap 0.248/0.248 >
Sierra 0.177/0.195. SCAPTURE is the mouse precision arm (0.694/0.672) as on human. SCAPTURE
replicate concordance 0.700/0.688 @100bp (22.2% exact) — between Sierra (0.79-0.83) and
PeakATail (0.73-0.75). All 6 tools × 2 datasets now have verified scores; nothing pending.
(scapture/replicate_concordance.tsv column order differs slightly from the sierra/peakatail
files — self-describing header, harvest by name not position.)

## Stage 2 — clip-seeded caller, first arm (testis mouse1, 2026-08-20)

| arm | n | P@100 | R_det | F1@100 | wall | RSS |
|---|---:|---:|---:|---:|---:|---:|
| shipped lambda_gradient | 45,921 | 0.369 | 0.198 | 0.258 | 2:06:48 | 14.0 GiB |
| **clip_seeded, both tiers** | 58,848 | 0.399 | 0.236 | **0.296** | **43:39** | **2.0 GiB** |
| clip_seeded tier-1 (38,690) | | 0.544 | 0.212 | 0.305 | | |
| clip_seeded tier-2 (20,158) | | 0.119 | 0.024 | 0.040 | | |
| *polyApipe (ref.)* | *89,527* | *0.401* | *0.250* | *0.308* | | |

Same denominator (126,686 sites), null ~0.015. Combined F1 +15% over shipped; tier-1 at parity with
polyApipe's F1 at markedly higher precision. PBMC arms (the gate) still running.
