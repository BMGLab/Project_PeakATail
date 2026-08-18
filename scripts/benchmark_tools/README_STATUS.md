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
