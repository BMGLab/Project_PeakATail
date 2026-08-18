# benchmark_tools — head-to-head PAS/APA caller benchmark harness

Protocol for benchmarking **PeakATail** against competitor single-cell PAS/APA
callers on shared data. Written 2026-08-18. Nothing in here runs a tool by
itself; each tool gets a wrapper script (`run_<tool>.sh`) that honours the
contract below, and `evaluate_vs_reference.py` scores every tool identically.

## 0. Machine rules (read first)

- `export LC_ALL=C` before ANY `sort` — this machine's `tr_TR` locale corrupts
  lexicographic chromosome order. Every script in this directory sets it itself.
- `/mnt/ssd2` is **read-only**. All outputs go under
  `/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/`.
- Tool isolation: R packages -> `~/R/bench-lib` (`.libPaths`), python tools ->
  a conda env `bench_<tool>` (conda at `/home/biolab/miniconda3/bin/conda`; no
  mamba) or a venv under
  `/mnt/ssd1/Projects/PeakATail_wd/tools/bench_envs/<tool>/`.
- No single command may run > 8 minutes; chunk installs, use `timeout`.
- Every tool gets an honest status file: `status/<tool>.md` — install commands
  that worked, version, smoke-test evidence, input/output formats, and the
  exact invocation for the dataset BAM. If an install stalls, record exactly
  how far it got and the next command.

## 1. One input BAM per dataset

Every tool consumes the **same BAM** per dataset. Tools that cannot start from
a BAM (FASTQ-only) must (a) say so in their wrapper header and `status/<tool>.md`,
and (b) derive FASTQ **from that same BAM** (`samtools fastq` keeping CB/UB in
read comments, or `bamtofastq` from 10x) — never from a different upstream —
so every tool sees identical reads.

### Dataset 1: `pbmc_10k_v3`

| item | value |
|---|---|
| BAM | `/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam` (44,198,457,238 B when complete; download tracked in `status/data_pbmc10k.md`) |
| provenance | CellRanger 3.0.0, reference `refdata-cellranger-GRCh38-3.0.0` (Ensembl 93 annotation), STAR 2.5.1b |
| chrom naming | **bare Ensembl** — `1..22, X, Y, MT` + accessioned scaffolds (`GL000225.1`, `KI270728.1`, …), 194 `@SQ` lines. Verified from the BAM header 2026-08-18. (This CORRECTS the guess in `status/data_pbmc10k.md` that the BAM was chr-prefixed.) |
| tags | `CB`/`UB` (corrected barcode/UMI; ~98% of reads carry `CB`), plus `CR/UR`, `GX/GN`, `RE`, `xf`. Tolerate missing tags. |
| read length | R2 = 91 bp (verified) |
| cells | 11,769 filtered barcodes in `pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz` |
| GTF | `/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf` (bare Ensembl names — matches this BAM directly) |

## 2. Standard result directory (the contract)

Each wrapper writes **exactly** this layout:

```
results/benchmark_tools/<dataset>/<tool>/
├── pas.bed            REQUIRED  BED6 point-mode PAS calls (see below)
├── counts.mtx         if the tool emits a cell x PAS matrix (MatrixMarket;
│                      sidecars allowed: counts_barcodes.tsv, counts_features.tsv)
├── runtime_mem.txt    REQUIRED  `/usr/bin/time -v` of the main tool invocation
├── tool_version.txt   REQUIRED  version string (+ git commit if built from source)
├── log                REQUIRED  full stdout+stderr of the run
└── run/               tool-native raw output, kept as-is for auditing
```

### pas.bed rules

- BED6: `chrom  start  end  name  score  strand`; **tab**-separated, no header.
- **Point coordinates**: each row is the single inferred cleavage base
  (`end == start + 1`). Tools that emit peak *intervals* collapse each peak to
  its 3'-most base by strand — the same convention as
  `scripts/manuscript_figures/null_control.py`:

  ```sh
  awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){s=$3-1;e=$3}else{s=$2;e=$2+1}print $1,s,e,$4,$5,$6}'
  ```

- **Ensembl chromosome names** (`1..22, X, Y, MT`, accessioned scaffolds).
  Tools that print `chr1/chrM` names pipe through the translation helper:
  `./normalize_chroms.py --to ensembl` (handles `chrM <-> MT` and UCSC-style
  scaffold names; see `--help` / `--selftest`).
- `score` = the tool's confidence / read support if it has one, else `0`.
- `strand` is required; strandless callers are noted in `status/<tool>.md` and
  scored unstranded *in addition to* the standard strand-matched mode, never
  instead of it.
- Sorted with `LC_ALL=C sort -k1,1 -k2,2n`.

## 3. Scoring

`./evaluate_vs_reference.py` scores any `pas.bed` the same way:
point-mode, **strand-matched** `bedtools closest -s -d -t first`, cutoffs
**10 / 25 / 50 / 100 / 200 bp**, with an optional shuffled-null control
(`bedtools shuffle -chrom -noOverlapping`, genome-wide and/or gene-body
restricted). The logic is adapted from
`scripts/manuscript_figures/null_control.py` — read that file's header for why
point-mode + strand matching + a matched null is the only honest geometry when
the reference is dense (18.4 M PolyASite entries put ~most of the genome
within 100 bp of a site; interval-mode precision is inflated by peak width).

Reference PAS sets (all bare-Ensembl naming, all on read-only ssd2):

| ref | path | note |
|---|---|---|
| PolyASite 3.0 full | `/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed` | 18.4 M rows; density caveat above — never report precision against this alone |
| PolyASite 3.0 hi-conf | `/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_hiconf_score0.1.bed` | preferred headline reference |
| PolyASite 2.0 | `/mnt/ssd2/Laugney_Aligned/refs/polyasite_2.0_GRCh38_sorted.bed` | for comparability with older papers |

Shuffle support files (bare-Ensembl, writable copies):
genome = `data/references/chrom.sizes.nochr.filt` (24 chroms),
gene bodies = `data/references/gene_end.bed`.

Example (not run automatically):

```sh
./evaluate_vs_reference.py results/benchmark_tools/pbmc_10k_v3/<tool>/pas.bed \
    --reference /mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_hiconf_score0.1.bed \
    --genome    data/references/chrom.sizes.nochr.filt \
    --null genome,genic --gene-bodies data/references/gene_end.bed \
    --recall --out results/benchmark_tools/pbmc_10k_v3/<tool>/eval_hiconf.tsv
```

## 4. Adding a tool — wrapper checklist

1. Install into an isolated env (rules in section 0); write `status/<tool>.md`.
2. Copy the shape of `run_peakatail.sh` (the reference wrapper). A wrapper must:
   - `set -euo pipefail`, `export LC_ALL=C`;
   - take the dataset BAM as its single data input (FASTQ derivation noted if unavoidable);
   - wrap the main tool command in `/usr/bin/time -v -o .../runtime_mem.txt`;
   - tee tool stdout+stderr to `log`;
   - write `tool_version.txt`;
   - convert native output -> `pas.bed` per section 2 (point collapse +
     `normalize_chroms.py --to ensembl` + `LC_ALL=C sort`);
   - copy/convert a counts matrix to `counts.mtx` when the tool makes one;
   - leave the tool's native output under `run/`.
3. Same machine, same BAM, threads pinned to the same count (`THREADS`, default 16)
   for every tool on a given dataset — runtimes are only comparable that way.

## 5. Files here

| file | role |
|---|---|
| `README.md` | this protocol |
| `normalize_chroms.py` | bidirectional `chr`/Ensembl chromosome-name mapper (`chrM <-> MT` incl.) |
| `evaluate_vs_reference.py` | uniform scorer: point-mode, strand-matched, cutoffs 10/25/50/100/200, shuffled-null option |
| `run_peakatail.sh` | reference wrapper for PeakATail (`ema run`) on pbmc_10k_v3 — written, **not yet executed** (BAM still downloading) |
| `status/` | per-tool install/status notes + dataset download status |
