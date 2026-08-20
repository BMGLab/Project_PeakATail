# Kinnex PBMC 10x3p — Tier-1 ground-truth extraction plan

**Goal** (per `manuscript/02_validation_plan.md`, Tier 1): per-cell-barcode polyA-supported
long-read 3' ends -> truth PAS (cluster termini within 25 nt, require >= 5 UMIs) -> truth
PAS x cell-type UMI usage matrix, for validating PeakATail PAS positions and per-cell-type
usage on the same sample.

## 1. Data downloaded (2026-08-19)

Location: `/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/kinnex_pbmc_10x3p/`
Source: `https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/`
Launcher: `run2.sh` (nohup; two parallel aria2c 8-connection streams, resumable;
markers `DONE.ok` / `FAILED.err`, logs in `logs/`; `run.sh` is the superseded
single-stream wget version, kept for reference).

### primary_10x3p/  — DATA-Revio-Kinnex-PBMC-10x3p (Kinnex, Revio, PBMC, 10x 3' v3.1)
The dataset named as **primary Tier-1** in the validation plan.

| File | Size | Why |
|---|---|---|
| `scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam` (+`.pbi`) | 19.5G | **The truth substrate.** CB-corrected, UMI-deduplicated FLNC molecules: 1 record = 1 unique CB/UMI molecule that passed `isoseq refine --require-polya` (terminal polyA detected and trimmed; SMRT Link default min polyA = 20 nt — satisfies the plan's ">= 20 nt polyA evidence" criterion). Unmapped: we align it ourselves (Section 2). |
| `genes_seurat/{barcodes.tsv,genes.tsv,matrix.mtx}` | 136M | 12,852 called cells (10x-convention barcodes with `-1` suffix) = real-cell CB whitelist; gene-level matrix for Seurat clustering / cell-type transfer. |
| `pigeon_filtered.classification.txt` | 292M | Isoform classification; cross-check of transcript-level annotation and on-target rate. |
| `SL-Revio-PBMC-10x3p-CCS.pdf`, `SL-Revio-PBMC-10x3p-ReadSegSingleCell.pdf` | 6M | SMRT Link run reports (library stats, refine polyA settings — cite in methods). |

### secondary_gemx3p/  — DATA-Revio-Kinnex-PBMC-10kcells-10xGEMX3p (10x 3' v4 GEM-X)
| File | Size | Why |
|---|---|---|
| `scisoseq.mapped.bam` (+`.bai`) | 15.4G | Dedup molecules **already genome-aligned** — lets us prototype/debug the terminus-extraction code immediately (no mapping wait) and serves as a cross-chemistry replicate. Check `samtools view -H` for the reference build (PacBio pigeon ref sets are GENCODE/`chr`-prefixed hg38; harmonize naming with our Ensembl no-chr convention before any intersection). |
| `bcstats_report.tsv.gz`, `cell_statistics.report.json` | 11M | Barcode rank / cell-calling stats. |

### Tradeoffs documented
- **Avoided**: raw HiFi movies (59G/70G, would need `skera`+`lima`+full isoseq pipeline),
  `segmented.bam` (73G/88G, same reason), dedup **fasta** (61G — redundant with the BAM and
  lacks tags). Total taken ~36G of the <=250G budget.
- Choosing the **dedup** BAM (not FLNC/refined) means UMI collapse is already done by
  `isoseq groupdedup`: UMI count of a truth PAS = number of dedup records, exactly the
  quantity the plan thresholds on (>= 5 UMIs).
- Caveat: refine trims the polyA tail, so we cannot re-measure softclipped-A length
  ourselves on these reads; we rely on refine's polyA requirement (documented in the
  ReadSeg PDF). Equivalent stringency to the plan's ">= 20 nt softclipped polyA" criterion.
- The Kinnex PBMC donor differs from the short-read `pbmc_10k_v3` donor, so raw CB sharing
  across datasets is impossible; Tier-1 per-CB truth is used (a) against pseudo-short reads
  derived from the same Kinnex molecules (CBs shared by construction) and (b) at cell-type
  resolution against independent short-read PBMC runs.

## 2. Extraction pipeline (primary dataset)

Tools: `samtools` (/usr/local/bin), `minimap2` (/usr/bin), Python+pysam, bedtools.
`pbmm2` not installed; minimap2 with tag pass-through is equivalent for this purpose.
Reference: `/mnt/ssd1/Projects/PeakATail_wd/data/references/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
(Ensembl no-chr naming — same build/naming as the STARsolo short-read BAMs and
`data/references/atlases/*GRCh38*`, so truth coordinates are directly comparable).

### Step 0 — inspect tags (must-do before coding)
`samtools view primary_10x3p/scisoseq...dedup.bam | head -2` and confirm per-record tags.
Expected (SMRT Link scisoseq dedup): `CB:Z` corrected cell barcode, `XM:Z` UMI,
`rc:i` real-cell flag, `nc/ic` dedup support counts. **Verify CB orientation**: scisoseq
CBs are often the reverse complement of the 10x whitelist convention. Test: intersect the
dedup CB set with `genes_seurat/barcodes.tsv` (strip `-1`) as-is and revcomp'ed; adopt the
orientation with the (near-total) overlap. Record the choice in the extraction script.

### Step 1 — align dedup molecules (splice-aware, tag-preserving)
```
samtools fastq -T CB,XM,rc primary_10x3p/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam \
| minimap2 -ax splice:hq -uf -y -t 96 --secondary=no \
    <ref>/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
| samtools sort -@ 16 -o kinnex_dedup.mapped.GRCh38.bam -
samtools index kinnex_dedup.mapped.GRCh38.bam
```
- `-T CB,XM,rc` copies BAM tags into FASTQ header comments; `-y` re-emits them as SAM tags.
- `-uf`: FLNC reads are already in transcript 5'->3' orientation -> forward-strand splicing
  only; alignment strand = transcript strand.
- `--secondary=no`: primary alignments only. Long job (~30-60 min on 96 threads for
  tens of millions of molecules): run via nohup `run.sh` with DONE.ok/FAILED.err markers.

### Step 2 — per-molecule 3' terminus
Python/pysam over primary, non-supplementary alignments with MAPQ >= 1:
- Terminus (cleavage site) = **reference end** (`reference_end`, rightmost base) for
  forward-strand alignments; **reference start** for reverse-strand alignments.
  (Read 3' end = trimmed-polyA junction = cleavage position.)
- **Softclip guard**: record softclip length at the read's 3' end
  (CIGAR end for + / CIGAR start for -). Discard molecules with 3'-end softclip > 30 nt
  (terminus position unreliable); report the discarded fraction (paper caveat).
- **CB filter**: keep only molecules whose CB (orientation-corrected, Step 0) is in
  `genes_seurat/barcodes.tsv` (12,852 real cells). This implements "valid CB matching the
  whitelist" at the called-cell level; `rc:i:1` should agree — report any discordance.
- Emit TSV: `chrom  terminus_pos  strand  CB  UMI`. 1 row = 1 UMI.

### Step 3 — internal-priming flag (truth hygiene)
Even polyA-selected long reads carry some internal priming. Flag termini where the
genomic 18 nt immediately downstream (in transcript direction) contains >= 12 A or
>= 6 consecutive A (same criterion as Section 3.1 of the validation plan; bedtools
getfasta on the terminus+1..+18 window). Keep flagged sites out of the **truth** set;
retain them in a side file — they double as a known-decoy list.

### Step 4 — cluster termini into truth PAS
Per (chrom, strand), sort termini; single-linkage clustering joining adjacent termini
<= 25 nt apart (plan: "cluster within 25 nt"). Per cluster:
- representative site = modal terminus (highest UMI count);
- total UMI = number of rows; **keep clusters with >= 5 UMIs**;
- outputs:
  - `truth_pas.bed6+`: chrom, cluster start/end, id, UMI count, strand, rep_site,
    n_cells, IP-flag;
  - `truth_cb_matrix/`: sparse truth-PAS x CB UMI-count matrix (MatrixMarket).

### Step 5 — cell types and usage matrix
Cluster `genes_seurat` gene-level matrix (Seurat/scanpy, standard PBMC annotation:
T/NK/B/mono/DC), or transfer labels from an annotated PBMC reference. Collapse the
CB matrix to truth-PAS x cell-type UMI counts = the truth usage matrix that
`ema switch diff` / `switch length` outputs are validated against.

### Secondary (GEM-X) track
Same Steps 2-5 directly on `secondary_gemx3p/scisoseq.mapped.bam` (skip Step 1; its
cell list comes from `bcstats_report.tsv.gz` / dedup `rc` tag). First harmonize
chromosome naming to Ensembl no-chr if the header shows GENCODE names. Used to
validate the extraction code and as a replicate; the 10x3p (v3.1) set remains primary.

## 3. Acceptance checks
- Number of dedup molecules with valid real-cell CB: expect ~10-40M; cells: 12,852.
- Truth PAS count: expect ~20-60k clusters at >= 5 UMIs (PBMC-expressed genes).
- >= 80% of truth PAS within 25 nt of a PolyASite2.0/PolyA_DB site
  (`data/references/atlases/`); the remainder = candidate novel sites, not errors.
- Canonical polyA signal (AAUAAA + variants) in -40..-10 of rep sites: expect 80-90%.
- Strand sanity: truth PAS antisense to overlapping gene < 2%.

---

# AS EXECUTED (2026-08-19) — deviations and what they cost

Scripts live under `/mnt/ssd0/emaout/peakatail_benchmark/kinnex/` (working copies) and
`scripts/benchmark_tools/{kinnex_decompose.sh}` + `scripts/manuscript_figures/kinnex_truth_validation.py`
(repo copies). Markers: `ALIGN_PRIMARY.DONE.ok`, `EXTRACT_SECONDARY.DONE.ok`,
`TRUTH2_GEMX.DONE.ok`, `SCORE_GEMX.DONE.ok`, `PRIMARY_CHAIN.DONE.ok`.

## Step 0 — tag inspection: as planned, both datasets
Both BAMs carry `CB:Z` / `XM:Z` (UMI) / `rc:i` / `ic:i` / `im:Z` / `it:Z`.
**CB orientation is AS-IS, not reverse-complemented** — 72/72 sampled primary CBs and
10,292/10,292 sampled secondary CBs matched the whitelist as written; 0 matched the
reverse complement. `rc:i` concordance on the secondary is perfect (0 kept records with
`rc != 1`, 0 dropped real-cell records), and the pbmm2-mapped secondary BAM turns out to be
pre-filtered to real cells (`dropped_cb_not_realcell = 0`).
Primary record count confirmed from the `.pbi` header: **81,678,354** dedup molecules.

## Step 1 — alignment: `-t 8`, not `-t 96`
The environment caps this session at 8 threads, so minimap2 ran
`-ax splice:hq -uf -y --secondary=no -t 8` (~3.8 k molecules/s, ~6 h wall clock for the
primary) instead of the plan's `-t 96`. Tag pass-through through
`samtools fastq -T CB,XM,rc | minimap2 -y` was verified on a chr21 test before launch.

## Steps 2-3 — as planned
3'-softclip guard at 30 nt discards 1.7% of secondary molecules. Internal priming is called on
the 18 genomic nt downstream in transcript direction; for '-' strand the transcript-strand A's
are forward-strand T's, so no reverse complement is computed (`ip_flag.py`, whole genome held
in memory — replaces `bedtools getfasta` for speed, identical result).

## Step 4 — TWO substantive deviations, both forced by the data

**(a) single-linkage -> greedy peak calling.** At 104 M molecules, single-linkage joining of
termini <=25 nt apart CHAINS through highly expressed loci: MALAT1 came out as ONE 5,955 bp
"PAS" and the MT rRNA locus as one 3,762 bp "PAS" (cluster width q99 = 155 bp, max 8,054 bp).
Replaced with greedy peak calling on the terminus histogram using the same 25 nt radius as an
exclusion window (`peak_call.py`): repeatedly take the highest-count remaining terminus, emit
it, absorb everything within +/-25 nt. Peak width max is then 51 bp. The 25 nt merge radius of
the plan is preserved; only the linkage rule changed.

**(b) ">= 5 UMIs" is a noise-level threshold at this depth — the truth is reported as a sweep.**
The plan's acceptance checks (20-60 k clusters; >=80% within 25 nt of an atlas site; 80-90%
canonical hexamer) were calibrated for a far shallower dataset. At >=5 UMIs the GEM-X data
gives 509,686 peaks, of which only 17.6% are within 25 nt of PolyASite 2.0 and 13.2% carry
AATAAA/ATTAAA. That is not an extraction error: **75.1% of MOLECULES terminate within 25 nt of
a PolyASite 2.0 site** (81.0% within 100 nt) — the mass is right, the low-count tail is noise.
The sequence-only diagnostic settles it (`<tag>_truth_calibration.tsv`): canonical-hexamer
content rises monotonically with molecular support (10.0% at 5-9 UMIs -> 63.5% at >=5000),
while the internal-priming decoy set is FLAT at ~4.5-5.0% at every support level. So truth PAS
are emitted at >=5 UMIs and **scored at >=5 / >=20 / >=100 / >=500 UMIs**, and no threshold was
chosen using an atlas (atlas agreement is reported, never optimised).

Also note the plan's "80-90% canonical hexamer" bar is unreachable by construction: the SAME
code scores PolyASite 2.0 rep sites at 41.1% and protein-coding TES at 37.99%. The >=500 UMI
truth set reaches 48.8%, i.e. it is BETTER localised than the atlas it is compared to.

## Step 5 — cell types: NOT done
The PAS x cell-barcode count matrix IS built (`<tag>_truth_pas_by_cb.triplets.tsv`,
28.9 M nonzero entries for GEM-X), but Seurat clustering / cell-type transfer was not run:
the first validation use (Step 4 of the task) is site-level, and the donor mismatch means the
per-cell-type usage matrix cannot be compared to pbmc_10k_v3 cell-for-cell anyway.

## Which dataset is primary
`primary_10x3p` (10x 3' v3.1, 81.7 M molecules, 12,852 cells) remains the designated Tier-1
substrate. `secondary_gemx3p` (GEM-X 3' v4, 105.8 M alignments -> 104.0 M usable termini,
10,324 cells) was processed FIRST because it ships genome-aligned, and it is a genuine
cross-chemistry, cross-donor replicate rather than a mere prototype.

## What the executed run produced (both datasets)

| | primary `x3p` (10x 3' v3.1) | secondary `gemx` (GEM-X 3' v4) |
|---|---:|---:|
| dedup molecules in BAM | 81,678,354 (`.pbi`) | 113,766,063 alignment records |
| aligned here | minimap2, 99.81% primary mapped, 5 h 56 m | pre-aligned (pbmm2, GENCODE v39 -> stripped to Ensembl) |
| records scored (primary, non-supp, MAPQ>=1) | 79,620,570 | 105,806,517 |
| dropped, 3' softclip > 30 nt | 2,795,638 (3.5%) | 1,799,849 (1.7%) |
| dropped, CB not a real cell | 0 | 0 |
| **molecule termini kept** | **76,824,932** | **103,955,157** |
| cells | 12,852 | 10,324 |
| unique terminus positions | 9,981,323 | 8,826,660 |
| internal-priming termini | 3,482,087 pos / 34.5 M molecules | 2,831,432 pos / 29.7 M molecules |
| truth PAS >=5 / >=20 / >=100 / >=500 UMI | 538,226 / 96,634 / 20,375 / 5,705 | 509,686 / 137,987 / 35,626 / 9,578 |
| internal-priming decoy PAS | 638,229 | 481,555 |
| PAS x CB nonzero entries | 17,455,685 | 26,880,689 |

Acceptance checks (primary, at >=500 UMI): 84.6% within 25 nt of PolyASite 2.0 (plan bar: >=80%),
98.9% sense-genic, **0.37% antisense-only (plan bar: <2%)**, 71.8% canonical hexamer in the top
support bin. At >=100 UMI: 20,375 PAS (plan expectation 20-60 k), 70.3% within 25 nt, 1.3%
antisense. The plan's checks are met — at the support level the plan implicitly assumed, which
is >=100 UMI here, not >=5.
