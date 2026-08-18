# Benchmark dataset 1: 10x pbmc_10k_v3 (download status)

Date started: 2026-08-18 23:01 +03
Target dir: `/mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3`
Download method: `nohup wget -c -P <dir> -o <dir>/<file>.log <URL> &` (resumable; re-running
the same command resumes a partial file).

## URLs

| File | URL | Size | Status |
|---|---|---|---|
| BAM | https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam | 44,198,457,238 B (41 GiB) | IN PROGRESS (PID 809842) |
| BAM index | https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam.bai | 17,316,008 B | COMPLETE (saved 17,316,008/17,316,008) |
| Filtered matrix | https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz | 94,334,700 B | COMPLETE (saved 94,334,700/94,334,700; `gzip -t` OK; contains filtered_feature_bc_matrix/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}) |
| Web summary | https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_web_summary.html | 4,606,473 B | COMPLETE (saved 4,606,473) |

PIDs at launch are recorded in `<dir>/download_pids.txt`
(BAM_PID=809842, BAI_PID=809843, MATRIX_PID=809844, WEBSUMMARY_PID=809845; started 23:01:12).

## Checksums / md5

10x Genomics publishes NO md5 checksums for the 3.0.0 pbmc_10k_v3 sample files — there is no
`md5sum` listing on the dataset page or on cf.10xgenomics.com. Integrity checks available instead:

1. Exact byte size: the server reports `Length: 44198457238` for the BAM (confirmed in the wget log).
   After completion verify `stat -c%s pbmc_10k_v3_possorted_genome_bam.bam` == 44198457238.
2. `samtools quickcheck -v pbmc_10k_v3_possorted_genome_bam.bam` (verifies BGZF EOF block + header).
3. wget itself only reports "saved [N/N]" when the byte count matches Content-Length.

## How to check progress

```
tail -2 /mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam.log
ls -l  /mnt/ssd1/Projects/PeakATail_wd/data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
```
The log ends with `... saved [44198457238/44198457238]` when done. Observed rate ~14 MB/s at start
(1.3 GB in the first 90 s) => ETA roughly 1 h if sustained. If the process dies, resume with the exact
same `nohup wget -c ...` command (see `download_pids.txt` dir); `-c` continues from the partial file.

## Expected BAM contents (for the benchmark harness)

- CellRanger 3.0.0 position-sorted BAM, GRCh38 (10x reference **GRCh38-3.0.0**, Ensembl 93 annotation).
- Per-read tags: `CB` (corrected cell barcode, e.g. `AAACCCA...-1`), `UB` (corrected UMI),
  plus raw `CR/UR`, quals `CY/UY`, gene tags `GX/GN`, region tag `RE`, and `xf`.
  Reads lacking a corrected barcode have no CB/UB tag — tools must tolerate missing tags.
- Chromosome names are **'chr'-prefixed** (`chr1 ... chr22, chrX, chrY, chrM`).
  IMPORTANT NAMING MISMATCH: PolyASite 2.0 (and most Ensembl-derived PAS atlases) use bare Ensembl
  names (`1 ... 22, X, Y, MT`). The benchmark harness MUST translate one side before intersecting
  (e.g. `sed 's/^/chr/; s/^chrMT/chrM/'` on the atlas BED, or strip `chr` from tool outputs).
  Also remember `chrM` (UCSC/10x) vs `MT` (Ensembl).
- Cell list: 11,769 cells in the filtered matrix (`filtered_feature_bc_matrix/barcodes.tsv.gz`);
  use these barcodes to restrict PAS calling to real cells.
- Locale note for this machine: `export LC_ALL=C` before any `sort` used by the harness
  (tr_TR locale breaks lexicographic chromosome sorting).

## Next step after BAM completes

```
export LC_ALL=C
samtools quickcheck -v pbmc_10k_v3_possorted_genome_bam.bam && echo BAM_OK
samtools view pbmc_10k_v3_possorted_genome_bam.bam | head -1   # eyeball CB/UB tags
```

## Dataset 2 (not downloaded yet)

Fetch script written (NOT run): `/mnt/ssd1/Projects/PeakATail_wd/scripts/benchmark_tools/fetch_gse104556.sh`
GSE104556 = SRP119327 = SRR6129050 (Mouse1_scRNAseq) + SRR6129051 (Mouse2_scRNAseq),
mouse testis, 10x Chromium 3' v2, ~199M read pairs each. Needs CellRanger + mm10 reference to
produce a CB/UB-tagged BAM (see script header).

## CORRECTION 2026-08-18 (verified from the partial BAM, harness agent)

The "Expected BAM contents" section above **guessed wrong on chrom naming**:
the header of the actual BAM shows **bare Ensembl names** (`SN:1 ... SN:22, X,
Y, MT` plus accessioned scaffolds `GL000225.1 / KI270728.1 ...`; 194 `@SQ`
lines; `@PG` confirms `refdata-cellranger-GRCh38-3.0.0`, STAR 2.5.1b). There
is NO `chr` prefix and mitochondrion is `MT`, so the BAM matches the Ensembl
GTF and the PolyASite references directly. Also verified: R2 length = 91 bp;
`CB` tag present on ~98% of the first 2000 reads. The translation helper
(`normalize_chroms.py`) is still needed for competitor tools that emit UCSC
names, just not for this BAM.
