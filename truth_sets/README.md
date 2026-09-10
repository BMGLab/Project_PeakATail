# truth_sets — the long-read reference sets the accuracy numbers are scored against

## Kinnex (PacBio, human PBMC)

Strand-aware 1-bp **point** BEDs, gzipped. These are the files `score_tool.py` consumes;
they are derived products, not raw data. The script that builds them from the vendor
downloads is in `scripts/benchmark_tools/` — see `DATA_ACCESSIONS.md` for the download
location of the source Kinnex libraries.

Two libraries:

- `x3p/`  — Kinnex 10x 3' library
- `gemx/` — Kinnex GEM-X library

Per library:

| file | what |
|---|---|
| `*_truth_t5.point.scorable.bed.gz` | truth sites at ≥5 supporting UMI — the primary threshold reported in the paper |
| `*_truth_t20.point.scorable.bed.gz` | ≥20 UMI, the stringent arm |
| `*_truth_t100.point.scorable.bed.gz`, `*_truth_t500.point.scorable.bed.gz` | higher thresholds, used for the stringency sweep |
| `*_decoy*.point.bed.gz` | the decoy set — long-read-supported positions that are *not* poly(A) termini, used for the IP-decoy read-out |

"scorable" means the set has already been restricted to the contigs the scorer accepts, so
that truth and call set share a denominator.

**Caveat carried from the manuscript.** The Kinnex donor is not the donor of the 10x PBMC
libraries, so these sets are an atlas-style, site-level truth, not cell-matched truth; and
the support counts are records, not de-duplicated molecules. Both limitations are stated in
the paper where the numbers are used.

The full uncompressed set (48 files, 576 MB, including intermediate and combined
truth+decoy files) is not distributed here; it is regenerable from the build script.
