# 26 — Second human donor: pre-registration of the generalisation test

**Written 2026-08-21, BEFORE any accuracy number from a second human donor exists.**
At the time of writing the only things that exist are (a) the HTTP headers and BAM *header* of the
candidate library, read to establish that it is technically usable at all, and (b) this file. The
caller has not been run on it. No PAS BED, no score TSV, no precision, no recall.

`manuscript/24_prime_preregistration.md` is owned by the concurrently running `peakAtail-prime`
job (last written 2026-08-21 23:34, process alive at the time this file was created), so this
pre-registration is a separate file rather than an appended section, as the PI instruction allowed.

---

## 1. Why this test exists

The manuscript's headline reliability claim — the pre-registered precision default reaches
**P@100 0.7062** on human PBMC ([19](19_final_gate_v2.md)) — rests on **one human donor and one
CellRanger BAM** (`pbmc_10k_v3`). Two mice ([19](19_final_gate_v2.md) §, testis mouse 1 and mouse 2)
provide biological replication in *mouse*, but nothing replicates the human arm. A reviewer is
entitled to ask whether 0.7062 is a property of the caller or a property of that one library.

This is therefore a **generalisation test**, not a tuning exercise. Its whole evidential value comes
from the fact that the answer is fixed here, before the number is seen.

---

## 2. The dataset

**10x Genomics `pbmc4k` — "4k PBMCs from a Healthy Donor".**

| property | value | how established |
|---|---|---|
| BAM URL | `https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_possorted_genome_bam.bam` | HTTP 200, `content-length: 33243581299` (33.24 GB), `last-modified: 2017-11-08` |
| chemistry | Single Cell 3' **v2** | `pbmc4k_web_summary.html`, Chemistry row |
| CellRanger | 2.1.0 | `pbmc4k_web_summary.html`, Cell Ranger Version row |
| reference | `GRCh38-1.2.0` | BAM `@PG`/`@CO` STAR command line (`--genomeDir /mnt/opt/refdata_cellranger/GRCh38-1.2.0/star`) |
| contig naming | Ensembl (`1`, `2`, … no `chr`), 194 `@SQ`, `SN:1 LN:248956422` | BAM header |
| estimated cells | 4,340 | web summary hero metric |
| mean reads/cell | 87,433 | web summary hero metric |
| read length (R2) | **98 bp** (200,000/200,000 first records) | `samtools view … \| awk '{print length($10)}'` on the partial download |
| cell-barcode length | 16 bp (`CB:Z:` + `-1` suffix) | first 50,000 records |
| UMI length | 10 bp (`UB:Z:`) | first 50,000 records |
| CB / UB tag presence | 197,083 / 199,972 of the first 200,000 records | as above |

**Why `pbmc4k` and not `pbmc8k`.** [04](04_datasets.md) lists both. `pbmc8k`'s
`possorted_genome_bam.bam` is **not hosted** — every URL variant probed this session returns HTTP
400/403 from the 10x CloudFront distribution, while its `filtered_gene_bc_matrices.tar.gz` returns
200. Only the count matrix is distributed, and a count matrix cannot feed a 3'-end peak caller.
`pbmc4k` is therefore the only second public human 10x 3' library with a CB/UB-tagged BAM available
without realignment. **This is recorded before the run, not chosen after seeing a result.**

**What this donor is and is not.** 10x publishes no donor identifier for `pbmc4k`, so it cannot be
*proved* to be a different individual from the `pbmc_10k_v3` donor. What is documented and verifiable
is that it is a **different library, prepared with a different chemistry (v2 vs v3), sequenced on a
different flow cell (`H53GNBCXY`), and processed with a different CellRanger version (2.1.0 vs
3.0.0)**, released a year apart. The paper must claim exactly that and no more: this is a
second-library / second-chemistry generalisation test, described as such. If a reviewer requires a
provably distinct donor, the honest answer is that no second public human 10x 3' library ships a
BAM with a donor identifier.

---

## 3. What is held fixed (the whole point)

The caller is the frozen merged tree `tools/pa-polya-run-9dfdefb3`
(`9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53`) — the same tree that produced every number in
[19](19_final_gate_v2.md). The run asserts the imported `ema` module path before starting, exactly as
`scripts/benchmark_tools/stage2_final_launch.sh` does.

The flag set is copied verbatim from
`results/benchmark_tools/pbmc_10k_v3/peakatail_clipseeded_final_v2_ipfilt/runtime_mem.txt`:

```
run --bam-dir <BAM> --gtf <HUM_GTF> --output <OUT> --threads <T>
    --barcode-tag CB --cb-len 16 --seq-len <L> --ignore-chro MT
    --plot-engine none --no-progress --peak-strategy clip_seeded
    --ip-filter --genome-fasta <HUM_FA> --ip-filter-mode filter
```

**Three arguments differ from the PBMC 10k v3 arm, and only three. Each is declared here, with its
justification, before the run:**

1. `--bam-dir`, `--output` — input and output paths. Not behaviour.
2. `--seq-len 98` instead of `91`. **This is a library descriptor, not a tuning knob.** In
   `ema/countmatrix/read.py` (frozen tree, lines ~100-110) `seq_len` normalises the read footprint:
   a record whose reference span *exceeds* `seq_len` is **discarded**, and a shorter span is
   *extended* to `seq_len`. Running a 98 bp library at `--seq-len 91` would throw away essentially
   every unspliced read; the parameter has one correct value per library and it is the read length.
   `pbmc_10k_v3` reads are 91 bp (100,000/100,000 first records, checked this session); `pbmc4k`
   reads are 98 bp. The precedent is already in the manuscript: the two testis mouse arms — v2
   chemistry, 98 bp — were run at `--seq-len 98` in the same `stage2_final_launch.sh` that ran PBMC
   at 91. **No other numeric parameter moves.** In particular `--default-threshold 5`,
   `--merge-len 100`, `--min-pas-prominence 5.0`, `--min-pas-spacing -1`, `--pas-gap 100`,
   `--max-pas 5`, `--smoothing-window 50`, `--floor-threshold 3`, `--ip-a-stretch 6`,
   `--ip-a-fraction 0.7` and the IP window all keep their PBMC values, which are the tool defaults.
   (`--min-pas-spacing -1` means the PAS merger's distance tier is inferred from the BAM's own
   median read length by `infer_median_read_length`; that is automatic and identical in kind to what
   the PBMC run did.)
3. `--threads` will be **8**, not 16, because other jobs own this box and the PI capped this task at
   12 threads. This is a resource parameter with **no effect on output**: `ema/countmatrix/
   chrom_parallel.py` partitions one job per (contig, strand) regardless of worker count and merges
   in the order the sequential run would have emitted (`+` strand first, contigs in BAM header
   order, then `-`), which its own docstring records as byte-identical to the sequential two-pass
   run, verified on the chr19+21 PBMC slice and pinned by `tests/test_chrom_parallel_identity.py`.
   Wall time and peak RSS will therefore not be comparable to the PBMC arm's 32:59 / 11.1 GB and
   will be reported as a run record only, never as a benchmark.

A `.bai` index must be built locally with `samtools index` because 10x does not distribute one for
this dataset (HTTP 403 on every `.bai` variant probed). The index is derived from the downloaded BAM
and changes no read.

**NOTHING WILL BE RETUNED.** No threshold, no filter, no tier rule, no post-hoc arm selection will be
adjusted after seeing a `pbmc4k` number. If the result is bad, the result is bad and it is reported
as bad. Any parameter change motivated by a `pbmc4k` score would void this test, and this paragraph
is the record that says so.

---

## 4. Scoring protocol

Identical scorer, identical references, identical nulls: `scripts/benchmark_tools/score_tool.py`,
invoked exactly as the human arms are invoked in `stage2_final_launch.sh`:

```
score_tool.py <pas.bed> <label> --outdir <out> --detected-atlas <donor2 pas2.in_detected_genes.bed>
```

Everything else stays at the scorer's human GRCh38 defaults, which is what makes the comparison
legitimate:

- **precision** vs `data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6` (569,005 rep sites)
  and `tes.protein_coding.GRCh38.99.bed6`, strand-matched `bedtools closest -s -d -t first`, at
  10/25/50/100/200 bp. **P@W = the fraction of calls with a strand-matched PolyASite 2.0
  representative site within W bp** — the standing definition, unchanged.
- **nulls** — 3 seeds of width-preserving `bedtools shuffle -chrom -noOverlapping -incl
  results/benchmark_tools/shared_refs/genebodies.merged.bed`, the same inclusion file and the same
  seeds 1/2/3 as every other arm.
- **recall (full)** vs the whole 569,005-site atlas — identical denominator to every other arm.
- **R_det** vs a detected-gene-restricted atlas **built for THIS library** by the recipe already
  recorded in `results/benchmark_tools/shared_refs_pbmc/README.md` and
  `data/references/atlases/README.md` §6, reused rather than re-invented:
  `cut -f5 <run/annotatedpas.bed> | sort -u > detected_genes.txt`, restrict
  `data/references/gene_end.bed` on column 4 to those gene ids, then
  `bedtools intersect -a polyasite2.GRCh38.96.rep_sites.bed6 -b <detected gene bodies> -s -u`.
  All sorts under `LC_ALL=C`.
- **F1** = harmonic mean of full-atlas precision and the stated recall flavour, as `score_tool.py`
  already computes it.

**One deviation from the PBMC denominator, declared now.** The PBMC `shared_refs_pbmc` denominator
was built from the *native* `lambda_gradient` run's `annotatedpas.bed`
(`results/benchmark_tools/pbmc_10k_v3/peakatail/run/`), because that run existed first. Only one
`pbmc4k` run is being made, the `clip_seeded --ip-filter` one, so its `annotatedpas.bed` is the only
source available. To keep this from being a silent confound, the size of the strategy effect will be
**measured on PBMC, where both runs exist**: the detected-gene set and restricted-atlas size will be
rebuilt from the PBMC `clip_seeded_final_v2_ipfilt` run and compared with the shipped 14,949-gene /
285,136-site denominator. That overlap will be reported alongside the `pbmc4k` numbers whatever it
says.

### Arms to be scored

Both, and both reported, no matter which looks better:

1. **the pre-registered precision default** — tier-1 ∩ IP-pass ∩ ≥2 clip molecules
   (`pas_PRESPEC_precision_default.bed`, produced by the same `awk -F'\t' '$5>=2'` rule as the
   PBMC arm).
2. **the ≥1-molecule arm** — `pas_tier1.bed`. Reported because
   [22](22_performance_roadmap.md)/[23](23_algorithm_roadmap.md) established that ≥2 molecules is a
   *reliability* choice and not the F1 optimum, and the paper must not let a reader infer otherwise.

`pas.bed`, `pas_tier2.bed` and `pas_tier1_ge2umi_POSTHOC.bed` are emitted and scored by the same
loop; they are secondary and labelled as such.

---

## 5. Cross-donor site concordance (declared now, computed after)

The human analogue of the cross-mouse reproducibility panel, which the paper currently cannot show.
Defined here so the window cannot be chosen to flatter the answer:

> Of this donor's **pre-registered-default** calls, the fraction with a **strand-matched** call in
> the `pbmc_10k_v3` pre-registered-default set within **25 bp** and within **100 bp**
> (`bedtools closest -s -d -t first`, points already reduced to the 3'-most base by strand), **and
> the reciprocal direction** (`pbmc_10k_v3` → `pbmc4k`). Both directions reported; neither is
> allowed to stand alone.

Because the two sets have different sizes, the two directions will differ and both numbers are
reported as-is. A **genic-shuffle null** for the concordance is computed the same way
`score_tool.py` does its nulls (3 seeds, `-incl genebodies.merged.bed`), so the concordance is read
against chance rather than in the absolute.

---

## 6. Success criterion — fixed now

> **The pre-registered gate `P@100 >= 0.50` must hold on `pbmc4k` for the pre-registered precision
> default (tier-1 ∩ IP-pass ∩ ≥2 clip molecules).**

This is the *same* gate, the *same* number, that `stage2_final_launch.sh` pre-registered for
`pbmc_10k_v3` and both testis mice, and that [19](19_final_gate_v2.md) records as passed on all
three (0.7062 / 0.7450 / 0.7572).

Reported but **not gated**: R_det, F1, P@10, P@25, the nulls, the ≥1-molecule arm, and the
cross-donor concordance. They are descriptive. Fixing a threshold on them now, with no prior, would
be inventing a criterion rather than testing one.

### What a failure would mean for the paper

Stated before the number exists, so it cannot be softened afterwards:

- **P@100 ≥ 0.50 on `pbmc4k`.** The reliability claim generalises beyond a single human library and
  the paper may say so, with the honest scope: *two human 10x 3' libraries (v2 and v3 chemistry) and
  two mice*. The single-donor weakness is closed to the extent one extra public library can close it.
- **P@100 < 0.50.** The gate has **failed on a second donor** and the paper cannot claim a
  generalised precision-first default. The correct response is *not* to retune: it is to report the
  failure, restate the headline claim as **library-specific to `pbmc_10k_v3`**, and treat the
  chemistry difference (v2, 98 bp, 10 bp UMI, CellRanger 2.1.0, shallower per-cell depth) as a
  candidate explanation to be *tested* in a separate, separately pre-registered analysis — not
  asserted as an excuse. A failure here is a genuine, publishable limitation of the method, and
  burying it would be worse for the paper than reporting it.
- **P@100 between 0.50 and materially below 0.7062.** The gate passes and is reported as passing,
  but the paper must additionally state that absolute precision is chemistry- and depth-dependent
  and quote both values side by side rather than the PBMC value alone.
- **Run fails / library proves technically unusable.** Reported as BLOCKED with the specific reason.
  No substitution of a different assay type or a non-10x library.

---

## 7. Provenance

- Frozen code: `tools/pa-polya-run-9dfdefb3` @ `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53`
- Run outputs: `/mnt/ssd0/emaout/peakatail_benchmark/donor2/`
- Reference GTF / FASTA: `/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf`,
  `Homo_sapiens.GRCh38.dna.primary_assembly.fa` — the same two files the PBMC arm used.
- Everything in this file is PROVISIONAL until a verifier passes it.

---

## A1 — pre-run data checks (2026-08-21 23:53, still BEFORE any caller run)

Recorded here because they are properties of the *input data*, established while the BAM was still
downloading and before the caller had been started. No accuracy number existed when these were
written. They are stated now so they cannot later be presented as post-hoc justification.

1. **The two human BAMs occupy the same coordinate space.** `diff` of the sorted `@SQ` lines of
   `pbmc_10k_v3_possorted_genome_bam.bam` and `pbmc4k_possorted_genome_bam.bam` is **empty**: the
   same 194 contigs with the same lengths. The CellRanger reference *packages* differ
   (`refdata-cellranger-GRCh38-3.0.0` vs `GRCh38-1.2.0`, from the two `@PG`/`@CO` STAR command
   lines) but the assembly and contig set do not, so scoring both against
   `polyasite2.GRCh38.96.rep_sites.bed6` and annotating both with `Homo_sapiens.GRCh38.99.gtf` is
   the *same* operation, not an approximation. Note that the PBMC arm already has this same mild
   annotation-vintage mismatch (BAM built on the 3.0.0 refdata, scored with Ensembl 99); the second
   donor is treated identically rather than specially.
2. **Poly(A) clip rate.** The frozen tree's own QC function
   (`ema.countmatrix.polya.check_clip_rate`, first 200,000 CB-tagged reads, `min_clip=6`,
   `min_purity=0.8`) measures **2.2565 %** on `pbmc_10k_v3`. The same measurement on `pbmc4k` will
   be taken once the BAM is complete (pysam refuses a truncated BGZF stream) and reported next to
   the scores as a **descriptive covariate**, because the clip-seeded strategy's evidence channel is
   exactly this signal. It is *not* part of the gate and no parameter may be changed on the strength
   of it. (The function's docstring quotes a "PBMC v3 reference: 1.15 %"; the value measured this
   session with this code on this BAM is 2.2565 %, so the measured value is the one quoted here and
   the docstring figure is treated as stale rather than authoritative.)
3. **No `.bai` is distributed** for `pbmc4k` (HTTP 403 on every index-name variant probed), so one is
   built locally with `samtools index`. `ema/countmatrix/chrom_parallel.py` requires an index to take
   its per-contig parallel path; without one the caller would fall back to the sequential two-pass
   scan, which the same module documents as producing byte-identical output. Either path is
   therefore acceptable; the indexed path is chosen purely for wall time.

---

# RESULT — 2026-08-22, run completed, gate outcome recorded

Written after the run, against the criteria fixed in §6 above. **PROVISIONAL until a verifier
passes it.** Nothing in §1–§6 was edited after the numbers appeared; this section only appends.

## R1. The run

| | value |
|---|---|
| BAM | `pbmc4k_possorted_genome_bam.bam`, 33,243,581,299 bytes (byte-exact vs `content-length`), md5 `52e86c3f1ff853bd9019e9b159ecb526`, `samtools quickcheck` OK |
| code | `/mnt/ssd1/Projects/PeakATail_wd/tools/pa-polya-run-9dfdefb3`, asserted at runtime — `code: …/tools/pa-polya-run-9dfdefb3/ema/strategies/clip_seeded.py` |
| start / end | 2026-08-22 00:07:45 → 00:23:26, `ema exit 0` |
| **wall time** | **15:40.87**, 627 % CPU (`--threads 8`) |
| **peak RSS** | **3,243,632 kB ≈ 3.24 GB** (PBMC 10k v3 arm: 32:59 / 11.12 GB at `--threads 16` — different thread count and a smaller BAM, so this is a run record, not a benchmark) |
| peak calling | 224 (contig, strand) jobs through 8 workers; largest `1+` with 37,349,233 mapped reads |
| outputs | `/mnt/ssd0/emaout/peakatail_benchmark/donor2/` (symlinked as `results/benchmark_tools/pbmc4k_donor2`), 3.7 GB |
| call counts | `pas.bed` 166,816 \| tier1 62,934 \| tier2 103,882 \| tier1 ∩ ≥2 mol **20,681** |

**The "nothing was retuned" claim is now a measurement, not a promise.** A field-by-field diff of the
two runs' `run_config.json` (`variables` + `filters` + `args` + `atlas_distance`) compares **108
recorded parameters**: **107 are identical and exactly one differs — `seqlen` 91 → 98**, the single
library parameter declared in §3. `--threads` is not captured in that block (it records `None` in
both runs) and is stated separately: 16 for PBMC, 8 here, output-irrelevant per §3.3. The caller's
*own* independent inference agrees with the declared value: `Auto-detected median read length 98 bp`.

Detected-gene denominator for this library, by the §4 recipe: **13,763 genes → 268,097 of 569,005
atlas sites** (md5 `3a44e78737f025f3e4de7034e0b114fa`). PBMC 10k v3: 14,949 → 285,136.

## R2. Scores — pre-registered default and the ≥1-molecule arm

Same scorer, same atlas, same nulls, same cutoffs. `n` is `n_matched` (scored calls), the quantity
[19](19_final_gate_v2.md) quotes; raw BED lines in brackets.

| dataset | arm | n | P@10 | P@25 | P@100 | R_det | F1_det | null P@100 (3 seeds) |
|---|---|---|---|---|---|---|---|---|
| **pbmc4k (donor 2)** | **pre-registered default** | **20,672** [20,681] | 0.6208 | 0.7614 | **0.8279** | 0.1104 | 0.1948 | 0.0216 / 0.0222 / 0.0217 |
| **pbmc4k (donor 2)** | ≥1 molecule (tier-1) | 62,911 [62,934] | 0.3419 | 0.4312 | 0.5042 | 0.1784 | 0.2635 | 0.0223 / 0.0227 / 0.0218 |
| pbmc_10k_v3 (donor 1) | pre-registered default | 46,524 | 0.5209 | 0.6389 | 0.7062 | 0.1754 | 0.2811 | 0.0217 / 0.0217 / 0.0217 |
| pbmc_10k_v3 (donor 1) | ≥1 molecule (tier-1) | 167,565 | 0.2298 | 0.2880 | 0.3520 | 0.2685 | 0.3046 | 0.0216 / 0.0219 / 0.0223 |

Secondary arms on pbmc4k (emitted and scored by the same loop): full `pas.bed` n=166,769,
P@100 0.2408, R_det 0.2173; tier-2 n=103,858, P@100 0.0813, R_det 0.0390. TES precision @100 for the
default is 0.5024 (PBMC 0.3307).

Nulls are ~0.022 on every arm and both donors — indistinguishable from the PBMC arm's, as they must
be, since the shuffle uses the same inclusion file. The default arm sits **38×** its own null.

**Independent recomputation.** P@10/P@25/P@100 for the pbmc4k default were recomputed from
`pas_PRESPEC_precision_default.bed` with a hand-written `bedtools closest -s -d -t first` pipeline
that does not call `score_tool.py`: **n=20,672, 0.6208 / 0.7614 / 0.8279 — identical to four decimal
places.** The arm file was also verified to be byte-identical to `awk -F'\t' '$5>=2' pas.bed`, i.e.
the same rule `stage2_final_launch.sh` applies to the PBMC arm, and `peak_filters:
ip_filter=True(mode=filter)` confirms the IP rule was active.

## R3. Cross-donor site concordance (§5 protocol, both directions, with nulls)

Computed from `score_tool.py`'s *own* point files and *own* shuffled nulls, so the point definition,
genome filter and null construction are identical to the scoring above.

| arm | direction | n_query | within 25 bp | within 100 bp | null @100 (3 seeds) |
|---|---|---|---|---|---|
| **default** | pbmc4k → pbmc_10k_v3 | 20,672 | **0.8132** | **0.8423** | 0.0024 / 0.0026 / 0.0020 |
| **default** | pbmc_10k_v3 → pbmc4k | 46,524 | **0.3613** | **0.3993** | 0.0011 / 0.0010 / 0.0011 |
| ≥1 molecule | pbmc4k → pbmc_10k_v3 | 62,911 | 0.4997 | 0.5755 | 0.0087 / 0.0084 / 0.0087 |
| ≥1 molecule | pbmc_10k_v3 → pbmc4k | 167,565 | 0.1878 | 0.2395 | 0.0033 / 0.0030 / 0.0032 |

**84.2 % of donor 2's default calls are reproduced in donor 1 within 100 bp — about 361× the genic-shuffle
null.** The reverse direction is 39.9 %, and the asymmetry is arithmetic, not disagreement: donor 1
emits 46,524 default calls to donor 2's 20,672, so most of donor 1's extra calls have no counterpart
in the smaller set. Both directions are reported, as pre-registered; neither is quoted alone. This is
the human analogue of the cross-mouse reproducibility panel the paper previously could not show.

## R4. Verdict against the §6 criterion

> **PASS. P@100 = 0.8279 on the pre-registered precision default, against a gate of ≥ 0.50.**

The gate now holds on **four independent libraries**: pbmc_10k_v3 0.7062, pbmc4k 0.8279, testis
mouse 1 0.7450, testis mouse 2 0.7572.

On the raw quantity the gate is defined over, this lands in §6's *first* bracket: 0.8279 exceeds
both the 0.50 gate and donor 1's 0.7062, so the "materially below 0.7062" bracket is not triggered.
What the paper must nonetheless say, because it is what the data show:

1. **Precision generalises; it does not improve.** At the shipped operating points donor 2 is
   uniformly higher at every cutoff and on every arm (default 0.8279 vs 0.7062; tier-1 0.5042 vs
   0.3520; full set 0.2408 vs 0.2044). But those are *different operating points*: the fixed rule
   emits 20,672 calls on donor 2 against 46,524 on donor 1. **At matched call count donor 1 is the
   more precise library, not donor 2** (see the matched-N control below). The defensible claim is
   that the precision-first behaviour is not an artefact of one donor — *not* that the tool is more
   precise on donor 2.
2. **Recall is lower on donor 2** (R_det 0.1104 vs 0.1754 on the default; 0.1784 vs 0.2685 at ≥1
   molecule), and the caller emits half as many calls (166,816 vs 334,005). pbmc4k is a 4,340-cell v2
   library against a ~10k-cell v3 library. **This is not a denominator artefact, and the direction
   proves it:** donor 2's detected-gene denominator is *smaller* (268,097 vs 285,136 sites), which
   inflates R_det, and donor 2's R_det is lower anyway. The operating point moved along the
   precision/recall trade-off toward precision; it did not degrade.
3. **The ≥2-molecule default is again not the F1 optimum** — on donor 2, ≥1 molecule gives F1_det
   0.2635 against the default's 0.1948. This replicates on a second human library the finding
   [22](22_performance_roadmap.md)/[23](23_algorithm_roadmap.md) established on PBMC and both mice,
   and strengthens the requirement that the paper present ≥2 molecules as a *reliability* choice
   rather than a maximum.
4. **The evidence channel was not weaker on donor 2, so it cannot be invoked as an excuse for the
   recall drop.** The frozen tree's own QC (`check_clip_rate`, first 200,000 CB reads) measures
   **3.6965 %** on pbmc4k versus **2.2565 %** on pbmc_10k_v3 — 1.64× *higher*, and both far above the
   0.3 % warn threshold.

### Matched-call-count control (added by the verifier, 2026-08-22)

Not pre-registered — added because [22](22_performance_roadmap.md) established for this project that
precision compared at unequal call counts is a property of the operating point, not of the tool, and
the same standard must apply when the comparison flatters us. Donor 1's default arm was cut down to
donor 2's call count by its own clip-molecule score (the tool's own ranking, no atlas involved), and
by molecule-threshold rules that need no tie-break at all:

| donor 1 default, reduced to ≈ donor 2's n | n | P@100 | R_det |
|---|---|---|---|
| top 20,672 by clip-molecule rank (coordinate tie-break) | 20,672 | **0.9049** | 0.1135 |
| top 20,672 by clip-molecule rank (reverse tie-break) | 20,672 | **0.9049** | 0.1134 |
| ≥ 4 molecules | 23,679 | 0.8856 | 0.1241 |
| ≥ 5 molecules | 20,320 | 0.9074 | 0.1122 |
| **donor 2 default (for comparison)** | **20,672** | **0.8279** | **0.1104** |

Every route to a matched call count puts donor 1 at P@100 0.886–0.907 and R_det 0.112–0.124 — i.e.
**donor 1 dominates donor 2 on both axes at matched N**, robustly to the tie-break (3,359 calls sit
at the boundary score of 4, so the tie-break was checked in both directions).

The correct reading is therefore: the pre-registered rule lands at a *more conservative point on the
trade-off curve* on the shallower v2 library, which raises its headline precision and lowers its
recall. The gate passes comfortably and the reliability behaviour reproduces — but donor 2 is not the
better library, and the paper must not say precision "improves" on it.

### Denominator-construction control (pre-declared in §4)

The §4 deviation was measured on PBMC, where both runs exist. Rebuilding the PBMC denominator from
the `clip_seeded_final_v2_ipfilt` run instead of the native `lambda_gradient` run gives **15,271
genes / 288,740 sites** versus the shipped **14,949 / 285,136** — a **+1.26 %** larger restricted
atlas, with **283,774 sites in common (99.5 % of the shipped set)** and 14,813 genes shared. The
construction choice is therefore worth ~1 % of denominator size, far too small to affect any
statement above, and it biases *against* the arm it is applied to (a larger denominator lowers
R_det). Reported as promised, whatever it said.

**Denominator swap (added by the verifier).** The stronger form of the same check is to score both
arms against *both* denominators, which removes denominator composition from the comparison entirely:

| default arm | vs donor 1 denom (285,136) | vs donor 2 denom (268,097) | vs clip-seeded-built PBMC denom (288,740) |
|---|---|---|---|
| donor 1 | **0.1754** | 0.1859 | 0.1735 |
| donor 2 | 0.1036 | **0.1104** | — |

The recall gap survives every denominator choice and is *wider*, not narrower, on any common
denominator (0.1036 vs 0.1754 on donor 1's; 0.1104 vs 0.1859 on donor 2's). §R4.2's conclusion is
confirmed and, if anything, understated: the recall drop on donor 2 is real and is not a denominator
artefact.

### Scope this licenses, and what it does not

The paper may now say the reliability claim holds on **two human 10x 3' libraries spanning two
chemistries (v2 and v3) and two CellRanger versions, plus two mice** — with the §2 caveat kept
verbatim: 10x publishes no donor identifier for `pbmc4k`, so this is demonstrably a second *library
and chemistry*, and only presumptively a second *individual*. It does not license "second donor"
without that qualification, and it does not license any claim about v2 vs v3 as such, which would
need its own pre-registered design with more than one library per chemistry.

### Caveats

- Everything here is PROVISIONAL until a verifier passes it.
- `n_ip_flagged` appears once per run in each log (pbmc4k 54,664; pbmc_10k_v3 110,000) and does not
  reconcile with the 68,855 recorded in the verified IP-rule comparison, so the two libraries' IP
  flag counts are **not** compared here; the field would need to be traced before use.
- Wall time and RSS are not comparable across the two arms (different thread counts, different BAM
  sizes) and are recorded only as a run log.

---

## V. Verifier record — 2026-08-22, adversarial verification

**Verdict: FIXED.** The gate result is trustworthy and the run is genuinely un-retuned. One
interpretive claim was overstated and has been corrected in place (§R4.1 + the matched-N control);
two controls were added that strengthen the builder's own conclusions.

Independently reproduced, not taken on report:

- **Provenance.** `tools/pa-polya-run-9dfdefb3` is at `9dfdefb3eb…` with a *clean* working tree;
  the run log's imported module path resolves inside that snapshot.
- **Parameter identity.** A flattened field-by-field diff of the two `run_config.json` files gives
  **114 fields, 110 identical, 4 different**: input BAM path, output path, timestamp, and
  `variables.seqlen` 91 → 98. Exactly one *knob* differs, as claimed. `read.py:106` confirms
  `span > seq_len → read dropped`, so 98 is a library descriptor, not a tuning choice; and
  `stage2_final_launch.sh` already ran both mouse arms at 98 from the same launcher that ran PBMC
  at 91. `args.seq_len` is `null` in both — the value is inferred, and the caller's own log agrees
  (`Auto-detected median read length 98 bp`).
- **Denominator.** Rebuilt from scratch: 13,763 genes / 268,097 sites, **md5-identical** to the
  builder's (`3ab0724a…`, `3a44e787…`). PBMC's shipped denominator also reproduces md5-identically
  (`bc4b56f2…`) from the source named in `shared_refs_pbmc/README.md`. The pre-declared control is a
  genuine like-for-like: its gene list is md5-identical to a rebuild from the *same run class*
  donor 2's denominator came from.
- **Scores.** Re-scored with a hand-written `bedtools closest` pipeline that never calls
  `score_tool.py`. All four arms match to 4 dp — donor 2 default 20,672 / 0.6208 / 0.7614 / 0.7995 /
  0.8279 / R_det 0.1104 / F1 0.1948; donor 2 tier-1 0.5042 / 0.1784 / 0.2635; donor 1 default
  0.7062 / 0.1754 / 0.2811; donor 1 tier-1 0.3520 / 0.2685 / 0.3046. My independently built point
  files are byte-identical to the scorer's, and `pas_PRESPEC_precision_default.bed` is byte-identical
  to `awk '$5>=2' pas.bed`.
- **Concordance.** Recomputed **both** directions from my own point files: 0.8132 / 0.8423 and
  0.3613 / 0.3993 — exact. The asymmetry is arithmetic as claimed: with only 20,672 donor-2 calls the
  ceiling for the reverse direction is 44.4 %, and 39.9 % is 90 % of that ceiling.
- **Clip rate.** Re-measured with the frozen tree's own `check_clip_rate`: **3.6965 %** vs
  **2.2565 %**, exact. The evidence channel is genuinely richer on donor 2.
- **BAM compatibility.** `@SQ` diff is empty (194 contigs). The header confirms a distinct library
  (flow cell `H53GNBCXY`, `GRCh38-1.2.0` refdata, STAR 2.5.1b).
- **Hygiene.** Nothing written to `/mnt/ssd2`; no `results/benchmark_tools/*final*` path modified;
  donor 1's reference arm and `shared_refs_pbmc/` carry their original mtimes; `--threads 8` is
  within the 12-thread cap; no commits (file is untracked).

**Gate: PASS, and it is the pre-existing gate.** `P@100 = 0.8279 ≥ 0.50`. The threshold was not
invented for this test — it is quoted verbatim from [19](19_final_gate_v2.md), which is committed and
already verified, so it demonstrably predates donor 2 and cannot have been fitted to 0.8279.

**One limitation the verifier cannot close.** Because no commit was made (correctly, per the standing
constraints), the filesystem records only the *last* write to this file (00:27:38, after the run), so
mtime alone cannot prove §1–§6 went unedited once the numbers existed. What *is* established: the
driver `run_donor2.sh`, frozen at 23:48:30 — 19 minutes before the caller started at 00:07:45 —
names this file and its 23:47 write time; the scoring scripts and the pre-declared denominator
control were all frozen 23:47–23:50; and the gate value itself is inherited from a committed
document. Closing it fully would need this file committed at write time.

**Corrections made by the verifier:** §R4.1 rewritten (precision *generalises*; it does not
*improve* — at matched call count donor 1 dominates on both axes), the matched-call-count control
added, and the denominator-swap table added to the §4 control. Nothing else in the RESULT sections
was altered; no number reported by the builder was found to be wrong.
