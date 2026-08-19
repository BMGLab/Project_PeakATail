# Caller fix plan — poly(A) evidence, blocking bugs, and the one re-run we pay for

**To:** Ebru (PI), Amir · **Date:** 2026-08-19 · **Branches:** manuscript `reorg-manuscript`, tool `fix/cellranger-input-compat`
**Decision requested:** approve Stage 1 Phase 2 (≈4 days of caller work + one ~3.5 h benchmark re-run), or take the fallback in §5.
**Companion docs:** `09_headtohead_results.md` (the negative result), `06_roadmap.md` (positioning), `github/issue7_fdr_calibration.md`, `github/issue8_cellranger_cb_suffix.md`.

---

## 1. Situation

The head-to-head put PeakATail **last among de novo callers on PBMC — F1@100 0.134 vs polyApipe 0.261**. We tested and
excluded every cheap explanation: it is not a threshold artefact (ranking the 277k calls by UMI count and scoring the top
22k/36k/106k gives precision 0.119/0.120/0.118, identical to the full set), not a coordinate offset (relaxing to 200 bp lifts
PeakATail 0.118→0.166 while polyApipe goes 0.380→0.410), and not distance-to-TES (median 14.1 kb for PeakATail vs 13.1 kb for
polyApipe on the same BAM). The root cause is now **measured, not inferred: PeakATail contains no read-level poly(A) evidence
anywhere, by construction.** `read_check()` (`ema/countmatrix/read.py:15-131`) discards the `AlignedSegment` at line 131; CIGAR
and SEQ are never inspected in the package (`grep -rn 'soft.?clip|polyA' ema/` returns nothing read-level), and the strategy
contract `find_pas(peak)` only ever sees coverage shape plus barcodes — there is no channel through which orthogonal evidence
could reach the caller. The cost of that, quantified on the PBMC BAM: **1.152%** of CB-bearing reads carry a non-templated
poly(A) soft clip (wrong-end control fires on 0.0125% → **92× specificity**), and **73.9%** of those clips fall within 100 bp of
a PolyASite 2.0 site; the same statistic for *any* read 3′ end is **13.9%** — statistically indistinguishable from PeakATail's
measured precision of **0.1175**. In one sentence for the PI: *a PeakATail call currently carries no more PAS information than a
randomly chosen read 3′ end.* **82.5% of its 277,164 PBMC calls have zero poly(A) read support within 100 bp.** Three orthogonal
correctness bugs (matrix keying, switch-length strand, CellRanger input) independently invalidate every count-dependent claim,
so Stage 0 has to happen whatever we decide about the caller.

---

## 2. THE DECIDING NUMBER

> ### P@100 = **0.375** for the poly(A)-supported subset of the *existing* calls (≥1 clip read within 100 bp, n = 48,506)
> versus **0.118** today — a **3.2×** improvement, level with polyApipe's 0.380.
> **The number is real and clearly better than 0.118. It is also not sufficient.**

**Why "not sufficient": filtering buys polyApipe's precision at half its recall.** Full grid (scored with `score_tool.py`'s exact
metric definitions, detected-gene denominator, `grid.tsv`):

| gate on existing calls | n | P@100 | R@100 | **F1@100** |
|---|---:|---:|---:|---:|
| none (shipped today) | 277,164 | 0.1175 | 0.1546 | **0.1335** |
| ≥1 clip read @100 bp | 48,506 | **0.3745** | 0.0979 | **0.1552** ← best pure filter |
| ≥2 clip reads | 22,373 | 0.5978 | 0.0773 | 0.1369 |
| ≥3 clip reads | 14,126 | 0.7673 | — | 0.1213 |
| ≥5 clip reads | 9,141 | 0.9038 | — | 0.1002 |
| ≥10 clip reads | 6,086 | 0.9525 | — | 0.0754 |
| *polyApipe (target)* | *120,916* | *0.380* | *0.199* | ***0.261*** |

**No pure-filter variant reaches F1 0.16.** Shipping the filter alone would produce a headline precision of 0.375 (or 0.90 at
≥5 reads) while F1 moves 0.134→0.155 — i.e. it would repeat exactly the overclaim the benchmark just exposed. **Do not ship
Phase 1 alone.**

**The number that justifies building: clip-*seeded* calling, not clip-*filtered* calling.** ~30 lines of untuned single-linkage
clustering of clip sites (25 bp, read-weighted mode, no coverage model, no internal-priming filter, no tuning), complete scan of
chr19+chr21:

| arm (chr19 + chr21, same scoring) | n | P@100 | R@100 | **F1@100** |
|---|---:|---:|---:|---:|
| **PeakATail as shipped** | 12,889 | 0.1614 | 0.1962 | **0.1771** |
| clip clusters, ≥1 read | 14,889 | 0.3218 | 0.3366 | **0.3291** |
| clip clusters, ≥2 reads | 5,678 | 0.5065 | 0.2339 | **0.3200** |
| clip clusters, ≥2 UMIs | 4,498 | 0.5887 | 0.2228 | **0.3232** |

**1.86× PeakATail's F1 on the same chromosomes, untuned, and above polyApipe (0.261) and scUTRquant (0.290) on that slice.**

**The ceiling, so nobody over-promises.** Of the 285,220 detected-gene atlas sites, only 28.77% carry ≥1 clip read within 100 bp
(19.86% ≥2, 15.83% ≥3). polyApipe's measured recall of **0.199** lands exactly on the ≥2-read ceiling of **0.1986** — an
independent confirmation that this measurement is calibrated. So a strictly clip-gated caller tops out at **R@100 ≈ 0.29,
F1 ≈ 0.42**. Any recall beyond that must come from a hybrid (clip-supported calls **plus** tagged coverage-only calls) reported
as **two explicit tiers**, or we are back to overclaiming.

**Where the evidence physically is** (exact counts, not extrapolated): on chr19, 87.7% of clip reads lie within ±150 bp of a
detected-gene atlas site but only **41.5%** lie near any PeakATail call. **~60% of the poly(A) evidence falls outside every
current peak window** — which is precisely why no filter can recover it and why Phase 2 must *seed* candidates from clips.

**Verdict: build it, as Phase 2. Phase 1 is scaffolding and a reviewer-facing evidence column, not the fix.**

---

## 3. Staged plan

| Stage | Content | Owner | Depends on | Cost |
|---|---|---|---|---|
| **0** | Blocking correctness bugs (matrix keying, switch-length strand, CellRanger branch is RED, no CI) | Amir (branch/CI/fixture), Claude (library fixes, under review) | — | 1–2 days, no re-run |
| **1** | Read-path evidence: poly(A) Phase 1 + Phase 2, UMI dedup + flag filtering, `read.py` coordinate defects | Claude → Amir review | Stage 0 green | 4–5 days |
| **2** | One benchmark re-run + re-score, all three arms concurrently | Claude | Stages 0+1 both landed | ~3.5 h wall + minutes |
| **3** | Re-test the three manuscript claims on re-keyed artifacts | Claude → Ebru | Stage 2 | 1–2 days |

**Budget exactly ONE full caller re-run.** The accuracy claim is scored from `pas.bed` **coordinates only**, so matrix-keying,
switch-length and FDR fixes do not invalidate F1@100; conversely novelty (AMI) and statistics both consume the count matrix, so
they *are* invalidated by any read-path or keying change. Landing Stages 0+1 together and re-running once is the whole
scheduling trick.

### Stage 0 — bug fixes that block correctness (no re-run needed)

| # | Bug | Evidence | State |
|---|---|---|---|
| 0a | **Matrix mis-keying.** `make_dataframe()` returns a full-height matrix with a *compacted* `pas_ids`; `annotate()` indexes rows by **position** in that list. Row `pas_id-1` is correct; `pas_id_to_row[pas_id]` is correct only if no PAS row is empty. | **45,860 / 45,921 = 99.87%** of annotated PAS on testis mouse1 carry another PAS's counts. Depth correlation (Spearman, chr19, n=1478): **0.016 shipped → 0.718 re-keyed**. Regression bisected to `04e0b3a` (2026-03-23): **everything produced after that date is affected.** | Fix + `tests/test_matrix_pas_id_row_alignment.py` written and validated in the working tree (suite 39F/871P vs 39F/869P on HEAD; +2, no regressions). **Uncommitted.** |
| 0b | **Switch-length strand inversion**, residual paths. `per_gene` was fixed at HEAD (`02f4153`), but `_build_gene_fallback_map` (`runner.py:121`) still hardcodes the `(gene,"_gene_",0,1,1)` sentinel used by `--isoform-agg per_isoform --utr-unmatched gene`. | Shipped AT2 table: `per_gene` pre-fix **2299/2299 = 100%** of minus-strand genes inverted (50.1% of all PDUI genes); `_gene_` fallback rows **1667/1667 = 100%** inverted, still on HEAD. Also: `switch length` has no `--pasbed` flag and silently reverts to plus-strand convention when the walk-up misses. | Diagnosed; `tests/test_switch_length_strand_and_counts.py` written (2 pass, 3 `xfail(strict)`). Fixes **not applied**. |
| 0c | **PDUI computed on TF-IDF weights, not counts.** `build_count_dfs` reads only `.X`, which `leiden_tfidf` overwrites in place; no `layers`, no `raw` anywhere in the package. | **100.0%** of non-zero `proximal_reads`/`distal_reads` in shipped tables are non-integer. On informative cells (both ends non-zero): r = **0.47**, **33.4%** cross 0.5; **10.6%** of cluster-pair ΔPDUI signs flip; top-200 switching-gene lists overlap only **64.7%**. | Fix known: stash `layers['counts']` before normalisation (`ad.concat` preserves layers — verified). Backfill possible from `06_preprocessing/*/preprocessed.h5ad` (int64, identical obs/var — verified). |
| 0d | **`fix/cellranger-input-compat` is RED.** 11 failed vs 7 on `develop` — 4 new failures all from `_MockRead` lacking `is_unmapped`; two fail *again* once that is added, because their RG values contain underscores the branch now rewrites. | Measured on both branches. Of the 7 pre-existing, 5 are real drift and 2 are environmental (no editable install). | **Do not merge or benchmark from this branch until green.** |
| 0e | **RG fix is a data-corrupting regression for merged BAMs.** `eb13529` *discards* any RG containing `_` instead of sanitising it — and `ema merge` (`samtools merge -r`) derives RG IDs from filenames. | Two samples merged from `sampleA_rep1.bam`/`sampleB_rep1.bam` collapse to one RG → identical CBs share a matrix column. | Replace with a bijective sanitiser (`rg.replace('_','-')` or a per-run RG→short-id table). |
| 0f | **No CI, no BAM fixture.** `.github/workflows/` has only `docs.yml`; PR #73 merged with no test gate. 8 of 9 skips are "no BAM available"; `.gitignore:15` (`*.bam`) guarantees they stay skipped; `test_region_fetch.py` points at `/home/user/PeakATail/test_run/chr22.bam`, a path from another machine. | — | Prototype fixture built and verified: 148 kB / 2072 reads, replays all three input bugs at the exact commits that introduced them. Needs `!tests/fixtures/*.bam` negation to be committable. |
| 0g | **Harness-side mis-pairing that survives the library fix.** `benchmark_headtohead.py:343-347` zips coordinate-sorted `pasbed.bed` names against numerically-ordered matrix rows; the `assert len(names)==len(sums)` passes, hiding it. | — | **FIXED 2026-08-19** (`depth_map()` re-keyed on `annotated_pas_ids.tsv`, length assert replaced by set equality; `tests/test_headtohead_depth_keying.py` 2F→0F). Blast radius measured: 100% of mouse PAS, 99.80% of PBMC PAS carried another PAS's depth. The top-N ranking arms were built by the same defect and are re-scored (0.120/0.118/0.118 vs the published 0.119/0.120/0.118; selection overlap only 7.89% at top-22k). Deltas tabulated in `05_figure_index.md` §8 under PENDING RE-VERIFICATION; §1 of THIS doc quotes the stale 0.119/0.120/0.118 and needs the same edit. |

**Definition of done:** suite green on `fix/cellranger-input-compat` (0 new failures vs `develop`; the 5 pre-existing drift failures fixed or explicitly xfailed with a tracking issue); `tests/fixtures/cellranger_pbmc_tiny.bam` + `.gitignore` negation + `test_cellranger_input_compat.py` committed; `ci.yml` (rebased `feat/release-engineering`, verified conflict-free) running on both `develop` and `develop→main` PRs; 0a/0b/0c fixes merged with their tests flipped from `xfail` to pass.
**Test that proves it:** the three xfail(strict) tests in `test_switch_length_strand_and_counts.py` become failures-if-unfixed; `test_matrix_pas_id_row_alignment.py` fails on HEAD and passes after; the fixture replay shows 845 barcode columns instead of 1.
**Stop signal:** if repairing the CellRanger branch turns out to need changes to `read_check`'s 5-tuple, stop and fold it into Stage 1 — the tuple is shared by the monolithic loop, the 3-stage pipeline (`peak_pipeline.py:151`) and the tile workers, and widening it in isolation produces silently divergent output between `--tiles`, `--pipeline` and the default path (the benchmark used the monolithic path, so a tile regression would not be caught by re-running the benchmark).

### Stage 1 — poly(A) evidence in the caller (+ the read-path fixes that share the re-run)

**Phase 1 — instrument and expose (1 day, no default behaviour change).**
New `ema/countmatrix/polya.py::clip_site(read, min_clip=6, min_purity=0.8)`; call it after the `read_check` success at
`peackcalling.py:302` behind `if polya_enabled:`; accumulate via a new `Peak.polya_counting(site, cb)` mirroring
`cb_position_counting`; score each emitted PAS by clip reads and distinct UMIs within `±polya_window` of its strand-aware 3′
base; write the count into **BED column 5 of `pasbed.bed`**, which `paswrite.py` currently hardcodes to `0` and every downstream
reader already names `score` and ignores — so the column is free and `annotatedpas.bed` inherits it automatically. CLI flags
alongside the existing `ip_*` block; `polya_mode='filter'` enforced at the existing `_apply_pas_filters()` seam
(`ema/main.py:364-500`). **Cost: +0.54 µs/read measured over 1.5M real reads ≈ +5% of peak-calling wall time, no second BAM
pass.** Do **not** widen `read_check`'s return tuple (see 0-stop-signal).

**Phase 2 — clip-*seeded* calling (3–4 days; this is the actual fix).**
New registered strategy `@register('polya_cluster')` so nothing existing regresses, plus one caller-level change: PAS candidates
must be seeded from clip-site clusters, not only coverage summits. (a) accumulate clip sites into the `Peak`; (b) `find_pas`
returns the read-weighted modal position of each 25 bp single-linkage cluster clearing `polya_min_reads` (**UMI-deduplicated
support scores best: P@100 0.676 vs 0.598 for raw reads at comparable n**); (c) coverage summits with no clip cluster are still
emitted, tagged `polya_support=0` — this is the second tier; (d) emit the cluster **mode** as the 3′ base rather than
`partition_peak_region`'s outer edge (**P@25 0.34→0.73 at ≥3 reads, free**); (e) split multi-PAS UTRs between clusters using the
coverage model — the one place PeakATail's coverage machinery genuinely beats polyApipe.

**Same wave, same re-run — do not defer:**
- **`read.py:124-127` drops every read whose reference span exceeds `--seq-len`**: 13.14% of reads in PeakATail's own peak
  windows, and 12.65% contain an N (junction) op — essentially every splice-junction read is silently discarded.
- **Short reads are padded** (`read_end = read_start + seq_len`), pushing a poly(A)-clipped read's 3′ end downstream by exactly
  its clip length (median 15 nt). The tool systematically moves the very reads that mark the cleavage site away from it.
- **No UMI dedup, no SAM-flag filtering anywhere in the package** (`UB`, `is_duplicate`, `is_secondary`, `is_supplementary`,
  `mapping_quality` appear in no source file). Measured on the PBMC BAM: 51.2% PCR duplicates in `1:150M-170M`, 44.5% secondary
  in `1:1-20M`, and 97.5% of secondaries carry a CB so all of them are ingested. A duplicate stack is exactly the sharp
  single-coordinate spike `Peak.pasfind` calls as a PAS. Downstream docstrings already claim "integer UMI counts"; they are raw
  read counts. **Cheapest large-effect accuracy fix on the list** — and it halves Fisher's effective n, feeding straight into
  issue #7.
- Turn on `--ip-filter` with `--genome-fasta` in the same pass: a genomic A-run downstream of a clip cluster is this evidence
  type's one systematic false-positive mode. **Blocker: no human GRCh38 FASTA is staged under `data/references` (mouse GRCm38
  only) — Amir, please stage one before Stage 2.**

**Definition of done:** `--peak-strategy polya_cluster` runs end-to-end on PBMC; BED column 5 carries per-PAS clip/UMI support;
default path byte-identical to today when `polya_evidence` is off; tile/pipeline/monolithic paths produce identical output.
**Test that proves it:** the chr19+chr21 slice reproduces F1@100 ≥ 0.30 through the real caller (the standalone prototype gives
0.329); UMI dedup is verified by a fixture with a known duplicate stack collapsing to one count.
**Stop signal:** if the in-caller chr19+21 F1 lands below **0.25** (i.e. materially under the standalone prototype's 0.329 and
below polyApipe's 0.261), stop before the full re-run — the coverage machinery is fighting the clip seeds, and the honest report
is the benchmark paper (§5), not another week of tuning.

### Stage 2 — re-run and re-score (one pass, ~3.5 h wall)

```bash
export LC_ALL=C
# all three arms CONCURRENTLY (box: 1005 GB RAM, 112 cores)
bash scripts/benchmark_tools/run_peakatail.sh                                  # PBMC: 3 h 27 m, 105 GB peak RSS
bash results/benchmark_tools/gse104556/peakatail/run.sh                        # testis m1 2 h 06 m / 14.7 GB, m2 1 h 26 m / 14.3 GB

# re-score (minutes)
python3 scripts/benchmark_tools/score_tool.py \
  results/benchmark_tools/pbmc_10k_v3/peakatail/pas.bed peakatail_polya \
  --detected-atlas results/benchmark_tools/shared_refs_pbmc/pas2.in_detected_genes.bed
python3 scripts/benchmark_tools/score_tool.py <mouse pas.bed> peakatail_polya_mouseN --outdir <tooldir> \
  --atlas data/references/atlases/polyasite2.GRCm38.96.rep_sites.bed6 \
  --tes data/references/atlases/tes.protein_coding.GRCm38.102.bed6 \
  --genome data/references/mouse/chrom.sizes.filt \
  --genebodies results/benchmark_tools/gse104556/shared_refs/genebodies.merged.bed \
  --detected-atlas results/benchmark_tools/gse104556/shared_refs/pas2.in_detected_genes.bed
```

Downstream stages replay from `filterdmatrix.mtx` (which is correctly keyed) without redoing peak calling.

**Definition of done — the publication gate. Ship the claim only if, on PBMC:**

| metric | gate | today | prototype predicts |
|---|---|---|---|
| **F1@100** (detected-gene atlas) | **> 0.261** (polyApipe) | 0.134 | ~0.33 naive, 0.40+ with coverage + IP filter |
| **P@100** | **≥ 0.38** | 0.118 | 0.32–0.59 depending on gate |
| clip-supported tier reported separately from coverage-only tier | mandatory | n/a | — |
| testis (STARsolo) arm | clip rate > 0.3% of CB reads, else report human-only | n/a | unmeasured — **validate first** |

**Stop signal:** F1@100 ≤ 0.261 on PBMC after Stage 1 → do not write an accuracy claim; go to §5. Peak RSS above ~150 GB on the
first run → stop and profile before launching the arms concurrently.

### Stage 3 — re-test the three claims

| Claim | Affected by | Action | Cost |
|---|---|---|---|
| **1. Accuracy (F1@100)** | Stage 1 only — it is coordinate-only, `pas.bed`, so keying/strand/FDR fixes do not touch it | Re-score from Stage 2 | minutes |
| **1b. "Ranking by expression does not enrich"** | 0a — `topN/pas_top*.bed` were ranked by row sums of a **mis-keyed** matrix | Re-rank on re-keyed counts, re-score | ~1 h |
| **2. FDR calibration** | 0a (its TREG h5ad is post-2026-03-23) **and** an independent bug: `fisher.py:129-130` `.astype(int)` on TF-IDF floats, so the test's effective n is a sum of IDF-weighted logs, not reads | Re-run nulls on re-keyed, counts-layer h5ad. **Start now — fully decoupled**: fisher 20 perms ≈ 11 min, nb_pairwise 20 perms ≈ 50 min, resumable via `DONE.ok` | ~1 h |
| **3. Novelty / AMI ablations** | 0a fully — the "REAL poly(A) sites vs non-PAS peaks" split partitions by labels that do not match the counts, so **0.32 vs 0.74 is not interpretable** | Re-run `pbmc_novelty_analysis.py` on re-keyed `clusters.h5ad`; also fix `benchmark_headtohead.py` `depth_map` (0g) | ~1 day |

**Stop signal:** if the AMI ablation does not survive re-keying (0.32 vs 0.74 collapses), the clustering-novelty leg of the new
positioning is gone too — escalate to Ebru immediately, because §5's fallback then rests on the benchmark methodology alone.

---

## 4. What to tell Amir now — should the #67 rerun wait?

**Yes. HOLD the SWITCH_CELLTYPE rerun.** The `main.nf` one-liner is correct and should be committed, but running it today buys a
result that must be thrown away — and worse, one that looks plausible.

| If Amir reruns… | What comes back |
|---|---|
| with the ema build that produced the 2026-08-17 outputs | **100% of minus-strand genes inverted** (2299/2299 measured in the shipped AT2 table; minus-strand = 50.1% of all PDUI genes). Every shortening call is a lengthening call. |
| with ema at current HEAD | `SWITCH_CELLTYPE_UTR` passes exactly `--isoform-agg per_isoform --utr-unmatched gene`, which still routes every UTR-unmatched PAS through the `rank=1` sentinel — **1667/1667 inverted** in the existing output. |
| either way | PDUI is computed on **TF-IDF weights, not counts**: r = 0.47 on informative cells, 33.4% cross 0.5, 10.6% of cluster-pair ΔPDUI signs flip, ~35% of the top-200 switching genes change. `--pdui-pseudocount 1.0` is calibrated for read counts and is being added to log-scaled weights. |

**Unblock path (cheap, ~half a day, no full re-run):**
1. Backfill `layers['counts']` into each `clusters.h5ad` from the sibling `06_preprocessing/<ds>/preprocessed.h5ad` — int64
   counts with byte-identical `obs_names`/`var_names`, verified on `GSM3516662-StageIA`. `ad.concat(join='outer',
   merge='first')` preserves layers, so it survives `switch combine`.
2. Land Stage 0b/0c fixes (kill the `_gene_` sentinel; add `--pasbed` to `switch length` and raise instead of silently ranking
   by input order; make `build_count_dfs` layer-aware and **raise** on a non-integral matrix behind
   `--allow-non-count-matrix`).
3. Then only **B2 celltyping + SWITCH_COMBINE + SWITCH_CELLTYPE** re-execute — not the whole sweep.

**Cheap permanent guard, worth adding to CI:** for every row of any `pdui_classic.tsv`, `strand '+' ⇒ proximal_start <
distal_start` and `'-' ⇒ proximal_start > distal_start`. That one awk-level check on real output is what produced the 2299/2299
evidence, and it catches the entire strand-inversion family.

**Amir's critical-path work that is NOT blocked and gates everything else:** Stage 0d/0e/0f — repair
`fix/cellranger-input-compat` (the `_MockRead.is_unmapped` attribute **and** the two underscore-bearing RG assertions, which fail
a second time otherwise), land the BAM fixture with the `.gitignore` negation, and push the rebased `ci.yml` with the 5 known-red
tests explicitly xfailed and tracked. **Do not push `ci.yml` before those 5 are handled** — a CI that is red on day one trains
the lab to ignore it.

---

## 5. Risks, and the honest fallback

| # | Risk | Number | Mitigation |
|---|---|---|---|
| R1 | **The recall ceiling is hard.** Only 28.8% of detected-gene atlas sites carry any clip read within 100 bp in this BAM. | R@100 ≤ 0.29, F1 ≤ ~0.42, ever | Report two tiers (clip-supported / coverage-only) separately from the first draft. Never quote the high-precision tier without its n and recall. |
| R2 | **Chemistry- and aligner-dependence.** 1.15% clip rate here is CellRanger 3.0.0 / 10x v3 / 91 bp R2. The testis arm is **STARsolo, 10x v2, 98 bp R2** and is unmeasured; any pipeline that trims poly(A) before alignment destroys the evidence entirely. | unmeasured on mouse | **Measure the clip rate on the GSE104556 BAMs before Stage 1 Phase 2 starts** (~30 min). Ship a loud startup warning when the observed clip rate is near zero rather than silently emitting an unsupported call set. |
| R3 | **Internal priming is the one systematic false positive of this evidence and is NOT covered by our measurements** — the 92× wrong-end control bounds alignment artefacts, not genomic A-runs. | untuned clip precision 0.32–0.51 | Pair with `--ip-filter`. Requires a human GRCh38 FASTA that is **not currently staged**. Measure the delta; do not assume it. |
| R4 | **Phase 1 shipped alone is a trap** — P@100 0.375 looks like a fix while F1 moves 0.134→0.155. | — | Phase 1 is annotate-only by default and is not a manuscript claim. |
| R5 | **Plumbing divergence.** `read_check`'s tuple is shared by monolithic / pipeline / tile paths; the benchmark used the monolithic path, so a tile-mode regression would not be caught by re-running the benchmark. | — | Add a three-path equality test on the fixture BAM in Stage 0. |
| R6 | **UMI dedup halves effective depth**, and `Peak.pasfind`'s hardcoded `max_height <= 20` reject will then drop many loci. | ~50% of ingested coverage in gene-dense regions | Re-tune that threshold in the same wave, or the yield collapse will be misattributed to the poly(A) filter. |
| R7 | **Every count artifact produced after 2026-03-23 is wrong and has no version stamp.** | 99.87% of rows on mouse1 | Delete rather than selectively trust; regenerate from `filterdmatrix.mtx`. |
| R8 | **Field velocity** — scTail (GB 2025), scPAISO (2025). | — | The fallback below is publishable now and does not decay. |

### The fallback: the benchmark / methods paper

**Trigger — take it if any of these fires:**
1. Stage 1 in-caller chr19+21 F1@100 **< 0.25** (stop before paying for the re-run), **or**
2. Stage 2 PBMC F1@100 **≤ 0.261** or P@100 **< 0.38**, **or**
3. clip rate on the STARsolo testis BAMs **< 0.3%** of CB reads *and* PBMC alone does not clear gate 2 (single-dataset accuracy
   claim is not defensible), **or**
4. the Stage 3 AMI ablation does not survive re-keying (then the clustering leg is gone as well, and the benchmark is all we
   have).

**What that paper is, and it is already ~65% collected:** six tools × two datasets under one identical fair matcher, with
gene-body-constrained 3-seed nulls, per-dataset detected-gene denominators, a replicate-reproducibility metric, the ~90–105 nt
cleavage-offset discovery, the finding that **evidence type — not peak-calling sophistication — is the discriminator**
(quantified here: clip reads 73.9% within 100 bp of a PAS vs 13.9% for arbitrary read 3′ ends), the demonstration that
per-dataset reporting is mandatory (PeakATail spans P 0.118→0.447 across datasets), the environment archaeology every tool
needed to run at all, and two independently reproduced correctness bugs found by adversarial re-analysis. PeakATail appears as
one tool among six, reported as measured — plus the peak-based clustering novelty as an orthogonal axis no competitor offers.
**That paper does not require the caller fix to be a success, and it is the version we can defend under review today.**

---

## Approvals needed

| Question | Who decides |
|---|---|
| Approve Stage 1 Phase 2 (~4 days + one 3.5 h re-run), or take §5 now? | **Ebru** |
| HOLD #67 until 0b/0c land + the counts backfill? (recommendation: **yes**) | **Ebru + Amir** |
| Stage 0d/0e/0f — branch repair, BAM fixture, CI — start immediately? | **Amir** |
| Stage 3 claim-2 FDR re-run — start now in parallel (it blocks nothing)? | **Ebru** |
