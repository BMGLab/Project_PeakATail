## Stage 1b + 1c: count `clip_seeded` tier-1 PAS from the cleavage site, and report clip support in molecules

Base: `develop` (at `6697b0d`, the merge of #92). Head: `feat/polya-evidence`, five commits on top of `f0370f7` — `7a54b01`, `002d976`, `ebb9fb3`, `0eebfe7`, `4efeb12`. 16 files, +2,358 / −150. Nothing in this PR moves a PAS coordinate; it fixes what the caller *counts* and what BED column 5 *means*.

## Summary

- **1b (`7a54b01`)** — tier-1 rows of the count matrix were being filled from the cluster's poly(A)-clipped reads only (~1% of the reads). They are now counted from the reads that pile up at the cleavage site: the suppressed coverage candidate's counts, partitioned between clusters at anchor midpoints, plus the read ends in a cleavage window (`--polya-count-window`, default `auto,25`) that belong to no candidate. No read is counted twice; tier-2 rows are untouched.
- **1c (`ebb9fb3`)** — BED column 5 now carries **distinct `(cell barcode, UMI)` molecules**, the unit the `--polya-min-*` gate always used. `--polya-min-reads` is renamed `--polya-min-umis` (old spelling accepted with a deprecation warning, CLI and YAML). Raw reads, `-F 3844` counts and the matrix-row read count move to a new `pas_support.tsv` sidecar; `pasbed.bed` stays BED6.
- **1c** also adds opt-in `--polya-clip-filter f3844` (default `none`, deliberately — it changes the call set) and closes a double count in the rare tier-1 clip-fallback path.
- Coordinates, ids, strand and tier tags are byte-identical to `develop` at the default gate on every acceptance run; only column 5 and the matrices change.
- Measured on GSE104556 mouse1: 1b restores 1,294/1,294 STARsolo cells (the clip-only matrix kept 995) at 1.053× the correctly keyed shipped gene-assigned mass; 1c changes the raw matrices by −0.0005% and leaves tier-1 scoring identical to 4 dp.

## Why

Three findings from the Stage-2 full-genome runs of the `clip_seeded` caller (#92):

1. **Tier-1 PAS were quantified by clip reads only.** `_flush_seeded` wrote `ClipSeeder`'s clip-only per-CB dict as every tier-1 row, so the annotated matrix carried ~10× fewer counts than the shipped caller and the `min_read` cell filter lost 299 of the 1,294 real STARsolo cells (23%) on mouse1. Every count-based downstream step (novelty, FDR, switch tests) on those runs was void.
2. **The score column and the gate disagreed on units.** `--polya-min-reads` documented "distinct molecules" and `ClipSeeder.flush()` did gate on molecules, but column 5 wrote raw clip *reads* — PCR duplicates and secondary alignments included (`read_check` applies no duplicate/secondary filter). On the PBMC run, 18,864 of the 79,751 "≥2-read" tier-1 sites (23.7%) were one molecule counted twice, and 47.7% of the exactly-2-read sites. A post-hoc "≥2 reads" sweep therefore did not describe what `--polya-min-reads 2` would emit, and its recall "above the ≥2-clip-read ceiling" (0.203 vs 0.1986) was the signature of duplicate counting, not a better caller (`manuscript/12_stage2_gate.md`, CORRECTION section).
3. **The clip fallback double counted.** A cluster whose candidate partition and cleavage window were both empty fell back to all of its clip reads, including ones whose `end1` lay past the midpoint inside the neighbour's window. Pinned at `13:43135248+` (pas#23403, mouse1): 37 clip reads past the clipped window's upper bound, row summed to 371 where candidate + remainder is 334. On mouse1 this touched 5,250 rows and ~8k of 140.7M counts (~0.006%).

## What changed, per commit

**`7a54b01` fix(caller): count clip_seeded tier-1 PAS from the read ends at the cleavage site** — 11 files, +1,030 / −55
- `ema/countmatrix/polya.py`: `ClipStream` records every accepted read end per chromosome/strand as two int32 arrays (`end1`, interned CB id; 8 bytes/read, sorted by construction, released at each chromosome flush, picklable for the pipeline finder → writer hand-off). `ClipSeeder.flush()` does the counting: a cluster that suppresses a coverage candidate takes that candidate's `cb_positions` restricted to its partition (`partition_peak_region` at anchor midpoints, `reconstruct_cb_dict` semantics on a compact `_CandidateSlice`; a non-positional strategy such as `original` is detected and its candidate goes whole to the nearest cluster), plus the read ends in `[site − seq_len, site + 25]` (transcript orientation) that lie in no candidate, midpoint-clipped to neighbours. Minus-strand anchors shifted by `seq_len` into `end1` space. Per-strand mass accounting logged as `clip_seeded counting: {...}`.
- `ema/countmatrix/peackcalling.py`, `peak_pipeline.py`, `tile_runner.py`, `ema/main.py`: the three peak-calling paths feed the same `ClipStream`; `JobSpec`/`run_tiled`/`run_pipeline` thread the new kwarg.
- `ema/cli/config_schema.py`: `--polya-count-window`.
- `tests/test_polya_tier1_counts.py` (new), `tests/test_polya_three_path_agreement.py` (now pins per-PAS counts, not only coordinates).
- `docs/cli/run.md`, `docs/strategies/peak-calling.md` (new step 6, "Counting"), `CHANGELOG.md`.

**`002d976` docs(changelog): measured acceptance numbers for the clip_seeded counting fix** — `CHANGELOG.md` only (+18): the mouse1 v2 numbers quoted under Validation.

**`ebb9fb3` fix(caller): count clip evidence in UMI-deduplicated molecules, not raw reads** — 16 files, +1,287 / −116
- `ema/countmatrix/polya.py`: `molecule_cb()` (drops the `RG_` prefix `read_check` prepends to the barcode), `clip_read_ok()` (== samtools `-F 3844`), both units tracked per site in the accumulator and the seeder; `--polya-clip-filter`; fallback path midpoint-clipped and candidate-excluded; `empty_tier1_rows` stat.
- `ema/countmatrix/paswrite.py`: `SUPPORT_COLUMNS`, `support_path_for()` — the `<bed>.support.tsv` sidecar.
- `ema/countmatrix/peak.py`, `peackcalling.py`, `peak_pipeline.py`, `tile_runner.py`, `ema/main.py`: sidecar written on all three paths; the tile merge re-keys it with the BED's own pasnumber remapping; single-BAM runs merge it to `<run>/pas_support.tsv`; CLI deprecation warning for `--polya-min-reads`.
- `ema/cli/config_schema.py`: `polya_min_reads` → `polya_min_umis` with `cli_aliases=("--polya-min-reads",)` and `legacy_alias="polya_min_reads"` (new `FieldSpec.cli_aliases`); `polya_clip_filter` with `choice=("none", "f3844")`.
- `tests/test_polya_clip_evidence_units.py` (new), `tests/test_polya_three_path_agreement.py` (pins the sidecar across paths), `tests/test_polya_clustering.py` (support-window score is molecules), `tests/test_polya_tier1_counts.py`.
- `.gitignore`: `peakatail_*.log`, `switch_combined/` (CLI run debris at the repo root).
- `docs/cli/run.md` (flag rows + "`pas_support.tsv` — where the raw counts went"), `docs/strategies/peak-calling.md`, `CHANGELOG.md`.

**`0eebfe7` docs(changelog): measured Stage-1c acceptance numbers (testis mouse1 v3)** — `CHANGELOG.md`, `ema/cli/config_schema.py`, `ema/main.py` (+40 / −5): the mouse1 v3 numbers; help/docstring text for `--polya-evidence`, `--polya-window` and `_apply_polya_gate` now says molecules (string literals only).

**`4efeb12` fix(config): warn on deprecated YAML key polya_min_reads (verifier repair, Stage 1c)** — `ema/cli/config_schema.py` (+9): the YAML alias now logs the same DEPRECATED warning the CLI alias does.

## Behaviour changes for users

- **`--polya-min-reads` → `--polya-min-umis`.** Default stays 1. The old flag and the old YAML key `polya_min_reads` still work and log `DEPRECATED`. The gate itself is unchanged (it always counted molecules), so no PAS is added or dropped at the default.
- **BED column 5 of `pasbed.bed` is distinct `(cell barcode, UMI)` molecules**, in both the `clip_seeded` tier-1 score and the per-PAS support annotation of the other strategies. A read with no `UB` tag counts as one molecule; PCR duplicates of one molecule count once; the read-group prefix is dropped from the molecule key (a BAM with one read group per lane no longer counts one molecule once per lane — +4.7% on the PBMC chr21 (+) slice, 13,100 vs 12,483 molecules over identical clusters). `score == 0` still means coverage-only, so `--polya-mode filter` behaves as before. Matrix cell identity is unchanged.
- **New `pas_support.tsv` sidecar** next to every caller BED (`<bed>.support.tsv`) and, for single-BAM runs, at `<run>/pas_support.tsv`: `pas_id`, `clip_reads`, `clip_umis` (== column 5), `clip_reads_f3844`, `clip_umis_f3844`, `window_reads` (reads counted into that row of the matrix), `tier`. It is a superset of the post-filter `pasbed.bed`; join on `pas_id`. No BED parser sees a new column.
- **`--polya-clip-filter {none,f3844}`, default `none`.** `f3844` restricts the molecule count — and therefore the gate, the score and the tier tag — to reads passing samtools `-F 3844`. It also drops every cluster whose evidence is entirely secondary/duplicate/qcfail alignments: −8.3% of chr19 (+) and −27.0% of chr21 (+) tier-1 clusters on PBMC. That is a call-set change and stays opt-in; the `*_f3844` sidecar columns let it be evaluated post hoc without a re-run. `read_check` is untouched (widening it would move every peak).
- **`--polya-count-window UP,DOWN`, default `auto,25`** (`auto` == `--seq-len`). Only affects tier-1 matrix counts, never the BED.
- **Tier-1 matrix rows now carry the cleavage-site pile-up**, so the post-filter `pasbed.bed` can grow: on mouse1 58,848 → 85,993 PAS, because ~27k tier-1 rows now pass `min_cells`; 120 one-to-two-read clusters drop out because their counts moved to cells below `min_read`.
- **A tier-1 fallback row may legitimately be empty** (its clip reads' ends belong to the neighbouring cluster); `empty_tier1_rows` in the per-strand counting log counts them. On mouse1: 562 of the 5,250 fallback rows (98 on `+`, all 464 on `−`).
- Memory: on the chr19+ PBMC slice the counting change peaks at 1.45 GB vs 1.08 GB before (5.5 GB when `Peak` objects were retained instead of `_CandidateSlice`). Full-genome footprint is a known problem — see follow-ups.

## Validation

**Tests** (`pytest`: 1215 passed / 2 failed at `ebb9fb3` vs 1189 / 2 at `002d976`; the 2 are the known environment failures in `tests/test_pyproject_install.py`).
- `tests/test_polya_tier1_counts.py` (new, 1b): a cluster inside a 200-read peak carries the peak, not ~2; the coverage strategy's mass at the locus is conserved; an isolated cluster with 30 in-window + 3 clip reads counts 33 on both strands; `--polya-count-window` narrows only the outside count; tier-2 rows byte-for-byte the coverage strategy's; the BED is invariant under any count window; every tier-1 row has counts and no mass is lost; midpoint partition between two clusters; no reads invented for a non-positional candidate; window count excludes candidate intervals and neighbours; the inside cluster also takes the peak's uncovered tail; minus-strand window shifted by `seq_len`; `parse_count_window` forms.
- `tests/test_polya_clip_evidence_units.py` (new, 1c, 26 collected): three reads of one molecule score 1 and report 3 reads; same UMI in two cells is two molecules; reads without `UB` count one each; `--polya-min-umis` rejects a duplicate stack; molecules keyed by barcode not read group; `clip_read_ok` is exactly `-F 3844` (parametrized over the flag bits); flagged reads excluded from the `_f3844` columns only; `f3844` scores, gates and drops an all-flagged cluster; both spellings and the YAML alias accepted; deprecated kwarg mapped to molecules; fallback excludes reads a neighbour already counted / keeps its own territory / ignores reads inside a candidate; sidecar written and joins the BED by `pas_id` end-to-end; `--polya-clip-filter` keeps coordinates end-to-end.
- `tests/test_polya_three_path_agreement.py`: monolithic, `--pipeline` and `--tiles` must agree on BED output, per-PAS counts (1b) and the support sidecar (1c), for both `lambda_gradient` and `clip_seeded`.

**Acceptance run, Stage 1b** — GSE104556 mouse1, `--seq-len 98`, 12 threads, 55 min wall, 26.9 GB peak RSS (`CHANGELOG.md`; verdict recorded in `manuscript/12_stage2_gate.md`):
- every caller BED (`posbed.bed`, `negbed.bed`, `01_peak_calling/*`, `annotatedpas.bed`) byte-identical to the Stage-2 run; the chr19+21 PBMC slice BEDs byte-identical to the Stage-2 BEDs;
- raw matrices 140.7M counts vs 132.9M shipped (1.059×); the suppressed candidates' partitions plus the tier-2 rows reproduce the shipped per-strand totals to the read (`+` 58,221,832 + 5,131,499 = 63,353,331; `−` 64,246,477 + 5,330,874 = 69,577,351); the +7.8M is the cleavage-window remainder;
- `min_read` keeps 10,339 cells including 1,294/1,294 STARsolo cells (the clip-only run kept 995);
- gene-assigned mass on the shipped run's cells: 1.053× the correctly keyed shipped `filterdmatrix.mtx` (per-gene Spearman 0.991, median per-gene ratio 1.039). The independent verifier returned PROBLEM on the pre-set rule ("within 0.8–1.2× of 41,944,353") because that figure is counts plus the MatrixMarket size line and the shipped `annotated_matrix.mtx` is itself mis-keyed (bug 0a, fixed in `51cbe67`); on the corrected baseline the verifier's own figures are 1.071× (shipped keyed by its pasnumbers), 1.059× (raw) and 1.074× (post-cell-filter). Accepted on that basis — "growth is retention/keying, not over-counting."

**Acceptance run, Stage 1c** — same dataset, code `ebb9fb3`, 12 threads, 54:18 wall, 26.9 GB peak RSS (`results/benchmark_tools/gse104556/peakatail_clipseeded_v3/mouse1`; `CHANGELOG.md`):
- every caller BED identical in columns 1–4 and 6 to the v2 run; only column 5 changes; `annotatedpas.bed` byte-identical (it carries no score column);
- clip reads 2,181,334 → 2,091,279 molecules (0.959×) over 148,715 tier-1 rows; 13.3% (+) / 13.2% (−) of tier-1 rows change value; 3,269 of the 62,159 "≥2 clip read" sites (5.3%) are a single molecule (23.7% on PBMC, where duplication is far higher);
- counting otherwise unchanged to the read (`reads_from_candidates`, `reads_from_window`, `reads_tier2`, `tier1`, `tier2`, `suppressed_candidates` identical); the raw matrices lose exactly the double-counted fallback reads, 140,744,000 → 140,743,255 (−745, −0.0005%); annotated matrix 52,703,320 → 52,703,105 (−215, −0.0004%; 85,816 vs 85,993 PAS — emptied fallback rows fall below `min_cells`); 10,339 cells, 1,294/1,294 STARsolo;
- tier-1 scoring unmoved: P@100 0.4992 / R_det 0.2991 / F1 0.3741, identical to 4 dp before and after (both tiers 0.4093 / 0.3231 / 0.3611);
- chr19+21 PBMC slice: BED columns 1–4 and 6 byte-identical to `002d976` on all four region BEDs; chr21 (+) column 5 sums 14,888 reads → 12,483 molecules. Post-hoc, labelled as such: ≥2 reads n 4,749 / P 0.5519 / R 0.2442 / F1 0.3386 vs ≥2 molecules n 3,664 / 0.6395 / 0.2245 / 0.3323 (≥1 unchanged at n 11,318 / 0.3591 / 0.3379 / 0.3482) — a read threshold and a molecule threshold are not the same operating point.
- `4efeb12` is the verifier's one repair: the YAML alias was accepted silently while the CLI alias warned.

**What is byte-identical vs what is not**, relative to `develop`, at the default gate:
- identical — BED columns 1–4 and 6 (coordinates, `pas_id`, strand), the tier tag (`score == 0` ↔ tier 2), the set of PAS the caller emits, `annotatedpas.bed`, tier-2 matrix rows;
- not identical — BED column 5 (reads → molecules), tier-1 matrix rows (cleavage-site counts), the post-`min_cells` row set of `pasbed.bed` and the annotated matrix, plus the new sidecar and the new `clip_seeded counting:` log line.

## Benchmark result at this commit (`4efeb125`, `manuscript/15_final_gate.md`, verifier verdict SOUND)

Scored with `scripts/benchmark_tools/score_tool.py` against PolyASite 2.0 at 100 bp, strand-matched, detected-gene denominator; every P@100 and R_det reproduced to six decimals by an independent `bedtools` pipeline; all arms ran from a frozen worktree at this HEAD, `LC_ALL=C`.

| output (PBMC 10k v3 unless noted) | n | P@100 | R_det | F1_det |
|---|---:|---:|---:|---:|
| **tier-1 ∩ IP-pass ∩ ≥2 molecules** (pre-registered "precision default") | 44,394 | **0.717** | 0.171 | 0.276 |
| same, testis mouse 1 / mouse 2 | 25,991 / 26,164 | **0.741 / 0.755** | 0.202 / 0.205 | 0.318 / 0.322 |
| tier-1 ≥1 molecule, IP filter (sensitivity arm) | 161,595 | 0.355 | 0.262 | 0.301 |
| tier-1 ≥1 molecule, no IP | 222,955 | 0.303 | 0.297 | 0.300 |
| both tiers, no IP (the tool's literal default output) | 402,765 | 0.193 | 0.336 | 0.246 |
| polyApipe (best other de novo) | 120,916 | 0.380 | 0.199 | 0.261 |
| SCAPTURE | 35,759 | 0.652 | 0.118 | 0.199 |
| PeakATail shipped (pre-#92) | 277,111 | 0.118 | 0.156 | 0.134 |

The pre-registered gate (P@100 ≥ 0.50 on PBMC and both mice) passes on all three call sets; genic-null P@100 is 0.013–0.023. No ≥1-molecule output clears the older 0.38 precision floor. Verifier's wording: the precision default is the most atlas-concordant de novo call set in the panel on both datasets, at polyApipe-or-lower recall (PBMC −14%, mouse −19%) and essentially tied F1 (+0.015 / +0.010) — the claim is precision, not overall accuracy. Note that `≥2 molecules` and the IP filter are run flags (`--polya-min-umis 2 --ip-filter --ip-filter-mode filter`), not defaults changed by this PR.

From `pas_support.tsv`: 72.5% of PBMC tier-1 sites are single-molecule (mouse ≈ 50%); of the ≥2-molecule PBMC sites, 90.2% keep ≥2 under `-F 3844` (mouse 94.8 / 95.6%), so `--polya-clip-filter f3844` would shrink the precision default by ~10% on PBMC.

Compute at this commit: PBMC `--threads 16`, no IP: 3:45:53 wall, **293.7 GB (280.1 GiB) peak RSS**, 143% CPU; IP arm 3:44:37, 239.5 GB; mice (`--threads 12`) 1:01 h, 23.1 / 22.4 GB.

## Known follow-ups (not in this PR)

Two separate PRs are coming; the rest are issues (draft: `manuscript/github/issue10_stage1d_ipfilter_memory.md`):

1. **PR — IP filter minus-strand window.** `ema/experimental/internal_priming.py` applies the forward window `[pos−10, pos+30)` on both strands, so on `−` it tests 30 bp upstream / 10 bp downstream in transcript orientation. A strand-symmetric rule flags 3.5% (PBMC) / 4.2%, 3.6% (mice) of surviving minus-strand sites; benchmark effect negligible (PBMC default P 0.7167 → 0.7182) but the rule is wrong.
2. **PR — memory / CPU.** 280 GiB on PBMC is 2.1× the Stage-2 run (140 GB) and 2.8× the shipped caller (105.6 GB), with ~1.4 cores used despite `--threads 16`; the 1b/1c bookkeeping (`ClipStream`, `ClipAccumulator`, candidate partitions) roughly doubled memory and peak calling is effectively serial. A 10k-cell PBMC BAM currently needs a ≥300 GB node.
3. **Unified BED column 5.** In multi-dataset runs `unified/multi_sample_merged.bed` column 5 is the *first* merged member's molecule count (`bedtools merge -o first`), not cohort support, and a unified PAS can merge a tier-1 member with a tier-2 member. Wanted: cohort-level `clip_umis_sum` / `any_tier1` per unified PAS (we currently recompute from the per-dataset `*.support.tsv` via `multi_sample_pas_mapping.tsv`).
4. **Per-site IP flag in `--ip-filter-mode filter`.** `annotatedpas.bed` carries no IP column in filter mode (only counts in `peak_filters_stats.json`); downstream "trusted novel PAS" needs one in both modes (or an `ip_flags.tsv`).
5. Smaller: `pasbed.bed` silently drops gene-assigned PAS with zero counts after the cell filter (0.45% PBMC, 3.7% mouse1) — document or emit with count 0; `ema/cli/run.py` logs "seq-len not set; defaulting to 150" although 91/98 was resolved.

## Checklist for reviewer

- [ ] On any `clip_seeded` run, `cut -f1-4,6 pasbed.bed` is identical between `develop` and this branch at default flags; only column 5 differs (and `annotatedpas.bed` is byte-identical).
- [ ] `ClipSeeder.flush()` in `ema/countmatrix/polya.py`: the three count sources (candidate partition, cleavage window, clip fallback) are mutually exclusive — each read end lands in at most one tier-1 row. Tests: `test_seeder_window_count_excludes_candidate_intervals_and_neighbours`, `test_clip_fallback_*`, `test_every_tier1_row_has_counts_and_mass_is_not_lost`.
- [ ] The non-positional-strategy detection (shares not summing to the candidate → candidate goes whole to the nearest cluster) is acceptable for `original`.
- [ ] `molecule_cb()` dropping the `RG_` prefix is the intended molecule identity for multi-read-group BAMs; matrix cell identity must stay `RG_` + barcode.
- [ ] `--polya-clip-filter` default `none` is the right call given −8.3% / −27.0% of tier-1 clusters under `f3844`; the sidecar's `*_f3844` columns make the alternative measurable post hoc.
- [ ] Pipeline path: `ClipStream` pickling in the finder → writer hand-off (`peak_pipeline.py`); tile path: sidecar re-keying follows the BED's pasnumber remapping (`tile_runner.py`). `tests/test_polya_three_path_agreement.py` covers both for both strategies.
- [ ] Deprecation: `--polya-min-reads 2` and YAML `polya_min_reads: 2` both warn and both set `polya_min_umis`; `tests/test_polya_clip_evidence_units.py::test_cli_accepts_both_spellings_and_yaml_alias`.
- [ ] `docs/cli/run.md`, `docs/strategies/peak-calling.md` (step 6 "Counting", hyperparameter table) and `CHANGELOG.md` read correctly; the CHANGELOG's two "Measured on GSE104556 mouse1" blocks are the 1b and 1c acceptance runs respectively.
- [ ] `pytest` locally: expect only the 2 known `tests/test_pyproject_install.py` environment failures.
- [ ] Acknowledge the memory footprint as a documented limitation until the follow-up PR lands.
