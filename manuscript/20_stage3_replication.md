# Stage 3 — Laughney patient-level replication (v2-code record, verified SOUND, 2026-08-22)

**Verifier verdict: SOUND.** This is the Stage-3 result of record. Code = frozen `9dfdefb3` (merge of
#93 + the #96 minus-strand internal-priming fix + #97), snapshot
`WD/tools/pa-polya-run-9dfdefb3` (clean worktree, all 17 per-library manifests and every `DONE.ok`
carry the same sha, exit 0). The chain repeats, on the fixed caller, exactly what the v1-code chain
did: same 17-library cohort in one unified PAS space, same pre-registered label-confirmation policy
(labels reused byte-identical — they are GEX-derived and code-independent), same `tier1_ge2` universe
rule, same Fisher cells-mode test with marker pre-selection off, 10 label-shuffle nulls per library,
patient as the unit of replication, MetBone excluded from the pre-registered primary. Universe,
replication and null control were each re-derived by the verifier from raw inputs along an
independent code path.

Pre-registration provenance (md5s recorded here, not only in the run manifests):
`scripts/stage3/UNIVERSE_POLICY.md` `1ac5f0a8cedb48021785f55615587fe8`,
`scripts/stage3/LABEL_POLICY.md` `ff787650a395f2bec37cecae528e931a`,
`scripts/reliability/replication_filter.py` (v0.2.0) `62d121ced8dbf241b7a4c875d3ed51d7`,
`labels/confirmed_labels.tsv` `89bc895007b9cb90c2ee236416b75065` (byte-identical to the v1 tree).

## Verified primary result (MetBone excluded: 15 GSMs = 12 patients)

Design A cohort run (one unified PAS space, 505,197 sites); universe = IP-pass ∧ tier-1 ∧ cohort
clip-molecule sum ≥2 → **80,464 PAS** (40,811 +, 39,653 −); 29,063 curated cells → **18,651 confirmed
(64.2%)** by the pre-registered label policy, 17,538 analysed across the 16 pair-bearing libraries
(17,373 in the 15 primary libraries); 59 cell-type pairs; test = Fisher cells-mode, no marker
pre-selection.

| config | (pair, feature) rows | replicated q<0.05, ≥2 patients, discordance veto | + \|Δprop\| ≥ 0.1 |
|---|---:|---:|---:|
| **PAS, K=2 (pre-registered)** | 2,128,711 | **15,942** (0.749%) | **15,212** |
| gene, K=2 | 482,196 | 10,415 | 10,104 |
| PAS, K=3 (sensitivity) | 2,128,711 | 5,438 | 5,242 |
| gene, K=3 | 482,196 | 3,855 | 3,770 |

The 15,942 replicated (pair, PAS) rows cover **5,951 distinct PAS in 2,883 genes across 47 pairs**
(the 15,212 that also clear the effect floor: 5,645 PAS / 2,841 genes / 47 pairs; K=3: 2,719 PAS /
1,596 genes / 35 pairs). Largest pairs: Epithelial_Tumor–T cell 3,031; Macrophage–T cell 2,774;
B cell–Epithelial_Tumor 1,555. Sensitivity including MetBone (16 GSMs = 13 patients):
16,010 / 15,282 over 2,129,233 rows — a 0.42% difference, and the verifier confirmed the whole
+68-PAS difference is confined to the single pair MetBone contributes
(Epithelial_Tumor_vs_Macrophage).

**Gene-level wording (verifier correction).** "10,415 gene-level switches replicate" would overstate
it: by construction of `collapse_to_gene` this is the count of **(cell-type pair, gene) combinations
containing at least one replicated PAS**, and it was verified to be *exactly* the (pair, gene) set of
the replicated PAS — not an independent gene-level test. Write it as **"10,415 (cell-type pair, gene)
combinations, covering 2,883 genes in 47 pairs."**

## What the #96 fix changed — v1-code comparison row (NOT the record)

The v1-code chain (frozen `4efeb125`, pre-Stage-1d, tree `stage3_laughney_v2/`) is retained only as
the labelled comparison below. Both rows are the same pre-registered primary (15 libraries /
12 patients), so they are like-for-like; the verifier re-parsed the v1 tree with the same parser to
confirm that.

| quantity | v1 code `4efeb125` | v2 code `9dfdefb` | change |
|---|---:|---:|---:|
| universe (PAS) | 74,954 | **80,464** | **+5,510 (+7.4%)** |
| tested (pair, PAS) rows | 1,977,134 | 2,128,711 | +151,577 (+7.7%) |
| distinct PAS tested | 49,874 | 53,531 | +3,657 |
| **replicated PAS, K=2** | **14,480** | **15,942** | **+1,462 (+10.1%)** |
| + \|Δprop\| ≥ 0.1 | 13,826 | 15,212 | +1,386 (+10.0%) |
| (pair, gene) combinations, K=2 | 9,579 | 10,415 | +836 (+8.7%) |
| PAS, K=3 | 4,937 | 5,438 | +501 (+10.1%) |
| distinct replicated PAS / genes / pairs | 5,395 / 2,687 / 47 | 5,951 / 2,883 / 47 | +556 / +196 / 0 |
| sensitivity (13 patients), K=2 | 14,541 | 16,010 | +1,469 (+10.1%) |
| replicated in any of 10 nulls | 0 | 0 | — |

**Where the universe growth comes from (corrected wording — the chain's own note was half the
story).** The change is entirely minus-strand and traces to #96: `posbed.bed`,
`multi_sample_merged.bed` and `multi_sample_pas_mapping.tsv` are byte-identical between the runs, and
so are `any_tier1` and `clip_umis_sum ≥ 2`; only `negbed.bed` changed. The corrected minus-strand
internal-priming window **admitted 20,874 sites and newly discarded 6,403** (negbed 176,952 →
191,423; `n_ip_flagged` 140,822 → 126,351, i.e. +14,471 ip_pass). At universe level that is
**net +5,510 = 6,391 PAS admitted and 881 newly flagged as internal priming**, all minus strand; the
plus-strand universe is unchanged at 40,811. Do **not** write "+7.4% because the fix stopped
discarding legitimate minus-strand sites" without the 881 losses.

**Where the extra replicated switches come from (verifier attribution — do not attribute the delta to
the universe alone).** The v2 replicated set is *not* a superset of v1: 14,123 (pair, PAS) are shared,
1,819 gained, 357 lost (net +1,462). Only **825 of the gains and 81 of the losses** involve PAS whose
universe membership changed (net +744). The other net +718 (994 gains, 276 losses) are PAS present in **both** universes whose
call moved because (i) the per-gene Fisher denominator changed when a sibling minus-strand PAS entered
or left the gene, and (ii) the per-pair BH family changed size — including 64 gained and 6 lost rows on
**plus-strand** PAS whose p-values are bit-identical between the two runs and which moved purely
through BH re-ranking. Of the 4,285 newly tested PAS, 3,974 are new universe admissions and 311 were
already in the v1 universe but only became testable once a restored sibling PAS made their gene
multi-PAS.

## Recommended claim wording (verifier's, adopted; numbers updated to v2)

> Across 12 patients (15 scRNA-seq libraries) of the Laughney lung-adenocarcinoma cohort, processed in
> a single cohort run so that all samples share one PAS identifier space, PeakATail called cell-type
> APA switches in 59 cell-type pairs over a pre-registered precision-first PAS universe (clip-supported,
> internal-priming-filtered, ≥2 clip molecules; 80,464 sites). Requiring a switch to be called (Fisher
> on cells, BH q<0.05) in the same direction in at least two independent patients, with any
> opposite-direction patient vetoing the call, 15,942 of 2,128,711 tested (pair, PAS) hypotheses
> replicate (0.75%), of which 15,212 also exceed an effect floor of |Δproportion| ≥ 0.1; these are
> 5,951 distinct PAS in 2,883 genes across 47 pairs, and 10,415 (cell-type pair, gene) combinations.
> Under ten patient-wise label-shuffle nulls run through the identical pipeline, no feature replicated
> in any combination. Requiring three patients retains 5,438 PAS-level switches.

## Null control — "none in 10 nulls", never an FDR

**0 replicated features in every one of 10 label-shuffle null combinations, in all four configs of
both arms (primary and sensitivity).** 58,859,840 null tests over the 15 primary libraries
(59,090,630 over all 16 label-confirmed libraries) produced **11 nominal q<0.05 hits**; the verifier
scanned all 2,100 null differential files exhaustively and found **no (pair, PAS) hit twice** — not
within a combination, not across combinations, not even by two libraries of the same patient — so
K≥2 replication under the null is 0 by construction, not by luck. MetBone contributes 0 null hits, so
the primary and sensitivity null totals are identical.

Ten index-aligned combinations resolve only to an empirical p = 1/11 = 0.0909 (recorded as
`empirical_p_q` in every per-pair row). **Report as "none in 10 nulls (empirical p ≤ 0.091, the
10-permutation floor)" — never as an FDR estimate.**

Calibration is conservative in every library: pooled null p<0.05 0.84–2.73% (rule ≤7%), p<0.01
0.09–0.47% (rule ≤1.5%), null q<0.05 rate ≤ 7×10⁻⁷ (rule ≤5%; max 6.8×10⁻⁷ in GSM3516663 on 3 hits in
4,389,688 tests — note the verifier's headline "≤3.3e-7" is the value for GSM3516672, not the cohort
maximum, which is recomputed here from the chain's own `null_summary` blocks).
**GSM3516672-StageIB again trips the pre-registered per-library count rule** (3 of 10 permutations
carry a q<0.05 hit; rule ≤25%, i.e. ≤2), on 3 hits in 9,096,302 null tests — BH discreteness in the
biggest library (28 pairs), not anti-conservativeness, since all four rate criteria pass comfortably.
Kept and disclosed exactly as in v1; the disclosure must travel with the number, because
"all libraries pass calibration" would be false.

Two-GSM patients verified to count once: for Macrophage_vs_T_cell the library-as-unit counterfactual
gives 2,863 vs 2,774 patient-as-unit — collapsing removes 89 PAS that would have "replicated" on two
libraries of a single patient (LX675 144, LX682 588, LX684 477 same-direction double calls). The
discordance veto also bites (18 and 10 PAS vetoed in the two audited pairs, including 2 where LX682's
two libraries call opposite directions).

## Disclosure issues (verified on v2; none change the counts)

1. **Top-gene names are unreliable in overlapping loci — the v1 problem persists, unchanged.** On v2,
   **12 of the top-30 gene rows (40%)** have a representative PAS
   that is *not* in the named gene's 3′ UTR, and **8 of 30 (27%)** sit in a *different* gene's 3′ UTR:
   AC073610.2→FKBP11, SENP3-EIF4A1→CD68 (×2), DGKZ→MDK, RHOA→GPX1 (×2), AC004922.1→ARPC1A,
   CD163L1→CD163. (**Two counts corrected here, 2026-08-22, by re-running the verifier's own
   `verify_v2chain/topgene_check.py`.** (i) The verifier record and the first draft of this note said
   "7 of 30" on the narrow reading, but the enumerated list above already holds **eight** rows; the
   script prints eight `other-gene 3' UTR` rows, so the narrow figure is 8/30 = 27%. (ii) Running the
   same script against the v1 tree gives the **same 12 of 30 (40%)** on the strict reading and 9 of 30
   (30%) on the narrow one, so the defect is **unchanged between v1 and v2**, not "slightly worse", and
   the "33% / 10 of 30" carried in the v1 note is superseded by that recount. **State which definition
   you use**, and quote 40% strict / 27% narrow for v2.) CD68, PTPRCAP, FKBP11, GPX1, MDK, STARD10 and ARPC1A still have
   **zero** assigned PAS cohort-wide, unchanged from v1; even "correct" rows are ambiguous (the CORO1B
   PAS also lies in PTPRCAP's 3′ UTR). Ambient-IG signature (IGLL5 in non-plasma pairs) was noted on v1 and was **not** re-checked on v2.
   Replication statistics are computed on PAS identifiers and are unaffected — but **no named
   top-gene list may be published** until PAS→gene re-assignment (tool issue #99).
   Genomic context of the 5,951 replicated PAS, recomputed against GRCh38.99 for Fig 6:
   **96.8% in a same-strand gene body, 61.5% in *any* same-strand 3′ UTR, 57.8% in the *assigned*
   gene's own 3′ UTR.** **Label correction:** the v1 note's "62.5% of replicated PAS in own-gene
   3′ UTRs" was mislabelled — recomputing v1 with an explicit definition gives 62.6% for any
   same-strand 3′ UTR and only 58.6% for the assigned gene's own 3′ UTR. Quote both figures, labelled.
2. **The MetBone exclusion stands (pre-registered) but its stated premise was wrong** — reproduced on
   v2 data. MetBone carries 306,202 clip molecules, indistinguishable from GSM3516662-StageIA's
   306,068 which the caller did *not* flag, and 8,753 of the 28,983 PAS it retains after the
   universe/tier filter (**30.2%**) carry ≥2 of its *own* clip molecules. (**Denominator label
   corrected 2026-08-22:** the verifier called 28,983 "its tested PAS"; 28,983 is
   `n_pas_after_tier_filter`, and the number MetBone actually contributes a test for is **23,079** —
   the rest are the only PAS in their gene, so no Fisher table exists. On that stricter denominator
   6,392 of 23,079 (27.7%) carry ≥2 of MetBone's own clip molecules. Either denominator refutes
   "no clip evidence".)
   The 0.015% clip rate came from the caller sampling only the first 200,000 CB
   reads of a coordinate-sorted BAM (head of chr1); the same head-sampling warning fired on **five**
   libraries (GSM3516663, GSM3516664, GSM3516671, GSM3516672, GSM3516673), three of which are among
   the highest-yield in the cohort. **Never repeat the "no clip evidence" / "0.015% clip rate"
   justification as fact.** Primary and sensitivity agree to 0.42%. Tool issue #99 §2 (the
   head-sampling clip-rate estimator; #99 §1 is the PAS→gene assignment defect of disclosure 1).
3. **Reporting hazard in the chain's own outputs.** `per_pair_*.tsv` and `REPORT.md` carry a column
   `exp_false_frac_q = 0.0` / `expected_false_replicated_fraction 0.0`. Read literally that asserts
   FDR = 0, which the pre-registration forbids (10-permutation resolution floor). **This column must
   be dropped or relabelled before any of these tables is used as a manuscript source.**
   Cosmetic provenance defects, unchanged from v1 and harmless to every number: each per-library
   switch `run_manifest.json` records a per-sample `ema run` flag string ending `--threads 8` for a
   per-sample run that never happens in cohort mode (the cohort caller used `--threads 16
   --peak-workers 16`), and the `low_clip_rate_warning` block records `clip_rate_pct 0.015` with no
   in-file caveat that the estimator samples only the first 200,000 CB reads.

Files: `/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/replication/{primary_noMetBone,sensitivity_all}/`,
`REPORT.md`, `PRIMARY_vs_SENSITIVITY.txt`; verifier artefacts (independent universe rebuild,
independent replication re-derivation, exhaustive null scan `null_hits.tsv`, v1↔v2 set comparison) in
`.../stage3_laughney_v3/verify_v2chain/`. Figure: [Fig 6 `laughney_switches`](figures/laughney_switches.caption.md),
regenerated from this tree with `LAUGHNEY_SWITCHES_VERSION=v2_code` (now the script default; `v1_code`
still reproduces the comparison row).
The v1-code tree `/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/` is retained read-only.

## Naming convention (use these exact phrases in the paper)

- **The Laughney cohort** = **17 libraries from 14 patients** (three patients contributed a tumour and a
  normal library). This is the cohort that was peak-called in one run and is what "17 samples/datasets"
  in [13](13_reliability_positioning.md) §2 and [19](19_final_gate_v2.md) §5 refers to.
- **The replication primary** = **15 libraries from 12 patients** — the cohort minus the pre-registered
  MetBone exclusion (one patient, one library) and minus GSM3516671, which retains only one confirmed
  cell type and therefore yields no testable pair.
- Null-test totals are scope-dependent: **59.1 M** across all 16 label-confirmed libraries, **58.9 M**
  across the 15 libraries of the primary. Quote the primary figure with the primary result.
- The labels are code-independent and were reused unchanged: `n_cells_h5ad` rose by only 13 cells
  cohort-wide (7 libraries, extra cells admitted by the min-PAS-per-cell filter now that more
  minus-strand PAS exist), while matched confirmed cells (18,651) and analysed cells (17,538) are
  **identical to v1 in every library**.
