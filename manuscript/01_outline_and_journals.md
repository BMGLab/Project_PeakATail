# PeakATail Manuscript Skeleton — Q1 Methods + Application Paper

**Working draft v0.2 (2026-08-21) — re-positioned per `13_reliability_positioning.md` (PI decision: reliability first).** Grounded in PeakATail `develop` post-PR #92 (clip-seeded calling merged, verified SOUND on the chr19+21 slice; Stage 1c UMI-dedup clip evidence in flight; repo `tools/PeakATail`, docs bmglab.github.io/PeakATail). The clustering-novelty claim is **dropped**: the ablation showed PAS-count clustering is 3'-end expression signal re-encoded (see Supplementary S4). The paper now claims four things — **trustworthy PAS detection (low false positives), reliable cell-type-specific APA switches, 3'UTR length change, and trusted de novo PAS** — each tied to a gate pre-registered in `13_reliability_positioning.md` before the final Stage-2/3 run exists.

**Number policy for this file.** Only adversarially verified numbers appear as results (sources: `07_curated_benchmark_report.md`, `09_headtohead_results.md`, `05_figure_index.md` §1–8, `12_stage2_gate.md` Stage-1b verdict). Anything that depends on the final clip-seeded + UMI-dedup run, on re-keyed count matrices (bug 0a), or on the quarantined spermatogenesis/Kinnex figures (`11_verifier_corrections.md`) is written as **[final run pending]**. The retired 0.9986 precision and the stale Stage-2 gate numbers are never quoted as results; the shipped-caller head-to-head numbers (F1 0.134 PBMC) are quoted only as the *baseline* the fix was measured against.

---

## (a) Ranked candidate journals

The ranking is unchanged; the fit argument is re-justified for a reliability/benchmark-centred methods paper whose headline is *how much of a single-cell PAS call set can be trusted, and how to make it trustworthy* — not a biological discovery and not a clustering novelty.

| Rank | Journal | Fit rationale (reliability/benchmark framing) | Format constraints (methods+application paper) |
|---|---|---|---|
| 1 | **Genome Biology** (IF ~12) | Still the best fit, and the fit has *improved*: the Method article type explicitly rewards rigorous benchmarking with nulls and honest comparison, and GB published APAeval — the community's own APA benchmark — so a paper whose core is "pre-registered gates, gene-body-shuffled nulls, long-read truth, replicate reproducibility, and a disclosed negative result that motivated a caller redesign" speaks their language. Sierra and scDaPars precedents remain. Reviewers there will accept a tool that is not the F1 leader if its precision-first tier, its calibration story and its provenance are demonstrably sound. | Method article: no hard word limit, ~6–8 main figures; structured abstract (Background/Results/Conclusions, ≤350 words — our ≤200 w draft fits); open access APC (~$4k); code in a DOI-archived repository (Zenodo) under an OSI license (MIT ✓); mandatory "Availability of data and materials". |
| 2 | **Nature Communications** (IF ~15) | Fit has *weakened*: without a headline biological discovery the paper is a methods/benchmark paper, which NatComms takes only when the benchmark itself changes field practice. Our strongest NatComms argument is now the field-level finding (evidence type, not peak-calling sophistication, discriminates callers; atlas-only benchmarking flatters internal-priming artefacts) plus the Laughney per-cell-type switch application **[pending Amir's rerun]**. Attempt only if Fig 6 yields a replicated, orthogonally supported program. | Article: ~5,000 words (excl. Methods), ≤10 display items, unstructured abstract ≤150 words (trim the draft), Methods unlimited; Reporting Summary + Code Availability mandatory; APC ~$6.8k. |
| 3 | **Genome Research** (IF ~6–7) | Methods/Resource track; excellent fit for accuracy-centred work with nulls and long-read truth (PolyA-seq/Derti lineage is a GR paper). Natural second submission after GB with no restructuring. | ~55,000 characters incl. references; abstract ≤250 words; typically ≤8 figures; software free to academics, data public before review. |
| 4 | **Nucleic Acids Research** (IF ~14–16) | PolyASite and PolyA_DB are NAR papers — a benchmark built on curated PolyASite 2.0 point sites, with the circularity caveats handled explicitly and a trusted-novel definition validated *independently* of the atlas, will be read by exactly the right reviewers. Risk unchanged: NAR methods papers often want a resource/webserver flavour. | No rigid length limit (~9 typeset pages); abstract ≤200 words + graphical abstract; open access mandatory (APC ~$5k). Fallback: **NAR Genomics & Bioinformatics**. |
| 5 | **Briefings in Bioinformatics** (IF ~7–9) | Fit has *improved*: BiB's centre of gravity is comprehensive tool comparison, and Fig 3 (six tools × two datasets × two truth sets, one scorer, three-seed nulls, reproducibility metric) is now one of the paper's two strongest assets. If the final-run gates are cleared only narrowly, this is the safest high-visibility home. | ~6,000–7,000 words; abstract ≤250 words; "Key Points" box (3–5 bullets); figures ~6–8; code availability mandatory. |
| 6 | **Bioinformatics** (OUP) (IF ~5–6) | Reliable venue (SCAPE, scDAPA). Length limit forces the application and half the benchmark into supplements; viable if we descope to "precision-first caller + calibration". | Original Paper ~7 pages (~5,000 words all-inclusive); structured abstract (Motivation/Results/Availability); reviewers install software. |
| 7 | **PLoS Computational Biology** (IF ~4) | Methods section; no length pressure; its open-science ethos fits the pre-registration, full disclosure of the failed gate and the FDR-failure-and-fix story. Safety net. | No length limit; abstract ≤300 words + lay author summary; unlimited figures; archived open-source software; APC ~$3k. |

**Recommended strategy:** write to Genome Biology's Method format (structured abstract, 6 main figures, pre-registered gates cited in Results). It degrades to Genome Research and BiB without restructuring. Escalate to Nature Communications only if Fig 6 delivers a replicated per-cell-type APA program with orthogonal support.

---

## (b) Title candidates

1. **PeakATail: trusted poly(A)-site detection and reliable cell-type APA switches from standard 3' scRNA-seq** *(safe, descriptive; leads with the two primary claims)*
2. **Precision-first polyadenylation-site calling in single cells: poly(A)-tail read evidence, calibrated testing and replication filters in PeakATail** *(mechanistic; names the three components that buy reliability)*
3. **How many single-cell poly(A) sites can be trusted? A null-referenced, long-read-validated benchmark and a precision-first caller** *(question title; benchmark-centred — best for BiB/GB if the gates clear narrowly)*
4. **Low-false-positive poly(A)-site discovery and replicated 3'UTR switches in single cells with PeakATail** *(compact; fits NatComms/NAR title length)*
5. **Reliable alternative polyadenylation analysis from 3' scRNA-seq: PeakATail's clip-seeded calling benchmarked against curated atlases and long-read truth** *(explicit about the truth sets; longest)*

Retired: every earlier candidate that led with "clustering", "cell-state discovery from PAS usage alone", "expression-independent", or "beyond gene counts".

---

## (c) Abstract draft (≤200 words; bracketed = final run pending)

> Alternative polyadenylation (APA) remodels 3'UTRs across cell types, and the poly(A)-site (PAS) signal is present in every 3' scRNA-seq dataset — but single-cell PAS calls are rarely benchmarked against independent truth, and we show that on deep 10x data coverage-shape callers, including our own earlier version, emit call sets in which most sites lack atlas or long-read support. We present PeakATail, an aligner-agnostic framework that seeds PAS from non-templated poly(A) soft-clipped reads, reports clip-supported and coverage-only sites as explicit tiers with per-site molecule counts, filters internal priming, and tests cell-type APA switches with calibrated, replication-filtered statistics. Against curated PolyASite 2.0 point sites with gene-body-shuffled nulls and against [104 M] PacBio Kinnex poly(A)-verified molecules, PeakATail's default tier reaches precision [0.XX] on PBMC and [0.XX/0.XX] on two testis replicates (pre-registered gate ≥ 0.50), at recall [0.XX]; a six-tool benchmark shows precision, recall and replicate reproducibility trade off and that evidence type, not peak-calling sophistication, separates callers. PeakATail recovers the spermatocyte→round→elongating spermatid 3'UTR shortening gradient [ρ], yields atlas-novel PAS of which [XX%] are long-read-confirmed, and maps replicated per-cell-type switches across [17] lung-cancer samples. Open source (MIT): github.com/BMGLab/PeakATail.

*(≈195 words. Trim to ≤150 for NatComms by cutting the benchmark clause. Every bracket is filled from the final Stage-2/3 run only — never from the stale run.)*

---

## (d) Figure-by-figure plan

### Fig 1 — Tool overview: clip-seeded, tiered PAS calling and the reliability pipeline
- **A.** Pipeline flow: tagged BAM(s) → per-strand read ingestion with UMI dedup and SAM-flag filtering (−F 3844 semantics) → **clip-seeded PAS calling** (poly(A) soft-clip clusters seed candidates; coverage summits without clip support are emitted as a second tier) → internal-priming filter (genome FASTA) → PAS×cell matrix → tiered gene-end/3'UTR annotation → `switch diff / length / trend / geneview`. (Base: `pipeline_flow_report.png`, to be redrawn.)
- **B.** Evidence schematic: why coverage shape alone is insufficient on deep data (clip reads sit ~90–105 nt past the coverage peak end — the cleavage-offset finding, verified in `07_curated_benchmark_report.md` §3) and how clip seeding recovers evidence that falls outside every coverage peak window.
- **C.** The two output tiers and the default: tier-1 = clip-supported with ≥2 distinct molecules and IP-pass (default output); tier-2 = coverage-only, written to a separate `low_confidence` file; the ≥1-molecule output as the sensitivity arm. State here that the default was pre-registered in `13_reliability_positioning.md` before the final run.
- **D.** Worked gene example: one 3'UTR, proximal/distal PAS, two cell types, with clip-support counts per site (real `switch geneview` track, regenerated from the final run).

### Fig 2 — Trustworthy PAS detection: accuracy vs curated atlas and long-read truth, with nulls, two tiers
- **A.** Precision@100 bp vs curated PolyASite 2.0 representative sites (strand-matched, point mode) for default tier, ≥1-molecule output and tier-2, on pbmc_10k_v3 and both GSE104556 testis mice, with the 3-seed gene-body-shuffled null (≈0.02) drawn on every bar. Gate lines: **P ≥ 0.50** (default output, pre-registered) and **P ≥ 0.38 / F1 > 0.261** (original two-sided gate, reported for the ≥1-molecule output). **[final run pending]**
- **B.** Precision–recall curves across windows (10/25/50/100/200 bp) with recall on the detected-gene denominator (285,136 PBMC / 126,686 testis sites); recall ceilings stated (only ~29% of detected-gene atlas sites carry any clip read in this BAM).
- **C.** Atlas-independent truth: precision vs PacBio Kinnex PBMC long-read 3' ends (≥5/≥20/≥100/≥500 record-support tiers; donor-matched and -mismatched), with the genuine / internal-priming / unsupported decomposition per tier. **[quarantined figure; regenerate per `11_verifier_corrections.md`; report record-count support or dedup first]**
- **D.** Sequence architecture (verified): AATAAA|ATTAAA within +0..+100 nt of the coverage-peak end 37.1% vs 14.7% null; modal hexamer +75 nt; A-fraction crest 0.426 at +98 nt — and the same panel re-anchored on clip-seeded cleavage sites, where the canonical −40..−5 window should now carry the signal **[final run pending]**.
- **E.** Internal-priming filter effect: A-rich downstream tracts 6.3% → 3.0% (null 2.2%) on the shipped arm (verified); fraction of tier-1 calls removed by IP filtering on the final run **[pending]**.
- **F.** Gate table inset: each pre-registered gate, the value, pass/fail — no arm omitted.

### Fig 3 — Head-to-head trade surface and replicate reproducibility
- **A.** Six tools (PeakATail default + sensitivity arm, polyApipe, SCAPTURE, Sierra, scAPAtrap; scUTRquant shown but unranked as annotation-based), two datasets, one scorer (`score_tool.py`), per-dataset detected-gene denominators, 3-seed nulls. Shipped-caller baseline row kept for transparency: PBMC P 0.118 / R 0.156 / F1 0.134 (last among de novo tools); testis F1 0.258/0.260 (second). polyApipe reference PBMC 0.380 / 0.199 / 0.261. Final-run PeakATail rows **[pending]**.
- **B.** Trade surface: precision vs recall with call count as marker size (21k–787k, 37-fold range) — the figure argues that single-number rankings are misleading; per-dataset reporting is mandatory (shipped PeakATail precision spans 0.118 PBMC → 0.369 testis → 0.447 Laughney).
- **C.** Replicate reproducibility@100 bp on the two testis mice: scAPAtrap 0.79–0.88, Sierra 0.79–0.83, PeakATail (shipped) 0.73–0.75, polyApipe 0.49–0.50; depth stratification as mechanism (62% of polyApipe's mouse calls are depth-1 and reproduce 0.31 vs 0.81 for depth ≥2). Honest scope: PeakATail's reproducibility advantage is **against polyApipe only**; scAPAtrap's lead is measured after its own `reducePeaks` depth cleaning. Final-run PeakATail **[pending]**.
- **D.** Long-read re-ranking: Spearman ρ between atlas and Kinnex precision orderings (0.90 at 7 of 8 arms, 0.80 at one — never say "every"); Sierra's atlas precision was ~half internal priming, a caution for atlas-only benchmarking. **[quarantined; regenerate]**
- **E.** Compute: runtime and peak RSS per tool (PeakATail is the heaviest in the panel; state it) and the `reannotate` re-branching cost argument.
- **F.** Fairness checks as a small table: all tools ≥95% gene-proximal (verified), so PeakATail's `max_gene_distance` restriction does not bias the comparison; gene filter is F1-neutral on raw output.

### Fig 4 — Reliable cell-type APA switches: calibration and the replication filter
- **A.** The failure (verified, shipped tests on mis-keyed matrices): label-shuffle null p-value histograms — Fisher 38.6% and NB-pairwise 21.3% of null tests p < 0.05; 20/20 null runs report q < 0.05 hits. This is reported as a finding, not hidden.
- **B.** The fix: re-calibration on correctly keyed count matrices with `--count-mode cells` (de-pseudoreplicated) and NB with shrinkage; null p < 0.05 fraction after fix **[final run pending]**. If still off nominal, the default becomes **permutation-calibrated q-values** (label-shuffle null, ≥20 permutations) — decision rule pre-registered in `13_reliability_positioning.md` §2.
- **C.** Replication filter: a switch is reported only if called in the same direction in ≥2 independent samples (testis 2 mice; Laughney 17 samples); number of switches before/after, and the expected false-positive rate of the filter from between-donor same-cell-type contrasts (the specificity null of `02_validation_plan.md`). **[pending]**
- **D.** Effect-size floor |ΔPDUI| ≥ 0.1 alongside q; volcano with both thresholds; canonical positive controls (IGHM plasma-vs-B, proliferation-linked global shortening) with correct direction **[pending]**.

### Fig 5 — 3'UTR length biology and trusted de novo PAS
- **A–C. Spermatogenesis positive control (GSE104556, both mice).** Per-cell 3'UTR length index across spermatocyte → round spermatid → elongating spermatid: **SPC>RS>ES is the claim** (ρ −0.73 to −0.86 under every labelling tried on the quarantined run); the spermatogonia first step is fragile (inverts under marker-argmax labels, rests on 64/68 cells) and goes to supplement. Per-gene replicate ρ (0.874 on the quarantined run), 6/6 literature genes in the expected direction, 200-permutation null; the intronic/other-exon peaks that flip the sign shown explicitly, and the `wul` vs `wdi` index discrepancy resolved before the atlas arm is used. **[all numbers final run pending; regenerate per `11_verifier_corrections.md`]**
- **D–F. Trusted de novo PAS (pre-registered definition, `13_reliability_positioning.md` §3):** clip-supported ≥2 molecules AND not IP-flagged AND canonical hexamer within −40..−5 nt AND ≥100 bp from any PolyASite 2.0 / PolyA_DB site. Validation: fraction with a Kinnex long-read 3' end within 25 bp — **gate ≥ 70%** against a shuffled-position null; genomic-region and motif-class composition of trusted-novel sites; a geneview vignette of one long-read-confirmed novel site. **[final run pending]**

### Fig 6 — Application: replicated per-cell-type APA switches in lung cancer (Laughney GSE123904, 17 samples)
- **A.** Cell-type annotation (GEX-derived labels; PeakATail does not claim to define cell types) and per-cell-type PAS coverage.
- **B.** Per-cell-type switches across Normal → Stage I → Metastasis with the replication filter (≥2 samples, same direction) and calibrated q; counts per cell type.
- **C.** TAM/myeloid 3'UTR programs as the worked case; direction calls from `switch length`; pathway enrichment of replicated switches.
- **D.** Geneview vignettes for 2–3 headline genes; orthogonal support (bulk 3'-seq or atlas) for at least one.
- **Status: entirely pending Amir's #67 SWITCH_CELLTYPE rerun on re-keyed, counts-layer matrices with the fixed `main.nf` and `--count-mode cells`.** If the rerun slips, Fig 6 folds into Fig 4 as the replication-filter demonstration and the paper ships as methods + benchmark.

### Supplementary figures
- S1. Dataset summaries (cells, depth, chemistry, QC; clip rate per BAM — CellRanger vs STARsolo, since any pipeline that trims poly(A) before alignment destroys the evidence).
- S2. Peak-call QC: width, per-gene PAS counts, tier composition, molecule-count distributions; cleavage-offset analysis on the shipped caller (verified) vs clip-seeded cleavage placement.
- S3. Strategy comparison in-framework (`original`, `lambda_poisson`, `sierra_iterative`, `lambda_gradient`, `clip_seeded`) on identical regions — the verified result that the old dump-based atlas could not rank them (spread 0.004) and the curated-atlas spread (0.355–0.492) that can.
- S4. **PAS-profile clustering as an optional utility (the dropped claim, stated honestly).** TF-IDF/LSI/Leiden on PAS counts recovers GEX cell types (PBMC AMI 0.708 / ARI 0.502; Laughney median AMI 0.662 / ARI 0.463, 17/17 samples, verified). **Ablation: collapsing the 275,370 sites to 14,891 gene totals gives AMI 0.698 — the site resolution adds nothing, so this is 3'-end expression signal re-encoded, not isoform choice.** The within-gene usage-fraction spaces (the isoform-only test) were computed on mis-keyed matrices and are re-run on re-keyed ones **[pending]**; the claim does not depend on their outcome. Offered as a convenience when no GEX labels exist, nothing more.
- S5. Statistical validation extended: `count_mode="cells"` vs `"reads"`; `nb_multi` `sample_split`; min-cells sensitivity; permutation-calibration runtime.
- S6. **The FDR-failure-and-fix story as a methodological contribution.** (i) Null-calibration design (label-shuffle, 20 perms, resumable); (ii) the three causes found — read-level pseudoreplication, marker double-dip, NB dispersion floor at 1e-4 carrying ~half the false hits; (iii) the independent matrix mis-keying bug (99.87% of annotated PAS carried another PAS's counts on testis mouse1; bisected to 2026-03-23) that voided every count-based result and the `test_matrix_pas_id_row_alignment` regression test; (iv) the strand inversion in `switch length` (100% of minus-strand genes) and its CI guard; (v) calibration after the fix. Framed as "what it takes to make single-cell APA statistics trustworthy", with a checklist other tools can apply.
- S7. Internal-priming characterization: genomic A-tract definitions, flag rate per dataset (3.8% on test data, verified live and off by default in v0.2.0; final-run rates pending), long-read-measured IP fraction per tool **[regenerate]**.
- S8. Benchmark reproducibility: environment archaeology and bug fixes every compared tool needed; `score_tool.py` contract; null construction; truth-set size dependence of F1 (why F1 is quoted only with a named truth set).
- S9. Gate history: the original two-sided gate, the stale Stage-2 run's failure, the post-hoc threshold sweep examined and disclosed (numbers not reported as results), the pre-registered default that followed, and the final-run outcome. Full disclosure that the default was set after seeing the stale run but before the final run.
- S10. Compute scaling detail; `collapse`/`combine` multi-library workflows; `reannotate` branching.
- S11. Spermatogonia step and the `wul`/`wdi` index reconciliation (from Fig 5 corrections).

### Supplementary tables
- T1. Full 19-tool feature comparison (from `combined_polya_scrna_methods.csv`), with the evidence-type column (coverage / poly(A) clip / sequence model / catalog) added.
- T2. Per-tool, per-dataset, per-window precision/recall/F1 with null values (`benchmark_consolidated.tsv`).
- T3. All pre-registered gates with values and pass/fail.
- T4. Replicated differential APA results per contrast (findings_long export), with per-sample direction calls.
- T5. Trusted-novel PAS list with long-read support.
- T6. PDUI/length scores per gene × cell type (length_long export).
- T7. Software versions, parameters, frozen commit hashes, `run_config.json` / `run_manifest.json` digests for every run (runs must freeze their code tree — a lesson from Stage 2).

---

## (e) Section-by-section outline

### Introduction (~5 paragraphs)
1. APA biology: 3'UTR isoform choice controls stability/localization/translation; remodeled in proliferation, cancer, immune activation.
2. 3' scRNA-seq reads pile near poly(A) sites — the signal is in every 10x dataset, but deep data also produces abundant coverage peaks that are *not* PAS, and single-cell PAS call sets are rarely validated against independent truth.
3. Existing single-cell APA tools (Sierra, scAPAtrap, SCAPTURE, polyApipe, scUTRquant, SCAPE, scDaPars, scTail, …) differ mainly in evidence type (coverage shape / poly(A) soft-clips / sequence models / catalogs); published benchmarks rarely use nulls, per-dataset denominators or long-read truth; differential tests are typically not calibrated and not replicated.
4. Gap → PeakATail: clip-seeded, tiered, IP-filtered PAS calling with per-site molecule support; calibrated and replication-filtered switch testing; an atlas-independent definition of trusted novel PAS; all within one aligner-agnostic, provenance-tracked framework with pre-registered gates.
5. Contributions summary: (i) precision-first PAS detection validated on two truth sets with nulls; (ii) a six-tool benchmark that reframes the field's accuracy claims; (iii) calibrated, replicated cell-type switches; (iv) UTR-length biology positive control; (v) trusted de novo PAS; (vi) open source with full disclosure of the negative result that drove the redesign.

### Results

**R1. PeakATail: clip-seeded, tiered PAS calling in an end-to-end reliability pipeline** *(descriptive; Fig 1)*
- **Key claim:** PeakATail turns any CB/UB-tagged BAM into tiered, evidence-annotated PAS calls and calibrated single-cell APA tests through one modular pipeline.
- **Exact evidence:** architecture description; the clip-seeding algorithm (25-bp single-linkage clusters of poly(A) soft-clip sites, read-weighted mode as cleavage position, UMI-deduplicated support, coverage summits without clips as tier-2); IP filter; UMI dedup / flag filtering; `reannotate` branching; provenance (ledgers, content-hashed manifests, frozen commit). Verified engineering evidence: on the chr19+21 slice the clip-seeded caller lifted F1@100 from 0.177 (shipped) to 0.291 both tiers / 0.348 tier-1, independently re-scored (Stage 1, verified SOUND). No gate applies here; no full-run number is quoted.

**R2. PeakATail's default output is a low-false-positive PAS set against curated-atlas and long-read truth** *(Fig 2)*
- **Key claim:** the pre-registered default (clip-supported, ≥2 distinct molecules, IP-pass) achieves **P@100 ≥ 0.50 on pbmc_10k_v3 and P ≥ 0.50 on both testis mice** against curated PolyASite 2.0 point sites (gate, `13_reliability_positioning.md` §1), ≥20× the gene-body-shuffled null, with recall reported and not gated; the ≥1-molecule output is reported against the original two-sided gate (F1 > 0.261 AND P ≥ 0.38). **[final run pending]**
- **Exact evidence needed:** `score_tool.py` on the final-run `pas.bed` for default tier, ≥1-molecule output, tier-2, per dataset, windows 10–200 bp, 3-seed nulls, detected-gene denominators; Kinnex long-read precision per support tier with genuine/IP/unsupported decomposition (record-count support stated or deduplicated); motif architecture re-anchored at clip-seeded cleavage sites; IP-filter delta. Verified context that stays: curated-atlas precision of the shipped arms 0.36–0.49 (~20× null at 100 bp), the ~90–105 nt cleavage-offset, IP filter 6.3% → 3.0%. **Honesty rules:** always report both tiers with n; never quote tier-1 precision without its recall; F1 only with the truth set named; recall ceilings stated.

**R3. A six-tool benchmark shows a trade surface, not a ranking — and that evidence type discriminates** *(Fig 3)*
- **Key claim:** across six tools and two datasets under one scorer, precision, recall, call count and replicate reproducibility trade off; PeakATail's final default sits on the high-precision side of that surface and its ≥1-molecule arm on the sensitivity side **[final run pending]**; the shipped coverage-only caller ranked last on PBMC (F1 0.134 vs polyApipe 0.261), and the gap was evidence type, not threshold or offset (verified).
- **Exact evidence needed:** final-run PeakATail rows added to `benchmark_consolidated.tsv`; replicate concordance@100 on testis; depth stratification; long-read re-ranking (ρ 0.90 at 7/8 arms); fairness table (≥95% gene-proximal for all tools, verified); resource logs. **Honesty rules:** report where competitors win (SCAPTURE precision 0.652/0.694 among de novo tools; scAPAtrap and Sierra reproducibility above PeakATail's shipped 0.73–0.75); state that PeakATail's reproducibility advantage is against polyApipe only; state PeakATail is the heaviest tool in compute.

**R4. Cell-type-specific APA switches are reported only when calibrated and replicated** *(Fig 4)*
- **Key claim:** PeakATail's shipped tests were anti-conservative (Fisher 38.6% / NB 21.3% null p < 0.05, verified); after re-keying and de-pseudoreplication the null p < 0.05 fraction is ≈ 5% **[pending]** — otherwise permutation-calibrated q-values become the default — and a switch is reported only if replicated in ≥2 samples with |ΔPDUI| ≥ 0.1 (`13_reliability_positioning.md` §2).
- **Exact evidence needed:** null p-value histograms before/after on re-keyed counts-layer matrices; permutation-q pipeline; replication-filter yields and between-donor specificity null; canonical positive controls (IGHM plasma-vs-B intronic switch; proliferation-linked shortening) recovered with correct direction at calibrated q. **Blocking dependency:** counts-layer backfill + `--count-mode cells` CLI exposure (#86).

**R5. PeakATail recovers the spermatogenesis 3'UTR shortening gradient and yields long-read-confirmed novel PAS** *(Fig 5)*
- **Key claim (length):** 3'UTR length decreases SPC > RS > ES in both mice under every labelling, with per-gene replicate agreement and literature genes in the expected direction; the SPG step is reported as fragile in supplement. **[final run pending; preliminary ρ −0.73..−0.86 quarantined]**
- **Key claim (trusted novel):** ≥70% of trusted-novel PBMC PAS (pre-registered definition) have a Kinnex 3' end within 25 bp, against a shuffled-position null (gate, `13_reliability_positioning.md` §3). **[final run pending]**
- **Exact evidence needed:** regenerated spermatogenesis panels on re-keyed matrices with both "by peak count" and "by UMI signal" framings and the `wul`/`wdi` discrepancy resolved; trusted-novel set, Kinnex support fraction with null, composition and vignette.

**R6. Application: replicated per-cell-type APA switches across 17 lung-cancer samples** *(Fig 6; pending Amir's rerun)*
- **Key claim:** per cell type, a set of APA switches across Normal → Stage I → Metastasis that replicate across samples under the calibrated test and replication filter; a TAM/myeloid 3'UTR program as the worked example with orthogonal support for ≥1 gene.
- **Exact evidence needed:** #67 SWITCH_CELLTYPE rerun on re-keyed, counts-layer matrices with the fixed `main.nf` (unified merged bed), `--count-mode cells`, replication statistics, effect-size correlation across samples. **Fallback:** fold into Fig 4 if the rerun slips; the methods + benchmark paper stands without it.

### Discussion
- What "trusted" means operationally and why precision-first with explicit tiers is the right default for a field where most published call sets have never been scored against a null.
- The benchmark lessons: evidence type over peak-calling sophistication; per-dataset reporting; atlas-only validation flatters internal-priming artefacts; F1 is truth-set-size dependent.
- The calibration lessons: pseudoreplication, double-dipping, dispersion floors, and the value of permutation nulls and replication filters.
- Interpretation of spermatogenesis and Laughney findings against the APA literature (scoped to what replicated).
- **Limitations (state plainly):** the clip-evidence recall ceiling (~29% of detected-gene atlas sites carry any clip read in the PBMC BAM; chemistry/aligner dependent — poly(A)-trimmed pipelines lose the evidence entirely); PeakATail is the heaviest tool in compute; the precision-first default was pre-registered after a stale run failed the original gate (disclosed); the PAS-profile clustering is an expression proxy, not isoform-choice discovery; Laughney samples are donor-correlated; no Read1 cleavage mode; long-read truth support counted in records, not molecules, in the first version.
- Future: Read1/TVN-priming cleavage mode, sequence-model scoring as a third evidence channel, spatial.

### Methods
1. Datasets & preprocessing (STARsolo/CellRanger as documented; whitelists; GEO/atlas accessions; Kinnex FLNC processing).
2. Read ingestion: UMI dedup, SAM-flag filtering (−F 3844), clip-site detection (`min_clip`, `min_purity`), wrong-end specificity control.
3. PAS calling: clip-seeded clustering (25 bp single-linkage, read-weighted mode, `--polya-min-reads` in distinct molecules), coverage-only tier-2, multi-PAS splitting; the legacy strategies as registry baselines.
4. Internal-priming filter: genomic A-tract definitions, `annotate` vs `filter` modes.
5. Annotation: tiered gene-end/3'UTR assignment (`max_gene_distance`, `utr_multiplier`, TIER_3 handling); note the TIER label reflects gene assignment, not 3'-end proximity.
6. Output tiers and the pre-registered default; `low_confidence` file.
7. Differential APA: within-gene Fisher with `count_mode="cells"`; NB with shrinkage; permutation calibration; replication filter; effect-size floor; marker pre-filter handling to avoid double-dipping.
8. Length metrics: PDUI-classic (proximal/distal rule, strand-correct), proportion, entropy; structural direction; `switch trend`.
9. Trusted-novel PAS definition and long-read validation protocol.
10. Validation protocol: curated PolyASite 2.0 point sites (569,005; GRCh38.96), protein-coding TES, detected-gene denominators, windows, gene-body-shuffled nulls (3 seeds), precision/recall/F1 definitions, motif and IP analyses; atlas-free runs for scoring (no snap circularity).
11. Benchmarking protocol: tool versions, parameters, environment fixes, hardware, resource measurement, `score_tool.py`.
12. PAS-profile clustering (TF-IDF/LSI/Leiden) as an optional utility, with the gene-sum ablation.
13. Reproducibility: frozen commit per run, `run_config.json`, content-hashed `run_manifest.json`, drop ledgers + survivor invariant, `reannotate` branching; code/data availability (MIT, GitHub + Zenodo DOI, all manifests as supplement).

---

## (f) Positioning vs existing single-cell APA tools (honest, reliability-centred)

**Defensible differentiators (each tied to evidence in this repo):**
1. **Evidence-tiered output with per-site molecule support.** Clip-supported vs coverage-only tiers, UMI-deduplicated clip counts in the BED score column, IP flag per site — users can see *why* a site was called. No compared tool exposes this.
2. **Pre-registered precision-first default and published gates.** The default output and its gates were fixed in writing before the final run; the stale run's failure and the post-hoc sweep are disclosed (S9).
3. **Calibrated, replicated switch testing.** Null-calibration as a shipped check, permutation-q fallback, replication filter across samples, effect-size floor — a documented answer to pseudoreplication and double-dipping.
4. **Atlas-independent trusted-novel PAS** with long-read validation against a positional null.
5. **Benchmark methodology:** one scorer, curated point-site reference, per-dataset denominators, gene-body-shuffled nulls, replicate reproducibility, long-read re-ranking, the cleavage-offset finding, and the evidence-type result.
6. **Engineering/reproducibility layer:** aligner-agnostic input (STARsolo, CellRanger, Alevin-fry — three CellRanger-input bugs found and fixed, with a fixture in CI), `reannotate` branching, append-only drop ledgers with survivor invariant, content-hashed manifests, frozen commit per run.
7. **Optional utility:** PAS-profile clustering for label-free exploration when no GEX labels exist — explicitly an expression proxy (S4).

**What we must NOT claim:**
- **Accuracy leadership** beyond what the final-run gates show; never a single cross-dataset accuracy number; never F1 without the truth set named.
- **Reproducibility leadership** — the shipped caller's advantage is against polyApipe only; scAPAtrap and Sierra are above it.
- **Clustering novelty or "cell identity from PAS usage alone"** — the ablation shows gene totals carry the signal.
- **FDR-controlled discovery** from the shipped Fisher/NB q-values; only calibrated or permutation-calibrated q-values after the fix.
- **The retired 0.9986 precision**, anything against the 18.4M dump, or atlas snap-rate as precision (circularity).
- **Stale Stage-2 numbers** as results; **quarantined spermatogenesis/Kinnex numbers** until regenerated per `11_verifier_corrections.md`.
- **Monotone four-stage spermatogenesis shortening** — the SPG step is fragile; claim three stages.
- Read1 cleavage-site precision, sequence-model filtering, spatial support — roadmap items.

**Positioning sentence for the Introduction:** "Rather than asking how many poly(A) sites a single-cell caller can report, PeakATail asks how many can be trusted: it seeds sites from poly(A)-tail read evidence, labels every call with its evidence tier and molecule support, tests cell-type switches only where they are calibrated and replicated, and validates its novel sites against long-read truth — within a single aligner-agnostic, provenance-tracked framework whose acceptance gates were registered before the results existed."

---

## Pre-submission checklist (condensed)
- [ ] Stage 1c merged (UMI-dedup clip evidence, −F 3844, BED score = molecules); final Stage-2 run launched from a frozen commit with the pre-registered default.
- [ ] Fill every [final run pending] bracket from the final run only; attach gate table (T3) with pass/fail.
- [ ] Regenerate spermatogenesis and Kinnex figures on re-keyed matrices with the `11_verifier_corrections.md` corrections as spec; resolve `wul`/`wdi`.
- [ ] FDR recalibration on counts-layer, re-keyed matrices (`--count-mode cells`, NB shrinkage); decide permutation-q default before freezing Methods.
- [ ] Amir's #67 SWITCH_CELLTYPE rerun (fixed `main.nf`, counts layer, replication filter) — decides whether Fig 6 ships.
- [ ] Re-run the novelty usage-space ablation on re-keyed matrices for S4 (does not change the claim).
- [ ] Atlas-free scoring runs for Fig 2; Kinnex support converted to molecules or stated as records.
- [ ] Real-cohort provenance reconcile (`reconcile_summary.json` all_ok) on the final Laughney run.
- [ ] Literature sweep of 2025–26 single-cell PAS callers using clip evidence or long-read truth (positioning, not novelty, risk).
- [ ] bioRxiv preprint + Zenodo DOI + Docker image before journal submission.
