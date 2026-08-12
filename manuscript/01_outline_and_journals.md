# PeakATail Manuscript Skeleton — Q1 Methods + Application Paper

**Working draft v0.1 — grounded in PeakATail v0.2.0 (repo: `tools/PeakATail`, docs: bmglab.github.io/PeakATail).** Capability claims below were checked against README.md, CHANGELOG.md, ROADMAP.md, HANDOFF.md, and `docs/` (CLI reference, strategies, tutorials). Items the current release does **not** yet do are flagged explicitly — do not let them leak into the submitted text as capabilities.

---

## (a) Ranked candidate journals

| Rank | Journal | Fit rationale | Format constraints (methods+application paper) |
|---|---|---|---|
| 1 | **Genome Biology** (IF ~12) | Best precedent fit: Sierra (Patrick 2020) and scDaPars (2021) published here; APAeval benchmark also here. "Method" article type explicitly wants new method + rigorous benchmarking + a biological application — exactly our three-pillar structure. Values honest tool comparison and reproducibility (our provenance/manifest machinery is a selling point). | Method article: no hard word limit but expect ~6–8 main figures; structured abstract (Background/Results/Conclusions, ≤350 words — our 150–200 w draft fits easily); open access APC (~$4k); code must be in a DOI-archived repository (Zenodo) under an OSI license (MIT ✓); mandatory "Availability of data and materials". |
| 2 | **Nature Communications** (IF ~15) | Takes methods papers when the application yields standalone biological insight — needs the TAM/PBMC application (Figs 5–6) to carry real novel biology, not just tool demos. Highest visibility; toughest bar: reviewers will demand the benchmarking AND a biological discovery. | Article: ~5,000 words (excluding Methods), up to 10 display items, unstructured abstract ≤150 words (our draft must be trimmed), Methods section effectively unlimited; Reporting Summary + Code Availability statement mandatory; APC ~$6.8k. |
| 3 | **Genome Research** (IF ~6–7) | "Methods/Resource" track; SCAPTURE-adjacent work appears here. Good fit if the emphasis is accuracy/benchmarking with a solid but not headline biological application. Slightly lower impact than GB with similar demands — natural second submission after GB. | Methods paper: ~55,000 characters incl. references; abstract ≤250 words; typically ≤8 figures; requires software freely available to academics, data in public repositories before review. |
| 4 | **Nucleic Acids Research** (IF ~14–16) | Strong IF; publishes computational method papers in the "Computational Biology" / "Methods" sections (PolyASite and PolyA_DB are NAR papers — reviewer pool knows PAS databases intimately, which cuts both ways: our PolyASite/PolyA_DB validation will be scrutinized by the people who built them). Fit is real but methods papers there often have a resource/webserver flavor we don't have. | Standard paper: no rigid length limit but concision expected (~9 typeset pages); abstract ≤200 words + mandatory graphical abstract; open access mandatory (APC ~$5k). Fallback within family: **NAR Genomics & Bioinformatics** (faster, lower IF ~4). |
| 5 | **Briefings in Bioinformatics** (IF ~7–9) | Published scAPAtrap and spvAPA; likes papers whose center of gravity is a comprehensive comparison against existing tools. If our benchmarking (Fig 3) ends up being the strongest asset and the biology thinner, this is the right home. | Problem-solving/method paper: ~6,000–7,000 words; abstract ≤250 words unstructured; "Key Points" box (3–5 bullets) required; figures ~6–8; code availability mandatory. |
| 6 | **Bioinformatics** (OUP) (IF ~5–6) | Published SCAPE, scDAPA. Reliable methods venue but strict length makes the three-pillar story (validation + benchmark + two applications) hard to fit — the application chapters would be compressed into supplements. Use if we descope to a leaner methods paper. | Original Paper: ~7 pages (~5,000 words all-inclusive) — tight; structured abstract (Motivation/Results/Availability); software must be freely available, installation-tested by reviewers; Application Note alternative (2 pages) is too small for this manuscript. |
| 7 | **PLoS Computational Biology** (IF ~4) | "Methods" section fits; no length pressure; strong open-science ethos matches our MIT/provenance story. Lower impact — safety net, not target. | Methods paper: no length limit; abstract ≤300 words + author summary (~200 words, lay); unlimited figures; software must be open source with archived version; APC ~$3k. |

**Recommended strategy:** write to Genome Biology's Method format (structured abstract, 6 main figures); it degrades gracefully to Genome Research and BiB without restructuring. Attempt Nature Communications only if Fig 5–6 produce a genuinely novel, orthogonally-supported biological finding (e.g., a TAM APA program invisible to gene expression).

---

## (b) Title candidates

1. **PeakATail: peak-level unsupervised clustering and differential analysis of alternative polyadenylation in single cells** *(safe, descriptive; leads with the novelty claim)*
2. **Cell-state discovery from polyadenylation site usage alone: PeakATail, a strategy-pluggable framework for single-cell APA analysis** *(leads with the conceptual claim; use only if Fig 4/Exp-2 shows PAS clustering adds populations)*
3. **PeakATail enables accurate poly(A)-site detection and expression-independent APA clustering from standard 3' scRNA-seq** *(claims accuracy + independence; needs Fig 2 numbers to back "accurate")*
4. **PeakATail: a modular toolkit for single-cell polyadenylation-site calling, clustering, and APA-switch testing, applied to the lung tumor microenvironment** *(methods+application framing for NatComms/GB)*
5. **Beyond gene counts: PeakATail resolves cell-type-specific 3'UTR isoform switching directly from single-cell poly(A) site peaks** *(punchier; risk: "beyond gene counts" over-promises unless Fig 4 delivers)*

---

## (c) Abstract draft (~180 words)

> Alternative polyadenylation (APA) diversifies 3'UTRs across cell types and is remodeled in cancer and immunity, yet standard single-cell RNA-seq analysis discards this signal by collapsing reads to gene counts. We present PeakATail, an aligner-agnostic framework that calls polyadenylation sites (PAS) directly from barcode- and UMI-tagged 3' scRNA-seq BAMs, builds PAS-by-cell count matrices, and — uniquely among current tools — clusters cells on PAS usage alone via TF-IDF/LSI and Leiden community detection, without gene-expression guidance or precomputed labels. Pluggable strategies for peak calling (including Poisson-significance and Sierra-style multi-PAS modes), reference-atlas coordinate unification, differential APA testing (Fisher and negative-binomial, FDR-controlled), 3'UTR-length metrics (PDUI, proportion, entropy), cross-dataset cluster matching, and gene-track visualization compose into one reproducible, provenance-tracked pipeline. PeakATail achieves [XX]% precision and [XX]% recall against PolyASite/PolyA_DB, matches or exceeds [Sierra/scAPAtrap/SCAPTURE] in benchmarks, and recovers established cell-type-specific APA switches in PBMCs. Applied to lung tumor microenvironment atlases, it reveals [tumor-associated-macrophage 3'UTR-shortening programs / finding TBD]. PeakATail is open source (MIT) at github.com/BMGLab/PeakATail.

*(Bracketed numbers/findings to be filled from the current rerun; trim to ≤150 words for NatComms.)*

---

## (d) Figure-by-figure plan

### Fig 1 — Tool overview and algorithm schematic
- **A.** Pipeline flow: tagged BAM(s) → per-strand read-level PAS calling → PAS×cell matrix → GTF gene-end/3'UTR annotation (tiered distance) → optional atlas snap to PolyASite → TF-IDF/LSI/Leiden clustering → `switch diff / length / trend / match / geneview`. (Base exists: `docs/assets/figures/pipeline_flow_report.png`.)
- **B.** Conceptual contrast: conventional route (gene counts → clusters → APA test on given labels) vs PeakATail (PAS counts → clusters directly). This panel carries the novelty claim.
- **C.** Peak-calling strategy registry explainer: `original` vs `lambda_poisson` (local-λ Poisson test) vs `sierra_iterative` (iterative subtraction, multi-PAS) vs `lambda_gradient` (smoothed-gradient candidates + per-summit Poisson p). (Base exists: `peak_strategy_explained.png`.)
- **D.** Worked gene example: one 3'UTR, two PAS, two cell populations using them differentially (real gene track from `switch geneview`, e.g. the existing `gene_ENSG00000103275.png`).

### Fig 2 — PAS-call validation against ground truth (sensitivity/specificity)
- **A.** Overlap of PeakATail PAS with PolyASite 2.0/3.0 and PolyA_DB v4 (upset/Venn; base exists: `peak_dual_db_validation.png`).
- **B.** Distance-to-nearest-known-PAS distribution (base exists: `peak_distance_hist.png`).
- **C.** Precision–recall across overlap windows (10/25/50/100/200 bp) and across the four peak-calling strategies; report F1. Targets from ROADMAP: precision ≥70%, recall ≥60%, F1 ≥0.65 (stretch: 80/70/0.75).
- **D.** Genomic-region distribution of calls (3'UTR / exonic / intronic / intergenic).
- **E.** Sequence-level validation: AAUAAA motif enrichment logo ±40 bp around called PAS vs shuffled controls; fraction of calls with A-rich genomic downstream stretches (internal-priming proxy) before/after filtering. **Honesty note:** `--ip-filter` is a documented no-op in v0.2.0 — either wire it before submission or perform IP assessment/filtering as an explicit post-hoc step in Methods; never describe it as an integrated pipeline stage.
- **F.** FDR calibration of the differential test: p-value histogram + QQ under null permutations; opt-in de-pseudoreplicated per-cell Fisher (`count_mode="cells"`) vs read-level default.

### Fig 3 — Benchmarking vs competing tools
- **A.** Head-to-head PAS detection (precision/recall/F1 vs PolyASite at matched windows) against Sierra, scAPAtrap, SCAPTURE (± SCAPE, scAPA) on shared datasets: mouse spermatogenesis GSE104556 (published peaks available), 10x PBMC, and the Laughney cohort. (Base exists: `peak_precision_comparison.png`, `peak_precision_v3.png`.)
- **B.** Differential-APA concordance: overlap of switch calls with tools that test APA (Sierra DTU, scAPAtrap+movAPA), plus recovery of a curated positive-control gene set.
- **C.** Compute: runtime and peak RSS vs cell count / BAM size (native from `resources.jsonl`, `tile_timings.json`; base exists: `resource_timeline.png`). Include `ema reannotate` speedup for parameter sweeps (branch downstream without re-peak-calling) — a workflow-cost argument no compared tool offers.
- **D.** Feature matrix table (may move to a Table): input compatibility, multi-PAS, peak significance, FDR, PAS-level clustering, atlas unification, cross-dataset matching, length metrics, trend testing, visualization — from `combined_polya_scrna_methods.csv` (19-tool survey). Every cell must be verifiable from the cited tool's docs.

### Fig 4 — PAS-level clustering recovers and extends known cell-type APA biology
- **A.** UMAPs: gene-expression clustering vs PeakATail PAS clustering on the same cells (PBMC + Laughney); ARI/NMI quantified (bases exist: `clustering_ari_ami_sweep.png`, `clustering_agreement_sweep.png`).
- **B.** Sankey/confusion of cluster correspondence (base exists: `match_sankey_report.png`); highlight populations split or merged only in PAS space. **Success criterion (from ROADMAP): ARI < 0.8 plus ≥1 biologically validated PAS-only population** — if this fails, the paper's framing falls back to "PAS clustering reproduces cell types from APA signal alone" (still novel, weaker).
- **C.** Recovery of established APA biology: proliferating/activated cells' global 3'UTR shortening (PDUI shift), known switch genes per lineage (e.g., immune-cell canon from the literature) — volcano + per-gene geneview tracks (bases: `diff_apa_volcano_report.png`, `volcano_0_vs_4.png`).
- **D.** 3'UTR-length landscape per cell type: PDUI-classic, proportion, Shannon-entropy distributions across clusters (bases: `pdui_distribution.png`, `proportion_distribution.png`); `switch trend` across an ordered axis (e.g., differentiation/stage) with slope + Spearman.

### Fig 5 — Application I: tumor-associated macrophages (Laughney GSE123904 + LuCA/Salcher subset)
- **A.** TAM identification and PAS-level substructure within myeloid compartment; cross-dataset canonical clusters via `switch match` (marker-overlap/MNN) between Laughney and LuCA subsets — replication across cohorts is the credibility anchor.
- **B.** Differential APA TAM vs monocytes / other myeloid states (`switch diff`, FDR-controlled; both Fisher and NB for robustness).
- **C.** Direction calls: shortening/lengthening per gene (structural direction from `length` output); pathway enrichment of switched genes.
- **D.** Gene-track vignettes (geneview, plotly/matplotlib) of the 2–3 headline genes; orthogonal support (bulk 3'-seq, long-read, or published atlas evidence) for at least one.

### Fig 6 — Application II: PBMC insight (or fold into Fig 5 if space demands)
- **A.** PAS-clustering of public 10x PBMC; annotation via markers.
- **B.** Cell-type-specific PAS switches among T/B/NK/monocyte lineages; concordance with published immune APA findings plus any novel switch.
- **C.** Cross-sample reproducibility (multiple PBMC runs; atlas-snap unified coordinates).

### Supplementary figures
- S1. Dataset summaries (cells, depth, chemistry, QC; UMAPs with reference annotations).
- S2. Peak-call QC: width, per-gene PAS counts, coverage distributions (`peak_qc_default.png`, `peak_counts.png`).
- S3. Strategy comparison in-framework: 4 peak callers on identical regions (counts, DB agreement, runtime).
- S4. Clustering robustness: resolution/n-neighbors/seed sweeps via `reannotate` branching; silhouette; subsampling stability.
- S5. Statistical validation extended: `count_mode="cells"` vs `"reads"`; `nb_multi` `sample_split`; min-cells sensitivity.
- S6. Atlas-snap behavior: snap-rate (reported as coverage, with the documented circularity caveat — snap-rate against the same atlas used for snapping is NOT precision), distance-to-summit distributions, with/without-atlas cluster stability.
- S7. Internal-priming characterization of called PAS (post-hoc), motif classes of novel (non-database) PAS + conservation.
- S8. Additional geneview vignettes; interactive HTML as supplementary files.
- S9. Compute scaling detail; `collapse`/`combine` multi-library workflows.

### Supplementary tables
- T1. Full 19-tool feature comparison (from `combined_polya_scrna_methods.csv`).
- T2. All differential APA results per contrast (findings_long export).
- T3. PDUI/length scores per gene × cluster (length_long export).
- T4. Validation metrics per dataset × strategy × window.
- T5. Software versions, parameters, `run_config.json` / `run_manifest.json` digests for every run (provenance ledgers make this nearly free — say so in Methods as a reproducibility feature).

---

## (e) Section-by-section outline

### Introduction (~5 paragraphs)
1. APA biology: 3'UTR isoform choice controls stability/localization/translation; remodeled in proliferation, cancer, immune activation.
2. 3' scRNA-seq reads pile at poly(A) sites — the signal is in every 10x dataset but standard pipelines collapse it to gene counts.
3. Existing single-cell APA tools (Sierra, scAPAtrap, SCAPTURE, SCAPE, scMAPA, scDaPars, scTail, …): all quantify APA *given* cell identities derived from gene expression or supplied annotations; none clusters cells on PAS usage alone; ecosystem is fragmented (separate tools for calling, testing, visualization) and inputs are often CellRanger-locked.
4. Gap → PeakATail: one aligner-agnostic, strategy-pluggable, provenance-tracked framework; expression-independent PAS clustering as the conceptual core.
5. Contributions summary: validated accuracy, benchmark, recovered biology, TME/PBMC application, open source.

### Results

**R1. PeakATail: an end-to-end, strategy-pluggable framework for single-cell APA** *(descriptive; Fig 1)*
- **Key claim:** PeakATail turns any CB/UB-tagged BAM into PAS-level single-cell analysis through one modular pipeline.
- **Evidence:** architecture description; four peak-calling strategies incl. Poisson-significance (`lambda_poisson`, `lambda_gradient`) and multi-PAS (`sierra_iterative`, `lambda_gradient`); TF-IDF/LSI/Leiden PAS clustering; atlas snapping; `reannotate` branching; provenance (drop ledgers, manifests, reconcile invariant). No experiment needed beyond a worked example run.

**R2. PeakATail calls PAS with high precision and recall against independent ground truth** *(Fig 2)*
- **Key claim:** ≥70% precision / ≥60% recall (window ≤50 bp) vs PolyASite + PolyA_DB, with canonical PAS sequence features and low internal-priming contamination.
- **Evidence needed:** run current rerun outputs (`pasbed.bed`) through the dual-DB overlap protocol at multiple windows; motif logo + IP-proxy analysis on genome FASTA; per-strategy breakdown; FDR calibration by permutation. **Blocking dependency:** validation must use an atlas-*free* run when scoring against PolyASite (else circular — HANDOFF D3 caveat); characterize novel PAS (motif/conservation) rather than counting them as errors.

**R3. PeakATail matches or exceeds existing tools at lower or comparable cost** *(Fig 3)*
- **Key claim:** detection accuracy ≥ Sierra/scAPAtrap/SCAPTURE on shared datasets; competitive runtime/memory; unique workflow economics via `reannotate`.
- **Evidence needed:** run all comparison tools on GSE104556 + PBMC with their default settings and identical references; identical scoring script for all; APAeval submission if timeline allows; resource logs from `resources.jsonl`. **Honesty rule:** report where competitors win too (e.g., if scAPAtrap precision is higher, say so and show recall/F1/cost trade-off).

**R4. Clustering on PAS usage alone recovers cell identity and exposes structure invisible to gene counts** *(Fig 4; the novelty-carrying section)*
- **Key claim (primary):** PAS-only clustering reproduces annotated cell types (ARI/NMI quantified) — cell identity is encoded in 3'-end choice.
- **Key claim (stretch, keep only if it survives):** ≥1 population or sub-state separable in PAS space but not in gene space, validated by orthogonal markers.
- **Evidence needed:** paired gene-based (Seurat/scanpy standard) vs PAS-based clustering on PBMC + Laughney; ARI/NMI + Sankey; marker-based identity of any PAS-only cluster; robustness sweeps (S4). Pre-register the ARI < 0.8 + validation criterion internally; if unmet, retitle section to the primary claim only.

**R5. PeakATail recovers established cell-type-specific APA programs** *(Fig 4C–D)*
- **Key claim:** known biology reproduced — global 3'UTR shortening in proliferative/activated states; canonical lineage-specific switch genes detected at FDR < 0.05 with correct direction.
- **Evidence needed:** curated positive-control gene list from literature (immune + cancer APA canon) assembled *before* looking at results; `switch diff` + `switch length` direction calls; geneview tracks for headline genes.

**R6. Application: APA remodeling in tumor-associated macrophages across two lung cancer cohorts** *(Fig 5)*
- **Key claim:** TAMs exhibit a reproducible APA program (e.g., 3'UTR shortening of [gene set TBD]) relative to blood monocytes/other myeloid states, replicated in Laughney and LuCA-subset cohorts.
- **Evidence needed:** myeloid-compartment PAS clustering; cross-cohort `switch match` canonical clusters; per-contrast diff + length with direction; replication statistics (effect-size correlation across cohorts); orthogonal support for ≥1 gene. This is the section that decides NatComms viability.

**R7. Application: PBMC cell-type APA atlas** *(Fig 6; can merge into R5/R6)*
- **Key claim:** a reusable per-cell-type PAS/PDUI reference across PBMC lineages, concordant across samples.
- **Evidence needed:** multi-sample PBMC runs, atlas-snapped unified coordinates, reproducibility metrics.

### Discussion
- What expression-independent PAS clustering buys, and when to use it vs gene clustering.
- Interpretation of TAM/PBMC findings against APA-in-cancer/immunity literature.
- **Limitations (state plainly):** 3'-tag data limits precision of cleavage-site localization (no Read1 mode yet); internal-priming filtering currently post-hoc, not an integrated stage; peak-level significance only in two of four strategies; NB strategies scale poorly over all cluster pairs (mitigated by `--cluster-pairs`, marker pre-filter); barcode correction delegated to the aligner; per-cell Fisher/NB improvements are opt-in pending default-flip characterization.
- Future: Read1 cleavage-site mode, integrated IP filter, deep-learning PAS validation, spatial.

### Methods
1. Datasets & preprocessing (STARsolo protocol as documented; whitelists; GEO/atlas accessions).
2. PAS calling: per-strand region building, the four strategies with hyperparameters and default values; PAS-summit merger; `pas_gap`.
3. Annotation: tiered gene-end/3'UTR assignment (`max_gene_distance`, `utr_multiplier`, TIER_3 handling).
4. Atlas mode: snap algorithm, `atlas_distance`, mapping provenance; when it is and isn't used (validation runs atlas-free).
5. Clustering: TF-IDF + LSI + Leiden (and libsize variant); depth-correlation component dropping; seeds.
6. Differential APA: within-gene Fisher; NB pairwise/multi; marker pre-filtering; FDR; `count_mode` and `sample_split` variants.
7. Length metrics: PDUI-classic (proximal/distal selection rule), proportion, Shannon entropy; structural direction call; `switch trend` model.
8. Cross-dataset matching: marker-overlap, Jaccard, MNN; canonical clusters.
9. Validation protocol: databases/versions, overlap windows, precision/recall/F1 definitions, motif and IP analyses, permutation nulls.
10. Benchmarking protocol: tool versions, parameters, hardware, resource measurement.
11. Reproducibility: `run_config.json`, `run_manifest.json` with content hashes, drop ledgers + survivor invariant, `reannotate` branching; code/data availability (MIT, GitHub + Zenodo DOI, all run manifests as supplement).

---

## (f) Novelty positioning vs existing single-cell APA tools (honest, from the docs)

**Defensible novel claims (each verified in repo docs):**
1. **Expression-independent PAS-level clustering.** TF-IDF/LSI + Leiden run directly on the PAS×cell matrix — cell states defined by 3'-end usage without gene-expression guidance or precomputed labels. ROADMAP records this as the confirmed first among 19 surveyed tools; must be re-verified against 2025–26 literature immediately before submission (stated ROADMAP risk).
2. **Strategy-pluggable architecture with in-framework baselines.** Registered, swappable strategies at every stage (4 peak callers incl. a Sierra-style reimplementation, 3 clustering methods, 3 diff tests, 3 length metrics, 3 matching methods) — enables like-for-like algorithm comparison inside one pipeline, which no compared tool offers.
3. **Aligner-agnostic input.** Any CB/UB-tagged BAM (STARsolo, CellRanger, Alevin-fry) — several competitors are CellRanger-output-locked.
4. **Cross-dataset machinery as first-class features:** PolyASite atlas snapping into a unified PAS coordinate space + `switch match` (marker-overlap/Jaccard/MNN) canonical clusters + `switch combine`/`collapse` for multi-library designs.
5. **`reannotate` branching:** re-run trim/filter/clustering variants from one peak-call at near-zero cost — makes parameter-sensitivity analysis (and honest reporting thereof) practical.
6. **`switch trend`:** ordered-covariate APA-length trend testing (slope + Spearman + direction) — beyond pairwise designs.
7. **Statistical care on pseudoreplication:** opt-in per-cell (de-pseudoreplicated) Fisher contingency and split-sample NB selection/inference — a documented answer to a known field-wide flaw in read-level testing.
8. **Engineering/reproducibility layer:** append-only PAS/cell drop ledgers with a survivor invariant, content-hashed run manifests, per-run resource telemetry — position as "auditable APA analysis", attractive to methods reviewers.

**What we must NOT claim (current v0.2.0 reality):**
- An **integrated internal-priming filter** — `--ip-filter`/`--genome-fasta`/`--annot-filter` are accepted but explicitly no-op with warnings (CHANGELOG 0.2.0). Wire it or present IP analysis as post-hoc Methods.
- Built-in **database-validation/benchmark subcommands** — `--validate-db`/`--benchmark` are likewise no-ops; validation runs via external scripts (fine, but describe honestly).
- **Read1-based cleavage-site precision**, **K-NN imputation**, **deep-learning PAS filtering**, **spatial support** — all roadmap items, not features. Mention only in Discussion/future work.
- **Peak-level FDR across all strategies** — Poisson significance exists in `lambda_poisson`/`lambda_gradient` only; `original` and `sierra_iterative` use height thresholds.
- Superiority claims not backed by Fig 3 numbers; and never present atlas snap-rate as precision (documented circularity caveat).

**Positioning sentence for the Introduction:** "Existing single-cell APA tools quantify polyadenylation *given* cell identities defined by gene expression; PeakATail inverts this, asking whether — and showing that — cell identity can be read from polyadenylation-site usage itself, within a single modular, aligner-agnostic, and fully provenance-tracked framework."

---

## Pre-submission checklist (from ROADMAP/HANDOFF, condensed)
- [ ] Wire or post-hoc-document internal-priming filtering (blocks Fig 2E claims).
- [ ] Atlas-free validation runs for Fig 2 (avoid circularity).
- [ ] Characterize `count_mode="cells"` / `sample_split=True` on the no-atlas rerun; decide defaults before freezing Methods.
- [ ] Real-cohort provenance reconcile (`reconcile_summary.json` all_ok) on the full Laughney rerun.
- [ ] Fresh literature sweep for any 2025–26 tool doing PAS-only clustering (novelty risk).
- [ ] bioRxiv preprint + Zenodo DOI + Docker image before journal submission.
