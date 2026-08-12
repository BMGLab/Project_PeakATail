# Benchmarking Component Design — PeakATail (`ema`) Manuscript

Prepared 2026-08-12. Facts below (publication venues/years, GitHub activity) were verified via web search and the GitHub API on this date. GitHub "last push" dates are exact API values.

---

## (a) Survey of single-cell APA / PAS detection tools

### Primary tools (call and/or quantify PAS from 3'-tag scRNA-seq)

| Tool | Year (venue) | PAS-calling approach | Quantification unit | Differential-APA test | Language | Actively maintained? (last GitHub push) |
|---|---|---|---|---|---|---|
| **scAPA** | 2019 (NAR; Shulman & Elkon) | De novo: cluster-pooled 3' read pileup peaks, wide peaks split by Gaussian mixture (mclust) | Peak × cell (and per-cluster) counts | Chi-square on proximal-PAS usage; proximal usage index | R + shell | **No** (2021-12) |
| **Sierra** | 2020 (Genome Biology; Patrick et al.) | De novo: spline fitting on per-gene 3' coverage to find peaks | Peak × cell UMI counts | DEXSeq-based differential usage (`DUTest`) | R | Low (2023-06) |
| **polyApipe / polyApiper** | 2019– (unpublished; Monash Bioinformatics Platform) | De novo: reads with non-templated poly(A) soft-clips grouped into poly(A) "ranges" | Peak(range) × cell UMI counts | Downstream in polyApiper (weitrix-based shift analysis) | Python + R | **Yes** (2025-02) |
| **scAPAtrap** | 2021 (Briefings in Bioinformatics; Wu et al.) | De novo: coverage peak finding + poly(A)-tail-anchored reads → near-nt resolution | PAS × cell counts | None built in (pairs with movAPA/scDAPA) | R | **Yes** (2025-10) |
| **SCAPTURE** | 2021 (Genome Biology; Li et al., YangLab) | De novo peak calling + **DeepPASS** deep-learning classifier to evaluate/filter PAS | PAS × cell UMI counts (featureCounts) | No native test (export to Seurat; usage-ratio tests) | Shell + Python | Low (2024-01) |
| **MAAPER** | 2021 (Genome Biology; Li, Zheng, Wang, Tian) | Annotation-based: probabilistic assignment of near-site reads to **PolyA_DB v3** PAS | Predicted PAS usage fractions per gene/condition | Built-in likelihood-based tests (RED/RLDu shift statistics) | R | **No** (2021-08) |
| **SAPAS** | 2021 (BMC Biology; Yang et al.) | De novo: 3'-tag cluster calling incl. novel sites, internal-priming filtering | PAS × cell counts | Cell-type-specific APA tests per paper (Fisher-type) | R | **No** (2021-04) |
| **scDaPars** | 2021 (Genome Research; Gao & Li) | No site calling: DaPars2 two-isoform coverage-ratio regression + NMF imputation | Gene-level PDUI per cell | Differential PDUI between clusters | R | **No** (2021-05) |
| **scMAPA** | 2022 (GigaScience; Bai et al.) | No de novo calling: 3'UTR coverage fit by quadratic programming per cluster BAM | Gene-level long/short isoform abundances | Logistic regression LRT for cluster-specific APA | R + Python | **No** (2021-11) |
| **SCAPE** | 2022 (NAR e66; Zhou et al., LuChenLab) | De novo: EM mixture model using cDNA **insert-size** + read position; PAS as Gaussians | PAS × cell (weighted) counts | Companion scripts (usage shifts); no formal native test | Python | Low (2024-03) |
| **scUTRquant** | 2024 (Nature Communications; Fansler, Mitschka & Mayr) | Annotation-based: kallisto/bustools against curated cleavage-site-augmented "UTRome" (200+ cell types; +40% sites vs GENCODE) | 3'UTR isoform × cell counts | **scUTRboot** bootstrap tests (WUI/LUI shifts) | Snakemake + R | **Yes** (2026-01) |
| **Infernape** | 2023 (Genome Research; Kang et al.) | De novo within-gene peak inference (change-point style) + internal-priming filtering | PAS × cell counts | Per-cluster usage tests (Fisher/chi-square) | R | Sporadic (2025-05) |
| **scraps** | 2022 preprint (bioRxiv; Fu, Gillen, Hesselberth et al.) | Read1-based: TVN-priming position from paired R1 → near-nucleotide priming-site calls; flags internal priming | Site × cell counts | DEXSeq via `scrapR` | Snakemake (Python/R) | **Yes** (2026-07); still unpublished |
| **scTail** | 2025 (Genome Biology; Hou & Huang, StatBiomed) | Read1-based: PAS from R1 cDNA 3' ends (STARsolo), pretrained sequence model filters false positives; R2 for quantification | PAS × cell counts | Built-in alternative-usage test | Python | **Yes** (2025-06) |
| **scPAISO** | 2025 preprint (bioRxiv, Aug 2025) | Read1-based: direct mRNA cleavage-site capture from R1; R2 anchored to PAS peaks | PAS-isoform × cell counts | Per paper (usage-ratio) | Python | New (preprint) |
| **scDeepAPA** | 2026 (Briefings in Bioinformatics) | Deep learning (CNN + Mamba state-space + BiLSTM) trained on PolyASite v3.0 for PAS identification | PAS × cell quantification | Framework includes functional/DE interpretation | Python | New |
| **SCINPAS** | 2023–2025 (Zavolan lab; workflow behind **PolyASite v3.0**, NAR 2025) | De novo from scRNA-seq incl. non-templated poly(A)-tail reads, dedicated IP filtering | PAS clusters (atlas-building) | n/a (atlas workflow) | Snakemake/Python | Maintained (Zavolan lab) |
| **PeakATail (`ema`)** — this work | 2026 (BMGLab) | De novo per-strand streaming peak calling on CB/UB-tagged BAMs; GTF 3'UTR-tiered annotation; atlas & internal-priming filter modes | PAS × cell UMI counts | `switch diff` (Fisher exact / NB regression, FDR), `switch length`, `switch match`, `switch trend` | Python | **Yes** (2026-08) |

Notes on the required columns: "Maintained?" = repository push within ~18 months of Aug 2026. scAPAtrap's continued pushes reflect the Wu (BMILAB) ecosystem (movAPA/vizAPA); Sierra and SCAPTURE remain installable but see little development; MAAPER/scDaPars/scMAPA/SAPAS/scAPA are effectively frozen.

### Adjacent / ecosystem tools (do not call PAS, or operate downstream — cite, don't benchmark)

| Tool | Year | Role |
|---|---|---|
| scDAPA | 2019 | Detects/visualizes dynamic APA from annotation, no PAS calling |
| movAPA / vizAPA | 2021 / 2024 | Downstream modeling & visualization of PAS matrices (Wu lab) |
| scAPAmod | 2022 | Profiles APA "modalities" (uni/bi-modal usage) in single cells |
| ReadZS | 2022 (Genome Biology) | Annotation-free read-distribution z-score; detects 3'-end shifts w/o PAS coordinates |
| spvAPA | 2025 (Briefings in Bioinformatics) | Supervised APA analysis for sc + spatial data |
| metaAPA | 2026 (Bioinformatics Advances) | **Consensus integration of PAS predictions across tools** — direct evidence that inter-tool discordance is a recognized problem |
| scPASU | 2026 (STAR Protocols) | Protocol wrapper for PAS-usage quantification from 3' scRNA-seq |
| scAPAdb / scAPAatlas | 2022 (NAR) | Databases of single-cell PAS/APA |
| Bulk-only (excluded) | — | DaPars2, APAlyzer, QAPA, TAPAS, REPAC, APAIQ, PolyAMiner-Bulk, LABRAT etc. |

---

## (b) Published benchmarks of single-cell APA tools

1. **NAR 2026 — "Benchmarking computational methods for identifying and quantifying polyadenylation sites from 3'-tag-based single-cell RNA-seq data"** (NAR 54(9):gkag490). The most comprehensive to date. **Tools (10):** de novo — scAPA, polyApipe, Sierra, scAPAtrap, SCAPE; annotation-based — MAAPER, SCAPTURE, Infernape, scUTRquant, scraps. **Data:** 9 simulated datasets (human, mouse, Arabidopsis) + 25 real samples over four protocols (10x, CEL-seq, Drop-seq, Microwell-seq); references: GENCODE, PolyA_DB 3, plant PAS DBs; matched bulk 3'-seq for validation. **Metrics:** precision/recall/accuracy at 0–100 nt distance windows; MAPE for quantification; chi-square on nucleotide composition; poly(A)-signal presence; DeepPASTA validation of novel sites; expression correlations; ARI/silhouette clustering; precision/recall for differential APA; runtime/memory. **Conclusions:** annotation-based methods more positionally accurate; polyApipe best sensitivity/accuracy balance among de novo; SCAPE, scUTRquant, SCAPTURE, polyApipe best quantifiers; overall recommendation **SCAPTURE + polyApipe**; most tools degrade badly off-10x (CEL-seq).
2. **bioRxiv 2024.10.15.618405 — "Benchmarking alternative polyadenylation detection in single-cell and spatial transcriptomes."** **Tools (9):** scAPA, scAPAtrap, Sierra, SAPAS, scDaPars, SCAPTURE, MAAPER, SCAPE, Infernape. **Data:** 7 3'-tag protocols; a protocol-aware **simulation pipeline** preserving peak-shape characteristics; **matched ONT long reads sharing cell barcodes/UMIs with the Illumina libraries as per-cell ground truth**. **Metrics:** PAS identification, quantification, DE-APA detection, sequencing-artifact (internal-priming) filtering, computational efficiency.
3. **bioRxiv 2024.11.29.626111 — "Guidelines for alternative polyadenylation identification tools using single-cell and spatial transcriptomics data."** Categorizes tools into three algorithmic classes (alignment-based being largest, with four sub-types); four sc/spatial datasets with **matched nanopore ground truth**; deep-dives Sierra, scAPAtrap, SCAPE. Findings: **no method wins everywhere**; SCAPE least sensitive to read-length/depth changes.
4. **NAR Genomics & Bioinformatics 2025 (PMC12076406)** — bulk-tool benchmark (APAlyzer, DaPars2, TAPAS, REPAC, APAIQ, PolyAMiner-Bulk, …) but methodologically a template we mirror: precision/sensitivity vs **PolyASite 2.0 + PolyA_DB 3 + GENCODE at 25–150 bp distance thresholds**, stratified by single- vs multi-PAS genes; found large inter-tool variation; deep-learning tools call many more sites with false-positive concerns; recommends multi-tool integration.
5. **metaAPA (Bioinformatics Advances, 2026)** — quantifies low overlap between per-tool PAS predictions and builds consensus sets; useful citation for motivating our precision-focused design and a possible consensus baseline.

Takeaways for our design: (i) distance-tolerance P/R against curated atlases is the accepted primary metric; (ii) matched long reads are the emerging gold standard and address atlas incompleteness; (iii) internal-priming behavior and protocol robustness are recognized differentiators; (iv) SCAPTURE, polyApipe, SCAPE are the standing performance references to beat; (v) no existing benchmark includes PeakATail — we must run the comparison ourselves with pre-registered metrics.

---

## (c) Recommended head-to-head panel (6 tools)

Selection criteria: runnable in 2026 on standard 10x BAMs/FASTQs, representative of each algorithmic family, and either top performers in published benchmarks or field-standard baselines.

| Tool | Family represented | Why include |
|---|---|---|
| **Sierra** | De novo coverage-peak (spline) | Most widely used/cited peak caller; native DEXSeq DU test makes it the key comparator for `ema switch diff`; still installable |
| **scAPAtrap** | De novo peak + poly(A)-read evidence | Nucleotide-resolution claims parallel to PeakATail; actively maintained (Oct 2025); deep-dived in the 2024 Guidelines benchmark |
| **SCAPE** | Model-based EM (insert-size mixture) | Top quantifier and most protocol-robust in NAR 2026 and Guidelines benchmarks; methodologically orthogonal to peak calling |
| **SCAPTURE** | De novo + deep-learning PAS filter | One of two overall winners in NAR 2026; its DeepPASS filter is the natural comparator for PeakATail's ip-filter mode |
| **polyApipe** | Poly(A)-soft-clip evidence | Second overall winner in NAR 2026; best sensitivity/accuracy balance among de novo tools; maintained (Feb 2025) |
| **scUTRquant** | Annotation-based pseudoalignment | Represents the annotation-driven paradigm (upper bound on positional accuracy, no novel-site discovery); published Nat Commun 2024; actively maintained; scUTRboot gives a principled DE comparator |

**Deliberate exclusions (state in manuscript):** scAPA, MAAPER, scDaPars, scMAPA, SAPAS — unmaintained since 2021 and/or not PAS-resolution (gene-level PDUI), already shown mid/low-tier in published benchmarks; Infernape — runnable but low adoption and overlapping methodology with included de novo tools (optional 7th if reviewers ask); scraps/scTail/scPAISO — require **read1-preserved FASTQs** (TVN-priming or R1 cDNA), which standard archived 10x data (incl. Laughney GSE123904 BAMs) lack; we will note them as complementary, not competing, inputs; scDeepAPA — published May 2026, include in survey and discussion, attempt to run only if trained models are released and stable. MAAPER may be added solely to the differential-APA concordance module (it is annotation-based and trivially runnable) if a 7th data point is desired.

---

## (d) Benchmark protocol for PeakATail

All runs pinned in versioned containers (Docker/Apptainer), orchestrated by Snakemake, same machine, fixed thread budget (8 threads), 3 random seeds where sampling is involved; configs and per-tool parameter files in the paper's repo. Default parameters for all tools unless a tool's own paper specifies 10x presets (record every deviation).

### D0. Datasets and ground-truth construction

- **Primary real data:** (1) 10x Genomics PBMC reference (e.g., 10k PBMC, v3 chemistry) — deeply characterized, used by nearly every prior tool paper; (2) Laughney et al. lung adenocarcinoma/metastasis (GSE123904) — our application dataset; (3) one benchmark dataset with **matched ONT long reads sharing CB/UMIs**, reused from bioRxiv 618405 / 626111 GEO deposits (retrieve accessions), or the public PacBio Kinnex single-cell PBMC set as long-read truth.
- **Ground-truth tiers (all lifted to GRCh38, strand-aware, ±collapsing of clustered atlas sites to representative site):**
  - **GT1 — PolyASite 2.0** (Herrmann et al., NAR 2020): orthogonal, bulk 3'-end-seq-derived. Primary atlas.
  - **GT2 — PolyA_DB v3** (Wang et al., NAR 2018; hg19 → liftOver): 3'READS-derived. Secondary atlas; intersection GT1∩GT2 = "high-confidence" set, union = "permissive" set.
  - **GT3 — matched long-read PAS:** ONT/Kinnex read 3' ends with soft-clipped poly(A) ≥ 10 nt, clustered within ±20 nt, supported by ≥3 UMIs; per-cell-type truth via shared barcodes.
  - **Caveat to state explicitly:** PolyASite **v3.0** (Moon et al., NAR 2025) is itself inferred from 10x scRNA-seq via SCINPAS — using it as primary truth would be circular for 10x-based callers; we use it only as a supplementary annotation, and this caveat is itself a point our Methods can make against naive benchmarking.
- **Recall denominator restriction:** recall is computed only over atlas PAS in **expressed** genes (gene ≥ 20 UMIs total and ≥ 10 cells in the dataset), within annotated or extended (≤5 kb, matching `--max-gene-distance`) 3'UTR/terminal-exon space. Report the denominator size; unexpressed-gene sites must not deflate recall.

### D1. PAS-level precision / recall vs atlases (distance-tolerance windows)

- Matching: strand-aware `bedtools closest` between each tool's called PAS (single representative coordinate per site; for interval callers like Sierra use the peak mode/3' boundary — record the rule per tool) and GT1/GT2/GT3.
- Windows: d ∈ {10, 25, 50, 100, 150, 200} nt.
- **Metrics:** Precision(d) = called PAS within d of an atlas site / all called; Recall(d) = expressed-gene atlas sites recovered / denominator; F1(d); headline **F1 at ±50 nt vs GT1**; positional accuracy = median |signed offset| to nearest atlas site (and offset bias direction); total sites called; % "novel" calls (no support at ±200 in GT1∪GT2), and of those, % validated by GT3 long reads and % carrying a canonical poly(A) hexamer (AAUAAA/AUUAAA + 12 minor variants) at −50..−10 nt; nucleotide-composition profile ±100 nt around called sites (expect A-rich downstream dip, U/GU-rich DSE).
- **Plots/tables:** Fig B1a — F1 vs window-size curves, one line per tool; Fig B1b — precision–recall scatter at ±50 nt (PeakATail highlighted, atlas-based tools flagged as having an inherent advantage); Fig B1c — signed-offset density per tool; Fig B1d — % hexamer-supported and % long-read-validated novel sites (stacked bars). Table B1 — full P/R/F1 grid (tool × window × atlas).

### D2. Sensitivity under read- and cell-downsampling

- **Read downsampling:** `samtools view -s` at 5/10/25/50/75/100% of reads, 3 seeds.
- **Cell downsampling:** random barcode subsets of 250 / 500 / 1,000 / 2,500 / 5,000 / all cells, 3 seeds.
- **Metrics per level:** (i) self-recall — fraction of the tool's own full-depth ±50 nt-matched-to-GT1 sites recovered (robustness); (ii) recall vs GT1 at ±50 nt; (iii) number of sites called (saturation curve); (iv) quantification stability — Spearman correlation of pseudobulk PAS-usage fractions vs full depth.
- **Plots:** Fig B2a — recall-vs-depth curves (reads); Fig B2b — recall-vs-cells curves; Fig B2c — usage-correlation heatmap (tool × depth). Mean ± SD over seeds.

### D3. Specificity via internal-priming (IP) negative controls

- **Negative region set:** genomic positions in expressed gene bodies with ≥6 consecutive genomic A (or ≥12 A in a 20-nt window) on the transcribed strand, located >200 nt from any GT1∪GT2 site — these are canonical IP traps. (This mirrors the emaout `negbed/negmatrix` machinery already in the pipeline.)
- **Metrics:** IP-call rate = % of a tool's called PAS falling in negative regions; A-rich unsupported fraction = % of calls with downstream (0..+20 nt) genomic A-fraction > 0.65 and no atlas/long-read support; per-tool **sensitivity-vs-IP-rate operating point**.
- **PeakATail-specific experiment:** run identical inputs with ip-filter mode off/on (branched cheaply via `ema reannotate` where applicable) → paired points with arrows showing the filter moves down in IP-rate at minimal recall cost; same comparison for SCAPTURE with/without DeepPASS filtering, scraps-style annotation of IP events noted in discussion.
- **Plots:** Fig B3 — recall (±50 nt, GT1) on y vs IP-call rate on x; one point per tool, arrows for filter on/off. Table B3 — negative-region counts and rates.

### D4. Runtime and peak memory

- `/usr/bin/time -v` (wall-clock, max RSS) per stage: peak calling / quantification / differential testing, separately; identical host (record CPU, RAM), 8 threads, warm cache.
- Scaling series: 5k / 10k / 25k cells (subsampled BAMs) and 100M / 250M / 500M reads.
- **PeakATail differentiator to report:** cost of a *re-analysis sweep* (10 parameter variants) — competitors must rerun end-to-end; `ema reannotate` branches from one peak-call, so report sweep totals as well as single-run cost.
- **Plots/tables:** Table B4 — wall-time and peak RSS per tool × stage × input size; Fig B4 — log–log runtime vs cells with per-tool slopes; inset: 10-variant sweep totals.

### D5. Quantification and clustering fidelity (supporting module)

- Pseudobulk PAS-usage fractions per cell type vs GT3 long-read usage: Pearson/Spearman r and MAPE (only sites matched at ±50 nt; restrict to multi-PAS genes).
- Cell-embedding utility: ARI and NMI of each tool's PAS-matrix Leiden clustering (PeakATail: TF-IDF+LSI pipeline) against expression-derived cell-type labels on PBMC.
- **Plots:** Fig B5a — per-tool r/MAPE bars; Fig B5b — ARI/NMI bars.

### D6. Differential-APA concordance across tools

- **Contrasts:** PBMC monocytes vs T cells (large, canonical); Laughney tumor epithelium vs normal epithelium (application-relevant); one simulated contrast with known DE-APA truth (reuse the protocol-aware simulator of bioRxiv 618405 or the NAR 2026 simulated datasets) for absolute precision/recall/AUPRC.
- **Design (i) — native end-to-end:** PeakATail `ema switch diff` (Fisher and NB variants) vs Sierra DUTest (DEXSeq), scUTRquant+scUTRboot (WUI/LUI), polyApiper shift test, SCAPE/SCAPTURE usage tests (or Wilcoxon on usage where no native test), gene-level significant sets at FDR 0.05.
- **Design (ii) — harmonized testing:** feed every tool's PAS×cell matrix into one common DEXSeq (and one common Fisher) framework → isolates quantification differences from statistical-test differences. This two-design structure is the key novelty over prior benchmarks.
- **Metrics:** pairwise Jaccard of significant gene sets; Spearman correlation of gene-level effect sizes (ΔPAS-usage / Δ weighted 3'UTR length) on shared genes; on simulation — precision/recall/AUPRC vs planted DE-APA; **positive-control recovery**: % recovery of (a) canonical proliferation-associated 3'UTR shortening (Sandberg 2008; Mayr & Bartel 2009 gene sets) in tumor-vs-normal, and (b) PBMC cell-type APA events reported by ≥2 independent prior single-cell APA papers (scAPA, SCAPTURE, MAAPER, scUTRquant analyses) — pre-registered gene list frozen before running.
- **Cross-check with `ema switch length`:** direction-of-effect agreement (global shortening in tumor/proliferative compartments) reported as a sign-concordance percentage.
- **Plots/tables:** Fig B6a — UpSet plot of significant APA genes across tools; Fig B6b — effect-size correlation heatmap; Fig B6c — PR curves on simulated truth; Fig B6d — positive-control recovery bars; Table B6 — per-tool counts, FDR settings, versions.

### Headline deliverables mapped to manuscript

| Item | Content |
|---|---|
| Main Fig. "Benchmark" (2 rows) | B1b PR scatter + B1a F1-window curves; B3 specificity plot; B4 runtime/memory; B6a UpSet |
| Supp. Figs | offset densities, hexamer/nucleotide profiles, downsampling curves, clustering ARI, simulation PR |
| Supp. Tables | full P/R/F1 grid; tool versions/parameters/commands; runtime table; DE gene lists; positive-control list (pre-registered) |
| Repro artifact | Snakemake benchmark repo + containers + config, cited in Data/Code Availability |

**Success criteria to claim in the paper (pre-registered):** PeakATail F1(±50 nt, GT1) ≥ best de novo competitor; median |offset| ≤ 50 nt; IP-call rate with ip-filter ≤ SCAPTURE-with-DeepPASS; graceful recall decay (≥80% self-recall at 25% reads); runtime/memory within 2× of fastest peak-based tool with a ≥5× advantage on parameter sweeps via `reannotate`; DE-APA recovery of positive controls ≥ every competitor at matched FDR.

---

## Sources

- [Benchmarking computational methods for identifying and quantifying polyadenylation sites from 3'-tag-based single-cell RNA-seq data — NAR 2026](https://academic.oup.com/nar/article/54/9/gkag490/8676205)
- [Benchmarking alternative polyadenylation detection in single-cell and spatial transcriptomes — bioRxiv 2024.10.15.618405](https://www.biorxiv.org/content/10.1101/2024.10.15.618405v1)
- [Guidelines for alternative polyadenylation identification tools using single-cell and spatial transcriptomics data — bioRxiv 2024.11.29.626111](https://www.biorxiv.org/content/10.1101/2024.11.29.626111v1.full)
- [Benchmarking of methods that identify alternative polyadenylation events in single-/multiple-polyadenylation site genes — NAR Genomics & Bioinformatics 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12076406/)
- [metaAPA: integration of PolyA site predictions — Bioinformatics Advances 2026](https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag147/8694770)
- [scTail — Genome Biology 2025](https://link.springer.com/article/10.1186/s13059-025-03710-7) · [GitHub](https://github.com/StatBiomed/scTail) · [preprint](https://www.biorxiv.org/content/10.1101/2024.07.05.602174v1.full)
- [scUTRquant — Fansler et al., Nature Communications 2024 (PubMed 38744866)](https://pubmed.ncbi.nlm.nih.gov/38744866/) · [GitHub Mayrlab/scUTRquant](https://github.com/Mayrlab/scUTRquant)
- [Infernape — Genome Research 2023](https://genome.cshlp.org/content/33/10/1774) · [PubMed 37907328](https://pubmed.ncbi.nlm.nih.gov/37907328/)
- [SCAPE — NAR 2022 e66](https://academic.oup.com/nar/article/50/11/e66/6548409) · [GitHub LuChenLab/SCAPE](https://github.com/LuChenLab/SCAPE)
- [scMAPA — GigaScience 2022](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac033/6576244)
- [SAPAS — BMC Biology 2021](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01076-3) · [GitHub YY-TMU/SAPAS](https://github.com/YY-TMU/SAPAS)
- [scDaPars — Genome Research 2021](https://genome.cshlp.org/content/31/10/1856)
- [scAPAtrap — Briefings in Bioinformatics 2021](https://academic.oup.com/bib/article/22/4/bbaa273/5952304)
- [SCAPTURE — Genome Biology 2021 (PMC8353616)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8353616/)
- [scraps — bioRxiv 2022](https://www.biorxiv.org/content/10.1101/2022.08.22.504859v1) · [Hesselberth Lab publications](https://hesselberthlab.org/publications/)
- [scPAISO — bioRxiv 2025.08.20.669565](https://www.biorxiv.org/content/10.1101/2025.08.20.669565v1.full)
- [scDeepAPA — Briefings in Bioinformatics 2026](https://academic.oup.com/bib/article/27/3/bbag339/8715718)
- [spvAPA — Briefings in Bioinformatics 2025 (PMC11724721)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11724721/)
- [scPASU — STAR Protocols 2026](https://www.sciencedirect.com/science/article/pii/S2666166726001978)
- [PolyASite 2.0 — NAR 2020](https://academic.oup.com/nar/article/48/D1/D174/5588346)
- [PolyASite v3.0 — NAR 2025](https://academic.oup.com/nar/article/53/D1/D197/7893321) · [PubMed 39530237](https://pubmed.ncbi.nlm.nih.gov/39530237/)
- [polyApipe — GitHub MonashBioinformaticsPlatform/polyApipe](https://github.com/MonashBioinformaticsPlatform/polyApipe)
- [scAPAdb — NAR 2022](https://academic.oup.com/nar/article/50/D1/D365/6368523) · [scAPAmod (PMC9329739)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9329739)
- GitHub maintenance dates: GitHub REST API (`repos/{owner}/{repo}.pushed_at`), queried 2026-08-12.
