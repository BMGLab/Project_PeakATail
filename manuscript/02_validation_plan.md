# PeakATail Validation Strategy
## Establishing Correctness, Sensitivity, and Specificity of PAS/APA Calls for the Methods Manuscript

**Scope.** This document specifies the validation experiments that will support Claim 1 of the manuscript: PeakATail's PAS detections and APA switch calls are correct, sensitive, and specific, with performance quantified against ground truth and benchmarked against competing tools (Sierra, scAPA, scAPAtrap, SCAPE, SCAPTURE, polyApipe, MAAPER). Every experiment below runs identically on PeakATail and comparators. Two recent benchmark studies define the field's expectations and should be mined for harness code and datasets: the NAR benchmark of 3'-tag scRNA-seq PAS tools (NAR 2026, gkag490) and "Guidelines for alternative polyadenylation identification tools using single-cell and spatial transcriptomics data" (bioRxiv 2024.11.29.626111), which used four paired short-read/long-read datasets.

---

## 1. Ground-truth options, ranked by strength

No single ground truth is sufficient: atlases give **positions** but not sample-specific usage; bulk 3'-end-seq gives sample-matched **quantitative** truth but no cell resolution; matched long reads give per-cell truth but with 3'-end softclipping/coverage caveats; simulation gives **perfect** truth (including differential events) but imperfect realism. We use all four tiers, each answering a different question.

### Tier 1 (strongest): matched single-cell long-read data — same library, same cell barcodes

The decisive experiment: long reads from the *same 10x cDNA library* (or same sample) define, per cell barcode, which 3' end each molecule used. This validates both PAS positions and per-cell-type usage proportions — the exact quantity `ema switch diff` / `switch length` operate on.

| Dataset | What it is | Access | Use |
|---|---|---|---|
| **PacBio Kinnex single-cell PBMC, 10x 3'** (`DATA-Revio-Kinnex-PBMC-10x3p`; a 5' counterpart also exists) | Kinnex (commercialized MAS-ISO-seq; Al'Khafaji et al., *Nat Biotechnol* 42:582, 2024) HiFi reads on a 10x Chromium 3' PBMC library, CB/UMI-tagged | Public: `downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/` (verified to exist, Aug 2026) | **Primary validation dataset.** Run PeakATail on the short-read/Illumina side (or the pseudo-short reads); derive truth PAS + per-barcode usage from HiFi 3' ends. Doubles as the PBMC application dataset (Goal 3) — validation and application on one sample. |
| **FLT-seq / FLAMES data** (Tian et al., *Genome Biol* 2021; scmixology1/2, 5 cell lines) and **BLAZE data** (You et al., *Genome Biol* 2023) | Matched 10x short-read + ONT PromethION long-read single-cell | GEO/ENA, public | Cross-platform replication (ONT vs PacBio truth); cell-line mixtures give clean discrete "cell types" for switch validation |
| **SCOTCH benchmark datasets** (bioRxiv 2024.04.29.590597) | 10x + ONT (R9/R10) and 10x + PacBio, including PBMC samples | Public per paper's data availability | Additional ONT PBMC truth; reuse their accessions rather than generating new data |

**Truth extraction protocol.** Keep only long reads with (a) soft-clipped terminal poly(A) ≥ 20 nt (evidence of genuine polyadenylated 3' end, not internal fragment), (b) valid CB matching the short-read whitelist. Cluster long-read 3' termini within 25 nt into truth PAS; require ≥ 5 UMIs. Truth usage matrix = truth-PAS × cell-type UMI counts (cell types transferred from short-read clustering, shared barcodes).

*Caveat to state in the paper:* long reads under-sample long transcripts and have their own internal-priming rate; therefore we require poly(A)-tail evidence and treat long-read truth as high-precision but not exhaustive (sensitivity is computed against it; PAS absent from it are not automatically false).

### Tier 2: matched bulk 3'-end-sequencing

Orthogonal chemistry (no 10x priming artifacts), quantitative, but pseudo-bulk only.

| Dataset | Contents | Access | Use |
|---|---|---|---|
| **Derti et al. 2012 PolyA-seq atlas** (*Genome Res* 22:1173) | PAS maps for 24 matched tissues, 5 mammals; 158,533 known + 280,857 novel human sites | **GEO GSE30198**, SRA SRA039286; also a UCSC track (`polyASeqSites`) | Tissue-level truth for lung (Laughney pseudo-bulk) and blood-adjacent tissues |
| **Gruber et al. 2014** (*Nat Commun* 5:5465) | A-seq 3'-end-seq of **naive and activated CD4+ T cells, human and mouse** | GEO (per paper) | Direct quantitative truth for the T-cell activation contrast used as a positive control (Section 2) — pseudo-bulk PeakATail T-cell clusters against it |
| **Lianoglou et al. 2013 3'-seq** (*Genes Dev* 27:2380) | 3'-seq of human tissues and cell lines incl. immune lines | GEO (locate accession from paper) | Secondary tissue-level truth |
| *(Optional, strongest possible)* In-house QuantSeq-REV / 3'READS on an aliquot of the same PBMC sample | Sample-matched bulk 3' ends | To generate (~1 lane) | Only needed if reviewers demand sample-matched orthogonal chemistry; Tier 1 Kinnex largely covers this |

**Comparison:** pseudo-bulk PeakATail PAS (sum over cells) vs bulk sites: positional match (Section 4 tolerances) and usage-proportion correlation per gene.

### Tier 3: PAS atlases (positional truth, genome-wide)

| Atlas | Basis | Notes |
|---|---|---|
| **PolyASite 2.0** (Herrmann et al., *NAR* 48:D174, 2020; polyasite.unibas.ch) | Consolidated human/mouse/worm 3'-end-seq, uniform QC, **internal-priming sites flagged**, clustered with representative site per cluster | **Primary atlas.** The IP-flagged set is also a ready-made known-false-positive catalog (Section 3) |
| **PolyA_DB v3** (Wang et al., *NAR* 46:D315, 2018; polya-db.org/v3) and **v4** (3'READS+ human/mouse, long-read-validated) | 3'READS — chemically resistant to internal priming | Independent chemistry from much of PolyASite; use v4 where hg38 coordinates are needed |
| **PolyASite v3.0** (*NAR* 53:D197, 2025) and **scAPAdb** (*NAR* 50:D365) | Inferred from **scRNA-seq** | **Do not use as truth** — circular with the method under test. Cite only for annotation context |

**High-confidence truth set** = PolyASite 2.0 ∩ PolyA_DB (±25 nt), restricted to expressed genes in the evaluated sample. **Supported set** = union. Report precision against the union (a call in either atlas is not false) and sensitivity against the intersection (only jointly-supported sites are demanded). This two-set convention prevents atlas incompleteness from inflating FDR estimates and atlas noise from inflating sensitivity demands.

### Tier 4: read simulation with known PAS (the only truth for switch-level FDR)

Simulation is the only tier where differential events are known by construction, so it carries the switch-detection sensitivity/FDR-calibration experiments.

- **Base simulator: scReadSim** (Yan & Li, *Nat Commun* 14:7482, 2023; github.com/JSB-UCLA/scReadSim) — generates synthetic BAM/FASTQ mimicking a real 10x BAM at read-sequence and read-count level, with user-specified ground truth. Anchor it on our own STARsolo PBMC BAM so coverage shape, CB/UMI structure, and depth are realistic.
- **APA-specific layer (custom, on top):**
  1. Sample truth PAS per gene from PolyASite 2.0 (1–4 PAS/gene, matching the empirical PAS-per-gene distribution).
  2. Place read 5' ends upstream of each PAS using the **empirical peak-offset distribution** learned from real data (distance from read start to cleavage site in genes with a single unambiguous atlas PAS) — this is what makes 10x 3' simulation realistic, since reads sit in a smeared window upstream of the cleavage site, not at it.
  3. Assign per-cell-type PAS mixing proportions; create two synthetic "cell types" and **spike APA switches** into a known gene list with controlled effect sizes ΔPDUI ∈ {0.05, 0.1, 0.2, 0.3, 0.5} and cluster sizes ∈ {50, 100, 250, 500, 1000} cells.
  4. **Inject internal-priming decoys**: generate a controlled fraction (5–20%) of reads primed at genomic A-rich stretches (≥ 6 consecutive A, or ≥ 12 A in 18 nt) inside gene bodies — labeled false sites for specificity scoring.
  5. Sweep depth (10k–50k reads/cell) and add an ambient-RNA fraction.
- Report all simulator parameters + seeds; release the simulation code with the paper (reviewers will ask).

---

## 2. Known-biology positive controls

Literature-proven APA events with defined direction. Convention: **proximal shift = 3'UTR shortening = PDUI decreases** (PDUI = distal PAS usage fraction, `ema switch length --strategy classic`). Each event passes only if the switch is detected at q < 0.05 in the expected direction (`switch diff` + directional check via `switch length`), and is inspected with `switch geneview`.

### 2.1 Event checklist table

| # | Gene | Contrast (cell types) | Expected shift | Mechanism / note | Reference | Testable in our data? |
|---|---|---|---|---|---|---|
| P1 | *Global signature* (transcriptome-wide PDUI) | Proliferating (MKI67+ cycling: tumor cells, cycling T cells) vs quiescent counterparts | Global shortening: median ΔPDUI < 0; significantly more shortened than lengthened genes (sign test p < 0.01) | Proliferation-linked APA | Sandberg et al., *Science* 320:1643 (2008); Mayr & Bartel, *Cell* 138:673 (2009) | Yes — Laughney tumor vs normal epithelium; cycling vs non-cycling T cells |
| P2 | **HIP2/UBE2K** (Hip2) | Activated/proliferating T cells vs naive T cells | Proximal shift (loss of long 3'UTR carrying miR-21/miR-155 sites) | The validated exemplar gene of Sandberg 2008 (luciferase-confirmed) | Sandberg et al. 2008 | Yes — PBMC T-cell subsets; mouse ortholog if mouse data used |
| P3 | **IGHM** | Plasma cells vs naive/memory B cells | Plasma cells: proximal, intronic μs PAS (secreted IgM); B cells: distal μm PAS (membrane exons M1/M2) | CstF-64/ELL2-driven; the single best-established APA switch in biology; large effect size, unambiguous direction | Takagaki et al., *Cell* 87:941 (1996); Takagaki & Manley, *Mol Cell* 2:761 (1998); Martincic/Milcarek ELL2, *Nat Immunol* 2009 | Yes — PBMC B/plasma cells; Laughney tumor-infiltrating plasma cells. **Headline positive control** |
| P4 | **CD47** | Activated immune cells / cells with high surface CD47 vs resting; tumor vs normal | Two robust annotated PAS (short-UTR SU and long-UTR LU isoforms) must both be detected; differential SU/LU usage across clusters, proximal shift in proliferative clusters | Long-UTR isoform scaffolds HuR/SET → plasma-membrane CD47 ("don't eat me"); TAM-relevant | Berkovits & Mayr, *Nature* 522:363 (2015) | Yes — detection of both PAS is the hard criterion; direction exploratory |
| P5 | **NUDT21/CFIm25 target set: CCND1, DICER1, TIMP2, PAK1, MECP2** | Clusters/cells with low vs high NUDT21 expression (e.g., tumor vs normal epithelium) | Low NUDT21 → proximal shift of targets; test as negative correlation between cluster-level NUDT21 expression and target PDUI across clusters (Spearman, expect ρ < 0) | CFIm25 loss shortens these targets (qRT-PCR validated for CCND1/DICER1/TIMP2; PAK1 in GBM; MECP2 dosage) | Masamha et al., *Nature* 510:412 (2014); Chu et al., *Oncogene* 2019 (PAK1); Gennarino et al., *eLife* 2015 (MECP2) | Yes — regulator-target coherence test, a stronger-than-anecdote control |
| P6 | *Innate-immunity shortening signature* | Activated/inflammatory monocytes-macrophages (LPS/infection-like state, M1-polarized TAMs) vs resting monocytes | Global proximal shift in immune-response genes; gene set = significantly shortened genes from Pai et al. supplementary tables (use their list, do not curate ad hoc); GSEA-style enrichment of that set among our shortened genes | Conserved 3'UTR shortening upon macrophage infection; miRNA-escape mechanism | Pai et al., *PLoS Genet* 12:e1006338 (2016); Frontiers Immunol 14:1182525 (2023) for M1/M2 polarization | Yes — TAM application dataset (Goal 3) and PBMC monocytes; links validation to the TAM biology chapter |
| P7 | *T-cell activation, quantitative* | Naive vs activated CD4+ T cells | Per-gene ΔPDUI from PeakATail correlates with ΔPAS-usage measured by A-seq (expect r > 0.5 on shared genes) | Sample-type-matched bulk 3'-end-seq truth (Tier 2) | Gruber et al., *Nat Commun* 5:5465 (2014) — human and mouse naive/activated CD4+ T cells | Yes — PBMC |
| P8 | **PAX3** (mouse) | Diaphragm quiescent satellite cells vs limb satellite cells / activated SCs | Diaphragm QSCs: proximal PAS (short 3'UTR escaping miR-206); limb: distal PAS retaining miR-206 site | miRNA-escape by APA in muscle stem cells | Boutet et al., *Cell Stem Cell* 10:327 (2012) | **Only if a mouse muscle scRNA-seq dataset is added** (e.g., Tabula Muris limb muscle). Optional; include as external-dataset control or drop |

### 2.2 Scoring the checklist

- Primary endpoint: **fraction of applicable checklist events recovered with correct direction at q < 0.05**. Pre-register the applicable set (P1–P7 for human immune/lung data) before the final run; report the same checklist for every benchmarked competitor tool — this becomes a main-text figure/table.
- For each event show a `switch geneview` panel (per-cluster PAS read distributions) plus, where Tier 1 data covers the gene, the long-read 3' ends overlaid — biology + orthogonal evidence in one figure.
- P1/P6 are signature-level (protect against single-gene dropout); P2–P5, P7 are gene-level.

---

## 3. Negative controls and specificity

### 3.1 Internal-priming (IP) artifacts

The dominant 10x 3' failure mode: oligo(dT) priming on genomic/intronic A-stretches creates fake PAS.

- **Sequence-based FP labeling:** flag any called PAS with ≥ 6 consecutive A or ≥ 12 A within the 18 nt genomic window immediately downstream of the inferred cleavage site (standard PolyA_DB/PolyASite-style criterion), *and* no atlas support and no long-read poly(A)-tail support. Report **IP rate = flagged calls / all calls**, before and after `ema` ip-filter mode — the before/after delta quantifies the filter's value (ablation figure).
- **Known-false catalog:** PolyASite 2.0's flagged putative internal-priming sites = ready-made decoy set. Metric: fraction of decoys that PeakATail calls as PAS (target ≈ 0 after filtering).
- **Poly(A)-signal enrichment (positive sequence diagnostic):** fraction of called PAS with AAUAAA or one of the 12–18 canonical variants in the −40…−10 nt window. Atlas-supported true sites show 80–90%; a call set whose novel sites drop far below its known sites indicates FP contamination. Also plot the single-nucleotide frequency profile ±100 nt around called cleavage sites: true PAS show the canonical A-rich(signal)/U-rich(DSE) profile; mispriming shows a downstream genomic A-spike.
- **Simulation decoys** (Tier 4, step 4): known-label IP reads → exact FP attribution.

### 3.2 Single-PAS housekeeping genes (known-negative gene set)

- Define empirically, not from folklore: genes with **exactly one** PAS cluster in *both* PolyASite 2.0 and PolyA_DB (and one 3' end in long-read truth where covered), expressed in ≥ 30% of cells. Expect n in the hundreds.
- Metrics: (a) **spurious-multiplicity rate** = fraction of these genes where PeakATail calls ≥ 2 PAS; (b) **spurious-switch rate** = fraction appearing as significant in `switch diff` at q < 0.05 (should be ≈ 0 since one PAS cannot switch).
- Report the same rates for competitor tools — this is where peak-caller over-segmentation shows up.

### 3.3 Statistical calibration by label permutation

- **Cluster-label permutation:** shuffle cell→cluster assignments (preserving cluster sizes), re-run `switch diff`, ≥ 100 permutations. Expected: uniform p-values (QQ plot, genomic inflation λ ≈ 1) and empirical FDR at q < 0.05 near zero. This tests the Fisher/NB machinery + multiple-testing stack end-to-end.
- **Within-cluster split:** randomly bisect one large homogeneous cluster into two pseudo-clusters; any significant switch is a false positive. Closer to the real null than full permutation (preserves within-cluster structure).
- **Replicate concordance null/positive pair:** same annotated cell type across different patients (Laughney) or PBMC donors — expect *few* switches between donors' matched cell types (specificity) but *reproducible* PDUI per gene (per-gene PDUI correlation across donors, r ≥ 0.8 for genes with ≥ 20 UMIs).
- **Ambient-RNA stress test:** repeat permutations on clusters with divergent total counts to confirm depth differences alone do not create switch calls.

### 3.4 Technical negative controls

- **Strand sanity:** fraction of called PAS antisense to their assigned gene (should be < 1–2%; PeakATail calls per-strand, so this audits the gene-assignment step, `--max-gene-distance`).
- **Parameter-sweep stability via `ema reannotate`:** OFAT sweep over gene-distance/filter/clustering knobs from one base peak-call; positive-control events (Section 2) must survive across the reasonable parameter range — guards against a tuned-to-the-demo method.

---

## 4. Metric definitions and pass/fail criteria

### 4.1 Definitions (fixed before the final run)

**Matching rule.** A called PAS matches a truth PAS if within **±50 nt** (primary tolerance; report a sweep at ±20/±50/±100/±150 as robustness — 10x 3' peak smearing makes single-nt matching meaningless, and the field's benchmarks use tolerances in this range). One-to-one greedy matching by distance; a truth site can validate only one call (extra calls at the same site count as FPs → penalizes over-segmentation).

**Evaluation universe.** Restrict to **expressed genes**: ≥ 20 pseudo-bulk UMIs in the evaluated dataset (sensitivity must not be punished for unexpressed truth sites); state this in Methods.

| Metric | Definition | Truth tier |
|---|---|---|
| Sensitivity (recall) | TP / (TP + FN); FN = truth PAS in expressed genes with no matched call | Tiers 1, 2, 3(∩), 4 |
| Precision | TP / (TP + FP); FP = call with no match in the *union* truth set | Tiers 1, 3(∪), 4 |
| Empirical FDR (detection) | 1 − precision | same |
| F1 | harmonic mean; headline benchmark number + PR curve over peak-score threshold | same |
| Positional error | median |called − truth| distance for TPs | Tiers 1, 2, 4 |
| Quantification accuracy | Pearson/Spearman between called per-PAS usage proportions and truth proportions (per gene, ≥ 2 PAS, ≥ 20 UMIs); plus per-cell-type PDUI correlation | Tiers 1, 2, 4 |
| Switch sensitivity | fraction of spiked switches detected at q < 0.05, as a function of ΔPDUI × cluster size (power surface figure) | Tier 4 |
| Switch FDR calibration | empirical FDR among switch calls vs nominal q, from spiked truth (Tier 4) and permutations (Section 3.3); calibration plot nominal-vs-empirical | Tier 4 + permutation |
| Direction accuracy | fraction of detected true switches with correct shortening/lengthening sign | Tiers 1, 4; checklist |
| IP rate / decoy recall / spurious-multiplicity / spurious-switch | as defined in Section 3 | Section 3 sets |

### 4.2 Pass/fail table (pre-registered acceptance criteria)

| # | Experiment | Metric | Pass threshold | Rationale |
|---|---|---|---|---|
| V1 | PBMC Kinnex 10x 3': PAS detection vs long-read truth | Precision / Sensitivity @ ±50 nt | ≥ 0.80 / ≥ 0.70 (truth sites with ≥ 10 long-read UMIs) | Competitive with best reported tools (SCAPE/scAPAtrap range) |
| V2 | Same, quantification | Per-gene usage-proportion correlation | median r ≥ 0.80; per-cell-type PDUI r ≥ 0.75 | Matches SCAPE-class quantification quality |
| V3 | Atlas comparison (all datasets) | Precision vs PolyASite∪PolyA_DB; sensitivity vs ∩ | ≥ 0.75 / ≥ 0.65 | Atlas incompleteness bounded by two-set convention |
| V4 | Bulk 3'-end-seq (Derti lung; Gruber T cells) | Pseudo-bulk usage correlation | r ≥ 0.7 on shared genes | Orthogonal-chemistry quantitative agreement |
| V5 | Simulation: detection | F1 @ ±50 nt across depth sweep | ≥ 0.85 at ≥ 20k reads/cell; graceful degradation curve reported | Controlled-truth ceiling performance |
| V6 | Simulation: switch power | Sensitivity for ΔPDUI ≥ 0.2, ≥ 250 cells/cluster | ≥ 0.80 | Effect sizes typical of real biology (IGHM-class events ≫ 0.2) |
| V7 | Simulation + permutation: FDR calibration | Empirical FDR at nominal q = 0.05 | ≤ 0.075 (≤ 1.5× nominal); QQ λ ∈ [0.9, 1.1] | Statistical validity of switch-diff |
| V8 | IP artifacts | Decoy/IP-flagged retention after ip-filter | ≤ 5% of injected/flagged decoys called; IP rate delta reported | Specificity against the dominant artifact |
| V9 | Single-PAS genes | Spurious-multiplicity rate; spurious-switch rate | ≤ 10%; ≤ 1% | Over-segmentation control |
| V10 | Positive-control checklist (P1–P7) | Events recovered, correct direction, q < 0.05 | ≥ 80% of applicable events; IGHM (P3) mandatory | Biological correctness |
| V11 | Donor/patient replicates | Cross-replicate PDUI correlation; cross-replicate switch rate | r ≥ 0.8; matched-cell-type switches < 5% of between-cell-type switches | Reproducibility + specificity |
| V12 | Benchmark ranking | F1 (V1, V5) and checklist recovery (V10) vs Sierra, scAPA, scAPAtrap, SCAPE, SCAPTURE, polyApipe | PeakATail in top 2 on ≥ 2 of 3 axes; any axis where it is not, discussed honestly | Q1-journal comparative bar |

**Failure handling (pre-registered):** a failed criterion triggers (i) diagnosis via the `ema reannotate` parameter sweep, (ii) one documented parameter revision applied uniformly across all experiments, (iii) re-evaluation; no per-experiment tuning. Thresholds, tolerances, and the checklist are frozen before the final benchmark run and stated as such in Methods.

### 4.3 Figure plan mapped to claims

1. **Fig. V-A:** Schema of the four-tier truth design.
2. **Fig. V-B:** PBMC Kinnex — PR curves + usage-correlation scatter, PeakATail vs competitors (V1/V2/V12).
3. **Fig. V-C:** Simulation — F1 vs depth; switch-power surface (ΔPDUI × cells); FDR calibration plot (V5–V7).
4. **Fig. V-D:** Specificity panel — nucleotide profiles around called PAS, IP rates before/after filter, single-PAS gene rates, permutation QQ (V8, V9, V7).
5. **Fig. V-E:** Positive-control gallery — IGHM geneview with long-read overlay, HIP2/UBE2K, CD47, NUDT21-target coherence, global-shortening signatures (V10).

---

## Key verified references and resources

- PolyASite 2.0: Herrmann et al., *NAR* 48:D174 (2020) — https://polyasite.unibas.ch (note: PolyASite v3.0, *NAR* 53:D197, is scRNA-seq-derived → excluded from truth for circularity)
- PolyA_DB v3: Wang et al., *NAR* 46:D315 (2018) — http://www.polya-db.org/v3; PolyA_DB v4 (3'READS+, long-read-validated)
- PolyA-seq atlas: Derti et al., *Genome Res* 22:1173 (2012) — GEO **GSE30198** / SRA SRA039286
- Gruber et al., *Nat Commun* 5:5465 (2014) — A-seq, naive/activated human + mouse CD4+ T cells
- Kinnex/MAS-ISO-seq: Al'Khafaji et al., *Nat Biotechnol* 42:582 (2024); public PBMC 10x-3' Kinnex data: https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/
- FLT-seq/FLAMES: Tian et al., *Genome Biol* 22:310 (2021); BLAZE: You et al., *Genome Biol* 24 (2023); SCOTCH: bioRxiv 2024.04.29.590597
- scReadSim: Yan & Li, *Nat Commun* 14:7482 (2023) — https://github.com/JSB-UCLA/scReadSim
- Benchmarks to align with: *NAR* 54:gkag490 (2026); bioRxiv 2024.11.29.626111
- Biology: Sandberg et al., *Science* 320:1643 (2008); Mayr & Bartel, *Cell* 138:673 (2009); Takagaki et al., *Cell* 87:941 (1996); Takagaki & Manley, *Mol Cell* 2:761 (1998); Berkovits & Mayr, *Nature* 522:363 (2015); Masamha et al., *Nature* 510:412 (2014); Gennarino et al., *eLife* (2015); Boutet et al., *Cell Stem Cell* 10:327 (2012); Pai et al., *PLoS Genet* 12:e1006338 (2016)
- Application data: Laughney et al., *Nat Med* (2020) — GEO **GSE123904** (verified)
