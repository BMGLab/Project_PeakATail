<!--
sec_front.md — title, structured abstract (Genome Biology Method track), key points.
Drafted 2026-09-02 under 21_paper_architecture.md (claim ladder, must-not-claim list) and the
binding language rules. Every number is quoted from a verified evidence document; the source map
is in the PROVENANCE comment at the end of this file. No number here is new.
-->

# PeakATail: precision-first poly(A)-site calling and calibrated alternative polyadenylation analysis in single-cell RNA-seq

<!-- TITLE REVISED 2026-09-02 per PI: the interrogative title was too long and too informal.
     Previous: "How much of a single-cell poly(A)-site call set can be trusted? Precision-first calling,
     calibrated switch testing and pre-registered evaluation in PeakATail".
     Runner-up kept for the record: "Precision-first poly(A)-site calling and calibrated APA switch
     testing in single cells with PeakATail". Pre-registered evaluation moved to the abstract/key points. -->


<!--
RUNNER-UP TITLE (21 §2 rank 2): "Precision-first poly(A)-site calling and calibrated cell-type APA
switch testing from 3′ single-cell RNA-seq". Safe and descriptive, but it does not signal the
benchmark/pre-registration half of the paper's value and omits the tool name; 21 rates it
"forgettable". The chosen title is 21's rank-1 candidate verbatim — the only one where the failed
pre-registered gate is an asset rather than an embarrassment, and it makes no accuracy-leadership
promise.
-->

Amir Amiri Tabat^1^, Yasin Kaymaz^2,\*^

^1^ Department of Computer Engineering, Faculty of Computer and Informatics Sciences, Ege University, İzmir, Türkiye
^2^ Department of Bioengineering, Faculty of Engineering, Ege University, İzmir, Türkiye

\* Corresponding author: Yasin Kaymaz, yasin.kaymaz@ege.edu.tr (ORCID: 0000-0002-9725-7536)

<!-- AUTHORS SET 2026-09-02 per PI instruction (verbatim): Amir Amiri Tabat (Ege University Faculty of
     Computer and Informatics Sciences, Department of Computer Engineering); Yasin Kaymaz* (Ege University
     Faculty of Engineering, Bioengineering Department); * Corresponding Author. Corresponding e-mail and Y.K. ORCID supplied by the PI 2026-09-02.
     A.A.T. ORCID: optional, not yet supplied — [[PLACEHOLDER: A.A.T. ORCID if available]]. -->

## Abstract

<!-- ABSTRACT-BEGIN (244 words incl. the three section labels, `wc -w` after stripping bold markers, LC_ALL=C, 2026-09-02; limit ≤ 250) -->

**Background:** Alternative polyadenylation (APA) remodels 3′UTRs across cell types, yet single-cell poly(A)-site (PAS) calls are rarely scored against independent truth and differential tests are rarely calibrated.

**Results:** PeakATail seeds PAS from non-templated poly(A)-tail reads, tiers them by molecule-counted clip support, filters internal priming, and adds a calibrated, replication-filtered switch test. Its pre-registered default reaches atlas-agreement precision (curated PolyASite 2.0, 100 bp) of 0.706 (PBMC), 0.745/0.757 (two mice) and, un-tuned, 0.828 (second human donor)—gate P ≥ 0.50 passed on all four libraries, gene-body-shuffled nulls 0.013–0.022—the most atlas-concordant de novo call set on both benchmark datasets. Detected-gene recall (R_det 0.110–0.208) is below polyApipe's at this conservative point, though at matched call count PeakATail leads on recall at every budget tested and precision above 20,000 calls. 76.5% of default PBMC sites match a poly(A)-verified, donor-mismatched Kinnex long-read 3′ end within 25 bp (≥5 records). Only one of six switch-test configurations—Fisher on cells, no marker pre-selection—controls the false-discovery rate (3.0% label-shuffled p < 0.05; shipped defaults 13.0–24.7%). In 15 libraries from 12 lung-cancer patients, 15,942 of 2,128,711 switch hypotheses replicated across patients—none in ten label-shuffle nulls (empirical p ≤ 0.091, the ten-permutation floor); in testis, 31.4%/30.5% of depth-guarded genes shorten 3′UTRs monotonically along spermatogenesis in both mice (z = 10.5/6.4). A pre-registered "trusted novel PAS" definition reached 48.5% against a 70% long-read target—a reported negative result. A 10k-cell library runs in 12.53 GB peak memory.

**Conclusions:** Pre-registered gates, calibrated tests and null-disciplined replication make single-cell APA auditable—failures included. Code: github.com/BMGLab/PeakATail.

<!-- ABSTRACT-END -->

**Keywords:** alternative polyadenylation; poly(A) sites; single-cell RNA-seq; 3′ end sequencing; benchmarking; false-discovery-rate calibration; pre-registration; replication

## Key points

<!-- One bullet per load-bearing claim (21 §5); the negatives appear as findings per 21 §6 rule (iv). -->

- **Evidence type, not peak-calling sophistication, separates single-cell PAS callers.** Re-seeding our own coverage-only caller on non-templated poly(A)-tail reads moved atlas-agreement precision at 100 bp from 0.118 to 0.7062 on the same PBMC BAM, and the six-tool panel ordering tracks the evidence available to each caller, not its algorithm class.

- **The pre-registered precision default** (clip-supported, ≥2 distinct molecules, internal-priming-filtered—a reliability choice, not the F1 optimum) **is the most atlas-concordant de novo call set on both benchmark datasets** (atlas-agreement P@100 0.7062 PBMC, 0.7450/0.7572 testis mice) **and passes its gate un-tuned on a second human donor** (pbmc4k 0.8279; 84.2% of donor-2 default calls reproduce in donor 1 within 100 bp). Among de novo tools at matched call count, PeakATail leads polyApipe on recall at every budget tested and on precision at every budget above 20,000 calls.

- **Only one of six differential-test configurations controls the false-discovery rate** (Fisher, cells mode, marker pre-selection off: 3.0% of label-shuffled tests p < 0.05); the defaults we originally shipped are anti-conservative (13.0–24.7%). We report this failure of our own tool as a finding, together with the label-shuffle calibration harness that exposed it.

- **Switches are reported only when they replicate across patients:** on the precision-first PAS universe, 15,942 of 2,128,711 tested (cell-type pair, PAS) hypotheses replicate across 15 libraries from 12 patients, with none replicated in ten label-shuffle nulls (empirical p ≤ 0.091, the ten-permutation floor); the testis control replicates across mice (6,070/9,480/11,219 same-direction PAS per stage pair, 0 in all 15 null pairings, per-gene effect ρ 0.641, n = 923).

- **The failures are results.** The pre-registered "trusted novel PAS" definition missed its 70% long-read target (48.5%), and no site in this paper is called trusted-novel; PAS-profile clustering recovers cell types no better than gene-level totals (AMI 0.708 vs 0.698), so that claim is withdrawn; a 16-gene literature panel of spermatogenesis 3′UTR shortening does not reproduce; and an atlas-trained scoring model that gained 16 precision points on the atlas lost 8–11% on every long-read axis—curated-atlas benchmarking rewards atlas-shaped priors. The pre-registration timeline that makes these outcomes auditable is Supplementary Fig. S11.

<!--
PROVENANCE (number → verified source; nothing below is quotable text, it is an audit map):
- 0.706/0.7062 (PBMC default P@100, arm: precision default), 0.745/0.7450, 0.757/0.7572 (mice),
  R_det 0.1754/0.2048/0.2080, gene-body-shuffled nulls 0.013–0.022 (0.0127–0.0217; pbmc4k seeds
  0.0216/0.0222/0.0217 round to 0.022), gate P ≥ 0.50 PASS: 19_final_gate_v2.md §1–2 via 21 L1/C1;
  pbmc4k nulls from 26 §"pbmc4k" table.
- 0.828/0.8279, R_det 0.1104, "gate holds on four independent libraries", 84.2% cross-donor
  concordance within 100 bp (donor 2 → donor 1): 26_second_donor_preregistration.md.
- R_det range 0.110–0.208 = pbmc4k 0.1104 … mouse 2 0.2080 (26 + 19).
- matched-call-count sentence: licensed wording of 25_competitive_position.md §8.1 ("recall at every
  budget tested", "precision at every budget above 20,000 calls"); single-point recall comparison
  never quoted without it (25 §8.2).
- 76.5% (0.7647 at ≥5 support, 25 bp, x3p truth; arm: full PBMC precision default), donor-mismatched:
  16_trusted_novel_kinnex.md §v2 via 21 L2/C2. "≥5 records" wording per 21 §7 item 8 (the truth-set
  support threshold is an alignment-record count, not a de-duplicated UMI count — do not write bare
  "≥5 UMI" without that caveat).
- calibration: 3.0% null p < 0.05 (Fisher cells-mode, marker pre-selection off), shipped/candidate
  arms 13.0–24.7% (13.0 cells-Fisher+top-200 markers, 20.3 reads-Fisher, 24.7 nb_pairwise):
  14_switch_calibration_v2.md via 21 L3/C3 and §6.
- Laughney: 15,942 / 2,128,711 replicate (15,212 over the |Δproportion| ≥ 0.1 floor — floored count
  not in the abstract for word budget, carried in Results), replication primary 15 libraries from
  12 patients (cohort: 17 libraries from 14 patients), 0 in each of 10 label-shuffle nulls, wording
  "none … (empirical p ≤ 0.091, the ten-permutation floor)" per 21 §7 item 2: 20_stage3_replication.md
  §v2 (verified SOUND).
- spermatogenesis v2 (the record): 31.4%/30.5% monotone shortening, z = 10.5/6.4, replication
  6,070/9,480/11,219, 0 in all 15 null pairings, ρ 0.641 (n = 923): 18_spermatogenesis_final.md
  v1→v2 table + "The claim, in the exact form that survives on v2".
- trusted-novel negative: 48.5% (0.4853) vs 70% target: 16 §v2 via 21 C11.
- compute: 12.53 GB peak RSS (PBMC default arm) — quoted without wall time so no concurrency
  disclosure is triggered (21 §7 item 5): 19_final_gate_v2.md §3.
- 0.118 → 0.7062 same-BAM before/after: 09_headtohead_results.md + 19 §1 via 21 C5.
- AMI 0.708 vs 0.698 clustering ablation, 16-gene literature panel: 21 §6 / 18 v2 negatives.
- atlas-trained scoring model "+16 atlas precision points, −8–11% on every long-read axis":
  23_algorithm_roadmap.md §4 item 4 (the atlas-vs-long-read decomposition).
- Supplementary figure numbering frozen per figures/FIGURES_MANIFEST.md: S11 = gatehistory.
Deliberately absent: any F1 claim; any "trusted novel" positive claim; clustering novelty; named
switch genes (issue #99); wall times (concurrency rule); 1.152% clip rate (corrected to 0.573%);
0.9986; licence statement (no LICENSE file yet — PI decision, per 21 §11 correction 5).
-->
