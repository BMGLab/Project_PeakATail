# NUMBERS_AUDIT — traceability audit of the manuscript text (2026-08-21)

**Scope audited:** `manuscript/01_outline_and_journals.md` (outline v0.4) and all five
`manuscript/figures/*.caption.md` (`final_benchmark`, `trade_reproducibility`,
`trusted_novel_funnel`, `spermatogenesis_final`, `laughney_switches`).

**Reference set (the only documents quotable as results):** `13_reliability_positioning.md`,
`14_switch_calibration_v2.md`, `15_final_gate.md` (v1 record), `16_trusted_novel_kinnex.md`
(incl. §v2), `18_spermatogenesis_final.md`, `19_final_gate_v2.md`,
`20_stage3_replication.md` — plus, where those documents cite one, the TSV itself.

**Method.** Every numeric assertion in scope was classified as
`traced` (verbatim or an exact rounding of a value in a reference doc or in a TSV that doc cites),
`stale` (a v1 value where a verified v2 value exists),
`inconsistent` (the same quantity quoted with two different values inside the manuscript),
`repo-source` (traced to a repo write-up **outside** the reference set — `05_figure_index.md`,
`07_curated_benchmark_report.md`, `09_headtohead_results.md`, `10_caller_fix_plan.md`,
`11_verifier_corrections.md`, `12_stage2_gate.md` — not re-verified by this audit), or
`untraceable`.
TSVs were read directly and, where a value was definitional rather than quoted, recomputed here from
primary data (BED/TSV) with the computation stated. `LC_ALL=C` throughout. The Stage-3 v2 tree
`/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3/` was **not touched**; the v1 tree
`stage3_laughney_v2/` was read only through the figure TSVs and `20`.

**Verdict:** **829 numeric values examined — 795 result numbers plus 34 editorial journal-format figures.** No fabricated number was found and no result number
was untraceable. Fixed: one stale v1 value (the "+0.015 / +0.010" F1 gap), three internal collisions
(the 72.1 / 72.2 / 72.5 single-molecule share, the 65–180× / 65–179× enrichment, and one wrong source
attribution), one forbidden-phrasing hit, two stale cross-document provenance statements, one wrong
self-reported word count and one figure-layout overlap. Listed for the PI: 4 further collisions,
13 `[final run pending]` placeholders and 6 claim-level conflicts with the disclosure sections of
`16` / `18` / `20`.

---

## 1. Summary by verdict

| verdict | values | notes |
|---|---:|---|
| traced (reference doc or a TSV it cites) | 705 | includes values recomputed here from the primary BED/TSV to confirm a definition |
| untraceable | 1 | row A14, the abstract's self-reported "199 words" — **fixed to 197** |
| repo-source (05 / 07 / 09 / 10 / 11 / 12) | 53 | legitimate per the outline's own number policy; outside this audit's re-verification mandate |
| mixed traced / repo-source in one row | 3 | rows A57, A64, A88 |
| inconsistent (collision inside the manuscript) | 23 | rows A16, A48, A70, A75, A93, B6, B48 — 3 fixed, 4 left for the PI |
| stale v1 | 7 | A54 (the "+0.015 / +0.010" F1 gap), B27 (16 cited as v1-only), B40 (a quote of 18 that 18 no longer contains) — all **fixed** |
| forbidden phrasing (values themselves correct) | 2 | B3 — **fixed** |
| contradicted by its evidence document | 1 | A72 ("6/6 literature genes") — disclosure added, claim left for the PI |
| editorial (journal IFs, word limits, APCs) | 34 | row A94; not result numbers |
| **total examined** | **829** | 795 of them result numbers |

**Untraceable:** exactly one — the abstract's self-reported length "199 words by `wc -w`" (measured
197, or 198 with the `>` marker). **Fixed.** No result number in scope was untraceable.

---

## 2. Table A — `manuscript/01_outline_and_journals.md`

`n` = number of individual numeric values covered by the row.

| # | quantity (as written) | location | source checked | n | verdict |
|---|---|---|---|---:|---|
| A1 | PBMC P@100 0.7167 → 0.7062; mice 0.7414 → 0.7450, 0.7554 → 0.7572 | header ¶ | `19` §1 table; `VERIFIED_v2.md` §1 | 6 | traced |
| A2 | every prediction hit to ≤ 0.001 | header ¶ | `VERIFIED_v2.md` "Movement vs the PR predictions" (residuals ≤ 0.00006) | 1 | traced |
| A3 | trusted-novel 49% vs 70% target on v2 | header ¶ | `16` §v2 (0.4853); `13` §3 (target 0.70) | 2 | traced |
| A4 | shipped-caller baseline P 0.118 / F1 0.134 (PBMC) | number policy ¶ | `15` §3 | 2 | traced |
| A5 | atlas-agreement P@100 0.706–0.757 vs 0.50 | journals, BiB row | `19` §1 | 3 | traced |
| A6 | nulls 0.013–0.022 | abstract | `VERIFIED_v2.md` §1 (0.0127–0.0217) | 2 | traced |
| A7 | P 0.706 PBMC, 0.745 / 0.757 mice; gate ≥ 0.50 | abstract | `19` §1 | 4 | traced |
| A8 | R_det 0.18–0.21 | abstract | `19` §1 (0.1754 / 0.2048 / 0.2080) | 2 | traced |
| A9 | 76% of PBMC sites within 25 bp of a Kinnex 3′ end at ≥5 UMI | abstract | `16` §v2 / `validation.tsv` (0.764659) | 3 | traced |
| A10 | trusted-novel definition missed 70% target (49%) | abstract | `16` §v2 | 2 | traced |
| A11 | Kinnex concordance 0.7647; trusted-novel 0.4853 vs 0.70 | abstract note | `validate_incl_genebodies/validation.tsv` | 3 | traced |
| A12 | 3.0% null p < 0.05 ("calibrated") | abstract note | `14` arm B0 | 1 | traced |
| A13 | full-atlas recall 0.0889 / 0.0977 / 0.0995 | abstract note | `VERIFIED_v2.md` §1 (R_full column) | 3 | traced |
| A14 | "(199 words by `wc -w`)" | abstract note | measured: 198 with the `>` prefix, **197** without | 1 | **untraceable → fixed to 197** |
| A15 | Fig 2 panel-d values 0.0558 / 0.3520 / 0.7062 | Fig 2 block | `VERIFIED_v2.md` (PBMC IP tier-2 0.0558; IP tier-1 0.3520; default 0.706195) | 3 | traced |
| A16 | 72.1% of tier-1 sites single-molecule | Fig 2 block, R1, R2 | `19` §4 (verifier recount); recomputed here: 1 − 62,110 / 222,955 = 0.7214 on the **no-IP** arm; the IP arm gives 0.7223 and `15` §4's v1 IP arm 0.7253 | 1 | **inconsistent (72.1 / 72.2 / 72.5) → fixed: 72.1 kept, definition stated everywhere** |
| A17 | default P 0.7062 (n 46,524), nulls 0.0217/0.0217/0.0217 | Fig 2A | `VERIFIED_v2.md` §1 | 5 | traced |
| A18 | mouse 1 0.7450 (n 26,255), nulls 0.0149/0.0129/0.0139 | Fig 2A | `VERIFIED_v2.md` §1 | 5 | traced |
| A19 | mouse 2 0.7572 (n 26,526), nulls 0.0137/0.0127/0.0150 | Fig 2A | `VERIFIED_v2.md` §1 | 5 | traced |
| A20 | null 33–55× below | Fig 2A | `19` §1 | 2 | traced |
| A21 | tier-1 ≥1 mol: 0.3520 (n 167,565) IP / 0.3032 (n 222,955) no-IP; mice 0.5686 / 0.5880 | Fig 2A | `VERIFIED_v2.md` §2 + remaining-outputs table | 6 | traced |
| A22 | tier-2 PBMC 0.0568 / 0.0558; mice 0.1191 / 0.1142 | Fig 2A | `VERIFIED_v2.md` remaining outputs | 4 | traced |
| A23 | both tiers PBMC 0.1932 / 0.2044; mice 0.4488 / 0.4539 | Fig 2A | `VERIFIED_v2.md` §2 + remaining outputs | 4 | traced |
| A24 | gate lines P ≥ 0.38, F1_det > 0.261; F1_det 0.3001 / 0.3046 | Fig 2A / 2F | `VERIFIED_v2.md` §2; `10` §5 for the gate | 4 | traced |
| A25 | `pas_tier1_ge2mol_noIP_POSTHOC.bed` n 62,110, P 0.591 | Fig 2A, positioning (f) | `VERIFIED_v2.md` (62,110 / 0.5907) | 2 | traced |
| A26 | denominators 285,136 / 126,686 | Fig 2B, R2, Methods 10 | `VERIFIED_v2.md` §1 | 2 | traced |
| A27 | R_det 0.1754 / 0.2048 / 0.2080; F1_det 0.2811 / 0.3213 / 0.3264 | Fig 2B, R2 | `VERIFIED_v2.md` §1 | 6 | traced |
| A28 | PBMC default across windows P 0.521 / 0.639 / 0.674 / 0.706 / 0.737 | Fig 2B | `score_pbmc_final_v2_ipfilt__PRESPEC_precision_default.tsv` (0.520871 … 0.737447) | 5 | traced |
| A29 | PBMC default across windows R_det 0.084 / 0.105 / 0.135 / 0.175 / 0.220 | Fig 2B | same TSV (0.083764 … 0.219818) | 5 | traced |
| A30 | tier-1 ≥1 mol R_det: 0.2685 / 0.2970; mice 0.2806 / 0.2826 | Fig 2B | `VERIFIED_v2.md` | 4 | traced |
| A31 | full-atlas recall 0.0889 (569,005) / 0.0977 / 0.0995 | Fig 2B | `VERIFIED_v2.md` §1 | 4 | traced |
| A32 | "~29% of detected-gene atlas sites carry any clip read" | Fig 2B, Limitations | `10_caller_fix_plan.md` §"The ceiling" (28.77% of 285,220) | 1 | repo-source |
| A33 | `pasbed.bed` drops 0.45% PBMC / 3.7% mouse 1 | Fig 2B, Methods 6 | `15` §6 item 3 | 2 | traced |
| A34 | Kinnex default 0.7647 / 0.5584 / 0.2807 / 0.0976 | Fig 2C, R2 | `validation.tsv` `CAL_default_output_all` | 4 | traced |
| A35 | atlas-known hexpass n 28,453: 0.8941 / 0.7072 / 0.3885 / 0.1415 | Fig 2C | `validation.tsv` `CAL_atlas_known_hexpass` | 5 | traced |
| A36 | trusted-novel n 7,259: 0.4853 / 0.2393 / 0.0657 / 0.0118 | Fig 2C, Fig 5D–F, R5 | `validation.tsv` `08_trusted_novel` | 5 | traced |
| A37 | gene-body-shuffled null ≈ 0.0074 at ≥5 UMI (10 seeds) | Fig 2C, Fig 5D–F | `validation.tsv` null_mean 0.007411 | 2 | traced |
| A38 | GEM-X 0.5154; pooled 0.5816 (0.6688 @50 bp, 0.7143 @100 bp) | Fig 2C, Fig 5D–F | `verifier_crosscheck/RESULTS.txt` | 4 | traced |
| A39 | atlas-known ~0.89–0.90 under the same truths | Fig 2C | `16` §v2 (0.8941 v2, 0.9023 v1) | 2 | traced |
| A40 | IP-decoy proximity 27.0% vs 7.6% (13.0% full default) | Fig 2C, Fig 5D–F, S7, Limitations | `distance_crosscheck.tsv` (0.2699 / 0.0763 / 0.1300) | 3 | traced |
| A41 | Kinnex termini upstream: 1,228 at −5..0 vs 210 at 0..5 (v1, flagged as v1) | Fig 2C | `16` (v1 body) | 2 | traced |
| A42 | hexamer/A-fraction architecture 37.1% vs 14.7% null, modal +75 nt, crest 0.426 at +98 nt | Fig 2D | `07_curated_benchmark_report.md` | 5 | repo-source |
| A43 | 76.7% hexamer (35,712 / 46,524); +0.046 (0.7647 → 0.8110) | Fig 2D, Fig 5D–F, R2 | `call_primary/funnel.tsv` (35,712, 0.767274) and `validation.tsv` (0.811044) | 6 | traced |
| A44 | hexamer-fail 0.4874 vs trusted-novel 0.4853 | Fig 2D, Fig 5D–F | `validation.tsv` `CAL_atlas_novel_hexfail` (0.487407) | 2 | traced |
| A45 | IP filter on shipped arm 6.3% → 3.0% (null 2.2%) | Fig 2E, R2 | `07_curated_benchmark_report.md` | 3 | repo-source |
| A46 | IP filter buys +5.3 pp P for −1.2 pp R_det (mouse 1); v1 "+4.9 / −1.4" superseded | Fig 2E, S7, R2 | `19` §1 (explicitly supersedes) | 4 | traced |
| A47 | #96 impact: '+' 0, '−' 5,599 removed / 11,217 kept, residual 0; n 44,394 → 46,524 | Fig 2E, S7 | `VERIFIED_v2.md` "#96 IP-fix impact" | 6 | traced |
| A48 | mouse full-v2 sets restore 560 / 589 sites | Fig 2E | `VERIFIED_v2.md` prose says 560 / 589; its own arithmetic gives 588 for mouse 2 (25,938 + 588 = 26,526) | 2 | **inconsistent inside the source doc — left as written, listed for the PI** |
| A49 | v2 IP re-derivation flags 36–661 sites (≤3% change) | Fig 2E, Fig 5D–F | `call_primary_iprecomp_{p10_p30,m10_p29}/funnel.tsv` (36 and 661) | 3 | traced |
| A50 | competitor PBMC rows (scUTRquant 40,519 / 0.793 / 0.178 / 0.290; SCAPTURE 35,759 / 0.652 / 0.118 / 0.199; polyApipe 120,916 / 0.380 / 0.199 / 0.261; Sierra 106,170 / 0.256 / 0.135 / 0.177; scAPAtrap 787,138 / 0.100 / 0.300 / 0.150; shipped 277,111 / 0.118 / 0.156 / 0.134) | Fig 3A | `15` §3 (competitors were not re-run for v2 — stated in the text) | 24 | traced |
| A51 | competitor mouse rows (8 tools × P / R_det / F1_det × 2 mice) | Fig 3A | `15` §3 | 42 | traced |
| A52 | v2 PeakATail rows in Fig 3A (default, tier-1 IP, tier-1 no-IP; mice) | Fig 3A | `VERIFIED_v2.md` | 18 | traced |
| A53 | recall gap −12% PBMC, −18% / −17% mouse; F1 gap +0.020, +0.013 / +0.014 | Fig 3A, R3, positioning (f) | computed from `VERIFIED_v2.md` + `15` §3 polyApipe rows (−11.9 / −18.1 / −17.1; +0.0201 / +0.0133 / +0.0144) | 6 | traced |
| A54 | "+0.015 / +0.010 over polyApipe" | **R3 honesty rules** | v1 wording of `15` §3; v2 values are +0.020 / +0.013–0.014 | 2 | **stale v1 → fixed** |
| A55 | scTail R1 = 28 bp (BLOCKED) | Fig 3A | `15` §3 | 1 | traced |
| A56 | call counts 21k–787k, 37-fold | Fig 3B | `05_figure_index.md` | 3 | repo-source |
| A57 | shipped precision 0.118 PBMC → 0.369 testis → 0.447 Laughney | Fig 3B | `15` §3 (0.118, 0.369) and `09_headtohead_results.md` (0.447) | 3 | repo-source / traced |
| A58 | replicate reproducibility scAPAtrap 0.79–0.88, Sierra 0.79–0.83, PeakATail shipped 0.73–0.75, polyApipe 0.49–0.50; 62% depth-1, 0.31 vs 0.81 | Fig 3C, R3 | `05_figure_index.md` / `09_headtohead_results.md`; the v2 default value (0.776–0.779) now exists in `trade_reproducibility_reproducibility.tsv` but the bracket is unfilled | 11 | repo-source |
| A59 | long-read re-ranking ρ 0.90 at 7 of 8 arms, 0.80 at one | Fig 3D, R3 | `11_verifier_corrections.md` | 3 | repo-source |
| A60 | v2 compute: 34:37 / 12.53 GB / 785%; 32:59 / 11.12 / 818%; 14:10 / 3.07; 16:00 / 3.65; 27m43s | Fig 3E, R3, Limitations | `19` §3 and `VERIFIED_v2.md` compute table | 11 | traced |
| A61 | 23.5× / 21.5× / 7.5× / 6.1× less peak RSS; v1 3:45:53 / 293.7 GB at ~1.4 cores | Fig 3E | `VERIFIED_v2.md` compute table; `15` §5 | 7 | traced |
| A62 | panel context polyApipe 3:26 h, Sierra 2:43 h, scUTRquant 0:30 h | Fig 3E | `15` §5 | 3 | traced |
| A63 | all tools ≥95% gene-proximal | Fig 3F, R3 | `05_figure_index.md` | 1 | repo-source |
| A64 | mis-keyed-matrix nulls Fisher 38.6% / NB 21.3%; 20/20 null runs with q<0.05 | Fig 4A, R4, S6 | `13` §2 table gives 38.6 / 21.3; the "20/20" for the mis-keyed run is not in the reference set (it is the v1 `fdr_calibration` figure) | 3 | traced (2) / repo-source (1) |
| A65 | correctly-keyed shipped defaults 20.3% / 13.0% / 24.7%; 19–20/20 null runs | Fig 4A, R4, S6 | `14` arms table | 4 | traced |
| A66 | 1,230 cells (SPC 370 / RS 504 / ES 356) × 76,711 PAS; 20 perms; 6 arms × 21 runs; 0 failures | Fig 4B | `14` "Input" | 8 | traced |
| A67 | B0: 3.0% null p<0.05, 0.00% null q<0.05, 0/20 runs, 0/60 BH families; 24.8% of null p exactly 1; 2.94–3.13%; 2.6–4.0% | Fig 4B, R4 | `14` arms table + per-pair breakdown + verifier corrections | 9 | traced |
| A68 | marker double-dip 17.4%; denominator 629 / 6,453; reads-mode 9.6% and ~1,500 hits; NB floor 1e-4, 5.1%, 67%, min p 4e-182 | Fig 4B | `14` "Mechanism" 1–4 | 9 | traced |
| A69 | permutation-q survivors 186/307, 640/668, 1,733/1,754 | Fig 4B | `14` "Permutation-calibrated q" | 6 | traced |
| A70 | replication filter ≥2 samples; testis 2 mice; Laughney 17 samples; \|ΔPDUI\| ≥ 0.1 | Fig 4C/D, R4 | `13` §2 (17 GSMs); **`20` primary is 12 patients / 15 libraries** | 4 | **inconsistent — listed for the PI (Stage-3)** |
| A71 | spermatogenesis quarantined ρ −0.73 to −0.86; SPG 64/68 cells; per-gene replicate ρ 0.874; 200-permutation null | Fig 5A–C, R5 | `11_verifier_corrections.md` (quarantined run, so labelled) | 5 | repo-source |
| A72 | "6/6 literature genes in the expected direction" | Fig 5A–C | **`18` says the 16-gene panel does not reproduce (4 / 5 / 5 / 2) and that the 6/6 came from a different index and gene set** | 1 | **contradicted by `18` → disclosure added, claim left for the PI** |
| A73 | trusted-novel funnel 46,524 → 46,524 → 35,712 (76.7%) → 7,259 (15.6%) → 4,329 (9.3%); 3,578 genes | Fig 5D–F, T5 | `call_primary/funnel.tsv`; 3,578 in `REPORT_PROVISIONAL.md` §3 | 8 | traced |
| A74 | trusted-novel 0.4853 [0.474–0.497]; strong 0.5024; ≥1-mol sensitivity n 40,077 → 0.2894 | Fig 5D–F | `validation.tsv` | 6 | traced |
| A75 | enrichment "~65–180×" | Fig 5D–F, Discussion | `validation.tsv` enrichment 65.48 (t5) / 179.07 (t20); the caption says 65–179× | 2 | **inconsistent → fixed to 65–179×** |
| A76 | composition 69.4% intronic / 14.9% 3′UTR; atlas-known 69.9% 3′UTR | Fig 5D–F, R5 | `feature_breakdown.tsv` (0.6939 / 0.1489 / 0.6989) | 3 | traced |
| A77 | corrected filter flags 110,000 vs 119,518 raw peaks | Fig 5D–F, S7 | `16` §v2 | 2 | traced |
| A78 | post-hoc strata 3′UTR 0.769 / ≥5 mol 0.692–0.728 / 2-mol 0.401, 0.144 | Fig 5D–F, R5 | `stratified_concordance_trusted_novel.tsv` | 5 | traced |
| A79 | v1 comparison 52% / 6,628 sites | Fig 5D–F, R5, T5, positioning (f) | `16` §v2 table (0.5211; 6,628) | 2 | traced |
| A80 | clustering AMI 0.708 / ARI 0.502; Laughney median AMI 0.662 / ARI 0.463, 17/17 | S4 | `results/figures/manuscript/pbmc_novelty.tsv` (0.708) and `05_figure_index.md` | 5 | repo-source |
| A81 | ablation 275,370 sites → 14,891 gene totals, AMI 0.698 | S4 | `pbmc_novelty.tsv` (n_features 275370 / 14891; 0.6978) | 3 | traced (figure TSV) |
| A82 | 99.87% of annotated PAS carried another PAS's counts; bisected 2026-03-23 | S6 | `10_caller_fix_plan.md` / `09_headtohead_results.md` | 1 | repo-source |
| A83 | IP flag rate 3.8% on test data | S7 | `08_ipfilter_live_check.md` | 1 | repo-source |
| A84 | strategy spread 0.004; curated-atlas spread 0.355–0.492 | S3 | `07_curated_benchmark_report.md` | 3 | repo-source |
| A85 | S9 v1 outcome 0.717 / 0.741 / 0.755 (labelled v1); v2 0.7062 / 0.7450 / 0.7572; P 0.3032–0.3520 < 0.38 | S9 | `15` §1 and `19` §1 | 9 | traced |
| A86 | pre-registration `0e27b1a` 01:18:59, arms 02:59:36 | S9, Fig 1C, Methods 6 | `15` "Pre-registration timing" | 2 | traced |
| A87 | T5 7,259 (v2) / 6,628 (v1) | T5 | `16` §v2 | 2 | traced |
| A88 | chr19+21 slice F1 0.177 → 0.291 / 0.348; suite 1,277 passed | R1 | `06_roadmap.md` (slice) ; `19` header (1,277) | 4 | repo-source (3) / traced (1) |
| A89 | 90.2% keep ≥2 molecules under `-F 3844` (mouse 94.8 / 95.6%) | R2 | `15` §4 (labelled "v1 measurement, not re-derived in v2") | 3 | traced |
| A90 | curated-atlas precision of shipped arms 0.36–0.49, ~20× null | R2 | `07_curated_benchmark_report.md` | 3 | repo-source |
| A91 | Methods: PolyASite 2.0 569,005 sites; IP rule −10..+30 / ≥6 A / ≥70% A; Kinnex rule +1..+18 / ≥12 of 18 | Methods 4, 10 | `VERIFIED_v2.md` (569,005); `16` (both rules) | 8 | traced |
| A92 | Methods: clip clustering 25 bp single-linkage | Methods 3, R1 | `13` §1 / tool description | 1 | repo-source |
| A93 | Fig 6 header "17 samples"; R6 "across 17 lung-cancer samples"; abstract "[17]" | Fig 6, R6, abstract | `20` primary = **15 libraries / 12 patients**; 17 is the raw cohort in `13` addendum | 3 | **inconsistent — listed for the PI (Stage-3)** |
| A94 | journals/format numbers (IFs, word limits, APCs, figure counts) | §(a) | editorial, not results | 34 | n/a (not result numbers) |

**Table A subtotal: 457 values (423 result numbers + 34 editorial).**

---

## 3. Table B — the five `.caption.md` files

| # | quantity | caption | source checked | n | verdict |
|---|---|---|---|---:|---|
| B1 | 285,136 / 126,686; full-atlas 0.089 / 0.098 / 0.100 (569,005; 301,006) | final_benchmark | `VERIFIED_v2.md` §1 | 7 | traced |
| B2 | gate PASS 0.706 / 0.745 / 0.757; v1 0.717 / 0.741 / 0.755; −1.1 pp | final_benchmark | `19` §1 | 7 | traced |
| B3 | "+5.3 pp precision for −1.2 pp recall" | final_benchmark | correct values (`19` §1) but **bare "precision"/"recall"**, which the outline's honesty rules forbid | 2 | **forbidden phrasing → fixed to "+5.3 pp atlas-agreement precision for −1.2 pp R_det"** |
| B4 | compute 12.5 GB / ~28–35 min (was 294 GB / 3:46 h) | final_benchmark | `19` §3 | 4 | traced |
| B5 | recall −12% PBMC / −18% mouse; F1_det lead +0.020 / +0.014 | final_benchmark | computed in-script from the plotted TSV; the mouse figures are the two-mouse means (the outline quotes them per mouse, −18% / −17% and +0.013 / +0.014) | 4 | traced (granularity difference only; noted, not changed) |
| B6 | panel d 0.056 / 0.352 / 0.706; 72.2% single-molecule (121,041 / 167,565); mouse 50.3% / 48.8% | final_benchmark | `final_benchmark_tiers.tsv`; recomputed from the arm BEDs | 8 | **inconsistent with the text's 72.1% → fixed by stating both definitions** |
| B7 | panel-a sweep P 0.706 → 0.941 (PBMC), R_det 0.175 → 0.083; ≥10 mol 0.941 / 0.885 / 0.892 at 0.083 / 0.098 / 0.102; n 9,449–167,565 | trade_reproducibility | `trade_reproducibility.tsv` (all 15 sweep rows) | 12 | traced |
| B8 | defaults 0.7062 / 0.1754 (46,524); 0.7450 / 0.2048; 0.7572 / 0.2080 | trade_reproducibility | `19` §1 + the same TSV (reproduces to 4 dp) | 7 | traced |
| B9 | P@10 / P@100 ratios 0.74 default, 0.77 polyApipe, 0.46 / 0.33 / 0.24; scUTRquant 0.72 | trade_reproducibility | per-tool score TSVs (default 0.520871 / 0.706195 = 0.7375) | 6 | traced |
| B10 | reproducibility 0.776–0.779 @100 bp, 0.729–0.737 @25 bp; Jaccard 0.551 / 0.604; ≥1-mol 0.636–0.651; +14.3 / +12.5 pts | trade_reproducibility | `trade_reproducibility_reproducibility.tsv` | 10 | traced |
| B11 | chance ≤0.010, dotted line 0.0095, span 0.0017–0.0094, default 0.0019 vs scAPAtrap 0.0094 | trade_reproducibility | same TSV `chance_mean` column | 6 | traced |
| B12 | Sierra 0.790 / 0.834; scAPAtrap 0.884 / 0.793 | trade_reproducibility | same TSV | 4 | traced |
| B13 | ρ 0.655 (n 917; null 0.001 ± 0.032); 6,049 / 9,469 / 11,263; 0 in all 15 null pairings | trade_reproducibility | `18` "Switch test + cross-mouse replication" | 8 | traced |
| B14 | compute v1 3:45:53 / 293.7 GB → v2 34:37.49 / 12.53 GB; 23.5×; audit rows 3:44:37 / 239.5 → 32:59.38 / 11.12; uncontended 27:43; scAPAtrap 4.07 h / 12:58:50 | trade_reproducibility | `trade_reproducibility_compute.tsv`, `15` §5, `19` §3 | 12 | traced |
| B15 | funnel 46,524 input; 7,259 (15.6%, 3,578 genes; 4,329 strong) | trusted_novel_funnel | `call_primary/funnel.tsv`; `REPORT_PROVISIONAL.md` §3 | 6 | traced |
| B16 | 0.485 [0.474–0.497] (3,523 / 7,259); 0.239 / 0.066 / 0.012 | trusted_novel_funnel | `validation.tsv` | 7 | traced |
| B17 | GEM-X 0.515; pooled 0.582 / 0.669 / 0.714 | trusted_novel_funnel | `verifier_crosscheck/RESULTS.txt` | 4 | traced |
| B18 | funnel-stage concordances 0.765 → 0.765 → 0.811 (+0.046) → 0.485 (−0.326) → 0.502 (+0.017) | trusted_novel_funnel | `validation.tsv` stage rows | 8 | traced |
| B19 | panel-b sets 28,453 (0.894 / 0.707 / 0.389 / 0.142); 46,524 (0.765 / 0.558 / 0.281 / 0.098); 6,432 (0.487); 40,077 (0.289); null 0.007 (range 0.006–0.010) | trusted_novel_funnel | `validation.tsv` | 15 | traced |
| B20 | feature class TN 14.9 / 7.5 / 69.4 / 8.2; atlas-known 69.9 / 5.0 / 18.0 / 7.1 | trusted_novel_funnel | `feature_breakdown.tsv` | 8 | traced |
| B21 | decoy 27.0% (1,959 / 7,259), strong 16.8%, default 13.0%, atlas-known 7.6% (2,171 / 28,453) | trusted_novel_funnel | `distance_crosscheck.tsv` | 8 | traced |
| B22 | panel d 3′UTR 0.77 / 0.55 (1,081); other exon 0.61 / 0.38 (545); intron 0.41 / 0.16 (5,037); intergenic 0.47 / 0.20 (596); support 2 / 3–4 / ≥5 bins incl. the 5–9 (0.692 / 0.491, 672) and 10+ (0.728 / 0.625, 629) split | trusted_novel_funnel | `stratified_concordance_trusted_novel.tsv` | 22 | traced |
| B23 | v1 → v2: 0.521 on 6,628, decoy 24.5%; Δ −0.036 / +0.025; reference-set decoy shift 6.6 → 7.6, 11.2 → 13.0 | trusted_novel_funnel | `16` §v2 and the v1 `distance_crosscheck.tsv` (0.2453 / 0.0660 / 0.1124) | 9 | traced |
| B24 | 0.41 gap; ≈0.39 at ≥100 UMI; decoy-proximal 0.39 vs 0.52; 65–179× | trusted_novel_funnel | `validation.tsv`, `stratified_*.tsv` (0.3885 / 0.3885 / 0.5211 / 65.5–179.1) | 6 | traced |
| B25 | IP re-derivation 36 and 661 sites, ≤3% | trusted_novel_funnel | `call_primary_iprecomp_*/funnel.tsv` | 3 | traced |
| B26 | 46,544-site default, 20 sites off-contig, 99.96% | trusted_novel_funnel | `call_primary/funnel.tsv` | 3 | traced |
| B27 | "16 (v1 record)" in the Sources line | trusted_novel_funnel | `16` now carries a **§v2** section that is the v2 record for every number on this figure | 0 | **stale provenance → fixed** |
| B28 | 1,294 / 1,364 labelled cells; 98.8% agreement | spermatogenesis_final | `14` (1,294); `spermatogenesis_final_reference_lines.tsv` (0.9878 / 0.9877) | 3 | traced |
| B29 | 1,553 / 1,194 guarded genes of ~10,400; 31.5% / 30.2% (489 / 361); nulls 17.5% / 21.0%; z 10.5 / 5.8 | spermatogenesis_final | `18` claim paragraph | 11 | traced |
| B30 | lengthening 0.209 / 0.248 (z 3.1 / 2.1); excess 489 vs 325 (p 9.9e-9), 361 vs 296 (p 0.0125); 1 of 40 draws | spermatogenesis_final | `18` cautions; `spermatogenesis_final_reference_lines.tsv` | 12 | traced |
| B31 | Cliff's δ SPC–ES 0.56 / 0.61; SPC–RS 0.28 / 0.28; RS–ES 0.41 / 0.46; null ranges | spermatogenesis_final | `18`; `spermatogenesis_final_percell_resid.tsv` | 10 | traced |
| B32 | per-cell residual medians +0.0017 → −0.0036 → −0.0115 and +0.0028 → −0.0010 → −0.0084; 296 / 326.5 / 116.5 genes | spermatogenesis_final | figure TSVs | 9 | traced |
| B33 | 11,263 / 6,049 / 9,469 of 36,826 / 35,644 / 35,372; sign agreement 99.7–99.8%; 3.2–5.1×; 0 in 15 null pairings | spermatogenesis_final | `18`; `spermatogenesis_final_replication.tsv` | 12 | traced |
| B34 | TRUE hits 50,415 / 52,043; null 3.2% / 3.2%; 0 q<0.05 in 30 BH families | spermatogenesis_final | `18`; reference-lines TSV (0.03168 / 0.03163) | 5 | traced |
| B35 | ρ 0.655, r 0.853, null 0.001 ± 0.032, 200 perms, max \|ρ\| 0.110; 338 / 278 / 244 / 57 of 917; 167 (18.2%) vs 9.5%; binomial 0.017 | spermatogenesis_final | `18`; `gene_replication_summary.json` via reference-lines TSV | 15 | traced |
| B36 | per-gene medians 0.7867 / 0.7944 / 0.7320 and 0.8780 / 0.8863 / 0.8787; means 0.593 / 0.576 / 0.556 and 0.635 / 0.620 / 0.602 | spermatogenesis_final | reference-lines TSV | 12 | traced |
| B37 | UMI-weighted index 0.589 / 0.578 / 0.638 and 0.641 / 0.636 / 0.672; Cliff −0.95 / −0.86; protamine ~26% | spermatogenesis_final | `18` negative finding 1; figure TSVs | 9 | traced |
| B38 | literature panel: 16 measurable, 4 / 5 / 5 / 2 (named genes) | spermatogenesis_final | `18` negative finding 2 | 5 | traced |
| B39 | cell-weighted arm 214 / 197 genes; Wilcoxon shuffle best 8.2e-4 / 2.4e-5; SPG 64 / 68 | spermatogenesis_final | `REPORT_PROVISIONAL.md` §2.3b; reference-lines TSV | 6 | traced |
| B40 | "`manuscript/18` says … permutation p < 0.005"; Mann–Whitney 3.0e-39 / 2.6e-34; null range 1/21 ≈ 0.048 | spermatogenesis_final | `18` **no longer says that** — it was corrected to the null-range wording in commit `4c190b6`, the same commit that added this caption | 5 | **stale cross-document quote → fixed** |
| B41 | 12 patients / 15 libraries; 59 pairs; 74,954-PAS universe | laughney_switches | `20` primary block; `laughney_switches_cohort.tsv` | 4 | traced |
| B42 | 14,480 of 1,977,134 (0.73%); 13,826; 5,395 PAS / 2,687 genes / 47 pairs; K=3 4,937 | laughney_switches | `20` table; cohort TSV | 8 | traced |
| B43 | funnel 62,243 → 14,643 → 14,480 (163 = 129 + 34) → 13,826; null 1,888,790 tested, mean 1.0 (range 0–3), 0 later | laughney_switches | `laughney_switches_funnel.tsv` + cohort TSV | 11 | traced |
| B44 | consensus \|Δprop\| median 0.34 / mean 0.36 | laughney_switches | cohort TSV (0.341803 / 0.364821) | 2 | traced |
| B45 | panel e 58.6 / 4.0 / 13.8 / 20.4 / 3.2; roll-ups 96.8 / 76.3 / 62.6; 216 (6.4%) of 3,375; 497 (9.2%) | laughney_switches | `laughney_switches_context_summary.tsv` | 13 | traced |
| B46 | `20` disclosure-1 roll-ups 62.5% and 96.8%, with the own-gene wording correction to 58.6% | laughney_switches | `20` disclosure 1 vs the recomputation | 3 | traced (the caption states the correction explicitly) |
| B47 | 33% of top-30 gene rows misnamed; named examples; MetBone 0.015% vs 306k molecules; 13-patient sensitivity 14,541 / 13,888 (0.4%) | laughney_switches | `20` disclosures 1–2 | 8 | traced |
| B48 | "54.4M null tests across the 15 libraries → 10 nominal q<0.05" | laughney_switches | cohort TSV `null_tests_total` = 54,433,203; **`20` says "54.6M null tests cohort-wide"** | 3 | **inconsistent (54.4M vs 54.6M) — the caption self-scopes to the 15 primary libraries; listed for the PI to confirm the two are different scopes** |
| B49 | empirical p ≤ 0.091 (1/11 floor) | laughney_switches, trusted_novel_funnel | `20`; `13` addendum item 10 | 2 | traced |
| B50 | pre-registration times / codes on every caption (`4efeb125`, `9dfdefb`, 01:19, 02:59, 16:14) | all captions | `15`, `19` | 8 | traced |

**Table B subtotal: 372 values.**

**Total examined: 829 (795 result numbers).**

---

## 4. Stale v1 numbers where a verified v2 value exists

| # | text | was | should be | status |
|---|---|---|---|---|
| S-1 | R3 honesty rules, "never claim F1 superiority (…)" | +0.015 / +0.010 (`15` §3) | +0.020 PBMC / +0.013–0.014 mouse (`19` §1 vs `15` §3 polyApipe) | **fixed** |

Everything else that quotes a v1 value does so **with an explicit v1 label** and is correct as such:
`15` §3 competitor rows (competitors were not re-run on `9dfdefb`), the v1 compute pair, the v1 Kinnex
values (0.5211 / 6,628 / 24.5%), the v1 signed-distance counts (1,228 / 210), the 90.2 / 94.8 / 95.6%
`-F 3844` measurement, S9's v1 gate outcome, and the Fig 3C reproducibility rows.

---

## 5. Internal inconsistencies (same quantity, two values)

| # | quantity | values in the manuscript | resolution |
|---|---|---|---|
| I-1 | PBMC tier-1 single-molecule share | 72.1% (outline ×3) / 72.2% (final_benchmark caption + figure) / 72.5% (`15` §4) | **fixed.** The verifier's 72.1% is `1 − 62,110 / 222,955` on the **no-IP** tier-1 output (recomputed here from `peakatail_clipseeded_final_v2/pas_tier1.bed`: 160,894 single-molecule rows of 223,026 = 0.7214, byte-identical v1/v2). 72.2% is the same ratio on the **IP** arm (`1 − 46,524 / 167,565`), which is what panel d draws; 72.5% was the v1 IP arm. The verifier's value is kept in the text and the definition is now stated in all four places; the caption states both with their arms named. |
| I-2 | trusted-novel enrichment over null | "~65–180×" (outline ×2) vs "65–179×" (caption) | **fixed** to 65–179× (TSV: 65.48× at ≥5 UMI, 179.07× at ≥20). |
| I-3 | mouse-2 sites restored by #96 | 560 / 589 (outline, following `VERIFIED_v2.md` prose) vs 588 implied by the same document's arithmetic (25,938 + 588 = 26,526) | **left as written** — the collision is inside the verified source document. `VERIFIED_v2.md` should be corrected first. |
| I-4 | Laughney cohort size | "17 samples" (outline Fig 6 / R6 / abstract `[17]`) vs "12 patients (15 libraries)" (`20` primary, laughney caption) | **left as written** — Stage-3 placeholder; changing it changes what R6 asserts. For the PI. |
| I-5 | Laughney null test count | 54.4M (caption, `= 54,433,203` recomputed) vs 54.6M (`20`) | **left as written** — the caption self-scopes ("across the 15 libraries"); `20` says "cohort-wide". Confirm they are different scopes. |
| I-6 | recall / F1 gap granularity | outline per mouse (−18% / −17%; +0.013 / +0.014) vs caption two-mouse mean (−18%; +0.014) | **left as written** — both are correct at their stated granularity; the caption's are computed in-script from the plotted TSV. |

---

## 6. Forbidden phrasing

| check | result |
|---|---|
| bare "precision" for atlas agreement | **1 hit**: `final_benchmark.caption.md` (and the on-figure footer) "the IP filter buys +5.3 pp **precision** for −1.2 pp **recall**". **Fixed** to "+5.3 pp atlas-agreement precision for −1.2 pp R_det" in the caption, the figure footer and the generating script (whose own NAMING RULES already forbade it). Remaining uses are all qualified ("atlas-agreement precision", "curated-atlas precision", "atlas precision", "precision default", "precision-first") or are the honesty rule itself. |
| "beats" / "outperforms" / "superior" on F1 or accuracy | **none as a claim.** The only occurrences are prohibitions ("never claim F1 superiority", "**Accuracy leadership** … must NOT claim", "not the F1 leader"). The approved comparative wording is the verifier's "the most atlas-concordant de novo call set in the panel", followed immediately by the recall cost and "the F1 lead is not a meaningful margin — the claim is precision, not overall accuracy". |
| positive "trusted novel PAS" claim | **none.** Every occurrence is framed as a negative result; the positioning section lists "'Trusted novel PAS' as a positive claim" under "What we must NOT claim"; the funnel caption opens with "negative result" and forbids quoting the 50 bp / 100 bp pooled numbers as meeting the target. |
| clustering-novelty claim | **none.** Dropped in the header, retired in the title list, and S4 states the gene-sum ablation. |
| the retired 0.9986 quoted as a result | **absent from the whole scope.** (It appears only in `manuscript/figures/README.md`, and there only as the description of the retired `null_control` figure.) |
| unverified Kinnex molecule counts ("104 M") | **absent from the whole scope.** ⚠ It does survive in a reference document — `13_reliability_positioning.md` §3, "(we have 104 M poly(A)-verified molecules)". `13` is a pre-registration, not a results document, but the phrase should be struck there before any draft quotes it. Listed for the PI (out of this audit's edit scope). |
| "records" vs "molecules" for Kinnex support | correct everywhere: support is stated in UMIs, and the Limitations paragraph explicitly retires the record-counting first-version figure. |

---

## 7. `[final run pending]` placeholders (13 in the outline; none in the captions)

| # | line | placeholder | which job fills it |
|---|---|---|---|
| P-1 | number policy ¶ | "[final run pending]" (the meta-statement listing what is still bracketed) | now partly obsolete: `18` and `20` have landed. Provenance note added; no bracket filled. |
| P-2 | journals, NatComms row | "[pending Amir's rerun]" | Amir's #67 SWITCH_CELLTYPE rerun (Fig 6) |
| P-3 | abstract | "[ρ]" | the spermatogenesis v2 re-run (after #96/#97). **Note: `18`'s surviving claim is a monotone-gene-fraction statement, not a ρ** — the bracket may need to change shape, not just be filled. |
| P-4 | abstract | "[17]" | Stage-3 Laughney. `20` (v1 code) reports 12 patients / 15 libraries; the running `stage3_laughney_v3` chain will restate it on `9dfdefb`. |
| P-5 | Fig 3C | "[final run pending — not in `15_final_gate.md`; compute on the two mouse `pas_PRESPEC_precision_default.bed` files …]" | **already computed**: `trade_reproducibility_reproducibility.tsv` gives 0.776–0.779 @100 bp (0.729–0.737 @25 bp) for the v2 default. Filling it is a PI call — the value's source is a figure TSV, not yet a verified write-up. |
| P-6 | Fig 3D | "[quarantined; regenerate]" (long-read re-ranking) | not regenerated for v2; the trade caption says so explicitly under "Not drawn". |
| P-7 | Fig 4B | "[final run pending — Stage 3]" (cross-dataset confirmation on Laughney) | `stage3_laughney_v3` chain (running) |
| P-8 | Fig 4C | "[pending]" (replication-filter yields, between-donor specificity) | `stage3_laughney_v3` chain |
| P-9 | Fig 4D | "[pending]" (effect-size floor, volcano, positive controls) | `stage3_laughney_v3` chain |
| P-10 | Fig 5A–C | "[all numbers final run pending; regenerate per `11_verifier_corrections.md`]" | the spermatogenesis v2 re-run. `18` already supplies the v1-code (`4efeb125`) values; the disclosure is now attached. |
| P-11 | S4 | "[pending]" (within-gene usage-fraction spaces on re-keyed matrices) | novelty-ablation re-run; the outline states the claim does not depend on it |
| P-12 | S7 | "[regenerate; Stage 3]" (long-read-measured IP fraction per tool) | Stage-3 / Kinnex figure regeneration |
| P-13 | R3 / R4 / R5 | "[final run pending]" repeats of P-5, P-7 and P-10 | as above |

---

## 8. Claims whose evidence document says the opposite

| # | claim in the manuscript | what the evidence document says | action |
|---|---|---|---|
| C-1 | Fig 5A–C plans "6/6 literature genes in the expected direction" | `18` negative finding 2: of 16 measurable genes **4 shorten in both mice, 5 lengthen in both, 5 discordant, 2 uninformative**; the quarantined 6/6 "came from a different metric and gene set"; "the literature-panel leg of Fig 5 is dropped or shown with these numbers" | **disclosure added to the bullet** (the plan text itself is a claim change — for the PI) |
| C-2 | Fig 5A–C plans a per-cell 3′UTR length index gradient (ρ) | `18` / caption: the **per-gene medians are not monotone** and the **across-gene UMI-weighted index reverses at RS→ES**; the surviving claim is the monotone-gene *fraction* | disclosure added; the claim reshape is for the PI |
| C-3 | Fig 5A–C says the `wul` vs `wdi` discrepancy is "resolved before the atlas arm is used" | caption: "The discrepancy is reproduced here, **not resolved**"; `18` documents it as a name collision | disclosure added |
| C-4 | Fig 6C/D plan "TAM/myeloid … worked case" and "Geneview vignettes for 2–3 headline genes" | `20` disclosure 1: **no named top-switch list may be published** until PAS→gene re-assignment (33% of the top-30 rows name a spanning/readthrough model; 7 named genes have zero assigned PAS cohort-wide) | **for the PI** — Stage-3 section, and the fix is a claim change |
| C-5 | Fig 4C plans "the expected false-positive rate of the filter" | `20`: the null resolves only to **empirical p ≤ 0.091** (10-permutation floor) — "report as 'none in 10 nulls', **never** as an FDR estimate" | **for the PI** — the wording "expected false-positive rate" invites exactly the forbidden reading |
| C-6 | Limitations: "the switch-test calibration rests on one mouse and one tissue until the Laughney replication lands" | `18`: arm B0 **re-validates on the final caller in both testis mice**; `20`: the Laughney patient-level replication **has landed** (v1 code, v2 chain running) | **for the PI** — updating it changes what the limitation asserts |

Also carried forward but already stated correctly in the text: the MetBone exclusion premise was wrong
(`20` disclosure 2 — the caption says so), the Kinnex donor mismatch does not explain the atlas-novel
gap (`16`), and the trusted-novel residual is the IP **rule's** leniency, not the strand bug (`16` §v2).

---

## 9. Edits made (all in this pass)

| file | edit |
|---|---|
| `manuscript/01_outline_and_journals.md` | R3 honesty rules: stale v1 F1 gap "+0.015 / +0.010" → "v2: +0.020 PBMC / +0.013–0.014 mouse … the v1 wording … is superseded" |
| `manuscript/01_outline_and_journals.md` | Fig 2 figure block, R1 and R2: 72.1% single-molecule share now carries its definition (no-IP arm, 1 − 62,110 / 222,955) and names the IP-arm (72.2%) and superseded v1 (72.5%) values |
| `manuscript/01_outline_and_journals.md` | R2: molecule-unit sentence no longer attributes 72.1% to `pas_support.tsv` (that file, tier-1 rows, gives 73.3% because it is pre-IP-filter); `pas_support.tsv` is kept as the source of the `-F 3844` figures, which is what it actually supports |
| `manuscript/01_outline_and_journals.md` | enrichment "~65–180×" → "65–179×" in all three places it appears (Fig 5D–F, R5, Discussion), with "(65.5× at ≥5 UMI, 179.1× at ≥20; `validation.tsv`)" spelled out at first use — matches the funnel caption |
| `manuscript/01_outline_and_journals.md` | abstract note "(199 words by `wc -w`)" → 197 (measured; 198 including the `>` blockquote marker) |
| `manuscript/01_outline_and_journals.md` | number-policy ¶: Stage-3 provenance note added — `18` is "still the `4efeb125` record" (v2 re-run not yet run) and `20` is "a v1-code record … with the Stage-3 v2 chain on `9dfdefb` running now"; **no bracket filled** |
| `manuscript/01_outline_and_journals.md` | Fig 5A–C: `18` disclosure attached (literature panel 4 / 5 / 5 / 2, the 6/6 provenance, the UMI-weighted reversal, `wul`/`wdi` not resolved) |
| `manuscript/figures/final_benchmark.caption.md` | "+5.3 pp precision for −1.2 pp recall" → "+5.3 pp atlas-agreement precision for −1.2 pp R_det" |
| `manuscript/figures/final_benchmark.caption.md` | panel d: definition sentence added reconciling 72.2% (IP arm, drawn) with the verifier's 72.1% (no-IP arm) and `15` §4's 72.5% |
| `manuscript/figures/final_benchmark.{png,pdf}` | regenerated from the corrected script (edge check 0 ink; PDF fonts CID TrueType only, no Type 3); panels a–d are **pixel-identical** to the previous version except the footer text, plus one layout fix: panel b's three-line annotation was right-anchored at 0.90 instead of 0.985 so it no longer overlaps the "F1 0.1" isoline label |
| `scripts/manuscript_figures/final_benchmark.py` | `ip_trade` strings (v2 and v1) requalified; per-version `single_mol_note` added and used in the sidecar caption; docstring records the definition nit; panel-b annotation anchor 0.985 → 0.900 |
| `manuscript/figures/trade_reproducibility.caption.md` + `scripts/manuscript_figures/trade_reproducibility.py` | the `18` citation "(v1 code)" → "(**still the `4efeb125` record — its v2 re-run after #96/#97 has not happened**)"; the on-figure footer was left untouched, so the figure was not regenerated |
| `manuscript/figures/spermatogenesis_final.caption.md` | the note "…`manuscript/18` says … 'permutation p < 0.005'" now reads as a historical correction: `18` was changed to the null-range wording in commit `4c190b6`, the same commit that added this caption |
| `manuscript/figures/trusted_novel_funnel.caption.md` | Sources line: `16` cited as "(v1 record)" → "**§v2** (the v2 record for every number on this figure) … the pre-§v2 body of 16 is the v1 record" |

Nothing was invented; the only values introduced that were not already in a reference document are
`160,894 / 223,026` and `1 − 62,110 / 222,955`, both recomputed here from the primary BEDs and the
`VERIFIED_v2.md` table, and both stated as such in the text.

---

## 10. Left for the PI (a fix would change what the paper asserts)

1. **Fig 6 / R6 / abstract "17 samples"** vs `20`'s pre-registered primary of **12 patients / 15
   libraries** (MetBone excluded by pre-registration, plus label-confirmation drop-out). Decide the
   headline unit before Fig 6 is written.
2. **Fig 5A–C's "6/6 literature genes"** and the ρ-shaped abstract bracket — `18` retires both.
3. **Fig 6C/D named-gene vignettes** — `20` forbids a named top-switch list until PAS→gene
   re-assignment (issue #99).
4. **Fig 4C's "expected false-positive rate of the filter"** — `20` allows only "none in 10 nulls,
   empirical p ≤ 0.091", never an FDR.
5. **Limitations' "until the Laughney replication lands"** — it has landed (v1 code).
6. **Fig 3C's replicate-concordance bracket** — the value (0.776–0.779) exists in a figure TSV but not
   yet in a verified write-up.
7. **`VERIFIED_v2.md` internal slip**: "restore 560 / 589 sites" vs its own `25,938 + 588 = 26,526`.
8. **`13` §3's "104 M poly(A)-verified molecules"** — an unverified Kinnex molecule count still sitting
   in a pre-registration document; strike it before any draft quotes it.
9. **`20` "54.6M null tests cohort-wide"** vs the recomputed 54,433,203 for the 15 primary libraries —
   confirm these are different scopes and say which one the paper quotes.
