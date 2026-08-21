# Roadmap to submission

**Target: Genome Biology (Method article), submission early November 2026.** Fallback chain that
needs no restructuring: Genome Research → Briefings in Bioinformatics / NAR.

**Status 2026-08-21 — PI DECISION: reliability-first positioning (13_reliability_positioning.md).**
Clustering-novelty claim DROPPED. Claims = low-FP PAS detection, reliable cell-type APA switches, UTR
length, trusted de novo PAS. Precision-first default output pre-registered BEFORE the final run
(tier-1 ∩ IP-pass ∩ ≥2 molecules; gate P ≥ 0.50 both datasets). Stage 3 = reliability program
(recalibration + replication filter + trusted-novel definition validated on Kinnex).

**Status 2026-08-19 — POSITIONING CHANGED.** The evidence base is now large and largely negative on
accuracy: PeakATail ranks last among de novo tools on public PBMC (F1 0.134 vs polyApipe 0.261) and
mid-pack on testis, and neither differential test is FDR-calibrated. Accuracy-leadership and
FDR-controlled-discovery framings are both retired. What remains defensible: (1) peak-based clustering
novelty, (2) replicate reproducibility at high recall (0.73–0.75 vs polyApipe 0.49–0.50), (3) the
benchmark methodology itself (nulls, per-dataset denominators, reproducibility metric, the ~90–105 nt
cleavage-offset discovery). Readiness: engineering ~85%, evidence ~65%, text ~20% — but the evidence
now argues for a different paper than originally planned. **Decision pending with Ebru/PI: reframe
around clustering + reproducibility, or pause and improve the caller first (add poly(A) soft-clip
evidence, the axis that separates the winners).**


### Stage 3 status (2026-08-21)
- **Switch recalibration v2 — DONE, verified SOUND** ([14_switch_calibration_v2.md](14_switch_calibration_v2.md)):
  only Fisher `--count-mode cells` with marker pre-selection OFF controls FDR (conservative: 3.0% null
  p<0.05, 0/20 null runs with a q<0.05 hit). Every shipped-default arm (top-200 marker pre-selection)
  is anti-conservative; nb_pairwise is tail-inflated by its dispersion floor. Decision: ship
  cells-mode Fisher, `--marker-top-n 0`, + replication filter + effect floor; permutation-calibrated
  q mandatory whenever pre-selection or nb is used. Tool issue drafted for Amir
  (github/issue9_switch_marker_preselection.md).
- **Reliability tooling — DONE** (`scripts/reliability/`, 39 unit tests): `replication_filter.py`
  (≥2 samples same direction, discordance veto, effect floor, null-based false-replication rate) and
  `trusted_novel_pas.py` (5-stage pre-registered funnel + Kinnex `validate` vs shuffled null,
  `meets_target_0.70`). Numbers deferred to the final run.
- **Outline v0.2** repositioned to reliability-first; clustering demoted to supplementary utility.
- **Stage 1c — DONE, verified FIXED** (feat/polya-evidence 0eebfe7 + verifier patch): BED col5 = distinct
  (CB,UMI) molecules; `--polya-min-reads` → `--polya-min-umis`; `pas_support.tsv` sidecar (reads /
  molecules / −F 3844 subsets / window reads per PAS); fallback double-count fixed (−0.0005%);
  coordinates + tier membership byte-identical at defaults. `--polya-clip-filter f3844` implemented,
  DEFAULT OFF (call-set change; −12.6% mouse tier-1 rows, to be scored post hoc from the sidecar).
  Mouse1 post-hoc ≥2-molecule tier-1: P 0.692 / R 0.216 (post hoc; the pre-registered test is the final run).
- **FINAL Stage-2 run — DONE, verified SOUND** ([15_final_gate.md](15_final_gate.md)): pre-registered precision
  default P@100 0.717 (PBMC) / 0.741 / 0.755 (mice) → gate P≥0.50 PASS on all three; R_det 0.17–0.20
  (≈ polyApipe or below), F1 tied with polyApipe. Original two-sided gate still FAILS on P for every
  ≥1-molecule output (reported as sensitivity arm). Compute limitation: 280 GiB / 3:46 h on PBMC.
  Stage-1d defects for Amir: IP-filter minus-strand window, memory, output bookkeeping
  (github/issue10_stage1d_ipfilter_memory.md).
- **Kinnex trusted-novel validation — DONE, verified FIXED** ([16_trusted_novel_kinnex.md](16_trusted_novel_kinnex.md)):
  pre-registered definition yields 6,628 trusted-novel sites but only **52% (not 70%)** within 25 bp of a
  long-read 3′ end; enrichment 70–250× over null; residual error = intronic sites + internal priming the
  current IP rule misses (24.5% coincide with Kinnex IP decoys). Reported as a negative result; no
  "trusted novel" claim; stricter v2 needs pre-registration + held-out validation. Stage-1d: `--ip-rule kinnex`.
- **Stage 1d (PI directive 2026-08-21 08:00):** fix IP-filter minus-strand window (branch fix/ip-filter-strand)
  and memory/CPU (perf/clip-memory); PR each to Amir; after merge pull the new version and RE-RUN the
  benchmark + Stage 3 from a fresh frozen snapshot. PR #93 (Stage 1b+1c) MERGED into develop 2026-08-21 (Amir approved); issues #94 (switch marker pre-selection) and #95
  (Stage 1d) open. Publishing from this session is blocked by the permission classifier → publish scripts in manuscript/github/.
- **Stage 1d-A IP-filter strand fix — DONE, verified FIXED, ready for PR** (branch fix/ip-filter-strand, 3 commits;
  19 new tests; '+' byte-identical; measured: corrected rule ≈ re-run would give PBMC default P@100 0.717 → ~0.706
  with +2,100 sites / +0.5 pp recall; mice 0.741 → 0.746, 0.755 → 0.759; all gates still pass). Publish:
  `bash manuscript/github/publish_stage1d_ipfilter.sh` (after publish_stage1bc.sh).
- **Stage 1d-B memory/CPU — DONE, verified FIXED, ready for PR** (branch perf/clip-memory, 8 commits, +42 tests):
  root cause of the 280 GiB was the DENSE TF-IDF in clustering (4 float64 copies of cells × PAS), not the caller;
  sparse TF-IDF is bit-identical; peak calling now per (contig, strand) in parallel with a deterministic merge;
  read_check reordered; CB filter streamed. PBMC 10k: **293.7 GB → 12.5 GB, 3h46m → 27m43s** (mouse1 23 GB → 3.7 GB,
  61 → 9 min), all outputs byte-identical to the final run. Publish: `bash manuscript/github/publish_stage1d_perf.sh`.
  Not done (follow-ups): --tiles/--pipeline remain broken (separate issue), single-pass both strands, chunked mmread.
- Stage-3 Laughney cohort run — DONE 14:24 (9h06m, 13.6 GB; frozen 4efeb125): 17 datasets, unified PAS space 505,197
  sites, every dataset's var pas_id maps into it, obs prefixes/counts verified (verify_cohort.json); cohort IP-flag
  rate 27.9%. Switch queue (17 GSM × 1 true + 10 null, 3 concurrent) started automatically; replication follows.
- **Spermatogenesis control on final caller — DONE, verified FIXED** ([18_spermatogenesis_final.md](18_spermatogenesis_final.md)):
  SPC>RS>ES survives as a per-gene equal-weight claim (monotone-shortening 31.5%/30.2% vs null 17.5%/21.0%;
  Cliff SPC vs ES 0.56/0.61); arm B0 calibration re-validates; cross-mouse replication 6.0-11.3k PAS/pair,
  99.7% sign agreement, 0 in 15 null pairings. Negatives: UMI-weighted index reverses at RS→ES (protamine
  ceiling); 16-gene literature panel does not reproduce (dropped); per_isoform degenerate-pair defect →
  issue11 draft. Switch queue DONE 15:47 (17/17 GSMs); replication driver running.
- **PRs #96 + #97 MERGED by Amir 2026-08-21 ~16:10** → develop tip 9dfdefb (IP-strand fix + perf, on 1b/1c).
  **v2 re-run LAUNCHED 16:14** from frozen tools/pa-polya-run-9dfdefb3 (suite: 1,277 passed, only the 2 known env
  failures), OUT_TAG=final_v2, same gates. Expect: PBMC default P@100 ~0.706 (from the measured IP-fix impact),
  mice ~0.746/0.759; runtime ~30 min PBMC.
- **v2 gate — VERIFIED FIXED** ([19_final_gate_v2.md](19_final_gate_v2.md)): PRESPEC P@100 0.7062/0.7450/0.7572,
  gate PASS all three; #97 identity confirmed byte-for-byte; #96 impact residual 0; PBMC 12.5 GB / ~28-35 min.
  Paper numbers move to v2.
- **Stage-3 v1 chain COMPLETE 17:08** (replication done; PROVISIONAL): patient-unit K2 PAS 14,541 replicated
  (13,888 with floor) over 1.98M (pair,feature); gene level 9,618; K3 4,976/3,559; **0 replicated in all 10 null
  combos at every level**. CAVEAT: MetBone was INCLUDED — the pre-registered primary (13 addendum 8) excludes it;
  verifier is running the excluded primary + full checks.
- **Stage-3 replication — VERIFIED FIXED** ([20_stage3_replication.md](20_stage3_replication.md)): primary
  (MetBone excluded, 12 patients): 14,480 replicated (pair,PAS) at K=2 (13,826 with effect floor) = 5,395 PAS /
  2,687 genes / 47 pairs; gene level 9,579; K=3 4,937; **0 replicated in all 10 null combos**. v1-code record —
  regenerated by Stage-3 v2. Disclosures: top-gene names unreliable in overlapping loci (issue #99: CD68 etc.);
  MetBone premise wrong (clip-rate sampling bug, same issue); null floor 1/11.
- **16 v2 section added (verified)**: trusted-novel negative robust to the fixes (0.485; decoy 27%).
- Fig 2 + funnel + outline moved to verified v2 numbers (v0.4; VERSIONS switch keeps v1 renders reproducible).
- **Stage-3 v2 chain LAUNCHED 17:55** (root /mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3, code 9dfdefb,
  ~9× faster peak calling; smoke test isolates #96: minus-strand IP removal 21.9% → 18.4%, plus strand byte-identical;
  both replication arms automated; ETA complete ~01:00 2026-08-22). Then: Stage-3 v2 verifier → 20 update → figures.
- Amir queue: issues #94, #95 (§1 fixed by #96, §2 by #97; §3 + --ip-rule kinnex open), #98, #99. (`stage2_final_launch.sh`,
  gate P ≥ 0.50 on the pre-registered default) → Stage-3 numbers (replication on 17 Laughney
  patients + 2 mice; trusted-novel on Kinnex).

## Owners

| Who | Role |
|---|---|
| Ebru (@ebrukocakaya) | Manuscript lead, analyses, decisions on framing |
| Amir (@TRextabat) | Tool development, sweep pipeline on ssd2, PR reviews |
| Claude | Analyses, figures, benchmark harness, drafts, PR implementation — all human-reviewed |

## Phases

### P0 — this week (unblock + de-risk)
- [ ] **Switch rerun** — apply the one-line `main.nf` fix (→ `unified/multi_sample_merged.bed`, NOT
      the 14%-coverage cohort bed) + `nextflow -resume`. *Owner: Amir. Everything biological waits on this.*
- [x] Null-control analysis proving the 18.4M-reference benchmark indefensible *(done, verified)*
- [x] **Curated re-benchmark** vs PolyASite 2.0 point-mode/strand-matched + null *(done, verified — see 07_curated_benchmark_report.md)*
- [x] **Motif validation** — atlas-independent accuracy *(done, verified; discovered the ~90–105 nt cleavage-offset)*
- [x] **GitHub issues + PR published** (Ebru ran publish_github.sh): issues #67–#72 filed, PR #73 open
      with @TRextabat as reviewer, branch pushed.
- [ ] **Publish the follow-ups**: `manuscript/github/publish_followups.sh` — issue 7 (neither diff test
      is FDR-calibrated), issue 8 (three CellRanger-input bugs), comments on #67/#69. Also still unstaged:
      an issue for the head-to-head negative result + the soft-clip-evidence proposal it implies.
- [x] **Decision checkpoint — TRIGGERED (2026-08-19).** Point-mode precision vs curated PolyASite 2.0
      is 0.36–0.49 @100bp (~20× null) — under the 70% bar; recall is arithmetically capped at 0.06–0.15
      by set size, so the ≥60% bar is unreachable *by design*. Decisions:
      (a) **reframe**: lead with clustering novelty + enrichment-over-null + canonical motif architecture,
          not absolute atlas accuracy; (b) the old ROADMAP bars must be restated (issue for Amir/PI);
      (c) the **cleavage-offset correction** (issue 6) is the concrete path to recovering precision —
          AATAAA mode at +75 nt implies cleavage ~+90–105 nt past the reported peak end;
      (d) the head-to-head matters even more: all tools get scored under the identical fair matcher,
          so *relative* standing is what the paper argues.

### P1 — weeks 1–3 (competitive evidence)
- [x] **Head-to-head COMPLETE (2026-08-19) — see 09_headtohead_results.md. RESULT IS NEGATIVE for
      accuracy claims:** PeakATail ranks LAST among de novo tools on pbmc_10k_v3 (F1 0.134 vs polyApipe
      0.261) and mid-pack on testis (F1 0.258-0.260, 2nd of 4). Not a threshold artifact (top-N ranking
      gives identical precision) nor a cleavage-offset artifact (@200bp still 0.166 vs 0.410). The
      discriminator is EVIDENCE TYPE: winners use poly(A) soft-clips (polyApipe) or sequence models
      (SCAPTURE) — PeakATail calls from coverage shape alone. Strategic consequence: drop any accuracy-
      leadership framing; lead with clustering novelty + reproducibility (PeakATail 0.73-0.75 vs
      polyApipe 0.49-0.50). Tool fix worth doing: add soft-clip polyA evidence to peak calling.
      Panel completed: PeakATail, Sierra, scAPAtrap, SCAPTURE, polyApipe, scUTRquant on both datasets;
      scTail excluded by chemistry (needs R1-preserved libraries; 10x v3 R1 is CB+UMI only). Every tool
      needed environment archaeology or bug fixes to run at all — reportable in its own right.
- [x] **Three CellRanger-input bugs found in PeakATail** (commits ad59cdd, 6df8eed, 2e7fc0d): GEM-group
      CB suffix, unmapped reads with reference_end None, underscore-bearing RG corrupting the sample_CB
      composite. All were invisible to the STARsolo-only test suite. Issue 8 asks for a CellRanger
      regression BAM in CI.
- [x] **GSE104556 testis** built (STARsolo, GRCm38 r102; both mice, CB/UB-tagged) and fully benchmarked
- [x] **FDR calibration** — DONE for fisher, verdict NEGATIVE: 20/20 null runs report q<0.05 hits
      (mean 42/run vs TRUE 40); causes = read pseudoreplication + marker double-dip. Issue 7 filed;
      nb_pairwise calibration DONE — also anti-conservative (21.3% null p<0.05); D4 cells-mode needs CLI exposure. Manuscript must not use
      fisher q-values as-is.
- [x] **FDR calibration (nb_pairwise)** — DONE (2026-08-19), verdict NEGATIVE too: 20/20 null runs report q<0.05 hits (mean 350/run vs TRUE 316; 21.3% null p<0.05, 12.0% null q<0.05, min null p 2.5e-36 — milder than fisher but still ~4x nominal; driver: per-PAS plug-in dispersion clipped at 1e-4 floor in 21% of null tests carrying ~half the false hits, no shrinkage, marker double-dip). nb_pairwise q-values also unusable at face value → manuscript significance needs D4 cells-mode (CLI exposure) or permutation-calibrated thresholds; figure: manuscript/figures/fdr_calibration.{png,pdf}.
- [x] `--ip-filter` verified LIVE (annotate flags 3.8%, filter drops exactly flagged, default=off) —
      resolves the lg_ip_off duplicate as expected behavior; docs rows to delete (comment staged for #69)
- [x] Point-mode benchmark PR #73 merged; Stage-0 PRs #81-#84 merged; Stage-1 PR #92 merged.
- [ ] **#67 SWITCH_CELLTYPE rerun is now UNBLOCKED** (#82 merged): Amir backfills `layers['counts']`
      from `06_preprocessing/*/preprocessed.h5ad`, passes `--pasbed` + `--counts-layer`, reruns only
      B2+combine+switch (half a day). Use `--count-mode cells` (#86) for fisher.

### P1b — IN FLIGHT (2026-08-19, parallel)
- [~] **PBMC novelty experiment** — does PAS-only clustering recover GEX cell types on public data
      (replicating Laughney AMI 0.662 / ARI 0.463), and do peak-space-only populations survive confound
      checks? THE decisive analysis for the paper's new positioning.
- [~] **Consolidated benchmark figure** + single source-of-truth table + depth-stratified reproducibility
- [~] **Spermatogenesis positive control** (testis, both mice) — textbook 3'UTR shortening gradient
- [~] **Kinnex long-read Tier-1 truth** — 81.7M poly(A)-verified molecules; atlas-independent re-ranking
      of every tool. If the unfavourable ranking flips here, that changes the accuracy story.
- [~] **SCAPTURE testis mouse2** rerun (first attempt destroyed by the ssd1 disk-full event; now on ssd0)

### Stage 1 — DONE (2026-08-20, verified SOUND)
- [x] **Clip-seeded calling built, verified SOUND, and MERGED as PR #92** (develop 6697b0d, 2026-08-20
      18:20). **CI green on develop post-merge (3.11 + 3.12, run 32402648097)** — the gate from #84 works and the merged caller passes it. chr19+21 slice, identical scoring: shipped F1 0.177 → clip_seeded both-tiers **0.291**
      (1.65×), clip-supported tier-1 **0.348** (P 0.359 / R 0.338, at the evidence recall ceiling).
      Verifier independently re-scored the BEDs, hand-checked clip_site on 365,606 real reads
      (0 mismatches), and confirmed three-path agreement. +71 tests, zero regressions.
- [ ] **Stage 1b (opened 2026-08-20): tier-1 quantification.** Stage 2 testis arms show clip_seeded
      counts tier-1 PAS from clip reads only -> 10x less matrix mass, 23% of real cells lost to min_read.
      Coordinates/accuracy unaffected. Fix = count all read ends in a cleavage window; acceptance = mass
      within 0.8-1.2x shipped, >=98% STARsolo cell recovery, pas.bed byte-identical. Then regenerate
      Stage-2 matrices (second re-run; disclosed).
- [x] **Stage 2 PBMC GATE: NOT CLEARED (two verifier passes, 12_stage2_gate.md).** Default output F1
      0.290 / P 0.308 fails the P≥0.38 floor. Post-hoc sweeps examined and disclosed: read-count ≥2 was
      duplicate-inflated (23.7% single molecules); molecule-based ≥2 UMIs gives P 0.593 / R 0.189 /
      F1 0.286 (numeric pass, recall below polyApipe 0.199); matched-recall precision 0.51 vs 0.38.
- [ ] **Stage 1c: UMI-dedup clip evidence** (-F 3844 semantics; BED score = molecules, reads as annotation;
      align --polya-min-reads unit). Then re-run Stage 2 once, with the default set BEFORE the run.
- [x] ~~Stage 2 LAUNCHED (2026-08-20)~~ from the rebased branch (f0370f7 on develop a3ddf4e, byte-identical
      slice output vs the verified build). Four concurrent arms via scripts/benchmark_tools/stage2_launch.sh:
      PBMC clip_seeded, PBMC clip_seeded + --ip-filter filter (human FASTA), testis m1, testis m2. Each
      self-scores (both tiers + combined) on its dataset's denominator and writes DONE.ok/FAILED.err.
      Gate unchanged: F1@100 > 0.261 on full PBMC. Pre-registered: tier-1 P is not final until the
      internal-priming filter runs on human; tiers reported separately, always.

### P2 — weeks 3–6 (ground truth + the make-or-break experiment)
- [~] **Kinnex PBMC long-read** Tier-1 truth: 36GB downloading (dedup FLNC BAM + GEM-X mapped BAM);
      extraction plan in scripts/benchmark_tools/kinnex_truth_plan.md
- [ ] Known-biology positive controls: T-cell activation genes, plasma-cell IGHM switch (checklist in
      `02_validation_plan.md`)
- [ ] **Novelty experiment**: a cell population/state visible in PAS space that GEX clustering misses,
      orthogonally validated — this decides whether the title leads with clustering
- [~] Release engineering (issue #70): branch feat/release-engineering READY (CI yml, validated
      CITATION.cff, README 12/12 commands; suite 925 pass / 5 pre-existing fails documented) —
      publish block staged in publish_github.sh; PyPI/Zenodo remain

### P3 — weeks 6–10 (biology + figure freeze)
- [ ] Laughney per-celltype APA switches: stage/metastasis programs per cell type (from P0 rerun)
- [ ] TAM 3′UTR programs across Normal → StageI → Met; geneview panels for headline genes
- [ ] Figure freeze: 6 main figures per `01_outline_and_journals.md`, all through adversarial verification
- [ ] Methods drafted against the *verified* docs (issues 3/5 resolved)

### P4 — weeks 10–14 (writing + submission)
- [ ] Full draft, internal review round with Amir + PI
- [ ] bioRxiv preprint, then Genome Biology submission

## Collaboration conventions (agreed workflow)

- **Tool changes** (`BMGLab/PeakATail`): feature branches off `develop`, conventional-commit titles,
  PR to `develop`, at least one human review (Amir reviews Claude/Ebru branches; Ebru reviews Amir's
  manuscript-relevant ones). Label `manuscript` on anything blocking the paper.
- **Sweep pipeline / ssd2** is Amir's area — others read-only; requests go through issues, never
  direct edits.
- **Manuscript workspace** (`PeakATail_wd`, branch `reorg-manuscript`): figures are only committed
  with their re-runnable script + `.tsv` audit file; every figure passes an adversarial verification
  pass before it enters `manuscript/figures/`.
- **GitHub posting from the BioLab box** currently authenticates as @yasinkaymaz — publish via the
  reviewed script, or re-auth (`gh auth login`) to post as yourself.
- Machine quirk that corrupts analyses silently: **always `export LC_ALL=C`** before sort/bedtools.

## Standing risks

1. **Field velocity**: scTail (GB 2025), scPAISO (2025) — each month costs novelty margin.
2. **Bar risk**: if curated-atlas precision lands under 70%, pivot framing to clustering novelty (P0 checkpoint).
3. **Single-dataset risk**: until PBMC + testis run, every result is Laughney-only.
4. **The lg_ip_off duplicate** means the IP-filter story currently rests on one arm.
