# Roadmap to submission

**Target: Genome Biology (Method article), submission early November 2026.** Fallback chain that
needs no restructuring: Genome Research → Briefings in Bioinformatics / NAR. Assessment of
2026-08-13: ~40% ready overall — engineering ~80%, evidence ~30%, text ~20%.

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
- [ ] **Publish GitHub issues + PR** — review `manuscript/github/`, then run `manuscript/github/publish_github.sh` *(Ebru)*
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
- [~] Head-to-head on **pbmc_10k_v3** — RUNNING (2026-08-19): PeakATail, Sierra (txdbmaker shim), scAPAtrap,
      SCAPTURE, polyApipe, scUTRquant all live; scTail inapplicable by chemistry (needs R1-preserved
      libraries; 10x v3 R1 = CB+UMI only) — documented for the comparison table.
- [~] **GSE104556 testis** build RUNNING (STARsolo path, mouse GRCm38 r102; ~6-9h)
- [x] **FDR calibration** — DONE for fisher, verdict NEGATIVE: 20/20 null runs report q<0.05 hits
      (mean 42/run vs TRUE 40); causes = read pseudoreplication + marker double-dip. Issue 7 filed;
      nb_pairwise calibration DONE — also anti-conservative (21.3% null p<0.05); D4 cells-mode needs CLI exposure. Manuscript must not use
      fisher q-values as-is.
- [x] **FDR calibration (nb_pairwise)** — DONE (2026-08-19), verdict NEGATIVE too: 20/20 null runs report q<0.05 hits (mean 350/run vs TRUE 316; 21.3% null p<0.05, 12.0% null q<0.05, min null p 2.5e-36 — milder than fisher but still ~4x nominal; driver: per-PAS plug-in dispersion clipped at 1e-4 floor in 21% of null tests carrying ~half the false hits, no shrinkage, marker double-dip). nb_pairwise q-values also unusable at face value → manuscript significance needs D4 cells-mode (CLI exposure) or permutation-calibrated thresholds; figure: manuscript/figures/fdr_calibration.{png,pdf}.
- [x] `--ip-filter` verified LIVE (annotate flags 3.8%, filter drops exactly flagged, default=off) —
      resolves the lg_ip_off duplicate as expected behavior; docs rows to delete (comment staged for #69)
- [ ] Merge the point-mode benchmark PR after Amir's review

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
