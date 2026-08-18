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
- [ ] **Curated re-benchmark** vs PolyASite 2.0 point-mode/strand-matched + null *(Claude — running)*
- [ ] **Motif validation** (AAUAAA enrichment, internal-priming rate) — atlas-independent accuracy *(Claude — running)*
- [ ] **Publish GitHub issues + PR** — review `manuscript/github/`, then run `manuscript/github/publish_github.sh` *(Ebru)*
- [ ] Decision checkpoint: does point-mode precision vs curated atlas clear the ROADMAP bar (≥70%)?
      If not, reframe the paper around the clustering novelty before more benchmarking is sunk.

### P1 — weeks 1–3 (competitive evidence)
- [ ] Head-to-head on **pbmc_10k_v3**: PeakATail vs Sierra, scAPAtrap, SCAPTURE, polyApipe, scTail
      (+ scUTRquant if build is tractable). Installs + harness in progress; 44 GB BAM downloading.
- [ ] Add **GSE104556 testis** (spermatogenesis gradient = benchmark + positive control in one)
- [ ] **FDR calibration** of `switch diff` via permuted cell labels (needs P0 switch rerun)
- [ ] Verify `--ip-filter` live; resolve the lg_ip_off duplicate-arm question (issue 3)
- [ ] Merge the point-mode benchmark PR after Amir's review

### P2 — weeks 3–6 (ground truth + the make-or-break experiment)
- [ ] **Kinnex PBMC long-read** Tier-1 truth: download, extract poly(A)-tailed 3′ ends per barcode,
      score PeakATail sensitivity/specificity against it
- [ ] Known-biology positive controls: T-cell activation genes, plasma-cell IGHM switch (checklist in
      `02_validation_plan.md`)
- [ ] **Novelty experiment**: a cell population/state visible in PAS space that GEX clustering misses,
      orthogonally validated — this decides whether the title leads with clustering
- [ ] Release engineering (issue 4): CI, PyPI, CITATION.cff, Zenodo DOI on `v0.3.0`

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
