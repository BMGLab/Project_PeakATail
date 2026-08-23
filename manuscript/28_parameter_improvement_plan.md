# Improvement plan from Amir's parameter reference (2026-08-23)

**Source:** `PeakATail_Parameters_Reference.pdf` (26 pp, generated from the `RunConfig` schema, so the
defaults are exactly what the pipeline runs). It documents knobs our measurement programme never
touched, and one of them bears directly on the finding we treated as a hard ceiling.

## 0. What it confirms (no action beyond hygiene)

- `--ip-filter-mode` defaults to **`annotate`**, which keeps every PAS and only records a flag. Our
  benchmark launcher passes `--ip-filter-mode filter` explicitly, so **every published number is
  safe**; the `peakAtail-prime` branch's new default was silently a no-op for exactly this reason,
  caught independently by the accuracy verifier and fixed in `6954082`.
- `--polya-min-umis` default 1, `--polya-min-clip` 6, `--polya-min-purity` 0.8 — as measured.
- **`--min-pas-prominence` does nothing on the default `lambda_gradient` strategy.** Check whether the
  retired `parameter_sweep` figure swept it; if so, that arm was a no-op and must be labelled.
- Use `--polya-mode filter` (the tool's own tier-1 selection) instead of our post-hoc `awk '$5>0'`.

## 1. The four leads, ranked

### Lead A — `--max-pas 5` is a hard ceiling per peak (highest value)
Amir states it plainly: *"A gene with more than five true PAS can never be fully recovered no matter
what else you change."* Our decomposition ([25](25_competitive_position.md) §4) found **79.96% of
missed atlas sites have no peak of any kind within 100 bp**, and the atlas places ~20 sites in an
average detected human gene while the default emits ~4 calls per gene. If one coverage peak spans a
3′UTR holding eight real sites, five is all we can ever emit. **This is a candidate-generator ceiling
we imposed, not one the data imposed** — and if it moves, our headline "the ceiling is the evidence in
the library" claim needs revising. Test: sweep `--max-pas` 5 → 10 → 20 → unlimited.

### Lead B — the ~90–105 nt coverage offset, applied to tier-2 (most likely to pay)
The reference explains that coverage-based peak 3′ ends *"stop ~90–105 nt short of the true cleavage
site because 10x R2 coverage runs out before the poly(A) junction"* — the rationale for
`--auto-cleavage-offset`. Clip-anchored tier-1 calls do not suffer this; **coverage-only tier-2 calls
do**, which would explain their 0.057 precision as a *systematic offset* rather than noise. We killed
"tier-2 rescue" after testing hexamer and internal-priming gating; we never tested offset correction.
166,355 currently worthless sites are at stake. Test: `--auto-cleavage-offset` and a fixed
`--cleavage-offset` sweep, scored on tier-2 alone.

### Lead C — `--pas-gap 100` and `--min-pas-spacing` merge neighbouring sites
Real PAS are frequently 30–80 bp apart; a 100 bp minimum gap merges them. Caveat from our own
measurement: the scorer's nearest-neighbour matching already credits both members of a merged doublet
(33% of current recall is many-to-one), so **splitting may not move benchmark recall** — but it does
change per-PAS counts, which feed every switch test. Evaluate on resolution (P@10/P@25) and on
per-PAS quantification, not on recall alone.

### Lead D — `--lambda-fold-change 2.0` is the actual candidate-generation threshold
The reference's own tuning map lists it first under "find more PAS". Lowering it creates more
candidate peaks, attacking the same peakless class as Lead A. Test jointly with A.

### Also: the cell-filter seam
`--min-read 1500`, `--min-cells 3`, `--min-pas-per-cell 50` explain why **31.2% of tier-1 clip
clusters never reach `pasbed.bed`** (measured, [23](23_algorithm_roadmap.md)). Test an exemption for
high-support clusters.

## 2. Pre-registered acceptance criteria (fixed BEFORE any sweep number exists)

A parameter change may become a `peakAtail-prime` default only if, on **both species**:
1. atlas-agreement precision at matched call count falls by **≤ 0.005**;
2. detected-gene recall rises by **≥ 0.010**;
3. **the long-read (Kinnex) truth moves in the SAME direction** — this is the anti-atlas-prior guard,
   added because a transferable scoring model already fooled the atlas while losing 15% of
   long-read-verified sites ([24](24_prime_preregistration.md), prime TASK D);
4. no more than a 2× increase in wall time or 1.5× in peak RSS.
Anything that passes 1–2 but fails 3 is reported as a negative result and ships behind a flag, off.

## 3. Execution order

1. **Slice sweep** (PBMC chr19+21 and mouse1 chr18+19) over Leads A–D, factorial where cheap,
   scored against atlas AND Kinnex, with matched-call-count comparisons.
2. **Full-BAM confirmation** of any arm that clears §2 on the slice.
3. **Prime update**: implement/default the winners; everything else ships flagged-off with its numbers.
4. **Four-way re-benchmark** (shipped → v1 → v2 → prime) and the version-progression figure.

## 4. To reconcile with Amir

His page-2 benchmark table (shipped 0.161/0.196/0.177; clip_seeded tier-1 only 0.359/0.338/0.348)
disagrees with our verified numbers (tier-1 ≥1 molecule 0.3520/0.2685/0.3046 — [19](19_final_gate_v2.md)).
Most likely a different recall denominator (full atlas vs detected-gene restricted). One set must win
before either is quoted.
