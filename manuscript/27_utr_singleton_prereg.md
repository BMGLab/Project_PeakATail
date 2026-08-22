# 27 — Pre-registration: 3′UTR promotion of single-clip-molecule tier-1 calls

**Status: PRE-REGISTRATION, written 2026-08-21, BEFORE any confirmatory evaluation.** Nothing in this
file is a result. It fixes, in advance, the exact rule, the datasets, the metrics, the decision
thresholds and the falsifiers, so that the promotion *can* ship later without the project having to
claim that a post-hoc observation was a finding.

Governed by [24](24_prime_preregistration.md) §3.1 (adoption criterion), §3.3 (Rule T — no threshold
may be tuned on the metric it is reported against) and §3.4 (post-hoc results are exploratory until
re-pre-registered). Discovery record: `results/perf_gap/D3_head_to_head/` §2.3,
`results/perf_gap/VERIFY/README.md` correction 2, roadmap entry [22](22_performance_roadmap.md) §4.2,
engineering draft `manuscript/github/issue20_threeutr_singleton_promotion.md`. Competitive context:
[25](25_competitive_position.md) §4 (the promotion can move at most ~+0.04 R_det; the ceiling is
elsewhere).

---

## 0. The one thing that must be read first

**Every number this idea currently has is post-hoc, and none of it can serve as evidence for the rule.**
The promotion was found by searching the *same* PBMC and mouse-1 data that the Stage-2 gate was scored
on — a 69-arm rule sweep (`results/perf_gap/D2_rule_sweep/`), a 17-candidate filter search
(`D3_head_to_head/` §3.1) and a promotion scan (`D3_head_to_head/` §2.3), plus the verifier's own
addition of the 3′UTR stratum. It is the best-measured candidate in that search, which is precisely why
it is the one most exposed to selection.

**Worse, and this is the reason this document exists:** the promotion is a *deterministic
post-processing function of already-frozen call sets*. Re-running it on PBMC 10k v3 or testis mouse 1
would reproduce §5's numbers exactly and would add **zero** information. A "confirmatory run" on the
discovery data would be theatre. Therefore:

> **The deciding evidence for this pre-registration is (i) testis mouse 2, on which this rule has never
> been evaluated, and (ii) any dataset not used in the discovery search. PBMC 10k v3 and testis mouse 1
> are the DISCOVERY SET; their numbers are reported for transparency and are labelled
> `DISCOVERY SET — not evidence` wherever they appear.**

---

## 1. The exact rule

### 1.1 Inputs, fixed now

| element | fixed definition |
|---|---|
| base call set | the shipped pre-registered precision default: clip-seeded **tier-1** ∩ **IP-pass** (`--ip-filter --ip-filter-mode filter`; `ema/experimental/internal_priming.py::ip_window`, transcript-relative **−9..+30** (40 nt), flagged if ≥6 consecutive A or A-fraction ≥ 0.70) ∩ **≥2 distinct clip molecules** |
| promotion candidates | tier-1 ∩ IP-pass calls with **exactly 1** distinct clip molecule (`clip_umis` = BED column 5 = 1). The IP filter remains a **hard veto**: an IP-flagged call is never promoted, whatever its annotation. |
| call position | the strand-aware **3′-most base**: `end − 1` on `+`, `start` on `−` — the same point the scorer uses |
| 3′UTR test | that base lies inside an annotated `three_prime_utr` feature of a gene on the **same strand**, under the hierarchical class **3′UTR > other exon > intron > intergenic** |
| annotation build | **Ensembl 99 (GRCh38)** for human, **Ensembl 102 (GRCm38)** for mouse — the builds already used by the benchmark. The build is recorded in the run manifest; changing it is an amendment (§7). |
| hexamer test | at least one of the **12 canonical hexamers** (`scripts/reliability/trusted_novel_pas.py` list) occurring **fully** within transcript-relative `r ∈ [−40, −5]`, computed on one symmetric ±45 bp `bedtools getfasta -s` window (the [16](16_trusted_novel_kinnex.md) rule) |

### 1.2 The two arms, and which one is the default candidate

> **R2 — PRIMARY (the only arm that may become a default).**
> `output = (tier-1 ∧ IP-pass ∧ molecules ≥ 2) ∪ (tier-1 ∧ IP-pass ∧ molecules = 1 ∧ 3′UTR ∧ hexamer_any12)`

> **R1 — SECONDARY, pre-declared as a labelled higher-recall arm and NOT a default candidate.**
> `output = (tier-1 ∧ IP-pass ∧ molecules ≥ 2) ∪ (tier-1 ∧ IP-pass ∧ molecules = 1 ∧ 3′UTR)`

**Why the hexamer-gated arm is primary, decided now.** The stated success criterion is *no precision
loss with a recall gain*. On the discovery set R1 loses precision on mouse 1 (§5), so pre-registering
R1 as the default candidate would be pre-registering a rule the available data already says will fail
its own criterion. R2 is the arm that can satisfy it. R1 is nevertheless retained, because its recall
gain is larger and a labelled sensitivity arm has value ([25](25_competitive_position.md) §6) — but it
is declared **now**, in advance, as an arm that may only ever be *reported*, never adopted as the
default. Choosing between R1 and R2 *after* seeing the confirmatory numbers is forbidden; the order is
fixed here.

**Free parameters, all fixed here (Rule T).** Molecule count for promotion = exactly 1. Hexamer set =
the 12-hexamer list, window `[−40, −5]`. Annotation feature = `three_prime_utr`, same strand.
No numeric threshold is introduced by this rule, so there is nothing left to tune. **If any element of
§1.1 or §1.2 is changed after the confirmatory numbers are seen, the result is void and the amendment
must be recorded under §7 with the pre-amendment numbers preserved.**

### 1.3 Shipping form

Behind a **flag**, default OFF, until §3 passes: `--promote-utr-singletons` (name provisional; the
engineering spec is `manuscript/github/issue20_threeutr_singleton_promotion.md`). Promoted calls carry
an explicit per-site flag column in the output so that any downstream analysis can exclude them.

---

## 2. Datasets and scoring protocol

Unchanged from [24](24_prime_preregistration.md) §2 and [19](19_final_gate_v2.md): one scorer
(`scripts/benchmark_tools/score_tool.py`), `P@100` against PolyASite 2.0 representative sites
(569,005 human / 301,006 mouse, strand-matched, `bedtools closest -s -d -t first`), `R_det@100` against
the detected-gene-restricted atlas (**285,136** human / **126,686** mouse), `F1_det` their harmonic mean.

| dataset | role | why |
|---|---|---|
| **testis mouse 2** | **DECIDING** | the rule has never been evaluated on it — verified: `results/perf_gap/D3_head_to_head/tsv/08*.tsv` covers PBMC and mouse 1 only, and no promotion arm for mouse 2 exists anywhere in `results/perf_gap/` |
| PBMC 10k v3 | DISCOVERY SET — not evidence | the rule was found here; re-running reproduces §5 exactly |
| testis mouse 1 | DISCOVERY SET — not evidence | idem |
| any dataset outside the discovery search (e.g. a Laughney cohort sample, Stage 3) | **supporting, deciding if available** | genuinely held out; but no atlas-independent truth and a different tissue, so reported with that caveat |

**How weak "held out" is here, stated in advance.** Mouse 2 is a *replicate of mouse 1 from the same
study, tissue, chemistry and atlas* — the manuscript already carries that caveat ([01](01_outline_and_journals.md)
honesty rules; [19](19_final_gate_v2.md)). It is held out only in the sense that this rule has never
been run on it, and the two mice are known to track each other closely on related rules (the
cross-replicate arms of `D2_arm_sweep.tsv` group 7 reproduce between the two mice to within 0.005 in
P@100 and 0.001 in R_det: `M1_crossrep_w25` 0.7346 / 0.2237 vs `M2_crossrep_w25` 0.7351 / 0.2237, and
`M1_union_sup2_or_1crossrep25` 0.7099 / 0.2394 vs `M2_...` 0.7144 / 0.2404). **A pass on mouse 2
is therefore weak, largely-anticipated evidence and must be reported in exactly those words: "confirmed
on one held-out replicate of the discovery mouse, same study and chemistry."** It is *not* independent
replication, and a pass on it does not license dropping the flag-default-OFF requirement of §1.3 without
a genuinely independent dataset. If a fourth, non-discovery dataset with a usable detected-gene
denominator becomes available before submission (a second donor — see [26](26_second_donor_preregistration.md)
— or a Stage-3 cohort sample), it is added to the deciding set; that decision is made **now**, not after
seeing mouse 2.

---

## 3. Success criteria — all of them, on the deciding data

Let `ΔP = P@100(promoted) − P@100(v2 default)`, `ΔR = R_det@100(promoted) − R_det@100(v2 default)`,
`ΔF1 = F1_det(promoted) − F1_det(v2 default)`, per dataset. The caller is deterministic, so these are
exact quantities, not estimates.

> **R2 PASSES and may become the default iff ALL of C1–C5 hold.**

**C1 — no material precision loss.** `ΔP ≥ −0.005` on **every** dataset evaluated (deciding *and*
discovery). The `−0.005` is the project's standing materiality allowance ([24](24_prime_preregistration.md)
§3.1), not a noise allowance. Report the strictness label: **STRICT** if `ΔP ≥ 0` everywhere.

**C2 — a real recall gain.** `ΔR ≥ +0.010` on every dataset evaluated.

**C3 — F1 does not regress.** `ΔF1 > 0` on every dataset evaluated.

**C4 — it survives the atlas-independent (long-read) check.** PBMC only, since the Kinnex panel is
human-only. All four must hold, against the v2 default's measured values (`D2_rule_sweep/kinnex/`,
`D3_head_to_head/tsv/11_*.tsv`, `CHECK22` §2):

| | bar, fixed now | v2 default's value |
|---|---|---|
| C4a combined-set Kinnex x3p ≥5 UMI within 25 bp | **≥ 0.7347** (default − 0.03) | 0.7647 |
| C4b **the promoted calls' own** concordance | **≥ 0.5300** — the midpoint between the singleton pool (0.2952) and the default (0.7647): a promoted call must look **more like a call than like a discard** | — (pool 0.2952) |
| C4c combined-set Kinnex **internal-priming decoy** proximity within 25 bp | **≤ 0.1400** (default + 0.01); *higher is worse* | 0.1300 |
| C4d combined-set **hard-FP rate** (no atlas site ≤100 bp **and** no Kinnex x3p ≥5 UMI 3′ end ≤25 bp) | **≤ 0.1510** (no worse than the default) | 0.1510 |

**Disclosure required by honesty, not by the protocol:** C4a, C4b and C4d are set with knowledge of the
discovery values (§5), and PBMC/x3p *is* the discovery truth. They are therefore **weak** criteria — they
can falsify the rule but a pass on them is not independent confirmation. C4c is the one long-read
quantity that has **not** been measured for either arm, so it is a genuine pre-registered test. The
GEM-X panel and the ≥20-UMI threshold are reported alongside as robustness, with no bar attached.

**C5 — de novo integrity.** (i) No atlas information may enter the rule at any point — the 3′UTR test
reads the GTF only. (ii) Because the rule reads an annotation, any sentence in the paper calling
PeakATail's output "de novo" must either exclude promoted calls or be restated; the promoted-call flag
(§1.3) exists to make that mechanical. (iii) The residual circularity the verifier itself flagged —
annotated 3′UTRs derive from transcript ends, so a weak prior is shared with long-read 3′ ends — is
declared here and must be stated wherever the rule's numbers appear. It is far weaker than the
motif/PolyASite circularity that sank the hexamer arm, but it is not zero.

**FAIL ⇒ the rule does not become the default.** It stays behind the flag, R2 and R1 are reported as
labelled exploratory arms only, and the miss is recorded as a negative result in the same way
[16](16_trusted_novel_kinnex.md)'s trusted-novel result was.

### 3.1 Falsifiers, declared in advance

* **Non-transfer.** If on mouse 2 the promoted subset's own `P@100` falls more than 0.05 below the
  mouse-2 default's `P@100`, declare non-transfer and stop, regardless of the combined numbers.
* **Recall bought from noise.** If C2 passes while C4a or C4d fails, the gain is not a gain — report as
  such and stop. This is the failure mode `results/algo_headroom/A1_end_pileup/` measured for the
  end-position channel and the one [24](24_prime_preregistration.md) §3.2(2) exists to catch.
* **Annotation dependence.** If repeating the mouse arms under a different Ensembl build moves `ΔP` or
  `ΔR` by more than 0.01, the rule is annotation-fragile; report it and keep the flag off by default.

---

## 4. Reporting requirements

Whatever the outcome, report: per dataset `n`, `P@100`, `R_det@100`, `F1_det`, the promoted subset's own
`n` and `P@100`; PBMC Kinnex x3p t5/t20 and GEM-X t5 concordance, decoy@25 and hard-FP rate for the
combined set **and** for the promoted subset alone; the strictness label from C1; and the
DISCOVERY/DECIDING label on every row. Report R1 in the same table, always labelled *"sensitivity arm,
not a default candidate"*. Never report a promoted arm against a competitor without the matched-N
context ([25](25_competitive_position.md) §9).

---

## 5. The post-hoc record — NOT EVIDENCE, recorded so the rule cannot be quietly changed

Source: `results/perf_gap/D3_head_to_head/tsv/08_promote_singlemolecule_pbmc.tsv`,
`08b_…mouse1.tsv`, `11_promotion_atlas_independent_check_pbmc.tsv`; verified in
`results/perf_gap/VERIFY/README.md` (correction 2 and "what reproduced exactly").
**These numbers were obtained on the data the rule was selected on. They are listed here only so that a
later reader can confirm the pre-registered rule is the rule that was explored — they do not and cannot
support adoption.**

| arm | dataset | n | P@100 | R_det@100 | F1_det | ΔP | ΔR |
|---|---|---:|---:|---:|---:|---:|---:|
| v2 default | PBMC `DISCOVERY` | 46,524 | 0.7062 | 0.1754 | 0.2811 | — | — |
| **R2** (3′UTR ∧ hex) | PBMC `DISCOVERY` | 54,220 | 0.7194 | 0.2015 | 0.3148 | **+0.0132** | **+0.0260** |
| R1 (3′UTR) | PBMC `DISCOVERY` | 61,232 | 0.7034 | 0.2150 | 0.3294 | **−0.0028** | +0.0396 |
| v2 default | mouse 1 `DISCOVERY` | 26,255 | 0.7450 | 0.2048 | 0.3213 | — | — |
| **R2** | mouse 1 `DISCOVERY` | 30,974 | 0.7549 | 0.2402 | 0.3644 | **+0.0099** | **+0.0354** |
| R1 | mouse 1 `DISCOVERY` | 34,865 | 0.7359 | 0.2545 | 0.3782 | **−0.0091** | +0.0497 |
| **any arm** | **mouse 2** | — | — | — | — | — | **never evaluated** |

Promoted subsets alone, and the atlas-independent view (PBMC):

| set | n | P@100 (atlas) | Kinnex x3p t5 ≤25 bp | GEM-X t5 ≤25 bp | hard-FP rate |
|---|---:|---:|---:|---:|---:|
| default (≥2 mol) | 46,524 | 0.7062 | 0.7647 | 0.7907 | 0.1510 |
| all singletons (the pool) | 121,041 | 0.2158 | 0.2952 | 0.2957 | 0.5906 |
| promoted by **R2** | 7,696 | 0.7994 | **0.5635** | 0.5997 | 0.1171 |
| promoted by R1 | 14,708 | 0.6947 | **0.5143** | 0.5623 | 0.1850 |
| combined **R2** | 54,220 | 0.7194 | 0.7361 | 0.7636 | 0.1462 |
| combined R1 | 61,232 | 0.7034 | 0.7045 | 0.7358 | 0.1592 |
| mouse 1 promoted by R2 / R1 | 4,719 / 8,610 | 0.8099 / 0.7080 | no mouse long-read truth | — | — |

**Read against §3, the discovery set already says: R1 fails C1 on mouse 1 (−0.0091) and fails C4b
(0.5143 < 0.5300); R2 passes C1–C3 on both discovery datasets and passes C4a/C4b/C4d.** That is exactly
why R2 is primary and why nothing here counts as evidence: the criteria were written knowing this.
The rule's whole remaining claim to credibility rests on **mouse 2** and on C4c — and §2 records how
little a mouse-2 pass is worth on its own.

**Context that limits the prize** ([25](25_competitive_position.md) §4): classes (a)–(d) are only
20.04 % of the PBMC recall gap and 15.47 % of the mouse gap; 79.96 % / 84.53 % of missed atlas sites have
no peak of any kind. This promotion can buy at most ~+0.04 R_det. It is a scoring improvement, not an
answer to the ceiling.

---

## 6. Relationship to the other singleton routes

| route | status | why it is not this |
|---|---|---|
| hexamer union (`≥2 OR 1-mol strong hexamer`) | **rejected** | the gain is largely atlas-circular: +0.078 atlas but only +0.007 long-read over a size-matched random top-up (`D2` §5) |
| cohort / cross-replicate corroboration | deferred, separate pre-registration | needs a multi-sample mode; non-circular but has no long-read confirmation (`D2` §2(4); issue17) |
| a calibrated per-site score | [24](24_prime_preregistration.md) / prime | may recover the same sites; if prime ships a score, this rule must be re-evaluated **against** it, not stacked on it without measurement |
| tier-2 rescue, relaxed clip detection, end-position calling, deconvolution, cleavage-offset correction, `--ip-rule kinnex` | **measured dead** | do not resurrect ([22](22_performance_roadmap.md) §4.6; [16](16_trusted_novel_kinnex.md) IP-rule note) |

---

## 7. Amendments

Any change to §1, §2 or §3 after a confirmatory number has been seen voids the pre-registration for that
run. Amendments are appended here, dated, with the reason and with the pre-amendment text preserved
verbatim — the convention used in [24](24_prime_preregistration.md) §6.

*No amendments.*

---

*This document is PROVISIONAL until a verifier passes on it. It contains no result; it contains a
promise about how a result will be judged.*
