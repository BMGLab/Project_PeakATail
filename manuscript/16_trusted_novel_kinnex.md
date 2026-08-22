# Trusted de novo PAS — pre-registered definition vs Kinnex long-read truth (verified FIXED, 2026-08-21)

**Outcome: the pre-registered definition does NOT meet its target.** It was set in
[13](13_reliability_positioning.md) §3 before any number existed: trusted-novel = clip-supported ∧ ≥2
molecules ∧ not internal-priming ∧ canonical hexamer in −40..−5 ∧ ≥100 bp from any atlas site; target
≥70% within 25 bp (strand-matched) of a poly(A)-verified Kinnex long-read 3′ end. Input = the final
PBMC precision default (IP arm, code `4efeb125`, 44,413 sites). Verifier reproduced every count
exactly; one README fix (IP-flag column lives in `03_gtf_annotation/default/annotatedpas.bed:10`, not
the 9-column top-level file) and one robustness table were added.

## Funnel (verified)

44,394 on primary contigs → clip-supported 44,394 → ≥2 molecules 44,394 → not-IP 44,394 (the run was
already IP-filtered; re-deriving IP from the genome would flag 389–736, ≤3% change) → hexamer (any
of 12) **34,192** (77.0%) → ≥100 bp from PolyASite 2.0 **6,628 trusted-novel** (14.9% of input) →
AATAAA/ATTAAA **4,037 strong** (9.1%). 3,350 genes. PolyA_DB could not be used (hg19 only).

## Long-read concordance (≤25 bp, same strand; Kinnex x3p truth; Wilson 95% CI)

| set | n | ≥5 UMI | ≥20 | ≥100 | ≥500 | target 0.70 |
|---|---:|---:|---:|---:|---:|---|
| **trusted-novel (pre-registered)** | 6,628 | **0.521** [0.509–0.533] | 0.261 | 0.073 | 0.013 | **not met** |
| trusted-novel, strong hexamer | 4,037 | 0.525 | 0.256 | 0.074 | 0.014 | not met |
| atlas-known complement, hexamer-pass (calibration) | 27,564 | 0.902 | 0.715 | 0.393 | 0.144 | — |
| full precision default | 44,394 | 0.786 | 0.576 | 0.289 | 0.101 | — |
| atlas-novel, hexamer-FAIL (diagnostic) | 5,972 | 0.528 | 0.281 | 0.078 | 0.016 | — |
| sensitivity: ≥1 molecule trusted-novel | 37,438 | 0.309 | 0.103 | 0.017 | 0.003 | — |

Gene-body-shuffled null ≈ 0.007 at ≥5 UMI (10 seeds) → enrichment ~70–250×, empirical p at the 1/11
floor: the sites are far from random, but the target is an absolute fraction. Robustness (verifier):
GEM-X truth 0.554; x3p+GEM-X pooled 0.623 at 25 bp (0.698 at 50 bp, 0.739 at 100 bp) — the
pre-registered 25-bp metric is missed under every truth choice. Donor mismatch cannot explain the gap:
the atlas-known complement scores 0.90–0.94 under the same truths.

## What the data say about *why* (diagnostic, not a re-definition)

- The atlas-novelty stage costs the most concordance (0.828 → 0.521); the hexamer stage adds +0.04 on
  the whole set but **nothing among atlas-novel sites** (hexamer-fail 0.528 vs trusted-novel 0.521).
- 68.7% of trusted-novel sites are **intronic** (atlas-known complement: 68.6% in 3′ UTRs); only 15.6%
  are in annotated 3′ UTRs.
- **24.5% of trusted-novel sites lie within 25 bp of a Kinnex internal-priming decoy terminus** (6.6%
  for the atlas-known reference) — this measurement **stands**. ~~The caller's IP rule (−10..+30, ≥6 A or
  ≥70% A) is looser than Kinnex's (+1..+18, ≥12 of 18 A) — the residual false-positive class is internal
  priming the filter does not catch. → Stage-1d follow-up (issue 10 addendum): a `--ip-rule kinnex`
  option.~~ **INVERTED — CORRECTED 2026-08-21, see the "IP-rule note" at the foot of this file.** The
  shipped rule is the *stricter* of the two by every measure of what it removes from a call set, the
  correct shipped window is **−9..+30 (40 nt)**, and `--ip-rule kinnex` is measured dead.
- Stratification (exploratory only): 3′-UTR trusted-novel sites 0.78 at ≥5 UMI / 0.56 at ≥20;
  ≥5-molecule sites 0.73–0.76 / 0.53–0.66; 2-molecule sites 0.43 / 0.16.
- Signed distances are asymmetric: Kinnex termini sit a few bp upstream of PeakATail cleavage points
  (1,228 hits at −5..0 vs 210 at 0..5) — consistent with the clip-site convention; no effect at 25 bp.

## What the paper says (reliability-first)

1. PeakATail's atlas-novel calls are strongly enriched for real cleavage sites (70–250× over null) but
   **the pre-registered "trusted de novo PAS" definition did not reach the 70% long-read target
   (52%)**; we report this as a negative result with the funnel figure.
2. The residual error is dominated by intronic sites and by internal priming the current filter misses.
   A stricter definition (3′-UTR-restricted and/or ≥5 molecules) reaches
   the target in exploratory stratification, but **was chosen after seeing the data**;
   *(the "and/or a Kinnex-style IP rule" option that stood here was removed 2026-08-21 — measured dead,
   see the IP-rule note below)*; it may be
   pre-registered as v2 and must be validated on data not used here (e.g. after the Stage-1d re-run, on
   a held-out truth set) before the paper calls any de novo site "trusted".
3. Until then, de novo sites are reported with their support (molecules), hexamer tier, feature class
   and IP status, and the paper makes no "trusted novel" claim.

Files: `results/reliability/trusted_novel_final_pbmc/` (REPORT_PROVISIONAL.md §1–9 incl. verifier
robustness table; call_primary/, call_sens_min1/, validate_*/). Tool: `scripts/reliability/trusted_novel_pas.py`.

## v2 section (2026-08-21, code `9dfdefb` — after the Stage-1d fixes; verified FIXED)

Re-running the pre-registered trusted-novel pipeline unchanged on the v2 PBMC precision default
(46,544 sites, code `9dfdefb3`) yields **7,259 trusted-novel PAS** (15.6%; 4,329 strong) whose Kinnex
x3p concordance at the pre-registered ≤25 bp / ≥5 UMI metric is **0.4853 [0.474–0.497]** — *below* the
v1 value of 0.5211 — with `meets_target_0.70 = 0` at every truth threshold and every robustness variant
(GEM-X 0.5154, pooled 0.5816 @25 bp) moving down, not up. An independent adversarial verifier re-derived
the input identity, the whole funnel, the concordances, the shuffled null and the decoy fractions from
the primary inputs and reproduced every number exactly, and confirmed that no threshold, window, seed,
atlas, truth set or `--incl` space differs from v1 — **so the negative result is not an artefact of the
Stage-1d IP bugs.**

| | v1 (4efeb125) | v2 (9dfdefb) |
|---|---:|---:|
| input (on contigs) | 44,394 | 46,524 |
| trusted-novel / strong | 6,628 / 4,037 | 7,259 / 4,329 |
| concordance ≥5 UMI | 0.5211 | 0.4853 |
| IP-decoy proximity | 24.5% | **27.0%** |
| atlas-known hexpass concordance | 0.9023 | 0.8941 |

Mechanism: the corrected strand-aware filter flags *fewer* raw peaks (110,000 vs 119,518), so more
A-rich-downstream sites survive. ~~the residual failure is the IP **rule's** leniency (issue #95
addendum: `--ip-rule kinnex`), not the strand bug.~~ **CORRECTED 2026-08-21:** the v1→v2 *direction* of
this movement stands, but the attribution to "rule leniency, fixable by adopting the Kinnex rule" does
not — see the IP-rule note below. What the residual decoy proximity shows is that a class of A-rich
loci survives **both** rules, not that the shipped rule is the more permissive one. The definitional conclusion in §"What the paper
says" is unchanged and now robust to the fix. Structural note for Methods: on the precision-default
input the clip/≥2-molecule/not-IP criteria remove 0 sites by construction — the trusted-novel metric
is carried entirely by hexamer ∧ atlas-novelty. Files:
`results/reliability/trusted_novel_final_v2_pbmc/` (REPORT_PROVISIONAL.md status VERIFIED;
`verifier_crosscheck/ADVERSARIAL_VERIFY_v2.txt`).


---

## IP-rule note — the "looser than Kinnex's" premise is INVERTED (dated correction, 2026-08-21)

**What this file said.** That the caller's internal-priming rule is *looser* than the Kinnex long-read
decoy rule, and that the residual internal-priming false positives should therefore be attacked with a
`--ip-rule kinnex` option (issue 10 §4 addendum; also echoed in [19](19_final_gate_v2.md) §4).

**What was measured** (`results/perf_gap/D2_rule_sweep/` §2(5) and `D2_arm_sweep.tsv` group 4;
adversarially verified, `results/perf_gap/VERIFY/README.md` correction 1, which re-derived the caller's
own per-peak flag on **402,860 / 402,860** called peaks with zero discrepancies):

| | shipped rule (`ip_window`, transcript-relative **−9..+30**, 40 nt; ≥6 consecutive A or A-fraction ≥0.70) | Kinnex decoy rule (**+1..+18**; ≥12 A of 18 or ≥6 consecutive A) |
|---|---|---|
| peaks flagged genome-wide | **68,855 of 402,860 (17.09 %)** | **12,186 (3.02 %)** |
| applied to the no-IP ≥2-molecule arm (n 62,110, P@100 0.5907) | removes **15,586** → n **46,524**, P **0.7062** *(this arm **is** the shipped precision default)* | removes **900** → n 61,210, P **0.5944** |
| Kinnex decoy proximity of the survivors | 0.3122 → **0.1300** | 0.3122 → 0.3040 |
| added on top of the shipped rule | — | 31 of 46,524 sites; P 0.7062 → 0.7065 (**no-op**) |

**The shipped rule is 5.65× stricter** (68,855 / 12,186), and at the level of what each rule actually
*removes from a call set* it is the far more aggressive of the two: it buys +11.6 pp of precision where
the Kinnex rule buys +0.4 pp. So the premise is inverted, and its two consequences are withdrawn:

1. **`--ip-rule kinnex` is measured dead** — a no-op as an addition, a large regression as a replacement
   (P 0.7062 → 0.5944, decoy proximity 0.130 → 0.304). Withdrawn in `manuscript/github/issue10_stage1d_ipfilter_memory.md` §4
   and in [19](19_final_gate_v2.md) §4. Do not resurrect it.
2. **The residual decoy proximity is not a rule-leniency artefact.** 27.0 % of v2 trusted-novel sites
   near a decoy terminus is a real, reproduced measurement; what is wrong is the diagnosis. The sites
   survive both genomic A-richness rules, so the fix — if there is one — is not another genomic window.

**Three bookkeeping corrections that come with it.** (a) The window quoted in this file as "−10..+30" is
the 41-nt off-by-one that the verifier found and fixed; the caller's window is **−9..+30 (40 nt)**,
anchored on BED `end` for `+` and BED `start` for `−`. (b) Do not quote "~17×" from
`results/perf_gap/VERIFY/README.md` — that ratio is an arithmetic slip in the verifier's own prose;
68,855 / 12,186 = **5.65** ([22](22_performance_roadmap.md) §6 item 5). (c) The decoy proximity of the
shipped rule's survivors is **0.1300**, the value in the verifier-corrected `D2_arm_sweep.tsv` row
`PN_sup2_callerrule`; the "~0.120" that still appears in `results/perf_gap/D2_rule_sweep/README.md` §2(5)
prose is the pre-correction 41-nt-window value (0.1204) and must not be quoted.

**What is unchanged.** Every concordance number, the funnel, the null, the decoy fractions and the
negative result itself (the pre-registered definition reaches 0.485 on v2 / 0.521 on v1 against a 0.70
target) — all reproduced exactly by the verifier. This note changes one *diagnosis*, not one measurement.
Competitive-framing context: [25](25_competitive_position.md) §7 (the default's residual error is
low-support intronic peaks without a poly(A) signal, **not** internal priming).
