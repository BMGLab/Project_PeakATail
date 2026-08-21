# Final Stage-2 run — verified gate results (code `4efeb125`, run 2026-08-21 02:59–06:46)

**Verifier verdict: SOUND.** Every number below was reproduced by an independent verifier: the
arm score TSVs are byte-identical to a re-run of `score_tool.py`, and an independent
`bedtools closest -s -d` / `bedtools window -sm` pipeline reproduced every P@100 and R_det@100 to
six decimals (max discrepancy 0.000000). Provenance: all four arms ran from the frozen worktree
`tools/pa-polya-run-4efeb125` (HEAD `4efeb1252e6d…`, clean), flags per arm confirmed from
`run_config.json`, `LC_ALL=C` in effect (all BEDs pass `sort -c` under C locale).

**Pre-registration timing (verified from git):** the precision-first default and its gate were
committed in `0e27b1a` at 2026-08-21 **01:18:59**; the arms started at **02:59:36**. The ≥2-molecule
threshold was chosen after seeing a stale Stage-2 sweep ([12](12_stage2_gate.md) CORRECTION) and
pre-registered before the final run — the paper discloses both facts.

## 1. Pre-registered precision-first default — gate P@100 ≥ 0.50 on PBMC AND both mice

Default output = clip-seeded tier-1 ∩ IP-pass (`--ip-filter --ip-filter-mode filter`) ∩ ≥2 distinct
clip molecules.

| dataset | n scored | **P@100** | R_det@100 (denominator) | F1_det@100 | null_genic P@100 (3 seeds) | gate |
|---|---:|---:|---:|---:|---|---|
| PBMC 10k v3 (IP arm) | 44,394 | **0.7167** | 0.1707 (285,136) | 0.2757 | 0.0226 / 0.0218 / 0.0218 | **PASS** |
| testis mouse 1 | 25,991 | **0.7414** | 0.2021 (126,686) | 0.3176 | 0.0146 / 0.0135 / 0.0143 | **PASS** |
| testis mouse 2 | 26,164 | **0.7554** | 0.2049 (126,686) | 0.3223 | 0.0154 / 0.0134 / 0.0142 | **PASS** |

Gate **passes on all three call sets**. Null bars are 30–50× below the real values. Full-atlas recall
(not gated, must be shown alongside): 0.0865 (PBMC), 0.0963 / 0.0980 (mice).

## 2. Original two-sided gate (F1_det > 0.261 AND P@100 ≥ 0.38, PBMC ≥1-molecule output)

| PBMC output | n | P@100 | R_det | F1_det | verdict |
|---|---:|---:|---:|---:|---|
| default arm (no IP), both tiers = [12]'s "default output" | 402,765 | 0.1932 | 0.3363 | 0.2455 | FAIL (F1 and P) |
| default arm (no IP), tier-1 ≥1 mol | 222,955 | 0.3032 | 0.2970 | 0.3001 | F1 pass / **P FAIL** |
| IP arm, tier-1 ≥1 mol (IP-on variant) | 161,595 | 0.3546 | 0.2619 | 0.3013 | F1 pass / **P FAIL** |
| IP arm, both tiers | 328,302 | 0.2030 | 0.2979 | 0.2415 | FAIL |

Unchanged from [12](12_stage2_gate.md): **no ≥1-molecule output clears the 0.38 precision floor.**
The ≥1-molecule tier-1 outputs are reported as the *sensitivity arm* (highest de novo F1, 0.300–0.301).

Other rows (verified): PBMC tier-2 P 0.057 / F1 0.047 (no IP) and 0.056 / 0.044 (IP); mouse both-tiers
P 0.443 / 0.448, tier-1 ≥1 mol P 0.561 / 0.581 (F1 0.371 / 0.377), tier-2 P 0.119 / 0.114.
**Mislabel warning:** the launcher also wrote a `pas_PRESPEC_precision_default.bed` in the *non-IP*
PBMC default arm (n 62,110, P 0.5907, R 0.1911, F1 0.2888). It is NOT the pre-registered default
(no IP filter); it has been renamed `pas_tier1_ge2mol_noIP_POSTHOC.bed` and the launcher fixed.

## 3. Head-to-head (same scorer, 100 bp strand-matched, same denominators)

PBMC 10k v3 — n / P@100 / R_det / F1_det. scUTRquant is catalog-based (shown, not ranked with de novo);
scTail has no number (R1 = 28 bp; BLOCKED).

| tool | n | P@100 | R_det | F1_det |
|---|---:|---:|---:|---:|
| scUTRquant* (catalog) | 40,519 | 0.793 | 0.178 | 0.290 |
| **PeakATail precision default** | 44,394 | **0.717** | 0.171 | 0.276 |
| SCAPTURE | 35,759 | 0.652 | 0.118 | 0.199 |
| polyApipe | 120,916 | 0.380 | 0.199 | 0.261 |
| PeakATail tier-1 ≥1 mol, IP | 161,595 | 0.355 | 0.262 | 0.301 |
| PeakATail tier-1 ≥1 mol, no IP | 222,955 | 0.303 | 0.297 | 0.300 |
| Sierra | 106,170 | 0.256 | 0.135 | 0.177 |
| PeakATail shipped (pre-fix) | 277,111 | 0.118 | 0.156 | 0.134 |
| scAPAtrap | 787,138 | 0.100 | 0.300 | 0.150 |

Mouse testis m1 / m2 — P@100 / R_det / F1_det: scUTRquant* 0.782/0.784 · 0.259/0.259 · 0.389/0.389;
**PeakATail precision default 0.741/0.755 · 0.202/0.205 · 0.318/0.322**; SCAPTURE 0.694/0.672 ·
0.155/0.147 · 0.253/0.241; PeakATail tier-1 ≥1 mol IP 0.561/0.581 · 0.278/0.279 · 0.371/0.377;
polyApipe 0.401/0.412 · 0.250/0.251 · 0.308/0.312; shipped 0.369/0.368 · 0.198/0.201 · 0.258/0.260;
scAPAtrap 0.250/0.239 · 0.245/0.259 · 0.248/0.248; Sierra 0.509/0.581 · 0.107/0.117 · 0.177/0.195.

**Honest one-liner (verifier's wording):** the precision default is the most atlas-concordant de novo
call set in the panel on both datasets, at polyApipe-or-lower recall (PBMC −14%, mouse −19% vs
polyApipe) and essentially tied F1 (+0.015 / +0.010). The F1 lead is not a meaningful margin; the
claim is precision, not overall accuracy. The IP filter alone buys +4.9 pp precision for −1.4 pp
recall (mouse 1, ≥2-mol, no-IP 0.692 → IP 0.741).

## 4. Molecule unit (from `pas_support.tsv`)

72.5% of PBMC tier-1 sites are single-molecule and are discarded by the default (mouse ≈ 50%). Of
the ≥2-molecule PBMC sites, 90.2% keep ≥2 molecules under a `-F 3844` alignment filter (mouse
94.8 / 95.6%), so the stricter variant would shrink the default by ~10% on PBMC.

## 5. Compute (limitation — must be stated; FIXED on perf/clip-memory, see addendum below)

| arm | wall | peak RSS | CPU |
|---|---|---|---|
| PBMC default (no IP), `--threads 16` | 3:45:53 | 293.7 GB (280.1 GiB) | 143% |
| PBMC IP arm | 3:44:37 | 239.5 GB (228.4 GiB) | 146% |
| mouse 1 / mouse 2 (IP, `--threads 12`) | 1:01:09 / 1:01:48 | 23.1 / 22.4 GB | ~140% |

2.1× the Stage-2 run (140 GB) and 2.8× the shipped caller (105.6 GB) on PBMC; uses ~1.4 cores despite
`--threads 16`. Heaviest tool in the panel by a wide margin (polyApipe 3:26 h, Sierra 2:43 h,
scUTRquant 0:30 h, none near 100 GB). A 10k-cell PBMC BAM needs a ≥300 GB node. The
[10](10_caller_fix_plan.md) stop signal (RSS > 150 GB → profile first) was exceeded; profiling is
a Stage-1d item for Amir.

**Addendum 2026-08-21 12:30 — compute fixed, outputs unchanged.** Profiling (results/perf/) showed the peak RSS is
not in peak calling but in `_tfidf_signac_method1` (clustering): four dense float64 copies of the cells × PAS matrix
(4 × 23,303 × 390,493 × 8 B = 291 GB; measured 287.8 GB in isolation). The 1b/1c "doubling" was indirect (2.07×
larger matrix area). Branch `perf/clip-memory` (verified FIXED, byte-identical on the full PBMC, mouse1 and slice
runs — 18/18 files incl. pasbed.bed, pas_support.tsv, matrices, clusters): sparse TF-IDF, per-(contig, strand)
parallel peak calling with deterministic merge, read_check reorder, streamed CB filter → PBMC **12.45 GB / 27m43s**
(process tree 19.3 GB), mouse1 3.68 GB / 9m03s. The paper reports the fixed footprint once the PR is merged and the
v2 re-run confirms it; the numbers in §1–3 are unaffected by construction.

## 6. Tool defects surfaced (Stage 1d; issue draft [github/issue10_stage1d_ipfilter_memory.md](github/issue10_stage1d_ipfilter_memory.md))

1. **IP filter minus-strand window.** `ema/experimental/internal_priming.py` applies the same
   forward-coordinate window [pos−10, pos+30) to both strands, so on "−" it tests 30 bp *upstream* /
   10 bp downstream in transcript orientation. A strand-symmetric re-application flags 3.5–4.2% of
   surviving minus-strand sites; benchmark effect negligible (PBMC default P 0.7167 → 0.7182 if
   removed), but the "IP-pass" criterion of the trusted-novel definition must recompute IP
   strand-correctly rather than trust the run's flag.
2. In `filter` mode `annotatedpas.bed` carries no per-site IP flag (only counts in
   `peak_filters_stats.json`); an explicit flag file is needed downstream.
3. `pasbed.bed` silently drops gene-assigned PAS with zero counts after the cell filter
   (0.45% PBMC, 3.7% mouse 1) — the scored set is the post-cell-filter set; say so in Methods.
4. Cosmetic: "seq-len not set; defaulting to 150" logged although 91/98 was resolved.

## 7. Overstatement guardrails for the paper

- Call it **atlas-agreement precision at 100 bp**, not precision against ground truth; atlas-novel
  true sites count as false positives here. The Kinnex long-read panel is the atlas-independent check
  (running on the PBMC default output; see 16 when it lands).
- Always show full-atlas recall next to detected-gene recall.
- Single PBMC donor / single CellRanger BAM; the two mice are one study and chemistry. The gate was met
  with a large margin but on three call sets from two datasets.
- Disclose: default threshold pre-registered before the final run, chosen after a stale run.
