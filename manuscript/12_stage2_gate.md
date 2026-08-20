# Stage 2 gate — full PBMC re-run of the clip-seeded caller (2026-08-20)

Pre-registered gate (10_caller_fix_plan.md §Stage 2): **F1@100 > 0.261 (polyApipe) on full pbmc_10k_v3**,
detected-gene denominator (285,136 sites), tiers reported separately. Numbers below are self-scored by
the launcher; an adversarial verification is running and must land before any of this is quoted.

| arm (pbmc_10k_v3) | n called | P@100 | R_det | **F1@100** | gate |
|---|---:|---:|---:|---:|:--:|
| shipped lambda_gradient | 277,111 | 0.118 | 0.156 | 0.134 | — |
| **clip_seeded, default output (both tiers)** | 377,236 | 0.188 | 0.314 | **0.235** | **FAIL** |
| clip_seeded **tier-1** (clip-supported) | 197,410 | 0.308 | 0.274 | **0.290** | PASS |
| clip_seeded tier-2 (coverage-only) | 179,826 | 0.057 | 0.040 | 0.047 | — |
| clip_seeded + IP filter, tier-1 | 142,400 | 0.362 | 0.241 | **0.290** | PASS |
| *polyApipe (best de novo, default output)* | *120,916* | *0.380* | *0.199* | *0.261* | *ref* |

Runtime 3:03 h, peak RSS 140 GB (plain) / 127 GB (IP) — still the heaviest tool in the panel.

## Honest reading

1. **As shipped by default, the fix does not clear the gate**: 0.235 < 0.261 (though 1.75× the old caller).
   The coverage-only tier (F1 0.047, P 0.057 vs genic null ~0.022) drags the combined number down.
2. **The clip-supported tier clears it**: 0.290 > 0.261 — and it is *not* the "Phase 1 trap" the plan
   warned about (a precision filter that sacrifices recall): tier-1 has **higher recall than polyApipe**
   (0.274 vs 0.199) on a **larger** call set (197k vs 121k). Adding the IP filter trades recall for
   precision at identical F1 (P 0.362 ≈ polyApipe's 0.380).
3. Therefore the decision is a **product decision, not a statistical one**: make tier-1 the default
   output (tier-2 written to a separate, clearly labelled low-confidence file) and the default output
   passes; keep both tiers merged and it fails. Either way the manuscript must disclose that the default
   was chosen *after* seeing these numbers, and must always report both tiers with n.
4. Testis replicates the direction (F1 0.258/0.260 → 0.296/0.307; tier-1 0.305/0.320).

## What is NOT resolved by this run

- Count matrices from these runs are invalid (Stage 1b tier-1 quantification fix in progress); all
  count-based re-tests (novelty, FDR) wait for the matrix regeneration.
- Resource footprint grew (105 → 140 GB RSS); acceptable for the box, a caveat for users.
