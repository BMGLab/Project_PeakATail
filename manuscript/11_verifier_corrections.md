# Verifier corrections — spermatogenesis & Kinnex arms (2026-08-20)

Both adversarial verifiers returned **PROBLEM**: load-bearing numbers reproduce, but headline prose
overstates in ways that must not reach the manuscript. Corrections below are binding on any use of
`spermatogenesis_control.*` and `kinnex_truth_validation.*`.

## Spermatogenesis positive control — corrected claims

| As reported | Corrected |
|---|---|
| "Monotone shortening across all four stages, both mice" | **SPC>RS>ES is rock-solid** under every labelling tried (ρ −0.73..−0.86). The **SPG first step is fragile**: holds under shipped GEX labels (p≈1e-9) but INVERTS under independent marker-argmax labels; it rests on 64/68 cells. Claim the three-stage gradient; present SPG as supplementary. |
| "Trend lives ONLY in atlas/3'UTR-supported PAS; unsupported majority poisons it" | True by **peak count** (36–37% supported) but not by **signal**: 3'UTR peaks carry **70–75% of UMIs**; atlas-supported carry 65–67%. The sign-flippers are specifically **intronic/other-exon peaks** (relative position ~0.3 = upstream); read-through peaks past the 3'UTR **preserve** the trend. Both framings must appear. |
| (undisclosed) | The shipped trend TSV carries a second index, **wul (bp from proximal PAS), which contradicts the atlas-arm headline** (wul ρ +0.52 where wdi is −0.74/−0.48) and is weaker and non-monotone on the 3'UTR arm. **Unexplained — must be resolved or the panel C atlas arm removed before manuscript use.** |
| "resolves ≥2 3'UTR PAS for only ~1,000 genes" | Undercount ~2×: 1,876 / 2,190 genes at ≥50 UMI; ~1,000 is post-pseudobulk-depth-filter. |
| "peak intervals up to ~1.3 kb" | 1.3 kb is ~p97; true max 6.4 / 5.3 kb (median 325/302 bp). |
| mouse1 depth-matched cell "NaN, underpowered" | Verifier's own thinning yields ρ = **−0.756** on 1,189 cells — fill it in. |

Everything else replicated exactly (repair numbers, per-gene ρ 0.874, 6/6 literature genes, tool-own
PDUI trend, 200-perm null).

## Kinnex long-read truth — corrected claims

| As reported | Corrected |
|---|---|
| "PeakATail LAST of six in both long-read datasets" | Last **only at ≥5/≥20-UMI truth**; at ≥100/≥500 UMI PeakATail is **5th and scAPAtrap last**, both donors. The two are indistinguishable at the bottom; say that. |
| "precision ordering replicates, ρ=0.90 at every stringency" | 0.90 at 7 of 8 arms, 0.80 at one; donor-donor 1.00 at 3 of 4, 0.90 at one. Direction unchanged; drop "every". |
| "IP decoy flat at 3.6–5.2% everywhere" | Flat through 4,999 UMI, **rises to 14.5% in the 138 peaks >5,000** (the figure shows it; prose didn't). |
| "GEM-X replicates within ~4 points" | Within ~2 pts for the four upper tools; **PeakATail and scAPAtrap shift 5–12 pts** (same direction, ranks unchanged). |
| top-N self-rank decomposition | **Computed on the forbidden mis-keyed topN sets** (7.9% member overlap with correct ranking). Verifier re-ran on re-keyed sets: genuine-fraction 22.7/22.7/22.6% — conclusion unchanged, but only the re-keyed numbers may be quoted. |
| truth-set "UMI" support | Actually **alignment-record counts**: ~10.5% of records repeat a (CB,UMI) at the same terminus (isoseq groupdedup splits indel-variant records). Support thresholds are therefore slightly optimistic; state as record-count support or dedup before v2. |

## Status

Both figures remain on disk but are **quarantined from manuscript use** until regenerated with these
corrections — which will happen anyway in Stage 3 (both arms are re-run on re-keyed matrices after
Stage 0/1 merge). The corrections themselves are now the spec for that regeneration.
