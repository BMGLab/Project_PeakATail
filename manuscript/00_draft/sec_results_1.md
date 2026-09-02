# Results — Part 1: the caller and the benchmark (Figs 1–3)

<!-- Drafted 2026-09-02 from verified sources only: 19 §1–§6, 23 §3/§4, 25 §2–§8 (head-to-head
     one-liner is §8.1 with edits preserving every qualifier; single-point rule §8.2/§9; recall
     budget §8.4/§8.5), 26 §R1–R4/V, fig1/fig2/fig3 caption sidecars, figS1/figS8 sidecars,
     FIGURES_MANIFEST.md. Figure numbering is the frozen scheme. Every number below appears in one
     of those records. S5/S6 are PENDING and nothing is cited from them. -->

## Direct poly(A) evidence: a scarce, specific channel and one pre-registered operating point (Fig 1)

PeakATail calls poly(A) sites (PAS) from the one element of a droplet 3′-end library that observes
cleavage directly: the non-templated poly(A) tail carried in a read's 3′-side terminal soft clip. A
read is accepted when it is mapped on the strand of the pass, carries a 16-nt cell barcode and
aligns over at most the library read length; it enters the poly(A) channel when its 3′-side terminal
soft clip is ≥6 nt and ≥80% A (+) / T (−), with a ≥6-nt run flush to the alignment edge (Fig 1a;
per-library parameters in Fig S1). In the measured prototype the flush condition alone buys a
92-fold specificity against wrong-end clips — a control on alignment artefacts, not on genomic
A-runs. Support is counted in distinct (cell barcode, UMI) molecules, never in reads. The channel is scarce by construction: genome-wide, 0.573% of accepted
cell-barcoded reads carry a qualifying clip (3,195,067 of 557,564,408 on the PBMC 10k v3 library; an
earlier 1.152% figure was a head-sampling artefact of the QC estimator and is reported only as a
corrected estimate). The channel is destroyed outright by any pipeline that trims poly(A)
before alignment.

Fig 1c shows what the scarcity buys. On PBMC 10k v3 (CellRanger 3.0.0 [[CITE: CellRanger paper]]),
402,860 candidate peaks enter the funnel; the internal-priming filter removes 68,855 (17.09%); the
survivors split into tier 1 (clip-supported; 167,629) and tier 2 (coverage-only summits; 166,376);
and the pre-registered default keeps tier-1, internal-priming-pass calls with ≥2 distinct clip
molecules (46,544 rows; 46,524 scored). Atlas-agreement precision at 100 bp — agreement with a
strand-matched PolyASite 2.0 representative site [[CITE: PolyASite 2.0 paper]], not ground truth —
tracks evidence type, not peak geometry: tier 2 scores 0.0558 (n 166,355), tier 1 at ≥1 molecule
0.3520 (n 167,565, IP arm), and the default 0.7062, against a gene-body-shuffled null of 0.0217
(peak-level QC in Fig S2).

The internal-priming filter is measured on its own output rather than assumed (Fig 1d). On the PBMC
≥2-molecule arm it removes 15,588 of 62,132 calls (25.1%), moving atlas-agreement P@100 from 0.5907
(no-IP arm) to 0.7062 while costing detected-gene recall (R_det) 0.1911 → 0.1754; the removed
calls' own atlas-agreement precision is 0.246, and 99.97% of removals are triggered by the
≥6-consecutive-A rule. The default is one point on a measured curve (the sweep is Fig 3a), and only
the ≥2-molecule point was pre-registered.

## Accuracy against the field, under a gate committed before the run (Fig 2)

Every comparison in Fig 2 uses one scorer and identical denominators for every tool.
Atlas-agreement precision (P@100) is the fraction of called 1-bp points with a strand-matched
PolyASite 2.0 representative site within 100 bp (atlases of 569,005 sites on GRCh38 and 301,006 on
GRCm38); an atlas-novel true site counts as a false positive. Detected-gene recall (R_det) is the
fraction of the atlas restricted to genes detected in each dataset that is recovered within 100 bp —
denominators 285,136 sites (PBMC 10k v3) and 126,686 (GSE104556 testis; STARsolo
[[CITE: STARsolo paper]]; datasets, Fig S1). The precision-first default and its gate,
P@100 ≥ 0.50, were committed before the benchmark arms ran (2026-08-21: gate 01:19 in commit
0e27b1a, arms 02:59, v2 re-run 16:14) and passed on all three: 0.7062 (PBMC, n 46,524), 0.7450
(mouse 1, n 26,255) and 0.7572 (mouse 2, n 26,526), each 33–55× its three-seed gene-body-shuffled
null (null design, Fig S3). The −1.1 pp PBMC move from the v1 code is the corrected minus-strand
internal-priming window, reconstructed key-for-key (version lineage, Fig S12).

Against the de novo field — polyApipe [[CITE: polyApipe paper]], Sierra [[CITE: Sierra paper]],
scAPAtrap [[CITE: scAPAtrap paper]] and SCAPTURE [[CITE: SCAPTURE paper]], with catalog-based
scUTRquant [[CITE: scUTRquant paper]] shown but not ranked — call-set sizes span 35,759 to 787,138
(22×), and both precision and recall move with call budget, so a single-point comparison measures
the operating point, not the tool. (scTail [[CITE: scTail paper]] was not runnable on this BAM,
whose R1 is 28 bp; SCAPTURE's mouse-2 run was scored at site level, P@100 0.672.) At its pre-registered
operating point PeakATail is the most atlas-concordant de novo call set in the panel on both
datasets (atlas-agreement P@100 0.7062 PBMC; 0.7450 / 0.7572 mice). Its recall at that point (R_det
0.1754 / 0.2048 / 0.2080) is below polyApipe's, but that comparison sets 46,524 of our calls against
120,916 of polyApipe's: at matched call count PeakATail leads polyApipe on recall at every call
budget tested and on precision at every budget above 20,000 calls, on both datasets — at polyApipe's
own N of 120,916, 0.4036 / 0.2336 against 0.3800 / 0.1988 — and on both mice the ≥1-molecule arm
exceeds polyApipe on precision, recall and F1 simultaneously (0.5686 / 0.2806 / 0.3758 and
0.5880 / 0.2826 / 0.3818 against 0.4005 / 0.2499 / 0.3078 and 0.4119 / 0.2508 / 0.3118). The default
is a deliberately conservative point on a curve that dominates polyApipe's operating point, not the
tool's frontier. Three bounds travel with this claim: it is restricted to de novo tools
(catalog-based scUTRquant is above PeakATail on both axes at several matched budgets, and dropping
the restriction makes the claim false); the precision half carries its lower bound — below
N ≈ 15,000 polyApipe's atlas-agreement precision is the higher of the two — while the recall half
holds at every N; and SCAPTURE ships no non-circular ranking column, so it cannot be truncated to a
matched N.

The default is also corroborated by truth that owes nothing to the atlas: 0.7647 of its 46,524 PBMC
sites (35,575; Wilson 95% CI 0.7608–0.7685) lie within 25 bp, strand-matched, of a poly(A)-verified
Kinnex [[CITE: Kinnex paper]] x3p long-read 3′ end supported at ≥5 UMI, against a gene-body-shuffled
null of 0.0072 (10 seeds; 106.5×); polyApipe's concordance on the same truth is 0.3895, and the
atlas-known, hexamer-pass complement scores 0.8941, calibrating the ceiling of the metric. Two
caveats travel with every Kinnex number: the long-read truth is from different donors than the
short-read library (donor-mismatched), and its "≥5 UMI" support is an un-de-duplicated
alignment-record count, so support thresholds are slightly optimistic.

Fig 2d decomposes where the default's stringency acts: 72.2% of tier-1 sites on the PBMC IP arm are
single-molecule (121,041 of 167,565; 72.1% on the no-IP arm; mice 50.3% / 48.8%), so the
≥2-molecule rule does most of the selection; Fig 3 prices what it costs.

## The trade surface: resolution, reproducibility across libraries, and compute (Fig 3)

Sweeping the molecule-support threshold on the IP-filtered v2 arms from ≥1 to ≥10 moves
atlas-agreement P@100 from 0.352 to 0.941 on PBMC while R_det falls 0.268 → 0.083; the mice reach
0.885 / 0.892 at ≥10 (Fig 3a). Only the ≥2-molecule point was pre-registered; the others are
descriptive, and ≥5 or ≥10 must not be read as a recommendation. The ≥2-molecule default is
not the F1 optimum on any dataset — the ≥1-molecule arm reaches F1_det 0.3046 (PBMC) and
0.3758 / 0.3818 (mice) against the default's 0.2811 and 0.3213 / 0.3264. The default is chosen for
reliability, not for F1: it is the arm whose calls a biologist can act on without re-validating each
site, and we report both arms as labelled operating points.

Tightening the matching window to 10 bp does not change the de novo precision ordering (Fig 3b). The
default retains P@10 / P@100 = 0.74 and polyApipe 0.77, against 0.46 (SCAPTURE), 0.33 (Sierra) and
0.24 (scAPAtrap): the precision lead is not an artefact of a loose window; catalog-based scUTRquant
(retention 0.72) sits above the default at every window, again shown, not ranked.

The two testis mice serve as biological replicates (Fig 3c): 0.776–0.779 of one mouse's default
calls have a strand-matched call in the other within 100 bp (0.729–0.737 at 25 bp), against a chance
level ≤0.010, and the ≥2-molecule threshold buys 14.3 / 12.5 points of replicate agreement over the
≥1-molecule arm (0.636–0.651). The honest scope: the v2 default is above polyApipe and PeakATail's own
shipped caller but below Sierra (0.790 / 0.834) and scAPAtrap (0.884 / 0.793 — a set its own
reducePeaks depth-cleaning step has already filtered), and the competitor rows are v1-era
measurements. PBMC is a single donor and a
single BAM, so the human replicate had to be created.

We therefore pre-registered a second-library generalisation test before any number from it existed
(Fig S4): the frozen v2 tree, run unchanged on pbmc4k — a 10x 3′ v2-chemistry, CellRanger 2.1.0,
98-bp library — with 107 of 108 recorded parameters identical and only the read-length descriptor
moving (91 → 98). The unchanged gate passed — atlas-agreement P@100 0.8279 (n 20,672) — and now
holds on four independent libraries (0.7062 / 0.8279 / 0.7450 / 0.7572). One caveat travels with that number: precision generalises, it does not improve.
The fixed rule simply lands at a more conservative point on the shallower v2 library — donor 1's
default, cut to 20,672 calls by the tool's own clip-molecule ranking, scores 0.886–0.907 and is
above donor 2 on both axes — and donor 2's R_det is lower (0.1104 vs 0.1754; denominators 268,097
and 285,136 sites; the gap survives every denominator choice).
Cross-donor concordance is the human analogue of Fig 3c: 84.2% of donor 2's default calls are
reproduced in donor 1 within 100 bp (~361× the genic-shuffle null); the reverse direction is 39.9%
against its 44.4% arithmetic ceiling — donor 1 emits 46,524 default calls to donor 2's 20,672 — and
both directions are reported, neither alone. 10x publishes no donor identifier for pbmc4k: this is
demonstrably a second library, chemistry and CellRanger version, only presumptively a second
individual.

The same code revision that produced the v2 numbers removed the tool's hardware barrier (Fig 3d;
detail in Fig S8). On the PBMC BAM, peak RSS fell 293.7 GB → 12.53 GB
(23.5×) — the ≥300 GB node requirement is gone — and wall time 3:45:53 → 34:37 on the like-for-like
pair (the clip-seeded arm without the internal-priming filter); the v2 wall time was measured with
all four benchmark arms running concurrently (~28–35 min; uncontended single run, 27:43). At cohort
scale — 17 libraries from 14 patients, ~224 GB of BAM — the
identical configuration fell from 9:06:26 to 1:06:26 (8.2×; v2 run with --peak-workers 16) at
identical output, 505,197 unified PAS from both codes, with peak RSS 13.63 → 13.53 GB.

**The recall budget, measured.** We measured how much of the recall gap any configuration of this
caller could ever recover. Of the 235,111
detected-gene atlas sites the PBMC default misses, 79.96% have no peak of any kind within 100 bp
(84.53% in mouse 1), so at most 20.04% / 15.47% of the gap is reachable by any threshold, filter or
tier policy; keeping every peak the caller ever produced reaches R_det 0.3406 (PBMC) and 0.3278
(mouse 1) — the hard ceiling for anything downstream of peak calling. The denominator is itself an
upper bound on what any tool can find in one library: only 44.6% of the detected-gene atlas is
corroborated by matched-tissue Kinnex long reads at ≥5 UMI (the union count, 127,188 of 285,136) —
long reads from a different donor, so "not corroborated" is not "not real" — and against the
x3p-corroborated subset (107,260 sites) the same default recalls 0.4372, not 0.1754. Those are facts
about two different subsets and must not be divided across each other. The atlas places 20.0 sites
in an average detected human gene where the default emits 4.0: absolute recall against such a
denominator should not be read as sensitivity.
