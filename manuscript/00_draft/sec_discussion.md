<!-- DRAFT: Discussion + Conclusions. Status: DRAFT for PI review, 2026-09-02.
Sources quoted: 14 (calibration), 16 §v2 (trusted-novel/Kinnex), 19 §1 (accuracy record), 20 §v2
(replication record), 23 (algorithm roadmap, corrected clip rate, scoring-model decomposition),
25 (competitive position — §8 wording preserved with every qualifier), 26 (second donor),
27 (singleton pre-registration), figS11_gatehistory.caption.md (timeline examples),
github/PUBLISHED.md (issues #98/#99). Figure refs use the frozen numbering. No number here is new. -->

# Discussion

What a precision-first default buys a biologist is trust per call. At its pre-registered
operating point — clip-supported tier-1 sites, internal-priming-filtered, ≥2 distinct clip
molecules — PeakATail's output carries atlas-agreement precision at 100 bp of 0.7062 (PBMC 10k v3)
and 0.7450/0.7572 (two testis mice) against curated PolyASite 2.0 representative sites
[[CITE: PolyASite 2.0 paper]], and 76.5% of the PBMC default's sites lie within 25 bp of a
poly(A)-verified Kinnex [[CITE: Kinnex paper]] long-read 3′ end at ≥5 UMI (Fig. 2). Under the
two-truth hard false-positive definition (no atlas site within 100 bp, no ≥5-UMI long-read end within
25 bp), 15.10% of the default's PBMC calls fail both truths, against 49.50% of polyApipe's
[[CITE: polyApipe paper]] (46,524 versus 120,916 calls; the matched-budget comparison follows). The
≥2-molecule default is not the F1 optimum on either dataset — against the same detected-gene
PolyASite denominators the ≥1-molecule sensitivity arm reaches F1_det 0.3046 (PBMC) and
0.3758/0.3818 (mice) versus the default's 0.2811 and 0.3213/0.3264. The default is chosen for
reliability, not for F1 — it is the arm whose calls a biologist can act on without re-validating
each site — and both arms ship as labelled operating points.

The cost is recall. At the default point, detected-gene recall (R_det; denominators 285,136
detected-gene atlas sites for PBMC, 126,686 for mouse) is 0.1754/0.2048/0.2080, below polyApipe's at
its own point — but that comparison sets 46,524 of our calls against 120,916 of polyApipe's, and both
axes move with call budget. At matched call count, among de novo tools, PeakATail leads polyApipe on
R_det at every budget tested and on atlas-agreement precision at every budget from 20,320 calls
upward, on both datasets (as-shipped rankings; the precision half is tie-sensitive inside
20,320 ≤ N ≤ 35,759 and reverses below N ≈ 15,000). At
polyApipe's own N of 120,916 the comparison is 0.4036/0.2336 against 0.3800/0.1988, and on both mice
the ≥1-molecule arm exceeds polyApipe on atlas-agreement precision, R_det and F1_det simultaneously
(Fig. 3). The default is a deliberately conservative point on a curve that dominates polyApipe's
operating point, not the tool's frontier. The de novo restriction is load-bearing: catalog-based
scUTRquant [[CITE: scUTRquant paper]] sits above every de novo arm at several matched budgets and is
shown unranked.

Both halves of this trade run into measured ceilings. The evidence channel is small: 0.573% of
accepted cell-barcoded reads in the PBMC 10k v3 BAM carry a qualifying poly(A) soft clip (correcting
our earlier documented 1.152%). 79.96% of the detected-gene atlas sites the PBMC default
misses have no peak of any kind within 100 bp (84.53% in mouse 1), so at most 20.04%/15.47% of the
recall gap is reachable by any threshold, filter or tier policy, and the hard ceiling for anything
downstream of peak calling is R_det 0.3406/0.3278. The denominator is itself an upper bound: only
44.6% of the detected-gene PBMC atlas is corroborated by Kinnex long reads at ≥5 UMI (the union of
both chemistries; the long reads are from a different donor, so "not corroborated" is not
"not real"), and — a separate fact — against the x3p-corroborated subset the same default
recalls 0.4372. Absolute recall against a pan-tissue atlas should not be read as sensitivity.

The negative results carry the field-level lessons. First, curated-atlas benchmarking
rewards atlas-shaped priors: an atlas-trained per-site scoring model gained 16 points of
atlas-agreement precision over the shipped ≥2-molecule rule while losing 8–11% on every long-read
axis, and the sites it swapped in were ~4× better by the atlas and ~2× worse by long reads. A method
fitted to the truth set improves against the truth set while degrading against the biology; every
gate here therefore pairs atlas agreement with an atlas-independent long-read check, and the field's
benchmarks should too. Second, pre-registration repeatedly caught our own overclaims; two examples
from Fig. S11's timeline: the sweep that informed the ≥2-molecule threshold was created twelve
minutes after the stale Stage-2 run failed the original two-sided gate on every arm — we withdrew our
earlier "pre-identified in the plan" description and disclose it as post hoc, informing a choice,
never cited as a result; and the adoption criterion for the revised defaults returned FAIL at
Δ = 0.000000 because it was drafted against a default configuration that never existed — the outputs
were byte-identical, an error only a pre-written criterion could expose. The same
discipline produced the paper's other negatives: only one of six switch-test configurations controls
the false-discovery rate, and the pre-registered trusted-novel definition reached 48.5% long-read
concordance against its 70% target — reported as a failure (Fig. S7).

Limitations. The human results rest on one donor at depth; the second library (pbmc4k) — pre-registered
default gate PASS at atlas-agreement P@100 0.8279, cross-donor site concordance 84.2% within 100 bp
(descriptive, not gated; Fig. S4) — carries no published donor identifier and is presumptively, not
provably, a second individual. The two testis mice come from one study and one chemistry. The Kinnex
truth is donor-mismatched and site-level: it corroborates positions, not donor-specific usage. Two
defects are open publicly: PAS-to-gene assignment in overlapping loci is under revision (issue #99),
so cohort switches are keyed by PAS identifier and no ranked gene list appears here; and the
per_isoform switch mode has a degenerate-pair defect (issue #98) and underlies no result here.

Future work follows the measured roadmap. Cross-sample "cohort borrowing" is the one recall
mechanism with verified headroom: mouse-1 single-clip-molecule tier-1 sites corroborated by a mouse-2
default call within 10 bp score atlas-agreement P@100 0.6683 versus 0.4144 for expression-matched
singletons — but internal priming replicates across samples too, the effect is unmeasured on human
data, and it ships only behind its own pre-registration. The 3′UTR singleton promotion is
pre-registered — rule, thresholds and falsifiers fixed, flag default OFF — and untested: its deciding
evidence is reserved to data unused in its discovery, and its ceiling is ~+0.04 R_det. Re-ranking is
nearly exhausted — a long-read-trained score gains +14.3% relative R_det at matched atlas-agreement
precision on the one library it was fitted on, but adopting it would replace the pre-registered
default and so requires its own pre-registration and multi-dataset runs; the next real gain lies in
the candidate generator or in the chemistry, not in thresholds.

# Conclusions

PeakATail calls poly(A) sites from direct poly(A)-tail read evidence and reports, for every claim,
the gate it was required to pass, the null it was tested against, and whether it failed. Under that
discipline its precision-first default is the most atlas-concordant de novo call set in the benchmark
on both datasets, corroborated by donor-independent long reads, and its calibrated,
replication-filtered switch test yields cell-type APA switches that recur across patients
(replication primary: 15 libraries from 12 patients) with none in 10 label-shuffle nulls (empirical
p ≤ 0.091, the 10-permutation floor). Just as deliberately, it reports a failed pre-registered gate,
a withdrawn clustering claim, anti-conservative shipped defaults and a literature panel that does not
reproduce. We offer the method — and the evaluation discipline it was built to survive — as the
contribution.
