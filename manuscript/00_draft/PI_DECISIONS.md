<!--
PI_DECISIONS.md — decisions and new-work items arising from the number audit (NUMBERS_LEDGER.tsv)
and the simulated referee report (REVIEW_SIMULATION.md), assembled during the 2026-09-02 fix pass.
Everything the evidence base could support was fixed directly in MAIN.md / SUPPLEMENTARY.md (see
STATUS.md for the list); this file holds only what needs new computation, new data, or a PI call.
Numbering follows REVIEW_SIMULATION.md (M# = major, m# = minor).
-->

# PI decisions and new-work queue

## A. Required before submission (reviewer-flagged as blocking; all need new work)

### A1. Panel-wide, de-duplicated long-read concordance (M1 + M8 — the review's single highest-value change)
**Problem.** The headline "most atlas-concordant de novo call set" rests on atlas agreement — the metric Fig S12 proves gameable — and the atlas-independent Kinnex read-out exists for only two tools (PeakATail 0.7647, polyApipe 0.3895). Also, the truth set counts un-deduplicated alignment records (~10.5% CB/UMI repeats), which is exactly what S5 is pending on, and the x3p/GEM-X datasets have no accession in Methods (a [[CITE]] slot for the data citation was added in the fix pass; the accession must come from the truth-build record).
**Options.**
1. De-duplicate the truth set (the S5 blocker), recompute the four inheriting numbers (0.7647, 0.3895, 15.10%/49.50%, S7's 48.5%), and extend concordance to Sierra, scAPAtrap, SCAPTURE and scUTRquant at shipped and matched call counts, presented as a co-primary axis in Fig 2 — this also completes S5 and answers the review's top request in full.
2. De-duplicate and recompute the existing numbers only; defer the panel extension with an explicit limitation sentence.
3. Submit as-is, with the un-deduplicated caveat that already travels with every Kinnex number.
**Recommendation:** Option 1. The reviewer is right that the paper's own S12 argument demands it, and the existing 0.7647-vs-0.3895 fragment suggests the claim strengthens. One de-duplication unblocks M1, M8 and S5 simultaneously. Fill the Kinnex data accession into the new [[CITE]] slot from the truth-build script's record.

### A2. Complete Supplementary Table T7 and resolve the S5/S6 slots (M10)
**Problem.** T7 (competitor versions/parameters/run commands) vouches for every competitor claim and does not exist; S5/S6 are "in preparation". Competitor version numbers are not in the quotable evidence base and must be read from the run logs (never reconstructed from memory).
**Options.** (1) Collate T7 from the verified run registry (the SUPPLEMENTARY note says content exists and needs only collation); build S5 via A1; build S6's shuffled-position null or drop the slot and renumber in one pass. (2) Drop S5/S6 now and renumber.
**Recommendation:** Collate T7 immediately (no science, only collation). Tie S5 to A1. For S6, decide by freeze date: the draft cites neither, so dropping costs nothing textually.

### A3. Make the companion repository public / mint the Zenodo DOI (M9vi)
Reviewers will expect working access at review time. The Zenodo placeholder already says the repository must be public first. This is a PI/maintainer action (known blockers are recorded in the project memory); no draft change needed beyond filling the placeholder.

## B. Strongly requested by the reviewer (new compute; feasible; PI to schedule or decline with a response letter)

### B1. More cohort null permutations (M3)
Ten permutations floor the empirical p at 0.091 under a zero-replication headline. The reviewer suggests ~100 permutations (their ~5-day estimate uses the 1:06 h cohort peak-calling time; the null pipeline re-runs only the switch/replication stage, so the true cost is likely lower and should be measured on one permutation first), or an analytical bound from the per-permutation single-patient call distribution (mean 1.1, already in Fig 6a).
**Recommendation:** Run 90 more permutations if the per-permutation cost measures ≤ ~2 h; otherwise add the analytical-bound supplement. Either satisfies the concern; the empirical route is cleaner to defend.

### B2. Competitor test through the calibration harness (M4)
Run Sierra's DEXSeq-based test (and optionally one more) through the identical label-permutation harness on testis mouse 1. Without it, the draft now states explicitly (fix pass) that other tools' calibration is an untested hypothesis. A second calibration dataset (e.g. PBMC labels, smaller effects) is the same reviewer's secondary request.
**Recommendation:** Do the Sierra/DEXSeq run — it is one tool, one dataset, 21 runs, and it converts the paper's field-level framing from scoped statement to demonstrated result whichever way it comes out. The second calibration dataset is worthwhile but severable.

### B3. pbmc_10k_v3 downsampled to pbmc4k depth (M5)
Separates depth from library/donor effects in the cross-donor asymmetry; ~half an hour by the paper's own timings, plus scoring.
**Recommendation:** Run it; report in Fig S4's audit TSV / legend. Low cost, closes a clean confound.

### B4. 200-shuffle null for the spermatogenesis PDUI arm (M6a)
Already a named pending item in the fig5 sidecar ("not covered by this run"). The fix pass added the 1/21-floor caution; 200 shuffles would make the z-tail claim empirical.
**Recommendation:** Run when the testis tree is next touched; not blocking if B1 is done (the cohort is where the zero-replication headline lives).

### B5. Manual PAS-to-gene check for the 16 literature-panel genes (M6c)
The fix pass added the caveat that mis-assignment (issue #99) is not excluded for individual panel verdicts. The reviewer calls the manual check "a one-afternoon" job.
**Recommendation:** Do it; whatever it shows, the panel paragraph then names the alternative explanation *tested*, matching the paper's own standard.

### B6. Fig 6e / gene-keyed robustness on the non-overlapping-locus subset (M7)
The fix pass marked 2,883 genes / 10,415 combinations / 57.8% as assignment-dependent. A recompute on the non-overlapping-locus subset would quantify robustness.
**Recommendation:** Cheap and worth doing together with B5; otherwise the added caution stands.

### B7. Tool-independent R_det denominator rebuild (M2 residual)
The fix pass documented the provenance (PeakATail `lambda_gradient` run) and the measured +1.26% construction sensitivity that biases *against* PeakATail. The reviewer's full request — rebuild from the CellRanger/STARsolo gene-count matrix at a stated expression threshold and show the ordering unchanged — is a small script plus a re-score.
**Recommendation:** Do it as a supplement row; the measured control already in Methods makes a surprise very unlikely.

## C. PI calls on wording/decisions (no new compute)

### C1. Abstract R_det sentence — applied, please confirm
The audit flagged that "R_det 0.110–0.208 ... below polyApipe's" folded pbmc4k (never run for polyApipe) into the comparison. The fix pass split it: "0.175–0.208 on the three benchmarked libraries — below polyApipe's — and 0.110 on the second donor", with denominators named. Confirm the abstract wording (it is ~15 words longer; the abstract stays well under GB limits).

### C2. Matched-N precision bound now reads "from 20,320 calls upward" everywhere
The draft carried both "above 20,000 calls" (25 §8.1's compression) and "from 20,320 upward" (25 §8.7's mandated form). The fix pass standardised on 20,320 in the abstract, key points, Results one-liner and Fig 2 legend — an edit 25 licenses ("an edit that preserves every qualifier"). Confirm.

### C3. Curated Laughney label provenance (M9iii)
Methods says cells carry "curated immune labels" confirmed by CellTypist; the provenance of the curation (Laughney atlas metadata vs in-house) is not in the quotable evidence base. One sentence from the analysis record (LABEL_POLICY.md / labels/confirmed_labels.tsv history) is needed.

### C4. Testis stage-label marker panel list (M9ii)
The panel genes behind the marker-argmax stage labels are not in the quotable evidence base; Methods now points to "panel and assignment script in the companion repository". Either leave the pointer (fine once A3 makes the repo public) or paste the panel into Methods from the script.

### C5. Binomial CIs on the four gate P@100 values (m7)
The Kinnex number has a Wilson CI; the gate values (0.7062/0.8279/0.7450/0.7572) do not. Computing CIs is trivial but produces numbers not in the evidence base, so it was not done in the fix pass. Recommend: compute Wilson CIs in the gate-table script (so they land in an audit TSV) and add to Fig 2/S4.

### C6. pbmc4k donor-identity genotype check (m8)
SNP genotype concordance between the two public BAMs would settle "presumptively a second individual" either way. Recommend doing it — the reviewer notes either outcome is informative; if declined, the current wording already survives.

### C7. Fig 3c v1-era competitor rows marked in-panel (m9)
Needs figure regeneration (greyed rows or asterisk) — a figure-script change, out of scope for the text-only fix pass. Recommend at next figure rebuild; the caption already discloses it.

### C8. `--polya-min-umis` default vs the paper's ≥2-molecule output (m14)
Users run the tool, not the paper. Options: change the shipped default to 2 (a tool release decision for Amir — interacts with PR #100's no-op guarantee), or add one Methods sentence saying why the caller stays at 1 (tiered output philosophy; the default *output* is a downstream filter). Recommend the one-sentence justification now, tool change on its own pre-registration.

### C9. Key points / abstract overlap (m13)
The matched-N sentence appears in both. GB allows it; trimming buys words for A1's new panel text. PI's stylistic call — no edit made.

### C10. Title
Resolved during the fix pass window: the PI set the title and author list directly (see MAIN.md comments, 2026-09-02). Reviewer m12 is thereby addressed. A.A.T. ORCID placeholder remains.

## D. Resolved in the fix pass (for the record — no action needed)
Audit PI-item 1 (Fig S12 over-pointing): the Results paragraph heading no longer cites Fig S12; the figure is cited only for the byte-for-byte version identity it actually carries. Review structure item 3 (two-truth composite debuting in Discussion): moved to Results with 25 §8.6's wording; Discussion keeps a back-reference. M2/M5c/M6b/M9i/M9iv/M9v/M11/M12 and minors 1, 2, 4, 5, 6, 10, 11 were addressed in-draft; see STATUS.md.
