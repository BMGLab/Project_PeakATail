<!--
cover_letter.md — draft cover letter to the Genome Biology editors, written in the voice of the
corresponding author (Yasin Kaymaz). Drafted 2026-09-02 against manuscript/00_draft/MAIN.md
(abstract, Background, Discussion, Conclusions) and manuscript/00_draft/REVIEW_SIMULATION.md, so
the letter meets the simulated referee's obvious objections in advance rather than after review.

Every quantity below is quoted from MAIN.md's verified abstract/key points. Nothing here is a new
number, and no claim is stronger than the licensed wording in the draft ("the most atlas-concordant
de novo call set in the benchmark", never bare "precision" or "outperforms").

Genome Biology's own cover-letter requirements (verified 2026-09-02 on the live submission
guidelines) are: why the manuscript should be published in Genome Biology; any issues relating to
journal policies; a declaration of potential competing interests; confirmation that all authors have
approved the manuscript for submission; confirmation that it has not been published or submitted
elsewhere. All five are covered below.

LENGTH: 744 words including the letterhead and sign-off (~600 of prose), which sets to about one
page at 11 pt with a compact letterhead. If it must be shorter, the first thing to cut is the closing
sentence of the negative-results paragraph ("In the same spirit … donor-mismatched") — it pre-empts
referee concerns the manuscript already states, so nothing load-bearing is lost.

NO DATE LINE: none was invented. Add the submission date above "To the Editors" when sending.

ADVERSARIAL VERIFICATION, 2026-09-03. Every factual sentence was re-checked against MAIN.md and
against live records. All quantities reconcile: 0.706 / 0.745 / 0.757 / 0.828, 76.5%, 15,942 of
2,128,711, ten shuffles at p <= 0.091, 3.0% vs 13.0-24.7%, 48.5% vs the 70% target, "16 points of
atlas-agreement precision ... losing 8-11%", and commits 4efeb125 / 9dfdefb. The ethics sentence
matches MAIN.md's "Ethics approval" section. All ten GitHub references were resolved against the
public API: PRs #92, #93, #96, #97 are merged and #100 is open; issues #94, #95, #98, #99, #101 are
open. (Note: manuscript/github/PUBLISHED.md does NOT list PR #92 — the PR is real and merged, so
the letter is right and PUBLISHED.md is simply incomplete.) Forbidden-phrasing scan is clean: no
"outperforms", no priority claim; every "first" is either the coined descriptor "precision-first"
or ordinary enumeration, and every "precision" is qualified as "atlas-agreement precision".

TWO WORDING FIXES were applied by that pass:
  * "In a 12-patient lung-adenocarcinoma cohort" -> "In 15 libraries from 12 lung-adenocarcinoma
    patients", to match MAIN.md's own phrasing exactly.
  * "a 16-gene literature panel" -> "a curated literature panel". THIS IS A MANUSCRIPT DEFECT, not
    just a letter one. The curated panel is 27 genes (hard-coded in
    scripts/manuscript_figures/spermatogenesis_control_steps/step7_genes.py); 16 is the number
    MEASURABLE in both mice under the depth guard. results/stage3_spermatogenesis_v2/
    REPORT_PROVISIONAL.md states it correctly ("The 27-gene curated panel ... has 16 genes
    measurable in both mice"), and so do MAIN.md's Results body and the Fig 5 caption ("of the 16
    genes ... measurable in both mice"). But MAIN.md's Key points and Discussion compress this to
    "a 16-gene literature panel", which misstates the panel's size. The letter now avoids the
    count entirely so it cannot be wrong; MAIN.md still needs the fix. See CITATION_GAPS.md item G8.

THREE PLACEHOLDERS REMAIN, all PI calls: suggested reviewers, opposed reviewers, and the competing-
interests wording (which must match the manuscript's Competing interests section verbatim). The
Zenodo DOI slot inherits the manuscript's placeholder and is blocked on the companion repository
being made public; do NOT send this letter with a repository-availability sentence that is not yet
true. See SUBMISSION_CHECKLIST.md in this directory for the full pre-send list.
-->

**Yasin Kaymaz**
Department of Bioengineering, Faculty of Engineering
Ege University, İzmir, Türkiye
yasin.kaymaz@ege.edu.tr · ORCID 0000-0002-9725-7536

To the Editors
*Genome Biology*

Dear Editors,

We submit **"PeakATail: precision-first poly(A)-site calling and calibrated alternative polyadenylation analysis in single-cell RNA-seq"** for consideration as a Method article.

A substantial tool ecosystem now calls poly(A) sites (PAS) from droplet 3′ single-cell data and ranks cell-type differences in their usage. What it lacks is a way of knowing how much of a call set — or of the switch list computed from it — deserves belief. PeakATail seeds PAS from direct poly(A) evidence, the non-templated soft clips that record cleavage on an individual molecule, tiers sites by molecule-counted support, filters internal priming, and couples them to a switch test whose false-discovery-rate behaviour we measured under label permutation rather than assumed. Its pre-registered default is the most atlas-concordant de novo call set in our benchmark on both datasets (atlas-agreement precision 0.706 PBMC, 0.745/0.757 in two mice, 0.828 un-tuned on a second donor), and 76.5% of its PBMC calls match a donor-mismatched long-read 3′ end within 25 bp. In 15 libraries from 12 lung-adenocarcinoma patients, 15,942 of 2,128,711 switch hypotheses replicate across patients, with none replicated in ten label-shuffle nulls.

This belongs in the Method track for two reasons. The first is what a reader can do with the tool: six tools run side by side on the same BAMs, on one machine, through one scorer, with shuffled-coordinate nulls under every headline metric and matched-call-count curves rather than a single operating point — so the ordering reflects the callers, not their default call budgets. The second is what a reader can take away without the tool. The evaluation is reusable: gates registered before the final run and auditable against git timestamps; a label-permutation harness that measures whether a differential test controls its own error rate; replication across patients, not libraries, as the condition for reporting a switch. Applied to our own defaults, that harness found only one of six switch-test configurations controls the FDR (3.0% of label-shuffled tests p < 0.05) against 13.0–24.7% for the configurations we had been shipping. We would rather the field had the harness than that we had a better-looking table.

The paper reports negative results, deliberately and in the main text. A pre-registered "trusted novel PAS" definition reached 48.5% long-read concordance against its 70% target and failed; no site here is called trusted-novel. A clustering claim was withdrawn when its own ablation showed per-gene totals do as well; a curated literature panel of spermatogenesis 3′UTR shortening does not reproduce under our PAS-to-gene assignment. An atlas-trained scoring model gained 16 points of atlas agreement while losing 8–11% against long reads — a caution about atlas-only benchmarking that applies to our own headline metric as much as to anyone else's. These are findings: a pre-registration that can only be confirmed is decoration, and each failure marks where the measurements stop being trustworthy. In the same spirit we are explicit about where our evidence is thinner — recall at the precision-first operating point is below polyApipe's, the cohort null resolves only to its ten-permutation floor (empirical p ≤ 0.091, reported as such and never as an FDR), and the long-read truth is donor-mismatched.

The development record is public: the pull requests and issues behind this work, including tool defects that were open when these results were produced and that limit what we claim, are in the tool repository github.com/BMGLab/PeakATail (PRs #92, #93, #96, #97, #100; issues #94, #95, #98, #99, #101), and the manuscript pins the two frozen commits behind its numbers, `4efeb125` (v1) and `9dfdefb` (v2, the record for every reported value). All analysis code, every figure script with the audit table of its plotted values, and the pre-registration documents with their amendments are in the companion repository github.com/BMGLab/Project_PeakATail, archived at [[PLACEHOLDER: Zenodo DOI — companion repository must be public before a DOI can be minted]].

The work is original, has not been published, and is not under consideration elsewhere. All authors have read and approved the manuscript for submission. All data analysed are previously published or publicly distributed, so no ethics approval was required. Competing interests: [[PLACEHOLDER: competing interests declaration — PI; must match the manuscript's Competing interests section verbatim]].

Suggested reviewers: [[PLACEHOLDER: suggested reviewers — PI]].
Reviewers we ask you to exclude: [[PLACEHOLDER: opposed reviewers — PI]].

Thank you for considering the manuscript.

Yours sincerely,

**Yasin Kaymaz**, on behalf of Amir Amiri Tabat and myself
Corresponding author
Department of Bioengineering, Faculty of Engineering, Ege University, İzmir, Türkiye
