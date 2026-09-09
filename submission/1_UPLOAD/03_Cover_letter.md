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
