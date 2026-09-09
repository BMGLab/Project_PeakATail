# Open items — what still blocks submission
Ordered by who owns them. Nothing here is something I can decide.

## 1. Declarations — PI only
Six placeholders sit in the manuscript and the paper cannot be submitted with any of
them unresolved. They are marked in the draft as `[[PLACEHOLDER ...]]`:
- `[[PLACEHOLDER: A.A.T. ORCID if available]]`
- `[[PLACEHOLDER: Zenodo DOI pending — repository must be public before a DOI can be minted]]`
- `[[PLACEHOLDER: acknowledgements — PI decision]]`
- `[[PLACEHOLDER: author contributions (CRediT) — PI decision. Authors are now set: A.A.T. and Y.K. (corresponding). Sugges`
- `[[PLACEHOLDER: competing interests statement — PI decision]]`
- `[[PLACEHOLDER: funding — PI decision]]`

## 2. Citation gaps — PI only
Two citations cannot be filled without your call:
- `[[CITE: literature sources of the spermatogenesis 3′UTR-shortening gene panel]]`
- `[[CITE: the public Kinnex PBMC datasets used as long-read truth (x3p and GEM-X library preparations) — data citation with accession/identifier]]`

## 3. Licence — PI only
Three candidate licences are staged in `3_SUPPORTING/` and none is chosen. There is also
a conflict to resolve first: the tool repo carries **MIT on `develop`** and **GPL-3.0 on
`main`**. That must be settled before anything is archived, because the archived licence
is what the paper cites. See `3_SUPPORTING/LICENSE_DECISION.md`.

## 4. Release and DOI — blocked on Amir, then you
The manuscript names the command `peakatail`, which no released version provides yet.
The chain is: **PR #111 merged -> a release cut -> that tag archived to Zenodo -> the DOI
goes in the paper**. The repo is now public, so the Zenodo precondition is met.
If #111 will not merge in time, the paper must be reverted to `ema` rather than naming a
command that does not exist. See `manuscript/00_draft/PI_DECISIONS.md` item C11.

## 5. Supplementary figures S5 and S6 — a decision at freeze
Both are declared in the supplementary as **[In preparation]** and neither file exists.
The draft already states the plan: if they are not built at freeze, the slots are dropped
and S7-S12 renumber in a single pass. Ten supplementary figures exist and ship today.
Decide: build them, or drop and renumber.

## 6. Cost — unresolved, and it may change the target journal
Genome Biology is fully open access: there is no subscription route, so an APC is
unavoidable there without a discretionary waiver. Verified on the publishers' own pages:
Genome Research charges **$1,500 on every accepted paper** regardless of open access;
Bioinformatics went fully OA in 2023 and Briefings in Bioinformatics in 2024. The
TUBITAK-Springer Nature national agreement does list **Ege University**, but it covers
**hybrid journals only** and has **exhausted its annual quota**. Springer Nature waivers
require a low-income economy; Turkiye is upper-middle-income and qualifies for neither
the waiver nor the 50% discount.

Practical options: request a discretionary waiver at submission (BMC has one for authors
in financial need, and it must be asked for at submission, not after); or target **PLOS
Computational Biology**, whose fee assistance is means-tested on circumstances rather than
nationality. Ask the Ege library what current agreements cover — they will know faster
than the web does.

---

Everything else is done: text, figures, typography, references, package structure.
