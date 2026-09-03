# Availability of data and materials — drop-in draft

Prepared 2026-09-02 for the Genome Biology Method-track submission.

**This file does not modify the manuscript.** `manuscript/00_draft/MAIN.md` is
off-limits to this task (legend sync pending), so Part A below is written to be
pasted over the existing `## Availability of data and materials` section by
whoever holds the draft. Parts B–E explain and justify every difference from
what is in the draft today, so nothing has to be taken on trust.

Every accession, identifier and URL in Part A is copied from a verified project
record; Part E names the record for each one. **Nothing is invented.** Where the
project record has no accession, Part A says so in the prose rather than
supplying a plausible-looking one.

---

## Part A — the statement, ready to paste

> ### Availability of data and materials
>
> **Tool code.** PeakATail (command-line name `ema`) is developed openly at
> https://github.com/BMGLab/PeakATail and is released under the
> [[PLACEHOLDER: licence — PI decision; see submission/LICENSE_DECISION.md]]
> licence. The development record behind this paper is public: pull requests #92
> (clip-seeded calling), #93 (tier-1 quantification and molecule counting), #96
> (minus-strand internal-priming window fix), #97 (memory and CPU) and #100 (the
> `peakAtail-prime` defaults change, exploratory with respect to this paper's
> gates), and issues #94 (switch-test marker pre-selection), #95
> (internal-priming filter and memory), #98 (`per_isoform` degenerate pairs),
> #99 (PAS-to-gene assignment in overlapping loci and the clip-rate sampling
> estimator) and #101 (dynamic-threshold crash and documentation defects). The
> frozen commits used in this paper are
> `4efeb1252e6d7c7d51b252b443b5be5947ae000f` (v1) and
> `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` (v2, the record for all reported
> numbers); both are reachable on the repository's default branch, `develop`. An
> archived snapshot of the v2 code is deposited at
> [[PLACEHOLDER: Zenodo DOI for the tool release]].
>
> **Analysis code and manuscript record.** The benchmark scorer
> (`score_tool.py`), the replication filter (`replication_filter.py` v0.2.0),
> the long-read truth-set build scripts, the label-confirmation and universe
> policy files, every figure-generating script with its audit TSVs and
> script-written captions, the pre-registration documents with their appended
> amendments, the per-tool environment exports, and the run manifests are in the
> companion repository https://github.com/BMGLab/Project_PeakATail, released
> under the same licence. An archived snapshot is deposited at
> [[PLACEHOLDER: Zenodo DOI for the companion repository]]. The figure source
> tables in `source_data/` are the plotted values themselves, so every number in
> every figure can be checked without access to the raw data or the pipeline.
>
> **Data.** No new sequencing data were generated; all sequencing data analysed
> here are previously published or publicly distributed.
>
> - **Human PBMC, donor 1** — the 10x Genomics public library `pbmc_10k_v3`
>   ("10k PBMCs from a Healthy Donor (v3 chemistry)"; Single Cell 3′ v3;
>   CellRanger 3.0.0 on the GRCh38-3.0.0 reference bundle), distributed at
>   https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
> - **Human PBMC, donor 2** — the 10x Genomics public library `pbmc4k` ("4k PBMCs
>   from a Healthy Donor"; Single Cell 3′ v2; CellRanger 2.1.0 on GRCh38-1.2.0),
>   distributed at
>   https://www.10xgenomics.com/datasets/4-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0
> - **Mouse testis** — GEO accession
>   [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556);
>   the two adult libraries analysed here are SRA runs SRR6129050 and SRR6129051,
>   realigned with STARsolo.
> - **Human lung adenocarcinoma cohort** — GEO accession
>   [GSE123904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123904)
>   (Laughney et al.); 17 libraries from 14 patients entered the analysis and 15
>   libraries from 12 patients form the pre-registered replication primary
>   (GSM3516664 is excluded by a pre-declared data-quality rule and GSM3516671
>   yields no testable cell-type pair).
> - **Long-read truth** — publicly released PacBio Kinnex single-cell RNA PBMC
>   datasets, downloaded on 2026-08-19 from
>   https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/ : the
>   sets `DATA-Revio-Kinnex-PBMC-10x3p` (the 10x 3′ v3.1 "x3p" preparation) and
>   `DATA-Revio-Kinnex-PBMC-10kcells-10xGEMX3p` (the GEM-X v4 preparation).
>   PacBio distributes these from its own download portal and assigns them no
>   INSDC accession. The poly(A)-verified truth and internal-priming decoy point
>   BEDs derived from them, and the script that builds them, are distributed with
>   the companion repository.
> - **Reference resources** — the PolyASite 2.0 representative-site atlas
>   (https://polyasite.unibas.ch), GRCh38 and GRCm38 builds; Ensembl 99 (GRCh38)
>   and Ensembl 102 (GRCm38) annotation.

---

## Part B — what changed from the draft, and why

Six changes. The draft's existing prose is otherwise preserved verbatim, because
it has already been through the number-trace audit.

1. **The licence is now named** (as a placeholder). *This is the one addition
   that is not optional.* Genome Biology's source-code policy requires code to be
   "released under a license complying with an Open Source Definition, as defined
   by the Open Source Initiative", and requires that to be stated in the
   manuscript. The current draft names no licence anywhere. See
   `submission/LICENSE_DECISION.md` — and note that the tool repository currently
   publishes **two contradictory licences** (MIT on `develop`, GPL-3.0 on `main`),
   which must be resolved before this sentence can be written truthfully.

2. **v1's commit is given in full** — `4efeb1252e6d7c7d51b252b443b5be5947ae000f`
   rather than the abbreviated `4efeb125`. The draft already gives v2 in full;
   this makes the pair symmetric. Abbreviations are not stable identifiers as a
   repository grows.

3. **"both are reachable on the repository's default branch, `develop`" added.**
   A reader who clones and lands on `main` will not find either commit — `main`
   is 322 commits behind and contains neither. Verified:
   `git merge-base --is-ancestor 9dfdefb origin/develop` → YES;
   `… origin/main` → NO. One clause prevents a reproducibility complaint.

4. **A Zenodo DOI slot was added for the tool as well as the companion repo.**
   Genome Biology asks authors to "archive the version used in the manuscript in
   a DOI-assigning repository". The version used is the tool at `9dfdefb`, so a
   companion-only DOI does not fully satisfy the policy. Both slots are
   placeholders and both are blocked behind the same work
   (`COMPANION_REPO_PLAN.md`).

5. **The Data paragraph became a list with retrieval routes.** The draft says
   "BAMs from the 10x public data portal" without naming the portal, and gives
   GSE104556 without the SRA runs. Reviewers of a Method paper check whether they
   can obtain the exact inputs. The 10x URLs, the SRR run ids and the Kinnex
   dataset directory names all come from the project's own download records
   (Part E), not from a search.

6. **The Kinnex entry now names the actual source, and says plainly that there
   is no accession.** This closes reviewer point M8 and fills the draft's open
   `[[CITE: the public Kinnex PBMC datasets used as long-read truth (x3p and
   GEM-X library preparations) — data citation with accession/identifier]]`. The
   companion repo's `DATA_ACCESSIONS.md` records only the portal URL, and no
   INSDC accession for these sets was found in any project record. **A URL plus
   an access date plus an explicit "no accession is assigned" is a correct data
   citation; a fabricated accession is a retraction risk.** If PacBio has since
   deposited these under an accession, replace the sentence — but verify it,
   do not assume it.

**One thing deliberately left alone:** the draft's "17 lung-adenocarcinoma
libraries of GEO accession GSE123904" is expanded to the audited "17 libraries
from 14 patients … 15 libraries from 12 patients" phrasing used everywhere else
in the manuscript. STATUS.md records this cohort phrasing as a checked
forbidden-phrasing invariant, so the two forms must not drift apart.

---

## Part C — placeholders in Part A, and who clears them

| Placeholder | Blocked on | Who |
|---|---|---|
| licence name | `LICENSE_DECISION.md`; a licence must be chosen **and** the tool repo's MIT/GPL branch conflict resolved | PI |
| Zenodo DOI, tool | tagging and archiving the tool at `9dfdefb` (`COMPANION_REPO_PLAN.md` §5 step 28) | PI |
| Zenodo DOI, companion | the whole `COMPANION_REPO_PLAN.md` critical path — repo must be public first | PI |

The other declarations in the draft (competing interests, funding, CRediT
contributions, acknowledgements, A.A.T.'s ORCID) are outside this task's scope
and remain as they are in `PI_DECISIONS.md`.

---

## Part D — optional "Availability and requirements" block

Genome Biology does **not** require the structured list that some BMC journals
ask of Software articles; its policy requires an OSI-compliant licence, a
DOI-assigning archive, and both stated in the manuscript, all of which Part A
does. Include this only if the editor asks. Every value is taken from the tool
repo's `pyproject.toml` and `README.md` at `9dfdefb`.

> **Project name:** PeakATail
> **Project home page:** https://github.com/BMGLab/PeakATail (documentation: https://bmglab.github.io/PeakATail/)
> **Archived version:** [[PLACEHOLDER: Zenodo DOI for the tool release]]
> **Operating system(s):** Platform independent (developed and benchmarked on Linux)
> **Programming language:** Python ≥ 3.11
> **Other requirements:** pysam ≥ 0.22, anndata ≥ 0.10, scanpy ≥ 1.10, pandas ≥ 2.1, scipy ≥ 1.11, numpy ≥ 1.26, pybedtools ≥ 0.10, statsmodels ≥ 0.14 (full pin: `pyproject.toml` and `uv.lock` at the frozen commit)
> **License:** [[PLACEHOLDER: licence — PI decision]]
> **Any restrictions to use by non-academics:** None

---

## Part E — provenance of every identifier in Part A

Nothing below was retrieved by search; each is a project record on this machine
or a read-only query whose command is shown.

| Item in Part A | Value | Source of record |
|---|---|---|
| Tool repo URL | github.com/BMGLab/PeakATail | MAIN.md Declarations; confirmed PUBLIC via `gh api repos/BMGLab/PeakATail` |
| PR / issue numbers | #92, #93, #96, #97, #100 / #94, #95, #98, #99, #101 | MAIN.md Declarations, cross-checked against `manuscript/github/PUBLISHED.md` |
| v1 commit | `4efeb1252e6d7c7d51b252b443b5be5947ae000f` | `git rev-parse 4efeb125` in `tools/PeakATail` |
| v2 commit | `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` | MAIN.md Methods; `git rev-parse 9dfdefb` |
| "on `develop`" | both are ancestors of `origin/develop` | `git merge-base --is-ancestor <c> origin/develop` → YES for both |
| Companion repo URL | github.com/BMGLab/Project_PeakATail | MAIN.md Declarations |
| `pbmc_10k_v3` URL | www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0 | `DATA_ACCESSIONS.md` §1 records the *legacy* `support.10xgenomics.com/…/3.0.0/pbmc_10k_v3` path; 10x has since migrated its portal, and the current landing page was confirmed 2026-09-03 |
| `pbmc4k` URL | www.10xgenomics.com/datasets/4-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0 | same migration; identified as `pbmc4k` by the recorded download URL `https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_possorted_genome_bam.bam` in the project record, and corroborated by the dataset page's own CellRanger summary — **4,340 estimated cells, mean 87,433 reads per cell**, matching MAIN.md Methods exactly |
| GSE104556 | GEO accession | MAIN.md Methods |
| SRR6129050 / SRR6129051 | SRA runs, mouse 1 / mouse 2 | companion repo `DATA_ACCESSIONS.md` §1 (retrieval script `fetch_gse104556.sh`) |
| GSE123904 | GEO accession | MAIN.md Methods |
| 17 libraries / 14 patients; 15 / 12 | cohort and replication primary | MAIN.md Methods (audited; NUMBERS_LEDGER.tsv) |
| GSM3516664, GSM3516671 | the two named exclusions | MAIN.md Methods |
| Kinnex portal | downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/ | `DATA_ACCESSIONS.md` §1 **and** the download script `data/benchmark/kinnex_pbmc_10x3p/run.sh` line 7 |
| Kinnex set names | `DATA-Revio-Kinnex-PBMC-10x3p`, `DATA-Revio-Kinnex-PBMC-10kcells-10xGEMX3p` | `run.sh` lines 23–24 / `run2.sh` lines 8–9 |
| Kinnex download date | 2026-08-19 | `data/benchmark/kinnex_pbmc_10x3p/DONE.ok` ("all files verified … 19 Ağu 2026 05:59:08 +03") and file mtimes |
| "no INSDC accession" | absence of record | no accession for these sets appears in `DATA_ACCESSIONS.md`, in MAIN.md, or in any download script; MAIN.md still carries the open `[[CITE]]` asking for one |
| PolyASite 2.0 | polyasite.unibas.ch | project record; URL confirmed to return HTTP 200 |
| Ensembl 99 / 102 | annotation builds | MAIN.md Methods |
| Part D requirements | Python ≥ 3.11 and the dependency floors | `tools/PeakATail/pyproject.toml` at `9dfdefb` |

### Note on the two 10x landing-page URLs

The `support.10xgenomics.com/single-cell-gene-expression/datasets/…` URL recorded
in `DATA_ACCESSIONS.md` **no longer serves the dataset pages** — 10x has migrated
its portal to `www.10xgenomics.com/datasets/…`. Part A therefore cites the
current pages, verified 2026-09-03; `DATA_ACCESSIONS.md` should be updated to
match (`COMPANION_REPO_PLAN.md` §3 lists it with the other manifest fixes).

Two things worth carrying forward:

- The stable, machine-readable fallbacks are the `cf.10xgenomics.com` file URLs
  already in the project record, e.g.
  `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam`.
  These have survived the portal migration and are what the download scripts
  used. Consider citing them alongside the landing pages.
- 10x publishes these datasets under **CC BY 4.0**. Nothing in this paper's use
  is inconsistent with that, but the licence is worth a line if a reviewer or
  the editor asks about redistribution — note that the companion repository
  redistributes only *derived* point BEDs and figure tables, never the 10x BAMs.

Both 10x URLs return HTTP 403/429 to a plain `curl` (bot protection), so they
cannot be verified from the analysis machine's shell. **Open both in a browser
before submission.**
