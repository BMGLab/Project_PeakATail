# Declarations

## Availability of data and materials

**Tool code.** PeakATail is developed openly at https://github.com/BMGLab/PeakATail. The development record behind this paper is public: pull requests #92 (clip-seeded calling), #93 (tier-1 quantification and molecule counting), #96 (minus-strand internal-priming window fix), #97 (memory/CPU) and #100 (the `peakAtail-prime` defaults change, exploratory with respect to this paper's gates), and issues #94 (switch-test marker pre-selection), #95 (internal-priming filter and memory), #98 (`per_isoform` degenerate pairs), #99 (PAS-to-gene assignment in overlapping loci and the clip-rate sampling estimator) and #101 (dynamic-threshold crash and documentation defects). The frozen commits used in this paper are `4efeb125` (v1) and `9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53` (v2).

**Analysis code and manuscript record.** The benchmark scorer (`score_tool.py`), the replication filter (`replication_filter.py` v0.2.0), the long-read truth-set build scripts, the label-confirmation and universe policy files, every figure-generating script with its audit TSVs and script-written captions, the pre-registration documents with their appended amendments, and the run manifests are in the companion repository https://github.com/BMGLab/Project_PeakATail. An archived snapshot will be deposited: [[PLACEHOLDER: Zenodo DOI pending — repository must be public before a DOI can be minted]].

**Data.** All sequencing data are previously published or publicly distributed: the 10x Genomics public libraries `pbmc_10k_v3` and `pbmc4k` (BAMs from the 10x public data portal), the mouse testis libraries of GEO accession GSE104556, and the 17 lung-adenocarcinoma libraries of GEO accession GSE123904. The accuracy reference is the PolyASite 2.0 atlas (GRCh38 and GRCm38 representative sites). The PacBio Kinnex long-read truth and decoy point BEDs (x3p and GEM-X) are distributed with the companion repository together with the script that builds them.

## Ethics approval and consent to participate

Not applicable. This study analyses only previously published, publicly available datasets; no new human or animal data were collected.

## Consent for publication

Not applicable.

## Competing interests

[[PLACEHOLDER: competing interests statement — PI decision]]

## Funding

[[PLACEHOLDER: funding — PI decision]]

## Authors' contributions

[[PLACEHOLDER: author contributions (CRediT) — PI decision; authors set 2026-09-02: A.A.T., Y.K. (corresponding)]]

## Acknowledgements

[[PLACEHOLDER: acknowledgements — PI decision]]
