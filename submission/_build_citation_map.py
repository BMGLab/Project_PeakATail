#!/usr/bin/env python3
import re, os, sys

WD = "/mnt/ssd1/Projects/PeakATail_wd"
FILES = ["manuscript/00_draft/MAIN.md", "manuscript/00_draft/SUPPLEMENTARY.md"]
OUT = os.path.join(WD, "submission/CITATION_MAP.tsv")

# token text (inside [[CITE: ...]]) -> (resolved, verification_url, note)
M = {
 "review of APA mechanisms and functional consequences": (
   "tian2017apa;mitschka2022apa",
   "https://doi.org/10.1038/nrm.2016.116 ; https://doi.org/10.1038/s41580-022-00507-5",
   "Two verified reviews offered. Tian & Manley 2017 is the general APA-mechanism review; Mitschka & Mayr 2022 is the one that covers stability/localization/translational output, which is what the sentence claims. Cite both or pick one - PI call, no gap."),
 "3'UTR shortening in proliferating and transformed cells": (
   "sandberg2008proliferating;mayr2009shortening",
   "https://doi.org/10.1126/science.1155390 ; https://doi.org/10.1016/j.cell.2009.06.016",
   "The sentence makes two claims: proliferating (Sandberg 2008) and transformed (Mayr & Bartel 2009). Both needed."),
 "APA remodelling during spermatogenesis": (
   "li2016spermatogenesis", "https://doi.org/10.1186/s12915-016-0229-6",
   "Genome-wide 3'UTR shortening across mouse spermatogenesis; direct match for 'remodel 3'-end usage wholesale'."),
 "internal-priming artefacts in oligo(dT)-primed sequencing": (
   "nam2002internalpriming", "https://doi.org/10.1073/pnas.092140899",
   "The founding demonstration of oligo(dT) internal priming at genomic A-rich stretches. Consider adding gruber2016polyasignals for the 3'-end-seq-era filtering practice."),
 "Sierra paper": ("patrick2020sierra", "https://doi.org/10.1186/s13059-020-02071-7", ""),
 "scAPAtrap paper": ("wu2021scapatrap", "https://doi.org/10.1093/bib/bbaa273", ""),
 "polyApipe paper": ("harrison2019polyapipe", "https://doi.org/10.7490/f1000research.1117076.1",
   "CAUTION - gap G2. polyApipe has no peer-reviewed algorithm paper. This is an F1000Research poster DOI (verified via Crossref) that the tool's own README points to; the poster title does not name polyApipe. PI must approve citing a poster for a benchmarked competitor."),
 "SCAPTURE paper": ("li2021scapture", "https://doi.org/10.1186/s13059-021-02437-5", ""),
 "scTail paper": ("hou2025sctail", "https://doi.org/10.1186/s13059-025-03710-7", ""),
 "scUTRquant paper": ("fansler2024scutrquant", "https://doi.org/10.1038/s41467-024-48254-9",
   "Tool name is not in the title; the verified abstract names scUTRquant as the pipeline this paper introduces."),
 "PolyASite 2.0 paper": ("herrmann2020polyasite", "https://doi.org/10.1093/nar/gkz918", ""),
 "recent benchmarks of single-cell APA tools": (
   "tian2025benchmark;zhao2024apaguidelines",
   "https://doi.org/10.1093/nargab/lqaf056 ; https://doi.org/10.1101/2024.11.29.626111",
   "Only one of the two is peer-reviewed (Tian 2025, NAR Genom Bioinform). Zhao & Rattray 2024 is still a bioRxiv preprint as of 2026-09-02 - re-check before submission. The sentence says 'published benchmarks' (plural), so at least the preprint status must be visible in the reference list."),
 "Kinnex paper": ("alkhafaji2024kinnex", "https://doi.org/10.1038/s41587-023-01815-7",
   "MAS-ISO-seq method paper; Kinnex is its commercialisation. There is no separate 'Kinnex' publication."),
 "CellRanger paper": ("zheng2017cellranger", "https://doi.org/10.1038/ncomms14049",
   "Standard citation for the 10x Chromium platform and Cell Ranger; 10x ships no separate Cell Ranger paper."),
 "STARsolo paper": ("kaminow2021starsolo;dobin2013star",
   "https://doi.org/10.1101/2021.05.05.442755 ; https://doi.org/10.1093/bioinformatics/bts635",
   "STARsolo is still a bioRxiv preprint (searched 2026-09-02, no journal version found). Pair it with the STAR aligner paper, as is standard practice."),
 "Benjamini-Hochberg FDR procedure": ("benjamini1995fdr", "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
   "MERGE PAIR 1 of 3: same target as the line-199 token 'Benjamini and Hochberg FDR'."),
 "Benjamini and Hochberg FDR": ("benjamini1995fdr", "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
   "MERGE PAIR 1 of 3: same target as the line-110 token 'Benjamini-Hochberg FDR procedure'."),
 "mouse spermatogenesis 10x scRNA-seq dataset GSE104556": (
   "lukassen2018testis;gse104556",
   "https://doi.org/10.1038/sdata.2018.192 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556",
   "MERGE PAIR 2 of 3 (with the line-165 token). Accession + source publication, as the brief requires. The GEO<->paper link was confirmed twice: the paper's PMC full text names GSE104556, and NCBI esummary for the series returns PMID 30204153."),
 "study describing the GSE104556 mouse testis 10x dataset": (
   "lukassen2018testis;gse104556",
   "https://doi.org/10.1038/sdata.2018.192 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556",
   "MERGE PAIR 2 of 3 (with the line-114 token). Note there is a second Lukassen et al. 2018 paper (Sci Rep 8:6521) analysing the same cells; the Scientific Data descriptor is the one that deposits GSE104556 and is the correct dataset citation."),
 "literature sources of the spermatogenesis 3'UTR-shortening gene panel": (
   "UNRESOLVED", "",
   "GAP G1 - the blocking one. The 27-gene panel is hard-coded at scripts/manuscript_figures/spermatogenesis_control_steps/step7_genes.py:70-71 with no source recorded anywhere in the repository. The PI must supply the paper(s) it was drawn from, or the word 'literature' must come out of the sentence. See CITATION_GAPS.md."),
 "Laughney et al. lung adenocarcinoma scRNA-seq cohort": (
   "laughney2020lungmet;gse123904",
   "https://doi.org/10.1038/s41591-019-0750-6 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123904",
   "MERGE PAIR 3 of 3 (with the line-165 token). GEO<->paper link confirmed: NCBI esummary for GSE123904 returns PMID 32042191."),
 "Laughney et al. lung adenocarcinoma single-cell atlas": (
   "laughney2020lungmet;gse123904",
   "https://doi.org/10.1038/s41591-019-0750-6 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123904",
   "MERGE PAIR 3 of 3 (with the line-126 token). This token sits next to the accession in Methods, so the @misc data citation belongs here at minimum."),
 "Signac TF-IDF paper": ("stuart2021signac", "https://doi.org/10.1038/s41592-021-01282-5",
   "An Author Correction exists (Nat Methods 2022;19:257) but does not replace the article."),
 "bedtools paper": ("quinlan2010bedtools", "https://doi.org/10.1093/bioinformatics/btq033", ""),
 "CellTypist paper": ("dominguezconde2022celltypist", "https://doi.org/10.1126/science.abl5197",
   "Tool name is not in the title; the verified abstract names CellTypist as the tool this paper introduces."),
 "the public Kinnex PBMC datasets used as long-read truth (x3p and GEM-X library preparations) - data citation with accession/identifier": (
   "PROVISIONAL:pacbio_kinnex_pbmc_x3p_gemx", "",
   "GAP G3. NOT RESOLVED - deliberately kept out of references.bib. The datasets are PacBio-hosted vendor downloads with no accession and no DOI: DATA-Revio-Kinnex-PBMC-10x3p and DATA-Revio-Kinnex-PBMC-10kcells-10xGEMX3p under https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/. Direct fetch of that host was REFUSED from this environment on 2026-09-02 (ECONNREFUSED 38.99.102.168:443 via WebFetch; curl timed out), so I could not see the listing myself. A draft entry and the evidence that does exist are in CITATION_GAPS.md."),
}

def normkey(s):
    # normalise the unicode punctuation the draft uses so the map dict can be plain ASCII
    return (s.replace("–", "-").replace("—", "-").replace("’", "'")
             .replace("′", "'").replace("‘", "'"))

rows = []
missing = []
pat = re.compile(r"\[\[CITE:\s*(.*?)\s*\]\]")
for rel in FILES:
    path = os.path.join(WD, rel)
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            for m in pat.finditer(line):
                tok = m.group(1)
                key = normkey(tok)
                if key not in M:
                    missing.append((rel, lineno, tok, key))
                    continue
                resolved, url, note = M[key]
                rows.append((rel, str(lineno), "[[CITE: %s]]" % tok, resolved, url, note))

if missing:
    for r in missing:
        print("UNMAPPED TOKEN:", r, file=sys.stderr)
    sys.exit(1)

header = [
 "# CITATION_MAP.tsv - one row per [[CITE]] token in the PeakATail draft. Built 2026-09-02.",
 "# Columns: file  line  token  resolved  verification_url  note",
 "# resolved vocabulary:",
 "#   <bibkey>            - verified entry present in submission/references.bib. Several keys separated",
 "#                         by ';' mean the token should carry more than one reference (e.g. accession +",
 "#                         source publication, or aligner + its solo mode).",
 "#   UNRESOLVED          - no verified target exists. Do NOT invent one; see CITATION_GAPS.md.",
 "#   PROVISIONAL:<key>   - a target that is identified but NOT verified to the standard used here, and",
 "#                         therefore deliberately absent from references.bib. Treat as unresolved until",
 "#                         the PI confirms it. See CITATION_GAPS.md.",
 "# verification_url is the canonical identifier of the record. The exact API URL that was queried, and",
 "# the access date, are recorded per entry in submission/references.bib.",
 "# Line numbers are against manuscript/00_draft/MAIN.md as of 2026-09-02 (327 lines, 90110 bytes).",
 "# SUPPLEMENTARY.md carries no [[CITE]] tokens.",
 "file\tline\ttoken\tresolved\tverification_url\tnote",
]
with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(header) + "\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

print("rows:", len(rows))
distinct = sorted({r[2] for r in rows})
print("distinct tokens:", len(distinct))
unres = [r for r in rows if r[3].startswith("UNRESOLVED") or r[3].startswith("PROVISIONAL")]
print("unresolved/provisional token instances:", len(unres))
print("distinct unresolved/provisional targets:", len({r[2] for r in unres}))
