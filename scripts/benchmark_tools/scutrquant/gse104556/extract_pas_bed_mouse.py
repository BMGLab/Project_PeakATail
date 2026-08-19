#!/usr/bin/env python
"""Derive scUTRquant PAS calls (pas.bed) for GSE104556 mouse testis, per mouse.

Mouse port of scripts/benchmark_tools/scutrquant/extract_pas_bed.py (pbmc_10k_v3);
identical derivation, run once per replicate:
 1. txs.mtx (bustools count --em, target utrome_mm10_v2) barcodes x merged
    isoforms; columns = txs.genes.txt = representative transcript_id of each
    <200 nt merge group.
 2. Cell filter = total UMIs across isoforms >= MIN_UMIS (500 = the run's
    config min_umis).
 3. Detected isoform = nonzero summed count in >= 1 passing cell.
 4. PAS = genomic 3' end of the representative transcript in the mm10 UTRome
    GTF ('+': GTF end, '-': GTF start) as 1-bp BED6.
 5. Chromosomes UCSC -> Ensembl (strip 'chr', chrM -> MT); mm10 == GRCm38.

scUTRquant is ANNOTATION-BASED (fixed mm10 UTRome catalog, Fansler et al. 2024);
it cannot discover novel PAS -- precision vs an annotation atlas is
near-tautological; headers carry the caveat.

Usage: python extract_pas_bed_mouse.py  (paths hard-coded below)
"""
import sys, re, datetime
import numpy as np
import pandas as pd

SQ   = "/mnt/ssd1/Projects/PeakATail_wd/tools/scUTRquant"
GTF  = f"{SQ}/extdata/targets/utrome_mm10_v2/utrome.e30.t5.gc25.pas3.f0.9999.w500.gtf"
BASE = "/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/gse104556/scutrquant"
MIN_UMIS = 500

# ---- GTF: representative transcript -> 3' end (parsed once) ----------------
pat = {k: re.compile(k + r' "([^"]+)"') for k in ("transcript_id", "gene_id", "gene_name")}
gtf = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[2] != "transcript":
        continue
    a = {k: (m.group(1) if (m := p.search(f[8])) else "NA") for k, p in pat.items()}
    gtf[a["transcript_id"]] = (f[0], int(f[3]), int(f[4]), f[6], a["gene_id"], a["gene_name"])
print(f"GTF: {len(gtf)} transcripts", file=sys.stderr)

def ens(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if c == "M" else c

def ckey(c):
    return (0, int(c)) if c.isdigit() else (1, {"X": "1", "Y": "2", "MT": "3"}.get(c, c))

for mouse in ("mouse1", "mouse2"):
    KDIR = f"{BASE}/{mouse}"
    with open(f"{KDIR}/txs.mtx") as fh:
        hdr = 0
        for line in fh:
            hdr += 1
            if line.startswith("%"):
                continue
            n_bx, n_tx, n_nz = (int(x) for x in line.split())
            break
    print(f"[{mouse}] mtx: {n_bx} barcodes x {n_tx} isoforms, {n_nz} nonzeros",
          file=sys.stderr)

    mtx = pd.read_csv(f"{KDIR}/txs.mtx", sep=" ", skiprows=hdr, header=None,
                      names=["bx", "tx", "ct"],
                      dtype={"bx": np.int32, "tx": np.int32, "ct": np.float64})
    assert len(mtx) == n_nz

    row_tot = np.bincount(mtx["bx"], weights=mtx["ct"], minlength=n_bx + 1)
    cells   = row_tot >= MIN_UMIS            # index 0 unused (mtx is 1-based)
    n_cells = int(cells[1:].sum())

    keep = cells[mtx["bx"].to_numpy()]
    tx_k, ct_k = mtx["tx"].to_numpy()[keep], mtx["ct"].to_numpy()[keep]
    umi_filt   = np.bincount(tx_k, weights=ct_k, minlength=n_tx + 1)
    ncell_filt = np.bincount(tx_k, minlength=n_tx + 1)
    umi_all    = np.bincount(mtx["tx"], weights=mtx["ct"], minlength=n_tx + 1)
    nbx_all    = np.bincount(mtx["tx"], minlength=n_tx + 1)

    tx_ids = [l.strip() for l in open(f"{KDIR}/txs.genes.txt")]
    assert len(tx_ids) == n_tx

    rows = []
    for i, tid in enumerate(tx_ids, start=1):
        if ncell_filt[i] == 0:
            continue
        chrom, start, end, strand, gid, gname = gtf[tid]
        pas = end if strand == "+" else start      # 1-based PAS position
        rows.append((ens(chrom), pas - 1, pas, tid, strand, gid, gname,
                     int(ncell_filt[i]), umi_filt[i], int(nbx_all[i]), umi_all[i]))

    n_det_filt = len(rows)
    n_det_all  = int((nbx_all[1:] > 0).sum())
    rows.sort(key=lambda r: (ckey(r[0]), r[1]))

    stamp = datetime.date.today().isoformat()
    hdr_lines = [
        f"# scUTRquant v0.5.1 PAS calls, GSE104556 {mouse} (mouse testis, 10xv2, "
        f"STARsolo BAM), generated {stamp}",
        "# ANNOTATION-BASED tool: sites are the fixed mm10 UTRome catalog (GENCODE 3'",
        "#  ends + cleavage sites, isoforms <200 nt apart merged; Fansler et al. 2024,",
        "#  Nat Commun 15:4050; target utrome_mm10_v2) filtered by detection in this",
        "#  dataset -- scUTRquant cannot discover novel PAS; precision vs an",
        "#  annotation atlas is near-tautological.",
        "# Derivation (identical to pbmc_10k_v3 pas.bed): txs.mtx (bustools count",
        "#  --em) columns = merge-group representative transcript_ids; cells =",
        "#  barcodes with >= 500 total UMIs (pipeline config min_umis); detected =",
        "#  nonzero count in >= 1 cell; PAS = genomic 3' end of representative",
        "#  transcript in UTRome GTF (+: end, -: start); chroms UCSC->Ensembl",
        "#  (chrM->MT); mm10 == GRCm38 coordinates.",
        "# Script: scripts/benchmark_tools/scutrquant/gse104556/extract_pas_bed_mouse.py",
        f"# n_isoforms_UTRome_catalog={n_tx} n_isoforms_detected={n_det_filt} (cells>=500 UMIs)",
        f"# n_isoforms_detected_any_barcode={n_det_all} n_cells={n_cells} (of {n_bx} barcodes)",
        "# BED6: chrom, start(0-based), end, name=transcript_id, score=min(1000,round(UMIs)), strand",
    ]
    with open(f"{KDIR}/pas.bed", "w") as fh:
        fh.write("\n".join(hdr_lines) + "\n")
        for r in rows:
            score = min(1000, int(round(r[8])))
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{score}\t{r[4]}\n")

    with open(f"{KDIR}/pas_counts.tsv", "w") as fh:
        fh.write("chrom\tpas_pos_1based\tstrand\ttranscript_id\tgene_id\tgene_name\t"
                 "n_cells_detected\ttotal_umis_cells\tn_barcodes_detected_all\ttotal_umis_all\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[2]}\t{r[4]}\t{r[3]}\t{r[5]}\t{r[6]}\t"
                     f"{r[7]}\t{r[8]:.2f}\t{r[9]}\t{r[10]:.2f}\n")

    with open(f"{KDIR}/pas_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in [("n_barcodes_total", n_bx), ("n_cells_min500umis", n_cells),
                     ("n_isoforms_utrome_catalog", n_tx),
                     ("n_isoforms_detected_cells", n_det_filt),
                     ("n_isoforms_detected_any_barcode", n_det_all),
                     ("total_umis_in_cells", f"{umi_filt[1:].sum():.0f}"),
                     ("min_umis_cell_filter", MIN_UMIS)]:
            fh.write(f"{k}\t{v}\n")
    print(f"[{mouse}] pas.bed: {n_det_filt} sites; cells={n_cells}; "
          f"all-barcode detected={n_det_all}", file=sys.stderr)
