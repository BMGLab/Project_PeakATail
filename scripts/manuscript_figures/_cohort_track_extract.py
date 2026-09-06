#!/usr/bin/env python3
"""Per-cell-type read coverage at replicated cell-type poly(A) switches (lung cohort).

The spermatogenesis companion (_switch_track_extract.py) draws an ordered
developmental axis in two mice.  Here the axis is cell identity, tumour epithelium
against T cells, and the replication unit is the patient: coverage is pooled over
the patients in which the switch replicated, and the per-patient usage that the
test actually scored is carried alongside so the pooling cannot hide a split.

Selection rules are the same two as the mouse figure and are asserted here:
exactly two replicated sites for the gene, and both inside the assigned gene with
no other annotated gene in the drawn window.

Output: results/figures/manuscript/cohort_tracks_{coverage.tsv.gz,meta.tsv,usage.tsv}
"""
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUT = WD / "results/figures/manuscript"
TREE = Path("/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v3")
BAMDIR = Path("/mnt/ssd2/Laugney_Aligned")
GTF = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf")
FLANK = 1500
C1, C2 = "Epithelial_Tumor", "T_cell"        # c1 / c2 exactly as the test named them
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

# 3 in the dominant direction + 1 counter-example.  The clean pool runs 49 to 3, so
# the counter-example is labelled as the rarity it is rather than presented as balance.
GENES = [
    dict(symbol="ZBTB38", gene="ENSG00000177311", chrom="3",  strand="+",
         prox=141442485, dist=141447567, direction="distal-up in tumour"),
    dict(symbol="YWHAB",  gene="ENSG00000166913", chrom="20", strand="+",
         prox=44906531,  dist=44908386,  direction="distal-up in tumour"),
    dict(symbol="CGGBP1", gene="ENSG00000163320", chrom="3",  strand="-",
         prox=88057326,  dist=88051949,  direction="distal-up in tumour"),
    dict(symbol="SPAG1",  gene="ENSG00000104450", chrom="8",  strand="+",
         prox=100241436, dist=100252087, direction="proximal-up in tumour"),
]


def assert_annotation_clean(genes):
    chroms = {g["chrom"] for g in genes}
    rows = []
    with open(GTF) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or f[0] not in chroms:
                continue
            rows.append((f[0], int(f[3]), int(f[4]),
                         re.search(r'gene_id "([^"]+)"', f[8]).group(1)))
    G = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_id"])
    for g in genes:
        s = [g["prox"], g["dist"]]
        lo, hi = min(s) - FLANK, max(s) + FLANK
        own = G[G.gene_id == g["gene"]]
        assert len(own) == 1, (g["symbol"], "gene not in the GTF")
        o = own.iloc[0]
        assert o.start <= min(s) and max(s) <= o.end, (g["symbol"], "site outside its gene")
        other = G[(G.chrom == g["chrom"]) & (G.end > lo) & (G.start < hi)
                  & (G.gene_id != g["gene"])]
        assert not len(other), (g["symbol"], f"{len(other)} other gene(s) in the window")
    print(f"annotation guard: {len(genes)} genes clean")


def main():
    assert_annotation_clean(GENES)
    srr = (pd.read_csv(WD / "data/laughney/20250507_laughney_metadata.csv")
             .rename(columns={"Unnamed: 0": "srr", "Sample.Name": "gsm_short"}))
    runs = srr.groupby("gsm_short").srr.apply(list).to_dict()

    gsms = sorted(p.parent.parent.parent.name for p in
                  TREE.glob("switch/*/true/differential/"
                            f"fisher_{C1}_vs_{C2}.tsv"))
    print(f"GSMs with a {C1} vs {C2} test: {len(gsms)}")

    # per-patient usage, straight from the tables the test itself wrote
    urows = []
    for gsm in gsms:
        d = pd.read_csv(TREE / f"switch/{gsm}/true/differential/fisher_{C1}_vs_{C2}.tsv",
                        sep="\t")
        for g in GENES:
            s = d[(d.gene_id == g["gene"]) & (d.start.isin([g["prox"], g["dist"]]))]
            for _, r in s.iterrows():
                site = "proximal" if r.start == g["prox"] else "distal"
                for cl, npas, ngene in ((C1, r.n_reads_pas_cluster1, r.n_reads_gene_cluster1),
                                        (C2, r.n_reads_pas_cluster2, r.n_reads_gene_cluster2)):
                    urows.append(dict(symbol=g["symbol"], gsm=gsm, celltype=cl, site=site,
                                      usage=npas / ngene if ngene else np.nan,
                                      n_reads_pas=int(npas), n_reads_gene=int(ngene),
                                      qvalue=r.qvalue))
    U = pd.DataFrame(urows)
    U.to_csv(OUT / "cohort_tracks_usage.tsv", sep="\t", index=False)
    print(f"usage rows: {len(U)}  genes covered: {U.symbol.nunique()}")

    # labels: barcode -> cell type, per GSM
    lab = {}
    ncell = {}
    for gsm in gsms:
        L = pd.read_csv(TREE / f"switch/{gsm}/labelled/{gsm}.labels.tsv",
                        sep="\t", index_col=0)
        d = {}
        for bc, ct in L.celltype.items():
            if ct in (C1, C2):
                d[bc.split("_")[-1]] = ct
        lab[gsm] = d
        ncell[gsm] = pd.Series(list(d.values())).value_counts().to_dict()

    rows, meta = [], []
    for g in GENES:
        lo, hi = min(g["prox"], g["dist"]) - FLANK, max(g["prox"], g["dist"]) + FLANK
        width = hi - lo
        # only patients whose test actually scored this gene contribute coverage
        contributing = sorted(set(U[U.symbol == g["symbol"]].gsm))
        acc = {C1: np.zeros(width, np.int64), C2: np.zeros(width, np.int64)}
        cells = {C1: 0, C2: 0}
        used = 0
        for gsm in contributing:
            cells[C1] += ncell[gsm].get(C1, 0)
            cells[C2] += ncell[gsm].get(C2, 0)
            for run in runs.get(gsm.split("-")[0], []):
                bam = BAMDIR / f"{run}_STAR" / f"{run}_Aligned.sortedByCoord.out.bam"
                if not bam.exists():
                    continue
                p = subprocess.Popen(["samtools", "view", "-F", "3844", str(bam),
                                      f"{g['chrom']}:{lo + 1}-{hi}"],
                                     stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
                for line in p.stdout:
                    f = line.split("\t", 11)
                    if len(f) < 12:
                        continue
                    cb = None
                    for tag in f[11].split("\t"):
                        if tag.startswith("CB:Z:"):
                            cb = tag[5:].strip()
                            break
                    ct = lab[gsm].get(cb) if cb else None
                    if ct is None:
                        continue
                    if ("-" if int(f[1]) & 16 else "+") != g["strand"]:
                        continue
                    pos = int(f[3]) - 1
                    a = acc[ct]
                    for n, op in CIG.findall(f[5]):
                        n = int(n)
                        if op in "M=X":
                            x0, x1 = max(pos - lo, 0), min(pos + n - lo, width)
                            if x1 > x0:
                                a[x0:x1] += 1
                            pos += n
                        elif op in "DN":
                            pos += n
                    used += 1
                p.stdout.close()
                p.wait()
        for ct, vec in acc.items():
            if vec.any():
                nz = np.nonzero(vec)[0]
                for i in range(nz[0], nz[-1] + 1):
                    rows.append((g["symbol"], ct, lo + i, int(vec[i])))
            meta.append(dict(symbol=g["symbol"], gene=g["gene"], celltype=ct,
                             n_cells=cells[ct], chrom=g["chrom"], strand=g["strand"],
                             lo=lo, hi=hi, prox=g["prox"], dist=g["dist"],
                             direction=g["direction"], n_patients=len(contributing)))
        print(f"  {g['symbol']:8s} {g['chrom']}:{lo}-{hi}  patients={len(contributing)}  "
              f"reads={used:,}  cells {C1}={cells[C1]:,} {C2}={cells[C2]:,}", flush=True)

    pd.DataFrame(rows, columns=["symbol", "celltype", "pos", "depth"]).to_csv(
        OUT / "cohort_tracks_coverage.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(meta).to_csv(OUT / "cohort_tracks_meta.tsv", sep="\t", index=False)
    print("wrote cohort_tracks_{coverage.tsv.gz,meta.tsv,usage.tsv}")


if __name__ == "__main__":
    main()
