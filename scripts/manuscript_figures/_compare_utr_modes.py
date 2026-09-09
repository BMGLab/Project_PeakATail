#!/usr/bin/env python3
"""What does --isoform-agg change about a switch call?

Every switch result in the methods paper was produced at the default
`--isoform-agg per_gene`, where each PAS is tested against the rest of its GENE.
The caller also offers:

  within_utr   each PAS against the other PAS sharing its 3'UTR isoform, i.e.
               tandem-UTR APA -- the event most readers picture
  between_utr  PAS collapsed to 3'UTR-level counts, testing differential 3'UTR
               PREFERENCE between groups (genes with >= 2 UTRs only)

The scope of the denominator is not a cosmetic setting: it decides what "a
switch" means. This quantifies the difference on the same cells, same labels,
same test, changing only that one flag.
"""
import sys
from pathlib import Path

import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
PER_GENE = WD / "results/stage3_spermatogenesis_v2/{m}/diff/true/differential"
ALT = Path(sys.argv[1] if len(sys.argv) > 1 else
           "/mnt/ssd2/claude-tmp/claude-1000/-mnt-ssd1-Projects-PeakATail-wd/"
           "93f9b511-1339-4cbc-9fdb-148fd4b2f08f/scratchpad/utrmode")
PAIRS = ["ES_vs_SPC", "ES_vs_RS", "RS_vs_SPC"]
Q = 0.05


def load(base: Path, m: str) -> pd.DataFrame:
    out = []
    for p in PAIRS:
        f = base / f"fisher_{p}.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        d["pair"] = p
        out.append(d)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def key(d):
    return set(zip(d.gene_id, d.start, d.pair))


rows = []
for m in ("mouse1", "mouse2"):
    sets = {}
    for mode, base in (("per_gene", Path(str(PER_GENE).format(m=m))),
                       ("within_utr", ALT / m / "within_utr" / "differential"),
                       ("between_utr", ALT / m / "between_utr" / "differential")):
        d = load(base, m)
        if not len(d):
            print(f"  {m}/{mode}: no output")
            continue
        hits = d[d.qvalue < Q]
        sets[mode] = key(hits)
        rows.append(dict(mouse=m, mode=mode, tested=len(d), hits=len(hits),
                         genes=hits.gene_id.nunique(),
                         median_abs_dprop=round(hits.delta_proportion.abs().median(), 3)))
    if "per_gene" in sets and "within_utr" in sets:
        a, b = sets["per_gene"], sets["within_utr"]
        print(f"\n{m}: per_gene vs within_utr")
        print(f"  per_gene hits           {len(a):>7,}")
        print(f"  within_utr hits         {len(b):>7,}")
        print(f"  in both                 {len(a & b):>7,}  "
              f"({len(a & b) / max(len(a), 1):.0%} of per_gene)")
        print(f"  per_gene only           {len(a - b):>7,}")
        print(f"  within_utr only         {len(b - a):>7,}")

print("\n" + pd.DataFrame(rows).to_string(index=False))
pd.DataFrame(rows).to_csv(WD / "results/figures/manuscript/utr_mode_comparison.tsv",
                          sep="\t", index=False)
print("\nwrote results/figures/manuscript/utr_mode_comparison.tsv")
