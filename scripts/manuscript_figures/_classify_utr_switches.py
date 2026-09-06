"""Classify each two-site replicated switch by where its sites sit in the annotation.

Tandem 3'UTR APA means both poly(A) sites lie in the SAME 3'UTR of the SAME
transcript.  Anything else is a different event: two last exons, an intronic call,
or a site outside any annotated UTR.
"""
import re
import pandas as pd
import numpy as np

WD = "/mnt/ssd1/Projects/PeakATail_wd"
GTF = f"{WD}/data/references/mouse/Mus_musculus.GRCm38.102.gtf"
REPL = f"{WD}/results/stage3_spermatogenesis_v2/summary/switch_replicated_pas_true.tsv"

# ---- the strict pool: exactly two replicated sites per gene, opposite directions
d = pd.read_csv(REPL, sep="\t")
ns = d.groupby("gene_id").start_m1.nunique()
d2 = d[d.gene_id.isin(ns[ns == 2].index)]
rows = []
for g, s in d2.groupby("gene_id"):
    st = s.strand.iloc[0]
    sites = sorted(s.start_m1.unique(), reverse=(st == "-"))
    prox, dist = sites
    ok = True
    sgn = {}
    for p in sites:
        r = s[s.start_m1 == p]
        if not ((r.dprop_m1 > 0) == (r.dprop_m2 > 0)).all() or r.dprop_m1.gt(0).nunique() != 1:
            ok = False; break
        sgn[p] = r.dprop_m1.iloc[0] > 0
    if not ok or sgn[prox] == sgn[dist]:
        continue
    rows.append(dict(gene=g, chrom=str(s.chrom.iloc[0]), strand=st, prox=int(prox), dist=int(dist),
                     span=abs(int(prox) - int(dist)),
                     direction="shortening" if sgn[prox] else "lengthening",
                     max_abs_dprop=float(s.dprop_m1.abs().max()), qmax=float(s.q_m1.max())))
P = pd.DataFrame(rows)
print(f"two-site, opposite-direction, replicated in both mice: {len(P)} genes")

# ---- annotation: gene spans, and 3'UTR intervals per transcript
genes, utr = [], []
want = set(P.gene)
with open(GTF) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t", 9)
        if f[2] not in ("gene", "three_prime_utr"):
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        if f[2] == "gene":
            nm = re.search(r'gene_name "([^"]+)"', f[8])
            genes.append((f[0], int(f[3]), int(f[4]), f[6], gid, nm.group(1) if nm else gid))
        elif gid in want:
            tx = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            utr.append((gid, tx, f[0], int(f[3]), int(f[4])))
G = pd.DataFrame(genes, columns=["chrom","start","end","strand","gene_id","name"])
U = pd.DataFrame(utr, columns=["gene_id","tx","chrom","start","end"])
print(f"3'UTR intervals loaded for {U.gene_id.nunique()} of them")

# ---- rule 2: annotation-clean window
FLANK = 1500
def clean(c):
    lo, hi = min(c.prox, c.dist) - FLANK, max(c.prox, c.dist) + FLANK
    own = G[G.gene_id == c.gene]
    if not len(own): return False
    o = own.iloc[0]
    if not (o.start <= min(c.prox, c.dist) and max(c.prox, c.dist) <= o.end): return False
    return not len(G[(G.chrom == c.chrom) & (G.end > lo) & (G.start < hi) & (G.gene_id != c.gene)])
P["clean"] = P.apply(clean, axis=1)
C = P[P.clean].copy()
print(f"...also annotation-clean (rule 2): {len(C)} genes")

# ---- the new test: same 3'UTR of the same transcript?
def classify(c):
    u = U[U.gene_id == c.gene]
    if not len(u):
        return "no annotated 3'UTR", ""
    same_iv = u[(u.start <= min(c.prox, c.dist)) & (u.end >= max(c.prox, c.dist))]
    if len(same_iv):
        return "same 3'UTR interval", same_iv.tx.iloc[0]
    hit_p = set(u[(u.start <= c.prox) & (u.end >= c.prox)].tx)
    hit_d = set(u[(u.start <= c.dist) & (u.end >= c.dist)].tx)
    both = hit_p & hit_d
    if both:
        return "same transcript, spliced 3'UTR", sorted(both)[0]
    if hit_p and hit_d:
        return "different transcripts", ""
    if hit_p or hit_d:
        return "only one site in a 3'UTR", ""
    return "neither site in a 3'UTR", ""

res = C.apply(lambda c: pd.Series(classify(c), index=["utr_class", "tx"]), axis=1)
C = pd.concat([C, res], axis=1)
C = C.merge(G[["gene_id","name"]], left_on="gene", right_on="gene_id", how="left")
print("\nWhere do the two sites sit?")
print(C.utr_class.value_counts().to_string())
tand = C[C.utr_class == "same 3'UTR interval"]
print(f"\n=> genuine tandem 3'UTR APA: {len(tand)} genes "
      f"({len(tand)/len(C):.0%} of the annotation-clean set)")
print(f"   direction: {tand.direction.value_counts().to_dict()}")
print(f"   median span {tand.span.median():,.0f} bp, median max|dprop| {tand.max_abs_dprop.median():.2f}")
C.to_csv("/tmp/utr_classified.tsv", sep="\t", index=False)
print("\nOur four drawn genes:")
print(C[C.name.isin(["Rragc","Bpgm","App","Prdm15"])]
      [["name","direction","span","utr_class","tx","max_abs_dprop"]].to_string(index=False))
