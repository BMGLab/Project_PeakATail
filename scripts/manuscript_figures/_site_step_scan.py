"""Does a replicated switch site sit at a real transcript 3' end?

For each sampled site, compare mean read depth in the 200 bp INSIDE the transcript
against the 200 bp OUTSIDE it (strand-aware).  A genuine poly(A) site in 3'-tag data
shows a sharp drop; a call sitting mid-transcript does not.
"""
import re, subprocess, sys
import numpy as np, pandas as pd

WD = "/mnt/ssd1/Projects/PeakATail_wd"
BAM = f"{WD}/data/benchmark/gse104556/starsolo/Mouse1_scRNAseq/Aligned.sortedByCoord.out.bam"
LAB = f"{WD}/results/stage3_spermatogenesis_v2/mouse1/input/stage3_labels.tsv"
REPL = f"{WD}/results/stage3_spermatogenesis_v2/summary/switch_replicated_pas_true.tsv"
W, PAD, N = 200, 600, int(sys.argv[1]) if len(sys.argv) > 1 else 300
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

d = pd.read_csv(REPL, sep="\t")
u = d.drop_duplicates(["gene_id", "start_m1"]).copy()
rng = np.random.default_rng(20260906)
s = u.iloc[rng.choice(len(u), size=min(N, len(u)), replace=False)]
cb = set(pd.read_csv(LAB, sep="\t").cb)
print(f"{len(u):,} distinct replicated (gene, site) pairs; scanning {len(s)}", flush=True)

out = []
for k, (_, r) in enumerate(s.iterrows()):
    pos = int(r.start_m1); lo, hi = pos - PAD, pos + PAD
    v = np.zeros(hi - lo)
    p = subprocess.Popen(["samtools", "view", "-F", "3844", BAM,
                          f"{r.chrom}:{lo+1}-{hi}"], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 11)
        if len(f) < 12: continue
        b = None
        for t in f[11].split("\t"):
            if t.startswith("CB:Z:"): b = t[5:].strip(); break
        if b not in cb: continue
        if ("-" if int(f[1]) & 16 else "+") != r.strand: continue
        rp = int(f[3]) - 1
        for n_, op in CIG.findall(f[5]):
            n_ = int(n_)
            if op in "M=X":
                x0, x1 = max(rp-lo, 0), min(rp+n_-lo, len(v))
                if x1 > x0: v[x0:x1] += 1
                rp += n_
            elif op in "DN": rp += n_
    p.stdout.close(); p.wait()
    i = pos - lo
    ins, o = (v[i-W:i], v[i:i+W]) if r.strand == "+" else (v[i:i+W], v[i-W:i])
    a, b_ = ins.mean(), o.mean()
    out.append(dict(gene=r.gene_id, pos=pos, strand=r.strand, inside=a, outside=b_,
                    step=(a+.01)/(b_+.01), depth=a, dprop=abs(r.dprop_m1)))
    if (k+1) % 60 == 0: print(f"  {k+1}/{len(s)}", flush=True)

S = pd.DataFrame(out)
S.to_csv(f"{WD}/results/figures/manuscript/site_step_scan.tsv", sep="\t", index=False)
cov = S[S.inside >= 1.0]
print(f"\nscanned {len(S)} sites; {len(cov)} have usable inside-coverage (>= 1 read-depth mean)")
for lab, m in (("clear 3' end (step >= 3)", cov.step >= 3),
               ("weak step (1.5 - 3)", (cov.step >= 1.5) & (cov.step < 3)),
               ("no step (< 1.5)", cov.step < 1.5),
               ("inverted (< 0.8, more coverage outside)", cov.step < 0.8)):
    print(f"  {lab:42s} {m.sum():>4}/{len(cov)}  ({m.mean():5.1%})")
print("\nstep quantiles:", cov.step.quantile([.1,.25,.5,.75,.9]).round(1).to_dict())
