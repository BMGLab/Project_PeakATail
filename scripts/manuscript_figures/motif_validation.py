#!/usr/bin/env python3
"""
motif_validation.py -- ATLAS-INDEPENDENT sequence validation of called PAS.

null_control.py showed that precision against the 18.4 M-entry PolyASite dump
is mostly reference density.  This figure asks the question no atlas can rig:
does the GENOMIC SEQUENCE around the inferred cleavage base look like a
polyadenylation site?  Bona fide PAS carry AATAAA/ATTAAA ~10-30 nt upstream of
cleavage in ~70-80% of cases (Tian 2005, Beaudoing/Graber 2000); random
genomic hexamer expectation is a few percent.  No PolyASite input anywhere.

WHAT THE DATA SHOWED (and the figure now says)
    Under the prescribed anchor -- each peak's 3'-most base per strand taken
    as the cleavage site -- the canonical signal is almost absent (7.7% vs
    5.4% in the gene-body null).  Scanning +-150 nt around that anchor finds
    the AATAAA density hump at +40..+90 nt DOWNSTREAM of the peak end
    (mode ~+73), the A-rich signal region peaking ~+98, and the U(T)-rich
    downstream element past +100: the full canonical PAS architecture,
    displaced ~+90 nt.  The called peaks are genuine PAS neighbourhoods
    (AATAAA/ATTAAA within +0..+100 of the peak end: 37% real vs 15% null),
    but the peak 3' end systematically UNDERSHOOTS the cleavage site --
    consistent with 10x R2 coverage dying off before the poly(A) junction.
    An earlier asymmetric-window diagnostic put a phantom hump at -525; that
    was a minus-strand index bug (asymmetric windows break the r = i - FLANK
    mapping under getfasta -s revcomp).  Symmetric windows are immune, and
    are what this script uses everywhere.

SETS   real         runs/grid/lg_annotate/pasbed.bed, 3'-most base per strand
       ip_filtered  runs/grid/lg_ip_filter/pasbed.bed, same collapse
       null         3 seeds of `bedtools shuffle -chrom -noOverlapping
                    -incl <merged gene bodies>` of the real intervals (same
                    N_genic model as null_control.py), then the same collapse.

SEQ    one strand-aware 321-nt window per site via `bedtools getfasta -s`,
       BED [s-160, s+161) around anchor base s; being SYMMETRIC, sequence
       index i maps to transcript-relative position r = i - 160 on both
       strands (r < 0 upstream / 5', r = 0 anchor base, r > 0 downstream).

PANELS a  fraction of sites with AATAAA or ATTAAA (and, in the tsv, each of
          the 12 canonical hexamer variants + their union) fully inside the
          prescribed upstream window r = -40..-5 AND inside r = +0..+100
          where the signal actually sits
       b  positional profile of AATAAA starts, r = -100..+150, 1 bp
          (prescribed range -60..+10 is a subset; its "peak", -8 nt at 0.28%
          of sites, is barely above the null band -- the anchor is wrong)
       c  internal-priming proxy: >=6 consecutive genomic A OR >=70% A in
          r = +10..+30 -- lg_ip_filter halves it vs lg_annotate
       d  single-nucleotide composition, r = -100..+150 (prescribed -50..+50
          is a subset and misses the displaced architecture)

Genome: /home/sharedFolder/humanSTARindex primary assembly (Ensembl GRCh38.99,
chrom naming '1' not 'chr1', matches the pasbed).  Sites on scaffolds absent
from data/references/chrom.sizes.nochr.filt are dropped (same rule as
null_control.py); windows clipped by a chromosome end are dropped and counted.
Everything streams through bedtools/awk; only the ~110 k small windows are
held in memory.  TSV carries every number on the figure and every per-hexamer
count.  Run with LC_ALL=C (set internally for every subprocess).
"""

import os
import subprocess
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- manuscript deliverables: PNG + vector PDF into manuscript/figures/ ---
import os as _os
import matplotlib as _mpl
_mpl.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (journal requirement)
_mpl.rcParams['ps.fonttype'] = 42
FIGDIR = '/mnt/ssd1/Projects/PeakATail_wd/manuscript/figures'
_os.makedirs(FIGDIR, exist_ok=True)
def save_manuscript(fig, name, **kw):
    for ext in ('png', 'pdf'):
        p = _os.path.join(FIGDIR, f'{name}.{ext}')
        fig.savefig(p, **({'dpi': 300} if ext == 'png' else {}), **kw)
        print('wrote', p)


# ----------------------------------------------------------------------------
# paths / constants
# ----------------------------------------------------------------------------
RUNS = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs")
FASTA = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
GENOME = Path("/mnt/ssd1/Projects/PeakATail_wd/data/references/chrom.sizes.nochr.filt")
GENE_END = Path("/mnt/ssd1/Projects/PeakATail_wd/data/references/gene_end.bed")

OUTDIR = Path("/mnt/ssd1/Projects/PeakATail_wd/results/figures/manuscript")
NAME = "motif_validation"
WORK = Path(
    "/mnt/ssd2/claude-tmp/claude-1000/-mnt-ssd1-Projects-PeakATail-wd/"
    "93f9b511-1339-4cbc-9fdb-148fd4b2f08f/scratchpad/motif_run"
)

FLANK = 160                     # window = [s-160, s+161), r = i - FLANK
SIG_WINDOWS = {                 # hexamer must sit FULLY inside [lo, hi]
    "upstream_-40..-5": (-40, -5),      # prescribed: anchor = cleavage site
    "downstream_+0..+100": (0, 100),    # where the signal actually sits
}
PROF_LO, PROF_HI = -100, 150    # panel b hexamer START positions
IP_LO, IP_HI = 10, 30           # panel c downstream window
COMP_LO, COMP_HI = -100, 150    # panel d
N_SEEDS = 3
HEX_OFFSET = 21                 # canonical AATAAA start sits ~15-30 nt
                                # (typ. ~21) upstream of the cleavage site

HEX12 = ["AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", "CATAAA",
         "AATATA", "GATAAA", "AATGAA", "AAGAAA", "ACTAAA", "AATAGA"]

ENV = dict(os.environ, LC_ALL="C")  # a Turkish locale would corrupt sort order

# validated palette -- one colour per ENTITY, fixed across every panel
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
INK, MUTED, GRID, SURFACE = "#1B2429", "#5A6B73", "#D8E0E3", "#FFFFFF"
COL = {"real": C[0], "ip_filtered": C[1], "null": C[2]}
LABEL = {"real": "Real called PAS (lg_annotate)",
         "ip_filtered": "After internal-priming filter (lg_ip_filter)",
         "null": "Null: shuffled in gene bodies (3 seeds)"}
NTCOL = {"A": C[0], "C": C[4], "G": C[3], "T": C[2]}


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, env=ENV,
                          executable="/bin/bash", capture_output=True,
                          text=True).stdout


def log(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------
# 0. workspace + sanity checks
# ----------------------------------------------------------------------------
for d in (WORK, OUTDIR):
    d.mkdir(parents=True, exist_ok=True)

genome_rows = [l.split("\t") for l in GENOME.read_text().strip().split("\n")]
CHROMLEN = {r[0]: int(r[1]) for r in genome_rows}
fai = {l.split("\t")[0]: int(l.split("\t")[1])
       for l in (FASTA.parent / (FASTA.name + ".fai")).read_text().strip().split("\n")}
for c, n in CHROMLEN.items():
    assert c in fai and fai[c] == n, f"fasta/.fai disagrees with chrom sizes at {c}"
assert not any(c.startswith("chr") for c in fai), "fasta uses chr-prefixed names"
log(f"[check] {FASTA.name}: Ensembl naming ('1' not 'chr1'), .fai covers all "
    f"{len(CHROMLEN)} benchmark chroms with matching lengths")


# ----------------------------------------------------------------------------
# 1. point sets (3'-most base per strand) + null shuffles
# ----------------------------------------------------------------------------
def make_point(iv: Path, pt: Path):
    """Collapse each peak to its 3'-most base by strand = assumed cleavage site."""
    sh("awk -F'\\t' 'BEGIN{OFS=\"\\t\"}"
       "{if($6==\"+\"){s=$3-1;e=$3}else{s=$2;e=$2+1} print $1,s,e,$4,$5,$6}' "
       f"{iv} | sort -k1,1 -k2,2n > {pt}")


def prep_real(src: Path, tag: str):
    raw_n = int(sh_out(f"wc -l < {src}").strip())
    iv, pt = WORK / f"{tag}.iv.bed", WORK / f"{tag}.pt.bed"
    if not pt.exists():
        sh(f"awk -F'\\t' 'NR==FNR{{ok[$1]=1;next}} ok[$1]' {GENOME} {src} "
           f"| sort -k1,1 -k2,2n > {iv}")
        make_point(iv, pt)
    n = int(sh_out(f"wc -l < {pt}").strip())
    return dict(tag=tag, raw_n=raw_n, n=n, dropped=raw_n - n, iv=iv, pt=pt)


SETS = {}
SETS["real"] = [prep_real(RUNS / "grid/lg_annotate/pasbed.bed", "real_lg")]
SETS["ip_filtered"] = [prep_real(RUNS / "grid/lg_ip_filter/pasbed.bed", "ip_lg")]
for k in ("real", "ip_filtered"):
    m = SETS[k][0]
    log(f"[set] {k}: n={m['n']:,} ({m['dropped']} off-genome scaffold peaks dropped)")

GENIC = WORK / "genebodies.merged.bed"
if not GENIC.exists():
    sh(f"sort -k1,1 -k2,2n {GENE_END} | bedtools merge -i - > {GENIC}")

SETS["null"] = []
for seed in range(1, N_SEEDS + 1):
    iv, pt = WORK / f"null.s{seed}.iv.bed", WORK / f"null.s{seed}.pt.bed"
    if not pt.exists():
        sh(f"bedtools shuffle -i {SETS['real'][0]['iv']} -g {GENOME} -chrom "
           f"-noOverlapping -maxTries 5000 -seed {seed} -incl {GENIC} "
           f"2>/dev/null | sort -k1,1 -k2,2n > {iv}")
        make_point(iv, pt)
    n = int(sh_out(f"wc -l < {pt}").strip())
    SETS["null"].append(dict(tag=f"null_s{seed}", raw_n=n, n=n, dropped=0,
                             iv=iv, pt=pt))
    log(f"[set] null seed {seed}: n={n:,} (gene-body-constrained shuffle)")


# ----------------------------------------------------------------------------
# 2. strand-aware SYMMETRIC 321-nt windows -> sequences
#    Symmetry is load-bearing: with getfasta -s revcomp, only a window
#    symmetric about the anchor keeps r = i - FLANK true on both strands.
# ----------------------------------------------------------------------------
def get_seqs(meta):
    win = WORK / f"{meta['tag']}.win{FLANK}.bed"
    tab = WORK / f"{meta['tag']}.seq{FLANK}.tab"
    if not tab.exists():
        # drop windows clipped by a chromosome end (bounds from chrom sizes)
        sh(f"awk -F'\\t' -v F={FLANK} 'NR==FNR{{L[$1]=$2;next}} "
           f"$2-F>=0 && $3+F<=L[$1] {{print $1\"\\t\"$2-F\"\\t\"$3+F\"\\t\"$4\"\\t\"$5\"\\t\"$6}}' "
           f"{GENOME} {meta['pt']} > {win}")
        sh(f"bedtools getfasta -fi {FASTA} -bed {win} -s -tab -fo {tab}")
    n_win = int(sh_out(f"wc -l < {win}").strip())
    seqs = []
    with open(tab) as fh:
        for line in fh:
            seqs.append(line.rstrip("\n").split("\t")[1].upper())
    assert len(seqs) == n_win and all(len(s) == 2 * FLANK + 1 for s in seqs)
    meta["n_edge_dropped"] = meta["n"] - n_win
    meta["n_seq"] = n_win
    meta["n_with_N"] = sum("N" in s for s in seqs)
    log(f"[seq] {meta['tag']}: {n_win:,} windows ({meta['n_edge_dropped']} "
        f"clipped at chrom ends dropped, {meta['n_with_N']} contain N)")
    return seqs


for k in SETS:
    for meta in SETS[k]:
        meta["seqs"] = get_seqs(meta)


# ----------------------------------------------------------------------------
# 3. metrics
# ----------------------------------------------------------------------------
def i_of(r):                     # sequence index of relative position r
    return r + FLANK


rows = []


def record(**kw):
    rows.append(kw)


def analyse(meta, set_name, replicate):
    seqs, n = meta["seqs"], meta["n_seq"]
    win_slice = {w: (i_of(lo), i_of(hi) + 1) for w, (lo, hi) in SIG_WINDOWS.items()}
    hexhit = {w: {h: 0 for h in HEX12} for w in SIG_WINDOWS}
    top2 = {w: 0 for w in SIG_WINDOWS}
    any12 = {w: 0 for w in SIG_WINDOWS}
    n_prof = PROF_HI - PROF_LO + 1
    prof = np.zeros(n_prof, dtype=int)             # AATAAA start positions
    pa, pb = i_of(PROF_LO), i_of(PROF_HI)          # allowed start index range
    ia, ib = i_of(IP_LO), i_of(IP_HI) + 1          # panel c window slice
    ip_len = ib - ia
    n_run6 = n_frac70 = n_either = 0
    ca, cb = i_of(COMP_LO), i_of(COMP_HI) + 1      # panel d slice
    comp_len = cb - ca
    counts = {nt: np.zeros(comp_len, dtype=np.int64) for nt in "ACGT"}
    cmat = np.zeros((256, comp_len), dtype=np.int64)

    for s in seqs:
        for w, (a, b) in win_slice.items():
            sub = s[a:b]
            hits = [h for h in HEX12 if h in sub]
            for h in hits:
                hexhit[w][h] += 1
            if "AATAAA" in sub or "ATTAAA" in sub:
                top2[w] += 1
            if hits:
                any12[w] += 1
        j = s.find("AATAAA", pa)
        while j != -1 and j <= pb:
            prof[j - pa] += 1
            j = s.find("AATAAA", j + 1)
        d = s[ia:ib]
        r6 = "AAAAAA" in d
        f70 = d.count("A") / ip_len >= 0.70
        n_run6 += r6
        n_frac70 += f70
        n_either += (r6 or f70)
        row = np.frombuffer(s[ca:cb].encode(), dtype=np.uint8)
        np.add.at(cmat, (row, np.arange(comp_len)), 1)
    for nt in "ACGT":
        counts[nt] = cmat[ord(nt)]

    kw = dict(set=set_name, replicate=replicate, n=n)
    for item, v in (("n_input_bed", meta["raw_n"]),
                    ("n_scaffold_dropped", meta["dropped"]),
                    ("n_edge_dropped", meta["n_edge_dropped"]),
                    ("n_windows_scored", n),
                    ("n_with_N", meta["n_with_N"])):
        record(panel="meta", window="", item=item, position="", count=v,
               value=v, **kw)
    for w in SIG_WINDOWS:
        record(panel="signal_window", window=w, item="AATAAA_or_ATTAAA",
               position="", count=top2[w], value=top2[w] / n, **kw)
        record(panel="signal_window", window=w, item="any_canonical_12",
               position="", count=any12[w], value=any12[w] / n, **kw)
        for h in HEX12:
            record(panel="signal_window", window=w, item=f"hex_{h}",
                   position="", count=hexhit[w][h], value=hexhit[w][h] / n, **kw)
    for j, r in enumerate(range(PROF_LO, PROF_HI + 1)):
        record(panel="profile_AATAAA", window=f"{PROF_LO}..{PROF_HI}",
               item="AATAAA_start", position=r, count=int(prof[j]),
               value=prof[j] / n, **kw)
    for item, v in (("run_A6", n_run6), ("fracA_ge70", n_frac70),
                    ("either", n_either)):
        record(panel="internal_priming", window=f"+{IP_LO}..+{IP_HI}",
               item=item, position="", count=v, value=v / n, **kw)
    valid = sum(counts[nt] for nt in "ACGT")       # per-position non-N total
    for nt in "ACGT":
        for j, r in enumerate(range(COMP_LO, COMP_HI + 1)):
            record(panel="composition", window=f"{COMP_LO}..{COMP_HI}",
                   item=nt, position=r, count=int(counts[nt][j]),
                   value=counts[nt][j] / valid[j], **kw)
    log(f"[metrics] {set_name} rep{replicate}: "
        f"AATAAA|ATTAAA -40..-5 {top2['upstream_-40..-5']/n:.4f}, "
        f"+0..+100 {top2['downstream_+0..+100']/n:.4f}, "
        f"IP-proxy {n_either/n:.4f}")


for set_name, metas in SETS.items():
    for rep, meta in enumerate(metas, start=0 if set_name != "null" else 1):
        analyse(meta, set_name, rep)

df = pd.DataFrame(rows)
df.to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False, float_format="%.6f")
log(f"[out] {OUTDIR / (NAME + '.tsv')}  ({len(df)} rows)")


# ----------------------------------------------------------------------------
# 4. pull numbers for the figure
# ----------------------------------------------------------------------------
def vals(panel, set_name, item, position=None, window=None):
    s = df[(df.panel == panel) & (df.set == set_name) & (df.item == item)]
    if window is not None:
        s = s[s.window == window]
    if position is not None:
        s = s[s.position == position]
    return s["value"].values


def one(panel, set_name, item, window=None):
    v = vals(panel, set_name, item, window=window)
    assert len(v) == 1, (panel, set_name, item, window, len(v))
    return v[0]


def null3(panel, item, position=None, window=None):
    v = vals(panel, "null", item, position, window)
    assert len(v) == N_SEEDS
    return v.mean(), v.min(), v.max()


def profile(set_name, item, positions, panel="profile_AATAAA"):
    s = df[(df.panel == panel) & (df.set == set_name) & (df.item == item)]
    if set_name == "null":
        g = s.groupby("position")["value"]
        agg = pd.concat([g.mean(), g.min(), g.max()], axis=1)
        agg.columns = ["m", "lo", "hi"]
        agg = agg.reindex(positions)
        return agg["m"].values, agg["lo"].values, agg["hi"].values
    m = s.set_index("position")["value"].reindex(positions).values
    return m, None, None


WUP, WDN = "upstream_-40..-5", "downstream_+0..+100"
sig = {w: {k: one("signal_window", k, "AATAAA_or_ATTAAA", w)
           for k in ("real", "ip_filtered")} for w in SIG_WINDOWS}
sig_n = {w: null3("signal_window", "AATAAA_or_ATTAAA", window=w)
         for w in SIG_WINDOWS}
a12 = {w: {k: one("signal_window", k, "any_canonical_12", w)
           for k in ("real", "ip_filtered")} for w in SIG_WINDOWS}
a12_n = {w: null3("signal_window", "any_canonical_12", window=w)
         for w in SIG_WINDOWS}

POS = list(range(PROF_LO, PROF_HI + 1))
prof_real, _, _ = profile("real", "AATAAA_start", POS)
prof_ip, _, _ = profile("ip_filtered", "AATAAA_start", POS)
prof_nm, prof_nlo, prof_nhi = profile("null", "AATAAA_start", POS)
peak_r = POS[int(np.argmax(prof_real))]
# hump extent: contiguous region around the mode where the 11-nt-smoothed real
# profile exceeds twice the null's mean per-position rate
ker = np.ones(11) / 11
sm = np.convolve(prof_real, ker, mode="same")
thr = 2.0 * float(prof_nm.mean())
ix = POS.index(peak_r)
lo_i = ix
while lo_i > 0 and sm[lo_i - 1] >= thr:
    lo_i -= 1
hi_i = ix
while hi_i < len(POS) - 1 and sm[hi_i + 1] >= thr:
    hi_i += 1
hump_lo, hump_hi = POS[lo_i], POS[hi_i]
# the hump is visibly bimodal: report its two highest local maxima (>=10 nt
# apart, raw profile, within the enriched region)
loc = [(v, r) for r, v in zip(POS, prof_real)
       if hump_lo <= r <= hump_hi and v == max(
           prof_real[max(0, POS.index(r) - 5):POS.index(r) + 6])]
loc = sorted(set(loc), reverse=True)
sub_peaks = []
for v, r in loc:
    if all(abs(r - r2) >= 10 for _, r2 in sub_peaks):
        sub_peaks.append((v, r))
    if len(sub_peaks) == 2:
        break
sub_peaks = sorted(r for _, r in sub_peaks)
# prescribed-range "peak" (-60..+10), for the record
sub = [(r, v) for r, v in zip(POS, prof_real) if -60 <= r <= 10]
presc_peak_r, presc_peak_v = max(sub, key=lambda t: t[1])
presc_null_max = max(v for r, v in zip(POS, prof_nhi) if -60 <= r <= 10)

ip_real = one("internal_priming", "real", "either")
ip_ip = one("internal_priming", "ip_filtered", "either")
ip_null = null3("internal_priming", "either")
ip_drop_rel = (ip_real - ip_ip) / ip_real

CPOS = list(range(COMP_LO, COMP_HI + 1))
comp_real = {nt: profile("real", nt, CPOS, "composition")[0] for nt in "ACGT"}
compA_nm, compA_nlo, compA_nhi = profile("null", "A", CPOS, "composition")
a_max_r = CPOS[int(np.argmax(comp_real["A"]))]
a_max_v = float(comp_real["A"].max())

N_REAL = SETS["real"][0]["n_seq"]
N_IP = SETS["ip_filtered"][0]["n_seq"]
N_NULL = SETS["null"][0]["n_seq"]
IMPL_LO, IMPL_HI = peak_r + 15, peak_r + 30      # implied cleavage-site range
                                                 # (canonical spacing off the mode)


# ----------------------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 7.8, "axes.titlesize": 8.6,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 6.8,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.linewidth": 0.8, "font.family": "DejaVu Sans",
})
fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.9))
ax_a, ax_b = axes[0]
ax_c, ax_d = axes[1]


def style(ax, ylab, xlab):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)


GROUPS = ["real", "ip_filtered", "null"]

# --- a. canonical signal fraction, prescribed vs shifted window -------------
xs, w = np.arange(2), 0.24
wlabs = [f"\u2212{-SIG_WINDOWS[WUP][0]}..\u2212{-SIG_WINDOWS[WUP][1]} nt "
         f"upstream of peak end\n(canonical, if peak end were\nthe cleavage site)",
         f"+{SIG_WINDOWS[WDN][0]}..+{SIG_WINDOWS[WDN][1]} nt downstream\n"
         f"(where the signal\nactually sits)"]
for j, g in enumerate(GROUPS):
    if g == "null":
        vv = [sig_n[w_][0] for w_ in (WUP, WDN)]
        lo = [sig_n[w_][1] for w_ in (WUP, WDN)]
        hi = [sig_n[w_][2] for w_ in (WUP, WDN)]
    else:
        vv = [sig[w_][g] for w_ in (WUP, WDN)]
    pos = xs + (j - 1) * (w + 0.035)
    ax_a.bar(pos, vv, width=w, color=COL[g], zorder=3, edgecolor=SURFACE,
             linewidth=0.6, label=LABEL[g])
    if g == "null":
        ax_a.errorbar(pos, vv, yerr=[np.array(vv) - lo, np.array(hi) - vv],
                      fmt="none", ecolor=INK, elinewidth=0.8, capsize=2, zorder=4)
    for p, v in zip(pos, vv):
        ax_a.text(p, v + 0.013, f"{100*v:.0f}%", ha="center", va="bottom",
                  fontsize=6.8, color=COL[g], fontweight="bold")
ax_a.axhspan(0.70, 0.80, color="#E9EEF1", zorder=0)
ax_a.text(0.985, 0.755, "bona fide PAS carry AATAAA/ATTAAA in ~70\u201380%\n"
          "of cases (Tian / Graber atlases); random genomic\n"
          "hexamer expectation is a few % \u2014 the null bars\n"
          "measure it per window directly",
          transform=ax_a.get_yaxis_transform(), fontsize=6.2, color=MUTED,
          va="center", ha="right", linespacing=1.35)
ax_a.set_xticks(xs)
ax_a.set_xticklabels(wlabs, fontsize=6.6)
ax_a.set_ylim(0, 0.92)
style(ax_a, "Fraction of sites with AATAAA or ATTAAA", "")
ax_a.set_title(f"a   Upstream of the peak end the signal is near-null "
               f"({100*sig[WUP]['real']:.0f}% vs {100*sig_n[WUP][0]:.0f}%);\n"
               f"downstream it is real but partial: {100*sig[WDN]['real']:.0f}% "
               f"vs {100*sig_n[WDN][0]:.0f}% null")
ax_a.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False,
            handlelength=1.2, borderpad=0.2, labelspacing=0.3)

# --- b. positional profile of AATAAA ----------------------------------------
Xp = np.array(POS)
ax_b.axvspan(SIG_WINDOWS[WUP][0], SIG_WINDOWS[WUP][1], color="#E9EEF1", zorder=0)
ax_b.fill_between(Xp, 100 * prof_nlo, 100 * prof_nhi, color=COL["null"],
                  alpha=0.30, lw=0)
ax_b.plot(Xp, 100 * prof_nm, color=COL["null"], lw=1.0)
ax_b.plot(Xp, 100 * prof_ip, color=COL["ip_filtered"], lw=1.1, alpha=0.9)
ax_b.plot(Xp, 100 * prof_real, color=COL["real"], lw=1.6, zorder=3)
ax_b.axvline(0, color=MUTED, lw=0.8, ls=":")
ymax = 100 * prof_real.max()
ax_b.annotate(f"mode {peak_r:+d} nt; enriched >2\u00d7 null\n"
              f"across {hump_lo:+d}..{hump_hi:+d} \u2192 implied cleavage\n"
              f"\u2248 +{IMPL_LO}..+{IMPL_HI} nt past the peak end",
              (peak_r, ymax), xytext=(-95, ymax * 0.86), fontsize=6.6,
              color=COL["real"], fontweight="bold", linespacing=1.3,
              arrowprops=dict(arrowstyle="-", color=COL["real"], lw=0.7))
ax_b.text(-22, ymax * 0.31, "shaded: expected AATAAA\nposition if peak end were\n"
          "the cleavage site \u2014 real \u2248 null",
          fontsize=6.0, color=MUTED, ha="center", va="bottom", linespacing=1.3)
ax_b.text(3, ymax * 0.60, "called peak\n3\u2032 end", fontsize=6.2, color=MUTED,
          ha="left")
ax_b.annotate("Real called PAS", (PROF_HI - 2, ymax * 0.62), ha="right",
              color=COL["real"], fontsize=6.8, fontweight="bold")
ax_b.annotate("IP-filtered (tracks real)", (PROF_HI - 2, ymax * 0.54),
              ha="right", color=COL["ip_filtered"], fontsize=6.8,
              fontweight="bold")
ax_b.annotate("Null band (3 seeds)", (PROF_HI - 2, ymax * 0.46), ha="right",
              color=COL["null"], fontsize=6.8, fontweight="bold")
ax_b.set_xlim(PROF_LO, PROF_HI)
style(ax_b, "% of sites with an AATAAA starting here",
      "AATAAA start position relative to peak 3\u2032 end (nt)")
ax_b.set_title(f"b   The AATAAA hump sits downstream of the peak end "
               f"(mode {peak_r:+d} nt) \u2014\nthe called peak stops short "
               f"of the cleavage site")

# --- c. internal-priming proxy ----------------------------------------------
crit = [("run_A6", "\u22656 consecutive A"), ("fracA_ge70", "\u226570% A"),
        ("either", "either (IP proxy)")]
xs3, w3 = np.arange(3), 0.24
for j, g in enumerate(GROUPS):
    if g == "null":
        vv = [null3("internal_priming", c)[0] for c, _ in crit]
        lo = [null3("internal_priming", c)[1] for c, _ in crit]
        hi = [null3("internal_priming", c)[2] for c, _ in crit]
    else:
        vv = [one("internal_priming", g, c) for c, _ in crit]
    pos = xs3 + (j - 1) * (w3 + 0.03)
    ax_c.bar(pos, vv, width=w3, color=COL[g], zorder=3, edgecolor=SURFACE,
             linewidth=0.6, label=LABEL[g])
    if g == "null":
        ax_c.errorbar(pos, vv, yerr=[np.array(vv) - lo, np.array(hi) - vv],
                      fmt="none", ecolor=INK, elinewidth=0.8, capsize=2, zorder=4)
    for p, v in zip(pos, vv):
        ax_c.text(p, v + 0.002, f"{100*v:.1f}%", ha="center", va="bottom",
                  fontsize=6.4, color=COL[g], fontweight="bold")
ax_c.set_xticks(xs3)
ax_c.set_xticklabels([c[1] for c in crit], fontsize=7.2)
ax_c.set_ylim(0, max(ip_real, ip_null[2]) * 1.5)
style(ax_c, f"Fraction flagged at +{IP_LO}..+{IP_HI} nt downstream", "")
verb = "halves" if 0.4 <= ip_drop_rel <= 0.6 else (
    "removes most of" if ip_drop_rel > 0.6 else
    ("reduces" if ip_drop_rel >= 0.05 else "barely changes"))
ax_c.set_title(f"c   The IP filter {verb} A-rich downstream tracts:\n"
               f"{100*ip_real:.1f}% \u2192 {100*ip_ip:.1f}% of sites "
               f"({100*ip_drop_rel:.0f}% relative reduction; null "
               f"{100*ip_null[0]:.1f}%)")
ax_c.legend(loc="upper left", frameon=False, handlelength=1.2, borderpad=0.2,
            labelspacing=0.3)
ax_c.text(0.985, 0.985, "Genomic A-tracts just downstream of an apparent\n"
          "cleavage base are the classic internal-priming\n"
          "signature (oligo-dT annealing inside the transcript).",
          transform=ax_c.transAxes, fontsize=6.2, color=MUTED, va="top",
          ha="right", linespacing=1.35)

# --- d. nucleotide composition ----------------------------------------------
Xc = np.array(CPOS)
ax_d.fill_between(Xc, compA_nlo, compA_nhi, color=COL["null"], alpha=0.22, lw=0)
ax_d.plot(Xc, compA_nm, color=COL["null"], lw=1.0, ls="--")
ends = sorted("ACGT", key=lambda nt: comp_real[nt][-1])
ypos, floor = {}, -1.0
for nt in ends:                       # push labels apart bottom-up
    y = max(float(comp_real[nt][-1]), floor)
    ypos[nt], floor = y, y + 0.028
for nt in "ACGT":
    ax_d.plot(Xc, comp_real[nt], color=NTCOL[nt], lw=1.4, zorder=3)
    ax_d.annotate(nt, (Xc[-1] + 2.5, ypos[nt]), color=NTCOL[nt],
                  fontsize=7.4, fontweight="bold", va="center",
                  annotation_clip=False)
ax_d.axvline(0, color=MUTED, lw=0.8, ls=":")
ax_d.annotate(f"A climbs to {100*a_max_v:.0f}% then falls off a cliff at "
              f"{a_max_r:+d} nt:\nthe modal cleavage position \u2014 genomic\n"
              f"A-richness ends where the poly(A) tail\nwould begin "
              f"(matches panel b: +{IMPL_LO}..+{IMPL_HI})",
              (a_max_r, a_max_v), xytext=(-95, 0.50), fontsize=6.6,
              color=NTCOL["A"], fontweight="bold", linespacing=1.3,
              arrowprops=dict(arrowstyle="-", color=NTCOL["A"], lw=0.7))
ax_d.annotate("T (U in RNA) overtakes past +100:\ndownstream U-rich element",
              (118, comp_real["T"][CPOS.index(118)]),
              xytext=(38, 0.13), fontsize=6.6, color=NTCOL["T"],
              fontweight="bold", linespacing=1.3,
              arrowprops=dict(arrowstyle="-", color=NTCOL["T"], lw=0.7))
ax_d.annotate("null A (gene bodies, dashed;\nband = 3 seeds)",
              (-60, compA_nm[CPOS.index(-60)] + 0.025), fontsize=6.2,
              color=COL["null"], ha="center", va="bottom", linespacing=1.3)
ax_d.text(3, 0.03, "called peak 3\u2032 end", fontsize=6.2, color=MUTED,
          ha="left")
ax_d.set_xlim(COMP_LO, COMP_HI)
ax_d.set_ylim(0, 0.60)
style(ax_d, "Per-base fraction across real called PAS",
      "Position relative to peak 3\u2032 end (nt)")
ax_d.set_title("d   Base composition shows the whole PAS architecture \u2014\n"
               "A-rich signal then U-rich element \u2014 shifted downstream "
               "of the peak end")

fig.suptitle("Sequence, with zero atlas input: called peaks flank genuine "
             "poly(A) sites,\nbut their 3\u2032 ends stop "
             f"~{IMPL_LO}\u2013{IMPL_HI} nt short of the cleavage site",
             x=0.006, ha="left", fontsize=11, fontweight="bold", color=INK,
             y=0.998, linespacing=1.25)
FOOTNOTE = (
    f"Anchor r=0 is each peak's 3\u2032-most base per strand (pasbed.bed carries no finer cleavage estimate); windows are "
    f"strand-aware genomic sequence, bedtools getfasta -s on symmetric \u00b1{FLANK} nt windows (GRCh38 primary assembly, "
    f"Ensembl.99 naming). Sets: real = lg_annotate (n={N_REAL:,}; {SETS['real'][0]['dropped']} scaffold + "
    f"{SETS['real'][0]['n_edge_dropped']} chrom-edge sites dropped), IP-filtered = lg_ip_filter (n={N_IP:,}), null = real "
    f"intervals shuffled within merged gene bodies on the same chromosome (bedtools shuffle -chrom -noOverlapping -incl, "
    f"3 seeds, same n and widths), all collapsed the same way. Hexamers must sit fully inside the stated window; with all 12 "
    f"canonical variants the fractions are {100*a12[WUP]['real']:.0f}% / {100*a12_n[WUP][0]:.0f}% (real/null, upstream "
    f"window) and {100*a12[WDN]['real']:.0f}% / {100*a12_n[WDN][0]:.0f}% (+0..+100; single-hexamer noise is high in a "
    f"101-nt window). The two windows differ in width, so their null expectations differ \u2014 compare each bar with its own "
    f"null. Implied cleavage offset assumes the canonical AATAAA-to-cleavage spacing of 15\u201330 nt; the displacement is "
    f"consistent with 10x R2 coverage ending before the poly(A) junction. In the originally prescribed profile range "
    f"(-60..+10) the real maximum ({presc_peak_r:+d} nt, {100*presc_peak_v:.2f}%) barely clears the null band "
    f"({100*presc_null_max:.2f}%); the downstream hump is bimodal (sub-peaks ~{sub_peaks[0]:+d} and ~{sub_peaks[-1]:+d} nt), "
    f"suggesting at least two offset populations. Even at +0..+100 the AATAAA/ATTAAA rate "
    f"({100*sig[WDN]['real']:.0f}%) stays below the 70\u201380% of curated PAS \u2014 offset variability and residual "
    f"false calls are not separable here. PolyASite is not used anywhere in this figure. All numbers: {NAME}.tsv.")
fig.text(0.006, 0.026, textwrap.fill(FOOTNOTE, 189),
         fontsize=6.3, color=MUTED, va="bottom", ha="left", linespacing=1.5)
fig.tight_layout(rect=[0, 0.125, 1, 0.930])
fig.savefig(OUTDIR / f"{NAME}.png", dpi=300, facecolor=SURFACE)
save_manuscript(fig, NAME, facecolor=SURFACE)
log(f"[out] {OUTDIR / (NAME + '.png')}")


# ----------------------------------------------------------------------------
# 6. console summary
# ----------------------------------------------------------------------------
log("\n================ HEADLINE ================")
for w in (WUP, WDN):
    log(f"AATAAA|ATTAAA {w:22s}: real {sig[w]['real']:.4f} | ip "
        f"{sig[w]['ip_filtered']:.4f} | null {sig_n[w][0]:.4f} "
        f"[{sig_n[w][1]:.4f}-{sig_n[w][2]:.4f}]")
    log(f"any-of-12     {w:22s}: real {a12[w]['real']:.4f} | ip "
        f"{a12[w]['ip_filtered']:.4f} | null {a12_n[w][0]:.4f}")
log(f"AATAAA profile: real enriched(>2x null, 11-nt smooth) {hump_lo:+d}..{hump_hi:+d} nt, mode {peak_r:+d} "
    f"({100*prof_real.max():.2f}%/site/nt) => implied cleavage +{IMPL_LO}..+{IMPL_HI} nt "
    f"past peak end; prescribed -60..+10 'peak' {presc_peak_r:+d} at "
    f"{100*presc_peak_v:.2f}% vs null max {100*presc_null_max:.2f}%")
log(f"AATAAA hump bimodal: sub-peaks at {sub_peaks[0]:+d} and {sub_peaks[-1]:+d} nt")
log(f"internal-priming proxy (+{IP_LO}..+{IP_HI}): real {ip_real:.4f} | "
    f"ip_filtered {ip_ip:.4f} ({100*ip_drop_rel:.1f}% rel. reduction) | "
    f"null {ip_null[0]:.4f} [{ip_null[1]:.4f}-{ip_null[2]:.4f}]")
log(f"composition: real A max {a_max_v:.3f} at {a_max_r:+d} nt (null A "
    f"~{compA_nm.mean():.3f} flat)")
hex_top = sorted(((one('signal_window', 'real', f'hex_{h}', WDN), h)
                  for h in HEX12), reverse=True)[:4]
log("top real hexamers in +0..+100: "
    + ", ".join(f"{h} {v:.3f}" for v, h in hex_top))
log("==========================================")
