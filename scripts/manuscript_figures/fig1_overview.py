#!/usr/bin/env python3
"""
fig1_overview.py -- manuscript Fig 1 ("fig1_overview"): what PeakATail does,
drawn on the real data it was measured on, and what its pre-registered default
output is.

WHY THIS SCRIPT EXISTS.  Fig 1 is the only main figure with no asset
(21_paper_architecture.md gap 3) and it is the figure every reader looks at
first.  It is therefore built the same way every other manuscript figure is
built: from files, by a script, with every printed number either READ from a
verified artefact or COMPUTED here by a method stated on the figure.  Nothing
on this figure is drawn "for illustration".

WHAT IS REAL ON THIS FIGURE
  * Panel b is one real locus of the PBMC 10k v3 BAM.  The reads are real BAM
    records (CIGAR, soft-clipped sequence, CB/UB tags); the sequence is real
    GRCh38 primary-assembly sequence; the clip positions, the cluster, the
    read-weighted mode and the molecule count are recomputed here from those
    records by a from-scratch re-implementation of ema/countmatrix/polya.py's
    clip_site() + read.py's read_check(), and the molecule count it produces is
    asserted equal to column 5 of the caller's own BED.
  * The lower track of panel b is a real >=2-molecule call the internal-priming
    filter REMOVED, drawn at the same scale and from the same kind of records.
  * Panel d's A-fraction profile is computed here over every kept and every
    removed call (46,544 and 15,588); the precision of the removed set is computed here with
    an independent bedtools pipeline whose value on the kept set reproduces the
    scorer's 0.7062 to 4 decimals (0.7060 -- the 20 off-genome peaks the scorer
    drops).
  * Every P@100 / R_det / n on panels c and e is read from a score_tool.py TSV
    (the single scoring path) or from the Fig 3 trade TSV.

SOURCES OF TRUTH (nothing numeric is typed by hand except where cited)
  19_final_gate_v2.md (FIXED)        pre-registered default 46,524 / 0.7062 /
                                     0.1754; tier-1 IP 0.3520 / 0.2685; nulls
  results/benchmark_tools/.../score_*.tsv   every P, R_det, F1, n on c and e
  results/figures/manuscript/fig3_tradeoff.tsv  the molecule sweep (Fig 3, stem fig3_tradeoff)
  22_performance_roadmap.md          IP rule flags 68,855 / 402,860 peaks (17.09%)
                                     and lifts P 0.5907 -> 0.7062 (shipped rule
                                     vs --ip-rule kinnex, verified negative)
  23_algorithm_roadmap.md section 2 + results/algo_headroom/VERIFY section 5
                                     genome-wide qualifying clip rate 0.5730%
                                     (3,195,067 / 557,564,408 accepted CB reads);
                                     this SUPERSEDES the 1.152% of 10 section R1
  14_switch_calibration_v2.md (SOUND)  the one calibrated switch-test arm
  20_stage3_replication.md           the replication RULE only (its counts are
                                     v1 code and are Fig 5's, not Fig 1's)
  the frozen worktree tools/pa-polya-run-9dfdefb3 and the run's own
  run_config.json                    every parameter printed on the figure

OUTPUTS
  manuscript/figures/fig1_overview.{png,pdf}      (PNG 600 dpi, PDF fonttype 42)
  manuscript/figures/fig1_overview.caption.md     (sidecar: '## Legend' = the
                                                   journal legend, single source
                                                   of the caption; '## Provenance')
  results/figures/manuscript/fig1_overview.tsv            every printed quantity
  results/figures/manuscript/fig1_overview_reads.tsv      every drawn read
  results/figures/manuscript/fig1_overview_sequence.tsv   every drawn base
  results/figures/manuscript/fig1_overview_ipprofile.tsv  panel d profile
  results/figures/manuscript/fig1_overview_trade.tsv      panel e points
  results/figures/manuscript/fig1_overview_work/          cached intermediates

Run:  export LC_ALL=C; python3 scripts/manuscript_figures/fig1_overview.py
      FIG1_REFRESH=1 forces the cached BAM / FASTA / bedtools work to be redone.
"""
import os
import re
import subprocess
import textwrap
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 7.0
matplotlib.rcParams["axes.titlesize"] = 8
matplotlib.rcParams["axes.labelsize"] = 7.2
matplotlib.rcParams["legend.fontsize"] = 6.2

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
H = BT / "pbmc_10k_v3"
IP_ARM = H / "peakatail_clipseeded_final_v2_ipfilt"      # the paper's arm (v2, code 9dfdefb)
NOIP_ARM = H / "peakatail_clipseeded_final_v2"           # same run without --ip-filter
ATLAS = WD / "data/references/atlases/polyasite2.GRCh38.96.rep_sites.bed6"
BAM = WD / "data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam"
FASTA = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
GTF = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf")   # the STAR index annotation
OUTDIR = WD / "results/figures/manuscript"
FIGDIR = WD / "manuscript/figures"
WORK = OUTDIR / "fig1_overview_work"
NAME = "fig1_overview"
for d in (OUTDIR, FIGDIR, WORK):
    d.mkdir(parents=True, exist_ok=True)
REFRESH = os.environ.get("FIG1_REFRESH", "") == "1"

CUTOFF = 100.0
CODE = "9dfdefb"          # frozen worktree tools/pa-polya-run-9dfdefb3; 19 header
RUN_STAMP = "2026-08-21 16:14-16:48"
PREREG = "default + gate committed 2026-08-21 01:18:59 (13 §1, commit 0e27b1a); arms 02:59:36; v2 re-run 16:14"

# Okabe-Ito colourblind-safe palette (manuscript/figures/README.md convention).
INK, MUTED, GRID, PAPER = "#1B2429", "#5A6B73", "#D8E0E3", "#F4F7F8"
BLUE = "#0072B2"     # aligned read / tier 1 / PeakATail
VERM = "#D55E00"     # non-templated poly(A) clip
GREEN = "#009E73"    # kept / pre-registered default output
PURPLE = "#CC79A7"   # internal priming / removed
ORANGE = "#E69F00"   # hexamer
SKY = "#56B4E9"      # tier 2 / coverage-only
GREY = "#999999"

# caller parameters, read from the run's own resolved config below and asserted
PARAMS_EXPECTED = dict(seqlen=91, cb_len=16, barcode_tag="CB", strategy="clip_seeded",
                       polya_min_clip=6, polya_min_purity=0.8, polya_seed_window=25,
                       polya_min_umis=1, polya_clip_filter="none",
                       ip_window_left=10, ip_window_right=30, ip_a_stretch=6, ip_a_fraction=0.7)

# genome-wide clip rate: 23_algorithm_roadmap.md §2 ([V] §5), which supersedes the
# 1.152% of 10_caller_fix_plan.md §R2 (a head-sampling artefact of the QC estimator).
CLIP_READS_GW, CB_READS_GW, CLIP_RATE_GW = 3_195_067, 557_564_408, 0.005730
# 22_performance_roadmap.md §6 / §354 (verified negative on --ip-rule kinnex)
IP_FLAGGED_GW, IP_PEAKS_GW = 68_855, 402_860
# 14_switch_calibration_v2.md (SOUND), arm B0
FDR_NULL_P05, FDR_NULL_RUNS, FDR_NULL_FAMILIES = 0.030, 20, 60

rows_tsv = []          # every printed quantity -> fig1_overview.tsv


def record(panel, quantity, value, source, note=""):
    rows_tsv.append(dict(panel=panel, quantity=quantity, value=value, source=source, note=note))
    return value


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, text=True,
                          capture_output=True, executable="/bin/bash").stdout


# ---------------------------------------------------------------------------
# 1. score TSVs -- the single scoring path (scripts/benchmark_tools/score_tool.py)
# ---------------------------------------------------------------------------
def read_arm(tsv: Path) -> dict:
    df = pd.read_csv(tsv, sep="\t")
    at = df[df.cutoff_bp == CUTOFF]

    def one(panel, series, reference, col="value", replicate=None):
        q = at[(at.panel == panel) & (at.series == series) & (at.reference == reference)]
        if replicate is not None:
            q = q[q.replicate == replicate]
        assert len(q) == 1, (tsv, panel, series, reference, replicate, len(q))
        return float(q[col].iloc[0])

    nulls = [one("precision", "null_genic", "atlas_full", replicate=s) for s in (1, 2, 3)]
    return dict(n=int(one("precision", "real", "atlas_full", col="n_query")),
                P=one("precision", "real", "atlas_full"),
                R_det=one("recall", "real", "atlas_detected"),
                F1_det=one("f1", "real", "atlas_detected"),
                null_P=float(np.mean(nulls)),
                source=str(tsv.relative_to(WD)))


ARMS = {
    "tier2_ip":   IP_ARM / "score_pbmc_final_v2_ipfilt__tier2.tsv",
    "tier1_ip":   IP_ARM / "score_pbmc_final_v2_ipfilt__tier1.tsv",
    "default":    IP_ARM / "score_pbmc_final_v2_ipfilt__PRESPEC_precision_default.tsv",
    "both_ip":    IP_ARM / "score_pbmc_final_v2_ipfilt.tsv",
    "tier1_noip": NOIP_ARM / "score_pbmc_final_v2__tier1.tsv",
    "ge2_noip":   NOIP_ARM / "score_pbmc_final_v2__tier1_ge2umi_POSTHOC.tsv",
    "tier2_noip": NOIP_ARM / "score_pbmc_final_v2__tier2.tsv",
}
S = {k: read_arm(v) for k, v in ARMS.items()}

# --- consistency with the verified gate table (19 §1/§2, verifier verdict FIXED)
assert (S["default"]["n"], round(S["default"]["P"], 4), round(S["default"]["R_det"], 4)) == (46524, 0.7062, 0.1754)
assert (S["tier1_ip"]["n"], round(S["tier1_ip"]["P"], 4), round(S["tier1_ip"]["R_det"], 4)) == (167565, 0.3520, 0.2685)
assert round(S["tier2_ip"]["P"], 4) == 0.0558 and round(S["tier2_noip"]["P"], 4) == 0.0568
assert (S["ge2_noip"]["n"], round(S["ge2_noip"]["P"], 4)) == (62110, 0.5907)
for k in ("default", "tier1_ip", "tier2_ip"):
    record("c", f"{k} P@100", S[k]["P"], S[k]["source"])
    record("c", f"{k} n scored", S[k]["n"], S[k]["source"])
record("c", "no-IP >=2 molecules P@100", S["ge2_noip"]["P"], S["ge2_noip"]["source"],
       "diagnostic arm only -- pas_tier1_ge2mol_noIP is NOT the default (21 §7)")

# --- call-set sizes straight off the BEDs (7-64 lines above the scored n: the
#     scorer drops non-primary contigs; 21 §10 item 1)
BEDS = {
    "peaks_noip": NOIP_ARM / "pas.bed",
    "peaks_ip": IP_ARM / "pas.bed",
    "tier1_ip": IP_ARM / "pas_tier1.bed",
    "tier2_ip": IP_ARM / "pas_tier2.bed",
    "tier1_noip": NOIP_ARM / "pas_tier1.bed",
    "tier2_noip": NOIP_ARM / "pas_tier2.bed",
    "default": IP_ARM / "pas_PRESPEC_precision_default.bed",
    "ge2_noip": NOIP_ARM / "pas_tier1_ge2umi_POSTHOC.bed",
}
NB = {k: int(sh(f"wc -l < {v}").strip()) for k, v in BEDS.items()}
assert NB["peaks_noip"] == IP_PEAKS_GW, NB          # 402,860 == 22's denominator
assert NB["peaks_noip"] - NB["peaks_ip"] == IP_FLAGGED_GW  # 68,855 flagged == 22
assert NB["tier1_ip"] + NB["tier2_ip"] == NB["peaks_ip"]
for k, v in NB.items():
    record("c", f"BED rows {k}", v, str(BEDS[k].relative_to(WD)))

# --- the run's own resolved parameters (so the figure cannot drift from the run)
import json
cfg = json.load(open(IP_ARM / "run/run_config.json"))
P_RUN = dict(seqlen=cfg["variables"]["seqlen"], cb_len=cfg["variables"]["cb_len"],
             barcode_tag=cfg["variables"]["barcode_tag"],
             **{k: cfg["args"][k] for k in
                ("strategy", "polya_min_clip", "polya_min_purity", "polya_seed_window",
                 "polya_min_umis", "polya_clip_filter", "ip_window_left", "ip_window_right",
                 "ip_a_stretch", "ip_a_fraction")})
assert P_RUN == PARAMS_EXPECTED, P_RUN
SEQLEN = P_RUN["seqlen"]
for k, v in P_RUN.items():
    record("all", f"run parameter {k}", v, "results/.../peakatail_clipseeded_final_v2_ipfilt/run/run_config.json")

# ---------------------------------------------------------------------------
# 2. panel e: the molecule sweep, from the Fig 3 TSV (one scorer, same denominators)
# ---------------------------------------------------------------------------
TRADE_SRC = OUTDIR / "fig3_tradeoff.tsv"
tr = pd.read_csv(TRADE_SRC, sep="\t")
tr = tr[(tr.panel == "a") & (tr.cutoff_bp == 100)].copy()
tr["min_molecules"] = tr["min_molecules"].astype(int)
trade = tr[["dataset", "min_molecules", "pre_registered", "n", "atlas_agreement_precision",
            "recall_detected_genes", "F1_detected_genes"]].sort_values(["dataset", "min_molecules"])
assert set(trade.dataset) == {"pbmc", "mouse1", "mouse2"}
_p = trade[(trade.dataset == "pbmc") & (trade.min_molecules == 2)].iloc[0]
assert (int(_p.n), round(_p.atlas_agreement_precision, 4)) == (46524, 0.7062)
_p1 = trade[(trade.dataset == "pbmc") & (trade.min_molecules == 1)].iloc[0]
assert round(_p1.F1_detected_genes, 4) == 0.3046 and round(_p.F1_detected_genes, 4) == 0.2811

# ---------------------------------------------------------------------------
# 3. what the internal-priming filter removes  (computed here)
# ---------------------------------------------------------------------------
REMOVED = WORK / "ip_removed.bed"
if REFRESH or not REMOVED.exists():
    sh(f"sort -k1,1 -k2,2n {BEDS['ge2_noip']} > {WORK}/noip_ge2.sorted.bed")
    sh(f"sort -k1,1 -k2,2n {BEDS['default']} > {WORK}/default.sorted.bed")
    sh(f"bedtools intersect -a {WORK}/noip_ge2.sorted.bed -b {WORK}/default.sorted.bed -s -v > {REMOVED}")
N_REMOVED = int(sh(f"wc -l < {REMOVED}").strip())
assert N_REMOVED == NB["ge2_noip"] - NB["default"] == 15588, N_REMOVED   # strict subset

PRECFILE = WORK / "ip_precision.tsv"
if REFRESH or not PRECFILE.exists():
    sh(f"sort -k1,1 -k2,2n {ATLAS} > {WORK}/atlas.sorted.bed")
    out = []
    for tag, bed in (("kept_default", WORK / "default.sorted.bed"),
                     ("removed_by_ip", REMOVED),
                     ("noip_ge2", WORK / "noip_ge2.sorted.bed")):
        n = int(sh(f"wc -l < {bed}").strip())
        m = int(sh(f"sort -k1,1 -k2,2n {bed} | bedtools window -a - -b {WORK}/atlas.sorted.bed "
                   f"-w 100 -sm -u | wc -l").strip())
        out.append(f"{tag}\t{n}\t{m}\t{m / n:.6f}")
    PRECFILE.write_text("set\tn\tn_matched_100bp\tprecision\n" + "\n".join(out) + "\n")
ipprec = pd.read_csv(PRECFILE, sep="\t").set_index("set")
P_REMOVED = float(ipprec.loc["removed_by_ip", "precision"])
# the same independent pipeline on the kept set must reproduce the scorer's 0.7062
assert abs(ipprec.loc["kept_default", "precision"] - S["default"]["P"]) < 0.001
assert abs(ipprec.loc["noip_ge2", "precision"] - S["ge2_noip"]["P"]) < 0.001
record("d", "precision@100 of the calls the IP filter removes", P_REMOVED,
       "computed here: bedtools window -w 100 -sm -u vs PolyASite 2.0 rep sites")
record("d", "precision@100 of the kept default, same independent pipeline",
       float(ipprec.loc["kept_default", "precision"]),
       "computed here; reproduces score_tool.py's 0.7062 (the 20 off-genome peaks it drops)")

# ---------------------------------------------------------------------------
# 4. genomic sequence helpers (samtools faidx; no pysam on this box)
# ---------------------------------------------------------------------------
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def faidx(regions, tag):
    """Fetch many regions at once; returns a list of upper-case sequences."""
    rf = WORK / f"regions_{tag}.txt"
    rf.write_text("\n".join(regions) + "\n")
    out = subprocess.run(["samtools", "faidx", str(FASTA), "-r", str(rf)],
                         capture_output=True, text=True, check=True).stdout
    seqs, buf, started = [], [], False
    for line in out.splitlines():
        if line.startswith(">"):
            if started:
                seqs.append("".join(buf))
            buf, started = [], True
        else:
            buf.append(line.strip())
    seqs.append("".join(buf))
    assert len(seqs) == len(regions), (len(seqs), len(regions))
    return [s.upper() for s in seqs]


def a_profile(bed: Path, tag: str, flank=50):
    """Mean A-fraction by transcript-oriented offset from the called position.

    Offset 0 is the last aligned base (the BED position); the caller's
    internal-priming window is offsets -9..+30 (ip_window() uses BED end on '+'
    and BED start on '-', i.e. pas_pos = site + 1 in 0-based coordinates).
    """
    cache = WORK / f"aprofile_{tag}.npz"
    if not REFRESH and cache.exists():
        z = np.load(cache)
        return z["prof"], int(z["n"])
    rows = [l.split("\t") for l in bed.read_text().splitlines()]
    regs, strands = [], []
    for r in rows:
        chrom, s, e, st = r[0], int(r[1]), int(r[2]), r[5]
        pos1 = e if st == "+" else s + 1        # 1-based coordinate of the called base
        if pos1 - flank < 1:
            continue
        regs.append(f"{chrom}:{pos1 - flank}-{pos1 + flank}")
        strands.append(st)
    seqs = faidx(regs, tag)
    n = 2 * flank + 1
    acc, cnt = np.zeros(n), 0
    for sq, st in zip(seqs, strands):
        if len(sq) != n:
            continue
        if st == "-":
            sq = sq.translate(RC)[::-1]
        acc += np.frombuffer(sq.encode(), dtype=np.uint8) == ord("A")
        cnt += 1
    prof = acc / cnt
    np.savez(cache, prof=prof, n=cnt)
    return prof, cnt


PROF_KEPT, N_PROF_KEPT = a_profile(BEDS["default"], "kept")
PROF_REM, N_PROF_REM = a_profile(REMOVED, "removed")
OFF = np.arange(-50, 51)


def ip_triggers(bed: Path, tag: str, flank=50):
    """Which of the two IP rules fires, and where the flagged A-run sits."""
    cache = WORK / f"iptriggers_{tag}.tsv"
    if not REFRESH and cache.exists():
        return pd.read_csv(cache, sep="\t").iloc[0].to_dict()
    rows = [l.split("\t") for l in bed.read_text().splitlines()]
    regs, strands = [], []
    for r in rows:
        chrom, s, e, st = r[0], int(r[1]), int(r[2]), r[5]
        pos1 = e if st == "+" else s + 1
        if pos1 - flank < 1:
            continue
        regs.append(f"{chrom}:{pos1 - flank}-{pos1 + flank}")
        strands.append(st)
    seqs = faidx(regs, tag + "_trig")
    a6 = frac = both = n = 0
    ends = []
    for sq, st in zip(seqs, strands):
        if len(sq) != 2 * flank + 1:
            continue
        if st == "-":
            sq = sq.translate(RC)[::-1]
        win = sq[flank - 9:flank + 31]           # offsets -9..+30 == the caller's window
        has6 = "AAAAAA" in win
        hasf = win.count("A") / len(win) >= P_RUN["ip_a_fraction"]
        a6 += has6
        frac += hasf
        both += has6 and hasf
        n += 1
        if has6:
            m = max(re.finditer(r"A{6,}", win), key=lambda mm: len(mm.group()))
            ends.append(m.end() - 1 - 9)
    ends = np.array(ends)
    d = dict(n=n, a6_rule=a6, a_fraction_rule=frac, both=both,
             a6_only=a6 - both, fraction_only=frac - both,
             run_ends_at_or_before_call=float((ends <= 0).mean()),
             median_run_end_offset=float(np.median(ends)))
    pd.DataFrame([d]).to_csv(cache, sep="\t", index=False)
    return d


TRIG = ip_triggers(REMOVED, "removed")
assert TRIG["n"] == N_REMOVED - 0 or TRIG["n"] <= N_REMOVED
record("d", "removed calls flagged by the >=6-A run rule", TRIG["a6_rule"] / TRIG["n"],
       "computed here on the caller's window: offsets -9..+30 from the called base "
       "(--ip-window-left 10 / --ip-window-right 30 are measured from the BED end = call + 1)")
record("d", "removed calls flagged by the >=70% A-fraction rule", TRIG["a_fraction_rule"] / TRIG["n"],
       "computed here")
record("d", "flagged A-run ends at or before the call", TRIG["run_ends_at_or_before_call"],
       "computed here")

# ---------------------------------------------------------------------------
# 5. one real locus: reads straight out of the BAM
#    from-scratch re-implementation of read.py:read_check + polya.py:clip_site
# ---------------------------------------------------------------------------
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_span(cigar):
    return sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")


def clip_site(flag, cigar, seq, start1, min_clip, min_purity):
    """ema/countmatrix/polya.py:clip_site, re-implemented (see module docstring)."""
    ops = CIG.findall(cigar)
    if not ops:
        return None, 0, ""
    if flag & 16:                                    # reverse: clip at the START, T-run
        n, op = int(ops[0][0]), ops[0][1]
        if op != "S" or n < min_clip:
            return None, 0, ""
        clip = seq[:n]
        if clip.count("T") / n < min_purity:
            return None, 0, ""
        if len(clip) - len(clip.rstrip("T")) < min_clip:   # run flush to the alignment edge
            return None, 0, ""
        return start1 - 1, n, clip
    n, op = int(ops[-1][0]), ops[-1][1]              # forward: clip at the END, A-run
    if op != "S" or n < min_clip:
        return None, 0, ""
    clip = seq[-n:]
    if clip.count("A") / n < min_purity:
        return None, 0, ""
    if len(clip) - len(clip.lstrip("A")) < min_clip:
        return None, 0, ""
    return start1 - 1 + ref_span(cigar) - 1, n, clip


def locus_reads(chrom, site0, strand, flank, tag):
    """Every read read_check() would accept on this strand near the site."""
    cache = WORK / f"reads_{tag}.tsv"
    if not REFRESH and cache.exists():
        return pd.read_csv(cache, sep="\t", keep_default_na=False)
    reg = f"{chrom}:{max(1, site0 + 1 - flank)}-{site0 + 1 + flank}"
    out = subprocess.run(["samtools", "view", str(BAM), reg],
                         capture_output=True, text=True, check=True).stdout
    recs = []
    for line in out.splitlines():
        f = line.split("\t")
        flag, cigar, seq, start1 = int(f[1]), f[5], f[9], int(f[3])
        if flag & 4:                                     # unmapped
            continue
        if (strand == "+") == bool(flag & 16):           # wrong strand for this pass
            continue
        cb = ub = ""
        for t in f[11:]:
            if t.startswith("CB:Z:"):
                cb = t[5:]
            elif t.startswith("UB:Z:"):
                ub = t[5:]
        if not cb:
            continue
        if len(cb) != P_RUN["cb_len"]:                   # strip the CellRanger GEM suffix
            dash = cb.rfind("-")
            if dash == P_RUN["cb_len"] and cb[dash + 1:].isdigit():
                cb = cb[:dash]
            if len(cb) != P_RUN["cb_len"]:
                continue
        span = ref_span(cigar)
        if span > SEQLEN:                                # read.py's span rule
            continue
        cs, clen, cseq = clip_site(flag, cigar, seq, start1,
                                   P_RUN["polya_min_clip"], P_RUN["polya_min_purity"])
        recs.append(dict(start0=start1 - 1, span=span, cigar=cigar, cb=cb, ub=ub,
                         clip_site=-1 if cs is None else cs, clip_len=clen, clip_seq=cseq,
                         f3844=int((flag & 3844) == 0)))
    df = pd.DataFrame(recs)
    df.to_csv(cache, sep="\t", index=False)
    return df


LOCI = {
    # the kept example: a pre-registered-default call, THAP3 3' end.  Chosen by
    # stated criteria (see the caption): '+' strand, 10-40 molecules, isolated
    # (nearest other call in this arm 608 bp away), a strand-matched atlas
    # representative site 1 bp away, a canonical AATAAA in -40..-5, IP-pass.
    "kept": dict(chrom="1", site0=6633561, strand="+", pas_id=817,
                 label="1:6,633,562 (+)  ·  THAP3 3′ end  ·  tier 1, 31 molecules, IP-pass, in the default",
                 bed=BEDS["default"]),
    # the removed example: a >=2-molecule tier-1 call the IP filter deleted.
    # Criteria: '+' strand, >=8 molecules, no strand-matched atlas site within
    # 200 bp (nearest is 16,389 bp), a genomic A-run ending at the call.
    "removed": dict(chrom="16", site0=22493756, strand="+", pas_id=112359,
                    label="16:22,493,757 (+)  ·  19 molecules  ·  flagged as internal priming and removed",
                    bed=REMOVED),
}
for key, L in LOCI.items():
    row = sh(f"awk -F'\\t' '$1==\"{L['chrom']}\" && $2=={L['site0']}' {L['bed']}").strip().split("\t")
    assert row and int(row[3]) == L["pas_id"], (key, row)
    L["molecules_bed"] = int(row[4])
    L["reads"] = locus_reads(L["chrom"], L["site0"], L["strand"], 300, key)
    cl = L["reads"][L["reads"].clip_site >= 0]
    cl = cl[(cl.clip_site - L["site0"]).abs() <= P_RUN["polya_seed_window"]]
    L["clipped"] = cl
    L["molecules_recomputed"] = len(set(zip(cl.cb, cl.ub)))
    # the re-implementation must reproduce the caller's own support column
    assert L["molecules_recomputed"] == L["molecules_bed"], (key, L["molecules_recomputed"], L["molecules_bed"])
    L["n_usable"] = len(L["reads"])
    L["seq"] = faidx([f"{L['chrom']}:{L['site0'] + 1 - 80}-{L['site0'] + 1 + 60}"], key + "_seq")[0]
    L["seq_lo"] = -80                                   # index 0 of L["seq"] is offset -80
    record("b", f"{key} locus: usable reads in +/-300 bp", L["n_usable"], "computed here from the BAM")
    record("b", f"{key} locus: qualifying clip reads at the site", len(cl), "computed here from the BAM")
    record("b", f"{key} locus: distinct (CB,UB) molecules", L["molecules_recomputed"],
           "computed here; equals column 5 of the caller's BED")

# ---------------------------------------------------------------------------
# 5b. how many calls satisfy the criteria each drawn locus was selected from
#     (so the examples are quantified, not asserted)
# ---------------------------------------------------------------------------
def selection_counts():
    cache = WORK / "selection_counts.tsv"
    if not REFRESH and cache.exists():
        return pd.read_csv(cache, sep="\t").iloc[0].to_dict()
    # --- kept-locus criteria, over the whole default call set
    sh(f"sort -k1,1 -k2,2n {BEDS['peaks_ip']} > {WORK}/ip_all.sorted.bed")
    iso = sh(f"bedtools closest -a {WORK}/default.sorted.bed -b {WORK}/ip_all.sorted.bed -io -d -t first")
    atl = sh(f"bedtools closest -a {WORK}/default.sorted.bed -b {WORK}/atlas.sorted.bed -s -d -t first")
    keep = []
    for li, la in zip(iso.splitlines(), atl.splitlines()):
        fi, fa = li.split("\t"), la.split("\t")
        assert (fi[0], fi[1]) == (fa[0], fa[1])
        if int(fi[4]) >= 10 and int(fi[-1]) > 300 and int(fa[-1]) <= 25:
            keep.append((fi[0], int(fi[1]), int(fi[2]), fi[5]))
    regs, strands = [], []
    for chrom, s, e, st in keep:
        pos1 = e if st == "+" else s + 1
        if pos1 - 80 < 1:
            continue
        regs.append(f"{chrom}:{pos1 - 80}-{pos1 + 60}")
        strands.append(st)
    n_kept = 0
    for sq, st in zip(faidx(regs, "selkept"), strands):
        if len(sq) != 141:
            continue
        if st == "-":
            sq = sq.translate(RC)[::-1]
        if "AATAAA" in sq[80 - 40:80 - 4] and "AAAAAA" not in sq[80 - 9:80 + 31]:
            n_kept += 1
    # --- removed-locus criteria, over every call the IP filter removed
    atl_r = sh(f"bedtools closest -a {REMOVED} -b {WORK}/atlas.sorted.bed -s -d -t first")
    cand = []
    for line in atl_r.splitlines():
        f = line.split("\t")
        if int(f[4]) >= 8 and int(f[-1]) > 200:
            cand.append((f[0], int(f[1]), int(f[2]), f[5]))
    regs, strands = [], []
    for chrom, s, e, st in cand:
        pos1 = e if st == "+" else s + 1
        if pos1 - 80 < 1:
            continue
        regs.append(f"{chrom}:{pos1 - 80}-{pos1 + 60}")
        strands.append(st)
    n_rem = 0
    for sq, st in zip(faidx(regs, "selrem"), strands):
        if len(sq) != 141:
            continue
        if st == "-":
            sq = sq.translate(RC)[::-1]
        win = sq[80 - 9:80 + 31]
        runs = [mm for mm in re.finditer(r"A{8,}", win) if mm.end() - 1 - 9 <= 0 and mm.end() - 1 - 9 >= -2]
        if runs:
            n_rem += 1
    d = dict(kept_criteria_n=n_kept, kept_denominator=NB["default"],
             removed_criteria_n=n_rem, removed_denominator=N_REMOVED)
    pd.DataFrame([d]).to_csv(cache, sep="\t", index=False)
    return d


# gene context of each drawn locus, from the same annotation the run used
GENE_BED = WD / "data/references/gene_end.bed"


def tx_three_prime(chrom, site0, strand, gene, tag):
    """Nearest ANNOTATED TRANSCRIPT 3' end of *gene* (Ensembl 99 -- the GTF the
    STAR index was built from), measured from the called base.

    This is measured, not asserted: the gene BODY of a multi-transcript gene can
    end kilobases past a real internal/proximal PAS, so "the gene's 3' end" is
    not a safe label for a call that merely lies inside the gene (caveat 4)."""
    cache = WORK / f"tx3p_{tag}.tsv"
    if REFRESH or not cache.exists():
        raw = sh(f"grep -F 'gene_name \"{gene}\"' {GTF} | "
                 f"awk -F'\\t' '$1==\"{chrom}\" && $3==\"transcript\" && $7==\"{strand}\"'")
        rows = []
        for line in raw.splitlines():
            f = line.split("\t")
            three1 = int(f[4]) if strand == "+" else int(f[3])
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            rows.append(dict(transcript=m.group(1), three_prime_1based=three1,
                             distance=abs(three1 - (site0 + 1))))
        pd.DataFrame(rows).sort_values("distance").to_csv(cache, sep="\t", index=False)
    d = pd.read_csv(cache, sep="\t")
    assert len(d), (gene, chrom)
    return d.iloc[0].to_dict()


def atlas_distance(chrom, site0, strand, tag):
    """bp to the nearest STRAND-MATCHED PolyASite 2.0 representative site."""
    cache = WORK / f"atlasdist_{tag}.txt"
    if REFRESH or not cache.exists():
        if not (WORK / "atlas.sorted.bed").exists():
            sh(f"sort -k1,1 -k2,2n {ATLAS} > {WORK}/atlas.sorted.bed")
        b = WORK / f"locus_{tag}.bed"
        b.write_text(f"{chrom}\t{site0}\t{site0 + 1}\t{tag}\t0\t{strand}\n")
        out = sh(f"bedtools closest -a {b} -b {WORK}/atlas.sorted.bed -s -d -t first")
        cache.write_text(out.strip().split("\t")[-1] + "\n")
    return int(cache.read_text().strip())


for key, L in LOCI.items():
    hits = sh(f"awk -F'\\t' '$1==\"{L['chrom']}\" && $2<={L['site0']} && $3>{L['site0']} "
              f"&& $6==\"{L['strand']}\"' {GENE_BED}").strip().splitlines()
    L["genes_same_strand"] = [r.split("\t")[4] for r in hits if r]
    record("b", f"{key} locus: same-strand gene bodies overlapping it",
           ",".join(L["genes_same_strand"]) or "none", str(GENE_BED.relative_to(WD)))
    L["gene"] = L["genes_same_strand"][0]
    L["gene_body_end"] = int(sh(f"awk -F'\\t' '$5==\"{L['gene']}\"{{print ($6==\"+\")?$3:$2+1}}' "
                                f"{GENE_BED}").strip().splitlines()[0])
    L["tx3p"] = tx_three_prime(L["chrom"], L["site0"], L["strand"], L["gene"], key)
    L["atlas_bp"] = atlas_distance(L["chrom"], L["site0"], L["strand"], key)
    record("b", f"{key} locus: nearest annotated transcript 3' end of {L['gene']} (bp)",
           int(L["tx3p"]["distance"]),
           f"computed here from {GTF} ({L['tx3p']['transcript']})")
    record("b", f"{key} locus: distance to the {L['gene']} GENE-body 3' end (bp)",
           abs(L["gene_body_end"] - (L["site0"] + 1)), str(GENE_BED.relative_to(WD)),
           "the gene body ends far past the call: 'gene 3' end' is not a safe label")
    record("b", f"{key} locus: nearest strand-matched atlas representative site (bp)",
           L["atlas_bp"], "computed here: bedtools closest -s -d against PolyASite 2.0")

# labels are BUILT from the measured facts, never typed
LOCI["kept"]["label"] = (
    f"{LOCI['kept']['chrom']}:{LOCI['kept']['site0'] + 1:,} ({LOCI['kept']['strand']})  ·  in "
    f"{LOCI['kept']['gene']}  ·  tier 1, {LOCI['kept']['molecules_bed']} "
    f"molecules, IP-pass, in the default")
LOCI["removed"]["label"] = (
    f"{LOCI['removed']['chrom']}:{LOCI['removed']['site0'] + 1:,} ({LOCI['removed']['strand']})"
    f"  ·  {LOCI['removed']['molecules_bed']} molecules  ·  flagged as internal priming and removed")

SEL = selection_counts()
record("b", "default calls satisfying the kept-example criteria", SEL["kept_criteria_n"],
       "computed here: >=10 molecules, nearest other call >300 bp, atlas site <=25 bp, "
       "canonical AATAAA in -40..-5, IP-pass")
record("b", "removed calls satisfying the removed-example criteria", SEL["removed_criteria_n"],
       "computed here: >=8 molecules, no strand-matched atlas site within 200 bp, "
       "flagged A-run >=8 nt ending at -2..0")

# clip-site composition of the kept locus (this is the clustering input)
KEPT = LOCI["kept"]
clip_counts = (KEPT["clipped"].clip_site - KEPT["site0"]).value_counts().sort_index()
MODE_OFF = int(clip_counts.idxmax())
assert MODE_OFF == 0, clip_counts       # the caller's read-weighted mode is the BED position
SPAN_CLIP = int(clip_counts.index.max() - clip_counts.index.min())
record("b", "kept locus: distinct clip positions in the cluster", len(clip_counts), "computed here")
record("b", "kept locus: span of the cluster (bp)", SPAN_CLIP, "computed here")

# coverage on the caller's convention: every accepted read extended to seq_len
def coverage(L, lo, hi):
    cov = np.zeros(hi - lo)
    for s in L["reads"].start0.values:
        a, b = max(int(s), lo), min(int(s) + SEQLEN, hi)
        if b > a:
            cov[a - lo:b - lo] += 1
    return cov


COV_LO, COV_HI = KEPT["site0"] - 300, KEPT["site0"] + 120
COV = coverage(KEPT, COV_LO, COV_HI)
COV_X = np.arange(COV_LO, COV_HI) - KEPT["site0"]
COV_MAX_OFF = int(COV_X[int(COV.argmax())])
record("b", "kept locus: coverage maximum (reads, seq_len-extended)", float(COV.max()), "computed here")
record("b", "kept locus: offset of the coverage maximum (bp)", COV_MAX_OFF, "computed here")
COV_AT_SITE = float(COV[KEPT["site0"] - COV_LO])
_cov_site = COV_AT_SITE
record("b", "kept locus: coverage at the cleavage site", COV_AT_SITE, "computed here")

# the count column the switch test actually uses (pas_support.tsv, written by the run)
sup = pd.read_csv(IP_ARM / "run/pas_support.tsv", sep="\t").set_index("pas_id")
SUP = sup.loc[KEPT["pas_id"]]
assert int(SUP.clip_reads) == len(KEPT["clipped"]) and int(SUP.clip_umis) == KEPT["molecules_bed"]
WINDOW_READS = int(SUP.window_reads)
record("b", "kept locus: reads counted into the matrix column", WINDOW_READS,
       "run/pas_support.tsv (window_reads)",
       "NOT simply the cleavage window: polya.py flush() gives the cluster the coverage candidate "
       "it suppresses (_assign_candidate) plus count_ends() over [site - seq_len, site + 25] minus "
       "those candidates, clipped at midpoints to neighbouring PAS")
# the cleavage window ALONE does not produce that number: in polya.py flush() a
# tier-1 cluster also takes the per-CB counts of every coverage candidate it
# suppresses (_assign_candidate), and only the residue of [site-seq_len,
# site+25] outside those candidates and inside its midpoint territory is added
# by count_ends().  Measured here so panel a/f cannot overstate the definition.
_ends_ext = KEPT["reads"].start0.values + np.maximum(KEPT["reads"].span.values, SEQLEN)
WINDOW_ENDS_ONLY = int(((_ends_ext >= KEPT["site0"] - SEQLEN)
                        & (_ends_ext <= KEPT["site0"] + 25)).sum())
record("b", "kept locus: accepted read ends inside [site - seq_len, site + 25] alone",
       WINDOW_ENDS_ONLY, "computed here from the BAM",
       f"vs window_reads {WINDOW_READS}: the cleavage window alone accounts for at most "
       f"{WINDOW_ENDS_ONLY / WINDOW_READS:.0%} of the count column; the rest is the coverage "
       f"candidate the cluster suppresses (polya.py flush)")

# hexamer position in the kept locus (canonical AATAAA in -40..-5, transcript orientation)
_up = KEPT["seq"][80 - 40:80 - 4]
_k = _up.rfind("AATAAA")
assert _k >= 0
HEX_OFF = _k - 40
record("b", "kept locus: canonical AATAAA start offset", HEX_OFF, "computed here from GRCh38")

# the removed locus: the genomic A-run the rule fires on
REM = LOCI["removed"]
_win = REM["seq"][80 - 9:80 + 31]
_m = max(re.finditer(r"A{6,}", _win), key=lambda mm: len(mm.group()))
REM_RUN = (len(_m.group()), _m.start() - 9, _m.end() - 1 - 9)
record("b", "removed locus: flagged A-run (length, start, end offsets)", str(REM_RUN),
       "computed here from GRCh38")

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 6. figure
#
# Page-realistic geometry: the canvas is 7.35 x 8.30 in (~the width the figure
# occupies on a journal page), so the point sizes here are the point sizes in
# print.  The former on-figure footer caption moved to the sidecar Legend at
# the 2026-09-02 submission pass and the canvas height shrank by the freed
# space (10.40 -> 8.30).  The width stays 7.35 in: panel a's box text and the
# base-resolution sequence panels do not tolerate narrowing to 180 mm without
# crowding.  Nothing is placed by eye: every panel is an
# axes positioned in INCHES from the top-left corner, every stack of text is
# spaced by lineh() (which turns a font size into axes-fraction for that
# panel's own height), and every wrapped string is wrapped to chars(), which
# turns the panel's own width into a character count.  Layout is therefore
# stable when a number changes length.
# ---------------------------------------------------------------------------
FIG_W, FIG_H = 7.35, 8.30
ML, MR = 0.42, 0.14
CW = FIG_W - ML - MR
fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor("white")


def axbox(x0, top, w, h, label=None):
    """Axes placed by inches from the top-left corner of the canvas."""
    ax = fig.add_axes([x0 / FIG_W, (FIG_H - top - h) / FIG_H, w / FIG_W, h / FIG_H], label=label)
    ax.set_facecolor("none")
    ax._h_in = h
    ax._w_in = w
    return ax


def plain(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2.5, labelsize=5.8)
    ax.set_axisbelow(True)
    return ax


def lineh(ax, pt, mult=1.30):
    """One line of pt-size text, in axes-fraction units of this axes."""
    return pt * mult / (ax._h_in * 72.0)


def chars(w_in, pt):
    """How many characters of pt-size DejaVu Sans fit in w_in inches."""
    return max(20, int(w_in * 72.0 / (pt * 0.545)))


def panel_tag(ax, letter, title, dy=0.012):
    ax.text(0.0, 1.0 + dy, letter, transform=ax.transAxes, fontsize=9.0, fontweight="bold",
            color=INK, ha="left", va="bottom")
    ax.text(0.034, 1.0 + dy, title, transform=ax.transAxes, fontsize=7.2, color=INK,
            ha="left", va="bottom")


def rbox(ax, x, y, w, h, fc="white", ec=MUTED, lw=0.8, z=3, r=0.010):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z,
                                transform=ax.transData, mutation_aspect=1))


def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=0.9, z=5, ms=5.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
                                 linewidth=lw, color=color, zorder=z, shrinkA=0, shrinkB=0))


# ===========================================================================
# panel a -- ingestion: what makes a read usable, and the channel it seeds from
# ===========================================================================
axA = plain(axbox(ML, 0.24, CW, 0.82, "a"))
panel_tag(axA, "a", "Read ingestion — what the caller accepts, and the small, highly specific channel it seeds from")
LA = lineh(axA, 5.5)
BOX_T, BOX_B = 0.95, 0.42
stages = [
    (0.000, 0.200, BLUE, "10x BAM",
     ["PBMC 10k v3 · CellRanger 3.0.0", "3′ v3 · 91 bp R2 · CB + UB tags"], None),
    (0.240, 0.240, BLUE, "usable read",
     ["mapped · strand of this pass", f"CB present, {P_RUN['cb_len']} nt · span ≤ {SEQLEN} nt"],
     f"{CB_READS_GW:,} accepted CB reads"),
    (0.520, 0.265, VERM, "qualifying poly(A) clip",
     [f"3′-side soft clip ≥ {P_RUN['polya_min_clip']} nt, ≥ {P_RUN['polya_min_purity']:.0%} A (+) / T (−)",
      f"with a ≥ {P_RUN['polya_min_clip']}-nt run flush to the edge"],
     f"{CLIP_READS_GW:,} = {CLIP_RATE_GW:.3%} of them  (23 §2)"),
    (0.825, 0.175, GREEN, "cleavage site + support",
     ["site = last aligned base", "support = distinct (CB, UB)"],
     "duplicates share a key"),
]
for x, w, col, head, subs, note in stages:
    rbox(axA, x, BOX_B, w, BOX_T - BOX_B, ec=col, lw=1.0)
    axA.text(x + w / 2, BOX_T - 0.05 - LA * 0.5, head, ha="center", va="center",
             fontsize=6.2, color=col, fontweight="bold")
    for i, s in enumerate(subs):
        axA.text(x + w / 2, BOX_T - 0.05 - LA * (1.95 + 1.15 * i), s, ha="center", va="center",
                 fontsize=5.4, color=INK)
    if note:
        axA.text(x + w / 2, BOX_B - 0.05, note, ha="center", va="top", fontsize=5.4,
                 color=VERM if col == VERM else MUTED)
for x0, x1 in ((0.200, 0.240), (0.480, 0.520), (0.785, 0.825)):
    arrow(axA, x0 + 0.004, (BOX_T + BOX_B) / 2, x1 - 0.004, (BOX_T + BOX_B) / 2, lw=1.0, ms=6)
# (the coverage-channel / count-column note that stood here moved to the
#  sidecar Legend, panel-a paragraph, at the 2026-09-02 submission pass)

# ===========================================================================
# panel b -- one real locus, to scale
# ===========================================================================
BX0, BX1 = -45.5, 32.5              # base-level window, offsets from the called base

# ---- b0: context strip ----------------------------------------------------
axB0 = style(axbox(ML, 1.14, CW, 0.48, "b0"))
axB0.set_xlim(-300, 120)
axB0.fill_between(COV_X, 0, COV, color=SKY, alpha=0.30, lw=0, zorder=2)
axB0.plot(COV_X, COV, color=SKY, lw=0.7, zorder=3)
axB0.axvline(0, color=VERM, lw=0.9, ls=(0, (3, 2)), zorder=4)
axB0.set_ylim(0, COV.max() * 1.50)
axB0.set_yticks([0, 1000, 2000])
axB0.set_ylabel("reads", fontsize=5.8, color=MUTED, labelpad=1)
axB0.set_xlabel("nt from the called cleavage site", fontsize=5.8, color=MUTED, labelpad=0.5)
panel_tag(axB0, "b", f"One real locus, drawn to scale — {KEPT['label']}", dy=0.22)
axB0.add_patch(Rectangle((BX0, 0), BX1 - BX0, COV.max() * 1.50, facecolor=ORANGE, alpha=0.13,
                         edgecolor="none", zorder=1))
# (the transcript-3'-end / gene-body-end distances and the "a coverage summit
#  is not one" qualification moved to the sidecar Legend, panel-b paragraph)
axB0.text(-296, COV.max() * 1.47,
          f"coverage of the {KEPT['n_usable']:,} accepted reads within ±300 bp, each extended to {SEQLEN} nt",
          fontsize=5.3, color=MUTED, ha="left", va="top")
axB0.annotate(f"local maximum {int(COV.max()):,} reads at {COV_MAX_OFF:+d} bp",
              xy=(COV_MAX_OFF, COV.max()), xytext=(-105, COV.max() * 1.18),
              fontsize=5.3, color=MUTED, ha="left", va="center",
              arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
axB0.annotate(f"{int(_cov_site)} reads at the cleavage site",
              xy=(3, _cov_site), xytext=(10, COV.max() * 1.30),
              fontsize=5.3, color=INK, ha="left", va="top", linespacing=1.3,
              arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
axB0.text(BX1 - 1, COV.max() * 0.06, "zoom", fontsize=5.2, color=ORANGE, ha="right", va="bottom")

# ---- b1: base-level zoom of the kept call ---------------------------------
B1_TOP, B1_H = 1.84, 1.44
axB1 = plain(axbox(ML, B1_TOP, CW, B1_H, "b1"))
axB1.set_xlim(BX0, BX1)

# zoom connectors from the context strip (figure coordinates)
for xa, xb in ((BX0, 0.0), (BX1, 1.0)):
    fig.add_artist(Line2D([(ML + CW * (xa + 300) / 420) / FIG_W, (ML + CW * xb) / FIG_W],
                          [(FIG_H - 1.62) / FIG_H, (FIG_H - B1_TOP) / FIG_H],
                          transform=fig.transFigure, color=ORANGE, lw=0.6, alpha=0.45, zorder=0))


def draw_sequence(ax, L, y, h, fs=5.2, lo=BX0, hi=BX1):
    """Real GRCh38 sequence; the A's are picked out because the A's are the point."""
    seq, s0 = L["seq"], L["seq_lo"]
    for off in range(int(np.ceil(lo)), int(np.floor(hi)) + 1):
        i = off - s0
        if not (0 <= i < len(seq)):
            continue
        b = seq[i]
        if b == "A":
            ax.add_patch(Rectangle((off - 0.5, y), 1.0, h, facecolor=VERM, alpha=0.18,
                                   edgecolor="none", zorder=2))
        ax.text(off, y + h / 2, b, ha="center", va="center", fontsize=fs,
                family="DejaVu Sans Mono", color=VERM if b == "A" else MUTED, zorder=3)


def draw_ruler(ax, y, h, fs=5.0):
    for off in range(-40, 31, 10):
        ax.add_line(Line2D([off, off], [y + h, y], color=GRID, lw=0.6, zorder=2))
        ax.text(off, y - 0.008, f"{off:+d}", ha="center", va="top", fontsize=fs, color=MUTED)


def draw_read(ax, x_start, x_end, clip_len, clip_seq, y, h, letters_from=12):
    """One BAM record, to scale: aligned bases from its own start to its OWN last
    aligned base (*x_end*, which is where clip_site() put this read's clip and is
    NOT always the cluster's called base), then the soft clip."""
    a = max(x_start, BX0)
    edge = x_end + 0.5                          # right edge of this read's last aligned base
    ax.add_patch(Rectangle((a - 0.5, y), edge - (a - 0.5), h, facecolor=BLUE, alpha=0.62,
                           edgecolor="none", zorder=3))
    cl = min(int(clip_len), int(BX1 - edge))
    if cl > 0:
        ax.add_patch(Rectangle((edge, y), cl, h, facecolor=VERM, edgecolor="none", alpha=0.92, zorder=4))
        if clip_len >= letters_from:
            ax.text(edge + cl / 2, y + h / 2, clip_seq[:cl], ha="center", va="center",
                    fontsize=5.0, family="DejaVu Sans Mono", color="white", zorder=5)
        if clip_len > cl:                       # clip runs past the drawn window
            ax.text(edge + cl + 0.4, y + h / 2, "»", ha="left", va="center", fontsize=5.0,
                    color=VERM, zorder=5)
    if x_start < BX0:
        ax.text(BX0 + 0.4, y + h / 2, "«", ha="left", va="center", fontsize=5.0, color="white", zorder=5)


IP_LO = -P_RUN["ip_window_left"] + 1 - 0.5      # offsets -9 .. +30 (see a_profile docstring)
IP_HI = P_RUN["ip_window_right"] + 0.5

# --- reads: a spread over the real clip lengths, plus reads with no clip
cl = KEPT["clipped"].drop_duplicates(subset=["cb", "ub"]).sort_values(["clip_len", "start0"])
sel = np.unique(np.linspace(0, len(cl) - 1, 6).round().astype(int))
show = cl.iloc[sel].sort_values("clip_len", ascending=False)
un = KEPT["reads"][KEPT["reads"].clip_site < 0].copy()
un["end_off"] = un.start0 + un.span - 1 - KEPT["site0"]
un = un[(un.end_off > BX0 + 8) & (un.end_off < -6)].sort_values("end_off")
un_show = un.iloc[np.unique(np.linspace(0, len(un) - 1, 2).round().astype(int))]

drawn_reads = []
axB1.text(BX0 + 0.5, 1.005, f"{len(show)} of the {len(KEPT['clipped'])} qualifying clip reads here "
          f"({KEPT['molecules_recomputed']} distinct molecules)", fontsize=5.6, color=BLUE,
          ha="left", va="top")
axB1.text(1.4, 1.005, "cleavage site", fontsize=5.6, color=VERM, ha="left", va="top")

Y_TOP, RH, PITCH = 0.905, 0.038, 0.054
y = Y_TOP
for _, r in show.iterrows():
    a = r.start0 - KEPT["site0"]
    draw_read(axB1, a, int(r.clip_site - KEPT["site0"]), int(r.clip_len), r.clip_seq, y, RH)
    drawn_reads.append(dict(locus="kept", start_offset=int(a), aligned_span=int(r.span),
                            clip_len=int(r.clip_len), clip_seq=r.clip_seq, cb=r.cb, ub=r.ub,
                            cigar=r.cigar))
    y -= PITCH
for j, (_, r) in enumerate(un_show.iterrows()):
    a = r.start0 - KEPT["site0"]
    axB1.add_patch(Rectangle((max(a, BX0) - 0.5, y), r.end_off - max(a, BX0) + 1, RH,
                             facecolor=GREY, alpha=0.40, edgecolor="none", zorder=3))
    if a < BX0:
        axB1.text(BX0 + 0.4, y + RH / 2, "«", ha="left", va="center", fontsize=5.0,
                  color="white", zorder=5)
    if j == len(un_show) - 1:                    # label inside the bar: no room beside it
        axB1.text((BX0 + r.end_off) / 2, y + RH / 2, "no qualifying clip", ha="center",
                  va="center", fontsize=5.0, color=INK, zorder=5)
    drawn_reads.append(dict(locus="kept", start_offset=int(a), aligned_span=int(r.span),
                            clip_len=0, clip_seq="", cb=r.cb, ub=r.ub, cigar=r.cigar))
    y -= PITCH

# --- clip-site stems: the clustering input
Y_STEM = y - 0.032
STEM_MAX = 0.078
for off, k in clip_counts.items():
    axB1.add_line(Line2D([off, off], [Y_STEM, Y_STEM + STEM_MAX * k / clip_counts.max()],
                         color=VERM, lw=1.4, zorder=4, solid_capstyle="butt"))
    axB1.text(off, Y_STEM + STEM_MAX * k / clip_counts.max() + 0.006, str(int(k)),
              ha="center", va="bottom", fontsize=5.0, color=VERM)
_b0 = clip_counts.index.min() - 0.5
_b1 = clip_counts.index.max() + 0.5      # the cluster's OWN extent: single linkage
                                          # bounds the GAP between neighbours, not the span
axB1.add_line(Line2D([_b0, _b1], [Y_STEM - 0.022, Y_STEM - 0.022], color=INK, lw=0.8, zorder=4))
for xx in (_b0, _b1):
    axB1.add_line(Line2D([xx, xx], [Y_STEM - 0.032, Y_STEM - 0.012], color=INK, lw=0.8, zorder=4))
axB1.text(BX1 - 0.5, Y_STEM + STEM_MAX * 0.45, "clip reads\nper position", fontsize=5.3,
          color=VERM, ha="right", va="center", linespacing=1.3)
# (the single-linkage sentence moved to the sidecar Legend; a short label
#  naming the bracket stays on the image)
axB1.text(BX0 + 0.5, Y_STEM - 0.042,
          f"cluster (single linkage, gap ≤ {P_RUN['polya_seed_window']} bp): "
          f"{len(clip_counts)} positions over {SPAN_CLIP} bp — mode = the call",
          ha="left", va="top", fontsize=5.3, color=INK, linespacing=1.35)

# --- sequence, hexamer, internal-priming window, ruler
Y_SEQ = Y_STEM - 0.305
H_SEQ = 0.062
axB1.add_patch(Rectangle((IP_LO, Y_SEQ - 0.070), IP_HI - IP_LO, Y_STEM - Y_SEQ + 0.185,
                         facecolor=PURPLE, alpha=0.09, edgecolor=PURPLE, lw=0.6,
                         ls=(0, (2, 2)), zorder=1))
axB1.text(IP_HI - 0.5, Y_STEM + 0.118, f"internal-priming window "
          f"({IP_LO + 0.5:+.0f} … {IP_HI - 0.5:+.0f} nt from the call = {int(IP_HI - IP_LO)} nt of genomic sequence)",
          ha="right", va="bottom", fontsize=5.3, color=PURPLE)
draw_sequence(axB1, KEPT, Y_SEQ, H_SEQ)
axB1.add_patch(Rectangle((HEX_OFF - 0.5, Y_SEQ - 0.008), 6, H_SEQ + 0.016, facecolor="none",
                         edgecolor=ORANGE, lw=1.0, zorder=5))
axB1.text(HEX_OFF + 3, Y_SEQ + H_SEQ + 0.016, f"hexamer AATAAA ({HEX_OFF:+d})",
          ha="center", va="bottom", fontsize=5.2, color=ORANGE)
draw_ruler(axB1, Y_SEQ - 0.040, 0.018)
axB1.axvline(0.5, color=VERM, lw=0.9, ls=(0, (3, 2)), zorder=6)
_kw = KEPT["seq"][80 - 9:80 + 31]
axB1.text(BX0 + 0.5, Y_SEQ - 0.098,
          f"strand-matched atlas representative site {KEPT['atlas_bp']} bp away",
          ha="left", va="top", fontsize=5.5, color=MUTED)
axB1.text(BX1 - 0.5, Y_SEQ - 0.098,
          f"longest genomic A-run in the window "
          f"{max(len(m.group()) for m in re.finditer('A+', _kw))} nt (flagged at ≥ {P_RUN['ip_a_stretch']}) · "
          f"A-fraction {_kw.count('A') / 40:.0%} (flagged at ≥ {P_RUN['ip_a_fraction']:.0%})  →  IP-PASS, kept",
          ha="right", va="top", fontsize=5.5, color=GREEN)

# ---- legend strip ---------------------------------------------------------
axLg = plain(axbox(ML, 3.34, CW, 0.12, "blegend"))
leg = [Line2D([], [], color=BLUE, lw=4, alpha=.62, label="aligned bases"),
       Line2D([], [], color=VERM, lw=4, label="non-templated poly(A) soft clip"),
       Line2D([], [], color=GREY, lw=4, alpha=.40, label="read with no qualifying clip"),
       Line2D([], [], color=ORANGE, lw=1.1, label="canonical hexamer"),
       Line2D([], [], color=PURPLE, lw=1.1, ls=(0, (2, 2)), label="internal-priming window"),
       Line2D([], [], marker="s", color=VERM, alpha=0.35, lw=0, ms=4, label="genomic A")]
axLg.legend(handles=leg, loc="center", ncol=6, frameon=False, handlelength=1.4,
            columnspacing=1.0, handletextpad=0.35, fontsize=5.3)

# ---- b2: the same geometry for a call the filter removed ------------------
axB2 = plain(axbox(ML, 3.58, CW, 0.96, "b2"))
axB2.set_xlim(BX0, BX1)
axB2.text(BX0 + 0.5, 1.10, REM["label"], fontsize=6.2, color=PURPLE, ha="left", va="top",
          fontweight="bold")

cl2 = REM["clipped"].drop_duplicates(subset=["cb", "ub"]).sort_values("clip_len", ascending=False)
sel2 = np.unique(np.linspace(0, len(cl2) - 1, 4).round().astype(int))
cl2 = cl2.iloc[sel2].sort_values("clip_len", ascending=False)
# ("read-level evidence indistinguishable from the call above" moved to the
#  sidecar Legend, panel-b paragraph)
axB2.text(BX0 + 0.5, 0.985, f"{len(sel2)} of {len(REM['clipped'])} clip reads here",
          fontsize=5.3, color=INK, ha="left", va="top")
Y2_TOP, RH2, PITCH2 = 0.855, 0.058, 0.076
y = Y2_TOP
for _, r in cl2.iterrows():
    a = r.start0 - REM["site0"]
    draw_read(axB2, a, int(r.clip_site - REM["site0"]), int(r.clip_len), r.clip_seq, y, RH2)
    drawn_reads.append(dict(locus="removed", start_offset=int(a), aligned_span=int(r.span),
                            clip_len=int(r.clip_len), clip_seq=r.clip_seq, cb=r.cb, ub=r.ub,
                            cigar=r.cigar))
    y -= PITCH2
Y_SEQ2 = y - 0.120
axB2.add_patch(Rectangle((IP_LO, Y_SEQ2 - 0.075), IP_HI - IP_LO, Y2_TOP - Y_SEQ2 + 0.155,
                         facecolor=PURPLE, alpha=0.09, edgecolor=PURPLE, lw=0.6,
                         ls=(0, (2, 2)), zorder=1))
axB2.axvline(0.5, color=VERM, lw=0.9, ls=(0, (3, 2)), zorder=6)
draw_sequence(axB2, REM, Y_SEQ2, 0.090)
axB2.add_patch(Rectangle((REM_RUN[1] - 0.5, Y_SEQ2 - 0.010), REM_RUN[0], 0.110, facecolor="none",
                         edgecolor=PURPLE, lw=1.1, zorder=5))
axB2.text(REM_RUN[1] + REM_RUN[0] / 2, Y_SEQ2 + 0.105,
          f"genomic A-run of {REM_RUN[0]} nt ending at the call — where oligo-dT primes",
          ha="center", va="bottom", fontsize=5.3, color=PURPLE)
draw_ruler(axB2, Y_SEQ2 - 0.055, 0.024)
axB2.text(BX1 - 0.5, Y_SEQ2 - 0.130,
          f"≥ {P_RUN['ip_a_stretch']} consecutive genomic A in the window  →  FLAGGED, removed; "
          f"nearest strand-matched atlas site {REM['atlas_bp']:,} bp away, nearest {REM['gene']} "
          f"transcript 3′ end {int(REM['tx3p']['distance']):,} bp",
          ha="right", va="top", fontsize=5.5, color=PURPLE)

# ===========================================================================
# panel c -- the two tiers, the filter, the threshold
# ===========================================================================
CW_C = 3.50
axC = plain(axbox(ML, 4.72, CW_C, 1.62, "c"))
panel_tag(axC, "c", "From peaks to the pre-registered default (PBMC 10k v3)")
LC = lineh(axC, 6.0)
BH = 3.6 * LC


def node(x, w, ytop, col, title, n, P=None, note=None, big=7.6):
    rbox(axC, x, ytop - BH, w, BH, ec=col, lw=1.0)
    axC.text(x + w / 2, ytop - LC * 0.75, title, ha="center", va="center", fontsize=5.8,
             color=col, fontweight="bold")
    axC.text(x + w / 2, ytop - LC * 1.80, n, ha="center", va="center", fontsize=big, color=INK)
    if P is not None:                       # a to-scale precision bar inside the node
        bw = w * 0.48
        bx = x + 0.02
        by = ytop - LC * 3.10
        axC.add_patch(Rectangle((bx, by), bw, LC * 0.45, facecolor=PAPER, edgecolor=GRID, lw=0.4,
                                zorder=4))
        axC.add_patch(Rectangle((bx, by), bw * P, LC * 0.45, facecolor=col, alpha=0.85,
                                edgecolor="none", zorder=5))
        axC.add_line(Line2D([bx + bw * S["default"]["null_P"]] * 2, [by, by + LC * 0.45],
                            color=MUTED, lw=0.7, zorder=6))
        axC.text(bx + bw + 0.012, by + LC * 0.22, f"P@100 {P:.3f}", ha="left", va="center",
                 fontsize=5.2, color=MUTED)
    elif note:
        axC.text(x + w / 2, ytop - LC * 2.85, note, ha="center", va="center", fontsize=5.3, color=MUTED)


node(0.02, 0.72, 1.00, MUTED, "peak candidates", f"{NB['peaks_noip']:,}",
     note="clip clusters + coverage summits")
axC.text(0.245, 1.00 - BH - LC * 0.50, f"internal-priming filter: −{IP_FLAGGED_GW:,} "
         f"({IP_FLAGGED_GW / IP_PEAKS_GW:.1%})", ha="left", va="center", fontsize=5.2, color=PURPLE)
Y2 = 1.00 - BH - LC * 1.00
node(0.02, 0.36, Y2, BLUE, "tier 1 · clip-supported", f"{NB['tier1_ip']:,}", P=S["tier1_ip"]["P"], big=7.2)
node(0.54, 0.36, Y2, SKY, "tier 2 · coverage-only", f"{NB['tier2_ip']:,}", P=S["tier2_ip"]["P"], big=7.2)
Y3 = Y2 - BH - LC * 1.00
node(0.02, 0.36, Y3, GREEN, "≥ 2 clip molecules", f"{NB['default']:,}", P=S["default"]["P"], big=7.2)
arrow(axC, 0.20, 1.00 - BH, 0.20, Y2 + 0.004, ms=5.5)
arrow(axC, 0.72, 1.00 - BH, 0.72, Y2 + 0.004, ms=5.5)
arrow(axC, 0.20, Y2 - BH, 0.20, Y3 + 0.004, ms=5.5)
axC.text(0.42, Y2 - BH - LC * 0.30, f"discards {NB['tier1_ip'] - NB['default']:,} single-molecule\n"
         f"sites ({1 - NB['default'] / NB['tier1_ip']:.1%} of tier 1, IP arm)",
         ha="left", va="top", fontsize=5.3, color=MUTED, linespacing=1.35)
YB = Y3 - BH - LC * 0.95
CC = chars(CW_C, 5.4)
# (the gate-PASS sentence, the null-ratio sentence, the tier-2-file note and
#  the BED-rows-vs-scored-n note moved to the sidecar Legend, panel-c
#  paragraph; the definition of the default and a short label for the grey
#  null tick stay on the image)
for i, (txt, col, fs, bold) in enumerate([
        ("= the pre-registered default: tier 1 ∩ IP-pass ∩ ≥ 2 distinct clip molecules", GREEN, 5.6, True),
        (f"grey tick = 3-seed gene-body-shuffled null P@100 ({S['default']['null_P']:.4f})", MUTED, 5.4, False)]):
    axC.text(0.0, YB - LC * 1.10 * i, textwrap.fill(txt, CC), ha="left", va="top", fontsize=fs,
             color=col, fontweight="bold" if bold else "normal")

# ===========================================================================
# panel d -- the internal-priming filter, measured on the whole call set
# ===========================================================================
DX, DW = ML + 3.78, CW - 3.78
axD = style(axbox(DX, 4.88, DW, 0.95, "d"))
panel_tag(axD, "d", "What the filter removes, measured", dy=0.30)
axD.add_patch(Rectangle((IP_LO, 0), IP_HI - IP_LO, 1.05, facecolor=PURPLE, alpha=0.08,
                        edgecolor="none", zorder=1))
axD.plot(OFF, PROF_KEPT, color=GREEN, lw=1.0, zorder=4, label=f"kept ({N_PROF_KEPT:,})")
axD.plot(OFF, PROF_REM, color=PURPLE, lw=1.0, zorder=4, label=f"removed ({N_PROF_REM:,})")
axD.axvline(0, color=VERM, lw=0.8, ls=(0, (3, 2)), zorder=3)
axD.set_xlim(-50, 50)
axD.set_ylim(0, 1.05)
axD.set_xticks([-50, -25, 0, 25, 50])
axD.set_yticks([0, 0.5, 1.0])
axD.set_xlabel("nt from the call", fontsize=5.8, color=MUTED, labelpad=0.5)
axD.set_ylabel("genomic A-fraction", fontsize=5.8, color=MUTED, labelpad=1)
axD.grid(color=GRID, lw=0.5, alpha=0.7)
axD.legend(loc="upper left", frameon=False, fontsize=5.2, handlelength=1.1, borderpad=0.05,
           labelspacing=0.2)
axD.annotate("hexamer", xy=(HEX_OFF + 3, PROF_KEPT[50 + HEX_OFF + 3]), xytext=(-47, 0.70),
             fontsize=5.2, color=ORANGE, ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.5, color=ORANGE))
axD.annotate("genomic A-tract\nending at the call", xy=(-3, 0.95), xytext=(9, 1.02),
             fontsize=5.2, color=PURPLE, ha="left", va="top", linespacing=1.3,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=PURPLE))

# (the three panel-d statistics sentences that stood in a text block below the
#  profile plot moved to the sidecar Legend, panel-d paragraph, verbatim in
#  substance: removed/kept counts and the P/R_det moves, the removed set's own
#  precision vs the kept set's, and the trigger breakdown)

# ===========================================================================
# panel e -- the molecule threshold is a trade surface; one point was gated
# ===========================================================================
axE = style(axbox(ML, 6.60, 2.48, 1.32, "e"))
panel_tag(axE, "e", "One point on a trade surface")
DS = [("pbmc", BLUE, "PBMC 10k v3"), ("mouse1", GREY, "testis mouse 1"), ("mouse2", GREY, "testis mouse 2")]
for ds, col, lab in DS:
    q = trade[trade.dataset == ds].sort_values("min_molecules")
    axE.plot(q.recall_detected_genes, q.atlas_agreement_precision, "-o", color=col, lw=0.9,
             ms=2.4, mfc="white", mew=0.9, alpha=1.0 if ds == "pbmc" else 0.55,
             zorder=4 if ds == "pbmc" else 3, label=lab)
    if ds == "pbmc":
        for _, r in q.iterrows():
            dx, dy = (-17, 1) if int(r.min_molecules) == 2 else (4, 3.5)
            axE.annotate(f"≥{int(r.min_molecules)}", (r.recall_detected_genes, r.atlas_agreement_precision),
                         textcoords="offset points", xytext=(dx, dy), fontsize=5.0, color=MUTED)
_d = trade[(trade.dataset == "pbmc") & (trade.min_molecules == 2)].iloc[0]
axE.plot(_d.recall_detected_genes, _d.atlas_agreement_precision, "o", ms=8.0, mfc="none",
         mec=GREEN, mew=1.3, zorder=5)
axE.annotate("the only pre-registered point", (_d.recall_detected_genes, _d.atlas_agreement_precision),
             textcoords="offset points", xytext=(8, -10), fontsize=5.4, color=GREEN,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=GREEN))
axE.set_xlabel("detected-gene recall R_det @100 bp", fontsize=5.8, color=MUTED, labelpad=0.5)
axE.set_ylabel("atlas-agreement P @100 bp", fontsize=5.8, color=MUTED, labelpad=1)
axE.set_xlim(0.055, 0.325)
axE.set_ylim(0.30, 1.03)     # top headroom shrank with the F1 note's move to the legend
axE.set_yticks([0.4, 0.6, 0.8, 1.0])
axE.grid(color=GRID, lw=0.5, alpha=0.7)
axE.legend(loc="lower left", frameon=False, fontsize=5.2, handlelength=1.2, borderpad=0.05,
           labelspacing=0.2)
# (the F1-honesty note — F1_det is higher at >=1 molecule than at the default,
#  a reliability choice, not the F1 optimum — moved to the sidecar Legend,
#  panel-e paragraph)

# ===========================================================================
# panel f -- what the default call set is used for
# ===========================================================================
FX, FW = ML + 2.88, CW - 2.88
axF = plain(axbox(FX, 6.60, FW, 1.32, "f"))
panel_tag(axF, "f", "What the paper does with the default call set")
LF = lineh(axF, 5.3)
CF = chars(FW * 0.94, 5.2)
STEPS = [
    (BLUE, "PAS × cell count matrix",
     [f"reads, not clips: the coverage peak this PAS claims, plus accepted read ends in",
      f"[site − {SEQLEN}, site + 25] outside it, split at midpoints between neighbouring PAS"]),
    (BLUE, "PAS usage per cell type",
     ["within-gene proportions, per-cell-type PAS counts, 3′UTR length"]),
    # (the "other five arms fail the pre-registered calibration rule" clause and
    #  the "yields are Fig 5's numbers" note moved to the sidecar Legend,
    #  panel-f paragraph)
    (GREEN, "calibrated switch test (14)",
     ["Fisher · --count-mode cells · marker pre-selection OFF",
      f"{FDR_NULL_P05:.1%} of label-shuffled tests p < 0.05; 0/{FDR_NULL_RUNS} null runs with a "
      f"q < 0.05 hit"]),
    (GREEN, "replication across biological units (20)",
     ["q < 0.05 in ≥ 2 patients, same direction, opposite-direction veto, |Δproportion| ≥ 0.1"]),
]
yy = 1.0
for step_i, (col, head, subs) in enumerate(STEPS):
    wrapped = [textwrap.fill(s, CF) for s in subs]
    nlines = sum(w.count("\n") + 1 for w in wrapped)
    h = LF * (1.15 + 1.05 * nlines)
    rbox(axF, 0.01, yy - h, 0.98, h, ec=col, lw=0.9)
    axF.text(0.026, yy - LF * 0.70, head, fontsize=5.8, color=col, fontweight="bold",
             ha="left", va="center")
    yline = yy - LF * 1.20
    for w in wrapped:
        axF.text(0.026, yline, w, fontsize=5.2, color=INK, ha="left", va="top", linespacing=1.3)
        yline -= LF * 1.05 * (w.count("\n") + 1)
    yy -= h
    if step_i < len(STEPS) - 1:          # arrows join boxes; none dangles after the last
        arrow(axF, 0.06, yy, 0.06, yy - LF * 0.45, ms=5.0)
        yy -= LF * 0.60

# ===========================================================================
# legend (sidecar only -- the on-figure footer was retired at the 2026-09-02
# submission pass; the sidecar '## Legend' is the single source of the caption)
# ===========================================================================
caption = (
    f"Figure 1 | What PeakATail does, drawn on the data it was measured on (code {CODE}; PBMC 10k v3 arm run "
    f"{RUN_STAMP}). "
    f"(a) A read is used when it is mapped, on the strand of the pass, carries a {P_RUN['cb_len']}-nt CB tag and "
    f"aligns over \u2264 {SEQLEN} nt; it enters the poly(A) channel when its 3\u2032-side terminal soft clip is "
    f"\u2265 {P_RUN['polya_min_clip']} nt and \u2265 {P_RUN['polya_min_purity']:.0%} A (+) / T (\u2212) with a "
    f"\u2265 {P_RUN['polya_min_clip']}-nt run flush to the alignment edge \u2014 {CLIP_RATE_GW:.3%} of accepted "
    f"CB reads genome-wide ({CLIP_READS_GW:,} / {CB_READS_GW:,}; 23 \u00a72, superseding the 1.152% of 10 "
    f"\u00a7R2). Support is distinct (CB, UB) molecules. "
    f"(b) Real BAM records and real GRCh38 sequence; the cluster, its read-weighted mode and its molecule count "
    f"are recomputed here from those records and reproduce column 5 of the caller\u2019s BED. Beneath it, at the "
    f"same scale, a real \u22652-molecule call the internal-priming filter removed. "
    f"(c) Sizes and precisions from one scorer (19 \u00a71/\u00a72). "
    f"(d) Profiles, trigger breakdown and the removed set\u2019s precision are computed by this script; the same "
    f"pipeline returns {float(ipprec.loc['kept_default', 'precision']):.4f} on the kept set, reproducing the "
    f"scorer\u2019s {S['default']['P']:.4f}. "
    f"(e) Fig 3\u2019s molecule sweep; only \u22652 was pre-registered. "
    f"(f) Only the FDR-controlling configuration is used, and only replicated switches are reported. "
    f"Caveats: atlas-agreement precision is agreement with a curated atlas, not ground truth; one donor, one "
    f"chemistry, and pipelines that trim poly(A) before alignment destroy this evidence (10 \u00a7R2); the "
    f"panel-b loci are examples chosen by stated criteria, not summaries; --polya-min-umis is 1 at the caller "
    f"and \u22652 molecules is the pre-registered OUTPUT. "
)

# ---------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------
for ext, kw in (("png", dict(dpi=600)), ("pdf", {})):
    p = FIGDIR / f"{NAME}.{ext}"
    fig.savefig(p, **kw)
    print("wrote", p)
fig.savefig(OUTDIR / f"{NAME}.png", dpi=600)

pd.DataFrame(rows_tsv).to_csv(OUTDIR / f"{NAME}.tsv", sep="\t", index=False)
pd.DataFrame(drawn_reads).to_csv(OUTDIR / f"{NAME}_reads.tsv", sep="\t", index=False)
seq_rows = []
for key, L in LOCI.items():
    for off in range(int(np.ceil(BX0)), int(np.floor(BX1)) + 1):
        i = off - L["seq_lo"]
        if 0 <= i < len(L["seq"]):
            seq_rows.append(dict(locus=key, chrom=L["chrom"], offset=off,
                                 pos_1based=L["site0"] + 1 + off, base=L["seq"][i]))
pd.DataFrame(seq_rows).to_csv(OUTDIR / f"{NAME}_sequence.tsv", sep="\t", index=False)
pd.DataFrame(dict(offset=OFF, A_fraction_kept=PROF_KEPT, A_fraction_removed=PROF_REM,
                  n_kept=N_PROF_KEPT, n_removed=N_PROF_REM)).to_csv(
    OUTDIR / f"{NAME}_ipprofile.tsv", sep="\t", index=False, float_format="%.6f")
trade.to_csv(OUTDIR / f"{NAME}_trade.tsv", sep="\t", index=False, float_format="%.6f")
print("wrote", OUTDIR / f"{NAME}.tsv", "and 4 companion TSVs")

# ---------------------------------------------------------------------------
# sidecar caption + index paragraph
# ---------------------------------------------------------------------------
cap_md = f"""# Fig 1 — `fig1_overview` caption (generated by `scripts/manuscript_figures/fig1_overview.py`)

## Legend

{caption}

### Panel-by-panel

**Panel a — ingestion.** The four acceptance conditions are `ema/countmatrix/read.py:read_check` as run
(`strategy={P_RUN['strategy']}`, `--seq-len {SEQLEN}`, `--barcode-tag {P_RUN['barcode_tag']}`, `cb_len
{P_RUN['cb_len']}`); the clip test is `ema/countmatrix/polya.py:clip_site` with `--polya-min-clip
{P_RUN['polya_min_clip']}` and `--polya-min-purity {P_RUN['polya_min_purity']}`, including the run-flush-to-the-boundary
condition that buys the 92× wrong-end specificity of the measured prototype (`10` §R3, which also bounds that
control: it constrains alignment artefacts, not genomic A-runs — which is what panel d is about). Support is distinct (CB, UB) molecules
(`--polya-min-umis {P_RUN['polya_min_umis']}` at the caller, with the ≥2-molecule output filtered downstream);
`--polya-clip-filter {P_RUN['polya_clip_filter']}`, so the −F 3844 subset is counted in parallel and written to
`pas_support.tsv` rather than used as the gate. The genome-wide clip rate **{CLIP_RATE_GW:.4%}**
({CLIP_READS_GW:,} / {CB_READS_GW:,} accepted CB reads) is `23_algorithm_roadmap.md` §2 / `results/algo_headroom/VERIFY`
§5 and **supersedes the 1.152% published in `10_caller_fix_plan.md`**, which was a head-sampling artefact of the
QC estimator. The same usable reads, extended to {SEQLEN} nt, are the coverage channel (tier 2, panel c) and the
count column the switch test uses (panel f): {WINDOW_READS:,} reads counted onto the panel-b PAS — the coverage
peak it claims plus read ends in [site − {SEQLEN}, site + 25], of which {WINDOW_ENDS_ONLY:,} — against
{len(KEPT['clipped'])} clip reads.

**Panel b — one real locus.** {KEPT['label']}, PAS id {KEPT['pas_id']} of the v2 IP arm. Context strip: coverage
of the {KEPT['n_usable']:,} accepted reads within ±300 bp, each extended to {SEQLEN} nt as the caller's coverage
path does; the local maximum is {int(COV.max()):,} reads at {COV_MAX_OFF:+d} bp while the cleavage site itself
carries {int(_cov_site)} — the reason a coverage summit is not a cleavage estimate (peak 3′ ends stop ~90–105 nt
short of it, `07` §3). Zoom: {len(show)} of the {len(KEPT['clipped'])} qualifying clip reads
({KEPT['molecules_recomputed']} distinct molecules) drawn to scale with their real soft-clipped sequence,
{len(un_show)} reads with no qualifying clip for contrast, the {len(clip_counts)} clip positions of the cluster with their read
counts, the cluster's own {SPAN_CLIP} bp extent (bracket; single linkage joins clip positions ≤
{P_RUN['polya_seed_window']} bp apart into one cluster whose read-weighted mode is the call — the parameter bounds the
GAP between consecutive positions, not the span, so the bracket is the cluster and not the parameter), the
canonical AATAAA at {HEX_OFF:+d}, and the internal-priming window (offsets
{IP_LO + 0.5:+.0f}…{IP_HI - 0.5:+.0f} from the called base). The lower track is {REM['label']} (PAS id {REM['pas_id']}), a ≥2-molecule tier-1 call
whose read-level evidence is indistinguishable from the call above and which the filter removes because a genomic
A-run of {REM_RUN[0]} nt ends at the call; its nearest strand-matched atlas representative site is
{REM['atlas_bp']:,} bp away (computed here). Every read is drawn from its own start to its OWN last aligned base,
so a drawn read whose clip sits away from the cluster's mode is drawn away from it: the reads and the
clip-position stems beneath them are the same records.
Gene context, from the annotation the run used (`data/references/gene_end.bed`): exactly one same-strand gene
body covers the kept call ({','.join(KEPT['genes_same_strand'])}), so its gene label is unambiguous on strand at
this locus — which is not a general claim, see caveat 4; the removed call also sits inside a same-strand gene
body ({','.join(REM['genes_same_strand'])}), which is the point, since being inside a gene is not evidence of a
PAS. **Neither call is labelled by its gene's 3′ end, because that label would be wrong**: measured here against
the Ensembl 99 GTF the STAR index was built from, the kept call is
{int(KEPT['tx3p']['distance'])} bp from the 3′ end of transcript {KEPT['tx3p']['transcript']} but
{abs(KEPT['gene_body_end'] - (KEPT['site0'] + 1)):,} bp from the {KEPT['gene']} **gene-body** end, and the removed
call is {int(REM['tx3p']['distance']):,} bp from the nearest {REM['gene']} transcript 3′ end — i.e. the kept
example is a transcript 3′ end and the removed example is not near one at all. Neither locus was picked by eye: both were drawn from a scripted filter over the whole call set, recorded in
this script and re-executed on every run. The kept example comes from calls with ≥10 molecules, no other call of
this arm within 300 bp, a strand-matched PolyASite 2.0 representative site within 25 bp, a canonical AATAAA in
−40…−5 and no ≥6-A run in the internal-priming window — **{SEL['kept_criteria_n']:,} of the
{SEL['kept_denominator']:,} default calls satisfy it**. The removed example comes from removals with ≥8
molecules, no strand-matched atlas site within 200 bp and a flagged A-run of ≥8 nt ending at −2…0 —
**{SEL['removed_criteria_n']:,} of the {SEL['removed_denominator']:,} removals satisfy it**. One survivor of each
filter was then taken for drawability ('+' strand so the drawing reads 5′→3′ left to right, and few enough reads
to draw). **They are examples, not summaries.** The molecule counts drawn
here were recomputed from the BAM by a from-scratch re-implementation of `read_check` + `clip_site` and are equal
to column 5 of the caller's own BED ({KEPT['molecules_bed']} and {REM['molecules_bed']}).

**Panel c — the funnel.** {NB['peaks_noip']:,} peak candidates → the internal-priming filter removes
{IP_FLAGGED_GW:,} ({IP_FLAGGED_GW / IP_PEAKS_GW:.2%}; `22` §6, the verified `--ip-rule kinnex` negative) →
tier 1 {NB['tier1_ip']:,} (clip-supported) and tier 2 {NB['tier2_ip']:,} (coverage summits no cluster claimed,
BED score 0) → ≥2 distinct clip molecules {NB['default']:,} = the pre-registered default. Scored precisions
(one scorer, `score_tool.py`, 100 bp point matching, strand-matched, curated PolyASite 2.0 representative sites):
tier 2 {S['tier2_ip']['P']:.4f} (n {S['tier2_ip']['n']:,}), tier 1 ≥1 molecule {S['tier1_ip']['P']:.4f}
(n {S['tier1_ip']['n']:,}), default {S['default']['P']:.4f} (n {S['default']['n']:,}), 3-seed gene-body-shuffled
null {S['default']['null_P']:.4f}. The pre-registered gate P@100 ≥ 0.50 is a PASS on all three arms —
{S['default']['P']:.4f} PBMC, 0.7450 / 0.7572 mice (`19` §1) — and the default sits
{S['default']['P'] / S['default']['null_P']:.0f}× above its shuffled null (grey tick in the in-node precision
bars). Tier 2 goes to its own file: it is not in the default and is not switch-tested here. The BED row counts
are 7–64 above the scored n because the scorer drops
non-primary contigs (`21` §10 item 1).

**Panel d — the internal-priming filter, measured here.** Mean genomic A-fraction by transcript-oriented offset
over all {N_PROF_KEPT:,} kept default calls and all {N_PROF_REM:,} calls the filter removed from the ≥2-molecule
arm. The kept curve carries the hexamer bump at ≈{HEX_OFF:+d}; the removed curve is the internal-priming
signature, and it lies **upstream of and at the call**, not downstream of it: the flagged A-run ends at or before
the called base in {TRIG['run_ends_at_or_before_call']:.1%} of removed calls, because the aligner runs into the
genomic A tract and clips where it ends, so the window's upstream arm (offsets −9…0 relative to the called
base; `--ip-window-left 10` counts from the BED end, which is the called base + 1) is what does the work even
though `ema/experimental/internal_priming.py` documents the rule as targeting a tract *downstream* of the
cleavage site. {TRIG['a6_rule'] / TRIG['n']:.2%} of removed calls are flagged by the
≥{P_RUN['ip_a_stretch']}-consecutive-A rule and only {TRIG['a_fraction_rule'] / TRIG['n']:.2%} by the
≥{P_RUN['ip_a_fraction']:.0%} A-fraction rule — the fraction rule is nearly inert as configured. The filter removes
{N_REMOVED:,} of {NB['ge2_noip']:,} ≥2-molecule calls ({N_REMOVED / NB['ge2_noip']:.1%}; BED rows — the scorer
sees {S['ge2_noip']['n'] - S['default']['n']:,} of {S['ge2_noip']['n']:,} on primary contigs, which is
19 §4's number), moving P@100
{S['ge2_noip']['P']:.4f} → {S['default']['P']:.4f} and R_det {S['ge2_noip']['R_det']:.4f} →
{S['default']['R_det']:.4f}; the removed calls' own atlas-agreement precision is **{P_REMOVED:.4f}** (computed
here with `bedtools window -w 100 -sm -u` against the same reference; the identical pipeline returns
{float(ipprec.loc['kept_default', 'precision']):.4f} on the kept set, reproducing the scorer's
{S['default']['P']:.4f}). Note that {1 - P_REMOVED:.1%} of what the filter removes has no atlas site within
100 bp but {P_REMOVED:.1%} does: the filter costs recall as well as buying precision. The non-IP ≥2-molecule
file is a diagnostic arm only and is never the default (`21` §7).

**Panel e — the trade surface.** The molecule sweep of Fig 3 (`results/figures/manuscript/fig3_tradeoff.tsv`,
same scorer and denominators): PBMC ≥1→≥10 molecules moves P@100 {_p1.atlas_agreement_precision:.3f}→
{trade[(trade.dataset == 'pbmc') & (trade.min_molecules == 10)].iloc[0].atlas_agreement_precision:.3f} while R_det
falls {_p1.recall_detected_genes:.3f}→
{trade[(trade.dataset == 'pbmc') & (trade.min_molecules == 10)].iloc[0].recall_detected_genes:.3f}. Only the ≥2
point was pre-registered ({PREREG}); the others are descriptive and **≥5 or ≥10 must never be presented as a
recommendation** (`21` §7 item 9). F1_det is higher at ≥1 molecule ({_p1.F1_detected_genes:.4f}) than at the
default ({_p.F1_detected_genes:.4f}) on PBMC, and likewise on both mice — the default is a reliability choice and
the figure says so rather than letting a reader infer that {_p.F1_detected_genes:.4f} is the tool's best.

**Panel f — downstream.** The count column is not the clip evidence. In `polya.py` `flush()` a tier-1 cluster
takes (i) the per-CB counts of every coverage candidate it suppresses (`_assign_candidate`) and (ii) the accepted
read ends in its cleavage window [site − {SEQLEN}, site + 25] (`--polya-count-window auto,25`) that lie outside
those candidates and inside its midpoint territory (`count_ends`, `_clip_to_neighbours`). At the panel-b site that
is {WINDOW_READS:,} reads against {len(KEPT['clipped'])} clip reads — and only {WINDOW_ENDS_ONLY:,} accepted read
ends fall in [site − {SEQLEN}, site + 25] at all, so **the window alone does not reproduce the count column**
({WINDOW_ENDS_ONLY / WINDOW_READS:.0%} of it) and must not be given as its definition. Of six differential-APA configurations only Fisher with `--count-mode cells`
and marker pre-selection off controls the false-discovery rate ({FDR_NULL_P05:.1%} of label-shuffled tests
p < 0.05, 0/{FDR_NULL_RUNS} null runs with any q < 0.05 hit, 0/{FDR_NULL_FAMILIES} BH families; `14`, SOUND);
switches are reported only when they replicate in ≥2 patients in the same direction with an opposite-direction
veto and |Δproportion| ≥ 0.1 (`20`). **No replication or switch count appears on this figure**: those numbers are
Fig 5's and are being regenerated on the Stage-3 v2 chain.

### Caveats that travel with this figure

1. Atlas-agreement precision is agreement with a curated atlas, not ground truth: atlas-novel true sites count as
   false positives, and the recall denominator is the atlas restricted to genes detected in the dataset.
2. The clip rate and every locus shown are CellRanger 3.0.0 / 10x 3′ v3 / 91 bp R2 on one PBMC donor. Any
   pipeline that trims poly(A) before alignment destroys this evidence channel entirely (`10` §R2), and the
   testis arm is STARsolo / 10x v2 / 98 bp R2.
3. The two loci in panel b are examples drawn from scripted filters that {SEL['kept_criteria_n']:,} and
   {SEL['removed_criteria_n']:,} calls respectively satisfy. They illustrate the mechanism; they are not evidence
   about the call set, which is what panels c–e are for.
4. The tier label reflects clip support, not 3′-end proximity or gene assignment; PAS→gene assignment in
   overlapping loci is under revision (`09` §5, issue #99), which is why no gene-level claim is made here. The
   panel-b labels name the nearest annotated **transcript** 3′ end and its distance, never "the gene's 3′ end":
   for the kept locus those differ by {abs(KEPT['gene_body_end'] - (KEPT['site0'] + 1)) - int(KEPT['tx3p']['distance']):,} bp.
5. `--polya-min-umis` is 1 at the caller; the ≥2-molecule default is the pre-registered *output*, filtered from
   the tier-1 file. Do not describe it as the caller's built-in default.
6. Panel d's profile and trigger breakdown are computed by this script, not quoted from a verified document; the
   method is stated on the figure and in this caption, and the same pipeline reproduces the scorer's precision on
   the kept set as a control.

## Provenance

- Script: `scripts/manuscript_figures/fig1_overview.py` (writes the figure, this sidecar and every audit TSV).
- Caller code {CODE} (frozen worktree `tools/pa-polya-run-9dfdefb3`); PBMC 10k v3 arm run {RUN_STAMP};
  pre-registration: {PREREG}.
- Every printed value, with its method and source: `results/figures/manuscript/{NAME}.tsv`, plus the
  `{NAME}_reads.tsv`, `{NAME}_sequence.tsv`, `{NAME}_ipprofile.tsv` and `{NAME}_trade.tsv` companions and the
  cached intermediates in `results/figures/manuscript/{NAME}_work/`.
- Sources of truth: `19_final_gate_v2.md` §1/§2 (verifier verdict FIXED), `22_performance_roadmap.md` §6,
  `23_algorithm_roadmap.md` §2 + `results/algo_headroom/VERIFY` §5, `14_switch_calibration_v2.md` (SOUND),
  `20_stage3_replication.md` (rule only), and the run's own `run_config.json` / `pas_support.tsv`.

### Index paragraph (for `05_figure_index.md` — add by hand, this script does not edit it)

**fig1_overview — Fig 1, what the method is, on real data ({RUN_STAMP[:10]}, built from verified artefacts
plus quantities computed in-script, each with its method recorded in `{NAME}.tsv`).** Six panels: (a) ingestion and the acceptance rules, with the corrected
genome-wide poly(A)-clip rate {CLIP_RATE_GW:.4%} ({CLIP_READS_GW:,} / {CB_READS_GW:,} accepted CB reads; `23` §2 —
this supersedes 1.152%); (b) one real PBMC locus drawn to scale — real reads, real soft clips, real GRCh38
sequence, the single-linkage cluster ({len(clip_counts)} clip positions spanning {SPAN_CLIP} bp, joined by the
{P_RUN['polya_seed_window']} bp gap rule) and its read-weighted mode, the hexamer at
{HEX_OFF:+d}, the internal-priming window — and beneath it a real ≥2-molecule call the filter removed;
(c) the funnel {NB['peaks_noip']:,} peaks → {IP_FLAGGED_GW:,} internal-priming removals → tier 1
{NB['tier1_ip']:,} / tier 2 {NB['tier2_ip']:,} → default {NB['default']:,}, with P@100
{S['tier2_ip']['P']:.4f} / {S['tier1_ip']['P']:.4f} / {S['default']['P']:.4f} against a
{S['default']['null_P']:.4f} shuffled null; (d) the filter measured — A-fraction profiles of kept vs removed
calls, {TRIG['a6_rule'] / TRIG['n']:.2%} of removals triggered by the ≥6-A run rule, the run ending at or before
the call in {TRIG['run_ends_at_or_before_call']:.0%}, and removed-call precision {P_REMOVED:.3f} vs kept
{S['default']['P']:.3f}; (e) the molecule trade surface with the single pre-registered point ringed and the
F1 honesty note; (f) the downstream chain, calibrated test and replication rule, with no Fig 5 counts.
**Caveats that travel with it:** atlas-agreement precision is not ground truth; single donor, single chemistry,
and poly(A)-trimming pipelines destroy the evidence (`10` §R2); the panel-b loci are examples, not summaries;
`--polya-min-umis` is 1 at the caller and ≥2 is the pre-registered *output*; panel d is computed in-script with
its method stated and a control that reproduces the scorer. Full caption:
`figures/fig1_overview.caption.md`; sources `19_final_gate_v2.md` §1/§2, `22_performance_roadmap.md` §6,
`23_algorithm_roadmap.md` §2, `14_switch_calibration_v2.md`, `20_stage3_replication.md`, and the run's own
`run_config.json` / `pas_support.tsv`.
"""
(FIGDIR / f"{NAME}.caption.md").write_text(cap_md)
print("wrote", FIGDIR / f"{NAME}.caption.md")

# ---------------------------------------------------------------------------
# edge check: nothing may be clipped; the outer 8 px of the PNG must be blank
# ---------------------------------------------------------------------------
try:
    from PIL import Image
    im = np.asarray(Image.open(FIGDIR / f"{NAME}.png").convert("L"))
    edge = 8
    border = np.concatenate([im[:edge, :].ravel(), im[-edge:, :].ravel(),
                             im[:, :edge].ravel(), im[:, -edge:].ravel()])
    ink = int((border < 250).sum())
    print(f"edge check: {im.shape[1]}x{im.shape[0]} px, ink pixels in the outer {edge} px = {ink}")
    assert ink == 0, "INK IN THE OUTER 8 PX -- something is clipped"
except ImportError:
    print("edge check skipped (no PIL)")

print(f"fig1_overview: default n={S['default']['n']:,} P@100={S['default']['P']:.4f} "
      f"R_det={S['default']['R_det']:.4f}; clip rate {CLIP_RATE_GW:.4%}; "
      f"IP removes {N_REMOVED:,} at P {P_REMOVED:.3f}")
