#!/usr/bin/env python3
"""
trusted_novel_pas.py -- the pre-registered "trusted de novo (atlas-novel) PAS" definition
(manuscript/13_reliability_positioning.md, section 3) and its atlas-independent validation.

Two sub-commands:

  call      apply the definition to a PAS BED, write trusted_novel.bed + a per-criterion
            funnel + one cleavage-point BED per funnel stage (inputs for `validate`)
  validate  fraction of a PAS set with a Kinnex long-read truth 3' end within 25 bp,
            against a width/chrom/strand-matched shuffled null (bedtools shuffle -chrom)

DEFINITION (13 section 3).  A called PAS is *trusted novel* iff ALL of:
  1. clip-supported        present in the clip-supported tier (ids of --clip-bed)
  2. support >= 2          BED column 5.  After Stage 1c that column is DISTINCT MOLECULES.
                           Older (stale) outputs carry READS there: pass --support-unit reads;
                           accepted with a loud warning -- the quantity is then NOT the
                           pre-registered one and the run is a mechanics dry-run only.
  3. not internal-priming  either the tool's annotatedpas.bed flag (--ip-from-annot) or
                           recomputed from the genome: >= 6 consecutive A or >= 70 % A in the
                           transcript-strand window +10..+30 nt downstream of the cleavage
                           site (same window as motif_validation.py panel c).
  4. canonical hexamer     AATAAA/ATTAAA ("strong") or any of the 12 canonical variants
                           ("weak" tier, a superset) lying FULLY inside -40..-5 nt upstream
                           of the cleavage site on the transcript strand.
  5. atlas-novel           >= 100 bp from every site of the curated atlas(es), strand-matched
                           (--atlas-ignore-strand for the stricter strand-agnostic variant).

trusted_novel.bed = sites passing 1-5 with the weak-or-strong hexamer tier (column 9 says
which); trusted_novel.strong.bed = the AATAAA/ATTAAA subset.

SEQUENCE WINDOWS.  One `bedtools getfasta -s -bedOut` call on a SYMMETRIC window
[c-W, c+W+1) around the cleavage base c.  Symmetric windows are the only layout for which
sequence index i maps to transcript-relative position r = i - W on BOTH strands; an
asymmetric window silently mis-registers the minus strand after getfasta's reverse
complement (the phantom -525 hump documented in scripts/manuscript_figures/motif_validation.py).
r < 0 is upstream (inside the transcript body), r > 0 is downstream of the poly(A) junction.

COORDINATES.  Cleavage base c (0-based) = end-1 on '+', start on '-'.  For the point BEDs the
pipeline writes (start = end-1) both rules give the same base.  Distances are bedtools
`closest -d` distances; for two 1-bp points that is exactly |c_query - c_ref|.

Only the standard library + bedtools (on PATH) are needed.  Every sort is done in Python on
ASCII strings (identical to LC_ALL=C order) and every subprocess runs with LC_ALL=C.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ENV = dict(os.environ, LC_ALL="C")

HEX_STRONG: Tuple[str, ...] = ("AATAAA", "ATTAAA")
# The 12 canonical poly(A) signal hexamers (Beaudoing 2000), in the order used by
# scripts/manuscript_figures/motif_validation.py.
HEX12: Tuple[str, ...] = ("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", "CATAAA",
                          "AATATA", "GATAAA", "AATGAA", "AAGAAA", "ACTAAA", "AATAGA")
TRUE_WORDS = {"true", "t", "1", "yes", "y"}
FALSE_WORDS = {"false", "f", "0", "no", "n"}
MD5_MAX_BYTES = 256 * 1024 * 1024


# ----------------------------------------------------------------------------------------
# basic I/O
# ----------------------------------------------------------------------------------------
@dataclass
class Site:
    uid: int            # row index in the input; the key used everywhere internally
    chrom: str
    start: int
    end: int
    name: str
    support: float
    strand: str

    @property
    def cleavage(self) -> int:
        """0-based cleavage base: 3'-most base of the feature on its own strand."""
        return self.end - 1 if self.strand == "+" else self.start


def warn(msg: str) -> None:
    sys.stderr.write(f"[trusted_novel_pas] WARNING: {msg}\n")
    sys.stderr.flush()


def info(msg: str) -> None:
    sys.stderr.write(f"[trusted_novel_pas] {msg}\n")
    sys.stderr.flush()


def read_bed6(path: Path) -> List[Site]:
    """Read a BED6(+) file into Site records.  Extra columns are ignored.

    Column 5 is parsed as a number (support); '.' becomes 0.  Strand must be + or -.
    """
    sites: List[Site] = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                raise ValueError(f"{path}:{ln}: expected >= 6 tab-separated columns, got {len(f)}")
            if f[5] not in ("+", "-"):
                raise ValueError(f"{path}:{ln}: strand must be '+' or '-', got {f[5]!r}")
            start, end = int(f[1]), int(f[2])
            if end <= start or start < 0:
                raise ValueError(f"{path}:{ln}: malformed interval {start}-{end}")
            sup = 0.0 if f[4] in (".", "") else float(f[4])
            sites.append(Site(len(sites), f[0], start, end, f[3], sup, f[5]))
    return sites


def read_fai(fasta: Path) -> Dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(
            f"{fai} not found; index the genome first (samtools faidx {fasta})")
    out: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def read_first_column(path: Path) -> List[str]:
    vals: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                vals.append(line.split("\t")[0].split()[0])
    return vals


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                f = line.split()
                out[f[0]] = int(f[1])
    return out


def md5_or_size(path: Path) -> str:
    p = Path(path)
    if not p.exists():
        return "missing"
    if p.stat().st_size > MD5_MAX_BYTES:
        return f"size={p.stat().st_size};mtime={int(p.stat().st_mtime)}"
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: Sequence[str], stdout_path: Optional[Path] = None, stderr_log: Optional[Path] = None) -> None:
    """Run a subprocess with LC_ALL=C; raise with the captured stderr on failure."""
    out_fh = open(stdout_path, "w") if stdout_path else subprocess.DEVNULL
    try:
        proc = subprocess.run(list(cmd), env=ENV, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    finally:
        if stdout_path:
            out_fh.close()
    if stderr_log and proc.stderr:
        with open(stderr_log, "a") as fh:
            fh.write(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")


def bedtools_version() -> str:
    try:
        return subprocess.run(["bedtools", "--version"], env=ENV, capture_output=True,
                              text=True).stdout.strip()
    except FileNotFoundError:
        raise SystemExit("bedtools not found on PATH")


def write_point_bed(sites: Iterable[Site], path: Path, name_of=lambda s: str(s.uid)) -> int:
    """Write sorted cleavage-point BED6 (chrom, c, c+1, name, support, strand)."""
    rows = sorted(((s.chrom, s.cleavage, s.cleavage + 1, name_of(s), s.support, s.strand)
                   for s in sites), key=lambda r: (r[0], r[1], r[2]))
    with open(path, "w") as fh:
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{fmt_num(r[4])}\t{r[5]}\n")
    return len(rows)


def fmt_num(x: float) -> str:
    return str(int(x)) if float(x).is_integer() else f"{x:g}"


def sorted_bed6_union(paths: Sequence[Path], out: Path, keep_chroms: Optional[set] = None) -> int:
    """Concatenate BED6+ files, keep the first 6 columns, sort (chrom, start, end), dedupe."""
    rows = set()
    for p in paths:
        with open(p) as fh:
            for ln, line in enumerate(fh, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    raise ValueError(f"{p}:{ln}: expected >= 6 columns")
                if f[5] not in ("+", "-"):
                    raise ValueError(f"{p}:{ln}: strand must be '+'/'-' (got {f[5]!r})")
                if keep_chroms is not None and f[0] not in keep_chroms:
                    continue
                rows.add((f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]))
    srt = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[5]))
    with open(out, "w") as fh:
        for r in srt:
            fh.write("\t".join(map(str, r)) + "\n")
    return len(srt)


# ----------------------------------------------------------------------------------------
# sequence criteria
# ----------------------------------------------------------------------------------------
def fetch_windows(sites: Sequence[Site], fasta: Path, contig_len: Dict[str, int], W: int,
                  workdir: Path, log: Optional[Path] = None) -> Dict[int, str]:
    """Strand-aware symmetric windows [c-W, c+W+1) via bedtools getfasta -s.

    Returns {uid: SEQ} (upper-case, transcript strand, length 2W+1); sites whose window
    would leave the contig, or whose contig is absent from the FASTA, are omitted
    (getfasta hard-errors on negative starts and silently skips over-end features, so they
    are filtered here and reported by the caller as 'sequence_testable = False').
    """
    bed = workdir / "windows.bed"
    n_written = 0
    with open(bed, "w") as fh:
        for s in sites:
            L = contig_len.get(s.chrom)
            if L is None:
                continue
            a, b = s.cleavage - W, s.cleavage + W + 1
            if a < 0 or b > L:
                continue
            fh.write(f"{s.chrom}\t{a}\t{b}\t{s.uid}\t0\t{s.strand}\n")
            n_written += 1
    seqs: Dict[int, str] = {}
    if n_written == 0:
        return seqs
    out = workdir / "windows.seq.tsv"
    run(["bedtools", "getfasta", "-fi", str(fasta), "-bed", str(bed), "-s", "-bedOut"],
        stdout_path=out, stderr_log=log)
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seq = f[-1].upper()
            if len(seq) != 2 * W + 1:
                raise RuntimeError(f"getfasta returned a {len(seq)}-nt window for uid {f[3]} "
                                   f"(expected {2 * W + 1}); symmetric-window invariant broken")
            seqs[int(f[3])] = seq
    return seqs


def rel_slice(seq: str, W: int, lo: int, hi: int) -> str:
    """Sub-sequence covering transcript-relative positions lo..hi inclusive (r = i - W)."""
    return seq[lo + W: hi + W + 1]


def hexamer_call(seq: str, W: int, lo: int = -40, hi: int = -5) -> Tuple[bool, bool, List[str]]:
    """(strong, any12, hits): hexamers lying FULLY inside positions lo..hi of the window."""
    sub = rel_slice(seq, W, lo, hi)
    hits = [h for h in HEX12 if h in sub]
    strong = any(h in HEX_STRONG for h in hits)
    return strong, bool(hits), hits


def ip_call(seq: str, W: int, lo: int = 10, hi: int = 30, a_stretch: int = 6,
            a_frac: float = 0.7) -> bool:
    """Internal-priming flag on the transcript-strand downstream window lo..hi (inclusive)."""
    d = rel_slice(seq, W, lo, hi)
    if not d:
        return False
    return ("A" * a_stretch in d) or (d.count("A") / len(d) >= a_frac)


def parse_flag(word: str) -> Optional[bool]:
    w = word.strip().lower()
    if w in TRUE_WORDS:
        return True
    if w in FALSE_WORDS:
        return False
    return None


def read_annot_flags(annot_bed: Path, column: int) -> Dict[str, Optional[bool]]:
    """{name (col 4): flag} from an annotatedpas.bed-style file; unparsable -> None."""
    flags: Dict[str, Optional[bool]] = {}
    with open(annot_bed) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < column:
                raise ValueError(f"{annot_bed}:{ln}: has {len(f)} columns, --ip-from-annot asks for {column}")
            flags[f[3]] = parse_flag(f[column - 1])
    return flags


# ----------------------------------------------------------------------------------------
# distances
# ----------------------------------------------------------------------------------------
def closest_distances(query_bed: Path, ref_bed: Path, workdir: Path, tag: str,
                      strand_aware: bool = True, log: Optional[Path] = None) -> Dict[str, Optional[int]]:
    """{query name: bedtools closest -d distance} (None when no reference on that chrom/strand).

    Both inputs must be sorted (chrom, start).  Ties resolve with -t first; the distance is
    the same for every tied record so the choice is irrelevant here.
    """
    out = workdir / f"closest.{tag}.tsv"
    cmd = ["bedtools", "closest", "-a", str(query_bed), "-b", str(ref_bed), "-d", "-t", "first"]
    if strand_aware:
        cmd.append("-s")
    run(cmd, stdout_path=out, stderr_log=log)
    dist: Dict[str, Optional[int]] = {}
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = int(f[-1])
            dist[f[3]] = None if (d < 0 or f[6] == ".") else d
    return dist


def wilson(k: int, n: int, z: float = 1.959964) -> Tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return centre - half, centre + half


# ----------------------------------------------------------------------------------------
# sub-command: call
# ----------------------------------------------------------------------------------------
def check_atlas_path(p: Path, allow_unsafe: bool) -> None:
    if "hg19" in p.name.lower() and not allow_unsafe:
        raise SystemExit(
            f"refusing atlas {p}: hg19 coordinates cannot be intersected with GRCh38 calls "
            f"(data/references/atlases/README.md); lift over first or pass --allow-unsafe-atlas")


def cmd_call(a: argparse.Namespace) -> int:
    t0 = time.time()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = Path(a.workdir) if a.workdir else outdir / ".work"
    workdir.mkdir(parents=True, exist_ok=True)
    log = workdir / "bedtools.stderr.log"
    if log.exists():
        log.unlink()
    stages_dir = outdir / "stages"
    stages_dir.mkdir(exist_ok=True)
    notes: List[str] = []

    if a.support_unit == "reads":
        msg = ("--support-unit reads: BED column 5 holds READS, not distinct molecules. The "
               "pre-registered criterion is >= 2 DISTINCT MOLECULES (Stage 1c). Treat this run "
               "as a mechanics dry-run; its numbers are NOT the pre-registered quantity.")
        warn(msg)
        notes.append(msg)

    fasta = Path(a.genome)
    contig_len = read_fai(fasta)
    sites = read_bed6(Path(a.pas_bed))
    info(f"read {len(sites)} PAS from {a.pas_bed}")

    # --- stage predicates ---------------------------------------------------------------
    keep_chroms: Optional[set] = None
    if a.restrict_chroms:
        keep_chroms = set(read_first_column(Path(a.restrict_chroms)))

    if a.clip_bed:
        clip_ids = set()
        with open(a.clip_bed) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    clip_ids.add(line.rstrip("\n").split("\t")[3])
        clip = {s.uid: (s.name in clip_ids) for s in sites}
        info(f"clip-supported ids: {len(clip_ids)} from {a.clip_bed}; "
             f"{sum(clip.values())} of {len(sites)} input PAS are clip-supported")
    else:
        msg = ("no --clip-bed given: assuming --pas-bed is ALREADY the clip-supported tier "
               "(every site passes criterion 1)")
        warn(msg)
        notes.append(msg)
        clip = {s.uid: True for s in sites}

    support_ok = {s.uid: (s.support >= a.min_support) for s in sites}

    W = max(abs(a.hex_lo), abs(a.hex_hi), abs(a.ip_lo), abs(a.ip_hi))
    seqs = fetch_windows(sites, fasta, contig_len, W, workdir, log)
    info(f"sequence windows (W={W}, symmetric, strand-aware): {len(seqs)} of {len(sites)} sites testable")

    annot_flags: Dict[str, Optional[bool]] = {}
    if a.ip_from_annot:
        ap, col = a.ip_from_annot.rsplit(":", 1)
        annot_flags = read_annot_flags(Path(ap), int(col))
        n_known = sum(v is not None for v in annot_flags.values())
        info(f"internal-priming flags from {ap} column {col}: {n_known} parsable of {len(annot_flags)}")

    ip_flag: Dict[int, Optional[bool]] = {}
    ip_source: Dict[int, str] = {}
    hex_strong: Dict[int, bool] = {}
    hex_any: Dict[int, bool] = {}
    hex_hits: Dict[int, List[str]] = {}
    for s in sites:
        seq = seqs.get(s.uid)
        flag = annot_flags.get(s.name) if annot_flags else None
        if flag is not None:
            ip_flag[s.uid], ip_source[s.uid] = flag, "annot"
        elif seq is not None:
            ip_flag[s.uid] = ip_call(seq, W, a.ip_lo, a.ip_hi, a.ip_a_stretch, a.ip_a_frac)
            ip_source[s.uid] = "recomputed"
        else:
            ip_flag[s.uid], ip_source[s.uid] = None, "untestable"
        if seq is not None:
            st, anyh, hits = hexamer_call(seq, W, a.hex_lo, a.hex_hi)
        else:
            st, anyh, hits = False, False, []
        hex_strong[s.uid], hex_any[s.uid], hex_hits[s.uid] = st, anyh, hits

    atlas_paths = [Path(p) for p in a.atlas]
    for p in atlas_paths:
        check_atlas_path(p, a.allow_unsafe_atlas)
    atlas_bed = workdir / "atlas.union.bed6"
    n_atlas = sorted_bed6_union(atlas_paths, atlas_bed)
    info(f"atlas union: {n_atlas} sites from {len(atlas_paths)} file(s)")
    q_bed = workdir / "query.points.bed"
    write_point_bed(sites, q_bed)
    dist = closest_distances(q_bed, atlas_bed, workdir, "atlas",
                             strand_aware=not a.atlas_ignore_strand, log=log)
    atlas_dist: Dict[int, Optional[int]] = {s.uid: dist.get(str(s.uid)) for s in sites}
    atlas_novel = {u: (d is None or d >= a.atlas_min_dist) for u, d in atlas_dist.items()}

    # --- the funnel (pre-registered order) ----------------------------------------------
    stages: List[Tuple[str, Dict[int, bool]]] = []
    if keep_chroms is not None:
        stages.append(("on_listed_contigs", {s.uid: s.chrom in keep_chroms for s in sites}))
    stages += [
        ("clip_supported", clip),
        (f"support_ge{fmt_num(a.min_support)}_{a.support_unit}", support_ok),
        ("sequence_testable", {s.uid: s.uid in seqs for s in sites}),
        ("not_internal_priming", {u: (v is False) for u, v in ip_flag.items()}),
        (f"hexamer_any12_{a.hex_lo}..{a.hex_hi}", hex_any),
        (f"atlas_novel_ge{a.atlas_min_dist}bp", atlas_novel),
    ]
    final_weak = "trusted_novel"
    final_strong = "trusted_novel_strong"

    alive = {s.uid for s in sites}
    fail_stage: Dict[int, str] = {}
    funnel_rows: List[dict] = []
    n_in = len(sites)
    funnel_rows.append({"stage": "input", "n_pass": n_in, "n_fail_at_stage": 0,
                        "frac_of_input": 1.0, "frac_of_previous": 1.0})
    stage_sets: List[Tuple[str, set]] = [("00_input", set(alive))]
    prev = n_in
    for k, (name, pred) in enumerate(stages, 1):
        nxt = {u for u in alive if pred.get(u, False)}
        for u in alive - nxt:
            fail_stage[u] = name
        funnel_rows.append({"stage": name, "n_pass": len(nxt), "n_fail_at_stage": len(alive) - len(nxt),
                            "frac_of_input": len(nxt) / n_in if n_in else float("nan"),
                            "frac_of_previous": len(nxt) / prev if prev else float("nan")})
        stage_sets.append((f"{k:02d}_{name}", set(nxt)))
        alive, prev = nxt, len(nxt)
    trusted = set(alive)
    funnel_rows.append({"stage": final_weak, "n_pass": len(trusted), "n_fail_at_stage": 0,
                        "frac_of_input": len(trusted) / n_in if n_in else float("nan"),
                        "frac_of_previous": 1.0})
    strong = {u for u in trusted if hex_strong[u]}
    for u in trusted - strong:
        fail_stage[u] = "hexamer_strong_only"
    funnel_rows.append({"stage": "hexamer_strong_AATAAA_ATTAAA", "n_pass": len(strong),
                        "n_fail_at_stage": len(trusted) - len(strong),
                        "frac_of_input": len(strong) / n_in if n_in else float("nan"),
                        "frac_of_previous": len(strong) / len(trusted) if trusted else float("nan")})
    funnel_rows.append({"stage": final_strong, "n_pass": len(strong), "n_fail_at_stage": 0,
                        "frac_of_input": len(strong) / n_in if n_in else float("nan"),
                        "frac_of_previous": 1.0})
    stage_sets.append((f"{len(stages) + 1:02d}_{final_weak}", trusted))
    stage_sets.append((f"{len(stages) + 2:02d}_{final_strong}", strong))

    # marginal (each criterion on its own, over the whole input)
    marginal = []
    for name, pred in stages:
        n_pass = sum(1 for s in sites if pred.get(s.uid, False))
        marginal.append({"criterion": name, "n_pass_of_input": n_pass, "n_fail_of_input": n_in - n_pass,
                         "frac_pass_of_input": n_pass / n_in if n_in else float("nan")})
    n_strong_all = sum(hex_strong.values())
    marginal.append({"criterion": "hexamer_strong_AATAAA_ATTAAA", "n_pass_of_input": n_strong_all,
                     "n_fail_of_input": n_in - n_strong_all,
                     "frac_pass_of_input": n_strong_all / n_in if n_in else float("nan")})
    n_untest = sum(1 for v in ip_source.values() if v == "untestable")
    marginal.append({"criterion": "untestable_sequence_window", "n_pass_of_input": n_in - n_untest,
                     "n_fail_of_input": n_untest,
                     "frac_pass_of_input": (n_in - n_untest) / n_in if n_in else float("nan")})

    # --- outputs --------------------------------------------------------------------------
    by_uid = {s.uid: s for s in sites}
    with open(outdir / "funnel.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["stage", "n_pass", "n_fail_at_stage", "frac_of_input",
                                           "frac_of_previous"], delimiter="\t")
        w.writeheader()
        for r in funnel_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(outdir / "funnel_marginal.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["criterion", "n_pass_of_input", "n_fail_of_input",
                                           "frac_pass_of_input"], delimiter="\t")
        w.writeheader()
        for r in marginal:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})

    with open(outdir / "site_annotations.tsv", "w") as fh:
        fh.write("\t".join(["uid", "name", "chrom", "start", "end", "strand", "cleavage0", "support",
                            "support_unit", "clip_supported", "support_pass", "sequence_testable",
                            "ip_flag", "ip_source", "hex_strong", "hex_any12", "hex_hits",
                            "atlas_dist_bp", "atlas_novel", "trusted_novel", "trusted_novel_strong",
                            "fail_stage"]) + "\n")
        for s in sites:
            fh.write("\t".join(map(str, [
                s.uid, s.name, s.chrom, s.start, s.end, s.strand, s.cleavage, fmt_num(s.support),
                a.support_unit, int(clip[s.uid]), int(support_ok[s.uid]), int(s.uid in seqs),
                "" if ip_flag[s.uid] is None else int(ip_flag[s.uid]), ip_source[s.uid],
                int(hex_strong[s.uid]), int(hex_any[s.uid]), ",".join(hex_hits[s.uid]) or ".",
                "" if atlas_dist[s.uid] is None else atlas_dist[s.uid], int(atlas_novel[s.uid]),
                int(s.uid in trusted), int(s.uid in strong), fail_stage.get(s.uid, "pass"),
            ])) + "\n")

    def write_trusted(uids: set, path: Path) -> None:
        rows = sorted((by_uid[u] for u in uids), key=lambda s: (s.chrom, s.cleavage, s.strand))
        with open(path, "w") as fh:
            fh.write("#chrom\tcleavage0\tcleavage0+1\tname\tsupport\tstrand\torig_start\torig_end\t"
                     "hex_tier\thex_hits\tip_flag\tatlas_dist_bp\n")
            for s in rows:
                fh.write("\t".join(map(str, [
                    s.chrom, s.cleavage, s.cleavage + 1, s.name, fmt_num(s.support), s.strand,
                    s.start, s.end, "strong" if hex_strong[s.uid] else "weak",
                    ",".join(hex_hits[s.uid]),
                    int(bool(ip_flag[s.uid])),
                    "NA" if atlas_dist[s.uid] is None else atlas_dist[s.uid]])) + "\n")

    write_trusted(trusted, outdir / "trusted_novel.bed")
    write_trusted(strong, outdir / "trusted_novel.strong.bed")
    for label, uids in stage_sets:
        write_point_bed((by_uid[u] for u in uids), stages_dir / f"{label}.bed", name_of=lambda s: s.name)

    manifest = {
        "tool": "trusted_novel_pas.py call",
        "argv": sys.argv,
        "started": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0)),
        "elapsed_s": round(time.time() - t0, 1),
        "bedtools": bedtools_version(),
        "python": sys.version.split()[0],
        "params": {
            "support_unit": a.support_unit, "min_support": a.min_support,
            "ip_window": [a.ip_lo, a.ip_hi], "ip_a_stretch": a.ip_a_stretch, "ip_a_frac": a.ip_a_frac,
            "ip_from_annot": a.ip_from_annot, "hex_window": [a.hex_lo, a.hex_hi],
            "hex_strong": list(HEX_STRONG), "hex12": list(HEX12),
            "atlas_min_dist": a.atlas_min_dist, "atlas_strand_aware": not a.atlas_ignore_strand,
            "symmetric_window_halfwidth_W": W,
        },
        "inputs": {
            "pas_bed": {"path": str(a.pas_bed), "md5": md5_or_size(Path(a.pas_bed))},
            "clip_bed": {"path": a.clip_bed, "md5": md5_or_size(Path(a.clip_bed))} if a.clip_bed else None,
            "genome": {"path": str(fasta), "md5": md5_or_size(fasta)},
            "atlas": [{"path": str(p), "md5": md5_or_size(p)} for p in atlas_paths],
            "restrict_chroms": a.restrict_chroms,
        },
        "funnel": funnel_rows,
        "n_trusted_novel": len(trusted),
        "n_trusted_novel_strong": len(strong),
        "warnings": notes,
        "stale_input_dry_run": a.support_unit == "reads",
    }
    with open(outdir / "run_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"# trusted_novel_pas call  (support unit = {a.support_unit}"
          f"{'  ** STALE / DRY-RUN: NOT the pre-registered quantity **' if a.support_unit == 'reads' else ''})")
    print("stage\tn_pass\tn_fail_at_stage\tfrac_of_input")
    for r in funnel_rows:
        print(f"{r['stage']}\t{r['n_pass']}\t{r['n_fail_at_stage']}\t{r['frac_of_input']:.4f}")
    print(f"# wrote {outdir / 'trusted_novel.bed'} ({len(trusted)} sites; strong subset {len(strong)})")
    if not a.keep_work:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


# ----------------------------------------------------------------------------------------
# sub-command: validate
# ----------------------------------------------------------------------------------------
def parse_labeled(items: Sequence[str]) -> List[Tuple[str, Path]]:
    out: List[Tuple[str, Path]] = []
    seen = set()
    for it in items:
        if "=" in it:
            label, p = it.split("=", 1)
        else:
            p = it
            label = Path(it).stem
        if label in seen:
            raise SystemExit(f"duplicate label {label!r}")
        seen.add(label)
        if not Path(p).exists():
            raise SystemExit(f"missing file: {p}")
        out.append((label, Path(p)))
    return out


def hit_fraction(dist: Dict[str, Optional[int]], names: Iterable[str], window: int) -> Tuple[int, int]:
    n = k = 0
    for nm in names:
        n += 1
        d = dist.get(nm)
        if d is not None and d <= window:
            k += 1
    return k, n


def cmd_validate(a: argparse.Namespace) -> int:
    t0 = time.time()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = Path(a.workdir) if a.workdir else outdir / ".work"
    workdir.mkdir(parents=True, exist_ok=True)
    log = workdir / "bedtools.stderr.log"
    if log.exists():
        log.unlink()
    sizes = read_chrom_sizes(Path(a.chrom_sizes))
    keep = set(sizes)
    queries = parse_labeled(a.query)
    truths = parse_labeled(a.truth)
    strand_aware = not a.ignore_strand
    seeds = [a.seed_base + i for i in range(a.seeds)]

    # truth sets restricted to the scorable contigs
    truth_beds: Dict[str, Path] = {}
    truth_n: Dict[str, int] = {}
    for label, p in truths:
        tb = workdir / f"truth.{label}.bed"
        truth_n[label] = sorted_bed6_union([p], tb, keep_chroms=keep)
        truth_beds[label] = tb
        info(f"truth {label}: {truth_n[label]} sites on listed contigs ({p})")

    rows: List[dict] = []
    seed_rows: List[dict] = []
    for qlabel, qpath in queries:
        sites = read_bed6(qpath)
        kept = [s for s in sites if s.chrom in keep]
        n_dropped = len(sites) - len(kept)
        if n_dropped:
            info(f"query {qlabel}: dropped {n_dropped} sites on contigs absent from {a.chrom_sizes}")
        # unique internal names so that duplicate input names cannot collide in closest output
        qbed = workdir / f"q.{qlabel}.bed"
        write_point_bed(kept, qbed)
        names = [str(s.uid) for s in kept]
        if not kept:
            warn(f"query {qlabel} is empty after contig restriction; skipping")
            continue
        # null sets: shuffle once per seed, reuse for every truth
        null_beds: Dict[int, Path] = {}
        for sd in seeds:
            raw = workdir / f"null.{qlabel}.s{sd}.raw.bed"
            srt = workdir / f"null.{qlabel}.s{sd}.bed"
            cmd = ["bedtools", "shuffle", "-i", str(qbed), "-g", str(a.chrom_sizes), "-chrom",
                   "-noOverlapping", "-seed", str(sd)]
            if a.incl:
                cmd += ["-incl", str(a.incl)]
            run(cmd, stdout_path=raw, stderr_log=log)
            sorted_bed6_union([raw], srt)
            null_beds[sd] = srt
        for tlabel, tbed in truth_beds.items():
            if truth_n[tlabel] == 0:
                warn(f"truth {tlabel} has no sites on listed contigs; skipping")
                continue
            d_obs = closest_distances(qbed, tbed, workdir, f"{qlabel}.{tlabel}.obs", strand_aware, log)
            k, n = hit_fraction(d_obs, names, a.window)
            frac = k / n
            lo, hi = wilson(k, n)
            null_fracs: List[float] = []
            for sd in seeds:
                d_null = closest_distances(null_beds[sd], tbed, workdir,
                                           f"{qlabel}.{tlabel}.null{sd}", strand_aware, log)
                kn, nn = hit_fraction(d_null, names, a.window)
                nf = kn / nn if nn else float("nan")
                null_fracs.append(nf)
                seed_rows.append({"query": qlabel, "truth": tlabel, "seed": sd, "n": nn, "n_hit": kn,
                                  "frac_hit": nf})
            m = sum(null_fracs) / len(null_fracs) if null_fracs else float("nan")
            sd_ = (math.sqrt(sum((x - m) ** 2 for x in null_fracs) / (len(null_fracs) - 1))
                   if len(null_fracs) > 1 else float("nan"))
            n_ge = sum(1 for x in null_fracs if x >= frac)
            rows.append({
                "query": qlabel, "truth": tlabel, "n_query": n, "n_truth": truth_n[tlabel],
                "window_bp": a.window, "strand_aware": int(strand_aware),
                "n_hit": k, "frac_hit": frac, "wilson95_lo": lo, "wilson95_hi": hi,
                "null_seeds": len(null_fracs), "null_mean": m, "null_sd": sd_,
                "null_min": min(null_fracs) if null_fracs else float("nan"),
                "null_max": max(null_fracs) if null_fracs else float("nan"),
                "enrichment": (frac / m) if m and m > 0 else float("inf") if frac > 0 else float("nan"),
                "z_vs_null": ((frac - m) / sd_) if sd_ and sd_ > 0 else float("nan"),
                "empirical_p": (1 + n_ge) / (len(null_fracs) + 1) if null_fracs else float("nan"),
                "meets_target_0.70": int(frac >= a.target),
            })
            info(f"{qlabel} vs {tlabel}: {k}/{n} = {frac:.4f} within {a.window} bp "
                 f"(null mean {m:.4f}, {len(null_fracs)} seeds)")

    cols = ["query", "truth", "n_query", "n_truth", "window_bp", "strand_aware", "n_hit", "frac_hit",
            "wilson95_lo", "wilson95_hi", "null_seeds", "null_mean", "null_sd", "null_min", "null_max",
            "enrichment", "z_vs_null", "empirical_p", "meets_target_0.70"]
    with open(outdir / "validation.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(outdir / "validation_null_seeds.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["query", "truth", "seed", "n", "n_hit", "frac_hit"], delimiter="\t")
        w.writeheader()
        for r in seed_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    manifest = {
        "tool": "trusted_novel_pas.py validate",
        "argv": sys.argv,
        "started": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0)),
        "elapsed_s": round(time.time() - t0, 1),
        "bedtools": bedtools_version(),
        "params": {"window_bp": a.window, "strand_aware": strand_aware, "seeds": seeds,
                   "incl": a.incl, "chrom_sizes": a.chrom_sizes, "target": a.target,
                   "null_model": "bedtools shuffle -chrom -noOverlapping (width, chrom and strand "
                                 "of every record preserved)" + (" restricted to -incl" if a.incl else "")},
        "queries": [{"label": l, "path": str(p), "md5": md5_or_size(p)} for l, p in queries],
        "truths": [{"label": l, "path": str(p), "md5": md5_or_size(p)} for l, p in truths],
        "stale_input_dry_run": bool(a.stale_dry_run),
    }
    with open(outdir / "run_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print("# trusted_novel_pas validate" + ("  ** STALE / DRY-RUN **" if a.stale_dry_run else ""))
    print("query\ttruth\tn_query\tn_hit\tfrac_hit\tnull_mean\tnull_sd\tenrichment\tempirical_p")
    for r in rows:
        print(f"{r['query']}\t{r['truth']}\t{r['n_query']}\t{r['n_hit']}\t{r['frac_hit']:.4f}\t"
              f"{r['null_mean']:.4f}\t{r['null_sd']:.4f}\t{r['enrichment']:.2f}\t{r['empirical_p']:.3f}")
    if not a.keep_work:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


# ----------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="trusted_novel_pas.py",
                                 description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("call", help="apply the trusted-novel definition; write BED + funnel")
    c.add_argument("--pas-bed", required=True, help="called PAS, BED6+ (col5 = support; 3' end = cleavage)")
    c.add_argument("--support-unit", required=True, choices=["molecules", "reads"],
                   help="what BED col 5 counts. 'reads' (pre-Stage-1c outputs) is accepted with a warning")
    c.add_argument("--min-support", type=float, default=2, help="criterion 2 threshold (default 2)")
    c.add_argument("--clip-bed", help="BED whose col-4 ids define the clip-supported tier (criterion 1); "
                                      "omit only if --pas-bed is already clip-only")
    c.add_argument("--genome", required=True, help="genome FASTA (.fai must exist) -- hexamer + IP windows")
    c.add_argument("--ip-from-annot", metavar="ANNOT_BED:COL",
                   help="take the internal-priming flag from column COL (1-based) of an annotatedpas.bed; "
                        "sites with an unparsable value fall back to recomputation")
    c.add_argument("--ip-lo", type=int, default=10, help="IP window start, transcript-relative (default +10)")
    c.add_argument("--ip-hi", type=int, default=30, help="IP window end, inclusive (default +30)")
    c.add_argument("--ip-a-stretch", type=int, default=6)
    c.add_argument("--ip-a-frac", type=float, default=0.7)
    c.add_argument("--hex-lo", type=int, default=-40, help="hexamer window start (default -40)")
    c.add_argument("--hex-hi", type=int, default=-5, help="hexamer window end, inclusive (default -5)")
    c.add_argument("--atlas", action="append", required=True,
                   help="curated atlas point BED6 (repeatable; union is used)")
    c.add_argument("--atlas-min-dist", type=int, default=100, help="criterion 5 distance (default 100)")
    c.add_argument("--atlas-ignore-strand", action="store_true",
                   help="measure atlas distance regardless of strand (default: strand-matched)")
    c.add_argument("--allow-unsafe-atlas", action="store_true", help="permit an atlas path containing 'hg19'")
    c.add_argument("--restrict-chroms", help="keep only contigs listed in col 1 of this file (e.g. chrom.sizes)")
    c.add_argument("--outdir", required=True)
    c.add_argument("--workdir", help="scratch dir (default <outdir>/.work)")
    c.add_argument("--keep-work", action="store_true")
    c.set_defaults(func=cmd_call)

    v = sub.add_parser("validate", help="Kinnex long-read concordance vs shuffled null")
    v.add_argument("--query", action="append", required=True, metavar="LABEL=BED",
                   help="PAS set(s) to score (repeatable); e.g. trusted_novel.bed or a stages/*.bed")
    v.add_argument("--truth", action="append", required=True, metavar="LABEL=BED",
                   help="Kinnex truth point BED6 (repeatable; e.g. t5=..., t20=..., t100=...)")
    v.add_argument("--chrom-sizes", required=True, help="two-column chrom sizes (shuffle space + contig filter)")
    v.add_argument("--incl", help="optional BED restricting where shuffled sites may land (e.g. gene bodies)")
    v.add_argument("--window", type=int, default=25, help="hit if truth 3' end within this many bp (default 25)")
    v.add_argument("--seeds", type=int, default=10, help="number of shuffled null replicates (default 10)")
    v.add_argument("--seed-base", type=int, default=1)
    v.add_argument("--ignore-strand", action="store_true", help="match truth regardless of strand")
    v.add_argument("--target", type=float, default=0.70, help="pre-registered target fraction (13 section 3)")
    v.add_argument("--stale-dry-run", action="store_true", help="mark the manifest/stdout as a dry run")
    v.add_argument("--outdir", required=True)
    v.add_argument("--workdir")
    v.add_argument("--keep-work", action="store_true")
    v.set_defaults(func=cmd_validate)
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    a = build_parser().parse_args(argv)
    if a.cmd == "call" and a.hex_lo > a.hex_hi:
        raise SystemExit("--hex-lo must be <= --hex-hi")
    if a.cmd == "call" and a.ip_lo > a.ip_hi:
        raise SystemExit("--ip-lo must be <= --ip-hi")
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
