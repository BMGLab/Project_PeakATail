#!/usr/bin/env python3
"""Cross-check the caller's emitted pas_support.tsv features against the
OFFLINE feature table the scoring analysis actually used.

Two independent comparisons, both keyed on (contig, cleavage, strand) so the
two runs' different pas_id numbering never enters:

  A. vs ``results/algo_headroom/A3_scoring_model/data/candidates.tsv.gz`` --
     the stored table produced by a separate genome-wide run and a separate
     implementation (bedtools getfasta -s over a symmetric +/-120 nt window,
     A3's seq_features.py).  Sequence features are a pure function of the
     coordinate and the genome, so every site present in BOTH tables must
     agree exactly.  Context features are a function of the CANDIDATE SET,
     which differs between a 2-contig slice run and a genome-wide run, so
     they are reported separately and disagreement there is expected.

  B. vs a fresh run of A3's own ``seq_features.py`` over ``bedtools getfasta
     -s`` windows for a random sample of the slice run's sites -- an
     end-to-end recomputation that does not depend on the stored table.

Usage:
  check_features.py --run RUNDIR --a3 candidates.tsv.gz [--fasta FA]
                    [--contigs 19,21] [--sample N] [--out report.tsv]
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import subprocess
import sys
import tempfile

# tool column -> A3 column (identical semantics, sometimes a different name)
SEQ_MAP = {
    "seq_ok": "seq_ok",
    "ip_tool_flag": "ip_tool_flag",
    "ip_tool_afrac": "ip_tool_afrac",
    "ip_tool_arun": "ip_tool_arun",
    "a_count_d18": "a_count_p1_18",
    "a_frac_d30": "a_frac_d30",
    "a_run_d30": "a_run_d30",
    "kin_ip_flag": "kin_ip_flag",
    "hex_strong": "hex_strong",
    "hex_any12": "hex_any12",
    "hex_n_types": "hex_n_types",
    "hex_best_off": "hex_best_off",
    "hex_strong_off": "hex_strong_off",
}
CTX_MAP = {c: c for c in (
    "d_prev_cand", "d_next_cand", "n_cand_100", "n_cand_500",
    "mol_500_sum", "is_local_mol_max", "mol_frac_local")}
# v2 columns worth joining as a sanity check that the two runs agree at all
V2_MAP = {c: c for c in ("clip_reads", "clip_umis", "window_reads", "tier")}


def _num(x):
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return x


def load_run(run_dir: str, contigs: set[str] | None):
    """{(chrom, cleavage, strand): {col: value}} from pasbed.bed + pas_support.tsv."""
    sup = {}
    with open(os.path.join(run_dir, "pas_support.tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sup[f[0]] = dict(zip(cols, f))
    out = {}
    n_bed = 0
    for name in ("pasbed.bed",):
        path = os.path.join(run_dir, name)
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                chrom, start, end, pid, _score, strand = p[:6]
                if contigs and chrom not in contigs:
                    continue
                n_bed += 1
                row = sup.get(pid)
                if row is None:
                    continue
                cleav = int(end) - 1 if strand == "+" else int(start)
                out[(chrom, cleav, strand)] = row
    return out, n_bed, len(sup)


def load_a3(path: str, contigs: set[str] | None):
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = dict(zip(cols, f))
            if contigs and r["chrom"] not in contigs:
                continue
            out[(r["chrom"], int(r["cleavage"]), r["strand"])] = r
    return out


def compare(run, a3, mapping, label, report):
    keys = sorted(set(run) & set(a3))
    stats = {}
    for tcol, acol in mapping.items():
        n = agree = 0
        examples = []
        for k in keys:
            tv, av = run[k].get(tcol), a3[k].get(acol)
            if tv is None or av is None:
                continue
            n += 1
            if _num(tv) == _num(av):
                agree += 1
            elif len(examples) < 3:
                examples.append(f"{k}: tool={tv} offline={av}")
        stats[tcol] = (n, agree)
        rate = agree / n if n else float("nan")
        report.append(f"{label}\t{tcol}\t{acol}\t{n}\t{agree}\t{rate:.6f}\t"
                      + ("; ".join(examples) if examples else ""))
    return len(keys), stats


def sample_recompute(run, fasta, a3code, n_sample, seed, report):
    """B: re-derive the sequence features with A3's own script + bedtools."""
    keys = sorted(run)
    random.Random(seed).shuffle(keys)
    keys = keys[:n_sample]
    W = 120
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, "win.bed")
        with open(bed, "w") as fh:
            for i, (chrom, c, strand) in enumerate(keys):
                if c - W < 0:
                    continue
                fh.write(f"{chrom}\t{c - W}\t{c + W + 1}\t{i}\t0\t{strand}\n")
        seqtsv = os.path.join(td, "win.seq.tsv")
        with open(seqtsv, "w") as out:
            subprocess.run(["bedtools", "getfasta", "-s", "-bedOut",
                            "-fi", fasta, "-bed", bed],
                           check=True, stdout=out)
        feat = os.path.join(td, "seq_feat.tsv")
        subprocess.run([sys.executable, a3code, seqtsv, feat], check=True)
        ref = {}
        with open(feat) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                ref[int(f[0])] = dict(zip(cols, f))
    n_tot = 0
    for tcol, acol in SEQ_MAP.items():
        n = agree = 0
        examples = []
        for i, k in enumerate(keys):
            r = ref.get(i)
            if r is None:
                continue
            tv, av = run[k].get(tcol), r.get(acol)
            if tv is None or av is None:
                continue
            n += 1
            if _num(tv) == _num(av):
                agree += 1
            elif len(examples) < 3:
                examples.append(f"{k}: tool={tv} bedtools+A3={av}")
        n_tot = max(n_tot, n)
        rate = agree / n if n else float("nan")
        report.append(f"B_recompute\t{tcol}\t{acol}\t{n}\t{agree}\t{rate:.6f}\t"
                      + ("; ".join(examples) if examples else ""))
    return n_tot


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--a3", required=True)
    ap.add_argument("--fasta")
    ap.add_argument("--a3-code")
    ap.add_argument("--contigs", default="19,21")
    ap.add_argument("--sample", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260822)
    ap.add_argument("--out")
    a = ap.parse_args()

    contigs = set(a.contigs.split(",")) if a.contigs else None
    run, n_bed, n_sup = load_run(a.run, contigs)
    a3 = load_a3(a.a3, contigs)
    report = ["comparison\ttool_column\toffline_column\tn\tn_agree\tagreement\texamples"]
    n_common, _ = compare(run, a3, {**V2_MAP, **SEQ_MAP}, "A_stored", report)
    compare(run, a3, CTX_MAP, "A_stored_context", report)
    n_recomp = 0
    if a.fasta and a.a3_code:
        n_recomp = sample_recompute(run, a.fasta, a.a3_code, a.sample,
                                    a.seed, report)
    head = (f"# run={a.run}\n# a3={a.a3}\n# contigs={sorted(contigs or [])}\n"
            f"# run_pas_in_bed={n_bed} run_sidecar_rows={n_sup} "
            f"run_keyed={len(run)} offline_keyed={len(a3)} "
            f"common_coordinates={n_common} recomputed_sample={n_recomp}")
    text = head + "\n" + "\n".join(report) + "\n"
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
