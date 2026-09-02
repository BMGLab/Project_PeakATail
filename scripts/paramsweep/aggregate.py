#!/usr/bin/env python3
"""scripts/paramsweep/aggregate.py -- score every arm x view and emit one tidy TSV.

Nothing re-implements a metric: every P/R/F1/null number is read out of
score_tool.py's own TSV (produced through scripts/prime/score_pas.sh, which is
symlinked into scripts/paramsweep/ so only the output roots differ), and the
Kinnex columns are read off that wrapper's own stdout.

Views per arm
  default  pas_tier1_ge2mol.bed   tier-1 AND IP-pass AND >=2 molecules (the v2 default arm)
  tier2    pas_tier2.bed          coverage-only calls, scored alone
  all      pas.bed                everything the arm emits
  matched  default, truncated to the baseline's call count by descending BED
           score (clip molecules), ties broken by genomic order -- the
           "an arm that just emits more calls looks better for free" control
  cand     every raw candidate the peak caller emitted (run/peakcalling/*.bed),
           before the IP veto, the cell filter and gene assignment -- the
           Lead-A ceiling set
"""
from __future__ import annotations
import os, re, subprocess, sys
from pathlib import Path
import pandas as pd

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
OUTROOT = Path("/mnt/ssd0/emaout/peakatail_benchmark/paramsweep")
SCORE_SH = WD / "scripts/paramsweep/score_pas.sh"
ENV = dict(os.environ, LC_ALL="C")
CONTIGS = {"human": "19,21", "mouse": "18,19"}


def run_score(bed: Path, label: str, species: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["bash", str(SCORE_SH), str(bed), label, "--species", species,
           "--slice-contigs", CONTIGS[species], "--outdir", str(outdir)]
    p = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    if p.returncode != 0:
        sys.stderr.write("SCORE FAILED %s: %s\n%s\n" % (label, p.stdout[-2000:], p.stderr[-2000:]))
        return None
    tsv = outdir / ("score_%s.tsv" % label)
    d = pd.read_csv(tsv, sep="\t")

    def g(panel, ref, cut, series="real"):
        r = d[(d.panel == panel) & (d.series == series) & (d.reference == ref)
              & (d.cutoff_bp == cut)]
        return float(r.value.iloc[0]) if len(r) else float("nan")

    meta = d[d.panel == "meta"]
    row = {
        "n_scored": int(meta.n_matched.iloc[0]) if len(meta) else -1,
        "n_raw": int(meta.n_query.iloc[0]) if len(meta) else -1,
        "P@10": g("precision", "atlas_full", 10),
        "P@25": g("precision", "atlas_full", 25),
        "P@50": g("precision", "atlas_full", 50),
        "P@100": g("precision", "atlas_full", 100),
        "R_det@100": g("recall", "atlas_detected", 100),
        "F1_det@100": g("f1", "atlas_detected", 100),
        "R_full@100": g("recall", "atlas_full", 100),
    }
    nulls = d[(d.panel == "precision") & (d.series == "null_genic")
              & (d.cutoff_bp == 100)].sort_values("replicate")
    row["null_P@100_mean"] = float(nulls.value.mean()) if len(nulls) else float("nan")
    for i, (_, r) in enumerate(nulls.iterrows(), 1):
        row["null_P@100_s%d" % i] = float(r.value)
    for key, pat in (("kin_t5@25", r"Kinnex t5\s+P@25\s+([0-9.]+)"),
                     ("kin_t20@25", r"Kinnex t20 P@25\s+([0-9.]+)"),
                     ("kin_decoy@25", r"Kinnex decoy@25\s+([0-9.]+)")):
        m = re.search(pat, p.stdout)
        row[key] = float(m.group(1)) if m else float("nan")
    return row


def top_n_by_score(src: Path, dst: Path, n: int):
    """Truncate a scored BED to its n highest-scoring calls (BED col 5 = clip
    molecules), ties broken by genomic order, then re-sort for bedtools."""
    rows = []
    with open(src) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            rows.append((float(p[4]), p))
    rows.sort(key=lambda t: (-t[0], t[1][0], int(t[1][1])))
    keep = [p for _, p in rows[:n]]
    keep.sort(key=lambda p: (p[0], int(p[1])))
    with open(dst, "w") as fh:
        for p in keep:
            fh.write("\t".join(p) + "\n")
    return len(keep)


def candidate_bed(run: Path, dst: Path):
    """Every raw candidate the peak caller emitted, as scorer-ready 3' points."""
    pk = run / "peakcalling"
    src = sorted(f for f in os.listdir(pk) if f.endswith((".pos.bed", ".neg.bed")))
    tmp = str(dst) + ".raw"
    n = 0
    with open(tmp, "w") as out:
        for f in src:
            with open(pk / f) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 6:
                        continue
                    s, e, st = int(p[1]), int(p[2]), p[5]
                    b = (e - 1, e) if st == "+" else (s, s + 1)
                    out.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (p[0], b[0], b[1], p[3], p[4], st))
                    n += 1
    subprocess.run("sort -k1,1 -k2,2n %s > %s && rm -f %s" % (tmp, dst, tmp),
                   shell=True, check=True, env=ENV)
    return n


def main():
    arms = [("pbmc_base", "human", "base", ""), ("m1_base", "mouse", "base", "")]
    seen = {"pbmc_base", "m1_base"}
    for g in ("GRID.tsv", "GRID_extra.tsv"):
        grid = WD / "scripts/paramsweep" / g
        if not grid.exists():
            continue
        with open(grid) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                label, species, block, flags = line.rstrip("\n").split("\t")
                if label in seen:
                    continue
                seen.add(label)
                arms.append((label, species, block, flags))

    only = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    base_n = {}
    out_rows = []
    for label, species, block, flags in arms:
        if only and label not in only:
            continue
        d = OUTROOT / label
        if not (d / "DONE.ok").exists():
            sys.stderr.write("[skip] %s not finished\n" % label)
            continue
        sc = d / "score"
        sc.mkdir(exist_ok=True)
        views = {"default": d / "pas_tier1_ge2mol.bed",
                 "tier2": d / "pas_tier2.bed",
                 "all": d / "pas.bed"}
        cand = sc / "candidates.bed"
        if not cand.exists():
            candidate_bed(d / "run", cand)
        views["cand"] = cand
        bn = base_n.get(species)
        if bn is None and label.endswith("_base"):
            with open(views["default"]) as fh:
                bn = base_n[species] = sum(1 for _ in fh)
        if bn is not None:
            with open(views["default"]) as fh:
                nn = sum(1 for _ in fh)
            if nn > bn:
                m = sc / "default_matched.bed"
                top_n_by_score(views["default"], m, bn)
                views["matched"] = m
        for view, bed in views.items():
            if not bed.exists() or bed.stat().st_size == 0:
                sys.stderr.write("[empty] %s %s\n" % (label, view))
                continue
            lab = "%s__%s" % (label, view)
            r = run_score(bed, lab, species, sc)
            if r is None:
                continue
            r.update(arm=label, species=species, block=block, view=view, flags=flags)
            out_rows.append(r)
            print("%-26s %-8s n=%-7d P@100=%.4f R_det=%.4f F1=%.4f kin_t5=%.4f"
                  % (label, view, r["n_scored"], r["P@100"], r["R_det@100"],
                     r["F1_det@100"], r["kin_t5@25"]), flush=True)

    df = pd.DataFrame(out_rows)
    cols = ["arm", "species", "block", "view", "n_scored", "n_raw",
            "P@10", "P@25", "P@50", "P@100", "R_det@100", "F1_det@100",
            "R_full@100", "null_P@100_mean", "null_P@100_s1", "null_P@100_s2",
            "null_P@100_s3", "kin_t5@25", "kin_t20@25", "kin_decoy@25", "flags"]
    df = df[[c for c in cols if c in df.columns]]
    dst = WD / "results/paramsweep/arms_all.tsv"
    if dst.exists() and only:
        old = pd.read_csv(dst, sep="\t")
        key = ["arm", "view"]
        old = old[~old.set_index(key).index.isin(df.set_index(key).index)]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(dst, sep="\t", index=False, float_format="%.6f")
    print("\nwrote %s (%d rows)" % (dst, len(df)))


if __name__ == "__main__":
    main()
