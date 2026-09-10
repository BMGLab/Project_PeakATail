#!/usr/bin/env python3
"""Collate Supplementary Table T7 from the verified run registry on disk.

T7 is the reproducibility record of the benchmark: for every competitor, on every
dataset, the exact version, the dependency versions, the invocation and the
run-log caveat. Nothing here is typed from memory -- versions come from each run's
own `tool_version.txt`, written by the run script at run time, and invocations from
the `run.sh` / driver actually executed.

  python3 scripts/benchmark_tools/collate_table_T7.py          # print markdown
  python3 scripts/benchmark_tools/collate_table_T7.py --tsv    # also write the TSV
"""
import re
import sys
from pathlib import Path

WD = Path("/mnt/ssd1/Projects/PeakATail_wd")
BT = WD / "results/benchmark_tools"
OUT_TSV = WD / "results/figures/manuscript/table_T7_competitor_runs.tsv"

DATASETS = [("pbmc_10k_v3", "PBMC 10k v3 (human, 10x 3' v3, CellRanger 3.0.0)"),
            ("gse104556", "GSE104556 testis (mouse, 10x 3' v2, STARsolo)")]

TOOLS = ["polyapipe", "scapture", "sierra", "scapatrap", "scutrquant"]
NAME = {"polyapipe": "polyApipe", "scapture": "SCAPTURE", "sierra": "Sierra",
        "scapatrap": "scAPAtrap", "scutrquant": "scUTRquant"}

# Parameters that decide the call set, read from each dataset's driver.
PARAMS = {
    ("polyapipe", "pbmc_10k_v3"): "`--cell_barcode_tag CB --umi_tag UB -t 16`",
    ("polyapipe", "gse104556"):   "`--cell_barcode_tag CB --umi_tag UB -t 16`",
    ("scapture", "pbmc_10k_v3"):  "`-m PAScall --species human -l 91 -p 16`; DeepPASS on CPU",
    ("scapture", "gse104556"):    "`-m PAScall --species mouse -l 98 -p 16`; annotation prebuilt from the run's own GTF",
    ("sierra", "pbmc_10k_v3"):    "`FindPeaks`/`CountPeaks`, `ncores = 16`, GRCh38.99 GTF",
    ("sierra", "gse104556"):      "`FindPeaks`/`CountPeaks`, `ncores = 16`, GRCm38.102 GTF",
    ("scapatrap", "pbmc_10k_v3"): "`readlength 91`; CellRanger filtered barcodes (11,769); `findUniqueMap` -> `dedupByPos` -> per-strand `findPeaksByStrand`",
    ("scapatrap", "gse104556"):   "`readlength 98`; STARsolo filtered barcodes; same four-stage chain",
    # the human run's own record states target and tech only; min_umis / strand are NOT
    # recorded for it, so they are not claimed here (they are, for the mouse run).
    ("scutrquant", "pbmc_10k_v3"): "`target utrome_hg38_v1`, `tech 10xv3`; min_umis and strand not recorded in this run's own version record",
    ("scutrquant", "gse104556"):   "`target utrome_mm10_v2`, `tech 10xv2`, `--fr-stranded`, `min_umis 500`",
}

# Run-log caveats, quoted from the verified compute record (manuscript/00_draft
# SUPPLEMENTARY Fig S8 legend) so the table and the figure cannot drift apart.
CAVEATS = {
    ("polyapipe", "pbmc_10k_v3"):  "Ran clean.",
    ("sierra", "pbmc_10k_v3"):     "Ran clean.",
    ("scapatrap", "pbmc_10k_v3"):  "Wall time 4:04:29 is a RESUMED run after an OOM-killed first attempt; a from-scratch run is 12:58:50 total machine time. Disclosed, never summed.",
    ("scutrquant", "pbmc_10k_v3"): "Wall time 30:45 SUMS two attempts (salvage rerun after exit 1).",
    ("scapture", "pbmc_10k_v3"):   "Wall time 12:14 sums its two mandatory stages; the annotation prebuild is untimed.",
    ("scapture", "gse104556"):     "Mouse-2 site-level run completed and was scored (P@100 0.672); only its per-cell quantification step failed. That run is not invalid.",
}
DEFAULT_CAVEAT = "No retry or skipped stage recorded in the run log."

NOT_RUN = [("scTail (Hou and Huang, 2025)", "both",
            "Not runnable on either BAM: the method reads the cleavage position from the "
            "cDNA end of read 1, and read 1 is 28 bp in these libraries. No number is "
            "reported for it anywhere.")]


def version_block(ds: str, tool: str) -> str:
    f = BT / ds / tool / "tool_version.txt"
    if not f.exists():
        return "**not recorded**"
    lines = [l.strip().replace("|", "/") for l in f.read_text().splitlines() if l.strip()]
    keep = [l for l in lines
            if re.search(r"\d+\.\d+|version|commit", l, re.I)
            and not l.startswith(("Usage", "-h/", "Module", "Description", "target:"))]
    return "; ".join(keep[:5]) if keep else lines[0]


def rows():
    for ds, ds_label in DATASETS:
        for tool in TOOLS:
            d = BT / ds / tool
            if not d.exists():
                continue
            yield {
                "tool": NAME[tool],
                "dataset": ds_label,
                "version": version_block(ds, tool),
                "params": PARAMS.get((tool, ds), "see run script"),
                "script": f"`results/benchmark_tools/{ds}/{tool}/run.sh`",
                "caveat": CAVEATS.get((tool, ds), DEFAULT_CAVEAT),
            }


def main():
    rs = list(rows())
    md = ["**Supplementary Table T7 | Competitor tool versions, parameters, run commands and "
          "run-log caveats.** Every version string is read from the `tool_version.txt` that "
          "each run wrote at run time; every invocation is the `run.sh` or driver actually "
          "executed. Both are in the analysis repository at the paths given. Caveats are "
          "carried verbatim from the verified compute record (Fig S8).", "",
          "| Tool | Dataset | Version and key dependencies | Parameters that set the call set | Run script | Run-log caveat |",
          "|---|---|---|---|---|---|"]
    for r in rs:
        md.append(f"| {r['tool']} | {r['dataset'].split(' (')[0]} | {r['version']} | "
                  f"{r['params']} | {r['script']} | {r['caveat']} |")
    for name, _, why in NOT_RUN:
        md.append(f"| {name} | both | not run | not applicable | not applicable | {why} |")
    md += ["", "All competitor runs used one scorer (`scripts/benchmark_tools/score_tool.py`) and "
           "the per-dataset detected-gene denominators given in Methods. PeakATail's own arms are "
           "recorded separately at the frozen commits named in Methods."]
    text = "\n".join(md)
    print(text)
    if "--tsv" in sys.argv:
        OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
        with open(OUT_TSV, "w") as fh:
            fh.write("tool\tdataset\tversion\tparameters\trun_script\tcaveat\n")
            for r in rs:
                fh.write("\t".join([r["tool"], r["dataset"], r["version"].replace("\t", " "),
                                    r["params"], r["script"], r["caveat"]]) + "\n")
        print(f"\nwrote {OUT_TSV}", file=sys.stderr)
    return text


if __name__ == "__main__":
    main()
