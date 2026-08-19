#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
rm -f $K/TRUTH2_GEMX.DONE.ok $K/TRUTH2_GEMX.FAILED.err
trap 'echo "FAILED rc=$?" > $K/TRUTH2_GEMX.FAILED.err' ERR
$K/build_truth2.sh gemx $K/scratch/gemx_bt /mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/kinnex_truth/gemx
touch $K/TRUTH2_GEMX.DONE.ok
