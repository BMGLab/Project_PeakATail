#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
K=/mnt/ssd0/emaout/peakatail_benchmark/kinnex
G=/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/kinnex_truth/gemx
R=/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/kinnex_truth
rm -f $K/SCORE_GEMX.DONE.ok $K/SCORE_GEMX.FAILED.err
trap 'echo "FAILED rc=$?" > $K/SCORE_GEMX.FAILED.err' ERR
$K/score_vs_kinnex.sh $G/gemx_truth_t5.point.bed              gemxt5    $R/gemxt5
$K/score_vs_kinnex.sh $G/gemx_truth_t20.point.bed             gemxt20   $R/gemxt20
$K/score_vs_kinnex.sh $G/gemx_truth_t100.point.bed            gemxt100  $R/gemxt100
$K/score_vs_kinnex.sh $G/gemx_truth_t500.point.bed            gemxt500  $R/gemxt500
$K/score_vs_kinnex.sh $G/gemx_truth_t5_plus_decoy.point.bed   gemxt5ip  $R/gemxt5ip
touch $K/SCORE_GEMX.DONE.ok
