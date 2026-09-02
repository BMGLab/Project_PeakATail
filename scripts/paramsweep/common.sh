#!/usr/bin/env bash
# scripts/paramsweep/common.sh -- the prime dev harness, with the output roots
# redirected to the paramsweep sandbox.  Nothing is re-implemented: this sources
# scripts/prime/common.sh and overrides exactly two variables, so run_slice.sh /
# score_pas.sh (symlinked into this directory) behave identically except for
# where they write.
source /mnt/ssd1/Projects/PeakATail_wd/scripts/prime/common.sh
OUTROOT=/mnt/ssd0/emaout/peakatail_benchmark/paramsweep
ARTROOT=$WD/results/paramsweep
