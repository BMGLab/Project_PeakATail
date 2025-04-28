#!/bin/bash

inputbam=$1
outprefix=$2

samtools $inputbam

conda activate deeptools

bamCoverage -b $inputbam -o "${outfile}".bw --binSize 100 --normalizeUsing CPM


