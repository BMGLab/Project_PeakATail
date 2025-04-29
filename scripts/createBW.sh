#!/bin/bash

inputbam=$1
outprefix=$2

# Check if arguments are provided
if [ -z "$inputbam" ] || [ -z "$outprefix" ]; then
  echo "Usage: $0 <input.bam> <output_prefix>"
  exit 1
fi

# Check if BAM index exists
if [ ! -f "${inputbam}.bai" ]; then
  echo "BAM index not found. Creating index..."
  samtools index "$inputbam"
fi


# Example samtools command (optional, remove if not needed)
# Just to check the BAM file, you might want something like:
samtools quickcheck "$inputbam"
if [ $? -ne 0 ]; then
  echo "samtools quickcheck failed. BAM file might be corrupt."
  exit 1
fi

# Activate conda environment properly
source $(conda info --base)/etc/profile.d/conda.sh
conda activate deeptools

# Run bamCoverage
bamCoverage -b "$inputbam" -o "${outprefix}.bw" --binSize 100 --normalizeUsing CPM
