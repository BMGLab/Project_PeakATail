# Define your input BAM
INPUT_BAM="input.bam"

samtools index $INPUT_BAM

# Get a list of chromosomes in the BAM file
CHRS=($(samtools idxstats "$INPUT_BAM" | cut -f1 | grep -v '*' ))

# Select one random chromosome
RANDOM_CHR=${CHRS[$RANDOM % ${#CHRS[@]}]}

# Output BAM filename
OUTPUT_BAM="${RANDOM_CHR}.bam"

# Extract reads mapped to the selected chromosome
samtools view -b "$INPUT_BAM" "$RANDOM_CHR" -o "$OUTPUT_BAM"

# Optional: index the output BAM
samtools index "$OUTPUT_BAM"

echo "Random chromosome selected: $RANDOM_CHR"
echo "Output written to: $OUTPUT_BAM"
