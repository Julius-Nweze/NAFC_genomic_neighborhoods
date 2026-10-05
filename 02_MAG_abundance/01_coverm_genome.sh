#!/bin/bash
# Title          : 01_coverm_genome.sh
# Description    : MAG relative abundance per sample with CoverM genome
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./01_coverm_genome.sh

set -euo pipefail

# Set the working directories
WKDIR="/path/to/your/Exp3_Metagenome"
TRIM_DIR="$WKDIR/trimmed_reads"
MAG_DIR="$WKDIR/MAGs"
OUTDIR="$WKDIR/CoverM"

mkdir -p "$OUTDIR"

# Paired read files of all samples (R1 R2 R1 R2 ...)
COUPLED=()
for r1 in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
    COUPLED+=("$r1" "${r1/_R1_trimmed/_R2_trimmed}")
done

# Relative abundance and TPM of every MAG in every sample
coverm genome \
    --coupled "${COUPLED[@]}" \
    --genome-fasta-directory "$MAG_DIR" \
    --genome-fasta-extension fa \
    -p minimap2-sr \
    -m relative_abundance tpm \
    --threads 10 \
    -o "$OUTDIR/MAG_relative_abundance.CoverM.tsv"
