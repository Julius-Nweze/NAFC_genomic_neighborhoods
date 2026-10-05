#!/bin/bash
# Title          : 03_anvio_contigs_database.sh
# Description    : Reformat co-assembled contigs (>=1 kbp), build the anvi'o contigs database, run HMMs, COG and SCG-taxonomy annotation, and report assembly statistics
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 03_anvio_contigs_database.sh

#SBATCH --account=<your-account>
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=60G
#SBATCH --job-name=anvio_database
#SBATCH --output=anvioDB_%j.out
#SBATCH --error=anvioDB_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --export=ALL

# Stop on errors and undefined variables
set -euo pipefail

### ---- 1. Load required modules ----
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 diamond/0.9.36 gcc/9.3.0 blast+/2.10.1

### ---- 2. Activate anvi'o virtual environment ----
source /path/to/your/anvio/bin/activate

### ---- 3. Set working directories ----
WKDIR="/path/to/your/Exp3_Metagenome"
INPUT_FASTA="$WKDIR/contigs/Contigs.fasta"
OUTPUT_CONTIGS="$WKDIR/Contigs"
OUTPUT_DB="$WKDIR/Database"

# Create output directories
mkdir -p "$OUTPUT_CONTIGS" "$OUTPUT_DB"

# Check input file exists
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA not found: $INPUT_FASTA" >&2
    exit 1
fi

### ---- 4. Reformat FASTA (keep contigs >= 1,000 bp, simplify names) ----
anvi-script-reformat-fasta \
    "$INPUT_FASTA" \
    -o "$OUTPUT_CONTIGS/Contigs.fa" \
    -l 1000 \
    --simplify-names \
    --report "$OUTPUT_CONTIGS/Contigs_name_conversions.txt" \
    --seq-type NT

### ---- 5. Create contigs database (Prodigal gene calls) ----
anvi-gen-contigs-database \
    -f "$OUTPUT_CONTIGS/Contigs.fa" \
    -o "$OUTPUT_DB/Contigs.db"

### ---- 6. Run HMMs (single-copy core genes) ----
anvi-run-hmms \
    -c "$OUTPUT_DB/Contigs.db"

### ---- 7. Annotate open reading frames with NCBI COGs ----
# One-time setup: anvi-setup-ncbi-cogs --just-do-it
anvi-run-ncbi-cogs \
    -c "$OUTPUT_DB/Contigs.db"

### ---- 8. Populate the contigs database with SCG taxonomy ----
# One-time setup: anvi-setup-scg-taxonomy
anvi-run-scg-taxonomy \
    -c "$OUTPUT_DB/Contigs.db"

### ---- 9. Estimate taxonomy in metagenome mode ----
anvi-estimate-scg-taxonomy \
    -c "$OUTPUT_DB/Contigs.db" \
    --metagenome-mode

### ---- 10. Generate stats table ----
anvi-display-contigs-stats \
    "$OUTPUT_DB/Contigs.db" \
    --report-as-text \
    -o "$OUTPUT_DB/Contigs_contigs_stats.csv"

echo "Anvi'o contigs database workflow completed successfully."
