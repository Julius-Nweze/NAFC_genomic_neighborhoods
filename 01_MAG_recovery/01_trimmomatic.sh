#!/bin/bash
# Title          : 01_trimmomatic.sh
# Description    : Quality-filter paired-end metagenomic reads with Trimmomatic
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 01_trimmomatic.sh

#SBATCH --account=<your-account>
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name=trimmomatic
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=trimmomatic_%j.out
#SBATCH --error=trimmomatic_%j.err

set -euo pipefail

# Load modules
module load StdEnv/2020 trimmomatic/0.39

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
RAW_DIR="$WKDIR/raw_reads"
TRIM_DIR="$WKDIR/trimmed_reads"

# Adapter folder supplied with the Trimmomatic module
adap="$EBROOTTRIMMOMATIC/adapters"

mkdir -p "$TRIM_DIR"

# Run Trimmomatic on the R1 and R2 files of each sample
for r1 in "$RAW_DIR"/*_R1.fastq.gz; do
    base=$(basename "$r1" _R1.fastq.gz)
    r2="$RAW_DIR/${base}_R2.fastq.gz"
    echo "Trimming: $base"

    java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" PE -threads "$SLURM_CPUS_PER_TASK" "$r1" "$r2" \
        "$TRIM_DIR/${base}_R1_trimmed.fastq.gz" "$TRIM_DIR/${base}_R1_unpaired.fastq.gz" \
        "$TRIM_DIR/${base}_R2_trimmed.fastq.gz" "$TRIM_DIR/${base}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$adap/NexteraPE-PE.fa":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Trimming completed successfully."
