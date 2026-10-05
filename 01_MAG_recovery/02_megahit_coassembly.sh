#!/bin/bash
# Title          : 02_megahit_coassembly.sh
# Description    : Co-assemble the trimmed reads of all samples with MEGAHIT
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 02_megahit_coassembly.sh

#SBATCH --account=<your-account>
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --job-name=megahit
#SBATCH --mem=500G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=megahit_%j.out
#SBATCH --error=megahit_%j.err

set -euo pipefail

# Load modules
module load StdEnv/2020 python/3.10 megahit/1.2.9

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
TRIM_DIR="$WKDIR/trimmed_reads"

# Comma-separated lists of the trimmed R1 and R2 files
R1s=$(ls "$TRIM_DIR"/*_R1_trimmed.fastq.gz | python3 -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
R2s=$(ls "$TRIM_DIR"/*_R2_trimmed.fastq.gz | python3 -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
echo "$R1s"
echo "$R2s"

# Co-assemble to produce a single contigs file
megahit -1 "$R1s" -2 "$R2s" -o "$WKDIR/ASSEMBLY" -t "$SLURM_CPUS_PER_TASK"

# Contigs used as input by 03_anvio_contigs_database.sh
mkdir -p "$WKDIR/contigs"
cp "$WKDIR/ASSEMBLY/final.contigs.fa" "$WKDIR/contigs/Contigs.fasta"

echo "Co-assembly completed successfully."
