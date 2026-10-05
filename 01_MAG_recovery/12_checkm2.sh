#!/bin/bash
# Title          : 12_checkm2.sh
# Description    : Estimate MAG completeness and contamination with CheckM2
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 12_checkm2.sh

#SBATCH --account=<your-account>
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --job-name=checkm2
#SBATCH --output=checkm2_%j.out
#SBATCH --error=checkm2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --export=ALL

# Load conda and activate checkm2 environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate checkm2

# Define input and output directories
input="/path/to/your/Exp3_Metagenome/Merged_metabat2_summary/MAGs1"
output_pred="/path/to/your/Exp3_Metagenome/Merged_metabat2_summary/MAGs1/checkm2"

# Requires the CheckM2 database (one-time setup: checkm2 database --download)
mkdir -p "$output_pred"

# Run CheckM2 prediction
checkm2 predict -i "$input" -o "$output_pred" -x fasta --force
