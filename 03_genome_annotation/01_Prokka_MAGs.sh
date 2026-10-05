#!/bin/bash
# Title          : 01_Prokka_MAGs.sh
# Description    : Prokka annotation of the plant-root MAGs
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./01_Prokka_MAGs.sh

set -euo pipefail

# Activate the environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate prokka

# Working directory containing MAGs/ (one FASTA file per MAG)
WKDIR="/path/to/your/Gene_analysis_mags/MAGs"
cd "$WKDIR"

mkdir -p Prokka/MAGs

# Prokka: rapid prokaryotic genome annotation
for f in MAGs/*.fa.fasta; do
    name=$(basename "$f")
    echo "Running Prokka on $name"
    prokka --outdir "Prokka/MAGs/$name" "$f" --cpus 10
done

echo "Prokka annotation of MAGs complete."
