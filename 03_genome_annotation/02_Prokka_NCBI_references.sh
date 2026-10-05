#!/bin/bash
# Title          : 02_Prokka_NCBI_references.sh
# Description    : Prokka annotation of the NCBI comparator genomes
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./02_Prokka_NCBI_references.sh

set -euo pipefail

# Activate the environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate prokka

# Working directory containing Reformatted/ (one FASTA file per comparator genome)
WKDIR="/path/to/your/Gene_analysis_mags/Ref"
cd "$WKDIR"

mkdir -p Prokka/Reformatted

# Prokka: rapid prokaryotic genome annotation
for f in Reformatted/*.fa; do
    name=$(basename "$f")
    echo "Running Prokka on $name"
    prokka --outdir "Prokka/Reformatted/$name" "$f" --cpus 2
done

echo "Prokka annotation of comparator genomes complete."
