#!/bin/bash
# Title          : 05_metabat2_binning.sh
# Description    : Bin the reformatted contigs with MetaBAT2
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 05_metabat2_binning.sh

#SBATCH --account=<your-account>
#SBATCH --time=00-05:00:00
#SBATCH --cpus-per-task=20
#SBATCH --job-name=binning
#SBATCH --mem=40G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=metabat2_%j.out
#SBATCH --error=metabat2_%j.err

set -euo pipefail

### ---- 1. Load required modules ----
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 gcc/9.3.0 metabat/2.14

### ---- 2. Activate anvi'o environment ----
source /path/to/your/anvio/bin/activate

### ---- 3. Set working directories ----
WKDIR="/path/to/your/Exp3_Metagenome"
OUTDIR="$WKDIR/Metabat2"
CONTIGS="$WKDIR/Contigs/Contigs.fa"

mkdir -p "$OUTDIR"

### ---- 4. Check for required input ----
if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Contigs file not found at $CONTIGS" >&2
    exit 1
fi

### ---- 5. Run MetaBAT2 ----
echo "Starting MetaBAT2 binning..."
metabat2 -i "$CONTIGS" -o "$OUTDIR/Exp3_metabat2" -t "$SLURM_CPUS_PER_TASK"

echo "MetaBAT2 binning completed successfully."
