#!/bin/bash
# Title          : 13_gtdbtk.sh
# Description    : Classify MAGs with GTDB-Tk (classify_wf)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 13_gtdbtk.sh

#SBATCH --account=<your-account>
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name=gtdbtkDB
#SBATCH --mem=40G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=gtdbtk_%j.out
#SBATCH --error=gtdbtk_%j.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Define input and output directories
INPUT=/path/to/your/Exp3_Metagenome/MAGs
OUTPUT=/path/to/your/Exp3_Metagenome/MAGs/GTDB

# Make sure output directory exists
mkdir -p $OUTPUT

# Run GTDB-Tk classification
gtdbtk classify_wf \
  --genome_dir $INPUT \
  --out_dir $OUTPUT \
  --extension fa \
  --cpus $SLURM_CPUS_PER_TASK \
  --skip_ani_screen
