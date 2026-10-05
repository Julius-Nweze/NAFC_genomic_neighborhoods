#!/bin/bash
# Title          : 14_MAG_phylogenomic_tree.sh
# Description    : MAG-only phylogenomic tree from 39 concatenated ribosomal single-copy genes (Bacteria_71)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 14_MAG_phylogenomic_tree.sh

#SBATCH --account=<your-account>
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name=phylogenomics
#SBATCH --mem=40G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=phylogenomics_%j.out
#SBATCH --error=phylogenomics_%j.err

set -euo pipefail

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
CONTIGS_DB="$WKDIR/Database/Contigs.db"
PROFILE_DB="$WKDIR/Merged_profile/PROFILE.db"
COLLECTION="MAGs_Metabat2"
OUTDIR="$WKDIR/Phylogenomics"

mkdir -p "$OUTDIR"

# Ribosomal proteins from Bacteria_71 used for the tree
GENES="Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_S7,Ribosomal_L9_C,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,ribosomal_L24,Ribosomal_L27,Ribosomal_L27A,Ribosomal_L28,Ribosomal_L29,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S8,Ribosomal_S9,Ribosomal_S10,Ribosomal_S11,Ribosom_S12_S23,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S19"

# Get the concatenated, aligned amino-acid sequences of these genes for every MAG in the collection
anvi-get-sequences-for-hmm-hits \
    -c "$CONTIGS_DB" \
    -p "$PROFILE_DB" \
    -C "$COLLECTION" \
    -o "$OUTDIR/MAGs-seqs-for-phylogenomics.fa" \
    --hmm-source Bacteria_71 \
    --gene-names "$GENES" \
    --return-best-hit \
    --get-aa-sequence \
    --concatenate-genes

# Compute the phylogenomic tree
anvi-gen-phylogenomic-tree \
    -f "$OUTDIR/MAGs-seqs-for-phylogenomics.fa" \
    -o "$OUTDIR/MAGs-phylogenomic-tree.txt"

echo "Phylogenomic tree completed successfully."
