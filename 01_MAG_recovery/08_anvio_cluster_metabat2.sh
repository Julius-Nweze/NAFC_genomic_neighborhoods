#!/bin/bash
# Title          : 08_anvio_cluster_metabat2.sh
# Description    : Bin contigs with MetaBAT2 through anvi'o
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 08_anvio_cluster_metabat2.sh

#SBATCH --account=<your-account>
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=30
#SBATCH --job-name=anvio_cluster
#SBATCH --mem=126G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=cluster_%j.out
#SBATCH --error=cluster_%j.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 gcc/9.3.0 metabat/2.14 

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
MERGED_DIR="$WKDIR/Merged_profile"
DATABASE="$WKDIR/Database"
CONTIGS_DB="$DATABASE/Contigs.db"
PROFILE_DB="$MERGED_DIR/PROFILE.db"

# Check for required files
if [[ ! -f "$PROFILE_DB" ]]; then
    echo "ERROR: PROFILE.db not found: $PROFILE_DB" >&2
    exit 1
fi

if [[ ! -f "$CONTIGS_DB" ]]; then
    echo "ERROR: Contigs database not found: $CONTIGS_DB" >&2
    exit 1
fi

# Run anvi-cluster-contigs
echo "Clustering contigs with MetaBAT2..."

anvi-cluster-contigs \
    -p "$PROFILE_DB" \
    -c "$CONTIGS_DB" \
    -C Bins1 \
    --driver metabat2 \
    --just-do-it

echo "Clustering completed."
