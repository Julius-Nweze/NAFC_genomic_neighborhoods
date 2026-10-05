#!/bin/bash
# Title          : 10_anvio_export_collection.sh
# Description    : Export the manually refined MAG collection
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 10_anvio_export_collection.sh

#SBATCH --account=<your-account>
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name=anvio_export_collection
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=export_collection_%j.out
#SBATCH --error=export_collection_%j.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 bowtie2/2.5.1 samtools/1.16

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
MAPPINGDIR="$WKDIR/Mapping"
DATABASE="$WKDIR/Database"
CONTIGS_DB="$DATABASE/Contigs.db"
MERGED_DIR="$WKDIR/Merged_profile"

# Check if merged profile exists
if [[ ! -f "$MERGED_DIR/PROFILE.db" ]]; then
    echo "ERROR: PROFILE.db not found in $MERGED_DIR" >&2
    exit 1
fi

# Run anvi-export-collection
anvi-export-collection \
    -p "$MERGED_DIR/PROFILE.db" \
    -C "MAGs_Metabat2" \
    -O "MAGs_Metabat2"

echo "Collection 'MAGs_Metabat2' exported successfully."
