#!/bin/bash
# Title          : 11_anvio_import_collection.sh
# Description    : Import the refined MAG collection with bin information
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 11_anvio_import_collection.sh

#SBATCH --account=<your-account>
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=anvio_import_collection
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=import_collection_%j.out
#SBATCH --error=import_collection_%j.err

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

# Files
COLLECTION_FILE="$WKDIR/MAGs_Metabat2.txt"
BINS_INFO_FILE="$WKDIR/MAGs_Metabat2-info.txt"

# Check files exist
if [[ ! -f "$MERGED_DIR/PROFILE.db" ]]; then
    echo "ERROR: PROFILE.db not found in $MERGED_DIR" >&2
    exit 1
fi

if [[ ! -f "$CONTIGS_DB" ]]; then
    echo "ERROR: Contigs.db not found in $DATABASE" >&2
    exit 1
fi

if [[ ! -f "$COLLECTION_FILE" ]]; then
    echo "ERROR: Collection file $COLLECTION_FILE not found" >&2
    exit 1
fi

if [[ ! -f "$BINS_INFO_FILE" ]]; then
    echo "ERROR: Bins info file $BINS_INFO_FILE not found" >&2
    exit 1
fi

# Run anvi-import-collection
anvi-import-collection \
    "$COLLECTION_FILE" \
    --bins-info "$BINS_INFO_FILE" \
    -p "$MERGED_DIR/PROFILE.db" \
    -c "$CONTIGS_DB" \
    -C "MAGs_Metabat2"

echo "Collection 'MAGs_Metabat2' imported successfully."
