#!/bin/bash
# Title          : 07_anvio_merge.sh
# Description    : Merge the per-sample anvi'o profiles
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 07_anvio_merge.sh

#SBATCH --account=<your-account>
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=30
#SBATCH --job-name=anvio_merge
#SBATCH --mem=256G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=merge_%j.out
#SBATCH --error=merge_%j.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
MAPPINGDIR="$WKDIR/Mapping"
DATABASE="$WKDIR/Database"
CONTIGS_DB="$DATABASE/Contigs.db"
MERGED_DIR="$WKDIR/Merged_profile"

# Check if contigs database exists
if [[ ! -f "$CONTIGS_DB" ]]; then
    echo "ERROR: Contigs database not found: $CONTIGS_DB" >&2
    exit 1
fi

# Find all PROFILE.db files
PROFILE_DBS=("$MAPPINGDIR"/*_PROFILE/PROFILE.db)

if [[ ${#PROFILE_DBS[@]} -eq 0 ]]; then
    echo "No PROFILE.db files found in $MAPPINGDIR" >&2
    exit 1
fi

# Run anvi-merge
echo "Merging profiles..."
anvi-merge "${PROFILE_DBS[@]}" -o "$MERGED_DIR" -c "$CONTIGS_DB" --skip-hierarchical-clustering

echo "Profiles merged successfully into $MERGED_DIR"
