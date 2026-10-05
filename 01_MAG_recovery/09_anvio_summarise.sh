#!/bin/bash
# Title          : 09_anvio_summarise.sh
# Description    : Summarise the MetaBAT2 bins
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 09_anvio_summarise.sh

#SBATCH --account=<your-account>
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=30
#SBATCH --job-name=anvio_summarize
#SBATCH --mem=126G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=summarize_%j.out
#SBATCH --error=summarize_%j.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
DATABASE="$WKDIR/Database"
MERGED_DIR="$WKDIR/Merged_profile"
PROFILE_DB="$MERGED_DIR/PROFILE.db"
OUTPUT_DIR="$WKDIR/Merged_metabat2_summary"
CONTOUR="Bins1"

# Check if merged profile exists
if [[ ! -f "$PROFILE_DB" ]]; then
    echo "ERROR: Merged profile database not found: $PROFILE_DB" >&2
    exit 1
fi

# Run anvi-summarize
anvi-summarize -p "$PROFILE_DB" -c "$DATABASE/Contigs.db" -o "$OUTPUT_DIR" -C "$CONTOUR"

echo "Anvi'o summary completed successfully."
