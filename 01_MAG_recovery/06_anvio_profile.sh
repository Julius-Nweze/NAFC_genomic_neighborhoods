#!/bin/bash
# Title          : 06_anvio_profile.sh
# Description    : Profile each sample's BAM file against the contigs database (one SLURM array task per sample)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch --array=1-118 06_anvio_profile.sh

#SBATCH --account=<your-account>
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name=anvio_profile
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=profile_%A_%a.out
#SBATCH --error=profile_%A_%a.err

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 bowtie2/2.5.1 samtools/1.16

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
MAPPINGDIR="$WKDIR/Mapping"
DATABASE="$WKDIR/Database"
CONTIGS_DB="$DATABASE/Contigs.db"

# Check if contigs database exists
if [[ ! -f "$CONTIGS_DB" ]]; then
    echo "ERROR: Contigs database not found: $CONTIGS_DB" >&2
    exit 1
fi

# Get BAM files into an array
BAM_FILES=("$MAPPINGDIR"/*.bam)

# Pick the BAM file corresponding to this array job
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
BAM="${BAM_FILES[$INDEX]}"

if [[ -z "$BAM" ]]; then
    echo "No BAM file found for task ID $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID processing BAM: $BAM"

# Run anvi-profile
anvi-profile -i "$BAM" -c "$CONTIGS_DB" -T "$SLURM_CPUS_PER_TASK" \
             --output-dir "$MAPPINGDIR/$(basename "$BAM" .bam)_PROFILE"

echo "Task $SLURM_ARRAY_TASK_ID finished successfully."
