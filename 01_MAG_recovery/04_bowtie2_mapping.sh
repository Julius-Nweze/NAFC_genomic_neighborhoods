#!/bin/bash
# Title          : 04_bowtie2_mapping.sh
# Description    : Map each sample's trimmed reads to the reformatted contigs with Bowtie2 and prepare BAM files for anvi'o
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : sbatch 04_bowtie2_mapping.sh (run after 03_anvio_contigs_database.sh)

#SBATCH --account=<your-account>
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=bowtie2_mapping
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --output=mapping_%j.out
#SBATCH --error=mapping_%j.err

set -euo pipefail

# Load modules
module load StdEnv/2020 python/3.10 scipy-stack/2023a prodigal/2.6.3 hmmer/3.2.1 bowtie2/2.5.1 samtools/1.16

# Activate anvi'o virtual environment
source /path/to/your/anvio/bin/activate

# Directories
WKDIR="/path/to/your/Exp3_Metagenome"
TRIM_DIR="$WKDIR/trimmed_reads"
CONTIGS="$WKDIR/Contigs/Contigs.fa"
MAPDIR="$WKDIR/Mapping"

# Check that the reformatted contigs exist
if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Contigs file not found: $CONTIGS" >&2
    exit 1
fi

mkdir -p "$MAPDIR"

# Build the Bowtie2 index for the contigs
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" "$CONTIGS" "$MAPDIR/contigs"

# Map each sample and convert to a sorted, indexed BAM file for anvi'o
for r1 in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
    sample=$(basename "$r1" _R1_trimmed.fastq.gz)
    r2="$TRIM_DIR/${sample}_R2_trimmed.fastq.gz"
    echo "Mapping: $sample"

    bowtie2 -x "$MAPDIR/contigs" -q -1 "$r1" -2 "$r2" --no-unal -p "$SLURM_CPUS_PER_TASK" -S "$MAPDIR/$sample.sam"
    samtools view -F 4 -bS "$MAPDIR/$sample.sam" > "$MAPDIR/$sample-RAW.bam"
    anvi-init-bam "$MAPDIR/$sample-RAW.bam" -o "$MAPDIR/$sample.bam"
    rm "$MAPDIR/$sample.sam" "$MAPDIR/$sample-RAW.bam"
done

# Mapping is done, and the Bowtie2 index is no longer needed
rm "$MAPDIR"/*.bt2

echo "Mapping completed successfully."
