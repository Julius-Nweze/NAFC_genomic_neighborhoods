#!/bin/bash
# Title          : 04_Prokka_OSPW_GROW_genomes.sh
# Description    : Prokka annotation of the GROW OSPW-associated genomes and export of CDS coordinates to TSV
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./04_Prokka_OSPW_GROW_genomes.sh

set -euo pipefail

# Activate the environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate prokka

# Folder containing the GROW_*_consensus.fasta genomes
BASE_DIR="/path/to/your/GROW_genomes"
CPUS="${CPUS:-11}"

shopt -s nullglob
GENOMES=("$BASE_DIR"/GROW_*_consensus.fasta)

if [ "${#GENOMES[@]}" -eq 0 ]; then
    echo "No GROW_*_consensus.fasta files found in $BASE_DIR" >&2
    exit 1
fi

for FILE in "${GENOMES[@]}"; do
    BASENAME="$(basename "$FILE" .fasta)"
    OUTPUT_DIR="$BASE_DIR/${BASENAME}_Prokka"
    GFF_FILE="$OUTPUT_DIR/${BASENAME}.gff"
    OUTPUT_TSV="$OUTPUT_DIR/${BASENAME}_cds.tsv"

    if [ -f "$GFF_FILE" ] && [ -f "$OUTPUT_TSV" ]; then
        echo "Skipping $BASENAME because $OUTPUT_DIR already exists."
        continue
    fi

    echo "Running Prokka on $BASENAME"
    prokka \
        --force \
        --outdir "$OUTPUT_DIR" \
        --prefix "$BASENAME" \
        --cpus "$CPUS" \
        "$FILE"

    # Output gene information and positions from the GFF file
    if [ -f "$GFF_FILE" ]; then
        {
            printf "Seqid\tStart\tEnd\tStrand\tID\teC_number\tName\tdb_xref\tgene\tinference\n"
            awk '
            BEGIN {
                FS = "\t";
                OFS = "\t";
            }
            $3 == "CDS" {
                split($9, attributes, ";");
                id = ""; eC_number = ""; name = ""; db_xref = ""; gene = ""; inference = "";
                for (i in attributes) {
                    split(attributes[i], kv, "=");
                    if (kv[1] == "ID") id = kv[2];
                    if (kv[1] == "eC_number") eC_number = kv[2];
                    if (kv[1] == "Name") name = kv[2];
                    if (kv[1] == "db_xref") db_xref = kv[2];
                    if (kv[1] == "gene") gene = kv[2];
                    if (kv[1] == "inference") inference = kv[2];
                }
                print $1, $4, $5, $7, id, eC_number, name, db_xref, gene, inference
            }' "$GFF_FILE"
        } > "$OUTPUT_TSV"
    fi
done

echo "Prokka annotation of GROW genomes complete."
