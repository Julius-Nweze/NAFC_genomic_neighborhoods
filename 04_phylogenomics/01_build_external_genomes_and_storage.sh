#!/bin/bash
# Title          : 01_build_external_genomes_and_storage.sh
# Description    : Build the anvi'o external-genomes table and genomes storage
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./01_build_external_genomes_and_storage.sh

set -euo pipefail

BASE_DIR="/path/to/your/Ref_genomes/Tree"
CONTIGS_DB_DIR="$BASE_DIR/Contigs_DB"
ANALYSIS_DIR="$BASE_DIR/Analysis"
EXTERNAL_GENOMES_TSV="$ANALYSIS_DIR/external-genomes.tsv"
GENOMES_STORAGE_DB="$ANALYSIS_DIR/Tree-GENOMES.db"
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="anvio-8"

if [[ ! -d "$CONTIGS_DB_DIR" ]]; then
    echo "Missing contigs database directory: $CONTIGS_DB_DIR" >&2
    exit 1
fi

if [[ ! -f "$CONDA_SH" ]]; then
    echo "Missing conda activation script: $CONDA_SH" >&2
    exit 1
fi

mkdir -p "$ANALYSIS_DIR"

mapfile -t DBS < <(find "$CONTIGS_DB_DIR" -maxdepth 1 -type f -name '*.db' | sort)

if [[ ${#DBS[@]} -eq 0 ]]; then
    echo "No .db files found in $CONTIGS_DB_DIR" >&2
    exit 1
fi

TMP_TSV=$(mktemp)
TMP_NAMES=$(mktemp)
trap 'rm -f "$TMP_TSV" "$TMP_NAMES"' EXIT

printf 'name\tcontigs_db_path\n' > "$TMP_TSV"

for DB in "${DBS[@]}"; do
    RAW_NAME=$(basename "$DB" .db)
    SAFE_NAME=$(printf '%s\n' "$RAW_NAME" | sed 's/[^A-Za-z0-9_]/_/g')

    if [[ ! "$SAFE_NAME" =~ ^[A-Za-z] ]]; then
        SAFE_NAME="g_${SAFE_NAME}"
    fi

    printf '%s\n' "$SAFE_NAME" >> "$TMP_NAMES"
    printf '%s\t%s\n' "$SAFE_NAME" "$DB" >> "$TMP_TSV"
done

SANITIZED_DUPLICATES=$(sort "$TMP_NAMES" | uniq -d || true)

if [[ -n "$SANITIZED_DUPLICATES" ]]; then
    echo "Duplicate genome names detected after anvi'o-safe name conversion:" >&2
    printf '%s\n' "$SANITIZED_DUPLICATES" >&2
    exit 1
fi

mv "$TMP_TSV" "$EXTERNAL_GENOMES_TSV"
trap - EXIT
rm -f "$TMP_NAMES"

if [[ -f "$GENOMES_STORAGE_DB" ]]; then
    BACKUP_PATH="${GENOMES_STORAGE_DB}.bak.$(date +%Y%m%d_%H%M%S)"
    mv "$GENOMES_STORAGE_DB" "$BACKUP_PATH"
    echo "Backed up existing genomes storage to: $BACKUP_PATH"
fi

source "$CONDA_SH"
conda activate "$CONDA_ENV"

anvi-gen-genomes-storage \
    -e "$EXTERNAL_GENOMES_TSV" \
    -o "$GENOMES_STORAGE_DB"

echo
echo "External genomes file: $EXTERNAL_GENOMES_TSV"
echo "Genomes storage DB:   $GENOMES_STORAGE_DB"
echo "Genome count:         ${#DBS[@]}"
