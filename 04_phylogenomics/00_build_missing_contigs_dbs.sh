#!/bin/bash
# Title          : 00_build_missing_contigs_dbs.sh
# Description    : Build anvi'o contigs databases for genomes that do not yet have one
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/14
# Usage          : ./00_build_missing_contigs_dbs.sh

set -euo pipefail

BASE_DIR="/path/to/your/Ref_genomes"
GENOMES_DIR="$BASE_DIR/Genomes"
TREE_DIR="$BASE_DIR/Tree"
CONTIGS_DB_DIR="$TREE_DIR/Contigs_DB"
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="anvio-8"
THREADS="${THREADS:-4}"
MAX_JOBS="${MAX_JOBS:-4}"

mkdir -p "$CONTIGS_DB_DIR"

if [[ ! -f "$CONDA_SH" ]]; then
    echo "Missing conda activation script: $CONDA_SH" >&2
    exit 1
fi

normalize_name() {
    local base="$1"
    case "$base" in
        *.medaka|*.consensus)
            printf '%s\n' "$base"
            ;;
        *.*)
            printf '%s\n' "${base%.*}"
            ;;
        *)
            printf '%s\n' "$base"
            ;;
    esac
}

shopt -s nullglob
GENOME_FILES=("$GENOMES_DIR"/*.fasta)

if [[ ${#GENOME_FILES[@]} -eq 0 ]]; then
    echo "No genome FASTA files found in $GENOMES_DIR" >&2
    exit 1
fi

MISSING=()
for genome_path in "${GENOME_FILES[@]}"; do
    genome_base="$(basename "$genome_path" .fasta)"
    db_base="$(normalize_name "$genome_base")"
    db_path="$CONTIGS_DB_DIR/${db_base}.db"

    if [[ ! -f "$db_path" ]]; then
        MISSING+=("$genome_path")
    fi
done

if [[ ${#MISSING[@]} -eq 0 ]]; then
    echo "No missing contigs databases found."
    exit 0
fi

build_one() {
    local genome_path="$1"
    local genome_base
    local db_base
    local db_path
    local tmp_fasta

    genome_base="$(basename "$genome_path" .fasta)"
    db_base="$(normalize_name "$genome_base")"
    db_path="$CONTIGS_DB_DIR/${db_base}.db"
    tmp_fasta="$(mktemp --suffix=.fasta)"

    awk -v prefix="$db_base" '
        /^>/ {
            header = substr($0, 2)
            split(header, parts, /[[:space:]]+/)
            id = parts[1]
            gsub(/[^A-Za-z0-9_]/, "_", id)
            if (id !~ /[A-Za-z]/) {
                id = prefix "_" id
            }
            print ">" id
            next
        }
        { print }
    ' "$genome_path" > "$tmp_fasta"

    rm -f "$db_path"

    anvi-gen-contigs-database \
        -f "$tmp_fasta" \
        -o "$db_path" \
        -n "$db_base" \
        --num-threads "$THREADS"

    rm -f "$tmp_fasta"
}

source "$CONDA_SH"
conda activate "$CONDA_ENV"

count=0
running=0
for genome_path in "${MISSING[@]}"; do
    genome_base="$(basename "$genome_path" .fasta)"
    count=$((count + 1))
    echo "[$count/${#MISSING[@]}] Building contigs database for $genome_base"

    build_one "$genome_path" &
    running=$((running + 1))

    if (( running >= MAX_JOBS )); then
        wait -n
        running=$((running - 1))
    fi
done

while (( running > 0 )); do
    wait -n
    running=$((running - 1))
done

echo
echo "Generated ${#MISSING[@]} contigs databases in $CONTIGS_DB_DIR"
