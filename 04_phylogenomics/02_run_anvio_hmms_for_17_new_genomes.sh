#!/bin/bash
# Title          : 02_run_anvio_hmms_for_17_new_genomes.sh
# Description    : Run anvi'o HMM searches (Bacteria_71) on newly added contigs databases
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/14
# Usage          : ./02_run_anvio_hmms_for_17_new_genomes.sh

set -euo pipefail

BASE_DIR="/path/to/your/Ref_genomes/Tree"
CONTIGS_DB_DIR="$BASE_DIR/Contigs_DB"
ANALYSIS_DIR="$BASE_DIR/Analysis"
LOG_DIR="$ANALYSIS_DIR/logs/anvi_run_hmms"
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="anvio-8"
THREADS="${THREADS:-4}"
HMM_SOURCE="${HMM_SOURCE:-Bacteria_71}"

mkdir -p "$ANALYSIS_DIR" "$LOG_DIR"

if [[ ! -f "$CONDA_SH" ]]; then
    echo "Missing conda activation script: $CONDA_SH" >&2
    exit 1
fi

source "$CONDA_SH"
conda activate "$CONDA_ENV"

mapfile -t DBS < <(
    find "$CONTIGS_DB_DIR" -maxdepth 1 -type f -name "*.db" -print0 |
        sort -z |
        while IFS= read -r -d '' db_path; do
            if ! sqlite3 "$db_path" \
                "select 1 from hmm_hits_info where source='$HMM_SOURCE' limit 1;" |
                grep -q 1; then
                basename "$db_path"
            fi
        done
)

if [[ "${#DBS[@]}" -eq 0 ]]; then
    echo "All contigs databases already have HMM source: $HMM_SOURCE"
    exit 0
fi

echo "Found ${#DBS[@]} contigs databases missing HMM source: $HMM_SOURCE"
echo

count=0

for db_name in "${DBS[@]}"; do
    db_path="$CONTIGS_DB_DIR/$db_name"
    log_path="$LOG_DIR/${db_name%.db}.log"

    if [[ ! -f "$db_path" ]]; then
        echo "Missing contigs database: $db_path" >&2
        exit 1
    fi

    count=$((count + 1))
    echo "[$count/${#DBS[@]}] Running anvi-run-hmms on $db_name"

    anvi-run-hmms \
        -c "$db_path" \
        -T "$THREADS" \
        --just-do-it \
        > "$log_path" 2>&1

    echo "[$count/${#DBS[@]}] Completed $db_name"
done

echo
echo "Processed ${#DBS[@]} contigs databases for HMM source: $HMM_SOURCE."
echo "Logs directory: $LOG_DIR"
