#!/bin/bash
# Title          : 03_build_bacteria71_phylogenomic_tree.sh
# Description    : Concatenate Bacteria_71 marker proteins, align with FAMSA and build a FastTree phylogeny
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./03_build_bacteria71_phylogenomic_tree.sh

set -euo pipefail

BASE_DIR="/path/to/your/Ref_genomes/Tree"
ANALYSIS_DIR="$BASE_DIR/Analysis"
EXTERNAL_GENOMES_TSV="$ANALYSIS_DIR/external-genomes.tsv"
OUT_DIR="$ANALYSIS_DIR/phylogenomics_bacteria71"
ALIGNMENT_FA="$OUT_DIR/Bacteria_71_concatenated_aa.fa"
PARTITIONS_TXT="$OUT_DIR/Bacteria_71_partitions.nexus"
TREE_NWK="$OUT_DIR/Bacteria_71_fasttree.nwk"
WRAPPER_SCRIPT="$ANALYSIS_DIR/anvio_famsa_stdio_wrapper.sh"
WRAPPER_BIN_DIR="$OUT_DIR/bin"
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="anvio-8"
ALIGNER="${ALIGNER:-famsa}"
TREE_PROGRAM="${TREE_PROGRAM:-fasttree}"

mkdir -p "$OUT_DIR"

if [[ ! -f "$EXTERNAL_GENOMES_TSV" ]]; then
    echo "Missing external genomes file: $EXTERNAL_GENOMES_TSV" >&2
    exit 1
fi

if [[ ! -f "$CONDA_SH" ]]; then
    echo "Missing conda activation script: $CONDA_SH" >&2
    exit 1
fi

source "$CONDA_SH"
conda activate "$CONDA_ENV"

if [[ "$ALIGNER" == "famsa" ]]; then
    if [[ ! -x "$WRAPPER_SCRIPT" ]]; then
        echo "Missing executable FAMSA wrapper: $WRAPPER_SCRIPT" >&2
        exit 1
    fi

    mkdir -p "$WRAPPER_BIN_DIR"
    ln -sf "$WRAPPER_SCRIPT" "$WRAPPER_BIN_DIR/famsa"
    export PATH="$WRAPPER_BIN_DIR:$PATH"
fi

anvi-get-sequences-for-hmm-hits \
    -e "$EXTERNAL_GENOMES_TSV" \
    --hmm-sources Bacteria_71 \
    --get-aa-sequences \
    --concatenate-genes \
    --partition-file "$PARTITIONS_TXT" \
    --align-with "$ALIGNER" \
    --return-best-hit \
    --unique-genes \
    --just-do-it \
    -o "$ALIGNMENT_FA"

anvi-gen-phylogenomic-tree \
    -f "$ALIGNMENT_FA" \
    -o "$TREE_NWK" \
    --program "$TREE_PROGRAM"

echo
echo "Alignment:  $ALIGNMENT_FA"
echo "Partitions: $PARTITIONS_TXT"
echo "Tree:       $TREE_NWK"
