#!/bin/bash
# Title          : 03_Prokka_OSPW_genomes.sh
# Description    : Prokka annotation of the OSPW-associated genomes (parallel, COLLA locus tags)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./03_Prokka_OSPW_genomes.sh <dataset_root> [parallel_jobs=3] [cpus_per_job=4]

set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <dataset_root> [parallel_jobs=3] [cpus_per_job=4]" >&2
  exit 2
fi

dataset_root=$1
parallel_jobs=${2:-3}
cpus_per_job=${3:-4}
conda_env=${PROKKA_CONDA_ENV:-prokka}

if [ ! -d "$dataset_root" ]; then
  echo "Dataset root not found: $dataset_root" >&2
  exit 1
fi

find "$dataset_root" -mindepth 1 -maxdepth 1 -type d \
  ! -name metadata \
  ! -name NA_cluster_review \
  ! -name MAGs_cluster_review \
  -print0 |
  xargs -0 -n 1 -P "$parallel_jobs" bash -c '
    set -euo pipefail

    genome_dir=$1
    conda_env=$2
    cpus_per_job=$3
    sample=$(basename "$genome_dir")

    fasta=$(
      {
        find "$genome_dir/raw" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) 2>/dev/null
        find "$genome_dir" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) 2>/dev/null
      } | sort | head -n 1
    )
    [ -n "$fasta" ] || exit 0

    safe=$(printf "%s" "$sample" | tr -c "A-Za-z0-9" "_" | sed "s/_*$//")
    prefix="PROKKA_${safe}"
    tag=$(printf "COLLA%s" "$sample" | tr "[:lower:]" "[:upper:]" | tr -cd "A-Z0-9" | cut -c1-16)

    if [ -f "$genome_dir/prokka/${prefix}.gff" ] && [ -f "$genome_dir/prokka/${prefix}.tsv" ]; then
      printf "Skipping %s; complete Prokka output exists\n" "$sample"
      exit 0
    fi

    printf "Running %s with %s CPUs\n" "$sample" "$cpus_per_job"
    conda run -n "$conda_env" prokka \
      --force \
      --outdir "$genome_dir/prokka" \
      --prefix "$prefix" \
      --locustag "$tag" \
      --cpus "$cpus_per_job" \
      --kingdom Bacteria \
      "$fasta" > "$genome_dir/prokka_run.log" 2>&1
  ' bash {} "$conda_env" "$cpus_per_job"
