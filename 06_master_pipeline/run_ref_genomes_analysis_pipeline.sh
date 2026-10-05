#!/bin/bash
# Title          : run_ref_genomes_analysis_pipeline.sh
# Description    : Master pipeline: BLASTP panel screen, identity/coverage filtering, best-hit resolution, iTOL datasets and neighborhood reconstruction
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/21
# Usage          : ./run_ref_genomes_analysis_pipeline.sh (optional env vars: MIN_PIDENT, MIN_QCOVS, SKIP_BLAST_SEARCH, BLAST_THREADS)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

ANALYSIS_ROOT="/path/to/your/Ref_genomes/Analysis"
GENOME_DIR="$ANALYSIS_ROOT/Genome_proteins"
QUERY_FILE="$ANALYSIS_ROOT/KO/KO_proteins.fasta"

DB_DIR="$ANALYSIS_ROOT/Database"
TXT_DIR="$ANALYSIS_ROOT/TXT"
GENE_HIT_DIR="$ANALYSIS_ROOT/Genes"
RESULT_DIR="$ANALYSIS_ROOT/Result"
ITOL_DIR="$RESULT_DIR/ITOL"
SUMMARY_DIR="$ITOL_DIR/gene_count_summary"
TO_USE_DIR="$ITOL_DIR/To_Use_for_ITOL"
READY_DIR="$TO_USE_DIR/Ready_To_Plot_ITOL"
SHARED_DIR="$READY_DIR/00_full_tree_shared"
BUNDLES_DIR="$READY_DIR/01_category_bundles"
ZERO_DIR="$READY_DIR/99_zero_besthit_no_tree"
AROMATICS_TMP_DIR="$ITOL_DIR/aromatics_split_tmp"
PLASTICS_TMP_DIR="$ITOL_DIR/plastics_split_tmp"
OPERON_CLUSTERS_DIR="$RESULT_DIR/Operon_clusters"

REF_ROOT="/path/to/your/Ref_genomes"
GENE_INFO_SOURCE="$REF_ROOT/Result/Gene_info.xlsx"
METADATA_SOURCE="$REF_ROOT/Result/NCBI_genome_info.xlsx"
TREE_SOURCE="$REF_ROOT/Result/ITOL/Bacteria_71_fasttree.nwk"
SUMMARIZE_SCRIPT="$REPO_ROOT/shared_dependencies/summarize_ref_genome_gene_counts_itol.py"
AROMATICS_SCRIPT="$REPO_ROOT/shared_dependencies/split_aromatics_itol_by_subcategory.py"
OPERON_CLUSTER_SCRIPT="$REPO_ROOT/07_neighborhood_reconstruction/build_nafc_gene_panel_operon_clusters.py"
OPERON_CLUSTER_SVG_SCRIPT="$REPO_ROOT/07_neighborhood_reconstruction/build_nafc_gene_panel_cluster_svgs.py"

THREADS="${BLAST_THREADS:-$(nproc 2>/dev/null || echo 1)}"
MIN_PIDENT="${MIN_PIDENT:-30}"
MIN_QCOVS="${MIN_QCOVS:-50}"
DB_FASTA="$DB_DIR/Genomes_MAGs.faa"
RAW_HITS="$TXT_DIR/KO_proteins.output.txt"
UNIQUE_HITS="$TXT_DIR/KO_proteins.output.unique.txt"
EXTRACTED_HITS="$GENE_HIT_DIR/KO_proteins.faa"
RESULT_RAW="$RESULT_DIR/KO_proteins.output.txt"
RESULT_UNIQUE="$RESULT_DIR/KO_proteins.output.unique.txt"
RESULT_SCORED="$RESULT_DIR/KO_proteins.output.scored.tsv"
TREE_RESULT="$ITOL_DIR/Bacteria_71_fasttree.nwk"
TREE_UPLOAD_COPY="$TO_USE_DIR/NCBI_MAGs_Colla_Bacteria_71_fasttree.nwk"
CHECKLIST="$TO_USE_DIR/CHECKLIST_verified.txt"

echo "Creating analysis output directories..."
mkdir -p \
    "$DB_DIR" \
    "$TXT_DIR" \
    "$GENE_HIT_DIR" \
    "$RESULT_DIR" \
    "$SUMMARY_DIR" \
    "$TO_USE_DIR" \
    "$SHARED_DIR" \
    "$BUNDLES_DIR" \
    "$ZERO_DIR"

if [[ ! -f "$QUERY_FILE" ]]; then
    echo "Missing query file: $QUERY_FILE" >&2
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "Missing genome protein directory: $GENOME_DIR" >&2
    exit 1
fi

SKIP_BLAST_SEARCH="${SKIP_BLAST_SEARCH:-0}"

if [[ "$SKIP_BLAST_SEARCH" == "1" ]]; then
    echo "SKIP_BLAST_SEARCH=1: reusing existing BLAST DB and hit files, rerunning scoring/iTOL only ..."
    for existing in "$RAW_HITS" "$UNIQUE_HITS" "$RESULT_SCORED"; do
        if [[ ! -s "$existing" ]]; then
            echo "SKIP_BLAST_SEARCH=1 but required file is missing/empty: $existing" >&2
            exit 1
        fi
    done
else
    echo "Combining staged genome proteins into $DB_FASTA ..."
    LC_ALL=C find "$GENOME_DIR" -maxdepth 1 -type f -name '*.faa' | sort | while read -r genome_faa; do
        cat "$genome_faa"
    done > "$DB_FASTA"

    echo "Building BLAST database..."
    makeblastdb -in "$DB_FASTA" -dbtype prot -parse_seqids -out "$DB_FASTA" >/dev/null

    echo "Running raw BLAST hits at e-value 1e-25 ..."
    blastp \
        -query "$QUERY_FILE" \
        -db "$DB_FASTA" \
        -out "$RAW_HITS" \
        -evalue 1e-25 \
        -outfmt '6 qseqid sseqid' \
        -num_threads "$THREADS"

    echo "Extracting unique hit proteins..."
    if [[ -s "$RAW_HITS" ]]; then
        awk '{print $2}' "$RAW_HITS" | sort -u > "$UNIQUE_HITS"
        blastdbcmd -db "$DB_FASTA" -entry_batch "$UNIQUE_HITS" -out "$EXTRACTED_HITS" -outfmt '%f'
    else
        : > "$UNIQUE_HITS"
        : > "$EXTRACTED_HITS"
    fi

    echo "Running scored BLAST hits at e-value 1e-25 ..."
    blastp \
        -query "$QUERY_FILE" \
        -db "$DB_FASTA" \
        -out "$RESULT_SCORED" \
        -evalue 1e-25 \
        -outfmt '6 qseqid sseqid pident qcovs length evalue bitscore' \
        -num_threads "$THREADS"
fi

echo "Copying core inputs into Analysis/Result ..."
cp -f "$RAW_HITS" "$RESULT_RAW"
cp -f "$UNIQUE_HITS" "$RESULT_UNIQUE"
cp -f "$GENE_INFO_SOURCE" "$RESULT_DIR/Gene_info.xlsx"
cp -f "$METADATA_SOURCE" "$RESULT_DIR/NCBI_genome_info.xlsx"
cp -f "$TREE_SOURCE" "$TREE_RESULT"

echo "Running summary and iTOL generation (min-pident=$MIN_PIDENT, min-qcovs=$MIN_QCOVS) ..."
python3 "$SUMMARIZE_SCRIPT" \
    --gene-info-xlsx "$RESULT_DIR/Gene_info.xlsx" \
    --ko-output "$RESULT_RAW" \
    --ko-scored-output "$RESULT_SCORED" \
    --metadata-xlsx "$RESULT_DIR/NCBI_genome_info.xlsx" \
    --tree "$TREE_RESULT" \
    --output-dir "$SUMMARY_DIR" \
    --min-pident "$MIN_PIDENT" \
    --min-qcovs "$MIN_QCOVS"

echo "Copying main iTOL outputs ..."
cp -f "$TREE_RESULT" "$TREE_UPLOAD_COPY"
cp -f "$TREE_RESULT" "$SHARED_DIR/Bacteria_71_fasttree.nwk"
cp -f "$TREE_UPLOAD_COPY" "$SHARED_DIR/NCBI_MAGs_Colla_Bacteria_71_fasttree.nwk"
cp -f "$SUMMARY_DIR/itol_phylum_annotation.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_genus_annotation.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_phylum_colorstrip.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_genus_colorstrip.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_source_colorstrip.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_completeness_piechart.txt" "$TO_USE_DIR/" || true
cp -f "$SUMMARY_DIR/itol_genome_name_labels.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_category_counts_besthit.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_subcategory_counts_besthit.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_gene_counts_besthit.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_gene_counts_besthit_by_category.tsv" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_category_counts_strict.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_subcategory_counts_strict.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_gene_counts_strict.txt" "$TO_USE_DIR/"
cp -f "$SUMMARY_DIR/itol_source_legend.tsv" "$TO_USE_DIR/" || true
cp -f "$SUMMARY_DIR/itol_phylum_legend.tsv" "$TO_USE_DIR/" || true
cp -f "$SUMMARY_DIR/itol_genus_legend.tsv" "$TO_USE_DIR/" || true

cp -f "$SUMMARY_DIR/itol_phylum_annotation.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_genus_annotation.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_phylum_colorstrip.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_genus_colorstrip.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_source_colorstrip.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_completeness_piechart.txt" "$SHARED_DIR/" || true
cp -f "$SUMMARY_DIR/itol_genome_name_labels.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_category_counts_besthit.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_subcategory_counts_besthit.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_gene_counts_besthit.txt" "$SHARED_DIR/"
cp -f "$SUMMARY_DIR/itol_gene_counts_besthit_by_category.tsv" "$SHARED_DIR/"

mkdir -p "$TO_USE_DIR/itol_gene_counts_besthit_by_category"
mkdir -p "$TO_USE_DIR/itol_gene_counts_besthit_by_category_trees"
mkdir -p "$TO_USE_DIR/itol_gene_counts_besthit_by_category_annotations"
cp -rf "$SUMMARY_DIR/itol_gene_counts_besthit_by_category/." "$TO_USE_DIR/itol_gene_counts_besthit_by_category/"
cp -rf "$SUMMARY_DIR/itol_gene_counts_besthit_by_category_trees/." "$TO_USE_DIR/itol_gene_counts_besthit_by_category_trees/"
cp -rf "$SUMMARY_DIR/itol_gene_counts_besthit_by_category_annotations/." "$TO_USE_DIR/itol_gene_counts_besthit_by_category_annotations/"

echo "Creating category-ready bundle folders ..."
export SUMMARY_DIR TO_USE_DIR SHARED_DIR BUNDLES_DIR ZERO_DIR TREE_RESULT
python3 - <<'PY'
import csv
import os
import shutil
from pathlib import Path

def slugify(label: str) -> str:
    import re
    slug = re.sub(r"[^A-Za-z0-9]+", "_", label.strip().lower()).strip("_")
    return slug or "unassigned"

summary_dir = Path(os.environ["SUMMARY_DIR"])
to_use_dir = Path(os.environ["TO_USE_DIR"])
bundles_dir = Path(os.environ["BUNDLES_DIR"])
zero_dir = Path(os.environ["ZERO_DIR"])

manifest_path = summary_dir / "itol_gene_counts_besthit_by_category.tsv"
with manifest_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for index, row in enumerate(reader, start=1):
        category = row["Category"]
        slug = slugify(category)
        heatmap = row["HeatmapFileName"]
        tree_file = row.get("TreeFileName", "")
        bundle_name = f"{index:02d}_{slug}"
        if tree_file:
            bundle_dir = bundles_dir / bundle_name
            bundle_dir.mkdir(parents=True, exist_ok=True)
            file_names = [
                ("itol_gene_counts_besthit_by_category", heatmap),
                ("itol_gene_counts_besthit_by_category_trees", tree_file),
                ("itol_gene_counts_besthit_by_category_annotations", row["PhylumAnnotationFileName"]),
                ("itol_gene_counts_besthit_by_category_annotations", row["GenusAnnotationFileName"]),
                ("itol_gene_counts_besthit_by_category_annotations", row["SourceStripFileName"]),
                ("itol_gene_counts_besthit_by_category_annotations", row["GenomeNameLabelsFileName"]),
                ("itol_gene_counts_besthit_by_category_annotations", row.get("CompletenessPiechartFileName", "")),
            ]
            compat_files = [
                f"itol_phylum_colorstrip__{slug}.txt",
                f"itol_genus_colorstrip__{slug}.txt",
            ]
            for folder, name in file_names:
                if not name:
                    continue
                source = to_use_dir / folder / name
                if source.exists():
                    shutil.copy2(source, bundle_dir / name)
            for name in compat_files:
                source = to_use_dir / "itol_gene_counts_besthit_by_category_annotations" / name
                if source.exists():
                    shutil.copy2(source, bundle_dir / name)
        else:
            source = to_use_dir / "itol_gene_counts_besthit_by_category" / heatmap
            if source.exists():
                zero_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, zero_dir / heatmap)
PY

echo "Running Aromatics split helper (min-pident=$MIN_PIDENT, min-qcovs=$MIN_QCOVS) ..."
python3 "$AROMATICS_SCRIPT" \
    --gene-info-xlsx "$RESULT_DIR/Gene_info.xlsx" \
    --ko-output "$RESULT_RAW" \
    --ko-scored-output "$RESULT_SCORED" \
    --metadata-xlsx "$RESULT_DIR/NCBI_genome_info.xlsx" \
    --tree "$TREE_RESULT" \
    --output-dir "$AROMATICS_TMP_DIR" \
    --min-pident "$MIN_PIDENT" \
    --min-qcovs "$MIN_QCOVS"

cp -f "$AROMATICS_TMP_DIR/itol_gene_counts_besthit__aromatics_split_by_subcategory.tsv" "$TO_USE_DIR/"
cp -f "$AROMATICS_TMP_DIR/README_aromatics_split.txt" "$TO_USE_DIR/"
cp -f "$AROMATICS_TMP_DIR/itol_gene_counts_besthit__aromatics_split_by_subcategory.tsv" "$SUMMARY_DIR/"
cp -f "$AROMATICS_TMP_DIR/README_aromatics_split.txt" "$SUMMARY_DIR/"
cp -f "$AROMATICS_TMP_DIR/itol_gene_counts_besthit__aromatics_split_by_subcategory.tsv" "$BUNDLES_DIR/"
cp -f "$AROMATICS_TMP_DIR/README_aromatics_split.txt" "$BUNDLES_DIR/"
for bundle_name in 04_aromatics_1 04_aromatics_2 04_aromatics_3 04_aromatics_4 04_aromatics_5; do
    mkdir -p "$BUNDLES_DIR/$bundle_name"
    cp -f "$AROMATICS_TMP_DIR/bundles/$bundle_name/"* "$BUNDLES_DIR/$bundle_name/"
done

echo "Running Plastics split helper (min-pident=$MIN_PIDENT, min-qcovs=$MIN_QCOVS) ..."
python3 "$AROMATICS_SCRIPT" \
    --gene-info-xlsx "$RESULT_DIR/Gene_info.xlsx" \
    --ko-output "$RESULT_RAW" \
    --ko-scored-output "$RESULT_SCORED" \
    --metadata-xlsx "$RESULT_DIR/NCBI_genome_info.xlsx" \
    --tree "$TREE_RESULT" \
    --category "Plastics" \
    --groups 3 \
    --bundle-prefix 10 \
    --output-dir "$PLASTICS_TMP_DIR" \
    --min-pident "$MIN_PIDENT" \
    --min-qcovs "$MIN_QCOVS"

cp -f "$PLASTICS_TMP_DIR/itol_gene_counts_besthit__plastics_split_by_subcategory.tsv" "$TO_USE_DIR/"
cp -f "$PLASTICS_TMP_DIR/README_plastics_split.txt" "$TO_USE_DIR/"
cp -f "$PLASTICS_TMP_DIR/itol_gene_counts_besthit__plastics_split_by_subcategory.tsv" "$SUMMARY_DIR/"
cp -f "$PLASTICS_TMP_DIR/README_plastics_split.txt" "$SUMMARY_DIR/"
cp -f "$PLASTICS_TMP_DIR/itol_gene_counts_besthit__plastics_split_by_subcategory.tsv" "$BUNDLES_DIR/"
cp -f "$PLASTICS_TMP_DIR/README_plastics_split.txt" "$BUNDLES_DIR/"
for bundle_name in 10_plastics_1 10_plastics_2 10_plastics_3; do
    mkdir -p "$BUNDLES_DIR/$bundle_name"
    cp -f "$PLASTICS_TMP_DIR/bundles/$bundle_name/"* "$BUNDLES_DIR/$bundle_name/"
done

echo "Running reference-panel operon/genomic-cluster context check ..."
python3 "$OPERON_CLUSTER_SCRIPT" \
    --besthit-table "$SUMMARY_DIR/besthit_assigned_hit_proteins.tsv" \
    --output-dir "$OPERON_CLUSTERS_DIR"

export OPERON_CLUSTERS_DIR
python3 - <<'PY'
import csv
import os
from pathlib import Path

operon_dir = Path(os.environ["OPERON_CLUSTERS_DIR"])
summary_path = operon_dir / "NAFC_gene_panel_neighborhood_summary.tsv"
with summary_path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
top = [r for r in rows if r["Cluster strength"] == "multi_hit" and int(r["Panel seed hit count"]) >= 3]
top.sort(key=lambda r: int(r["Panel seed hit count"]), reverse=True)
fields = list(rows[0].keys()) if rows else []
out_path = operon_dir / "NAFC_gene_panel_top_multi_hit_neighborhoods.tsv"
with out_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(top)
print(f"Top multi-hit shortlist (>=3 panel genes): {len(top)} clusters -> {out_path}")
PY

echo "Copying operon/cluster context summary into the iTOL bundle root ..."
cp -f "$OPERON_CLUSTERS_DIR/README_operon_cluster_summary.txt" "$TO_USE_DIR/"
cp -f "$OPERON_CLUSTERS_DIR/README_operon_cluster_summary.txt" "$SHARED_DIR/"
cp -f "$OPERON_CLUSTERS_DIR/NAFC_gene_panel_top_multi_hit_neighborhoods.tsv" "$TO_USE_DIR/"

cat > "$CHECKLIST" <<EOF
Run settings

Analysis rerun root: $ANALYSIS_ROOT
BLAST query file: $QUERY_FILE
BLAST database FASTA: $DB_FASTA
BLAST e-value cutoff: 1e-25
Scored output columns: qseqid sseqid pident qcovs length evalue bitscore

Identity/coverage floor:
- A hit only counts toward gene/category/heatmap totals if pident >= $MIN_PIDENT and qcovs >= $MIN_QCOVS,
  in addition to the e-value <= 1e-25 cutoff from the BLASTP run itself.
- Rationale: e-value alone does not guard against low-identity cross-hits within large enzyme
  superfamilies (P450s, Rieske oxygenases, dehydrogenases) that dominate these categories.
- Per-run filtered-row counts are recorded in gene_count_summary/README_summary.txt.
- Override via MIN_PIDENT / MIN_QCOVS environment variables before rerunning this script.

Use phylum and genus taxonomy files as TREE_COLORS annotations, not as DATASET_COLORSTRIP files.
The preferred files are:
- itol_phylum_annotation.txt
- itol_genus_annotation.txt
- per-category itol_phylum_annotation__<category>.txt
- per-category itol_genus_annotation__<category>.txt

Source remains a DATASET_COLORSTRIP file.
Genome names remain a LABELS file.

Phylum colors use a fixed hex palette (PHYLUM_COLORS in summarize_ref_genome_gene_counts_itol.py).
Phyla not in that palette get an auto-generated color that does not collide with it.
Genus is auto-colored.

Genome completeness is a DATASET_PIECHART file (itol_completeness_piechart.txt,
and itol_completeness_piechart__<category>.txt per category bundle): an external 2-slice
pie per genome (radius 1, red=Completeness%, green=100-Completeness%), sourced from the
"Completeness (%)" column of NCBI_genome_info.xlsx. Use it to sanity-check gene-absence calls on
genomes below ~70% CheckM completeness before drawing presence/absence conclusions.

Aromatics split bundles:
- Bundles 04_aromatics_1 to 04_aromatics_5.
- Manifest: itol_gene_counts_besthit__aromatics_split_by_subcategory.tsv
- Bundle root: Ready_To_Plot_ITOL/01_category_bundles/

Plastics split bundles:
- Bundles 10_plastics_1 to 10_plastics_3 (split_aromatics_itol_by_subcategory.py invoked with
  --category Plastics --groups 3 --bundle-prefix 10).
- Manifest: itol_gene_counts_besthit__plastics_split_by_subcategory.tsv
- Bundle root: Ready_To_Plot_ITOL/01_category_bundles/

Candidate genomic neighborhoods:
- Groups panel gene hits into candidate neighborhoods (same contig, gene-order gap <=12,
  coordinate gap <=25 kb) and classifies every gene with the shared enzyme_class() vocabulary.
- Script: build_nafc_gene_panel_operon_clusters.py
- Full tables (per-cluster summary, per-gene members, top >=3-hit shortlist) live in
  Analysis/Result/Operon_clusters/ - not copied here in full because they are large
  (tens of MB) and are not iTOL dataset files.
- README_operon_cluster_summary.txt (run-level stats) and NAFC_gene_panel_top_multi_hit_neighborhoods.tsv
  (shortlist of multi-gene clusters) ARE copied here, into
  To_Use_for_ITOL/ and this SHARED_DIR, so anyone starting from the iTOL bundle can see the
  operon-context headline numbers and shortlist without leaving this folder tree.
- A cluster's "Panel categories" column matches the category names used for the bundle
  folders above (e.g. Aromatics, Plastics), so the shortlist can be cross-referenced
  directly against 01_category_bundles/.
- Both build_nafc_gene_panel_operon_clusters.py and build_nafc_gene_panel_cluster_svgs.py
  (the SVG gene-map renderer for the >=3-hit shortlist; not run automatically on every
  rerun since it writes ~2,000 files) are copied into Analysis/Result/Operon_clusters/
  alongside their outputs, matching how SUMMARIZE_SCRIPT/AROMATICS_SCRIPT are copied into
  gene_count_summary/. Rendered SVGs + a browsable index.html live in
  Analysis/Result/Operon_clusters/svgs/.

Category heatmaps preserve local workbook row order from Gene_info.xlsx, not first global gene occurrence.
Within a subcategory, repeated workbook rows are collapsed to one gene column.
A gene still repeats across different subcategory blocks when Gene_info.xlsx places it in more than one subcategory.
EOF
cp -f "$CHECKLIST" "$SHARED_DIR/CHECKLIST_verified.txt"

cp -f "$SUMMARIZE_SCRIPT" "$SUMMARY_DIR/"
cp -f "$AROMATICS_SCRIPT" "$SUMMARY_DIR/"
cp -f "$OPERON_CLUSTER_SCRIPT" "$OPERON_CLUSTERS_DIR/"
cp -f "$OPERON_CLUSTER_SVG_SCRIPT" "$OPERON_CLUSTERS_DIR/"

cat > "$RESULT_DIR/README_analysis_rerun.txt" <<EOF
Ref_genomes analysis run completed.

Inputs:
- Query FASTA: $QUERY_FILE
- Genome proteins: $GENOME_DIR
- Gene info workbook: $RESULT_DIR/Gene_info.xlsx
- Genome metadata workbook: $RESULT_DIR/NCBI_genome_info.xlsx
- Tree: $TREE_RESULT

Core outputs:
- Raw BLAST hits: $RESULT_RAW
- Unique BLAST hit subjects: $RESULT_UNIQUE
- Extracted hit proteins: $EXTRACTED_HITS
- Scored BLAST hits: $RESULT_SCORED
- Summary/iTOL outputs: $SUMMARY_DIR
- Upload-ready bundle root: $READY_DIR

Gene-count filtering:
- Identity/coverage floor: pident >= $MIN_PIDENT, qcovs >= $MIN_QCOVS (on top of e-value <= 1e-25).
- See gene_count_summary/README_summary.txt for how many rows this dropped.
- Genome completeness is now included as an iTOL DATASET_PIECHART track (whole-tree and per-category),
  matching the style of Analysis/Epibolus/3_ITOL_MAG_Completion_E.pulchripes.txt.
- Phylum ring colors now use the fixed palette from Summary_Bins3/Plot/MAG_abundance.Rmd.
EOF

echo "Analysis rerun complete."
echo "Raw hits: $RESULT_RAW"
echo "Scored hits: $RESULT_SCORED"
echo "Summary dir: $SUMMARY_DIR"
echo "Ready-to-plot iTOL bundles: $READY_DIR"
