# NAFC genomic-neighborhood analysis scripts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23165100.svg)](https://doi.org/10.5281/zenodo.23165100)

Analysis code for the article "Wetland-root microbial genomes contain recurrent multi-gene neighborhoods
with potential relevance to naphthenic acid fractional compound transformation" (Nweze et al.).

- **Author:** Julius Eyiuche Nweze
- **Contact:** julipeale2001@gmail.com
- **Repository:** https://github.com/Julius-Nweze/NAFC_genomic_neighborhoods
- **DOI:** [10.5281/zenodo.23165100](https://doi.org/10.5281/zenodo.23165100)

This repository contains the analysis scripts and a browsable gene-map atlas. Raw genomes, annotations, and
sequencing reads are deposited separately (see the Data Availability section of the article).

## Gene-map atlas

**Browse online:** https://julius-nweze.github.io/NAFC_genomic_neighborhoods/

- `docs/Architecture_gene_maps/`: one gene-map atlas (SVG) and neighborhood table (TSV) per recurring
  gene-set signature, pooling every qualifying genome.
- `docs/svgs_by_genome/`: one combined gene map per genome, showing every candidate neighborhood with at
  least three independent panel hits (288 genomes, 2,095 neighborhoods).

## Scripts: folder layout and run order

Folders are numbered in the order they were used. Within each folder, run the scripts in the order listed.
`shared_dependencies/` is not a step. It holds modules imported by steps 06 to 08.

| Step | Folder | Scripts, in run order |
|---|---|---|
| 01 | `01_MAG_recovery/` | `01_trimmomatic.sh` → `02_megahit_coassembly.sh` → `03_anvio_contigs_database.sh` (contigs ≥1 kbp, contigs database, HMMs) → `04_bowtie2_mapping.sh` → `05_metabat2_binning.sh` → `06_anvio_profile.sh` (SLURM array, one task per sample) → `07_anvio_merge.sh` → `08_anvio_cluster_metabat2.sh` → `09_anvio_summarise.sh` → manual bin refinement in anvi'o → `10_anvio_export_collection.sh` → `11_anvio_import_collection.sh` → `12_checkm2.sh` → `13_gtdbtk.sh` → `14_MAG_phylogenomic_tree.sh` (39 ribosomal genes). |
| 02 | `02_MAG_abundance/` | `01_coverm_genome.sh` (CoverM relative abundance of every MAG per sample) → `02_MAG_abundance.Rmd` (joins abundance with sample metadata, anvi'o bin summaries and GTDB-Tk novelty; plots unmapped reads, MAG counts, phylum- and family-level composition, and genus- and species-level novelty). |
| 03 | `03_genome_annotation/` | `01_Prokka_MAGs.sh` (291 MAGs), `02_Prokka_NCBI_references.sh` (47 comparator genomes), `03_Prokka_OSPW_genomes.sh` and `04_Prokka_OSPW_GROW_genomes.sh` (OSPW-associated genomes). Prokka GFF/TSV coordinates and protein FASTA are used by all later steps. |
| 04 | `04_phylogenomics/` | `00_build_missing_contigs_dbs.sh` → `01_build_external_genomes_and_storage.sh` → `02_run_anvio_hmms_for_17_new_genomes.sh` → `03_build_bacteria71_phylogenomic_tree.sh` (calls `anvio_famsa_stdio_wrapper.sh`). anvi'o v8, Bacteria_71 HMMs, FAMSA alignment, FastTree tree used by iTOL. |
| 05 | `05_reference_panel/` | `fix_KO_proteins_fasta_panel.py`. Reconciles the curated 667-protein panel (`KO_proteins.fasta`) with `Gene_info.xlsx`. |
| 06 | `06_master_pipeline/` | `run_ref_genomes_analysis_pipeline.sh`. BLASTP screen of all 365 genomes, identity/coverage filtering, best-hit resolution, iTOL datasets, aromatics/plastics subcategory split, then step 07. Environment variables: `MIN_PIDENT` (default 30), `MIN_QCOVS` (default 50), `SKIP_BLAST_SEARCH=1` (reuse existing BLAST output), `BLAST_THREADS`. |
| 07 | `07_neighborhood_reconstruction/` | `build_na_gene_panel_operon_clusters.py` → `build_na_gene_panel_cluster_svgs.py`. Called by step 06. Groups best hits into candidate neighborhoods (same contig, gap ≤12 genes and ≤25 kb), expands each to all intervening CDSs, and renders gene maps. |
| 08 | `08_iTOL_gene_count_datasets/` | `category_bundles_corrected/build_corrected_category_itol_heatmaps.py` → `build_category_itol_bundles_min_distinct_gene_pct.py` → `rescope_bundle_shared_annotations.py`, then `build_transport_alkane_alkene_genome_hit_summary.py` → `build_transport_alkane_alkene_gene_taxonomy_table.py` → `build_all_categories_gene_taxonomy_tables.py`. |
| 09 | `09_architecture_signature_atlas/` | `build_architecture_gene_map_figures.py` (imports `render_missing_cluster_svgs.py` and `panel_style.py`) → `build_signature_qualified_neighborhoods_table.py` → `build_curated_figures.py` → `build_atlas_index.py` → `build_svgs_by_genome_index.py`. |
| 10 | `10_annotation_confirmation/` | `export_reference_panel_cluster_reblast_targets.py` → BLASTP against Swiss-Prot (`-outfmt '6 qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'`) → `parse_swissprot_blast_generic.py` → `review_colla_supported_updates.py` → `followup_colla_nohits_local.py` (local IS/AMR BLASTP and HAMAP follow-up of Swiss-Prot no-hits). |
| 11 | `11_treatment_comparison/` | `build_row_ospw_mag_gene_content_comparison.py`. ROW-versus-OSPW MAG detection- and abundance-based comparison. |
| 12 | `12_figure_generation/` | `build_operon_category_summary_chart.py`, `build_figure4_landscape_chart.py` → `build_figure4_composite.py`, `build_cluster_priority_summary_chart.py`, `build_annotation_confirmation_summary.py`. Builds figure panels from tables produced in steps 06 to 11. |
| – | `shared_dependencies/` | `generate_ref_cluster_outputs.py` (shared `enzyme_class()` vocabulary), `summarize_ref_genome_gene_counts_itol.py`, `split_aromatics_itol_by_subcategory.py`. |

## Dependencies

- Trimmomatic, MEGAHIT, Bowtie2, SAMtools, anvi'o v8, MetaBAT2, CheckM2, GTDB-Tk (step 01)
- CoverM (step 02)
- R with `tidyverse`, `readODS`, `reshape2`, `rmarkdown` (step 02)
- Prokka (step 03)
- Python ≥3.10 with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`, `biopython`
- BLAST+ 2.17.0, anvi'o v8, FAMSA v2, FastTree v2, bash

## Inputs (not included)

- Protein FASTA and Prokka annotations (GFF/TSV) for all 365 genomes (291 MAGs, 27 OSPW-associated genomes, 47 comparator genomes)
- `KO_proteins.fasta` (curated panel), `Gene_info.xlsx`, `NCBI_genome_info.xlsx`
- MAG abundance tables (CoverM) for the treatment comparison

## Paths

Input and output locations are set near the top of each script as `/path/to/your/...`. Replace these with the
location of your own data. Cluster scripts also use `<your-account>` and `<your-email>` in their `#SBATCH`
lines. Script-to-script references resolve relative to this folder.

## How to cite

Nweze JE. 2026. NAFC genomic-neighborhood analysis scripts (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23165100

© 2026 Julius Eyiuche Nweze
