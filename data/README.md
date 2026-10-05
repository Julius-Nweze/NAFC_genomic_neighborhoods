# Data files

| Folder | Contents |
|---|---|
| `reference_panel/` | `KO_proteins.fasta`: the curated panel of 667 reference proteins in 22 functional categories |
| `blastp_besthit/` | `besthit_assigned_hit_proteins.tsv`: one row per genome protein (28,490) assigned to its best panel query after identity (>=30%) and coverage (>=50%) filtering; `besthit_resolved_ambiguous_hit_proteins.tsv`: proteins that matched several panel entries and how each was resolved |
| `neighborhoods/` | `NAFC_gene_panel_neighborhood_summary.tsv`: one row per candidate neighborhood (18,029); `NAFC_gene_panel_neighborhood_members.tsv`: every gene in each expanded neighborhood span with its functional class; `NAFC_gene_panel_top_multi_hit_neighborhoods.tsv`: neighborhoods with >=3 independent panel hits (2,095) |
| `itol/` | Unified FastTree phylogenomic tree of 364 genomes (Exp3_MAG_2 excluded); `Bacteria_71_fasttree.nwk` and `NCBI_MAGs_Colla_Bacteria_71_fasttree.nwk` are identical copies under the two names the scripts use and iTOL annotation datasets (phylum, genus and genome-source colour strips, genome labels, completeness pie charts, gene and category heatmaps); `category_heatmaps/` holds one heatmap dataset per functional category |
| `prokka/` | Prokka annotation of all 365 genomes: `<genome>.gff.gz` (feature coordinates, without the appended genome sequence), `<genome>.tsv.gz` (feature table) and `<genome>.faa.gz` (protein sequences) |

## Genome identifiers

- `Exp3_MAG_<n>`: plant-root metagenome-assembled genomes (291).
- `GROW_<n>_consensus`, `<n>.medaka`, `L<n>.medaka`: OSPW-associated genomes (27).
- `GCF_...`, `CP...`, `AP...`, `LMAZ...`: NCBI reference genomes (47). The best-hit table writes
  these without the version suffix (e.g. `GCF_000267545` for `GCF_000267545.1`).

In the `Tier` and source columns, `MAG` = plant-root MAG, `Colla` = OSPW-associated genome, and
`Ref` = NCBI reference genome.
