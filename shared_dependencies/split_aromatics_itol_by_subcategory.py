#!/usr/bin/env python3
# Title          : split_aromatics_itol_by_subcategory.py
# Description    : Split Aromatics and Plastics iTOL datasets into subcategory panels (called by the master pipeline)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/20
# Usage          : python3 split_aromatics_itol_by_subcategory.py [options]

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

from Bio import Phylo

import summarize_ref_genome_gene_counts_itol as base


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split a category's iTOL gene heatmap into smaller subcategory-based "
            "bundles with matching pruned trees and annotation files. Originally written "
            "for Aromatics (the default --category), but works for any category in "
            "Gene_info.xlsx via --category/--groups/--bundle-prefix."
        )
    )
    parser.add_argument(
        "--gene-info-xlsx",
        default="/path/to/your/Ref_genomes/Result/Gene_info.xlsx",
    )
    parser.add_argument(
        "--ko-output",
        default="/path/to/your/Ref_genomes/Result/KO_proteins.output.txt",
    )
    parser.add_argument(
        "--ko-scored-output",
        default="/path/to/your/Ref_genomes/Result/KO_proteins.output.scored.tsv",
    )
    parser.add_argument(
        "--metadata-xlsx",
        default="/path/to/your/Ref_genomes/Result/NCBI_genome_info.xlsx",
    )
    parser.add_argument(
        "--tree",
        default="/path/to/your/Ref_genomes/Result/ITOL/Bacteria_71_fasttree.nwk",
    )
    parser.add_argument(
        "--category",
        default="Aromatics",
        help="Category to split.",
    )
    parser.add_argument(
        "--groups",
        type=int,
        default=5,
        help="Number of subcategory groups to create.",
    )
    parser.add_argument(
        "--bundle-prefix",
        default="04",
        help=(
            "Numeric prefix matching this category's position in "
            "itol_gene_counts_besthit_by_category.tsv / the 01_category_bundles folder "
            "(e.g. '10' for Plastics, which sits at 10_plastics). Used to name split "
            "bundle folders as '<prefix>_<category_slug>_<group>' consistently with the "
            "single-bundle folder naming."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="/tmp/aromatics_itol_split",
        help="Temporary output directory for the split bundles.",
    )
    parser.add_argument(
        "--min-pident",
        type=float,
        default=30.0,
        help="Minimum percent identity floor, passed through to process_scored_ko_output.",
    )
    parser.add_argument(
        "--min-qcovs",
        type=float,
        default=50.0,
        help="Minimum query coverage floor, passed through to process_scored_ko_output.",
    )
    parser.add_argument(
        "--filtered-min-distinct-fraction",
        type=float,
        default=0.0,
        help=(
            "Optional fraction (0-1) used to create an extra filtered copy of the "
            "subcategory bundles that keeps only genomes with enough distinct genes "
            "in that bundle. The threshold is ceil(distinct bundle gene count * "
            "fraction). Set to 0 to disable this extra output."
        ),
    )
    return parser.parse_args()


def load_counts(args: argparse.Namespace):
    _header, raw_rows = base.load_xlsx_rows(Path(args.gene_info_xlsx).resolve())
    gene_info_rows = [{key: base.normalize_text(value) for key, value in row.items()} for row in raw_rows]
    (
        gene_order,
        gene_meta,
        by_uniprot,
        by_gene,
        _category_gene_row_order,
        _subcategory_gene_row_order,
    ) = base.read_gene_info(Path(args.gene_info_xlsx).resolve())
    tree_path = Path(args.tree).resolve()
    tree_labels = base.read_tree_labels(tree_path)
    metadata_header, metadata_lookup = base.read_metadata(Path(args.metadata_xlsx).resolve())

    missing_tree_metadata = [label for label in tree_labels if label not in metadata_lookup]
    if missing_tree_metadata:
        raise ValueError(
            "Tree labels missing from metadata workbook: "
            + ", ".join(missing_tree_metadata[:20])
        )

    (
        query_cache,
        hit_to_genes,
        hit_to_queries,
        hit_to_genome,
        genome_hit_counter,
        _raw_total_rows,
        _scored_total_rows,
        besthit_by_hit,
        _filtered_low_identity_rows,
        _filtered_low_coverage_rows,
    ) = base.process_scored_ko_output(
        Path(args.ko_output).resolve(),
        Path(args.ko_scored_output).resolve(),
        by_uniprot,
        by_gene,
        min_pident=args.min_pident,
        min_qcovs=args.min_qcovs,
    )

    allowed_tree_labels = set(tree_labels)
    filtered_hit_ids = {
        hit_id for hit_id, genome_id in hit_to_genome.items() if genome_id in allowed_tree_labels
    }
    hit_to_genes = {
        hit_id: genes for hit_id, genes in hit_to_genes.items() if hit_id in filtered_hit_ids
    }
    hit_to_queries = {
        hit_id: queries for hit_id, queries in hit_to_queries.items() if hit_id in filtered_hit_ids
    }
    besthit_by_hit = {
        hit_id: best_row for hit_id, best_row in besthit_by_hit.items() if hit_id in filtered_hit_ids
    }
    hit_to_genome = {
        hit_id: genome_id for hit_id, genome_id in hit_to_genome.items() if hit_id in filtered_hit_ids
    }
    genome_hit_counter = Counter(
        {
            genome_id: count
            for genome_id, count in genome_hit_counter.items()
            if genome_id in allowed_tree_labels
        }
    )

    counts = base.count_hits(
        tree_labels,
        gene_order,
        gene_meta,
        hit_to_genes,
        hit_to_genome,
        hit_to_queries,
        query_cache,
        besthit_by_hit,
    )

    phylum_colors = base.make_phylum_color_map(
        [
            base.normalize_text(metadata_lookup[genome].get("Phylum", "")) or "Unassigned"
            for genome in tree_labels
        ],
    )
    genus_colors = base.make_color_map(
        [
            base.normalize_text(metadata_lookup[genome].get("Genus", "")) or "Unassigned"
            for genome in tree_labels
        ],
        prefer_high_contrast=True,
    )

    source_colors = None
    if "Source" in metadata_header:
        source_colors = base.make_color_map(
            [
                base.normalize_text(metadata_lookup[genome].get("Source", "")) or "Unassigned"
                for genome in tree_labels
            ]
        )

    completeness_available = "Completeness (%)" in metadata_header

    return (
        tree_path,
        tree_labels,
        metadata_lookup,
        gene_meta,
        gene_info_rows,
        counts,
        phylum_colors,
        genus_colors,
        source_colors,
        completeness_available,
    )


def get_category_subcategory_maps(
    category: str,
    gene_info_rows: list[dict[str, str]],
):
    subcategory_order: list[str] = []
    subcategory_to_genes: dict[str, set[str]] = defaultdict(set)
    category_gene_rows: list[str] = []
    seen_category_subcategory_gene: set[tuple[str, str, str]] = set()

    for row in gene_info_rows:
        gene = base.normalize_text(row.get("Gene", ""))
        if not gene:
            continue
        categories = base.split_multi_value(base.normalize_text(row.get("Category", "")))
        if category not in categories:
            continue
        subcategories = base.split_multi_value(
            base.normalize_text(row.get("Sub category", ""))
        ) or ["Unassigned"]
        for subcategory in subcategories:
            key = (category, subcategory, gene)
            if key not in seen_category_subcategory_gene:
                category_gene_rows.append(gene)
                seen_category_subcategory_gene.add(key)
            if subcategory not in subcategory_order:
                subcategory_order.append(subcategory)
            subcategory_to_genes[subcategory].add(gene)

    if not category_gene_rows:
        raise ValueError(f"No genes found for category {category}")

    return subcategory_order, subcategory_to_genes, category_gene_rows


def get_group_gene_rows(
    category: str,
    subcategories: list[str],
    gene_info_rows: list[dict[str, str]],
) -> list[str]:
    allowed_subcategories = set(subcategories)
    group_gene_rows: list[str] = []
    seen_subcategory_gene: set[tuple[str, str]] = set()
    for row in gene_info_rows:
        gene = base.normalize_text(row.get("Gene", ""))
        if not gene:
            continue
        categories = base.split_multi_value(base.normalize_text(row.get("Category", "")))
        if category not in categories:
            continue
        row_subcategories = base.split_multi_value(base.normalize_text(row.get("Sub category", ""))) or [
            "Unassigned"
        ]
        for subcategory in row_subcategories:
            if subcategory not in allowed_subcategories:
                continue
            key = (subcategory, gene)
            if key not in seen_subcategory_gene:
                group_gene_rows.append(gene)
                seen_subcategory_gene.add(key)
    return group_gene_rows


def partition_subcategories(
    subcategory_order: list[str],
    subcategory_to_genes: dict[str, set[str]],
    group_count: int,
) -> list[dict[str, object]]:
    if group_count < 1:
        raise ValueError("group_count must be at least 1")
    if group_count > len(subcategory_order):
        group_count = len(subcategory_order)

    size = len(subcategory_order)
    segment_gene_counts = [[0] * (size + 1) for _ in range(size)]
    for left in range(size):
        genes: set[str] = set()
        for right in range(left, size):
            genes.update(subcategory_to_genes[subcategory_order[right]])
            segment_gene_counts[left][right + 1] = len(genes)

    all_genes: set[str] = set()
    for subcategory in subcategory_order:
        all_genes.update(subcategory_to_genes[subcategory])
    target = len(all_genes) / group_count

    dp = [[math.inf] * (size + 1) for _ in range(group_count + 1)]
    previous = [[-1] * (size + 1) for _ in range(group_count + 1)]
    dp[0][0] = 0.0

    for groups_used in range(1, group_count + 1):
        for right in range(1, size + 1):
            for left in range(groups_used - 1, right):
                gene_count = segment_gene_counts[left][right]
                score = dp[groups_used - 1][left] + (gene_count - target) ** 2
                if score < dp[groups_used][right]:
                    dp[groups_used][right] = score
                    previous[groups_used][right] = left

    ranges = []
    right = size
    for groups_used in range(group_count, 0, -1):
        left = previous[groups_used][right]
        if left < 0:
            raise ValueError("Failed to partition subcategories")
        ranges.append((left, right))
        right = left
    ranges.reverse()

    partitions = []
    for index, (left, right) in enumerate(ranges, start=1):
        subcategories = subcategory_order[left:right]
        genes: set[str] = set()
        for subcategory in subcategories:
            genes.update(subcategory_to_genes[subcategory])
        partitions.append(
            {
                "group_index": index,
                "subcategories": subcategories,
                "gene_count": len(genes),
            }
        )
    return partitions


def prune_tree(tree_text: str, keep_labels: set[str], output_path: Path) -> None:
    tree = Phylo.read(StringIO(tree_text), "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep_labels:
            tree.prune(target=terminal)
    Phylo.write(tree, output_path, "newick")


def format_percent_slug(value: float) -> str:
    percent = value * 100.0
    if abs(percent - round(percent)) < 1e-9:
        return f"{int(round(percent))}pct"
    return f"{f'{percent:.3f}'.rstrip('0').rstrip('.').replace('.', 'p')}pct"


def write_bundle(
    bundle_dir: Path,
    slug: str,
    category: str,
    subcategories: list[str],
    tree_path: Path,
    tree_text: str,
    tree_labels: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    besthit_gene_counts: dict[str, Counter[str]],
    category_genes: list[str],
    phylum_colors: dict[str, str],
    genus_colors: dict[str, str],
    source_colors: dict[str, str] | None,
    min_distinct_fraction: float = 0.0,
    completeness_available: bool = False,
) -> dict[str, str]:
    bundle_dir.mkdir(parents=True, exist_ok=True)

    distinct_genes = base.ordered_unique(category_genes)
    genomes_with_any_hits = [
        genome
        for genome in tree_labels
        if any(besthit_gene_counts.get(genome, Counter()).get(gene, 0) > 0 for gene in distinct_genes)
    ]
    min_distinct_genes_required = 1
    cutoff_tag = ""
    if min_distinct_fraction > 0:
        min_distinct_genes_required = max(1, math.ceil(len(distinct_genes) * min_distinct_fraction))
        cutoff_tag = f"min{min_distinct_genes_required}of{len(distinct_genes)}"
        genomes_with_hits = [
            genome
            for genome in tree_labels
            if sum(
                1
                for gene in distinct_genes
                if besthit_gene_counts.get(genome, Counter()).get(gene, 0) > 0
            )
            >= min_distinct_genes_required
        ]
    else:
        genomes_with_hits = genomes_with_any_hits
    if not genomes_with_hits:
        raise ValueError(f"No genomes with hits for {slug}")

    file_slug = slug if not cutoff_tag else f"{slug}__{cutoff_tag}"
    heatmap_name = f"itol_gene_counts_besthit__{file_slug}.txt"
    tree_name = f"{tree_path.stem}__{file_slug}.nwk"
    phylum_annotation_name = f"itol_phylum_annotation__{file_slug}.txt"
    genus_annotation_name = f"itol_genus_annotation__{file_slug}.txt"
    phylum_compat_name = f"itol_phylum_colorstrip__{file_slug}.txt"
    genus_compat_name = f"itol_genus_colorstrip__{file_slug}.txt"
    source_name = f"itol_source_colorstrip__{file_slug}.txt"
    labels_name = f"itol_genome_name_labels__{file_slug}.txt"
    completeness_name = f"itol_completeness_piechart__{file_slug}.txt"

    base.write_itol_heatmap(
        bundle_dir / heatmap_name,
        (
            f"gene_counts_besthit__{category}__{' | '.join(subcategories)}"
            if not cutoff_tag
            else f"gene_counts_besthit__{category}__{' | '.join(subcategories)}__{cutoff_tag}"
        ),
        genomes_with_hits,
        category_genes if not cutoff_tag else distinct_genes,
        besthit_gene_counts,
        show_labels=True,
        fixed_max_value=base.CATEGORY_HEATMAP_SHARED_MAX,
    )

    tree_file = bundle_dir / tree_name
    prune_tree(tree_text, set(genomes_with_hits), tree_file)

    base.write_itol_tree_colors(
        bundle_dir / phylum_annotation_name,
        tree_file,
        metadata_lookup,
        "Phylum",
        phylum_colors,
    )
    base.write_itol_tree_colors(
        bundle_dir / genus_annotation_name,
        tree_file,
        metadata_lookup,
        "Genus",
        genus_colors,
    )
    base.write_itol_tree_colors(
        bundle_dir / phylum_compat_name,
        tree_file,
        metadata_lookup,
        "Phylum",
        phylum_colors,
    )
    base.write_itol_tree_colors(
        bundle_dir / genus_compat_name,
        tree_file,
        metadata_lookup,
        "Genus",
        genus_colors,
    )

    if source_colors is not None:
        base.write_itol_colorstrip(
            bundle_dir / source_name,
            f"Genome source __ {category} __ group {slug}",
            genomes_with_hits,
            metadata_lookup,
            "Source",
            strip_width=28,
            color_map=source_colors,
        )

    base.write_itol_labels(
        bundle_dir / labels_name,
        genomes_with_hits,
        metadata_lookup,
        "Genome name",
    )

    if completeness_available:
        base.write_itol_piechart(
            bundle_dir / completeness_name,
            f"Genome completeness (%) __ {category} __ group {slug}",
            genomes_with_hits,
            metadata_lookup,
            "Completeness (%)",
        )

    return {
        "Group": slug,
        "SubcategoryCount": str(len(subcategories)),
        "Subcategories": "; ".join(subcategories),
        "GeneCount": str(len(category_genes)),
        "DistinctGeneCount": str(len(distinct_genes)),
        "MinDistinctGenesRequired": str(min_distinct_genes_required),
        "GenomeCountWithAnyHits": str(len(genomes_with_any_hits)),
        "GenomeCount": str(len(genomes_with_hits)),
        "CutoffTag": cutoff_tag,
        "TreeFileName": tree_name,
        "HeatmapFileName": heatmap_name,
        "PhylumAnnotationFileName": phylum_annotation_name,
        "GenusAnnotationFileName": genus_annotation_name,
        "SourceStripFileName": source_name if source_colors is not None else "",
        "GenomeNameLabelsFileName": labels_name,
        "CompletenessPiechartFileName": completeness_name if completeness_available else "",
        "BundleFolder": bundle_dir.name,
    }


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.filtered_min_distinct_fraction <= 1.0:
        raise ValueError("--filtered-min-distinct-fraction must be between 0 and 1")
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    (
        tree_path,
        tree_labels,
        metadata_lookup,
        gene_meta,
        gene_info_rows,
        counts,
        phylum_colors,
        genus_colors,
        source_colors,
        completeness_available,
    ) = load_counts(args)

    subcategory_order, subcategory_to_genes, _category_gene_rows = get_category_subcategory_maps(
        args.category,
        gene_info_rows,
    )
    partitions = partition_subcategories(subcategory_order, subcategory_to_genes, args.groups)

    category_slug = base.slugify_label(args.category)
    tree_text = tree_path.read_text(encoding="utf-8")
    manifest_rows = []
    for partition in partitions:
        group_index = partition["group_index"]
        subcategories = partition["subcategories"]
        slug = f"{category_slug}_{group_index}"
        bundle_name = f"{args.bundle_prefix}_{category_slug}_{group_index}"
        genes = get_group_gene_rows(args.category, subcategories, gene_info_rows)
        manifest_rows.append(
            write_bundle(
                output_dir / "bundles" / bundle_name,
                slug,
                args.category,
                subcategories,
                tree_path,
                tree_text,
                tree_labels,
                metadata_lookup,
                counts["besthit_gene_counts"],
                genes,
                phylum_colors,
                genus_colors,
                source_colors,
                completeness_available=completeness_available,
            )
        )

    manifest_name = f"itol_gene_counts_besthit__{category_slug}_split_by_subcategory.tsv"
    base.write_rows(
        output_dir / manifest_name,
        manifest_rows,
        [
            "Group",
            "SubcategoryCount",
            "Subcategories",
            "GeneCount",
            "DistinctGeneCount",
            "MinDistinctGenesRequired",
            "GenomeCountWithAnyHits",
            "GenomeCount",
            "CutoffTag",
            "TreeFileName",
            "HeatmapFileName",
            "PhylumAnnotationFileName",
            "GenusAnnotationFileName",
            "SourceStripFileName",
            "GenomeNameLabelsFileName",
            "CompletenessPiechartFileName",
            "BundleFolder",
        ],
    )

    readme_lines = [
        f"Category: {args.category}",
        f"Requested groups: {args.groups}",
        f"Created groups: {len(manifest_rows)}",
        "",
        "Group summary:",
    ]
    for row in manifest_rows:
        readme_lines.append(
            f"- {row['Group']}: {row['GeneCount']} genes, {row['GenomeCount']} genomes, "
            f"{row['SubcategoryCount']} subcategories"
        )
        readme_lines.append(f"  Subcategories: {row['Subcategories']}")
    readme_name = f"README_{category_slug}_split.txt"
    (output_dir / readme_name).write_text(
        "\n".join(readme_lines) + "\n",
        encoding="utf-8",
    )

    filtered_manifest_name = ""
    filtered_readme_name = ""
    filtered_bundle_root_name = ""
    if args.filtered_min_distinct_fraction > 0:
        filtered_suffix = f"min{format_percent_slug(args.filtered_min_distinct_fraction)}_distinct_genes"
        filtered_bundle_root = output_dir / f"bundles_{filtered_suffix}"
        filtered_bundle_root.mkdir(parents=True, exist_ok=True)

        filtered_manifest_rows = []
        filtered_skipped_groups = []
        for partition in partitions:
            group_index = partition["group_index"]
            subcategories = partition["subcategories"]
            slug = f"{category_slug}_{group_index}"
            bundle_name = f"{args.bundle_prefix}_{category_slug}_{group_index}"
            genes = get_group_gene_rows(args.category, subcategories, gene_info_rows)
            try:
                filtered_manifest_rows.append(
                    write_bundle(
                        filtered_bundle_root / bundle_name,
                        slug,
                        args.category,
                        subcategories,
                        tree_path,
                        tree_text,
                        tree_labels,
                        metadata_lookup,
                        counts["besthit_gene_counts"],
                        genes,
                        phylum_colors,
                        genus_colors,
                        source_colors,
                        min_distinct_fraction=args.filtered_min_distinct_fraction,
                        completeness_available=completeness_available,
                    )
                )
            except ValueError as exc:
                if "No genomes with hits" not in str(exc):
                    raise
                filtered_skipped_groups.append(
                    {
                        "Group": slug,
                        "Subcategories": "; ".join(subcategories),
                    }
                )

        filtered_manifest_name = (
            f"itol_gene_counts_besthit__{category_slug}_split_by_subcategory_{filtered_suffix}.tsv"
        )
        base.write_rows(
            output_dir / filtered_manifest_name,
            filtered_manifest_rows,
            [
                "Group",
                "SubcategoryCount",
                "Subcategories",
                "GeneCount",
                "DistinctGeneCount",
                "MinDistinctGenesRequired",
                "GenomeCountWithAnyHits",
                "GenomeCount",
                "CutoffTag",
                "TreeFileName",
                "HeatmapFileName",
                "PhylumAnnotationFileName",
                "GenusAnnotationFileName",
                "SourceStripFileName",
                "GenomeNameLabelsFileName",
                "CompletenessPiechartFileName",
                "BundleFolder",
            ],
        )

        filtered_readme_lines = [
            f"Category: {args.category}",
            f"Distinct-gene filter fraction: {args.filtered_min_distinct_fraction}",
            (
                "Rule: keep a genome only when its distinct genes present in the "
                "bundle are >= ceil(distinct bundle gene count * fraction)."
            ),
            f"Bundle root: {filtered_bundle_root.name}",
            f"Retained filtered groups: {len(filtered_manifest_rows)}",
            f"Skipped filtered groups: {len(filtered_skipped_groups)}",
            "",
            "Filtered group summary:",
        ]
        for row in filtered_manifest_rows:
            filtered_readme_lines.append(
                f"- {row['Group']}: kept {row['GenomeCount']} of {row['GenomeCountWithAnyHits']} genomes, "
                f"cutoff {row['CutoffTag']}"
            )
            filtered_readme_lines.append(f"  Subcategories: {row['Subcategories']}")
        if filtered_skipped_groups:
            filtered_readme_lines.extend(["", "Skipped filtered groups (no genomes passed the cutoff):"])
            for row in filtered_skipped_groups:
                filtered_readme_lines.append(f"- {row['Group']}: {row['Subcategories']}")
        filtered_readme_name = f"README_{category_slug}_split_{filtered_suffix}.txt"
        (output_dir / filtered_readme_name).write_text(
            "\n".join(filtered_readme_lines) + "\n",
            encoding="utf-8",
        )
        filtered_bundle_root_name = filtered_bundle_root.name

    print(f"Output directory: {output_dir}")
    print(f"Manifest: {manifest_name}")
    for row in manifest_rows:
        print(row["BundleFolder"])
    if filtered_manifest_name:
        print(f"Filtered manifest: {filtered_manifest_name}")
        print(f"Filtered bundle root: {filtered_bundle_root_name}")
        print(f"Filtered README: {filtered_readme_name}")


if __name__ == "__main__":
    main()
