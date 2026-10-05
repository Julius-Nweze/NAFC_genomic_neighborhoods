#!/usr/bin/env python3
# Title          : rescope_bundle_shared_annotations.py
# Description    : Regenerate each bundle's shared iTOL annotation files to match its pruned tree
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/07
# Usage          : python3 rescope_bundle_shared_annotations.py [options]

"""Regenerate each bundle's shared iTOL annotation files to match its pruned tree."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

from Bio import Phylo

SUMMARIZE_SCRIPT = Path(__file__).resolve().parent.parent / "shared_dependencies" / "summarize_ref_genome_gene_counts_itol.py"
DEFAULT_METADATA_XLSX = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/NCBI_genome_info.xlsx"
)
DEFAULT_PHYLUM_LEGEND = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/itol_phylum_legend.tsv"
)
DEFAULT_SOURCE_LEGEND = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/itol_source_legend.tsv"
)
DEFAULT_BUNDLE_ROOTS = [
    Path(
        "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
        "gene_count_summary/category_bundles_corrected"
    ),
    Path(
        "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
        "gene_count_summary/category_bundles_min30pct_distinct_genes"
    ),
    Path(
        "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
        "gene_count_summary/category_bundles_min10pct_distinct_genes"
    ),
]


def load_legend(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        return {row[0]: row[1] for row in reader if len(row) >= 2}


def import_summarize_module():
    spec = importlib.util.spec_from_file_location("summarize_ref_genome_gene_counts_itol", SUMMARIZE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def find_bundle_tree(folder: Path) -> Path | None:
    candidates = sorted(folder.glob("Bacteria_71_fasttree*.nwk"))
    return candidates[0] if candidates else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-xlsx", default=str(DEFAULT_METADATA_XLSX))
    parser.add_argument("--phylum-legend", default=str(DEFAULT_PHYLUM_LEGEND))
    parser.add_argument("--source-legend", default=str(DEFAULT_SOURCE_LEGEND))
    parser.add_argument(
        "--bundle-root", action="append", dest="bundle_roots",
        help="A category_bundles_*/ directory to rescope. Repeatable. Defaults to all known bundle sets.",
    )
    args = parser.parse_args()

    bundle_roots = [Path(p) for p in args.bundle_roots] if args.bundle_roots else DEFAULT_BUNDLE_ROOTS

    mod = import_summarize_module()
    _header, metadata_lookup = mod.read_metadata(Path(args.metadata_xlsx))
    phylum_colors = load_legend(Path(args.phylum_legend))
    source_colors = load_legend(Path(args.source_legend))

    total_bundles = 0
    for bundle_root in bundle_roots:
        if not bundle_root.exists():
            print(f"SKIP (missing): {bundle_root}")
            continue
        for folder in sorted(p for p in bundle_root.iterdir() if p.is_dir()):
            tree_path = find_bundle_tree(folder)
            if tree_path is None:
                continue
            genomes = [t.name for t in Phylo.read(tree_path, "newick").get_terminals() if t.name]

            mod.write_itol_tree_colors(
                folder / "itol_phylum_colorstrip.txt",
                tree_path,
                metadata_lookup,
                "Phylum",
                phylum_colors,
            )
            mod.write_itol_colorstrip(
                folder / "itol_source_colorstrip.txt",
                "Genome source",
                genomes,
                metadata_lookup,
                "Source",
                strip_width=28,
                color_map=source_colors,
            )
            mod.write_itol_labels(
                folder / "itol_genome_name_labels.txt",
                genomes,
                metadata_lookup,
                "Genome name",
            )
            mod.write_itol_piechart(
                folder / "itol_completeness_piechart.txt",
                "Genome completeness (%)",
                genomes,
                metadata_lookup,
                "Completeness (%)",
            )
            print(f"  {bundle_root.name}/{folder.name}: rescoped 4 shared files to {len(genomes)} genomes")
            total_bundles += 1

    print(f"Rescoped shared annotation files in {total_bundles} bundle folders across {len(bundle_roots)} bundle sets")


if __name__ == "__main__":
    main()
