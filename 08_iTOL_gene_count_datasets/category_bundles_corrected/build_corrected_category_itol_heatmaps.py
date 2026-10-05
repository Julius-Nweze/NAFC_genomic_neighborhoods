#!/usr/bin/env python3
# Title          : build_corrected_category_itol_heatmaps.py
# Description    : Build per-category iTOL heatmap datasets with one shared color scale
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/07
# Usage          : python3 build_corrected_category_itol_heatmaps.py [options]

"""Build per-category iTOL heatmap datasets with one shared color scale."""
from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from io import StringIO
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from Bio import Phylo

XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

DEFAULT_BESTHIT = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/besthit_assigned_hit_proteins.tsv"
)
DEFAULT_GENE_INFO = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/Gene_info.xlsx"
)
DEFAULT_OUTPUT_DIR = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/itol_gene_counts_besthit_by_category_corrected"
)
DEFAULT_BUNDLES_DIR = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/category_bundles_corrected"
)

# Category display order used for numbering plotting-ready bundle folders.
# Aromatics/Plastics use their split groups here (legibility, see docstring) in
# place of the single whole-category file, which still exists in --output-dir
# for analysis use (e.g. the Table S9 top-genome/core-gene computation).
CATEGORY_PLOT_ORDER = [
    "alkanes", "alkenes", "aromatics_1", "aromatics_2", "aromatics_3",
    "aromatic_carboxylate", "phenylacetate",
    "cyclohexanecarboxylate", "cyclohexylacetate", "naphthalene", "catechol_ortho",
    "catechol_meta", "protocatechuate", "oxalate", "benzoyl_coa", "tannin", "benzoate",
    "gallate", "gentisate", "pyrogallol", "hydroquinone", "hydroxyquinol",
    "plastics_1", "plastics_2", "plastics_3",
    "transportation", "beta_oxidation",
]

# Full, unpruned genome tree. Each bundle folder gets its own category-pruned
# copy of this (see write_pruned_bundle_tree) rather than a verbatim copy --
# otherwise a narrow category like Benzoyl-CoA (19 hit genomes out of 364)
# would still plot against all 364 leaves, most with a zero/blank row.
FULL_TREE_PATH = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/00_full_tree_shared/Bacteria_71_fasttree.nwk"
)

# Shared annotation files copied verbatim into every per-category bundle folder
# alongside that category's own heatmap and pruned tree, so each folder is a
# self-contained iTOL upload set. These are genome-keyed datasets, not tied to
# tree topology, so iTOL simply ignores rows without a matching pruned-tree leaf.
SHARED_BUNDLE_FILES = [
    Path("/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_phylum_colorstrip.txt"),
    Path("/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_source_colorstrip.txt"),
    Path("/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_genome_name_labels.txt"),
    Path("/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_completeness_piechart.txt"),
]

HEATMAP_COLOR = "#2166ac"
COLOR_MIN = "#f7fbff"
COLOR_MAX = "#08306b"
COLOR_NAN = "#d9d9d9"

# Fresh balanced split of the CURRENT Aromatics category (24 sub-categories,
# 187 genes), contiguous by sub-category name, ~3 roughly equal groups.
AROMATICS_SPLIT = {
    "aromatics_1": [
        "1-Methylnaphthalene/2-Methylnaphthalene", "2-Aminophenol", "2-Nitrobenzoate",
        "4-Hydroxyphenylacetate", "Aniline", "Anthranilate", "Benzene", "Biphenyl",
    ],
    "aromatics_2": [
        "Carbazol", "Chlorobenzen", "Chlorocatechol", "Dibenzofuran", "Fluoren",
        "Homogentisate", "Homoprotocatechuate", "Naphthalene", "Phenol", "Phthalate",
        "Pyrene", "Salicylate",
    ],
    "aromatics_3": ["Styrene", "Terephthalate", "Toluene", "p-Cumate"],
}

# Original subcategory grouping for Plastics, unaffected by the category
# restructuring (see module docstring) - reused as-is against current data.
PLASTICS_SPLIT = {
    "plastics_1": [
        "3PET/PCL/PLA", "LDPE", "NR", "Nylon", "O-PVA/PVA",
        "P3HP/P4HB/PEA/PES/PHA/PHB", "P3HV/PHA/PHBV", "PBAT",
        "PBAT/PBS/PBSA/PCL/PET/PHA/PHB/PLA", "PBS", "PBS/PBSA", "PBS/PBSA/PCL",
        "PBS/PBSA/PCL/PES/PHA/PHB/PLA", "PBS/PBSA/PCL/PES/PHA/PLA",
        "PBS/PCL/PET/PHA/PHB/PLA/PU", "PBS/PCL/PHA/PHB/PLA", "PBSA",
        "PBSA/PCL/PES/PHA/PHBV/PLA", "PBSA/PCL/PLA", "PBSA/PLA",
    ],
    "plastics_2": [
        "PCL", "PCL/PES/PET", "PCL/PET", "PCL/PET/PU", "PCL/PHA/PHBV/PHPV",
        "PCL/PU", "PE", "PEG", "PET", "PET/PLA", "PET/PU",
    ],
    "plastics_3": [
        "PHA", "PHA/PHB", "PHA/PHB/PHBV", "PHA/PHB/PPL", "PHA/PHO", "PHB", "PLA",
        "PMCL", "PU", "PVA",
    ],
}


def load_gene_subcategories(path: Path, category: str) -> dict[str, list[str]]:
    """Return {gene: [subcategory, ...]} restricted to one Category, using the
    same raw-XML reader as load_gene_categories (see that function for why)."""
    with ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", XLSX_NS):
                shared_strings.append(
                    "".join(node.text or "" for node in si.iterfind(".//a:t", XLSX_NS))
                )
        sheet_root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        sheet_data = sheet_root.find("a:sheetData", XLSX_NS)
        rows: list[list[str]] = []
        for row in sheet_data.findall("a:row", XLSX_NS):
            values: list[str] = []
            for cell in row.findall("a:c", XLSX_NS):
                cell_type = cell.get("t")
                if cell_type == "inlineStr":
                    is_el = cell.find("a:is", XLSX_NS)
                    text = "".join(n.text or "" for n in is_el.iterfind(".//a:t", XLSX_NS)) if is_el is not None else ""
                elif cell_type == "s":
                    v_el = cell.find("a:v", XLSX_NS)
                    text = shared_strings[int(v_el.text)] if v_el is not None else ""
                else:
                    v_el = cell.find("a:v", XLSX_NS)
                    text = v_el.text if v_el is not None else ""
                values.append(text or "")
            rows.append(values)

    header = rows[0]
    cat_idx = header.index("Category")
    sub_idx = header.index("Sub category")
    gene_idx = header.index("Gene")
    result: dict[str, list[str]] = defaultdict(list)
    for row in rows[1:]:
        if len(row) <= max(cat_idx, sub_idx, gene_idx):
            continue
        if row[cat_idx].strip() != category:
            continue
        gene, sub = row[gene_idx].strip(), row[sub_idx].strip()
        if gene and sub:
            result[gene].append(sub)
    return dict(result)


def load_gene_categories(path: Path) -> dict[str, str]:
    """Read the Gene_info.xlsx 'Code' sheet's Category/Gene columns directly
    from the raw sheet XML (handles both shared-string and inline-string
    cells, both of which occur in Gene_info.xlsx)."""
    with ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", XLSX_NS):
                shared_strings.append(
                    "".join(node.text or "" for node in si.iterfind(".//a:t", XLSX_NS))
                )

        # "Code" sheet is expected to be the first sheet in this workbook.
        sheet_root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        sheet_data = sheet_root.find("a:sheetData", XLSX_NS)
        rows: list[list[str]] = []
        for row in sheet_data.findall("a:row", XLSX_NS):
            values: list[str] = []
            for cell in row.findall("a:c", XLSX_NS):
                cell_type = cell.get("t")
                if cell_type == "inlineStr":
                    is_el = cell.find("a:is", XLSX_NS)
                    text = "".join(n.text or "" for n in is_el.iterfind(".//a:t", XLSX_NS)) if is_el is not None else ""
                elif cell_type == "s":
                    v_el = cell.find("a:v", XLSX_NS)
                    text = shared_strings[int(v_el.text)] if v_el is not None else ""
                else:
                    v_el = cell.find("a:v", XLSX_NS)
                    text = v_el.text if v_el is not None else ""
                values.append(text or "")
            rows.append(values)

    header = rows[0]
    cat_idx = header.index("Category")
    gene_idx = header.index("Gene")
    gene_cat: dict[str, str] = {}
    for row in rows[1:]:
        if len(row) <= max(cat_idx, gene_idx):
            continue
        gene, cat = row[gene_idx].strip(), row[cat_idx].strip()
        if gene and cat:
            gene_cat[gene] = cat
    return gene_cat


def load_genome_gene_counts(path: Path) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["genome_id"]][row["chosen_gene"]] += 1
    return counts


def slug(category: str) -> str:
    return category.lower().replace(" ", "_").replace("-", "_")


def write_category_heatmap(
    out_path: Path,
    category: str,
    genes: list[str],
    genome_gene_counts: dict[str, dict[str, int]],
    global_min: int,
    global_max: int,
) -> int:
    genes = sorted(genes)
    lines = [
        "DATASET_HEATMAP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\tgene_counts_besthit__{category}",
        f"COLOR\t{HEATMAP_COLOR}",
        "FIELD_LABELS\t" + "\t".join(genes),
        f"COLOR_MIN\t{COLOR_MIN}",
        f"COLOR_MAX\t{COLOR_MAX}",
        f"COLOR_NAN\t{COLOR_NAN}",
        "AUTO_LEGEND\t1",
        "SHOW_INTERNAL\t0",
        "SHOW_LABELS\t1",
        "LABEL_ROTATION\t90",
        "STRIP_WIDTH\t18",
        "MARGIN\t20",
        "BORDER_WIDTH\t0",
        f"USER_MIN_VALUE\t{global_min}",
        f"USER_MAX_VALUE\t{global_max}",
        "DATA",
    ]
    n_rows = 0
    for genome in sorted(genome_gene_counts):
        row_counts = [genome_gene_counts[genome].get(g, 0) for g in genes]
        if any(row_counts):
            lines.append(genome + "\t" + "\t".join(str(v) for v in row_counts))
            n_rows += 1
    out_path.write_text("\n".join(lines) + "\n")
    return n_rows


def write_pruned_bundle_tree(
    out_path: Path,
    full_tree_text: str,
    genes: list[str],
    genome_gene_counts: dict[str, dict[str, int]],
) -> int:
    """Prune the full tree down to genomes with at least one hit among `genes`
    (the same genome set that ends up as a row in that category's heatmap)."""
    hit_genomes = {
        genome
        for genome, gene_counts in genome_gene_counts.items()
        if any(gene_counts.get(g, 0) for g in genes)
    }
    tree = Phylo.read(StringIO(full_tree_text), "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in hit_genomes:
            tree.prune(target=terminal)
    Phylo.write(tree, out_path, "newick")
    return len(hit_genomes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--besthit-table", default=str(DEFAULT_BESTHIT))
    parser.add_argument("--gene-info-xlsx", default=str(DEFAULT_GENE_INFO))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument(
        "--bundles-dir",
        default=str(DEFAULT_BUNDLES_DIR),
        help="Where to write one self-contained plotting folder per category "
        "(heatmap + tree + shared annotation files). Pass an empty string to skip.",
    )
    parser.add_argument(
        "--global-min-value",
        type=int,
        default=1,
        help="USER_MIN_VALUE applied to every category (default 1).",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_cat = load_gene_categories(Path(args.gene_info_xlsx))
    genome_gene_counts = load_genome_gene_counts(Path(args.besthit_table))

    category_genes: dict[str, list[str]] = defaultdict(list)
    for gene, cat in gene_cat.items():
        category_genes[cat].append(gene)

    global_max = max(
        count
        for gene_counts in genome_gene_counts.values()
        for count in gene_counts.values()
    )

    print(f"Global scale applied to all categories: {args.global_min_value}-{global_max}")
    written: dict[str, Path] = {}
    written_genes: dict[str, list[str]] = {}
    for category, genes in sorted(category_genes.items()):
        out_path = out_dir / f"itol_gene_counts_besthit__{slug(category)}.txt"
        n_rows = write_category_heatmap(
            out_path, category, genes, genome_gene_counts, args.global_min_value, global_max
        )
        written[slug(category)] = out_path
        written_genes[slug(category)] = genes
        print(f"{category}: {n_rows} genome rows -> {out_path}")

    # Legibility splits for the two largest categories (see module docstring).
    for parent_category, split in (("Aromatics", AROMATICS_SPLIT), ("Plastics", PLASTICS_SPLIT)):
        gene_subs = load_gene_subcategories(Path(args.gene_info_xlsx), parent_category)
        for group_name, subcats in split.items():
            subcats_set = set(subcats)
            genes = sorted({g for g, subs in gene_subs.items() if subcats_set & set(subs)})
            out_path = out_dir / f"itol_gene_counts_besthit__{group_name}.txt"
            n_rows = write_category_heatmap(
                out_path, group_name, genes, genome_gene_counts, args.global_min_value, global_max
            )
            written[group_name] = out_path
            written_genes[group_name] = genes
            print(f"{group_name} ({len(genes)} genes, subset of {parent_category}): {n_rows} genome rows -> {out_path}")

    if args.bundles_dir:
        bundles_dir = Path(args.bundles_dir)
        bundles_dir.mkdir(parents=True, exist_ok=True)
        # Whole-category Aromatics/Plastics are excluded from plotting bundles --
        # their split groups (already in CATEGORY_PLOT_ORDER) supersede them for
        # legibility. The whole-category files remain in --output-dir for analysis.
        superseded = {"aromatics", "plastics"}
        ordered = [s for s in CATEGORY_PLOT_ORDER if s in written] + [
            s for s in written if s not in CATEGORY_PLOT_ORDER and s not in superseded
        ]
        full_tree_text = FULL_TREE_PATH.read_text(encoding="utf-8")
        for i, cat_slug in enumerate(ordered, start=1):
            folder = bundles_dir / f"{i:02d}_{cat_slug}"
            folder.mkdir(parents=True, exist_ok=True)
            shutil.copy2(written[cat_slug], folder / written[cat_slug].name)
            for shared_path in SHARED_BUNDLE_FILES:
                if shared_path.exists():
                    shutil.copy2(shared_path, folder / shared_path.name)
            # Tree file name states the category (see FULL_TREE_PATH docstring note) so
            # each bundle's tree is distinguishable once uploaded to iTOL - a bare
            # "Bacteria_71_fasttree.nwk" repeated across all 26 bundles left no way to
            # tell which uploaded tree belonged to which category.
            tree_name = f"Bacteria_71_fasttree__{cat_slug}.nwk"
            n_leaves = write_pruned_bundle_tree(
                folder / tree_name,
                full_tree_text,
                written_genes[cat_slug],
                genome_gene_counts,
            )
            print(f"  {cat_slug}: pruned tree to {n_leaves} hit genomes -> {folder / tree_name}")
        print(f"Wrote {len(ordered)} self-contained plotting bundles -> {bundles_dir}")


if __name__ == "__main__":
    main()
