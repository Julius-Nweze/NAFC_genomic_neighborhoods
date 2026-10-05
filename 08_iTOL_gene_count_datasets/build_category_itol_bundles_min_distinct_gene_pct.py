#!/usr/bin/env python3
# Title          : build_category_itol_bundles_min_distinct_gene_pct.py
# Description    : Build per-category iTOL bundles keeping genomes above a minimum fraction of distinct category genes
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/11
# Usage          : python3 build_category_itol_bundles_min_distinct_gene_pct.py [options]

"""Build per-category iTOL bundles keeping genomes above a minimum fraction of distinct category genes."""
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
FULL_TREE_PATH = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/00_full_tree_shared/Bacteria_71_fasttree.nwk"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary"
)

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

# Same category display order / subcategory splits as
# category_bundles_corrected/build_corrected_category_itol_heatmaps.py, kept in
# sync deliberately so the two bundle sets cover the identical 26 groups.
CATEGORY_PLOT_ORDER = [
    "alkanes", "alkenes", "aromatics_1", "aromatics_2", "aromatics_3",
    "aromatic_carboxylate", "phenylacetate",
    "cyclohexanecarboxylate", "cyclohexylacetate", "naphthalene", "catechol_ortho",
    "catechol_meta", "protocatechuate", "oxalate", "benzoyl_coa", "tannin", "benzoate",
    "gallate", "gentisate", "pyrogallol", "hydroquinone", "hydroxyquinol",
    "plastics_1", "plastics_2", "plastics_3",
    "transportation", "beta_oxidation",
]

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


def _read_xlsx_rows(path: Path) -> list[list[str]]:
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
    return rows


def load_gene_categories(path: Path) -> dict[str, str]:
    rows = _read_xlsx_rows(path)
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


def load_gene_subcategories(path: Path, category: str) -> dict[str, list[str]]:
    rows = _read_xlsx_rows(path)
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


def load_genome_gene_counts(path: Path) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["genome_id"]][row["chosen_gene"]] += 1
    return counts


def slug(category: str) -> str:
    return category.lower().replace(" ", "_").replace("-", "_")


def kept_genomes(
    genes: list[str],
    genome_gene_counts: dict[str, dict[str, int]],
    min_fraction: float,
) -> dict[str, int]:
    """Return {genome: distinct_genes_present} for genomes clearing the
    distinct-gene-fraction cutoff (>=1 hit counts as present regardless of
    copy number)."""
    n_genes = len(genes)
    kept: dict[str, int] = {}
    for genome, counts in genome_gene_counts.items():
        distinct_present = sum(1 for g in genes if counts.get(g, 0) > 0)
        if distinct_present == 0:
            continue
        if distinct_present / n_genes >= min_fraction - 1e-9:
            kept[genome] = distinct_present
    return kept


def cap_to_top_genomes(
    kept: dict[str, int],
    genes: list[str],
    genome_gene_counts: dict[str, dict[str, int]],
    max_genomes: int,
) -> dict[str, int]:
    """If more than max_genomes genomes clear the fraction cutoff, keep only the
    max_genomes with the highest distinct-gene count (readability cap for
    trees that get merged into composite figures). Ties broken by
    total gene-copy count, then genome ID, for determinism."""
    if len(kept) <= max_genomes:
        return kept

    def rank_key(genome: str) -> tuple[int, int, str]:
        total_copies = sum(genome_gene_counts[genome].get(g, 0) for g in genes)
        return (-kept[genome], -total_copies, genome)

    top_genomes = sorted(kept, key=rank_key)[:max_genomes]
    return {genome: kept[genome] for genome in top_genomes}


def order_genes_by_subcategory(
    gene_subs: dict[str, list[str]],
    subcats_order: list[str],
) -> list[str]:
    """Order genes by their subcategory's position in subcats_order (so iTOL heatmap
    columns group genes by subcategory instead of a flat alphabetical gene sort),
    breaking ties alphabetically by gene name within a subcategory. A gene tagged
    with multiple subcategories is placed at the earliest (first-listed) one."""
    subcat_rank = {name: i for i, name in enumerate(subcats_order)}
    subcats_set = set(subcats_order)
    genes = {g for g, subs in gene_subs.items() if subcats_set & set(subs)}

    def sort_key(gene: str) -> tuple[int, str]:
        ranks = [subcat_rank[s] for s in gene_subs[gene] if s in subcat_rank]
        return (min(ranks) if ranks else len(subcats_order), gene)

    return sorted(genes, key=sort_key)


def write_heatmap(
    out_path: Path,
    category_label: str,
    genes: list[str],
    genome_gene_counts: dict[str, dict[str, int]],
    kept: dict[str, int],
    global_min: int,
    global_max: int,
    field_order: list[str] | None = None,
) -> None:
    genes = field_order if field_order is not None else sorted(genes)
    lines = [
        "DATASET_HEATMAP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\tgene_counts_besthit__{category_label}",
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
    for genome in sorted(kept):
        row_counts = [genome_gene_counts[genome].get(g, 0) for g in genes]
        lines.append(genome + "\t" + "\t".join(str(v) for v in row_counts))
    out_path.write_text("\n".join(lines) + "\n")


def write_pruned_tree(out_path: Path, full_tree_text: str, keep_names: set[str]) -> None:
    tree = Phylo.read(StringIO(full_tree_text), "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep_names:
            tree.prune(target=terminal)
    Phylo.write(tree, out_path, "newick")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--besthit-table", default=str(DEFAULT_BESTHIT))
    parser.add_argument("--gene-info-xlsx", default=str(DEFAULT_GENE_INFO))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument(
        "--min-fraction", type=float, default=0.30,
        help="Minimum distinct-genes-present / genes-in-category fraction required to keep a genome (default 0.30).",
    )
    parser.add_argument(
        "--global-min-value", type=int, default=1,
        help="USER_MIN_VALUE applied to every category heatmap (default 1).",
    )
    parser.add_argument(
        "--category-fraction", action="append", default=[],
        help="Per-category cutoff override as slug=fraction, e.g. alkanes=0.15. Repeatable. "
        "Overrides --min-fraction for just that category; folder/file names reflect its own cutoff.",
    )
    parser.add_argument(
        "--only", default=None,
        help="Comma-separated category slugs to (re)generate, leaving every other existing bundle folder "
        "untouched. Defaults to all categories in CATEGORY_PLOT_ORDER.",
    )
    parser.add_argument(
        "--max-genomes", type=int, default=None,
        help="If a category's cutoff-passing genome count exceeds this, keep only the top N by distinct-gene "
        "count (readability cap for trees merged into composite figures, e.g. Figure 1a-d). Default: no cap.",
    )
    args = parser.parse_args()

    fraction_overrides: dict[str, float] = {}
    for item in args.category_fraction:
        cat_slug_part, _, frac_part = item.partition("=")
        fraction_overrides[cat_slug_part.strip()] = float(frac_part)

    only_slugs = {s.strip() for s in args.only.split(",")} if args.only else None

    pct_label = f"min{round(args.min_fraction * 100):g}pct"
    bundles_dir = Path(args.output_root) / f"category_bundles_{pct_label}_distinct_genes"
    bundles_dir.mkdir(parents=True, exist_ok=True)

    gene_cat = load_gene_categories(Path(args.gene_info_xlsx))
    genome_gene_counts = load_genome_gene_counts(Path(args.besthit_table))

    category_genes: dict[str, list[str]] = defaultdict(list)
    for gene, cat in gene_cat.items():
        category_genes[cat].append(gene)

    slug_genes: dict[str, list[str]] = {slug(cat): genes for cat, genes in category_genes.items()}
    field_order_by_group: dict[str, list[str]] = {}

    for parent_category, split in (("Aromatics", AROMATICS_SPLIT), ("Plastics", PLASTICS_SPLIT)):
        gene_subs = load_gene_subcategories(Path(args.gene_info_xlsx), parent_category)
        for group_name, subcats in split.items():
            subcats_set = set(subcats)
            genes = sorted({g for g, subs in gene_subs.items() if subcats_set & set(subs)})
            slug_genes[group_name] = genes
            # iTOL heatmap columns grouped by subcategory for Aromatics and Plastics, instead of
            # the flat alphabetical gene order used for every other category.
            field_order_by_group[group_name] = order_genes_by_subcategory(gene_subs, subcats)

    global_max = max(
        count
        for gene_counts in genome_gene_counts.values()
        for count in gene_counts.values()
    )

    full_tree_text = FULL_TREE_PATH.read_text(encoding="utf-8")

    print(f"Default cutoff: keep genomes with >= {args.min_fraction:.0%} of a category's distinct genes present")
    if fraction_overrides:
        print(f"Per-category overrides: {fraction_overrides}")
    if only_slugs:
        print(f"Only (re)generating: {sorted(only_slugs)} - every other existing bundle folder left untouched")
    if args.max_genomes:
        print(f"Readability cap: at most {args.max_genomes} genomes per category, ranked by distinct-gene count")
    print(f"Global heatmap scale: {args.global_min_value}-{global_max}")

    manifest_path = bundles_dir / f"README_{pct_label}_distinct_genes_cutoff.tsv"
    manifest_by_category: dict[str, dict[str, str]] = {}
    if manifest_path.exists():
        with manifest_path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row.setdefault("cutoff_fraction", str(args.min_fraction))
                row.setdefault("genomes_before_cap", row.get("genomes_kept", ""))
                row.setdefault("capped", "no")
                manifest_by_category[row["category"]] = row

    for i, cat_slug in enumerate(CATEGORY_PLOT_ORDER, start=1):
        if only_slugs is not None and cat_slug not in only_slugs:
            continue
        genes = slug_genes.get(cat_slug)
        if not genes:
            print(f"  SKIP {cat_slug}: no genes found")
            continue
        fraction = fraction_overrides.get(cat_slug, args.min_fraction)
        cat_pct_label = f"min{round(fraction * 100):g}pct"
        kept = kept_genomes(genes, genome_gene_counts, fraction)
        min_required = __import__("math").ceil(len(genes) * fraction - 1e-9)
        n_before_cap = len(kept)
        if args.max_genomes:
            kept = cap_to_top_genomes(kept, genes, genome_gene_counts, args.max_genomes)
        capped = args.max_genomes is not None and n_before_cap > len(kept)

        # If this category was previously built at a different cutoff (e.g. re-running
        # with --only after a --category-fraction override), remove its old folder so a
        # stale, differently-cut bundle doesn't linger alongside the new one.
        for stale in bundles_dir.glob(f"{i:02d}_{cat_slug}_min*pct"):
            if stale.name != f"{i:02d}_{cat_slug}_{cat_pct_label}":
                shutil.rmtree(stale)
                print(f"  removed stale folder from a previous cutoff: {stale}")

        if not kept:
            print(
                f"  SKIP {cat_slug}: 0 genomes have >= {min_required} of {len(genes)} "
                f"distinct genes (>= {fraction:.0%}) - no tree/heatmap to build"
            )
            manifest_by_category[cat_slug] = {
                "category": cat_slug,
                "genes_in_panel": str(len(genes)),
                "cutoff_fraction": str(fraction),
                "min_distinct_genes_required": str(min_required),
                "genomes_kept": "0",
                "folder": "(skipped - no genome meets cutoff)",
            }
            continue

        folder = bundles_dir / f"{i:02d}_{cat_slug}_{cat_pct_label}"
        folder.mkdir(parents=True, exist_ok=True)

        heatmap_path = folder / f"itol_gene_counts_besthit__{cat_slug}_{cat_pct_label}.txt"
        write_heatmap(
            heatmap_path, f"{cat_slug}_{cat_pct_label}", genes, genome_gene_counts, kept,
            args.global_min_value, global_max,
            field_order=field_order_by_group.get(cat_slug),
        )

        # Tree file name states both the category and the cutoff so each uploaded tree is identifiable.
        tree_path = folder / f"Bacteria_71_fasttree__{cat_slug}_{cat_pct_label}.nwk"
        write_pruned_tree(tree_path, full_tree_text, set(kept))

        for shared_path in SHARED_BUNDLE_FILES:
            if shared_path.exists():
                shutil.copy2(shared_path, folder / shared_path.name)

        cap_note = f", capped from {n_before_cap} (top {len(kept)} by gene count)" if capped else ""
        print(
            f"  {cat_slug}: {len(genes)} genes in panel, "
            f"{len(kept)} genomes kept (>= {min_required} genes, {fraction:.0%}{cap_note}) -> {folder}"
        )
        manifest_by_category[cat_slug] = {
            "category": cat_slug,
            "genes_in_panel": str(len(genes)),
            "cutoff_fraction": str(fraction),
            "min_distinct_genes_required": str(min_required),
            "genomes_before_cap": str(n_before_cap),
            "genomes_kept": str(len(kept)),
            "capped": "yes" if capped else "no",
            "folder": str(folder),
        }

    manifest_rows = [manifest_by_category[c] for c in CATEGORY_PLOT_ORDER if c in manifest_by_category]
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["category", "genes_in_panel", "cutoff_fraction", "min_distinct_genes_required",
                        "genomes_before_cap", "genomes_kept", "capped", "folder"],
            delimiter="\t",
            extrasaction="ignore",
            restval="",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Wrote/updated {len(manifest_rows)} rows in manifest -> {manifest_path}")


if __name__ == "__main__":
    main()
