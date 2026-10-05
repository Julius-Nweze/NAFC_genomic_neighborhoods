#!/usr/bin/env python3
# Title          : build_transport_alkane_alkene_genome_hit_summary.py
# Description    : Per-genome hit-count summary for the Transportation, Alkanes and Alkenes categories
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/17
# Usage          : python3 build_transport_alkane_alkene_genome_hit_summary.py

from __future__ import annotations

import csv
import importlib.util
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUMMARY_SCRIPT = HERE.parent / "shared_dependencies" / "summarize_ref_genome_gene_counts_itol.py"
ANALYSIS_ROOT = Path("/path/to/your/Ref_genomes/Analysis")
METADATA_XLSX_CANDIDATES = [
    HERE.parents[1] / "NCBI_genome_info.xlsx",  # .../Result/NCBI_genome_info.xlsx
    ANALYSIS_ROOT / "NCBI_genome_info.xlsx",
]

TREE_PATH = HERE / "00_full_tree_shared" / "NCBI_MAGs_Colla_Bacteria_71_fasttree.nwk"

HEATMAP_DIR = HERE / "itol_gene_counts_besthit_by_category_corrected"
CATEGORY_FILES = {
    "Transportation": HEATMAP_DIR / "itol_gene_counts_besthit__transportation.txt",
    "Alkanes": HEATMAP_DIR / "itol_gene_counts_besthit__alkanes.txt",
    "Alkenes": HEATMAP_DIR / "itol_gene_counts_besthit__alkenes.txt",
}

OUTPUT_PATH = HERE / "transport_alkane_alkene_genome_hit_summary_all_364.tsv"


def load_summary_module():
    spec = importlib.util.spec_from_file_location("summarize_ref_genome_gene_counts_itol", SUMMARY_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def resolve_metadata_path() -> Path:
    for candidate in METADATA_XLSX_CANDIDATES:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("NCBI_genome_info.xlsx not found at any expected location")


def load_all_tree_tips(path: Path) -> list[str]:
    text = path.read_text()
    return re.findall(r"[\(,]([A-Za-z0-9_\.\-]+):", text)


def parse_heatmap(path: Path) -> tuple[list[str], dict[str, list[int]]]:
    genes: list[str] = []
    rows: dict[str, list[int]] = {}
    in_data = False
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("FIELD_LABELS"):
                genes = line.split("\t")[1:]
                continue
            if line == "DATA":
                in_data = True
                continue
            if in_data:
                parts = line.split("\t")
                genome_id = parts[0]
                counts = [int(x) if x not in ("", "-") else 0 for x in parts[1:]]
                rows[genome_id] = counts
    if not genes:
        raise ValueError(f"No FIELD_LABELS found in {path}")
    return genes, rows


def main() -> None:
    module = load_summary_module()
    metadata_path = resolve_metadata_path()
    _, lookup = module.read_metadata(metadata_path)

    all_genomes = load_all_tree_tips(TREE_PATH)
    assert len(all_genomes) == 364, f"expected 364 tree tips, found {len(all_genomes)}"

    out_rows = []
    missing_taxonomy: set[str] = set()

    for category, path in CATEGORY_FILES.items():
        genes, hit_rows = parse_heatmap(path)
        panel_size = len(genes)
        for genome_id in all_genomes:
            counts = hit_rows.get(genome_id, [0] * panel_size)
            distinct_genes = sum(1 for c in counts if c > 0)
            total_copies = sum(counts)
            genes_present = ";".join(g for g, c in zip(genes, counts) if c > 0)

            meta = lookup.get(genome_id) or lookup.get(module.canonicalize_genome_id(genome_id))
            if meta is None:
                missing_taxonomy.add(genome_id)
                meta = {}

            out_rows.append(
                {
                    "Category": category,
                    "Genome_ID": genome_id,
                    "Genome_name": meta.get("Genome name", ""),
                    "Phylum": meta.get("Phylum", ""),
                    "Class": meta.get("Class", ""),
                    "Order": meta.get("Order", ""),
                    "Family": meta.get("Family", ""),
                    "Genus": meta.get("Genus", ""),
                    "Species": meta.get("Species", ""),
                    "Source": meta.get("Source", ""),
                    "Distinct_genes_present": distinct_genes,
                    "Genes_in_panel": panel_size,
                    "Pct_distinct_genes": round(100 * distinct_genes / panel_size, 1),
                    "Total_copy_number": total_copies,
                    "Genes_present": genes_present,
                }
            )

    out_rows.sort(key=lambda r: (r["Category"], -r["Distinct_genes_present"], -r["Total_copy_number"], r["Genome_ID"]))

    with OUTPUT_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "Category",
                "Genome_ID",
                "Genome_name",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "Source",
                "Distinct_genes_present",
                "Genes_in_panel",
                "Pct_distinct_genes",
                "Total_copy_number",
                "Genes_present",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows ({len(all_genomes)} genomes x {len(CATEGORY_FILES)} categories) to {OUTPUT_PATH}")
    if missing_taxonomy:
        print(f"WARNING: {len(missing_taxonomy)} genome IDs had no taxonomy match: {sorted(missing_taxonomy)[:10]}...")


if __name__ == "__main__":
    main()
