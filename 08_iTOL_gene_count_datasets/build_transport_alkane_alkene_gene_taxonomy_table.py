#!/usr/bin/env python3
# Title          : build_transport_alkane_alkene_gene_taxonomy_table.py
# Description    : Long-format table of Transportation, Alkanes and Alkenes gene hits per genome with taxonomy
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/17
# Usage          : python3 build_transport_alkane_alkene_gene_taxonomy_table.py

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUMMARY_SCRIPT = HERE.parent / "shared_dependencies" / "summarize_ref_genome_gene_counts_itol.py"
METADATA_XLSX = Path("/path/to/your/Ref_genomes/Analysis/NCBI_genome_info.xlsx")  # fallback
METADATA_XLSX_ALT = HERE.parents[1] / "NCBI_genome_info.xlsx"  # .../Result/NCBI_genome_info.xlsx

HEATMAP_DIR = HERE / "itol_gene_counts_besthit_by_category_corrected"
CATEGORY_FILES = {
    "Transportation": HEATMAP_DIR / "itol_gene_counts_besthit__transportation.txt",
    "Alkanes": HEATMAP_DIR / "itol_gene_counts_besthit__alkanes.txt",
    "Alkenes": HEATMAP_DIR / "itol_gene_counts_besthit__alkenes.txt",
}

OUTPUT_PATH = HERE / "transport_alkane_alkene_gene_hits_by_genome_taxonomy.tsv"


def load_summary_module():
    spec = importlib.util.spec_from_file_location("summarize_ref_genome_gene_counts_itol", SUMMARY_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def resolve_metadata_path() -> Path:
    for candidate in (METADATA_XLSX_ALT, METADATA_XLSX):
        if candidate.exists():
            return candidate
    raise FileNotFoundError("NCBI_genome_info.xlsx not found at either expected location")


def parse_heatmap(path: Path) -> tuple[list[str], list[tuple[str, list[int]]]]:
    genes: list[str] = []
    rows: list[tuple[str, list[int]]] = []
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
                rows.append((genome_id, counts))
    if not genes:
        raise ValueError(f"No FIELD_LABELS found in {path}")
    return genes, rows


def main() -> None:
    module = load_summary_module()
    metadata_path = resolve_metadata_path()
    _, lookup = module.read_metadata(metadata_path)

    out_rows = []
    missing_taxonomy: set[str] = set()

    for category, path in CATEGORY_FILES.items():
        genes, rows = parse_heatmap(path)
        for genome_id, counts in rows:
            meta = lookup.get(genome_id) or lookup.get(module.canonicalize_genome_id(genome_id))
            if meta is None:
                missing_taxonomy.add(genome_id)
                meta = {}
            for gene, count in zip(genes, counts):
                if count <= 0:
                    continue
                out_rows.append(
                    {
                        "Category": category,
                        "Gene": gene,
                        "Genome_ID": genome_id,
                        "Genome_name": meta.get("Genome name", ""),
                        "Phylum": meta.get("Phylum", ""),
                        "Class": meta.get("Class", ""),
                        "Order": meta.get("Order", ""),
                        "Family": meta.get("Family", ""),
                        "Genus": meta.get("Genus", ""),
                        "Species": meta.get("Species", ""),
                        "Source": meta.get("Source", ""),
                        "Copy_number": count,
                    }
                )

    out_rows.sort(key=lambda r: (r["Category"], r["Gene"], -r["Copy_number"], r["Genome_ID"]))

    with OUTPUT_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "Category",
                "Gene",
                "Genome_ID",
                "Genome_name",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "Source",
                "Copy_number",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows to {OUTPUT_PATH}")
    if missing_taxonomy:
        print(f"WARNING: {len(missing_taxonomy)} genome IDs had no taxonomy match: {sorted(missing_taxonomy)[:10]}...")


if __name__ == "__main__":
    main()
