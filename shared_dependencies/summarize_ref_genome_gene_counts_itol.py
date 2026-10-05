#!/usr/bin/env python3
# Title          : summarize_ref_genome_gene_counts_itol.py
# Description    : Per-genome, per-category best-hit gene counts and iTOL datasets (called by the master pipeline)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/20
# Usage          : python3 summarize_ref_genome_gene_counts_itol.py [options]

from __future__ import annotations

import argparse
import colorsys
import csv
from io import StringIO
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable
from zipfile import ZipFile
from xml.etree import ElementTree as ET


XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize per-genome gene copy numbers from KO BLASTP output, "
            "compute averages at the gene/category/sub-category levels, "
            "and generate iTOL heatmap datasets aligned to an existing tree."
        )
    )
    parser.add_argument(
        "--gene-info-xlsx",
        default="/path/to/your/Ref_genomes/Result/Gene_info.xlsx",
        help="Gene annotation workbook.",
    )
    parser.add_argument(
        "--ko-output",
        default="/path/to/your/Ref_genomes/Result/KO_proteins.output.txt",
        help="Raw two-column KO BLASTP output.",
    )
    parser.add_argument(
        "--ko-scored-output",
        default="/path/to/your/Ref_genomes/Result/KO_proteins.output.scored.tsv",
        help=(
            "BLASTP output with score columns used for best-hit assignment. "
            "Expected columns: qseqid sseqid pident qcovs length evalue bitscore"
        ),
    )
    parser.add_argument(
        "--metadata-xlsx",
        default="/path/to/your/Ref_genomes/Result/NCBI_genome_info.xlsx",
        help="Genome metadata workbook.",
    )
    parser.add_argument(
        "--tree",
        default="/path/to/your/Ref_genomes/Result/ITOL/Bacteria_71_fasttree.nwk",
        help="Newick tree whose tip labels define the desired genome order for iTOL.",
    )
    parser.add_argument(
        "--output-dir",
        default="/tmp/ref_genome_gene_count_itol",
        help="Directory for summary tables and iTOL datasets.",
    )
    parser.add_argument(
        "--min-pident",
        type=float,
        default=30.0,
        help=(
            "Minimum percent identity (0-100) required for a scored BLASTP row to "
            "count as a gene hit, applied in addition to the e-value cutoff already "
            "baked into the BLASTP run. Rows below this are dropped before best-hit "
            "resolution. Set to 0 to disable."
        ),
    )
    parser.add_argument(
        "--min-qcovs",
        type=float,
        default=50.0,
        help=(
            "Minimum query coverage percent (0-100) required for a scored BLASTP row "
            "to count as a gene hit, applied alongside --min-pident. Set to 0 to disable."
        ),
    )
    parser.add_argument(
        "--filtered-category-min-distinct-fraction",
        type=float,
        default=0.0,
        help=(
            "Optional fraction (0-1) used to generate an extra set of simplified "
            "per-category iTOL bundles that keep only genomes with enough distinct "
            "genes in that category. The threshold is ceil(distinct_category_gene_count "
            "* fraction), so 0.3 keeps genomes with at least 3 distinct genes in a "
            "10-gene category. Set to 0 to disable this extra output."
        ),
    )
    return parser.parse_args()


def normalize_text(value: str) -> str:
    return value.strip() if isinstance(value, str) else value


def canonicalize_genome_id(genome_id: str) -> str:
    genome_id = normalize_text(genome_id)

    if genome_id.endswith(".medaka"):
        stem = genome_id[: -len(".medaka")].replace(".", "_")
        if stem and stem[0].isdigit():
            return f"g_{stem}_medaka"
        return f"{stem}_medaka"

    if re.match(r"^(GCF|GCA)_", genome_id):
        return genome_id.rsplit(".", 1)[0]

    if re.match(r"^[A-Z][A-Z0-9_]*\.\d+$", genome_id):
        return genome_id.rsplit(".", 1)[0]

    return genome_id


def col_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    total = 0
    for ch in letters:
        total = total * 26 + (ord(ch.upper()) - 64)
    return total - 1


def trim_trailing_empty(values: list[str]) -> list[str]:
    while values and values[-1] == "":
        values.pop()
    return values


def load_xlsx_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            shared_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in shared_root.findall("a:si", XLSX_NS):
                shared_strings.append(
                    "".join(node.text or "" for node in si.iterfind(".//a:t", XLSX_NS))
                )

        sheet_root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        sheet_data = sheet_root.find("a:sheetData", XLSX_NS)
        if sheet_data is None:
            raise ValueError(f"sheetData not found in {path}")

        parsed_rows: list[dict[int, str]] = []
        max_index = 0
        for row in sheet_data.findall("a:row", XLSX_NS):
            parsed: dict[int, str] = {}
            for cell in row.findall("a:c", XLSX_NS):
                idx = col_index(cell.attrib.get("r", "A1"))
                cell_type = cell.attrib.get("t")
                if cell_type == "inlineStr":
                    is_node = cell.find("a:is", XLSX_NS)
                    value = "".join(t.text or "" for t in is_node.iterfind(".//a:t", XLSX_NS)) if is_node is not None else ""
                else:
                    value_node = cell.find("a:v", XLSX_NS)
                    if value_node is None:
                        value = ""
                    else:
                        raw = value_node.text or ""
                        if cell_type == "s":
                            value = shared_strings[int(raw)]
                        else:
                            value = raw
                parsed[idx] = normalize_text(value)
                max_index = max(max_index, idx)
            parsed_rows.append(parsed)

    table: list[list[str]] = []
    for parsed in parsed_rows:
        row = [parsed.get(i, "") for i in range(max_index + 1)]
        table.append(trim_trailing_empty(row))

    if not table:
        raise ValueError(f"{path} is empty")

    header = table[0]
    rows = []
    for row in table[1:]:
        padded = row + [""] * max(0, len(header) - len(row))
        rows.append({header[i]: padded[i] for i in range(len(header))})
    return header, rows


def read_tree_labels(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8").strip()
    labels = re.findall(r"(?<=[(,])([^():,]+):", text)
    if not labels:
        raise ValueError(f"No tip labels found in {path}")
    if len(labels) != len(set(labels)):
        raise ValueError(f"Duplicate tip labels found in {path}")
    return labels


def read_gene_info(path: Path):
    header, rows = load_xlsx_rows(path)
    required = ["Category", "Sub category", "Gene", "Gene name", "KO", "UniProt", "EC"]
    missing = [field for field in required if field not in header]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")

    unique_gene_order: list[str] = []
    seen_gene = set()
    gene_meta_sets: dict[str, set[tuple[str, str, str, str, str]]] = defaultdict(set)
    query_to_rows_by_uniprot: dict[str, list[dict[str, str]]] = defaultdict(list)
    query_to_rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    category_gene_row_order: dict[str, list[str]] = defaultdict(list)
    subcategory_gene_row_order: dict[str, list[str]] = defaultdict(list)
    seen_category_subcategory_gene: set[tuple[str, str, str]] = set()
    seen_subcategory_gene: set[tuple[str, str]] = set()

    for row in rows:
        clean = {key: normalize_text(value) for key, value in row.items()}
        gene_code = clean["Gene"]
        if not gene_code:
            continue
        subcategories = split_multi_value(clean["Sub category"]) or [""]
        for category in split_multi_value(clean["Category"]):
            for subcategory in subcategories:
                category_key = (category, subcategory, gene_code)
                if category_key not in seen_category_subcategory_gene:
                    category_gene_row_order[category].append(gene_code)
                    seen_category_subcategory_gene.add(category_key)
        for subcategory in subcategories:
            subcategory_key = (subcategory, gene_code)
            if subcategory_key not in seen_subcategory_gene:
                subcategory_gene_row_order[subcategory].append(gene_code)
                seen_subcategory_gene.add(subcategory_key)
        if gene_code not in seen_gene:
            seen_gene.add(gene_code)
            unique_gene_order.append(gene_code)
        gene_meta_sets[gene_code].add(
            (
                clean["Category"],
                clean["Sub category"],
                clean["Gene name"],
                clean["KO"],
                clean["EC"],
            )
        )
        if clean["UniProt"]:
            query_to_rows_by_uniprot[clean["UniProt"]].append(clean)
        query_to_rows_by_gene[gene_code].append(clean)

    gene_meta: dict[str, dict[str, str]] = {}
    for gene_code, meta_set in gene_meta_sets.items():
        categories = sorted({item[0] for item in meta_set if item[0]})
        sub_categories = sorted({item[1] for item in meta_set if item[1]})
        gene_names = sorted({item[2] for item in meta_set if item[2]})
        kos = sorted({item[3] for item in meta_set if item[3]})
        ecs = sorted({item[4] for item in meta_set if item[4]})
        gene_meta[gene_code] = {
            "Category": "; ".join(categories),
            "Sub category": "; ".join(sub_categories),
            "Gene name": "; ".join(gene_names),
            "KO": "; ".join(kos),
            "EC": "; ".join(ecs),
        }

    return (
        unique_gene_order,
        gene_meta,
        query_to_rows_by_uniprot,
        query_to_rows_by_gene,
        category_gene_row_order,
        subcategory_gene_row_order,
    )


def read_metadata(path: Path):
    header, rows = load_xlsx_rows(path)
    required = ["Assembly Accession", "Genome name", "Phylum", "Class", "Order", "Family", "Genus"]
    missing = [field for field in required if field not in header]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")

    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        clean = {key: normalize_text(value) for key, value in row.items()}
        keys = {
            clean.get("Assembly Accession", ""),
            clean.get("Genome name", ""),
            clean.get("Organism Name", ""),
            clean.get("Species", ""),
        }
        for key in keys:
            if key:
                lookup.setdefault(key, clean)
                lookup.setdefault(canonicalize_genome_id(key), clean)
    return header, lookup


def resolve_query(
    query_id: str,
    by_uniprot: dict[str, list[dict[str, str]]],
    by_gene: dict[str, list[dict[str, str]]],
) -> dict[str, str | list[str]]:
    tokens = [token.strip() for token in query_id.split("|") if token.strip()]

    matched_rows: list[dict[str, str]] = []
    match_mode = ""
    for token in tokens:
        matched_rows.extend(by_uniprot.get(token, []))
    if matched_rows:
        match_mode = "uniprot_token"
    else:
        for token in tokens:
            matched_rows.extend(by_gene.get(token, []))
        if matched_rows:
            match_mode = "gene_token"

    if not matched_rows:
        return {
            "query_id": query_id,
            "match_mode": "unmapped",
            "genes": [],
            "categories": [],
            "sub_categories": [],
        }

    genes = sorted({normalize_text(row["Gene"]) for row in matched_rows if normalize_text(row["Gene"])})
    categories = sorted(
        {normalize_text(row["Category"]) for row in matched_rows if normalize_text(row["Category"])}
    )
    sub_categories = sorted(
        {
            normalize_text(row["Sub category"])
            for row in matched_rows
            if normalize_text(row["Sub category"])
        }
    )
    return {
        "query_id": query_id,
        "match_mode": match_mode,
        "genes": genes,
        "categories": categories,
        "sub_categories": sub_categories,
    }


def process_scored_ko_output(
    raw_path: Path,
    scored_path: Path,
    by_uniprot: dict[str, list[dict[str, str]]],
    by_gene: dict[str, list[dict[str, str]]],
    min_pident: float = 0.0,
    min_qcovs: float = 0.0,
):
    query_cache: dict[str, dict[str, str | list[str]]] = {}
    hit_to_genes: dict[str, set[str]] = defaultdict(set)
    hit_to_queries: dict[str, set[str]] = defaultdict(set)
    hit_to_genome: dict[str, str] = {}
    genome_hit_counter: Counter[str] = Counter()
    raw_total_rows = 0
    scored_total_rows = 0
    besthit_by_hit: dict[str, dict[str, str | int | float]] = {}
    filtered_low_identity_rows = 0
    filtered_low_coverage_rows = 0

    # The raw two-column output carries no pident/qcovs, so the identity/coverage
    # floor can only be evaluated from the scored output. Pre-compute which
    # (query_id, hit_id) pairs survive the floor and use that set to gate both
    # the raw-derived counters below and the best-hit pass further down.
    passing_pairs: set[tuple[str, str]] = set()
    with scored_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            query_id = parts[0].strip()
            hit_id = parts[1].strip()
            pident = float(parts[2])
            qcovs = float(parts[3])
            if pident < min_pident:
                filtered_low_identity_rows += 1
                continue
            if qcovs < min_qcovs:
                filtered_low_coverage_rows += 1
                continue
            passing_pairs.add((query_id, hit_id))

    with raw_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            raw_total_rows += 1
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            query_id = parts[0].strip()
            hit_id = parts[1].strip()
            if (query_id, hit_id) not in passing_pairs:
                continue
            genome_id = canonicalize_genome_id(hit_id.rsplit("__", 1)[-1])
            hit_to_genome[hit_id] = genome_id
            hit_to_queries[hit_id].add(query_id)
            genome_hit_counter[genome_id] += 1

            if query_id not in query_cache:
                query_cache[query_id] = resolve_query(query_id, by_uniprot, by_gene)

            for gene in query_cache[query_id]["genes"]:
                hit_to_genes[hit_id].add(gene)

    def better_besthit(candidate, incumbent) -> bool:
        if incumbent is None:
            return True
        candidate_key = (
            float(candidate["bitscore"]),
            float("inf") if float(candidate["evalue"]) == 0 else -math.log10(float(candidate["evalue"])),
            float(candidate["pident"]),
            float(candidate["qcovs"]),
            int(candidate["length"]),
            str(candidate["query_id"]),
        )
        incumbent_key = (
            float(incumbent["bitscore"]),
            float("inf") if float(incumbent["evalue"]) == 0 else -math.log10(float(incumbent["evalue"])),
            float(incumbent["pident"]),
            float(incumbent["qcovs"]),
            int(incumbent["length"]),
            str(incumbent["query_id"]),
        )
        return candidate_key > incumbent_key

    with scored_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            query_id = parts[0].strip()
            hit_id = parts[1].strip()
            if (query_id, hit_id) not in passing_pairs:
                continue
            scored_total_rows += 1
            if query_id not in query_cache:
                query_cache[query_id] = resolve_query(query_id, by_uniprot, by_gene)
            mapping = query_cache[query_id]
            candidate = {
                "query_id": query_id,
                "hit_protein_id": hit_id,
                "genome_id": hit_to_genome.get(hit_id, canonicalize_genome_id(hit_id.rsplit("__", 1)[-1])),
                "pident": float(parts[2]),
                "qcovs": float(parts[3]),
                "length": int(float(parts[4])),
                "evalue": float(parts[5]),
                "bitscore": float(parts[6]),
                "mapped_gene_count": len(mapping["genes"]),
                "mapped_genes": list(mapping["genes"]),
                "mapped_categories": list(mapping["categories"]),
                "mapped_sub_categories": list(mapping["sub_categories"]),
                "match_mode": str(mapping["match_mode"]),
            }
            incumbent = besthit_by_hit.get(hit_id)
            if better_besthit(candidate, incumbent):
                besthit_by_hit[hit_id] = candidate

    return (
        query_cache,
        hit_to_genes,
        hit_to_queries,
        hit_to_genome,
        genome_hit_counter,
        raw_total_rows,
        scored_total_rows,
        besthit_by_hit,
        filtered_low_identity_rows,
        filtered_low_coverage_rows,
    )


def ordered_unique(items: Iterable[str]) -> list[str]:
    seen = set()
    ordered = []
    for item in items:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def split_multi_value(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def slugify_label(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return slug or "unassigned"


def rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    return "#" + "".join(f"{round(max(0, min(1, channel)) * 255):02x}" for channel in rgb)


def hex_to_rgb01(color: str) -> tuple[float, float, float]:
    color = color.lstrip("#")
    return tuple(int(color[index : index + 2], 16) / 255.0 for index in (0, 2, 4))


def srgb_to_linear(channel: float) -> float:
    if channel <= 0.04045:
        return channel / 12.92
    return ((channel + 0.055) / 1.055) ** 2.4


def rgb01_to_lab(rgb: tuple[float, float, float]) -> tuple[float, float, float]:
    r, g, b = (srgb_to_linear(channel) for channel in rgb)
    x = (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / 0.95047
    y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b
    z = (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / 1.08883

    def f(value: float) -> float:
        if value > 0.008856:
            return value ** (1 / 3)
        return (7.787 * value) + (16 / 116)

    fx, fy, fz = f(x), f(y), f(z)
    return ((116 * fy) - 16, 500 * (fx - fy), 200 * (fy - fz))


def color_distance_sq(color_a: tuple[float, float, float], color_b: tuple[float, float, float]) -> float:
    return sum((left - right) ** 2 for left, right in zip(color_a, color_b))


QUALITATIVE_TAXONOMY_PALETTE = [
    "#1f77b4",
    "#d62728",
    "#ff7f0e",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#17becf",
    "#bcbd22",
    "#7f7f7f",
    "#2ca02c",
    "#393b79",
    "#ad494a",
    "#c49c94",
    "#8c6d31",
    "#843c39",
    "#7b4173",
    "#3182bd",
    "#e6550d",
    "#f1ce63",
    "#756bb1",
    "#636363",
    "#dd1c77",
    "#6baed6",
    "#fd8d3c",
]


def build_high_contrast_palette(size: int, unassigned_color: str) -> list[str]:
    seed_colors = QUALITATIVE_TAXONOMY_PALETTE[: min(size, len(QUALITATIVE_TAXONOMY_PALETTE))]
    if size <= len(seed_colors):
        return seed_colors

    candidate_hexes = []
    for hue_steps in (24, 36):
        for hue_index in range(hue_steps):
            hue = hue_index / hue_steps
            for saturation in (0.45, 0.60, 0.75, 0.90):
                for lightness in (0.34, 0.45, 0.56, 0.67):
                    rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
                    candidate_hexes.append(rgb_to_hex(rgb))

    deduped_candidates = []
    seen = set(seed_colors)
    seen.add(unassigned_color)
    for color in candidate_hexes:
        if color not in seen:
            deduped_candidates.append(color)
            seen.add(color)

    candidate_pool = [(color, rgb01_to_lab(hex_to_rgb01(color))) for color in deduped_candidates]
    selected = list(seed_colors)
    selected_labs = [rgb01_to_lab(hex_to_rgb01(color)) for color in selected]

    while len(selected) < size and candidate_pool:
        best_index = 0
        best_distance = -1.0
        for index, (_, candidate_lab) in enumerate(candidate_pool):
            min_distance = min(color_distance_sq(candidate_lab, selected_lab) for selected_lab in selected_labs)
            if min_distance > best_distance:
                best_distance = min_distance
                best_index = index
        chosen_color, chosen_lab = candidate_pool.pop(best_index)
        selected.append(chosen_color)
        selected_labs.append(chosen_lab)
    return selected[:size]


def make_color_map(
    values: list[str],
    unassigned_color: str = "#bdbdbd",
    prefer_qualitative: bool = False,
    prefer_high_contrast: bool = False,
) -> dict[str, str]:
    uniq = sorted(set(values))
    color_map: dict[str, str] = {}
    assigned_values = [value for value in uniq if value != "Unassigned"]
    if prefer_qualitative and len(assigned_values) <= len(QUALITATIVE_TAXONOMY_PALETTE):
        for value, color in zip(assigned_values, QUALITATIVE_TAXONOMY_PALETTE):
            color_map[value] = color
    elif prefer_high_contrast:
        palette = build_high_contrast_palette(len(assigned_values), unassigned_color)
        for value, color in zip(assigned_values, palette):
            color_map[value] = color
    else:
        for index, value in enumerate(assigned_values):
            hue = ((index * 137.508) % 360) / 360.0
            saturation = 0.55 + 0.15 * (index % 3)
            lightness = 0.48 + 0.06 * (index % 2)
            rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
            color_map[value] = rgb_to_hex(rgb)
    if "Unassigned" in uniq:
        color_map["Unassigned"] = unassigned_color
    return color_map


# Fixed phylum palette, shared with the MAG abundance figures so phylum colors match across figures.
PHYLUM_COLORS = {
    "Acidobacteriota": "#1f77b4",
    "Actinomycetota": "#ff7f0e",
    "Armatimonadota": "#2ca02c",
    "Ascomycota": "#ffb6c1",
    "Bacteroidota": "#9467bd",
    "Bdellovibrionota": "#EEC900",
    "Cyanobacteriota": "#8c564b",
    "Deinococcota": "#b22222",
    "Fibrobacterota": "#ee82ee",
    "Gemmatimonadota": "#17becf",
    "Myxococcota": "#228b22",
    "Oomycota": "#a0522d",
    "Planctomycetota": "#c49c94",
    "Pseudomonadota": "#8B0A50",
    "Spirochaetota": "#98df8a",
    "Verrucomicrobiota": "#dbdb8d",
    "Arthropoda": "#ff4500",
    "Deferribacterota": "#ff69b4",
    "Desulfobacterota": "#ffd700",
    "Bacillota": "#d62728",
}

# Category-panel ITOL figures are presented side by side, so they use one shared
# legend range rather than per-panel autoscaling.
CATEGORY_HEATMAP_SHARED_MAX = 24


def make_phylum_color_map(values: list[str], unassigned_color: str = "#bdbdbd") -> dict[str, str]:
    uniq = sorted(set(values))
    color_map: dict[str, str] = {}
    unmapped = [value for value in uniq if value != "Unassigned" and value not in PHYLUM_COLORS]
    for value in uniq:
        if value in PHYLUM_COLORS:
            color_map[value] = PHYLUM_COLORS[value]
    if unmapped:
        reserved = set(PHYLUM_COLORS.values())
        reserved.add(unassigned_color)
        pool_size = len(unmapped)
        palette: list[str] = []
        while len(palette) < len(unmapped):
            pool_size += len(reserved) + 8
            candidates = build_high_contrast_palette(pool_size, unassigned_color)
            palette = [color for color in candidates if color not in reserved]
        for value, color in zip(unmapped, palette):
            color_map[value] = color
    if "Unassigned" in uniq:
        color_map["Unassigned"] = unassigned_color
    return color_map


def build_level_order(gene_order: list[str], gene_meta: dict[str, dict[str, str]], field: str) -> list[str]:
    ordered = []
    for gene in gene_order:
        values = split_multi_value(normalize_text(gene_meta[gene].get(field, "")))
        ordered.extend(values)
    return ordered_unique(ordered)


def count_hits(
    genomes: list[str],
    gene_order: list[str],
    gene_meta: dict[str, dict[str, str]],
    hit_to_genes: dict[str, set[str]],
    hit_to_genome: dict[str, str],
    hit_to_queries: dict[str, set[str]],
    query_cache: dict[str, dict[str, str | list[str]]],
    besthit_by_hit: dict[str, dict[str, str | int | float]],
):
    inclusive_counts: dict[str, Counter[str]] = {genome: Counter() for genome in genomes}
    strict_counts: dict[str, Counter[str]] = {genome: Counter() for genome in genomes}
    besthit_counts: dict[str, Counter[str]] = {genome: Counter() for genome in genomes}

    ambiguous_hits = []
    strict_hit_total = 0
    inclusive_hit_total = 0
    besthit_hit_total = 0
    besthit_unmapped_total = 0

    for hit_id, genes in hit_to_genes.items():
        genome = hit_to_genome[hit_id]
        if genome not in inclusive_counts:
            inclusive_counts[genome] = Counter()
            strict_counts[genome] = Counter()
            besthit_counts[genome] = Counter()
        for gene in genes:
            inclusive_counts[genome][gene] += 1
            inclusive_hit_total += 1

        if len(genes) == 1:
            gene = next(iter(genes))
            strict_counts[genome][gene] += 1
            strict_hit_total += 1
        elif len(genes) > 1:
            categories = sorted({gene_meta[gene]["Category"] for gene in genes if gene in gene_meta})
            sub_categories = sorted(
                {gene_meta[gene]["Sub category"] for gene in genes if gene in gene_meta}
            )
            ambiguous_hits.append(
                {
                    "hit_protein_id": hit_id,
                    "genome_id": genome,
                    "mapped_gene_count": str(len(genes)),
                    "mapped_genes": "; ".join(sorted(genes)),
                    "mapped_categories": "; ".join(cat for cat in categories if cat),
                    "mapped_sub_categories": "; ".join(sub for sub in sub_categories if sub),
                }
            )

    besthit_assignment_rows = []
    besthit_resolved_ambiguous_rows = []
    for hit_id, candidate_genes_set in hit_to_genes.items():
        best = besthit_by_hit.get(hit_id)
        if best is None and len(candidate_genes_set) == 1:
            chosen_query_id = sorted(hit_to_queries.get(hit_id, []))[0] if hit_to_queries.get(hit_id) else ""
            mapping = query_cache.get(
                chosen_query_id,
                {"genes": sorted(candidate_genes_set), "categories": [], "sub_categories": [], "match_mode": ""},
            )
            best = {
                "query_id": chosen_query_id,
                "hit_protein_id": hit_id,
                "genome_id": hit_to_genome[hit_id],
                "pident": "",
                "qcovs": "",
                "length": "",
                "evalue": "",
                "bitscore": "",
                "mapped_gene_count": len(mapping["genes"]),
                "mapped_genes": list(mapping["genes"]),
                "mapped_categories": list(mapping["categories"]),
                "mapped_sub_categories": list(mapping["sub_categories"]),
                "match_mode": str(mapping["match_mode"]),
            }
        elif best is None:
            best = {
                "query_id": "",
                "hit_protein_id": hit_id,
                "genome_id": hit_to_genome[hit_id],
                "pident": "",
                "qcovs": "",
                "length": "",
                "evalue": "",
                "bitscore": "",
                "mapped_gene_count": 0,
                "mapped_genes": [],
                "mapped_categories": [],
                "mapped_sub_categories": [],
                "match_mode": "unresolved_multigene",
            }

        genome = str(best["genome_id"])
        chosen_genes = list(best["mapped_genes"])
        chosen_gene = chosen_genes[0] if len(chosen_genes) == 1 else ""
        if chosen_gene:
            besthit_counts[genome][chosen_gene] += 1
            besthit_hit_total += 1
        else:
            besthit_unmapped_total += 1
        candidate_genes = sorted(candidate_genes_set)
        candidate_categories = sorted(
            {gene_meta[gene]["Category"] for gene in candidate_genes if gene in gene_meta and gene_meta[gene]["Category"]}
        )
        candidate_subcategories = sorted(
            {
                gene_meta[gene]["Sub category"]
                for gene in candidate_genes
                if gene in gene_meta and gene_meta[gene]["Sub category"]
            }
        )
        row = {
            "hit_protein_id": hit_id,
            "genome_id": genome,
            "chosen_query_id": str(best["query_id"]),
            "chosen_gene": chosen_gene,
            "chosen_category": "; ".join(best["mapped_categories"]),
            "chosen_sub_category": "; ".join(best["mapped_sub_categories"]),
            "chosen_match_mode": str(best["match_mode"]),
            "bitscore": (
                f"{float(best['bitscore']):.6f}" if str(best["bitscore"]) not in {"", "None"} else ""
            ),
            "evalue": (
                f"{float(best['evalue']):.6g}" if str(best["evalue"]) not in {"", "None"} else ""
            ),
            "pident": (
                f"{float(best['pident']):.6f}" if str(best["pident"]) not in {"", "None"} else ""
            ),
            "qcovs": (
                f"{float(best['qcovs']):.6f}" if str(best["qcovs"]) not in {"", "None"} else ""
            ),
            "alignment_length": (
                str(int(float(best["length"]))) if str(best["length"]) not in {"", "None"} else ""
            ),
            "candidate_gene_count": str(len(candidate_genes)),
            "candidate_genes": "; ".join(candidate_genes),
            "candidate_categories": "; ".join(candidate_categories),
            "candidate_sub_categories": "; ".join(candidate_subcategories),
            "chosen_query_mapped": "yes" if chosen_gene else "no",
        }
        besthit_assignment_rows.append(row)
        if len(candidate_genes) > 1:
            besthit_resolved_ambiguous_rows.append(row)

    category_order = build_level_order(gene_order, gene_meta, "Category")
    subcategory_order = build_level_order(gene_order, gene_meta, "Sub category")

    def aggregate(
        source: dict[str, Counter[str]],
        field: str,
        ordered_fields: list[str],
    ) -> dict[str, Counter[str]]:
        aggregated: dict[str, Counter[str]] = {genome: Counter() for genome in source}
        for genome, counter in source.items():
            for gene, value in counter.items():
                labels = split_multi_value(normalize_text(gene_meta[gene].get(field, "")))
                for label in labels:
                    aggregated[genome][label] += value
        for genome in list(aggregated):
            for label in ordered_fields:
                aggregated[genome][label] += 0
        return aggregated

    return {
        "gene_order": gene_order,
        "category_order": category_order,
        "subcategory_order": subcategory_order,
        "besthit_gene_counts": besthit_counts,
        "inclusive_gene_counts": inclusive_counts,
        "strict_gene_counts": strict_counts,
        "besthit_category_counts": aggregate(besthit_counts, "Category", category_order),
        "inclusive_category_counts": aggregate(inclusive_counts, "Category", category_order),
        "strict_category_counts": aggregate(strict_counts, "Category", category_order),
        "besthit_subcategory_counts": aggregate(besthit_counts, "Sub category", subcategory_order),
        "inclusive_subcategory_counts": aggregate(
            inclusive_counts, "Sub category", subcategory_order
        ),
        "strict_subcategory_counts": aggregate(strict_counts, "Sub category", subcategory_order),
        "ambiguous_hits": ambiguous_hits,
        "besthit_assignment_rows": besthit_assignment_rows,
        "besthit_resolved_ambiguous_rows": besthit_resolved_ambiguous_rows,
        "besthit_hit_total": besthit_hit_total,
        "besthit_unmapped_total": besthit_unmapped_total,
        "inclusive_hit_total": inclusive_hit_total,
        "strict_hit_total": strict_hit_total,
    }


def write_matrix(
    path: Path,
    genomes: list[str],
    labels: list[str],
    counts: dict[str, Counter[str]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", *labels])
        for genome in genomes:
            writer.writerow([genome, *[counts.get(genome, Counter()).get(label, 0) for label in labels]])


def build_summary_rows(
    labels: list[str],
    genomes: list[str],
    besthit_counts: dict[str, Counter[str]],
    strict_counts: dict[str, Counter[str]],
    inclusive_counts: dict[str, Counter[str]],
    label_meta: dict[str, dict[str, str]],
    level_name: str,
) -> list[dict[str, str]]:
    rows = []
    for label in labels:
        besthit_values = [besthit_counts.get(genome, Counter()).get(label, 0) for genome in genomes]
        strict_values = [strict_counts.get(genome, Counter()).get(label, 0) for genome in genomes]
        inclusive_values = [inclusive_counts.get(genome, Counter()).get(label, 0) for genome in genomes]
        besthit_present = [value for value in besthit_values if value > 0]
        strict_present = [value for value in strict_values if value > 0]
        inclusive_present = [value for value in inclusive_values if value > 0]
        meta = label_meta.get(label, {})
        rows.append(
            {
                "Level": level_name,
                "Label": label,
                "Category": meta.get("Category", ""),
                "Sub_category": meta.get("Sub category", ""),
                "Gene_name": meta.get("Gene name", ""),
                "KO": meta.get("KO", ""),
                "EC": meta.get("EC", ""),
                "Genome_count": str(len(genomes)),
                "Besthit_total": str(sum(besthit_values)),
                "Besthit_genomes_with_hits": str(sum(1 for value in besthit_values if value > 0)),
                "Besthit_average_per_genome": f"{mean(besthit_values):.6f}",
                "Besthit_average_if_present": f"{mean(besthit_present):.6f}" if besthit_present else "0",
                "Besthit_max_in_single_genome": str(max(besthit_values) if besthit_values else 0),
                "Strict_total": str(sum(strict_values)),
                "Strict_genomes_with_hits": str(sum(1 for value in strict_values if value > 0)),
                "Strict_average_per_genome": f"{mean(strict_values):.6f}",
                "Strict_average_if_present": f"{mean(strict_present):.6f}" if strict_present else "0",
                "Strict_max_in_single_genome": str(max(strict_values) if strict_values else 0),
                "Inclusive_total": str(sum(inclusive_values)),
                "Inclusive_genomes_with_hits": str(sum(1 for value in inclusive_values if value > 0)),
                "Inclusive_average_per_genome": f"{mean(inclusive_values):.6f}",
                "Inclusive_average_if_present": (
                    f"{mean(inclusive_present):.6f}" if inclusive_present else "0"
                ),
                "Inclusive_max_in_single_genome": str(max(inclusive_values) if inclusive_values else 0),
            }
        )
    return rows


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_itol_heatmap(
    path: Path,
    label: str,
    genomes: list[str],
    fields: list[str],
    counts: dict[str, Counter[str]],
    show_labels: bool,
    fixed_max_value: int | None = None,
) -> None:
    max_value = max(
        (counts.get(genome, Counter()).get(field, 0) for genome in genomes for field in fields),
        default=0,
    )
    if fixed_max_value is not None:
        max_value = fixed_max_value
    header = [
        "DATASET_HEATMAP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{label}",
        "COLOR\t#2166ac",
        "FIELD_LABELS\t" + "\t".join(fields),
        "COLOR_MIN\t#f7fbff",
        "COLOR_MAX\t#08306b",
        "COLOR_NAN\t#d9d9d9",
        "AUTO_LEGEND\t1",
        "SHOW_INTERNAL\t0",
        f"SHOW_LABELS\t{1 if show_labels else 0}",
        "LABEL_ROTATION\t90",
        "STRIP_WIDTH\t18" if show_labels else "STRIP_WIDTH\t10",
        "MARGIN\t20",
        "BORDER_WIDTH\t0",
        "USER_MIN_VALUE\t0",
        f"USER_MAX_VALUE\t{max_value}",
        "DATA",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(header))
        handle.write("\n")
        for genome in genomes:
            values = [str(counts.get(genome, Counter()).get(field, 0)) for field in fields]
            handle.write(genome)
            handle.write("\t")
            handle.write("\t".join(values))
            handle.write("\n")


def write_itol_colorstrip(
    path: Path,
    dataset_label: str,
    genomes: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    rank_field: str,
    strip_width: int,
    label_field: str | None = None,
    color_map: dict[str, str] | None = None,
) -> dict[str, str]:
    values = []
    for genome in genomes:
        value = normalize_text(metadata_lookup[genome].get(rank_field, "")) or "Unassigned"
        values.append(value)
    if color_map is None:
        color_map = make_color_map(
            values,
            prefer_qualitative=(rank_field == "Phylum"),
            prefer_high_contrast=(rank_field == "Genus"),
        )

    legend_values = sorted(set(values))
    header = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{dataset_label}",
        "COLOR\t#444444",
        f"STRIP_WIDTH\t{strip_width}",
        "MARGIN\t5",
        "SHOW_INTERNAL\t0",
        f"LEGEND_TITLE\t{rank_field}",
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_values)),
        "LEGEND_COLORS\t" + "\t".join(color_map[value] for value in legend_values),
        "LEGEND_LABELS\t" + "\t".join(legend_values),
        "DATA",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(header))
        handle.write("\n")
        for genome in genomes:
            value = normalize_text(metadata_lookup[genome].get(rank_field, "")) or "Unassigned"
            label_value = value
            if label_field:
                label_value = normalize_text(metadata_lookup[genome].get(label_field, "")) or value
            handle.write(f"{genome}\t{color_map[value]}\t{label_value}\n")
    return color_map


def write_itol_piechart(
    path: Path,
    dataset_label: str,
    genomes: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    value_field: str,
    field_labels: tuple[str, str] = ("Completeness", "Incomplete"),
    field_colors: tuple[str, str] = ("#ff0000", "#00ff00"),
    radius: float = 1,
) -> None:
    # A 2-slice external
    # pie chart per leaf (value, 100 - value), COMMA-separated, fixed radius 1.
    header = [
        "DATASET_PIECHART",
        "SEPARATOR COMMA",
        f"DATASET_LABEL,{dataset_label}",
        f"COLOR,{field_colors[0]}",
        f"FIELD_COLORS,{field_colors[0]},{field_colors[1]}",
        f"FIELD_LABELS,{field_labels[0]},{field_labels[1]}",
        "DATA",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(header))
        handle.write("\n")
        for genome in genomes:
            raw_value = normalize_text(metadata_lookup[genome].get(value_field, ""))
            try:
                completeness = float(raw_value)
            except (TypeError, ValueError):
                completeness = 0.0
            incomplete = max(0.0, 100.0 - completeness)
            handle.write(f"{genome},-1,{radius:g},{completeness:g},{incomplete:g}\n")


def write_itol_tree_colors(
    path: Path,
    tree_path: Path,
    metadata_lookup: dict[str, dict[str, str]],
    rank_field: str,
    color_map: dict[str, str],
) -> None:
    from Bio import Phylo

    tree = Phylo.read(tree_path, "newick")
    leaf_values = {
        terminal.name: normalize_text(metadata_lookup[terminal.name].get(rank_field, "")) or "Unassigned"
        for terminal in tree.get_terminals()
    }

    annotations: list[tuple[str, str, str]] = []

    def descend(clade, parent_value: str | None = None) -> None:
        terminals = clade.get_terminals()
        terminal_names = [terminal.name for terminal in terminals if terminal.name]
        if not terminal_names:
            return

        values = {leaf_values[name] for name in terminal_names}
        if len(values) == 1:
            value = next(iter(values))
            if value != parent_value:
                if len(terminal_names) == 1:
                    node_id = terminal_names[0]
                else:
                    node_id = f"{terminal_names[0]}|{terminal_names[-1]}"
                annotations.append((node_id, color_map[value], value))
            return

        for child in clade.clades:
            descend(child)

    descend(tree.root)

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("TREE_COLORS\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATA\n")
        for node_id, color, value in annotations:
            handle.write(f"{node_id}\trange\t{color}\t{value}\n")


def write_taxonomy_legend(path: Path, rank_field: str, color_map: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([rank_field, "ColorHex"])
        for label in sorted(color_map):
            writer.writerow([label, color_map[label]])


def write_category_gene_heatmaps(
    output_dir: Path,
    genomes: list[str],
    gene_order: list[str],
    gene_meta: dict[str, dict[str, str]],
    category_gene_row_order: dict[str, list[str]],
    gene_counts: dict[str, Counter[str]],
    category_order: list[str],
) -> list[dict[str, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for category in category_order:
        category_genes = list(category_gene_row_order.get(category, []))
        if not category_genes:
            category_genes = [
                gene
                for gene in gene_order
                if category in split_multi_value(normalize_text(gene_meta[gene].get("Category", "")))
            ]
        if not category_genes:
            continue
        genomes_with_hits = [
            genome
            for genome in genomes
            if any(gene_counts.get(genome, Counter()).get(gene, 0) > 0 for gene in category_genes)
        ]
        file_name = f"itol_gene_counts_besthit__{slugify_label(category)}.txt"
        write_itol_heatmap(
            output_dir / file_name,
            f"gene_counts_besthit__{category}",
            genomes_with_hits,
            category_genes,
            gene_counts,
            show_labels=True,
        )
        manifest_rows.append(
            {
                "Category": category,
                "HeatmapFileName": file_name,
                "GeneCount": str(len(category_genes)),
            }
        )
    return manifest_rows


def write_category_pruned_trees(
    output_dir: Path,
    tree_path: Path,
    genomes: list[str],
    category_counts: dict[str, Counter[str]],
    category_order: list[str],
) -> list[dict[str, str]]:
    from Bio import Phylo

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    tree_text = tree_path.read_text(encoding="utf-8")
    tree_stem = tree_path.stem
    for category in category_order:
        genomes_with_hits = [
            genome
            for genome in genomes
            if category_counts.get(genome, Counter()).get(category, 0) > 0
        ]
        if not genomes_with_hits:
            continue
        keep_labels = set(genomes_with_hits)
        tree = Phylo.read(StringIO(tree_text), "newick")
        for terminal in list(tree.get_terminals()):
            if terminal.name not in keep_labels:
                tree.prune(target=terminal)
        file_name = f"{tree_stem}__{slugify_label(category)}.nwk"
        Phylo.write(tree, output_dir / file_name, "newick")
        manifest_rows.append(
            {
                "Category": category,
                "TreeFileName": file_name,
                "GenomeCount": str(len(genomes_with_hits)),
            }
        )
    return manifest_rows


def write_category_annotation_files(
    output_dir: Path,
    tree_dir: Path,
    tree_stem: str,
    genomes: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    category_counts: dict[str, Counter[str]],
    category_order: list[str],
    phylum_colors: dict[str, str],
    genus_colors: dict[str, str],
    source_colors: dict[str, str] | None,
    completeness_available: bool = False,
) -> list[dict[str, str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for category in category_order:
        genomes_with_hits = [
            genome
            for genome in genomes
            if category_counts.get(genome, Counter()).get(category, 0) > 0
        ]
        if not genomes_with_hits:
            continue
        slug = slugify_label(category)
        tree_file = tree_dir / f"{tree_stem}__{slug}.nwk"

        phylum_file = f"itol_phylum_annotation__{slug}.txt"
        write_itol_tree_colors(
            output_dir / phylum_file,
            tree_file,
            metadata_lookup,
            "Phylum",
            phylum_colors,
        )
        write_itol_tree_colors(
            output_dir / f"itol_phylum_colorstrip__{slug}.txt",
            tree_file,
            metadata_lookup,
            "Phylum",
            phylum_colors,
        )

        genus_file = f"itol_genus_annotation__{slug}.txt"
        write_itol_tree_colors(
            output_dir / genus_file,
            tree_file,
            metadata_lookup,
            "Genus",
            genus_colors,
        )
        write_itol_tree_colors(
            output_dir / f"itol_genus_colorstrip__{slug}.txt",
            tree_file,
            metadata_lookup,
            "Genus",
            genus_colors,
        )

        source_file = ""
        if source_colors is not None:
            source_file = f"itol_source_colorstrip__{slug}.txt"
            write_itol_colorstrip(
                output_dir / source_file,
                f"Genome source __ {category}",
                genomes_with_hits,
                metadata_lookup,
                "Source",
                strip_width=28,
                color_map=source_colors,
            )

        label_file = f"itol_genome_name_labels__{slug}.txt"
        write_itol_labels(
            output_dir / label_file,
            genomes_with_hits,
            metadata_lookup,
            "Genome name",
        )

        completeness_file = ""
        if completeness_available:
            completeness_file = f"itol_completeness_piechart__{slug}.txt"
            write_itol_piechart(
                output_dir / completeness_file,
                f"Genome completeness (%) __ {category}",
                genomes_with_hits,
                metadata_lookup,
                "Completeness (%)",
            )

        manifest_rows.append(
            {
                "Category": category,
                "PhylumAnnotationFileName": phylum_file,
                "GenusAnnotationFileName": genus_file,
                "SourceStripFileName": source_file,
                "GenomeNameLabelsFileName": label_file,
                "CompletenessPiechartFileName": completeness_file,
            }
        )
    return manifest_rows


def format_percent_slug(value: float) -> str:
    percent = value * 100.0
    if abs(percent - round(percent)) < 1e-9:
        return f"{int(round(percent))}pct"
    return f"{f'{percent:.3f}'.rstrip('0').rstrip('.').replace('.', 'p')}pct"


def build_filtered_category_distinct_gene_specs(
    genomes: list[str],
    gene_order: list[str],
    gene_meta: dict[str, dict[str, str]],
    category_gene_row_order: dict[str, list[str]],
    gene_counts: dict[str, Counter[str]],
    category_order: list[str],
    min_fraction: float,
) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []
    for category in category_order:
        category_genes = ordered_unique(category_gene_row_order.get(category, []))
        if not category_genes:
            category_genes = ordered_unique(
                [
                    gene
                    for gene in gene_order
                    if category in split_multi_value(normalize_text(gene_meta[gene].get("Category", "")))
                ]
            )
        if not category_genes:
            continue

        distinct_gene_count = len(category_genes)
        min_distinct_genes_required = max(1, math.ceil(distinct_gene_count * min_fraction))
        genomes_with_any_hits = []
        genomes_passing_filter = []
        for genome in genomes:
            represented_distinct_genes = sum(
                1
                for gene in category_genes
                if gene_counts.get(genome, Counter()).get(gene, 0) > 0
            )
            if represented_distinct_genes > 0:
                genomes_with_any_hits.append(genome)
            if represented_distinct_genes >= min_distinct_genes_required:
                genomes_passing_filter.append(genome)

        specs.append(
            {
                "category": category,
                "category_slug": slugify_label(category),
                "category_genes": category_genes,
                "distinct_gene_count": distinct_gene_count,
                "min_distinct_genes_required": min_distinct_genes_required,
                "genomes_with_any_hits": genomes_with_any_hits,
                "genomes_passing_filter": genomes_passing_filter,
                "cutoff_tag": f"min{min_distinct_genes_required}of{distinct_gene_count}",
            }
        )
    return specs


def write_filtered_category_distinct_gene_bundles(
    output_dir: Path,
    tree_path: Path,
    genomes: list[str],
    gene_order: list[str],
    gene_meta: dict[str, dict[str, str]],
    category_gene_row_order: dict[str, list[str]],
    gene_counts: dict[str, Counter[str]],
    category_order: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    phylum_colors: dict[str, str],
    genus_colors: dict[str, str],
    source_colors: dict[str, str] | None,
    min_fraction: float,
    completeness_available: bool = False,
) -> dict[str, object]:
    from Bio import Phylo

    if min_fraction <= 0:
        raise ValueError("min_fraction must be > 0 for filtered category bundle export")

    fraction_slug = format_percent_slug(min_fraction)
    bundle_suffix = f"min{fraction_slug}_distinct_genes"
    heatmap_dir = output_dir / f"itol_gene_counts_besthit_by_category_{bundle_suffix}"
    tree_dir = output_dir / f"itol_gene_counts_besthit_by_category_{bundle_suffix}_trees"
    annotation_dir = output_dir / f"itol_gene_counts_besthit_by_category_{bundle_suffix}_annotations"
    manifest_name = f"itol_gene_counts_besthit_by_category_{bundle_suffix}.tsv"
    readme_name = f"README_itol_gene_counts_besthit_by_category_{bundle_suffix}.txt"

    heatmap_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)

    tree_text = tree_path.read_text(encoding="utf-8")
    tree_stem = tree_path.stem
    manifest_rows: list[dict[str, str]] = []
    empty_specs: list[dict[str, object]] = []

    for spec in build_filtered_category_distinct_gene_specs(
        genomes,
        gene_order,
        gene_meta,
        category_gene_row_order,
        gene_counts,
        category_order,
        min_fraction,
    ):
        genomes_passing_filter = list(spec["genomes_passing_filter"])
        if not genomes_passing_filter:
            empty_specs.append(spec)
            continue

        category = str(spec["category"])
        category_slug = str(spec["category_slug"])
        category_genes = list(spec["category_genes"])
        distinct_gene_count = int(spec["distinct_gene_count"])
        min_distinct_genes_required = int(spec["min_distinct_genes_required"])
        genomes_with_any_hits = list(spec["genomes_with_any_hits"])
        cutoff_tag = str(spec["cutoff_tag"])
        file_suffix = f"{category_slug}__{cutoff_tag}"

        heatmap_name = f"itol_gene_counts_besthit__{file_suffix}.txt"
        tree_name = f"{tree_stem}__{file_suffix}.nwk"
        phylum_annotation_name = f"itol_phylum_annotation__{file_suffix}.txt"
        genus_annotation_name = f"itol_genus_annotation__{file_suffix}.txt"
        phylum_compat_name = f"itol_phylum_colorstrip__{file_suffix}.txt"
        genus_compat_name = f"itol_genus_colorstrip__{file_suffix}.txt"
        source_name = f"itol_source_colorstrip__{file_suffix}.txt"
        labels_name = f"itol_genome_name_labels__{file_suffix}.txt"
        completeness_name = f"itol_completeness_piechart__{file_suffix}.txt"

        write_itol_heatmap(
            heatmap_dir / heatmap_name,
            f"gene_counts_besthit__{category}__{cutoff_tag}",
            genomes_passing_filter,
            category_genes,
            gene_counts,
            show_labels=True,
            fixed_max_value=CATEGORY_HEATMAP_SHARED_MAX,
        )

        tree = Phylo.read(StringIO(tree_text), "newick")
        keep_labels = set(genomes_passing_filter)
        for terminal in list(tree.get_terminals()):
            if terminal.name not in keep_labels:
                tree.prune(target=terminal)
        tree_file = tree_dir / tree_name
        Phylo.write(tree, tree_file, "newick")

        write_itol_tree_colors(
            annotation_dir / phylum_annotation_name,
            tree_file,
            metadata_lookup,
            "Phylum",
            phylum_colors,
        )
        write_itol_tree_colors(
            annotation_dir / phylum_compat_name,
            tree_file,
            metadata_lookup,
            "Phylum",
            phylum_colors,
        )
        write_itol_tree_colors(
            annotation_dir / genus_annotation_name,
            tree_file,
            metadata_lookup,
            "Genus",
            genus_colors,
        )
        write_itol_tree_colors(
            annotation_dir / genus_compat_name,
            tree_file,
            metadata_lookup,
            "Genus",
            genus_colors,
        )

        source_file = ""
        if source_colors is not None:
            source_file = source_name
            write_itol_colorstrip(
                annotation_dir / source_name,
                f"Genome source __ {category} __ {cutoff_tag}",
                genomes_passing_filter,
                metadata_lookup,
                "Source",
                strip_width=28,
                color_map=source_colors,
            )

        write_itol_labels(
            annotation_dir / labels_name,
            genomes_passing_filter,
            metadata_lookup,
            "Genome name",
        )

        completeness_file = ""
        if completeness_available:
            completeness_file = completeness_name
            write_itol_piechart(
                annotation_dir / completeness_name,
                f"Genome completeness (%) __ {category} __ {cutoff_tag}",
                genomes_passing_filter,
                metadata_lookup,
                "Completeness (%)",
            )

        manifest_rows.append(
            {
                "Category": category,
                "CategorySlug": category_slug,
                "DistinctGeneCount": str(distinct_gene_count),
                "MinDistinctGenesRequired": str(min_distinct_genes_required),
                "MinDistinctGeneFraction": f"{min_fraction:.6f}",
                "GenomeCountWithAnyHits": str(len(genomes_with_any_hits)),
                "GenomeCount": str(len(genomes_passing_filter)),
                "CutoffTag": cutoff_tag,
                "HeatmapFileName": heatmap_name,
                "TreeFileName": tree_name,
                "PhylumAnnotationFileName": phylum_annotation_name,
                "GenusAnnotationFileName": genus_annotation_name,
                "SourceStripFileName": source_file,
                "GenomeNameLabelsFileName": labels_name,
                "CompletenessPiechartFileName": completeness_file,
            }
        )

    write_rows(
        output_dir / manifest_name,
        manifest_rows,
        [
            "Category",
            "CategorySlug",
            "DistinctGeneCount",
            "MinDistinctGenesRequired",
            "MinDistinctGeneFraction",
            "GenomeCountWithAnyHits",
            "GenomeCount",
            "CutoffTag",
            "HeatmapFileName",
            "TreeFileName",
            "PhylumAnnotationFileName",
            "GenusAnnotationFileName",
            "SourceStripFileName",
            "GenomeNameLabelsFileName",
            "CompletenessPiechartFileName",
        ],
    )

    readme_lines = [
        f"Distinct-gene filter fraction: {min_fraction}",
        (
            "Rule: keep a genome only when the number of distinct genes present in a "
            "category is >= ceil(distinct category gene count * fraction)."
        ),
        f"Bundle suffix: {bundle_suffix}",
        f"Categories with at least one retained genome: {len(manifest_rows)}",
        f"Categories with zero retained genomes: {len(empty_specs)}",
        "",
        "Created bundles:",
    ]
    for row in manifest_rows:
        readme_lines.append(
            f"- {row['Category']}: kept {row['GenomeCount']} of {row['GenomeCountWithAnyHits']} genomes "
            f"using cutoff {row['CutoffTag']}"
        )
    if empty_specs:
        readme_lines.extend(["", "Categories with zero retained genomes:"])
        for spec in empty_specs:
            readme_lines.append(
                f"- {spec['category']}: cutoff {spec['cutoff_tag']} removed all "
                f"{len(spec['genomes_with_any_hits'])} genomes with any distinct-gene hits"
            )
    (output_dir / readme_name).write_text("\n".join(readme_lines) + "\n", encoding="utf-8")

    return {
        "bundle_suffix": bundle_suffix,
        "heatmap_dir_name": heatmap_dir.name,
        "tree_dir_name": tree_dir.name,
        "annotation_dir_name": annotation_dir.name,
        "manifest_name": manifest_name,
        "readme_name": readme_name,
        "bundle_count": len(manifest_rows),
        "empty_category_count": len(empty_specs),
    }


def write_itol_labels(
    path: Path,
    genomes: list[str],
    metadata_lookup: dict[str, dict[str, str]],
    label_field: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("LABELS\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATA\n")
        for genome in genomes:
            label = normalize_text(metadata_lookup[genome].get(label_field, "")) or genome
            handle.write(f"{genome}\t{label}\n")


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.filtered_category_min_distinct_fraction <= 1.0:
        raise ValueError("--filtered-category-min-distinct-fraction must be between 0 and 1")
    gene_info_path = Path(args.gene_info_xlsx).resolve()
    ko_output_path = Path(args.ko_output).resolve()
    ko_scored_output_path = Path(args.ko_scored_output).resolve()
    metadata_path = Path(args.metadata_xlsx).resolve()
    tree_path = Path(args.tree).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    (
        gene_order,
        gene_meta,
        by_uniprot,
        by_gene,
        category_gene_row_order,
        _subcategory_gene_row_order,
    ) = read_gene_info(gene_info_path)
    tree_labels = read_tree_labels(tree_path)
    metadata_header, metadata_lookup = read_metadata(metadata_path)

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
        raw_total_rows,
        scored_total_rows,
        besthit_by_hit,
        filtered_low_identity_rows,
        filtered_low_coverage_rows,
    ) = process_scored_ko_output(
        ko_output_path,
        ko_scored_output_path,
        by_uniprot,
        by_gene,
        min_pident=args.min_pident,
        min_qcovs=args.min_qcovs,
    )

    ko_genomes = sorted(set(hit_to_genome.values()))
    excluded_tree_genomes = sorted(set(ko_genomes).difference(tree_labels))
    excluded_hit_proteins = 0
    excluded_raw_hit_rows = 0
    if excluded_tree_genomes:
        allowed_tree_labels = set(tree_labels)
        excluded_raw_hit_rows = sum(
            count for genome, count in genome_hit_counter.items() if genome not in allowed_tree_labels
        )
        original_mapped_hit_ids = set(hit_to_genes)
        filtered_hit_ids = {
            hit_id for hit_id, genome_id in hit_to_genome.items() if genome_id in allowed_tree_labels
        }
        excluded_hit_proteins = sum(
            1 for hit_id in original_mapped_hit_ids if hit_to_genome.get(hit_id) not in allowed_tree_labels
        )
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

    counts = count_hits(
        tree_labels,
        gene_order,
        gene_meta,
        hit_to_genes,
        hit_to_genome,
        hit_to_queries,
        query_cache,
        besthit_by_hit,
    )

    gene_summary_meta = {
        gene: {
            "Category": gene_meta[gene].get("Category", ""),
            "Sub category": gene_meta[gene].get("Sub category", ""),
            "Gene name": gene_meta[gene].get("Gene name", ""),
            "KO": gene_meta[gene].get("KO", ""),
            "EC": gene_meta[gene].get("EC", ""),
        }
        for gene in counts["gene_order"]
    }
    category_summary_meta = {
        category: {"Category": category, "Sub category": "", "Gene name": "", "KO": "", "EC": ""}
        for category in counts["category_order"]
    }
    subcategory_to_categories: dict[str, set[str]] = defaultdict(set)
    for gene in counts["gene_order"]:
        categories = split_multi_value(gene_meta[gene].get("Category", ""))
        subcategories = split_multi_value(gene_meta[gene].get("Sub category", ""))
        for sub in subcategories:
            subcategory_to_categories[sub].update(categories)
    subcategory_summary_meta = {
        sub: {
            "Category": "; ".join(sorted(subcategory_to_categories[sub])),
            "Sub category": sub,
            "Gene name": "",
            "KO": "",
            "EC": "",
        }
        for sub in counts["subcategory_order"]
    }

    gene_summary_rows = build_summary_rows(
        counts["gene_order"],
        tree_labels,
        counts["besthit_gene_counts"],
        counts["strict_gene_counts"],
        counts["inclusive_gene_counts"],
        gene_summary_meta,
        "Gene",
    )
    category_summary_rows = build_summary_rows(
        counts["category_order"],
        tree_labels,
        counts["besthit_category_counts"],
        counts["strict_category_counts"],
        counts["inclusive_category_counts"],
        category_summary_meta,
        "Category",
    )
    subcategory_summary_rows = build_summary_rows(
        counts["subcategory_order"],
        tree_labels,
        counts["besthit_subcategory_counts"],
        counts["strict_subcategory_counts"],
        counts["inclusive_subcategory_counts"],
        subcategory_summary_meta,
        "Sub_category",
    )

    write_matrix(
        output_dir / "gene_counts_besthit.tsv",
        tree_labels,
        counts["gene_order"],
        counts["besthit_gene_counts"],
    )
    write_matrix(
        output_dir / "gene_counts_strict.tsv",
        tree_labels,
        counts["gene_order"],
        counts["strict_gene_counts"],
    )
    write_matrix(
        output_dir / "gene_counts_inclusive.tsv",
        tree_labels,
        counts["gene_order"],
        counts["inclusive_gene_counts"],
    )
    write_matrix(
        output_dir / "category_counts_besthit.tsv",
        tree_labels,
        counts["category_order"],
        counts["besthit_category_counts"],
    )
    write_matrix(
        output_dir / "category_counts_strict.tsv",
        tree_labels,
        counts["category_order"],
        counts["strict_category_counts"],
    )
    write_matrix(
        output_dir / "category_counts_inclusive.tsv",
        tree_labels,
        counts["category_order"],
        counts["inclusive_category_counts"],
    )
    write_matrix(
        output_dir / "subcategory_counts_besthit.tsv",
        tree_labels,
        counts["subcategory_order"],
        counts["besthit_subcategory_counts"],
    )
    write_matrix(
        output_dir / "subcategory_counts_strict.tsv",
        tree_labels,
        counts["subcategory_order"],
        counts["strict_subcategory_counts"],
    )
    write_matrix(
        output_dir / "subcategory_counts_inclusive.tsv",
        tree_labels,
        counts["subcategory_order"],
        counts["inclusive_subcategory_counts"],
    )

    summary_fields = [
        "Level",
        "Label",
        "Category",
        "Sub_category",
        "Gene_name",
        "KO",
        "EC",
        "Genome_count",
        "Besthit_total",
        "Besthit_genomes_with_hits",
        "Besthit_average_per_genome",
        "Besthit_average_if_present",
        "Besthit_max_in_single_genome",
        "Strict_total",
        "Strict_genomes_with_hits",
        "Strict_average_per_genome",
        "Strict_average_if_present",
        "Strict_max_in_single_genome",
        "Inclusive_total",
        "Inclusive_genomes_with_hits",
        "Inclusive_average_per_genome",
        "Inclusive_average_if_present",
        "Inclusive_max_in_single_genome",
    ]
    write_rows(output_dir / "gene_average_summary.tsv", gene_summary_rows, summary_fields)
    write_rows(output_dir / "category_average_summary.tsv", category_summary_rows, summary_fields)
    write_rows(
        output_dir / "subcategory_average_summary.tsv",
        subcategory_summary_rows,
        summary_fields,
    )

    query_mapping_rows = []
    for query_id in sorted(query_cache):
        mapping = query_cache[query_id]
        query_mapping_rows.append(
            {
                "query_id": query_id,
                "match_mode": str(mapping["match_mode"]),
                "mapped_gene_count": str(len(mapping["genes"])),
                "mapped_genes": "; ".join(mapping["genes"]),
                "mapped_categories": "; ".join(mapping["categories"]),
                "mapped_sub_categories": "; ".join(mapping["sub_categories"]),
            }
        )
    write_rows(
        output_dir / "query_mapping_summary.tsv",
        query_mapping_rows,
        [
            "query_id",
            "match_mode",
            "mapped_gene_count",
            "mapped_genes",
            "mapped_categories",
            "mapped_sub_categories",
        ],
    )

    unmapped_queries = [row for row in query_mapping_rows if row["match_mode"] == "unmapped"]
    write_rows(
        output_dir / "unmapped_queries.tsv",
        unmapped_queries,
        [
            "query_id",
            "match_mode",
            "mapped_gene_count",
            "mapped_genes",
            "mapped_categories",
            "mapped_sub_categories",
        ],
    )

    write_rows(
        output_dir / "besthit_assigned_hit_proteins.tsv",
        counts["besthit_assignment_rows"],
        [
            "hit_protein_id",
            "genome_id",
            "chosen_query_id",
            "chosen_gene",
            "chosen_category",
            "chosen_sub_category",
            "chosen_match_mode",
            "bitscore",
            "evalue",
            "pident",
            "qcovs",
            "alignment_length",
            "candidate_gene_count",
            "candidate_genes",
            "candidate_categories",
            "candidate_sub_categories",
            "chosen_query_mapped",
        ],
    )
    write_rows(
        output_dir / "besthit_resolved_ambiguous_hit_proteins.tsv",
        counts["besthit_resolved_ambiguous_rows"],
        [
            "hit_protein_id",
            "genome_id",
            "chosen_query_id",
            "chosen_gene",
            "chosen_category",
            "chosen_sub_category",
            "chosen_match_mode",
            "bitscore",
            "evalue",
            "pident",
            "qcovs",
            "alignment_length",
            "candidate_gene_count",
            "candidate_genes",
            "candidate_categories",
            "candidate_sub_categories",
            "chosen_query_mapped",
        ],
    )

    ambiguous_hit_rows = []
    for row in counts["ambiguous_hits"]:
        row = dict(row)
        row["query_count"] = str(len(hit_to_queries[row["hit_protein_id"]]))
        row["queries"] = "; ".join(sorted(hit_to_queries[row["hit_protein_id"]]))
        ambiguous_hit_rows.append(row)
    write_rows(
        output_dir / "ambiguous_hit_proteins.tsv",
        ambiguous_hit_rows,
        [
            "hit_protein_id",
            "genome_id",
            "mapped_gene_count",
            "mapped_genes",
            "mapped_categories",
            "mapped_sub_categories",
            "query_count",
            "queries",
        ],
    )

    genome_metadata_rows = []
    for genome in tree_labels:
        row = {field: metadata_lookup[genome].get(field, "") for field in metadata_header}
        row["Tree_label"] = genome
        row["KO_raw_hit_rows"] = str(genome_hit_counter.get(genome, 0))
        row["Besthit_unique_gene_hits"] = str(sum(counts["besthit_gene_counts"][genome].values()))
        row["Strict_unique_gene_hits"] = str(sum(counts["strict_gene_counts"][genome].values()))
        row["Inclusive_unique_gene_hits"] = str(
            sum(counts["inclusive_gene_counts"][genome].values())
        )
        genome_metadata_rows.append(row)
    write_rows(
        output_dir / "genome_metadata_tree_order.tsv",
        genome_metadata_rows,
        [
            "Tree_label",
            "KO_raw_hit_rows",
            "Besthit_unique_gene_hits",
            "Strict_unique_gene_hits",
            "Inclusive_unique_gene_hits",
        ]
        + metadata_header,
    )

    write_itol_heatmap(
        output_dir / "itol_gene_counts_besthit.txt",
        "gene_counts_besthit",
        tree_labels,
        counts["gene_order"],
        counts["besthit_gene_counts"],
        show_labels=False,
    )
    write_itol_heatmap(
        output_dir / "itol_gene_counts_strict.txt",
        "gene_counts_strict",
        tree_labels,
        counts["gene_order"],
        counts["strict_gene_counts"],
        show_labels=False,
    )
    write_itol_heatmap(
        output_dir / "itol_gene_counts_inclusive.txt",
        "gene_counts_inclusive",
        tree_labels,
        counts["gene_order"],
        counts["inclusive_gene_counts"],
        show_labels=False,
    )
    write_itol_heatmap(
        output_dir / "itol_subcategory_counts_besthit.txt",
        "subcategory_counts_besthit",
        tree_labels,
        counts["subcategory_order"],
        counts["besthit_subcategory_counts"],
        show_labels=True,
    )
    write_itol_heatmap(
        output_dir / "itol_subcategory_counts_strict.txt",
        "subcategory_counts_strict",
        tree_labels,
        counts["subcategory_order"],
        counts["strict_subcategory_counts"],
        show_labels=True,
    )
    write_itol_heatmap(
        output_dir / "itol_category_counts_besthit.txt",
        "category_counts_besthit",
        tree_labels,
        counts["category_order"],
        counts["besthit_category_counts"],
        show_labels=True,
    )
    write_itol_heatmap(
        output_dir / "itol_category_counts_strict.txt",
        "category_counts_strict",
        tree_labels,
        counts["category_order"],
        counts["strict_category_counts"],
        show_labels=True,
    )

    phylum_colors = make_phylum_color_map(
        [
            normalize_text(metadata_lookup[genome].get("Phylum", "")) or "Unassigned"
            for genome in tree_labels
        ],
    )
    genus_colors = make_color_map(
        [
            normalize_text(metadata_lookup[genome].get("Genus", "")) or "Unassigned"
            for genome in tree_labels
        ],
        prefer_high_contrast=True,
    )
    write_itol_tree_colors(
        output_dir / "itol_phylum_annotation.txt",
        tree_path,
        metadata_lookup,
        "Phylum",
        phylum_colors,
    )
    write_itol_tree_colors(
        output_dir / "itol_phylum_colorstrip.txt",
        tree_path,
        metadata_lookup,
        "Phylum",
        phylum_colors,
    )
    write_itol_tree_colors(
        output_dir / "itol_genus_annotation.txt",
        tree_path,
        metadata_lookup,
        "Genus",
        genus_colors,
    )
    write_itol_tree_colors(
        output_dir / "itol_genus_colorstrip.txt",
        tree_path,
        metadata_lookup,
        "Genus",
        genus_colors,
    )
    source_strip_available = "Source" in metadata_header
    if source_strip_available:
        source_colors = write_itol_colorstrip(
            output_dir / "itol_source_colorstrip.txt",
            "Genome source",
            tree_labels,
            metadata_lookup,
            "Source",
            strip_width=28,
        )
    write_taxonomy_legend(output_dir / "itol_phylum_legend.tsv", "Phylum", phylum_colors)
    write_taxonomy_legend(output_dir / "itol_genus_legend.tsv", "Genus", genus_colors)
    if source_strip_available:
        write_taxonomy_legend(output_dir / "itol_source_legend.tsv", "Source", source_colors)
    completeness_available = "Completeness (%)" in metadata_header
    if completeness_available:
        write_itol_piechart(
            output_dir / "itol_completeness_piechart.txt",
            "Genome completeness (%)",
            tree_labels,
            metadata_lookup,
            "Completeness (%)",
        )
    write_itol_labels(
        output_dir / "itol_genome_name_labels.txt",
        tree_labels,
        metadata_lookup,
        "Genome name",
    )
    category_gene_dir = output_dir / "itol_gene_counts_besthit_by_category"
    category_gene_manifest = write_category_gene_heatmaps(
        category_gene_dir,
        tree_labels,
        counts["gene_order"],
        gene_meta,
        category_gene_row_order,
        counts["besthit_gene_counts"],
        counts["category_order"],
    )
    category_tree_dir = output_dir / "itol_gene_counts_besthit_by_category_trees"
    category_tree_manifest = write_category_pruned_trees(
        category_tree_dir,
        tree_path,
        tree_labels,
        counts["besthit_category_counts"],
        counts["category_order"],
    )
    category_tree_by_name = {row["Category"]: row for row in category_tree_manifest}
    category_annotation_dir = output_dir / "itol_gene_counts_besthit_by_category_annotations"
    category_annotation_manifest = write_category_annotation_files(
        category_annotation_dir,
        category_tree_dir,
        tree_path.stem,
        tree_labels,
        metadata_lookup,
        counts["besthit_category_counts"],
        counts["category_order"],
        phylum_colors,
        genus_colors,
        source_colors if source_strip_available else None,
        completeness_available,
    )
    category_annotation_by_name = {row["Category"]: row for row in category_annotation_manifest}
    category_manifest_rows = []
    for row in category_gene_manifest:
        merged = dict(row)
        merged.update(category_tree_by_name.get(row["Category"], {}))
        merged.update(category_annotation_by_name.get(row["Category"], {}))
        category_manifest_rows.append(merged)
    write_rows(
        output_dir / "itol_gene_counts_besthit_by_category.tsv",
        category_manifest_rows,
        [
            "Category",
            "HeatmapFileName",
            "GeneCount",
            "TreeFileName",
            "GenomeCount",
            "PhylumAnnotationFileName",
            "GenusAnnotationFileName",
            "SourceStripFileName",
            "GenomeNameLabelsFileName",
            "CompletenessPiechartFileName",
        ],
    )

    filtered_category_bundle_info = None
    if args.filtered_category_min_distinct_fraction > 0:
        filtered_category_bundle_info = write_filtered_category_distinct_gene_bundles(
            output_dir,
            tree_path,
            tree_labels,
            counts["gene_order"],
            gene_meta,
            category_gene_row_order,
            counts["besthit_gene_counts"],
            counts["category_order"],
            metadata_lookup,
            phylum_colors,
            genus_colors,
            source_colors if source_strip_available else None,
            args.filtered_category_min_distinct_fraction,
            completeness_available,
        )

    summary_lines = [
        f"Tree genomes: {len(tree_labels)}",
        f"KO raw rows: {raw_total_rows}",
        f"Identity/coverage floor applied: min_pident={args.min_pident}, min_qcovs={args.min_qcovs}",
        f"Rows dropped for pident < {args.min_pident}: {filtered_low_identity_rows}",
        f"Rows dropped for qcovs < {args.min_qcovs}: {filtered_low_coverage_rows}",
        f"KO scored rows used for best-hit resolution (post identity/coverage filter): {scored_total_rows}",
        f"Unique queries in KO output: {len(query_cache)}",
        f"Mapped queries: {sum(1 for row in query_cache.values() if row['genes'])}",
        f"Unmapped queries: {sum(1 for row in query_cache.values() if not row['genes'])}",
        f"Unique hit proteins with >=1 mapped gene: {len(hit_to_genes)}",
        f"Excluded genomes absent from tree: {len(excluded_tree_genomes)}",
        (
            "Excluded genome labels: " + ", ".join(excluded_tree_genomes)
            if excluded_tree_genomes
            else "Excluded genome labels: None"
        ),
        f"Excluded KO raw rows from pruned genomes: {excluded_raw_hit_rows}",
        f"Excluded hit proteins from pruned genomes: {excluded_hit_proteins}",
        f"Best-hit assigned proteins: {counts['besthit_hit_total']}",
        f"Best-hit proteins whose chosen query is unmapped: {counts['besthit_unmapped_total']}",
        f"Strict unique hit proteins (1 mapped gene): {sum(1 for genes in hit_to_genes.values() if len(genes) == 1)}",
        f"Ambiguous hit proteins (>1 mapped gene): {sum(1 for genes in hit_to_genes.values() if len(genes) > 1)}",
        f"Best-hit gene-hit total: {counts['besthit_hit_total']}",
        f"Strict gene-hit total: {counts['strict_hit_total']}",
        f"Inclusive gene-hit total: {counts['inclusive_hit_total']}",
        f"Gene columns: {len(counts['gene_order'])}",
        f"Category columns: {len(counts['category_order'])}",
        f"Sub-category columns: {len(counts['subcategory_order'])}",
        "",
        "Recommended iTOL uploads:",
        f"- Tree: {tree_path}",
        f"- Phylum annotation: {output_dir / 'itol_phylum_annotation.txt'}",
        f"- Genus annotation: {output_dir / 'itol_genus_annotation.txt'}",
    ]
    if source_strip_available:
        summary_lines.append(f"- Source color strip: {output_dir / 'itol_source_colorstrip.txt'}")
    if completeness_available:
        summary_lines.append(
            f"- Genome completeness pie chart: {output_dir / 'itol_completeness_piechart.txt'}"
        )
    summary_lines.extend([
        f"- Genome-name labels: {output_dir / 'itol_genome_name_labels.txt'}",
        f"- Category heatmap: {output_dir / 'itol_category_counts_besthit.txt'}",
        f"- Sub-category heatmap: {output_dir / 'itol_subcategory_counts_besthit.txt'}",
        f"- Gene heatmap: {output_dir / 'itol_gene_counts_besthit.txt'}",
        f"- Gene heatmaps split by category: {output_dir / 'itol_gene_counts_besthit_by_category'}",
        f"- Category-pruned trees: {output_dir / 'itol_gene_counts_besthit_by_category_trees'}",
        f"- Category-pruned annotation files: {output_dir / 'itol_gene_counts_besthit_by_category_annotations'}",
        f"- Category heatmap manifest: {output_dir / 'itol_gene_counts_besthit_by_category.tsv'}",
        "",
        "Interpretation note:",
        "- Besthit counts keep sole-gene proteins as-is and assign ambiguous proteins to the single top-scoring KO query from the targeted scored BLAST table.",
        "- Strict counts keep only proteins that mapped to exactly one gene definition.",
        "- Inclusive counts let one protein contribute to every gene it matched and are therefore more permissive.",
    ])
    if filtered_category_bundle_info is not None:
        summary_lines.extend(
            [
                "",
                "Distinct-gene filtered category bundles:",
                (
                    f"- Fraction: {args.filtered_category_min_distinct_fraction} "
                    "(keep genomes with >= ceil(distinct category genes * fraction))"
                ),
                f"- Filtered heatmaps: {output_dir / filtered_category_bundle_info['heatmap_dir_name']}",
                f"- Filtered trees: {output_dir / filtered_category_bundle_info['tree_dir_name']}",
                f"- Filtered annotations: {output_dir / filtered_category_bundle_info['annotation_dir_name']}",
                f"- Filtered manifest: {output_dir / filtered_category_bundle_info['manifest_name']}",
                f"- Filtered README: {output_dir / filtered_category_bundle_info['readme_name']}",
            ]
        )
    (output_dir / "README_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"Output directory: {output_dir}")
    output_names = [
        "gene_average_summary.tsv",
        "category_average_summary.tsv",
        "subcategory_average_summary.tsv",
        "itol_phylum_annotation.txt",
        "itol_genus_annotation.txt",
        "itol_genome_name_labels.txt",
        "itol_category_counts_besthit.txt",
        "itol_subcategory_counts_besthit.txt",
        "itol_gene_counts_besthit.txt",
        "itol_gene_counts_besthit_by_category.tsv",
        "itol_category_counts_strict.txt",
        "itol_subcategory_counts_strict.txt",
        "itol_gene_counts_strict.txt",
        "README_summary.txt",
    ]
    if filtered_category_bundle_info is not None:
        output_names.extend(
            [
                str(filtered_category_bundle_info["manifest_name"]),
                str(filtered_category_bundle_info["readme_name"]),
            ]
        )
    if source_strip_available:
        output_names.insert(5, "itol_source_colorstrip.txt")
    for name in output_names:
        print(name)


if __name__ == "__main__":
    main()
