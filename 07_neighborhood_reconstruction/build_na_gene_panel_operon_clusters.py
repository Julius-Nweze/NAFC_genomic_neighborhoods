#!/usr/bin/env python3
# Title          : build_na_gene_panel_operon_clusters.py
# Description    : Group best-hit panel genes into candidate genomic neighborhoods from Prokka coordinates (same contig, gap <=12 genes and <=25 kb)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/22
# Usage          : python3 build_na_gene_panel_operon_clusters.py [options]

"""Group best-hit panel genes into candidate genomic neighborhoods from Prokka coordinates (same contig, gap <=12 genes and <=25 kb)."""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

REF_CLUSTER_SCRIPT_DIR = Path(__file__).resolve().parent.parent / "shared_dependencies"
sys.path.insert(0, str(REF_CLUSTER_SCRIPT_DIR))
from generate_ref_cluster_outputs import (  # type: ignore  # noqa: E402
    enzyme_class,
    parse_gff_attributes,
    read_tsv,
    write_tsv,
)

ROW_GAP_MAX = 12
COORD_GAP_MAX = 25000

PROKKA_ROOT = Path("/path/to/your/Ref_genomes/Prokka")
PROKKA_FLAT_ROOT = Path("/path/to/your/Ref_genomes/Prokka/Prokka")
BESTHIT_TABLE = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/ITOL/"
    "gene_count_summary/besthit_assigned_hit_proteins.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters"
)


def resolve_prokka_folder(stem: str) -> Path | None:
    # Prokka/Prokka/ is a flat, non-tree-structured copy of the same Prokka output
    # (confirmed byte-identical to Prokka/{Ref,MAGs,Colla}/ across all 365 genomes);
    # it is checked first.
    # Ref/ and MAGs/ keep Prokka output directly in the genome folder; Colla/
    # nests it one level deeper under a "prokka" subfolder. Some genomes (mostly
    # Colla) use an older "<genome>_Prokka" naming convention instead.
    candidates = [
        PROKKA_FLAT_ROOT / stem,
        PROKKA_FLAT_ROOT / f"{stem}.fa",
        PROKKA_FLAT_ROOT / f"{stem}_Prokka",
        PROKKA_ROOT / "Ref" / f"{stem}.fa",
        PROKKA_ROOT / "MAGs" / stem,
        PROKKA_ROOT / "Colla" / stem / "prokka",
        PROKKA_ROOT / "Colla" / stem,
        PROKKA_ROOT / "Colla" / f"{stem}_Prokka",
        PROKKA_ROOT / "Ref" / f"{stem}_Prokka",
        PROKKA_ROOT / "MAGs" / f"{stem}_Prokka",
    ]
    for candidate in candidates:
        if candidate.is_dir() and (
            list(candidate.glob("PROKKA_*.gff")) or list(candidate.glob("*.gff"))
        ):
            return candidate
    return None


def locate_gff(folder: Path) -> Path:
    matches = sorted(folder.glob("PROKKA_*.gff")) or sorted(folder.glob("*.gff"))
    return matches[0]


def locate_annotation_tsv(folder: Path) -> Path | None:
    candidates = []
    for path in sorted(folder.glob("*.tsv")):
        name = path.name
        if name.startswith("NA_cluster_") or name.endswith("_cluster_assessment.tsv"):
            continue
        if name.endswith("_gene_maps.tsv") or name.endswith("_gene_refinement_layer.tsv"):
            continue
        if name.endswith("_reblast_targets.tsv"):
            continue
        if name.endswith("_cds.tsv"):
            continue
        candidates.append(path)
    if not candidates:
        return None
    prokka = [path for path in candidates if path.name.startswith("PROKKA_")]
    return (prokka or candidates)[0]


def read_gene_records(folder: Path) -> dict[str, dict[str, object]]:
    """Return locus_tag -> {contig, start, end, strand, row_index, gene, product, EC, COG}."""
    gff_path = locate_gff(folder)
    coords: list[dict[str, object]] = []
    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_gff_attributes(parts[8])
            locus = attrs.get("locus_tag") or attrs.get("ID")
            if not locus:
                continue
            coords.append(
                {
                    "locus_tag": locus,
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            )
    coords.sort(key=lambda row: (row["contig"], row["start"]))
    for index, row in enumerate(coords):
        row["row_index"] = index

    ann_path = locate_annotation_tsv(folder)
    ann_by_locus: dict[str, dict[str, str]] = {}
    if ann_path is not None:
        for row in read_tsv(ann_path):
            locus = row.get("locus_tag", "")
            if locus:
                ann_by_locus[locus] = row

    records: dict[str, dict[str, object]] = {}
    for row in coords:
        locus = row["locus_tag"]
        ann = ann_by_locus.get(locus, {})
        records[locus] = {
            **row,
            "gene": ann.get("gene", ""),
            "product": ann.get("product", ""),
            "EC_number": ann.get("EC_number", ""),
            "COG": ann.get("COG", ""),
        }
    return records


def group_seed_clusters(seed_rows: list[dict[str, object]]) -> list[list[dict[str, object]]]:
    ordered = sorted(seed_rows, key=lambda row: (row["contig"], row["start"]))
    clusters: list[list[dict[str, object]]] = []
    current = [ordered[0]]
    for row in ordered[1:]:
        prev = current[-1]
        same_contig = row["contig"] == prev["contig"]
        row_gap = row["row_index"] - prev["row_index"]
        coord_gap = abs(row["start"] - prev["end"])
        if same_contig and row_gap <= ROW_GAP_MAX and coord_gap <= COORD_GAP_MAX:
            current.append(row)
        else:
            clusters.append(current)
            current = [row]
    clusters.append(current)
    return clusters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--besthit-table", default=str(BESTHIT_TABLE))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    seed_rows_by_genome: dict[str, list[dict[str, object]]] = defaultdict(list)
    with open(args.besthit_table, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hit_id = row["hit_protein_id"]
            locus_tag, sep, stem = hit_id.rpartition("__")
            if not sep:
                continue
            seed_rows_by_genome[stem].append({**row, "locus_tag": locus_tag})

    genomes_without_prokka: list[dict[str, str]] = []
    genomes_with_locus_mismatch: list[dict[str, str]] = []
    cluster_summary_rows: list[dict[str, str]] = []
    member_rows: list[dict[str, str]] = []

    category_total_counter: Counter = Counter()
    category_multi_counter: Counter = Counter()

    cluster_counter = 0
    processed_genomes = 0

    for stem in sorted(seed_rows_by_genome):
        seed_rows = seed_rows_by_genome[stem]
        folder = resolve_prokka_folder(stem)
        if folder is None:
            genomes_without_prokka.append({"Genome": stem, "SeedHitCount": str(len(seed_rows))})
            continue

        gene_records = read_gene_records(folder)

        seeds_with_coords: list[dict[str, object]] = []
        unmatched = 0
        for row in seed_rows:
            rec = gene_records.get(row["locus_tag"])
            if rec is None:
                unmatched += 1
                continue
            seeds_with_coords.append({**rec, "seed_row": row})
        if unmatched:
            genomes_with_locus_mismatch.append(
                {
                    "Genome": stem,
                    "UnmatchedSeedHits": str(unmatched),
                    "TotalSeedHits": str(len(seed_rows)),
                }
            )
        if not seeds_with_coords:
            continue

        processed_genomes += 1
        clusters = group_seed_clusters(seeds_with_coords)

        for cluster in clusters:
            cluster_counter += 1
            cluster_id = f"{stem}_panelcluster_{cluster_counter}"
            contig = cluster[0]["contig"]
            span_start = min(int(gene["start"]) for gene in cluster)
            span_end = max(int(gene["end"]) for gene in cluster)

            span_members = [
                rec
                for rec in gene_records.values()
                if rec["contig"] == contig and rec["start"] >= span_start and rec["end"] <= span_end
            ]
            span_members.sort(key=lambda row: (row["start"], row["end"]))

            seed_locus_tags = {gene["locus_tag"] for gene in cluster}
            seed_by_locus = {gene["locus_tag"]: gene["seed_row"] for gene in cluster}
            strength = "multi_hit" if len(seed_locus_tags) >= 2 else "single_hit"

            enzyme_counts: Counter = Counter()
            seed_genes: list[str] = []
            seed_categories: list[str] = []
            seed_subcategories: list[str] = []

            for rec in span_members:
                is_seed = rec["locus_tag"] in seed_locus_tags
                cls = enzyme_class(rec.get("gene", ""), rec.get("product", ""))
                enzyme_counts[cls] += 1
                seed_row = seed_by_locus.get(rec["locus_tag"])
                member_rows.append(
                    {
                        "Cluster ID": cluster_id,
                        "Genome": stem,
                        "Locus tag": rec["locus_tag"],
                        "Contig": contig,
                        "Start": str(rec["start"]),
                        "End": str(rec["end"]),
                        "Strand": rec.get("strand", ""),
                        "Gene": rec.get("gene", ""),
                        "Product": rec.get("product", ""),
                        "Enzyme class": cls,
                        "Is panel seed hit": "Yes" if is_seed else "No",
                        "Panel gene": seed_row.get("chosen_gene", "") if seed_row else "",
                        "Panel category": seed_row.get("chosen_category", "") if seed_row else "",
                        "Panel sub-category": seed_row.get("chosen_sub_category", "") if seed_row else "",
                        "Panel percent identity": seed_row.get("pident", "") if seed_row else "",
                        "Panel query coverage": seed_row.get("qcovs", "") if seed_row else "",
                        "Panel bit score": seed_row.get("bitscore", "") if seed_row else "",
                    }
                )
                if is_seed and seed_row:
                    seed_genes.append(seed_row.get("chosen_gene", ""))
                    categories = [c for c in seed_row.get("chosen_category", "").split("; ") if c]
                    seed_categories.extend(categories)
                    seed_subcategories.append(seed_row.get("chosen_sub_category", ""))
                    for cat in categories:
                        category_total_counter[cat] += 1
                        if strength == "multi_hit":
                            category_multi_counter[cat] += 1

            enzyme_summary = "; ".join(f"{cls}:{count}" for cls, count in enzyme_counts.most_common())
            cluster_summary_rows.append(
                {
                    "Cluster ID": cluster_id,
                    "Genome": stem,
                    "Cluster strength": strength,
                    "Contig": contig,
                    "Start": str(span_start),
                    "End": str(span_end),
                    "Span (bp)": str(span_end - span_start),
                    "Panel seed hit count": str(len(seed_locus_tags)),
                    "Total CDS in span": str(len(span_members)),
                    "Panel genes": "; ".join(sorted(set(g for g in seed_genes if g))),
                    "Panel categories": "; ".join(sorted(set(seed_categories))),
                    "Panel sub-categories": "; ".join(sorted(set(c for c in seed_subcategories if c))),
                    "Enzyme class counts": enzyme_summary,
                }
            )

    summary_fields = [
        "Cluster ID", "Genome", "Cluster strength", "Contig", "Start", "End", "Span (bp)",
        "Panel seed hit count", "Total CDS in span", "Panel genes", "Panel categories",
        "Panel sub-categories", "Enzyme class counts",
    ]
    member_fields = [
        "Cluster ID", "Genome", "Locus tag", "Contig", "Start", "End", "Strand", "Gene", "Product",
        "Enzyme class", "Is panel seed hit", "Panel gene", "Panel category", "Panel sub-category",
        "Panel percent identity", "Panel query coverage", "Panel bit score",
    ]

    write_tsv(output_dir / "NA_gene_panel_clusters_summary.tsv", summary_fields, cluster_summary_rows)
    write_tsv(output_dir / "NA_gene_panel_cluster_members.tsv", member_fields, member_rows)
    write_tsv(
        output_dir / "genomes_without_prokka_annotation.tsv",
        ["Genome", "SeedHitCount"],
        genomes_without_prokka,
    )
    if genomes_with_locus_mismatch:
        write_tsv(
            output_dir / "genomes_with_unmatched_seed_hits.tsv",
            ["Genome", "UnmatchedSeedHits", "TotalSeedHits"],
            genomes_with_locus_mismatch,
        )

    total_seed_hits = sum(len(v) for v in seed_rows_by_genome.values())
    hits_without_prokka = sum(int(row["SeedHitCount"]) for row in genomes_without_prokka)
    multi_hit_clusters = [row for row in cluster_summary_rows if row["Cluster strength"] == "multi_hit"]
    single_hit_clusters = [row for row in cluster_summary_rows if row["Cluster strength"] == "single_hit"]
    seed_hits_in_multi = sum(int(row["Panel seed hit count"]) for row in multi_hit_clusters)
    seed_hits_in_single = sum(int(row["Panel seed hit count"]) for row in single_hit_clusters)

    lines = [
        f"Total panel seed hits (best-hit table): {total_seed_hits}",
        f"Genomes with at least one seed hit: {len(seed_rows_by_genome)}",
        f"Genomes processed (Prokka annotation found): {processed_genomes}",
        f"Genomes skipped (no Prokka annotation): {len(genomes_without_prokka)}",
        f"Seed hits affected by missing Prokka annotation: {hits_without_prokka}",
        "",
        f"Clusters found: {len(cluster_summary_rows)}",
        f"  Multi-hit clusters (>=2 panel genes co-occurring): {len(multi_hit_clusters)}",
        f"  Single-hit loci (1 panel gene + expanded genomic context): {len(single_hit_clusters)}",
        f"Seed hits inside multi-hit clusters: {seed_hits_in_multi}",
        f"Seed hits inside single-hit loci: {seed_hits_in_single}",
        "",
        "Per-category seed-hit distribution (counted per seed hit, not per cluster):",
    ]
    for cat, total in category_total_counter.most_common():
        multi = category_multi_counter.get(cat, 0)
        pct = (100.0 * multi / total) if total else 0.0
        lines.append(f"  {cat}: {total} seed hits, {multi} ({pct:.1f}%) in multi-hit clusters")

    if genomes_without_prokka:
        lines.append("")
        lines.append("Genomes skipped for lacking Prokka annotation:")
        for row in genomes_without_prokka:
            lines.append(f"  {row['Genome']} ({row['SeedHitCount']} seed hits not assessed)")

    (output_dir / "README_operon_cluster_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\nOutputs written to {output_dir}")


if __name__ == "__main__":
    main()
