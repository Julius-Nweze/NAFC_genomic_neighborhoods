#!/usr/bin/env python3
# Title          : export_reference_panel_cluster_reblast_targets.py
# Description    : Export unnamed, generic or hypothetical neighborhood genes as Swiss-Prot re-BLAST targets, per genome-source tier
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/31
# Usage          : python3 export_reference_panel_cluster_reblast_targets.py [options]

"""Export unnamed, generic or hypothetical neighborhood genes as Swiss-Prot re-BLAST targets, per genome-source tier."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

GENERIC_PRODUCT_RE = re.compile(
    r"^(hypothetical protein|uncharacterized protein|putative protein|domain[- ]containing protein|"
    r"duf\d+.*|conserved protein)$",
    re.IGNORECASE,
)


def load_source_map(path: Path) -> dict[str, str]:
    m: dict[str, str] = {}
    started = False
    for line in path.read_text().splitlines():
        if line.strip() == "DATA":
            started = True
            continue
        if started and line.strip():
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                m[parts[0]] = parts[2]
    return m


def normalize_genome_id(g: str) -> str:
    g = g.strip()
    g = re.sub(r"\.\d+$", "", g)
    return g.replace(".", "_")


def is_generic(gene: str, product: str) -> bool:
    gene = (gene or "").strip()
    product = (product or "").strip()
    if not gene:
        return True
    return bool(GENERIC_PRODUCT_RE.match(product))


def refined_name_for(product: str) -> str:
    product = (product or "hypothetical protein").strip()
    return product.lower().replace(" ", "_")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--members",
        default="/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters/NAFC_gene_panel_neighborhood_members.tsv",
    )
    parser.add_argument(
        "--source-colorstrip",
        default="/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_source_colorstrip.txt",
    )
    parser.add_argument(
        "--blast-db",
        default="/path/to/your/Ref_genomes/Analysis/Database/Genomes_MAGs.faa",
    )
    parser.add_argument(
        "--out-dir",
        default="/path/to/your/Ref_genomes/Analysis/Result/Swissprot_reexport",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    src_raw = load_source_map(Path(args.source_colorstrip))
    norm_src: dict[str, str] = {}
    for k, v in src_raw.items():
        norm_src[normalize_genome_id(re.sub(r"^g_", "", k))] = v
        norm_src[normalize_genome_id(k)] = v

    tier_label = {
        "Plant roots": "MAGs_batch",
        "Genome bank": "Ref_batch",
        "Soil-containing OSPW": "Colla_batch",
    }

    rows_by_tier: dict[str, list[dict[str, str]]] = {label: [] for label in tier_label.values()}

    with Path(args.members).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["Is panel seed hit"] == "Yes":
                continue
            gene = row.get("Gene", "")
            product = row.get("Product", "")
            if not is_generic(gene, product):
                continue
            src = norm_src.get(normalize_genome_id(row["Genome"]))
            label = tier_label.get(src)
            if label is None:
                continue
            flags = "no_gene;generic_or_hypothetical;unresolved_function" if not gene.strip() else "generic_or_hypothetical"
            rows_by_tier[label].append(
                {
                    "Genome": row["Genome"],
                    "Cluster ID": row["Cluster ID"],
                    "Locus tag": row["Locus tag"],
                    "Original gene": gene,
                    "Refined name": refined_name_for(product),
                    "Original product": product,
                    "Enzyme class": row.get("Enzyme class", ""),
                    "Ambiguity flags": flags,
                    "Confidence": "low",
                }
            )

    for label, rows in rows_by_tier.items():
        # de-duplicate by locus tag (a gene can appear in >1 cluster as expanded context)
        seen: dict[str, dict[str, str]] = {}
        for r in rows:
            seen.setdefault(r["Locus tag"], r)
        rows = list(seen.values())

        tsv_path = out_dir / f"{label}_reblast_targets.tsv"
        with tsv_path.open("w", newline="", encoding="utf-8") as handle:
            fieldnames = ["Genome", "Cluster ID", "Locus tag", "Original gene", "Refined name",
                          "Original product", "Enzyme class", "Ambiguity flags", "Confidence"]
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        entry_batch = out_dir / f"{label}_entry_batch.txt"
        with entry_batch.open("w") as handle:
            for r in rows:
                handle.write(f"{r['Locus tag']}__{r['Genome']}\n")

        extracted = out_dir / f"{label}_reblast_targets_raw.faa"
        subprocess.run(
            [
                "blastdbcmd", "-db", args.blast_db,
                "-entry_batch", str(entry_batch),
                "-out", str(extracted),
                "-outfmt", "%f",
            ],
            check=True,
        )

        meta_by_key = {f"{r['Locus tag']}__{r['Genome']}": r for r in rows}
        faa_path = out_dir / f"{label}_reblast_targets.faa"
        with extracted.open() as fin, faa_path.open("w") as fout:
            header = None
            seq_lines: list[str] = []

            def flush():
                if header is None:
                    return
                key = header.split()[0][1:]
                meta = meta_by_key.get(key)
                if meta is None:
                    return
                locus = meta["Locus tag"]
                fout.write(
                    f">{locus} cluster={meta['Cluster ID']} refined={meta['Refined name']} "
                    f"flags={meta['Ambiguity flags']} product={meta['Original product'].replace(' ', '_')}\n"
                )
                fout.write("\n".join(seq_lines) + "\n")

            for line in fin:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    header = line
                    seq_lines = []
                else:
                    seq_lines.append(line)
            flush()

        print(f"{label}: {len(rows)} targets exported -> {tsv_path}, {faa_path}")


if __name__ == "__main__":
    main()
