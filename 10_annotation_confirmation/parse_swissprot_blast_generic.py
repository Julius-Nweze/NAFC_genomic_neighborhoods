#!/usr/bin/env python3
# Title          : parse_swissprot_blast_generic.py
# Description    : Parse Swiss-Prot BLASTP results and assign annotation-confirmation actions
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/24
# Usage          : python3 parse_swissprot_blast_generic.py --targets <tsv> --blast <tsv> --output-dir <dir> --prefix <name>

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


BLAST_FIELDS = [
    "qseqid",
    "sacc",
    "stitle",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
]

META_FIELDS = [
    "Genome",
    "Cluster ID",
    "Locus tag",
    "Original gene",
    "Refined name",
    "Original product",
    "Enzyme class",
    "Ambiguity flags",
    "Confidence",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in rows])


def parse_blast(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if len(parts) != len(BLAST_FIELDS):
                continue
            rows.append(dict(zip(BLAST_FIELDS, parts)))
    return rows


def protein_name(stitle: str) -> str:
    return re.split(r"\sOS=", stitle or "", maxsplit=1)[0].strip()


def gene_name(stitle: str) -> str:
    match = re.search(r"\sGN=([A-Za-z0-9_.-]+)", stitle or "")
    return match.group(1) if match else ""


def is_uncharacterized(name: str) -> bool:
    return bool(re.search(r"\b(uncharacterized|hypothetical|unknown function|putative protein)\b", name or "", re.I))


def as_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def classify_hit(row: dict[str, str]) -> str:
    pname = protein_name(row["stitle"])
    pident = as_float(row["pident"])
    qcov = as_float(row["qcovs"])
    evalue = as_float(row["evalue"])
    if is_uncharacterized(pname):
        return "low_uncharacterized_hit"
    if evalue <= 1e-50 and pident >= 45 and qcov >= 75:
        return "high"
    if evalue <= 1e-20 and pident >= 35 and qcov >= 60:
        return "medium"
    if evalue <= 1e-5 and pident >= 25 and qcov >= 40:
        return "low"
    return "very_low"


def is_missing_or_generic(meta: dict[str, str]) -> bool:
    current_gene = (meta.get("Original gene") or "").strip()
    refined = (meta.get("Refined name") or "").strip().lower()
    flags = meta.get("Ambiguity flags", "")
    if not current_gene or current_gene.upper() == "HP":
        return True
    if "no_gene" in flags or "generic_or_hypothetical" in flags:
        return True
    return refined in {"hypothetical/accessory protein", "hypothetical protein", "uncharacterized protein"}


def suggested_action(meta: dict[str, str], hit: dict[str, str] | None, confidence: str) -> str:
    if hit is None:
        return "retain_flag_for_TrEMBL_or_HMM_review"
    if confidence in {"low_uncharacterized_hit", "very_low"}:
        return "do_not_rename_uncharacterized_or_weak_hit"
    if confidence == "low":
        return "review_low_confidence_hit"
    if gene_name(hit["stitle"]) and is_missing_or_generic(meta):
        return "accept_blast_supported_name"
    return "accept_blast_supported_product_only"


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse focused Swiss-Prot BLAST hits for legacy Stage 2 re-BLAST targets.")
    parser.add_argument("--targets", type=Path, required=True, help="*_batch_reblast_targets.tsv file.")
    parser.add_argument("--blast", type=Path, required=True, help="BLASTP output in tabular format 6 with 14 fields.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for parsed Swiss-Prot outputs.")
    parser.add_argument("--prefix", required=True, help="Dataset prefix, for example Colla_batch or Ref_batch.")
    parser.add_argument("--database-note", default="UniProt Swiss-Prot local BLAST database.", help="Free-text database note for the summary markdown.")
    args = parser.parse_args()

    targets = read_tsv(args.targets)
    meta_by_locus = {row["Locus tag"]: row for row in targets}
    hits = parse_blast(args.blast)
    hits_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in hits:
        hits_by_query[hit["qseqid"]].append(hit)
    for query_hits in hits_by_query.values():
        query_hits.sort(key=lambda row: (as_float(row["evalue"]), -as_float(row["bitscore"]), -as_float(row["qcovs"])))

    all_hit_rows: list[dict[str, str]] = []
    top_hit_rows: list[dict[str, str]] = []
    suggestion_rows: list[dict[str, str]] = []
    supported_rows: list[dict[str, str]] = []

    hit_fields = META_FIELDS + BLAST_FIELDS + ["blast_protein_name", "blast_gene_name", "top_hit_confidence"]
    suggestion_fields = META_FIELDS + [
        "blast_hit",
        "top_accession",
        "top_protein_name",
        "top_gene_name",
        "pident",
        "qcovs",
        "evalue",
        "bitscore",
        "top_hit_confidence",
        "suggested_gene",
        "suggested_product",
        "suggested_action",
    ]

    for hit in hits:
        meta = meta_by_locus.get(hit["qseqid"], {})
        confidence = classify_hit(hit)
        row = {**meta, **hit}
        row["blast_protein_name"] = protein_name(hit["stitle"])
        row["blast_gene_name"] = gene_name(hit["stitle"])
        row["top_hit_confidence"] = confidence
        all_hit_rows.append(row)

    for target in targets:
        locus = target["Locus tag"]
        top = hits_by_query.get(locus, [None])[0]
        if top is None:
            confidence = "no_hit"
            action = suggested_action(target, None, confidence)
            row = {**target}
            row.update(
                {
                    "blast_hit": "No",
                    "top_accession": "",
                    "top_protein_name": "",
                    "top_gene_name": "",
                    "pident": "",
                    "qcovs": "",
                    "evalue": "",
                    "bitscore": "",
                    "top_hit_confidence": confidence,
                    "suggested_gene": target.get("Refined name", ""),
                    "suggested_product": target.get("Original product", ""),
                    "suggested_action": action,
                }
            )
        else:
            confidence = classify_hit(top)
            action = suggested_action(target, top, confidence)
            blast_gene = gene_name(top["stitle"])
            blast_product = protein_name(top["stitle"])
            row = {**target}
            row.update(
                {
                    "blast_hit": "Yes",
                    "top_accession": top["sacc"],
                    "top_protein_name": blast_product,
                    "top_gene_name": blast_gene,
                    "pident": top["pident"],
                    "qcovs": top["qcovs"],
                    "evalue": top["evalue"],
                    "bitscore": top["bitscore"],
                    "top_hit_confidence": confidence,
                    "suggested_gene": blast_gene if action == "accept_blast_supported_name" else target.get("Refined name", ""),
                    "suggested_product": blast_product if action in {"accept_blast_supported_name", "accept_blast_supported_product_only", "review_low_confidence_hit"} else target.get("Original product", ""),
                    "suggested_action": action,
                }
            )
            top_row = {**target, **top}
            top_row["blast_protein_name"] = blast_product
            top_row["blast_gene_name"] = blast_gene
            top_row["top_hit_confidence"] = confidence
            top_hit_rows.append(top_row)
        suggestion_rows.append(row)
        if row["suggested_action"] in {"accept_blast_supported_name", "accept_blast_supported_product_only", "review_low_confidence_hit"}:
            supported_rows.append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / f"{args.prefix}_swissprot_blast_all_hits.tsv", hit_fields, all_hit_rows)
    write_tsv(args.output_dir / f"{args.prefix}_swissprot_blast_top_hits.tsv", hit_fields, top_hit_rows)
    write_tsv(args.output_dir / f"{args.prefix}_swissprot_annotation_suggestions.tsv", suggestion_fields, suggestion_rows)
    write_tsv(args.output_dir / f"{args.prefix}_swissprot_supported_updates.tsv", suggestion_fields, supported_rows)

    confidence_counts = Counter(row["top_hit_confidence"] for row in top_hit_rows)
    action_counts = Counter(row["suggested_action"] for row in suggestion_rows)
    summary = [
        f"# {args.prefix} Swiss-Prot BLAST confirmation summary",
        "",
        f"Database: {args.database_note}",
        "",
        f"- Query proteins: {len(targets)}",
        f"- BLAST hit lines: {len(hits)}",
        f"- Queries with at least one Swiss-Prot hit at e <= 1e-5: {len(hits_by_query)}",
        f"- Queries without Swiss-Prot hit: {len(targets) - len(hits_by_query)}",
        "",
        "## Top-hit confidence counts",
        "",
    ]
    for key, count in confidence_counts.most_common():
        summary.append(f"- {key}: {count}")
    summary.extend(["", "## Suggested action counts", ""])
    for key, count in action_counts.most_common():
        summary.append(f"- {key}: {count}")
    (args.output_dir / f"{args.prefix}_swissprot_blast_summary.md").write_text("\n".join(summary) + "\n")

    print("hits", len(hits))
    print("queries_with_hits", len(hits_by_query))
    print("suggestions", len(suggestion_rows))
    print("accepted_or_review", len(supported_rows))
    print("top confidence", dict(confidence_counts))
    print("actions", dict(action_counts))


if __name__ == "__main__":
    main()
