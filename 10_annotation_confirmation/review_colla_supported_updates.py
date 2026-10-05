#!/usr/bin/env python3
# Title          : review_colla_supported_updates.py
# Description    : Review accepted gene-name updates for symbol stability
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/24
# Usage          : python3 review_colla_supported_updates.py --suggestions <tsv> --out-dir <dir>

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


PATHWAY_GENE_PATTERNS = [
    r"^pca[A-Z0-9a-z]*$",
    r"^cat[A-Z0-9a-z]*$",
    r"^mhp[A-Z0-9a-z]*$",
    r"^paa[A-Z0-9a-z]*$",
    r"^fad[A-Z0-9a-z]*$",
    r"^etf[A-Z0-9a-z]*$",
    r"^phh[A-Z0-9a-z]*$",
    r"^chn[A-Z0-9a-z]*$",
    r"^acmsd$",
    r"^qhn[A-Z0-9a-z]*$",
    r"^yoaI$",
    r"^hglS$",
    r"^ipd[A-Z0-9a-z]*$",
    r"^echA\d+[A-Za-z]*$",
    r"^ltp\d*[A-Za-z]*$",
    r"^chsH\d*[A-Za-z]*$",
    r"^gbcA$",
    r"^folA$",
    r"^rocR$",
    r"^pdeB$",
    r"^mucR$",
    r"^nrfD$",
    r"^tfpB$",
    r"^moxR\d*[A-Za-z]*$",
    r"^comD$",
    r"^ndhG$",
    r"^pcl\d+[A-Za-z]*$",
    r"^htdZ$",
    r"^tsdA$",
    r"^trhO$",
]

GENERIC_PRODUCT_TERMS = [
    "upf",
    "outer membrane",
    "binding protein",
    "chaperone-like",
    "family protein",
    "domain-containing",
    "protein ",
]

MOBILE_TERMS = [
    "transposase",
    "insertion sequence",
    "ins",
    "tnp",
]

SYSTEMATIC_GENE_PATTERNS = [
    r"^[A-Z]{2,}\d{3,}[A-Za-z]*$",
    r"^[A-Z][A-Za-z0-9]+_\d+[A-Za-z]*$",
    r"^[A-Z]{2,}[A-Za-z0-9]*_\d+[A-Za-z]*$",
    r"^[A-Z][a-z]+\d+[A-Za-z]*$",
    r"^[A-Z]{2,}\d+[A-Za-z]*_[A-Z0-9]+$",
]

EUK_STYLE_GENE_PATTERNS = [
    r"^[A-Z][A-Z0-9]+[a-z][A-Za-z0-9]*$",
    r"^[A-Z][a-z]+[A-Z0-9][A-Za-z0-9]*$",
]

REVIEW_FIELDS = [
    "Genome",
    "Cluster ID",
    "Locus tag",
    "Original gene",
    "Refined name",
    "Original product",
    "Enzyme class",
    "Ambiguity flags",
    "Confidence",
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
    "manual_review_action",
    "curated_gene_name",
    "curated_product_name",
    "manual_review_reason",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in rows])


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def matches_any(value: str, patterns: list[str]) -> bool:
    return any(re.match(pattern, value or "") for pattern in patterns)


def contains_any(value: str, terms: list[str]) -> bool:
    value_l = (value or "").lower()
    return any(term in value_l for term in terms)


def is_pathway_gene(gene: str) -> bool:
    return matches_any(gene, PATHWAY_GENE_PATTERNS)


def is_systematic_gene(gene: str) -> bool:
    if matches_any(gene, SYSTEMATIC_GENE_PATTERNS):
        return True
    if gene.startswith(("PA", "VF_", "Pfl", "ABS", "TCP", "PD_", "SE_", "SH", "HI_", "MW", "SAR", "HP_")):
        return True
    return False


def is_y_gene(gene: str) -> bool:
    return bool(re.match(r"^y[a-z]{2,}[A-Za-z0-9]*$", gene or ""))


def is_euk_style_gene(gene: str) -> bool:
    return matches_any(gene, EUK_STYLE_GENE_PATTERNS)


def is_mobile(row: dict[str, str]) -> bool:
    return contains_any(row.get("suggested_gene", ""), MOBILE_TERMS) or contains_any(row.get("suggested_product", ""), MOBILE_TERMS)


def clean_product_name(product: str) -> str:
    product = (product or "").strip()
    product = re.sub(r"^(Probable|Putative)\s+", "", product, flags=re.I)
    return product


def review_row(row: dict[str, str]) -> tuple[str, str, str, str]:
    gene = (row.get("suggested_gene") or "").strip()
    product = clean_product_name(row.get("suggested_product", ""))
    top_conf = (row.get("top_hit_confidence") or "").strip().lower()
    enzyme_class = (row.get("Enzyme class") or "").strip()
    pident = as_float(row.get("pident", "0"))
    qcovs = as_float(row.get("qcovs", "0"))
    bitscore = as_float(row.get("bitscore", "0"))

    if is_mobile(row):
        return (
            "retain_existing_mobile_or_accessory_annotation",
            "",
            product or row.get("Refined name", ""),
            "Mobile-element-like Swiss-Prot hit; keep accessory/mobile annotation and do not promote a gene symbol.",
        )

    if contains_any(product, ["upf", "uncharacterized"]) or is_y_gene(gene):
        return (
            "downgrade_to_product_only",
            "",
            product or row.get("Refined name", ""),
            "Hit supports a product description but the proposed gene symbol is generic, provisional, or not pathway-informative.",
        )

    if is_systematic_gene(gene) or is_euk_style_gene(gene):
        return (
            "downgrade_to_product_only",
            "",
            product or row.get("Refined name", ""),
            "Top hit supports a function, but the proposed gene name is a locus-tag-like or non-standard ortholog label rather than a stable cross-genome symbol.",
        )

    if contains_any(product, GENERIC_PRODUCT_TERMS):
        return (
            "downgrade_to_product_only",
            "",
            product or row.get("Refined name", ""),
            "Swiss-Prot hit is informative at the product level, but the description remains too generic for a confident gene-symbol rename.",
        )

    if is_pathway_gene(gene) and qcovs >= 80 and bitscore >= 120 and (top_conf == "high" or pident >= 40):
        return (
            "accept_gene_symbol_and_product",
            gene,
            product or row.get("Refined name", ""),
            "High-support hit to a pathway-relevant or biochemically interpretable gene family with a stable short gene symbol.",
        )

    if enzyme_class in {
        "CoA ligase/synthetase",
        "Dehydrogenase/oxidoreductase",
        "Hydratase/dehydratase/isomerase",
        "Oxygenase/dioxygenase",
        "Hydrolase/decarboxylase/lyase",
        "Electron transfer/redox",
        "Carboxylase",
    } and qcovs >= 85 and bitscore >= 150 and pident >= 40 and top_conf in {"high", "medium"}:
        return (
            "downgrade_to_product_only",
            "",
            product or row.get("Refined name", ""),
            "Functional support is strong enough to refine the product name, but the proposed gene symbol is not sufficiently standardized for cross-genome use.",
        )

    return (
        "needs_case_by_case_followup",
        "",
        product or row.get("Refined name", ""),
        "Suggestion remains plausible but still needs additional context, orthology, or domain support before a gene-level rename.",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Manual-review layer for Colla Swiss-Prot-supported gene-name suggestions.")
    parser.add_argument("--suggestions", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="Colla_batch")
    args = parser.parse_args()

    rows = read_tsv(args.suggestions)
    accepted = [row for row in rows if row.get("suggested_action") == "accept_blast_supported_name"]

    reviewed_rows: list[dict[str, str]] = []
    counts = Counter()
    for row in accepted:
        action, curated_gene, curated_product, reason = review_row(row)
        counts[action] += 1
        out = dict(row)
        out["manual_review_action"] = action
        out["curated_gene_name"] = curated_gene
        out["curated_product_name"] = curated_product
        out["manual_review_reason"] = reason
        reviewed_rows.append(out)

    reviewed_rows.sort(key=lambda r: (r["manual_review_action"], r["Genome"], r["Cluster ID"], r["Locus tag"]))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    reviewed_tsv = args.out_dir / f"{args.prefix}_manual_reviewed_gene_renames.tsv"
    write_tsv(reviewed_tsv, REVIEW_FIELDS, reviewed_rows)

    summary = args.out_dir / f"{args.prefix}_manual_reviewed_gene_renames_summary.md"
    with summary.open("w") as handle:
        handle.write(f"# {args.prefix} manual review of Swiss-Prot gene-name suggestions\n\n")
        handle.write(f"- Reviewed raw name-suggestion rows: {len(accepted)}\n")
        for key, value in sorted(counts.items()):
            handle.write(f"- {key}: {value}\n")
        handle.write("\n## Review logic\n\n")
        handle.write("- Mobile-element-like hits were not promoted to gene symbols.\n")
        handle.write("- Locus-tag-like, systematic, provisional `y*`, and eukaryotic-style ortholog labels were downgraded to product-only support.\n")
        handle.write("- Stable pathway-relevant gene symbols with strong Swiss-Prot support were retained as gene-level renames.\n")
        handle.write("- Ambiguous residual cases were left for additional follow-up rather than overcalled.\n")

    print(f"Reviewed {len(accepted)} raw gene-name suggestions")
    for key, value in sorted(counts.items()):
        print(f"{key}\t{value}")
    print(f"Wrote {reviewed_tsv}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
