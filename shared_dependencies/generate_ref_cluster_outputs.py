#!/usr/bin/env python3
# Title          : generate_ref_cluster_outputs.py
# Description    : Shared enzyme_class() functional vocabulary and SVG gene-map helpers (imported module)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/23
# Usage          : python3 generate_ref_cluster_outputs.py (imported by other scripts)

from __future__ import annotations

import csv
import html
import math
import re
import textwrap
from collections import Counter, defaultdict
from pathlib import Path


# CP010516.1/CP010517.1 (two replicons of Cupriavidus gilardii CR3) and CP063454.1/CP063455.1
# are kept as separate entries, matching the 47-genome comparator set.
SKIP_GENOMES: set[str] = set()

SUMMARY_FIELDS = [
    "Genome",
    "Cluster ID",
    "Gene IDs",
    "Annotation (refined)",
    "Relevance (Yes/No)",
    "Functional Role",
    "Evidence",
    "Category",
]

REFINEMENT_FIELDS = [
    "Genome",
    "Cluster ID",
    "Locus tag",
    "Original gene",
    "Refined name",
    "Original product",
    "Refined function",
    "Enzyme class",
    "EC",
    "COG",
    "Ambiguity flags",
    "Relevant rule hit",
    "Inference source",
    "Confidence",
    "Contig",
    "Start",
    "End",
    "Strand",
]

GENERIC_PAT = re.compile(
    r"hypothetical|DUF|domain-containing|unknown|uncharacterized|^HP$|\bHP\b|family protein",
    re.I,
)

COLOR = {
    "Oxygenase/dioxygenase": "#d73027",
    "CoA ligase/synthetase": "#1f78b4",
    "Carboxylase": "#6a3d9a",
    "Dehydrogenase/oxidoreductase": "#33a02c",
    "Hydratase/dehydratase/isomerase": "#ff7f00",
    "Hydrolase/decarboxylase/lyase": "#b15928",
    "Thiolase/transferase": "#a6cee3",
    "Transport/efflux": "#cab2d6",
    "Regulator": "#fdbf6f",
    "Electron transfer/redox": "#fb9a99",
    "Mobile/stress/accessory": "#b2df8a",
    "Hypothetical/accessory": "#d9d9d9",
    "Accessory/other": "#cccccc",
}

CATEGORY_COLOR = {
    "High Priority": "#b2182b",
    "Medium Priority": "#ef8a62",
    "Low Priority": "#878787",
}

REFERENCES = {
    "R1": (
        "Wasi et al. 2023, microbial degradation of naphthenic acids: metabolic and genomic insights",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10710301/",
    ),
    "R3": (
        "Harwood & Parales 1996, beta-ketoadipate pathway review",
        "https://pubmed.ncbi.nlm.nih.gov/8905091/",
    ),
    "R4": (
        "Teufel et al. 2010, phenylacetate catabolic pathway",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC2922514/",
    ),
    "R5": (
        "Diao et al. 2013, cyclohexanecarboxyl-CoA dehydrogenases",
        "https://pubmed.ncbi.nlm.nih.gov/23667239/",
    ),
    "R7": (
        "Diaz et al. 1998, hca cluster for 3-phenylpropionic acid catabolism",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC107259/",
    ),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in rows])


def parse_rule_counts(value: str) -> Counter:
    counts = Counter()
    for item in (value or "").split(";"):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            key, val = item.split(":", 1)
            try:
                counts[key.strip()] += int(val)
            except ValueError:
                counts[key.strip()] += 1
        else:
            counts[item] += 1
    return counts


def cluster_number(source_cluster_id: str) -> int:
    match = re.search(r"_cluster_(\d+)$", source_cluster_id)
    if not match:
        raise ValueError(f"Cannot parse cluster number from {source_cluster_id!r}")
    return int(match.group(1))


def output_cluster_id(genome: str, source_cluster_id: str) -> str:
    return f"{genome.replace('.fa', '')}-C{cluster_number(source_cluster_id):02d}"


def source_cluster_id(genome: str, output_id: str) -> str:
    match = re.search(r"-C(\d+)$", output_id)
    if not match:
        raise ValueError(f"Cannot parse output cluster id {output_id!r}")
    return f"{genome}_cluster_{int(match.group(1))}"


def locate_annotation_tsv(folder: Path) -> Path:
    candidates = []
    for path in sorted(folder.glob("*.tsv")):
        name = path.name
        if name.startswith("NA_cluster_") or name.endswith("_cluster_assessment.tsv"):
            continue
        if name.endswith("_gene_maps.tsv") or name.endswith("_gene_refinement_layer.tsv"):
            continue
        candidates.append(path)
    prokka = [p for p in candidates if p.name.startswith("PROKKA_")]
    return (prokka or candidates)[0]


def locate_faa(folder: Path) -> Path:
    faa = sorted(folder.glob("*.faa"))
    if not faa:
        raise FileNotFoundError(f"No .faa file in {folder}")
    prokka = [p for p in faa if p.name.startswith("PROKKA_")]
    return (prokka or faa)[0]


def locate_coord_table(folder: Path) -> Path:
    ods = sorted(folder.glob("PROKKA_*.ods"))
    if not ods:
        ods = sorted(folder.glob("*.ods"))
    if not ods:
        raise FileNotFoundError(f"No Prokka coordinate .ods/TSV file in {folder}")
    # The Prokka .ods files in this dataset are tab-separated text.
    non_curated = [p for p in ods if not p.name.endswith("_2.ods")]
    return (non_curated or ods)[0]


def parse_gff_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = val
    return attrs


def read_coords(folder: Path) -> list[dict[str, str]]:
    coord_path = locate_coord_table(folder)
    try:
        with coord_path.open("rb") as handle:
            prefix = handle.read(4)
        if not prefix.startswith(b"PK"):
            return read_tsv(coord_path)
    except UnicodeDecodeError:
        pass

    gff_files = sorted(folder.glob("PROKKA_*.gff")) or sorted(folder.glob("*.gff"))
    if not gff_files:
        raise FileNotFoundError(f"No text coordinate table or GFF found in {folder}")
    rows: list[dict[str, str]] = []
    with gff_files[0].open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_gff_attributes(parts[8])
            locus = attrs.get("ID") or attrs.get("locus_tag")
            if not locus:
                continue
            rows.append(
                {
                    "Seqid": parts[0],
                    "Start": parts[3],
                    "End": parts[4],
                    "Strand": parts[6],
                    "ID": locus,
                    "inference": attrs.get("inference", ""),
                }
            )
    return rows


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current:
        seqs[current] = "".join(chunks)
    return seqs


def tokens(gene: str) -> list[str]:
    return [token for token in re.split(r"[^a-z0-9]+", (gene or "").lower()) if token]


def starts_any(gene_tokens: list[str], *prefixes: str) -> bool:
    return any(any(token.startswith(prefix) for prefix in prefixes) for token in gene_tokens)


def has_text(text: str, *patterns: str) -> bool:
    return any(re.search(pattern, text or "", re.I) for pattern in patterns)


def enzyme_class(gene: str, product: str) -> str:
    product = product or ""
    gene_tokens = tokens(gene)
    if has_text(
        product,
        r"transcriptional regulator",
        r"\bregulator\b",
        r"\brepressor\b",
        r"\bactivator\b",
        r"sensor histidine kinase",
        r"response regulator",
    ) or starts_any(
        gene_tokens,
        "kdgr",
        "paax",
        "kstr",
        "ytr",
        "xynr",
        "hcar",
        "kipr",
        "ydjf",
        "tipa",
        "acor",
        "lutr",
        "walr",
        "srr",
        "yiaj",
        "bagr",
        "teni",
        "beti",
        "iclr",
        "gntr",
        "lysr",
        "tetr",
        "marr",
        "xre",
    ):
        return "Regulator"
    if has_text(
        product,
        r"transport",
        r"transporter",
        r"porin",
        r"symporter",
        r"permease",
        r"efflux",
        r"solute-binding",
        r"substrate-binding",
        r"\bMFS\b",
    ) or starts_any(gene_tokens, "actp", "satp", "pcak", "benk", "bene", "nicp", "nict", "tctc", "aaea", "aaeb", "aaex", "oprm", "oprd", "exut", "ydco"):
        return "Transport/efflux"
    if has_text(product, r"electron transfer flavoprotein", r"quinone oxidoreductase", r"ferredoxin", r"iron-sulfur", r"flavodoxin", r"rubredoxin") or starts_any(gene_tokens, "etfa", "etfb", "ndh", "nuo", "nqo", "fadf", "yfhp", "alke", "alkt"):
        return "Electron transfer/redox"
    if has_text(product, r"coa ligase", r"coa synthetase", r"phenylacetate-coenzyme a ligase", r"fatty-acid--coa ligase", r"acyl-coa synthetase") or starts_any(gene_tokens, "fadd", "lcfb", "bcla", "paak", "fcs"):
        return "CoA ligase/synthetase"
    if has_text(product, r"oxygenase", r"dioxygenase", r"hydroxylase", r"monooxygenase") or starts_any(gene_tokens, "hpcb", "nagx", "hpab", "alkb", "vana", "bena", "benb", "benc", "pcag", "pcah", "cata", "hmga", "ksha", "kshe", "naga", "ant"):
        return "Oxygenase/dioxygenase"
    if has_text(product, r"carboxylase", r"carboxyltransferase") or starts_any(gene_tokens, "mccb", "acca", "atuc", "atuf"):
        return "Carboxylase"
    if has_text(product, r"dehydrogenase", r"oxidoreductase", r"\breductase\b", r"\bSDR\b") or starts_any(gene_tokens, "gbs", "betb", "amn", "aldh", "carc", "acda", "mmgc", "bcd", "fade", "dmdc", "bac", "hbd", "adh", "yajo", "pao", "vdh", "hca"):
        return "Dehydrogenase/oxidoreductase"
    if has_text(product, r"hydratase", r"dehydratase", r"isomerase", r"crotonase") or starts_any(gene_tokens, "fadb", "echa", "crt", "paag", "paaf", "fadn", "maia", "nagl", "menb", "phaj"):
        return "Hydratase/dehydratase/isomerase"
    if has_text(product, r"thiolase", r"coa-transferase", r"coenzyme a transferase", r"acyltransferase") or starts_any(gene_tokens, "fada", "paaj", "catj", "gcta", "yfde", "uctc", "pcai", "pcaj", "thla"):
        return "Thiolase/transferase"
    if has_text(product, r"hydrolase", r"thioesterase", r"decarboxylase", r"lyase", r"aldolase", r"lactonase") or starts_any(gene_tokens, "paai", "paaz", "menc", "pcac", "catd", "pcad", "faha", "yhaa", "menh"):
        return "Hydrolase/decarboxylase/lyase"
    if has_text(product, r"transposase", r"integrase", r"stress protein", r"SOS response"):
        return "Mobile/stress/accessory"
    if GENERIC_PAT.search(((gene or "") + " " + product).strip()):
        return "Hypothetical/accessory"
    return "Accessory/other"


def accession_from_text(value: str) -> str:
    match = re.search(r"UniProtKB:([A-Z0-9]+)", value or "")
    return match.group(1) if match else ""


def product_label(product: str) -> str:
    product = (product or "").strip()
    if not product:
        return ""
    replacements = [
        (r"Long-chain-fatty-acid--CoA ligase FadD13", "fadD13"),
        (r"Long-chain-fatty-acid--CoA ligase FadD15", "fadD15"),
        (r"Long-chain-fatty-acid--CoA ligase", "fadD/lcfB-like"),
        (r"Phenylacetate-coenzyme A ligase", "paaK"),
        (r"4-hydroxyphenylacetate 3-monooxygenase", "hpaB-like"),
        (r"3,4-dihydroxyphenylacetate 2,3-dioxygenase", "hpcB-like"),
        (r"Acetaldehyde dehydrogenase", "aldB-like"),
        (r"4-hydroxy-2-oxovalerate aldolase", "mhpE-like aldolase"),
        (r"Acyl-CoA dehydrogenase FadE34", "fadE34"),
        (r"Acyl-CoA dehydrogenase FadE26", "fadE26"),
        (r"Acyl-CoA dehydrogenase FadE17", "fadE17"),
        (r"Crotonyl-CoA hydratase", "crt-like"),
        (r"4-hydroxybenzoyl-CoA thioesterase", "4hbt-like"),
        (r"Methylmalonyl-CoA carboxyltransferase", "mccB-like"),
        (r"NADH:quinone reductase", "NADH-quinone reductase"),
        (r"Cyclic di-GMP phosphodiesterase", "cyclic-di-GMP phosphodiesterase"),
        (r"Solute-binding protein", "substrate-binding protein"),
        (r"Putative universal stress protein", "universal stress protein"),
        (r"hypothetical protein", "hypothetical/accessory protein"),
    ]
    for pattern, label in replacements:
        if re.search(pattern, product, re.I):
            return label
    return product.split("/")[0][:55]


def load_review_maps(review_dir: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    review: dict[tuple[str, str, str], dict[str, str]] = {}

    def put(row: dict[str, str], gene: str, product: str, confidence: str, source: str) -> None:
        genome = row.get("genome") or row.get("Genome")
        cluster = row.get("cluster_id") or row.get("Cluster ID")
        locus = row.get("locus_tag") or row.get("Locus tag")
        if not genome or not cluster or not locus:
            return
        key = (genome, cluster, locus)
        current = review.get(key)
        priority = {"high": 3, "medium": 2, "low": 1, "local_only": 1}.get((confidence or "").lower(), 0)
        old_priority = {"high": 3, "medium": 2, "low": 1, "local_only": 1}.get((current or {}).get("confidence", "").lower(), -1)
        if current and old_priority > priority:
            return
        review[key] = {
            "gene": gene or "",
            "product": product or "",
            "confidence": confidence or "",
            "source": source,
            "accession": row.get("inferred_uniprot_accession", ""),
            "entry": row.get("uniprot_entry_name", ""),
            "reviewed": row.get("uniprot_reviewed", ""),
            "status": row.get("verification_status", row.get("review_decision", "")),
        }

    for name in [
        "UniProt_manual_review_assignments.tsv",
        "UniProt_manual_review_decisions.tsv",
    ]:
        for row in read_tsv(review_dir / name):
            gene = row.get("reviewed_gene") or row.get("proposed_gene") or row.get("corrected_gene") or ""
            product = row.get("reviewed_product") or row.get("proposed_product") or row.get("corrected_product") or ""
            confidence = row.get("review_confidence") or row.get("proposal_confidence") or row.get("correction_confidence") or ""
            if gene or product:
                put(row, gene, product, confidence, name)

    for name in [
        "UniProt_annotation_corrections.tsv",
        "UniProt_possible_misnamed_gene_review.tsv",
    ]:
        for row in read_tsv(review_dir / name):
            gene = row.get("corrected_gene") or row.get("proposed_gene") or ""
            product = row.get("corrected_product") or row.get("proposed_product") or ""
            confidence = row.get("correction_confidence") or row.get("proposal_confidence") or ""
            if gene or product:
                put(row, gene, product, confidence, name)

    for name in [
        "UniProt_manual_curation_round2.tsv",
        "UniProt_manual_curation_round1.tsv",
        "UniProt_unnamed_gene_review.tsv",
        "UniProt_cluster_gene_verification_all.tsv",
        "UniProt_low_priority_gene_verification.tsv",
    ]:
        for row in read_tsv(review_dir / name):
            gene = row.get("proposed_gene") or row.get("corrected_gene") or row.get("reviewed_gene") or ""
            product = row.get("proposed_product") or row.get("corrected_product") or row.get("reviewed_product") or row.get("uniprot_protein_name") or ""
            confidence = row.get("proposal_confidence") or row.get("correction_confidence") or row.get("review_confidence") or ""
            if gene or product:
                put(row, gene, product, confidence, name)

    return review


def refined_from_review(
    genome: str,
    cluster_id: str,
    locus: str,
    original_gene: str,
    product: str,
    review: dict[tuple[str, str, str], dict[str, str]],
) -> tuple[str, str, str, str]:
    row = review.get((genome, cluster_id, locus))
    if row:
        name = row["gene"] or original_gene or product_label(row["product"] or product) or locus
        product_out = row["product"] or product
        source = f"{row['source']}"
        if row.get("accession"):
            source += f"; UniProtKB:{row['accession']}"
        if row.get("entry"):
            source += f"; entry:{row['entry']}"
        if row.get("status"):
            source += f"; status:{row['status']}"
        return name, product_out, row.get("confidence", ""), source
    if original_gene and original_gene.upper() != "HP":
        return original_gene, product, "high", "local Prokka gene/product"
    label = product_label(product)
    if label:
        confidence = "low" if label == "hypothetical/accessory protein" else "medium"
        return label, product, confidence, "product-derived local annotation"
    return locus, product, "low", "no local product annotation"


def ambiguity_flags(original_gene: str, product: str, refined_name: str) -> str:
    flags: list[str] = []
    if not original_gene or original_gene.upper() == "HP":
        flags.append("no_gene")
    if GENERIC_PAT.search(((original_gene or "") + " " + (product or "")).strip()):
        flags.append("generic_or_hypothetical")
    if refined_name == "hypothetical/accessory protein":
        flags.append("unresolved_function")
    return ";".join(dict.fromkeys(flags)) if flags else "none"


def refined_function(cls: str, product: str, candidate: dict[str, str] | None) -> str:
    if candidate and candidate.get("potential function in the cluster"):
        return candidate["potential function in the cluster"]
    if cls == "Hypothetical/accessory":
        return "Ambiguous accessory gene inside the candidate span; no confident function from local annotation."
    mapping = {
        "CoA ligase/synthetase": "Activates carboxylic-acid substrates as CoA thioesters for downstream beta-oxidation or ring processing.",
        "Oxygenase/dioxygenase": "Introduces oxygen into aromatic or hydroxylated substrates; likely ring activation/cleavage entry step.",
        "Dehydrogenase/oxidoreductase": "Oxidizes alcohol/aldehyde/acyl-CoA intermediates and supports beta-oxidation-like flux.",
        "Hydratase/dehydratase/isomerase": "Processes enoyl-CoA or ring-opened intermediates during beta-oxidation-like metabolism.",
        "Thiolase/transferase": "CoA transfer, thiolysis, or CoA recycling in lower aromatic/fatty-acid pathways.",
        "Hydrolase/decarboxylase/lyase": "Hydrolysis, ring-opened intermediate processing, or thioester control.",
        "Transport/efflux": "Imports carboxylates/aromatic acids or exports toxic intermediates.",
        "Regulator": "Controls inducible carboxylate/aromatic degradation genes.",
        "Electron transfer/redox": "Supplies redox/electron-transfer capacity for anaerobic or oxygenase-linked metabolism.",
        "Carboxylase": "Handles branched acyl-CoA or propionyl/methylmalonyl-CoA intermediates.",
        "Mobile/stress/accessory": "Mobile, stress, or accessory function near the candidate cluster.",
    }
    return mapping.get(cls, "Accessory gene; indirect or unknown pathway role.")


def relevance_and_category(summary: dict[str, str]) -> tuple[str, str]:
    counts = parse_rule_counts(summary.get("matched_rule_counts", ""))
    direct_keys = {"paa_operon", "beta_oxidation", "coa_activation", "aromatic_or_cyclic"}
    direct_score = sum(counts[key] for key in direct_keys)
    only_support = direct_score == 0 and any(counts.values())
    if only_support:
        return "No", "Low Priority"
    confidence = summary.get("confidence", "").lower()
    if confidence == "high":
        return "Yes", "High Priority"
    if confidence == "medium":
        return "Yes", "Medium Priority"
    if direct_score >= 6 or (counts.get("paa_operon", 0) >= 2 and counts.get("coa_activation", 0) >= 1):
        return "Yes", "Medium Priority"
    return "Yes", "Low Priority"


def functional_role(summary: dict[str, str], category: str) -> str:
    role = summary.get("overall_role", "Candidate NA-like degradation cluster")
    counts = parse_rule_counts(summary.get("matched_rule_counts", ""))
    pieces = []
    if counts.get("paa_operon"):
        pieces.append("phenylacetate/CoA-dependent aromatic-acid processing")
    if counts.get("aromatic_or_cyclic"):
        pieces.append("aromatic or cyclic scaffold modification")
    if counts.get("coa_activation"):
        pieces.append("carboxylate CoA activation")
    if counts.get("beta_oxidation"):
        pieces.append("beta-oxidation-like acyl-CoA processing")
    if counts.get("transport"):
        pieces.append("substrate transport")
    if counts.get("electron_transfer"):
        pieces.append("redox/electron-transfer support")
    if counts.get("regulation"):
        pieces.append("local transcriptional regulation")
    if not pieces:
        return f"{role}; low-specificity support only."
    qualifier = "Clear" if category == "High Priority" else "Partial" if category == "Medium Priority" else "Weak/indirect"
    return f"{qualifier} {role}: " + "; ".join(pieces) + "."


def evidence_for(summary: dict[str, str], relevance: str) -> str:
    if relevance == "No":
        return "No direct degradation evidence; supporting respiratory/redox/transport genes only."
    role = summary.get("overall_role", "").lower()
    refs = ["R1"]
    if "phenylacetate" in role:
        refs.append("R4")
    if "benzoate" in role or "hydroxybenzoate" in role or "cinnamate" in role:
        refs.append("R3")
        refs.append("R7")
    if "alicyclic" in role or "cycloalkyl" in role or "steroid" in role:
        refs.append("R5")
    return "; ".join(dict.fromkeys(refs))


def cluster_records(
    genome: str,
    summary: dict[str, str],
    annotation: dict[str, dict[str, str]],
    coords: list[dict[str, str]],
    candidates: dict[str, dict[str, str]],
    review: dict[tuple[str, str, str], dict[str, str]],
) -> list[dict[str, str]]:
    source_id = summary["cluster_id"]
    contig = summary["contig"]
    start = int(summary["start"])
    end = int(summary["end"])
    records: list[dict[str, str]] = []
    for coord in coords:
        if coord.get("Seqid") != contig or not coord.get("Start", "").isdigit():
            continue
        coord_start = int(coord["Start"])
        coord_end = int(coord["End"])
        if coord_start < start or coord_end > end:
            continue
        locus = coord["ID"]
        ann = annotation.get(locus, {})
        original_gene = ann.get("gene", "")
        product = ann.get("product", "")
        refined_name, refined_product, review_conf, review_source = refined_from_review(
            genome, source_id, locus, original_gene, product, review
        )
        cls = enzyme_class(refined_name if refined_name != locus else original_gene, refined_product or product)
        flags = ambiguity_flags(original_gene, product, refined_name)
        candidate = candidates.get(locus)
        inference_parts = [review_source, "NA_cluster_candidates rule hit" if candidate else "intervening CDS in expanded span"]
        accession = accession_from_text(coord.get("inference", ""))
        if accession and "UniProtKB:" not in "; ".join(inference_parts):
            inference_parts.append(f"UniProtKB:{accession}")
        confidence = review_conf or "high"
        if cls == "Hypothetical/accessory":
            confidence = "low"
        elif flags != "none" and confidence == "high":
            confidence = "medium"
        records.append(
            {
                "Genome": genome,
                "Cluster ID": output_cluster_id(genome, source_id),
                "Locus tag": locus,
                "Original gene": original_gene,
                "Refined name": refined_name,
                "Original product": product,
                "Refined function": refined_function(cls, refined_product or product, candidate),
                "Enzyme class": cls,
                "EC": ann.get("EC_number", ""),
                "COG": ann.get("COG", ""),
                "Ambiguity flags": flags,
                "Relevant rule hit": "Yes" if candidate else "No",
                "Inference source": "; ".join(part for part in inference_parts if part),
                "Confidence": confidence,
                "Contig": contig,
                "Start": str(coord_start),
                "End": str(coord_end),
                "Strand": coord.get("Strand", ""),
            }
        )
    return sorted(records, key=lambda row: (int(row["Start"]), int(row["End"])))


def cluster_annotation(summary: dict[str, str], records: list[dict[str, str]]) -> str:
    counts = Counter(row["Enzyme class"] for row in records)
    order = [
        "Oxygenase/dioxygenase",
        "CoA ligase/synthetase",
        "Carboxylase",
        "Dehydrogenase/oxidoreductase",
        "Hydratase/dehydratase/isomerase",
        "Hydrolase/decarboxylase/lyase",
        "Thiolase/transferase",
        "Transport/efflux",
        "Regulator",
        "Electron transfer/redox",
        "Mobile/stress/accessory",
        "Hypothetical/accessory",
        "Accessory/other",
    ]
    parts = [f"{key} x{counts[key]}" for key in order if counts.get(key)]
    key_genes = [
        f"{row['Refined name']}:{row['Enzyme class']}"
        for row in records
        if row["Relevant rule hit"] == "Yes"
    ]
    if len(key_genes) > 18:
        key_genes = key_genes[:18] + ["..."]
    flagged = [row["Refined name"] for row in records if row["Ambiguity flags"] != "none"]
    text = (
        f"Expanded from candidate span {summary['contig']}:{summary['start']}-{summary['end']}; "
        f"refined role classes: {', '.join(parts)}. "
        f"Key pathway genes: {'; '.join(key_genes)}"
    )
    if flagged:
        text += "; no-name/hypothetical/generic genes needing caution: "
        text += ", ".join(flagged[:18])
        if len(flagged) > 18:
            text += ", ..."
    return text


def gene_map(records: list[dict[str, str]]) -> str:
    return " -> ".join(
        f"{row['Refined name']}({row['Strand'] or '?'}; {row['Enzyme class'].replace('/', '-')})"
        for row in records
    )


def write_markdown(folder: Path, prefix: str, clusters: list[dict[str, str]], records_by_cluster: dict[str, list[dict[str, str]]]) -> None:
    def esc(value: str) -> str:
        return str(value).replace("|", "/").replace("\n", "<br>")

    path = folder / f"{prefix}_cluster_assessment.md"
    with path.open("w") as handle:
        handle.write(f"# {prefix} naphthenic-acid-like degradation cluster assessment\n\n")
        handle.write(
            "Input: Prokka outputs plus `NA_cluster_candidates.tsv` and `NA_cluster_summary.tsv`. "
            "Candidate spans were expanded to include all CDS between start/end coordinates so unnamed and hypothetical genes inside each operon/cluster are retained.\n\n"
        )
        handle.write(
            "Annotation-refinement note: this pass uses Prokka product/EC/COG calls, embedded Prokka UniProtKB inference accessions, "
            "the `NA_cluster_review` UniProt/manual curation tables where available, candidate-rule context, and pathway-context reasoning. "
            "True hypothetical proteins without supporting curation remain flagged as low-confidence rather than over-renamed.\n\n"
        )
        handle.write("## Strict Summary Table\n\n")
        handle.write("| " + " | ".join(SUMMARY_FIELDS) + " |\n")
        handle.write("| " + " | ".join(["---"] * len(SUMMARY_FIELDS)) + " |\n")
        for cluster in clusters:
            handle.write("| " + " | ".join(esc(cluster[field]) for field in SUMMARY_FIELDS) + " |\n")
        handle.write("\n## Gene Maps For Relevant Clusters\n\n")
        handle.write("Format: `refined-name(strand; enzyme class)`. The `->` separator follows genomic order.\n\n")
        for cluster in clusters:
            if cluster["Relevance (Yes/No)"] != "Yes":
                continue
            records = records_by_cluster[cluster["Cluster ID"]]
            handle.write(f"### {cluster['Cluster ID']} - {cluster['Category']}\n\n")
            handle.write(f"Functional role: {cluster['Functional Role']}\n\n")
            handle.write("```text\n")
            handle.write("\n".join(textwrap.wrap(gene_map(records), width=150, break_long_words=False, break_on_hyphens=False)))
            handle.write("\n```\n\n")
            handle.write("Per-gene role notes:\n\n")
            for record in records:
                flag = f"; FLAG {record['Ambiguity flags']}" if record["Ambiguity flags"] != "none" else ""
                rule = "; candidate-rule hit" if record["Relevant rule hit"] == "Yes" else "; intervening CDS"
                handle.write(
                    f"- {record['Locus tag']} ({record['Refined name']}): {record['Enzyme class']}; "
                    f"product={record['Original product'] or 'missing'}; function={record['Refined function']}; "
                    f"EC={record['EC'] or 'NA'}; COG={record['COG'] or 'NA'}{rule}{flag}\n"
                )
            handle.write("\n")
        handle.write("## Non-Relevant / Low-Confidence Notes\n\n")
        for cluster in clusters:
            if cluster["Relevance (Yes/No)"] == "No":
                handle.write(f"- {cluster['Cluster ID']} ({cluster['Category']}): {cluster['Functional Role']}\n")
        handle.write("\n## Literature Key\n\n")
        for key, (title, url) in REFERENCES.items():
            handle.write(f"- {key}: {title}. {url}\n")


def svg_text(x: float, y: float, text: str, size: int = 11, weight: str = "400", anchor: str = "start", fill: str = "#222") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{html.escape(str(text))}</text>'
    )


def wrap_label(text: str, width: int) -> list[str]:
    return textwrap.wrap(" ".join(str(text).split()), width=width, break_long_words=False, break_on_hyphens=False)


def arrow_points(x: float, y: float, width: float, height: float, strand: str) -> str:
    head = min(max(width * 0.33, 7), 18)
    if width < 18:
        head = width * 0.45
    if strand == "-":
        points = [(x + width, y), (x + head, y), (x, y + height / 2), (x + head, y + height), (x + width, y + height)]
    else:
        points = [(x, y), (x + width - head, y), (x + width, y + height / 2), (x + width - head, y + height), (x, y + height)]
    return " ".join(f"{px:.1f},{py:.1f}" for px, py in points)


def render_cluster_svg(cluster: dict[str, str], records: list[dict[str, str]], y: int, total_width: int = 1700) -> tuple[list[str], int]:
    left, label_width, right = 32, 370, 36
    plot_x = left + label_width
    plot_width = total_width - plot_x - right
    row_height, gene_height = 112, 22
    start = min(int(row["Start"]) for row in records)
    end = max(int(row["End"]) for row in records)
    scale = plot_width / max(end - start + 1, 1)
    elements = [
        f'<rect x="18" y="{y - 22}" width="{total_width - 36}" height="{row_height - 8}" rx="6" fill="#ffffff" stroke="#dddddd"/>',
        f'<rect x="18" y="{y - 22}" width="7" height="{row_height - 8}" fill="{CATEGORY_COLOR.get(cluster["Category"], "#666")}"/>',
        svg_text(left, y, f"{cluster['Cluster ID']} ({cluster['Category']})", 14, "700"),
        svg_text(left, y + 18, f"{cluster['Relevance (Yes/No)']} relevance", 11, "700", fill=CATEGORY_COLOR.get(cluster["Category"], "#666")),
    ]
    for idx, line in enumerate(wrap_label(cluster["Functional Role"], 52)[:4]):
        elements.append(svg_text(left, y + 35 + idx * 13, line, 10, fill="#444"))
    base = y + 43
    elements.append(f'<line x1="{plot_x}" y1="{base + gene_height / 2}" x2="{plot_x + plot_width}" y2="{base + gene_height / 2}" stroke="#d0d0d0"/>')
    for row in records:
        row_start = int(row["Start"])
        row_end = int(row["End"])
        gene_x = plot_x + (row_start - start) * scale
        gene_width = max((row_end - row_start + 1) * scale, 18)
        if gene_x + gene_width > plot_x + plot_width:
            gene_x = plot_x + plot_width - gene_width
        fill = COLOR.get(row["Enzyme class"], "#ccc")
        stroke = "#111" if row["Ambiguity flags"] != "none" else "#555"
        title = f"{row['Locus tag']} | {row['Refined name']} | {row['Enzyme class']} | {row['Original product']}"
        elements.append(
            f'<polygon points="{arrow_points(gene_x, base, gene_width, gene_height, row["Strand"])}" fill="{fill}" stroke="{stroke}" stroke-width="0.8">'
            f"<title>{html.escape(title)}</title></polygon>"
        )
        label = row["Refined name"]
        if len(label) > 22:
            label = label[:19] + "..."
        tx = gene_x + gene_width / 2
        ty = base + gene_height + 13
        if len(records) > 12:
            elements.append(
                f'<text x="{tx:.1f}" y="{ty:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="9" text-anchor="end" '
                f'transform="rotate(-45 {tx:.1f} {ty:.1f})" fill="#222">{html.escape(label)}</text>'
            )
        else:
            elements.append(svg_text(tx, ty, label, 9, anchor="middle"))
    elements.append(svg_text(plot_x, y + 10, f"{start:,} bp", 9, fill="#666"))
    elements.append(svg_text(plot_x + plot_width, y + 10, f"{end:,} bp", 9, anchor="end", fill="#666"))
    return elements, row_height


def write_svg(path: Path, genome: str, clusters: list[dict[str, str]], records_by_cluster: dict[str, list[dict[str, str]]]) -> None:
    total_width = 1700
    top = 100
    row_height = 112
    legend_rows = math.ceil(len(COLOR) / 4)
    height = top + len(clusters) * row_height + 56 + legend_rows * 24
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{height}" viewBox="0 0 {total_width} {height}">',
        '<rect width="100%" height="100%" fill="#f8f8f6"/>',
        svg_text(24, 34, f"{genome}: cluster/operon gene maps", 22, "700"),
        svg_text(24, 56, f"{len(clusters)} expanded candidate clusters from Prokka annotations and NA candidate tables", 12, fill="#555"),
        svg_text(24, 74, "Arrow direction shows CDS strand. Color shows inferred enzyme/function class. Dark outlines mark no-name/hypothetical/generic genes.", 11, fill="#555"),
    ]
    y = top
    for cluster in clusters:
        parts, consumed = render_cluster_svg(cluster, records_by_cluster[cluster["Cluster ID"]], y, total_width)
        elements.extend(parts)
        y += consumed
    legend_y = y + 20
    elements.append(svg_text(24, legend_y, "Legend", 13, "700"))
    x0, y0, cell_width = 90, legend_y - 12, 380
    for idx, (label, color) in enumerate(COLOR.items()):
        col = idx % 4
        row = idx // 4
        x = x0 + col * cell_width
        yy = y0 + row * 24
        elements.append(f'<rect x="{x}" y="{yy}" width="18" height="13" fill="{color}" stroke="#555" stroke-width="0.6"/>')
        elements.append(svg_text(x + 25, yy + 11, label, 10, fill="#333"))
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n")


def process_genome(folder: Path, review: dict[tuple[str, str, str], dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    genome = folder.name
    prefix = genome.replace(".fa", "")
    summary_rows = read_tsv(folder / "NA_cluster_summary.tsv")
    candidate_rows = read_tsv(folder / "NA_cluster_candidates.tsv")
    annotation_rows = read_tsv(locate_annotation_tsv(folder))
    coord_rows = read_coords(folder)
    seqs = read_fasta(locate_faa(folder))

    annotation = {row["locus_tag"]: row for row in annotation_rows if row.get("locus_tag")}
    candidates = {row["locus_tag"]: row for row in candidate_rows if row.get("locus_tag")}

    cluster_rows: list[dict[str, str]] = []
    refinement_rows: list[dict[str, str]] = []
    records_by_cluster: dict[str, list[dict[str, str]]] = {}
    reblast_rows: list[dict[str, str]] = []

    for summary in sorted(summary_rows, key=lambda row: cluster_number(row["cluster_id"])):
        relevance, category = relevance_and_category(summary)
        role = functional_role(summary, category)
        evidence = evidence_for(summary, relevance)
        records = cluster_records(genome, summary, annotation, coord_rows, candidates, review)
        cluster_id = output_cluster_id(genome, summary["cluster_id"])
        records_by_cluster[cluster_id] = records
        refinement_rows.extend(records)
        cluster = {
            "Genome": genome,
            "Cluster ID": cluster_id,
            "Gene IDs": "; ".join(row["Locus tag"] for row in records),
            "Annotation (refined)": cluster_annotation(summary, records),
            "Relevance (Yes/No)": relevance,
            "Functional Role": role,
            "Evidence": evidence,
            "Category": category,
        }
        cluster_rows.append(cluster)

    write_tsv(folder / f"{prefix}_cluster_assessment.tsv", SUMMARY_FIELDS, cluster_rows)
    write_tsv(folder / f"{prefix}_gene_refinement_layer.tsv", REFINEMENT_FIELDS, refinement_rows)
    map_fields = ["Genome", "Cluster ID", "Category", "Functional Role", "Gene Map"]
    map_rows = [
        {
            "Genome": genome,
            "Cluster ID": cluster["Cluster ID"],
            "Category": cluster["Category"],
            "Functional Role": cluster["Functional Role"],
            "Gene Map": gene_map(records_by_cluster[cluster["Cluster ID"]]),
        }
        for cluster in cluster_rows
        if cluster["Relevance (Yes/No)"] == "Yes"
    ]
    write_tsv(folder / f"{prefix}_gene_maps.tsv", map_fields, map_rows)
    write_markdown(folder, prefix, cluster_rows, records_by_cluster)

    svg_dir = folder / "cluster_operon_svgs"
    svg_dir.mkdir(exist_ok=True)
    write_svg(svg_dir / f"{prefix}_clusters_operons.svg", genome, cluster_rows, records_by_cluster)
    (svg_dir / "index.html").write_text(
        f'<!doctype html><meta charset="utf-8"><title>{html.escape(prefix)} cluster SVG</title>'
        f'<h1>{html.escape(prefix)} cluster/operon SVG</h1><ul><li><a href="{html.escape(prefix)}_clusters_operons.svg">'
        f'{html.escape(prefix)}_clusters_operons.svg</a></li></ul>'
    )

    with (folder / f"{prefix}_reblast_targets.faa").open("w") as fasta:
        for row in refinement_rows:
            if row["Ambiguity flags"] == "none":
                continue
            reblast_rows.append(
                {
                    "Genome": row["Genome"],
                    "Cluster ID": row["Cluster ID"],
                    "Locus tag": row["Locus tag"],
                    "Original gene": row["Original gene"],
                    "Refined name": row["Refined name"],
                    "Original product": row["Original product"],
                    "Enzyme class": row["Enzyme class"],
                    "Ambiguity flags": row["Ambiguity flags"],
                    "Confidence": row["Confidence"],
                }
            )
            seq = seqs.get(row["Locus tag"], "")
            if not seq:
                continue
            header = (
                f">{row['Locus tag']} cluster={row['Cluster ID']} "
                f"refined={row['Refined name'].replace(' ', '_')} "
                f"flags={row['Ambiguity flags']} product={row['Original product'].replace(' ', '_')}"
            )
            fasta.write(header + "\n")
            for idx in range(0, len(seq), 60):
                fasta.write(seq[idx : idx + 60] + "\n")
    write_tsv(
        folder / f"{prefix}_reblast_targets.tsv",
        ["Genome", "Cluster ID", "Locus tag", "Original gene", "Refined name", "Original product", "Enzyme class", "Ambiguity flags", "Confidence"],
        reblast_rows,
    )

    return cluster_rows, refinement_rows, map_rows


def main() -> None:
    root = Path(".").resolve()
    review = load_review_maps(root / "NA_cluster_review")
    genomes = [
        folder
        for folder in sorted(root.iterdir())
        if folder.is_dir()
        and folder.name != "NA_cluster_review"
        and folder.name not in SKIP_GENOMES
        and (folder / "NA_cluster_summary.tsv").exists()
        and (folder / "NA_cluster_candidates.tsv").exists()
    ]

    all_clusters: list[dict[str, str]] = []
    all_refinements: list[dict[str, str]] = []
    all_maps: list[dict[str, str]] = []
    per_genome_stats: list[dict[str, str]] = []

    for folder in genomes:
        clusters, refinements, maps = process_genome(folder, review)
        all_clusters.extend(clusters)
        all_refinements.extend(refinements)
        all_maps.extend(maps)
        per_genome_stats.append(
            {
                "Genome": folder.name,
                "Clusters": str(len(clusters)),
                "Relevant clusters": str(sum(row["Relevance (Yes/No)"] == "Yes" for row in clusters)),
                "High Priority": str(sum(row["Category"] == "High Priority" for row in clusters)),
                "Medium Priority": str(sum(row["Category"] == "Medium Priority" for row in clusters)),
                "Low Priority": str(sum(row["Category"] == "Low Priority" for row in clusters)),
                "Expanded genes": str(len(refinements)),
                "Reblast targets": str(sum(row["Ambiguity flags"] != "none" for row in refinements)),
            }
        )

    write_tsv(root / "Ref_batch_cluster_assessment.tsv", SUMMARY_FIELDS, all_clusters)
    write_tsv(root / "Ref_batch_gene_refinement_layer.tsv", REFINEMENT_FIELDS, all_refinements)
    write_tsv(root / "Ref_batch_gene_maps.tsv", ["Genome", "Cluster ID", "Category", "Functional Role", "Gene Map"], all_maps)
    write_tsv(
        root / "Ref_batch_run_summary.tsv",
        ["Genome", "Clusters", "Relevant clusters", "High Priority", "Medium Priority", "Low Priority", "Expanded genes", "Reblast targets"],
        per_genome_stats,
    )

    with (root / "Ref_batch_reblast_targets.tsv").open("w", newline="") as handle:
        fields = ["Genome", "Cluster ID", "Locus tag", "Original gene", "Refined name", "Original product", "Enzyme class", "Ambiguity flags", "Confidence"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in all_refinements:
            if row["Ambiguity flags"] != "none":
                writer.writerow({field: row.get(field, "") for field in fields})

    svg_index_dir = root / "Ref_cluster_operon_svgs"
    svg_index_dir.mkdir(exist_ok=True)
    links = []
    for folder in genomes:
        prefix = folder.name.replace(".fa", "")
        src = folder / "cluster_operon_svgs" / f"{prefix}_clusters_operons.svg"
        if src.exists():
            dest = svg_index_dir / src.name
            dest.write_text(src.read_text())
            links.append(dest.name)
    (svg_index_dir / "index.html").write_text(
        "<!doctype html><meta charset=\"utf-8\"><title>Ref cluster operon SVGs</title>"
        "<h1>Ref cluster/operon SVGs</h1><ul>"
        + "\n".join(f'<li><a href="{html.escape(name)}">{html.escape(name)}</a></li>' for name in links)
        + "</ul>"
    )

    print(f"Processed genomes: {len(genomes)}")
    print(f"Clusters: {len(all_clusters)}")
    print(f"Relevant clusters: {sum(row['Relevance (Yes/No)'] == 'Yes' for row in all_clusters)}")
    print(f"Expanded genes: {len(all_refinements)}")
    print(f"Reblast targets: {sum(row['Ambiguity flags'] != 'none' for row in all_refinements)}")
    print("Category counts:", dict(Counter(row["Category"] for row in all_clusters)))


if __name__ == "__main__":
    main()
