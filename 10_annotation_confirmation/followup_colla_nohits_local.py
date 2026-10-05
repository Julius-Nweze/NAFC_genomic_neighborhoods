#!/usr/bin/env python3
# Title          : followup_colla_nohits_local.py
# Description    : Follow up Swiss-Prot no-hit proteins with local IS, AMR and HAMAP searches
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/24
# Usage          : python3 followup_colla_nohits_local.py --suggestions <tsv> --targets-faa <faa> --out-dir <dir> --hamap-hmm <hmm> --is-db <db> --amr-db <db>

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import Counter
from pathlib import Path


SUGGESTION_FIELDS = [
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


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def parse_fasta(path: Path) -> dict[str, tuple[str, str]]:
    seqs: dict[str, tuple[str, str]] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    seqs[header.split()[0]] = (header, "".join(seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            seqs[header.split()[0]] = (header, "".join(seq_lines))
    return seqs


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def parse_blast6(path: Path) -> dict[str, dict[str, str]]:
    top: dict[str, dict[str, str]] = {}
    if not path.exists():
        return top
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) != 8:
                continue
            qseqid, sacc, stitle, pident, length, evalue, bitscore, qcovs = row
            current = {
                "qseqid": qseqid,
                "sacc": sacc,
                "stitle": stitle,
                "pident": pident,
                "length": length,
                "evalue": evalue,
                "bitscore": bitscore,
                "qcovs": qcovs,
            }
            prev = top.get(qseqid)
            if prev is None or (as_float(evalue), -as_float(bitscore), -as_float(qcovs)) < (
                as_float(prev["evalue"]),
                -as_float(prev["bitscore"]),
                -as_float(prev["qcovs"]),
            ):
                top[qseqid] = current
    return top


def parse_hmmscan_tblout(path: Path) -> dict[str, dict[str, str]]:
    top: dict[str, dict[str, str]] = {}
    if not path.exists():
        return top
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            target_name = parts[0]
            target_accession = parts[1]
            query_name = parts[2]
            query_accession = parts[3]
            full_evalue = parts[4]
            full_score = parts[5]
            full_bias = parts[6]
            best_dom_evalue = parts[7]
            row = {
                "query_name": query_name,
                "target_name": target_name,
                "target_accession": target_accession,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "full_bias": full_bias,
                "best_domain_evalue": best_dom_evalue,
            }
            prev = top.get(query_name)
            if prev is None or (as_float(full_evalue), -as_float(full_score)) < (
                as_float(prev["full_evalue"]),
                -as_float(prev["full_score"]),
            ):
                top[query_name] = row
    return top


def decide_followup(is_hit: dict[str, str] | None, amr_hit: dict[str, str] | None, hamap_hit: dict[str, str] | None) -> tuple[str, str]:
    if is_hit and as_float(is_hit["bitscore"]) >= 80 and as_float(is_hit["qcovs"]) >= 70:
        return ("IS_supported_mobile_element", is_hit["stitle"])
    if amr_hit and as_float(amr_hit["bitscore"]) >= 80 and as_float(amr_hit["qcovs"]) >= 70:
        return ("AMR_supported_function", amr_hit["stitle"])
    if hamap_hit and as_float(hamap_hit["full_score"]) >= 40 and as_float(hamap_hit["full_evalue"]) <= 1e-5:
        return ("HAMAP_supported_family", hamap_hit["target_name"])
    if hamap_hit and as_float(hamap_hit["full_score"]) >= 25 and as_float(hamap_hit["full_evalue"]) <= 1e-3:
        return ("weak_HAMAP_family_signal", hamap_hit["target_name"])
    return ("still_unresolved_after_local_followup", "")


def main() -> None:
    parser = argparse.ArgumentParser(description="Local follow-up for Colla Swiss-Prot no-hit targets using IS/AMR BLAST and HAMAP HMM support.")
    parser.add_argument("--suggestions", type=Path, required=True)
    parser.add_argument("--targets-faa", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="Colla_batch")
    parser.add_argument("--hamap-hmm", type=Path, required=True)
    parser.add_argument("--is-db", type=Path, required=True)
    parser.add_argument("--amr-db", type=Path, required=True)
    parser.add_argument("--threads", default="8")
    parser.add_argument("--force", action="store_true", help="Rerun BLAST/HMM searches even if raw output files already exist.")
    args = parser.parse_args()

    rows = read_tsv(args.suggestions)
    nohit_rows = [row for row in rows if row.get("suggested_action") == "retain_flag_for_TrEMBL_or_HMM_review"]
    nohit_ids = {row["Locus tag"] for row in nohit_rows}

    seqs = parse_fasta(args.targets_faa)
    records = [seqs[qid] for qid in nohit_ids if qid in seqs]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    nohit_faa = args.out_dir / f"{args.prefix}_nohit_targets.faa"
    write_fasta(nohit_faa, records)

    is_out = args.out_dir / f"{args.prefix}_nohit_vs_is.blastp.tsv"
    amr_out = args.out_dir / f"{args.prefix}_nohit_vs_amr.blastp.tsv"
    hamap_tbl = args.out_dir / f"{args.prefix}_nohit_vs_hamap.tbl"

    if records and (args.force or not is_out.exists() or is_out.stat().st_size == 0):
        run(
            [
                "blastp",
                "-query",
                str(nohit_faa),
                "-db",
                str(args.is_db),
                "-evalue",
                "1e-5",
                "-max_target_seqs",
                "5",
                "-num_threads",
                str(args.threads),
                "-outfmt",
                "6 qseqid sacc stitle pident length evalue bitscore qcovs",
                "-out",
                str(is_out),
            ]
        )
    if records and (args.force or not amr_out.exists() or amr_out.stat().st_size == 0):
        run(
            [
                "blastp",
                "-query",
                str(nohit_faa),
                "-db",
                str(args.amr_db),
                "-evalue",
                "1e-5",
                "-max_target_seqs",
                "5",
                "-num_threads",
                str(args.threads),
                "-outfmt",
                "6 qseqid sacc stitle pident length evalue bitscore qcovs",
                "-out",
                str(amr_out),
            ]
        )
    if records and (args.force or not hamap_tbl.exists() or hamap_tbl.stat().st_size == 0):
        run(
            [
                "hmmscan",
                "--cpu",
                str(args.threads),
                "--tblout",
                str(hamap_tbl),
                "-o",
                "/dev/null",
                str(args.hamap_hmm),
                str(nohit_faa),
            ]
        )

    is_hits = parse_blast6(is_out)
    amr_hits = parse_blast6(amr_out)
    hamap_hits = parse_hmmscan_tblout(hamap_tbl)

    out_rows: list[dict[str, str]] = []
    counts = Counter()
    for row in nohit_rows:
        locus = row["Locus tag"]
        is_hit = is_hits.get(locus)
        amr_hit = amr_hits.get(locus)
        hamap_hit = hamap_hits.get(locus)
        category, label = decide_followup(is_hit, amr_hit, hamap_hit)
        counts[category] += 1
        out = {field: row.get(field, "") for field in SUGGESTION_FIELDS}
        out["followup_category"] = category
        out["recommended_product_or_family"] = label
        out["is_top_hit"] = is_hit["stitle"] if is_hit else ""
        out["is_pident"] = is_hit["pident"] if is_hit else ""
        out["is_qcovs"] = is_hit["qcovs"] if is_hit else ""
        out["is_bitscore"] = is_hit["bitscore"] if is_hit else ""
        out["amr_top_hit"] = amr_hit["stitle"] if amr_hit else ""
        out["amr_pident"] = amr_hit["pident"] if amr_hit else ""
        out["amr_qcovs"] = amr_hit["qcovs"] if amr_hit else ""
        out["amr_bitscore"] = amr_hit["bitscore"] if amr_hit else ""
        out["hamap_top_model"] = hamap_hit["target_name"] if hamap_hit else ""
        out["hamap_evalue"] = hamap_hit["full_evalue"] if hamap_hit else ""
        out["hamap_score"] = hamap_hit["full_score"] if hamap_hit else ""
        out_rows.append(out)

    out_rows.sort(key=lambda r: (r["followup_category"], r["Genome"], r["Cluster ID"], r["Locus tag"]))

    fields = SUGGESTION_FIELDS + [
        "followup_category",
        "recommended_product_or_family",
        "is_top_hit",
        "is_pident",
        "is_qcovs",
        "is_bitscore",
        "amr_top_hit",
        "amr_pident",
        "amr_qcovs",
        "amr_bitscore",
        "hamap_top_model",
        "hamap_evalue",
        "hamap_score",
    ]
    followup_tsv = args.out_dir / f"{args.prefix}_nohit_followup_local_evidence.tsv"
    write_tsv(followup_tsv, fields, out_rows)

    summary = args.out_dir / f"{args.prefix}_nohit_followup_local_evidence_summary.md"
    with summary.open("w") as handle:
        handle.write(f"# {args.prefix} local follow-up of Swiss-Prot no-hit targets\n\n")
        handle.write(f"- No-hit targets reviewed: {len(nohit_rows)}\n")
        handle.write(f"- FASTA records extracted: {len(records)}\n")
        handle.write("- Local follow-up resources: Prokka `IS` BLAST DB, Prokka `AMR` BLAST DB, and `HAMAP.hmm`\n")
        handle.write("- Note: no local UniRef or TrEMBL database was available in this environment, so broader UniRef/TrEMBL confirmation remains pending if those databases are installed later.\n\n")
        for key, value in sorted(counts.items()):
            handle.write(f"- {key}: {value}\n")

    print(f"No-hit targets reviewed: {len(nohit_rows)}")
    for key, value in sorted(counts.items()):
        print(f"{key}\t{value}")
    print(f"Wrote {followup_tsv}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
