#!/usr/bin/env python3
# Title          : fix_KO_proteins_fasta_panel.py
# Description    : Reconcile the curated reference-panel FASTA (KO_proteins.fasta) with the 667-gene list in Gene_info.xlsx
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/31
# Usage          : python3 fix_KO_proteins_fasta_panel.py

from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

import openpyxl

GENE_INFO = Path("/path/to/your/Ref_genomes/Gene_info.xlsx")
FASTA_PATHS = [
    Path("/path/to/your/Ref_genomes/KO/KO_proteins.fasta"),
    Path("/path/to/your/Ref_genomes/Analysis/KO/KO_proteins.fasta"),
]
REPORT_PATH = Path(
    "/path/to/your/Ref_genomes/KO/KO_proteins_panel_correction_report_2026-07-31.tsv"
)
UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def load_gene_rows():
    wb = openpyxl.load_workbook(GENE_INFO, read_only=True)
    ws = wb["Code"]
    rows = list(ws.iter_rows(values_only=True))[1:]
    return [
        {
            "idx": i,
            "cat": cat,
            "gene": (gene or "").strip(),
            "name": name,
            "ko": (ko or "").strip(),
            "uniprot": (uni or "").strip(),
            "ec": ec,
        }
        for i, (cat, gene, name, ko, uni, ec) in enumerate(rows)
    ]


def load_fasta(path):
    seqs = []
    header, seq_lines = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq_lines)))
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line)
    if header is not None:
        seqs.append((header, "".join(seq_lines)))
    parsed = []
    for h, s in seqs:
        toks = [t.strip() for t in re.split(r"[|]", h)]
        gn = re.search(r"GN=(\S+)", h)
        if gn:
            toks.append(gn.group(1))
        parsed.append({"header": h, "seq": s, "tokens": toks, "consumed": False})
    return parsed


def build_indexes(fasta_parsed):
    exact_index, gene_index = {}, {}
    for fp in fasta_parsed:
        for t in fp["tokens"]:
            exact_index.setdefault(t, []).append(fp)
            tclean = re.split(r"[\s_]", t)[0].lower()
            if 1 < len(tclean) <= 20:
                gene_index.setdefault(tclean, []).append(fp)
    return exact_index, gene_index


def match_rows(gene_rows, fasta_parsed):
    exact_index, gene_index = build_indexes(fasta_parsed)
    row_to_seq = {}
    for gr in gene_rows:  # pass 1: exact token match (any identifier type)
        if gr["uniprot"] and gr["uniprot"] in exact_index:
            cands = [c for c in exact_index[gr["uniprot"]] if not c["consumed"]]
            if cands:
                cands[0]["consumed"] = True
                row_to_seq[gr["idx"]] = cands[0]
    unresolved = []
    for gr in gene_rows:  # pass 2: gene-symbol fallback, consume-once
        if gr["idx"] in row_to_seq:
            continue
        cands = [c for c in gene_index.get(gr["gene"].lower(), []) if not c["consumed"]]
        if cands:
            cands[0]["consumed"] = True
            row_to_seq[gr["idx"]] = cands[0]
        else:
            unresolved.append(gr)
    return row_to_seq, unresolved


def fetch_uniprot_or_uniparc(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=15) as resp:
        text = resp.read().decode()
    if text.startswith(">"):
        return text.strip()
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=15) as resp:
        j = json.loads(resp.read().decode())
    upi = j.get("extraAttributes", {}).get("uniParcId")
    if not upi:
        raise RuntimeError(f"{acc}: no UniProt sequence and no UniParc fallback")
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniparc/{upi}.fasta", timeout=15) as resp:
        text2 = resp.read().decode()
    lines = text2.strip().split("\n")
    return f">{acc}|{lines[0][1:]}\n" + "\n".join(lines[1:])


def main():
    gene_rows = load_gene_rows()
    fasta_parsed = load_fasta(FASTA_PATHS[0])
    row_to_seq, unresolved = match_rows(gene_rows, fasta_parsed)

    fetched = {}
    for gr in unresolved:
        acc = gr["uniprot"]
        if not (acc and UNIPROT_RE.match(acc)):
            raise RuntimeError(f"row {gr} has no fetchable accession - resolve manually")
        fetched[acc] = fetch_uniprot_or_uniparc(acc)
        time.sleep(0.12)

    out_records = []
    for gr in gene_rows:
        if gr["idx"] in row_to_seq:
            fp = row_to_seq[gr["idx"]]
            out_records.append((gr, fp["header"], fp["seq"], "kept_from_existing_fasta"))
        else:
            lines = fetched[gr["uniprot"]].strip().split("\n")
            out_records.append((gr, lines[0][1:], "".join(lines[1:]), "fetched_from_uniprot_or_uniparc"))

    assert len(out_records) == 667, f"expected 667 records, got {len(out_records)}"
    assert len({h for _, h, _, _ in out_records}) == 667, "duplicate header assigned to two rows"

    for path in FASTA_PATHS:
        with open(path, "w") as f:
            for _, header, seq, _ in out_records:
                f.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
        print(f"wrote 667 records to {path}")

    removed = [fp for fp in fasta_parsed if not fp["consumed"]]
    with open(REPORT_PATH, "w") as f:
        f.write("action\tcategory\tgene\tuniprot\tsource_or_note\n")
        for gr, header, _, source in out_records:
            f.write(f"kept_or_added\t{gr['cat']}\t{gr['gene']}\t{gr['uniprot']}\t{source}\n")
        for fp in removed:
            f.write(f"removed_not_in_gene_info\t\t\t\t{fp['header']}\n")
    print(f"report written: {REPORT_PATH} ({len(removed)} sequences removed)")


if __name__ == "__main__":
    main()
