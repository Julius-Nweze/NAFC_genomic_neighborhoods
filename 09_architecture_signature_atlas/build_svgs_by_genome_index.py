#!/usr/bin/env python3
# Title          : build_svgs_by_genome_index.py
# Description    : Write a browsable index.html of per-genome combined gene maps
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/25
# Usage          : python3 build_svgs_by_genome_index.py

import csv
import html
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHORTLIST_TSV = HERE.parent / "NA_gene_panel_top_multi_hit_clusters.tsv"

REF_BARE = {"AP042450.1", "CP010516.1", "CP010517.1", "CP063454.1", "CP063455.1", "LMAZ01000003.1"}


def tier(genome: str) -> str:
    if genome.startswith("Exp3_MAG_"):
        return "MAG"
    if genome.startswith("GCF_") or genome in REF_BARE:
        return "Ref"
    return "Colla"


def main() -> None:
    cluster_count = Counter()
    with SHORTLIST_TSV.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cluster_count[row["Genome"]] += 1

    genomes = sorted((p.stem, p.name) for p in HERE.glob("*.svg"))

    rows = [
        "<!doctype html><meta charset=\"utf-8\">",
        "<title>Per-genome candidate-cluster gene-map atlas</title>",
        "<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#fafafa}"
        "table{border-collapse:collapse;width:100%}th,td{padding:6px 10px;border-bottom:1px solid #ddd;"
        "text-align:left;font-size:13px}th{background:#eee;position:sticky;top:0}"
        "tr:hover{background:#f0f6ff}a{color:#1f4e8c;text-decoration:none}a:hover{text-decoration:underline}"
        ".tag{display:inline-block;padding:1px 8px;border-radius:10px;font-size:11px;color:#fff}"
        ".MAG{background:#1b9e77}.Ref{background:#d95f02}.Colla{background:#7570b3}</style>",
        "<h1>Per-genome candidate-cluster gene-map atlas</h1>",
        f"<p>{len(genomes)} genomes, each combined SVG showing every candidate cluster for that "
        "genome with &ge;3 independent panel gene hits (the same high-stringency shortlist used "
        "throughout this analysis, 2,095 clusters total), any functional category. This is the "
        "exhaustive per-genome view; for clusters matching a specific recurring gene-set signature "
        "across genomes instead, see <a href=\"../Architecture_gene_maps/index.html\">"
        "Architecture_gene_maps/index.html</a>.</p>",
        "<table><tr><th>Genome</th><th>Tier</th><th>Shortlisted clusters (&ge;3 seed hits)</th>"
        "<th>Gene-map atlas</th></tr>",
    ]
    for stem, fname in genomes:
        t = tier(stem)
        n = cluster_count.get(stem, 0)
        rows.append(
            "<tr>"
            f"<td>{html.escape(stem)}</td>"
            f'<td><span class="tag {t}">{t}</span></td>'
            f"<td>{n}</td>"
            f'<td><a href="{html.escape(fname)}" target="_blank">{html.escape(fname)}</a></td>'
            "</tr>"
        )
    rows.append("</table>")
    (HERE / "index.html").write_text("\n".join(rows), encoding="utf-8")
    print(f"Wrote {HERE / 'index.html'} ({len(genomes)} genomes)")


if __name__ == "__main__":
    main()
