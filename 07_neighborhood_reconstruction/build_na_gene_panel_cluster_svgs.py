#!/usr/bin/env python3
# Title          : build_na_gene_panel_cluster_svgs.py
# Description    : Render one SVG gene map per multi-hit candidate neighborhood and write an index.html
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/22
# Usage          : python3 build_na_gene_panel_cluster_svgs.py

from __future__ import annotations

import csv
import html
import math
import sys
from collections import defaultdict
from pathlib import Path

REF_CLUSTER_SCRIPT_DIR = Path(__file__).resolve().parent.parent / "shared_dependencies"
sys.path.insert(0, str(REF_CLUSTER_SCRIPT_DIR))
from generate_ref_cluster_outputs import (  # type: ignore  # noqa: E402
    COLOR,
    arrow_points,
    svg_text,
    wrap_label,
)

STRENGTH_COLOR = {
    "multi_hit": "#b2182b",
    "single_hit": "#878787",
}

OPERON_DIR = Path("/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters")
SHORTLIST_PATH = OPERON_DIR / "NA_gene_panel_top_multi_hit_clusters.tsv"
MEMBERS_PATH = OPERON_DIR / "NA_gene_panel_cluster_members.tsv"
SVG_DIR = OPERON_DIR / "svgs"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def render_cluster_block(cluster: dict[str, str], records: list[dict[str, str]], y: int, total_width: int = 1700) -> tuple[list[str], int]:
    """Return (svg element strings, vertical space consumed) for one cluster, anchored at y."""
    left, label_width, right = 32, 420, 36
    plot_x = left + label_width
    plot_width = total_width - plot_x - right
    gene_height = 22
    top = y
    base = top + 62
    block_height = 62 + gene_height + 34

    start = min(int(row["Start"]) for row in records)
    end = max(int(row["End"]) for row in records)
    scale = plot_width / max(end - start + 1, 1)

    strength = cluster["Cluster strength"]
    strength_label = "Multi-hit cluster (>=2 independent panel gene hits)" if strength == "multi_hit" else "Single-hit locus"

    elements = [
        f'<rect x="18" y="{top - 18}" width="{total_width - 36}" height="{block_height}" rx="6" fill="#ffffff" stroke="#dddddd"/>',
        f'<rect x="18" y="{top - 18}" width="7" height="{block_height}" fill="{STRENGTH_COLOR.get(strength, "#666")}"/>',
        svg_text(left, top, f"{cluster['Cluster ID']}", 16, "700"),
        svg_text(left, top + 18, f"{cluster['Genome']} | {cluster['Contig']} | {strength_label}", 11, "700", fill=STRENGTH_COLOR.get(strength, "#666")),
    ]
    detail_lines = wrap_label(f"Panel genes: {cluster['Panel genes']}", 68)[:2]
    detail_lines += wrap_label(f"Panel categories: {cluster['Panel categories']}", 68)[:2]
    for idx, line in enumerate(detail_lines):
        elements.append(svg_text(left, top + 34 + idx * 12, line, 9, fill="#444"))

    elements.append(
        f'<line x1="{plot_x}" y1="{base + gene_height / 2}" x2="{plot_x + plot_width}" y2="{base + gene_height / 2}" stroke="#d0d0d0"/>'
    )
    for row in records:
        row_start = int(row["Start"])
        row_end = int(row["End"])
        gene_x = plot_x + (row_start - start) * scale
        gene_width = max((row_end - row_start + 1) * scale, 18)
        if gene_x + gene_width > plot_x + plot_width:
            gene_x = plot_x + plot_width - gene_width
        is_seed = row["Is panel seed hit"] == "Yes"
        fill = COLOR.get(row["Enzyme class"], "#ccc")
        stroke = "#111" if is_seed else "#888"
        stroke_width = "2.2" if is_seed else "0.8"
        gene_label = row["Gene"] or row["Locus tag"]
        title_bits = [row["Locus tag"], gene_label, row["Enzyme class"], row["Product"]]
        if is_seed:
            title_bits.append(f"PANEL HIT: {row['Panel gene']} ({row['Panel category']} / {row['Panel sub-category']})")
            title_bits.append(f"pident={row['Panel percent identity']} qcovs={row['Panel query coverage']}")
        title = " | ".join(bit for bit in title_bits if bit)
        elements.append(
            f'<polygon points="{arrow_points(gene_x, base, gene_width, gene_height, row["Strand"])}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"><title>{html.escape(title)}</title></polygon>'
        )
        label = f"{gene_label} ({row['Panel gene']})" if is_seed and row["Panel gene"] else gene_label
        if len(label) > 24:
            label = label[:21] + "..."
        tx = gene_x + gene_width / 2
        ty = base + gene_height + 13
        weight = "700" if is_seed else "400"
        if len(records) > 12:
            elements.append(
                f'<text x="{tx:.1f}" y="{ty:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="9" '
                f'font-weight="{weight}" text-anchor="end" transform="rotate(-45 {tx:.1f} {ty:.1f})" '
                f'fill="{"#111" if is_seed else "#555"}">{html.escape(label)}</text>'
            )
        else:
            elements.append(svg_text(tx, ty, label, 9, weight, anchor="middle", fill="#111" if is_seed else "#555"))

    elements.append(svg_text(plot_x, top + 10, f"{start:,} bp", 9, fill="#666"))
    elements.append(svg_text(plot_x + plot_width, top + 10, f"{end:,} bp", 9, anchor="end", fill="#666"))

    return elements, block_height + 24


def wrap_svg_document(inner_elements: list[str], total_width: int, height: int) -> str:
    doc = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{height}" viewBox="0 0 {total_width} {height}">',
        '<rect width="100%" height="100%" fill="#f8f8f6"/>',
    ]
    doc.extend(inner_elements)
    doc.append("</svg>")
    return "\n".join(doc) + "\n"


def legend_block(y: int, total_width: int = 1700) -> tuple[list[str], int]:
    legend_rows = math.ceil(len(COLOR) / 4)
    elements = [
        svg_text(24, y, "Legend (bold outline = reference-panel seed hit; thin outline = expanded genomic context)", 10, "700"),
    ]
    x0, y0, cell_width = 24, y + 10, 400
    for idx, (label, color) in enumerate(COLOR.items()):
        col = idx % 4
        row = idx // 4
        x = x0 + col * cell_width
        yy = y0 + row * 22
        elements.append(f'<rect x="{x}" y="{yy}" width="16" height="12" fill="{color}" stroke="#555" stroke-width="0.6"/>')
        elements.append(svg_text(x + 22, yy + 10, label, 9, fill="#333"))
    return elements, 10 + legend_rows * 22 + 10


def render_cluster_svg(cluster: dict[str, str], records: list[dict[str, str]], total_width: int = 1700) -> str:
    elements, consumed = render_cluster_block(cluster, records, y=24, total_width=total_width)
    legend_elements, legend_height = legend_block(24 + consumed, total_width=total_width)
    elements.extend(legend_elements)
    return wrap_svg_document(elements, total_width, 24 + consumed + legend_height)


def render_genome_svg(genome: str, clusters: list[dict[str, str]], members_by_cluster: dict[str, list[dict[str, str]]], total_width: int = 1700) -> str:
    elements = [
        svg_text(24, 34, f"{genome}: reference-panel operon/cluster gene maps", 22, "700"),
        svg_text(24, 56, f"{len(clusters)} shortlisted clusters (>=3 independent panel gene hits) for this genome", 12, fill="#555"),
        svg_text(24, 74, "Arrow direction shows CDS strand. Color shows enzyme class. Bold outline = reference-panel seed hit.", 11, fill="#555"),
    ]
    y = 100
    for cluster in clusters:
        records = sorted(members_by_cluster[cluster["Cluster ID"]], key=lambda row: (int(row["Start"]), int(row["End"])))
        block_elements, consumed = render_cluster_block(cluster, records, y=y, total_width=total_width)
        elements.extend(block_elements)
        y += consumed
    legend_elements, legend_height = legend_block(y, total_width=total_width)
    elements.extend(legend_elements)
    return wrap_svg_document(elements, total_width, y + legend_height)


def slugify(cluster_id: str) -> str:
    return cluster_id.replace("/", "_")


def main() -> None:
    SVG_DIR.mkdir(parents=True, exist_ok=True)
    genome_dir = SVG_DIR.parent / "svgs_by_genome"
    genome_dir.mkdir(parents=True, exist_ok=True)

    shortlist = read_tsv(SHORTLIST_PATH)
    wanted_ids = {row["Cluster ID"] for row in shortlist}
    print(f"Rendering SVGs for {len(wanted_ids)} shortlisted clusters ...")

    members_by_cluster: dict[str, list[dict[str, str]]] = defaultdict(list)
    with MEMBERS_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["Cluster ID"] in wanted_ids:
                members_by_cluster[row["Cluster ID"]].append(row)

    index_rows = []
    clusters_by_genome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for rank, cluster in enumerate(shortlist, start=1):
        cluster_id = cluster["Cluster ID"]
        records = sorted(members_by_cluster[cluster_id], key=lambda row: (int(row["Start"]), int(row["End"])))
        svg_text_content = render_cluster_svg(cluster, records)
        file_name = f"{slugify(cluster_id)}.svg"
        (SVG_DIR / file_name).write_text(svg_text_content, encoding="utf-8")
        index_rows.append((rank, cluster, file_name))
        clusters_by_genome[cluster["Genome"]].append(cluster)
        if rank % 200 == 0:
            print(f"  ... {rank}/{len(shortlist)}")

    print(f"Rendering {len(clusters_by_genome)} per-genome combined SVGs (all shortlisted clusters for that genome, one file) ...")
    genome_file_by_name: dict[str, str] = {}
    for genome, clusters in clusters_by_genome.items():
        clusters_sorted = sorted(clusters, key=lambda row: (row["Contig"], int(row["Start"])))
        svg_text_content = render_genome_svg(genome, clusters_sorted, members_by_cluster)
        file_name = f"{slugify(genome)}.svg"
        (genome_dir / file_name).write_text(svg_text_content, encoding="utf-8")
        genome_file_by_name[genome] = file_name

    index_lines = [
        "<!doctype html><meta charset=\"utf-8\">",
        "<title>Reference-panel operon/cluster gene maps</title>",
        "<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#fafafa}"
        "table{border-collapse:collapse;width:100%}th,td{padding:6px 10px;border-bottom:1px solid #ddd;"
        "text-align:left;font-size:13px}th{background:#eee;position:sticky;top:0}"
        "tr:hover{background:#f0f6ff}a{color:#1f4e8c;text-decoration:none}a:hover{text-decoration:underline}</style>",
        "<h1>Reference-panel operon/cluster gene maps</h1>",
        f"<p>{len(index_rows)} multi-hit clusters with &ge;3 independent reference-panel gene hits across "
        f"{len(clusters_by_genome)} genomes, sorted by seed-hit count. Bold-outlined arrows in each SVG are "
        "the panel seed hits; thin-outlined arrows are expanded genomic context (any CDS in the same span). "
        "The \"Genome (all clusters)\" column links to one combined SVG per genome "
        "(../svgs_by_genome/) with every shortlisted cluster for that genome stacked in one file.</p>",
        "<table><tr><th>#</th><th>Cluster</th><th>Genome</th><th>Genome (all clusters)</th><th>Panel seed hits</th>"
        "<th>Total CDS</th><th>Panel genes</th><th>Panel categories</th></tr>",
    ]
    for rank, cluster, file_name in index_rows:
        genome = cluster["Genome"]
        genome_file = genome_file_by_name.get(genome, "")
        genome_link = (
            f'<a href="../svgs_by_genome/{html.escape(genome_file)}" target="_blank">{html.escape(genome)} ({len(clusters_by_genome[genome])})</a>'
            if genome_file else html.escape(genome)
        )
        index_lines.append(
            "<tr>"
            f"<td>{rank}</td>"
            f'<td><a href="{html.escape(file_name)}" target="_blank">{html.escape(cluster["Cluster ID"])}</a></td>'
            f"<td>{html.escape(genome)}</td>"
            f"<td>{genome_link}</td>"
            f"<td>{html.escape(cluster['Panel seed hit count'])}</td>"
            f"<td>{html.escape(cluster['Total CDS in span'])}</td>"
            f"<td>{html.escape(cluster['Panel genes'])}</td>"
            f"<td>{html.escape(cluster['Panel categories'])}</td>"
            "</tr>"
        )
    index_lines.append("</table>")
    (SVG_DIR / "index.html").write_text("\n".join(index_lines), encoding="utf-8")

    print(f"Wrote {len(index_rows)} per-cluster SVGs -> {SVG_DIR}")
    print(f"Wrote {len(clusters_by_genome)} per-genome combined SVGs -> {genome_dir}")
    print(f"index.html -> {SVG_DIR / 'index.html'}")


if __name__ == "__main__":
    main()
