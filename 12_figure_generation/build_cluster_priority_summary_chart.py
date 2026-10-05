#!/usr/bin/env python3
# Title          : build_cluster_priority_summary_chart.py
# Description    : Chart candidate-neighborhood strength (single, multi-hit, three-plus-hit) by genome-source tier
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/31
# Usage          : python3 build_cluster_priority_summary_chart.py [options]

"""Chart candidate-neighborhood strength (single, multi-hit, three-plus-hit) by genome-source tier."""
import argparse
import csv
import re
from pathlib import Path

COLOR_HIGH = "#b2182b"
COLOR_MEDIUM = "#ef8a62"
COLOR_LOW = "#878787"

MAX_BAR_PX = 480.0
BAR_X = 220.0
BAR_HEIGHT = 62
ROW_GAP = 142
FIRST_ROW_Y = 95
LEFT_LABEL_X = 40
CANVAS_WIDTH = 1420


def esc(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def svg_text(x, y, text, size, weight=400, anchor="start", fill="#222"):
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(text)}</text>'
    )


def load_source_map(path):
    m = {}
    with open(path) as f:
        started = False
        for line in f:
            if line.strip() == "DATA":
                started = True
                continue
            if started and line.strip():
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    m[parts[0]] = parts[2]
    return m


def normalize_genome_id(g):
    g = g.strip()
    g = re.sub(r"\.\d+$", "", g)
    g = g.replace(".", "_")
    return g


def compute_counts(clusters_tsv, source_colorstrip_txt, wanted_source):
    src_raw = load_source_map(source_colorstrip_txt)
    norm_src = {}
    for k, v in src_raw.items():
        norm_src[normalize_genome_id(re.sub(r"^g_", "", k))] = v
        norm_src[normalize_genome_id(k)] = v

    high = medium = low = 0
    with open(clusters_tsv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            src = norm_src.get(normalize_genome_id(row["Genome"]))
            if src != wanted_source:
                continue
            hits = int(row["Panel seed hit count"])
            strength = row["Cluster strength"]
            if strength == "multi_hit" and hits >= 3:
                high += 1
            elif strength == "multi_hit":
                medium += 1
            else:
                low += 1
    return {"high": high, "medium": medium, "low_relevant": low}


def build_svg(datasets, max_total):
    scale = MAX_BAR_PX / max_total
    segments = [
        ("high", COLOR_HIGH, "#fff"),
        ("medium", COLOR_MEDIUM, "#fff"),
        ("low_relevant", COLOR_LOW, "#fff"),
    ]
    elements = []
    elements.append(svg_text(40, 45, "Candidate cluster strength", 28, 700))
    elements.append(
        svg_text(
            40, 72,
            "Expanded cluster spans retain all CDS inside each candidate operon/genomic island.",
            14, 400, fill="#555",
        )
    )

    row_y = FIRST_ROW_Y
    for ds in datasets:
        counts = ds["counts"]
        total = sum(counts.values())
        label_y = row_y + BAR_HEIGHT - 24
        label = ds["label"]
        if len(label) > 20:
            breakpoints = [i for i, c in enumerate(label) if c in " -"]
            mid = min(breakpoints, key=lambda i: abs(i - len(label) / 2))
            line1, line2 = label[: mid + (1 if label[mid] == "-" else 0)], label[mid + 1 :]
            elements.append(svg_text(LEFT_LABEL_X, label_y - 10, line1.strip(), 13, 700))
            elements.append(svg_text(LEFT_LABEL_X, label_y + 9, line2.strip(), 13, 700))
        else:
            elements.append(svg_text(LEFT_LABEL_X, label_y, label, 18, 700))

        x = BAR_X
        for key, color, text_color in segments:
            count = counts[key]
            w = count * scale
            if w > 0:
                elements.append(
                    f'<rect x="{x:.1f}" y="{row_y}" width="{w:.1f}" height="{BAR_HEIGHT}" '
                    f'fill="{color}" stroke="#fff" stroke-width="1"/>'
                )
                if w > 24:
                    elements.append(
                        svg_text(x + w / 2, label_y, str(count), 15, 700, anchor="middle", fill=text_color)
                    )
            x += w

        elements.append(svg_text(825, label_y, f"{total:,} clusters", 16, 700))
        row_y += ROW_GAP

    legend_y = row_y - ROW_GAP + BAR_HEIGHT + 26
    elements.append(svg_text(LEFT_LABEL_X, legend_y, "Cluster strength", 16, 700))
    swatch_y = legend_y + BAR_HEIGHT - 48
    labels = [
        ("high", COLOR_HIGH, ">=3 panel seed hits"),
        ("medium", COLOR_MEDIUM, "2 panel seed hits (multi-hit)"),
        ("low_relevant", COLOR_LOW, "Single-hit locus"),
    ]
    positions = [220, 410, 660]
    for (key, color, label), sx in zip(labels, positions):
        elements.append(
            f'<rect x="{sx}" y="{swatch_y}" width="18" height="18" fill="{color}" stroke="#777" stroke-width="0.5"/>'
        )
        elements.append(svg_text(sx + 26, swatch_y + 14, label, 14, 400))

    footnote_y = swatch_y + 40
    footnote_parts = []
    for ds in datasets:
        c = ds["counts"]
        total = sum(c.values())
        footnote_parts.append(f"{ds['label']}: {total:,} clusters, {c['high']:,} with >=3 panel hits")
    elements.append(
        svg_text(LEFT_LABEL_X, footnote_y, ". ".join(footnote_parts) + ".", 13, 400, fill="#555")
    )

    height = footnote_y + 30
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS_WIDTH}" height="{height}" '
        f'viewBox="0 0 {CANVAS_WIDTH} {height}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(elements) + "\n</svg>\n"
    )


def write_data_tsv(path, datasets):
    label_lookup = {
        "high": ">=3 panel seed hits",
        "medium": "2 panel seed hits (multi-hit)",
        "low_relevant": "Single-hit locus",
    }
    with open(path, "w") as handle:
        handle.write("dataset\tclass\tcount\n")
        for ds in datasets:
            for key in ("high", "medium", "low_relevant"):
                handle.write(f"{ds['label']}\t{label_lookup[key]}\t{ds['counts'][key]}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clusters-tsv",
        default="/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters/NAFC_gene_panel_neighborhood_summary.tsv",
    )
    parser.add_argument(
        "--source-colorstrip",
        default="/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/itol_source_colorstrip.txt",
    )
    parser.add_argument(
        "--out-dir",
        default="/path/to/your/Figures",
    )
    args = parser.parse_args()

    datasets = [
        {"label": "MAGs", "counts": compute_counts(args.clusters_tsv, args.source_colorstrip, "Plant roots")},
        {"label": "NCBI references", "counts": compute_counts(args.clusters_tsv, args.source_colorstrip, "Genome bank")},
        {"label": "OSPW-containing-mesocosm isolates", "counts": compute_counts(args.clusters_tsv, args.source_colorstrip, "Soil-containing OSPW")},
    ]
    max_total = max(sum(ds["counts"].values()) for ds in datasets)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cluster_priority_summary.svg").write_text(build_svg(datasets, max_total))
    write_data_tsv(out_dir / "cluster_priority_summary_data.tsv", datasets)

    for ds in datasets:
        c = ds["counts"]
        total = sum(c.values())
        print(f"{ds['label']}: {total} clusters, "
              f"{c['high']} with >=3 hits, {c['medium']} with 2 hits, {c['low_relevant']} single-hit")


if __name__ == "__main__":
    main()
