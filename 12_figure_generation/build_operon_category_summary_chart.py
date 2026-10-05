#!/usr/bin/env python3
# Title          : build_operon_category_summary_chart.py
# Description    : Chart the percent of panel hits in multi-hit neighborhoods per functional category
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/22
# Usage          : python3 build_operon_category_summary_chart.py

from __future__ import annotations

import csv
import html
from collections import Counter
from pathlib import Path

OPERON_DIR = Path("/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters")
SUMMARY_PATH = OPERON_DIR / "NAFC_gene_panel_neighborhood_summary.tsv"
MEMBERS_PATH = OPERON_DIR / "NAFC_gene_panel_neighborhood_members.tsv"
OUT_SVG = OPERON_DIR / "figure_category_multihit_summary.svg"
OUT_TSV = OPERON_DIR / "figure_category_multihit_summary_data.tsv"

COLOR_MULTI = "#2a78d6"
COLOR_SINGLE = "#eb6834"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#8a8980"
GRID = "#dedcd3"


def svg_text(x: float, y: float, text: str, size: int = 12, weight: str = "400", anchor: str = "start", fill: str = TEXT_PRIMARY) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{html.escape(str(text))}</text>'
    )


def main() -> None:
    with SUMMARY_PATH.open(newline="", encoding="utf-8") as handle:
        strength_by_cluster = {row["Cluster ID"]: row["Cluster strength"] for row in csv.DictReader(handle, delimiter="\t")}

    total_counter: Counter = Counter()
    multi_counter: Counter = Counter()
    with MEMBERS_PATH.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Is panel seed hit"] != "Yes":
                continue
            categories = [c for c in row["Panel category"].split("; ") if c]
            if not categories:
                continue
            strength = strength_by_cluster.get(row["Cluster ID"], "single_hit")
            for category in categories:
                total_counter[category] += 1
                if strength == "multi_hit":
                    multi_counter[category] += 1

    rows = []
    for category, total in total_counter.items():
        multi = multi_counter.get(category, 0)
        rows.append((category, total, multi, 100.0 * multi / total if total else 0.0))
    rows.sort(key=lambda r: r[3], reverse=True)

    # --- layout ---
    width = 980
    left_label_w = 190
    right_margin = 90
    bar_x0 = left_label_w
    bar_max_w = width - left_label_w - right_margin
    bar_h = 22
    row_h = 40
    top = 92
    n = len(rows)
    height = top + n * row_h + 70

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Arial, Helvetica, sans-serif">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 30, "Reference-panel gene hits are predominantly organized into multi-gene clusters", 16, "700"),
        svg_text(24, 50, "Percent of panel seed hits per category found inside a multi-hit genomic cluster (≥2 independent panel genes co-occurring)", 11.5, fill=TEXT_SECONDARY),
    ]

    # gridlines at 0/25/50/75/100%
    grid_y0 = top - 6
    grid_y1 = top + n * row_h - (row_h - bar_h) - 6
    for frac in (0, 25, 50, 75, 100):
        gx = bar_x0 + bar_max_w * frac / 100
        elements.append(f'<line x1="{gx:.1f}" y1="{grid_y0}" x2="{gx:.1f}" y2="{grid_y1 + 10}" stroke="{GRID}" stroke-width="1"/>')
        elements.append(svg_text(gx, grid_y0 - 8, f"{frac}%", 9, anchor="middle", fill=TEXT_MUTED))

    y = top
    for category, total, multi, pct in rows:
        single = total - multi
        gap = 2
        multi_w = bar_max_w * (multi / total)
        single_w = max(bar_max_w * (single / total) - gap, 0)
        elements.append(svg_text(left_label_w - 12, y + bar_h * 0.72, category, 12, "600", anchor="end", fill=TEXT_PRIMARY))
        elements.append(svg_text(left_label_w - 12, y + bar_h + 11, f"n={total:,} panel hits", 9, anchor="end", fill=TEXT_MUTED))

        # square data-end where segments meet the gap; rounded data-end on the bar's outer (right) end
        elements.append(
            f'<rect x="{bar_x0:.1f}" y="{y}" width="{max(multi_w - gap, 0):.1f}" height="{bar_h}" fill="{COLOR_MULTI}"/>'
        )
        seg2_x = bar_x0 + multi_w
        elements.append(
            f'<path d="M {seg2_x:.1f} {y} h {max(single_w - 4, 0):.1f} q 4 0 4 4 v {bar_h - 8} q 0 4 -4 4 h -{max(single_w - 4, 0):.1f} z" fill="{COLOR_SINGLE}"/>'
        )
        elements.append(svg_text(bar_x0 + bar_max_w + 10, y + bar_h * 0.72, f"{pct:.1f}%", 12, "700", fill=TEXT_PRIMARY))
        y += row_h

    # legend
    legend_y = top + n * row_h + 26
    elements.append(f'<rect x="24" y="{legend_y - 12}" width="14" height="14" rx="2" fill="{COLOR_MULTI}"/>')
    elements.append(svg_text(44, legend_y, "Multi-hit cluster (≥2 independent panel gene hits)", 11, fill=TEXT_SECONDARY))
    elements.append(f'<rect x="380" y="{legend_y - 12}" width="14" height="14" rx="2" fill="{COLOR_SINGLE}"/>')
    elements.append(svg_text(400, legend_y, "Single-hit locus", 11, fill=TEXT_SECONDARY))

    elements.append("</svg>")
    OUT_SVG.write_text("\n".join(elements) + "\n", encoding="utf-8")

    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "total_panel_seed_hits", "multi_hit_seed_hits", "percent_multi_hit"])
        for category, total, multi, pct in rows:
            writer.writerow([category, total, multi, round(pct, 1)])

    print(f"Wrote {OUT_SVG}")
    print(f"Wrote {OUT_TSV}")
    for category, total, multi, pct in rows:
        print(f"  {category}: {total} seed hits, {multi} multi-hit ({pct:.1f}%)")


if __name__ == "__main__":
    main()
