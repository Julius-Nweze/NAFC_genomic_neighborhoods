#!/usr/bin/env python3
# Title          : build_figure4_landscape_chart.py
# Description    : Chart candidate-neighborhood counts and multi-hit fractions by genome-source tier (Fig. 4a)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/24
# Usage          : python3 build_figure4_landscape_chart.py

import html
from pathlib import Path

OUT_SVG = Path(__file__).resolve().parent / "figure4_tier_landscape.svg"

COLOR_TOTAL = "#c9c7bc"
COLOR_MULTI = "#2a78d6"
COLOR_THREE = "#1a4a86"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#8a8980"
GRID = "#dedcd3"

# (tier label, total, multi-hit, multi-hit %, three-plus, three-plus %)
ROWS = [
    ("NCBI reference genomes", 3010, 1231, 40.9, 652, 21.7),
    ("OSPW-containing-mesocosm isolates", 1614, 575, 35.6, 274, 17.0),
    ("Plant-root MAGs", 13405, 3204, 23.9, 1169, 8.7),
]


def svg_text(x, y, text, size=12, weight="400", anchor="start", fill=TEXT_PRIMARY):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
            f'fill="{fill}">{html.escape(str(text))}</text>')


def main():
    width = 980
    left_label_w = 260
    right_margin = 90
    bar_x0 = left_label_w
    bar_max_w = width - left_label_w - right_margin
    bar_h = 16
    group_h = 3 * (bar_h + 4) + 22
    top = 92
    n = len(ROWS)
    height = top + n * group_h + 60

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Arial, Helvetica, sans-serif">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 30, "Candidate-neighborhood counts by genome-source tier", 16, "700"),
        svg_text(24, 50, "Total candidate neighborhoods, and the percent reaching multi-hit "
                 "(≥2 seed hits) or three-plus-hit (≥3) thresholds", 11.5, fill=TEXT_SECONDARY),
    ]

    grid_y0 = top - 6
    grid_y1 = top + n * group_h - 14
    for frac in (0, 25, 50, 75, 100):
        gx = bar_x0 + bar_max_w * frac / 100
        elements.append(f'<line x1="{gx:.1f}" y1="{grid_y0}" x2="{gx:.1f}" y2="{grid_y1 + 10}" '
                         f'stroke="{GRID}" stroke-width="1"/>')
        elements.append(svg_text(gx, grid_y0 - 8, f"{frac}%", 9, anchor="middle", fill=TEXT_MUTED))

    y = top
    for label, total, multi, multi_pct, three, three_pct in ROWS:
        elements.append(svg_text(left_label_w - 12, y + bar_h * 0.8, label, 12, "600", anchor="end"))
        elements.append(svg_text(left_label_w - 12, y + bar_h + 11, f"n={total:,} total neighborhoods",
                                  9, anchor="end", fill=TEXT_MUTED))
        bars = [
            ("All neighborhoods (100%)", 100.0, COLOR_TOTAL, f"{total:,}"),
            (f"Multi-hit (≥2 seed hits)", multi_pct, COLOR_MULTI, f"{multi:,} ({multi_pct:.1f}%)"),
            (f"Three-plus-hit (≥3 seed hits)", three_pct, COLOR_THREE, f"{three:,} ({three_pct:.1f}%)"),
        ]
        by = y + 20
        for sub_label, pct, color, text in bars:
            w = bar_max_w * pct / 100
            elements.append(f'<rect x="{bar_x0:.1f}" y="{by}" width="{w:.1f}" height="{bar_h}" '
                             f'rx="2" fill="{color}"/>')
            elements.append(svg_text(bar_x0 + w + 8, by + bar_h * 0.78, text, 10, "600"))
            by += bar_h + 4
        y += group_h

    legend_y = top + n * group_h + 30
    lx = 24
    for label, color in [("All neighborhoods", COLOR_TOTAL), ("Multi-hit (≥2)", COLOR_MULTI),
                          ("Three-plus-hit (≥3)", COLOR_THREE)]:
        elements.append(f'<rect x="{lx}" y="{legend_y - 12}" width="14" height="14" rx="2" fill="{color}"/>')
        elements.append(svg_text(lx + 20, legend_y, label, 11, fill=TEXT_SECONDARY))
        lx += 20 + 10 * len(label) + 30

    elements.append("</svg>")
    OUT_SVG.write_text("\n".join(elements) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
