#!/usr/bin/env python3
# Title          : build_figure4_composite.py
# Description    : Combine the two Fig. 4 panels into one SVG (run build_figure4_landscape_chart.py first)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/24
# Usage          : python3 build_figure4_composite.py

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PANEL_A = HERE / "figure4_tier_landscape.svg"
PANEL_B = HERE / "figure_category_multihit_summary.svg"
OUT_DIR = HERE.parent / "Figures" / "Main_Figures" / "Figure_4_cluster_landscape"

TITLE = "Fig. 4 [pending number]. Candidate-cluster landscape across genome-source tiers and functional categories"


def inner_svg(path: Path) -> str:
    text = path.read_text()
    start = text.index(">", text.index("<svg")) + 1
    end = text.rindex("</svg>")
    return text[start:end]


def dims(path: Path) -> tuple[int, int]:
    text = path.read_text()
    w = int(re.search(r'width="(\d+)', text).group(1))
    h = int(re.search(r'height="(\d+)', text).group(1))
    return w, h


def main() -> None:
    a_w, a_h = dims(PANEL_A)
    b_w, b_h = dims(PANEL_B)
    title_w_estimate = 11 * len(TITLE)  # rough px width at font-size 20, bold
    total_w = max(a_w, b_w, title_w_estimate) + 40

    y = 60
    parts = [
        f'<text x="20" y="30" font-family="Arial, Helvetica, sans-serif" font-size="20" '
        f'font-weight="700" fill="#111">{TITLE}</text>',
        f'<text x="20" y="{y - 10}" font-family="Arial, Helvetica, sans-serif" font-size="14" '
        f'font-weight="700" fill="#111">A</text>',
        f'<svg x="10" y="{y}" width="{a_w}" height="{a_h}" viewBox="0 0 {a_w} {a_h}">'
        f'{inner_svg(PANEL_A)}</svg>',
    ]
    y2 = y + a_h + 30
    parts.append(f'<text x="20" y="{y2 - 10}" font-family="Arial, Helvetica, sans-serif" '
                  f'font-size="14" font-weight="700" fill="#111">B</text>')
    parts.append(f'<svg x="10" y="{y2}" width="{b_w}" height="{b_h}" viewBox="0 0 {b_w} {b_h}">'
                  f'{inner_svg(PANEL_B)}</svg>')
    total_h = y2 + b_h + 20

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" '
           f'viewBox="0 0 {total_w} {total_h}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
           + "".join(parts) + "</svg>\n")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "composite.svg").write_text(svg)
    print(f"Wrote {OUT_DIR / 'composite.svg'} ({total_w}x{total_h})")


if __name__ == "__main__":
    main()
