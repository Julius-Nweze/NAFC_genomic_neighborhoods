#!/usr/bin/env python3
# Title          : panel_style.py
# Description    : Shared per-panel header style for gene-map composites (imported module)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/25
# Usage          : python3 panel_style.py (imported by other scripts; not run directly)

import html

PANEL_W = 1700
HEADER_ROW_H = 30      # per-panel: source/genome/category text row
BAR_H = 5               # red divider bar
BP_ROW_H = 20           # start/end bp row
ARROWS_TOP = 80         # y in the source per-cluster SVG where the arrow diagram begins
ARROWS_BOTTOM = 166     # y where the legend starts (verified constant across all source panels)
ARROWS_H = ARROWS_BOTTOM - ARROWS_TOP

PANEL_HEADER_H = HEADER_ROW_H + BAR_H + BP_ROW_H  # total header block height per panel
PANEL_TOTAL_H = PANEL_HEADER_H + ARROWS_H          # total height of one panel, header + arrows

BAR_COLOR = "#b2182b"
SOURCE_CANVAS_H = 274  # full height of every original per-cluster source SVG
LEGEND_H = SOURCE_CANVAS_H - ARROWS_BOTTOM  # = 108, legend block height

TIER_FULL_LABEL = {
    "MAG": "Plant-root MAGs",
    "Ref": "NCBI reference genomes",
    "Colla": "OSPW-containing-mesocosm isolate genomes",
}

COLUMN_HEADER_H = 34  # once-per-composite "Source | Genome name | Category" title row


def column_header_svg(y: int) -> str:
    return (
        f'<text x="20" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="20" '
        f'font-weight="700" fill="#111">Source</text>\n'
        f'<text x="{PANEL_W * 0.42:.0f}" y="{y}" text-anchor="middle" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="700" '
        f'fill="#111">Genome name</text>\n'
        f'<text x="{PANEL_W - 20}" y="{y}" text-anchor="end" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="700" '
        f'fill="#111">Category</text>\n'
    )


def split_off_header(content: str) -> str:
    """Discard the original per-cluster SVG's header block (cluster-ID title, genome/strength
    line, panel-genes list, panel-categories list): everything up to and including the
    horizontal baseline line the gene arrows sit on. The new header is rebuilt from the cluster
    TSV row instead of parsed back out of this text."""
    marker = 'stroke="#d0d0d0"'
    idx = content.index(marker)
    cut = content.rindex("<line", 0, idx)
    return content[cut:]


def split_off_legend(content: str) -> tuple[str, str]:
    """Physically cut the legend block out of one panel's content (viewBox-only clipping is not
    reliably respected across SVG renderers)."""
    marker = ">Legend (bold outline"
    idx = content.index(marker)
    cut = content.rindex("<text", 0, idx)
    return content[:cut], content[cut:]


def arrows_and_legend(content: str) -> tuple[str, str]:
    """Header stripped off the front, legend split off the back. Returns (arrows_only,
    legend_fragment): arrows_only is just the gene-arrow diagram + gene labels + start/end bp
    text; legend_fragment is kept once per composite and discarded for every other panel."""
    no_header = split_off_header(content)
    no_legend, legend_rest = split_off_legend(no_header)
    return no_legend, legend_rest


def panel_svg(row: dict, tier_code: str, y_top: float, arrows_content: str,
              prefix: str = "") -> str:
    """One panel: red bar + Source/Genome name/Category row, then the cropped arrow diagram
    (embedded at its original y-coordinates via a viewBox offset, so nothing needs rewriting).
    `prefix` (e.g. "A. " or "F. Gentisate: ") is prepended to the Genome name column, for the
    curated composites' lettered panels."""
    genome = row["Genome"]
    contig = row["Contig"]
    cats = row["Panel categories"]
    source_label = TIER_FULL_LABEL[tier_code]
    bar_y = y_top + HEADER_ROW_H
    arrows_y = bar_y + BAR_H
    return (
        f'<rect x="0" y="{bar_y:.1f}" width="{PANEL_W}" height="{BAR_H}" fill="{BAR_COLOR}"/>\n'
        f'<text x="0" y="{y_top + 16:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="11" fill="#444">{html.escape(source_label)}</text>\n'
        f'<text x="{PANEL_W * 0.42:.0f}" y="{y_top + 16:.1f}" text-anchor="middle" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="13" font-weight="700" '
        f'fill="{BAR_COLOR}">{html.escape(prefix)}{html.escape(genome)} | '
        f'{html.escape(contig)}</text>\n'
        f'<text x="{PANEL_W}" y="{y_top + 16:.1f}" text-anchor="end" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="9" fill="#444">Panel categories: '
        f'{html.escape(cats)}</text>\n'
        f'<svg x="0" y="{arrows_y:.1f}" width="{PANEL_W}" height="{ARROWS_H}" '
        f'viewBox="0 {ARROWS_TOP} {PANEL_W} {ARROWS_H}">{arrows_content}</svg>\n'
    )


def legend_svg(legend_fragment: str, y_top: float) -> str:
    return (
        f'<svg x="0" y="{y_top:.1f}" width="{PANEL_W}" height="{LEGEND_H}" '
        f'viewBox="0 {ARROWS_BOTTOM} {PANEL_W} {LEGEND_H}">{legend_fragment}</svg>\n'
    )
