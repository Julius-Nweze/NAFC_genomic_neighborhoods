#!/usr/bin/env python3
# Title          : build_annotation_confirmation_summary.py
# Description    : Chart annotation-confirmation outcomes by genome-source tier
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/01
# Usage          : python3 build_annotation_confirmation_summary.py [options]

"""Chart annotation-confirmation outcomes by genome-source tier."""
import argparse
from pathlib import Path

ACTIONS = [
    ("no_hit", "Swiss-Prot no hit", "#c7c7c7"),
    ("review", "Review low-confidence hit", "#fdae61"),
    ("accept_name", "Accept name update", "#1b9e77"),
    ("accept_product", "Accept product only", "#66c2a5"),
    ("do_not_rename", "Do not rename weak/uncharacterized", "#7570b3"),
]

DATASETS = [
    {
        "label": "MAGs",
        "counts": {"no_hit": 2785, "review": 1629, "accept_name": 38, "accept_product": 480, "do_not_rename": 745},
    },
    {
        "label": "NCBI references",
        "counts": {"no_hit": 1282, "review": 882, "accept_name": 37, "accept_product": 335, "do_not_rename": 483},
    },
    {
        "label": "OSPW-containing-mesocosm isolates",
        "counts": {"no_hit": 483, "review": 280, "accept_name": 10, "accept_product": 141, "do_not_rename": 232},
    },
]

MAX_BAR_PX = 650.0
BAR_X = 220.0
BAR_HEIGHT = 62
ROW_GAP = 147
FIRST_ROW_Y = 100
LEFT_LABEL_X = 40
CANVAS_WIDTH = 1100


def esc(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def svg_text(x, y, text, size, weight=400, anchor="start", fill="#222"):
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{fill}">{esc(text)}</text>'
    )


def build_svg(datasets, max_total):
    scale = MAX_BAR_PX / max_total
    elements = []
    elements.append(svg_text(40, 45, "Focused annotation confirmation of caution targets", 28, 700))
    elements.append(
        svg_text(
            40,
            72,
            "No-hit proteins remain candidates for broader TrEMBL, UniRef, domain, or HMM confirmation.",
            14,
            400,
            fill="#555",
        )
    )

    row_y = FIRST_ROW_Y
    for ds in datasets:
        total = sum(ds["counts"].values())
        label_y = row_y + BAR_HEIGHT - 24
        label = ds["label"]
        if len(label) > 20:
            # split near the midpoint, preferring a space or hyphen so neither line overflows
            breakpoints = [i for i, c in enumerate(label) if c in " -"]
            mid = min(breakpoints, key=lambda i: abs(i - len(label) / 2))
            line1, line2 = label[: mid + (1 if label[mid] == "-" else 0)], label[mid + 1 :]
            elements.append(svg_text(LEFT_LABEL_X, label_y - 10, line1.strip(), 13, 700))
            elements.append(svg_text(LEFT_LABEL_X, label_y + 9, line2.strip(), 13, 700))
        else:
            elements.append(svg_text(LEFT_LABEL_X, label_y, label, 18, 700))

        x = BAR_X
        no_hit_count = ds["counts"]["no_hit"]
        no_hit_w = no_hit_count * scale
        elements.append(
            f'<rect x="{x:.1f}" y="{row_y}" width="{no_hit_w:.1f}" height="{BAR_HEIGHT}" '
            f'fill="#c7c7c7" stroke="#fff" stroke-width="1"/>'
        )
        elements.append(
            svg_text(x + no_hit_w / 2, label_y, str(no_hit_count), 14, 700, anchor="middle")
        )
        x += no_hit_w
        for key, _, color in ACTIONS[1:]:
            w = ds["counts"][key] * scale
            if w > 0:
                elements.append(
                    f'<rect x="{x:.1f}" y="{row_y}" width="{w:.1f}" height="{BAR_HEIGHT}" '
                    f'fill="{color}" stroke="#fff" stroke-width="1"/>'
                )
            x += w

        elements.append(svg_text(895, label_y, f"{total:,} targets", 16, 700))
        row_y += ROW_GAP

    legend_y = row_y - ROW_GAP + BAR_HEIGHT + 28
    elements.append(svg_text(LEFT_LABEL_X, legend_y, "Suggested action", 16, 700))
    swatch_row1_y = legend_y + BAR_HEIGHT - 48
    swatch_row2_y = swatch_row1_y + 28
    positions = [
        (LEFT_LABEL_X, swatch_row1_y),
        (270, swatch_row1_y),
        (560, swatch_row1_y),
        (LEFT_LABEL_X, swatch_row2_y),
        (280, swatch_row2_y),
    ]
    for (key, label, color), (sx, sy) in zip(ACTIONS, positions):
        elements.append(
            f'<rect x="{sx}" y="{sy}" width="18" height="18" fill="{color}" '
            f'stroke="#777" stroke-width="0.5"/>'
        )
        elements.append(svg_text(sx + 26, sy + 14, label, 13, 400))

    footnote_y = swatch_row2_y + 40
    elements.append(
        svg_text(
            LEFT_LABEL_X,
            footnote_y,
            "Swiss-Prot is a conservative curated database; sparse hits are expected for novel environmental genomes.",
            13,
            400,
            fill="#555",
        )
    )

    height = footnote_y + 20
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS_WIDTH}" height="{height}" '
        f'viewBox="0 0 {CANVAS_WIDTH} {height}">\n'
        f'<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(elements)
        + "\n</svg>\n"
    )
    return svg


def write_data_tsv(path, datasets):
    label_lookup = {key: label for key, label, _ in ACTIONS}
    with open(path, "w") as handle:
        handle.write("dataset\taction\tcount\n")
        for ds in datasets:
            for key, _, _ in ACTIONS:
                handle.write(f"{ds['label']}\t{label_lookup[key]}\t{ds['counts'][key]}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        default="/path/to/your/Figures",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    max_total = max(sum(ds["counts"].values()) for ds in DATASETS)
    svg = build_svg(DATASETS, max_total)
    (out_dir / "annotation_confirmation_summary.svg").write_text(svg)
    write_data_tsv(out_dir / "annotation_confirmation_summary_data.tsv", DATASETS)
    print(f"Wrote {out_dir / 'annotation_confirmation_summary.svg'}")
    print(f"Wrote {out_dir / 'annotation_confirmation_summary_data.tsv'}")


if __name__ == "__main__":
    main()
