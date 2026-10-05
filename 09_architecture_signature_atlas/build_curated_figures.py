#!/usr/bin/env python3
# Title          : build_curated_figures.py
# Description    : Build curated gene-map figure panels from the architecture atlases
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/25
# Usage          : python3 build_curated_figures.py

import csv
import html
from pathlib import Path

from panel_style import (
    PANEL_W, COLUMN_HEADER_H, PANEL_TOTAL_H, LEGEND_H, column_header_svg, arrows_and_legend,
    panel_svg as render_panel, legend_svg,
)

OPERON_DIR = Path(__file__).resolve().parent
CLUSTERS_TSV = OPERON_DIR / "NAFC_gene_panel_neighborhood_summary.tsv"
SVG_DIR = OPERON_DIR / "svgs"
BAR_CHART_SVG = OPERON_DIR / "figure_category_multihit_summary.svg"
MAIN_FIG_DIR = OPERON_DIR.parent / "Figures" / "Main_Figures"
SUPP_FIG_DIR = OPERON_DIR.parent / "Figures" / "Supplementary_Figures"

REF_BARE = {"AP042450.1", "CP010516.1", "CP010517.1", "CP063454.1", "CP063455.1", "LMAZ01000003.1"}
GAP = 10
COLUMN_HEADER_Y = 72   # below the composite's own title (y=26) and subtitle (y=46)
HEADER_H = COLUMN_HEADER_Y + 10


def tier(genome: str) -> str:
    if genome.startswith("Exp3_MAG_"):
        return "MAG"
    if genome.startswith("GCF_") or genome in REF_BARE:
        return "Ref"
    return "Colla"


def load_cluster_rows():
    with open(CLUSTERS_TSV) as f:
        return {row["Cluster ID"]: row for row in csv.DictReader(f, delimiter="\t")}


CLUSTERS = load_cluster_rows()

# (cluster_id, one-line reason this panel was picked)
FIG5_PANELS = [
    ("GCF_029269295.1_panelcluster_16141",
     "intact catA/catB/catC/catR operon, all 4 genes seed hits, already cited in Results"),
    ("Exp3_MAG_68_panelcluster_12479",
     "catA/catB/catC + benA/benB/benD, plant-root MAG, already cited in Results"),
    ("GROW_115_consensus_panelcluster_16861",
     "catA/catB/catC/catR + benzoate genes at one locus, OSPW isolate"),
]
FIG6_PANELS = [
    ("GCF_000470885.1_panelcluster_14061",
     "bad-ali/chc, Rhodococcus aetherivorans BCP1 (literature-characterized, 5/5 genes)"),
    ("GCF_001190925.1_panelcluster_14753",
     "bad-ali/chc, Aromatoleum/Azoarcus sp. CIB (literature-characterized, 5/5 genes)"),
    ("Exp3_MAG_162_panelcluster_3892",
     "bad-ali/chc, plant-root MAG Curvibacter sp. (5/5 genes)"),
    ("GROW_115_consensus_panelcluster_16852",
     "bad-ali/chc, OSPW isolate Acinetobacter bohemicus (4/5 genes)"),
    ("GCF_000736435.1_panelcluster_14676",
     "phenylacetate, Rhodococcus opacus R7 (literature-characterized; not a bad-ali/chc match, see caption)"),
    ("CP010516.1_panelcluster_177",
     "CR3-type beta-oxidation/CoA core, Cupriavidus gilardii CR3 (the literature source locus)"),
]
S15_PANELS = [
    ("cyclohexanecarboxylate_core", "GROW_123_consensus_panelcluster_16989"),
    ("catechol_meta_core", "GCF_040170545.1_panelcluster_16253"),
    ("aromatics_core", "GCF_000470885.1_panelcluster_14061"),
    ("gallate_core", "GROW_68_consensus_panelcluster_17408"),
    ("hydroxyquinol_core", "GCF_000736435.1_panelcluster_14694"),
    ("gentisate_core", "Exp3_MAG_129_panelcluster_1988"),
    ("cyclohexylacetate_core", "GCF_000012925.1_panelcluster_13900"),
]
S16_PANELS = [
    ("alkanes_core", "GCF_026428335.1_panelcluster_15500"),
    ("alkenes_core", "GCF_000736435.1_panelcluster_14564"),
    ("benzoyl_coa_core", "GCF_002009335.2_panelcluster_14829"),
    ("naphthalene_core", "Exp3_MAG_129_panelcluster_1988"),
    ("oxalate_core", "GCF_000736435.1_panelcluster_14670"),
    ("plastics_core", "GCF_029269195.1_panelcluster_15761"),
    ("pyrogallol_core", "Exp3_MAG_210_panelcluster_5887"),
    ("tannin_core", "GROW_73_consensus_panelcluster_17596"),
    ("transportation_core", "Exp3_MAG_56_panelcluster_11962"),
]


def inner_svg_content(path: Path) -> str:
    text = path.read_text()
    start = text.index(">", text.index("<svg")) + 1
    end = text.rindex("</svg>")
    return text[start:end]


def letter(i: int) -> str:
    return chr(ord("A") + i)


def build_composite(out_dir: Path, title: str, subtitle: str, panels: list[tuple], category_label=False,
                     lead_bar_chart=False) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    body = [column_header_svg(COLUMN_HEADER_Y)]
    manifest_rows = []
    y = HEADER_H
    total_w = PANEL_W + 40
    idx = 0

    if lead_bar_chart:
        bar_svg_raw = BAR_CHART_SVG.read_text()
        bw = int(bar_svg_raw.split('width="')[1].split('"')[0])
        bh = int(bar_svg_raw.split('height="')[1].split('"')[0])
        scale = PANEL_W / bw
        disp_h = bh * scale
        body.append(f'<text x="20" y="{y + 14}" font-family="Arial, Helvetica, sans-serif" '
                    f'font-size="16" font-weight="700" fill="#111">A</text>')
        body.append(f'<svg x="50" y="{y}" width="{PANEL_W - 30}" height="{disp_h:.0f}" '
                    f'viewBox="0 0 {bw} {bh}">{inner_svg_content(BAR_CHART_SVG)}</svg>')
        y += disp_h + GAP + 10
        manifest_rows.append(("A", "figure_category_multihit_summary.svg", "-", "-", "-",
                               "percent of panel seed hits per category found in a multi-hit cluster, all 22 categories"))
        idx = 1

    legend_fragment = None
    for entry in panels:
        if category_label:
            cat_name, cluster_id = entry
            reason = f"top seed-hit-count example of the {cat_name} signature"
        else:
            cluster_id, reason = entry
            cat_name = ""
        row = CLUSTERS[cluster_id]
        genome = row["Genome"]
        gtier = tier(genome)
        panel_letter = letter(idx)
        prefix = f"{panel_letter}. {cat_name}: " if cat_name else f"{panel_letter}. "
        arrows_content, legend_rest = arrows_and_legend(
            inner_svg_content(SVG_DIR / f"{cluster_id}.svg"))
        if legend_fragment is None:
            legend_fragment = legend_rest
        body.append(render_panel(row, gtier, y, arrows_content, prefix=prefix))
        manifest_rows.append((panel_letter, cluster_id, genome, gtier, row["Panel seed hit count"], reason))
        y += PANEL_TOTAL_H + GAP
        idx += 1

    if legend_fragment is not None:
        body.append(legend_svg(legend_fragment, y))
        y += LEGEND_H + GAP

    total_h = y + 10
    header = (
        f'<rect width="100%" height="100%" fill="#f3f3f0"/>\n'
        f'<text x="20" y="26" font-family="Arial, Helvetica, sans-serif" font-size="20" '
        f'font-weight="700" fill="#111">{html.escape(title)}</text>\n'
        f'<text x="20" y="46" font-family="Arial, Helvetica, sans-serif" font-size="12" '
        f'fill="#555">{html.escape(subtitle)}</text>\n'
    )
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" '
           f'viewBox="0 0 {total_w} {total_h}">\n' + header + "".join(body) + "</svg>\n")
    (out_dir / "composite.svg").write_text(svg)

    with open(out_dir / "manifest.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Panel", "Cluster ID / source", "Genome", "Tier", "Seed hits", "Why this panel"])
        for r in manifest_rows:
            w.writerow(r)
    print(f"Wrote {out_dir / 'composite.svg'} ({idx} gene-map panels{' + bar chart' if lead_bar_chart else ''})")


def main():
    build_composite(
        MAIN_FIG_DIR / "Figure_5_aromatic_ring_cleavage_maps",
        "Fig. 5 [pending number]. Catechol ortho-cleavage and protocatechuate gene maps",
        "Panel A: all-22-category multi-hit-clustering summary. Panels B-D: representative loci, one per genome-source tier.",
        FIG5_PANELS, lead_bar_chart=True,
    )
    build_composite(
        MAIN_FIG_DIR / "Figure_6_literature_linked_maps",
        "Fig. 6 [pending number]. Literature-linked and NAFC-relevant gene maps",
        "bad-ali/chc (A-D, all tiers), phenylacetate (E), and CR3-type beta-oxidation/CoA core (F).",
        FIG6_PANELS,
    )
    build_composite(
        SUPP_FIG_DIR / "Figure_S15_coherent_data_driven_maps",
        "Fig. S15 [pending number]. Coherent data-driven recurring neighborhoods",
        "One representative example (highest seed-hit count) per category. Not literature-cited operons.",
        S15_PANELS, category_label=True,
    )
    build_composite(
        SUPP_FIG_DIR / "Figure_S16_thin_contrast_maps",
        "Fig. S16 [pending number]. Limited-co-occurrence and non-degradative contrast panels",
        "One representative example per category, shown as a contrast to Fig. 5/6/S15, not as architectures.",
        S16_PANELS, category_label=True,
    )


if __name__ == "__main__":
    main()
