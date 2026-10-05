#!/usr/bin/env python3
# Title          : build_atlas_index.py
# Description    : Write a browsable index.html for the architecture gene-map atlases
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/27
# Usage          : python3 build_atlas_index.py

import csv
import html
from pathlib import Path

HERE = Path(__file__).resolve().parent

ARCHS = [
    ("bad_ali_chc", "bad-ali/chc (literature-anchored)"),
    ("catechol_ortho", "Catechol ortho-cleavage (literature-anchored)"),
    ("phenylacetate", "Phenylacetate (literature-anchored)"),
    ("protocatechuate", "Protocatechuate (literature-anchored)"),
    ("beta_ox_core", "CR3-type beta-oxidation/CoA core (literature-anchored)"),
    ("beta_oxidation_full", "Full Beta-oxidation category, same 18-gene panel as Fig. 3a iTOL heatmap (data-driven)"),
    ("cyclohexanecarboxylate_core", "Cyclohexanecarboxylate (data-driven, coherent)"),
    ("catechol_meta_core", "Catechol-meta (data-driven, coherent)"),
    ("aromatics_core", "Aromatics (data-driven, coherent)"),
    ("gallate_core", "Gallate (data-driven, coherent)"),
    ("hydroxyquinol_core", "Hydroxyquinol (data-driven, coherent)"),
    ("gentisate_core", "Gentisate (data-driven, coherent)"),
    ("cyclohexylacetate_core", "Cyclohexylacetate (data-driven, coherent)"),
    ("alkanes_core", "Alkanes (data-driven, thin/contrast)"),
    ("alkenes_core", "Alkenes (data-driven, thin/contrast)"),
    ("benzoyl_coa_core", "Benzoyl-CoA (data-driven, thin/contrast)"),
    ("naphthalene_core", "Naphthalene (data-driven, thin/contrast)"),
    ("oxalate_core", "Oxalate (data-driven, thin/contrast)"),
    ("plastics_core", "Plastics (data-driven, thin/contrast)"),
    ("pyrogallol_core", "Pyrogallol (data-driven, thin/contrast)"),
    ("tannin_core", "Tannin (data-driven, thin/contrast)"),
    ("transportation_core", "Transportation/Mla (data-driven, contextual, non-degradative)"),
    # Full-panel companions: every Gene_info.xlsx gene in the category rather than a 2-5 gene subset.
    ("alkanes_full", "Alkanes, full 38-gene panel (data-driven)"),
    ("alkenes_full", "Alkenes, full 21-gene panel (data-driven)"),
    ("aromatics_full", "Aromatics, full 149-gene panel (data-driven; mixes sub-pathways, see README)"),
    ("aromatics_subcat1_full", "Aromatics subcategory 1 (S8a grouping), 60 genes (data-driven)"),
    ("aromatics_subcat2_full", "Aromatics subcategory 2 (S8b grouping), 46 genes (data-driven)"),
    ("aromatics_subcat3_full", "Aromatics subcategory 3 (S8c grouping), 43 genes (data-driven)"),
    ("benzoate_full", "Benzoate, full 8-gene panel (data-driven)"),
    ("benzoyl_coa_full", "Benzoyl-CoA, full 9-gene panel (data-driven)"),
    ("catechol_meta_full", "Catechol-meta, full 12-gene panel (data-driven)"),
    ("catechol_ortho_full", "Catechol-ortho, full 11-gene panel (data-driven)"),
    ("cyclohexanecarboxylate_full", "Cyclohexanecarboxylate, full 14-gene panel (data-driven)"),
    ("cyclohexylacetate_full", "Cyclohexylacetate, full 3-gene panel (data-driven)"),
    ("gallate_full", "Gallate, full 7-gene panel (data-driven)"),
    ("gentisate_full", "Gentisate, full 4-gene panel (data-driven)"),
    ("naphthalene_full", "Naphthalene, full 13-gene panel (data-driven)"),
    ("oxalate_full", "Oxalate, full 10-gene panel (data-driven)"),
    ("plastics_full", "Plastics, full 99-gene panel (data-driven)"),
    ("protocatechuate_full", "Protocatechuate, full 10-gene panel (data-driven)"),
    ("pyrogallol_full", "Pyrogallol, full 4-gene panel (data-driven)"),
    ("tannin_full", "Tannin, full 7-gene panel (data-driven)"),
    ("transportation_full", "Transportation, full 13-gene panel (data-driven)"),
]

rows = ["<!doctype html><meta charset=\"utf-8\">",
        "<title>NAFC candidate-architecture gene-map atlases</title>",
        "<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#fafafa}"
        "table{border-collapse:collapse;width:100%}th,td{padding:6px 10px;border-bottom:1px solid #ddd;"
        "text-align:left;font-size:13px}th{background:#eee}a{color:#1f4e8c;text-decoration:none}"
        "a:hover{text-decoration:underline}</style>",
        "<h1>NAFC candidate-architecture gene-map atlases</h1>",
        "<p>Full atlas per architecture: every genome/locus qualifying for that signature, pooled "
        "across all three genome-source tiers. Fig. 5, Fig. 6, Fig. S15, and "
        "Fig. S16 each show only a curated representative subset "
        "(<code>Result/Figures/Main_Figures/</code> and <code>Supplementary_Figures/</code>); "
        "this index is the underlying supplementary/data-availability resource. Most categories "
        "have two entries: a narrow hand-picked/literature \"core\" signature (a handful of "
        "specific genes) and a \"full\" companion using every official Gene_info.xlsx panel gene "
        "for that category instead &mdash; these are genuinely different "
        "gene sets, not a relabeling; see README.md. See "
        "<a href=\"README.md\">README.md</a> for methodology, thresholds, and the "
        "literature-anchored vs. data-driven distinction.</p>",
        "<table><tr><th>Architecture</th><th>Category</th><th>Panels</th><th>Genomes</th>"
        "<th>Full atlas</th><th>Manifest</th></tr>"]

for name, label in ARCHS:
    manifest = HERE / f"{name}_manifest.tsv"
    with manifest.open() as f:
        data = list(csv.DictReader(f, delimiter="\t"))
    genomes = len(set(r["Genome"] for r in data))
    rows.append(
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f"<td>{html.escape(label)}</td>"
        f"<td>{len(data)}</td>"
        f"<td>{genomes}</td>"
        f'<td><a href="{name}.svg" target="_blank">{name}.svg</a></td>'
        f'<td><a href="{name}_manifest.tsv">{name}_manifest.tsv</a></td>'
        "</tr>"
    )
rows.append("</table>")
(HERE / "index.html").write_text("\n".join(rows), encoding="utf-8")
print(f"Wrote {HERE / 'index.html'}")
