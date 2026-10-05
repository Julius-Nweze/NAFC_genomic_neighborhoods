#!/usr/bin/env python3
# Title          : build_architecture_gene_map_figures.py
# Description    : Compile candidate-neighborhood SVGs into one gene-map atlas per architecture signature
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/27
# Usage          : python3 build_architecture_gene_map_figures.py

import csv
import html
from pathlib import Path

from panel_style import (
    PANEL_W as STYLE_PANEL_W, COLUMN_HEADER_H, PANEL_TOTAL_H, LEGEND_H as STYLE_LEGEND_H,
    column_header_svg, arrows_and_legend, panel_svg, legend_svg,
)

OPERON_DIR = Path(__file__).resolve().parent
CLUSTERS_TSV = OPERON_DIR / "NA_gene_panel_clusters_summary.tsv"
SVG_DIR = OPERON_DIR / "svgs"
OUT_DIR = OPERON_DIR / "Architecture_gene_maps"

REF_BARE = {"AP042450.1", "CP010516.1", "CP010517.1", "CP063454.1", "CP063455.1", "LMAZ01000003.1"}

ARCHITECTURES = {
    # --- literature-anchored signatures ---
    "bad_ali_chc": (["badH", "badI", "badK", "aliA", "aliB"], 4,
                     "bad-ali/chc architecture (>=4 of badH/badI/badK/aliA/aliB)"),
    "catechol_ortho": (["catA", "catB", "catC", "catR"], 3,
                        "Catechol ortho-cleavage operon (>=3 of catA/catB/catC/catR)"),
    "phenylacetate": (["paaA", "paaB", "paaC", "paaD", "paaE", "paaF", "paaG", "paaH", "paaI",
                        "paaJ", "paaK", "paaX", "paaY", "paaZ"], 8,
                       "Phenylacetate architecture (>=8 of the 14-gene paa panel)"),
    "protocatechuate": (["pcaD", "pcaF", "pcaI", "pcaJ"], 3,
                         "Protocatechuate branch (>=3 of pcaD/pcaF/pcaI/pcaJ)"),
    # hmgL belongs to the Cyclohexylacetate category, so this signature (from the Cupriavidus
    # gilardii CR3 badI-acdA-atoB-hmgL locus) is not equivalent to the full 18-gene
    # Beta-oxidation category.
    "beta_ox_core": (["acdA", "atoB", "badI", "hmgL"], 3,
                      "CR3-type beta-oxidation/CoA core (>=3 of acdA/atoB/badI/hmgL; not the full "
                      "Beta-oxidation category - see README)"),
    # Full 18-gene Beta-oxidation panel. Data-driven threshold: >=3 of the 18 genes co-located
    # in one candidate neighborhood.
    "beta_oxidation_full": (["GMET_RS16560", "acdA", "acdB", "acsA", "aliA", "aliB", "atoB",
                              "badH", "badI", "crt", "echA", "echA1", "fadB", "fadB2", "fadD",
                              "fadE", "fadJ", "fadK"], 3,
                             "Full Beta-oxidation category, same 18-gene panel as the Fig. 3a "
                             "iTOL heatmap (>=3 of 18 co-located; DATA-DRIVEN, not a cited operon)"),
    # Full-panel signatures for the remaining categories (Gene_info.xlsx gene lists), with
    # data-driven thresholds. Not built for Hydroquinone (too rare), Hydroxyquinol (same as
    # hydroxyquinol_core), Beta-oxidation (above) or Phenylacetate (already requires 8 of 14 genes).
    "alkanes_full": (["BVMO", "CHMO", "CYP153", "CYP450", "ahpC", "ahpF", "alkB", "alkH_ald", "alkJ_adh", "alkL", "alkMa", "alkMb", "alkN_mcp", "alkS", "almA", "blc", "bmoX", "bmoY", "bmoZ", "chnA", "chnB", "chnC", "chnD", "chnE", "ladA", "mmoB", "mmoC", "mmoX", "mmoY", "mmoZ", "pmoA1", "prmA", "prmB", "prmC", "prmD", "rdx", "rdxR", "ssuD"], 3,
                'Full Alkanes category (38 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "alkenes_full": (["etnE", "isoA", "isoB", "isoC", "isoD", "isoE", "isoF", "isoH", "isoI", "mpdB", "mpdC", "xamoA", "xamoB", "xamoC", "xamoD", "xamoE", "xamoF", "xecA", "xecC", "xecD", "xecE1"], 2,
                'Full Alkenes category (21 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "aromatics_full": (["GI:23664434", "abmA", "abmG1", "abmG2", "amnA", "amnB", "amnC", "amnD", "andAa", "andAb", "andAc", "andAd", "antA", "antB", "antC", "atdA1", "bedA", "bedB", "bedC1", "bedC2", "benA", "benB", "benC", "benD", "boxA", "boxB", "boxC", "bphA", "bphA1", "bphA2", "bphA3", "bphAa", "bphAb", "bphAc", "bphAd", "bphB", "bphC", "bphD", "bphE", "bphF", "bphG", "bphI", "bphJ", "bphX2", "bphX3", "carAa", "carAc", "carC", "clcA", "cmtAa", "cmtAb", "cmtAc", "cmtAd", "cmtB", "dbfB", "diox1", "diox2", "dmpK", "dmpL", "dmpM", "dmpN", "dmpO", "dmpP", "dpgC", "dszC", "flnB", "hmgA", "hpaB", "hpaC", "hpaH", "hpcB", "hpcC", "hpcD", "hpcE", "hpcG", "hpcH", "mdlD", "nagAa", "nagAb", "nagAc", "nagAd", "nagB", "nagG", "nagH", "nahG", "nahW", "nbaA", "nidA", "nidB", "ntnC", "padAa1", "padAb1", "padAc1", "padAd1", "padB2", "padC2", "pchA", "pchC", "pchF", "phdI", "phdJ", "phdK", "phnB", "pht2", "pht3", "pht4", "pht5", "pobA", "pqsA", "pydA", "pydB", "rauF", "salA", "sdc", "sdgA", "sdgB", "sdgC", "styA", "styB", "tcbAb", "tcbC", "tfdC", "tmoA", "tmoB", "tmoC", "tmoD", "tmoE", "tmoF", "todA", "todB", "todD", "todE", "todF", "tomA1", "tomA2", "tomA3", "tomA4", "tomA5", "tphA1I", "tphA2I", "tphA3I", "tphB", "xylB", "xylC", "xylL", "xylM", "xylX", "xylY", "xylZ"], 5,
                'Full Aromatics category (149 genes; >=5 co-located; DATA-DRIVEN, not a cited '
                'operon; a 149-gene umbrella spanning several unrelated sub-pathways, e.g. '
                'tomA1-tomA5 toluene monooxygenase clusters never co-occur with the separate '
                'hpcB/C/E/G/H homoprotocatechuate clusters that also qualify - see README)'),
    # Aromatics split into the three subcategory groupings used for the iTOL heatmap panels,
    # so each atlas is internally coherent.
    "aromatics_subcat1_full": (["GI:23664434", "abmA", "abmG1", "abmG2", "amnA", "amnB", "amnC", "amnD", "andAa", "andAb", "andAc", "andAd", "antA", "antB", "antC", "atdA1", "bedA", "bedB", "bedC1", "bedC2", "benA", "benB", "benC", "benD", "boxA", "boxB", "boxC", "bphA", "bphA1", "bphA2", "bphA3", "bphAa", "bphAb", "bphAc", "bphAd", "bphB", "bphC", "bphD", "bphE", "bphF", "bphG", "bphI", "bphJ", "bphX2", "bphX3", "hpaB", "hpaC", "hpaH", "nagAa", "nagAb", "nagAc", "nagAd", "nagB", "nbaA", "pqsA", "rauF", "xylL", "xylX", "xylY", "xylZ"], 4,
                'Aromatics subcategory (S8a grouping, 60 genes; >=4 co-located; DATA-DRIVEN, not a cited operon): 1-/2-Methylnaphthalene, 2-Aminophenol, 2-Nitrobenzoate, 4-Hydroxyphenylacetate, Aniline, Anthranilate, Benzene, Biphenyl - same grouping as Fig. S8a. Still spans several sub-pathways: pairwise check found boxB (Benzene/benzoyl-CoA box pathway) does not co-occur with the dominant xylZ/xylL/benA/xylY/benB (benzoate/toluate) cluster - residual heterogeneity within this subcategory, much reduced from flat aromatics_full but not fully eliminated'),
    "aromatics_subcat2_full": (["carAa", "carAc", "carC", "clcA", "dbfB", "dmpK", "dmpL", "dmpM", "dmpN", "dmpO", "dmpP", "dszC", "flnB", "hmgA", "hpcB", "hpcC", "hpcD", "hpcE", "hpcG", "hpcH", "nagG", "nagH", "nahG", "nahW", "nidA", "nidB", "padAa1", "padAb1", "padAc1", "padAd1", "padB2", "padC2", "pht2", "pht3", "pht4", "pht5", "pydA", "pydB", "salA", "sdc", "sdgA", "sdgB", "sdgC", "tcbAb", "tcbC", "tfdC"], 3,
                'Aromatics subcategory (S8b grouping, 46 genes; >=3 co-located; DATA-DRIVEN, not a cited operon): Carbazol, Chlorobenzene, Chlorocatechol, Dibenzofuran, Fluorene, Homogentisate, Homoprotocatechuate, Naphthalene, Phenol, Phthalate, Pyrene, Salicylate - same grouping as Fig. S8b'),
    "aromatics_subcat3_full": (["cmtAa", "cmtAb", "cmtAc", "cmtAd", "cmtB", "diox1", "diox2", "dpgC", "mdlD", "ntnC", "pchA", "pchC", "pchF", "phdI", "phdJ", "phdK", "phnB", "pobA", "styA", "styB", "tmoA", "tmoB", "tmoC", "tmoD", "tmoE", "tmoF", "todA", "todB", "todD", "todE", "todF", "tomA1", "tomA2", "tomA3", "tomA4", "tomA5", "tphA1I", "tphA2I", "tphA3I", "tphB", "xylB", "xylC", "xylM"], 3,
                'Aromatics subcategory (S8c grouping, 43 genes; >=3 co-located; DATA-DRIVEN, not a cited operon): Styrene, Terephthalate, Toluene, p-Cumate - same grouping as Fig. S8c'),

    "benzoate_full": (["BenA-xylX", "benA", "benB", "benB-xylY", "benC", "benC-xylZ", "benD", "benD-xylL"], 3,
                'Full Benzoate category (8 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "benzoyl_coa_full": (["bamB", "bamC", "bcrA", "bcrB", "bcrC", "bcrD", "dch", "had", "oah"], 2,
                'Full Benzoyl-CoA category (9 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "catechol_meta_full": (["catE", "mhpD", "mhpE", "mhpF", "xylE", "xylF", "xylG", "xylH", "xylI", "xylJ", "xylK", "xylQ"], 3,
                'Full Catechol-meta category (12 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "catechol_ortho_full": (["catA", "catB", "catC", "catR", "fadA", "fadI", "pcaD", "pcaF", "pcaI", "pcaJ", "pcaL"], 4,
                'Full Catechol-ortho category (11 genes; >=4 co-located; DATA-DRIVEN, not a cited operon)'),
    "cyclohexanecarboxylate_full": (["aabB", "aabC", "aabD", "badK", "chcA", "chcAa", "chcAb", "chcAc", "chcB1", "chcB2", "chcC1", "chcC2", "pimC/aabA", "pimD"], 3,
                'Full Cyclohexanecarboxylate category (14 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "cyclohexylacetate_full": (["atoA", "atoD", "hmgL"], 2,
                'Full Cyclohexylacetate category (3 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "gallate_full": (["galA", "galB", "galC", "galD", "galP", "galR", "galT"], 3,
                'Full Gallate category (7 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "gentisate_full": (["mhbD", "mhbH", "mhbl", "mobA"], 3,
                'Full Gentisate category (4 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "naphthalene_full": (["nagAa", "nagAb", "nagAc", "nagAd", "nagB", "nagC", "nagD", "nagE", "nagF", "nagG", "nagH", "nagR", "nahG"], 3,
                'Full Naphthalene category (13 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),
    "oxalate_full": (["AAE3", "frc", "oorA", "oorB", "oorD", "oxc", "oxdD", "oxlT", "uctC", "yfdE"], 2,
                'Full Oxalate category (10 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "plastics_full": (["AAW51_2473", "ABO_1197", "ABO_2449", "AFUA_4G03560", "ALC24_4107", "AvCA6_03910", "BAY15_3292", "BSQ33_03270", "CLE1", "CLOSTHATH_07193", "DLH", "FB4_0780", "GEN0105", "HRbin29_00073", "ISF6_0224", "ISF6_4831", "NFA_28320", "OA86_11720", "PLAase2", "PLAase3", "PROK", "PaCLE1", "RBTH_03267", "RPA1511", "Rru_A1969", "SAMN04487908_10753", "SUBB_LEDLE", "SUBS_LEDLE", "Tcur_0390", "Tcur_1278", "WG66_16246", "amdA", "bta2", "chi25", "cut", "cut1", "cut190", "cut2", "cut3", "cutA", "cutL", "cut_1", "est", "est1", "est2", "estA", "estC", "fkbU", "hbd", "lac", "latA", "lcp", "lip", "lip1", "lipA", "lipB", "lipIAF5-2", "mnp1", "nylA", "nylB", "nylC", "oph", "pDS-PS", "pbath", "pbsA", "pclh", "pegA", "pegC", "peth", "phaY2", "phaZ", "phaZ1", "phaZ2", "phaZ3", "phaZ4", "phaZ5", "phaZ6", "phaZ7", "phaZCac", "phaZPfu", "phaZbd", "phaZc", "phaZd", "phaZpst", "plaA", "plaM4", "plaM5", "plaM7", "plaM9", "pnbA", "pudA", "pueA", "pueB", "pulA", "pvaA", "pvaB", "roxA", "roxB", "subC"], 2,
                'Full Plastics category (99 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "protocatechuate_full": (["pcaB", "pcaC", "pcaD", "pcaF", "pcaG", "pcaH", "pcaI", "pcaJ", "pobA", "pobB"], 5,
                'Full Protocatechuate category (10 genes; >=5 co-located; DATA-DRIVEN, not a cited operon)'),
    "pyrogallol_full": (["DPGC", "PGRcl", "athL", "bthL"], 2,
                'Full Pyrogallol category (4 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "tannin_full": (["lpdB", "lpdC", "lpdD", "tanA", "tanB", "tanL", "ubiX"], 2,
                'Full Tannin category (7 genes; >=2 co-located; DATA-DRIVEN, not a cited operon)'),
    "transportation_full": (["fadL", "mlaB", "mlaD", "mlaE", "mlaF", "ompW", "opdK", "pcaK", "pqiA", "pqiB", "pqiC", "yebS", "yebT"], 3,
                'Full Transportation category (13 genes; >=3 co-located; DATA-DRIVEN, not a cited operon)'),

    # --- data-driven signatures: top co-occurring genes within each category's multi-hit
    # neighborhoods (>=2 panel genes). Default threshold >=2 of the top 4 genes, adjusted per
    # category. Hydroquinone is excluded (its 2 genes co-occur in only 1 genome).
    "alkanes_core": (["alkJ_adh", "alkN_mcp", "ahpF", "BVMO"], 2,
                      "Alkanes core (>=2 of alkJ_adh/alkN_mcp/ahpF/BVMO)"),
    "alkenes_core": (["mpdC", "xecD", "xecE1", "mpdB"], 2,
                      "Alkenes core (>=2 of mpdC/xecD/xecE1/mpdB)"),
    "aromatics_core": (["xylB", "GI:23664434", "hpcH", "pchA"], 2,
                        "Aromatics core (>=2 of xylB/GI:23664434/hpcH/pchA)"),
    "benzoyl_coa_core": (["bamB", "bamC"], 2,
                          "Benzoyl-CoA core (bamB+bamC pair; data-driven, bcrD and oah excluded "
                          "because they rarely co-occur with bamB/bamC)"),
    "catechol_meta_core": (["mhpD", "xylG", "mhpF", "xylI"], 2,
                            "Catechol-meta core (>=2 of mhpD/xylG/mhpF/xylI)"),
    "cyclohexanecarboxylate_core": (["pimC/aabA", "pimD", "chcB2", "chcB1"], 3,
                                     "Cyclohexanecarboxylate core (>=3 of pimC/aabA/pimD/chcB2/chcB1)"),
    "cyclohexylacetate_core": (["atoA", "hmgL", "atoD"], 2,
                                "Cyclohexylacetate core (>=2 of the 3-gene atoA/hmgL/atoD panel)"),
    "gallate_core": (["galT", "galC", "galD", "galP"], 2,
                      "Gallate core (>=2 of galT/galC/galD/galP)"),
    "gentisate_core": (["mhbl", "mhbH", "mhbD"], 2,
                        "Gentisate core (>=2 of mhbl/mhbH/mhbD; data-driven, mobA excluded "
                        "because it never co-occurs with these three)"),
    "hydroxyquinol_core": (["chqB", "macA"], 2,
                            "Hydroxyquinol (both genes in the 2-gene chqB/macA panel)"),
    "naphthalene_core": (["nagR", "nagF", "nagAa", "nagG"], 2,
                          "Naphthalene core (>=2 of nagR/nagF/nagAa/nagG)"),
    "oxalate_core": (["uctC", "yfdE", "AAE3", "oxlT"], 2,
                      "Oxalate core (>=2 of uctC/yfdE/AAE3/oxlT)"),
    "plastics_core": (["est", "hbd", "pnbA", "pegC"], 2,
                       "Plastics core (>=2 of est/hbd/pnbA/pegC)"),
    "pyrogallol_core": (["PGRcl", "bthL", "athL", "DPGC"], 2,
                         "Pyrogallol core (>=2 of PGRcl/bthL/athL/DPGC)"),
    "tannin_core": (["lpdB", "ubiX", "lpdC", "tanL"], 2,
                     "Tannin core (>=2 of lpdB/ubiX/lpdC/tanL)"),
    "transportation_core": (["mlaF", "mlaE", "mlaD"], 3,
                             "Mla lipid-transport complex (all 3 of mlaF/mlaE/mlaD; data-driven, "
                             "opdK excluded because it never co-occurs with these three; broad "
                             "transport machinery, not a degradative architecture)"),
}

PANEL_W = STYLE_PANEL_W
GAP = 10
# Per-panel header (Source / Genome name / Category) from panel_style, with one shared
# column-title row and one shared legend per composite.
HEADER_H = COLUMN_HEADER_H + 20


def tier(genome: str) -> str:
    if genome.startswith("Exp3_MAG_"):
        return "MAG"
    if genome.startswith("GCF_") or genome in REF_BARE:
        return "Ref"
    return "Colla"


def load_rows():
    rows = []
    with open(CLUSTERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["genes"] = set(g.strip() for g in row["Panel genes"].split(";") if g.strip())
            row["seedn"] = int(row["Panel seed hit count"])
            rows.append(row)
    return rows


def inner_svg_content(path: Path) -> str:
    text = path.read_text()
    start = text.index(">", text.index("<svg")) + 1
    end = text.rindex("</svg>")
    return text[start:end]


def build_figure(name: str, geneset: list[str], minhit: int, caption: str, rows: list[dict]) -> None:
    members = []
    for row in rows:
        n = len(row["genes"] & set(geneset))
        if n >= minhit:
            members.append((row, n))
    members.sort(key=lambda x: (tier(x[0]["Genome"]), x[0]["Genome"], -x[1]))

    from render_missing_cluster_svgs import ensure_rendered
    failed = ensure_rendered({row["Cluster ID"] for row, n in members})
    if failed:
        print(f"  WARNING: {len(failed)} clusters could not be rendered (not found in source "
              f"tables): {failed}")

    OUT_DIR.mkdir(exist_ok=True)
    manifest_path = OUT_DIR / f"{name}_manifest.tsv"
    with open(manifest_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Cluster ID", "Genome", "Tier", "Genes matched", "Of", "Seed hit count"])
        for row, n in members:
            w.writerow([row["Cluster ID"], row["Genome"], tier(row["Genome"]), n, len(geneset),
                        row["seedn"]])

    panel_y = HEADER_H
    body_parts = [column_header_svg(panel_y - 6)]
    panel_y += 4
    missing = []
    legend_fragment = None
    for row, n in members:
        cluster_id = row["Cluster ID"]
        genome = row["Genome"]
        gtier = tier(genome)
        src = SVG_DIR / f"{cluster_id}.svg"
        if not src.exists():
            missing.append(cluster_id)
            continue
        arrows_content, legend_rest = arrows_and_legend(inner_svg_content(src))
        if legend_fragment is None:
            legend_fragment = legend_rest
        body_parts.append(panel_svg(row, gtier, panel_y, arrows_content))
        panel_y += PANEL_TOTAL_H + GAP

    if legend_fragment is not None:
        body_parts.append(legend_svg(legend_fragment, panel_y))
        panel_y += STYLE_LEGEND_H + GAP

    total_h = panel_y + 10
    total_w = PANEL_W + 40
    header = (
        f'<rect width="100%" height="100%" fill="#f3f3f0"/>\n'
        f'<text x="20" y="26" font-family="Arial, Helvetica, sans-serif" font-size="20" '
        f'font-weight="700" fill="#111">{html.escape(caption)}</text>\n'
    )
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" '
           f'viewBox="0 0 {total_w} {total_h}">\n' + header + "".join(body_parts) + "</svg>\n")
    (OUT_DIR / f"{name}.svg").write_text(svg)
    print(f"{name}: {len(members)} panels ({len(set(r['Genome'] for r, n in members))} genomes), "
          f"missing SVGs: {len(missing)} -> {OUT_DIR / (name + '.svg')}")
    if missing:
        print("  MISSING:", missing)


def main():
    rows = load_rows()
    for name, (geneset, minhit, caption) in ARCHITECTURES.items():
        build_figure(name, geneset, minhit, caption, rows)


if __name__ == "__main__":
    main()
