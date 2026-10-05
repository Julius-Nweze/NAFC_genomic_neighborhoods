#!/usr/bin/env python3
# Title          : build_row_ospw_mag_gene_content_comparison.py
# Description    : Compare candidate gene and neighborhood content of MAGs detected under ROW versus OSPW
# Author         : Julius Eyiuche Nweze
# Date           : 2026/07/31
# Usage          : python3 build_row_ospw_mag_gene_content_comparison.py [options]

"""Compare candidate gene and neighborhood content of MAGs detected under ROW versus OSPW."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_ABUNDANCE = "/path/to/your/MAGs_abundance_metadata_MG.csv"
DEFAULT_BESTHIT = "/path/to/your/Ref_genomes/Analysis/Result/ITOL/gene_count_summary/besthit_assigned_hit_proteins.tsv"
DEFAULT_CLUSTERS = "/path/to/your/Ref_genomes/Analysis/Result/Operon_clusters/NAFC_gene_panel_neighborhood_summary.tsv"


def classify_group(row):
    if row["detected_ROW_only"]:
        return "ROW_only"
    if row["detected_OSPW_only"]:
        return "OSPW_only"
    if row["detected_both"]:
        return "both"
    return "neither"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abundance", default=DEFAULT_ABUNDANCE)
    parser.add_argument("--besthit", default=DEFAULT_BESTHIT)
    parser.add_argument("--clusters", default=DEFAULT_CLUSTERS)
    parser.add_argument(
        "--out-dir",
        default="/path/to/your/Ref_genomes/Analysis/Result/ROW_vs_OSPW_MAG_comparison",
    )
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_all = pd.read_csv(args.abundance)
    # Exp3_MAG_2 is the non-microbial, arthropod-derived (Chironomus riparius) MAG excluded from
    # every other panel-search/comparative-tree analysis in this pipeline (see NCBI_genome_info.xlsx,
    # Phylum=Arthropoda). It was never searched against the reference protein panel, so it always has
    # zero panel hits by construction, not by biology, and must be excluded here too for consistency.
    df_all = df_all[df_all["MAGs"] != "Exp3_MAG_2"].copy()
    df = df_all[df_all["Time"] != "D-0"].copy()

    besthit = pd.read_csv(args.besthit, sep="\t")
    mag_hits = besthit[besthit["genome_id"].str.startswith("Exp3_MAG_")]
    panel_hits = mag_hits.groupby("genome_id").size().rename("panel_hits")
    panel_categories = mag_hits.groupby("genome_id")["chosen_category"].nunique().rename("panel_categories")

    clusters = pd.read_csv(args.clusters, sep="\t")
    # Multi-hit = >=2 independently found panel seed hits in the cluster span; three-plus-hit = >=3.
    clusters["is_relevant"] = clusters["Panel seed hit count"] >= 2
    clusters["is_high"] = clusters["Panel seed hit count"] >= 3
    cluster_summary = (
        clusters.groupby("Genome")
        .agg(total_clusters=("Cluster ID", "count"), relevant_clusters=("is_relevant", "sum"), high_priority_clusters=("is_high", "sum"))
        .reset_index()
        .rename(columns={"Genome": "MAG"})
    )

    # ---- Analysis 1+2: per-MAG detection groups ----
    g = (
        df.groupby(["MAGs", "Water type"])["Relative_abundance"]
        .agg(mean_abund="mean", detection_rate=lambda x: (x > 0).mean())
        .reset_index()
    )
    piv_mean = g.pivot(index="MAGs", columns="Water type", values="mean_abund").fillna(0)
    piv_det = g.pivot(index="MAGs", columns="Water type", values="detection_rate").fillna(0)
    piv_mean.columns = [f"mean_{c}" for c in piv_mean.columns]
    piv_det.columns = [f"det_{c}" for c in piv_det.columns]
    wt = piv_mean.join(piv_det)
    wt["detected_ROW_only"] = (wt["det_ROW"] > 0) & (wt["det_OSPW"] == 0)
    wt["detected_OSPW_only"] = (wt["det_OSPW"] > 0) & (wt["det_ROW"] == 0)
    wt["detected_both"] = (wt["det_ROW"] > 0) & (wt["det_OSPW"] > 0)
    wt["group"] = wt.apply(classify_group, axis=1)
    wt = wt.reset_index().rename(columns={"MAGs": "MAG"})

    merged = (
        wt.merge(panel_hits, left_on="MAG", right_index=True, how="left")
        .merge(panel_categories, left_on="MAG", right_index=True, how="left")
        .merge(cluster_summary, on="MAG", how="left")
    )
    for col in ["panel_hits", "panel_categories", "total_clusters", "relevant_clusters", "high_priority_clusters"]:
        merged[col] = merged[col].fillna(0)
    merged.to_csv(out_dir / "MAG_ROW_vs_OSPW_gene_content_per_MAG.tsv", sep="\t", index=False)

    row_only = merged[merged["group"] == "ROW_only"]
    ospw_only = merged[merged["group"] == "OSPW_only"]
    # Restrict "both" to gene-carrying MAGs (panel_hits > 0), matching the Analysis-3 gene-repertoire
    # universe below. ROW_only/OSPW_only are already 100% gene-carrying by construction (see
    # pct_with_any_panel_hit), but "both" (raw water-type detection) includes one MAG with zero panel
    # hits; leaving it in would silently drag a structurally-zero data point into every gene/cluster
    # statistic computed over this group (Mann-Whitney vs both, both Spearman correlations, Fig. S14A-B).
    both = merged[(merged["group"] == "both") & (merged["panel_hits"] > 0)].copy()
    both["log2fc_OSPW_over_ROW"] = np.log2((both["mean_OSPW"] + 1e-6) / (both["mean_ROW"] + 1e-6))

    summary_rows = []
    for name, grp in [("ROW_only", row_only), ("OSPW_only", ospw_only), ("both", both)]:
        summary_rows.append(
            {
                "group": name,
                "n_MAGs": len(grp),
                "pct_with_any_panel_hit": round((grp["panel_hits"] > 0).mean() * 100, 1),
                "mean_panel_hits": round(grp["panel_hits"].mean(), 2),
                "median_panel_hits": grp["panel_hits"].median(),
                "mean_panel_categories_of_22": round(grp["panel_categories"].mean(), 2),
                "mean_relevant_clusters": round(grp["relevant_clusters"].mean(), 2),
                "mean_high_priority_clusters": round(grp["high_priority_clusters"].mean(), 2),
                "pct_with_high_priority_cluster": round((grp["high_priority_clusters"] > 0).mean() * 100, 1),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "MAG_ROW_vs_OSPW_gene_content_summary.tsv", sep="\t", index=False)

    u_hits_row_ospw, p_hits_row_ospw = stats.mannwhitneyu(row_only["panel_hits"], ospw_only["panel_hits"], alternative="two-sided")
    u_clust_row_ospw, p_clust_row_ospw = stats.mannwhitneyu(row_only["relevant_clusters"], ospw_only["relevant_clusters"], alternative="two-sided")
    u_hits_row_both, p_hits_row_both = stats.mannwhitneyu(row_only["panel_hits"], both["panel_hits"], alternative="two-sided")
    u_clust_row_both, p_clust_row_both = stats.mannwhitneyu(row_only["relevant_clusters"], both["relevant_clusters"], alternative="two-sided")
    r_hits, p_r_hits = stats.spearmanr(both["log2fc_OSPW_over_ROW"], both["panel_hits"])
    r_clust, p_r_clust = stats.spearmanr(both["log2fc_OSPW_over_ROW"], both["relevant_clusters"])

    # ---- Analysis 3: gene-level repertoire overlap ----
    det_by_mag_water = df.groupby(["MAGs", "Water type"])["Relative_abundance"].apply(lambda x: (x > 0).any())
    det_piv = det_by_mag_water.reset_index().pivot(index="MAGs", columns="Water type", values="Relative_abundance").fillna(False)
    mags_with_hits = set(mag_hits["genome_id"].unique())
    row_gene_carrying = set(det_piv[det_piv.get("ROW", False) == True].index) & mags_with_hits
    ospw_gene_carrying = set(det_piv[det_piv.get("OSPW", False) == True].index) & mags_with_hits
    row_genes = set(mag_hits[mag_hits["genome_id"].isin(row_gene_carrying)]["chosen_gene"].unique())
    ospw_genes = set(mag_hits[mag_hits["genome_id"].isin(ospw_gene_carrying)]["chosen_gene"].unique())

    repertoire_lines = [
        "# Gene-level repertoire overlap (post-exposure samples, D6-D89)",
        "",
        f"Gene-carrying MAGs detected under ROW: {len(row_gene_carrying)}",
        f"Gene-carrying MAGs detected under OSPW: {len(ospw_gene_carrying)}",
        f"Shared gene-carrying MAGs: {len(row_gene_carrying & ospw_gene_carrying)}",
        f"Distinct reference-panel genes recovered under ROW: {len(row_genes)}",
        f"Distinct reference-panel genes recovered under OSPW: {len(ospw_genes)}",
        f"Shared distinct genes: {len(row_genes & ospw_genes)}",
        f"OSPW-only genes: {sorted(ospw_genes - row_genes)}",
        f"ROW-only genes: {sorted(row_genes - ospw_genes)}",
    ]
    (out_dir / "MAG_ROW_vs_OSPW_gene_repertoire_overlap.txt").write_text("\n".join(repertoire_lines) + "\n")

    # ---- Analysis 4: abundance-weighted community gene-copy burden ----
    df2 = df.merge(panel_hits.rename("n_hits"), left_on="MAGs", right_index=True, how="left")
    df2["n_hits"] = df2["n_hits"].fillna(0)
    df2["weighted_burden"] = df2["Relative_abundance"] * df2["n_hits"]
    sample_burden = df2.groupby(["SampleID", "Water type", "Plant species"])["weighted_burden"].sum().reset_index()
    burden_by_water = sample_burden.groupby("Water type")["weighted_burden"].mean()
    burden_by_plant_water = sample_burden.groupby(["Plant species", "Water type"])["weighted_burden"].mean()

    det_per_sample = df2[df2["n_hits"] > 0].copy()
    det_per_sample["detected"] = det_per_sample["Relative_abundance"] > 0
    mag_count_per_sample = det_per_sample.groupby(["SampleID", "Water type"])["detected"].sum().reset_index()
    mean_mag_count_by_water = mag_count_per_sample.groupby("Water type")["detected"].mean()

    burden_lines = [
        "# Abundance-weighted community gene-copy burden (post-exposure samples, D6-D89)",
        "",
        "Mean abundance-weighted panel-hit-gene-copy burden per sample, by water type:",
        burden_by_water.to_string(),
        "",
        f"Percent change OSPW vs ROW: {(burden_by_water['OSPW'] / burden_by_water['ROW'] - 1) * 100:.1f}%",
        "",
        "Mean burden per sample, by plant species x water type:",
        burden_by_plant_water.to_string(),
        "",
        "Mean detected gene-carrying MAGs per sample, by water type:",
        mean_mag_count_by_water.to_string(),
    ]
    (out_dir / "MAG_ROW_vs_OSPW_abundance_weighted_burden.txt").write_text("\n".join(burden_lines) + "\n")
    sample_burden.to_csv(out_dir / "MAG_ROW_vs_OSPW_sample_burden.tsv", sep="\t", index=False)

    stats_lines = [
        "# ROW vs. OSPW MAG gene-content comparison - statistical tests (post-exposure samples, D6-D89)",
        "",
        f"Mann-Whitney U, ROW-only vs OSPW-only panel hits: U={u_hits_row_ospw:.1f}, p={p_hits_row_ospw:.3f} (n={len(row_only)},{len(ospw_only)})",
        f"Mann-Whitney U, ROW-only vs OSPW-only relevant clusters: U={u_clust_row_ospw:.1f}, p={p_clust_row_ospw:.3f}",
        f"Mann-Whitney U, ROW-only vs both panel hits: U={u_hits_row_both:.1f}, p={p_hits_row_both:.3f} (n={len(row_only)},{len(both)})",
        f"Mann-Whitney U, ROW-only vs both relevant clusters: U={u_clust_row_both:.1f}, p={p_clust_row_both:.3f}",
        f"Spearman, within shared MAGs (n={len(both)}): OSPW/ROW log2 enrichment vs panel hits: rho={r_hits:.3f}, p={p_r_hits:.3f}",
        f"Spearman, within shared MAGs (n={len(both)}): OSPW/ROW log2 enrichment vs relevant clusters: rho={r_clust:.3f}, p={p_r_clust:.3f}",
    ]
    (out_dir / "MAG_ROW_vs_OSPW_statistical_tests.txt").write_text("\n".join(stats_lines) + "\n")
    print("\n".join(stats_lines))
    print()
    print(summary.to_string(index=False))
    print()
    print("\n".join(repertoire_lines))
    print()
    print("\n".join(burden_lines))

    # ---- Figure: 2x2 panels ----
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    labels = [f"ROW-only\n(n={len(row_only)})", f"OSPW-only\n(n={len(ospw_only)})", f"Both\n(n={len(both)})"]
    colors = ["#4393c3", "#d6604d", "#878787"]
    data = [row_only["panel_hits"], ospw_only["panel_hits"], both["panel_hits"]]
    bp = axes[0, 0].boxplot(data, tick_labels=labels, patch_artist=True, showfliers=True, widths=0.55)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    for median in bp["medians"]:
        median.set_color("#222222")
    axes[0, 0].set_ylabel("Reference-panel best-hit genes per MAG")
    axes[0, 0].set_title("A. Gene content by water-treatment\ndetection group", fontsize=11)
    axes[0, 0].spines[["top", "right"]].set_visible(False)

    axes[0, 1].scatter(both["log2fc_OSPW_over_ROW"], both["panel_hits"], s=18, alpha=0.55, color="#4d4d4d", edgecolors="none")
    axes[0, 1].set_xlabel("log2(mean OSPW abundance / mean ROW abundance)")
    axes[0, 1].set_ylabel("Reference-panel best-hit genes per MAG")
    axes[0, 1].set_title(f"B. Within shared MAGs (n={len(both)}): no relationship\nto OSPW/ROW enrichment (rho={r_hits:.2f}, p={p_r_hits:.2f})", fontsize=10.5)
    axes[0, 1].axvline(0, color="#bbbbbb", linewidth=1, linestyle="--", zorder=0)
    axes[0, 1].spines[["top", "right"]].set_visible(False)

    venn_vals = [len(row_genes - ospw_genes), len(row_genes & ospw_genes), len(ospw_genes - row_genes)]
    venn_labels = [f"ROW only\n({venn_vals[0]} genes)", f"Shared\n({venn_vals[1]} genes)", f"OSPW only\n({venn_vals[2]} genes)"]
    axes[1, 0].bar(venn_labels, venn_vals, color=["#4393c3", "#878787", "#d6604d"], alpha=0.75)
    axes[1, 0].set_ylabel("Distinct reference-panel gene names")
    axes[1, 0].set_title("C. Gene-level repertoire overlap", fontsize=11)
    axes[1, 0].spines[["top", "right"]].set_visible(False)
    for i, v in enumerate(venn_vals):
        axes[1, 0].text(i, v + 5, str(v), ha="center", fontsize=10)

    plant_water = burden_by_plant_water.reset_index()
    plants = ["Typha", "Juncus"]
    x = np.arange(len(plants))
    width = 0.35
    row_vals = [plant_water[(plant_water["Plant species"] == p) & (plant_water["Water type"] == "ROW")]["weighted_burden"].iloc[0] for p in plants]
    ospw_vals = [plant_water[(plant_water["Plant species"] == p) & (plant_water["Water type"] == "OSPW")]["weighted_burden"].iloc[0] for p in plants]
    axes[1, 1].bar(x - width / 2, row_vals, width, label="ROW", color="#4393c3", alpha=0.75)
    axes[1, 1].bar(x + width / 2, ospw_vals, width, label="OSPW", color="#d6604d", alpha=0.75)
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(plants)
    axes[1, 1].set_ylabel("Mean abundance-weighted gene-copy\nburden per sample")
    axes[1, 1].set_title("D. Community-level burden diverges\nby plant species under OSPW", fontsize=11)
    axes[1, 1].legend(frameon=False, fontsize=9)
    axes[1, 1].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_dir / "MAG_ROW_vs_OSPW_gene_content_figure.png", dpi=300)
    fig.savefig(out_dir / "MAG_ROW_vs_OSPW_gene_content_figure.pdf")
    fig.savefig(out_dir / "MAG_ROW_vs_OSPW_gene_content_figure.svg")
    print(f"\nWrote {out_dir / 'MAG_ROW_vs_OSPW_gene_content_figure.png'}")
    print(f"Wrote {out_dir / 'MAG_ROW_vs_OSPW_gene_content_figure.svg'}")


if __name__ == "__main__":
    main()
