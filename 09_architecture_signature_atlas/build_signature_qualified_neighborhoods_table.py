#!/usr/bin/env python3
# Title          : build_signature_qualified_neighborhoods_table.py
# Description    : Per-genome counts of signature-qualified candidate neighborhoods
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/28
# Usage          : python3 build_signature_qualified_neighborhoods_table.py

import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

LITERATURE = ["bad_ali_chc", "catechol_ortho", "phenylacetate", "protocatechuate", "beta_ox_core"]
CORE = [
    "cyclohexanecarboxylate_core", "catechol_meta_core", "aromatics_core", "gallate_core",
    "hydroxyquinol_core", "gentisate_core", "cyclohexylacetate_core", "alkanes_core",
    "alkenes_core", "benzoyl_coa_core", "naphthalene_core", "oxalate_core", "plastics_core",
    "pyrogallol_core", "tannin_core", "transportation_core",
]
ARCHITECTURES_21 = LITERATURE + CORE
assert len(ARCHITECTURES_21) == 21

TIER_LABEL = {"MAG": "Plant-root MAG", "Ref": "NCBI reference genome", "Colla": "OSPW-associated comparison genome"}

# genome -> {"tier": ..., "clusters": {cluster_id: set(architecture_names)}}
genomes = defaultdict(lambda: {"tier": None, "clusters": defaultdict(set)})

for arch in ARCHITECTURES_21:
    manifest = HERE / f"{arch}_manifest.tsv"
    with manifest.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = row["Genome"]
            genomes[g]["tier"] = row["Tier"]
            genomes[g]["clusters"][row["Cluster ID"]].add(arch)

rows = []
for genome, data in genomes.items():
    n = len(data["clusters"])
    if n == 0:
        continue
    archs_hit = sorted({a for archs in data["clusters"].values() for a in archs})
    rows.append({
        "Genome": genome,
        "Genome-source tier": TIER_LABEL[data["tier"]],
        "Signature-qualified candidate neighborhoods": n,
        "Matching architectures": "; ".join(archs_hit),
    })

rows.sort(key=lambda r: (-r["Signature-qualified candidate neighborhoods"], r["Genome"]))

out = HERE / "signature_qualified_neighborhoods_by_genome.tsv"
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["Genome", "Genome-source tier",
                                       "Signature-qualified candidate neighborhoods",
                                       "Matching architectures"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# Tier-level summary, matching the Results-text sentence structure
by_tier = defaultdict(list)
for r in rows:
    by_tier[r["Genome-source tier"]].append(r["Signature-qualified candidate neighborhoods"])

print(f"Wrote {out} ({len(rows)} genomes with >=1 signature-qualified neighborhood)")
print()
for tier_label, counts in by_tier.items():
    print(f"{tier_label}: {len(counts)} genomes carrying >=1 signature-qualified neighborhood "
          f"(mean {sum(counts)/len(counts):.2f}, median {sorted(counts)[len(counts)//2] if len(counts)%2 else (sorted(counts)[len(counts)//2-1]+sorted(counts)[len(counts)//2])/2})")
print()
print("Top 10 overall:")
for r in rows[:10]:
    print(f"  {r['Genome']} ({r['Genome-source tier']}): {r['Signature-qualified candidate neighborhoods']}")
