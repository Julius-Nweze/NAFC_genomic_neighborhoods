#!/usr/bin/env python3
# Title          : render_missing_cluster_svgs.py
# Description    : Render gene-map SVGs for neighborhoods selected by a signature but not yet rendered
# Author         : Julius Eyiuche Nweze
# Date           : 2026/08/21
# Usage          : python3 render_missing_cluster_svgs.py

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

OPERON_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(OPERON_DIR))
sys.path.insert(0, str(OPERON_DIR.parent / "07_neighborhood_reconstruction"))
from build_nafc_gene_panel_cluster_svgs import render_cluster_svg, slugify  # noqa: E402

CLUSTERS_TSV = OPERON_DIR / "NAFC_gene_panel_neighborhood_summary.tsv"
MEMBERS_TSV = OPERON_DIR / "NAFC_gene_panel_neighborhood_members.tsv"
SVG_DIR = OPERON_DIR / "svgs"

_cluster_rows_cache = None
_members_cache = None


def _load():
    global _cluster_rows_cache, _members_cache
    if _cluster_rows_cache is not None:
        return
    with open(CLUSTERS_TSV) as f:
        _cluster_rows_cache = {row["Cluster ID"]: row for row in csv.DictReader(f, delimiter="\t")}
    _members_cache = defaultdict(list)
    with open(MEMBERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            _members_cache[row["Cluster ID"]].append(row)


def ensure_rendered(cluster_ids: set[str]) -> list[str]:
    """Render any cluster_ids not already present in svgs/. Returns list of cluster IDs that
    could not be rendered (not found in the source tables)."""
    SVG_DIR.mkdir(exist_ok=True)
    to_render = [cid for cid in cluster_ids if not (SVG_DIR / f"{slugify(cid)}.svg").exists()]
    if not to_render:
        return []
    _load()
    failed = []
    rendered = 0
    for cid in to_render:
        cluster = _cluster_rows_cache.get(cid)
        records = _members_cache.get(cid)
        if not cluster or not records:
            failed.append(cid)
            continue
        records = sorted(records, key=lambda row: (int(row["Start"]), int(row["End"])))
        svg = render_cluster_svg(cluster, records)
        (SVG_DIR / f"{slugify(cid)}.svg").write_text(svg, encoding="utf-8")
        rendered += 1
    print(f"render_missing_cluster_svgs: rendered {rendered} new SVGs, {len(failed)} failed")
    return failed


if __name__ == "__main__":
    ids = {line.strip() for line in sys.stdin if line.strip()}
    failed = ensure_rendered(ids)
    if failed:
        print("FAILED:", failed)
