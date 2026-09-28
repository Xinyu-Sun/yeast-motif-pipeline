#!/usr/bin/env python3
"""Plot motif/domain dot-matrix figures from a Clustal Omega Newick tree.

This script builds SVG figures from:
- a guide tree (Newick)
- raw SGD domain hits
- protein metadata

It supports:
- all-protein and category-specific figures (Pol I / Pol II / Pol III)
- shared-domain filtering by protein count
- optional exclusion of MobiDBLite hits
- optional PNG/PDF export
- layout and style tuning through CLI arguments and/or JSON config
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TREE = REPO_ROOT / "data" / "input" / "clustalo-all48.phylotree"
DEFAULT_HITS = REPO_ROOT / "data" / "output" / "polymerase_tf_domain_hits_raw.csv"
DEFAULT_PROTEINS = REPO_ROOT / "data" / "output" / "polymerase_tf_proteins.csv"
DEFAULT_OUT = REPO_ROOT / "data" / "output" / "all48_motif_figure.svg"
DEFAULT_ORDER_CSV = REPO_ROOT / "data" / "output" / "all48_motif_order.csv"
SVG_DPI = 96.0

DEFAULT_CATEGORY_COLORS = {
    "Pol I": "#1f77b4",
    "Pol I axial": "#084594",
    "Pol I peri-axial": "#6baed6",
    "Pol II": "#d62728",
    "Pol III": "#2ca02c",
}
UNMAPPED_MOTIF_COLOR = "#8a8f98"
POL_I_SUBGROUP_LABELS = {
    "original_set": "Pol I axial",
    "original": "Pol I axial",
    "transcription_factor": "Pol I axial",
    "tf": "Pol I axial",
    "periaxial": "Pol I peri-axial",
    "peri_axial": "Pol I peri-axial",
    "paraaxial": "Pol I peri-axial",
    "axial": "Pol I axial",
}
DEFAULT_LEGEND_ORDER = ["Pol I axial", "Pol I peri-axial", "Pol I", "Pol II", "Pol III"]
DEFAULT_SUBTITLES = {
    "tree": "Ordered by Clustal Omega sequence similarity (neighbour-joining)",
    "motifs": "Ordered by shared motifs (Jaccard distance, average linkage)",
}
ORDER_CAPTIONS = {
    "tree": (
        "Columns follow a neighbour-joining tree built from a Clustal Omega alignment. When most proteins "
        "share under 20% sequence identity, the dendrogram is a display order, not evidence of "
        "evolutionary relationships."
    ),
    "motifs": (
        "Columns follow average-linkage clustering of the plotted motifs: the distance between two proteins "
        "is 1 minus the fraction of their motifs they share (Jaccard). Proteins with no plotted motif share "
        "nothing with any protein and are placed last. The order does not depend on input order."
    ),
}
LEGEND_BOX_WIDTH = 150.0

DEFAULT_STYLE = {
    "margin_left": 172.0,
    "margin_right": 16.0,
    "margin_top": 0.0,
    "margin_bottom": 12.0,
    "cell_width": 13.0,
    "cell_height": 14.0,
    "column_label_height": 94.0,
    "dendrogram_gap": 3.0,
    "dendrogram_height": 76.0,
    "dendrogram_bottom_padding": 8.0,
    "font_family": "Arial, Helvetica, sans-serif",
    "title_font_size": 19.0,
    "subtitle_font_size": 14.0,
    "protein_label_font_size": 8.7,
    "domain_label_font_size": 9.5,
    "legend_title_font_size": 10.5,
    "legend_font_size": 9.5,
    "title_baseline_offset": -3.0,
    "title_subtitle_gap": 18.0,
    "subtitle_matrix_gap": 8.0,
    "title_line_half_length": 88.0,
    "title_line_offset": 9.0,
    "title_line_width": 0.0,
    "grid_color": "#666666",
    "grid_width": 0.45,
    "frame_color": "#000000",
    "frame_width": 0.8,
    "tree_stroke_color": "#000000",
    "tree_stroke_width": 1.1,
    "dot_fill": None,
    "dot_stroke": None,
    "dot_radius": 3.55,
    "dot_stroke_width": 0.8,
    "background_color": "#ffffff",
    "legend_position": "auto",
    "legend_x": None,
    "legend_y": None,
    "legend_swatch_size": 9.0,
    "legend_row_gap": 12.0,
    "legend_title": "Category",
    "category_colors": DEFAULT_CATEGORY_COLORS,
    "unmapped_motif_color": UNMAPPED_MOTIF_COLOR,
    "auto_margin_left": True,
    "min_matrix_width": 156.0,
    "min_figure_width": 520.0,
}


@dataclass
class Node:
    name: Optional[str] = None
    length: float = 0.0
    children: List["Node"] = field(default_factory=list)
    height: float = 0.0
    x: float = 0.0
    y: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def parse(self) -> Node:
        node = self._parse_subtree()
        self._skip_ws()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self._skip_ws()
        if self.i != len(self.text):
            raise ValueError(f"Unexpected trailing Newick content at position {self.i}")
        return node

    def _parse_subtree(self) -> Node:
        self._skip_ws()
        if self._peek() == "(":
            self.i += 1
            children = []
            while True:
                children.append(self._parse_subtree())
                self._skip_ws()
                ch = self._peek()
                if ch == ",":
                    self.i += 1
                    continue
                if ch == ")":
                    self.i += 1
                    break
                raise ValueError(f"Expected ',' or ')' in Newick at position {self.i}, got {ch!r}")
            name = self._parse_name()
            length = self._parse_length()
            return Node(name=name or None, length=length, children=children)

        name = self._parse_name()
        length = self._parse_length()
        return Node(name=name or None, length=length)

    def _parse_name(self) -> str:
        self._skip_ws()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start:self.i].strip()

    def _parse_length(self) -> float:
        self._skip_ws()
        if self._peek() != ":":
            return 0.0
        self.i += 1
        self._skip_ws()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        raw = self.text[start:self.i].strip()
        return float(raw) if raw else 0.0

    def _skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _peek(self) -> str:
        return self.text[self.i] if self.i < len(self.text) else ""


def iter_leaves(node: Node) -> Iterable[Node]:
    if node.is_leaf:
        yield node
    else:
        for child in node.children:
            yield from iter_leaves(child)


def normalize_tree_label(label: str) -> str:
    label = label.strip()
    if label.startswith("sp|") or label.startswith("tr|"):
        parts = label.split("|")
        if len(parts) >= 3:
            return parts[2].split()[0].strip()
    return label


def compute_heights(node: Node) -> float:
    if node.is_leaf:
        node.height = 0.0
        return 0.0
    child_heights = [compute_heights(child) + child.length for child in node.children]
    node.height = max(child_heights) if child_heights else 0.0
    return node.height


def assign_positions(node: Node, leaf_x: Dict[str, float], tree_top: float, tree_height_px: float, max_height: float) -> None:
    if node.is_leaf:
        key = normalize_tree_label(node.name or "")
        if key not in leaf_x:
            raise KeyError(f"Leaf {node.name!r} missing x position")
        node.x = leaf_x[key]
        node.y = tree_top
        return

    for child in node.children:
        assign_positions(child, leaf_x, tree_top, tree_height_px, max_height)

    node.x = sum(child.x for child in node.children) / len(node.children)
    scale = 0.0 if max_height == 0 else tree_height_px / max_height
    node.y = tree_top + node.height * scale


def iter_tree_segments(node: Node) -> Iterable[tuple[float, float, float, float]]:
    if node.is_leaf:
        return

    child_x = [child.x for child in node.children]
    yield min(child_x), node.y, max(child_x), node.y
    for child in node.children:
        yield child.x, child.y, child.x, node.y
        yield from iter_tree_segments(child)


def draw_tree_segments(node: Node, stroke: str, stroke_width: float) -> List[str]:
    return [
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{stroke}" stroke-width="{stroke_width:.2f}" />'
        for x1, y1, x2, y2 in iter_tree_segments(node)
    ]


def prune_tree(node: Node, keep_labels: set[str]) -> Optional[Node]:
    if node.is_leaf:
        key = normalize_tree_label(node.name or "")
        if key in keep_labels:
            return Node(name=node.name, length=node.length)
        return None

    kept_children: List[Node] = []
    for child in node.children:
        pruned = prune_tree(child, keep_labels)
        if pruned is not None:
            kept_children.append(pruned)

    if not kept_children:
        return None

    if len(kept_children) == 1:
        only = kept_children[0]
        only.length += node.length
        return only

    return Node(name=node.name, length=node.length, children=kept_children)


def parse_header_identifier(header: str) -> Optional[str]:
    header = header.strip()
    if not header.startswith(">"):
        return None
    parts = header[1:].split("|")
    if len(parts) >= 3:
        return parts[2].split()[0].strip()
    token = header[1:].split()[0].strip()
    return token or None


def alias_forms(text: str) -> List[str]:
    text = (text or "").strip()
    if not text:
        return []
    compact = re.sub(r"[^A-Za-z0-9]+", "", text).lower()
    out = {text, text.upper(), text.lower(), compact}
    return [value for value in out if value]


def read_protein_metadata(path: Path) -> Dict[str, dict]:
    mapping: Dict[str, dict] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            header = row.get("fasta_header", "") or ""
            aliases = []
            identifier = parse_header_identifier(header)
            if identifier:
                aliases.extend(alias_forms(identifier))
            aliases.extend(alias_forms(row.get("standard_name", "") or ""))
            aliases.extend(alias_forms(row.get("protein", "") or ""))
            aliases.extend(alias_forms(row.get("y_name", "") or ""))
            token = header.lstrip(">").split()[0] if header else ""
            aliases.extend(alias_forms(token))
            for alias in aliases:
                mapping[alias] = row
    return mapping


def resolve_protein_row(label: str, protein_meta: Dict[str, dict]) -> Optional[dict]:
    key = normalize_tree_label(label)
    for alias in alias_forms(key):
        if alias in protein_meta:
            return protein_meta[alias]
    return None


def display_group(row: dict) -> str:
    category = (row.get("category") or "").strip()
    subgroup = (row.get("pol_i_subgroup") or "").strip().lower().replace(" ", "_").replace("-", "_")
    if category == "Pol I" and subgroup:
        return POL_I_SUBGROUP_LABELS.get(subgroup, f"Pol I {subgroup.replace('_', ' ')}")
    return category


def ordered_display_groups(rows: List[dict]) -> List[str]:
    present = {display_group(row) for row in rows if display_group(row)}
    ordered = [group for group in DEFAULT_LEGEND_ORDER if group in present]
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def group_color(row: dict, category_colors: dict[str, str]) -> str:
    return category_colors.get(display_group(row), category_colors.get(row.get("category", ""), "#000000"))


def dot_colors(row: dict, style: dict, category_colors: dict[str, str]) -> tuple[str, str]:
    fill = str(style["dot_fill"] or group_color(row, category_colors))
    stroke = str(style["dot_stroke"] or fill)
    return fill, stroke


def read_domain_hits(
    hits_path: Path,
    proteins_path: Path,
    protein_meta: Dict[str, dict],
    tree_labels: List[str],
    category: Optional[str],
    min_proteins: int,
    exclude_mobidblite: bool,
) -> tuple[list[str], dict[str, dict], list[dict], list[str]]:
    with proteins_path.open(newline="", encoding="utf-8-sig") as handle:
        protein_rows = list(csv.DictReader(handle))

    selected_rows = [row for row in protein_rows if not category or row["category"] == category]
    selected_y_names = {row["y_name"] for row in selected_rows}

    tree_rows = []
    unresolved_tree = []
    for label in tree_labels:
        row = resolve_protein_row(label, protein_meta)
        if not row:
            unresolved_tree.append(label)
            continue
        if row["y_name"] in selected_y_names:
            tree_rows.append(row)

    if unresolved_tree:
        raise SystemExit(f"Unresolved tree labels: {unresolved_tree}")

    tree_y_names = {row["y_name"] for row in tree_rows}
    missing_from_tree = [
        row["standard_name"] or row["protein"] or row["y_name"]
        for row in selected_rows
        if row["y_name"] not in tree_y_names
    ]

    domains_by_y: Dict[str, set[str]] = defaultdict(set)
    with hits_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["y_name"] not in tree_y_names:
                continue
            if category and row["category"] != category:
                continue

            raw_key = (row.get("raw_domain_key") or "").strip()
            source = (row.get("source") or "").strip()
            if not raw_key:
                continue
            if exclude_mobidblite and (source == "MobiDBLite" or "MobiDBLite" in raw_key):
                continue
            domains_by_y[row["y_name"]].add(raw_key)

    domain_presence = Counter()
    for domain_set in domains_by_y.values():
        for domain in domain_set:
            domain_presence[domain] += 1

    domain_keys = sorted([domain for domain, count in domain_presence.items() if count >= min_proteins])
    matrix_by_y: Dict[str, dict] = {}
    for row in tree_rows:
        out = {
            "category": row["category"],
            "protein": row["protein"],
            "y_name": row["y_name"],
            "standard_name": row.get("standard_name", ""),
            "pol_i_subgroup": row.get("pol_i_subgroup", ""),
        }
        present = domains_by_y.get(row["y_name"], set())
        for domain in domain_keys:
            out[domain] = 1 if domain in present else 0
        matrix_by_y[row["y_name"]] = out

    return domain_keys, matrix_by_y, tree_rows, missing_from_tree


def is_harmonized_interpro_mapping(row: dict) -> bool:
    accession = (row.get("harmonized_accession") or "").strip()
    return bool(
        (row.get("harmonization_status") or "").strip() == "mapped"
        and (row.get("harmonized_key") or "").strip().startswith("InterPro::")
        and accession
        and ";" not in accession
    )


def read_deduped_domain_hits(
    hits_path: Path,
    proteins_path: Path,
    protein_meta: Dict[str, dict],
    tree_labels: List[str],
    category: Optional[str],
    min_proteins: int,
    exclude_mobidblite: bool,
) -> tuple[list[str], dict[str, dict], list[dict], list[str], dict[str, dict]]:
    with proteins_path.open(newline="", encoding="utf-8-sig") as handle:
        protein_rows = list(csv.DictReader(handle))

    selected_rows = [row for row in protein_rows if not category or row["category"] == category]
    selected_y_names = {row["y_name"] for row in selected_rows}

    tree_rows = []
    unresolved_tree = []
    for label in tree_labels:
        row = resolve_protein_row(label, protein_meta)
        if not row:
            unresolved_tree.append(label)
            continue
        if row["y_name"] in selected_y_names:
            tree_rows.append(row)

    if unresolved_tree:
        raise SystemExit(f"Unresolved tree labels: {unresolved_tree}")

    tree_y_names = {row["y_name"] for row in tree_rows}
    missing_from_tree = [
        row["standard_name"] or row["protein"] or row["y_name"]
        for row in selected_rows
        if row["y_name"] not in tree_y_names
    ]

    domains_by_y: Dict[str, set[str]] = defaultdict(set)
    key_meta: dict[str, dict] = {}
    support_raw_keys: Dict[str, set[str]] = defaultdict(set)
    support_sources: Dict[str, set[str]] = defaultdict(set)

    with hits_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["y_name"] not in tree_y_names:
                continue
            if category and row["category"] != category:
                continue

            raw_key = (row.get("raw_domain_key") or "").strip()
            source = (row.get("source") or "").strip()
            if not raw_key:
                continue
            if exclude_mobidblite and (source == "MobiDBLite" or "MobiDBLite" in raw_key):
                continue

            mapped = is_harmonized_interpro_mapping(row)
            if mapped:
                key = (row.get("harmonized_key") or "").strip()
                label = (row.get("harmonized_name") or "").strip() or key
                accession = (row.get("harmonized_accession") or "").strip()
            else:
                key = f"Unmapped::{raw_key}"
                label = (row.get("domain_description") or "").strip() or raw_key
                accession = raw_key
            if label == "-":
                label = raw_key

            domains_by_y[row["y_name"]].add(key)
            support_raw_keys[key].add(raw_key)
            if source:
                support_sources[key].add(source)
            key_meta.setdefault(
                key,
                {
                    "status": "mapped" if mapped else "unmapped",
                    "label": label,
                    "accession": accession,
                },
            )

    domain_presence = Counter()
    for domain_set in domains_by_y.values():
        for domain in domain_set:
            domain_presence[domain] += 1

    def sort_key(key: str) -> tuple[int, str, str]:
        meta = key_meta[key]
        # Sort by the label drawn on the y-axis so unmapped (grey) rows read alphabetically too.
        return (0 if meta["status"] == "mapped" else 1, short_deduped_label(key, meta).lower(), key.lower())

    domain_keys = sorted([domain for domain, count in domain_presence.items() if count >= min_proteins], key=sort_key)
    for key in domain_keys:
        key_meta[key]["protein_count"] = domain_presence[key]
        key_meta[key]["support_raw_keys"] = sorted(support_raw_keys[key])
        key_meta[key]["support_sources"] = sorted(support_sources[key])

    matrix_by_y: Dict[str, dict] = {}
    for row in tree_rows:
        out = {
            "category": row["category"],
            "protein": row["protein"],
            "y_name": row["y_name"],
            "standard_name": row.get("standard_name", ""),
            "pol_i_subgroup": row.get("pol_i_subgroup", ""),
        }
        present = domains_by_y.get(row["y_name"], set())
        for domain in domain_keys:
            out[domain] = 1 if domain in present else 0
        matrix_by_y[row["y_name"]] = out

    return domain_keys, matrix_by_y, tree_rows, missing_from_tree, key_meta


def short_domain_label(raw: str) -> str:
    return raw.split("::", 1)[1] if "::" in raw else raw


def short_deduped_label(key: str, meta: dict) -> str:
    label = str(meta.get("label") or key)
    if meta.get("status") == "mapped":
        accession = str(meta.get("accession") or key.replace("InterPro::", ""))
        return f"{label} ({accession})"
    return short_domain_label(key.replace("Unmapped::", ""))


def domain_display_label(domain: str, domain_meta: Optional[dict[str, dict]] = None) -> str:
    meta = (domain_meta or {}).get(domain)
    return short_deduped_label(domain, meta) if meta else short_domain_label(domain)


def estimate_svg_text_width(text: str, font_size: float) -> float:
    return len(str(text)) * font_size * 0.58


def escape(text: str) -> str:
    return html.escape(text, quote=True)


def has_hit(row: dict, domain: str) -> bool:
    return str(row.get(domain, "0") or "0").strip() not in {"0", "", "0.0"}


def motif_cluster_tree(rows: List[dict], labels: List[str], domain_keys: List[str]) -> Node:
    """Average-linkage (UPGMA) tree on the Jaccard distance between plotted motif sets.

    Distance = 1 - shared motifs / motifs in either protein. A protein with no
    plotted motif is at distance 1 from every protein, including other empty ones.
    Ties are broken by protein label, and children are ordered with motif-carrying
    and larger clusters first, so the result does not depend on input order.
    """
    motif_sets = [frozenset(domain for domain in domain_keys if has_hit(row, domain)) for row in rows]

    def jaccard_distance(a: frozenset, b: frozenset) -> float:
        union = a | b
        return 1.0 - len(a & b) / len(union) if union else 1.0

    clusters: Dict[int, dict] = {
        i: {"node": Node(name=labels[i]), "size": 1, "height": 0.0, "key": labels[i], "with_motifs": int(bool(motif_sets[i]))}
        for i in range(len(rows))
    }
    dist: Dict[tuple, float] = {}
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            dist[(i, j)] = jaccard_distance(motif_sets[i], motif_sets[j])

    def d(a: int, b: int) -> float:
        return dist[(a, b) if a < b else (b, a)]

    def display_order(cluster: dict) -> tuple:
        return (-cluster["with_motifs"], -cluster["size"], cluster["key"])

    next_id = len(rows)
    while len(clusters) > 1:
        ids = sorted(clusters)
        best = min(
            ((a, b) for idx, a in enumerate(ids) for b in ids[idx + 1 :]),
            key=lambda pair: (
                round(d(*pair), 12),
                min(clusters[pair[0]]["key"], clusters[pair[1]]["key"]),
                max(clusters[pair[0]]["key"], clusters[pair[1]]["key"]),
            ),
        )
        a, b = best
        height = d(a, b)
        if height >= 1.0 - 1e-12:
            # The remaining groups share no motif at all, so any pairing would be arbitrary.
            # Join them in one node: motif groups first (largest first), then proteins without motifs.
            remaining = sorted(clusters.values(), key=display_order)
            for child in remaining:
                child["node"].length = 1.0 - child["height"]
            return Node(children=[child["node"] for child in remaining])
        first, second = sorted((clusters[a], clusters[b]), key=display_order)
        for child in (first, second):
            child["node"].length = height - child["height"]
        merged = {
            "node": Node(children=[first["node"], second["node"]]),
            "size": first["size"] + second["size"],
            "height": height,
            "key": min(first["key"], second["key"]),
            "with_motifs": first["with_motifs"] + second["with_motifs"],
        }
        for other in ids:
            if other in (a, b):
                continue
            dist[(other, next_id)] = (clusters[a]["size"] * d(a, other) + clusters[b]["size"] * d(b, other)) / merged["size"]
        del clusters[a], clusters[b]
        clusters[next_id] = merged
        next_id += 1
    return next(iter(clusters.values()))["node"]


def protein_labels_in_table(path: Path, category: Optional[str]) -> List[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            row["y_name"]
            for row in csv.DictReader(handle)
            if (row.get("y_name") or "").strip() and (not category or row.get("category") == category)
        ]


def prepare_plot_data(args: argparse.Namespace) -> dict:
    """Select proteins and motifs, then order the columns by the tree file or by shared motifs."""
    protein_meta = read_protein_metadata(args.proteins)
    if args.order_by == "tree":
        tree: Optional[Node] = NewickParser(args.tree.read_text(encoding="utf-8")).parse()
        raw_labels = [leaf.name or "" for leaf in iter_leaves(tree)]
    else:
        tree = None
        raw_labels = protein_labels_in_table(args.proteins, args.category)
    labels = [normalize_tree_label(label) for label in raw_labels]

    reader = read_deduped_domain_hits if args.dedupe_harmonized else read_domain_hits
    result = reader(
        hits_path=args.hits,
        proteins_path=args.proteins,
        protein_meta=protein_meta,
        tree_labels=labels,
        category=args.category,
        min_proteins=args.min_proteins,
        exclude_mobidblite=args.exclude_mobidblite,
    )
    domain_keys, matrix_by_y, missing_from_tree = result[0], result[1], result[3]
    domain_meta = result[4] if args.dedupe_harmonized else None

    rows: List[dict] = []
    kept: List[str] = []
    for raw_label, label in zip(raw_labels, labels):
        row = resolve_protein_row(raw_label, protein_meta) or resolve_protein_row(label, protein_meta)
        if not row or (args.category and row["category"] != args.category):
            continue
        y_name = row["y_name"]
        if y_name not in matrix_by_y:
            continue
        if domain_keys and args.hide_empty_columns and not any(has_hit(matrix_by_y[y_name], domain) for domain in domain_keys):
            continue
        rows.append(matrix_by_y[y_name])
        kept.append(label)
    if not rows:
        raise SystemExit("No proteins selected for plotting. Check --category and input files.")

    if tree is None:
        tree = motif_cluster_tree(rows, kept, domain_keys)
        position = {label: idx for idx, label in enumerate(kept)}
        kept = [normalize_tree_label(leaf.name or "") for leaf in iter_leaves(tree)]
        rows = [rows[position[label]] for label in kept]

    pruned = prune_tree(tree, set(kept))
    if pruned is None:
        raise SystemExit("Could not prune tree to selected proteins.")
    return {
        "tree": pruned,
        "labels": kept,
        "rows": rows,
        "domain_keys": domain_keys,
        "domain_meta": domain_meta,
        "missing_from_tree": missing_from_tree,
    }


def write_order_csv(path: Path, ordered_rows: List[dict], tree_labels: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_label",
                "y_name",
                "standard_name",
                "protein",
                "category",
                "pol_i_subgroup",
                "display_group",
            ],
        )
        writer.writeheader()
        for label, row in zip(tree_labels, ordered_rows):
            writer.writerow(
                {
                    "tree_label": label,
                    "y_name": row["y_name"],
                    "standard_name": row.get("standard_name", ""),
                    "protein": row.get("protein", ""),
                    "category": row.get("category", ""),
                    "pol_i_subgroup": row.get("pol_i_subgroup", ""),
                    "display_group": display_group(row),
                }
            )


def merge_style(base: dict, patch: dict) -> dict:
    merged = dict(base)
    for key, value in patch.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            nested = dict(merged[key])
            nested.update(value)
            merged[key] = nested
        else:
            merged[key] = value
    return merged


def load_style_config(path: Optional[Path]) -> dict:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return json.load(handle)


def manifest_path_for(figure_path: Path) -> Path:
    return figure_path.with_name(figure_path.stem + ".plot.json")


def portable_relpath(path: Optional[Path], start: Path) -> Optional[str]:
    """Return a forward-slash path relative to ``start`` (absolute if on another drive)."""
    if path is None:
        return None
    try:
        return Path(os.path.relpath(path.resolve(), start.resolve())).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def write_plot_manifest(
    args: argparse.Namespace,
    renderer: str,
    ordered_rows: List[dict],
    domain_keys: List[str],
    missing_from_tree: List[str],
) -> Path:
    """Write a small JSON sidecar describing one figure, used by build_output_index.py."""
    path = manifest_path_for(args.out)
    base_dir = path.parent
    dot_count = sum(1 for row in ordered_rows for domain in domain_keys if has_hit(row, domain))
    manifest = {
        "schema": "yeast-motif-pipeline/plot-manifest/v1",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "renderer": renderer,
        "title": args.title,
        "subtitle": args.subtitle,
        "category": args.category or "All",
        "order_by": args.order_by,
        "caption": ORDER_CAPTIONS[args.order_by],
        "options": {
            "min_proteins": args.min_proteins,
            "exclude_mobidblite": args.exclude_mobidblite,
            "dedupe_harmonized": args.dedupe_harmonized,
            "hide_empty_columns": args.hide_empty_columns,
            "export_dpi": args.export_dpi if args.png else None,
        },
        # Inputs are recorded as given on the command line (normally relative to the repository root).
        "inputs": {
            "tree": args.tree.as_posix() if args.order_by == "tree" else None,
            "hits": args.hits.as_posix(),
            "proteins": args.proteins.as_posix(),
        },
        "outputs": {
            "svg": portable_relpath(args.out, base_dir),
            "png": portable_relpath(args.png, base_dir),
            "pdf": portable_relpath(args.pdf, base_dir),
            "order_csv": portable_relpath(args.order_csv, base_dir),
        },
        "counts": {
            "proteins_plotted": len(ordered_rows),
            "motifs_plotted": len(domain_keys),
            "dots_plotted": dot_count,
            "display_groups": dict(Counter(display_group(row) for row in ordered_rows)),
        },
        "missing_from_tree": list(missing_from_tree),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path


def parse_category_color_overrides(overrides: List[str]) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for raw in overrides:
        if "=" not in raw:
            raise SystemExit(f"Invalid --category-color value {raw!r}. Use format 'Pol I=#1f77b4'.")
        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not value:
            raise SystemExit(f"Missing color for category {key!r}")
        parsed[key] = value
    return parsed


def top_legend_y(style: dict, width: float, legend_x: float) -> float:
    """Keep a top legend in the top margin unless it would overlap the centered title block.

    Moving it below the subtitle only when needed keeps it clear of the rotated
    protein labels on wide figures.
    """
    title_font_size = float(style["title_font_size"])
    subtitle_font_size = float(style["subtitle_font_size"])
    half_text_w = 12.0 + max(
        estimate_svg_text_width(str(style.get("_current_title") or ""), title_font_size),
        estimate_svg_text_width(str(style.get("_current_subtitle") or ""), subtitle_font_size),
    ) / 2.0
    if legend_x + LEGEND_BOX_WIDTH < width / 2.0 - half_text_w or legend_x > width / 2.0 + half_text_w:
        return 8.0
    subtitle_y = (
        float(style["margin_top"])
        + title_font_size
        + float(style.get("title_baseline_offset", 0.0))
        + float(style["title_subtitle_gap"])
    )
    return subtitle_y + 15.0


def determine_legend_position(
    style: dict,
    width: float,
    matrix_left: float,
    matrix_w: float,
    matrix_top: float,
    matrix_h: float,
    tree_top: float,
    tree_h: float,
) -> tuple[float, float]:
    position = str(style["legend_position"]).lower()
    if position == "custom":
        if style.get("legend_x") is None or style.get("legend_y") is None:
            raise SystemExit("legend_position=custom requires legend_x and legend_y")
        return float(style["legend_x"]), float(style["legend_y"])

    if position == "auto":
        position = "top-left" if width < 420 else "top-right"

    if position == "top-left":
        return 12.0, top_legend_y(style, width, 12.0)
    if position == "top-right":
        legend_x = max(width - 165.0, matrix_left + matrix_w - 95.0)
        return legend_x, top_legend_y(style, width, legend_x)
    if position == "bottom-left":
        return 12.0, tree_top + tree_h + 2.0
    if position == "bottom-right":
        return max(width - 165.0, matrix_left + matrix_w - 95.0), tree_top + tree_h + 2.0

    raise SystemExit(f"Unsupported legend_position: {style['legend_position']!r}")


def prepare_layout(
    tree: Node,
    tree_labels: List[str],
    domain_keys: List[str],
    style: dict,
    domain_meta: Optional[dict[str, dict]] = None,
) -> dict:
    margin_left = float(style["margin_left"])
    margin_right = float(style["margin_right"])
    margin_top = float(style["margin_top"])
    margin_bottom = float(style["margin_bottom"])
    cell_w = float(style["cell_width"])
    cell_h = float(style["cell_height"])
    col_label_h = float(style["column_label_height"])
    dendro_gap = float(style["dendrogram_gap"])
    dendro_h = float(style["dendrogram_height"])
    dendro_bottom_pad = float(style["dendrogram_bottom_padding"])
    min_matrix_w = max(0.0, float(style.get("min_matrix_width", 0.0)))
    min_figure_w = max(0.0, float(style.get("min_figure_width", 0.0)))

    title_font_size = float(style["title_font_size"])
    title_baseline_offset = float(style.get("title_baseline_offset", 0.0))
    title_subtitle_gap = float(style["title_subtitle_gap"])
    subtitle_matrix_gap = float(style["subtitle_matrix_gap"])
    if style.get("auto_margin_left", True):
        label_font_size = float(style["domain_label_font_size"])
        longest_label_width = max(
            [estimate_svg_text_width(domain_display_label(domain, domain_meta), label_font_size) for domain in domain_keys],
            default=0.0,
        )
        margin_left = max(margin_left, longest_label_width + 18.0)

    label_count = max(len(tree_labels), 1)
    matrix_w = max(label_count * cell_w, min_matrix_w)
    cell_w = matrix_w / label_count
    matrix_h = max(len(domain_keys), 1) * cell_h

    title_y = margin_top + title_font_size + title_baseline_offset
    subtitle_y = title_y + title_subtitle_gap
    matrix_top = subtitle_y + subtitle_matrix_gap + col_label_h
    matrix_left = margin_left

    tree_top = matrix_top + matrix_h + dendro_gap
    text_min_width = max(
        estimate_svg_text_width(str(style.get("legend_title", "")), float(style["legend_title_font_size"])) + 34.0,
        estimate_svg_text_width("", title_font_size),
    )
    width = margin_left + matrix_w + margin_right
    required_width = max(
        min_figure_w,
        estimate_svg_text_width("", title_font_size),
        estimate_svg_text_width("", float(style["subtitle_font_size"])),
        text_min_width,
    )
    if title := style.get("_current_title"):
        required_width = max(required_width, estimate_svg_text_width(str(title), title_font_size) + 28.0)
    if subtitle := style.get("_current_subtitle"):
        required_width = max(required_width, estimate_svg_text_width(str(subtitle), float(style["subtitle_font_size"])) + 28.0)
    if width < required_width:
        margin_right += required_width - width
        width = required_width
    height = tree_top + dendro_h + margin_bottom

    leaf_x = {label: matrix_left + (idx + 0.5) * cell_w for idx, label in enumerate(tree_labels)}
    max_height = compute_heights(tree)
    assign_positions(
        tree,
        leaf_x=leaf_x,
        tree_top=tree_top,
        tree_height_px=max(0.0, dendro_h - dendro_bottom_pad),
        max_height=max_height,
    )

    return {
        "margin_left": margin_left,
        "margin_right": margin_right,
        "margin_top": margin_top,
        "margin_bottom": margin_bottom,
        "cell_w": cell_w,
        "cell_h": cell_h,
        "col_label_h": col_label_h,
        "dendro_gap": dendro_gap,
        "dendro_h": dendro_h,
        "dendro_bottom_pad": dendro_bottom_pad,
        "title_font_size": title_font_size,
        "subtitle_font_size": float(style["subtitle_font_size"]),
        "title_baseline_offset": title_baseline_offset,
        "title_subtitle_gap": title_subtitle_gap,
        "subtitle_matrix_gap": subtitle_matrix_gap,
        "matrix_w": matrix_w,
        "matrix_h": matrix_h,
        "matrix_top": matrix_top,
        "matrix_left": matrix_left,
        "tree_top": tree_top,
        "width": width,
        "height": height,
        "title_x": width / 2.0,
        "title_y": title_y,
        "subtitle_y": subtitle_y,
        "label_y": matrix_top - 2.0,
    }


def style_with_text_context(style: dict, title: str, subtitle: str) -> dict:
    contextual_style = dict(style)
    contextual_style["_current_title"] = title
    contextual_style["_current_subtitle"] = subtitle
    return contextual_style


def build_svg(
    tree: Node,
    tree_labels: List[str],
    domain_keys: List[str],
    ordered_rows: List[dict],
    out: Path,
    title: str,
    subtitle: str,
    show_legend: bool,
    style: dict,
    domain_meta: Optional[dict[str, dict]] = None,
) -> None:
    style = style_with_text_context(style, title, subtitle)
    layout = prepare_layout(tree, tree_labels, domain_keys, style, domain_meta=domain_meta)
    margin_left = layout["margin_left"]
    cell_w = layout["cell_w"]
    cell_h = layout["cell_h"]
    dendro_h = layout["dendro_h"]
    title_font_size = layout["title_font_size"]
    subtitle_font_size = layout["subtitle_font_size"]
    matrix_w = layout["matrix_w"]
    matrix_h = layout["matrix_h"]
    matrix_top = layout["matrix_top"]
    matrix_left = layout["matrix_left"]
    tree_top = layout["tree_top"]
    width = layout["width"]
    height = layout["height"]
    title_x = layout["title_x"]
    title_y = layout["title_y"]
    subtitle_y = layout["subtitle_y"]
    label_y = layout["label_y"]

    category_colors = dict(DEFAULT_CATEGORY_COLORS)
    category_colors.update(style.get("category_colors", {}))
    unmapped_color = str(style.get("unmapped_motif_color", UNMAPPED_MOTIF_COLOR))

    pieces: List[str] = []
    pieces.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.2f}" height="{height:.2f}" viewBox="0 0 {width:.2f} {height:.2f}">')
    pieces.append(f'<rect width="100%" height="100%" fill="{escape(str(style["background_color"]))}" />')

    pieces.append(
        f'<text x="{title_x:.2f}" y="{title_y:.2f}" text-anchor="middle" font-family="{escape(str(style["font_family"]))}" '
        f'font-size="{title_font_size:.2f}" font-weight="700">{escape(title)}</text>'
    )

    line_half = float(style["title_line_half_length"])
    line_y = title_y + float(style["title_line_offset"])
    line_w = float(style["title_line_width"])
    if line_w > 0:
        pieces.append(
            f'<line x1="{title_x - line_half:.2f}" y1="{line_y:.2f}" x2="{title_x + line_half:.2f}" y2="{line_y:.2f}" '
            f'stroke="{escape(str(style["frame_color"]))}" stroke-width="{line_w:.2f}" />'
        )

    pieces.append(
        f'<text x="{title_x:.2f}" y="{subtitle_y:.2f}" text-anchor="middle" font-family="{escape(str(style["font_family"]))}" '
        f'font-size="{subtitle_font_size:.2f}">{escape(subtitle)}</text>'
    )

    if show_legend:
        legend_x, legend_y = determine_legend_position(
            style=style,
            width=width,
            matrix_left=matrix_left,
            matrix_w=matrix_w,
            matrix_top=matrix_top,
            matrix_h=matrix_h,
            tree_top=tree_top,
            tree_h=dendro_h,
        )
        pieces.append(
            f'<text x="{legend_x:.2f}" y="{legend_y + 2:.2f}" font-family="{escape(str(style["font_family"]))}" '
            f'font-size="{float(style["legend_title_font_size"]):.2f}" font-weight="700">{escape(str(style["legend_title"]))}</text>'
        )
        swatch_size = float(style["legend_swatch_size"])
        row_gap = float(style["legend_row_gap"])
        for idx, category in enumerate(ordered_display_groups(ordered_rows)):
            y = legend_y + 12.0 + idx * row_gap
            color = category_colors.get(category, "#000000")
            pieces.append(
                f'<rect x="{legend_x:.2f}" y="{y - swatch_size + 2:.2f}" width="{swatch_size:.2f}" height="{swatch_size:.2f}" '
                f'fill="{escape(color)}" stroke="{escape(str(style["frame_color"]))}" stroke-width="0.4" />'
            )
            pieces.append(
                f'<text x="{legend_x + swatch_size + 5:.2f}" y="{y + 1:.2f}" font-family="{escape(str(style["font_family"]))}" '
                f'font-size="{float(style["legend_font_size"]):.2f}">{escape(category)}</text>'
            )

    for idx, row in enumerate(ordered_rows):
        short = row.get("standard_name") or row.get("protein") or row.get("y_name")
        x = matrix_left + (idx + 0.5) * cell_w
        fill = group_color(row, category_colors)
        pieces.append(
            f'<text x="0" y="0" transform="translate({x:.2f},{label_y:.2f}) rotate(-90)" text-anchor="start" '
            f'font-family="{escape(str(style["font_family"]))}" font-size="{float(style["protein_label_font_size"]):.2f}" '
            f'dominant-baseline="middle" fill="{escape(fill)}">{escape(str(short))}</text>'
        )

    grid_color = escape(str(style["grid_color"]))
    grid_w = float(style["grid_width"])
    for idx in range(len(tree_labels) + 1):
        x = matrix_left + idx * cell_w
        pieces.append(
            f'<line x1="{x:.2f}" y1="{matrix_top:.2f}" x2="{x:.2f}" y2="{matrix_top + matrix_h:.2f}" '
            f'stroke="{grid_color}" stroke-width="{grid_w:.2f}" />'
        )

    dot_radius = float(style["dot_radius"])
    dot_stroke_w = float(style["dot_stroke_width"])
    for row_idx, domain in enumerate(domain_keys):
        cy = matrix_top + (row_idx + 0.5) * cell_h
        meta = (domain_meta or {}).get(domain)
        label = domain_display_label(domain, domain_meta)
        label_fill = unmapped_color if meta and meta.get("status") == "unmapped" else "#000000"
        pieces.append(
            f'<text x="{matrix_left - 8:.2f}" y="{cy + 3:.2f}" text-anchor="end" '
            f'font-family="{escape(str(style["font_family"]))}" font-size="{float(style["domain_label_font_size"]):.2f}" '
            f'fill="{escape(label_fill)}">{escape(label)}</text>'
        )
        for col_idx, row in enumerate(ordered_rows):
            raw = str(row.get(domain, "0") or "0").strip()
            if raw in {"0", "", "0.0"}:
                continue
            cx = matrix_left + (col_idx + 0.5) * cell_w
            if meta and meta.get("status") == "unmapped":
                dot_fill = dot_stroke = unmapped_color
                mapped_note = "unmapped source motif; kept separate"
            else:
                dot_fill, dot_stroke = dot_colors(row, style, category_colors)
                mapped_note = "InterPro-deduped motif" if meta else "raw motif"
            support = ""
            if meta:
                support_values = meta.get("support_raw_keys", [])[:5]
                support = ", ".join(support_values)
                if len(meta.get("support_raw_keys", [])) > 5:
                    support += ", ..."
            tooltip = f'{label} in {row.get("standard_name") or row.get("protein") or row.get("y_name")}'
            if meta:
                tooltip += f" | {mapped_note}"
            if support:
                tooltip += f" | source IDs: {support}"
            pieces.append(
                f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{dot_radius:.2f}" fill="{escape(dot_fill)}" stroke="{escape(dot_stroke)}" stroke-width="{dot_stroke_w:.2f}">'
                f'<title>{escape(tooltip)}</title>'
                f'</circle>'
            )

    pieces.extend(
        draw_tree_segments(
            tree,
            stroke=escape(str(style["tree_stroke_color"])),
            stroke_width=float(style["tree_stroke_width"]),
        )
    )

    pieces.append(
        f'<rect x="{matrix_left:.2f}" y="{matrix_top:.2f}" width="{matrix_w:.2f}" height="{matrix_h:.2f}" fill="none" '
        f'stroke="{escape(str(style["frame_color"]))}" stroke-width="{float(style["frame_width"]):.2f}" />'
    )
    pieces.append("</svg>")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(pieces), encoding="utf-8")


def parse_color(value: str) -> tuple[int, int, int]:
    text = str(value).strip().lower()
    named = {
        "black": "#000000",
        "white": "#ffffff",
        "red": "#ff0000",
        "green": "#008000",
        "blue": "#0000ff",
        "transparent": "#ffffff",
    }
    text = named.get(text, text)
    if text.startswith("#") and len(text) == 4:
        text = "#" + "".join(ch * 2 for ch in text[1:])
    if not (text.startswith("#") and len(text) == 7):
        raise SystemExit(f"PNG/PDF export only supports hex colors like #178b1d, got {value!r}")
    try:
        return int(text[1:3], 16), int(text[3:5], 16), int(text[5:7], 16)
    except ValueError as exc:
        raise SystemExit(f"Invalid color for PNG/PDF export: {value!r}") from exc


def load_pillow():
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise SystemExit(
            "PNG export requires Pillow. Install it with `python3 -m pip install pillow`, "
            "or omit --png and use the SVG output."
        ) from exc
    return Image, ImageDraw, ImageFont


def load_reportlab():
    try:
        from reportlab.lib.colors import Color
        from reportlab.pdfgen import canvas
    except ImportError as exc:
        raise SystemExit(
            "PDF export requires ReportLab. Install it with `python3 -m pip install reportlab`, "
            "or omit --pdf and use the SVG output."
        ) from exc
    return Color, canvas


def scaled_font(image_font, size: float):
    try:
        return image_font.truetype("Arial.ttf", int(round(size)))
    except OSError:
        try:
            return image_font.truetype("DejaVuSans.ttf", int(round(size)))
        except OSError:
            return image_font.load_default()


def draw_png_text(draw, image_font, xy: tuple[float, float], text: str, size: float, fill, anchor: str = "la") -> None:
    draw.text(xy, text, font=scaled_font(image_font, size), fill=fill, anchor=anchor)


def draw_rotated_png_text(
    image_module,
    image_draw_module,
    image,
    image_draw,
    image_font,
    x: float,
    y: float,
    text: str,
    size: float,
    fill,
) -> None:
    font = scaled_font(image_font, size)
    bbox = image_draw.textbbox((0, 0), text, font=font)
    width = max(1, bbox[2] - bbox[0] + 4)
    height = max(1, bbox[3] - bbox[1] + 4)
    text_image = image_module.new("RGBA", (width, height), (255, 255, 255, 0))
    text_draw = image_draw_module.Draw(text_image)
    text_draw.text((2 - bbox[0], 2 - bbox[1]), text, font=font, fill=fill)
    rotated = text_image.rotate(90, expand=True)
    image.alpha_composite(rotated, (int(round(x)), int(round(y - rotated.height))))


def build_png(
    tree: Node,
    tree_labels: List[str],
    domain_keys: List[str],
    ordered_rows: List[dict],
    out: Path,
    title: str,
    subtitle: str,
    show_legend: bool,
    style: dict,
    dpi: float,
    domain_meta: Optional[dict[str, dict]] = None,
) -> None:
    if dpi <= 0:
        raise SystemExit("--export-dpi must be > 0")
    Image, ImageDraw, ImageFont = load_pillow()
    style = style_with_text_context(style, title, subtitle)
    layout = prepare_layout(tree, tree_labels, domain_keys, style, domain_meta=domain_meta)
    scale = dpi / SVG_DPI

    def sx(value: float) -> float:
        return value * scale

    category_colors = dict(DEFAULT_CATEGORY_COLORS)
    category_colors.update(style.get("category_colors", {}))
    unmapped_color = str(style.get("unmapped_motif_color", UNMAPPED_MOTIF_COLOR))
    background = parse_color(str(style["background_color"]))
    image = Image.new("RGBA", (int(round(sx(layout["width"]))), int(round(sx(layout["height"])))), (*background, 255))
    draw = ImageDraw.Draw(image)

    def color(value: str):
        return (*parse_color(value), 255)

    frame_color = color(str(style["frame_color"]))
    grid_color = color(str(style["grid_color"]))
    tree_color = color(str(style["tree_stroke_color"]))

    title_x = sx(layout["title_x"])
    title_y = sx(layout["title_y"])
    draw_png_text(draw, ImageFont, (title_x, title_y), title, sx(layout["title_font_size"]), (0, 0, 0, 255), anchor="mm")
    line_half = sx(float(style["title_line_half_length"]))
    line_y = sx(layout["title_y"] + float(style["title_line_offset"]))
    title_line_width = float(style["title_line_width"])
    if title_line_width > 0:
        draw.line((title_x - line_half, line_y, title_x + line_half, line_y), fill=frame_color, width=max(1, int(round(sx(title_line_width)))))
    draw_png_text(
        draw,
        ImageFont,
        (title_x, sx(layout["subtitle_y"])),
        subtitle,
        sx(layout["subtitle_font_size"]),
        (0, 0, 0, 255),
        anchor="mm",
    )

    matrix_left = layout["matrix_left"]
    matrix_top = layout["matrix_top"]
    matrix_w = layout["matrix_w"]
    matrix_h = layout["matrix_h"]
    cell_w = layout["cell_w"]
    cell_h = layout["cell_h"]

    if show_legend:
        legend_x, legend_y = determine_legend_position(
            style=style,
            width=layout["width"],
            matrix_left=matrix_left,
            matrix_w=matrix_w,
            matrix_top=matrix_top,
            matrix_h=matrix_h,
            tree_top=layout["tree_top"],
            tree_h=layout["dendro_h"],
        )
        draw_png_text(
            draw,
            ImageFont,
            (sx(legend_x), sx(legend_y + 2)),
            str(style["legend_title"]),
            sx(float(style["legend_title_font_size"])),
            (0, 0, 0, 255),
        )
        swatch_size = float(style["legend_swatch_size"])
        row_gap = float(style["legend_row_gap"])
        for idx, category in enumerate(ordered_display_groups(ordered_rows)):
            y = legend_y + 12.0 + idx * row_gap
            fill = color(category_colors.get(category, "#000000"))
            draw.rectangle(
                (sx(legend_x), sx(y - swatch_size + 2), sx(legend_x + swatch_size), sx(y + 2)),
                fill=fill,
                outline=frame_color,
                width=max(1, int(round(sx(0.4)))),
            )
            draw_png_text(
                draw,
                ImageFont,
                (sx(legend_x + swatch_size + 5), sx(y + 1)),
                category,
                sx(float(style["legend_font_size"])),
                (0, 0, 0, 255),
            )

    label_y = layout["label_y"]
    for idx, row in enumerate(ordered_rows):
        short = row.get("standard_name") or row.get("protein") or row.get("y_name")
        x = matrix_left + (idx + 0.5) * cell_w
        fill = color(group_color(row, category_colors))
        draw_rotated_png_text(
            Image,
            ImageDraw,
            image,
            draw,
            ImageFont,
            sx(x),
            sx(label_y),
            str(short),
            sx(float(style["protein_label_font_size"])),
            fill,
        )

    grid_width = max(1, int(round(sx(float(style["grid_width"])))))
    for idx in range(len(tree_labels) + 1):
        x = matrix_left + idx * cell_w
        draw.line((sx(x), sx(matrix_top), sx(x), sx(matrix_top + matrix_h)), fill=grid_color, width=grid_width)

    dot_radius = sx(float(style["dot_radius"]))
    dot_width = max(1, int(round(sx(float(style["dot_stroke_width"])))))
    for row_idx, domain in enumerate(domain_keys):
        cy = matrix_top + (row_idx + 0.5) * cell_h
        meta = (domain_meta or {}).get(domain)
        label = domain_display_label(domain, domain_meta)
        label_color = color(unmapped_color) if meta and meta.get("status") == "unmapped" else (0, 0, 0, 255)
        draw_png_text(
            draw,
            ImageFont,
            (sx(matrix_left - 8), sx(cy + 3)),
            label,
            sx(float(style["domain_label_font_size"])),
            label_color,
            anchor="ra",
        )
        for col_idx, row in enumerate(ordered_rows):
            raw = str(row.get(domain, "0") or "0").strip()
            if raw in {"0", "", "0.0"}:
                continue
            cx = matrix_left + (col_idx + 0.5) * cell_w
            if meta and meta.get("status") == "unmapped":
                dot_fill = dot_stroke = unmapped_color
            else:
                dot_fill, dot_stroke = dot_colors(row, style, category_colors)
            draw.ellipse(
                (sx(cx) - dot_radius, sx(cy) - dot_radius, sx(cx) + dot_radius, sx(cy) + dot_radius),
                fill=color(dot_fill),
                outline=color(dot_stroke),
                width=dot_width,
            )

    tree_width = max(1, int(round(sx(float(style["tree_stroke_width"])))))
    for x1, y1, x2, y2 in iter_tree_segments(tree):
        draw.line((sx(x1), sx(y1), sx(x2), sx(y2)), fill=tree_color, width=tree_width)

    frame_width = max(1, int(round(sx(float(style["frame_width"])))))
    draw.rectangle(
        (sx(matrix_left), sx(matrix_top), sx(matrix_left + matrix_w), sx(matrix_top + matrix_h)),
        outline=frame_color,
        width=frame_width,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(out, dpi=(dpi, dpi))


def build_pdf(
    tree: Node,
    tree_labels: List[str],
    domain_keys: List[str],
    ordered_rows: List[dict],
    out: Path,
    title: str,
    subtitle: str,
    show_legend: bool,
    style: dict,
    domain_meta: Optional[dict[str, dict]] = None,
) -> None:
    Color, canvas = load_reportlab()
    style = style_with_text_context(style, title, subtitle)
    layout = prepare_layout(tree, tree_labels, domain_keys, style, domain_meta=domain_meta)

    def pdf_y(value: float) -> float:
        return layout["height"] - value

    def set_color(c, value: str) -> None:
        red, green, blue = parse_color(value)
        c.setFillColor(Color(red / 255.0, green / 255.0, blue / 255.0))
        c.setStrokeColor(Color(red / 255.0, green / 255.0, blue / 255.0))

    category_colors = dict(DEFAULT_CATEGORY_COLORS)
    category_colors.update(style.get("category_colors", {}))
    unmapped_color = str(style.get("unmapped_motif_color", UNMAPPED_MOTIF_COLOR))

    out.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(out), pagesize=(layout["width"], layout["height"]))
    set_color(c, str(style["background_color"]))
    c.rect(0, 0, layout["width"], layout["height"], fill=1, stroke=0)

    set_color(c, "#000000")
    c.setFont("Helvetica-Bold", layout["title_font_size"])
    c.drawCentredString(layout["title_x"], pdf_y(layout["title_y"]), title)
    set_color(c, str(style["frame_color"]))
    title_line_width = float(style["title_line_width"])
    if title_line_width > 0:
        c.setLineWidth(title_line_width)
        line_half = float(style["title_line_half_length"])
        line_y = layout["title_y"] + float(style["title_line_offset"])
        c.line(layout["title_x"] - line_half, pdf_y(line_y), layout["title_x"] + line_half, pdf_y(line_y))
    set_color(c, "#000000")
    c.setFont("Helvetica", layout["subtitle_font_size"])
    c.drawCentredString(layout["title_x"], pdf_y(layout["subtitle_y"]), subtitle)

    matrix_left = layout["matrix_left"]
    matrix_top = layout["matrix_top"]
    matrix_w = layout["matrix_w"]
    matrix_h = layout["matrix_h"]
    cell_w = layout["cell_w"]
    cell_h = layout["cell_h"]

    if show_legend:
        legend_x, legend_y = determine_legend_position(
            style=style,
            width=layout["width"],
            matrix_left=matrix_left,
            matrix_w=matrix_w,
            matrix_top=matrix_top,
            matrix_h=matrix_h,
            tree_top=layout["tree_top"],
            tree_h=layout["dendro_h"],
        )
        set_color(c, "#000000")
        c.setFont("Helvetica-Bold", float(style["legend_title_font_size"]))
        c.drawString(legend_x, pdf_y(legend_y + 2), str(style["legend_title"]))
        swatch_size = float(style["legend_swatch_size"])
        row_gap = float(style["legend_row_gap"])
        for idx, category in enumerate(ordered_display_groups(ordered_rows)):
            y = legend_y + 12.0 + idx * row_gap
            set_color(c, category_colors.get(category, "#000000"))
            c.rect(legend_x, pdf_y(y + 2), swatch_size, swatch_size, fill=1, stroke=1)
            set_color(c, "#000000")
            c.setFont("Helvetica", float(style["legend_font_size"]))
            c.drawString(legend_x + swatch_size + 5, pdf_y(y + 1), category)

    for idx, row in enumerate(ordered_rows):
        short = row.get("standard_name") or row.get("protein") or row.get("y_name")
        x = matrix_left + (idx + 0.5) * cell_w
        set_color(c, group_color(row, category_colors))
        c.saveState()
        c.translate(x, pdf_y(layout["label_y"]))
        c.rotate(90)
        c.setFont("Helvetica", float(style["protein_label_font_size"]))
        c.drawString(0, 0, str(short))
        c.restoreState()

    set_color(c, str(style["grid_color"]))
    c.setLineWidth(float(style["grid_width"]))
    for idx in range(len(tree_labels) + 1):
        x = matrix_left + idx * cell_w
        c.line(x, pdf_y(matrix_top), x, pdf_y(matrix_top + matrix_h))

    dot_radius = float(style["dot_radius"])
    for row_idx, domain in enumerate(domain_keys):
        cy = matrix_top + (row_idx + 0.5) * cell_h
        meta = (domain_meta or {}).get(domain)
        label = domain_display_label(domain, domain_meta)
        set_color(c, unmapped_color if meta and meta.get("status") == "unmapped" else "#000000")
        c.setFont("Helvetica", float(style["domain_label_font_size"]))
        c.drawRightString(matrix_left - 8, pdf_y(cy + 3), label)
        for col_idx, row in enumerate(ordered_rows):
            raw = str(row.get(domain, "0") or "0").strip()
            if raw in {"0", "", "0.0"}:
                continue
            cx = matrix_left + (col_idx + 0.5) * cell_w
            if meta and meta.get("status") == "unmapped":
                fill = stroke = unmapped_color
            else:
                fill, stroke = dot_colors(row, style, category_colors)
            red, green, blue = parse_color(fill)
            c.setFillColor(Color(red / 255.0, green / 255.0, blue / 255.0))
            red, green, blue = parse_color(stroke)
            c.setStrokeColor(Color(red / 255.0, green / 255.0, blue / 255.0))
            c.setLineWidth(float(style["dot_stroke_width"]))
            c.circle(cx, pdf_y(cy), dot_radius, fill=1, stroke=1)

    set_color(c, str(style["tree_stroke_color"]))
    c.setLineWidth(float(style["tree_stroke_width"]))
    for x1, y1, x2, y2 in iter_tree_segments(tree):
        c.line(x1, pdf_y(y1), x2, pdf_y(y2))

    set_color(c, str(style["frame_color"]))
    c.setLineWidth(float(style["frame_width"]))
    c.rect(matrix_left, pdf_y(matrix_top + matrix_h), matrix_w, matrix_h, fill=0, stroke=1)
    c.showPage()
    c.save()


def export_png_from_svg(svg_path: Path, png_path: Path, dpi: float) -> bool:
    if dpi <= 0:
        raise SystemExit("--export-dpi must be > 0")

    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        svg_text = svg_path.read_text(encoding="utf-8")
        match = re.search(r'<svg[^>]*\bwidth="([0-9.]+)"[^>]*\bheight="([0-9.]+)"', svg_text)
        if match:
            width = max(1, int(round(float(match.group(1)) * dpi / SVG_DPI)))
            height = max(1, int(round(float(match.group(2)) * dpi / SVG_DPI)))
            png_path.parent.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(
                [rsvg, "--format", "png", "--width", str(width), "--height", str(height), "--output", str(png_path), str(svg_path)],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0 and png_path.exists():
                try:
                    Image, _ImageDraw, _ImageFont = load_pillow()
                    image = Image.open(png_path)
                    image.save(png_path, dpi=(dpi, dpi))
                except SystemExit:
                    pass
                return True

    sips = shutil.which("sips")
    if not sips:
        return False

    Image, _ImageDraw, _ImageFont = load_pillow()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_png = Path(tmpdir) / "figure.png"
        result = subprocess.run(
            [sips, "-s", "format", "png", str(svg_path), "--out", str(tmp_png)],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0 or not tmp_png.exists():
            return False

        rendered = Image.open(tmp_png).convert("RGBA")
        white = Image.new("RGBA", rendered.size, (255, 255, 255, 255))
        white.alpha_composite(rendered)
        image = white.convert("RGB")
        if dpi != SVG_DPI:
            scale = dpi / SVG_DPI
            width = max(1, int(round(image.width * scale)))
            height = max(1, int(round(image.height * scale)))
            image = image.resize((width, height), Image.Resampling.LANCZOS)
        image.save(png_path, dpi=(dpi, dpi))
    return True


def export_pdf_from_svg(svg_path: Path, pdf_path: Path) -> bool:
    # rsvg-convert writes a vector PDF that matches the SVG exactly. sips is not
    # used here because it rasterizes SVG input.
    rsvg = shutil.which("rsvg-convert")
    if not rsvg:
        return False
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [rsvg, "--format", "pdf", "--output", str(pdf_path), str(svg_path)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.returncode == 0 and pdf_path.exists()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE)
    parser.add_argument("--hits", type=Path, default=DEFAULT_HITS)
    parser.add_argument("--proteins", type=Path, default=DEFAULT_PROTEINS)
    parser.add_argument("--ptm-intersections", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--show-ptm-overlay", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--order-csv", type=Path, default=DEFAULT_ORDER_CSV)
    parser.add_argument("--png", type=Path, help="Optional high-resolution PNG export path. Requires Pillow.")
    parser.add_argument("--pdf", type=Path, help="Optional vector PDF export path. Requires ReportLab.")
    parser.add_argument("--export-dpi", type=float, default=600.0, help="DPI for PNG export")

    parser.add_argument("--title", help="Default: 'All proteins', or '<category> proteins' with --category.")
    parser.add_argument(
        "--subtitle",
        help="Default describes the column order, e.g. 'Ordered by Clustal Omega sequence similarity (neighbour-joining)'.",
    )
    parser.add_argument(
        "--order-by",
        choices=["tree", "motifs"],
        default="tree",
        help=(
            "Column order. tree: leaf order of --tree (Clustal Omega). motifs: average-linkage clustering on the "
            "fraction of plotted motifs two proteins share; --tree is not used."
        ),
    )
    parser.add_argument("--category", choices=["Pol I", "Pol II", "Pol III"])
    parser.add_argument("--min-proteins", type=int, default=2)
    parser.add_argument(
        "--dedupe-harmonized",
        action="store_true",
        help="Collapse mapped motif/domain hits to InterPro rows; keep unmapped source motifs as grey rows.",
    )
    parser.add_argument(
        "--hide-empty-columns",
        action="store_true",
        help="Remove protein columns that have no plotted motif/domain hits after filtering.",
    )

    parser.add_argument(
        "--exclude-mobidblite",
        dest="exclude_mobidblite",
        action="store_true",
        default=True,
        help="Exclude MobiDBLite domains (default)",
    )
    parser.add_argument(
        "--include-mobidblite",
        dest="exclude_mobidblite",
        action="store_false",
        help="Include MobiDBLite domains",
    )

    parser.add_argument("--style-config", type=Path, help="Optional JSON file with style/layout overrides")
    parser.add_argument("--no-legend", action="store_true")

    parser.add_argument("--font-family")
    parser.add_argument("--margin-left", type=float)
    parser.add_argument("--margin-right", type=float)
    parser.add_argument("--margin-top", type=float)
    parser.add_argument("--margin-bottom", type=float)
    parser.add_argument("--cell-width", type=float)
    parser.add_argument("--cell-height", type=float)
    parser.add_argument("--column-label-height", type=float)
    parser.add_argument("--dendrogram-gap", type=float)
    parser.add_argument("--dendrogram-height", type=float)

    parser.add_argument("--title-font-size", type=float)
    parser.add_argument("--subtitle-font-size", type=float)
    parser.add_argument("--protein-label-font-size", type=float)
    parser.add_argument("--domain-label-font-size", type=float)
    parser.add_argument("--title-subtitle-gap", type=float)
    parser.add_argument("--subtitle-matrix-gap", type=float)

    parser.add_argument("--grid-color")
    parser.add_argument("--dot-fill")
    parser.add_argument("--dot-stroke")
    parser.add_argument("--unmapped-motif-color")
    parser.add_argument("--background-color")
    parser.add_argument("--category-color", action="append", default=[], help="Repeatable: 'Pol I=#1f77b4'")

    parser.add_argument(
        "--legend-position",
        choices=["auto", "top-left", "top-right", "bottom-left", "bottom-right", "custom"],
    )
    parser.add_argument("--legend-x", type=float)
    parser.add_argument("--legend-y", type=float)

    parser.add_argument("--write-default-style-config", type=Path, help="Write default style JSON and exit")
    args = parser.parse_args()
    if args.title is None:
        args.title = f"{args.category} proteins" if args.category else "All proteins"
    if args.subtitle is None:
        args.subtitle = DEFAULT_SUBTITLES[args.order_by]
    return args


def build_style(args: argparse.Namespace) -> dict:
    style = merge_style(DEFAULT_STYLE, {})
    if args.style_config:
        style = merge_style(style, load_style_config(args.style_config))

    mapping = {
        "font_family": args.font_family,
        "margin_left": args.margin_left,
        "margin_right": args.margin_right,
        "margin_top": args.margin_top,
        "margin_bottom": args.margin_bottom,
        "cell_width": args.cell_width,
        "cell_height": args.cell_height,
        "column_label_height": args.column_label_height,
        "dendrogram_gap": args.dendrogram_gap,
        "dendrogram_height": args.dendrogram_height,
        "title_font_size": args.title_font_size,
        "subtitle_font_size": args.subtitle_font_size,
        "protein_label_font_size": args.protein_label_font_size,
        "domain_label_font_size": args.domain_label_font_size,
        "title_subtitle_gap": args.title_subtitle_gap,
        "subtitle_matrix_gap": args.subtitle_matrix_gap,
        "grid_color": args.grid_color,
        "dot_fill": args.dot_fill,
        "dot_stroke": args.dot_stroke,
        "unmapped_motif_color": args.unmapped_motif_color,
        "background_color": args.background_color,
        "legend_position": args.legend_position,
        "legend_x": args.legend_x,
        "legend_y": args.legend_y,
    }
    for key, value in mapping.items():
        if value is not None:
            style[key] = value

    if args.category_color:
        parsed = parse_category_color_overrides(args.category_color)
        colors = dict(style.get("category_colors", {}))
        colors.update(parsed)
        style["category_colors"] = colors

    return style


def validate_style(style: dict) -> None:
    nonnegative_keys = [
        "margin_top",
    ]
    positive_keys = [
        "margin_left",
        "margin_right",
        "margin_bottom",
        "cell_width",
        "cell_height",
        "column_label_height",
        "dendrogram_height",
        "title_font_size",
        "subtitle_font_size",
        "protein_label_font_size",
        "domain_label_font_size",
    ]
    for key in nonnegative_keys:
        value = float(style[key])
        if value < 0:
            raise SystemExit(f"Style value {key} must be >= 0, got {value}")
    for key in positive_keys:
        value = float(style[key])
        if value <= 0:
            raise SystemExit(f"Style value {key} must be > 0, got {value}")


def main() -> int:
    args = parse_args()

    if args.write_default_style_config:
        args.write_default_style_config.parent.mkdir(parents=True, exist_ok=True)
        args.write_default_style_config.write_text(json.dumps(DEFAULT_STYLE, indent=2), encoding="utf-8")
        print(f"Wrote default style config: {args.write_default_style_config}")
        return 0

    if args.min_proteins < 1:
        raise SystemExit("--min-proteins must be >= 1")

    style = build_style(args)
    validate_style(style)

    data = prepare_plot_data(args)
    pruned_tree = data["tree"]
    resolved_tree_labels = data["labels"]
    ordered_rows = data["rows"]
    domain_keys = data["domain_keys"]
    domain_meta = data["domain_meta"]
    missing_from_tree = data["missing_from_tree"]

    write_order_csv(args.order_csv, ordered_rows, resolved_tree_labels)
    build_svg(
        tree=pruned_tree,
        tree_labels=resolved_tree_labels,
        domain_keys=domain_keys,
        ordered_rows=ordered_rows,
        out=args.out,
        title=args.title,
        subtitle=args.subtitle,
        show_legend=not args.no_legend,
        style=style,
        domain_meta=domain_meta,
    )
    if args.png:
        exported = export_png_from_svg(args.out, args.png, args.export_dpi)
        if not exported:
            build_png(
                tree=pruned_tree,
                tree_labels=resolved_tree_labels,
                domain_keys=domain_keys,
                ordered_rows=ordered_rows,
                out=args.png,
                title=args.title,
                subtitle=args.subtitle,
                show_legend=not args.no_legend,
                style=style,
                dpi=args.export_dpi,
                domain_meta=domain_meta,
            )
    if args.pdf and not export_pdf_from_svg(args.out, args.pdf):
        build_pdf(
            tree=pruned_tree,
            tree_labels=resolved_tree_labels,
            domain_keys=domain_keys,
            ordered_rows=ordered_rows,
            out=args.pdf,
            title=args.title,
            subtitle=args.subtitle,
            show_legend=not args.no_legend,
            style=style,
            domain_meta=domain_meta,
        )
    manifest = write_plot_manifest(args, "native", ordered_rows, domain_keys, missing_from_tree)

    print(f"Wrote figure: {args.out}")
    if args.png:
        print(f"Wrote PNG: {args.png}")
    if args.pdf:
        print(f"Wrote PDF: {args.pdf}")
    print(f"Wrote order table: {args.order_csv}")
    print(f"Wrote figure manifest: {manifest}")
    print(f"Proteins plotted: {len(ordered_rows)}")
    print(f"Shared domains plotted: {len(domain_keys)}")
    print(f"Excluded MobiDBLite: {args.exclude_mobidblite}")
    print(f"Min proteins threshold: {args.min_proteins}")
    print(f"Deduped harmonized motifs: {args.dedupe_harmonized}")
    print(f"Hidden empty protein columns: {args.hide_empty_columns}")
    print(f"Column order: {args.order_by} ({args.subtitle})")
    if missing_from_tree:
        print("Missing from tree:", ", ".join(missing_from_tree))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
