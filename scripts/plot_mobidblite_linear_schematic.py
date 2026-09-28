#!/usr/bin/env python3
"""Render separate linear protein schematics for domains and MobiDBLite disorder.

This script intentionally produces its own figure and tables. It does not modify
or feed the motif/domain matrix plot.
"""

from __future__ import annotations

import argparse
import colorsys
import csv
import html
from collections import defaultdict
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROTEINS = REPO_ROOT / "data" / "output" / "polymerase_tf_proteins.csv"
DEFAULT_HITS = REPO_ROOT / "data" / "output" / "polymerase_tf_domain_hits_raw.csv"
DEFAULT_PTMS = REPO_ROOT / "data" / "output" / "polymerase_tf_ptm_sites.csv"
DEFAULT_OUT = REPO_ROOT / "data" / "output" / "all48_mobidblite_linear_schematic.svg"
DEFAULT_REGIONS_OUT = REPO_ROOT / "data" / "output" / "all48_mobidblite_linear_regions.csv"
DEFAULT_SUMMARY_OUT = REPO_ROOT / "data" / "output" / "all48_mobidblite_linear_summary.csv"
DEFAULT_MARKDOWN_OUT = REPO_ROOT / "data" / "output" / "all48_mobidblite_linear_summary.md"

# Domain/motif blocks use identity colors, not source colors. The exact label is
# available through SVG hover tooltips, so the color only helps distinguish nearby
# features visually.
DISORDER_COLOR = "#5B5B5B"
PTM_COLORS = {
    "phosphorylation": "#D62728",
    "sumoylation": "#9467BD",
    "ubiquitination": "#2CA02C",
    "acetylation": "#17BECF",
    "methylation": "#E377C2",
    "succinylation": "#8C564B",
    "other": "#7F7F7F",
}


def clean_int(value: str | int | None) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def protein_length(row: dict) -> int:
    explicit = clean_int(row.get("length") or row.get("protein_length"))
    if explicit:
        return explicit
    sequence = (row.get("amino_acid_sequence") or row.get("sequence") or "").strip()
    return len(sequence.replace(" ", "").replace("\n", ""))


def short_label(row: dict) -> str:
    raw_key = (row.get("raw_domain_key") or "").strip()
    if raw_key and "::" in raw_key:
        return raw_key.split("::", 1)[1]
    if raw_key:
        return raw_key
    return (row.get("domain_description") or row.get("source") or "region").strip()


def is_mobidblite(row: dict) -> bool:
    source = (row.get("source") or "").strip()
    raw_key = (row.get("raw_domain_key") or "").strip()
    return source == "MobiDBLite" or "MobiDBLite" in raw_key


def read_proteins(path: Path, category: str | None = None) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    out = []
    for row in rows:
        if category and row.get("category") != category:
            continue
        length = protein_length(row)
        if length <= 0:
            continue
        copied = dict(row)
        copied["protein_length"] = length
        out.append(copied)
    return out


def build_regions(hits_path: Path, proteins: Iterable[dict]) -> list[dict]:
    protein_by_y = {row["y_name"]: row for row in proteins}
    regions = []
    with hits_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            y_name = row.get("y_name", "")
            protein = protein_by_y.get(y_name)
            if not protein:
                continue

            start = clean_int(row.get("start"))
            end = clean_int(row.get("end"))
            if start is None or end is None:
                continue
            if start > end:
                start, end = end, start

            length = int(protein["protein_length"])
            if start > length:
                # Coordinates from a different isoform/sequence version; clamping would invent a 1-aa region.
                continue
            start = max(1, min(start, length))
            end = max(1, min(end, length))
            if end < start:
                continue

            feature_type = "disorder" if is_mobidblite(row) else "domain"
            region_len = end - start + 1
            regions.append(
                {
                    "category": protein.get("category", ""),
                    "protein": protein.get("protein", ""),
                    "standard_name": protein.get("standard_name", ""),
                    "y_name": y_name,
                    "accession": protein.get("accession") or protein.get("uniprot_id", ""),
                    "protein_length": length,
                    "feature_type": feature_type,
                    "source": row.get("source", ""),
                    "raw_domain_key": row.get("raw_domain_key", ""),
                    "harmonized_key": row.get("harmonized_key", ""),
                    "label": "MobiDBLite disorder" if feature_type == "disorder" else short_label(row),
                    "description": row.get("domain_description", ""),
                    "start": start,
                    "end": end,
                    "region_length": region_len,
                    "fraction_of_protein": round(region_len / length, 5),
                }
            )
    regions.sort(key=lambda r: (r["category"], r["protein"], r["feature_type"] != "domain", r["start"], r["end"], r["label"]))
    return regions


def merged_coverage(intervals: Iterable[tuple[int, int]]) -> int:
    sorted_intervals = sorted((min(a, b), max(a, b)) for a, b in intervals if a and b)
    if not sorted_intervals:
        return 0
    total = 0
    current_start, current_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    total += current_end - current_start + 1
    return total


def summarize_regions(proteins: list[dict], regions: list[dict]) -> list[dict]:
    by_y: dict[str, list[dict]] = defaultdict(list)
    for region in regions:
        by_y[region["y_name"]].append(region)

    summaries = []
    for protein in proteins:
        y_name = protein["y_name"]
        length = int(protein["protein_length"])
        protein_regions = by_y.get(y_name, [])
        domain_regions = [r for r in protein_regions if r["feature_type"] == "domain"]
        disorder_regions = [r for r in protein_regions if r["feature_type"] == "disorder"]
        domain_coverage = merged_coverage((int(r["start"]), int(r["end"])) for r in domain_regions)
        disorder_coverage = merged_coverage((int(r["start"]), int(r["end"])) for r in disorder_regions)
        summaries.append(
            {
                "category": protein.get("category", ""),
                "protein": protein.get("protein", ""),
                "standard_name": protein.get("standard_name", ""),
                "y_name": y_name,
                "accession": protein.get("accession") or protein.get("uniprot_id", ""),
                "protein_length": length,
                "domain_region_count": len(domain_regions),
                "disorder_region_count": len(disorder_regions),
                "domain_coverage_aa": domain_coverage,
                "domain_coverage_fraction": round(domain_coverage / length, 5),
                "disorder_coverage_aa": disorder_coverage,
                "disorder_coverage_fraction": round(disorder_coverage / length, 5),
            }
        )
    return summaries


def assign_tracks(regions: list[dict]) -> dict[int, int]:
    tracks_end: list[int] = []
    assignment: dict[int, int] = {}
    for idx, region in sorted(enumerate(regions), key=lambda item: (int(item[1]["start"]), int(item[1]["end"]))):
        start = int(region["start"])
        end = int(region["end"])
        for track, last_end in enumerate(tracks_end):
            if start > last_end + 4:
                tracks_end[track] = end
                assignment[idx] = track
                break
        else:
            assignment[idx] = len(tracks_end)
            tracks_end.append(end)
    return assignment


def stable_hash(text: str) -> int:
    value = 0
    for ch in text:
        value = (value * 131 + ord(ch)) % 1_000_003
    return value


def domain_color(label: str, raw_key: str) -> str:
    identity = label or raw_key or "domain"
    hue = (stable_hash(identity) % 360) / 360.0
    red, green, blue = colorsys.hls_to_rgb(hue, 0.54, 0.66)
    return f"#{int(red * 255):02X}{int(green * 255):02X}{int(blue * 255):02X}"


def normalize_ptm_family(value: str | None) -> str:
    text = (value or "").strip().lower()
    return text if text in PTM_COLORS else "other"


def read_ptms(path: Path | None, proteins: Iterable[dict]) -> list[dict]:
    if path is None or not path.is_file():
        return []
    protein_by_y = {row["y_name"]: row for row in proteins}
    ptms = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            y_name = row.get("y_name", "")
            protein = protein_by_y.get(y_name)
            if not protein:
                continue
            site = clean_int(row.get("site_index"))
            if site is None:
                continue
            length = int(protein["protein_length"])
            if site < 1 or site > length:
                continue
            family = normalize_ptm_family(row.get("ptm_family"))
            copied = dict(row)
            copied["site_index"] = site
            copied["ptm_family"] = family
            ptms.append(copied)
    ptms.sort(key=lambda row: (row["y_name"], row["site_index"], row["ptm_family"]))
    return ptms


def render_svg(proteins: list[dict], regions: list[dict], ptms: list[dict], out: Path, title: str) -> None:
    by_y: dict[str, list[dict]] = defaultdict(list)
    for region in regions:
        by_y[region["y_name"]].append(region)
    ptms_by_y: dict[str, list[dict]] = defaultdict(list)
    for ptm in ptms:
        ptms_by_y[ptm["y_name"]].append(ptm)

    row_layout = []
    current_y = 104
    for protein in proteins:
        protein_regions = by_y.get(protein["y_name"], [])
        domains = [r for r in protein_regions if r["feature_type"] == "domain"]
        disorders = [r for r in protein_regions if r["feature_type"] == "disorder"]
        protein_ptms = ptms_by_y.get(protein["y_name"], [])
        domain_tracks = assign_tracks(domains)
        disorder_tracks = assign_tracks(disorders)
        domain_track_count = max(domain_tracks.values(), default=-1) + 1
        disorder_track_count = max(disorder_tracks.values(), default=-1) + 1
        top_pad = max(24, 17 + domain_track_count * 8)
        bottom_pad = max(32, 22 + disorder_track_count * 8)
        backbone_y = current_y + top_pad
        row_layout.append((protein, backbone_y, domains, disorders, protein_ptms, domain_tracks, disorder_tracks))
        current_y = backbone_y + bottom_pad

    max_len = max(int(row["protein_length"]) for row in proteins)
    left = 175
    right = 42
    top = 104
    axis_w = 820
    width = left + axis_w + right
    height = current_y + 64

    def sx(position: int) -> float:
        return left + ((position - 1) / max(max_len - 1, 1)) * axis_w

    tick_step = 250 if max_len > 1000 else 100
    ticks = list(range(0, max_len + tick_step, tick_step))
    ticks[0] = 1

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.small{font-size:10px}.label{font-size:11px}.title{font-size:18px;font-weight:700}.muted{fill:#5f6368}.axis{stroke:#9aa0a6;stroke-width:1}.domain{stroke:#ffffff;stroke-width:.6}.disorder{stroke:#7a3900;stroke-width:.4;opacity:.72}",
        "</style>",
        f'<text x="{left}" y="30" class="title">{html.escape(title)}</text>',
        f'<text x="{left}" y="50" class="small muted">Domains above each protein; MobiDBLite disorder below; known SGD PTMs on the backbone. Hover in HTML for exact labels.</text>',
    ]

    axis_y = top - 18
    parts.append(f'<line x1="{left}" y1="{axis_y}" x2="{left + axis_w}" y2="{axis_y}" class="axis"/>')
    for tick in ticks:
        position = min(max(tick, 1), max_len)
        x = sx(position)
        parts.append(f'<line x1="{x:.2f}" y1="{axis_y - 4}" x2="{x:.2f}" y2="{axis_y + 4}" class="axis"/>')
        parts.append(f'<text x="{x:.2f}" y="{axis_y - 8}" class="small muted" text-anchor="middle">{position}</text>')

    for protein, y, domains, disorders, protein_ptms, domain_tracks, disorder_tracks in row_layout:
        name = protein.get("protein") or protein.get("standard_name") or protein.get("y_name")
        category = protein.get("category", "")
        length = int(protein["protein_length"])
        y_name = protein["y_name"]
        parts.append(f'<text x="18" y="{y + 4}" class="label">{html.escape(name)}</text>')
        parts.append(f'<text x="18" y="{y + 18}" class="small muted">{html.escape(category)} | {html.escape(y_name)} | {length} aa</text>')
        parts.append(f'<line x1="{left}" y1="{y}" x2="{sx(length):.2f}" y2="{y}" stroke="#3c4043" stroke-width="2.2" stroke-linecap="round"/>')

        for ptm in protein_ptms:
            x = sx(int(ptm["site_index"]))
            family = normalize_ptm_family(ptm.get("ptm_family"))
            color = PTM_COLORS[family]
            residue = ptm.get("site_residue", "")
            protein_name = protein.get("protein") or protein.get("standard_name") or protein.get("y_name")
            tooltip = f'{protein_name} | {family}: {residue}{ptm["site_index"]}' if residue else f'{protein_name} | {family}: {ptm["site_index"]}'
            parts.append(f'<line x1="{x:.2f}" y1="{y - 7:.2f}" x2="{x:.2f}" y2="{y + 7:.2f}" stroke="{color}" stroke-width="1.25"><title>{html.escape(tooltip)}</title></line>')
            parts.append(f'<circle cx="{x:.2f}" cy="{y - 9:.2f}" r="2.3" fill="{color}" stroke="#ffffff" stroke-width="0.5"><title>{html.escape(tooltip)}</title></circle>')

        for idx, region in enumerate(domains):
            x1 = sx(int(region["start"]))
            x2 = sx(int(region["end"]))
            track = domain_tracks[idx]
            bar_y = y - 17 - track * 8
            color = domain_color(region.get("label", ""), region["raw_domain_key"])
            source = region.get("source", "domain")
            protein_name = protein.get("protein") or protein.get("standard_name") or protein.get("y_name")
            tooltip = f'{protein_name} | {region["label"]} ({source}): {region["start"]}-{region["end"]} aa'
            parts.append(f'<rect class="domain" x="{x1:.2f}" y="{bar_y:.2f}" width="{max(x2 - x1, 2):.2f}" height="7" rx="1.5" fill="{color}"><title>{html.escape(tooltip)}</title></rect>')

        for idx, region in enumerate(disorders):
            x1 = sx(int(region["start"]))
            x2 = sx(int(region["end"]))
            track = disorder_tracks[idx]
            bar_y = y + 6 + track * 8
            protein_name = protein.get("protein") or protein.get("standard_name") or protein.get("y_name")
            tooltip = f'{protein_name} | MobiDBLite disorder: {region["start"]}-{region["end"]} aa'
            parts.append(f'<rect class="disorder" x="{x1:.2f}" y="{bar_y:.2f}" width="{max(x2 - x1, 2):.2f}" height="7" rx="1.5" fill="{DISORDER_COLOR}"><title>{html.escape(tooltip)}</title></rect>')

    legend_y = height - 42
    parts.append(f'<text x="{left}" y="{legend_y - 28}" class="small muted">Color legend: each rainbow domain block color represents a distinct motif/domain label, so repeated colors mark the same annotation type across proteins.</text>')
    parts.append(f'<text x="{left}" y="{legend_y - 15}" class="small muted">Gray blocks mark MobiDBLite-predicted disorder below the backbone; colored circles/vertical ticks mark known SGD PTM residue sites. Hover in HTML for exact labels and coordinates.</text>')
    legend_items = [
        ("domain / motif identity", domain_color("example domain", "example domain")),
        ("MobiDBLite disorder", DISORDER_COLOR),
        ("phosphorylation", PTM_COLORS["phosphorylation"]),
        ("sumoylation", PTM_COLORS["sumoylation"]),
        ("ubiquitination", PTM_COLORS["ubiquitination"]),
        ("methylation", PTM_COLORS["methylation"]),
        ("acetylation/other", PTM_COLORS["acetylation"]),
    ]
    x_cursor = left
    y_cursor = legend_y
    for label, color in legend_items:
        if x_cursor > left + 650:
            x_cursor = left
            y_cursor += 16
        if label in PTM_COLORS or label == "acetylation/other":
            parts.append(f'<circle cx="{x_cursor + 6}" cy="{y_cursor - 4}" r="3" fill="{color}"/>')
        else:
            parts.append(f'<rect x="{x_cursor}" y="{y_cursor - 9}" width="16" height="7" fill="{color}" class="domain"/>')
        parts.append(f'<text x="{x_cursor + 22}" y="{y_cursor}" class="small">{html.escape(label)}</text>')
        x_cursor += max(92, len(label) * 6 + 34)
    parts.append("</svg>")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, summaries: list[dict], regions: list[dict], title: str, ptm_count: int = 0) -> None:
    domain_count = sum(1 for row in regions if row["feature_type"] == "domain")
    disorder_count = sum(1 for row in regions if row["feature_type"] == "disorder")
    disorder_proteins = sum(1 for row in summaries if int(row["disorder_region_count"]) > 0)
    top_disorder = sorted(summaries, key=lambda row: float(row["disorder_coverage_fraction"]), reverse=True)[:10]
    lines = [
        f"# {title}",
        "",
        "This is a separate MobiDBLite-style linear protein schematic output. It is not part of the motif/domain matrix plot.",
        "",
        "## Totals",
        "",
        f"- Proteins: {len(summaries)}",
        f"- Domain regions: {domain_count}",
        f"- MobiDBLite disorder regions: {disorder_count}",
        f"- Proteins with MobiDBLite disorder: {disorder_proteins}",
        "",
        "## Highest MobiDBLite disorder coverage",
        "",
    ]
    for row in top_disorder:
        percent = 100 * float(row["disorder_coverage_fraction"])
        lines.append(f"- {row['protein']} ({row['y_name']}): {row['disorder_coverage_aa']} aa ({percent:.1f}%)")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Domains and MobiDBLite disorder are drawn on a protein-length amino-acid axis.",
            (
                f"- {ptm_count} known SGD PTM sites are marked on the protein backbone. PTMs are never added to the matrix plots."
                if ptm_count
                else "- No PTM sites were drawn on this schematic. PTMs are never added to the matrix plots."
            ),
            "- Coverage values merge overlapping intervals within each feature class before calculating amino-acid coverage.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, default=DEFAULT_PROTEINS)
    parser.add_argument("--hits", type=Path, default=DEFAULT_HITS)
    parser.add_argument(
        "--ptm-sites",
        type=Path,
        help=(
            "Known PTM site table to mark on the backbone. Default: polymerase_tf_ptm_sites.csv next to --proteins "
            "(same run). Pass --ptm-sites \"\" to omit PTMs."
        ),
    )
    parser.add_argument("--category", help="Optional category filter, e.g. 'Pol II'")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="SVG output path")
    parser.add_argument("--regions-csv", type=Path, default=DEFAULT_REGIONS_OUT)
    parser.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY_OUT)
    parser.add_argument("--summary-md", type=Path, default=DEFAULT_MARKDOWN_OUT)
    parser.add_argument("--title", default="All 48 Protein TFs: Domains and MobiDBLite Disorder")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    proteins = read_proteins(args.proteins, category=args.category)
    if not proteins:
        raise SystemExit("No proteins with sequence lengths were found.")

    regions = build_regions(args.hits, proteins)
    ptm_path = args.ptm_sites if args.ptm_sites is not None else args.proteins.parent / DEFAULT_PTMS.name
    ptms = read_ptms(ptm_path, proteins)
    summaries = summarize_regions(proteins, regions)

    region_fields = [
        "category",
        "protein",
        "standard_name",
        "y_name",
        "accession",
        "protein_length",
        "feature_type",
        "source",
        "raw_domain_key",
        "harmonized_key",
        "label",
        "description",
        "start",
        "end",
        "region_length",
        "fraction_of_protein",
    ]
    summary_fields = [
        "category",
        "protein",
        "standard_name",
        "y_name",
        "accession",
        "protein_length",
        "domain_region_count",
        "disorder_region_count",
        "domain_coverage_aa",
        "domain_coverage_fraction",
        "disorder_coverage_aa",
        "disorder_coverage_fraction",
    ]

    render_svg(proteins, regions, ptms, args.out, args.title)
    write_csv(args.regions_csv, regions, region_fields)
    write_csv(args.summary_csv, summaries, summary_fields)
    write_markdown(args.summary_md, summaries, regions, args.title, ptm_count=len(ptms))

    print(f"Proteins plotted: {len(proteins)}")
    print(f"Regions written: {len(regions)}")
    print(f"PTM markers drawn: {len(ptms)}")
    print(f"SVG: {args.out}")
    print(f"Region table: {args.regions_csv}")
    print(f"Summary table: {args.summary_csv}")
    print(f"Summary markdown: {args.summary_md}")


if __name__ == "__main__":
    main()
