#!/usr/bin/env python3
"""Build a standalone HTML index of the figures and tables in a pipeline output folder.

The index is a single HTML file (no server, no internet, no extra packages)
written into the output folder. Open it by double-clicking it on macOS or
Windows. It groups motif dot-matrix figures into All proteins / Pol I / Pol II /
Pol III sections, links every SVG/PNG/PDF/CSV next to each preview, and shows
each file's full path on the viewer's own machine with a copy button.

Links are relative, so the folder can be moved, zipped, or shared between
macOS and Windows without breaking them. Use --embed-previews to inline the
preview images as well, so the page still shows the figures when the HTML file
is sent on its own (the file links then need the folder next to it).

Figures written by plot_clustalo_motif_figure*.py carry a ``<name>.plot.json``
sidecar that records their category, options and counts. Figures without one
(older outputs) are classified from their file names.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from fnmatch import fnmatch
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import quote

FIGURE_EXTS = (".svg", ".png", ".pdf")
PREVIEW_EXTS = (".png", ".svg")
EMBED_EXTS = (".svg", ".png")
TABLE_EXTS = {".csv", ".tsv", ".md", ".json", ".xlsx", ".txt", ".html", ".zip", ".ndjson"}
SKIP_DIRS = {"cache", "__pycache__", "mplconfig"}
SKIP_FILES = {".ds_store", "thumbs.db", "desktop.ini"}
MAX_CSV_BYTES_TO_COUNT = 50 * 1024 * 1024

SECTION_ORDER = ["All", "Pol I", "Pol II", "Pol III", "Other"]
SECTION_TITLES = {
    "All": "All proteins",
    "Pol I": "Pol I",
    "Pol II": "Pol II",
    "Pol III": "Pol III",
    "Other": "Other figures",
}
GROUP_COLORS = {
    "Pol I": "#1f77b4",
    "Pol I axial": "#084594",
    "Pol I peri-axial": "#6baed6",
    "Pol II": "#d62728",
    "Pol III": "#2ca02c",
}
CATEGORY_PATTERN = re.compile(r"(?:^|[^a-z0-9])pol[_\-. ]?(iii|ii|i|3|2|1)(?=$|[^a-z0-9])")
CATEGORY_TOKENS = {"i": "Pol I", "1": "Pol I", "ii": "Pol II", "2": "Pol II", "iii": "Pol III", "3": "Pol III"}

# (glob pattern, description). First match wins, so specific patterns go first.
TABLE_DESCRIPTIONS = [
    ("run_summary.json", "Extraction run summary: input file, category counts and domain totals."),
    ("polymerase_tf_proteins.csv", "Proteins looked up in SGD, with IDs, group labels and sequences."),
    ("polymerase_tf_domain_hits_raw.csv", "One row per domain/motif hit, with its InterPro mapping. Input for the plots."),
    ("polymerase_tf_domain_hits.csv", "Same rows as polymerase_tf_domain_hits_raw.csv (kept for older scripts)."),
    ("polymerase_tf_domain_hits_harmonized.csv", "Hits collapsed to harmonized (InterPro) domain keys."),
    ("polymerase_tf_domain_source_summary.csv", "Hit counts and InterPro mapping rates per source database."),
    ("polymerase_tf_domain_matrix_harmonized_counts_all.csv", "Protein x harmonized-domain hit counts, all domains."),
    ("polymerase_tf_domain_matrix_harmonized_counts_min*.csv", "Protein x harmonized-domain hit counts, shared domains only."),
    ("polymerase_tf_domain_matrix_harmonized_all.csv", "Protein x harmonized-domain presence (0/1), all domains."),
    ("polymerase_tf_domain_matrix_harmonized_min*.csv", "Protein x harmonized-domain presence (0/1), shared domains only."),
    ("polymerase_tf_domain_matrix_counts_all.csv", "Protein x raw-domain hit counts, all domains."),
    ("polymerase_tf_domain_matrix_counts_min*.csv", "Protein x raw-domain hit counts, shared domains only."),
    ("polymerase_tf_domain_matrix_all.csv", "Protein x raw-domain presence (0/1), all domains."),
    ("polymerase_tf_domain_matrix_min*.csv", "Protein x raw-domain presence (0/1), shared domains only."),
    ("polymerase_tf_ptm_sites_raw.csv", "SGD known PTM evidence records (one row per record)."),
    ("polymerase_tf_ptm_sites.csv", "SGD known PTM sites, one row per protein, position and PTM type."),
    ("polymerase_tf_ptm_domain_intersections.csv", "Known PTM sites classified as inside, near, or outside domains."),
    ("polymerase_tf_ptm_summary_by_protein.csv", "Known PTM counts per protein, by PTM family."),
    ("polymerase_tf_ptm_predictions_mtprompt.csv", "Optional MTPrompt-PTM predictions (hypotheses, not known sites)."),
    ("polymerase_tf_ptm_candidate_sites.csv", "Known and predicted PTM candidate sites with evidence tiers."),
    ("sgd_failures.csv", "SGD lookups that failed. Should be empty."),
    ("interpro_failures.csv", "InterPro lookups that failed. Should be empty."),
    ("*shared_domain_category_summary.csv", "Shared raw domains with per-category protein counts."),
    ("*mobidblite_linear_regions.csv", "Domain and MobiDBLite disorder intervals used in the linear schematic."),
    ("*mobidblite_linear_summary.csv", "Per-protein length, domain coverage and disorder coverage."),
    ("*mobidblite_linear_summary.md", "Text summary of the linear schematic."),
    ("*_order.csv", "Column order for a figure."),
    ("*_columns.csv", "Column order for a figure."),
    ("*.xlsx", "Excel workbook."),
    ("*.html", "HTML report."),
    ("*.md", "Markdown notes."),
    ("*.zip", "Zip archive."),
]
FAILURE_TABLES = {"sgd_failures.csv", "interpro_failures.csv"}


@dataclass
class Figure:
    stem: str
    files: Dict[str, Path] = field(default_factory=dict)
    manifest: Optional[dict] = None
    order_csv: Optional[Path] = None
    section: str = "Other"
    kind: str = "Figure"
    tags: List[str] = field(default_factory=list)
    title: str = ""
    subtitle: str = ""
    caption: str = ""
    counts: Dict[str, object] = field(default_factory=dict)
    groups: Dict[str, int] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


def rel_href(path: Path, html_dir: Path) -> str:
    """URL-encoded link from the HTML file to ``path``; forward slashes work on every OS."""
    try:
        rel = Path(os.path.relpath(path.resolve(), html_dir.resolve())).as_posix()
    except ValueError:  # Windows: different drive letters
        return path.resolve().as_uri()
    return "/".join(quote(part) for part in rel.split("/"))


def rel_display(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{num_bytes} B"


def modified(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")


def read_json(path: Path) -> Optional[dict]:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def csv_shape(path: Path) -> Optional[tuple]:
    if path.stat().st_size > MAX_CSV_BYTES_TO_COUNT:
        return None
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    try:
        with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, None)
            if header is None:
                return (0, 0)
            rows = sum(1 for row in reader if any(cell.strip() for cell in row))
    except (OSError, csv.Error):
        return None
    return (rows, len(header))


def describe_table(name: str) -> str:
    lowered = name.lower()
    for pattern, description in TABLE_DESCRIPTIONS:
        if fnmatch(lowered, pattern):
            return description
    return ""


def category_from_name(stem: str) -> Optional[str]:
    match = CATEGORY_PATTERN.search(stem.lower())
    return CATEGORY_TOKENS[match.group(1)] if match else None


def group_counts_from_order_csv(path: Path) -> Dict[str, int]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return {}
    counts: Counter = Counter()
    for row in rows:
        group = (row.get("display_group") or row.get("category") or "").strip()
        if group:
            counts[group] += 1
    return dict(counts)


def collect_files(outdir: Path, recursive: bool, html_out: Path) -> tuple:
    files: List[Path] = []
    subfolders: List[Path] = []
    for entry in sorted(outdir.iterdir(), key=lambda p: p.name.lower()):
        if entry.name.startswith(".") or entry.name.lower() in SKIP_FILES:
            continue
        if entry.is_dir():
            subfolders.append(entry)
            if recursive and entry.name.lower() not in SKIP_DIRS:
                for child in sorted(entry.rglob("*"), key=lambda p: p.as_posix().lower()):
                    parts = child.relative_to(outdir).parts
                    if any(part.startswith(".") or part.lower() in SKIP_DIRS for part in parts):
                        continue
                    if child.is_file() and child.name.lower() not in SKIP_FILES:
                        files.append(child)
            continue
        if entry.is_file():
            files.append(entry)
    html_resolved = html_out.resolve()
    files = [path for path in files if path.resolve() != html_resolved]
    return files, subfolders


def order_csv_candidates(stem_path: Path) -> List[Path]:
    name = stem_path.name
    candidates = [stem_path.with_name(name + "_columns.csv"), stem_path.with_name(name + "_order.csv")]
    if "motif_figure" in name:
        candidates.append(stem_path.with_name(name.replace("motif_figure", "motif_order") + ".csv"))
    return candidates


def build_figures(files: List[Path]) -> tuple:
    by_stem: Dict[Path, Dict[str, Path]] = defaultdict(dict)
    for path in files:
        if path.suffix.lower() in FIGURE_EXTS:
            by_stem[path.with_suffix("")][path.suffix.lower()] = path
    file_set = {path.resolve() for path in files}
    consumed: set = set()

    figures: List[Figure] = []
    for stem_path, formats in sorted(by_stem.items(), key=lambda item: item[0].as_posix().lower()):
        figure = Figure(stem=stem_path.name, files=formats)
        manifest_path = stem_path.with_name(stem_path.name + ".plot.json")
        if manifest_path.resolve() in file_set:
            figure.manifest = read_json(manifest_path)
            consumed.add(manifest_path.resolve())

        if figure.manifest:
            order_rel = (figure.manifest.get("outputs") or {}).get("order_csv")
            if order_rel:
                candidate = (manifest_path.parent / order_rel)
                if candidate.exists():
                    figure.order_csv = candidate
        if figure.order_csv is None:
            for candidate in order_csv_candidates(stem_path):
                if candidate.resolve() in file_set:
                    figure.order_csv = candidate
                    break
        if figure.order_csv is not None:
            consumed.add(figure.order_csv.resolve())

        classify_figure(figure)
        figures.append(figure)
    return figures, consumed


def classify_figure(figure: Figure) -> None:
    stem = figure.stem.lower()
    manifest = figure.manifest
    if manifest:
        figure.kind = "Motif dot matrix"
        category = str(manifest.get("category") or "All")
        figure.section = category if category in SECTION_ORDER else "Other"
        figure.title = str(manifest.get("title") or "")
        figure.subtitle = str(manifest.get("subtitle") or "")
        figure.caption = str(manifest.get("caption") or "")
        order_by = manifest.get("order_by") or "tree"
        figure.tags.append("Ordered by shared motifs" if order_by == "motifs" else "Ordered by Clustal Omega tree")
        options = manifest.get("options") or {}
        if options.get("dedupe_harmonized"):
            figure.tags.append("InterPro-deduped")
        else:
            figure.tags.append("Raw motif IDs")
        if options.get("hide_empty_columns"):
            figure.tags.append("Empty columns hidden")
        if options.get("min_proteins") is not None:
            figure.tags.append(f"Shared by >= {options['min_proteins']} proteins")
        if options.get("exclude_mobidblite") is False:
            figure.tags.append("MobiDBLite included")
        renderer = manifest.get("renderer")
        if renderer:
            figure.tags.append(f"{renderer} renderer")
        counts = manifest.get("counts") or {}
        figure.counts = {
            "proteins": counts.get("proteins_plotted"),
            "motifs": counts.get("motifs_plotted"),
            "dots": counts.get("dots_plotted"),
        }
        figure.groups = dict(counts.get("display_groups") or {})
        missing = manifest.get("missing_from_tree") or []
        if missing:
            figure.warnings.append(f"{len(missing)} protein(s) missing from the tree: {', '.join(map(str, missing))}")
        if figure.counts.get("motifs") == 0:
            figure.warnings.append("No shared motifs at this threshold; the matrix is empty.")
        return

    category = category_from_name(stem)
    if "schematic" in stem:
        figure.kind = "Linear schematic"
        figure.section = "Other"
    elif "motif" in stem and ("figure" in stem or "matrix" in stem) or figure.order_csv is not None:
        figure.kind = "Motif dot matrix"
        figure.section = category or "All"
    else:
        figure.section = "Other"
    if "dedup" in stem:
        figure.tags.append("InterPro-deduped")
    if "nonempty" in stem or "non_empty" in stem:
        figure.tags.append("Empty columns hidden")
    if figure.kind == "Motif dot matrix":
        figure.tags.append("details from file name")
    if figure.order_csv is not None:
        figure.groups = group_counts_from_order_csv(figure.order_csv)
        if figure.groups:
            figure.counts = {"proteins": sum(figure.groups.values())}


def preview_source(figure: Figure, html_dir: Path, embed: bool) -> Optional[str]:
    if embed:
        for ext in EMBED_EXTS:
            path = figure.files.get(ext)
            if path is None:
                continue
            mime = "image/svg+xml" if ext == ".svg" else "image/png"
            payload = base64.b64encode(path.read_bytes()).decode("ascii")
            return f"data:{mime};base64,{payload}"
        return None
    for ext in PREVIEW_EXTS:
        path = figure.files.get(ext)
        if path is not None:
            return rel_href(path, html_dir)
    return None


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def path_row(label: str, path: Path, html_dir: Path, root: Path) -> str:
    href = rel_href(path, html_dir)
    shown = rel_display(path, root)
    return (
        f'<li><span class="path-label">{esc(label)}</span>'
        f'<code class="path" data-href="{esc(href)}" title="{esc(shown)}">{esc(shown)}</code>'
        f'<button type="button" class="copy" data-href="{esc(href)}">Copy path</button></li>'
    )


def group_chips(groups: Dict[str, int]) -> str:
    if not groups:
        return ""
    chips = []
    for group, count in sorted(groups.items(), key=lambda item: group_sort_key(item[0])):
        color = GROUP_COLORS.get(group, "#777777")
        chips.append(
            f'<span class="chip"><span class="swatch" style="background:{esc(color)}"></span>'
            f"{esc(group)} <b>{esc(count)}</b></span>"
        )
    return f'<div class="chips">{"".join(chips)}</div>'


def group_sort_key(group: str) -> tuple:
    order = ["Pol I axial", "Pol I peri-axial", "Pol I", "Pol II", "Pol III"]
    return (order.index(group) if group in order else len(order), group)


def render_figure_card(figure: Figure, html_dir: Path, root: Path, embed: bool) -> str:
    preview = preview_source(figure, html_dir, embed)
    open_target = figure.files.get(".png") or figure.files.get(".svg") or figure.files.get(".pdf")
    open_href = rel_href(open_target, html_dir) if open_target is not None else "#"
    heading = figure.title or figure.stem
    if preview:
        preview_html = (
            f'<a class="preview" href="{esc(open_href)}" target="_blank" rel="noopener" title="Open full size">'
            f'<img src="{esc(preview)}" alt="{esc(heading)}" loading="lazy"></a>'
        )
    else:
        preview_html = f'<a class="preview empty" href="{esc(open_href)}" target="_blank" rel="noopener">PDF only - click to open</a>'

    links = []
    for ext, label in ((".svg", "SVG"), (".png", "PNG"), (".pdf", "PDF")):
        if ext in figure.files:
            links.append(f'<a class="btn" href="{esc(rel_href(figure.files[ext], html_dir))}" target="_blank" rel="noopener">{label}</a>')
    if figure.order_csv is not None:
        links.append(f'<a class="btn ghost" href="{esc(rel_href(figure.order_csv, html_dir))}" target="_blank" rel="noopener">Column order CSV</a>')

    counts = []
    for key in ("proteins", "motifs", "dots"):
        value = figure.counts.get(key)
        if value is not None:
            counts.append(f"{value} {key}")
    counts_html = f'<p class="counts">{esc(" · ".join(counts))}</p>' if counts else ""
    subtitle_html = f'<p class="sub">{esc(figure.subtitle)}</p>' if figure.subtitle else ""
    caption_html = f'<p class="caption">{esc(figure.caption)}</p>' if figure.caption else ""
    tags_html = "".join(f"<li>{esc(tag)}</li>" for tag in [figure.kind, *figure.tags])
    warnings_html = "".join(f'<p class="warn">{esc(message)}</p>' for message in figure.warnings)

    path_rows = [path_row(ext.lstrip(".").upper(), path, html_dir, root) for ext, path in sorted(figure.files.items())]
    if figure.order_csv is not None:
        path_rows.append(path_row("Order", figure.order_csv, html_dir, root))
    newest = max(figure.files.values(), key=lambda path: path.stat().st_mtime)
    search = " ".join([figure.stem, heading, figure.subtitle, figure.section, *figure.tags, *figure.groups]).lower()

    return f"""
<article class="card" data-search="{esc(search)}">
  {preview_html}
  <div class="card-body">
    <h3>{esc(heading)}</h3>
    {subtitle_html}
    <p class="file-name">{esc(figure.stem)}</p>
    <ul class="tags">{tags_html}</ul>
    {counts_html}
    {group_chips(figure.groups)}
    {warnings_html}
    {caption_html}
    <div class="links">{''.join(links)}</div>
    <details>
      <summary>File locations <span class="muted">(updated {esc(modified(newest))})</span></summary>
      <ul class="paths">{''.join(path_rows)}</ul>
    </details>
  </div>
</article>"""


def render_run_summary(summary: dict) -> str:
    rows = []
    if summary.get("input_file"):
        rows.append(("Input", f"<code>{esc(summary['input_file'])}</code>"))
    category_counts = summary.get("category_counts") or {}
    if category_counts:
        total = sum(int(value) for value in category_counts.values())
        detail = ", ".join(f"{esc(key)}: {esc(value)}" for key, value in sorted(category_counts.items()))
        rows.append(("Proteins", f"{total} ({detail})"))
    subgroup_counts = summary.get("pol_i_subgroup_counts") or {}
    if subgroup_counts:
        rows.append(("Pol I subgroups", ", ".join(f"{esc(key)}: {esc(value)}" for key, value in sorted(subgroup_counts.items()))))
    for label, key in (("Raw domains", "raw"), ("Harmonized (InterPro) domains", "harmonized")):
        block = summary.get(key) or {}
        if block:
            rows.append(
                (
                    label,
                    f"{esc(block.get('domain_hit_count', '?'))} hits, {esc(block.get('unique_domains', '?'))} unique, "
                    f"{esc(block.get('unique_domains_min_count', '?'))} shared by >= {esc(block.get('min_proteins_per_domain', '?'))} proteins",
                )
            )
    if not rows:
        return ""
    body = "".join(f"<tr><th>{esc(label)}</th><td>{value}</td></tr>" for label, value in rows)
    return f'<section id="run-summary"><h2>Run summary</h2><table class="kv">{body}</table></section>'


def render_tables(tables: List[Path], html_dir: Path, root: Path) -> str:
    if not tables:
        return ""
    rows = []
    for path in tables:
        shape = csv_shape(path) if path.suffix.lower() in {".csv", ".tsv"} else None
        shape_text = f"{shape[0]} rows × {shape[1]} cols" if shape else ""
        href = rel_href(path, html_dir)
        shown = rel_display(path, root)
        rows.append(
            f'<tr data-search="{esc(shown.lower())}">'
            f'<td><a href="{esc(href)}" target="_blank" rel="noopener">{esc(shown)}</a></td>'
            f"<td>{esc(describe_table(path.name))}</td>"
            f'<td class="num">{esc(shape_text)}</td>'
            f'<td class="num">{esc(human_size(path.stat().st_size))}</td>'
            f'<td class="num">{esc(modified(path))}</td>'
            f'<td><button type="button" class="copy" data-href="{esc(href)}">Copy path</button></td></tr>'
        )
    return (
        '<section id="tables"><h2>Tables and other files</h2>'
        '<div class="table-wrap"><table class="files"><thead><tr><th>File</th><th>What it is</th>'
        '<th class="num">Rows × cols</th><th class="num">Size</th><th class="num">Updated</th><th></th></tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></div></section>'
    )


def render_subfolders(subfolders: List[Path], html_dir: Path, recursive: bool) -> str:
    if not subfolders:
        return ""
    items = []
    for folder in subfolders:
        note = ""
        if folder.name.lower() == "cache":
            note = "API response cache (reused on the next run)"
        elif recursive and folder.name.lower() not in SKIP_DIRS:
            note = "Included above"
        nested_index = folder / "index.html"
        extra = ""
        if nested_index.exists():
            extra = f' · <a href="{esc(rel_href(nested_index, html_dir))}">its index</a>'
        items.append(
            f'<li><a href="{esc(rel_href(folder, html_dir))}/" target="_blank" rel="noopener">{esc(folder.name)}/</a>'
            f'{extra} <span class="muted">{esc(note)}</span></li>'
        )
    hint = "" if recursive else '<p class="muted">Not scanned. Rerun with <code>--recursive</code> to include their files.</p>'
    return f'<section id="folders"><h2>Subfolders</h2>{hint}<ul class="folders">{"".join(items)}</ul></section>'


def collect_warnings(figures: List[Figure], tables: List[Path]) -> List[str]:
    warnings = []
    for path in tables:
        if path.name.lower() in FAILURE_TABLES:
            shape = csv_shape(path)
            if shape and shape[0] > 0:
                warnings.append(f"{path.name} lists {shape[0]} failed lookup(s). Those proteins are missing domain data.")
    for figure in figures:
        for message in figure.warnings:
            warnings.append(f"{figure.stem}: {message}")
    return warnings


def render_html(
    outdir: Path,
    html_out: Path,
    title: str,
    recursive: bool,
    embed: bool,
) -> tuple:
    html_dir = html_out.parent
    files, subfolders = collect_files(outdir, recursive, html_out)
    figures, consumed = build_figures(files)
    consumed.update(path.resolve() for figure in figures for path in figure.files.values())
    tables = [
        path
        for path in files
        if path.resolve() not in consumed and (path.suffix.lower() in TABLE_EXTS or path.suffix.lower() not in FIGURE_EXTS)
    ]

    summary_path = outdir / "run_summary.json"
    summary = read_json(summary_path) if summary_path.exists() else None
    warnings = collect_warnings(figures, tables)

    by_section: Dict[str, List[Figure]] = defaultdict(list)
    for figure in figures:
        by_section[figure.section].append(figure)
    for section_figures in by_section.values():
        section_figures.sort(key=lambda fig: (fig.kind != "Motif dot matrix", fig.stem.lower()))

    nav = []
    sections_html = []
    for section in SECTION_ORDER:
        section_figures = by_section.get(section, [])
        if not section_figures:
            continue
        anchor = "fig-" + re.sub(r"[^a-z0-9]+", "-", section.lower()).strip("-")
        nav.append(f'<a href="#{anchor}">{esc(SECTION_TITLES[section])} <span class="count">{len(section_figures)}</span></a>')
        cards = "".join(render_figure_card(figure, html_dir, outdir, embed) for figure in section_figures)
        sections_html.append(
            f'<section id="{anchor}" class="fig-section"><h2>{esc(SECTION_TITLES[section])}</h2>'
            f'<div class="grid">{cards}</div></section>'
        )
    if summary:
        nav.insert(0, '<a href="#run-summary">Run summary</a>')
    if tables:
        nav.append(f'<a href="#tables">Tables <span class="count">{len(tables)}</span></a>')
    if subfolders:
        nav.append(f'<a href="#folders">Subfolders <span class="count">{len(subfolders)}</span></a>')

    warnings_html = ""
    if warnings:
        items = "".join(f"<li>{esc(message)}</li>" for message in warnings)
        warnings_html = f'<section id="warnings" class="warnings"><h2>Check these</h2><ul>{items}</ul></section>'

    folder_href = rel_href(outdir, html_dir)
    folder_href = folder_href + "/" if not folder_href.endswith("/") else folder_href
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    embed_note = (
        "Previews are embedded, so this page shows the figures even if it is sent on its own. "
        "The file links still need the output folder next to it."
        if embed
        else "Keep this file inside the output folder. Links are relative, so the folder can be moved, "
        "zipped, or shared between macOS and Windows."
    )

    page = PAGE_TEMPLATE.format(
        title=esc(title),
        generated=esc(generated),
        n_figures=len(figures),
        n_tables=len(tables),
        folder_href=esc(folder_href),
        folder_display=esc(outdir.resolve().as_posix()),
        note=esc(embed_note),
        nav="".join(nav),
        warnings=warnings_html,
        run_summary=render_run_summary(summary) if summary else "",
        sections="".join(sections_html) or '<p class="muted">No figures found in this folder.</p>',
        tables=render_tables(tables, html_dir, outdir),
        folders=render_subfolders(subfolders, html_dir, recursive),
        css=PAGE_CSS,
        js=PAGE_JS,
    )
    html_out.parent.mkdir(parents=True, exist_ok=True)
    html_out.write_text(page, encoding="utf-8")
    return len(figures), len(tables), warnings


PAGE_CSS = """
:root { --bg:#f7f7f5; --panel:#ffffff; --ink:#1d1d1f; --muted:#6b6b70; --line:#e2e2df; --accent:#0b5cad;
  --warn-bg:#fff6e0; --warn-line:#e7b54a; --chip:#f0f0ee; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#161618; --panel:#202023; --ink:#ececef; --muted:#a0a0a8; --line:#34343a; --accent:#6aa9ff;
    --warn-bg:#3a3020; --warn-line:#b38a2e; --chip:#2b2b30; }
}
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
header, main { max-width: 1400px; margin: 0 auto; padding: 0 20px; }
header { padding-top: 28px; }
h1 { font-size: 26px; margin: 0 0 4px; }
h2 { font-size: 19px; margin: 34px 0 12px; }
h3 { font-size: 15px; margin: 0; }
a { color: var(--accent); }
code { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12.5px; word-break: break-all; }
.muted { color: var(--muted); font-weight: normal; }
.meta { color: var(--muted); margin: 0 0 10px; }
.folder { display:flex; flex-wrap:wrap; gap:8px; align-items:center; background:var(--panel); border:1px solid var(--line);
  border-radius:8px; padding:8px 12px; }
.note { color: var(--muted); font-size: 13.5px; margin: 8px 0 0; }
nav { position: sticky; top: 0; z-index: 5; background: var(--bg); border-bottom: 1px solid var(--line);
  display:flex; flex-wrap:wrap; gap:6px 16px; align-items:center; padding: 10px 20px; margin-top: 16px; }
nav a { text-decoration: none; font-weight: 600; }
nav .count { color: var(--muted); font-weight: normal; }
nav input { margin-left: auto; min-width: 220px; padding: 6px 10px; border-radius: 6px; border: 1px solid var(--line);
  background: var(--panel); color: var(--ink); font: inherit; }
.grid { display:grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 16px; }
.card { background: var(--panel); border: 1px solid var(--line); border-radius: 10px; overflow: hidden; display:flex; flex-direction:column; }
.preview { display:flex; align-items:center; justify-content:center; height: 240px; background:#ffffff; border-bottom:1px solid var(--line); }
.preview img { max-width: 100%; max-height: 100%; object-fit: contain; }
.preview.empty { color:#555; text-decoration:none; font-size: 13px; }
.card-body { padding: 12px 14px 14px; display:flex; flex-direction:column; gap: 6px; }
.sub { margin:0; color: var(--muted); font-size: 13.5px; }
.file-name { margin:0; font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12px; color: var(--muted); word-break: break-all; }
.tags { list-style:none; margin: 2px 0 0; padding:0; display:flex; flex-wrap:wrap; gap:4px; }
.tags li { background: var(--chip); border-radius: 4px; padding: 1px 7px; font-size: 12px; }
.counts { margin: 0; font-size: 13.5px; }
.caption { margin: 0; font-size: 12.5px; color: var(--muted); line-height: 1.45; }
.chips { display:flex; flex-wrap:wrap; gap:4px 10px; font-size: 12.5px; }
.chip { display:inline-flex; align-items:center; gap:5px; }
.swatch { width: 10px; height: 10px; border-radius: 2px; display:inline-block; }
.links { display:flex; flex-wrap:wrap; gap:6px; margin-top: 4px; }
.btn { display:inline-block; padding: 4px 12px; border-radius: 6px; background: var(--accent); color: #fff; text-decoration:none;
  font-weight: 600; font-size: 13px; }
.btn.ghost { background: transparent; color: var(--accent); border: 1px solid var(--accent); }
details { font-size: 13px; }
summary { cursor: pointer; color: var(--muted); }
.paths { list-style:none; padding:0; margin: 6px 0 0; display:flex; flex-direction:column; gap:6px; }
.paths li { display:grid; grid-template-columns: 44px 1fr auto; gap: 6px; align-items:start; }
.path-label { font-weight: 600; font-size: 12px; color: var(--muted); padding-top: 1px; }
button.copy { font: inherit; font-size: 12px; padding: 2px 8px; border-radius: 5px; border:1px solid var(--line);
  background: var(--panel); color: var(--ink); cursor: pointer; white-space: nowrap; }
button.copy.done { border-color: #2e8b57; color: #2e8b57; }
.warn { margin: 0; font-size: 13px; background: var(--warn-bg); border-left: 3px solid var(--warn-line); padding: 4px 8px; }
.warnings { background: var(--warn-bg); border: 1px solid var(--warn-line); border-radius: 10px; padding: 2px 16px 8px; margin-top: 20px; }
.warnings h2 { margin-top: 12px; }
table { border-collapse: collapse; width: 100%; background: var(--panel); }
.table-wrap { overflow-x: auto; border: 1px solid var(--line); border-radius: 10px; }
th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--line); vertical-align: top; font-size: 13.5px; }
th { font-weight: 600; }
td.num, th.num { text-align: right; white-space: nowrap; color: var(--muted); }
table.kv { border: 1px solid var(--line); border-radius: 10px; max-width: 900px; }
table.kv th { width: 230px; color: var(--muted); font-weight: 600; }
.folders { padding-left: 18px; }
footer { max-width: 1400px; margin: 40px auto; padding: 0 20px 40px; color: var(--muted); font-size: 13px; }
.hidden { display: none !important; }
"""

PAGE_JS = r"""
// Turn a link relative to this page into a full path on the viewer's machine.
function localPath(href) {
  var url;
  try { url = new URL(href, window.location.href); } catch (e) { return href; }
  if (url.protocol !== 'file:') { return url.href; }
  var path = decodeURIComponent(url.pathname);
  if (/^\/[A-Za-z]:\//.test(path)) {               // Windows drive: /C:/Users/... -> C:\Users\...
    return path.slice(1).replace(/\//g, '\\').replace(/\\$/, '');
  }
  if (url.host) {                                   // Windows network share: \\server\share\...
    return ('\\\\' + url.host + path).replace(/\//g, '\\').replace(/\\$/, '');
  }
  return path.length > 1 ? path.replace(/\/$/, '') : path;  // macOS / Linux
}
function copyText(text, button) {
  function done() { button.textContent = 'Copied'; button.classList.add('done');
    setTimeout(function () { button.textContent = 'Copy path'; button.classList.remove('done'); }, 1500); }
  function fallback() {
    var area = document.createElement('textarea');
    area.value = text; document.body.appendChild(area); area.select();
    try { document.execCommand('copy'); done(); } catch (e) { window.prompt('Copy this path:', text); }
    document.body.removeChild(area);
  }
  if (navigator.clipboard && window.isSecureContext) {
    navigator.clipboard.writeText(text).then(done, fallback);
  } else { fallback(); }
}
document.querySelectorAll('code.path[data-href]').forEach(function (el) {
  el.textContent = localPath(el.getAttribute('data-href'));
});
document.querySelectorAll('button.copy[data-href]').forEach(function (button) {
  button.addEventListener('click', function () { copyText(localPath(button.getAttribute('data-href')), button); });
});
var filter = document.getElementById('filter');
if (filter) {
  filter.addEventListener('input', function () {
    var needle = filter.value.trim().toLowerCase();
    document.querySelectorAll('[data-search]').forEach(function (el) {
      el.classList.toggle('hidden', needle !== '' && el.getAttribute('data-search').indexOf(needle) === -1);
    });
    document.querySelectorAll('section.fig-section').forEach(function (section) {
      var visible = section.querySelectorAll('.card:not(.hidden)').length;
      section.classList.toggle('hidden', visible === 0);
    });
  });
}
"""

PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <p class="meta">Generated {generated} · {n_figures} figures · {n_tables} tables and other files</p>
  <div class="folder">
    <strong>Folder</strong>
    <code class="path" data-href="{folder_href}">{folder_display}</code>
    <button type="button" class="copy" data-href="{folder_href}">Copy path</button>
    <a href="{folder_href}" target="_blank" rel="noopener">Open folder</a>
  </div>
  <p class="note">{note} Click a preview to open it full size. Use "File locations" on a figure to see and copy each file's full path.</p>
</header>
<nav>{nav}<input id="filter" type="search" placeholder="Filter figures and files" aria-label="Filter figures and files"></nav>
<main>
{warnings}
{run_summary}
{sections}
{tables}
{folders}
</main>
<footer>Built by <code>scripts/build_output_index.py</code>. Rerun it after new plots are added to refresh this page.</footer>
<script>{js}</script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, required=True, help="Output folder to summarize.")
    parser.add_argument("--html-out", type=Path, help="Where to write the index (default: <outdir>/index.html).")
    parser.add_argument("--title", help="Page title (default: 'Output summary: <folder name>').")
    parser.add_argument("--recursive", action="store_true", help="Also include files in subfolders (cache/ is always skipped).")
    parser.add_argument(
        "--embed-previews",
        action="store_true",
        help="Inline the preview images so the page shows figures even when sent without the folder (larger file).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outdir = args.outdir
    if not outdir.is_dir():
        raise SystemExit(f"Output folder not found: {outdir}")
    html_out = args.html_out or (outdir / "index.html")
    title = args.title or f"Output summary: {outdir.resolve().name}"
    n_figures, n_tables, warnings = render_html(outdir, html_out, title, args.recursive, args.embed_previews)
    print(f"Wrote output index: {html_out}")
    print(f"Open it in a browser: {html_out.resolve().as_uri()}")
    print(f"Figures: {n_figures}  Tables/other files: {n_tables}  Warnings: {len(warnings)}")
    for message in warnings:
        print(f"  - {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
