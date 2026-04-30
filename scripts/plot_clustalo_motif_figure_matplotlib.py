#!/usr/bin/env python3
"""Matplotlib renderer for motif/domain dot-matrix figures.

This script keeps the native SVG writer intact and provides a standard
plotting-library implementation using the same input files and most of the same
CLI/style options as ``plot_clustalo_motif_figure.py``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

import plot_clustalo_motif_figure as base

PT_PER_SVG_PX = 72.0 / base.SVG_DPI


def load_matplotlib():
    try:
        import matplotlib.pyplot as plt
        from matplotlib import patches
    except ImportError as exc:
        raise SystemExit(
            "The Matplotlib renderer requires matplotlib. "
            "Install it with `python3 -m pip install matplotlib`."
        ) from exc
    return plt, patches


def mpl_font_family(style: dict) -> list[str]:
    families = [part.strip().strip("'\"") for part in str(style["font_family"]).split(",") if part.strip()]
    if "DejaVu Sans" not in families:
        families.append("DejaVu Sans")
    return families


def pt(value: float | str) -> float:
    return float(value) * PT_PER_SVG_PX


def draw_matplotlib_figure(
    tree: base.Node,
    tree_labels: List[str],
    domain_keys: List[str],
    ordered_rows: List[dict],
    title: str,
    subtitle: str,
    show_legend: bool,
    style: dict,
):
    plt, patches = load_matplotlib()
    layout = base.prepare_layout(tree, tree_labels, domain_keys, style)

    fig = plt.figure(
        figsize=(layout["width"] / base.SVG_DPI, layout["height"] / base.SVG_DPI),
        dpi=base.SVG_DPI,
        facecolor=str(style["background_color"]),
    )
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, layout["width"])
    ax.set_ylim(layout["height"], 0)
    ax.set_axis_off()
    ax.set_facecolor(str(style["background_color"]))

    category_colors = dict(base.DEFAULT_CATEGORY_COLORS)
    category_colors.update(style.get("category_colors", {}))
    font_family = mpl_font_family(style)

    ax.text(
        layout["title_x"],
        layout["title_y"],
        title,
        ha="center",
        va="baseline",
        fontsize=pt(style["title_font_size"]),
        fontweight="bold",
        family=font_family,
        color="black",
    )
    line_y = layout["title_y"] + float(style["title_line_offset"])
    line_half = float(style["title_line_half_length"])
    ax.plot(
        [layout["title_x"] - line_half, layout["title_x"] + line_half],
        [line_y, line_y],
        color=str(style["frame_color"]),
        linewidth=pt(style["title_line_width"]),
    )
    ax.text(
        layout["title_x"],
        layout["subtitle_y"],
        subtitle,
        ha="center",
        va="baseline",
        fontsize=pt(style["subtitle_font_size"]),
        family=font_family,
        color="black",
    )

    matrix_left = layout["matrix_left"]
    matrix_top = layout["matrix_top"]
    matrix_w = layout["matrix_w"]
    matrix_h = layout["matrix_h"]
    cell_w = layout["cell_w"]
    cell_h = layout["cell_h"]

    if show_legend:
        legend_x, legend_y = base.determine_legend_position(
            style=style,
            width=layout["width"],
            matrix_left=matrix_left,
            matrix_w=matrix_w,
            matrix_top=matrix_top,
            matrix_h=matrix_h,
            tree_top=layout["tree_top"],
            tree_h=layout["dendro_h"],
        )
        ax.text(
            legend_x,
            legend_y + 2,
            str(style["legend_title"]),
            ha="left",
            va="baseline",
            fontsize=pt(style["legend_title_font_size"]),
            fontweight="bold",
            family=font_family,
            color="black",
        )
        swatch_size = float(style["legend_swatch_size"])
        row_gap = float(style["legend_row_gap"])
        for idx, category in enumerate(["Pol I", "Pol II", "Pol III"]):
            y = legend_y + 12.0 + idx * row_gap
            ax.add_patch(
                patches.Rectangle(
                    (legend_x, y - swatch_size + 2),
                    swatch_size,
                    swatch_size,
                    facecolor=category_colors.get(category, "#000000"),
                    edgecolor=str(style["frame_color"]),
                    linewidth=pt(0.4),
                )
            )
            ax.text(
                legend_x + swatch_size + 5,
                y + 1,
                category,
                ha="left",
                va="baseline",
                fontsize=pt(style["legend_font_size"]),
                family=font_family,
                color="black",
            )

    for idx, row in enumerate(ordered_rows):
        short = row.get("standard_name") or row.get("protein") or row.get("y_name")
        x = matrix_left + (idx + 0.5) * cell_w
        ax.text(
            x,
            layout["label_y"],
            str(short),
            ha="left",
            va="center",
            rotation=90,
            rotation_mode="anchor",
            fontsize=pt(style["protein_label_font_size"]),
            family=font_family,
            color=category_colors.get(row.get("category", ""), "#000000"),
        )

    for idx in range(len(tree_labels) + 1):
        x = matrix_left + idx * cell_w
        ax.plot(
            [x, x],
            [matrix_top, matrix_top + matrix_h],
            color=str(style["grid_color"]),
            linewidth=pt(style["grid_width"]),
        )

    for row_idx, domain in enumerate(domain_keys):
        cy = matrix_top + (row_idx + 0.5) * cell_h
        label = base.short_domain_label(domain)
        ax.text(
            matrix_left - 8,
            cy + 3,
            label,
            ha="right",
            va="center",
            fontsize=pt(style["domain_label_font_size"]),
            family=font_family,
            color="black",
        )
        for col_idx, row in enumerate(ordered_rows):
            raw = str(row.get(domain, "0") or "0").strip()
            if raw in {"0", "", "0.0"}:
                continue
            cx = matrix_left + (col_idx + 0.5) * cell_w
            ax.add_patch(
                patches.Circle(
                    (cx, cy),
                    radius=float(style["dot_radius"]),
                    facecolor=str(style["dot_fill"]),
                    edgecolor=str(style["dot_stroke"]),
                    linewidth=pt(style["dot_stroke_width"]),
                )
            )

    for x1, y1, x2, y2 in base.iter_tree_segments(tree):
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=str(style["tree_stroke_color"]),
            linewidth=pt(style["tree_stroke_width"]),
        )

    ax.add_patch(
        patches.Rectangle(
            (matrix_left, matrix_top),
            matrix_w,
            matrix_h,
            fill=False,
            edgecolor=str(style["frame_color"]),
            linewidth=pt(style["frame_width"]),
        )
    )

    return fig


def save_figure(fig, path: Path, dpi: float | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches=None, pad_inches=0)


def main() -> int:
    args = base.parse_args()

    if args.write_default_style_config:
        args.write_default_style_config.parent.mkdir(parents=True, exist_ok=True)
        args.write_default_style_config.write_text(json.dumps(base.DEFAULT_STYLE, indent=2))
        print(f"Wrote default style config: {args.write_default_style_config}")
        return 0

    if args.min_proteins < 1:
        raise SystemExit("--min-proteins must be >= 1")

    style = base.build_style(args)
    base.validate_style(style)

    tree = base.NewickParser(args.tree.read_text()).parse()
    raw_tree_labels = [leaf.name or "" for leaf in base.iter_leaves(tree)]
    tree_labels = [base.normalize_tree_label(label) for label in raw_tree_labels]

    protein_meta = base.read_protein_metadata(args.proteins)
    domain_keys, matrix_by_y, _tree_rows, missing_from_tree = base.read_domain_hits(
        hits_path=args.hits,
        proteins_path=args.proteins,
        protein_meta=protein_meta,
        tree_labels=tree_labels,
        category=args.category,
        min_proteins=args.min_proteins,
        exclude_mobidblite=args.exclude_mobidblite,
    )

    ordered_rows: List[dict] = []
    resolved_tree_labels: List[str] = []
    for raw_label, norm_label in zip(raw_tree_labels, tree_labels):
        row = base.resolve_protein_row(raw_label, protein_meta) or base.resolve_protein_row(norm_label, protein_meta)
        if not row:
            continue
        if args.category and row["category"] != args.category:
            continue

        y_name = row["y_name"]
        if y_name not in matrix_by_y:
            continue
        ordered_rows.append(matrix_by_y[y_name])
        resolved_tree_labels.append(norm_label)

    if not ordered_rows:
        raise SystemExit("No proteins selected for plotting. Check --category and input files.")

    keep_labels = set(resolved_tree_labels)
    pruned_tree = base.prune_tree(tree, keep_labels)
    if pruned_tree is None:
        raise SystemExit("Could not prune tree to selected proteins.")

    base.write_order_csv(args.order_csv, ordered_rows, resolved_tree_labels)
    fig = draw_matplotlib_figure(
        tree=pruned_tree,
        tree_labels=resolved_tree_labels,
        domain_keys=domain_keys,
        ordered_rows=ordered_rows,
        title=args.title,
        subtitle=args.subtitle,
        show_legend=not args.no_legend,
        style=style,
    )

    save_figure(fig, args.out)
    if args.png:
        save_figure(fig, args.png, dpi=args.export_dpi)
    if args.pdf:
        save_figure(fig, args.pdf)

    print(f"Wrote Matplotlib figure: {args.out}")
    if args.png:
        print(f"Wrote Matplotlib PNG: {args.png}")
    if args.pdf:
        print(f"Wrote Matplotlib PDF: {args.pdf}")
    print(f"Wrote order table: {args.order_csv}")
    print(f"Proteins plotted: {len(ordered_rows)}")
    print(f"Shared domains plotted: {len(domain_keys)}")
    print(f"Excluded MobiDBLite: {args.exclude_mobidblite}")
    print(f"Min proteins threshold: {args.min_proteins}")
    if missing_from_tree:
        print("Missing from tree:", ", ".join(missing_from_tree))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
