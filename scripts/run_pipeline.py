#!/usr/bin/env python3
"""Run the standard pipeline for one protein CSV and finish with an HTML index.

One command, the same on macOS and Windows:

1. extract SGD/InterPro domains and SGD PTMs (skip with --skip-extract)
2. summarize shared domains by category
3. draw motif dot matrices for all proteins and for each category present
   (Pol I / Pol II / Pol III). Each comes with raw motif IDs and InterPro-deduped
   (empty columns hidden), and in two column orders: the Clustal Omega tree and
   clustering on shared motifs. Files are named
   <prefix>_<scope>_<raw|deduped>_by_<tree|motifs>.svg/.png/.pdf
4. draw the MobiDBLite linear schematic
5. write <outdir>/index.html

The tree (``--tree``) is only needed for tree-ordered figures. Make one with
``scripts/run_clustalo_tree.py``, or pass ``--order-by motifs`` to skip it.
Every step runs the same script documented in the README, and each command is
printed before it runs.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
CATEGORY_ORDER = ["Pol I", "Pol II", "Pol III"]


def run(step: str, args: list) -> None:
    cmd = [sys.executable, str(SCRIPTS / args[0]), *[str(arg) for arg in args[1:]]]
    print(f"\n== {step}\n$ " + " ".join(f'"{part}"' if " " in part else part for part in cmd), flush=True)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"Step failed ({step}); stopping. Fix the error above and rerun.")


def slug(text: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in text.lower()).strip("_")


def categories_present(proteins_csv: Path) -> list:
    with proteins_csv.open(newline="", encoding="utf-8-sig") as handle:
        present = {(row.get("category") or "").strip() for row in csv.DictReader(handle)}
    return [category for category in CATEGORY_ORDER if category in present]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, help="Protein CSV (see README 'Inputs'). Required unless --skip-extract.")
    parser.add_argument("--tree", type=Path, help="Newick tree for tree-ordered figures (not needed with --order-by motifs).")
    parser.add_argument(
        "--order-by",
        choices=["both", "tree", "motifs"],
        default="both",
        help="Column order(s) to draw: the Clustal Omega tree, clustering on shared motifs, or both (default).",
    )
    parser.add_argument("--outdir", type=Path, required=True, help="Folder for all outputs of this run.")
    parser.add_argument("--prefix", help="File-name prefix for figures (default: the output folder name).")
    parser.add_argument("--title", default="All proteins", help="Title of the all-protein figures.")
    parser.add_argument("--min-proteins", type=int, default=2, help="Plot motifs shared by at least this many proteins.")
    parser.add_argument("--export-dpi", type=float, default=600.0, help="PNG resolution.")
    parser.add_argument("--renderer", choices=["matplotlib", "native"], default="matplotlib")
    parser.add_argument("--skip-extract", action="store_true", help="Reuse the extraction tables already in --outdir.")
    parser.add_argument("--embed-previews", action="store_true", help="Pass --embed-previews to the HTML index.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outdir = args.outdir
    prefix = args.prefix or slug(outdir.resolve().name)
    orders = ["tree", "motifs"] if args.order_by == "both" else [args.order_by]
    if "tree" in orders and (args.tree is None or not args.tree.is_file()):
        raise SystemExit(f"Tree file not found: {args.tree}. Pass --tree, or use --order-by motifs.")

    if not args.skip_extract:
        if args.csv is None:
            raise SystemExit("--csv is required unless --skip-extract is given")
        run("Extract domains and PTMs", ["extract_sgd_domains.py", "--csv", args.csv, "--outdir", outdir,
                                         "--min-proteins-per-domain", args.min_proteins])
    proteins = outdir / "polymerase_tf_proteins.csv"
    hits = outdir / "polymerase_tf_domain_hits_raw.csv"
    for required in (proteins, hits):
        if not required.is_file():
            raise SystemExit(f"Missing {required}. Run without --skip-extract first.")

    run("Shared-domain summary", ["summarize_shared_domains.py", "--hits", hits,
                                  "--out", outdir / f"{prefix}_shared_domain_category_summary.csv",
                                  "--min-proteins", args.min_proteins])

    plot_script = "plot_clustalo_motif_figure_matplotlib.py" if args.renderer == "matplotlib" else "plot_clustalo_motif_figure.py"
    scopes = [(None, "all", args.title)] + [(category, slug(category), f"{category} proteins") for category in categories_present(proteins)]
    variants = [
        ("raw", [], f"shared motifs (>= {args.min_proteins} proteins)"),
        ("deduped", ["--dedupe-harmonized", "--hide-empty-columns"], "InterPro-deduped shared motifs"),
    ]
    for category, scope_slug, scope_title in scopes:
        for variant, flags, what in variants:
            for order in orders:
                stem = outdir / f"{prefix}_{scope_slug}_{variant}_by_{order}"
                title = f"{scope_title}: {what}"
                cmd = [plot_script, "--order-by", order, "--hits", hits, "--proteins", proteins,
                       "--out", stem.with_suffix(".svg"), "--png", stem.with_suffix(".png"), "--pdf", stem.with_suffix(".pdf"),
                       "--order-csv", stem.with_name(stem.name + "_columns.csv"),
                       "--title", title, "--min-proteins", args.min_proteins, "--export-dpi", args.export_dpi, *flags]
                if order == "tree":
                    cmd += ["--tree", args.tree]
                if category:
                    cmd += ["--category", category]
                run(f"Dot matrix: {title} (by {order})", cmd)

    run("MobiDBLite linear schematic", ["plot_mobidblite_linear_schematic.py", "--proteins", proteins, "--hits", hits,
                                        "--out", outdir / f"{prefix}_mobidblite_linear_schematic.svg",
                                        "--regions-csv", outdir / f"{prefix}_mobidblite_linear_regions.csv",
                                        "--summary-csv", outdir / f"{prefix}_mobidblite_linear_summary.csv",
                                        "--summary-md", outdir / f"{prefix}_mobidblite_linear_summary.md",
                                        "--title", f"{args.title}: domains and MobiDBLite disorder"])

    index_cmd = ["build_output_index.py", "--outdir", outdir, "--title", f"Output summary: {prefix}"]
    if args.embed_previews:
        index_cmd.append("--embed-previews")
    run("HTML output index", index_cmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
