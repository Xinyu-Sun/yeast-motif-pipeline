import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import build_output_index as index
import plot_clustalo_motif_figure as base


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class OutputIndexTests(unittest.TestCase):
    def test_category_from_file_name(self):
        cases = {
            "example_motif_figure_ebi_pol_i_deduped": "Pol I",
            "example_motif_figure_ebi_pol_ii": "Pol II",
            "example_motif_figure_ebi_pol_iii_nonempty_columns": "Pol III",
            "pol2_motif_figure": "Pol II",
            "polymerase_tf_domain_matrix_all": None,
            "all48_motif_figure": None,
        }
        for stem, expected in cases.items():
            self.assertEqual(index.category_from_name(stem), expected, stem)

    def test_index_groups_figures_links_relatively_and_flags_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "run with spaces"
            outdir.mkdir()
            for ext in (".svg", ".png", ".pdf"):
                (outdir / f"example_motif_figure_ebi_pol_ii{ext}").write_bytes(b"x")
            (outdir / "Pol I figure.svg").write_text("<svg/>", encoding="utf-8")
            (outdir / "Pol I figure.plot.json").write_text(
                json.dumps(
                    {
                        "category": "Pol I",
                        "title": "Pol I TFs",
                        "options": {"dedupe_harmonized": True, "hide_empty_columns": True, "min_proteins": 2},
                        "outputs": {"order_csv": "Pol I figure_order.csv"},
                        "counts": {"proteins_plotted": 2, "motifs_plotted": 3, "display_groups": {"Pol I axial": 2}},
                        "missing_from_tree": ["ABC1"],
                    }
                ),
                encoding="utf-8",
            )
            write_csv(outdir / "Pol I figure_order.csv", [{"y_name": "Y1", "display_group": "Pol I axial"}], ["y_name", "display_group"])
            write_csv(outdir / "example_motif_order_ebi_pol_ii.csv", [{"y_name": "Y2", "category": "Pol II"}], ["y_name", "category"])
            write_csv(outdir / "sgd_failures.csv", [{"y_name": "Y9", "error": "timeout"}], ["y_name", "error"])
            write_csv(outdir / "polymerase_tf_proteins.csv", [{"y_name": "Y1"}], ["y_name"])
            (outdir / "cache").mkdir()

            html_out = outdir / "index.html"
            n_figures, n_tables, warnings = index.render_html(outdir, html_out, "Test run", False, False)
            page = html_out.read_text(encoding="utf-8")

            self.assertEqual(n_figures, 2)
            self.assertIn('id="fig-pol-i"', page)
            self.assertIn('id="fig-pol-ii"', page)
            self.assertIn('href="Pol%20I%20figure.svg"', page)
            self.assertIn('href="Pol%20I%20figure_order.csv"', page)
            self.assertIn('href="example_motif_order_ebi_pol_ii.csv"', page)
            self.assertIn("InterPro-deduped", page)
            self.assertIn("2 proteins · 3 motifs", page)
            hrefs = re.findall(r'(?:href|src|data-href)="([^"#]+)"', page)
            self.assertTrue(hrefs)
            for href in hrefs:
                self.assertFalse(href.startswith(("/", "file:")), href)
                self.assertNotIn("\\", href)
            self.assertEqual(n_tables, 2)  # order CSVs and the manifest are attached to their figures
            self.assertTrue(any("sgd_failures.csv lists 1" in message for message in warnings))
            self.assertTrue(any("ABC1" in message for message in warnings))

            # Rebuilding must not list the previous index.html as a table.
            index.render_html(outdir, html_out, "Test run", False, False)
            self.assertNotIn('href="index.html"', html_out.read_text(encoding="utf-8"))

    def test_embed_previews_inlines_images(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            (outdir / "all_motif_figure.svg").write_text("<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8")
            html_out = outdir / "index.html"
            index.render_html(outdir, html_out, "Embedded", False, True)
            self.assertIn('src="data:image/svg+xml;base64,', html_out.read_text(encoding="utf-8"))


class PlotManifestTests(unittest.TestCase):
    def test_manifest_records_options_counts_and_relative_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            args = argparse.Namespace(
                out=outdir / "pol_i.svg",
                png=outdir / "pol_i.png",
                pdf=None,
                order_csv=outdir / "pol_i_order.csv",
                tree=Path("data/input/tree.phylotree"),
                hits=Path("data/output/hits.csv"),
                proteins=Path("data/output/proteins.csv"),
                title="Pol I TFs",
                subtitle="Deduped",
                category="Pol I",
                order_by="tree",
                min_proteins=2,
                exclude_mobidblite=True,
                dedupe_harmonized=True,
                hide_empty_columns=False,
                export_dpi=600.0,
            )
            rows = [
                {"y_name": "Y1", "category": "Pol I", "pol_i_subgroup": "axial", "D1": 1, "D2": 0},
                {"y_name": "Y2", "category": "Pol I", "pol_i_subgroup": "periaxial", "D1": 1, "D2": 1},
            ]
            path = base.write_plot_manifest(args, "matplotlib", rows, ["D1", "D2"], [])
            manifest = json.loads(path.read_text(encoding="utf-8"))

            self.assertEqual(path.name, "pol_i.plot.json")
            self.assertEqual(manifest["category"], "Pol I")
            self.assertEqual(manifest["outputs"]["png"], "pol_i.png")
            self.assertIsNone(manifest["outputs"]["pdf"])
            self.assertEqual(manifest["inputs"]["tree"], "data/input/tree.phylotree")
            self.assertEqual(manifest["counts"]["dots_plotted"], 3)
            self.assertEqual(manifest["counts"]["display_groups"], {"Pol I axial": 1, "Pol I peri-axial": 1})
            self.assertEqual(manifest["order_by"], "tree")
            self.assertIn("display order", manifest["caption"])


class LinearSchematicCliTests(unittest.TestCase):
    def test_cli_writes_region_table_with_harmonized_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            write_csv(
                tmpdir / "polymerase_tf_proteins.csv",
                [{"category": "Pol II", "protein": "Example", "y_name": "YAL001C", "amino_acid_sequence": "M" * 100}],
                ["category", "protein", "y_name", "amino_acid_sequence"],
            )
            write_csv(
                tmpdir / "hits.csv",
                [
                    {"y_name": "YAL001C", "source": "Pfam", "raw_domain_key": "Pfam::PF1", "harmonized_key": "InterPro::IPR1", "start": "10", "end": "30"},
                    {"y_name": "YAL001C", "source": "Pfam", "raw_domain_key": "Pfam::PF2", "harmonized_key": "", "start": "150", "end": "200"},
                ],
                ["y_name", "source", "raw_domain_key", "harmonized_key", "start", "end"],
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "plot_mobidblite_linear_schematic.py"),
                    "--proteins", str(tmpdir / "polymerase_tf_proteins.csv"),
                    "--hits", str(tmpdir / "hits.csv"),
                    "--out", str(tmpdir / "s.svg"),
                    "--regions-csv", str(tmpdir / "regions.csv"),
                    "--summary-csv", str(tmpdir / "summary.csv"),
                    "--summary-md", str(tmpdir / "summary.md"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (tmpdir / "regions.csv").open(newline="", encoding="utf-8") as handle:
                regions = list(csv.DictReader(handle))
            # The out-of-range hit (150-200 on a 100-aa protein) is dropped, not clamped to a 1-aa region.
            self.assertEqual([row["harmonized_key"] for row in regions], ["InterPro::IPR1"])
            self.assertTrue((tmpdir / "summary.md").exists())


if __name__ == "__main__":
    unittest.main()
