import inspect
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import plot_clustalo_motif_figure_matplotlib as mpl_plot
import plot_clustalo_motif_figure as svg_plot


class MatrixRendererNoPtmTests(unittest.TestCase):
    def test_ptm_overlay_flags_are_not_documented_for_matrix_renderer(self):
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "plot_clustalo_motif_figure_matplotlib.py"), "--help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        self.assertNotIn("--show-ptm-overlay", result.stdout)
        self.assertNotIn("--ptm-intersections", result.stdout)

    def test_matplotlib_draw_function_has_no_ptm_overlay_input(self):
        signature = inspect.signature(mpl_plot.draw_matplotlib_figure)
        self.assertNotIn("ptm_overlay", signature.parameters)

    def test_pol_i_subgroup_display_groups_are_ordered_for_legend(self):
        rows = [
            {"category": "Pol III"},
            {"category": "Pol I", "pol_i_subgroup": "axial"},
            {"category": "Pol I", "pol_i_subgroup": "transcription_factor"},
            {"category": "Pol II"},
            {"category": "Pol I", "pol_i_subgroup": "periaxial"},
        ]

        self.assertEqual(svg_plot.display_group(rows[1]), "Pol I axial")
        self.assertEqual(
            svg_plot.ordered_display_groups(rows),
            ["Pol I axial", "Pol I peri-axial", "Pol II", "Pol III"],
        )

    def test_default_dot_colors_follow_display_group(self):
        style = dict(svg_plot.DEFAULT_STYLE)
        category_colors = dict(svg_plot.DEFAULT_CATEGORY_COLORS)
        fill, stroke = svg_plot.dot_colors(
            {"category": "Pol I", "pol_i_subgroup": "periaxial"},
            style,
            category_colors,
        )

        self.assertEqual(fill, category_colors["Pol I peri-axial"])
        self.assertEqual(stroke, category_colors["Pol I peri-axial"])

    def test_deduped_matrix_collapses_mapped_interpro_hits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            proteins = tmp / "proteins.csv"
            hits = tmp / "hits.csv"
            proteins.write_text(
                "category,protein,y_name,standard_name,pol_i_subgroup,fasta_header\n"
                "Pol II,A,YAL001C,A1,,>sp|P1|YAL001C protein A\n"
                "Pol III,B,YBL001C,B1,,>sp|P2|YBL001C protein B\n"
            )
            hits.write_text(
                "category,y_name,raw_domain_key,source,harmonization_status,harmonized_key,harmonized_accession,harmonized_name,domain_description\n"
                "Pol II,YAL001C,SGD::alpha,SGD,mapped,InterPro::IPR000001,IPR000001,Shared domain,Alpha source\n"
                "Pol II,YAL001C,Pfam::alpha,InterPro,mapped,InterPro::IPR000001,IPR000001,Shared domain,Alpha pfam\n"
                "Pol III,YBL001C,SGD::alpha2,SGD,mapped,InterPro::IPR000001,IPR000001,Shared domain,Alpha source 2\n"
            )
            meta = svg_plot.read_protein_metadata(proteins)

            domain_keys, matrix_by_y, _tree_rows, _missing, key_meta = svg_plot.read_deduped_domain_hits(
                hits_path=hits,
                proteins_path=proteins,
                protein_meta=meta,
                tree_labels=["YAL001C", "YBL001C"],
                category=None,
                min_proteins=2,
                exclude_mobidblite=True,
            )

            self.assertEqual(domain_keys, ["InterPro::IPR000001"])
            self.assertEqual(matrix_by_y["YAL001C"]["InterPro::IPR000001"], 1)
            self.assertEqual(matrix_by_y["YBL001C"]["InterPro::IPR000001"], 1)
            self.assertEqual(key_meta["InterPro::IPR000001"]["status"], "mapped")

    def test_svg_protein_labels_use_column_center_baseline(self):
        text = inspect.getsource(svg_plot.build_svg)

        self.assertIn('transform="translate({x:.2f},{label_y:.2f}) rotate(-90)"', text)
        self.assertIn('dominant-baseline="middle"', text)

    def test_auto_margin_expands_for_long_deduped_labels(self):
        tree = svg_plot.NewickParser("(YAL001C:0.1,YBL001C:0.2);").parse()
        style = dict(svg_plot.DEFAULT_STYLE)
        style["margin_left"] = 120
        domain_meta = {
            "InterPro::IPR999999": {
                "status": "mapped",
                "label": "Very long InterPro motif label that should remain fully visible",
                "accession": "IPR999999",
            }
        }

        layout = svg_plot.prepare_layout(
            tree,
            ["YAL001C", "YBL001C"],
            ["InterPro::IPR999999"],
            style,
            domain_meta=domain_meta,
        )

        self.assertGreater(layout["margin_left"], 300)

    def test_sparse_category_layout_keeps_readable_canvas_and_matrix(self):
        tree = svg_plot.NewickParser("(TFC6:0.1,TFC8:0.2);").parse()
        style = svg_plot.style_with_text_context(
            dict(svg_plot.DEFAULT_STYLE),
            "Pol III proteins",
            "EBI Clustal Omega tree; empty protein columns hidden",
        )

        layout = svg_plot.prepare_layout(tree, ["TFC6", "TFC8"], ["A", "B"], style)

        self.assertGreaterEqual(layout["width"], svg_plot.DEFAULT_STYLE["min_figure_width"])
        self.assertGreaterEqual(layout["matrix_w"], svg_plot.DEFAULT_STYLE["min_matrix_width"])
        self.assertGreater(layout["cell_w"], svg_plot.DEFAULT_STYLE["cell_width"])

    def test_no_shared_motifs_with_hide_empty_columns_keeps_plot_structure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            proteins = tmp / "proteins.csv"
            hits = tmp / "hits.csv"
            tree = tmp / "tree.phylotree"
            out = tmp / "out.svg"
            order = tmp / "order.csv"
            proteins.write_text(
                "category,protein,y_name,standard_name,pol_i_subgroup,fasta_header\n"
                "Pol II,A,YAL001C,A1,,>sp|P1|YAL001C protein A\n"
                "Pol II,B,YBL001C,B1,,>sp|P2|YBL001C protein B\n"
            )
            hits.write_text(
                "category,y_name,raw_domain_key,source,harmonization_status,harmonized_key,harmonized_accession,harmonized_name,domain_description\n"
                "Pol II,YAL001C,SGD::singleton,SGD,,,,,Singleton domain\n"
            )
            tree.write_text("(YAL001C:0.1,YBL001C:0.2);")

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "plot_clustalo_motif_figure.py"),
                    "--tree",
                    str(tree),
                    "--hits",
                    str(hits),
                    "--proteins",
                    str(proteins),
                    "--out",
                    str(out),
                    "--order-csv",
                    str(order),
                    "--min-proteins",
                    "3",
                    "--hide-empty-columns",
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            self.assertIn("Proteins plotted: 2", result.stdout)
            self.assertIn("Shared domains plotted: 0", result.stdout)
            svg = out.read_text()
            self.assertIn('height="14.00"', svg)
            self.assertIn('y="54.00"', svg)


if __name__ == "__main__":
    unittest.main()
