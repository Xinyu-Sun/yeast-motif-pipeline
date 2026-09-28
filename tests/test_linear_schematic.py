import csv
import tempfile
import unittest
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import plot_mobidblite_linear_schematic as schematic


class LinearSchematicTests(unittest.TestCase):
    def write_csv(self, path, rows, fieldnames):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def test_build_regions_splits_domains_and_mobidblite(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            hits = tmpdir / "hits.csv"
            proteins = [
                {
                    "category": "Pol II",
                    "protein": "Example",
                    "standard_name": "EXA",
                    "y_name": "YAL001C",
                    "accession": "P00000",
                    "protein_length": 100,
                }
            ]
            self.write_csv(
                hits,
                [
                    {
                        "category": "Pol II",
                        "protein": "Example",
                        "y_name": "YAL001C",
                        "source": "Pfam",
                        "raw_domain_key": "Pfam::PF00001",
                        "domain_description": "Example domain",
                        "start": "10",
                        "end": "30",
                    },
                    {
                        "category": "Pol II",
                        "protein": "Example",
                        "y_name": "YAL001C",
                        "source": "MobiDBLite",
                        "raw_domain_key": "MobiDBLite::MobiDBLite",
                        "domain_description": "consensus disorder prediction",
                        "start": "80",
                        "end": "120",
                    },
                ],
                ["category", "protein", "y_name", "source", "raw_domain_key", "domain_description", "start", "end"],
            )

            regions = schematic.build_regions(hits, proteins)

            self.assertEqual([row["feature_type"] for row in regions], ["domain", "disorder"])
            self.assertEqual(regions[0]["label"], "PF00001")
            self.assertEqual(regions[1]["end"], 100)
            self.assertEqual(regions[1]["region_length"], 21)

    def test_merged_coverage_and_summary(self):
        proteins = [{"category": "Pol I", "protein": "Example", "standard_name": "EXA", "y_name": "YAL001C", "protein_length": 100}]
        regions = [
            {"y_name": "YAL001C", "feature_type": "domain", "start": 1, "end": 20},
            {"y_name": "YAL001C", "feature_type": "domain", "start": 15, "end": 30},
            {"y_name": "YAL001C", "feature_type": "disorder", "start": 60, "end": 70},
        ]

        summaries = schematic.summarize_regions(proteins, regions)

        self.assertEqual(summaries[0]["domain_coverage_aa"], 30)
        self.assertEqual(summaries[0]["disorder_coverage_aa"], 11)
        self.assertEqual(summaries[0]["domain_region_count"], 2)

    def test_render_svg_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "schematic.svg"
            proteins = [{"category": "Pol I", "protein": "Example", "standard_name": "EXA", "y_name": "YAL001C", "protein_length": 100}]
            regions = [
                {
                    "y_name": "YAL001C",
                    "feature_type": "disorder",
                    "raw_domain_key": "MobiDBLite::MobiDBLite",
                    "label": "MobiDBLite disorder",
                    "start": 60,
                    "end": 70,
                }
            ]

            schematic.render_svg(proteins, regions, [], out, "Test")

            text = out.read_text(encoding="utf-8")
            self.assertIn("<svg", text)
            self.assertIn("MobiDBLite disorder", text)


if __name__ == "__main__":
    unittest.main()
