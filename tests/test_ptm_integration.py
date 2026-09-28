import csv
import tempfile
import unittest
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import extract_sgd_domains as extract
import run_mtprompt_ptm as mtprompt


class PtmExtractionTests(unittest.TestCase):
    def write_csv(self, path, rows, fieldnames):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def test_normalize_ptm_family(self):
        cases = {
            "phosphorylated residue": "phosphorylation",
            "ubiquitinated lysine": "ubiquitination",
            "SUMOylated residue": "sumoylation",
            "acetylated residue": "acetylation",
            "methylated residue": "methylation",
            "succinylated residue": "succinylation",
            "palmitoylated residue": "palmitoylation",
            "glycosylated residue": "other",
        }
        for text, expected in cases.items():
            self.assertEqual(extract.normalize_ptm_family(text), expected)

    def test_flatten_and_collapse_sgd_ptm_records(self):
        protein = {"category": "Pol II", "protein": "Example", "y_name": "YAL001C", "accession": "P00000"}
        protein_row = {
            **protein,
            "sgdid": "S000000001",
            "standard_name": "TFC3",
            "systematic_name": "YAL001C",
            "uniprot_id": "P00000",
        }
        records = [
            {
                "site_index": 10,
                "site_residue": "S",
                "type": "phosphorylated residue",
                "id": 1,
                "source": {"display_name": "SGD"},
                "reference": {"display_name": "A (2020)", "pubmed_id": 111, "link": "/reference/S1"},
                "locus": {"display_name": "TFC3", "format_name": "YAL001C"},
            },
            {
                "site_index": 10,
                "site_residue": "S",
                "type": "phosphorylated residue",
                "id": 2,
                "source": {"display_name": "SGD"},
                "reference": {"display_name": "B (2021)", "pubmed_id": 222, "link": "/reference/S2"},
                "locus": {"display_name": "TFC3", "format_name": "YAL001C"},
            },
        ]
        raw = extract.flatten_sgd_ptm_records(protein, protein_row, records)
        self.assertEqual(len(raw), 2)
        self.assertEqual(raw[0]["ptm_family"], "phosphorylation")
        collapsed = extract.collapse_ptm_evidence(raw)
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(collapsed[0]["evidence_record_count"], 2)
        self.assertEqual(collapsed[0]["ptm_evidence_ids"], "1; 2")
        self.assertEqual(collapsed[0]["pubmed_ids"], "111;222")

    def test_domain_relations(self):
        domains = [
            {"raw_domain_key": "Pfam::A", "start": "20", "end": "40"},
            {"raw_domain_key": "Pfam::B", "start": "80", "end": "100"},
        ]
        relation, matches = extract.classify_ptm_domain_relation(25, domains, 10)
        self.assertEqual(relation, "inside_domain")
        self.assertEqual(matches[0]["raw_domain_key"], "Pfam::A")

        relation, matches = extract.classify_ptm_domain_relation(70, domains, 10)
        self.assertEqual(relation, "near_domain_boundary")
        self.assertEqual(matches[0]["raw_domain_key"], "Pfam::B")

        relation, matches = extract.classify_ptm_domain_relation(55, domains, 10)
        self.assertEqual(relation, "outside_domains")
        self.assertEqual(matches, [])

        relation, matches = extract.classify_ptm_domain_relation(55, [], 10)
        self.assertEqual(relation, "no_domains")
        self.assertEqual(matches, [])

    def test_csv_metadata_is_preserved_in_domain_matrices(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "proteins.csv"
            self.write_csv(
                path,
                [
                    {
                        "category": "Pol I",
                        "protein": "Example",
                        "y_name": "YAL001C",
                        "pol_i_subgroup": "axial",
                    }
                ],
                ["category", "protein", "y_name", "pol_i_subgroup"],
            )

            entries = extract.read_csv_entries(path)
            self.assertEqual(entries[0]["pol_i_subgroup"], "axial")

            raw_all, _raw_filtered, _raw_counts, _raw_filtered_counts, domains, _filtered = extract.build_domain_matrices(
                [{"y_name": "YAL001C", "raw_domain_key": "Pfam::A"}],
                entries,
                key_field="raw_domain_key",
                min_count=1,
            )

            self.assertEqual(domains, ["Pfam::A"])
            self.assertEqual(raw_all[0]["pol_i_subgroup"], "axial")
            self.assertEqual(raw_all[0]["Pfam::A"], 1)


class MtpromptTests(unittest.TestCase):
    def write_csv(self, path, rows, fieldnames):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def test_prediction_parsing_rejects_residue_mismatch_and_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            proteins = [
                {
                    "category": "Pol II",
                    "protein": "Example",
                    "y_name": "YAL001C",
                    "standard_name": "TFC3",
                    "uniprot_id": "P00000",
                    "amino_acid_sequence": "MSTK",
                }
            ]
            prediction_dir = tmpdir / "predictions"
            prediction_dir.mkdir()
            self.write_csv(
                prediction_dir / "Phosphorylation_S.csv",
                [
                    {"protein_id": "YAL001C", "position": "2", "score": "0.9"},
                    {"protein_id": "YAL001C", "position": "3", "residue": "S", "score": "0.8"},
                ],
                ["protein_id", "position", "residue", "score"],
            )

            predictions = mtprompt.normalize_prediction_rows(prediction_dir, proteins)
            self.assertEqual(len(predictions), 1)
            self.assertEqual(predictions[0]["site_index"], 2)

            known = [
                {
                    "y_name": "YAL001C",
                    "site_index": "2",
                    "site_residue": "S",
                    "ptm_family": "phosphorylation",
                    "ptm_evidence_ids": "1; 2",
                }
            ]
            candidates = mtprompt.build_candidate_sites(predictions, known)
            self.assertEqual(candidates[0]["evidence_tier"], "predicted_known_overlap")
            self.assertEqual(candidates[0]["known_sgd_ptm_evidence_ids"], "1; 2")

    def test_candidate_table_keeps_known_sgd_without_prediction(self):
        candidates = mtprompt.build_candidate_sites(
            predictions=[],
            known_ptms=[
                {
                    "category": "Pol II",
                    "protein": "Example",
                    "y_name": "YAL001C",
                    "standard_name": "TFC3",
                    "site_index": "2",
                    "site_residue": "S",
                    "ptm_family": "phosphorylation",
                    "ptm_evidence_ids": "1",
                }
            ],
        )
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["evidence_tier"], "known_sgd")
        self.assertEqual(candidates[0]["prediction_count"], 0)


if __name__ == "__main__":
    unittest.main()
