import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import build_pol1_extended_input as pol1


def uniprot_hit(gene, locus, synonyms=()):
    return {
        "primaryAccession": "P00001",
        "uniProtkbId": f"{gene}_YEAST",
        "genes": [
            {
                "geneName": {"value": gene},
                "orderedLocusNames": [{"value": locus}],
                "synonyms": [{"value": name} for name in synonyms],
            }
        ],
        "sequence": {"value": "MKV"},
    }


class UniProtLookupTests(unittest.TestCase):
    def build(self, cache_dir, query, hit, workbook_y_name="YAL001C"):
        hits = hit if isinstance(hit, list) else [hit]
        (cache_dir / pol1.cache_name(query)).write_text(json.dumps({"results": hits}), encoding="utf-8")
        source = {"protein": query, "workbook_y_name": workbook_y_name, "source_row": "3", "length_aa": ""}
        return pol1.build_pol_i_row(source, cache_dir=cache_dir, offline=True, timeout=1.0)

    def test_matching_gene_is_used_and_locus_correction_is_noted(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = self.build(Path(tmp), "NOP10", uniprot_hit("NOP10", "YHR072W-A"), workbook_y_name="YDL208W")
        self.assertEqual(row["y_name"], "YHR072W-A")
        self.assertIn("corrected to YHR072W-A", row["source_note"])

    def test_unique_alias_match_is_accepted_and_noted(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = self.build(Path(tmp), "SOF3", uniprot_hit("NOP10", "YHR072W-A", synonyms=["SOF3"]))
        self.assertEqual(row["protein"], "NOP10")
        self.assertIn("matched UniProt alias of NOP10", row["source_note"])

    def test_main_name_match_wins_over_an_earlier_alias_match(self):
        hits = [uniprot_hit("OTHER1", "YBR001W", synonyms=["ABC1"]), uniprot_hit("ABC1", "YAL002W")]
        with tempfile.TemporaryDirectory() as tmp:
            row = self.build(Path(tmp), "ABC1", hits, workbook_y_name="YAL002W")
        self.assertEqual((row["protein"], row["y_name"], row["source_note"]), ("ABC1", "YAL002W", ""))

    def test_ambiguous_alias_is_refused(self):
        hits = [uniprot_hit("GENE1", "YBR001W", synonyms=["ABC1"]), uniprot_hit("GENE2", "YBR002W", synonyms=["ABC1"])]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as caught:
                self.build(Path(tmp), "ABC1", hits)
        self.assertIn("several genes list it as an alias", str(caught.exception))

    def test_hit_for_a_different_gene_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as caught:
                self.build(Path(tmp), "ABC1", uniprot_hit("XYZ2", "YBR999W"))
        self.assertIn("no hit lists it as a gene name", str(caught.exception))

    def test_cached_lookups_for_the_lab_example_pass_the_check(self):
        cache = REPO_ROOT / "data" / "output" / "all63" / "cache" / "uniprot"  # lab example files, not in git
        cached = sorted(cache.glob("*.json")) if cache.is_dir() else []
        if not cached:
            self.skipTest("UniProt cache for the lab example is not present")
        for path in cached:
            source = {"protein": path.stem, "workbook_y_name": "", "source_row": "?", "length_aa": ""}
            row = pol1.build_pol_i_row(source, cache_dir=cache, offline=True, timeout=1.0)
            self.assertTrue(row["amino_acid_sequence"], path.stem)


if __name__ == "__main__":
    unittest.main()
