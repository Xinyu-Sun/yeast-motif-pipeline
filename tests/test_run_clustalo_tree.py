import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import run_clustalo_tree


class RunClustaloTreeTests(unittest.TestCase):
    def test_ebi_result_identifier_selection_matches_job_dispatcher_outputs(self):
        result_types = {
            "aln-clustal_num": "Alignment",
            "phylotree": "Phylogenetic tree",
            "tree": "Guide Tree",
            "pim": "Percent identity matrix",
        }

        self.assertEqual(run_clustalo_tree.choose_result_identifier(result_types, "alignment"), "aln-clustal_num")
        self.assertEqual(run_clustalo_tree.choose_result_identifier(result_types, "phylo_tree"), "phylotree")
        self.assertEqual(run_clustalo_tree.choose_result_identifier(result_types, "guide_tree"), "tree")
        self.assertEqual(run_clustalo_tree.choose_result_identifier(result_types, "distmat"), "pim")

    def test_phylotree_never_falls_back_to_the_guide_tree(self):
        result_types = {"aln-clustal_num": "Alignment", "tree": "Guide Tree"}
        self.assertIsNone(run_clustalo_tree.choose_result_identifier(result_types, "phylo_tree"))

    def test_identity_ignores_gap_columns_and_marks_pairs_without_overlap(self):
        alignment = [("A", "MKV--"), ("B", "MRVLE"), ("C", "---LE")]
        identity, distance = run_clustalo_tree.identity_and_distance(alignment)
        self.assertAlmostEqual(identity[0][1], 200.0 / 3)  # 2 of 3 shared ungapped columns identical
        self.assertAlmostEqual(distance[0][1], 1.0 / 3)
        self.assertIsNone(identity[0][2])  # A and C share no ungapped column
        self.assertEqual(distance[0][2], 0.0)  # ClustalW/EBI convention
        text = run_clustalo_tree.format_percent_identity_matrix(["A", "B", "C"], identity)
        self.assertIn("     1: A   100.00   66.67    -nan", text)

    def test_neighbor_joining_recovers_an_additive_tree_in_clustalw_format(self):
        # Tree ((A:1,B:2):3, C:4, (D:5,E:6):7); distances are path lengths.
        lengths = {"A": 1, "B": 2, "C": 4, "D": 5, "E": 6}
        stem = {"A": 3, "B": 3, "C": 0, "D": 7, "E": 7}
        names = list("ABCDE")
        distance = [
            [0 if x == y else lengths[x] + lengths[y] + (0 if stem[x] == stem[y] and x in "ABDE" else stem[x] + stem[y]) for y in names]
            for x in names
        ]
        newick = run_clustalo_tree.neighbor_joining_newick(names, distance)
        self.assertEqual(newick, "(\n(\nA:1.00000,\nB:2.00000)\n:3.00000,\nC:4.00000,\n(\nD:5.00000,\nE:6.00000)\n:7.00000);\n")

    @unittest.skipUnless(shutil.which("clustalo"), "Clustal Omega is not installed")
    def test_local_backend_reproduces_the_lab_ebi_tree(self):
        csv_path = REPO_ROOT / "data" / "input" / "proteins.all48_plus_pol1_axial_periaxial.csv"
        ebi_tree = REPO_ROOT / "data" / "input" / "clustalo-all48-plus-pol1-axial-periaxial.ebi.phylotree"
        if not (csv_path.exists() and ebi_tree.exists()):
            self.skipTest("lab example input files are not present")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_clustalo_tree.py"),
                    "--proteins", str(csv_path),
                    "--fasta", str(tmpdir / "in.fasta"),
                    "--alignment-out", str(tmpdir / "aln.clu"),
                    "--out", str(tmpdir / "tree.phylotree"),
                    "--distmat-out", str(tmpdir / "pim.distmat"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((tmpdir / "tree.phylotree").read_bytes(), ebi_tree.read_bytes())

    def test_pol_i_subgroup_is_written_to_fasta_description(self):
        records = run_clustalo_tree.fasta_records(
            [
                {
                    "y_name": "YOR000W",
                    "standard_name": "Example1",
                    "amino_acid_sequence": "M K*Z",
                    "pol_i_subgroup": "axial",
                }
            ]
        )

        self.assertEqual(records[0], ">YOR000W Example1 axial")
        self.assertEqual(records[1], "MK*Z")


if __name__ == "__main__":
    unittest.main()
