import argparse
import csv
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import plot_clustalo_motif_figure as base

DOMAINS = ["x", "y", "z"]
PROFILES = {"P1": {"x", "y"}, "P2": {"z"}, "P3": {"x", "y"}, "P4": set(), "P5": {"x"}}


def rows_for(labels):
    return [{"y_name": label, **{d: int(d in PROFILES[label]) for d in DOMAINS}} for label in labels]


def leaf_names(node):
    return [leaf.name for leaf in base.iter_leaves(node)]


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class MotifClusterTreeTests(unittest.TestCase):
    def test_worked_example(self):
        # d(P1,P3) = 0 (same motifs); d(P1P3, P5) = 1 - 1/2 = 0.5; P2 and P4 share nothing with anyone.
        tree = base.motif_cluster_tree(rows_for(PROFILES), list(PROFILES), DOMAINS)
        self.assertEqual(leaf_names(tree), ["P1", "P3", "P5", "P2", "P4"])
        self.assertAlmostEqual(base.compute_heights(tree), 1.0)
        motif_group = tree.children[0]
        self.assertEqual(leaf_names(motif_group), ["P1", "P3", "P5"])
        self.assertAlmostEqual(motif_group.height, 0.5)
        self.assertAlmostEqual(motif_group.children[0].height, 0.0)  # P1 + P3 identical
        self.assertEqual(len(tree.children), 3)  # unrelated groups join in one node, empty protein last

    def test_order_does_not_depend_on_input_order(self):
        expected = leaf_names(base.motif_cluster_tree(rows_for(PROFILES), list(PROFILES), DOMAINS))
        rng = random.Random(0)
        for _ in range(20):
            labels = list(PROFILES)
            rng.shuffle(labels)
            self.assertEqual(leaf_names(base.motif_cluster_tree(rows_for(labels), labels, DOMAINS)), expected)

    def test_single_protein(self):
        self.assertEqual(leaf_names(base.motif_cluster_tree(rows_for(["P1"]), ["P1"], DOMAINS)), ["P1"])


class PreparePlotDataTests(unittest.TestCase):
    def test_motif_order_needs_no_tree_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            proteins = tmpdir / "proteins.csv"
            hits = tmpdir / "hits.csv"
            write_csv(
                proteins,
                [{"category": "Pol II", "protein": label, "y_name": label, "standard_name": label} for label in PROFILES],
                ["category", "protein", "y_name", "standard_name"],
            )
            write_csv(
                hits,
                [
                    {"category": "Pol II", "y_name": label, "raw_domain_key": f"Pfam::{d}", "source": "Pfam"}
                    for label, motifs in PROFILES.items()
                    for d in sorted(motifs)
                ],
                ["category", "y_name", "raw_domain_key", "source"],
            )
            args = argparse.Namespace(
                proteins=proteins, hits=hits, tree=tmpdir / "missing.nwk", order_by="motifs", category=None,
                min_proteins=1, exclude_mobidblite=True, dedupe_harmonized=False, hide_empty_columns=False,
            )
            data = base.prepare_plot_data(args)
        self.assertEqual(data["labels"], ["P1", "P3", "P5", "P2", "P4"])
        self.assertEqual([row["y_name"] for row in data["rows"]], data["labels"])
        self.assertEqual(data["missing_from_tree"], [])

    def test_defaults_describe_the_column_order(self):
        with mock.patch.object(sys, "argv", ["plot"]):
            args = base.parse_args()
        self.assertEqual((args.title, args.order_by), ("All proteins", "tree"))
        self.assertEqual(args.subtitle, "Ordered by Clustal Omega sequence similarity (neighbour-joining)")
        with mock.patch.object(sys, "argv", ["plot", "--order-by", "motifs", "--category", "Pol I"]):
            args = base.parse_args()
        self.assertEqual(args.title, "Pol I proteins")
        self.assertEqual(args.subtitle, "Ordered by shared motifs (Jaccard distance, average linkage)")

    def test_lab_example_keeps_identical_motif_sets_together(self):
        run_dir = REPO_ROOT / "data" / "output" / "all63"  # created from the lab example files, not in git
        if not (run_dir / "polymerase_tf_domain_hits_raw.csv").exists():
            self.skipTest("lab example run outputs are not present")
        for dedupe in (False, True):
            args = argparse.Namespace(
                proteins=run_dir / "polymerase_tf_proteins.csv", hits=run_dir / "polymerase_tf_domain_hits_raw.csv",
                tree=None, order_by="motifs", category=None, min_proteins=2, exclude_mobidblite=True,
                dedupe_harmonized=dedupe, hide_empty_columns=dedupe,
            )
            data = base.prepare_plot_data(args)
            profiles = [frozenset(d for d in data["domain_keys"] if base.has_hit(row, d)) for row in data["rows"]]
            for profile in set(p for p in profiles if p):
                positions = [idx for idx, p in enumerate(profiles) if p == profile]
                self.assertEqual(positions, list(range(positions[0], positions[-1] + 1)))


class RunnerTests(unittest.TestCase):
    def test_runner_requires_tree_only_for_tree_order(self):
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "run_pipeline.py"), "--outdir", "unused", "--skip-extract"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--order-by motifs", result.stderr)


if __name__ == "__main__":
    unittest.main()
