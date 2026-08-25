"""Tests for the independent published S3 graph-data replay."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import validate_published_s3 as validation  # noqa: E402

_ARTIFACT = Path(__file__).resolve().parents[1] / "data" / "published_s3_validation.json"


class TestGapParser(unittest.TestCase):
    def test_two_gap_records_are_parsed_without_executing_gap(self):
        payload = b"""gama:=[ rec(
adjacencies := [ [ 2 ], [ 1, 3 ], [ 2 ] ],
autGroup := Group( [ () ] ) ), rec(
adjacencies := [ [ 2, 3 ], [ 1, 3 ], [ 1, 2 ] ],
autGroup := Group( [ () ] ) ) ];
"""
        graphs = validation.parse_gap_adjacencies(payload)
        self.assertEqual(graphs, [
            (1 << 1, (1 << 0) | (1 << 2), 1 << 1),
            ((1 << 1) | (1 << 2), (1 << 0) | (1 << 2),
             (1 << 0) | (1 << 1)),
        ])

    def test_unbalanced_payload_is_rejected(self):
        with self.assertRaises(validation.ValidationViolation):
            validation.parse_gap_adjacencies(b"adjacencies := [ [ 2 ]")

    def test_published_involution_has_fixed_five_cycle_type(self):
        involution = validation._specified_involution()
        self.assertEqual(sum(involution[v] == v for v in range(45)), 5)
        self.assertEqual(sum(v < involution[v] for v in range(45)), 20)
        self.assertTrue(all(involution[involution[v]] == v for v in range(45)))


class TestPublishedArtifact(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = validation.load_document(_ARTIFACT.read_text())

    def test_source_integrity_and_record_count(self):
        self.assertEqual(self.document["source"], {
            "url": "https://www.math.uniri.hr/~mmaksimovic/srg45.txt",
            "bytes": 1470604,
            "sha256": "a53366d919f5d29876d07031a65c1bd722d05e88b59b39a04e9988c9162ebec0",
            "format": "GAP_grape_records",
        })
        self.assertEqual(self.document["records"], 288)

    def test_support_census_dependency_is_hash_bound(self):
        self.assertEqual(self.document["support_census_dependency"], {
            "relative_path": "r55/data/involution_f5_support_census.json",
            "bytes": 503066,
            "sha256":
                "520e2cb453cc2efee7925636af018316acca3e375a8107886c5e9c90c3727aaa",
            "disposition":
                "F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS",
            "representatives": 844,
        })

    def test_every_record_has_both_ramsey_obstructions(self):
        self.assertEqual(self.document["validation"], {
            "srg_45_22_10_11": 288,
            "distinct_labelled_adjacency_records": 288,
            "specified_order3_automorphism": 288,
            "order3_fixed_points": 9,
            "order3_cycles": 12,
            "specified_involution_automorphism": 288,
            "involution_fixed_points": 5,
            "involution_transpositions": 20,
            "fixed_graph_C5": 288,
            "balanced_involution_support": 288,
            "signed_seidel_completion_equations": 288,
            "labelled_involution_support_signatures": 38,
            "dihedral_involution_support_types": 26,
            "support_census_membership_verified": True,
            "support_types_in_balanced_census": 26,
            "involution_coverage_sha256":
                "6f8533cf6d35cf0e967da50be94e63b9ebbc9e67c10eb3189fa326dcdae96ac9",
            "with_K5": 288,
            "with_I5": 288,
            "with_both": 288,
            "witness_coverage_sha256":
                "18fa4784de5d2cf61b20b340c25bf66b511c6677119c506177d50d5b0106fec2",
        })

    def test_scope_is_empirical_only(self):
        self.assertEqual(self.document["scope"], {
            "role":
                "source_verified_positive_controls_for_order3_and_fixed5_"
                "involution_signed_completion",
            "catalog_completeness_claimed": False,
            "novelty_claimed": False,
            "general_ramsey_bound_claimed": False,
        })

    def test_duplicate_json_members_are_rejected(self):
        payload = _ARTIFACT.read_text()
        needle = '"general_ramsey_bound_claimed": false'
        replacement = (
            '"general_ramsey_bound_claimed": true, '
            '"general_ramsey_bound_claimed": false'
        )
        self.assertIn(needle, payload)
        with self.assertRaises(validation.ValidationViolation):
            validation.load_document(payload.replace(needle, replacement, 1))

    def test_exact_tree_comparison_rejects_json_type_confusion(self):
        false_as_zero = copy.deepcopy(self.document)
        false_as_zero["scope"]["general_ramsey_bound_claimed"] = 0
        count_as_float = copy.deepcopy(self.document)
        count_as_float["records"] = 288.0
        self.assertFalse(validation._exact_tree_equal(
            self.document, false_as_zero
        ))
        self.assertFalse(validation._exact_tree_equal(
            self.document, count_as_float
        ))


if __name__ == "__main__":
    unittest.main()
