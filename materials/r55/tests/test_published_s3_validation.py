"""Tests for the independent published S3 graph-data replay."""

from __future__ import annotations

import json
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


class TestPublishedArtifact(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(_ARTIFACT.read_text())

    def test_source_integrity_and_record_count(self):
        self.assertEqual(self.document["source"], {
            "url": "https://www.math.uniri.hr/~mmaksimovic/srg45.txt",
            "bytes": 1470604,
            "sha256": "a53366d919f5d29876d07031a65c1bd722d05e88b59b39a04e9988c9162ebec0",
            "format": "GAP_grape_records",
        })
        self.assertEqual(self.document["records"], 288)

    def test_every_record_has_both_ramsey_obstructions(self):
        self.assertEqual(self.document["validation"], {
            "srg_45_22_10_11": 288,
            "distinct_labelled_adjacency_records": 288,
            "specified_order3_automorphism": 288,
            "order3_fixed_points": 9,
            "order3_cycles": 12,
            "with_K5": 288,
            "with_I5": 288,
            "with_both": 288,
            "witness_coverage_sha256":
                "18fa4784de5d2cf61b20b340c25bf66b511c6677119c506177d50d5b0106fec2",
        })

    def test_scope_is_empirical_only(self):
        self.assertEqual(self.document["scope"], {
            "role": "independent_empirical_corroboration_of_order3_fixed9_exclusion",
            "catalog_completeness_claimed": False,
            "novelty_claimed": False,
            "general_ramsey_bound_claimed": False,
        })


if __name__ == "__main__":
    unittest.main()
