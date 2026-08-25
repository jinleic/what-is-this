"""Tests for the exact involution frontier of ``srg(45,22,10,11)``.

The theorem is SRG-only: an involution can fix only 1, 5, 9, or 13
vertices.  It does not exclude those four cases and changes no bound on
``R(5,5)``.
"""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import check_involution_srg_frontier as check  # noqa: E402
import involution_srg_frontier as frontier  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_ARTIFACT = _ROOT / "data" / "involution_srg_frontier.json"


class TestSpectralReduction(unittest.TestCase):
    def test_fixed_count_candidates_are_exact(self):
        record = frontier.fixed_count_reduction()
        self.assertEqual(
            record["spectral_candidates"],
            [1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41],
        )
        self.assertEqual(record["half_incidence_candidates"], [1, 5, 9, 13, 17, 21])
        self.assertEqual(record["surviving_fixed_counts"], [1, 5, 9, 13])

    def test_anti_invariant_row_has_eleven_half_incidence_blocks(self):
        record = frontier.anti_invariant_identity()
        self.assertEqual(record["matrix_equation"], "C^2+C=11I")
        self.assertEqual(record["off_diagonal_values"], [-1, 0, 1])
        self.assertEqual(record["half_incidence_blocks_per_row"], 11)
        self.assertEqual(record["edge_pair_orbits"], "c/2")


class TestSmallAntiInvariantEnumeration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = frontier.enumerate_small_anti_obstructions()

    def test_fixed_twenty_one_and_seventeen_are_closed(self):
        self.assertEqual(self.record["closed_fixed_counts"], [17, 21])
        self.assertEqual(self.record["c12"]["solutions"], 0)
        self.assertEqual(self.record["c14"]["solutions"], 0)


    def test_zero_neighborhood_parity_already_closes_both_cases(self):
        self.assertEqual(
            self.record["c12"]["zero_neighborhood_parity_violations"], 30
        )
        self.assertTrue(all(
            item["zero_neighborhood_parity_violations"] > 0
            for item in self.record["c14"]["zero_graph_shapes"]
        ))


    def test_c14_zero_graph_classification_is_complete(self):
        records = self.record["c14"]["zero_graph_shapes"]
        self.assertEqual(len(records), 12)
        self.assertEqual({item["solutions"] for item in records}, {0})
        self.assertEqual(
            {item["cross_vertices_per_part"] for item in records},
            {0, 2, 3, 4, 7},
        )
        self.assertTrue(all(item["search_sha256"] for item in records))


class TestRamseyInterface(unittest.TestCase):
    def test_fixed_five_graph_is_uniquely_c5(self):
        record = frontier.fixed_five_graph_enumeration()
        self.assertEqual(record["labelled_solutions"], 12)
        self.assertEqual(record["isomorphism_types"], 1)
        self.assertEqual(record["unique_type"], "C5")

    def test_fixed_five_support_design_survives(self):
        record = frontier.fixed_five_support_design()
        self.assertEqual(len(record["edge_pair_support_masks"]), 10)
        self.assertEqual(len(record["nonedge_pair_support_masks"]), 10)
        self.assertEqual(record["replication"], [10] * 5)
        self.assertEqual(record["pair_multiplicity"], 5)
        self.assertEqual(record["disposition"], "EXACT_SUPPORT_DESIGN_SURVIVES")


class TestArtifactAndChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = frontier.run_analysis()

    def test_disposition_and_scope(self):
        self.assertEqual(
            self.document["disposition"],
            "INVOLUTION_FIXED_COUNTS_REDUCED_TO_1_5_9_13_FOR_SRG",
        )
        self.assertEqual(self.document["claim"], {
            "graph_family": "strongly_regular",
            "parameters": [45, 22, 10, 11],
            "involution_fixed_counts": [1, 5, 9, 13],
            "ramsey_hypothesis_used": False,
            "general_ramsey_bound_claimed": False,
        })
        self.assertEqual(
            self.document["ramsey_frontier"]["disposition"],
            "OPEN_RAMSEY_COMPLETION_FOR_FIXED_1_5_9_13",
        )

    def test_committed_artifact_matches_producer(self):
        self.assertEqual(json.loads(_ARTIFACT.read_text()), self.document)

    def test_independent_checker_accepts_control(self):
        self.assertEqual(check.verify_document(self.document), {
            "closed_fixed_counts": [17, 21],
            "surviving_fixed_counts": [1, 5, 9, 13],
            "c14_zero_graph_shapes": 12,
            "ramsey_status": "OPEN_RAMSEY_COMPLETION_FOR_FIXED_1_5_9_13",
        })


    def test_search_hashes_bind_their_domains(self):
        records = self.document["anti_invariant_search"]["c14"][
            "zero_graph_shapes"
        ]
        self.assertEqual(len({record["search_sha256"] for record in records}), 12)
        self.assertTrue(all(
            record["search_domain"]["zero_edges"] == record["zero_edges"]
            for record in records
        ))

    def test_exact_json_types_are_required(self):
        mutations = []
        bad = copy.deepcopy(self.document)
        bad["schema_version"] = True
        mutations.append(bad)
        bad = copy.deepcopy(self.document)
        bad["schema_version"] = 1.0
        mutations.append(bad)
        bad = copy.deepcopy(self.document)
        bad["anti_invariant_search"]["c12"]["solutions"] = False
        mutations.append(bad)
        for mutation in mutations:
            with self.assertRaises(check.CheckViolation):
                check.verify_document(mutation)

    def test_duplicate_json_members_are_rejected(self):
        payload = _ARTIFACT.read_text()
        needle = '"general_ramsey_bound_claimed": false'
        duplicate = (
            '"general_ramsey_bound_claimed": true, '
            '"general_ramsey_bound_claimed": false'
        )
        self.assertIn(needle, payload)
        with self.assertRaises(check.CheckViolation):
            check.load_document(payload.replace(needle, duplicate, 1))

    def test_mutating_returned_documents_cannot_poison_caches(self):
        produced = frontier.run_analysis()
        produced["anti_invariant_search"]["c14"]["solutions"] = 1
        produced["ramsey_frontier"]["fixed_5"]["support_design"][
            "pair_multiplicity"
        ] = 4
        pristine = frontier.run_analysis()
        self.assertEqual(pristine["anti_invariant_search"]["c14"]["solutions"], 0)
        self.assertEqual(
            pristine["ramsey_frontier"]["fixed_5"]["support_design"][
                "pair_multiplicity"
            ],
            5,
        )
        expected = check._expected_document()
        expected["schema_version"] = False
        check.verify_document(pristine)


    def test_corrupted_c14_solution_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["anti_invariant_search"]["c14"]["zero_graph_shapes"][0]["solutions"] = 1
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_search_digest_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["anti_invariant_search"]["c14"]["zero_graph_shapes"][0][
            "search_sha256"
        ] = "0" * 64
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_surviving_fixed_count_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["claim"]["involution_fixed_counts"].append(17)
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_fixed_five_design_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["ramsey_frontier"]["fixed_5"]["support_design"][
            "pair_multiplicity"
        ] = 4
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_false_general_ramsey_claim_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["claim"]["general_ramsey_bound_claimed"] = True
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)


if __name__ == "__main__":
    unittest.main()
