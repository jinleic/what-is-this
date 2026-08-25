"""Tests for the exact R(3,3) support filter in the fixed-five branch.

The filter is necessary for a Ramsey-good completion.  It does not establish
signed completion or any general R(5,5) bound.
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

import check_involution_f5_ramsey_filter as check  # noqa: E402
import involution_f5_ramsey_filter as frontier  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_ARTIFACT = _ROOT / "data" / "involution_f5_ramsey_filter.json"
_SUPPORT_ARTIFACT = _ROOT / "data" / "involution_f5_support_census.json"


class TestForcedRamseyIntersection(unittest.TestCase):
    def test_opposite_fixed_windows_are_declared(self):
        self.assertEqual(
            frontier.OPPOSITE_WINDOWS,
            (
                ((0, 1), (2, 4)),
                ((1, 2), (0, 3)),
                ((2, 3), (1, 4)),
                ((3, 4), (0, 2)),
                ((0, 4), (1, 3)),
            ),
        )

    def test_six_vertices_force_a_clique_or_independent_set(self):
        self.assertTrue(check.verify_r33_six_vertex_lemma())

    def test_four_vertex_relation_domains_are_exact(self):
        self.assertEqual(
            frontier.four_vertex_allowed_w_values(0, 0), (-1, 0))
        self.assertEqual(
            frontier.four_vertex_allowed_w_values(0, 1), (0,))
        self.assertEqual(
            frontier.four_vertex_allowed_w_values(1, 0), (0,))
        self.assertEqual(
            frontier.four_vertex_allowed_w_values(1, 1), (0, 1))
        for left_internal in (0, 1):
            for right_internal in (0, 1):
                self.assertEqual(
                    frontier.four_vertex_allowed_w_values(
                        left_internal, right_internal),
                    check.enumerate_four_vertex_allowed_w_values(
                        left_internal, right_internal),
                )


class TestExactRamseyFilterCensus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = frontier.run_analysis()
        cls.census = cls.document["census"]

    def test_exact_rejected_and_surviving_counts(self):
        self.assertEqual(self.census["input_d5_orbits"], 844)
        self.assertEqual(self.census["input_labelled_supports"], 7872)
        self.assertEqual(self.census["rejected_d5_orbits"], 139)
        self.assertEqual(self.census["rejected_labelled_supports"], 1245)
        self.assertEqual(self.census["surviving_d5_orbits"], 705)
        self.assertEqual(self.census["surviving_labelled_supports"], 6627)
        self.assertEqual(
            self.census["surviving_orbit_size_histogram"],
            {"1": 2, "5": 81, "10": 622},
        )
        self.assertEqual(
            self.census["rejected_orbit_size_histogram"],
            {"1": 0, "5": 29, "10": 110},
        )

    def test_every_candidate_has_five_opposite_window_counts(self):
        rejected = 0
        survivors = 0
        for candidate in self.census["candidates"]:
            counts = candidate["opposite_window_multiplicities"]
            self.assertEqual(len(counts), 5)
            self.assertTrue(all(type(value) is int for value in counts))
            if candidate["status"] == "REJECTED_FORCED_R33_SIX_SET":
                rejected += 1
                self.assertEqual(max(counts), 3)
                self.assertTrue(candidate["violating_window_indices"])
                self.assertEqual(candidate["relation_restrictions"], [])
            else:
                survivors += 1
                self.assertEqual(
                    candidate["status"], "SURVIVES_R33_SUPPORT_FILTER")
                self.assertLessEqual(max(counts), 2)
                self.assertEqual(candidate["violating_window_indices"], [])
        self.assertEqual((rejected, survivors), (139, 705))

    def test_relation_restrictions_have_exact_domains(self):
        support = json.loads(_SUPPORT_ARTIFACT.read_text())
        representatives = support["support_census"]["representatives"]
        for candidate in self.census["candidates"]:
            if candidate["status"] != "SURVIVES_R33_SUPPORT_FILTER":
                continue
            representative = representatives[candidate["source_index"]]
            masks, internal = frontier.expand_representative(
                support["support_census"], representative)
            expected = frontier.r33_relation_restrictions(masks, internal)
            self.assertEqual(candidate["relation_restrictions"], expected)

    def test_filter_is_invariant_under_the_full_dihedral_action(self):
        support = json.loads(_SUPPORT_ARTIFACT.read_text())
        for representative in support["support_census"]["representatives"]:
            masks, _ = frontier.expand_representative(
                support["support_census"], representative)
            status = max(frontier.opposite_window_multiplicities(masks)) <= 2
            for permutation in frontier.DIHEDRAL_PERMUTATIONS:
                image = tuple(frontier.permute_mask(mask, permutation)
                              for mask in masks)
                self.assertEqual(
                    max(frontier.opposite_window_multiplicities(image)) <= 2,
                    status,
                )

    def test_new_frontier_does_not_claim_signed_completion(self):
        claim = self.document["claim"]
        self.assertEqual(
            claim["signed_completion_status"],
            "OPEN_SIGNED_COMPLETION_FOR_705_R33_SURVIVING_SUPPORT_ORBITS",
        )
        self.assertTrue(claim["strongly_regular_fixed_five_branch_only"])
        self.assertFalse(claim["general_ramsey_bound_claimed"])


class TestSignedConsequences(unittest.TestCase):
    def test_sign_degree_table_and_global_counts(self):
        model = frontier.signed_frontier_model()
        self.assertEqual(
            model["w_sign_degrees_by_support_size"],
            {
                "0": {"W=+2": 3, "W=-2": 5},
                "1": {"W=+2": 3, "W=-2": 5},
                "2": {"W=+2": 4, "W=-2": 4},
                "3": {"W=+2": 4, "W=-2": 4},
                "4": {"W=+2": 5, "W=-2": 3},
                "5": {"W=+2": 5, "W=-2": 3},
            },
        )
        self.assertEqual(
            model["global_unordered_relation_counts"],
            {"W=+2": 40, "W=-2": 40, "T-supported": 110},
        )
        self.assertEqual(model["w_action_on_one"], "W 1=-R^T 1")
        self.assertEqual(model["w_action_on_support_sum"], "W R^T 1=-5 1")


class TestArtifactAndIndependentChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.produced = frontier.run_analysis()
        cls.document = check.load_document(_ARTIFACT.read_bytes())

    def test_committed_artifact_matches_producer(self):
        self.assertEqual(self.document, self.produced)

    def test_independent_checker_accepts_artifact(self):
        self.assertEqual(
            check.verify_document(self.document),
            {
                "rejected_d5_orbits": 139,
                "surviving_d5_orbits": 705,
                "surviving_labelled_supports": 6627,
                "signed_completion_status":
                    "OPEN_SIGNED_COMPLETION_FOR_705_R33_SURVIVING_SUPPORT_ORBITS",
                "general_ramsey_bound_claimed": False,
            },
        )

    def test_exact_json_types_are_required(self):
        for path, value in (
            (("schema_version",), True),
            (("schema_version",), 1.0),
            (("census", "surviving_d5_orbits"), 705.0),
            (("claim", "general_ramsey_bound_claimed"), 0),
        ):
            bad = copy.deepcopy(self.document)
            target = bad
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.assertRaises(check.CheckViolation):
                check.verify_document(bad)

    def test_duplicate_json_members_are_rejected(self):
        payload = _ARTIFACT.read_text()
        needle = '"general_ramsey_bound_claimed": false'
        replacement = (
            '"general_ramsey_bound_claimed": true, '
            '"general_ramsey_bound_claimed": false'
        )
        self.assertIn(needle, payload)
        with self.assertRaises(check.CheckViolation):
            check.load_document(payload.replace(needle, replacement, 1))

    def test_dependency_hash_is_enforced(self):
        bad = copy.deepcopy(self.document)
        bad["input_dependency"]["sha256"] = "0" * 64
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_noncanonical_dependency_path_is_rejected(self):
        noncanonical = Path("/tmp/noncanonical-involution-f5-support.json")
        with self.assertRaises(frontier.FrontierViolation):
            frontier.load_support_artifact(noncanonical)
        with self.assertRaises(check.CheckViolation):
            check.verify_document(self.document, noncanonical)

    def test_corrupted_candidate_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["census"]["candidates"][0][
            "opposite_window_multiplicities"][0] += 1
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_false_completion_or_ramsey_claim_is_rejected(self):
        for key, value in (
            ("signed_completion_status", "COMPLETE"),
            ("general_ramsey_bound_claimed", True),
        ):
            bad = copy.deepcopy(self.document)
            bad["claim"][key] = value
            with self.assertRaises(check.CheckViolation):
                check.verify_document(bad)

    def test_mutation_cannot_poison_cached_document(self):
        first = frontier.run_analysis()
        first["census"]["surviving_d5_orbits"] = 1
        second = frontier.run_analysis()
        self.assertEqual(second["census"]["surviving_d5_orbits"], 705)
        check.verify_document(second)


if __name__ == "__main__":
    unittest.main()
