"""Tests for the exact fixed-five involution balanced-support census.

The census is support-only.  Signed SRG completion and every Ramsey-number
claim remain open.
"""

from __future__ import annotations

import copy
import json
import itertools
import sys
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import check_involution_f5_support_census as check  # noqa: E402
import involution_f5_support_census as census  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_ARTIFACT = _ROOT / "data" / "involution_f5_support_census.json"


class TestWalshReduction(unittest.TestCase):
    def test_rotation_design_partner_is_reconstructed(self):
        edge_supports = (3, 6, 12, 15, 17, 23, 24, 27, 29, 30)
        expected_nonedge = (1, 2, 4, 8, 11, 13, 16, 21, 22, 26)
        counts = census.counts_from_supports(edge_supports, census.EVEN_MASKS)
        partner = census.walsh_partner_odd_counts(counts)
        self.assertIsNotNone(partner)
        self.assertEqual(
            census.supports_from_counts(partner, census.ODD_MASKS),
            expected_nonedge,
        )

    def test_invalid_even_multiset_has_no_partner(self):
        counts = [0] * 16
        counts[0] = 10
        self.assertIsNone(census.walsh_partner_odd_counts(tuple(counts)))


class TestExactCensus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = census.support_census()

    def test_labelled_and_dihedral_counts(self):
        self.assertEqual(self.record["balanced_labelled_designs"], 7872)
        self.assertEqual(self.record["dihedral_orbits"], 844)
        self.assertEqual(
            self.record["orbit_size_histogram"],
            {"1": 2, "5": 110, "10": 732},
        )
        self.assertEqual(len(self.record["representatives"]), 844)

    def test_burnside_audit(self):
        self.assertEqual(
            self.record["burnside_fixed_designs"],
            {
                "identity": 7872,
                "nonidentity_rotations_each": 2,
                "reflections_each": 112,
            },
        )
        fixed_sum = 7872 + 4 * 2 + 5 * 112
        self.assertEqual(fixed_sum // 10, 844)

    def test_every_representative_is_balanced(self):
        for representative in self.record["representatives"]:
            check.check_representative(representative)

    def test_explicit_positive_control_is_present(self):
        control = self.record["positive_control"]
        representative = self.record["representatives"][
            control["canonical_representative_index"]
        ]
        self.assertEqual(control["dihedral_orbit_size"], 1)
        self.assertEqual(representative["orbit_size"], 1)


class TestSignedCompletionBoundary(unittest.TestCase):
    def test_exact_equations_are_recorded_without_claiming_completion(self):
        model = census.signed_completion_model()
        self.assertEqual(model["invariant_square"], "W^2+2R^T R=45I-2J")
        self.assertEqual(model["anti_invariant_square"], "T^2=45I")
        self.assertEqual(model["off_diagonal_compatibility"], "W_ij^2+T_ij^2=4")
        self.assertEqual(
            model["anti_invariant_diagonal"],
            "T_ii=+1 on edge pairs and -1 on nonedge pairs",
        )
        self.assertIn("c_ij=-T_ij/2", model["relation_recovery"])
        self.assertEqual(
            model["status"],
            "OPEN_SIGNED_COMPLETION_FOR_844_BALANCED_SUPPORT_ORBITS",
        )
        self.assertFalse(model["ramsey_bound_claimed"])

    def test_true_anti_seidel_sign_recovers_parallel_and_crossed_bits(self):
        for parallel, crossed in itertools.product((0, 1), repeat=2):
            invariant = 2 - 2 * (parallel + crossed)
            anti = 2 * (crossed - parallel)
            block_sum = 1 - invariant // 2
            block_difference = -anti // 2
            self.assertEqual(
                (
                    (block_sum + block_difference) // 2,
                    (block_sum - block_difference) // 2,
                ),
                (parallel, crossed),
            )


class TestArtifactAndIndependentChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = census.run_analysis()

    def test_committed_artifact_matches_producer(self):
        self.assertEqual(json.loads(_ARTIFACT.read_text()), self.document)

    def test_independent_checker_accepts_control(self):
        self.assertEqual(
            check.verify_document(self.document),
            {
                "balanced_labelled_designs": 7872,
                "dihedral_orbits": 844,
                "signed_completion_status":
                    "OPEN_SIGNED_COMPLETION_FOR_844_BALANCED_SUPPORT_ORBITS",
                "ramsey_bound_claimed": False,
            },
        )

    def test_exact_json_types_are_required(self):
        for path, value in (
            (("schema_version",), True),
            (("schema_version",), 1.0),
            (("support_census", "dihedral_orbits"), False),
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

    def test_corrupted_count_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["support_census"]["balanced_labelled_designs"] = 7871
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_representative_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["support_census"]["representatives"][0]["edge_even_counts"][0] += 1
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_false_ramsey_claim_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["claim"]["general_ramsey_bound_claimed"] = True
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_mutation_cannot_poison_cached_documents(self):
        first = census.run_analysis()
        first["support_census"]["dihedral_orbits"] = 1
        first["signed_completion"]["ramsey_bound_claimed"] = True
        second = census.run_analysis()
        self.assertEqual(second["support_census"]["dihedral_orbits"], 844)
        self.assertFalse(second["signed_completion"]["ramsey_bound_claimed"])
        check.verify_document(second)


if __name__ == "__main__":
    unittest.main()
