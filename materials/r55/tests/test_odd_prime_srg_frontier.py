"""Tests for the odd-prime automorphism frontier of srg(45,22,10,11).

The theorem is confined to Ramsey-good strongly regular graphs: their
full automorphism group is a 2-group.  It changes no bound on R(5,5).
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

import check_odd_prime_srg_frontier as check  # noqa: E402
import odd_prime_srg_frontier as frontier  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_ARTIFACT = _ROOT / "data" / "odd_prime_srg_frontier.json"
_ORDER3_ARTIFACT = _ROOT / "data" / "order3_srg_frontier.json"


class TestOrderFive(unittest.TestCase):
    def test_trace_forces_every_moving_orbit_to_be_c5(self):
        self.assertEqual(frontier.order5_trace_record(), {
            "quotient_trace": "2*c",
            "internal_degree_counts_equation": "n0=n4",
            "ramsey_allowed_internal_degrees": [2],
            "moving_orbit_graph": "C5",
        })

    def test_moment_excludes_every_fixed_count_at_least_ten_without_ramsey(self):
        records = frontier.order5_large_fixed_records()
        self.assertEqual(
            [(record["fixed_points"], record["orbit_5_cycles"],
              record["minimum_cauchy_gap"]) for record in records],
            [(10, 7, 15), (15, 6, 44), (20, 5, 60), (25, 4, 64),
             (30, 3, 64), (35, 2, 64), (40, 1, 64)],
        )
        for record in records:
            self.assertEqual(
                [item["internal_degree"] for item in record["internal_degree_records"]],
                [0, 2, 4],
            )
            self.assertTrue(all(
                gap > 0
                for item in record["internal_degree_records"]
                for gap in item["cauchy_gaps"]
            ))

    def test_fixed_free_quotient_search_is_complete_and_empty(self):
        records = frontier.enumerate_order5_free_quotients()
        self.assertEqual([record["solutions"] for record in records], [0, 0, 0])
        self.assertEqual(
            [record["tested_by_row"] for record in records],
            [
                {"1": 735, "2": 15120, "3": 9680},
                {"1": 735, "2": 18420, "3": 6680},
                {"1": 105, "2": 2520, "3": 2760},
            ],
        )
        self.assertEqual(
            [record["surviving_by_row"] for record in records],
            [
                {"1": 72, "2": 132},
                {"1": 64, "2": 96},
                {"1": 8, "2": 48},
            ],
        )

    def test_fixed_five_case_has_a_solver_free_ramsey_bridge(self):
        record = frontier.order5_fixed_five_obstruction()
        self.assertEqual(record["fixed_graph"], "C5")
        self.assertEqual(record["incidence_row_size"], 4)
        self.assertEqual(record["incidence_pair_multiplicity"], 2)
        self.assertEqual(record["moving_vertices_in_bridge"], 10)
        self.assertEqual(record["ramsey_numbers"], {"R(3,4)": 9, "R(4,3)": 9})
        self.assertTrue(record["all_triple_multiplicities_one_is_impossible"])


class TestOrderSeven(unittest.TestCase):
    def test_fixed_point_candidates_are_all_excluded(self):
        record = frontier.order7_obstructions()
        self.assertEqual(record["spectral_candidates"], [
            {"fixed_points": 3, "orbit_7_cycles": 6},
            {"fixed_points": 17, "orbit_7_cycles": 4},
            {"fixed_points": 31, "orbit_7_cycles": 2},
        ])
        self.assertEqual(record["fixed_17"]["maximum_support_family_size"], 6)
        self.assertEqual(record["fixed_31"]["support_group_sizes"], [16, 14])
        self.assertEqual(record["fixed_31"]["rayleigh_fourfold_gap"], -11130)


class TestOrderEleven(unittest.TestCase):
    def test_fixed_one_quotient_is_unique_up_to_swap(self):
        matrix = frontier.order11_canonical_quotient()
        self.assertEqual(matrix, [
            [4, 6, 4, 7],
            [6, 4, 7, 4],
            [4, 7, 6, 5],
            [7, 4, 5, 6],
        ])
        self.assertTrue(frontier.verify_order11_quotient(matrix))

    def test_supported_two_orbit_block_is_exhausted(self):
        coverage = frontier.order11_block_coverage()
        self.assertEqual(coverage["configurations"], 46200)
        self.assertEqual(coverage["k4_primary_witnesses"], 44000)
        self.assertEqual(coverage["i5_fallback_witnesses"], 2200)
        self.assertEqual(
            coverage["coverage_sha256"],
            "b888c63a7b77419c3e6df970acd3759391da073a3357f3318f11f0611d5538a6",
        )

    def test_fixed_twenty_three_support_pigeonhole_forces_i5(self):
        record = frontier.order11_fixed_23_obstruction()
        self.assertEqual(record["largest_single_support_group_lower_bound"], 11)
        self.assertEqual(record["forced_independent_set_lower_bound"], 11)


class TestLargePrimes(unittest.TestCase):
    def test_every_prime_at_least_thirteen_is_excluded(self):
        record = frontier.large_prime_obstructions()
        self.assertEqual(
            [(item["prime"], item["fixed_points"], item["fixed_degree"])
             for item in record["handshake_cases"]],
            [(13, 19, 9), (17, 11, 5), (19, 7, 3)],
        )
        self.assertEqual(record["cycle_parity_cutoff_prime"], 23)
        self.assertEqual(record["largest_relevant_prime"], 43)


class TestArtifactAndChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.order3_document = json.loads(_ORDER3_ARTIFACT.read_text())
        cls.document = frontier.run_analysis(_ORDER3_ARTIFACT)

    def test_disposition_and_scope(self):
        self.assertEqual(
            self.document["disposition"],
            "AUTOMORPHISM_GROUP_IS_A_2_GROUP_FOR_RAMSEY_SRG",
        )
        self.assertEqual(self.document["claim"], {
            "graph_family": "ramsey_good_strongly_regular",
            "parameters": [45, 22, 10, 11],
            "conclusion": "automorphism_group_order_is_a_power_of_two",
            "excluded_odd_prime_orders": [3, 5, 7, 11],
            "excluded_odd_prime_orders_at_least": 13,
            "general_ramsey_bound_claimed": False,
        })
        self.assertEqual(
            self.document["order_5"]["fixed_10_plus"]["disposition"],
            "IMPOSSIBLE_ORBIT_DIAGONAL_MOMENT_FOR_ANY_SRG",
        )
        self.assertFalse(
            self.document["order_5"]["fixed_10_plus"]["ramsey_hypothesis_used"])

    def test_committed_artifact_matches_producer(self):
        self.assertEqual(json.loads(_ARTIFACT.read_text()), self.document)

    def test_independent_checker_accepts_control(self):
        summary = check.verify_document(self.document, self.order3_document)
        self.assertEqual(summary, {
            "excluded_odd_primes": "all",
            "order5_free_quotient_solutions": 0,
            "order11_configurations": 46200,
            "automorphism_group": "2-group",
        })

    def test_checker_rederives_proof_critical_finite_bounds(self):
        self.assertEqual(check.independent_proof_summary(), {
            "order5_large_minimum_gap": 15,
            "order5_fixed5_block_designs": 60,
            "order5_fixed5_bridged_designs": 60,
            "order7_fixed17_maximum_support_family": 6,
            "order11_fixed23_support_group_lower_bound": 11,
            "prime13_plus_excluded": [13, 17, 19, 23, 29, 31, 37, 41, 43],
        })

    def test_corrupted_order_five_search_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["order_5"]["fixed_0"]["matrix_search"][0]["solutions"] = 1
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, self.order3_document)

    def test_corrupted_order_five_moment_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["order_5"]["fixed_10_plus"]["records"][0][
            "internal_degree_records"][0]["minimum_cauchy_gap"] = 38
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, self.order3_document)

    def test_corrupted_order_seven_rayleigh_gap_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["order_7"]["fixed_31"]["rayleigh_fourfold_gap"] = -11129
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, self.order3_document)

    def test_corrupted_order_eleven_coverage_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["order_11"]["fixed_1"]["block_coverage"]["configurations"] = 46199
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, self.order3_document)

    def test_false_general_bound_claim_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["claim"]["general_ramsey_bound_claimed"] = True
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, self.order3_document)


if __name__ == "__main__":
    unittest.main()
