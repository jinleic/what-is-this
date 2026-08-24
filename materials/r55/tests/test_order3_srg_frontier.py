"""Tests for the order-three automorphism frontier of srg(45,22,10,11).

The claimed theorem is narrow but complete within its lane: a Ramsey(5,5)-good
srg(45,22,10,11) has no automorphism of order three.  It says nothing about
non-strongly-regular Ramsey graphs and changes no bound on R(5,5).
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

import check_order3_srg_frontier as check  # noqa: E402
import order3_srg_frontier as frontier  # noqa: E402

_ARTIFACT = Path(__file__).resolve().parents[1] / "data" / "order3_srg_frontier.json"


class TestFixedGraph(unittest.TestCase):
    def test_constructed_graph_is_srg_9_4_1_2(self):
        graph = frontier.construct_fixed_graph()
        self.assertEqual(frontier.srg_parameters(graph), (9, 4, 1, 2))
        self.assertEqual(len(frontier.independent_triples(graph)), 6)
        self.assertEqual(len(frontier.triangles(graph)), 6)

    def test_local_reduction_has_one_isomorphism_type(self):
        candidates = frontier.enumerate_fixed_graph_candidates()
        self.assertEqual(len(candidates), 8)
        reference = frontier.construct_fixed_graph()
        self.assertTrue(all(frontier.is_isomorphic(reference, graph)
                            for graph in candidates))

    def test_degree_one_and_seven_residue_obstructions_are_complete(self):
        self.assertEqual(frontier.degree_extreme_obstructions(), {
            "degree_one_max_common_fixed": 1,
            "nonedge_required_residue": 2,
            "possible_common_fixed_counts": [0, 1],
            "residue_hits": [],
            "degree_seven_complement_degree": 1,
        })


class TestIncidenceFrame(unittest.TestCase):
    def test_ramsey_columns_are_unique(self):
        graph = frontier.construct_fixed_graph()
        incidence = frontier.build_incidence(graph)
        self.assertEqual(frontier.row_sums(incidence), (6,) * 9)
        self.assertEqual(frontier.column_sums(incidence), (3,) * 6 + (6,) * 6)
        self.assertEqual(frontier.gram(incidence), frontier.target_gram())
        self.assertEqual(frontier.incidence_multiplicity_solutions(graph),
                         [((1,) * 6, (1,) * 6)])

    def test_columns_have_the_ramsey_meaning(self):
        graph = frontier.construct_fixed_graph()
        incidence = frontier.build_incidence(graph)
        for column in range(6):
            members = frontier.column_members(incidence, column)
            self.assertTrue(frontier.is_independent(graph, members))
        for column in range(6, 12):
            nonmembers = tuple(v for v in range(9)
                               if not incidence[v][column])
            self.assertTrue(frontier.is_clique(graph, nonmembers))


class TestQuotientReduction(unittest.TestCase):
    def test_linear_family_is_complete(self):
        system = frontier.quotient_linear_system()
        self.assertEqual(system["variables"], 78)
        self.assertEqual(system["rank"], 77)
        self.assertEqual(system["augmented_rank"], 77)
        self.assertTrue(frontier.verify_affine_family(system))

    def test_exact_quotient_candidates_are_tau_one_and_two(self):
        candidates = frontier.quotient_candidates()
        self.assertEqual([record["tau"] for record in candidates], [1, 2])
        for record in candidates:
            matrix = record["matrix"]
            self.assertTrue(frontier.verify_quotient_matrix(matrix))
            self.assertEqual(tuple(matrix[i][i] for i in range(12)),
                             (2,) * 6 + (0,) * 6)

    def test_both_candidates_force_the_same_independent_five(self):
        witness = frontier.forced_independent_five()
        self.assertEqual(witness, {
            "fixed_vertices": [2, 3, 6],
            "orbit_pair": [0, 3],
            "orbit_positions": [0, 1],
            "normalized_phase": 0,
        })
        for record in frontier.quotient_candidates():
            self.assertTrue(frontier.verify_forced_independent_five(
                record["matrix"], witness))


class TestLargeFixedCounts(unittest.TestCase):
    def test_trace_forces_half_the_nonfixed_orbits_to_be_triangles(self):
        records = frontier.large_fixed_count_obstructions()
        self.assertEqual(
            [(record["fixed_points"], record["orbit_3_cycles"],
              record["triangle_orbits"]) for record in records],
            [(15, 10, 5), (21, 8, 4), (27, 6, 3), (33, 4, 2), (39, 2, 1)],
        )

    def test_triangle_moment_excludes_every_large_fixed_count_without_ramsey(self):
        records = frontier.large_fixed_count_obstructions()
        self.assertEqual(
            [record["minimum_cauchy_gap"] for record in records],
            [16, 44, 66, 88, 110],
        )
        for record in records:
            self.assertEqual(record["fixed_support_values"], list(range(10)))
            self.assertTrue(all(gap > 0 for gap in record["cauchy_gaps"]))

    def test_moment_cut_does_not_claim_the_fixed_nine_case(self):
        boundary = frontier.triangle_orbit_moment(12, 4)
        self.assertEqual(boundary, {
            "fixed_support": 4,
            "other_orbits": 11,
            "off_orbit_degree_sum": 16,
            "off_orbit_square_sum": 26,
            "cauchy_gap": -30,
        })
        self.assertLessEqual(boundary["cauchy_gap"], 0)


class TestArtifactAndChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = frontier.run_analysis()

    def test_disposition_and_scope(self):
        self.assertEqual(self.document["schema_version"], 2)
        self.assertEqual(
            self.document["disposition"],
            "ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG",
        )
        self.assertEqual(self.document["claim"], {
            "graph_family": "ramsey_good_strongly_regular",
            "parameters": [45, 22, 10, 11],
            "automorphism_order": 3,
            "conclusion": "no_automorphism_of_order_three",
            "excluded_spectral_fixed_point_counts": [3, 9, 15, 21, 27, 33, 39],
            "general_ramsey_bound_claimed": False,
        })
        self.assertEqual(
            self.document["fixed_15_plus"]["disposition"],
            "IMPOSSIBLE_TRIANGLE_ORBIT_MOMENT_FOR_ANY_SRG",
        )
        self.assertFalse(self.document["fixed_15_plus"]["ramsey_hypothesis_used"])

    def test_committed_artifact_matches_producer(self):
        self.assertEqual(json.loads(_ARTIFACT.read_text()), self.document)

    def test_independent_checker_accepts_control(self):
        summary = check.verify_document(self.document)
        self.assertEqual(summary["quotient_candidates"], 2)
        self.assertEqual(summary["forced_independent_sets"], 2)
        self.assertEqual(summary["large_fixed_counts"], 5)
        self.assertEqual(summary["minimum_cauchy_gap"], 16)

    def test_corrupted_tau_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["fixed_9"]["quotient_candidates"][0]["tau"] = 0
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_large_fixed_count_moment_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["fixed_15_plus"]["records"][0]["minimum_cauchy_gap"] = 15
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_large_fixed_count_trace_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["fixed_15_plus"]["records"][0]["triangle_orbits"] = 4
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_degree_residue_obstruction_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["fixed_9"]["degree_extreme_obstructions"]["residue_hits"] = [1]
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)


    def test_false_general_bound_claim_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["claim"]["general_ramsey_bound_claimed"] = True
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)

    def test_corrupted_witness_is_rejected(self):
        bad = copy.deepcopy(self.document)
        bad["fixed_9"]["forced_independent_five"]["fixed_vertices"] = [0, 1, 2]
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad)


if __name__ == "__main__":
    unittest.main()
