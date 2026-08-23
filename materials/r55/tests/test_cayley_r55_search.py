"""Exact tests for the exhaustive order-45 Cayley Ramsey search.

The theorem under test is deliberately narrow: no undirected Cayley graph on
45 vertices avoids both a clique and an independent set of order five.  It is
not an upper bound on R(5,5), because a Ramsey graph need not be Cayley.
"""

from __future__ import annotations

import copy
import itertools
import json
import sys
import unittest
from functools import lru_cache
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import cayley_r55_search as search  # noqa: E402
import check_cayley_r55 as check  # noqa: E402

_ARTIFACT = Path(__file__).resolve().parents[1] / "data" / "cayley_r55_45.json"


EXPECTED = {
    "Z45": {
        0: (1, 1),
        2: (22, 22),
        4: (231, 230),
        6: (1540, 1520),
        8: (7315, 7102),
        10: (26334, 24788),
        12: (74613, 66066),
        14: (170544, 133366),
        16: (319770, 194744),
        18: (497420, 186764),
        20: (646646, 101549),
        22: (352716, 13457),
    },
    "Z15xZ3": {
        0: (1, 1),
        2: (22, 22),
        4: (231, 230),
        6: (1540, 1520),
        8: (7315, 7084),
        10: (26334, 24488),
        12: (74613, 63864),
        14: (170544, 124036),
        16: (319770, 170222),
        18: (497420, 149264),
        20: (646646, 72668),
        22: (352716, 8820),
    },
}

EXPECTED_DIGESTS = {
    "Z45": "ad52c249c7019a228214629d1adc28b309fe9e859ca99a9575a3d3dcee491665",
    "Z15xZ3": "856f050fc1f7728d299d7202c8d8d6a104fc319cd1ae54b729647bc0063d136d",
}

EXPECTED_CNF = {
    "Z45": (12749, 25498,
            "c84ade0feee37ae76bf76e4b1f8c1ad3f31dc32c36bfc43bea87b3025e766382"),
    "Z15xZ3": (12650, 25300,
               "2337e59cbc005447125d7455534ab7db955ac7dbdcf0f3374107c53a5302f004"),
}


@lru_cache(maxsize=1)
def _document():
    return search.run_search()


def _brute_has_clique(group, connection, size):
    vertices = range(group.order)
    for subset in itertools.combinations(vertices, size):
        if all(group.sub[u][v] in connection
               for u, v in itertools.combinations(subset, 2)):
            return True
    return False


def _connection_sets(group):
    pairs = group.inverse_pairs
    for chosen in range(1 << len(pairs)):
        yield frozenset(
            x
            for bit, pair in enumerate(pairs)
            if chosen >> bit & 1
            for x in pair
        )


class TestGroupClassification(unittest.TestCase):
    def test_there_are_exactly_two_group_models(self):
        groups = search.groups_of_order_45()
        self.assertEqual([(g.name, g.moduli) for g in groups],
                         [("Z45", (45,)), ("Z15xZ3", (15, 3))])
        for group in groups:
            self.assertEqual(group.order, 45)
            self.assertEqual(len(group.inverse_pairs), 22)
            self.assertEqual(
                {x for pair in group.inverse_pairs for x in pair},
                set(range(1, 45)),
            )
            for x in range(45):
                self.assertEqual(group.add[x][group.neg[x]], 0)
                self.assertEqual(group.neg[group.neg[x]], x)

    def test_complement_representatives_are_exact(self):
        for group in search.groups_of_order_45():
            counts = {degree: 0 for degree in range(0, 23, 2)}
            total = 0
            for low, high in search.complement_representatives(group):
                self.assertEqual(low & high, 0)
                self.assertEqual(low | high, group.nonzero_mask)
                self.assertLessEqual(low.bit_count(), high.bit_count())
                counts[low.bit_count()] += 1
                total += 1
            self.assertEqual(total, 1 << 21)
            for degree in range(0, 22, 2):
                self.assertEqual(counts[degree],
                                 search.comb(22, degree // 2))
            self.assertEqual(counts[22], search.comb(22, 11) // 2)


class TestTranslatedCliqueCriterion(unittest.TestCase):
    def test_fast_criterion_matches_full_graph_brute_force(self):
        # Every inverse-closed connection set in four small odd abelian groups.
        for moduli in ((5,), (7,), (9,), (3, 3)):
            group = search.make_group("x", moduli)
            for connection in _connection_sets(group):
                mask = sum(1 << x for x in connection)
                witness = search.first_cayley_k5(group, mask)
                brute = _brute_has_clique(group, connection, 5)
                self.assertEqual(witness is not None, brute,
                                 (moduli, sorted(connection)))
                if witness is not None:
                    self.assertEqual(witness[0], 0)
                    self.assertEqual(len(witness), 5)
                    self.assertEqual(len(set(witness)), 5)
                    self.assertTrue(all(
                        group.sub[u][v] in connection
                        for u, v in itertools.combinations(witness, 2)
                    ))

    def test_disjoint_kernel_returns_the_same_canonical_witness(self):
        for moduli in ((5,), (7,), (9,), (3, 3)):
            producer_group = search.make_group("x", moduli)
            checker_group = check._audit_group("x", moduli)
            for connection in _connection_sets(producer_group):
                mask = sum(1 << x for x in connection)
                self.assertEqual(
                    search.first_cayley_k5(producer_group, mask),
                    check._first_k5(checker_group, mask),
                )

    def test_compact_cnf_is_equivalent_on_every_z9_connection_set(self):
        group = search.make_group("Z9", (9,))
        patterns = search.ramsey_cnf_patterns(group)
        clauses = search.ramsey_cnf_clauses(patterns)
        owner = {x: index + 1
                 for index, pair in enumerate(group.inverse_pairs)
                 for x in pair}
        for connection in _connection_sets(group):
            assignment = {owner[x] for x in connection}
            cnf_satisfied = all(any(
                (literal > 0 and literal in assignment)
                or (literal < 0 and -literal not in assignment)
                for literal in clause
            ) for clause in clauses)
            mask = sum(1 << x for x in connection)
            complement = group.nonzero_mask ^ mask
            ramsey = (
                search.first_cayley_k5(group, mask) is None
                and search.first_cayley_k5(group, complement) is None
            )
            self.assertEqual(cnf_satisfied, ramsey)

    def test_positive_and_negative_paths(self):
        group = search.make_group("Z9", (9,))
        complete = group.nonzero_mask
        empty = 0
        self.assertIsNotNone(search.first_cayley_k5(group, complete))
        self.assertIsNone(search.first_cayley_k5(group, empty))
        self.assertIsNone(search.first_cayley_k5(group, 0b110))


class TestFullExhaustion(unittest.TestCase):
    def test_every_cayley_complement_class_is_rejected(self):
        document = _document()
        self.assertEqual(document["schema_version"], 2)
        self.assertEqual(document["disposition"],
                         "NO_CAYLEY_RAMSEY_5_5_45")
        self.assertEqual(document["total_complement_classes"], 1 << 22)
        self.assertEqual(document["total_connection_sets"], 1 << 23)
        self.assertEqual(document["ramsey_graphs"], 0)
        self.assertEqual(document["claim"], {
            "graph_family": "undirected_cayley",
            "order": 45,
            "forbidden_clique_order": 5,
            "forbidden_independent_set_order": 5,
            "general_ramsey_bound_claimed": False,
        })

        observed = {}
        for result in document["groups"]:
            observed[result["name"]] = {
                row["minimum_degree"]: (
                    row["complement_classes"],
                    row["low_side_k5_free"],
                )
                for row in result["degree_classes"]
            }
            self.assertEqual(result["complement_classes"], 1 << 21)
            self.assertEqual(result["ramsey_graphs"], 0)
            self.assertEqual(result["conference_pds"], 0)
            self.assertEqual(
                (result["normalized_k5_patterns"], result["ramsey_cnf_clauses"],
                 result["ramsey_cnf_sha256"]),
                EXPECTED_CNF[result["name"]],
            )
        self.assertEqual(observed, EXPECTED)

    def test_committed_artifact_matches_producer_exactly(self):
        self.assertEqual(json.loads(_ARTIFACT.read_text()), _document())


    def test_every_rejection_carries_a_valid_witness_commitment(self):
        document = _document()
        for group in document["groups"]:
            digest = group["coverage_sha256"]
            self.assertEqual(digest, EXPECTED_DIGESTS[group["name"]])
            self.assertEqual(group["witnessed_classes"], 1 << 21)


class TestIndependentChecker(unittest.TestCase):
    def test_compact_truth_table_accepts_control(self):
        summary = check.verify_document(_document(), resweep=False)
        self.assertEqual(summary["complement_classes"], 1 << 22)
        self.assertEqual(summary["ramsey_graphs"], 0)
        self.assertEqual(summary["truth_table_survivors"], 0)
        self.assertFalse(summary["witness_replay_verified"])

    def test_truth_table_kernel_against_literal_enumeration(self):
        clauses = ((1, 2), (-1, 2), (1, -2))
        survivors, used = check._truth_table_survivors(clauses, 2)
        brute = 0
        for assignment in itertools.product((False, True), repeat=2):
            if all(any(
                assignment[literal - 1] if literal > 0
                else not assignment[-literal - 1]
                for literal in clause
            ) for clause in clauses):
                brute += 1
        self.assertEqual((survivors, used), (brute, len(clauses)))
        self.assertEqual(brute, 1)

        unsat = clauses + ((-1, -2),)
        self.assertEqual(check._truth_table_survivors(unsat, 2)[0], 0)

    def test_mutated_cnf_commitment_is_rejected(self):
        bad = copy.deepcopy(_document())
        bad["groups"][0]["ramsey_cnf_sha256"] = "0" * 64
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, resweep=False)

    def test_false_general_ramsey_claim_is_rejected(self):
        bad = copy.deepcopy(_document())
        bad["claim"]["general_ramsey_bound_claimed"] = True
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, resweep=False)

    def test_mutated_disposition_is_rejected_before_resweep(self):
        bad = copy.deepcopy(_document())
        bad["disposition"] = "R55_SOLVED"
        with self.assertRaises(check.CheckViolation):
            check.verify_document(bad, resweep=False)


if __name__ == "__main__":
    unittest.main()
