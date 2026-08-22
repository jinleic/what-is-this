"""Tests for the level-2 pair lift and the identity it retires.

The load-bearing test here is `test_invalid_ez_row_is_the_trap`.  A misread of
the z-stratum convention produces a row that looks like a breakthrough and is
not sound; that test pins both halves of the trap so it cannot be silently
reintroduced.
"""

from __future__ import annotations

import itertools
import json
import random
import sys
import unittest
from fractions import Fraction
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pair_lift_mixed as lift                                  # noqa: E402

_ARTIFACT = Path(__file__).resolve().parents[1] / "data" / "engstrom_identity.json"


def _states():
    return lift.load_states(_ARTIFACT)


def _brute_force_facts(adj, n):
    """Independent, slow reference: sets and combinations, no bitmask tricks."""
    nbrs = [{u for u in range(n) if adj[v] >> u & 1} for v in range(n)]
    deg = [len(nbrs[v]) for v in range(n)]
    edges = {frozenset((u, v)) for u in range(n) for v in nbrs[u]}
    e = len(edges)
    inside = lambda members: sum(                               # noqa: E731
        1 for u, v in itertools.combinations(sorted(members), 2) if v in nbrs[u])
    per_vertex = []
    for v in range(n):
        dv = set(range(n)) - nbrs[v] - {v}
        per_vertex.append((sum(deg[u] for u in nbrs[v]), inside(nbrs[v]), inside(dv)))
    return e, deg, per_vertex


class TestIdentity(unittest.TestCase):
    def test_identity_against_independent_brute_force(self):
        """(P) and (S) verified against a set-based reference on all n<=6."""
        checked = 0
        for n in range(2, 7):
            pairs = list(itertools.combinations(range(n), 2))
            for mask in range(1 << len(pairs)):
                adj = [0] * n
                for bit, (i, j) in enumerate(pairs):
                    if mask >> bit & 1:
                        adj[i] |= 1 << j
                        adj[j] |= 1 << i
                e, deg, per_vertex = _brute_force_facts(adj, n)
                for nbr_deg, ex, ed in per_vertex:
                    self.assertEqual(nbr_deg, e + ex - ed)      # identity (P)
                self.assertEqual(                               # identity (S)
                    sum(d * d for d in deg),
                    n * e + sum(p[1] for p in per_vertex)
                    - sum(p[2] for p in per_vertex))
                checked += 1
        self.assertEqual(checked, 2 + 8 + 64 + 1024 + 32768)

    def test_module_sweep_agrees(self):
        stats = lift.verify_identity_exhaustive(max_order=6, random_orders=(11,),
                                                per_order=5)
        self.assertEqual(stats["exhaustive_graphs"], 2 + 8 + 64 + 1024 + 32768)
        self.assertEqual(stats["random_graphs"], 5)

    def test_identity_detects_a_broken_graph(self):
        """The checker must fail on data that violates the identity."""
        adj = [0b110, 0b101, 0b011]
        self.assertTrue(lift._identity_holds(adj, 3))
        adj_broken = [0b110, 0b101, 0b001]                      # asymmetric
        self.assertFalse(lift._identity_holds(adj_broken, 3))


class TestRowCollapse(unittest.TestCase):
    def test_row_equals_excess_balance_on_every_state(self):
        states = _states()
        self.assertEqual(lift.assert_row_is_excess_balance(states), 3215)

    def test_row_adds_no_rank(self):
        states = _states()
        self.assertEqual(lift.exact_rank(states), 8)
        self.assertEqual(lift.exact_rank(states, lift.neighbour_degree_row), 8)

    def test_invalid_ez_row_is_the_trap(self):
        """Pin both halves: the trap row is independent AND unsound.

        Independent, so it looks like a new constraint and appears to cut every
        route hard.  Unsound, because it substitutes the recorded e_z (edges of
        the COMPLEMENT on D(v)) for e(G[D(v)]) = e_y.  Anyone who re-derives it
        must fail this test rather than publish the bound.
        """
        states = _states()
        self.assertEqual(lift.exact_rank(states, lift.invalid_ez_row), 9)
        differs = [s for s in states
                   if lift.invalid_ez_row(s) != s["excess_balance"]]
        self.assertTrue(differs, "trap row must not coincide with the m=2 row")
        # the sound row and the trap row disagree precisely by the e_z/e_y swap
        for state in states[:200]:
            m = 44 - state["d"]
            swap = 2 * (lift.c2(m) - 2 * (lift.E45[m] - state["b"]))
            self.assertEqual(lift.neighbour_degree_row(state),
                             lift.invalid_ez_row(state) + swap)

    def test_mis_keyed_row_is_not_even_an_identity(self):
        """The cheapest possible falsifier, and the one that should run first.

        The sound row sums to zero over the vertices of ANY graph, because it is
        identity (S).  The mis-keyed row does not: on the path $P_4$ it sums to
        $-4$, and it is nonzero on most small graphs.  A single call to this
        would have refuted the bad row before any LP was written.
        """
        def sums(adj, n):
            deg = [bin(a).count("1") for a in adj]
            full = (1 << n) - 1
            sound = mis_keyed = 0
            for v in range(n):
                nv = adj[v]
                dv = full & ~nv & ~(1 << v)
                d = deg[v]
                m = n - 1 - d
                ex = lift._edges_within(adj, nv)
                ed = lift._edges_within(adj, dv)
                sound += 2 * d * d - n * d - 2 * ex + 2 * ed
                mis_keyed += 2 * d * d - n * d - 2 * ex + 2 * (lift.c2(m) - ed)
            return sound, mis_keyed

        path4 = [0b0010, 0b0101, 0b1010, 0b0100]
        self.assertEqual(sums(path4, 4), (0, -4))

        nonzero = 0
        total = 0
        for n in range(3, 7):
            pairs = list(itertools.combinations(range(n), 2))
            for mask in range(1 << len(pairs)):
                adj = [0] * n
                for bit, (i, j) in enumerate(pairs):
                    if mask >> bit & 1:
                        adj[i] |= 1 << j
                        adj[j] |= 1 << i
                sound, mis_keyed = sums(adj, n)
                self.assertEqual(sound, 0, "identity (S) must hold on every graph")
                total += 1
                if mis_keyed != 0:
                    nonzero += 1
        self.assertEqual(total, 8 + 64 + 1024 + 32768)
        self.assertEqual(nonzero, 4 + 54 + 938 + 30652)

    def test_identity_is_owned_by_theorem_2(self):
        """SSOT: identity (S) is already Theorem 2 of the structural notes."""
        note = (Path(__file__).resolve().parents[1] / "notes"
                / "structural_constraints_2026-08-15.md")
        text = note.read_text()
        self.assertIn("Theorem 2 (excess identity)", text)
        self.assertIn("e(F_v^-) − e(F_v^+)", text)

    def test_ex_and_ey_match_recorded_convention(self):
        states = _states()
        for state in states[:400]:
            m = 44 - state["d"]
            self.assertEqual(lift.e_x(state), lift.E45[state["d"]] - state["a"])
            self.assertEqual(lift.e_y(state),
                             lift.c2(m) - (lift.E45[m] - state["b"]))
            self.assertEqual(lift.phi(state),
                             lift.e_x(state) - lift.e_y(state))


class TestCodegreeFacts(unittest.TestCase):
    def test_adjacent_bound_is_the_r35_bound(self):
        self.assertEqual(lift.lam_adjacent_max(22, 22), 13)
        self.assertEqual(lift.lam_adjacent_max(20, 20), 13)

    def test_nonadjacent_bound_is_an_upper_bound(self):
        """Sign check.  Reading this as a lower bound fakes a refutation."""
        self.assertEqual(lift.lam_nonadjacent_max(22, 22), 14)
        self.assertEqual(lift.lam_nonadjacent_max(20, 20), 10)
        self.assertEqual(lift.lam_nonadjacent_max(24, 24), 18)
        for d in lift.STATE_DEGREES:
            for dp in lift.STATE_DEGREES:
                self.assertLessEqual(lift.lam_nonadjacent_min(d, dp),
                                     lift.lam_nonadjacent_max(d, dp))

    def test_codegree_counting_on_a_concrete_graph(self):
        """sum over ordered non-adjacent pairs of lam = 2 (sum C(d,2) - 3T)."""
        n = 45
        adj = [0] * n
        for v in range(n):
            for off in range(1, 12):
                u = (v + off) % n
                adj[v] |= 1 << u
                adj[u] |= 1 << v
        deg = [bin(a).count("1") for a in adj]
        triangles = sum(
            bin(adj[u] & adj[v]).count("1")
            for u in range(n) for v in range(u + 1, n) if adj[u] >> v & 1) // 3
        lam = lambda u, v: bin(adj[u] & adj[v]).count("1")      # noqa: E731
        ordered_adj = sum(lam(u, v) for u in range(n) for v in range(n)
                          if u != v and adj[u] >> v & 1)
        ordered_non = sum(lam(u, v) for u in range(n) for v in range(n)
                          if u != v and not adj[u] >> v & 1)
        self.assertEqual(ordered_adj, 6 * triangles)
        self.assertEqual(ordered_non,
                         2 * (sum(lift.c2(d) for d in deg) - 3 * triangles))


class TestLift(unittest.TestCase):
    def test_exact_bounds_and_verdict(self):
        states = _states()
        bounds = lift.lift_bounds(states)
        self.assertEqual(bounds.pop("_e_window"), lift.EXPECTED_E_WINDOW)
        self.assertEqual(bounds, dict(lift.EXPECTED_LIFT))

    def test_no_route_is_accepted(self):
        for route, value in lift.EXPECTED_LIFT.items():
            edge, op = lift.ACCEPTANCE[route]
            accepted = value <= edge if op == "<=" else value < edge
            self.assertFalse(accepted, f"{route} must stay rejected")

    def test_lift_never_exceeds_the_frozen_bound(self):
        """The lift restricts the frozen programme, so it cannot get worse."""
        for route, value in lift.EXPECTED_LIFT.items():
            self.assertLessEqual(value, lift.FROZEN[route])

    def test_degree20_improves_by_exactly_one_over_349(self):
        self.assertEqual(lift.FROZEN["degree20_count"]
                         - lift.EXPECTED_LIFT["degree20_count"],
                         Fraction(1, 349))

    def test_dual_certificate_is_verified_not_trusted(self):
        """A tampered objective must not keep the recorded bound."""
        states = _states()
        value = lift.exact_dual_bound(states, "degree20_count", 500)
        self.assertIsNotNone(value)
        self.assertLessEqual(value, Fraction(41))

    def test_programme_is_infeasible_outside_the_e_window(self):
        states = _states()
        low, high = lift.EXPECTED_E_WINDOW
        self.assertIsNone(lift.exact_dual_bound(states, "total_deficiency", low - 1))
        self.assertIsNone(lift.exact_dual_bound(states, "total_deficiency", high + 1))


class TestNoBoundClaimed(unittest.TestCase):
    def test_module_claims_no_ramsey_bound(self):
        text = (_SRC / "pair_lift_mixed.py").read_text()
        self.assertIn("NO BOUND ON R(5,5) IS", text)
        self.assertIn("MIXED_PAIR_LIFT_CLOSED", text)

    def test_artifact_is_not_mutated(self):
        before = _ARTIFACT.read_bytes()
        lift.load_states(_ARTIFACT)
        lift.exact_rank(_states(), lift.neighbour_degree_row)
        self.assertEqual(_ARTIFACT.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
