#!/usr/bin/env python3
"""Exact mixed (Engström K=4) identity kernel tests — TDD Task 1.

Every observable contract below is specified in
notes/engstrom_identity_mixed_plan_2026-08-18.md, Task 1. The fast
counters are cross-validated against a naive 4-subset enumeration on all
labeled graphs through six vertices (33,867 graphs), the corrected-math
mutation pins are asserted on their named witness graphs, the three pc3
closed forms are enforced equal exhaustively, and the identity residual is
replayed on all 656 published R(5,5,42) graphs and complements plus the
n=49 type-combo fixture.
"""

import itertools
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line  # noqa: E402
from subgraph_identities import (  # noqa: E402
    complement,
    edge_count,
    induced_p3_count,
    induced_subgraph,
    triangle_count,
)
from m4_subgraph_identities import (  # noqa: E402
    diamond_count,
    independent_quad_count,
    induced_k3k1_count,
    k4_count,
)
from mixed_subgraph_identities import (  # noqa: E402  (RED before Task 1)
    MotifSweep,
    c4_count,
    claw_count,
    complement_p3_count,
    k2_2k1_count,
    mixed_residual,
    mixed_vertex_row,
    motif_sweep,
    p3_k1_count,
    p4_count,
    paw_count,
    two_k2_count,
)


def graph_from_mask(n, mask):
    """Convert a vertex-pair bitmask to an adjacency list of bitmasks."""
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


# ------------------------------------------------------------ naive refs ---

_IDX = {}
for _n in range(1, 9):
    _IDX[_n] = list(itertools.combinations(range(_n), 4))


def _is_edge(adj, u, v):
    return (adj[u] >> v) & 1


_PB = {(i, j): bit for bit, (i, j) in enumerate(itertools.combinations(range(4), 2))}


def _canon(n_adj_tuple):
    """Canonical form of an adjacency tuple on 4 labeled vertices."""
    best = None
    for perm in itertools.permutations(range(4)):
        pat = 0
        for i in range(4):
            for j in range(i + 1, 4):
                if (n_adj_tuple[perm[i]] >> perm[j]) & 1:
                    pat |= 1 << _PB[(i, j)]
        if best is None or pat < best:
            best = pat
    return best


# Memoized canonical pattern for every 6-bit 4-vertex graph pattern.
_CANONPAT = {}
for _p in range(1 << 6):
    _adjp = [0] * 4
    for (u, v), bit in _PB.items():
        if (_p >> bit) & 1:
            _adjp[u] |= 1 << v
            _adjp[v] |= 1 << u
    _CANONPAT[_p] = _canon(tuple(_adjp))


# The eleven induced 4-vertex class patterns, built from explicit edge lists
# on {0,1,2,3} and canonicalized (relabeling-invariant 6-bit patterns).
def _class_pattern(edges):
    a = [0] * 4
    for (i, j) in edges:
        a[i] |= 1 << j
        a[j] |= 1 << i
    return _canon(tuple(a))


_K4 = _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
_DIA = _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)])   # K4 - e
_PAW = _class_pattern([(0, 1), (0, 2), (0, 3), (1, 2)])         # paw
_C4 = _class_pattern([(0, 1), (0, 2), (1, 3), (2, 3)])          # C4
_CLAW = _class_pattern([(0, 1), (0, 2), (0, 3)])                  # K1,3
_P4 = _class_pattern([(0, 1), (1, 2), (2, 3)])                    # P4
_K3K1 = _class_pattern([(0, 1), (0, 2), (1, 2)])                  # K3 + K1
_P3K1 = _class_pattern([(0, 1), (1, 2)])                            # P3 + K1
_K22 = _class_pattern([(0, 1), (2, 3)])                             # 2K2
_K2_2K1 = _class_pattern([(0, 1)])                                   # K2 + 2K1
_E4 = _class_pattern([])                                                  # empty


def _canon_pat_atom(adj, quad):
    """Canonical class of the 4-subset graph via the memoized table."""
    pat = 0
    for (i, j), bit in _PB.items():
        if _is_edge(adj, quad[i], quad[j]):
            pat |= 1 << bit
    return _CANONPAT[pat]


def _class_count(adj, cls_pattern):
    return sum(1 for q in _IDX[len(adj)] if _canon_pat_atom(adj, q) == cls_pattern)


def naive_p4(adj):
    return _class_count(adj, _P4)


def naive_c4(adj):
    return _class_count(adj, _C4)


def naive_claw(adj):
    return _class_count(adj, _CLAW)


def naive_paw(adj):
    return _class_count(adj, _PAW)


def naive_two_k2(adj):
    return _class_count(adj, _K22)


def naive_k2_2k1(adj):
    return _class_count(adj, _K2_2K1)


def naive_p3_k1(adj):
    return _class_count(adj, _P3K1)


def _nc(n):
    return itertools.combinations(range(n), 2)


def _pairs_inside_mask(adj, mask):
    e = 0
    mm = mask
    while mm:
        bit = mm & -mm
        v = bit.bit_length() - 1
        e += (mm & adj[v]).bit_count()
        mm ^= bit
    return e


def pc3_form_a(adj):
    """s(P3, complement): sum over edges of (nu - d_u - d_v + c_uv)."""
    n = len(adj)
    degs = [r.bit_count() for r in adj]
    total = 0
    for u, v in _nc(n):
        if not _is_edge(adj, u, v):
            continue
        c = (adj[u] & adj[v]).bit_count()
        total += n - degs[u] - degs[v] + c
    return total


def pc3_form_a_without_cuv(adj):
    n = len(adj)
    degs = [r.bit_count() for r in adj]
    return sum(n - degs[u] - degs[v] for u, v in _nc(n) if _is_edge(adj, u, v))


def pc3_form_c(adj):
    """sum_v C(nu-1-d_v, 2) - 3*i3, with i3 = independent triples."""
    n = len(adj)
    degs = [r.bit_count() for r in adj]
    i3 = triangle_count(complement(adj))
    return sum((n - 1 - d) * (n - 2 - d) // 2 for d in degs) - 3 * i3


def _c3(d):
    return d * (d - 1) * (d - 2) // 6


def mutated_residual(adj, mutation):
    """Residual with one coefficient of the mixed row perturbed."""
    n = len(adj)
    total = 0
    for v in range(n):
        d = (adj[v]).bit_count()
        x = motif_sweep(
            induced_subgraph(
                adj, [u for u in range(n) if (adj[v] >> u) & 1]
            )
        )
        z = motif_sweep(
            complement(
                induced_subgraph(
                    adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
                )
            )
        )
        row = mixed_vertex_row(n, d, x, z)
        if mutation == "c4_coef":
            row += x.c4  # 72 -> 73 on s(C4, X)
        elif mutation == "cross_coef":
            # p_3: -2(n-2) -> -2(n-1): delta -2*d*s(K2,Y), exact (Y=comp Z)
            ey = (n - 1 - d) * (n - 2 - d) // 2 - z.edges
            row -= 2 * d * ey
        elif mutation == "p2_c4_coef":
            row += z.two_k2  # -8 -> -7 on s(C4, Y) = two_k2(Z)
        total += row
    return total


class MixedIdentityTests(unittest.TestCase):
    def test_fast_counters_match_naive_through_six_vertices(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(c4_count(adj), naive_c4(adj), (n, mask))
                self.assertEqual(claw_count(adj), naive_claw(adj), (n, mask))
                self.assertEqual(paw_count(adj), naive_paw(adj), (n, mask))
                self.assertEqual(two_k2_count(adj), naive_two_k2(adj), (n, mask))
                self.assertEqual(k2_2k1_count(adj), naive_k2_2k1(adj), (n, mask))
                self.assertEqual(p3_k1_count(adj), naive_p3_k1(adj), (n, mask))
                self.assertEqual(p4_count(adj), naive_p4(adj), (n, mask))
                self.assertEqual(
                    complement_p3_count(adj), naive_p3(complement(adj)), (n, mask)
                )

    def test_complement_pair_map_exhaustive_through_six_vertices(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                cj = complement(adj)
                self.assertEqual(naive_p3(cj), complement_p3_count(adj), (n, mask))
                self.assertEqual(naive_c4(cj), two_k2_count(adj), (n, mask))
                self.assertEqual(naive_claw(cj), induced_k3k1_count(adj), (n, mask))
                self.assertEqual(naive_paw(cj), p3_k1_count(adj), (n, mask))
                self.assertEqual(naive_diamond(cj), k2_2k1_count(adj), (n, mask))
                self.assertEqual(naive_k4(cj), independent_quad_count(adj), (n, mask))
                self.assertEqual(naive_p4(cj), p4_count(adj), (n, mask))

    def test_pc3_three_closed_forms_agree_exhaustive(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                a = pc3_form_a(adj)
                b = (
                    edge_count(adj) * (len(adj) - 2)
                    - 2 * induced_p3_count(adj)
                    - 3 * triangle_count(adj)
                )
                c = pc3_form_c(adj)
                self.assertEqual(a, b, (n, mask))
                self.assertEqual(b, c, (n, mask))
        k3 = graph_from_mask(3, 0b111)
        self.assertEqual(pc3_form_a(k3), 0)
        self.assertEqual(pc3_form_a_without_cuv(k3), -3)

    def test_motif_sweep_matches_direct_counters_through_six_vertices(self):
        # Guards the single-pass reuse chain (p3_k1 via complement-paw with
        # s(diamond, comp) = k2_2k1 and s(K4, comp) = i4) against the direct
        # counters; pc3 via the motif3-only form B.
        for n in range(4, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                sw = motif_sweep(adj)
                self.assertEqual(sw.order, n)
                self.assertEqual(sw.e4, naive_i4(adj))
                self.assertEqual(sw.k2_2k1, k2_2k1_count(adj))
                self.assertEqual(sw.two_k2, two_k2_count(adj))
                self.assertEqual(sw.p3_k1, p3_k1_count(adj))
                self.assertEqual(sw.p4, p4_count(adj))
                self.assertEqual(sw.claw, claw_count(adj))
                self.assertEqual(sw.k3k1, induced_k3k1_count(adj))
                self.assertEqual(sw.c4, c4_count(adj))
                self.assertEqual(sw.paw, paw_count(adj))
                self.assertEqual(sw.diamond, diamond_count(adj))
                self.assertEqual(sw.k4, k4_count(adj))
                self.assertEqual(sw.pc3, complement_p3_count(adj))

    def test_mixed_residual_exhaustive_through_six_vertices(self):
        total = 0
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(mixed_residual(graph_from_mask(n, mask)), 0, (n, mask))
                total += 1
        self.assertEqual(total, 33867)

    def test_spot_random_seven_eight_vertices(self):
        rng = random.Random(20260818)
        for _ in range(40):
            n = rng.choice([7, 8])
            adj = [0] * n
            for u, v in itertools.combinations(range(n), 2):
                if rng.random() < rng.choice([0.3, 0.5, 0.7]):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
            self.assertEqual(mixed_residual(adj), 0)

    def test_mutated_coefficients_are_detected(self):
        for mutation in ("c4_coef", "cross_coef", "p2_c4_coef"):
            worst = max(
                abs(mutated_residual(graph_from_mask(n, mask), mutation))
                for n in range(4, 7)
                for mask in range(1 << (n * (n - 1) // 2))
            )
            self.assertGreater(worst, 0, mutation)

    def test_closed_form_corrections_are_detected(self):
        diamond = graph_from_mask(4, 0b111110)  # K4 minus edge (0,1)
        self.assertEqual(sum(_c3(d) for d in (3, 3, 2, 2)), 2)
        self.assertEqual(claw_count(diamond), 0)  # = 2 - paw - 2*dia - 4*K4
        two_k2 = graph_from_mask(4, 0b010010)  # edges: (0,2), (1,3)
        self.assertEqual(k2_2k1_count(two_k2), 0)  # bare C(2,2)*2 = 2 - 2*1
        k4 = graph_from_mask(4, (1 << 6) - 1)
        self.assertEqual(paw_count(k4), 0)  # = 12 - 0*dia - 12*K4

    def test_c4_incidence_relation_exhaustive(self):
        # For every non-edge uv with cn = common-neighbor set of u and v,
        # C(|cn|,2) - e(G[cn]) counts the induced C4s whose two non-edges
        # include uv: a non-adjacent pair inside cn completes a C4, an
        # adjacent pair closes a diamond (which has only its ONE non-edge, hence
        # contributes nothing here).  Each C4 is counted exactly twice:
        #   sum_{uv not in E} [C(c_uv,2) - e(G[cn(u,v)])] = 2*C4.
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                s = 0
                for u, v in _nc(n):
                    if _is_edge(adj, u, v):
                        continue
                    cn = adj[u] & adj[v]
                    c = cn.bit_count()
                    s += c * (c - 1) // 2 - _pairs_inside_mask(adj, cn)
                self.assertEqual(s, 2 * c4_count(adj), (n, mask))

    def test_all_656_known_r55_graphs_have_zero_mixed_residual(self):
        published = []
        with (ROOT / "data" / "r55_42some.g6").open(encoding="ascii") as fh:
            for line in fh:
                if line.strip():
                    _, adj = parse_graph6_line(line)
                    published.append(adj)
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(mixed_residual(adj), 0)
            self.assertEqual(mixed_residual(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

    def test_n49_type_combo_replay(self):
        graphs = []
        with (ROOT / "data" / "r45_24.g6").open(encoding="ascii") as fh:
            for line in fh:
                if not line.strip():
                    continue
                _, adj = parse_graph6_line(line)
                if sum(row.bit_count() for row in adj) // 2 == 132:
                    graphs.append(adj)
        self.assertEqual(len(graphs), 2)
        for x in graphs:  # 11-regular; Thm 3.2 anchor
            sw = motif_sweep(x)
            self.assertEqual(sw.triangles, 176)
            self.assertEqual(sw.diamond, 792)
            self.assertEqual(sw.paw, 1584)  # new deliberate pin on X1/X2
            self.assertEqual(sw.claw, 792)  # = sum C(11,3)*24 - paw - 2*diamond
            self.assertEqual(sw.k3k1, 528)
        values = sorted(
            mixed_vertex_row(49, 24, motif_sweep(xi), motif_sweep(zj))
            for xi in graphs
            for zj in graphs
        )
        self.assertEqual(values, [0, 144, 288, 432])
        self.assertEqual(sum(v == 0 for v in values), 1)  # unique zero corner


# Naive p3 count of an induced graph (already used above).
def naive_p3(adj):
    return induced_p3_count(adj)


def naive_i4(adj):
    return independent_quad_count(adj)


def naive_diamond(adj):
    return diamond_count(adj)


def naive_k4(adj):
    return k4_count(adj)


if __name__ == "__main__":
    unittest.main()
