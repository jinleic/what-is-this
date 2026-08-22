import itertools
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line
from subgraph_identities import complement, induced_subgraph
from m4_subgraph_identities import (  # expected to fail before Task 1 implementation
    diamond_count, independent_quad_count, induced_k3k1_count, k4_count,
    motif4, m4_residual_specialized, m4_residual_universal,
    m4_vertex_row_scaled,
)


def graph_from_mask(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


# ---------------------------------------------------------------------------
# Test-only naive references: 4-subset enumeration, independent of the
# production common-neighborhood algebra.
# ---------------------------------------------------------------------------


def _quartet_edges(adj, quad):
    return sum(
        1 for u, v in itertools.combinations(quad, 2) if (adj[u] >> v) & 1
    )


def _quartet_triangles(adj, quad):
    return sum(
        1
        for a, b, c in itertools.combinations(quad, 3)
        if (adj[a] >> b) & 1 and (adj[a] >> c) & 1 and (adj[b] >> c) & 1
    )


def naive_diamond(adj):
    """Induced K4-minus-an-edge: the quartet spans exactly five edges."""
    return sum(
        1 for q in itertools.combinations(range(len(adj)), 4)
        if _quartet_edges(adj, q) == 5
    )


def naive_i4(adj):
    """Induced independent quartet: zero edges."""
    return sum(
        1 for q in itertools.combinations(range(len(adj)), 4)
        if _quartet_edges(adj, q) == 0
    )


def naive_k3k1(adj):
    """Induced K3+K1: three edges that form exactly one triangle."""
    return sum(
        1 for q in itertools.combinations(range(len(adj)), 4)
        if _quartet_edges(adj, q) == 3 and _quartet_triangles(adj, q) == 1
    )


def naive_k4(adj):
    """Induced clique quartet: six edges."""
    return sum(
        1 for q in itertools.combinations(range(len(adj)), 4)
        if _quartet_edges(adj, q) == 6
    )


def _naive_triangles(adj):
    n = len(adj)
    return sum(
        1
        for a, b, c in itertools.combinations(range(n), 3)
        if (adj[a] >> b) & 1 and (adj[a] >> c) & 1 and (adj[b] >> c) & 1
    )


def mutated_universal_residual(adj, mutation):
    """Universal residual with one coefficient perturbed (test-only).

    Recomputed purely from naive enumeration, so the mutation scan checks
    the identity constants independently of the production kernel:
    sum_v [3(n+2-2d_v) t(N_v) + c_diamond diamond(N_v) - 6 s(K3+K1, N_v)]
        + c_k4 sum_v K4(N_v) - c_k4 sum_v K4(D_v),
    with K4(D_v) = i4(complement(D_v)); c_diamond = 4 -> 3 under
    "diamond_coef" and c_k4 = 12 -> 11 under "k4_correction".
    """
    if mutation == "diamond_coef":
        diamond_coef = 3
    elif mutation == "k4_correction":
        diamond_coef = 4
    else:
        raise ValueError("unknown mutation")
    k4_correction = 11 if mutation == "k4_correction" else 12
    n = len(adj)
    total = 0
    for v in range(n):
        nbhd = induced_subgraph(
            adj, [u for u in range(n) if (adj[v] >> u) & 1]
        )
        dual = induced_subgraph(
            adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
        )
        total += 3 * (n + 2 - 2 * len(nbhd)) * _naive_triangles(nbhd)
        total += diamond_coef * naive_diamond(nbhd)
        total -= 6 * naive_k3k1(nbhd)
        total += k4_correction * naive_k4(nbhd)
        total -= k4_correction * naive_i4(complement(dual))
    return total


def _one_triangle_neighborhood(d):
    """Order-d neighborhood carrying exactly one triangle and no other
    4-vertex motif: triangle {0, 1, 2} plus d-3 pendants on vertex 0."""
    nbhd = [0] * d
    nbhd[0] = (1 << 1) | (1 << 2) | (((1 << (d - 3)) - 1) << 3)
    nbhd[1] = (1 << 0) | (1 << 2)
    nbhd[2] = (1 << 0) | (1 << 1)
    for w in range(3, d):
        nbhd[w] = 1 << 0
    return nbhd


class M4IdentityTests(unittest.TestCase):
    def test_fast_counters_match_naive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(diamond_count(adj), naive_diamond(adj))
                self.assertEqual(induced_k3k1_count(adj), naive_k3k1(adj))
                self.assertEqual(independent_quad_count(adj), naive_i4(adj))
                self.assertEqual(k4_count(adj), naive_k4(adj))

    def test_m4_universal_identity_exhaustive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(m4_residual_universal(graph_from_mask(n, mask)), 0)

    def test_specialized_identity_on_k4free_neighborhood_graphs(self):
        tested = 0  # plan-drafting reference count over all labeled n<=6: 33,694
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                if any(k4_count(induced_subgraph(
                        adj, [u for u in range(n) if (adj[v] >> u) & 1]))
                       for v in range(n)):
                    continue
                self.assertEqual(m4_residual_specialized(adj), 0)
                tested += 1
        self.assertGreater(tested, 30000)

    def test_specialized_raises_on_k4_neighborhood(self):
        k5 = graph_from_mask(5, (1 << 10) - 1)
        with self.assertRaises(ValueError):
            m4_residual_specialized(k5)

    def test_mutated_coefficients_are_detected(self):
        # Each perturbation is caught separately: diamond coefficient 4 -> 3
        # and K4-correction 12 -> 11 must each yield a nonzero residual on at
        # least one labeled graph on <= 6 vertices.
        for mutation in ("diamond_coef", "k4_correction"):
            worst = max(
                abs(mutated_universal_residual(
                    graph_from_mask(n, mask), mutation))
                for n in range(4, 7)
                for mask in range(1 << (n * (n - 1) // 2))
            )
            self.assertGreater(worst, 0, mutation)
        # Hand-derived witnesses (each sums to exactly zero unperturbed):
        # K5: every neighborhood is a K4, total 5*(3*(7-8)*4 + c) = 5c - 60.
        self.assertEqual(
            mutated_universal_residual(
                graph_from_mask(5, (1 << 10) - 1), "k4_correction"),
            -5,
        )
        # K5 minus one edge: vertices 2, 3, 4 each see exactly one diamond,
        # total 6 + 3*(coef - 6) = 3*coef - 12.
        self.assertEqual(
            mutated_universal_residual(
                graph_from_mask(5, (1 << 10) - 2), "diamond_coef"),
            -3,
        )

    def test_mr97_n49_ground_truth_replay(self):
        graphs = []
        with (ROOT / "data" / "r45_24.g6").open(encoding="ascii") as fh:
            for line in fh:
                if not line.strip():
                    continue
                _, adj = parse_graph6_line(line)
                if sum(row.bit_count() for row in adj) // 2 == 132:
                    graphs.append(adj)
        self.assertEqual(len(graphs), 2)
        for x in graphs:
            f = motif4(x)
            self.assertEqual(f.triangles, 176)
            self.assertEqual(f.diamonds, 792)
            self.assertEqual(f.k3k1, 528)
        self.assertEqual(sorted(independent_quad_count(x) for x in graphs),
                         [138, 144])
        for x in graphs:
            row = 3 * (49 + 2 - 2 * 24) * 176 + 4 * 792 - 6 * 528
            self.assertEqual(row, 1584)          # = 12 * 132, forced i4 mean
            self.assertGreater(independent_quad_count(x), 132)

    def test_all_656_known_r55_graphs_have_zero_specialized_residual(self):
        published = []
        with (ROOT / "data" / "r55_42some.g6").open(encoding="ascii") as fh:
            for line in fh:
                if line.strip():
                    _, adj = parse_graph6_line(line)
                    published.append(adj)
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(m4_residual_specialized(adj), 0)
            self.assertEqual(m4_residual_specialized(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

    def test_spot_random_seven_eight_vertices(self):
        rng = random.Random(20260817)
        for _ in range(40):
            n = rng.choice([7, 8])
            adj = [0] * n
            for u, v in itertools.combinations(range(n), 2):
                if rng.random() < rng.choice([0.3, 0.5, 0.7]):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
            self.assertEqual(m4_residual_universal(adj), 0)

    def test_hand_built_motif_values(self):
        # K4 minus one edge: exactly one induced diamond, two triangles.
        diamond4 = graph_from_mask(4, 62)
        self.assertEqual(diamond_count(diamond4), 1)
        self.assertEqual(motif4(diamond4), (4, 2, 1, 0, 0, 0))
        # Empty 4-graph: exactly one independent quad.
        empty4 = [0] * 4
        self.assertEqual(independent_quad_count(empty4), 1)
        self.assertEqual(motif4(empty4), (4, 0, 0, 0, 1, 0))
        # K3 plus an isolated vertex: exactly one induced K3+K1.
        k3k1 = graph_from_mask(4, 11)
        self.assertEqual(induced_k3k1_count(k3k1), 1)
        self.assertEqual(motif4(k3k1), (4, 1, 0, 1, 0, 0))
        # Paw (triangle with pendant): a 4-edge quartet that is neither a
        # diamond nor a K3+K1.
        paw = graph_from_mask(4, 15)
        self.assertEqual(motif4(paw), (4, 1, 0, 0, 0, 0))
        # K4 itself: four triangles, one clique quartet.
        self.assertEqual(motif4(graph_from_mask(4, 63)), (4, 4, 0, 0, 0, 1))
        # C5 built explicitly so bit ordering cannot hide a test error:
        # triangle-free, independence number 2, so every 4-motive is absent.
        c5 = [0] * 5
        for u, v in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0)):
            c5[u] |= 1 << v
            c5[v] |= 1 << u
        self.assertEqual(motif4(c5), (5, 0, 0, 0, 0, 0))

    def test_m4_vertex_row_scaled_sign_table(self):
        # The n=45 usable rows: 3*(47-2d) in {21, 15, 9, 3, -3} for
        # d in {20..24} -- note the sign flip at d=24.  The neighborhood
        # carries exactly one triangle and no other 4-vertex motif; the
        # empty dual makes the complement-side independent-quad term zero.
        for d, expected in zip(range(20, 25), (21, 15, 9, 3, -3)):
            nbhd = _one_triangle_neighborhood(d)
            dual = [0] * (45 - 1 - d)
            self.assertEqual(m4_vertex_row_scaled(45, nbhd, dual), expected)
        with self.assertRaises(ValueError):
            m4_vertex_row_scaled(45, [0] * 20, [0] * 25)  # 45 != 44

    def test_complement_counter_consistency(self):
        # Independent quads of G are clique quartets of the complement and
        # vice versa, through the two distinct fast code paths.
        for n in range(1, 6):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(
                    independent_quad_count(adj), k4_count(complement(adj))
                )
                self.assertEqual(
                    k4_count(adj), independent_quad_count(complement(adj))
                )


if __name__ == "__main__":
    unittest.main()
