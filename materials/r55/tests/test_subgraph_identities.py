import itertools
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from subgraph_identities import (  # expected to fail before Task 1 implementation
    complement, independent_triple_count, induced_p3_count, induced_subgraph,
    m3_residual_scaled, motif3, triangle_count,
)


def graph_from_mask(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


class M3IdentityTests(unittest.TestCase):
    def test_motif_counts(self):
        # Build C5 explicitly so graph-bit ordering cannot hide a test error.
        c5 = [0] * 5
        for u, v in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0)):
            c5[u] |= 1 << v
            c5[v] |= 1 << u
        self.assertEqual(motif3(c5), (5, 5, 0, 5, 0))
        empty4 = [0] * 4
        self.assertEqual(independent_triple_count(empty4), 4)

    def test_complement_is_involution(self):
        for n in range(1, 6):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(complement(complement(adj)), adj)

    def test_m3_identity_exhaustive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(m3_residual_scaled(graph_from_mask(n, mask)), 0)

    def test_mutated_triangle_coefficient_is_detected(self):
        k4 = graph_from_mask(4, (1 << 6) - 1)
        self.assertEqual(m3_residual_scaled(k4), 0)
        # The production identity uses coefficient 6. K4 detects replacing it by 5.
        n = 4
        mutated = 0
        for v in range(n):
            neighborhood_vertices = [
                u for u in range(n) if (k4[v] >> u) & 1
            ]
            dual_vertices = [
                u for u in range(n)
                if u != v and not ((k4[v] >> u) & 1)
            ]
            neighborhood = induced_subgraph(k4, neighborhood_vertices)
            dual = induced_subgraph(k4, dual_vertices)
            f = motif3(neighborhood)
            mutated += 3 * triangle_count(dual) - (
                (n + 3 - 3 * len(neighborhood)) * f.edges
                + 3 * f.induced_p3 + 5 * f.triangles
            )
        self.assertNotEqual(mutated, 0)


if __name__ == "__main__":
    unittest.main()
