#!/usr/bin/env python3
"""Exact local m=3 subgraph identity for the R(5,5) campaign.

For a graph G on n vertices write G_v^+ for the subgraph induced on the
neighborhood of v and G_v^- for the subgraph induced on the non-neighborhood
of v (v itself excluded from both). All counts below are counts of *induced*
subgraphs: e is edges, t is triangles, p_3 is induced paths on three vertices.

The identity is stored after multiplication by three so that every quantity
is an exact integer:

    3 sum_v t(G_v^-)
        = sum_v [ (n + 3 - 3 d_v) e(G_v^+) + 3 p_3(G_v^+) + 6 t(G_v^+) ].

Independent counting derivation (unscaled)
------------------------------------------
Divide the displayed identity by three:

    sum_v t(G_v^-)
        = sum_v [ ((n + 3)/3 - d_v) e(G_v^+) + p_3(G_v^+) + 2 t(G_v^+) ].

Left side. A triangle in the non-neighborhood of v is a triple {a,b,c}
spanning a triangle with v adjacent to none of a, b, c, i.e. an induced
K_3 + K_1 (disjoint union) with v as the isolated vertex. Each induced
K_3 + K_1 has a unique isolated vertex, so the left side counts induced
K_3 + K_1 subgraphs of G exactly once each.

First right-hand term. An edge of G_v^+ is a triangle of G through v, so
summing over v regroups this term by triangles: each triangle T = abc is
visited once per vertex v in T, contributing

    sum_{v in T} ((n + 3)/3 - d_v) = (n + 3) - sum_{v in T} d_v.

Let r_x be the number of neighbors of x inside T. Since the edges leaving T
number sum_{v in T} (d_v - 2), the same value is reached from outside:

    sum_{x not in T} (1 - r_x) = (n - 3) - (sum_{v in T} d_v - 6)
                              = (n + 3) - sum_{v in T} d_v,

so the first term equals sum_{triangle T} sum_{x not in T} (1 - r_x). Every
summand is now read off the four-vertex subgraph induced on T + {x} by r_x:

    r_x = 0  ->  K_3 + K_1  ->  +1  (its single triangle, so +1 per K_3 + K_1)
    r_x = 1  ->  paw        ->   0  (its single triangle, so 0 per paw)
    r_x = 2  ->  diamond    ->  -1  (2 triangles per diamond, so -2 each)
    r_x = 3  ->  K_4        ->  -2  (4 triangles per K_4, so -8 each)

Correction terms. p_3(G_v^+) is nonzero exactly when the neighborhood of v
contains an induced path, which on four vertices means {v} plus an induced
path: summed over v this contributes +2 per diamond (the two degree-3
vertices of a diamond each see an induced P_3) and 0 per K_4 and paw.
2 t(G_v^+) contributes 2 per (v, triangle in N(v)) pair, i.e. +8 per K_4
(all four vertices see a triangle) and 0 per diamond, paw and K_3 + K_1.

Adding the three groups, diamonds cancel (-2 + 2 + 0 = 0), K_4 cancels
(-8 + 0 + 8 = 0), paws cancel (0 + 0 + 0 = 0), and only the +1 per induced
K_3 + K_1 survives: exactly the unscaled left side. Multiplying back by
three yields the stored integer identity.

`tests/test_subgraph_identities.py` re-verifies the identity by exhaustive
enumeration of all labeled graphs on up to six vertices, a mechanically
independent check of this derivation.
"""

from collections import namedtuple

from check_ramsey import popcount

Motif3 = namedtuple(
    "Motif3", "order edges triangles induced_p3 independent_triples"
)


def complement(adj):
    """Complement graph, same vertex order, no loops."""
    n = len(adj)
    full = (1 << n) - 1
    return [full & ~adj[v] & ~(1 << v) for v in range(n)]


def induced_subgraph(adj, vertices):
    """Subgraph induced on `vertices`, relabeled 0..len(vertices)-1 in order."""
    out = [0] * len(vertices)
    for i, u in enumerate(vertices):
        for j in range(i + 1, len(vertices)):
            v = vertices[j]
            if (adj[u] >> v) & 1:
                out[i] |= 1 << j
                out[j] |= 1 << i
    return out


def edge_count(adj):
    return sum(popcount(row) for row in adj) // 2


def triangle_count(adj):
    total = 0
    for u in range(len(adj)):
        later = adj[u] & ~((1 << (u + 1)) - 1)
        while later:
            bit = later & -later
            v = bit.bit_length() - 1
            total += popcount(adj[u] & adj[v])
            later ^= bit
    return total // 3


def induced_p3_count(adj):
    """Induced paths on three vertices: wedges minus the three per triangle."""
    triangles = triangle_count(adj)
    wedges = sum(popcount(row) * (popcount(row) - 1) // 2 for row in adj)
    return wedges - 3 * triangles


def independent_triple_count(adj):
    return triangle_count(complement(adj))


def motif3(adj):
    """All three-vertex induced counts of `adj`, plus its order."""
    return Motif3(
        len(adj), edge_count(adj), triangle_count(adj), induced_p3_count(adj),
        independent_triple_count(adj),
    )


def m3_vertex_scaled(n, neighborhood, dual):
    """Scaled residual of the identity at one vertex of an n-vertex graph.

    `neighborhood` is G_v^+ and `dual` is G_v^-, both as induced subgraphs.
    """
    if len(neighborhood) + len(dual) != n - 1:
        raise ValueError("neighborhood and dual orders must sum to n-1")
    d = len(neighborhood)
    f = motif3(neighborhood)
    return 3 * triangle_count(dual) - (
        (n + 3 - 3 * d) * f.edges + 3 * f.induced_p3 + 6 * f.triangles
    )


def m3_residual_scaled(adj):
    """Total scaled residual; zero for every graph by the identity above."""
    n = len(adj)
    return sum(
        m3_vertex_scaled(
            n,
            induced_subgraph(adj, [u for u in range(n) if (adj[v] >> u) & 1]),
            induced_subgraph(
                adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
            ),
        )
        for v in range(n)
    )
