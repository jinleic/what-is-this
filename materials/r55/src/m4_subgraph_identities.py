#!/usr/bin/env python3
"""Exact local m=4 subgraph identity for the R(5,5) campaign.

For a graph G on n vertices write N_v for the subgraph induced on the
neighbors of v and D_v for the subgraph induced on the non-neighbors of
v (v itself excluded from both), d_v for the degree of v. All counts
below are counts of *induced* subgraphs: t is triangles, "diamond" is
K4 minus an edge, K3+K1 the disjoint union of a triangle and an isolated
vertex, i4 the number of independent vertex quartets, K4(J) the number
of induced 4-cliques of J, and cn(S) the set of common neighbors of a
vertex set S.

Two identities are stored, both after multiplication by twelve so that
every quantity is an exact integer.

Universal (every graph):

    sum_v [ 3 (n + 2 - 2 d_v) t(N_v) + 4 diamond(N_v)
            - 6 s(K3+K1, N_v) ]
        + 12 sum_v K4(N_v) - 12 sum_v K4(D_v)  =  0.

Specialized to graphs whose every neighborhood is K4-free -- in
particular every (5,5)-graph, since a K4 in N_v would complete a K5
with v -- the K4(N_v) term vanishes and

    sum_v h_v = 0,   h_v = 3 (n + 2 - 2 d_v) t(N_v) + 4 diamond(N_v)
                        - 6 s(K3+K1, N_v) - 12 i4(complement(D_v)),

where K4(D_v) = i4(complement(D_v)): independent quartets of the
complemented dual are exactly the clique quartets of the dual. At n = 45
the usable coefficient 3 (47 - 2 d) takes the values 21, 15, 9, 3, -3
for d = 20..24 -- note the sign flip at d = 24.

Independent counting derivation (the ledger)
--------------------------------------------
L1 (dual side).  A K4 inside D_v is a 4-clique none of whose vertices
is adjacent to v: an induced K4+K1 with v as its isolated vertex. Every
induced K4+K1 has exactly one isolated vertex (its other four vertices
span all six edges), so

    12 sum_v K4(D_v) = 12 * #(induced K4+K1),

each pattern counted once, by its isolated vertex.

L2 (triangle mass).  A triangle of N_v is a triangle T all of whose
vertices v is adjacent to, i.e. T + {v} is a K4 through v; conversely
each K4 S contributes the triangle S - {v} to t(N_v) from each of its
four vertices. Hence sum_v t(N_v) = 4 k4(G). (A paw contributes
nothing: its pendant is adjacent to only one vertex of the triangle, so
the pendant's neighborhood never contains the whole triangle.)

L3 (incidences).  For J in {diamond, K3+K1, K4}, a copy of J inside N_v
is a copy of J in G, not containing v, all of whose vertices v is
adjacent to. So

    sum_v J(N_v) = sum_{S ~= J} |cn(S)|,

an incidence count between induced 4-sets and outside common neighbors.

L4 (regrouping).  Split the degree coefficient

    3 (n + 2 - 2 d_v) = 3 (n - 4) + 6 (3 - d_v)

and regroup both t(N_v) terms by the K4 they extend to (L2), using
sum_{v in S} d_v = 12 + sum_{x not in S} r_x for a K4 S, where r_x
counts the edges from an outside vertex x into S:

    sum_v 3 (n - 4) t(N_v) = sum_{S ~= K4} sum_{x not in S} 12,
    sum_v 6 (3 - d_v) t(N_v)
        = 72 k4(G) - 6 sum_{S ~= K4} sum_{v in S} d_v
        = sum_{S ~= K4} sum_{x not in S} (-6 r_x).

Together with L1 and L3, every term of the universal identity is a
per-incidence weight carried by a 5-vertex set F = S + {x}, read off
the induced pattern of F - x and deg_F(x) = r_x:

    F - x = K4:        12 - 6 deg_F(x) + 12 [deg = 4] - 12 [deg = 0]
                       (from L4, from +12 sum K4(N_v) via L3, from L1)
    F - x = diamond:   +4 if deg_F(x) = 4    (4 sum diamond(N_v), L3)
    F - x = K3 + K1:   -6 if deg_F(x) = 4    (-6 sum s(K3+K1), L3)

The ledger closes pattern by pattern over induced 5-sets:

    F = K5:          each of the 5 vertices gives 12 - 24 + 12 = 0.
                     The K5-type incidence is carried exactly by the
                     +12 sum_v K4(N_v) correction term.
    F = K5 - e:      the two degree-3 vertices give -6 each (their
                     F - x is a K4); the three degree-4 vertices give
                     +4 each (their F - x is the diamond across the
                     missing edge); -12 + 12 = 0.
    F = K4 + K1:     the isolated vertex gives 12 - 12 = 0 (the L1
                     cancellation); the clique vertices see a K3+K1
                     but have degree 3, not 4.
    F = K4 + pendant: the pendant gives +6 (K4, degree 1); the
                     attachment gives -6 (K3+K1, degree 4); the two
                     remaining vertices see paws.
    F = K4 + 2-nbr:  the 2-neighbor gives 12 - 12 = 0; the two
                     opposite clique vertices see diamonds but have
                     degree 3, not 4; the two adjacent ones see paws.

Every 5-set carrying any incidence is one of these: it either contains
a K4 (K4 plus an r-neighbor, r = 0..4, with r = 3 being K5 - e and
r = 4 being K5), or contains a diamond or K3+K1 with a common neighbor
(diamond plus common neighbor = K5 - e; K3+K1 plus common neighbor =
K4 + pendant). All nets are zero, so the universal residual vanishes
on every graph; on K4-free neighborhoods the +12 term is empty and the
specialized row h_v above remains. `tests/test_m4_subgraph_identities.py`
re-proves the closure mechanically: exhaustive enumeration of all
labeled graphs on up to six vertices, the two coefficient-mutation
witnesses, all 656 published R(5,5,42) graphs and complements, and the
MR97 n = 49 replay below.

Exact fast counters
-------------------
diamond(G) = sum over non-edges uv of e(G[A_u & A_v]).  An induced
diamond has a unique non-edge, its degree-2 pair, and for such uv the
quartet {u, v, w, x} spans five edges exactly when w, x is an adjacent
pair of common neighbors of both u and v (a non-adjacent pair spans a
C4 instead). Each diamond is counted once, by its non-edge.

i4(G) = (1/6) sum over non-edges uv of [ C(c, 2) - e(G[W]) ], where W
is the common non-neighborhood of u and v with u, v excluded and
c = |W|. Every pair {a, b} in W is non-adjacent to both u and v, so
{u, v, a, b} is an independent quartet exactly when ab is a non-edge;
each independent quartet contains six non-edges and is counted six
times, and the divisibility by six is asserted. The edge subtraction
is REQUIRED: without it the sum equals 6 i4 + X, where X counts the
induced one-edge quartets (on K2+2K1 it already gives 1 instead of 0);
the research shortcut C(c-2, 2) is provably inexact and is not used.

s(K3+K1, G) = sum over triangles T = {a, b, c} of the number of common
non-neighbors of T outside T. The fourth vertex is non-adjacent to all
of a, b, c, so the quartet spans exactly the triangle; every induced
K3+K1 is counted once, by its triangle.

k4(G) = i4(complement(G)): clique quartets are the independent
quartets of the complement.

Ground-truth anchor (McKay-Radziszowski 1997, Thm 3.2, n = 49,
https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf; replayed on
data/r45_24.g6): exactly two of the 352,366 R(4,5,24) graphs have 132
edges; both have t = 176, diamond = 792, s(K3+K1) = 528, and i4 in
{138, 144}. The n = 49, d = 24 row 9*176 + 4*792 - 6*528 = 1584 =
12*132 forces per-vertex mean i4 = 132, while both admissible
neighborhood types have i4 >= 138 > 132: the identity alone refutes
every 24-regular hypothetical Ramsey(5,5,49) configuration.
"""

from collections import namedtuple

from check_ramsey import popcount
from subgraph_identities import complement, induced_subgraph, triangle_count

Motif4 = namedtuple(
    "Motif4", "order triangles diamonds k3k1 independent_quads cliques"
)


def _pairs_inside(adj, mask):
    """Induced edge count inside the vertex bitmask, via per-vertex rows.

    Triangular counting: each vertex contributes its edges to the strictly
    REMAINING mask members (mask shrinks as it is consumed), so every edge
    inside the mask is counted exactly once.
    """
    edges = 0
    while mask:
        bit = mask & -mask
        v = bit.bit_length() - 1
        edges += popcount(mask & adj[v])
        mask ^= bit
    return edges


def diamond_count(adj):
    """Induced K4-minus-an-edge count, via each diamond's unique non-edge."""
    n = len(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                continue
            total += _pairs_inside(adj, adj[u] & adj[v])
    return total


def independent_quad_count(adj):
    """Independent quartet count i4, via common non-neighborhoods.

    (1/6) * sum over non-edges uv of [C(c, 2) - e(G[W])], W the common
    non-neighborhood of u and v, u and v excluded; the edge subtraction
    is REQUIRED (see the module docstring for the failure mode).
    """
    n = len(adj)
    dual = complement(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                continue
            w = dual[u] & dual[v] & ~((1 << u) | (1 << v))
            c = popcount(w)
            total += c * (c - 1) // 2 - _pairs_inside(adj, w)
    assert total % 6 == 0
    return total // 6


def k4_count(adj):
    """Induced clique quartets: the independent quartets of the complement."""
    return independent_quad_count(complement(adj))


def induced_k3k1_count(adj):
    """Induced K3+K1 count, via triangles and their common non-neighbors."""
    n = len(adj)
    dual = complement(adj)
    total = 0
    for a in range(n):
        later = adj[a] & ~((1 << (a + 1)) - 1)
        while later:
            b1 = later & -later
            b = b1.bit_length() - 1
            tri = adj[a] & adj[b] & ~((1 << (b + 1)) - 1)
            while tri:
                b2 = tri & -tri
                c = b2.bit_length() - 1
                total += popcount(
                    dual[a] & dual[b] & dual[c]
                    & ~((1 << a) | (1 << b) | (1 << c))
                )
                tri ^= b2
            later ^= b1
    return total


def _neighborhood(adj, v):
    """N_v: the subgraph induced on the neighbors of v."""
    return induced_subgraph(
        adj, [u for u in range(len(adj)) if (adj[v] >> u) & 1]
    )


def _dual(adj, v):
    """D_v: the subgraph induced on the non-neighbors of v, v excluded."""
    n = len(adj)
    return induced_subgraph(
        adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
    )


def motif4(adj):
    """All four-vertex induced motif counts of `adj`, plus its order."""
    return Motif4(
        len(adj),
        triangle_count(adj),
        diamond_count(adj),
        induced_k3k1_count(adj),
        independent_quad_count(adj),
        k4_count(adj),
    )


def m4_vertex_row_scaled(n, neighborhood, dual):
    """Scaled specialized row h_v of one vertex (see module docstring)."""
    if len(neighborhood) + len(dual) != n - 1:
        raise ValueError("neighborhood and dual orders must sum to n-1")
    d = len(neighborhood)
    return (
        3 * (n + 2 - 2 * d) * triangle_count(neighborhood)
        + 4 * diamond_count(neighborhood)
        - 6 * induced_k3k1_count(neighborhood)
        - 12 * independent_quad_count(complement(dual))
    )


def m4_residual_specialized(adj):
    """sum_v h_v; zero for every graph with K4-free neighborhoods.

    Raises ValueError when any neighborhood contains a K4: the
    specialized form drops the +12 sum_v K4(N_v) correction term of the
    universal identity and is valid only under that hypothesis.
    """
    n = len(adj)
    total = 0
    for v in range(n):
        nbhd = _neighborhood(adj, v)
        if k4_count(nbhd):
            raise ValueError("specialized identity needs K4-free neighborhoods")
        total += m4_vertex_row_scaled(n, nbhd, _dual(adj, v))
    return total


def m4_residual_universal(adj):
    """Total scaled residual of the universal identity; zero for every graph.

    Per vertex: the specialized row h_v (which already carries the
    -12 K4(D_v) side through i4(complement(D_v)) = K4(D_v)) plus the
    +12 K4(N_v) correction, i.e. exactly the displayed universal row.
    """
    n = len(adj)
    total = 0
    for v in range(n):
        nbhd = _neighborhood(adj, v)
        total += m4_vertex_row_scaled(n, nbhd, _dual(adj, v))
        total += 12 * k4_count(nbhd)  # the K5-type correction, L3
    return total
