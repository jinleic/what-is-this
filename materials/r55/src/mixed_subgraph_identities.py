#!/usr/bin/env python3
"""Exact local Engström K=4 "mixed" subgraph identity for the R(5,5) campaign.

For a graph G on n vertices write X = G_v^+ for the subgraph induced on
the neighborhood of v and Z = complement(G_v^-) for the complement of the
subgraph induced on the non-neighborhood of v (v itself excluded from both),
d_v for the degree of v, and s(J, H) for the number of induced copies of
J in H. Engström's K=4 identity (arXiv:1002.4304, Theorem
"(Conjecture 1 in McKay and Radziszowski)"; VERBATIM coefficients in
notes/engstrom_identity_extraction_2026-08-18.md section 1) is the
universal per-vertex sum

    sum_v [ p_1(G_v^+) + p_2(G_v^-) + p_3(G_v^+, G_v^-) ]  =  0,

with Y = G_v^- (so Z = complement(Y)) and

    p_1(X) = n(n-3) s(K1,X) - (n^2+2n-6) s(K1,X)^2
              + 3n s(K1,X)^3 - 2 s(K1,X)^4
              + 2(n^2+n-8) s(K2,X) - 12 s(K2,X)^2
              - 12(n-1) s(K1,X) s(K2,X) + 12 s(K1,X)^2 s(K2,X)
              + 72 s(C4,X) + 12(n-2) s(K3,X) + 24 s(K1,3,X)
              + 24 s(P4,X) + 24 s(T31,X) + 12(n+2) s(P3,X)
              - 24 s(K1,X) s(P3,X) + 32 s(T32,X),
    p_2(Y) = 4 s(K2,Y)^2 - 12 s(K1,3,Y) - 8 s(C4,Y)
              - 8 s(T31,Y) - 24 s(T32,Y) + 2(n-8) s(P3,Y),
    p_3(X,Y) = 4 s(K1,X) s(P3,Y) - 2(n-2) s(K1,X) s(K2,Y)
              + 4 s(K1,X)^2 s(K2,Y),

where T31 is the paw and T32 the diamond. Unlike the m=3/m=4 rows this
identity needs NO K4-free hypothesis and NO s(K4,*) counter: it holds for
every graph, so the kernel stores exactly one (universal) residual form.

Independent counting derivation (the ledger)
-------------------------------------------
L1 (incidence).  For a 4-vertex pattern J, a copy of J inside N_v is a
copy of J in G, not containing v, all of whose vertices v is adjacent to;
a copy of J inside D_v is a copy of J in G, not containing v, none of
whose vertices v is adjacent to. So with cn(S), cnn(S) the common
neighbors / common non-neighbors of a vertex set S:

    sum_v s(J, N_v) = sum_{S ~= J} |cn(S)|,
    sum_v s(J, D_v) = sum_{S ~= J} |cnn(S)|.

Since Z = complement(D_v) is exactly the complement of the dual, the Y-side
counts s(J, Y) = s(J, D_v) are read off Z through the complement-pair
map of section 5 of the extraction note (each pinned in the tests):

    s(K2,  Y) = e(Z-bar)  = C(mu,2) - e(Z),
    s(P3,  Y) = pc3(Z),
    s(C4,  Y) = 2K2(Z),
    s(K1,3,Y) = s(K3+K1, Z),
    s(T31, Y) = s(P3+K1, Z),
    s(T32, Y) = s(K2+2K1, Z),

where mu = |Y| = |Z| = n-1-d_v.

L2 (cross products).  s(K1,X) = d_v and the p_3 monomials carry
degree-weighted sums:

    sum_v d_v^k s(J, D_v) = sum_{S ~= J} sum_{w in cnn(S)} d_w^k,

k in {1,2}, and dually sum_v d_v^k s(J, N_v) = sum_{S~=J} sum_{w in cn(S)} d_w^k.

L3 (same-side products).  A product s(J1,X) s(J2,X) expands as
sum over (S1 ~= J1, S2 ~= J2) of the number of vertices common to
cn(S1) and cn(S2); regroups by the induced type on S1 u S2.  The K=4
bound keeps every such union inside four vertices, so the regrouping closes on
the eleven 4-vertex patterns.

L4 (closure).  The 25 monomials (16 in p_1, 6 in p_2, 3 in p_3)
regroup through L1-L3 so that every per-incidence weight is carried by a
5-vertex set F = S + {x} (with the roles of cn and cnn, and the
degree weights, folded in); reading the induced pattern of F - x and
deg_F(x) from the four-vertex type table shows every net is zero.  This is
the pattern-by-pattern closure claimed by Engström's Theorem 2.18 /
Corollary 2.19; the exhaustive labeled-graph test below (L5) is the
second, mechanically independent check.

L5 (mechanical).  `tests/test_mixed_subgraph_identities.py` verifies the
residual is identically zero on all 33,867 labeled graphs n <= 6, plus
40 random graphs at n in {7, 8}, the 656 published R(5,5,42) graphs
and complements, and the n=49 type-combo multiset [0, 144, 288, 432]
(unique zero corner) on the four (X,Z) combos of the two 132-edge
R(4,5,24) census graphs.  Each fast closed form below is also
cross-validated against a naive 4-subset enumeration on all labeled graphs
n <= 6, and the corrected-math mutation pins are asserted on their witnesses.

Corrected-math rule (v3-typed).  The exact counter forms supersede six
provably WRONG candidate short forms (each pinned as a detected mutation in the
tests):  claw needs BOTH -2*diamond and -4*K4 (the 4-vertex diamond
has sum C(d,3) = 2 but claw = 0);  s(K1,3, complement) is
s(K3+K1, ·), not paw;  s(C4, complement) is 2K2(·), C4 is not
self-complementary;  pc3 needs the +c_uv term (witness K3: -3 vs 0);
K2+2K1 needs the -2*(2K2) subtraction (witness 2K2);  paw needs the
universal -12*K4 term (witness K4: 12 vs 0).  On the swept R(4,5)
catalogs k4 == 0 by Ramseyhood, asserted during streaming, never assumed.
"""

from collections import namedtuple

from check_ramsey import popcount
from subgraph_identities import (
    complement,
    induced_subgraph,
    induced_p3_count,
    motif3,
    triangle_count,
)
from m4_subgraph_identities import (
    _dual,
    _neighborhood,
    _pairs_inside,
    diamond_count,
    independent_quad_count,
    induced_k3k1_count,
    k4_count,
    motif4,
)

MotifSweep = namedtuple(
    "MotifSweep",
    "order edges wedges triangles induced_p3 independent_triples "
    "e4 k2_2k1 two_k2 p3_k1 p4 claw k3k1 c4 paw diamond k4 pc3",
)


def _c2(x):
    return x * (x - 1) // 2


def _c3(x):
    return x * (x - 1) * (x - 2) // 6


def _c4(x):
    return x * (x - 1) * (x - 2) * (x - 3) // 24


def paw_count(adj):
    """Induced paws (T31): sum_v (d_v - 2) t_v - 4*diamond - 12*K4.

    t_v counts triangles through v, one per triangle via the edges inside
    N_v (_pairs_inside).  The -4*diamond and -12*K4 corrections are
    universal; on Ramsey catalogs k4 == 0 (asserted at stream time).
    """
    n = len(adj)
    tri_through = [_pairs_inside(adj, adj[v]) for v in range(n)]
    total = sum((popcount(adj[v]) - 2) * tv for v, tv in enumerate(tri_through))
    return total - 4 * diamond_count(adj) - 12 * k4_count(adj)


def claw_count(adj):
    """Induced K_{1,3}: sum_v C(d_v,3) - paw - 2*diamond - 4*K4."""
    return sum(_c3(popcount(r)) for r in adj) - (
        paw_count(adj) + 2 * diamond_count(adj) + 4 * k4_count(adj)
    )


def c4_count(adj):
    """Induced C4: (1/2) * sum over non-edges uv of [C(c_uv,2) - e(G[cn])],
    cn = the common-neighbor set of u and v (any non-edge of cn completes a
    C4 with uv).  Each C4 is counted exactly twice, once per non-edge.  The
    e(G[cn]) subtraction removes the diamond completion (an adjacent pair inside cn
    makes uvwx a diamond, not a C4).  A threshold form ("non-edges with
    exactly two common neighbors") is WRONG: a non-edge can have c > 2 common
    neighbors and complete C(c,2) - e(G[cn]) distinct C4s.
    """
    n = len(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                continue
            cn = adj[u] & adj[v]
            c = popcount(cn)
            total += c * (c - 1) // 2 - _pairs_inside(adj, cn)
    return total // 2


def two_k2_count(adj):
    """Induced 2K2: half of sum over edges uv of edges inside the common
    NON-neighborhood of u and v (u, v excluded).  Each 2K2 is counted
    once from each of its two edges."""
    n = len(adj)
    full = (1 << n) - 1
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if not ((adj[u] >> v) & 1):
                continue
            mask = full & ~adj[u] & ~adj[v] & ~((1 << u) | (1 << v))
            total += _pairs_inside(adj, mask)
    return total // 2


def k2_2k1_count(adj):
    """Induced K2 + 2K1: sum over edges uv of C(nu - |N[u] u N[v]|, 2)
    minus 2*(2K2).  The subtraction is REQUIRED: each 2K2 is counted
    twice here (once from each edge's outside pair) and has no isolated
    vertex relative to the edge."""
    n = len(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if not ((adj[u] >> v) & 1):
                continue
            closed = popcount(adj[u] | adj[v])  # |N[u] u N[v]|, no +2: u,v in it
            outside = n - closed
            total += _c2(outside)
    return total - 2 * two_k2_count(adj)


def p3_k1_count(adj):
    """Induced P3 + K1: sum over x of induced_P3(G - N[x]).  Each P3+K1
    is counted once, by its isolated vertex x."""
    n = len(adj)
    total = 0
    for x in range(n):
        keep = [u for u in range(n) if u != x and not ((adj[x] >> u) & 1)]
        total += induced_p3_count(induced_subgraph(adj, keep))
    return total


def p4_count(adj):
    """Induced P4 on four vertices: C(nu,4) minus the other ten
    induced 4-vertex classes (the eleven classes partition all 4-subsets)."""
    n = len(adj)
    other = (
        independent_quad_count(adj)
        + k2_2k1_count(adj)
        + two_k2_count(adj)
        + p3_k1_count(adj)
        + induced_k3k1_count(adj)
        + claw_count(adj)
        + c4_count(adj)
        + paw_count(adj)
        + diamond_count(adj)
        + k4_count(adj)
    )
    return _c4(n) - other


def complement_p3_count(adj):
    """s(P3, complement(adj)); form B: e(nu-2) - 2*p3 - 3*t.

    Provably equivalent to form A sum_{uv in E}(nu - d_u - d_v + c_uv)
    and form C sum_v C(nu-1-d_v, 2) - 3*i3; all three are enforced
    equal exhaustively in the tests."""
    n = len(adj)
    m3 = motif3(adj)
    return m3.edges * (n - 2) - 2 * m3.induced_p3 - 3 * m3.triangles


def _counts4(adj):
    """The eleven raw induced 4-vertex classes of `adj` in ONE pass.

    Every sub-quantity is computed exactly once (the naive "derive p4 from
    the other ten" architecture re-prices diamond/paw/K4/C4/paw several
    times).  The complement of `adj` is materialized once and reused by
    k4 (= i4 of the complement) and p3_k1 (= paw of the complement), whose
    diamond/K4 corrections then come for free from the already-computed
    complement-pair identities s(diamond, comp Z) = k2_2k1(Z) and
    s(K4, comp Z) = i4(Z) (both machine-verified exhaustively n <= 6):

        p3_k1(Z) = [sum_v (d_comp(v)-2) t_v(comp)] - 4*k2_2k1(Z) - 12*i4(Z).
    """
    n = len(adj)
    cj = complement(adj)
    dia = diamond_count(adj)
    i4v = independent_quad_count(adj)
    k4v = independent_quad_count(cj)
    k3k1v = induced_k3k1_count(adj)
    t2k2 = two_k2_count(adj)
    c4v = c4_count(adj)
    deg = [popcount(r) for r in adj]
    pawv = (sum((d - 2) * _pairs_inside(adj, adj[v])
                for v, d in enumerate(deg))
            - 4 * dia - 12 * k4v)
    clawv = sum(_c3(d) for d in deg) - pawv - 2 * dia - 4 * k4v
    closed_sum = 0
    for u in range(n):
        for v in range(u + 1, n):
            if not ((adj[u] >> v) & 1):
                continue
            outside = n - popcount(adj[u] | adj[v])
            closed_sum += _c2(outside)
    k2_2k1v = closed_sum - 2 * t2k2
    p3k1v = (sum((popcount(cj[v]) - 2) * _pairs_inside(cj, cj[v])
                 for v in range(n))
             - 4 * k2_2k1v - 12 * i4v)
    p4v = _c4(n) - (i4v + k2_2k1v + t2k2 + p3k1v + k3k1v + clawv
                    + c4v + pawv + dia + k4v)
    return {
        "e4": i4v, "k2_2k1": k2_2k1v, "two_k2": t2k2, "p3_k1": p3k1v,
        "p4": p4v, "claw": clawv, "k3k1": k3k1v, "c4": c4v,
        "paw": pawv, "diamond": dia, "k4": k4v,
    }


def motif_sweep(adj):
    """All fifteen aggregates used by the mixed row plus raw 3-motifs, in one
    pass over a graph recording no intermediate structures."""
    m3 = motif3(adj)
    c4map = _counts4(adj)
    return MotifSweep(
        order=len(adj),
        edges=m3.edges,
        wedges=sum(_c2(popcount(r)) for r in adj),
        triangles=m3.triangles,
        induced_p3=m3.induced_p3,
        independent_triples=m3.independent_triples,
        e4=c4map["e4"],
        k2_2k1=c4map["k2_2k1"],
        two_k2=c4map["two_k2"],
        p3_k1=c4map["p3_k1"],
        p4=c4map["p4"],
        claw=c4map["claw"],
        k3k1=c4map["k3k1"],
        c4=c4map["c4"],
        paw=c4map["paw"],
        diamond=c4map["diamond"],
        k4=c4map["k4"],
        pc3=m3.edges * (len(adj) - 2) - 2 * m3.induced_p3 - 3 * m3.triangles,
    )


def mixed_vertex_row(n, d, x, z):
    """Engström F(X, Z) at one vertex, exact integer.

    X = G_v^+ of order d; Z = complement(G_v^-) of order m = n-1-d.
    The Y-side (Y = G_v^- = complement(Z)) is single-sourced through the
    complement-pair map, so only Z's sweep is passed in.
    """
    if x.order != d or x.order + z.order != n - 1:
        raise ValueError("neighborhood and complemented-dual orders must sum to n-1")
    ey = _c2(n - 1 - d) - z.edges  # s(K2, Y), Y = complement(Z)
    exact = (
        n * (n - 3) * d
        - (n * n + 2 * n - 6) * d * d
        + 3 * n * d ** 3
        - 2 * d ** 4
        + 2 * (n * n + n - 8) * x.edges
        - 12 * x.edges * x.edges
        - 12 * (n - 1) * d * x.edges
        + 12 * d * d * x.edges
        + 4 * ey * ey
        - 2 * (n - 2) * d * ey
        + 4 * d * d * ey
    )
    x_motifs = (
        72 * x.c4
        + 12 * (n - 2) * x.triangles
        + 24 * x.claw
        + 24 * x.p4
        + 24 * x.paw
        + (12 * (n + 2) - 24 * d) * x.induced_p3
        + 32 * x.diamond
    )
    z_motifs = (
        (2 * (n - 8) + 4 * d) * z.pc3
        - 12 * z.k3k1
        - 8 * z.two_k2
        - 8 * z.p3_k1
        - 24 * z.k2_2k1
    )
    return exact + x_motifs + z_motifs


def mixed_residual(adj):
    """sum_v F_v for an n-vertex graph; exactly zero for every graph
    (universal identity, no K4-free hypothesis needed)."""
    n = len(adj)
    total = 0
    for v in range(n):
        d_v = popcount(adj[v])
        total += mixed_vertex_row(
            n,
            d_v,
            motif_sweep(_neighborhood(adj, v)),
            motif_sweep(complement(_dual(adj, v))),
        )
    return total
