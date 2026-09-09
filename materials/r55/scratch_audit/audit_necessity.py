"""Necessity audit: for every graph in the published strata, decompose at
EVERY vertex (frame claims) and at every max-degree vertex (constraint claims)
and verify every spec claim holds. Any violation = spec would miss real graphs.
"""
import os
import sys
from itertools import combinations
from math import ceil

from g6lib import (parse_g6, ecount, has_clique, has_indep, indep_sets,
                   cliques, induced, load_g6_file)

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'r45extreme')

STRATA = [  # (n, e, expected iso classes)
    (12, 48, 1), (13, 53, 2), (13, 52, 10), (14, 60, 1), (15, 66, 1),
    (16, 72, 5), (16, 71, 138), (17, 79, 1), (17, 78, 86), (21, 107, 31),
]

fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL:", msg)


def audit_graph(tag, n, adj, e_expected):
    # 0. sanity: the published graph really is in R(4,5,n,e)
    e = ecount(n, adj)
    if e != e_expected:
        fail(f"{tag}: edge count {e} != {e_expected}")
    if has_clique(n, adj, 4):
        fail(f"{tag}: published graph has a K4")
    if has_indep(n, adj, 5):
        fail(f"{tag}: published graph has an I5")

    degs = [bin(a).count("1") for a in adj]
    Delta = max(degs)

    # spec Delta-range claim
    lo, hi = ceil(2 * e / n), min(13, n - 1)
    if not (lo <= Delta <= hi):
        fail(f"{tag}: Delta={Delta} outside spec range [{lo},{hi}]")

    # frame claims: for ANY vertex v
    for v in range(n):
        d = degs[v]
        Hverts = [u for u in range(n) if (adj[v] >> u) & 1]
        Kverts = [u for u in range(n) if u != v and not (adj[v] >> u) & 1]
        assert len(Hverts) == d and len(Kverts) == n - 1 - d
        hn, hadj = induced(adj, Hverts)
        kn, kadj = induced(adj, Kverts)
        if has_clique(hn, hadj, 3):
            fail(f"{tag} v={v}: H has a triangle")
        if has_indep(hn, hadj, 5):
            fail(f"{tag} v={v}: H has an I5")
        if has_clique(kn, kadj, 4):
            fail(f"{tag} v={v}: K has a K4")
        if has_indep(kn, kadj, 4):
            fail(f"{tag} v={v}: K has an I4")
        if d > 13:
            fail(f"{tag} v={v}: d={d} > 13")
        if n - 1 - d > 17:
            fail(f"{tag} v={v}: |K|={n-1-d} > 17")

    # constraint claims at every max-degree vertex
    for v in range(n):
        if degs[v] != Delta:
            continue
        d = Delta
        Hverts = [u for u in range(n) if (adj[v] >> u) & 1]
        Kverts = [u for u in range(n) if u != v and not (adj[v] >> u) & 1]
        hn, hadj = induced(adj, Hverts)
        kn, kadj = induced(adj, Kverts)
        hpos = {u: i for i, u in enumerate(Hverts)}
        # cone: S[k] = bitmask over H-positions
        S = []
        for kv in Kverts:
            m = 0
            for u in Hverts:
                if (adj[kv] >> u) & 1:
                    m |= 1 << hpos[u]
            S.append(m)

        cone = sum(bin(s).count("1") for s in S)
        eH, eK = ecount(hn, hadj), ecount(kn, kadj)
        # edge budget identity + range
        if cone != e_expected - d - eH - eK:
            fail(f"{tag} v={v}: cone {cone} != e-d-eH-eK "
                 f"{e_expected - d - eH - eK}")
        if not (0 <= cone <= d * (n - 1 - d)):
            fail(f"{tag} v={v}: cone {cone} outside [0, d*(n-1-d)]")

        # degree-cap arithmetic (identity AND <= d)
        conecnt = [0] * hn
        for s in S:
            for i in range(hn):
                if (s >> i) & 1:
                    conecnt[i] += 1
        for i, u in enumerate(Hverts):
            degH = bin(hadj[i]).count("1")
            lhs = 1 + degH + conecnt[i]
            if lhs != degs[u]:
                fail(f"{tag} v={v}: deg identity h={u}: 1+{degH}+{conecnt[i]}"
                     f" != {degs[u]}")
            if lhs > d:
                fail(f"{tag} v={v}: deg cap violated at h={u}: {lhs} > {d}")
        for j, u in enumerate(Kverts):
            degK = bin(kadj[j]).count("1")
            lhs = degK + bin(S[j]).count("1")
            if lhs != degs[u]:
                fail(f"{tag} v={v}: deg identity k={u}: {degK}+|S| != {degs[u]}")
            if lhs > d:
                fail(f"{tag} v={v}: deg cap violated at k={u}: {lhs} > {d}")

        # precompute H independent sets
        I2 = [frozenset(t) for t in indep_sets(hn, hadj, 2)]
        I3 = [frozenset(t) for t in indep_sets(hn, hadj, 3)]
        I4 = [frozenset(t) for t in indep_sets(hn, hadj, 4)]
        Hedges = [(a, b) for a in range(hn) for b in range(a + 1, hn)
                  if (hadj[a] >> b) & 1]

        def mask2set(m):
            return frozenset(i for i in range(hn) if (m >> i) & 1)

        # (2,2): every K-edge: S_k1 & S_k2 independent in H
        for a in range(kn):
            for b in range(a + 1, kn):
                if (kadj[a] >> b) & 1:
                    common = S[a] & S[b]
                    for (x, y) in Hedges:
                        if (common >> x) & 1 and (common >> y) & 1:
                            fail(f"{tag} v={v}: (2,2) violated K-edge "
                                 f"({a},{b}) H-edge ({x},{y})")
        # (1,3): every K-triangle: triple intersection empty
        for t in cliques(kn, kadj, 3):
            if S[t[0]] & S[t[1]] & S[t[2]]:
                fail(f"{tag} v={v}: (1,3) violated K-triangle {t}")
        # (4,1): every S_k hits every I4(H)
        for j in range(kn):
            sk = mask2set(S[j])
            for q in I4:
                if not (sk & q):
                    fail(f"{tag} v={v}: (4,1) violated k={j} I4={q}")
        # (3,2): every K-independent pair: union hits every I3(H)
        for a in range(kn):
            for b in range(a + 1, kn):
                if not (kadj[a] >> b) & 1:
                    un = mask2set(S[a] | S[b])
                    for t in I3:
                        if not (un & t):
                            fail(f"{tag} v={v}: (3,2) violated pair ({a},{b})"
                                 f" I3={t}")
        # (2,3): every K-independent triple: union hits every I2(H)
        for t in indep_sets(kn, kadj, 3):
            un = mask2set(S[t[0]] | S[t[1]] | S[t[2]])
            for p in I2:
                if not (un & p):
                    fail(f"{tag} v={v}: (2,3) violated triple {t} I2={p}")


def main():
    total = 0
    for n, e, expected in STRATA:
        path = f"{DATA}/r45{n}.{e}.g6"
        graphs = load_g6_file(path)
        if len(graphs) != expected:
            fail(f"{path}: {len(graphs)} lines, expected {expected}")
        for gi, (gn, adj) in enumerate(graphs):
            if gn != n:
                fail(f"{path}#{gi}: n={gn} != {n}")
            audit_graph(f"r45{n}.{e}#{gi}", gn, adj, e)
            total += 1
        print(f"done r45{n}.{e}: {len(graphs)} graphs")
    print(f"\naudited {total} graphs; {len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
