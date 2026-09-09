"""For every published graph in selected strata and EVERY max-degree vertex:
H = G[N(v)] must be isomorphic to some line of r35_d.g6 and
K = G[rest] isomorphic to some line of r44_q.g6.
Own backtracking isomorphism test (invariant prefilter + refinement).
"""
import os
import sys
from itertools import combinations
from g6lib import load_g6_file, induced, ecount

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL:", msg)


def invariant(n, adj):
    degs = sorted(bin(a).count("1") for a in adj)
    tri = 0
    for a in range(n):
        for b in range(a + 1, n):
            if (adj[a] >> b) & 1:
                tri += bin(adj[a] & adj[b]).count("1")
    return (n, tuple(degs), tri)  # tri = 3*#triangles


def vertex_inv(n, adj, v):
    d = bin(adj[v]).count("1")
    nbr_degs = tuple(sorted(bin(adj[u]).count("1")
                            for u in range(n) if (adj[v] >> u) & 1))
    return (d, nbr_degs)


def isomorphic(n1, a1, n2, a2):
    if n1 != n2:
        return False
    if invariant(n1, a1) != invariant(n2, a2):
        return False
    n = n1
    inv1 = [vertex_inv(n, a1, v) for v in range(n)]
    inv2 = [vertex_inv(n, a2, v) for v in range(n)]
    if sorted(inv1) != sorted(inv2):
        return False
    # backtracking: map vertices of g1 (ordered by rarity of invariant) to g2
    order = sorted(range(n), key=lambda v: (inv1.count(inv1[v]), v))
    mapping = [-1] * n     # g1 -> g2
    used = [False] * n

    def bt(i):
        if i == n:
            return True
        v = order[i]
        for w in range(n):
            if used[w] or inv1[v] != inv2[w]:
                continue
            ok = True
            for j in range(i):
                u = order[j]
                if ((a1[v] >> u) & 1) != ((a2[w] >> mapping[u]) & 1):
                    ok = False
                    break
            if ok:
                mapping[v] = w
                used[w] = True
                if bt(i + 1):
                    return True
                used[w] = False
                mapping[v] = -1
        return False

    return bt(0)


def main():
    strata = [(12, 48), (13, 53), (13, 52), (16, 71), (17, 78), (21, 107)]
    r35 = {}
    r44 = {}
    checked = 0
    for n, e in strata:
        graphs = load_g6_file(f"{DATA}/r45extreme/r45{n}.{e}.g6")
        for gi, (gn, adj) in enumerate(graphs):
            degs = [bin(a).count("1") for a in adj]
            Delta = max(degs)
            for v in range(gn):
                if degs[v] != Delta:
                    continue
                d = Delta
                q = gn - 1 - d
                Hverts = [u for u in range(gn) if (adj[v] >> u) & 1]
                Kverts = [u for u in range(gn)
                          if u != v and not (adj[v] >> u) & 1]
                hn, hadj = induced(adj, Hverts)
                kn, kadj = induced(adj, Kverts)
                if d not in r35:
                    r35[d] = load_g6_file(f"{DATA}/r35_{d}.g6")
                if q not in r44:
                    r44[q] = (load_g6_file(f"{DATA}/r44_{q}.g6")
                              if q > 0 else [(0, [])])
                hinv = invariant(hn, hadj)
                if not any(isomorphic(hn, hadj, cn, cadj)
                           for cn, cadj in r35[d]
                           if invariant(cn, cadj) == hinv):
                    fail(f"r45{n}.{e}#{gi} v={v}: H (d={d}) NOT in r35_{d}")
                kinv = invariant(kn, kadj)
                if q == 0:
                    ok = (kn == 0)
                else:
                    ok = any(isomorphic(kn, kadj, cn, cadj)
                             for cn, cadj in r44[q]
                             if invariant(cn, cadj) == kinv)
                if not ok:
                    fail(f"r45{n}.{e}#{gi} v={v}: K (q={q}) NOT in r44_{q}")
                checked += 1
        print(f"r45{n}.{e}: membership checked")
    print(f"\n{checked} (v, H, K) decompositions checked; {len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
