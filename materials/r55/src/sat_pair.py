#!/usr/bin/env python3
"""SAT path for census gluing pairs that exceed the DFS node budget.

Builds the spec's constraint set as CNF (variables x[h][k] = cone edge) with
cardinality constraints (total = C_t, per-vertex degree windows) and
enumerates ALL solutions with blocking clauses over the x-variables (aux vars
from cardinality encodings are not blocked, so each cone appears exactly
once). Every solution is re-verified with the trusted checker before output.

Called by run_stratum_c.py with a todo list; also usable standalone.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import (parse_graph6_line, encode_graph6,
                          check_ramsey_graph, popcount)  # noqa: E402

from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
from pysat.solvers import Cadical195  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

_HCACHE = {}


def h_structs(d, h_idx, adj):
    key = (d, h_idx)
    if key in _HCACHE:
        return _HCACHE[key]
    full = (1 << d) - 1
    closed = [adj[v] | (1 << v) for v in range(d)]
    mis = bytearray(1 << d)
    for m in range(1, 1 << d):
        v = (m & -m).bit_length() - 1
        a = mis[m & (m - 1)]
        b = 1 + mis[m & ~closed[v]]
        mis[m] = b if b > a else a
    isets = {2: [], 3: [], 4: []}
    for m in range(1 << d):
        pc = popcount(m)
        if pc in isets and mis[m] == pc:  # m is independent
            isets[pc].append([i for i in range(d) if m >> i & 1])
    edges = [(a, b) for a in range(d) for b in range(a + 1, d)
             if adj[a] >> b & 1]
    _HCACHE[key] = (edges, isets)
    return _HCACHE[key]


def solve_pair(n, e, d, H_adj, K_adj, h_idx, dmin, stats=None):
    """Enumerate all valid cones for (H, K); return list of glued adj lists."""
    q = len(K_adj)
    eH = sum(popcount(a) for a in H_adj) // 2
    eK = sum(popcount(a) for a in K_adj) // 2
    c_t = e - d - eH - eK
    if c_t < 0 or c_t > d * q:
        return []
    h_edges, h_isets = h_structs(d, h_idx, H_adj)
    k_deg = [popcount(a) for a in K_adj]
    caps = [d - k_deg[k] for k in range(q)]
    lov = [max(0, dmin - k_deg[k]) for k in range(q)]
    caph = [d - 1 - popcount(H_adj[h]) for h in range(d)]
    need = [max(0, dmin - 1 - popcount(H_adj[h])) for h in range(d)]
    if any(c < 0 for c in caps) or any(need[h] > caph[h] for h in range(d)):
        return []

    def x(h, k):
        return h * q + k + 1

    pool = IDPool(start_from=d * q + 1)
    cls = []
    k_edges, k_ipairs, k_tris, k_itris = [], [], [], []
    for a in range(q):
        for b in range(a + 1, q):
            (k_edges if K_adj[a] >> b & 1 else k_ipairs).append((a, b))
            for c in range(b + 1, q):
                ab, ac, bc = (K_adj[a] >> b & 1, K_adj[a] >> c & 1,
                              K_adj[b] >> c & 1)
                if ab and ac and bc:
                    k_tris.append((a, b, c))
                elif not (ab or ac or bc):
                    k_itris.append((a, b, c))
    for (h1, h2) in h_edges:                      # (2,2)
        for (k1, k2) in k_edges:
            cls.append([-x(h1, k1), -x(h1, k2), -x(h2, k1), -x(h2, k2)])
    for (k1, k2, k3) in k_tris:                   # (1,3)
        for h in range(d):
            cls.append([-x(h, k1), -x(h, k2), -x(h, k3)])
    for quad in h_isets[4]:                       # (4,1)
        for k in range(q):
            cls.append([x(h, k) for h in quad])
    for tri in h_isets[3]:                        # (3,2)
        for (k1, k2) in k_ipairs:
            cls.append([x(h, ki) for h in tri for ki in (k1, k2)])
    for pair in h_isets[2]:                       # (2,3)
        for (k1, k2, k3) in k_itris:
            cls.append([x(h, ki) for h in pair for ki in (k1, k2, k3)])
    allx = [x(h, k) for h in range(d) for k in range(q)]
    cls.extend(CardEnc.equals(allx, c_t, vpool=pool,
                              encoding=EncType.seqcounter).clauses)
    for k in range(q):
        col = [x(h, k) for h in range(d)]
        if caps[k] < d:
            cls.extend(CardEnc.atmost(col, caps[k], vpool=pool,
                                      encoding=EncType.seqcounter).clauses)
        if lov[k] > 0:
            cls.extend(CardEnc.atleast(col, lov[k], vpool=pool,
                                       encoding=EncType.seqcounter).clauses)
    for h in range(d):
        row = [x(h, k) for k in range(q)]
        if caph[h] < q:
            cls.extend(CardEnc.atmost(row, caph[h], vpool=pool,
                                      encoding=EncType.seqcounter).clauses)
        if need[h] > 0:
            cls.extend(CardEnc.atleast(row, need[h], vpool=pool,
                                       encoding=EncType.seqcounter).clauses)
    out = []
    with Cadical195(bootstrap_with=cls) as solver:
        while solver.solve():
            model = solver.get_model()
            pos = {v for v in model if 0 < v <= d * q}
            # build glued graph: 0 = v, 1..d = H, d+1.. = K
            adj = [0] * n
            for h in range(d):
                adj[0] |= 1 << (1 + h)
                adj[1 + h] |= 1
                a = H_adj[h]
                while a:
                    b = a & -a
                    adj[1 + h] |= 1 << (1 + (b.bit_length() - 1))
                    a ^= b
            for k in range(q):
                a = K_adj[k]
                while a:
                    b = a & -a
                    adj[1 + d + k] |= 1 << (1 + d + (b.bit_length() - 1))
                    a ^= b
            for h in range(d):
                for k in range(q):
                    if x(h, k) in pos:
                        adj[1 + h] |= 1 << (1 + d + k)
                        adj[1 + d + k] |= 1 << (1 + h)
            assert sum(popcount(a) for a in adj) // 2 == e
            assert check_ramsey_graph(n, adj, 4, 5) == (True, True), \
                "SAT solution failed trusted verification"
            out.append(encode_graph6(n, adj))
            solver.add_clause([-x(h, k) if x(h, k) in pos else x(h, k)
                               for h in range(d) for k in range(q)])
        if stats is not None:
            acc = solver.accum_stats()
            for key in ("conflicts", "decisions", "propagations"):
                stats[key] = stats.get(key, 0) + acc.get(key, 0)
    return out


_CATS = {}


def catalog(kind, order):
    key = (kind, order)
    if key not in _CATS:
        if order == 0:
            _CATS[key] = [()]
        else:
            out = []
            for line in open(os.path.join(DATA, f"{kind}_{order}.g6")):
                if line.strip():
                    nn, adj = parse_graph6_line(line)
                    out.append(tuple(adj))
            _CATS[key] = out
    return _CATS[key]


def run_todo(n, e, todos, dmin):
    """todos: list of (d, h_idx, k_idx). Returns (g6lines, counts rows)."""
    lines, rows = [], []
    for d, h_idx, k_idx in todos:
        H = catalog("r35", d)[h_idx]
        K = catalog("r44", n - 1 - d)[k_idx]
        sols = solve_pair(n, e, d, H, K, h_idx, dmin)
        for s in sols:
            lines.append(s)
        if sols:
            rows.append((d, h_idx, k_idx, len(sols)))
    return lines, rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--e", type=int, required=True)
    ap.add_argument("--dmin", type=int, default=0)
    ap.add_argument("--pairs", required=True,
                    help="comma-separated d:h:k triples or @file with lines d,h,k")
    args = ap.parse_args()
    if args.pairs.startswith("@"):
        todos = [tuple(int(x) for x in ln.split(",")[:3])
                 for ln in open(args.pairs[1:]) if ln.strip()]
    else:
        todos = [tuple(int(x) for x in p.split(":")) for p in args.pairs.split(",")]
    lines, rows = run_todo(args.n, args.e, todos, args.dmin)
    for ln in lines:
        print(ln)
    for r in rows:
        print("#", ",".join(map(str, r)), file=sys.stderr)
