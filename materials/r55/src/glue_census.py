#!/usr/bin/env python3
"""Implementation A of the census gluing spec (notes/census_gluing_spec.md).

Enumerates all Ramsey(4,5,n)-graphs with exactly e edges by gluing the
complete R(3,5,d) catalogs (neighborhood of a maximum-degree vertex v) to the
complete R(4,4,n-1-d) catalogs (non-neighborhood), over all valid bipartite
cones, for d in [ceil(2e/n), min(13, n-1)] with all degrees capped at d.

Output contract per the spec: vertex 0 = v, 1..d = H in catalog order,
d+1..n-1 = K in catalog order; raw .g6 (duplicates allowed) + counts.csv
rows (d, h_idx, k_idx, n_solutions). Every emitted graph is re-verified by
the independent checker in check_ramsey.py before being written.
"""

import argparse
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import (parse_graph6_line, encode_graph6,
                          check_ramsey_graph, popcount)  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load_catalog(kind, order):
    """kind in {'r35','r44'}; returns list of adj-bitmask tuples."""
    if order == 0:
        return [()]  # the empty graph
    path = os.path.join(DATA, f"{kind}_{order}.g6")
    out = []
    for line in open(path):
        if line.strip():
            n, adj = parse_graph6_line(line)
            assert n == order
            out.append(tuple(adj))
    return out


class HSide:
    """Per-H precomputation: MIS-size DP, hitting/independence tables,
    candidate cone-neighborhoods grouped by popcount."""

    def __init__(self, adj):
        d = len(adj)
        self.d = d
        self.adj = adj
        full = (1 << d) - 1
        self.full = full
        closed = [adj[v] | (1 << v) for v in range(d)]
        mis = bytearray(1 << d)
        hasedge = bytearray(1 << d)
        for m in range(1, 1 << d):
            v = (m & -m).bit_length() - 1
            rest = m & (m - 1)
            a = mis[rest]
            b = 1 + mis[m & ~closed[v]]
            mis[m] = b if b > a else a
            hasedge[m] = 1 if (hasedge[rest] or (adj[v] & rest)) else 0
        self.mis = mis
        # independence of a subset of V(H) (for the (2,2) K4 constraint)
        self.indep_ok = bytearray(1 if not hasedge[m] else 0
                                  for m in range(1 << d))
        # hitting tables: hit_r_ok[m] == 1 iff m intersects every independent
        # r-subset of H, i.e. mis(V \ m) <= r-1
        self.hit4_ok = bytearray(1 if mis[full & ~m] <= 3 else 0
                                 for m in range(1 << d))
        self.hit3_ok = bytearray(1 if mis[full & ~m] <= 2 else 0
                                 for m in range(1 << d))
        self.hit2_ok = bytearray(1 if mis[full & ~m] <= 1 else 0
                                 for m in range(1 << d))
        # candidate S_k masks (satisfy the (4,1) constraint), by popcount
        self.by_pc = [[] for _ in range(d + 1)]
        for m in range(1 << d):
            if self.hit4_ok[m]:
                self.by_pc[popcount(m)].append(m)
        self.minp = next((p for p in range(d + 1) if self.by_pc[p]), None)
        self.edges = sum(popcount(a) for a in adj) // 2
        self.degs = [popcount(a) for a in adj]
        # cap on cone-degree of h: deg_G(h) = 1 + deg_H(h) + used(h) <= d
        self.caph = [d - 1 - self.degs[h] for h in range(d)]


class KSide:
    """Per-K precomputation: edges/triangles/independent pairs and triples."""

    def __init__(self, adj):
        q = len(adj)
        self.q = q
        self.adj = adj
        self.edges = sum(popcount(a) for a in adj) // 2
        self.degs = [popcount(a) for a in adj]
        self.epairs, self.ipairs = [], []
        self.triangles, self.itriples = [], []
        for a in range(q):
            for b in range(a + 1, q):
                if adj[a] >> b & 1:
                    self.epairs.append((a, b))
                else:
                    self.ipairs.append((a, b))
                for c in range(b + 1, q):
                    ab = adj[a] >> b & 1
                    ac = adj[a] >> c & 1
                    bc = adj[b] >> c & 1
                    if ab and ac and bc:
                        self.triangles.append((a, b, c))
                    elif not (ab or ac or bc):
                        self.itriples.append((a, b, c))


def solve_pair(H, K, d, e_target, emit):
    """Enumerate all valid cones for one (H, K) pair; call emit(S_list) for
    each solution, where S_list[k] is the cone mask of catalog K-vertex k.
    Returns number of solutions."""
    q = K.q
    c_target = e_target - d - H.edges - K.edges
    if c_target < 0 or c_target > d * q:
        return 0
    if q == 0:
        if c_target == 0:
            emit([])
            return 1
        return 0
    if H.minp is None:
        return 0  # no S_k can hit all I4s of H (cannot happen for d<=13)
    caps = [d - K.degs[k] for k in range(q)]
    if min(caps) < H.minp:
        return 0
    # DFS over K vertices, tightest cap first
    perm = sorted(range(q), key=lambda k: caps[k])
    pos = {k: i for i, k in enumerate(perm)}
    # static budget bounds
    sufmin = [0] * (q + 1)
    sufmax = [0] * (q + 1)
    for i in range(q - 1, -1, -1):
        sufmin[i] = sufmin[i + 1] + H.minp
        sufmax[i] = sufmax[i + 1] + caps[perm[i]]
    if not (sufmin[0] <= c_target <= sufmax[0]):
        return 0
    # constraint groups: each keyed by the latest DFS position it involves
    checks = [[] for _ in range(q)]
    for a, b in K.epairs:
        i, j = pos[a], pos[b]
        checks[max(i, j)].append((0, min(i, j)))          # (2,2): indep cap
    for a, b in K.ipairs:
        i, j = pos[a], pos[b]
        checks[max(i, j)].append((1, min(i, j)))          # (3,2): hit3 on union
    for a, b, c in K.triangles:
        i, j, l = sorted((pos[a], pos[b], pos[c]))
        checks[l].append((2, i, j))                       # (1,3): empty triple int.
    for a, b, c in K.itriples:
        i, j, l = sorted((pos[a], pos[b], pos[c]))
        checks[l].append((3, i, j))                       # (2,3): hit2 on union
    indep_ok, hit3_ok, hit2_ok = H.indep_ok, H.hit3_ok, H.hit2_ok
    by_pc, caph = H.by_pc, H.caph
    d_h = H.d
    used = [0] * d_h
    S = [0] * q
    nsol = 0

    def dfs(i, remaining, fullmask):
        nonlocal nsol
        if i == q:
            out = [0] * q
            for j in range(q):
                out[perm[j]] = S[j]
            emit(out)
            nsol += 1
            return
        k = perm[i]
        hi = min(caps[k], remaining - sufmin[i + 1])
        lo = max(H.minp, remaining - sufmax[i + 1])
        my_checks = checks[i]
        for p in range(hi, lo - 1, -1):
            for m in by_pc[p]:
                if m & fullmask:
                    continue
                ok = True
                for chk in my_checks:
                    t = chk[0]
                    if t == 0:
                        if not indep_ok[S[chk[1]] & m]:
                            ok = False
                            break
                    elif t == 1:
                        if not hit3_ok[S[chk[1]] | m]:
                            ok = False
                            break
                    elif t == 2:
                        if S[chk[1]] & S[chk[2]] & m:
                            ok = False
                            break
                    else:
                        if not hit2_ok[S[chk[1]] | S[chk[2]] | m]:
                            ok = False
                            break
                if not ok:
                    continue
                S[i] = m
                newfull = fullmask
                mm = m
                bad = False
                while mm:
                    b = mm & -mm
                    h = b.bit_length() - 1
                    used[h] += 1
                    if used[h] > caph[h]:
                        bad = True  # exceeded cap: undo and reject
                    elif used[h] == caph[h]:
                        newfull |= b
                    mm ^= b
                if not bad:
                    dfs(i + 1, remaining - p, newfull)
                mm = m
                while mm:
                    b = mm & -mm
                    used[b.bit_length() - 1] -= 1
                    mm ^= b
        S[i] = 0

    init_full = 0
    for h in range(d_h):
        if caph[h] == 0:
            init_full |= 1 << h
    dfs(0, c_target, init_full)
    return nsol


# ------------------------------------------------------------ worker glue ---

_G = {}


def _init(n, e, dlist):
    _G["n"], _G["e"] = n, e
    _G["K"] = {}
    for d in dlist:
        q = n - 1 - d
        _G["K"][d] = [KSide(adj) for adj in load_catalog("r44", q)]
    _G["H"] = {d: load_catalog("r35", d) for d in dlist}


def _run_h(task):
    d, h_idx = task
    n, e = _G["n"], _G["e"]
    H = HSide(_G["H"][d][h_idx])
    lines, counts = [], []
    for k_idx, K in enumerate(_G["K"][d]):
        sols = []
        cnt = solve_pair(H, K, d, e, sols.append)
        if not cnt:
            continue
        counts.append((d, h_idx, k_idx, cnt))
        q = K.q
        for out in sols:
            adj = [0] * n
            for h in range(d):
                adj[0] |= 1 << (1 + h)
                adj[1 + h] |= 1
            for h in range(d):
                a = H.adj[h]
                while a:
                    b = a & -a
                    adj[1 + h] |= 1 << (1 + (b.bit_length() - 1))
                    a ^= b
            for k in range(q):
                a = K.adj[k]
                while a:
                    b = a & -a
                    adj[1 + d + k] |= 1 << (1 + d + (b.bit_length() - 1))
                    a ^= b
                m = out[k]
                while m:
                    b = m & -m
                    h = b.bit_length() - 1
                    adj[1 + d + k] |= 1 << (1 + h)
                    adj[1 + h] |= 1 << (1 + d + k)
                    m ^= b
            # trusted final check (safety net)
            assert sum(popcount(a) for a in adj) // 2 == e
            assert max(popcount(a) for a in adj) == d
            assert check_ramsey_graph(n, adj, 4, 5) == (True, True), \
                f"constraint bug at d={d} h={h_idx} k={k_idx}"
            lines.append(encode_graph6(n, adj))
    return lines, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--e", type=int, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--outdir", default=os.path.join(DATA, "..", "outA"))
    args = ap.parse_args()
    n, e = args.n, args.e
    d_lo = max((2 * e + n - 1) // n, 0)
    dlist = [d for d in range(d_lo, min(13, n - 1) + 1) if n - 1 - d <= 17]
    t0 = time.time()
    hcats = {d: load_catalog("r35", d) for d in dlist}
    tasks = [(d, i) for d in dlist for i in range(len(hcats[d]))]
    os.makedirs(args.outdir, exist_ok=True)
    g6_path = os.path.join(args.outdir, f"r45_{n}_{e}.g6")
    csv_path = os.path.join(args.outdir, f"r45_{n}_{e}.counts.csv")
    total = 0
    with Pool(args.workers, initializer=_init, initargs=(n, e, dlist)) as pool, \
         open(g6_path, "w") as fg, open(csv_path, "w") as fc:
        fc.write("d,h_idx,k_idx,n_solutions\n")
        allcounts = []
        done = 0
        for lines, counts in pool.imap_unordered(_run_h, tasks, chunksize=1):
            for ln in lines:
                fg.write(ln + "\n")
            total += len(lines)
            allcounts.extend(counts)
            done += 1
            if done % 10 == 0 or done == len(tasks):
                print(f"  {done}/{len(tasks)} H-tasks, {total} solutions, "
                      f"{time.time()-t0:.0f}s", file=sys.stderr, flush=True)
        for row in sorted(allcounts):
            fc.write(",".join(map(str, row)) + "\n")
    dt = time.time() - t0
    print(f"n={n} e={e}: d in {dlist}, {len(tasks)} H-tasks, "
          f"{total} raw solutions, {dt:.1f}s -> {g6_path}")


if __name__ == "__main__":
    main()
