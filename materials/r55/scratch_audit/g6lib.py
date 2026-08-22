"""Minimal independent graph6 parser + bitmask graph utilities.

Written from the graph6 format description from scratch for the audit;
deliberately does NOT import anything from src/.
"""
from itertools import combinations


def parse_g6(line):
    """Return (n, adj) where adj[i] is an int bitmask of neighbors of i."""
    s = line.strip()
    if s.startswith(">>graph6<<"):
        s = s[10:]
    data = [ord(c) - 63 for c in s]
    assert all(0 <= x <= 63 for x in data), "bad g6 char"
    if data[0] <= 62:
        n = data[0]
        bits_start = 1
    else:
        # 126 -> long form (not needed for n<=63 but implement 3-byte form)
        assert data[1] <= 62
        n = (data[1] << 12) | (data[2] << 6) | data[3]
        bits_start = 4
    # bit stream, column-major over upper triangle (j from 1..n-1, i from 0..j-1)
    bits = []
    for x in data[bits_start:]:
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    need = n * (n - 1) // 2
    assert len(bits) >= need, "g6 too short"
    adj = [0] * n
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            idx += 1
    return n, adj


def edges(n, adj):
    return [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i] >> j & 1]


def ecount(n, adj):
    return sum(bin(a).count("1") for a in adj) // 2


def has_clique(n, adj, k):
    """True if graph contains a k-clique."""
    for combo in combinations(range(n), k):
        ok = True
        for a, b in combinations(combo, 2):
            if not (adj[a] >> b) & 1:
                ok = False
                break
        if ok:
            return True
    return False


def has_indep(n, adj, k):
    for combo in combinations(range(n), k):
        ok = True
        for a, b in combinations(combo, 2):
            if (adj[a] >> b) & 1:
                ok = False
                break
        if ok:
            return True
    return False


def indep_sets(n, adj, k):
    """Yield all independent k-sets (as tuples)."""
    for combo in combinations(range(n), k):
        ok = True
        for a, b in combinations(combo, 2):
            if (adj[a] >> b) & 1:
                ok = False
                break
        if ok:
            yield combo


def cliques(n, adj, k):
    for combo in combinations(range(n), k):
        ok = True
        for a, b in combinations(combo, 2):
            if not (adj[a] >> b) & 1:
                ok = False
                break
        if ok:
            yield combo


def induced(adj, verts):
    """Induced subgraph on list verts; returns (m, adj2) with relabeling by position."""
    m = len(verts)
    pos = {v: i for i, v in enumerate(verts)}
    adj2 = [0] * m
    for i, v in enumerate(verts):
        for j, w in enumerate(verts):
            if i < j and (adj[v] >> w) & 1:
                adj2[i] |= 1 << j
                adj2[j] |= 1 << i
    return m, adj2


def load_g6_file(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(parse_g6(line))
    return out
