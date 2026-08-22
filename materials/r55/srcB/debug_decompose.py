#!/usr/bin/env python3
"""Decompose the published (12,48) graph at a max-degree vertex and test
whether the constraint frame accepts the true (H, K, S_k) triple."""
import itertools, sys

def parse_g6(line):
    n = ord(line[0]) - 63
    adj = [0] * n
    bits = []
    for ch in line[1:]:
        c = ord(ch) - 63
        bits.extend((c >> (5 - t)) & 1 for t in range(6))
    p = 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            p += 1
    return n, adj

line = open(sys.argv[1]).readline().strip()
n, adj = parse_g6(line)
e = sum(bin(a).count("1") for a in adj) // 2
degs = [bin(a).count("1") for a in adj]
print("n", n, "e", e, "degs", sorted(degs))

v = max(range(n), key=lambda x: degs[x])
d = degs[v]
Hv = [u for u in range(n) if adj[v] >> u & 1]
Kv = [u for u in range(n) if u != v and not (adj[v] >> u & 1)]
q = len(Kv)
print("v", v, "d", d, "q", q)

def sub(verts):
    m = len(verts)
    a = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if adj[verts[i]] >> verts[j] & 1:
                a[i] |= 1 << j
                a[j] |= 1 << i
    return a

aH = sub(Hv)
aK = sub(Kv)
eH = sum(bin(x).count("1") for x in aH) // 2
eK = sum(bin(x).count("1") for x in aK) // 2
S = []
for k in range(q):
    s = 0
    for i in range(d):
        if adj[Kv[k]] >> Hv[i] & 1:
            s |= 1 << i
    S.append(s)
cone = sum(bin(s).count("1") for s in S)
print("eH", eH, "eK", eK, "cone", cone, "= e-d-eH-eK?", e - d - eH - eK)

def indep_sets(a, m, size):
    return [c for c in itertools.combinations(range(m), size)
            if all(not (a[x] >> y & 1) for x, y in itertools.combinations(c, 2))]

I2 = indep_sets(aH, d, 2); I3 = indep_sets(aH, d, 3); I4 = indep_sets(aH, d, 4)
print("H: #I2", len(I2), "#I3", len(I3), "#I4", len(I4),
      "triangles", len([c for c in itertools.combinations(range(d), 3)
                        if all(aH[x] >> y & 1 for x, y in itertools.combinations(c, 2))]))

def hits(s, fam):
    return all(any(s >> x & 1 for x in c) for c in fam)

ok = True
for k in range(q):
    if not hits(S[k], I4):
        print("FAIL (4,1) at k", k); ok = False
    if bin(S[k]).count("1") + bin(aK[k]).count("1") > d:
        print("FAIL degcap K at", k); ok = False
for h in range(d):
    dg = 1 + bin(aH[h]).count("1") + sum(S[k] >> h & 1 for k in range(q))
    if dg > d:
        print("FAIL degcap H at", h, dg); ok = False
for k1, k2 in itertools.combinations(range(q), 2):
    if aK[k1] >> k2 & 1:
        inter = S[k1] & S[k2]
        bad = any((inter >> x & 1) and (inter >> y & 1) for x, y in
                  ((x, y) for x in range(d) for y in range(x + 1, d) if aH[x] >> y & 1))
        if bad:
            print("FAIL (2,2) at", k1, k2); ok = False
    else:
        if not hits(S[k1] | S[k2], I3):
            print("FAIL (3,2) at", k1, k2); ok = False
for k1, k2, k3 in itertools.combinations(range(q), 3):
    a12, a13, a23 = aK[k1] >> k2 & 1, aK[k1] >> k3 & 1, aK[k2] >> k3 & 1
    if a12 and a13 and a23:
        if S[k1] & S[k2] & S[k3]:
            print("FAIL (1,3) at", k1, k2, k3); ok = False
    if not a12 and not a13 and not a23:
        if not hits(S[k1] | S[k2] | S[k3], I2):
            print("FAIL (2,3) at", k1, k2, k3); ok = False
print("frame constraints all pass:", ok)
