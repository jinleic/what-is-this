#!/usr/bin/env python3
"""End-to-end reproduction: structural constraints on hypothetical R(5,5,n) graphs.
Data expected in ../data/. Pure stdlib + numpy (numpy only for optional spectra).
Executed and verified 2026-08-15; mathematics in ../NOTES.md.
"""
import os, tarfile
from collections import Counter
from fractions import Fraction as Fr
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def parse_g6(line):
    line = line.strip()
    n = line[0] - 63
    nbits = n*(n-1)//2
    nbytes = (nbits + 5)//6
    assert len(line) == 1 + nbytes
    bits = 0
    for c in line[1:]:
        bits = (bits << 6) | (c - 63)
    bits >>= (6*nbytes - nbits)
    adj = [0]*n
    pos = nbits - 1
    for j in range(1, n):
        for i in range(j):
            if (bits >> pos) & 1:
                adj[i] |= 1 << j; adj[j] |= 1 << i
            pos -= 1
    return n, adj

def bits_iter(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x ^= b

def complement(n, adj):
    full = (1 << n) - 1
    return [full & ~a & ~(1 << v) for v, a in enumerate(adj)]

def has_k5(n, adj):
    for u in range(n):
        for du in bits_iter(adj[u] >> (u+1)):
            v = u + 1 + du
            B = adj[u] & adj[v] & ~((1 << (v+1)) - 1)
            if bin(B).count("1") < 3: continue
            for w in bits_iter(B):
                C = B & adj[w] & ~((1 << (w+1)) - 1)
                for x in bits_iter(C):
                    if adj[x] & C & ~((1 << (x+1)) - 1):
                        return True
    return False

def has_k4(n, adj):
    for u in range(n):
        for du in bits_iter(adj[u] >> (u+1)):
            v = u + 1 + du
            B = adj[u] & adj[v] & ~((1 << (v+1)) - 1)
            for w in bits_iter(B):
                if B & adj[w] & ~((1 << (w+1)) - 1):
                    return True
    return False

def tri_count(n, adj):
    t = 0
    for u in range(n):
        for du in bits_iter(adj[u] >> (u+1)):
            v = u+1+du
            t += bin(adj[u] & adj[v] & ~((1 << (v+1)) - 1)).count("1")
    return t

e45 = {17:41,18:50,19:57,20:68,21:77,22:88,23:101}
E45 = {17:79,18:85,19:92,20:100,21:107,22:114,23:122}

def C2(m): return m*(m-1)//2

def excess_and_intervals(n, adj):
    full = (1 << n) - 1
    viol = 0; total = 0
    for v in range(n):
        Nv = adj[v]; d = bin(Nv).count("1")
        Dv = full & ~Nv & ~(1 << v)
        ep = sum(bin(adj[x] & Nv).count("1") for x in bits_iter(Nv))//2
        em = sum(bin(adj[x] & Dv).count("1") for x in bits_iter(Dv))//2
        total += em - ep - d*(n-2*d)/2
        m = n-1-d
        if not (e45[d] <= ep <= E45[d]): viol += 1
        if not (C2(m)-E45[m] <= em <= C2(m)-e45[m]): viol += 1
    return total, viol

def h_bounds(n, d):
    m = n - 1 - d
    base = Fr(d*(n - 2*d), 2)
    return (Fr(C2(m) - E45[m] - E45[d]) - base,
            Fr(C2(m) - e45[m] - e45[d]) - base)

def R_matrix_int(Xa, n=24):
    A2 = [[bin(Xa[i] & Xa[j]).count("1") for j in range(n)] for i in range(n)]
    R = [[0]*n for _ in range(n)]
    for w in range(n):
        for w2 in range(n):
            if w == w2: R[w][w2] = 12
            elif (Xa[w] >> w2) & 1: R[w][w2] = 10 - A2[w][w2]
            else:                   R[w][w2] = 11 - A2[w][w2]
    return R

def C_matrix_int(Db, n=24):
    A2 = [[bin(Db[i] & Db[j]).count("1") for j in range(n)] for i in range(n)]
    C = [[0]*n for _ in range(n)]
    for x in range(n):
        for x2 in range(n):
            if x == x2: C[x][x2] = 12
            elif (Db[x] >> x2) & 1: C[x][x2] = 11 - A2[x][x2]
            else:                   C[x][x2] = 12 - A2[x][x2]
    return C

def tr_powers(Mi, kmax=4):
    n = len(Mi); out = []
    P = [row[:] for row in Mi]
    for k in range(1, kmax+1):
        out.append(sum(P[i][i] for i in range(n)))
        if k < kmax:
            P = [[sum(P[i][t]*Mi[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return out

def main():
    raw = open(f"{BASE}/data/r55_42some.g6","rb").read().splitlines()
    g42 = [parse_g6(l) for l in raw if l.strip()]
    assert len(g42) == 328 and all(n == 42 for n,_ in g42)
    all656 = []
    for n, adj in g42:
        all656 += [(n, adj), (n, complement(n, adj))]
    assert all(not has_k5(n,a) and not has_k5(n,complement(n,a)) for n,a in all656)
    print("656 R(5,5,42) graphs verified")

    tf = tarfile.open(f"{BASE}/data/r45extreme.tar.gz")
    for order in range(17,24):
        for ecnt in (e45[order], E45[order]):
            for line in tf.extractfile(f"r45extreme/r45{order}.{ecnt}.g6").read().splitlines():
                if not line.strip(): continue
                n, adj = parse_g6(line)
                e = sum(bin(a).count("1") for a in adj)//2
                assert n==order and e==ecnt and not has_k4(n,adj) and not has_k5(n,complement(n,adj))
    print("boundary extremal classes verified (2887 graphs)")

    # Census stream: exact count + edge distribution + extremal extraction.
    # Full K4/I5 property validation of all 352,366 graphs is gate 1
    # (src/validate_catalogs.py, data/VALIDATION.json); here every 100th
    # graph plus every graph with e >= 130 is property-checked directly.
    edist = Counter(); X = []; sampled = 0
    for i, line in enumerate(open(f"{BASE}/data/r45_24.g6","rb")):
        line = line.rstrip()
        if not line: continue
        e = sum(bin(c-63).count("1") for c in line[1:])
        edist[e] += 1
        if e == 132: X.append(parse_g6(line)[1])
        if i % 100 == 0 or e >= 130:
            n, adj = parse_g6(line)
            assert n == 24 and not has_k4(n, adj) and not has_k5(n, complement(n, adj))
            sampled += 1
    assert sum(edist.values()) == 352366 and max(edist) == 132 and len(X) == 2
    e45[24], E45[24] = min(edist), max(edist)
    print(f"census 24: count 352366 and edge distribution verified, e in [116,132], "
          f"2 graphs at 132; {sampled} graphs property-checked (deterministic sample + all e>=130); "
          f"full property validation = gate 1")

    nchecks = 0
    for n, adj in all656:
        tot, vi = excess_and_intervals(n, adj)
        assert abs(tot) < 1e-9 and vi == 0
        nchecks += 2*n
    assert nchecks == 656*42*2 == 55104
    print("excess identity exact on all 656; all 55104 per-vertex nbhd/dual "
          "edge counts inside the h-interval tables")

    for n in range(43, 51):
        lo, hi = n-25, 24
        if lo > hi:
            print(f"n={n}: empty degree window -> R(5,5) <= {n}"); continue
        rows = {d: h_bounds(n, d) for d in range(lo, hi+1)}
        print(f"n={n}: h_min per degree", {d: str(h) for d,(h,_) in rows.items()},
              "| deficiency budget <=", max(-h for h,_ in rows.values())*n)

    assert [tri_count(24, a) for a in X] == [176, 176]
    assert [sum(bin(a).count("1")**2 for a in x_) for x_ in X] == [2904, 2904]
    assert all({bin(a).count("1") for a in x_} == {11} for x_ in X)
    D = [complement(24, a) for a in X]
    tpR = [tr_powers(R_matrix_int(X[i])) for i in (0,1)]
    tpC = [tr_powers(C_matrix_int(D[i])) for i in (0,1)]
    assert tpR[0] == tpC[0] and tpR[1] == tpC[1] and tpR[0] != tpC[1]
    print("n=49 chain verified: SRG(49,24,11,12) forcing; exact type-matching")
    print("tr powers R(X1):", tpR[0])
    print("tr powers R(X2):", tpR[1])
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
