"""Equivalence audit v2 (feasible in pure Python).

Part A: EXHAUSTIVE over all cone assignments for all catalog (H,K) pairs
        with d*q <= 12 (d,q in {1..4}x{1..4} plus d in {1,2},q up to 6 etc.)
Part B: boundary mutation testing around REAL decompositions taken from the
        published strata: flip 1-3 random cone bits, check
        constraints(H,K,S) == Ramsey(glue) every time. This probes exactly
        the boundary where a wrong/missing constraint would show up.
"""
import random
import sys
from itertools import combinations

from g6lib import (parse_g6, has_clique, has_indep, indep_sets, cliques,
                   induced, load_g6_file, ecount)

DATA = "/Users/jinleic/jinleic-workspace/math/r55/data"
random.seed(20260813)
fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL:", msg)


def constraints_ok(hn, hadj, kn, kadj, S, pre):
    Hedges, I2m, I3m, I4m, ktri, kipair, kitri, kedge = pre
    for a, b in kedge:                      # (2,2)
        c = S[a] & S[b]
        for em in Hedges:
            if (c & em) == em:
                return False
    for a, b, c in ktri:                    # (1,3)
        if S[a] & S[b] & S[c]:
            return False
    for j in range(kn):                     # (4,1)
        for q in I4m:
            if not (S[j] & q):
                return False
    for a, b in kipair:                     # (3,2)
        un = S[a] | S[b]
        for t in I3m:
            if not (un & t):
                return False
    for a, b, c in kitri:                   # (2,3)
        un = S[a] | S[b] | S[c]
        for p in I2m:
            if not (un & p):
                return False
    return True


def precompute(hn, hadj, kn, kadj):
    Hedges = [(1 << a) | (1 << b) for a in range(hn) for b in range(a + 1, hn)
              if (hadj[a] >> b) & 1]
    I2m = [sum(1 << i for i in p) for p in indep_sets(hn, hadj, 2)]
    I3m = [sum(1 << i for i in t) for t in indep_sets(hn, hadj, 3)]
    I4m = [sum(1 << i for i in q) for q in indep_sets(hn, hadj, 4)]
    ktri = list(cliques(kn, kadj, 3))
    kipair = [(a, b) for a in range(kn) for b in range(a + 1, kn)
              if not (kadj[a] >> b) & 1]
    kitri = list(indep_sets(kn, kadj, 3))
    kedge = [(a, b) for a in range(kn) for b in range(a + 1, kn)
             if (kadj[a] >> b) & 1]
    return Hedges, I2m, I3m, I4m, ktri, kipair, kitri, kedge


def glue(hn, hadj, kn, kadj, S):
    n = 1 + hn + kn
    adj = [0] * n
    for i in range(hn):
        adj[0] |= 1 << (1 + i)
        adj[1 + i] |= 1
    for a in range(hn):
        for b in range(a + 1, hn):
            if (hadj[a] >> b) & 1:
                adj[1 + a] |= 1 << (1 + b)
                adj[1 + b] |= 1 << (1 + a)
    for a in range(kn):
        for b in range(a + 1, kn):
            if (kadj[a] >> b) & 1:
                adj[1 + hn + a] |= 1 << (1 + hn + b)
                adj[1 + hn + b] |= 1 << (1 + hn + a)
    for j in range(kn):
        for i in range(hn):
            if (S[j] >> i) & 1:
                adj[1 + hn + j] |= 1 << (1 + i)
                adj[1 + i] |= 1 << (1 + hn + j)
    return n, adj


def ramsey45(n, adj):
    return not has_clique(n, adj, 4) and not has_indep(n, adj, 5)


def check_one(tag, hn, hadj, kn, kadj, S, pre):
    c = constraints_ok(hn, hadj, kn, kadj, S, pre)
    n, adj = glue(hn, hadj, kn, kadj, S)
    r = ramsey45(n, adj)
    if c != r:
        fail(f"{tag}: constraints={c} but ramsey={r} S={S}")
    return c


def main():
    # ---- Part A: exhaustive small pairs ----
    total = sat = 0
    for d in (1, 2, 3, 4):
        r35 = load_g6_file(f"{DATA}/r35_{d}.g6")
        for q in (1, 2, 3, 4):
            if d * q > 12:
                continue
            r44 = load_g6_file(f"{DATA}/r44_{q}.g6")
            for hi, (hn, hadj) in enumerate(r35):
                for ki, (kn, kadj) in enumerate(r44):
                    pre = precompute(hn, hadj, kn, kadj)
                    for code in range(1 << (hn * kn)):
                        S = [(code >> (j * hn)) & ((1 << hn) - 1)
                             for j in range(kn)]
                        if check_one(f"d{d}h{hi}q{q}k{ki}", hn, hadj, kn,
                                     kadj, S, pre):
                            sat += 1
                        total += 1
                        if len(fails) > 10:
                            print("too many failures, aborting")
                            sys.exit(1)
        print(f"exhaustive d={d} done (cum {total}, sat {sat})")

    # ---- Part B: mutation around real decompositions ----
    strata = [(13, 52), (16, 71), (17, 78), (21, 107)]
    mtotal = 0
    for n, e in strata:
        graphs = load_g6_file(f"{DATA}/r45extreme/r45{n}.{e}.g6")
        for gi, (gn, adj) in enumerate(graphs):
            degs = [bin(a).count("1") for a in adj]
            Delta = max(degs)
            v = degs.index(Delta)
            Hverts = [u for u in range(gn) if (adj[v] >> u) & 1]
            Kverts = [u for u in range(gn) if u != v and not (adj[v] >> u) & 1]
            hn, hadj = induced(adj, Hverts)
            kn, kadj = induced(adj, Kverts)
            hpos = {u: i for i, u in enumerate(Hverts)}
            S0 = []
            for kv in Kverts:
                m = 0
                for u in Hverts:
                    if (adj[kv] >> u) & 1:
                        m |= 1 << hpos[u]
                S0.append(m)
            pre = precompute(hn, hadj, kn, kadj)
            # sanity: the real cone must satisfy both sides
            if not check_one(f"real r45{n}.{e}#{gi}", hn, hadj, kn, kadj,
                             list(S0), pre):
                fail(f"real cone of r45{n}.{e}#{gi} does not satisfy "
                     f"constraints")
            # mutations: flip 1..3 random bits, 400 mutants per graph
            for _ in range(400):
                S = list(S0)
                for _ in range(random.randint(1, 3)):
                    j = random.randrange(kn)
                    i = random.randrange(hn)
                    S[j] ^= 1 << i
                check_one(f"mut r45{n}.{e}#{gi}", hn, hadj, kn, kadj, S, pre)
                mtotal += 1
                if len(fails) > 10:
                    print("too many failures, aborting")
                    sys.exit(1)
        print(f"mutations r45{n}.{e} done (cum {mtotal})")

    print(f"\nPart A: {total} exhaustive assignments ({sat} sat); "
          f"Part B: {mtotal} mutants. {len(fails)} mismatches")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
