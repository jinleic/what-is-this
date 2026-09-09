"""Soundness/completeness audit of the constraint set as an *equivalence*.

For (H,K) pairs drawn from the real catalogs, enumerate cone assignments
(exhaustively when 2^(d*q) is small, randomly otherwise) and check:

    constraints(H,K,S) == (glued graph has no K4 and no I5)

for EVERY assignment. Any direction failing breaks the spec:
 - constraint holds but glue not Ramsey  -> unsound (spec emits bad graphs)
 - glue Ramsey but constraint fails     -> incomplete (spec misses graphs) FATAL
"""
import os
import random
import sys
from itertools import combinations

from g6lib import (parse_g6, has_clique, has_indep, indep_sets, cliques,
                   load_g6_file)

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
random.seed(20260813)

fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL:", msg)


def constraints_ok(hn, hadj, kn, kadj, S):
    """The spec's constraint set, verbatim."""
    Hedges = [(a, b) for a in range(hn) for b in range(a + 1, hn)
              if (hadj[a] >> b) & 1]
    # (2,2)
    for a in range(kn):
        for b in range(a + 1, kn):
            if (kadj[a] >> b) & 1:
                c = S[a] & S[b]
                if any((c >> x) & 1 and (c >> y) & 1 for x, y in Hedges):
                    return False
    # (1,3)
    for t in cliques(kn, kadj, 3):
        if S[t[0]] & S[t[1]] & S[t[2]]:
            return False
    # (4,1)
    I4masks = [sum(1 << i for i in q) for q in indep_sets(hn, hadj, 4)]
    for j in range(kn):
        for q in I4masks:
            if not (S[j] & q):
                return False
    # (3,2)
    I3masks = [sum(1 << i for i in t) for t in indep_sets(hn, hadj, 3)]
    for a in range(kn):
        for b in range(a + 1, kn):
            if not (kadj[a] >> b) & 1:
                un = S[a] | S[b]
                for t in I3masks:
                    if not (un & t):
                        return False
    # (2,3)
    I2masks = [sum(1 << i for i in p) for p in indep_sets(hn, hadj, 2)]
    for t in indep_sets(kn, kadj, 3):
        un = S[t[0]] | S[t[1]] | S[t[2]]
        for p in I2masks:
            if not (un & p):
                return False
    return True


def glue(hn, hadj, kn, kadj, S):
    """Build G: vertex 0=v, 1..hn=H, hn+1..hn+kn=K."""
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


def run_pair(tag, hn, hadj, kn, kadj, exhaustive_limit=17, rand_samples=4000):
    nbits = hn * kn
    checked = sat = 0
    if nbits <= exhaustive_limit:
        space = range(1 << nbits)
    else:
        space = (random.getrandbits(nbits) for _ in range(rand_samples))
    for code in space:
        S = [(code >> (j * hn)) & ((1 << hn) - 1) for j in range(kn)]
        c = constraints_ok(hn, hadj, kn, kadj, S)
        n, adj = glue(hn, hadj, kn, kadj, S)
        r = ramsey45(n, adj)
        checked += 1
        if c != r:
            fail(f"{tag}: constraints={c} ramsey={r} S={S}")
            if len(fails) > 5:
                return checked, sat
        if c:
            sat += 1
    return checked, sat


def main():
    r35 = {d: load_g6_file(f"{DATA}/r35_{d}.g6") for d in (3, 4, 5)}
    r44 = {q: load_g6_file(f"{DATA}/r44_{q}.g6") for q in (3, 4, 5)}

    total = totsat = 0
    # exhaustive on all (H,K) pairs with d*q <= 17
    for d in (3, 4, 5):
        for q in (3, 4, 5):
            if d * q > 17:
                continue
            for hi, (hn, hadj) in enumerate(r35[d]):
                for ki, (kn, kadj) in enumerate(r44[q]):
                    c, s = run_pair(f"d{d}h{hi}q{q}k{ki}", hn, hadj, kn, kadj)
                    total += c
                    totsat += s
        print(f"exhaustive d={d} done (cum {total} assignments, {totsat} sat)")

    # random sampling on bigger pairs (d=5,q=4/5; d=4,q=5)
    for d, q in ((4, 5), (5, 4), (5, 5)):
        for hi, (hn, hadj) in enumerate(r35[d][:10]):
            for ki, (kn, kadj) in enumerate(r44[q][:10]):
                c, s = run_pair(f"R d{d}h{hi}q{q}k{ki}", hn, hadj, kn, kadj,
                                exhaustive_limit=0, rand_samples=800)
                total += c
                totsat += s
        print(f"random d={d} q={q} done (cum {total}, sat {totsat})")

    print(f"\nchecked {total} cone assignments, {totsat} satisfied both;"
          f" {len(fails)} mismatches")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
