"""Equivalence audit v3: fast.

Part A: EXHAUSTIVE cone enumeration for all catalog (H,K) pairs with
        d*q <= 12, comparing constraint-set verdict vs direct Ramsey check.
Part B: mutation testing around real decompositions with an INCREMENTAL
        Ramsey check (a violating K4/I5 in the mutant must contain both
        endpoints of some flipped v-pair, since the base graph is Ramsey).
"""
import os
import random
import sys
from itertools import combinations

from g6lib import (load_g6_file, has_clique, has_indep, indep_sets, cliques,
                   induced)

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
random.seed(20260813)
fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL:", msg, flush=True)


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


def constraints_ok(kn, S, pre):
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
        sj = S[j]
        for q in I4m:
            if not (sj & q):
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


def main():
    # ---------------- Part A ----------------
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
                    mask = (1 << hn) - 1
                    for code in range(1 << (hn * kn)):
                        S = [(code >> (j * hn)) & mask for j in range(kn)]
                        c = constraints_ok(kn, S, pre)
                        n, adj = glue(hn, hadj, kn, kadj, S)
                        r = (not has_clique(n, adj, 4)
                             and not has_indep(n, adj, 5))
                        if c != r:
                            fail(f"A d{d}h{hi}q{q}k{ki}: constr={c} "
                                 f"ramsey={r} S={S}")
                            if len(fails) > 10:
                                sys.exit(1)
                        if c:
                            sat += 1
                        total += 1
        print(f"Part A d={d} done (cum {total}, sat {sat})", flush=True)

    # ---------------- Part B ----------------
    strata = [(13, 52), (16, 71), (17, 78), (21, 107)]
    mtotal = 0
    for n0, e0 in strata:
        graphs = load_g6_file(f"{DATA}/r45extreme/r45{n0}.{e0}.g6")
        for gi, (gn, adj0) in enumerate(graphs):
            degs = [bin(a).count("1") for a in adj0]
            Delta = max(degs)
            v = degs.index(Delta)
            Hverts = [u for u in range(gn) if (adj0[v] >> u) & 1]
            Kverts = [u for u in range(gn)
                      if u != v and not (adj0[v] >> u) & 1]
            hn, hadj = induced(adj0, Hverts)
            kn, kadj = induced(adj0, Kverts)
            hpos = {u: i for i, u in enumerate(Hverts)}
            S0 = []
            for kv in Kverts:
                m = 0
                for u in Hverts:
                    if (adj0[kv] >> u) & 1:
                        m |= 1 << hpos[u]
                S0.append(m)
            pre = precompute(hn, hadj, kn, kadj)
            if not constraints_ok(kn, list(S0), pre):
                fail(f"B real r45{n0}.{e0}#{gi}: real cone fails constraints")
            n, adjbase = glue(hn, hadj, kn, kadj, S0)
            # (glued base graph is Ramsey: necessity audit already proved the
            # published graph is; glue() rebuilds an isomorphic copy)
            others = list(range(n))
            for _ in range(600):
                S = list(S0)
                flips = set()
                for _ in range(random.randint(1, 3)):
                    j = random.randrange(kn)
                    i = random.randrange(hn)
                    S[j] ^= 1 << i
                    p = (1 + i, 1 + hn + j)   # glued labels (h, k)
                    if p in flips:
                        flips.discard(p)      # net flip cancelled
                    else:
                        flips.add(p)
                # build mutant adjacency
                adj = list(adjbase)
                for (a, b) in flips:
                    adj[a] ^= 1 << b
                    adj[b] ^= 1 << a
                # incremental Ramsey check
                viol = False
                for (a, b) in flips:
                    if viol:
                        break
                    rest = [x for x in others if x != a and x != b]
                    if (adj[a] >> b) & 1:
                        # possible new K4 through edge (a,b)
                        common = [x for x in rest
                                  if (adj[a] >> x) & 1 and (adj[b] >> x) & 1]
                        for x, y in combinations(common, 2):
                            if (adj[x] >> y) & 1:
                                viol = True
                                break
                    else:
                        # possible new I5 through non-edge (a,b)
                        noncom = [x for x in rest
                                  if not (adj[a] >> x) & 1
                                  and not (adj[b] >> x) & 1]
                        for x, y, z in combinations(noncom, 3):
                            if not ((adj[x] >> y) & 1 or (adj[x] >> z) & 1
                                    or (adj[y] >> z) & 1):
                                viol = True
                                break
                r = not viol
                c = constraints_ok(kn, S, pre)
                if c != r:
                    fail(f"B mut r45{n0}.{e0}#{gi}: constr={c} ramsey={r} "
                         f"flips={sorted(flips)}")
                    if len(fails) > 10:
                        sys.exit(1)
                mtotal += 1
        print(f"Part B r45{n0}.{e0} done (cum mutants {mtotal})", flush=True)

    print(f"\nPart A: {total} exhaustive ({sat} sat); Part B: {mtotal} "
          f"mutants. {len(fails)} mismatches", flush=True)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
