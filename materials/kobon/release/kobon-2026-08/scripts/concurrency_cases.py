#!/usr/bin/env python3
"""n=14 / target-54 degeneracy signatures and concurrency cube builder.

Theorems used (math/kobon/report.md + CrossingBudgetMath polished transcript;
the multipoint extension is under independent audit, so every UNSAT result is
"discovery UNSAT"):

  segment/capacity budget (crossing-refined multipoint capacity theorem):
      3T + C  <=  n(n-2) - 2Q - sum_p k_p(k_p-4)
      =>  C + 2Q - 3*N3 + 0*N4 + 5*N5 + 12*N6 + 21*N7 <= 6      (B1)

  face budget (audited bounded-face penalty, engine.add_face_bound):
      T + C + Q + sum_p C(k_p-1,2) <= C(n-1,2)
      =>  Q + N3 + 3*N4 + 6*N5 + 10*N6 + 15*N7 <= 24            (B2)

  pair account (elementary):
      sum_p C(k_p,2) <= C(n,2) - Q
      =>  3*N3 + 6*N4 + 10*N5 + 15*N6 + 21*N7 <= 91 - Q         (B3)

  exact distinct-finite-point count:
      pts = C(n,2) - Q - sum_p (C(k_p,2)-1)
          = C(n,2) - Q - (2*N3 + 5*N4 + 9*N5 + 14*N6 + 20*N7)
  (NOTE: the assignment's formula "C(n,2)-Q-sum_p C(k_p-1,2)" is the bounded
  FACE count F_b = 1-n+sum_p(k_p-1) of the arrangement, not the point count;
  the point-count reduction per k-fold point is C(k,2)-1 = C(k-1,2)+(k-2).)

  k >= 8 is infeasible: a k=8 point costs face budget C(7,2)=21 and segment
  budget k(k-4)=32, hence needs 3*N3 >= 32-6 -> N3 >= 9, giving face total
  >= 21+9 = 30 > 24; k >= 9 exceeds face slack alone (C(k-1,2) >= 28 > 24).

Cube covering argument (relabeling/nesting):
  A cube fixes an exact parallel-class pattern (P-units for every pair) plus a
  positive C-unit for ONE distinguished concurrent triple in a given orbit of
  the pattern, and leaves all other concurrency FREE.  Every arrangement with
  >= 1 multipoint has a concurrent triple (any k >= 4 point yields C-true
  triples), so the union over parallel patterns x orbits covers every
  concurrency signature.  Richer signatures nest into narrower cubes (dropping
  negative C-units is monotone), so no negative C-units are needed for
  soundness/coverage.

CLI:
  enum           print the complete admissible signature lattice + tables
  orbits         print parallel patterns and triple-orbit cube recipes
  build NAME     build the cube CNF (writes scratch/kobon/n14/NAME.cnf)
"""

from __future__ import annotations

import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

import os
N = 14
TARGET = int(os.environ.get("N14_TARGET", "54"))
SEG_SLACK = N * (N - 2) - 3 * TARGET          # T=54: 6; T=55: 3
FACE_SLACK = comb(N - 1, 2) - TARGET          # T=54: 24; T=55: 23

# ----------------------------- exact arithmetic -----------------------------

def seg_penalty(k):   # k_p(k_p-4) per k-fold point
    return k * (k - 4)


def face_penalty(k):
    return comb(k - 1, 2)


def pair_cost(k):
    return comb(k, 2)


def merge_cost(k):    # distinct finite points lost per k-fold point
    return comb(k, 2) - 1


def admissible(Q, ns):
    """ns: multiplicity->count for k>=3.  Returns (ok, info)."""
    seg = 2 * Q + sum(c * seg_penalty(k) for k, c in ns.items())
    cmax = SEG_SLACK - seg
    face = Q + sum(c * face_penalty(k) for k, c in ns.items())
    pairs = sum(c * pair_cost(k) for k, c in ns.items())
    pts = comb(N, 2) - Q - sum(c * merge_cost(k) for k, c in ns.items())
    info = dict(C_max=cmax, face=face, face_slack=FACE_SLACK - face,
                pairs=pairs, pairs_slack=comb(N, 2) - Q - pairs, pts=pts)
    ok = (cmax >= 0 and face <= FACE_SLACK and pairs <= comb(N, 2) - Q
          and pts >= 1)
    return ok, info


# --------------------------- parallel class patterns -------------------------

def parallel_patterns(Q, n=N):
    """Multisets (descending tuples) of class sizes s>=2, sum C(s_i,2)=Q,
    total lines <= n."""
    if Q == 0:
        return [()]
    out = []

    def rec(q_left, max_s, sizes, lines_used):
        if q_left == 0:
            if sizes:
                out.append(tuple(sizes))
            return
        for s in range(min(max_s, n - lines_used), 1, -1):
            c = comb(s, 2)
            if 0 < c <= q_left:
                rec(q_left - c, s, sizes + [s], lines_used + s)

    rec(Q, n, [], 0)
    return sorted(set(out), key=lambda t: (-len(t), t))


def triple_orbits(pattern, n=N):
    """Orbits of one distinguished concurrent 3-set: <= 1 line per class.

    Orbit key = sorted multiset of the chosen class sizes (one line each);
    remaining 3-r lines are singletons (need enough lines)."""
    lines_used = sum(pattern)
    single_avail = n - lines_used
    orbits = []
    for r in range(0, min(3, len(pattern)) + 1):
        if r < 3 and single_avail < 3 - r:
            continue
        chosen = set()
        sizes = sorted(set(pattern))
        counts = {}

        def rec(i, cur):
            if len(cur) == r:
                chosen.add(tuple(sorted(cur)))
                return
            if i == len(sizes):
                return
            s = sizes[i]
            rec(i + 1, cur)
            if counts.get(s, 0) < pattern.count(s):
                counts[s] = counts.get(s, 0) + 1
                cur.append(s)
                rec(i, cur)
                cur.pop()
                counts[s] -= 1

        rec(0, [])
        orbits.extend(sorted(chosen))
    return orbits


def orbit_units(pattern, chosen_sizes, n=N):
    """Labels: classes first in pattern order, then singleton lines.  Triple =
    one line from each chosen class + singletons.  Returns (classes, triple)."""
    classes = []
    nxt = 0
    for s in pattern:
        classes.append(list(range(nxt, nxt + s)))
        nxt += s
    singletons = list(range(nxt, n))
    triple = []
    used = [False] * len(classes)
    for s in sorted(chosen_sizes):
        for i, cls in enumerate(classes):
            if not used[i] and len(cls) == s:
                used[i] = True
                triple.append(cls[0])
                break
        else:
            raise ValueError("orbit not realizable in pattern")
    for x in singletons:
        if len(triple) < 3:
            triple.append(x)
    assert len(triple) == 3
    return classes, tuple(sorted(triple))


# ------------------------------- enumeration --------------------------------

def enumerate_all():
    """Complete admissible (Q, {N3..N7}) lattice.  Returns list of
    (Q, ns_dict, info)."""
    rows = []
    for Q in range(0, FACE_SLACK + 1):        # face: Q <= 24
        for N3 in range(0, FACE_SLACK - Q + 1):
            for N4 in range(0, (FACE_SLACK - Q - N3) // 3 + 1):
                for N5 in range(0, 5):
                    for N6 in range(0, 3):
                        for N7 in range(0, 2):
                            ns = {k: v for k, v in
                                  ((3, N3), (4, N4), (5, N5), (6, N6),
                                   (7, N7)) if v}
                            ok, info = admissible(Q, ns)
                            if ok:
                                rows.append((Q, ns, info))
    return rows


def summarize(rows):
    print(f"n={N} T={TARGET}: segment slack {SEG_SLACK}, face slack "
          f"{FACE_SLACK}")
    print(f"TOTAL admissible signatures: {len(rows)}")
    by_bulk = {}
    for Q, ns, info in rows:
        bulk = tuple(sorted((k, c) for k, c in ns.items() if k != 3))
        by_bulk.setdefault(bulk, []).append((Q, ns.get(3, 0), info))
    for bulk, lst in sorted(by_bulk.items()):
        qs = sorted(set(Q for Q, _, _ in lst))
        n3s = sorted(set(n3 for _, n3, _ in lst))
        print(f"bulk {bulk or '()'}: {len(lst)} sigs; "
              f"Q in {qs[0]}..{qs[-1]}; N3 in {n3s[0]}..{n3s[-1]}")
    return by_bulk


# ------------------------------- cube builder --------------------------------

CUBES = {
    # name: (parallel classes, [concurrency triple units])
    # classes=None marks the true all-degeneracy monolith (no unit clauses).
    "n14monolith": (None, []),
    "n14c-q0tp1": ([], [(0, 1, 2)]),
    "n14c-q1tp1d": ([[0, 1]], [(2, 3, 4)]),
    "n14c-q1tp1s": ([[0, 1]], [(0, 2, 3)]),
    "n14c-q0tp2d": ([], [(0, 1, 2), (3, 4, 5)]),
    "n14c-q0tp2s": ([], [(0, 1, 2), (0, 3, 4)]),
    "n14c-q0k4": ([], [(0, 1, 2), (0, 1, 3)]),
}


def capacity_list(classes, triples, n=N):
    """cap_r = n-2-q_r-sum_{p on r}(k_p-4) for the designated signature
    (exactly the planted multipoints; triples sharing a pair cluster into a
    higher-multiplicity point)."""
    q = [0] * n
    for cls in classes:
        for r in cls:
            q[r] = len(cls) - 1
    pts = []
    for t in triples:
        for p in pts:
            if len(set(p) & set(t)) >= 2:
                for x in t:
                    if x not in p:
                        p.append(x)
                break
        else:
            pts.append(list(t))
    caps = []
    for r in range(n):
        c = n - 2 - q[r]
        for p in pts:
            if r in p:
                c -= (len(p) - 4)
        caps.append(c)
    return pts, caps


def build_cube(name):
    sys.path.insert(0, str(ROOT / "math" / "kobon"))
    import engine

    if name not in CUBES:
        raise SystemExit(f"unknown cube {name}; choices: {sorted(CUBES)}")
    classes, triples = CUBES[name]
    started = time.time()
    cnf, pool = engine.build_model(N, TARGET)
    print(f"build_model {time.time()-started:.1f}s vars={cnf.nv} "
          f"clauses={len(cnf.clauses)}", flush=True)
    t0 = time.time()
    crossing_lits = engine.add_triangle_crossing_indicators(cnf, pool, N)
    engine.add_face_bound(cnf, pool, N, TARGET,
                          crossing_lits=crossing_lits, per_line=True)
    engine.add_exact_selection(cnf, pool, N, TARGET)
    print(f"constraints {time.time()-t0:.1f}s vars={cnf.nv} "
          f"clauses={len(cnf.clauses)}", flush=True)

    P = lambda i, j: pool.id(("P",) + tuple(sorted((i, j))))
    Cv = lambda t: pool.id(("C",) + tuple(sorted(t)))
    pre_units = len(cnf.clauses)
    if classes is None:
        # true all-degeneracy monolith: leave every P(i,j) and C(t) FREE;
        # the regression below asserts zero unit clauses were emitted.
        par_pairs = set()
    else:
        par_pairs = {tuple(sorted(p)) for cls in classes for p in
                     combinations(cls, 2)}
        for i, j in combinations(range(N), 2):
            cnf.append([-P(i, j) if (i, j) in par_pairs else P(i, j)])
        for t in triples:
            cnf.append([Cv(t)])
    if classes is None:
        assert len(cnf.clauses) == pre_units, (
            f"monolith regression: {len(cnf.clauses)-pre_units} unit clauses "
            "emitted for classes=None")

    out = OUT / f"{name}.cnf"
    cnf.to_file(out)
    print(f"cube={name} vars={cnf.nv} clauses={len(cnf.clauses)} "
          f"bytes={out.stat().st_size} total={time.time()-started:.1f}s",
          flush=True)
    if classes is not None:
        pts, caps = capacity_list(classes, triples)
        print(f"classes={classes} triples={triples} multipoints={pts}")
        print(f"caps={caps}")
    else:
        print("monolith: all P/C free (zero unit clauses)")


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    cmd = sys.argv[1]
    if cmd == "enum":
        summarize(enumerate_all())
    elif cmd == "orbits":
        for Q in range(0, 7):
            pats = parallel_patterns(Q)
            print(f"Q={Q}: {len(pats)} parallel patterns")
            for pat in pats:
                for orb in triple_orbits(pat):
                    print("   pat", pat, "orbit", orb, "->",
                          orbit_units(pat, orb))
    elif cmd == "build":
        build_cube(sys.argv[2])
    else:
        raise SystemExit(f"bad command {cmd}")


if __name__ == "__main__":
    main()
