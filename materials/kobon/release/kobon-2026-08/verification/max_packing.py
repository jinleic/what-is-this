#!/usr/bin/env python3
"""Exact maximum interior-disjoint triangle packing of a FIXED arrangement
(broad Kobon convention: selected triangles may be crossed by other lines;
only pairwise interior-disjointness is required).

For an arrangement A of n lines this computes the independence number of the
"overlap graph" whose vertices are all nondegenerate triples of A and whose
edges join triples with intersecting open interiors. That number is exactly
the Kobon count of A in the broad convention, so it is a *lower bound
certificate* for K_gen(n) whenever it exceeds the published record.

All geometry is exact rational (campaign engine predicates); the packing
optimum is certified by SAT: a K-family is a model, and UNSAT at K+1
certifies optimality for this arrangement.

Usage:
  max_packing.py FILE.json [--from K] [--time-budget SECONDS]

Input JSON accepts three schemas:
  {"lines_frac": [[A,B,C], ...]}   (A x + B y = C)
  {"lines":      [[m,b], ...]}     (y = m x + b)
  {"ms": [...], "bs": [...]}
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402
from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import CNF, IDPool  # noqa: E402
from pysat.solvers import Solver  # noqa: E402


def load_lines(path):
    d = json.loads(Path(path).read_text())
    if "lines_frac" in d:
        L = [tuple(Fr(str(x)) for x in t) for t in d["lines_frac"]]
        if any(b == 0 for (_, b, _) in L):
            raise SystemExit("vertical line: slope-intercept engine cannot model it")
        ms = [-a / b for (a, b, _) in L]
        bs = [c / b for (_, b, c) in L]
    elif "lines" in d:
        ms = [Fr(str(p[0])) for p in d["lines"]]
        bs = [Fr(str(p[1])) for p in d["lines"]]
    elif "ms" in d and "bs" in d:
        ms = [Fr(str(x)) for x in d["ms"]]
        bs = [Fr(str(x)) for x in d["bs"]]
    else:
        raise SystemExit(f"unrecognized schema: {sorted(d)}")
    tri = d.get("triangles")
    known = d.get("count")
    if known is None:
        known = tri if isinstance(tri, int) else (len(tri) if tri else None)
    return ms, bs, known, d


def build_graph(ms, bs):
    n = len(ms)
    Xp = engine.crossings(n, ms, bs)
    cands = [t for t in combinations(range(n), 3) if engine.tri_ok(Xp, t)]
    H = {t: tuple(engine.hom(v) for v in engine.tri_verts(Xp, t)) for t in cands}
    conflicts = []
    for a, b in combinations(cands, 2):
        if not engine.interiors_disjoint_h(H[a], H[b]):
            conflicts.append((a, b))
    return cands, conflicts


def max_packing(cands, conflicts, start, budget):
    idx = {t: i + 1 for i, t in enumerate(cands)}
    base = CNF()
    for a, b in conflicts:
        base.append([-idx[a], -idx[b]])
    best, best_sel = start - 1, None
    K = start
    t0 = time.time()
    while True:
        if time.time() - t0 > budget:
            print(json.dumps({"stage": "budget-exhausted", "last_proved": best,
                              "next_untested": K}), flush=True)
            return best, best_sel, False
        cnf = CNF(from_clauses=base.clauses)
        pool = IDPool(start_from=len(cands) + 1)
        card = CardEnc.atleast(lits=list(range(1, len(cands) + 1)), bound=K,
                              vpool=pool, encoding=EncType.totalizer)
        cnf.extend(card.clauses)
        with Solver(name="cadical195", bootstrap_with=cnf.clauses) as s:
            sat = s.solve()
            if sat:
                model = set(l for l in s.get_model() if 0 < l <= len(cands))
                best, best_sel = K, [list(t) for t in cands if idx[t] in model]
                print(json.dumps({"stage": "SAT", "K": K,
                                  "family": len(best_sel)}), flush=True)
                K += 1
                continue
        print(json.dumps({"stage": "UNSAT", "K": K, "optimum": best}), flush=True)
        return best, best_sel, True


def main():
    path = sys.argv[1]
    start = None
    budget = 900.0
    args = sys.argv[2:]
    for i, a in enumerate(args):
        if a == "--from":
            start = int(args[i + 1])
        if a == "--time-budget":
            budget = float(args[i + 1])
    ms, bs, known, _ = load_lines(path)
    n = len(ms)
    t0 = time.time()
    cands, conflicts = build_graph(ms, bs)
    print(json.dumps({"stage": "graph", "n": n, "candidate_triangles": len(cands),
                      "conflict_pairs": len(conflicts), "known_count": known,
                      "seconds": round(time.time() - t0, 1)}), flush=True)
    if start is None:
        start = (known or 1)
    opt, sel, proved = max_packing(cands, conflicts, start, budget)
    out = {"file": path, "n": n, "known_count": known, "packing_optimum": opt,
           "optimality_proved": proved, "selection": sel}
    print(json.dumps({"stage": "RESULT", "n": n, "known": known,
                      "max_packing": opt, "proved_optimal": proved}), flush=True)
    outp = Path(path).with_suffix(".maxpacking.json")
    outp.write_text(json.dumps(out, indent=1))
    print(json.dumps({"stage": "written", "path": str(outp)}), flush=True)


if __name__ == "__main__":
    main()
