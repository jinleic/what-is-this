#!/usr/bin/env python3
"""Canonical Q=1, no-concurrency triangular-face SAT encoding.

For a simple affine pseudoline arrangement with one parallel pair, a pairwise
crossing triple bounds a triangular face iff its two vertices on each of the
three supporting lines are consecutive in that line's crossing order.  This
script adds those direct gap clauses to the audited A1--A7 order relaxation in
``engine.build_model(..., rules=())``.  No sidedness or pairwise-disjointness
variables are needed: distinct arrangement faces have disjoint interiors.

The unique parallel class is fixed to {0,1} WLOG by choosing the cut in the
circular slope order immediately before that class.  Every finite triple is
fixed nonconcurrent.  UNSAT is sound for the complete Q=1/no-concurrency branch;
SAT remains a combinatorial candidate until straight-line realization.

Usage: q1_gap_faces.py N T [--out FILE] [--solve] [--solver cadical195]
                           [--symbreak]
"""
from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def build(n, target, symbreak=False):
    cnf, pool = engine.build_model(n, target, rules=())
    pair = lambda i, j: pool.id(("P",) + tuple(sorted((i, j))))
    concurrent = lambda t: pool.id(("C",) + tuple(sorted(t)))
    order = lambda r, i, j: pool.id(("X", r, i, j))
    selected = lambda t: pool.id(("S",) + tuple(sorted(t)))

    for i, j in combinations(range(n), 2):
        cnf.append([-pair(i, j) if (i, j) == (0, 1) else pair(i, j)])
    for triple in combinations(range(n), 3):
        cnf.append([-concurrent(triple)])
    if symbreak:
        # P(0,1)=false makes engine.add_symbreak's guarded (0,1,2) clause
        # tautological.  Triple (0,2,3) is fixed crossing/nonconcurrent here,
        # and the global b -> -b involution flips this order bit.
        cnf.append([order(0, 2, 3)])

    # If t is selected, no outside crossing may lie strictly between its two
    # vertices on any side line.  The two clauses cover both endpoint orders.
    for triple in combinations(range(n), 3):
        s = selected(triple)
        for side in triple:
            a, b = (x for x in triple if x != side)
            for outside in range(n):
                if outside in triple:
                    continue
                cnf.append([
                    -s, -order(side, a, outside),
                    -order(side, outside, b),
                ])
                cnf.append([
                    -s, -order(side, b, outside),
                    -order(side, outside, a),
                ])

    q_min = [1, 1] + [0] * (n - 2)
    engine.add_no_concurrency_capacities(
        cnf, pool, n, q_min, target=target)
    engine.add_exact_selection(cnf, pool, n, target)
    return cnf, pool


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    n, target = int(sys.argv[1]), int(sys.argv[2])
    output = None
    solver = "cadical195"
    solve = "--solve" in sys.argv
    symbreak = "--symbreak" in sys.argv
    for i, arg in enumerate(sys.argv):
        if arg == "--out":
            output = Path(sys.argv[i + 1])
        elif arg == "--solver":
            solver = sys.argv[i + 1]

    started = time.time()
    cnf, pool = build(n, target, symbreak=symbreak)
    print(json.dumps({
        "stage": "built",
        "formulation": "q1-gap-faces-symbreak" if symbreak else "q1-gap-faces",
        "n": n,
        "target": target,
        "vars": cnf.nv,
        "clauses": len(cnf.clauses),
        "seconds": round(time.time() - started, 3),
        "rss_gb": round(rss_gb(), 3),
    }), flush=True)
    if output is not None:
        cnf.to_file(str(output))
        print(json.dumps({
            "stage": "written", "path": str(output),
            "bytes": output.stat().st_size,
        }), flush=True)
    if not solve:
        return 0

    from pysat.solvers import Solver
    solve_started = time.time()
    with Solver(name=solver, bootstrap_with=cnf.clauses) as sat:
        satisfiable = sat.solve()
        stats = sat.accum_stats()
        model = sat.get_model() if satisfiable else None
    result = {
        "stage": "VERDICT",
        "n": n,
        "target": target,
        "result": "SAT" if satisfiable else "UNSAT",
        "seconds": round(time.time() - solve_started, 3),
        "stats": stats,
    }
    if model is not None:
        result.update(engine.extract(model, pool, n))
    print(json.dumps(result), flush=True)
    return 10 if satisfiable else 20


if __name__ == "__main__":
    sys.exit(main())
