#!/usr/bin/env python3
"""Universal direct-gap encoding for triangular faces of line arrangements.

A nondegenerate triple bounds a triangular arrangement face iff no outside
line intersects the relative interior of any of its three sides.  On side r
with endpoint lines a,b, an outside line u intersects that open side exactly
when its crossing on r is strictly between P_ra and P_rb.  The two forbidden
order patterns are encoded directly in X; a tie at an endpoint (concurrency)
is intentionally allowed.  Parallels and all finite multiple points remain
free.

Distinct arrangement faces automatically have disjoint interiors, so the
broad R1--R4 polygon-overlap rules and the SA/SB, B, IN and TC auxiliaries are
unnecessary.  The exact face budget and per-line capacity inequalities remain
as accelerators.  Every added clause is necessary for every real arrangement;
therefore UNSAT with a checked proof establishes K(n) < T unconditionally.
SAT remains a combinatorial candidate until exact straight-line realization.

Usage:  gap_faces.py N T [--out FILE] [--solve] [--solver cadical195]
                        [--global-only] [--sector-bounds]
                        [--sub-bound SIZE:UPPER]

``--global-only`` omits the per-line capacity counters while retaining the
exact bounded-face budget.  It is a weaker but sound low-memory variant.
``--sector-bounds`` adds local face-sector constraints at simple and multiple
vertices.  Each repeatable ``--sub-bound k:u`` adds the hereditary constraint
that every ``k``-line subarrangement supports at most ``u`` selected faces.
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


def add_direct_face_clauses(cnf, pool, n):
    order = lambda r, i, j: pool.id(("X", r, i, j))
    selected = lambda t: pool.id(("S",) + tuple(sorted(t)))
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


def build(
        n, target, per_line=True, subarrangement_bounds=None,
        sector_bounds=False):
    cnf, pool = engine.build_model(n, target, rules=())
    add_direct_face_clauses(cnf, pool, n)
    if sector_bounds:
        engine.add_selected_face_sector_bounds(cnf, pool, n)
    engine.add_face_bound(cnf, pool, n, target, per_line=per_line)
    engine.add_exact_selection(cnf, pool, n, target)
    if subarrangement_bounds:
        engine.add_subarrangement_face_bounds(
            cnf, pool, n, subarrangement_bounds, exact_target=target)
    return cnf, pool


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    n, target = int(sys.argv[1]), int(sys.argv[2])
    output = None
    solver = "cadical195"
    solve = "--solve" in sys.argv
    global_only = "--global-only" in sys.argv
    sector_bounds = "--sector-bounds" in sys.argv
    subarrangement_bounds = {}
    for i, arg in enumerate(sys.argv):
        if arg == "--out":
            output = Path(sys.argv[i + 1])
        elif arg == "--solver":
            solver = sys.argv[i + 1]
        elif arg == "--sub-bound":
            size, upper = map(int, sys.argv[i + 1].split(":"))
            subarrangement_bounds[size] = upper
    started = time.time()
    cnf, pool = build(
        n, target, per_line=not global_only,
        subarrangement_bounds=subarrangement_bounds,
        sector_bounds=sector_bounds)
    print(json.dumps({
        "stage": "built",
        "formulation": (
            "direct-gap-faces-global-only"
            if global_only else "direct-gap-faces"
        ) + ("+sectors" if sector_bounds else "") + (
            "+subarrangement" if subarrangement_bounds else ""),
        "sector_bounds": sector_bounds,
        "subarrangement_bounds": subarrangement_bounds,
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
