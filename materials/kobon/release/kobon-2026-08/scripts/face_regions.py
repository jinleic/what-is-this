#!/usr/bin/env python3
"""Stronger faces-only Kobon CNF, sound by the convention-collapse theorem.

Unlike ``faces_only.py`` (whose active proof runs must remain reproducible),
this formulation uses ``build_model(..., rules=..., 'FACE')`` directly.  A
selected triple is constrained both by vertex sidedness and by consecutive
intersection order on each side.  It therefore needs no TC crossing variables.
The extra side-order clauses are redundant for realizable line arrangements,
but strengthen the SAT relaxation substantially.

UNSAT with a checked DRAT proof implies K(n) < T for arbitrary straight-line
arrangements, including parallels and multiple points.  SAT is only a
combinatorial candidate and still requires exact straight-line realization.

Usage:
  face_regions.py N T [--out FILE] [--solve] [--solver cadical195]
                      [--distilled] [--q1-no-concurrency]

``--distilled`` drops R1--R4.  This is still sound for UNSAT because
distinct arrangement faces are automatically interior-disjoint; it is a
weaker relaxation and is provided as a low-memory solver-diversity variant.

``--q1-no-concurrency`` fixes the unique parallel pair to the canonical
slope class {0,1}, fixes every other pair crossing, and forbids finite
multiple points.  With exactly one parallel class, rotating the slope circle
and relabelling puts that class first WLOG.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402


RULES = ("R1", "R2", "R3", "R4", "FACE")
DISTILLED_RULES = ("FACE",)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def build(n, target, distilled=False):
    rules = DISTILLED_RULES if distilled else RULES
    cnf, pool = engine.build_model(n, target, rules=rules)
    engine.add_face_bound(cnf, pool, n, target, per_line=True)
    engine.add_exact_selection(cnf, pool, n, target)
    return cnf, pool


def add_q1_no_concurrency(cnf, pool, n):
    pair = lambda i, j: pool.id(("P",) + tuple(sorted((i, j))))
    concurrent = lambda t: pool.id(("C",) + tuple(sorted(t)))
    for i in range(n):
        for j in range(i + 1, n):
            cnf.append([-pair(i, j) if (i, j) == (0, 1) else pair(i, j)])
    from itertools import combinations
    for triple in combinations(range(n), 3):
        cnf.append([-concurrent(triple)])


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    n, target = int(sys.argv[1]), int(sys.argv[2])
    out = None
    solver = "cadical195"
    solve = "--solve" in sys.argv
    distilled = "--distilled" in sys.argv
    q1_no_concurrency = "--q1-no-concurrency" in sys.argv
    for i, arg in enumerate(sys.argv):
        if arg == "--out":
            out = Path(sys.argv[i + 1])
        elif arg == "--solver":
            solver = sys.argv[i + 1]

    started = time.time()
    cnf, pool = build(n, target, distilled=distilled)
    if q1_no_concurrency:
        add_q1_no_concurrency(cnf, pool, n)
    print(json.dumps({
        "stage": "built",
        "formulation": (
            ("face-regions-distilled" if distilled else "face-regions")
            + ("-q1-no-concurrency" if q1_no_concurrency else "")
        ),
        "n": n,
        "target": target,
        "vars": cnf.nv,
        "clauses": len(cnf.clauses),
        "seconds": round(time.time() - started, 3),
        "rss_gb": round(rss_gb(), 3),
    }), flush=True)
    if out is not None:
        cnf.to_file(str(out))
        print(json.dumps({
            "stage": "written", "path": str(out), "bytes": out.stat().st_size
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
