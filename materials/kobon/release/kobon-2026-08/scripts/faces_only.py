#!/usr/bin/env python3
"""Crossing-free (faces-only) monolith: sound by the collapse theorem.

Theorem (collapse).  For every arrangement A, the maximum interior-disjoint
family of triangles has the same size as the family of triangular faces.
Proof: corner descent replaces any crossed member by a triangular face
inside it (see scratch/kobon/collapse_verify.py for the machine check).

Consequence used here: "some arrangement of n lines carries T pairwise
interior-disjoint triangles" is equivalent to "some arrangement of n lines
carries T triangles none of which is crossed".  Forcing every TC literal
false is therefore a WLOG, and it removes the entire crossing degree of
freedom from the search while keeping all parallel/concurrency patterns free.

UNSAT  =>  K(n) < T for the classical Kobon problem.

Usage:  faces_only.py N T [--solve] [--out PATH] [--solver cadical195]
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


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def build(n, T):
    cnf, pool = engine.build_model(n, T)
    tc = engine.add_triangle_crossing_indicators(cnf, pool, n)
    for lit in tc:                       # WLOG: no selected triangle crossed
        cnf.append([-lit])
    engine.add_face_bound(cnf, pool, n, T, crossing_lits=tc, per_line=True)
    engine.add_exact_selection(cnf, pool, n, T)
    return cnf, pool, tc


def main() -> int:
    n, T = int(sys.argv[1]), int(sys.argv[2])
    out, solver, do_solve = None, "cadical195", "--solve" in sys.argv
    for i, a in enumerate(sys.argv):
        if a == "--out":
            out = Path(sys.argv[i + 1])
        if a == "--solver":
            solver = sys.argv[i + 1]
    t0 = time.time()
    print(json.dumps({"stage": "build", "n": n, "T": T}), flush=True)
    cnf, pool, tc = build(n, T)
    print(json.dumps({"stage": "built", "n": n, "T": T, "vars": cnf.nv,
                      "clauses": len(cnf.clauses), "tc_forced_false": len(tc),
                      "seconds": round(time.time() - t0, 1),
                      "rss_gb": round(rss_gb(), 2)}), flush=True)
    if out:
        cnf.to_file(str(out))
        print(json.dumps({"stage": "cnf", "path": str(out),
                          "bytes": out.stat().st_size}), flush=True)
    if not do_solve:
        return 0
    from pysat.solvers import Solver
    t1 = time.time()
    with Solver(name=solver, bootstrap_with=cnf.clauses) as s:
        sat = s.solve()
        secs = round(time.time() - t1, 1)
        if not sat:
            print(json.dumps({"stage": "VERDICT", "n": n, "T": T,
                              "result": "UNSAT", "K_upper": T - 1,
                              "seconds": secs,
                              "stats": s.accum_stats()}), flush=True)
            return 20
        model = s.get_model()
    info = engine.extract(model, pool, n)
    res = {"stage": "VERDICT", "n": n, "T": T, "result": "SAT",
           "seconds": secs, "parallel_pairs": info["par"],
           "concurrent_triples": info["conc"], "selected": info["sel"]}
    print(json.dumps(res), flush=True)
    p = ROOT / "scratch" / "kobon" / "discoveries" / f"facesonly_n{n}_t{T}.json"
    p.write_text(json.dumps(res, indent=1))
    print(json.dumps({"stage": "written", "path": str(p)}), flush=True)
    return 10


if __name__ == "__main__":
    sys.exit(main())
