#!/usr/bin/env python3
"""Broad-convention ladder prober.

For a pair (n, T) build the all-degeneracy monolith (every parallel/concurrency
pattern free; per-line capacity + face budget + exact selection) and decide it.
On SAT, dump the combinatorial model so it can be straightened and verified
exactly; on UNSAT report the verdict (discovery-only until DRAT-verified).

Usage:  broad_ladder.py N T [--solver cadical195] [--dump]
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
from pysat.solvers import Solver  # noqa: E402


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def main() -> int:
    n, T = int(sys.argv[1]), int(sys.argv[2])
    solver = "cadical195"
    for i, a in enumerate(sys.argv):
        if a == "--solver":
            solver = sys.argv[i + 1]
    t0 = time.time()
    print(json.dumps({"stage": "build", "n": n, "T": T}), flush=True)
    cnf, pool = engine.build_model(n, T)
    tc = engine.add_triangle_crossing_indicators(cnf, pool, n)
    engine.add_face_bound(cnf, pool, n, T, crossing_lits=tc, per_line=True)
    engine.add_exact_selection(cnf, pool, n, T)
    print(json.dumps({"stage": "built", "vars": cnf.nv,
                      "clauses": len(cnf.clauses),
                      "seconds": round(time.time() - t0, 1),
                      "rss_gb": round(rss_gb(), 2)}), flush=True)
    t1 = time.time()
    with Solver(name=solver, bootstrap_with=cnf.clauses) as s:
        sat = s.solve()
        secs = round(time.time() - t1, 1)
        if not sat:
            print(json.dumps({"stage": "VERDICT", "n": n, "T": T,
                              "result": "UNSAT-discovery", "seconds": secs,
                              "stats": s.accum_stats()}), flush=True)
            return 20
        model = s.get_model()
    info = engine.extract(model, pool, n)
    out = {"stage": "VERDICT", "n": n, "T": T, "result": "SAT", "seconds": secs,
           "parallel_pairs": info["par"], "concurrent_triples": info["conc"],
           "selected": info["sel"]}
    print(json.dumps(out), flush=True)
    p = ROOT / "scratch" / "kobon" / "discoveries" / f"ladder_n{n}_t{T}_model.json"
    p.write_text(json.dumps(out, indent=1))
    print(json.dumps({"stage": "written", "path": str(p)}), flush=True)
    return 10


if __name__ == "__main__":
    sys.exit(main())
