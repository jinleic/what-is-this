#!/usr/bin/env python3
"""Do crossings ever occur in an OPTIMAL broad-convention family?

For (n, T) build the all-degeneracy monolith with exact selection of T
triangles and additionally require at least one line-through-open-interior
incidence (an exact TC literal).  T is meant to be the certified value
Kgen(n).

UNSAT  =>  every optimal family at that n is crossing-free.  A selected
           triangle whose open interior meets no line IS a triangular face
           of the arrangement, so the family is a family of triangular
           faces and Kgen(n) = K(n) there: the crossing credit in the
           capacity ceiling is not realised at the optimum.
SAT    =>  an optimal family with a crossed triangle exists at pseudoline
           level; dump it for straightening.

Usage:  crossing_forced.py N T [--solver cadical195] [--dump-cnf PATH]
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
    solver, dump = "cadical195", None
    for i, a in enumerate(sys.argv):
        if a == "--solver":
            solver = sys.argv[i + 1]
        if a == "--dump-cnf":
            dump = Path(sys.argv[i + 1])
    t0 = time.time()
    print(json.dumps({"stage": "build", "n": n, "T": T}), flush=True)
    cnf, pool = engine.build_model(n, T)
    tc = engine.add_triangle_crossing_indicators(cnf, pool, n)
    engine.add_face_bound(cnf, pool, n, T, crossing_lits=tc, per_line=True)
    engine.add_exact_selection(cnf, pool, n, T)
    cnf.append(list(tc))          # at least one crossed selected triangle
    print(json.dumps({"stage": "built", "vars": cnf.nv,
                      "clauses": len(cnf.clauses), "tc_lits": len(tc),
                      "seconds": round(time.time() - t0, 1),
                      "rss_gb": round(rss_gb(), 2)}), flush=True)
    if dump:
        cnf.to_file(str(dump))
        print(json.dumps({"stage": "cnf", "path": str(dump)}), flush=True)
    t1 = time.time()
    with Solver(name=solver, bootstrap_with=cnf.clauses) as s:
        sat = s.solve()
        secs = round(time.time() - t1, 1)
        if not sat:
            print(json.dumps({"stage": "VERDICT", "n": n, "T": T,
                              "result": "UNSAT-no-optimal-crossing",
                              "seconds": secs,
                              "stats": s.accum_stats()}), flush=True)
            return 20
        model = s.get_model()
    info = engine.extract(model, pool, n)
    tvals = set(l for l in model if l > 0)
    inc = sorted(k[1:] for k, v in pool.obj2id.items()
                 if k[0] == "TC" and v in tvals)
    out = {"stage": "VERDICT", "n": n, "T": T, "result": "SAT",
           "seconds": secs, "parallel_pairs": info["par"],
           "concurrent_triples": info["conc"], "selected": info["sel"],
           "crossing_incidences": inc}
    print(json.dumps(out), flush=True)
    p = ROOT / "scratch" / "kobon" / "discoveries" / f"crossopt_n{n}_t{T}.json"
    p.write_text(json.dumps(out, indent=1))
    print(json.dumps({"stage": "written", "path": str(p)}), flush=True)
    return 10


if __name__ == "__main__":
    sys.exit(main())
