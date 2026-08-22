"""EXP-005: which PBB codes admit a one-ancilla syndrome-extraction schedule?

A translation-invariant schedule assigns a time slot to each *monomial
direction*; there are  |A|+|B|+|C|+|D| (block 1) + |A|+|B| (block 2)
directions.  Once T reaches that count, every relative ordering of the
directions is realisable, so feasibility at T = #directions decides
feasibility at every larger T.  We therefore probe

    T = lb, lb+1, ..., #directions + 2

and treat "INFEASIBLE at all of them" as a proof that no translation-invariant
schedule exists.  For codes that fail we additionally probe the *general*
(non-symmetric) model, which is a strictly larger search space.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import pbb_supports_and_orbits  # noqa: E402
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound)
from qec_research.codes.bicycle import PBBSpec  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def probe(r: dict, tl: float, workers: int) -> dict:
    t0 = time.time()
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    code, sup, orb = pbb_supports_and_orbits(spec)
    n = code.n
    lb = depth_lower_bound(sup, n)
    ndir = len(set(orb.values()))
    best = None
    trace = []
    for T in range(lb, ndir + 3):
        res = cpsat_schedule(sup, n, T=T, time_limit_s=tl, workers=workers,
                             symmetry_orbits=orb)
        ok = res.slot is not None and res.verification["valid"]
        trace.append((T, res.status, ok))
        if ok:
            best = T
            break
    all_infeasible = all(s == "INFEASIBLE" for _, s, _ in trace)
    return {
        "code_id": r.get("code_id"), "n": r["n"], "k": r["k"], "d": r["d"],
        "ell": r["ell"], "m": r["m"],
        "nA": len(r["A_terms"]), "nB": len(r["B_terms"]),
        "nC": len(r["C_terms"] or []), "nD": len(r["D_terms"] or []),
        "lower_bound": lb, "num_directions": ndir,
        "orbit_depth": best,
        "orbit_schedulable": best is not None,
        "orbit_proven_unschedulable": (best is None and all_infeasible),
        "trace": trace,
        "wall_time_s": round(time.time() - t0, 1),
    }


def _job(a):
    r, tl, w = a
    try:
        return probe(r, tl, w)
    except Exception as e:
        return {"code_id": r.get("code_id"), "n": r.get("n"), "error": repr(e)}


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 180
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 9
    rows = [json.loads(l) for l in open(CATALOG) if json.loads(l)["n"] <= max_n]
    print(f"probing {len(rows)} PBB codes with n <= {max_n}", flush=True)
    out = []
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = [ex.submit(_job, (r, tl, 3)) for r in rows]
        for i, f in enumerate(as_completed(futs)):
            rec = f.result()
            out.append(rec)
            if i % 25 == 0:
                print(f"  {i+1}/{len(rows)}", flush=True)
    (OUT / f"exp005_schedulability_n{max_n}.json").write_text(json.dumps(out, indent=2))

    ok = [x for x in out if x.get("orbit_schedulable")]
    bad = [x for x in out if x.get("orbit_proven_unschedulable")]
    und = [x for x in out if not x.get("orbit_schedulable") and not x.get("orbit_proven_unschedulable")]
    print(f"\nschedulable={len(ok)}  proven-unschedulable={len(bad)}  undecided={len(und)}")

    print("\n--- correlation with perturbation weight parity ---")
    tab = defaultdict(lambda: [0, 0])
    for x in out:
        if "error" in x:
            continue
        w = x["nC"] + x["nD"]
        tab[w][0 if x.get("orbit_schedulable") else 1] += 1
    for w in sorted(tab):
        s, u = tab[w]
        print(f"  |C|+|D| = {w}:  schedulable {s:4d}   unschedulable {u:4d}")

    print("\n--- depth distribution among schedulable ---")
    print("  ", dict(sorted(Counter(x["orbit_depth"] for x in ok).items())))
