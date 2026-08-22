"""EXP-004: certified syndrome-extraction depth for BB vs PBB.

For each code we
  1. compute the combinatorial lower bound  lb = max(max check weight, max qubit degree);
  2. ask CP-SAT for a translation-invariant schedule at T = lb, lb+1, ... and
     record whether each infeasibility is *proven* (status INFEASIBLE) rather
     than a timeout;
  3. verify the resulting schedule independently (ancilla/qubit disjointness
     plus the anticommuting-overlap parity rule);
  4. build the Stim circuit and confirm zero detection events at p = 0.

Output feeds the end-to-end resource comparison.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec, build_memory_circuit, schedule_stats)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers, verify_schedule)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "raw"
OUT.mkdir(parents=True, exist_ok=True)


def certified_min_depth(sup, n, orb, *, max_extra: int = 8, tl: float = 120.0,
                        workers: int = 12) -> dict:
    lb = depth_lower_bound(sup, n)
    trace = []
    for T in range(lb, lb + max_extra + 1):
        r = cpsat_schedule(sup, n, T=T, time_limit_s=tl, workers=workers,
                           symmetry_orbits=orb)
        trace.append({"T": T, "status": r.status,
                      "valid": bool(r.verification and r.verification["valid"])})
        if r.slot is not None and r.verification["valid"]:
            proven = all(t["status"] == "INFEASIBLE" for t in trace[:-1])
            return {"lower_bound_combinatorial": lb, "depth": T, "slot": r.slot,
                    "depth_is_certified_min_for_orbit_schedules": proven,
                    "trace": trace, "verification": r.verification}
    return {"lower_bound_combinatorial": lb, "depth": None, "slot": None,
            "depth_is_certified_min_for_orbit_schedules": False, "trace": trace}


def evaluate(code, sup, orb, label: str, tl: float) -> dict:
    t0 = time.time()
    n = code.n
    st = schedule_stats(code.H)
    res = certified_min_depth(sup, n, orb, tl=tl)
    rec = {"label": label, **st,
           "lower_bound": res["lower_bound_combinatorial"],
           "depth": res["depth"],
           "depth_certified_min": res["depth_is_certified_min_for_orbit_schedules"],
           "depth_search_trace": res["trace"]}
    if res["slot"] is not None:
        obs = pure_z_logical_basis(code)
        circ, meta = build_memory_circuit(CircuitSpec(
            H=code.H, observables=obs, rounds=4, p=0.0, basis="Z",
            layers=slots_to_layers(res["slot"], res["depth"])))
        d, o = circ.compile_detector_sampler().sample(4000, separate_observables=True)
        rec.update(noiseless_detector_firings=int(d.sum()),
                   noiseless_observable_flips=int(o.sum()),
                   num_observables=meta["num_observables"],
                   num_detectors=meta["num_detectors"],
                   schedule_verified=res["verification"]["valid"])
        rec["slot"] = [[int(c), int(q), int(t)] for (c, q), t in res["slot"].items()]
    rec["wall_time_s"] = round(time.time() - t0, 1)
    return rec


if __name__ == "__main__":
    tl = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    out = []

    for nm in ["[[72,12,6]]", "[[108,8,10]]", "[[144,12,12]]", "[[288,12,18]]"]:
        code, sup, orb = bb_supports_and_orbits(BRAVYI_BB[nm])
        rec = evaluate(code, sup, orb, f"BB {nm}", tl)
        rec["family"] = "CSS-BB"
        out.append(rec)
        print(f"{rec['label']:22s} lb={rec['lower_bound']} depth={rec['depth']} "
              f"certified={rec['depth_certified_min']} 2q={rec['total_two_qubit_gates']} "
              f"maxw={rec['max_check_weight']} mixed={rec['num_mixed_checks']} "
              f"det0={rec.get('noiseless_detector_firings')}", flush=True)

    rows = [json.loads(l) for l in open(CATALOG)]
    # survey PBB check-weight distribution first (cheap, no solver)
    wsurvey = Counter()
    for r in rows:
        w = len(r["A_terms"]) + len(r["B_terms"]) + len(r["C_terms"] or []) + len(r["D_terms"] or [])
        wsurvey[(r["n"], w)] += 1
    print("\nPBB catalogue (n, |A|+|B|+|C|+|D|) -> count:",
          dict(sorted(wsurvey.items())), flush=True)

    # pick, for each (n,k,d) headline family, the entry with the SMALLEST
    # perturbation weight -- the most favourable case for the PBB construction
    best: dict[tuple, dict] = {}
    for r in rows:
        key = (r["n"], r["k"], r["d"])
        w = len(r["C_terms"] or []) + len(r["D_terms"] or [])
        if key not in best or w < best[key]["_w"]:
            best[key] = {**r, "_w": w}

    targets = [k for k in best if k[0] <= 144]
    targets.sort(key=lambda t: (-t[2] * t[1] / t[0], t[0]))
    for key in targets[:12]:
        r = best[key]
        spec = PBBSpec(r["ell"], r["m"],
                       [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                       [tuple(t) for t in (r["C_terms"] or [])],
                       [tuple(t) for t in (r["D_terms"] or [])])
        code, sup, orb = pbb_supports_and_orbits(spec)
        lbl = f"PBB [[{r['n']},{r['k']},{r['d']}]] {r.get('code_id')}"
        rec = evaluate(code, sup, orb, lbl, tl)
        rec["family"] = "nonCSS-PBB"
        rec["catalog"] = {k: r[k] for k in ("code_id", "ell", "m", "n", "k", "d",
                                            "trust_level", "A_terms", "B_terms",
                                            "C_terms", "D_terms")}
        out.append(rec)
        print(f"{rec['label']:34s} lb={rec['lower_bound']} depth={rec['depth']} "
              f"certified={rec['depth_certified_min']} 2q={rec['total_two_qubit_gates']} "
              f"maxw={rec['max_check_weight']} mixed={rec['num_mixed_checks']} "
              f"det0={rec.get('noiseless_detector_firings')}", flush=True)

    (OUT / "exp004_circuit_cost_survey.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/raw/exp004_circuit_cost_survey.json")
