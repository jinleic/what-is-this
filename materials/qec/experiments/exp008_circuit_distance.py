"""EXP-008: circuit distance of the BB and PBB syndrome-extraction circuits.

Stim's ``search_for_undetectable_logical_errors`` returns a set of elementary
DEM errors whose combined effect flips an observable while firing no detector.
Its size is an UPPER BOUND on the circuit distance; it is exact only when the
search is unrestricted, which we cannot afford at n = 144.  Every number here
is therefore reported as an upper bound with the search restrictions stated.

We also report the graphlike distance (``shortest_graphlike_error``) where the
decomposition succeeds, which is a different, also-upper-bound quantity.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import CircuitSpec, build_memory_circuit  # noqa: E402
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "certificates"
OUT.mkdir(parents=True, exist_ok=True)


def build(code, sup, orb, rounds, p):
    lb = depth_lower_bound(sup, code.n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, code.n, T=T, time_limit_s=120, workers=12,
                           symmetry_orbits=orb)
        if r.slot is not None and r.verification["valid"]:
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=rounds,
                p=p, basis="Z", layers=slots_to_layers(r.slot, T)))
            meta["schedule_depth"] = T
            return circ, meta
    return None, None


def probe(circ: stim.Circuit, cap: int) -> dict:
    out: dict = {}
    t0 = time.time()
    try:
        err = circ.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=cap,
            dont_explore_edges_with_degree_above=cap,
            dont_explore_edges_increasing_symptom_degree=False,
            canonicalize_circuit_errors=True,
        )
        out["undetectable_logical_error_size"] = len(err)
        out["is_upper_bound"] = True
        out["search_cap"] = cap
    except Exception as e:
        out["undetectable_logical_error_size"] = None
        out["search_error"] = str(e)[:300]
    out["search_wall_s"] = round(time.time() - t0, 1)
    t1 = time.time()
    try:
        g = circ.shortest_graphlike_error(ignore_ungraphlike_errors=True)
        out["graphlike_error_size"] = len(g)
    except Exception as e:
        out["graphlike_error_size"] = None
        out["graphlike_error"] = str(e)[:200]
    out["graphlike_wall_s"] = round(time.time() - t1, 1)
    return out


def main() -> None:
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    out = []

    for nm in ["[[72,12,6]]", "[[144,12,12]]"]:
        code, sup, orb = bb_supports_and_orbits(BRAVYI_BB[nm])
        circ, meta = build(code, sup, orb, rounds, 1e-3)
        rec = {"label": f"CSS-BB {nm}", "family": "CSS-BB", "rounds": rounds,
               "depth": meta["schedule_depth"],
               "code_distance": 6 if "72" in nm else 12, **probe(circ, cap)}
        out.append(rec)
        print(f"{rec['label']:32s} depth={rec['depth']} d_code={rec['code_distance']} "
              f"circuit_d<={rec['undetectable_logical_error_size']} "
              f"graphlike<={rec['graphlike_error_size']} "
              f"[{rec['search_wall_s']}s]", flush=True)

    rows = [json.loads(l) for l in open(CATALOG)]
    sched = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    depth_by_id = {x["code_id"]: x["orbit_depth"] for x in sched}
    for nkd in [(72, 12, 6), (144, 12, 12)]:
        fam = [r for r in rows if (r["n"], r["k"], r["d"]) == nkd
               and depth_by_id.get(r.get("code_id")) is not None]
        if not fam:
            continue
        fam.sort(key=lambda r: depth_by_id[r["code_id"]])
        best = fam[0]
        spec = PBBSpec(best["ell"], best["m"],
                       [tuple(t) for t in best["A_terms"]], [tuple(t) for t in best["B_terms"]],
                       [tuple(t) for t in best["C_terms"]], [tuple(t) for t in best["D_terms"]])
        code, sup, orb = pbb_supports_and_orbits(spec)
        circ, meta = build(code, sup, orb, rounds, 1e-3)
        rec = {"label": f"nonCSS-PBB [[{nkd[0]},{nkd[1]},{nkd[2]}]] {best['code_id']}",
               "family": "nonCSS-PBB", "rounds": rounds,
               "depth": meta["schedule_depth"], "code_distance": nkd[2],
               "catalog_id": best.get("code_id"), **probe(circ, cap)}
        out.append(rec)
        print(f"{rec['label']:32s} depth={rec['depth']} d_code={rec['code_distance']} "
              f"circuit_d<={rec['undetectable_logical_error_size']} "
              f"graphlike<={rec['graphlike_error_size']} "
              f"[{rec['search_wall_s']}s]", flush=True)

    (OUT / f"exp008_circuit_distance_r{rounds}_cap{cap}.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote results/certificates/exp008_circuit_distance_r{rounds}_cap{cap}.json")


if __name__ == "__main__":
    main()
