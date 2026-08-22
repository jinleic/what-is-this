"""EXP-007: matched circuit-level benchmark, CSS Gross vs non-CSS PBB.

Everything is matched except the code+circuit under test:
  * identical uniform circuit-level depolarising noise model and strength p
  * identical number of syndrome rounds
  * identical decoder (BP min-sum + OSD on the *undecomposed* hypergraph DEM,
    same iteration count, same OSD order, same scaling factor) -- so X/Z/Y
    correlations are preserved for the non-CSS code
  * identical observable convention: all k pure-Z logicals; a shot fails if
    any observable is mispredicted (block logical error)
  * identical shot budget and stopping rule
  * identical worker count, on one machine, for every timing number

Resource accounting includes data qubits, ancillas, two-qubit gate count and
the certified schedule depth, so the comparison is space-time complete.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import CircuitSpec, build_memory_circuit  # noqa: E402
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.decoders.parallel_harness import collect  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "raw"
OUT.mkdir(parents=True, exist_ok=True)

MAX_ITER, OSD_ORDER, OSD_METHOD = 30, 0, "osd0"


def find_schedule(sup, n, orb, tl=120.0):
    lb = depth_lower_bound(sup, n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, n, T=T, time_limit_s=tl, workers=12, symmetry_orbits=orb)
        if r.slot is not None and r.verification["valid"]:
            return T, r
    return None, None


def build(code, sup, orb, rounds, p, tl=120.0):
    T, r = find_schedule(sup, code.n, orb, tl)
    if T is None:
        return None, None
    obs = pure_z_logical_basis(code)
    circ, meta = build_memory_circuit(CircuitSpec(
        H=code.H, observables=obs, rounds=rounds, p=p, basis="Z",
        layers=slots_to_layers(r.slot, T)))
    meta["schedule_depth"] = T
    return circ, meta


def main() -> None:
    max_shots = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    ps = [float(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2
                             else ["0.001", "0.002", "0.003"])]
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 24
    max_errors = int(sys.argv[5]) if len(sys.argv) > 5 else 300

    targets = []
    code, sup, orb = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    targets.append(("CSS-BB [[144,12,12]] Gross", code, sup, orb, None))

    rows = [json.loads(l) for l in open(CATALOG)]
    fam = [r for r in rows if r["n"] == 144 and r["k"] == 12 and r["d"] == 12]
    sched = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    depth_by_id = {x["code_id"]: x["orbit_depth"] for x in sched}
    fam = [r for r in fam if depth_by_id.get(r.get("code_id")) is not None]
    fam.sort(key=lambda r: depth_by_id[r["code_id"]])
    for best in fam[:1]:                        # most favourable PBB member
        spec = PBBSpec(best["ell"], best["m"],
                       [tuple(t) for t in best["A_terms"]], [tuple(t) for t in best["B_terms"]],
                       [tuple(t) for t in best["C_terms"]], [tuple(t) for t in best["D_terms"]])
        pc, psu, po = pbb_supports_and_orbits(spec)
        targets.append((f"nonCSS-PBB [[144,12,12]] {best['code_id']}", pc, psu, po, best))

    env = {
        "date": "2026-08-11",
        "python": platform.python_version(), "platform": platform.platform(),
        "stim": stim.__version__, "workers": workers,
        "rounds": rounds, "max_shots": max_shots, "max_errors": max_errors,
        "p_values": ps,
        "decoder": f"BP min-sum ({MAX_ITER} iter, scale 0.625) + {OSD_METHOD} "
                   f"order {OSD_ORDER}, undecomposed hypergraph DEM",
        "noise": "uniform circuit-level depolarising: DEPOLARIZE2 after every 2q gate, "
                 "DEPOLARIZE1 on every idle location, X/Z_ERROR on reset and "
                 "before every measurement",
        "failure_criterion": "block logical error: any of the k=12 pure-Z observables mispredicted",
    }
    results = {"env": env, "codes": []}

    for label, c, s, o, src in targets:
        circ0, meta = build(c, s, o, rounds, 0.0)
        if circ0 is None:
            print(f"{label}: NO VALID SCHEDULE", flush=True)
            continue
        d0, o0 = circ0.compile_detector_sampler().sample(4000, separate_observables=True)
        entry = {"label": label, "meta": {k: v for k, v in meta.items() if k != "slot"},
                 "noiseless_detector_firings": int(d0.sum()),
                 "noiseless_observable_flips": int(o0.sum()),
                 "catalog": src, "points": []}
        print(f"\n=== {label} ===", flush=True)
        print(f"  data={meta['n_data']} anc={meta['n_ancilla']} depth={meta['schedule_depth']} "
              f"2q/round={meta['total_two_qubit_gates']} maxw={meta['max_check_weight']} "
              f"mixed={meta['num_mixed_checks']} obs={meta['num_observables']} "
              f"noiseless_dets={int(d0.sum())}", flush=True)
        for p in ps:
            circ, _ = build(c, s, o, rounds, p)
            r = collect(circ, max_shots=max_shots, max_errors=max_errors,
                        workers=workers, chunk=200, max_iter=MAX_ITER,
                        osd_order=OSD_ORDER, osd_method=OSD_METHOD)
            d = r.to_dict()
            d["p"] = p
            ler = r.ler_per_shot
            d["block_ler_per_round"] = 1 - (1 - ler) ** (1 / rounds) if ler < 1 else 1.0
            d["ci95_lo_per_round"] = 1 - (1 - r.ci95_lo) ** (1 / rounds)
            d["ci95_hi_per_round"] = 1 - (1 - r.ci95_hi) ** (1 / rounds)
            entry["points"].append(d)
            print(f"  p={p:<7g} shots={r.shots:6d} fails={r.failures:5d} "
                  f"LER/shot={ler:.4e} [{r.ci95_lo:.2e},{r.ci95_hi:.2e}] "
                  f"LER/round={d['block_ler_per_round']:.4e} "
                  f"p50={r.latency_p50_ms:.1f}ms p99={r.latency_p99_ms:.1f}ms "
                  f"wall={r.wall_s}s", flush=True)
        results["codes"].append(entry)

    (OUT / "exp007_matched_benchmark.json").write_text(json.dumps(results, indent=2, default=str))
    print("\nwrote results/raw/exp007_matched_benchmark.json")


if __name__ == "__main__":
    main()
