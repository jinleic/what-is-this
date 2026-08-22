"""EXP-014: why mixed checks cost more -- the hook-error mechanism.

For a pure-X check measured with an ancilla control, a fault on the ancilla
midway through the gate sequence propagates only X errors onto the remaining
data targets.  Those are all in one sector and are seen by the Z-checks.

For a *mixed* check, the same fault propagates X onto the remaining X-targets
AND Z onto the remaining Z-targets, producing a single-fault data error with
support in both sectors.  Two measurements make this concrete:

(A) Hook spectrum.  For every check and every gate position, propagate a single
    ancilla fault forward and record the symplectic weight and the Pauli
    composition of the resulting data error.  Purely combinatorial; no
    simulation, no solver.

(B) Circuit distance.  Stim's undetectable-logical-error search on the actual
    circuit, at n = 72 where it terminates.  Reported as an UPPER BOUND with
    the search restrictions stated.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec, build_memory_circuit, generator_supports)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def hook_spectrum(H: np.ndarray, slot: dict[tuple[int, int], int]) -> dict:
    """Symplectic weight and Pauli composition of every single-ancilla hook.

    An ancilla fault occurring after the gates in time slots <= t has already
    delivered those Paulis; the *remaining* targets receive the propagated
    error.  We enumerate every (check, cut point) pair.
    """
    sup = generator_supports(H)
    n = H.shape[1] // 2
    weights: Counter = Counter()
    mixed_hooks = 0
    total = 0
    worst = 0
    for s in sup:
        order = sorted(s.paulis, key=lambda j: slot[(s.index, j)])
        for cut in range(1, len(order)):          # fault between gate cut-1 and cut
            rest = order[cut:]
            xs = sum(1 for j in rest if s.paulis[j] in ("X", "Y"))
            zs = sum(1 for j in rest if s.paulis[j] in ("Z", "Y"))
            w = len(rest)
            weights[w] += 1
            total += 1
            worst = max(worst, w)
            if xs and zs:
                mixed_hooks += 1
    return {
        "num_hooks": total,
        "num_two_sector_hooks": mixed_hooks,
        "frac_two_sector_hooks": mixed_hooks / total if total else 0.0,
        "max_hook_weight": worst,
        "mean_hook_weight": (sum(w * c for w, c in weights.items()) / total) if total else 0.0,
        "hook_weight_histogram": dict(sorted(weights.items())),
    }


def build(code, sup, orb, rounds, p):
    lb = depth_lower_bound(sup, code.n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, code.n, T=T, time_limit_s=120, workers=6,
                           symmetry_orbits=orb)
        if r.slot is not None and r.verification["valid"]:
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=rounds,
                p=p, basis="Z", layers=slots_to_layers(r.slot, T)))
            meta["schedule_depth"] = T
            return circ, meta, r.slot
    raise RuntimeError("no valid schedule")


def main() -> None:
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    do_search = (sys.argv[3].lower() != "no") if len(sys.argv) > 3 else True

    out = []
    code, sup, orb = bb_supports_and_orbits(BRAVYI_BB["[[72,12,6]]"])
    circ, meta, slot = build(code, sup, orb, rounds, 1e-3)
    rec = {"label": "CSS-BB [[72,12,6]]", "family": "CSS-BB", "code_distance": 6,
           "depth": meta["schedule_depth"],
           "two_qubit_gates_per_round": meta["total_two_qubit_gates"],
           "max_check_weight": meta["max_check_weight"],
           "num_mixed_checks": meta["num_mixed_checks"],
           **hook_spectrum(code.H, slot)}
    circuits = [(rec, circ)]

    rows = [json.loads(l) for l in open(CATALOG)]
    sch = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    dby = {x["code_id"]: x["orbit_depth"] for x in sch}
    fam = [r for r in rows if (r["n"], r["k"], r["d"]) == (72, 12, 6)
           and dby.get(r.get("code_id")) is not None]
    fam.sort(key=lambda r: dby[r["code_id"]])
    if fam:
        b = fam[0]
        spec = PBBSpec(b["ell"], b["m"], [tuple(t) for t in b["A_terms"]],
                       [tuple(t) for t in b["B_terms"]], [tuple(t) for t in b["C_terms"]],
                       [tuple(t) for t in b["D_terms"]])
        pc, ps, po = pbb_supports_and_orbits(spec)
        pcirc, pmeta, pslot = build(pc, ps, po, rounds, 1e-3)
        prec = {"label": f"nonCSS-PBB [[72,12,6]] {b['code_id']}", "family": "nonCSS-PBB",
                "code_distance": 6, "depth": pmeta["schedule_depth"],
                "two_qubit_gates_per_round": pmeta["total_two_qubit_gates"],
                "max_check_weight": pmeta["max_check_weight"],
                "num_mixed_checks": pmeta["num_mixed_checks"],
                **hook_spectrum(pc.H, pslot)}
        circuits.append((prec, pcirc))

    for rec, circ in circuits:
        print(f"{rec['label']:34s} depth={rec['depth']} 2q={rec['two_qubit_gates_per_round']} "
              f"mixed_checks={rec['num_mixed_checks']}", flush=True)
        print(f"    hooks={rec['num_hooks']} two-sector={rec['num_two_sector_hooks']} "
              f"({100*rec['frac_two_sector_hooks']:.1f}%) "
              f"max_wt={rec['max_hook_weight']} mean_wt={rec['mean_hook_weight']:.2f}", flush=True)
        if do_search:
            t0 = time.time()
            try:
                err = circ.search_for_undetectable_logical_errors(
                    dont_explore_detection_event_sets_with_size_above=cap,
                    dont_explore_edges_with_degree_above=cap,
                    dont_explore_edges_increasing_symptom_degree=False,
                    canonicalize_circuit_errors=True)
                rec["circuit_distance_upper_bound"] = len(err)
                rec["search_cap"] = cap
            except Exception as e:
                rec["circuit_distance_upper_bound"] = None
                rec["search_error"] = str(e)[:200]
            rec["search_wall_s"] = round(time.time() - t0, 1)
            print(f"    circuit distance <= {rec.get('circuit_distance_upper_bound')} "
                  f"(cap {cap}, code distance {rec['code_distance']}) "
                  f"[{rec['search_wall_s']}s]", flush=True)
        out.append(rec)

    (OUT / "exp014_hook_mechanism.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/processed/exp014_hook_mechanism.json")


if __name__ == "__main__":
    main()
