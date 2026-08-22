"""EXP-009: is translation-invariance the obstruction, or is the code itself?

Theorem C4 says 83/318 catalogue PBB codes admit no *translation-invariant*
one-ancilla schedule.  That restriction could be the whole story.  Here we take
the smallest such codes and re-run the scheduler on the FULL, unrestricted
model (one independent time variable per (check, qubit) edge), which is a
strictly larger search space.

  * general model INFEASIBLE  -> the obstruction is intrinsic to the code, and
    Theorem C4 upgrades to an unconditional no-go for one-ancilla extraction.
  * general model FEASIBLE    -> translation-invariance was the binding
    restriction; the no-go must be stated for that schedule class only, and the
    achievable depth of the general schedule is the honest cost figure.

Either outcome is reportable; we must not assume the first.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import CircuitSpec, build_memory_circuit  # noqa: E402
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers)
from qec_research.codes.bicycle import PBBSpec  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
SCHED = ROOT / "results" / "processed" / "exp005_schedulability_n180.json"
OUT = ROOT / "results" / "processed"


def main() -> None:
    how_many = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600.0
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    extra = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    sched = {x["code_id"]: x for x in json.load(open(SCHED)) if x.get("code_id")}
    rows = {r.get("code_id"): r for r in map(json.loads, open(CATALOG))}
    bad = [cid for cid, x in sched.items() if x.get("orbit_proven_unschedulable")]
    bad.sort(key=lambda cid: (sched[cid]["n"], sched[cid]["nC"] + sched[cid]["nD"]))
    print(f"{len(bad)} orbit-unschedulable codes; probing the {how_many} smallest "
          f"with the unrestricted model", flush=True)

    out = []
    for cid in bad[:how_many]:
        r, s = rows[cid], sched[cid]
        spec = PBBSpec(r["ell"], r["m"],
                       [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                       [tuple(t) for t in (r["C_terms"] or [])],
                       [tuple(t) for t in (r["D_terms"] or [])])
        code, sup, orb = pbb_supports_and_orbits(spec)
        lb = depth_lower_bound(sup, code.n)
        nedges = sum(x.weight for x in sup)
        rec = {"code_id": cid, "n": r["n"], "k": r["k"], "d": r["d"],
               "nC": s["nC"], "nD": s["nD"], "lower_bound": lb,
               "num_edges": nedges, "orbit_unschedulable": True, "trace": []}
        found = None
        for T in range(lb, lb + extra + 1):
            t0 = time.time()
            res = cpsat_schedule(sup, code.n, T=T, time_limit_s=tl, workers=workers)
            ok = res.slot is not None and res.verification["valid"]
            rec["trace"].append({"T": T, "status": res.status, "valid": ok,
                                 "s": round(time.time() - t0, 1)})
            print(f"  {cid} n={r['n']} T={T}: {res.status} valid={ok} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            if ok:
                found = (T, res)
                break
        rec["general_depth"] = found[0] if found else None
        rec["general_schedulable"] = found is not None
        rec["general_proven_unschedulable"] = (
            found is None and all(t["status"] == "INFEASIBLE" for t in rec["trace"]))
        if found:
            T, res = found
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=4, p=0.0,
                basis="Z", layers=slots_to_layers(res.slot, T)))
            d0, o0 = circ.compile_detector_sampler().sample(4000, separate_observables=True)
            rec["noiseless_detector_firings"] = int(d0.sum())
            rec["noiseless_observable_flips"] = int(o0.sum())
            print(f"    -> general schedule at depth {T}; noiseless detectors "
                  f"{int(d0.sum())}", flush=True)
        out.append(rec)

    (OUT / "exp009_general_schedule_nogo.json").write_text(json.dumps(out, indent=2))
    ns = sum(1 for x in out if x["general_proven_unschedulable"])
    ok = sum(1 for x in out if x["general_schedulable"])
    print(f"\nSUMMARY over {len(out)} probed codes: "
          f"general-schedulable={ok}  general-proven-unschedulable={ns}")
    if ok:
        print("=> translation-invariance was the binding restriction for at least "
              "one code; Theorem C4 must remain scoped to that schedule class.")
    elif ns == len(out):
        print("=> obstruction is intrinsic on every probed code; the no-go is "
              "stronger than the translation-invariant class (on this sample).")


if __name__ == "__main__":
    main()
