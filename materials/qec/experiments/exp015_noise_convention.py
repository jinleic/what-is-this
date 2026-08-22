"""EXP-015: reconcile our absolute LER with the artifact's noise convention.

Our matched benchmark applies DEPOLARIZE1(p) at *every* idle location, data and
ancilla alike.  The Bravyi artifact (and the independent reproduction built
from it) applies idle noise only at scheduled **data**-idle locations.  That is
a strictly weaker noise model, so our absolute logical error rates must be
higher.

This script quantifies the gap by rebuilding the *same* Gross circuit under
both conventions and measuring both, so the discrepancy is explained rather
than left dangling.  The CSS-vs-PBB comparison is unaffected either way,
because both codes are treated identically within each convention -- which this
script also checks by running the PBB code under both.
"""

from __future__ import annotations

import json
import sys
import time
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
from qec_research.decoders.parallel_harness import collect  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "raw"


def strip_ancilla_idle_noise(circ: stim.Circuit, n_data: int) -> stim.Circuit:
    """Return a copy in which DEPOLARIZE1 targets on ancillas are removed.

    This converts our 'idle noise everywhere' convention into the artifact's
    'idle noise on data qubits only' convention, changing nothing else.
    """
    out = stim.Circuit()
    for inst in circ:
        if inst.name == "DEPOLARIZE1":
            keep = [t.value for t in inst.targets_copy() if t.value < n_data]
            if keep:
                out.append("DEPOLARIZE1", keep, inst.gate_args_copy())
        else:
            out.append(inst)
    return out


def build(code, sup, orb, rounds, p):
    lb = depth_lower_bound(sup, code.n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, code.n, T=T, time_limit_s=120, workers=8,
                           symmetry_orbits=orb)
        if r.slot is not None and r.verification["valid"]:
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=rounds,
                p=p, basis="Z", layers=slots_to_layers(r.slot, T)))
            meta["schedule_depth"] = T
            return circ, meta
    raise RuntimeError("no schedule")


def main() -> None:
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    p = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 24

    targets = []
    c, s, o = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    targets.append(("CSS-BB [[144,12,12]] Gross", c, s, o))
    rows = [json.loads(l) for l in open(CATALOG)]
    sch = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    dby = {x["code_id"]: x["orbit_depth"] for x in sch}
    fam = [r for r in rows if (r["n"], r["k"], r["d"]) == (144, 12, 12)
           and dby.get(r.get("code_id")) is not None]
    fam.sort(key=lambda r: dby[r["code_id"]])
    b = fam[0]
    spec = PBBSpec(b["ell"], b["m"], [tuple(t) for t in b["A_terms"]],
                   [tuple(t) for t in b["B_terms"]], [tuple(t) for t in b["C_terms"]],
                   [tuple(t) for t in b["D_terms"]])
    pc, ps, po = pbb_supports_and_orbits(spec)
    targets.append((f"nonCSS-PBB [[144,12,12]] {b['code_id']}", pc, ps, po))

    res = {"p": p, "rounds": rounds, "shots": shots, "codes": []}
    for label, code, sup, orb in targets:
        circ_all, meta = build(code, sup, orb, rounds, p)
        circ_data = strip_ancilla_idle_noise(circ_all, code.n)
        d0 = circ_data.compile_detector_sampler().sample(0)  # sanity: compiles
        entry = {"label": label, "depth": meta["schedule_depth"],
                 "dem_errors_all_idles": circ_all.detector_error_model(
                     decompose_errors=False).num_errors,
                 "dem_errors_data_idles_only": circ_data.detector_error_model(
                     decompose_errors=False).num_errors, "runs": []}
        for tag, cc in (("all_idles_noisy (ours)", circ_all),
                        ("data_idles_only (artifact convention)", circ_data)):
            r = collect(cc, max_shots=shots, workers=workers, chunk=200)
            entry["runs"].append({"convention": tag, **{k: v for k, v in r.to_dict().items()
                                                        if k not in ("latency_p50_ms",
                                                                     "latency_p95_ms",
                                                                     "latency_p99_ms")}})
            print(f"{label:38s} {tag:40s} shots={r.shots} fails={r.failures} "
                  f"LER={r.ler_per_shot:.4e} CI=[{r.ci95_lo:.2e},{r.ci95_hi:.2e}]",
                  flush=True)
        a, bb_ = entry["runs"][0]["ler_per_shot"], entry["runs"][1]["ler_per_shot"]
        entry["ratio_all_over_data_only"] = (a / bb_) if bb_ else None
        res["codes"].append(entry)
        print(f"    dem errors: {entry['dem_errors_all_idles']} (ours) vs "
              f"{entry['dem_errors_data_idles_only']} (artifact convention); "
              f"LER ratio {entry['ratio_all_over_data_only']}", flush=True)

    if len(res["codes"]) == 2:
        g, q = res["codes"]
        for i, tag in enumerate(("ours", "artifact convention")):
            gl, ql = g["runs"][i]["ler_per_shot"], q["runs"][i]["ler_per_shot"]
            print(f"\nPBB/CSS LER ratio under {tag}: {ql/gl:.2f}" if gl else "")
            res.setdefault("pbb_over_css_ratio", {})[tag] = (ql / gl) if gl else None

    (OUT / "exp015_noise_convention.json").write_text(json.dumps(res, indent=2, default=str))
    print("\nwrote results/raw/exp015_noise_convention.json")


if __name__ == "__main__":
    main()
