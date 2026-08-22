"""EXP-017: is the decoder slowdown caused by BP failing to converge?

The isolated latency measurement (EXP-013) showed the PBB circuit costs ~10x
more per shot to decode while its detector error model is only 1.24x larger.
The natural explanation is that belief propagation converges less often on the
mixed-check Tanner graph, so more shots fall through to the ordered-statistics
stage.  That was an assertion; this measures it.

`ldpc.BpOsdDecoder` exposes `.converge` after each call, so the BP success rate
is directly observable.  We also record per-shot latency split by whether BP
converged, which is what turns the explanation into evidence.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis)
from qec_research.circuits.mixed_stabilizer import CircuitSpec, build_memory_circuit  # noqa: E402
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, depth_lower_bound, slots_to_layers)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.decoders.bposd_dem import dem_to_matrices  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "raw"

MAX_ITER, OSD_ORDER, OSD_METHOD = 30, 0, "osd0"


def build(code, sup, orb, rounds, p):
    lb = depth_lower_bound(sup, code.n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, code.n, T=T, time_limit_s=180, symmetry_orbits=orb,
                           random_seed=0)
        if r.slot is not None and r.verification["valid"]:
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=rounds,
                p=p, basis="Z", layers=slots_to_layers(r.slot, T)))
            meta["schedule_depth"] = T
            return circ, meta
    raise RuntimeError("no schedule")


def measure(circ, shots, seed):
    import scipy.sparse as sp
    from ldpc import BpOsdDecoder

    dem = circ.detector_error_model(decompose_errors=False)
    M = dem_to_matrices(dem)
    pr = np.clip(M.priors, 1e-12, 1 - 1e-12)
    dec = BpOsdDecoder(sp.csr_matrix(M.H), error_channel=list(pr), max_iter=MAX_ITER,
                       bp_method="ms", ms_scaling_factor=0.625,
                       osd_method=OSD_METHOD, osd_order=OSD_ORDER)
    smp = circ.compile_detector_sampler(seed=seed)
    det, obs = smp.sample(shots, separate_observables=True)
    det = det.astype(np.uint8)
    conv, lat = [], []
    for i in range(shots):
        t = time.perf_counter()
        dec.decode(det[i])
        lat.append((time.perf_counter() - t) * 1e3)
        conv.append(bool(getattr(dec, "converge", True)))
    conv = np.array(conv)
    lat = np.array(lat)
    ndet = det.sum(axis=1)
    out = {
        "shots": shots, "dem_errors": dem.num_errors, "dem_detectors": dem.num_detectors,
        "bp_converged_frac": float(conv.mean()),
        "mean_detection_events_per_shot": float(ndet.mean()),
        "latency_mean_ms": float(lat.mean()),
        "latency_p50_ms": float(np.percentile(lat, 50)),
        "latency_p99_ms": float(np.percentile(lat, 99)),
    }
    if conv.any():
        out["latency_mean_ms_bp_converged"] = float(lat[conv].mean())
    if (~conv).any():
        out["latency_mean_ms_bp_failed"] = float(lat[~conv].mean())
    return out


def main() -> None:
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    p = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002

    print(f"load average: {os.getloadavg()}", flush=True)
    res = {"env": {"shots": shots, "rounds": rounds, "p": p,
                   "decoder": f"BP min-sum {MAX_ITER} iter + {OSD_METHOD} order {OSD_ORDER}",
                   "load_at_start": os.getloadavg()}, "codes": []}

    c, s, o = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    circ, meta = build(c, s, o, rounds, p)
    r = measure(circ, shots, 7)
    res["codes"].append({"label": "CSS-BB [[144,12,12]] Gross",
                         "depth": meta["schedule_depth"], **r})

    rows = [json.loads(l) for l in open(CATALOG)]
    sch = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    dby = {x["code_id"]: x["orbit_depth"] for x in sch}
    fam = [x for x in rows if (x["n"], x["k"], x["d"]) == (144, 12, 12)
           and dby.get(x.get("code_id")) is not None]
    fam.sort(key=lambda x: dby[x["code_id"]])
    b = fam[0]
    spec = PBBSpec(b["ell"], b["m"], [tuple(t) for t in b["A_terms"]],
                   [tuple(t) for t in b["B_terms"]], [tuple(t) for t in b["C_terms"]],
                   [tuple(t) for t in b["D_terms"]])
    pc, ps, po = pbb_supports_and_orbits(spec)
    pcirc, pmeta = build(pc, ps, po, rounds, p)
    pr_ = measure(pcirc, shots, 7)
    res["codes"].append({"label": f"nonCSS-PBB [[144,12,12]] {b['code_id']}",
                         "depth": pmeta["schedule_depth"], **pr_})

    for x in res["codes"]:
        print(f"{x['label']:38s} BP converged {100*x['bp_converged_frac']:5.1f}% "
              f"| mean lat {x['latency_mean_ms']:8.1f} ms "
              f"(converged {x.get('latency_mean_ms_bp_converged', float('nan')):.1f}, "
              f"failed {x.get('latency_mean_ms_bp_failed', float('nan')):.1f}) "
              f"| dets/shot {x['mean_detection_events_per_shot']:.1f} "
              f"| dem_err {x['dem_errors']}", flush=True)
    a, bq = res["codes"]
    res["ratios"] = {
        "bp_failure_rate": ((1 - bq["bp_converged_frac"]) /
                            (1 - a["bp_converged_frac"])) if a["bp_converged_frac"] < 1 else None,
        "mean_latency": bq["latency_mean_ms"] / a["latency_mean_ms"],
        "dem_errors": bq["dem_errors"] / a["dem_errors"],
        "detection_events": (bq["mean_detection_events_per_shot"]
                             / a["mean_detection_events_per_shot"]),
    }
    print("\nPBB/CSS ratios:", {k: (round(v, 3) if v is not None else None)
                                for k, v in res["ratios"].items()})
    (OUT / "exp017_bp_convergence.json").write_text(json.dumps(res, indent=2, default=str))
    print("wrote results/raw/exp017_bp_convergence.json")


if __name__ == "__main__":
    main()
