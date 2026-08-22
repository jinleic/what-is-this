"""EXP-013: decoder latency, measured on an idle machine.

Supersedes the quarantined timing numbers of exp007 (see failed_routes FR-004).
Protocol:
  * refuse to start unless the machine load is below a threshold;
  * single worker process, one BLAS/OMP thread, pinned decoder settings;
  * discard a warm-up batch before timing;
  * identical shot count, identical decoder object lifetime for both codes;
  * interleave the two codes in alternating blocks so any residual drift hits
    both equally;
  * report p50/p95/p99 and throughput, plus the observed load average.
"""

from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

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
        r = cpsat_schedule(sup, code.n, T=T, time_limit_s=120, workers=4,
                           symmetry_orbits=orb)
        if r.slot is not None and r.verification["valid"]:
            circ, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=pure_z_logical_basis(code), rounds=rounds,
                p=p, basis="Z", layers=slots_to_layers(r.slot, T)))
            meta["schedule_depth"] = T
            return circ, meta
    raise RuntimeError("no schedule")


def make_decoder(circ):
    import scipy.sparse as sp
    from ldpc import BpOsdDecoder

    dem = circ.detector_error_model(decompose_errors=False)
    M = dem_to_matrices(dem)
    pr = np.clip(M.priors, 1e-12, 1 - 1e-12)
    dec = BpOsdDecoder(sp.csr_matrix(M.H), error_channel=list(pr),
                       max_iter=MAX_ITER, bp_method="ms", ms_scaling_factor=0.625,
                       osd_method=OSD_METHOD, osd_order=OSD_ORDER)
    return dec, M, dem


def timed(dec, M, circ, shots, seed):
    smp = circ.compile_detector_sampler(seed=seed)
    det, obs = smp.sample(shots, separate_observables=True)
    det = det.astype(np.uint8)
    lat = []
    fails = 0
    for i in range(shots):
        t = time.perf_counter()
        corr = dec.decode(det[i])
        lat.append((time.perf_counter() - t) * 1e3)
        if np.any(((M.L @ corr) % 2).astype(np.uint8) != obs[i].astype(np.uint8)):
            fails += 1
    return lat, fails


def main() -> None:
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    p = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002
    blocks = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    max_load = float(sys.argv[5]) if len(sys.argv) > 5 else 4.0

    load1 = os.getloadavg()[0]
    if load1 > max_load:
        print(f"REFUSING TO MEASURE: 1-minute load average is {load1:.1f} > {max_load}. "
              f"Wait for the machine to go idle.")
        sys.exit(2)
    print(f"load average at start: {os.getloadavg()}  ncpu={os.cpu_count()}", flush=True)

    code, sup, orb = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    circ_css, meta_css = build(code, sup, orb, rounds, p)

    rows = [json.loads(l) for l in open(CATALOG)]
    fam = [r for r in rows if r["n"] == 144 and r["k"] == 12 and r["d"] == 12]
    sch = json.load(open(ROOT / "results" / "processed" / "exp005_schedulability_n180.json"))
    dby = {x["code_id"]: x["orbit_depth"] for x in sch}
    fam = [r for r in fam if dby.get(r.get("code_id")) is not None]
    fam.sort(key=lambda r: dby[r["code_id"]])
    b = fam[0]
    spec = PBBSpec(b["ell"], b["m"], [tuple(t) for t in b["A_terms"]],
                   [tuple(t) for t in b["B_terms"]], [tuple(t) for t in b["C_terms"]],
                   [tuple(t) for t in b["D_terms"]])
    pcode, psup, porb = pbb_supports_and_orbits(spec)
    circ_pbb, meta_pbb = build(pcode, psup, porb, rounds, p)

    entries = [("CSS-BB [[144,12,12]] Gross", circ_css, meta_css, None),
               (f"nonCSS-PBB [[144,12,12]] {b['code_id']}", circ_pbb, meta_pbb, b["code_id"])]
    decs = {}
    for label, circ, meta, _ in entries:
        dec, M, dem = make_decoder(circ)
        decs[label] = (dec, M, circ, dem)
        timed(dec, M, circ, 40, seed=1)          # warm-up, discarded
    print("warm-up complete", flush=True)

    acc: dict[str, list[float]] = {label: [] for label, *_ in entries}
    per_block = max(1, shots // blocks)
    for bkt in range(blocks):                     # interleave the two codes
        for label, *_ in entries:
            dec, M, circ, dem = decs[label]
            lat, _ = timed(dec, M, circ, per_block, seed=1000 + bkt)
            acc[label].extend(lat)

    load_end = os.getloadavg()
    res = {"env": {"shots_per_code": sum(len(v) for v in acc.values()) // 2,
                   "rounds": rounds, "p": p, "blocks": blocks,
                   "workers": 1, "omp_threads": os.environ.get("OMP_NUM_THREADS"),
                   "decoder": f"BP min-sum {MAX_ITER} iter + {OSD_METHOD} order {OSD_ORDER}",
                   "load_at_start": load1, "load_at_end": load_end,
                   "ncpu": os.cpu_count(), "isolated": True},
           "codes": []}
    for label, circ, meta, cid in entries:
        a = np.asarray(acc[label])
        dem = decs[label][3]
        rec = {"label": label, "catalog_id": cid,
               "depth": meta["schedule_depth"],
               "two_qubit_gates_per_round": meta["total_two_qubit_gates"],
               "max_check_weight": meta["max_check_weight"],
               "num_mixed_checks": meta["num_mixed_checks"],
               "dem_errors": dem.num_errors, "dem_detectors": dem.num_detectors,
               "shots_timed": int(a.size),
               "latency_mean_ms": float(a.mean()),
               "latency_p50_ms": float(np.percentile(a, 50)),
               "latency_p95_ms": float(np.percentile(a, 95)),
               "latency_p99_ms": float(np.percentile(a, 99)),
               "throughput_shots_per_s_1core": float(1000.0 / a.mean())}
        res["codes"].append(rec)
        print(f"{label:38s} depth={rec['depth']} dem_err={rec['dem_errors']:6d} "
              f"p50={rec['latency_p50_ms']:7.1f}ms p95={rec['latency_p95_ms']:8.1f}ms "
              f"p99={rec['latency_p99_ms']:8.1f}ms  {rec['throughput_shots_per_s_1core']:.2f} shots/s/core",
              flush=True)
    a0 = res["codes"][0]
    a1 = res["codes"][1]
    res["ratios"] = {k: a1[k] / a0[k] for k in
                     ("latency_p50_ms", "latency_p95_ms", "latency_p99_ms", "dem_errors")}
    print("\nPBB / CSS ratios:", {k: round(v, 2) for k, v in res["ratios"].items()})
    print(f"load average at end: {load_end}")
    (OUT / "exp013_isolated_latency.json").write_text(json.dumps(res, indent=2))
    print("wrote results/raw/exp013_isolated_latency.json")


if __name__ == "__main__":
    main()
