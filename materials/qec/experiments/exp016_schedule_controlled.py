"""EXP-016: schedule-controlled matched benchmark (supersedes EXP-007).

Why this exists
---------------
EXP-007 re-solved the scheduling model on every call, unseeded, so the circuit
was not held fixed across noise points.  Schedule choice moves the logical
error rate by ~2x for the *same* code, so that sweep is confounded
(`notes/failed_routes.md` FR-007).

Design
------
The translation-invariant schedules at the minimum depth form a finite set that
CP-SAT can enumerate **exhaustively** (proved by status OPTIMAL, no truncation):
8496 for the Gross code at depth 7, 9968 for the PBB code at depth 8.  We
therefore

  1. enumerate the complete set,
  2. draw a **uniform random sample** of S schedules with a fixed RNG seed --
     unlike varying a solver seed, this is a genuine uniform sample,
  3. persist the full slot map of every sampled schedule (not just a hash: a
     solver upgrade would make a seed unrecoverable),
  4. build one circuit per (code, schedule) and reuse that exact slot map for
     both the p=0 validity check and the noisy run,
  5. report per-schedule counts with **Bonferroni-corrected simultaneous**
     Clopper-Pearson intervals at familywise 95 %, and a t-interval on the
     schedule-averaged LER.

Claims this licenses
--------------------
* *Primary*: the schedule-averaged LER of the two codes, with an interval that
  includes between-schedule variance.
* *Sampled separation*: every sampled PBB simultaneous lower bound exceeds
  every sampled CSS simultaneous upper bound.  This is a statement about the
  sample, made with familywise-corrected intervals; it is NOT a claim about the
  global best or worst schedule, which we do not measure.
"""

from __future__ import annotations

import hashlib
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
    cpsat_schedule, depth_lower_bound, enumerate_orbit_schedules,
    slots_to_layers, verify_schedule)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.decoders.bposd_dem import clopper_pearson  # noqa: E402
from qec_research.decoders.parallel_harness import collect  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "raw"
OUT.mkdir(parents=True, exist_ok=True)

MAX_ITER, OSD_ORDER, OSD_METHOD = 30, 0, "osd0"


def slot_hash(slot: dict) -> str:
    b = ",".join(f"{c}:{q}:{t}" for (c, q), t in sorted(slot.items())).encode()
    return hashlib.sha256(b).hexdigest()[:16]


def min_depth(sup, n, orb, tl=180.0) -> int | None:
    lb = depth_lower_bound(sup, n)
    ndir = len(set(orb.values()))
    for T in range(lb, ndir + 3):
        r = cpsat_schedule(sup, n, T=T, time_limit_s=tl, symmetry_orbits=orb,
                           random_seed=0)
        if r.slot is not None and r.verification["valid"]:
            return T
    return None


def main() -> None:
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    p = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 26
    S = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    rng_seed = int(sys.argv[6]) if len(sys.argv) > 6 else 20260811

    targets = []
    c, s, o = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    targets.append(("CSS-BB [[144,12,12]] Gross", c, s, o, None))
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
    targets.append((f"nonCSS-PBB [[144,12,12]] {b['code_id']}", pc, ps, po, b))

    # familywise 95%: 2 codes x S schedules simultaneous intervals
    alpha_fw = 0.05
    alpha_each = alpha_fw / (2 * S)

    res = {"env": {"date": "2026-08-11", "python": platform.python_version(),
                   "stim": stim.__version__, "rounds": rounds, "p": p,
                   "shots_per_schedule": shots, "workers": workers,
                   "schedules_sampled_per_code": S, "sampling_rng_seed": rng_seed,
                   "familywise_alpha": alpha_fw, "per_interval_alpha": alpha_each,
                   "decoder": f"BP min-sum {MAX_ITER} iter scale 0.625 + {OSD_METHOD} "
                              f"order {OSD_ORDER}, undecomposed hypergraph DEM",
                   "protocol": "complete enumeration of minimum-depth translation-invariant "
                               "schedules, then uniform random sample; full slot maps persisted; "
                               "the same slot map is used for the p=0 check and the noisy run"},
           "codes": []}

    rng = np.random.default_rng(rng_seed)
    for label, code, sup, orb, src in targets:
        T = min_depth(sup, code.n, orb)
        t0 = time.time()
        allsch, complete = enumerate_orbit_schedules(sup, code.n, T, orb,
                                                     limit=200000, time_limit_s=600)
        print(f"\n=== {label} : depth {T}, {len(allsch)} valid schedules "
              f"(enumeration complete={complete}) [{time.time()-t0:.0f}s] ===", flush=True)
        idx = rng.choice(len(allsch), size=min(S, len(allsch)), replace=False)
        entry = {"label": label, "catalog": src, "depth": T,
                 "total_valid_schedules": len(allsch),
                 "enumeration_complete": bool(complete),
                 "sampled_indices": [int(i) for i in idx], "runs": []}
        obs = pure_z_logical_basis(code)
        for i in idx:
            slot = allsch[int(i)]
            assert verify_schedule(sup, slot)["valid"]
            layers = slots_to_layers(slot, T)
            circ0, meta = build_memory_circuit(CircuitSpec(
                H=code.H, observables=obs, rounds=rounds, p=0.0, basis="Z", layers=layers))
            d0, o0 = circ0.compile_detector_sampler().sample(4000, separate_observables=True)
            assert int(d0.sum()) == 0 and int(o0.sum()) == 0, "invalid circuit"
            circ, _ = build_memory_circuit(CircuitSpec(
                H=code.H, observables=obs, rounds=rounds, p=p, basis="Z", layers=layers))
            r = collect(circ, max_shots=shots, workers=workers, chunk=200)
            lo, hi = clopper_pearson(r.failures, r.shots, alpha=alpha_each)
            rec = {"schedule_index": int(i), "schedule_hash": slot_hash(slot),
                   "slot_map": [[int(cc), int(qq), int(tt)] for (cc, qq), tt in sorted(slot.items())],
                   "depth": T, "two_qubit_gates_per_round": meta["total_two_qubit_gates"],
                   "max_check_weight": meta["max_check_weight"],
                   "num_mixed_checks": meta["num_mixed_checks"],
                   "noiseless_detector_firings": 0, "noiseless_observable_flips": 0,
                   "shots": r.shots, "failures": r.failures,
                   "ler_per_shot": r.ler_per_shot,
                   "ci_simultaneous_lo": lo, "ci_simultaneous_hi": hi,
                   "dem_errors": r.dem_errors, "wall_s": r.wall_s}
            entry["runs"].append(rec)
            print(f"  sched#{i:<5d} {rec['schedule_hash']} shots={r.shots} "
                  f"fails={r.failures:5d} LER={r.ler_per_shot:.4e} "
                  f"simCI=[{lo:.2e},{hi:.2e}] dem={r.dem_errors}", flush=True)
        lers = np.array([x["ler_per_shot"] for x in entry["runs"]])
        entry["ler_mean_over_schedules"] = float(lers.mean())
        entry["ler_sd_over_schedules"] = float(lers.std(ddof=1)) if len(lers) > 1 else 0.0
        if len(lers) > 1:
            from scipy import stats
            tcrit = float(stats.t.ppf(0.975, len(lers) - 1))
            half = tcrit * lers.std(ddof=1) / np.sqrt(len(lers))
            entry["ler_mean_ci95"] = [float(lers.mean() - half), float(lers.mean() + half)]
        entry["sampled_min_ler"] = float(lers.min())
        entry["sampled_max_ler"] = float(lers.max())
        entry["sampled_spread_ratio"] = float(lers.max() / lers.min()) if lers.min() else None
        entry["max_simultaneous_hi"] = max(x["ci_simultaneous_hi"] for x in entry["runs"])
        entry["min_simultaneous_lo"] = min(x["ci_simultaneous_lo"] for x in entry["runs"])
        print(f"  -> schedule-averaged LER {entry['ler_mean_over_schedules']:.3e} "
              f"CI95 {entry.get('ler_mean_ci95')}  sampled spread "
              f"{entry['sampled_spread_ratio']:.2f}x", flush=True)
        res["codes"].append(entry)

    if len(res["codes"]) == 2:
        css, pbb = res["codes"]
        sep = pbb["min_simultaneous_lo"] > css["max_simultaneous_hi"]
        point_sep = (
            min(r["ler_per_shot"] for r in pbb["runs"])
            > max(r["ler_per_shot"] for r in css["runs"])
        )
        res["verdict"] = {
            "css_max_simultaneous_upper_bound": css["max_simultaneous_hi"],
            "pbb_min_simultaneous_lower_bound": pbb["min_simultaneous_lo"],
            # FR-013 lesson: field names are contracts.  The familywise interval
            # test and the point-estimate rank separation are different claims
            # and get different, unambiguous names.
            "every_pair_separated_at_familywise_95": bool(sep),
            "every_sampled_pbb_point_estimate_exceeds_every_css": bool(point_sep),
            "familywise_confidence": 1 - alpha_fw,
            "schedule_averaged_ratio": (pbb["ler_mean_over_schedules"]
                                        / css["ler_mean_over_schedules"]),
            "scope": ("statement about the uniform sample of TRANSLATION-INVARIANT "
                      "minimum-depth schedules actually measured, with familywise-"
                      "corrected intervals; NOT a claim about the global "
                      "best or worst schedule, nor about non-TI schedules"),
        }
        print(f"\nVERDICT (familywise 95%): every sampled PBB lower bound "
              f"({pbb['min_simultaneous_lo']:.3e}) > every sampled CSS upper bound "
              f"({css['max_simultaneous_hi']:.3e}) ? {sep}")
        print(f"point-estimate rank separation (all PBB > all CSS)? {point_sep}")
        print(f"schedule-averaged LER ratio PBB/CSS = "
              f"{res['verdict']['schedule_averaged_ratio']:.2f}")

    (OUT / "exp016_schedule_controlled.json").write_text(json.dumps(res, indent=2, default=str))
    print("\nwrote results/raw/exp016_schedule_controlled.json")


if __name__ == "__main__":
    main()
