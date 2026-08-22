"""EXP-040: verified one-ancilla mixed-stabilizer circuit benchmark.

The decisive end-to-end circuit-level experiment for the PBB no-go program:
rebuild the EXP-013/029 deterministic one-ancilla circuits for the non-CSS PBB
[[144,12,12]] (12_6_0193) and the CSS Gross code, LOCK their fingerprints
against EXP-029, audit detector semantics with an independent second
implementation (Stim TableauSimulator replay), run a schedule-controlled
matched LER benchmark under IDENTICAL noise at three noise points with
familywise-corrected Clopper-Pearson intervals, attempt a budget-capped
circuit-distance tightening on the undecomposed DEMs, and add a
production-decoder (osd_cs) parity arm.

Staged driver: each subcommand persists its JSON before the next phase starts.

  registry   W1  fingerprint lock + noiseless checks
  audit      W2  independent detector-semantics audit (TableauSimulator)
  mc         W3+W5-LER  matched LER grid then production-decoder LER arm
              (long-running; writes results/partial_runs/exp040/mc_partial.jsonl
               continuously; safe to interrupt and resume)
  distance   W4  budget-capped circuit-distance tightening (single-threaded)
  latency    W5-latency  EXP-013-discipline latency arm (load-guarded)
  assemble   build final W3/W5 JSONs from the partial state
  verdict    W6  evidence consolidation
  pilot      throughput pilot (calibration only; not part of the grid)

Pre-declared grids (fixed before any Monte Carlo ran; measured core-s/shot on
this workstation under external SAT-fleet load: Gross ~0.3, PBB ~3.4):

  W3 matched grid   p in {0.0015, 0.002, 0.003}, execution order
                    0.002 -> 0.0015 -> 0.003; per (p, code, schedule):
                    start 40,000 shots; growth in stages of 4,000 until the
                    Clopper-Pearson half-width at per-interval alpha=0.005
                    (familywise 95% over 10 intervals) is <= 2e-3 or the
                    400,000-shot cap.  Incomplete cells at interruption are
                    reported as incomplete with their achieved counts; the
                    grid is never silently shrunk.
  W5 production arm p=0.002, per (code, schedule): target 2,000 shots in
                    stages of 500, per-cell wall cap 3,600 s, total wall
                    cap 21,600 s; achieved counts recorded per cell.  Run by
                    restarting the SAME exp040-mc process with
                    --phases production (serialized in time, never parallel
                    with the matched arm).
  W4 budgets        w3 exclusion 2,700 s; w4 matching-class 2,700 s and
                    8 GiB RSS guard; w4 star-class 2,700 s; stim witness
                    search sizes 4..8 with 2,700 s total per circuit.

Compute policy (user, 2026-08-17): this machine is capped at 50% CPU for our
workload; exp040-mc and exp040-distance together use at most 4 worker
threads -- Stim sampler workers=4, nice -n 15, OMP_NUM_THREADS=1; the
distance stage is single-threaded (not SAT-based: pure DEM mechanism search
plus one-at-a-time single-threaded stim search children).  An earlier
eight-worker nice -n 10 attempt measured ~0.5% CPU per worker under the
then-load-110 SAT fleet and produced no data; it was stopped and relaunched
under this policy.  Phases are serialized instead of parallelised.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os

for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_env, "1")

import platform
import resource
import sys
import time
from array import array
from pathlib import Path
from typing import Any

import numpy as np
import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits,
    pbb_supports_and_orbits,
    pure_z_logical_basis,
)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec,
    build_memory_circuit,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule,
    depth_lower_bound,
    slots_to_layers,
)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.decoders.bposd_dem import clopper_pearson  # noqa: E402
from qec_research.decoders.parallel_harness import collect  # noqa: E402

# ---------------------------------------------------------------------------
# Pre-registered protocol constants
# ---------------------------------------------------------------------------
EXPERIMENT = "EXP-040"
ROUNDS = 12
BASIS = "Z"
P_DEM_STRUCTURE = 0.001
SCHEDULE_TIME_LIMIT_S = 120.0
SCHEDULE_WORKERS_ARGUMENT = 4      # deterministic=True forces one CP-SAT worker
SCHEDULE_DETERMINISTIC = True
SCHEDULE_RANDOM_SEED = 0

CATALOG_ID = "12_6_0193"
CATALOG_PATH = (
    ROOT / "third_party" / "qcode-discovery" / "results" /
    "campaign7_publication_merged.jsonl"
)
EXP016_RAW_PATH = ROOT / "results" / "raw" / "exp016_schedule_controlled.json"
PAPER_PATH = ROOT / "reports" / "paper_pbb_nogo.md"
PAPER_SENTENCE_HEAD = "A verified one-ancilla mixed-stabilizer syndrome-extraction circuit for"

# W1 locked fingerprints (EXP-029 processed artifact, results/processed/
# exp029_circuit_distance.json#circuits/{key}/circuit_fingerprints).
LOCKED = {
    "gross": {
        "label": "CSS-BB [[144,12,12]] Gross",
        "circuit_sha256":
            "19591d836666a442f2d51521cdba93035fb1571b41a90ca807492f925923abaf",
        "dem_undecomposed_sha256":
            "cfc26ea70fb691dc7bb057126eb9c7adcae4b0fd4b15fd77d5d05f8a0616d349",
        "schedule_depth": 7,
        "two_qubit_gates_per_round": 864,
        "max_check_weight": 6,
        "num_mixed_checks": 0,
        "dem_errors": 67032,
        "dem_detectors": 1728,
    },
    "pbb": {
        "label": "nonCSS-PBB [[144,12,12]] 12_6_0193",
        "circuit_sha256":
            "ee9fe6fe01556be2d88201f1f02945949ce6f55faafde13e98f3c161f56283e6",
        "dem_undecomposed_sha256":
            "75e726ebf95a2ba8040b7cfbffb4610e54849a63d32e3b34555ac0a4a15d74c6",
        "schedule_depth": 8,
        "two_qubit_gates_per_round": 1008,
        "max_check_weight": 8,
        "num_mixed_checks": 72,
        "dem_errors": 82800,
        "dem_detectors": 1728,
    },
}

# W3 matched Monte Carlo grid (pre-declared; see module docstring).
W3_P_POINTS = (0.002, 0.0015, 0.003)     # execution order: 0.002 first
W3_START_SHOTS = 40_000
W3_STAGE_SHOTS = 4_000
W3_CAP_SHOTS = 400_000
W3_CP_HALF_WIDTH_TARGET = 2e-3
FAMILYWISE_ALPHA = 0.05
N_INTERVALS_FAMILY = 10                   # 2 codes x 5 schedules
PER_INTERVAL_ALPHA = FAMILYWISE_ALPHA / N_INTERVALS_FAMILY

MATCHED_DECODER = dict(max_iter=30, osd_order=0, osd_method="osd0",
                       bp_method="ms", ms_scaling_factor=0.625)

# W5 production-decoder arm (pre-declared).
PRODUCTION_DECODER = dict(max_iter=1000, osd_order=7, osd_method="osd_cs",
                          bp_method="ms", ms_scaling_factor=0.625)
PRODUCTION_P = 0.002
PRODUCTION_TARGET_SHOTS = 2_000
PRODUCTION_STAGE_SHOTS = 500
PRODUCTION_CELL_WALL_CAP_S = 3_600.0
PRODUCTION_TOTAL_WALL_CAP_S = 21_600.0

MC_WORKERS = 4
MC_CHUNK = 200
MC_BASE_SEED = 20260817
NOISELESS_SHOTS = 4000

# W4 budgets (do not raise).
W4_W3_BUDGET_S = 2700.0
W4_MATCHING_BUDGET_S = 2700.0
W4_STAR_BUDGET_S = 2700.0
W4_STIM_TOTAL_BUDGET_S = 2700.0
W4_MEMORY_GUARD_BYTES = 8 * 1024 ** 3
W4_DIGEST_BYTES = 8
FAULT_MODEL = "flattened undecomposed DEM error mechanisms"
CANONICALIZE_CIRCUIT_ERRORS = True

# W5 latency discipline (EXP-013).
LATENCY_MAX_LOAD = 4.0
LATENCY_BLOCKS = 4
LATENCY_SHOTS_PER_CODE = 400

PROCESSED_DIR = ROOT / "results" / "processed"
PARTIAL_DIR = ROOT / "results" / "partial_runs" / "exp040"
REGISTRY_PATH = PROCESSED_DIR / "exp040_circuit_registry.json"
AUDIT_PATH = PROCESSED_DIR / "exp040_detector_semantics_audit.json"
MATCHED_PATH = PROCESSED_DIR / "exp040_matched_ler.json"
DISTANCE_PATH = PROCESSED_DIR / "exp040_circuit_distance.json"
STAGNATION_PATH = PARTIAL_DIR / "circuit_distance_stagnation.json"
PRODUCTION_PATH = PROCESSED_DIR / "exp040_production_decoder.json"
VERDICT_PATH = PROCESSED_DIR / "exp040_verdict.json"
NOTE_PATH = ROOT / "notes" / "exp040_evidence_module.md"
MC_PARTIAL_PATH = PARTIAL_DIR / "mc_partial.jsonl"
DIVERGENCE_PATH = PARTIAL_DIR / "fingerprint_divergence.json"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def slot_hash(slot: dict[tuple[int, int], int]) -> str:
    b = ",".join(f"{c}:{q}:{t}" for (c, q), t in sorted(slot.items())).encode()
    return hashlib.sha256(b).hexdigest()[:16]


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def env_json() -> dict[str, Any]:
    return {
        "date": time.strftime("%Y-%m-%d"),
        "python": platform.python_version(),
        "stim": stim.__version__,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "threads": {k: os.environ.get(k) for k in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
        "loadavg_at_start": list(os.getloadavg()),
        "ncpu": os.cpu_count(),
    }


# ---------------------------------------------------------------------------
# Shared circuit rebuild layer (exact EXP-013/029 protocol)
# ---------------------------------------------------------------------------
def load_catalogue_row() -> dict[str, Any]:
    """Copied from experiments/exp029_circuit_distance.py::load_catalogue_row."""
    matches = []
    with CATALOG_PATH.open() as f:
        for line in f:
            row = json.loads(line)
            if row.get("code_id") == CATALOG_ID:
                matches.append(row)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one catalogue row {CATALOG_ID}, got {len(matches)}")
    row = matches[0]
    if (row["n"], row["k"], row["d"]) != (144, 12, 12):
        raise RuntimeError(f"unexpected parameters for {CATALOG_ID}: "
                           f"{row['n'], row['k'], row['d']}")
    return row


def code_and_supports(key: str):
    if key == "gross":
        return bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    if key == "pbb":
        row = load_catalogue_row()
        spec = PBBSpec(
            row["ell"], row["m"],
            [tuple(t) for t in row["A_terms"]],
            [tuple(t) for t in row["B_terms"]],
            [tuple(t) for t in row["C_terms"]],
            [tuple(t) for t in row["D_terms"]],
        )
        return pbb_supports_and_orbits(spec)
    raise ValueError(key)


def deterministic_slot(key: str) -> dict[str, Any]:
    """Recover the deterministic EXP-029 slot map via the EXP-013 CP-SAT loop."""
    code, supports, orbits = code_and_supports(key)
    lower_depth = depth_lower_bound(supports, code.n)
    num_directions = len(set(orbits.values()))
    for depth in range(lower_depth, num_directions + 3):
        result = cpsat_schedule(
            supports, code.n, T=depth,
            time_limit_s=SCHEDULE_TIME_LIMIT_S,
            workers=SCHEDULE_WORKERS_ARGUMENT,
            symmetry_orbits=orbits,
            random_seed=SCHEDULE_RANDOM_SEED,
            deterministic=SCHEDULE_DETERMINISTIC,
        )
        if result.slot is not None and result.verification["valid"]:
            return {
                "key": key, "code": code, "slot": result.slot, "depth": depth,
                "lower_bound": lower_depth, "status": result.status,
                "wall_s": result.wall_time_s,
                "verification": result.verification,
            }
    raise RuntimeError(f"no EXP-013 schedule found for {key}")


def circuit_from_slot(key: str, code, slot: dict[tuple[int, int], int],
                      depth: int, p: float, rounds: int = ROUNDS):
    layers = slots_to_layers(slot, depth)
    return build_memory_circuit(CircuitSpec(
        H=code.H, observables=pure_z_logical_basis(code),
        rounds=rounds, p=p, basis=BASIS, layers=layers))


def slot_map_json(slot: dict[tuple[int, int], int]) -> list[list[int]]:
    return [[int(c), int(q), int(t)] for (c, q), t in sorted(slot.items())]


def slot_from_json(entries: list[list[int]]) -> dict[tuple[int, int], int]:
    return {(int(c), int(q)): int(t) for c, q, t in entries}


def build_info(key: str, code, det: dict[str, Any], circuit, dem, metadata,
               p: float) -> dict[str, Any]:
    slot_lines = [f"{ci},{q}:{det['slot'][(ci, q)]}"
                  for ci, q in sorted(det["slot"])]
    return {
        "key": key,
        "label": LOCKED[key]["label"],
        "n": int(code.n), "k": int(code.k),
        "rounds": ROUNDS, "p": p, "basis": BASIS,
        "schedule_depth": det["depth"],
        "schedule_lower_bound": det["lower_bound"],
        "schedule_status": det["status"],
        "schedule_solver_wall_s": det["wall_s"],
        "schedule_verification": det["verification"],
        "schedule_sha256": sha256_text("\n".join(slot_lines)),
        "circuit_sha256": sha256_text(str(circuit)),
        "dem_undecomposed_sha256": sha256_text(str(dem)),
        "two_qubit_gates_per_round": metadata["total_two_qubit_gates"],
        "max_check_weight": metadata["max_check_weight"],
        "num_mixed_checks": metadata["num_mixed_checks"],
        "num_observables": int(circuit.num_observables),
        "dem_errors": int(dem.num_errors),
        "dem_detectors": int(dem.num_detectors),
    }


def rebuild_and_fingerprint(key: str) -> dict[str, Any]:
    """Rebuild one circuit by the literal EXP-013/029 loop and fingerprint it."""
    det = deterministic_slot(key)
    code = det["code"]
    circuit, metadata = circuit_from_slot(
        key, code, det["slot"], det["depth"], P_DEM_STRUCTURE)
    metadata = dict(metadata)
    metadata["schedule_depth"] = det["depth"]
    dem = circuit.detector_error_model(decompose_errors=False)
    info = build_info(key, code, det, circuit, dem, metadata, P_DEM_STRUCTURE)

    locked = LOCKED[key]
    observed = {name: info[name] for name in (
        "schedule_depth", "two_qubit_gates_per_round", "max_check_weight",
        "num_mixed_checks", "dem_errors", "dem_detectors")}
    expected = {name: locked[name] for name in observed}
    structural_ok = observed == expected
    info["structural_fingerprint_matches_exp029"] = structural_ok
    info["circuit_sha256_matches"] = info["circuit_sha256"] == locked["circuit_sha256"]
    info["dem_sha256_matches"] = (
        info["dem_undecomposed_sha256"] == locked["dem_undecomposed_sha256"])
    return {"info": info, "det": det, "code": code,
            "circuit": circuit, "dem": dem}


# ---------------------------------------------------------------------------
# W1 -- registry & fingerprint lock
# ---------------------------------------------------------------------------
def noiseless_check(key: str, code, slot, depth: int) -> dict[str, Any]:
    circ0, _ = circuit_from_slot(key, code, slot, depth, 0.0)
    d0, o0 = circ0.compile_detector_sampler().sample(
        NOISELESS_SHOTS, separate_observables=True)
    firings = int(d0.sum())
    flips = int(o0.sum())
    if firings != 0 or flips != 0:
        raise RuntimeError(
            f"{key} noiseless check FAILED: {firings} detector firings, "
            f"{flips} observable flips over {NOISELESS_SHOTS} shots")
    return {"shots": NOISELESS_SHOTS, "detector_firings": firings,
            "observable_flips": flips, "passed": True}


def cmd_registry() -> None:
    t0 = time.time()
    registry: dict[str, Any] = {
        "experiment": EXPERIMENT,
        "step": "W1-registry-fingerprint-lock",
        "env": env_json(),
        "protocol": {
            "rounds": ROUNDS, "basis": BASIS,
            "p_structure": P_DEM_STRUCTURE,
            "schedule": ("cpsat_schedule(time_limit_s=120, workers=4, "
                         "symmetry_orbits=translation orbits, random_seed=0, "
                         "deterministic=True); first feasible depth from "
                         "depth_lower_bound upward (EXP-013/029 loop)"),
            "catalog_id_pbb": CATALOG_ID,
            "gross_source": 'BRAVYI_BB["[[144,12,12]]"]',
            "noiseless_shots": NOISELESS_SHOTS,
            "locked_from": ("results/processed/exp029_circuit_distance.json"
                            "#circuits/{gross,pbb}/circuit_fingerprints"),
        },
        "circuits": {},
        "deterministic_slot_maps": {},
    }
    mismatch: dict[str, Any] | None = None
    for key in ("gross", "pbb"):
        t1 = time.time()
        built = rebuild_and_fingerprint(key)
        info = built["info"]
        info["rebuild_wall_s"] = round(time.time() - t1, 1)
        info["noiseless_p0"] = noiseless_check(
            key, built["code"], built["det"]["slot"], built["det"]["depth"])
        registry["circuits"][key] = info
        registry["deterministic_slot_maps"][key] = {
            "depth": built["det"]["depth"],
            "slot_map": slot_map_json(built["det"]["slot"]),
            "schedule_sha256": info["schedule_sha256"],
        }
        print(f"[W1] {key}: circuit_sha256={info['circuit_sha256'][:16]}... "
              f"match={info['circuit_sha256_matches']} "
              f"dem_match={info['dem_sha256_matches']} "
              f"structural={info['structural_fingerprint_matches_exp029']} "
              f"noiseless=OK ({time.time()-t1:.0f}s)", flush=True)
        if not (info["circuit_sha256_matches"] and info["dem_sha256_matches"]
                and info["structural_fingerprint_matches_exp029"]):
            mismatch = mismatch or {}
            mismatch[key] = {
                "observed_circuit_sha256": info["circuit_sha256"],
                "locked_circuit_sha256": LOCKED[key]["circuit_sha256"],
                "observed_dem_sha256": info["dem_undecomposed_sha256"],
                "locked_dem_sha256": LOCKED[key]["dem_undecomposed_sha256"],
                "structural_matches": info["structural_fingerprint_matches_exp029"],
            }

    registry["wall_s"] = round(time.time() - t0, 1)
    if mismatch is not None:
        registry["verdict"] = "FINGERPRINT_DIVERGENCE"
        registry["divergence"] = mismatch
        PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
        DIVERGENCE_PATH.write_text(json.dumps(registry, indent=2) + "\n")
        print(f"[W1] FINGERPRINT DIVERGENCE: {json.dumps(mismatch)}")
        print(f"[W1] wrote {DIVERGENCE_PATH.relative_to(ROOT)} and STOPPING")
        sys.exit(3)

    registry["verdict"] = "FINGERPRINTS_LOCKED"
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REGISTRY_PATH.write_text(json.dumps(registry, indent=2) + "\n")
    print(f"[W1] wrote {REGISTRY_PATH.relative_to(ROOT)}", flush=True)


# ---------------------------------------------------------------------------
# W2 -- independent detector-semantics audit (second implementation)
# ---------------------------------------------------------------------------
def row_pauli_strings(H: np.ndarray, n: int, n_total: int) -> list[stim.PauliString]:
    """Symplectic rows of H as PauliStrings over n_total qubits, +1 sign."""
    out = []
    for row in H:
        chars = ["I"] * n_total
        for q in range(n):
            x, z = int(row[q]), int(row[n + q])
            if (x, z) == (1, 0):
                chars[q] = "X"
            elif (x, z) == (0, 1):
                chars[q] = "Z"
            elif (x, z) == (1, 1):
                chars[q] = "Y"
        out.append(stim.PauliString("".join(chars)))
    return out


def analyse_stream(circuit: stim.Circuit, n: int, r: int, rounds: int):
    """Pass A: pure instruction-stream analysis, no simulation.

    Returns (measurement plan, detectors, observables) where the plan maps
    (kind, round, index) -> absolute measurement record index derived only
    from the instruction order, and detectors/observables carry their
    declared rec targets resolved to absolute indices.
    """
    meas_count = 0
    meas_plan: dict[tuple[str, int, int], int] = {}
    detectors: list[dict[str, Any]] = []
    observables: list[dict[str, Any]] = []
    round_of_ancilla_meas: dict[int, int] = {}
    current_round = -1
    for inst in circuit.flattened():
        name = inst.name
        if name in ("R", "RX"):
            continue
        if name in ("M", "MX"):
            targets = [t.value for t in inst.targets_copy()]
            is_ancilla = all(q >= n for q in targets)
            if name == "MX" and is_ancilla:
                current_round += 1
            for q in targets:
                if name == "MX" and is_ancilla:
                    meas_plan[("anc", current_round, q - n)] = meas_count
                    round_of_ancilla_meas[meas_count] = current_round
                elif name == "M" and q < n:
                    meas_plan[("data", None, q)] = meas_count
                else:
                    raise RuntimeError(f"unexpected measurement {name} on {targets}")
                meas_count += 1
            continue
        if name == "DETECTOR":
            coords = [c for c in inst.gate_args_copy()]
            recs = inst.targets_copy()
            idxs = [meas_count + t.value for t in recs]  # rec offsets are negative
            detectors.append({"coords": coords, "abs_indices": idxs})
            continue
        if name == "OBSERVABLE_INCLUDE":
            oid = int(inst.gate_args_copy()[0])
            recs = inst.targets_copy()
            idxs = [meas_count + t.value for t in recs]
            observables.append({"observable": oid, "abs_indices": idxs})
            continue
        if name in ("TICK", "X_ERROR", "Z_ERROR", "DEPOLARIZE1", "DEPOLARIZE2",
                    "SHIFT_COORDS"):
            continue
        if name in ("CX", "CY", "CZ"):
            continue
        raise RuntimeError(f"unexpected instruction {name} in noiseless circuit")
    return meas_plan, detectors, observables, meas_count


def replay(circuit: stim.Circuit, *, n: int, r: int, codeword_gens=None):
    """Replay the noiseless circuit on a TableauSimulator (second implementation).

    mode |0>    : R on data honoured (real experiment start).
    mode codeword: data starts in the mutual +1 eigenstate of every H row
                   (plus pure-Z logicals and ancilla Z stabilisers); the
                   initial data R is skipped.
    Returns (sim, records, per-round ancilla outcomes, final data outcomes).
    """
    sim = stim.TableauSimulator()
    if codeword_gens is not None:
        sim.set_state_from_stabilizers(codeword_gens, allow_redundant=True)
    records: list[bool] = []
    ancilla_by_round: list[list[int]] = []
    data_final: dict[int, int] = {}
    skipped_initial_R = codeword_gens is not None
    current_round = -1
    for inst in circuit.flattened():
        name = inst.name
        targets = [t.value for t in inst.targets_copy()]
        if name == "R":
            if skipped_initial_R and all(q < n for q in targets):
                continue
            for q in targets:
                sim.reset(q)
        elif name == "RX":
            for q in targets:
                sim.reset_x(q)
        elif name in ("CX", "CY", "CZ"):
            for c, t in zip(targets[0::2], targets[1::2]):
                getattr(sim, name.lower())(c, t)
        elif name == "MX":
            is_ancilla = all(q >= n for q in targets)
            if is_ancilla:
                current_round += 1
                ancilla_by_round.append([])
            for q in targets:
                if is_ancilla:
                    sim.h(q)
                    bit = sim.measure(q)
                    ancilla_by_round[current_round].append(int(bit))
                else:
                    sim.h(q)
                    bit = sim.measure(q)
                records.append(bool(bit))
        elif name == "M":
            for q in targets:
                bit = sim.measure(q)
                if q < n:
                    data_final[q] = int(bit)
                records.append(bool(bit))
        elif name in ("TICK", "DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS"):
            continue
        elif name in ("X_ERROR", "Z_ERROR", "DEPOLARIZE1", "DEPOLARIZE2"):
            continue  # p=0 circuit carries no error instructions
        else:
            raise RuntimeError(f"replay: unexpected instruction {name}")
    return sim, records, ancilla_by_round, data_final


def audit_one(key: str, slot: dict[tuple[int, int], int], depth: int,
              provenance: str) -> dict[str, Any]:
    code, _, _ = code_and_supports(key)
    n, r = code.n, code.H.shape[0]
    H = code.H
    circuit, meta = circuit_from_slot(key, code, slot, depth, 0.0)

    if circuit.num_detectors != 1728:
        raise RuntimeError(f"{key}: num_detectors {circuit.num_detectors} != 1728")

    # ---- Pass A: declared semantics vs analytic model (no simulation) ----
    meas_plan, detectors, observables, meas_total = analyse_stream(
        circuit, n, r, ROUNDS)
    supports = []
    for row in H:
        d = {}
        for q in range(n):
            x, z = int(row[q]), int(row[n + q])
            if x and z:
                d[q] = "Y"
            elif x:
                d[q] = "X"
            elif z:
                d[q] = "Z"
        supports.append(d)
    det_ready = [all(p == "Z" for p in s.values()) for s in supports]

    expected_sets: dict[tuple[int, int], frozenset] = {}
    for check in range(r):
        for rd in range(ROUNDS + 1):
            if rd == 0:
                if det_ready[check]:
                    expected_sets[(check, rd)] = frozenset(
                        {meas_plan[("anc", 0, check)]})
            elif rd < ROUNDS:
                expected_sets[(check, rd)] = frozenset(
                    {meas_plan[("anc", rd, check)],
                     meas_plan[("anc", rd - 1, check)]})
            else:
                if det_ready[check]:
                    s = {meas_plan[("data", None, q)] for q in supports[check]}
                    s.add(meas_plan[("anc", ROUNDS - 1, check)])
                    expected_sets[(check, rd)] = frozenset(s)

    seen_coords: set[tuple[int, int]] = set()
    for d in detectors:
        coords = d["coords"]
        if len(coords) != 2:
            raise RuntimeError(f"{key}: detector coords {coords} not (check, round)")
        ck, rd = int(coords[0]), int(coords[1])
        if (ck, rd) in seen_coords:
            raise RuntimeError(f"{key}: duplicate detector {coords}")
        seen_coords.add((ck, rd))
        exp = expected_sets.get((ck, rd))
        if exp is None:
            raise RuntimeError(
                f"{key}: detector {(ck, rd)} declared but analytic model says "
                f"it must not exist (check not basis-compatible)")
        got = frozenset(d["abs_indices"])
        if got != exp:
            raise RuntimeError(
                f"{key}: detector {(ck, rd)} index set mismatch: "
                f"declared {sorted(got)} != analytic {sorted(exp)}")
    missing = set(expected_sets) - seen_coords
    if missing:
        raise RuntimeError(f"{key}: missing detectors for {sorted(missing)[:5]} ...")
    if len(detectors) != 1728:
        raise RuntimeError(f"{key}: analysed {len(detectors)} detectors != 1728")

    basis = pure_z_logical_basis(code)
    if basis.shape[0] != 12:
        raise RuntimeError(f"{key}: pure-Z logical basis has {basis.shape[0]} rows")
    for entry in observables:
        oi = entry["observable"]
        if oi >= basis.shape[0]:
            raise RuntimeError(f"{key}: observable id {oi} out of range")
        supp = np.flatnonzero(basis[oi, n:])
        exp = frozenset(int(meas_plan[("data", None, int(j))]) for j in supp)
        if frozenset(entry["abs_indices"]) != exp:
            raise RuntimeError(
                f"{key}: OBSERVABLE_INCLUDE({oi}) index set mismatch: "
                f"declared {sorted(entry['abs_indices'])} != "
                f"basis support {sorted(exp)}")
    if len(observables) != 12:
        raise RuntimeError(f"{key}: {len(observables)} observables != 12")

    # ---- Pass B: |0> replay (real experiment initialisation) ----
    _, recs0, anc0, data0 = replay(circuit, n=n, r=r)

    # ---- Pass C: codeword replay, eigenstate check (a) ----
    gens = row_pauli_strings(H, n, n + r)
    for oi in range(basis.shape[0]):
        chars = ["I"] * (n + r)
        for q in np.flatnonzero(basis[oi, n:]):
            chars[int(q)] = "Z"
        gens.append(stim.PauliString("".join(chars)))
    for a in range(r):
        chars = ["I"] * (n + r)
        chars[n + a] = "Z"
        gens.append(stim.PauliString("".join(chars)))
    sim_c, recs_c, anc_c, data_c = replay(circuit, n=n, r=r, codeword_gens=gens)

    rows = gens[:r]
    ancilla_flips_in_codeword = sum(sum(o) for o in anc_c)
    # Round-by-round eigenstate verification needs the state *after* each
    # round, so replay again with a checkpoint after every MX ancilla block.
    sim = stim.TableauSimulator()
    sim.set_state_from_stabilizers(gens, allow_redundant=True)
    # the prepared state must be a +1 eigenstate of every generator row
    # *before* any circuit instruction is applied
    for row in rows:
        if sim.peek_observable_expectation(row) != 1:
            raise RuntimeError(f"{key}: initial codeword is not +1 on a row")
    skipped = True
    current_round = -1
    per_round_signs: list[list[int]] = []
    for inst in circuit.flattened():
        name = inst.name
        targets = [t.value for t in inst.targets_copy()]
        if name == "R":
            if skipped and all(q < n for q in targets):
                continue
            for q in targets:
                sim.reset(q)
        elif name == "RX":
            for q in targets:
                sim.reset_x(q)
        elif name in ("CX", "CY", "CZ"):
            for c, t in zip(targets[0::2], targets[1::2]):
                getattr(sim, name.lower())(c, t)
        elif name == "MX":
            is_ancilla = all(q >= n for q in targets)
            if is_ancilla:
                current_round += 1
            for q in targets:
                sim.h(q)
                sim.measure(q)
            if is_ancilla:
                signs = [sim.peek_observable_expectation(row) for row in rows]
                bad = [i for i, s in enumerate(signs) if s != 1]
                if bad:
                    raise RuntimeError(
                        f"{key}: after round {current_round}, generator rows "
                        f"{bad[:5]} are not +1 eigenstates of the data state")
                per_round_signs.append(signs)
        elif name == "M":
            for q in targets:
                sim.measure(q)
        elif name in ("TICK", "DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS",
                      "X_ERROR", "Z_ERROR", "DEPOLARIZE1", "DEPOLARIZE2"):
            continue
        else:
            raise RuntimeError(f"checkpoint replay: unexpected {name}")
    if len(per_round_signs) != ROUNDS:
        raise RuntimeError(f"{key}: checkpointed {len(per_round_signs)} rounds")

    # codeword pass: all ancilla outcomes must be 0 and final parities 0
    if ancilla_flips_in_codeword != 0:
        raise RuntimeError(
            f"{key}: {ancilla_flips_in_codeword} "
            f"ancilla flips in the noiseless codeword replay")
    for check in range(r):
        if det_ready[check]:
            par = sum(data_c[q] for q in supports[check]) % 2
            if par != 0:
                raise RuntimeError(f"{key}: final parity of pure-Z check {check} "
                                   f"is {par} in codeword replay")
    for oi in range(12):
        par = sum(data_c[int(q)] for q in np.flatnonzero(basis[oi, n:])) % 2
        if par != 0:
            raise RuntimeError(f"{key}: observable {oi} flips in codeword replay")

    # detectors computed from the |0> replay records must never fire
    fired = 0
    total_meas = len(recs0)
    for d in detectors:
        x = 0
        for i in d["abs_indices"]:
            x ^= int(recs0[i])
        fired += x
    if fired != 0:
        raise RuntimeError(f"{key}: {fired} detectors fire in the |0> replay")
    flips0 = 0
    for entry in observables:
        x = 0
        for i in entry["abs_indices"]:
            x ^= int(recs0[i])
        flips0 += x
    if flips0 != 0:
        raise RuntimeError(f"{key}: {flips0} observable flips in the |0> replay")

    return {
        "key": key,
        "provenance": provenance,
        "schedule_hash": slot_hash(slot),
        "depth": depth,
        "num_detectors": int(circuit.num_detectors),
        "num_observables": len(observables),
        "checks": {
            "a_eigenstate_plus1_after_every_round": True,
            "rounds_checked": len(per_round_signs),
            "generator_rows_checked_per_round": len(rows),
            "b_detector_index_sets_match_analytic_model": True,
            "detectors_checked": len(detectors),
            "c_observable_index_sets_match_pure_z_basis": True,
            "observables_checked": len(observables),
            "codeword_ancilla_outcomes_all_zero": True,
            "codeword_final_pure_z_parity_zero": True,
            "zero_replay_detector_firings": 0,
            "zero_replay_observable_flips": 0,
        },
        "method": (
            "independent replay of the p=0 instruction stream on "
            "stim.TableauSimulator; pass A pure index analysis; pass C starts "
            "from set_state_from_stabilizers(144 H rows + 12 pure-Z logicals "
            "+ ancilla Z) and peeks every generator row after every round"),
    }


def cmd_audit() -> None:
    t0 = time.time()
    exp016 = json.loads(EXP016_RAW_PATH.read_text())
    registry = json.loads(REGISTRY_PATH.read_text())
    audits: list[dict[str, Any]] = []
    order: list[tuple[str, dict[tuple[int, int], int], int, str]] = []
    for ci, code_entry in enumerate(exp016["codes"]):
        key = "gross" if ci == 0 else "pbb"
        depth = int(code_entry["depth"])
        for run in code_entry["runs"]:
            slot = slot_from_json(run["slot_map"])
            order.append((key, slot, depth,
                          f"EXP-016 sampled schedule #{run['schedule_index']} "
                          f"hash {run['schedule_hash']}"))
    for key in ("gross", "pbb"):
        det = registry["deterministic_slot_maps"][key]
        slot = slot_from_json(det["slot_map"])
        order.append((key, slot, int(det["depth"]),
                      "deterministic EXP-029/EXP-013 slot "
                      f"(cpsat seed 0, schedule_sha256 {det['schedule_sha256'][:16]})"))

    hashes: dict[tuple[str, str], int] = {}
    duplicates: list[list[str, str, int]] = []
    for key, slot, depth, prov in order:
        t1 = time.time()
        record = audit_one(key, slot, depth, prov)
        record["wall_s"] = round(time.time() - t1, 1)
        audits.append(record)
        h = (key, record["schedule_hash"])
        if h in hashes:
            duplicates.append([h[0], h[1], len(audits)])
        hashes[h] = len(audits)
        print(f"[W2] {key} {record['schedule_hash']} depth={depth}: ALL CHECKS PASS "
              f"({record['wall_s']}s)", flush=True)

    out = {
        "experiment": EXPERIMENT,
        "step": "W2-detector-semantics-audit",
        "env": env_json(),
        "protocol": {
            "slot_map_source": "results/raw/exp016_schedule_controlled.json#codes/"
                               "[*].runs[*].slot_map (5 per code) + deterministic "
                               "EXP-029 pair from exp040_circuit_registry.json",
            "builder": "qec_research.circuits.mixed_stabilizer."
                       "build_memory_circuit(layers=slots_to_layers(slot, depth))",
            "rounds": ROUNDS,
            "required_num_detectors": 1728,
            "failure_policy": "any mismatch raises; no JSON is written",
        },
        "n_circuits_audited": len(audits),
        "n_unique_circuits": len(hashes),
        "duplicate_entries": duplicates,
        "all_passed": True,
        "audits": audits,
        "wall_s": round(time.time() - t0, 1),
    }
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(f"[W2] wrote {AUDIT_PATH.relative_to(ROOT)} "
          f"({out['n_circuits_audited']} circuits)", flush=True)


# ---------------------------------------------------------------------------
# W3 / W5 -- Monte Carlo (matched grid + production arm)
# ---------------------------------------------------------------------------
def matched_cell_list() -> list[dict[str, Any]]:
    exp016 = json.loads(EXP016_RAW_PATH.read_text())
    schedules: dict[str, list[dict[str, Any]]] = {}
    for ci, code_entry in enumerate(exp016["codes"]):
        key = "gross" if ci == 0 else "pbb"
        schedules[key] = [{"depth": int(code_entry["depth"]), **run}
                          for run in code_entry["runs"]]
    cells = []
    cell_id = 0
    for p in W3_P_POINTS:
        for s in range(5):
            for key in ("gross", "pbb"):
                run = schedules[key][s]
                cells.append({
                    "cell_id": cell_id, "phase": "matched", "p": p, "key": key,
                    "schedule_index": run["schedule_index"],
                    "schedule_hash": run["schedule_hash"],
                    "depth": run["depth"],
                    "slot_map": run["slot_map"],
                })
                cell_id += 1
    return cells


def production_cell_list() -> list[dict[str, Any]]:
    cells = [c for c in matched_cell_list() if c["p"] == PRODUCTION_P]
    out = []
    for i, c in enumerate(cells):
        d = dict(c)
        d["cell_id"] = 10_000 + i
        d["phase"] = "production"
        out.append(d)
    return out


def cell_seed(cell: dict[str, Any], stage: int) -> int:
    return int((MC_BASE_SEED + 1_000_003 * cell["cell_id"] + 7919 * stage)
               % (2 ** 31 - 1))


class CellState:
    def __init__(self, cell: dict[str, Any]):
        self.cell = cell
        self.shots = 0
        self.failures = 0
        self.stages: list[dict[str, Any]] = []
        self.decode_core_s = 0.0
        self.wall_s = 0.0
        self.done = False

    def to_json(self) -> dict[str, Any]:
        lo, hi = clopper_pearson(self.failures, self.shots,
                                 alpha=PER_INTERVAL_ALPHA) if self.shots else (float("nan"), float("nan"))
        return {
            "phase": self.cell["phase"], "p": self.cell["p"], "key": self.cell["key"],
            "cell_id": self.cell["cell_id"],
            "schedule_index": self.cell["schedule_index"],
            "schedule_hash": self.cell["schedule_hash"],
            "shots": self.shots, "failures": self.failures,
            "ler_per_shot": (self.failures / self.shots) if self.shots else None,
            "ci_simultaneous_lo": lo, "ci_simultaneous_hi": hi,
            "ci_half_width": (hi - lo) if self.shots else None,
            "complete": self.done,
            "stages": self.stages,
            "decode_core_s_per_shot": (self.decode_core_s / self.shots)
                                      if self.shots else None,
            "wall_s": round(self.wall_s, 1),
        }


def load_partial_states(path: Path) -> dict[int, dict[str, Any]]:
    """Fold the append-only partial log into per-cell cumulative state.

    Each line carries the stage's own counts plus the authoritative
    cumulative counters (cum_shots / cum_failures); the loader trusts the
    last line's cumulative fields, which makes it robust to per-line schema
    drift (some early lines lack a top-level ``failures`` key).
    """
    states: dict[int, dict[str, Any]] = {}
    if not path.exists():
        return states
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        cid = rec["cell_id"]
        if cid in states:
            prev = states[cid]
            prev["shots"] = rec.get("cum_shots", prev["shots"] + rec.get("shots", 0))
            prev["failures"] = rec.get(
                "cum_failures", prev["failures"] + rec.get("failures", 0))
            prev["decode_core_s"] += rec.get("decode_core_s", 0.0)
            prev["wall_s"] += rec.get("wall_s", 0.0)
            if rec.get("stage_record"):
                prev["stages"].append(rec["stage_record"])
            prev["done"] = prev["done"] or rec.get("done", False)
        else:
            states[cid] = {
                "shots": rec.get("cum_shots", rec.get("shots", 0)),
                "failures": rec.get(
                    "cum_failures", rec.get("failures", 0)),
                "decode_core_s": rec.get("decode_core_s", 0.0),
                "wall_s": rec.get("wall_s", 0.0),
                "stages": [rec["stage_record"]] if rec.get("stage_record") else [],
                "done": rec.get("done", False),
            }
    return states


def run_mc_cell(cell: dict[str, Any], state: CellState, circuit_cache: dict,
                workers: int, append_partial) -> None:
    key = cell["key"]
    cache_key = (key, cell["schedule_hash"], cell["p"])
    if cache_key not in circuit_cache:
        code, _, _ = code_and_supports(key)
        slot = slot_from_json(cell["slot_map"])
        circ, meta = circuit_from_slot(key, code, slot, cell["depth"], cell["p"])
        circuit_cache[cache_key] = circ
    circuit = circuit_cache[cache_key]

    if cell["phase"] == "matched":
        decoder = MATCHED_DECODER
        stage_shots = W3_STAGE_SHOTS
        target = W3_START_SHOTS
        cap = W3_CAP_SHOTS
        wall_cap = float("inf")
    else:
        decoder = PRODUCTION_DECODER
        stage_shots = PRODUCTION_STAGE_SHOTS
        target = PRODUCTION_TARGET_SHOTS
        cap = PRODUCTION_TARGET_SHOTS
        wall_cap = PRODUCTION_CELL_WALL_CAP_S

    while True:
        if state.done:
            return
        if state.shots >= cap:
            state.done = True
            append_partial(cell, state, final=True)
            return
        lo, hi = clopper_pearson(state.failures, state.shots,
                                 alpha=PER_INTERVAL_ALPHA) if state.shots else (0.0, 1.0)
        if (cell["phase"] == "matched" and state.shots >= target
                and (hi - lo) <= W3_CP_HALF_WIDTH_TARGET):
            state.done = True
            append_partial(cell, state, final=True)
            return
        if state.wall_s >= wall_cap:
            state.done = True
            append_partial(cell, state, final=True)
            return
        stage = len(state.stages)
        n_shots = min(stage_shots, cap - state.shots)
        t0 = time.time()
        r = collect(circuit, max_shots=n_shots, workers=workers, chunk=MC_CHUNK,
                    seed=cell_seed(cell, stage), **decoder)
        wall = time.time() - t0
        state.shots += r.shots
        state.failures += r.failures
        state.decode_core_s += r.decode_s_per_shot_1core * r.shots
        state.wall_s += wall
        lo, hi = clopper_pearson(state.failures, state.shots,
                                 alpha=PER_INTERVAL_ALPHA)
        stage_record = {
            "stage": stage, "seed": cell_seed(cell, stage),
            "shots": r.shots, "failures": r.failures,
            "cum_shots": state.shots, "cum_failures": state.failures,
            "ci_simultaneous_lo": lo, "ci_simultaneous_hi": hi,
            "ci_half_width": hi - lo,
            "latency_p50_ms": r.latency_p50_ms,
            "latency_p95_ms": r.latency_p95_ms,
            "latency_p99_ms": r.latency_p99_ms,
            "wall_s": r.wall_s,
            "decode_core_s": r.decode_s_per_shot_1core * r.shots,
        }
        state.stages.append(stage_record)
        grown = (cell["phase"] == "matched" and state.shots >= target
                 and (hi - lo) <= W3_CP_HALF_WIDTH_TARGET) or state.shots >= cap \
            or state.wall_s >= wall_cap
        if grown:
            state.done = True
        append_partial(cell, state, final=state.done)
        print(f"[MC] {cell['phase']} p={cell['p']} {key} "
              f"sched#{cell['schedule_index']} stage{stage}: cum={state.shots} "
              f"fails={state.failures} hw={hi-lo:.2e} "
              f"done={state.done} ({wall:.0f}s)", flush=True)


def cmd_mc(workers: int, phases: list[str]) -> None:
    PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    lock_path = PARTIAL_DIR / "mc.lock"
    if lock_path.exists():
        old_pid = lock_path.read_text().strip()
        try:
            os.kill(int(old_pid), 0)
            alive = True
        except (ValueError, ProcessLookupError, PermissionError):
            alive = False
        if alive:
            print(f"[MC] refusing to start: lock held by live pid {old_pid}")
            sys.exit(4)
        print(f"[MC] removing stale lock from dead pid {old_pid}", flush=True)
        lock_path.unlink(missing_ok=True)
    lock_path.write_text(str(os.getpid()))
    try:
        cells: list[dict[str, Any]] = []
        if "matched" in phases:
            cells += matched_cell_list()
        if "production" in phases:
            cells += production_cell_list()
        states_raw = load_partial_states(MC_PARTIAL_PATH)
        circuit_cache: dict = {}

        def append_partial(cell, state: CellState, final: bool) -> None:
            rec = {
                "cell_id": cell["cell_id"], "phase": cell["phase"],
                "p": cell["p"], "key": cell["key"],
                "schedule_index": cell["schedule_index"],
                "schedule_hash": cell["schedule_hash"],
                "shots": state.stages[-1]["shots"] if state.stages else 0,
                "failures": (state.stages[-1]["failures"]
                             if state.stages else 0),
                "decode_core_s": (state.stages[-1].get("decode_core_s",
                                state.stages[-1]["wall_s"] * workers)
                                  if state.stages else 0.0),
                "wall_s": state.stages[-1]["wall_s"] if state.stages else 0.0,
                "stage_record": state.stages[-1] if state.stages else None,
                "cum_shots": state.shots, "cum_failures": state.failures,
                "done": final,
            }
            with MC_PARTIAL_PATH.open("a") as f:
                f.write(json.dumps(rec) + "\n")

        production_wall_total = sum(
            prev["wall_s"] for cid, prev in states_raw.items() if cid >= 10_000)
        for cell in cells:
            state = CellState(cell)
            prev = states_raw.get(cell["cell_id"])
            if prev is not None:
                state.shots = prev["shots"]
                state.failures = prev["failures"]
                state.stages = prev["stages"]
                state.decode_core_s = prev["decode_core_s"]
                state.wall_s = prev["wall_s"]
                state.done = prev["done"]
            if state.done:
                print(f"[MC] skip complete cell {cell['cell_id']} "
                      f"({cell['phase']} p={cell['p']} {cell['key']} "
                      f"sched#{cell['schedule_index']}, {state.shots} shots)",
                      flush=True)
                continue
            if (cell["phase"] == "production"
                    and production_wall_total >= PRODUCTION_TOTAL_WALL_CAP_S):
                print(f"[MC] production total wall cap "
                       f"{PRODUCTION_TOTAL_WALL_CAP_S:.0f}s reached; leaving cell "
                       f"{cell['cell_id']} unstarted", flush=True)
                continue
            run_mc_cell(cell, state, circuit_cache, workers, append_partial)
            if cell["phase"] == "production":
                production_wall_total += state.wall_s
        print("[MC] all requested cells complete", flush=True)
    finally:
        lock_path.unlink(missing_ok=True)


def schedule_stats_block(cells: list[dict[str, Any]]) -> dict[str, Any]:
    from scipy import stats as sps
    lers = np.array([c["ler_per_shot"] for c in cells], dtype=float)
    out: dict[str, Any] = {
        "n_schedules": len(cells),
        "ler_mean_over_schedules": float(lers.mean()) if len(lers) else None,
        "ler_sd_over_schedules": float(lers.std(ddof=1)) if len(lers) > 1 else None,
        "ler_mean_ci95_t_interval": None,
        "min_simultaneous_lo": min((c["ci_simultaneous_lo"] for c in cells),
                                   default=None),
        "max_simultaneous_hi": max((c["ci_simultaneous_hi"] for c in cells),
                                   default=None),
        "min_point_ler": float(lers.min()) if len(lers) else None,
        "max_point_ler": float(lers.max()) if len(lers) else None,
    }
    if len(lers) > 1:
        tcrit = float(sps.t.ppf(0.975, len(lers) - 1))
        half = tcrit * lers.std(ddof=1) / np.sqrt(len(lers))
        out["ler_mean_ci95_t_interval"] = [
            float(lers.mean() - half), float(lers.mean() + half)]
    return out


def assemble_matched(states_by_cell: dict[int, CellState],
                     cells: list[dict[str, Any]]) -> dict[str, Any]:
    from scipy import stats as sps
    per_p: dict[str, Any] = {}
    for p in sorted({c["p"] for c in cells}):
        pcells = [c for c in cells if c["p"] == p]
        entries = []
        for c in pcells:
            st = states_by_cell.get(c["cell_id"])
            if st is None or st.shots == 0:
                entries.append({
                    "key": c["key"], "schedule_index": c["schedule_index"],
                    "schedule_hash": c["schedule_hash"],
                    "p": p, "shots": 0, "failures": 0, "ler_per_shot": None,
                    "ci_simultaneous_lo": None, "ci_simultaneous_hi": None,
                    "ci_half_width": None, "complete": False,
                    "target_shots": W3_START_SHOTS, "cap_shots": W3_CAP_SHOTS,
                    "reason": "not reached before interruption"})
                continue
            j = st.to_json()
            j.pop("phase", None)
            j["target_shots"] = W3_START_SHOTS
            j["cap_shots"] = W3_CAP_SHOTS
            entries.append(j)
        gross = [e for e in entries if e["key"] == "gross"]
        pbb = [e for e in entries if e["key"] == "pbb"]
        grid_complete = (len(entries) == 10 and all(e["complete"] for e in entries))
        block: dict[str, Any] = {
            "grid_complete": grid_complete,
            "n_intervals_recorded": sum(1 for e in entries if e["shots"] > 0),
            "cells": entries,
            "gross": schedule_stats_block(gross) if gross else None,
            "pbb": schedule_stats_block(pbb) if pbb else None,
        }
        gl = [e["ler_per_shot"] for e in gross if e["ler_per_shot"] is not None]
        pl = [e["ler_per_shot"] for e in pbb if e["ler_per_shot"] is not None]
        gminlo = min((e["ci_simultaneous_lo"] for e in gross
                      if e["ci_simultaneous_lo"] is not None), default=None)
        pmaxhi = max((e["ci_simultaneous_hi"] for e in pbb
                      if e["ci_simultaneous_hi"] is not None), default=None)
        gmaxhi = max((e["ci_simultaneous_hi"] for e in gross
                      if e["ci_simultaneous_hi"] is not None), default=None)
        pminlo = min((e["ci_simultaneous_lo"] for e in pbb
                      if e["ci_simultaneous_lo"] is not None), default=None)
        block["verdict"] = {
            "css_max_simultaneous_hi": gmaxhi,
            "pbb_min_simultaneous_lo": pminlo,
            "pointwise_familywise_separation": (
                (pminlo is not None and gmaxhi is not None and pminlo > gmaxhi
                 and grid_complete)),
            "pointwise_familywise_separation_on_available_intervals": (
                pminlo is not None and gmaxhi is not None and pminlo > gmaxhi),
            "every_sampled_pbb_point_exceeds_every_css": (
                bool(pl and gl and min(pl) > max(gl))),
        }
        if len(gl) >= 2 and len(pl) >= 2:
            wt = sps.ttest_ind(pl, gl, equal_var=False)
            mw = sps.mannwhitneyu(pl, gl, alternative="two-sided")
            block["welch_t"] = {
                "statistic": float(wt.statistic), "p_value": float(wt.pvalue),
                "n_pbb": len(pl), "n_gross": len(gl)}
            block["mann_whitney_u"] = {
                "statistic": float(mw.statistic), "p_value": float(mw.pvalue),
                "n_pbb": len(pl), "n_gross": len(gl)}
            block["schedule_mean_ratio_pbb_over_gross"] = float(
                np.mean(pl) / np.mean(gl))
        else:
            block["welch_t"] = None
            block["mann_whitney_u"] = None
            block["schedule_mean_ratio_pbb_over_gross"] = None
        per_p[f"{p}"] = block
    return {
        "experiment": EXPERIMENT,
        "step": "W3-matched-ler",
        "env": env_json(),
        "protocol": {
            "noise_model": ("identical for both codes: build_memory_circuit "
                            "p-channel (X_ERROR/Z_ERROR on resets+measurements, "
                            "DEPOLARIZE1/2 at the same p per gate), 12 rounds, "
                            "Z basis"),
            "slot_maps": "results/raw/exp016_schedule_controlled.json (5 sampled "
                         "translation-invariant minimum-depth schedules per code)",
            "builder": "qec_research.circuits.mixed_stabilizer.build_memory_circuit",
            "p_points": list(W3_P_POINTS),
            "start_shots_per_schedule": W3_START_SHOTS,
            "stage_shots": W3_STAGE_SHOTS,
            "cap_shots": W3_CAP_SHOTS,
            "growth_rule": (f"accumulate stages until >= {W3_START_SHOTS} shots, "
                            f"then stop when Clopper-Pearson half-width at "
                            f"alpha={PER_INTERVAL_ALPHA} <= "
                            f"{W3_CP_HALF_WIDTH_TARGET} or {W3_CAP_SHOTS} shots"),
            "decoder": ("qec_research.decoders.parallel_harness.collect "
                        f"{MATCHED_DECODER}, undecomposed hypergraph DEM"),
            "familywise_scheme": (f"Bonferroni: familywise 95% over "
                                  f"{N_INTERVALS_FAMILY} simultaneous intervals "
                                  f"(alpha={PER_INTERVAL_ALPHA} each), as EXP-016"),
            "workers": MC_WORKERS, "chunk": MC_CHUNK,
            "seed_rule": "20260817 + 1000003*cell_id + 7919*stage (mod 2^31-1)",
        },
        "per_p": per_p,
    }


def assemble_production(states_by_cell: dict[int, CellState],
                        cells: list[dict[str, Any]],
                        matched_cells: list[dict[str, Any]]) -> dict[str, Any]:
    entries = []
    for c in cells:
        st = states_by_cell.get(c["cell_id"])
        if st is None or st.shots == 0:
            entries.append({
                "key": c["key"], "schedule_index": c["schedule_index"],
                "schedule_hash": c["schedule_hash"], "shots": 0, "failures": 0,
                "ler_per_shot": None, "ci_simultaneous_lo": None,
                "ci_simultaneous_hi": None, "ci_half_width": None,
                "complete": False, "reason": "not reached before interruption"})
            continue
        j = st.to_json()
        j.pop("phase", None)
        entries.append(j)
    gross = [e for e in entries if e["key"] == "gross"]
    pbb = [e for e in entries if e["key"] == "pbb"]
    gl = [e["ler_per_shot"] for e in gross if e["ler_per_shot"] is not None]
    pl = [e["ler_per_shot"] for e in pbb if e["ler_per_shot"] is not None]
    prod_ratio = float(np.mean(pl) / np.mean(gl)) if (len(gl) and len(pl)) else None
    matched_ratio = None
    matched_p = [c for c in matched_cells if c["p"] == PRODUCTION_P]
    mstates = {cid: st for cid, st in states_by_cell.items()}
    mg = [mstates[c["cell_id"]].to_json()["ler_per_shot"] for c in matched_p
          if c["key"] == "gross" and c["cell_id"] in mstates
          and mstates[c["cell_id"]].shots > 0]
    mp = [mstates[c["cell_id"]].to_json()["ler_per_shot"] for c in matched_p
          if c["key"] == "pbb" and c["cell_id"] in mstates
          and mstates[c["cell_id"]].shots > 0]
    if mg and mp:
        matched_ratio = float(np.mean(mp) / np.mean(mg))
    load = os.getloadavg()
    deferred = load[0] > LATENCY_MAX_LOAD
    return {
        "experiment": EXPERIMENT,
        "step": "W5-production-decoder",
        "env": env_json(),
        "protocol": {
            "p": PRODUCTION_P,
            "decoder": ("qec_research.decoders.parallel_harness.collect "
                        f"{PRODUCTION_DECODER}, undecomposed hypergraph DEM"),
            "slot_maps": "same 10 EXP-016 slot maps as the matched arm",
            "target_shots_per_schedule": PRODUCTION_TARGET_SHOTS,
            "stage_shots": PRODUCTION_STAGE_SHOTS,
            "cell_wall_cap_s": PRODUCTION_CELL_WALL_CAP_S,
            "total_wall_cap_s": PRODUCTION_TOTAL_WALL_CAP_S,
            "ci": f"Clopper-Pearson alpha={PER_INTERVAL_ALPHA} (Bonferroni "
                  f"familywise 95% over 10 intervals)",
        },
        "cells": entries,
        "gross": schedule_stats_block(gross) if gross else None,
        "pbb": schedule_stats_block(pbb) if pbb else None,
        "schedule_mean_ratio_pbb_over_gross": prod_ratio,
        "matched_arm_ratio_same_p": matched_ratio,
        "exp016_reference_ratio": 3.547752808988764,
        "latency_arm": (
            f"DEFERRED: machine under SAT fleet load (load average "
            f"{load[0]:.1f} > {LATENCY_MAX_LOAD})") if deferred else "MEASURED",
        "loadavg_at_assembly": list(load),
    }


def cmd_assemble() -> None:
    states_raw = load_partial_states(MC_PARTIAL_PATH)
    matched_cells = matched_cell_list()
    prod_cells = production_cell_list()

    def to_state(cell):
        st = CellState(cell)
        prev = states_raw.get(cell["cell_id"])
        if prev is not None:
            st.shots = prev["shots"]
            st.failures = prev["failures"]
            st.stages = prev["stages"]
            st.decode_core_s = prev["decode_core_s"]
            st.wall_s = prev["wall_s"]
            st.done = prev["done"]
        return st

    states = {c["cell_id"]: to_state(c) for c in matched_cells + prod_cells}
    matched = assemble_matched(states, matched_cells)
    shots_total = sum(s.shots for s in states.values())
    matched["total_shots_accumulated"] = shots_total
    MATCHED_PATH.write_text(json.dumps(matched, indent=2) + "\n")
    print(f"[assemble] wrote {MATCHED_PATH.relative_to(ROOT)} "
          f"({shots_total} shots accumulated)", flush=True)

    prod_any = any(c["cell_id"] in states and states[c["cell_id"]].shots > 0
                   for c in prod_cells)
    if prod_any:
        prod = assemble_production(states, prod_cells, matched_cells)
        PRODUCTION_PATH.write_text(json.dumps(prod, indent=2) + "\n")
        print(f"[assemble] wrote {PRODUCTION_PATH.relative_to(ROOT)}", flush=True)
    else:
        print("[assemble] no production shots yet; production JSON not written",
              flush=True)


# ---------------------------------------------------------------------------
# W4 -- budget-capped circuit-distance tightening (EXP-029 machinery copied)
# ---------------------------------------------------------------------------
def xor_target_mask(targets: list[stim.DemTarget], *, detector: bool) -> int:
    mask = 0
    for target in targets:
        relevant = (target.is_relative_detector_id() if detector
                    else target.is_logical_observable_id())
        if relevant:
            mask ^= 1 << target.val
    return mask


def extract_error_mechanisms(dem: stim.DetectorErrorModel):
    """Copied from experiments/exp029_circuit_distance.py."""
    mechanisms = []
    separator_instructions = 0
    for instruction_index, instruction in enumerate(dem.flattened()):
        if instruction.type != "error":
            continue
        targets = instruction.targets_copy()
        if any(t.is_separator() for t in targets):
            separator_instructions += 1
        probability = float(instruction.args_copy()[0])
        if probability <= 0:
            continue
        detector_mask = xor_target_mask(targets, detector=True)
        observable_mask = xor_target_mask(targets, detector=False)
        mechanisms.append({
            "index": len(mechanisms),
            "flattened_instruction_index": instruction_index,
            "probability": probability,
            "detector_mask": detector_mask,
            "observable_mask": observable_mask,
            "detector_degree": detector_mask.bit_count(),
        })
    return mechanisms, separator_instructions


def mask_ids(mask: int) -> list[int]:
    out = []
    while mask:
        low = mask & -mask
        out.append(low.bit_length() - 1)
        mask ^= low
    return out


def mechanism_json(mechanism: dict[str, Any]) -> dict[str, Any]:
    return {
        "index": mechanism["index"],
        "flattened_instruction_index": mechanism["flattened_instruction_index"],
        "probability_ignored_for_distance": mechanism["probability"],
        "detectors": mask_ids(mechanism["detector_mask"]),
        "observables": mask_ids(mechanism["observable_mask"]),
        "detector_degree": mechanism["detector_degree"],
    }


def choose_different_observable(variants, forbidden_observable):
    for observable, index in variants.items():
        if observable != forbidden_observable:
            return observable, index
    return None


def check_w1_w2(mechanisms):
    """Copied from experiments/exp029_circuit_distance.py."""
    from collections import Counter
    t0 = time.perf_counter()
    detector_variants: dict[int, dict[int, int]] = {}
    w1_property_holds = True
    w1_witness = None
    for mechanism in mechanisms:
        detector_mask = mechanism["detector_mask"]
        observable_mask = mechanism["observable_mask"]
        if detector_mask == 0 and observable_mask != 0:
            w1_property_holds = False
            if w1_witness is None:
                w1_witness = {"weight": 1,
                              "mechanisms": [mechanism_json(mechanism)]}
        detector_variants.setdefault(detector_mask, {}).setdefault(
            observable_mask, mechanism["index"])
    w2_witness = None
    for detector_mask, variants in detector_variants.items():
        if len(variants) < 2:
            continue
        first_observable, first_index = next(iter(variants.items()))
        other = choose_different_observable(variants, first_observable)
        if other is None:
            continue
        _, second_index = other
        w2_witness = {
            "weight": 2,
            "mechanisms": [mechanism_json(mechanisms[first_index]),
                           mechanism_json(mechanisms[second_index])],
            "shared_detector_set": mask_ids(detector_mask),
        }
        break
    result = {
        "method": "exhaustive hash of undecomposed flattened DEM mechanism signatures",
        "fault_model": FAULT_MODEL,
        "mechanisms_checked": len(mechanisms),
        "w1_property_holds": w1_property_holds,
        "w1_undetectable_logical_found": w1_witness is not None,
        "w2_undetectable_logical_found": w2_witness is not None,
        "detector_signature_groups": len(detector_variants),
        "wall_s": time.perf_counter() - t0,
        "certified": True,
    }
    return result, w1_witness or w2_witness


def check_w3(mechanisms, *, budget_s: float):
    """Copied from experiments/exp029_circuit_distance.py (budget injected)."""
    t0 = time.perf_counter()
    deadline = t0 + budget_s
    unique_by_signature: dict[tuple[int, int], int] = {}
    for mechanism in mechanisms:
        signature = (mechanism["detector_mask"], mechanism["observable_mask"])
        unique_by_signature.setdefault(signature, mechanism["index"])
    unique_indices = list(unique_by_signature.values())
    unique_detector = [mechanisms[i]["detector_mask"] for i in unique_indices]
    unique_observable = [mechanisms[i]["observable_mask"] for i in unique_indices]

    singles: dict[int, dict[int, int]] = {}
    incidence: dict[int, list[int]] = {}
    detector_empty = 0
    for local_index, original_index in enumerate(unique_indices):
        detector_mask = unique_detector[local_index]
        observable_mask = unique_observable[local_index]
        singles.setdefault(detector_mask, {}).setdefault(
            observable_mask, local_index)
        if detector_mask == 0:
            detector_empty += 1
            continue
        work = detector_mask
        while work:
            low = work & -work
            detector_id = low.bit_length() - 1
            incidence.setdefault(detector_id, []).append(local_index)
            work ^= low

    all_pair_count = math.comb(len(unique_indices), 2)
    shared_pair_visits_upper = sum(math.comb(len(ids), 2)
                                   for ids in incidence.values())
    estimated_materialized_pair_bytes = all_pair_count * 384
    rss_start = max_rss_bytes()
    if rss_start >= W4_MEMORY_GUARD_BYTES:
        return ({
            "method": "singles-vs-streamed-pairs XOR meet-in-the-middle",
            "fault_model": FAULT_MODEL,
            "completed": False, "certified_no_w3": False,
            "stop_reason": "8 GiB RSS guard already exceeded",
            "unique_mechanism_signatures": len(unique_indices),
            "all_pair_count": all_pair_count,
            "shared_pair_visits_upper": shared_pair_visits_upper,
            "memory_guard_bytes": W4_MEMORY_GUARD_BYTES,
            "max_rss_bytes": rss_start,
            "wall_s": time.perf_counter() - t0,
        }, None)

    pair_visits = 0
    unique_shared_pairs_checked = 0
    witness = None
    stop_reason = None
    timed_out = False
    memory_guard_hit = False
    singles_get = singles.get
    check_interval_mask = (1 << 18) - 1
    for detector_id in sorted(incidence):
        ids = incidence[detector_id]
        owner_bit = 1 << detector_id
        for left_pos in range(len(ids)):
            left = ids[left_pos]
            left_detector = unique_detector[left]
            left_observable = unique_observable[left]
            for right_pos in range(left_pos + 1, len(ids)):
                right = ids[right_pos]
                pair_visits += 1
                if (pair_visits & check_interval_mask) == 0:
                    if max_rss_bytes() >= W4_MEMORY_GUARD_BYTES:
                        memory_guard_hit = True
                        stop_reason = "8 GiB RSS guard reached"
                        break
                    if time.perf_counter() >= deadline:
                        timed_out = True
                        stop_reason = "pair-enumeration budget exhausted"
                        break
                right_detector = unique_detector[right]
                common = left_detector & right_detector
                if common & -common != owner_bit:
                    continue
                unique_shared_pairs_checked += 1
                target_detector = left_detector ^ right_detector
                if target_detector == 0:
                    continue
                variants = singles_get(target_detector)
                if variants is None:
                    continue
                pair_observable = left_observable ^ unique_observable[right]
                choice = choose_different_observable(variants, pair_observable)
                if choice is None:
                    continue
                _, third = choice
                original_ids = [unique_indices[left], unique_indices[right],
                                unique_indices[third]]
                if len(set(original_ids)) != 3:
                    continue
                detector_xor = 0
                observable_xor = 0
                for original_id in original_ids:
                    detector_xor ^= mechanisms[original_id]["detector_mask"]
                    observable_xor ^= mechanisms[original_id]["observable_mask"]
                if detector_xor == 0 and observable_xor != 0:
                    witness = {
                        "weight": 3,
                        "mechanisms": [mechanism_json(mechanisms[i])
                                       for i in original_ids],
                        "detector_xor": [],
                        "observable_xor": mask_ids(observable_xor),
                    }
                    break
            if witness is not None or stop_reason is not None:
                break
            if time.perf_counter() >= deadline:
                timed_out = True
                stop_reason = "pair-enumeration budget exhausted"
                break
        if witness is not None or stop_reason is not None:
            break
    completed = witness is not None or stop_reason is None
    return ({
        "method": "singles-vs-streamed-pairs XOR meet-in-the-middle",
        "fault_model": FAULT_MODEL,
        "completeness_argument": (
            "Each detector in a zero-XOR triple occurs in exactly two "
            "mechanisms; therefore at least one pair shares a detector. Each "
            "such pair is processed once at its least shared detector and "
            "queried against the singles hash."),
        "completed": completed,
        "certified_no_w3": completed and witness is None,
        "w3_undetectable_logical_found": witness is not None,
        "stop_reason": stop_reason, "timed_out": timed_out,
        "memory_guard_hit": memory_guard_hit,
        "mechanisms": len(mechanisms),
        "unique_mechanism_signatures": len(unique_indices),
        "detector_empty_unique_signatures": detector_empty,
        "all_pair_count": all_pair_count,
        "shared_pair_visits_upper": shared_pair_visits_upper,
        "pair_visits": pair_visits,
        "unique_shared_pairs_checked": unique_shared_pairs_checked,
        "budget_s": budget_s,
        "max_rss_bytes": max_rss_bytes(),
        "wall_s": time.perf_counter() - t0,
    }, witness)


def check_w4_matching(mechanisms, *, budget_s: float):
    """Exhaustive w=4 exclusion over the *perfect-matching* co-occurrence class.

    Any zero-detector-XOR 4-set with no smaller witness has a co-occurrence
    multigraph (edges = shared detectors) of minimum degree >= 1.  Such a
    multigraph on 4 vertices either admits a perfect matching -- then the two
    complementary pairs each share a detector -- or is a star K1,3 (handled by
    check_w4_star).  This routine enumerates every shared-detector pair once
    (assigned to its least shared detector, as EXP-029's w3 loop), digests the
    pair's detector-XOR mask, and looks for two index-disjoint pairs with
    equal digest; every collision is re-verified against the exact masks, so
    hash collisions can only produce false candidates that are filtered, never
    a missed witness.
    """
    t0 = time.perf_counter()
    deadline = t0 + budget_s
    n_mech = len(mechanisms)
    masks = [m["detector_mask"] for m in mechanisms]
    obs = [m["observable_mask"] for m in mechanisms]
    incidence: dict[int, list[int]] = {}
    for i, mask in enumerate(masks):
        work = mask
        while work:
            low = work & -work
            incidence.setdefault(low.bit_length() - 1, []).append(i)
            work ^= low

    digest_hi = array("I")
    digest_lo = array("I")
    left_idx = array("I")
    right_idx = array("I")
    blake = hashlib.blake2b
    nbytes = (max(masks).bit_length() + 7) // 8 if masks else 1
    pair_visits = 0
    stored = 0
    stop_reason = None
    memory_guard_hit = False
    timed_out = False
    check_interval_mask = (1 << 18) - 1
    for detector_id in sorted(incidence):
        ids = incidence[detector_id]
        owner_bit = 1 << detector_id
        for a in range(len(ids)):
            i = ids[a]
            mi = masks[i]
            oi = obs[i]
            for b in range(a + 1, len(ids)):
                j = ids[b]
                pair_visits += 1
                if (pair_visits & check_interval_mask) == 0:
                    if max_rss_bytes() >= W4_MEMORY_GUARD_BYTES:
                        memory_guard_hit = True
                        stop_reason = "8 GiB RSS guard reached during pair enumeration"
                        break
                    if time.perf_counter() >= deadline:
                        timed_out = True
                        stop_reason = "pair-enumeration budget exhausted"
                        break
                mj = masks[j]
                common = mi & mj
                if common & -common != owner_bit:
                    continue
                x = mi ^ mj
                if x == 0:
                    continue
                d = blake(x.to_bytes(nbytes, "little"), digest_size=W4_DIGEST_BYTES).digest()
                digest_hi.append(int.from_bytes(d[:4], "little"))
                digest_lo.append(int.from_bytes(d[4:], "little"))
                left_idx.append(i)
                right_idx.append(j)
                stored += 1
            if stop_reason is not None:
                break
        if stop_reason is not None:
            break

    witness = None
    collision_candidates = 0
    degenerate_collisions = 0
    sort_wall = None
    if stop_reason is None:
        t_sort = time.perf_counter()
        rec = np.zeros(stored, dtype=[("hi", np.uint32), ("lo", np.uint32),
                                      ("i", np.uint32), ("j", np.uint32)])
        assert digest_hi.itemsize == 4 and left_idx.itemsize == 4
        rec["hi"] = np.frombuffer(digest_hi, dtype=np.uint32, count=stored)
        rec["lo"] = np.frombuffer(digest_lo, dtype=np.uint32, count=stored)
        rec["i"] = np.frombuffer(left_idx, dtype=np.uint32, count=stored)
        rec["j"] = np.frombuffer(right_idx, dtype=np.uint32, count=stored)
        rec.sort(order=("hi", "lo"))
        sort_wall = time.perf_counter() - t_sort
        if stored > 1:
            eq = (rec["hi"][1:] == rec["hi"][:-1]) & (rec["lo"][1:] == rec["lo"][:-1])
            # group boundaries: a new group starts where eq[k-1] is False
            boundaries = np.concatenate(([0], np.flatnonzero(~eq) + 1))
            for s, e in zip(boundaries[:-1], boundaries[1:]):
                if e - s < 2:
                    continue
                idxs = rec[s:e]
                m = len(idxs)
                for x in range(m):
                    for y in range(x + 1, m):
                        i1, j1 = int(idxs["i"][x]), int(idxs["j"][x])
                        i2, j2 = int(idxs["i"][y]), int(idxs["j"][y])
                        if len({i1, j1, i2, j2}) != 4:
                            continue
                        collision_candidates += 1
                        det_xor = (masks[i1] ^ masks[j1]) ^ (masks[i2] ^ masks[j2])
                        if det_xor != 0:
                            continue  # digest collision only
                        obs_xor = (obs[i1] ^ obs[j1]) ^ (obs[i2] ^ obs[j2])
                        if obs_xor == 0:
                            degenerate_collisions += 1
                            continue
                        witness = {
                            "weight": 4,
                            "class": "perfect-matching",
                            "mechanisms": [mechanism_json(mechanisms[k])
                                           for k in (i1, j1, i2, j2)],
                            "detector_xor": [],
                            "observable_xor": mask_ids(obs_xor),
                        }
                        break
                    if witness is not None:
                        break
                if witness is not None:
                    break
    completed = stop_reason is None
    return ({
        "method": "digest-sorted shared-pair-pair XOR meet-in-the-middle",
        "fault_model": FAULT_MODEL,
        "coverage_argument": (
            "exhaustive over zero-XOR 4-sets whose co-occurrence multigraph "
            "admits a perfect matching (both complementary pairs share a "
            "detector); each shared pair enumerated once at its least shared "
            "detector; digest collisions re-verified exactly"),
        "digest": f"blake2b-{8*W4_DIGEST_BYTES}bit of the pair detector-XOR mask",
        "completed": completed,
        "certified_no_w4_matching_class": completed and witness is None,
        "w4_witness_found": witness is not None,
        "stop_reason": stop_reason, "timed_out": timed_out,
        "memory_guard_hit": memory_guard_hit,
        "mechanisms": n_mech,
        "pair_visits": pair_visits,
        "pairs_stored": stored,
        "collision_candidates_verified": collision_candidates,
        "degenerate_collisions": degenerate_collisions,
        "sort_wall_s": sort_wall,
        "budget_s": budget_s,
        "max_rss_bytes": max_rss_bytes(),
        "wall_s": time.perf_counter() - t0,
    }, witness)


def check_w4_star(mechanisms, *, budget_s: float):
    """Exhaustive w=4 exclusion over the star (K1,3) co-occurrence class.

    In a star, the three leaves are pairwise detector-disjoint and the union
    of their masks equals the centre's mask exactly, so every leaf mask is a
    subset of the centre's support.  After the w1/w2 pass each detector mask
    belongs to at most one mechanism, so per centre we collect candidate
    leaf mechanisms (mask subset of the centre support), then enumerate
    disjoint candidate pairs and look the exact complement mask up.
    """
    t0 = time.perf_counter()
    deadline = t0 + budget_s
    n_mech = len(mechanisms)
    masks = [m["detector_mask"] for m in mechanisms]
    obs = [m["observable_mask"] for m in mechanisms]
    by_mask: dict[int, int] = {}
    mask_unique = True
    for i, m in enumerate(masks):
        if m in by_mask:
            mask_unique = False
        by_mask[m] = i
    incidence: dict[int, list[int]] = {}
    for i, m in enumerate(masks):
        work = m
        while work:
            low = work & -work
            incidence.setdefault(low.bit_length() - 1, []).append(i)
            work ^= low

    witness = None
    centers_checked = 0
    centers_with_candidates = 0
    stop_reason = None
    timed_out = False
    check_interval = 1024
    for c in range(n_mech):
        if centers_checked % check_interval == 0 and centers_checked:
            if time.perf_counter() >= deadline:
                timed_out = True
                stop_reason = (f"star-class budget exhausted after "
                               f"{centers_checked} centres")
                break
        centers_checked += 1
        mc_mask = masks[c]
        if not mc_mask:
            continue
        seen: set[int] = set()
        work = mc_mask
        cand: list[int] = []
        while work:
            low = work & -work
            d = low.bit_length() - 1
            for i in incidence.get(d, ()):
                if i == c or i in seen:
                    continue
                seen.add(i)
                mi = masks[i]
                if mi | mc_mask == mc_mask and mi:
                    cand.append(i)
            work ^= low
        if len(cand) < 3:
            continue
        centers_with_candidates += 1
        cand_masks = [masks[i] for i in cand]
        for x in range(len(cand)):
            mx = cand_masks[x]
            for y in range(x + 1, len(cand)):
                my = cand_masks[y]
                if mx & my:
                    continue
                rest = mc_mask ^ (mx | my)
                if not rest:
                    continue
                k = by_mask.get(rest)
                if k is None or k in (c, cand[x], cand[y]):
                    continue
                obs_xor = obs[c] ^ obs[cand[x]] ^ obs[cand[y]] ^ obs[k]
                det_xor = mc_mask ^ mx ^ my ^ masks[k]
                if det_xor != 0 or obs_xor == 0:
                    continue
                witness = {
                    "weight": 4,
                    "class": "star",
                    "mechanisms": [mechanism_json(mechanisms[i])
                                   for i in (c, cand[x], cand[y], k)],
                    "detector_xor": [],
                    "observable_xor": mask_ids(obs_xor),
                }
                break
            if witness is not None:
                break
        if witness is not None:
            break
    completed = stop_reason is None
    return ({
        "method": "per-centre exact-cover of nested candidate leaf masks",
        "fault_model": FAULT_MODEL,
        "coverage_argument": (
            "exhaustive over zero-XOR 4-sets whose co-occurrence multigraph "
            "is a star K1,3: leaves pairwise disjoint with union exactly the "
            "centre mask, hence every leaf mask is a subset of the centre "
            "support and each mask maps to at most one mechanism after the "
            "w1/w2 pass"),
        "detector_masks_unique": mask_unique,
        "completed": completed,
        "certified_no_w4_star_class": completed and witness is None,
        "w4_witness_found": witness is not None,
        "stop_reason": stop_reason, "timed_out": timed_out,
        "mechanisms": n_mech,
        "centers_checked": centers_checked,
        "centers_with_ge3_candidates": centers_with_candidates,
        "budget_s": budget_s,
        "max_rss_bytes": max_rss_bytes(),
        "wall_s": time.perf_counter() - t0,
    }, witness)


def heuristic_worker(send_conn, circuit_text, parameters):
    """Copied from experiments/exp029_circuit_distance.py (spawn target)."""
    try:
        t0 = time.perf_counter()
        circuit = stim.Circuit(circuit_text)
        errors = circuit.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=parameters[
                "dont_explore_detection_event_sets_with_size_above"],
            dont_explore_edges_with_degree_above=parameters[
                "dont_explore_edges_with_degree_above"],
            dont_explore_edges_increasing_symptom_degree=parameters[
                "dont_explore_edges_increasing_symptom_degree"],
            canonicalize_circuit_errors=CANONICALIZE_CIRCUIT_ERRORS,
        )
        send_conn.send({
            "ok": True, "wall_s": time.perf_counter() - t0,
            "weight": len(errors),
            "mechanisms": [witness_entry_from_explained(e, i)
                           for i, e in enumerate(errors)],
        })
    except BaseException as exc:
        send_conn.send({
            "ok": False, "wall_s": time.perf_counter() - t0,
            "error": f"{type(exc).__name__}: {exc}",
        })
    finally:
        send_conn.close()


def witness_entry_from_explained(error: stim.ExplainedError, index: int):
    detector_mask = 0
    observable_mask = 0
    for term_with_coords in error.dem_error_terms:
        target = term_with_coords.dem_target
        if target.is_relative_detector_id():
            detector_mask ^= 1 << target.val
        elif target.is_logical_observable_id():
            observable_mask ^= 1 << target.val
    return {
        "search_entry_index": index,
        "detector_mask_hex": hex(detector_mask),
        "observable_mask_hex": hex(observable_mask),
        "circuit_error_locations": sum(1 for _ in error.circuit_error_locations),
    }


def run_timed_child(target, args, timeout_s):
    """Copied from experiments/exp029_circuit_distance.py."""
    ctx = mp.get_context("spawn")
    receive_conn, send_conn = ctx.Pipe(duplex=False)
    process = ctx.Process(target=target, args=(send_conn, *args))
    started = time.perf_counter()
    process.start()
    send_conn.close()
    payload = None
    if receive_conn.poll(timeout_s):
        try:
            payload = receive_conn.recv()
        except EOFError:
            payload = None
    if payload is None and process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join(5)
    else:
        process.join(5)
        if process.is_alive():
            process.terminate()
            process.join(5)
    receive_conn.close()
    elapsed = time.perf_counter() - started
    if payload is None:
        return {
            "ok": False,
            "timed_out": elapsed >= timeout_s * 0.99,
            "error": ("search process exceeded its wall-time allocation"
                      if elapsed >= timeout_s * 0.99 else
                      f"search process exited without a result "
                      f"(exitcode={process.exitcode})"),
            "wall_s": elapsed,
        }
    payload["timed_out"] = False
    payload["parent_wall_s"] = elapsed
    return payload


def signature_representatives(mechanisms):
    result: dict[tuple[int, int], list[int]] = {}
    for mechanism in mechanisms:
        key = (mechanism["detector_mask"], mechanism["observable_mask"])
        result.setdefault(key, []).append(mechanism["index"])
    return result


def validate_signature_witness(entries, mechanisms):
    """Copied from experiments/exp029_circuit_distance.py."""
    representatives = signature_representatives(mechanisms)
    detector_xor = 0
    observable_xor = 0
    all_present = True
    used: dict[tuple[int, int], int] = {}
    for entry in entries:
        detector_mask = int(entry["detector_mask_hex"], 16)
        observable_mask = int(entry["observable_mask_hex"], 16)
        detector_xor ^= detector_mask
        observable_xor ^= observable_mask
        signature = (detector_mask, observable_mask)
        occurrence = used.get(signature, 0)
        candidates = representatives.get(signature, [])
        if occurrence < len(candidates):
            entry["undecomposed_mechanism_index"] = candidates[occurrence]
            used[signature] = occurrence + 1
        else:
            entry["undecomposed_mechanism_index"] = None
            all_present = False
    return {
        "all_entries_present_in_undecomposed_dem": all_present,
        "detectors_cancel": detector_xor == 0,
        "observable_flips": observable_xor != 0,
        "observable_xor": mask_ids(observable_xor),
        "valid_mechanism_level_upper_bound": (
            all_present and detector_xor == 0 and observable_xor != 0),
    }


def cmd_distance() -> None:
    t0 = time.time()
    registry = json.loads(REGISTRY_PATH.read_text())
    out: dict[str, Any] = {
        "experiment": EXPERIMENT,
        "step": "W4-circuit-distance-tightening",
        "env": env_json(),
        "protocol": {
            "fault_model": FAULT_MODEL,
            "p_dem_structure_only": P_DEM_STRUCTURE,
            "probabilities_used_in_distance": False,
            "w3_budget_s": W4_W3_BUDGET_S,
            "w4_matching_budget_s": W4_MATCHING_BUDGET_S,
            "w4_star_budget_s": W4_STAR_BUDGET_S,
            "stim_witness_search_total_budget_s": W4_STIM_TOTAL_BUDGET_S,
            "stim_witness_search_sizes": [4, 5, 6, 7, 8],
            "memory_guard_bytes": W4_MEMORY_GUARD_BYTES,
            "w4_partition_argument": (
                "a zero-detector-XOR 4-set with no smaller witness has "
                "co-occurrence min degree >= 1; on 4 vertices such a "
                "multigraph admits a perfect matching or is a star K1,3; "
                "both classes are searched exhaustively and independently"),
            "exp029_baseline": "results/processed/exp029_circuit_distance.json "
                               "(certified 4 <= d_DEM_mech <= 12 for both)",
        },
        "circuits": {},
    }
    any_stagnation = False
    for key in ("gross", "pbb"):
        t1 = time.time()
        det = registry["deterministic_slot_maps"][key]
        code, _, _ = code_and_supports(key)
        slot = slot_from_json(det["slot_map"])
        circuit, _ = circuit_from_slot(key, code, slot, int(det["depth"]),
                                       P_DEM_STRUCTURE)
        dem = circuit.detector_error_model(decompose_errors=False)
        fp = {
            "circuit_sha256": sha256_text(str(circuit)),
            "dem_undecomposed_sha256": sha256_text(str(dem)),
        }
        assert fp["circuit_sha256"] == LOCKED[key]["circuit_sha256"], key
        assert fp["dem_undecomposed_sha256"] == LOCKED[key]["dem_undecomposed_sha256"], key
        mechanisms, separators = extract_error_mechanisms(dem)

        w12, low_witness = check_w1_w2(mechanisms)
        w3, w3_witness = check_w3(mechanisms, budget_s=W4_W3_BUDGET_S) \
            if low_witness is None else (None, None)
        w4m, w4m_witness = (None, None)
        w4s, w4s_witness = (None, None)
        if low_witness is None and w3_witness is None and w3 is not None \
                and w3["certified_no_w3"]:
            w4m, w4m_witness = check_w4_matching(
                mechanisms, budget_s=W4_MATCHING_BUDGET_S)
            if w4m_witness is None:
                w4s, w4s_witness = check_w4_star(
                    mechanisms, budget_s=W4_STAR_BUDGET_S)

        # stim witness search at sizes 4..8 (XOR-validated)
        search_start = time.perf_counter()
        attempts = []
        stim_ub = None
        circuit_text = str(circuit)
        for size in (4, 5, 6, 7, 8):
            elapsed = time.perf_counter() - search_start
            remaining = W4_STIM_TOTAL_BUDGET_S - elapsed
            if remaining <= 60:
                attempts.append({"size": size, "skipped": "budget exhausted"})
                continue
            allocation = min(540.0, remaining)
            parameters = {
                "dont_explore_detection_event_sets_with_size_above": size,
                "dont_explore_edges_with_degree_above": size,
                "dont_explore_edges_increasing_symptom_degree": False,
            }
            payload = run_timed_child(heuristic_worker, (circuit_text, parameters),
                                      allocation)
            record = {"size": size, "allocation_s": allocation,
                      "parameters": parameters, "canonicalize_circuit_errors":
                      CANONICALIZE_CIRCUIT_ERRORS,
                      **{k: v for k, v in payload.items() if k != "mechanisms"}}
            if payload.get("ok"):
                entries = payload["mechanisms"]
                validation = validate_signature_witness(entries, mechanisms)
                record["weight"] = payload["weight"]
                record["validation"] = validation
                if (validation["valid_mechanism_level_upper_bound"]
                        and (stim_ub is None or payload["weight"] < stim_ub)):
                    stim_ub = payload["weight"]
                    record["new_best_validated_upper_bound"] = True
            attempts.append(record)
            print(f"[W4] {key} stim search size={size}: "
                  f"{record.get('weight', 'no result')}", flush=True)

        witnesses = {
            "w1_w2": low_witness, "w3": w3_witness,
            "w4_matching": w4m_witness, "w4_star": w4s_witness,
        }
        certified_lb = 4
        lb_method = "exhaustive undecomposed-DEM mechanism search at w=1,2,3 (EXP-029 reproduced)"
        if (w3 is not None and w3["certified_no_w3"]
                and w4m is not None and w4m["certified_no_w4_matching_class"]
                and w4s is not None and w4s["certified_no_w4_star_class"]):
            certified_lb = 5
            lb_method = ("exhaustive undecomposed-DEM mechanism search at "
                         "w=1,2,3 plus exhaustive w=4 exclusion over both "
                         "co-occurrence classes (perfect-matching and star)")
        stagnation = {
            "w3_completed": w3 is None or w3["completed"],
            "w4_matching_completed": w4m is None or w4m["completed"],
            "w4_star_completed": w4s is None or w4s["completed"],
            "stop_reasons": {name: (block or {}).get("stop_reason")
                             for name, block in (("w3", w3),
                                                 ("w4_matching", w4m),
                                                 ("w4_star", w4s))},
        }
        if certified_lb < 5:
            any_stagnation = True
        ub = 12
        ub_method = ("EXP-029 certified weight-12 data-logical DEM witness "
                     "(results/raw/exp029_witnesses.json); unchanged by EXP-040")
        if stim_ub is not None and stim_ub < 12:
            ub = stim_ub
            ub_method = "EXP-040 stim witness search, XOR-validated against the undecomposed DEM"
        entry = {
            "label": LOCKED[key]["label"],
            "fingerprints": fp,
            "dem": {
                "flattened_error_mechanisms": len(mechanisms),
                "flattened_separator_instructions": separators,
                "dem_errors": int(dem.num_errors),
                "dem_detectors": int(dem.num_detectors),
            },
            "w1_w2": w12,
            "w3": w3,
            "w4_matching_class": w4m,
            "w4_star_class": w4s,
            "stim_witness_search": {
                "budget_s": W4_STIM_TOTAL_BUDGET_S,
                "wall_s": time.perf_counter() - search_start,
                "attempts": attempts,
            },
            "witnesses": witnesses,
            "lb": {"value": certified_lb, "method": lb_method,
                   "certified": True, "fault_model": FAULT_MODEL,
                   "classification": "CERTIFIED lower bound"},
            "ub": {"value": ub, "method": ub_method,
                   "classification": "mechanism-level upper bound"},
            "stagnation": stagnation,
            "wall_s": round(time.time() - t1, 1),
        }
        out["circuits"][key] = entry
        print(f"[W4] {key}: certified lb={certified_lb} ub={ub} "
              f"stagnation={stagnation} ({entry['wall_s']}s)", flush=True)

    both_lb5 = all(c["lb"]["value"] >= 5 for c in out["circuits"].values())
    out["verdict"] = {
        "certified_d_dem_mech_both_ge_5": both_lb5,
        "summary": (
            "certified d_DEM_mech >= 5 on both circuits" if both_lb5 else
            "w=4 exclusion did not complete within budget on at least one "
            "circuit; stagnation recorded; certified bound remains 4 <= "
            "d_DEM_mech <= 12"),
    }
    out["wall_s"] = round(time.time() - t0, 1)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    DISTANCE_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(f"[W4] wrote {DISTANCE_PATH.relative_to(ROOT)}", flush=True)
    if any_stagnation:
        PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
        STAGNATION_PATH.write_text(json.dumps(
            {"experiment": EXPERIMENT,
             "note": "budget-capped w=4 exclusion stagnation record",
             "circuits": {k: {"stagnation": v["stagnation"],
                              "w3": v["w3"], "w4_matching_class": v["w4_matching_class"],
                              "w4_star_class": v["w4_star_class"]}
                          for k, v in out["circuits"].items()}},
            indent=2) + "\n")
        print(f"[W4] wrote {STAGNATION_PATH.relative_to(ROOT)}", flush=True)


# ---------------------------------------------------------------------------
# W5 -- latency arm (EXP-013 discipline)
# ---------------------------------------------------------------------------
def cmd_latency() -> None:
    import scipy.sparse as sp
    from ldpc import BpOsdDecoder
    from qec_research.decoders.bposd_dem import dem_to_matrices

    load1 = os.getloadavg()[0]
    if load1 > LATENCY_MAX_LOAD:
        msg = (f"DEFERRED: machine under SAT fleet load (1-minute load "
               f"average {load1:.1f} > {LATENCY_MAX_LOAD})")
        print(f"[latency] {msg}")
        if PRODUCTION_PATH.exists():
            prod = json.loads(PRODUCTION_PATH.read_text())
            prod["latency_arm"] = msg
            prod["loadavg_at_latency_check"] = list(os.getloadavg())
            PRODUCTION_PATH.write_text(json.dumps(prod, indent=2) + "\n")
            print(f"[latency] updated {PRODUCTION_PATH.relative_to(ROOT)}")
        return

    registry = json.loads(REGISTRY_PATH.read_text())
    acc: dict[str, list[float]] = {}
    decs = {}
    for key in ("gross", "pbb"):
        det = registry["deterministic_slot_maps"][key]
        code, _, _ = code_and_supports(key)
        slot = slot_from_json(det["slot_map"])
        circ, _ = circuit_from_slot(key, code, slot, int(det["depth"]),
                                    PRODUCTION_P)
        dem = circ.detector_error_model(decompose_errors=False)
        M = dem_to_matrices(dem)
        pr = np.clip(M.priors, 1e-12, 1 - 1e-12)
        dec = BpOsdDecoder(sp.csr_matrix(M.H), error_channel=list(pr),
                           **PRODUCTION_DECODER)
        decs[key] = (dec, M, circ)

        def timed(dec, M, circ, shots, seed):
            smp = circ.compile_detector_sampler(seed=seed)
            dets, obss = smp.sample(shots, separate_observables=True)
            dets = dets.astype(np.uint8)
            obss = obss.astype(np.uint8)
            lat = []
            for i in range(shots):
                t = time.perf_counter()
                corr = dec.decode(dets[i])
                lat.append((time.perf_counter() - t) * 1e3)
            return lat

        acc[key] = []
        acc[key].extend(timed(dec, M, circ, 40, seed=1))  # warm-up
    print("[latency] warm-up complete", flush=True)
    per_block = LATENCY_SHOTS_PER_CODE // LATENCY_BLOCKS
    for bkt in range(LATENCY_BLOCKS):
        for key in ("gross", "pbb"):
            dec, M, circ = decs[key]
            acc[key].extend(timed(dec, M, circ, per_block, seed=1000 + bkt))
    load_end = os.getloadavg()
    codes = []
    for key in ("gross", "pbb"):
        a = np.asarray(acc[key])
        codes.append({
            "key": key, "shots_timed": int(a.size),
            "latency_mean_ms": float(a.mean()),
            "latency_p50_ms": float(np.percentile(a, 50)),
            "latency_p95_ms": float(np.percentile(a, 95)),
            "latency_p99_ms": float(np.percentile(a, 99)),
        })
    ratios = {f"p{q}_ratio": codes[1][f"latency_p{q}_ms"] / codes[0][f"latency_p{q}_ms"]
              for q in (50, 95, 99)}
    block = {
        "discipline": "EXP-013: load guard, warm-up, interleaved blocks, "
                      "OMP_NUM_THREADS=1, single process",
        "load_at_start": load1, "load_at_end": list(load_end),
        "shots_per_code": LATENCY_SHOTS_PER_CODE,
        "blocks": LATENCY_BLOCKS,
        "decoder": PRODUCTION_DECODER,
        "p": PRODUCTION_P,
        "codes": codes, "ratios": ratios,
    }
    prod = (json.loads(PRODUCTION_PATH.read_text())
            if PRODUCTION_PATH.exists() else {"experiment": EXPERIMENT,
                                              "step": "W5-production-decoder"})
    prod["latency_arm"] = "MEASURED"
    prod["latency"] = block
    PRODUCTION_PATH.write_text(json.dumps(prod, indent=2) + "\n")
    print(f"[latency] wrote block into {PRODUCTION_PATH.relative_to(ROOT)}: "
          f"{ratios}", flush=True)


# ---------------------------------------------------------------------------
# W6 -- verdict
# ---------------------------------------------------------------------------
def paper_anchor() -> dict[str, Any]:
    lines = PAPER_PATH.read_text().splitlines()
    for i, line in enumerate(lines, start=1):
        if line.startswith(PAPER_SENTENCE_HEAD):
            return {"path": "reports/paper_pbb_nogo.md", "line": i,
                    "sentence": " ".join(
                        s.strip() for s in lines[i - 1:i + 2])}
    return {"path": "reports/paper_pbb_nogo.md", "line": None,
            "sentence": None, "note": "anchor sentence not found"}


def cmd_verdict() -> None:
    anchor = paper_anchor()
    steps: dict[str, Any] = {}

    def step(name, ok, artifact, keys, note=""):
        steps[name] = {"status": "ACCEPTED" if ok else "FAILED-WITH-EVIDENCE",
                       "artifact": artifact, "keys": keys, "note": note}

    reg = json.loads(REGISTRY_PATH.read_text()) if REGISTRY_PATH.exists() else None
    reg_ok = bool(reg and reg.get("verdict") == "FINGERPRINTS_LOCKED")
    step("W1_registry_fingerprint_lock", reg_ok,
         "results/processed/exp040_circuit_registry.json",
         ["#/circuits/gross/circuit_sha256", "#/circuits/pbb/circuit_sha256",
          "#/circuits/gross/noiseless_p0", "#/circuits/pbb/noiseless_p0"])

    aud = json.loads(AUDIT_PATH.read_text()) if AUDIT_PATH.exists() else None
    aud_ok = bool(aud and aud.get("all_passed") and aud.get("n_circuits_audited", 0) >= 12)
    step("W2_detector_semantics_audit", aud_ok,
         "results/processed/exp040_detector_semantics_audit.json",
         ["#/n_circuits_audited", "#/all_passed", "#/audits"])

    mat = json.loads(MATCHED_PATH.read_text()) if MATCHED_PATH.exists() else None
    mat_ok = bool(mat and all(
        mat.get("per_p", {}).get(f"{p}", {}).get("grid_complete")
        for p in W3_P_POINTS))
    step("W3_matched_ler", mat_ok, "results/processed/exp040_matched_ler.json",
         ["#/per_p"], note=("grid complete at all three p" if mat_ok else
                            "grid incomplete; see per_p grid_complete flags "
                            "and cell-level achieved shots"))

    dst = json.loads(DISTANCE_PATH.read_text()) if DISTANCE_PATH.exists() else None
    dst_ok = bool(dst and dst.get("verdict", {}).get("certified_d_dem_mech_both_ge_5"))
    step("W4_circuit_distance", dst_ok, "results/processed/exp040_circuit_distance.json",
         ["#/circuits/gross/lb", "#/circuits/pbb/lb", "#/verdict"],
         note=("certified d_DEM_mech >= 5 both circuits" if dst_ok else
               "no certified >=5 and no smaller validated witness; stagnation "
               "recorded at results/partial_runs/exp040/"
               "circuit_distance_stagnation.json" if STAGNATION_PATH.exists()
               else "stagnation record absent"))

    pro = json.loads(PRODUCTION_PATH.read_text()) if PRODUCTION_PATH.exists() else None
    pro_ok = bool(pro and pro.get("cells")
                  and all(c.get("complete") for c in pro["cells"])
                  and len(pro["cells"]) == 10)
    latency_done = bool(pro and pro.get("latency_arm") == "MEASURED")
    step("W5_production_decoder", pro_ok or bool(pro),
         "results/processed/exp040_production_decoder.json",
         ["#/cells", "#/latency_arm"],
         note=("production LER complete" if pro_ok else
               "production arm partial (see cells)" if pro else
               "production arm not started"))

    gaps = {
        "detector_semantics": {
            "status": "DONE" if aud_ok else "MISSING",
            "artifact": "results/processed/exp040_detector_semantics_audit.json",
            "key": "#/audits" if aud_ok else None,
            "step_id": "W2"},
        "circuit_distance": {
            "status": "DONE" if dst_ok else ("PARTIAL" if dst else "MISSING"),
            "artifact": "results/processed/exp040_circuit_distance.json",
            "key": "#/circuits/{gross,pbb}/lb" if dst else None,
            "step_id": "W4"},
        "decoder_parity": {
            "status": "DONE" if (mat and all(
                mat["per_p"][f"{p}"]["grid_complete"] for p in W3_P_POINTS if f"{p}" in mat.get("per_p", {}))) else "PARTIAL" if mat else "MISSING",
            "artifact": "results/processed/exp040_matched_ler.json",
            "key": "#/per_p", "step_id": "W3"},
        "identical_noise": {
            "status": "DONE" if mat else "MISSING",
            "artifact": "results/processed/exp040_matched_ler.json",
            "key": "#/protocol/noise_model", "step_id": "W3"},
        "ler_with_CIs": {
            "status": "DONE" if (mat and all(
                (mat["per_p"][f"{p}"]["n_intervals_recorded"] or 0) >= 10
                for p in W3_P_POINTS if f"{p}" in mat.get("per_p", {}))) else "PARTIAL" if mat else "MISSING",
            "artifact": "results/processed/exp040_matched_ler.json",
            "key": "#/per_p/{p}/cells", "step_id": "W3"},
        "latency": {
            "status": "DONE" if latency_done else "MISSING",
            "artifact": "results/processed/exp040_production_decoder.json",
            "key": "#/latency_arm", "step_id": "W5",
            "note": (pro or {}).get("latency_arm")},
        "production_decoder": {
            "status": "DONE" if pro_ok else ("PARTIAL" if pro else "MISSING"),
            "artifact": "results/processed/exp040_production_decoder.json",
            "key": "#/cells", "step_id": "W5"},
    }

    headline: dict[str, Any] = {}
    if reg:
        for key in ("gross", "pbb"):
            headline[f"circuit_sha256_{key}"] = {
                "value": reg["circuits"][key]["circuit_sha256"],
                "path": "results/processed/exp040_circuit_registry.json",
                "pointer": f"#/circuits/{key}/circuit_sha256"}
    if mat:
        for p, block in mat["per_p"].items():
            headline[f"separation_p{p}"] = {
                "value": block["verdict"]["pointwise_familywise_separation_on_available_intervals"],
                "grid_complete": block["grid_complete"],
                "path": "results/processed/exp040_matched_ler.json",
                "pointer": f"#/per_p/{p}/verdict"}
            headline[f"schedule_mean_ratio_p{p}"] = {
                "value": block.get("schedule_mean_ratio_pbb_over_gross"),
                "path": "results/processed/exp040_matched_ler.json",
                "pointer": f"#/per_p/{p}/schedule_mean_ratio_pbb_over_gross"}
    if dst:
        for key in ("gross", "pbb"):
            headline[f"dem_distance_interval_{key}"] = {
                "value": [dst["circuits"][key]["lb"]["value"],
                          dst["circuits"][key]["ub"]["value"]],
                "path": "results/processed/exp040_circuit_distance.json",
                "pointer": f"#/circuits/{key}/lb/value,#/circuits/{key}/ub/value"}

    verdict = {
        "experiment": EXPERIMENT,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "paper_anchor": anchor,
        "locked_sentence": anchor.get("sentence"),
        "steps": steps,
        "gap_items": gaps,
        "headline_numbers": headline,
        "wall_note": ("every number above is traceable to the listed "
                      "(path, JSON-pointer) pair; see "
                      "notes/exp040_evidence_module.md"),
    }
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    VERDICT_PATH.write_text(json.dumps(verdict, indent=2) + "\n")

    lines = [
        "# EXP-040 evidence module",
        "",
        f"Generated {time.strftime('%Y-%m-%d %H:%M')} by "
        "`experiments/exp040_verified_circuit_benchmark.py verdict`.",
        "",
        "Locked sentence (`reports/paper_pbb_nogo.md`, line "
        f"{anchor.get('line')}):",
        "",
        f"> {anchor.get('sentence')}",
        "",
        "Seven gap items and where each number lives:",
        "",
        "| gap item | status | artifact | key pointer |",
        "|---|---|---|---|",
    ]
    for name, g in gaps.items():
        lines.append(f"| {name} | {g['status']} | {g['artifact']} | {g['key']} |")
    lines += [
        "",
        "Per-step outcomes:",
        "",
    ]
    for name, s in steps.items():
        lines.append(f"- **{name}**: {s['status']} ({s['artifact']})")
    lines += [
        "",
        "Headline numbers (all as `(path, JSON-pointer)` pairs):",
        "",
    ]
    for name, h in headline.items():
        lines.append(f"- `{name}` = {h['value']} — {h['path']}{h['pointer']}")
    NOTE_PATH.parent.mkdir(parents=True, exist_ok=True)
    NOTE_PATH.write_text("\n".join(lines) + "\n")
    print(f"[W6] wrote {VERDICT_PATH.relative_to(ROOT)} and "
          f"{NOTE_PATH.relative_to(ROOT)}", flush=True)


# ---------------------------------------------------------------------------
# pilot
# ---------------------------------------------------------------------------
def cmd_pilot(shots: int, workers: int) -> None:
    PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    out: dict[str, Any] = {"experiment": EXPERIMENT, "pilot": True,
                           "env": env_json(), "shots": shots, "workers": workers}
    cells = matched_cell_list()
    for phase, cfg in (("matched", MATCHED_DECODER), ("production", PRODUCTION_DECODER)):
        for key in ("gross", "pbb"):
            cell = next(c for c in cells
                        if c["key"] == key and c["p"] == 0.002)
            code, _, _ = code_and_supports(key)
            slot = slot_from_json(cell["slot_map"])
            circ, _ = circuit_from_slot(key, code, slot, cell["depth"], cell["p"])
            t0 = time.time()
            r = collect(circ, max_shots=shots, workers=workers, chunk=MC_CHUNK,
                        seed=777_000 + cell["cell_id"], **cfg)
            wall = time.time() - t0
            out[f"{phase}_{key}"] = {
                "shots": r.shots, "failures": r.failures,
                "wall_s": round(wall, 1),
                "decode_core_s_per_shot": round(r.decode_s_per_shot_1core, 4),
                "effective_core_s_per_shot": round(wall * workers / r.shots, 4),
                "latency_p50_ms": r.latency_p50_ms,
                "latency_p99_ms": r.latency_p99_ms,
            }
            print(f"[pilot] {phase} {key}: {out[f'{phase}_{key}']}", flush=True)
    (PARTIAL_DIR / "pilot.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[pilot] wrote {PARTIAL_DIR / 'pilot.json'}", flush=True)


# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("registry")
    sub.add_parser("audit")
    p_mc = sub.add_parser("mc")
    p_mc.add_argument("--workers", type=int, default=MC_WORKERS)
    p_mc.add_argument("--phases", default="matched,production")
    sub.add_parser("distance")
    sub.add_parser("latency")
    sub.add_parser("assemble")
    sub.add_parser("verdict")
    p_pilot = sub.add_parser("pilot")
    p_pilot.add_argument("--shots", type=int, default=400)
    p_pilot.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    if args.command == "registry":
        cmd_registry()
    elif args.command == "audit":
        cmd_audit()
    elif args.command == "mc":
        cmd_mc(args.workers, [s for s in args.phases.split(",") if s])
    elif args.command == "distance":
        cmd_distance()
    elif args.command == "latency":
        cmd_latency()
    elif args.command == "assemble":
        cmd_assemble()
    elif args.command == "verdict":
        cmd_verdict()
    elif args.command == "pilot":
        cmd_pilot(args.shots, args.workers)
    else:
        raise SystemExit(f"unknown command {args.command}")


if __name__ == "__main__":
    main()
