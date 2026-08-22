"""EXP-024: two-ancilla cat extraction for PBB ``12_6_0193``.

The verdict depth is the number of *all* two-qubit-gate layers in a round.
Cat preparation and merge CNOTs therefore count, even when they overlap data
couplings belonging to other checks.  A separate data-only scheduling problem
is solved to expose (but not confuse with the verdict) the apparent 8 -> 7
improvement obtained by omitting that plumbing.

Cat policy
----------
The experiment uses an unverified two-qubit cat.  Every weight-8 mixed check is
split 4+4 across the Bell-pair controls; there is no third verification/flag
ancilla.  A fault on either cat ancilla can consequently propagate to as many
as four data qubits.  The undecomposed Stim detector error model retains these
correlated hyperedges, but this experiment is not a fault-tolerance proof.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any, NoReturn

# Pre-registered protocol constants.  The production artifact is accepted only
# when these exact values are used.
PHYSICAL_ERROR_RATE = 0.002
ROUNDS = 12
SHOTS = 20_000
BENCHMARK_SEEDS = (20_260_812, 20_260_813, 20_260_814)
SHOTS_BY_SEED = (6_667, 6_667, 6_666)
NOISELESS_SHOTS = 4_096
NOISELESS_SEED = 20_260_812
MAX_MONTE_CARLO_WORKERS = 6
MAX_SOLVER_WORKERS = 8
INTEGRATED_DEPTH7_TIME_LIMIT_S = 60.0
INTEGRATED_DEPTH8_BUDGET_S = 60.0
CONSTRUCTIVE_DEPTH = 9
DATA_ONLY_DEPTH = 7
MAX_ITER = 30
BP_METHOD = "ms"
MS_SCALING_FACTOR = 0.625
OSD_METHOD = "osd0"
OSD_ORDER = 0
CAT_POLICY = "unverified_two_qubit_cat_no_flag"
TARGET_CODE_ID = "12_6_0193"
BASELINE_SCHEDULE_POLICY = "first persisted uniformly sampled EXP-016 schedule"

# One C++ decoder is run by each process.  Do not permit hidden BLAS fan-out.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np  # noqa: E402
import stim  # noqa: E402
from ortools.sat.python import cp_model  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits,
    pbb_supports_and_orbits,
    pure_z_logical_basis,
)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec,
    PauliSupport,
    build_memory_circuit,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    anticommuting_overlaps,
    depth_lower_bound,
    parity_defects,
    slots_to_layers,
    verify_schedule,
)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.decoders.bposd_dem import clopper_pearson  # noqa: E402
from qec_research.decoders.parallel_harness import collect  # noqa: E402
from qec_research.artifacts import canonical_route  # noqa: E402

CATALOGUE = (
    ROOT
    / "third_party"
    / "qcode-discovery"
    / "results"
    / "campaign7_publication_merged.jsonl"
)
EXP016_RESULT = ROOT / "results" / "raw" / "exp016_schedule_controlled.json"
OUTPUT_NAME = "exp024_cat_extraction.json"
REPORT_NAME = "exp024_cat_extraction.md"
OUTPUT = ROOT / "results" / "processed" / OUTPUT_NAME
REPORT = ROOT / "notes" / "agent_reports" / REPORT_NAME
STRUCTURAL_OUTPUT = (
    ROOT / "results" / "processed" / "exp024_cat_extraction_structural.json"
)
STRUCTURAL_REPORT = (
    ROOT / "notes" / "agent_reports" / "exp024_cat_extraction_structural.md"
)
STRUCTURAL_QUARANTINE = (
    ROOT / "results" / "quarantine" / "exp024_cat_extraction_structural-failing.json"
)
FEASIBILITY_AMENDMENT = (
    ROOT / "results" / "partial_runs" / "exp024_benchmark_feasibility_amendment.json"
)
V2_GATE = ROOT / "results" / "partial_runs" / "exp024_benchmark_v2_harness_gate.json"
V2_CHECKPOINT_DIR = ROOT / "results" / "partial_runs" / "exp024_benchmark_v2"
EXP026_HARNESS = ROOT / "experiments" / "exp026_decoder_codesign.py"
EXP026_HARNESS_TEST = ROOT / "tests" / "test_exp026_decoder_codesign.py"
EXP026_HARNESS_SHA256 = "257e0b1827b4e18b5d70fabd67751922eadc482cfa5fb36e70a0961d38a8d9f2"
EXP026_TEST_SHA256 = "04fac68276d3d6ffcd950bda03a6a2ac8000be51d7f95fe3bdc6c6d659f36efe"
EXP016_SHA256 = "49f70b8c1d8bbbe6e8a85296d74e93f156760013e40ed0c60207a7d3f207cdb3"

_PAULI_GATE = {"X": "CX", "Y": "CY", "Z": "CZ"}


def _versions() -> dict[str, str]:
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ("numpy", "stim", "scipy", "ldpc", "ortools"):
        try:
            out[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            out[name] = "not-installed"
    return out

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def _atomic_write_json(path: Path, value: Any) -> None:
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True))

def benchmark_v2_protocol() -> dict[str, Any]:
    amendment = json.loads(FEASIBILITY_AMENDMENT.read_text())
    return amendment["v2_frozen_protocol"]


def benchmark_v2_protocol_sha256() -> str:
    payload = json.dumps(
        benchmark_v2_protocol(), sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def require_benchmark_v2_execution_gate() -> dict[str, Any]:
    """Refuse production until the shared supervised harness is attested."""
    if not V2_GATE.exists():
        raise RuntimeError(
            "EXP-024 benchmark v2 is frozen but blocked: missing supervised "
            f"harness gate {V2_GATE}"
        )
    gate = json.loads(V2_GATE.read_text())
    expected = {
        "protocol_sha256": benchmark_v2_protocol_sha256(),
        "exp026_harness_sha256": EXP026_HARNESS_SHA256,
        "exp026_test_sha256": EXP026_TEST_SHA256,
    }
    if any(gate.get(key) != value for key, value in expected.items()):
        raise RuntimeError("EXP-024 benchmark v2 gate binding mismatch")
    required = (
        "timeout_restart_verified",
        "nonmultiple_resume_verified",
        "binding_drift_rejection_verified",
        "full_coverage_route_verified",
        "atomic_arm_checkpoint_verified",
        "development_isolation_verified",
        "production_execution_authorized",
    )
    missing = [key for key in required if gate.get(key) is not True]
    if missing:
        raise RuntimeError(f"EXP-024 benchmark v2 gate is incomplete: {missing}")
    if _sha256(EXP026_HARNESS) != EXP026_HARNESS_SHA256:
        raise RuntimeError("shared EXP-026 harness changed after verification")
    if _sha256(EXP026_HARNESS_TEST) != EXP026_TEST_SHA256:
        raise RuntimeError("shared EXP-026 harness tests changed after verification")
    return gate


def run_benchmark_v2_production() -> NoReturn:
    """Fail closed until the frozen supervised v2 driver is implemented."""
    require_benchmark_v2_execution_gate()
    raise RuntimeError(
        "EXP-024 benchmark v2 production is blocked: the supervised, "
        "deadline-bounded, resumable driver is not implemented; the legacy "
        "_benchmark/collect path is retired and cannot write canonical output"
    )


def benchmark_v2_binding(
    phase: str,
    arm: str,
    circuit_sha256: str,
    seed: int,
    sample_sha256: str,
    expected_shots: int,
) -> dict[str, Any]:
    if phase not in ("development", "heldout"):
        raise ValueError(phase)
    if arm not in ("cat_pbb", "one_ancilla_pbb", "gross"):
        raise ValueError(arm)
    return {
        "protocol_sha256": benchmark_v2_protocol_sha256(),
        "phase": phase,
        "arm": arm,
        "circuit_sha256": circuit_sha256,
        "seed": int(seed),
        "sample_sha256": sample_sha256,
        "expected_shots": int(expected_shots),
    }


def _validate_v2_records(
    binding: dict[str, Any], records: list[dict[str, Any]]
) -> None:
    shot_ids = [int(record["shot_id"]) for record in records]
    if shot_ids != list(range(len(records))):
        raise RuntimeError("v2 checkpoint is not a unique contiguous shot prefix")
    if len(records) > int(binding["expected_shots"]):
        raise RuntimeError("v2 checkpoint exceeds frozen expected shots")


def save_benchmark_v2_arm_checkpoint(
    binding: dict[str, Any], records: list[dict[str, Any]]
) -> Path:
    """Atomically persist one bound arm prefix; safe to resume after interruption."""
    _validate_v2_records(binding, records)
    path = (
        V2_CHECKPOINT_DIR
        / f"{binding['phase']}-{binding['arm']}.json"
    )
    _atomic_write_json(
        path,
        {
            "status": (
                "complete"
                if len(records) == int(binding["expected_shots"])
                else "partial"
            ),
            "artifact_route": canonical_route(
                clean=True,
                full_coverage=False,
            ),
            "binding": binding,
            "records": records,
        },
    )
    return path


def load_benchmark_v2_arm_checkpoint(
    binding: dict[str, Any],
) -> list[dict[str, Any]]:
    path = V2_CHECKPOINT_DIR / f"{binding['phase']}-{binding['arm']}.json"
    if not path.exists():
        return []
    checkpoint = json.loads(path.read_text())
    if checkpoint.get("binding") != binding:
        raise RuntimeError("v2 checkpoint protocol/circuit/sample binding mismatch")
    records = checkpoint.get("records", [])
    _validate_v2_records(binding, records)
    return records


def _load_target() -> tuple[dict[str, Any], PBBSpec]:
    row = None
    with CATALOGUE.open() as f:
        for line in f:
            candidate = json.loads(line)
            if candidate.get("code_id") == TARGET_CODE_ID:
                row = candidate
                break
    if row is None:
        raise RuntimeError(f"catalogue member {TARGET_CODE_ID!r} not found")
    if (row["n"], row["k"], row["d"]) != (144, 12, 12):
        raise RuntimeError("target catalogue parameters changed")
    spec = PBBSpec(
        row["ell"],
        row["m"],
        row["A_terms"],
        row["B_terms"],
        row["C_terms"],
        row["D_terms"],
    )
    return row, spec


def _exp016_slot(label_prefix: str) -> tuple[dict[tuple[int, int], int], dict[str, Any]]:
    prior = json.loads(EXP016_RESULT.read_text())
    entry = next(x for x in prior["codes"] if x["label"].startswith(label_prefix))
    run = entry["runs"][0]
    slot = {(int(c), int(q)): int(t) for c, q, t in run["slot_map"]}
    source = {
        "experiment": "EXP-016",
        "selection_policy": BASELINE_SCHEDULE_POLICY,
        "schedule_class": "translation-invariant (orbit)",
        "schedule_index": int(run["schedule_index"]),
        "schedule_hash": run["schedule_hash"],
        "depth": int(entry["depth"]),
        "enumeration_complete_within_class": bool(entry["enumeration_complete"]),
        "total_valid_schedules_within_class": int(entry["total_valid_schedules"]),
    }
    return slot, source


def _add_semantic_parity_constraints(
    model: cp_model.CpModel,
    supports: list[PauliSupport],
    data_time: dict[tuple[int, int], cp_model.IntVar],
) -> None:
    """Impose the controlled-Pauli even-crossing measurement criterion."""
    for (a, b), overlap in anticommuting_overlaps(supports).items():
        before = []
        for q in overlap:
            bit = model.new_bool_var(f"before_{a}_{b}_{q}")
            model.add(data_time[(a, q)] < data_time[(b, q)]).only_enforce_if(bit)
            model.add(data_time[(a, q)] > data_time[(b, q)]).only_enforce_if(
                bit.negated()
            )
            before.append(bit)
        half_count = model.new_int_var(0, len(before) // 2, f"parity_{a}_{b}")
        model.add(sum(before) == 2 * half_count)


def _solve_unrestricted_integrated(
    supports: list[PauliSupport],
    n: int,
    depth: int,
    time_limit_s: float,
) -> dict[str, Any]:
    """Search all 4+4 splits and all schedules at one integrated depth.

    The only symmetry break fixes the lower-index support edge to physical cat
    half 0.  Exchanging the two Bell controls is an exact symmetry, so this does
    not remove a physical extraction circuit.
    """
    started = time.time()
    model = cp_model.CpModel()
    edges = [(s.index, q) for s in supports for q in s.paulis]
    data_time = {
        edge: model.new_int_var(0, depth - 1, f"data_{edge[0]}_{edge[1]}")
        for edge in edges
    }
    perturbed = [s for s in supports if s.weight == 8]
    prep_time = {
        s.index: model.new_int_var(0, depth - 1, f"prep_{s.index}")
        for s in perturbed
    }
    merge_time = {
        s.index: model.new_int_var(0, depth - 1, f"merge_{s.index}")
        for s in perturbed
    }
    half_var: dict[tuple[int, int], cp_model.IntVar] = {}

    for s in supports:
        qs = sorted(s.paulis)
        if s.weight != 8:
            model.add_all_different([data_time[(s.index, q)] for q in qs])
            continue
        bits = []
        for q in qs:
            bit = model.new_bool_var(f"half_{s.index}_{q}")
            half_var[(s.index, q)] = bit
            bits.append(bit)
            model.add(prep_time[s.index] < data_time[(s.index, q)])
            model.add(data_time[(s.index, q)] < merge_time[s.index])
        model.add(sum(bits) == 4)
        model.add(bits[0] == 0)
        for q1, q2 in itertools.combinations(qs, 2):
            same_half = model.new_bool_var(f"same_half_{s.index}_{q1}_{q2}")
            model.add(
                half_var[(s.index, q1)] == half_var[(s.index, q2)]
            ).only_enforce_if(same_half)
            model.add(
                half_var[(s.index, q1)] != half_var[(s.index, q2)]
            ).only_enforce_if(same_half.negated())
            model.add(
                data_time[(s.index, q1)] != data_time[(s.index, q2)]
            ).only_enforce_if(same_half)

    for q in range(n):
        touching = [data_time[(s.index, q)] for s in supports if q in s.paulis]
        model.add_all_different(touching)
    _add_semantic_parity_constraints(model, supports, data_time)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = MAX_SOLVER_WORKERS
    solver.parameters.random_seed = BENCHMARK_SEEDS[0]
    status = solver.solve(model)
    status_name = solver.status_name(status)
    result: dict[str, Any] = {
        "depth": depth,
        "status": status_name,
        "scope": "unrestricted over per-check 4+4 partitions and non-translation-invariant slots",
        "time_limit_s": time_limit_s,
        "solver_workers": MAX_SOLVER_WORKERS,
        "wall_s": time.time() - started,
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result["schedule"] = {
            "depth": depth,
            "data_slot": {edge: int(solver.value(var)) for edge, var in data_time.items()},
            "half": {
                edge: int(solver.value(var)) for edge, var in half_var.items()
            },
            "prep_slot": {
                c: int(solver.value(var)) for c, var in prep_time.items()
            },
            "merge_slot": {
                c: int(solver.value(var)) for c, var in merge_time.items()
            },
        }
    return result


def _solve_invariant(
    supports: list[PauliSupport],
    n: int,
    orbits: dict[tuple[int, int], int],
    depth: int,
    include_plumbing: bool,
) -> dict[str, Any]:
    """Solve the deterministic translation-invariant 0..3 / 4..7 split."""
    started = time.time()
    model = cp_model.CpModel()
    orbit_ids = sorted(set(orbits.values()))
    orbit_time = {
        o: model.new_int_var(0, depth - 1, f"orbit_{o}") for o in orbit_ids
    }
    data_time = {edge: orbit_time[orbit] for edge, orbit in orbits.items()}
    mixed_orbits = sorted(
        {
            orbits[(s.index, q)]
            for s in supports
            if s.weight == 8
            for q in s.paulis
        }
    )
    if len(mixed_orbits) != 8:
        raise RuntimeError(f"expected 8 perturbed-check orbits, got {mixed_orbits}")
    primary_orbits = set(mixed_orbits[:4])
    half = {
        (s.index, q): int(orbits[(s.index, q)] not in primary_orbits)
        for s in supports
        if s.weight == 8
        for q in s.paulis
    }

    prep = merge = None
    if include_plumbing:
        prep = model.new_int_var(0, depth - 1, "cat_prep")
        merge = model.new_int_var(0, depth - 1, "cat_merge")

    for s in supports:
        if s.weight == 8:
            for h in (0, 1):
                values = [
                    data_time[(s.index, q)]
                    for q in s.paulis
                    if half[(s.index, q)] == h
                ]
                if include_plumbing:
                    values = [prep, merge, *values]
                model.add_all_different(values)
            if include_plumbing:
                for q in s.paulis:
                    model.add(prep < data_time[(s.index, q)])
                    model.add(data_time[(s.index, q)] < merge)
        else:
            model.add_all_different(
                [data_time[(s.index, q)] for q in s.paulis]
            )
    for q in range(n):
        model.add_all_different(
            [data_time[(s.index, q)] for s in supports if q in s.paulis]
        )
    _add_semantic_parity_constraints(model, supports, data_time)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 120.0
    solver.parameters.num_workers = 1
    solver.parameters.random_seed = BENCHMARK_SEEDS[0]
    status = solver.solve(model)
    status_name = solver.status_name(status)
    result: dict[str, Any] = {
        "depth": depth,
        "status": status_name,
        "scope": "translation-invariant fixed orbit split",
        "primary_orbits": sorted(primary_orbits),
        "secondary_orbits": sorted(set(mixed_orbits) - primary_orbits),
        "solver_workers": 1,
        "wall_s": time.time() - started,
    }
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return result
    slot = {edge: int(solver.value(var)) for edge, var in data_time.items()}
    schedule = {
        "depth": depth,
        "data_slot": slot,
        "half": half,
        "prep_slot": {},
        "merge_slot": {},
    }
    if include_plumbing:
        p = int(solver.value(prep))
        m = int(solver.value(merge))
        schedule["prep_slot"] = {
            s.index: p for s in supports if s.weight == 8
        }
        schedule["merge_slot"] = {
            s.index: m for s in supports if s.weight == 8
        }
        result["shared_prep_slot"] = p
        result["shared_merge_slot"] = m
    result["orbit_slot"] = {
        str(o): int(solver.value(var)) for o, var in orbit_time.items()
    }
    result["schedule"] = schedule
    return result


def _cat_ancillas(
    n: int, supports: list[PauliSupport]
) -> tuple[dict[int, int], dict[int, int], int]:
    primary = {s.index: n + s.index for s in supports}
    perturbed = [s.index for s in supports if s.weight == 8]
    secondary = {c: n + len(supports) + i for i, c in enumerate(perturbed)}
    return primary, secondary, n + len(supports) + len(perturbed)


def _verify_cat_schedule(
    supports: list[PauliSupport], n: int, schedule: dict[str, Any]
) -> dict[str, Any]:
    data_slot = schedule["data_slot"]
    half = schedule["half"]
    prep_slot = schedule["prep_slot"]
    merge_slot = schedule["merge_slot"]
    depth = schedule["depth"]
    primary, secondary, total_qubits = _cat_ancillas(n, supports)
    expected_edges = {(s.index, q) for s in supports for q in s.paulis}
    half_counts: dict[int, list[int]] = {}
    layer_ops = [dict(prep=0, data=0, merge=0) for _ in range(depth)]
    touched: list[set[int]] = [set() for _ in range(depth)]
    collisions: list[tuple[int, int, str]] = []

    def touch(t: int, qubits: list[int], kind: str) -> None:
        if not 0 <= t < depth:
            collisions.append((t, -1, f"out_of_range_{kind}"))
            return
        for q in qubits:
            if q in touched[t]:
                collisions.append((t, q, kind))
            touched[t].add(q)
        layer_ops[t][kind] += 1

    for s in supports:
        if s.weight == 8:
            counts = [0, 0]
            touch(prep_slot[s.index], [primary[s.index], secondary[s.index]], "prep")
            touch(merge_slot[s.index], [primary[s.index], secondary[s.index]], "merge")
            for q in s.paulis:
                h = int(half[(s.index, q)])
                if h not in (0, 1):
                    collisions.append((-1, q, "invalid_half"))
                    continue
                counts[h] += 1
                anc = primary[s.index] if h == 0 else secondary[s.index]
                touch(data_slot[(s.index, q)], [anc, q], "data")
                if not (
                    prep_slot[s.index]
                    < data_slot[(s.index, q)]
                    < merge_slot[s.index]
                ):
                    collisions.append((data_slot[(s.index, q)], anc, "precedence"))
            half_counts[s.index] = counts
        else:
            for q in s.paulis:
                touch(data_slot[(s.index, q)], [primary[s.index], q], "data")

    defects = parity_defects(supports, data_slot)
    valid = (
        set(data_slot) == expected_edges
        and not collisions
        and not defects
        and all(counts == [4, 4] for counts in half_counts.values())
    )
    return {
        "valid": valid,
        "covers_all_edges": set(data_slot) == expected_edges,
        "physical_qubit_conflict_free": not collisions,
        "collision_examples": collisions[:10],
        "all_weight8_splits_are_4_plus_4": all(
            counts == [4, 4] for counts in half_counts.values()
        ),
        "num_weight8_checks": len(half_counts),
        "semantic_parity_defects": len(defects),
        "semantic_parity_defect_examples": defects[:10],
        "layer_operation_counts": layer_ops,
        "num_qubits_total": total_qubits,
    }


def _schedule_hash(schedule: dict[str, Any]) -> str:
    rows = []
    for (c, q), t in sorted(schedule["data_slot"].items()):
        rows.append(f"D:{c}:{q}:{schedule['half'].get((c, q), 0)}:{t}")
    for c, t in sorted(schedule["prep_slot"].items()):
        rows.append(f"P:{c}:{t}")
    for c, t in sorted(schedule["merge_slot"].items()):
        rows.append(f"M:{c}:{t}")
    return hashlib.sha256("|".join(rows).encode()).hexdigest()[:16]


def _serialise_schedule(schedule: dict[str, Any]) -> dict[str, Any]:
    return {
        "depth": int(schedule["depth"]),
        "schedule_hash": _schedule_hash(schedule),
        "data_slot": [
            [int(c), int(q), int(schedule["half"].get((c, q), 0)), int(t)]
            for (c, q), t in sorted(schedule["data_slot"].items())
        ],
        "prep_slot": [
            [int(c), int(t)] for c, t in sorted(schedule["prep_slot"].items())
        ],
        "merge_slot": [
            [int(c), int(t)] for c, t in sorted(schedule["merge_slot"].items())
        ],
    }


def build_cat_memory_circuit(
    H: np.ndarray,
    observables: np.ndarray,
    rounds: int,
    p: float,
    schedule: dict[str, Any],
) -> tuple[stim.Circuit, dict[str, Any]]:
    """Build the unverified 2-cat memory circuit without changing src/."""
    H = np.asarray(H, dtype=np.uint8) & 1
    observables = np.asarray(observables, dtype=np.uint8) & 1
    n = H.shape[1] // 2
    from qec_research.circuits.mixed_stabilizer import generator_supports

    supports = generator_supports(H)
    verification = _verify_cat_schedule(supports, n, schedule)
    if not verification["valid"]:
        raise RuntimeError(f"invalid cat schedule: {verification}")
    primary, secondary, total_qubits = _cat_ancillas(n, supports)
    all_qubits = list(range(total_qubits))
    perturbed = [s.index for s in supports if s.weight == 8]
    unperturbed = [s.index for s in supports if s.weight != 8]
    cat_ancillas = [q for c in perturbed for q in (primary[c], secondary[c])]
    unperturbed_ancillas = [primary[c] for c in unperturbed]
    depth = int(schedule["depth"])

    c = stim.Circuit()
    c.append("R", list(range(n)))
    if p:
        c.append("X_ERROR", list(range(n)), p)

    measurement_count = 0
    measurement_index: list[dict[int, int]] = []
    det_ready = [s.is_pure_z for s in supports]

    for rd in range(rounds):
        # One single-qubit preparation layer.  R + H has the same preparation
        # flip channel as RX: X after R becomes Z after H.
        c.append("R", cat_ancillas)
        if p:
            c.append("X_ERROR", cat_ancillas, p)
        c.append("H", [primary[ci] for ci in perturbed])
        c.append("RX", unperturbed_ancillas)
        if p:
            c.append("Z_ERROR", unperturbed_ancillas, p)
        c.append("TICK")

        for layer in range(depth):
            by_gate: dict[str, list[int]] = {}
            noisy_pairs: list[int] = []
            busy: set[int] = set()

            def add_two_qubit(gate: str, a: int, b: int) -> None:
                if a in busy or b in busy:
                    raise RuntimeError(f"layer {layer} reuses physical qubit")
                busy.add(a)
                busy.add(b)
                by_gate.setdefault(gate, []).extend([a, b])
                noisy_pairs.extend([a, b])

            for ci, t in schedule["prep_slot"].items():
                if t == layer:
                    add_two_qubit("CX", primary[ci], secondary[ci])
            for s in supports:
                for q, pauli in s.paulis.items():
                    if schedule["data_slot"][(s.index, q)] != layer:
                        continue
                    h = schedule["half"].get((s.index, q), 0)
                    ancilla = primary[s.index] if h == 0 else secondary[s.index]
                    add_two_qubit(_PAULI_GATE[pauli], ancilla, q)
            for ci, t in schedule["merge_slot"].items():
                if t == layer:
                    add_two_qubit("CX", primary[ci], secondary[ci])

            for gate in sorted(by_gate):
                c.append(gate, by_gate[gate])
            if p:
                if noisy_pairs:
                    c.append("DEPOLARIZE2", noisy_pairs, p)
                idle = [q for q in all_qubits if q not in busy]
                if idle:
                    c.append("DEPOLARIZE1", idle, p)
            c.append("TICK")

        if p:
            c.append("Z_ERROR", [primary[i] for i in range(len(supports))], p)
        c.append("MX", [primary[i] for i in range(len(supports))])
        this_round = {}
        for i in range(len(supports)):
            this_round[i] = measurement_count
            measurement_count += 1
        measurement_index.append(this_round)

        for i in range(len(supports)):
            if rd == 0:
                if det_ready[i]:
                    off = measurement_count - measurement_index[0][i]
                    c.append("DETECTOR", [stim.target_rec(-off)], [i, rd])
            else:
                current = measurement_count - measurement_index[rd][i]
                previous = measurement_count - measurement_index[rd - 1][i]
                c.append(
                    "DETECTOR",
                    [stim.target_rec(-current), stim.target_rec(-previous)],
                    [i, rd],
                )
        c.append("TICK")

    if p:
        c.append("X_ERROR", list(range(n)), p)
    c.append("M", list(range(n)))
    data_measurement_base = measurement_count
    measurement_count += n

    for i, s in enumerate(supports):
        if not det_ready[i]:
            continue
        targets = [
            stim.target_rec(-(measurement_count - (data_measurement_base + q)))
            for q in s.paulis
        ]
        targets.append(
            stim.target_rec(
                -(measurement_count - measurement_index[rounds - 1][i])
            )
        )
        c.append("DETECTOR", targets, [i, rounds])

    num_observables = 0
    for obs in observables:
        if obs[:n].any():
            continue
        targets = [
            stim.target_rec(-(measurement_count - (data_measurement_base + int(q))))
            for q in np.flatnonzero(obs[n:])
        ]
        c.append("OBSERVABLE_INCLUDE", targets, num_observables)
        num_observables += 1

    data_gates = sum(s.weight for s in supports)
    cat_checks = len(perturbed)
    meta = {
        "n_data": n,
        "n_ancilla": len(supports) + cat_checks,
        "n_qubits_total": total_qubits,
        "num_checks": len(supports),
        "num_weight8_cat_checks": cat_checks,
        "num_unperturbed_one_ancilla_checks": len(unperturbed),
        "data_two_qubit_gates_per_round": data_gates,
        "cat_prep_cnot_per_round": cat_checks,
        "cat_merge_cnot_per_round": cat_checks,
        "total_two_qubit_gates_per_round": data_gates + 2 * cat_checks,
        "total_two_qubit_depth_per_round": depth,
        "single_qubit_prep_layers_per_round": 1,
        "measurement_layers_per_round": 1,
        "full_operation_layer_depth_per_round": depth + 2,
        "rounds": rounds,
        "p": p,
        "basis": "Z",
        "num_observables": num_observables,
        "num_detectors": c.num_detectors,
        "schedule_verification": verification,
    }
    return c, meta


def _observable_definitions(circuit: stim.Circuit) -> list[str]:
    return [
        str(instruction)
        for instruction in circuit
        if instruction.name == "OBSERVABLE_INCLUDE"
    ]


def _noiseless_evidence(
    circuit: stim.Circuit, seed: int
) -> dict[str, int]:
    detector, observable = circuit.compile_detector_sampler(seed=seed).sample(
        NOISELESS_SHOTS, separate_observables=True
    )
    return {
        "shots": NOISELESS_SHOTS,
        "seed": seed,
        "detector_firings": int(detector.sum()),
        "observable_flips": int(observable.sum()),
        "num_detectors": int(detector.shape[1]),
        "num_observables": int(observable.shape[1]),
    }


def _benchmark(
    label: str, circuit: stim.Circuit, workers: int
) -> dict[str, Any]:
    started = time.time()
    subruns = []
    total_shots = 0
    total_failures = 0
    for seed, shots in zip(BENCHMARK_SEEDS, SHOTS_BY_SEED):
        run = collect(
            circuit,
            max_shots=shots,
            workers=workers,
            chunk=200,
            seed=seed,
            max_iter=MAX_ITER,
            osd_order=OSD_ORDER,
            osd_method=OSD_METHOD,
            bp_method=BP_METHOD,
            ms_scaling_factor=MS_SCALING_FACTOR,
        )
        record = run.to_dict()
        record["seed"] = seed
        if run.failures == 0:
            record["ler_per_shot"] = None
            record["zero_failure_statement"] = (
                f"0/{run.shots}; 95% upper bound {run.ci95_hi:.8g}, not zero"
            )
        subruns.append(record)
        total_shots += run.shots
        total_failures += run.failures
    lo, hi = clopper_pearson(total_failures, total_shots)
    rate = total_failures / total_shots if total_failures else None
    return {
        "label": label,
        "shots": total_shots,
        "failures": total_failures,
        "ler_per_shot": rate,
        "ci95_clopper_pearson": [lo, hi],
        "subruns": subruns,
        "wall_s": time.time() - started,
        "workers": workers,
        "decoder": {
            "max_iter": MAX_ITER,
            "bp_method": BP_METHOD,
            "ms_scaling_factor": MS_SCALING_FACTOR,
            "osd_method": OSD_METHOD,
            "osd_order": OSD_ORDER,
            "dem": "undecomposed hypergraph",
        },
    }


def _resource_row(
    label: str,
    n_data: int,
    n_ancilla: int,
    data_gates: int,
    plumbing_gates: int,
    total_depth: int,
    data_only_depth: int,
) -> dict[str, Any]:
    total_qubits = n_data + n_ancilla
    full_depth = total_depth + 2
    return {
        "label": label,
        "data_qubits": n_data,
        "ancilla_qubits": n_ancilla,
        "total_physical_qubits": total_qubits,
        "data_two_qubit_gates_per_round": data_gates,
        "cat_plumbing_two_qubit_gates_per_round": plumbing_gates,
        "total_two_qubit_gates_per_round": data_gates + plumbing_gates,
        "data_cnot_only_depth": data_only_depth,
        "total_two_qubit_gate_layer_depth_per_round": total_depth,
        "total_two_qubit_depth_status": "achieved upper bound",
        "single_qubit_prep_layers_per_round": 1,
        "measurement_layers_per_round": 1,
        "full_operation_layer_depth_per_round": full_depth,
        "full_operation_depth_status": "achieved upper bound",
        "qubits_x_total_two_qubit_depth": total_qubits * total_depth,
        "qubits_x_full_operation_depth": total_qubits * full_depth,
    }


def _format_ci(row: dict[str, Any]) -> str:
    lo, hi = row["ci95_clopper_pearson"]
    if row["failures"] == 0:
        estimate = f"< {hi:.4e} (95% upper bound)"
    else:
        estimate = f"{row['ler_per_shot']:.4e}"
    return f"{estimate} [{lo:.4e}, {hi:.4e}]"


def _report_text(result: dict[str, Any]) -> str:
    rows = {x["label"]: x for x in result["benchmark"]["rows"]}
    cat = rows["PBB two-ancilla unverified cat"]
    pbb = rows["PBB one ancilla (EXP-016 schedule)"]
    css = rows["CSS Gross one ancilla (EXP-016 schedule)"]
    resources = {x["label"]: x for x in result["resource_table"]}
    cat_resource = resources["PBB two-ancilla unverified cat"]
    pbb_resource = resources["PBB one ancilla"]
    css_resource = resources["CSS Gross one ancilla"]
    search = result["schedule_search"]
    sem = result["semantics"]
    details = result["verdict_details"]
    ratio = details["cat_to_one_ancilla_pbb_ler_ratio"]
    ratio_text = (
        f"{ratio:.4f}" if ratio is not None
        else "undefined because at least one design had zero observed failures"
    )

    report = f"""{result['verdict']} — {details['reason']}

# EXP-024: two-ancilla cat extraction for `{TARGET_CODE_ID}`

## Protocol and policy

- Cat policy: **unverified two-qubit cat, no flag ancilla**.  Each of the 72
  weight-8 mixed checks uses two syndrome ancillas, with exactly four data
  couplings per half.  The Bell pair is prepared by H + CNOT, uncomputed by one
  merge CNOT, and the primary ancilla is measured in X.
- Correlated-error caveat: one cat-ancilla fault can propagate to as many as
  four data qubits.  The decoder consumes the undecomposed Stim hypergraph DEM,
  but this is not a full fault-tolerance proof.
- Noise: p={PHYSICAL_ERROR_RATE}, {ROUNDS} rounds; `DEPOLARIZE2(p)` after every
  two-qubit gate, `DEPOLARIZE1(p)` on every data or ancilla idle in every
  two-qubit layer, and preparation/measurement flips at p.  This is the
  EXP-016/mixed-stabilizer convention, extended to the cat plumbing layers.
- Decoder: BP min-sum, {MAX_ITER} iterations, scale {MS_SCALING_FACTOR},
  {OSD_METHOD} order {OSD_ORDER}, undecomposed DEM.
- Each design received {SHOTS} shots with the identical seed allocation
  {list(zip(BENCHMARK_SEEDS, SHOTS_BY_SEED))}; intervals are exact two-sided
  Clopper-Pearson 95%.

## Depth result

The verdict metric counts **all two-qubit gates**, including cat-preparation
and merge CNOTs.  The unrestricted CP-SAT model at total depth 7 returned
`{search['integrated_depth7']['status']}` in
{search['integrated_depth7']['wall_s']:.3f} s.  This search allows arbitrary
per-check 4+4 partitions and arbitrary non-translation-invariant slots, while
enforcing physical-qubit conflicts, prep < data < merge, and the controlled-
Pauli even-crossing semantics.  Thus depth <= 7 is ruled out for this specified
scheme **exactly (INFEASIBLE)**.

A data-gates-only model does reach depth {search['data_only']['depth']} with
status `{search['data_only']['status']}` (exact(OPTIMAL), since every data qubit
has degree 7).  That is not a physical round depth: putting the cat plumbing
back changes the answer.  A verified translation-invariant construction was
found at total two-qubit depth {search['constructed']['depth']} (`{search['constructed']['status']}`),
so the integrated optimum is bracketed 8..9.  The depth-8 unrestricted search
was budgeted at {INTEGRATED_DEPTH8_BUDGET_S:.0f} s and returned
`{search['integrated_depth8_budgeted']['status']}`; no exhaustiveness claim is
made for depth 8.  The achieved full operation-layer depth is
{cat_resource['full_operation_layer_depth_per_round']} after separately
bookkeeping one H/reset/prep layer and one measurement layer, just as for the
baselines.

**Conclusion on Q4:** the apparent data-only 8 -> 7 recovery is a bookkeeping
artifact.  This two-cat construction does not recover the Gross-code total
two-qubit depth of 7.

## Measurement semantics

- Cat schedule structural verifier: `{sem['cat_schedule_structural_valid']}`;
  physical-qubit conflicts = {sem['cat_schedule_collision_count']}, semantic
  parity defects = {sem['cat_schedule_semantic_parity_defects']}.
- Noiseless cat circuit: **{sem['cat_noiseless']['detector_firings']} detector
  firings in {sem['cat_noiseless']['shots']} shots** and
  {sem['cat_noiseless']['observable_flips']} observable flips.
- Cat vs one-ancilla PBB observable definitions identical:
  `{sem['pbb_observable_definitions_identical']}`
  ({sem['cat_noiseless']['num_observables']} observables; SHA-256
  `{sem['observable_definition_sha256']}`).

## Matched logical-error benchmark

| design | shots | fails | LER per shot (95% CP CI) |
|---|---:|---:|---:|
| Cat PBB | {cat['shots']} | {cat['failures']} | {_format_ci(cat)} |
| One-ancilla PBB | {pbb['shots']} | {pbb['failures']} | {_format_ci(pbb)} |
| CSS Gross | {css['shots']} | {css['failures']} | {_format_ci(css)} |

Observed cat/one-ancilla-PBB LER ratio:
{ratio_text}.  This is a point-estimate
ratio from one fixed schedule per design, not a schedule-averaged claim.

## Resource accounting (per round)

| design | ancillas | total qubits | data 2q gates | cat prep+merge CNOTs | all 2q gates | data-only depth | all-2q depth | full depth | qubits x all-2q depth | qubits x full depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cat PBB | {cat_resource['ancilla_qubits']} | {cat_resource['total_physical_qubits']} | {cat_resource['data_two_qubit_gates_per_round']} | {cat_resource['cat_plumbing_two_qubit_gates_per_round']} | {cat_resource['total_two_qubit_gates_per_round']} | {cat_resource['data_cnot_only_depth']} | {cat_resource['total_two_qubit_gate_layer_depth_per_round']} | {cat_resource['full_operation_layer_depth_per_round']} | {cat_resource['qubits_x_total_two_qubit_depth']} | {cat_resource['qubits_x_full_operation_depth']} |
| One-ancilla PBB | {pbb_resource['ancilla_qubits']} | {pbb_resource['total_physical_qubits']} | {pbb_resource['data_two_qubit_gates_per_round']} | 0 | {pbb_resource['total_two_qubit_gates_per_round']} | {pbb_resource['data_cnot_only_depth']} | {pbb_resource['total_two_qubit_gate_layer_depth_per_round']} | {pbb_resource['full_operation_layer_depth_per_round']} | {pbb_resource['qubits_x_total_two_qubit_depth']} | {pbb_resource['qubits_x_full_operation_depth']} |
| CSS Gross | {css_resource['ancilla_qubits']} | {css_resource['total_physical_qubits']} | {css_resource['data_two_qubit_gates_per_round']} | 0 | {css_resource['total_two_qubit_gates_per_round']} | {css_resource['data_cnot_only_depth']} | {css_resource['total_two_qubit_gate_layer_depth_per_round']} | {css_resource['full_operation_layer_depth_per_round']} | {css_resource['qubits_x_total_two_qubit_depth']} | {css_resource['qubits_x_full_operation_depth']} |

Relative to one-ancilla PBB, the achieved cat circuit uses
{details['cat_to_one_ancilla_pbb_space_time_ratio']:.4f}x the
qubits x all-two-qubit-depth space-time and
{details['cat_to_one_ancilla_pbb_two_qubit_gate_ratio']:.4f}x the two-qubit
gates per round.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/python experiments/exp024_cat_extraction.py --workers 6
```

Wall time: {result['wall_s']:.1f} s on {result['environment']['platform']}.
Monte Carlo was capped at {result['protocol']['workers']} workers; shared-machine
load can affect wall time but not counts or seeded samples.

## Caveats

1. The cat is unverified; no third flag ancilla or pre-interaction verification
   detector is present.  Flag-assisted variants remain untested.
2. This is one fixed schedule for each design.  EXP-016 showed material
   between-schedule LER variation, so no best/worst or schedule-averaged claim
   is licensed.
3. Depth 7 is exactly infeasible for the specified H+CNOT / 4+4 / merge-CNOT
   extraction model.  Depth 8 remains open after a budgeted search; achieved
   depth 9 is an upper bound, not an optimum claim.
4. Exact catalogue distance 12 is inherited from the verified catalogue; this
   experiment makes no new code-distance or full fault-tolerance claim.
"""
    return report

def _structural_report_text(result: dict[str, Any]) -> str:
    search = result["schedule_search"]
    semantics = result["semantics"]
    resources = {row["label"]: row for row in result["resource_table"]}
    cat = resources["PBB two-ancilla unverified cat"]
    pbb = resources["PBB one ancilla"]
    gross = resources["CSS Gross one ancilla"]
    return f"""NEGATIVE_STRUCTURAL — the specified unverified 4+4 two-ancilla scheme cannot recover Gross's seven-layer all-two-qubit round.

# EXP-024 structural result: two-ancilla cat extraction

## Scope and policy

This result concerns catalogue member `{TARGET_CODE_ID}` with the exact
catalogue polynomials persisted in the JSON artifact.  Every weight-8
perturbed check uses an **unverified** two-qubit cat: H plus a cat-preparation
CNOT, four controlled-Pauli data couplings from each half, one merge CNOT, and
X measurement of the primary ancilla.  The 72 weight-6 checks remain
one-ancilla measurements.

The verdict metric counts every two-qubit gate layer, including cat
preparation and merge.  Single-qubit preparation/H and measurement are
bookkept separately.  This policy is not fault tolerant: one fault on a cat
ancilla can fan out to four data qubits.  No flag or cat-verification detector
is present.

## Exact and budgeted depth evidence

- Integrated all-two-qubit depth 7: **{search['integrated_depth7']['status']}**,
  exact for the specified scheme.  The unrestricted model permits a separate
  arbitrary 4+4 partition for every check and non-translation-invariant data,
  preparation, and merge slots; it enforces physical-qubit disjointness,
  prep < data < merge, and the controlled-Pauli even-crossing criterion.
- Integrated depth 8: **{search['integrated_depth8_budgeted']['status']}** after
  the preregistered {search['integrated_depth8_budgeted']['time_limit_s']:.0f} s
  budget.  This is budgeted inconclusive evidence, not an infeasibility proof.
- Achieved construction: total all-two-qubit depth
  **{cat['total_two_qubit_gate_layer_depth_per_round']}**, a valid
  translation-invariant upper bound.
- Data-coupling-only model: depth **{search['data_only']['depth']}**
  ({search['data_only']['status']}).  Excluding cat plumbing creates an
  apparent 8 -> 7 improvement, but that is a bookkeeping artifact and is not
  the verdict depth.
- Integrated optimum bracket for this scheme: **[8, 9]**.

The one-ancilla PBB depth 8 label is class-free exact: weight 8 is a lower
bound and the persisted schedule is a matching witness.  The Gross depth 7
schedule and the EXP-016 schedule populations are translation-invariant; our
local artifact certifies optimality within that class, while the unrestricted
depth-6 exclusion is external (ASC, arXiv:2603.21499).

## Measurement semantics

- Structural schedule valid: `{semantics['cat_schedule_structural_valid']}`.
- Physical-qubit collisions:
  {semantics['cat_schedule_collision_count']}.
- Controlled-Pauli semantic parity defects:
  {semantics['cat_schedule_semantic_parity_defects']}.
- Noiseless sampling: **{semantics['cat_noiseless']['detector_firings']}**
  detector firings and **{semantics['cat_noiseless']['observable_flips']}**
  observable flips in {semantics['cat_noiseless']['shots']} shots.
- Cat and one-ancilla PBB observable definitions identical:
  `{semantics['pbb_observable_definitions_identical']}`
  (SHA-256 `{semantics['observable_definition_sha256']}`).

## Exact resource accounting per round

| design | data | ancillas | total qubits | data 2q | cat prep+merge | all 2q gates | data-only depth | all-2q depth | full operation depth | qubits x all-2q depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| two-ancilla PBB | {cat['data_qubits']} | {cat['ancilla_qubits']} | {cat['total_physical_qubits']} | {cat['data_two_qubit_gates_per_round']} | {cat['cat_plumbing_two_qubit_gates_per_round']} | {cat['total_two_qubit_gates_per_round']} | {cat['data_cnot_only_depth']} | {cat['total_two_qubit_gate_layer_depth_per_round']} | {cat['full_operation_layer_depth_per_round']} | {cat['qubits_x_total_two_qubit_depth']} |
| one-ancilla PBB | {pbb['data_qubits']} | {pbb['ancilla_qubits']} | {pbb['total_physical_qubits']} | {pbb['data_two_qubit_gates_per_round']} | 0 | {pbb['total_two_qubit_gates_per_round']} | {pbb['data_cnot_only_depth']} | {pbb['total_two_qubit_gate_layer_depth_per_round']} | {pbb['full_operation_layer_depth_per_round']} | {pbb['qubits_x_total_two_qubit_depth']} |
| Gross | {gross['data_qubits']} | {gross['ancilla_qubits']} | {gross['total_physical_qubits']} | {gross['data_two_qubit_gates_per_round']} | 0 | {gross['total_two_qubit_gates_per_round']} | {gross['data_cnot_only_depth']} | {gross['total_two_qubit_gate_layer_depth_per_round']} | {gross['full_operation_layer_depth_per_round']} | {gross['qubits_x_total_two_qubit_depth']} |

## Verdict and limitations

**NEGATIVE_STRUCTURAL.**  This specified scheme cannot recover Gross's
seven-layer round: depth 7 is exactly infeasible and the achieved physical
construction has depth 9.  This is not a best-multi-ancilla-scheme claim, a
logical-error claim, or a fault-tolerance proof.  Flag-assisted and verified
cat variants remain open.

The 3-way BP+OSD benchmark is not reported here.  Its v1 feasibility failure
and frozen v2 protocol are recorded in
`results/partial_runs/exp024_benchmark_feasibility_amendment.json`; the
zero-byte canonical benchmark sentinel remains untouched.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \\
NUMEXPR_NUM_THREADS=1 PYTHONPATH=src \\
.venv/bin/python experiments/exp024_cat_extraction.py --structural-only
```

Wall time: {result['wall_s']:.3f} s on
{result['environment']['platform']}.  Depth-8 status is budget dependent;
no exhaustiveness claim is made for that search.
"""


def _validate_structural_acceptance(result: dict[str, Any]) -> None:
    search = result["schedule_search"]
    semantics = result["semantics"]
    resources = {row["label"]: row for row in result["resource_table"]}
    cat = resources["PBB two-ancilla unverified cat"]
    required = {
        "target": result["catalogue_member"]["code_id"] == TARGET_CODE_ID,
        "catalogue_parameters": (
            result["catalogue_member"]["n"],
            result["catalogue_member"]["k"],
            result["catalogue_member"]["d"],
        ) == (144, 12, 12),
        "depth7_exact_infeasible": (
            search["integrated_depth7"]["status"] == "INFEASIBLE"
            and search["depth7_claim"] == "exact(INFEASIBLE)"
        ),
        "depth8_budgeted_unknown": (
            search["integrated_depth8_budgeted"]["status"] == "UNKNOWN"
        ),
        "depth9_constructed": (
            search["constructed"]["depth"] == 9
            and search["constructed"]["status"] in ("OPTIMAL", "FEASIBLE")
        ),
        "data_only_depth7": (
            search["data_only"]["depth"] == 7
            and search["data_only"]["status"] in ("OPTIMAL", "FEASIBLE")
        ),
        "schedule_valid": search["actual_schedule_verification"]["valid"],
        "noiseless_semantics": (
            semantics["cat_noiseless"]["shots"] == NOISELESS_SHOTS
            and semantics["cat_noiseless"]["detector_firings"] == 0
            and semantics["cat_noiseless"]["observable_flips"] == 0
        ),
        "observables_identical": semantics[
            "pbb_observable_definitions_identical"
        ],
        "resources": (
            cat["data_qubits"],
            cat["ancilla_qubits"],
            cat["total_physical_qubits"],
            cat["total_two_qubit_gates_per_round"],
            cat["total_two_qubit_gate_layer_depth_per_round"],
        ) == (144, 216, 360, 1152, 9),
    }
    failed = [name for name, passed in required.items() if not passed]
    result["acceptance_checks"] = required
    if failed:
        raise RuntimeError(f"structural acceptance checks failed: {failed}")


def run_structural() -> dict[str, Any]:
    benchmark_hash_before = _sha256(OUTPUT) if OUTPUT.exists() else None
    started = time.time()
    try:
        result = run(workers=1, validate_only=True)
        result["experiment"] = "EXP-024-structural"
        result["scope"] = (
            "specified unverified two-ancilla 4+4 extraction; structural "
            "depth, semantics, and resources only"
        )
        result["verdict"] = "NEGATIVE_STRUCTURAL"
        result["verdict_details"] = {
            "reason": (
                "unrestricted integrated depth 7 is exactly infeasible; "
                "the valid physical construction has depth 9 and the "
                "data-only depth-7 result is a bookkeeping artifact"
            ),
            "integrated_optimum_bracket": [8, 9],
            "logical_error_claim": False,
            "fault_tolerance_claim": False,
            "best_multi_ancilla_scheme_claim": False,
        }
        result["benchmark"] = {
            "status": "NOT_RUN_IN_STRUCTURAL_PROTOCOL",
            "canonical_benchmark_path": str(OUTPUT.relative_to(ROOT)),
            "canonical_benchmark_sha256_before": benchmark_hash_before,
        }
        result["protocol"]["protocol_kind"] = "structural_only"
        result["protocol"]["full_scope"] = True
        result["protocol"]["benchmark_shots"] = 0
        result["wall_s"] = time.time() - started
        _validate_structural_acceptance(result)
        if canonical_route(clean=True, full_coverage=True) != "canonical":
            raise RuntimeError("structural canonical route rejected full clean run")
        benchmark_hash_after = _sha256(OUTPUT) if OUTPUT.exists() else None
        result["benchmark"]["canonical_benchmark_sha256_after"] = benchmark_hash_after
        result["benchmark"]["canonical_benchmark_untouched"] = (
            benchmark_hash_before == benchmark_hash_after
        )
        if not result["benchmark"]["canonical_benchmark_untouched"]:
            raise RuntimeError("structural run touched the benchmark artifact")
        _atomic_write_json(STRUCTURAL_OUTPUT, result)
        _atomic_write_text(STRUCTURAL_REPORT, _structural_report_text(result))
        return result
    except BaseException as error:
        failure = {
            "experiment": "EXP-024-structural",
            "status": "QUARANTINED_FAILING",
            "error_type": type(error).__name__,
            "error": str(error),
            "wall_s": time.time() - started,
            "artifact_route": canonical_route(clean=False, full_coverage=False),
        }
        _atomic_write_json(STRUCTURAL_QUARANTINE, failure)
        raise


def run(workers: int, validate_only: bool = False) -> dict[str, Any]:
    if workers < 1 or workers > MAX_MONTE_CARLO_WORKERS:
        raise ValueError(f"workers must be in [1, {MAX_MONTE_CARLO_WORKERS}]")
    if sum(SHOTS_BY_SEED) != SHOTS:
        raise RuntimeError("pre-registered shot allocation does not sum to SHOTS")
    if not validate_only:
        # Fail before schedule solving or circuit construction.  The legacy
        # benchmark remains below only for forensic reproducibility.
        run_benchmark_v2_production()
    started = time.time()
    row, pbb_spec = _load_target()
    pbb_code, pbb_supports, pbb_orbits = pbb_supports_and_orbits(pbb_spec)
    css_code, css_supports, css_orbits = bb_supports_and_orbits(
        BRAVYI_BB["[[144,12,12]]"]
    )
    if [s.weight for s in pbb_supports].count(8) != 72:
        raise RuntimeError("expected 72 weight-8 perturbed checks")

    depth7 = _solve_unrestricted_integrated(
        pbb_supports,
        pbb_code.n,
        7,
        INTEGRATED_DEPTH7_TIME_LIMIT_S,
    )
    if depth7["status"] != "INFEASIBLE":
        raise RuntimeError(
            f"pre-registered depth-7 gate expected INFEASIBLE, got {depth7['status']}"
        )
    depth8 = _solve_unrestricted_integrated(
        pbb_supports,
        pbb_code.n,
        8,
        INTEGRATED_DEPTH8_BUDGET_S,
    )
    data_only = _solve_invariant(
        pbb_supports,
        pbb_code.n,
        pbb_orbits,
        DATA_ONLY_DEPTH,
        include_plumbing=False,
    )
    if data_only["status"] not in ("OPTIMAL", "FEASIBLE"):
        raise RuntimeError(f"data-only depth-7 solve failed: {data_only['status']}")

    if depth8.get("schedule") is not None:
        cat_schedule = depth8["schedule"]
        constructed = {
            key: value for key, value in depth8.items() if key != "schedule"
        }
        constructed["source"] = "unrestricted depth-8 search"
    else:
        construction = _solve_invariant(
            pbb_supports,
            pbb_code.n,
            pbb_orbits,
            CONSTRUCTIVE_DEPTH,
            include_plumbing=True,
        )
        if construction["status"] not in ("OPTIMAL", "FEASIBLE"):
            raise RuntimeError(
                f"constructive depth-{CONSTRUCTIVE_DEPTH} solve failed: "
                f"{construction['status']}"
            )
        cat_schedule = construction.pop("schedule")
        constructed = construction
        constructed["source"] = "translation-invariant fixed orbit split"

    cat_schedule_check = _verify_cat_schedule(
        pbb_supports, pbb_code.n, cat_schedule
    )
    if not cat_schedule_check["valid"]:
        raise RuntimeError(f"cat schedule verification failed: {cat_schedule_check}")

    pbb_slot, pbb_schedule_source = _exp016_slot(
        f"nonCSS-PBB [[144,12,12]] {TARGET_CODE_ID}"
    )
    css_slot, css_schedule_source = _exp016_slot("CSS-BB [[144,12,12]] Gross")
    if not verify_schedule(pbb_supports, pbb_slot)["valid"]:
        raise RuntimeError("persisted EXP-016 PBB schedule failed verification")
    if not verify_schedule(css_supports, css_slot)["valid"]:
        raise RuntimeError("persisted EXP-016 CSS schedule failed verification")

    observables_pbb = pure_z_logical_basis(pbb_code)
    observables_css = pure_z_logical_basis(css_code)
    cat0, cat_meta0 = build_cat_memory_circuit(
        pbb_code.H, observables_pbb, ROUNDS, 0.0, cat_schedule
    )
    pbb0, pbb_meta0 = build_memory_circuit(
        CircuitSpec(
            H=pbb_code.H,
            observables=observables_pbb,
            rounds=ROUNDS,
            p=0.0,
            basis="Z",
            layers=slots_to_layers(pbb_slot, pbb_schedule_source["depth"]),
        )
    )
    css0, css_meta0 = build_memory_circuit(
        CircuitSpec(
            H=css_code.H,
            observables=observables_css,
            rounds=ROUNDS,
            p=0.0,
            basis="Z",
            layers=slots_to_layers(css_slot, css_schedule_source["depth"]),
        )
    )
    cat_noiseless = _noiseless_evidence(cat0, NOISELESS_SEED)
    pbb_noiseless = _noiseless_evidence(pbb0, NOISELESS_SEED + 1)
    css_noiseless = _noiseless_evidence(css0, NOISELESS_SEED + 2)
    cat_observable_defs = _observable_definitions(cat0)
    pbb_observable_defs = _observable_definitions(pbb0)
    observable_defs_identical = cat_observable_defs == pbb_observable_defs
    observable_hash = hashlib.sha256(
        "\n".join(cat_observable_defs).encode()
    ).hexdigest()
    if cat_noiseless["detector_firings"] or cat_noiseless["observable_flips"]:
        raise RuntimeError(f"cat noiseless semantics failed: {cat_noiseless}")
    if not observable_defs_identical:
        raise RuntimeError("cat/PBB observable definitions differ")

    actual_total_depth = int(cat_schedule["depth"])
    resource_table = [
        _resource_row(
            "PBB two-ancilla unverified cat",
            pbb_code.n,
            cat_meta0["n_ancilla"],
            cat_meta0["data_two_qubit_gates_per_round"],
            cat_meta0["cat_prep_cnot_per_round"]
            + cat_meta0["cat_merge_cnot_per_round"],
            actual_total_depth,
            DATA_ONLY_DEPTH,
        ),
        _resource_row(
            "PBB one ancilla",
            pbb_code.n,
            pbb_meta0["n_ancilla"],
            pbb_meta0["two_qubit_gates_per_round"],
            0,
            pbb_schedule_source["depth"],
            pbb_schedule_source["depth"],
        ),
        _resource_row(
            "CSS Gross one ancilla",
            css_code.n,
            css_meta0["n_ancilla"],
            css_meta0["two_qubit_gates_per_round"],
            0,
            css_schedule_source["depth"],
            css_schedule_source["depth"],
        ),
    ]
    # PBB depth 8 is class-free exact: check weight gives the lower bound and
    # this TI witness supplies the matching upper bound.  Gross depth 7 is
    # certified minimal only inside the enumerated TI class by our artifacts;
    # ASC arXiv:2603.21499 supplies the unrestricted depth-6 exclusion.
    resource_table[1]["total_two_qubit_depth_status"] = (
        "exact (class-free; weight-8 lower bound plus valid witness)"
    )
    resource_table[2]["total_two_qubit_depth_status"] = (
        "exact within translation-invariant class; unrestricted depth-6 "
        "exclusion is external (ASC arXiv:2603.21499)"
    )
    if actual_total_depth == 8:
        resource_table[0]["total_two_qubit_depth_status"] = (
            "exact within specified unverified 4+4 cat-extraction model"
        )

    protocol = {
        "p": PHYSICAL_ERROR_RATE,
        "rounds": ROUNDS,
        "shots_per_design": SHOTS,
        "benchmark_seeds": list(BENCHMARK_SEEDS),
        "shots_by_seed": list(SHOTS_BY_SEED),
        "noiseless_shots": NOISELESS_SHOTS,
        "noiseless_seed": NOISELESS_SEED,
        "workers": workers,
        "max_workers_cap": MAX_MONTE_CARLO_WORKERS,
        "cat_policy": CAT_POLICY,
        "baseline_schedule_class": (
            "translation-invariant (orbit); selected from EXP-016's complete "
            "enumeration within that class"
        ),
        "cat_split": "exactly 4+4 for every weight-8 check",
        "depth_convention": (
            "all two-qubit gates count; cat prep/merge CNOTs may overlap data "
            "couplings on disjoint physical qubits; single-qubit prep/H and "
            "measurement are separately bookkept"
        ),
        "noise_model": (
            "DEPOLARIZE2(p) after every two-qubit gate; DEPOLARIZE1(p) on "
            "every idle data and ancilla in each two-qubit layer; prep and "
            "measurement flips p"
        ),
        "decoder": {
            "max_iter": MAX_ITER,
            "bp_method": BP_METHOD,
            "ms_scaling_factor": MS_SCALING_FACTOR,
            "osd_method": OSD_METHOD,
            "osd_order": OSD_ORDER,
            "dem": "undecomposed hypergraph",
        },
        "confidence_interval": "two-sided Clopper-Pearson 95%",
    }

    result: dict[str, Any] = {
        "experiment": "EXP-024",
        "environment": _versions(),
        "protocol": protocol,
        "catalogue_member": {
            key: row[key]
            for key in (
                "code_id",
                "ell",
                "m",
                "A_terms",
                "B_terms",
                "C_terms",
                "D_terms",
                "n",
                "k",
                "d",
                "d_is_exact",
            )
            if key in row
        },
        "schedule_search": {
            "integrated_depth7": {
                key: value for key, value in depth7.items() if key != "schedule"
            },
            "integrated_depth8_budgeted": {
                key: value for key, value in depth8.items() if key != "schedule"
            },
            "data_only": {
                key: value for key, value in data_only.items() if key != "schedule"
            },
            "constructed": constructed,
            "integrated_optimum_bracket": [8, actual_total_depth],
            "depth7_claim": "exact(INFEASIBLE)",
            "data_only_depth_claim": "exact(OPTIMAL)",
            "achieved_total_depth_claim": (
                "exact(OPTIMAL)" if actual_total_depth == 8 else "upper bound"
            ),
            "actual_schedule": _serialise_schedule(cat_schedule),
            "actual_schedule_verification": cat_schedule_check,
            "pbb_baseline_schedule": pbb_schedule_source,
            "css_baseline_schedule": css_schedule_source,
        },
        "semantics": {
            "cat_schedule_structural_valid": cat_schedule_check["valid"],
            "cat_schedule_collision_count": len(
                cat_schedule_check["collision_examples"]
            ),
            "cat_schedule_semantic_parity_defects": cat_schedule_check[
                "semantic_parity_defects"
            ],
            "cat_noiseless": cat_noiseless,
            "pbb_noiseless": pbb_noiseless,
            "css_noiseless": css_noiseless,
            "pbb_observable_definitions_identical": observable_defs_identical,
            "observable_definitions": cat_observable_defs,
            "observable_definition_sha256": observable_hash,
            "observable_matrix_sha256": hashlib.sha256(
                np.ascontiguousarray(observables_pbb).tobytes()
            ).hexdigest(),
        },
        "resource_table": resource_table,
        "benchmark": {"rows": []},
        "verdict": "INCONCLUSIVE" if validate_only else None,
        "verdict_details": {
            "reason": "validation-only run; Monte Carlo not executed"
            if validate_only
            else None
        },
        "wall_s": None,
    }

    if validate_only:
        result["wall_s"] = time.time() - started
        route = canonical_route(clean=True, full_coverage=False)
        partial = ROOT / "results" / route / "exp024_cat_extraction-validate_only.json"
        result["artifact_route"] = route
        result["artifact_path"] = str(partial.relative_to(ROOT))
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.write_text(json.dumps(result, indent=2))
        return result

    cat_circuit, _ = build_cat_memory_circuit(
        pbb_code.H,
        observables_pbb,
        ROUNDS,
        PHYSICAL_ERROR_RATE,
        cat_schedule,
    )
    pbb_circuit, _ = build_memory_circuit(
        CircuitSpec(
            H=pbb_code.H,
            observables=observables_pbb,
            rounds=ROUNDS,
            p=PHYSICAL_ERROR_RATE,
            basis="Z",
            layers=slots_to_layers(pbb_slot, pbb_schedule_source["depth"]),
        )
    )
    css_circuit, _ = build_memory_circuit(
        CircuitSpec(
            H=css_code.H,
            observables=observables_css,
            rounds=ROUNDS,
            p=PHYSICAL_ERROR_RATE,
            basis="Z",
            layers=slots_to_layers(css_slot, css_schedule_source["depth"]),
        )
    )

    benchmark_rows = [
        _benchmark("PBB two-ancilla unverified cat", cat_circuit, workers),
        _benchmark("PBB one ancilla (EXP-016 schedule)", pbb_circuit, workers),
        _benchmark("CSS Gross one ancilla (EXP-016 schedule)", css_circuit, workers),
    ]
    result["benchmark"]["rows"] = benchmark_rows
    by_label = {x["label"]: x for x in benchmark_rows}
    cat_bench = by_label["PBB two-ancilla unverified cat"]
    pbb_bench = by_label["PBB one ancilla (EXP-016 schedule)"]
    resource_by_label = {x["label"]: x for x in resource_table}
    cat_resource = resource_by_label["PBB two-ancilla unverified cat"]
    pbb_resource = resource_by_label["PBB one ancilla"]

    semantics_ok = (
        cat_noiseless["detector_firings"] == 0
        and cat_noiseless["observable_flips"] == 0
        and observable_defs_identical
        and cat_schedule_check["valid"]
    )
    depth_recovered = actual_total_depth <= 7
    cat_rate = cat_bench["ler_per_shot"]
    pbb_rate = pbb_bench["ler_per_shot"]
    ler_not_worse = (
        cat_rate is not None and pbb_rate is not None and cat_rate <= pbb_rate
    )
    space_time_not_worse = (
        cat_resource["qubits_x_total_two_qubit_depth"]
        <= pbb_resource["qubits_x_total_two_qubit_depth"]
    )
    if not semantics_ok:
        verdict = "INCONCLUSIVE"
        reason = "measurement semantics could not be verified"
    elif depth_recovered and ler_not_worse and space_time_not_worse:
        verdict = "BREAKTHROUGH_CANDIDATE"
        reason = "total two-qubit depth, LER, and space-time all match or improve on one-ancilla PBB"
    elif depth_recovered:
        verdict = "POSITIVE"
        reason = "total two-qubit depth recovered, but LER or space-time is worse"
    else:
        verdict = "NEGATIVE"
        reason = (
            "unrestricted CP-SAT proves total two-qubit depth 7 infeasible; "
            f"the verified construction uses depth {actual_total_depth}, so the "
            "data-only depth-7 appearance is a bookkeeping artifact"
        )

    result["verdict"] = verdict
    result["verdict_details"] = {
        "reason": reason,
        "semantics_verified": semantics_ok,
        "total_two_qubit_depth_recovered": depth_recovered,
        "ler_not_worse_than_one_ancilla_pbb_point_estimate": ler_not_worse,
        "space_time_not_worse_than_one_ancilla_pbb": space_time_not_worse,
        "cat_to_one_ancilla_pbb_ler_ratio": (
            cat_rate / pbb_rate
            if cat_rate is not None and pbb_rate not in (None, 0.0)
            else None
        ),
        "cat_to_one_ancilla_pbb_space_time_ratio": (
            cat_resource["qubits_x_total_two_qubit_depth"]
            / pbb_resource["qubits_x_total_two_qubit_depth"]
        ),
        "cat_to_one_ancilla_pbb_two_qubit_gate_ratio": (
            cat_resource["total_two_qubit_gates_per_round"]
            / pbb_resource["total_two_qubit_gates_per_round"]
        ),
    }
    result["wall_s"] = time.time() - started
    route = canonical_route(clean=True, full_coverage=True)
    if route != "canonical":
        raise RuntimeError(f"full clean run unexpectedly routed to {route}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    output_tmp.write_text(json.dumps(result, indent=2))
    output_tmp.replace(OUTPUT)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    report_tmp = REPORT.with_suffix(REPORT.suffix + ".tmp")
    report_tmp.write_text(_report_text(result))
    report_tmp.replace(REPORT)
    return result


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=MAX_MONTE_CARLO_WORKERS)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="solve schedules and run 4096-shot noiseless checks, but skip Monte Carlo",
    )
    parser.add_argument(
        "--structural-only",
        action="store_true",
        help=(
            "run the complete structural protocol and atomically write the "
            "separate structural artifact/report; never touch benchmark output"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    if args.structural_only and args.validate_only:
        raise SystemExit("--structural-only and --validate-only are mutually exclusive")
    if args.structural_only:
        outcome = run_structural()
        print(
            f"{outcome['verdict']}: wrote {STRUCTURAL_OUTPUT.relative_to(ROOT)} "
            f"and {STRUCTURAL_REPORT.relative_to(ROOT)} in "
            f"{outcome['wall_s']:.1f}s"
        )
    else:
        outcome = run(workers=args.workers, validate_only=args.validate_only)
        if args.validate_only:
            print(json.dumps(outcome, indent=2))
        else:
            print(
                f"{outcome['verdict']}: wrote {OUTPUT.relative_to(ROOT)} and "
                f"{REPORT.relative_to(ROOT)} in {outcome['wall_s']:.1f}s"
            )
