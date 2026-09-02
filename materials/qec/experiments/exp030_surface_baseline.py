"""EXP-030: real rotated-surface-code circuit baseline for the [[144,12,12]] codes.

The surface-code points are simulated with Stim's generated circuits and decoded
with PyMatching.  Gross/PBB results are imported, not rerun, from EXP-013 and
EXP-016.  The generated surface-code noise is deliberately described as a
mapped baseline because its noise locations differ from the custom qLDPC
circuit convention.
"""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import pymatching
from scipy.stats import beta
import stim

# ---------------------------------------------------------------------------
# Pre-registered protocol constants.  Do not change these after collecting data.
# ---------------------------------------------------------------------------
EXPERIMENT = "EXP-030"
DISTANCE = 12
ROUNDS = 12
PHYSICAL_ERROR_RATES = (0.001, 0.002, 0.003)
MEMORY_TASKS = (
    ("Z", "surface_code:rotated_memory_z"),
    ("X", "surface_code:rotated_memory_x"),
)
SHOTS_PER_POINT = 200_000
BATCH_SHOTS = 20_000
BASE_SEED = 20_260_812
CONFIDENCE = 0.95
MATCHED_LOGICAL_QUBITS = 12
VERDICT = "POSITIVE"
GENERATOR_NOISE_ARGUMENTS = {
    "after_clifford_depolarization": "p",
    "before_round_data_depolarization": "p",
    "before_measure_flip_probability": "p",
    "after_reset_flip_probability": "p",
}
POINT_ORDER = tuple(
    (p, basis, task)
    for p in PHYSICAL_ERROR_RATES
    for basis, task in MEMORY_TASKS
)
POINT_SEEDS = tuple(BASE_SEED + i for i in range(len(POINT_ORDER)))

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "processed" / "exp030_surface_baseline.json"
REPORT = ROOT / "notes" / "agent_reports" / "exp030_surface_baseline.md"
EXP013_SOURCE = ROOT / "results" / "raw" / "exp013_isolated_latency.json"
EXP016_SOURCE = ROOT / "results" / "raw" / "exp016_schedule_controlled.json"

CONVENTION_DELTAS = (
    {
        "topic": "data idle noise",
        "stim_generated": (
            "before_round_data_depolarization=p applies one DEPOLARIZE1(p) to every "
            "data qubit at the start of each syndrome round, independent of whether "
            "that qubit is idle in a particular CX layer"
        ),
        "custom_qldpc": (
            "DEPOLARIZE1(p) is applied separately at every two-qubit schedule layer "
            "to each data qubit idle in that layer; a data qubit participating in the "
            "layer instead receives the two-qubit gate's DEPOLARIZE2(p)"
        ),
    },
    {
        "topic": "ancilla idle noise",
        "stim_generated": (
            "there is no general ancilla-idle channel; ancillas receive noise after "
            "reset, before measurement, after an H when present, and after a CX when active"
        ),
        "custom_qldpc": (
            "every ancilla idle in each two-qubit schedule layer receives "
            "DEPOLARIZE1(p)"
        ),
    },
    {
        "topic": "single-qubit Clifford noise",
        "stim_generated": (
            "after_clifford_depolarization=p adds DEPOLARIZE1(p) after both H layers "
            "on the X-check ancillas"
        ),
        "custom_qldpc": (
            "ancillas are prepared and measured directly with RX/MX and no separate "
            "post-H noise location exists"
        ),
    },
    {
        "topic": "two-qubit operations",
        "stim_generated": (
            "a four-layer nearest-neighbour CX schedule is used and DEPOLARIZE2(p) "
            "follows every CX"
        ),
        "custom_qldpc": (
            "seven (Gross) or eight (PBB) controlled-Pauli layers are used; "
            "DEPOLARIZE2(p) follows every CX/CY/CZ operation"
        ),
    },
    {
        "topic": "ancilla reset, basis rotation, and measurement ordering",
        "stim_generated": (
            "ancillas are reset with R; X-check ancillas are rotated by H; all "
            "ancillas are measured and reset together with MR.  The measurement flip "
            "is X_ERROR(p) immediately before MR, and the next reset flip is "
            "X_ERROR(p) immediately after MR"
        ),
        "custom_qldpc": (
            "all check ancillas are freshly prepared with RX then Z_ERROR(p), interact "
            "as controls, and are measured separately with Z_ERROR(p) then MX; the next "
            "round performs a separate RX instead of reusing an MR reset"
        ),
    },
    {
        "topic": "round-boundary and terminal ordering",
        "stim_generated": (
            "the data-round DEPOLARIZE1 follows the leading TICK and the reset error; "
            "MR also resets ancillas after the final measured round and therefore emits "
            "a terminal after-reset error on ancillas that are never used again"
        ),
        "custom_qldpc": (
            "there is no once-per-round data channel or combined MR boundary; idle "
            "channels occur inside interaction layers, and no unused reset follows the "
            "last ancilla measurement"
        ),
    },
    {
        "topic": "final data measurement",
        "stim_generated": (
            "before_measure_flip_probability=p inserts the basis-opposite Pauli just "
            "before the final M/MX data readout"
        ),
        "custom_qldpc": (
            "the same type of basis-opposite probability-p flip is inserted before "
            "the final data readout; this location is aligned even though the preceding "
            "round noise is not"
        ),
    },
    {
        "topic": "decoder and DEM mapping",
        "stim_generated": (
            "the DEM is graphlike-decomposed and decoded by PyMatching"
        ),
        "custom_qldpc": (
            "EXP-016 used BP+OSD on an undecomposed hypergraph DEM"
        ),
    },
    {
        "topic": "logical-block mapping",
        "stim_generated": (
            "each circuit contains one logical qubit; twelve statistically independent "
            "patches are mapped to k=12 using 1-(1-r)^12 rather than simulated jointly"
        ),
        "custom_qldpc": (
            "each sampled shot contains the actual twelve-logical-qubit code block and "
            "fails when any decoded logical observable is wrong"
        ),
    },
    {
        "topic": "depth metric",
        "stim_generated": (
            "depth/round is the literal TICK count, including the two H boundaries"
        ),
        "custom_qldpc": (
            "the imported depth is the number of two-qubit controlled-Pauli schedule "
            "layers and excludes reset/measurement basis-operation slices"
        ),
    },
)


def clopper_pearson(
    failures: int, shots: int, confidence: float = CONFIDENCE
) -> tuple[float, float]:
    """Exact equal-tailed binomial confidence interval (copied from EXP-010)."""
    if shots <= 0 or failures < 0 or failures > shots:
        raise ValueError("require 0 <= failures <= shots and shots > 0")
    alpha = 1.0 - confidence
    lower = 0.0 if failures == 0 else float(
        beta.ppf(alpha / 2, failures, shots - failures + 1)
    )
    upper = 1.0 if failures == shots else float(
        beta.ppf(1 - alpha / 2, failures + 1, shots - failures)
    )
    return lower, upper


def per_round_rate(per_shot_rate: float, rounds: int) -> float:
    """Copy of the exact EXP-010 block-to-round convention."""
    return float(-np.expm1(np.log1p(-per_shot_rate) / rounds)) if per_shot_rate < 1 else 1.0


def independent_block_rate(rate: float, copies: int) -> float:
    """Probability that at least one of ``copies`` independent blocks fails."""
    return float(-np.expm1(copies * np.log1p(-rate))) if rate < 1 else 1.0


def generated_circuit(task: str, p: float) -> stim.Circuit:
    return stim.Circuit.generated(
        task,
        distance=DISTANCE,
        rounds=ROUNDS,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def circuit_resources(circuit: stim.Circuit) -> dict[str, Any]:
    """Count active physical qubits and per-round gates from the circuit itself."""
    flattened = list(circuit.flattened())
    coordinates = circuit.get_final_qubit_coordinates()
    active_qubits = set(coordinates)

    final_data: set[int] | None = None
    for instruction in flattened:
        if instruction.name in {"M", "MX"}:
            final_data = {
                int(target.value)
                for target in instruction.targets_copy()
                if target.is_qubit_target
            }
    if final_data is None:
        raise AssertionError("generated circuit has no final data measurement")
    if not final_data <= active_qubits:
        raise AssertionError("final measurement contains a qubit without QUBIT_COORDS")
    ancillas = active_qubits - final_data

    ticks = 0
    cx_pairs = 0
    depolarize2_pairs = 0
    h_targets = 0
    mr_targets = 0
    d1_data_targets = 0
    d1_ancilla_targets = 0
    unexpected_two_qubit_gates: set[str] = set()
    known_two_qubit_gates = {"CX", "CY", "CZ", "SWAP", "ISWAP", "SQRT_XX", "SQRT_YY", "SQRT_ZZ"}
    for instruction in flattened:
        name = instruction.name
        qubits = [
            int(target.value)
            for target in instruction.targets_copy()
            if target.is_qubit_target
        ]
        if name == "TICK":
            ticks += 1
        elif name == "CX":
            cx_pairs += len(qubits) // 2
        elif name in known_two_qubit_gates:
            unexpected_two_qubit_gates.add(name)
        elif name == "DEPOLARIZE2":
            depolarize2_pairs += len(qubits) // 2
        elif name == "H":
            h_targets += len(qubits)
        elif name == "MR":
            mr_targets += len(qubits)
        elif name == "DEPOLARIZE1":
            d1_data_targets += sum(q in final_data for q in qubits)
            d1_ancilla_targets += sum(q in ancillas for q in qubits)

    if unexpected_two_qubit_gates:
        raise AssertionError(f"unexpected two-qubit gates: {sorted(unexpected_two_qubit_gates)}")
    for name, count in {
        "TICK": ticks,
        "CX": cx_pairs,
        "DEPOLARIZE2": depolarize2_pairs,
        "H": h_targets,
        "MR": mr_targets,
        "data DEPOLARIZE1": d1_data_targets,
        "ancilla DEPOLARIZE1": d1_ancilla_targets,
    }.items():
        if count % ROUNDS:
            raise AssertionError(f"{name} count {count} is not divisible by {ROUNDS}")

    data_count = len(final_data)
    ancilla_count = len(ancillas)
    active_count = len(active_qubits)
    if data_count + ancilla_count != active_count:
        raise AssertionError("active-qubit partition is inconsistent")
    return {
        "per_patch": {
            "data_qubits": data_count,
            "ancilla_qubits": ancilla_count,
            "total_active_qubits": active_count,
            "stim_qubit_index_span": circuit.num_qubits,
            "depth_ticks_per_round": ticks // ROUNDS,
            "depth_status": (
                "exact(OPTIMAL) for counting the fixed generated circuit: literal "
                f"{ticks} TICKs/{ROUNDS} rounds; not a global schedule-optimization claim"
            ),
            "two_qubit_gates_per_round": cx_pairs // ROUNDS,
            "two_qubit_gate_count_status": "exact circuit count",
        },
        "matched_k12_parallel_patches": {
            "logical_qubits": MATCHED_LOGICAL_QUBITS,
            "data_qubits": MATCHED_LOGICAL_QUBITS * data_count,
            "ancilla_qubits": MATCHED_LOGICAL_QUBITS * ancilla_count,
            "total_active_qubits": MATCHED_LOGICAL_QUBITS * active_count,
            "depth_ticks_per_round": ticks // ROUNDS,
            "two_qubit_gates_per_round": MATCHED_LOGICAL_QUBITS * (cx_pairs // ROUNDS),
            "scaling_assumption": "12 independent patches operated in parallel",
        },
        "noise_instruction_audit_per_round": {
            "data_depolarize1_target_locations": d1_data_targets // ROUNDS,
            "ancilla_depolarize1_target_locations": d1_ancilla_targets // ROUNDS,
            "h_target_locations": h_targets // ROUNDS,
            "cx_target_pairs": cx_pairs // ROUNDS,
            "depolarize2_target_pairs": depolarize2_pairs // ROUNDS,
            "mr_target_locations": mr_targets // ROUNDS,
        },
    }


def simulate_point(index: int, p: float, basis: str, task: str) -> dict[str, Any]:
    seed = POINT_SEEDS[index]
    started = time.perf_counter()
    circuit = generated_circuit(task, p)
    resources = circuit_resources(circuit)
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler(seed=seed)

    failures = 0
    completed = 0
    while completed < SHOTS_PER_POINT:
        batch = min(BATCH_SHOTS, SHOTS_PER_POINT - completed)
        detector_samples, observable_samples = sampler.sample(
            batch, separate_observables=True
        )
        predicted_observables = matching.decode_batch(detector_samples)
        failures += int(
            np.count_nonzero(
                np.any(predicted_observables != observable_samples, axis=1)
            )
        )
        completed += batch

    empirical = failures / completed
    ci_lower, ci_upper = clopper_pearson(failures, completed)
    if failures == 0:
        reported_single_block = ci_upper
        report_type = "clopper_pearson_95_upper_bound"
    else:
        reported_single_block = empirical
        report_type = "empirical_point_estimate"

    mapped_block = independent_block_rate(
        reported_single_block, MATCHED_LOGICAL_QUBITS
    )
    mapped_block_ci = {
        "lower": independent_block_rate(ci_lower, MATCHED_LOGICAL_QUBITS),
        "upper": independent_block_rate(ci_upper, MATCHED_LOGICAL_QUBITS),
    }
    mapped_round = per_round_rate(mapped_block, ROUNDS)
    mapped_round_ci = {
        "lower": per_round_rate(mapped_block_ci["lower"], ROUNDS),
        "upper": per_round_rate(mapped_block_ci["upper"], ROUNDS),
    }
    return {
        "point_index": index,
        "p": p,
        "basis": basis,
        "stim_task": task,
        "seed": seed,
        "shots": completed,
        "logical_failures": failures,
        "rate_report_type": report_type,
        "single_patch_12_round_block": {
            "empirical_ler": empirical,
            "reported_ler": reported_single_block,
            "clopper_pearson_95_ci": {"lower": ci_lower, "upper": ci_upper},
        },
        "matched_k12_12_round_block": {
            "reported_ler": mapped_block,
            "clopper_pearson_95_ci": mapped_block_ci,
            "mapping": "1-(1-single_patch_block_ler)^12",
        },
        "matched_k12_ler_per_round": {
            "reported_ler": mapped_round,
            "clopper_pearson_95_ci": mapped_round_ci,
            "conversion": "-expm1(log1p(-block_ler)/12)",
        },
        "detector_error_model": {
            "decompose_errors": True,
            "detectors": dem.num_detectors,
            "observables": dem.num_observables,
            "error_instructions": dem.num_errors,
            "matching_edges": matching.num_edges,
        },
        "resources": resources,
        "wall_s": time.perf_counter() - started,
    }


def imported_qldpc_results() -> list[dict[str, Any]]:
    exp013 = json.loads(EXP013_SOURCE.read_text())
    exp016 = json.loads(EXP016_SOURCE.read_text())
    if exp013["env"]["rounds"] != ROUNDS or exp013["env"]["p"] != 0.002:
        raise AssertionError("EXP-013 source is not the expected 12-round p=0.002 artifact")
    if exp016["env"]["rounds"] != ROUNDS or exp016["env"]["p"] != 0.002:
        raise AssertionError("EXP-016 source is not the expected 12-round p=0.002 artifact")

    resource_by_label = {entry["label"]: entry for entry in exp013["codes"]}
    imported: list[dict[str, Any]] = []
    for code in exp016["codes"]:
        label = code["label"]
        source_resource = resource_by_label[label]
        qubit_counts: set[tuple[int, int]] = set()
        for run in code["runs"]:
            checks = {int(item[0]) for item in run["slot_map"]}
            data = {int(item[1]) for item in run["slot_map"]}
            qubit_counts.add((len(data), len(checks)))
        if len(qubit_counts) != 1:
            raise AssertionError(f"inconsistent persisted slot-map qubit counts for {label}")
        data_qubits, ancilla_qubits = qubit_counts.pop()

        block_ler = float(code["ler_mean_over_schedules"])
        block_ci = [float(x) for x in code["ler_mean_ci95"]]
        total_shots = sum(int(run["shots"]) for run in code["runs"])
        total_failures = sum(int(run["failures"]) for run in code["runs"])
        if not np.isclose(total_failures / total_shots, block_ler):
            raise AssertionError(f"schedule mean and pooled point estimate differ for {label}")
        per_round_ci = {
            "lower": per_round_rate(block_ci[0], ROUNDS),
            "upper": per_round_rate(block_ci[1], ROUNDS),
        }
        is_gross = "Gross" in label
        imported.append(
            {
                "code": "Gross [[144,12,12]]" if is_gross else "PBB 12_6_0193 [[144,12,12]]",
                "source_label": label,
                "basis": "Z",
                "p": 0.002,
                "rounds": ROUNDS,
                "distance": DISTANCE,
                "distance_status": (
                    "independently certified exact (two-sided solver proof plus witness)"
                    if is_gross else
                    "catalogue metadata marks exact; independent full symplectic "
                    "distance certificate pending"
                ),
                "resources": {
                    "data_qubits": data_qubits,
                    "ancilla_qubits": ancilla_qubits,
                    "total_qubits": data_qubits + ancilla_qubits,
                    "depth_per_round": int(source_resource["depth"]),
                    "depth_status": (
                        "exact within translation-invariant class; unrestricted "
                        "depth-6 exclusion is external (ASC arXiv:2603.21499)"
                        if is_gross else
                        "class-free exact: weight-8 lower bound plus valid depth-8 "
                        "witness; sampled schedule population is translation-invariant"
                    ),
                    "two_qubit_gates_per_round": int(
                        source_resource["two_qubit_gates_per_round"]
                    ),
                },
                "schedule_sample": {
                    "schedule_class": "translation-invariant (orbit)",
                    "sampling_seed": int(exp016["env"]["sampling_rng_seed"]),
                    "schedules": len(code["runs"]),
                    "total_valid_minimum_depth_schedules_within_class": int(
                        code["total_valid_schedules"]
                    ),
                    "enumeration_complete_within_class": bool(code["enumeration_complete"]),
                    "shots": total_shots,
                    "failures": total_failures,
                    "per_schedule": [
                        {
                            "schedule_index": int(run["schedule_index"]),
                            "shots": int(run["shots"]),
                            "failures": int(run["failures"]),
                            "block_ler": float(run["ler_per_shot"]),
                            "bonferroni_simultaneous_ci": {
                                "lower": float(run["ci_simultaneous_lo"]),
                                "upper": float(run["ci_simultaneous_hi"]),
                            },
                            "wall_s": float(run["wall_s"]),
                        }
                        for run in code["runs"]
                    ],
                    "source_decode_wall_s": sum(
                        float(run["wall_s"]) for run in code["runs"]
                    ),
                },
                "block_ler": {
                    "estimate": block_ler,
                    "ci95": {"lower": block_ci[0], "upper": block_ci[1]},
                    "ci_kind": (
                        "two-sided t interval over five uniformly sampled schedule "
                        "point estimates; includes measured between-schedule variation"
                    ),
                },
                "ler_per_round": {
                    "estimate": per_round_rate(block_ler, ROUNDS),
                    "ci95": per_round_ci,
                    "conversion": "-expm1(log1p(-block_ler)/12)",
                },
                "sources": {
                    "resources": "results/raw/exp013_isolated_latency.json",
                    "qubit_counts_and_statistics": "results/raw/exp016_schedule_controlled.json",
                },
            }
        )
    if len(imported) != 2:
        raise AssertionError("expected exactly Gross and PBB in EXP-016")
    return imported


def rate_text(estimate: float, interval: dict[str, float], upper_bound: bool = False) -> str:
    if upper_bound:
        return f"<={estimate:.6g} (CP95 upper; [{interval['lower']:.6g}, {interval['upper']:.6g}])"
    return f"{estimate:.6g} [{interval['lower']:.6g}, {interval['upper']:.6g}]"


def make_report(result: dict[str, Any]) -> str:
    qldpc = {entry["code"]: entry for entry in result["qldpc_imported"]}
    points = {(entry["p"], entry["basis"]): entry for entry in result["surface"]["points"]}
    surface_resources = result["surface"]["resources"]["matched_k12_parallel_patches"]
    gross = qldpc["Gross [[144,12,12]]"]
    pbb = qldpc["PBB 12_6_0193 [[144,12,12]]"]

    lines = [
        "POSITIVE — Both real Stim/PyMatching d=12 memory circuits were sampled at all three pre-registered noise points, closing the p=0.002 matched-k=12 Gross/PBB/surface table with the noise-convention deltas explicit.",
        "",
        "# EXP-030 real circuit-level rotated-surface-code baseline",
        "",
        "## Verdict and headline",
        "",
    ]
    z002 = points[(0.002, "Z")]
    x002 = points[(0.002, "X")]
    z_round = z002["matched_k12_ler_per_round"]
    x_round = x002["matched_k12_ler_per_round"]
    gross_total = gross["resources"]["total_qubits"]
    surface_total = surface_resources["total_active_qubits"]
    lines.extend(
        [
            f"At matched k=12, twelve parallel d=12 rotated patches use **{surface_total} active physical qubits** versus **{gross_total}** for Gross: **{surface_total - gross_total} more ({surface_total / gross_total:.3f}x)**.",
            (
                "At p=0.002 their mapped k=12 LER/round is "
                f"Z-memory {rate_text(z_round['reported_ler'], z_round['clopper_pearson_95_ci'], z002['logical_failures'] == 0)} "
                f"and X-memory {rate_text(x_round['reported_ler'], x_round['clopper_pearson_95_ci'], x002['logical_failures'] == 0)}."
            ),
            "This is a **mapped baseline, not an identical-noise comparison**; the complete placement differences are listed below.",
            "",
            "## Complete three-way Pareto table at p=0.002, 12 rounds, k=12",
            "",
            "Surface resources and two-qubit gate counts are multiplied by 12 independent patches; patch depth is unchanged because the patches run in parallel. Gross/PBB statistics are imported without rerunning them.",
            "",
            "| code / memory | data | ancilla | total | depth/round | 2q gates/round | shots / fails | 12-round block LER (95% CI) | LER/round (95% CI) | CI type |",
            "|---|---:|---:|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for entry in (gross, pbb):
        r = entry["resources"]
        s = entry["schedule_sample"]
        lines.append(
            f"| {entry['code']} / Z | {r['data_qubits']} | {r['ancilla_qubits']} | {r['total_qubits']} | "
            f"{r['depth_per_round']}* | {r['two_qubit_gates_per_round']} | "
            f"{s['shots']} / {s['failures']} | "
            f"{entry['block_ler']['estimate']:.6g} [{entry['block_ler']['ci95']['lower']:.6g}, {entry['block_ler']['ci95']['upper']:.6g}] | "
            f"{entry['ler_per_round']['estimate']:.6g} [{entry['ler_per_round']['ci95']['lower']:.6g}, {entry['ler_per_round']['ci95']['upper']:.6g}] | schedule-mean t CI |"
        )
    for basis, entry in (("Z", z002), ("X", x002)):
        block = entry["matched_k12_12_round_block"]
        per_round = entry["matched_k12_ler_per_round"]
        upper = entry["logical_failures"] == 0
        lines.append(
            f"| 12 x rotated surface d=12 / {basis} | {surface_resources['data_qubits']} | "
            f"{surface_resources['ancilla_qubits']} | {surface_resources['total_active_qubits']} | "
            f"{surface_resources['depth_ticks_per_round']} exact fixed-circuit TICK count* | "
            f"{surface_resources['two_qubit_gates_per_round']} | {entry['shots']} / {entry['logical_failures']} | "
            f"{rate_text(block['reported_ler'], block['clopper_pearson_95_ci'], upper)} | "
            f"{rate_text(per_round['reported_ler'], per_round['clopper_pearson_95_ci'], upper)} | CP95 on one patch, monotone-mapped x12 |"
        )
    lines.extend(
        [
            "",
            "* Gross depth 7 is exact within the translation-invariant class; ASC arXiv:2603.21499 supplies the unrestricted depth-6 exclusion. PBB depth 8 is class-free exact here (weight-8 lower bound plus a valid depth-8 witness), while the five-schedule statistical population remains translation-invariant. Surface depth is the literal fixed-circuit count (84 TICKs / 12), not a global optimum. qLDPC depth counts only two-qubit layers; surface TICK depth includes H boundaries.",
            "",
            "## Surface sweep: every simulated point",
            "",
            "| p | memory | seed | shots | fails | one-patch 12-round block | mapped k=12 12-round block LER (CP95) | mapped k=12 LER/round (CP95) | wall (s) |",
            "|---:|:---:|---:|---:|---:|---|---|---|---:|",
        ]
    )
    for entry in result["surface"]["points"]:
        single = entry["single_patch_12_round_block"]
        block = entry["matched_k12_12_round_block"]
        per_round = entry["matched_k12_ler_per_round"]
        upper = entry["logical_failures"] == 0
        lines.append(
            f"| {entry['p']:.3f} | {entry['basis']} | {entry['seed']} | {entry['shots']} | "
            f"{entry['logical_failures']} | {rate_text(single['reported_ler'], single['clopper_pearson_95_ci'], upper)} | "
            f"{rate_text(block['reported_ler'], block['clopper_pearson_95_ci'], upper)} | "
            f"{rate_text(per_round['reported_ler'], per_round['clopper_pearson_95_ci'], upper)} | {entry['wall_s']:.2f} |"
        )
    lines.extend(
        [
            "",
            "Zero-failure points are never reported as zero: the table substitutes the 95% Clopper–Pearson upper confidence bound and labels it as an upper bound. Block-to-round conversion exactly copies EXP-010: `-expm1(log1p(-block_ler) / 12)`. The k=12 surface block first applies `1-(1-single_patch_block_ler)^12`; interval endpoints are transformed monotonically in both steps.",
            "",
            "Gross and PBB have no p=0.001 or p=0.003 measurements in the permitted source artifacts, so those cells are intentionally unavailable rather than fabricated. Their complete imported comparison point is p=0.002; the surface sweep covers all requested p values in both bases.",
            "",
            "## Circuit-derived surface resources",
            "",
        ]
    )
    per_patch = result["surface"]["resources"]["per_patch"]
    audit = result["surface"]["resources"]["noise_instruction_audit_per_round"]
    lines.extend(
        [
            f"The generated circuit contains {per_patch['data_qubits']} final-data-measurement targets and {per_patch['ancilla_qubits']} other active QUBIT_COORDS targets, hence {per_patch['total_active_qubits']} active physical qubits per patch. Stim's `num_qubits={per_patch['stim_qubit_index_span']}` is only the sparse index span (maximum index plus one), not the active physical-qubit count.",
            f"Exact fixed-circuit counts per round are {per_patch['depth_ticks_per_round']} TICKs, {per_patch['two_qubit_gates_per_round']} CX pairs, {audit['data_depolarize1_target_locations']} data DEPOLARIZE1 target locations, {audit['ancilla_depolarize1_target_locations']} ancilla DEPOLARIZE1 target locations, and {audit['h_target_locations']} H target locations. The equality of ancilla DEPOLARIZE1 and H target counts reflects post-H noise, not ancilla-idle noise.",
            "",
            "## Convention deltas — mapped baseline, not identical noise",
            "",
            "The same numeric p therefore does **not** denote the same set or number of noisy circuit locations.",
            "",
        ]
    )
    for i, delta in enumerate(result["convention_deltas"], start=1):
        lines.append(
            f"{i}. **{delta['topic']}.** Stim generated: {delta['stim_generated']}. Custom qLDPC: {delta['custom_qldpc']}."
        )
    lines.extend(
        [
            "",
            "## Imported evidence and provenance",
            "",
            "- `results/raw/exp013_isolated_latency.json`: Gross/PBB depth and two-qubit gates per round. No latency number is used in this baseline.",
            "- `results/raw/exp016_schedule_controlled.json`: persisted slot maps (used to count 144 data indices and 144 check/ancilla indices), p=0.002 shots/failures, schedule-mean block LER, and schedule-mean t intervals. The five schedules per code were a fixed-seed uniform sample from completely enumerated minimum-depth translation-invariant schedules.",
            "- Neither imported artifact was rerun or modified.",
            "- EXP-030 performed no code-distance search. Gross distance 12 has an independent exact certificate. PBB distance 12 is catalogue metadata pending an independent full-symplectic certificate. Surface `d=12` is Stim's nominal generated-code parameter; no independent physical-location circuit-distance certification was performed.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "cd /Users/jinleic/jinleic-workspace/math/qec",
            "PYTHONPATH=src .venv/bin/python experiments/exp030_surface_baseline.py",
            "```",
            "",
            f"Pre-registered point order is p outermost and Z then X within p; seeds are 20260812+i. Total new shots: {sum(p['shots'] for p in result['surface']['points'])}. Experiment wall time: {result['wall_s']:.2f} s.",
            f"Environment: Python {result['environment']['python']}, Stim {result['environment']['stim']}, PyMatching {result['environment']['pymatching']}, NumPy {result['environment']['numpy']}, SciPy {result['environment']['scipy']}; {result['environment']['platform']}.",
            f"Shared-machine caveat: one process ran batched PyMatching decoding while other agents could be active; load average moved from {result['environment']['load_at_start']} to {result['environment']['load_at_end']}. Wall times are operational metadata, not latency benchmarks.",
            "",
            "## Caveats",
            "",
            "- The surface k=12 rate and resource row assumes twelve independent, identical patches operated in parallel. It is an analytic monotone mapping of one-patch samples, not a jointly sampled 12-patch circuit.",
            "- Z-memory is the direct basis match to EXP-016's qLDPC measurement; X-memory is reported separately and is not averaged into it.",
            "- Gross/PBB confidence intervals are schedule-mean t intervals that include measured schedule-to-schedule variation. Surface intervals are binomial Clopper–Pearson intervals conditional on Stim's fixed generated schedule. They answer related but not identical uncertainty questions.",
            "- PyMatching and BP+OSD are different decoders, and the two DEM treatments differ. Decoder quality is part of this mapped end-to-end baseline.",
            "- No claim of matched target logical error, identical noise-location count, identical geometry/connectivity, or globally optimized wall-clock cycle duration is made.",
        ]
    )
    return "\n".join(lines) + "\n"


def package_version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def main() -> None:
    started_clock = time.perf_counter()
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    load_at_start = os.getloadavg()

    surface_points = []
    shared_resources: dict[str, Any] | None = None
    for index, (p, basis, task) in enumerate(POINT_ORDER):
        point = simulate_point(index, p, basis, task)
        if shared_resources is None:
            shared_resources = point["resources"]
        elif point["resources"] != shared_resources:
            raise AssertionError("surface resource counts changed across p or memory basis")
        surface_points.append(point)
        print(
            f"p={p:.3f} basis={basis} seed={point['seed']} shots={point['shots']} "
            f"fails={point['logical_failures']} mapped-k12-round="
            f"{point['matched_k12_ler_per_round']['reported_ler']:.8g} "
            f"wall={point['wall_s']:.2f}s",
            flush=True,
        )
    assert shared_resources is not None

    imported = imported_qldpc_results()
    completed_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    wall_s = time.perf_counter() - started_clock
    resource_rows = [
        {
            "code": entry["code"],
            **entry["resources"],
            "scaling": "native k=12 block",
        }
        for entry in imported
    ]
    resource_rows.append(
        {
            "code": "12 x rotated surface code d=12",
            "data_qubits": shared_resources["matched_k12_parallel_patches"]["data_qubits"],
            "ancilla_qubits": shared_resources["matched_k12_parallel_patches"]["ancilla_qubits"],
            "total_qubits": shared_resources["matched_k12_parallel_patches"]["total_active_qubits"],
            "depth_per_round": shared_resources["matched_k12_parallel_patches"]["depth_ticks_per_round"],
            "depth_status": shared_resources["per_patch"]["depth_status"],
            "two_qubit_gates_per_round": shared_resources["matched_k12_parallel_patches"]["two_qubit_gates_per_round"],
            "scaling": "12 independent patches in parallel",
        }
    )

    surface_by_key = {(point["p"], point["basis"]): point for point in surface_points}
    rate_matrix = []
    for p in PHYSICAL_ERROR_RATES:
        row: dict[str, Any] = {"p": p}
        for entry in imported:
            row[entry["code"]] = (
                {
                    "block_ler": entry["block_ler"],
                    "ler_per_round": entry["ler_per_round"],
                }
                if p == 0.002
                else None
            )
        for basis, _ in MEMORY_TASKS:
            point = surface_by_key[(p, basis)]
            row[f"surface_memory_{basis.lower()}"] = {
                "shots": point["shots"],
                "failures": point["logical_failures"],
                "block_ler": point["matched_k12_12_round_block"],
                "ler_per_round": point["matched_k12_ler_per_round"],
            }
        rate_matrix.append(row)

    gross_total = imported[0]["resources"]["total_qubits"]
    surface_total = shared_resources["matched_k12_parallel_patches"]["total_active_qubits"]
    result = {
        "experiment": EXPERIMENT,
        "description": (
            "Real Stim/PyMatching rotated-surface-code d=12 sweep and imported "
            "Gross/PBB p=0.002 Pareto comparison"
        ),
        "started_utc": started_utc,
        "completed_utc": completed_utc,
        "protocol": {
            "pre_registered": True,
            "distance": DISTANCE,
            "distance_status": (
                "nominal generated surface-code distance parameter; no independent "
                "physical-location circuit-distance certification"
            ),
            "rounds": ROUNDS,
            "physical_error_rates": list(PHYSICAL_ERROR_RATES),
            "memory_tasks": [
                {"basis": basis, "stim_task": task} for basis, task in MEMORY_TASKS
            ],
            "shots_per_p_basis_point": SHOTS_PER_POINT,
            "total_new_shots": SHOTS_PER_POINT * len(POINT_ORDER),
            "batch_shots": BATCH_SHOTS,
            "base_seed": BASE_SEED,
            "seed_rule": "20260812+i in p-outer, Z-then-X point order",
            "point_seeds": [
                {
                    "i": i,
                    "p": p,
                    "basis": basis,
                    "stim_task": task,
                    "seed": POINT_SEEDS[i],
                }
                for i, (p, basis, task) in enumerate(POINT_ORDER)
            ],
            "confidence_interval": "two-sided equal-tailed Clopper-Pearson, 95%",
            "zero_failure_reporting": (
                "use the Clopper-Pearson 95% upper bound, never zero"
            ),
            "generator_noise_arguments": GENERATOR_NOISE_ARGUMENTS,
            "decoder": (
                "pymatching.Matching.from_detector_error_model after "
                "stim detector_error_model(decompose_errors=True)"
            ),
            "logical_failure": "decoded logical observable differs from sampled observable",
            "matched_k_mapping": "12 statistically independent surface patches",
            "block_to_round_formula": "-expm1(log1p(-block_ler)/12)",
        },
        "sources": {
            "qldpc_resources": "results/raw/exp013_isolated_latency.json",
            "qldpc_statistics_and_slot_maps": "results/raw/exp016_schedule_controlled.json",
            "sources_rerun": False,
        },
        "surface": {
            "resources": shared_resources,
            "points": surface_points,
        },
        "qldpc_imported": imported,
        "pareto_resources": resource_rows,
        "rate_matrix": rate_matrix,
        "convention_deltas": list(CONVENTION_DELTAS),
        "verdict": {
            "classification": VERDICT,
            "table_complete_at_p": 0.002,
            "both_surface_memory_bases_complete_at_all_requested_p": True,
            "gross_pbb_unavailable_p": [0.001, 0.003],
            "surface_vs_gross_total_qubits": {
                "surface_matched_k12": surface_total,
                "gross": gross_total,
                "difference": surface_total - gross_total,
                "ratio": surface_total / gross_total,
                "surface_uses": "more",
            },
            "comparability": "mapped baseline; generated and custom noise placements differ",
        },
        "environment": {
            "python": platform.python_version(),
            "stim": stim.__version__,
            "pymatching": package_version("pymatching"),
            "numpy": np.__version__,
            "scipy": package_version("scipy"),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "workers": 1,
            "load_at_start": list(load_at_start),
            "load_at_end": list(os.getloadavg()),
            "shared_load_caveat": (
                "other agents may have been active; wall_s is not a latency benchmark"
            ),
        },
        "wall_s": wall_s,
        "machine_readable_verdict": VERDICT,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    REPORT.write_text(make_report(result))
    print(f"wrote {OUTPUT.relative_to(ROOT)}", flush=True)
    print(f"wrote {REPORT.relative_to(ROOT)}", flush=True)
    print(f"verdict={VERDICT} wall={wall_s:.2f}s", flush=True)


if __name__ == "__main__":
    main()
