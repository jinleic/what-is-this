"""EXP-010: Stim circuit-level BB memory baseline and validation.

This experiment builds the depth-7 syndrome-extraction circuit, validates its
Stim detector model and noiseless behavior, searches for a circuit-distance
upper bound, then decodes circuit samples with BP-OSD over the detector error
model's fault variables.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import multiprocessing as mp
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

# Each process runs one C++ BP decoder.  Prevent hidden BLAS/OpenMP fan-out.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import ldpc  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse  # noqa: E402
from scipy.stats import beta  # noqa: E402
import stim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ldpc.ckt_noise.dem_matrices import detector_error_model_to_check_matrices  # noqa: E402
from qec_research.circuits.bb_syndrome import (  # noqa: E402
    DEPTH7_SCHEDULE,
    build_bb_memory_circuit,
)
from qec_research.codes.bicycle import BRAVYI_BB, build_bb  # noqa: E402

OUTPUT = ROOT / "results" / "raw" / "exp010_bb_circuit_baseline.json"
PHYSICAL_ERROR_RATES = (0.001, 0.002, 0.003, 0.005)
CODE_ROUNDS = {"[[72,12,6]]": 6, "[[144,12,12]]": 12}
BASE_SEED = 10_082_023_008_079_015
NOISE_MODEL_P_FOR_VALIDATION = 0.001
NOISELESS_VALIDATION_SEED = 230_807_915
NOISELESS_VALIDATION_SHOTS = 10_000

DECODER_PARAMETERS = {
    "implementation": "ldpc.BpOsdDecoder",
    "bp_method": "ms",
    "max_iter": 200,
    "ms_scaling_factor": 0.0,
    "schedule": "parallel",
    "omp_thread_count": 1,
    "osd_method": "osd_0",
    "osd_order": 0,
}

ARTIFACT_EVIDENCE = {
    "schedule": [
        "third_party/BivariateBicycleCodes/decoder_setup.py:45-51",
        "third_party/BivariateBicycleCodes/decoder_setup.py:208-264",
    ],
    "direction_mapping": "third_party/BivariateBicycleCodes/decoder_setup.py:181-205",
    "noise_strengths": "third_party/BivariateBicycleCodes/decoder_setup.py:35-40",
    "noise_sampling": "third_party/BivariateBicycleCodes/decoder_run.py:84-190",
    "decoder_parameters": "third_party/BivariateBicycleCodes/decoder_run.py:67-72",
}

# Globals are populated before forking so large sparse matrices are inherited
# copy-on-write instead of serialized once per worker.
_WORKER_CIRCUIT_TEXT: str | None = None
_WORKER_CHECK_MATRIX: scipy.sparse.csc_matrix | None = None
_WORKER_OBSERVABLES_MATRIX: scipy.sparse.csc_matrix | None = None
_WORKER_PRIORS: np.ndarray | None = None


def clopper_pearson(failures: int, shots: int, confidence: float = 0.95) -> tuple[float, float]:
    """Exact equal-tailed binomial confidence interval."""
    if shots <= 0 or failures < 0 or failures > shots:
        raise ValueError("require 0 <= failures <= shots and shots > 0")
    alpha = 1.0 - confidence
    lower = 0.0 if failures == 0 else float(beta.ppf(alpha / 2, failures, shots - failures + 1))
    upper = 1.0 if failures == shots else float(
        beta.ppf(1 - alpha / 2, failures + 1, shots - failures)
    )
    return lower, upper


def per_round_rate(per_shot_rate: float, rounds: int) -> float:
    return float(-np.expm1(np.log1p(-per_shot_rate) / rounds)) if per_shot_rate < 1 else 1.0


def _package_versions() -> dict[str, str]:
    distributions = (
        "stim",
        "sinter",
        "pymatching",
        "ldpc",
        "qldpc",
        "numpy",
        "scipy",
        "networkx",
        "python-sat",
        "ortools",
    )
    versions: dict[str, str] = {}
    for distribution in distributions:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "not installed"
    return versions


def _circuit_for(code: str, p: float) -> stim.Circuit:
    hx, hz = build_bb(BRAVYI_BB[code])
    return build_bb_memory_circuit(hx, hz, rounds=CODE_ROUNDS[code], p=p, basis="Z")


def validate_circuits(run_distance_search: bool = True) -> dict[str, Any]:
    """Execute acceptance gates (a), (b), (c), and (d)."""
    per_code: dict[str, Any] = {}
    all_passed = True
    for code, rounds in CODE_ROUNDS.items():
        hx, _ = build_bb(BRAVYI_BB[code])
        n = hx.shape[1]

        noisy = _circuit_for(code, NOISE_MODEL_P_FOR_VALIDATION)
        dem = noisy.detector_error_model(decompose_errors=False)
        gate_a = {
            "passed": True,
            "p": NOISE_MODEL_P_FOR_VALIDATION,
            "detectors": noisy.num_detectors,
            "observables": noisy.num_observables,
            "dem_error_instructions": dem.num_errors,
        }

        noiseless = _circuit_for(code, 0.0)
        detector_samples, observable_samples = noiseless.compile_detector_sampler(
            seed=NOISELESS_VALIDATION_SEED
        ).sample(NOISELESS_VALIDATION_SHOTS, separate_observables=True)
        detection_event_count = int(np.count_nonzero(detector_samples))
        logical_flip_count = int(np.count_nonzero(observable_samples))
        detection_event_shots = int(np.count_nonzero(np.any(detector_samples, axis=1)))
        logical_failure_shots = int(np.count_nonzero(np.any(observable_samples, axis=1)))
        gate_b_passed = detection_event_count == 0 and logical_flip_count == 0
        gate_b = {
            "passed": gate_b_passed,
            "p": 0.0,
            "shots": NOISELESS_VALIDATION_SHOTS,
            "seed": NOISELESS_VALIDATION_SEED,
            "detection_events": detection_event_count,
            "shots_with_detection_events": detection_event_shots,
            "logical_observable_flips": logical_flip_count,
            "shots_with_logical_failure": logical_failure_shots,
        }

        gate_d_passed = noisy.num_qubits == 2 * n
        gate_d = {
            "passed": gate_d_passed,
            "code_data_qubits": n,
            "stim_num_qubits": noisy.num_qubits,
            "expected_num_qubits": 2 * n,
            "x_ancillas": n // 2,
            "z_ancillas": n // 2,
        }
        all_passed &= gate_b_passed and gate_d_passed
        per_code[code] = {"rounds": rounds, "gate_a_dem_build": gate_a,
                          "gate_b_noiseless": gate_b, "gate_d_qubit_count": gate_d}
        print(
            f"{code}: DEM built ({dem.num_errors} errors); noiseless "
            f"events={detection_event_count}, flips={logical_flip_count}; "
            f"qubits={noisy.num_qubits}=2*{n}",
            flush=True,
        )

    gate_c: dict[str, Any] = {"code": "[[72,12,6]]", "executed": run_distance_search}
    if run_distance_search:
        circuit = _circuit_for("[[72,12,6]]", NOISE_MODEL_P_FOR_VALIDATION)
        try:
            graphlike = circuit.shortest_graphlike_error()
            gate_c.update({
                "passed": True,
                "method": "stim.Circuit.shortest_graphlike_error",
                "fault_set_size": len(graphlike),
                "is_exact": True,
                "bound": "exact within the graphlike subset",
                "fault_set": [str(error) for error in graphlike],
            })
        except ValueError as error:
            search_parameters = {
                "dont_explore_detection_event_sets_with_size_above": 6,
                "dont_explore_edges_with_degree_above": 6,
                "dont_explore_edges_increasing_symptom_degree": True,
                "canonicalize_circuit_errors": True,
            }
            fault_set = circuit.search_for_undetectable_logical_errors(**search_parameters)
            gate_c.update({
                "passed": len(fault_set) > 0,
                "method": "stim.Circuit.search_for_undetectable_logical_errors",
                "graphlike_search_result": str(error),
                "search_parameters": search_parameters,
                "fault_set_size": len(fault_set),
                "is_exact": False,
                "bound": "heuristic upper bound on circuit distance",
                "fault_set": [str(item) for item in fault_set],
            })
        all_passed &= bool(gate_c["passed"])
        print(
            f"[[72,12,6]] circuit-distance search: {gate_c['fault_set_size']} faults "
            f"({gate_c['bound']})",
            flush=True,
        )
    gate_c["passed"] = bool(gate_c.get("passed", False)) if run_distance_search else None
    return {
        "all_executed_gates_passed": all_passed if run_distance_search else None,
        "per_code": per_code,
        "gate_c_circuit_distance": gate_c,
    }


def _decode_worker(task: tuple[int, int]) -> dict[str, int | float]:
    shots, seed = task
    if (_WORKER_CIRCUIT_TEXT is None or _WORKER_CHECK_MATRIX is None
            or _WORKER_OBSERVABLES_MATRIX is None or _WORKER_PRIORS is None):
        raise RuntimeError("decoder worker globals were not initialized")
    circuit = stim.Circuit(_WORKER_CIRCUIT_TEXT)
    detector_samples, observable_samples = circuit.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True
    )
    decoder = ldpc.BpOsdDecoder(
        _WORKER_CHECK_MATRIX,
        error_channel=_WORKER_PRIORS.tolist(),
        max_iter=DECODER_PARAMETERS["max_iter"],
        bp_method=DECODER_PARAMETERS["bp_method"],
        ms_scaling_factor=DECODER_PARAMETERS["ms_scaling_factor"],
        schedule=DECODER_PARAMETERS["schedule"],
        omp_thread_count=DECODER_PARAMETERS["omp_thread_count"],
        osd_method=DECODER_PARAMETERS["osd_method"],
        osd_order=DECODER_PARAMETERS["osd_order"],
    )
    failures = 0
    bp_converged = 0
    detector_event_shots = 0
    started = time.perf_counter()
    for detector_sample, actual_observables in zip(
        detector_samples, observable_samples, strict=True
    ):
        detector_event_shots += bool(np.any(detector_sample))
        correction = decoder.decode(detector_sample)
        bp_converged += bool(decoder.converge)
        predicted_observables = np.asarray(
            _WORKER_OBSERVABLES_MATRIX @ correction, dtype=np.uint8
        ).reshape(-1) & 1
        failures += bool(np.any(predicted_observables != actual_observables))
    return {
        "shots": shots,
        "failures": failures,
        "bp_converged_shots": bp_converged,
        "shots_with_detection_events": detector_event_shots,
        "seed": seed,
        "wall_time_s": time.perf_counter() - started,
    }


def _split_work(shots: int, workers: int, point_seed: int) -> list[tuple[int, int]]:
    workers = min(workers, shots)
    quotient, remainder = divmod(shots, workers)
    return [
        (quotient + (worker < remainder), point_seed + worker)
        for worker in range(workers)
    ]


def benchmark_point(code: str, p: float, shots: int, workers: int, point_seed: int) -> dict[str, Any]:
    """Sample and BP-OSD decode one (code, p) point."""
    global _WORKER_CIRCUIT_TEXT, _WORKER_CHECK_MATRIX
    global _WORKER_OBSERVABLES_MATRIX, _WORKER_PRIORS

    circuit = _circuit_for(code, p)
    dem = circuit.detector_error_model(decompose_errors=False)
    matrices = detector_error_model_to_check_matrices(
        dem, allow_undecomposed_hyperedges=True
    )
    _WORKER_CIRCUIT_TEXT = str(circuit)
    _WORKER_CHECK_MATRIX = matrices.check_matrix.tocsc()
    _WORKER_OBSERVABLES_MATRIX = matrices.observables_matrix.tocsc()
    _WORKER_PRIORS = np.asarray(matrices.priors, dtype=np.float64)

    tasks = _split_work(shots, workers, point_seed)
    started = time.perf_counter()
    if len(tasks) == 1:
        chunks = [_decode_worker(tasks[0])]
        process_start_method = "single-process"
    else:
        if "fork" not in mp.get_all_start_methods():
            raise RuntimeError("parallel execution requires the fork multiprocessing start method")
        process_start_method = "fork"
        with mp.get_context("fork").Pool(processes=len(tasks)) as pool:
            chunks = pool.map(_decode_worker, tasks)
    wall_time = time.perf_counter() - started

    completed_shots = sum(int(chunk["shots"]) for chunk in chunks)
    failures = sum(int(chunk["failures"]) for chunk in chunks)
    if completed_shots != shots:
        raise AssertionError(f"requested {shots} shots but workers returned {completed_shots}")
    ler = failures / shots
    ci_lower, ci_upper = clopper_pearson(failures, shots)
    rounds = CODE_ROUNDS[code]
    if failures == 0:
        reported_ler = ci_upper
        rate_report_type = "clopper_pearson_95_upper_bound"
    else:
        reported_ler = ler
        rate_report_type = "empirical_point_estimate"
    record = {
        "code": code,
        "n": BRAVYI_BB[code].n,
        "k": 12,
        "distance": rounds,
        "basis": "Z",
        "p": p,
        "rounds": rounds,
        "shots": shots,
        "logical_failures": failures,
        "logical_error_rate_per_shot": reported_ler,
        "logical_error_rate_per_round": per_round_rate(reported_ler, rounds),
        "logical_error_rate_report_type": rate_report_type,
        "clopper_pearson_95_ci_per_shot": {"lower": ci_lower, "upper": ci_upper},
        "clopper_pearson_95_ci_per_round": {
            "lower": per_round_rate(ci_lower, rounds),
            "upper": per_round_rate(ci_upper, rounds),
        },
        "bp_converged_shots": sum(int(chunk["bp_converged_shots"]) for chunk in chunks),
        "shots_with_detection_events": sum(
            int(chunk["shots_with_detection_events"]) for chunk in chunks
        ),
        "detector_error_model": {
            "decompose_errors": False,
            "detectors": dem.num_detectors,
            "observables": dem.num_observables,
            "error_instructions": dem.num_errors,
            "check_matrix_shape": list(_WORKER_CHECK_MATRIX.shape),
            "check_matrix_nonzeros": int(_WORKER_CHECK_MATRIX.nnz),
        },
        "decoder": DECODER_PARAMETERS.copy(),
        "worker_count": len(tasks),
        "worker_omp_threads": DECODER_PARAMETERS["omp_thread_count"],
        "process_start_method": process_start_method,
        "stim_seeds": [int(chunk["seed"]) for chunk in chunks],
        "worker_shots": [int(chunk["shots"]) for chunk in chunks],
        "wall_time_s": wall_time,
    }
    rate_label = "LER_CP95_upper" if failures == 0 else "LER"
    print(
        f"{code} p={p:.4g} rounds={rounds} shots={shots} failures={failures} "
        f"{rate_label}={reported_ler:.8g} CP95=[{ci_lower:.8g}, {ci_upper:.8g}] "
        f"time={wall_time:.1f}s",
        flush=True,
    )
    return record


def run(
    *,
    codes: list[str],
    probabilities: list[float],
    shots: int,
    workers: int,
    output: Path,
    run_validation: bool,
    run_distance_search: bool,
) -> dict[str, Any]:
    if shots < 1:
        raise ValueError("shots must be positive")
    if workers < 1:
        raise ValueError("workers must be positive")
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    validation = validate_circuits(run_distance_search) if run_validation else {
        "executed": False
    }
    benchmarks: list[dict[str, Any]] = []
    for code_index, code in enumerate(codes):
        for probability_index, p in enumerate(probabilities):
            point_seed = BASE_SEED + code_index * 1_000_000 + probability_index * 10_000
            benchmarks.append(benchmark_point(code, p, shots, workers, point_seed))

    result = {
        "experiment": "EXP-010",
        "description": "Stim circuit-level Z-memory baseline for CSS bivariate-bicycle codes",
        "started_utc": started_utc,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "artifact_commit": "fa77e3333d3ec44c79d8f914dd24c040d1da471b",
        "artifact_evidence": ARTIFACT_EVIDENCE,
        "recovered_schedule": {
            "depth_7_layers_xz": [list(layer) for layer in DEPTH7_SCHEDULE],
            "sX": [layer[0] for layer in DEPTH7_SCHEDULE],
            "sZ": [layer[1] for layer in DEPTH7_SCHEDULE],
            "timeline_depth_including_final_measure_reset_slice": 8,
        },
        "noise_model": {
            "description": (
                "Uniform strength p: DEPOLARIZE2 after CNOT, DEPOLARIZE1 on each "
                "scheduled data idle, basis-opposite Pauli error after reset and before measurement"
            ),
            "two_qubit": "DEPOLARIZE2(p) after CX",
            "idle": "DEPOLARIZE1(p)",
            "z_reset": "R then X_ERROR(p)",
            "x_reset": "RX then Z_ERROR(p)",
            "z_measurement": "X_ERROR(p) then MZ",
            "x_measurement": "Z_ERROR(p) then MX",
        },
        "validation": validation,
        "benchmarks": benchmarks,
        "confidence_interval": "two-sided equal-tailed Clopper-Pearson, 95%",
        "logical_failure_definition": (
            "at least one of the 12 decoded pure-Z logical observables differs from Stim's sampled frame"
        ),
        "software_versions": _package_versions(),
        "hardware": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "os_cpu_count": os.cpu_count(),
            "decoder_processes_requested": workers,
            "decoder_omp_threads_per_process": 1,
            "thread_environment": {
                name: os.environ[name]
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                )
            },
        },
        "seeding": {
            "base_seed": BASE_SEED,
            "noiseless_validation_seed": NOISELESS_VALIDATION_SEED,
            "point_rule": "base + code_index*1000000 + probability_index*10000 + worker_index",
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {output}", flush=True)
    return result


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=20_000, help="shots per code/p point")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--code", action="append", choices=tuple(CODE_ROUNDS), dest="codes")
    parser.add_argument("--p", action="append", type=float, dest="probabilities")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--skip-validation", action="store_true")
    parser.add_argument("--skip-distance-search", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = _parse_args()
    run(
        codes=arguments.codes or list(CODE_ROUNDS),
        probabilities=arguments.probabilities or list(PHYSICAL_ERROR_RATES),
        shots=arguments.shots,
        workers=arguments.workers,
        output=arguments.output,
        run_validation=not arguments.skip_validation,
        run_distance_search=not arguments.skip_distance_search,
    )
