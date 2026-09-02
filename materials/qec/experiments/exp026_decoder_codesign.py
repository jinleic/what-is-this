"""EXP-026 v2: finite decoder co-design with paired samples and hard deadlines.

Protocol v1 (29 configurations x 2 codes x 20,000 shots) was retired before
production data after 0/58 rows completed and a conservative >216
baseline-equivalent core-hour lower bound.  Its smoke/guard artifacts remain in
results/partial_runs as feasibility evidence.  V2 freezes five candidates,
paired held-out samples, supervised 30 s decode deadlines, and isolated timing.
"""
from __future__ import annotations

import os
for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import platform
import queue
import statistics
import tempfile
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.sparse as sp
import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits, pure_z_logical_basis,
)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec, build_memory_circuit, generator_supports,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule, slots_to_layers, verify_schedule,
)
from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec, BRAVYI_BB, PBBSpec, build_bb,
)
from qec_research.decoders.bposd_dem import (  # noqa: E402
    clopper_pearson, dem_to_matrices,
)

# --------------------------- frozen protocol v2 ---------------------------
P = 0.002
ROUNDS = 12
PBB_CODE_ID = "12_6_0193"
DEADLINE_S = 30.0
CHECKPOINT_EVERY = 25
DEV_SHOTS = 200
DEV_SEED = 20_260_813
HELDOUT_SHOTS = 5_000
HELDOUT_SEED_BASE = 20_270_813
TIMING_WARMUP = 50
TIMING_SHOTS = 500
TIMING_SEED_BASE = 20_280_813
TIMING_BLOCK = 25
BOOTSTRAP_REPS = 10_000
BOOTSTRAP_SEED = 20_290_813
NOISELESS_SHOTS = 4_000
MAX_IDLE_LOAD = 4.0
CANDIDATES = (
    {"id": "baseline_bposd", "kind": "bposd", "max_iter": 30,
     "alpha": 0.625, "osd_method": "osd0", "osd_order": 0},
    {"id": "shadow_two_stage", "kind": "shadow", "max_iter": 30,
     "alpha": 0.625, "osd_method": "osd0", "osd_order": 0},
    {"id": "bplsd_b4", "kind": "bplsd", "max_iter": 30,
     "alpha": 0.625, "bits_per_step": 4, "lsd_method": "lsd_0", "lsd_order": 0},
    {"id": "bplsd_b16", "kind": "bplsd", "max_iter": 30,
     "alpha": 0.625, "bits_per_step": 16, "lsd_method": "lsd_0", "lsd_order": 0},
    {"id": "bplsd_b64", "kind": "bplsd", "max_iter": 30,
     "alpha": 0.625, "bits_per_step": 64, "lsd_method": "lsd_0", "lsd_order": 0},
)
CODE_KEYS = ("gross", "pbb")
EXPECTED_DEPTH = {"gross": 7, "pbb": 8}

CATALOGUE = ROOT / "third_party/qcode-discovery/results/campaign7_publication_merged.jsonl"
PARTIAL = ROOT / "results/partial_runs/exp026_decoder_codesign-v2-partial.json"
QUARANTINE = ROOT / "results/quarantine/exp026_decoder_codesign-v2-failing.json"
RAW = ROOT / "results/raw/exp026_decoder_codesign.json"
PROCESSED = ROOT / "results/processed/exp026_decoder_codesign.json"
REPORT = ROOT / "notes/agent_reports/exp026_decoder_codesign.md"
V1_GUARD = ROOT / "results/partial_runs/exp026_decoder_codesign-artifact-guard.json"
V1_SMOKE = ROOT / "results/partial_runs/exp026_decoder_codesign-smoke.json"


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise


def write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True))


def protocol() -> dict[str, Any]:
    return {
        "version": 2,
        "revision_reason": {
            "v1_production_rows": 0,
            "v1_declared_rows": 58,
            "baseline_equivalent_core_hour_lower_bound": 216,
            "revision_precedes_production_data": True,
            "v1_smoke_artifact": str(V1_SMOKE.relative_to(ROOT)),
            "v1_guard_artifact": str(V1_GUARD.relative_to(ROOT)),
        },
        "p": P, "rounds": ROUNDS, "deadline_s": DEADLINE_S,
        "checkpoint_every_shots_max": CHECKPOINT_EVERY,
        "dev": {
            "code": "pbb", "shots": DEV_SHOTS, "seed": DEV_SEED,
            "paired": True, "claims_allowed": False,
            "prior_development_note": (
                "Earlier first-stream 10-shot smoke observed b1/b4/b16; b1 was "
                "omitted from v2 for slowness. This is development data only."
            ),
            "all_candidates_proceed_unless": "fatal API/correction failure",
        },
        "heldout": {
            "shots_per_code": HELDOUT_SHOTS,
            "seed_formula": "20270813 + code_index", "paired": True,
            "operational_failure": (
                "logical mismatch OR invalid syndrome correction OR decoder timeout"
            ),
            "early_stopping": False,
            "isolation_gate": "same load/process gate as timing; passed per session",
            "resume_provenance": (
                "each resumed accuracy session requires a newly passed gate snapshot"
            ),
        },
        "timing": {
            "warmup_discarded": TIMING_WARMUP, "shots_per_code": TIMING_SHOTS,
            "seed_formula": "20280813 + code_index", "paired": True,
            "candidate_order": "seeded balanced rotations in 25-shot blocks",
            "deadline_s": DEADLINE_S, "bootstrap_reps": BOOTSTRAP_REPS,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "resumability": (
                "intentionally non-resumable: any interrupted timing pass is "
                "discarded and the entire isolated two-code pass reruns"
            ),
            "isolation_gate": {
                "load_max": "absolute 1-minute load average <= 4.0",
                "continuous_recheck": "startup and at least every <=25-shot block",
                "heavy_processes_forbidden": (
                    "any non-self Python qec experiment process or known solver job"
                ),
            },
        },
        "candidates": list(CANDIDATES),
        "verdict": (
            "PBB challenger joint improvement iff lower operational-failure point "
            "estimate, one-sided paired McNemar superiority with Holm correction "
            "across four challengers at alpha .05, and bootstrap upper 95% bound "
            "on paired p95 latency ratio <1; otherwise NEGATIVE within finite v2 set"
        ),
        "thread_limits": {name: 1 for name in (
            "OMP", "OPENBLAS", "MKL", "VECLIB", "NUMEXPR"
        )},
    }


def protocol_hash() -> str:
    return hashlib.sha256(json.dumps(protocol(), sort_keys=True).encode()).hexdigest()


def load_member() -> dict[str, Any]:
    rows = [json.loads(line) for line in CATALOGUE.open()]
    matches = [row for row in rows if row.get("code_id") == PBB_CODE_ID]
    if len(matches) != 1:
        raise RuntimeError(f"expected one {PBB_CODE_ID}, found {len(matches)}")
    return matches[0]


def _pbb_spec(row: dict[str, Any]) -> PBBSpec:
    return PBBSpec(
        int(row["ell"]), int(row["m"]),
        [tuple(x) for x in row["A_terms"]], [tuple(x) for x in row["B_terms"]],
        [tuple(x) for x in row["C_terms"]], [tuple(x) for x in row["D_terms"]],
        name=PBB_CODE_ID,
    )


def _schedule(code, supports, orbits, depth: int):
    result = cpsat_schedule(
        supports, code.n, T=depth, time_limit_s=180, workers=1,
        symmetry_orbits=orbits, random_seed=0, deterministic=True,
    )
    if result.slot is None or not result.verification["valid"]:
        raise RuntimeError(f"no valid fixed depth-{depth} schedule: {result.status}")
    if not verify_schedule(supports, result.slot)["valid"]:
        raise RuntimeError("independent schedule verification failed")
    return result


def _pure_rows(circuit: stim.Circuit, code) -> np.ndarray:
    pure = {s.index for s in generator_supports(code.H) if s.is_pure_z}
    coords = circuit.get_detector_coordinates()
    rows = [d for d in range(circuit.num_detectors) if int(coords[d][0]) in pure]
    return np.asarray(rows, dtype=np.int64)


def build_circuits() -> dict[str, dict[str, Any]]:
    gross_spec = BRAVYI_BB["[[144,12,12]]"]
    gross = bb_supports_and_orbits(gross_spec)
    row = load_member()
    pspec = _pbb_spec(row)
    pbb = pbb_supports_and_orbits(pspec)
    entries = {
        "gross": ("CSS-BB [[144,12,12]] Gross", gross_spec, *gross),
        "pbb": (f"nonCSS-PBB [[144,12,12]] {PBB_CODE_ID}", pspec, *pbb),
    }
    out = {}
    for ci, (key, (label, spec, code, supports, orbits)) in enumerate(entries.items()):
        sch = _schedule(code, supports, orbits, EXPECTED_DEPTH[key])
        layers = slots_to_layers(sch.slot, EXPECTED_DEPTH[key])
        observables = pure_z_logical_basis(code)
        c0, _ = build_memory_circuit(CircuitSpec(
            H=code.H, observables=observables, rounds=ROUNDS,
            p=0, basis="Z", layers=layers,
        ))
        d0, o0 = c0.compile_detector_sampler(seed=DEV_SEED + 50_000 + ci).sample(
            NOISELESS_SHOTS, separate_observables=True,
        )
        if d0.any() or o0.any():
            raise RuntimeError(f"{key} noiseless validation failed")
        circuit, meta = build_memory_circuit(CircuitSpec(
            H=code.H, observables=observables, rounds=ROUNDS,
            p=P, basis="Z", layers=layers,
        ))
        dem = circuit.detector_error_model(decompose_errors=False)
        slot_payload = [[int(a), int(b), int(t)] for (a, b), t in sorted(sch.slot.items())]
        rec = {
            "label": label, "depth": EXPECTED_DEPTH[key],
            "schedule_status": sch.status,
            "schedule_hash": hashlib.sha256(json.dumps(slot_payload).encode()).hexdigest()[:16],
            "slot_map": slot_payload,
            "schedule_class": (
                "translation-invariant witness; exact only within TI class; unrestricted "
                "depth-6 exclusion is external ASC" if key == "gross" else
                "class-free exact from max-check-weight lower bound plus depth-8 witness"
            ),
            "schedule_provenance": (
                "EXP-013 deterministic CP-SAT witness, symmetry orbits, seed 0"
            ),
            "dem_errors": dem.num_errors, "dem_detectors": dem.num_detectors,
            "circuit_sha256": hashlib.sha256(str(circuit).encode()).hexdigest(),
            "noiseless_shots": NOISELESS_SHOTS,
            "noiseless_detector_firings": int(d0.sum()),
            "noiseless_observable_flips": int(o0.sum()),
            "meta": meta,
        }
        if key == "pbb":
            parent_hx, parent_hz = build_bb(BBSpec(
                pspec.ell, pspec.m, list(pspec.A), list(pspec.B)
            ))
            dim = pspec.ell * pspec.m
            rec["shadow_parent_checks"] = {
                "top_x_equals_parent_hx": bool(np.array_equal(code.H[:dim, :code.n], parent_hx)),
                "bottom_z_equals_parent_hz": bool(np.array_equal(code.H[dim:, code.n:], parent_hz)),
            }
            if not all(rec["shadow_parent_checks"].values()):
                raise RuntimeError("PBB parent shadow equality failed")
        out[key] = {"circuit": circuit, "pure_rows": _pure_rows(circuit, code), "record": rec}
    return out


def decoder_matrices(circuit: stim.Circuit, pure_rows: np.ndarray) -> dict[str, Any]:
    mats = dem_to_matrices(circuit.detector_error_model(decompose_errors=False))
    h = sp.csr_matrix(mats.H)
    priors = np.clip(mats.priors, 1e-12, 1 - 1e-12)
    projected = h[pure_rows, :].tocsr()
    active = np.flatnonzero(np.asarray(projected.getnnz(axis=0)).reshape(-1))
    return {
        "H": h, "L": mats.L, "priors": priors, "pure_rows": pure_rows,
        "stage_H": projected[:, active].tocsr(), "active": active,
        "H_active": h[:, active].tocsr(), "L_active": mats.L[:, active],
        "stage_priors": priors[active],
    }


def make_decoder(data: dict[str, Any], candidate: dict[str, Any], stage: bool = False):
    matrix = data["stage_H"] if stage else data["H"]
    priors = data["stage_priors"] if stage else data["priors"]
    common = dict(
        error_channel=priors.tolist(), max_iter=candidate["max_iter"],
        bp_method="ms", ms_scaling_factor=candidate["alpha"],
        schedule="parallel", omp_thread_count=1,
    )
    if candidate["kind"] == "bplsd":
        from ldpc import BpLsdDecoder
        return BpLsdDecoder(
            matrix, **common, bits_per_step=candidate["bits_per_step"],
            lsd_method=candidate["lsd_method"], lsd_order=candidate["lsd_order"],
        )
    from ldpc import BpOsdDecoder
    return BpOsdDecoder(
        matrix, **common, osd_method=candidate["osd_method"],
        osd_order=candidate["osd_order"],
    )


def parity(matrix, vector: np.ndarray) -> np.ndarray:
    return np.asarray(matrix @ vector, dtype=np.uint8).reshape(-1) & 1


def decode_once(data: dict[str, Any], candidate: dict[str, Any], decoders, syndrome):
    if candidate["kind"] != "shadow":
        corr = decoders[0].decode(syndrome)
        return parity(data["L"], corr), parity(data["H"], corr)
    full, stage = decoders
    ca = stage.decode(syndrome[data["pure_rows"]])
    residual = syndrome ^ parity(data["H_active"], ca)
    prediction = parity(data["L_active"], ca)
    if residual.any():
        cb = full.decode(residual)
        prediction ^= parity(data["L"], cb)
        residual ^= parity(data["H"], cb)
    return prediction, syndrome ^ residual


def _worker(command_q, reply_q, circuit_text: str, pure_rows, candidate):
    circuit = stim.Circuit(circuit_text)
    data = decoder_matrices(circuit, np.asarray(pure_rows, dtype=np.int64))
    decoders = (make_decoder(data, candidate),)
    if candidate["kind"] == "shadow":
        decoders = (decoders[0], make_decoder(data, candidate, stage=True))
    reply_q.put({"kind": "ready"})
    while True:
        message = command_q.get()
        if message is None:
            return
        shot_id, syndrome = message
        started = time.perf_counter()
        try:
            prediction, reproduced = decode_once(data, candidate, decoders, syndrome)
            reply_q.put({
                "kind": "result", "shot_id": shot_id,
                "prediction": prediction.tolist(),
                "valid": bool(np.array_equal(reproduced, syndrome)),
                "latency_s": time.perf_counter() - started,
            })
        except BaseException as error:
            reply_q.put({
                "kind": "fatal", "shot_id": shot_id,
                "error_type": type(error).__name__, "error": str(error),
            })


@dataclass
class DecodeResult:
    prediction: np.ndarray | None
    valid: bool
    timeout: bool
    fatal: str | None
    latency_s: float
    restarted: bool


class SupervisedDecoder:
    """Persistent decoder worker; hard per-call deadline and restart on timeout."""
    def __init__(self, circuit: stim.Circuit, pure_rows: np.ndarray, candidate: dict[str, Any],
                 deadline_s: float = DEADLINE_S):
        self.circuit_text = str(circuit)
        self.pure_rows = pure_rows.tolist()
        self.candidate = candidate
        self.deadline_s = deadline_s
        self.startup_deadline_s = max(DEADLINE_S, deadline_s)
        self.context = mp.get_context("spawn")
        self.process = None
        self.command_q = None
        self.reply_q = None
        self.restarts = 0
        self._start()

    def _start(self) -> None:
        self.command_q = self.context.Queue()
        self.reply_q = self.context.Queue()
        self.process = self.context.Process(
            target=_worker,
            args=(self.command_q, self.reply_q, self.circuit_text,
                  self.pure_rows, self.candidate),
        )
        self.process.start()
        try:
            ready = self.reply_q.get(timeout=self.startup_deadline_s)
        except queue.Empty as error:
            self._terminate()
            raise TimeoutError("decoder initialization exceeded deadline") from error
        if ready.get("kind") != "ready":
            self._terminate()
            raise RuntimeError(f"decoder failed to initialize: {ready}")

    def _terminate(self) -> None:
        if self.process is not None and self.process.is_alive():
            self.process.terminate()
            self.process.join(timeout=5)
            if self.process.is_alive():
                self.process.kill()
                self.process.join(timeout=5)

    def decode(self, shot_id: int, syndrome: np.ndarray) -> DecodeResult:
        self.command_q.put((shot_id, syndrome.astype(np.uint8)))
        started = time.perf_counter()
        try:
            reply = self.reply_q.get(timeout=self.deadline_s)
        except queue.Empty:
            elapsed = time.perf_counter() - started
            self._terminate()
            self.restarts += 1
            self._start()
            return DecodeResult(None, False, True, None, elapsed, True)
        if reply.get("shot_id") != shot_id:
            raise RuntimeError("supervisor received unpaired shot id")
        if reply["kind"] == "fatal":
            return DecodeResult(None, False, False, reply["error"],
                                time.perf_counter() - started, False)
        return DecodeResult(
            np.asarray(reply["prediction"], dtype=np.uint8), bool(reply["valid"]),
            False, None, float(reply["latency_s"]), False,
        )

    def close(self) -> None:
        if self.process is not None and self.process.is_alive():
            self.command_q.put(None)
            self.process.join(timeout=5)
        self._terminate()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def sample_pairs(circuit: stim.Circuit, shots: int, seed: int):
    det, obs = circuit.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True,
    )
    return det.astype(np.uint8), obs.astype(np.uint8)


def validate_record_prefix(records: list[dict[str, Any]], expected_shots: int) -> int:
    shot_ids = [int(record["shot_id"]) for record in records]
    if shot_ids != list(range(len(records))):
        raise RuntimeError("checkpoint records are not a unique contiguous shot prefix")
    if len(records) > expected_shots:
        raise RuntimeError("checkpoint contains more shots than the frozen sample")
    return len(records)


def sample_sha256(det: np.ndarray, obs: np.ndarray) -> str:
    return hashlib.sha256(det.tobytes() + obs.tobytes()).hexdigest()


def checkpoint_binding(seed: int, det: np.ndarray, obs: np.ndarray,
                       circuit_sha256: str) -> dict[str, Any]:
    return {
        "protocol_sha256": protocol_hash(), "circuit_sha256": circuit_sha256,
        "seed": seed, "paired_sample_sha256": sample_sha256(det, obs),
        "expected_shots": int(len(det)),
    }


def validate_binding(entry: dict[str, Any], binding: dict[str, Any]) -> None:
    stored = entry.get("binding")
    if stored is not None and stored != binding:
        raise RuntimeError("checkpoint protocol/circuit/sample binding mismatch")
    validate_record_prefix(entry.get("records", []), binding["expected_shots"])


def evaluate_candidate(circuit, pure_rows, candidate, det, obs,
                       checkpoint: Callable | None = None,
                       initial_records: list[dict[str, Any]] | None = None,
                       supervisor_factory=SupervisedDecoder):
    records = list(initial_records or [])
    start = validate_record_prefix(records, len(det))
    with supervisor_factory(circuit, pure_rows, candidate) as supervisor:
        for i in range(start, len(det)):
            syndrome, actual = det[i], obs[i]
            decoded = supervisor.decode(i, syndrome)
            mismatch = decoded.prediction is None or bool(np.any(decoded.prediction != actual))
            operational = mismatch or not decoded.valid or decoded.timeout
            records.append({
                "shot_id": i, "logical_mismatch": mismatch,
                "invalid_correction": not decoded.valid, "timeout": decoded.timeout,
                "operational_failure": operational, "fatal": decoded.fatal,
                "latency_s": decoded.latency_s, "worker_restarted": decoded.restarted,
            })
            if checkpoint and (
                len(records) % CHECKPOINT_EVERY == 0 or len(records) == len(det)
            ):
                checkpoint(records)
            if decoded.fatal:
                break
    return records


def counts(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    failures = sum(r["operational_failure"] for r in records)
    lo, hi = clopper_pearson(failures, n) if n else (0.0, 1.0)
    return {
        "shots": n, "operational_failures": failures,
        "operational_failure_rate": failures / n if n else None,
        "ci95": [lo, hi],
        "logical_mismatches": sum(r["logical_mismatch"] for r in records),
        "invalid_corrections": sum(r["invalid_correction"] for r in records),
        "timeouts": sum(r["timeout"] for r in records),
        "worker_restarts": sum(r["worker_restarted"] for r in records),
        "fatal_errors": sum(r["fatal"] is not None for r in records),
    }


def one_sided_mcnemar(baseline: list[bool], challenger: list[bool]) -> dict[str, Any]:
    """Exact H0 equal discordant probability; alternative challenger fails less."""
    from scipy.stats import binomtest
    b_only = sum(b and not c for b, c in zip(baseline, challenger, strict=True))
    c_only = sum(c and not b for b, c in zip(baseline, challenger, strict=True))
    discordant = b_only + c_only
    pvalue = 1.0 if discordant == 0 else float(
        binomtest(b_only, discordant, 0.5, alternative="greater").pvalue
    )
    return {"baseline_only_failures": b_only, "challenger_only_failures": c_only,
            "discordant": discordant, "pvalue_one_sided": pvalue}


def holm(pvalues: dict[str, float], alpha: float = 0.05) -> dict[str, dict[str, Any]]:
    ordered = sorted(pvalues.items(), key=lambda item: (item[1], item[0]))
    m = len(ordered)
    rejected_prefix = True
    adjusted_running = 0.0
    out = {}
    for rank, (name, pvalue) in enumerate(ordered, start=1):
        threshold = alpha / (m - rank + 1)
        rejected_prefix = rejected_prefix and pvalue <= threshold
        adjusted_running = max(adjusted_running, (m - rank + 1) * pvalue)
        out[name] = {
            "raw_p": pvalue, "holm_rank": rank, "threshold": threshold,
            "holm_adjusted_p": min(1.0, adjusted_running),
            "rejected": bool(rejected_prefix),
        }
    return out


def balanced_orders(candidate_ids: list[str], blocks: int, seed: int) -> list[list[str]]:
    rng = np.random.default_rng(seed)
    base = list(candidate_ids)
    rng.shuffle(base)
    orders = []
    for block in range(blocks):
        shift = block % len(base)
        order = base[shift:] + base[:shift]
        if (block // len(base)) % 2:
            order = list(reversed(order))
        orders.append(order)
    return orders


def bootstrap_p95_ratio(base: np.ndarray, challenger: np.ndarray, seed: int,
                        reps: int = BOOTSTRAP_REPS) -> dict[str, Any]:
    if len(base) != len(challenger):
        raise ValueError("paired latency arrays must have equal length")
    if len(base) == 0:
        return {"ratio": None, "ci95": [None, None], "reps": reps}
    rng = np.random.default_rng(seed)
    ratios = np.empty(reps, dtype=float)
    n = len(base)
    for i in range(reps):
        idx = rng.integers(0, n, size=n)
        ratios[i] = np.percentile(challenger[idx], 95) / np.percentile(base[idx], 95)
    ratio = float(np.percentile(challenger, 95) / np.percentile(base, 95))
    return {"ratio": ratio, "ci95": [float(np.percentile(ratios, 2.5)),
                                      float(np.percentile(ratios, 97.5))], "reps": reps}

def heavy_processes() -> list[str]:
    import subprocess
    rows = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,command="], check=True, text=True,
        capture_output=True,
    ).stdout.splitlines()
    processes: dict[int, tuple[int, str]] = {}
    for row in rows:
        fields = row.strip().split(maxsplit=2)
        if len(fields) != 3:
            continue
        try:
            processes[int(fields[0])] = (int(fields[1]), fields[2])
        except ValueError:
            continue
    self_pid = os.getpid()
    descendants = {self_pid}
    changed = True
    while changed:
        changed = False
        for pid, (ppid, _) in processes.items():
            if ppid in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    heavy = []
    for pid, (_, command) in processes.items():
        if pid in descendants:
            continue
        lower = command.lower()
        is_qec_experiment = (
            "python" in lower
            and ("experiments/exp" in lower or (
                ("qec-codesign" in lower or "math/qec" in lower)
                and "pytest" not in lower
                and "oh-my-pi" not in lower and "omp" not in lower
            ))
        )
        is_known_solver = any(token in lower for token in (
            "kissat", "cadical", "cryptominisat", "glucose", "lingeling",
            "exp024_", "exp025_", "schedule_opt", "cat_extraction",
        ))
        if is_qec_experiment or is_known_solver:
            heavy.append(f"{pid}:{command}")
    return sorted(heavy)


def timing_gate() -> dict[str, Any]:
    load = list(os.getloadavg())
    ncpu = os.cpu_count() or 1
    heavy = heavy_processes()
    return {
        "load_average": load, "cpu_count": ncpu,
        "max_load": MAX_IDLE_LOAD, "heavy_processes": heavy,
        "passed": load[0] <= MAX_IDLE_LOAD and not heavy,
    }


def state_template(circuits) -> dict[str, Any]:
    return {
        "status": "running", "artifact_route": "partial_runs",
        "protocol": protocol(), "protocol_sha256": protocol_hash(),
        "environment": {
            "python": platform.python_version(), "platform": platform.platform(),
            "stim": stim.__version__, "ldpc": importlib.metadata.version("ldpc"),
            "cpu_count": os.cpu_count(),
        },
        "circuits": {key: value["record"] for key, value in circuits.items()},
        "development": {}, "heldout": {}, "timing": {}, "wall_s": 0.0,
        "heldout_gate_sessions": [], "quarantined_sessions": [],
    }


def save_partial(state):
    if canonical_route(True, False) != "partial_runs":
        raise RuntimeError("partial artifact router failed")
    write_json(PARTIAL, state)


def run_dev(state, circuits):
    det, obs = sample_pairs(circuits["pbb"]["circuit"], DEV_SHOTS, DEV_SEED)
    binding = checkpoint_binding(
        DEV_SEED, det, obs, circuits["pbb"]["record"]["circuit_sha256"]
    )
    for candidate in CANDIDATES:
        cid = candidate["id"]
        entry = state["development"].get(cid, {})
        validate_binding(entry, binding)
        initial = entry.get("records", [])
        if len(initial) == DEV_SHOTS:
            continue
        def checkpoint(records, cid=cid):
            state["development"][cid] = {
                "binding": binding, "records": records, "counts": counts(records),
            }
            save_partial(state)
        records = evaluate_candidate(
            circuits["pbb"]["circuit"], circuits["pbb"]["pure_rows"], candidate,
            det, obs, checkpoint, initial_records=initial,
        )
        if any(r["fatal"] for r in records):
            raise RuntimeError(f"fatal development failure for {cid}")




def require_claiming_gate(state: dict[str, Any], phase: str) -> dict[str, Any]:
    gate = {**timing_gate(), "time_ns": time.time_ns()}
    snapshot = {
        **gate, "phase": phase,
        "session_index": len(state.setdefault("heldout_gate_sessions", [])),
        "start_gate": gate, "rechecks": [], "aborted": not gate["passed"],
        "claiming": gate["passed"],
    }
    state["heldout_gate_sessions"].append(snapshot)
    if not gate["passed"]:
        state.setdefault("quarantined_sessions", []).append(dict(snapshot))
    save_partial(state)
    if not gate["passed"]:
        raise RuntimeError(f"{phase} isolation gate refused run: {snapshot}")
    return snapshot


def recheck_claiming_gate(state: dict[str, Any], session_gate: dict[str, Any],
                          phase: str) -> dict[str, Any]:
    check = {**timing_gate(), "time_ns": time.time_ns()}
    session_gate.setdefault("rechecks", []).append(check)
    if check["passed"]:
        return check
    session_gate["claiming"] = False
    session_gate["aborted"] = True
    session_gate["aborted_reason"] = check
    state.setdefault("quarantined_sessions", []).append(dict(session_gate))
    save_partial(state)
    raise RuntimeError(f"{phase} continuous isolation gate failed: {check}")


def run_heldout(state, circuits):
    session_gate = require_claiming_gate(state, "heldout")
    for code_index, code_key in enumerate(CODE_KEYS):
        seed = HELDOUT_SEED_BASE + code_index
        det, obs = sample_pairs(circuits[code_key]["circuit"], HELDOUT_SHOTS, seed)
        binding = checkpoint_binding(
            seed, det, obs, circuits[code_key]["record"]["circuit_sha256"]
        )
        state["heldout"].setdefault(code_key, {})
        for candidate in CANDIDATES:
            cid = candidate["id"]
            entry = state["heldout"][code_key].get(cid, {})
            validate_binding(entry, binding)
            initial = entry.get("records", [])
            if len(initial) == HELDOUT_SHOTS:
                continue
            segments = list(entry.get("gate_segments", []))
            session_index = session_gate["session_index"]
            cursor = {"shot": len(initial), "gate": dict(session_gate["start_gate"])}
            def checkpoint(records, cid=cid, segments=segments, cursor=cursor):
                end_gate = recheck_claiming_gate(state, session_gate, "heldout")
                segments.append({
                    "session_index": session_gate["session_index"],
                    "start_shot": cursor["shot"],
                    "end_shot_exclusive": len(records),
                    "start_gate": cursor["gate"],
                    "end_gate": end_gate,
                    "passed": True,
                })
                cursor["shot"] = len(records)
                cursor["gate"] = end_gate
                state["heldout"][code_key][cid] = {
                    "binding": binding, "records": records, "counts": counts(records),
                    "gate_segments": segments,
                }
                save_partial(state)
            evaluate_candidate(
                circuits[code_key]["circuit"], circuits[code_key]["pure_rows"],
                candidate, det, obs, checkpoint, initial_records=initial,
            )

def recheck_timing_gate(state: dict[str, Any], code_key: str,
                        shot_boundary: int) -> dict[str, Any]:
    snapshot = {
        **timing_gate(), "code": code_key, "shot_boundary": shot_boundary,
        "time_ns": time.time_ns(),
    }
    state["timing_run"].setdefault("rechecks", []).append(snapshot)
    if snapshot["passed"]:
        return snapshot
    state["timing_run"]["aborted"] = True
    state["timing_run"]["aborted_reason"] = snapshot
    state.setdefault("quarantined_sessions", []).append(dict(state["timing_run"]))
    state["timing"] = {}
    save_partial(state)
    raise RuntimeError(f"timing continuous isolation gate failed: {snapshot}")


def run_timing(state, circuits):
    gate = {**timing_gate(), "time_ns": time.time_ns()}
    state["timing_gate"] = gate
    if not gate["passed"]:
        state.setdefault("quarantined_sessions", []).append({
            **gate, "phase": "timing", "aborted": True, "time_ns": time.time_ns(),
        })
    save_partial(state)
    if not gate["passed"]:
        raise RuntimeError(f"timing isolation gate refused run: {gate}")
    # Intentionally non-resumable for drift safety: no prior timing records are
    # reused. A complete pass has one run id, one gate snapshot, and all 5x2 rows.
    state["timing"] = {}
    run_id = f"pid-{os.getpid()}-ns-{time.time_ns()}"
    state["timing_run"] = {
        "id": run_id, "intentionally_non_resumable": True,
        "interrupted_pass_policy": "discard entire timing pass and rerun in isolation",
        "gate": gate,
    }
    save_partial(state)
    for code_index, code_key in enumerate(CODE_KEYS):
        seed = TIMING_SEED_BASE + code_index
        det, obs = sample_pairs(circuits[code_key]["circuit"],
                                TIMING_WARMUP + TIMING_SHOTS, seed)
        binding = checkpoint_binding(
            seed, det, obs, circuits[code_key]["record"]["circuit_sha256"]
        )
        supervisors = {
            c["id"]: SupervisedDecoder(circuits[code_key]["circuit"],
                                        circuits[code_key]["pure_rows"], c)
            for c in CANDIDATES
        }
        candidates = {c["id"]: c for c in CANDIDATES}
        ids = list(candidates)
        blocks = (TIMING_WARMUP + TIMING_SHOTS + TIMING_BLOCK - 1) // TIMING_BLOCK
        orders = balanced_orders(ids, blocks, seed)
        records = {cid: [] for cid in ids}
        gate_segments = []
        block_start_gate = None
        block_start_shot = 0
        try:
            for shot_id, (syndrome, actual) in enumerate(zip(det, obs, strict=True)):
                order = orders[shot_id // TIMING_BLOCK]
                if shot_id % TIMING_BLOCK == 0:
                    block_start_shot = shot_id
                    block_start_gate = recheck_timing_gate(state, code_key, shot_id)
                for cid in order:
                    decoded = supervisors[cid].decode(shot_id, syndrome)
                    records[cid].append({
                        "shot_id": shot_id, "warmup": shot_id < TIMING_WARMUP,
                        "latency_s": min(decoded.latency_s, DEADLINE_S),
                        "timeout": decoded.timeout, "invalid_correction": not decoded.valid,
                        "logical_mismatch": decoded.prediction is None or bool(
                            np.any(decoded.prediction != actual)),
                        "worker_restarted": decoded.restarted,
                    })
                if (
                    (shot_id + 1) % TIMING_BLOCK == 0
                    or shot_id + 1 == len(det)
                ):
                    end_gate = recheck_timing_gate(state, code_key, shot_id + 1)
                    gate_segments.append({
                        "start_shot": block_start_shot,
                        "end_shot_exclusive": shot_id + 1,
                        "start_gate": block_start_gate,
                        "end_gate": end_gate,
                        "passed": True,
                    })
        finally:
            for supervisor in supervisors.values():
                supervisor.close()
        state["timing"][code_key] = {
            "binding": binding, "timing_run_id": run_id,
            "balanced_orders": orders, "records": records,
            "gate_segments": gate_segments,
        }
        save_partial(state)


def gate_segments_cover(
    segments: list[dict[str, Any]],
    expected_shots: int,
    max_span: int,
    sessions: list[dict[str, Any]] | None = None,
    phase: str | None = None,
) -> bool:
    if not segments:
        return False
    cursor = 0
    for segment in segments:
        start = segment.get("start_shot")
        end = segment.get("end_shot_exclusive")
        valid = (
            segment.get("passed") is True
            and start == cursor
            and isinstance(end, int)
            and cursor < end <= expected_shots
            and end - cursor <= max_span
            and valid_isolation_snapshot(segment.get("start_gate", {}))
            and valid_isolation_snapshot(segment.get("end_gate", {}))
        )
        if sessions is not None:
            index = segment.get("session_index")
            valid_session = isinstance(index, int) and 0 <= index < len(sessions)
            if valid_session:
                session = sessions[index]
                linked_gates = [
                    session.get("start_gate", {}),
                    *session.get("rechecks", []),
                ]
                valid_session = (
                    session.get("phase") == phase
                    and segment.get("start_gate") in linked_gates
                    and segment.get("end_gate") in linked_gates
                )
            valid = valid and valid_session
        if not valid:
            return False
        cursor = end
    return cursor == expected_shots


def expected_binding(state: dict[str, Any], code_key: str, seed: int,
                     shots: int, sample_hash: str) -> dict[str, Any]:
    return {
        "protocol_sha256": protocol_hash(),
        "circuit_sha256": state["circuits"][code_key]["circuit_sha256"],
        "seed": seed, "paired_sample_sha256": sample_hash,
        "expected_shots": shots,
    }


def valid_isolation_snapshot(snapshot: dict[str, Any]) -> bool:
    load = snapshot.get("load_average", [])
    return bool(
        snapshot.get("passed")
        and snapshot.get("max_load") == MAX_IDLE_LOAD
        and load and load[0] <= MAX_IDLE_LOAD
        and not snapshot.get("heavy_processes")
    )


def validate_full_coverage(state: dict[str, Any], circuits: dict[str, Any] | None = None) -> None:
    expected_ids = {candidate["id"] for candidate in CANDIDATES}
    if circuits is not None:
        for code_key in CODE_KEYS:
            if (
                state["circuits"][code_key]["circuit_sha256"]
                != circuits[code_key]["record"]["circuit_sha256"]
            ):
                raise RuntimeError(f"stored circuit hash differs from regenerated {code_key}")
    if set(state.get("development", {})) != expected_ids:
        raise RuntimeError("development candidate coverage incomplete")
    dev_bindings = [entry.get("binding", {}) for entry in state["development"].values()]
    if len({json.dumps(binding, sort_keys=True) for binding in dev_bindings}) != 1:
        raise RuntimeError("development candidates do not share one frozen sample binding")
    dev_expected = expected_binding(
        state, "pbb", DEV_SEED, DEV_SHOTS,
        dev_bindings[0].get("paired_sample_sha256", ""),
    )
    if dev_bindings[0] != dev_expected:
        raise RuntimeError("development binding violates frozen contract")
    for cid, entry in state["development"].items():
        validate_record_prefix(entry.get("records", []), DEV_SHOTS)
        if len(entry["records"]) != DEV_SHOTS or entry["binding"] != dev_expected:
            raise RuntimeError(f"development coverage/binding incomplete for {cid}")
    sessions = state.get("heldout_gate_sessions", [])
    if not sessions:
        raise RuntimeError("heldout gate provenance missing")
    for code_index, code_key in enumerate(CODE_KEYS):
        entries = state.get("heldout", {}).get(code_key, {})
        if set(entries) != expected_ids:
            raise RuntimeError(f"heldout candidate coverage incomplete for {code_key}")
        bindings = [entry.get("binding", {}) for entry in entries.values()]
        if len({json.dumps(binding, sort_keys=True) for binding in bindings}) != 1:
            raise RuntimeError(f"heldout candidates do not share one sample for {code_key}")
        held_expected = expected_binding(
            state, code_key, HELDOUT_SEED_BASE + code_index, HELDOUT_SHOTS,
            bindings[0].get("paired_sample_sha256", ""),
        )
        if bindings[0] != held_expected:
            raise RuntimeError(f"heldout frozen binding violated for {code_key}")
        for cid, entry in entries.items():
            segments = entry.get("gate_segments", [])
            provenance_valid = gate_segments_cover(
                segments, HELDOUT_SHOTS, CHECKPOINT_EVERY, sessions, "heldout"
            )
            if (
                len(entry["records"]) != HELDOUT_SHOTS
                or entry["binding"] != held_expected
                or not provenance_valid
            ):
                raise RuntimeError(f"heldout coverage/provenance incomplete for {code_key}/{cid}")
    run_id = state.get("timing_run", {}).get("id")
    if not run_id:
        raise RuntimeError("complete isolated timing run missing")
    timing_run = state["timing_run"]
    if (
        timing_run.get("aborted")
        or not valid_isolation_snapshot(timing_run.get("gate", {}))
        or not timing_run.get("rechecks")
        or not all(valid_isolation_snapshot(check)
                   for check in timing_run["rechecks"])
    ):
        raise RuntimeError("timing continuous isolation provenance incomplete")
    for code_index, code_key in enumerate(CODE_KEYS):
        entry = state.get("timing", {}).get(code_key, {})
        if entry.get("timing_run_id") != run_id:
            raise RuntimeError("timing rows do not share one isolated pass")
        if set(entry.get("records", {})) != expected_ids:
            raise RuntimeError(f"timing candidate coverage incomplete for {code_key}")
        binding = entry.get("binding", {})
        timing_expected = expected_binding(
            state, code_key, TIMING_SEED_BASE + code_index,
            TIMING_WARMUP + TIMING_SHOTS,
            binding.get("paired_sample_sha256", ""),
        )
        if binding != timing_expected:
            raise RuntimeError(f"timing frozen binding violated for {code_key}")
        if not gate_segments_cover(
            entry.get("gate_segments", []),
            TIMING_WARMUP + TIMING_SHOTS,
            TIMING_BLOCK,
        ):
            raise RuntimeError(f"timing gate-segment coverage incomplete for {code_key}")
        for cid, records in entry["records"].items():
            validate_record_prefix(records, TIMING_WARMUP + TIMING_SHOTS)
            if len(records) != TIMING_WARMUP + TIMING_SHOTS:
                raise RuntimeError(f"timing coverage incomplete for {code_key}/{cid}")
    if not valid_isolation_snapshot(state.get("timing_gate", {})):
        raise RuntimeError("timing isolation gate did not pass")


def report_text(result: dict[str, Any]) -> str:
    verdict = result["verdict"]
    reason = (
        f"{len(verdict['joint_winners'])} challenger(s) met all joint gates"
        if verdict["classification"] == "POSITIVE"
        else "no challenger met accuracy significance and paired latency gates"
    )
    lines = [
        f"{verdict['classification']} — {reason} within the finite EXP-026 v2 set.",
        "", "# EXP-026 decoder co-design v2", "", "## Protocol revision", "",
        "V1 produced **0/58 production rows** before revision. Its conservative "
        "baseline-equivalent lower bound exceeded **216 core-hours**; revision "
        "therefore preceded production data. V1 smoke and canonical-write guard "
        "remain in `results/partial_runs/`.", "", "## Evidence", "",
    ]
    for code_key in CODE_KEYS:
        lines += [
            f"### {code_key}", "",
            "| decoder | failures/shots | rate [95% CP] | timeouts | invalid | p50/p95/p99 ms |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for candidate in CANDIDATES:
            cid = candidate["id"]
            acc, tim = result["heldout"][code_key][cid], result["timing"][code_key][cid]
            lines.append(
                f"| `{cid}` | {acc['operational_failures']}/{acc['shots']} | "
                f"{acc['operational_failure_rate']:.6g} "
                f"[{acc['ci95'][0]:.6g},{acc['ci95'][1]:.6g}] | {acc['timeouts']} | "
                f"{acc['invalid_corrections']} | {tim['p50_ms']:.4g}/"
                f"{tim['p95_ms']:.4g}/{tim['p99_ms']:.4g} |"
            )
        lines.append("")
    lines += [
        "## PBB joint gates", "",
        "| challenger | McNemar p | Holm reject | p95 ratio 95% CI |",
        "|---|---:|---:|---:|",
    ]
    for candidate in CANDIDATES[1:]:
        cid = candidate["id"]
        ratio = result["timing"]["pbb"][cid]["paired_p95_ratio"]
        lines.append(
            f"| `{cid}` | {result['mcnemar'][cid]['pvalue_one_sided']:.6g} | "
            f"{result['holm'][cid]['rejected']} | {ratio['ratio']:.4g} "
            f"[{ratio['ci95'][0]:.4g},{ratio['ci95'][1]:.4g}] |"
        )
    lines += [
        "", "## Reproduction (resumes validated partial checkpoints by default)", "",
        "```bash",
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
        "VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=src "
        ".venv/bin/python experiments/exp026_decoder_codesign.py --phase all",
        "```", "",
        "`--no-resume` requests an explicit clean rerun. Held-out gate sessions and "
        "the one-pass timing isolation gate are stored in the raw artifact.",
        f" Timing isolation gate: `{result['timing_gate']}`.",
        f"Total wall: {result['wall_s']:.3f} s. Exact schedule claim scopes and "
        "full slot maps are in the processed JSON.",
    ]
    return "\n".join(lines) + "\n"


def replace_canonical_bundle(payloads: dict[Path, str],
                             replacer: Callable[[Path, Path], None] = os.replace) -> None:
    staged: list[tuple[Path, Path]] = []
    backups: dict[Path, Path | None] = {}
    replaced: list[Path] = []
    try:
        for path, content in payloads.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.stage.", dir=path.parent)
            temp_path = Path(temp_name)
            with os.fdopen(fd, "w") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            staged.append((temp_path, path))
            if path.exists():
                backup_fd, backup_name = tempfile.mkstemp(
                    prefix=f".{path.name}.backup.", dir=path.parent
                )
                backup = Path(backup_name)
                with os.fdopen(backup_fd, "wb") as handle:
                    handle.write(path.read_bytes())
                    handle.flush()
                    os.fsync(handle.fileno())
                backups[path] = backup
            else:
                backups[path] = None
        # Rendering, staging, and rollback backups all complete before replacement.
        for temp_path, path in staged:
            replacer(temp_path, path)
            replaced.append(path)
    except BaseException:
        for path in reversed(replaced):
            backup = backups[path]
            if backup is None:
                path.unlink(missing_ok=True)
            else:
                os.replace(backup, path)
                backups[path] = None
        raise
    finally:
        for temp_path, _ in staged:
            temp_path.unlink(missing_ok=True)
        for backup in backups.values():
            if backup is not None:
                backup.unlink(missing_ok=True)


def finalize(state, circuits):
    validate_full_coverage(state, circuits)
    analysis = {"heldout": {}, "timing": {}, "mcnemar": {}, "holm": {}, "verdict": {}}
    for code_key in CODE_KEYS:
        analysis["heldout"][code_key] = {
            cid: entry["counts"] for cid, entry in state["heldout"][code_key].items()
        }
        analysis["timing"][code_key] = {}
        for cid, recs in state["timing"][code_key]["records"].items():
            timed = [r for r in recs if not r["warmup"]]
            lat = np.asarray([r["latency_s"] for r in timed]) * 1000
            analysis["timing"][code_key][cid] = {
                "timed_shots": len(timed), "timeouts": sum(r["timeout"] for r in timed),
                "p50_ms": float(np.percentile(lat, 50)),
                "p95_ms": float(np.percentile(lat, 95)),
                "p99_ms": float(np.percentile(lat, 99)),
            }
    base_records = state["heldout"]["pbb"]["baseline_bposd"]["records"]
    base_fail = [r["operational_failure"] for r in base_records]
    pvalues = {}
    for index, candidate in enumerate(CANDIDATES[1:], start=1):
        cid = candidate["id"]
        ch_fail = [r["operational_failure"] for r in state["heldout"]["pbb"][cid]["records"]]
        result = one_sided_mcnemar(base_fail, ch_fail)
        analysis["mcnemar"][cid] = result
        pvalues[cid] = result["pvalue_one_sided"]
    analysis["holm"] = holm(pvalues)
    base_timed = [r for r in state["timing"]["pbb"]["records"]["baseline_bposd"]
                  if not r["warmup"]]
    base_lat = np.asarray([r["latency_s"] for r in base_timed])
    winners = []
    for index, candidate in enumerate(CANDIDATES[1:], start=1):
        cid = candidate["id"]
        ch_timed = [r for r in state["timing"]["pbb"]["records"][cid]
                    if not r["warmup"]]
        ch_lat = np.asarray([r["latency_s"] for r in ch_timed])
        ratio = bootstrap_p95_ratio(base_lat, ch_lat, BOOTSTRAP_SEED + index)
        analysis["timing"]["pbb"][cid]["paired_p95_ratio"] = ratio
        lower_rate = (analysis["heldout"]["pbb"][cid]["operational_failure_rate"] <
                      analysis["heldout"]["pbb"]["baseline_bposd"]["operational_failure_rate"])
        joint = lower_rate and analysis["holm"][cid]["rejected"] and ratio["ci95"][1] < 1
        if joint:
            winners.append(cid)
    analysis["verdict"] = {
        "classification": "POSITIVE" if winners else "NEGATIVE",
        "scope": "finite EXP-026 v2 candidate set on PBB",
        "joint_winners": winners,
    }
    state["analysis"] = analysis
    state["status"] = "complete"
    state["artifact_route"] = canonical_route(True, True)
    if state["artifact_route"] != "canonical":
        raise RuntimeError("complete artifact did not route canonical")
    processed = {
        "experiment": "EXP-026", "protocol": state["protocol"],
        "protocol_sha256": state["protocol_sha256"], "circuits": state["circuits"],
        "development": {cid: entry["counts"] for cid, entry in state["development"].items()},
        **analysis, "timing_gate": state["timing_gate"],
        "heldout_gate_sessions": state["heldout_gate_sessions"],
        "wall_s": state["wall_s"],
        "machine_readable_verdict": analysis["verdict"]["classification"],
    }
    replace_canonical_bundle({
        RAW: json.dumps(state, indent=2, sort_keys=True),
        PROCESSED: json.dumps(processed, indent=2, sort_keys=True),
        REPORT: report_text(processed),
    })


def write_report(result):
    atomic_write_text(REPORT, report_text(result))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("all", "dev", "accuracy", "timing", "smoke"),
                        default="all")
    parser.add_argument("--no-resume", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    started = time.perf_counter()
    circuits = build_circuits()
    if args.no_resume or not PARTIAL.exists():
        state = state_template(circuits)
    else:
        state = json.loads(PARTIAL.read_text())
        if state["protocol_sha256"] != protocol_hash():
            raise RuntimeError("partial checkpoint protocol mismatch")
    if args.phase == "smoke":
        det, obs = sample_pairs(circuits["pbb"]["circuit"], 1, DEV_SEED)
        for candidate in CANDIDATES:
            records = evaluate_candidate(circuits["pbb"]["circuit"],
                                         circuits["pbb"]["pure_rows"], candidate,
                                         det, obs)
            if records[0]["fatal"]:
                raise RuntimeError(records[0]["fatal"])
        state["status"] = "smoke_complete"
        state["wall_s"] = time.perf_counter() - started
        save_partial(state)
        return
    if args.phase in ("all", "dev"):
        run_dev(state, circuits)
    if args.phase in ("all", "accuracy"):
        run_heldout(state, circuits)
    if args.phase in ("all", "timing"):
        run_timing(state, circuits)
    state["wall_s"] += time.perf_counter() - started
    if args.phase == "all":
        finalize(state, circuits)
    else:
        state["status"] = f"partial_{args.phase}"
        save_partial(state)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        write_json(QUARANTINE, {
            "experiment": "EXP-026", "artifact_route": canonical_route(False, False),
            "error_type": type(error).__name__, "error": str(error),
            "protocol_sha256": protocol_hash(), "time_ns": time.time_ns(),
        })
        raise
