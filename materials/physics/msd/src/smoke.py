"""Bounded Gate-A construction, static, postselection, and pilot smoke run."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import uuid

import time
import stim

try:
    from .code_832 import LX, Z_STABS
    from .encoder import ENCODER_LAYERS, PLUS_INPUTS
except ImportError:  # pragma: no cover - exercised by direct CLI invocation.
    from code_832 import LX, Z_STABS
    from encoder import ENCODER_LAYERS, PLUS_INPUTS


try:
    from .estimate import sample_gate_a
    from .gate_a import build_gate_a_circuit, circuit_metadata
except ImportError:  # pragma: no cover - exercised by direct CLI invocation.
    from estimate import sample_gate_a
    from gate_a import build_gate_a_circuit, circuit_metadata


ROOT = Path(__file__).resolve().parents[1]


def _encoder_stabilizer_check(*, shots: int, seed: int) -> dict[str, object]:
    circuit = stim.Circuit()
    circuit.append("R", list(range(8)))
    circuit.append("RX", list(PLUS_INPUTS))
    for layer in ENCODER_LAYERS:
        circuit.append("CX", [q for pair in layer for q in pair])

    def pauli(kind: str, support: set[int]) -> stim.PauliString:
        text = ["_"] * 8
        for q in support:
            text[q] = kind
        return stim.PauliString("".join(text))

    checks = [pauli("X", set(range(8)))]
    checks.extend(pauli("X", support) for support in LX)
    checks.append(pauli("Z", set(range(8))))
    checks.extend(pauli("Z", support) for support in Z_STABS[1:])
    for check in checks:
        circuit.append("MPP", check)
    samples = circuit.compile_sampler(seed=seed).sample(shots)
    nonzero = int(samples.any(axis=1).sum())
    if nonzero:
        raise AssertionError(f"encoder stabilizer violation in {nonzero} shots")
    return {
        "shots": shots,
        "checks": len(checks),
        "nonzero_shots": nonzero,
        "status": "PASS",
        "evidence": "REPRODUCED",
    }


def run_smoke(*, pilot_shots: int = 1000, seed: int = 20260627) -> dict[str, object]:
    """Run deterministic checks plus one bounded noisy pilot."""

    encoder_static = _encoder_stabilizer_check(shots=256, seed=seed)

    ideal_circuit = build_gate_a_circuit(0.0)
    ideal_detectors, ideal_observables = ideal_circuit.compile_detector_sampler(
        seed=seed
    ).sample(256, separate_observables=True)
    ideal_detector_nonzero = int(ideal_detectors.any(axis=1).sum())
    ideal_observable_nonzero = int(ideal_observables.any(axis=1).sum())
    if ideal_detector_nonzero or ideal_observable_nonzero:
        raise AssertionError(
            "ideal circuit is not deterministic: "
            f"detectors={ideal_detector_nonzero}, observables={ideal_observable_nonzero}"
        )

    # A deliberately stronger noise probe verifies that detector postselection
    # is active, independently of whether the 1e-3 pilot happens to reject.
    probe = sample_gate_a(0.1, 512, seed=seed + 1)
    if probe.rejected <= 0:
        raise AssertionError("postselection probe produced no rejected shots")

    start = time.perf_counter()
    pilot = sample_gate_a(1e-3, pilot_shots, seed=seed + 2)
    elapsed = time.perf_counter() - start
    pilot_dict = pilot.as_dict()
    pilot_dict["wall_seconds"] = elapsed
    pilot_dict["wall_seconds_per_shot"] = elapsed / pilot_shots

    ideal_text = str(ideal_circuit)
    noisy_circuit = build_gate_a_circuit(1e-3)
    noisy_text = str(noisy_circuit)
    return {
        "schema": "msd.gate_a.smoke.v1",
        "evidence": "REPRODUCED",
        "run": {
            "script": "src/smoke.py",
            "requested_priority": "nice -n 10",
            "single_core": True,
            "pilot_shots": pilot_shots,
            "seed": seed,
        },
        "encoder_static": encoder_static,
        "ideal_static": {
            "shots": 256,
            "num_qubits": ideal_circuit.num_qubits,
            "num_ticks": ideal_circuit.num_ticks,
            "num_detectors": ideal_circuit.num_detectors,
            "num_observables": ideal_circuit.num_observables,
            "detector_nonzero_shots": ideal_detector_nonzero,
            "observable_nonzero_shots": ideal_observable_nonzero,
            "status": "PASS",
        },
        "postselection_probe": {
            **probe.as_dict(),
            "status": "PASS",
            "note": "p=0.1 probe only checks that explicit detectors reject faults",
        },
        "pilot": {
            **pilot_dict,
            "evidence": "NUMERICAL",
        },
        "circuit": {
            **circuit_metadata(1e-3),
            "ideal_sha256": hashlib.sha256(ideal_text.encode()).hexdigest(),
            "noisy_sha256": hashlib.sha256(noisy_text.encode()).hexdigest(),
        },
        "artifact_files": ["summary.json", "circuit_p0.stim", "circuit_p1e-3.stim"],
    }


def _default_artifact_dir(summary: dict[str, object]) -> Path:
    digest = str(summary["circuit"]["noisy_sha256"])[:12]  # type: ignore[index]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_uuid = uuid.uuid4().hex[:8]
    return ROOT / "campaigns" / f"{stamp}_{run_uuid}_{digest}"



def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-shots", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260627)
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=None,
        help="write an immutable smoke snapshot here (default: campaigns/<UTC>_<uuid>_<hash>)",
    )
    args = parser.parse_args()
    if args.pilot_shots <= 0:
        parser.error("--pilot-shots must be positive")

    summary = run_smoke(pilot_shots=args.pilot_shots, seed=args.seed)
    artifact_dir = args.artifact_dir or _default_artifact_dir(summary)
    artifact_dir.mkdir(parents=True, exist_ok=False)
    (artifact_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (artifact_dir / "circuit_p0.stim").write_text(str(build_gate_a_circuit(0.0)))
    (artifact_dir / "circuit_p1e-3.stim").write_text(str(build_gate_a_circuit(1e-3)))
    print(json.dumps(summary, indent=2))
    print(f"artifact_dir={artifact_dir}")


if __name__ == "__main__":
    main()
