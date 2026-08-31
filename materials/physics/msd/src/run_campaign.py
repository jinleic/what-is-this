"""Run a finite Gate-A p-grid and freeze the exact circuits/results."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import math

import json
from pathlib import Path
import time
import uuid

try:
    from .estimate import sample_gate_a, summarize
    from .gate_a import build_gate_a_circuit, circuit_metadata
except ImportError:  # pragma: no cover - exercised by direct CLI invocation.
    from estimate import sample_gate_a, summarize
    from gate_a import build_gate_a_circuit, circuit_metadata


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_P_VALUES = (5e-4, 1e-3, 3e-3, 1e-2)


def _token(p: float) -> str:
    return f"{p:.0e}".replace("+", "p").replace("-", "m")


def run_campaign(
    p_values: tuple[float, ...],
    shots: int,
    *,
    seed: int = 20260627,
    label: str = "gate_a",
) -> tuple[dict[str, object], dict[float, str]]:
    """Execute exactly one finite sample batch for each requested p."""

    if shots <= 0:
        raise ValueError(f"shots must be positive, got {shots!r}")
    if not p_values:
        raise ValueError("at least one physical error rate is required")
    if not label:
        raise ValueError("label must be non-empty")

    circuits = {p: build_gate_a_circuit(p) for p in p_values}
    estimates = []
    circuit_text = {}
    point_wall_seconds: dict[str, float] = {}
    for index, p in enumerate(p_values):
        circuit = circuits[p]
        circuit_text[p] = str(circuit)
        start = time.perf_counter()
        estimate = sample_gate_a(p, shots, seed=seed + index)
        point_wall_seconds[_token(p)] = time.perf_counter() - start
        estimates.append(estimate)

    result = summarize(estimates)
    fit = result["quadratic_fit"]
    coefficient = fit.get("coefficient_c") if isinstance(fit, dict) else None
    if coefficient is None:
        verdict = "INCONCLUSIVE"
    elif 150.0 <= float(coefficient) <= 600.0:
        verdict = "PASS"
    else:
        verdict = "REFUTE_FOR_THIS_RECONSTRUCTION"
    prime_result: dict[str, object] | None = None
    if label == "gate_a_prime":
        prime_row = next(
            (
                row
                for row in estimates
                if math.isclose(row.physical_error_rate, 1e-4, abs_tol=1e-15)
            ),
            None,
        )
        if prime_row is None or prime_row.logical_error_rate is None:
            prime_result = {
                "verdict": "INCONCLUSIVE",
                "reason": "missing p=1e-4 row or zero accepted shots",
                "evidence": "NUMERICAL",
            }
        else:
            low_p = prime_row.logical_error_rate
            c_eff = low_p / (1e-4 * 1e-4)
            prime_result = {
                "verdict": (
                    "REPRODUCTION_COMPATIBLE_AT_LOW_P"
                    if 1.5e-6 <= low_p <= 6e-6
                    else "SURROGATE_REFUTATION_AT_LOW_P"
                ),
                "p": 1e-4,
                "target_p_l": 3e-6,
                "compatibility_interval_p_l": [1.5e-6, 6e-6],
                "p_l": low_p,
                "p_l_wilson95": prime_row.as_dict()["logical_error_rate_wilson95"],
                "c_eff": c_eff,
                "evidence": "NUMERICAL",
            }
    summary_verdict = (
        prime_result["verdict"] if prime_result is not None else verdict
    )
    fit_verdict = "DIAGNOSTIC_ONLY" if label == "gate_a_prime" else verdict



    circuit_hash_material = "\n".join(circuit_text[p] for p in p_values).encode()
    circuit_hash = hashlib.sha256(circuit_hash_material).hexdigest()
    summary = {
        "schema": "msd.gate_a.campaign.v1",
        "evidence": "NUMERICAL",
        "run": {
            "script": "src/run_campaign.py",
            "requested_priority": "nice -n 10",
            "single_core": True,
            "shots_per_point": shots,
            "seed": seed,
            "stop_condition": f"completed exactly {len(p_values)} p points x {shots} shots",
            "label": label,

            "p_values": list(p_values),
            "wall_seconds_by_point": point_wall_seconds,
        },
        "paper_headline": {
            "reported_coefficient_c": 300.0,
            "pass_interval": [150.0, 600.0],
            "source": "arXiv:2605.21867 Fig. 10",
            "evidence": "REPORTED",
        },
        "quadratic_fit_verdict": fit_verdict,
        "gate_a_verdict": summary_verdict,
        "gate_a_prime": prime_result,
        "comparison_scope": "The verdict applies only to this Clifford/effective-output reconstruction; the paper does not name a decoder in its text.",
        "circuit": {
            "sha256": circuit_hash,
            "per_point": {
                _token(p): {
                    **circuit_metadata(p),
                    "sha256": hashlib.sha256(circuit_text[p].encode()).hexdigest(),
                    "stim_filename": f"circuit_p{_token(p)}.stim",
                }
                for p in p_values
            },
        },
        **result,
    }
    files = {p: f"circuit_p{_token(p)}.stim" for p in p_values}
    summary["artifact_files"] = ["summary.json", *files.values()]
    return summary, circuit_text


def _default_artifact_dir(summary: dict[str, object]) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_uuid = uuid.uuid4().hex[:8]
    digest = str(summary["circuit"]["sha256"])[:12]  # type: ignore[index]
    return ROOT / "campaigns" / f"{stamp}_{run_uuid}_{digest}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=100_000)
    parser.add_argument(
        "--p-values",
        default=",".join(str(p) for p in DEFAULT_P_VALUES),
        help="comma-separated physical error rates",
    )
    parser.add_argument(
        "--label", default="gate_a", choices=("gate_a", "gate_a_prime")
    )
    parser.add_argument("--seed", type=int, default=20260627)
    parser.add_argument("--artifact-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.shots <= 0:
        parser.error("--shots must be positive")
    try:
        p_values = tuple(
            float(part.strip()) for part in args.p_values.split(",") if part.strip()
        )
    except ValueError as exc:
        parser.error(f"invalid --p-values: {exc}")
    if any(not 0.0 <= p <= 1.0 for p in p_values):
        parser.error("all p values must be in [0, 1]")

    summary, circuit_text = run_campaign(
        p_values, args.shots, seed=args.seed, label=args.label
    )
    artifact_dir = args.artifact_dir or _default_artifact_dir(summary)
    artifact_dir.mkdir(parents=True, exist_ok=False)
    (artifact_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for p, text in circuit_text.items():
        (artifact_dir / f"circuit_p{_token(p)}.stim").write_text(text)
    print(json.dumps(summary, indent=2))
    print(f"artifact_dir={artifact_dir}")


if __name__ == "__main__":
    main()
