"""Sampling driver for the distance-carrying output patch (msd/src/expanded.py).

Reuses ``estimate.Estimate`` rows and ``fit_quadratic`` so results are
directly comparable with the frozen Gate-A campaign.  No directories besides
``--scratch-root`` (default ``msd/scratch``) are written; campaigns/ is
explicitly forbidden by task constraints.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    from .estimate import Estimate, _wilson_interval, fit_quadratic
    from .expanded import build_expanded_circuit
except ImportError:  # direct script invocation
    from estimate import Estimate, _wilson_interval, fit_quadratic
    from expanded import build_expanded_circuit


def _token(p: float) -> str:
    return f"{p:.0e}".replace("+", "p").replace("-", "m")



def _single_call_counts(
    circuit: "stim.Circuit", shots: int, seed: int
) -> tuple[int, int, np.ndarray]:
    """One big ``sample`` call per point.

    Chunked sampling through one continued detector sampler was verified NOT
    bit-equal to a single call at equal seed (p=1e-3, 2000 shots, 4x500,
    verified 2026-08-30), so per-point counts come from a single sample.
    Returns (rejected, accepted_failures_total, per_observable_counts).
    """

    sampler = circuit.compile_detector_sampler(seed=seed)
    detectors, observables = sampler.sample(shots, separate_observables=True)
    accepted_mask = ~np.any(detectors, axis=1)
    rejected = int(np.count_nonzero(np.any(detectors, axis=1)))
    failure_mask = accepted_mask & np.any(observables, axis=1)
    failures = int(np.count_nonzero(failure_mask))
    per_obs = np.array(
        [
            int(np.count_nonzero(accepted_mask & observables[:, k]))
            for k in range(observables.shape[1])
        ],
        dtype=np.int64,
    )
    return rejected, failures, per_obs



def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--p", dest="p_values", action="append", type=float, required=True,
        help="physical error rate; repeatable",
    )
    parser.add_argument("--shots", type=int, required=True)
    parser.add_argument("--seed", type=int, default=20260627)
    parser.add_argument("--distance", type=int, default=3)
    parser.add_argument(
        "--scratch-root", type=Path,
        default=Path(__file__).resolve().parents[1] / "scratch",
    )
    args = parser.parse_args()

    p_values: list[float] = list(dict.fromkeys(args.p_values))
    if not p_values:
        raise SystemExit("at least one --p is required")
    if args.shots <= 0:
        raise SystemExit("--shots must be positive")
    if args.distance < 3 or args.distance % 2 == 0:
        raise SystemExit("--distance must be an odd integer >= 3")

    scratch_root = args.scratch_root
    scratch_root.mkdir(parents=True, exist_ok=True)

    estimates: list[Estimate] = []
    per_point: list[dict[str, object]] = []
    for index, p in enumerate(p_values):
        per_point_seed = args.seed + index
        circuit = build_expanded_circuit(p, output_distance=args.distance)
        det_rejected, failures, per_obs = _single_call_counts(
            circuit, args.shots, per_point_seed
        )
        accepted = args.shots - det_rejected
        if accepted <= 0:
            raise SystemExit(f"no accepted shots at p={p}; cannot estimate")
        per_point.append(
            {
                "p": p,
                "seed": per_point_seed,
                "shots": args.shots,
                "accepted": accepted,
                "rejected": det_rejected,
                "accepted_failures_total": failures,
                "per_observable_failures": {
                    str(k): int(per_obs[k]) for k in range(len(per_obs))
                },
                "per_observable_wilson95": {
                    str(k): list(
                        _wilson_interval(int(per_obs[k]), accepted)
                    ) for k in range(len(per_obs))
                },
                "logical_error_rate": failures / accepted,
            }
        )
        estimates.append(
            Estimate(
                physical_error_rate=p,
                shots=args.shots,
                accepted=accepted,
                rejected=det_rejected,
                logical_failures=failures,
                acceptance_rate=accepted / args.shots,
                logical_error_rate=failures / accepted,
            )
        )

    fit = fit_quadratic(estimates)
    summary = {
        "evidence": "NUMERICAL",
        "module": "msd/src/expanded.py",
        "distance": args.distance,
        "shots_per_point": args.shots,
        "seed_base": args.seed,
        "p_values": p_values,
        "per_point": per_point,
        "estimates": [row.as_dict() for row in estimates],
        "quadratic_fit": fit,
    }
    token = "s".join(_token(p) for p in p_values)
    out_name = (
        f"expanded_d{args.distance}_{token}_shots{args.shots}"
        f"_seed{args.seed}.json"
    )
    out_path = scratch_root / out_name
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out_path}")
    for row, point in zip(estimates, per_point):
        print(
            f"p={row.physical_error_rate:g} accepted={point['accepted']}"
            f" failures={point['accepted_failures_total']}"
            f" per_obs={point['per_observable_failures']}"
            f" pL={row.logical_error_rate:.4g}"
            f" c_eff={point['logical_error_rate']/ (p*p):.4g}"
        )
    print(f"quadratic_fit: {fit}")


if __name__ == "__main__":
    main()
