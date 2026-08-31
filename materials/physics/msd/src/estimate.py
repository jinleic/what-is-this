"""Sampling and fit helpers for the Gate-A Clifford surrogate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

import numpy as np

try:
    from .gate_a import build_gate_a_circuit
except ImportError:  # pragma: no cover - exercised by direct CLI invocation.
    from gate_a import build_gate_a_circuit


@dataclass(frozen=True)
class Estimate:
    physical_error_rate: float
    shots: int
    accepted: int
    rejected: int
    logical_failures: int
    acceptance_rate: float
    logical_error_rate: float | None
    evidence: str = "NUMERICAL"

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        if self.logical_error_rate is None:
            data["logical_error_rate_wilson95"] = None
            data["coefficient_c_eff"] = None
        else:
            data["logical_error_rate_wilson95"] = list(
                _wilson_interval(self.logical_failures, self.accepted)
            )
            data["coefficient_c_eff"] = self.logical_error_rate / (
                self.physical_error_rate * self.physical_error_rate
            )
        return data


def _wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0:
        return (math.nan, math.nan)
    phat = successes / trials
    denominator = 1.0 + z * z / trials
    center = (phat + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt(phat * (1.0 - phat) / trials + z * z / (4.0 * trials * trials)) / denominator
    return (max(0.0, center - radius), min(1.0, center + radius))


def sample_gate_a(
    physical_error_rate: float,
    shots: int,
    *,
    seed: int | None = 0,
) -> Estimate:
    """Sample accepted/failing shots from the detector-defined circuit.

    ``accepted`` means every verification/syndrome detector is zero.  A
    logical failure is an accepted shot with at least one of the three
    teleportation-frame observables set.  No decoder beyond that explicit
    detector postselection is assumed.
    """

    if shots <= 0:
        raise ValueError(f"shots must be positive, got {shots!r}")
    circuit = build_gate_a_circuit(float(physical_error_rate))
    sampler = circuit.compile_detector_sampler(seed=seed)
    detectors, observables = sampler.sample(shots, separate_observables=True)
    accepted_mask = ~np.any(detectors, axis=1)
    failure_mask = accepted_mask & np.any(observables, axis=1)
    accepted = int(np.count_nonzero(accepted_mask))
    failures = int(np.count_nonzero(failure_mask))
    rejected = int(shots - accepted)
    logical_error_rate = failures / accepted if accepted else None
    return Estimate(
        physical_error_rate=float(physical_error_rate),
        shots=int(shots),
        accepted=accepted,
        rejected=rejected,
        logical_failures=failures,
        acceptance_rate=accepted / shots,
        logical_error_rate=logical_error_rate,
    )


def fit_quadratic(estimates: Sequence[Estimate]) -> dict[str, object]:
    """Fit ``log(p_L) = slope*log(p) + intercept`` for positive estimates."""

    usable = [
        row
        for row in estimates
        if row.logical_error_rate is not None and row.logical_error_rate > 0.0
    ]
    if len(usable) < 2:
        return {
            "status": "INSUFFICIENT_POSITIVE_POINTS",
            "points_used": len(usable),
            "slope": None,
            "coefficient_c": None,
            "evidence": "NUMERICAL",
        }
    x = np.log([row.physical_error_rate for row in usable])
    y = np.log([row.logical_error_rate for row in usable])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    residual = y - predicted
    coefficient = float(math.exp(intercept))
    return {
        "status": "FIT",
        "points_used": len(usable),
        "slope": float(slope),
        "coefficient_c": coefficient,
        "rmse_log": float(math.sqrt(float(np.mean(residual * residual)))),
        "evidence": "NUMERICAL",
    }


def summarize(estimates: Iterable[Estimate]) -> dict[str, object]:
    rows = list(estimates)
    return {
        "estimates": [row.as_dict() for row in rows],
        "quadratic_fit": fit_quadratic(rows),
        "acceptance_ci_note": "Wilson 95% intervals can be recomputed from accepted/shots; not used in fit",
    }
