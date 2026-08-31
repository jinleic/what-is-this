"""Bootstrap confidence intervals for logical error rates.

Standard nonparametric bootstrap over shots (resample errors with
replacement via binomial resampling of the failure count) and a
Wilson-score interval helper. 95% CI default (pre_statement.md).
"""
from __future__ import annotations

import math

import numpy as np


def wilson_interval(failures: int, shots: int, z: float = 1.959963985) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if shots == 0:
        return (0.0, 1.0)
    ph = failures / shots
    denom = 1 + z * z / shots
    center = (ph + z * z / (2 * shots)) / denom
    half = z * math.sqrt(ph * (1 - ph) / shots + z * z / (4 * shots * shots)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def bootstrap_ci_binomial(
    failures: int,
    shots: int,
    n_boot: int = 2000,
    ci: float = 0.95,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Percentile bootstrap CI for a binomial proportion.

    Deterministic given rng; use a seed recorded in the campaign ledger.
    """
    rng = rng or np.random.default_rng(0)
    phat = failures / shots if shots else 0.0
    draws = rng.binomial(shots, phat, size=n_boot) / shots
    lo_q = (1 - ci) / 2 * 100
    hi_q = (1 + ci) / 2 * 100
    return tuple(np.percentile(draws, [lo_q, hi_q]))  # type: ignore[return-value]


def _percentile_with_infinity(values: np.ndarray, q: float) -> float:
    """NumPy's linear percentile, with finite-to-infinity interpolation = infinity."""
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    if ordered.size == 0:
        raise ValueError("percentile requires at least one value")
    position = (ordered.size - 1) * q / 100.0
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    hi = float(ordered[upper])
    if math.isinf(hi):
        return hi
    lo = float(ordered[lower])
    return lo + (position - lower) * (hi - lo)


def ratio_ci(
    num_failures: int,
    num_shots: int,
    den_failures: int,
    den_shots: int,
    n_boot: int = 2000,
    ci: float = 0.95,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Parametric-bootstrap CI for a ratio of independent binomial rates.

    Positive-numerator draws with a zero denominator are retained as
    ``+inf``; silently dropping them truncates the upper tail exactly in the
    sparse-failure regime. Draws with both counts zero have no defined ratio
    and are excluded.
    """
    rng = rng or np.random.default_rng(0)
    a = num_failures / num_shots if num_shots else 0.0
    b = den_failures / den_shots if den_shots else 0.0
    num_draws = rng.binomial(num_shots, a, size=n_boot)
    den_draws = rng.binomial(den_shots, b, size=n_boot)
    defined = (num_draws != 0) | (den_draws != 0)
    num_draws = num_draws[defined]
    den_draws = den_draws[defined]
    if num_draws.size == 0:
        return (0.0, float("inf"))

    ratios = np.full(num_draws.shape, float("inf"), dtype=np.float64)
    nonzero_den = den_draws != 0
    ratios[nonzero_den] = (
        num_draws[nonzero_den] * den_shots
        / (den_draws[nonzero_den] * max(num_shots, 1))
    )
    lo_q = (1 - ci) / 2 * 100
    hi_q = (1 + ci) / 2 * 100
    return (
        _percentile_with_infinity(ratios, lo_q),
        _percentile_with_infinity(ratios, hi_q),
    )


if __name__ == "__main__":
    print(wilson_interval(5, 1000))
    print(bootstrap_ci_binomial(5, 1000, n_boot=500, rng=np.random.default_rng(7)))
    print(ratio_ci(5, 1000, 50, 1000, n_boot=500, rng=np.random.default_rng(7)))
