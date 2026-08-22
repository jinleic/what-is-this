"""Standalone acceptance checks for the certified critical-coupling bounds."""

from __future__ import annotations

import os

import mpmath as mp

from ising.rigorous_bounds import (
    OEIS_A001412,
    connective_constant_upper,
    connective_constant_upper_interval,
    enumerate_saw_counts,
    fetch_oeis_a001412,
    infrared_kc_upper,
    infrared_kc_upper_interval,
    saw_kc_lower,
    saw_kc_lower_interval,
    square_lattice_kc,
    square_lattice_kc_interval,
    watson_closed_form,
    watson_i3_interval,
    watson_integral_numeric,
)


PRECISION = 80
mp.mp.dps = PRECISION
BENCHMARK_KC = mp.mpf("0.221654626")
MAX_SAW_STEPS = 14


def main() -> None:
    workers = min(8, os.cpu_count() or 1)
    counts = enumerate_saw_counts(MAX_SAW_STEPS, workers=workers)
    expected = OEIS_A001412[: MAX_SAW_STEPS + 1]
    live_expected = fetch_oeis_a001412(MAX_SAW_STEPS)
    assert counts == expected == live_expected, (
        "backtracking SAW counts disagree with live OEIS A001412 b-file",
        counts,
        live_expected,
    )
    print(f"SAW counts c_0..c_{MAX_SAW_STEPS} match live OEIS A001412: PASS")

    mu_upper = connective_constant_upper(counts[-1], MAX_SAW_STEPS, PRECISION)
    lower = saw_kc_lower(counts[-1], MAX_SAW_STEPS, PRECISION)
    square_upper = square_lattice_kc(PRECISION)
    infrared_upper = infrared_kc_upper(3, PRECISION)
    mu_lo, mu_hi = (
        mp.mpf(x)
        for x in connective_constant_upper_interval(
            counts[-1], MAX_SAW_STEPS, PRECISION
        )
    )
    lower_lo, lower_hi = (
        mp.mpf(x)
        for x in saw_kc_lower_interval(counts[-1], MAX_SAW_STEPS, PRECISION)
    )
    square_lo, square_hi = (
        mp.mpf(x) for x in square_lattice_kc_interval(PRECISION)
    )
    infrared_lo, infrared_hi = (
        mp.mpf(x) for x in infrared_kc_upper_interval(PRECISION)
    )
    assert mu_lo <= mu_upper <= mu_hi
    assert lower_lo <= lower <= lower_hi
    assert square_lo <= square_upper <= square_hi
    assert infrared_lo <= infrared_upper <= infrared_hi
    assert lower < BENCHMARK_KC, (lower, BENCHMARK_KC)
    assert square_upper > BENCHMARK_KC, (square_upper, BENCHMARK_KC)
    assert infrared_upper > BENCHMARK_KC, (infrared_upper, BENCHMARK_KC)
    assert infrared_upper < square_upper
    print(
        "bound-side sanity checks: PASS",
        f"({mp.nstr(lower, 18)} < Kc < {mp.nstr(infrared_upper, 18)})",
    )

    integral = watson_integral_numeric(3, PRECISION)
    corrected_closed_form = watson_closed_form(PRECISION, denominator=32) / 3
    stated_closed_form = watson_closed_form(PRECISION, denominator=4) / 3
    with mp.workdps(15):
        low_context_upper = infrared_kc_upper(3, PRECISION)
    assert abs(low_context_upper - corrected_closed_form / 2) < mp.mpf("1e-70")
    assert abs(integral - corrected_closed_form) < mp.mpf("1e-60")
    assert abs(stated_closed_form / integral - 8) < mp.mpf("1e-60")
    lo, hi = (mp.mpf(x) for x in watson_i3_interval(PRECISION))
    assert lo <= corrected_closed_form <= hi, (lo, corrected_closed_form, hi)
    assert hi - lo < mp.mpf("1e-70")
    print("Watson integral and corrected gamma closed form agree beyond 25 digits: PASS")
    print("problem-stated 1/(4*pi^3) gamma formula rejected (factor 8): PASS")

    assert mp.isinf(watson_integral_numeric(2, PRECISION))
    upper4 = infrared_kc_upper(4, PRECISION)
    upper5 = infrared_kc_upper(5, PRECISION)
    assert infrared_upper > upper4 > upper5 > 0
    print("d=2 divergence and d=3,4,5 normalization trend: PASS")

    assert mu_upper > 1
    print("ALL rigorous-bounds checks: PASS")


if __name__ == "__main__":
    main()
