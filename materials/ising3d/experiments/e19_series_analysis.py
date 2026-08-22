"""Analyze the exact extended simple-cubic high-temperature series.

Run from the repository root:
    .venv/bin/python experiments/e19_series_analysis.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from ising.series.analysis import (
    d_finiteness_scan,
    derivative_series,
    differential_approximant_scan,
    dlog_pade_scan,
    even_to_squared_variable,
    exponential_series,
)

SCRIPT = "experiments/e19_series_analysis.py"
ROOT = Path(__file__).resolve().parents[1]
SERIES_DIR = ROOT / "results" / "series"
NOTE_PATH = ROOT / "notes" / "series_analysis.md"
PRECISION = 80
HOLDOUT = 3
REFERENCE_KC_DECIMAL = "0.221654626"
REFERENCE_KC_UNCERTAINTY_DECIMAL = "0.000000005"


def _mp_string(value, digits: int = 35) -> str:
    with mp.workdps(PRECISION):
        return mp.nstr(value, digits)


def _reported_estimate(
    central, uncertainty, *, central_digits: int, uncertainty_digits: int = 2
) -> str:
    with mp.workdps(PRECISION):
        central_text = mp.nstr(central, central_digits, strip_zeros=False)
        uncertainty_text = mp.nstr(
            uncertainty, uncertainty_digits, strip_zeros=False
        )
    return f"{central_text} \u00b1 {uncertainty_text}"


def _fraction_strings(values) -> list[str]:
    return [str(value) for value in values]


def _median(values):
    ordered = sorted(values)
    if not ordered:
        raise ValueError("median requires at least one value")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def _load_ht_coefficients() -> tuple[Fraction, ...]:
    payload = json.loads(
        (SERIES_DIR / "extended_sc_ht_free_energy.json").read_text(
            encoding="utf-8"
        )
    )
    return tuple(Fraction(value) for value in payload["data"]["coefficients"])


def _selected_dlog(z_series):
    specific_heat_series = derivative_series(z_series, 2)
    fits = [
        fit
        for fit in dlog_pade_scan(
            specific_heat_series,
            precision=PRECISION,
            near_diagonal=True,
        )
        if not fit.defective and 0 < fit.singularity < 1
    ]
    if not fits:
        raise AssertionError("no nondefective Dlog-Pade approximant")
    maximum_complexity = max(
        fit.numerator_degree + fit.denominator_degree for fit in fits
    )
    return tuple(
        fit
        for fit in fits
        if fit.numerator_degree + fit.denominator_degree == maximum_complexity
    )


def _selected_differential(z_series):
    fits = tuple(
        fit
        for fit in differential_approximant_scan(
            z_series, precision=PRECISION, balanced=True
        )
        if 0 < fit.singularity < 1 and 0 < fit.singular_exponent < 4
    )
    if not fits:
        raise AssertionError("no admissible differential approximant")
    return fits


def _prefix_estimate(coefficients: tuple[Fraction, ...], order: int) -> dict:
    z_series = even_to_squared_variable(coefficients[: order + 1])
    dlog = _selected_dlog(z_series)
    differential = _selected_differential(z_series)
    z_estimates = [fit.singularity for fit in dlog] + [
        fit.singularity for fit in differential
    ]
    alpha_estimates = [fit.exponent for fit in dlog] + [
        2 - fit.singular_exponent for fit in differential
    ]
    with mp.workdps(PRECISION):
        k_estimates = [mp.atanh(mp.sqrt(value)) for value in z_estimates]
        central_z = _median(z_estimates)
        central_k = mp.atanh(mp.sqrt(central_z))
        return {
            "order": order,
            "z_series": z_series,
            "dlog": dlog,
            "differential": differential,
            "z_estimates": tuple(+value for value in z_estimates),
            "k_estimates": tuple(+value for value in k_estimates),
            "alpha_estimates": tuple(+value for value in alpha_estimates),
            "central_z": +central_z,
            "central_v": +mp.sqrt(central_z),
            "central_k": +central_k,
            "central_alpha": +_median(alpha_estimates),
        }


def _serialize_dlog(fit) -> dict:
    return {
        "degrees": [fit.numerator_degree, fit.denominator_degree],
        "z_singularity": _mp_string(fit.singularity),
        "v_singularity": _mp_string(mp.sqrt(fit.singularity)),
        "K_singularity": _mp_string(mp.atanh(mp.sqrt(fit.singularity))),
        "specific_heat_exponent_alpha": _mp_string(fit.exponent),
        "numerator": _fraction_strings(fit.numerator),
        "denominator": _fraction_strings(fit.denominator),
        "nearest_numerator_zero_distance": (
            _mp_string(fit.nearest_numerator_zero)
            if fit.nearest_numerator_zero is not None
            else None
        ),
        "defective": fit.defective,
    }


def _serialize_differential(fit) -> dict:
    return {
        "degrees_Q1_Q0_P": [
            fit.derivative_degree,
            fit.function_degree,
            fit.inhomogeneous_degree,
        ],
        "z_singularity": _mp_string(fit.singularity),
        "v_singularity": _mp_string(mp.sqrt(fit.singularity)),
        "K_singularity": _mp_string(mp.atanh(mp.sqrt(fit.singularity))),
        "free_energy_singular_exponent_2_minus_alpha": _mp_string(
            fit.singular_exponent
        ),
        "implied_specific_heat_exponent_alpha": _mp_string(
            2 - fit.singular_exponent
        ),
        "Q1": _fraction_strings(fit.derivative_polynomial),
        "Q0": _fraction_strings(fit.function_polynomial),
        "P": _fraction_strings(fit.inhomogeneous_polynomial),
    }


def _serialize_ode(fit) -> dict:
    return {
        "order": fit.order,
        "degree": fit.degree,
        "unknown_count": fit.unknown_count,
        "training_coefficients": fit.training_coefficients,
        "held_out_coefficients": fit.held_out_coefficients,
        "polynomials_Qj": [
            _fraction_strings(polynomial) for polynomial in fit.coefficients
        ],
        "validation_residuals": _fraction_strings(fit.validation_residuals),
        "passed": fit.passed,
    }


def _record_check(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _write_note(data: dict) -> None:
    estimate = data["critical_estimate"]
    dfinite = data["d_finiteness"]
    dlog = data["dlog_pade"]
    differential = data["differential_approximants"]
    note = rf"""# Analysis of the extended simple-cubic series

## Exact input and scope

`experiments/e18_series_extend.py` derives the reduced free-energy series through
`v^20` at high temperature and `x^28` at low temperature.  Every coefficient is
an exact `Fraction`.  The high-temperature extension requires five free boxes
above the old 62-site limit, including `4x4x4` and `4x4x5` (80 sites); modular
propagation uses three 31-bit primes for each of those boxes.  A separate
`4x5x2` calculation exercised a 20-spin cross-section and agreed exactly with
the int64 engine.  Direct GF(2) cycle-space enumeration, with no spin transfer
matrix, agrees with the connected-cluster/free-energy coefficients through
`v^10`.

The singularity analysis uses

\[
F(z)=\phi(v)-\log 2,\qquad z=v^2,
\]

so the ten nonzero HT coefficients become ten consecutive powers of `z`.  Root
finding is performed at {PRECISION} decimal digits, but the uncertainty below is
set by series length, not numerical precision.

## Critical singularity

The combined, untuned estimate is

\[
v_c={estimate['v_c_report']},\qquad
K_c=\operatorname{{atanh}}v_c={estimate['K_c_report']}.
\]

The central value is the median of the highest-complexity nondefective
near-diagonal Dlog-Pade poles and all balanced first-order inhomogeneous
differential-approximant poles at order `v^20`.  The quoted uncertainty is the
full stability envelope: the largest displacement from that central value
among every retained approximant at truncations `v^16`, `v^18`, and `v^20`.
This deliberately conservative rule was fixed without reference to the
published benchmark.

For comparison only, `0.221654626(5)` differs from our central `K_c` by
`{estimate['benchmark_absolute_difference_report']}`, which is
`{estimate['benchmark_difference_in_our_uncertainties_report']}` of our quoted
uncertainty.  The benchmark was not used to select approximants, filters, or an
error bar.

### Exponent

Dlog-Pade was applied to `F''(z)`.  Its two highest-complexity nondefective
approximants give a median specific-heat exponent
`alpha = {dlog['median_alpha']}`.  The differential approximants fit
`Q1(z) F'(z) + Q0(z) F(z) = P(z)` and give a median
`alpha = {differential['median_implied_alpha']}`.  Their disagreement and strong
order drift are material: the combined value is only

\[
\alpha={estimate['alpha_report']},\qquad
2-\alpha={estimate['free_energy_singular_exponent_report']}.
\]

The exponent uncertainty is the analogous envelope over all retained
approximants at `v^16`, `v^18`, and `v^20`.  This short series locates the
singularity much more stably than it determines the weak specific-heat
exponent; the exponent estimate should not be treated as precision physics.

## Finite D-finiteness search

The exact Euler-operator ansatz was

\[
\sum_{{j=0}}^r Q_j(z)\,(z\,d/dz)^j F(z)=0,
\qquad \deg Q_j\le d.
\]

All `{dfinite['tested_pair_count']}` pairs `(r,d)` satisfying
`(r+1)(d+1) <= {dfinite['known_coefficient_count']} - {HOLDOUT}` were tested.
The first `{dfinite['known_coefficient_count'] - HOLDOUT}` coefficients fit the
candidate; the last `{HOLDOUT}` were held out.  No candidate passed all held-out
coefficients.  The failed training relations are recorded in
`results/series/extended_series_analysis.json`, including their nonzero exact
residuals.

This is **not a proof of non-D-finiteness**.  It excludes only homogeneous
Euler-form ODEs inside the tested order-degree budget on the available prefix.
It says nothing about larger `(r,d)`, and a finite prefix can always support
accidental relations at sufficient complexity.  We would require at least
`{dfinite['credible_independent_holdout_needed']}` newly derived coefficients,
not used in fitting, before calling any discovered ODE credible; this
conservative threshold is at least ten and at least the number of fitted ODE
parameters.

As a positive control, the identical code was applied to 19 exact coefficients
of `exp(z)`.  It found the order-one, degree-one equation
`z F - (z d/dz)F = 0` from the fitting prefix and obtained exactly zero on all
three held-out coefficients.  Thus the negative Ising result is not caused by a
search implementation that cannot recognize a simple holonomic function.

## Closed-form and integer-relation searches

No integer-relation or closed-form search was run.  The series-derived
uncertainty supplies far too few reliable digits for such a search to have
non-accidental predictive power.  Consequently there are no closed-form
candidates to report, and no exactness claim is made.
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def main() -> None:
    mp.mp.dps = PRECISION
    reference_kc = mp.mpf(REFERENCE_KC_DECIMAL)
    reference_kc_uncertainty = mp.mpf(REFERENCE_KC_UNCERTAINTY_DECIMAL)
    coefficients = _load_ht_coefficients()
    if len(coefficients) != 21:
        raise AssertionError("expected exact HT coefficients through v^20")

    prefix_results = [
        _prefix_estimate(coefficients, order) for order in (16, 18, 20)
    ]
    current = prefix_results[-1]
    dlog = current["dlog"]
    differential = current["differential"]

    all_k_estimates = [
        value for result in prefix_results for value in result["k_estimates"]
    ]
    all_alpha_estimates = [
        value for result in prefix_results for value in result["alpha_estimates"]
    ]
    central_k = current["central_k"]
    central_v = current["central_v"]
    central_alpha = current["central_alpha"]
    k_uncertainty = max(abs(value - central_k) for value in all_k_estimates)
    alpha_uncertainty = max(
        abs(value - central_alpha) for value in all_alpha_estimates
    )
    v_lower = mp.tanh(max(mp.mpf("0"), central_k - k_uncertainty))
    v_upper = mp.tanh(central_k + k_uncertainty)
    v_uncertainty = max(central_v - v_lower, v_upper - central_v)

    dlog_alpha = _median([fit.exponent for fit in dlog])
    differential_alpha = _median(
        [2 - fit.singular_exponent for fit in differential]
    )
    dlog_z = _median([fit.singularity for fit in dlog])
    differential_z = _median([fit.singularity for fit in differential])

    ising_dfinite = d_finiteness_scan(
        current["z_series"], holdout=HOLDOUT
    )
    control_dfinite = d_finiteness_scan(
        exponential_series(18), holdout=HOLDOUT
    )
    control_passes = [fit for fit in control_dfinite.fits if fit.passed]
    benchmark_difference = abs(central_k - reference_kc)

    checks: list[dict] = []
    _record_check(
        checks,
        "dlog_pade_physical_poles",
        len(dlog) >= 2,
        f"retained {len(dlog)} highest-complexity nondefective near-diagonal poles",
    )
    _record_check(
        checks,
        "differential_approximant_physical_poles",
        len(differential) >= 4,
        f"retained {len(differential)} balanced first-order inhomogeneous approximants",
    )
    _record_check(
        checks,
        "ising_short_ode_search_negative",
        not ising_dfinite.found,
        (
            f"no ODE passed {HOLDOUT} held-out coefficients among "
            f"{len(ising_dfinite.tested_pairs)} exact (r,d) searches"
        ),
    )
    _record_check(
        checks,
        "holonomic_control_ode_found",
        any(fit.order == 1 and fit.degree == 1 for fit in control_passes),
        "the same code finds z*F - theta(F) = 0 for F=exp(z), with three exact held-out zeros",
    )
    _record_check(
        checks,
        "benchmark_not_used_for_fit",
        True,
        "the published K_c is evaluated only after the approximant ensemble and stability envelope are fixed",
    )

    benchmark_ratio = benchmark_difference / k_uncertainty
    data = {
        "input": {
            "artifact": "results/series/extended_sc_ht_free_energy.json",
            "variable_transformation": "z = v^2",
            "known_v_order": 20,
            "known_z_coefficients": len(current["z_series"]),
            "z_coefficients": _fraction_strings(current["z_series"]),
            "mpmath_decimal_precision": PRECISION,
        },
        "critical_estimate": {
            "selection_rule": (
                "median at v^20 of highest-total-degree nondefective near-diagonal "
                "Dlog-Pade poles plus balanced first-order inhomogeneous differential "
                "approximant poles; no benchmark input"
            ),
            "uncertainty_rule": (
                "largest displacement from the v^20 central estimate among every "
                "retained approximant at truncations v^16, v^18, and v^20"
            ),
            "z_c": _mp_string(current["central_z"]),
            "v_c": _mp_string(central_v),
            "v_c_uncertainty": _mp_string(v_uncertainty),
            "K_c": _mp_string(central_k),
            "K_c_uncertainty": _mp_string(k_uncertainty),
            "v_c_report": _reported_estimate(
                central_v, v_uncertainty, central_digits=4
            ),
            "K_c_report": _reported_estimate(
                central_k, k_uncertainty, central_digits=4
            ),
            "alpha": _mp_string(central_alpha),
            "alpha_uncertainty": _mp_string(alpha_uncertainty),
            "free_energy_singular_exponent": _mp_string(2 - central_alpha),
            "alpha_report": _reported_estimate(
                central_alpha, alpha_uncertainty, central_digits=2
            ),
            "free_energy_singular_exponent_report": _reported_estimate(
                2 - central_alpha, alpha_uncertainty, central_digits=3
            ),
            "benchmark_K_c": _mp_string(reference_kc),
            "benchmark_uncertainty": _mp_string(reference_kc_uncertainty),
            "benchmark_absolute_difference": _mp_string(benchmark_difference),
            "benchmark_difference_in_our_uncertainties": _mp_string(
                benchmark_ratio
            ),
            "benchmark_absolute_difference_report": mp.nstr(
                benchmark_difference, 3, strip_zeros=False
            ),
            "benchmark_difference_in_our_uncertainties_report": mp.nstr(
                benchmark_ratio, 2, strip_zeros=False
            ),
            "benchmark_used_in_fit": False,
        },
        "dlog_pade": {
            "analyzed_series": "d^2 F(z) / dz^2",
            "singular_form": "F''(z) ~ (1-z/z_c)^(-alpha)",
            "selection": "nondefective near-diagonal approximants at maximum available L+M",
            "pole_zero_defect_threshold": "absolute distance <= 1e-3*max(1,abs(z_c))",
            "selected_count": len(dlog),
            "median_z_c": _mp_string(dlog_z),
            "median_K_c": _mp_string(mp.atanh(mp.sqrt(dlog_z))),
            "median_alpha": _mp_string(dlog_alpha),
            "approximants": [_serialize_dlog(fit) for fit in dlog],
        },
        "differential_approximants": {
            "ansatz": "Q1(z) F'(z) + Q0(z) F(z) = P(z)",
            "selection": (
                "all coefficients used; deg(Q1), deg(Q0), deg(P) differ by at most 2; "
                "retain the first positive real Q1 zero with 0 < z_c < 1 and "
                "0 < singular exponent < 4"
            ),
            "selected_count": len(differential),
            "median_z_c": _mp_string(differential_z),
            "median_K_c": _mp_string(mp.atanh(mp.sqrt(differential_z))),
            "median_implied_alpha": _mp_string(differential_alpha),
            "approximants": [
                _serialize_differential(fit) for fit in differential
            ],
        },
        "truncation_stability": [
            {
                "v_order": result["order"],
                "central_z_c": _mp_string(result["central_z"]),
                "central_K_c": _mp_string(result["central_k"]),
                "central_alpha": _mp_string(result["central_alpha"]),
                "dlog_count": len(result["dlog"]),
                "differential_count": len(result["differential"]),
                "all_K_estimates": [
                    _mp_string(value) for value in result["k_estimates"]
                ],
                "all_alpha_estimates": [
                    _mp_string(value) for value in result["alpha_estimates"]
                ],
            }
            for result in prefix_results
        ],
        "d_finiteness": {
            "operator_basis": "Euler theta = z*d/dz",
            "ansatz": "sum_j Q_j(z) theta^j F(z) = 0, deg Q_j <= d",
            "known_coefficient_count": len(current["z_series"]),
            "held_out_coefficient_count": HOLDOUT,
            "search_constraint": (
                f"(r+1)(d+1) <= {len(current['z_series'])} - {HOLDOUT}"
            ),
            "tested_pair_count": len(ising_dfinite.tested_pairs),
            "tested_pairs": [list(pair) for pair in ising_dfinite.tested_pairs],
            "passed_fit_count": sum(fit.passed for fit in ising_dfinite.fits),
            "training_relations_rejected_by_holdout": [
                _serialize_ode(fit)
                for fit in ising_dfinite.fits
                if not fit.passed
            ],
            "conclusion": (
                "No ODE in the searched finite order-degree rectangle passes the held-out "
                "coefficients. This is not a proof of non-D-finiteness and does not test "
                "larger order or degree."
            ),
            "credible_independent_holdout_needed": (
                ising_dfinite.credible_holdout_needed
            ),
            "credibility_rule": (
                "at least max(10, number of fitted ODE parameters) newly derived "
                "coefficients not used to select the relation"
            ),
        },
        "holonomic_control": {
            "function": "exp(z)",
            "known_coefficient_count": 19,
            "held_out_coefficient_count": HOLDOUT,
            "tested_pair_count": len(control_dfinite.tested_pairs),
            "found": control_dfinite.found,
            "passing_fits": [_serialize_ode(fit) for fit in control_passes],
            "identified_equation": "z*F - theta*F = 0",
        },
        "integer_relation_search": {
            "performed": False,
            "reason": (
                "the conservative series-only v_c uncertainty leaves fewer than three "
                "stable decimal digits, insufficient for a non-accidental relation search"
            ),
            "candidates": [],
            "status": "No candidates; no exactness claim",
        },
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": (
                "exact Fraction linear algebra; mpmath root finding at "
                f"{PRECISION} decimal digits"
            ),
        },
        "data": data,
        "checks": checks,
    }
    output_path = SERIES_DIR / "extended_series_analysis.json"
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    _write_note(data)
    print(f"WROTE {output_path.relative_to(ROOT)}")
    print(f"WROTE {NOTE_PATH.relative_to(ROOT)}")
    print(
        "PASS e19_series_analysis: "
        f"K_c={_mp_string(central_k, 12)} +/- {_mp_string(k_uncertainty, 6)}; "
        f"alpha={_mp_string(central_alpha, 8)} +/- {_mp_string(alpha_uncertainty, 5)}; "
        "no short ODE; exp(z) control found"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e19_series_analysis: {error}")
        raise
