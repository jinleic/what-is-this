"""Compute and certify elementary rigorous bounds on the simple-cubic Ising Kc.

Run from the repository root:
    .venv/bin/python experiments/e05_kc_bounds.py
"""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
import json
from datetime import datetime, timezone
import os
from pathlib import Path
import time

import mpmath as mp

from ising.rigorous_bounds import (
    OEIS_A001412,
    OEIS_A001412_URL,
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


SCRIPT = "experiments/e05_kc_bounds.py"
RESULT_PATH = Path(__file__).resolve().parents[1] / "results" / "bounds" / "kc_bounds.json"
PRECISION = 80
mp.mp.dps = PRECISION
MAX_SAW_STEPS = 14
WORKERS = min(8, os.cpu_count() or 1)
BENCHMARK_KC = mp.mpf("0.221654626")
CERTIFIED_DECIMAL_PLACES = 40


def _record_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _mp_text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def _outward_decimal(value: str, places: int, *, upper: bool) -> str:
    with localcontext() as context:
        context.prec = max(PRECISION + 20, places + 20)
        quantum = Decimal(1).scaleb(-places)
        rounding = ROUND_CEILING if upper else ROUND_FLOOR
        return format(Decimal(value).quantize(quantum, rounding=rounding), "f")


def _submultiplicative_through(counts: tuple[int, ...]) -> bool:
    return all(
        counts[left + right] <= counts[left] * counts[right]
        for left in range(len(counts))
        for right in range(len(counts) - left)
    )


def main() -> None:
    mp.mp.dps = PRECISION
    checks: list[dict[str, object]] = []

    start = time.perf_counter()
    counts = enumerate_saw_counts(MAX_SAW_STEPS, workers=WORKERS)
    enumeration_seconds = time.perf_counter() - start
    oeis_counts = fetch_oeis_a001412(MAX_SAW_STEPS)
    _record_check(
        checks,
        "saw_counts_match_live_oeis_a001412",
        counts == oeis_counts == OEIS_A001412[: MAX_SAW_STEPS + 1],
        f"exact backtracking c_0..c_{MAX_SAW_STEPS} agree with {OEIS_A001412_URL}",
    )
    _record_check(
        checks,
        "finite_submultiplicativity_check",
        _submultiplicative_through(counts),
        "c_(m+n) <= c_m*c_n for every available pair m+n <= 14",
    )

    count_n = counts[MAX_SAW_STEPS]
    mu_upper = connective_constant_upper(count_n, MAX_SAW_STEPS, PRECISION)
    mu_interval = connective_constant_upper_interval(
        count_n, MAX_SAW_STEPS, PRECISION
    )
    lower = saw_kc_lower(count_n, MAX_SAW_STEPS, PRECISION)
    lower_interval = saw_kc_lower_interval(count_n, MAX_SAW_STEPS, PRECISION)
    certified_lower_decimal = _outward_decimal(
        lower_interval[0], CERTIFIED_DECIMAL_PLACES, upper=False
    )
    _record_check(
        checks,
        "saw_lower_bound_benchmark_side",
        lower < BENCHMARK_KC,
        f"{_mp_text(lower, 30)} < 0.221654626 (benchmark used only for falsification)",
    )

    square_upper = square_lattice_kc(PRECISION)
    square_interval = square_lattice_kc_interval(PRECISION)
    certified_square_decimal = _outward_decimal(
        square_interval[1], CERTIFIED_DECIMAL_PLACES, upper=True
    )
    _record_check(
        checks,
        "square_layer_upper_bound_benchmark_side",
        square_upper > BENCHMARK_KC,
        f"{_mp_text(square_upper, 30)} > 0.221654626",
    )

    i3_numeric = watson_integral_numeric(3, PRECISION)
    w3_numeric = 3 * i3_numeric
    w3_corrected = watson_closed_form(PRECISION, denominator=32)
    w3_formula_as_stated = watson_closed_form(PRECISION, denominator=4)
    i3_interval = watson_i3_interval(PRECISION)
    i3_lo, i3_hi = (mp.mpf(endpoint) for endpoint in i3_interval)
    _record_check(
        checks,
        "watson_numeric_vs_corrected_closed_form",
        abs(w3_numeric - w3_corrected) < mp.mpf("1e-60"),
        "Bessel-Laplace quadrature agrees with the Glasser-Zucker 1/(32*pi^3) gamma formula beyond 60 digits",
    )
    _record_check(
        checks,
        "watson_interval_encloses_closed_form",
        i3_lo <= w3_corrected / 3 <= i3_hi
        and i3_hi - i3_lo < mp.mpf("1e-70"),
        f"directed-rounding interval width < 1e-70 at {PRECISION} decimal digits",
    )
    _record_check(
        checks,
        "problem_stated_watson_formula_factor_error",
        abs(w3_formula_as_stated / w3_corrected - 8) < mp.mpf("1e-70"),
        "the supplied 1/(4*pi^3) expression is exactly 8 times W_sc; the literature formula has 1/(32*pi^3)",
    )

    infrared_upper = infrared_kc_upper(3, PRECISION)
    infrared_interval = infrared_kc_upper_interval(PRECISION)
    certified_infrared_decimal = _outward_decimal(
        infrared_interval[1], CERTIFIED_DECIMAL_PLACES, upper=True
    )
    _record_check(
        checks,
        "infrared_upper_bound_benchmark_side",
        infrared_upper > BENCHMARK_KC,
        f"I_3/2 = {_mp_text(infrared_upper, 30)} > 0.221654626; no factor-of-two normalization failure",
    )
    _record_check(
        checks,
        "infrared_bound_improves_layer_bound",
        infrared_upper < square_upper,
        f"{_mp_text(infrared_upper, 20)} < {_mp_text(square_upper, 20)}",
    )

    i2 = watson_integral_numeric(2, PRECISION)
    i4 = watson_integral_numeric(4, PRECISION)
    i5 = watson_integral_numeric(5, PRECISION)
    _record_check(
        checks,
        "infrared_dimension_cross_checks",
        mp.isinf(i2) and i3_numeric > i4 > i5 > 0,
        "I_2 diverges (bound vacuous), while I_3 > I_4 > I_5 > 0",
    )

    _record_check(
        checks,
        "final_interval_benchmark_sanity",
        mp.mpf(certified_lower_decimal)
        < BENCHMARK_KC
        < mp.mpf(certified_infrared_decimal),
        f"benchmark 0.221654626 lies inside [{certified_lower_decimal}, {certified_infrared_decimal}]",
    )

    data = {
        "model": "nearest-neighbor ferromagnetic Ising model on the simple cubic lattice",
        "coupling_convention": "K = beta*J and v = tanh(K)",
        "arithmetic": {
            "combinatorial": "exact Python integers",
            "mpmath_version": mp.__version__,
            "transcendental_point_values": f"mpmath with mp.dps={PRECISION}",
            "certified_intervals": f"mpmath.iv directed rounding with iv.dps={PRECISION}",
            "reported_certified_decimal_places": CERTIFIED_DECIMAL_PLACES,
        },
        "lower_bound_saw": {
            "theorem": "<sigma_0 sigma_x> <= sum_{self-avoiding 0-to-x paths gamma} v^|gamma|",
            "logical_chain": [
                "c_(m+n) <= c_m*c_n",
                "mu = inf_n c_n^(1/n) <= c_14^(1/14)",
                "v_c >= 1/mu >= c_14^(-1/14)",
                "K_c >= atanh(c_14^(-1/14))",
            ],
            "enumeration": {
                "method": "visited-set backtracking, quotiented by all 48 signed coordinate permutations at depth 6",
                "workers": WORKERS,
                "elapsed_seconds_numerical_metadata": enumeration_seconds,
                "max_steps": MAX_SAW_STEPS,
                "counts": {str(index): value for index, value in enumerate(counts)},
                "oeis_cross_check_url": OEIS_A001412_URL,
                "oeis_live_match": True,
            },
            "mu_upper_exact": f"{count_n}^(1/{MAX_SAW_STEPS})",
            "mu_upper_interval": list(mu_interval),
            "mu_upper_point_numerical": _mp_text(mu_upper),
            "kc_lower_exact": f"atanh({count_n}^(-1/{MAX_SAW_STEPS}))",
            "kc_lower_interval": list(lower_interval),
            "kc_lower_point_numerical": _mp_text(lower),
            "certified_decimal_lower": certified_lower_decimal,
        },
        "upper_bound_square_layers": {
            "theorem": "GKS-II monotonicity under addition of nonnegative vertical bonds",
            "logical_chain": [
                "decoupled square layers order for K > Kc_square",
                "adding ferromagnetic vertical bonds cannot decrease correlations or plus magnetization",
                "Kc_simple_cubic <= Kc_square",
            ],
            "exact": "log(1+sqrt(2))/2",
            "interval": list(square_interval),
            "point_numerical": _mp_text(square_upper),
            "certified_decimal_upper": certified_square_decimal,
        },
        "upper_bound_infrared": {
            "normalization": {
                "ordered_pair_interaction": "j(+/-e_i)=J/2, so -sum_(x,y) j(y-x)sigma_x sigma_y = -J sum_<x,y> sigma_x sigma_y",
                "dispersion": "E(k)=sum_z j(z)*(1-cos(k.z))=J*(3-cos(k1)-cos(k2)-cos(k3))",
                "infrared_bound": "Ghat(k) <= 1/(2*beta*E(k)) = 1/(2*K*(3-sum_i cos(ki))) for k != 0",
                "sum_rule": "1 = m_LRO^2 + integral Ghat_regular",
                "conclusion": "Kc_simple_cubic <= I_3/2",
            },
            "watson": {
                "I3_definition": "(2*pi)^(-3) integral_[-pi,pi]^3 dk/(3-cos(k1)-cos(k2)-cos(k3))",
                "Wsc_definition": "3*I3",
                "I3_bessel_quadrature_numerical": _mp_text(i3_numeric),
                "I3_certified_interval": list(i3_interval),
                "Wsc_bessel_quadrature_numerical": _mp_text(w3_numeric),
                "Wsc_corrected_gamma_formula": "sqrt(6)*Gamma(1/24)*Gamma(5/24)*Gamma(7/24)*Gamma(11/24)/(32*pi^3)",
                "Wsc_corrected_gamma_numerical": _mp_text(w3_corrected),
                "formula_in_task": "sqrt(6)*Gamma(1/24)*Gamma(5/24)*Gamma(7/24)*Gamma(11/24)/(4*pi^3)",
                "formula_in_task_numerical": _mp_text(w3_formula_as_stated),
                "formula_in_task_ratio_to_Wsc": _mp_text(
                    w3_formula_as_stated / w3_corrected
                ),
                "discrepancy": "the task-stated formula is too large by exactly a factor of 8",
            },
            "kc_upper_exact": "I_3/2 = W_sc/6",
            "kc_upper_interval": list(infrared_interval),
            "kc_upper_point_numerical": _mp_text(infrared_upper),
            "certified_decimal_upper": certified_infrared_decimal,
            "dimension_cross_check": {
                "d2": "diverges; infrared magnetization criterion is vacuous",
                "d3_I_d": _mp_text(i3_numeric),
                "d3_I_d_over_2": _mp_text(i3_numeric / 2),
                "d4_I_d": _mp_text(i4),
                "d4_I_d_over_2": _mp_text(i4 / 2),
                "d5_I_d": _mp_text(i5),
                "d5_I_d_over_2": _mp_text(i5 / 2),
            },
        },
        "final_certified_interval": {
            "exact": f"[atanh({count_n}^(-1/{MAX_SAW_STEPS})), W_sc/6]",
            "decimal_outward_rounded": [
                certified_lower_decimal,
                certified_infrared_decimal,
            ],
            "benchmark_falsification_target": "0.221654626(5)",
            "benchmark_role": "comparison only; no constant was fitted to it",
        },
        "rigor_classification": {
            "rigorous_theorems": [
                "high-temperature SAW domination of ferromagnetic Ising correlations",
                "submultiplicativity of rooted SAW counts and Fekete's lemma",
                "GKS-II monotonicity under addition of ferromagnetic bonds",
                "Onsager square-lattice critical coupling",
                "Frohlich-Simon-Spencer Gaussian domination / infrared bound",
                "Glasser-Zucker corrected gamma evaluation of the simple-cubic Watson integral",
            ],
            "computer_assisted_rigorous": [
                "exact integer backtracking counts c_0 through c_14",
                "directed-rounding enclosures for the algebraic/transcendental bound constants",
            ],
            "numerical_only": [
                "Bessel-Laplace quadrature point values",
                "d=4 and d=5 trend values",
                "wall-clock enumeration time",
                "comparison with the published Kc benchmark",
            ],
        },
    }

    payload = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": PRECISION,
        },
        "data": data,
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(RESULT_PATH.parents[2])}")
    print(
        "PASS e05_kc_bounds: certified "
        f"{certified_lower_decimal} <= Kc <= {certified_infrared_decimal}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e05_kc_bounds: {error}")
        raise
