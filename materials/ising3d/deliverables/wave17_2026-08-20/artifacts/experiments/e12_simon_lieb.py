"""Compute certified Simon--Lieb finite-box lower bounds on simple-cubic Kc.

Run from the repository root:
    .venv/bin/python experiments/e12_simon_lieb.py
"""

from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp

from ising.rigorous_bounds.simon_lieb import (
    atanh_rational_interval,
    boundary_multiplicities,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
    kappa_float_from_v,
    kappa_mpf_from_v,
    mpmath_relative_roundoff_bound,
    solve_box_root_float,
)


SCRIPT = "experiments/e12_simon_lieb.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "simon_lieb_bounds.json"
SAW_RESULT_PATH = ROOT / "results" / "bounds" / "kc_bounds.json"
BOXES = (
    (1, 1, 1),
    (2, 2, 2),
    (3, 3, 3),
    (4, 4, 4),
    (4, 4, 8),
    (4, 4, 16),
    (4, 4, 32),
)
EXACT_RESIDUAL_BOXES = frozenset(BOXES[:4])
V_DECIMAL_PLACES = 14
K_DECIMAL_PLACES = 14
LOW_DPS = 50
HIGH_DPS = 80
BENCHMARK_KC = mp.mpf("0.221654626")


def _record_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _mp_text(value: object, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def _floor_decimal(value: str, places: int) -> str:
    with localcontext() as context:
        context.prec = max(120, places + 30)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _guarded_lower_rational(root: float) -> tuple[int, int, str]:
    scale = 10**V_DECIMAL_PLACES
    with localcontext() as context:
        context.prec = 50
        scaled = Decimal(repr(root)) * scale
        numerator = int(scaled.to_integral_value(rounding=ROUND_FLOOR)) - 1
    if numerator <= 0:
        raise ValueError("root is too small for guarded decimal rounding")
    text = f"{numerator // scale}.{numerator % scale:0{V_DECIMAL_PLACES}d}"
    return numerator, scale, text


def _integer_fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    byte_count = max(1, (magnitude.bit_length() + 7) // 8)
    encoded = magnitude.to_bytes(byte_count, byteorder="big", signed=False)
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(encoded).hexdigest(),
    }


def _coefficient_fingerprint(coefficients: tuple[int, ...]) -> str:
    encoded = ",".join(str(value) for value in coefficients).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _load_saw_bound() -> tuple[str, str]:
    payload = json.loads(SAW_RESULT_PATH.read_text(encoding="utf-8"))
    value = payload["data"]["lower_bound_saw"]["certified_decimal_lower"]
    return str(value), str(SAW_RESULT_PATH.relative_to(ROOT))


def _scaled_polynomial(coefficients: tuple[int, ...], p: int, q: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * p**power * q ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def main() -> None:
    checks: list[dict[str, object]] = []

    _record_check(
        checks,
        "singleton_exact_criterion",
        exact_criterion_residual((1, 1, 1), 1, 6) == 0
        and exact_criterion_residual((1, 1, 1), 1, 7) < 0
        and exact_criterion_residual((1, 1, 1), 1, 5) > 0,
        "P=1 and Q=6, so kappa_{\u007b0\u007d}=6v and v_root=1/6 exactly",
    )

    enumeration_start = time.perf_counter()
    cube3_p, cube3_q = enumerate_pq_polynomials((3, 3, 3))
    enumeration_seconds = time.perf_counter() - enumeration_start
    cycle_count = 1 << (54 - 27 + 1)
    _record_check(
        checks,
        "cube3_enumeration_invariants",
        len(cube3_p) == len(cube3_q) == 55
        and cube3_p[0] == 1
        and cube3_q[0] == 0
        and all(value >= 0 for value in cube3_p)
        and all(value >= 0 for value in cube3_q)
        and sum(cube3_p) == cycle_count
        and sum(cube3_q) == 54 * cycle_count,
        "2^26 fixed-origin configurations give exact nonnegative P,Q; P(1)=2^28 and Q(1)=54*2^28",
    )
    cross_p, cross_q = exact_pq_at_rational((3, 3, 3), 1, 5)
    _record_check(
        checks,
        "cube3_enumeration_vs_transfer",
        cross_p == _scaled_polynomial(cube3_p, 1, 5)
        and cross_q == _scaled_polynomial(cube3_q, 1, 5),
        "independent exact spin enumeration and parity transfer agree at v=1/5",
    )

    box_results: list[dict[str, object]] = []
    for shape in BOXES:
        root_start = time.perf_counter()
        v_root, k_root = solve_box_root_float(shape)
        root_seconds = time.perf_counter() - root_start
        p_lower, q_lower, v_lower_text = _guarded_lower_rational(v_root)
        k_interval = atanh_rational_interval(p_lower, q_lower, dps=HIGH_DPS)
        certified_k_decimal = _floor_decimal(k_interval[0], K_DECIMAL_PLACES)
        origin = tuple((side - 1) // 2 for side in shape)
        boundary_edges = sum(boundary_multiplicities(shape))
        certification: dict[str, object]

        if shape in EXACT_RESIDUAL_BOXES:
            exact_start = time.perf_counter()
            residual = exact_criterion_residual(shape, p_lower, q_lower)
            exact_seconds = time.perf_counter() - exact_start
            _record_check(
                checks,
                f"exact_residual_{'x'.join(map(str, shape))}",
                residual < 0,
                f"exact integer sign(vQ-P)=- at v={v_lower_text}",
            )
            certification = {
                "method": "exact Python-integer parity transfer",
                "criterion_residual": _integer_fingerprint(residual),
                "elapsed_seconds_numerical_metadata": exact_seconds,
            }
        else:
            low_start = time.perf_counter()
            kappa_low_precision = kappa_mpf_from_v(
                shape, v_lower_text, dps=LOW_DPS
            )
            low_seconds = time.perf_counter() - low_start
            high_start = time.perf_counter()
            kappa_high_precision = kappa_mpf_from_v(
                shape, v_lower_text, dps=HIGH_DPS
            )
            high_seconds = time.perf_counter() - high_start
            relative_low = mpmath_relative_roundoff_bound(shape, LOW_DPS)
            relative_high = mpmath_relative_roundoff_bound(shape, HIGH_DPS)
            with mp.workdps(HIGH_DPS + 20):
                agreement = abs(kappa_low_precision - kappa_high_precision)
                agreement_allowance = (
                    4
                    * max(abs(kappa_low_precision), abs(kappa_high_precision))
                    * relative_low
                )
                rigorous_kappa_upper = kappa_high_precision / (1 - relative_high)
                criterion_margin = 1 - rigorous_kappa_upper
            _record_check(
                checks,
                f"dual_precision_{'x'.join(map(str, shape))}",
                agreement <= agreement_allowance,
                "mp.dps=50 and 80 agree within the conservative positive-recurrence roundoff allowance",
            )
            _record_check(
                checks,
                f"safe_criterion_{'x'.join(map(str, shape))}",
                rigorous_kappa_upper < 1,
                f"roundoff-enlarged kappa upper endpoint is below 1 by {_mp_text(criterion_margin, 12)}",
            )
            certification = {
                "method": "positive-term mpmath parity transfer with analytic forward-error allowance",
                "mpmath_dps": [LOW_DPS, HIGH_DPS],
                "kappa_dps_50": _mp_text(kappa_low_precision),
                "kappa_dps_80": _mp_text(kappa_high_precision),
                "absolute_precision_difference": _mp_text(agreement),
                "relative_roundoff_allowance_dps_50": _mp_text(relative_low),
                "relative_roundoff_allowance_dps_80": _mp_text(relative_high),
                "roundoff_enlarged_kappa_upper": _mp_text(rigorous_kappa_upper),
                "certified_margin_below_one": _mp_text(criterion_margin),
                "elapsed_seconds_numerical_metadata": {
                    "dps_50": low_seconds,
                    "dps_80": high_seconds,
                },
            }

        box_results.append(
            {
                "shape": list(shape),
                "origin_index_coordinates": list(origin),
                "sites": shape[0] * shape[1] * shape[2],
                "internal_bonds": (shape[0] - 1) * shape[1] * shape[2]
                + shape[0] * (shape[1] - 1) * shape[2]
                + shape[0] * shape[1] * (shape[2] - 1),
                "crossing_boundary_bonds": boundary_edges,
                "transfer_cross_section_sites": shape[0] * shape[1],
                "transfer_states": 1 << (shape[0] * shape[1]),
                "numerical_root": {
                    "v": format(v_root, ".17g"),
                    "K": format(k_root, ".17g"),
                    "method": "binary64 positive-term bisection; locator only",
                    "elapsed_seconds_numerical_metadata": root_seconds,
                },
                "certified_endpoint": {
                    "v_exact": f"{p_lower}/{q_lower}",
                    "v_decimal": v_lower_text,
                    "rounding": f"floor to {V_DECIMAL_PLACES} places, then subtract one unit in the last place",
                    "criterion": "v*Q(v)-P(v) < 0",
                    "K_exact": f"atanh({p_lower}/{q_lower})",
                    "K_directed_rounding_interval": list(k_interval),
                    "reported_decimal_lower": certified_k_decimal,
                },
                "certification": certification,
            }
        )

    certified_k_values = [
        mp.mpf(entry["certified_endpoint"]["reported_decimal_lower"])
        for entry in box_results
    ]
    _record_check(
        checks,
        "certified_sequence_monotone",
        all(
            right >= left
            for left, right in zip(certified_k_values, certified_k_values[1:])
        ),
        "certified K lower endpoints are non-decreasing over 1^3,2^3,3^3,4^3,4x4x8,4x4x16,4x4x32",
    )
    final_bound = certified_k_values[-1]
    _record_check(
        checks,
        "final_bound_benchmark_side",
        final_bound < BENCHMARK_KC,
        f"{_mp_text(final_bound, 20)} < 0.221654626 (benchmark is comparison only)",
    )

    saw_bound_text, saw_provenance = _load_saw_bound()
    saw_bound = mp.mpf(saw_bound_text)
    _record_check(
        checks,
        "comparison_with_saw_bound",
        saw_bound > final_bound,
        f"SAW bound {saw_bound_text} is stronger than Simon--Lieb {box_results[-1]['certified_endpoint']['reported_decimal_lower']}",
    )

    data = {
        "model": "nearest-neighbor ferromagnetic Ising model on the simple-cubic lattice",
        "coupling_convention": "K=beta*J, v=tanh(K), and the boundary-edge factor is tanh(K)",
        "criterion": {
            "P": "sum over even edge subsets F of B of v^|F|",
            "Q": "sum_x q_B(x) sum_{partial F={0,x}} v^|F|",
            "kappa": "v*Q(v)/P(v)",
            "safe_condition": "kappa_B(K)<1 implies exponential decay and K<Kc",
        },
        "arithmetic": {
            "cube_3_enumeration": "exact int64 histograms and Python-integer polynomial conversion",
            "small_box_certificates": "exact Python integers at rational v",
            "large_box_certificates": "positive-term mpmath recurrence at 50 and 80 decimal digits",
            "large_box_forward_error": "gamma_n bound with depth 32*(E+N+1), unit roundoff inflated to 10^(-(dps-5)), and an additional factor 8",
            "K_conversion": f"mpmath.iv directed rounding at {HIGH_DPS} decimal digits",
        },
        "cube3_exact_enumeration": {
            "shape": [3, 3, 3],
            "free_boundary_conditions": True,
            "configurations_explicitly_enumerated": 1 << 26,
            "global_spin_flip_multiplicity": 2,
            "elapsed_seconds_numerical_metadata": enumeration_seconds,
            "P_coefficients": list(cube3_p),
            "Q_coefficients": list(cube3_q),
            "P_sha256_comma_decimal": _coefficient_fingerprint(cube3_p),
            "Q_sha256_comma_decimal": _coefficient_fingerprint(cube3_q),
        },
        "boxes": box_results,
        "best_certified_simon_lieb_lower_bound": {
            "shape": list(BOXES[-1]),
            "exact": box_results[-1]["certified_endpoint"]["K_exact"],
            "decimal_lower": box_results[-1]["certified_endpoint"][
                "reported_decimal_lower"
            ],
        },
        "comparisons": {
            "saw_connective_constant_lower_bound": saw_bound_text,
            "saw_provenance": saw_provenance,
            "stronger_rigorous_lower_bound": "SAW/connective-constant",
            "published_benchmark_falsification_target": "0.221654626",
            "benchmark_role": "comparison only; not used to choose or fit any endpoint",
        },
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": {
                "mpmath_dps": [LOW_DPS, HIGH_DPS],
                "interval_dps": HIGH_DPS,
                "certified_v_decimal_places": V_DECIMAL_PLACES,
                "reported_K_decimal_places": K_DECIMAL_PLACES,
            },
        },
        "data": data,
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    print(
        "PASS e12_simon_lieb: certified Kc >= "
        f"{box_results[-1]['certified_endpoint']['reported_decimal_lower']} "
        f"from B={BOXES[-1]} (SAW bound remains stronger)"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e12_simon_lieb: {error}")
        raise
