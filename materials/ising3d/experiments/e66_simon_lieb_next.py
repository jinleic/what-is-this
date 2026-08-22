"""Exact next-step Simon--Lieb finite-set certificates and obstruction audit.

This experiment deliberately does not use the benchmark critical coupling to
select a geometry or endpoint.  It enumerates a fixed size-ordered family of
open boxes, certifies every displayed rational value with Python integers, and
records why that exact family cannot improve the current rigorous lower bound.

Run from the repository root:
    .venv/bin/python experiments/e66_simon_lieb_next.py
"""

from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp

from ising.rigorous_bounds.simon_lieb import (
    atanh_rational_interval,
    boundary_multiplicities,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
)

SCRIPT = "experiments/e66_simon_lieb_next.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "simon_lieb_next.json"
INCUMBENT_PATH = ROOT / "results" / "bounds" / "improved_bounds.json"
DPS = 90
REPORT_PLACES = 40

# Geometry search is independent of any critical-coupling benchmark.  These
# are all sorted triples 2 <= a <= b <= c <= 5 with transfer cross-section
# a*b <= 16, ordered by site count and then lexicographically.  The four
# longest prisms extend the saturated 4x4 family by site count alone.
SEARCH_SHAPES = tuple(
    sorted(
        {
            (a, b, c)
            for a in range(2, 5)
            for b in range(a, 5)
            for c in range(b, 6)
            if a * b <= 16
        }
        | {(4, 4, 8), (4, 4, 12), (4, 4, 16), (4, 4, 24)},
        key=lambda shape: (math.prod(shape), shape),
    )
)
# Decimal rationals are chosen below the numerical roots, never from K_c.
# They are intentionally modest-denominator exact points so a clean rerun is
# fast while still exposing the finite-family frontier.
CERTIFICATE_POINTS = {
    (2, 2, 2): (1813, 10_000),
    (2, 2, 3): (1856, 10_000),
    (2, 2, 4): (1866, 10_000),
    (2, 3, 3): (1901, 10_000),
    (2, 2, 5): (1875, 10_000),
    (2, 3, 4): (1912, 10_000),
    (3, 3, 3): (1947, 10_000),
    (2, 3, 5): (1923, 10_000),
    (2, 4, 4): (1925, 10_000),
    (2, 4, 5): (1936, 10_000),
    (3, 3, 4): (1961, 10_000),
    (3, 3, 5): (1972, 10_000),
    (3, 4, 4): (1975, 10_000),
    (3, 4, 5): (1987, 10_000),
    (4, 4, 4): (1991, 10_000),
    (4, 4, 5): (2005, 10_000),
    (4, 4, 8): (2017, 10_000),
    (4, 4, 12): (2022, 10_000),
    (4, 4, 16): (2023, 10_000),
    (4, 4, 24): (2023, 10_000),
}
# exact signs below are obstruction certificates for all tested shapes: none
# can prove the desired improvement by the scalar Simon--Lieb criterion.
CHALLENGE = (23, 110)
CONTROL_SHAPE = (2, 3, 2)
CONTROL_POINT = (2, 9)


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(encoded).hexdigest(),
    }


def _floor_decimal(value: str, places: int = REPORT_PLACES) -> str:
    with localcontext() as context:
        context.prec = max(DPS + 20, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _scaled_polynomial(coefficients: tuple[int, ...], p: int, q: int) -> int:
    degree = len(coefficients) - 1
    return sum(c * p**i * q ** (degree - i) for i, c in enumerate(coefficients))


def _load_incumbent() -> tuple[str, str, str]:
    payload = json.loads(INCUMBENT_PATH.read_text(encoding="utf-8"))
    candidate = next(
        item
        for item in payload["data"]["candidates"]
        if item["id"] == "published_exact_saw_c36"
    )
    certificate = candidate["certificate"]
    upper_v = str(certificate["v_exact_rational_bracket"][1])
    numerator, denominator = upper_v.split("/", maxsplit=1)
    return (
        str(payload["data"]["final_certified_interval"]["decimal_outward_rounded"][0]),
        numerator,
        denominator,
    )


def main() -> None:
    checks: list[dict[str, object]] = []
    incumbent_k, incumbent_v_upper_num, incumbent_v_den = _load_incumbent()
    incumbent_v_upper = Decimal(incumbent_v_upper_num) / Decimal(incumbent_v_den)
    challenge_p, challenge_q = CHALLENGE
    challenge_interval = atanh_rational_interval(challenge_p, challenge_q, dps=DPS)
    challenge_k_lower = _floor_decimal(challenge_interval[0])
    _check(
        checks,
        "challenge_strictly_exceeds_incumbent",
        Decimal(challenge_p) / challenge_q > incumbent_v_upper
        and Decimal(challenge_k_lower) > Decimal(incumbent_k),
        f"atanh(23/110) >= {challenge_k_lower} > incumbent {incumbent_k}",
    )

    # Independent route 1: exhaustive spin enumeration produces the complete
    # integer P,Q polynomials on 2x2x2.  Route 2 is the parity transfer at p/q.
    control_p_poly, control_q_poly = enumerate_pq_polynomials((2, 2, 2))
    transfer_p, transfer_q = exact_pq_at_rational((2, 2, 2), 1, 5)
    _check(
        checks,
        "control_polynomial_enumeration_vs_parity_transfer",
        transfer_p == _scaled_polynomial(control_p_poly, 1, 5)
        and transfer_q == _scaled_polynomial(control_q_poly, 1, 5),
        "complete spin enumeration and exact parity transfer agree coefficient evaluation at t=1/5",
    )
    cp, cq = CONTROL_POINT
    control_full_p, control_full_q = exact_pq_at_rational(CONTROL_SHAPE, cp, cq)
    control_residual = exact_criterion_residual(CONTROL_SHAPE, cp, cq)
    _check(
        checks,
        "control_full_values_vs_residual_stream",
        control_residual == cp * control_full_q - cq * control_full_p,
        "independent full P,Q streams and direct residual stream agree exactly on 2x3x2",
    )

    records: list[dict[str, object]] = []
    best: dict[str, object] | None = None
    for shape in SEARCH_SHAPES:
        if shape not in CERTIFICATE_POINTS:
            continue
        p, q = CERTIFICATE_POINTS[shape]
        safe_residual = exact_criterion_residual(shape, p, q)
        challenge_residual = exact_criterion_residual(shape, challenge_p, challenge_q)
        interval = atanh_rational_interval(p, q, dps=DPS)
        record = {
            "shape": list(shape),
            "sites": math.prod(shape),
            "internal_bonds": (shape[0] - 1) * shape[1] * shape[2]
            + shape[0] * (shape[1] - 1) * shape[2]
            + shape[0] * shape[1] * (shape[2] - 1),
            "crossing_bonds": sum(boundary_multiplicities(shape)),
            "transfer_cross_section_sites": shape[0] * shape[1],
            "safe_exact_t": f"{p}/{q}",
            "safe_K_interval": list(interval),
            "safe_K_decimal_lower": _floor_decimal(interval[0]),
            "safe_residual_tQ_minus_P": _fingerprint(safe_residual),
            "challenge_exact_t": f"{challenge_p}/{challenge_q}",
            "challenge_residual_tQ_minus_P": _fingerprint(challenge_residual),
            "challenge_passes": challenge_residual < 0,
        }
        _check(
            checks,
            f"safe_exact_{'x'.join(map(str, shape))}",
            safe_residual < 0,
            f"exact integer sign(tQ-P) is negative at t={p}/{q}",
        )
        _check(
            checks,
            f"challenge_rejected_{'x'.join(map(str, shape))}",
            challenge_residual >= 0,
            "exact integer sign(tQ-P) is nonnegative at t=23/110",
        )
        records.append(record)
        if best is None or Decimal(record["safe_K_decimal_lower"]) > Decimal(
            best["safe_K_decimal_lower"]
        ):
            best = record

    assert best is not None
    _check(
        checks,
        "search_order_is_size_then_lexicographic",
        list(SEARCH_SHAPES) == sorted(SEARCH_SHAPES, key=lambda s: (math.prod(s), s)),
        "geometry selection is fixed by size and lexicographic order, not by a benchmark",
    )
    _check(
        checks,
        "finite_family_does_not_improve_incumbent",
        Decimal(best["safe_K_decimal_lower"]) < Decimal(incumbent_k),
        f"best tested exact endpoint {best['safe_K_decimal_lower']} remains below {incumbent_k}",
    )

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "Python integers for every criterion sign; exhaustive integer polynomial control; mpmath.iv directed rounding",
            "mpmath_iv_dps": DPS,
        },
        "data": {
            "classification": "[UNRESOLVED] No Simon--Lieb improvement over the exact c_36 incumbent was found in the completed exact family",
            "model": "nearest-neighbour ferromagnetic Ising model on the simple-cubic lattice; K=beta*J and t=tanh(K)",
            "theorem": {
                "finite_set_coefficient": "phi_S(t)=t*sum_(x in S) q_S(x)<sigma_0 sigma_x>_(S,free)",
                "high_temperature_form": "phi_S(t)=t*Q_S(t)/P_S(t)",
                "safe_condition": "phi_S(t)<1 implies exponential decay, hence atanh(t)<K_c",
            },
            "search_policy": {
                "statement": "[COMPUTATION] fixed size-ordered open-box family; no benchmark enters geometry or endpoint selection",
                "shapes": [list(shape) for shape in SEARCH_SHAPES if shape in CERTIFICATE_POINTS],
                "largest_exact_transfer_cross_section_sites": max(
                    shape[0] * shape[1] for shape in CERTIFICATE_POINTS
                ),
            },
            "independent_exact_controls": {
                "enumerated_shape": [2, 2, 2],
                "P_coefficients": list(control_p_poly),
                "Q_coefficients": list(control_q_poly),
                "second_control_shape": list(CONTROL_SHAPE),
                "second_control_exact_t": f"{cp}/{cq}",
                "full_P_fingerprint": _fingerprint(control_full_p),
                "full_Q_fingerprint": _fingerprint(control_full_q),
                "residual_fingerprint": _fingerprint(control_residual),
            },
            "challenge_beyond_incumbent": {
                "exact_t": f"{challenge_p}/{challenge_q}",
                "K_directed_interval": list(challenge_interval),
                "K_decimal_lower": challenge_k_lower,
                "incumbent_K_decimal_lower": incumbent_k,
                "outcome": "[COMPUTATION] every tested shape has exact phi_S(23/110)>=1",
            },
            "certificates": records,
            "best_completed_exact_certificate": best,
            "resource_obstruction": {
                "statement": "[COMPUTATION] Every completed cross-section<=16 certificate fails before the incumbent; a wider exact transfer is required for this route",
                "state_scaling": "parity transfer has 2^(a*b) states for an a-by-b cross-section",
                "not_a_global_no_go": "[UNRESOLVED] failure of this finite family does not rule out larger boxes, diamonds, cylinders, or another exact contraction",
            },
            "final_lower_bound": {
                "status": "[THEOREM from EXTERNAL EXACT COMPUTATION] incumbent unchanged",
                "decimal_lower": incumbent_k,
                "improves_incumbent": False,
            },
        },
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    print(
        "PASS e66_simon_lieb_next: exact finite-family obstruction; "
        f"incumbent remains K_c >= {incumbent_k}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e66_simon_lieb_next: {error}")
        raise
