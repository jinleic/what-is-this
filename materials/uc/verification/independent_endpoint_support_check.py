#!/usr/bin/env python3
"""Independent arithmetic check of the endpoint support and exact q lift."""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_boundary_layer import gap_mp

ctx.prec = max(ctx.prec, 320)

ENDPOINT_REPORT = HERE / "results/liu9-endpoint-support.json"
OUTPUT = HERE / "results/liu9-endpoint-support-independent.json"


def _digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _arbf(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _h_point(value: arb) -> arb:
    if value == 0 or value == 1:
        return arb(0)
    return -(value * value.log() + (1 - value) * (1 - value).log())


def _h_increasing(interval: arb) -> arb:
    if interval.lower() < 0 or interval.upper() > arb(1) / 2:
        raise ValueError("independent entropy range expects [0,1/2]")
    return _h_point(arb(interval.lower())).union(_h_point(arb(interval.upper())))


def _pi(left: arb, right: arb) -> arb:
    return left * right * (1 + (1 - left) * (1 - right))


def _kernel(left: arb, right: arb, beta: arb) -> arb:
    product = left * right
    protocol = _pi(left, right)
    return (1 - beta) * _h_increasing(product) + beta * _h_increasing(protocol)


def _mu_lower(radius: Fraction) -> arb:
    r = _arbf(radius)
    return -((1 - r) * (1 - r).log()) / r


def _cross_max(x_lo: arb, x_hi: arb, y_lo: arb) -> arb:
    z_lo, z_hi = y_lo * x_lo, y_lo * x_hi

    def g(value: arb) -> arb:
        return value * ((1 - value) / value).log()

    def derivative(value: arb) -> arb:
        return ((1 - value) / value).log() - 1 / (1 - value)

    if derivative(z_hi) > 0:
        maximum = g(z_hi)
    elif derivative(z_lo) < 0:
        maximum = g(z_lo)
    else:
        raise AssertionError("independent blocker ranges should miss the critical point")
    return maximum / y_lo


def independent_blocker_lower(parameters) -> arb:
    # Necessary cuts derived without the production contractor:
    # a1(1-b1)<=1-m and max mean at b5=255/256 is 5053/8192<m.
    deficit = arb(1) - parameters.mean
    a1_lo = _arbf(Fraction(13, 32))
    a1_hi = deficit / _arbf(Fraction(15, 16))
    a1 = arb(a1_lo.lower()).union(arb(a1_hi.upper()))
    a2 = arb(0).union(_arbf(Fraction(1, 32)))
    a3_lo = 1 - a1_hi - _arbf(Fraction(1, 32))
    a3_hi = 1 - a1_lo
    a3 = arb(a3_lo.lower()).union(arb(a3_hi.upper()))
    b1_lo = 1 - deficit / a1_lo
    b1 = arb(b1_lo.lower()).union(_arbf(Fraction(1, 16)))
    b3 = _arbf(Fraction(13, 32)).union(_arbf(Fraction(7, 16)))
    y_lo = _arbf(Fraction(255, 256))
    s_upper = Fraction(1, 256)
    beta = parameters.beta
    lam = _arbf(Fraction(-1, 64))

    # Exact b5=1 face grouping.
    h1, h3 = _h_increasing(b1), _h_increasing(b3)
    face_gap = (
        (2 * a3 - 1) * (a1 * h1 + a2 * h3)
        + a1**2 * _kernel(b1, b1, beta)
        + 2 * a1 * a2 * _kernel(b1, b3, beta)
        + a2**2 * _kernel(b3, b3, beta)
    )
    face_mean = a1 * b1 + a2 * b3 + a3
    face_lower = (face_gap + lam * (face_mean - parameters.mean)).lower()

    # Independent s*log(s) combination for lowering b5 from one.
    s = _arbf(s_upper)
    a_one = _arbf(2 - s_upper).union(arb(2))
    a_two_lo_fraction = (
        2 - 2 * s_upper + 2 * s_upper**2 - s_upper**3)
    a_two = _arbf(a_two_lo_fraction).union(arb(2))
    abar = (1 - beta) * a_one + beta * a_two
    coefficient = a3 - a3**2 * abar
    r1 = s_upper * (2 - s_upper)
    r2 = s_upper * a_two_lo_fraction
    mu1 = _mu_lower(r1).union(arb(1))
    mu2 = _mu_lower(r2).union(arb(1))
    mus = _mu_lower(s_upper).union(arb(1))
    analytic = a3**2 * (
        (1 - beta) * a_one * (mu1 - a_one.log())
        + beta * a_two * (mu2 - a_two.log())
    ) - a3 * mus
    cross_rate = 2 * a3.upper() * (
        a1.upper() * _cross_max(
            arb(b1.lower()), arb(b1.upper()), y_lo)
        + a2.upper() * _cross_max(
            arb(b3.lower()), arb(b3.upper()), y_lo)
    )
    shifted = analytic - cross_rate - lam * a3
    if not coefficient.upper() < 0:
        raise AssertionError("independent endpoint coefficient lost its sign")
    inside = coefficient.upper() * s.log() + shifted.lower()
    remainder = s * inside if inside < 0 else arb(0)
    return face_lower + remainder


def _h_mp(value: mpmath.mpf) -> mpmath.mpf:
    if value == 0 or value == 1:
        return mpmath.mpf(0)
    return -value * mpmath.log(value) - (1 - value) * mpmath.log(1 - value)


def _q_coefficients(values, beta):
    a1, a2, _, *supports = values
    masses = (a1, a2, 1 - a1 - a2)
    p0, p1 = supports[:3], supports[3:]

    def energy(left, right, protocol=False):
        total = mpmath.mpf(0)
        for i in range(3):
            for j in range(3):
                argument = left[i] * right[j]
                if protocol:
                    argument *= 1 + (1 - left[i]) * (1 - right[j])
                total += masses[i] * masses[j] * _h_mp(argument)
        return total

    def entropy(points):
        return sum(masses[i] * _h_mp(points[i]) for i in range(3))

    u0, u1 = energy(p0, p0), energy(p1, p1)
    cross = energy(p0, p1)
    w0, w1 = energy(p0, p0, True), energy(p1, p1, True)
    v0, v1 = entropy(p0), entropy(p1)
    return (
        (1 - beta) * u1 + beta * w1 - v1,
        2 * (1 - beta) * (cross - u1) + beta * (w0 - w1) - (v0 - v1),
        (1 - beta) * (u0 - 2 * cross + u1),
    )


def main() -> int:
    parameters_mp = solve_equation_parameters(100)
    parameters = certify_equation_parameters(parameters_mp)
    independent_lower = independent_blocker_lower(parameters)
    if not independent_lower > 0:
        raise AssertionError("independent endpoint lower is not positive")

    production = json.loads(ENDPOINT_REPORT.read_text())
    claimed = production.pop("report_sha256")
    if _digest(production) != claimed or production["claim_status"] != "PROVED":
        raise AssertionError("production endpoint report failed authentication")

    with mpmath.workdps(100):
        base = [
            mpmath.mpf("0.40625"), mpmath.mpf("0.001"), mpmath.mpf(0),
            mpmath.mpf("0.2"), mpmath.mpf("0.4"), mpmath.mpf("0.8"),
            mpmath.mpf("0.0625"), mpmath.mpf("0.42"), mpmath.mpf("0.999"),
        ]
        g0, g1, g2 = _q_coefficients(tuple(base), parameters_mp.beta)
        maximum_residual = mpmath.mpf(0)
        for q in (mpmath.mpf("0.5"), mpmath.mpf("0.9375"),
                  mpmath.mpf("0.99"), mpmath.mpf(1)):
            base[2] = q
            exact = gap_mp(tuple(base), parameters_mp.beta)
            radius = 1 - q
            residual = abs(exact - (g0 + radius * g1 + radius**2 * g2))
            maximum_residual = max(maximum_residual, residual)
        if maximum_residual > mpmath.mpf("1e-90"):
            raise AssertionError("independent exact-q identity failed")

    report = {
        "tool": "independent_endpoint_support_check.py",
        "claim_status": "PASS",
        "production_report_sha256": claimed,
        "independent_blocker_lower": str(independent_lower),
        "independent_blocker_lower_float": float(independent_lower.lower()),
        "q_polynomial_max_residual": mpmath.nstr(maximum_residual, 20),
        "q_polynomial_points": 4,
        "imports_production_endpoint_module": False,
        "imports_production_q_lift_module": False,
        "uses_sanctioned_gap_mp": True,
    }
    report["report_sha256"] = _digest(report)
    OUTPUT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("INDEPENDENT LIU H2 ENDPOINT SUPPORT CHECK")
    print("PASS blocker lower %.12g; exact-q residual %s"
          % (report["independent_blocker_lower_float"],
             report["q_polynomial_max_residual"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
