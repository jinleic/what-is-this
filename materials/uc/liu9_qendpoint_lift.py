#!/usr/bin/env python3
"""Exact inward-q lift and global endpoint-complement certifier for Liu H2.

For r=1-q, Liu's quotient-free raw gap is exactly

    G(r) = G(0) + r G1 + r^2 G2.

No entropy argument depends on q.  This module encloses those three coefficients
with Arb, combines the q=1 endpoint support theorem from ``liu9_qone_face``,
and certifies the full feasible layer through nonpositive mean multipliers.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_objective import (
    _ArbOps, _a3_feasible_enclosure, _semantic_bounds, arb_fraction,
    arb_hull, evaluate_arb,
)
from liu9_pilot import (
    ONE, ZERO, GradientUnavailable, _box_arbs, mean_corner_range,
)
from liu9_qone_face import (
    LAMBDAS, _cross_slope_ceiling, _endpoint_log_domain, _float_lower,
    _float_upper, _mu_lower, canonical_qone_box, contract_qone_mean,
    qone_feasible_shift_lower, qone_shift_lower, zero_mass_face,
    zero_mass_total_upper,
)
from liu9_survivors import _box_from_json, _load_report
from liu9_tube import (
    Box, _bisect, complement_box_bound, distance_squared_box,
)
from liu9_size_biased import centered_kernel_range, kernel_range

ctx.prec = max(ctx.prec, 320)

SOURCE_DEFAULT = (
    HERE / "verification/results/liu9-complement-residual-rho1-1728-100k.json")
OUTPUT_DEFAULT = HERE / "verification/results/liu9-qendpoint-global.json"
RHO = Fraction(1, 1728)
GLOBAL_LAMBDAS = tuple(Fraction(value) for value in (
    "0", "-0.1", "-0.2", "-0.4", "-0.8", "-0.9", "-1", "-1.2",
))
_OPS = _ArbOps(monotone_corners=True)
_BASE_CACHE: dict[tuple[Any, ...], tuple[arb, str]] = {}
_MASS_FACE_CACHE: dict[
    tuple[Any, ...], tuple[arb, dict[str, Any]]
] = {}
_ACTIVE_TARGET_CACHE: dict[tuple[Any, ...], arb] = {}
QONE_KERNEL_PATH = (
    HERE / "verification/results/liu9-size-biased-qone-kernel.json")
QONE_KERNEL_SHA256 = (
    "eb645792526578d115da8356f7d9378812164273900a60607691a8f32661f19a")
_QONE_KERNEL_VERIFIED = False


def _verify_qone_kernel() -> None:
    global _QONE_KERNEL_VERIFIED
    if _QONE_KERNEL_VERIFIED:
        return
    report = json.loads(QONE_KERNEL_PATH.read_text())
    claimed = report.pop("report_sha256", None)
    if (
        report.get("claim_status") != "PROVED"
        or claimed != QONE_KERNEL_SHA256
        or _canonical_digest(report) != claimed
    ):
        raise AssertionError("size-biased q-one kernel theorem is unavailable")
    _QONE_KERNEL_VERIFIED = True




def _cached_base(
    active: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
    feasible: bool,
) -> tuple[arb, str]:
    key = (
        tuple(active[index] for index in (0, 1, 6, 7, 8)),
        str(target_mean), str(beta), str(lam), feasible,
    )
    cached = _BASE_CACHE.get(key)
    if cached is not None:
        return cached
    if feasible:
        value = qone_feasible_shift_lower(
            active, target_mean, beta, lam)
    else:
        value = qone_shift_lower(active, target_mean, beta, lam)
    if len(_BASE_CACHE) >= 200000:
        _BASE_CACHE.clear()
    _BASE_CACHE[key] = value
    return value




def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _unit_interval(value: arb) -> arb:
    lo, hi = _semantic_bounds(value, arb(0), arb(1))
    return arb_hull(lo, hi)


def _masses(box: Box) -> tuple[arb, arb, arb]:
    values = _box_arbs(box)
    a1, a2 = _unit_interval(values[0]), _unit_interval(values[1])
    return a1, a2, _a3_feasible_enclosure(a1, a2)


def _entropy(point: arb) -> arb:
    return _OPS.entropy(point)


def _product_kernel(left: arb, right: arb) -> arb:
    return _OPS.entropy(_OPS.xy_argument(left, right))


def _protocol_kernel(left: arb, right: arb) -> arb:
    return _OPS.entropy(_OPS.pi_argument(left, right))


def _linear(masses: Sequence[arb], points: Sequence[arb]) -> arb:
    return sum(
        (mass * _entropy(point) for mass, point in zip(masses, points)),
        arb(0),
    )


def _bilinear(
    masses: Sequence[arb],
    left: Sequence[arb],
    right: Sequence[arb],
    kernel,
) -> arb:
    return sum((
        masses[i] * masses[j] * kernel(left[i], right[j])
        for i in range(3) for j in range(3)
    ), arb(0))


def _mean(masses: Sequence[arb], points: Sequence[arb]) -> arb:
    return sum((mass * point for mass, point in zip(masses, points)), arb(0))


def q_gap_coefficients(box: Box, beta: arb) -> tuple[arb, arb, arb, arb, arb]:
    """Outward enclosures of G0,G1,G2,M0,M1 for r=1-q.

    The identity is algebraic: EHXY is quadratic in the mixture weights while
    EHPI and EHX are affine.  All entropy arguments are q-independent.
    """
    values = _box_arbs(box)
    masses = _masses(box)
    p0 = values[3:6]
    p1 = values[6:9]

    u0 = _bilinear(masses, p0, p0, _product_kernel)
    u1 = _bilinear(masses, p1, p1, _product_kernel)
    cross = _bilinear(masses, p0, p1, _product_kernel)
    w0 = _bilinear(masses, p0, p0, _protocol_kernel)
    w1 = _bilinear(masses, p1, p1, _protocol_kernel)
    v0 = _linear(masses, p0)
    v1 = _linear(masses, p1)

    g0 = (1 - beta) * u1 + beta * w1 - v1
    g1 = 2 * (1 - beta) * (cross - u1) + beta * (w0 - w1) - (v0 - v1)
    g2 = (1 - beta) * (u0 - 2 * cross + u1)
    return g0, g1, g2, _mean(masses, p0), _mean(masses, p1)


def _mu_point(value: arb) -> arb:
    if value == 0:
        return arb(1)
    if value == 1:
        return arb(0)
    return -((1 - value) * (1 - value).log()) / value


def _mu_range(value: arb) -> arb:
    lo, hi = arb(value.lower()), arb(value.upper())
    return _mu_point(hi).union(_mu_point(lo))


def _size_biased_gap_lower(
    weights: Sequence[arb],
    points: Sequence[arb],
    mean: arb,
    beta: arb,
    nonnegative_kernel: bool = False,
    kernel_mean: Optional[arb] = None,
) -> arb:
    """Exact q-one gap in its cancellation-free size-biased kernel form."""
    if not mean.lower() > 0:
        raise ValueError("size-biased form needs a positive mean floor")
    evaluation_mean = mean if kernel_mean is None else kernel_mean
    numerator = arb(0)
    for i, x in enumerate(points):
        for j, y in enumerate(points):
            pair, _ = kernel_range(x, y, evaluation_mean, beta)
            try:
                centered = centered_kernel_range(
                    x, y, evaluation_mean, beta)
                if centered.lower() > pair.lower():
                    pair = centered
            except GradientUnavailable:
                pass
            if nonnegative_kernel and pair.lower() < 0:
                if pair.upper() < 0:
                    raise AssertionError("q-one kernel enclosure misses zero")
                pair = arb(0).union(pair.upper())
            numerator += weights[i] * weights[j] * pair
    return (numerator / mean).lower()


def _direct_qone_gap(
    masses: Sequence[arb],
    points: Sequence[arb],
    beta: arb,
) -> arb:
    product = _bilinear(masses, points, points, _product_kernel)
    protocol = _bilinear(masses, points, points, _protocol_kernel)
    entropy = _linear(masses, points)
    return (1 - beta) * product + beta * protocol - entropy


def _component_mean_range(box: Box, q_value: float) -> arb:
    component = list(box)
    component[2] = (q_value, q_value)
    value = mean_corner_range(tuple(component))
    if value is None:
        raise ValueError("mass box misses the simplex")
    return value


def _q_one_minus_q_range(q_lo: float, q_hi: float) -> arb:
    def value(q: float) -> arb:
        point = arb(q)
        return point * (1 - point)
    endpoints = (value(q_lo), value(q_hi))
    lower = min((item.lower() for item in endpoints), key=float)
    if q_lo <= 0.5 <= q_hi:
        upper = arb(1) / 4
    else:
        upper = max((item.upper() for item in endpoints), key=float)
    return arb(lower).union(arb(upper))


def _qone_decomposition_lower(box: Box, parameters) -> arb:
    """Two exact q-one decompositions of the full gap; take their best lower."""
    values = _box_arbs(box)
    masses = _masses(box)
    q = _unit_interval(values[2])
    qbar = 1 - q
    p0 = tuple(_unit_interval(value) for value in values[3:6])
    p1 = tuple(_unit_interval(value) for value in values[6:9])
    global_weights = tuple(
        _unit_interval(qbar * mass) for mass in masses
    ) + tuple(_unit_interval(q * mass) for mass in masses)
    global_points = p0 + p1
    global_mean = mean_corner_range(box)
    if global_mean is None:
        return arb(0)
    mean_lo = max(global_mean.lower(), parameters.mean.lower(), key=float)
    feasible_mean = arb(mean_lo).union(global_mean.upper())
    overall_interval = _size_biased_gap_lower(
        global_weights, global_points, feasible_mean, parameters.beta,
        nonnegative_kernel=True, kernel_mean=parameters.mean)
    overall = max(arb(0), arb(overall_interval), key=float)

    mean0 = _component_mean_range(box, 0.0)
    mean1 = _component_mean_range(box, 1.0)
    component0_proved = bool(mean0.lower() >= parameters.mean.upper())
    component1_proved = bool(mean1.lower() >= parameters.mean.upper())
    if mean0.lower() > 0:
        f0 = _size_biased_gap_lower(
            masses, p0, mean0, parameters.beta,
            nonnegative_kernel=component0_proved,
            kernel_mean=parameters.mean if component0_proved else None)
    else:
        f0 = _direct_qone_gap(masses, p0, parameters.beta).lower()
    if component0_proved:
        f0 = max(arb(0), arb(f0), key=float)
    if mean1.lower() > 0:
        f1 = _size_biased_gap_lower(
            masses, p1, mean1, parameters.beta,
            nonnegative_kernel=component1_proved,
            kernel_mean=parameters.mean if component1_proved else None)
    else:
        f1 = _direct_qone_gap(masses, p1, parameters.beta).lower()
    if component1_proved:
        f1 = max(arb(0), arb(f1), key=float)

    product00 = _bilinear(masses, p0, p0, _product_kernel)
    product11 = _bilinear(masses, p1, p1, _product_kernel)
    product01 = _bilinear(masses, p0, p1, _product_kernel)
    protocol00 = _bilinear(masses, p0, p0, _protocol_kernel)
    protocol11 = _bilinear(masses, p1, p1, _protocol_kernel)
    protocol01 = _bilinear(masses, p0, p1, _protocol_kernel)
    signed_product = product00 - 2 * product01 + product11
    signed_protocol = protocol00 - 2 * protocol01 + protocol11
    mixing = _q_one_minus_q_range(box[2][0], box[2][1])

    from_global = (
        overall + parameters.beta * mixing * signed_protocol).lower()
    weighted_components = qbar * f0 + q * f1
    from_components = (
        weighted_components
        - (1 - parameters.beta) * mixing * signed_product).lower()
    return max(from_global, from_components, key=float)


def _quadratic_remainder_lower(
    q_lo: float,
    q_hi: float,
    linear: arb,
    quadratic: arb,
) -> arb:
    """Minimize the common-r lower polynomial instead of decorrelating r,r²."""
    r_lo = arb(1) - arb(q_hi)
    r_hi = arb(1) - arb(q_lo)
    slope = arb(linear.lower())
    curvature = arb(quadratic.lower())

    def value(radius: arb) -> arb:
        return slope * radius + curvature * radius**2

    endpoints = (value(r_lo), value(r_hi))
    radius = r_lo.union(r_hi)
    fallback = (radius * linear + radius**2 * quadratic).lower()
    if curvature <= 0:
        return min((item.lower() for item in endpoints), key=float)
    if not curvature > 0:
        return fallback
    derivative_lo = slope + 2 * curvature * r_lo
    derivative_hi = slope + 2 * curvature * r_hi
    if derivative_lo >= 0:
        return endpoints[0].lower()
    if derivative_hi <= 0:
        return endpoints[1].lower()
    if not (derivative_lo < 0 and derivative_hi > 0):
        return fallback
    vertex = -(slope**2) / (4 * curvature)
    return vertex.lower()


def qendpoint_shift_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> tuple[arb, dict[str, Any]]:
    """Lower-bound gap+lambda*(mean-target) using the exact q polynomial."""
    q_lo, q_hi = box[2]
    if q_lo <= 0.0 or q_hi > 1.0:
        raise ValueError("q lift requires an interval strictly above q=0")
    if lam > 0:
        raise ValueError("feasible-set multiplier must be nonpositive")

    active = canonical_qone_box(
        tuple(box[index] for index in (0, 1, 6, 7, 8)))
    base, base_method = qone_shift_lower(
        active, target_mean, beta, lam)
    _, linear, quadratic, mean0, mean1 = q_gap_coefficients(box, beta)
    shifted_linear = linear + arb_fraction(lam) * (mean0 - mean1)
    remainder_lower = _quadratic_remainder_lower(
        q_lo, q_hi, shifted_linear, quadratic)
    lower = base + remainder_lower
    return lower, {
        "base_method": base_method,
        "base_lower": float(base),
        "linear_lower": float(shifted_linear.lower()),
        "linear_upper": float(shifted_linear.upper()),
        "quadratic_lower": float(quadratic.lower()),
        "quadratic_upper": float(quadratic.upper()),
        "r_lower": 1.0 - q_hi,
        "r_upper": 1.0 - q_lo,
        "lambda": str(lam),
    }


def _uniform_zero_margin(target_mean: arb, beta: arb) -> arb:
    """Simultaneous zero-support reduction margin at y0=1/32."""
    y0 = arb_fraction(Fraction(1, 32))
    two_one_minus_beta = 2 * (1 - beta)
    lam_uniform = two_one_minus_beta * (target_mean - y0) - 1
    penalty = two_one_minus_beta * (-(1 - y0).log())
    return lam_uniform * (-y0.log() + 1) - penalty


def _zero_reduced_box(box: Box) -> Optional[Box]:
    reduced = list(box)
    changed = False
    for index in (3, 4, 5, 6, 7, 8):
        if box[index][1] <= 1.0 / 32.0 and box[index] != (0.0, 0.0):
            reduced[index] = (0.0, 0.0)
            changed = True
    return tuple(reduced) if changed else None

def _single_mass_high_endpoint_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> arb:
    """Endpoint-factored F_lambda when a3=1 and P0's support tends to one."""
    if box[0] != (0.0, 0.0) or box[1] != (0.0, 0.0):
        raise ValueError("single-mass endpoint bound needs a1=a2=0")
    if box[5][1] != 1.0 or not _endpoint_log_domain(box[5][0]):
        raise ValueError("P0 high support must end at one")
    if box[8][1] > 0.5:
        raise ValueError("P1 partner must lie in [0,1/2]")

    q = arb(box[2][0]).union(arb(box[2][1]))
    qbar = 1 - q
    v = arb(box[8][0]).union(arb(box[8][1]))
    hv = _entropy(v)
    face_gap = q * (
        (2 * qbar * (1 - beta) - 1) * hv
        + q * (1 - beta) * _product_kernel(v, v)
        + beta * _protocol_kernel(v, v)
    )
    face_mean = qbar + q * v
    face_shift = (
        face_gap
        + arb_fraction(lam) * (face_mean - target_mean)
    )

    s_range = arb(1) - arb(box[5][0])
    if s_range.upper() <= 0:
        return face_shift.lower()
    s = arb(s_range.upper())
    a_one = (arb(2) - s).union(arb(2))
    a_two_lo = 2 - 2 * s + 2 * s**2 - s**3
    a_two = arb(a_two_lo.lower()).union(arb(2))
    log_coefficient = (
        qbar - qbar**2 * (1 - beta) * a_one
        - beta * qbar * a_two
    )
    r_one_upper = s * (2 - s)
    r_two_upper = s * a_two_lo
    mu_one = _mu_lower(r_one_upper).union(arb(1))
    mu_two = _mu_lower(r_two_upper).union(arb(1))
    mu_s = _mu_lower(s).union(arb(1))
    analytic = (
        qbar**2 * (1 - beta) * a_one * (mu_one - a_one.log())
        + beta * qbar * a_two * (mu_two - a_two.log())
        - qbar * mu_s
    )
    cross_rate = (
        2 * qbar * q * (1 - beta)
        * _cross_slope_ceiling(box[8][0], box[8][1], box[5][0])
    )
    shifted_linear = analytic - cross_rate - arb_fraction(lam) * qbar
    log_s_upper = s.log()
    if log_coefficient.upper() < 0:
        inside = (
            log_coefficient.upper() * log_s_upper
            + shifted_linear.lower()
        )
        remainder = s * inside if inside < 0 else arb(0)
    else:
        logarithmic = log_coefficient.upper() * s * log_s_upper
        linear = (
            s * shifted_linear.lower()
            if shifted_linear.lower() < 0 else arb(0)
        )
        remainder = logarithmic + linear
    return face_shift.lower() + remainder



def _cached_mass_face_lower(
    face: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> tuple[arb, dict[str, Any]]:
    key = (
        face[2], face[5], face[8],
        str(target_mean), str(beta), str(lam),
    )
    cached = _MASS_FACE_CACHE.get(key)
    if cached is not None:
        return cached
    lower, details = qendpoint_shift_lower(
        face, target_mean, beta, lam)
    if (
        face[5][1] == 1.0
        and _endpoint_log_domain(face[5][0])
        and face[8][1] <= 0.5
    ):
        factored = _single_mass_high_endpoint_lower(
            face, target_mean, beta, lam)
        if factored > lower:
            lower = factored
            details = {
                **details,
                "base_method": "single-mass-high-endpoint",
            }
    value = lower, details
    if len(_MASS_FACE_CACHE) >= 200000:
        _MASS_FACE_CACHE.clear()
    _MASS_FACE_CACHE[key] = value
    return value


def _refined_zero_mass_penalty(
    box: Box,
    face: Box,
    beta: arb,
    lam: Fraction,
) -> arb:
    """Face-aware one-sided mass loss, sharper than the universal bound."""
    s = arb(zero_mass_total_upper(box).upper())
    terms = evaluate_arb(_box_arbs(face), beta, monotone_corners=True)
    numerator_upper = terms.numerator.upper()
    entropy_lower = max(arb(0), terms.ehx.lower(), key=float)
    return (
        (2 * s - s**2) * numerator_upper
        + s * (arb(2).log() - entropy_lower)
        + abs(arb_fraction(lam)) * s
    )


def _implied_active_mean_target(box: Box, target_mean: arb) -> arb:
    """Necessary P1-mean floor from full feasibility and P0's box maximum."""
    key = (
        box[0], box[1], box[2], box[3], box[4], box[5],
        str(target_mean),
    )
    cached = _ACTIVE_TARGET_CACHE.get(key)
    if cached is not None:
        return cached
    p0_box = list(box)
    p0_box[2] = (0.0, 0.0)
    p0_mean = mean_corner_range(tuple(p0_box))
    if p0_mean is None:
        return arb(1)
    upper0 = arb(p0_mean.upper())
    q_lo, q_hi = box[2]
    values = [
        upper0 + (target_mean - upper0) / arb(q_lo),
        upper0 + (target_mean - upper0) / arb(q_hi),
    ]
    lower = min((value.lower() for value in values), key=float)
    answer = arb(lower)
    if len(_ACTIVE_TARGET_CACHE) >= 100000:
        _ACTIVE_TARGET_CACHE.clear()
    _ACTIVE_TARGET_CACHE[key] = answer
    return answer


def best_qendpoint_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
) -> tuple[arb, str, dict[str, Any]]:
    q_lo, q_hi = box[2]
    active = canonical_qone_box(
        tuple(box[index] for index in (0, 1, 6, 7, 8)))
    _, linear, quadratic, mean0, mean1 = q_gap_coefficients(box, beta)
    candidates = []

    for lam in GLOBAL_LAMBDAS:
        base, base_method = _cached_base(
            active, target_mean, beta, lam, feasible=False)
        shifted_linear = linear + arb_fraction(lam) * (mean0 - mean1)
        remainder_lower = _quadratic_remainder_lower(
            q_lo, q_hi, shifted_linear, quadratic)
        details = {
            "base_method": base_method,
            "base_lower": float(base),
            "linear_lower": float(shifted_linear.lower()),
            "linear_upper": float(shifted_linear.upper()),
            "quadratic_lower": float(quadratic.lower()),
            "quadratic_upper": float(quadratic.upper()),
            "r_lower": 1.0 - q_hi,
            "r_upper": 1.0 - q_lo,
            "lambda": str(lam),
            "active_target": None,
        }
        candidates.append((
            base + remainder_lower,
            f"qendpoint(lambda={lam})",
            details,
        ))

    # Full feasibility and M0<=1 imply
    # M1 >= 1-(1-target)/q >= 1-(1-target)/q_lo.
    # This second family can therefore use the sharper feasible q-one support
    # reductions without assuming that M1 reaches the original target.
    # This is deliberately not an expansion of a shifted full-q functional:
    # on the implied active-feasible set, lambda<=0 gives
    # G0+lambda*(M1-active_target) <= G0 pointwise.  The remainder is therefore
    # the raw exact polynomial; no r*lambda*(M0-M1) term is needed, regardless
    # of the sign of M0-M1.
    if q_lo >= 0.984375:
        active_target = _implied_active_mean_target(box, target_mean)
        contracted_active, _ = contract_qone_mean(active, active_target)
        if contracted_active is None:
            contracted_active = active
        raw_remainder_lower = _quadratic_remainder_lower(
            q_lo, q_hi, linear, quadratic)
        for lam in LAMBDAS[:5]:
            base, base_method = _cached_base(
                contracted_active, active_target, beta, lam, feasible=True)
            details = {
                "base_method": base_method,
                "base_lower": float(base),
                "linear_lower": float(linear.lower()),
                "linear_upper": float(linear.upper()),
                "quadratic_lower": float(quadratic.lower()),
                "quadratic_upper": float(quadratic.upper()),
                "r_lower": 1.0 - q_hi,
                "r_upper": 1.0 - q_lo,
                "lambda": str(lam),
                "active_target": str(active_target),
            }
            candidates.append((
                base + raw_remainder_lower,
                f"qactive(lambda={lam})",
                details,
            ))
    zero_reduced = _zero_reduced_box(box) if q_lo >= 0.875 else None
    zero_margin = _uniform_zero_margin(target_mean, beta)
    if zero_reduced is not None and zero_margin > 0:

        for lam in GLOBAL_LAMBDAS:
            if zero_margin + arb_fraction(lam) < 0:
                continue
            lower, details = qendpoint_shift_lower(
                zero_reduced, target_mean, beta, lam)
            details["zero_margin"] = str(zero_margin)
            candidates.append((
                lower,
                f"qzero(lambda={lam})",
                details,
            ))
    if zero_mass_total_upper(box).upper() <= arb(1) / 16:
        face = zero_mass_face(box)
        for lam in GLOBAL_LAMBDAS:
            lower, details = _cached_mass_face_lower(
                face, target_mean, beta, lam)
            lower -= _refined_zero_mass_penalty(
                box, face, beta, lam)
            candidates.append((
                lower,
                f"mass-zero(lambda={lam})",
                details,
            ))
    return max(candidates, key=lambda item: float(item[0]))
def _paired_mean_kernel_lower(box: Box, parameters) -> arb:
    """Pairwise cancellation bound: sum D_m = m*N-M*H <= m*gap."""
    values = _box_arbs(box)
    masses = _masses(box)
    q = _unit_interval(values[2])
    qbar = 1 - q
    p0 = tuple(_unit_interval(value) for value in values[3:6])
    p1 = tuple(_unit_interval(value) for value in values[6:9])
    state_mean = tuple(
        qbar * p0[i] + q * p1[i] for i in range(3))
    state_entropy = tuple(
        qbar * _entropy(p0[i]) + q * _entropy(p1[i])
        for i in range(3))
    total = arb(0)
    for i in range(3):
        for j in range(3):
            product_energy = (
                qbar**2 * _product_kernel(p0[i], p0[j])
                + qbar * q * _product_kernel(p0[i], p1[j])
                + q * qbar * _product_kernel(p1[i], p0[j])
                + q**2 * _product_kernel(p1[i], p1[j])
            )
            protocol_energy = (
                qbar * _protocol_kernel(p0[i], p0[j])
                + q * _protocol_kernel(p1[i], p1[j])
            )
            pair_energy = (
                (1 - parameters.beta) * product_energy
                + parameters.beta * protocol_energy
            )
            pair = (
                parameters.mean * pair_energy
                - (
                    state_mean[j] * state_entropy[i]
                    + state_mean[i] * state_entropy[j]
                ) / 2
            )
            total += masses[i] * masses[j] * pair
    return (total / parameters.mean).lower()



def _swap_components(box: Box) -> Box:
    """Simultaneous (q,P0,P1)->(1-q,P1,P0) symmetry."""
    q_lo, q_hi = box[2]
    return (
        box[0], box[1], (
            max(0.0, _float_lower(arb(1) - arb(q_hi))),
            min(1.0, _float_upper(arb(1) - arb(q_lo))),
        ),
        box[6], box[7], box[8], box[3], box[4], box[5],
    )


def _full_bound(box: Box, parameters) -> tuple[arb, arb, str, dict[str, Any]]:
    ordinary = complement_box_bound(box, parameters, GLOBAL_LAMBDAS)
    gap_lower = ordinary.gap_lower
    objective_lower = ordinary.objective_lower
    method = ordinary.method
    details: dict[str, Any] = {}
    _verify_qone_kernel()
    support_width = max(hi - lo for lo, hi in box[3:])
    if support_width <= 1.0 / 16.0:
        decomposition = _qone_decomposition_lower(box, parameters)
        if decomposition > gap_lower:
            gap_lower = decomposition
            method = "qone-decomposition"
    guidance_lower = -math.inf
    orientations = []
    if box[2][0] > 0.0:
        orientations.append((box, False))
    if box[2][1] < 1.0:
        orientations.append((_swap_components(box), True))
    for oriented, swapped in orientations:
        endpoint, endpoint_method, endpoint_details = best_qendpoint_lower(
            oriented, parameters.mean, parameters.beta)
        endpoint_details = {**endpoint_details, "swapped": swapped}
        if float(endpoint) > guidance_lower:
            guidance_lower = float(endpoint)
            details = endpoint_details
        if endpoint > gap_lower:
            gap_lower = endpoint
            method = (
                f"qswap({endpoint_method})" if swapped else endpoint_method)
            if gap_lower >= ZERO:
                objective_lower = max(objective_lower, ONE, key=float)
    return gap_lower, objective_lower, method, details


def _bisect_coordinate(box: Box, index: int) -> tuple[Box, Box]:
    lo, hi = box[index]
    mid = (lo + hi) / 2.0
    if not lo < mid < hi:
        candidates = [
            candidate for candidate, (lower, upper) in enumerate(box)
            if lower < (lower + upper) / 2.0 < upper
        ]
        if not candidates:
            raise ArithmeticError(f"unresolved point box: {box!r}")
        index = max(
            candidates, key=lambda candidate: box[candidate][1] - box[candidate][0])
        lo, hi = box[index]
        mid = (lo + hi) / 2.0
    left, right = list(box), list(box)
    left[index] = (lo, mid)
    right[index] = (mid, hi)
    return tuple(left), tuple(right)


def _lift_bisect(
    box: Box,
    details: dict[str, Any],
    method: str,
) -> tuple[Box, Box]:
    """Split the failed variables named by the certificate, never its gauges."""
    if box[2] == (0.0, 1.0):
        return _bisect_coordinate(box, 2)
    if (
        method.startswith("mass-zero")
        or zero_mass_total_upper(box).upper() <= 3 * arb(1) / 64
    ):
        relevant = (0, 1, 2, 5, 8)
        index = max(relevant, key=lambda item: box[item][1] - box[item][0])
        if box[index][1] > box[index][0]:
            return _bisect_coordinate(box, index)
        raise ArithmeticError("unresolved exact zero-mass face")
    if details:
        if details["base_lower"] <= 0.0:
            swapped = bool(details.get("swapped"))
            support = (3, 4, 5) if swapped else (6, 7, 8)
            inactive = (6, 7, 8) if swapped else (3, 4, 5)
            active = (0, 1, *support)
            other = (2, *inactive)
            active_index = max(
                active, key=lambda item: box[item][1] - box[item][0])
            other_index = max(
                other, key=lambda item: box[item][1] - box[item][0])
            active_width = box[active_index][1] - box[active_index][0]
            other_width = box[other_index][1] - box[other_index][0]
            return _bisect_coordinate(
                box,
                active_index if active_width > other_width / 8 else other_index,
            )
        q_width = box[2][1] - box[2][0]
        other_width = max(
            box[index][1] - box[index][0]
            for index in (0, 1, 3, 4, 5, 6, 7, 8)
        )
        if q_width > other_width / 8.0:
            return _bisect_coordinate(box, 2)
    widths = tuple(hi - lo for lo, hi in box)
    return _bisect_coordinate(
        box, max(range(len(box)), key=lambda index: widths[index]))


def _initial_boxes(source: Path) -> tuple[list[Box], dict[str, Any]]:
    report, metadata = _load_report(source)
    boxes = [_box_from_json(saved["box"]) for saved in report["residual"]]
    return boxes, {
        "path": str(source.relative_to(HERE)) if source.is_relative_to(HERE) else str(source),
        "report_sha256": metadata["report_sha256"],
        "residual_boxes": len(boxes),
    }


def _active_group_key(box: Box) -> tuple[tuple[float, float], ...]:
    return tuple(box[index] for index in (0, 1, 6, 7, 8))


def run(
    source: Path,
    budget: int,
    seconds: float,
    checkpoint: Optional[Path],
    slice_index: int = 0,
    slice_count: int = 1,
) -> dict[str, Any]:
    all_initial, source_meta = _initial_boxes(source)
    if slice_count <= 0 or not 0 <= slice_index < slice_count:
        raise ValueError("slice requires 0 <= index < count")
    active_groups = sorted({_active_group_key(box) for box in all_initial})
    selected_groups = set(active_groups[slice_index::slice_count])
    initial = [
        box for box in all_initial if _active_group_key(box) in selected_groups
    ]
    source_meta.update({
        "active_groups_total": len(active_groups),
        "active_groups_selected": len(selected_groups),
        "slice_index": slice_index,
        "slice_count": slice_count,
        "slice_source_boxes": len(initial),
    })
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    rho_squared = arb_fraction(RHO * RHO)
    heap: list[tuple[int, int, float, int, int, Box]] = []
    serial = 0
    for box in initial:
        endpoint_priority = 0 if box[2][1] == 1.0 else 1
        heapq.heappush(
            heap, (0, endpoint_priority, -10.0, serial, 0, box))
        serial += 1

    started = time.monotonic()
    processed = split = mean_infeasible = tube_cleared = 0
    gap_cleared = objective_cleared = endpoint_cleared = 0
    maximum_depth = 0
    method_counts: dict[str, int] = {}

    while heap and processed < budget and time.monotonic() - started < seconds:
        _, _, _, _, depth, box = heapq.heappop(heap)
        processed += 1
        maximum_depth = max(maximum_depth, depth)
        mean_range = mean_corner_range(box)
        if mean_range is None or mean_range.upper() < parameters.mean.lower():
            mean_infeasible += 1
            continue
        distance = distance_squared_box(box, parameters)
        if distance.upper() < rho_squared:
            tube_cleared += 1
            continue
        gap_lower, objective_lower, method, details = _full_bound(
            box, parameters)
        method_counts[method] = method_counts.get(method, 0) + 1
        if gap_lower >= ZERO:
            gap_cleared += 1
            endpoint_cleared += int(method.startswith("qendpoint"))
            continue
        if objective_lower >= ONE:
            objective_cleared += 1
            continue

        left, right = _lift_bisect(box, details, method)
        split += 1
        priority = float(gap_lower)
        for child in (left, right):
            serial += 1
            # Endpoint boxes first: they contain the only non-smooth face.
            endpoint_priority = 0 if child[2][1] == 1.0 else 1
            heapq.heappush(heap, (
                -(depth + 1), endpoint_priority, priority, serial,
                depth + 1, child,
            ))

        if checkpoint and processed % 10000 == 0:
            checkpoint.write_text(json.dumps({
                "processed": processed,
                "frontier": len(heap),
                "split": split,
                "mean_infeasible": mean_infeasible,
                "tube_cleared": tube_cleared,
                "gap_cleared": gap_cleared,
                "endpoint_cleared": endpoint_cleared,
                "maximum_depth": maximum_depth,
                "current_depth": depth,
                "current_method": method,
                "current_gap_lower": float(gap_lower),
                "current_box": [[float(lo), float(hi)] for lo, hi in box],
                "elapsed_seconds": time.monotonic() - started,
            }, indent=1, sort_keys=True) + "\n")

    # Pending boxes themselves form a conservative residual cover.  Re-running
    # every expensive endpoint enclosure cannot promote an incomplete run.
    residual = [{
        "box": [[float(lo), float(hi)] for lo, hi in box],
        "depth": depth,
    } for _, _, _, _, depth, box in sorted(heap)]

    elapsed = time.monotonic() - started
    report: dict[str, Any] = {
        "tool": "liu9_qendpoint_lift.py",
        "claim_status": "PROVED" if not residual else "NUMERICAL",
        "source": source_meta,
        "rho": str(RHO),
        "run": {
            "processed": processed,
            "split": split,
            "mean_infeasible": mean_infeasible,
            "tube_cleared": tube_cleared,
            "gap_cleared": gap_cleared,
            "objective_cleared": objective_cleared,
            "endpoint_cleared": endpoint_cleared,
            "residual_count": len(residual),
            "maximum_depth": maximum_depth,
            "elapsed_seconds": elapsed,
            "stop_reason": "complete" if not residual else (
                "box_budget" if processed >= budget else "time_cap"),
            "method_counts": method_counts,
        },
        "identity": (
            "gap(q)=G0+(1-q)G1+(1-q)^2G2 exactly; "
            "mean(q)=M1+(1-q)(M0-M1)"
        ),
        "scope": (
            "the source residual cover, including q=1, with the proved "
            "rho=1/1701 tube used conservatively only at rho=1/1728"
        ),
        "residual": residual,
        "limitations": [] if not residual else [
            "An unfinished cover is not a proof of Hypothesis 2.",
        ],
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    parser.add_argument("--budget", type=int, default=500001)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--slice-index", type=int, default=0)
    parser.add_argument("--slice-count", type=int, default=1)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = run(
        args.source, args.budget, args.seconds, args.checkpoint,
        args.slice_index, args.slice_count)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    data = report["run"]
    print("LIU H2 EXACT Q-ENDPOINT LIFT")
    print("%s: processed=%d gap=%d endpoint=%d tube=%d mean=%d residual=%d depth=%d elapsed=%.2fs"
          % (report["claim_status"], data["processed"], data["gap_cleared"],
             data["endpoint_cleared"], data["tube_cleared"],
             data["mean_infeasible"], data["residual_count"],
             data["maximum_depth"], data["elapsed_seconds"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if not report["residual"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
