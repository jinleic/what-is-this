#!/usr/bin/env python3
"""Specialized q=1-face branch-and-bound for Liu Hypothesis 2.

The rho=1/4096 complement frontier contains 48,996 boxes whose q interval
reaches one but not zero.  Their three P0 support coordinates are inactive
exactly at q=1; retaining those gauges duplicates only 230 distinct active
five-dimensional boxes by a factor above 200 and disables the centered gradient
bound whenever a gauge interval touches an entropy endpoint.

This module verifies the source report digest, projects every such box to
(a1,a2,b1,b3,b5), canonicalizes q=1 and the inactive supports to 1/2, then runs
a deterministic Arb branch-and-bound splitting ONLY those five active
coordinates.  Every discard is rigorous.  An unfinished frontier is labelled
NUMERICAL and is not a proof of the face or of Hypothesis 2.

The raw gap and objective ratio are kept separate.  Boundary laws with EHX=0
may have raw gap zero while the ratio is undefined; no ratio value is invented.
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
    _ArbOps, _a3_feasible_enclosure, arb_fraction, evaluate_arb, mean_value,
)
from liu9_pilot import (
    ONE, ZERO, GradientUnavailable, _box_arbs, _gap_mean_gradient,
    mean_corner_range,
)
from liu9_survivors import _box_from_json, _load_report
from liu9_tube import Box, _centered_shift_enclosure, complement_box_bound

ctx.prec = max(ctx.prec, 320)

SOURCE_DEFAULT = (
    HERE / "verification/results/liu9-complement-residual-rho1-1728-100k.json")
OUTPUT_DEFAULT = HERE / "verification/results/liu9-qone-face.json"
ACTIVE = (0, 1, 6, 7, 8)
ACTIVE_NAMES = ("a1", "a2", "b1", "b3", "b5")
LAMBDAS = tuple(Fraction(v) for v in (
    "0", "-0.015625", "-0.03125", "-0.046875", "-0.06249",
    "-0.1", "-0.2", "-0.4", "-0.8", "-0.9", "-1", "-1.2",
))
ZERO_LAYER = Fraction(1, 32)
_QONE_MASS_CACHE: dict[
    tuple[Any, ...], tuple[arb, str]
] = {}
HALF = arb_fraction(Fraction(1, 2))
_OPS = _ArbOps(monotone_corners=True)


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _active_key(box: Box) -> tuple[tuple[float, float], ...]:
    return tuple(tuple(float(v) for v in box[i]) for i in ACTIVE)


def canonical_qone_box(key: Sequence[Sequence[float]]) -> Box:
    a1, a2, b1, b3, b5 = (tuple(pair) for pair in key)
    gauge = (0.5, 0.5)
    return (a1, a2, (1.0, 1.0), gauge, gauge, gauge, b1, b3, b5)


def active_widths(box: Box) -> tuple[float, ...]:
    return tuple(box[index][1] - box[index][0] for index in ACTIVE)


def bisect_active(box: Box) -> tuple[Box, Box]:
    widths = active_widths(box)
    widest = max(widths)
    active_position = next(i for i, width in enumerate(widths) if width == widest)
    coordinate = ACTIVE[active_position]
    lo, hi = box[coordinate]
    mid = (lo + hi) / 2.0
    left, right = list(box), list(box)
    left[coordinate] = (lo, mid)
    right[coordinate] = (mid, hi)
    return tuple(left), tuple(right)


def active_compact(box: Box) -> list[list[float]]:
    return [[float(box[i][0]), float(box[i][1])] for i in ACTIVE]


def _float_lower(value: arb) -> float:
    return math.nextafter(float(value.lower()), -math.inf)


def _float_upper(value: arb) -> float:
    return math.nextafter(float(value.upper()), math.inf)


def _simplex_vertices_fraction(box: Box) -> tuple[tuple[Fraction, Fraction], ...]:
    l1, u1 = (Fraction.from_float(value) for value in box[0])
    l2, u2 = (Fraction.from_float(value) for value in box[1])
    candidates = {
        (x, y) for x in (l1, u1) for y in (l2, u2)
        if x + y <= 1
    }
    for x in (l1, u1):
        y = 1 - x
        if l2 <= y <= u2:
            candidates.add((x, y))
    for y in (l2, u2):
        x = 1 - y
        if l1 <= x <= u1:
            candidates.add((x, y))
    return tuple(sorted(candidates))




def contract_qone_mean(box: Box, target_mean: arb) -> tuple[Optional[Box], int]:
    """Contract q=1 by the exact nonnegative mean deficit, always outward."""
    answer = list(box)
    contractions = 0
    deficit = arb(1) - arb(target_mean.lower())

    for _ in range(4):
        changed = False
        a1_lo, a1_hi = answer[0]
        a2_lo, a2_hi = answer[1]
        b1_lo, b1_hi = answer[6]
        b3_lo, b3_hi = answer[7]

        new_a1_hi = min(
            a1_hi, _float_upper(arb(1) - arb(a2_lo)))
        if b1_hi < 1.0:
            new_a1_hi = min(new_a1_hi, _float_upper(
                deficit / (arb(1) - arb(b1_hi))))
        new_a2_hi = min(
            a2_hi, _float_upper(arb(1) - arb(a1_lo)))
        if b3_hi < 1.0:
            new_a2_hi = min(new_a2_hi, _float_upper(
                deficit / (arb(1) - arb(b3_hi))))
        if new_a1_hi < a1_lo or new_a2_hi < a2_lo:
            return None, contractions
        if new_a1_hi < a1_hi:
            answer[0] = (a1_lo, new_a1_hi)
            contractions += 1
            changed = True
        if new_a2_hi < a2_hi:
            answer[1] = (a2_lo, new_a2_hi)
            contractions += 1
            changed = True

        a1_lo, _ = answer[0]
        a2_lo, _ = answer[1]
        if a1_lo > 0.0:
            new_b1_lo = max(b1_lo, _float_lower(
                arb(1) - deficit / arb(a1_lo)))
            if new_b1_lo > b1_hi:
                return None, contractions
            if new_b1_lo > b1_lo:
                answer[6] = (new_b1_lo, b1_hi)
                contractions += 1
                changed = True
        if a2_lo > 0.0:
            new_b3_lo = max(b3_lo, _float_lower(
                arb(1) - deficit / arb(a2_lo)))
            if new_b3_lo > b3_hi:
                return None, contractions
            if new_b3_lo > b3_lo:
                answer[7] = (new_b3_lo, b3_hi)
                contractions += 1
                changed = True

        # The admissible s=1-b5 is linear-fractional in (a1,a2), so its
        # maximum is at an exact-rational vertex of the clipped mass rectangle.
        b1_hi = answer[6][1]
        b3_hi = answer[7][1]
        s_max = arb(0)
        unconstrained = False
        vertices = _simplex_vertices_fraction(tuple(answer))
        if not vertices:
            return None, contractions
        for a1, a2 in vertices:
            a3 = 1 - a1 - a2
            remaining = (
                deficit
                - arb_fraction(a1) * (arb(1) - arb(b1_hi))
                - arb_fraction(a2) * (arb(1) - arb(b3_hi))
            )
            if a3 == 0:
                if remaining.upper() >= 0:
                    unconstrained = True
                    break
                continue
            ratio = remaining / arb_fraction(a3)
            if ratio.upper() > s_max.upper():
                s_max = ratio
        if not unconstrained:
            if s_max.upper() < 0:
                s_max = arb(0)
            new_b5_lo = max(
                answer[8][0], _float_lower(arb(1) - s_max.upper()))
            if new_b5_lo > answer[8][1]:
                return None, contractions
            if new_b5_lo > answer[8][0]:
                answer[8] = (new_b5_lo, answer[8][1])
                contractions += 1
                changed = True

        if not changed:
            break

    contracted = tuple(answer)
    if mean_corner_range(contracted) is None:
        return None, contractions
    return contracted, contractions


def _kernel_arb(x: arb, y: arb, beta: arb) -> arb:
    return (
        (1 - beta) * _OPS.entropy(_OPS.xy_argument(x, y))
        + beta * _OPS.entropy(_OPS.pi_argument(x, y))
    )


def _g_product(value: arb) -> arb:
    """z*h'(z), extended continuously by zero at z=0."""
    if value == 0:
        return arb(0)
    return value * ((1 - value) / value).log()


def _cross_slope_ceiling(x_lo: float, x_hi: float,
                         y_lo: float) -> arb:
    """Upper-bound x*h'(yx) for x in [x_lo,x_hi], y>=y_lo."""
    if x_hi <= 0.0:
        return arb(0)
    if x_lo < 0.0 or x_hi > 0.5 or not 0.0 < y_lo <= 1.0:
        raise ValueError("cross-slope domain must satisfy 0<=x<=1/2, 0<y<=1")
    y = arb(y_lo)
    z_lo = y * arb(x_lo)
    z_hi = y * arb(x_hi)
    if not z_hi > 0:
        return arb(0)

    def derivative(z: arb) -> arb:
        return ((1 - z) / z).log() - 1 / (1 - z)

    if derivative(z_hi) > 0:
        maximum = _g_product(z_hi)
    elif z_lo > 0 and derivative(z_lo) < 0:
        maximum = _g_product(z_lo)
    else:
        # g(z)=z log((1-z)/z) is strictly concave.  Its unique critical
        # point is bracketed here by derivative signs, so an interval
        # evaluation on this rational cell encloses the global maximum.
        critical = (
            arb_fraction(Fraction(2178, 10000))
            .union(arb_fraction(Fraction(2179, 10000)))
        )
        if not (derivative(arb_fraction(Fraction(2178, 10000))) > 0
                and derivative(arb_fraction(Fraction(2179, 10000))) < 0):
            raise AssertionError("product-entropy critical bracket failed")
        maximum = _g_product(critical)
    return maximum / arb(y_lo)


def _mu_lower(radius: float | arb) -> arb:
    """Lower endpoint of mu(r)=-(1-r)log(1-r)/r on [0,radius]."""
    r = arb(radius)
    if r <= 0:
        return arb(1)
    return -((1 - r) * (1 - r).log()) / r


def _endpoint_log_domain(lower: float) -> bool:
    return bool(arb(1) - arb(lower) <= (-arb(1)).exp())


def _endpoint_remainder_lower(
    box: Box,
    beta: arb,
    lam: Fraction,
) -> arb:
    """Lower-bound F_lambda(b5)-F_lambda(b5=1), including b5=1.

    Here ``F_lambda = gap + lambda*(mean-target)``.  The two low supports
    lie in [0,1/2], the high support is ``b5=1-s`` with ``s<=1/2``.
    Entropy singularities are combined analytically before interval
    evaluation; no interval ever evaluates log(0).
    """
    a1 = arb(box[0][0]).union(arb(box[0][1]))
    a2 = arb(box[1][0]).union(arb(box[1][1]))
    a3 = _a3_feasible_enclosure(a1, a2)
    s_range = arb(1) - arb(box[8][0])
    if s_range.upper() <= 0:
        return arb(0)
    if s_range.upper() > arb(1) / 2:
        raise ValueError("endpoint remainder needs b5 in [1/2,1]")
    s = arb(s_range.upper())
    y_lower = box[8][0]

    cross_rate = 2 * a3.upper() * (
        a1.upper() * _cross_slope_ceiling(
            box[6][0], box[6][1], y_lower)
        + a2.upper() * _cross_slope_ceiling(
            box[7][0], box[7][1], y_lower)
    )

    a_one = (arb(2) - s).union(arb(2))
    a_two_lo = 2 - 2 * s + 2 * s**2 - s**3
    a_two = arb(a_two_lo.lower()).union(arb(2))
    abar = (1 - beta) * a_one + beta * a_two
    log_coefficient = a3 - a3**2 * abar

    r_one_upper = s * (2 - s)
    r_two_upper = s * a_two_lo
    mu_one = _mu_lower(r_one_upper).union(arb(1))
    mu_two = _mu_lower(r_two_upper).union(arb(1))
    mu_s = _mu_lower(s).union(arb(1))
    analytic = a3**2 * (
        (1 - beta) * a_one * (mu_one - a_one.log())
        + beta * a_two * (mu_two - a_two.log())
    ) - a3 * mu_s
    shifted_linear = (
        analytic - cross_rate - arb_fraction(lam) * a3
    )

    log_s_upper = s.log()
    if log_coefficient.upper() < 0:
        inside = (
            log_coefficient.upper() * log_s_upper
            + shifted_linear.lower()
        )
        return s * inside if inside < 0 else arb(0)

    # For nonnegative or sign-indefinite coefficient, s*log(s) is decreasing
    # on [0,s_upper] (all endpoint boxes use s_upper<=1/2; boxes above 1/e
    # are not sent here).  The largest coefficient gives the lower product.
    if not s <= (-arb(1)).exp():
        raise ValueError("positive log coefficient needs s<=1/e")
    logarithmic = log_coefficient.upper() * s * log_s_upper
    linear = (
        s * shifted_linear.lower()
        if shifted_linear.lower() < 0 else arb(0)
    )
    return logarithmic + linear


def qone_endpoint_shift_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> arb:
    """Endpoint-aware feasible-set lower bound for the q=1 raw gap."""
    if box[2] != (1.0, 1.0):
        raise ValueError("q-one endpoint bound requires q=1 exactly")
    if box[8][1] != 1.0:
        raise ValueError("endpoint bound requires a b5 interval ending at 1")
    if max(box[6][1], box[7][1]) > 0.5 or box[8][0] < 0.5:
        raise ValueError("endpoint ordering b1,b3<=1/2<=b5 is required")
    if lam > 0:
        raise ValueError("feasible-set multiplier must be nonpositive")

    a1 = arb(box[0][0]).union(arb(box[0][1]))
    a2 = arb(box[1][0]).union(arb(box[1][1]))
    a3 = _a3_feasible_enclosure(a1, a2)
    b1 = arb(box[6][0]).union(arb(box[6][1]))
    b3 = arb(box[7][0]).union(arb(box[7][1]))
    h1, h3 = _OPS.entropy(b1), _OPS.entropy(b3)
    low_entropy = a1 * h1 + a2 * h3
    face_gap = (
        (2 * a3 - 1) * low_entropy
        + a1**2 * _kernel_arb(b1, b1, beta)
        + 2 * a1 * a2 * _kernel_arb(b1, b3, beta)
        + a2**2 * _kernel_arb(b3, b3, beta)
    )
    face_mean = a1 * b1 + a2 * b3 + a3
    face_shift = (
        face_gap
        + arb_fraction(lam) * (face_mean - target_mean)
    )
    if lam == 0 and a3.lower() >= HALF:
        # Exact algebra: every displayed term in face_gap is nonnegative.
        face_lower = arb(0)
    else:
        face_lower = face_shift.lower()
    return face_lower + _endpoint_remainder_lower(box, beta, lam)


def best_qone_endpoint_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
) -> tuple[arb, str]:
    candidates = [
        (qone_endpoint_shift_lower(box, target_mean, beta, lam),
         f"endpoint(lambda={lam})")
        for lam in LAMBDAS
    ]
    return max(candidates, key=lambda item: float(item[0]))


def reduce_simplex_face(box: Box) -> tuple[Box, bool]:
    """Canonicalize cells whose simplex lower endpoints already sum to one.

    If a1>=l1, a2>=l2, a1+a2<=1 and l1+l2=1, feasibility forces
    a1=l1, a2=l2 and a3=0 exactly.  The third active support b5 is then a gauge.
    """
    l1, _ = box[0]
    l2, _ = box[1]
    if Fraction.from_float(l1) + Fraction.from_float(l2) != 1:
        return box, False
    answer = list(box)
    answer[0] = (l1, l1)
    answer[1] = (l2, l2)
    answer[8] = (0.5, 0.5)
    return tuple(answer), True


def reduce_zero_layer(box: Box) -> tuple[Box, int]:
    """Reduce active supports wholly in [0,1/32] to the proved zero face."""
    answer = list(box)
    reduced = 0
    for index in (6, 7, 8):
        lo, hi = answer[index]
        if lo == 0.0 and 0.0 < hi <= float(ZERO_LAYER):
            answer[index] = (0.0, 0.0)
            reduced += 1
    return tuple(answer), reduced


def initial_active_boxes(path: Path) -> tuple[list[Box], dict[str, Any]]:
    source, metadata = _load_report(path)
    keys = set()
    source_count = 0
    q_intervals: dict[str, int] = {}
    for saved in source["residual"]:
        box = _box_from_json(saved["box"])
        qlo, qhi = box[2]
        if float(qhi) != 1.0 or float(qlo) <= 0.0:
            continue
        source_count += 1
        qkey = str([float(qlo), float(qhi)])
        q_intervals[qkey] = q_intervals.get(qkey, 0) + 1
        keys.add(_active_key(box))
    boxes = [canonical_qone_box(key) for key in sorted(keys)]
    return boxes, {
        "source": metadata,
        "q_at_one_source_boxes": source_count,
        "unique_active_boxes": len(boxes),
        "duplication_factor": source_count / len(boxes),
        "q_intervals": q_intervals,
    }


def feasible(box: Box, target_mean: arb) -> bool:
    mean_range = mean_corner_range(box)
    return (mean_range is not None
            and mean_range.upper() >= target_mean.lower())


def qone_shift_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> tuple[arb, str]:
    """Lower-bound gap+lambda*(mean-target) without assuming feasibility."""
    if box[2] != (1.0, 1.0):
        raise ValueError("shifted q-one bound requires q=1 exactly")
    values = _box_arbs(box)
    terms = evaluate_arb(values, beta, monotone_corners=True)
    direct = (
        terms.numerator - terms.ehx
        + arb_fraction(lam) * (mean_value(values) - target_mean)
    )
    candidates = [(direct.lower(), f"shift-direct(lambda={lam})")]
    if (
        max(upper - lower for lower, upper in box[3:]) <= 0.25
        and all(
            lower == upper or (lower > 0.0 and upper < 1.0)
            for lower, upper in box[3:]
        )
    ):
        try:
            centered = _centered_shift_enclosure(
                box, target_mean, beta, arb_fraction(lam),
                _gap_mean_gradient(box, beta))
            candidates.append(
                (centered.lower(), f"shift-center(lambda={lam})"))
        except GradientUnavailable:
            pass
    if (
        box[8][1] == 1.0
        and _endpoint_log_domain(box[8][0])
        and max(box[6][1], box[7][1]) <= 0.5
    ):
        candidates.append((
            qone_endpoint_shift_lower(box, target_mean, beta, lam),
            f"shift-endpoint(lambda={lam})",
        ))
    return max(candidates, key=lambda item: float(item[0]))


def _single_zero_margin(target_mean: arb, beta: arb) -> arb:
    """Sharp one-coordinate boundary-layer margin at y0=1/16."""
    y0 = arb_fraction(Fraction(1, 16))
    two_one_minus_beta = 2 * (1 - beta)
    lam_feasible = two_one_minus_beta * target_mean - 1
    log_inverse = -y0.log()
    penalty = two_one_minus_beta * (-(1 - y0).log())
    return (
        lam_feasible * (log_inverse + 1)
        - two_one_minus_beta * y0 * (log_inverse / 2 + arb(3) / 4)
        - penalty
    )


def qone_feasible_shift_lower(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> tuple[arb, str]:
    """Shifted lower bound on the feasible subset, composing one zero layer."""
    candidates = [qone_shift_lower(box, target_mean, beta, lam)]
    margin = _single_zero_margin(target_mean, beta)
    if margin > 0 and margin + arb_fraction(lam) >= 0:
        for index in (6, 7):
            if box[index][1] <= 1.0 / 16.0:
                reduced = list(box)
                reduced[index] = (0.0, 0.0)
                lower, method = qone_shift_lower(
                    tuple(reduced), target_mean, beta, lam)
                candidates.append((
                    lower,
                    f"single-zero({index},{method})",
                ))
    return max(candidates, key=lambda item: float(item[0]))


def zero_mass_face(box: Box) -> Box:
    """Move masses a1,a2 to a3 and canonicalize their four support gauges."""
    answer = list(box)
    answer[0] = (0.0, 0.0)
    answer[1] = (0.0, 0.0)
    for index in (3, 4, 6, 7):
        answer[index] = (0.5, 0.5)
    return tuple(answer)


def zero_mass_total_upper(box: Box) -> arb:
    return arb(box[0][1]) + arb(box[1][1])


def zero_mass_shift_penalty(box: Box, lam: Fraction) -> arb:
    """One-sided loss for replacing (a1,a2,a3) by (0,0,1).

    Put s=a1+a2.  In each component, and hence in the global mixture, the
    original law is (1-s)Q+sR while the face law is Q.  For any kernel in
    [0,log(2)], its quadratic expectation can fall by at most
    (2s-s^2)log(2) <= 2s log(2); the entropy moment can rise by at most
    s log(2), and the mean changes by at most s.  Therefore F_lambda can fall
    by at most (3 log(2)+|lambda|)s.
    """
    s_upper = arb(zero_mass_total_upper(box).upper())
    return (
        3 * arb(2).log() + abs(arb_fraction(lam))
    ) * s_upper


def _cached_qone_mass_face(
    face: Box,
    target_mean: arb,
    beta: arb,
    lam: Fraction,
) -> tuple[arb, str]:
    key = (face[8], str(target_mean), str(beta), str(lam))
    cached = _QONE_MASS_CACHE.get(key)
    if cached is not None:
        return cached
    value = qone_shift_lower(face, target_mean, beta, lam)
    if len(_QONE_MASS_CACHE) >= 100000:
        _QONE_MASS_CACHE.clear()
    _QONE_MASS_CACHE[key] = value
    return value


def certified_qone_bound(box: Box, parameters) -> tuple[arb, arb, str]:
    """Best rigorous bound currently available on the feasible part of box."""
    ordinary = complement_box_bound(box, parameters, LAMBDAS)
    gap_lower = ordinary.gap_lower
    objective_lower = ordinary.objective_lower
    method = ordinary.method
    if min(box[6][1], box[7][1]) <= 1.0 / 16.0:
        for lam in LAMBDAS[:5]:
            shifted, shifted_method = qone_feasible_shift_lower(
                box, parameters.mean, parameters.beta, lam)
            if shifted > gap_lower:
                gap_lower = shifted
                method = shifted_method
    elif (
        box[8][1] == 1.0
        and _endpoint_log_domain(box[8][0])
        and max(box[6][1], box[7][1]) <= 0.5
    ):
        endpoint_lower, endpoint_method = best_qone_endpoint_lower(
            box, parameters.mean, parameters.beta)
        if endpoint_lower > gap_lower:
            gap_lower = endpoint_lower
            method = endpoint_method
    if zero_mass_total_upper(box).upper() <= arb(1) / 16:
        face = zero_mass_face(box)
        for lam in LAMBDAS:
            face_lower, face_method = _cached_qone_mass_face(
                face, parameters.mean, parameters.beta, lam)
            mass_lower = face_lower - zero_mass_shift_penalty(box, lam)
            if mass_lower > gap_lower:
                gap_lower = mass_lower
                method = f"mass-zero({face_method})"
    if gap_lower >= ZERO:
        objective_lower = max(objective_lower, ONE, key=float)
    return gap_lower, objective_lower, method


def run(
    source: Path,
    budget: int,
    seconds: float,
    checkpoint: Optional[Path],
    slice_index: int = 0,
    slice_count: int = 1,
) -> dict[str, Any]:
    all_initial, source_meta = initial_active_boxes(source)
    if slice_count <= 0 or not 0 <= slice_index < slice_count:
        raise ValueError("slice requires 0 <= index < count")
    initial = all_initial[slice_index::slice_count]
    source_meta.update({
        "slice_index": slice_index,
        "slice_count": slice_count,
        "slice_active_boxes": len(initial),
    })
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    heap: list[tuple[int, float, int, int, Box]] = []
    serial = 0
    initial_keys = set()
    for original in initial:
        box, _ = reduce_simplex_face(original)
        if box in initial_keys:
            continue
        initial_keys.add(box)
        heapq.heappush(heap, (0, -10.0, serial, 0, box))
        serial += 1

    started = time.monotonic()
    processed = split = mean_infeasible = gap_cleared = objective_cleared = 0
    mean_contractions = 0
    # Zeroing a positive support lowers the mean and is not a gauge.  The
    # boundary-layer gap inequality alone therefore does not authorize the old
    # q-face feasibility discard after that reduction.
    zero_layer_reductions = 0
    simplex_face_reductions = sum(
        reduce_simplex_face(box)[1] for box in initial)
    method_counts: dict[str, int] = {}
    maximum_depth = 0

    while heap and processed < budget and time.monotonic() - started < seconds:
        _, _, _, depth, box = heapq.heappop(heap)
        processed += 1
        maximum_depth = max(maximum_depth, depth)
        box, contractions = contract_qone_mean(box, parameters.mean)
        mean_contractions += contractions
        if box is None or not feasible(box, parameters.mean):
            mean_infeasible += 1
            continue
        gap_lower, objective_lower, method = certified_qone_bound(
            box, parameters)
        method_counts[method] = method_counts.get(method, 0) + 1
        if gap_lower >= ZERO:
            gap_cleared += 1
            continue
        if objective_lower >= ONE:
            objective_cleared += 1
            continue
        children = bisect_active(box)
        split += 1
        priority = float(gap_lower)
        for child in children:
            child, face = reduce_simplex_face(child)
            simplex_face_reductions += int(face)
            heapq.heappush(
                heap, (-(depth + 1), priority, serial, depth + 1, child))
            serial += 1

        if checkpoint and processed % 10000 == 0:
            payload = {
                "processed": processed,
                "frontier": len(heap),
                "split": split,
                "mean_infeasible": mean_infeasible,
                "gap_cleared": gap_cleared,
                "objective_cleared": objective_cleared,
                "maximum_depth": maximum_depth,
                "mean_contractions": mean_contractions,
                "elapsed_seconds": time.monotonic() - started,
            }
            checkpoint.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")

    # An unfinished heap is already a conservative residual cover.  Rebinding
    # every pending box would duplicate the expensive endpoint proof and cannot
    # strengthen a NUMERICAL result; a completed proof has an empty heap.
    residual = [{
        "active_box": active_compact(box),
        "depth": depth,
        "maximum_width": max(active_widths(box)),
    } for _, _, _, depth, box in sorted(heap)]

    elapsed = time.monotonic() - started
    stop = "complete" if not residual else (
        "box_budget" if processed >= budget else "time_cap")
    report: dict[str, Any] = {
        "tool": "liu9_qone_face.py",
        "claim_status": "PROVED" if not residual else "NUMERICAL",
        "source": source_meta,
        "run": {
            "processed": processed,
            "split": split,
            "mean_infeasible": mean_infeasible,
            "gap_cleared": gap_cleared,
            "objective_cleared": objective_cleared,
            "residual_count": len(residual),
            "maximum_depth": maximum_depth,
            "elapsed_seconds": elapsed,
            "stop_reason": stop,
            "method_counts": method_counts,
            "mean_contractions": mean_contractions,
            "zero_layer_reductions": zero_layer_reductions,
            "simplex_face_reductions": simplex_face_reductions,
        },
        "residual": residual,
        "scope": (
            "exact q=1 face; inactive P0 gauges removed; exact mean-deficit "
            "contractor and endpoint-factored support bound composed"
        ),
        "component_reports": {},
        "limitations": (
            [] if not residual else [
                "An unfinished q=1 frontier is not a proof of the face.",
                "Finite intervals q<1 require a separate inward-q remainder bound.",
            ]),
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    parser.add_argument("--budget", type=int, default=200001)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--slice-index", type=int, default=0)
    parser.add_argument("--slice-count", type=int, default=1)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = run(
        args.source, args.budget, args.seconds, args.checkpoint,
        args.slice_index, args.slice_count)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    run_data = report["run"]
    print("LIU H2 Q=1 ACTIVE-FACE CERTIFIER")
    print("MACHINE VERIFIED [source]: %d q-at-one boxes -> %d unique active boxes "
          "(%.3fx duplicate)." % (
              report["source"]["q_at_one_source_boxes"],
              report["source"]["unique_active_boxes"],
              report["source"]["duplication_factor"]))
    print("%s [run]: processed=%d, gap-cleared=%d, objective-cleared=%d, "
          "mean-infeasible=%d, residual=%d, depth=%d, stop=%s, elapsed=%.2fs."
          % (report["claim_status"], run_data["processed"],
             run_data["gap_cleared"], run_data["objective_cleared"],
             run_data["mean_infeasible"], run_data["residual_count"],
             run_data["maximum_depth"], run_data["stop_reason"],
             run_data["elapsed_seconds"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if not report["residual"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
