#!/usr/bin/env python3
"""Feasibility measurement for a tube-plus-complement proof of Liu H2.

This is deliberately a pilot, not a proof trace.  Individual interval box
classifications are rigorous, but an unfinished branch-and-bound and every
cost extrapolation are labelled NUMERICAL or CONJECTURED.  The mathematical
parameters are obtained from Liu's defining equations, never from the rounded
printed decimal for c'.

Run from the repository root with

    math/.venv/bin/python math/uc/liu9_tube.py
"""

from __future__ import annotations

import argparse
import heapq
import math
import os
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Optional, Sequence

# This workstation is shared with live certification workers.  These variables
# must be set before importing modules which may load NumPy/SciPy or a BLAS.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    certify_curvatures,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import (  # noqa: E402
    _semantic_bounds,
    arb_fraction,
    arb_hull,
    evaluate_arb,
)
from liu9_pilot import (  # noqa: E402
    Box,
    GradientUnavailable,
    ONE,
    SPLIT_TIE_ORDER,
    ZERO,
    _box_arbs,
    _gap_mean_gradient,
    _raw_gap_at,
    _simplex_vertices,
    mean_corner_range,
)

ctx.prec = max(ctx.prec, 320)

DEFAULT_RHOS = ("0.1", "0.03", "0.01", "0.003")
DEFAULT_LAMBDAS = ("0", "-0.4", "-0.8", "-0.9", "-1.0", "-1.2")
DEFAULT_BOX_BUDGET = 75_001
DEFAULT_SECONDS_PER_RHO = 100.0
DEFAULT_WALL_SECONDS = 420.0


@dataclass(frozen=True)
class LocalAttempt:
    rho: Fraction
    smooth_kappa: arb
    third_derivative: arb
    tube_certified: bool
    boundary_margin: Optional[arb]
    blocker: str


@dataclass(frozen=True)
class BoxBound:
    gap_lower: arb
    objective_lower: arb
    method: str
    centered_usable: bool

    @property
    def key(self) -> float:
        return float(self.objective_lower)


@dataclass
class ComplementStats:
    rho: Fraction
    processed: int = 0
    objective_cleared: int = 0
    mean_infeasible: int = 0
    tube_excluded: int = 0
    split: int = 0
    residual: int = 0
    fully_outside_evaluated: int = 0
    straddling_evaluated: int = 0
    centered_uses: int = 0
    maximum_depth: int = 0
    weakest_cleared_gap: Optional[arb] = None
    cover_gap_lower: float = -math.inf
    elapsed: float = 0.0
    stop_reason: str = "not_started"
    clearing_rate: float = 0.0
    resolved_leaf_rate: float = 0.0
    extrapolated_boxes: float = math.inf
    extrapolated_days_28: float = math.inf
    throughput: float = 0.0


# ---------------------------------------------------------------------------
# Equality quotient and a gauge-invariant distance.


def distance_squared_mp(values: Sequence[mpmath.mpf], parameters: MPParameters) -> mpmath.mpf:
    """The exact symbolic distance formula evaluated with mpmath scalars."""
    if len(values) != 9:
        raise ValueError("expected nine coordinates")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    a3 = 1 - a1 - a2
    masses = (a1, a2, a3)

    def component(points: Sequence[mpmath.mpf]) -> mpmath.mpf:
        mean = sum(mass * point for mass, point in zip(masses, points))
        support_defect = sum(
            mass * (point * (point - parameters.x)) ** 2
            for mass, point in zip(masses, points)
        )
        return (mean - parameters.mean) ** 2 + support_defect

    return (1 - q) * component((b0, b2, b4)) + q * component((b1, b3, b5))


def _arb_interval(lower: float, upper: float) -> arb:
    answer = arb(lower).union(arb(upper))
    if not answer.is_finite():
        raise ArithmeticError("non-finite interval construction")
    return answer


def _square_range(value: arb) -> arb:
    """Exact range of z^2 over one real Arb interval."""
    lower, upper = value.lower(), value.upper()
    upper_square = max((lower * lower, upper * upper), key=float)
    if value.contains(0):
        lower_square = arb(0)
    else:
        lower_square = min((lower * lower, upper * upper), key=float)
    answer = arb_hull(lower_square, upper_square)
    if not answer.is_finite():
        raise ArithmeticError("invalid square enclosure")
    return answer


def _affine_simplex_range(box: Box, coefficients: Sequence[arb]) -> arb:
    """Range of sum a_i*c_i with independent c_i intervals and clipped a box."""
    vertices = _simplex_vertices(box)
    if not vertices:
        raise ValueError("box misses the mass simplex")
    bounds = tuple(_semantic_bounds(coefficient, arb(0), None)
                   for coefficient in coefficients)
    lower_values: list[arb] = []
    upper_values: list[arb] = []
    for a1, a2 in vertices:
        masses = (arb(a1), arb(a2), arb(1.0 - a1 - a2))
        lower_values.append(sum(
            (mass * bound[0] for mass, bound in zip(masses, bounds)),
            arb(0),
        ))
        upper_values.append(sum(
            (mass * bound[1] for mass, bound in zip(masses, bounds)),
            arb(0),
        ))
    lower = min((value.lower() for value in lower_values), key=float)
    upper = max((value.upper() for value in upper_values), key=float)
    return arb_hull(lower, upper)


def _component_distance_range(
    box: Box,
    support_indices: Sequence[int],
    x: arb,
    target_mean: arb,
) -> arb:
    points = tuple(_arb_interval(*box[index]) for index in support_indices)
    mean_range = _affine_simplex_range(box, points)
    mean_defect = _square_range(mean_range - target_mean)
    support_coefficients = tuple(
        _square_range(point * (point - x)) for point in points
    )
    support_defect = _affine_simplex_range(box, support_coefficients)
    answer = mean_defect + support_defect
    lower, upper = _semantic_bounds(answer, arb(0), None)
    return arb_hull(lower, upper)


def distance_squared_box(box: Box, parameters: ArbParameters) -> arb:
    """Certified enclosure of dist^2 on a box intersected with the simplex."""
    d0 = _component_distance_range(box, (3, 4, 5), parameters.x, parameters.mean)
    d1 = _component_distance_range(box, (6, 7, 8), parameters.x, parameters.mean)
    q_lower, q_upper = (arb(box[2][0]), arb(box[2][1]))
    lower_candidates = (
        (1 - q_lower) * d0.lower() + q_lower * d1.lower(),
        (1 - q_upper) * d0.lower() + q_upper * d1.lower(),
    )
    upper_candidates = (
        (1 - q_lower) * d0.upper() + q_lower * d1.upper(),
        (1 - q_upper) * d0.upper() + q_upper * d1.upper(),
    )
    lower = min(lower_candidates, key=float).lower()
    upper = max(upper_candidates, key=float).upper()
    answer = arb_hull(lower, upper)
    if not answer.is_finite():
        raise ArithmeticError("invalid distance enclosure")
    semantic_lower, semantic_upper = _semantic_bounds(answer, arb(0), None)
    return arb_hull(semantic_lower, semantic_upper)


def canonical_values(
    parameters: MPParameters, q: mpmath.mpf, zero_split: mpmath.mpf
) -> tuple[mpmath.mpf, ...]:
    """A canonical section of the repeated-zero-atom representation gauge."""
    return (
        parameters.p,
        zero_split,
        q,
        parameters.x,
        mpmath.mpf(0),
        mpmath.mpf(0),
        parameters.x,
        mpmath.mpf(0),
        mpmath.mpf(0),
    )


def verify_distance_samples(mp_parameters: MPParameters, arb_parameters: ArbParameters) -> None:
    tolerance = mpmath.mpf("1e-70")
    samples: list[tuple[str, tuple[mpmath.mpf, ...], bool]] = []
    zero_mass = 1 - mp_parameters.p
    for q_text, split_fraction in (("0", "0.2"), ("0.5", "0.5"), ("1", "0.8")):
        q = mpmath.mpf(q_text)
        split = zero_mass * mpmath.mpf(split_fraction)
        samples.append((f"canonical q={q_text}", canonical_values(mp_parameters, q, split), True))

    # At q=0 the inactive law is a zero-mixture-weight gauge and must disappear.
    endpoint = list(canonical_values(mp_parameters, mpmath.mpf(0), zero_mass / 2))
    endpoint[6:] = (mpmath.mpf("0.17"), mpmath.mpf("0.63"), mpmath.mpf("0.91"))
    samples.append(("q=0 arbitrary inactive P1", tuple(endpoint), True))

    paired = list(canonical_values(mp_parameters, mpmath.mpf("0.5"), zero_mass / 2))
    paired[3] += mpmath.mpf("0.03")
    paired[6] -= mpmath.mpf("0.03")
    samples.append(("paired support displacement", tuple(paired), False))

    generic = (
        mpmath.mpf("0.65"), mpmath.mpf("0.2"), mpmath.mpf("0.4"),
        mpmath.mpf("0.9"), mpmath.mpf("0.1"), mpmath.mpf("0.3"),
        mpmath.mpf("0.8"), mpmath.mpf("0.2"), mpmath.mpf("0.4"),
    )
    samples.append(("generic off-manifold point", generic, False))

    for name, values, should_vanish in samples:
        squared = distance_squared_mp(values, mp_parameters)
        if should_vanish and abs(squared) > tolerance:
            raise AssertionError(f"distance failed to vanish at {name}: {squared}")
        if not should_vanish and not squared > 0:
            raise AssertionError(f"distance failed to separate {name}: {squared}")
        print(
            "NUMERICAL [mpmath distance sample]: %s; dist^2=%s; expected=%s"
            % (
                name,
                mpmath.nstr(squared, 12),
                "zero in quotient" if should_vanish else "positive",
            )
        )

    # A separated exact-dyadic point also gets a rigorous Arb check.
    off_box: Box = (
        (0.625, 0.625), (0.25, 0.25), (0.5, 0.5),
        (0.875, 0.875), (0.125, 0.125), (0.25, 0.25),
        (0.75, 0.75), (0.25, 0.25), (0.5, 0.5),
    )
    off_enclosure = distance_squared_box(off_box, arb_parameters)
    if not off_enclosure > 0:
        raise AssertionError("Arb did not separate the off-manifold sample")
    print(f"PROVED [Arb sample separation]: off-manifold dist^2 in {off_enclosure} > 0.")


# ---------------------------------------------------------------------------
# Local expansion measurement.


def _nonfinite_third_entropy_enclosure(rho: Fraction) -> arb:
    """Attempt h''' on the zero-support part of the requested raw tube."""
    support = arb(0).union(arb_fraction(rho))
    denominator = support**2 * (1 - support) ** 2
    return (1 - 2 * support) / denominator


def measure_local_attempts(
    rhos: Sequence[Fraction], parameters: ArbParameters
) -> tuple[LocalAttempt, ...]:
    curvatures = certify_curvatures(parameters)

    # On the active-mean mode, dist^2=p*x^2*ds^2+o(ds^2).  On the paired
    # component mode, dist^2=q(1-q)*(p^2+p*x^2)*dd^2+o(dd^2).
    # Since the objective expansion uses one half of the Hessian, these are
    # the two exact quadratic comparison ratios for the displayed smooth modes.
    face_ratio = curvatures.face / (2 * parameters.p * parameters.x**2)
    split_ratio = curvatures.split_unit / (
        2 * (parameters.p**2 + parameters.p * parameters.x**2)
    )
    ratio_lower = min(face_ratio.lower(), split_ratio.lower(), key=float)
    ratio_upper = min(face_ratio.upper(), split_ratio.upper(), key=float)
    smooth_kappa = arb_hull(ratio_lower, ratio_upper)
    if not smooth_kappa > arb_fraction(Fraction(3, 5)):
        raise AssertionError("smooth quadratic comparison did not certify kappa>3/5")

    attempts: list[LocalAttempt] = []
    for rho in rhos:
        third = _nonfinite_third_entropy_enclosure(rho)
        if third.is_finite():
            raise AssertionError("h''' enclosure touching zero unexpectedly finite")
        attempts.append(LocalAttempt(
            rho=rho,
            smooth_kappa=smooth_kappa,
            third_derivative=third,
            tube_certified=False,
            boundary_margin=None,
            blocker=(
                "raw h'''([0,rho]) is non-finite; the proved y*log(1/y) barrier "
                "and endpoint D(Q)>=0 have no quantitative uniform q-remainder constant"
            ),
        ))
    return tuple(attempts)


# ---------------------------------------------------------------------------
# Counterfactual complement branch-and-bound.


def _centered_shift_enclosure(
    box: Box,
    target_mean: arb,
    beta: arb,
    lam: arb,
    gradients: Optional[tuple[tuple[arb, ...], tuple[arb, ...]]] = None,
) -> arb:
    if gradients is None:
        gradients = _gap_mean_gradient(box, beta)
    gap_gradient, mean_gradient = gradients
    centers = tuple((lower + upper) / 2.0 for lower, upper in box)
    center_arbs = tuple(arb(value) for value in centers)
    gap_center, mean_center = _raw_gap_at(center_arbs, beta)
    answer = gap_center + lam * (mean_center - target_mean)
    for index, ((lower, upper), center) in enumerate(zip(box, centers)):
        displacement = arb(lower - center).union(arb(upper - center))
        answer += (gap_gradient[index] + lam * mean_gradient[index]) * displacement
    if not answer.is_finite():
        raise GradientUnavailable("non-finite centered enclosure")
    return answer


def _best_lower(candidates: Iterable[tuple[arb, str]]) -> tuple[arb, str]:
    choices = list(candidates)
    if not choices:
        raise ValueError("no bound candidates")
    return max(choices, key=lambda item: float(item[0]))


def complement_box_bound(
    box: Box,
    parameters: ArbParameters,
    lambdas: Sequence[Fraction],
) -> BoxBound:
    terms = evaluate_arb(_box_arbs(box), parameters.beta, monotone_corners=True)
    direct_gap = terms.numerator - terms.ehx
    candidates: list[tuple[arb, str]] = [(direct_gap.lower(), "direct")]
    centered_usable = False

    if (
        max(upper - lower for lower, upper in box[3:]) <= 0.25
        and all(lower > 0.0 and upper < 1.0 for lower, upper in box[3:])
    ):
        try:
            gradients = _gap_mean_gradient(box, parameters.beta)
            centered_usable = True
            for lam_fraction in lambdas:
                enclosure = _centered_shift_enclosure(
                    box,
                    parameters.mean,
                    parameters.beta,
                    arb_fraction(lam_fraction),
                    gradients,
                )
                candidates.append((enclosure.lower(), f"center(lambda={lam_fraction})"))
        except GradientUnavailable:
            centered_usable = False

    gap_lower, method = _best_lower(candidates)
    numerator_lower, _ = _semantic_bounds(terms.numerator, ZERO, None)
    denominator_lower, denominator_upper = _semantic_bounds(terms.ehx, ZERO, None)
    if denominator_upper > ZERO:
        direct_ratio = numerator_lower / denominator_upper
    else:
        direct_ratio = ZERO
    objective_candidates = [direct_ratio.lower()]
    if gap_lower >= ZERO:
        objective_candidates.append(ONE)
    elif denominator_lower > ZERO:
        objective_candidates.append((ONE + gap_lower / denominator_lower).lower())
    objective_lower = max(objective_candidates, key=float)
    if not gap_lower.is_finite() or not objective_lower.is_finite():
        raise ArithmeticError("non-finite objective box bound")
    return BoxBound(gap_lower, objective_lower, method, centered_usable)


def _widest_coordinate(box: Box) -> int:
    widths = tuple(upper - lower for lower, upper in box)
    widest = max(widths)
    return next(index for index in SPLIT_TIE_ORDER if widths[index] == widest)


def _bisect(box: Box) -> tuple[Box, Box]:
    coordinate = _widest_coordinate(box)
    lower, upper = box[coordinate]
    midpoint = (lower + upper) / 2.0
    left, right = list(box), list(box)
    left[coordinate] = (lower, midpoint)
    right[coordinate] = (midpoint, upper)
    return tuple(left), tuple(right)


def run_complement(
    rho: Fraction,
    parameters: ArbParameters,
    lambdas: Sequence[Fraction],
    box_budget: int,
    seconds: float,
    global_deadline: float,
) -> ComplementStats:
    stats = ComplementStats(rho=rho)
    rho_squared = arb_fraction(rho * rho)
    root: Box = tuple((0.0, 1.0) for _ in range(9))
    heap: list[tuple[int, float, float, int, int, Box, BoxBound, arb]] = []
    serial = 0
    started = time.monotonic()
    deadline = min(started + seconds, global_deadline)

    def classify(box: Box, depth: int) -> None:
        nonlocal serial
        stats.processed += 1
        stats.maximum_depth = max(stats.maximum_depth, depth)
        mean_range = mean_corner_range(box)
        if mean_range is None or mean_range.upper() < parameters.mean.lower():
            stats.mean_infeasible += 1
            return

        distance_range = distance_squared_box(box, parameters)
        if distance_range.upper() < rho_squared:
            stats.tube_excluded += 1
            return

        straddles_tube = not distance_range.lower() >= rho_squared
        if straddles_tube:
            stats.straddling_evaluated += 1
        else:
            stats.fully_outside_evaluated += 1

        bound = complement_box_bound(box, parameters, lambdas)
        stats.centered_uses += int(bound.centered_usable)
        if bound.gap_lower >= ZERO or bound.objective_lower >= ONE:
            stats.objective_cleared += 1
            if (
                stats.weakest_cleared_gap is None
                or bound.gap_lower < stats.weakest_cleared_gap
            ):
                stats.weakest_cleared_gap = bound.gap_lower
            return

        serial += 1
        # The hard boxes are those near equality.  Resolve possible tube
        # membership before spending the budget refining boxes already proved
        # wholly outside; among straddlers, take the one closest to fitting.
        geometry_priority = (
            float(distance_range.upper() - rho_squared)
            if straddles_tube else bound.key
        )
        heapq.heappush(
            heap,
            (0 if straddles_tube else 1, geometry_priority, bound.key,
             serial, depth, box, bound, distance_range),
        )

    classify(root, 0)
    stats.stop_reason = "complete"
    while heap:
        if stats.processed + 2 > box_budget:
            stats.stop_reason = "box_budget"
            break
        if time.monotonic() >= deadline:
            stats.stop_reason = (
                "global_wall_clock" if deadline == global_deadline else "rho_wall_clock"
            )
            break
        _, _, _, _, depth, box, _, _ = heapq.heappop(heap)
        stats.split += 1
        left, right = _bisect(box)
        classify(left, depth + 1)
        classify(right, depth + 1)

    stats.elapsed = time.monotonic() - started
    stats.residual = len(heap)
    stats.throughput = stats.processed / max(stats.elapsed, 1e-12)
    if heap:
        stats.cover_gap_lower = min(float(item[6].gap_lower) for item in heap)
    else:
        stats.cover_gap_lower = 0.0

    # The clearing rate deliberately excludes mean-infeasible and inside-tube
    # leaves: it asks how often the objective itself cleared a surviving
    # complement leaf.  Treating every residual as genuinely feasible is an
    # optimistic convention for extrapolation.
    objective_denominator = stats.objective_cleared + stats.residual
    stats.clearing_rate = (
        stats.objective_cleared / objective_denominator
        if objective_denominator else 1.0
    )
    leaves = (
        stats.objective_cleared
        + stats.mean_infeasible
        + stats.tube_excluded
        + stats.residual
    )
    resolved = leaves - stats.residual
    stats.resolved_leaf_rate = resolved / leaves if leaves else 1.0
    if stats.residual == 0:
        stats.extrapolated_boxes = float(stats.processed)
    elif stats.clearing_rate > 0:
        stats.extrapolated_boxes = stats.processed / stats.clearing_rate
    else:
        stats.extrapolated_boxes = math.inf
    if math.isfinite(stats.extrapolated_boxes) and stats.throughput > 0:
        stats.extrapolated_days_28 = (
            stats.extrapolated_boxes / (28 * stats.throughput * 86_400)
        )
    return stats


# ---------------------------------------------------------------------------
# Reporting.


def _fraction_text(value: Fraction) -> str:
    return format(float(value), ".6g")


def _count_text(value: float) -> str:
    return "infinity" if not math.isfinite(value) else f"{value:.3g}"


def print_geometry(
    mp_parameters: MPParameters,
    arb_parameters: ArbParameters,
    dps: int,
) -> None:
    print("1. EQUALITY QUOTIENT AND DISTANCE")
    print(f"PROVED [320-bit Arb equation bracket]: x in {arb_parameters.x}")
    print(f"PROVED [320-bit Arb equation propagation]: p in {arb_parameters.p}")
    print(f"PROVED [320-bit Arb equation propagation]: m=p*x=1-c' in {arb_parameters.mean}")
    print("PROVED [explicit canonical section]: for q in [0,1] and r in [0,1-p],")
    print("   (a1,a2,q,b0,b2,b4,b1,b3,b5)=(p,r,q,x,0,0,x,0,0).")
    print("   Changing r only splits the zero atom.  Simultaneous atom permutations,")
    print("   coincident-support splits, and arbitrary zero-mass supports are quotient")
    print("   gauges.  At q=0 (respectively q=1), all coordinates of inactive P1")
    print("   (respectively P0) are additionally zero-mixture-weight gauges.")
    print("PROVED [definition]: for P with masses a_i and supports y_i, put")
    print("   delta(P)^2=(sum_i a_i*y_i-m)^2 + sum_i a_i*[y_i(y_i-x)]^2,")
    print("   dist^2=(1-q)*delta(P0)^2+q*delta(P1)^2, and dist=sqrt(dist^2).")
    print("PROVED [nonnegative-summand algebra]: delta(P)=0 iff every positive-mass")
    print("   support lies in {0,x} and its total mass at x is m/x=p, i.e. P=P*.")
    print("   Therefore dist=0 exactly on the diagonal equality family after the")
    print("   stated representation and endpoint inactive-law gauges; it is positive")
    print("   off that quotient.  This is invariant under atom permutations.")
    verify_distance_samples(mp_parameters, arb_parameters)
    print(f"NUMERICAL [mpmath precision]: distance samples used {dps} decimal digits.")
    print()


def print_local_analysis(
    attempts: Sequence[LocalAttempt],
    parameters: ArbParameters,
) -> None:
    curvatures = certify_curvatures(parameters)
    print("2. LOCAL TUBE ATTEMPT")
    print(f"PROVED [Arb, imported binding calculation]: A in {curvatures.face}")
    print(f"PROVED [Arb, imported binding calculation]: C in {curvatures.split_unit}")
    print("PROVED [exact reduced coordinates]: the Hessian is")
    print("   diag(0,A,C*q*(1-q),0) in (q,s,d,r); the q and r nulls are gauges.")
    print("PROVED [distance comparison on the two smooth displayed modes]:")
    print("   face ratio=A/(2*p*x^2), split ratio=C/[2*(p^2+p*x^2)], so")
    print(f"   their common quadratic kappa ceiling is in {attempts[0].smooth_kappa} > 3/5.")
    print("PROVED [liu9_endpoint.py, prior certificate]: D(Q)>=0 at q=0 and q=1,")
    print("   with equality exactly at Q=P* modulo its coordinate gauges.  This handles")
    print("   the sign of the endpoint first variation, but does not provide a uniform")
    print("   positive lower coefficient or a q-uniform higher-order remainder.")
    print("PROVED [Arb obstruction to the requested single cubic Taylor bound]: every")
    print("   raw-coordinate tube contains positive-mass supports tending to zero, where")
    print("   h'''(y)=(1-2y)/(y^2(1-y)^2) is unbounded.  The positive local term is")
    print("   B*y*log(1/y), not a C^3 remainder; B was Arb-certified positive in")
    print("   liu9_binding.py.  A separate quantitative boundary-layer lemma would be")
    print("   required to turn that sign information into kappa*dist^2-C3*dist^3.")
    for attempt in attempts:
        print(
            "PROVED [Arb tube attempt, rho=%s]: h'''([0,rho])=%s; finite=%s; "
            "tube inequality certified=NO; guaranteed boundary margin=unavailable."
            % (
                _fraction_text(attempt.rho),
                attempt.third_derivative,
                attempt.third_derivative.is_finite(),
            )
        )
    print("CONJECTURED [local-proof diagnosis]: the inequality itself may still be true;")
    print("   what fails is the proposed one-piece finite-C3 certification from the")
    print("   currently proved Hessian, insertion-potential, and endpoint-sign data.")
    print()


def print_complement_results(
    results: Sequence[ComplementStats],
    box_budget: int,
    seconds_per_rho: float,
) -> None:
    print("3. COUNTERFACTUAL COMPLEMENT COST")
    print("NUMERICAL [measurement scope]: no rho passed the local tube step, so these")
    print("   runs measure the complement filter counterfactually; they are not pieces")
    print("   of a complete certificate.  Every individual mean/tube/objective discard")
    print("   is nevertheless certified by exact dyadic geometry plus Arb.")
    print(
        "NUMERICAL [fixed budget]: at most %d evaluated boxes or %.1f seconds per rho, one core."
        % (box_budget, seconds_per_rho)
    )
    print("NUMERICAL [rate definition]: clear-rate=objective-cleared/(objective-cleared+residual);")
    print("   mean-infeasible and wholly-inside-tube leaves are reported separately.")
    print("NUMERICAL [complement table]: rho | processed | obj-clear | mean-infeasible | tube-out | residual | split | max-depth | clear-rate | resolved-leaf-rate | boxes/s | optimistic-total | ideal-28-core-days | cover-gap-lb | seconds | stop")
    for stats in results:
        print(
            "NUMERICAL [measured]: %s | %d | %d | %d | %d | %d | %d | %d | %.6g | %.6g | %.3g | %s | %s | %+.6g | %.3f | %s"
            % (
                _fraction_text(stats.rho),
                stats.processed,
                stats.objective_cleared,
                stats.mean_infeasible,
                stats.tube_excluded,
                stats.residual,
                stats.split,
                stats.maximum_depth,
                stats.clearing_rate,
                stats.resolved_leaf_rate,
                stats.throughput,
                _count_text(stats.extrapolated_boxes),
                _count_text(stats.extrapolated_days_28),
                stats.cover_gap_lower,
                stats.elapsed,
                stats.stop_reason,
            )
        )
        print(
            "NUMERICAL [rho=%s branch detail]: fully-outside evaluated=%d; tube-straddling evaluated=%d; centered-MVT uses=%d; weakest cleared gap=%s."
            % (
                _fraction_text(stats.rho),
                stats.fully_outside_evaluated,
                stats.straddling_evaluated,
                stats.centered_uses,
                "none" if stats.weakest_cleared_gap is None else str(stats.weakest_cleared_gap),
            )
        )
    print("CONJECTURED [extrapolation convention]: optimistic-total=processed/clear-rate")
    print("   and ideal days assume the measured one-core throughput, perfect 28-core")
    print("   scaling, and a stationary clearing rate.  A growing residual frontier makes")
    print("   this an optimistic feasibility indicator, not a certified upper bound.")
    print()


def print_verdict(
    attempts: Sequence[LocalAttempt],
    results: Sequence[ComplementStats],
) -> None:
    by_rho = {result.rho: result for result in results}
    print("4. RHO TRADE-OFF AND VERDICT")
    print("rho | certified kappa | tube inequality | boundary margin | complement clear-rate | extrapolated boxes | feasible on 28 cores in days")
    for attempt in attempts:
        result = by_rho[attempt.rho]
        print(
            "%s | -- (smooth-mode ceiling %s) | NO | -- | %.6g | %s | NO"
            % (
                _fraction_text(attempt.rho),
                attempt.smooth_kappa,
                result.clearing_rate,
                _count_text(result.extrapolated_boxes),
            )
        )

    best = max(
        results,
        key=lambda item: (item.clearing_rate, -item.extrapolated_boxes),
    )
    print(
        "NUMERICAL [best measured complement radius]: rho=%s, clear-rate=%.6g, optimistic total=%s boxes, idealized 28-core time=%s days."
        % (
            _fraction_text(best.rho),
            best.clearing_rate,
            _count_text(best.extrapolated_boxes),
            _count_text(best.extrapolated_days_28),
        )
    )
    print("CONJECTURED [plain feasibility verdict]: INFEASIBLE AS SPECIFIED.  No tested")
    print("   rho supplies both halves: every local finite-C3 attempt is blocked by the")
    print("   zero-support entropy singularity and the absence of a quantitative uniform")
    print("   endpoint remainder, while the measured complement runs remain unfinished.")
    print("CONJECTURED [what remains]: prove a piecewise local lemma combining (i) the")
    print("   smooth A/C Hessian modes, (ii) an explicit y*log(1/y) boundary layer, and")
    print("   (iii) a q-weighted quantitative consequence of D(Q)>=0.  Only then can the")
    print("   best complement rho be promoted from a counterfactual timing point to a")
    print("   certificate radius; the current script does not overclaim that theorem.")


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rhos", nargs="+", default=DEFAULT_RHOS)
    parser.add_argument("--box-budget", type=int, default=DEFAULT_BOX_BUDGET)
    parser.add_argument("--seconds-per-rho", type=float, default=DEFAULT_SECONDS_PER_RHO)
    parser.add_argument("--wall-seconds", type=float, default=DEFAULT_WALL_SECONDS)
    parser.add_argument("--prec", type=int, default=320)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--lambdas", nargs="+", default=DEFAULT_LAMBDAS)
    parser.add_argument("--skip-complement", action="store_true")
    args = parser.parse_args(argv)

    if args.box_budget < 101:
        parser.error("--box-budget must be at least 101")
    if args.seconds_per_rho <= 0 or not 1 <= args.wall_seconds <= 840:
        parser.error("positive per-rho time and --wall-seconds in [1,840] are required")
    if args.prec < 320 or args.dps < 70:
        parser.error("--prec >= 320 and --dps >= 70 are required")

    rhos = tuple(Fraction(value) for value in args.rhos)
    if not rhos or any(not 0 < rho < 1 for rho in rhos):
        parser.error("every rho must lie strictly between zero and one")
    lambdas = tuple(Fraction(value) for value in args.lambdas)
    if any(value > 0 for value in lambdas):
        parser.error("all feasibility-shift lambdas must be <= 0")

    ctx.prec = args.prec
    mpmath.mp.dps = max(mpmath.mp.dps, args.dps + 30)
    started = time.monotonic()
    global_deadline = started + args.wall_seconds

    print("LIU HYPOTHESIS-2 TUBE-PLUS-COMPLEMENT FEASIBILITY PILOT")
    print("Status labels are PROVED, NUMERICAL, or CONJECTURED.")
    print("PROVED [runtime configuration]: BLAS/OpenMP use one core; no background jobs;")
    print(
        "   Arb precision=%d bits; per-rho budget=%d; per-rho cap=%.1fs; global cap=%.1fs."
        % (args.prec, args.box_budget, args.seconds_per_rho, args.wall_seconds)
    )
    print()

    mp_parameters = solve_equation_parameters(args.dps + 20)
    arb_parameters = certify_equation_parameters(mp_parameters)
    print_geometry(mp_parameters, arb_parameters, args.dps)
    attempts = measure_local_attempts(rhos, arb_parameters)
    print_local_analysis(attempts, arb_parameters)

    if args.skip_complement:
        print("NUMERICAL [explicit CLI choice]: complement measurement skipped.")
        return

    results: list[ComplementStats] = []
    for rho in rhos:
        if time.monotonic() >= global_deadline:
            stats = ComplementStats(rho=rho)
            stats.stop_reason = "global_wall_clock_before_start"
            results.append(stats)
            continue
        results.append(run_complement(
            rho,
            arb_parameters,
            lambdas,
            args.box_budget,
            args.seconds_per_rho,
            global_deadline,
        ))
    print_complement_results(results, args.box_budget, args.seconds_per_rho)
    print_verdict(attempts, results)


if __name__ == "__main__":
    main()
