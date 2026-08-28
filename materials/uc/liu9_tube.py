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
                "raw h'''([0,rho]) is non-finite, so no one-piece cubic bound "
                "exists; the y*log(1/y) layer is now quantitative "
                "(liu9_boundary_layer.py), but endpoint D(Q)>=0 still has no "
                "quantitative uniform q-remainder constant"
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
    print("   UNITS: that ratio is (Phi-1)/dist^2.  Since gap = EHX*(Phi-1) exactly and")
    print("   EHX(P*)=p*h(x)=0.5526667300..., the same modes cap the raw-gap ratio")
    print("   gap/dist^2 at 0.38712508776..., which is the number the boundary-layer")
    print("   and mirror-layer estimates have to beat.")
    print("PROVED [liu9_endpoint.py, prior certificate]: D(Q)>=0 at q=0 and q=1,")
    print("   with equality exactly at Q=P* modulo its coordinate gauges.  This handles")
    print("   the sign of the endpoint first variation, but does not provide a uniform")
    print("   positive lower coefficient or a q-uniform higher-order remainder.")
    print("PROVED [Arb obstruction to the requested single cubic Taylor bound]: every")
    print("   raw-coordinate tube contains positive-mass supports tending to zero, where")
    print("   h'''(y)=(1-2y)/(y^2(1-y)^2) is unbounded.  The positive local term is")
    print("   B*y*log(1/y), not a C^3 remainder; B was Arb-certified positive in")
    print("   liu9_binding.py.  That sign information is now quantitative:")
    print("PROVED [liu9_boundary_layer.py, 2026-08-27]: with y0=1/32 and")
    print("   Lam=2(1-beta)(mean-y0)-1, K=2(1-beta)log(1/(1-y0)), every")
    print("   mean-feasible point obeys gap(V)>=gap(V|layer:=0)+m0*sum w_j b_j")
    print("   with m0=Lam*(log(1/y0)+1)-K>0 Arb-certified.  The layer is")
    print("   therefore reduced to the face b_j=0, where |h'''|<=1022.94 is")
    print("   finite for every surviving support in [y0,1-y0].")
    print("PROVED [liu9_mirror_layer.py, 2026-08-28]: that face is only C^3 for")
    print("   supports in [y0,1-y0].  The mirror stratum b_j->1 carries the same")
    print("   singularity with leading coefficient 1-2(1-beta)W1-2*beta*A1, which")
    print("   the mean does not sign; a tube does, via W1 <= rho^2/g(1-t0).  The")
    print("   reduction raises the mirror atoms to b_j=1, so it also preserves")
    print("   mean-feasibility.  Ingredients (ii) and (iv) are supplied; (i) and")
    print("   (iii) are not.")
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
    print("      (the ceiling below is in (Phi-1)/dist^2 units; multiply by 0.5526667300... for raw gap)")
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
    print("   rho supplies both halves: no one-piece finite-C3 local bound exists, and")
    print("   the measured complement runs remain unfinished.")
    print("CONJECTURED [what remains]: prove a piecewise local lemma combining (i) the")
    print("   smooth A/C Hessian modes, (ii) an explicit y*log(1/y) boundary layer,")
    print("   (iii) a q-weighted quantitative consequence of D(Q)>=0, and (iv) the")
    print("   mirror layer b_j->1 that (ii) exposed.")
    print("PROVED [liu9_boundary_layer.py, 2026-08-27]: ingredient (ii).  At y0=1/32 the")
    print("   layer admits kappa up to 0.873 for component-mean deviations up to 0.1.")
    print("PROVED [liu9_mirror_layer.py, 2026-08-28]: ingredient (iv).  The leading")
    print("   coefficient 1-2(1-beta)W1-2*beta*A1 is unsigned in general, but a tube")
    print("   forces W1 <= rho^2/g(1-t0) with g(b)=[b(b-x)]^2, so it is positive below")
    print("   an explicit critical radius; the reduction then raises the mirror atoms")
    print("   to the exact face b_j=1, which also preserves mean-feasibility.  It is")
    print("   certified at (rho,t0)=(1/10,1/64) with kappa 0.393693 and REFUTED at")
    print("   (1/10,1/32), so the threshold is sharp in t0.")
    print("PROVED [units, 2026-08-28]: gap = EHX*(Phi-1) exactly, and EHX(P*)=p*h(x)")
    print("   = 0.5526667300...  The ceiling 0.70046750915... printed in section 2 is")
    print("   in (Phi-1)/dist^2 units; the raw-gap ceiling that (ii) and (iv) must beat")
    print("   is 0.38712508776...  Both do, at their certified parameters.")
    print("PROVED [liu9_survivors.py, 2026-08-27]: every one of the 79 residual")
    print("   complement boxes is covered by the zero-face hypothesis")
    print("   2(1-beta)*mean-1 >= 0.11105875229... > 0 (70 boxes) or is mirror-only and")
    print("   entirely outside the tube (9 boxes).  The measured frontier of section 3")
    print("   is therefore exactly the union of the two boundary strata.")
    print("PROVED [liu9_smooth_chart.py + liu9_chart_centered.py, 2026-08-28]:")
    print("   ingredient (i) on the active-mean chart.  Both gap and dist^2 vanish to")
    print("   second order at the chart centre, so Taylor's MEAN-VALUE form gives")
    print("   f(v)=(1/2)v^T Hess f(xi) v exactly and h''' is never needed -- the very")
    print("   thing that defeated the one-piece cubic bound.  The right comparison is")
    print("   the pencil Hess gap - kappa*Hess dist^2 being PSD, in which the q(1-q)")
    print("   factors cancel, so the ceiling is q-free.  A two-variable interval jet")
    print("   reproduces H*diag(A,C q(1-q)) and diag(2px^2,2(p^2+px^2)q(1-q)) to")
    print("   1.9e-67, independently confirming the section 2 curvatures.  Certified:")
    print("   gap >= kappa*dist^2 on |s-x|,|d| <= 1/256 for every q in [1/4,3/4], with")
    print("   raw-gap kappa >= 4119063/33554432 = 0.1227576.  A centered third-order")
    print("   form cut the enclosure inflation constant from 811.29 to 72.20.")
    print("PROVED [liu9_chart_cover.py, 2026-08-28]: that pencil statement now holds")
    print("   on |s-x|,|d| <= 1/32 for every q in [1/4096,4095/4096], at the same")
    print("   kappa, over 393216 exhaustively abutting cells; weakest determinant")
    print("   margin 2.48e-05 at (s-x,d,q)=(15/512,-15/512,1/4096).  Two premises")
    print("   behind the old plan were wrong.  First, a cover DOES suffice: PSD is a")
    print("   pointwise property, and the box is convex and contains the centre, so")
    print("   the mean-value step gets its segments from geometry once rather than")
    print("   per cell -- no radial quadrature is needed.  Second, the radius was not")
    print("   the binding gap; the q range was.  Ingredient (iii) needs the chart down")
    print("   to q=1/4096, and a uniform q grid stops certifying below q=1/64 because")
    print("   gap.hdd and dist.hdd both carry an exact factor q(1-q) -- the ratio")
    print("   m22/(q(1-q)) is 0.69512040 at q=1/4,1/64,1/128,1/1024,1/4096 alike -- so")
    print("   the true entry vanishes linearly while a fixed-width cell's enclosure")
    print("   error does not.  Octave cells of relative width 1/64 fix it.")
    print("CONDITIONAL [liu9_chart_cover.py]: that PSD statement is Step A.  Turning")
    print("   it into gap >= kappa*dist^2 is Step B, and Step B does NOT give an")
    print("   exact zero: Liu's x is a numerically determined root of (87)-(90), so")
    print("   the centre carries that residual.  The certified form is")
    print("   gap - kappa*dist^2 >= -7.361e-69 on the box, and it is an all-q")
    print("   enclosure rather than a sampled bound -- the centre value and d/ds")
    print("   come from a q-free jet, and d/dd is exactly zero for every q because")
    print("   each term carries a q-free multiplier times (1-q)(-q)+q(1-q), whose")
    print("   coefficients all vanish in exact Fraction arithmetic.  The deficit is")
    print("   inherited unchanged from the radius 1/256 certificate.")
    print("PROVED [liu9_qdegenerate.py, 2026-08-28]: ingredient (iii) on the pure-d")
    print("   endpoint chart.  Inner core |d|<=1/32 has q-uniform kappa 1/3; the")
    print("   endpoint annulus 1/32<=|d|<=1/4 with min(q,1-q)<=1/4096 has kappa 1/20,")
    print("   from H*D(Q_d) >= (1/4)delta(Q_d)^2 over 384 exact dyadic cells.  The")
    print("   seam is PROVED under rho^2 <= p^2 eps_sm^2 q_*(1-q_*), giving")
    print("   rho=1/4096, so no pure-d tube point falls between the two q regimes.")
    print("   Only the simultaneous swap (q,P0,P1)->(1-q,P1,P0) is a symmetry; q->1-q")
    print("   alone is REFUTED, changing the gap by 2.155e-02.")
    print("REFUTED [liu9_ninevar.py, 2026-08-28]: the chart result does NOT extend to")
    print("   an unrestricted ambient nine-variable neighbourhood.  At")
    print("   (y,eps,q)=(1/32,1/1024,1/2) the split (p-eps,eps,q,x,y,0,x,y,0) gives")
    print("   raw gap -5.5338e-04 and pencil -5.5348e-04; exact Fraction arithmetic")
    print("   puts the linear coefficient below -1/2 and the finite pencil below")
    print("   -1/2048.  The leading order is eps with no eps*log(1/eps) term.")
    print("CONDITIONAL: that counterexample is mean-INFEASIBLE (mean-target")
    print("   -6.4408e-04), so it does not refute Liu's Hypothesis 2, which lives in")
    print("   the half-space mean >= p*x.  The exact mean-preserving control at the")
    print("   same y has linear coefficient above 1/64, and over a 65-point y scan")
    print("   the mean-preserving family is positive at every point (0 negative)")
    print("   while the ambient split is negative at 46.  Its weakest coefficient,")
    print("   7.865e-08, sits at y=707/1024 ~ x, where the inserted atom merges with")
    print("   the existing one and the perturbation degenerates.")
    print("OPEN, and now a single named gap.  A feasible-half-space nine-variable")
    print("   theorem is not proved -- only the chart slice is, plus positive")
    print("   evidence on one transverse family.  Until that is closed no tube")
    print("   radius follows, because the tube is a nine-variable neighbourhood and")
    print("   the chart is a four-parameter slice of it.  No tube radius is")
    print("   certified and this script does not overclaim one.")


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
