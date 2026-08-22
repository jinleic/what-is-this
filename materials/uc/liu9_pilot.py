#!/usr/bin/env python3
"""Feasibility pilot for interval branch-and-bound on Liu's 9-D objective.

This is deliberately not a proof-trace producer.  It measures how a direct
nine-dimensional Arb cover behaves at rational targets below Liu's reported
constant.  Every box decision is interval-certified, but an unfinished run is
only a feasibility experiment, not a global certificate.
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
from typing import Any, Iterable, Optional, Sequence

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_objective import (  # noqa: E402
    LIU_BETA_DECIMAL,
    LIU_C_DECIMAL,
    _ArbOps,
    _ScalarOps,
    _formula,
    _semantic_bounds,
    arb_fraction,
    arb_hull,
    as_fraction,
    evaluate_arb,
    h_arb,
    mean_value,
    stationary_parameters,
)

ctx.prec = max(ctx.prec, 128)

Box = tuple[tuple[float, float], ...]
ZERO = arb(0)
ONE = arb(1)
C_STAR_DECIMAL = "0.3823455333667027"

DEFAULT_TARGETS = (
    "0.3827090", "0.38268", "0.38260", "0.38250", "0.38240",
)
# Valid because M >= 1-c and lambda <= 0 imply G+lambda(M-(1-c)) <= G.
DEFAULT_LAMBDAS = ("0", "-0.2", "-0.4", "-0.6", "-0.8", "-0.9", "-1.0", "-1.2")
# Widest-coordinate ties are resolved in this order; support coordinates expose
# mean-infeasible boxes before the more expensive weight subdivision.
SPLIT_TIE_ORDER = (3, 4, 5, 6, 7, 8, 0, 1, 2)


class GradientUnavailable(RuntimeError):
    pass


@dataclass(frozen=True)
class BoundResult:
    gap_lower: arb
    objective_lower: arb
    direct_gap: arb
    denominator: arb
    numerator: arb
    method: str
    centered_usable: bool

    @property
    def gap_key(self) -> float:
        return float(self.gap_lower)

    @property
    def objective_key(self) -> float:
        return float(self.objective_lower)


@dataclass(frozen=True)
class RunStats:
    target: Fraction
    beta: Fraction
    processed: int
    cleared: int
    infeasible: int
    split: int
    residual: int
    elapsed: float
    stop_reason: str
    cover_objective_lower: float
    cover_gap_lower: float
    largest_objective_lower: float
    largest_gap_lower: float
    centered_uses: int
    candidate_objective_lower: float
    candidate_objective_upper: float
    candidate_margin: float


@dataclass(frozen=True)
class Jet:
    """First-order interval jet used only for rigorous mean-value forms."""

    value: arb
    gradient: tuple[arb, ...]

    @classmethod
    def constant(cls, value: Any, dimension: int) -> "Jet":
        return cls(arb(value), (ZERO,) * dimension)

    def _coerce(self, other: Any) -> "Jet":
        if isinstance(other, Jet):
            if len(other.gradient) != len(self.gradient):
                raise ValueError("jet dimensions differ")
            return other
        return Jet.constant(other, len(self.gradient))

    def __add__(self, other: Any) -> "Jet":
        other = self._coerce(other)
        return Jet(self.value + other.value,
                   tuple(x + y for x, y in zip(self.gradient, other.gradient)))

    __radd__ = __add__

    def __neg__(self) -> "Jet":
        return Jet(-self.value, tuple(-entry for entry in self.gradient))

    def __sub__(self, other: Any) -> "Jet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "Jet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "Jet":
        other = self._coerce(other)
        return Jet(
            self.value * other.value,
            tuple(x * other.value + self.value * y
                  for x, y in zip(self.gradient, other.gradient)),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "Jet":
        other = self._coerce(other)
        if not (other.value.lower() > 0 or other.value.upper() < 0):
            raise GradientUnavailable("jet division by interval containing zero")
        denominator = other.value * other.value
        return Jet(
            self.value / other.value,
            tuple((x * other.value - self.value * y) / denominator
                  for x, y in zip(self.gradient, other.gradient)),
        )


class _JetOps(_ScalarOps):
    def __init__(self, dimension: int):
        self.dimension = dimension
        super().__init__(self._entropy, Jet.constant(1, dimension))
        self._corner_ops = _ArbOps(True)

    def quotient(self, numerator: Jet, denominator: Jet) -> None:
        # The pilot centers the gap N-D, not the quotient.  Avoid needlessly
        # rejecting boxes whose entropy denominator contains zero.
        return None

    def xy_argument(self, x: Jet, y: Jet) -> Jet:
        result = x * y
        return Jet(self._corner_ops.xy_argument(x.value, y.value),
                   result.gradient)

    def pi_argument(self, x: Jet, y: Jet) -> Jet:
        one = self.one
        result = x * y + x * (one - x) * y * (one - y)
        return Jet(self._corner_ops.pi_argument(x.value, y.value),
                   result.gradient)

    @staticmethod
    def _entropy(argument: Jet) -> Jet:
        lo, hi = _semantic_bounds(argument.value, ZERO, ONE)
        if not (lo > ZERO and hi < ONE):
            raise GradientUnavailable("entropy derivative is unbounded at 0 or 1")
        # h'(z)=log((1-z)/z), decreasing on (0,1).
        derivative = arb_hull(
            ((ONE - hi) / hi).log(),
            ((ONE - lo) / lo).log(),
        )
        if not derivative.is_finite():
            raise GradientUnavailable("non-finite entropy derivative")
        return Jet(h_arb(argument.value),
                   tuple(derivative * entry for entry in argument.gradient))


class _NoQuotientArbOps(_ScalarOps):
    def __init__(self):
        super().__init__(h_arb, ONE)

    def quotient(self, numerator: arb, denominator: arb) -> None:
        return None


def _box_arbs(box: Box) -> tuple[arb, ...]:
    return tuple(arb(lo).union(arb(hi)) for lo, hi in box)


def _box_midpoints(box: Box) -> tuple[float, ...]:
    return tuple((lo + hi) / 2.0 for lo, hi in box)


def _simplex_vertices(box: Box) -> tuple[tuple[float, float], ...]:
    l1, u1 = box[0]
    l2, u2 = box[1]
    if l1 + l2 > 1.0:
        return ()
    candidates = {
        (x, y) for x in (l1, u1) for y in (l2, u2)
        if x + y <= 1.0
    }
    for x in (l1, u1):
        y = 1.0 - x
        if l2 <= y <= u2:
            candidates.add((x, y))
    for y in (l2, u2):
        x = 1.0 - y
        if l1 <= x <= u1:
            candidates.add((x, y))
    return tuple(sorted(candidates))


def _mean_at_corner(a1: float, a2: float, q: float,
                    points: Sequence[float]) -> arb:
    values = tuple(arb(value) for value in (a1, a2, q, *points))
    answer = mean_value(values)
    if not answer.is_finite():
        raise ValueError("non-finite corner mean")
    return answer


def mean_corner_range(box: Box) -> Optional[arb]:
    """Exact range of the multilinear mean on box intersect simplex.

    The mean is coordinatewise increasing in every support coordinate.  At
    fixed support endpoints it is affine in q and in (a1,a2), so its extrema
    occur at q endpoints and vertices of the clipped simplex rectangle.
    """
    vertices = _simplex_vertices(box)
    if not vertices:
        return None
    lows = tuple(interval[0] for interval in box[3:])
    highs = tuple(interval[1] for interval in box[3:])
    values = (
        _mean_at_corner(a1, a2, q, points)
        for points in (lows, highs)
        for q in box[2]
        for a1, a2 in vertices
    )
    result = next(values)
    for value in values:
        result = result.union(value)
    if not result.is_finite():
        raise ValueError("non-finite mean range")
    return result


def _raw_gap_at(values: Sequence[arb], beta: arb) -> tuple[arb, arb]:
    terms = _formula(values, beta, _NoQuotientArbOps())
    gap = terms.numerator - terms.ehx
    mean = mean_value(values)
    if not gap.is_finite() or not mean.is_finite():
        raise ValueError("non-finite centered point value")
    return gap, mean


def _gap_mean_gradient(box: Box, beta: arb) -> tuple[tuple[arb, ...],
                                                     tuple[arb, ...]]:
    dimension = len(box)
    jets = []
    for index, value in enumerate(_box_arbs(box)):
        gradient = [ZERO] * dimension
        gradient[index] = ONE
        jets.append(Jet(value, tuple(gradient)))
    ops = _JetOps(dimension)
    terms = _formula(tuple(jets), beta, ops)
    gap = terms.numerator - terms.ehx

    a1, a2, q, b0, b2, b4, b1, b3, b5 = jets
    one = ops.one
    a3 = one - a1 - a2
    mean = ((one - q) * (a1 * b0 + a2 * b2 + a3 * b4)
            + q * (a1 * b1 + a2 * b3 + a3 * b5))
    if any(not entry.is_finite() for entry in gap.gradient + mean.gradient):
        raise GradientUnavailable("non-finite interval gradient")
    return gap.gradient, mean.gradient


def centered_shift_enclosure(box: Box, target: Fraction, beta: arb,
                              lam: arb,
                              gradients: Optional[tuple[tuple[arb, ...],
                                                        tuple[arb, ...]]] = None) -> arb:
    """MVT enclosure of G+lambda*(M-(1-c)) over the whole box."""
    if gradients is None:
        gradients = _gap_mean_gradient(box, beta)
    gap_gradient, mean_gradient = gradients
    centers = _box_midpoints(box)
    center_arbs = tuple(arb(value) for value in centers)
    gap_center, mean_center = _raw_gap_at(center_arbs, beta)
    threshold = ONE - arb_fraction(target)
    answer = gap_center + lam * (mean_center - threshold)
    for index, ((lo, hi), center) in enumerate(zip(box, centers)):
        displacement = arb(lo - center).union(arb(hi - center))
        derivative = gap_gradient[index] + lam * mean_gradient[index]
        answer += derivative * displacement
    if not answer.is_finite():
        raise GradientUnavailable("non-finite centered enclosure")
    return answer


def _max_lower(candidates: Iterable[tuple[arb, str]]) -> tuple[arb, str]:
    candidates = list(candidates)
    if not candidates:
        raise ValueError("no lower-bound candidates")
    # Scheduling comparisons need not be certified: every candidate is itself
    # a certified lower bound, so selecting either overlapping value is sound.
    return max(candidates, key=lambda item: float(item[0]))


def box_bound(box: Box, target: Fraction, beta_fraction: Fraction,
              lambdas: Sequence[Fraction]) -> BoundResult:
    beta = arb_fraction(beta_fraction)
    terms = evaluate_arb(_box_arbs(box), beta, monotone_corners=True)
    direct_gap = terms.numerator - terms.ehx
    candidates: list[tuple[arb, str]] = [(direct_gap.lower(), "direct")]
    centered_usable = False

    # Full gradients require all entropy arguments to stay strictly interior.
    # Weight/mixture coordinates may still be wider: entropy singularities
    # occur only in the six support coordinates.
    if (max(hi - lo for lo, hi in box[3:]) <= 0.25
            and all(lo > 0.0 and hi < 1.0 for lo, hi in box[3:])):
        try:
            gradients = _gap_mean_gradient(box, beta)
            centered_usable = True
            for lam_fraction in lambdas:
                lam = arb_fraction(lam_fraction)
                enclosure = centered_shift_enclosure(
                    box, target, beta, lam, gradients,
                )
                candidates.append((enclosure.lower(),
                                   f"center(lambda={lam_fraction})"))
        except GradientUnavailable:
            centered_usable = False

    gap_lower, method = _max_lower(candidates)
    numerator_lo, _ = _semantic_bounds(terms.numerator, ZERO, None)
    denominator_lo, denominator_hi = _semantic_bounds(terms.ehx, ZERO, None)
    if denominator_hi > ZERO:
        direct_ratio = numerator_lo / denominator_hi
    else:
        # A zero-entropy point has no defined quotient in frankl5.m.  Zero is a
        # conservative lower bound for limiting positive-denominator points.
        direct_ratio = ZERO
    objective_candidates = [direct_ratio.lower()]
    if gap_lower >= ZERO:
        objective_candidates.append(ONE)
    elif denominator_lo > ZERO:
        objective_candidates.append((ONE + gap_lower / denominator_lo).lower())
    objective_lower = max(objective_candidates, key=float)
    if not all(value.is_finite() for value in (
            gap_lower, objective_lower, direct_gap, terms.ehx,
            terms.numerator)):
        raise ValueError("NaN/non-finite box bound")
    return BoundResult(
        gap_lower=gap_lower,
        objective_lower=objective_lower,
        direct_gap=direct_gap,
        denominator=terms.ehx,
        numerator=terms.numerator,
        method=method,
        centered_usable=centered_usable,
    )


def _widest_coordinate(box: Box) -> int:
    widths = tuple(hi - lo for lo, hi in box)
    widest = max(widths)
    return next(index for index in SPLIT_TIE_ORDER if widths[index] == widest)


def _bisect(box: Box) -> tuple[Box, Box]:
    coordinate = _widest_coordinate(box)
    lo, hi = box[coordinate]
    midpoint = (lo + hi) / 2.0
    left, right = list(box), list(box)
    left[coordinate] = (lo, midpoint)
    right[coordinate] = (midpoint, hi)
    return tuple(left), tuple(right)


def _reduced_candidate(target: Fraction, beta: Fraction,
                       dps: int = 80) -> tuple[arb, float]:
    """Numerically minimize the observed two-atom face, then enclose its value."""
    with mpmath.workdps(dps):
        threshold = mpmath.mpf(target.numerator) / target.denominator
        threshold = 1 - threshold
        beta_mp = mpmath.mpf(beta.numerator) / beta.denominator

        def h(z: mpmath.mpf) -> mpmath.mpf:
            if z == 0 or z == 1:
                return mpmath.mpf(0)
            return -(z * mpmath.log(z) + (1 - z) * mpmath.log(1 - z))

        def face(x: mpmath.mpf) -> mpmath.mpf:
            p = threshold / x
            xy = x * x
            pi = xy + (x * (1 - x)) ** 2
            return p * ((1 - beta_mp) * h(xy) + beta_mp * h(pi)) / h(x)

        x0 = stationary_parameters(dps)["x"]
        derivative = lambda z: mpmath.diff(face, z)
        try:
            x = mpmath.findroot(derivative, x0)
        except (ValueError, ZeroDivisionError):
            x = x0
        if not threshold < x < 1:
            x = x0
        p = threshold / x
        values_mp = (p, 1 - p, mpmath.mpf(0), x,
                     mpmath.mpf(0), mpmath.mpf(0),
                     mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(0))
        values_arb = tuple(arb(mpmath.nstr(value, dps)) for value in values_mp)
        terms = evaluate_arb(values_arb, arb_fraction(beta))
        if terms.objective is None or not terms.objective.is_finite():
            raise ValueError("candidate objective is non-finite")
        face_value = float(face(x))
        return terms.objective, face_value - 1.0


def run_target(target: Fraction, beta: Fraction,
               box_budget: int, seconds: float,
               lambdas: Sequence[Fraction]) -> RunStats:
    root: Box = tuple((0.0, 1.0) for _ in range(9))
    threshold = ONE - arb_fraction(target)
    started = time.monotonic()
    serial = 0
    heap: list[tuple[float, int, Box, BoundResult]] = []
    processed = cleared = infeasible = split = centered_uses = 0
    largest_objective = 0.0
    largest_gap = -math.inf

    def classify(box: Box) -> None:
        nonlocal serial, processed, cleared, infeasible, centered_uses
        nonlocal largest_objective, largest_gap
        mean_range = mean_corner_range(box)
        processed += 1
        if mean_range is None or mean_range.upper() < threshold:
            infeasible += 1
            return
        bound = box_bound(box, target, beta, lambdas)
        centered_uses += int(bound.centered_usable)
        largest_objective = max(largest_objective, bound.objective_key)
        largest_gap = max(largest_gap, bound.gap_key)
        if bound.gap_lower >= ZERO or bound.objective_lower >= ONE:
            cleared += 1
            return
        serial += 1
        heapq.heappush(heap, (bound.objective_key, serial, box, bound))

    classify(root)
    stop_reason = "complete"
    while heap:
        elapsed = time.monotonic() - started
        if processed + 2 > box_budget:
            stop_reason = "box_budget"
            break
        if elapsed >= seconds:
            stop_reason = "wall_clock"
            break
        _, _, box, _ = heapq.heappop(heap)
        left, right = _bisect(box)
        split += 1
        classify(left)
        classify(right)

    elapsed = time.monotonic() - started
    residual = len(heap)
    if heap:
        cover_objective = min(item[3].objective_key for item in heap)
        cover_gap = min(item[3].gap_key for item in heap)
    else:
        cover_objective = 1.0
        cover_gap = 0.0
    candidate_interval, candidate_margin = _reduced_candidate(target, beta)
    return RunStats(
        target=target,
        beta=beta,
        processed=processed,
        cleared=cleared,
        infeasible=infeasible,
        split=split,
        residual=residual,
        elapsed=elapsed,
        stop_reason=stop_reason,
        cover_objective_lower=cover_objective,
        cover_gap_lower=cover_gap,
        largest_objective_lower=largest_objective,
        largest_gap_lower=largest_gap,
        centered_uses=centered_uses,
        candidate_objective_lower=float(candidate_interval.lower()),
        candidate_objective_upper=float(candidate_interval.upper()),
        candidate_margin=candidate_margin,
    )


def _width(value: arb) -> float:
    return float(value.upper() - value.lower())


def print_pathology_diagnostics(beta: Fraction,
                                target: Fraction) -> tuple[float, float]:
    root: Box = tuple((0.0, 1.0) for _ in range(9))
    root_terms = evaluate_arb(_box_arbs(root), arb_fraction(beta), True)
    naive_quotient = root_terms.numerator / root_terms.ehx
    print("PROVED [Arb finite check]: the root-box denominator contains zero="
          f"{root_terms.ehx.contains(0)}; certified ehx enclosure={root_terms.ehx}.")
    print("PROVED [Arb finite check]: naive root quotient is finite="
          f"{naive_quotient.is_finite()}; the evaluator withholds division and certifies N-D instead.")
    print("PROVED [NaN discipline by construction]: no square root occurs; entropy endpoints are handled before logarithms, and every accepted Arb objective/gradient bound passes is_finite().")

    # A fully interior, partially mean-feasible box makes all compared rules
    # applicable on exactly the same geometry.
    benchmark: Box = (
        (0.86, 0.92), (0.07, 0.13), (0.01, 0.05),
        (0.66, 0.72), (0.02, 0.08), (0.10, 0.18),
        (0.54, 0.62), (0.25, 0.35), (0.76, 0.84),
    )
    values = _box_arbs(benchmark)
    natural = evaluate_arb(values, arb_fraction(beta), False)
    corner = evaluate_arb(values, arb_fraction(beta), True)
    natural_gap = natural.numerator - natural.ehx
    corner_gap = corner.numerator - corner.ehx
    beta_free = evaluate_arb(values, arb(0).union(arb(1)), True)
    beta_free_gap = beta_free.numerator - beta_free.ehx
    mean_natural = mean_value(values)
    mean_corner = mean_corner_range(benchmark)
    if mean_corner is None:
        raise AssertionError("benchmark unexpectedly misses simplex")

    gradients = _gap_mean_gradient(benchmark, arb_fraction(beta))
    centered = centered_shift_enclosure(
        benchmark, target, arb_fraction(beta), ZERO, gradients,
    )
    kkt_lambda = arb_fraction(Fraction("-0.9"))
    shifted = centered_shift_enclosure(
        benchmark, target, arb_fraction(beta), kkt_lambda, gradients,
    )
    print("NUMERICAL [same-box Arb widths]: natural gap=%.9g; monotone-kernel-corner gap=%.9g; centered gap=%.9g; KKT-shifted(lambda=-0.9) centered functional=%.9g."
          % (_width(natural_gap), _width(corner_gap),
             _width(centered), _width(shifted)))
    print("NUMERICAL [same-box Arb widths]: beta left as an interval variable in [0,1] gives gap=%.9g; eliminating it from the 9-D box by fixing frankl5.m's reported rational gives gap=%.9g."
          % (_width(beta_free_gap), _width(corner_gap)))
    print("NUMERICAL [same-box Arb widths]: natural mean=%.9g; exact monotone-corner mean over the simplex-clipped box=%.9g."
          % (_width(mean_natural), _width(mean_corner)))
    print("PROVED [calculus + Arb]: xy and xy+xy(1-x)(1-y) are coordinatewise nondecreasing on [0,1]^2, so their corner ranges are valid two-sided enclosures.")
    print("PROVED [feasibility shift]: each lambda<=0 bound encloses G+lambda*(M-(1-c))<=G on M>=1-c; no NaN or non-finite centered enclosure was accepted.")
    return _width(shifted), max(hi - lo for lo, hi in benchmark)



def print_table(stats: Sequence[RunStats], box_budget: int,
                seconds: float) -> None:
    print("NUMERICAL [measured one-core pilot]: fixed budget per target: "
          f"{box_budget} evaluated boxes or {seconds:g} seconds, whichever first.")
    print("NUMERICAL [trend table]: c | processed | cleared | infeasible | split | residual | centered | cover_obj_lb | cover_gap_lb | largest_obj_lb | candidate_margin | candidate_obj_interval | seconds | stop")
    for item in stats:
        target_text = format(float(item.target), ".15g")
        candidate = "[%.12g,%.12g]" % (
            item.candidate_objective_lower, item.candidate_objective_upper)
        print("NUMERICAL [measured]: %s | %d | %d | %d | %d | %d | %d | %.9g | %+.9g | %.9g | %+.9g | %s | %.3f | %s"
              % (target_text, item.processed, item.cleared, item.infeasible,
                 item.split, item.residual, item.centered_uses,
                 item.cover_objective_lower, item.cover_gap_lower,
                 item.largest_objective_lower, item.candidate_margin,
                 candidate, item.elapsed, item.stop_reason))
    margins = [item.candidate_margin for item in stats]
    positive = [value for value in margins if value > 0.0]
    if positive:
        scale = max(positive) / min(positive)
        print("NUMERICAL [trend interpretation]: the observed two-atom-face margin grows by a factor %.6g across the target retreat; the direct 9-D cover cleared %d boxes in total at this resolution."
              % (scale, sum(item.cleared for item in stats)))



def print_verdict(stats: Sequence[RunStats], diagnostic_shift_width: float,
                  diagnostic_box_width: float) -> None:
    # The largest retreat (smallest c) is the easiest and therefore the honest
    # best case for extrapolating this direct rule set.
    easiest = min(stats, key=lambda item: float(item.target))
    resolved = easiest.cleared + easiest.infeasible
    leaves = resolved + easiest.residual
    resolved_fraction = resolved / leaves if leaves else 1.0
    if resolved_fraction <= 0.0:
        linear_floor = math.inf
    elif resolved_fraction >= 1.0:
        linear_floor = float(easiest.processed)
    else:
        linear_floor = easiest.processed / resolved_fraction
    linear_text = ("infinite from zero observed progress"
                   if not math.isfinite(linear_floor)
                   else f">={linear_floor:.3g}")
    # Optimistically assume the KKT-shifted MVT width becomes quadratic under
    # isotropic refinement.  This is much kinder than the observed boundary
    # behavior and therefore should not be read as an upper bound.
    margin = max(easiest.candidate_margin, sys.float_info.min)
    required_width = diagnostic_box_width * math.sqrt(
        margin / diagnostic_shift_width)
    isotropic_cells = (1.0 / required_width) ** 9
    print("CONJECTURED [box-count extrapolation from easiest measured target]: at c=%s the observed-volume floor is %s boxes; optimistic quadratic KKT scaling from the measured same-box width gives resolution %.3g and about %.3g isotropic 9-D cells."
          % (format(float(easiest.target), ".15g"), linear_text,
             required_width, isotropic_cells))
    print("CONJECTURED [feasibility verdict]: attack the structure first, not the raw nine-dimensional box.  A viable campaign needs certified symmetry/face reduction (q=0 and a two-atom P0, or an equivalent analytic elimination), plus sink-safe mixed centered forms; monotone corners and fixed beta alone do not remove the dimensional explosion.")


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--beta", default=LIU_BETA_DECIMAL,
                        help="exact decimal rational; default is frankl5.m Beta")
    parser.add_argument("--box-budget", type=int, default=50001)
    parser.add_argument("--seconds", type=float, default=60.0,
                        help="wall-clock cap per target")
    parser.add_argument("--prec", type=int, default=128,
                        help="Arb midpoint precision in bits")
    parser.add_argument("--lambdas", nargs="+", default=DEFAULT_LAMBDAS)
    parser.add_argument("--diagnostics-only", action="store_true")
    args = parser.parse_args(argv)
    if args.box_budget < 1 or args.seconds <= 0 or args.prec < 80:
        parser.error("positive budgets and --prec >= 80 are required")
    ctx.prec = args.prec
    targets = tuple(as_fraction(value) for value in args.targets)
    beta = as_fraction(args.beta)
    lambdas = tuple(as_fraction(value) for value in args.lambdas)
    if any(value > 0 for value in lambdas):
        parser.error("all feasibility-shift lambdas must be <= 0")

    source_constant = as_fraction(LIU_C_DECIMAL)
    c_star = as_fraction(C_STAR_DECIMAL)
    for target in targets:
        relation = "<" if target < source_constant else ">="
        star_relation = ">" if target > c_star else "<="
        print("PROVED [integer comparison of terminating decimals]: target %s %s Liu's reported decimal %s and %s c*=%s."
              % (format(float(target), ".15g"), relation, LIU_C_DECIMAL,
                 star_relation, C_STAR_DECIMAL))
    diagnostic_shift_width, diagnostic_box_width = \
        print_pathology_diagnostics(beta, targets[-1])
    if args.diagnostics_only:
        return
    results = [
        run_target(target, beta, args.box_budget, args.seconds, lambdas)
        for target in targets
    ]
    print_table(results, args.box_budget, args.seconds)
    print_verdict(results, diagnostic_shift_width, diagnostic_box_width)


if __name__ == "__main__":
    main()
