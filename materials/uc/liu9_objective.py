#!/usr/bin/env python3
"""Liu's nine-parameter objective (arXiv:2306.08824v1, Section V-B).

The source of truth transcribed here is
https://jingbol.web.illinois.edu/frankl5.m.  The single ``_formula`` routine is
used by the binary64, mpmath, and Arb front ends.  Entropies are in nats, as in
the Matlab source; the quotient is independent of the logarithm base.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Callable, Iterable, Optional, Sequence

import mpmath
from flint import arb, ctx

ctx.prec = max(ctx.prec, 160)

VARIABLE_NAMES = ("a1", "a2", "q", "b0", "b2", "b4", "b1", "b3", "b5")
LIU_C_DECIMAL = "0.382709087918741"
LIU_BETA_DECIMAL = "0.100052559862974"
LIU_X_DECIMAL = "0.690787593924988"
# The paper prints ...5457; solving its defining equations gives ...5466.
LIU_P_PAPER_DECIMAL = "0.893604513905457"
LIU_P_EQUATION_DECIMAL = "0.893604513905466"

SOURCE_OBJECTIVE_LINE = (
    "obj = ((1-beta)*ehxy+beta*ehpi)/ehx;"
)
SOURCE_ENTROPY_LINE = (
    "hxy=-xy.*log(xy)-(1-xy).*log(1-xy);"
)
SOURCE_BETA_LINE = "beta=Beta(j);"
SOURCE_MEAN_LINE = (
    "meancons =qb*(x(1)*x(4)+x(2)*x(5)+(1-x(1)-x(2))*x(6)) "
    "+ x(3)*(x(1)*x(7)+x(2)*x(8)+(1-x(1)-x(2))*x(9)) >=1-c;"
)


@dataclass(frozen=True)
class ObjectiveTerms:
    """The three expectations and their assembled objective."""

    ehxy: Any
    ehpi: Any
    ehx: Any
    numerator: Any
    objective: Optional[Any]


def as_fraction(value: str | int | Fraction) -> Fraction:
    """Parse a CLI decimal as the exact rational it denotes."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    return Fraction(value)


def arb_fraction(value: Fraction | str | int) -> arb:
    """An Arb enclosure of an exact rational."""
    value = as_fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def arb_hull(lo: arb, hi: arb) -> arb:
    """Outward-rounded Arb hull."""
    return arb(lo).union(arb(hi))


def _semantic_bounds(value: arb, lo: Optional[arb] = None,
                     hi: Optional[arb] = None) -> tuple[arb, arb]:
    """Endpoint balls clipped using a known semantic range."""
    if not value.is_finite():
        raise ValueError("non-finite Arb input")
    lower, upper = value.lower(), value.upper()
    if lo is not None and lower < lo:
        lower = arb(lo)
    if hi is not None and upper > hi:
        upper = arb(hi)
    if lower > upper:
        raise ValueError("Arb enclosure misses its semantic range")
    return lower, upper


def _h_float(x: float) -> float:
    if x == 0.0 or x == 1.0:
        return 0.0
    if not 0.0 < x < 1.0:
        raise ValueError(f"entropy argument outside [0,1]: {x!r}")
    return -(x * math.log(x) + (1.0 - x) * math.log1p(-x))


def _h_mp(x: mpmath.mpf) -> mpmath.mpf:
    if x == 0 or x == 1:
        return mpmath.mpf(0)
    if not 0 < x < 1:
        raise ValueError(f"entropy argument outside [0,1]: {x!r}")
    return -(x * mpmath.log(x) + (1 - x) * mpmath.log(1 - x))


def _h_arb_point(x: arb) -> arb:
    """Binary entropy at an endpoint ball known to lie in [0,1]."""
    if x <= 0 or x >= 1:
        return arb(0)
    ans = -(x * x.log() + (arb(1) - x) * (arb(1) - x).log())
    if not ans.is_finite():
        raise ValueError("non-finite entropy endpoint")
    return ans


def h_arb(value: arb) -> arb:
    """Exact range enclosure of binary entropy over an Arb ball.

    The input has the semantics of a probability.  Natural interval products
    can protrude slightly outside [0,1], so endpoint balls are clipped to that
    known domain before any logarithm is taken.  No logarithm of a ball
    containing zero is attempted.
    """
    zero, one, half = arb(0), arb(1), arb(1) / 2
    lo, hi = _semantic_bounds(arb(value), zero, one)
    hlo = _h_arb_point(lo) if lo > zero else zero
    hhi = _h_arb_point(hi) if hi < one else zero
    out = arb_hull(hlo, hhi)
    if not (hi <= half or lo >= half):
        out = out.union(arb(2).log())
    if not out.is_finite():
        raise ValueError("non-finite entropy enclosure")
    return out


def _nonnegative_product_arb(*values: arb) -> arb:
    """Product enclosure when every factor is semantically nonnegative."""
    lo, hi = arb(1), arb(1)
    for value in values:
        vlo, vhi = _semantic_bounds(arb(value), arb(0), None)
        lo *= vlo
        hi *= vhi
    return arb_hull(lo, hi)


def _nonnegative_sum_arb(values: Iterable[arb]) -> arb:
    lo, hi = arb(0), arb(0)
    for value in values:
        vlo, vhi = _semantic_bounds(value, arb(0), None)
        lo += vlo
        hi += vhi
    return arb_hull(lo, hi)


class _ScalarOps:
    def __init__(self, entropy: Callable[[Any], Any], one: Any):
        self.entropy = entropy
        self.one = one

    def nonnegative_product(self, *values: Any) -> Any:
        out = self.one
        for value in values:
            out *= value
        return out

    def nonnegative_sum(self, values: Iterable[Any]) -> Any:
        return sum(values, self.one - self.one)
    def quotient(self, numerator: Any, denominator: Any) -> Optional[Any]:
        return numerator / denominator


    def xy_argument(self, x: Any, y: Any) -> Any:
        return x * y

    def pi_argument(self, x: Any, y: Any) -> Any:
        return x * y + x * (self.one - x) * y * (self.one - y)


class _ArbOps(_ScalarOps):
    def __init__(self, monotone_corners: bool):
        super().__init__(h_arb, arb(1))
        self.monotone_corners = monotone_corners

    def nonnegative_product(self, *values: arb) -> arb:
        return _nonnegative_product_arb(*values)

    def nonnegative_sum(self, values: Iterable[arb]) -> arb:
        return _nonnegative_sum_arb(values)
    def quotient(self, numerator: arb, denominator: arb) -> Optional[arb]:
        dlo, _ = _semantic_bounds(denominator, arb(0), None)
        if not (dlo > 0):
            return None
        answer = numerator / denominator
        return answer if answer.is_finite() else None


    @staticmethod
    def _corner_range(function: Callable[[arb, arb], arb], x: arb,
                      y: arb) -> arb:
        xlo, xhi = _semantic_bounds(x, arb(0), arb(1))
        ylo, yhi = _semantic_bounds(y, arb(0), arb(1))
        return arb_hull(function(xlo, ylo), function(xhi, yhi))

    def xy_argument(self, x: arb, y: arb) -> arb:
        if not self.monotone_corners:
            return x * y
        return self._corner_range(lambda u, v: u * v, x, y)

    def pi_argument(self, x: arb, y: arb) -> arb:
        def kernel(u: arb, v: arb) -> arb:
            return u * v + u * (arb(1) - u) * v * (arb(1) - v)

        if not self.monotone_corners:
            return kernel(x, y)
        # d/du [uv + u(1-u)v(1-v)] >= 0 on [0,1]^2, and likewise
        # for v, so the two same-corner evaluations give the exact range.
        return self._corner_range(kernel, x, y)


def _formula(values: Sequence[Any], beta: Any, ops: _ScalarOps,
             a3_override: Optional[Any] = None) -> ObjectiveTerms:
    """Single exact formula used by all three arithmetic front ends."""
    if len(values) != 9:
        raise ValueError("Liu's objective needs exactly nine variables")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    one = ops.one
    a3 = one - a1 - a2 if a3_override is None else a3_override
    masses = (a1, a2, a3)
    p0 = (b0, b2, b4)
    p1 = (b1, b3, b5)
    qbar = one - q
    support = p0 + p1
    weights = tuple(
        ops.nonnegative_product(component, mass)
        for component in (qbar, q) for mass in masses
    )

    ehxy = ops.nonnegative_sum(
        ops.nonnegative_product(
            weights[i], weights[j],
            ops.entropy(ops.xy_argument(support[i], support[j])),
        )
        for i in range(6) for j in range(6)
    )

    def component_pi(points: Sequence[Any]) -> Any:
        return ops.nonnegative_sum(
            ops.nonnegative_product(
                masses[i], masses[j],
                ops.entropy(ops.pi_argument(points[i], points[j])),
            )
            for i in range(3) for j in range(3)
        )

    ehpi = ops.nonnegative_sum((
        ops.nonnegative_product(qbar, component_pi(p0)),
        ops.nonnegative_product(q, component_pi(p1)),
    ))
    ehx = ops.nonnegative_sum(
        ops.nonnegative_product(weights[i], ops.entropy(support[i]))
        for i in range(6)
    )
    numerator = ops.nonnegative_sum((
        ops.nonnegative_product(one - beta, ehxy),
        ops.nonnegative_product(beta, ehpi),
    ))
    objective = ops.quotient(numerator, ehx)
    return ObjectiveTerms(ehxy, ehpi, ehx, numerator, objective)


def _check_point(values: Sequence[Any], beta: Any) -> None:
    if len(values) != 9:
        raise ValueError("expected nine variables")
    a1, a2, q, *points = values
    if not (0 <= a1 <= 1 and 0 <= a2 <= 1 and a1 + a2 <= 1):
        raise ValueError("(a1,a2,1-a1-a2) is not on the simplex")
    if not 0 <= q <= 1 or any(not 0 <= point <= 1 for point in points):
        raise ValueError("q and b0,...,b5 must lie in [0,1]")
    if not 0 <= beta <= 1:
        raise ValueError("beta must lie in [0,1]")


def evaluate_float64(values: Sequence[float], beta: float) -> ObjectiveTerms:
    """Evaluate the exact transcription using IEEE binary64 arithmetic."""
    values = tuple(float(value) for value in values)
    beta = float(beta)
    _check_point(values, beta)
    return _formula(values, beta, _ScalarOps(_h_float, 1.0))


def evaluate_mpmath(values: Sequence[Any], beta: Any,
                    dps: int = 80) -> ObjectiveTerms:
    """Evaluate the exact transcription with arbitrary precision."""
    with mpmath.workdps(dps):
        converted = tuple(mpmath.mpf(str(value)) for value in values)
        beta_mp = mpmath.mpf(str(beta))
        _check_point(converted, beta_mp)
        terms = _formula(converted, beta_mp, _ScalarOps(_h_mp, mpmath.mpf(1)))
        # Unary plus rounds persistent results to the requested precision.
        return ObjectiveTerms(*(+value for value in (
            terms.ehxy, terms.ehpi, terms.ehx, terms.numerator,
            terms.objective,
        )))


def _a3_feasible_enclosure(a1: arb, a2: arb) -> arb:
    raw = arb(1) - a1 - a2
    lo, hi = _semantic_bounds(raw, arb(0), arb(1))
    return arb_hull(lo, hi)


def evaluate_arb(values: Sequence[arb], beta: arb,
                 monotone_corners: bool = True) -> ObjectiveTerms:
    """Certified two-sided Arb interval extension of the shared formula.

    ``values`` may be points or balls.  The call is for the feasible portion of
    the box; consequently the derived a3 range is intersected with [0,1].
    """
    if len(values) != 9:
        raise ValueError("expected nine Arb values")
    converted = tuple(arb(value) for value in values)
    if any(not value.is_finite() for value in converted) or not arb(beta).is_finite():
        raise ValueError("all Arb inputs must be finite")
    a3 = _a3_feasible_enclosure(converted[0], converted[1])
    terms = _formula(converted, arb(beta), _ArbOps(monotone_corners), a3)
    if any(value is not None and not value.is_finite() for value in (
            terms.ehxy, terms.ehpi, terms.ehx, terms.numerator,
            terms.objective)):
        raise ValueError("formula produced a non-finite Arb enclosure")
    return terms


def mean_value(values: Sequence[Any]) -> Any:
    """The left side of frankl5.m's sole nonlinear constraint."""
    if len(values) != 9:
        raise ValueError("expected nine variables")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    one = type(q)(1) if not isinstance(q, arb) else arb(1)
    a3 = one - a1 - a2
    return ((one - q) * (a1 * b0 + a2 * b2 + a3 * b4)
            + q * (a1 * b1 + a2 * b3 + a3 * b5))


def stationary_parameters(dps: int = 100) -> dict[str, mpmath.mpf]:
    """Solve Liu's displayed structural equations and constrained stationarity."""
    with mpmath.workdps(dps):
        one = mpmath.mpf(1)
        h = _h_mp
        hp = lambda z: mpmath.log((one - z) / z)
        x = mpmath.findroot(
            lambda z: z * z + z * z * (one + (one - z) ** 2) - one,
            mpmath.mpf(LIU_X_DECIMAL),
        )
        p = h(x) / h(x * x)
        c = one - p * x
        u = x * x
        t = u + (x * (one - x)) ** 2
        du = 2 * x
        dt = 2 * x + 2 * x * (one - x) ** 2 - 2 * x * x * (one - x)
        a_prime = hp(u) * du
        b_prime = hp(t) * dt
        beta = ((hp(x) + h(x) / x) / p - a_prime) / (b_prime - a_prime)
        return {name: +value for name, value in (
            ("x", x), ("p", p), ("c", c), ("beta", beta),
        )}


def structural_point(parameters: dict[str, Any]) -> tuple[Any, ...]:
    """A nine-vector representing P0=p*delta_x+(1-p)*delta_0, q=0."""
    p, x = parameters["p"], parameters["x"]
    zero = p - p
    one = zero + 1
    return (p, one - p, zero, x, zero, zero, zero, zero, zero)


def _arb_point_from_mp(value: mpmath.mpf, digits: int = 100) -> arb:
    return arb(mpmath.nstr(value, digits))


def _random_feasible_points(count: int, c: float,
                            seed: int) -> list[tuple[float, ...]]:
    rng = random.Random(seed)
    threshold = 1.0 - c
    answer: list[tuple[float, ...]] = []
    while len(answer) < count:
        raw = [rng.expovariate(1.0) for _ in range(3)]
        total = sum(raw)
        masses = [value / total for value in raw]
        q = rng.uniform(0.01, 0.99)
        points = [rng.uniform(0.01, 0.99) for _ in range(6)]
        values = (masses[0], masses[1], q, *points)
        if mean_value(values) >= threshold:
            answer.append(values)
    return answer


def print_feasible_set() -> None:
    print("PROVED [verbatim frankl5.m constraints]: feasible set for fixed rational c and fixed beta:")
    print("PROVED [bounds + ca12]: 0 <= a1,a2,q,b0,b2,b4,b1,b3,b5 <= 1; a1+a2 <= 1; a3=1-a1-a2.")
    print("PROVED [meancons]: (1-q)[a1*b0+a2*b2+a3*b4] + q[a1*b1+a2*b3+a3*b5] >= 1-c.")
    print(f"PROVED [source line]: `{SOURCE_MEAN_LINE}`")
    print("PROVED [absence of constraints in frankl5.m]: there are no ordering constraints; permutations of the three paired atoms are symmetries only.")
    print(f"PROVED [source line]: fixed beta is external to the nine variables: `{SOURCE_BETA_LINE}`")
    print(f"PROVED [source line]: objective transcription: `{SOURCE_OBJECTIVE_LINE}`")
    print(f"PROVED [source line]: boundary ambiguity: `{SOURCE_ENTROPY_LINE}` uses 0*log(0), while no explicit ehx>0 constraint or endpoint convention is stated.")


def validate(random_points: int = 20, seed: int = 230608824,
             dps: int = 100) -> None:
    print_feasible_set()
    params = stationary_parameters(dps)
    with mpmath.workdps(dps):
        point_mp = structural_point(params)
        beta_mp = mpmath.mpf(LIU_BETA_DECIMAL)
        p_reproduction_delta = abs(
            params["p"] - mpmath.mpf(LIU_P_EQUATION_DECIMAL))
        beta_source_delta = abs(
            params["beta"] - mpmath.mpf(LIU_BETA_DECIMAL))
    point_float = tuple(float(value) for value in point_mp)
    f64 = evaluate_float64(point_float, float(beta_mp))
    high = evaluate_mpmath(
        tuple(mpmath.nstr(value, dps) for value in point_mp),
        mpmath.nstr(beta_mp, dps), dps,
    )
    point_arb = tuple(_arb_point_from_mp(value, dps) for value in point_mp)
    interval = evaluate_arb(point_arb, _arb_point_from_mp(beta_mp, dps))

    print("NUMERICAL [100-digit solve of (90),(91) and constrained stationarity]: "
          f"x*={mpmath.nstr(params['x'], 30)}")
    print("NUMERICAL [same solve]: p*=" + mpmath.nstr(params["p"], 30))
    print("NUMERICAL [frankl5.m reported value used by all evaluators]: beta*=" + LIU_BETA_DECIMAL)
    print("NUMERICAL [same solve]: c'=1-p*x=" + mpmath.nstr(params["c"], 30))
    print("NUMERICAL [paper-vs-equation transcription check]: paper prints p*="
          + LIU_P_PAPER_DECIMAL + "; its defining equations give "
          + mpmath.nstr(params["p"], 18)
          + ", differing from the independently reproduced decimal "
          + LIU_P_EQUATION_DECIMAL + " by "
          + mpmath.nstr(p_reproduction_delta, 6)
          + ".")
    print("NUMERICAL [stationarity cross-check on the transcribed formula]: beta="
          + mpmath.nstr(params["beta"], 30)
          + "; it differs from frankl5.m by "
          + mpmath.nstr(beta_source_delta, 6)
          + ", consistent with the source's numerical optimizer rather than an exact identity.")
    print("NUMERICAL [float64 shared formula at structural point]: objective="
          f"{f64.objective:.17g}, ehxy={f64.ehxy:.17g}, ehpi={f64.ehpi:.17g}, ehx={f64.ehx:.17g}")
    print("NUMERICAL [mpmath shared formula at structural point]: objective="
          + mpmath.nstr(high.objective, 50))
    print("PROVED [Arb evaluation of the exact decimal enclosure]: objective in "
          + str(interval.objective))
    if not interval.objective.contains(_arb_point_from_mp(high.objective, dps)):
        raise AssertionError("Arb optimum enclosure does not contain mpmath value")
    if abs(f64.objective - float(high.objective)) > 5e-14:
        raise AssertionError("float64 and mpmath disagree at structural point")

    points = _random_feasible_points(random_points, float(params["c"]), seed)
    worst_float = 0.0
    widest_arb = 0.0
    for index, values in enumerate(points):
        f = evaluate_float64(values, float(beta_mp))
        m = evaluate_mpmath(tuple(repr(value) for value in values),
                            mpmath.nstr(beta_mp, dps), dps)
        av = tuple(arb(repr(value)) for value in values)
        a = evaluate_arb(av, _arb_point_from_mp(beta_mp, dps))
        probe = _arb_point_from_mp(m.objective, dps)
        if not a.objective.contains(probe):
            raise AssertionError(f"Arb enclosure misses random point {index}")
        worst_float = max(worst_float, abs(f.objective - float(m.objective)))
        widest_arb = max(widest_arb,
                         float(a.objective.upper() - a.objective.lower()))
        if not all(value.is_finite() for value in (
                a.ehxy, a.ehpi, a.ehx, a.numerator, a.objective)):
            raise AssertionError(f"non-finite Arb result at random point {index}")
    print(f"NUMERICAL [deterministic random feasible-point cross-check]: {random_points} points, max |float64-mpmath|={worst_float:.3e}.")
    print(f"PROVED [Arb containment at rationalized random inputs]: {random_points}/{random_points} mpmath approximations enclosed; all intervals finite; max point-enclosure width={widest_arb:.3e}.")


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-points", type=int, default=20)
    parser.add_argument("--seed", type=int, default=230608824)
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args(argv)
    validate(args.random_points, args.seed, args.dps)


if __name__ == "__main__":
    main()
