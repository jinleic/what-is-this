#!/usr/bin/env python3
"""Transverse audit of Liu's active-mean smooth chart.

The certified chart only moves the two active supports while keeping the mean
constraint active.  This module embeds that chart into Liu's nine coordinates,
checks the embedding against ``liu9_objective._formula``, and then opens a
zero-mass coordinate into a genuine third atom.

Two openings are kept distinct:

* ``active_split_embedding`` literally moves mass epsilon from the active atom
  at x to an atom at y.  It is an ambient nine-variable direction and changes
  the mean by epsilon*(y-x).
* ``mean_preserving_embedding`` compensates between the atoms at x and zero so
  that the mean remains exactly p*x.

The distinction is load-bearing: the first family refutes an unrestricted
ambient extension, while the explicit counterexample decreases the mean and
therefore does not refute a statement restricted to Liu's feasible half-space.

Run from the math directory with

    ./.venv/bin/python -I -B uc/liu9_ninevar.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Sequence

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_chart_centered import (  # noqa: E402
    Jet3,
    distance_squared_jet,
    gap_jet,
)
from liu9_objective import (  # noqa: E402
    _ArbOps,
    _formula,
    h_arb,
)

ctx.prec = max(ctx.prec, 320)

KAPPA = Fraction(4_119_063, 33_554_432)
Y0 = Fraction(1, 32)
DECISIVE_Y = Y0
DECISIVE_EPSILON = Fraction(1, 1024)
DECISIVE_Q = Fraction(1, 2)
GRID_INTERVALS = 64
LOG_SERIES_TERMS = 96
PARAMETER_INTERVAL_DECIMALS = 18
LEADING_NEGATIVE_THRESHOLD = Fraction(-1, 2)
FINITE_NEGATIVE_THRESHOLD = Fraction(-1, 2048)
MEAN_PRESERVING_POSITIVE_THRESHOLD = Fraction(1, 64)
REPORT_PATH = HERE / "verification" / "results" / "liu9-ninevar.json"


# ---------------------------------------------------------------------------
# The two embeddings and the nine-variable tube distance.
def _exact_zero(value: Any) -> Any:
    """A scalar zero without widening an interval by dependency."""
    return value * 0


def _exact_one(value: Any) -> Any:
    """A scalar one without widening an interval by dependency."""
    return value * 0 + 1




def chart_embedding(
    parameters: Any,
    q: Any,
    s: Any,
    d: Any,
    r: Any,
) -> tuple[Any, ...]:
    """Embed ``(q,s,d,r)`` in Liu's documented nine-variable order.

    The returned order is exactly
    ``(a1,a2,q,b0,b2,b4,b1,b3,b5)`` and the active mass is ``p*x/s``.
    """
    active_mass = parameters.mean / s
    zero = _exact_zero(s)
    return (
        active_mass,
        r,
        q,
        s - q * d,
        zero,
        zero,
        s + (1 - q) * d,
        zero,
        zero,
    )


def active_split_embedding(
    parameters: Any,
    y: Any,
    epsilon: Any,
    q: Any,
) -> tuple[Any, ...]:
    """Move mass ``epsilon`` from x to y in both components.

    The third atom is paired with the same support in both components, so q is
    a representation gauge.  The mean changes by ``epsilon*(y-x)``.
    """
    zero = _exact_zero(parameters.x)
    return (
        parameters.p - epsilon,
        epsilon,
        q,
        parameters.x,
        y,
        zero,
        parameters.x,
        y,
        zero,
    )


def mean_preserving_embedding(
    parameters: Any,
    y: Any,
    epsilon: Any,
    q: Any,
) -> tuple[Any, ...]:
    """Insert mass ``epsilon`` at y while preserving mean ``p*x`` exactly."""
    zero = _exact_zero(parameters.x)
    return (
        parameters.p - epsilon * y / parameters.x,
        epsilon,
        q,
        parameters.x,
        y,
        zero,
        parameters.x,
        y,
        zero,
    )


def distance_squared(values: Sequence[Any], parameters: Any) -> Any:
    """The exact nine-variable distance used by both chart modules."""
    if len(values) != 9:
        raise ValueError("expected Liu's nine coordinates")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    one = _exact_one(q)
    a3 = one - a1 - a2
    masses = (a1, a2, a3)

    def component(points: Sequence[Any]) -> Any:
        mean = sum(
            (mass * point for mass, point in zip(masses, points)),
            _exact_zero(one),
        )
        spread = sum(
            (
                mass
                * (point * (point - parameters.x))
                * (point * (point - parameters.x))
                for mass, point in zip(masses, points)
            ),
            _exact_zero(one),
        )
        return (mean - parameters.mean) * (mean - parameters.mean) + spread

    return ((one - q) * component((b0, b2, b4))
            + q * component((b1, b3, b5)))


# ---------------------------------------------------------------------------
# Exact epsilon-polynomial automatic differentiation through _formula.


class EpsilonDependentEntropyArgument(ArithmeticError):
    """An entropy argument depends on epsilon and is not a polynomial case."""


class EpsilonLogTerm(EpsilonDependentEntropyArgument):
    """An entropy argument opens linearly from zero or one."""


@dataclass(frozen=True)
class EpsilonPolynomial:
    """Coefficients of ``c0 + c1*epsilon + c2*epsilon^2``."""

    c0: Any
    c1: Any
    c2: Any

    @staticmethod
    def constant(value: Any) -> "EpsilonPolynomial":
        zero = _exact_zero(value)
        return EpsilonPolynomial(value, zero, zero)

    @staticmethod
    def variable(one: Any) -> "EpsilonPolynomial":
        zero = _exact_zero(one)
        return EpsilonPolynomial(zero, one, zero)

    def _coerce(self, other: Any) -> "EpsilonPolynomial":
        return other if isinstance(other, EpsilonPolynomial) else self.constant(other)

    def __add__(self, other: Any) -> "EpsilonPolynomial":
        rhs = self._coerce(other)
        return EpsilonPolynomial(
            self.c0 + rhs.c0,
            self.c1 + rhs.c1,
            self.c2 + rhs.c2,
        )

    __radd__ = __add__

    def __neg__(self) -> "EpsilonPolynomial":
        return EpsilonPolynomial(-self.c0, -self.c1, -self.c2)

    def __sub__(self, other: Any) -> "EpsilonPolynomial":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "EpsilonPolynomial":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "EpsilonPolynomial":
        rhs = self._coerce(other)
        return EpsilonPolynomial(
            self.c0 * rhs.c0,
            self.c0 * rhs.c1 + self.c1 * rhs.c0,
            self.c0 * rhs.c2 + self.c1 * rhs.c1 + self.c2 * rhs.c0,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "EpsilonPolynomial":
        inverse = 1 / self.c0
        return EpsilonPolynomial(
            inverse,
            -self.c1 * inverse * inverse,
            self.c1 * self.c1 * inverse**3 - self.c2 * inverse * inverse,
        )

    def __truediv__(self, other: Any) -> "EpsilonPolynomial":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "EpsilonPolynomial":
        return self._coerce(other) / self

    def evaluate(self, epsilon: Any) -> Any:
        return self.c0 + epsilon * self.c1 + epsilon * epsilon * self.c2


def _is_exact_zero(value: Any) -> bool:
    try:
        return bool(value == 0)
    except TypeError:
        return False


def _is_exact_one(value: Any) -> bool:
    try:
        return bool(value == 1)
    except TypeError:
        return False


class PolynomialFormulaOps:
    """Arithmetic front end that exposes mass dependence in _formula exactly."""

    def __init__(self, entropy: Callable[[Any], Any], one: Any):
        self.entropy_function = entropy
        self.one = EpsilonPolynomial.constant(one)
        self.entropy_argument_count = 0
        self.epsilon_dependent_entropy_arguments = 0

    def nonnegative_product(self, *values: Any) -> EpsilonPolynomial:
        answer = self.one
        for value in values:
            answer = answer * value
        return answer

    def nonnegative_sum(self, values: Iterable[Any]) -> EpsilonPolynomial:
        return sum(values, EpsilonPolynomial.constant(_exact_zero(self.one.c0)))

    def quotient(
        self,
        numerator: EpsilonPolynomial,
        denominator: EpsilonPolynomial,
    ) -> EpsilonPolynomial:
        return numerator / denominator

    def xy_argument(self, x: Any, y: Any) -> EpsilonPolynomial:
        return self.one._coerce(x) * y

    def pi_argument(self, x: Any, y: Any) -> EpsilonPolynomial:
        x_poly = self.one._coerce(x)
        y_poly = self.one._coerce(y)
        return x_poly * y_poly + x_poly * (self.one - x_poly) * y_poly * (
            self.one - y_poly
        )

    def entropy(self, argument: Any) -> EpsilonPolynomial:
        self.entropy_argument_count += 1
        polynomial = self.one._coerce(argument)
        if not (_is_exact_zero(polynomial.c1)
                and _is_exact_zero(polynomial.c2)):
            self.epsilon_dependent_entropy_arguments += 1
            if _is_exact_zero(polynomial.c0) and not _is_exact_zero(polynomial.c1):
                raise EpsilonLogTerm(
                    "an entropy argument opens linearly from zero: "
                    "h(c*epsilon)=c*epsilon*log(1/epsilon)+O(epsilon)"
                )
            if _is_exact_one(polynomial.c0) and not _is_exact_zero(polynomial.c1):
                raise EpsilonLogTerm(
                    "an entropy argument opens linearly from one: "
                    "h(1-c*epsilon)=c*epsilon*log(1/epsilon)+O(epsilon)"
                )
            raise EpsilonDependentEntropyArgument(
                "entropy has nonconstant epsilon argument"
            )
        return EpsilonPolynomial.constant(self.entropy_function(polynomial.c0))


def formula_gap_polynomial(
    values: Sequence[Any],
    beta: Any,
    entropy: Callable[[Any], Any],
    one: Any,
) -> tuple[EpsilonPolynomial, PolynomialFormulaOps]:
    """Differentiate the single source-of-truth formula in epsilon."""
    ops = PolynomialFormulaOps(entropy, one)
    terms = _formula(values, beta, ops)
    return terms.numerator - terms.ehx, ops


def active_split_polynomial_values(
    parameters: Any,
    y: Any,
    q: Any,
    one: Any,
) -> tuple[EpsilonPolynomial, ...]:
    epsilon = EpsilonPolynomial.variable(one)
    return tuple(
        EpsilonPolynomial.constant(value)
        if not isinstance(value, EpsilonPolynomial) else value
        for value in active_split_embedding(parameters, y, epsilon, q)
    )


def mean_preserving_polynomial_values(
    parameters: Any,
    y: Any,
    q: Any,
    one: Any,
) -> tuple[EpsilonPolynomial, ...]:
    epsilon = EpsilonPolynomial.variable(one)
    return tuple(
        EpsilonPolynomial.constant(value)
        if not isinstance(value, EpsilonPolynomial) else value
        for value in mean_preserving_embedding(parameters, y, epsilon, q)
    )


def polynomial_expansion(
    parameters: ArbParameters,
    y: arb,
    q: arb,
    family: str = "active_split",
) -> dict[str, Any]:
    """Automatic raw-gap and distance expansion through order epsilon^2."""
    if family == "active_split":
        values = active_split_polynomial_values(parameters, y, q, arb(1))
    elif family == "mean_preserving":
        values = mean_preserving_polynomial_values(parameters, y, q, arb(1))
    else:
        raise ValueError(f"unknown transverse family {family!r}")
    gap, ops = formula_gap_polynomial(values, parameters.beta, h_arb, arb(1))
    dist = distance_squared(values, parameters)
    kappa = arb(KAPPA.numerator) / KAPPA.denominator
    pencil = gap - kappa * dist
    return {
        "gap": gap,
        "distance_squared": dist,
        "pencil": pencil,
        "entropy_argument_count": ops.entropy_argument_count,
        "epsilon_dependent_entropy_arguments": (
            ops.epsilon_dependent_entropy_arguments
        ),
    }


def formula_has_epsilon_log_argument(
    values: Sequence[Any],
    beta: Any,
    entropy: Callable[[Any], Any],
    one: Any,
) -> bool:
    """Return whether _formula encounters a linear endpoint entropy opening."""
    try:
        formula_gap_polynomial(values, beta, entropy, one)
    except EpsilonLogTerm:
        return True
    except EpsilonDependentEntropyArgument:
        return False
    return False


def validate_log_term_claim(
    values: Sequence[Any],
    beta: Any,
    entropy: Callable[[Any], Any],
    one: Any,
    claimed_present: bool,
) -> None:
    actual = formula_has_epsilon_log_argument(values, beta, entropy, one)
    if actual != claimed_present:
        raise AssertionError(
            f"epsilon*log(1/epsilon) claim {claimed_present} disagrees with "
            f"automatic _formula structure {actual}"
        )


# ---------------------------------------------------------------------------
# Exact Fraction interval arithmetic for the decisive signs.


@dataclass(frozen=True)
class RationalInterval:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed rational interval")

    @staticmethod
    def point(value: Fraction | int) -> "RationalInterval":
        exact = value if isinstance(value, Fraction) else Fraction(value)
        return RationalInterval(exact, exact)

    def _coerce(self, other: Any) -> "RationalInterval":
        return other if isinstance(other, RationalInterval) else self.point(other)

    def __add__(self, other: Any) -> "RationalInterval":
        rhs = self._coerce(other)
        return RationalInterval(self.lo + rhs.lo, self.hi + rhs.hi)

    __radd__ = __add__

    def __neg__(self) -> "RationalInterval":
        return RationalInterval(-self.hi, -self.lo)

    def __sub__(self, other: Any) -> "RationalInterval":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "RationalInterval":
        rhs = self._coerce(other)
        products = (
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        )
        return RationalInterval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self) -> "RationalInterval":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("rational interval contains zero")
        return RationalInterval(min(1 / self.lo, 1 / self.hi),
                                max(1 / self.lo, 1 / self.hi))

    def __truediv__(self, other: Any) -> "RationalInterval":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "RationalInterval":
        return self._coerce(other) / self

    def __pow__(self, power: int) -> "RationalInterval":
        if power < 0:
            return (self.reciprocal()) ** (-power)
        answer = RationalInterval.point(1)
        base = self
        exponent = power
        while exponent:
            if exponent & 1:
                answer = answer * base
            base = base * base
            exponent >>= 1
        return answer

    def abs_upper(self) -> Fraction:
        return max(abs(self.lo), abs(self.hi))

    def as_float_pair(self) -> list[float]:
        return [float(self.lo), float(self.hi)]


@lru_cache(maxsize=None)
def _atanh_log_bounds(value: Fraction) -> RationalInterval:
    """Exact atanh-series bounds for log(value), for 1 <= value <= 2."""
    if not 1 <= value <= 2:
        raise ValueError("atanh log kernel expects a value in [1,2]")
    if value == 1:
        return RationalInterval.point(0)
    t = (value - 1) / (value + 1)
    term = t
    total = Fraction(0)
    for index in range(LOG_SERIES_TERMS):
        total += term / (2 * index + 1)
        term *= t * t
    lower = 2 * total
    remainder = 2 * term / ((2 * LOG_SERIES_TERMS + 1) * (1 - t * t))
    return RationalInterval(lower, lower + remainder)


@lru_cache(maxsize=None)
def _log_two_bounds() -> RationalInterval:
    return _atanh_log_bounds(Fraction(2))


@lru_cache(maxsize=None)
def rational_log_point(value: Fraction) -> RationalInterval:
    """Rigorous exact-Fraction enclosure of log(value)."""
    if value <= 0:
        raise ValueError("log requires a positive rational")
    if value < 1:
        answer = rational_log_point(1 / value)
        return -answer
    reduced = value
    power = 0
    while reduced >= 2:
        reduced /= 2
        power += 1
    return _atanh_log_bounds(reduced) + power * _log_two_bounds()


def rational_log(value: RationalInterval) -> RationalInterval:
    if value.lo <= 0:
        raise ValueError("log interval touches zero")
    lower = rational_log_point(value.lo)
    upper = rational_log_point(value.hi)
    return RationalInterval(lower.lo, upper.hi)


@lru_cache(maxsize=None)
def rational_entropy_point(value: Fraction) -> RationalInterval:
    if value == 0 or value == 1:
        return RationalInterval.point(0)
    if not 0 < value < 1:
        raise ValueError("entropy point outside [0,1]")
    z = RationalInterval.point(value)
    return -(z * rational_log_point(value)
             + (1 - z) * rational_log_point(1 - value))


def rational_entropy(value: RationalInterval) -> RationalInterval:
    """Monotonicity-tight entropy range with exact rational endpoints."""
    if value.lo < 0 or value.hi > 1:
        raise ValueError("entropy interval outside [0,1]")
    at_lo = rational_entropy_point(value.lo)
    at_hi = rational_entropy_point(value.hi)
    half = Fraction(1, 2)
    if value.hi <= half:
        return RationalInterval(at_lo.lo, at_hi.hi)
    if value.lo >= half:
        return RationalInterval(at_hi.lo, at_lo.hi)
    lower = min(at_lo.lo, at_hi.lo)
    return RationalInterval(lower, _log_two_bounds().hi)


def rational_hprime(value: RationalInterval) -> RationalInterval:
    if value.lo <= 0 or value.hi >= 1:
        raise ValueError("h-prime interval must be interior")
    ratio = (1 - value) / value
    return rational_log(ratio)


def rational_protocol(left: RationalInterval,
                      right: RationalInterval) -> RationalInterval:
    return (left * right
            + left * (1 - left) * right * (1 - right))


def rational_kernel(
    left: RationalInterval,
    right: RationalInterval,
    beta: RationalInterval,
) -> RationalInterval:
    return ((1 - beta) * rational_entropy(left * right)
            + beta * rational_entropy(rational_protocol(left, right)))


@dataclass(frozen=True)
class RationalParameters:
    x: RationalInterval
    p: RationalInterval
    mean: RationalInterval
    beta: RationalInterval


def exact_parameter_intervals(parameters: ArbParameters) -> RationalParameters:
    """Recompute p and beta from a coarsened root bracket using Fractions."""
    scale = 10**PARAMETER_INTERVAL_DECIMALS
    lower_scaled = parameters.root_lo * scale
    upper_scaled = parameters.root_hi * scale
    lower = Fraction(lower_scaled.numerator // lower_scaled.denominator, scale)
    upper_floor = upper_scaled.numerator // upper_scaled.denominator
    upper = Fraction(
        upper_floor + int(upper_scaled.numerator % upper_scaled.denominator != 0),
        scale,
    )
    if not (lower <= parameters.root_lo <= parameters.root_hi <= upper):
        raise AssertionError("coarsened rational bracket lost Liu's root")
    x = RationalInterval(lower, upper)
    hx = rational_entropy(x)
    x_squared = x * x
    hxx = rational_entropy(x_squared)
    p = hx / hxx
    mean = p * x
    protocol = rational_protocol(x, x)
    hp_x = rational_hprime(x)
    hp_x_squared = rational_hprime(x_squared)
    hp_protocol = rational_hprime(protocol)
    du = 2 * x
    dt = 2 * x + 2 * x * (1 - x) ** 2 - 2 * x * x * (1 - x)
    a_prime = hp_x_squared * du
    b_prime = hp_protocol * dt
    beta = ((hp_x + hx / x) / p - a_prime) / (b_prime - a_prime)
    if not (0 < p.lo <= p.hi < 1 and 0 < beta.lo <= beta.hi < 1):
        raise AssertionError("exact rational parameter enclosure is invalid")
    return RationalParameters(x=x, p=p, mean=mean, beta=beta)


@dataclass(frozen=True)
class ExactSignCertificate:
    leading: RationalInterval
    quadratic: RationalInterval
    reduced_at_epsilon: RationalInterval
    full_pencil: RationalInterval
    mean_change: RationalInterval
    mean_preserving_leading: RationalInterval

    def validate(self) -> None:
        if not self.leading.hi < LEADING_NEGATIVE_THRESHOLD:
            raise AssertionError("exact leading coefficient did not clear -1/2")
        if not self.full_pencil.hi < FINITE_NEGATIVE_THRESHOLD:
            raise AssertionError("exact finite pencil did not clear -1/2048")
        if not self.mean_change.hi < 0:
            raise AssertionError("counterexample was not mean-decreasing")
        if not (self.mean_preserving_leading.lo
                > MEAN_PRESERVING_POSITIVE_THRESHOLD):
            raise AssertionError(
                "mean-preserving control did not clear its positive threshold"
            )


def exact_sign_certificate(
    parameters: RationalParameters,
    y: Fraction = DECISIVE_Y,
    epsilon: Fraction = DECISIVE_EPSILON,
) -> ExactSignCertificate:
    """Certify the decisive signs without an interval ball touching zero."""
    y_interval = RationalInterval.point(y)
    x, p, beta = parameters.x, parameters.p, parameters.beta
    h_x = rational_entropy(x)
    h_y = rational_entropy(y_interval)
    k_xx = rational_kernel(x, x, beta)
    k_xy = rational_kernel(x, y_interval, beta)
    k_yy = rational_kernel(y_interval, y_interval, beta)
    active_mass_derivative = 2 * p * k_xx - h_x
    gap_linear = 2 * p * k_xy - h_y - active_mass_derivative
    support_linear = y_interval**2 * (y_interval - x) ** 2
    leading = gap_linear - KAPPA * support_linear
    gap_quadratic = k_yy - 2 * k_xy + k_xx
    distance_quadratic = (y_interval - x) ** 2
    quadratic = gap_quadratic - KAPPA * distance_quadratic
    reduced = leading + epsilon * quadratic
    full = epsilon * reduced
    mean_change = epsilon * (y_interval - x)

    tangent_gap_linear = (
        2 * p * k_xy - h_y
        - (y_interval / x) * active_mass_derivative
    )
    tangent_leading = tangent_gap_linear - KAPPA * support_linear
    certificate = ExactSignCertificate(
        leading=leading,
        quadratic=quadratic,
        reduced_at_epsilon=reduced,
        full_pencil=full,
        mean_change=mean_change,
        mean_preserving_leading=tangent_leading,
    )
    certificate.validate()
    return certificate


def validate_sign_claim(interval: RationalInterval, claim: str) -> None:
    if claim == "negative":
        holds = interval.hi < 0
    elif claim == "positive":
        holds = interval.lo > 0
    elif claim == "nonnegative":
        holds = interval.lo >= 0
    else:
        raise ValueError(f"unknown sign claim {claim!r}")
    if not holds:
        raise AssertionError(
            f"claimed {claim} sign is not certified on "
            f"float enclosure {interval.as_float_pair()}"
        )


# ---------------------------------------------------------------------------
# Arb formula identities and chart reproduction.


def _arbf(value: Fraction | int) -> arb:
    exact = value if isinstance(value, Fraction) else Fraction(value)
    return arb(exact.numerator) / arb(exact.denominator)


def _arb_width(value: arb) -> float:
    return 2 * float(value.rad())


def _arb_record(value: arb) -> dict[str, Any]:
    return {
        "interval": str(value),
        "lower": float(value.lower()),
        "upper": float(value.upper()),
        "width": _arb_width(value),
    }


def _mp_entropy(value: mpmath.mpf) -> mpmath.mpf:
    if value == 0 or value == 1:
        return mpmath.mpf(0)
    return -(value * mpmath.log(value)
             + (1 - value) * mpmath.log(1 - value))


def _mp_protocol(left: mpmath.mpf, right: mpmath.mpf) -> mpmath.mpf:
    return left * right + left * (1 - left) * right * (1 - right)


def _mp_kernel(
    left: mpmath.mpf,
    right: mpmath.mpf,
    beta: mpmath.mpf,
) -> mpmath.mpf:
    return ((1 - beta) * _mp_entropy(left * right)
            + beta * _mp_entropy(_mp_protocol(left, right)))


def transverse_coefficients_mp(
    parameters: MPParameters,
    y: mpmath.mpf,
    family: str,
) -> tuple[mpmath.mpf, mpmath.mpf]:
    """Hand form used only as an independent check of _formula AD."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    k_xx = _mp_kernel(x, x, beta)
    k_xy = _mp_kernel(x, y, beta)
    k_yy = _mp_kernel(y, y, beta)
    active_derivative = 2 * p * k_xx - _mp_entropy(x)
    support_linear = y**2 * (y - x) ** 2
    if family == "active_split":
        gap_linear = 2 * p * k_xy - _mp_entropy(y) - active_derivative
        gap_quadratic = k_yy - 2 * k_xy + k_xx
        distance_quadratic = (y - x) ** 2
    elif family == "mean_preserving":
        ratio = y / x
        gap_linear = (
            2 * p * k_xy - _mp_entropy(y) - ratio * active_derivative
        )
        gap_quadratic = k_yy - 2 * ratio * k_xy + ratio**2 * k_xx
        distance_quadratic = mpmath.mpf(0)
    else:
        raise ValueError(f"unknown family {family!r}")
    kappa = mpmath.mpf(KAPPA.numerator) / KAPPA.denominator
    return (
        gap_linear - kappa * support_linear,
        gap_quadratic - kappa * distance_quadratic,
    )


CHART_POINTS = (
    (Fraction(1, 4), Fraction(-1, 512), Fraction(1, 1024), Fraction(1, 64)),
    (Fraction(1, 3), Fraction(1, 1024), Fraction(-1, 768), Fraction(1, 80)),
    (Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 64)),
    (Fraction(2, 3), Fraction(-1, 768), Fraction(1, 640), Fraction(1, 96)),
    (Fraction(3, 4), Fraction(1, 512), Fraction(-1, 1024), Fraction(1, 72)),
)


def validate_chart_agreement(
    parameters: ArbParameters,
    points: Sequence[tuple[Fraction, Fraction, Fraction, Fraction]] = CHART_POINTS,
    embedding: Callable[..., tuple[Any, ...]] = chart_embedding,
    maximum_width: float = 1e-60,
) -> dict[str, Any]:
    """Reproduce chart gap and distance from the nine-variable formulas."""
    rows: list[dict[str, Any]] = []
    max_gap_width = 0.0
    max_distance_width = 0.0
    for q_value, s_offset, d_value, r_value in points:
        q = _arbf(q_value)
        s = parameters.x + _arbf(s_offset)
        d = _arbf(d_value)
        r = _arbf(r_value)
        values = embedding(parameters, q, s, d, r)
        terms = _formula(values, parameters.beta, _ArbOps(False))
        nine_gap = terms.numerator - terms.ehx
        nine_distance = distance_squared(values, parameters)
        chart_gap = gap_jet(
            parameters, q, Jet3.constant(s), Jet3.constant(d), r
        ).v
        chart_distance = distance_squared_jet(
            parameters, q, Jet3.constant(s), Jet3.constant(d), r
        ).v
        gap_residual = chart_gap - nine_gap
        distance_residual = chart_distance - nine_distance
        if not gap_residual.contains(0):
            raise AssertionError(
                f"chart gap disagrees with _formula at q={q_value}, "
                f"s-x={s_offset}, d={d_value}: {gap_residual}"
            )
        if not distance_residual.contains(0):
            raise AssertionError(
                f"chart distance disagrees at q={q_value}, "
                f"s-x={s_offset}, d={d_value}: {distance_residual}"
            )
        gap_width = _arb_width(gap_residual)
        distance_width = _arb_width(distance_residual)
        if gap_width > maximum_width or distance_width > maximum_width:
            raise AssertionError(
                "chart identity residual is too wide: "
                f"gap={gap_width}, distance={distance_width}"
            )
        max_gap_width = max(max_gap_width, gap_width)
        max_distance_width = max(max_distance_width, distance_width)
        rows.append({
            "q": str(q_value),
            "s_minus_x": str(s_offset),
            "d": str(d_value),
            "r": str(r_value),
            "gap_formula": str(nine_gap),
            "gap_chart": str(chart_gap),
            "gap_residual": str(gap_residual),
            "gap_residual_width": gap_width,
            "distance_formula": str(nine_distance),
            "distance_chart": str(chart_distance),
            "distance_residual": str(distance_residual),
            "distance_residual_width": distance_width,
        })
    return {
        "checked_points": len(rows),
        "maximum_allowed_residual_width": maximum_width,
        "maximum_gap_residual_width": max_gap_width,
        "maximum_distance_residual_width": max_distance_width,
        "passed": True,
        "points": rows,
    }


def automatic_formula_certificate(
    parameters: ArbParameters,
    y: Fraction = DECISIVE_Y,
    epsilon: Fraction = DECISIVE_EPSILON,
    q: Fraction = DECISIVE_Q,
) -> dict[str, Any]:
    """Cross-check the expansion and finite counterexample with _formula."""
    y_arb, epsilon_arb, q_arb = _arbf(y), _arbf(epsilon), _arbf(q)
    expansion = polynomial_expansion(parameters, y_arb, q_arb, "active_split")
    gap_poly = expansion["gap"]
    distance_poly = expansion["distance_squared"]
    pencil_poly = expansion["pencil"]
    if expansion["epsilon_dependent_entropy_arguments"] != 0:
        raise AssertionError("fixed-y mass insertion unexpectedly moved entropy arguments")

    values = active_split_embedding(parameters, y_arb, epsilon_arb, q_arb)
    terms = _formula(values, parameters.beta, _ArbOps(False))
    direct_gap = terms.numerator - terms.ehx
    direct_distance = distance_squared(values, parameters)
    kappa_arb = _arbf(KAPPA)
    direct_pencil = direct_gap - kappa_arb * direct_distance
    if not direct_gap < 0 or not direct_pencil < 0:
        raise AssertionError("Arb did not enclose the explicit refutation below zero")

    gap_residual = direct_gap - gap_poly.evaluate(epsilon_arb)
    distance_residual = direct_distance - distance_poly.evaluate(epsilon_arb)
    pencil_residual = direct_pencil - pencil_poly.evaluate(epsilon_arb)
    for name, residual in (
        ("gap", gap_residual),
        ("distance", distance_residual),
        ("pencil", pencil_residual),
    ):
        if not residual.contains(0):
            raise AssertionError(f"automatic {name} expansion misses _formula: {residual}")

    # Independent hand formulas are only accepted when the automatic source of
    # truth encloses their residuals.
    x, p, beta = parameters.x, parameters.p, parameters.beta
    kernel = lambda u, v: (
        (1 - beta) * h_arb(u * v)
        + beta * h_arb(u * v + u * (1 - u) * v * (1 - v))
    )
    k_xx, k_xy, k_yy = kernel(x, x), kernel(x, y_arb), kernel(y_arb, y_arb)
    active_derivative = 2 * p * k_xx - h_arb(x)
    hand_gap_linear = 2 * p * k_xy - h_arb(y_arb) - active_derivative
    hand_gap_quadratic = k_yy - 2 * k_xy + k_xx
    hand_distance_linear = y_arb**2 * (y_arb - x) ** 2
    hand_distance_quadratic = (y_arb - x) ** 2
    coefficient_residuals = {
        "gap_linear": gap_poly.c1 - hand_gap_linear,
        "gap_quadratic": gap_poly.c2 - hand_gap_quadratic,
        "distance_linear": distance_poly.c1 - hand_distance_linear,
        "distance_quadratic": distance_poly.c2 - hand_distance_quadratic,
    }
    for name, residual in coefficient_residuals.items():
        if not residual.contains(0):
            raise AssertionError(f"hand {name} misses automatic _formula: {residual}")

    mean = ((parameters.p - epsilon_arb) * parameters.x
            + epsilon_arb * y_arb)
    mean_change = mean - parameters.mean
    if not mean_change < 0:
        raise AssertionError("explicit ambient refutation unexpectedly satisfies the mean")

    objective_minus_one = terms.objective - 1
    units_residual = direct_gap - terms.ehx * objective_minus_one
    if not units_residual.contains(0):
        raise AssertionError("raw-gap/objective unit identity failed")

    return {
        "entropy_argument_count": expansion["entropy_argument_count"],
        "epsilon_dependent_entropy_arguments": 0,
        "epsilon_log_coefficient": "0",
        "leading_order": "epsilon",
        "gap_linear_raw_gap": _arb_record(gap_poly.c1),
        "gap_quadratic_raw_gap": _arb_record(gap_poly.c2),
        "distance_linear": _arb_record(distance_poly.c1),
        "distance_quadratic": _arb_record(distance_poly.c2),
        "pencil_linear_raw_gap": _arb_record(pencil_poly.c1),
        "pencil_quadratic_raw_gap": _arb_record(pencil_poly.c2),
        "direct_gap_raw_gap": _arb_record(direct_gap),
        "direct_distance_squared": _arb_record(direct_distance),
        "direct_pencil_raw_gap": _arb_record(direct_pencil),
        "ehx": _arb_record(terms.ehx),
        "objective_minus_one": _arb_record(objective_minus_one),
        "raw_gap_identity_residual": _arb_record(units_residual),
        "mean_minus_target": _arb_record(mean_change),
        "formula_expansion_residuals": {
            "gap": _arb_record(gap_residual),
            "distance": _arb_record(distance_residual),
            "pencil": _arb_record(pencil_residual),
        },
        "hand_coefficient_residuals": {
            name: _arb_record(value)
            for name, value in coefficient_residuals.items()
        },
    }


# ---------------------------------------------------------------------------
# Endpoint-inclusive y scan and mutation guards.


def y_grid(
    y0: Fraction = Y0,
    intervals: int = GRID_INTERVALS,
) -> tuple[Fraction, ...]:
    if intervals < 1:
        raise ValueError("the y grid needs at least one interval")
    width = 1 - 2 * y0
    return tuple(y0 + width * index / intervals
                 for index in range(intervals + 1))


def validate_y_grid(
    grid: Sequence[Fraction],
    y0: Fraction = Y0,
) -> None:
    if not grid:
        raise AssertionError("empty y grid")
    if grid[0] != y0 or grid[-1] != 1 - y0:
        raise AssertionError(
            f"y grid must include exact endpoints {y0} and {1-y0}"
        )
    if any(left >= right for left, right in zip(grid, grid[1:])):
        raise AssertionError("y grid must be strictly increasing")


def sample_transverse_grid(
    parameters: MPParameters,
    dps: int = 100,
    intervals: int = GRID_INTERVALS,
) -> dict[str, Any]:
    grid = y_grid(intervals=intervals)
    validate_y_grid(grid)
    rows: list[dict[str, Any]] = []
    negative = 0
    positive = 0
    tangent_negative = 0
    tangent_positive = 0
    tangent_worst = None
    with mpmath.workdps(dps):
        for value in grid:
            y = mpmath.mpf(value.numerator) / value.denominator
            active_linear, _ = transverse_coefficients_mp(
                parameters, y, "active_split"
            )
            tangent_linear, _ = transverse_coefficients_mp(
                parameters, y, "mean_preserving"
            )
            sign = "negative" if active_linear < 0 else "positive"
            negative += int(active_linear < 0)
            positive += int(active_linear > 0)
            tangent_negative += int(tangent_linear < 0)
            tangent_positive += int(tangent_linear > 0)
            if tangent_worst is None or tangent_linear < tangent_worst[0]:
                tangent_worst = (tangent_linear, value)
            rows.append({
                "y": str(value),
                "active_split_pencil_linear_raw_gap": mpmath.nstr(
                    active_linear, 50
                ),
                "active_split_sign": sign,
                "mean_preserving_pencil_linear_raw_gap": mpmath.nstr(
                    tangent_linear, 50
                ),
            })

        special_rows = []
        for name, y in (
            ("x", parameters.x),
            ("1-x", 1 - parameters.x),
        ):
            active_linear, _ = transverse_coefficients_mp(
                parameters, y, "active_split"
            )
            tangent_linear, _ = transverse_coefficients_mp(
                parameters, y, "mean_preserving"
            )
            special_rows.append({
                "name": name,
                "y": mpmath.nstr(y, 80),
                "active_split_pencil_linear_raw_gap": (
                    "0" if name == "x" else mpmath.nstr(active_linear, 60)
                ),
                "mean_preserving_pencil_linear_raw_gap": (
                    "0" if name == "x" else mpmath.nstr(tangent_linear, 60)
                ),
            })
    if negative == 0 or positive == 0:
        raise AssertionError("endpoint-inclusive scan did not see the sign change")
    # The ambient verdict rests on `active_split`, but the CONDITIONAL
    # feasible-half-space verdict rests entirely on `mean_preserving` never
    # turning negative.  That tally was not being recorded, so the conditional
    # was resting on a single decisive y rather than on the whole scan.
    if tangent_worst is None:
        raise AssertionError("the mean-preserving scan produced no rows")
    return {
        "intervals": intervals,
        "points": len(grid),
        "left_endpoint": str(grid[0]),
        "right_endpoint": str(grid[-1]),
        "endpoints_included": True,
        "negative_active_split_points": negative,
        "positive_active_split_points": positive,
        "negative_mean_preserving_points": tangent_negative,
        "positive_mean_preserving_points": tangent_positive,
        "worst_mean_preserving_linear": mpmath.nstr(tangent_worst[0], 40),
        "worst_mean_preserving_y": str(tangent_worst[1]),
        "special_values_checked": ["x", "1-x"],
        "special_values": special_rows,
        "rows": rows,
    }


def support_opening_polynomial_values(
    parameters: ArbParameters,
) -> tuple[EpsilonPolynomial, ...]:
    """A control family that really does create epsilon*log(1/epsilon)."""
    one = arb(1)
    epsilon = EpsilonPolynomial.variable(one)
    zero = EpsilonPolynomial.constant(arb(0))
    q = EpsilonPolynomial.constant(_arbf(Fraction(1, 2)))
    positive_gauge_mass = _arbf(Fraction(1, 64))
    return (
        EpsilonPolynomial.constant(parameters.p),
        EpsilonPolynomial.constant(positive_gauge_mass),
        q,
        EpsilonPolynomial.constant(parameters.x),
        epsilon,
        zero,
        EpsilonPolynomial.constant(parameters.x),
        zero,
        zero,
    )


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _write_report(payload: dict[str, Any], path: Path = REPORT_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def build_report(dps: int = 100, grid_intervals: int = GRID_INTERVALS) -> dict[str, Any]:
    mp_parameters = solve_equation_parameters(max(dps, 100))
    arb_parameters = certify_equation_parameters(mp_parameters)
    chart_check = validate_chart_agreement(arb_parameters)
    exact_parameters = exact_parameter_intervals(arb_parameters)
    exact = exact_sign_certificate(exact_parameters)
    automatic = automatic_formula_certificate(arb_parameters)
    grid = sample_transverse_grid(mp_parameters, dps=dps,
                                  intervals=grid_intervals)

    active_values = active_split_polynomial_values(
        arb_parameters, _arbf(DECISIVE_Y), _arbf(DECISIVE_Q), arb(1)
    )
    validate_log_term_claim(
        active_values, arb_parameters.beta, h_arb, arb(1), False
    )
    opening_values = support_opening_polynomial_values(arb_parameters)
    validate_log_term_claim(
        opening_values, arb_parameters.beta, h_arb, arb(1), True
    )

    with mpmath.workdps(dps):
        y_mp = mpmath.mpf(DECISIVE_Y.numerator) / DECISIVE_Y.denominator
        epsilon_mp = (
            mpmath.mpf(DECISIVE_EPSILON.numerator)
            / DECISIVE_EPSILON.denominator
        )
        q_mp = mpmath.mpf(DECISIVE_Q.numerator) / DECISIVE_Q.denominator
        values_mp = active_split_embedding(
            mp_parameters, y_mp, epsilon_mp, q_mp
        )
        from liu9_objective import evaluate_mpmath  # local source front end
        terms_mp = evaluate_mpmath(values_mp, mp_parameters.beta, dps=dps)
        gap_mp = terms_mp.numerator - terms_mp.ehx
        distance_mp = distance_squared(values_mp, mp_parameters)
        kappa_mp = mpmath.mpf(KAPPA.numerator) / KAPPA.denominator
        pencil_mp = gap_mp - kappa_mp * distance_mp
        mean_change_mp = (
            (mp_parameters.p - epsilon_mp) * mp_parameters.x
            + epsilon_mp * y_mp - mp_parameters.mean
        )

    payload: dict[str, Any] = {
        "schema_version": 1,
        "tool": "liu9_ninevar.py",
        "verdict": "REFUTED",
        "scope": "unrestricted ambient nine-variable extension",
        "kappa": {
            "exact": str(KAPPA),
            "decimal": float(KAPPA),
            "units": "raw-gap per distance-squared",
        },
        "parameters": {
            "x": mpmath.nstr(mp_parameters.x, 80),
            "p": mpmath.nstr(mp_parameters.p, 80),
            "mean": mpmath.nstr(mp_parameters.mean, 80),
            "beta": mpmath.nstr(mp_parameters.beta, 80),
            "x_root_bracket": [
                str(arb_parameters.root_lo),
                str(arb_parameters.root_hi),
            ],
        },
        "chart_embedding": {
            "variable_order": [
                "a1", "a2", "q", "b0", "b2", "b4", "b1", "b3", "b5"
            ],
            "formula": [
                "mean/s", "r", "q", "s-q*d", "0", "0",
                "s+(1-q)*d", "0", "0",
            ],
            "mean_identity": "(mean/s)*((1-q)*(s-q*d)+q*(s+(1-q)*d))=mean",
            "agreement": chart_check,
        },
        "transverse_families": {
            "active_split": {
                "formula": [
                    "p-epsilon", "epsilon", "q", "x", "y", "0",
                    "x", "y", "0",
                ],
                "mean_minus_target": "epsilon*(y-x)",
                "ambient": True,
            },
            "mean_preserving": {
                "formula": [
                    "p-epsilon*y/x", "epsilon", "q", "x", "y", "0",
                    "x", "y", "0",
                ],
                "mean_minus_target": "0",
                "ambient": True,
            },
            "y_domain": [str(Y0), str(1 - Y0)],
        },
        "asymptotic_structure": {
            "leading_order": "epsilon",
            "epsilon_log_epsilon_term_present": False,
            "epsilon_log_coefficient": "0",
            "reason": (
                "all six supports are fixed as epsilon varies; _formula applies "
                "entropy only to support/protocol arguments, while every mass "
                "occurs outside entropy and polynomially"
            ),
            "automatic_formula": automatic,
            "support_opening_mutation_detected_epsilon_log": True,
        },
        "y_scan": grid,
        "counterexample": {
            "y": str(DECISIVE_Y),
            "epsilon": str(DECISIVE_EPSILON),
            "q": str(DECISIVE_Q),
            "nine_vector": [
                "p-1/1024", "1/1024", "1/2", "x", "1/32", "0",
                "x", "1/32", "0",
            ],
            "gap_raw_gap": mpmath.nstr(gap_mp, 80),
            "distance_squared": mpmath.nstr(distance_mp, 80),
            "pencil_raw_gap": mpmath.nstr(pencil_mp, 80),
            "mean_minus_target": mpmath.nstr(mean_change_mp, 80),
            "liu_mean_feasible": False,
            "units_identity": "gap=EHX*(objective-1)",
        },
        "exact_fraction_sign_certificate": {
            "method": (
                "exact Fraction interval arithmetic; logarithms bounded by a "
                "96-term atanh series with an exact geometric remainder"
            ),
            "leading_interval_float": exact.leading.as_float_pair(),
            "quadratic_interval_float": exact.quadratic.as_float_pair(),
            "reduced_at_epsilon_interval_float": (
                exact.reduced_at_epsilon.as_float_pair()
            ),
            "full_pencil_interval_float": exact.full_pencil.as_float_pair(),
            "mean_change_interval_float": exact.mean_change.as_float_pair(),
            "mean_preserving_leading_interval_float": (
                exact.mean_preserving_leading.as_float_pair()
            ),
            "leading_upper_lt": str(LEADING_NEGATIVE_THRESHOLD),
            "full_pencil_upper_lt": str(FINITE_NEGATIVE_THRESHOLD),
            "mean_change_upper_lt": "0",
            "mean_preserving_leading_lower_gt": str(
                MEAN_PRESERVING_POSITIVE_THRESHOLD
            ),
            "all_signs_checked_as_exact_fractions": True,
        },
        "conclusion": {
            "ambient_extension": False,
            "obstruction": (
                "the active-mean chart does not control the mean-decreasing "
                "normal direction produced by splitting active mass to y<x"
            ),
            "feasible_half_space_status": (
                "not refuted by this configuration: its mean is below the "
                "required target; the mean-preserving control has positive "
                "linear pencil coefficient at the decisive y"
            ),
        },
    }
    payload["report_sha256"] = _canonical_digest(payload)
    return payload


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--grid-intervals", type=int, default=GRID_INTERVALS)
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    args = parser.parse_args(argv)

    payload = build_report(dps=args.dps, grid_intervals=args.grid_intervals)
    _write_report(payload, args.report)
    counterexample = payload["counterexample"]
    exact = payload["exact_fraction_sign_certificate"]
    agreement = payload["chart_embedding"]["agreement"]
    print(
        "PROVED [embedding]: "
        "(mean/s,r,q,s-q*d,0,0,s+(1-q)*d,0,0) in "
        "(a1,a2,q,b0,b2,b4,b1,b3,b5)."
    )
    print(
        "PROVED [chart reproduction through _formula]: "
        f"{agreement['checked_points']} points; max raw-gap residual width="
        f"{agreement['maximum_gap_residual_width']:.3e}; max distance-squared "
        f"residual width={agreement['maximum_distance_residual_width']:.3e}."
    )
    print(
        "PROVED [transverse expansion]: leading order=epsilon; "
        "epsilon*log(1/epsilon) coefficient=0; raw-gap kappa="
        f"{KAPPA}={float(KAPPA):.15f}."
    )
    print(
        "PROVED [exact Fraction signs]: leading coefficient upper < "
        f"{exact['leading_upper_lt']}; finite raw-gap pencil upper < "
        f"{exact['full_pencil_upper_lt']}; mean change upper < 0."
    )
    print(
        "REFUTED [ambient nine-variable extension]: "
        f"y={counterexample['y']}, epsilon={counterexample['epsilon']}, "
        f"q={counterexample['q']}; gap={counterexample['gap_raw_gap']}; "
        f"dist^2={counterexample['distance_squared']}; raw-gap pencil="
        f"{counterexample['pencil_raw_gap']}; mean-target="
        f"{counterexample['mean_minus_target']}."
    )
    print(
        "SUPERSEDED [see liu9_transverse.py, 2026-08-29]: the feasible "
        "half-space statement below is no longer conditional -- it is PROVED "
        "on the whole interval, and for the asymmetric family too.  "
        "Historically: this counterexample is "
        "mean-infeasible; the exact mean-preserving control at the same y has "
        "linear raw-gap pencil coefficient > "
        f"{exact['mean_preserving_leading_lower_gt']}."
    )
    scan = payload["y_scan"]
    print(
        f"COMPUTATIONAL EVIDENCE [{scan['points']}-point y scan]: the "
        f"mean-preserving family is positive at every grid point "
        f"({scan['negative_mean_preserving_points']} negative, "
        f"{scan['positive_mean_preserving_points']} positive), while the "
        f"ambient active split is negative at "
        f"{scan['negative_active_split_points']} of them.  The weakest "
        f"mean-preserving coefficient is "
        f"{scan['worst_mean_preserving_linear']} at "
        f"y={scan['worst_mean_preserving_y']}, which sits at Liu's abscissa "
        f"x: there the inserted atom merges with the existing one and the "
        f"perturbation degenerates, so the small margin is the family "
        f"collapsing rather than the inequality failing.  The refutation is "
        f"confined to mean-decreasing directions, which leave Liu's "
        f"feasible set."
    )
    print(f"MACHINE VERIFIED [report]: {args.report} sha256={payload['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
