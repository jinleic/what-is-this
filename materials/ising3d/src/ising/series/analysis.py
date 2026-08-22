"""Exact-coefficient tools for Padé, differential, and D-finiteness analyses.

All fitting linear algebra is performed over :class:`fractions.Fraction` (via
SymPy's exact rational matrices).  Floating-point arithmetic is used only to
locate roots of the fitted exact polynomials, under an explicit ``mpmath``
precision.

The D-finiteness scan is deliberately a *finite search*, not a decision
procedure.  It tests every Euler-form homogeneous ODE

    sum_{j=0}^r Q_j(z) (z d/dz)^j F(z) = 0,

with ``deg(Q_j) <= d`` and the requested parameter budget.  Three or more
coefficients are withheld from fitting.  Failure means only that no ODE in the
searched rectangle explains the available coefficients; it is not evidence
that an ODE of larger order or degree does not exist, and is not a proof of
non-D-finiteness.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial, gcd, lcm
from typing import Iterable, Sequence

import mpmath as mp
from sympy import Matrix, Rational

__all__ = [
    "DFinitenessReport",
    "DifferentialApproximant",
    "DlogPadeApproximant",
    "ODEFit",
    "PadeApproximant",
    "derivative_series",
    "d_finiteness_scan",
    "differential_approximant_scan",
    "dlog_pade_scan",
    "even_to_squared_variable",
    "exponential_series",
    "first_order_differential_approximant",
    "guess_linear_ode",
    "logarithmic_derivative",
    "pade_approximant",
]

RationalSeries = tuple[Fraction, ...]


def _fractions(values: Sequence[int | Fraction]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def _from_sympy(value) -> Fraction:
    rational = Rational(value)
    return Fraction(int(rational.p), int(rational.q))


def _primitive_vector(values: Iterable) -> tuple[Fraction, ...]:
    fractions = tuple(_from_sympy(value) for value in values)
    denominator = 1
    for value in fractions:
        denominator = lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in fractions]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    if divisor:
        integers = [value // divisor for value in integers]
    for value in integers:
        if value:
            if value < 0:
                integers = [-entry for entry in integers]
            break
    return tuple(Fraction(value) for value in integers)


def _dot(row: Sequence[Fraction], vector: Sequence[Fraction]) -> Fraction:
    return sum((left * right for left, right in zip(row, vector)), Fraction(0))


def _poly_eval(coefficients: Sequence[Fraction], value: mp.mpf | mp.mpc):
    result = mp.mpf("0")
    for coefficient in reversed(coefficients):
        result = result * value + mp.mpf(coefficient.numerator) / coefficient.denominator
    return result


def _poly_derivative(coefficients: Sequence[Fraction]) -> RationalSeries:
    return tuple(index * value for index, value in enumerate(coefficients[1:], 1))


def _polynomial_roots(coefficients: Sequence[Fraction], precision: int):
    polynomial = list(coefficients)
    while len(polynomial) > 1 and not polynomial[-1]:
        polynomial.pop()
    degree = len(polynomial) - 1
    if degree < 1:
        return ()
    with mp.workdps(precision):
        descending = [
            mp.mpf(value.numerator) / value.denominator for value in reversed(polynomial)
        ]
        if degree == 1:
            return (-descending[1] / descending[0],)
        try:
            roots = mp.polyroots(descending, maxsteps=1000, error=False)
        except (mp.libmp.libhyper.NoConvergence, ValueError, ZeroDivisionError):
            return ()
        return tuple(+root for root in roots)


def _smallest_positive_real_root(
    coefficients: Sequence[Fraction], precision: int
) -> mp.mpf | None:
    with mp.workdps(precision):
        tolerance = mp.power(10, -max(12, precision // 3))
        candidates = []
        for root in _polynomial_roots(coefficients, precision):
            if abs(mp.im(root)) <= tolerance * max(mp.mpf(1), abs(root)):
                real_root = mp.re(root)
                if real_root > tolerance:
                    candidates.append(real_root)
        if not candidates:
            return None
        return +min(candidates)


def derivative_series(
    coefficients: Sequence[int | Fraction], times: int = 1
) -> RationalSeries:
    """Differentiate a truncated formal series ``times`` times exactly."""

    result = _fractions(coefficients)
    times = int(times)
    if times < 0:
        raise ValueError("times must be non-negative")
    for _ in range(times):
        result = tuple(index * value for index, value in enumerate(result[1:], 1))
    return result


def even_to_squared_variable(
    coefficients: Sequence[int | Fraction],
) -> RationalSeries:
    """Map an even series ``sum a[2n] v^(2n)`` to ``sum a[2n] z^n``."""

    series = _fractions(coefficients)
    if any(series[index] for index in range(1, len(series), 2)):
        raise ValueError("series contains a nonzero odd coefficient")
    return tuple(series[index] for index in range(0, len(series), 2))


def logarithmic_derivative(
    coefficients: Sequence[int | Fraction],
) -> RationalSeries:
    """Return the exact known terms of ``F'(z)/F(z)``."""

    series = _fractions(coefficients)
    if not series or not series[0]:
        raise ValueError("logarithmic derivative requires a nonzero constant term")
    derivative = derivative_series(series)
    quotient: list[Fraction] = []
    for degree, numerator in enumerate(derivative):
        convolution = sum(
            (
                quotient[index] * series[degree - index]
                for index in range(degree)
            ),
            Fraction(0),
        )
        quotient.append((numerator - convolution) / series[0])
    return tuple(quotient)


@dataclass(frozen=True)
class PadeApproximant:
    numerator_degree: int
    denominator_degree: int
    numerator: RationalSeries
    denominator: RationalSeries


def pade_approximant(
    coefficients: Sequence[int | Fraction],
    numerator_degree: int,
    denominator_degree: int,
) -> PadeApproximant | None:
    """Fit an exact ``[L/M]`` Padé approximant to a formal series."""

    series = _fractions(coefficients)
    numerator_degree = int(numerator_degree)
    denominator_degree = int(denominator_degree)
    if numerator_degree < 0 or denominator_degree < 0:
        raise ValueError("Padé degrees must be non-negative")
    required = numerator_degree + denominator_degree + 1
    if len(series) < required:
        raise ValueError(f"Padé [{numerator_degree}/{denominator_degree}] needs {required} terms")
    if denominator_degree == 0:
        return PadeApproximant(
            numerator_degree,
            0,
            tuple(series[: numerator_degree + 1]),
            (Fraction(1),),
        )

    rows = []
    right_hand_side = []
    for degree in range(
        numerator_degree + 1,
        numerator_degree + denominator_degree + 1,
    ):
        rows.append(
            [
                Rational(series[degree - shift].numerator, series[degree - shift].denominator)
                for shift in range(1, denominator_degree + 1)
            ]
        )
        value = -series[degree]
        right_hand_side.append(Rational(value.numerator, value.denominator))
    matrix = Matrix(rows)
    if matrix.rank() != denominator_degree:
        return None
    solution = matrix.inv() * Matrix(right_hand_side)
    denominator = (Fraction(1),) + tuple(_from_sympy(value) for value in solution)
    numerator = []
    for degree in range(numerator_degree + 1):
        numerator.append(
            sum(
                (
                    denominator[shift] * series[degree - shift]
                    for shift in range(min(degree, denominator_degree) + 1)
                ),
                Fraction(0),
            )
        )
    return PadeApproximant(
        numerator_degree,
        denominator_degree,
        tuple(numerator),
        denominator,
    )


@dataclass(frozen=True)
class DlogPadeApproximant:
    numerator_degree: int
    denominator_degree: int
    numerator: RationalSeries
    denominator: RationalSeries
    singularity: mp.mpf
    exponent: mp.mpf
    nearest_numerator_zero: mp.mpf | None
    defective: bool


def _nearest_root_distance(
    polynomial: Sequence[Fraction], point: mp.mpf, precision: int
) -> mp.mpf | None:
    roots = _polynomial_roots(polynomial, precision)
    if not roots:
        return None
    with mp.workdps(precision):
        return +min(abs(root - point) for root in roots)


def dlog_pade_scan(
    coefficients: Sequence[int | Fraction],
    *,
    precision: int = 80,
    near_diagonal: bool = False,
) -> tuple[DlogPadeApproximant, ...]:
    """Fit all available Dlog--Padé approximants and extract positive poles.

    The returned ``exponent`` is ``-Res(F'/F, z_c)``.  Thus it is positive for
    a divergence ``F ~ (1-z/z_c)^(-exponent)``.
    """

    precision = int(precision)
    if precision < 30:
        raise ValueError("precision must be at least 30 decimal digits")
    dlog = logarithmic_derivative(coefficients)
    results = []
    for denominator_degree in range(1, len(dlog)):
        for numerator_degree in range(0, len(dlog) - denominator_degree):
            if near_diagonal and abs(numerator_degree - denominator_degree) > 1:
                continue
            approximant = pade_approximant(
                dlog, numerator_degree, denominator_degree
            )
            if approximant is None:
                continue
            singularity = _smallest_positive_real_root(
                approximant.denominator, precision
            )
            if singularity is None:
                continue
            with mp.workdps(precision):
                derivative = _poly_derivative(approximant.denominator)
                denominator_slope = _poly_eval(derivative, singularity)
                if not denominator_slope:
                    continue
                residue = _poly_eval(approximant.numerator, singularity) / denominator_slope
                exponent = -residue
                zero_distance = _nearest_root_distance(
                    approximant.numerator, singularity, precision
                )
                cancellation_scale = mp.mpf("1e-3") * max(
                    mp.mpf(1), abs(singularity)
                )
                defective = zero_distance is not None and zero_distance <= cancellation_scale
                results.append(
                    DlogPadeApproximant(
                        numerator_degree=approximant.numerator_degree,
                        denominator_degree=approximant.denominator_degree,
                        numerator=approximant.numerator,
                        denominator=approximant.denominator,
                        singularity=+singularity,
                        exponent=+exponent,
                        nearest_numerator_zero=(+zero_distance if zero_distance is not None else None),
                        defective=bool(defective),
                    )
                )
    return tuple(results)


@dataclass(frozen=True)
class DifferentialApproximant:
    derivative_degree: int
    function_degree: int
    inhomogeneous_degree: int
    derivative_polynomial: RationalSeries
    function_polynomial: RationalSeries
    inhomogeneous_polynomial: RationalSeries
    singularity: mp.mpf
    singular_exponent: mp.mpf


def first_order_differential_approximant(
    coefficients: Sequence[int | Fraction],
    derivative_degree: int,
    function_degree: int,
    inhomogeneous_degree: int,
    *,
    precision: int = 80,
) -> DifferentialApproximant | None:
    """Fit ``Q1 F' + Q0 F = P`` and return its first positive singularity."""

    series = _fractions(coefficients)
    derivative_degree = int(derivative_degree)
    function_degree = int(function_degree)
    inhomogeneous_degree = int(inhomogeneous_degree)
    if min(derivative_degree, function_degree, inhomogeneous_degree) < 0:
        raise ValueError("differential-approximant degrees must be non-negative")
    equation_count = len(series) - 1
    unknown_count = (
        derivative_degree + function_degree + inhomogeneous_degree + 3
    )
    if unknown_count > equation_count + 1:
        raise ValueError("not enough series coefficients for this differential approximant")

    rows: list[list[Rational]] = []
    for degree in range(equation_count):
        row: list[Fraction] = []
        for shift in range(derivative_degree + 1):
            source = degree - shift + 1
            row.append(
                source * series[source] if 0 <= source < len(series) else Fraction(0)
            )
        for shift in range(function_degree + 1):
            source = degree - shift
            row.append(series[source] if source >= 0 else Fraction(0))
        for shift in range(inhomogeneous_degree + 1):
            row.append(Fraction(-1 if degree == shift else 0))
        rows.append(
            [Rational(value.numerator, value.denominator) for value in row]
        )
    nullspace = Matrix(rows).nullspace()
    if len(nullspace) != 1:
        return None
    vector = _primitive_vector(nullspace[0])
    split_one = derivative_degree + 1
    split_two = split_one + function_degree + 1
    derivative_polynomial = vector[:split_one]
    function_polynomial = vector[split_one:split_two]
    inhomogeneous_polynomial = vector[split_two:]
    if not any(derivative_polynomial):
        return None
    singularity = _smallest_positive_real_root(derivative_polynomial, precision)
    if singularity is None:
        return None
    with mp.workdps(precision):
        slope = _poly_eval(_poly_derivative(derivative_polynomial), singularity)
        if not slope:
            return None
        singular_exponent = -_poly_eval(function_polynomial, singularity) / slope
        return DifferentialApproximant(
            derivative_degree=derivative_degree,
            function_degree=function_degree,
            inhomogeneous_degree=inhomogeneous_degree,
            derivative_polynomial=derivative_polynomial,
            function_polynomial=function_polynomial,
            inhomogeneous_polynomial=inhomogeneous_polynomial,
            singularity=+singularity,
            singular_exponent=+singular_exponent,
        )


def differential_approximant_scan(
    coefficients: Sequence[int | Fraction],
    *,
    precision: int = 80,
    balanced: bool = True,
) -> tuple[DifferentialApproximant, ...]:
    """Fit first-order inhomogeneous approximants using all known coefficients."""

    series = _fractions(coefficients)
    target_unknowns = len(series)
    results = []
    for derivative_degree in range(target_unknowns - 2):
        for function_degree in range(target_unknowns - 2):
            inhomogeneous_degree = (
                target_unknowns - derivative_degree - function_degree - 3
            )
            if inhomogeneous_degree < 0:
                continue
            if balanced and (
                max(derivative_degree, function_degree, inhomogeneous_degree)
                - min(derivative_degree, function_degree, inhomogeneous_degree)
                > 2
            ):
                continue
            fit = first_order_differential_approximant(
                series,
                derivative_degree,
                function_degree,
                inhomogeneous_degree,
                precision=precision,
            )
            if fit is not None:
                results.append(fit)
    return tuple(results)


@dataclass(frozen=True)
class ODEFit:
    order: int
    degree: int
    coefficients: tuple[RationalSeries, ...]
    training_coefficients: int
    held_out_coefficients: int
    validation_residuals: RationalSeries

    @property
    def passed(self) -> bool:
        return all(residual == 0 for residual in self.validation_residuals)

    @property
    def unknown_count(self) -> int:
        return (self.order + 1) * (self.degree + 1)


def _euler_ode_row(
    series: RationalSeries, coefficient_index: int, order: int, degree: int
) -> RationalSeries:
    row = []
    for derivative_order in range(order + 1):
        for polynomial_degree in range(degree + 1):
            source = coefficient_index - polynomial_degree
            if source < 0:
                row.append(Fraction(0))
            else:
                row.append(series[source] * (source**derivative_order))
    return tuple(row)


def guess_linear_ode(
    coefficients: Sequence[int | Fraction],
    order: int,
    degree: int,
    *,
    holdout: int = 3,
) -> ODEFit | None:
    """Fit one Euler-form ODE on a prefix and validate it on held-out terms.

    A result is returned only when the training equations determine a unique
    one-dimensional nullspace.  This prevents held-out equations from being
    used to select a relation from an underdetermined family.
    """

    series = _fractions(coefficients)
    order = int(order)
    degree = int(degree)
    holdout = int(holdout)
    if order < 0 or degree < 0:
        raise ValueError("ODE order and degree must be non-negative")
    if holdout < 1 or holdout >= len(series):
        raise ValueError("holdout must leave at least one fitting coefficient")
    training_count = len(series) - holdout
    unknown_count = (order + 1) * (degree + 1)
    if unknown_count > training_count:
        raise ValueError("ODE has more unknowns than fitting coefficients")

    training_rows = [
        _euler_ode_row(series, index, order, degree)
        for index in range(training_count)
    ]
    matrix = Matrix(
        [
            [Rational(value.numerator, value.denominator) for value in row]
            for row in training_rows
        ]
    )
    nullspace = matrix.nullspace()
    if len(nullspace) != 1:
        return None
    vector = _primitive_vector(nullspace[0])
    residuals = tuple(
        _dot(_euler_ode_row(series, index, order, degree), vector)
        for index in range(training_count, len(series))
    )
    polynomials = tuple(
        tuple(
            vector[
                derivative_order * (degree + 1)
                + polynomial_degree
            ]
            for polynomial_degree in range(degree + 1)
        )
        for derivative_order in range(order + 1)
    )
    return ODEFit(
        order=order,
        degree=degree,
        coefficients=polynomials,
        training_coefficients=training_count,
        held_out_coefficients=holdout,
        validation_residuals=residuals,
    )


@dataclass(frozen=True)
class DFinitenessReport:
    known_coefficients: int
    holdout: int
    tested_pairs: tuple[tuple[int, int], ...]
    fits: tuple[ODEFit, ...]

    @property
    def found(self) -> bool:
        return any(fit.passed for fit in self.fits)

    @property
    def credible_holdout_needed(self) -> int:
        """Conservative future-coefficient target independent of this short fit."""

        largest_fit = max((fit.unknown_count for fit in self.fits), default=0)
        return max(10, largest_fit)


def d_finiteness_scan(
    coefficients: Sequence[int | Fraction],
    *,
    holdout: int = 3,
) -> DFinitenessReport:
    """Search every ``(r,d)`` with ``(r+1)(d+1) <= N-holdout`` exactly."""

    series = _fractions(coefficients)
    holdout = int(holdout)
    if holdout < 1 or len(series) <= holdout:
        raise ValueError("not enough coefficients for the requested holdout")
    budget = len(series) - holdout
    tested = []
    fits = []
    for order in range(budget):
        for degree in range(budget):
            if (order + 1) * (degree + 1) > budget:
                continue
            tested.append((order, degree))
            fit = guess_linear_ode(
                series, order, degree, holdout=holdout
            )
            if fit is not None:
                fits.append(fit)
    return DFinitenessReport(
        known_coefficients=len(series),
        holdout=holdout,
        tested_pairs=tuple(tested),
        fits=tuple(fits),
    )


def exponential_series(order: int) -> RationalSeries:
    """Exact coefficients of the holonomic control ``exp(z)`` through ``z**order``."""

    order = int(order)
    if order < 0:
        raise ValueError("order must be non-negative")
    return tuple(Fraction(1, factorial(degree)) for degree in range(order + 1))
