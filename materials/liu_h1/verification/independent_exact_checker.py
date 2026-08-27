#!/usr/bin/env python3
"""CAS-free exact checks for the algebraic core of Liu Hypothesis 1.

This checker is intentionally independent of SymPy, FLINT/Arb, NumPy, and the
canonical implementation in ``math/uc``.  It uses only integer and
``fractions.Fraction`` arithmetic.  A tiny multivariate-polynomial type checks
identities by coefficient comparison; a Leibniz determinant is adequate for
the two 4-by-4 determinants used here.

What this file checks exactly:

* the polynomial part of the entropy-kernel reduction;
* the coefficient matrix and characteristic polynomial of D_r;
* the cleared-denominator Lorentz factorization;
* exact ordered-field certificates used for the Descartes and cubic arguments;
* uniform denominator lower-bound identities at the singular corner; and
* the rational 2-point counterexample to the discarded Cauchy-kernel route.

Explicit limit: this is not a formal proof assistant.  It does not mechanize
Taylor's theorem, logarithm differentiation, the Schur product theorem,
Descartes' rule, convergence of Gram series, Fubini, or the signed-measure
extension.  Those human-audited steps and their hypotheses are recorded in
``LIU_H1/AUDIT.md``.

Run from the ``math`` directory:

    python3 -I -B LIU_H1/verification/independent_exact_checker.py
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from typing import Iterable, Sequence

# Exponent order for every polynomial in this checker.
VARIABLES = ("r", "s", "t", "lambda")
VARIABLE_COUNT = len(VARIABLES)
Exponent = tuple[int, ...]


class Polynomial:
    """Sparse exact polynomial over Q in (r, s, t, lambda)."""

    def __init__(self, terms: dict[Exponent, Fraction] | None = None) -> None:
        cleaned: dict[Exponent, Fraction] = {}
        for exponent, coefficient in (terms or {}).items():
            if len(exponent) != VARIABLE_COUNT or any(value < 0 for value in exponent):
                raise ValueError("invalid exponent")
            exact = Fraction(coefficient)
            if exact:
                cleaned[tuple(exponent)] = exact
        self.terms = cleaned

    @classmethod
    def constant(cls, value: int | Fraction) -> "Polynomial":
        exact = Fraction(value)
        return cls({(0,) * VARIABLE_COUNT: exact} if exact else {})

    @classmethod
    def variable(cls, index: int) -> "Polynomial":
        exponent = [0] * VARIABLE_COUNT
        exponent[index] = 1
        return cls({tuple(exponent): Fraction(1)})

    @staticmethod
    def coerce(value: object) -> "Polynomial":
        if isinstance(value, Polynomial):
            return value
        if isinstance(value, (int, Fraction)):
            return Polynomial.constant(value)
        raise TypeError(f"cannot coerce {type(value)!r} to Polynomial")

    def __add__(self, other: object) -> "Polynomial":
        rhs = self.coerce(other)
        terms = dict(self.terms)
        for exponent, coefficient in rhs.terms.items():
            terms[exponent] = terms.get(exponent, Fraction(0)) + coefficient
            if not terms[exponent]:
                del terms[exponent]
        return Polynomial(terms)

    def __radd__(self, other: object) -> "Polynomial":
        return self + other

    def __neg__(self) -> "Polynomial":
        return Polynomial({exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Polynomial":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "Polynomial":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "Polynomial":
        rhs = self.coerce(other)
        terms: dict[Exponent, Fraction] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in rhs.terms.items():
                exponent = tuple(
                    left + right for left, right in zip(left_exponent, right_exponent)
                )
                terms[exponent] = (
                    terms.get(exponent, Fraction(0))
                    + left_coefficient * right_coefficient
                )
        return Polynomial(terms)

    def __rmul__(self, other: object) -> "Polynomial":
        return self * other

    def __pow__(self, power: int) -> "Polynomial":
        if power < 0:
            raise ValueError("polynomial powers must be nonnegative")
        result = Polynomial.constant(1)
        base = self
        remaining = power
        while remaining:
            if remaining & 1:
                result = result * base
            base = base * base
            remaining //= 2
        return result

    def __eq__(self, other: object) -> bool:
        try:
            return self.terms == self.coerce(other).terms
        except TypeError:
            return False

    def derivative(self, variable: int) -> "Polynomial":
        terms: dict[Exponent, Fraction] = {}
        for exponent, coefficient in self.terms.items():
            degree = exponent[variable]
            if not degree:
                continue
            reduced = list(exponent)
            reduced[variable] -= 1
            terms[tuple(reduced)] = coefficient * degree
        return Polynomial(terms)

    def evaluate(self, values: Sequence[int | Fraction]) -> Fraction:
        if len(values) != VARIABLE_COUNT:
            raise ValueError("one value is required for each variable")
        exact_values = tuple(Fraction(value) for value in values)
        return sum(
            coefficient
            * product(value**power for value, power in zip(exact_values, exponent))
            for exponent, coefficient in self.terms.items()
        )


def product(values: Iterable[object]) -> object:
    result: object = 1
    for value in values:
        result = result * value
    return result


def permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: Sequence[Sequence[object]]) -> Polynomial:
    """Exact Leibniz determinant; matrices in this checker have order at most 4."""

    dimension = len(matrix)
    if any(len(row) != dimension for row in matrix):
        raise ValueError("determinant needs a square matrix")
    total = Polynomial.constant(0)
    for permutation in permutations(range(dimension)):
        term = Polynomial.constant(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            term *= Polynomial.coerce(matrix[row][column])
        total += term
    return total


def sign_variations(signs: Sequence[int]) -> int:
    nonzero = [sign for sign in signs if sign]
    if any(sign not in (-1, 1) for sign in nonzero):
        raise ValueError("signs must be -1, 0, or 1")
    return sum(left != right for left, right in zip(nonzero, nonzero[1:]))


def require(condition: bool, statement: str) -> None:
    if not condition:
        raise AssertionError(statement)
    print(f"MACHINE-VERIFIED: {statement}")


def main() -> None:
    r = Polynomial.variable(0)
    s = Polynomial.variable(1)
    t = Polynomial.variable(2)
    eigenvalue = Polynomial.variable(3)

    # Exact reduction identities.  The logarithmic factors are not represented;
    # these identities verify the polynomial/rank-one bookkeeping around them.
    u = 1 - s
    v = 1 - t
    a_s = s * u
    a_t = t * v
    A = u * v
    B = s * t
    z = A + a_s * a_t
    require(z == A * (1 + B), "z=uv+a(s)a(t)=A(1+B)")
    require(z - A * B == A, "the pointwise entropy-to-R remainder uses z-AB=A")

    projection_basis = (
        (1, 0, 0),
        (0, 1, 0),
        (0, 1, -1),
    )
    require(
        determinant(projection_basis) == -1,
        "{1,s,s(1-s)} spans every polynomial of degree at most two",
    )

    # Coefficient matrix of D_r=(1+r st)[1-(1-s)(1-t)(1+r st)].
    D = (1 + r * B) * (1 - A * (1 + r * B))
    M: list[list[object]] = [
        [0, 1, 0, 0],
        [1, -r - 1, 2 * r, 0],
        [0, 2 * r, -r * (r + 2), r**2],
        [0, 0, r**2, -(r**2)],
    ]
    s_basis = (Polynomial.constant(1), s, s**2, s**3)
    t_basis = (Polynomial.constant(1), t, t**2, t**3)
    reconstructed_D = sum(
        Polynomial.coerce(M[row][column]) * s_basis[row] * t_basis[column]
        for row in range(4)
        for column in range(4)
    )
    require(reconstructed_D == D, "the displayed 4-by-4 matrix represents D_r(s,t)")
    require(determinant(M) == -2 * r**3, "det M(r)=-2r^3")

    characteristic = determinant(
        [
            [
                (eigenvalue if row == column else 0) - Polynomial.coerce(M[row][column])
                for column in range(4)
            ]
            for row in range(4)
        ]
    )
    c1 = (r + 1) * (2 * r + 1)
    c2 = 4 * r**3 + 2 * r - 1
    c3 = -2 * r * (r**3 - r**2 + r + 1)
    c4 = -2 * r**3
    expected_characteristic = (
        eigenvalue**4
        + c1 * eigenvalue**3
        + c2 * eigenvalue**2
        + c3 * eigenvalue
        + c4
    )
    require(
        characteristic == expected_characteristic,
        "det(lambda I-M) has the five claimed coefficient polynomials",
    )

    positive_factor = r**2 - r + 1
    require(
        positive_factor == (r - Fraction(1, 2)) ** 2 + Fraction(3, 4),
        "r^2-r+1=(r-1/2)^2+3/4",
    )
    require(
        all(sign_variations((1, 1, middle, -1, -1)) == 1 for middle in (-1, 0, 1)),
        "each abstract sign pattern (+,+,+/0/-,-,-) has one variation",
    )

    # Clear the harmless positive denominators (r+1)(r+2) in the explicit
    # Lorentz coordinates.  This avoids radicals and rational functions.
    z1_s = 1 + 2 * r * s**2
    z1_t = 1 + 2 * r * t**2
    n0_s = 2 * r * s**2 - (r + 1) * s + 1
    n0_t = 2 * r * t**2 - (r + 1) * t + 1
    q_s = s**2 * (r + 2 - r * s)
    q_t = t**2 * (r + 2 - r * t)
    cleared_lorentz = (
        (r + 2) * (z1_s * z1_t - n0_s * n0_t)
        - r * (r + 1) * q_s * q_t
        - 2 * r**2 * (r + 1) * s**3 * t**3
    )
    require(
        (r + 1) * (r + 2) * D == cleared_lorentz,
        "the cleared-denominator Lorentz factorization equals D_r, including r=0",
    )

    # The cubic derivative has a manifestly negative exact form.  Together with
    # p(1)=1 this is the ordered-field certificate p(s)>=1 on [0,1].
    p = 2 - 2 * s + 2 * s**2 - s**3
    require(
        -p.derivative(1) == 3 * (s - Fraction(2, 3)) ** 2 + Fraction(2, 3),
        "-p'(s)=3(s-2/3)^2+2/3",
    )
    require(
        p.evaluate((0, 0, 0, 0)) == 2 and p.evaluate((0, 1, 0, 0)) == 1,
        "p(0)=2 and p(1)=1",
    )

    diagonal_bracket = 1 - (1 - s) ** 2 * (1 + r * s**2)
    diagonal_D = (1 + r * s**2) * diagonal_bracket
    require(
        diagonal_bracket - s * p == (1 - r) * (1 - s) ** 2 * s**2,
        "1-(1-s)^2(1+rs^2)-s p(s)=(1-r)(1-s)^2s^2",
    )
    require(
        diagonal_D - s * p
        == (1 - r) * (1 - s) ** 2 * s**2 + r * s**2 * diagonal_bracket,
        "D_r(s,s)-s p(s) has the exact displayed factor decomposition",
    )

    # This stronger exact bound supplies the endpoint estimate omitted from the
    # symbolic canonical checker.  If H=1-A(1+rB), then H>=s and H>=t on the
    # unit cube, hence D_r=(1+rB)H>=max(s,t).  Therefore AB^2/D_r tends to zero
    # uniformly in r at (0,0).
    denominator_gap = 1 - A * (1 + r * B)
    require(
        denominator_gap - s == (1 - s) * t * (1 - r * s * (1 - t)),
        "1-A(1+rB)-s=(1-s)t[1-rs(1-t)]",
    )
    require(
        denominator_gap - t == (1 - t) * s * (1 - r * t * (1 - s)),
        "1-A(1+rB)-t=(1-t)s[1-rt(1-s)]",
    )

    theta = Fraction(1, 2)
    points = (Fraction(1, 4), Fraction(3, 4))
    cauchy = tuple(
        tuple(Fraction(1, 1) / (1 + theta * left * right) for right in points)
        for left in points
    )
    cauchy_determinant = cauchy[0][0] * cauchy[1][1] - cauchy[0][1] ** 2
    require(
        cauchy_determinant == Fraction(-131072, 1657425),
        "the 2-point Gram determinant of 1/(1+st/2) is -131072/1657425",
    )

    print(
        "HUMAN-AUDITED: ordered-field signs turn the checked identities into the "
        "stated interval inequalities."
    )
    print(
        "OPEN: this checker deliberately leaves calculus, kernel-limit, and "
        "measure-theoretic theorems to the explicit audit."
    )


if __name__ == "__main__":
    main()
