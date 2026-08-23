#!/usr/bin/env python3
"""Geometric reduction for the vertical irreducibility problem.

Main unconditional result (all identities replayed exactly below):

    Descent theorem.  For every Z0 in Qbar with Z0 not in {0,1}, put
    r0=(1-Z0)/Z0^2.  If the plane pencil member S_{r0}(a,y) is absolutely
    irreducible, then H(a,Z0,b) is irreducible in Qbar[a,b]; when Z0 is
    rational it is therefore irreducible in Q(a)[b].

    Noether system.  Absolute irreducibility of S_{r0} holds iff an explicit
    797x483 integer matrix M(r0)=M0+r0*M1 has full column rank, i.e. iff the
    integer polynomial N(r)=gcd of the 483-minors does not vanish at r0.

So the last non-Schinzel obligation is reduced to: N((1-z^3)/z^6)!=0 for the
cell's z.  That final nonvanishing for every rational z in the cell family is
NOT proved here; finite scans below are labeled EVIDENCE and never promoted.

Replay from the workspace root:

    nice -n 19 python3 math/h10q/l22_fiber_geometry.py
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l22_fiber_geometry.jsonl"
REPORT = Path("/tmp/l22_fiber_geometry.md")
CERTIFICATE_PRIME = 1009
SCAN_PRIMES = (101, 103)
PACE_EVERY = 200
PACE_SECONDS = 0.005

Monomial = tuple[int, ...]
Polynomial = dict[Monomial, F]
Univariate = dict[int, F]


class Pacer:
    """Yield briefly inside the finite matrix loops."""

    def __init__(self) -> None:
        self.iterations = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.iterations += 1
        if self.iterations % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


# ----------------------------------------------------------------------------
# exact multivariate polynomial arithmetic (dict monomial -> Fraction)
# ----------------------------------------------------------------------------

def constant(dimension: int, value: int | F = 1) -> Polynomial:
    coefficient = F(value)
    return {(0,) * dimension: coefficient} if coefficient else {}


def variable(dimension: int, index: int) -> Polynomial:
    exponent = [0] * dimension
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def add(*polynomials: Polynomial) -> Polynomial:
    result: defaultdict[Monomial, F] = defaultdict(F)
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] += coefficient
    return {monomial: coefficient for monomial, coefficient in result.items()
            if coefficient}


def scale(polynomial: Polynomial, scalar: int | F) -> Polynomial:
    scalar_f = F(scalar)
    return {monomial: scalar_f * coefficient
            for monomial, coefficient in polynomial.items()
            if scalar_f * coefficient}


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: defaultdict[Monomial, F] = defaultdict(F)
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(x + y for x, y in zip(left_monomial,
                                                    right_monomial))
            result[monomial] += left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in result.items()
            if coefficient}


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    dimension = len(next(iter(polynomial))) if polynomial else 1
    result = constant(dimension)
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        remaining >>= 1
        if remaining:
            factor = multiply(factor, factor)
    return result


def derivative(polynomial: Polynomial, index: int) -> Polynomial:
    result: Polynomial = {}
    for monomial, coefficient in polynomial.items():
        exponent = monomial[index]
        if exponent:
            derived = list(monomial)
            derived[index] -= 1
            result[tuple(derived)] = coefficient * exponent
    return result


def total_degree(polynomial: Polynomial) -> int:
    return max(sum(monomial) for monomial in polynomial)


def reduce_mod_A(polynomial: Polynomial, a_index: int = 0) -> Polynomial:
    """Reduce modulo A=1+4a^2 by substituting a^2 -> -1/4."""
    result: defaultdict[Monomial, F] = defaultdict(F)
    for monomial, coefficient in polynomial.items():
        exponent = monomial[a_index]
        reduced = list(monomial)
        reduced[a_index] = exponent % 2
        result[tuple(reduced)] += coefficient * F(-1, 4) ** (exponent // 2)
    return {monomial: coefficient for monomial, coefficient in result.items()
            if coefficient}


# ----------------------------------------------------------------------------
# the pencil S_r = F0 + r*G in Q[a,y] and its (a,y,Z) companions
# ----------------------------------------------------------------------------

def integral_pencil() -> tuple[Polynomial, Polynomial, Polynomial,
                               Polynomial, Polynomial, Polynomial]:
    """Return A,V,B,E,G,F0 with S_r=F0+r*G=E+(r-a^2)*G in Q[a,y]."""
    one = constant(2)
    a = variable(2, 0)
    y = variable(2, 1)
    a_minus_one = add(a, scale(one, -1))
    A = add(one, scale(power(a, 2), 4))
    V = add(multiply(A, power(y, 2)), scale(one, 4))
    B = scale(multiply(A, power(a_minus_one, 2)), 2)
    E = add(
        scale(
            multiply(
                multiply(
                    multiply(power(a, 8), A),
                    power(a_minus_one, 4),
                ),
                power(V, 2),
            ),
            64,
        ),
        scale(multiply(power(a, 4), power(add(V, scale(B, -1)), 4)), -1),
    )
    G = scale(
        multiply(
            multiply(multiply(power(A, 3), power(a_minus_one, 4)), y),
            power(V, 2),
        ),
        4,
    )
    F0 = add(E, scale(multiply(power(a, 2), G), -1))
    return A, V, B, E, G, F0


def norm_identity_record() -> dict[str, object]:
    """Verify 2H=2X^2+A*L*Y^2 and the cleared parametrization B^4*R=A*T."""
    one = constant(3)
    a = variable(3, 0)
    Z = variable(3, 1)
    b = variable(3, 2)
    A = add(one, scale(power(a, 2), 4))
    D = add(one, scale(Z, -1),
            scale(multiply(power(a, 2), power(Z, 2)), -1))
    L = add(scale(one, 2),
            scale(multiply(multiply(A, power(add(a, scale(one, -1)), 2)), b),
                  -1))
    Ng = add(scale(multiply(power(a, 4), power(b, 2)), 16),
             scale(multiply(A, power(add(b, scale(one, -1)), 4)), -1))
    X = multiply(multiply(power(a, 4), power(Z, 2)), Ng)
    Y = scale(multiply(multiply(A, D), power(b, 2)), 2)
    H = add(
        multiply(multiply(power(a, 8), power(Z, 4)), power(Ng, 2)),
        scale(multiply(multiply(power(A, 3), power(D, 2)), power(b, 4)), 4),
        scale(
            multiply(
                multiply(
                    multiply(power(A, 4), power(add(a, scale(one, -1)), 2)),
                    power(D, 2),
                ),
                power(b, 5),
            ),
            -2,
        ),
    )
    assert scale(H, 2) == add(scale(power(X, 2), 2),
                              multiply(multiply(A, L), power(Y, 2)))

    # (a,y,Z) companions and the denominator-cleared R.
    one3 = constant(3)
    a3 = variable(3, 0)
    y3 = variable(3, 1)
    Z3 = variable(3, 2)
    am1 = add(a3, scale(one3, -1))
    A3 = add(one3, scale(power(a3, 2), 4))
    D3 = add(one3, scale(Z3, -1),
             scale(multiply(power(a3, 2), power(Z3, 2)), -1))
    V3 = add(multiply(A3, power(y3, 2)), scale(one3, 4))
    B3 = scale(multiply(A3, power(am1, 2)), 2)
    E3 = add(
        scale(
            multiply(
                multiply(multiply(power(a3, 8), A3), power(am1, 4)),
                power(V3, 2),
            ),
            64,
        ),
        scale(multiply(power(a3, 4),
                       power(add(V3, scale(B3, -1)), 4)), -1),
    )
    G3 = scale(
        multiply(
            multiply(multiply(power(A3, 3), power(am1, 4)), y3),
            power(V3, 2),
        ),
        4,
    )

    cleared_X = multiply(
        multiply(power(a3, 4), power(Z3, 2)),
        add(
            scale(multiply(multiply(power(a3, 4), power(V3, 2)),
                           power(B3, 2)), 16),
            scale(multiply(A3, power(add(V3, scale(B3, -1)), 4)), -1),
        ),
    )
    cleared_tY = multiply(
        multiply(multiply(power(A3, 2), D3), y3),
        multiply(power(V3, 2), power(B3, 2)),
    )
    cleared_R = add(cleared_X, cleared_tY)
    T = add(multiply(power(Z3, 2), E3), multiply(D3, G3))
    assert cleared_R == multiply(A3, T)

    return {
        "type": "norm-birational-reduction",
        "label": "PROVED",
        "norm_field": "K'=K(t), K=Qbar(Z)(a,b), t^2=-A*L/2",
        "norm_element": "R=X+tY, X=a^4*Z^2*N_g, Y=2*A*D*b^2",
        "identity": "H=X^2+(A*L/2)*Y^2=Norm_{K'/K}(R)",
        "rational_coordinates": "t=A*y/2, V=A*y^2+4, B=2*A*(a-1)^2, b=V/B",
        "involution_identity": "t^2+A*L(a,V/B)/2=0 holds identically",
        "cleared_identity": "B^4*(X+tY)|_{b=V/B} = A*T, T=Z^2*E+D*G",
        "parameter": "r=(1-Z)/Z^2, D=Z^2*(r-a^2), T=Z^2*S_r",
        "pencil": (
            "S_r=E+(r-a^2)*G, E=64*a^8*A*(a-1)^4*V^2-a^4*(V-B)^4, "
            "G=4*A^3*(a-1)^4*y*V^2"
        ),
    }


def descent_theorem_record() -> dict[str, object]:
    """Verify the product identity and every fact used by the descent proof.

    Identity (*):  sum_j H_j(a,Z)*V^j*B^(8-j) = A^2*(Z^4*E^2-D^2*G^2)
    in Q[a,y,Z], where H_j are the b-coefficients of H.
    """
    one = constant(3)
    a = variable(3, 0)
    y = variable(3, 1)
    Z = variable(3, 2)
    am1 = add(a, scale(one, -1))
    A = add(one, scale(power(a, 2), 4))
    D = add(one, scale(Z, -1),
            scale(multiply(power(a, 2), power(Z, 2)), -1))
    V = add(multiply(A, power(y, 2)), scale(one, 4))
    B = scale(multiply(A, power(am1, 2)), 2)
    E = add(
        scale(
            multiply(
                multiply(multiply(power(a, 8), A), power(am1, 4)),
                power(V, 2),
            ),
            64,
        ),
        scale(multiply(power(a, 4), power(add(V, scale(B, -1)), 4)), -1),
    )
    G = scale(
        multiply(multiply(multiply(power(A, 3), power(am1, 4)), y),
                 power(V, 2)),
        4,
    )

    # H in variables (a,Z,b), then b-coefficients as (a,y,Z) polynomials.
    oneH = constant(3)
    aH = variable(3, 0)
    ZH = variable(3, 1)
    bH = variable(3, 2)
    AH = add(oneH, scale(power(aH, 2), 4))
    DH = add(oneH, scale(ZH, -1),
             scale(multiply(power(aH, 2), power(ZH, 2)), -1))
    NgH = add(scale(multiply(power(aH, 4), power(bH, 2)), 16),
              scale(multiply(AH, power(add(bH, scale(oneH, -1)), 4)), -1))
    HH = add(
        multiply(multiply(power(aH, 8), power(ZH, 4)), power(NgH, 2)),
        scale(multiply(multiply(power(AH, 3), power(DH, 2)), power(bH, 4)),
              4),
        scale(
            multiply(
                multiply(
                    multiply(power(AH, 4),
                             power(add(aH, scale(oneH, -1)), 2)),
                    power(DH, 2),
                ),
                power(bH, 5),
            ),
            -2,
        ),
    )
    H_coefficients: dict[int, Polynomial] = {j: {} for j in range(9)}
    for (a_exp, Z_exp, b_exp), coefficient in HH.items():
        H_coefficients[b_exp][(a_exp, 0, Z_exp)] = coefficient

    # (fact 1) deg_b H = 8 with H_8 = a^8*Z^4*A^2.
    assert H_coefficients[8] == multiply(multiply(power(a, 8), power(Z, 4)),
                                         power(A, 2))
    # (fact 2) H(a,Z,0) = a^8*Z^4*A^2, so content(H) | a^8*Z^4*A^2 in Qbar[a].
    assert H_coefficients[0] == multiply(multiply(power(a, 8), power(Z, 4)),
                                         power(A, 2))
    # (fact 3) H(0,Z,b) = 2*(1-Z)^2*b^4*(2-b): the factor a is excluded
    # from the content whenever Z0 != 1.
    H_at_a0 = {(Z_exp, b_exp): coefficient
               for (a_exp, Z_exp, b_exp), coefficient in HH.items()
               if a_exp == 0}
    one2 = constant(2)
    Z2 = variable(2, 0)
    b2 = variable(2, 1)
    expected_a0 = scale(
        multiply(
            multiply(power(add(one2, scale(Z2, -1)), 2), power(b2, 4)),
            add(scale(one2, 2), scale(b2, -1)),
        ),
        2,
    )
    assert H_at_a0 == expected_a0
    # (fact 4) H mod A = Z^4*b^4/256: the factors of A are excluded from the
    # content whenever Z0 != 0.
    assert reduce_mod_A(HH) == {(0, 4, 4): F(1, 256)}

    # (*) the product identity.
    left = constant(3, 0)
    for b_exp in range(9):
        left = add(left,
                   multiply(H_coefficients[b_exp],
                            multiply(power(V, b_exp), power(B, 8 - b_exp))))
    right = multiply(
        power(A, 2),
        add(multiply(power(Z, 4), power(E, 2)),
            scale(multiply(power(D, 2), power(G, 2)), -1)),
    )
    assert left == right

    # (fact 5) evenness bookkeeping: E is even in y, G is odd in y.
    assert all(monomial[1] % 2 == 0 for monomial in E)
    assert all(monomial[1] % 2 == 1 for monomial in G)
    # (fact 6) top form of S_r is 4096*a^18*y^4 for every r
    # (deg E=22 > max(deg a^2*G, deg G)+... : deg G=19, so deg S_r=22).
    assert E[(18, 4, 0)] == 4096
    assert total_degree(E) == 22 and total_degree(G) == 19

    return {
        "type": "fixed-fibre-descent-theorem",
        "label": "PROVED",
        "product_identity": (
            "sum_j H_j(a,Z)*V^j*B^(8-j) = A^2*(Z^4*E^2-D^2*G^2)"
            " = A^2*T(a,y,Z)*T(a,-y,Z)"
        ),
        "identity_terms": len(left),
        "H_b_degree": 8,
        "H_leading": "a^8*Z^4*A^2",
        "H_at_b0": "a^8*Z^4*A^2",
        "H_at_a0": "2*(1-Z)^2*b^4*(2-b)",
        "H_mod_A": "Z^4*b^4/256",
        "statement": (
            "for Z0 in Qbar-{0,1} and r0=(1-Z0)/Z0^2: S_{r0} absolutely "
            "irreducible => H(a,Z0,b) irreducible in Qbar[a,b] => (Z0 "
            "rational) H(a,Z0,b) irreducible in Q(a)[b]"
        ),
        "proof_ingredients": [
            "content(H_Z0)=1 in Qbar[a] from facts 2-4",
            "pullback h_i(a,y)=B^(d_i)*H_i(a,V/B) is even in y with "
            "y-degree exactly 2*d_i (leading c_{d_i}(a)*A^(d_i))",
            "h_1*h_2=A^2*Z0^4*S_{r0}(a,y)*S_{r0}(a,-y) by (*)",
            "S_{r0}(a,y) and S_{r0}(a,-y) are non-associate irreducibles "
            "(odd part (r0-a^2)*G nonzero, even part nonzero)",
            "unique factorization in Qbar[a,y] forces one h_i into Qbar[a], "
            "contradicting y-degree 2*d_i>=2",
            "Gauss: primitive+irreducible in Qbar[a][b] => irreducible in "
            "Qbar(a)[b] => irreducible in Q(a)[b]",
        ],
        "cells": (
            "every cell z is a nonzero rational with v_w(z)>=1, so "
            "Z=z^3 avoids {0,1} and the theorem applies with "
            "r0=(1-z^3)/z^6"
        ),
    }


# ----------------------------------------------------------------------------
# generic geometric integrality of the pencil
# ----------------------------------------------------------------------------

def substitute_a(polynomial: Polynomial, value: int | F) -> Univariate:
    value_f = F(value)
    result: defaultdict[int, F] = defaultdict(F)
    for (a_exponent, y_exponent), coefficient in polynomial.items():
        result[y_exponent] += coefficient * value_f ** a_exponent
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def substitute_diagonal(polynomial: Polynomial) -> Univariate:
    result: defaultdict[int, F] = defaultdict(F)
    for (a_exponent, y_exponent), coefficient in polynomial.items():
        result[a_exponent + y_exponent] += coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def univariate_list(polynomial: Univariate) -> list[F]:
    if not polynomial:
        return [F(0)]
    return [polynomial.get(exponent, F(0))
            for exponent in range(max(polynomial) + 1)]


def trim(polynomial: list[F]) -> list[F]:
    result = list(polynomial)
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def divide_univariate(dividend: list[F],
                      divisor: list[F]) -> tuple[list[F], list[F]]:
    remainder = trim(dividend)
    divisor_t = trim(divisor)
    if len(divisor_t) == 1 and not divisor_t[0]:
        raise ZeroDivisionError("zero polynomial")
    quotient = [F(0)] * max(1, len(remainder) - len(divisor_t) + 1)
    while len(remainder) >= len(divisor_t) and any(remainder):
        offset = len(remainder) - len(divisor_t)
        coefficient = remainder[-1] / divisor_t[-1]
        quotient[offset] = coefficient
        for index, value in enumerate(divisor_t):
            remainder[offset + index] -= coefficient * value
        remainder = trim(remainder)
    return trim(quotient), remainder


def gcd_univariate(left: list[F], right: list[F]) -> list[F]:
    first = trim(left)
    second = trim(right)
    while not (len(second) == 1 and not second[0]):
        _, remainder = divide_univariate(first, second)
        first, second = second, remainder
    leading = first[-1]
    return [coefficient / leading for coefficient in first]


def restriction_record(E: Polynomial, G: Polynomial,
                       F0: Polynomial) -> dict[str, object]:
    """Coprimality of (F0,G) and two restrictions proving non-compositeness."""
    # gcd(F0,G)=1: G=4*A^3*(a-1)^4*y*V^2 has known irreducible factors, and
    # each is excluded from F0 by an exact reduction.
    assert reduce_mod_A(F0) == {(0, 0): F(-16)}
    at_a1: defaultdict[int, F] = defaultdict(F)
    for (a_exponent, y_exponent), coefficient in F0.items():
        at_a1[y_exponent] += coefficient
    one2 = constant(2)
    y2 = variable(2, 1)
    expected_a1 = scale(power(add(scale(power(y2, 2), 5),
                                  scale(one2, 4)), 4), -1)
    assert ({(0, exponent): coefficient
             for exponent, coefficient in at_a1.items() if coefficient}
            == expected_a1)
    at_y0 = {a_exponent: coefficient
             for (a_exponent, y_exponent), coefficient in F0.items()
             if y_exponent == 0}
    assert max(at_y0) == 20 and at_y0[20] == -4096 and min(at_y0) == 4
    # mod V: F0 = 64a^8*A*(a-1)^4*V^2 - a^4*(V-B)^4 - a^2*G is congruent to
    # -a^4*B^4, a nonzero y-free polynomial, and every nonzero multiple of
    # V=A*y^2+4 has y-degree >= 2.  So V does not divide F0.  (formal)

    vertical_f = substitute_a(F0, 2)
    vertical_g = substitute_a(G, 2)
    vertical_gcd = gcd_univariate(univariate_list(vertical_f),
                                  univariate_list(vertical_g))
    assert vertical_gcd == [F(1)]
    assert (min(vertical_f), max(vertical_f)) == (0, 8)
    assert (min(vertical_g), max(vertical_g)) == (1, 5)
    assert vertical_f[0] == -8503552
    assert vertical_f[8] == -1336336

    diagonal_f = substitute_diagonal(F0)
    diagonal_g = substitute_diagonal(G)
    diagonal_gcd = gcd_univariate(univariate_list(diagonal_f),
                                  univariate_list(diagonal_g))
    assert diagonal_gcd == [F(0), F(1)]
    assert (min(diagonal_f), max(diagonal_f)) == (3, 22)
    assert (min(diagonal_g), max(diagonal_g)) == (1, 19)
    assert diagonal_f[3] == -64 and diagonal_f[22] == 4096
    assert diagonal_g[1] == 64 and diagonal_g[19] == 4096

    return {
        "type": "generic-geometric-integrality",
        "label": "PROVED",
        "pencil": "S_r=F0+r*G",
        "coprimality": {
            "F0_mod_A": "-16",
            "F0_at_a1": "-(5*y^2+4)^4",
            "F0_at_y0": "degree 20, leading -4096, a-valuation 4",
            "F0_mod_V": "-a^4*B^4 (formal; V has y-degree 2)",
        },
        "rational_function": "psi=-F0/G in Qbar(a,y)",
        "restriction_degrees": {
            "a=2": {"degree_after_gcd": 8, "gcd": "1"},
            "y=a": {"degree_after_gcd": 21, "gcd": "a"},
        },
        "noncomposite_reason": "an outer factor's degree divides gcd(8,21)=1",
        "bertini_krull_conclusion": (
            "S_r is irreducible in overline(Qbar(r))[a,y] for transcendental r"
        ),
        "base_change": (
            "r=(1-Z)/Z^2 identifies overline(Qbar(r))=overline(Qbar(Z)); "
            "the generic-Z fibre polynomial Z^2*S_r=T is geometrically "
            "integral"
        ),
    }


# ----------------------------------------------------------------------------
# Ruppert/Noether system for the fixed spectrum
# ----------------------------------------------------------------------------

def homogenize(polynomial: Polynomial, degree: int) -> Polynomial:
    result: Polynomial = {}
    for (a_exponent, y_exponent), coefficient in polynomial.items():
        homogenizing_exponent = degree - a_exponent - y_exponent
        if homogenizing_exponent < 0:
            raise ValueError("homogenization degree too small")
        result[(a_exponent, y_exponent, homogenizing_exponent)] = coefficient
    return result


def ruppert_domain_basis(degree: int) -> list[tuple[Polynomial, Polynomial]]:
    """Basis of pairs (U,W), homogeneous degree d-1, with c dividing aU+yW."""
    basis: list[tuple[Polynomial, Polynomial]] = []
    for a_exponent in range(degree):
        for y_exponent in range(degree - a_exponent):
            c_exponent = degree - 1 - a_exponent - y_exponent
            if c_exponent >= 1:
                monomial = {(a_exponent, y_exponent, c_exponent): F(1)}
                basis.append((monomial, {}))
                basis.append(({}, monomial))
    for a_exponent in range(degree - 1):
        y_exponent = degree - 2 - a_exponent
        basis.append((
            {(a_exponent, y_exponent + 1, 0): F(1)},
            {(a_exponent + 1, y_exponent, 0): F(-1)},
        ))
    return basis


def ruppert_column(f: Polynomial, g: Polynomial, h: Polynomial) -> Polynomial:
    """One column of the map (U,W) -> [f(U_y - W_a) - U f_y + W f_a]/c."""
    numerator = add(
        multiply(f, add(derivative(g, 1), scale(derivative(h, 0), -1))),
        scale(multiply(g, derivative(f, 1)), -1),
        multiply(h, derivative(f, 0)),
    )
    if any(monomial[2] < 1 for monomial in numerator):
        raise AssertionError("Ruppert numerator is not divisible by c")
    return {(a_exponent, y_exponent, c_exponent - 1): coefficient
            for (a_exponent, y_exponent, c_exponent), coefficient
            in numerator.items()}


def build_ruppert_matrix(F0: Polynomial, G: Polynomial,
                          pacer: Pacer) -> tuple[
                              dict[Monomial, dict[int, tuple[int, int]]],
                              dict[str, object],
                          ]:
    degree = 22
    f0_h = homogenize(F0, degree)
    g_h = homogenize(G, degree)
    basis = ruppert_domain_basis(degree)
    assert len(basis) == degree * degree - 1 == 483

    rows: defaultdict[Monomial, dict[int, tuple[int, int]]] = defaultdict(dict)
    digest = hashlib.sha256()
    nonzero_f0 = 0
    nonzero_g = 0
    for column_index, (basis_g, basis_h) in enumerate(basis):
        column_f0 = ruppert_column(f0_h, basis_g, basis_h)
        column_g = ruppert_column(g_h, basis_g, basis_h)
        nonzero_f0 += len(column_f0)
        nonzero_g += len(column_g)
        for part, column in enumerate((column_f0, column_g)):
            for monomial, coefficient in sorted(column.items()):
                assert coefficient.denominator == 1
                digest.update(
                    f"{column_index}|{part}|{monomial}|"
                    f"{coefficient.numerator}\n".encode()
                )
        for monomial in set(column_f0) | set(column_g):
            rows[monomial][column_index] = (
                int(column_f0.get(monomial, 0)),
                int(column_g.get(monomial, 0)),
            )
        pacer.tick()

    assert all(sum(monomial) == 41 for monomial in rows)
    nominal_rows = math.comb(43, 2)
    assert nominal_rows == 903
    stats: dict[str, object] = {
        "type": "ruppert-noether-system",
        "label": "PROVED",
        "homogeneous_degree": degree,
        "matrix": "M(r)=M0+r*M1 with integer entries",
        "nominal_shape": [nominal_rows, len(basis)],
        "nonzero_row_shape": [len(rows), len(basis)],
        "M0_nonzero_entries": nonzero_f0,
        "M1_nonzero_entries": nonzero_g,
        "canonical_sha256": digest.hexdigest(),
        "criterion": (
            "for fixed r0: homogeneous S_{r0} (degree 22, top form "
            "4096*a^18*y^4 for all r0) is irreducible over Qbar iff "
            "rank M(r0)=483; equivalently iff N(r0)!=0 where N in Z[r] is "
            "the gcd of the 483-minors (Ruppert; Buse-Cheze Thm 8)"
        ),
        "noether_polynomial": (
            "N(r) integer polynomial, deg<=483; N!=0 by the generic "
            "theorem; N(0)=0; N not computed here"
        ),
    }
    return dict(rows), stats


def reduce_rows_mod(rows: dict[Monomial, dict[int, tuple[int, int]]],
                    prime: int) -> list[list[tuple[int, int, int]]]:
    """Per-row lists of (column, M0 entry, M1 entry) reduced mod prime."""
    reduced: list[list[tuple[int, int, int]]] = []
    for monomial in sorted(rows):
        entries = []
        for column, (constant_part, linear_part) in rows[monomial].items():
            c0 = constant_part % prime
            c1 = linear_part % prime
            if c0 or c1:
                entries.append((column, c0, c1))
        if entries:
            reduced.append(entries)
    return reduced


def rank_at(rows_mod: list[list[tuple[int, int, int]]], parameter: int,
            prime: int, pacer: Pacer) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for entries in rows_mod:
        row: dict[int, int] = {}
        for column, c0, c1 in entries:
            value = (c0 + parameter * c1) % prime
            if value:
                row[column] = value
        while row:
            pivot_column = min(row)
            if pivot_column not in pivots:
                inverse = pow(row[pivot_column], -1, prime)
                row = {column: value * inverse % prime
                       for column, value in row.items()
                       if value * inverse % prime}
                pivots[pivot_column] = row
                break
            pivot = pivots[pivot_column]
            multiplier = row[pivot_column]
            for column, value in pivot.items():
                reduced = (row.get(column, 0) - multiplier * value) % prime
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
        pacer.tick()
        if len(pivots) == 483:
            break
    return len(pivots)


def modular_fraction(value: F, prime: int) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    if not denominator:
        raise ZeroDivisionError("parameter denominator vanishes modulo prime")
    return numerator * pow(denominator, -1, prime) % prime


def certificate_record(rows_mod: list[list[tuple[int, int, int]]],
                       parameter: F, description: str,
                       pacer: Pacer) -> dict[str, object]:
    parameter_mod = modular_fraction(parameter, CERTIFICATE_PRIME)
    rank = rank_at(rows_mod, parameter_mod, CERTIFICATE_PRIME, pacer)
    assert rank == 483
    return {
        "type": "fixed-fibre-ruppert-certificate",
        "label": "PROVED",
        "description": description,
        "r": str(parameter),
        "prime": CERTIFICATE_PRIME,
        "r_mod_prime": parameter_mod,
        "rank": rank,
        "columns": 483,
        "char0_conclusion": (
            "a maximal integer minor is nonzero mod 1009, hence nonzero in "
            "Q; N(r0)!=0; S_{r0} is absolutely irreducible; by the descent "
            "theorem the vertical fibre is irreducible in Q(a)[b]"
        ),
    }


def degeneration_record(rows_mod: list[list[tuple[int, int, int]]],
                        F0: Polynomial, pacer: Pacer) -> dict[str, object]:
    """The r=0 member: exactly the a^2 content degenerates."""
    a_valuation = min(a_exponent for (a_exponent, _) in F0)
    assert a_valuation == 2
    rank0 = rank_at(rows_mod, 0, CERTIFICATE_PRIME, pacer)
    assert rank0 == 480
    return {
        "type": "r0-degeneration",
        "label": "PROVED",
        "r": "0",
        "structure": "F0=a^2*U with U(0,y)!=0, deg U=20",
        "rank_mod_1009": rank0,
        "kernel_dimension_char0": 3,
        "argument": (
            "rank_Q>=rank_F1009=480 gives dim ker<=3; the Ruppert dimension "
            "formula for f=a^2*prod(u_i^{e_i}) gives dim ker=s-1+"
            "binom(3+sum d_i(e_i-1),2)>=3 with equality iff s=1 and U "
            "squarefree; hence U is absolutely irreducible and the only "
            "degeneration at r=0 (i.e. Z=1) is the visible a^2 content"
        ),
        "interpretation": (
            "r=a^2 is where D vanishes: D=Z^2*(r-a^2); at r=0 that locus "
            "degenerates onto a=0 doubly; Z=1 is excluded from every cell"
        ),
    }


def scan_record(rows: dict[Monomial, dict[int, tuple[int, int]]],
                prime: int, pacer: Pacer) -> dict[str, object]:
    """Full residue scan of rank drops mod prime.  EVIDENCE only."""
    rows_mod = reduce_rows_mod(rows, prime)
    drops = []
    for residue in range(prime):
        rank = rank_at(rows_mod, residue, prime, pacer)
        if rank < 483:
            drops.append([residue, rank])
    return {
        "type": "spectrum-residue-scan",
        "label": "EVIDENCE",
        "prime": prime,
        "rank_drops": drops,
        "reading": (
            "every rational root u/v of N with gcd(u,v)=1 and prime not "
            "dividing v reduces to a rank-drop residue; drops == {0} means "
            "any nonzero rational root must have numerator divisible by "
            "the prime"
        ),
        "not_a_theorem": (
            "a finite scan cannot prove N has no nonzero rational root; "
            "irrational spectrum points may also avoid every scanned prime"
        ),
    }


def distinctions_record() -> dict[str, object]:
    return {
        "type": "field-distinction",
        "label": "PROVED",
        "statements": {
            "generic_arithmetic": (
                "irreducibility in Q(a,Z)[b] (already known for H/P)"
            ),
            "generic_geometric": (
                "irreducibility in overline(Q(Z))[a,b], equivalently after "
                "every constant-field extension"
            ),
            "fixed_bivariate_arithmetic": (
                "irreducibility of H(a,Z0,b) in Q(a)[b]"
            ),
            "fixed_bivariate_absolute": (
                "irreducibility of the primitive H_Z0 in Qbar[a,b]"
            ),
            "cover": (
                "fixed bivariate absolute irreducibility equals geometric "
                "connectedness of the a-line cover; the descent theorem "
                "transports it from the pencil member S_{r0}"
            ),
            "univariate": (
                "H(a0,Z0,b) in Q[b] may be arithmetically irreducible, but "
                "a degree-8 univariate always splits over Qbar"
            ),
        },
        "nonimplication": (
            "generic arithmetic irreducibility implies neither generic "
            "geometric irreducibility nor any fixed-fibre assertion; "
            "b^2-a^2-r is generically geometrically irreducible yet its "
            "r=0 member factors"
        ),
    }


def open_record() -> dict[str, object]:
    return {
        "type": "uniform-cell-verdict",
        "label": "OPEN",
        "proved": [
            "descent: N(r_z)!=0 => vertical fibre irreducible in Q(a)[b], "
            "r_z=(1-z^3)/z^6",
            "the generic member of the r-pencil is geometrically integral, "
            "so N!=0",
            "the r=0 degeneration is exactly the a^2 content",
        ],
        "not_proved": (
            "N(r_z)!=0 for every nonzero rational z with v_w(z)>=1"
        ),
        "missing_lemma": (
            "compute the integer polynomial N(r) (gcd of the 483-minors of "
            "M0+r*M1) and prove N(r)=c*r^m; failing that, prove that "
            "z^(6*deg N)*N((1-z^3)/z^6) has no rational root in the cell "
            "family"
        ),
        "why_finiteness_is_insufficient": (
            "Bertini-Noether finiteness of the spectrum leaves finitely "
            "many bad r; nothing yet excludes a rational bad r of the "
            "special shape (1-z^3)/z^6"
        ),
        "finite_checks_are_uniform_evidence": False,
        "refusals_are_evidence": False,
    }


def write_report(records: list[dict[str, object]]) -> None:
    matrix = next(row for row in records
                  if row["type"] == "ruppert-noether-system")
    descent = next(row for row in records
                   if row["type"] == "fixed-fibre-descent-theorem")
    certificates = [row for row in records
                    if row["type"] == "fixed-fibre-ruppert-certificate"]
    scans = [row for row in records
             if row["type"] == "spectrum-residue-scan"]
    degeneration = next(row for row in records
                        if row["type"] == "r0-degeneration")
    open_row = next(row for row in records
                    if row["type"] == "uniform-cell-verdict")

    certificate_lines = [
        f"- `{row['description']}`: `r={row['r']}`, rank `{row['rank']}/483` "
        f"modulo `{row['prime']}` -- a characteristic-zero **PROVED** "
        "certificate for that fixed fibre."
        for row in certificates
    ]
    scan_lines = [
        f"- `p={row['prime']}`: rank drops exactly at residues "
        f"`{[d[0] for d in row['rank_drops']]}` "
        f"(ranks `{[d[1] for d in row['rank_drops']]}`).  **EVIDENCE**."
        for row in scans
    ]
    lines = [
        "# L22: fibre geometry and the vertical spectrum",
        "",
        "## Status",
        "",
        "- **PROVED** (descent theorem): for `Z0 not in {0,1}`, absolute "
        "irreducibility of the pencil member `S_{r0}`, `r0=(1-Z0)/Z0^2`, "
        "implies `H(a,Z0,b)` irreducible in `Q(a)[b]`.",
        "- **PROVED**: the generic member of the pencil is geometrically "
        "integral, so the Noether polynomial `N` is not identically zero.",
        "- **PROVED**: `S_{r0}` absolutely irreducible `<=> N(r0)!=0`, with "
        "`N` the gcd of the `483`-minors of an explicit integer matrix "
        "`M0+r*M1`.",
        "- **PROVED**: at `r=0` (i.e. `Z=1`, excluded from all cells) the "
        "only degeneration is the visible `a^2` content.",
        "- **EVIDENCE**: full residue scans modulo 101 and 103 find rank "
        "drops only at `r=0`.",
        "- **OPEN**: `N((1-z^3)/z^6)!=0` for every cell `z`.",
        "",
        "No finite scan is promoted to a theorem; refusals are never "
        "evidence.",
        "",
        "## 1. Four different irreducibility statements",
        "",
        "Put `k=Qbar`.  The proved statement `P` irreducible in `Q(a,Z)[b]` "
        "is *arithmetic generic* irreducibility.  *Geometric generic* means "
        "irreducible in `overline(Q(Z))[a,b]`.  For fixed `Z_0`, *arithmetic "
        "vertical* means irreducible in `Q(a)[b]`, and *absolute vertical* "
        "means the primitive part of `H(a,Z_0,b)` is irreducible in "
        "`k[a,b]`; the latter equals geometric connectedness of the cover "
        "of the `a`-line.  After also fixing `a=a_0`, a degree-8 univariate "
        "can be irreducible in `Q[b]` while splitting completely over `k`.  "
        "These levels must not be conflated: `b^2-a^2-r` is geometrically "
        "irreducible for generic `r` yet factors at `r=0`, so a generic "
        "theorem never settles all vertical fibres by itself.",
        "",
        "## 2. Norm identity and the plane-pencil model",
        "",
        "On the constructed branch, with",
        "",
        "```text",
        "A=1+4a^2,  L=2-A(a-1)^2 b,  D=1-Z-a^2 Z^2,",
        "X=a^4 Z^2 N_g,  Y=2 A D b^2,",
        "```",
        "",
        "exact expansion gives the norm form `H=X^2+(A L/2)Y^2`: over "
        "`K'=K(t)`, `t^2=-A L/2`, `K=k(Z)(a,b)`, one has "
        "`H=Norm_{K'/K}(X+tY)`.  The cover is rationalized by",
        "",
        "```text",
        "t=A y/2,  V=A y^2+4,  B=2A(a-1)^2,  b=V/B,",
        "```",
        "",
        "under which `t^2+A L(a,V/B)/2=0` holds *identically*, and with "
        "`r=(1-Z)/Z^2` (so `D=Z^2(r-a^2)`):",
        "",
        "```text",
        "B^4 (X+tY)|_{b=V/B} = A T,   T=Z^2 E + D G = Z^2 S_r,",
        "S_r = E+(r-a^2)G,",
        "E = 64a^8 A(a-1)^4 V^2 - a^4 (V-B)^4,   G = 4A^3(a-1)^4 y V^2.",
        "```",
        "",
        "`E` is even and `G` is odd in `y`; `deg S_r=22` with top form "
        "`4096 a^18 y^4` for every `r`.",
        "",
        "## 3. The descent theorem (PROVED)",
        "",
        "**Theorem.** Let `Z_0 in k`, `Z_0 not in {0,1}`, "
        "`r_0=(1-Z_0)/Z_0^2`.  If `S_{r_0}` is irreducible in `k[a,y]`, "
        "then `H(a,Z_0,b)` is irreducible in `k[a,b]`, hence in `k(a)[b]`, "
        "hence in `Q(a)[b]` when `Z_0` is rational.",
        "",
        "*Proof.*  Write `H=sum_j H_j(a,Z) b^j`, `deg_b H=8`.  The script "
        "verifies the exact identity",
        "",
        "```text",
        "(*)   sum_j H_j(a,Z) V^j B^(8-j) = A^2 (Z^4 E^2 - D^2 G^2)",
        "                                 = A^2 T(a,y,Z) T(a,-y,Z).",
        "```",
        "",
        "Content: `H_0=H(a,Z_0,0)=a^8 Z_0^4 A^2`, so any content prime of "
        "`H_{Z_0}` in `k[a]` divides `a` or a factor of `A`; but "
        "`H(0,Z_0,b)=2(1-Z_0)^2 b^4 (2-b)!=0` (uses `Z_0!=1`) and "
        "`H mod A=Z_0^4 b^4/256!=0` (uses `Z_0!=0`).  So `H_{Z_0}` is "
        "primitive over `k[a]` and any nontrivial factorization "
        "`H_{Z_0}=H_1 H_2` has `deg_b H_i=d_i>=1`, `d_1+d_2=8`.",
        "",
        "Pull back: `h_i(a,y):=B^{d_i} H_i(a,V/B)` is a polynomial, *even* "
        "in `y` (both `V` and `B` are), with exact `y`-degree `2d_i>=2` "
        "(leading term `c_{d_i}(a) A^{d_i} y^{2d_i}`).  Identity `(*)` "
        "gives",
        "",
        "```text",
        "h_1 h_2 = A^2 Z_0^4 S_{r_0}(a,y) S_{r_0}(a,-y).",
        "```",
        "",
        "`S_{r_0}(a,y)` and `S_{r_0}(a,-y)` are irreducible (hypothesis) "
        "and non-associate: the odd part `(r_0-a^2)G` and the even part are "
        "both nonzero, so `S(a,-y)=cS(a,y)` would force `c=1` and `c=-1` "
        "simultaneously.  In the UFD `k[a,y]`, applying `y->-y` to the "
        "factorization of the even `h_i` shows `S(a,y)` and `S(a,-y)` occur "
        "in `h_i` with *equal* multiplicities.  Since each appears exactly "
        "once in the product, exactly one `h_i` contains both and the other "
        "contains neither; the latter divides `A^2 Z_0^4`, so it lies in "
        "`k[a]` -- contradicting `y`-degree `2d_i>=2`.  So `H_{Z_0}` is "
        "irreducible in `k[a,b]`; Gauss (primitivity) upgrades to "
        "`k(a)[b]`, and `Q(a)[b]`-irreducibility follows because a "
        "factorization over the subfield `Q(a)` would persist over `k(a)`. "
        " QED.",
        "",
        f"The verified identity `(*)` has `{descent['identity_terms']}` "
        "terms on each side.  Every cell has a nonzero rational `z` with "
        "`v_w(z)>=1`, so `Z=z^3` avoids `{0,1}` and the theorem applies "
        "with `r_z=(1-z^3)/z^6`.",
        "",
        "## 4. Generic geometric integrality (PROVED)",
        "",
        "`gcd(F_0,G)=1` where `F_0=E-a^2G`: the irreducible factors of "
        "`G=4A^3(a-1)^4yV^2` are excluded from `F_0` by four exact "
        "reductions (`F_0 mod A=-16`, `F_0(1,y)=-(5y^2+4)^4`, "
        "`F_0(a,0)` nonzero, `F_0 mod V=-a^4B^4`).  Let `psi=-F_0/G`.  If "
        "`psi=u o h` with `deg u=n>1`, then `n` divides the degree of the "
        "restriction of `psi` to any line where it is nonconstant.  Exact "
        "computation gives degree `8` on `a=2` (gcd `1`) and degree `21` on "
        "`y=a` (gcd exactly `a`), and `gcd(8,21)=1`.  So `psi` is "
        "non-composite, and Bertini--Krull for the coprime pencil "
        "`mu F_0+lambda G` (char 0) makes the spectrum finite and the "
        "generic member irreducible over `overline(Qbar(r))`.  The "
        "substitution `r=(1-Z)/Z^2` is an isomorphism of `Qbar(r)` onto a "
        "subfield of `Qbar(Z)` over which `Qbar(Z)` is algebraic, so "
        "`overline(Qbar(r))=overline(Qbar(Z))` and the generic-`Z` fibre "
        "is geometrically integral.  (Independently, each PROVED "
        "certificate below exhibits an absolutely irreducible member, "
        "which already forces non-compositeness.)",
        "",
        "## 5. Noether forms for the fixed spectrum (PROVED)",
        "",
        "Homogenize `S_r` to `f_r` of degree `d=22` in `(a,y,c)`; since the "
        "top form is `4096a^18y^4` for every `r`, no fibre loses degree and "
        "`c` never divides `f_r`.  Ruppert's characteristic-zero criterion "
        "(Buse--Cheze, Theorem 8): `f` of degree `d` is irreducible over "
        "`k` iff the linear map",
        "",
        "```text",
        "(U,W) -> [f(U_y-W_a)-U f_y+W f_a]/c",
        "```",
        "",
        "on pairs of degree-`(d-1)` forms with `c | aU+yW` has zero "
        "kernel; the domain has dimension `d^2-1=483`.  In monomial bases "
        "the matrix is affine in `r`:",
        "",
        f"```text\nM(r)=M0+r*M1, nominal shape {matrix['nominal_shape']}, "
        f"nonzero rows {matrix['nonzero_row_shape'][0]},\n"
        f"nnz(M0)={matrix['M0_nonzero_entries']}, "
        f"nnz(M1)={matrix['M1_nonzero_entries']}, integer entries,\n"
        f"sha256={matrix['canonical_sha256']}\n```",
        "",
        "The `483`-minors are integer polynomials in `r` of degree at most "
        "`483`; their gcd `N(r) in Z[r]` is the Noether polynomial: "
        "**`S_{r_0}` is absolutely reducible iff `N(r_0)=0`**.  `N!=0` by "
        "section 4.  Combining with section 3:",
        "",
        "```text",
        "For a cell z (so Z=z^3 not in {0,1}):",
        "N((1-z^3)/z^6)!=0  =>  H(a,z^3,b) irreducible in Q(a)[b].",
        "```",
        "",
        "### The `r=0` member",
        "",
        f"`F_0` has `a`-valuation exactly 2, and rank `M(0)` mod 1009 is "
        f"`{degeneration['rank_mod_1009']}`, so the characteristic-zero "
        "kernel has dimension exactly 3.  The Ruppert dimension formula "
        "for `f=a^2 prod u_i^{e_i}` reads "
        "`dim ker = s-1+binom(3+sum d_i(e_i-1),2) >= 3`, with equality iff "
        "`s=1` and `U=F_0/a^2` is squarefree.  Hence `U` is absolutely "
        "irreducible: at `r=0` -- that is `Z=1`, which no cell contains -- "
        "the *only* degeneration is the `a^2` content (the collapsed "
        "`D=Z^2(r-a^2)` locus).",
        "",
        "## 6. Exact finite corroboration",
        "",
        *certificate_lines,
        *scan_lines,
        "",
        "A full modular column rank at fixed `r_0` proves `N(r_0)!=0` in "
        "characteristic zero (a nonzero minor mod `p` lifts).  The residue "
        "scans say: any nonzero rational root `u/v` of `N` (lowest terms) "
        "must have `101*103 | u` or `101*103 | v`-type obstructions -- "
        "strong structure evidence for the spectrum lemma below, but not a "
        "proof.",
        "",
        "## 7. No-go and the exact next lemma",
        "",
        "The uniform verdict stays **OPEN**.  What is missing is precisely:",
        "",
        "```text",
        "Spectrum lemma (target):  N(r)=c*r^m, c in Z-{0}, m>=1.",
        "```",
        "",
        "Since every cell has `r_z=(1-z^3)/z^6!=0`, the lemma would close "
        "every cell unconditionally through the descent theorem.  If the "
        "computed `N` has more factors, the fallback is the finitely-many "
        "root check of `z^(6 deg N) N((1-z^3)/z^6)` over the cell family.  "
        "Until then, generic integrality, Bertini--Noether finiteness, and "
        "the scans do not remove the vertical premise: a thin/finite locus "
        "can still contain rational points of the special shape.",
        "",
        "## Reference",
        "",
        "L. Buse and G. Cheze, *On the total order of reducibility of a "
        "pencil of algebraic plane curves*, arXiv:0812.4706 -- Bertini--"
        "Krull equivalence (introduction) and the homogeneous Ruppert-map "
        "kernel criterion with the multiplicity-aware dimension formula "
        "(Theorem 8).",
        "",
        f"Machine verdict: **{open_row['label']}**.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    started = time.monotonic()
    pacer = Pacer()
    records: list[dict[str, object]] = [distinctions_record()]
    records.append(norm_identity_record())
    records.append(descent_theorem_record())

    A, V, B, E, G, F0 = integral_pencil()
    assert total_degree(F0) == 22
    assert total_degree(G) == 19
    records.append(restriction_record(E, G, F0))

    rows, matrix_record = build_ruppert_matrix(F0, G, pacer)
    records.append(matrix_record)

    certificate_rows = reduce_rows_mod(rows, CERTIFICATE_PRIME)
    records.append(degeneration_record(certificate_rows, F0, pacer))
    records.append(certificate_record(certificate_rows, F(1),
                                      "pencil parameter r=1", pacer))
    Z0 = F(27)
    r0 = (1 - Z0) / Z0**2
    assert r0 == F(-26, 729)
    records.append(certificate_record(certificate_rows, r0,
                                      "rational-cube fibre Z=27=3^3", pacer))

    for prime in SCAN_PRIMES:
        record = scan_record(rows, prime, pacer)
        assert record["rank_drops"] == [[0, 480]]
        records.append(record)

    records.append(open_record())

    elapsed = time.monotonic() - started
    records.append({
        "type": "summary",
        "label": "OPEN",
        "descent_theorem": "PROVED",
        "generic_geometric": "PROVED",
        "fixed_spectrum_noether_system": "PROVED",
        "r0_degeneration": "PROVED",
        "residue_scans": "EVIDENCE",
        "all_cells": "OPEN",
        "records": len(records) + 1,
        "pacer_iterations": pacer.iterations,
        "pacer_sleeps": pacer.sleeps,
        "wall_seconds": round(elapsed, 6),
    })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "\n".join(json.dumps(record, sort_keys=True, separators=(",", ":"),
                             default=str)
                  for record in records) + "\n",
        encoding="utf-8",
    )
    write_report(records)
    print(
        "L22 fibre geometry: descent theorem PROVED; generic geometric "
        "PROVED; spectrum reduced to N(r)!=0; scans EVIDENCE; uniform "
        f"verdict OPEN; wall={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
