#!/usr/bin/env python3
"""Exact square-class families and obstructions on the canonical branch.

Replay from the workspace repository root:

    nice -n 19 python3 math/h10q/l23_squareclass.py

The script is stdlib-only.  Symbolic identities and local obstructions are
labelled PROVED.  Bounded searches are discovery evidence only; an empty scan
is never used as a negative conclusion.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h10q  # noqa: E402 -- repository exact-arithmetic authority


OUT = HERE / "data" / "l23_squareclass.jsonl"
REPORT = Path("/tmp/l23_squareclass.md")
PACE_EVERY = 250
PACE_SECONDS = 0.003
SQUARECLASS_MENU = (-15, -7, -3, -1, 1, 3, 5, 7, 15)
SEARCH_HEIGHT = 9
IRRED_PRIME_LIMIT = 60


class Pacer:
    def __init__(self) -> None:
        self.steps = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.steps += 1
        if self.steps % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def is_square_fraction(value: F) -> bool:
    value = F(value)
    return (
        value >= 0
        and math.isqrt(value.numerator) ** 2 == value.numerator
        and math.isqrt(value.denominator) ** 2 == value.denominator
    )


def poly_add(left: list[F], right: list[F]) -> list[F]:
    answer = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] += coefficient
    for index, coefficient in enumerate(right):
        answer[index] += coefficient
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_scale(poly: list[F], scalar: F | int) -> list[F]:
    return [F(scalar) * coefficient for coefficient in poly]


def poly_mul(left: list[F], right: list[F]) -> list[F]:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_pow(poly: list[F], exponent: int) -> list[F]:
    assert exponent >= 0
    answer = [F(1)]
    base = poly
    while exponent:
        if exponent & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return answer


def poly_eval(poly: list[F], value: F | int) -> F:
    answer = F(0)
    value = F(value)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


# Exact sparse Laurent polynomials.  Negative exponents are used only for b^-1
# in reciprocal trace identities.
VARIABLES = ("a", "A", "Z", "D", "s", "b", "rho", "X", "Y", "L")
NVAR = len(VARIABLES)
Monomial = tuple[int, ...]
Sparse = dict[Monomial, F]
ZERO_MONOMIAL = (0,) * NVAR


def sp_clean(poly: Sparse) -> Sparse:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def sp_const(value: F | int) -> Sparse:
    value = F(value)
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def sp_var(name: str, exponent: int = 1) -> Sparse:
    powers = [0] * NVAR
    powers[VARIABLES.index(name)] = exponent
    return {tuple(powers): F(1)}


def sp_add(*polys: Sparse) -> Sparse:
    answer: Sparse = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            answer[monomial] = answer.get(monomial, F(0)) + coefficient
    return sp_clean(answer)


def sp_scale(poly: Sparse, scalar: F | int) -> Sparse:
    scalar = F(scalar)
    return sp_clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def sp_neg(poly: Sparse) -> Sparse:
    return sp_scale(poly, -1)


def sp_sub(left: Sparse, right: Sparse) -> Sparse:
    return sp_add(left, sp_neg(right))


def sp_mul(*polys: Sparse) -> Sparse:
    answer = sp_const(1)
    for poly in polys:
        product: Sparse = {}
        for left_monomial, left_coefficient in answer.items():
            for right_monomial, right_coefficient in poly.items():
                monomial = tuple(x + y for x, y in zip(left_monomial, right_monomial))
                product[monomial] = product.get(monomial, F(0)) + left_coefficient * right_coefficient
        answer = sp_clean(product)
    return answer


def sp_pow(poly: Sparse, exponent: int) -> Sparse:
    assert exponent >= 0
    answer = sp_const(1)
    base = poly
    while exponent:
        if exponent & 1:
            answer = sp_mul(answer, base)
        base = sp_mul(base, base)
        exponent >>= 1
    return answer


def symbolic_branch_identities() -> dict[str, Any]:
    a, A, Z, D, s, b, rho, X, Y, L = (sp_var(name) for name in VARIABLES)
    one = sp_const(1)
    b_minus_one = sp_add(b, sp_const(-1))
    Ng = sp_add(
        sp_scale(sp_mul(sp_pow(a, 4), sp_pow(b, 2)), 16),
        sp_neg(sp_mul(A, sp_pow(b_minus_one, 4))),
    )
    defined_L = sp_sub(one, sp_scale(sp_mul(A, sp_pow(s, 2), b), 2))
    defined_X = sp_mul(sp_pow(a, 4), sp_pow(Z, 2), Ng)
    defined_Y = sp_scale(sp_mul(A, D, sp_pow(b, 2)), 2)
    H = sp_add(sp_pow(defined_X, 2), sp_mul(A, defined_L, sp_pow(defined_Y, 2)))
    AP = sp_add(
        sp_scale(sp_mul(sp_pow(D, 2), sp_pow(A, 3), sp_pow(b, 4)), 16),
        sp_scale(sp_mul(sp_pow(a, 8), sp_pow(Z, 4), sp_pow(Ng, 2)), 4),
        sp_scale(sp_mul(sp_pow(A, 4), sp_pow(s, 2), sp_pow(D, 2), sp_pow(b, 5)), -32),
    )
    assert not sp_sub(AP, sp_scale(H, 4))

    formal_H = sp_add(sp_pow(X, 2), sp_mul(A, L, sp_pow(Y, 2)))
    quartic_product = sp_mul(
        sp_sub(X, sp_mul(A, rho, Y)),
        sp_add(X, sp_mul(A, rho, Y)),
    )
    factor_remainder = sp_mul(A, sp_add(L, sp_mul(A, sp_pow(rho, 2))), sp_pow(Y, 2))
    assert not sp_sub(sp_sub(formal_H, quartic_product), factor_remainder)

    # Reciprocal quartic trace identity on L=-A*rho^2.
    u_minus_two = sp_add(b, sp_var("b", -1), sp_const(-2))
    reciprocal_core = sp_add(
        sp_scale(sp_mul(sp_pow(a, 4), sp_pow(Z, 2), sp_pow(a, 4)), 16),
        sp_neg(sp_mul(sp_pow(a, 4), sp_pow(Z, 2), A, sp_pow(u_minus_two, 2))),
    )
    for sign in (-1, 1):
        factor = sp_add(defined_X, sp_scale(sp_mul(A, rho, defined_Y), sign))
        trace = sp_add(
            reciprocal_core,
            sp_scale(sp_mul(sp_pow(A, 2), rho, D), 2 * sign),
        )
        assert not sp_sub(factor, sp_mul(sp_pow(b, 2), trace))

    return {
        "type": "symbolic-branch-identities",
        "label": "PROVED",
        "canonical_kernel_terms": len(H),
        "cleared_P_terms": len(AP),
        "identities": [
            "A*P=4*H",
            "H=X^2+A*L*Y^2",
            "H-(X-A*rho*Y)(X+A*rho*Y)=A(L+A*rho^2)Y^2",
            "F_eps/b^2=a^4 Z^2(16a^4-A(u-2)^2)+eps*2A^2 rho D",
        ],
        "quartic_trace_discriminant": "8*A*a^4*Z^2*(8*a^8*Z^2+eps*A^2*rho*D)",
    }


# A second sparse ring for the exact quartic discriminant determinant.
DMonomial = tuple[int, int, int]  # powers of q, C, E
DPoly = dict[DMonomial, int]
DZERO = (0, 0, 0)


def dp_clean(poly: DPoly) -> DPoly:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def dp_const(value: int) -> DPoly:
    return {} if value == 0 else {DZERO: value}


def dp_var(index: int) -> DPoly:
    powers = [0, 0, 0]
    powers[index] = 1
    return {tuple(powers): 1}


def dp_add(*polys: DPoly) -> DPoly:
    answer: DPoly = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            answer[monomial] = answer.get(monomial, 0) + coefficient
    return dp_clean(answer)


def dp_scale(poly: DPoly, scalar: int) -> DPoly:
    return dp_clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def dp_mul(*polys: DPoly) -> DPoly:
    answer = dp_const(1)
    for poly in polys:
        product: DPoly = {}
        for left_monomial, left_coefficient in answer.items():
            for right_monomial, right_coefficient in poly.items():
                monomial = tuple(x + y for x, y in zip(left_monomial, right_monomial))
                product[monomial] = product.get(monomial, 0) + left_coefficient * right_coefficient
        answer = dp_clean(product)
    return answer


def dp_pow(poly: DPoly, exponent: int) -> DPoly:
    answer = dp_const(1)
    for _ in range(exponent):
        answer = dp_mul(answer, poly)
    return answer


def dp_det(matrix: list[list[DPoly]]) -> DPoly:
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    answer: DPoly = {}
    for permutation in itertools.permutations(range(size)):
        product = dp_const(1)
        inversions = 0
        for row, column in enumerate(permutation):
            entry = matrix[row][column]
            if not entry:
                product = {}
                break
            product = dp_mul(product, entry)
            inversions += sum(permutation[prior] > column for prior in range(row))
        if product:
            answer = dp_add(answer, dp_scale(product, -1 if inversions % 2 else 1))
        PACER.tick()
    return answer


def sylvester_matrix(left: list[DPoly], right: list[DPoly]) -> list[list[DPoly]]:
    m, n = len(left) - 1, len(right) - 1
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    size = m + n
    zero = dp_const(0)
    matrix: list[list[DPoly]] = []
    for shift in range(n):
        matrix.append([zero] * shift + left_descending + [zero] * (n - 1 - shift))
    for shift in range(m):
        matrix.append([zero] * shift + right_descending + [zero] * (m - 1 - shift))
    assert len(matrix) == size and all(len(row) == size for row in matrix)
    return matrix


def symbolic_discriminant() -> dict[str, Any]:
    q, C, E = dp_var(0), dp_var(1), dp_var(2)
    one = dp_const(1)
    coefficients = [
        C,
        dp_scale(C, -2),
        dp_mul(C, dp_add(one, dp_scale(q, -2))),
        dp_scale(dp_mul(q, C), 2),
        dp_add(dp_mul(dp_pow(q, 2), C), E),
    ]
    derivative = [dp_scale(coefficients[index], index) for index in range(1, 5)]
    resultant = dp_det(sylvester_matrix(coefficients, derivative))
    A_square = dp_pow(dp_add(one, dp_scale(q, 4)), 2)
    discriminant = dp_scale(
        dp_mul(dp_pow(C, 3), dp_pow(E, 2), dp_add(dp_mul(A_square, C), dp_scale(E, 16))),
        16,
    )
    expected_resultant = dp_mul(coefficients[-1], discriminant)
    assert resultant == expected_resultant
    return {
        "type": "quartic-discriminant",
        "label": "PROVED",
        "polynomial": "C(1-Z-a^2 Z^2)^2+E Z^4",
        "factorization": "Disc_Z=16*C^3*E^2*(A^2*C+16*E)",
        "sylvester_dimension": 7,
        "resultant_terms": len(resultant),
        "genus_for_Phi_pair": 1,
    }


def canonical_data(a: F, Z: F, b: F) -> tuple[F, F, F, F, F, list[F], F]:
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    D = 1 - Z - a * a * Z * Z
    delta = -4 * a**4 / A
    coefficients = h10q._l10_P(a, Z, D, A, delta, s)
    value = poly_eval(coefficients, b)
    return A, s, D, delta, value, coefficients, 2 * b


def dyadic_replay() -> dict[str, Any]:
    samples = (
        (F(-3), F(1), F(-1, 37)),
        (F(-1), F(27), F(3, 5)),
        (F(1), F(1, 2), F(1)),
        (F(3), F(2), F(19, 37)),
        (F(5), F(1, 4), F(-1, 101)),
        (F(7), F(4), F(11, 197)),
    )
    rows = []
    for a, Z, b in samples:
        assert h10q.vp(a, 2) == 0 and h10q.vp(b, 2) == 0
        A, s, D, _, P_value, _, d = canonical_data(a, Z, b)
        m = h10q.vp(Z, 2)
        expected = 4 + 4 * min(m, 0)
        assert h10q.vp(P_value, 2) == expected
        assert h10q.hilbert(A, d, 2) == -1
        rows.append(
            {
                "a": frac_text(a), "Z": frac_text(Z), "b": frac_text(b),
                "v2_Z": m, "v2_P": expected,
                "hilbert_A_2b_at_2": -1,
            }
        )
    return {
        "type": "dyadic-obstructions",
        "label": "PROVED",
        "formula": "v2(P)=4+4*min(v2(Z),0)",
        "consequences": [
            "P=2b*y^2 has no Phi-admissible point",
            "P=-2b*y^2 has no Phi-admissible point",
            "P=A*Norm_Q(sqrt(2b))/Q forces the target Hilbert symbol -1 at 2",
        ],
        "sample_replays": rows,
    }


def reciprocal_symbolic() -> tuple[dict[str, Any], list[F]]:
    Z = sp_var("Z")
    b = sp_var("b")
    one = sp_const(1)
    D = sp_add(one, sp_neg(Z), sp_neg(sp_pow(Z, 2)))
    Ng = sp_add(
        sp_scale(sp_pow(b, 2), 16),
        sp_scale(sp_pow(sp_add(b, sp_const(-1)), 4), -5),
    )
    P = sp_add(
        sp_scale(sp_mul(sp_pow(D, 2), sp_pow(b, 4)), 400),
        sp_scale(sp_mul(sp_pow(Z, 4), sp_pow(Ng, 2)), F(4, 5)),
    )
    u_minus_two = sp_add(b, sp_var("b", -1), sp_const(-2))
    T = sp_add(
        sp_scale(sp_pow(D, 2), 400),
        sp_scale(
            sp_mul(sp_pow(Z, 4), sp_pow(sp_add(sp_const(16), sp_scale(sp_pow(u_minus_two, 2), -5)), 2)),
            F(4, 5),
        ),
    )
    assert not sp_sub(P, sp_mul(sp_pow(b, 4), T))
    T_at_two = sp_add(
        sp_scale(sp_pow(D, 2), 400),
        sp_scale(sp_pow(Z, 4), F(1024, 5)),
    )
    T_at_minus_two = sp_add(
        sp_scale(sp_pow(D, 2), 400),
        sp_scale(sp_pow(Z, 4), F(16384, 5)),
    )
    norm_factor_64 = sp_add(
        sp_scale(sp_pow(D, 2), 125),
        sp_scale(sp_pow(Z, 4), 64),
    )
    norm_factor_1024 = sp_add(
        sp_scale(sp_pow(D, 2), 125),
        sp_scale(sp_pow(Z, 4), 1024),
    )
    assert not sp_sub(sp_scale(T_at_two, F(4, 5)), sp_scale(norm_factor_64, F(64, 25)))
    assert not sp_sub(
        sp_scale(T_at_minus_two, F(4, 5)),
        sp_scale(norm_factor_1024, F(64, 25)),
    )

    P_at_b_one = T_at_two
    comparison_square = sp_scale(sp_pow(Z, 2), F(84, 5))
    factors = sp_mul(
        sp_add(sp_scale(Z, 2), sp_const(5)),
        sp_add(sp_scale(Z, 4), sp_const(5)),
        sp_add(sp_scale(Z, 7), sp_const(-5)),
        sp_add(sp_scale(Z, 9), sp_const(-5)),
    )
    assert not sp_sub(
        sp_sub(P_at_b_one, sp_pow(comparison_square, 2)),
        sp_scale(factors, F(16, 25)),
    )
    comparison_square_two = sp_scale(sp_pow(Z, 2), F(324, 5))
    factors_two = sp_mul(
        sp_add(sp_scale(Z, 8), sp_const(5)),
        sp_add(sp_scale(Z, 13), sp_const(-5)),
        sp_add(sp_scale(sp_pow(Z, 2), 54), sp_scale(Z, -25), sp_const(25)),
    )
    assert not sp_sub(
        sp_sub(P_at_b_one, sp_pow(comparison_square_two, 2)),
        sp_scale(factors_two, F(-16, 25)),
    )

    point_specs = [
        (F(-5, 2), F(84, 5)), (F(-5, 4), F(84, 5)),
        (F(5, 7), F(84, 5)), (F(5, 9), F(84, 5)),
        (F(-5, 8), F(324, 5)), (F(5, 13), F(324, 5)),
    ]
    roots = [root for root, _ in point_specs]
    point_rows = []
    for root, multiplier in point_specs:
        _, _, _, _, value, _, d = canonical_data(F(1), root, F(1))
        y = multiplier * root * root
        assert value == y * y and not h10q.ramified(value, d)
        point_rows.append({"Z": frac_text(root), "y": frac_text(y), "Z_is_rational_cube": False})

    parity_rows = []
    for p_parity, q_parity in ((0, 1), (1, 0), (1, 1)):
        d_mod_8 = (q_parity * q_parity - p_parity * q_parity - p_parity * p_parity) % 2
        assert d_mod_8 == 1
        row = {"p_mod_2": p_parity, "q_mod_2": q_parity}
        for coefficient in (64, 1024):
            residue = (125 + coefficient * p_parity**4) % 8
            assert residue == 5
            row[f"125d2_plus_{coefficient}p4_mod_8"] = residue
        parity_rows.append(row)

    return (
        {
            "type": "reciprocal-capell",
            "label": "PROVED",
            "scope": "a=1, every rational nonzero Z",
            "trace_identity": "P=b^4*T(u), T=400D^2+(4/5)Z^4(16-5(u-2)^2)^2",
            "square_norm_obstructions": [
                "125D^2+64Z^4",
                "125D^2+1024Z^4",
            ],
            "norm_identities": [
                "(4/5)T(2)=(64/25)(125D^2+64Z^4)",
                "(4/5)T(-2)=(64/25)(125D^2+1024Z^4)",
            ],
            "conclusion": "2*beta is not a square in Q[beta]/(P)",
            "parity_replay": parity_rows,
            "off_cell_square_points": point_rows,
            "point_factor_identities": [
                "P(1,Z,1)-(84Z^2/5)^2=(16/25)(2Z+5)(4Z+5)(7Z-5)(9Z-5)",
                "P(1,Z,1)-(324Z^2/5)^2=-(16/25)(8Z+5)(13Z-5)(54Z^2-25Z+25)",
            ],
        },
        roots,
    )


def structured_family_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    # L is a square and Phi-compatible, but this does not force splitting.
    a, Z = F(-3), F(1)
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    u = (a - 1) / 4
    b = -1 / (A * u * u)
    L = 1 - 2 * A * s * s * b
    _, _, _, _, P_value, _, d = canonical_data(a, Z, b)
    assert L == 9 and h10q.vp(b, 2) == 0
    assert h10q.hilbert(P_value, d, 37) == -1
    rows.append(
        {
            "type": "L-square-counterexample", "label": "PROVED",
            "a": "-3", "Z": "1", "b": "-1/37", "L": "9",
            "bad_place": 37, "hilbert": -1,
            "conclusion": "L square is Phi-compatible but is not a splitting identity",
        }
    )

    # Parent factor family: exact Phi point and a forced bad prime.
    a, rho, Z = F(3), F(1), F(1)
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    b = (1 + A * rho * rho) / (2 * A * s * s)
    L = 1 - 2 * A * s * s * b
    _, _, D, _, P_value, _, d = canonical_data(a, Z, b)
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    X = a**4 * Z * Z * Ng
    Y = 2 * A * D * b * b
    assert L == -A * rho * rho
    assert (A * P_value / 4) == (X - A * rho * Y) * (X + A * rho * Y)
    assert h10q.vp(b, 2) == 0 and h10q.hilbert(P_value, d, 19) == -1
    rows.append(
        {
            "type": "L-minus-A-square-obstruction", "label": "PROVED",
            "a": "3", "rho": "1", "Z": "1", "b": "19/37",
            "L": "-37", "forced_prime": 19, "hilbert": -1,
            "factorization": "H=(X-A*rho*Y)(X+A*rho*Y)",
        }
    )

    # Exact aligned exception: factorization alone neither forces failure nor success.
    a, rho, Z = F(-1), F(1), F(27)
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    b = (1 + A * rho * rho) / (2 * A * s * s)
    _, _, D, _, P_value, _, d = canonical_data(a, Z, b)
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    X = a**4 * Z * Z * Ng
    Y = 2 * A * D * b * b
    assert 1 - 2 * A * s * s * b == -A * rho * rho
    assert A * P_value / 4 == (X - A * rho * Y) * (X + A * rho * Y)
    ramification = h10q.ramified(P_value, d)
    assert not ramification
    rows.append(
        {
            "type": "L-minus-A-square-aligned-specialization",
            "label": "PROVED_SPECIALIZATION",
            "a": "-1", "rho": "1", "Z": "27", "b": "3/5",
            "ramified_places": ramification,
            "scope": "one actual cube Z=3^3; not an identity in Z",
        }
    )
    return rows


def rational_search(known_roots: list[F]) -> dict[str, Any]:
    z_values = sorted(
        {
            F(numerator, denominator)
            for denominator in range(1, SEARCH_HEIGHT + 1)
            for numerator in range(-SEARCH_HEIGHT, SEARCH_HEIGHT + 1)
            if numerator and math.gcd(abs(numerator), denominator) == 1
        }
    )
    r_values = (F(1, 3), F(1), F(3))
    a_values = (F(-3), F(-1), F(1), F(3), F(5))
    hits: list[dict[str, str]] = []
    checked = 0
    cube_checked = 0
    cube_hits: list[dict[str, str]] = []
    for a in a_values:
        for c in SQUARECLASS_MENU:
            for r in r_values:
                b = F(c) * r * r
                if h10q.vp(b, 2) != 0:
                    continue
                for Z in z_values:
                    _, _, _, _, P_value, _, _ = canonical_data(a, Z, b)
                    checked += 1
                    if P_value and is_square_fraction(P_value):
                        hits.append({"condition": "P square", "a": frac_text(a), "Z": frac_text(Z), "b": frac_text(b)})
                    if P_value and is_square_fraction(P_value / (2 * b)):
                        hits.append({"condition": "P=2b square", "a": frac_text(a), "Z": frac_text(Z), "b": frac_text(b)})
                    if P_value and is_square_fraction(-P_value / (2 * b)):
                        hits.append({"condition": "P=-2b square", "a": frac_text(a), "Z": frac_text(Z), "b": frac_text(b)})
                    PACER.tick()
                for z in (F(-3), F(-2), F(-1), F(1, 2), F(1), F(2), F(3)):
                    Z = z**3
                    _, _, _, _, P_value, _, _ = canonical_data(a, Z, b)
                    cube_checked += 1
                    if P_value and is_square_fraction(P_value):
                        cube_hits.append({"a": frac_text(a), "z": frac_text(z), "b": frac_text(b)})
                    PACER.tick()
    expected = {(F(1), root, F(1)) for root in known_roots if root in z_values}
    observed = {(F(row["a"]), F(row["Z"]), F(row["b"])) for row in hits if row["condition"] == "P square"}
    assert expected <= observed
    assert not [row for row in hits if row["condition"] in ("P=2b square", "P=-2b square")]
    return {
        "type": "bounded-rational-search", "label": "EVIDENCE",
        "box": {
            "a": [frac_text(value) for value in a_values],
            "squareclasses_c": list(SQUARECLASS_MENU),
            "r": [frac_text(value) for value in r_values],
            "Z_height": SEARCH_HEIGHT,
        },
        "inputs_checked": checked,
        "hits": hits,
        "cube_inputs_checked": cube_checked,
        "cube_P_square_hits": cube_hits,
        "interpretation": "identity discovery/counterexample search only; no-hit subsets are not negative evidence",
    }


# Small generic finite-field engine for degree-16 Capell-composition evidence.
def int_primitive(poly: list[F]) -> list[int]:
    denominator = 1
    for coefficient in poly:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in poly]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, abs(coefficient))
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return integers


def fp_trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def fp_sub(left: list[int], right: list[int], prime: int) -> list[int]:
    answer = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] = (answer[index] + coefficient) % prime
    for index, coefficient in enumerate(right):
        answer[index] = (answer[index] - coefficient) % prime
    return fp_trim(answer)


def fp_divmod(dividend: list[int], divisor: list[int], prime: int) -> tuple[list[int], list[int]]:
    remainder = fp_trim([coefficient % prime for coefficient in dividend])
    divisor = fp_trim([coefficient % prime for coefficient in divisor])
    assert divisor != [0]
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    inverse = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (remainder[index + shift] - coefficient * value) % prime
        fp_trim(remainder)
    return fp_trim(quotient), fp_trim(remainder)


def fp_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left, right = fp_trim(left[:]), fp_trim(right[:])
    while right != [0]:
        _, remainder = fp_divmod(left, right, prime)
        left, right = right, remainder
    inverse = pow(left[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in left]


def fp_mul_mod(left: list[int], right: list[int], modulus: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            product[i + j] = (product[i + j] + x * y) % prime
    _, remainder = fp_divmod(fp_trim(product), modulus, prime)
    return remainder


def fp_pow_mod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    answer = [1]
    while exponent:
        if exponent & 1:
            answer = fp_mul_mod(answer, base, modulus, prime)
        base = fp_mul_mod(base, base, modulus, prime)
        exponent >>= 1
    return answer


def fp_irreducible(poly: list[int], prime: int) -> bool:
    reduced = fp_trim([coefficient % prime for coefficient in poly])
    degree = len(reduced) - 1
    if degree <= 0 or reduced[-1] == 0:
        return False
    inverse = pow(reduced[-1], -1, prime)
    modulus = [(coefficient * inverse) % prime for coefficient in reduced]
    x = [0, 1]
    frobenius = x
    halfway = None
    for step in range(1, degree + 1):
        frobenius = fp_pow_mod(frobenius, prime, modulus, prime)
        if step == degree // 2:
            halfway = frobenius
    assert halfway is not None
    if fp_sub(frobenius, x, prime) != [0]:
        return False
    return len(fp_gcd(modulus, fp_sub(halfway, x, prime), prime)) == 1

def fp_factor_degrees(poly: list[int], prime: int) -> list[int] | None:
    """Distinct-degree factorization; None means bad/non-squarefree reduction."""
    reduced = fp_trim([coefficient % prime for coefficient in poly])
    if len(reduced) <= 1 or reduced[-1] == 0:
        return None
    inverse = pow(reduced[-1], -1, prime)
    remaining = [(coefficient * inverse) % prime for coefficient in reduced]
    derivative = [
        (degree * coefficient) % prime
        for degree, coefficient in enumerate(remaining[1:], start=1)
    ]
    if len(fp_gcd(remaining, fp_trim(derivative), prime)) != 1:
        return None
    degrees: list[int] = []
    degree = 1
    x = [0, 1]
    while 2 * degree <= len(remaining) - 1:
        frobenius = fp_pow_mod(x, prime**degree, remaining, prime)
        block = fp_gcd(remaining, fp_sub(frobenius, x, prime), prime)
        block_degree = len(block) - 1
        if block_degree:
            assert block_degree % degree == 0
            degrees.extend([degree] * (block_degree // degree))
            quotient, remainder = fp_divmod(remaining, block, prime)
            assert remainder == [0]
            remaining = quotient
        degree += 1
    if len(remaining) > 1:
        degrees.append(len(remaining) - 1)
    assert sum(degrees) == len(reduced) - 1
    return sorted(degrees)


def capell_composition(P: list[F]) -> list[int]:
    composition = [F(0)] * 17
    for degree, coefficient in enumerate(P):
        composition[2 * degree] = coefficient / (2**degree)
    return int_primitive(composition)


def capell_evidence() -> dict[str, Any]:
    samples = (
        (F(-3), F(1, 2)),
        (F(-1), F(1)),
        (F(1), F(2)),
        (F(3), F(3, 2)),
        (F(5), F(4)),
    )
    rows = []
    for a, Z in samples:
        b_probe = F(1)
        _, _, _, _, _, P, _ = canonical_data(a, Z, b_probe)
        composition = capell_composition(P)
        certificate_prime = None
        certificate_degrees = None
        for prime in h10q.primerange(3, IRRED_PRIME_LIMIT + 1):
            PACER.tick()
            degrees = fp_factor_degrees(composition, prime)
            if degrees is not None and max(degrees) > 8:
                certificate_prime = prime
                certificate_degrees = degrees
                assert (degrees == [16]) == fp_irreducible(composition, prime)
                break
        rows.append(
            {
                "a": frac_text(a), "Z": frac_text(Z),
                "v2_Z": h10q.vp(Z, 2),
                "capell_exclusion_prime": certificate_prime,
                "factor_degrees_mod_prime": certificate_degrees,
            }
        )
    return {
        "type": "capell-composition-scan", "label": "EVIDENCE",
        "composition": "primitive integer scalar multiple of P(u^2/2)",
        "rows": rows,
        "prime_limit": IRRED_PRIME_LIMIT,
        "scope": "a factor degree above 8 excludes a rational 8+8 factor at that specialization; missing certificates are refusals, not negative evidence",
    }


def build_report(payload: list[dict[str, Any]], elapsed: float) -> str:
    search = next(row for row in payload if row["type"] == "bounded-rational-search")
    capell = next(row for row in payload if row["type"] == "capell-composition-scan")
    reciprocal = next(row for row in payload if row["type"] == "reciprocal-capell")
    lines = [
        "# L23 square-class families on the canonical branch",
        "",
        "## Status and scope",
        "",
        "- **PROVED:** on every `Phi`-admissible canonical specialization, `v2(P)` is divisible by four. Hence neither `P=2b*square` nor `P=-2b*square` can occur.",
        "- **PROVED:** `P=A*Norm_{Q(sqrt(2b))/Q}` cannot be a splitting mechanism on `Phi`, because it forces `(P,2b)_2=(A,2b)_2=-1`.",
        "- **PROVED:** the tempting `L=-A*rho^2` locus really enters `Phi` and factors `H` into reciprocal quartics, but it does not force all-plus signs; a broad coprime subfamily has a forced `p=3 mod 4` obstruction.",
        "- **PROVED:** on the reciprocal slice `a=1`, the algebraic-root/Capell condition `2*beta in K^2` is impossible for every rational nonzero `Z`.",
        "- **OPEN:** the full `(a,Z)` Capell square locus away from `a=1`. The modular rows below are only evidence.",
        "",
        "Throughout,",
        "",
        "```text",
        "A=1+4a^2,  s=(a-1)/2,  D=1-Z-a^2Z^2,",
        "Ng=16a^4b^2-A(b-1)^4,  L=1-2As^2b,",
        "X=a^4Z^2Ng,  Y=2ADb^2,",
        "H=(A/4)P=X^2+A L Y^2,  d=2b.",
        "```",
        "",
        "`Phi` means `a,b` are 2-adic units, equivalently `s in Z_2` and `v2(b)=0`. Finite searches are labelled **EVIDENCE** and are never used as negative evidence.",
        "",
        "## 1. Dyadic eliminations (PROVED)",
        "",
        "For a 2-adic unit `b`, `b-1 in 2 Z_2`, hence `Ng in 16 Z_2`. Put `m=v2(Z)`. Exact dominance gives",
        "",
        "```text",
        "v2(D)=0                  (m>=0),",
        "v2(D)=2m                 (m<0),",
        "v2(X^2)>=4m+8,",
        "v2(A L Y^2)=2+4 min(m,0).",
        "```",
        "",
        "The second summand is uniquely minimal, so",
        "",
        "```text",
        "v2(P)=4+4 min(m,0) == 0 (mod 4).",
        "```",
        "",
        "But `v2(plus_or_minus 2b*y^2)=1 mod 2`. Therefore both proposed curves `P=2b*y^2` and `P=-2b*y^2` have **no Phi-admissible rational point**. This is an exact local no-go, not a scan.",
        "",
        "For the broader norm proposal, if `P=A*N` with `N` a nonzero norm from `Q(sqrt(d))`, then",
        "",
        "```text",
        "(P,d)_2=(A,d)_2=-1.",
        "```",
        "For fixed data the norm equation is a genus-zero conic `R^2-dS^2=P/A`. Whenever that conic has a point, the desired quaternion has dyadic symbol `-1`; the norm condition forces failure rather than splitting.",
        "",
        "Indeed `A=5 mod 8`, `d=2*(odd unit)`, and the exact dyadic Hilbert formula has exponent `(A^2-1)/8=1 mod 2`; the mixed unit term vanishes because `(A-1)/2` is even. Thus the `P=A*norm` condition is itself incompatible with the desired split.",
        "",
        "## 2. The exact square curves (PROVED)",
        "",
        "For fixed `(a,Z)`, `y^2=P(b)` is generically a degree-eight hyperelliptic curve of genus `3`; `Y^2=plus_or_minus 2bP(b)` is generically degree nine of genus `4`. Section 1 eliminates the latter two Phi loci entirely.",
        "",
        "For fixed `(a,b)`, write",
        "",
        "```text",
        "P(Z)=C(1-Z-a^2Z^2)^2+E Z^4,",
        "C=16A^2b^4L,  E=(4a^8/A)Ng^2.",
        "```",
        "",
        "A direct `7 x 7` Sylvester determinant factors symbolically as",
        "",
        "```text",
        "Disc_Z(P)=16 C^3 E^2 (A^2 C+16E),",
        "A^2 C+16E=(16/A)(A^5b^4L+4a^8Ng^2).",
        "```",
        "",
        "On `Phi`, `C!=0`; `Ng=0` would make the 2-adic nonsquare `A` a rational square; the `Z^4` coefficient is nonzero by the same dyadic dominance; and the last discriminant parenthesis is a 2-adic unit plus a multiple of four. Hence the quartic has degree four and its discriminant never vanishes. The `P`-square locus in `Z` is therefore a smooth genus-one curve for every fixed Phi pair: the hoped-for genus-zero degeneration is **proved absent**.",
        "",
        "There is an exact off-cell split locus at `a=b=1`:",
        "",
        "```text",
        reciprocal["point_factor_identities"][0],
        reciprocal["point_factor_identities"][1],
        "```",
        "",
        "The first identity gives `Z=-5/2,-5/4,5/7,5/9`; the second gives `Z=-5/8,5/13`. At all six points `P` is the displayed square and the quaternion splits. None is a rational cube, so these are algebraic-family points and counterexamples to blanket emptiness, **not** members for the actual substitution `Z=z^3`.",
        "",
        "## 3. `L` square, twice-square, and fixed `b` squareclasses (PROVED)",
        "",
        "For `s!=0`, `L=ell^2` gives the rational parameter `b=(1-ell^2)/(2As^2)`. If `e=v2(s)`, Phi compatibility is exactly",
        "",
        "```text",
        "e>=1, ell in Z_2^x, and v2(1-ell^2)=1+2e.",
        "```",
        "",
        "For example, `a=1+4u` with `u` odd and `b=-1/(A u^2)` has `L=9` and lies in Phi. It does **not** force splitting: at `(a,Z,b)=(-3,1,-1/37)`, the exact symbol `(P,2b)_37=-1`.",
        "",
        "By contrast, `L=2ell^2` is never Phi-compatible. If `ell` is integral its numerator `1-2ell^2` is odd; if `v2(ell)<0`, subtracting `v2(2As^2)` still leaves an even nonzero valuation for `b`. In neither case is `v2(b)=0`.",
        "",
        "For the fixed odd squareclass menu `c in " + str(SQUARECLASS_MENU) + "`, put `b=c r^2`. The simultaneous condition `L=ell^2` is the genus-zero conic",
        "",
        "```text",
        "ell^2+2As^2 c r^2=1,",
        "ell=(1-2As^2c t^2)/(1+2As^2c t^2),",
        "r=2t/(1+2As^2c t^2).",
        "```",
        "",
        "For `s!=0` it meets Phi iff `v2(s)>=1`; then the two parameter ranges `v2(t)=-1` and `v2(t)=-2v2(s)` give unit `r`. For `s` odd it has no unit-`r` point. On the separate reciprocal slice `s=0`, `L=1` automatically and any unit `r` is allowed. This is an exact elimination of the squareclass condition, but, as the explicit `a=-3` row shows, it is not a splitting identity.",
        "",
        "## 4. The Phi-compatible factor family `L=-A*rho^2` (PROVED)",
        "",
        "Here",
        "",
        "```text",
        "b=(1+A rho^2)/(2As^2),",
        "H=(X-A rho Y)(X+A rho Y).",
        "```",
        "",
        "If `s,rho` are 2-adic units, then `1+A rho^2=6 mod 8`, so its valuation is exactly one and `b` is a unit. This genuinely crosses the L8c 2-adic wall. Each factor is a reciprocal quartic. With `u=b+b^-1`,",
        "",
        "```text",
        "F_eps/b^2=a^4Z^2(16a^4-A(u-2)^2)+eps*2A^2rho D,",
        "disc_u(F_eps/b^2)=8Aa^4Z^2(8a^8Z^2+eps*A^2rho D).",
        "```",
        "",
        "Thus the reciprocal-quadratic further-factor locus is the explicit genus-zero conic",
        "",
        "```text",
        "v^2=2A(8a^8Z^2+eps*A^2rho(1-Z-a^2Z^2)).",
        "```",
        "",
        "Factorization does not force all-plus signs. For odd integral `a=3 mod 4` and odd integral `rho`, put `N=(1+A rho^2)/2`. Then `N=3 mod 4`, so some `p=3 mod 4` divides `N` oddly. If `gcd(N,s)=1` and `Z` is a `p`-unit at one such prime, then `v_p(b)` is odd and",
        "",
        "```text",
        "P(b)=4a^8AZ^4 mod p,",
        "(A|p)=(-1|p)=-1,",
        "(P,2b)_p=-1.",
        "```",
        "",
        "The exact row `(a,rho,Z,b)=(3,1,1,19/37)` realizes the forced prime `p=19`. Scope is essential: `(a,rho,Z,b)=(-1,1,27,3/5)` is an exact globally split actual-cell specialization. Therefore this locus may align at selected cells, but it is not an identity killing all value-prime signs.",
        "",
        "## 5. Reciprocal Capell locus (PROVED empty)",
        "",
        "Set `a=1`, so `A=5,s=0`, and let `beta` be a root of the irreducible reciprocal octic. With `u=beta+beta^-1`,",
        "",
        "```text",
        reciprocal["trace_identity"],
        "```",
        "",
        "Let `K=Q(beta)` and `E=Q(u)`. If `y^2=2beta` in `K`, the involution `beta -> beta^-1` sends `y` to `plus_or_minus 2/y`. Accordingly one of",
        "",
        "```text",
        "(y+2/y)^2=2(u+2),",
        "(y-2/y)^2=2(u-2)",
        "```",
        "",
        "is a square in `E`. Taking norms gives the necessary rational-square conditions",
        "",
        "```text",
        "125D^2+1024Z^4 is a square, or",
        "125D^2+64Z^4 is a square.",
        "```",
        "Each displayed condition is itself a smooth genus-one quartic in `Z`: its discriminant is the nonzero specialization `16 C^3 E^2(25C+16E)` with `(C,E)=(125,1024)` or `(125,64)`.",
        "",
        "Write `Z=p/q` in lowest terms and `d0=q^2-pq-p^2`. In all three primitive parity classes `d0` is odd, while both cleared numerators",
        "",
        "```text",
        "125d0^2+1024p^4,  125d0^2+64p^4",
        "```",
        "",
        "are `5 mod 8`. Neither is a square. Hence `2beta` is **never** a square in `K` on the entire reciprocal slice. This concerns the algebraic root `beta`, not a rational Phi value and not L8c.",
        "",
        "## 6. Bounded discovery evidence",
        "",
        f"The squareclass scan checked {search['inputs_checked']} exact `(a,c,r,Z)` rows and {search['cube_inputs_checked']} cube rows. It found {len(search['hits'])} recorded structured hits and {len(search['cube_P_square_hits'])} `P`-square cube hits. Empty subsets are not promoted to theorems; Sections 1-5 contain the proofs.",
        "",
        "Degree-16 modular certificates for `P(u^2/2)`:",
        "",
        "| a | Z | v2(Z) | certificate prime | factor degrees mod p |",
        "|---:|---:|---:|---:|---|",
    ]
    for row in capell["rows"]:
        certificate = row["capell_exclusion_prime"]
        degrees = row["factor_degrees_mod_prime"]
        degree_text = ",".join(str(degree) for degree in degrees) if degrees else "REFUSAL"
        lines.append(
            f"| {row['a']} | {row['Z']} | {row['v2_Z']} | "
            f"{certificate if certificate is not None else 'REFUSAL'} | {degree_text} |"
        )
    lines.extend(
        [
            "",
            "A factor of degree above eight modulo a good prime cannot come from either member of a rational `8+8` factorization, so each displayed certificate rules out the Capell factor at that specialization. `REFUSAL` means only that no prime below the frozen bound certified it.",
            "",
            "## Machine artifacts",
            "",
            f"- JSONL: `{OUT}`",
            f"- Replay: `nice -n 19 python3 {Path(__file__).as_posix()}`",
            f"- Elapsed: `{elapsed:.3f}` seconds; pacer sleeps: `{PACER.sleeps}`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    started = time.monotonic()
    payload: list[dict[str, Any]] = []
    payload.append(
        {
            "type": "metadata", "label": "PROVED/EVIDENCE/OPEN",
            "branch": "tau_dagger=(1+2a^2)/(1+4a^2)",
            "squareclass_menu": list(SQUARECLASS_MENU),
            "stdlib_only": True,
        }
    )
    payload.append(symbolic_branch_identities())
    payload.append(symbolic_discriminant())
    payload.append(dyadic_replay())
    reciprocal, roots = reciprocal_symbolic()
    payload.append(reciprocal)
    payload.extend(structured_family_rows())
    payload.append(rational_search(roots))
    payload.append(capell_evidence())
    elapsed = time.monotonic() - started
    payload.append(
        {
            "type": "summary", "label": "PROVED",
            "proved_new_no_go": [
                "P=plus_or_minus 2b square is Phi-empty",
                "P=A times a norm from Q(sqrt(2b)) forces nonsplitting at 2",
                "reciprocal a=1 Capell locus is empty",
                "P-square quartic in Z is never genus-zero on Phi",
            ],
            "open": "full Capell square locus away from a=1",
            "elapsed_seconds": elapsed,
            "pacer_sleeps": PACER.sleeps,
        }
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in payload))
    REPORT.write_text(build_report(payload, elapsed))
    print("L23 squareclass replay: PASS")
    print(f"  records={len(payload)} output={OUT}")
    print(f"  report={REPORT}")
    print("  PROVED: dyadic P=±2b square no-go; P=A*norm no-go")
    print("  PROVED: reciprocal a=1 Capell locus empty; L=-A*rho^2 scoped obstruction")
    print("  OPEN: full Capell locus away from a=1")


if __name__ == "__main__":
    main()
