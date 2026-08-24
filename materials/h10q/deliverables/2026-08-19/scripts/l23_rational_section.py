#!/usr/bin/env python3
"""Exact rational-section analysis for the tau_dagger tied conic bundle.

This script is stdlib-only.  It proves identities and valuation no-go results;
it does not promote finite searches to theorems.  The first admissible direct
factor ansatz is normalized to a genus-4 hyperelliptic curve.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import reduce
import json
from math import gcd
from pathlib import Path
import time
from typing import Iterable


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l23_rational_section.jsonl"
REPORT = Path("/tmp/l23_rational_section.md")

# Sparse exact polynomials use the variables (a, Z, b, q, y, r).
VARIABLES = ("a", "Z", "b", "q", "y", "r")
A_I, Z_I, B_I, Q_I, Y_I, R_I = range(len(VARIABLES))
Monomial = tuple[int, ...]
Polynomial = dict[Monomial, F]


def constant(value: int | F = 1, dimension: int = len(VARIABLES)) -> Polynomial:
    coefficient = F(value)
    return {(0,) * dimension: coefficient} if coefficient else {}


def variable(index: int, dimension: int = len(VARIABLES)) -> Polynomial:
    exponent = [0] * dimension
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def add(*polynomials: Polynomial) -> Polynomial:
    answer: defaultdict[Monomial, F] = defaultdict(F)
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            answer[monomial] += coefficient
    return {monomial: coefficient for monomial, coefficient in answer.items()
            if coefficient}


def scale(polynomial: Polynomial, scalar: int | F) -> Polynomial:
    scalar_f = F(scalar)
    return {monomial: scalar_f * coefficient
            for monomial, coefficient in polynomial.items()
            if scalar_f * coefficient}


def multiply(*polynomials: Polynomial) -> Polynomial:
    if not polynomials:
        raise ValueError("multiply needs at least one polynomial")
    answer = constant(1, len(next(iter(polynomials[0]))))
    for polynomial in polynomials:
        product: defaultdict[Monomial, F] = defaultdict(F)
        for left_monomial, left_coefficient in answer.items():
            for right_monomial, right_coefficient in polynomial.items():
                monomial = tuple(x + y for x, y in zip(left_monomial,
                                                        right_monomial))
                product[monomial] += left_coefficient * right_coefficient
        answer = {monomial: coefficient
                  for monomial, coefficient in product.items() if coefficient}
    return answer


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("negative polynomial exponent")
    answer = constant(1, len(next(iter(polynomial), (0,) * len(VARIABLES))))
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, factor)
        factor = multiply(factor, factor)
        remaining //= 2
    return answer


def degree(polynomial: Polynomial, index: int) -> int:
    return max((monomial[index] for monomial in polynomial), default=-1)


def coefficient(polynomial: Polynomial, index: int,
                target_degree: int) -> Polynomial:
    answer: Polynomial = {}
    for monomial, value in polynomial.items():
        if monomial[index] != target_degree:
            continue
        reduced = list(monomial)
        reduced[index] = 0
        answer[tuple(reduced)] = value
    return answer


def weighted_initial(polynomial: Polynomial,
                     weights: dict[int, int]) -> tuple[int, Polynomial]:
    values = {
        monomial: sum(monomial[index] * weight
                      for index, weight in weights.items())
        for monomial in polynomial
    }
    minimum = min(values.values())
    return minimum, {
        monomial: coefficient_value
        for monomial, coefficient_value in polynomial.items()
        if values[monomial] == minimum
    }


def specialize_univariate(polynomial: Polynomial, index: int,
                          values: dict[int, int | F]) -> list[F]:
    coefficients: defaultdict[int, F] = defaultdict(F)
    for monomial, coefficient_value in polynomial.items():
        value = coefficient_value
        for variable_index, variable_value in values.items():
            value *= F(variable_value) ** monomial[variable_index]
        for variable_index, exponent in enumerate(monomial):
            if variable_index == index or variable_index in values:
                continue
            if exponent:
                raise AssertionError("unspecialized variable in univariate extraction")
        coefficients[monomial[index]] += value
    if not coefficients:
        return [F(0)]
    return trim([coefficients[exponent]
                 for exponent in range(max(coefficients) + 1)])


def trim(polynomial: list[F]) -> list[F]:
    answer = list(polynomial)
    while len(answer) > 1 and not answer[-1]:
        answer.pop()
    return answer


def univariate_derivative(polynomial: list[F]) -> list[F]:
    if len(polynomial) <= 1:
        return [F(0)]
    return trim([F(exponent) * polynomial[exponent]
                 for exponent in range(1, len(polynomial))])


def univariate_divmod(dividend: list[F], divisor: list[F]) -> tuple[list[F], list[F]]:
    remainder = trim(dividend)
    divisor_t = trim(divisor)
    if divisor_t == [0]:
        raise ZeroDivisionError
    quotient = [F(0)] * max(1, len(remainder) - len(divisor_t) + 1)
    while remainder != [0] and len(remainder) >= len(divisor_t):
        shift = len(remainder) - len(divisor_t)
        multiplier = remainder[-1] / divisor_t[-1]
        quotient[shift] = multiplier
        for index, coefficient_value in enumerate(divisor_t):
            remainder[index + shift] -= multiplier * coefficient_value
        remainder = trim(remainder)
    return trim(quotient), remainder


def univariate_gcd(left: list[F], right: list[F]) -> list[F]:
    first = trim(left)
    second = trim(right)
    while second != [0]:
        _, remainder = univariate_divmod(first, second)
        first, second = second, remainder
    if first == [0]:
        return first
    return [coefficient_value / first[-1] for coefficient_value in first]


def primitive_integer_coefficients(polynomial: list[F]) -> list[int]:
    denominator = 1
    for coefficient_value in polynomial:
        denominator = denominator * coefficient_value.denominator // gcd(
            denominator, coefficient_value.denominator)
    integers = [coefficient_value.numerator *
                (denominator // coefficient_value.denominator)
                for coefficient_value in polynomial]
    content = reduce(gcd, (abs(value) for value in integers if value), 0)
    if content:
        integers = [value // content for value in integers]
    if integers and integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def frac_text(value: F | int) -> str:
    value_f = F(value)
    return (str(value_f.numerator) if value_f.denominator == 1
            else str(value_f))


def build_bundle() -> dict[str, Polynomial]:
    one = constant()
    a = variable(A_I)
    Z = variable(Z_I)
    b = variable(B_I)
    y = variable(Y_I)
    r = variable(R_I)

    A = add(one, scale(power(a, 2), 4))
    s = scale(add(a, scale(one, -1)), F(1, 2))
    D = add(one, scale(Z, -1),
            scale(multiply(power(a, 2), power(Z, 2)), -1))
    Ng = add(scale(multiply(power(a, 4), power(b, 2)), 16),
             scale(multiply(A, power(add(b, scale(one, -1)), 4)), -1))
    X = multiply(power(a, 4), power(Z, 2), Ng)
    Y = scale(multiply(A, D, power(b, 2)), 2)
    L = add(one, scale(multiply(A, power(s, 2), b), -2))
    H = add(power(X, 2), multiply(A, L, power(Y, 2)))

    U = multiply(power(a, 2), y)
    V = scale(r, 2)
    denominator = multiply(power(A, 3), power(D, 2), power(b, 4))
    norm_form = add(power(U, 2), scale(multiply(b, power(V, 2)), -2))
    cleared_norm_equation = add(multiply(denominator, norm_form), scale(H, -1))

    original_cleared = add(
        scale(multiply(power(a, 8), power(Z, 4), power(Ng, 2)), -4),
        scale(multiply(power(a, 4), denominator, power(y, 2)), 4),
        scale(multiply(power(A, 4), power(s, 2), power(D, 2), power(b, 5)), 32),
        scale(multiply(power(A, 3), power(D, 2), power(b, 5), power(r, 2)), -32),
        scale(denominator, -16),
    )
    assert original_cleared == scale(cleared_norm_equation, 4)

    endpoint = multiply(power(a, 8), power(A, 2), power(Z, 4))
    assert degree(H, B_I) == 8
    assert coefficient(H, B_I, 0) == endpoint
    assert coefficient(H, B_I, 8) == endpoint

    left_weight, left_initial = weighted_initial(H, {Z_I: 1, B_I: 1})
    expected_left = add(
        multiply(power(a, 8), power(A, 2), power(Z, 4)),
        scale(multiply(power(A, 3), power(b, 4)), 4),
    )
    assert left_weight == 4 and left_initial == expected_left

    horizontal_weight, horizontal_initial = weighted_initial(H, {Z_I: 1})
    expected_horizontal = scale(
        multiply(power(A, 3), power(b, 4), L), 4)
    assert horizontal_weight == 0 and horizontal_initial == expected_horizontal

    return {
        "one": one, "a": a, "Z": Z, "b": b, "y": y, "r": r,
        "A": A, "s": s, "D": D, "Ng": Ng, "X": X, "Y": Y,
        "L": L, "H": H, "denominator": denominator,
        "cleared_norm_equation": cleared_norm_equation,
        "left_initial": left_initial,
        "horizontal_initial": horizontal_initial,
    }


def normalization_record(bundle: dict[str, Polynomial]) -> dict[str, object]:
    return {
        "type": "normalization",
        "label": "PROVED",
        "proof_mode": "exact sparse polynomial identity over Q",
        "definitions": {
            "A": "1+4*a^2",
            "D": "1-Z-a^2*Z^2",
            "s": "(a-1)/2",
            "Ng": "16*a^4*b^2-A*(b-1)^4",
            "X": "a^4*Z^2*Ng",
            "Y": "2*A*D*b^2",
            "L": "1-2*A*s^2*b",
            "H": "X^2+A*L*Y^2",
        },
        "exact_conic_normalization": (
            "A^3*D^2*b^4*((a^2*y)^2-2*b*(2*r)^2)=H"
        ),
        "original_equation_multiplier": "A^3*D^2*b^4",
        "original_equation_clears_to": (
            "4*(A^3*D^2*b^4*((a^2*y)^2-2*b*(2*r)^2)-H)=0"
        ),
        "H_degree_in_b": degree(bundle["H"], B_I),
        "H_equal_endpoints": "a^8*A^2*Z^4 at b^0 and b^8",
    }


def generic_fiber_record() -> dict[str, object]:
    return {
        "type": "generic-b-fiber-no-section",
        "label": "PROVED",
        "field": "k=Q(a,Z), generic base k(b)",
        "equation": "U^2-2*b*V^2=h(b), h=H/(A^3*D^2*b^4)",
        "b_zero_data": {
            "valuation_h": -4,
            "leading_unit": "(a^4*Z^2/D)^2/A",
            "square_class": "1/A",
        },
        "valuation_proof": (
            "At b=0, ord(U^2) is even and ord(2*b*V^2) is odd, so the "
            "two leading valuations cannot cancel.  Since ord(h)=-4, ord(U)=-2 "
            "and the residue equation would make 1/A a square in Q(a,Z).  "
            "The divisor A=1+4*a^2 occurs to odd order, so this is impossible."
        ),
        "b_infinity_data": {
            "valuation_h": -4,
            "leading_unit": "(a^4*Z^2/D)^2/A",
            "reason": "the b^0 and b^8 coefficients of H are equal",
        },
        "conclusion": (
            "The generic conic over Q(a,Z)(b) has no rational point.  In "
            "particular no Laurent-polynomial, polynomial, or rational-function "
            "section with b as the bundle coordinate exists; every degree-one "
            "Mobius pullback retains a rational simple zero/pole and the same obstruction."
        ),
        "scope": (
            "This does not rule out a ramified nonlinear pullback b=b(t), nor an "
            "isolated rational b in a fixed cell."
        ),
    }


def split_pullback_record() -> dict[str, object]:
    return {
        "type": "split-quadratic-pullback",
        "label": "PROVED",
        "parameter": "t",
        "substitution": "b=t^2/2",
        "section": {
            "h": "H/(A^3*D^2*b^4) evaluated at b=t^2/2",
            "y": "(h+1)/(2*a^2)",
            "r": "(1-h)/(4*t)",
        },
        "identity": (
            "((h+1)/2)^2-t^2*((1-h)/(2*t))^2=h"
        ),
        "two_adic_wall": "v2(b)=2*v2(t)-1 is odd for every rational t!=0",
        "admissible_members": 0,
        "conclusion": (
            "The first quadratic line construction gives a genuine uniform "
            "section only off the Phi condition v2(b)=0."
        ),
    }


def direct_factor_records(bundle: dict[str, Polynomial]) -> tuple[
        list[dict[str, object]], dict[int, Polynomial]]:
    one = bundle["one"]
    a = bundle["a"]
    Z = bundle["Z"]
    q = variable(Q_I)
    A = bundle["A"]
    D = bundle["D"]
    am1 = add(a, scale(one, -1))
    T = add(one, multiply(A, power(q, 2)))
    R = multiply(A, power(am1, 2))
    U0 = add(scale(T, 2), scale(R, -1))
    K = add(
        scale(multiply(power(a, 4), A, power(am1, 4), power(T, 2)), 64),
        scale(power(U0, 4), -1),
    )

    # b=2T/R.  This is the cleared numerator identity for Ng after substitution.
    ng_numerator = add(
        scale(multiply(power(a, 4), power(T, 2), power(R, 2)), 64),
        scale(multiply(A, power(U0, 4)), -1),
    )
    assert ng_numerator == multiply(A, K)
    assert add(one, scale(T, -1)) == scale(multiply(A, power(q, 2)), -1)

    # Verify the complete H factorization after substitution, over denominator R^8.
    X_numerator = multiply(power(a, 4), power(Z, 2), A, K)
    Y_numerator = scale(multiply(A, D, power(T, 2)), 8)
    H_numerator = add(
        power(X_numerator, 2),
        multiply(A, scale(multiply(A, power(q, 2)), -1),
                 power(Y_numerator, 2), power(R, 4)),
    )
    minus_factor = add(X_numerator,
                       scale(multiply(A, q, Y_numerator, power(R, 2)), -1))
    plus_factor = add(X_numerator,
                      multiply(A, q, Y_numerator, power(R, 2)))
    assert H_numerator == multiply(minus_factor, plus_factor)

    rows: list[dict[str, object]] = []
    W_by_sign: dict[int, Polynomial] = {}
    for sign in (1, -1):
        C = scale(multiply(power(A, 3), power(am1, 4), q, power(T, 2)),
                  sign * 8)
        F_direct = add(multiply(power(a, 4), power(Z, 2), K),
                       scale(multiply(C, D), -1))
        expected_cleared_factor = add(
            X_numerator,
            scale(multiply(A, q, Y_numerator, power(R, 2)), -sign),
        )
        assert expected_cleared_factor == multiply(A, F_direct)

        quadratic_lead = add(multiply(power(a, 4), K),
                             multiply(power(a, 2), C))
        expanded_quadratic = add(
            multiply(quadratic_lead, power(Z, 2)), multiply(C, Z), scale(C, -1))
        assert F_direct == expanded_quadratic

        V = add(
            scale(multiply(power(a, 8), A, power(am1, 4), power(T, 2)), 64),
            scale(multiply(power(A, 4), power(am1, 4), q, power(T, 2)),
                  sign * 2),
            scale(multiply(power(a, 4), power(U0, 4)), -1),
        )
        discriminant = add(power(C, 2),
                           scale(multiply(quadratic_lead, C), 4))
        hyperelliptic_rhs = scale(multiply(A, q, V), sign * 2)
        square_factor = scale(multiply(A, power(am1, 2), T), 4)
        assert discriminant == multiply(power(square_factor, 2),
                                        hyperelliptic_rhs)
        W = multiply(q, V)
        assert degree(F_direct, Q_I) == 8
        assert degree(F_direct, Z_I) == 2
        assert max(monomial[Q_I] + monomial[Z_I]
                   for monomial in F_direct) == 10
        assert degree(W, Q_I) == 9

        W_at_3 = specialize_univariate(W, Q_I, {A_I: 3})
        gcd_certificate = univariate_gcd(
            W_at_3, univariate_derivative(W_at_3))
        assert gcd_certificate == [F(1)]
        primitive = primitive_integer_coefficients(W_at_3)
        W_by_sign[sign] = W

        rows.append({
            "type": "admissible-factor-residual",
            "label": "PROVED",
            "sign": sign,
            "ansatz": {
                "L": "-A*q^2",
                "b": "(1+A*q^2)/(2*A*s^2)=2*(1+A*q^2)/(A*(a-1)^2)",
                "factorization": "H=(X-A*q*Y)*(X+A*q*Y)",
                "zero_factor_equation": (
                    "a^4*Z^2*K-sign*8*A^3*(a-1)^4*q*T^2*D=0"
                ),
                "T": "1+A*q^2",
                "U0": "2*T-A*(a-1)^2",
                "K": "64*a^4*A*(a-1)^4*T^2-U0^4",
                "tied_point": "y=r=0 whenever the selected zero factor holds",
            },
            "two_adic_admissibility": (
                "If a=3 mod 4 and q is a 2-adic unit, then s is a unit, "
                "A=5 mod 8, T=6 mod 8, and v2(b)=1-1=0."
            ),
            "fixed_cell_degree_in_q": 8,
            "curve_over_Q(a)": {
                "coordinates": ["q", "Z"],
                "bidegree_in_P1xP1": [8, 2],
                "affine_total_degree": 10,
                "quadratic_in_Z": (
                    "(a^4*K+a^2*C)*Z^2+C*Z-C=0, "
                    "C=sign*8*A^3*(a-1)^4*q*T^2"
                ),
                "birational_hyperelliptic_model": "w^2=2*sign*A*q*V(q)",
                "V": (
                    "64*a^8*A*(a-1)^4*T^2+sign*2*A^4*(a-1)^4*q*T^2"
                    "-a^4*U0^4"
                ),
                "birational_w": (
                    "w=(2*(a^4*K+a^2*C)*Z+C)/(4*A*(a-1)^2*T)"
                ),
                "branch_polynomial_degree": 9,
                "geometric_genus": 4,
            },
            "squarefree_certificate": {
                "specialization": "a=3",
                "primitive_coefficients_ascending": primitive,
                "gcd_with_derivative": [1],
                "logical_use": (
                    "This exact specialization makes the generic discriminant "
                    "nonzero; it is a proof of generic squarefreeness, not a point hunt."
                ),
            },
            "uniform_section_conclusion": (
                "The normalization is a geometrically integral genus-4 curve, "
                "so this ansatz has no rational parametrization/section over the Z-line."
            ),
            "per_cell_scope": (
                "For fixed rational (a,Z), the degree-8 equation may still have "
                "isolated rational q; genus 4 does not decide those points.  Such "
                "a point only gives an admissible tied-conic b after the stated "
                "2-adic hypotheses, and it is not automatically a member of the "
                "preselected aligned congruence class."
            ),
        })
    return rows, W_by_sign


def capell_record(bundle: dict[str, Polynomial]) -> dict[str, object]:
    return {
        "type": "capell-square-criterion",
        "label": "PROVED",
        "generic_function_field": {
            "field": "Q(a,Z), beta a root of the already-proved irreducible P(beta)",
            "Z_adic_newton_edge": [[0, 4], [4, 0]],
            "residual_after_beta=Z*u": "4*A*(a^8+4*A*u^4)",
            "separability": "derivative 64*A^2*u^3; gcd is 1 in characteristic 0",
            "valuation": "an unramified prime has ord(beta)=1, hence ord(2*beta)=1",
            "conclusion": (
                "2*beta is not a square in Q(a,Z)[beta].  By Capell, "
                "P(u^2/2) is generically irreducible of degree 16, not an 8x8 split."
            ),
        },
        "actual_cell_target_prime": {
            "hypotheses": (
                "w odd, e=v_w(z)>=1, Z=z^3, (A|w)=-1, w does not divide a*A, D=1 mod w"
            ),
            "left_edge": [[0, "12*e"], [4, 0]],
            "left_residual": (
                "4*A*(a^8*(Z/w^(3e))^4+4*A*u^4), separable"
            ),
            "odd_e_conclusion": (
                "Four local roots have ord(beta)=3e; if e is odd then 2*beta "
                "has odd valuation and is not a square."
            ),
            "horizontal_edge_if_w_not_divide_s": [[4, 0], [5, 0]],
            "horizontal_residual": "4*A^3*b^4*(1-2*A*s^2*b)",
            "simple_root": "beta=1/(2*A*s^2)",
            "unit_square_class": "2*beta=1/(A*s^2), a nonsquare because (A|w)=-1",
            "all_e_conclusion_if_w_not_divide_s": True,
            "class_choice_refinement": (
                "There are (w-(-1|w))/2 >= 2 residues a with (A|w)=-1 "
                "(the polynomial A has 1+(-1|w) zero residues), and at most "
                "one nonsquare residue has a=1 mod w; CRT can therefore impose "
                "a!=1 mod w, i.e. w does not divide s."
            ),
        },
        "open_edge": (
            "For a pre-existing fixed class with w|s and even v_w(z), these "
            "two target-prime arguments alone are inconclusive.  No negative "
            "claim is made for that specialization."
        ),
        "relevance": (
            "The untwisted algebraic-root square criterion cannot be a uniform "
            "escape on the refined L20 class choice.  Algebraic beta is not the "
            "rational Phi variable b."
        ),
    }


def constant_twist_record() -> dict[str, object]:
    # Cross-multiplied parametrization identity in Q[a,theta,j].
    dimension = 3
    a = variable(0, dimension)
    theta = variable(1, dimension)
    j = variable(2, dimension)
    one = constant(1, dimension)
    A = add(one, scale(power(a, 2), 4))
    jt2 = multiply(j, power(theta, 2))
    assert add(power(add(jt2, scale(A, -1)), 2),
               scale(power(add(jt2, A), 2), -1)) == scale(
                   multiply(A, j, power(theta, 2)), -4)

    return {
        "type": "constant-twist-norm-no-go",
        "label": "PROVED",
        "setup": "2*b=j*q^2 with constant j in Q^x",
        "endpoint_necessity": (
            "At q=0 the pulled-back norm is U^2-j*(q*V)^2=h and "
            "ord_q(h)=-8.  Its leading residue has square class 1/A, so "
            "1/A, equivalently A, must be a norm from Q(a,Z)(sqrt(j))."
        ),
        "divisor_proof": (
            "At the irreducible divisor A=1+4*a^2, the residue field is Q(i).  "
            "A constant j is a square in Q(i) iff j is a rational square or "
            "minus a rational square.  If the divisor is inert, every norm has "
            "even A-valuation, while A has valuation one."
        ),
        "exact_constant_twists": (
            "For nonsquare j, A is a norm iff j=-1 modulo Q^x2; explicitly "
            "A=Norm_Q(i)(1+2*a*i).  Square j is the split case."
        ),
        "two_adic_wall": (
            "The only permitted constant square classes j=+1 or -1 have even "
            "v2(j), so v2(b)=v2(j)+2*v2(q)-1 is odd.  A twist with odd "
            "v2(j) would make b admissible but fails the necessary A-norm condition."
        ),
        "exact_factor_subroute": {
            "extra_condition": "A*L=-j*k^2",
            "norm_identity": "H=X^2-j*(k*Y)^2",
            "matching_conic": "k^2-A^2*s^2*q^2=-A/j",
            "parametrization": {
                "k": "(theta-A/(j*theta))/2",
                "q": "-(theta+A/(j*theta))/(2*A*s)",
            },
        },
        "conclusion": (
            "No constant-squareclass quadratic pullback can carry a rational "
            "section and satisfy v2(b)=0."
        ),
        "scope": (
            "This does not rule out nonconstant twists, higher-degree pullbacks, "
            "or isolated rational points."
        ),
    }


def two_adic_newton_record() -> dict[str, object]:
    return {
        "type": "two-adic-capell-newton-correction",
        "label": "PROVED",
        "hypotheses": "a odd, m=v2(Z)>=1, v2(s)=0, P=(4/A)*H",
        "coefficient_valuations_b0_to_b8": (
            "[4m+2,4m+5,4m+4,4m+5,4,5,4m+4,4m+5,4m+2]"
        ),
        "newton_hulls": {
            "m=1": [[0, 6], [4, 4], [8, 6]],
            "m>=2": [[0, "4m+2"], [4, 4], [5, 5], [8, "4m+2"]],
        },
        "correction": (
            "The b^5 perturbation -32*A^3*s^2*D^2 has valuation 5 in P "
            "and creates the extra vertex for m>=2."
        ),
        "m_nonpositive": (
            "The lower hull is horizontal and the first residual is (b+1)^8 "
            "mod 2; higher Newton analysis is required."
        ),
        "scope": (
            "Repeated residuals and possible higher ramification prevent a "
            "universal 2-adic nonsquare claim from these first polygons alone."
        ),
    }


def build_report(records: list[dict[str, object]], elapsed: float) -> str:
    residuals = [record for record in records
                 if record["type"] == "admissible-factor-residual"]
    coefficients = "\n".join(
        f"- sign {record['sign']:+d}: `{record['squarefree_certificate']['primitive_coefficients_ascending']}`"
        for record in residuals
    )
    return "\n".join([
        "# L23 rational-section analysis",
        "",
        "## Status",
        "",
        "**OPEN overall; exact no-go/reduction proved.**  No admissible uniform",
        "rational section was found.  The generic conic over `Q(a,Z)(b)` has no",
        "point, the first split quadratic pullback violates `v2(b)=0`, and the",
        "first genuinely admissible zero-factor construction reduces exactly to",
        "a geometrically integral genus-4 hyperelliptic curve.",
        "",
        "No finite point hunt is used as negative evidence.",
        "",
        "## 1. Exact normalization",
        "",
        "Put",
        "",
        "```text",
        "A=1+4a^2, D=1-Z-a^2Z^2, s=(a-1)/2,",
        "Ng=16a^4b^2-A(b-1)^4,",
        "X=a^4Z^2Ng, Y=2ADb^2, L=1-2As^2b,",
        "H=X^2+A L Y^2.",
        "```",
        "",
        "Exact denominator clearing gives",
        "",
        "```text",
        "A^3 D^2 b^4 ((a^2 y)^2-2b(2r)^2)=H.",
        "```",
        "",
        "This is exactly the equation in the assignment, multiplied by",
        "`A^3 D^2 b^4/4`.  The sparse replay proves the identity term by term.",
        "",
        "## 2. Degenerate endpoints: no section over the b-line",
        "",
        "At both `b=0` and `b=infinity`, `H` has endpoint coefficient",
        "`a^8 A^2 Z^4`.  Thus for `h=H/(A^3D^2b^4)` the endpoint unit is",
        "",
        "```text",
        "(a^4 Z^2/D)^2 / A.",
        "```",
        "",
        "At `b=0`, the valuations of `U^2` and `2bV^2` have opposite parity.",
        "They cannot cancel.  Matching `ord(h)=-4` forces the residue of `U^2`",
        "to have square class `1/A`, impossible in `Q(a,Z)` because the prime",
        "divisor `A=1+4a^2` has odd valuation.  The equal endpoint gives the",
        "same obstruction at infinity.  Therefore there is no point on the",
        "generic conic over `Q(a,Z)(b)`, hence no Laurent/polynomial/rational",
        "section using `b` as parameter and no degree-one Mobius pullback.",
        "",
        "This is a generic-fibre statement, not a denial of specialized points.",
        "",
        "## 3. First ramified line and its 2-adic wall",
        "",
        "The pullback `b=t^2/2` splits the norm form and gives",
        "",
        "```text",
        "y=(h+1)/(2a^2),   r=(1-h)/(4t).",
        "```",
        "",
        "It is an exact uniform section, but `v2(b)=2v2(t)-1` is always odd.",
        "So it has no admissible member.",
        "",
        "## 4. First admissible direct factor and residual curve",
        "",
        "Impose `L=-Aq^2`.  Then",
        "",
        "```text",
        "b=(1+Aq^2)/(2As^2),",
        "H=(X-AqY)(X+AqY).",
        "```",
        "",
        "For `a=3 mod 4` and 2-adic-unit `q`, this has `v2(b)=0`.  Setting",
        "one factor to zero makes `H=0`, so `y=r=0` is a tied-conic point.",
        "With `T=1+Aq^2`, `U0=2T-A(a-1)^2`, and",
        "`K=64a^4A(a-1)^4T^2-U0^4`, the two signs give",
        "",
        "```text",
        "a^4 Z^2 K - sign*8 A^3(a-1)^4 q T^2(1-Z-a^2Z^2)=0.",
        "```",
        "",
        "For fixed `(a,Z)` this has degree 8 in `q`.  For fixed `a` it is",
        "bidegree `(8,2)` in `(q,Z)` (affine total degree 10).  Completing the",
        "quadratic in `Z` gives the birational model",
        "",
        "```text",
        "w^2 = 2 sign A q V(q),",
        "V=64a^8A(a-1)^4T^2 + sign*2A^4(a-1)^4qT^2 - a^4U0^4.",
        "```",
        "",
        "The branch polynomial `qV` has degree 9.  Exact Euclidean gcd at",
        "`a=3` is 1 for both signs, proving the generic discriminant is nonzero",
        "(this specialization is a proof certificate, not a rational-point hunt).",
        "The smooth projective normalization therefore has 10 branch points",
        "including infinity, and Riemann-Hurwitz gives genus 4.",
        "",
        "Primitive `a=3` branch certificates (ascending coefficients):",
        "",
        coefficients,
        "",
        "A geometrically integral genus-4 curve has no rational parametrization,",
        "so this ansatz supplies no uniform section over the `Z`-line.  It does",
        "not rule out isolated rational `q` for a fixed cell.  Even such a point",
        "must separately lie in the preselected aligned congruence class; the",
        "factor equation alone proves only the tied-conic identity and `v2(b)=0`.",
        "",
        "## 5. Capell algebraic-root criterion",
        "",
        "Let `beta` be the algebraic root of the irreducible octic.  At the",
        "generic `Z`-adic place the first Newton edge is `(0,4)-(4,0)`, with",
        "separable residual `4A(a^8+4Au^4)`.  It yields an unramified prime",
        "with `ord(beta)=1`, hence `2beta` is not a square.  Capell therefore",
        "makes `P(u^2/2)` generically irreducible of degree 16, rather than an",
        "8-by-8 split.",
        "",
        "For an actual cell, write `e=v_w(z)`.  The target-prime left edge kills",
        "the square criterion when `e` is odd.  If `w` does not divide `s`, the",
        "horizontal edge has the simple root `beta=1/(2As^2)`, so exactly",
        "`2beta=1/(As^2)`, of nonsquare class `1/A`; this kills every `e`.",
        "L20's character choice can impose `a!=1 mod w` (hence `w` does not divide `s`) because",
        "at least two residues make `A` nonsquare.  For a pre-existing class with",
        "`w|s` and even `e`, this local argument is explicitly OPEN.",
        "",
        "The independent 2-adic first polygon must be scoped more carefully.",
        "For `m=v2(Z)>=1` and 2-adic-unit `s`, the coefficients of",
        "`P=(4/A)H` have valuations",
        "",
        "```text",
        "[4m+2,4m+5,4m+4,4m+5,4,5,4m+4,4m+5,4m+2].",
        "```",
        "",
        "Thus `m=1` has hull `(0,6)-(4,4)-(8,6)`, while `m>=2` has the",
        "extra perturbation vertex `(5,5)`.  For `m<=0` the first residual is",
        "the repeated `(b+1)^8` modulo 2.  These first polygons alone do not",
        "prove a 2-adic nonsquare statement; higher ramification is OPEN.",
        "",
        "The algebraic root `beta` is not the rational Phi variable `b`.",
        "",
        "## 6. Constant-squareclass quadratic pullbacks",
        "",
        "For `2b=jq^2`, the pulled-back norm field is the constant quadratic",
        "`Q(a,Z)(sqrt(j))`.  At `q=0`, the leading coefficient of `h` has",
        "square class `1/A`; hence `A` must be a norm from that constant field.",
        "At the divisor `A=0` the residue field is `Q(i)`, so this is possible",
        "only when `j` is a square or minus a square.  Both classes have even",
        "`v2(j)`, making `v2(b)=v2(j)+2v2(q)-1` odd.  Therefore no",
        "constant-squareclass quadratic pullback can be both section-bearing",
        "and admissible.  The tempting extra match `AL=-jk^2` is a rational",
        "subroute of this argument, not an escape.",
        "",
        "## 7. Scope ledger",
        "",
        "- **PROVED:** exact normalization and equal endpoints.",
        "- **PROVED:** no generic `Q(a,Z)(b)` point; all Laurent/rational",
        "  ansatzes for `y,r` over the generic `b`-line fail.",
        "- **PROVED:** split quadratic section exists but is 2-adically inadmissible.",
        "- **PROVED:** the admissible `L=-Aq^2`, `X=sign*AqY` route reduces to",
        "  a degree-8 fixed-cell equation / genus-4 hyperelliptic curve.",
        "- **PROVED:** generic and refined-cell Capell untwisted square criterion fails.",
        "- **PROVED:** every constant-squareclass quadratic pullback hits the",
        "  scalar-norm/2-adic wall.",
        "- **OPEN:** rational points on the genus-4 residual for every cell, and",
        "  any different nonlinear/nonconstant-twist section of the total surface.",
        "",
        "A rational section over `Q(a,Z)` would be one uniform identity.  A",
        "per-cell point may depend arithmetically on the fixed rational `(a,Z)`",
        "and is not supplied or excluded by the genus computation.  Nor does a",
        "point on the residual automatically satisfy the fixed aligned-class",
        "congruences required by the member-existence step.",
        "",
        f"Replay elapsed: {elapsed:.3f} seconds.",
        "",
    ])


def main() -> int:
    started = time.perf_counter()
    bundle = build_bundle()
    normalization = normalization_record(bundle)
    generic = generic_fiber_record()
    split = split_pullback_record()
    residuals, _ = direct_factor_records(bundle)
    capell = capell_record(bundle)
    twist = constant_twist_record()
    two_adic = two_adic_newton_record()
    elapsed = time.perf_counter() - started

    summary = {
        "type": "summary",
        "label": "OPEN",
        "uniform_admissible_section": None,
        "exact_reduction": "two genus-4 hyperelliptic curves, signs +/-",
        "fixed_cell_residual_degree_in_q": 8,
        "generic_residual_genus": 4,
        "schinzel_removed": False,
        "finite_point_hunts": 0,
        "proved_no_go_routes": [
            "generic b-line/Laurent/rational section",
            "split quadratic pullback under v2(b)=0",
            "all constant-squareclass quadratic pullbacks under v2(b)=0",
            "generic and refined-cell untwisted Capell square criterion",
        ],
        "remaining_open": [
            "per-cell rational points on the genus-4 residual",
            "intersection of those points with each preselected aligned class",
            "nonconstant-twist or different nonlinear sections of the total surface",
            "Capell specialization with w|s and even v_w(z) for a pre-existing class",
        ],
        "wall_seconds": elapsed,
    }
    meta = {
        "type": "meta",
        "label": "PROVED",
        "schema": "l23-rational-section-v1",
        "source_script": "math/h10q/l23_rational_section.py",
        "report": "/tmp/l23_rational_section.md",
        "stdlib_only": True,
        "refusals_are_negative_evidence": False,
        "finite_scans_imply_uniform_theorem": False,
    }
    records: list[dict[str, object]] = [
        meta, normalization, generic, split, *residuals, capell, twist,
        two_adic, summary,
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records) + "\n", encoding="utf-8")
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print(
        "L23 rational section: OPEN; generic b-section ruled out; "
        "admissible direct factor -> genus 4 (degree 8 per cell); "
        f"Capell/twist no-gos proved; {elapsed:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
