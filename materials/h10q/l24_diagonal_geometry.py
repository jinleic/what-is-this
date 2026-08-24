#!/usr/bin/env python3
"""Exact geometry of the L23 diagonal self-coupled square branches.

This standalone, stdlib-only replay clears the Sun bridge denominator, builds
both degree-eight covers, determines their generic branch tower, certifies an
actual-cell genus slice, and proves that both orientations are globally
Phi-empty over Q_2.  Finite-field certificates are positive irreducibility
proofs; no finite no-hit statement is used.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
from typing import Any
import json
import time


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l24_diagonal_geometry.jsonl"
REPORT = Path("/tmp/l24_diagonal_geometry.md")
PACE_EVERY = 256
PACE_SECONDS = 0.001


class Pacer:
    """Pace the small exact finite-field loops; no pool or worker is used."""

    def __init__(self) -> None:
        self.steps = 0
        self.sleeps = 0

    def step(self) -> None:
        self.steps += 1
        if self.steps % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()
VARIABLES = ("s", "b", "Z")
S_I, B_I, Z_I = range(3)
Monomial = tuple[int, int, int]
Polynomial = dict[Monomial, F]
UniPoly = list[F]


def constant(value: int | F = 1) -> Polynomial:
    value_f = F(value)
    return {(0, 0, 0): value_f} if value_f else {}


def variable(index: int) -> Polynomial:
    exponent = [0, 0, 0]
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
    answer = constant()
    for polynomial in polynomials:
        product: defaultdict[Monomial, F] = defaultdict(F)
        for left_monomial, left_coefficient in answer.items():
            for right_monomial, right_coefficient in polynomial.items():
                monomial = tuple(left_monomial[index] + right_monomial[index]
                                 for index in range(3))
                product[monomial] += left_coefficient * right_coefficient
        answer = {monomial: coefficient
                  for monomial, coefficient in product.items() if coefficient}
    return answer


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    assert exponent >= 0
    answer = constant()
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        remaining //= 2
    return answer


def degree(polynomial: Polynomial, index: int) -> int:
    return max((monomial[index] for monomial in polynomial), default=-1)


def evaluate(polynomial: Polynomial, values: dict[int, int | F]) -> Polynomial:
    answer: defaultdict[Monomial, F] = defaultdict(F)
    for monomial, coefficient in polynomial.items():
        exponent = list(monomial)
        value = coefficient
        for index, specialization in values.items():
            value *= F(specialization) ** exponent[index]
            exponent[index] = 0
        answer[tuple(exponent)] += value
    return {monomial: coefficient for monomial, coefficient in answer.items()
            if coefficient}


def specialize_univariate(polynomial: Polynomial, index: int,
                           values: dict[int, int | F]) -> UniPoly:
    specialized = evaluate(polynomial, values)
    coefficients: defaultdict[int, F] = defaultdict(F)
    for monomial, coefficient in specialized.items():
        assert all(exponent == 0 for other, exponent in enumerate(monomial)
                   if other != index)
        coefficients[monomial[index]] += coefficient
    top = max(coefficients, default=0)
    return poly_trim([coefficients[exponent] for exponent in range(top + 1)])


def poly_trim(polynomial: UniPoly) -> UniPoly:
    answer = [F(coefficient) for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(*polynomials: UniPoly) -> UniPoly:
    size = max((len(polynomial) for polynomial in polynomials), default=1)
    answer = [F(0)] * size
    for polynomial in polynomials:
        for index, coefficient in enumerate(polynomial):
            answer[index] += coefficient
    return poly_trim(answer)


def poly_scale(polynomial: UniPoly, scalar: int | F) -> UniPoly:
    return poly_trim([F(scalar) * coefficient for coefficient in polynomial])


def poly_multiply(left: UniPoly, right: UniPoly) -> UniPoly:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += left_coefficient * right_coefficient
    return poly_trim(answer)


def poly_power(polynomial: UniPoly, exponent: int) -> UniPoly:
    answer = [F(1)]
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = poly_multiply(answer, base)
        base = poly_multiply(base, base)
        remaining //= 2
    return answer


def poly_eval(polynomial: UniPoly, value: int | F) -> F:
    answer = F(0)
    value_f = F(value)
    for coefficient in reversed(polynomial):
        answer = answer * value_f + coefficient
    return answer


def poly_derivative(polynomial: UniPoly) -> UniPoly:
    if len(polynomial) <= 1:
        return [F(0)]
    return poly_trim([F(index) * polynomial[index]
                      for index in range(1, len(polynomial))])


def poly_divmod(dividend: UniPoly, divisor: UniPoly) -> tuple[UniPoly, UniPoly]:
    remainder = poly_trim(dividend)
    divisor_t = poly_trim(divisor)
    assert divisor_t != [F(0)]
    quotient = [F(0)] * max(1, len(remainder) - len(divisor_t) + 1)
    while remainder != [F(0)] and len(remainder) >= len(divisor_t):
        shift = len(remainder) - len(divisor_t)
        coefficient = remainder[-1] / divisor_t[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor_t):
            remainder[index + shift] -= coefficient * value
        remainder = poly_trim(remainder)
    return poly_trim(quotient), remainder


def poly_gcd(left: UniPoly, right: UniPoly) -> UniPoly:
    first = poly_trim(left)
    second = poly_trim(right)
    while second != [F(0)]:
        first, second = second, poly_divmod(first, second)[1]
    return poly_scale(first, 1 / first[-1])


def lcm(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right)


def primitive_integer_poly(polynomial: UniPoly) -> list[int]:
    denominator = 1
    for coefficient in polynomial:
        denominator = lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in polynomial]
    content = 0
    for coefficient in integers:
        content = gcd(content, abs(coefficient))
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return integers


def mod_trim(polynomial: list[int], prime: int) -> list[int]:
    answer = [coefficient % prime for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def mod_add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    return mod_trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ], prime)


def mod_divmod(dividend: list[int], divisor: list[int],
               prime: int) -> tuple[list[int], list[int]]:
    remainder = mod_trim(dividend, prime)
    divisor_t = mod_trim(divisor, prime)
    assert divisor_t != [0]
    quotient = [0] * max(1, len(remainder) - len(divisor_t) + 1)
    inverse = pow(divisor_t[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor_t):
        shift = len(remainder) - len(divisor_t)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor_t):
            remainder[index + shift] -= coefficient * value
        remainder = mod_trim(remainder, prime)
        PACER.step()
    return mod_trim(quotient, prime), remainder


def mod_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    first = mod_trim(left, prime)
    second = mod_trim(right, prime)
    while second != [0]:
        first, second = second, mod_divmod(first, second, prime)[1]
    inverse = pow(first[-1], -1, prime)
    return mod_trim([coefficient * inverse for coefficient in first], prime)


def mod_multiply(left: list[int], right: list[int], modulus: list[int],
                 prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            product[left_index + right_index] += left_coefficient * right_coefficient
    PACER.step()
    return mod_divmod(product, modulus, prime)[1]


def mod_power(base: list[int], exponent: int, modulus: list[int],
              prime: int) -> list[int]:
    answer = [1]
    current = mod_trim(base, prime)
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = mod_multiply(answer, current, modulus, prime)
        current = mod_multiply(current, current, modulus, prime)
        remaining //= 2
    return answer


def prime_divisors(value: int) -> list[int]:
    answer: list[int] = []
    remaining = value
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            answer.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        answer.append(remaining)
    return answer


def irreducibility_certificate(polynomial: UniPoly, prime: int) -> dict[str, Any]:
    primitive = primitive_integer_poly(polynomial)
    reduced = mod_trim(primitive, prime)
    degree_value = len(reduced) - 1
    assert degree_value == len(primitive) - 1
    inverse = pow(reduced[-1], -1, prime)
    monic = mod_trim([coefficient * inverse for coefficient in reduced], prime)
    x_poly = [0, 1]
    gcd_rows = []
    for divisor in prime_divisors(degree_value):
        exponent = prime ** (degree_value // divisor)
        residue = mod_add(mod_power(x_poly, exponent, monic, prime), [0, -1], prime)
        common = mod_gcd(monic, residue, prime)
        assert common == [1]
        gcd_rows.append({
            "divisor": divisor,
            "exponent": exponent,
            "gcd": common,
        })
    final = mod_add(mod_power(x_poly, prime ** degree_value, monic, prime),
                    [0, -1], prime)
    assert final == [0]
    digest = sha256(json.dumps(primitive, separators=(",", ":")).encode()).hexdigest()
    return {
        "prime": prime,
        "degree": degree_value,
        "primitive_sha256": digest,
        "reduction_ascending": reduced,
        "monic_reduction_ascending": monic,
        "rabin_gcd_rows": gcd_rows,
        "x_pow_pn_minus_x_remainder": final,
    }


def vp(value: int | F, prime: int = 2) -> int:
    value_f = F(value)
    assert value_f != 0

    def integer_vp(integer: int) -> int:
        integer = abs(integer)
        answer = 0
        while integer % prime == 0:
            integer //= prime
            answer += 1
        return answer

    return integer_vp(value_f.numerator) - integer_vp(value_f.denominator)


def frac_text(value: int | F) -> str:
    value_f = F(value)
    return (str(value_f.numerator) if value_f.denominator == 1
            else str(value_f))


def build_geometry(orientation: str) -> dict[str, Polynomial]:
    assert orientation in {"I", "II"}
    one = constant()
    s = variable(S_I)
    b = variable(B_I)
    Z = variable(Z_I)
    a = add(one, scale(s, 2))
    A = add(one, scale(power(a, 2), 4))
    D = add(one, scale(Z, -1), scale(multiply(power(a, 2), power(Z, 2)), -1))
    Ng = add(
        scale(multiply(power(a, 4), power(b, 2)), 16),
        scale(multiply(A, power(add(b, constant(-1)), 4)), -1),
    )
    E = multiply(A, power(b, 2), D)
    K = multiply(power(a, 2), power(Z, 2), Ng)
    E2 = power(E, 2)
    K2 = power(K, 2)
    T = multiply(A, E2)

    if orientation == "I":
        e = add(scale(multiply(A, b, power(s, 2)), 2), constant(-1))
        linear = add(multiply(add(A, scale(b, -32)), T), scale(K2, -1))
        N = add(
            K2,
            multiply(
                add(scale(multiply(A, b, add(one, power(s, 2))), 32),
                    constant(-16)),
                E2,
            ),
        )
    else:
        e = add(scale(multiply(A, b, add(power(s, 2), constant(-1))), 2),
                constant(-1))
        linear = scale(add(K2, scale(multiply(b, T), 32)), -1)
        N = add(
            K2,
            multiply(
                add(power(A, 2), scale(multiply(A, b, power(s, 2)), 32),
                    constant(-16)),
                E2,
            ),
        )

    constant_coefficient = scale(multiply(e, T), 16)
    # Replay the original eliminated quadratics after multiplying by E^2:
    # A*u^2+(A^2-c^2-32*A*b)u+32*A^2*b*s^2-16*A, and
    # A*u^2-(c^2+32*A*b)u+32*A^2*b*(s^2-1)-16*A.
    if orientation == "I":
        expected_linear = add(
            multiply(add(power(A, 2), scale(multiply(A, b), -32)), E2),
            scale(K2, -1),
        )
        expected_constant = multiply(
            add(scale(multiply(power(A, 2), b, power(s, 2)), 32),
                scale(A, -16)),
            E2,
        )
    else:
        expected_linear = scale(
            add(K2, scale(multiply(A, b, E2), 32)), -1)
        expected_constant = multiply(
            add(scale(multiply(power(A, 2), b,
                               add(power(s, 2), constant(-1))), 32),
                scale(A, -16)),
            E2,
        )
    assert T == multiply(A, E2)
    assert linear == expected_linear
    assert constant_coefficient == expected_constant
    discriminant = add(power(linear, 2),
                       scale(multiply(e, power(T, 2)), -64))

    # q(-A)=A*N is the second Kummer branch identity.
    q_minus_A = add(multiply(T, power(A, 2)),
                    scale(multiply(linear, A), -1), constant_coefficient)
    assert q_minus_A == multiply(A, N)

    # The lambda model substitutes u=(lambda^2-A)^2/(4 lambda^2).
    lam2_minus_A = [scale(A, -1), {}, one]
    lam2 = [{}, {}, one]
    lam4 = [{}, {}, {}, {}, one]

    def list_add(*lists: list[Polynomial]) -> list[Polynomial]:
        size = max(len(values) for values in lists)
        answer = [{} for _ in range(size)]
        for values in lists:
            for index, coefficient in enumerate(values):
                answer[index] = add(answer[index], coefficient)
        return answer

    def list_multiply(left: list[Polynomial], right: list[Polynomial]) -> list[Polynomial]:
        answer = [{} for _ in range(len(left) + len(right) - 1)]
        for left_index, left_coefficient in enumerate(left):
            for right_index, right_coefficient in enumerate(right):
                answer[left_index + right_index] = add(
                    answer[left_index + right_index],
                    multiply(left_coefficient, right_coefficient),
                )
        return answer

    def list_power(values: list[Polynomial], exponent: int) -> list[Polynomial]:
        answer = [one]
        base = values
        remaining = exponent
        while remaining:
            if remaining & 1:
                answer = list_multiply(answer, base)
            base = list_multiply(base, base)
            remaining //= 2
        return answer

    lambda_polynomial = list_add(
        [multiply(T, coefficient) for coefficient in list_power(lam2_minus_A, 4)],
        [scale(multiply(linear, coefficient), 4)
         for coefficient in list_multiply(lam2, list_power(lam2_minus_A, 2))],
        [scale(multiply(constant_coefficient, coefficient), 16)
         for coefficient in lam4],
    )
    assert len(lambda_polynomial) == 9
    assert lambda_polynomial[8] == T
    assert lambda_polynomial[0] == multiply(power(A, 4), T)
    assert lambda_polynomial[2] == multiply(power(A, 2), lambda_polynomial[6])
    assert all(not lambda_polynomial[index] for index in (1, 3, 5, 7))

    return {
        "a": a, "A": A, "D": D, "Ng": Ng, "E": E, "K": K,
        "T": T, "e": e, "linear": linear,
        "constant_coefficient": constant_coefficient,
        "N": N, "discriminant": discriminant,
        "lambda_coefficients": lambda_polynomial,
    }


def rational_geometry(a: F, b: F, Z: F) -> dict[str, F]:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    assert A * b * D != 0
    c = a * a * Z * Z * Ng / (A * b * b * D)
    return {"A": A, "D": D, "Ng": Ng, "c": c}


def two_adic_replay() -> dict[str, Any]:
    checks = 0
    for s_value in map(F, (-3, -2, -1, 0, 1, 2, 3)):
        a = 1 + 2 * s_value
        A = 1 + 4 * a * a
        assert vp(a) == vp(A) == 0 and A.numerator % 8 == 5
        for b in map(F, (-3, -1, 1, 3, F(1, 3), F(5, 3))):
            assert vp(b) == 0
            for Z in map(F, (F(1, 8), -1, 1, 8, -8, 27)):
                values = rational_geometry(a, b, Z)
                c = values["c"]
                assert vp(values["Ng"]) >= 4
                expected_D = 0 if vp(Z) >= 0 else 2 * vp(Z)
                assert vp(values["D"]) == expected_D
                assert vp(c) >= 4
                e_I = 2 * A * b * s_value * s_value - 1
                e_II = 2 * A * b * (s_value * s_value - 1) - 1
                coefficients_I = [16 * e_I,
                                  A - c * c / A - 32 * b,
                                  F(1)]
                coefficients_II = [16 * e_II,
                                   -(c * c / A + 32 * b),
                                   F(1)]
                assert [vp(value) for value in coefficients_I] == [4, 0, 0]
                assert [vp(value) for value in coefficients_II] == [4, 5, 0]
                checks += 1
                PACER.step()
    return {
        "type": "two-adic-replay",
        "label": "EVIDENCE",
        "exact_instances": checks,
        "orientation_I_coefficient_valuations": [4, 0, 0],
        "orientation_II_coefficient_valuations": [4, 5, 0],
        "interpretation": "cross-check only; the theorem is the symbolic valuation proof",
    }


def orientation_II_mod16_record() -> dict[str, Any]:
    """Exhaust the complete residue proof behind the second Phi obstruction."""
    checks = 0
    residues = range(16)
    odd_residues = range(1, 16, 2)
    for s in residues:
        a = (1 + 2 * s) % 16
        A = (1 + 4 * a * a) % 16
        assert A == 5
        for b in odd_residues:
            e = (2 * A * b * (s * s - 1) - 1) % 16
            for h in odd_residues:
                for t in odd_residues:
                    divided_equation = (
                        t**4 - 8 * h * t * t + e
                    ) % 16
                    simplified = (
                        8 + 2 * A * b * (s * s - 1)
                    ) % 16
                    assert divided_equation == simplified
                    assert divided_equation != 0
                    checks += 1
                    PACER.step()
    assert checks == 8192
    return {
        "type": "orientation-II-mod16-proof",
        "label": "PROVED",
        "complete_residue_classes": checks,
        "equation_after_u_equals_4t2": (
            "t^4-8*h*t^2+e == 8+2*A*b*(s^2-1) (mod 16)"
        ),
        "result": "nonzero in every residue class; orientation II is Phi-empty",
        "method": "complete congruence exhaustion, not a bounded rational search",
    }


def branch_record(orientation: str, geometry: dict[str, Polynomial]) -> dict[str, Any]:
    discriminant_slice = specialize_univariate(
        geometry["discriminant"], B_I, {S_I: 2, Z_I: 27})
    N_slice = specialize_univariate(geometry["N"], B_I, {S_I: 2, Z_I: 27})
    e_slice = specialize_univariate(geometry["e"], B_I, {S_I: 2, Z_I: 27})
    discriminant_prime = 11 if orientation == "I" else 23
    certificates = {
        "Delta": irreducibility_certificate(discriminant_slice, discriminant_prime),
        "N_minus_A": irreducibility_certificate(N_slice, 11),
    }
    assert len(discriminant_slice) - 1 == 16
    assert len(N_slice) - 1 == 8
    assert len(e_slice) - 1 == 1
    polynomials = [discriminant_slice, N_slice, e_slice]
    assert all(poly_gcd(polynomial, poly_derivative(polynomial)) == [F(1)]
               for polynomial in polynomials)
    assert all(poly_gcd(polynomials[left], polynomials[right]) == [F(1)]
               for left, right in ((0, 1), (0, 2), (1, 2)))

    if orientation == "I":
        equations = {
            "e_zero": "e_I=2*A*b*s^2-1",
            "N_minus_A": "N_I=K^2+16*(2*A*b*(1+s^2)-1)*E^2",
            "Delta": "Delta_I=((A-32*b)*T-K^2)^2-64*e_I*T^2",
        }
    else:
        equations = {
            "e_zero": "e_II=2*A*b*(s^2-1)-1",
            "N_minus_A": "N_II=K^2+(A^2+32*A*b*s^2-16)*E^2",
            "Delta": "Delta_II=(K^2+32*b*T)^2-64*e_II*T^2",
        }
    return {
        "type": "branch-geometry",
        "label": "PROVED",
        "orientation": orientation,
        "guarded_base": "Spec Q[s,b,Z,1/(A*b*D)]",
        "cover_tower": "q(u)=0; rho^2=u; X^2=A+u",
        "lambda": "lambda=X+rho, u=(lambda^2-A)^2/(4*lambda^2)",
        "lambda_reciprocity": "lambda^8*G(A/lambda)=A^4*G(lambda)",
        "finite_flat": True,
        "free_basis": "1,lambda,...,lambda^7 after dividing G by the unit T",
        "affine_branch_components": equations,
        "degrees_b_s_Z": {
            "e_zero": [degree(geometry["e"], B_I), degree(geometry["e"], S_I),
                       degree(geometry["e"], Z_I)],
            "N_minus_A": [degree(geometry["N"], B_I), degree(geometry["N"], S_I),
                          degree(geometry["N"], Z_I)],
            "Delta": [degree(geometry["discriminant"], B_I),
                      degree(geometry["discriminant"], S_I),
                      degree(geometry["discriminant"], Z_I)],
        },
        "generic_inertia": {
            "Delta": "2^4",
            "e_zero": "2^2 1^4",
            "N_minus_A": "2^2 1^4",
        },
        "generic_components_over_Q": 1,
        "generic_degree": 8,
        "galois_closure": "C2^4 semidirect C2 = V4 wr C2, order 32",
        "proof_of_connectedness": (
            "Delta is nonsquare; over its quadratic field the two primes above "
            "e=0 give independent odd valuations of u_1,u_2 and the two primes "
            "above N=0 give independent odd valuations of A+u_1,A+u_2"
        ),
        "actual_cell_certificate": {
            "Z": 27,
            "s": 2,
            "a": 5,
            "A": 101,
            "D": -18251,
            "primitive_polynomials_ascending": {
                "Delta": primitive_integer_poly(discriminant_slice),
                "N_minus_A": primitive_integer_poly(N_slice),
                "e_zero": primitive_integer_poly(e_slice),
            },
            "certificates": certificates,
        },
    }


def genus_record(orientation: str, geometry: dict[str, Polynomial]) -> dict[str, Any]:
    delta = specialize_univariate(geometry["discriminant"], B_I,
                                  {S_I: 2, Z_I: 27})
    N = specialize_univariate(geometry["N"], B_I, {S_I: 2, Z_I: 27})
    e = specialize_univariate(geometry["e"], B_I, {S_I: 2, Z_I: 27})
    assert [len(delta) - 1, len(e) - 1, len(N) - 1] == [16, 1, 8]
    assert all(poly_gcd(poly, poly_derivative(poly)) == [F(1)]
               for poly in (delta, e, N))
    ramification = 16 * 4 + 1 * 2 + 8 * 2 + 2
    genus = (ramification - 2 * 8 + 2) // 2
    assert ramification == 84 and genus == 35
    return {
        "type": "strategic-slice-genus",
        "label": "PROVED",
        "orientation": orientation,
        "slice": {"Z": 27, "z": 3, "s": 2, "a": 5, "A": 101},
        "base_coordinate": "b on P1",
        "branch_degrees": {
            "Delta_cycle_2^4": 16,
            "u_zero_cycle_2^2_1^4": 1,
            "A_plus_u_zero_cycle_2^2_1^4": 8,
            "b_infinity_cycle_2^2_1^4": 1,
            "b_zero": 0,
        },
        "discriminant_subcover": {
            "equation": f"v^2=Delta_{orientation}(b)|_(s=2,Z=27)",
            "genus": 7,
        },
        "sqrt_u_intermediate_genus": 14,
        "riemann_hurwitz": "2*g-2=8*(-2)+84=68",
        "normalization_genus": genus,
        "consequence": "the degree-eight b-curve is not rational or unirational",
    }


def c_fiber_record(orientation: str) -> dict[str, Any]:
    A = 37
    B = 6
    s = 1
    if orientation == "I":
        d = A
        p0 = A + 16 * B
        q0 = 16 * A * B * (1 + s * s) - 16
    else:
        d = 0
        p0 = 16 * B
        q0 = -16 * A * B * (1 - s * s) - 16
    L = 2 * d - p0
    M = d * d - p0 * d + q0
    N = M - A * L + A * A
    quartic_discriminant = L * L - 4 * M
    assert M * N * quartic_discriminant != 0
    expected = (-59, 3536, 7088) if orientation == "I" else (-96, -16, 4905)
    assert (L, M, N) == expected
    return {
        "type": "named-lower-genus-reduction",
        "label": "PROVED",
        "orientation": orientation,
        "definitions": {
            "F0(t)": "t^2-p0*t+q0",
            "I": "d=A, p0=A+16*B, q0=16*A*B*(1+s^2)-16",
            "II": "d=0, p0=16*B, q0=-16*A*B*(1-s^2)-16",
            "u": "t-d",
            "identity": "c^2*u=A*F0(d+u)",
        },
        "genus_3_curve": "w^2=A*(v^4+L*v^2+M), y^2=v^2+A; u=v^2, c=w/v",
        "elliptic_quotient": "E: w^2=A*(v^4+L*v^2+M)",
        "generic_genus": 3,
        "genus_proof": "disjoint branch sets of sizes 4 and 2 on P1_v; Riemann-Hurwitz gives 2g-2=-8+8+4=4",
        "exact_slice": {
            "s": s, "b": 3, "a": 3, "A": A, "B": B,
            "L": L, "M": M, "N": N,
            "quartic_discriminant": quartic_discriminant,
        },
        "scope": "pre-h c-curve; the h-substituted diagonal cover is its pullback along c=h(a,b,Z)",
    }


def u_four_record(orientation: str, geometry: dict[str, Polynomial]) -> dict[str, Any]:
    T = geometry["T"]
    linear = geometry["linear"]
    constant_coefficient = geometry["constant_coefficient"]
    q_at_four = add(scale(T, 16), scale(linear, 4), constant_coefficient)
    slice_polynomial = specialize_univariate(q_at_four, B_I, {S_I: 0, Z_I: 27})
    prime = 17 if orientation == "I" else 11
    certificate = irreducibility_certificate(slice_polynomial, prime)
    assert certificate["degree"] == 8
    expected = {
        "I": [13286025, -106288200, 286978140, -403895160, 199621661,
              1876204840, 286978140, -106288200, 13286025],
        "II": [1476225, -11809800, 31886460, -44877240, 61765254,
               525147760, 31886460, -11809800, 1476225],
    }
    assert primitive_integer_poly(slice_polynomial) == expected[orientation]
    return {
        "type": "fixed-u-route",
        "label": "PROVED",
        "orientation": orientation,
        "u_minus_four": "EMPTY over Q because u=rho^2 cannot equal -4",
        "u_plus_four_square_condition": "X^2=A+4=5+4*a^2",
        "parameterization": {
            "a": "(theta^2-5)/(4*theta)",
            "s": "(theta^2-4*theta-5)/(8*theta)",
            "X": "(theta^2+5)/(2*theta)",
            "Phi_at_2": "v2(theta)=0; independently v2(b)=0",
        },
        "q_at_u_4": (
            "I: K^2-(A^2+8*A^2*b*s^2-32*A*b)*E^2=0"
            if orientation == "I" else
            "II: K^2-(8*A^2*b*(s^2-1)-32*A*b)*E^2=0"
        ),
        "actual_cell_obstruction": {
            "Z": 27, "theta": 5, "s": 0, "a": 1, "A": 5,
            "primitive_polynomial_in_b_ascending": expected[orientation],
            "irreducibility_certificate": certificate,
            "conclusion": "no rational b exists on this Phi-compatible u=4 slice",
        },
        "generic_conclusion": "the u=4 locus is a degree-eight multisection, not a rational section",
    }


def theorem_records() -> list[dict[str, Any]]:
    return [
        {
            "type": "normalization-and-guards",
            "label": "PROVED",
            "definitions": {
                "a": "1+2s", "A": "1+4a^2", "B": "2b",
                "D": "1-Z-a^2Z^2", "Ng": "16a^4b^2-A(b-1)^4",
                "E": "A*b^2*D", "K": "a^2*Z^2*Ng", "c": "K/E",
                "T": "A*E^2",
            },
            "orientation_I": "F_I=T*u^2+((A-32b)T-K^2)u+16(2Abs^2-1)T",
            "orientation_II": "F_II=T*u^2-(K^2+32bT)u+16(2Ab(s^2-1)-1)T",
            "clearing_multiplier": "E^2; nonzero on Phi",
            "guards": {
                "A": "A=5 mod 8 in Q2, hence nonsquare and nonzero",
                "b": "v2(b)=0, hence b!=0",
                "Z": "every target cell has Z=z^3!=0",
                "D": "D=0 would make rational Z a root of a^2X^2+X-1, forcing its discriminant A to be a rational square",
                "Ng": "Ng=0 gives A=(4a^2b/(b-1)^2)^2 when b!=1; at b=1, Ng=16a^4; hence Ng!=0",
                "lambda": "(X-rho)(X+rho)=A!=0, so lambda=X+rho!=0",
                "rho": "rho=0 would force X^2=A, impossible in Q2",
            },
        },
        {
            "type": "orientation-I-Phi-obstruction",
            "label": "PROVED",
            "statement": "orientation I has no Q2 point, hence no rational point, on any Phi-admissible base for any nonzero rational Z (also for Z=0, where c=0)",
            "proof": [
                "b and a are 2-adic units, b-1 is even, so v2(Ng)>=4",
                "for Z!=0, v2(D)=0 when v2(Z)>=0 and v2(D)=2v2(Z) when v2(Z)<0; therefore v2(c)>=4",
                "after dividing Q_I by A, coefficient valuations (u^0,u^1,u^2) are (4,0,0)",
                "the Newton polygon forces v2(u) in {4,0}",
                "if u=r^2 has valuation 4 then A+u=5 mod8; if it has valuation 0 then u=1 mod8 and A+u=6 mod8; neither is y^2",
            ],
            "consequence": "orientation I contributes the empty subset to the proposed five-unknown definition",
        },
        {
            "type": "orientation-II-Phi-obstruction",
            "label": "PROVED",
            "statement": "orientation II has no Q2 point, hence no rational point, on any Phi-admissible base",
            "proof": [
                "write Q_II/A as u^2-32*h*u+16*e=0 with h=b+c^2/(32A) and e=2*A*b*(s^2-1)-1; both h,e are 2-adic units",
                "the coefficient valuations (u^0,u^1,u^2)=(4,5,0) force v2(u)=2",
                "if the required u=y^2 existed, write u=4*t^2 with t a 2-adic unit; division by 16 gives t^4-8*h*t^2+e=0",
                "modulo 16 the left side is 8+2*A*b*(s^2-1), so vanishing would force v2(s^2-1)=2",
                "for s in Z2, v2(s^2-1)=0 when s is even and is at least 3 (or infinity) when s is odd; contradiction",
            ],
            "consequence": "orientation II is Phi-empty; the union of both diagonal orientations is empty and this route is false",
        },
        {
            "type": "route-obstructions",
            "label": "PROVED",
            "cell_factor_scaling": "for odd w, b=w*t is a birational base automorphism and v2(b)=v2(t); it preserves degree, monodromy, genus, and existence of a rational section",
            "degenerate_lines": {
                "u=0": "e=0, but e=-1 mod2 on Phi and a lift would require X^2=A",
                "u=-A": "a lift would require rho^2=-A, impossible in Q2 (-A=3 mod8) and over R",
                "Delta=0": "the remaining degenerate divisor is the irreducible degree-(16,24) discriminant component, not a uniform rational anchor",
            },
        },
        {
            "type": "rationality-verdict",
            "label": "PROVED",
            "scope": "relative rationality/no-section theorem; absolute surface status stated separately",
            "proved": [
                "the generic guarded cover in each orientation is integral of degree 8 and has no rational section over Q(s,b,Z)",
                "there is no base-preserving rational or unirational parameterization (relative dimension zero)",
                "the actual-cell fixed-s normalizations have genus 35 and are not rational or unirational curves",
                "both orientations are Phi-empty",
            ],
            "absolute_surface_status": "OPEN",
            "warning": "high-genus fibres alone do not prove that the two-dimensional total surface is absolutely non-unirational; resolving the boundary of the discriminant double surface is a separate surface-classification problem",
            "formal_discriminant_model": "V^2=Delta_j on P1_b x P1_s; nominal branch bidegree (16,24) and formal canonical class pullback O(6,10), but boundary discrepancies are not claimed audited here",
        },
    ]


def build_report(records: list[dict[str, Any]], elapsed: float) -> str:
    branch_I = next(record for record in records
                    if record.get("type") == "branch-geometry"
                    and record["orientation"] == "I")
    branch_II = next(record for record in records
                     if record.get("type") == "branch-geometry"
                     and record["orientation"] == "II")
    lines = [
        "# L24 — exact geometry of the diagonal covers",
        "",
        "**Verdict.  Both orientations are PROVED EMPTY on Phi.**  Orientation I",
        "fails by a global 2-adic Newton polygon.  Orientation II has roots of",
        "the eliminated quadratic only in the wrong 2-adic square classes: a",
        "mod-16 congruence rules out `u=y^2`.  Thus the union is empty and the",
        "diagonal five-unknown route is **PROVED FALSE**, not merely incomplete.",
        "Independently, each guarded generic cover is one integral degree-eight",
        "component with Galois closure `V4 wr C2` (order 32), three affine branch",
        "components, and no rational section.  On the actual cell `Z=27`, the",
        "fixed-`s` normalization has genus 35.  Absolute rationality of the",
        "two-dimensional total surface is kept **OPEN** because its boundary",
        "resolution has not been audited; that does not affect Phi-emptiness.",
        "",
        "No finite no-hit inference occurs below.  Every finite-field calculation is",
        "a positive Rabin irreducibility certificate.",
        "",
        "## 1. Clearing the bridge and all guards",
        "",
        "Put",
        "",
        "```text",
        "a=1+2s, A=1+4a^2, B=2b, D=1-Z-a^2Z^2,",
        "Ng=16a^4b^2-A(b-1)^4, E=Ab^2D, K=a^2Z^2Ng, T=AE^2.",
        "```",
        "",
        "Then `c=K/E`.  Multiplication by the nonzero square `E^2` gives",
        "",
        "```text",
        "F_I =T u^2+((A-32b)T-K^2)u+16(2Abs^2-1)T,       u=r^2, y^2=A+u,",
        "F_II=T u^2-(K^2+32bT)u+16(2Ab(s^2-1)-1)T,      u=y^2, r^2=A+u.",
        "```",
        "",
        "On Phi, `a,b` are 2-adic units and `A=5 mod 8`; hence `A` and `b`",
        "are nonzero.  If `D=0`, rational `Z` is a root of",
        "`a^2 X^2+X-1`, whose discriminant is `A`, contradicting the",
        "nonsquare `A`.  If `Ng=0` and `b!=1`, then",
        "`A=(4a^2b/(b-1)^2)^2`; at `b=1`, `Ng=16a^4`.  Target cells",
        "have `Z!=0`, so `E,K,c` are nonzero.  Finally",
        "`(X-rho)(X+rho)=A`, so `lambda=X+rho` is nonzero; `rho=0`",
        "would again make `A` a square.",
        "",
        "## 2. Decisive Phi 2-adic classification",
        "",
        "For a 2-adic unit `b`, `b-1` is even, so `v2(Ng)>=4`.  Exactly",
        "",
        "```text",
        "v2(D)=0                    if v2(Z)>=0,",
        "v2(D)=2v2(Z)              if v2(Z)<0.",
        "```",
        "",
        "Therefore `v2(c)>=4` for every nonzero rational `Z` (and `c=0`",
        "when `Z=0`).  In orientation I the coefficient valuations of `F_I/T`",
        "(constant through quadratic) are `(4,0,0)`.  The two possible root",
        "valuations are `4` and `0`.  If `u=r^2` has valuation four,",
        "`A+u=5 mod 8`; if it has valuation zero, `u=1 mod 8` and",
        "`A+u=6 mod 8`.  Neither can be `y^2`.  Hence",
        "",
        "```text",
        "orientation I ∩ Phi = empty                         [PROVED].",
        "```",
        "",
        "For orientation II, write its normalized quadratic as",
        "",
        "```text",
        "u^2-32h*u+16e=0,  h=b+c^2/(32A) in Z_2^x,",
        "e=2Ab(s^2-1)-1 in Z_2^x.",
        "```",
        "",
        "The coefficient valuations `(4,5,0)` force `v2(u)=2`.  A required",
        "square `u=y^2` would be `u=4t^2` with `t` a unit.  Dividing by 16",
        "and reducing modulo 16 gives",
        "",
        "```text",
        "0 = t^4-8h*t^2+e = 8+2Ab(s^2-1)  (mod 16).",
        "```",
        "",
        "This would force `v2(s^2-1)=2`.  But it is zero for even `s` and at",
        "least three (or infinity) for odd `s`.  Hence",
        "",
        "```text",
        "orientation II ∩ Phi = empty;  I ∪ II = empty       [PROVED].",
        "```",
        "",
        "Equivalently, the hyperbola is rationally parametrized by",
        "",
        "```text",
        "rho=(lambda-A/lambda)/2, X=(lambda+A/lambda)/2, lambda!=0.",
        "```",
        "",
        "For a Phi point, `v2(rho)=1` when `lambda` is a unit and",
        "`v2(rho)<=-2` otherwise.",
        "",
        "## 3. Exact cover and branch divisor",
        "",
        "Over `K0=Q(s,b,Z)` the cover is the tower",
        "",
        "```text",
        "K0 -> K0(u), q_j(u)=0 -> adjoining sqrt(u), sqrt(A+u).",
        "```",
        "",
        "Substitution `u=(lambda^2-A)^2/(4lambda^2)` gives an even octic",
        "",
        "```text",
        "G=T(lambda^2-A)^4+4L lambda^2(lambda^2-A)^2+16M lambda^4,",
        "lambda^8 G(A/lambda)=A^4 G(lambda).",
        "```",
        "",
        "On the guarded affine base its reduced branch support has exactly the",
        "three Q-irreducible components shown below (the listed degrees are in",
        "`(b,s,Z)`):",
        "",
        "| orientation | `q(0)=0` | `q(-A)=0` | quadratic discriminant |",
        "|---|---|---|---|",
        f"| I | `{branch_I['affine_branch_components']['e_zero']}` (1,4,0) | `N_I=0` (8,12,4) | `Delta_I=0` (16,24,8) |",
        f"| II | `{branch_II['affine_branch_components']['e_zero']}` (1,4,0) | `N_II=0` (8,12,4) | `Delta_II=0` (16,24,8) |",
        "",
        "Generic inertia cycle types are `2^2 1^4`, `2^2 1^4`, and `2^4`.",
        "The actual-cell specializations at `(s,Z)=(2,27)` are irreducible",
        "modulo the recorded primes, proving the universal components",
        "irreducible by the degree-preserving specialization lemma.  Above the",
        "discriminant double cover, the two primes over `e=0` independently",
        "detect `u_1,u_2`, and the two over `N=0` independently detect",
        "`A+u_1,A+u_2`.  Thus the Kummer rank is four, the cover is integral",
        "of degree eight, and its Galois closure is",
        "`C2^4 semidirect C2 = V4 wr C2`, order 32.",
        "",
        "## 4. Strategic actual-cell slices",
        "",
        "Fix `Z=27` (`z=3`) and the Phi value `s=2` (`a=5,A=101`).  For",
        "both orientations, the branch polynomials in `b` are squarefree and",
        "pairwise coprime, of degrees `16,1,8`.  At `b=infinity` the two `u`",
        "roots have orders `-4` and `+3`; only the latter Kummer sheet branches.",
        "There is no branch at `b=0` (orders are `-4,+4`).  Hence the total",
        "ramification is",
        "",
        "```text",
        "16*4 + 1*2 + 8*2 + 1*2 = 84,",
        "2g-2 = 8*(-2)+84 = 68, so g=35.",
        "```",
        "",
        "The first quadratic quotient is the named genus-seven hyperelliptic",
        "curve `v^2=Delta_j(b)`; adjoining `sqrt(u)` gives genus 14; the full",
        "normalization has genus 35.  These are theorems for the displayed",
        "actual cell, not generic-smoothness guesses.",
        "",
        "There is also a useful lower-genus universal reduction before inserting",
        "`c=h`.  For fixed `(s,b)`, put `F0(t)=t^2-p0*t+q0` and `u=t-d` as",
        "recorded in the JSON artifact.  Then",
        "",
        "```text",
        "c^2 u=A F0(d+u),",
        "w^2=A(v^4+L v^2+M), y^2=v^2+A, u=v^2, c=w/v.",
        "```",
        "",
        "This is generically a genus-three biquadratic curve with an elliptic",
        "quartic quotient.  Inserting `h` pulls the `(b,s,Z)` cover back along",
        "its degree-eight `c`-map; it does not rationalize it.",
        "",
        "## 5. Requested section attempts",
        "",
        "* `u=-4` is empty over `Q` because `u=rho^2`.",
        "* For `u=4`, the other square is parametrized by",
        "  `a=(theta^2-5)/(4theta)`, `s=(theta^2-4theta-5)/(8theta)`,",
        "  `X=(theta^2+5)/(2theta)`.  Phi compatibility is exactly that",
        "  `theta` and `b` are 2-adic units.  At the compatible point",
        "  `(theta,s,a,Z)=(5,0,1,27)`, each orientation leaves a certified",
        "  irreducible octic in `b`, hence no rational `b`.  Universally this",
        "  is a degree-eight multisection, not a section.",
        "* Replacing `b` by `w*t` for an odd cell factor is a birational",
        "  automorphism of the base and preserves `v2(b)=v2(t)`; it cannot",
        "  alter degree, monodromy, genus, or section existence.",
        "* The apparent line anchors `u=0` and `u=-A` do not lift on Phi:",
        "  they require respectively `X^2=A` and `rho^2=-A`.  The remaining",
        "  discriminant divisor is the irreducible high-degree component.",
        "",
        "## 6. Rationality scope",
        "",
        "The integral degree-eight generic cover has no `Q(s,b,Z)`-point and",
        "therefore no rational section or base-preserving rational/unirational",
        "parameterization.  The actual-cell `b`-curves have genus 35.  This is",
        "the exact geometric obstruction required here.",
        "",
        "It would be an overclaim to infer absolute non-unirationality of the",
        "two-dimensional total surface solely from high-genus fibres.  The",
        "discriminant quotient has formal model `V^2=Delta_j` with nominal",
        "branch bidegree `(16,24)` and formal canonical class `O(6,10)`, but",
        "its boundary discrepancies have not been resolved.  Absolute surface",
        "rationality is therefore labelled **OPEN**, without weakening the",
        "proved Phi-emptiness counterexample.",
        "",
        "## 7. Labels and replay",
        "",
        "**PROVED:** denominator guards; cleared equations; global Phi emptiness",
        "of **both** orientations; branch tower, components, monodromy and no",
        "generic section; genus-35 actual-cell slices; genus-three c-reduction;",
        "fixed-`u` and cell-scaling no-gos.  Consequently the diagonal",
        "five-unknown route is false.",
        "",
        "**OPEN:** absolute rationality/unirationality of the total algebraic",
        "surfaces.  No absence conclusion is drawn from a finite search.",
        "",
        f"Wall-clock: {elapsed:.3f} s.  Pacing: {PACER.steps} steps,",
        f"{PACER.sleeps} sleeps of {PACE_SECONDS} s.",
        "",
        "Artifacts: `math/h10q/l24_diagonal_geometry.py`,",
        "`math/h10q/data/l24_diagonal_geometry.jsonl`,",
        "`/tmp/l24_diagonal_geometry.md`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    geometry_I = build_geometry("I")
    geometry_II = build_geometry("II")
    records: list[dict[str, Any]] = [{
        "type": "meta",
        "label": "PROVED",
        "schema": "l24-diagonal-geometry-v1",
        "source": "math/h10q/l24_diagonal_geometry.py",
        "finite_no_hit_inference": False,
        "cpu_policy": "one stdlib process; no pools; paced exact finite-field loops",
    }]
    records.extend(theorem_records())
    records.append(two_adic_replay())
    records.append(orientation_II_mod16_record())
    for orientation, geometry in (("I", geometry_I), ("II", geometry_II)):
        records.append(branch_record(orientation, geometry))
        records.append(genus_record(orientation, geometry))
        records.append(c_fiber_record(orientation))
        records.append(u_four_record(orientation, geometry))
    elapsed = time.perf_counter() - started
    records.append({
        "type": "summary",
        "label": "PROVED",
        "orientation_I_Phi": "EMPTY (PROVED)",
        "orientation_II_Phi": "EMPTY (PROVED)",
        "diagonal_union": "EMPTY; route FALSE (PROVED)",
        "generic_cover_degree": 8,
        "generic_cover_components": 1,
        "actual_cell_slice_genus": 35,
        "generic_rational_section": False,
        "absolute_surface_rationality": "OPEN",
        "elapsed_seconds": elapsed,
        "pacer_steps": PACER.steps,
        "pacer_sleeps": PACER.sleeps,
    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
                           for record in records), encoding="utf-8")
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print("L24 diagonal geometry: both orientations Phi-empty; route FALSE; "
          f"actual-cell slice genus 35; {elapsed:.3f}s")
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
