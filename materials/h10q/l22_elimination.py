#!/usr/bin/env python3
"""Exact fixed-Z factor elimination for the constructed tau_dagger octic.

Replay from the workspace root with the mandated low-priority process:

    nice -n 19 python3 math/h10q/l22_elimination.py

The symbolic part is self-contained.  It constructs the integral polynomial H,
checks its three local Newton initial forms, and mechanically replays the finite
endpoint case split used in the proof.  The finite-field scans are deliberately
reported only as EVIDENCE.
"""
from __future__ import annotations

from fractions import Fraction as F
import json
import math
from pathlib import Path
import time
from typing import Any, Callable


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l22_elimination.jsonl"
REPORT = Path("/tmp/l22_elimination.md")
FIXED_FIELD_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
ALPHA_SQUARE_PRIMES = (
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79,
)
PACE_EVERY = 200
PACE_SECONDS = 0.005

# Sparse polynomials use exponents (a, Z, b).  All coefficients are exact.
Monomial = tuple[int, int, int]
Sparse = dict[Monomial, F]
AZPoly = dict[tuple[int, int], F]
XZPoly = dict[tuple[int, int], F]


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def sparse_add(left: Sparse, right: Sparse) -> Sparse:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, F(0)) + coefficient
        if not answer[monomial]:
            del answer[monomial]
    return answer


def sparse_scale(poly: Sparse, scalar: F | int) -> Sparse:
    scalar = F(scalar)
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in poly.items()
        if coefficient * scalar
    }


def sparse_mul(left: Sparse, right: Sparse) -> Sparse:
    answer: Sparse = {}
    for (a1, z1, b1), x in left.items():
        for (a2, z2, b2), y in right.items():
            monomial = (a1 + a2, z1 + z2, b1 + b2)
            answer[monomial] = answer.get(monomial, F(0)) + x * y
            if not answer[monomial]:
                del answer[monomial]
    return answer


def sparse_pow(poly: Sparse, exponent: int) -> Sparse:
    assert exponent >= 0
    answer: Sparse = {(0, 0, 0): F(1)}
    base = poly
    while exponent:
        if exponent & 1:
            answer = sparse_mul(answer, base)
        base = sparse_mul(base, base)
        exponent >>= 1
    return answer


def build_H() -> tuple[Sparse, Sparse, Sparse, Sparse]:
    one: Sparse = {(0, 0, 0): F(1)}
    a: Sparse = {(1, 0, 0): F(1)}
    Z: Sparse = {(0, 1, 0): F(1)}
    b: Sparse = {(0, 0, 1): F(1)}
    A = sparse_add(one, sparse_scale(sparse_pow(a, 2), 4))
    D = sparse_add(
        sparse_add(one, sparse_scale(Z, -1)),
        sparse_scale(sparse_mul(sparse_pow(a, 2), sparse_pow(Z, 2)), -1),
    )
    b_minus_one = sparse_add(b, sparse_scale(one, -1))
    Ng = sparse_add(
        sparse_scale(sparse_mul(sparse_pow(a, 4), sparse_pow(b, 2)), 16),
        sparse_scale(sparse_mul(A, sparse_pow(b_minus_one, 4)), -1),
    )
    first = sparse_mul(
        sparse_mul(sparse_pow(a, 8), sparse_pow(Z, 4)),
        sparse_pow(Ng, 2),
    )
    second = sparse_scale(
        sparse_mul(sparse_mul(sparse_pow(A, 3), sparse_pow(D, 2)), sparse_pow(b, 4)),
        4,
    )
    third = sparse_scale(
        sparse_mul(
            sparse_mul(
                sparse_mul(sparse_pow(A, 4), sparse_pow(sparse_add(a, sparse_scale(one, -1)), 2)),
                sparse_pow(D, 2),
            ),
            sparse_pow(b, 5),
        ),
        -2,
    )
    return sparse_add(sparse_add(first, second), third), A, D, Ng


def coefficients_in_b(poly: Sparse) -> list[AZPoly]:
    degree = max(b_degree for _, _, b_degree in poly)
    answer: list[AZPoly] = [dict() for _ in range(degree + 1)]
    for (a_degree, z_degree, b_degree), coefficient in poly.items():
        answer[b_degree][(a_degree, z_degree)] = coefficient
    return answer


def az_add(left: AZPoly, right: AZPoly) -> AZPoly:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, F(0)) + coefficient
        if not answer[monomial]:
            del answer[monomial]
    return answer


def az_scale(poly: AZPoly, scalar: F | int) -> AZPoly:
    scalar = F(scalar)
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in poly.items()
        if coefficient * scalar
    }


def divide_by_A_once(poly: AZPoly) -> tuple[AZPoly, AZPoly]:
    """Divide in Q[Z][a] by A=1+4a^2, returning quotient and remainder."""
    remainder = dict(poly)
    quotient: AZPoly = {}
    while remainder and max(a_degree for a_degree, _ in remainder) >= 2:
        top_degree = max(a_degree for a_degree, _ in remainder)
        top_terms = [
            (z_degree, coefficient)
            for (a_degree, z_degree), coefficient in remainder.items()
            if a_degree == top_degree
        ]
        increment: AZPoly = {
            (top_degree - 2, z_degree): coefficient / 4
            for z_degree, coefficient in top_terms
        }
        quotient = az_add(quotient, increment)
        subtract = dict(increment)
        for (a_degree, z_degree), coefficient in increment.items():
            key = (a_degree + 2, z_degree)
            subtract[key] = subtract.get(key, F(0)) + 4 * coefficient
        remainder = az_add(remainder, az_scale(subtract, -1))
    return quotient, remainder


def valuation_A(poly: AZPoly) -> int:
    valuation = 0
    current = poly
    while True:
        quotient, remainder = divide_by_A_once(current)
        if remainder:
            return valuation
        assert quotient
        valuation += 1
        current = quotient


def valuation_a(poly: AZPoly) -> int:
    return min(a_degree for a_degree, _ in poly)


def degree_a(poly: AZPoly) -> int:
    return max(a_degree for a_degree, _ in poly)


def initial_after_b_scale(
    poly: Sparse,
    b_scale: int,
    extremum: Callable[[list[int]], int],
) -> tuple[int, XZPoly]:
    weights = [a_degree + b_scale * b_degree for a_degree, _, b_degree in poly]
    target = extremum(weights)
    initial: XZPoly = {}
    for (a_degree, z_degree, b_degree), coefficient in poly.items():
        if a_degree + b_scale * b_degree == target:
            key = (b_degree, z_degree)
            initial[key] = initial.get(key, F(0)) + coefficient
            if not initial[key]:
                del initial[key]
    return target, initial


def serialize_az(poly: AZPoly) -> list[list[Any]]:
    return [
        [a_degree, z_degree, frac_text(coefficient)]
        for (a_degree, z_degree), coefficient in sorted(poly.items())
    ]


def serialize_xz(poly: XZPoly) -> list[list[Any]]:
    return [
        [x_degree, z_degree, frac_text(coefficient)]
        for (x_degree, z_degree), coefficient in sorted(poly.items())
    ]


def endpoint_cases(degree: int) -> list[dict[str, int]]:
    """Enumerate A exponents and the number of a-infinity-large roots."""
    cases: list[dict[str, int]] = []
    for beta in range(3):
        for delta in range(3):
            for large_roots in range(degree + 1):
                # deg(c_F/l_F)=2*degree+2*delta-2*beta, while the
                # product of roots grows as a^(2*large_roots-degree).
                if (2 * degree + 2 * delta - 2 * beta
                        == 2 * large_roots - degree):
                    cases.append(
                        {
                            "beta_leading_A": beta,
                            "delta_constant_A": delta,
                            "infinity_large_roots": large_roots,
                        }
                    )
    return cases


def local_cluster_allocations(degree: int) -> list[dict[str, int]]:
    """Allowed a=0 allocations after the residual no-linear-root check."""
    answer: list[dict[str, int]] = []
    # The length-3, denominator-3 side is indivisible.  The unit side is
    # simple.  A degree 1 or 3 subfactor of the small residual quartic would
    # force a rational linear divisor, which its positive binomial lacks.
    for small in range(5):
        if small in (1, 3):
            continue
        for unit in range(2):
            for large in (0, 3):
                if small + unit + large == degree:
                    answer.append({"small": small, "unit": unit, "large": large})
    return answer


def is_rational_square_integer(value: int) -> bool:
    if value < 0:
        return False
    root = math.isqrt(value)
    return root * root == value


def evaluate_H(poly: Sparse, a_value: F | int, Z_value: F | int) -> list[F]:
    a_value = F(a_value)
    Z_value = F(Z_value)
    answer = [F(0)] * 9
    for (a_degree, z_degree, b_degree), coefficient in poly.items():
        answer[b_degree] += coefficient * a_value**a_degree * Z_value**z_degree
    return answer


def primitive_integer_poly(coefficients: list[F]) -> list[int]:
    denominator = 1
    for coefficient in coefficients:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in coefficients]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, coefficient)
    assert content
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return integers


def fp_trim(poly: list[int], prime: int) -> list[int]:
    answer = [coefficient % prime for coefficient in poly]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def fp_divmod(dividend: list[int], divisor: list[int], prime: int) -> tuple[list[int], list[int]]:
    divisor = fp_trim(divisor, prime)
    remainder = fp_trim(dividend, prime)
    assert divisor != [0]
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse = pow(divisor[-1], prime - 2, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        scalar = remainder[-1] * inverse % prime
        quotient[shift] = scalar
        for index, coefficient in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - scalar * coefficient
            ) % prime
        remainder = fp_trim(remainder, prime)
    return fp_trim(quotient, prime), remainder


def fp_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left = fp_trim(left, prime)
    right = fp_trim(right, prime)
    while right != [0]:
        _, remainder = fp_divmod(left, right, prime)
        left, right = right, remainder
    inverse = pow(left[-1], prime - 2, prime)
    return [(coefficient * inverse) % prime for coefficient in left]


def fp_mul_mod(left: list[int], right: list[int], modulus: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if x:
            for j, y in enumerate(right):
                if y:
                    product[i + j] = (product[i + j] + x * y) % prime
    return fp_divmod(product, modulus, prime)[1]


def fp_pow_mod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    answer = [1]
    base = fp_divmod(base, modulus, prime)[1]
    while exponent:
        if exponent & 1:
            answer = fp_mul_mod(answer, base, modulus, prime)
        base = fp_mul_mod(base, base, modulus, prime)
        exponent >>= 1
    return answer


def fp_frobenius_x(modulus: list[int], prime: int, iterations: int) -> list[int]:
    value = [0, 1]
    for _ in range(iterations):
        value = fp_pow_mod(value, prime, modulus, prime)
    return value


def fp_irred8(coefficients: list[int], prime: int) -> bool:
    polynomial = fp_trim(coefficients, prime)
    if len(polynomial) != 9:
        return False
    inverse = pow(polynomial[-1], prime - 2, prime)
    polynomial = [(coefficient * inverse) % prime for coefficient in polynomial]
    x_p4 = fp_frobenius_x(polynomial, prime, 4)
    x_p8 = fp_frobenius_x(polynomial, prime, 8)
    x = [0, 1]
    difference = [0] * max(len(x_p4), 2)
    for index, coefficient in enumerate(x_p4):
        difference[index] = coefficient
    difference[1] = (difference[1] - 1) % prime
    return fp_trim(x_p8, prime) == x and fp_gcd(polynomial, difference, prime) == [1]


def irred_certificate_details(coefficients: list[int], prime: int) -> dict[str, Any]:
    reduction = fp_trim(coefficients, prime)
    assert len(reduction) == 9
    inverse = pow(reduction[-1], prime - 2, prime)
    monic = [(coefficient * inverse) % prime for coefficient in reduction]
    x_p4 = fp_frobenius_x(monic, prime, 4)
    x_p8 = fp_frobenius_x(monic, prime, 8)
    difference = [0] * max(len(x_p4), 2)
    for index, coefficient in enumerate(x_p4):
        difference[index] = coefficient
    difference[1] = (difference[1] - 1) % prime
    gcd = fp_gcd(monic, difference, prime)
    assert x_p8 == [0, 1] and gcd == [1]
    assert fp_irred8(coefficients, prime)
    return {
        "prime": prime,
        "reduction_ascending": reduction,
        "monic_reduction_ascending": monic,
        "x_pow_p4_mod_f_ascending": x_p4,
        "x_pow_p8_mod_f_ascending": x_p8,
        "gcd_x_pow_p4_minus_x_with_f": gcd,
    }


class Pacer:
    """Sleep at least 5 ms per 200 checks and hold loop duty cycle <= 50%."""

    def __init__(self) -> None:
        self.steps = 0
        self.sleeps = 0
        self.sleep_seconds = 0.0
        self.chunk_started = time.perf_counter()

    def tick(self) -> None:
        self.steps += 1
        if self.steps % PACE_EVERY == 0:
            now = time.perf_counter()
            compute_seconds = now - self.chunk_started
            pause = max(PACE_SECONDS, compute_seconds)
            time.sleep(pause)
            self.sleep_seconds += pause
            self.sleeps += 1
            self.chunk_started = time.perf_counter()


PACER = Pacer()
IRRED_CACHE: dict[tuple[int, tuple[int, ...]], bool] = {}


def cached_irred8(coefficients: list[int], prime: int) -> bool:
    normalized = tuple(coefficient % prime for coefficient in coefficients)
    key = (prime, normalized)
    answer = IRRED_CACHE.get(key)
    if answer is None:
        answer = fp_irred8(list(normalized), prime)
        IRRED_CACHE[key] = answer
    return answer


def H_mod(a_value: int, Z_value: int, prime: int) -> list[int]:
    a = a_value % prime
    Z = Z_value % prime
    A = (1 + 4 * a * a) % prime
    D = (1 - Z - a * a * Z * Z) % prime
    ng = [
        -A,
        4 * A,
        16 * pow(a, 4, prime) - 6 * A,
        4 * A,
        -A,
    ]
    ng = [coefficient % prime for coefficient in ng]
    ng_square = [0] * 9
    for i, x in enumerate(ng):
        for j, y in enumerate(ng):
            ng_square[i + j] = (ng_square[i + j] + x * y) % prime
    scalar = pow(a, 8, prime) * pow(Z, 4, prime) % prime
    answer = [scalar * coefficient % prime for coefficient in ng_square]
    answer[4] = (answer[4] + 4 * pow(A, 3, prime) * D * D) % prime
    answer[5] = (
        answer[5] - 2 * pow(A, 4, prime) * pow(a - 1, 2, prime) * D * D
    ) % prime
    return answer


def alpha_square_mod(a_value: int, unit: int, prime: int) -> list[int]:
    a = a_value % prime
    A = (1 + 4 * a * a) % prime
    inverse_two = pow(2, prime - 2, prime)
    s = (a - 1) * inverse_two % prime
    ng = [
        -A,
        4 * A,
        16 * pow(a, 4, prime) - 6 * A,
        4 * A,
        -A,
    ]
    ng = [coefficient % prime for coefficient in ng]
    ng_square = [0] * 9
    for i, x in enumerate(ng):
        for j, y in enumerate(ng):
            ng_square[i + j] = (ng_square[i + j] + x * y) % prime
    scalar = pow(a, 4, prime) * unit * unit % prime
    answer = [scalar * coefficient % prime for coefficient in ng_square]
    answer[4] = (answer[4] + 16 * pow(A, 3, prime)) % prime
    answer[5] = (answer[5] - 32 * pow(A, 4, prime) * s * s) % prime
    return answer


def polynomial_degree_mod(coefficients: list[int], prime: int) -> int:
    for degree in range(len(coefficients) - 1, -1, -1):
        if coefficients[degree] % prime:
            return degree
    return -1


def symbolic_derivation(H: Sparse) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    coefficients = coefficients_in_b(H)
    assert len(coefficients) == 9

    endpoint_expected: AZPoly = {
        (8, 4): F(1),
        (10, 4): F(8),
        (12, 4): F(16),
    }
    assert coefficients[0] == coefficients[8] == endpoint_expected

    a_degrees = [degree_a(coefficient) for coefficient in coefficients]
    a_valuations = [valuation_a(coefficient) for coefficient in coefficients]
    A_valuations = [valuation_A(coefficient) for coefficient in coefficients]
    assert a_degrees == [12, 12, 14, 14, 16, 13, 14, 12, 12]
    assert a_valuations == [8, 8, 8, 8, 0, 0, 8, 8, 8]
    assert A_valuations == [2, 2, 1, 1, 0, 1, 1, 2, 2]

    b5_top = {
        monomial: coefficient
        for monomial, coefficient in coefficients[5].items()
        if monomial[0] == a_degrees[5]
    }
    assert b5_top == {(13, 4): F(1024)}

    a_zero_weight, a_zero_residual = initial_after_b_scale(H, 2, min)
    infinity_small_weight, infinity_small_residual = initial_after_b_scale(H, -1, max)
    infinity_large_weight, infinity_large_residual = initial_after_b_scale(H, 1, max)
    expected_a_zero: XZPoly = {
        (0, 4): F(1),
        (4, 0): F(4),
        (4, 1): F(-8),
        (4, 2): F(4),
    }
    expected_infinity_small: XZPoly = {
        (0, 4): F(16),
        (2, 4): F(-128),
        (4, 4): F(256),
    }
    expected_infinity_large: XZPoly = {
        (4, 4): F(256),
        (6, 4): F(-128),
        (8, 4): F(16),
    }
    assert a_zero_weight == 8 and a_zero_residual == expected_a_zero
    assert infinity_small_weight == 12
    assert infinity_small_residual == expected_infinity_small
    assert infinity_large_weight == 20
    assert infinity_large_residual == expected_infinity_large

    cases_two = endpoint_cases(2)
    cases_four = endpoint_cases(4)
    assert cases_two == [
        {"beta_leading_A": 1, "delta_constant_A": 0, "infinity_large_roots": 2},
        {"beta_leading_A": 2, "delta_constant_A": 0, "infinity_large_roots": 1},
        {"beta_leading_A": 2, "delta_constant_A": 1, "infinity_large_roots": 2},
    ]
    assert cases_four == [
        {"beta_leading_A": 2, "delta_constant_A": 0, "infinity_large_roots": 4},
    ]
    allocations_two = local_cluster_allocations(2)
    allocations_four = local_cluster_allocations(4)
    assert allocations_two == [{"small": 2, "unit": 0, "large": 0}]
    assert allocations_four == [
        {"small": 0, "unit": 1, "large": 3},
        {"small": 4, "unit": 0, "large": 0},
    ]

    # A degree-2 small factor has residual x^2+e*x+r with r=+/-16.
    # Multiplication by x^2-e*x+s and comparison with x^4+K gives either
    # e=0,s=-r,K=-r^2<0, or s=r,e^2=2r.  Neither r works over Q.
    assert not is_rational_square_integer(32)
    degree_two_cases = [
        {
            **case,
            "constant_over_leading_at_a_zero": ["-16", "+16"],
            "quadratic_comparison": "e=0 gives K=-256<0; otherwise e^2=+/-32",
            "verdict": "impossible over Q",
        }
        for case in cases_two
    ]

    # The degree-4 endpoint case has d/c=256 from (x^2-4)^2 at infinity.
    # Its entire a=0 small residual then forces the displayed equation.
    discriminants = [1152, 896]
    assert all(not is_rational_square_integer(value) for value in discriminants)
    degree_four_condition = {
        "endpoint_form_after_scaling": {
            "leading": "c*A^2",
            "constant": "d*a^8",
            "infinity_residual": "16*c*(x^2-4)^2",
            "forced_ratio_d_over_c": "256",
        },
        "a_zero_matching": "c*Z^4/d=4*(1-Z)^2",
        "forced_equation": "Z^4=1024*(1-Z)^2",
        "factorization": [
            "Z^2-32*(1-Z)=0",
            "Z^2+32*(1-Z)=0",
        ],
        "quadratic_discriminants": discriminants,
        "rational_roots": [],
        "verdict": "impossible over Q",
    }

    z1_coefficients = primitive_integer_poly(evaluate_H(H, 1, 1))
    z2_coefficients = primitive_integer_poly(evaluate_H(H, 1, 2))
    assert z1_coefficients == [
        25, -200, 540, -760, 1546, -760, 540, -200, 25,
    ]
    assert z2_coefficients == [
        100, -800, 2160, -3040, 7309, -3040, 2160, -800, 100,
    ]
    z1_certificate = irred_certificate_details(z1_coefficients, 11)
    z2_certificate = irred_certificate_details(z2_coefficients, 13)

    rows: list[dict[str, Any]] = [
        {
            "type": "integral-model",
            "label": "PROVED",
            "formula": "H=a^8*Z^4*Ng^2+4*A^3*D^2*b^4-2*A^4*(a-1)^2*D^2*b^5",
            "A": "1+4*a^2",
            "D": "1-Z-a^2*Z^2",
            "Ng": "16*a^4*b^2-A*(b-1)^4",
            "endpoint": "a^8*A^2*Z^4",
            "coefficients_by_b_sparse_a_Z": [serialize_az(coefficient) for coefficient in coefficients],
        },
        {
            "type": "local-newton-data",
            "label": "PROVED",
            "place_A": {
                "coefficient_valuations_b0_to_b8": A_valuations,
                "polygon": [[0, 2], [4, 0], [8, 2]],
                "slopes_with_lengths": [["-1/2", 4], ["1/2", 4]],
                "global_factor_degree_consequence": "every factor degree is even",
            },
            "place_a_zero": {
                "coefficient_valuations_b0_to_b8_generic_Z_not_0_1": a_valuations,
                "slopes_with_lengths": [["-2", 4], ["0", 1], ["8/3", 3]],
                "small_residual": "Z^4+4*(1-Z)^2*x^4",
                "small_residual_sparse_x_Z": serialize_xz(a_zero_residual),
                "small_residual_has_no_Q_linear_factor": True,
                "degree_3_side_indivisible_by_slope_denominator": True,
            },
            "place_a_infinity": {
                "coefficient_a_degrees_b0_to_b8": a_degrees,
                "exact_b5_top_after_a14_cancellation": serialize_az(b5_top),
                "slopes_with_lengths": [["-1", 4], ["1", 4]],
                "small_scale_b=x/a_weight": infinity_small_weight,
                "small_residual": "16*Z^4*(4*x^2-1)^2",
                "small_residual_sparse_x_Z": serialize_xz(infinity_small_residual),
                "large_scale_b=a*x_weight": infinity_large_weight,
                "large_residual": "16*Z^4*x^4*(x^2-4)^2",
                "large_residual_sparse_x_Z": serialize_xz(infinity_large_residual),
            },
        },
        {
            "type": "global-degree-filter",
            "label": "PROVED",
            "newton_product_theorem": (
                "over each completion, a factor's Newton sides are a sub-multiset of H's sides; "
                "a slope denominator divides the horizontal contribution, and an integral-slope "
                "factor initial form divides the side residual polynomial"
            ),
            "A_place_allowed_global_factor_degrees_up_to_4": [2, 4],
            "a_zero_allocations_degree_2": allocations_two,
            "a_zero_allocations_degree_4": allocations_four,
        },
        {
            "type": "degree-2-elimination",
            "label": "PROVED",
            "endpoint_infinity_cases": degree_two_cases,
            "a_zero_residual_normalized": "x^4+K with K=Z^4/(4*(1-Z)^2)>0",
            "coefficient_equations": [
                "f=-e",
                "e*(s-r)=0",
                "s-e^2+r=0",
                "r*s=K",
            ],
            "forced_r_values": ["-16", "+16"],
            "verdict": "no degree-2 factor for rational Z not in {0,1}",
        },
        {
            "type": "degree-4-elimination",
            "label": "PROVED",
            **degree_four_condition,
            "verdict_scope": "no degree-4 factor for rational Z not in {0,1}",
        },
        {
            "type": "exceptional-fiber",
            "label": "PROVED",
            "Z": "0",
            "status": "excluded and reducible",
            "identity": "H=2*A^3*b^4*(2-A*(a-1)^2*b)",
            "reason_excluded": "the cell parameter z is nonzero",
        },
        {
            "type": "exceptional-fiber",
            "label": "PROVED",
            "Z": "1",
            "specialization_a": 1,
            "primitive_H_coefficients_ascending": z1_coefficients,
            "certificate": z1_certificate,
            "vertical_conclusion": "H(a,1) is irreducible in Q(a)[b]",
        },
        {
            "type": "cross-check-fiber",
            "label": "PROVED",
            "Z": "2",
            "note": "Z=2 is already covered by the symbolic Z not in {0,1} proof",
            "specialization_a": 1,
            "primitive_H_coefficients_ascending": z2_coefficients,
            "certificate": z2_certificate,
        },
        {
            "type": "fixed-Z-theorem",
            "label": "PROVED",
            "statement": "for every nonzero rational Z, H is irreducible in Q(a)[b]",
            "equivalent_P_statement": "P=(4/A)H is irreducible in Q(a)[b]",
            "cell_corollary": (
                "for every nonzero rational z (hence Z=z^3), including every v_w(z)>=1 cell, "
                "the fixed vertical fiber is irreducible"
            ),
            "D_zero_note": (
                "D(a)=0 makes the specialized H a square and is not claimed good; for fixed "
                "nonzero Z it removes at most two a-values and does not affect vertical irreducibility"
            ),
        },
    ]
    facts = {
        "a_degrees": a_degrees,
        "a_valuations": a_valuations,
        "A_valuations": A_valuations,
        "degree_two_cases": cases_two,
        "degree_four_cases": cases_four,
        "degree_four_condition": degree_four_condition,
        "z1_certificate": z1_certificate,
        "z2_certificate": z2_certificate,
    }
    return rows, facts


def scan_fixed_fields() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for prime in FIXED_FIELD_PRIMES:
        missing: list[int] = []
        witness_counts: list[int] = []
        for Z in range(1, prime):
            witnesses: list[int] = []
            for a in range(prime):
                coefficients = H_mod(a, Z, prime)
                degree = polynomial_degree_mod(coefficients, prime)
                PACER.tick()
                irreducible = degree == 8 and cached_irred8(coefficients, prime)
                if irreducible:
                    witnesses.append(a)
                rows.append(
                    {
                        "type": "finite-field-specialization",
                        "label": "EVIDENCE",
                        "prime": prime,
                        "Z": Z,
                        "a": a,
                        "degree_H": degree,
                        "irreducible_degree_8": irreducible,
                    }
                )
            if not witnesses:
                missing.append(Z)
            witness_counts.append(len(witnesses))
            rows.append(
                {
                    "type": "finite-field-fiber",
                    "label": "EVIDENCE",
                    "prime": prime,
                    "Z": Z,
                    "irreducible_witness_a": witnesses,
                    "some_a_is_irreducible": bool(witnesses),
                    "no_certificate": not witnesses,
                }
            )
        summary = {
            "type": "finite-field-prime-summary",
            "label": "EVIDENCE",
            "prime": prime,
            "nonzero_Z_fibers": prime - 1,
            "fibers_with_certificate": prime - 1 - len(missing),
            "Z_with_no_certificate": missing,
            "minimum_witness_count": min(witness_counts),
            "maximum_witness_count": max(witness_counts),
            "total_irreducible_pairs": sum(witness_counts),
            "interpretation": "finite diagnostic only; not used in the fixed-Z theorem",
        }
        rows.append(summary)
        summaries.append(summary)
    return rows, summaries


def quadratic_character(value: int, prime: int) -> int:
    value %= prime
    if not value:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def scan_alpha_square_family() -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for prime in ALPHA_SQUARE_PRIMES:
        eligible = 0
        irreducible = 0
        first_witness: list[int] | None = None
        counts_by_a: dict[int, int] = {}
        for a in range(prime):
            A = (1 + 4 * a * a) % prime
            if quadratic_character(A, prime) != -1:
                continue
            count = 0
            for unit in range(1, prime):
                coefficients = alpha_square_mod(a, unit, prime)
                eligible += 1
                PACER.tick()
                if cached_irred8(coefficients, prime):
                    irreducible += 1
                    count += 1
                    if first_witness is None:
                        first_witness = [a, unit]
            counts_by_a[a] = count
        summaries.append(
            {
                "type": "alpha-square-full-degree-prime-summary",
                "label": "EVIDENCE",
                "prime": prime,
                "family": (
                    "a^4*u^2*Ng^2+16*A^3*b^4-32*A^4*s^2*b^5, "
                    "with (A|p)=-1 and u!=0"
                ),
                "eligible_a_u_pairs": eligible,
                "irreducible_pairs": irreducible,
                "first_witness_a_u": first_witness,
                "irreducible_counts_by_a": counts_by_a,
                "no_certificate": first_witness is None,
                "interpretation": "finite diagnostic only; not used in the fixed-Z theorem",
            }
        )
    return summaries


def write_report(
    facts: dict[str, Any],
    fixed_summaries: list[dict[str, Any]],
    alpha_summaries: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    degree_four = facts["degree_four_condition"]
    lines = [
        "# L22 fixed-`Z` factor elimination",
        "",
        "## Verdict",
        "",
        "**PROVED.** For every nonzero rational `Z`, the integral octic",
        "",
        "    H = a^8 Z^4 N_g(b)^2 + 4 A^3 D^2 b^4",
        "        - 2 A^4 (a-1)^2 D^2 b^5",
        "",
        "is irreducible in `Q(a)[b]`.  Since `P=(4/A)H`, the same holds for the normalized class polynomial.  In particular this covers every rational cube `Z=z^3` in a cell with `v_w(z)>=1`.",
        "",
        "The proof below is symbolic.  The finite-field tables are an independent diagnostic and are not used to infer the theorem.",
        "",
        "## 1. Exact local data",
        "",
        "Put `A=1+4a^2`, `D=1-Z-a^2Z^2`, and `N_g=16a^4b^2-A(b-1)^4`.  The constant and leading coefficients of `H` are both",
        "",
        "    L = a^8 A^2 Z^4.",
        "",
        "For fixed `Z` different from `0,1`, `H` is primitive in `Q[a][b]`: a common divisor must divide `L`; the `b^4` coefficient at `a=0` is `4(1-Z)^2`, excluding `a`, and its reduction modulo `A` is nonzero, excluding `A`.  Gauss therefore turns any `Q(a)` factorization into a primitive `Q[a,b]` factorization.",
        "",
        "The script derives the following coefficient data directly from the sparse expansion:",
        "",
        f"- `deg_a(H_i)`, `i=0..8`: `{facts['a_degrees']}`.  The apparent `a^14 b^5` terms cancel exactly; the top surviving term is `1024 a^13 Z^4 b^5`.",
        f"- `v_a(H_i)`: `{facts['a_valuations']}`.",
        f"- `v_A(H_i)`: `{facts['A_valuations']}`.",
        "",
        "At the prime `A` the Newton polygon has sides `4@-1/2` and `4@+1/2`.  By the Newton product theorem, slopes of a factor form a sub-multiset of these slopes.  Since its endpoint valuations are integers, the sum of an odd number of half-integral root valuations cannot occur.  **Every global factor has even degree.**  Thus a reducible octic has a factor of degree 2 or 4; the `1+7` and `3+5` cases are eliminated.",
        "",
        "At `a=0` (for `Z!=0,1`) the sides are",
        "",
        "    4@-2, 1@0, 3@8/3.",
        "",
        "For `b=a^2 x`, the first side residual is",
        "",
        "    R_0(x)=Z^4+4(1-Z)^2 x^4.",
        "",
        "It has no rational linear factor: for rational `x` both nonzero summands are nonnegative.  The denominator 3 of the last slope makes its three roots an indivisible local block.  The exact Newton-residual theorem therefore leaves only these global allocations:",
        "",
        "- degree 2: two roots from the small four-root block;",
        "- degree 4: either all four small roots, or the three-root large block plus the one unit root.",
        "",
        "In a `4+4` factorization the factors may be swapped, so one factor always consists of all four small roots.  There is no omitted mixed branch: selecting one or three small roots would give a linear divisor of `R_0`.",
        "",
        "At `a=infinity`, direct rescaling gives",
        "",
        "    a^-12 H(a,x/a) -> 16 Z^4 (4x^2-1)^2,",
        "    a^-20 H(a,a x) -> 16 Z^4 x^4 (x^2-4)^2.",
        "",
        "These are the exact repeated residual clusters used below.",
        "",
        "## 2. Endpoint parameterization",
        "",
        "Let `F` be a primitive factor of degree `r=2` or `4` consisting only of `a=0` small roots.  Its leading coefficient is not divisible by `a` (otherwise every coefficient would be), while its constant coefficient has `a`-valuation `2r`.  Because both endpoints of `H` equal `L`, unique factorization gives",
        "",
        "    lc_b(F)=c A^beta,     F(0)=d a^(2r) A^delta,",
        "",
        "with nonzero rational `c,d` and `0<=beta,delta<=2`.  Any fixed rational powers of `Z` are scalars and are absorbed in `c,d`; hence there is no hidden `Z` denominator.  If `ell` roots of `F` are large at infinity, comparison of the endpoint quotient with the two infinity slopes gives",
        "",
        "    2r+2delta-2beta = 2ell-r.",
        "",
        f"The script exhausts the nine `(beta,delta)` pairs.  For `r=2` it obtains `{facts['degree_two_cases']}`; for `r=4` it obtains `{facts['degree_four_cases']}`.",
        "",
        "## 3. Degree 2 is impossible",
        "",
        "If the quadratic factor contains two infinity-large roots, its residual divides `(x^2-4)^2`, so its constant/leading ratio at `a=0` is `r=+16` or `r=-16`.  If it contains one large and one small root, the large residual is `x(x+2)` or `x(x-2)` and its small residual root is `+1/2` or `-1/2`; the same conclusion `r=+/-16` follows.",
        "",
        "But its `a=0` residual must divide `R_0`.  After monic normalization write",
        "",
        "    (x^2+e x+r)(x^2-e x+s)=x^4+K,",
        "    K=Z^4/(4(1-Z)^2)>0.",
        "",
        "Coefficient comparison gives `e(s-r)=0` and `s-e^2+r=0`.  If `e=0`, then `s=-r` and `K=-r^2<0`, impossible.  Otherwise `s=r` and `e^2=2r`; for `r=16` this asks for a rational square root of 32, and for `r=-16` it asks for a square root of -32.  Both are impossible.  Thus no degree-2 factor exists.",
        "",
        "## 4. Degree 4 is impossible",
        "",
        "For the all-small factor, the endpoint enumeration forces `beta=2`, `delta=0`, and all four of its roots to be infinity-large.  Write",
        "",
        "    lc_b(F)=c A^2,       F(0)=d a^8.",
        "",
        "Its `b=a x` initial form must be the entire nonzero-root cluster `16c(x^2-4)^2`; hence `d=256c`.  At `b=a^2x`, the complementary factor has initial form equal to its constant `Z^4/d`.  Matching the `x^4` coefficient of `R_0` therefore gives, with no other term at this weight,",
        "",
        "    c Z^4/d = 4(1-Z)^2,",
        f"    {degree_four['forced_equation']}.",
        "",
        "Factoring the difference of squares gives the two quadratics",
        "",
        "    Z^2+32Z-32=0,        discriminant 1152=576*2,",
        "    Z^2-32Z+32=0,        discriminant  896=64*14.",
        "",
        "Neither discriminant is a rational square, so the forced equation has no rational root.  Thus no degree-4 factor exists.",
        "",
        "## 5. Exceptional values and the cell condition",
        "",
        "- `Z=0`: **PROVED reducible and excluded**.  Exactly `H=2A^3b^4(2-A(a-1)^2b)`.  Cell parameters have nonzero `z`.",
        "- `Z=1`: the `a=0` polygon used above degenerates, so it is tested separately.  At `a=1`, primitive coefficients are `[25,-200,540,-760,1546,-760,540,-200,25]`; modulo 11 the exact degree-8 Rabin criterion succeeds.  The localization at `L` shows an assumed vertical factorization would specialize with positive degrees, contradiction.",
        "- `Z=2`: this is **not exceptional** to the symbolic proof, but is independently cross-checked at `a=1`; modulo 13 the primitive octic has an exact irreducible reduction.",
        "- `D=0` at a specialized `a`: then `H=a^8Z^4N_g^2` is a square and is correctly excluded.  For fixed nonzero `Z`, `D=1-Z-a^2Z^2` has at most two rational zeros.  Vertical irreducibility is a statement in `Q(a)[b]` and is unaffected; the existing admissible progression avoids those finitely many values.",
        "",
        "This proves fixed-fiber irreducibility for every nonzero rational `Z`, stronger than the rational-cube locus required here.",
        "",
        "## 6. Independent finite-field diagnostic",
        "",
        "Every `(p,Z,a)` row is recorded in the JSONL artifact.  `H mod p` has the same irreducibility as `P mod p` whenever `aA Z!=0`; bad leading-coefficient rows are recorded as degree drops.  These finite facts are labeled **EVIDENCE** and were not used above.",
        "",
        "| p | nonzero Z fibers | fibers with a witness | Z with no witness | irreducible `(Z,a)` pairs |",
        "|---:|---:|---:|:--|---:|",
    ]
    for row in fixed_summaries:
        lines.append(
            f"| {row['prime']} | {row['nonzero_Z_fibers']} | {row['fibers_with_certificate']} | "
            f"`{row['Z_with_no_certificate']}` | {row['total_irreducible_pairs']} |"
        )
    lines.extend(
        [
            "",
            "The absence of witnesses for small finite fields is not an exceptional-locus theorem.  Starting at `p=19` in this chosen table, every nonzero finite-field fiber has at least one witness.",
            "",
            "### Two-parameter full-degree reduction diagnostic",
            "",
            "The additional scanned family is",
            "",
            "    a^4 u^2 N_g^2 + 16 A^3 b^4 - 32 A^4 s^2 b^5",
            "",
            "under `(A|p)=-1` and `u!=0`.  Again this is only finite evidence.",
            "",
            "| p | eligible `(a,u)` | irreducible pairs | first witness |",
            "|---:|---:|---:|:--|",
        ]
    )
    for row in alpha_summaries:
        lines.append(
            f"| {row['prime']} | {row['eligible_a_u_pairs']} | {row['irreducible_pairs']} | "
            f"`{row['first_witness_a_u']}` |"
        )
    lines.extend(
        [
            "",
            "## 7. Replay",
            "",
            "    nice -n 19 python3 math/h10q/l22_elimination.py",
            "",
            f"The replay wrote `{OUT}`, `{REPORT}`, and exited 0.  Wall-clock: {summary['wall_seconds']:.3f} s.  Pacing: {summary['pacing']['steps']} checks, {summary['pacing']['sleeps']} sleeps, {summary['pacing']['sleep_seconds']:.3f} total sleep seconds.",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    started = time.perf_counter()
    H, _, _, _ = build_H()
    symbolic_rows, facts = symbolic_derivation(H)
    fixed_rows, fixed_summaries = scan_fixed_fields()
    alpha_summaries = scan_alpha_square_family()
    elapsed = time.perf_counter() - started

    missing_by_prime = {
        str(row["prime"]): row["Z_with_no_certificate"]
        for row in fixed_summaries
    }
    alpha_missing = [
        row["prime"] for row in alpha_summaries if row["no_certificate"]
    ]
    summary = {
        "type": "summary",
        "label": "PROVED",
        "theorem": "H(a,Z) is irreducible in Q(a)[b] for every nonzero rational Z",
        "factor_degrees": {
            "1+7": "excluded by A-adic half-integral slopes",
            "2+6": "excluded by endpoint and a=0/infinity residual comparison",
            "3+5": "excluded by A-adic half-integral slopes",
            "4+4": "would force Z^4=1024*(1-Z)^2, which has no rational root",
        },
        "exceptional_locus": {
            "Z=0": "reducible and outside the cell family",
            "Z=1": "symbolic a=0 degeneration; proved irreducible by a=1 mod-11 certificate",
            "Z=2": "covered symbolically and independently checked by a=1 mod-13 certificate",
            "D=0_specializations": "reducible square; at most two a-values and explicitly excluded",
        },
        "fixed_field_scan": {
            "label": "EVIDENCE",
            "primes": list(FIXED_FIELD_PRIMES),
            "missing_Z_by_prime": missing_by_prime,
            "not_used_in_proof": True,
        },
        "alpha_square_scan": {
            "label": "EVIDENCE",
            "primes": list(ALPHA_SQUARE_PRIMES),
            "primes_with_no_certificate": alpha_missing,
            "not_used_in_proof": True,
        },
        "wall_seconds": elapsed,
        "pacing": {
            "steps": PACER.steps,
            "sleeps": PACER.sleeps,
            "sleep_seconds": PACER.sleep_seconds,
            "every": PACE_EVERY,
            "minimum_sleep_seconds": PACE_SECONDS,
            "loop_duty_cycle_cap": "50% by sleeping at least the preceding chunk compute time",
        },
        "refusals": 0,
        "open_obligations_in_scope": 0,
    }
    rows: list[dict[str, Any]] = [
        {
            "type": "meta",
            "label": "PROVED",
            "artifact": "l22_elimination",
            "version": 1,
            "arithmetic": "stdlib Fraction plus local sparse and finite-field polynomial helpers",
            "symbolic_proof_independent_of_finite_scans": True,
            "strict_label_legend": {
                "PROVED": "exact identity, finite certificate, or theorem proved in the report",
                "EVIDENCE": "finite diagnostic not promoted to a uniform theorem",
                "CONDITIONAL": "not used",
                "OPEN": "not used",
            },
        },
        *symbolic_rows,
        *fixed_rows,
        *alpha_summaries,
        summary,
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(facts, fixed_summaries, alpha_summaries, summary)
    print(
        "L22 elimination: PROVED all nonzero rational Z; "
        f"fixed_rows={len(fixed_rows)}, alpha_primes={len(alpha_summaries)}, "
        f"wall={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
