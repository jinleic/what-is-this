#!/usr/bin/env python3
"""Exact multivariable obstructions on the canonical tau_dagger branch.

This replay tests whether moving ``a`` together with the construction prime can
lower the prime-value problem. It uses only exact standard-library arithmetic.
Finite-field computations certify one displayed diagonal; they are not a
search-based uniform theorem.

Replay from the repository root:

    nice -n 19 python3 math/h10q/l23_multivar.py
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l23_multivar.jsonl"
REPORT = Path("/tmp/l23_multivar.md")
PACE_EVERY = 200
PACE_SECONDS = 0.002


class Pacer:
    """Keep finite certificate searches at deliberately low duty cycle."""

    def __init__(self) -> None:
        self.iterations = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.iterations += 1
        if self.iterations % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def trim(poly: list) -> list:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(left: list, right: list) -> list:
    zero = (left[0] if left else right[0]) * 0
    answer = [zero] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return trim(answer)


def poly_scale(poly: list, scalar: F | int) -> list:
    return trim([scalar * value for value in poly])


def poly_mul(left: list, right: list) -> list:
    answer = [left[0] * right[0] * 0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if not x:
            continue
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return trim(answer)


def poly_pow(poly: list, exponent: int) -> list:
    assert exponent >= 0
    answer = [poly[0] * 0 + 1]
    base = list(poly)
    while exponent:
        if exponent & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return answer


def poly_eval(poly: list, value: F | int) -> F | int:
    answer = poly[0] * 0
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def poly_degree(poly: list) -> int:
    return len(trim(poly)) - 1


def poly_divmod_q(left: list, right: list) -> tuple[list[F], list[F]]:
    remainder = trim([F(value) for value in left])
    divisor = trim([F(value) for value in right])
    assert divisor != [0]
    if len(remainder) < len(divisor):
        return [F(0)], remainder
    quotient = [F(0)] * (len(remainder) - len(divisor) + 1)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = trim(remainder)
    return trim(quotient), remainder


def poly_gcd_q(left: list, right: list) -> list[F]:
    a = trim([F(value) for value in left])
    b = trim([F(value) for value in right])
    while b != [0]:
        _, remainder = poly_divmod_q(a, b)
        a, b = b, remainder
    return [value / a[-1] for value in a]


def primitive_integer_poly(poly: list) -> tuple[list[int], F]:
    values = [F(value) for value in trim(poly)]
    denominator = math.lcm(*(value.denominator for value in values))
    integers = [
        value.numerator * (denominator // value.denominator)
        for value in values
    ]
    content = math.gcd(*(abs(value) for value in integers))
    assert content > 0
    integers = [value // content for value in integers]
    scale = F(content, denominator)
    if integers[-1] < 0:
        integers = [-value for value in integers]
        scale = -scale
    assert [scale * value for value in integers] == values
    return integers, scale


def fixed_divisor(poly: list[int]) -> int:
    """Exact fixed divisor: degree+1 consecutive values suffice."""
    answer = 0
    for value in range(poly_degree(poly) + 1):
        PACER.tick()
        answer = math.gcd(answer, abs(int(poly_eval(poly, value))))
    return answer


def tuple_fixed_divisor(polys: list[list[int]]) -> int:
    """Exact fixed divisor of the product of the displayed polynomials."""
    degree = sum(poly_degree(poly) for poly in polys)
    answer = 0
    for value in range(degree + 1):
        PACER.tick()
        product = 1
        for poly in polys:
            product *= int(poly_eval(poly, value))
        answer = math.gcd(answer, abs(product))
    return answer


def vp(value: F | int, prime: int) -> int:
    value = F(value)
    assert value
    numerator = abs(value.numerator)
    denominator = value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def mod_fraction(value: F | int, prime: int) -> int:
    value = F(value)
    return value.numerator * pow(value.denominator, -1, prime) % prime


def quadratic_character(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    result = pow(value, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


def canonical_polynomials(a: list, b: list, Z: F | int) -> dict[str, list]:
    """Construct A,D,Ng,H with H=(A/4)P on tau_dagger.

    H is polynomial even when a moves. The actual tied square class is
    P=4H/A, so the denominator A must not be silently discarded.
    """
    one = [a[0] * 0 + 1]
    A = poly_add(one, poly_scale(poly_pow(a, 2), 4))
    D = poly_add(
        poly_scale(one, 1 - F(Z)),
        poly_scale(poly_pow(a, 2), -F(Z) ** 2),
    )
    b_minus_one = poly_add(b, poly_scale(one, -1))
    Ng = poly_add(
        poly_scale(poly_mul(poly_pow(a, 4), poly_pow(b, 2)), 16),
        poly_scale(poly_mul(A, poly_pow(b_minus_one, 4)), -1),
    )
    H = poly_add(
        poly_add(
            poly_scale(poly_mul(poly_pow(a, 8), poly_pow(Ng, 2)), F(Z) ** 4),
            poly_scale(
                poly_mul(poly_mul(poly_pow(A, 3), poly_pow(D, 2)), poly_pow(b, 4)),
                4,
            ),
        ),
        poly_scale(
            poly_mul(
                poly_mul(
                    poly_mul(poly_pow(A, 4), poly_pow(poly_add(a, [-1]), 2)),
                    poly_pow(D, 2),
                ),
                poly_pow(b, 5),
            ),
            -2,
        ),
    )
    return {"A": A, "D": D, "Ng": Ng, "H": H}


def canonical_P_scalar(a: F, Z: F) -> list[F]:
    """Ascending coefficients of canonical P(b), exactly."""
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    b = [F(0), F(1)]
    n = [-A, 4 * A, 16 * a**4 - 6 * A, 4 * A, -A]
    first = poly_scale(poly_pow(n, 2), 4 * a**8 * Z**4 / A)
    second = [F(0)] * 4 + [16 * A * A * D * D]
    third = [F(0)] * 5 + [-32 * A**3 * s * s * D * D]
    answer = [F(value) for value in poly_add(poly_add(first, second), third)]
    direct = canonical_polynomials([a], b, Z)
    reconstructed = [4 * coefficient / A for coefficient in direct["H"]]
    assert answer == reconstructed
    return answer


def reciprocal_lift(trace: list[F]) -> list[F]:
    """Return b^4*trace(b+b^-1), with exponents 0 through 8."""
    answer = [F(0)] * 9
    for degree, coefficient in enumerate(trace):
        for chosen in range(degree + 1):
            exponent = 4 + degree - 2 * chosen
            answer[exponent] += coefficient * math.comb(degree, chosen)
    return trim(answer)


def trace_a_one(Z: F) -> tuple[list[F], list[F]]:
    """Monic F=P/lc and its reciprocal trace polynomial at a=1."""
    D = 1 - Z - Z * Z
    n_trace = [F(-4, 5), F(4), F(-1)]
    U = 20 * D * D / Z**4
    trace = poly_add(poly_pow(n_trace, 2), [U])
    P = canonical_P_scalar(F(1), Z)
    monic = [coefficient / P[-1] for coefficient in P]
    assert monic == reciprocal_lift(trace)
    return monic, trace


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        PACER.tick()
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def primes_up_to(limit: int) -> Iterable[int]:
    for candidate in range(2, limit + 1):
        PACER.tick()
        if is_prime(candidate):
            yield candidate


def ff_trim(poly: list[int], prime: int) -> list[int]:
    answer = [value % prime for value in poly]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def ff_divmod(left: list[int], right: list[int], prime: int) -> tuple[list[int], list[int]]:
    remainder = ff_trim(left, prime)
    divisor = ff_trim(right, prime)
    assert divisor != [0]
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % prime
        remainder = ff_trim(remainder, prime)
    return ff_trim(quotient, prime), remainder


def ff_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    a, b = ff_trim(left, prime), ff_trim(right, prime)
    while b != [0]:
        _, remainder = ff_divmod(a, b, prime)
        a, b = b, remainder
    inverse = pow(a[-1], -1, prime)
    return [(value * inverse) % prime for value in a]


def ff_mul_mod(left: list[int], right: list[int], modulus: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            product[i + j] = (product[i + j] + x * y) % prime
    return ff_divmod(product, modulus, prime)[1]


def ff_pow_mod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    answer = [1]
    power = ff_divmod(base, modulus, prime)[1]
    while exponent:
        if exponent & 1:
            answer = ff_mul_mod(answer, power, modulus, prime)
        power = ff_mul_mod(power, power, modulus, prime)
        exponent >>= 1
    return answer


def prime_divisors(value: int) -> list[int]:
    answer: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        PACER.tick()
        if value % divisor == 0:
            answer.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        answer.append(value)
    return answer


def irreducible_mod(poly: list[int], prime: int) -> bool:
    reduced = ff_trim(poly, prime)
    degree = poly_degree(poly)
    if poly_degree(reduced) != degree:
        return False
    inverse = pow(reduced[-1], -1, prime)
    modulus = [(value * inverse) % prime for value in reduced]
    x = [0, 1]
    checks = {degree // divisor for divisor in prime_divisors(degree)}
    frobenius = x
    for index in range(1, degree + 1):
        frobenius = ff_pow_mod(frobenius, prime, modulus, prime)
        difference = list(frobenius)
        if len(difference) < 2:
            difference.extend([0] * (2 - len(difference)))
        difference[1] = (difference[1] - 1) % prime
        difference = ff_trim(difference, prime)
        if index in checks and poly_degree(ff_gcd(modulus, difference, prime)) > 0:
            return False
    return ff_trim(difference, prime) == [0]


def no_linear_quadratic_factor_mod(poly: list[int], prime: int) -> bool:
    reduced = ff_trim(poly, prime)
    if poly_degree(reduced) != poly_degree(poly):
        return False
    inverse = pow(reduced[-1], -1, prime)
    modulus = [(value * inverse) % prime for value in reduced]
    xp2 = ff_pow_mod([0, 1], prime * prime, modulus, prime)
    if len(xp2) < 2:
        xp2.extend([0] * (2 - len(xp2)))
    xp2[1] = (xp2[1] - 1) % prime
    return poly_degree(ff_gcd(modulus, ff_trim(xp2, prime), prime)) == 0


def certificate_prime(poly: list[int], predicate, limit: int = 300) -> int | None:
    for prime in primes_up_to(limit):
        PACER.tick()
        if predicate(poly, prime):
            return prime
    return None


def degree_rows() -> list[dict]:
    """Mechanically check sharp cases in the general degree calculation."""
    q = [F(0), F(1)]
    cases = [
        ("fixed_a", [F(3)], poly_scale(q, 3), 0, 8, 4),
        ("linear_zero_constant", q, q, 1, 20, 6),
        ("linear_generic", [F(1), F(2)], poly_scale(q, 3), 1, 20, 6),
        ("linear_lead_cancel", [F(3, 2), F(1, 2)], q, 1, 19, 5),
        ("linear_double_cancel", [F(-1), F(1, 2)], q, 1, 19, 4),
        ("quadratic", [F(1), F(0), F(1)], q, 2, 36, 10),
        ("cubic", [F(1), F(0), F(0), F(1)], q, 3, 52, 14),
    ]
    rows: list[dict] = [
        {
            "type": "polynomial_degree_theorem",
            "label": "PROVED",
            "scope": "actual cell Z!=0,1; fixed kappa!=0; b=kappa*Q; and nonconstant a(Q) of degree r with leading coefficient c",
            "raw_degree_Ng": {
                "r>=2": "4*r+2",
                "r=1, 4*c^2!=kappa^2": 6,
                "r=1, 4*c^2=kappa^2": "<=5",
            },
            "raw_degree_H": {
                "r>=2": "16*r+4",
                "r=1, 4*c^2!=kappa^2": 20,
                "r=1, 4*c^2=kappa^2": 19,
            },
            "raw_same_squareclass_polynomial": "A*H, because (A*H)/P=(A/2)^2",
            "raw_degree_AH": {
                "r>=2": "18*r+4",
                "r=1, 4*c^2!=kappa^2": 22,
                "r=1, 4*c^2=kappa^2": 21,
            },
            "a_at_Q0_nonzero": "Q does not divide H; the displayed a(0)!=0 branches have raw degrees 21 or larger, but no universal squarefreeness/minimality is asserted",
            "a_at_Q0_zero": {
                "exact_valuation": "ord_Q(H)=4, since the Q^4 coefficient is 4*kappa^4*(1-Z)^2",
                "same_squareclass_representative": "A*(H/Q^4)",
                "degree_r>=2": "18*r",
                "degree_r=1_admissible": 18,
                "why_no_linear_cancellation": "v2(c)=v2(kappa)=0 makes 4*c^2=kappa^2 impossible",
            },
            "coprimality_identity": "H=256*Z^4*a^16*kappa^4*Q^4 mod A",
            "coprimality": "gcd(A,H)=1 because gcd(A,a)=gcd(A,Q)=1",
            "conclusion": "raw A*H is not always squarefree or minimal; the forced-square-reduced a(0)=0 representative has degree 18*r and the certified linear example has exact squarefree factor degrees 2+16",
        }
    ]
    for name, a, b, degree_a, expected_H, expected_Ng in cases:
        PACER.tick()
        data = canonical_polynomials(a, b, 27)
        actual_H = poly_degree(data["H"])
        actual_Ng = poly_degree(data["Ng"])
        assert actual_H == expected_H
        assert actual_Ng == expected_Ng
        rows.append(
            {
                "type": "degree_case",
                "label": "PROVED",
                "name": name,
                "degree_a_in_Q": degree_a,
                "degree_b_in_Q": poly_degree(b),
                "degree_Ng_in_Q": actual_Ng,
                "degree_H_numerator_in_Q": actual_H,
                "degree_A_factor_in_Q": poly_degree(data["A"]),
                "degree_same_squareclass_G_equals_AH": poly_degree(
                    poly_mul(data["A"], data["H"])
                ),
            }
        )
    return rows


def capell_rows() -> list[dict]:
    rows: list[dict] = []

    trace_checks = []
    for Z in (F(27), F(1, 2), F(-3, 5)):
        PACER.tick()
        _monic, trace = trace_a_one(Z)
        D = 1 - Z - Z * Z
        factor_64 = 125 * D * D + 64 * Z**4
        factor_1024 = 125 * D * D + 1024 * Z**4
        norm_minus = 16 * poly_eval(trace, F(2))
        norm_plus = 16 * poly_eval(trace, F(-2))
        square_multiplier = F(64, 25 * Z**4)
        assert norm_minus == square_multiplier * factor_64
        assert norm_plus == square_multiplier * factor_1024
        trace_checks.append(
            {
                "Z": frac_text(Z),
                "factor_64": frac_text(factor_64),
                "factor_1024": frac_text(factor_1024),
            }
        )
    parity = []
    for p_mod_2, q_mod_2 in ((0, 1), (1, 0), (1, 1)):
        d_mod_2 = (q_mod_2 - p_mod_2 * q_mod_2 - p_mod_2) % 2
        assert d_mod_2 == 1
        residues = {
            str(coefficient): (125 * d_mod_2**2 + coefficient * p_mod_2**4) % 8
            for coefficient in (64, 1024)
        }
        assert set(residues.values()) == {5}
        parity.append(
            {
                "p_mod_2": p_mod_2,
                "q_mod_2": q_mod_2,
                "d_mod_2": d_mod_2,
                "numerators_mod_8": residues,
            }
        )
    rows.append(
        {
            "type": "capell_reciprocal",
            "label": "PROVED",
            "scope": "a=1 and every rational Z!=0",
            "criterion": "if 2*b were square in K, one of 2*(t-2), 2*(t+2) would be square in the reciprocal quartic subfield",
            "norm_squareclasses": [
                "125*(1-Z-Z^2)^2+64*Z^4",
                "125*(1-Z-Z^2)^2+1024*Z^4",
            ],
            "parity_obstruction": "for Z=p/q reduced, d=q^2-p*q-p^2 is odd and both cleared numerators are 5 mod 8",
            "square_locus": "EMPTY",
            "mechanical_trace_checks": trace_checks,
            "primitive_parity_classes": parity,
        }
    )

    # On every actual cell Z=z^3.  The left target-prime Newton edge has
    # endpoints (0, 12*e_z) and (4, 0).  Its residual is separable, so its
    # local factors are unramified and ord_w(b)=3*e_z.
    edge_w = 3
    edge_z = F(3)
    edge_Z = edge_z**3
    edge_P = canonical_P_scalar(F(5), edge_Z)
    edge_values = [vp(value, edge_w) for value in edge_P]
    assert edge_values == [12, 12, 16, 12, 0, 0, 16, 12, 12]
    assert all(
        edge_values[index] > 12 - 3 * index
        for index in (1, 2, 3)
    )
    edge_constant = mod_fraction(edge_P[0] / edge_w**12, edge_w)
    edge_leading = mod_fraction(edge_P[4], edge_w)
    edge_residual = [edge_constant, 0, 0, 0, edge_leading]
    edge_derivative = [0, 0, 0, 4 * edge_leading]
    assert ff_gcd(edge_residual, edge_derivative, edge_w) == [1]
    rows.append(
        {
            "type": "capell_target_newton_edge",
            "label": "PROVED",
            "scope": "actual cell Z=z^3 with odd target w and w-adic units a,A,D",
            "edge": "(0,12*e_z)-(4,0), where e_z=v_w(z)",
            "root_valuation": "v_w(b)=3*e_z on four roots",
            "residual": "c0+c4*y^4 with c0*c4 nonzero; separable because w is odd",
            "ramification_index": 1,
            "obstruction": "if e_z is odd, ord_w(2*b)=3*e_z is odd",
            "square_locus_when_e_z_odd": "EMPTY",
            "mechanical_check": {
                "w": edge_w,
                "z": frac_text(edge_z),
                "coefficient_valuations": edge_values,
                "residual": edge_residual,
            },
        }
    )

    w = 3
    a = F(5)
    Z = F(27)
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    assert quadratic_character(A.numerator, w) == -1
    assert mod_fraction(s, w) != 0 and mod_fraction(Z, w) == 0
    P = canonical_P_scalar(a, Z)
    reduced = [mod_fraction(value, w) for value in P]
    expected = [0] * 9
    expected[4] = 16 * mod_fraction(A * A, w) % w
    expected[5] = -32 * mod_fraction(A**3 * s * s, w) % w
    assert reduced == expected
    root = pow(mod_fraction(2 * A * s * s, w), -1, w)
    derivative = (expected[4] * 4 * pow(root, 3, w) + expected[5] * 5 * pow(root, 4, w)) % w
    assert poly_eval(reduced, root) % w == 0 and derivative != 0
    assert quadratic_character(2 * root, w) == -1
    rows.append(
        {
            "type": "capell_aligned_bad_root",
            "label": "PROVED",
            "scope": "odd target prime w with v_w(Z)>0; a,A,D,s are w-adic units; and Legendre(A,w)=-1",
            "reduction": "P(b)=16*A^2*b^4*(1-2*A*s^2*b) mod w",
            "simple_root": "r=(2*A*s^2)^(-1)",
            "bad_unit": "2*r=(A*s^2)^(-1) has Legendre symbol -1",
            "square_locus": "EMPTY in this aligned stratum",
            "sharp_diagonal_check": {
                "w": w,
                "a_mod_w": mod_fraction(a, w),
                "Z_mod_w": mod_fraction(Z, w),
                "root": root,
                "derivative": derivative,
                "Legendre_2r": quadratic_character(2 * root, w),
            },
        }
    )
    rows.append(
        {
            "type": "capell_scope",
            "label": "OPEN",
            "proved_empty_strata": [
                "a=1, all rational Z!=0",
                "actual cells with odd v_w(z), for every w-unit a",
                "aligned choices with s=(a-1)/2 nonzero mod target w",
            ],
            "class_selection_consequence": "when v_w(z) is even, L20 may refine its character choice to a!=1 mod w, placing the selected admissible a in the simple-bad-root stratum",
            "selected_protocol_status": "EMPTY Capell square locus after that permitted a-refinement",
            "unresolved_stratum": "as a statement about every possible a: a!=1, even v_w(z), and w dividing s",
            "warning": "Norm_K/Q(2*b)=256 is already a square, so endpoint/norm arguments alone do not obstruct Capell factorization",
        }
    )
    return rows


def natural_locus_rows() -> list[dict]:
    return [
        {
            "type": "natural_locus",
            "label": "PROVED",
            "locus": "degree_b(P)<8",
            "equation": "4*a^8*(1+4*a^2)*Z^4=0",
            "admissible_result": "EMPTY because odd a and a cell both require a*Z!=0",
        },
        {
            "type": "natural_locus",
            "label": "PROVED",
            "locus": "Ng=0",
            "equation": "16*a^4*b^2=A*(b-1)^4",
            "consequence": "A=(4*a^2*b/(b-1)^2)^2",
            "admissible_result": "EMPTY: for v2(a)=0, A=1+4*a^2 is 5 mod 8 and is not a rational square",
        },
        {
            "type": "natural_locus",
            "label": "PROVED",
            "locus": "L=1-2*A*s^2*b=0",
            "consequence": "b=1/(2*A*s^2) and v2(b)=-1-2*v2(s)",
            "admissible_result": "EMPTY under v2(b)=0",
        },
        {
            "type": "natural_locus",
            "label": "PROVED",
            "locus": "D=0 or Z=0",
            "factorization": "D=0 makes P=(4*a^8*Z^4/A)*Ng^2; Z=0 makes P=16*A^2*D^2*b^4*L",
            "admissible_result": "EXCLUDED, not a reduction: D is the denominator of c and Z=0 is not a target cell",
        },
        {
            "type": "norm_match_wall",
            "label": "PROVED",
            "identity": "H=X^2+A*L*Y^2",
            "required_squareclass": "-2*A*b*L must be a rational square to make this a norm from Q(sqrt(2*b))",
            "obstruction": "v2(-2*A*b*L)=1 for v2(a)=v2(b)=0",
            "extra_parameter_count": 1,
            "extra_parameter": "a putative square root; it cannot be treated as free and the parity wall leaves no point",
        },
        {
            "type": "constant_twist_audit",
            "label": "PROVED",
            "twisted_match": "a norm from Q(sqrt(2*j*b)) requires -2*A*j*b*L square",
            "two_adic_condition": "v2(j) must be odd",
            "desired_symbol_identity": "(P,2*b)=(P,2*j*b)*(P,j)",
            "conclusion": "the twisted norm identity alone is insufficient: a nonsquare twist leaves the separate global obligation (P,j)=1",
            "parameter_count": "j is one additional datum unless fixed in advance; either way its Hilbert condition must be independently proved",
        },
    ]


def sharp_diagonal_row() -> dict:
    """Certify the sharp raw degree on one a(0)!=0 diagonal."""
    Q = [7, 12]
    a = [5, 6]                 # a=(Q+3)/2
    b = list(Q)
    Z = 27
    data = canonical_polynomials(a, b, Z)
    H, A = data["H"], data["A"]
    G_squareclass = poly_mul(A, H)
    assert poly_degree(data["Ng"]) == 5
    assert poly_degree(H) == 19 and poly_degree(A) == 2
    assert poly_degree(G_squareclass) == 21
    assert poly_gcd_q(H, A) == [F(1)]

    H_integer, H_scale = primitive_integer_poly(H)
    A_integer, A_scale = primitive_integer_poly(A)
    Q_integer, Q_scale = primitive_integer_poly(Q)
    H_fixed = fixed_divisor(H_integer)
    A_fixed = fixed_divisor(A_integer)
    Q_fixed = fixed_divisor(Q_integer)
    product_fixed = tuple_fixed_divisor([Q_integer, A_integer, H_integer])
    assert H_scale == -12 and A_scale == Q_scale == 1
    # Scaling every value by the same integer scales its exact fixed divisor.
    H_actual_fixed = abs(H_scale.numerator) * H_fixed
    product_actual_fixed = abs(H_scale.numerator) * product_fixed
    assert H_actual_fixed == product_actual_fixed == 12

    low_factor_prime = certificate_prime(H_integer, no_linear_quadratic_factor_mod)
    full_irred_prime = certificate_prime(H_integer, irreducible_mod)
    A_discriminant = A[1] * A[1] - 4 * A[0] * A[2]
    assert A_discriminant == -576
    assert low_factor_prime is not None and full_irred_prime is not None
    assert all(int(poly_eval(Q, t)) % 2 for t in range(6))
    assert all(int(poly_eval(a, t)) % 2 for t in range(6))
    assert int(poly_eval(a, 0)) % 3 == 2
    assert int(poly_eval(poly_scale(poly_add(a, [-1]), F(1, 2)), 0)) % 3 == 2

    return {
        "type": "sharp_polynomial_diagonal",
        "label": "PROVED",
        "cell": {"w": 3, "z": 3, "Z": 27},
        "parameter_contract": {
            "Q(t)": "7+12*t",
            "b(t)": "Q(t)",
            "a(t)": "(Q(t)+3)/2=5+6*t",
            "new_free_parameters": 0,
            "identification": "t is exactly the existing moving-prime coordinate",
        },
        "admissibility": {
            "v2_a": 0,
            "v2_b": 0,
            "A_character_mod_3": -1,
            "s_nonzero_mod_3": True,
        },
        "degrees": {
            "Ng": poly_degree(data["Ng"]),
            "H_reduced_numerator": poly_degree(H),
            "A_squareclass_factor": poly_degree(A),
            "same_squareclass_G_equals_AH": poly_degree(G_squareclass),
            "Q": 1,
        },
        "coprime_H_A": True,
        "primitive_scales": {
            "H": frac_text(H_scale),
            "A": frac_text(A_scale),
            "Q": frac_text(Q_scale),
            "frozen_H_scale_squareclass": -3,
        },
        "fixed_divisors": {
            "H_actual": H_actual_fixed,
            "product_Q_A_H_actual": product_actual_fixed,
            "H_primitive": H_fixed,
            "A_primitive": A_fixed,
            "Q": Q_fixed,
            "product_Q_A_H_primitive": product_fixed,
            "method": "gcd of values at 0 through the exact product degree",
        },
        "factor_certificates": {
            "A_quadratic_discriminant": A_discriminant,
            "A_irreducible_over_Q": True,
            "H_no_Q_linear_or_quadratic_factor_mod_prime": low_factor_prime,
            "H_full_irreducibility_mod_prime": full_irred_prime,
            "H_full_irreducibility_status": "PROVED" if full_irred_prime is not None else "OPEN (not needed for the degree obstruction)",
            "G_factor_degrees": [2, 19],
        },
        "effective_prime_value_complexity": "P and A*H have the same square class; G=A*H factors into a coprime irreducible quadratic and irreducible degree-19 polynomial, in addition to linear Q",
        "capell_status": "PROVED impossible by the w=3 simple bad root",
    }


def zero_constant_diagonal_row() -> dict:
    """Certify the forced Q^4 square on the admissible a(0)=0 branch."""
    Q = [7, 12]
    a = list(Q)                  # as a polynomial in abstract Q: a(Q)=Q
    b = list(Q)
    Z = 27
    data = canonical_polynomials(a, b, Z)
    H, A = data["H"], data["A"]
    Q4 = poly_pow(Q, 4)
    H_reduced, remainder = poly_divmod_q(H, Q4)
    assert remainder == [F(0)]
    _, next_remainder = poly_divmod_q(H_reduced, Q)
    assert next_remainder != [F(0)]
    assert poly_degree(H) == 20
    assert poly_degree(H_reduced) == 16
    assert poly_degree(poly_mul(A, H)) == 22
    assert poly_degree(poly_mul(A, H_reduced)) == 18
    assert poly_gcd_q(A, H_reduced) == [F(1)]

    H_integer, H_scale = primitive_integer_poly(H_reduced)
    A_integer, A_scale = primitive_integer_poly(A)
    Q_integer, Q_scale = primitive_integer_poly(Q)
    assert H_scale == 4 and A_scale == Q_scale == 1
    H_fixed = fixed_divisor(H_integer)
    A_fixed = fixed_divisor(A_integer)
    Q_fixed = fixed_divisor(Q_integer)
    product_fixed = tuple_fixed_divisor([Q_integer, A_integer, H_integer])
    assert H_fixed == A_fixed == Q_fixed == product_fixed == 1
    full_irred_prime = certificate_prime(H_integer, irreducible_mod, 500)
    assert full_irred_prime == 43
    A_discriminant = A[1] * A[1] - 4 * A[0] * A[2]
    assert A_discriminant == -2304

    assert all(int(poly_eval(Q, t)) % 2 for t in range(6))
    assert quadratic_character(int(poly_eval(A, 0)), 3) == -1
    s = poly_scale(poly_add(a, [-1]), F(1, 2))
    assert int(poly_eval(s, 0)) % 3 == 0

    return {
        "type": "zero_constant_polynomial_diagonal",
        "label": "PROVED",
        "cell": {"w": 3, "z": 3, "Z": 27, "v_w_z": 1},
        "parameter_contract": {
            "Q(t)": "7+12*t",
            "a_as_polynomial_in_Q": "a(Q)=Q, so a(0)=0",
            "a(t)": "7+12*t",
            "b(t)": "Q(t)",
            "new_free_parameters": 0,
        },
        "admissibility": {
            "v2_a": 0,
            "v2_b": 0,
            "A_character_mod_3": -1,
            "s_mod_3": 0,
            "capell_obstruction": "odd v_3(z) target Newton edge",
        },
        "forced_square": {
            "exact_Q_valuation_H": 4,
            "identity_source": "the Q^4 coefficient is 4*kappa^4*(1-Z)^2",
            "raw_same_squareclass_polynomial": "A*H",
            "raw_degree": 22,
            "reduced_same_squareclass_polynomial": "A*(H/Q^4)",
            "reduced_degree": 18,
        },
        "primitive_scales": {
            "H_over_Q4": frac_text(H_scale),
            "square_scale_removed": 4,
            "A": frac_text(A_scale),
            "Q": frac_text(Q_scale),
        },
        "fixed_divisors": {
            "H_over_Q4_actual": 4,
            "product_Q_A_H_over_Q4_actual": 4,
            "H_over_Q4_primitive": H_fixed,
            "A_primitive": A_fixed,
            "Q": Q_fixed,
            "primitive_tuple": product_fixed,
        },
        "factor_certificates": {
            "A_quadratic_discriminant": A_discriminant,
            "A_irreducible_over_Q": True,
            "H_over_Q4_primitive_degree": 16,
            "H_over_Q4_irreducible_mod_prime": full_irred_prime,
            "reduced_squareclass_factor_degrees": [2, 16],
            "reduced_representative_squarefree": True,
        },
        "conclusion": "the exact squarefree representative has degree 18, still strictly above the fixed-a octic degree 8",
    }


def newton_audit_row() -> dict:
    """Record the b5 cancellation that limits a tempting 2-adic shortcut."""
    a, Z = F(3), F(8)
    valuations = [vp(value, 2) for value in canonical_P_scalar(a, Z)]
    assert valuations == [14, 17, 16, 17, 4, 5, 16, 17, 14]
    return {
        "type": "newton_audit",
        "label": "PROVED",
        "parameters": {"a": 3, "Z": 8, "m=v2(Z)": 3, "e=v2((a-1)/2)": 0},
        "coefficient_v2_P_b0_through_b8": valuations,
        "conclusion": "the b^5 perturbation is a lower-hull point, so the three-vertex polygon (0,14)-(4,4)-(8,14) is false here",
        "general_p5_warning": "for m>=1 the two b^5 summands have valuations 5+4*m and 5+2*e; equality can also cancel",
        "m_le_0_warning": "the horizontal hull has residual (b+1)^8 mod 2, so root units alone do not prove odd normalized valuations",
        "uniform_2_adic_capell_verdict": "OPEN",
    }


def render_report(rows: list[dict]) -> str:
    sharp = next(row for row in rows if row["type"] == "sharp_polynomial_diagonal")
    zero = next(row for row in rows if row["type"] == "zero_constant_polynomial_diagonal")
    newton = next(row for row in rows if row["type"] == "newton_audit")
    degree_cases = [row for row in rows if row["type"] == "degree_case"]
    case_lines = [
        f"| `{row['name']}` | {row['degree_a_in_Q']} | {row['degree_Ng_in_Q']} | {row['degree_H_numerator_in_Q']} | {row['degree_A_factor_in_Q']} | {row['degree_same_squareclass_G_equals_AH']} |"
        for row in degree_cases
    ]
    fixed = sharp["fixed_divisors"]
    certs = sharp["factor_certificates"]
    zero_fixed = zero["fixed_divisors"]
    zero_certs = zero["factor_certificates"]
    return "\n".join(
        [
            "# L23 multivariable prime-complexity audit",
            "",
            "**Final status:** **PROVED obstruction for the tested families.** No unconditional member-existence theorem is obtained. Raw `A*H` has degree at least 21 on the tested `a(0)!=0` branches, but it is not asserted squarefree or minimal. On the missing admissible `a(0)=0` branch, exactly `Q^4||H`; removing that forced square gives a certified squarefree degree-18 representative with irreducible factor degrees 2 and 16, still above the fixed-a octic. **OPEN as a theorem for every possible `a`:** the residual Capell stratum `a!=1`, even `v_w(z)`, and `w|(a-1)/2`.",
            "",
            "## 1. Canonical identities and parameter accounting",
            "",
            "On `tau_dagger`, put `A=1+4a^2`, `D=1-Z-a^2Z^2`, `s=(a-1)/2`, `Ng=16a^4b^2-A(b-1)^4`, and `L=1-2As^2b`. The exact polynomial numerator is",
            "",
            "```text",
            "H=(A/4)P=a^8 Z^4 Ng^2+4A^3D^2b^4-2A^4(a-1)^2D^2b^5",
            "             =X^2+A*L*Y^2.",
            "```",
            "",
            "When `a` moves, the tied square class is still `P=4H/A`; dropping the now-variable denominator `A` changes the quaternion. Multiplying by the square `A^2/4` shows that the polynomial `G_move=A*H` represents exactly the same square class. Every curve below either uses the existing moving-prime coordinate `Q` as its sole parameter or explicitly counts an added square-root parameter. A nonconstant `a(Q)` also leaves the already-proved fixed aligned class, whose modulus freezes `a`; the diagonal is therefore an exact complexity test, not a reuse of that class theorem.",
            "",
            "## 2. Capell/Fable square criterion",
            "",
            "For irreducible `P`, `2b` is a square in `K=Q[b]/(P)` exactly when the monic normalization of `P(u^2/2)` splits as `R(u)R(-u)` with both factors of degree 8. Its norm is automatically `Norm_K/Q(2b)=256`; this is not an obstruction.",
            "",
            "### Reciprocal stratum `a=1`: empty",
            "",
            "Set `t=b+b^-1`. The monic octic is `b^4 T(t)`, where",
            "",
            "```text",
            "T(t)=(-t^2+4t-4/5)^2+20(1-Z-Z^2)^2/Z^4.",
            "```",
            "",
            "If `u^2=2b` in `K`, the involution `b -> b^-1` sends `u` to an element whose product with `u` is `+2` or `-2`. Hence one of `2(t+2)` and `2(t-2)` is a square in the quartic trace field. Taking norms forces one of",
            "",
            "```text",
            "125(1-Z-Z^2)^2+64Z^4,",
            "125(1-Z-Z^2)^2+1024Z^4",
            "```",
            "",
            "to be a rational square. For reduced `Z=p/q`, `d=q^2-pq-p^2` is odd in all three primitive parity classes, and each cleared numerator is `5 mod 8`. Both are nonsquares. Thus the square locus is **PROVED EMPTY for every rational `Z!=0`**.",
            "",
            "### Actual-cell target Newton edge",
            "",
            "Write `Z=z^3` and `e_z=v_w(z)>0`. For `w`-unit `a,A,D`, the left lower edge of `P` is `(0,12e_z)-(4,0)`. Its residual has the form `c0+c4*y^4` with nonzero endpoints and is separable because `w` is odd. Thus its local factors are unramified and four roots have `v_w(b)=3e_z`. If `e_z` is odd, `ord_w(2b)` is odd and the Capell square locus is **PROVED EMPTY**, without any assumption on `s`.",
            "",
            "### Nonreciprocal aligned stratum with `w` not dividing `s`: empty",
            "",
            "At the target prime, `Z=0 mod w`, `D=1 mod w`, and",
            "",
            "```text",
            "P(b)=16A^2b^4(1-2As^2b) (mod w).",
            "```",
            "",
            "If `s` is a unit, `r=(2As^2)^-1` is a simple root. Its bad unit is `2r=(As^2)^-1`, with character `(A|w)=-1`. Hensel gives a degree-one completion in which `2b` is nonsquare, so Capell factorization is impossible. This is **PROVED** and applies identically to the sharp diagonal at `w=3`. When `v_w(z)` is even, L20 can refine its character choice to `a!=1 mod w`, so the selected protocol again has empty Capell square locus.",
            "",
            "**OPEN only as a statement about every possible choice of `a`:** `a!=1`, even `v_w(z)`, and `w|s`. No conclusion is inferred there from a finite scan or from the incomplete 2-adic Newton shortcut.",
            "",
            "## 3. Natural factor, norm, and twist loci",
            "",
            "- `deg_b(P)<8` requires `aZ=0`, forbidden by oddness and the nonzero target cell.",
            "- `Ng=0` forces `A=(4a^2b/(b-1)^2)^2`, while odd `a` gives `A=5 mod 8` in `Q_2`.",
            "- `L=0` forces `v2(b)=-1-2v2(s)`, not zero.",
            "- `D=0` makes `P` a square but destroys the denominator defining `c`; `Z=0` is not a target cell.",
            "- Matching `H=X^2+ALY^2` to a norm from `Q(sqrt(2b))` requires `-2AbL` square, but its 2-adic valuation is exactly 1. The putative square root is one added parameter, not free branch data.",
            "- A constant twist `j` requires odd `v2(j)` to evade that wall, hence is nonsquare; then `(P,2b)=(P,2jb)(P,j)`. The twisted norm identity alone is insufficient because `(P,j)=1` remains a separate global obligation.",
            "",
            "## 4. Polynomial `a(Q)`: raw degree versus squareclass degree",
            "",
            "Fix actual cell data `Z!=0,1` and write `b=kappa*Q` with fixed nonzero `kappa`. If `a(Q)` has degree `r>=1` and leading coefficient `c`, exact leading-term comparison gives",
            "",
            "```text",
            "r>=2: deg Ng=4r+2, deg H=16r+4, deg(AH)=18r+4;",
            "r=1 and 4c^2!=kappa^2: deg Ng=6, deg H=20, deg(AH)=22;",
            "r=1 and 4c^2=kappa^2: deg Ng<=5, deg H=19, deg(AH)=21.",
            "```",
            "",
            "For the displayed `a(0)!=0` branches, `Q` does not divide `H`; these are exact **raw** degrees only. No claim that `A*H` is always squarefree or a minimal squareclass representative is made.",
            "",
            "If instead `a(0)=0`, then `A(0)=1`, `D(0)=1-Z`, and the `4A^3D^2b^4` term proves",
            "",
            "```text",
            "ord_Q(H)=4, with leading quotient coefficient 4*kappa^4*(1-Z)^2.",
            "```",
            "",
            "Thus `A*(H/Q^4)` represents the same square class. Its degree is exactly `18r` for `r>=2`; for admissible linear `a=cQ`, oddness gives `v2(c)=v2(kappa)=0`, so the exceptional equation `4c^2=kappa^2` is impossible and the degree is 18. This corrects the raw degree-21 minimum claim.",
            "",
            "The identity `H = 256 Z^4 a^16 kappa^4 Q^4 (mod A)` still proves `gcd(A,H)=1`; it does not remove the separate square factor `Q^4` inside `H`.",
            "",
            "| case | deg a | deg Ng | raw deg H | deg A | raw deg A*H |",
            "|---|---:|---:|---:|---:|---:|",
            *case_lines,
            "",
            "## 5. Certified `a(0)=0` branch: exact degree 18",
            "",
            "Use the same cell and prime coordinate but set",
            "",
            "```text",
            "Q(t)=7+12t,  a(Q)=Q,  b(Q)=Q,  Z=27.",
            "```",
            "",
            "Here `a,b` are odd, `(A|3)=-1`, and `v_3(z)=1`, so the target Newton edge supplies the Capell obstruction even though `3|s`. Exact division gives `H=Q^4 H_0`, `deg(H_0)=16`, and raw `deg(AH)=22`; the same-squareclass representative `A H_0` has degree 18. Moreover `H_0=4H_{0,prim}`, and 4 is a square.",
            "",
            "```text",
            f"actual: fd(H_0)={zero_fixed['H_over_Q4_actual']}, fd(Q*A*H_0)={zero_fixed['product_Q_A_H_over_Q4_actual']};",
            f"primitive: fd(Q)={zero_fixed['Q']}, fd(A)={zero_fixed['A_primitive']}, fd(H_0,prim)={zero_fixed['H_over_Q4_primitive']}, fd(tuple)={zero_fixed['primitive_tuple']}.",
            "```",
            "",
            f"`A` has discriminant `{zero_certs['A_quadratic_discriminant']}`; `H_0,prim` is irreducible modulo `{zero_certs['H_over_Q4_irreducible_mod_prime']}`. They are coprime, so the degree-18 representative is squarefree with exact factor degrees `{zero_certs['reduced_squareclass_factor_degrees']}`. This is the smallest certified moving-`a` representative in this audit, not a universal minimum over all rational substitutions.",
            "",
            "## 6. Certified `a(0)!=0` sharp raw-degree branch",
            "",
            "Use the cell `w=3,z=3,Z=27` and identify the sole parameter with the existing prime coordinate:",
            "",
            "```text",
            "Q(t)=7+12t,  b(t)=Q(t),  a(t)=(Q(t)+3)/2=5+6t.",
            "```",
            "",
            "For every integer `t`, `a,b` are odd, `A` is a nonresidue mod 3, and `s` is nonzero mod 3. Exact construction gives `deg Ng=5`, `deg H=19`, `deg A=2`, `deg G=A*H=21`, and `gcd(H,A)=1`. Because `a` moves, this is not one fixed aligned class; it is the requested test using the same six-variable formula and no extra coordinate. Degree-plus-one value gcds give",
            "",
            "```text",
            f"actual: fd(H)={fixed['H_actual']}, fd(Q*A*H)={fixed['product_Q_A_H_actual']};",
            f"after extracting H=-12*H_prim: fd(Q)={fixed['Q']}, fd(A)={fixed['A_primitive']}, fd(H_prim)={fixed['H_primitive']}, fd(Q*A*H_prim)={fixed['product_Q_A_H_primitive']}.",
            "```",
            "",
            f"`A` has discriminant `{certs['A_quadratic_discriminant']}` and is irreducible. Reduction of `H_prim` modulo `{certs['H_no_Q_linear_or_quadratic_factor_mod_prime']}` has no factor of degree 1 or 2; reduction modulo `{certs['H_full_irreducibility_mod_prime']}` is irreducible. Thus `G=A*H=-12*A*H_prim` has frozen constant squareclass `-3` and exact nonconstant factor degrees `{certs['G_factor_degrees']}`, not a tuple of only linear/quadratic forms.",
            "",
            "## 7. Newton-polygon scope correction",
            "",
            f"For `a=3,Z=8`, the exact coefficient valuations are `{newton['coefficient_v2_P_b0_through_b8']}`. The `b^5` point `(5,5)` lies below the segment from `(4,4)` to `(8,14)`. In general its two summands have valuations `5+4m` and `5+2e`; equality may cancel. For `m<=0`, the horizontal residual is `(b+1)^8`. Consequently the tempting universal three-vertex/odd-valuation proof is **OPEN**, not used above.",
            "",
            "## Conclusion",
            "",
            "**PROVED obstruction for the tested family:** joint polynomial motion of `a` does not beat the fixed degree-8 octic. On the `a(0)=0` branch, removing the exact forced square `Q^4` lowers raw degree 22 to a certified squarefree degree-18 representative with factor degrees 2 and 16. The separate `a(0)!=0` diagonal has raw degree 21 and certified factor degrees 2 and 19; that claim is scoped only to that diagonal. The natural square/norm/degree-drop loci are empty or invalid. The Capell route is ruled out for `a=1`, for actual cells with odd `v_w(z)`, and for the aligned `s mod w !=0` stratum; the permitted L20 refinement therefore rules it out for the selected protocol. **No tested locus reaches a known unconditional theorem sufficient for the tied split; a generic bounded-almost-prime statement would not control the individual Hilbert signs.**",
            "",
        ]
    )


def main() -> int:
    started = time.perf_counter()
    rows: list[dict] = [
        {
            "type": "meta",
            "artifact": "l23_multivar",
            "version": 5,
            "labels": {
                "symbolic_identities": "PROVED",
                "tested_family_obstructions": "PROVED",
                "finite_field_certificates": "PROVED for the displayed polynomial",
                "unresolved_capell_stratum": "OPEN",
                "unconditional_member_existence": "OPEN",
            },
            "branch": "tau_dagger=(1+2*a^2)/(1+4*a^2)",
            "constraints": ["v2(a)=0", "v2(b)=0", "Z!=0", "D!=0"],
            "authority_boundary": "finite scans and refusals are never negative evidence",
        }
    ]
    rows.extend(capell_rows())
    rows.extend(natural_locus_rows())
    rows.extend(degree_rows())
    rows.append(sharp_diagonal_row())
    rows.append(zero_constant_diagonal_row())
    rows.append(newton_audit_row())
    rows.append(
        {
            "type": "summary",
            "label": "PROVED",
            "proved_result": "exact obstruction/degree lower bound for the tested multivariable families",
            "raw_AH_degree_on_certified_a_at_Q0_nonzero_diagonal": 21,
            "forced_square_reduced_degree_on_a_at_Q0_zero_diagonal": 18,
            "smallest_certified_squarefree_representative_degree": 18,
            "factor_degrees_on_a_at_Q0_zero_diagonal": [2, 16],
            "factor_degrees_on_sharp_a_at_Q0_nonzero_diagonal": [2, 19],
            "unconditional_theorem_reduction": False,
            "capell_proved_empty": [
                "a=1 for all rational nonzero Z",
                "actual cells with odd v_w(z)",
                "aligned s nonzero mod target w",
            ],
            "capell_selected_protocol": "EMPTY after the permitted a!=1 mod w refinement when v_w(z) is even",
            "capell_open": "every-a statement on a!=1, even v_w(z), and target w dividing s",
            "extra_free_parameters_used": 0,
            "wall_seconds": round(time.perf_counter() - started, 6),
            "pacer_iterations": PACER.iterations,
            "pacer_sleeps": PACER.sleeps,
        }
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    encoded = "\n".join(json.dumps(row, sort_keys=True, separators=(",", ":")) for row in rows) + "\n"
    OUT.write_text(encoded)
    REPORT.write_text(render_report(rows))

    replayed = [json.loads(line) for line in OUT.read_text().splitlines()]
    assert replayed == rows
    assert REPORT.read_text().startswith("# L23 multivariable prime-complexity audit\n")
    summary = rows[-1]
    print(
        "L23 multivar: "
        f"{len(rows)} rows; smallest certified moving-a squarefree degree "
        f"{summary['smallest_certified_squarefree_representative_degree']}; "
        "tested-family obstruction PROVED; residual Capell stratum OPEN"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
