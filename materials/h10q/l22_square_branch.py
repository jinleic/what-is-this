#!/usr/bin/env python3
"""Uniform alpha-square branch, irreducibility, and local-alignment replay.

Replay from the workspace repository root:

    nice -n 19 python3 math/h10q/l22_square_branch.py

The algebraic and generalized-class theorem is symbolic.  It is not eligible
for the existing six-witness chain because its HIT parameter ranges through
the infinite L11c branch family.  Finite rows are labelled EVIDENCE.
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h10q  # noqa: E402 (the repository's exact-arithmetic authority)


OUT = HERE / "data" / "l22_square_branch.jsonl"
REPORT = Path("/tmp/l22_square_branch.md")
PACE_EVERY = 200
PACE_SECONDS = 0.005
IRRED_PRIME_LIMIT = 1000
SAMPLE_CELLS = (
    (3, F(6)),
    (5, F(5)),
    (7, F(7, 2)),
    (11, F(-11, 3)),
    (59, F(-177, 7)),
    (89, F(89, 8)),
    (179, F(358)),
)
CHAIN_OBSTRUCTION = (
    "fixed-finite-branch compatibility: prove that a predeclared finite L11c "
    "branch menu contains a locally aligned irreducible specialization for "
    "every cell, without introducing the HIT parameter as a seventh witness"
)
L11C_CITATION = (
    "THEOREMS.md, L11c: a finite square-branch subfamily may be adjoined; "
    "an infinite one may not because its parameter would be an extra witness"
)


class Pacer:
    """Keep finite searches below the requested duty cycle."""

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


def place_text(place: int | str) -> str:
    return "inf" if place == h10q.OO else str(place)


def gcd_many(values: list[int]) -> int:
    answer = 0
    for value in values:
        answer = math.gcd(answer, abs(value))
    return answer


# A tiny exact sparse Laurent-polynomial ring.  Negative exponents are used
# only for A^{-1} and lambda^{-1}; this is enough to replay every displayed
# branch identity without a CAS dependency.
VARIABLES = ("a", "A", "Z", "D", "r", "s", "b", "lambda")
NVAR = len(VARIABLES)
Monomial = tuple[int, ...]
Laurent = dict[Monomial, F]
ZERO_MONOMIAL = (0,) * NVAR


def lp_clean(poly: Laurent) -> Laurent:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def lp_const(value: F | int) -> Laurent:
    value = F(value)
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def lp_monomial(coefficient: F | int = 1, **powers: int) -> Laurent:
    exponent = [0] * NVAR
    for variable, power in powers.items():
        exponent[VARIABLES.index(variable)] = power
    coefficient = F(coefficient)
    return {} if coefficient == 0 else {tuple(exponent): coefficient}


def lp_add(*polys: Laurent) -> Laurent:
    answer: Laurent = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            answer[monomial] = answer.get(monomial, F(0)) + coefficient
    return lp_clean(answer)


def lp_neg(poly: Laurent) -> Laurent:
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def lp_sub(left: Laurent, right: Laurent) -> Laurent:
    return lp_add(left, lp_neg(right))


def lp_scale(poly: Laurent, scalar: F | int) -> Laurent:
    scalar = F(scalar)
    return lp_clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def lp_mul(*polys: Laurent) -> Laurent:
    answer = lp_const(1)
    for poly in polys:
        product: Laurent = {}
        for left_monomial, left_coefficient in answer.items():
            for right_monomial, right_coefficient in poly.items():
                monomial = tuple(
                    left_monomial[index] + right_monomial[index]
                    for index in range(NVAR)
                )
                product[monomial] = product.get(monomial, F(0)) + (
                    left_coefficient * right_coefficient
                )
        answer = lp_clean(product)
    return answer


def lp_pow(poly: Laurent, exponent: int) -> Laurent:
    assert exponent >= 0
    answer = lp_const(1)
    base = poly
    power = exponent
    while power:
        if power & 1:
            answer = lp_mul(answer, base)
        base = lp_mul(base, base)
        power >>= 1
    return answer


def lp_replace_monomial(poly: Laurent, source: str, target: str | None) -> Laurent:
    """Replace source by target, or by 1 when target is None."""
    source_index = VARIABLES.index(source)
    target_index = VARIABLES.index(target) if target is not None else None
    answer: Laurent = {}
    for monomial, coefficient in poly.items():
        exponent = list(monomial)
        moved = exponent[source_index]
        exponent[source_index] = 0
        if target_index is not None:
            exponent[target_index] += moved
        key = tuple(exponent)
        answer[key] = answer.get(key, F(0)) + coefficient
    return lp_clean(answer)


def symbolic_replay() -> dict[str, Any]:
    a = lp_monomial(a=1)
    A = lp_monomial(A=1)
    Z = lp_monomial(Z=1)
    D = lp_monomial(D=1)
    r = lp_monomial(r=1)
    s = lp_monomial(s=1)
    b = lp_monomial(b=1)
    lam = lp_monomial(**{"lambda": 1})
    A_inverse = lp_monomial(A=-1)
    lam_inverse = lp_monomial(**{"lambda": -1})
    one = lp_const(1)

    X_parameter = lp_scale(lp_add(lam, lp_mul(A, lam_inverse)), F(1, 2))
    r_parameter = lp_scale(lp_sub(lam, lp_mul(A, lam_inverse)), F(1, 2))
    tau_parameter = lp_mul(X_parameter, A_inverse)
    delta_parameter = lp_sub(one, lp_mul(A, lp_pow(tau_parameter, 2)))
    assert not lp_sub(lp_sub(lp_pow(X_parameter, 2), lp_pow(r_parameter, 2)), A)
    assert not lp_add(delta_parameter, lp_mul(lp_pow(r_parameter, 2), A_inverse))
    assert not lp_sub(lp_neg(lp_mul(A, delta_parameter)), lp_pow(r_parameter, 2))

    X_at_A = lp_replace_monomial(X_parameter, "lambda", "A")
    r_at_A = lp_replace_monomial(r_parameter, "lambda", "A")
    tau_at_A = lp_replace_monomial(tau_parameter, "lambda", "A")
    assert not lp_sub(lp_scale(X_at_A, 2), lp_add(A, one))
    assert not lp_sub(lp_scale(r_at_A, 2), lp_sub(A, one))
    assert not lp_sub(lp_scale(lp_mul(A, tau_at_A), 2), lp_add(A, one))

    X_at_one = lp_replace_monomial(X_parameter, "lambda", None)
    r_at_one = lp_replace_monomial(r_parameter, "lambda", None)
    tau_at_one = lp_replace_monomial(tau_parameter, "lambda", None)
    assert not lp_sub(lp_scale(X_at_one, 2), lp_add(A, one))
    assert not lp_add(lp_scale(r_at_one, 2), lp_sub(A, one))
    assert not lp_sub(lp_scale(lp_mul(A, tau_at_one), 2), lp_add(A, one))

    b_minus_one = lp_sub(b, one)
    Ng = lp_sub(
        lp_scale(lp_mul(lp_pow(a, 4), lp_pow(b, 2)), 16),
        lp_mul(A, lp_pow(b_minus_one, 4)),
    )
    E = lp_sub(one, lp_scale(lp_mul(A, lp_pow(s, 2), b), 2))
    P = lp_add(
        lp_scale(lp_mul(lp_pow(D, 2), lp_pow(A, 2), lp_pow(b, 4)), 16),
        lp_mul(
            lp_pow(r, 2),
            A_inverse,
            lp_pow(a, 4),
            lp_pow(Z, 4),
            lp_pow(Ng, 2),
        ),
        lp_scale(lp_mul(lp_pow(A, 3), lp_pow(s, 2), lp_pow(D, 2), lp_pow(b, 5)), -32),
    )
    square_one = lp_pow(lp_mul(r, lp_pow(a, 2), lp_pow(Z, 2), Ng), 2)
    square_two = lp_pow(lp_scale(lp_mul(D, A, lp_pow(b, 2)), 4), 2)
    norm_rhs = lp_add(square_one, lp_mul(A, square_two, E))
    assert not lp_sub(lp_mul(A, P), norm_rhs)

    B = lp_mul(A, square_two, E)
    C = lp_pow(lp_mul(lp_pow(a, 2), lp_pow(Z, 2), Ng), 2)
    assert not lp_sub(lp_mul(A, P), lp_add(B, lp_mul(lp_pow(r, 2), C)))

    return {
        "kind": "symbolic-replay",
        "label": "PROVED",
        "parameter_conic_identity": True,
        "delta_identity": True,
        "alpha_identity": True,
        "canonical_lambda_A": {
            "X": "(A+1)/2=1+2*a^2",
            "r": "(A-1)/2=2*a^2",
            "tau": "(A+1)/(2*A)=(1+2*a^2)/A",
        },
        "canonical_lambda_1": {
            "X": "(A+1)/2=1+2*a^2",
            "r": "-(A-1)/2=-2*a^2",
            "tau": "(A+1)/(2*A)=(1+2*a^2)/A",
        },
        "norm_identity": True,
        "pencil_identity": True,
        "norm_lhs_terms": len(lp_mul(A, P)),
        "norm_rhs_terms": len(norm_rhs),
    }


Polynomial = list[F]


def poly_trim(poly: Polynomial) -> Polynomial:
    answer = [F(value) for value in poly]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [F(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return poly_trim(answer)


def poly_scale(poly: Polynomial, scalar: F | int) -> Polynomial:
    return poly_trim([F(scalar) * value for value in poly])


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return poly_trim(answer)


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    answer = [F(1)]
    base = poly
    power = exponent
    while power:
        if power & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        power >>= 1
    return answer


def poly_eval(poly: Polynomial | list[int], value: F | int) -> F:
    answer = F(0)
    value = F(value)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def poly_compose_affine(poly: Polynomial, origin: int, step: int) -> Polynomial:
    answer = [F(0)]
    affine = [F(origin), F(step)]
    power = [F(1)]
    for coefficient in poly:
        answer = poly_add(answer, poly_scale(power, coefficient))
        power = poly_mul(power, affine)
    return poly_trim(answer)


def primitive_integer_poly(poly: Polynomial) -> tuple[list[int], F]:
    denominator = 1
    for coefficient in poly:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in poly]
    content = gcd_many(integers)
    assert content > 0
    integers = [value // content for value in integers]
    scale = F(content, denominator)
    if integers[-1] < 0:
        integers = [-value for value in integers]
        scale = -scale
    assert all(F(scale) * integers[index] == poly[index] for index in range(len(poly)))
    return integers, scale


def polynomial_data(a: int, z: F, tau: F) -> dict[str, Any]:
    A = F(1 + 4 * a * a)
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = F(a - 1, 2)
    delta = 1 - A * tau * tau
    alpha = -delta * A
    P = h10q._l10_P(F(a), Z, D, A, delta, s)
    return {
        "A": A,
        "Z": Z,
        "D": D,
        "s": s,
        "delta": delta,
        "alpha": alpha,
        "P": P,
    }


def ng_polynomial(a: int, A: int) -> Polynomial:
    return poly_add(
        poly_scale(poly_pow([F(-1), F(1)], 4), -A),
        [F(0), F(0), F(16 * a**4)],
    )


def replay_numeric_norm(a: int, z: F, lam: F) -> tuple[dict[str, Any], dict[str, Any]]:
    A = F(1 + 4 * a * a)
    X = (lam + A / lam) / 2
    r = (lam - A / lam) / 2
    tau = X / A
    data = polynomial_data(a, z, tau)
    Z, D, s, P = data["Z"], data["D"], data["s"], data["P"]
    assert data["delta"] == -(r * r) / A
    assert data["alpha"] == r * r
    assert X * X - r * r == A

    Ng = ng_polynomial(a, int(A))
    first_square = poly_pow(poly_scale(Ng, r * a * a * Z * Z), 2)
    second_square = [F(0), F(0), 4 * D * A]
    E = [F(1), -2 * A * s * s]
    rhs = poly_add(first_square, poly_mul(poly_scale(poly_pow(second_square, 2), A), E))
    assert poly_scale(P, A) == rhs

    return data, {
        "X": X,
        "r": r,
        "tau": tau,
        "norm_identity": True,
        "conic_identity": True,
        "delta_identity": True,
        "alpha_identity": True,
    }


def quadratic_character(value: F | int, prime: int) -> int:
    value = F(value)
    if h10q.vp(value, prime) != 0:
        return 0
    residue = h10q.unit_mod(value, prime)
    symbol = pow(residue, (prime - 1) // 2, prime)
    return 1 if symbol == 1 else -1


def choose_odd_a(prime: int) -> tuple[int, dict[str, int]]:
    assert prime % 2 == 1 and h10q._is_prime(prime)
    character_sum = 0
    nonsquare_residues = 0
    zero_residues = 0
    for residue in range(prime):
        PACER.tick()
        value = (1 + 4 * residue * residue) % prime
        if value == 0:
            character = 0
            zero_residues += 1
        else:
            character = 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1
            nonsquare_residues += character == -1
        character_sum += character
    assert character_sum == -1
    assert nonsquare_residues > 0
    for a in range(1, 4 * prime, 2):
        PACER.tick()
        if quadratic_character(1 + 4 * a * a, prime) == -1:
            return a, {
                "character_sum": character_sum,
                "nonsquare_residues": nonsquare_residues,
                "zero_residues": zero_residues,
            }
    raise AssertionError((prime, "no odd lift of a nonsquare residue"))


def irreducibility_certificate(P: Polynomial) -> tuple[int, list[int]]:
    primitive, _ = primitive_integer_poly(P)
    assert len(primitive) == 9
    tried: list[int] = []
    for prime in h10q.primerange(2, IRRED_PRIME_LIMIT + 1):
        PACER.tick()
        tried.append(prime)
        if h10q._l13_irred8(primitive, prime):
            return prime, tried
    raise AssertionError("no degree-8 Frobenius certificate <= 1000")


def choose_lambda(a: int, z: F, w: int) -> tuple[int, F, dict[str, Any], dict[str, Any], int, int]:
    A = 1 + 4 * a * a
    candidates = [value for n in range(1, 25) for value in (n, -n)]
    certificates_tried = 0
    for k in candidates:
        PACER.tick()
        lam = F(A + 8 * w * k)
        data, identities = replay_numeric_norm(a, z, lam)
        try:
            certificate, tried = irreducibility_certificate(data["P"])
        except AssertionError:
            certificates_tried += len(tuple(h10q.primerange(2, IRRED_PRIME_LIMIT + 1)))
            continue
        certificates_tried += len(tried)
        return k, lam, data, identities, certificate, certificates_tried
    raise AssertionError((w, z, "no finite specialization certificate in search box"))


def finite_support(*values: F | int) -> set[int]:
    support: set[int] = set()
    for value in values:
        support.update(h10q._l10_supp(F(value)))
    return support


def crt_pair(left: int, left_modulus: int, right: int, right_modulus: int) -> tuple[int, int]:
    assert math.gcd(left_modulus, right_modulus) == 1
    multiplier = ((right - left) * pow(left_modulus, -1, right_modulus)) % right_modulus
    modulus = left_modulus * right_modulus
    return (left + left_modulus * multiplier) % modulus, modulus


def q1_progression(S: set[int], w: int) -> tuple[int, int]:
    residue, modulus = 1, 8
    for prime in sorted(S - {2}):
        target = 1 if prime == w else pow(2 * w, -1, prime)
        residue, modulus = crt_pair(residue, modulus, target, prime)
    assert residue % 8 == 1
    assert math.gcd(residue, modulus) == 1
    return residue, modulus


def local_quantities(a: int, z: F, r: F, b: int) -> tuple[F, F, F, F]:
    A = F(1 + 4 * a * a)
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = F(a - 1, 2)
    delta = -(r * r) / A
    alpha = r * r
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    c = F(a * a) * Z * Z * Ng / (A * b * b * D)
    M = 16 - delta * c * c - 32 * A * b * s * s
    P_value = poly_eval(h10q._l10_P(F(a), Z, D, A, delta, s), b)
    assert M == P_value / (b**4 * D * D * A * A)
    x = alpha * M
    d = 2 * alpha * b
    return x, d, M, P_value


def choose_q1(
    a: int,
    z: F,
    w: int,
    r: F,
    P: Polynomial,
    S: set[int],
) -> tuple[int, int, int, dict[str, int]]:
    residue, modulus = q1_progression(S, w)
    primality_attempts = 0
    for index in range(0, 20000):
        PACER.tick()
        q1 = residue + index * modulus
        if q1 < 3:
            continue
        primality_attempts += 1
        try:
            is_prime = h10q._is_prime(q1)
        except h10q.PrimalityBound:
            continue
        if not is_prime:
            continue
        b0 = w * q1
        if poly_eval(P, b0) == 0:
            continue
        assert all(
            quadratic_character(2 * w * q1, prime) == 1
            for prime in S - {2, w}
        )
        A = 1 + 4 * a * a
        assert quadratic_character(A, q1) == 1
        x, d, _, _ = local_quantities(a, z, r, b0)
        symbols = {place_text(place): h10q.hilbert(x, d, place) for place in sorted(S)}
        symbols[str(q1)] = h10q.hilbert(x, d, q1)
        symbols["inf"] = h10q.hilbert(x, d, h10q.OO)
        assert all(symbol == 1 for symbol in symbols.values())
        return q1, residue, modulus, symbols
    raise AssertionError((w, z, modulus, "no proven prime in CRT progression search"))


def class_modulus(P: Polynomial, b0: int, A: int, S: set[int]) -> tuple[int, dict[int, int]]:
    exponents: dict[int, int] = {}
    for prime in sorted(S):
        PACER.tick()
        exponents[prime] = h10q._l10_exponent(P, F(b0), prime)
    modulus = math.lcm(8, 4 * A, *(prime ** exponent for prime, exponent in exponents.items()))
    return modulus, exponents


def expected_content(P_value: F, S: set[int]) -> F:
    answer = F(1)
    for prime in S:
        answer *= F(prime) ** h10q.vp(P_value, prime)
    return answer


def hilbert_two_odd_unit(unit_mod_8: int, odd_part_mod_8: int) -> int:
    """Exact (u, 2*v)_2 formula for odd u and v, using residue lifts."""
    assert unit_mod_8 in (1, 3, 5, 7)
    assert odd_part_mod_8 in (1, 3, 5, 7)
    exponent = (
        ((unit_mod_8 - 1) // 2) * ((odd_part_mod_8 - 1) // 2)
        + (unit_mod_8 * unit_mod_8 - 1) // 8
    )
    return -1 if exponent % 2 else 1


def two_adic_row(a: int, b: int, M: F) -> dict[str, Any]:
    A = 1 + 4 * a * a
    s = (a - 1) // 2
    U = M / 16
    assert h10q.vp(U, 2) == 0
    residue = h10q.unit_mod(U, 8)
    expected = 1 if s % 2 == 0 else (1 - 2 * (A % 8) * (b % 8)) % 8
    assert residue == expected
    table = {
        u: {v: hilbert_two_odd_unit(u, v) for v in (1, 3, 5, 7)}
        for u in (1, 3, 5, 7)
    }
    assert table[residue][b % 8] == 1
    return {
        "s_parity": s % 2,
        "b_mod_8": b % 8,
        "M_over_16_mod_8": residue,
        "hilbert_symbol": 1,
    }


def replay_sample(w: int, z: F) -> tuple[dict[str, Any], dict[str, Any]]:
    assert h10q.vp(z, w) >= 1
    a, character_count = choose_odd_a(w)
    assert a % 2 == 1
    A = 1 + 4 * a * a
    assert quadratic_character(A, w) == -1
    k, lam, data, identities, certificate, certificate_attempts = choose_lambda(a, z, w)
    r = identities["r"]
    assert h10q.vp(r, 2) == 1
    assert h10q.vp(r, w) == 0
    assert data["alpha"] == r * r
    assert data["delta"] == -(r * r) / A

    S = {2, 3, 5, 7, w} | finite_support(
        a,
        A,
        r,
        data["delta"],
        data["s"],
        z,
        data["D"],
    )
    q1, q_residue, q_modulus, symbols = choose_q1(a, z, w, r, data["P"], S)
    b0 = w * q1
    x0, d0, M0, P0 = local_quantities(a, z, r, b0)
    assert all(h10q.hilbert(x0, d0, place) == 1 for place in S)
    assert h10q.hilbert(x0, d0, q1) == 1
    assert h10q.hilbert(x0, d0, h10q.OO) == 1

    m = h10q.vp(data["Z"], w)
    assert m >= 3
    selected_valuations = {
        "b4_term": 4,
        "r2_Z4_Ng2_term": 2 * h10q.vp(r, w) + 4 * m,
        "b5_term_lower_bound": None if data["s"] == 0 else 5 + 2 * h10q.vp(data["s"], w),
    }
    assert selected_valuations["r2_Z4_Ng2_term"] >= 12
    if selected_valuations["b5_term_lower_bound"] is not None:
        assert selected_valuations["b5_term_lower_bound"] >= 5
    assert h10q.vp(P0, w) == 4
    assert h10q.vp(M0, w) == 0
    assert quadratic_character(M0, w) == 1

    two_row = two_adic_row(a, b0, M0)
    assert h10q.hilbert(x0, d0, 2) == 1

    N, exponents = class_modulus(data["P"], b0, A, S)
    assert math.gcd(q1, N) == 1
    assert all(N % prime == 0 for prime in S)
    b1 = w * (q1 + N)
    x1, d1, _, _ = local_quantities(a, z, r, b1)
    frozen_replay = {
        place_text(place): h10q.hilbert(x1, d1, place) for place in sorted(S)
    }
    assert all(symbol == 1 for symbol in frozen_replay.values())

    F_poly = poly_compose_affine(data["P"], b0, w * N)
    G, content = primitive_integer_poly(F_poly)
    assert len(G) == 9
    assert G[-1] > 0
    assert content == expected_content(P0, S)
    assert finite_support(content) <= S
    product_values = [(q1 + N * t) * int(poly_eval(G, t)) for t in range(9)]
    fixed_divisor_gcd = gcd_many(product_values)
    assert fixed_divisor_gcd == 1
    assert G[0] != 0

    lambda_balanced = F(1, w ** (2 * m))
    _, balanced_identities = replay_numeric_norm(a, z, lambda_balanced)
    r_balanced = balanced_identities["r"]
    assert h10q.vp(r_balanced, w) == -2 * m
    xb, db, _, Pb = local_quantities(a, z, r_balanced, b0)
    balanced_valuations = {
        "b4_term": 4,
        "r2_Z4_Ng2_term": 2 * h10q.vp(r_balanced, w) + 4 * m,
        "b5_term_lower_bound": None if data["s"] == 0 else 5 + 2 * h10q.vp(data["s"], w),
    }
    assert balanced_valuations["r2_Z4_Ng2_term"] == 0
    assert h10q.vp(Pb, w) == 0
    assert h10q.vp(xb, w) % 2 == 0
    assert h10q.vp(db, w) % 2 == 1
    assert quadratic_character(h10q.unit_part(xb, w), w) == -1
    balanced_symbol = h10q.hilbert(xb, db, w)
    assert balanced_symbol == -1

    sample = {
        "kind": "sample-class",
        "label": "EVIDENCE",
        "w": w,
        "z": frac_text(z),
        "v_w_z": h10q.vp(z, w),
        "a": a,
        "A": A,
        "a_character_count": character_count,
        "lambda_progression": "A+8*w*k",
        "k": k,
        "lambda": frac_text(lam),
        "r": frac_text(r),
        "tau": frac_text(identities["tau"]),
        "delta": frac_text(data["delta"]),
        "alpha": frac_text(data["alpha"]),
        "D": frac_text(data["D"]),
        "identity_checks": {
            "conic": identities["conic_identity"],
            "delta": identities["delta_identity"],
            "alpha": identities["alpha_identity"],
            "norm": identities["norm_identity"],
        },
        "irreducibility_certificate_prime": certificate,
        "certificate_primes_tried": certificate_attempts,
        "S": sorted(S),
        "q1": q1,
        "q1_crt_residue": q_residue,
        "q1_crt_modulus": q_modulus,
        "base_symbols": symbols,
        "selected_w_valuations": selected_valuations,
        "two_adic_replay": two_row,
        "N": N,
        "taylor_exponents": {str(prime): exponent for prime, exponent in exponents.items()},
        "frozen_symbols_at_t_1": frozen_replay,
        "normalized_degree": len(G) - 1,
        "normalized_leading_coefficient": G[-1],
        "content": frac_text(content),
        "content_support": sorted(finite_support(content)),
        "fixed_divisor_gcd": fixed_divisor_gcd,
    }
    obstruction = {
        "kind": "balanced-w-obstruction",
        "label": "PROVED",
        "w": w,
        "z": frac_text(z),
        "m_v_w_Z": m,
        "lambda": frac_text(lambda_balanced),
        "v_w_r": h10q.vp(r_balanced, w),
        "term_valuations": balanced_valuations,
        "v_w_x_even": h10q.vp(xb, w),
        "v_w_d_odd": h10q.vp(db, w),
        "x_unit_character": quadratic_character(h10q.unit_part(xb, w), w),
        "hilbert_symbol": balanced_symbol,
    }
    return sample, obstruction


def theorem_rows(symbolic: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "kind": "parameterization-theorem",
            "label": "PROVED",
            "statement": (
                "Every rational point X^2-r^2=A is obtained uniquely from "
                "lambda=X+r!=0 by X=(lambda+A/lambda)/2 and "
                "r=(lambda-A/lambda)/2; tau=X/A, delta=-r^2/A, alpha=r^2."
            ),
            "symbolic_replay": symbolic["parameter_conic_identity"],
        },
        {
            "kind": "irreducibility-theorem",
            "label": "PROVED",
            "statement": (
                "For every nonzero integer a and nonzero rational Z with "
                "D=1-Z-a^2*Z^2 nonzero, the degree-8 square-branch P is "
                "irreducible in Q(lambda)[b]."
            ),
            "pencil": "A*P=B+rho*C, gcd(B,C)=1",
            "gcd_exceptions": "exactly a=0, Z=0, or D=0; none can occur in a cell",
            "first_pullback": "rho=r^2; -B/C=-A*(square)*L is nonsquare",
            "hyperbola_pullback": (
                "Q(lambda)=Q(r,sqrt(r^2+A)); neither E+A nor (E+A)/E "
                "is square, so no 4x4 splitting"
            ),
            "first_unproved_condition": None,
        },
        {
            "kind": "uniform-free-parameter-class-theorem",
            "label": "PROVED",
            "statement": (
                "For every odd prime w and rational z with v_w(z)>=1, an odd "
                "a, a HIT specialization lambda=A+8*w*k, a prime q1, and a "
                "modulus N produce an aligned class b=w*(q1+N*t) whose "
                "primitive normalized degree-8 G(t) is irreducible and whose "
                "Schinzel pair (q1+N*t,G(t)) has no fixed prime divisor."
            ),
            "scope": "algebraic/class existence with a cell-dependent free branch parameter",
            "six_count_applicability": "NOT APPLICABLE",
            "selection_order": ["a", "lambda by HIT", "finite S", "q1 by CRT/Dirichlet", "N"],
            "l19_l20_interface": {
                "branch_equations": "alpha=r^2 and delta=-r^2/A",
                "support_replacement": "S={2,3,5,7,w} union supp(a,A,r,delta,s,z,D)",
                "unchanged_after_replacement": [
                    "CRT/Dirichlet q1",
                    "Taylor square-class freezing",
                    "literal S-supported primitive content",
                    "gcd(q1,N)=1",
                    "S-unit normalization",
                    "fixed-prime-divisor proof",
                ],
            },
            "first_unproved_condition": CHAIN_OBSTRUCTION,
            "schinzel_prime_values": "CONDITIONAL for this off-chain selected class only",
        },
        {
            "kind": "six-count-chain-interface",
            "label": "OPEN",
            "verdict": "NOT APPLICABLE",
            "citation": L11C_CITATION,
            "reason": (
                "HIT selects k, hence lambda and tau, after seeing the cell. "
                "These values range through the infinite L11c family and k "
                "would be an additional witness."
            ),
            "first_unproved_condition": CHAIN_OBSTRUCTION,
            "conditional_record_updated": False,
        },
        {
            "kind": "balanced-route-theorem",
            "label": "PROVED",
            "statement": (
                "The tempting choice v_w(r)=-2*v_w(Z), used to expose a full "
                "degree reduction at w, forces (x,d)_w=(A|w)=-1 when "
                "v_w(b)=1 and therefore cannot be an aligned class."
            ),
        },
    ]


def audit_rows() -> list[dict[str, Any]]:
    entries = (
        (
            "choose a",
            "none",
            "a is an odd lift with (A|w)=-1; sum_a (1+4a^2|w)=-1",
            "PROVED",
        ),
        (
            "lambda specialization",
            "canonical value replaced",
            "HIT inside lambda=A+8*w*k; v_2(r)=1 and v_w(r)=0 for every k",
            "PROVED",
        ),
        (
            "w-symbol",
            "only its 2- and w-adic neighborhoods",
            "term valuations are (4, >=12, >=5), hence M is a w-adic square and x=r^2*M is square",
            "PROVED",
        ),
        (
            "2-symbol",
            "only v_2(r)=1",
            "a,b odd; v_2(c)>=4 and the exact U=M/16 residue table gives (U,2b)_2=+1",
            "PROVED",
        ),
        (
            "other frozen symbols",
            "none",
            "after lambda is fixed, put every fixed-data prime in S and impose 2*w*q1 square modulo p; then d is a p-adic square",
            "PROVED",
        ),
        (
            "moving q1/Q symbol",
            "none",
            "alpha=r^2 and delta=-r^2/A give squareclass(x)=A; reciprocity plus the CRT conditions gives (A|Q)=+1",
            "PROVED",
        ),
        (
            "infinite symbol",
            "none",
            "d=2*r^2*b is positive on positive prime values",
            "PROVED",
        ),
        (
            "Taylor freezing",
            "none",
            "h10q._l10_exponent is applied after lambda and S are frozen",
            "PROVED",
        ),
        (
            "content / normalization",
            "canonical r formerly hid support(a)",
            "the enlarged S contains supports of a,A,r,delta,s,z,D; the primitive content is S-supported",
            "PROVED",
        ),
        (
            "irreducibility",
            "canonical lambda was the old gap",
            "irreducibility over Q(lambda), followed by HIT and affine substitution",
            "PROVED",
        ),
        (
            "no fixed prime divisor",
            "none",
            "S contains 2,3,5,7; outside S a nonzero degree-8 polynomial cannot vanish on every nonzero residue",
            "PROVED",
        ),
        (
            "six-count branch interface",
            "tau_dagger is one predeclared branch",
            "cell-dependent HIT ranges through the infinite L11c family and costs an extra witness; a fixed finite menu theorem is absent",
            "OPEN",
        ),
    )
    return [
        {
            "kind": "local-audit",
            "label": label,
            "step": step,
            "canonical_r_dependency": dependency,
            "generalized_input": generalized,
        }
        for step, dependency, generalized, label in entries
    ]


def write_report(
    symbolic: dict[str, Any],
    theorem: list[dict[str, Any]],
    audits: list[dict[str, Any]],
    samples: list[dict[str, Any]],
    obstructions: list[dict[str, Any]],
    elapsed: float,
) -> None:
    lines = [
        "# L22 - the free alpha-square branch",
        "",
        "**Overall verdict - PROVED algebraic/class theorem; NOT APPLICABLE to the six-count.**  For every cell `(w,z)`, the free parameter produces an aligned admissible class whose primitive normalized octic is irreducible.  However, this does **not** close the last chain premise and does **not** update the conditional record: the HIT-selected `lambda` (equivalently the L11c branch parameter) depends on the cell and ranges through an infinite branch family.",
        "",
        "The chain-interface verdict is **OPEN**.  `THEOREMS.md` L11c explicitly permits adjoining a finite square-branch subfamily but forbids an infinite one because its parameter would require another witness.  The first unproved condition is fixed-finite-branch compatibility: a predeclared finite menu must be proved to contain a locally aligned irreducible specialization for every cell without adding a seventh witness.",
        "",
        "Finite rows in the JSON artifact are **EVIDENCE** only.  The off-chain uniform class theorem comes from the symbolic pencil/base-change proof, Hilbert irreducibility, quadratic reciprocity, and the exact local calculations below—not from the scan.",
        "",
        "## 1. All rational alpha-square points",
        "",
        "Put `A=1+4a^2` and `X=A*tau`.  The condition that `alpha=-A*(1-A*tau^2)=X^2-A` be a rational square says",
        "",
        "    X^2-r^2=A.",
        "",
        "Every rational point is parameterized, with no omission, by `lambda=X+r != 0`:",
        "",
        "    X=(lambda+A/lambda)/2,    r=(lambda-A/lambda)/2,",
        "    tau=X/A,                  delta=-r^2/A,",
        "    alpha=r^2.",
        "",
        "Indeed `(X-r)(X+r)=A`, so the inverse is `lambda=X+r`.  Conversely the displayed formula gives `X^2-r^2=A`.  Since `A` lies strictly between `(2|a|)^2` and `(2|a|+1)^2`, it is not a rational square; hence rational `lambda != 0` never gives `r=0`.  At `lambda=A`, `r=2a^2`; at `lambda=1`, `r=-2a^2`; both give the old `tau_dagger=(1+2a^2)/A`.  The sparse Laurent replay checked all of these cross-identities exactly.",
        "",
        "## 2. Replayed norm identity",
        "",
        "Write",
        "",
        "    Z=z^3,  D=1-Z-a^2*Z^2,  s=(a-1)/2,",
        "    Ng=16*a^4*b^2-A*(b-1)^4,  L=1-2*A*s^2*b.",
        "",
        "On the free square branch, the repository polynomial is",
        "",
        "    P=16*D^2*A^2*b^4+(r^2/A)*a^4*Z^4*Ng^2-32*A^3*s^2*D^2*b^5.",
        "",
        "Exact expansion in the local sparse Laurent ring proves",
        "",
        "    A*P=(r*a^2*Z^2*Ng)^2 + A*(4*D*A*b^2)^2*(1-2*A*s^2*b).",
        "",
        f"Both sides normalize to the same {symbolic['norm_lhs_terms']} nonzero Laurent monomials.  Every finite sample independently replayed the same identity through `h10q._l10_P`.",
        "",
        "## 3. Irreducibility before and after both pullbacks",
        "",
        "Fix nonzero integral `a` and rational `Z != 0`; put `D=1-Z-a^2 Z^2 != 0`.  Set",
        "",
        "    B=16*A^3*D^2*b^4*L,       C=(a^2*Z^2*Ng)^2.",
        "",
        "Then `A*P=B+rho*C`.  The two coefficients are coprime in `Q[b]`.  Certainly `C(0) != 0`, so the factor `b` is harmless.  If `s=0`, then `L=1`.  If `s!=0` and a zero beta of `L` were also a zero of `Ng`, then `beta != 1` and",
        "",
        "    A=(4*a^2*beta/(beta-1)^2)^2,",
        "",
        "contradicting that the integer `A=1+4a^2` is not a rational square.  Thus `gcd(B,C)=1`.  Conversely, `a=0` or `Z=0` makes `C=0`, and `D=0` makes `B=0`; these are exactly the pencil-gcd exceptions.  None occurs in a cell: the chosen `a` is nonzero, `v_w(z)>=1` makes `Z` nonzero, and `D=1 mod w`.  A primitive polynomial of degree one in `rho` can factor only through a common divisor of its two coefficients, so `B+rho*C` is irreducible over `Q(rho)[b]`.",
        "",
        "For the first pullback `rho=r^2`, the cover can split only if `E=-B/C` is a square in `Q(b)`.  But",
        "",
        "    E=-A*L*(4*A*D*b^2/(a^2*Z^2*Ng))^2.",
        "",
        "When `s!=0`, the simple zero of `L` makes this nonsquare.  When `s=0` its square class is the negative rational `-A`, again nonsquare.  Hence `P` is irreducible over `Q(r)[b]`.",
        "",
        "The hyperbola parameter gives `Q(lambda)=Q(r,sqrt(r^2+A))`.  In the first cover's function field `r^2=E`, so a second split would require `E+A` to be a square in `Q(b,sqrt(E))`.  For an element of `Q(b)`, this happens exactly when `E+A` or `(E+A)/E` is a square in `Q(b)`.  At `b=0`, the first is regular with value `A`, a nonsquare.  The second has a pole of order four whose leading coefficient is minus a rational square, so it too is nonsquare.  This rules out the possible `4 x 4` factorization and proves `P` irreducible in `Q(lambda)[b]` for every fixed admissible `(a,z)`.",
        "",
        "## 4. HIT in an explicit local neighborhood",
        "",
        "Choose `a` odd with `(A|w)=-1`.  Such an odd lift always exists because the elementary character sum `sum_u (1+4u^2|w)=-1` has a nonsquare term, and parity can be changed by adding `w`.  Apply Hilbert irreducibility after the invertible change",
        "",
        "    lambda=A+8*w*k.",
        "",
        "The resulting polynomial remains irreducible over `Q(k)`, so infinitely many integral `k` give an irreducible specialized octic.  This whole progression lies in the needed local neighborhood.  Since `A-1=4a^2` has 2-adic valuation two and the perturbation is divisible by eight, `v_2(r)=1`.  Modulo `w`, `r=(A-1)/2=2a^2` is a unit, so `v_w(r)=0`.  Thus HIT is done first; only then is the finite support of the selected `r` frozen.",
        "",
        "**Six-count warning.**  HIT chooses `k` only after `(w,z)` is fixed, and the resulting `tau=(lambda+A/lambda)/(2A)` ranges through the one-parameter square-branch family of `THEOREMS.md` L11c.  L11c allows only a predeclared finite subfamily in the existing formula; an infinite family needs its parameter as another witness.  Thus this specialization theorem is mathematically uniform but is not an admissible replacement inside the six-witness record.",
        "",
        "## 5. Exact local alignment and class existence",
        "",
        "Let `Z=z^3`, so `m=v_w(Z)>=3`; also `D` and `A` are `w`-units.  For `b=w*q` with `q` a `w`-unit, the three terms of `P` have valuations",
        "",
        "    4,        2*v_w(r)+4*m >= 12,        at least 5.",
        "",
        "This is the load-bearing cell calculation: the `b^4` term uniquely dominates.  Therefore `M=P/(b^4 D^2 A^2)` is in `16*(1+w Z_w)`, `x=r^2 M` is a square, and `(x,2r^2b)_w=+1`.  It would be false without `v_w(z)>=1`.",
        "",
        "At two, `a,b` are odd, `A=5 mod 8`, and `v_2(r)=1`.  For `c=a^2 Z^2 Ng/(A b^2 D)`, one has `v_2(c)>=4` for every rational `z`: if `v_2(Z)<0`, the `a^2Z^2` term uniquely controls `D`; otherwise `D` is odd.  Put `U=M/16`.  Modulo eight,",
        "",
        "    U=1                         if s is even,",
        "    U=1-2*A*b                   if s is odd.",
        "",
        "The four odd residue classes of `b` give `(U,2b)_2=+1` in both cases.  This is the exact 2-adic table; no appeal to vague local constancy is needed.",
        "",
        "After selecting `lambda`, enlarge the controlled set to",
        "",
        "    S={2,3,5,7,w} union supp(a,A,r,delta,s,z,D).",
        "",
        "For every odd `p in S\\{w}`, impose `(2*w*q1|p)=+1`.  These independent nonzero residue conditions and `q1=1 mod 8` define a reduced CRT class, hence Dirichlet supplies arbitrarily large prime `q1`; discard the finitely many roots of `P(w*q1)`.  Then `d=2r^2wq1` is a square at every such `p`, so every frozen symbol is `+1`.",
        "",
        "At the moving prime `q1` (and later at a prime `Q=q1+Nt`), `v_Q(c)=-2`; consequently `x` has even valuation and square class `A`, so the symbol is `(A|Q)`.  Every prime divisor of `A` is `1 mod 4`.  The CRT equations give",
        "",
        "    (A|q1)=(q1|A)=(2*w|A)=(2|A)*(w|A)=(-1)*(-1)=+1.",
        "",
        "Taking `N` divisible by `8`, `4A`, and every Taylor exponent `p^k` for `p in S` freezes both the controlled symbols and this moving character.  At infinity `d>0`.  Thus every named symbol is aligned.",
        "",
        "Finally set `b(t)=w*(q1+N*t)` and primitively normalize `P(b(t))=c*G(t)`.  The enlarged support makes `c` literally `S`-supported.  Affine substitution preserves irreducibility, so `G` is an irreducible positive-leading octic.  The linear polynomial and `G` have no common root because at the linear root `b=0` and `P(0)=r^2*A*a^4*Z^4 != 0`.  There is no fixed prime divisor: primes in `S` see units, while outside `S` a nonzero degree-8 polynomial cannot vanish on every nonzero residue unless `p-1<=8`; all such primes are already in `{2,3,5,7}`.  This proves the generalized free-parameter class theorem, but not its eligibility for the fixed six-count branch menu.",
        "",
        "### Explicit interface with L19/L20",
        "",
        "This is a clean generalization, not a citation of the canonical formula `delta=-4a^4/A`.  Two support implications used implicitly on `tau_dagger` are no longer available: `alpha=4a^4` formerly put every prime of `a` into `S`, and `delta=-4a^4/A` formerly put every prime of `A` into `S`.  The new definition of `S` inserts `a` and `A` explicitly (and inserts `r,delta,s,z,D` so all fixed coefficients and moving-prime unit guards are explicit).",
        "",
        "With that replacement, the L19/L20 interface is unchanged for exact reasons.  The CRT prime `q1` is a unit at every prime of `S`; `N` is supported on `S` and is divisible by the Taylor exponents and `4A`, so `gcd(q1,N)=1`.  The Taylor lemma makes `P(b(t))/P(b(0))` and `b(t)/b(0)` local squares, giving coefficientwise `S`-unit normalization and frozen symbols.  Outside `S`, the octic leading coefficient `r^2*A*a^4*Z^4*(wN)^8` is a unit, so the primitive content has no outside prime and is literally `S`-supported.  The fixed-divisor proof then uses only primitivity, degree eight, and `{2,3,5,7} subset S`.  The moving-prime computation uses only `alpha=r^2` and `delta=-r^2/A`.  Thus no step in the generalized class construction invokes `r=2a^2` or `delta=-4a^4/A` verbatim.",
        "",
        "This exact L19/L20 extension therefore proves the class theorem **after a branch has been selected**.  It does not repair the earlier L11c interface: selecting that branch by a cell-dependent HIT parameter is the separate OPEN condition.",
        "",
        "## 6. Sharp audit of the former canonical dependencies",
        "",
        "| Step | Dependence on numerical `r=2a^2` | General replacement | Label |",
        "|---|---|---|---|",
    ]
    for audit in audits:
        lines.append(
            f"| {audit['step']} | {audit['canonical_r_dependency']} | {audit['generalized_input']} | **{audit['label']}** |"
        )
    lines.extend(
        [
            "",
            f"**First unproved condition - OPEN:** {CHAIN_OBSTRUCTION}.  The algebraic irreducibility and generalized L19/L20 class mechanics are proved; the obstacle is their ineligibility for the existing fixed six-count branch interface.",
            "",
            "## 7. Why the balanced w-adic reduction is not the proof",
            "",
            "If instead one forces `v_w(r)=-2m` to expose all octic coefficients modulo `w`, then the three term valuations become `(4,0,>=5)`.  The norm-square term uniquely dominates, `x` has even valuation and unit square class `A`, while `d` has odd valuation.  Hence",
            "",
            "    (x,d)_w=(A|w)=-1.",
            "",
            "So the attractive full-degree reduction at the alignment prime is exactly incompatible with the orthodox `v_w(b)=1` class.  This is a **PROVED obstruction** to that route, not evidence against the HIT-in-a-neighborhood route.",
            "",
            "## 8. Exact replay evidence",
            "",
            f"The script replayed {len(samples)} varied cells, including `w=3,5,7,89,179` and all three signs of `v_2(Z)`.  Every row has: a finite-field Frobenius certificate for its selected specialized octic; the norm identity; the selected `(4,>=12,>=5)` valuation pattern; the exact 2-adic table; all frozen, moving, and infinite symbols; Taylor freezing at `t=1`; primitive content support; and fixed-divisor gcd one.  Each balanced comparison returned symbol `-1`.  These rows are labelled **EVIDENCE** and are not used to infer the theorem.",
            "",
            f"Labels: branch parameterization **PROVED**; norm identity **PROVED**; generic and pulled-back irreducibility **PROVED**; generalized free-parameter class existence **PROVED**; balanced-route obstruction **PROVED**; finite-menu/six-count compatibility **OPEN**.  Six-count applicability: **NOT APPLICABLE**.  The conditional record is unchanged.",
            "",
            f"Verification wall-clock: **{elapsed:.3f} s**. Pacing: {PACER.steps} finite-loop steps, {PACER.sleeps} sleeps of {PACE_SECONDS} s.",
            "",
            "Artifacts: `math/h10q/l22_square_branch.py`, `math/h10q/data/l22_square_branch.jsonl`, `/tmp/l22_square_branch.md`.",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    symbolic = symbolic_replay()
    theorem = theorem_rows(symbolic)
    audits = audit_rows()
    samples: list[dict[str, Any]] = []
    obstructions: list[dict[str, Any]] = []
    for w, z in SAMPLE_CELLS:
        sample, obstruction = replay_sample(w, z)
        samples.append(sample)
        obstructions.append(obstruction)
    elapsed = time.perf_counter() - started
    summary = {
        "kind": "summary",
        "label": "OPEN",
        "algebraic_class_theorem": "PROVED",
        "irreducibility": "PROVED",
        "local_alignment": "PROVED",
        "normalized_irreducible_class": "PROVED",
        "six_count_applicability": "NOT APPLICABLE",
        "closes_last_chain_gap": False,
        "first_unproved_condition": CHAIN_OBSTRUCTION,
        "l11c_citation": L11C_CITATION,
        "conditional_record_updated": False,
        "downstream_schinzel_prime_values": "CONDITIONAL for the off-chain class only",
        "sample_rows": len(samples),
        "sample_rows_label": "EVIDENCE",
        "balanced_obstructions": len(obstructions),
        "refusals_used_as_evidence": False,
        "wall_seconds": round(elapsed, 6),
        "pacing": {"steps": PACER.steps, "sleeps": PACER.sleeps, "sleep_seconds": PACE_SECONDS},
    }
    rows = [symbolic, *theorem, *audits, *samples, *obstructions, summary]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(symbolic, theorem, audits, samples, obstructions, elapsed)
    print(
        "L22 square branch: "
        f"samples={len(samples)}, balanced-obstructions={len(obstructions)}, "
        f"chain-open=1, wall={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
