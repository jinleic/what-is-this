#!/usr/bin/env python3
"""L23 zero-quantifier-cost branch/norm absorption audit.

This is a standalone stdlib-only exact-arithmetic replay.  It proves the
self-coupled square-branch identities and their finite-fibre count, audits the
places where a direct identification changes the L6 architecture, and records
positive finite-field certificates against the untwisted Capell escape.
Finite scans are never used as negative evidence.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
from typing import Any
import json
import math
import time


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l23_absorption.jsonl"
REPORT = Path("/tmp/l23_absorption.md")
PACE_EVERY = 400
PACE_SECONDS = 0.002
CERTIFICATE_LIMIT = 251


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
Polynomial = list[F]


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def vp(value: F | int, prime: int) -> int:
    value = F(value)
    if value == 0:
        return 10**9
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


def rational_sqrt(value: F | int) -> F | None:
    value = F(value)
    if value < 0:
        return None
    numerator = math.isqrt(value.numerator)
    denominator = math.isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return F(numerator, denominator)
    return None


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
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return poly_trim(answer)


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    answer = [F(1)]
    base = poly
    while exponent:
        if exponent & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return answer


def poly_eval(poly: Polynomial | list[int], value: F | int) -> F:
    answer = F(0)
    value = F(value)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def poly_derivative(poly: Polynomial | list[int]) -> Polynomial:
    return [F(index) * F(poly[index]) for index in range(1, len(poly))] or [F(0)]


def primitive_integer_poly(poly: Polynomial) -> list[int]:
    denominator = 1
    for coefficient in poly:
        denominator = math.lcm(denominator, F(coefficient).denominator)
    integers = [int(F(coefficient) * denominator) for coefficient in poly]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, abs(coefficient))
    assert content > 0
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return integers


def ng_polynomial(a: F, A: F) -> Polynomial:
    return poly_add(
        poly_scale(poly_pow([F(-1), F(1)], 4), -A),
        [F(0), F(0), 16 * a**4],
    )


def square_branch_polynomial(a_value: F | int, Z_value: F | int, rho_value: F | int) -> Polynomial:
    a = F(a_value)
    Z = F(Z_value)
    rho = F(rho_value)
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    Ng = ng_polynomial(a, A)
    first = [F(0)] * 4 + [16 * D * D * A * A]
    second = poly_scale(poly_mul(Ng, Ng), rho * rho * a**4 * Z**4 / A)
    third = [F(0)] * 5 + [-32 * A**3 * s * s * D * D]
    return poly_add(poly_add(first, second), third)


def canonical_polynomial(a_value: F | int, Z_value: F | int) -> Polynomial:
    a = F(a_value)
    return square_branch_polynomial(a, Z_value, 2 * a * a)


def ng_value(a: F, A: F, b: F) -> F:
    return 16 * a**4 * b * b - A * (b - 1) ** 4


def bridge_value(a_value: F | int, Z_value: F | int, b_value: F | int) -> F:
    a = F(a_value)
    Z = F(Z_value)
    b = F(b_value)
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    assert b != 0 and D != 0
    return a * a * Z * Z * ng_value(a, A, b) / (A * b * b * D)


def identity_replay() -> dict[str, Any]:
    checks = 0
    for a_int in (1, 3, 5, 7):
        a = F(a_int)
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        for Z in (F(-8), F(-1, 8), F(1), F(27)):
            D = 1 - Z - a * a * Z * Z
            assert D != 0
            for lam in (F(1), F(2), F(3, 2), A):
                X = (lam + A / lam) / 2
                rho = (lam - A / lam) / 2
                assert X * X - rho * rho == A
                assert 1 - A * (X / A) ** 2 == -rho * rho / A
                P = square_branch_polynomial(a, Z, rho)
                assert len(P) == 9
                endpoint = rho * rho * A * a**4 * Z**4
                assert P[0] == endpoint and P[8] == endpoint
                for b in (F(-3), F(3), F(5, 3)):
                    c = bridge_value(a, Z, b)
                    assert vp(c, 2) >= 4
                    B = 2 * b
                    M = 16 + rho * rho * c * c / A - 16 * A * B * s * s
                    assert poly_eval(P, b) == D * D * A * A * b**4 * M

                    # Orientation I: branch (X,rho) is the tied (y,r).
                    y, r = X, rho
                    C1 = y * y - r * r - A
                    C2 = -r * r * (c * c - A * y * y) - 16 * A * B * (r * r - A * s * s) - 16 * A
                    u = r * r
                    Q1 = A * u * u + (A * A - c * c - 16 * A * B) * u + 16 * A * A * B * s * s - 16 * A
                    assert C1 == 0 and C2 == Q1

                    # Orientation II: swap the hyperbola coordinates relative to tied y,r.
                    y2, r2 = rho, X
                    C1_swap = r2 * r2 - y2 * y2 - A
                    C2_swap = -y2 * y2 * (c * c - A * y2 * y2) - 16 * A * B * (r2 * r2 - A * s * s) - 16 * A
                    u2 = y2 * y2
                    Q2 = A * u2 * u2 - (c * c + 16 * A * B) * u2 + 16 * A * A * B * (s * s - 1) - 16 * A
                    assert C1_swap == 0 and C2_swap == Q2
                    checks += 1

            P_can = canonical_polynomial(a, Z)
            canonical_endpoint = 4 * A * a**8 * Z**4
            assert P_can[0] == P_can[8] == canonical_endpoint
            X_can = 1 + 2 * a * a
            correction_I = poly_mul(
                [F(0)] * 4 + [4 * a**4 * D * D * A * A],
                [X_can * X_can, F(-32)],
            )
            J_I = poly_add(P_can, poly_scale(correction_I, -1))
            correction_II = poly_mul(
                [F(0)] * 4 + [D * D * A * A],
                [16 * a**8, -32 * X_can * X_can],
            )
            J_II = poly_add(P_can, poly_scale(correction_II, -1))
            assert len(J_I) == len(J_II) == 9
            assert J_I[0] == J_I[8] == canonical_endpoint
            assert J_II[0] == J_II[8] == canonical_endpoint
            checks += 1

    return {
        "type": "identity-replay",
        "label": "PROVED",
        "exact_instances": checks,
        "hyperbola_identity": "X^2-rho^2=A",
        "delta_identity": "delta=-rho^2/A",
        "bridge_identity": "P=D^2*A^2*b^4*M",
        "orientation_I_eliminant": "A*u^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A",
        "orientation_II_eliminant": "A*u^2-(c^2+16AB)u+16A^2B(s^2-1)-16A",
        "canonical_endpoints": "P_0=P_8=4*A*a^8*Z^4",
    }


# Finite-field polynomial arithmetic.  These routines certify individual rows;
# absence of a certificate is an assertion failure, never a negative conclusion.
def primes_up_to(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return primes


PRIMES = primes_up_to(CERTIFICATE_LIMIT)


def mod_trim(poly: list[int], prime: int) -> list[int]:
    answer = [coefficient % prime for coefficient in poly]
    while answer and answer[-1] == 0:
        answer.pop()
    return answer


def mod_reduce(poly: list[int], modulus: list[int], prime: int) -> list[int]:
    answer = mod_trim(poly, prime)
    modulus = mod_trim(modulus, prime)
    inverse = pow(modulus[-1], prime - 2, prime)
    while len(answer) >= len(modulus):
        scalar = answer[-1] * inverse % prime
        shift = len(answer) - len(modulus)
        for index, coefficient in enumerate(modulus):
            answer[index + shift] = (answer[index + shift] - scalar * coefficient) % prime
        answer = mod_trim(answer, prime)
    return answer


def mod_mul(left: list[int], right: list[int], modulus: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            product[i + j] = (product[i + j] + left_value * right_value) % prime
    return mod_reduce(product, modulus, prime)


def mod_pow(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    answer = [1]
    base = mod_reduce(base, modulus, prime)
    while exponent:
        if exponent & 1:
            answer = mod_mul(answer, base, modulus, prime)
        base = mod_mul(base, base, modulus, prime)
        exponent >>= 1
    return answer


def mod_xpk(modulus: list[int], prime: int, exponent: int) -> list[int]:
    answer = [0, 1]
    for _ in range(exponent):
        answer = mod_pow(answer, prime, modulus, prime)
    return answer


def mod_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left = mod_trim(left, prime)
    right = mod_trim(right, prime)
    while right:
        left, right = right, mod_reduce(left, right, prime)
    return left


def irreducible_mod(poly: list[int], prime: int) -> bool:
    modulus = mod_trim(poly, prime)
    degree = len(modulus) - 1
    if degree <= 0 or len(modulus) != len(poly):
        return False
    if mod_xpk(modulus, prime, degree) != [0, 1]:
        return False
    divisors: set[int] = set()
    remaining = degree
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            divisors.add(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        divisors.add(remaining)
    for divisor in divisors:
        power = mod_xpk(modulus, prime, degree // divisor)
        length = max(2, len(power))
        difference = [0] * length
        for index, coefficient in enumerate(power):
            difference[index] = coefficient
        difference[1] = (difference[1] - 1) % prime
        if len(mod_gcd(difference, modulus, prime)) != 1:
            return False
    return True


def eval_mod(poly: list[int], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(poly):
        answer = (answer * value + coefficient) % prime
    return answer


def legendre(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0:
        return 0
    power = pow(residue, (prime - 1) // 2, prime)
    return 1 if power == 1 else -1


def irreducibility_certificate(poly: Polynomial) -> int:
    primitive = primitive_integer_poly(poly)
    for prime in PRIMES:
        PACER.tick()
        if irreducible_mod(primitive, prime):
            return prime
    raise AssertionError("no irreducibility certificate in the declared finite window")


def bad_simple_root_certificate(poly: Polynomial) -> tuple[int, int]:
    primitive = primitive_integer_poly(poly)
    derivative = [index * primitive[index] for index in range(1, len(primitive))]
    for prime in PRIMES:
        if primitive[-1] % prime == 0:
            continue
        for root in range(1, prime):
            PACER.tick()
            if eval_mod(primitive, root, prime) == 0 and eval_mod(derivative, root, prime) != 0:
                if legendre(2 * root, prime) == -1:
                    return prime, root
    raise AssertionError("no positive bad-root certificate in the declared finite window")


def capell_certificate_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    Z_values = (F(-27), F(-8), F(-1), F(-1, 8), F(1, 8), F(1), F(8), F(27))
    for a in (1, 3, 5, 7, 9):
        for Z in Z_values:
            P = canonical_polynomial(a, Z)
            A = 1 + 4 * a * a
            D = 1 - Z - a * a * Z * Z
            assert D != 0 and len(P) == 9
            assert P[0] == P[8] == 4 * A * a**8 * Z**4
            irred_prime = irreducibility_certificate(P)
            bad_prime, root = bad_simple_root_certificate(P)
            primitive = primitive_integer_poly(P)
            derivative = [index * primitive[index] for index in range(1, len(primitive))]
            assert irreducible_mod(primitive, irred_prime)
            assert eval_mod(primitive, root, bad_prime) == 0
            assert eval_mod(derivative, root, bad_prime) != 0
            assert legendre(2 * root, bad_prime) == -1
            rows.append(
                {
                    "type": "capell-certificate",
                    "label": "PROVED",
                    "a": a,
                    "Z": frac_text(Z),
                    "v2_Z": vp(Z, 2),
                    "degree_P": 8,
                    "P_irreducible_mod_prime": irred_prime,
                    "bad_simple_root": {"prime": bad_prime, "root": root, "legendre_2root": -1},
                    "conclusion": "2*beta is not a square in Q[beta]/(P)",
                    "scope": "this exact (a,Z) only; the finite table is not negative evidence about untested parameters",
                }
            )
    return rows


def choose_target_a(prime: int) -> int:
    for residue in range(prime):
        A_residue = (1 + 4 * residue * residue) % prime
        if residue != 1 and legendre(A_residue, prime) == -1:
            lift = residue if residue % 2 else residue + prime
            assert lift % 2 == 1 and lift % prime != 1
            return lift
    raise AssertionError((prime, "no s-unit nonsquare-A residue"))


def target_newton_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    character_audits: list[dict[str, int]] = []
    for prime in PRIMES:
        if prime == 2:
            continue
        nonsquare = 0
        nonsquare_away_from_one = 0
        zeros = 0
        for residue in range(prime):
            character = legendre(1 + 4 * residue * residue, prime)
            nonsquare += character == -1
            zeros += character == 0
            nonsquare_away_from_one += character == -1 and residue != 1
        expected = (prime + 1) // 2 if prime % 4 == 3 else (prime - 1) // 2
        assert nonsquare == expected and nonsquare_away_from_one > 0
        character_audits.append(
            {"prime": prime, "nonsquare_A_residues": nonsquare, "zeros": zeros, "away_from_a=1": nonsquare_away_from_one}
        )

    for prime in (3, 5, 7, 11, 13, 17, 19, 23):
        a = choose_target_a(prime)
        s = F(a - 1, 2)
        A = F(1 + 4 * a * a)
        assert vp(a, prime) == vp(A, prime) == vp(s, prime) == 0
        assert legendre(A.numerator, prime) == -1
        e = 1
        Z = F(prime ** (3 * e))
        D = 1 - Z - a * a * Z * Z
        P = canonical_polynomial(a, Z)
        valuations = [vp(coefficient, prime) for coefficient in P]
        assert valuations[0] == valuations[8] == 12 * e
        assert valuations[4] == valuations[5] == 0
        beta = pow((2 * (A.numerator % prime) * (s.numerator % prime) ** 2) % prime, prime - 2, prime)
        primitive = primitive_integer_poly(P)
        derivative = [index * primitive[index] for index in range(1, len(primitive))]
        assert eval_mod(primitive, beta, prime) == 0
        assert eval_mod(derivative, beta, prime) != 0
        assert legendre(2 * beta, prime) == -1
        rows.append(
            {
                "type": "target-newton",
                "label": "PROVED",
                "prime": prime,
                "e_vw_z": e,
                "a": a,
                "s": frac_text(s),
                "case": "w does not divide s",
                "newton_edge": [[4, 0], [5, 0]],
                "simple_residual_root": beta,
                "root_formula": "beta=1/(2*A*s^2) mod w",
                "legendre_2beta": -1,
                "conclusion": "2*beta is nonsquare in the w-adic linear completion",
            }
        )

    # The only first-order survivor: w|s and even v_w(z).  Odd e is killed by
    # the squarefree left residual edge; even e is deliberately left OPEN.
    prime, a = 3, 1
    for e in (1, 2):
        Z = F(prime ** (3 * e))
        P = canonical_polynomial(a, Z)
        valuations = [vp(coefficient, prime) for coefficient in P]
        assert valuations[0] == valuations[8] == 12 * e and valuations[4] == 0
        unit_constant = int((P[0] / prime ** (12 * e)) % prime)
        unit_b4 = int(P[4] % prime)
        assert unit_constant % prime != 0 and unit_b4 % prime != 0
        rows.append(
            {
                "type": "target-newton",
                "label": "PROVED" if e % 2 else "OPEN",
                "prime": prime,
                "e_vw_z": e,
                "a": a,
                "s": "0",
                "case": "w divides s",
                "left_edge": [[0, 12 * e], [4, 0]],
                "left_residual": f"{unit_constant % prime}+{unit_b4 % prime}*T^4 over F_{prime}",
                "left_residual_squarefree": True,
                "conclusion": (
                    "unramified roots have valuation 3*e, odd; 2*beta is nonsquare"
                    if e % 2
                    else "first-order Newton analysis does not decide the square class"
                ),
            }
        )

    return rows, {
        "type": "character-count",
        "label": "PROVED",
        "audited_primes": len(character_audits),
        "formula": "# {a mod w : (1+4a^2|w)=-1}=(w+1)/2 for w=3 mod 4, (w-1)/2 for w=1 mod 4",
        "all_have_choice_away_from_a=1": True,
        "first": character_audits[:5],
        "last": character_audits[-3:],
    }


def theorem_rows() -> list[dict[str, Any]]:
    return [
        {
            "type": "theorem",
            "name": "direct-lambda-identification",
            "label": "PROVED",
            "equation": "(v^2-A)^2(c^2-Ay^2)+64*A*v^2*(B(r^2-As^2)+1)=0, v in {y,r,s}",
            "count": {"base": 2, "Phi_rank": 3, "tied_rank": 2, "DDF_intersection": "3+2-1=4", "total": 6},
            "scope": "sound zero-cost square-branch substitution; lambda=s loses s=0, lambda=y/r makes P witness-dependent; no member theorem follows",
        },
        {
            "type": "theorem",
            "name": "diagonal-orientation-I",
            "label": "PROVED",
            "identification": "X=y, rho=r, lambda=y+r",
            "equations": [
                "C1=y^2-r^2-A=0",
                "C2=-r^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0",
            ],
            "single_polynomial_encoding": "C1^2+C2^2=0 over Q; no additional witness",
            "eliminant": "u=r^2: A*u^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A=0",
            "fiber_degree": 8,
            "finite_flat": True,
            "forces": "the branch-dependent quaternion (P_rho(b),2b) splits",
            "dyadic_obstruction": "PROVED empty on Phi: Q_I/A has root valuation 0 or 4, incompatible with u and A+u both squares when A=5 mod 8",
            "completeness": "PROVED empty on Phi",
        },
        {
            "type": "theorem",
            "name": "diagonal-orientation-II",
            "label": "PROVED",
            "identification": "X=r, rho=y, lambda=r+y",
            "equations": [
                "C1=r^2-y^2-A=0",
                "C2=-y^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0",
            ],
            "eliminant": "u=y^2: A*u^2-(c^2+16AB)u+16A^2B(s^2-1)-16A=0",
            "fiber_degree": 8,
            "finite_flat": True,
            "single_polynomial_encoding": "C1^2+C2^2=0 over Q; no additional witness",
            "target_reduction": "u in {4,-4}, with u a square and A+u a square; necessary residue condition only, and lifting may be singular",
            "forces": "the branch-dependent quaternion (P_rho(b),2b) splits",
            "dyadic_obstruction": "PROVED empty on Phi: any root has v2(u)=2; square u=4t^2 makes Q_II/(16A) nonzero mod 16",
            "completeness": "PROVED empty on Phi",
        },
        {
            "type": "theorem",
            "name": "finite-fiber-count",
            "label": "PROVED",
            "statement": "each diagonal block is the image of a finite flat degree-8 morphism over the A!=0 base, hence tied existential rank <=1",
            "conditional_count": {"base": 2, "Phi_rank": 3, "diagonal_rank_at_most": 1, "intersection_at_most": 3, "total_at_most": 5},
            "warning": "the finite-flat/rank theorem is algebraically valid, but its rational image on Phi is empty, so the <=5 count is vacuous",
        },
        {
            "type": "theorem",
            "name": "canonical-diagonal-slice",
            "label": "PROVED",
            "fixed_point": "rho=2a^2, X=1+2a^2",
            "base_equations": [
                "orientation I: J_I=P_can-4a^4D^2A^2b^4*((1+2a^2)^2-32b)=0",
                "orientation II: J_II=P_can-D^2A^2b^4*(16a^8-32b*(1+2a^2)^2)=0",
            ],
            "degrees": {"J_I": 8, "J_II": 8},
            "endpoints": "both J_0=J_8=4Aa^8Z^4 !=0",
            "conclusion": "recovering the canonical P imposes a nonzero octic equation in either orientation; neither lowers G nor supplies a class-wide section",
        },
        {
            "type": "theorem",
            "name": "rational-exact-norm-match-wall",
            "label": "PROVED",
            "equation": "2b=-A*(1-2As^2b)*q^2",
            "valuation": "v2(left)=1, v2(right)=2*v2(q)",
            "conclusion": "no q in Q exists on Phi; identifying q with y, r, s, or a Phi witness leaves the intersection empty",
        },
        {
            "type": "theorem",
            "name": "Phi-witness-coupling-count",
            "label": "PROVED",
            "known_bound": "2 base +3 Phi witnesses +2 tied witnesses =7",
            "conclusion": "lambda=xi for an internal Phi witness couples the two covers, so the DDF 3+2-1 intersection theorem no longer applies",
            "recovery_of_six": "OPEN without an explicit joint finite-fiber model; Phi rank 3 alone does not name a shareable coordinate",
        },
        {
            "type": "theorem",
            "name": "Capell-square-criterion",
            "label": "PROVED",
            "statement": "for irreducible octic P and beta=b mod P, 2beta square in Q(beta) iff P(U^2/2) is reducible (then 8+8)",
            "norm": "Norm(2beta)=2^8*P_0/P_8=256",
            "warning": "the square norm is necessary, not sufficient",
        },
        {
            "type": "theorem",
            "name": "target-local-Capell-locus",
            "label": "PROVED",
            "necessary_conditions": ["v_w(z) is even", "a=1 mod w (equivalently w|s)", "(5|w)=-1"],
            "scope": "necessary, not sufficient; applies when (A|w)=-1 on the canonical branch",
            "conclusion": "all odd-valuation cells and all s-unit choices are excluded; the even-valuation a=1 residue remains OPEN",
        },
        {
            "type": "theorem",
            "name": "fixed-constant-twist-no-uniformity",
            "label": "PROVED",
            "statement": "no fixed j in Q* makes j*2beta square for every target cell",
            "proof_key": "choose a split prime w outside supp(2j), set v_w(z)=1: w|s gives odd root valuation, while w∤s gives the horizontal root with class j/A and (j|w)=+1, (A|w)=-1",
            "norm_identity_scope": "the displayed H=X^2+ALY^2 cannot be a fixed-j norm via U=X,V=qY when s!=0 because L has a simple zero",
        },
    ]


def dyadic_diagonal_replay() -> dict[str, Any]:
    # Orientation I: a square unit is 1 mod 8; a square of valuation four
    # vanishes mod 8.  Neither can make A+u a square when A=5 mod 8.
    square_residues_mod_8 = {value * value % 8 for value in range(8)}
    assert square_residues_mod_8 == {0, 1, 4}
    assert (5 + 1) % 8 not in square_residues_mod_8
    assert 5 not in square_residues_mod_8

    # Orientation II after Newton forces u=4*t^2:
    # Q_II/(16A) = t^4 - 8*h*t^2 + e,
    # e=2*A*b*(s^2-1)-1.  Exhaust all relevant residues modulo 16.
    checked = 0
    for s in range(16):
        for b in range(1, 16, 2):
            for t in range(1, 16, 2):
                for h in range(1, 16, 2):
                    residue = (
                        t**4
                        - 8 * h * t * t
                        + 2 * 5 * b * (s * s - 1)
                        - 1
                    ) % 16
                    assert residue != 0
                    checked += 1
    assert checked == 8192
    return {
        "type": "dyadic-diagonal-no-go",
        "label": "PROVED",
        "v2_c_lower_bound": 4,
        "orientation_I": {
            "equation": "u^2+E*u+16*K=0 with E,K 2-adic units",
            "root_valuations": [0, 4],
            "squareclass_contradiction": "u square gives A+u=6 mod 8 or 5 mod 8",
        },
        "orientation_II": {
            "equation": "u^2-32*h*u+16*e=0 with h,e 2-adic units",
            "root_valuation": 2,
            "square_substitution": "u=4*t^2 => q/16=t^4-8*h*t^2+2*A*b*(s^2-1)-1 !=0 mod 16",
            "residue_tuples_checked": checked,
        },
        "conclusion": "both diagonal covers have empty rational image on Phi",
        "count_consequence": "the finite-flat degree-8 and tied-rank <=1 statements remain true, but the <=5 count is vacuous",
    }


def direct_target_orientation_I_obstruction() -> dict[str, Any]:
    prime = 3
    A = 2
    points = []
    for y in range(prime):
        for r in range(prime):
            C1 = (y * y - r * r - A) % prime
            C2_target = (A * r * r * y * y - 16 * A) % prime
            if C1 == 0 and C2_target == 0:
                points.append([y, r])
    assert points == []
    return {
        "type": "target-residue-obstruction",
        "label": "PROVED",
        "orientation": "X=y,rho=r",
        "prime": prime,
        "A_residue": A,
        "equations": ["y^2-r^2=A", "r^2*y^2=16"],
        "rational_points_mod_prime": points,
        "scope": "standard W1/L20 residue subroute v_w(b)>=1, hence B=c=0 mod w",
        "conclusion": "orientation I has no residue point in this w=3 subroute; this is scoped corroboration, while the separate Q_2 theorem proves Phi-emptiness",
    }


def write_report(rows: list[dict[str, Any]], elapsed: float) -> None:
    capell_rows = [row for row in rows if row.get("type") == "capell-certificate"]
    newton_rows = [row for row in rows if row.get("type") == "target-newton"]
    lines = [
        "# L23 - zero-cost parameter absorption audit",
        "",
        "**Overall verdict.**  Exact zero-cost self-coupling exists algebraically, but it cannot remove the member-existence hypothesis: both diagonal covers are **PROVED empty on `Phi` over `Q_2`**.  Before imposing `Phi`, identifying the square-branch hyperbola coordinates with the tied variables gives an exact branch point, forces the resulting branch-dependent quaternion to split, and makes each tied cover finite flat of degree eight.  Its tied rank is at most one and the formal intersection architecture would use at most five unknowns, but that count is vacuous because the rational image on `Phi` is empty.  The orientation-I `w=3`, `v_w(b)>=1` residue block remains only a scoped corroboration, not the global proof.  Recovering the fixed canonical octic also restores a nonzero degree-eight equation.",
        "",
        "The rational exact-norm section remains **PROVED empty** on `Phi` by the 2-adic wall.  The algebraic Capell escape is **PROVED impossible** on a large exact target-local locus and on every certified row below, but the residual even-valuation, `a=1 mod w`, `(5|w)=-1` locus remains **OPEN**.  No finite scan is treated as negative evidence.  Schinzel H is not removed.",
        "",
        "## 1. Baseline and notation",
        "",
        "Put",
        "",
        "```text",
        "a=1+2s,  A=1+4a^2,  B=2b,  Z=z^3,  D=1-Z-a^2Z^2,",
        "Ng=16a^4b^2-A(b-1)^4.",
        "```",
        "",
        "On the square branch use hyperbola coordinates",
        "",
        "```text",
        "X=(lambda+A/lambda)/2,  rho=(lambda-A/lambda)/2,",
        "X^2-rho^2=A,  tau=X/A,  delta=-rho^2/A,  alpha=rho^2.",
        "```",
        "",
        "The exact tied equation in the original tied variables `(y,r)` is",
        "",
        "```text",
        "delta(c^2-Ay^2)-16B(r^2-As^2)=16.                 (T)",
        "```",
        "",
        "The original accounting is `2 + (3+2-1) = 6`: base `(b,s)`, rank-three `Phi`, rank-two tied block, and the DDF intersection saving one.",
        "",
        "## 2. Direct lambda reuse: sound, but still rank two",
        "",
        "For `v=y`, `r`, or `s`, substituting `lambda=v` and clearing `4Av^2` gives exactly",
        "",
        "```text",
        "(v^2-A)^2(c^2-Ay^2)+64Av^2(B(r^2-As^2)+1)=0.     (D_v)",
        "```",
        "",
        "There is no hidden `v=0` rational branch on the admissible L6 base: at `v=0`, `(D_v)` says `c^2=Ay^2`; here `A` is nonsquare in `Q_2` and `c!=0`.  The latter follows from `Ng=0 => A=(4a^2b/(b-1)^2)^2`, impossible (with `b=1` checked separately).  Thus every rational solution is an honest square-branch solution.",
        "",
        "For `v=s`, the branch is fixed after the base is chosen, but its octic endpoints are `rho^2*A*a^4*Z^4 !=0`, so the degree stays eight.  For `v=y` or `r`, `rho` depends on a tied witness; there is no fixed univariate `G(t)` to which L19 can be applied.  Each `(D_v)` is still one equation in two tied variables, with generic fibre dimension one and tied rank at most two.  The count stays six.",
        "",
        "A finite disjunction is unrelated: multiplying finitely many fixed-branch equations reuses `(y,r)`, but it does not quantify an infinite `lambda`.  Calling the latter free because the former is free is invalid.",
        "",
        "## 3. Diagonal hyperbola absorption",
        "",
        "### Orientation I: `X=y`, `rho=r`",
        "",
        "Set `lambda=y+r` and impose",
        "",
        "```text",
        "C1 = y^2-r^2-A = 0,",
        "C2 = -r^2(c^2-Ay^2)-16AB(r^2-As^2)-16A = 0.      (I)",
        "```",
        "",
        "Because `(y-r)(y+r)=A!=0`, `lambda` is nonzero.  A rational solution cannot have `r=0`, since that would make `A=y^2`, contrary to `A` being nonsquare in `Q_2`.  Hence `delta=-r^2/A` is nonzero and `(I)` is exactly `(T)`, not a surrogate.",
        "",
        "With `u=r^2`, elimination of `y^2=A+u` gives",
        "",
        "```text",
        "A u^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A=0.          (Q_I)",
        "```",
        "The global obstruction is dyadic.  Write `Q_I/A=u^2+E*u+16K`, where `E=A-c^2/A-32b` and `K=2Abs^2-1`.  On `Phi`, `v_2(c)>=4`, `A=5 mod 8`, and `b` is a unit, so both `E` and `K` are units.  Newton valuations force any root to have `v_2(u)=0` or `4`.  But `u=r^2` and `A+u=y^2`: in the first case `u=1 mod 8` and `A+u=6 mod 8` has odd valuation; in the second `A+u=5 mod 8`.  Neither is a square.  Thus orientation I is **PROVED empty on `Phi`**.",
        "",
        "",
        "In the standard W1/L20 target subroute `v_w(b)>=1`, one has `B=c=0` in the residue field, so `u(A+u)=16`.  Together with `u=r^2`, `A+u=y^2`, this says `(ry)^2=16`.  At `w=3`, nonsquare `A=2`; the first equation forces `(y^2,r^2)=(0,1)`, contradicting `(ry)^2=1`.  Thus this **specific residue block is PROVED empty**.  This is not a global incompleteness theorem: when `v_w(b)=0`, `B` survives in the residue equations and points can exist.",
        "",
        "### Orientation II: `X=r`, `rho=y`",
        "",
        "The swapped exact system is",
        "",
        "```text",
        "C1' = r^2-y^2-A = 0,",
        "C2' = -y^2(c^2-Ay^2)-16AB(r^2-As^2)-16A = 0,",
        "A u^2-(c^2+16AB)u+16A^2B(s^2-1)-16A=0, u=y^2.  (Q_II)",
        "```",
        "",
        "It has the same exact soundness.  In the `B=c=0` target reduction, `u^2=16`, hence `u=+4` or `u=-4`.  Since `u=y^2` and `r^2=A+u`, the exact necessary residue disjunction is: `(u=4 and A+4 is a square)` or `(u=-4 is a square and A-4 is a square)`.  This condition is not sufficient for a rational or local point; the residue point can be singular and a Hensel lift still has to be proved.  The `w=3`, `u=4` residue point is therefore only a target-local observation; the global dyadic calculation below proves that it cannot lift to a rational `Phi` point.",
        "The global obstruction is again dyadic and stronger.  Put `L=c^2/A+32b=32h` and `e=2Ab(s^2-1)-1`; `h,e` are 2-adic units.  Then `Q_II/A=u^2-32hu+16e`.  Newton forces `v_2(u)=2`.  If `u=y^2`, write `u=4t^2` with `t` a unit.  Division by `16` gives",
        "",
        "```text",
        "t^4-8h t^2+2Ab(s^2-1)-1.",
        "```",
        "",
        "Modulo `16` this is `8+2Ab(s^2-1)`: it is `8` when `s` is odd and has valuation one when `s` is even.  It never vanishes.  Therefore orientation II is also **PROVED empty on `Phi`**.",
        "",
        "",
        "In either orientation, `rho^2` is a square and the displayed tied point proves that the branch-dependent quaternion `(P_rho(b),2b)` splits.  This is a theorem about a self-coupled branch, not the fixed canonical `P`.",
        "",
        "Over `Q`, the conjunction costs no variable: `C1^2+C2^2=0` has exactly the same rational solutions.  For the fibre argument below one uses the actual complete intersection `V(C1,C2)`, whose `Q`-point image is the same; one must not measure the geometric fibres of the sum-of-squares hypersurface.",
        "",
        "## 4. Finite fibres and exact quantifier accounting",
        "",
        "After localizing at the already-proved nonzero `A`, each system is equivalent to a monic quartic in the identified tied square variable (`r` in I, `y` in II) followed by one monic quadratic for the other variable.  The coordinate algebra is free with basis of size `4*2=8`; each complete-intersection cover is therefore finite flat of degree eight over the guarded base.",
        "",
        "By the DDF rank/essential-fibre-dimension theorem, the image of this finite cover has tied existential rank at most one (or rank zero if it is quantifier-free).  Consequently its intersection with rank-three `Phi` has rank at most three: `3+1-1=3` in the positive-rank case, and the elementary conjunction bound is three in the rank-zero case.  Adding `(b,s)` gives at most `2+3=5` unknowns.",
        "",
        "The finite-flat and rank statements are algebraically **PROVED**, but the nominal five-count is vacuous: both diagonal images have no rational point on `Phi`.  Neither can define any part of `Q\\Z` in the L6 architecture.",
        "",
        "If one freezes the diagonal point back to the canonical branch, `rho=2a^2` and `X=1+2a^2`, the two orientations impose",
        "",
        "```text",
        "J_I(b)  = P_can(b)-4a^4D^2A^2b^4((1+2a^2)^2-32b)=0,",
        "J_II(b) = P_can(b)-D^2A^2b^4(16a^8-32b(1+2a^2)^2)=0.",
        "```",
        "",
        "Each correction has degree at most five, while both polynomials have `J_0=J_8=4Aa^8Z^4!=0`; both remain nonzero octics.  Hence neither canonical slice lowers `G` nor gives a class-wide section.",
        "",
        "## 5. Why a Phi witness is not a free coordinate",
        "",
        "The input fact is only that `Phi` has existential rank at most three.  Naming one witness `xi` and setting `lambda=xi` couples the tied equation to the fibre of a particular rank-three presentation.  The projected set is no longer the intersection of a rank-three subset and an independent rank-two subset of the common `(b,s,z)` base, so DDF's `3+2-1` theorem does not apply.  The available direct count is `2+3+2=7`.  Recovering six would require an explicit joint finite-fibre model; it does not follow from the black-box rank bound and is **OPEN**.",
        "",
        "## 6. Rational exact norm matching is empty",
        "",
        "Let `L=1-2As^2b`.  On `Phi`, `v_2(a)=v_2(b)=0` and `v_2(s)>=0`, so `A` and `L` are 2-adic units (`L=1 mod 2`).  The tempting exact match",
        "",
        "```text",
        "2b=-A L q^2",
        "```",
        "",
        "has left valuation one and right valuation `2v_2(q)`, an even integer.  It has no rational solution, including `q=0`.  Reusing `y`, `r`, `s`, or any `Phi` witness for `q` cannot change this parity contradiction.  This proves the tested norm-section route empty, rather than merely failing to find examples.",
        "",
        "## 7. The algebraic Capell escape and its exact local locus",
        "",
        "For the canonical octic let `beta` be the class of `b` in `K=Q[b]/(P)`.  Capell's degree argument gives",
        "",
        "```text",
        "2 beta is a square in K  <=>  P(U^2/2) is reducible over Q.",
        "```",
        "",
        "In the square case the degree-sixteen composition is a product of two degree-eight conjugate factors.  The endpoint identity",
        "",
        "```text",
        "P_0=P_8=4Aa^8Z^4,  Norm_K/Q(2 beta)=2^8 P_0/P_8=256",
        "```",
        "",
        "shows only that the norm obstruction vanishes; it does not prove the element square.",
        "",
        "There is, however, an exact target-local restriction.  Let `e=v_w(z)>=1`, so `v_w(Z)=3e`, and choose `(A|w)=-1`.",
        "",
        "* If `w` does not divide `s`, the Newton polygon has a horizontal edge from `(4,0)` to `(5,0)`.  Its simple residual root is `beta=1/(2As^2)`, and `(2beta|w)=(1/A|w)=-1`.  Hensel gives a `Q_w` embedding at which `2beta` is nonsquare.",
        "* If `w|s` and `e` is odd, the left edge from `(0,12e)` to `(4,0)` has squarefree residual `c_0+c_4T^4`.  Its roots lie in unramified extensions and have valuation `3e`, which is odd; again `2beta` is nonsquare.",
        "",
        "Therefore any square locus must satisfy all three conditions",
        "",
        "```text",
        "v_w(z) even,  a=1 mod w (w|s),  and (5|w)=-1.",
        "```",
        "",
        "The last condition follows because then `A=5 mod w` while `(A|w)=-1`.  These conditions are necessary, not sufficient.  The even-valuation exceptional residue is **OPEN**.  The character count also proves that L20 can always choose a nonsquare-`A` residue away from `a=1`; that deliberately lands in the proved no-go case and cannot be misread as excluding the exceptional choice.",
        "",
        "A fixed constant twist cannot repair this uniformly.  For fixed `j`, quadratic Chebotarev supplies target primes `w` outside `supp(2j)` with `(j|w)=+1`.  On an `e=1` cell, `w|s` leaves the odd-valuation obstruction, while `w∤s` gives horizontal-root class `j/A`, also nonsquare.  Thus no fixed `j` makes `j*2beta` square in every cell.  Moreover a same-shape fixed-`j` norm identity `H=X^2-j(qY)^2` would require `AL=-jq^2`; for `s!=0`, `L` has a simple zero in `Q(b)`, so this is impossible.  Nonconstant twists require a new argument and leave the norm condition `(P(b),j)=1`; no closure is claimed.",
        "",
        "## 8. Exact certificate rows",
        "",
        f"The artifact contains **{len(capell_rows)}** parameterwise theorems, not a no-hit statistic.  Each row first proves its canonical octic irreducible by an irreducible reduction, then exhibits a different (or possibly equal) prime with a simple root `r` for which `(2r|p)=-1`.  Hensel embeds the octic field into `Q_p`, proving `2beta` nonsquare for that exact `(a,Z)`.  The rows cover `v_2(Z)<0`, `=0`, and `>0`; untested parameters remain unjudged.",
        "",
        f"The target audit has **{len(newton_rows)}** exact rows, including the horizontal-root theorem and both the proved odd-valuation and open even-valuation `w|s` cases.",
        "",
        "Labels: direct and diagonal identities **PROVED**; diagonal finite-flat degree-eight/rank count **PROVED but vacuous**; both diagonal rational images on `Phi` **PROVED empty over `Q_2`**; the orientation-I standard `w=3` W1/L20 residue block remains a scoped **PROVED** corroboration; rational norm match **PROVED empty**; Capell equivalence and target-local necessary locus **PROVED**; residual even-valuation square locus **OPEN**; finite parameter certificates **PROVED per row**; global Capell absence is not inferred.",
        "",
        f"Verification wall-clock: **{elapsed:.3f} s**.  Pacing: {PACER.steps} finite-loop steps, {PACER.sleeps} sleeps of {PACE_SECONDS} s.",
        "",
        "Artifacts: `math/h10q/l23_absorption.py`, `math/h10q/data/l23_absorption.jsonl`, `/tmp/l23_absorption.md`.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    rows: list[dict[str, Any]] = [
        {
            "type": "meta",
            "artifact": "l23_absorption",
            "labels": ["PROVED", "CONDITIONAL", "EVIDENCE", "OPEN"],
            "finite_scan_policy": "positive certificates only; absence/refusal is never negative evidence",
            "stdlib_only": True,
            "certificate_prime_limit": CERTIFICATE_LIMIT,
        },
        identity_replay(),
        dyadic_diagonal_replay(),
        *theorem_rows(),
        direct_target_orientation_I_obstruction(),
    ]
    target_rows, character_row = target_newton_rows()
    rows.append(character_row)
    rows.extend(target_rows)
    rows.extend(capell_certificate_rows())
    elapsed = time.perf_counter() - started
    summary = {
        "type": "summary",
        "label": "PROVED",
        "rows": len(rows) + 1,
        "capell_parameterwise_nonsquare_certificates": sum(row.get("type") == "capell-certificate" for row in rows),
        "target_newton_rows": sum(row.get("type") == "target-newton" for row in rows),
        "zero_cost_identifications_proved": ["lambda=y", "lambda=r", "lambda=s", "diagonal I", "diagonal II"],
        "diagonal_phi_image": "PROVED empty over Q_2; formal <=5 count is vacuous",
        "record_closure": "OPEN",
        "schinzel_removed": False,
        "elapsed_seconds": elapsed,
        "pacing": {"steps": PACER.steps, "sleeps": PACER.sleeps, "sleep_seconds": PACE_SECONDS},
    }
    rows.append(summary)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(rows, elapsed)
    print(
        f"l23_absorption: rows={len(rows)} capell={summary['capell_parameterwise_nonsquare_certificates']} "
        f"newton={summary['target_newton_rows']} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
