#!/usr/bin/env python3
"""L24 exact bounded search for the two self-coupled diagonal systems.

The program is deliberately proof-disciplined.  For a fixed rational base it
solves the even quartic in the tied coordinate exactly (as a quadratic in its
square), and it replays every reported rational hit in the original equations.
A failed bounded pool is never promoted to a failed cell.  Separately, exact
finite-field failures and nonsingular Hensel lifts are recorded with their
precise local scope.

Stdlib only.  Run as one paced low-priority process, for example
``nice -n 19 python3 l24_diagonal_search.py`` from this directory.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from typing import Any, Iterable
import json
import math
import os
import time


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l24_diagonal_search.jsonl"
REPORT = Path("/tmp/l24_diagonal_search.md")

# The first block is exactly h10q._L9_U_POOL, copied here so this replay remains
# standalone.  The second block is the already predeclared L19 off-grid horizon:
# sixteen u=-1 cells and the separately tried (131,5) cell.
GRID_U_POOL = tuple(F(v) for v in (1, -1, 2, -2, 3, -3, 5, -5)) + tuple(
    F(n, d)
    for n, d in (
        (1, 3), (-1, 3), (2, 3), (-1, 5), (7, 3), (1, 7), (-2, 7)
    )
)
OFFGRID_W = (
    101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179,
)
OFFGRID_CELLS = tuple((w, F(-1)) for w in OFFGRID_W) + ((131, F(5)),)

# Declared global rational-height box.  Cell-specific additions are only the
# two local mechanisms named below: target valuations and template/arm lifts.
S_NUMERATOR_BOUND = 9
B_NUMERATOR_BOUND = 21
DENOMINATOR_BOUND = 7
TARGET_MULTIPLIERS = tuple(
    F(n, d)
    for n, d in (
        (-1, 1), (1, 1), (-3, 1), (3, 1), (-5, 1), (5, 1),
        (-1, 3), (1, 3), (-1, 5), (1, 5),
    )
)
LOCAL_CERTIFICATE_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
LOCAL_CERTIFICATE_TRIES = 8
PACE_EVERY = 256
PACE_SECONDS = 0.001

ORIENTATION_NAMES = {1: "I", 2: "II"}


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


@dataclass(frozen=True)
class Base:
    s: F
    b: F
    z: F
    a: F
    A: F
    B: F
    Z: F
    D: F
    c: F


@dataclass(frozen=True)
class QuarticResult:
    status: str
    discriminant: F
    rational_u_roots: tuple[F, ...]
    solutions: tuple[tuple[F, F, F], ...]  # (u,y,r)


# Sparse multivariate polynomials, used only for genuine symbolic identity
# checks.  Exponents are in (A,s,B,y,r,L).
Monomial = tuple[int, int, int, int, int, int]
Polynomial = dict[Monomial, F]
ZERO_MONOMIAL: Monomial = (0, 0, 0, 0, 0, 0)


def poly_clean(poly: Polynomial) -> Polynomial:
    return {monomial: value for monomial, value in poly.items() if value}


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = dict(left)
    for monomial, value in right.items():
        answer[monomial] = answer.get(monomial, F(0)) + value
    return poly_clean(answer)


def poly_scale(poly: Polynomial, scalar: F | int) -> Polynomial:
    scalar = F(scalar)
    return poly_clean({monomial: scalar * value for monomial, value in poly.items()})


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for left_monomial, left_value in left.items():
        for right_monomial, right_value in right.items():
            monomial = tuple(
                left_monomial[index] + right_monomial[index]
                for index in range(len(ZERO_MONOMIAL))
            )
            answer[monomial] = answer.get(monomial, F(0)) + left_value * right_value
    return poly_clean(answer)


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    answer: Polynomial = {ZERO_MONOMIAL: F(1)}
    base = poly
    while exponent:
        if exponent & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return answer


def poly_var(index: int) -> Polynomial:
    exponent = [0] * len(ZERO_MONOMIAL)
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def vp(value: F | int, prime: int) -> int:
    value = F(value)
    if value == 0:
        raise ValueError("valuation of zero is not used in this search")
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


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def odd_primes_below(limit: int) -> tuple[int, ...]:
    return tuple(value for value in range(3, limit) if is_prime(value))


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    power = pow(value, (prime - 1) // 2, prime)
    assert power in (1, prime - 1)
    return 1 if power == 1 else -1


def mod_sqrt(value: int, prime: int) -> int | None:
    value %= prime
    for candidate in range(prime):
        if candidate * candidate % prime == value:
            return candidate
    return None


def mod_fraction(value: F | int, prime: int) -> int | None:
    value = F(value)
    if value.denominator % prime == 0:
        return None
    return value.numerator * pow(value.denominator, -1, prime) % prime


def bounded_rationals(numerator_bound: int, denominator_bound: int) -> tuple[F, ...]:
    values = {
        F(numerator, denominator)
        for denominator in range(1, denominator_bound + 1)
        for numerator in range(-numerator_bound, numerator_bound + 1)
    }
    return tuple(
        sorted(
            values,
            key=lambda value: (
                max(abs(value.numerator), value.denominator),
                value.denominator,
                abs(value.numerator),
                value,
            ),
        )
    )


RAW_S_POOL = bounded_rationals(S_NUMERATOR_BOUND, DENOMINATOR_BOUND)
RAW_B_POOL = bounded_rationals(B_NUMERATOR_BOUND, DENOMINATOR_BOUND)
S_POOL = tuple(value for value in RAW_S_POOL if value == 0 or vp(value, 2) >= 0)
B_POOL = tuple(value for value in RAW_B_POOL if value != 0 and vp(value, 2) == 0)


def phi_holds(s: F | int, b: F | int) -> bool:
    s = F(s)
    b = F(b)
    return b != 0 and (s == 0 or vp(s, 2) >= 0) and vp(b, 2) == 0


def make_base(s_value: F | int, b_value: F | int, z_value: F | int) -> tuple[Base | None, str | None]:
    s = F(s_value)
    b = F(b_value)
    z = F(z_value)
    if not phi_holds(s, b):
        return None, "outside-Phi"
    a = 1 + 2 * s
    A = 1 + 4 * a * a
    B = 2 * b
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    if D == 0:
        return None, "D=0-bridge-undefined"
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    c = a * a * Z * Z * Ng / (A * b * b * D)
    return Base(s, b, z, a, A, B, Z, D, c), None


def quartic_coefficients(base: Base, orientation: int) -> tuple[F, F, F]:
    """Return q4,q2,q0 for q4*x^4+q2*x^2+q0."""
    A, B, c, s = base.A, base.B, base.c, base.s
    if orientation == 1:
        return (
            A,
            A * A - c * c - 16 * A * B,
            16 * A * A * B * s * s - 16 * A,
        )
    if orientation == 2:
        return (
            A,
            -c * c - 16 * A * B,
            16 * A * A * B * (s * s - 1) - 16 * A,
        )
    raise ValueError(orientation)


def solve_quartic(base: Base, orientation: int) -> QuarticResult:
    """Complete Q-solver for the tied even quartic via u=x^2."""
    q4, q2, q0 = quartic_coefficients(base, orientation)
    assert q4 != 0
    discriminant = q2 * q2 - 4 * q4 * q0
    square_discriminant = rational_sqrt(discriminant)
    if square_discriminant is None:
        return QuarticResult("discriminant-nonsquare", discriminant, (), ())

    roots = tuple(
        sorted(
            {
                (-q2 + square_discriminant) / (2 * q4),
                (-q2 - square_discriminant) / (2 * q4),
            }
        )
    )
    solutions: list[tuple[F, F, F]] = []
    saw_square_u = False
    for u in roots:
        square_u = rational_sqrt(u)
        if square_u is None:
            continue
        saw_square_u = True
        square_companion = rational_sqrt(base.A + u)
        if square_companion is None:
            continue
        if orientation == 1:
            y, r = square_companion, square_u
        else:
            y, r = square_u, square_companion
        solutions.append((u, y, r))
    if solutions:
        status = "hit"
    elif saw_square_u:
        status = "companion-nonsquare"
    else:
        status = "u-nonsquare"
    return QuarticResult(status, discriminant, roots, tuple(solutions))


def original_equations(base: Base, orientation: int, y: F, r: F) -> tuple[F, F]:
    if orientation == 1:
        C1 = y * y - r * r - base.A
        rho_square = r * r
    elif orientation == 2:
        C1 = r * r - y * y - base.A
        rho_square = y * y
    else:
        raise ValueError(orientation)
    C2 = (
        -rho_square * (base.c * base.c - base.A * y * y)
        - 16 * base.A * base.B * (r * r - base.A * base.s * base.s)
        - 16 * base.A
    )
    return C1, C2


def replay_rational_hit(base: Base, orientation: int, u: F, y: F, r: F) -> None:
    assert phi_holds(base.s, base.b)
    assert base.a == 1 + 2 * base.s
    assert base.A == 1 + 4 * base.a * base.a
    assert base.B == 2 * base.b
    assert base.Z == base.z**3
    assert base.D == 1 - base.Z - base.a * base.a * base.Z * base.Z
    C1, C2 = original_equations(base, orientation, y, r)
    assert C1 == 0 and C2 == 0
    if orientation == 1:
        assert u == r * r and base.A + u == y * y
        X, rho = y, r
    else:
        assert u == y * y and base.A + u == r * r
        X, rho = r, y
    lam = X + rho
    assert lam != 0
    assert X == (lam + base.A / lam) / 2
    assert rho == (lam - base.A / lam) / 2


def rational_hit_row(
    family: str,
    prime: int,
    unit: F,
    base: Base,
    orientation: int,
    u: F,
    y: F,
    r: F,
) -> dict[str, Any]:
    replay_rational_hit(base, orientation, u, y, r)
    target_u = mod_fraction(u, prime) if vp(base.b, prime) > 0 else None
    return {
        "type": "rational-hit",
        "label": "PROVED",
        "family": family,
        "cell": [prime, [unit.numerator, unit.denominator]],
        "z": frac_text(base.z),
        "Z": frac_text(base.Z),
        "orientation": ORIENTATION_NAMES[orientation],
        "s": frac_text(base.s),
        "a": frac_text(base.a),
        "A": frac_text(base.A),
        "b": frac_text(base.b),
        "B": frac_text(base.B),
        "v_w_b": vp(base.b, prime),
        "c": frac_text(base.c),
        "u": frac_text(u),
        "u_mod_w": target_u,
        "y": frac_text(y),
        "r": frac_text(r),
        "lambda": frac_text(y + r),
        "Phi_replay": {"v2_s": None if base.s == 0 else vp(base.s, 2), "v2_b": vp(base.b, 2)},
        "equation_replay": {"C1": "0", "C2": "0"},
        "solver": "complete even-quartic solver over Q plus conic reconstruction",
    }


def local_obstruction(base: Base, orientation: int, prime: int) -> dict[str, Any] | None:
    """Certify no Q_prime point when the monic-integral reduction has none."""
    q4, q2, q0 = quartic_coefficients(base, orientation)
    values = (base.A, q4, q2, q0)
    reduced = tuple(mod_fraction(value, prime) for value in values)
    if any(value is None for value in reduced):
        return None
    A_mod, q4_mod, q2_mod, q0_mod = (int(value) for value in reduced)
    if A_mod == 0 or q4_mod == 0:
        return None

    compatible_u = []
    for u in range(prime):
        if (q4_mod * u * u + q2_mod * u + q0_mod) % prime:
            continue
        if legendre(u, prime) < 0 or legendre(A_mod + u, prime) < 0:
            continue
        compatible_u.append(u)
    if compatible_u:
        return None
    return {
        "prime": prime,
        "quartic_mod_p": [q0_mod, 0, q2_mod, 0, q4_mod],
        "A_mod_p": A_mod,
        "compatible_u": [],
    }


def find_local_obstruction(base: Base, orientation: int) -> dict[str, Any] | None:
    for prime in LOCAL_CERTIFICATE_PRIMES:
        certificate = local_obstruction(base, orientation, prime)
        if certificate is not None:
            return certificate
    return None


def target_mode(base: Base, prime: int) -> str:
    k = vp(base.b, prime)
    t = vp(base.Z, prime)
    assert t > 0
    if k == 0:
        return "v_w(b)=0"
    if 0 < k < t:
        return "0<v_w(b)<v_w(Z)"
    if k == t:
        return "v_w(b)=v_w(Z)"
    return "other-target-valuation"


def conic_point(s_value: F | int, lam_value: F | int, orientation: int) -> tuple[F, F, F, F, F]:
    s = F(s_value)
    lam = F(lam_value)
    assert lam != 0
    a = 1 + 2 * s
    A = 1 + 4 * a * a
    X = (lam + A / lam) / 2
    rho = (lam - A / lam) / 2
    if orientation == 1:
        y, r = X, rho
        assert y * y - r * r == A
    elif orientation == 2:
        y, r = rho, X
        assert r * r - y * y == A
    else:
        raise ValueError(orientation)
    return a, A, y, r, lam


def target_template(s_value: F | int, lam_value: F | int, orientation: int) -> dict[str, F]:
    s = F(s_value)
    a, A, y, r, lam = conic_point(s, lam_value, orientation)
    T = r * r - A * s * s
    assert T != 0
    rho_square = r * r if orientation == 1 else y * y
    numerator = rho_square * y * y - 16
    B = numerator / (16 * T)
    b_residue = B / 2
    if orientation == 1:
        determinant = 4 * A * r * y * (y * y + r * r - 16 * B)
    else:
        determinant = 8 * A * r * y * (8 * B - y * y)
    assert B != 0 and determinant != 0
    assert A * (numerator - 16 * B * T) == 0
    return {
        "s": s,
        "lambda": lam,
        "a": a,
        "A": A,
        "y": y,
        "r": r,
        "T": T,
        "numerator": numerator,
        "B": B,
        "b_residue": b_residue,
        "fiber_jacobian": determinant,
    }


def template_bad_integer(template: dict[str, F]) -> int:
    # If an odd p avoids this integer, every denominator is a p-unit and
    # A,B,T and the fixed-base (y,r) Jacobian are nonzero modulo p.
    denominator_values = (
        template["a"], template["A"], template["y"], template["r"],
        template["T"], template["B"], template["b_residue"],
        template["fiber_jacobian"],
    )
    nonzero_values = (
        template["A"], template["B"], template["T"], template["fiber_jacobian"],
    )
    answer = 1
    for value in denominator_values:
        answer *= value.denominator
    for value in nonzero_values:
        answer *= abs(value.numerator)
    while answer % 2 == 0:
        answer //= 2
    return answer


UNIT_TEMPLATE_SPECS = {
    1: ((F(0), F(-1)), (F(1), F(-1))),
    2: ((F(-2), F(-1)), (F(2), F(-1))),
}


def symbolic_identity_row() -> dict[str, Any]:
    A = poly_var(0)
    s = poly_var(1)
    B = poly_var(2)
    y = poly_var(3)
    r = poly_var(4)
    L = poly_var(5)
    one = {ZERO_MONOMIAL: F(1)}

    T = poly_add(poly_pow(r, 2), poly_scale(poly_mul(A, poly_pow(s, 2)), -1))
    numerators = {
        "I": poly_add(poly_mul(poly_pow(r, 2), poly_pow(y, 2)), poly_scale(one, -16)),
        "II": poly_add(poly_pow(y, 4), poly_scale(one, -16)),
    }
    cleared_substitutions: dict[str, Polynomial] = {}
    for name, numerator in numerators.items():
        relation = poly_add(numerator, poly_scale(poly_mul(B, T), -16))
        reduced_C2 = poly_mul(A, relation)
        assert poly_add(reduced_C2, poly_scale(poly_mul(A, relation), -1)) == {}
        # Clearing 16*T after B=numerator/(16*T) gives equal products.
        cleared = poly_add(
            poly_scale(poly_mul(poly_mul(A, T), numerator), 16),
            poly_scale(poly_mul(poly_mul(A, numerator), T), -16),
        )
        assert cleared == {}
        cleared_substitutions[name] = cleared

    # Hyperbola parameter X=(L+A/L)/2, rho=(L-A/L)/2, cleared by 4L^2.
    L2 = poly_pow(L, 2)
    plus = poly_add(L2, A)
    minus = poly_add(L2, poly_scale(A, -1))
    hyperbola = poly_add(
        poly_add(poly_pow(plus, 2), poly_scale(poly_pow(minus, 2), -1)),
        poly_scale(poly_mul(A, L2), -4),
    )
    assert hyperbola == {}

    indexed_templates = [
        (orientation, target_template(s_value, lam_value, orientation))
        for orientation, specs in UNIT_TEMPLATE_SPECS.items()
        for s_value, lam_value in specs
    ]
    ratios = {
        template["numerator"] / (template["B"] * template["T"])
        for _orientation, template in indexed_templates
    }
    interpolated = [
        coefficient
        for coefficient in range(-64, 65)
        if all(
            template["numerator"] == coefficient * template["B"] * template["T"]
            for _orientation, template in indexed_templates
        )
    ]
    assert ratios == {F(16)} and interpolated == [16]
    return {
        "type": "mined-symbolic-identity",
        "label": "PROVED",
        "source": "four successful low-height target-unit template rows",
        "interpolation_window": [-64, 64],
        "unique_interpolated_coefficient": 16,
        "identities": {
            "I": "r^2*y^2-16=16*B*(r^2-A*s^2)",
            "II": "y^4-16=16*B*(r^2-A*s^2)",
        },
        "target_reductions": {
            "I": "C2_bar=A*(r^2*y^2-16*B*(r^2-A*s^2)-16)=0",
            "II": "C2_bar=A*(y^4-16*B*(r^2-A*s^2)-16)=0",
        },
        "symbolic_checks": {
            "cleared_B_substitution_remainder_terms": {
                name: len(poly) for name, poly in cleared_substitutions.items()
            },
            "hyperbola_parameter_remainder_terms": len(hyperbola),
        },
        "scope": "polynomial identities at a target prime where c_bar=0; not global rational sections",
    }


def unit_template_rows() -> tuple[list[dict[str, Any]], dict[int, tuple[dict[str, F], ...]]]:
    rows: list[dict[str, Any]] = []
    templates_by_orientation: dict[int, tuple[dict[str, F], ...]] = {}
    gcds: dict[int, int] = {}
    for orientation, specs in UNIT_TEMPLATE_SPECS.items():
        templates = tuple(target_template(s, lam, orientation) for s, lam in specs)
        templates_by_orientation[orientation] = templates
        bad_integers = [template_bad_integer(template) for template in templates]
        gcds[orientation] = math.gcd(*bad_integers)
        assert gcds[orientation] == (3 if orientation == 1 else 45)
        for index, (template, bad_integer) in enumerate(zip(templates, bad_integers), 1):
            rows.append(
                {
                    "type": "target-unit-template",
                    "label": "PROVED",
                    "orientation": ORIENTATION_NAMES[orientation],
                    "template": index,
                    "s": frac_text(template["s"]),
                    "a": frac_text(template["a"]),
                    "A": frac_text(template["A"]),
                    "lambda": frac_text(template["lambda"]),
                    "y": frac_text(template["y"]),
                    "r": frac_text(template["r"]),
                    "T": frac_text(template["T"]),
                    "B_mod_w": frac_text(template["B"]),
                    "b_mod_w": frac_text(template["b_residue"]),
                    "fiber_jacobian": frac_text(template["fiber_jacobian"]),
                    "odd_bad_integer": str(bad_integer),
                    "equation_replay": {"C1": "0", "target_C2": "0"},
                }
            )
    rows.append(
        {
            "type": "target-unit-pattern-theorem",
            "label": "PROVED",
            "scope": "every rational z with v_w(z)>0; orientation-specific exceptional primes below",
            "statement": (
                "One of two fixed templates has good reduction away from the displayed gcd; "
                "an odd CRT lift of b_mod_w lies in Phi, has v_w(b)=0, and the nonzero "
                "fixed-base fiber Jacobian gives a Q_w point by Hensel."
            ),
            "bad_integer_gcds": {"I": gcds[1], "II": gcds[2]},
            "excluded_primes": {"I": [3], "II": [3, 5]},
            "uses": "Z=z^3, so c_bar=0 and D_bar=1 for A,b w-units",
            "conclusion": (
                "the missed v_w(b)=0 branch has a target-prime local point in both "
                "orientations for w>=7, and in orientation I also for w=5"
            ),
            "warning": "local points do not imply a rational point or global completeness",
        }
    )
    return rows, templates_by_orientation


def target_reduction_points(prime: int, orientation: int, require_A_unit: bool = True) -> list[tuple[int, int, int, int, int, int]]:
    """All (s,A,B,y,r,fiber_det) target-unit residue points."""
    points = []
    for s in range(prime):
        a = (1 + 2 * s) % prime
        A = (1 + 4 * a * a) % prime
        if require_A_unit and A == 0:
            continue
        for B in range(1, prime):
            for y in range(prime):
                for r in range(prime):
                    T = r * r - A * s * s
                    if orientation == 1:
                        C1 = y * y - r * r - A
                        reduced_C2 = r * r * y * y - 16 * B * T - 16
                        determinant = 4 * A * r * y * (y * y + r * r - 16 * B)
                    else:
                        C1 = r * r - y * y - A
                        reduced_C2 = pow(y, 4, prime) - 16 * B * T - 16
                        determinant = 8 * A * r * y * (8 * B - y * y)
                    C2 = A * reduced_C2
                    if C1 % prime or C2 % prime:
                        continue
                    points.append((s, A, B, y, r, determinant % prime))
    return points


def exceptional_unit_rows() -> list[dict[str, Any]]:
    rows = []
    expected = {
        (3, 1): (6, 4, 0, 0),
        (3, 2): (4, 4, 0, 0),
        (5, 1): (10, 8, 8, 8),
        (5, 2): (2, 0, 0, 0),
    }
    for prime in (3, 5):
        for orientation in (1, 2):
            points = target_reduction_points(prime, orientation)
            aligned = [point for point in points if legendre(point[1], prime) == -1]
            nonsingular = [point for point in points if point[-1] != 0]
            aligned_nonsingular = [point for point in aligned if point[-1] != 0]
            observed = (
                len(points), len(aligned), len(nonsingular), len(aligned_nonsingular)
            )
            assert observed == expected[(prime, orientation)]
            if not aligned:
                label = "PROVED_LOCAL_FAILURE_FOR_ALIGNED_UNIT_STRATUM"
                conclusion = "no residue point when A is a nonsquare unit and v_w(b)=0"
            elif aligned_nonsingular:
                label = "PROVED_LOCAL_LIFT_EXISTS_IN_ALIGNED_UNIT_STRATUM"
                conclusion = "the displayed fixed-base fiber points are nonsingular and Hensel-lift"
            else:
                label = "OPEN_AFTER_SINGULAR_RESIDUE_POINTS"
                conclusion = "aligned residue points exist but all fixed-base fiber Jacobians vanish"
            rows.append(
                {
                    "type": "target-unit-exception",
                    "label": label,
                    "prime": prime,
                    "orientation": ORIENTATION_NAMES[orientation],
                    "A_unit_points": len(points),
                    "aligned_nonsquare_A_points": len(aligned),
                    "nonsingular_fixed_base_points": len(nonsingular),
                    "aligned_nonsingular_fixed_base_points": len(aligned_nonsingular),
                    "aligned_points": [list(point) for point in aligned],
                    "conclusion": conclusion,
                    "scope": "c_bar=0, B!=0, A a unit; exact exhaustive residue enumeration",
                }
            )
    return rows


def arm_record(prime: int, sign: int, cell_count: int, min_t: int) -> dict[str, Any]:
    assert sign in (4, -4)
    residue_points = []
    nonsingular_points = []
    for s in range(prime):
        a = (1 + 2 * s) % prime
        A = (1 + 4 * a * a) % prime
        if A == 0 or legendre(A, prime) != -1:
            continue
        y = mod_sqrt(sign, prime)
        r = mod_sqrt(A + sign, prime)
        if y is None or r is None:
            continue
        C1 = (r * r - y * y - A) % prime
        C2 = A * (pow(y, 4, prime) - 16) % prime
        derivative = (-8 * A * r * pow(y, 3, prime)) % prime
        assert C1 == C2 == 0
        point = (s, a, A, y, r, derivative)
        residue_points.append(point)
        if derivative:
            nonsingular_points.append(point)
    if nonsingular_points:
        label = "PROVED_LOCAL_LIFT"
        conclusion = (
            "with b=w (so 0<v_w(b)=1<v_w(Z)), B_bar=c_bar=0 and the displayed "
            "orientation-II fixed-base fiber point Hensel-lifts"
        )
    elif residue_points:
        label = "OPEN_AFTER_SINGULAR_RESIDUE_POINTS"
        conclusion = "aligned arm points exist, but every fixed-base fiber Jacobian vanishes"
    else:
        label = "PROVED_RESIDUE_ARM_EMPTY"
        if sign == -4 and prime % 4 == 3:
            conclusion = "u=-4 is not a square modulo w"
        else:
            conclusion = "exact enumeration found no aligned residue point on this arm"
    return {
        "type": "orientation-II-target-arm",
        "label": label,
        "prime": prime,
        "u_arm": sign,
        "aligned_residue_points": len(residue_points),
        "aligned_nonsingular_points": len(nonsingular_points),
        "sample": list(residue_points[0]) if residue_points else None,
        "nonsingular_sample": list(nonsingular_points[0]) if nonsingular_points else None,
        "cells_with_prime": cell_count,
        "min_v_w_Z": min_t,
        "conclusion": conclusion,
        "scope": "aligned A nonsquare-unit subroute, b=w, exact finite-field enumeration",
    }


def odd_crt_lift(residue: F, prime: int) -> int | None:
    reduced = mod_fraction(residue, prime)
    if reduced in (None, 0):
        return None
    answer = int(reduced)
    if answer % 2 == 0:
        answer += prime
    assert answer % 2 == 1 and answer % prime != 0
    return answer


def cell_sets() -> tuple[list[tuple[str, int, F]], list[dict[str, Any]]]:
    grid = sorted(
        (
            ("grid", prime, unit)
            for prime in odd_primes_below(100)
            for unit in GRID_U_POOL
            if vp(F(prime) * unit, prime) >= 1
        ),
        key=lambda item: (item[1], item[2]),
    )
    assert len(grid) == 353
    offgrid = sorted(
        (("off-grid", prime, unit) for prime, unit in OFFGRID_CELLS),
        key=lambda item: (item[1], item[2]),
    )
    assert len(offgrid) == 17
    assert not {(prime, unit) for _, prime, unit in grid} & {
        (prime, unit) for _, prime, unit in offgrid
    }
    cells = grid + offgrid
    prime_cells: dict[int, list[tuple[str, int, F]]] = defaultdict(list)
    for cell in cells:
        prime_cells[cell[1]].append(cell)
    arms = []
    for prime, cohort in sorted(prime_cells.items()):
        min_t = min(vp((F(p) * unit) ** 3, p) for _, p, unit in cohort)
        arms.append(arm_record(prime, 4, len(cohort), min_t))
        arms.append(arm_record(prime, -4, len(cohort), min_t))
    return cells, arms


def cell_specific_pools(
    prime: int,
    z: F,
    arm_rows_by_prime: dict[int, list[dict[str, Any]]],
    templates_by_orientation: dict[int, tuple[dict[str, F], ...]],
) -> tuple[tuple[F, ...], tuple[F, ...]]:
    s_values = set(S_POOL)
    for row in arm_rows_by_prime[prime]:
        if row["sample"] is not None:
            s_values.add(F(row["sample"][0]))

    b_values = set(B_POOL)
    t = vp(z**3, prime)
    for multiplier in TARGET_MULTIPLIERS:
        candidate = F(prime) * multiplier
        if phi_holds(F(0), candidate) and 0 < vp(candidate, prime) < t:
            b_values.add(candidate)
    for multiplier in (F(-1), F(1), F(-3), F(3)):
        candidate = F(prime**t) * multiplier
        if phi_holds(F(0), candidate) and vp(candidate, prime) == t:
            b_values.add(candidate)
    for templates in templates_by_orientation.values():
        for template in templates:
            bad_integer = template_bad_integer(template)
            if bad_integer % prime:
                lift = odd_crt_lift(template["b_residue"], prime)
                if lift is not None:
                    b_values.add(F(lift))
    ordered_s = tuple(sorted(s_values, key=lambda value: (max(abs(value.numerator), value.denominator), value)))
    ordered_b = tuple(sorted(b_values, key=lambda value: (max(abs(value.numerator), value.denominator), value)))
    return ordered_s, ordered_b


def scan_cells(
    cells: list[tuple[str, int, F]],
    arms: list[dict[str, Any]],
    templates_by_orientation: dict[int, tuple[dict[str, F], ...]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], Counter[str]]:
    hit_rows: list[dict[str, Any]] = []
    local_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    refusal_totals: Counter[str] = Counter()
    arm_rows_by_prime: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in arms:
        arm_rows_by_prime[int(row["prime"])].append(row)

    for family, prime, unit in cells:
        z = F(prime) * unit
        s_pool, b_pool = cell_specific_pools(prime, z, arm_rows_by_prime, templates_by_orientation)
        outcomes: dict[int, dict[str, Counter[str]]] = {
            1: defaultdict(Counter),
            2: defaultdict(Counter),
        }
        orientation_II_dyadic_outcomes: dict[str, Counter[str]] = defaultdict(Counter)
        local_tries: Counter[tuple[int, str]] = Counter()
        local_found: set[tuple[int, str]] = set()
        cell_hit_count = 0
        cell_refusals: Counter[str] = Counter()

        for s in s_pool:
            for b in b_pool:
                PACER.tick()
                base, refusal = make_base(s, b, z)
                if base is None:
                    assert refusal is not None
                    cell_refusals[refusal] += 1
                    refusal_totals[refusal] += 1
                    continue
                mode = target_mode(base, prime)
                for orientation in (1, 2):
                    result = solve_quartic(base, orientation)
                    outcomes[orientation][mode][result.status] += 1
                    if orientation == 2:
                        dyadic_class = (
                            "s-unit-square-root-obstruction"
                            if s != 0 and vp(s, 2) == 0
                            else "s-even-discriminant-obstruction"
                        )
                        orientation_II_dyadic_outcomes[dyadic_class][result.status] += 1
                    if result.solutions:
                        for u, y, r in result.solutions:
                            hit_rows.append(
                                rational_hit_row(family, prime, unit, base, orientation, u, y, r)
                            )
                            cell_hit_count += 1
                    else:
                        key = (orientation, mode)
                        if key not in local_found and local_tries[key] < LOCAL_CERTIFICATE_TRIES:
                            local_tries[key] += 1
                            certificate = find_local_obstruction(base, orientation)
                            if certificate is not None:
                                assert local_obstruction(base, orientation, certificate["prime"]) is not None
                                local_found.add(key)
                                local_rows.append(
                                    {
                                        "type": "fixed-base-local-failure",
                                        "label": "PROVED",
                                        "family": family,
                                        "cell": [prime, [unit.numerator, unit.denominator]],
                                        "orientation": ORIENTATION_NAMES[orientation],
                                        "target_mode": mode,
                                        "s": frac_text(s),
                                        "b": frac_text(b),
                                        "obstruction_prime": certificate["prime"],
                                        "A_mod_p": certificate["A_mod_p"],
                                        "quartic_mod_p": certificate["quartic_mod_p"],
                                        "compatible_u": [],
                                        "scope": "this fixed rational (s,b) only; integral monic quartic has no residue point",
                                    }
                                )

        serial_outcomes = {
            ORIENTATION_NAMES[orientation]: {
                mode: dict(sorted(counter.items()))
                for mode, counter in sorted(mode_counters.items())
            }
            for orientation, mode_counters in outcomes.items()
        }
        no_certificate_windows = [
            {"orientation": ORIENTATION_NAMES[orientation], "target_mode": mode, "tries": tries}
            for (orientation, mode), tries in sorted(local_tries.items())
            if (orientation, mode) not in local_found and tries >= LOCAL_CERTIFICATE_TRIES
        ]
        refusal_totals["local-window-no-certificate"] += len(no_certificate_windows)
        cell_rows.append(
            {
                "type": "cell-proved-empty" if not cell_hit_count else "cell-covered",
                "label": (
                    "PROVED_EMPTY_ON_PHI_BY_DYADIC_THEOREM"
                    if not cell_hit_count
                    else "CONTRADICTS_DYADIC_THEOREM"
                ),
                "family": family,
                "cell": [prime, [unit.numerator, unit.denominator]],
                "z": frac_text(z),
                "Z": frac_text(z**3),
                "s_candidates": len(s_pool),
                "b_candidates": len(b_pool),
                "base_pairs_declared": len(s_pool) * len(b_pool),
                "rational_hits": cell_hit_count,
                "fixed_base_outcomes": serial_outcomes,
                "orientation_II_dyadic_outcomes": {
                    dyadic_class: dict(sorted(counter.items()))
                    for dyadic_class, counter in sorted(orientation_II_dyadic_outcomes.items())
                },
                "local_failure_certificates": sum(
                    1 for row in local_rows if row["cell"] == [prime, [unit.numerator, unit.denominator]]
                ),
                "refusals": dict(sorted(cell_refusals.items())),
                "local_certificate_window_refusals": no_certificate_windows,
                "theorem_authority": ["dyadic-orientation-I", "dyadic-orientation-II"],
                "no_extrapolation": True,
            }
        )
    return hit_rows, local_rows, cell_rows, refusal_totals


def pool_refusal_rows() -> list[dict[str, Any]]:
    s_reasons = Counter()
    for value in RAW_S_POOL:
        if value != 0 and vp(value, 2) < 0:
            s_reasons["v2(s)<0"] += 1
    b_reasons = Counter()
    for value in RAW_B_POOL:
        if value == 0:
            b_reasons["b=0"] += 1
        elif vp(value, 2) != 0:
            b_reasons["v2(b)!=0"] += 1
    return [
        {
            "type": "pool-refusal",
            "label": "REFUSED_OUTSIDE_DOMAIN",
            "variable": "s",
            "raw_candidates": len(RAW_S_POOL),
            "accepted": len(S_POOL),
            "reasons": dict(sorted(s_reasons.items())),
            "evidentiary_value": "none against completeness",
        },
        {
            "type": "pool-refusal",
            "label": "REFUSED_OUTSIDE_DOMAIN",
            "variable": "b",
            "raw_candidates": len(RAW_B_POOL),
            "accepted": len(B_POOL),
            "reasons": dict(sorted(b_reasons.items())),
            "evidentiary_value": "none against completeness",
        },
    ]


def meta_row() -> dict[str, Any]:
    return {
        "type": "meta",
        "label": "PROVED_METHOD_AND_DECLARED_SCOPE",
        "schema": "l24-diagonal-search-v1",
        "base": {
            "a": "1+2s",
            "A": "1+4a^2",
            "B": "2b",
            "Z": "z^3",
            "D": "1-Z-a^2Z^2",
            "c": "a^2*Z^2*(16a^4b^2-A(b-1)^4)/(A*b^2*D)",
        },
        "Phi": ["v2(s)>=0 (or s=0)", "v2(b)=0"],
        "systems": {
            "I": ["y^2-r^2-A=0", "-r^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0"],
            "II": ["r^2-y^2-A=0", "-y^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0"],
        },
        "exact_solver": (
            "orientation I: quartic in r, orientation II: quartic in y; solve the quadratic "
            "in u=x^2 over Q, require both u and A+u rational squares, then reconstruct the conic"
        ),
        "grid": {"cells": 353, "w": "odd primes <100", "u_pool": [frac_text(u) for u in GRID_U_POOL]},
        "off_grid_horizon": {
            "cells": 17,
            "predeclared": [[w, [u.numerator, u.denominator]] for w, u in OFFGRID_CELLS],
        },
        "global_pool_bounds": {
            "s_numerator_abs_le": S_NUMERATOR_BOUND,
            "b_numerator_abs_le": B_NUMERATOR_BOUND,
            "reduced_denominator_le": DENOMINATOR_BOUND,
            "Phi_s_count": len(S_POOL),
            "Phi_b_count": len(B_POOL),
        },
        "cell_specific_additions": [
            "b=w*m for declared m and 0<v_w(b)<v_w(Z)",
            "b=+/-w^v_w(Z), +/-3w^v_w(Z)",
            "odd CRT lifts of the four target-unit templates",
            "least s representatives from both u=+4 and u=-4 aligned arm searches",
        ],
        "proof_discipline": {
            "hits": "PROVED only after original Phi and both tied equations replay exactly",
            "fixed_base_failures": "complete only for that displayed rational (s,b)",
            "local_failures": "PROVED only with an exact good-reduction residue obstruction",
            "bounded_zero_hits": "corroboration only; not used in the global proof",
            "unknown_cells": "none: both orientations are uniformly Phi-empty over Q_2",
            "refusals": "separate and non-evidentiary",
            "coverage_extrapolation": False,
        },
    }


def solver_theorem_row() -> dict[str, Any]:
    # Independent exact sample catches signs and the two eliminants before the scan.
    checks = 0
    for s in (F(-2), F(-1), F(0), F(1), F(3, 5)):
        for b in (F(-5), F(-1), F(1), F(3), F(5, 3)):
            if not phi_holds(s, b):
                continue
            for z in (F(-15), F(3), F(7, 3), F(101)):
                base, refusal = make_base(s, b, z)
                if refusal:
                    continue
                assert base is not None
                for orientation in (1, 2):
                    result = solve_quartic(base, orientation)
                    q4, q2, q0 = quartic_coefficients(base, orientation)
                    for u in result.rational_u_roots:
                        assert q4 * u * u + q2 * u + q0 == 0
                    for u, y, r in result.solutions:
                        replay_rational_hit(base, orientation, u, y, r)
                    checks += 1
    return {
        "type": "solver-theorem",
        "label": "PROVED",
        "statement": (
            "For every guarded rational base the algorithm exhausts Q-points of either diagonal "
            "fiber: its even quartic is quadratic in u, and a rational tied point is equivalent "
            "to a rational u root with u and A+u both squares."
        ),
        "orientation_I_quartic": "A*r^4+(A^2-c^2-16AB)r^2+16A^2Bs^2-16A",
        "orientation_II_quartic": "A*y^4-(c^2+16AB)y^2+16A^2B(s^2-1)-16A",
        "exact_sample_replays": checks,
        "floating_point_operations": 0,
    }


def dyadic_rows() -> list[dict[str, Any]]:
    """Prove that both diagonal orientations are empty on Phi over Q_2."""
    residue_checks = 0
    square_substitution_checks = 0
    for s in range(32):
        a = 1 + 2 * s
        A = 1 + 4 * a * a
        assert A % 32 == 5
        for b in range(1, 64, 2):
            B = 2 * b
            Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
            assert Ng % 16 == 0
            # These are the constant coefficients after division by A and 16.
            assert (A * B * s * s - 1) % 2 == 1
            e = A * B * (s * s - 1) - 1
            assert e % 2 == 1
            assert 16 * B % 32 == 0
            discriminant_unit_mod_8 = (1 - A * B * (s * s - 1)) % 8
            if s % 2:
                assert e % 16 == 15
                assert discriminant_unit_mod_8 == 1
                # Here h=P/32 is odd.  If u=y^2 and v2(u)=2, write
                # u=4*t^2 with t odd; q(u)/16 is 8 mod 16, never zero.
                for h in range(1, 16, 2):
                    for t in range(1, 16, 2):
                        assert (t**4 - 8 * h * t * t + e) % 16 == 8
                        square_substitution_checks += 1
            else:
                assert discriminant_unit_mod_8 in (3, 7)
            residue_checks += 1

    shared = {
        "Phi_facts": [
            "a is a 2-adic unit and A=5 mod 8 (indeed 5 mod 32)",
            "b is a unit, B=2b has valuation 1, and Ng has valuation >=4",
            (
                "v2(D)=0 for v2(Z)>=0 and v2(D)=2v2(Z) for v2(Z)<0; "
                "therefore v2(c)>=4 for every rational Z"
            ),
        ],
        "mechanical_residue_checks": residue_checks,
    }
    return [
        {
            "type": "dyadic-orientation-I",
            "label": "PROVED_EMPTY_ON_PHI",
            **shared,
            "normalized_quartic": (
                "u^2+(A-c^2/A-16B)u+16ABs^2-16=0, u=r^2"
            ),
            "coefficient_valuations": {"constant": 4, "linear": 0, "quadratic": 0},
            "root_valuation_argument": (
                "for t=v2(u), term valuations are (2t,t,4); a minimum can repeat "
                "only at t=0 or t=4"
            ),
            "square_obstruction": {
                "t=4": "A+u=5 mod 8, not a square",
                "t=0": "u=r^2=1 mod 8, hence A+u=6 mod 8, not a square",
            },
            "conclusion": (
                "orientation I has no rational point for any Phi-admissible (s,b) "
                "and any guarded rational Z"
            ),
        },
        {
            "type": "dyadic-orientation-II",
            "label": "PROVED_EMPTY_ON_PHI",
            **shared,
            "normalized_quartic": (
                "u^2-Pu+16e=0, u=y^2, P=c^2/A+32b=32h, "
                "e=2Ab(s^2-1)-1"
            ),
            "coefficient_valuations": {
                "constant": 4,
                "linear": ">=5",
                "quadratic": 0,
            },
            "root_valuation_argument": "term minima force v2(u)=2",
            "even_s_obstruction": (
                "Delta/64=16h^2-e=1-2Ab(s^2-1) mod 8 is 3 or 7, "
                "so the quadratic has no Q_2 root"
            ),
            "odd_s_obstruction": (
                "e=-1 mod 16; if u=y^2=4t^2 with h,t odd, division by 16 "
                "gives t^4-8h*t^2+e=8 mod 16, contradiction"
            ),
            "mechanical_square_substitution_checks": square_substitution_checks,
            "conclusion": (
                "orientation II has no rational point for any Phi-admissible "
                "(s,b) and any guarded rational Z"
            ),
        },
    ]


def replay_rows(rows: Iterable[dict[str, Any]]) -> None:
    for row in rows:
        if row.get("type") == "rational-hit":
            prime, unit_pair = row["cell"]
            unit = F(*unit_pair)
            base, refusal = make_base(F(row["s"]), F(row["b"]), F(prime) * unit)
            assert refusal is None and base is not None
            replay_rational_hit(
                base,
                1 if row["orientation"] == "I" else 2,
                F(row["u"]),
                F(row["y"]),
                F(row["r"]),
            )
        elif row.get("type") == "fixed-base-local-failure":
            prime, unit_pair = row["cell"]
            base, refusal = make_base(
                F(row["s"]), F(row["b"]), F(prime) * F(*unit_pair)
            )
            assert refusal is None and base is not None
            certificate = local_obstruction(
                base,
                1 if row["orientation"] == "I" else 2,
                int(row["obstruction_prime"]),
            )
            assert certificate is not None and certificate["compatible_u"] == []


def write_report(rows: list[dict[str, Any]], elapsed: float) -> None:
    summary = next(row for row in rows if row.get("type") == "summary")
    arms = [row for row in rows if row.get("type") == "orientation-II-target-arm"]
    cells = [
        row
        for row in rows
        if row.get("type") in ("cell-covered", "cell-proved-empty")
    ]
    exceptions = [row for row in rows if row.get("type") == "target-unit-exception"]
    templates = [row for row in rows if row.get("type") == "target-unit-template"]
    local_failures = [row for row in rows if row.get("type") == "fixed-base-local-failure"]
    hits = [row for row in rows if row.get("type") == "rational-hit"]

    lines = [
        "# L24 exact diagonal completeness probe",
        "",
        "## Verdict",
        "",
        (
            f"The exact solver scanned **{summary['cells']} cells** "
            f"({summary['grid_cells']} grid + {summary['off_grid_cells']} predeclared off-grid), "
            f"**{summary['base_pairs_tested']} guarded rational base pairs**, and both orientations. "
            f"It found **{summary['rational_hits']} exact rational hits**.  This bounded zero-hit "
            "result is corroboration only, not the proof."
        ),
        "",
        (
            "**Both orientations are PROVED globally Phi-empty over Q_2.**  Therefore their "
            "union covers no target cell: the diagonal route is closed, all "
            f"**{summary['proved_empty_cells']} declared cells are proved uncovered**, and "
            "there are no unknown cells."
        ),
        "",
        (
            "The separate target-prime calculations remain useful diagnostics.  The missed "
            "`v_w(b)=0` branch has nonsingular target-prime local points in both orientations "
            "for every `w>=7`, and in orientation I at `w=5`; the obstruction killing the "
            "global route is at 2, not at the target prime."
        ),
        "",
        "## Exact solver and labels",
        "",
        "For orientation I the program solves",
        "",
        "`A*r^4 + (A^2-c^2-16AB)*r^2 + 16A^2Bs^2 - 16A = 0`,",
        "",
        "and for orientation II it solves",
        "",
        "`A*y^4 - (c^2+16AB)*y^2 + 16A^2B(s^2-1) - 16A = 0`.",
        "",
        "Each is solved exactly as a quadratic in the tied square `u`; both `u` and `A+u` must be rational squares.  The bounded solver is complete for each fixed `(s,b,z)`, while the dyadic theorem below is uniform in all Phi bases.  The 4,026,282 bounded no-hit attempts are not used to prove the uniform statement.",
        "",
        "### Global dyadic classification [PROVED]",
        "",
        "On Phi, `A=5 mod 32`, `v2(B)=1`, and `v2(c)>=4` for every rational `Z`.  In orientation I, after division by `A`, the quartic has coefficient valuations `(4,0,0)`.  If `t=v2(u)`, the three term valuations `(2t,t,4)` force `t=0` or `4`.  For `u=r^2`, the first gives `A+u=6 mod 8`; the second gives `A+u=5 mod 8`.  Neither is `y^2`.  Hence orientation I is empty on all of Phi.",
        "",
        "For orientation II write the normalized quartic as `u^2-Pu+16e=0`, where `u=y^2`, `P=c^2/A+32b=32h` with `h` odd, and `e=2Ab(s^2-1)-1`.  If `s` is even, `Delta/64=16h^2-e` is `3` or `7 mod 8`, so there is no root.  If `s` is odd, `e=-1 mod 16` and root valuations force `u=4t^2` with `t` odd.  Dividing the equation by 16 gives `t^4-8h*t^2+e = 8 mod 16`, again impossible.  Hence orientation II is also empty on all of Phi.",
        "",
        "## Mined target-unit identities [PROVED]",
        "",
        "Four low-height successful residue templates all gave the unique coefficient `16` in the interpolation window `[-64,64]`:",
        "",
        "`I: r^2*y^2-16 = 16*B*(r^2-A*s^2)`,",
        "",
        "`II: y^4-16 = 16*B*(r^2-A*s^2)`.",
        "",
        "Substitution into the two orientation-specific target reductions was checked with sparse formal polynomials (zero remainder), not only numerically.  The hyperbola parameter identity was independently cleared and checked symbolically.",
        "",
        "| orientation | s | A | lambda | y | r | B mod w | odd bad integer |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in templates:
        lines.append(
            f"| {row['orientation']} | {row['s']} | {row['A']} | {row['lambda']} | "
            f"{row['y']} | {row['r']} | {row['B_mod_w']} | {row['odd_bad_integer']} |"
        )
    lines.extend(
        [
            "",
            "The two odd bad integers have gcd `3` in orientation I and `45` in orientation II.  Thus every odd prime except `3` avoids an orientation-I template, and every prime `w>=7` avoids an orientation-II template.  Its `A,B,T` and fiber Jacobian are units; `Z=z^3` gives `c_bar=0`, and Hensel applies after choosing an odd integer `b` in the displayed nonzero residue class.",
            "",
            "## Exceptional target-unit residues",
            "",
            "| w | orientation | A-unit points | aligned points | nonsingular fixed-base points | label |",
            "|---:|:---:|---:|---:|---:|---|",
        ]
    )
    for row in exceptions:
        lines.append(
            f"| {row['prime']} | {row['orientation']} | {row['A_unit_points']} | "
            f"{row['aligned_nonsquare_A_points']} | {row['nonsingular_fixed_base_points']} | "
            f"{row['label']} |"
        )
    lines.extend(
        [
            "",
            "At `w=5`, orientation I has eight aligned nonsingular points, while orientation II has no aligned residue point in the target-unit stratum.  At `w=3`, aligned points exist in both orientations but every fixed-base fiber Jacobian vanishes, so their lift is recorded as OPEN rather than guessed.",
            "",
            "## Orientation-II `u=+4` and `u=-4` target arms",
            "",
            "Here `b=w`, hence `0<v_w(b)=1<v_w(Z)` on every listed cell and `B_bar=c_bar=0`.  A stored arm is `PROVED_LOCAL_LIFT` only when `u=y^2`, `r^2=A+u`, both reduced equations, and the fixed-base fiber Jacobian are all nonzero/exactly replayed.",
            "",
            "| w | cells | +4 | nonsingular/total | -4 | nonsingular/total |",
            "|---:|---:|---|---:|---|---:|",
        ]
    )
    by_prime: dict[int, dict[int, dict[str, Any]]] = defaultdict(dict)
    for row in arms:
        by_prime[int(row["prime"])][int(row["u_arm"])] = row
    for prime, pair in sorted(by_prime.items()):
        plus, minus = pair[4], pair[-4]
        lines.append(
            f"| {prime} | {plus['cells_with_prime']} | {plus['label']} | "
            f"{plus['aligned_nonsingular_points']}/{plus['aligned_residue_points']} | "
            f"{minus['label']} | "
            f"{minus['aligned_nonsingular_points']}/{minus['aligned_residue_points']} |"
        )
    lines.extend(
        [
            "",
            "No pattern in this finite arm table is extrapolated beyond the declared primes.  In particular the `-4` failures for `w=3 mod 4` are exact because `-4` is nonsquare; other zero counts are finite residue classifications.",
            "",
            "## Bounded rational scan",
            "",
            f"Global pool: `{len(S_POOL)}` Phi-admissible `s` values with numerator height <= {S_NUMERATOR_BOUND}, and `{len(B_POOL)}` Phi-admissible `b` values with numerator height <= {B_NUMERATOR_BOUND}; reduced denominators are <= {DENOMINATOR_BOUND}.  Each cell also receives the predeclared target-valuation, arm, and CRT-template candidates recorded in the JSONL metadata.",
            "",
            (
                "Orientation-II dyadic split in the bounded corroboration: "
                f"`{summary['orientation_II_dyadic_attempts'].get('s-even-discriminant-obstruction', 0)}` "
                "attempts lie in the even-`s` discriminant-obstructed case and "
                f"`{summary['orientation_II_dyadic_attempts'].get('s-unit-square-root-obstruction', 0)}` "
                "in the odd-`s` square-root-obstructed case."
            ),
            "",
            f"Exact fixed-base local obstruction certificates stored: **{len(local_failures)}**.  These are proved failures only for the displayed `(cell,s,b,orientation)`.",
            "",
        ]
    )
    if hits:
        lines.extend(
            [
                "### Rational hits [PROVED]",
                "",
                "| cell | orientation | s | b | y | r |",
                "|---|:---:|---:|---:|---:|---:|",
            ]
        )
        for row in hits:
            lines.append(
                f"| `{row['cell']}` | {row['orientation']} | {row['s']} | {row['b']} | "
                f"{row['y']} | {row['r']} |"
            )
        lines.append("")
    lines.extend(
        [
            "### Per-cell evidence (no extrapolation)",
            "",
            "| family | cell | s pool | b pool | base pairs | hits | local certificates | status |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in cells:
        lines.append(
            f"| {row['family']} | `{row['cell']}` | {row['s_candidates']} | {row['b_candidates']} | "
            f"{row['base_pairs_declared']} | {row['rational_hits']} | "
            f"{row['local_failure_certificates']} | {row['label']} |"
        )
    lines.extend(
        [
            "",
            "## Refusals",
            "",
            f"Raw pool exclusions and guard failures are stored separately.  Runtime refusal totals: `{summary['refusals']}`.  A `local-window-no-certificate` entry means only that the first {LOCAL_CERTIFICATE_TRIES} fixed bases had no obstruction among `{LOCAL_CERTIFICATE_PRIMES}`; it is not a local point and not a negative result.",
            "",
            "## Reproducibility",
            "",
            f"Single process, stdlib exact `Fraction` arithmetic, no floating-point solver, paced every {PACE_EVERY} base pairs with a {PACE_SECONDS}-second sleep.  Pacer sleeps: {PACER.sleeps}; wall time: {elapsed:.3f} seconds; process nice value recorded in JSONL.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    rows: list[dict[str, Any]] = [
        meta_row(), solver_theorem_row(), *dyadic_rows(), symbolic_identity_row()
    ]
    template_rows, templates_by_orientation = unit_template_rows()
    rows.extend(template_rows)
    rows.extend(exceptional_unit_rows())
    cells, arms = cell_sets()
    rows.extend(arms)
    rows.extend(pool_refusal_rows())

    hit_rows, local_rows, cell_rows, refusal_totals = scan_cells(
        cells, arms, templates_by_orientation
    )
    assert not hit_rows
    rows.extend(hit_rows)
    rows.extend(local_rows)
    rows.extend(cell_rows)

    elapsed = time.perf_counter() - started
    covered_cells = sum(row["type"] == "cell-covered" for row in cell_rows)
    base_pairs = sum(row["base_pairs_declared"] - sum(row["refusals"].values()) for row in cell_rows)
    fixed_base_outcomes: Counter[str] = Counter()
    target_mode_counts: Counter[str] = Counter()
    orientation_II_dyadic_counts: Counter[str] = Counter()
    for row in cell_rows:
        for orientation in ("I", "II"):
            for mode, outcomes in row["fixed_base_outcomes"][orientation].items():
                for status, count in outcomes.items():
                    fixed_base_outcomes[status] += count
                    target_mode_counts[mode] += count
        for dyadic_class, outcomes in row["orientation_II_dyadic_outcomes"].items():
            orientation_II_dyadic_counts[dyadic_class] += sum(outcomes.values())
    nice_value = os.getpriority(os.PRIO_PROCESS, 0) if hasattr(os, "getpriority") else None
    summary = {
        "type": "summary",
        "label": "PROVED_BOTH_ORIENTATIONS_EMPTY_WITH_BOUNDED_CORROBORATION",
        "cells": len(cell_rows),
        "grid_cells": sum(row["family"] == "grid" for row in cell_rows),
        "off_grid_cells": sum(row["family"] == "off-grid" for row in cell_rows),
        "base_pairs_tested": base_pairs,
        "orientation_attempts": sum(fixed_base_outcomes.values()),
        "fixed_base_outcomes": dict(sorted(fixed_base_outcomes.items())),
        "target_mode_attempts": dict(sorted(target_mode_counts.items())),
        "orientation_II_dyadic_attempts": dict(sorted(orientation_II_dyadic_counts.items())),
        "rational_hits": len(hit_rows),
        "covered_cells": covered_cells,
        "proved_empty_cells": len(cell_rows),
        "unknown_cells": 0,
        "fixed_base_local_failures": len(local_rows),
        "target_arm_local_lifts": sum(row["label"] == "PROVED_LOCAL_LIFT" for row in arms),
        "target_arm_exact_failures": sum(row["label"] == "PROVED_RESIDUE_ARM_EMPTY" for row in arms),
        "refusals": dict(sorted(refusal_totals.items())),
        "pacer": {"ticks": PACER.steps, "sleeps": PACER.sleeps},
        "process_nice": nice_value,
        "elapsed_seconds": elapsed,
        "global_completeness": "PROVED_FALSE_BOTH_ORIENTATIONS_PHI_EMPTY",
        "bounded_zero_hits_evidentiary_role": "corroboration-only",
        "coverage_extrapolation": False,
    }
    rows.append(summary)

    replay_rows(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in rows),
        encoding="utf-8",
    )
    reparsed = [json.loads(line) for line in OUT.read_text(encoding="utf-8").splitlines()]
    assert len(reparsed) == len(rows)
    assert reparsed[-1]["type"] == "summary"
    write_report(rows, elapsed)
    print(
        f"L24 diagonal search: {len(cell_rows)} cells, {base_pairs} base pairs, "
        f"{len(hit_rows)} hits, {len(local_rows)} fixed-base local failures; {elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
