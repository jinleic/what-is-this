#!/usr/bin/env python3
"""L30 constant-two quartic cover and rigidity/frontier reductions.

The positive result is a new lower-degree section of the L24--L27 square
branch: fix y=2 and r=sX+rho.  The bridge-specialized cover is even quartic
in the hyperbola parameter, not degree eight.  It is uniformly soluble at 2,
has a fibre-regular aligned target point for every odd w>=5, an exact
strong-Hensel fibre at w=3, and a guarded real stratum.  Global rationality
remains open; the selected w=3 fibre is not uniform in a.

The same producer records two independent scope controls: every rational L11c
square-branch tie has a b-adic generic norm-section obstruction, and the L26
trace-base degree-eight fixed-field bundle still carries a separate rational
reciprocal-lift conic.  Recent average-family theorems are not promoted to this
fixed family.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

from l25_scaled_coupling import (
    bridge_c,
    rational_sqrt,
    squareclass_2adic,
    vp,
)
from l26_reciprocal_tie import (
    bridge_data,
    evaluate,
    multiply,
    reciprocal_polynomials,
    reciprocal_tie,
)

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l30_quartic_frontier.jsonl"
REPORT = Path("/tmp/l30_quartic_frontier.md")


def legendre(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0:
        return 0
    answer = pow(residue, (prime - 1) // 2, prime)
    return -1 if answer == prime - 1 else answer


def primes_below(limit: int) -> list[int]:
    return [
        candidate
        for candidate in range(2, limit)
        if all(
            candidate % divisor
            for divisor in range(2, int(candidate**0.5) + 1)
        )
    ]


def unit_mod_8(value: F) -> int:
    valuation = vp(value, 2)
    unit = value / F(2) ** valuation
    assert unit.denominator % 2
    return unit.numerator * pow(unit.denominator, -1, 8) % 8


def hilbert_two(left: F, right: F) -> int:
    left_valuation = vp(left, 2)
    right_valuation = vp(right, 2)
    left_unit = unit_mod_8(left)
    right_unit = unit_mod_8(right)
    epsilon = lambda residue: ((residue - 1) // 2) % 2
    omega = lambda residue: ((residue * residue - 1) // 8) % 2
    exponent = (
        epsilon(left_unit) * epsilon(right_unit)
        + left_valuation * omega(right_unit)
        + right_valuation * omega(left_unit)
    )
    return -1 if exponent % 2 else 1


def constant_two_data(
    a: F, b: F, Z: F, lam: F
) -> tuple[F, F, F, F, F, F, F, F]:
    """Return A,c,X,rho,C1,C2,m,H for y=2 and r=sX+rho."""
    A = 1 + 4 * a * a
    B = 2 * b
    s = (a - 1) / 2
    c = bridge_c(a, Z, b)
    X = (lam + A / lam) / 2
    rho = (lam - A / lam) / 2
    y = F(2)
    r = s * X + rho
    C1 = X * X - rho * rho - A
    C2 = (
        -rho * rho * (c * c - A * y * y)
        - 16 * A * B * (r * r - A * s * s)
        - 16 * A
    )
    m = lam * lam
    q = (
        (s + 1) ** 2 * m * m
        - 2 * A * (s * s + 1) * m
        + (s - 1) ** 2 * A * A
    )
    H = -(c * c - 4 * A) * (m - A) ** 2 - 16 * A * B * q - 64 * A * m
    return A, c, X, rho, C1, C2, m, H


def constant_two_coefficients(
    a: F, b: F, Z: F
) -> tuple[F, F, F, F, F]:
    """Return A,c and ascending coefficients of H_2(m)."""
    A = 1 + 4 * a * a
    B = 2 * b
    s = (a - 1) / 2
    c = bridge_c(a, Z, b)
    E = 4 * A - c * c
    leading = E - 16 * A * B * (s + 1) ** 2
    linear = -2 * A * E + 32 * A * A * B * (s * s + 1) - 64 * A
    constant = A * A * E - 16 * A**3 * B * (s - 1) ** 2
    return A, c, constant, linear, leading


def constant_two_identity() -> dict:
    checks = 0
    for a in map(F, (1, 3, 5, -1, F(5, 3))):
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        for Z in (F(-8), F(-1, 8), F(27), F(1, 27)):
            D = 1 - Z - a * a * Z * Z
            if D == 0:
                continue
            for b in (F(-3), F(-1), F(3), F(5, 3)):
                if b in (0, 1):
                    continue
                Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
                if Ng == 0:
                    continue
                A0, c, constant, linear, leading = constant_two_coefficients(
                    a, b, Z
                )
                assert A0 == A
                B = 2 * b
                E = 4 * A - c * c
                J = (
                    4 * A
                    + 16
                    - c * c
                    - 16 * A * B * (s * s + 1)
                    + 16 * A * A * B * B * s * s
                )
                discriminant = linear * linear - 4 * leading * constant
                assert discriminant == 256 * A * A * J
                for lam in (F(1), F(2), F(3, 2), F(-3), A):
                    _, _, _, _, C1, C2, m, H = constant_two_data(
                        a, b, Z, lam
                    )
                    assert C1 == 0
                    assert H == 4 * m * C2
                    assert H == constant + linear * m + leading * m * m
                    checks += 1
    return {
        "type": "constant-two-quartic-identity",
        "label": "PROVED exact identity",
        "section": "y=2, r=s*X+rho, X^2-rho^2=A",
        "lambda_cover": "H_2(lambda^2)=0",
        "lambda_degree": 4,
        "m_degree": 2,
        "identity": "H_2(m)=4*m*C2",
        "q_polynomial": (
            "Q_s(m)=(s+1)^2*m^2-2*A*(s^2+1)*m+(s-1)^2*A^2"
        ),
        "H_polynomial": (
            "-(c^2-4*A)*(m-A)^2-16*A*B*Q_s(m)-64*A*m"
        ),
        "discriminant": "Disc_m(H_2)=256*A^2*J",
        "J": (
            "4*A+16-c^2-16*A*B*(s^2+1)+16*A^2*B^2*s^2"
        ),
        "exact_instances": checks,
    }


def dyadic_hensel() -> dict:
    rows = 0
    minimum_value_valuation = 10**9
    for a_residue in range(1, 64, 2):
        a = F(a_residue)
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        assert A % 32 == 5
        for b_residue in range(1, 32, 2):
            b = F(b_residue)
            for c_multiplier in range(4):
                c = F(16 * c_multiplier)
                B = 2 * b
                leading = 4 * A - c * c - 16 * A * B * (s + 1) ** 2
                assert vp(leading, 2) == 2
                e = c * c / A
                q_at_one = (
                    (s + 1) ** 2
                    - 2 * A * (s * s + 1)
                    + (s - 1) ** 2 * A * A
                )
                q_factor = -16 * a * a * (
                    s * (1 + 2 * a * a) - a * a * (s * s + 1)
                )
                assert q_at_one == q_factor
                h_at_one = (
                    (4 - e) * (1 - A) ** 2
                    - 32 * b * q_at_one
                    - 64
                )
                assert vp(h_at_one, 2) >= 9
                minimum_value_valuation = min(
                    minimum_value_valuation, vp(h_at_one, 2)
                )
                q_prime_at_one = 2 * (
                    (s + 1) ** 2 - A * (s * s + 1)
                )
                normalized_derivative = (
                    (4 - e) * (1 - A) / 16
                    - b * q_prime_at_one
                    - 2
                )
                assert vp(normalized_derivative, 2) == 0
                rows += 1
    return {
        "type": "constant-two-uniform-dyadic-hensel",
        "label": "PROVED on every Phi dyadic stratum",
        "substitution": "m=1+8*t",
        "normalized_function": "g(t)=H_2(1+8*t)/(256*A)",
        "value": "g(0) is even",
        "derivative": "g'(0) is a 2-adic unit",
        "leading_coefficient": "v2([m^2]H_2)=2 on Phi, so the cover is genuinely quartic",
        "reason": (
            "Q_s(1) is divisible by 16, odd a has a^4=1 mod 16, "
            "and c in 16*Z_2 contributes only above the decisive precision"
        ),
        "consequence": (
            "Hensel gives t in Z_2; m=1+8*t is a Q_2-square, so lambda exists"
        ),
        "finite_residue_rows": rows,
        "minimum_replayed_v2_H_over_A_at_m1": minimum_value_valuation,
    }


def target_rows(prime: int) -> list[dict]:
    rows = []
    for a in range(prime):
        A = (1 + 4 * a * a) % prime
        shifted = (A + 4) % prime
        if (
            A
            and shifted
            and legendre(A, prime) == -1
            and legendre(shifted, prime) == 1
        ):
            X = next(
                candidate
                for candidate in range(prime)
                if candidate * candidate % prime == shifted
            )
            rows.append({"a": a, "A": A, "X": X, "rho": 2 % prime})
    return rows


def target_character_theorem() -> dict:
    small = {}
    for prime in (5, 7, 11, 13, 17, 19):
        rows = target_rows(prime)
        assert rows
        sample = rows[0]
        determinant = 16 * sample["A"] * sample["X"] * sample["rho"]
        assert determinant % prime
        small[str(prime)] = {"count": len(rows), "sample": sample}

    # w=3 is singular in the first residue coordinates.  The total-space
    # row explains the geometry, and a=5 supplies an exact fixed fibre by
    # one-variable strong Hensel in m.
    prime = 3
    a, A, X, rho = 2, 2, 0, 1
    assert (X * X - rho * rho - A) % prime == 0
    assert (4 * A * (rho * rho - 4)) % prime == 0
    jacobian = ((-8 * a) % prime, (-2 * rho) % prime, 0, (8 * A * rho) % prime)
    determinant = (jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2]) % prime
    assert determinant

    soluble_a, soluble_b, soluble_z = F(5), F(3), F(3)
    soluble_Z = soluble_z**3
    soluble_A, soluble_c, soluble_constant, soluble_linear, soluble_leading = (
        constant_two_coefficients(soluble_a, soluble_b, soluble_Z)
    )
    value_at_one = soluble_leading + soluble_linear + soluble_constant
    derivative_at_one = 2 * soluble_leading + soluble_linear
    assert vp(soluble_c, 3) == 4
    assert vp(soluble_leading, 3) == 0
    assert vp(value_at_one, 3) == 3
    assert vp(derivative_at_one, 3) == 1
    assert vp(value_at_one, 3) > 2 * vp(derivative_at_one, 3)

    # Fixed-a solubility is not uniform: this nearby guarded fibre is empty.
    fixed_a, fixed_b, fixed_z = F(1), F(3), F(3)
    fixed_Z = fixed_z**3
    fixed_A, fixed_c, fixed_constant, fixed_linear, fixed_leading = (
        constant_two_coefficients(fixed_a, fixed_b, fixed_Z)
    )
    fixed_s = (fixed_a - 1) / 2
    fixed_B = 2 * fixed_b
    fixed_J = (
        4 * fixed_A
        + 16
        - fixed_c * fixed_c
        - 16 * fixed_A * fixed_B * (fixed_s * fixed_s + 1)
        + 16 * fixed_A * fixed_A * fixed_B * fixed_B * fixed_s * fixed_s
    )
    fixed_discriminant = (
        fixed_linear * fixed_linear - 4 * fixed_leading * fixed_constant
    )
    assert fixed_discriminant == 256 * fixed_A * fixed_A * fixed_J
    assert vp(fixed_c, 3) == 4
    assert vp(fixed_leading, 3) == 0
    assert vp(fixed_J, 3) == 1

    # The Weil lower bound is already positive after the two possible
    # zero-value corrections for every prime >= 23.
    for prime in primes_below(1000):
        if prime < 23:
            continue
        bound = (prime - 3 * math.sqrt(prime)) / 4 - 2
        assert bound > 0
        assert target_rows(prime)

    return {
        "type": "constant-two-target-character-theorem",
        "label": (
            "PROVED fibre-regular for every odd w>=5; "
            "w=3 has an exact strong-Hensel fibre"
        ),
        "target_reduction": (
            "c=0 and B=0 mod w give rho^2=4 and X^2=A+4"
        ),
        "desired_characters": (
            "chi(1+4*a^2)=-1 and chi(5+4*a^2)=+1"
        ),
        "large_prime_indicator": (
            "(1-chi(A))*(1+chi(A+4))/4"
        ),
        "large_prime_bound": (
            "N_w >= (w-3*sqrt(w))/4-2 > 0 for w>=23"
        ),
        "weil_input": (
            "|(sum_a chi((4*a^2+1)*(4*a^2+5)))| <= 3*sqrt(w)"
        ),
        "small_fibre_regular_rows": small,
        "w3": {
            "residue_point": {"a": a, "A": A, "X": X, "rho": rho},
            "total_jacobian_columns": ["a", "rho"],
            "total_jacobian_determinant_mod_3": determinant,
            "soluble_fixed_fibre": {
                "base": {"a": "5", "b": "3", "z": "3", "Z": "27"},
                "v3_c": vp(soluble_c, 3),
                "v3_leading": vp(soluble_leading, 3),
                "v3_H_at_m1": vp(value_at_one, 3),
                "v3_Hprime_at_m1": vp(derivative_at_one, 3),
                "consequence": "strong Hensel gives m in 1+9*Z_3, hence m is a square",
            },
            "empty_control_fibre": {
                "base": {"a": "1", "b": "3", "z": "3", "Z": "27"},
                "v3_c": vp(fixed_c, 3),
                "v3_leading": vp(fixed_leading, 3),
                "v3_J": vp(fixed_J, 3),
                "reason": "disc_m(H_2)=256*A^2*J has odd 3-adic valuation",
            },
            "scope": "target existence is fibrewise but not uniform in the selected a",
        },
        "exhaustive_crosscheck": "every odd prime below 1000",
    }


def real_certificate() -> dict:
    a, b, Z = F(1), F(-1), F(1, 8)
    A = 1 + 4 * a * a
    B = 2 * b
    s = (a - 1) / 2
    assert s == 0
    c = bridge_c(a, Z, b)
    denominator = 4 * A - c * c - 16 * A * B
    u = 16 * A / denominator
    assert denominator > 0 and u > 0 and A + u > 0
    assert u * denominator == 16 * A
    return {
        "type": "constant-two-real-certificate",
        "label": "PROVED on one guarded open real stratum",
        "base": {"a": str(a), "b": str(b), "Z": str(Z)},
        "c": str(c),
        "rho_squared": str(u),
        "X_squared": str(A + u),
        "positive_denominator": str(denominator),
        "reason": (
            "for s=0 the equation is rho^2*(4*A-c^2-16*A*B)=16*A"
        ),
    }


def bounded_quartic_probe() -> dict:
    values = sorted(
        {
            F(numerator, denominator)
            for denominator in range(1, 6)
            for numerator in range(-8, 9)
            if numerator
        }
    )
    checked = 0
    hits = []
    for a in map(F, (1, 3, 5, 7)):
        for b in values:
            if b in (0, 1):
                continue
            for z in values:
                Z = z**3
                A = 1 + 4 * a * a
                D = 1 - Z - a * a * Z * Z
                Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
                if D == 0 or Ng == 0:
                    continue
                _, _, constant, linear, leading = constant_two_coefficients(
                    a, b, Z
                )
                checked += 1
                if leading == 0:
                    continue
                square_discriminant = rational_sqrt(
                    linear * linear - 4 * leading * constant
                )
                if square_discriminant is None:
                    continue
                for m in (
                    (-linear + square_discriminant) / (2 * leading),
                    (-linear - square_discriminant) / (2 * leading),
                ):
                    lam = rational_sqrt(m)
                    if lam not in (None, 0):
                        hits.append(
                            {
                                "a": str(a),
                                "b": str(b),
                                "z": str(z),
                                "lambda": str(lam),
                            }
                        )
    assert checked == 13224
    assert hits == []
    return {
        "type": "constant-two-bounded-global-probe",
        "label": "EVIDENCE ONLY; zero hits is not a no-point theorem",
        "a_values": [1, 3, 5, 7],
        "rational_grid": "nonzero numerators |n|<=8, denominators<=5",
        "Z": "z^3",
        "guarded_rows": checked,
        "rational_lambda_hits": hits,
    }


def general_shear_cancellation() -> dict:
    checks = 0
    for a in map(F, (1, 3, 5, -1)):
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        for p, q, r, t in (
            (F(2), F(28), s, F(1)),
            (F(1), F(3), F(2), F(-1)),
            (F(5, 3), F(-2), F(1, 2), F(7, 3)),
        ):
            for lam in (F(1), F(2), F(3, 2), F(-3)):
                m = lam * lam
                X = (lam + A / lam) / 2
                rho = (lam - A / lam) / 2
                u = rho * rho
                V = X * rho
                coefficient_direct = 4 * m * m * (
                    -16 * A * (r * r + t * t) * u
                    - 32 * A * r * t * V
                    + 16 * A * A * (s * s - r * r)
                )
                H = (
                    (r + t) ** 2 * m * m
                    - 2 * A * (t * t - r * r + 2 * s * s) * m
                    + A * A * (r - t) ** 2
                )
                coefficient_factored = -16 * A * m * H
                assert coefficient_direct == coefficient_factored
                discriminant = (
                    (-2 * A * (t * t - r * r + 2 * s * s)) ** 2
                    - 4 * (r + t) ** 2 * A * A * (r - t) ** 2
                )
                assert discriminant == 16 * A * A * s * s * (
                    s * s + t * t - r * r
                )
                checks += 1
    return {
        "type": "general-linear-shear-B-cancellation-no-go",
        "label": "PROVED for B-independent shear entries and s!=0",
        "matrix": "[[p,q],[r,t]] acting on (X,rho)",
        "B_coefficient": "-16*A*B*m*H_L(m)",
        "H_L": (
            "(r+t)^2*m^2-2*A*(t^2-r^2+2*s^2)*m+A^2*(r-t)^2"
        ),
        "discriminant": "16*A^2*s^2*(s^2+t^2-r^2)",
        "obstruction": (
            "every nonzero rational root has m/A a rational square; "
            "m=lambda^2 and A=5 mod 8 make this impossible over Q_2"
        ),
        "degenerate_case": (
            "if r+t=0, the root is m=A*(r/s)^2 and is equally impossible"
        ),
        "loophole": "s=0 leaves a residual conic and is not excluded",
        "scope": (
            "b-dependent entries are not literal B coefficients; the separate "
            "b-adic norm-section theorem covers generic rational sections"
        ),
        "exact_instances": checks,
    }


def square_branch_rigidity() -> dict:
    checks = 0
    for a in map(F, (1, 3, 5, F(5, 3))):
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        for d in (F(1), F(2), a, 2 * a, F(3, 2)):
            ell = (A - d * d) / (2 * d)
            mu = (A + d * d) / (2 * d)
            delta = -ell * ell / A
            assert mu * mu - ell * ell == A
            assert 1 - A * (mu / A) ** 2 == delta
            for Z in (F(-8), F(1, 8), F(27)):
                D = 1 - Z - a * a * Z * Z
                if D == 0:
                    continue
                for b in (F(-3), F(-1), F(3), F(5, 3)):
                    if b in (0, 1):
                        continue
                    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
                    if Ng == 0:
                        continue
                    c = bridge_c(a, Z, b)
                    B = 2 * b
                    M = 16 + ell * ell * c * c / A - 16 * A * B * s * s
                    P = A * A * b**4 * D * D * M
                    expanded = (
                        16 * A * A * D * D * b**4
                        + ell * ell * a**4 * Z**4 * Ng * Ng / A
                        - 32 * A**3 * s * s * D * D * b**5
                    )
                    source_norm = (
                        (ell * a * a * Z * Z * Ng) ** 2
                        + A
                        * (1 - 2 * A * b * s * s)
                        * (4 * A * D * b * b) ** 2
                    )
                    assert P == expanded
                    assert A * P == source_norm
                    checks += 1
    return {
        "type": "square-branch-tie-rigidity",
        "label": "PROVED generic-section obstruction for every rational d(b)",
        "branch": {
            "ell_d": "(A-d^2)/(2*d)",
            "mu_d": "(A+d^2)/(2*d)",
            "delta_d": "-ell_d^2/A",
            "alpha_d": "ell_d^2",
        },
        "target_norm_field": "Q(sqrt(2*b)), independent of d",
        "M_d": "16+(ell_d^2/A)*c^2-16*A*B*s^2",
        "P_d": "A^2*b^4*D^2*M_d",
        "source_norm_identity": (
            "A*P_d=(ell_d*a^2*Z^2*N_g)^2+"
            "A*(1-2*A*b*s^2)*(4*A*D*b^2)^2"
        ),
        "b_adic_proof": (
            "over K0=Q(a,Z), k=ord_b(ell_d)<=0 for every d in K0(b)^x; "
            "P_d has even valuation 2*k and leading unit A times a square. "
            "A norm Y^2-2*b*R^2 of even valuation has square leading unit, "
            "because its two term valuations have opposite parity"
        ),
        "consequence": (
            "no rational d-tie, including poles or nonlinear witness formulas, "
            "gives a generic b-section; specialized members remain open"
        ),
        "branch_symmetry": "d ~ -d ~ A/d ~ -A/d gives the same delta, M, and P",
        "exact_instances": checks,
    }


def marginal_branch_counterexample() -> dict:
    # The b-adic theorem leaves d of 2-adic valuation +/-1 as the only
    # nonunit class not killed by immediate dominance.  The simplest
    # representative d=2*a has ell=1/(4*a), but is not uniformly dyadic.
    a, b, z = F(1), F(-15), F(-9)
    Z = z**3
    A = 1 + 4 * a * a
    B = 2 * b
    s = (a - 1) / 2
    d = 2 * a
    ell = (A - d * d) / (2 * d)
    assert ell == 1 / (4 * a)
    c = bridge_c(a, Z, b)
    M = 16 + ell * ell * c * c / A - 16 * A * B * s * s
    assert vp(a, 2) == vp(b, 2) == 0
    assert vp(c, 2) == 4
    assert vp(M, 2) == 5
    assert unit_mod_8(M) == 3 and unit_mod_8(B) == 1
    assert hilbert_two(M, B) == -1
    return {
        "type": "marginal-square-branch-counterexample",
        "label": "PROVED exact Phi-base obstruction; not a uniform no-member theorem",
        "branch": "d=2*a, ell=1/(4*a)",
        "base": {"a": str(a), "b": str(b), "z": str(z), "Z": str(Z)},
        "c": str(c),
        "M": str(M),
        "v2_c": vp(c, 2),
        "v2_M": vp(M, 2),
        "unit_M_mod_8": unit_mod_8(M),
        "unit_2b_mod_8": unit_mod_8(B),
        "dyadic_symbol": "(M,2*b)_2=-1",
        "consequence": (
            "the sole simple nonunit d-class surviving the coarse valuation "
            "screen is not uniformly dyadically admissible; isolated members remain open"
        ),
    }


def compose_quadratic(polynomial: list[F], constant: F, quadratic: F) -> list[F]:
    base = [constant, F(0), quadratic]
    answer = [F(0)]
    power = [F(1)]
    for coefficient in polynomial:
        if len(answer) < len(power):
            answer.extend([F(0)] * (len(power) - len(answer)))
        for index, value in enumerate(power):
            answer[index] += coefficient * value
        power = multiply(power, base)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def trace_base_route() -> dict:
    identity_rows = 0
    degree_rows = 0
    lift_rows = 0
    for a in map(F, (1, 3, 5, F(5, 3))):
        A = 1 + 4 * a * a
        assert not squareclass_2adic(A)
        for z in (F(1), F(3), F(1, 2), F(5, 3)):
            Z = z * z
            D = 1 - Z - a * a * Z * Z
            if D == 0:
                continue
            delta = -4 * a**4 / A
            P, T = reciprocal_polynomials(a, Z, delta)
            leading = T[-1]
            norm_two_theta_plus_two = 16 * evaluate(T, F(-2)) / leading
            assert squareclass_2adic(norm_two_theta_plus_two / A)
            assert not squareclass_2adic(norm_two_theta_plus_two)
            for b in (F(3), F(-3), F(5, 3), F(7, 5)):
                A0, D0, _, c = bridge_data(a, Z, b)
                eta = reciprocal_tie(Z, D, b)
                M = 16 - delta * c * c - 32 * A * b * eta * eta
                u = b + 1 / b
                assert (A0, D0) == (A, D)
                assert M == evaluate(T, u) / (D * D * A * A)
                assert (u + 2) / b == ((b + 1) / b) ** 2
                identity_rows += 1
            for w in (F(3), F(5), F(7), F(11)):
                g = compose_quadratic(T, F(-2), w)
                assert len(g) == 9
                assert g[-1] == A * (2 * a**4 * Z**2 * w * w) ** 2
                # If (theta+2)/w were square in K_2, its norm would be square;
                # L28 makes that norm class A.  Capell therefore gives degree 8.
                assert squareclass_2adic(
                    (norm_two_theta_plus_two / 16) / (w**4 * A)
                )
                degree_rows += 1
    for w in map(F, (3, 5, 7, 11)):
        for q in (F(1), F(2), F(3, 2), F(5, 3)):
            rho = q + 1 / (w * q)
            lam = w * q - 1 / q
            b = w * q * q
            u = w * rho * rho - 2
            assert lam * lam == w * (w * rho * rho - 4)
            assert u == b + 1 / b
            lift_rows += 1
    return {
        "type": "trace-base-fixed-field-route",
        "label": "PROVED exact descent and exact remaining lift obstruction",
        "hypotheses": (
            "the L26 square pullback Z=z^2; the L29 compressed pullback "
            "is also covered on its v2(Z)>=2 stratum"
        ),
        "descent": (
            "M=T(u)/(D*A)^2 and [2*b]=[2*(u+2)] for u=b+b^-1"
        ),
        "auxiliary_base_change": "u+2=w*rho^2",
        "auxiliary_conic": "U^2-2*w*V^2=T(w*rho^2-2)*W^2",
        "degree": 8,
        "leading_squareclass": "A",
        "irreducibility": (
            "L28 gives Norm((theta+2)/w) in A*Q_2^x2, hence nonsquare; "
            "Capell makes T(w*rho^2-2) irreducible over Q_2 and Q"
        ),
        "separate_reciprocal_lift": "lambda^2=w*(w*rho^2-4)",
        "lift_parameterization": (
            "rho=q+1/(w*q), lambda=w*q-1/q, b=w*q^2"
        ),
        "circularity_guard": (
            "parameterizing the lift first recovers exactly L29's b=w*q^2 "
            "fixed-field slice and its original parity problem"
        ),
        "Mestre_scope": (
            "the auxiliary has eight geometric simple poles and constant "
            "residue, but arXiv:2511.17213's direct criterion fails: its two "
            "diagonal orientations have leading squareclasses -2*w*A and A"
        ),
        "identity_rows": identity_rows,
        "degree_rows": degree_rows,
        "lift_rows": lift_rows,
    }


def recent_sources() -> dict:
    return {
        "type": "recent-source-routing",
        "label": "SOURCE-CHECKED; no average theorem promoted pointwise",
        "sources": [
            {
                "url": "https://arxiv.org/abs/2604.07047",
                "result": (
                    "Frei-Sofos prove a randomness law for an analytic Hilbert "
                    "symbol and a 100% Hasse principle for coefficient boxes"
                ),
                "first_mismatch": (
                    "Theorem 1.14 is L2 over the full binary-form coefficient "
                    "box; the H10 pair is a fixed form on a prime progression"
                ),
            },
            {
                "url": "https://arxiv.org/abs/2511.17213",
                "result": (
                    "Casarotti-Gammelgaard-Massarenti construct multisections "
                    "for degree-eight conic bundles and discuss Mestre's method"
                ),
                "first_mismatch": (
                    "the direct Mestre sufficient criterion fails on the exact "
                    "diagonal trace-base subfamily, and the rational lift remains"
                ),
            },
            {
                "url": "https://arxiv.org/abs/2607.16006",
                "result": (
                    "Hassett-Hernandez study moduli and monodromy of conic "
                    "surfaces over the projective line"
                ),
                "first_mismatch": (
                    "the paper supplies geometric parameter-space structure, "
                    "not a rational-point theorem for this fixed arithmetic fibre"
                ),
            },
            {
                "url": "https://arxiv.org/abs/2607.25287",
                "result": (
                    "Biswas-C-Ramachandran-Samanta kill Brauer-Manin "
                    "obstructions after suitable finite base extensions"
                ),
                "first_mismatch": (
                    "the needed point is over Q itself; extension-degree "
                    "divisibility does not descend a rational member"
                ),
            },
            {
                "url": "https://github.com/lena-ji/conicbundles",
                "result": "Sage/Macaulay2 checks for real conic-bundle geometry",
                "first_mismatch": (
                    "the scripts treat quartic discriminant curves and real "
                    "signature/torsor checks, not this fixed arithmetic norm family"
                ),
            },
        ],
        "fixed_family_first_estimate": (
            "for L23 bad-prime indicators eta_p(t), prove the two-large sector "
            "sum over z<=p1<p2, p1*p2>D_BV is o(X/(log X)^(3/2)); "
            "Frei-Sofos does not supply this prime-weighted line estimate"
        ),
    }


def audit_record() -> dict:
    return {
        "type": "independent-derivation-audit",
        "label": "REVISED THEN VALIDATED",
        "initial_objection": (
            "an empty fixed fibre at (a,b,z)=(1,3,3) was initially "
            "misread as refuting existential target selection at w=3"
        ),
        "resolution": (
            "the objection was retracted; the exact selected fibre "
            "(5,3,3) satisfies strong Hensel, while (1,3,3) remains "
            "a useful fixed-a nonuniformity control"
        ),
        "surviving_scope_fix": (
            "the trace-base Capell theorem is stated on L26's square "
            "pullback (or L29's compressed v2(Z)>=2 stratum), matching L28"
        ),
        "independently_rederived": [
            "H_2(m)=4*m*C2 and Disc_m(H_2)=256*A^2*J",
            "the uniform dyadic Hensel normalization and exact quartic degree",
            "the w>=5 character count and the w=3 strong-Hensel fibre",
            "the b-adic no-section theorem for every rational d(b)",
            "the general B-independent shear factorization",
            "the trace descent, Capell step, and separate reciprocal lift",
            "the d=2*a dyadic counterexample",
        ],
    }


def frontier_record() -> dict:
    return {
        "type": "strict-frontier",
        "label": "NO unconditional quantifier or H10/Q claim",
        "closed": [
            "an even quartic constant-y section replacing L27's degree-eight equation",
            "uniform dyadic solubility and fibre-regular aligned targets for every odd w>=5",
            "an exact strong-Hensel fibre at w=3, with fixed-a nonuniformity exposed",
            "generic rational-section rigidity for every rational square-branch d-tie",
            "the general B-independent linear-shear cancellation no-go",
            "the trace-base degree-eight fixed-field descent and its separate lift equation",
        ],
        "open": [
            "a rational root of the constant-two bridge-specialized quartic",
            "simultaneous compatibility at all remaining controlled places",
            "a globally good reciprocal member or fixed-family dispersion theorem",
            "removing classical Schinzel H or improving the conditional six-count",
            "Hilbert's tenth problem over Q",
        ],
        "rank_change": (
            "the algebraic cover degree falls from 8 to 4; no completeness or "
            "Diophantine-rank conclusion follows without a rational point theorem"
        ),
    }


def build_report(records: list[dict], elapsed: float) -> str:
    lines = [
        "# L30 quartic frontier replay",
        "",
        "The constant-two cover has degree four in lambda, is uniform at 2,",
        "fibre-regular at aligned targets w>=5, and has a w=3 Hensel fibre.",
        "Global rationality remains OPEN.",
        "",
    ]
    for record in records:
        lines.extend(
            [
                "## " + record["type"],
                "",
                "```json",
                json.dumps(record, indent=2, sort_keys=True),
                "```",
                "",
            ]
        )
    lines.append(f"elapsed_seconds: {elapsed:.3f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [
        constant_two_identity(),
        dyadic_hensel(),
        target_character_theorem(),
        real_certificate(),
        bounded_quartic_probe(),
        general_shear_cancellation(),
        square_branch_rigidity(),
        marginal_branch_counterexample(),
        trace_base_route(),
        recent_sources(),
        audit_record(),
        frontier_record(),
    ]
    elapsed = time.perf_counter() - started
    with OUT.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
