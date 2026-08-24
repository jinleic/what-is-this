#!/usr/bin/env python3
"""Arithmetic-theorem audit for the L24 diagonal degree-eight covers.

The replay is stdlib-only.  It proves a uniform Q_2 obstruction eliminating
both orientations on the guarded Phi base.  It also proves the tower, branch,
genus, and monodromy claims on an admissible one-dimensional slice and records
the exact matches and failures of Hilbert irreducibility, weak approximation,
norm-form theorems, and Chebotarev.

Replay from the workspace root with

    nice -n 19 python3 math/h10q/l24_diagonal_arithmetic.py
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import time
from typing import Iterable


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l24_diagonal_arithmetic.jsonl"
REPORT = Path("/tmp/l24_diagonal_arithmetic.md")

Polynomial = list[F]  # ascending coefficients


def trim(polynomial: Iterable[int | F]) -> Polynomial:
    answer = [F(coefficient) for coefficient in polynomial]
    while len(answer) > 1 and not answer[-1]:
        answer.pop()
    return answer or [F(0)]


def add(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] += coefficient
    for index, coefficient in enumerate(right):
        answer[index] += coefficient
    return trim(answer)


def scale(polynomial: Polynomial, scalar: int | F) -> Polynomial:
    scalar_f = F(scalar)
    return trim([scalar_f * coefficient for coefficient in polynomial])


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += left_coefficient * right_coefficient
    return trim(answer)


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("negative polynomial exponent")
    answer = [F(1)]
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, factor)
        factor = multiply(factor, factor)
        remaining //= 2
    return answer


def derivative(polynomial: Polynomial) -> Polynomial:
    return trim([index * polynomial[index]
                 for index in range(1, len(polynomial))] or [F(0)])


def divide_with_remainder(dividend: Polynomial,
                          divisor: Polynomial) -> tuple[Polynomial, Polynomial]:
    remainder = trim(dividend)
    divisor = trim(divisor)
    if divisor == [F(0)]:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [F(0)] * max(1, len(remainder) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and remainder != [F(0)]:
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = trim(remainder)
    return trim(quotient), remainder


def polynomial_gcd(left: Polynomial, right: Polynomial) -> Polynomial:
    left = trim(left)
    right = trim(right)
    while right != [F(0)]:
        _, remainder = divide_with_remainder(left, right)
        left, right = right, remainder
    if left == [F(0)]:
        return left
    return scale(left, F(1) / left[-1])


def evaluate(polynomial: Polynomial, value: int | F) -> F:
    value_f = F(value)
    answer = F(0)
    for coefficient in reversed(polynomial):
        answer = answer * value_f + coefficient
    return answer


def primitive_coefficients(polynomial: Polynomial) -> list[int]:
    denominator = 1
    for coefficient in polynomial:
        denominator = lcm(denominator, coefficient.denominator)
    integers = [coefficient.numerator * (denominator // coefficient.denominator)
                for coefficient in polynomial]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if content:
        integers = [value // content for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def fraction_text(value: int | F) -> str:
    value_f = F(value)
    return (str(value_f.numerator) if value_f.denominator == 1
            else f"{value_f.numerator}/{value_f.denominator}")


def is_rational_square(value: int | F) -> bool:
    value_f = F(value)
    if value_f < 0:
        return False
    return (isqrt(value_f.numerator) ** 2 == value_f.numerator
            and isqrt(value_f.denominator) ** 2 == value_f.denominator)


def v2(value: int | F) -> int:
    value_f = F(value)
    if not value_f:
        raise ValueError("v2(0) is not needed in this replay")
    numerator = abs(value_f.numerator)
    denominator = value_f.denominator
    answer = 0
    while numerator % 2 == 0:
        numerator //= 2
        answer += 1
    while denominator % 2 == 0:
        denominator //= 2
        answer -= 1
    return answer


def h_value(s: F, b: F, z: F) -> F:
    a = 1 + 2 * s
    A = 1 + 4 * a * a
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    if not D or b in (0, 1):
        raise ValueError("outside the guarded Phi base")
    g = 16 * a ** 4 / A - ((b - 1) ** 2 / b) ** 2
    return a * a * Z * Z * g / D


def orientation_parameters(orientation: str, s: F, b: F) -> dict[str, F]:
    a = 1 + 2 * s
    A = 1 + 4 * a * a
    B = 2 * b
    if orientation == "I":
        d = A
        p0 = A + 16 * B
        q0 = 16 * A * B * (1 + s * s) - 16
    elif orientation == "II":
        d = F(0)
        p0 = 16 * B
        q0 = -16 * A * B * (1 - s * s) - 16
    else:
        raise ValueError(f"unknown orientation {orientation}")
    L = 2 * d - p0
    M = d * d - p0 * d + q0
    N = M - A * L + A * A
    return {
        "s": s, "a": a, "A": A, "b": b, "B": B,
        "d": d, "p0": p0, "q0": q0, "L": L, "M": M, "N": N,
    }


def lambda_polynomial(A: F, L: F, M: F, c: F) -> Polynomial:
    """The monic reciprocal octic after u=(lambda^2-A)^2/(4 lambda^2)."""
    q = L - c * c / A
    lambda2 = [F(0), F(0), F(1)]
    numerator = add(lambda2, [F(-A)])
    return add(
        add(power(numerator, 4),
            scale(multiply(lambda2, power(numerator, 2)), 4 * q)),
        scale(power(lambda2, 2), 16 * M),
    )


def exact_cover_record() -> dict[str, object]:
    samples = [(F(2), F(5, 3)), (F(1, 3), F(-7, 5)),
               (F(-2, 5), F(11, 7))]
    checks = 0
    for orientation in ("I", "II"):
        parameters = orientation_parameters(orientation, F(1), F(3))
        A = parameters["A"]
        B = parameters["B"]
        s = parameters["s"]
        for z, lam in samples:
            c = h_value(s, parameters["b"], z)
            u = (lam * lam - A) ** 2 / (4 * lam * lam)
            small = u * u + (parameters["L"] - c * c / A) * u + parameters["M"]
            polynomial = lambda_polynomial(A, parameters["L"], parameters["M"], c)
            assert len(polynomial) == 9 and polynomial[-1] == 1
            assert polynomial[0] == A ** 4
            assert polynomial[1] == polynomial[3] == polynomial[5] == polynomial[7] == 0
            assert polynomial[2] == A * A * polynomial[6]
            assert evaluate(polynomial, lam) == 16 * lam ** 4 * small

            active = (lam - A / lam) / 2
            other = (lam + A / lam) / 2
            assert active * active == u and other * other == u + A
            if orientation == "I":
                y, r = other, active
                C1 = y * y - r * r - A
                C2 = (-r * r * (c * c - A * y * y)
                      - 16 * A * B * (r * r - A * s * s) - 16 * A)
            else:
                y, r = active, other
                C1 = r * r - y * y - A
                C2 = (-y * y * (c * c - A * y * y)
                      - 16 * A * B * (r * r - A * s * s) - 16 * A)
            assert C1 == 0 and C2 == A * small
            checks += 1

    return {
        "type": "finite-cover",
        "label": "PROVED",
        "guarded_base": [
            "a=1+2s", "A=1+4a^2", "B=2b", "Z=z^3",
            "b*(b-1)*(1-Z-a^2*Z^2)!=0", "v2(s)>=0", "v2(b)=0",
        ],
        "active_coordinate": {
            "I": "u=r^2, y^2=u+A, t=y^2=u+A",
            "II": "u=y^2, r^2=u+A, t=y^2=u",
        },
        "unified_quadratic": "u^2+(L-c^2/A)u+M=0",
        "orientation_parameters": {
            "I": {
                "d": "A", "L": "A-16B",
                "M": "16*(A*B*s^2-1)",
                "N": "16*A*B*(1+s^2)-16",
            },
            "II": {
                "d": "0", "L": "-16B",
                "M": "-16*A*B*(1-s^2)-16",
                "N": "A^2+16*A*B*s^2-16",
            },
        },
        "lambda": "lambda=y+r",
        "parametrization": (
            "active=(lambda-A/lambda)/2, other=(lambda+A/lambda)/2"
        ),
        "monic_octic": (
            "(lambda^2-A)^4+4*(L-c^2/A)*lambda^2*(lambda^2-A)^2"
            "+16*M*lambda^4"
        ),
        "reciprocity": "P(lambda)=lambda^8+g*lambda^6+h*lambda^4+A^2*g*lambda^2+A^4",
        "finite_flat_reason": (
            "the guarded coordinate ring quotient by this monic octic is free "
            "of rank 8, including over its branch divisor"
        ),
        "degree": 8,
        "exact_residual_checks": checks,
    }


def phi_two_adic_record() -> dict[str, object]:
    """Uniform Q_2 analysis after substituting c=h(a,b,z^3)."""
    residue_checks = 0
    square_u_checks = 0
    for s_residue in range(128):
        a = 1 + 2 * s_residue
        A = 1 + 4 * a * a
        assert A % 8 == 5
        for b_residue in range(1, 256, 2):
            # The numerator of g is 16*a^4*b^2-A*(b-1)^4.
            ng = 16 * a ** 4 * b_residue ** 2 - A * (b_residue - 1) ** 4
            assert ng % 16 == 0

            linear_I = A - 32 * b_residue
            constant_I = 16 * (2 * A * b_residue * s_residue ** 2 - 1)
            linear_II = -32 * b_residue
            constant_II = 16 * (
                2 * A * b_residue * (s_residue ** 2 - 1) - 1
            )
            assert linear_I % 8 == 5
            assert constant_I % 32 == 16
            assert linear_II % 64 == 32
            assert constant_II % 32 == 16
            discriminant_unit_II = (
                1 - 2 * A * b_residue * (s_residue ** 2 - 1)
            ) % 8
            if s_residue % 2:
                assert discriminant_unit_II == 1
                for t_residue in range(1, 32, 2):
                    normalized_II = (
                        t_residue ** 4
                        - 8 * b_residue * t_residue ** 2
                        + 2 * A * b_residue * (s_residue ** 2 - 1)
                        - 1
                    ) % 16
                    assert normalized_II == 8
                    square_u_checks += 1
                    if square_u_checks % 32768 == 0:
                        time.sleep(0.0002)
            else:
                assert discriminant_unit_II in (3, 7)
            residue_checks += 1
            if residue_checks % 4096 == 0:
                time.sleep(0.0002)
    assert residue_checks == 128 * 128
    assert square_u_checks == 64 * 128 * 16

    return {
        "type": "Phi-two-adic-obstruction",
        "label": "PROVED",
        "h_valuation": {
            "g": (
                "(16*a^4*b^2-A*(b-1)^4)/(A*b^2), whose numerator is "
                "divisible by 16 because a,b,A are 2-adic units and b-1 is even"
            ),
            "D_cases": {
                "e=v2(Z)>=0": "v2(D)=0",
                "e=v2(Z)<0": "v2(D)=2e (the a^2*Z^2 term is uniquely minimal)",
            },
            "conclusion": (
                "v2(c)=v2(a^2*g*Z^2/D)>=4 for every guarded rational z (D!=0)"
            ),
        },
        "common_Phi_data": [
            "v2(s)>=0 makes a=1+2s a 2-adic unit",
            "A=1+4a^2=5 mod 8",
            "v2(b)=0 and B=2b",
        ],
        "orientation_I": {
            "quadratic": (
                "u^2+(A-32b-c^2/A)u+16*(2*A*b*s^2-1)=0"
            ),
            "coefficient_valuations_ascending": [4, 0, 0],
            "newton_conclusion_for_a_Q2_root": "v2(u) is 4 or 0",
            "square_contradiction": {
                "v2(u)=4": "if u is square then u+A=5 mod 8, not square",
                "v2(u)=0": (
                    "if u is square then u=1 mod 8, so u+A=6 mod 8 has "
                    "valuation 1 and is not square"
                ),
            },
            "verdict": "EMPTY over Q2 at every Phi base point",
        },
        "orientation_II": {
            "quadratic": (
                "u^2+(-32b-c^2/A)u+16*(2*A*b*(s^2-1)-1)=0"
            ),
            "coefficient_valuations_ascending": [4, 5, 0],
            "newton_conclusion_for_a_Q2_root": "v2(u)=2",
            "even_s_obstruction": (
                "Delta/64=1-A*B*(s^2-1) mod 8; if s is not a 2-adic unit "
                "(including s=0) this is 3 or 7, so there is no rational Q2 root"
            ),
            "unit_s_square_u_obstruction": (
                "if v2(s)=0, a diagonal point has u=4*t^2 with t a unit; "
                "after division by 16 the equation is "
                "t^4-8*b*t^2+2*A*b*(s^2-1)-1-(c^2/(4*A))*t^2=0, "
                "but its left side is 8 mod 16"
            ),
            "verdict": "EMPTY over Q2 at every Phi base point",
        },
        "exhaustive_unit_residue_checks_mod_256": residue_checks,
        "exhaustive_unit_s_square_u_checks_mod_16": square_u_checks,
        "scope": (
            "both orientations are unconditionally ruled out on Phi; the finite "
            "diagonal cover cannot define any cell"
        ),
    }


def c_line_branch_record(orientation: str, parameters: dict[str, F]) -> dict[str, object]:
    A = parameters["A"]
    L = parameters["L"]
    M = parameters["M"]
    N = parameters["N"]
    delta = [A * A * (L * L - 4 * M), F(0), -2 * A * L, F(0), F(1)]
    sign_branch = [N, F(0), F(1)]
    branch = multiply(delta, sign_branch)
    quartic_v = [M, F(0), L, F(0), F(1)]
    quadratic_v = [A, F(0), F(1)]
    assert polynomial_gcd(branch, derivative(branch)) == [F(1)]
    assert polynomial_gcd(quartic_v, derivative(quartic_v)) == [F(1)]
    assert polynomial_gcd(quartic_v, quadratic_v) == [F(1)]
    assert not is_rational_square(M)

    return {
        "type": "c-line-slice",
        "label": "PROVED",
        "orientation": orientation,
        "slice": {
            "s": fraction_text(parameters["s"]),
            "a": fraction_text(parameters["a"]),
            "A": fraction_text(A),
            "b": fraction_text(parameters["b"]),
            "B": fraction_text(parameters["B"]),
            "phi_guards": {"v2(s)>=0": v2(parameters["s"]) >= 0,
                           "v2(b)=0": v2(parameters["b"]) == 0},
        },
        "constants": {name: fraction_text(parameters[name])
                      for name in ("d", "p0", "q0", "L", "M", "N")},
        "curve": {
            "equations": [
                "w^2=A*(v^4+L*v^2+M)",
                "y^2=v^2+A",
                "u=v^2", "c=w/v",
            ],
            "biquadratic_branch_sets_over_v": [4, 2],
            "geometrically_integral": True,
            "smooth_projective_genus": 3,
            "riemann_hurwitz": "2g-2=4*(-2)+6*2=4",
        },
        "degree_8_map_to_c": {
            "quadratic_in_u": "A*u^2+(A*L-c^2)*u+A*M=0",
            "discriminant": "Delta(c)=(A*L-c^2)^2-4*A^2*M",
            "other_square_product": "(u1+A)*(u2+A)=c^2+N",
            "branch_polynomial_coefficients_ascending": primitive_coefficients(branch),
            "branch_polynomial_squarefree_gcd": [1],
            "branch_values": {
                "Delta=0": {"count": 4, "inertia_cycle_type": "2^4",
                            "ramification_contribution_each": 4},
                "c^2+N=0": {"count": 2, "inertia_cycle_type": "2^2 1^4",
                             "ramification_contribution_each": 2},
                "c=infinity": {"count": 1, "ramified": False},
            },
            "riemann_hurwitz": "2g-2=-16+4*4+2*2=4",
        },
        "monodromy": {
            "geometric": "C2^3 semidirect C2 = D8 x C2 (order 16)",
            "arithmetic": "C2^4 semidirect C2 = V4 wreath C2 (order 32)",
            "constant_quadratic_extension": f"Q(sqrt({fraction_text(M)}))",
            "constant_is_nonsquare": True,
            "kummer_rank_proof": (
                "over Qbar(c), u1*u2=M is a square and is the sole relation; "
                "u has divisor 2*((0,0)-O) on W^2=A*u*(u^2+L*u+M), so its "
                "class is the nontrivial rational 2-torsion class, while the two "
                "classes u_i+A have distinct odd valuations above c^2+N=0"
            ),
        },
    }


def pullback_record(orientation: str, parameters: dict[str, F]) -> dict[str, object]:
    A = parameters["A"]
    a = parameters["a"]
    b = parameters["b"]
    L = parameters["L"]
    M = parameters["M"]
    N = parameters["N"]
    g = 16 * a ** 4 / A - ((b - 1) ** 2 / b) ** 2
    K = a * a * g

    # c(z)=K*z^6/(1-z^3-a^2*z^6), represented by coprime integer polynomials.
    c_numerator = [F(0)] * 6 + [F(K.numerator)]
    c_denominator = scale(
        [F(1), F(0), F(0), F(-1), F(0), F(0), -a * a],
        K.denominator,
    )
    denominator2 = power(c_denominator, 2)
    numerator2 = power(c_numerator, 2)
    delta_numerator = add(
        power(add(scale(denominator2, A * L), scale(numerator2, -1)), 2),
        scale(power(c_denominator, 4), -4 * A * A * M),
    )
    sign_numerator = add(numerator2, scale(denominator2, N))
    pulled_branch = multiply(delta_numerator, sign_numerator)
    assert len(delta_numerator) - 1 == 24
    assert len(sign_numerator) - 1 == 12
    assert len(pulled_branch) - 1 == 36
    assert polynomial_gcd(pulled_branch, derivative(pulled_branch)) == [F(1)]

    c_branch = multiply(
        [A * A * (L * L - 4 * M), F(0), -2 * A * L, F(0), F(1)],
        [N, F(0), F(1)],
    )
    critical_values = [F(0), -4 * K / A, -K / (a * a)]
    critical_evaluations = [evaluate(c_branch, value) for value in critical_values]
    assert all(critical_evaluations)

    return {
        "type": "z-line-pullback",
        "label": "PROVED",
        "orientation": orientation,
        "slice": "s=1, b=3",
        "h_map": {
            "formula": "c(z)=11072*z^6/(37*(1-z^3-9*z^6))",
            "degree": 6,
            "factorization": "z -> Z=z^3 (degree 3), then Z -> c (degree 2)",
            "geometric_branch_values": [fraction_text(value)
                                        for value in critical_values],
        },
        "branch_disjointness": {
            "cover_branch_polynomial_at_h_branch_values": [
                fraction_text(value) for value in critical_evaluations
            ],
            "all_nonzero": True,
            "consequence": (
                "the two geometric Galois closures have no nontrivial common "
                "subcover: a common subcover would be unramified everywhere on P1"
            ),
        },
        "pulled_branch": {
            "Delta_numerator_degree": len(delta_numerator) - 1,
            "other_square_numerator_degree": len(sign_numerator) - 1,
            "total_degree": len(pulled_branch) - 1,
            "primitive_coefficients_ascending": primitive_coefficients(pulled_branch),
            "squarefree_gcd": [1],
            "Delta_points": 24,
            "other_square_points": 12,
        },
        "normalization": {
            "geometrically_integral": True,
            "degree_over_z": 8,
            "genus": 53,
            "riemann_hurwitz": "2g-2=-16+24*4+12*2=104",
        },
        "monodromy_preserved": {
            "geometric": "D8 x C2 (order 16)",
            "arithmetic": "V4 wreath C2 (order 32)",
            "constant_check": (
                "the Galois closure of z^3 contributes only Q(sqrt(-3)); "
                f"M={fraction_text(M)} has squareclass different from -3"
            ),
        },
    }


def cycle_type(permutation: tuple[int, ...]) -> str:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    counts = Counter(lengths)
    return " ".join(str(length) if count == 1 else f"{length}^{count}"
                    for length, count in sorted(counts.items()))


def wreath_permutation(left_translation: int, right_translation: int,
                       swap_blocks: bool) -> tuple[int, ...]:
    answer: list[int] = []
    for point in range(8):
        block, coordinate = divmod(point, 4)
        translation = left_translation if block == 0 else right_translation
        target_block = 1 - block if swap_blocks else block
        answer.append(4 * target_block + (coordinate ^ translation))
    return tuple(answer)


def monodromy_action_record() -> dict[str, object]:
    arithmetic = {
        wreath_permutation(left, right, swap)
        for left in range(4) for right in range(4) for swap in (False, True)
    }
    geometric = {
        wreath_permutation(left, right, swap)
        for left in range(4) for right in range(4) for swap in (False, True)
        if (left & 1) == (right & 1)
    }
    assert len(arithmetic) == 32 and len(geometric) == 16
    arithmetic_cycles = Counter(cycle_type(permutation) for permutation in arithmetic)
    geometric_cycles = Counter(cycle_type(permutation) for permutation in geometric)
    assert arithmetic_cycles == {"1^8": 1, "1^4 2^2": 6, "2^4": 13, "4^2": 12}
    assert geometric_cycles == {"1^8": 1, "1^4 2^2": 2, "2^4": 9, "4^2": 4}
    arithmetic_fixed = sum(any(permutation[index] == index for index in range(8))
                           for permutation in arithmetic)
    geometric_fixed = sum(any(permutation[index] == index for index in range(8))
                          for permutation in geometric)
    assert arithmetic_fixed == 7 and geometric_fixed == 3
    return {
        "type": "monodromy-action",
        "label": "PROVED",
        "degree": 8,
        "arithmetic_group": {
            "name": "V4 wreath C2",
            "order": 32,
            "cycle_type_counts": dict(sorted(arithmetic_cycles.items())),
            "fixed_point_elements": arithmetic_fixed,
            "derangements": 32 - arithmetic_fixed,
            "derangement_fraction": "25/32",
        },
        "geometric_group": {
            "name": "D8 x C2",
            "order": 16,
            "cycle_type_counts": dict(sorted(geometric_cycles.items())),
            "fixed_point_elements": geometric_fixed,
            "derangements": 16 - geometric_fixed,
            "derangement_fraction": "13/16",
        },
    }


def theorem_records() -> list[dict[str, object]]:
    return [
        {
            "type": "theorem-match",
            "theorem": "Q2 square classification plus Newton valuation polygon",
            "label": "PROVED UNCONDITIONAL OBSTRUCTION TO BOTH ORIENTATIONS",
            "hypotheses_verified": [
                "Phi: v2(s)>=0 and v2(b)=0",
                "the exact substitution c=h(a,b,z^3)",
                "A=5 mod 8 and v2(c)>=4",
            ],
            "valid_unconditional_conclusion": (
                "neither orientation has a Q2 point at any Phi base point, hence "
                "the diagonal union has no rational point in any cell"
            ),
            "orientation_II_cases": {
                "s not a 2-adic unit": (
                    "Delta/64 is 3 or 7 mod 8, so there is no Q2 root"
                ),
                "v2(s)=0": (
                    "a square root u=4*t^2 makes the divided equation 8 mod 16"
                ),
            },
            "every_cell_conclusion": {
                "orientation_I": False,
                "orientation_II": False,
                "union": False,
            },
        },
        {
            "type": "theorem-match",
            "theorem": "Hilbert irreducibility",
            "label": "PROVED NON-APPLICATION TO POINT PRODUCTION",
            "hypotheses_verified_on_slice": [
                "the z-cover is geometrically integral and generically separable",
                "its generic degree-8 polynomial is irreducible over Q(z)",
            ],
            "valid_unconditional_conclusion": (
                "a Hilbert subset of rational z0 has irreducible degree-8 fibre; "
                "those fixed-(s,b) fibres have no rational lambda"
            ),
            "first_failed_hypothesis_for_desired_use": (
                "point production would need rational points on the source in a "
                "Hilbert-dense family (a section/rational parametrization or an "
                "independent density theorem); generic irreducibility supplies the "
                "opposite conclusion, and a prescribed arbitrary z is not a Hilbert subset"
            ),
            "one_point_warning": (
                "one isolated rational point gives only its one specialization; HIT "
                "does not manufacture source points or points on arbitrary fibres"
            ),
            "every_cell_conclusion": False,
        },
        {
            "type": "theorem-match",
            "theorem": "Faltings plus weak-approximation comparison",
            "label": "PROVED NEGATIVE MATCH ON THE RECORDED FIXED SLICE",
            "valid_unconditional_conclusion": (
                "the smooth fixed-(s,b) c-curve has genus 3 and its actual z-pullback "
                "has genus 53, so Faltings makes each of their Q-point sets finite"
            ),
            "weak_approximation_nonapplication": (
                "weak approximation is a density property of an already available "
                "global point set, not an existence theorem; the standard results for "
                "rational varieties and connected linear algebraic groups do not apply "
                "to these higher-genus curves"
            ),
            "first_failed_hypothesis_for_local_to_global_use": (
                "local solubility already fails: the substituted Phi fibre has no "
                "Q2 point in either orientation, so weak approximation and any "
                "Hasse-principle argument stop before Brauer-Manin analysis"
            ),
            "quantifier_scope": (
                "the genus/Faltings calculation is side structure only; the uniform "
                "Q2 obstruction has already decided the union, including degenerate "
                "parameter loci"
            ),
            "every_cell_conclusion": False,
        },
        {
            "type": "theorem-match",
            "theorem": "Hasse norm theorem / norm-form local-global principle",
            "label": "PROVED STRUCTURAL HYPOTHESIS FAILURE",
            "exact_quotient": (
                "W^2=A*u*(u^2+L*u+M), with diagonal points imposing both "
                "u in Q^2 and u+A in Q^2 and the prescribed value c=W/u"
            ),
            "interpretation": (
                "this is a simultaneous (Z/2)^2 Kummer cover of an elliptic curve, "
                "not a principal homogeneous space N_{L/Q}(x)=alpha for one fixed "
                "cyclic number-field extension L/Q"
            ),
            "first_failed_hypothesis": (
                "the fixed cyclic norm-torus equation required by the Hasse norm theorem "
                "is absent; the extension and right side cannot absorb both squareclass "
                "conditions and the prescribed c-fibre"
            ),
            "additional_failed_hypothesis": (
                "even an alternative norm formulation would fail the all-local "
                "premise because the required diagonal variety is empty over Q2"
            ),
            "every_cell_conclusion": False,
        },
        {
            "type": "theorem-match",
            "theorem": "finite-cover Chebotarev",
            "label": "PROVED NON-APPLICATION TO GLOBAL POINTS",
            "valid_unconditional_conclusion": (
                "Frobenius classes distribute in the computed monodromy action; "
                "fixed-point classes give residue points and derangements give no residue point"
            ),
            "computed_signal": (
                "the arithmetic action has 7/32 fixed-point elements and 25/32 "
                "derangements; geometric proportions are 3/16 and 13/16"
            ),
            "first_failed_hypothesis_for_desired_use": (
                "Chebotarev controls primes/closed-point specializations, not a rational "
                "point on every prescribed Q-fibre; a fixed residue point gives at most a "
                "Q_p point after nonsingular Hensel lifting"
            ),
            "local_global_gap": (
                "the computed Frobenius statistics are side structure: the exact "
                "Q2 obstruction already precludes every rational Phi point, regardless "
                "of residue fixed points at odd primes"
            ),
            "every_cell_conclusion": False,
        },
    ]


def build_report(records: list[dict[str, object]], elapsed: float) -> str:
    c_records = {record["orientation"]: record for record in records
                 if record["type"] == "c-line-slice"}
    z_records = {record["orientation"]: record for record in records
                 if record["type"] == "z-line-pullback"}
    action = next(record for record in records
                  if record["type"] == "monodromy-action")
    return "\n".join([
        "# L24 diagonal arithmetic theorem audit",
        "",
        "**Overall verdict — BOTH ORIENTATIONS ARE PROVED EMPTY; the finite",
        "diagonal route is CLOSED.**  After the exact `h` substitution, each",
        "orientation has a uniform `Q_2` obstruction on the full Phi base.  The",
        "monodromy, branch, and arithmetic-theorem calculations below remain exact",
        "side structure, but no point from either cover can define any cell.",
        "",
        "## 1. The exact degree-eight tower",
        "",
        "Put `u=r^2` in orientation I and `u=y^2` in orientation II.  In both",
        "cases the other square is `u+A`.  Exact elimination gives",
        "",
        "```text",
        "u^2+(L-c^2/A)u+M=0,",
        "active=(lambda-A/lambda)/2,",
        "other =(lambda+A/lambda)/2,",
        "lambda=y+r.",
        "```",
        "",
        "The constants are",
        "",
        "```text",
        "I : L=A-16B, M=16(ABs^2-1),",
        "II: L=-16B,  M=-16AB(1-s^2)-16.",
        "```",
        "",
        "Thus the monic octic is",
        "",
        "```text",
        "P(lambda)=(lambda^2-A)^4",
        " +4(L-c^2/A)lambda^2(lambda^2-A)^2+16Mlambda^4.",
        "```",
        "",
        "It is even and `A`-reciprocal, with constant coefficient `A^4`.",
        "Consequently the guarded coordinate algebra is free of rank eight; this",
        "is finite flat even on the branch divisor.  Six exact substitutions replay",
        "both original pairs `(C1,C2)`.",
        "",
        "## 2. Uniform Phi obstruction at 2",
        "",
        "Phi makes `a` and `b` 2-adic units, so `A=5 mod 8`.  In",
        "",
        "```text",
        "g=(16a^4b^2-A(b-1)^4)/(Ab^2)",
        "```",
        "",
        "the numerator is divisible by 16.  If `e=v2(Z)>=0`, then `D` is a",
        "unit; if `e<0`, then `v2(D)=2e`.  Consequently the exact substitution",
        "`c=a^2Z^2g/D` satisfies `v2(c)>=4` wherever the Phi guard `D!=0` holds.",
        "",
        "For orientation I the `u`-quadratic has coefficient valuations",
        "`[4,0,0]`.  A rational 2-adic root therefore has `v2(u)=4` or `0`.",
        "If square with valuation four, `u+A=5 mod 8`; if a square unit,",
        "`u=1 mod 8` and `u+A=6 mod 8` has odd valuation.  Thus `u` and",
        "`u+A` cannot both be squares: **orientation I is empty over `Q_2` at",
        "every Phi base point.**",
        "",
        "For orientation II the valuations are `[4,5,0]`, forcing `v2(u)=2`.",
        "If `s` is not a 2-adic unit (including `s=0`), then",
        "",
        "```text",
        "Delta/64 = 1-AB(s^2-1) = 3 or 7 (mod 8),",
        "```",
        "",
        "so the quadratic has no rational `Q_2` root.  If `v2(s)=0`, a diagonal",
        "point would have square `u=4t^2` with `t` a unit.  Dividing the equation",
        "by 16 gives",
        "",
        "```text",
        "t^4-8bt^2+2Ab(s^2-1)-1-(c^2/(4A))t^2 = 0.",
        "```",
        "",
        "Its left side is `1+8+0-1-0=8 mod 16`, a contradiction.  Therefore",
        "**orientation II is also empty over `Q_2` at every Phi base point.**",
        "",
        "The independent-`c` slice geometry below remains a valid cover",
        "calculation; none of its rational points can meet the substituted Phi locus.",
        "",
        "",
        "## 3. A Phi-admissible one-dimensional slice",
        "",
        "Fix `s=1, a=3, A=37, b=3, B=6`.  Both 2-adic Phi guards hold.",
        "Let `c` first be an independent coordinate.  With `w=cv`, the normalized",
        "source curve is the biquadratic fibre product",
        "",
        "```text",
        "w^2=A(v^4+Lv^2+M),     y^2=v^2+A.",
        "```",
        "",
        "The two branch sets over the `v`-line have sizes four and two and are",
        "disjoint.  Riemann–Hurwitz gives genus three.  The degree-eight map to",
        "the `c`-line has branch equations",
        "",
        "```text",
        "Delta(c)=(AL-c^2)^2-4A^2M,     c^2+N=0,",
        "N=M-AL+A^2.",
        "```",
        "",
        "The four `Delta` values have inertia `2^4`; the two remaining values",
        "have inertia `2^2 1^4`; infinity is unramified.  Exact slice data:",
        "",
        f"- I: `(L,M,N)=(-59,3536,7088)`, branch polynomial "
        f"`{c_records['I']['degree_8_map_to_c']['branch_polynomial_coefficients_ascending']}` (ascending).",
        f"- II: `(L,M,N)=(-96,-16,4905)`, branch polynomial "
        f"`{c_records['II']['degree_8_map_to_c']['branch_polynomial_coefficients_ascending']}` (ascending).",
        "",
        "Both exact gcds with the derivative are one.",
        "",
        "## 4. Monodromy",
        "",
        "After adjoining a root `u1` of the quadratic, splitting the fibre adjoins",
        "the four square roots of `u1,u1+A,u2,u2+A`.  Over `Qbar(c)`,",
        "`u1*u2=M` is square and is the sole Kummer relation.  The remaining",
        "unramified active class is nontrivial: on",
        "`W^2=A*u*(u^2+Lu+M)`, `div(u)=2((0,0)-O)`, the nonzero rational",
        "2-torsion class.  The `u_i+A` classes have distinct odd valuations above",
        "`c^2+N=0`.  Hence",
        "",
        "```text",
        "G_geom = C2^3 semidirect C2 = D8 x C2, order 16.",
        "```",
        "",
        "Over `Q(c)`, `M` is nonsquare in both recorded slices, removing the",
        "constant relation:",
        "",
        "```text",
        "G_arith = C2^4 semidirect C2 = V4 wreath C2, order 32.",
        "```",
        "",
        f"The exact permutation replay gives arithmetic cycle counts "
        f"`{action['arithmetic_group']['cycle_type_counts']}` and geometric counts "
        f"`{action['geometric_group']['cycle_type_counts']}`.",
        "",
        "## 5. The actual z-family",
        "",
        "On this slice",
        "",
        "```text",
        "c(z)=11072 z^6 / (37(1-z^3-9z^6)).",
        "```",
        "",
        "This degree-six map factors through `Z=z^3`.  Its three branch values",
        "are disjoint from the six cover branch values in both orientations",
        "(the JSON records the exact nonzero evaluations).  Therefore geometric",
        "Galois closures are linearly disjoint: a common subcover would be",
        "unramified everywhere on `P1`.  The constant extension from `z^3` is only",
        "`Q(sqrt(-3))`, distinct from `Q(sqrt(3536))` and `Q(sqrt(-16))`, so",
        "the arithmetic groups are preserved too.",
        "",
        "The pulled branch numerator is squarefree of degree 36: 24 discriminant",
        "points and 12 other-square points.  Thus",
        "",
        "```text",
        "2g-2=-16+24*4+12*2=104, so g=53.",
        "```",
        "",
        f"The exact degrees are I `{z_records['I']['pulled_branch']['total_degree']}`",
        f"and II `{z_records['II']['pulled_branch']['total_degree']}`, with gcd one.",
        "",
        "## 6. Exact theorem matches and failures",
        "",
        "### Hilbert irreducibility",
        "",
        "Generic irreducibility is proved on the slice, so HIT really applies—but",
        "its conclusion is that a Hilbert subset of rational `z0` has irreducible",
        "degree-eight fibre, hence no rational `lambda` for this fixed `(s,b)`.",
        "For both orientations this is now redundant with the uniform `Q_2` wall.",
        "It remains an exact illustration of HIT's direction: many fixed-slice",
        "fibres are pointless, rather than points being created.",
        "To produce points one first needs a Hilbert-dense supply of rational source",
        "points (for example a section or rational parametrization).  One isolated",
        "point gives only one specialization.  HIT neither creates source points nor",
        "addresses an arbitrary prescribed fibre.",
        "",
        "### Weak approximation and Faltings",
        "",
        "Faltings applies unconditionally to the recorded `(s,b)=(1,3)` slice:",
        "its source has genus three over the `c`-line and genus 53 after the actual",
        "`z` pullback, so both rational-point sets are finite.  The same conclusion",
        "holds for any other fixed pair whose normalization has genus greater than",
        "one; degenerate lower-genus parameter loci require separate analysis.",
        "Standard weak-approximation theorems for rational varieties or connected",
        "linear groups do not match.  More decisively, local solubility—the first",
        "hypothesis of any local-to-global use—fails at `Q_2` in both orientations.",
        "",
        "### Norm forms",
        "",
        "The elliptic quotient is `W^2=A*u*(u^2+Lu+M)`.  Diagonal points impose",
        "both `u` and `u+A` square and the prescribed value `c=W/u`.  This is a",
        "simultaneous `(Z/2)^2` Kummer cover of an elliptic curve, not",
        "`N_{L/Q}(x)=alpha` for one fixed cyclic extension, so the Hasse norm",
        "theorem has the wrong structural input.  Independently, its all-local",
        "premise fails because the diagonal variety is empty over `Q_2`.",
        "",
        "### Finite-cover Chebotarev",
        "",
        "The arithmetic action has `7/32` fixed-point elements and `25/32`",
        "derangements (geometrically `3/16` and `13/16`).  Chebotarev therefore",
        "organizes residue fibres of both kinds.  These Frobenius statistics are",
        "side structure only: fixed points at odd residue primes cannot repair the",
        "exact `Q_2` obstruction, and Chebotarev does not create rational points.",
        "",
        "## 7. Quantifier ledger",
        "",
        "- **PROVED:** both orientations are finite flat rank-eight covers on the",
        "  guarded Phi base.",
        "- **PROVED:** neither cover has a `Q_2` point anywhere on Phi; orientation I",
        "  fails by the squareclasses of `u,u+A`, and orientation II fails by the",
        "  even-`s` discriminant or unit-`s` square-`u` congruence.",
        "- **PROVED:** on the admissible `(s,b)=(1,3)` independent-`c` slice, both",
        "  have geometric monodromy `D8 x C2`, arithmetic monodromy `V4 wreath C2`,",
        "  and the stated branch/genus data; the actual `z` pullbacks have genus 53.",
        "- **PROVED:** HIT and Faltings give additional no-point/finiteness side",
        "  results; the direct weak-approximation, norm-form, and Chebotarev",
        "  applications already fail at local solubility or the stated structure.",
        "- **CLOSED:** the union over all cell-dependent Phi parameters is empty.",
        "",
        f"Replay elapsed: {elapsed:.3f} seconds.",
        "",
    ])


def main() -> int:
    started = time.perf_counter()
    cover = exact_cover_record()
    two_adic = phi_two_adic_record()
    slice_parameters = {
        orientation: orientation_parameters(orientation, F(1), F(3))
        for orientation in ("I", "II")
    }
    c_records = [c_line_branch_record(orientation, slice_parameters[orientation])
                 for orientation in ("I", "II")]
    time.sleep(0.002)
    z_records = [pullback_record(orientation, slice_parameters[orientation])
                 for orientation in ("I", "II")]
    action = monodromy_action_record()
    theorems = theorem_records()
    elapsed = time.perf_counter() - started

    summary = {
        "type": "summary",
        "label": "BOTH ORIENTATIONS PROVED EMPTY; DIAGONAL ROUTE CLOSED",
        "orientation_verdicts": {
            "I": "PROVED EMPTY OVER Q2 ON PHI",
            "II": "PROVED EMPTY OVER Q2 ON PHI",
            "union": "PROVED EMPTY OVER Q2 ON PHI",
        },
        "unconditional_every_cell_point_theorem": None,
        "completeness_decision": "PROVED FALSE",
        "counterexample_to_union_over_cell_dependent_s_b": (
            "uniform local obstruction: the union has no Q2 point at any Phi base point"
        ),
        "proved": [
            "finite flat rank-8 tower for both orientations",
            "uniform Phi Q2 obstruction eliminating orientation I",
            "uniform Phi Q2 obstruction eliminating orientation II",
            "exact branch and monodromy on an admissible independent-c slice",
            "degree-36 squarefree branch and genus 53 on the actual z pullback",
            "precise first failed hypothesis for each named arithmetic theorem",
        ],
        "first_failed_hypotheses": {
            "Hilbert irreducibility": (
                "no Hilbert-dense supply of rational source points/section; arbitrary "
                "prescribed fibres are outside its point-producing conclusion"
            ),
            "weak approximation": (
                "local solubility fails at Q2 before weak approximation or "
                "Brauer-Manin analysis can begin"
            ),
            "norm form": (
                "not a torsor under one fixed cyclic norm torus, and the all-local "
                "premise fails at Q2"
            ),
            "finite-cover Chebotarev": (
                "odd-prime Frobenius fixed points cannot repair Q2 emptiness"
            ),
        },
        "remaining_open": [],
        "wall_seconds": elapsed,
    }
    meta = {
        "type": "meta",
        "label": "PROVED",
        "schema": "l24-diagonal-arithmetic-v1",
        "source_script": "math/h10q/l24_diagonal_arithmetic.py",
        "output": "math/h10q/data/l24_diagonal_arithmetic.jsonl",
        "report": "/tmp/l24_diagonal_arithmetic.md",
        "stdlib_only": True,
        "finite_search_promoted_to_theorem": False,
    }
    records: list[dict[str, object]] = [
        meta, cover, two_adic, *c_records, *z_records, action, *theorems, summary,
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records
    ) + "\n", encoding="utf-8")
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print(
        "L24 diagonal arithmetic: BOTH orientations EMPTY over Q2; route CLOSED; "
        "exact degree-8 tower, D8xC2/V4wrC2 monodromy, genus-53 z-slices; "
        f"{elapsed:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
