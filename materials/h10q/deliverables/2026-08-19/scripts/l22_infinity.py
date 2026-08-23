#!/usr/bin/env python3
"""Exact higher-order Newton analysis of the constructed octic at a=infinity.

This is a stdlib-only symbolic replay. Sparse Laurent polynomials prove the
identities; truncated series are used only after the formal implicit-function
argument has identified the unique analytic centres. The finite-field row is
an exact representative replay and is never promoted to a uniform claim.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l22_infinity.jsonl"
REPORT = Path("/tmp/l22_infinity.md")

VARIABLES = ("a", "Z", "b", "t", "y", "x")
A_VAR, Z_VAR, B_VAR, T_VAR, Y_VAR, X_VAR = range(len(VARIABLES))
ZERO_EXPONENT = (0,) * len(VARIABLES)


class Poly:
    """Small exact sparse Laurent polynomial over Q."""

    def __init__(self, terms: dict[tuple[int, ...], F | int] | None = None) -> None:
        clean: dict[tuple[int, ...], F] = {}
        for exponent, coefficient in (terms or {}).items():
            assert len(exponent) == len(VARIABLES)
            value = F(coefficient)
            if value:
                clean[exponent] = clean.get(exponent, F(0)) + value
        self.terms = {exponent: coefficient for exponent, coefficient in clean.items()
                      if coefficient}

    @staticmethod
    def constant(value: F | int) -> "Poly":
        value = F(value)
        return Poly({ZERO_EXPONENT: value}) if value else Poly()

    @staticmethod
    def variable(index: int, exponent: int = 1) -> "Poly":
        powers = list(ZERO_EXPONENT)
        powers[index] = exponent
        return Poly({tuple(powers): F(1)})

    @staticmethod
    def coerce(value: "Poly" | F | int) -> "Poly":
        return value if isinstance(value, Poly) else Poly.constant(value)

    def __add__(self, other: "Poly" | F | int) -> "Poly":
        other = Poly.coerce(other)
        terms = dict(self.terms)
        for exponent, coefficient in other.terms.items():
            terms[exponent] = terms.get(exponent, F(0)) + coefficient
        return Poly(terms)

    __radd__ = __add__

    def __neg__(self) -> "Poly":
        return Poly({exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: "Poly" | F | int) -> "Poly":
        return self + (-Poly.coerce(other))

    def __rsub__(self, other: "Poly" | F | int) -> "Poly":
        return Poly.coerce(other) - self

    def __mul__(self, other: "Poly" | F | int) -> "Poly":
        other = Poly.coerce(other)
        terms: dict[tuple[int, ...], F] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in other.terms.items():
                exponent = tuple(left + right for left, right in
                                 zip(left_exponent, right_exponent))
                terms[exponent] = (terms.get(exponent, F(0))
                                   + left_coefficient * right_coefficient)
        return Poly(terms)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Poly":
        assert exponent >= 0
        answer = Poly.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power >>= 1
        return answer

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (Poly, int, F)):
            return False
        return self.terms == Poly.coerce(other).terms

    def substitute(self, replacements: dict[int, "Poly" | F | int]) -> "Poly":
        answer = Poly()
        for exponent, coefficient in self.terms.items():
            term = Poly.constant(coefficient)
            for index, power in enumerate(exponent):
                if not power:
                    continue
                assert power >= 0, "substitution source must be an ordinary polynomial"
                base = Poly.coerce(replacements.get(index, Poly.variable(index)))
                term = term * base**power
            answer = answer + term
        return answer

    def coefficient(self, variable: int, exponent: int) -> "Poly":
        terms: dict[tuple[int, ...], F] = {}
        for powers, coefficient in self.terms.items():
            if powers[variable] == exponent:
                reduced = list(powers)
                reduced[variable] = 0
                key = tuple(reduced)
                terms[key] = terms.get(key, F(0)) + coefficient
        return Poly(terms)

    def degree(self, variable: int) -> int:
        assert self.terms
        return max(exponent[variable] for exponent in self.terms)

    def order(self, variable: int) -> int:
        assert self.terms
        return min(exponent[variable] for exponent in self.terms)

    def leading_coefficient(self, variable: int) -> "Poly":
        return self.coefficient(variable, self.degree(variable))

    def initial_coefficient(self, variable: int) -> "Poly":
        return self.coefficient(variable, self.order(variable))

    def at_zero(self, variable: int) -> "Poly":
        return self.coefficient(variable, 0)

    def has_only_nonnegative_exponents(self, variable: int) -> bool:
        return all(exponent[variable] >= 0 for exponent in self.terms)


class Pacer:
    def __init__(self) -> None:
        self.iterations = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.iterations += 1
        if self.iterations % 200 == 0:
            time.sleep(0.005)
            self.sleeps += 1


PACER = Pacer()


def lower_hull(points: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            x1, y1 = hull[-2]
            x2, y2 = hull[-1]
            x3, y3 = point
            old_slope = F(y2 - y1, x2 - x1)
            new_slope = F(y3 - y2, x3 - x2)
            if old_slope >= new_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def series_mul(left: list[F], right: list[F], degree: int) -> list[F]:
    answer = [F(0)] * (degree + 1)
    for i, left_coefficient in enumerate(left[:degree + 1]):
        for j, right_coefficient in enumerate(right[:degree + 1 - i]):
            answer[i + j] += left_coefficient * right_coefficient
    return answer


def series_add(left: list[F], right: list[F], degree: int) -> list[F]:
    return [
        (left[index] if index < len(left) else F(0))
        + (right[index] if index < len(right) else F(0))
        for index in range(degree + 1)
    ]


def series_pow(value: list[F], exponent: int, degree: int) -> list[F]:
    answer = [F(1)] + [F(0)] * degree
    for _ in range(exponent):
        answer = series_mul(answer, value, degree)
    return answer


def series_inverse(value: list[F], degree: int) -> list[F]:
    assert value[0]
    answer = [F(1) / value[0]] + [F(0)] * degree
    for index in range(1, degree + 1):
        answer[index] = (-sum(value[j] * answer[index - j]
                              for j in range(1, index + 1)) / value[0])
    return answer


def n_left_series(y_series: list[F], degree: int) -> list[F]:
    t_y = [F(0)] + y_series[:degree]
    one_minus_t_y = series_add(
        [F(1)] + [F(0)] * degree,
        [-coefficient for coefficient in t_y],
        degree,
    )
    t2_plus_4 = [F(4), F(0), F(1)] + [F(0)] * max(0, degree - 2)
    return series_add(
        [16 * coefficient for coefficient in
         series_mul(y_series, y_series, degree)],
        [-coefficient for coefficient in
         series_mul(t2_plus_4, series_pow(one_minus_t_y, 4, degree), degree)],
        degree,
    )


def analytic_center(sign: int, degree: int) -> list[F]:
    """Unique r_sign with n_L(t,r_sign)=0 and r_sign(0)=sign/2."""
    assert sign in (-1, 1)
    answer = [F(sign, 2)] + [F(0)] * degree
    derivative_at_origin = F(16 * sign)
    for index in range(1, degree + 1):
        coefficient = n_left_series(answer, degree)[index]
        answer[index] -= coefficient / derivative_at_origin
        assert n_left_series(answer, degree)[index] == 0
    assert all(coefficient == 0 for coefficient in n_left_series(answer, degree))
    return answer


def series_poly(coefficients: list[F], variable: int = T_VAR) -> Poly:
    answer = Poly()
    for exponent, coefficient in enumerate(coefficients):
        answer += coefficient * Poly.variable(variable, exponent)
    return answer


def fraction_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def mod_fraction(value: F, prime: int) -> int:
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, prime - 2, prime) % prime


def univariate_mod(poly: Poly, variable: int, prime: int) -> dict[int, int]:
    answer: dict[int, int] = {}
    for exponent, coefficient in poly.terms.items():
        assert all(power == 0 for index, power in enumerate(exponent) if index != variable)
        degree = exponent[variable]
        answer[degree] = (answer.get(degree, 0) + mod_fraction(coefficient, prime)) % prime
    return {degree: coefficient for degree, coefficient in answer.items() if coefficient}


def roots_mod_quadratic(target_square: F, prime: int) -> list[int]:
    target = mod_fraction(target_square, prime)
    roots: list[int] = []
    for candidate in range(prime):
        PACER.tick()
        if candidate * candidate % prime == target:
            roots.append(candidate)
    return roots


def write_report(records: list[dict[str, object]]) -> None:
    lines = [
        "# L22: higher-order Newton polygons at `a = infinity`",
        "",
        "## Status",
        "",
        "- **PROVED:** the exact coefficient degrees and ordinary lower hull.",
        "- **PROVED:** the repeated left and right residuals and all four second-order clusters.",
        "- **PROVED:** over `Q(Z)((t))`, `Z != 0`, the octic has four irreducible quadratic local factors. Therefore a proper global factor can only have degree `2`, `4`, or `6`.",
        "- **PROVED:** the natural exact place `a=0` sharpens many rational fibers to a possible `4+4` split, but has an explicit infinite exceptional locus where it leaves degrees `2,4,6`.",
        "- **OPEN:** irreducibility for every fixed rational nonzero `Z`. Nothing here promotes a finite replay to that theorem.",
        "",
        "Throughout,",
        "",
        "```text",
        "A=1+4a^2,  D=1-Z-a^2 Z^2,",
        "N_g=16a^4 b^2-A(b-1)^4,",
        "H=a^8 Z^4 N_g^2+4A^3D^2b^4-2A^4(a-1)^2D^2b^5.",
        "```",
        "",
        "Multiplication by the nonzero scalar relating `H` and the normalized octic does not affect factorization in `b`.",
        "",
        "## L22-I1. Coefficient degrees and ordinary polygon (PROVED)",
        "",
        "Writing `N_g=sum n_i b^i`, its exact coefficient vector is",
        "",
        "```text",
        "(-A, 4A, 16a^4-6A, 4A, -A).",
        "```",
        "",
        "Squaring and adding the two central terms gives:",
        "",
        "| b-degree i | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        "| deg_a H_i | 12 | 12 | 14 | 14 | 16 | 13 | 14 | 12 | 12 |",
        "| lc_a(H_i)/Z^4 | 16 | -128 | -128 | 512 | 256 | 1024 | -128 | -128 | 16 |",
        "",
        "At `i=5`, the degree-14 term `512 Z^4 a^14` from `a^8Z^4N_g^2` cancels exactly with `-512 Z^4 a^14` from the perturbation; the next term is `1024 Z^4 a^13`. For `t=1/a`, the coefficient valuations are the negatives of these degrees, and the lower hull is",
        "",
        "```text",
        "(0,-12) -> (4,-16) -> (8,-12).",
        "```",
        "",
        "Its two length-four residuals are not squarefree, so first order is inconclusive.",
        "",
        "## L22-I2. Left chart and second-order edges (PROVED)",
        "",
        "Put `b=t*y` and `F_L=t^12 H`. Direct Laurent substitution gives",
        "",
        "```text",
        "n_L = 16y^2-(t^2+4)(ty-1)^4,",
        "d   = (1-Z)t^2-Z^2,",
        "F_L = Z^4 n_L^2",
        "      -2t^3(t^2+4)^4(1-t)^2 d^2 y^5",
        "      +4t^6(t^2+4)^3 d^2 y^4.",
        "```",
        "",
        "At `t=0`, `F_L(0,y)=16 Z^4(1-4y^2)^2`.",
        "",
        "For `sigma=+1,-1`, `partial_y n_L(0,sigma/2)=16 sigma`, so formal implicit-function theory gives a unique `r_sigma(t) in Q[[t]]` with `n_L(t,r_sigma)=0` and `r_sigma(0)=sigma/2`:",
        "",
        "```text",
        "r_+(t)= 1/2 - t/2 + 11t^2/16 - t^3 + O(t^4),",
        "r_-(t)=-1/2 - t/2 - 11t^2/16 - t^3 + O(t^4).",
        "```",
        "",
        "At `y=r_sigma+x`, the coefficients of `x^0,x^1,x^2` have valuations `(3,3,0)` and initial coefficients `Z^4*(-16 sigma,-160,256)`. Thus the new edge is `(0,3)->(2,0)`, slope `-3/2`. With `x=t^(3/2)v`,",
        "",
        "```text",
        "256 Z^4 v^2-16 sigma Z^4=0,  so v^2=sigma/16.",
        "```",
        "",
        "The small clusters are `b=t r_sigma(t)+t^(5/2)v+O(t^3)`. For `sigma=+1`, `v=+/-1/4`; for `sigma=-1`, `v=+/-i/4`.",
        "",
        "## L22-I3. Right chart and second-order edges (PROVED)",
        "",
        "Put `b=y/t` and `F_R=t^20 H`. Exact substitution gives",
        "",
        "```text",
        "n_R = 16y^2-(t^2+4)(y-t)^4,",
        "F_R = Z^4 n_R^2",
        "      -2t(t^2+4)^4(1-t)^2 d^2 y^5",
        "      +4t^6(t^2+4)^3 d^2 y^4.",
        "```",
        "",
        "Hence `F_R(0,y)=16Z^4 y^4(4-y^2)^2`. Also `n_R(t,1/y)=y^(-4)n_L(t,y)`, so the analytic centres are `R_sigma=1/r_sigma`:",
        "",
        "```text",
        "R_+(t)= 2+2t-3t^2/4+t^3/2+O(t^4),",
        "R_-(t)=-2+2t+3t^2/4+t^3/2+O(t^4).",
        "```",
        "",
        "At `y=R_sigma+x`, the `x^0,x^1,x^2` coefficients have valuations `(1,1,0)` and initials `Z^4*(-16384 sigma,-40960,4096)`. The new edge `(0,1)->(2,0)` has slope `-1/2`; with `x=t^(1/2)v`,",
        "",
        "```text",
        "4096 Z^4 v^2-16384 sigma Z^4=0,  so v^2=4 sigma.",
        "```",
        "",
        "The large clusters are `b=t^(-1)R_sigma(t)+t^(-1/2)v+O(1)`. For `sigma=+1`, `v=+/-2`; for `sigma=-1`, `v=+/-2i`.",
        "",
        "## L22-I4. Exact local factor degrees (PROVED)",
        "",
        "Let `K=Q(Z)((t))`, `Z != 0`. Coprime first-residual factors separate the four sign/size clusters by Hensel-Weierstrass preparation. Each cluster has length two and reduced slope denominator two. In the standard denominator-two coordinate, the edge residuals (up to nonzero scalar) are `16U-sigma` for `U=x^2/t^3` on the small side and `U-4sigma` for `U=x^2/t` on the large side. Equivalently, its Weierstrass quadratic has",
        "",
        "```text",
        "small: q_sigma=-(sigma/16)t^3+O(t^4),  disc=(sigma/4)t^3+O(t^4),",
        "large: q_sigma=-4sigma*t+O(t^2),       disc=16sigma*t+O(t^2).",
        "```",
        "",
        "Every discriminant has odd `t`-valuation, so each cluster is one irreducible quadratic over `K`. The exact local degree multiset is `{2,2,2,2}`. The `+` quadratics split over `K(sqrt(t))`; the `-` quadratics split over `K(sqrt(-t))`. Small and large clusters of a fixed sign therefore have the same, not incompatible, ramified square class. Because every initial coefficient used above is a nonzero multiple of `Z^4`, the same conclusion specializes to every fixed `Z in Q^*`.",
        "",
        "Every divisor is a product of whole local quadratics. Its degree is therefore one of `0,2,4,6,8`; the possible proper nonconstant degrees are `2,4,6`. This is a local permission statement, not an assertion that a proper factor exists.",
        "",
        "## L22-I5. Involution and natural finite place (PROVED obstruction)",
        "",
        "Exact calculation gives",
        "",
        "```text",
        "H(-a,b)-H(a,b)=-8a A^4D^2 b^5,",
        "H(a,b)-b^8H(a,1/b)=2A^4(a-1)^2D^2(b^3-b^5).",
        "```",
        "",
        "Also `gcd_b(H(a,b),H(-a,b))=1`: a common nonconstant divisor would divide their monomial difference, hence `b`, while the constant coefficient `a^8A^2Z^4` is nonzero. More generally, suppose `H(-a,b)=c b^8 H(a,k/b)` with `c,k in Q(a)^*`. The equal nonzero endpoints force `c=1` and `k^8=1`, hence `k=+/-1`; the nonzero `b^1` coefficient forces `k=1`, after which the displayed `b^3/b^5` defect makes the identity impossible. Thus `a -> -a` maps a factorization to one of a different polynomial and imposes no pairing of the local clusters.",
        "",
        "At `a=0`, for fixed rational `Z != 0,1`, the coefficient valuations are `(8,8,8,8,0,0,8,8,8)` and the hull is `(0,8)->(4,0)->(5,0)->(8,8)`. The clusters have degrees `4,1,3`; the small residual is",
        "",
        "```text",
        "R_Z(y)=Z^4+4(1-Z)^2 y^4.",
        "```",
        "",
        "Coefficient comparison proves, for rational `Z != 0,1`,",
        "",
        "```text",
        "R_Z reducible over Q  <=>  1-Z is a square or Z-1 is a square in Q.",
        "```",
        "",
        "Indeed, after division by `4(1-Z)^2`, this is `y^4+c^2`, `c=Z^2/(2(1-Z))`; it is reducible exactly when `2c` or `-2c` is a square. On either locus `Z=1+/-q^2`,",
        "",
        "```text",
        "Z^4+4q^4y^4",
        "=(2q^2y^2+2Zqy+Z^2)(2q^2y^2-2Zqy+Z^2).",
        "```",
        "",
        "For nonexceptional `Z`, combining zero with infinity leaves only a possible `4+4`. On the infinite exceptional locus the zero-place local degrees are `{2,2,1,3}`, so `2,4,6` still survive. At `Z=1`, the small residual is `1+4y^4=4(y^2+y+1/2)(y^2-y+1/2)`, the same unresolved pattern.",
        "",
        "Thus these exact places do not prove irreducibility for every rational nonzero `Z`. The sharp universal infinity obstruction is: **odd degrees are impossible, but degrees 2, 4, and 6 all survive the local degree test**. Away from the explicit zero-place exceptional locus, only quartic-plus-quartic survives the two-place test.",
        "",
        "## Representative finite-field/Laurent replay",
        "",
        "**PROVED for the stated replay only.** The script reduces the exact identities modulo `17`, takes `Z=27`, sets `t=s^2`, and substitutes one root of each second residual. It checks cancellation through the residual order in all four clusters. This is not evidence for a uniform irreducibility theorem.",
        "",
        "## Verdict",
        "",
        "**OPEN.** Higher-order analysis resolves the repeated residuals into four ramified irreducible quadratics rather than one degree-eight orbit. The sign symmetry is absent, and `a=0` has the displayed square exceptional locus. A new global argument is required to exclude the surviving even-degree factorizations uniformly.",
        "",
        f"Machine records: `{OUT}` ({len(records)} JSONL rows).",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()

    a = Poly.variable(A_VAR)
    Z = Poly.variable(Z_VAR)
    b = Poly.variable(B_VAR)
    t = Poly.variable(T_VAR)
    y = Poly.variable(Y_VAR)
    x = Poly.variable(X_VAR)

    A = 1 + 4 * a**2
    D = 1 - Z - a**2 * Z**2
    N_g = 16 * a**4 * b**2 - A * (b - 1)**4
    reciprocal_core = a**8 * Z**4 * N_g**2 + 4 * A**3 * D**2 * b**4
    perturbation = 2 * A**4 * (a - 1)**2 * D**2
    H = reciprocal_core - perturbation * b**5
    ng_coefficients = [N_g.coefficient(B_VAR, degree) for degree in range(5)]
    assert ng_coefficients == [-A, 4 * A, 16 * a**4 - 6 * A, 4 * A, -A]
    square_b5 = (a**8 * Z**4 * N_g**2).coefficient(B_VAR, 5)
    perturbation_b5 = (-perturbation * b**5).coefficient(B_VAR, 5)
    assert square_b5.degree(A_VAR) == perturbation_b5.degree(A_VAR) == 14
    assert square_b5.leading_coefficient(A_VAR) == 512 * Z**4
    assert perturbation_b5.leading_coefficient(A_VAR) == -512 * Z**4

    coefficients = [H.coefficient(B_VAR, degree) for degree in range(9)]
    assert H.degree(B_VAR) == 8
    assert sum((coefficient * b**degree for degree, coefficient in
                enumerate(coefficients)), Poly()) == H

    a_degrees = [coefficient.degree(A_VAR) for coefficient in coefficients]
    expected_degrees = [12, 12, 14, 14, 16, 13, 14, 12, 12]
    assert a_degrees == expected_degrees
    expected_leads = [16, -128, -128, 512, 256, 1024, -128, -128, 16]
    for coefficient, expected in zip(coefficients, expected_leads):
        assert coefficient.leading_coefficient(A_VAR) == expected * Z**4
    infinity_points = [(degree, -a_degree)
                       for degree, a_degree in enumerate(a_degrees)]
    infinity_hull = lower_hull(infinity_points)
    assert infinity_hull == [(0, -12), (4, -16), (8, -12)]

    d = (1 - Z) * t**2 - Z**2
    n_left = 16 * y**2 - (t**2 + 4) * (t * y - 1)**4
    F_left = (Z**4 * n_left**2
              - 2 * t**3 * (t**2 + 4)**4 * (1 - t)**2 * d**2 * y**5
              + 4 * t**6 * (t**2 + 4)**3 * d**2 * y**4)
    from_H_left = H.substitute({A_VAR: Poly.variable(T_VAR, -1),
                                B_VAR: t * y}) * t**12
    assert from_H_left == F_left
    assert F_left.has_only_nonnegative_exponents(T_VAR)
    left_first_residual = 16 * Z**4 * (1 - 4 * y**2)**2
    assert F_left.at_zero(T_VAR) == left_first_residual

    n_right = 16 * y**2 - (t**2 + 4) * (y - t)**4
    F_right = (Z**4 * n_right**2
               - 2 * t * (t**2 + 4)**4 * (1 - t)**2 * d**2 * y**5
               + 4 * t**6 * (t**2 + 4)**3 * d**2 * y**4)
    from_H_right = H.substitute({A_VAR: Poly.variable(T_VAR, -1),
                                 B_VAR: Poly.variable(T_VAR, -1) * y}) * t**20
    assert from_H_right == F_right
    assert F_right.has_only_nonnegative_exponents(T_VAR)
    right_first_residual = 16 * Z**4 * y**4 * (4 - y**2)**2
    assert F_right.at_zero(T_VAR) == right_first_residual
    reciprocal_n_identity = n_right.substitute({Y_VAR: Poly.variable(Y_VAR, -1)})
    assert reciprocal_n_identity * y**4 == n_left

    series_degree = 8
    centers: dict[int, list[F]] = {}
    reciprocal_centers: dict[int, list[F]] = {}
    expected_centers = {
        1: [F(1, 2), F(-1, 2), F(11, 16), F(-1)],
        -1: [F(-1, 2), F(-1, 2), F(-11, 16), F(-1)],
    }
    expected_reciprocals = {
        1: [F(2), F(2), F(-3, 4), F(1, 2)],
        -1: [F(-2), F(2), F(3, 4), F(1, 2)],
    }
    centered_data: dict[str, dict[str, object]] = {}
    for sign in (1, -1):
        center = analytic_center(sign, series_degree)
        reciprocal_center = series_inverse(center, series_degree)
        centers[sign] = center
        reciprocal_centers[sign] = reciprocal_center
        assert center[:4] == expected_centers[sign]
        assert reciprocal_center[:4] == expected_reciprocals[sign]

        center_poly = series_poly(center)
        reciprocal_poly = series_poly(reciprocal_center)
        left_centered = F_left.substitute({Y_VAR: center_poly + x})
        right_centered = F_right.substitute({Y_VAR: reciprocal_poly + x})
        left_coefficients = [left_centered.coefficient(X_VAR, degree)
                             for degree in range(3)]
        right_coefficients = [right_centered.coefficient(X_VAR, degree)
                              for degree in range(3)]
        left_orders = [coefficient.order(T_VAR) for coefficient in left_coefficients]
        right_orders = [coefficient.order(T_VAR) for coefficient in right_coefficients]
        assert left_orders == [3, 3, 0]
        assert right_orders == [1, 1, 0]
        left_initials = [coefficient.initial_coefficient(T_VAR)
                         for coefficient in left_coefficients]
        right_initials = [coefficient.initial_coefficient(T_VAR)
                          for coefficient in right_coefficients]
        assert left_initials == [-16 * sign * Z**4, -160 * Z**4, 256 * Z**4]
        assert right_initials == [-16384 * sign * Z**4,
                                  -40960 * Z**4, 4096 * Z**4]
        centered_data[str(sign)] = {
            "r_initial": [fraction_text(value) for value in center[:4]],
            "R_initial": [fraction_text(value) for value in reciprocal_center[:4]],
            "left_orders_x0_x1_x2": left_orders,
            "left_initials_over_Z4": [-16 * sign, -160, 256],
            "left_puiseux_residual": f"v^2-({fraction_text(F(sign, 16))})",
            "left_discriminant_lead": f"({fraction_text(F(sign, 4))})*t^3",
            "right_orders_x0_x1_x2": right_orders,
            "right_initials_over_Z4": [-16384 * sign, -40960, 4096],
            "right_puiseux_residual": f"v^2-({4 * sign})",
            "right_discriminant_lead": f"{16 * sign}*t",
        }

    H_negative_a = H.substitute({A_VAR: -a})
    assert H_negative_a - H == -8 * a * A**4 * D**2 * b**5
    reciprocal_H = H.substitute({B_VAR: Poly.variable(B_VAR, -1)}) * b**8
    assert H - reciprocal_H == perturbation * (b**3 - b**5)
    assert H.coefficient(B_VAR, 0) == a**8 * A**2 * Z**4
    assert H.coefficient(B_VAR, 1) != Poly()

    zero_orders = [coefficient.order(A_VAR) for coefficient in coefficients]
    assert zero_orders == [8, 8, 8, 8, 0, 0, 8, 8, 8]
    zero_hull = lower_hull(list(enumerate(zero_orders)))
    assert zero_hull == [(0, 8), (4, 0), (5, 0), (8, 8)]
    zero_small = (H.substitute({B_VAR: a**2 * y})
                  * Poly.variable(A_VAR, -8)).at_zero(A_VAR)
    expected_zero_small = Z**4 + 4 * (1 - Z)**2 * y**4
    assert zero_small == expected_zero_small
    zero_large = (H.substitute({A_VAR: t**3,
                                B_VAR: Poly.variable(T_VAR, -8) * y})
                  * t**40).at_zero(T_VAR)
    assert zero_large == y**5 * (Z**4 * y**3 - 2 * (1 - Z)**2)

    q = x
    for specialized_Z in (1 - q**2, 1 + q**2):
        residual = expected_zero_small.substitute({Z_VAR: specialized_Z})
        plus = (2 * q**2 * y**2 + 2 * specialized_Z * q * y
                + specialized_Z**2)
        minus = (2 * q**2 * y**2 - 2 * specialized_Z * q * y
                 + specialized_Z**2)
        assert residual == plus * minus

    H_Z1 = H.substitute({Z_VAR: 1})
    Z1_coefficients = [H_Z1.coefficient(B_VAR, degree) for degree in range(9)]
    Z1_orders = [coefficient.order(A_VAR) for coefficient in Z1_coefficients]
    assert Z1_orders == [8, 8, 8, 8, 4, 4, 8, 8, 8]
    assert lower_hull(list(enumerate(Z1_orders))) == [
        (0, 8), (4, 4), (5, 4), (8, 8)
    ]
    Z1_small = (H_Z1.substitute({B_VAR: a * y})
                * Poly.variable(A_VAR, -8)).at_zero(A_VAR)
    assert Z1_small == 1 + 4 * y**4
    assert Z1_small == 4 * (y**2 + y + F(1, 2)) * (y**2 - y + F(1, 2))

    prime = 17
    replay_Z = 27
    s = x
    replay_rows: list[dict[str, object]] = []
    residual_targets = {
        ("small", 1): F(1, 16),
        ("small", -1): F(-1, 16),
        ("large", 1): F(4),
        ("large", -1): F(-4),
    }
    expected_roots = {
        ("small", 1): [4, 13],
        ("small", -1): [1, 16],
        ("large", 1): [2, 15],
        ("large", -1): [8, 9],
    }
    for chart in ("small", "large"):
        for sign in (1, -1):
            roots = roots_mod_quadratic(residual_targets[(chart, sign)], prime)
            assert roots == expected_roots[(chart, sign)]
            chosen = roots[0]
            if chart == "small":
                centre_in_s = series_poly(centers[sign], X_VAR).substitute(
                    {X_VAR: s**2})
                replay = F_left.substitute({T_VAR: s**2,
                                            Y_VAR: centre_in_s + chosen * s**3,
                                            Z_VAR: replay_Z})
                residual_order = 6
                normalization = "s^24 H(s^-2, s^2*r_sigma(s^2)+v*s^5)"
            else:
                centre_in_s = series_poly(reciprocal_centers[sign], X_VAR).substitute(
                    {X_VAR: s**2})
                replay = F_right.substitute({T_VAR: s**2,
                                             Y_VAR: centre_in_s + chosen * s,
                                             Z_VAR: replay_Z})
                residual_order = 2
                normalization = "s^40 H(s^-2, s^-2*R_sigma(s^2)+v*s^-1)"
            reduced = univariate_mod(replay, X_VAR, prime)
            assert all(reduced.get(degree, 0) == 0
                       for degree in range(residual_order + 1))
            next_order = min(reduced) if reduced else None
            assert next_order is None or next_order > residual_order
            replay_rows.append({
                "chart": chart,
                "sign": sign,
                "prime": prime,
                "Z_integer": replay_Z,
                "Z_mod_prime": replay_Z % prime,
                "residual_target_square": fraction_text(
                    residual_targets[(chart, sign)]),
                "all_residual_roots_mod_prime": roots,
                "chosen_v": chosen,
                "normalization": normalization,
                "cancelled_through_s_order": residual_order,
                "next_nonzero_s_order_in_truncated_replay": next_order,
            })

    records: list[dict[str, object]] = [
        {
            "type": "coefficient-degrees-and-hull",
            "label": "PROVED",
            "proof_mode": "exact sparse symbolic identity",
            "coefficient_a_degrees_b0_to_b8": a_degrees,
            "leading_coefficients_over_Z4": expected_leads,
            "degree_14_cancellation_at_b5": {
                "Ng_square_contribution": 512,
                "perturbation_contribution": -512,
                "next_a13_coefficient_over_Z4": 1024,
            },
            "t_valuation_points": infinity_points,
            "ordinary_lower_hull": infinity_hull,
        },
        {
            "type": "first-residuals",
            "label": "PROVED",
            "proof_mode": "exact Laurent substitution",
            "left_chart": {
                "substitution": "a=t^-1,b=t*y,F_L=t^12*H",
                "residual": "16*Z^4*(1-4*y^2)^2",
                "repeated_roots": ["1/2", "-1/2"],
            },
            "right_chart": {
                "substitution": "a=t^-1,b=y/t,F_R=t^20*H",
                "residual": "16*Z^4*y^4*(4-y^2)^2",
                "repeated_nonzero_roots": ["2", "-2"],
            },
            "reciprocal_center_identity": "n_R(t,1/y)=y^-4*n_L(t,y)",
        },
        {
            "type": "higher-order-clusters",
            "label": "PROVED",
            "proof_mode": "formal implicit functions and exact higher Newton edges",
            "centers_and_edges": centered_data,
            "small_edge": {"vertices": [[0, 3], [2, 0]], "slope": "-3/2"},
            "large_edge": {"vertices": [[0, 1], [2, 0]], "slope": "-1/2"},
            "local_irreducible_factor_degrees": [2, 2, 2, 2],
            "possible_proper_global_factor_degrees": [2, 4, 6],
            "all_local_divisor_degrees_including_unit_and_full_polynomial": [0, 2, 4, 6, 8],
            "warning": "locally permitted does not mean globally realized",
        },
        {
            "type": "sign-and-zero-place-obstruction",
            "label": "PROVED",
            "proof_mode": "exact symbolic identities and local polygons",
            "sign_defect": "H(-a,b)-H(a,b)=-8*a*A^4*D^2*b^5",
            "reciprocal_defect": (
                "H(a,b)-b^8*H(a,1/b)="
                "2*A^4*(a-1)^2*D^2*(b^3-b^5)"
            ),
            "gcd_H_and_sign_conjugate": 1,
            "a0_hull_for_Z_not_1": zero_hull,
            "a0_local_degrees_nonexceptional": [4, 1, 3],
            "a0_small_residual": "Z^4+4*(1-Z)^2*y^4",
            "a0_residual_reducible_iff": (
                "1-Z is a rational square or Z-1 is a rational square"
            ),
            "a0_local_degrees_exceptional": [2, 2, 1, 3],
            "two_place_survivors_nonexceptional": [4],
            "two_place_survivors_exceptional": [2, 4, 6],
            "Z_equals_1_small_residual": "1+4*y^4",
        },
        {
            "type": "finite-field-laurent-replay",
            "label": "PROVED",
            "scope": "exact representative replay only",
            "uniform_implication": False,
            "rows": replay_rows,
            "pacing": {
                "finite_loop_iterations": PACER.iterations,
                "sleeps": PACER.sleeps,
                "sleep_seconds_every_200": "0.005",
            },
        },
        {
            "type": "verdict",
            "label": "OPEN",
            "uniform_irreducibility_for_every_rational_nonzero_Z": False,
            "proved_universal_infinity_obstruction": (
                "odd factor degrees are impossible; degrees 2,4,6 remain locally permitted"
            ),
            "proved_two_place_sharpening": (
                "outside Z=1 and the loci 1-Z square or Z-1 square, only 4+4 remains"
            ),
            "no_overclaim": (
                "no finite replay or local permission is promoted to a global theorem"
            ),
        },
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(record, sort_keys=True) + "\n"
                           for record in records), encoding="utf-8")
    write_report(records)

    elapsed = time.perf_counter() - started
    print(
        "L22 infinity: degrees and both Laurent identities proved; "
        "four local quadratics; global degrees 2/4/6 remain; "
        f"wrote {OUT} and {REPORT}; {elapsed:.3f}s"
    )


if __name__ == "__main__":
    main()
