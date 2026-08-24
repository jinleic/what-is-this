#!/usr/bin/env python3
"""L24 target-place analysis for the two diagonal self-couplings.

This standalone, standard-library-only replay separates the three possible
valuations of the bridge c at an odd target place, verifies the exact residue
character sums, and records the Hensel Jacobians.  Exhaustion is used only for
the finite fields not covered by the displayed Hasse bound.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
from typing import Any, Iterable
import json
import math
import time


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l24_diagonal_local.jsonl"
REPORT = Path("/tmp/l24_diagonal_local.md")
PACE_EVERY = 250
PACE_SECONDS = 0.002
INF = 10**9


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


class FiniteField:
    """Tiny exact F_p[x]/(modulus), with integers encoding coefficient tuples."""

    def __init__(self, p: int, modulus: tuple[int, ...], name: str) -> None:
        assert p >= 3 and all(p % d for d in range(2, math.isqrt(p) + 1))
        assert len(modulus) >= 2 and modulus[-1] % p == 1
        self.p = p
        self.modulus = tuple(value % p for value in modulus)
        self.degree = len(modulus) - 1
        self.q = p**self.degree
        self.name = name
        self.zero = 0
        self.one = 1

    def decode(self, value: int) -> list[int]:
        coefficients = []
        for _ in range(self.degree):
            coefficients.append(value % self.p)
            value //= self.p
        return coefficients

    def encode(self, coefficients: Iterable[int]) -> int:
        answer = 0
        place = 1
        values = list(coefficients)
        assert len(values) == self.degree
        for value in values:
            answer += (value % self.p) * place
            place *= self.p
        return answer

    def const(self, value: int) -> int:
        return value % self.p

    def add(self, left: int, right: int) -> int:
        a = self.decode(left)
        b = self.decode(right)
        return self.encode(x + y for x, y in zip(a, b))

    def neg(self, value: int) -> int:
        return self.encode(-x for x in self.decode(value))

    def sub(self, left: int, right: int) -> int:
        return self.add(left, self.neg(right))

    def mul(self, left: int, right: int) -> int:
        a = self.decode(left)
        b = self.decode(right)
        product = [0] * (2 * self.degree - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                product[i + j] = (product[i + j] + x * y) % self.p
        for degree in range(len(product) - 1, self.degree - 1, -1):
            coefficient = product[degree] % self.p
            if coefficient == 0:
                continue
            product[degree] = 0
            shift = degree - self.degree
            for j in range(self.degree):
                product[shift + j] = (
                    product[shift + j] - coefficient * self.modulus[j]
                ) % self.p
        return self.encode(product[: self.degree])

    def pow(self, value: int, exponent: int) -> int:
        assert exponent >= 0
        answer = self.one
        base = value
        while exponent:
            if exponent & 1:
                answer = self.mul(answer, base)
            base = self.mul(base, base)
            exponent >>= 1
        return answer

    def inv(self, value: int) -> int:
        assert value != 0
        return self.pow(value, self.q - 2)

    def div(self, numerator: int, denominator: int) -> int:
        return self.mul(numerator, self.inv(denominator))

    def scale(self, scalar: int, value: int) -> int:
        return self.mul(self.const(scalar), value)

    def square(self, value: int) -> int:
        return self.mul(value, value)

    def chi(self, value: int) -> int:
        if value == 0:
            return 0
        power = self.pow(value, (self.q - 1) // 2)
        if power == self.one:
            return 1
        assert power == self.const(-1)
        return -1

    def text(self, value: int) -> int | list[int]:
        coefficients = self.decode(value)
        return coefficients[0] if self.degree == 1 else coefficients


SMALL_FIELDS = [
    FiniteField(3, (0, 1), "F_3"),
    FiniteField(5, (0, 1), "F_5"),
    FiniteField(7, (0, 1), "F_7"),
    FiniteField(3, (1, 0, 1), "F_9=F_3[t]/(t^2+1)"),
    FiniteField(11, (0, 1), "F_11"),
    FiniteField(13, (0, 1), "F_13"),
]


def ff_sum(field: FiniteField, *values: int) -> int:
    answer = field.zero
    for value in values:
        answer = field.add(answer, value)
    return answer


def ff_product(field: FiniteField, *values: int) -> int:
    answer = field.one
    for value in values:
        answer = field.mul(answer, value)
    return answer


def matrix_rank(field: FiniteField, matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = field.inv(work[row][column])
        work[row] = [field.mul(value, inverse) for value in work[row]]
        for index in range(len(work)):
            if index == row or work[index][column] == 0:
                continue
            multiple = work[index][column]
            work[index] = [
                field.sub(value, field.mul(multiple, pivot_value))
                for value, pivot_value in zip(work[index], work[row])
            ]
        row += 1
        if row == len(work):
            break
    return row


def base_values(field: FiniteField, a: int) -> tuple[int, int]:
    A = field.add(field.one, field.scale(4, field.square(a)))
    s = field.div(field.sub(a, field.one), field.const(2))
    return A, s


def relation_row(field: FiniteField) -> list[int]:
    # Columns are (a,s,b,y,r), for E=a-1-2s.
    return [field.one, field.const(-2), field.zero, field.zero, field.zero]


def low_jacobian(
    field: FiniteField, orientation: str, a: int, y: int, r: int
) -> list[list[int]]:
    A, _ = base_values(field, a)
    y2, r2 = field.square(y), field.square(r)
    c1 = [
        field.scale(-8 if orientation == "I" else -8, a),
        field.zero,
        field.zero,
        field.scale(2 if orientation == "I" else -2, y),
        field.scale(-2 if orientation == "I" else 2, r),
    ]
    if orientation == "I":
        c2 = [
            field.scale(8, field.mul(a, field.sub(field.mul(r2, y2), field.const(16)))),
            field.zero,
            field.zero,
            field.scale(2, ff_product(field, A, r2, y)),
            field.scale(2, ff_product(field, A, r, y2)),
        ]
    else:
        y4 = field.square(y2)
        c2 = [
            field.scale(8, field.mul(a, field.sub(y4, field.const(16)))),
            field.zero,
            field.zero,
            field.scale(4, ff_product(field, A, y2, y)),
            field.zero,
        ]
    return [relation_row(field), c1, c2]


def unit_jacobian(
    field: FiniteField, orientation: str, a: int, b: int, y: int, r: int
) -> list[list[int]]:
    A, s = base_values(field, a)
    B = field.scale(2, b)
    y2, r2, s2 = field.square(y), field.square(r), field.square(s)
    c1 = [
        field.scale(-8, a),
        field.zero,
        field.zero,
        field.scale(2 if orientation == "I" else -2, y),
        field.scale(-2 if orientation == "I" else 2, r),
    ]
    core = field.sub(field.mul(A, s2), r2)
    if orientation == "I":
        da_inside = ff_sum(
            field,
            field.mul(r2, y2),
            field.neg(field.scale(16, field.mul(B, r2))),
            field.scale(32, ff_product(field, A, B, s2)),
            field.const(-16),
        )
        c2 = [
            field.scale(8, field.mul(a, da_inside)),
            field.scale(32, ff_product(field, A, A, B, s)),
            field.scale(32, field.mul(A, core)),
            field.scale(2, ff_product(field, A, r2, y)),
            field.scale(2, ff_product(field, A, r, field.sub(y2, field.scale(16, B)))),
        ]
    else:
        y4 = field.square(y2)
        da_inside = ff_sum(
            field,
            y4,
            field.neg(field.scale(16, field.mul(B, r2))),
            field.scale(32, ff_product(field, A, B, s2)),
            field.const(-16),
        )
        c2 = [
            field.scale(8, field.mul(a, da_inside)),
            field.scale(32, ff_product(field, A, A, B, s)),
            field.scale(32, field.mul(A, core)),
            field.scale(4, ff_product(field, A, y2, y)),
            field.scale(-32, ff_product(field, A, B, r)),
        ]
    return [relation_row(field), c1, c2]


def critical_jacobian(
    field: FiniteField, orientation: str, a: int, beta: int, y: int, r: int
) -> list[list[int]]:
    A, _ = base_values(field, a)
    y2, r2 = field.square(y), field.square(r)
    beta_inverse = field.inv(beta)
    c2_value = field.pow(field.mul(a, beta_inverse), 4)
    dc_da = field.div(field.scale(4, c2_value), a)
    dc_dbeta = field.neg(field.div(field.scale(4, c2_value), beta))
    c1 = [
        field.scale(-8, a),
        field.zero,
        field.zero,
        field.scale(2 if orientation == "I" else -2, y),
        field.scale(-2 if orientation == "I" else 2, r),
    ]
    if orientation == "I":
        c2 = [
            field.sub(
                field.scale(8, field.mul(a, field.sub(field.mul(r2, y2), field.const(16)))),
                field.mul(dc_da, r2),
            ),
            field.zero,
            field.neg(field.mul(dc_dbeta, r2)),
            field.scale(2, ff_product(field, A, r2, y)),
            field.scale(2, field.mul(r, field.sub(field.mul(A, y2), c2_value))),
        ]
    else:
        y4 = field.square(y2)
        c2 = [
            field.sub(
                field.scale(8, field.mul(a, field.sub(y4, field.const(16)))),
                field.mul(dc_da, y2),
            ),
            field.zero,
            field.neg(field.mul(dc_dbeta, y2)),
            field.scale(2, field.mul(y, field.sub(field.scale(2, field.mul(A, y2)), c2_value))),
            field.zero,
        ]
    return [relation_row(field), c1, c2]


def solution_text(
    field: FiniteField, a: int, b: int, y: int, r: int, rank: int
) -> dict[str, Any]:
    _, s = base_values(field, a)
    return {
        "a": field.text(a),
        "s": field.text(s),
        "b_or_beta": field.text(b),
        "y": field.text(y),
        "r": field.text(r),
        "total_jacobian_rank": rank,
    }


def low_counts(field: FiniteField) -> tuple[dict[str, int], dict[str, Any]]:
    counts = {"I": 0, "II": 0}
    regular_fibre = {"I": 0, "II": 0}
    total_regular = {"I": 0, "II": 0}
    witnesses: dict[str, Any] = {}
    for a in range(field.q):
        A, _ = base_values(field, a)
        if field.chi(A) != -1:
            continue
        for y in range(field.q):
            y2 = field.square(y)
            for r in range(field.q):
                PACER.tick()
                r2 = field.square(r)
                c1_i = field.sub(field.sub(y2, r2), A)
                c2_i = field.mul(A, field.sub(field.mul(r2, y2), field.const(16)))
                if c1_i == c2_i == 0:
                    counts["I"] += 1
                    rank = matrix_rank(field, low_jacobian(field, "I", a, y, r))
                    assert rank == 3
                    total_regular["I"] += 1
                    fibre_det = field.scale(
                        4,
                        ff_product(
                            field,
                            y,
                            r,
                            A,
                            field.add(A, field.scale(2, r2)),
                        ),
                    )
                    regular_fibre["I"] += fibre_det != 0
                    witnesses.setdefault("I", solution_text(field, a, 1, y, r, rank))
                c1_ii = field.sub(field.sub(r2, y2), A)
                c2_ii = field.mul(A, field.sub(field.square(y2), field.const(16)))
                if c1_ii == c2_ii == 0:
                    counts["II"] += 1
                    rank = matrix_rank(field, low_jacobian(field, "II", a, y, r))
                    assert rank == 3
                    total_regular["II"] += 1
                    fibre_det = field.scale(-8, ff_product(field, A, y2, y, r))
                    regular_fibre["II"] += fibre_det != 0
                    witnesses.setdefault("II", solution_text(field, a, 1, y, r, rank))
    formula_I = 0
    formula_II = 0
    epsilon = field.chi(field.const(-1))
    for a in range(field.q):
        A, _ = base_values(field, a)
        if field.chi(A) != -1:
            continue
        for u in range(field.q):
            if field.mul(u, field.add(A, u)) == field.const(16):
                formula_I += (1 + field.chi(u)) * (1 + field.chi(field.add(A, u)))
        formula_II += 2 * (1 + field.chi(field.add(A, field.const(4))))
        if epsilon == 1:
            formula_II += 2 * (1 + field.chi(field.sub(A, field.const(4))))
    assert counts == {"I": formula_I, "II": formula_II}
    return (
        {
            "normalized_unit_parts_beta": field.q - 1,
            "orientation_I_points": counts["I"] * (field.q - 1),
            "orientation_II_points": counts["II"] * (field.q - 1),
            "orientation_I_beta_suppressed_points": counts["I"],
            "orientation_II_beta_suppressed_points": counts["II"],
            "orientation_I_fibre_regular": regular_fibre["I"] * (field.q - 1),
            "orientation_II_fibre_regular": regular_fibre["II"] * (field.q - 1),
            "orientation_I_total_regular": total_regular["I"] * (field.q - 1),
            "orientation_II_total_regular": total_regular["II"] * (field.q - 1),
        },
        witnesses,
    )


def low_character_sum(field: FiniteField) -> dict[str, int]:
    epsilon = field.chi(field.const(-1))
    sum_A = 0
    sum_C = 0
    T = 0
    zeros_A = 0
    zeros_C = 0
    plus_square_or_zero = 0
    plus_nonzero_square = 0
    nonsquare_A = 0
    for a in range(field.q):
        A, _ = base_values(field, a)
        C = field.add(A, field.const(4))
        chi_A, chi_C = field.chi(A), field.chi(C)
        sum_A += chi_A
        sum_C += chi_C
        T += field.chi(field.mul(A, C))
        zeros_A += A == 0
        zeros_C += C == 0
        if chi_A == -1:
            nonsquare_A += 1
            plus_square_or_zero += chi_C >= 0
            plus_nonzero_square += chi_C == 1
    K_numerator = (
        field.q + sum_C - sum_A - T + zeros_C * (1 - epsilon) - 2 * zeros_A
    )
    M_numerator = (
        field.q + sum_C - sum_A - T + zeros_C * (epsilon - 1) - 2 * zeros_A
    )
    assert K_numerator % 4 == M_numerator % 4 == 0
    assert K_numerator // 4 == plus_square_or_zero
    assert M_numerator // 4 == plus_nonzero_square
    if field.p != 5:
        assert sum_A == sum_C == -1
        assert (T + 1) ** 2 <= 4 * field.q
    else:
        assert sum_A == -1 and sum_C == field.q - 1 and T == -2
        assert plus_square_or_zero == (field.q - 1) // 2
    return {
        "epsilon=chi(-1)": epsilon,
        "sum_chi(A)": sum_A,
        "sum_chi(A+4)": sum_C,
        "T=sum_chi(A(A+4))": T,
        "zeros_A": zeros_A,
        "zeros_A+4": zeros_C,
        "nonsquare_A": nonsquare_A,
        "K_plus_square_or_zero": plus_square_or_zero,
        "M_plus_nonzero_square": plus_nonzero_square,
        "K_formula_numerator": K_numerator,
        "M_formula_numerator": M_numerator,
    }


def unit_solution_count_formula(field: FiniteField, orientation: str) -> int:
    answer = 0
    for a in range(field.q):
        A, s = base_values(field, a)
        if field.chi(A) != -1:
            continue
        s2 = field.square(s)
        for u in range(field.q):
            weight = (1 + field.chi(u)) * (1 + field.chi(field.add(A, u)))
            if weight == 0:
                continue
            if orientation == "I":
                n = ff_sum(field, field.square(u), field.mul(A, u), field.const(-16))
                d = field.sub(field.mul(A, s2), u)
            else:
                n = field.sub(field.square(u), field.const(16))
                d = field.sub(field.mul(A, field.sub(s2, field.one)), u)
            if n == 0 and d == 0:
                multiplicity = field.q - 1
            elif n != 0 and d != 0:
                multiplicity = 1
            else:
                multiplicity = 0
            answer += weight * multiplicity
    return answer


def unit_counts(field: FiniteField) -> tuple[dict[str, int], dict[str, Any]]:
    counts = {"I": 0, "II": 0}
    regular = {"I": 0, "II": 0}
    witnesses: dict[str, Any] = {}
    for a in range(field.q):
        A, s = base_values(field, a)
        if field.chi(A) != -1:
            continue
        s2 = field.square(s)
        for b in range(1, field.q):
            B = field.scale(2, b)
            for y in range(field.q):
                y2 = field.square(y)
                for r in range(field.q):
                    PACER.tick()
                    r2 = field.square(r)
                    c1_i = field.sub(field.sub(y2, r2), A)
                    c2_i = ff_sum(
                        field,
                        field.mul(A, field.mul(r2, y2)),
                        field.neg(field.scale(16, ff_product(field, A, B, field.sub(r2, field.mul(A, s2))))),
                        field.scale(-16, A),
                    )
                    if c1_i == c2_i == 0:
                        counts["I"] += 1
                        rank = matrix_rank(field, unit_jacobian(field, "I", a, b, y, r))
                        regular["I"] += rank == 3
                        if rank == 3:
                            witnesses.setdefault("I", solution_text(field, a, b, y, r, rank))
                    c1_ii = field.sub(field.sub(r2, y2), A)
                    c2_ii = ff_sum(
                        field,
                        field.mul(A, field.square(y2)),
                        field.neg(field.scale(16, ff_product(field, A, B, field.sub(r2, field.mul(A, s2))))),
                        field.scale(-16, A),
                    )
                    if c1_ii == c2_ii == 0:
                        counts["II"] += 1
                        rank = matrix_rank(field, unit_jacobian(field, "II", a, b, y, r))
                        regular["II"] += rank == 3
                        if rank == 3:
                            witnesses.setdefault("II", solution_text(field, a, b, y, r, rank))
    assert counts["I"] == unit_solution_count_formula(field, "I")
    assert counts["II"] == unit_solution_count_formula(field, "II")
    assert regular["I"] or regular["II"]
    return (
        {
            "orientation_I_points": counts["I"],
            "orientation_II_points": counts["II"],
            "orientation_I_total_regular": regular["I"],
            "orientation_II_total_regular": regular["II"],
        },
        witnesses,
    )


def fourth_root_count(field: FiniteField, value: int) -> int:
    return sum(field.pow(g, 4) == value for g in range(1, field.q))


def critical_solution_count_formula(field: FiniteField, orientation: str) -> int:
    answer = 0
    for a in range(field.q):
        A, _ = base_values(field, a)
        if field.chi(A) != -1:
            continue
        for u in range(field.q):
            weight = (1 + field.chi(u)) * (1 + field.chi(field.add(A, u)))
            if weight == 0:
                continue
            assert u != 0
            if orientation == "I":
                numerator = field.mul(
                    A,
                    ff_sum(field, field.square(u), field.mul(A, u), field.const(-16)),
                )
            else:
                numerator = field.mul(A, field.sub(field.square(u), field.const(16)))
            target = field.div(numerator, u)
            answer += weight * fourth_root_count(field, target)
    return answer


def critical_counts(field: FiniteField) -> tuple[dict[str, int], dict[str, Any]]:
    counts = {"I": 0, "II": 0}
    regular = {"I": 0, "II": 0}
    witnesses: dict[str, Any] = {}
    for a in range(field.q):
        A, _ = base_values(field, a)
        if field.chi(A) != -1:
            continue
        for beta in range(1, field.q):
            c2_value = field.pow(field.div(a, beta), 4)
            for y in range(field.q):
                y2 = field.square(y)
                for r in range(field.q):
                    PACER.tick()
                    r2 = field.square(r)
                    c1_i = field.sub(field.sub(y2, r2), A)
                    c2_i = ff_sum(
                        field,
                        field.mul(A, field.mul(r2, y2)),
                        field.neg(field.mul(c2_value, r2)),
                        field.scale(-16, A),
                    )
                    if c1_i == c2_i == 0:
                        counts["I"] += 1
                        rank = matrix_rank(field, critical_jacobian(field, "I", a, beta, y, r))
                        regular["I"] += rank == 3
                        if rank == 3:
                            witnesses.setdefault("I", solution_text(field, a, beta, y, r, rank))
                    c1_ii = field.sub(field.sub(r2, y2), A)
                    c2_ii = ff_sum(
                        field,
                        field.mul(A, field.square(y2)),
                        field.neg(field.mul(c2_value, y2)),
                        field.scale(-16, A),
                    )
                    if c1_ii == c2_ii == 0:
                        counts["II"] += 1
                        rank = matrix_rank(field, critical_jacobian(field, "II", a, beta, y, r))
                        regular["II"] += rank == 3
                        if rank == 3:
                            witnesses.setdefault("II", solution_text(field, a, beta, y, r, rank))
    assert counts["I"] == critical_solution_count_formula(field, "I")
    assert counts["II"] == critical_solution_count_formula(field, "II")
    if field.q == 5:
        assert counts == {"I": 0, "II": 0}
    return (
        {
            "orientation_I_points": counts["I"],
            "orientation_II_points": counts["II"],
            "orientation_I_total_regular": regular["I"],
            "orientation_II_total_regular": regular["II"],
        },
        witnesses,
    )


def small_field_rows() -> list[dict[str, Any]]:
    expected = {
        3: {"K": 2, "M": 0, "low": (0, 4), "unit": (4, 4), "critical": (8, 0)},
        5: {"K": 2, "M": 2, "low": (0, 8), "unit": (8, 0), "critical": (0, 0)},
        7: {"K": 4, "M": 2, "low": (0, 12), "unit": (24, 10), "critical": (32, 16)},
        9: {"K": 0, "M": 0, "low": (0, 16), "unit": (32, 16), "critical": (64, 0)},
        11: {"K": 2, "M": 2, "low": (8, 8), "unit": (52, 52), "critical": (56, 48)},
        13: {"K": 2, "M": 2, "low": (32, 16), "unit": (40, 56), "critical": (32, 128)},
    }
    rows = []
    for field in SMALL_FIELDS:
        low, low_witnesses = low_counts(field)
        character = low_character_sum(field)
        unit, unit_witnesses = unit_counts(field)
        critical, critical_witnesses = critical_counts(field)
        assert low["orientation_II_total_regular"] > 0
        certificate = expected[field.q]
        assert (
            character["K_plus_square_or_zero"],
            character["M_plus_nonzero_square"],
        ) == (certificate["K"], certificate["M"])
        assert (
            low["orientation_I_beta_suppressed_points"],
            low["orientation_II_beta_suppressed_points"],
        ) == certificate["low"]
        assert (unit["orientation_I_points"], unit["orientation_II_points"]) == certificate["unit"]
        assert (
            critical["orientation_I_points"],
            critical["orientation_II_points"],
        ) == certificate["critical"]
        rows.append(
            {
                "type": "small-finite-field",
                "label": "PROVED",
                "field": field.name,
                "q": field.q,
                "scope": "exact exhaustive finite case required outside the general elementary/Hasse bounds",
                "low_positive_k<t": low,
                "low_positive_witnesses": low_witnesses,
                "character_sum": character,
                "unit_k=0": unit,
                "unit_witnesses": unit_witnesses,
                "critical_k=t": critical,
                "critical_witnesses": critical_witnesses,
            }
        )
    return rows


def vp(value: F | int, prime: int) -> int:
    value = F(value)
    if value == 0:
        return INF
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
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, -1, prime) % prime


def bridge_value(a: int, b: int, Z: int) -> F:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    assert b and D
    return F(a * a * Z * Z * Ng, A * b * b * D)


def valuation_rows() -> list[dict[str, Any]]:
    choices = {3: 5, 5: 7, 7: 1}
    rows = []
    for prime, a in choices.items():
        A = 1 + 4 * a * a
        assert pow(A % prime, (prime - 1) // 2, prime) == prime - 1
        for e in (1, 2):
            Z = prime ** (3 * e)
            t = 3 * e
            D = 1 - Z - a * a * Z * Z
            assert vp(D, prime) == 0
            cases = []
            for k in (0, 1, t, t + 1):
                b = prime**k
                c = bridge_value(a, b, Z)
                if k == 0:
                    assert vp(c, prime) >= 2 * t
                    expected = ">=2t"
                    residue = 0
                else:
                    assert vp(c, prime) == 2 * (t - k)
                    expected = "2(t-k)"
                    residue = mod_fraction(c, prime) if k == t else None
                    if k == t:
                        assert residue == (-a * a) % prime
                cases.append(
                    {
                        "k=v_w(b)": k,
                        "v_w(c)": vp(c, prime),
                        "formula": expected,
                        "c_residue_if_unit": residue,
                    }
                )
            rows.append(
                {
                    "type": "valuation-replay",
                    "label": "PROVED",
                    "w": prime,
                    "e=v_w(z)": e,
                    "t=v_w(Z)": t,
                    "a": a,
                    "A": A,
                    "D_unit": True,
                    "cases": cases,
                }
            )
    return rows


def critical_counterexample_row() -> dict[str, Any]:
    w = 5
    z = 5
    Z = z**3
    a = 7
    s = 3
    b = 5**3
    A = 1 + 4 * a * a
    c = bridge_value(a, b, Z)
    assert a == 1 + 2 * s
    assert vp(s, 2) >= 0 and vp(b, 2) == 0
    assert vp(Z, w) == vp(b, w) == 3 and vp(c, w) == 0
    assert A % w == 2 and mod_fraction(c * c, w) == 1
    # I: 2u^2+3u+3, discriminant 0, sole root 3 (a nonsquare).
    squares = {x * x % w for x in range(w)}
    roots_I = [u for u in range(w) if (2 * u * u + 3 * u + 3) % w == 0]
    roots_II = [u for u in range(w) if (2 * u * u + 4 * u + 3) % w == 0]
    assert roots_I == [3] and 3 not in squares and roots_II == []
    return {
        "type": "critical-stratum-obstruction",
        "label": "PROVED",
        "cell": {"w": w, "z": z, "Z": Z, "t=v_w(Z)": 3},
        "base": {"a": a, "s": s, "b": b, "k=v_w(b)": 3, "A": A},
        "Phi_compatibility_used": {"v2(s)>=0": True, "v2(b)=0": True},
        "bridge": {"v_w(c)": 0, "c^2_mod_w": 1},
        "orientation_I_reduction": "2*u^2+3*u+3=0; sole root u=3 is nonsquare",
        "orientation_II_reduction": "2*u^2+4*u+3=0; discriminant 2 is nonsquare",
        "conclusion": "neither orientation has a Q_5 point over this guarded base",
        "scope_warning": "the two displayed 2-adic conditions are necessary Phi compatibility, not a proof that this exact rational base is a global Phi member",
    }


def dyadic_supersession_row() -> dict[str, Any]:
    orientation_I_checks = 0
    orientation_II_checks = 0
    odd_residues = range(1, 16, 2)
    for s in range(16):
        a = 1 + 2 * s
        A = 1 + 4 * a * a
        assert A % 32 == 5
        for t in odd_residues:
            # The two possible I-root valuations are zero and four.
            assert (A + t * t) % 8 == 6
            assert (A + 16 * t * t) % 8 == 5
            orientation_I_checks += 2
            for b in odd_residues:
                for h in odd_residues:
                    # II after u=4t^2 and division by 16:
                    # t^4-8ht^2+[2Ab(s^2-1)-1].
                    residue = (
                        t**4 - 8 * h * t * t + 2 * A * b * (s * s - 1) - 1
                    ) % 16
                    assert residue != 0
                    orientation_II_checks += 1
                    PACER.tick()
    assert orientation_I_checks == 256
    assert orientation_II_checks == 8192
    return {
        "type": "dyadic-global-supersession",
        "label": "PROVED",
        "source": [
            "math/h10q/l24_diagonal_geometry.py",
            "math/h10q/data/l24_diagonal_geometry.jsonl",
            "/tmp/l24_diagonal_geometry.md",
        ],
        "bridge_bound": "on Phi, v2(c)>=4: v2(Ng)>=4 and v2(D)=0 for v2(Z)>=0, while v2(D)=2v2(Z) for v2(Z)<0",
        "orientation_I": "Q_I/A has coefficient valuations (4,0,0); roots have v2(u)=4 or 0, and A+u is respectively 5 or 6 mod 8, never a square",
        "orientation_II": "Q_II/A=u^2-32h*u+16e with h,e odd; roots have v2(u)=2. Writing u=4t^2 gives mod 16 the impossible congruence 8+2Ab(s^2-1)=0",
        "parity_split": "v2(s^2-1)=0 for s even and >=3 for s odd, so the last congruence has valuation 1 or residue 8 mod 16",
        "exact_modular_checks": {
            "orientation_I": orientation_I_checks,
            "orientation_II": orientation_II_checks,
        },
        "consequence": "both diagonal orientations are Phi-empty over Q_2; the odd target-place theorem cannot produce a global five-unknown definition",
    }


def theorem_rows() -> list[dict[str, Any]]:
    return [
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "bridge-valuation-trichotomy",
            "hypotheses": "odd target DVR; t=v(Z)>0; a,A,D units; (A|kappa)=-1; k=v(b)>=0",
            "result": {
                "k=0": "v(c)>=2t, hence Bbar!=0 and cbar=0",
                "0<k<t": "v(c)=2(t-k)>0, hence Bbar=cbar=0",
                "k=t": "cbar=-a_bar^2*zeta^2/beta^2 and cbar^2 is a nonzero fourth power",
                "k>t": "v(c)=2(t-k)<0",
            },
            "proof": "for k>0, Ng=16a^4b^2-A(b-1)^4 is congruent to -A; for k=0 its valuation is nonnegative",
        },
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "exact-eliminants-and-fibre-Jacobians",
            "orientation_I": {
                "u": "r^2",
                "squares": ["u", "A+u=y^2"],
                "Q": "A*u^2+(A^2-c^2-16*A*B)*u+16*A^2*B*s^2-16*A",
                "det_d(C1,C2)/d(y,r)": "4*y*r*(A*(y^2+r^2)-c^2-16*A*B)=4*y*r*Q_I'(u)",
            },
            "orientation_II": {
                "u": "y^2",
                "squares": ["u", "A+u=r^2"],
                "Q": "A*u^2-(c^2+16*A*B)*u+16*A^2*B*(s^2-1)-16*A",
                "det_d(C1',C2')/d(y,r)": "4*y*r*(c^2+16*A*B-2*A*u)=-4*y*r*Q_II'(u)",
            },
            "Hensel": "rank two in a chosen pair of variables over the residue field is sufficient; total-space rank may use s or b when the fibre ramifies",
        },
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "low-positive-uniform-local-solubility",
            "hypotheses": "odd residue field F_q; 0<k<t; A=1+4a^2 is a nonsquare",
            "reduction": {
                "I": "u*(A+u)=16",
                "II": "u^2=16, i.e. u=+4 or u=-4",
            },
            "exact_counts": {
                "I": "(q-1)*sum_{chi(A)=-1,u}(1+chi(u))(1+chi(A+u))*1_{u(A+u)=16}",
                "II": "2*(q-1)*sum_{chi(A)=-1}[(1+chi(A+4))+1_{chi(-1)=1}(1+chi(A-4))]",
                "normalization": "q-1 choices of beta=b/pi^k are included",
            },
            "orientation_II_plus_arm": "u=4, y=+-2, r^2=A+4",
            "character_sum": "K_q=[q+S_C-S_A-T+n_C(1-epsilon)-2n_A]/4 for A+4 square or zero; M_q=[q+S_C-S_A-T+n_C(epsilon-1)-2n_A]/4 for A+4 nonzero square",
            "curve": "T=sum_a chi((1+4a^2)(5+4a^2)); if char!=5, #E(F_q)=q+T+2 for E:v^2=(1+4x^2)(5+4x^2), so |T+1|<=2sqrt(q)",
            "finite_cases": [3, 5, 7, 9],
            "conclusion": "the plus arm exists for every odd q except F_9; F_9 has a fibre-regular minus-arm point. For prime w>=5 a nonzero-square A+4 may be chosen, giving fibrewise Hensel; w=3 is lifted on the s-chart.",
            "target_consequence": "t=3*v_w(z)>=3, so k=1 always lies in this theorem",
        },
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "unit-residue-system-and-count",
            "hypotheses": "k=0, so B=2b is nonzero and cbar=0",
            "orientation_I": "n_I=u^2+A*u-16, d_I=A*s^2-u, n_I+16*B*d_I=0",
            "orientation_II": "n_II=u^2-16, d_II=A*(s^2-1)-u, n_II+16*B*d_II=0",
            "exact_count": "N_O=sum_{chi(A)=-1,u}(1+chi(u))(1+chi(A+u))*L_q(n_O,d_O), where L=q-1 if n=d=0, L=1 if n*d!=0, and L=0 otherwise",
            "existence": "for q>=17, the orientation-II hyperbola has q-1 ordered points and at most 12 lie over n*d=0; q=3,5,7,9,11,13 are exhausted exactly",
            "Hensel": "when n*d!=0, dC2/db=32*A*d is a unit and one nonzero y/r derivative of C1 supplies rank two",
        },
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "critical-residue-system",
            "hypotheses": "k=t; write Z=pi^t*zeta and b=pi^t*beta",
            "bridge": "cbar=-g^2 with g=a*zeta/beta, hence cbar^2=g^4 and g ranges F_q^* with beta",
            "counts": {
                "I": "sum (1+chi(u))(1+chi(A+u))*nu_4(A*(u^2+A*u-16)/u)",
                "II": "sum (1+chi(u))(1+chi(A+u))*nu_4(A*(u^2-16)/u)",
                "nu_4": "#{g in F_q^*:g^4=x}",
            },
            "integrality": "unit leading and constant coefficients force every local root u to be a unit, so residue square tests are necessary",
            "obstruction": "q=5 has count zero in both orientations on the entire k=t, nonsquare-A residue stratum",
        },
        {
            "type": "theorem",
            "label": "PROVED",
            "name": "pole-range-obstruction",
            "hypotheses": "k>t and A is a nonsquare unit",
            "newton": "each eliminant has root valuations +4(k-t) and -4(k-t)",
            "proof": "a positive-valuation square u makes A+u reduce to nonsquare A; a negative-valuation root satisfies u/(c^2/A)=1 mod pi, so u has nonsquare class A^{-1}",
            "conclusion": "neither orientation has a local point",
        },
    ]


def render_report(rows: list[dict[str, Any]], elapsed: float) -> None:
    small = [row for row in rows if row["type"] == "small-finite-field"]
    table = []
    for row in small:
        low = row["low_positive_k<t"]
        unit = row["unit_k=0"]
        critical = row["critical_k=t"]
        chars = row["character_sum"]
        table.append(
            f"| {row['q']} | {chars['T=sum_chi(A(A+4))']} | {chars['K_plus_square_or_zero']} | "
            f"{chars['M_plus_nonzero_square']} | {low['orientation_I_points']} / {low['orientation_II_points']} | "
            f"{unit['orientation_I_points']} / {unit['orientation_II_points']} | "
            f"{critical['orientation_I_points']} / {critical['orientation_II_points']} |"
        )
    lines = [
        "# L24 - target-place solvability of the diagonal branches",
        "",
        "**Verdict.** The odd target-place problem itself is **PROVED uniformly soluble** on the valuation stratum used by the standard construction: if `e=v_w(z)>=1`, take `k=v_w(b)=1`; then `t=v_w(Z)=3e`, so `0<k<t`, and orientation II has a nonsingular total-space local point for every odd `w`. For `w>=5` one can choose a fibre-regular `u=4` point and keep the selected base fixed. At `w=3` the only `u=4` points have `r=0`; the fibre Jacobian drops rank, but the `s` derivative restores full rank and multivariable Hensel lifts the point. **This does not rescue the five-unknown route:** the independent dyadic audit, reproduced in Section 7, proves both orientations Phi-empty over `Q_2`.",
        "",
        "Not every target valuation of `b` works: the exact critical stratum `k=t` has a **PROVED common obstruction at `w=5`**, and the pole stratum `k>t` is **PROVED empty in both orientations** whenever `A` is the required nonsquare unit. The explicit rational critical base below satisfies the two stated 2-adic Phi-compatibility conditions, but global membership of that exact base in `Phi` is not asserted. More decisively, the uniform dyadic obstruction applies to every Phi base. Thus the odd target-place obstruction is removed for the existential standard choice `k=1`, while global diagonal completeness is **PROVED FALSE**.",
        "",
        "## 1. Setup and bridge valuations",
        "",
        "Let `K` be a complete discretely valued field of odd residue field `F_q`, uniformizer `pi`, and put",
        "",
        "```text",
        "a=1+2s, A=1+4a^2, B=2b, Z=z^3, D=1-Z-a^2Z^2,",
        "Ng=16a^4b^2-A(b-1)^4, c=a^2 Z^2 Ng/(A b^2 D).",
        "```",
        "",
        "At a target place `t=v(Z)>0`; the selected `A` is a nonsquare residue unit, hence `a,A,D` are units and `D=1 mod pi`. If `k=v(b)=0`, then `v(c)>=2t`, so `c=0` in the residue field. If `k>0`, then `Ng=-A mod pi` and therefore",
        "",
        "```text",
        "v(c)=2(t-k).",
        "```",
        "",
        "Consequently: `0<k<t` gives `B=c=0`; `k=t`, with `Z=pi^t zeta` and `b=pi^t beta`, gives `c=-a^2 zeta^2/beta^2 mod pi`, so `c^2` is a nonzero fourth power; and `k>t` makes `c` nonintegral. These are equalities, not heuristic leading terms.",
        "",
        "## 2. Exact eliminants and Jacobians",
        "",
        "For orientation I set `u=r^2`, so `y^2=A+u`. For orientation II set `u=y^2`, so `r^2=A+u`. In both cases `u` and `A+u` must be squares. Direct elimination gives",
        "",
        "```text",
        "Q_I  = A u^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A,",
        "Q_II = A u^2-(c^2+16AB)u+16A^2B(s^2-1)-16A.",
        "```",
        "",
        "For the two fibre variables `(y,r)`, the exact determinants are",
        "",
        "```text",
        "det J_I  =  4yr[A(y^2+r^2)-c^2-16AB] =  4yr Q_I'(u),",
        "det J_II =  4yr[c^2+16AB-2Au]        = -4yr Q_II'(u).",
        "```",
        "",
        "A nonzero determinant gives ordinary Hensel with the base fixed. At a ramified square cover the correct total-space test also includes `E=a-1-2s`; the replay forms the full `3 x 5` Jacobian in `(a,s,b,y,r)` (or the normalized unit `beta` on `k=t`) and checks rank three.",
        "",
        "## 3. The uniform low-positive theorem (`0<k<t`)",
        "",
        "Here `B=c=0`. Orientation I reduces to `u(A+u)=16`. Orientation II reduces to",
        "",
        "```text",
        "u^2=16,  hence u=+4 or u=-4.",
        "```",
        "The exact normalized `(a,s,beta,y,r)` counts (with `beta=b/pi^k`) are",
        "",
        "```text",
        "N_I =(q-1) sum_{chi(A)=-1,u}(1+chi(u))(1+chi(A+u)) 1_{u(A+u)=16},",
        "N_II=2(q-1) sum_{chi(A)=-1}[(1+chi(A+4))",
        "                    +1_{chi(-1)=1}(1+chi(A-4))].",
        "```",
        "",
        "",
        "The plus arm always has `y=+-2` and asks only that `A+4=5+4a^2` be a square, zero allowed. Write `chi` for the quadratic character, `epsilon=chi(-1)`, `n_A=#{A=0}`, `n_C=#{A+4=0}`, and",
        "",
        "```text",
        "S_A=sum chi(A), S_C=sum chi(A+4), T=sum chi(A(A+4)).",
        "K_q=#{a:chi(A)=-1, chi(A+4)>=0}",
        "   =[q+S_C-S_A-T+n_C(1-epsilon)-2n_A]/4.          (1)",
        "M_q=#{a:chi(A)=-1, chi(A+4)=1}",
        "   =[q+S_C-S_A-T+n_C(epsilon-1)-2n_A]/4.          (2)",
        "```",
        "",
        "For characteristic other than 5, `S_A=S_C=-1`. The smooth quartic genus-one curve",
        "",
        "```text",
        "E: v^2=(1+4x^2)(5+4x^2)",
        "```",
        "",
        "has two rational points at infinity, so `#E(F_q)=q+T+2` and Hasse gives `|T+1|<=2 sqrt(q)`. If `epsilon=-1`, (1) is already strictly positive. If `epsilon=+1`, it is positive for `q>9`. In characteristic 5 the quartic degenerates harmlessly and direct algebra gives `K_q=M_q=(q-1)/2`. Exact exhaustion gives `K_9=M_9=0`, but the `u=-4` arm over `F_9` has 16 fibre-regular `(a,s,y,r)` points and hence 128 normalized points after the eight choices of `beta`. The same table handles `F_3,F_5,F_7`. Thus orientation II has a full-rank point over every odd finite field. For prime fields `w>=5`, (2) and the same Hasse bound (plus exact `w=5,7`) give `M_w>0`, so the base can be fixed; only `w=3` needs the total-space chart.",
        "",
        "At every plus-arm point, `Q_II'(4)=8A` is nonzero. If `r!=0`, `det J_II` is nonzero. If `r=0`, the fibre determinant vanishes, but `a` is nonzero (since `A` is nonsquare), and the two rows have independent entries `dC1/ds=-16a` and `dC2/dy=4Ay^3` after substituting `a=1+2s`. This includes every singular `A+4=0` branch. For the actual target `t=3e`, choosing `k=1` proves the claimed uniform total-space local point.",
        "",
        "The minus arm is also exact: it requires `chi(-1)=+1` and `A-4` square. It is needed for the `F_9` extension-field exception, but not for any prime target place; both `u=+4` and `u=-4` are included in every finite replay count.",
        "",
        "## 4. Unit `b` (`k=0`)",
        "",
        "Now `B=2b` is nonzero and `c=0`. The two systems, linear in `B`, are",
        "",
        "```text",
        "I:  n_I+16B d_I=0, n_I=u^2+Au-16, d_I=As^2-u,",
        "II: n_II+16B d_II=0, n_II=u^2-16, d_II=A(s^2-1)-u.",
        "```",
        "",
        "Define `L_q(n,d)=q-1` if `n=d=0`, `1` if `nd!=0`, and `0` otherwise. The exact ordered `(a,b,s,y,r)` count is",
        "",
        "```text",
        "N_O=sum_{chi(A)=-1,u}(1+chi(u))(1+chi(A+u)) L_q(n_O,d_O).   (2)",
        "```",
        "",
        "For existence when `q>=17`, use orientation II. For any nonsquare `A`, the hyperbola `r^2-y^2=A` has exactly `q-1` ordered points. The two roots of `n_II` and the one root of `d_II` account for at most twelve ordered points, so some point has `n_II d_II!=0`; it determines a unique nonzero `B`. Moreover `dC2/db=32A d_II` is a unit and a `y` or `r` derivative of `C1` is nonzero, giving full Hensel rank. Formula (2), direct tuple enumeration, and the Jacobian agree for `q=3,5,7,9,11,13`. Orientation I supplies the small `q=5` unit point where orientation II has none.",
        "",
        "## 5. Critical and pole ranges",
        "",
        "At `k=t`, put `g=a zeta/beta`. Since `beta` varies over the units, so does `g`, and `c^2=g^4`. If `nu_4(x)=#{g in F_q^*:g^4=x}`, the exact counts are",
        "",
        "```text",
        "N_I^crit  = sum (1+chi(u))(1+chi(A+u)) nu_4(A(u^2+Au-16)/u),",
        "N_II^crit = sum (1+chi(u))(1+chi(A+u)) nu_4(A(u^2-16)/u).       (3)",
        "```",
        "",
        "The weights force `u!=0`, so (3) has no hidden division. The critical eliminants have integral coefficients with unit leading and constant terms; any local root `u` must therefore be a unit, making the residue square tests necessary. Over `F_5`, nonsquare `A=1+4a^2` forces `A=2`, while every nonzero fourth power is `c^2=1`. Orientation I becomes `2u^2+3u+3`; its sole root is the nonsquare `u=3`. Orientation II becomes `2u^2+4u+3`; its discriminant is the nonsquare `2`. Thus both orientations have zero points on this critical stratum. The exact rational base `z=5,a=7,s=3,b=125` has `t=k=3`, `v_5(c)=0`, and `c^2=1 mod 5`, so reduction proves no `Q_5` point over that base. It satisfies `v_2(s)>=0` and `v_2(b)=0`; no stronger global-Phi claim is made, and this is not a whole-cell obstruction because another valuation of `b` may be chosen.",
        "",
        "If `k>t`, put `d=k-t>0`. In either eliminant the constant and leading coefficients are units and the linear coefficient has valuation `-4d`, so the two Newton slopes give root valuations `+4d` and `-4d`. A positive-valuation square `u` would make `A+u` reduce to the nonsquare `A`. For the negative root, division by the two dominant terms gives `u/(c^2/A)=1 mod pi`; hence `u` has square class `A^{-1}`, again nonsquare. Both orientations are therefore empty throughout the pole range, including nonintegral possible hyperbola coordinates.",
        "",
        "## 6. Exact finite fields",
        "",
        "The table entries are exact ordered-point counts `I / II`. In the low-positive column they include all `q-1` normalized unit parts `beta=b/pi^k`; the unit and critical columns already enumerate `b` or `beta` directly. `K` allows the ramified value `A+4=0`; `M` requires a nonzero square and hence a fibre-regular plus arm.",
        "",
        "| q | T | K | M | low `0<k<t` | unit `k=0` | critical `k=t` |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        *table,
        "",
        "Every listed count is independently obtained both from the character formula and direct tuple equations; every recorded positive certificate is checked against the full Jacobian. The `q=5` critical zero is then proved symbolically above, not inferred from an unbounded no-hit search.",
        "",
        "## 7. Dyadic supersession: both global branches are empty",
        "",
        "The independent geometry artifact `math/h10q/l24_diagonal_geometry.py` established the uniform bridge bound `v_2(c)>=4`; this replay also checks the terminal congruences. Indeed `a,b,A` are 2-adic units, `v_2(Ng)>=4`, and",
        "",
        "```text",
        "v_2(D)=0       if v_2(Z)>=0,",
        "v_2(D)=2v_2(Z) if v_2(Z)<0,",
        "```",
        "",
        "so the displayed bound follows directly from the bridge formula. For orientation I, `Q_I/A` has coefficient valuations `(4,0,0)`, hence possible root valuations four and zero. If the square `u=r^2` has valuation four, `A+u=5 mod 8`; if it is a unit square, `A+u=6 mod 8`. Neither is a `Q_2` square.",
        "",
        "For orientation II write exactly",
        "",
        "```text",
        "Q_II/A = u^2-32h*u+16e,",
        "h=b+c^2/(32A) in Z_2^*,  e=2Ab(s^2-1)-1 in Z_2^*.",
        "```",
        "",
        "Its coefficient valuations `(4,5,0)` force `v_2(u)=2`. If the required `u=y^2` exists, write `u=4t^2` with `t` odd. Division by 16 and reduction modulo 16 gives",
        "",
        "```text",
        "0 = t^4-8h*t^2+e = 8+2Ab(s^2-1)  (mod 16).",
        "```",
        "",
        "For even `s`, `v_2(s^2-1)=0`, so the second term has valuation one; for odd `s`, `v_2(s^2-1)>=3`, so the residue is eight. Both are impossible. The script exhausts all 8,192 residue tuples `(s,b,t,h) mod 16` as an exact replay, not as evidence. Therefore both diagonal orientations are **PROVED Phi-empty**, and their union cannot define the target set.",
        "",
        "## 8. Scope",
        "",
        "Labels: valuation formulas, eliminants, Jacobians, character identities, the Hasse deduction, all finite exceptional cases, low-positive uniform local solubility, the `F_5` critical obstruction, the pole obstruction, and the dyadic global obstruction are **PROVED**. The odd-place theorem remains an exact local classification, but the dyadic theorem supersedes it for completeness: the proposed diagonal five-unknown route is **PROVED empty**, not open.",
        "",
        f"Verification wall-clock: **{elapsed:.3f} s**. Pacing: {PACER.steps} inner-loop steps and {PACER.sleeps} sleeps of {PACE_SECONDS} s.",
        "",
        "Artifacts: `math/h10q/l24_diagonal_local.py`, `math/h10q/data/l24_diagonal_local.jsonl`, `/tmp/l24_diagonal_local.md`.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    rows: list[dict[str, Any]] = [
        {
            "type": "meta",
            "artifact": "l24_diagonal_local",
            "labels": ["PROVED", "OPEN"],
            "stdlib_only": True,
            "finite_scan_policy": "only q=3,5,7,9,11,13, exactly the cases outside the stated general bounds",
            "target": "local solvability of both diagonal orientations",
        },
        *theorem_rows(),
        dyadic_supersession_row(),
        *valuation_rows(),
        *small_field_rows(),
        critical_counterexample_row(),
    ]
    elapsed = time.perf_counter() - started
    summary = {
        "type": "summary",
        "label": "PROVED",
        "rows": len(rows) + 1,
        "uniform_target_local_theorem": "orientation II is locally soluble for k=1<t=3*v_w(z), every odd target w",
        "unit_branch_uniform": True,
        "critical_stratum_common_obstruction": {"w": 5, "k=t": 3},
        "pole_range_both_empty": True,
        "global_diagonal_completeness": "PROVED EMPTY: both orientations have a uniform Phi-compatible Q_2 obstruction; the odd-place theorem is only local",
        "elapsed_seconds": elapsed,
        "pacing": {
            "steps": PACER.steps,
            "sleeps": PACER.sleeps,
            "sleep_seconds": PACE_SECONDS,
        },
    }
    rows.append(summary)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    render_report(rows, elapsed)
    print(
        f"l24_diagonal_local: rows={len(rows)} small_fields={len(SMALL_FIELDS)} "
        f"steps={PACER.steps} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
