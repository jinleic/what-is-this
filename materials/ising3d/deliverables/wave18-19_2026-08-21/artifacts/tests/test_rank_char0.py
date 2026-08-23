"""Standalone clean-room verification of the exact characteristic-zero T counts.

The verifier imports none of e166/e167/e168.  It rebuilds R_T by an explicit
Fraction triple product and clears the actual entry-denominator lcm.  It then
counts product values with companion, exterior-square, square, and Kronecker
matrices -- a different construction from the producer's resultants.

It also rebuilds from raw graph/parameter inputs the F_1000003 pair polynomials
for T disjoint union P_2 at all four named points and T disjoint union P_4 at
t=1/3.  Those one-sided characteristic-zero lower bounds are squeezed against
the independently recomputed exact core upper bounds.

Run from the repository root with:
    .venv/bin/python tests/test_rank_char0.py
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import resource
import sys
import time
from fractions import Fraction
from math import comb, gcd
from pathlib import Path

# Keep this standalone audit polite on the contended host.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "rank_char0.json"
T_BONDS: tuple[tuple[int, int], ...] = ((0, 1), (0, 2), (0, 3), (3, 4))
NAMED_T: tuple[Fraction, ...] = (
    Fraction(1, 3),
    Fraction(1, 4),
    Fraction(2, 5),
    Fraction(1, 5),
)
PRIME = 1_000_003
CPU_BUDGET_SECONDS = 1_200.0
RSS_LIMIT_BYTES = 4_000_000_000
Z = sp.Symbol("z")
PASSED: list[str] = []
FAILED: list[str] = []


class NonDecisiveBudget(RuntimeError):
    pass


class CpuBudget:
    def __init__(self, limit: float = CPU_BUDGET_SECONDS) -> None:
        self.limit = limit
        self.started = time.process_time()

    def check(self, stage: str) -> None:
        used = time.process_time() - self.started
        if used > self.limit:
            raise NonDecisiveBudget(
                f"NON-DECISIVE: process-CPU budget {self.limit:.1f}s expired "
                f"at {stage} after {used:.3f}s"
            )
        rss = max_rss_bytes()
        if rss > RSS_LIMIT_BYTES:
            raise NonDecisiveBudget(
                f"NON-DECISIVE: RSS budget {RSS_LIMIT_BYTES} exceeded at {stage}; "
                f"ru_maxrss={rss}"
            )


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(f"[{'PASS' if condition else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    (PASSED if condition else FAILED).append(name)


def text_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def physical_weight(t: Fraction) -> Fraction:
    return (1 + t * t) / (2 * t)


# ---------------------------------------------------------------- exact independent core

def build_rational_operator(t: Fraction, w: Fraction, budget: CpuBudget) -> list[list[Fraction]]:
    """Build R=PDP directly from its entries, without an integer Kronecker sweep."""

    dimension = 32
    p_matrix = [[Fraction(0)] * dimension for _ in range(dimension)]
    for left in range(dimension):
        if left % 4 == 0:
            budget.check("direct rational P")
        for right in range(dimension):
            p_matrix[left][right] = t ** ((left ^ right).bit_count())

    diagonal: list[Fraction] = []
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (4 - site)) & 1) for site in range(5)]
        aligned = sum(spins[left] == spins[right] for left, right in T_BONDS)
        diagonal.append(w**aligned)

    result = [[Fraction(0)] * dimension for _ in range(dimension)]
    for left in range(dimension):
        budget.check(f"direct rational PDP row {left}")
        for right in range(left, dimension):
            value = sum(
                (p_matrix[left][middle] * diagonal[middle] * p_matrix[middle][right]
                 for middle in range(dimension)),
                Fraction(0),
            )
            result[left][right] = result[right][left] = value
    return result


def clear_actual_lcm(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = scale * value.denominator // gcd(scale, value.denominator)
    integer = [[int(value * scale) for value in row] for row in matrix]
    assert all(value * scale == integer[i][j]
               for i, row in enumerate(matrix) for j, value in enumerate(row))
    return scale, integer


def companion_matrix(factor: sp.Poly) -> sp.Matrix:
    degree = factor.degree()
    coefficients = factor.all_coeffs()
    matrix = sp.zeros(degree)
    for row in range(1, degree):
        matrix[row, row - 1] = 1
    for row in range(degree):
        matrix[row, degree - 1] = -coefficients[degree - row]
    assert sp.Poly(matrix.charpoly(Z).as_expr(), Z, domain=sp.QQ) == factor
    return matrix


def exterior_square(matrix: sp.Matrix) -> sp.Matrix:
    pairs = [(left, right) for left in range(matrix.rows) for right in range(left + 1, matrix.rows)]
    result = sp.zeros(len(pairs))
    for column, (left, right) in enumerate(pairs):
        for row, (top, bottom) in enumerate(pairs):
            result[row, column] = (
                matrix[top, left] * matrix[bottom, right]
                - matrix[top, right] * matrix[bottom, left]
            )
    return result


def characteristic(poly_matrix: sp.Matrix) -> sp.Poly:
    return sp.Poly(poly_matrix.charpoly(Z).as_expr(), Z, domain=sp.QQ).monic()


def squarefree_lcm(polynomials: list[sp.Poly], budget: CpuBudget, stage: str) -> sp.Poly:
    union = sp.Poly(1, Z, domain=sp.QQ)
    for index, polynomial in enumerate(polynomials):
        budget.check(f"{stage} component {index}")
        squarefree = sp.Poly(polynomial.sqf_part(), Z, domain=sp.QQ).monic()
        overlap = sp.Poly(sp.gcd(union, squarefree), Z, domain=sp.QQ).monic()
        union = sp.Poly((union * squarefree).exquo(overlap), Z, domain=sp.QQ).monic()
    return union


def explicit_product_counts(
    factors: list[tuple[sp.Poly, int]],
    budget: CpuBudget,
) -> tuple[int, int]:
    companions = [companion_matrix(factor) for factor, _exponent in factors]
    r_components: list[sp.Poly] = []
    diagonal_only: list[sp.Poly] = []
    for index, ((factor, exponent), companion) in enumerate(zip(factors, companions)):
        budget.check(f"explicit factor representation {index}")
        if factor.degree() >= 2:
            r_components.append(characteristic(exterior_square(companion)))
        square_poly = characteristic(companion * companion)
        if exponent >= 2:
            r_components.append(square_poly)
        else:
            diagonal_only.append(square_poly)
    for left in range(len(companions)):
        for right in range(left + 1, len(companions)):
            budget.check(f"explicit Kronecker {left},{right}")
            r_components.append(characteristic(sp.kronecker_product(
                companions[left], companions[right]
            )))
    r_union = squarefree_lcm(r_components, budget, "explicit r lcm")
    all_union = squarefree_lcm([r_union] + diagonal_only, budget, "explicit r_all lcm")
    return r_union.degree(), all_union.degree()


def exact_core_at(t: Fraction, budget: CpuBudget) -> dict[str, object]:
    rational = build_rational_operator(t, physical_weight(t), budget)
    actual_scale, integer = clear_actual_lcm(rational)
    charpoly = sp.Poly(sp.Matrix(integer).charpoly(Z).as_expr(), Z, domain=sp.ZZ)
    content, raw_factors = sp.factor_list(charpoly)
    assert content == 1
    factors = [(sp.Poly(factor, Z, domain=sp.QQ).monic(), int(exponent))
               for factor, exponent in raw_factors]
    factors.sort(key=lambda item: (item[0].degree(), tuple(item[0].all_coeffs())))
    reconstruction = sp.Poly(1, Z, domain=sp.QQ)
    for factor, exponent in factors:
        assert factor.is_irreducible
        reconstruction *= factor**exponent
    assert reconstruction == sp.Poly(charpoly, Z, domain=sp.QQ)
    r_value, r_all_value = explicit_product_counts(factors, budget)
    return {
        "actual_lcm_scale": actual_scale,
        "factor_pattern": [(factor.degree(), exponent) for factor, exponent in factors],
        "r": r_value,
        "r_all": r_all_value,
        "pair_gcd_degree": 496 - r_value,
        "symmetric_gcd_degree": 528 - r_all_value,
    }


# ------------------------------------------------------- independent modular lower bounds

def modular_matrix(
    n: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    w: Fraction,
    p: int,
    budget: CpuBudget,
) -> np.ndarray:
    dimension = 1 << n
    p_matrix = np.zeros((dimension, dimension), dtype=np.int64)
    a, b = t.numerator, t.denominator
    for left in range(dimension):
        if left % 32 == 0:
            budget.check(f"modular P row {left}/{dimension}")
        for right in range(dimension):
            difference = left ^ right
            value = 1
            for site in range(n):
                value = value * (a if (difference >> (n - 1 - site)) & 1 else b) % p
            p_matrix[left, right] = value

    c, d = w.numerator, w.denominator
    diagonal = np.zeros(dimension, dtype=np.int64)
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        value = 1
        for left, right in bonds:
            value = value * (c if spins[left] == spins[right] else d) % p
        diagonal[state] = value
    return modmul((p_matrix * diagonal[np.newaxis, :]) % p, p_matrix, p)


def modmul(left: np.ndarray, right: np.ndarray, p: int) -> np.ndarray:
    assert left.shape[1] * (p - 1) ** 2 < 2**53
    product = np.asarray(left, dtype=np.float64) @ np.asarray(right, dtype=np.float64)
    return np.mod(product, float(p)).astype(np.int64)


def moddot(left: np.ndarray, right: np.ndarray, p: int) -> int:
    assert left.size == right.size
    assert left.size * (p - 1) ** 2 < 2**63 - 1
    return int(np.dot(left, right)) % p


def trim(poly: np.ndarray) -> np.ndarray:
    for index in range(min(poly.size, 64)):
        if poly[index]:
            return poly if index == 0 else poly[index:]
    nonzero = np.flatnonzero(poly)
    return np.zeros(1, dtype=np.int64) if nonzero.size == 0 else poly[nonzero[0]:]


def is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def monic(poly: np.ndarray, p: int) -> np.ndarray:
    result = trim(np.asarray(poly, dtype=np.int64) % p)
    if is_zero(result):
        return result
    lead = int(result[0])
    return result if lead == 1 else (result * pow(lead, -1, p)) % p


def remainder(
    dividend: np.ndarray,
    divisor: np.ndarray,
    p: int,
    budget: CpuBudget,
    stage: str,
) -> np.ndarray:
    assert int(divisor[0]) == 1
    work = (np.asarray(dividend, dtype=np.int64) % p).copy()
    divisor_degree = divisor.size - 1
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return trim(work)
    steps = dividend_degree - divisor_degree + 1
    for offset in range(steps):
        if offset % 4096 == 0:
            budget.check(f"{stage} remainder {offset}/{steps}")
        lead = int(work[offset])
        if lead:
            work[offset:offset + divisor_degree + 1] = (
                work[offset:offset + divisor_degree + 1] - lead * divisor
            ) % p
    return trim(work[steps:])


def poly_gcd(
    left: np.ndarray,
    right: np.ndarray,
    p: int,
    budget: CpuBudget,
    stage: str,
) -> np.ndarray:
    left, right = monic(left, p), monic(right, p)
    iteration = 0
    while not is_zero(right):
        budget.check(f"{stage} Euclid iteration {iteration}")
        left, right = right, monic(remainder(left, right, p, budget, stage), p)
        iteration += 1
    return left


def derivative(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    return np.ascontiguousarray(
        ((np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p)[::-1]
    )


def newton(
    power_sums: np.ndarray,
    degree: int,
    p: int,
    budget: CpuBudget,
    stage: str,
) -> np.ndarray:
    chunk = (2**53) // ((p - 1) ** 2)
    assert chunk >= 1
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    sums = np.ascontiguousarray(power_sums[:degree + 1], dtype=np.float64)
    window = np.zeros(degree + 1, dtype=np.float64)
    window[degree] = 1.0
    for order in range(1, degree + 1):
        if order % 256 == 0:
            budget.check(f"{stage} Newton {order}/{degree}")
        total = 0
        base = degree - order
        low = 1
        while low <= order:
            high = min(order, low + chunk - 1)
            total += int(np.dot(
                sums[low:high + 1],
                window[base + low:base + high + 1],
            )) % p
            low = high + 1
        value = (-(total % p) * pow(order, -1, p)) % p
        coefficients[order] = value
        window[degree - order] = float(value)
    return coefficients


def coefficient_digest_ascending(poly_descending: np.ndarray) -> str:
    payload = ",".join(str(int(value)) for value in poly_descending[::-1])
    return hashlib.sha256(payload.encode()).hexdigest()


def modular_pair_counts(
    n: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    w: Fraction,
    budget: CpuBudget,
    want_all: bool,
) -> dict[str, object]:
    matrix = modular_matrix(n, bonds, t, w, PRIME, budget)
    dimension = matrix.shape[0]
    pair_slots = comb(dimension, 2)
    all_slots = comb(dimension + 1, 2)
    needed = 2 * (all_slots if want_all else pair_slots)

    traces = np.zeros(needed + 1, dtype=np.int64)
    power = np.eye(dimension, dtype=np.int64)
    for order in range(1, dimension + 1):
        if order % 16 == 0:
            budget.check(f"modular matrix traces {order}/{dimension}")
        power = modmul(power, matrix, PRIME)
        traces[order] = int(power.trace()) % PRIME
    characteristic_poly = newton(
        traces[:dimension + 1], dimension, PRIME, budget, "base characteristic"
    )
    tail = np.ascontiguousarray(characteristic_poly[1:])
    for order in range(dimension + 1, needed + 1):
        if order % 2048 == 0:
            budget.check(f"trace recurrence {order}/{needed}")
        window = np.ascontiguousarray(traces[order - dimension:order][::-1])
        traces[order] = (-moddot(tail, window, PRIME)) % PRIME

    inverse_two = pow(2, -1, PRIME)

    def one_pair_polynomial(slots: int, sign: int, label: str) -> tuple[np.ndarray, int]:
        sums = np.zeros(slots + 1, dtype=np.int64)
        for order in range(1, slots + 1):
            if order % 2048 == 0:
                budget.check(f"{label} power sums {order}/{slots}")
            value = (int(traces[order]) * int(traces[order]) + sign * int(traces[2 * order])) % PRIME
            sums[order] = value * inverse_two % PRIME
        polynomial = newton(sums, slots, PRIME, budget, label)
        common = poly_gcd(
            polynomial,
            derivative(polynomial, PRIME),
            PRIME,
            budget,
            label,
        )
        return polynomial, slots - (common.size - 1)

    pair_poly, r_value = one_pair_polynomial(pair_slots, -1, "pair polynomial")
    result: dict[str, object] = {
        "r_modular_lower_bound": r_value,
        "pair_gcd_degree": pair_slots - r_value,
        "pair_poly_sha256": coefficient_digest_ascending(pair_poly),
    }
    if want_all:
        _all_poly, r_all_value = one_pair_polynomial(all_slots, +1, "symmetric polynomial")
        result["r_all_modular_lower_bound"] = r_all_value
        result["symmetric_gcd_degree"] = all_slots - r_all_value
    return result


def chain_bonds(m: int) -> list[tuple[int, int]]:
    return [(5 + index, 6 + index) for index in range(m - 1)]


def tensor_upper(m: int, r_core: int, r_all_core: int) -> int:
    return (3**m - 2**m) * r_all_core + 2**m * r_core


def gaussian_ceiling(n: int) -> int:
    return 3**n - 2**n

# ------------------------------------------------------------------------------- checks

def main() -> None:
    budget = CpuBudget()
    artifact = json.loads(ARTIFACT.read_text())
    data = artifact["data"]

    independent_core: dict[str, dict[str, object]] = {}
    expected_pattern = [(1, 2), (1, 2), (3, 1), (3, 1), (11, 1), (11, 1)]
    for t in NAMED_T:
        label = text_fraction(t)
        row = exact_core_at(t, budget)
        independent_core[label] = row
        ok(
            f"clean-room exact-Q core at t={label}",
            row["factor_pattern"] == expected_pattern
            and row["r"] == 417
            and row["r_all"] == 445
            and row["pair_gcd_degree"] == 79
            and row["symmetric_gcd_degree"] == 83,
            f"r={row['r']}, r_all={row['r_all']}, gcd degrees {row['pair_gcd_degree']}/{row['symmetric_gcd_degree']}",
        )

    artifact_points = {
        row["t"]: row for row in data["exact_characteristic_zero_points"]
    }
    ok(
        "artifact covers exactly all four named core rows",
        set(artifact_points) == {text_fraction(t) for t in NAMED_T},
        ", ".join(sorted(artifact_points)),
    )
    for t in NAMED_T:
        label = text_fraction(t)
        stored = artifact_points[label]["pair_product_certificate"]
        rebuilt = independent_core[label]
        ok(
            f"stored exact row agrees with clean-room row at t={label}",
            stored["r"] == rebuilt["r"]
            and stored["r_all"] == rebuilt["r_all"]
            and stored["r_two_sided_sandwich"] == {"lower": 417, "upper": 417}
            and stored["r_all_two_sided_sandwich"] == {"lower": 445, "upper": 445},
        )

    modular_m2: dict[str, dict[str, object]] = {}
    for t in NAMED_T:
        label = text_fraction(t)
        modular = modular_pair_counts(
            7,
            list(T_BONDS) + chain_bonds(2),
            t,
            physical_weight(t),
            budget,
            want_all=True,
        )
        modular_m2[label] = modular
        upper = tensor_upper(2, int(independent_core[label]["r"]), int(independent_core[label]["r_all"]))
        all_upper = 3**2 * int(independent_core[label]["r_all"])
        ok(
            f"independent U2 squeeze m=2 at t={label}",
            modular["r_modular_lower_bound"] == upper == 3893
            and modular["r_all_modular_lower_bound"] == all_upper == 4005,
            f"r >= {modular['r_modular_lower_bound']} and r <= {upper}; r_all >= {modular['r_all_modular_lower_bound']} and r_all <= {all_upper}",
        )

    t = Fraction(1, 3)
    modular_m4 = modular_pair_counts(
        9,
        list(T_BONDS) + chain_bonds(4),
        t,
        physical_weight(t),
        budget,
        want_all=False,
    )
    upper_m4 = tensor_upper(4, 417, 445)
    ok(
        "independent U2 squeeze m=4 at t=1/3",
        modular_m4["r_modular_lower_bound"] == upper_m4 == 35597
        and modular_m4["pair_gcd_degree"] == 95219,
        f"r >= {modular_m4['r_modular_lower_bound']} and r <= {upper_m4}; gcd degree {modular_m4['pair_gcd_degree']}",
    )

    formula_ok = True
    for m, exact_r in ((2, 3893), (4, 35597)):
        formula_ok &= (
            exact_r - gaussian_ceiling(5 + m)
            == 202 * 3**m + 4 * 2**m
        )
    ok(
        "exact restored excess formula at m=2,4",
        formula_ok,
        "excesses 1834 and 16426 equal 202*3^m+4*2^m",
    )

    theorem_rows = data["restored_U2_named_points"]
    theorem_keys = {(int(row["m"]), str(row["t"])) for row in theorem_rows}
    expected_keys = {(2, text_fraction(t)) for t in NAMED_T} | {(4, "1/3")}
    ok(
        "every universal named-point U2 row is present and independently covered",
        theorem_keys == expected_keys and len(theorem_rows) == len(expected_keys),
        str(sorted(theorem_keys)),
    )
    theorem_rows_ok = True
    for row in theorem_rows:
        key = (int(row["m"]), str(row["t"]))
        if key[0] == 2:
            lower = int(modular_m2[key[1]]["r_modular_lower_bound"])
            upper = tensor_upper(2, 417, 445)
        else:
            lower = int(modular_m4["r_modular_lower_bound"])
            upper = upper_m4
        theorem_rows_ok &= (
            row["two_sided_sandwich"] == {"lower": lower, "upper": upper}
            and row["exact_r"] == lower == upper
            and row["injectivity_INJ_m_certified"] is True
        )
    ok("all stored U2 sandwiches and INJ conclusions replay", theorem_rows_ok)

    p4_stored = next(row for row in theorem_rows if int(row["m"]) == 4)
    ok(
        "independent m=4 pair-polynomial digest matches the cited raw certificate",
        modular_m4["pair_poly_sha256"] == p4_stored["pair_polynomial_sha256_mod_prime"],
        str(modular_m4["pair_poly_sha256"]),
    )

    generic_rows = {int(row["m"]): row for row in data["F2_status"]["unconditional_generic_rows"]}
    generic_ok = set(generic_rows) == {2, 4}
    for m in (2, 4):
        generic_ok &= (
            generic_rows[m]["generic_lower_bound"]
            == gaussian_ceiling(5 + m) + 202 * 3**m + 4 * 2**m
            and generic_rows[m]["all_F2_hypotheses_discharged_for_this_m"] is True
        )
    ok("strengthened F2 rows follow from the two certified witnesses", generic_ok)

    ok(
        "producer recorded only passing checks",
        bool(data["checks"]) and all(bool(row["passed"]) for row in data["checks"]),
        f"{len(data['checks'])} producer checks",
    )
    ok(
        "verifier stayed below declared RSS cap",
        max_rss_bytes() < RSS_LIMIT_BYTES,
        f"peak RSS {max_rss_bytes()} < {RSS_LIMIT_BYTES}",
    )

    if FAILED:
        raise SystemExit(f"FAIL test_rank_char0: {len(FAILED)} failed of {len(PASSED) + len(FAILED)}")
    print(f"PASS test_rank_char0: {len(PASSED)} checks", flush=True)


if __name__ == "__main__":
    main()
