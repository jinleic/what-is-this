"""Exact characteristic-zero spectrum factorization for the five-site tree T.

This module is the first stage of the rank_char0 certificate.  It constructs an
integer scalar multiple of R_T(t,w) exactly, takes its characteristic polynomial
over Z, and factors it over Q.  No modular arithmetic or floating point enters.

Run from the repository root with:
    .venv/bin/python experiments/e166_rank_char0_spectrum.py
"""

from __future__ import annotations

import hashlib
import platform
import resource
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
T_BONDS: tuple[tuple[int, int], ...] = ((0, 1), (0, 2), (0, 3), (3, 4))
NAMED_T: tuple[Fraction, ...] = (
    Fraction(1, 3),
    Fraction(1, 4),
    Fraction(2, 5),
    Fraction(1, 5),
)
Z = sp.Symbol("z")
CPU_BUDGET_SECONDS = 1_200.0
RSS_LIMIT_BYTES = 4_000_000_000


class NonDecisiveBudget(RuntimeError):
    """A declared process-CPU or RSS budget expired without a mathematical verdict."""


@dataclass
class CpuBudget:
    limit_seconds: float = CPU_BUDGET_SECONDS
    started: float = 0.0

    def __post_init__(self) -> None:
        if not self.started:
            self.started = time.process_time()

    def check(self, stage: str) -> None:
        used = time.process_time() - self.started
        if used > self.limit_seconds:
            raise NonDecisiveBudget(
                f"NON-DECISIVE: process-CPU budget {self.limit_seconds:.1f}s "
                f"expired at {stage} after {used:.3f}s"
            )
        rss = max_rss_bytes()
        if rss > RSS_LIMIT_BYTES:
            raise NonDecisiveBudget(
                f"NON-DECISIVE: RSS budget {RSS_LIMIT_BYTES} bytes exceeded "
                f"at {stage} (ru_maxrss={rss})"
            )


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def physical_weight(t: Fraction) -> Fraction:
    if not t:
        raise ValueError("the physical curve excludes t=0")
    return (1 + t * t) / (2 * t)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _hash_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def matrix_digest(rows: Iterable[Iterable[int]]) -> str:
    payload = ";".join(",".join(str(int(value)) for value in row) for row in rows)
    return _hash_text(payload)


def polynomial_digest(poly: sp.Poly) -> str:
    """Hash monic rational coefficients without decimal integer conversion."""

    monic = sp.Poly(poly, Z, domain=sp.QQ).monic()
    digest = hashlib.sha256()
    for coefficient in monic.all_coeffs():
        numerator, denominator = (int(value) for value in sp.fraction(coefficient))
        for integer in (numerator, denominator):
            magnitude = abs(integer)
            raw = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
            digest.update(b"-" if integer < 0 else b"+")
            digest.update(len(raw).to_bytes(8, "big"))
            digest.update(raw)
    return digest.hexdigest()


def polynomial_coefficients(poly: sp.Poly) -> list[str]:
    monic = sp.Poly(poly, Z, domain=sp.QQ).monic()
    result: list[str] = []
    for coefficient in monic.all_coeffs():
        numerator, denominator = sp.fraction(coefficient)
        if denominator == 1:
            result.append(str(int(numerator)))
        else:
            result.append(f"{int(numerator)}/{int(denominator)}")
    return result


def build_integer_operator(
    t: Fraction,
    w: Fraction,
    budget: CpuBudget | None = None,
) -> tuple[int, list[list[int]]]:
    """Return M = den(t)^10 den(w)^4 R_T(t,w) as an exact integer matrix.

    In the computational spin basis, P(t)_{sigma,tau}=t^Hamming(sigma,tau) and
    D_T(w)_{sigma,sigma}=w^(number of aligned T-edges).  The returned matrix is
    P_int D_int P_int, where P_int has local factor [[b,a],[a,b]] for t=a/b,
    and an aligned/non-aligned edge contributes c/d to D before clearing d^4.
    """

    if budget is None:
        budget = CpuBudget()
    a, b = t.numerator, t.denominator
    c, d = w.numerator, w.denominator
    dimension = 1 << 5
    diagonal: list[int] = []
    for state in range(dimension):
        if state % 8 == 0:
            budget.check("integer diagonal")
        spins = [1 - 2 * ((state >> (4 - site)) & 1) for site in range(5)]
        value = 1
        for left, right in T_BONDS:
            value *= c if spins[left] == spins[right] else d
        diagonal.append(value)

    matrix = [[0] * dimension for _ in range(dimension)]
    for state, value in enumerate(diagonal):
        matrix[state][state] = value

    def apply_local_factors(rows: list[list[int]]) -> None:
        for site in range(5):
            bit = 1 << (4 - site)
            budget.check(f"Kronecker factor {site}")
            for state in range(dimension):
                if state & bit:
                    continue
                row_zero = rows[state]
                row_one = rows[state | bit]
                for column in range(dimension):
                    zero, one = row_zero[column], row_one[column]
                    row_zero[column] = b * zero + a * one
                    row_one[column] = a * zero + b * one

    apply_local_factors(matrix)
    matrix = [list(column) for column in zip(*matrix)]
    apply_local_factors(matrix)
    assert all(matrix[i][j] == matrix[j][i] for i in range(dimension) for j in range(dimension))
    scale = b**10 * d**4
    return scale, matrix


def factor_characteristic_polynomial(
    matrix: list[list[int]],
    budget: CpuBudget,
) -> tuple[sp.Poly, list[tuple[sp.Poly, int]], dict[str, object]]:
    budget.check("before characteristic polynomial")
    characteristic = sp.Poly(sp.Matrix(matrix).charpoly(Z).as_expr(), Z, domain=sp.ZZ)
    budget.check("after characteristic polynomial")
    content, raw_factors = sp.factor_list(characteristic)
    assert content == 1
    factors = [(sp.Poly(factor, Z, domain=sp.ZZ).monic(), int(exponent))
               for factor, exponent in raw_factors]
    factors.sort(key=lambda item: (item[0].degree(), tuple(item[0].all_coeffs())))

    reconstruction = sp.Poly(1, Z, domain=sp.QQ)
    factor_rows: list[dict[str, object]] = []
    for index, (factor, exponent) in enumerate(factors):
        budget.check(f"factor audit {index}")
        reconstruction *= factor**exponent
        factor_rows.append({
            "index": index,
            "degree": factor.degree(),
            "exponent": exponent,
            "coefficients_desc": polynomial_coefficients(factor),
            "sha256": polynomial_digest(factor),
            "irreducible_over_Q": bool(factor.is_irreducible),
        })
    reconstruction = sp.Poly(reconstruction, Z, domain=sp.QQ)
    assert reconstruction == sp.Poly(characteristic, Z, domain=sp.QQ)
    assert all(row["irreducible_over_Q"] for row in factor_rows)

    record: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "degree": characteristic.degree(),
        "sha256": polynomial_digest(characteristic),
        "factorization_reconstructs": True,
        "factor_degree_exponent_pattern": [
            {"degree": factor.degree(), "exponent": exponent} for factor, exponent in factors
        ],
        "distinct_eigenvalue_values": sum(factor.degree() for factor, _ in factors),
        "slot_multiplicity_excess": sum(factor.degree() * (exponent - 1)
                                         for factor, exponent in factors),
        "factors": factor_rows,
    }
    return characteristic, factors, record


def exact_spectrum_record(
    t: Fraction,
    budget: CpuBudget,
) -> tuple[dict[str, object], list[tuple[sp.Poly, int]]]:
    w = physical_weight(t)
    scale, matrix = build_integer_operator(t, w, budget)
    characteristic, factors, characteristic_record = factor_characteristic_polynomial(matrix, budget)
    record: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "t": fraction_text(t),
        "w": fraction_text(w),
        "graph_vertices": [0, 1, 2, 3, 4],
        "graph_edges": [list(edge) for edge in T_BONDS],
        "operator_dimension": len(matrix),
        "integer_scale_M_over_R": scale,
        "integer_matrix_sha256": matrix_digest(matrix),
        "matrix_symmetric": all(matrix[i][j] == matrix[j][i]
                                for i in range(len(matrix)) for j in range(len(matrix))),
        "characteristic_polynomial": characteristic_record,
    }
    assert characteristic.degree() == 32
    return record, factors


def main() -> None:
    budget = CpuBudget()
    record, factors = exact_spectrum_record(Fraction(1, 3), budget)
    pattern = [(factor.degree(), exponent) for factor, exponent in factors]
    checks = [
        record["operator_dimension"] == 32,
        record["matrix_symmetric"] is True,
        pattern == [(1, 2), (1, 2), (3, 1), (3, 1), (11, 1), (11, 1)],
        record["characteristic_polynomial"]["factorization_reconstructs"] is True,
        record["characteristic_polynomial"]["distinct_eigenvalue_values"] == 30,
    ]
    if not all(checks):
        raise SystemExit("FAIL e166_rank_char0_spectrum")
    print(
        "PASS e166_rank_char0_spectrum: exact-Q factor pattern "
        "(1^2,1^2,3,3,11,11)",
        flush=True,
    )


if __name__ == "__main__":
    main()
