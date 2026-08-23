"""[COMPUTATION] Exact characteristic-polynomial factors for small Ising layers.

The canonical rational representative is ``e38_gaussianity_certificate.build_R``.
For each graph at ``t=1/3`` this module clears the actual entry denominators,
divides the resulting integer matrix by its entry content, and factors its monic
characteristic polynomial over ``Z``.  The two degree-six factors selected for
Galois certification are chosen deterministically, not by numerical roots.

Run from the repository root with:
    .venv/bin/python experiments/e203_galois_charpoly.py
"""

from __future__ import annotations

import hashlib
import math
import os
import platform
import resource
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

# Keep accidental numerical-library imports single-threaded on the shared host.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from e38_gaussianity_certificate import build_R  # noqa: E402
from ising.transfer_matrix import layer_bonds  # noqa: E402

X = sp.Symbol("x")
T_SPECIALIZATION = Fraction(1, 3)
Q_SPECIALIZATION = Fraction(5, 3)
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2_000_000_000
TARGET_CASES = ("ising_2d_open_width6", "ising_3d_open_2x3")


class NonDecisiveBudget(RuntimeError):
    """[UNRESOLVED] The declared local resource budget expired before a verdict."""


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


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
                f"[UNRESOLVED] process-CPU budget {self.limit_seconds:.1f}s "
                f"expired at {stage} after {used:.3f}s"
            )
        rss = max_rss_bytes()
        if rss > RSS_LIMIT_BYTES:
            raise NonDecisiveBudget(
                f"[UNRESOLVED] RSS limit {RSS_LIMIT_BYTES} bytes exceeded at "
                f"{stage}; ru_maxrss={rss}"
            )


@dataclass(frozen=True)
class CaseSpec:
    key: str
    spatial_dimension: int
    cross_section: str
    n_sites: int
    bonds: tuple[tuple[int, int], ...]
    boundary: str
    interpretation: str


def _bonds(cross: tuple[int, ...], periodic: tuple[bool, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((int(left), int(right)) for left, right in layer_bonds(cross, periodic))


CASES: tuple[CaseSpec, ...] = (
    CaseSpec(
        key="ising_1d_single_spin",
        spatial_dimension=1,
        cross_section="single site",
        n_sites=1,
        bonds=(),
        boundary="no transverse bonds",
        interpretation="[COMPUTATION] one-dimensional Ising transfer control",
    ),
    CaseSpec(
        key="ising_2d_periodic_width4",
        spatial_dimension=2,
        cross_section="C4",
        n_sites=4,
        bonds=_bonds((4,), (True,)),
        boundary="periodic transverse chain",
        interpretation="[EXTERNAL] standard finite-width integrable 2D Ising control",
    ),
    CaseSpec(
        key="ising_3d_open_2x2",
        spatial_dimension=3,
        cross_section="2x2",
        n_sites=4,
        bonds=_bonds((2, 2), (False, False)),
        boundary="open square layer",
        interpretation="[COMPUTATION] smallest open square 3D layer; its graph is C4",
    ),
    CaseSpec(
        key="ising_2d_open_width6",
        spatial_dimension=2,
        cross_section="P6",
        n_sites=6,
        bonds=_bonds((6,), (False,)),
        boundary="open transverse chain",
        interpretation="[EXTERNAL] standard open-strip free-fermion 2D Ising control",
    ),
    CaseSpec(
        key="ising_3d_open_2x3",
        spatial_dimension=3,
        cross_section="2x3",
        n_sites=6,
        bonds=_bonds((2, 3), (False, False)),
        boundary="open rectangular layer",
        interpretation="[COMPUTATION] first open rectangular layer not isomorphic to a cycle",
    ),
)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def integer_matrix_digest(rows: Iterable[Iterable[int]]) -> str:
    payload = ";".join(",".join(str(int(value)) for value in row) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def polynomial_digest(poly: sp.Poly) -> str:
    exact = sp.Poly(poly, X, domain=sp.QQ).monic()
    payload = "\n".join(str(coefficient) for coefficient in exact.all_coeffs())
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def coefficients_as_text(poly: sp.Poly) -> list[str]:
    exact = sp.Poly(poly, X, domain=sp.QQ).monic()
    return [str(coefficient) for coefficient in exact.all_coeffs()]


def primitive_integer_operator(
    spec: CaseSpec,
    budget: CpuBudget,
) -> tuple[list[list[int]], dict[str, object]]:
    """[COMPUTATION] Return the primitive integer scalar multiple of e38's R."""

    budget.check(f"before build_R for {spec.key}")
    rational, q_value, parity = build_R(spec.n_sites, spec.bonds, T_SPECIALIZATION)
    if q_value != Q_SPECIALIZATION:
        raise ArithmeticError("unexpected physical-curve weight")

    denominator_lcm = 1
    for row in rational:
        for value in row:
            denominator_lcm = math.lcm(denominator_lcm, value.denominator)
    integer = [[int(value * denominator_lcm) for value in row] for row in rational]
    if not all(
        value * denominator_lcm == integer[row_index][column_index]
        for row_index, row in enumerate(rational)
        for column_index, value in enumerate(row)
    ):
        raise ArithmeticError("denominator clearing was not exact")

    entry_content = 0
    for row in integer:
        for value in row:
            entry_content = math.gcd(entry_content, abs(value))
    if entry_content == 0:
        raise ArithmeticError("zero transfer representative")
    primitive = [[value // entry_content for value in row] for row in integer]
    if math.gcd(*(abs(value) for row in primitive for value in row)) != 1:
        raise ArithmeticError("primitive matrix still has nontrivial entry content")
    dimension = 1 << spec.n_sites
    symmetric = all(
        primitive[left][right] == primitive[right][left]
        for left in range(dimension)
        for right in range(dimension)
    )
    budget.check(f"after primitive matrix for {spec.key}")
    metadata: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "q": fraction_text(q_value),
        "bond_parity": int(parity),
        "entry_denominator_lcm": denominator_lcm,
        "cleared_matrix_entry_content": entry_content,
        "primitive_scale_B_over_R": fraction_text(Fraction(denominator_lcm, entry_content)),
        "primitive_matrix_sha256": integer_matrix_digest(primitive),
        "primitive_matrix_symmetric": symmetric,
    }
    return primitive, metadata


def factor_characteristic_polynomial(
    matrix: list[list[int]],
    budget: CpuBudget,
    stage: str,
) -> tuple[list[tuple[sp.Poly, int]], dict[str, object]]:
    """[COMPUTATION] Factor ``det(xI-B)`` exactly over characteristic zero."""

    budget.check(f"before charpoly {stage}")
    characteristic = sp.Poly(sp.Matrix(matrix).charpoly(X).as_expr(), X, domain=sp.ZZ)
    budget.check(f"after charpoly {stage}")
    content, raw_factors = sp.factor_list(characteristic)
    if content != 1:
        raise ArithmeticError(f"nonunit characteristic content at {stage}: {content}")

    factors: list[tuple[sp.Poly, int]] = []
    for raw_factor, exponent in raw_factors:
        factor = sp.Poly(raw_factor, X, domain=sp.ZZ)
        if factor.LC() != 1:
            raise ArithmeticError("nonmonic factor of a monic integer characteristic polynomial")
        factors.append((factor, int(exponent)))
    factors.sort(
        key=lambda item: (
            item[0].degree(),
            tuple(int(coefficient) for coefficient in item[0].all_coeffs()),
            item[1],
        )
    )

    reconstruction = sp.Poly(1, X, domain=sp.ZZ)
    factor_rows: list[dict[str, object]] = []
    for index, (factor, exponent) in enumerate(factors):
        budget.check(f"factor audit {stage} index {index}")
        reconstruction *= factor**exponent
        irreducible = bool(factor.is_irreducible)
        factor_rows.append(
            {
                "claim_tag": "[COMPUTATION]",
                "index": index,
                "degree": factor.degree(),
                "exponent": exponent,
                "coefficients_desc": coefficients_as_text(factor),
                "sha256": polynomial_digest(factor),
                "irreducible_over_Q": irreducible,
            }
        )
    reconstruction_ok = reconstruction == characteristic
    all_irreducible = all(bool(row["irreducible_over_Q"]) for row in factor_rows)
    if not reconstruction_ok or not all_irreducible:
        raise ArithmeticError(f"exact factor audit failed at {stage}")

    record: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "degree": characteristic.degree(),
        "sha256": polynomial_digest(characteristic),
        "factorization_reconstructs": reconstruction_ok,
        "all_factors_irreducible_over_Q": all_irreducible,
        "factor_degree_exponent_pattern": [
            {"degree": factor.degree(), "exponent": exponent}
            for factor, exponent in factors
        ],
        "factors": factor_rows,
    }
    return factors, record


def _valuation(value: int, prime: int) -> int:
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def maximal_integral_root_scale(factor: sp.Poly) -> int:
    """[LEMMA] Largest integer s with ``s^k`` dividing coefficient k."""

    coefficients = [int(value) for value in factor.all_coeffs()]
    candidates: dict[int, int] = {}
    for index, coefficient in enumerate(coefficients[1:], start=1):
        if coefficient:
            candidates = {
                int(prime): int(exponent) // index
                for prime, exponent in sp.factorint(abs(coefficient)).items()
            }
            break
    if not candidates:
        return 1
    for prime in list(candidates):
        bound = candidates[prime]
        for index, coefficient in enumerate(coefficients[1:], start=1):
            if coefficient:
                bound = min(bound, _valuation(abs(coefficient), prime) // index)
        candidates[prime] = bound
    scale = math.prod(prime**exponent for prime, exponent in candidates.items() if exponent > 0)
    for index, coefficient in enumerate(coefficients):
        if coefficient % scale**index:
            raise ArithmeticError("computed root scale is not integral")
    return scale


def selected_sextic_certificate(
    factors: list[tuple[sp.Poly, int]],
    budget: CpuBudget,
    case_key: str,
) -> dict[str, object]:
    """[COMPUTATION] Select and normalize the lexicographically first sextic."""

    sextics = [(factor, exponent) for factor, exponent in factors if factor.degree() == 6]
    if not sextics:
        raise ArithmeticError(f"no sextic factor for {case_key}")
    sextics.sort(key=lambda item: tuple(int(value) for value in item[0].all_coeffs()))
    factor, exponent = sextics[0]
    scale = maximal_integral_root_scale(factor)
    raw_coefficients = [int(value) for value in factor.all_coeffs()]
    normalized_coefficients = [
        coefficient // scale**index for index, coefficient in enumerate(raw_coefficients)
    ]
    normalized = sp.Poly.from_list(normalized_coefficients, gens=X, domain=sp.ZZ)
    substituted = sp.Poly(factor.as_expr().subs(X, scale * X) / scale**6, X, domain=sp.QQ)
    normalization_ok = substituted == sp.Poly(normalized, X, domain=sp.QQ)
    if not normalization_ok:
        raise ArithmeticError("rational root normalization identity failed")

    budget.check(f"before discriminant {case_key}")
    discriminant = int(sp.discriminant(normalized, X))
    discriminant_factors = {
        int(prime): int(power) for prime, power in sp.factorint(abs(discriminant)).items()
    }
    factorization_reconstructs = (
        math.prod(prime**power for prime, power in discriminant_factors.items())
        == abs(discriminant)
    )
    factors_prime = all(bool(sp.isprime(prime)) for prime in discriminant_factors)
    discriminant_square = math.isqrt(abs(discriminant)) ** 2 == abs(discriminant)
    odd_primes = sorted(
        prime for prime, power in discriminant_factors.items() if power % 2
    )
    budget.check(f"after discriminant {case_key}")

    return {
        "claim_tag": "[COMPUTATION]",
        "selection_rule": (
            "[COMPUTATION] lexicographically smallest descending coefficient tuple "
            "among irreducible degree-six factors of the primitive matrix charpoly"
        ),
        "raw_factor_exponent_in_charpoly": exponent,
        "raw_factor_coefficients_desc": [str(value) for value in raw_coefficients],
        "raw_factor_sha256": polynomial_digest(factor),
        "rational_root_scale_s": scale,
        "normalization_identity": (
            "[LEMMA] g(y)=f_B(s*y)/s^6; rational root scaling leaves the "
            "splitting field over Q unchanged"
        ),
        "normalization_identity_checked": normalization_ok,
        "normalized_factor_coefficients_desc": [
            str(value) for value in normalized_coefficients
        ],
        "normalized_factor_sha256": polynomial_digest(normalized),
        "normalized_factor_irreducible_over_Q_from_charpoly_factorization": bool(
            normalized.is_irreducible
        ),
        "discriminant": str(discriminant),
        "discriminant_sign": 1 if discriminant > 0 else -1 if discriminant < 0 else 0,
        "discriminant_factorization": {
            str(prime): power for prime, power in sorted(discriminant_factors.items())
        },
        "discriminant_factorization_reconstructs": factorization_reconstructs,
        "discriminant_factors_all_prime": factors_prime,
        "discriminant_is_square": discriminant_square,
        "odd_exponent_discriminant_primes": odd_primes,
    }


def compute_case(spec: CaseSpec, budget: CpuBudget) -> tuple[dict[str, object], list[tuple[sp.Poly, int]]]:
    primitive, matrix_record = primitive_integer_operator(spec, budget)
    factors, characteristic_record = factor_characteristic_polynomial(
        primitive, budget, spec.key
    )
    record: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "key": spec.key,
        "spatial_dimension": spec.spatial_dimension,
        "cross_section": spec.cross_section,
        "n_sites": spec.n_sites,
        "operator_dimension": 1 << spec.n_sites,
        "bonds": [list(edge) for edge in spec.bonds],
        "boundary": spec.boundary,
        "interpretation": spec.interpretation,
        "t": fraction_text(T_SPECIALIZATION),
        **matrix_record,
        "characteristic_polynomial": characteristic_record,
    }
    if spec.key in TARGET_CASES:
        record["selected_sextic"] = selected_sextic_certificate(
            factors, budget, spec.key
        )
    return record, factors


def compute_all_cases(
    budget: CpuBudget | None = None,
) -> tuple[list[dict[str, object]], dict[str, list[tuple[sp.Poly, int]]]]:
    if budget is None:
        budget = CpuBudget()
    rows: list[dict[str, object]] = []
    factors_by_case: dict[str, list[tuple[sp.Poly, int]]] = {}
    for spec in CASES:
        budget.check(f"start case {spec.key}")
        row, factors = compute_case(spec, budget)
        rows.append(row)
        factors_by_case[spec.key] = factors
    return rows, factors_by_case


def main() -> None:
    budget = CpuBudget()
    rows, _factors = compute_all_cases(budget)
    by_key = {str(row["key"]): row for row in rows}
    expected_keys = {spec.key for spec in CASES}
    exact_audits = all(
        row["operator_dimension"] == row["characteristic_polynomial"]["degree"]
        and row["primitive_matrix_symmetric"] is True
        and row["characteristic_polynomial"]["factorization_reconstructs"] is True
        and row["characteristic_polynomial"]["all_factors_irreducible_over_Q"] is True
        for row in rows
    )
    square_pattern_equal = (
        by_key["ising_2d_periodic_width4"]["characteristic_polynomial"]["factors"]
        == by_key["ising_3d_open_2x2"]["characteristic_polynomial"]["factors"]
    )
    target_sextics = all(
        by_key[key]["selected_sextic"]["normalization_identity_checked"] is True
        and by_key[key]["selected_sextic"]["discriminant_is_square"] is False
        for key in TARGET_CASES
    )
    if set(by_key) != expected_keys or not exact_audits or not square_pattern_equal or not target_sextics:
        raise SystemExit("FAIL e203_galois_charpoly")
    print(
        "PASS e203_galois_charpoly: exact factors for 1D, 2D controls, "
        "and open 2x2/2x3 3D layers",
        flush=True,
    )


if __name__ == "__main__":
    main()
