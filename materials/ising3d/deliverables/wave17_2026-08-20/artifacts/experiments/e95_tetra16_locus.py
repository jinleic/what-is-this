"""Exact specialization-locus certificate for the cabled tetra16 RLLL system.

This experiment avoids a dense symbolic 128x128 determinant.  It evaluates the
four selected sector minors modulo enough exact primes, interpolates their
univariate determinant polynomials, and lifts the coefficients by CRT with a
recorded coefficient-height bound.  Two disjoint row selections give a
complementary-minor Bezout certificate for the remaining all-specialization
candidate locus.

Run from the repository root with::

    timeout 7200 .venv/bin/python experiments/e95_tetra16_locus.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import signal
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_add, gf_gcdex, gf_mul
from scipy.optimize import linear_sum_assignment

import e83_tetra16 as T16


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "integrability" / "tetra16.json"
OUTPUT = ROOT / "results" / "integrability" / "tetra16_locus.json"
Q = sp.Symbol("q")
SIZE = T16.SECTOR_SIZE
DENSE_SCALE_EXPONENT = 12 * SIZE
# Products in the modular eliminator never exceed prime**2 in one int64 lane.
PRIME_CEILING = 2_000_000_000
BEZOUT_PRIME_ATTEMPTS = 16
EXACT_GCD_WALL_SECONDS = int(os.environ.get("TETRA16_LOCUS_EXACT_GCD_WALL_SECONDS", "600"))
COMMON_FACTOR_WALL_SECONDS = int(
    os.environ.get("TETRA16_LOCUS_COMMON_FACTOR_WALL_SECONDS", "300")
)


class FactorizationTimeout(TimeoutError):
    """Raised only after a recorded bounded exact factorization search."""


@contextmanager
def _hard_timeout(seconds: int):
    previous_handler = signal.getsignal(signal.SIGALRM)

    def _raise_timeout(_signum, _frame):
        raise FactorizationTimeout(f"exact factorization exceeded {seconds} seconds")

    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


@dataclass(frozen=True)
class EntryData:
    """A nonzero Z[q] matrix entry, with coefficients in ascending degree."""

    row: int
    column: int
    coefficients: tuple[int, ...]

    @property
    def valuation(self) -> int:
        return next(index for index, coefficient in enumerate(self.coefficients) if coefficient)

    @property
    def degree(self) -> int:
        return len(self.coefficients) - 1

    @property
    def l1_norm(self) -> int:
        return sum(abs(coefficient) for coefficient in self.coefficients)


@dataclass(frozen=True)
class MinorData:
    """A fixed 128x128 sector minor and exact combinatorial degree data."""

    label: str
    sector: str
    rows: tuple[int, ...]
    entries: tuple[EntryData, ...]
    valuation_lower_bound: int
    determinant_degree_upper_bound: int
    shifted_degree_upper_bound: int
    height_bound: int
    min_assignment: tuple[int, ...]
    max_assignment: tuple[int, ...]


def _q_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _integer_poly_coefficients(value: sp.Expr) -> tuple[int, ...]:
    polynomial = sp.Poly(sp.expand(value), Q, domain=sp.ZZ)
    if polynomial.is_zero:
        return ()
    degree = int(polynomial.degree())
    coefficients = [0] * (degree + 1)
    for (power,), coefficient in polynomial.as_dict().items():
        coefficients[int(power)] = int(coefficient)
    return tuple(coefficients)


def _assignment_extremum(weights: np.ndarray, *, minimize: bool) -> tuple[int, tuple[int, ...]]:
    """Return an exact integer assignment extremum for a finite small cost matrix.

    The input costs are small exact integers (at most the entry degrees).  SciPy's
    assignment implementation is used only as a combinatorial optimizer; every
    selected edge and its integer sum are rechecked below.
    """

    chosen_rows, chosen_columns = linear_sum_assignment(weights)
    if len(chosen_rows) != SIZE or len(chosen_columns) != SIZE:
        raise ArithmeticError("degree-bound assignment did not cover every row and column")
    columns_by_row = [-1] * SIZE
    for row, column in zip(chosen_rows, chosen_columns, strict=True):
        columns_by_row[int(row)] = int(column)
    if any(column < 0 for column in columns_by_row):
        raise ArithmeticError("degree-bound assignment omitted a row")
    total = sum(int(weights[row, column]) for row, column in enumerate(columns_by_row))
    if minimize:
        return total, tuple(columns_by_row)
    return -total, tuple(columns_by_row)


def _minor_data(label: str, sector: str, rows: Sequence[int], forward: sp.Matrix, reverse: sp.Matrix) -> MinorData:
    """Build compact polynomial entries and tropical determinant bounds."""

    matrix = T16._symbolic_sector_minor(rows, forward, reverse, sector)
    if matrix.shape != (SIZE, SIZE):
        raise ArithmeticError("sector minor has the wrong shape")
    infinite = 10**6
    low = np.full((SIZE, SIZE), infinite, dtype=np.int64)
    negative_high = np.full((SIZE, SIZE), infinite, dtype=np.int64)
    entries: list[EntryData] = []
    row_l1_squares = [0] * SIZE
    for row in range(SIZE):
        for column in range(SIZE):
            value = matrix[row, column]
            coefficients = _integer_poly_coefficients(value)
            if not coefficients:
                continue
            entry = EntryData(row, column, coefficients)
            entries.append(entry)
            low[row, column] = entry.valuation
            negative_high[row, column] = -entry.degree
            row_l1_squares[row] += entry.l1_norm**2
    if not entries or any(value == 0 for value in row_l1_squares):
        raise ArithmeticError("selected minor has an empty support row")
    valuation, min_assignment = _assignment_extremum(low, minimize=True)
    degree, max_assignment = _assignment_extremum(negative_high, minimize=False)
    if any(low[row, column] >= infinite for row, column in enumerate(min_assignment)):
        raise ArithmeticError("minimum-degree matching uses a zero entry")
    if any(negative_high[row, column] >= infinite for row, column in enumerate(max_assignment)):
        raise ArithmeticError("maximum-degree matching uses a zero entry")
    if degree < valuation:
        raise ArithmeticError("tropical determinant bounds are inconsistent")
    # On |q|=1, every entry is bounded by its coefficient l1 norm.  Hadamard
    # then bounds |det M(q)| by prod_rows sqrt(sum_columns ||M_ij||_1^2).
    # Cauchy's coefficient inequality gives the same bound for every coefficient.
    squared_height_bound = math.prod(row_l1_squares)
    height_bound = math.isqrt(squared_height_bound)
    if height_bound * height_bound < squared_height_bound:
        height_bound += 1
    return MinorData(
        label=label,
        sector=sector,
        rows=tuple(int(row) for row in rows),
        entries=tuple(entries),
        valuation_lower_bound=valuation,
        determinant_degree_upper_bound=degree,
        shifted_degree_upper_bound=degree - valuation,
        height_bound=height_bound,
        min_assignment=min_assignment,
        max_assignment=max_assignment,
    )


def _alternate_rows(sector: str, forward: sp.Matrix, reverse: sp.Matrix, prime: int = 107) -> tuple[int, ...]:
    """Choose a complementary full-rank row set by reverse-order F_p elimination."""

    basis: list[list[int] | None] = [None] * SIZE
    selected: list[int] = []
    all_rows = list(T16._sector_row_indices(sector))
    for row_index in reversed(all_rows):
        output_state, input_state = divmod(row_index, T16.DIMENSION)
        vector = [
            int(value) % prime
            for value in T16.component_row_sector(
                output_state, input_state, forward, reverse, sector
            )
        ]
        for column, pivot in enumerate(basis):
            if pivot is None or vector[column] == 0:
                continue
            coefficient = vector[column]
            vector = [
                (value - coefficient * pivot_value) % prime
                for value, pivot_value in zip(vector, pivot, strict=True)
            ]
        pivot_column = next((column for column, value in enumerate(vector) if value), None)
        if pivot_column is None:
            continue
        inverse = pow(vector[pivot_column], -1, prime)
        basis[pivot_column] = [(value * inverse) % prime for value in vector]
        selected.append(row_index)
        if len(selected) == SIZE:
            break
    if len(selected) != SIZE:
        raise ArithmeticError(f"reverse row selection failed in {sector}")
    return tuple(selected)


def _modular_determinant(entries: Sequence[EntryData], node: int, prime: int) -> int:
    """Evaluate a sparse polynomial minor then take its exact determinant over F_p."""

    if not (prime < PRIME_CEILING and prime * prime < 2**63):
        raise ArithmeticError("chosen modular prime is unsafe for int64 elimination")
    x = node % prime
    powers = [1]
    maximum_degree = max(entry.degree for entry in entries)
    for _ in range(maximum_degree):
        powers.append((powers[-1] * x) % prime)
    matrix = np.zeros((SIZE, SIZE), dtype=np.int64)
    for entry in entries:
        value = sum(coefficient * powers[power] for power, coefficient in enumerate(entry.coefficients))
        matrix[entry.row, entry.column] = value % prime
    determinant = 1
    for column in range(SIZE):
        pivot = next((row for row in range(column, SIZE) if matrix[row, column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[[column, pivot]] = matrix[[pivot, column]]
            determinant = -determinant
        pivot_value = int(matrix[column, column])
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        multipliers = (matrix[column + 1 :, column] * inverse) % prime
        if len(multipliers):
            matrix[column + 1 :, column + 1 :] = (
                matrix[column + 1 :, column + 1 :]
                - multipliers[:, None] * matrix[column, column + 1 :]
            ) % prime
            matrix[column + 1 :, column] = 0
    return int(determinant % prime)


def _interpolate_mod_prime(minor: MinorData, prime: int) -> tuple[int, ...]:
    """Interpolate q^-v det(M(q)) in F_p[q] at 1,...,degree+1 exactly."""

    point_count = minor.shifted_degree_upper_bound + 1
    if prime <= point_count:
        raise ArithmeticError("interpolation prime does not separate the integer nodes")
    values = []
    for node in range(1, point_count + 1):
        determinant = _modular_determinant(minor.entries, node, prime)
        divisor = pow(node, minor.valuation_lower_bound, prime)
        values.append(determinant * pow(divisor, prime - 2, prime) % prime)
    # In-place Newton divided differences at equally spaced nodes 1,...,N.
    differences = values[:]
    for order in range(1, point_count):
        inverse = pow(order, prime - 2, prime)
        for index in range(point_count - 1, order - 1, -1):
            differences[index] = (differences[index] - differences[index - 1]) * inverse % prime
    # Convert Newton basis prod_{j<k}(q-(j+1)) to the monomial basis.
    coefficients = [0] * point_count
    basis = [1] + [0] * (point_count - 1)
    for order, coefficient in enumerate(differences):
        if coefficient:
            for power, basis_coefficient in enumerate(basis):
                if basis_coefficient:
                    coefficients[power] = (coefficients[power] + coefficient * basis_coefficient) % prime
        if order + 1 == point_count:
            continue
        next_basis = [0] * point_count
        node = order + 1
        for power in range(point_count - 1, -1, -1):
            value = basis[power]
            if value == 0:
                continue
            if power + 1 < point_count:
                next_basis[power + 1] = (next_basis[power + 1] + value) % prime
            next_basis[power] = (next_basis[power] - value * node) % prime
        basis = next_basis
    # Re-evaluate at all interpolation nodes; this catches any implementation error
    # before CRT turns it into a large false certificate.
    for node, value in enumerate(values, start=1):
        recovered = 0
        for coefficient in reversed(coefficients):
            recovered = (recovered * node + coefficient) % prime
        if recovered != value:
            raise ArithmeticError(f"modular interpolation self-check failed in {minor.label}")
    return tuple(int(value) for value in coefficients)


def _primes_for_height(height_bound: int) -> tuple[int, ...]:
    """Return enough distinct safe primes for unambiguous balanced CRT lifting."""

    modulus = 1
    candidate = PRIME_CEILING
    primes: list[int] = []
    target = 2 * height_bound
    while modulus <= target:
        prime = int(sp.prevprime(candidate))
        if prime <= 610:
            raise ArithmeticError("ran out of safe interpolation primes")
        primes.append(prime)
        modulus *= prime
        candidate = prime
    return tuple(primes)


def _crt_lift(modular_coefficients: Sequence[Sequence[int]], primes: Sequence[int], height_bound: int) -> tuple[tuple[int, ...], int]:
    """Lift coefficient vectors with an a-priori height proof, not heuristic CRT."""

    if not modular_coefficients or len(modular_coefficients) != len(primes):
        raise ArithmeticError("CRT inputs have incompatible lengths")
    count = len(modular_coefficients[0])
    if any(len(row) != count for row in modular_coefficients):
        raise ArithmeticError("CRT coefficient-vector lengths disagree")
    residues = [0] * count
    modulus = 1
    for prime, coefficients in zip(primes, modular_coefficients, strict=True):
        inverse = pow(modulus % prime, -1, prime)
        for index, target in enumerate(coefficients):
            correction = (int(target) - residues[index] % prime) * inverse % prime
            residues[index] += modulus * correction
        modulus *= prime
    if modulus <= 2 * height_bound:
        raise ArithmeticError("CRT modulus does not exceed the rigorous coefficient bound")
    balanced = tuple(value if 2 * value <= modulus else value - modulus for value in residues)
    return balanced, modulus


def _poly_from_shifted(coefficients: Sequence[int], valuation: int) -> sp.Poly:
    expression = sum(int(coefficient) * Q**power for power, coefficient in enumerate(coefficients))
    return sp.Poly(sp.expand(Q**valuation * expression), Q, domain=sp.ZZ)


def _evaluate_shifted(coefficients: Sequence[int], valuation: int, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + int(coefficient)
    return value**valuation * result


def _factor_record(polynomial: sp.Poly, *, timeout_seconds: int) -> dict:
    """Factor a Z[q] polynomial under a recorded wall and retain exact factor data."""

    started = time.perf_counter()
    try:
        with _hard_timeout(timeout_seconds):
            coefficient, factors = sp.factor_list(polynomial.as_expr(), Q)
        return {
            "status": "completed",
            "coefficient": str(coefficient),
            "factors": [
                {
                    "expression": str(factor),
                    "degree": int(sp.Poly(factor, Q).degree()),
                    "exponent": int(exponent),
                    "coefficients_low_to_high": [
                        str(value)
                        for value in reversed(sp.Poly(factor, Q, domain=sp.ZZ).all_coeffs())
                    ],
                }
                for factor, exponent in factors
            ],
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        }
    except FactorizationTimeout as exc:
        return {
            "status": "timed_out",
            "timeout_seconds": timeout_seconds,
            "detail": str(exc),
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        }
def _primitive_integer_poly(polynomial: sp.Poly) -> sp.Poly:
    """Normalize a nonzero rational univariate polynomial to primitive Z[q]."""

    rational = sp.Poly(polynomial, Q, domain=sp.QQ)
    denominators = [coefficient.q for coefficient in rational.all_coeffs()]
    multiple = int(sp.ilcm(*denominators)) if denominators else 1
    coefficients = [int(coefficient * multiple) for coefficient in rational.all_coeffs()]
    divisor = 0
    for coefficient in coefficients:
        divisor = math.gcd(divisor, abs(coefficient))
    coefficients = [coefficient // divisor for coefficient in coefficients]
    if coefficients[0] < 0:
        coefficients = [-coefficient for coefficient in coefficients]
    return sp.Poly.from_list(coefficients, gens=Q, domain=sp.ZZ)


def _exact_common_gcd_record(primary: sp.Poly, alternate: sp.Poly) -> tuple[sp.Poly | None, dict]:
    """Run a bounded exact Z[q] gcd computation after modular evidence is inconclusive."""

    started = time.perf_counter()
    try:
        with _hard_timeout(EXACT_GCD_WALL_SECONDS):
            common = _primitive_integer_poly(sp.gcd(primary, alternate))
        primary_quotient, primary_remainder = sp.div(primary, common, domain=sp.ZZ)
        alternate_quotient, alternate_remainder = sp.div(alternate, common, domain=sp.ZZ)
        if not primary_remainder.is_zero or not alternate_remainder.is_zero:
            raise ArithmeticError("exact gcd output does not divide both full minors")
        coefficients = [int(value) for value in reversed(common.all_coeffs())]
        return common, {
            "status": "completed",
            "method": "exact SymPy Z[q] gcd followed by exact Z[q] division checks",
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
            "degree": int(common.degree()),
            "coefficients_low_to_high": [str(value) for value in coefficients],
            "coefficients_sha256": hashlib.sha256(
                ",".join(str(value) for value in coefficients).encode()
            ).hexdigest(),
            "primary_quotient_degree": int(primary_quotient.degree()),
            "alternate_quotient_degree": int(alternate_quotient.degree()),
        }
    except FactorizationTimeout as exc:
        return None, {
            "status": "timed_out",
            "method": "exact SymPy Z[q] gcd",
            "timeout_seconds": EXACT_GCD_WALL_SECONDS,
            "detail": str(exc),
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        }






def _polynomial_record(minor: MinorData, coefficients: Sequence[int], modulus: int, primes: Sequence[int]) -> tuple[sp.Poly, dict]:
    polynomial = _poly_from_shifted(coefficients, minor.valuation_lower_bound)
    actual_degree = int(polynomial.degree())
    if actual_degree > minor.determinant_degree_upper_bound:
        raise ArithmeticError("CRT polynomial breaches its tropical degree bound")
    hash_input = ",".join(str(value) for value in coefficients).encode()
    return polynomial, {
        "label": minor.label,
        "sector": minor.sector,
        "rows": list(minor.rows),
        "shape": [SIZE, SIZE],
        "entry_count": len(minor.entries),
        "valuation_lower_bound": minor.valuation_lower_bound,
        "determinant_degree_upper_bound": minor.determinant_degree_upper_bound,
        "shifted_degree_upper_bound": minor.shifted_degree_upper_bound,
        "actual_degree": actual_degree,
        "actual_shifted_degree": max((index for index, value in enumerate(coefficients) if value), default=-1),
        "height_bound": str(minor.height_bound),
        "height_bound_bit_length": minor.height_bound.bit_length(),
        "crt_modulus": str(modulus),
        "crt_modulus_bit_length": modulus.bit_length(),
        "crt_prime_count": len(primes),
        "crt_primes": [str(prime) for prime in primes],
        "minimum_degree_matching_columns": list(minor.min_assignment),
        "maximum_degree_matching_columns": list(minor.max_assignment),
        "shifted_coefficients_low_to_high": [str(value) for value in coefficients],
        "shifted_coefficients_sha256": hashlib.sha256(hash_input).hexdigest(),
        "factorization": {"status": "not_attempted"},
    }


def _direct_minor_values(specifications: Sequence[MinorData], polynomials: dict[str, tuple[Sequence[int], int]]) -> list[dict]:
    """Fresh exact rational checks, independent of all modular interpolation nodes."""

    checks = []
    for q_value in (Fraction(-3, 1), Fraction(3, 7)):
        forward, reverse, _, _ = T16.cabled_products(q_value)
        denominator_scale = q_value.denominator**DENSE_SCALE_EXPONENT
        for minor in specifications:
            determinant = int(T16._sector_minor(minor.rows, forward, reverse, minor.sector).det(method="domain-ge"))
            coefficients, valuation = polynomials[minor.label]
            expected = _evaluate_shifted(coefficients, valuation, q_value) * denominator_scale
            if expected.denominator != 1 or determinant != expected.numerator:
                raise ArithmeticError(f"fresh exact determinant verification failed for {minor.label} at q={_q_text(q_value)}")
            checks.append(
                {
                    "label": minor.label,
                    "sector": minor.sector,
                    "q": _q_text(q_value),
                    "integer_scale": f"{q_value.denominator}^{DENSE_SCALE_EXPONENT}",
                    "direct_scaled_determinant": str(determinant),
                    "polynomial_times_scale": str(expected.numerator),
                    "passed": True,
                }
            )
    return checks


def _stored_control_checks(primary: dict[str, MinorData], polynomials: dict[str, tuple[Sequence[int], int]], source: dict) -> list[dict]:
    """Check the inherited q=1/4,1/3,2/5 full-rank determinants exactly."""

    checks = []
    for control in source["data"]["fixed_q_certificates"]:
        numerator, denominator = (
            map(int, control["q"].split("/")) if "/" in control["q"] else (int(control["q"]), 1)
        )
        q_value = Fraction(numerator, denominator)
        scale = denominator**DENSE_SCALE_EXPONENT
        for sector, minor in primary.items():
            coefficients, valuation = polynomials[minor.label]
            expected = _evaluate_shifted(coefficients, valuation, q_value) * scale
            stored = int(control["sectors"][sector]["determinant"])
            if expected.denominator != 1 or expected.numerator != stored or stored == 0:
                raise ArithmeticError(f"stored full-rank determinant does not match interpolation at q={control['q']}")
            checks.append(
                {
                    "q": control["q"],
                    "sector": sector,
                    "stored_scaled_determinant": str(stored),
                    "polynomial_times_scale": str(expected.numerator),
                    "passed": True,
                }
            )
    return checks


def _exact_sector_rank(q_value: Fraction, sector: str) -> int:
    """Exact Q-rank of every row in one 8192x128 parity sector."""

    forward, reverse, _, _ = T16.cabled_products(q_value)
    rows = [
        T16.component_row_sector(output_state, input_state, forward, reverse, sector)
        for output_state in range(T16.DIMENSION)
        for input_state in range(T16.DIMENSION)
        if T16._row_sector(output_state, input_state) == sector
    ]
    return int(sp.Matrix(rows).rank())


def _quadratic_root_full_rank_certificates(
    symbolic_forward: sp.Matrix, symbolic_reverse: sp.Matrix
) -> list[dict]:
    """Certify rank 256 at the two q²+1 roots by exact residue-field minors."""

    certificates = []
    for algebraic_point, residues in (
        ("i", ((5, 2), (13, 5))),
        ("-i", ((5, 3), (13, 8))),
    ):
        reductions = []
        for prime, residue in residues:
            forward, reverse, _, _ = T16.cabled_products(Fraction(residue, 1))
            sectors = {}
            for sector in T16.SECTOR_COLUMNS:
                rows, rank = T16._modular_independent_rows(
                    forward, reverse, sector, prime
                )
                if rank != SIZE:
                    raise ArithmeticError(
                        f"{algebraic_point} reduction mod {prime} has sector rank {rank}"
                    )
                minor = _minor_data(
                    f"{algebraic_point}_mod_{prime}_{sector}",
                    sector,
                    tuple(rows),
                    symbolic_forward,
                    symbolic_reverse,
                )
                determinant = _modular_determinant(minor.entries, residue, prime)
                if determinant == 0:
                    raise ArithmeticError(
                        f"{algebraic_point} selected sector minor vanished mod {prime}"
                    )
                sectors[sector] = {
                    "rows": list(rows),
                    "determinant_mod_prime": determinant,
                    "rank_lower_bound": SIZE,
                    "rank_upper_bound": SIZE,
                }
            reductions.append(
                {
                    "prime": prime,
                    "q_residue": residue,
                    "homomorphism": (
                        f"Z[i] -> F_{prime}, q={algebraic_point} -> {residue}"
                    ),
                    "sectors": sectors,
                }
            )
        certificates.append(
            {
                "claim_tag": "[THEOREM]",
                "q": algebraic_point,
                "minimal_polynomial": "q**2 + 1",
                "field": "Q(i)",
                "full_system_rank": 2 * SIZE,
                "sector_rank": SIZE,
                "method": (
                    "At each of two split-prime reductions, a recorded 128x128 "
                    "sector minor is nonzero. Thus rank is at least 128 over Q(i), "
                    "and the exact 128-column upper bound gives equality."
                ),
                "reductions": reductions,
            }
        )
    return certificates


def _rational_roots_from_factor_record(record: dict) -> tuple[Fraction, ...]:
    """Extract every rational root only when factorization completed exactly."""

    if record["status"] != "completed":
        return ()
    roots = set()
    for factor in record["factors"]:
        coefficients = [int(value) for value in factor["coefficients_low_to_high"]]
        if len(coefficients) != 2:
            continue
        constant, linear = coefficients
        roots.add(Fraction(-constant, linear))
    return tuple(sorted(roots))


def _rank_records(
    common_factor: sp.Poly,
    common_factorization: dict,
    *,
    resolved_algebraic_factors: set[str] = set(),
) -> tuple[list[dict], list[dict]]:
    """Directly decide rational roots and report any unhandled algebraic roots."""

    records = []
    unresolved = []
    if common_factorization["status"] != "completed":
        unresolved.append(
            {
                "claim_tag": "[UNRESOLVED]",
                "reason": "common-factor factorization timed out; rational roots were not enumerated",
            }
        )
        return records, unresolved
    for factor in common_factorization["factors"]:
        degree = int(factor["degree"])
        if degree > 1 and factor["expression"] not in resolved_algebraic_factors:
            unresolved.append(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "factor": factor["expression"],
                    "degree": degree,
                    "reason": "non-rational algebraic common-factor roots require number-field full-system rank elimination",
                }
            )
    for root in _rational_roots_from_factor_record(common_factorization):
        sector_ranks = {sector: _exact_sector_rank(root, sector) for sector in T16.SECTOR_COLUMNS}
        total_rank = sum(sector_ranks.values())
        if total_rank >= 2 * SIZE:
            raise ArithmeticError(f"common-factor root q={_q_text(root)} unexpectedly has full rank")
        records.append(
            {
                "q": _q_text(root),
                "claim_tag": "[COMPUTATION]",
                "sector_ranks": sector_ranks,
                "full_system_rank": total_rank,
                "nullity": 2 * SIZE - total_rank,
                "method": "exact SymPy Q-rank of every row in each 8192x128 parity sector",
            }
        )
    return records, unresolved


def _linear_multiplicity(polynomial: sp.Poly, root: int) -> int:
    """Return the exact multiplicity of q-root in a Z[q] polynomial."""

    divisor = sp.Poly(Q - root, Q, domain=sp.ZZ)
    work = polynomial
    multiplicity = 0
    while True:
        quotient, remainder = sp.div(work, divisor, domain=sp.ZZ)
        if not remainder.is_zero:
            return multiplicity
        work = quotient
        multiplicity += 1


def _linear_common_factor(primary: sp.Poly, alternate: sp.Poly) -> tuple[sp.Poly, dict[int, int]]:
    """Extract the candidate common factor supported at q=0,+1,-1 exactly."""

    multiplicities = {
        root: min(_linear_multiplicity(primary, root), _linear_multiplicity(alternate, root))
        for root in (0, 1, -1)
    }
    expression = sp.Integer(1)
    for root, multiplicity in multiplicities.items():
        expression *= (Q - root) ** multiplicity
    return sp.Poly(sp.expand(expression), Q, domain=sp.ZZ), multiplicities


def _linear_factor_record(multiplicities: dict[int, int]) -> dict:
    """Serialize the exact q, q-1, q+1 factorization without a dense factor search."""

    factors = []
    for root, expression, coefficients in (
        (0, "q", [0, 1]),
        (1, "q - 1", [-1, 1]),
        (-1, "q + 1", [1, 1]),
    ):
        if multiplicities[root]:
            factors.append(
                {
                    "expression": expression,
                    "degree": 1,
                    "exponent": multiplicities[root],
                    "coefficients_low_to_high": [str(value) for value in coefficients],
                }
            )
    return {
        "status": "completed",
        "method": "exact repeated Z[q] division by q, q-1, and q+1",
        "coefficient": "1",
        "factors": factors,
    }


def _modular_bezout_record(primary: sp.Poly, alternate: sp.Poly, common: sp.Poly) -> dict:
    """Seek a safe-prime Bézout certificate for the exact residual gcd."""

    reduced_primary, remainder_primary = sp.div(primary, common, domain=sp.ZZ)
    reduced_alternate, remainder_alternate = sp.div(alternate, common, domain=sp.ZZ)
    if not remainder_primary.is_zero or not remainder_alternate.is_zero:
        raise ArithmeticError("candidate common factor does not exactly divide both full minors")
    # galoistools keeps coefficients as dense Python integer lists in
    # high-to-low degree order and avoids dense symbolic-expression expansion.
    candidate = PRIME_CEILING
    rejected_primes = []
    final_prime = None
    final_gcd = None
    for _ in range(BEZOUT_PRIME_ATTEMPTS):
        prime = int(sp.prevprime(candidate))
        candidate = prime
        leading = int(reduced_primary.LC()) % prime
        if not leading:
            rejected_primes.append({"prime": prime, "reason": "reduced primary leading coefficient vanishes"})
            continue
        primary_mod = [int(value) % prime for value in reduced_primary.all_coeffs()]
        alternate_mod = [int(value) % prime for value in reduced_alternate.all_coeffs()]
        coefficient_a, coefficient_b, gcd = gf_gcdex(primary_mod, alternate_mod, prime, ZZ)
        if len(gcd) != 1:
            rejected_primes.append(
                {
                    "prime": prime,
                    "reason": "nonconstant residual gcd modulo prime",
                    "residual_gcd_degree": len(gcd) - 1,
                }
            )
            final_prime = prime
            final_gcd = gcd
            continue
        identity = gf_add(
            gf_mul(coefficient_a, primary_mod, prime, ZZ),
            gf_mul(coefficient_b, alternate_mod, prime, ZZ),
            prime,
            ZZ,
        )
        if identity != [1]:
            raise ArithmeticError("dense modular Bezout identity did not normalize to one")
        return {
            "status": "completed",
            "prime": prime,
            "reduced_primary_leading_coefficient_mod_prime": leading,
            "modular_residual_gcd": "1",
            "bezout_a_coefficients_low_to_high_mod_prime": [
                str(int(value)) for value in reversed(coefficient_a)
            ],
            "bezout_b_coefficients_low_to_high_mod_prime": [
                str(int(value)) for value in reversed(coefficient_b)
            ],
            "identity_verified_exactly_over_Fp": True,
            "rejected_modular_primes": rejected_primes,
            "logic": (
                "Any nonconstant primitive Q[q] common divisor of the reduced minors has "
                "leading coefficient dividing that of reduced_primary. The recorded leading "
                "coefficient is nonzero mod p, so such a divisor remains nonconstant mod p, "
                "contradicting the displayed F_p Bezout identity."
            ),
        }
    if final_prime is None or final_gcd is None:
        raise ArithmeticError("no safe prime remained after leading-coefficient exclusions")
    gcd_expression = sum(
        int(value) * Q ** (len(final_gcd) - 1 - power)
        for power, value in enumerate(final_gcd)
    )
    return {
        "status": "unresolved_nonconstant_modular_residual_gcd",
        "prime": final_prime,
        "modular_residual_gcd": str(gcd_expression),
        "modular_residual_gcd_degree": len(final_gcd) - 1,
        "modular_residual_gcd_coefficients_high_to_low_mod_prime": [
            str(int(value)) for value in final_gcd
        ],
        "rejected_modular_primes": rejected_primes,
        "logic": (
            "Every tested safe prime retained a nonconstant residual gcd; this is "
            "evidence, not an exact Q[q] gcd computation."
        ),
    }


def build_artifact() -> dict:
    """Produce exact interpolated minors and the complementary-minor locus certificate."""

    started = time.perf_counter()
    source = json.loads(SOURCE.read_text())
    fixed = source["data"]["fixed_q_certificates"][0]
    q_control = Fraction(1, 4)
    forward, reverse, _, _ = T16.cabled_products(q_control)
    symbolic_forward, symbolic_reverse = T16._symbolic_products(Q)
    primary_rows = {
        sector: tuple(int(row) for row in fixed["sectors"][sector]["rows"])
        for sector in T16.SECTOR_COLUMNS
    }
    alternate_rows = {
        sector: _alternate_rows(sector, forward, reverse)
        for sector in T16.SECTOR_COLUMNS
    }
    if any(set(primary_rows[sector]) & set(alternate_rows[sector]) for sector in T16.SECTOR_COLUMNS):
        raise ArithmeticError("reverse-order row selection did not produce a disjoint complementary minor")
    specifications = [
        _minor_data("primary_preserving", "parity_preserving", primary_rows["parity_preserving"], symbolic_forward, symbolic_reverse),
        _minor_data("primary_reversing", "parity_reversing", primary_rows["parity_reversing"], symbolic_forward, symbolic_reverse),
        _minor_data("alternate_preserving", "parity_preserving", alternate_rows["parity_preserving"], symbolic_forward, symbolic_reverse),
        _minor_data("alternate_reversing", "parity_reversing", alternate_rows["parity_reversing"], symbolic_forward, symbolic_reverse),
    ]
    height_bound = max(specification.height_bound for specification in specifications)
    primes = _primes_for_height(height_bound)
    modular_coefficients: dict[str, list[tuple[int, ...]]] = {specification.label: [] for specification in specifications}
    interpolation_started = time.perf_counter()
    for prime in primes:
        for specification in specifications:
            modular_coefficients[specification.label].append(_interpolate_mod_prime(specification, prime))
    interpolation_elapsed = time.perf_counter() - interpolation_started
    polynomials: dict[str, sp.Poly] = {}
    coefficient_data: dict[str, tuple[tuple[int, ...], int]] = {}
    records: dict[str, dict] = {}
    for specification in specifications:
        coefficients, modulus = _crt_lift(
            modular_coefficients[specification.label], primes, specification.height_bound
        )
        polynomial, record = _polynomial_record(specification, coefficients, modulus, primes)
        polynomials[specification.label] = polynomial
        coefficient_data[specification.label] = (coefficients, specification.valuation_lower_bound)
        records[specification.label] = record
    primary = {
        sector: next(specification for specification in specifications if specification.label == f"primary_{'preserving' if sector == 'parity_preserving' else 'reversing'}")
        for sector in T16.SECTOR_COLUMNS
    }
    direct_checks = _direct_minor_values(specifications, coefficient_data)
    control_checks = _stored_control_checks(primary, coefficient_data, source)
    primary_full = sp.Poly(
        polynomials["primary_preserving"].as_expr() * polynomials["primary_reversing"].as_expr(),
        Q,
        domain=sp.ZZ,
    )
    alternate_full = sp.Poly(
        polynomials["alternate_preserving"].as_expr() * polynomials["alternate_reversing"].as_expr(),
        Q,
        domain=sp.ZZ,
    )
    candidate_common, common_multiplicities = _linear_common_factor(primary_full, alternate_full)
    candidate_factorization = _linear_factor_record(common_multiplicities)
    initial_bezout = _modular_bezout_record(primary_full, alternate_full, candidate_common)
    exact_common: sp.Poly | None = None
    exact_gcd_record: dict
    final_bezout = initial_bezout
    common_factorization = candidate_factorization
    exact_common_proven = False
    unresolved: list[dict] = []
    if initial_bezout["status"] == "completed":
        exact_common = candidate_common
        exact_common_proven = True
        exact_gcd_record = {
            "status": "certified_by_safe_prime",
            "method": "candidate exact divisions plus safe-prime residual Bezout identity",
            "degree": int(candidate_common.degree()),
            "coefficients_low_to_high": [
                str(value) for value in reversed(candidate_common.all_coeffs())
            ],
        }
    else:
        exact_common, exact_gcd_record = _exact_common_gcd_record(primary_full, alternate_full)
        if exact_common is None:
            unresolved.append(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "reason": "bounded exact Z[q] gcd did not complete after every tested safe prime retained a nonconstant residual gcd",
                    "wall_record": exact_gcd_record,
                }
            )
        else:
            candidate_quotient, candidate_remainder = sp.div(
                exact_common, candidate_common, domain=sp.ZZ
            )
            if not candidate_remainder.is_zero:
                raise ArithmeticError("exact gcd omitted a certified linear common factor")
            final_bezout = _modular_bezout_record(
                primary_full, alternate_full, exact_common
            )
            if final_bezout["status"] != "completed":
                raise ArithmeticError(
                    "exact gcd completed but no safe-prime residual Bezout verifier was found"
                )
            exact_common_proven = True
            common_factorization = _factor_record(
                exact_common, timeout_seconds=COMMON_FACTOR_WALL_SECONDS
            )
            if common_factorization["status"] != "completed":
                unresolved.append(
                    {
                        "claim_tag": "[UNRESOLVED]",
                        "reason": "exact common gcd was reconstructed but its factorization timed out",
                        "wall_record": common_factorization,
                    }
                )
    algebraic_full_rank_certificates: list[dict] = []
    resolved_algebraic_factors: set[str] = set()
    if exact_common_proven and common_factorization["status"] == "completed":
        quadratic_factors = [
            factor
            for factor in common_factorization["factors"]
            if factor["expression"] == "q**2 + 1"
        ]
        if quadratic_factors:
            algebraic_full_rank_certificates = _quadratic_root_full_rank_certificates(
                symbolic_forward, symbolic_reverse
            )
            if len(algebraic_full_rank_certificates) != 2:
                raise ArithmeticError("did not obtain both q**2+1 root certificates")
            resolved_algebraic_factors.add("q**2 + 1")
    rank_factorization = (
        common_factorization
        if exact_common_proven and common_factorization["status"] == "completed"
        else candidate_factorization
    )
    rank_common = exact_common if rank_factorization is common_factorization else candidate_common
    rank_records, rank_unresolved = _rank_records(
        rank_common,
        rank_factorization,
        resolved_algebraic_factors=resolved_algebraic_factors,
    )
    unresolved.extend(rank_unresolved)
    if exact_common_proven and common_factorization["status"] == "completed":
        complete_locus = not unresolved
        status = (
            "[THEOREM] The full cabled 16384x256 system has rank 256 exactly for "
            "q not in {-1,0,1}. The q=0,+1,-1 ranks are computed over Q; the "
            "only other exact-gcd roots q=+i,-i have two-prime full-rank squeezes."
            if complete_locus
            else "[THEOREM] The full cabled 16384x256 system has rank 256 outside the "
            "roots of the exact complementary-minor gcd; [UNRESOLVED] some algebraic "
            "gcd roots are not directly ranked."
        )
    else:
        status = (
            "[COMPUTATION] Four exact sector-minor polynomials and the q=0,+1,-1 "
            "rank-drop certificates were obtained. [UNRESOLVED] the exact full "
            "complementary-minor gcd did not complete within its recorded wall."
        )
    checks = [
        {
            "name": "four sector determinant polynomials are exactly reconstructed by modular interpolation and height-bounded CRT",
            "passed": all(
                int(record["crt_modulus"]) > 2 * int(record["height_bound"])
                for record in records.values()
            ),
        },
        {
            "name": "fresh exact rational determinant evaluations agree with every interpolated sector polynomial",
            "passed": all(check["passed"] for check in direct_checks),
        },
        {
            "name": "stored q=1/4,1/3,2/5 full-rank sector determinants agree with primary interpolants",
            "passed": all(check["passed"] for check in control_checks),
        },
        {
            "name": "a safe-prime complementary-minor Bezout identity certifies the exact common factor",
            "passed": exact_common_proven
            and final_bezout["status"] == "completed"
            and final_bezout["modular_residual_gcd"] == "1"
            and final_bezout["identity_verified_exactly_over_Fp"],
        },
        {
            "name": "both q**2+1 roots have two-prime full-rank sector squeezes",
            "passed": (
                not resolved_algebraic_factors
                or len(algebraic_full_rank_certificates) == 2
                and all(
                    len(certificate["reductions"]) == 2
                    and all(
                        sector["rank_lower_bound"] == SIZE
                        and sector["rank_upper_bound"] == SIZE
                        for reduction in certificate["reductions"]
                        for sector in reduction["sectors"].values()
                    )
                    for certificate in algebraic_full_rank_certificates
                )
            ),
        },
    ]
    core_failures = [check["name"] for check in checks[:3] if not check["passed"]]
    if core_failures:
        raise ArithmeticError("tetra16 locus core certificate check failed: " + "; ".join(core_failures))
    return {
        "meta": {
            "provenance": "experiments/e95_tetra16_locus.py",
            "predecessor": "experiments/e83_tetra16.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "arithmetic": "Z[q] entries, exact F_p determinant evaluations, exact finite-field interpolation, CRT with a proved integer coefficient bound, and exact Q-rank at named rational specializations; NumPy int64 operations are range-bounded below 2^63.",
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        },
        "data": {
            "source_artifact": "results/integrability/tetra16.json",
            "selected_rows": {
                "primary": {sector: list(rows) for sector, rows in primary_rows.items()},
                "alternate_reverse_scan": {sector: list(rows) for sector, rows in alternate_rows.items()},
            },
            "sector_degree_bounds": {
                specification.label: {
                    "sector": specification.sector,
                    "entry_degree_upper_bound": max(entry.degree for entry in specification.entries),
                    "valuation_lower_bound": specification.valuation_lower_bound,
                    "determinant_degree_upper_bound": specification.determinant_degree_upper_bound,
                    "shifted_degree_upper_bound": specification.shifted_degree_upper_bound,
                    "degree_bound_logic": "Every nonzero determinant term is a perfect support matching. The minimum/maximum sums of entry valuations/degrees over matchings bound its valuation/degree.",
                }
                for specification in specifications
            },
            "interpolation": {
                "status": "completed",
                "method": "For each q^-v det(M(q)), evaluate at q=1,...,D+1 modulo each CRT prime; Newton-interpolate over F_p; lift every coefficient with a height-bound CRT.",
                "nodes": "integer q=1,...,D+1 separately for each shifted determinant",
                "per_point_exact_method": "modular Gaussian elimination of the evaluated 128x128 minor",
                "prime_ceiling": PRIME_CEILING,
                "int64_safety": "prime^2 < 2^63; every modular elimination product is reduced before multiplication by a row entry",
                "interpolation_elapsed_seconds": f"{interpolation_elapsed:.6f}",
                "degree_bound": max(specification.shifted_degree_upper_bound for specification in specifications),
                "per_point_wall_seconds": "not applicable: all interpolation points completed; fresh exact Bareiss checks are stored separately",
            },
            "minor_determinants": records,
            "fresh_evaluations": direct_checks,
            "stored_full_rank_cross_checks": control_checks,
            "complementary_minor_locus": {
                "primary_full_degree": int(primary_full.degree()),
                "alternate_full_degree": int(alternate_full.degree()),
                "candidate_linear_common_factor_degree": int(candidate_common.degree()),
                "candidate_linear_common_factor_coefficients_low_to_high": [
                    str(value) for value in reversed(candidate_common.all_coeffs())
                ],
                "candidate_linear_common_factorization": candidate_factorization,
                "initial_candidate_bezout": initial_bezout,
                "exact_common_gcd": exact_gcd_record,
                "exact_common_factor_coefficients_low_to_high": (
                    [str(value) for value in reversed(exact_common.all_coeffs())]
                    if exact_common_proven and exact_common is not None
                    else None
                ),
                "exact_common_factorization": (
                    common_factorization if exact_common_proven else None
                ),
                "bezout": final_bezout,
                "logic": (
                    "If the full system has rank below 256, both stored 256x256 "
                    "block-diagonal minors vanish. Their exact gcd must vanish; outside "
                    "its roots, at least one full minor is nonzero and the rank is 256. "
                    "The linear candidate is exact only after the displayed safe-prime "
                    "Bezout identity or a bounded exact-gcd computation plus its own "
                    "safe-prime residual identity."
                ),
            },
            "exceptional_specializations": rank_records,
            "algebraic_full_rank_specializations": algebraic_full_rank_certificates,
            "unresolved_algebraic_remainders": unresolved,
            "status": status,
        },
        "checks": checks,
    }


def main() -> None:
    artifact = build_artifact()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    for check in artifact["checks"]:
        print(f"{check['name']}: {'PASS' if check['passed'] else 'FAIL'}")
    print("written", OUTPUT.relative_to(ROOT))
    print("PASS" if all(check["passed"] for check in artifact["checks"]) else "UNRESOLVED")


if __name__ == "__main__":
    main()
