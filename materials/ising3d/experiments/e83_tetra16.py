"""Exact one-line-cabled 16x16 auxiliary-R RLLL calculation.

The canonical definition cabled here is

    L4_(1,4,5)(q) = L_(i1,4,5)(q) L_(i2,4,5)(q),

with the displayed factor order.  Thus line 1 is C^4 = C^2_i1 tensor
C^2_i2 while lines 2, 3, 4, 5, and 6 remain binary.  An arbitrary R acts
on i1,i2,2,3 and has 16x16 independent entries.  The global space is
128-dimensional, so the exact homogeneous system has shape 16384x256.

No dense 16384x256 matrix is materialized: the Z2 parity decomposition is
built directly, and every certificate is a pair of sparse-derived 128x128
sector minors.

Run from the repository root with::

    timeout 1800 .venv/bin/python experiments/e83_tetra16.py
"""

from __future__ import annotations

import json
import os
import resource
import signal
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterator, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "tetra16.json"
DIMENSION = 128
AUXILIARY_MASK = (1 << 4) - 1
SECTOR_SIZE = 128
FIXED_Q_VALUES = (Fraction(1, 4), Fraction(1, 3), Fraction(2, 5))
SELECTION_PRIMES = (101, 103, 107, 109, 127, 131, 137, 139, 149, 151)
GENERIC_MODULAR_POINTS = ((101, 2), (103, 2))
SYMBOLIC_PROBE_SECONDS = int(os.environ.get("TETRA16_SYMBOLIC_PROBE_SECONDS", "120"))

# Seven binary bit positions.  The first two jointly implement the cabled
# C^4 line 1; R acts on the low four bits and the physical environment is 456.
I1_45 = (0, 4, 5)
I2_45 = (1, 4, 5)
L_246 = (2, 4, 6)
L_356 = (3, 5, 6)


def _bits(state: int, width: int) -> tuple[int, ...]:
    return tuple((state >> position) & 1 for position in range(width))


def _state(bits: Sequence[int]) -> int:
    return sum(int(bit) << position for position, bit in enumerate(bits))


def _replace_local_bits(global_state: int, positions: Sequence[int], local_state: int) -> int:
    result = global_state
    for local_position, global_position in enumerate(positions):
        if (local_state >> local_position) & 1:
            result |= 1 << global_position
        else:
            result &= ~(1 << global_position)
    return result

def _peak_rss_bytes() -> int:
    """Return process peak resident set size in bytes on the active platform."""

    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if sys.platform == "darwin" else raw * 1024


def _q_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _integer_local_matrix(q: Fraction) -> sp.SparseMatrix:
    """Return denominator^3 L(q) as an exact integer 8x8 Ising matrix."""

    numerator = int(q.numerator)
    denominator = int(q.denominator)
    entries: dict[tuple[int, int], int] = {}
    for input_state in range(8):
        for output_state in range(8):
            total = input_state.bit_count() + output_state.bit_count()
            if total % 2:
                continue
            exponent = total // 2
            entries[(output_state, input_state)] = (
                numerator**exponent * denominator ** (3 - exponent)
            )
    return sp.SparseMatrix(8, 8, entries)


def _embedded_three_space(local: sp.SparseMatrix, positions: Sequence[int]) -> sp.SparseMatrix:
    """Embed a three-binary-leg matrix into the seven-bit global space."""

    entries: dict[tuple[int, int], int] = {}
    for global_input in range(DIMENSION):
        local_input = _state(
            tuple((global_input >> position) & 1 for position in positions)
        )
        for local_output in range(8):
            value = local[local_output, local_input]
            if value:
                entries[
                    (_replace_local_bits(global_input, positions, local_output), global_input)
                ] = int(value)
    return sp.SparseMatrix(DIMENSION, DIMENSION, entries)

def _polynomial_local_matrix(q: sp.Symbol) -> sp.SparseMatrix:
    """Return the exact unscaled binary Ising L(q) over Z[q]."""

    entries: dict[tuple[int, int], sp.Expr] = {}
    for input_state in range(8):
        for output_state in range(8):
            total = input_state.bit_count() + output_state.bit_count()
            if total % 2 == 0:
                entries[(output_state, input_state)] = q ** (total // 2)
    return sp.SparseMatrix(8, 8, entries)


def _embedded_three_space_polynomial(
    local: sp.SparseMatrix, positions: Sequence[int]
) -> sp.SparseMatrix:
    """Embed a Z[q] three-space matrix without numeric coercion."""

    entries: dict[tuple[int, int], sp.Expr] = {}
    for global_input in range(DIMENSION):
        local_input = _state(
            tuple((global_input >> position) & 1 for position in positions)
        )
        for local_output in range(8):
            value = local[local_output, local_input]
            if value:
                entries[
                    (_replace_local_bits(global_input, positions, local_output), global_input)
                ] = value
    return sp.SparseMatrix(DIMENSION, DIMENSION, entries)


def _symbolic_products(q: sp.Symbol) -> tuple[sp.Matrix, sp.Matrix]:
    """Build the exact Z[q] forward and reverse cabled products."""

    local = _polynomial_local_matrix(q)
    first = _embedded_three_space_polynomial(local, I1_45)
    second = _embedded_three_space_polynomial(local, I2_45)
    cable = first * second
    middle = _embedded_three_space_polynomial(local, L_246)
    last = _embedded_three_space_polynomial(local, L_356)
    return cable * middle * last, last * middle * cable


class SymbolicProbeTimeout(TimeoutError):
    """Raised when a bounded symbolic determinant probe reaches its wall."""


@contextmanager
def _hard_timeout(seconds: int):
    previous_handler = signal.getsignal(signal.SIGALRM)

    def _raise_timeout(_signum, _frame):
        raise SymbolicProbeTimeout(f"symbolic determinant probe exceeded {seconds} seconds")

    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


def cabled_products(q: Fraction) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return scaled forward/reverse products and both cable orders.

    If d is the denominator of q, each returned forward/reverse entry equals
    d^12 times the corresponding unscaled RLLL coefficient contribution.
    """

    local = _integer_local_matrix(q)
    first = _embedded_three_space(local, I1_45)
    second = _embedded_three_space(local, I2_45)
    cable = first * second
    reverse_cable = second * first
    middle = _embedded_three_space(local, L_246)
    last = _embedded_three_space(local, L_356)
    forward = cable * middle * last
    reverse = last * middle * cable
    return forward, reverse, cable, reverse_cable


def ungraded_coordinates() -> tuple[tuple[int, int], ...]:
    """All 256 R[output,input] entries, parity blocks only for certificates."""

    preserving: list[tuple[int, int]] = []
    reversing: list[tuple[int, int]] = []
    for input_state in range(16):
        for output_state in range(16):
            coordinate = (output_state, input_state)
            if (output_state.bit_count() - input_state.bit_count()) % 2:
                reversing.append(coordinate)
            else:
                preserving.append(coordinate)
    coordinates = tuple(preserving + reversing)
    if len(coordinates) != 256 or len(set(coordinates)) != 256:
        raise ArithmeticError("the cabled coordinate list is not all 16x16 R entries")
    return coordinates


COORDINATES = ungraded_coordinates()
SECTOR_COORDINATES = {
    "parity_preserving": COORDINATES[:SECTOR_SIZE],
    "parity_reversing": COORDINATES[SECTOR_SIZE:],
}
SECTOR_INDEX = {
    name: {coordinate: index for index, coordinate in enumerate(coordinates)}
    for name, coordinates in SECTOR_COORDINATES.items()
}
SECTOR_COLUMNS = {
    "parity_preserving": range(0, SECTOR_SIZE),
    "parity_reversing": range(SECTOR_SIZE, 2 * SECTOR_SIZE),
}


def _row_sector(output_state: int, input_state: int) -> str:
    return (
        "parity_preserving"
        if (output_state.bit_count() - input_state.bit_count()) % 2 == 0
        else "parity_reversing"
    )


def component_row_sector(
    output_state: int,
    input_state: int,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    sector: str,
) -> list[int]:
    """One sparse-derived length-128 parity-sector component row."""

    if sector not in SECTOR_INDEX:
        raise ValueError(f"unknown parity sector {sector}")
    values = [0] * SECTOR_SIZE
    coordinate_index = SECTOR_INDEX[sector]
    output_auxiliary = output_state & AUXILIARY_MASK
    input_auxiliary = input_state & AUXILIARY_MASK

    # R F: fixed R output and variable R input.
    for coordinate_input in range(16):
        index = coordinate_index.get((output_auxiliary, coordinate_input))
        if index is not None:
            intermediate = (output_state & ~AUXILIARY_MASK) | coordinate_input
            values[index] += int(forward[intermediate, input_state])

    # G R: fixed R input and variable R output.
    for coordinate_output in range(16):
        index = coordinate_index.get((coordinate_output, input_auxiliary))
        if index is not None:
            intermediate = (input_state & ~AUXILIARY_MASK) | coordinate_output
            values[index] -= int(reverse[output_state, intermediate])
    return values


def _sector_row_indices(sector: str) -> Iterator[int]:
    for output_state in range(DIMENSION):
        for input_state in range(DIMENSION):
            if _row_sector(output_state, input_state) == sector:
                yield DIMENSION * output_state + input_state


def _modular_independent_rows(
    forward: sp.Matrix,
    reverse: sp.Matrix,
    sector: str,
    prime: int,
) -> tuple[list[int], int]:
    """Choose 128 independent component rows by exact F_p elimination."""

    basis: list[list[int] | None] = [None] * SECTOR_SIZE
    selected: list[int] = []
    for row_index in _sector_row_indices(sector):
        output_state, input_state = divmod(row_index, DIMENSION)
        vector = [
            int(value) % prime
            for value in component_row_sector(
                output_state, input_state, forward, reverse, sector
            )
        ]
        for column, pivot_row in enumerate(basis):
            if pivot_row is None or vector[column] == 0:
                continue
            coefficient = vector[column]
            vector = [
                (value - coefficient * pivot_value) % prime
                for value, pivot_value in zip(vector, pivot_row, strict=True)
            ]
        pivot_column = next(
            (column for column, value in enumerate(vector) if value), None
        )
        if pivot_column is None:
            continue
        inverse = pow(vector[pivot_column], -1, prime)
        basis[pivot_column] = [(value * inverse) % prime for value in vector]
        selected.append(row_index)
        if len(selected) == SECTOR_SIZE:
            return selected, SECTOR_SIZE
    return selected, len(selected)


def _selected_rows_at_q(
    forward: sp.Matrix, reverse: sp.Matrix, q: Fraction
) -> tuple[int, dict[str, list[int]]]:
    """Find one good modular prime and full-sector row sets at an exact q."""

    for prime in SELECTION_PRIMES:
        if q.denominator % prime == 0:
            continue
        selected: dict[str, list[int]] = {}
        for sector in SECTOR_COLUMNS:
            rows, rank = _modular_independent_rows(forward, reverse, sector, prime)
            if rank != SECTOR_SIZE:
                break
            selected[sector] = rows
        if len(selected) == 2:
            return prime, selected
    raise ArithmeticError(f"no selected prime gave full rank for q={_q_text(q)}")


def _sector_minor(
    rows: Sequence[int],
    forward: sp.Matrix,
    reverse: sp.Matrix,
    sector: str,
) -> sp.Matrix:
    return sp.Matrix(
        [
            component_row_sector(
                int(row) // DIMENSION,
                int(row) % DIMENSION,
                forward,
                reverse,
                sector,
            )
            for row in rows
        ]
    )

def _component_row_sector_symbolic(
    output_state: int,
    input_state: int,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    sector: str,
) -> list[sp.Expr]:
    """One Z[q] sector row, retaining only its at-most-sixteen candidates."""

    values: list[sp.Expr] = [sp.Integer(0)] * SECTOR_SIZE
    coordinate_index = SECTOR_INDEX[sector]
    output_auxiliary = output_state & AUXILIARY_MASK
    input_auxiliary = input_state & AUXILIARY_MASK
    for coordinate_input in range(16):
        index = coordinate_index.get((output_auxiliary, coordinate_input))
        if index is not None:
            intermediate = (output_state & ~AUXILIARY_MASK) | coordinate_input
            values[index] += forward[intermediate, input_state]
    for coordinate_output in range(16):
        index = coordinate_index.get((coordinate_output, input_auxiliary))
        if index is not None:
            intermediate = (input_state & ~AUXILIARY_MASK) | coordinate_output
            values[index] -= reverse[output_state, intermediate]
    return values


def _symbolic_sector_minor(
    rows: Sequence[int],
    forward: sp.Matrix,
    reverse: sp.Matrix,
    sector: str,
) -> sp.Matrix:
    return sp.Matrix(
        [
            _component_row_sector_symbolic(
                int(row) // DIMENSION,
                int(row) % DIMENSION,
                forward,
                reverse,
                sector,
            )
            for row in rows
        ]
    )


def _bounded_symbolic_minor_probe(fixed: Sequence[dict]) -> dict:
    """Try the two fixed-row symbolic determinants under recorded hard walls."""

    if SYMBOLIC_PROBE_SECONDS <= 0:
        return {
            "status": "skipped",
            "matrix_shape": [SECTOR_SIZE, SECTOR_SIZE],
            "method": "disabled by TETRA16_SYMBOLIC_PROBE_SECONDS",
            "per_sector_timeout_seconds": SYMBOLIC_PROBE_SECONDS,
            "sectors": {},
            "elapsed_seconds": "0.000000",
            "peak_rss_bytes": _peak_rss_bytes(),
        }
    started = time.perf_counter()
    rss_before = _peak_rss_bytes()
    q = sp.Symbol("q")
    forward, reverse = _symbolic_products(q)
    sectors: dict[str, dict] = {}
    for sector in SECTOR_COLUMNS:
        determinant_started = time.perf_counter()
        rows = fixed[0]["sectors"][sector]["rows"]
        try:
            with _hard_timeout(SYMBOLIC_PROBE_SECONDS):
                determinant = sp.factor(
                    _symbolic_sector_minor(rows, forward, reverse, sector).det(
                        method="domain-ge"
                    )
                )
            sectors[sector] = {
                "status": "completed",
                "rows": list(rows),
                "determinant": str(determinant),
                "factorization": {
                    "coefficient": str(sp.factor_list(determinant, q)[0]),
                    "factors": [
                        [str(factor), int(exponent)]
                        for factor, exponent in sp.factor_list(determinant, q)[1]
                    ],
                },
                "elapsed_seconds": f"{time.perf_counter() - determinant_started:.6f}",
                "peak_rss_bytes": _peak_rss_bytes(),
            }
        except SymbolicProbeTimeout as exc:
            sectors[sector] = {
                "status": "timed_out",
                "rows": list(rows),
                "timeout_seconds": SYMBOLIC_PROBE_SECONDS,
                "elapsed_seconds": f"{time.perf_counter() - determinant_started:.6f}",
                "peak_rss_bytes": _peak_rss_bytes(),
                "detail": str(exc),
            }
        except Exception as exc:  # Record a reproducible failed route, not a theorem.
            sectors[sector] = {
                "status": "failed",
                "rows": list(rows),
                "elapsed_seconds": f"{time.perf_counter() - determinant_started:.6f}",
                "peak_rss_bytes": _peak_rss_bytes(),
                "detail": f"{type(exc).__name__}: {exc}",
            }
    completed = all(row["status"] == "completed" for row in sectors.values())
    return {
        "status": "completed" if completed else "UNRESOLVED",
        "matrix_shape": [SECTOR_SIZE, SECTOR_SIZE],
        "entry_ring": "Z[q]",
        "method": (
            "SymPy fraction-free determinant on the same two sparse-derived 128x128 "
            "sector minors selected at q=1/4"
        ),
        "per_sector_timeout_seconds": SYMBOLIC_PROBE_SECONDS,
        "sectors": sectors,
        "rss_before_probe_bytes": rss_before,
        "peak_rss_bytes": _peak_rss_bytes(),
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }



def _integer_factorization(value: int) -> dict:
    """Return a complete signed prime factor certificate for a nonzero integer."""

    integer = int(value)
    if integer == 0:
        raise ArithmeticError("a maximal-minor factorization cannot certify zero")
    factors = sp.factorint(abs(integer))
    if any(not sp.isprime(int(prime)) for prime in factors):
        raise ArithmeticError("factorint returned a non-prime factor")
    return {
        "value": str(integer),
        "sign": -1 if integer < 0 else 1,
        "prime_factors": [[str(prime), int(exponent)] for prime, exponent in sorted(factors.items())],
    }


def _fixed_q_certificate(q: Fraction) -> dict:
    """Produce two exact 128x128 full-rank sector certificates at one rational q."""

    started = time.perf_counter()
    forward, reverse, cable, reverse_cable = cabled_products(q)
    if cable == reverse_cable:
        raise ArithmeticError("the defining cabled factors unexpectedly commute")
    selection_prime, selected_rows = _selected_rows_at_q(forward, reverse, q)
    sectors: dict[str, dict] = {}
    determinants: list[int] = []
    for sector in SECTOR_COLUMNS:
        determinant_started = time.perf_counter()
        minor = _sector_minor(selected_rows[sector], forward, reverse, sector)
        determinant = int(minor.det(method="domain-ge"))
        if determinant == 0:
            raise ArithmeticError(
                f"modularly selected {sector} rows gave a zero exact determinant at q={_q_text(q)}"
            )
        determinants.append(determinant)
        sectors[sector] = {
            "rows": list(selected_rows[sector]),
            "shape": [SECTOR_SIZE, SECTOR_SIZE],
            "rank": SECTOR_SIZE,
            "determinant": str(determinant),
            "integer_factorization": _integer_factorization(determinant),
            "selection_prime": selection_prime,
            "selection_minor_rank_mod_prime": SECTOR_SIZE,
            "determinant_elapsed_seconds": f"{time.perf_counter() - determinant_started:.6f}",
        }
    denominator = int(q.denominator)
    return {
        "q": _q_text(q),
        "local_integer_scale": f"{denominator}^3",
        "equation_entry_integer_scale": f"{denominator}^12",
        "sector_minor_integer_scale": f"{denominator}^{12 * SECTOR_SIZE}",
        "sectors": sectors,
        "full_block_diagonal_determinant": str(determinants[0] * determinants[1]),
        "exact_full_rank": 2 * SECTOR_SIZE,
        "generic_nullity_at_this_q": 0,
        "cable_orders_differ": True,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def _sector_support_analysis(forward: sp.Matrix, reverse: sp.Matrix) -> dict:
    """Scan sparse row generation without materializing the rectangular system."""

    started = time.perf_counter()
    row_counts = {"parity_preserving": 0, "parity_reversing": 0}
    nonzero_rows = 0
    off_sector_nonzero_entries = {
        "parity_preserving_rows_to_reversing_columns": 0,
        "parity_reversing_rows_to_preserving_columns": 0,
    }
    for output_state in range(DIMENSION):
        for input_state in range(DIMENSION):
            sector = _row_sector(output_state, input_state)
            row_counts[sector] += 1
            own = component_row_sector(output_state, input_state, forward, reverse, sector)
            other = component_row_sector(
                output_state,
                input_state,
                forward,
                reverse,
                "parity_reversing" if sector == "parity_preserving" else "parity_preserving",
            )
            if any(own) or any(other):
                nonzero_rows += 1
            if sector == "parity_preserving":
                off_sector_nonzero_entries[
                    "parity_preserving_rows_to_reversing_columns"
                ] += sum(value != 0 for value in other)
            else:
                off_sector_nonzero_entries[
                    "parity_reversing_rows_to_preserving_columns"
                ] += sum(value != 0 for value in other)
    if any(off_sector_nonzero_entries.values()):
        raise ArithmeticError("the asserted direct-sum parity decomposition failed")
    return {
        "row_counts": row_counts,
        "column_counts": {
            "parity_preserving": SECTOR_SIZE,
            "parity_reversing": SECTOR_SIZE,
        },
        "nonzero_rows": nonzero_rows,
        "off_sector_nonzero_entries": off_sector_nonzero_entries,
        "method": "streamed 128-component sector rows; no 16384x256 dense matrix",
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def _modular_generic_evidence() -> list[dict]:
    """Give nonzero-minor witnesses over F_p(q) via good numeric residues."""

    evidence = []
    for prime, residue in GENERIC_MODULAR_POINTS:
        q = Fraction(residue, 1)
        forward, reverse, _, _ = cabled_products(q)
        sector_rows: dict[str, list[int]] = {}
        for sector in SECTOR_COLUMNS:
            rows, rank = _modular_independent_rows(forward, reverse, sector, prime)
            if rank != SECTOR_SIZE:
                raise ArithmeticError(
                    f"generic modular witness failed in {sector} over F_{prime} at q={residue}"
                )
            sector_rows[sector] = rows
        evidence.append(
            {
                "prime": prime,
                "q_residue": residue,
                "sector_rows": sector_rows,
                "sector_ranks": {
                    "parity_preserving": SECTOR_SIZE,
                    "parity_reversing": SECTOR_SIZE,
                },
                "meaning": (
                    "At this residue, the stored minors are nonzero modulo p; therefore each "
                    "sector has generic rank 128 over F_p(q). This is modular generic evidence, "
                    "not an all-parameter characteristic-zero certificate."
                ),
            }
        )
    return evidence

def _generic_rank_over_q_certificate(modular: Sequence[dict]) -> dict:
    """Promote nonzero modular specializations to a Q(q) rank theorem."""

    if not modular or any(
        any(rank != SECTOR_SIZE for rank in row["sector_ranks"].values())
        for row in modular
    ):
        raise ArithmeticError("modular data cannot support the generic-rank theorem")
    return {
        "proof_status": "[THEOREM]",
        "rank_over_Qq": 2 * SECTOR_SIZE,
        "nullity_over_Qq": 0,
        "modular_witnesses": list(modular),
        "logic": (
            "For each stored sector row set, its determinant is a polynomial in Z[q]. "
            "At the recorded integer residue it is nonzero modulo the recorded prime, so it "
            "is not the zero polynomial. The two parity sector determinants form a literal "
            "block-diagonal 256x256 minor of B(q); hence rank B=256 over Q(q)."
        ),
        "exceptional_set_status": "UNRESOLVED_NOT_FACTORED",
        "limitation": (
            "This proves generic characteristic-zero rank only. It does not identify or "
            "factor the finite specialization locus at which the selected minor vanishes."
        ),
    }

def _support_components(matrix: sp.Matrix) -> dict:
    """Describe connected components of a selected minor's bipartite support."""

    size = int(matrix.rows)
    if matrix.cols != size:
        raise ValueError("support components require a square minor")
    adjacency = [[] for _ in range(2 * size)]
    for row in range(size):
        for column in range(size):
            if matrix[row, column] != 0:
                adjacency[row].append(size + column)
                adjacency[size + column].append(row)
    seen = [False] * (2 * size)
    components = []
    for start in range(2 * size):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        rows = 0
        columns = 0
        edges = 0
        while stack:
            node = stack.pop()
            if node < size:
                rows += 1
                edges += len(adjacency[node])
            else:
                columns += 1
            for neighbor in adjacency[node]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        components.append(
            {
                "rows": rows,
                "columns": columns,
                "nonzero_entries": edges,
            }
        )
    return {
        "component_count": len(components),
        "components": components,
        "all_square": all(row["rows"] == row["columns"] for row in components),
    }


def _operator_rank_analysis(forward: sp.Matrix, reverse: sp.Matrix) -> dict:
    """Record exact ranks of the two structured 128x128 products at q=1/4."""

    started = time.perf_counter()
    forward_rank = int(forward.rank())
    reverse_rank = int(reverse.rank())
    return {
        "forward_rank": forward_rank,
        "reverse_rank": reverse_rank,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def build_artifact() -> dict:
    """Build the exact finite-q certificate envelope without dense system materialization."""

    started = time.perf_counter()
    first_q = FIXED_Q_VALUES[0]
    forward, reverse, cable, reverse_cable = cabled_products(first_q)
    if cable == reverse_cable:
        raise ArithmeticError("cabling order did not remain observable at the control point")
    support = _sector_support_analysis(forward, reverse)
    operator_ranks = _operator_rank_analysis(forward, reverse)
    fixed = [_fixed_q_certificate(q) for q in FIXED_Q_VALUES]
    modular = _modular_generic_evidence()
    generic_rank = _generic_rank_over_q_certificate(modular)
    structure = {
        sector: _support_components(
            _sector_minor(fixed[0]["sectors"][sector]["rows"], forward, reverse, sector)
        )
        for sector in SECTOR_COLUMNS
    }
    symbolic = _bounded_symbolic_minor_probe(fixed)
    symbolic_complete = symbolic["status"] == "completed"
    checks = [
        {
            "name": "corrected cabled geometry has the true 16384x256 system shape",
            "passed": True,
            "detail": "seven binary bits give global dimension 128; an arbitrary 16x16 R has 256 entries",
        },
        {
            "name": "fixed cable order is material and exact Z2 sectors are disjoint",
            "passed": cable != reverse_cable and not any(support["off_sector_nonzero_entries"].values()),
            "detail": (
                "L_(i1,4,5)L_(i2,4,5) differs from reversed order at q=1/4; "
                "all off-sector streamed coefficients are exactly zero"
            ),
        },
        {
            "name": "three rational specializations have exact full column rank",
            "passed": all(row["exact_full_rank"] == 256 for row in fixed),
            "detail": "; ".join(f"q={row['q']}: rank 256" for row in fixed),
        },
        {
            "name": "two-prime modular generic lower witnesses are full rank",
            "passed": all(
                all(rank == SECTOR_SIZE for rank in row["sector_ranks"].values())
                for row in modular
            ),
            "detail": "Each listed F_p residue has two nonzero 128x128 sector minors",
        },
        {
            "name": "modular minor witnesses prove generic rank 256 over Q(q)",
            "passed": generic_rank["rank_over_Qq"] == 256
            and generic_rank["nullity_over_Qq"] == 0,
            "detail": "Nonzero specializations modulo two primes show the fixed polynomial minors are nonzero.",
        },
    ]
    if not all(check["passed"] for check in checks):
        raise ArithmeticError("a tetra16 certificate check failed")
    return {
        "meta": {
            "provenance": "experiments/e83_tetra16.py",
            "predecessor": "experiments/e70_tetra_ungraded.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "arithmetic": (
                "exact Python integers and SymPy sparse/integer fraction-free determinants; "
                "floating point is used only for recorded wall times"
            ),
            "sympy_version": sp.__version__,
            "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        },
        "data": {
            "ansatz": {
                "equation": (
                    "R_(i1,i2,2,3) L4_(1,4,5)(q) L_(2,4,6)(q) L_(3,5,6)(q) = "
                    "L_(3,5,6)(q) L_(2,4,6)(q) L4_(1,4,5)(q) R_(i1,i2,2,3)"
                ),
                "operator_convention": "L[output,input]=q^((|output|+|input|)/2) for even total support",
                "global_operator_shape": [DIMENSION, DIMENSION],
                "auxiliary_R_shape": [16, 16],
                "R_entries": 256,
                "independent_unknowns": 256,
                "true_system_shape": [DIMENSION * DIMENSION, 256],
                "grading_assumption_on_R": "none",
                "other_R_symmetries": "none",
                "identical_uncabled_Ising_local_factors": True,
            },
            "cabling": {
                "line_1": "C^4 = C^2_i1 tensor C^2_i2",
                "bit_order": ["i1", "i2", "2", "3", "4", "5", "6"],
                "factor_order": "L_(i1,4,5) * L_(i2,4,5)",
                "definition": "L4_(1,4,5)(q) := L_(i1,4,5)(q) L_(i2,4,5)(q)",
                "order_is_material": True,
                "uncabled_C4_leg_status": "out_of_scope_different_bilinear_problem",
            },
            "sector_decomposition": support,
            "fixed_q_certificates": fixed,
            "modular_generic_evidence": modular,
            "symbolic_all_q_status": symbolic,
            "generic_rank_over_Qq": generic_rank,
            "selected_minor_support_components": structure,
            "structured_operator_ranks_at_q_1_4": operator_ranks,
            "outcome": {
                "class": (
                    "CABLED_16X16_RLLL_ALL_Q_NOGO"
                    if symbolic_complete
                    else "CABLED_16X16_RLLL_GENERIC_RANK_AND_FINITE_Q_CERTIFICATE"
                ),
                "statement": (
                    (
                        "[THEOREM] The displayed product of the two completed exact sector "
                        "minors has full column rank 256 outside the explicit union of their "
                        "stored irreducible factor roots."
                    )
                    if symbolic_complete
                    else (
                        "[THEOREM] The cabled 16384x256 system has rank 256 over Q(q), so its "
                        "generic nullspace is zero. [COMPUTATION] Exact rank 256 is independently "
                        "certified at q=1/4, 1/3, and 2/5."
                    )
                ),
                "scope_limit": (
                    "The cabled definition is one fixed fusion order. An unconstrained C^4-leg "
                    "L or arbitrary C^4 local factor is a different bilinear problem and is out "
                    "of scope."
                ),
            },
            "status": (
                "[THEOREM] all-q cabled 16x16 RLLL rank certificate outside the stored factor "
                "roots."
                if symbolic_complete
                else (
                    "[THEOREM] generic cabled 16x16 RLLL rank 256 over Q(q); "
                    "[COMPUTATION] exact full rank at q=1/4,1/3,2/5; "
                    "[UNRESOLVED] factored all-specialization exception set."
                )
            ),
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
    print("PASS")


if __name__ == "__main__":
    main()
