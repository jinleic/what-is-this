#!/usr/bin/env python3
"""Exact Hom-block closure for the open 2x6 Ising-ladder commutant.

For A=sum X_v and B=sum Z_u Z_v, split the spin space by row reflection,
column reflection, and global spin flip.  B is diagonal in each orbit-character
basis, so every Hom block is reduced before any elimination.  Prime-field
nullspaces give rank lower bounds; canonically aligned bases at two primes are
CRT-combined and rationally reconstructed before exact substitution in every
integer equation.  Those two inequalities certify each rational rank exactly
without a 2^24 raw-adjoint matrix.
"""
from __future__ import annotations

import argparse
import math
import hashlib
import itertools
import json
import platform
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "commutant_2x6.json"
SCRIPT = "experiments/e89_commutant_2x6.py"
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))
LABELS = tuple("".join(str(bit) for bit in character) for character in CHARACTERS)
PRIMARY_PRIME = 1_000_000_007
SECONDARY_PRIME = 2_147_483_647

SparseRow = dict[int, int]


@dataclass(frozen=True)
class Sector:
    character: tuple[int, int, int]
    representatives: tuple[int, ...]
    a_columns: tuple[dict[int, int], ...]
    a_rows: tuple[dict[int, int], ...]
    b_diagonal: tuple[int, ...]
    norms: tuple[int, ...]

    @property
    def dimension(self) -> int:
        return len(self.representatives)


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    print(f"{'PASS' if passed else 'FAIL'}: {name} ({detail})")
    if not passed:
        raise AssertionError(name)
    return {"name": name, "passed": True, "detail": detail}


def grid_edges(columns: int) -> tuple[tuple[int, int], ...]:
    return tuple((column, columns + column) for column in range(columns)) + tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(2)
        for column in range(columns - 1)
    )


def act_configuration(state: int, element: tuple[int, int, int], columns: int) -> int:
    row_reflection, column_reflection, spin_flip = element
    output = 0
    for row in range(2):
        for column in range(columns):
            old = row * columns + column
            new_row = 1 - row if row_reflection else row
            new_column = columns - 1 - column if column_reflection else column
            new = new_row * columns + new_column
            output |= ((state >> old) & 1) << new
    if spin_flip:
        output ^= (1 << (2 * columns)) - 1
    return output


def symmetry_sector(columns: int, character: tuple[int, int, int]) -> Sector:
    sites = 2 * columns
    seen: set[int] = set()
    representatives: list[int] = []
    orbit_vectors: list[dict[int, int]] = []
    for state in range(1 << sites):
        if state in seen:
            continue
        orbit = {act_configuration(state, element, columns) for element in CHARACTERS}
        seen.update(orbit)
        projected: dict[int, int] = {}
        for element in CHARACTERS:
            sign = -1 if sum(a * b for a, b in zip(character, element)) & 1 else 1
            image = act_configuration(state, element, columns)
            projected[image] = projected.get(image, 0) + sign
        projected = {state: value for state, value in projected.items() if value}
        if not projected:
            continue
        representative = min(projected)
        scale = projected[representative]
        vector = {state: value // scale for state, value in projected.items()}
        if not all(value in (-1, 1) for value in vector.values()):
            raise AssertionError("orbit-character normalization is not integral")
        representatives.append(representative)
        orbit_vectors.append(vector)

    dimension = len(representatives)
    a_columns: list[dict[int, int]] = [dict() for _ in range(dimension)]
    a_rows: list[dict[int, int]] = [dict() for _ in range(dimension)]
    b_diagonal: list[int] = []
    edges = grid_edges(columns)
    for source, vector in enumerate(orbit_vectors):
        image: dict[int, int] = {}
        for state, coefficient in vector.items():
            for site in range(sites):
                target = state ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for target, representative in enumerate(representatives):
            coefficient = image.get(representative, 0)
            if coefficient:
                a_columns[source][target] = coefficient
                a_rows[target][source] = coefficient
        representative = representatives[source]
        cut = sum(
            ((representative >> left) ^ (representative >> right)) & 1
            for left, right in edges
        )
        b_diagonal.append(len(edges) - 2 * cut)

    norms = tuple(sum(value * value for value in vector.values()) for vector in orbit_vectors)
    sector = Sector(
        character=character,
        representatives=tuple(representatives),
        a_columns=tuple(a_columns),
        a_rows=tuple(a_rows),
        b_diagonal=tuple(b_diagonal),
        norms=norms,
    )
    if any(
        norms[row] * coefficient != norms[column] * a_columns[row].get(column, 0)
        for column, entries in enumerate(a_columns)
        for row, coefficient in entries.items()
    ):
        raise AssertionError("restricted A is not self-adjoint in the exact orbit metric")
    return sector


def block_variables(target: Sector, source: Sector) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    variables = [
        (row, column)
        for row in range(target.dimension)
        for column in range(source.dimension)
        if target.b_diagonal[row] == source.b_diagonal[column]
    ]
    return variables, {pair: index for index, pair in enumerate(variables)}


def iter_block_equations(
    target: Sector, source: Sector, variable_index: dict[tuple[int, int], int]
) -> Iterator[SparseRow]:
    """Yield XA_source - A_target X after the exact diagonal B equation."""
    for row in range(target.dimension):
        target_a_row = target.a_rows[row]
        for column in range(source.dimension):
            equation: SparseRow = {}
            for middle, coefficient in source.a_columns[column].items():
                variable = variable_index.get((row, middle))
                if variable is not None:
                    equation[variable] = equation.get(variable, 0) + coefficient
            for middle, coefficient in target_a_row.items():
                variable = variable_index.get((middle, column))
                if variable is not None:
                    equation[variable] = equation.get(variable, 0) - coefficient
            equation = {index: value for index, value in equation.items() if value}
            if equation:
                yield equation


def modular_matrix(
    rows: Iterable[SparseRow], n_columns: int, prime: int
) -> tuple[sp.polys.matrices.DomainMatrix, dict[str, int]]:
    domain = sp.GF(prime)
    data: dict[int, dict[int, object]] = {}
    count = 0
    maximum_support = 0
    for source in rows:
        row = {
            column: domain.convert(value % prime)
            for column, value in source.items()
            if value % prime
        }
        if row:
            data[count] = row
            count += 1
            maximum_support = max(maximum_support, len(row))
    return sp.polys.matrices.DomainMatrix(data, (count, n_columns), domain), {
        "nonzero_A_equations": count,
        "maximum_A_equation_support": maximum_support,
    }


def modular_rank(
    target: Sector, source: Sector, variable_index: dict[tuple[int, int], int], n_columns: int, prime: int
) -> tuple[int, dict[str, int]]:
    matrix, stats = modular_matrix(
        iter_block_equations(target, source, variable_index), n_columns, prime
    )
    return int(matrix.rank()), stats


def modular_nullspace(
    target: Sector, source: Sector, variable_index: dict[tuple[int, int], int], n_columns: int, prime: int
) -> tuple[int, list[list[int]], dict[str, int]]:
    matrix, stats = modular_matrix(
        iter_block_equations(target, source, variable_index), n_columns, prime
    )
    rref, pivots = matrix.rref()
    nullspace = rref.nullspace_from_rref(pivots).to_Matrix()
    nullity = int(nullspace.rows)
    rank = len(pivots)
    basis = [
        [int(nullspace[row, column]) % prime for column in range(n_columns)]
        for row in range(nullity)
    ]
    return rank, basis, stats


def verify_exact_kernel(
    target: Sector,
    source: Sector,
    variable_index: dict[tuple[int, int], int],
    basis: Sequence[Sequence[Fraction]],
) -> bool:
    if not basis:
        return True
    for equation in iter_block_equations(target, source, variable_index):
        residuals = [Fraction(0) for _ in basis]
        for column, coefficient in equation.items():
            for vector_index, vector in enumerate(basis):
                residuals[vector_index] += coefficient * vector[column]
        if any(residuals):
            return False
    return True


def vector_digest(basis: Sequence[Sequence[Fraction]]) -> str:
    digest = hashlib.sha256()
    for row, vector in enumerate(basis):
        for column, value in enumerate(vector):
            if value:
                digest.update(
                    f"{row}:{column}:{value.numerator}:{value.denominator};".encode()
                )
    return digest.hexdigest()


def modular_anchor_certificate(basis: Sequence[Sequence[int]], prime: int) -> tuple[list[int], int]:
    """Return independent coordinate columns and their determinant modulo prime."""
    if not basis:
        return [], 1
    dimension = len(basis)
    pivots: dict[int, list[int]] = {}
    anchors: list[int] = []
    for column in range(len(basis[0])):
        vector = [basis[row][column] % prime for row in range(dimension)]
        while any(vector):
            lead = min(index for index, value in enumerate(vector) if value)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = [value * inverse % prime for value in vector]
                anchors.append(column)
                break
            factor = vector[lead]
            vector = [
                (value - factor * old_value) % prime
                for value, old_value in zip(vector, old)
            ]
        if len(anchors) == dimension:
            break
    if len(anchors) != dimension:
        raise AssertionError("modular nullspace rows are dependent")

    matrix = [[basis[row][column] % prime for column in anchors] for row in range(dimension)]
    determinant = 1
    for pivot in range(dimension):
        source = next((row for row in range(pivot, dimension) if matrix[row][pivot]), None)
        if source is None:
            raise AssertionError("selected anchor minor is singular")
        if source != pivot:
            matrix[pivot], matrix[source] = matrix[source], matrix[pivot]
            determinant = -determinant
        pivot_value = matrix[pivot][pivot]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        for row in range(pivot + 1, dimension):
            factor = matrix[row][pivot] * inverse % prime
            for column in range(pivot, dimension):
                matrix[row][column] = (matrix[row][column] - factor * matrix[pivot][column]) % prime
    return anchors, determinant % prime


def inverse_mod_square(matrix: Sequence[Sequence[int]], prime: int) -> list[list[int]]:
    size = len(matrix)
    augmented = [
        [value % prime for value in row]
        + [1 if row_index == column else 0 for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for pivot in range(size):
        source = next((row for row in range(pivot, size) if augmented[row][pivot]), None)
        if source is None:
            raise AssertionError("modular anchor minor is singular")
        if source != pivot:
            augmented[pivot], augmented[source] = augmented[source], augmented[pivot]
        inverse = pow(augmented[pivot][pivot], prime - 2, prime)
        augmented[pivot] = [value * inverse % prime for value in augmented[pivot]]
        for row in range(size):
            if row == pivot:
                continue
            factor = augmented[row][pivot]
            if factor:
                augmented[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(augmented[row], augmented[pivot])
                ]
    return [row[size:] for row in augmented]


def determinant_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    square = [[value % prime for value in row] for row in matrix]
    determinant = 1
    for pivot in range(len(square)):
        source = next((row for row in range(pivot, len(square)) if square[row][pivot]), None)
        if source is None:
            return 0
        if source != pivot:
            square[pivot], square[source] = square[source], square[pivot]
            determinant = -determinant
        pivot_value = square[pivot][pivot]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        for row in range(pivot + 1, len(square)):
            factor = square[row][pivot] * inverse % prime
            for column in range(pivot, len(square)):
                square[row][column] = (
                    square[row][column] - factor * square[pivot][column]
                ) % prime
    return determinant % prime


def normalize_modular_basis(
    basis: Sequence[Sequence[int]], anchors: Sequence[int], prime: int
) -> list[list[int]]:
    if not basis:
        return []
    anchor_matrix = [[row[column] % prime for column in anchors] for row in basis]
    inverse = inverse_mod_square(anchor_matrix, prime)
    return [
        [
            sum(inverse[row][middle] * basis[middle][column] for middle in range(len(basis)))
            % prime
            for column in range(len(basis[0]))
        ]
        for row in range(len(basis))
    ]


def crt_pair(first: int, first_modulus: int, second: int, second_modulus: int) -> int:
    correction = (
        (second - first)
        * pow(first_modulus % second_modulus, second_modulus - 2, second_modulus)
    ) % second_modulus
    return first + first_modulus * correction


def rational_reconstruct(residue: int, modulus: int) -> Fraction:
    """Recover the unique small rational represented by residue modulo modulus."""
    residue %= modulus
    if residue == 0:
        return Fraction(0)
    bound = math.isqrt(modulus // 2)
    previous_remainder, remainder = modulus, residue
    previous_denominator, denominator = 0, 1
    while abs(remainder) > bound:
        quotient = previous_remainder // remainder
        previous_remainder, remainder = remainder, previous_remainder - quotient * remainder
        previous_denominator, denominator = denominator, previous_denominator - quotient * denominator
    if denominator == 0 or abs(denominator) > bound:
        raise AssertionError("rational reconstruction exceeded the two-prime bound")
    if denominator < 0:
        remainder, denominator = -remainder, -denominator
    value = Fraction(remainder, denominator)
    if (value.denominator * residue - value.numerator) % modulus:
        raise AssertionError("rational reconstruction congruence failed")
    return value


def reconstruct_two_prime_basis(
    primary_basis: Sequence[Sequence[int]], secondary_basis: Sequence[Sequence[int]]
) -> tuple[list[list[Fraction]], list[int], int, int]:
    if len(primary_basis) != len(secondary_basis):
        raise AssertionError("two-prime nullities disagree")
    if not primary_basis:
        return [], [], 1, 1
    anchors, primary_determinant = modular_anchor_certificate(primary_basis, PRIMARY_PRIME)
    primary = normalize_modular_basis(primary_basis, anchors, PRIMARY_PRIME)
    secondary = normalize_modular_basis(secondary_basis, anchors, SECONDARY_PRIME)
    secondary_determinant = determinant_mod(
        [[row[column] for column in anchors] for row in secondary_basis],
        SECONDARY_PRIME,
    )
    if not secondary_determinant:
        raise AssertionError("primary nullspace anchor minor vanishes at the secondary prime")
    modulus = PRIMARY_PRIME * SECONDARY_PRIME
    basis = [
        [
            rational_reconstruct(
                crt_pair(
                    primary[row][column],
                    PRIMARY_PRIME,
                    secondary[row][column],
                    SECONDARY_PRIME,
                ),
                modulus,
            )
            for column in range(len(primary[0]))
        ]
        for row in range(len(primary))
    ]
    if any(
        basis[row][anchor] != (Fraction(1) if row == column else Fraction(0))
        for row in range(len(basis))
        for column, anchor in enumerate(anchors)
    ):
        raise AssertionError("two-prime basis lost its canonical anchor normalization")
    return basis, anchors, primary_determinant, secondary_determinant


def maximum_fraction_bits(basis: Sequence[Sequence[Fraction]]) -> int:
    return max(
        (
            max(abs(value.numerator).bit_length(), value.denominator.bit_length())
            for vector in basis
            for value in vector
            if value
        ),
        default=0,
    )


def fraction_echelon(rows: Iterable[dict[int, int | Fraction]]) -> dict[int, dict[int, Fraction]]:
    pivots: dict[int, dict[int, Fraction]] = {}
    for source in rows:
        row = {column: Fraction(value) for column, value in source.items() if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                pivot = row[lead]
                pivots[lead] = {column: value / pivot for column, value in row.items() if value}
                break
            factor = row[lead]
            for column, value in old.items():
                reduced = row.get(column, Fraction(0)) - factor * value
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return pivots


def fraction_nullspace(
    pivots: dict[int, dict[int, Fraction]], n_columns: int
) -> list[dict[int, Fraction]]:
    free = [column for column in range(n_columns) if column not in pivots]
    basis: list[dict[int, Fraction]] = []
    for free_column in free:
        vector: dict[int, Fraction] = {free_column: Fraction(1)}
        for pivot in sorted(pivots, reverse=True):
            value = -sum(
                coefficient * vector.get(column, Fraction(0))
                for column, coefficient in pivots[pivot].items()
                if column != pivot
            )
            if value:
                vector[pivot] = value
        basis.append(vector)
    return basis


def common_zero_dimensions(sectors: Sequence[Sector]) -> tuple[list[int], list[str]]:
    dimensions: list[int] = []
    digests: list[str] = []
    for sector in sectors:
        zero_columns = [column for column, energy in enumerate(sector.b_diagonal) if energy == 0]
        variable = {column: index for index, column in enumerate(zero_columns)}
        equations: list[SparseRow] = []
        for row in sector.a_rows:
            equation = {
                variable[column]: coefficient
                for column, coefficient in row.items()
                if column in variable
            }
            if equation:
                equations.append(equation)
        pivots = fraction_echelon(equations)
        basis = fraction_nullspace(pivots, len(zero_columns))
        if any(
            sum(Fraction(coefficient) * vector.get(column, Fraction(0)) for column, coefficient in equation.items())
            for equation in equations
            for vector in basis
        ):
            raise AssertionError("common-zero basis fails exact substitution")
        dimensions.append(len(basis))
        digest = hashlib.sha256()
        for vector_index, vector in enumerate(basis):
            for column, value in sorted(vector.items()):
                digest.update(
                    f"{vector_index}:{zero_columns[column]}:{value.numerator}:{value.denominator};".encode()
                )
        digests.append(digest.hexdigest())
    return dimensions, digests


def subtract_common_zero_block(
    table: Sequence[Sequence[int]], common_zero: Sequence[int]
) -> list[list[int]]:
    residual = [
        [
            table[target][source] - common_zero[target] * common_zero[source]
            for source in range(len(table))
        ]
        for target in range(len(table))
    ]
    if any(value < 0 for row in residual for value in row):
        raise AssertionError("common-zero Hom contribution exceeds an exact Hom block")
    return residual


def infer_control_type(closure: dict[str, object]) -> dict[str, object]:
    table = closure["hom_dimension_matrix_Q_rows_target_columns_source"]
    if not isinstance(table, list):
        raise AssertionError("control Hom table is malformed")
    expected = [
        [2, 0, 1, 0, 0, 0, 0, 0],
        [0, 3, 0, 0, 0, 0, 0, 0],
        [1, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 2],
    ]
    if table != expected:
        raise AssertionError("same-pipeline 2x5 control Hom table changed")
    if closure["common_A0_B0_dimensions_by_character"] != [0] * 8:
        raise AssertionError("2x5 has an unexpected common-zero space")
    multiplicities = [2, 2] + [1] * 12
    if sum(value * value for value in multiplicities) != 20:
        raise AssertionError("2x5 Wedderburn count failed")
    return {
        "claim_tag": "[COMPUTATION]",
        "commutant_dimension_Q_equals_dimension_C": 20,
        "wedderburn_type_over_C": "M_2(C)^2 direct-sum C^12",
        "centre_dimension_over_C": 14,
        "exact_Hom_table_reproduced": True,
        "type_derivation": (
            "The [000,010] residual table [[2,1],[1,2]] gives M_2(C) plus C^2; "
            "[100,110] gives M_2(C); the four isolated diagonal blocks of dimensions "
            "3,2,3,2 give C^10."
        ),
    }


def infer_2x6_type(closure: dict[str, object]) -> dict[str, object]:
    table = closure["hom_dimension_matrix_Q_rows_target_columns_source"]
    common_zero = closure["common_A0_B0_dimensions_by_character"]
    if not isinstance(table, list) or not isinstance(common_zero, list):
        raise AssertionError("2x6 closure payload is malformed")
    if common_zero != [1, 0, 2, 0, 2, 0, 2, 0]:
        raise AssertionError("2x6 common-zero dimensions changed")
    residual = subtract_common_zero_block(table, common_zero)
    expected_residual = [
        [2, 0, 1, 0, 0, 0, 0, 0],
        [0, 3, 0, 0, 0, 0, 0, 0],
        [1, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, 0, 3, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 3],
    ]
    if residual != expected_residual:
        raise AssertionError("2x6 residual Hom pattern changed")
    multiplicities = [7, 2, 2] + [1] * 14
    if sum(value * value for value in multiplicities) != sum(map(sum, table)):
        raise AssertionError("2x6 Wedderburn dimension count failed")
    return {
        "claim_tag": "[THEOREM]",
        "commutant_dimension_Q_equals_dimension_C": sum(map(sum, table)),
        "common_zero_dimension": sum(common_zero),
        "common_zero_matrix_block_over_C": "M_7(C)",
        "residual_Hom_table_after_M7_rows_target_columns_source": residual,
        "wedderburn_type_over_C": "M_7(C) direct-sum M_2(C)^2 direct-sum C^14",
        "centre_dimension_over_C": 17,
        "multiplicities_of_irreducible_C_modules": multiplicities,
        "structure_checks": {
            "common_zero_subspace_is_exactly_computed": True,
            "common_zero_subspace_is_proper_in_every_character_sector": True,
            "all_nonzero_Hom_dimensions_are_exact_rational_nullities": True,
            "sum_of_multiplicity_squares": sum(value * value for value in multiplicities),
        },
        "type_derivation": (
            "The exact common-zero reducing space W has seven copies of the common "
            "(A,B)=(0,0) line and contributes End(W)=M_7(C).  Its orthogonal complement "
            "has the displayed residual Hom table.  A finite-dimensional complex *-algebra "
            "of dimension two or three has only scalar simple factors.  Thus the "
            "[000,010] residual subtable [[2,1],[1,2]] is M_2(C) plus C^2, the "
            "[100,110] subtable [[1,1],[1,1]] is M_2(C), and the four isolated "
            "dimension-three diagonal blocks give C^12."
        ),
    }


def hom_closure(columns: int, *, progress: bool = False) -> dict[str, object]:
    """Certify every C2^3 Hom block over Q by one modular rank and exact lifts."""
    started = time.monotonic()
    sectors = [symmetry_sector(columns, character) for character in CHARACTERS]
    common_zero, common_zero_digests = common_zero_dimensions(sectors)
    table = [[0] * len(sectors) for _ in sectors]
    records: list[dict[str, object]] = []
    for target_index, target in enumerate(sectors):
        for source_index in range(target_index, len(sectors)):
            source = sectors[source_index]
            variables, variable_index = block_variables(target, source)
            block_started = time.monotonic()
            if progress:
                print(
                    f"BEGIN {LABELS[target_index]}<-{LABELS[source_index]} "
                    f"variables={len(variables)}",
                    flush=True,
                )
            rank_primary, modular_basis, stats = modular_nullspace(
                target, source, variable_index, len(variables), PRIMARY_PRIME
            )
            after_primary = time.monotonic()
            nullity = len(modular_basis)
            if rank_primary + nullity != len(variables):
                raise AssertionError("primary modular rank-nullity mismatch")
            if nullity:
                rank_secondary, secondary_basis, secondary_stats = modular_nullspace(
                    target, source, variable_index, len(variables), SECONDARY_PRIME
                )
                after_secondary = time.monotonic()
                exact_basis, anchors, primary_determinant, secondary_determinant = (
                    reconstruct_two_prime_basis(modular_basis, secondary_basis)
                )
                exact_kernel = verify_exact_kernel(target, source, variable_index, exact_basis)
                if not exact_kernel:
                    raise AssertionError(
                        "two-prime rational reconstruction is not exact in block "
                        f"{LABELS[target_index]}<-{LABELS[source_index]}"
                    )
                after_exact = time.monotonic()
            else:
                rank_secondary, secondary_stats = modular_rank(
                    target, source, variable_index, len(variables), SECONDARY_PRIME
                )
                after_secondary = time.monotonic()
                exact_basis = []
                anchors = []
                primary_determinant = 1
                secondary_determinant = 1
                after_exact = after_secondary
            if rank_primary != rank_secondary:
                raise AssertionError(
                    f"prime ranks differ in block {LABELS[target_index]}<-{LABELS[source_index]}"
                )
            table[target_index][source_index] = nullity
            table[source_index][target_index] = nullity
            records.append(
                {
                    "target_character": LABELS[target_index],
                    "source_character": LABELS[source_index],
                    "variables_after_B": len(variables),
                    **stats,
                    "rank_mod_primary": rank_primary,
                    "rank_mod_secondary": rank_secondary,
                    "nullity_Q": nullity,
                    "exact_rational_kernel_vectors": nullity,
                    "kernel_lift_sha256": vector_digest(exact_basis),
                    "kernel_maximum_coefficient_bits": maximum_fraction_bits(exact_basis),
                    "independence_anchor_columns": anchors,
                    "independence_anchor_determinant_mod_primary": primary_determinant,
                    "independence_anchor_determinant_mod_secondary": secondary_determinant,
                    "exact_substitution_of_reconstructed_kernel": True,
                    "secondary_prime_nonzero_A_equations": secondary_stats["nonzero_A_equations"],
                }
            )
            if progress:
                print(
                    f"END {LABELS[target_index]}<-{LABELS[source_index]} "
                    f"nullity={nullity} primary={after_primary - block_started:.3f}s "
                    f"secondary={after_secondary - after_primary:.3f}s "
                    f"exact={after_exact - after_secondary:.3f}s",
                    flush=True,
                )
    return {
        "shape": [2, columns],
        "character_order": list(LABELS),
        "sector_dimensions": [sector.dimension for sector in sectors],
        "full_End_dimension_before_symmetry_split": (1 << (2 * columns)) ** 2,
        "largest_Hom_variables_after_B": max(
            int(record["variables_after_B"]) for record in records
        ),
        "orbit_metric_hermiticity_checked_exactly": True,
        "common_A0_B0_dimensions_by_character": common_zero,
        "common_A0_B0_basis_sha256_by_character": common_zero_digests,
        "hom_dimension_matrix_Q_rows_target_columns_source": table,
        "commutant_dimension_Q_equals_dimension_C": sum(map(sum, table)),
        "upper_triangle_Hom_certificates": records,
        "certification_method": (
            "For each integer Hom system M, rank over the primary prime field is at most "
            "rank_Q(M).  The primary and secondary modular nullspaces are normalized at the "
            "same independent coordinate anchors, combined by CRT, rationally reconstructed, "
            "and substituted exactly into every integer equation.  The reconstructed vectors "
            "are anchor-independent, so their number n-rank_Fp(M) is a rational kernel lower "
            "bound.  The two inequalities force rank_Q(M)=rank_Fp(M) block by block."
        ),
        "primary_prime": PRIMARY_PRIME,
        "secondary_prime": SECONDARY_PRIME,
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def probe(columns: int) -> int:
    closure = hom_closure(columns, progress=True)
    print(
        json.dumps(
            {
                "sector_dimensions": closure["sector_dimensions"],
                "hom_table": closure["hom_dimension_matrix_Q_rows_target_columns_source"],
                "dimension": closure["commutant_dimension_Q_equals_dimension_C"],
                "elapsed_seconds": closure["elapsed_seconds"],
            },
            sort_keys=True,
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", type=int, metavar="L", help="run one complete 2xL closure")
    args = parser.parse_args()
    if args.probe is not None:
        return probe(args.probe)

    started = time.monotonic()
    prior_path = ROOT / "results" / "integrability" / "commutant_family.json"
    prior = json.loads(prior_path.read_text(encoding="utf-8"))
    control = hom_closure(5)
    control_type = infer_control_type(control)
    theorem_closure = hom_closure(6)
    theorem_type = infer_2x6_type(theorem_closure)
    prior_control = prior["data"]["exact_finite_cases"]["2x5"]
    prior_scout = prior["data"]["modular_scout_2x6"]
    lower_bound = prior_scout["exact_Q_lower_bound"]["commutant_dimension_at_least"]
    upper_bound = prior_scout["full_commutant_nullity_over_Fp"]

    checks = [
        check(
            "control_2x5_dimension",
            control["commutant_dimension_Q_equals_dimension_C"] == 20,
            "same two-prime/exact-substitution pipeline returns 20",
        ),
        check(
            "control_2x5_prior_table",
            control["hom_dimension_matrix_Q_rows_target_columns_source"]
            == prior_control["hom_dimension_matrix_Q_rows_target_columns_source"],
            "all 36 Hom dimensions reproduce the stored exact 2x5 table",
        ),
        check(
            "control_2x5_split",
            control_type["wedderburn_type_over_C"] == prior_control["wedderburn_type_over_C"],
            "Hom-table semisimplicity deduction reproduces M_2(C)^2 direct-sum C^12",
        ),
        check(
            "two_primes_agree_every_block",
            all(
                record["rank_mod_primary"] == record["rank_mod_secondary"]
                for record in theorem_closure["upper_triangle_Hom_certificates"]
            ),
            "both exact finite-field RREF ranks agree in all 36 upper-triangle blocks",
        ),
        check(
            "all_reconstructed_kernels_substitute_exactly",
            all(
                record["exact_substitution_of_reconstructed_kernel"]
                for record in theorem_closure["upper_triangle_Hom_certificates"]
            ),
            "every reconstructed rational Hom basis annihilates its integer equations",
        ),
        check(
            "stored_modular_scout_reproduced",
            theorem_closure["hom_dimension_matrix_Q_rows_target_columns_source"]
            == prior_scout["hom_dimension_matrix_over_Fp_rows_target_columns_source"],
            "exact-Q Hom table equals the prior F_1000000007 scout table",
        ),
        check(
            "common_zero_lower_bound_reproduced",
            theorem_type["common_zero_dimension"] == 7
            and theorem_type["common_zero_dimension"] ** 2 + 8 == lower_bound,
            "exact W has dimension 7 and gives the recorded lower bound 57",
        ),
        check(
            "gap_closed_against_stored_bounds",
            lower_bound
            <= theorem_type["commutant_dimension_Q_equals_dimension_C"]
            == upper_bound,
            f"stored interval [{lower_bound},{upper_bound}] closes at {upper_bound}",
        ),
        check(
            "wedderburn_dimension_count",
            theorem_type["structure_checks"]["sum_of_multiplicity_squares"]
            == theorem_type["commutant_dimension_Q_equals_dimension_C"],
            "49 + 4 + 4 + 14 = 71",
        ),
    ]
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "certifying_arithmetic": (
                "Python integers; fractions.Fraction; exact SymPy GF(p) RREF; "
                "two-prime CRT and exact rational reconstruction"
            ),
            "method": (
                "C2^3 orbit-character split; B-energy Hom variables; exact finite-field "
                "RREF lower ranks; two-prime canonical nullspace reconstruction; exact "
                "integer-equation substitution; reducing common-zero block and Hom-table "
                "Wedderburn deduction"
            ),
            "total_elapsed_seconds": round(time.monotonic() - started, 6),
        },
        "data": {
            "control_2x5": {**control, **control_type},
            "exact_2x6": {
                **theorem_closure,
                **theorem_type,
                "prior_bounds_cross_check": {
                    "prior_exact_lower_bound": lower_bound,
                    "prior_modular_upper_bound": upper_bound,
                    "closed_exact_value": theorem_type["commutant_dimension_Q_equals_dimension_C"],
                },
            },
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: wrote {RESULT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
