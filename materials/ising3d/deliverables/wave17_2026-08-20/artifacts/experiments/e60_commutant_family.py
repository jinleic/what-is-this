#!/usr/bin/env python3
"""Exact symmetry-module commutants for open 2 x L Ising ladders.

For A=sum_v X_v and B=sum_{uv in E} Z_u Z_v, the Hilbert space is first
split by row reflection, column reflection, and global spin flip.  Every Hom
block is then solved over Q after imposing the diagonal B equation
structurally.  Primitive invariant projectors are obtained inside the exact
commutant algebra of each character module, not from floating eigenvectors.

The exact cases are L=2,3,4,5.  L=6 is only a finite-field scout: its modular
rank is used as a lower bound on rational rank (hence an upper bound on the
rational commutant dimension), never as an equality over Q.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import random
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Iterator

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "commutant_family.json"
SCRIPT = "experiments/e60_commutant_family.py"
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))
CHARACTER_LABELS = tuple("".join(map(str, character)) for character in CHARACTERS)
PRIME_EXACT_CHECK = 2_147_483_647
PRIME_SCOUT = 1_000_000_007

SparseRow = dict[int, int]
RationalRow = dict[int, Fraction]
MatrixRows = list[dict[int, Fraction]]


def grid_edges(columns: int) -> tuple[tuple[int, int], ...]:
    return tuple((column, columns + column) for column in range(columns)) + tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(2)
        for column in range(columns - 1)
    )


def act_configuration(
    configuration: int, character_element: tuple[int, int, int], columns: int
) -> int:
    flip_rows, flip_columns, spin_flip = character_element
    output = 0
    for row in range(2):
        for column in range(columns):
            old = row * columns + column
            new_row = 1 - row if flip_rows else row
            new_column = columns - 1 - column if flip_columns else column
            new = new_row * columns + new_column
            if (configuration >> old) & 1:
                output |= 1 << new
    if spin_flip:
        output ^= (1 << (2 * columns)) - 1
    return output


@dataclass
class Sector:
    character: tuple[int, int, int]
    representatives: list[int]
    orbit_vectors: list[dict[int, int]]
    norms: list[int]
    a_columns: list[dict[int, int]]
    a_rows: list[dict[int, int]]
    b_diagonal: list[int]

    @property
    def dimension(self) -> int:
        return len(self.representatives)


def symmetry_sector(columns: int, character: tuple[int, int, int]) -> Sector:
    n = 2 * columns
    seen: set[int] = set()
    representatives: list[int] = []
    orbit_vectors: list[dict[int, int]] = []
    for configuration in range(1 << n):
        if configuration in seen:
            continue
        orbit = {
            act_configuration(configuration, element, columns)
            for element in CHARACTERS
        }
        seen.update(orbit)
        values: dict[int, int] = {}
        for element in CHARACTERS:
            sign = -1 if sum(x * y for x, y in zip(character, element)) & 1 else 1
            image = act_configuration(configuration, element, columns)
            values[image] = values.get(image, 0) + sign
        values = {state: value for state, value in values.items() if value}
        if not values:
            continue
        representative = min(values)
        scale = values[representative]
        vector = {state: value // scale for state, value in values.items()}
        if not all(value in (-1, 1) for value in vector.values()):
            raise AssertionError("orbit normalization ceased to be integral")
        representatives.append(representative)
        orbit_vectors.append(vector)

    dimension = len(representatives)
    a_columns: list[dict[int, int]] = [dict() for _ in range(dimension)]
    a_rows: list[dict[int, int]] = [dict() for _ in range(dimension)]
    b_diagonal: list[int] = []
    edges = grid_edges(columns)
    for source, vector in enumerate(orbit_vectors):
        image: dict[int, int] = {}
        for configuration, coefficient in vector.items():
            for site in range(n):
                target = configuration ^ (1 << site)
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
    norms = [sum(value * value for value in vector.values()) for vector in orbit_vectors]
    sector = Sector(
        character,
        representatives,
        orbit_vectors,
        norms,
        a_columns,
        a_rows,
        b_diagonal,
    )
    if any(
        norms[row] * coefficient
        != norms[column] * a_columns[row].get(column, 0)
        for column, entries in enumerate(a_columns)
        for row, coefficient in entries.items()
    ):
        raise AssertionError("restricted A is not self-adjoint in the orbit metric")
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
    """Yield X A_source - A_target X=0 after [X,B]=0 was imposed."""
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


def echelon_q(
    rows: Iterable[dict[int, int | Fraction]],
) -> tuple[dict[int, RationalRow], dict[str, int]]:
    pivots: dict[int, RationalRow] = {}
    row_count = 0
    maximum_bits = 0
    maximum_support = 0
    for source in rows:
        row_count += 1
        row = {column: Fraction(value) for column, value in source.items() if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                pivot = row[lead]
                row = {column: value / pivot for column, value in row.items() if value}
                pivots[lead] = row
                maximum_support = max(maximum_support, len(row))
                maximum_bits = max(
                    maximum_bits,
                    max(
                        max(abs(value.numerator).bit_length(), value.denominator.bit_length())
                        for value in row.values()
                    ),
                )
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = row.get(column, Fraction(0)) - coefficient * value
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return pivots, {
        "input_nonzero_rows": row_count,
        "rank": len(pivots),
        "maximum_coefficient_bits": maximum_bits,
        "maximum_pivot_row_support": maximum_support,
    }


def rank_mod(
    rows: Iterable[dict[int, int | Fraction]], prime: int
) -> tuple[int, dict[str, int]]:
    """Exact sparse finite-field rank via SymPy's DomainMatrix."""
    domain = sp.GF(prime)
    data: dict[int, dict[int, object]] = {}
    row_count = 0
    maximum_input_support = 0
    maximum_column = -1
    for source in rows:
        row: dict[int, int] = {}
        for column, raw in source.items():
            if isinstance(raw, int):
                residue = raw % prime
            else:
                value = Fraction(raw)
                residue = value.numerator * pow(value.denominator, prime - 2, prime) % prime
            if residue:
                row[column] = domain.convert(residue)
                maximum_column = max(maximum_column, column)
        if row:
            data[row_count] = row
            row_count += 1
            maximum_input_support = max(maximum_input_support, len(row))
    if not data:
        rank = 0
    else:
        matrix = sp.polys.matrices.DomainMatrix(
            data, (row_count, maximum_column + 1), domain
        )
        rank = int(matrix.rank())
    return rank, {
        "input_nonzero_rows": row_count,
        "maximum_input_row_support": maximum_input_support,
    }


def nullspace_q(pivots: dict[int, RationalRow], n_columns: int) -> tuple[list[RationalRow], list[int]]:
    free_columns = [column for column in range(n_columns) if column not in pivots]
    basis: list[RationalRow] = []
    for free in free_columns:
        vector: RationalRow = {free: Fraction(1)}
        for pivot in sorted(pivots, reverse=True):
            value = -sum(
                coefficient * vector.get(column, Fraction(0))
                for column, coefficient in pivots[pivot].items()
                if column != pivot
            )
            if value:
                vector[pivot] = value
        basis.append(vector)
    return basis, free_columns


def verify_kernel(
    equations: Iterable[SparseRow], basis: list[RationalRow]
) -> bool:
    return all(
        sum(Fraction(value) * vector.get(column, Fraction(0)) for column, value in equation.items()) == 0
        for equation in equations
        for vector in basis
    )


def rational_rank(vectors: Iterable[dict[int, Fraction]]) -> int:
    pivots, _ = echelon_q(vectors)
    return len(pivots)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def vector_digest(basis: list[RationalRow], variables: list[tuple[int, int]]) -> str:
    digest = hashlib.sha256()
    for basis_index, vector in enumerate(basis):
        for variable, value in sorted(vector.items()):
            row, column = variables[variable]
            digest.update(f"{basis_index}:{row}:{column}:{value.numerator}:{value.denominator};".encode())
    return digest.hexdigest()


def vector_to_rows(
    variables: list[tuple[int, int]], vector: RationalRow, n_rows: int
) -> MatrixRows:
    rows: MatrixRows = [dict() for _ in range(n_rows)]
    for variable, value in vector.items():
        row, column = variables[variable]
        if value:
            rows[row][column] = value
    return rows


def matrix_add(left: MatrixRows, right: MatrixRows, right_scale: Fraction = Fraction(1)) -> MatrixRows:
    output: MatrixRows = [dict(row) for row in left]
    for row, entries in enumerate(right):
        for column, value in entries.items():
            total = output[row].get(column, Fraction(0)) + right_scale * value
            if total:
                output[row][column] = total
            else:
                output[row].pop(column, None)
    return output


def matrix_scale(matrix: MatrixRows, scalar: Fraction) -> MatrixRows:
    if not scalar:
        return [dict() for _ in matrix]
    return [
        {column: scalar * value for column, value in row.items() if scalar * value}
        for row in matrix
    ]


def matrix_multiply(left: MatrixRows, right: MatrixRows) -> MatrixRows:
    output: MatrixRows = [dict() for _ in left]
    for row, left_entries in enumerate(left):
        result = output[row]
        for middle, left_value in left_entries.items():
            for column, right_value in right[middle].items():
                total = result.get(column, Fraction(0)) + left_value * right_value
                if total:
                    result[column] = total
                else:
                    result.pop(column, None)
    return output


def matrix_linear_combination(coefficients: list[Fraction], basis: list[MatrixRows]) -> MatrixRows:
    output: MatrixRows = [dict() for _ in basis[0]]
    for coefficient, matrix in zip(coefficients, basis):
        if coefficient:
            output = matrix_add(output, matrix, coefficient)
    return output


def matrix_equal(left: MatrixRows, right: MatrixRows) -> bool:
    return left == right


def matrix_trace(matrix: MatrixRows) -> Fraction:
    return sum((row.get(index, Fraction(0)) for index, row in enumerate(matrix)), Fraction(0))


def matrix_rank(matrix: MatrixRows) -> int:
    return rational_rank(matrix)


def matrix_digest(matrix: MatrixRows) -> str:
    digest = hashlib.sha256()
    for row, entries in enumerate(matrix):
        for column, value in sorted(entries.items()):
            digest.update(f"{row}:{column}:{value.numerator}:{value.denominator};".encode())
    return digest.hexdigest()


def matrix_self_adjoint(matrix: MatrixRows, norms: list[int]) -> bool:
    return all(
        Fraction(norms[row]) * matrix[row].get(column, Fraction(0))
        == Fraction(norms[column]) * matrix[column].get(row, Fraction(0))
        for row in range(len(matrix))
        for column in range(len(matrix))
    )


def matrix_adjoint(matrix: MatrixRows, norms: list[int]) -> MatrixRows:
    output: MatrixRows = [dict() for _ in matrix]
    for row, entries in enumerate(matrix):
        for column, value in entries.items():
            output[column][row] = value * Fraction(norms[row], norms[column])
    return output


def matrix_from_a(sector: Sector) -> MatrixRows:
    return [
        {column: Fraction(value) for column, value in entries.items()}
        for entries in sector.a_rows
    ]


def matrix_from_b(sector: Sector) -> MatrixRows:
    return [
        ({row: Fraction(value)} if value else {})
        for row, value in enumerate(sector.b_diagonal)
    ]


def identity_matrix(dimension: int) -> MatrixRows:
    return [{index: Fraction(1)} for index in range(dimension)]


def domain_nullspace_q(
    rows: Iterable[dict[int, int | Fraction]], n_columns: int
) -> tuple[int, list[RationalRow], list[int], dict[str, int]]:
    """Exact sparse QQ rank/nullspace, normalized at independent anchor columns."""
    domain = sp.QQ
    data: dict[int, dict[int, object]] = {}
    row_count = 0
    maximum_input_support = 0
    for source in rows:
        row = {
            column: domain.convert(value)
            for column, value in source.items()
            if value
        }
        if row:
            data[row_count] = row
            row_count += 1
            maximum_input_support = max(maximum_input_support, len(row))
    matrix = sp.polys.matrices.DomainMatrix(
        data, (row_count, n_columns), domain
    )
    rank = int(matrix.rank())
    nullspace_matrix = matrix.nullspace().to_Matrix()
    nullity = int(nullspace_matrix.rows)
    if rank + nullity != n_columns:
        raise AssertionError("exact QQ rank-nullity failed")
    if nullity == 0:
        return rank, [], [], {
            "input_nonzero_rows": row_count,
            "rank": rank,
            "maximum_coefficient_bits": 0,
            "maximum_basis_vector_support": 0,
        }
    _, anchor_tuple = nullspace_matrix.rref()
    anchors = list(anchor_tuple)
    if len(anchors) != nullity:
        raise AssertionError("nullspace rows are not independent")
    anchor_matrix = nullspace_matrix[:, anchors]
    normalized = anchor_matrix.inv() * nullspace_matrix
    basis: list[RationalRow] = []
    maximum_bits = 0
    maximum_support = 0
    for row_index in range(nullity):
        vector: RationalRow = {}
        for column in range(n_columns):
            raw = normalized[row_index, column]
            if raw:
                value = Fraction(int(raw.p), int(raw.q))
                vector[column] = value
                maximum_bits = max(
                    maximum_bits,
                    abs(value.numerator).bit_length(),
                    value.denominator.bit_length(),
                )
        basis.append(vector)
        maximum_support = max(maximum_support, len(vector))
    if any(
        basis[row].get(anchor, Fraction(0))
        != (Fraction(1) if row == column else Fraction(0))
        for row in range(nullity)
        for column, anchor in enumerate(anchors)
    ):
        raise AssertionError("anchor-normalized nullspace is not canonical")
    return rank, basis, anchors, {
        "input_nonzero_rows": row_count,
        "rank": rank,
        "maximum_coefficient_bits": maximum_bits,
        "maximum_basis_vector_support": maximum_support,
    }


@dataclass
class ExactBlock:
    target_index: int
    source_index: int
    variables: list[tuple[int, int]]
    variable_index: dict[tuple[int, int], int]
    rank_q: int
    basis: list[RationalRow]
    free_columns: list[int]
    matrices: list[MatrixRows]
    stats: dict[str, int]
    modular_rank: int
    modular_stats: dict[str, int]

    @property
    def hom_dimension(self) -> int:
        return len(self.basis)


def exact_block(target_index: int, source_index: int, target: Sector, source: Sector) -> ExactBlock:
    variables, variable_index = block_variables(target, source)
    rank_q, basis, free_columns, stats = domain_nullspace_q(
        iter_block_equations(target, source, variable_index), len(variables)
    )
    if rational_rank(basis) != len(basis):
        raise AssertionError("canonical nullspace basis is dependent")
    if not verify_kernel(iter_block_equations(target, source, variable_index), basis):
        raise AssertionError("canonical nullspace basis failed exact substitution")
    modular_rank, modular_stats = rank_mod(
        iter_block_equations(target, source, variable_index), PRIME_EXACT_CHECK
    )
    if modular_rank != rank_q:
        raise AssertionError("good-prime rank disagrees with exact rank")
    matrices = [vector_to_rows(variables, vector, target.dimension) for vector in basis]
    return ExactBlock(
        target_index,
        source_index,
        variables,
        variable_index,
        rank_q,
        basis,
        free_columns,
        matrices,
        stats,
        modular_rank,
        modular_stats,
    )


def coordinates_from_matrix(block: ExactBlock, matrix: MatrixRows) -> list[Fraction]:
    coordinates = []
    for free in block.free_columns:
        row, column = block.variables[free]
        coordinates.append(matrix[row].get(column, Fraction(0)))
    reconstructed = matrix_linear_combination(coordinates, block.matrices)
    if reconstructed != matrix:
        raise AssertionError("product did not close in canonical Hom basis")
    return coordinates


def algebra_multiply(
    left: list[Fraction], right: list[Fraction], block: ExactBlock
) -> list[Fraction]:
    left_matrix = matrix_linear_combination(left, block.matrices)
    right_matrix = matrix_linear_combination(right, block.matrices)
    return coordinates_from_matrix(block, matrix_multiply(left_matrix, right_matrix))


def primitive_projectors(sector: Sector, block: ExactBlock) -> tuple[list[dict[str, object]], list[MatrixRows], dict[str, object]]:
    """Split a character module into exact irreducible copies.

    A generic orbit-metric-self-adjoint element of the diagonal commutant is
    diagonalized in the small regular representation.  Its rational spectral
    projectors are accepted only when every exact corner P Comm P has
    dimension one.  This also handles repeated constituents inside one
    character sector (the 000 sector at L=4), where the commutant is
    noncommutative.
    """
    d = block.hom_dimension
    identity = identity_matrix(sector.dimension)
    identity_coordinates = coordinates_from_matrix(block, identity)

    commutative = True
    for left in range(d):
        for right in range(d):
            product = matrix_multiply(block.matrices[left], block.matrices[right])
            reverse = matrix_multiply(block.matrices[right], block.matrices[left])
            commutative &= product == reverse

    self_adjoint_coordinates = []
    for matrix in block.matrices:
        symmetrized = matrix_add(matrix, matrix_adjoint(matrix, sector.norms))
        self_adjoint_coordinates.append(coordinates_from_matrix(block, symmetrized))

    rng = random.Random(0xC011_6000 + sector.dimension * 17 + int("".join(map(str, sector.character)), 2))
    chosen_weights: list[int] | None = None
    chosen_eigenvalues: list[Fraction] | None = None
    chosen_regular_multiplicities: list[int] | None = None
    splitter_coordinates: list[Fraction] | None = None
    projector_rows: list[MatrixRows] | None = None
    for trial in range(1, 2049):
        if trial == 1:
            weights = list(range(1, d + 1))
        else:
            weights = [rng.randint(-17, 17) for _ in range(d)]
        splitter = [
            sum(
                Fraction(weights[generator]) * self_adjoint_coordinates[generator][coordinate]
                for generator in range(d)
            )
            for coordinate in range(d)
        ]
        if not any(splitter):
            continue
        splitter_matrix = matrix_linear_combination(splitter, block.matrices)
        if not matrix_self_adjoint(splitter_matrix, sector.norms):
            raise AssertionError("candidate splitter is not self-adjoint")
        regular_columns = []
        for right in range(d):
            unit = [Fraction(0)] * d
            unit[right] = Fraction(1)
            regular_columns.append(algebra_multiply(splitter, unit, block))
        regular = sp.Matrix(
            d,
            d,
            lambda row, column: sp.Rational(
                regular_columns[column][row].numerator,
                regular_columns[column][row].denominator,
            ),
        )
        eigenvalue_dict = regular.eigenvals()
        if not eigenvalue_dict or not all(value.is_Rational for value in eigenvalue_dict):
            continue
        eigenvalues = sorted(Fraction(int(value.p), int(value.q)) for value in eigenvalue_dict)
        candidate_projectors: list[MatrixRows] = []
        candidate_ok = True
        for eigenvalue in eigenvalues:
            coordinates = list(identity_coordinates)
            for other in eigenvalues:
                if other == eigenvalue:
                    continue
                factor = [
                    splitter[index] - other * identity_coordinates[index]
                    for index in range(d)
                ]
                coordinates = [
                    value / (eigenvalue - other)
                    for value in algebra_multiply(coordinates, factor, block)
                ]
            projector = matrix_linear_combination(coordinates, block.matrices)
            corner_vectors = []
            for matrix in block.matrices:
                corner = matrix_multiply(projector, matrix_multiply(matrix, projector))
                corner_coordinates = coordinates_from_matrix(block, corner)
                corner_vectors.append(
                    {index: value for index, value in enumerate(corner_coordinates) if value}
                )
            candidate_ok &= (
                matrix_multiply(projector, projector) == projector
                and matrix_self_adjoint(projector, sector.norms)
                and rational_rank(corner_vectors) == 1
            )
            candidate_projectors.append(projector)
        candidate_sum: MatrixRows = [dict() for _ in range(sector.dimension)]
        for projector in candidate_projectors:
            candidate_sum = matrix_add(candidate_sum, projector)
        candidate_ok &= candidate_sum == identity
        candidate_ok &= all(
            not any(matrix_multiply(candidate_projectors[left], candidate_projectors[right]))
            for left in range(len(candidate_projectors))
            for right in range(len(candidate_projectors))
            if left != right
        )
        if not candidate_ok:
            continue
        chosen_weights = weights
        chosen_eigenvalues = eigenvalues
        chosen_regular_multiplicities = [
            int(eigenvalue_dict[sp.Rational(value.numerator, value.denominator)])
            for value in eigenvalues
        ]
        splitter_coordinates = splitter
        projector_rows = candidate_projectors
        break
    if (
        chosen_weights is None
        or chosen_eigenvalues is None
        or chosen_regular_multiplicities is None
        or splitter_coordinates is None
        or projector_rows is None
    ):
        raise AssertionError("failed to find exact rational minimal invariant projectors")

    projector_records: list[dict[str, object]] = []
    a = matrix_from_a(sector)
    b = matrix_from_b(sector)
    a_squared = matrix_multiply(a, a)
    b_squared = matrix_multiply(b, b)
    for eigenvalue, projector in zip(chosen_eigenvalues, projector_rows):
        if matrix_multiply(projector, projector) != projector:
            raise AssertionError("spectral polynomial is not idempotent")
        if not matrix_self_adjoint(projector, sector.norms):
            raise AssertionError("primitive idempotent is not an orthogonal projector")
        rank = matrix_rank(projector)
        if matrix_trace(projector) != rank:
            raise AssertionError("projector trace and exact rank disagree")
        trace_a = matrix_trace(matrix_multiply(projector, a))
        trace_b = matrix_trace(matrix_multiply(projector, b))
        trace_a2 = matrix_trace(matrix_multiply(projector, a_squared))
        trace_b2 = matrix_trace(matrix_multiply(projector, b_squared))
        common_line = None
        if rank == 1:
            a_eigenvalue = trace_a
            b_eigenvalue = trace_b
            a_residual = matrix_add(
                matrix_multiply(a, projector), projector, -a_eigenvalue
            )
            b_residual = matrix_add(
                matrix_multiply(b, projector), projector, -b_eigenvalue
            )
            if any(a_residual) or any(b_residual):
                raise AssertionError("rank-one projector did not give a common eigenline")
            common_line = {
                "A_eigenvalue": fraction_text(a_eigenvalue),
                "B_eigenvalue": fraction_text(b_eigenvalue),
            }
        corner_vectors = []
        for matrix in block.matrices:
            corner = matrix_multiply(projector, matrix_multiply(matrix, projector))
            corner_coordinates = coordinates_from_matrix(block, corner)
            corner_vectors.append(
                {index: value for index, value in enumerate(corner_coordinates) if value}
            )
        corner_dimension = rational_rank(corner_vectors)
        if corner_dimension != 1:
            raise AssertionError("projector range is not irreducible")
        projector_records.append(
            {
                "splitter_eigenvalue": fraction_text(eigenvalue),
                "rank": rank,
                "trace_A": fraction_text(trace_a),
                "trace_B": fraction_text(trace_b),
                "trace_A2": fraction_text(trace_a2),
                "trace_B2": fraction_text(trace_b2),
                "common_eigenline": common_line,
                "projector_sha256": matrix_digest(projector),
                "exact_checks": {
                    "idempotent": True,
                    "orbit_metric_self_adjoint": True,
                    "rank_equals_trace": True,
                    "corner_endomorphism_dimension_Q": corner_dimension,
                },
            }
        )

    return projector_records, projector_rows, {
        "commutant_dimension_Q": d,
        "commutative_exactly": commutative,
        "canonical_basis_sha256": vector_digest(block.basis, block.variables),
        "splitter_weights_in_self_adjoint_generators": chosen_weights,
        "splitter_eigenvalues": [fraction_text(value) for value in chosen_eigenvalues],
        "regular_eigenvalue_multiplicities": chosen_regular_multiplicities,
        "projector_polynomial": "prod_{mu != lambda}(T-mu I)/(lambda-mu)",
        "complete_orthogonal_resolution_exactly": True,
        "reason_irreducible_over_C": (
            "Each exact rational projector P has dim_Q(P End_(A,B) P)=1. "
            "After scalar extension this corner is C, so semisimplicity makes its range irreducible over C."
        ),
    }


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def projected_hom_dimension(
    target_projector: MatrixRows,
    source_projector: MatrixRows,
    block: ExactBlock,
) -> tuple[int, int, str | None]:
    coordinate_vectors: list[RationalRow] = []
    witness: MatrixRows | None = None
    for matrix in block.matrices:
        projected = matrix_multiply(target_projector, matrix_multiply(matrix, source_projector))
        coordinates = coordinates_from_matrix(block, projected)
        vector = {index: value for index, value in enumerate(coordinates) if value}
        coordinate_vectors.append(vector)
        if witness is None and vector:
            witness = projected
    dimension = rational_rank(coordinate_vectors)
    witness_rank = matrix_rank(witness) if witness is not None else 0
    return dimension, witness_rank, matrix_digest(witness) if witness is not None else None


def site_permutation_action(configuration: int, permutation: tuple[int, ...], spin_flip: int) -> int:
    output = 0
    for old, new in enumerate(permutation):
        if (configuration >> old) & 1:
            output |= 1 << new
    if spin_flip:
        output ^= (1 << len(permutation)) - 1
    return output


def permutation_operator_rank(columns: int, full_square_group: bool) -> tuple[int, int]:
    n = 2 * columns
    edges = {tuple(sorted(edge)) for edge in grid_edges(columns)}
    if full_square_group:
        site_permutations = [
            permutation
            for permutation in itertools.permutations(range(n))
            if {
                tuple(sorted((permutation[left], permutation[right])))
                for left, right in edges
            }
            == edges
        ]
    else:
        site_permutations = []
        for row_flip, column_flip in itertools.product((0, 1), repeat=2):
            permutation = []
            for row in range(2):
                for column in range(columns):
                    new_row = 1 - row if row_flip else row
                    new_column = columns - 1 - column if column_flip else column
                    permutation.append(new_row * columns + new_column)
            site_permutations.append(tuple(permutation))
    dimension = 1 << n
    operators: list[RationalRow] = []
    for permutation in site_permutations:
        for spin_flip in (0, 1):
            operators.append(
                {
                    site_permutation_action(source, permutation, spin_flip) * dimension + source: Fraction(1)
                    for source in range(dimension)
                }
            )
    return len(site_permutations) * 2, rational_rank(operators)


def wedderburn_text(multiplicities: list[int], field: str) -> str:
    counts: dict[int, int] = {}
    for multiplicity in multiplicities:
        counts[multiplicity] = counts.get(multiplicity, 0) + 1
    pieces = []
    for multiplicity in sorted(counts, reverse=True):
        count = counts[multiplicity]
        block = field if multiplicity == 1 else f"M_{multiplicity}({field})"
        pieces.append(block if count == 1 else f"{block}^{count}")
    return " direct-sum ".join(pieces)


def analyze_exact(columns: int) -> dict[str, object]:
    started = time.monotonic()
    sectors = [symmetry_sector(columns, character) for character in CHARACTERS]
    common_zero_by_character = exact_common_zero_dimensions(sectors)
    blocks: dict[tuple[int, int], ExactBlock] = {}
    hom_table = [[0] * 8 for _ in range(8)]
    block_records = []
    for target_index in range(8):
        for source_index in range(target_index, 8):
            block = exact_block(
                target_index,
                source_index,
                sectors[target_index],
                sectors[source_index],
            )
            blocks[(target_index, source_index)] = block
            hom_table[target_index][source_index] = block.hom_dimension
            hom_table[source_index][target_index] = block.hom_dimension
            block_records.append(
                {
                    "target_character": CHARACTER_LABELS[target_index],
                    "source_character": CHARACTER_LABELS[source_index],
                    "variables_after_B": len(block.variables),
                    "nonzero_A_equations": block.stats["input_nonzero_rows"],
                    "rank_Q": block.stats["rank"],
                    "hom_dimension_Q": block.hom_dimension,
                    "canonical_basis_sha256": vector_digest(block.basis, block.variables),
                    "basis_verified_by_exact_substitution": True,
                    "maximum_coefficient_bits": block.stats["maximum_coefficient_bits"],
                    "modular_cross_check": {
                        "prime": PRIME_EXACT_CHECK,
                        "rank": block.modular_rank,
                        "interpretation": "cross-check only; the theorem uses the exact-Q row reduction and basis",
                    },
                }
            )

    projector_records_by_sector: list[list[dict[str, object]]] = []
    projectors_by_sector: list[list[MatrixRows]] = []
    sector_algebra_records = []
    nodes = []
    node_lookup: dict[tuple[int, int], int] = {}
    for sector_index, sector in enumerate(sectors):
        records, projectors, algebra_record = primitive_projectors(
            sector, blocks[(sector_index, sector_index)]
        )
        projector_records_by_sector.append(records)
        projectors_by_sector.append(projectors)
        for projector_index, record in enumerate(records):
            node_index = len(nodes)
            node_lookup[(sector_index, projector_index)] = node_index
            node = {
                "id": f"{CHARACTER_LABELS[sector_index]}:p{projector_index}",
                "character": CHARACTER_LABELS[sector_index],
                "projector_index": projector_index,
                **record,
            }
            nodes.append(node)
        sector_algebra_records.append(
            {
                "character": CHARACTER_LABELS[sector_index],
                "sector_dimension": sector.dimension,
                "orbit_representatives_sha256": hashlib.sha256(
                    ",".join(map(str, sector.representatives)).encode()
                ).hexdigest(),
                **algebra_record,
                "primitive_projectors": records,
            }
        )

    union = UnionFind(len(nodes))
    projected_records = []
    projected_sum = 0
    for target_index in range(8):
        for source_index in range(target_index, 8):
            block = blocks[(target_index, source_index)]
            block_sum = 0
            for target_projector_index, target_projector in enumerate(projectors_by_sector[target_index]):
                for source_projector_index, source_projector in enumerate(projectors_by_sector[source_index]):
                    dimension, witness_rank, witness_digest = projected_hom_dimension(
                        target_projector, source_projector, block
                    )
                    block_sum += dimension
                    if dimension:
                        if dimension != 1:
                            raise AssertionError("primitive constituent Hom dimension is not one")
                        left = node_lookup[(target_index, target_projector_index)]
                        right = node_lookup[(source_index, source_projector_index)]
                        union.union(left, right)
                        if left != right:
                            projected_records.append(
                                {
                                    "target": nodes[left]["id"],
                                    "source": nodes[right]["id"],
                                    "hom_dimension_Q": 1,
                                    "witness_rank": witness_rank,
                                    "witness_sha256": witness_digest,
                                    "intertwines_exactly": True,
                                }
                            )
            if block_sum != block.hom_dimension:
                raise AssertionError("primitive projector Hom decomposition lost dimension")
            projected_sum += block_sum if target_index == source_index else 2 * block_sum

    components_by_root: dict[int, list[int]] = {}
    for node_index in range(len(nodes)):
        components_by_root.setdefault(union.find(node_index), []).append(node_index)
    components = []
    for indices in sorted(components_by_root.values(), key=lambda values: [nodes[index]["id"] for index in values]):
        ranks = {int(nodes[index]["rank"]) for index in indices}
        if len(ranks) != 1:
            raise AssertionError("equivalent constituents have unequal dimensions")
        multiplicity = len(indices)
        components.append(
            {
                "multiplicity": multiplicity,
                "irreducible_module_dimension": ranks.pop(),
                "copies": [nodes[index]["id"] for index in indices],
                "wedderburn_block_over_Q": "Q" if multiplicity == 1 else f"M_{multiplicity}(Q)",
                "wedderburn_block_over_R": "R" if multiplicity == 1 else f"M_{multiplicity}(R)",
                "wedderburn_block_over_C": "C" if multiplicity == 1 else f"M_{multiplicity}(C)",
                "common_eigenline_component": all(nodes[index]["common_eigenline"] is not None for index in indices),
            }
        )
    multiplicities = [int(component["multiplicity"]) for component in components]
    commutant_dimension = sum(sum(row) for row in hom_table)
    if sum(value * value for value in multiplicities) != commutant_dimension:
        raise AssertionError("Wedderburn dimensions do not sum to the Hom table")
    if projected_sum != commutant_dimension:
        raise AssertionError("projected Hom dimensions do not sum to the full commutant")

    diagonal_dimension = sum(hom_table[index][index] for index in range(8))
    group_order, spatial_rank = permutation_operator_rank(columns, full_square_group=columns == 2)
    baseline_order, baseline_rank = permutation_operator_rank(columns, full_square_group=False)
    if baseline_order != 8 or baseline_rank != 8:
        raise AssertionError("represented C2^3 geometry algebra is not eight-dimensional")
    expected_group_order = 16 if columns == 2 else 8
    expected_spatial_rank = 13 if columns == 2 else 8
    if group_order != expected_group_order or spatial_rank != expected_spatial_rank:
        raise AssertionError("unexpected represented graph/spin symmetry dimension")

    return {
        "claim_tag": "[THEOREM]",
        "shape": [2, columns],
        "sites": 2 * columns,
        "character_order": list(CHARACTER_LABELS),
        "sector_dimensions": [sector.dimension for sector in sectors],
        "hom_dimension_matrix_Q_rows_target_columns_source": hom_table,
        "commutant_dimension_Q_equals_dimension_C": commutant_dimension,
        "exact_upper_triangle_block_certificates": block_records,
        "character_module_decomposition": sector_algebra_records,
        "irreducible_copies": nodes,
        "constituent_intertwiner_certificates": projected_records,
        "wedderburn_components": components,
        "wedderburn_type_over_Q": wedderburn_text(multiplicities, "Q"),
        "wedderburn_type_over_R": wedderburn_text(multiplicities, "R"),
        "wedderburn_type_over_C": wedderburn_text(multiplicities, "C"),
        "centre_dimension_over_Q_R_C": len(components),
        "structure_checks": {
            "all_minimal_projector_corners_have_dimension_one_over_Q": True,
            "all_primitive_projectors_exact_complete_orthogonal": True,
            "all_projected_Hom_dimensions_zero_or_one": True,
            "sum_multiplicity_squares": sum(value * value for value in multiplicities),
        },
        "common_A0_B0_dimensions_by_character": common_zero_by_character,
        "common_A0_B0_total_dimension": sum(common_zero_by_character),
        "geometric_versus_module_symmetries": {
            "baseline_C2_cubed_group_algebra_dimension": baseline_rank,
            "full_graph_automorphism_times_spin_group_order": group_order,
            "represented_spatial_group_algebra_dimension": spatial_rank,
            "full_commutant_dimension": commutant_dimension,
            "dimensions_beyond_baseline_geometry": commutant_dimension - baseline_rank,
            "dimensions_beyond_all_spatial_permutation_symmetries": commutant_dimension - spatial_rank,
            "within_character_dimensions_beyond_sector_identities": diagonal_dimension - baseline_rank,
            "extra_primitive_projector_directions": len(nodes) - baseline_rank,
            "intra_character_intertwiner_directions": diagonal_dimension - len(nodes),
            "off_character_intertwiner_dimensions": commutant_dimension - diagonal_dimension,
            "interpretation": (
                "The exact extra basis is exhausted by primitive module projectors and inter-character "
                "or intra-character intertwiners.  For L=2, the sixteen D4-times-spin group elements have a "
                "thirteen-dimensional represented algebra, five dimensions beyond C2^3; "
                "for L>=3 the full spatial image is already C2^3."
            ),
        },
        "common_eigenline_components": [
            component for component in components if component["common_eigenline_component"]
        ],
        "elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
    }


def exact_common_zero_dimensions(sectors: list[Sector]) -> list[int]:
    """Dimensions of ker(A) intersect ker(B) in every character sector."""
    dimensions = []
    for sector in sectors:
        zero_energy_columns = [
            column for column, energy in enumerate(sector.b_diagonal) if energy == 0
        ]
        variable = {
            column: index for index, column in enumerate(zero_energy_columns)
        }
        equations = []
        for row in sector.a_rows:
            equation = {
                variable[column]: coefficient
                for column, coefficient in row.items()
                if column in variable
            }
            if equation:
                equations.append(equation)
        pivots, _ = echelon_q(equations)
        basis, _ = nullspace_q(pivots, len(zero_energy_columns))
        if not verify_kernel(equations, basis):
            raise AssertionError("common-zero basis failed exact substitution")
        dimensions.append(len(basis))
    return dimensions


def analyze_modular_scout(columns: int) -> dict[str, object]:
    started = time.monotonic()
    sectors = [symmetry_sector(columns, character) for character in CHARACTERS]
    common_zero_by_character = exact_common_zero_dimensions(sectors)
    common_zero_dimension = sum(common_zero_by_character)
    every_character_has_complement = all(
        sector.dimension > zero_dimension
        for sector, zero_dimension in zip(sectors, common_zero_by_character)
    )
    if not every_character_has_complement:
        raise AssertionError("geometric algebra may intersect the common-zero endomorphism algebra")
    exact_lower_bound = common_zero_dimension * common_zero_dimension + 8
    table = [[0] * 8 for _ in range(8)]
    block_records = []
    total_variables = 0
    total_rank = 0
    for target_index in range(8):
        for source_index in range(target_index, 8):
            target = sectors[target_index]
            source = sectors[source_index]
            variables, variable_index = block_variables(target, source)
            rank, stats = rank_mod(
                iter_block_equations(target, source, variable_index), PRIME_SCOUT
            )
            nullity = len(variables) - rank
            table[target_index][source_index] = nullity
            table[source_index][target_index] = nullity
            factor = 1 if target_index == source_index else 2
            total_variables += factor * len(variables)
            total_rank += factor * rank
            block_records.append(
                {
                    "target_character": CHARACTER_LABELS[target_index],
                    "source_character": CHARACTER_LABELS[source_index],
                    "variables_after_B": len(variables),
                    "rank_mod_p": rank,
                    "nullity_mod_p": nullity,
                    "nonzero_A_equations": stats["input_nonzero_rows"],
                }
            )
    modular_nullity = total_variables - total_rank
    if modular_nullity < exact_lower_bound:
        raise AssertionError("modular upper bound contradicts exact rational lower bound")
    return {
        "claim_tag": "[COMPUTATION]",
        "shape": [2, columns],
        "sites": 2 * columns,
        "prime": PRIME_SCOUT,
        "sector_dimensions": [sector.dimension for sector in sectors],
        "hom_dimension_matrix_over_Fp_rows_target_columns_source": table,
        "full_commutant_nullity_over_Fp": modular_nullity,
        "rigorous_Q_statement": (
            f"rank over F_{PRIME_SCOUT} is a lower bound on rank over Q for the same integer blocks, "
            f"so dim_Q Comm(A,B) <= {modular_nullity}; equality is not asserted"
        ),
        "exact_Q_lower_bound": {
            "claim_tag": "[THEOREM]",
            "common_A0_B0_dimensions_by_character": common_zero_by_character,
            "common_A0_B0_total_dimension": common_zero_dimension,
            "endomorphism_algebra_dimension": common_zero_dimension * common_zero_dimension,
            "geometric_C2_cubed_dimension_with_zero_intersection": 8,
            "commutant_dimension_at_least": exact_lower_bound,
            "proof": (
                "A and B vanish on the exact rational common-zero space W, and Hermiticity makes W reducing, "
                "so End(W) embeds in the commutant.  Every character sector has a nonzero complement to W, "
                "hence no nonzero C2^3 character-projector combination is supported only on W; the two "
                "subspaces intersect trivially."
            ),
        },
        "status_over_Q": "[UNRESOLVED]",
        "upper_triangle_block_records": block_records,
        "elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
    }


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    started = time.monotonic()
    exact_cases = {f"2x{columns}": analyze_exact(columns) for columns in (2, 3, 4, 5)}
    scout = analyze_modular_scout(6)
    expected_dimensions = {"2x2": 27, "2x3": 12, "2x4": 41, "2x5": 20}
    expected_types = {
        "2x2": "M_4(C) direct-sum M_2(C)^2 direct-sum C^3",
        "2x3": "M_2(C) direct-sum C^8",
        "2x4": "M_5(C) direct-sum M_2(C) direct-sum C^12",
        "2x5": "M_2(C)^2 direct-sum C^12",
    }
    checks = [
        check(
            "exact_finite_dimensions",
            all(
                exact_cases[name]["commutant_dimension_Q_equals_dimension_C"] == dimension
                for name, dimension in expected_dimensions.items()
            ),
            "exact symmetry-character Hom sums are 27,12,41,20 for L=2,3,4,5",
        ),
        check(
            "two_by_three_structural_theorem",
            exact_cases["2x3"]["wedderburn_type_over_C"] == expected_types["2x3"]
            and exact_cases["2x3"]["centre_dimension_over_Q_R_C"] == 9,
            "primitive exact projectors and the unique repeated constituent give M2(C)+C^8",
        ),
        check(
            "all_exact_wedderburn_types",
            all(exact_cases[name]["wedderburn_type_over_C"] == expected for name, expected in expected_types.items()),
            "all four exact cases have the independently dimension-checked displayed Wedderburn type",
        ),
        check(
            "new_two_by_five_theorem",
            exact_cases["2x5"]["commutant_dimension_Q_equals_dimension_C"] == 20
            and exact_cases["2x5"]["wedderburn_type_over_C"] == "M_2(C)^2 direct-sum C^12"
            and exact_cases["2x5"]["centre_dimension_over_Q_R_C"] == 14,
            "2x5 has exact dimension 20 and type M2(C)^2+C^12",
        ),
        check(
            "geometry_separated_from_module_symmetry",
            all(
                case["geometric_versus_module_symmetries"]["baseline_C2_cubed_group_algebra_dimension"] == 8
                for case in exact_cases.values()
            )
            and exact_cases["2x2"]["geometric_versus_module_symmetries"]["represented_spatial_group_algebra_dimension"] == 13
            and all(
                exact_cases[f"2x{columns}"]["geometric_versus_module_symmetries"]["represented_spatial_group_algebra_dimension"] == 8
                for columns in (3, 4, 5)
            ),
            "C2^3 contributes eight dimensions; the square's sixteen D4-times-spin elements span dimension thirteen",
        ),
        check(
            "two_by_six_exact_lower_bound",
            scout["exact_Q_lower_bound"]["common_A0_B0_total_dimension"] == 7
            and scout["exact_Q_lower_bound"]["commutant_dimension_at_least"] == 57
            and scout["full_commutant_nullity_over_Fp"] == 71,
            "an exact seven-dimensional common A=B=0 space plus disjoint C2^3 geometry gives dim_Q Comm >=57",
        ),
        check(
            "modular_scout_is_only_upper_bound",
            scout["status_over_Q"] == "[UNRESOLVED]"
            and "<=" in scout["rigorous_Q_statement"],
            "the 2x6 prime-field result is explicitly retained only as a rational upper bound",
        ),
    ]
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python_version": platform.python_version(),
            "certifying_arithmetic": "integers, fractions.Fraction, exact SymPy QQ, and exact SymPy GF(p) only",
            "method": (
                "C2^3 orbit-character modules; exact sparse Q Hom kernels; exact multiplication tables; "
                "rational spectral idempotents in split character commutants; exact projected Hom graph"
            ),
            "total_elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
        },
        "data": {
            "definition": {
                "claim_tag": "[LEMMA]",
                "A": "sum_v X_v",
                "B": "sum_{uv in E(2xL open ladder)} Z_u Z_v",
                "commutant": "{T:[T,A]=[T,B]=0}",
                "character_group": "row reflection x column reflection x global spin flip = C2^3",
            },
            "exact_finite_cases": exact_cases,
            "modular_scout_2x6": scout,
            "finite_family_summary": {
                "claim_tag": "[THEOREM]",
                "dimensions_L_2_through_5": [27, 12, 41, 20],
                "wedderburn_types": {name: case["wedderburn_type_over_C"] for name, case in exact_cases.items()},
                "scope": "only L=2,3,4,5; no all-L dimension or Wedderburn formula is asserted",
            },
            "candidate_formula_audit": {
                "claim_tag": "[COMPUTATION]",
                "odd_constant_candidate_from_L3": "M2(C)+C^8 (dimension 12) for every odd L",
                "odd_constant_minimal_counterexample": "L=5: exact type M2(C)^2+C^12 and dimension 20",
                "two_point_parity_linear_fit": "4L for odd L=3,5 and 7L+13 for even L=2,4",
                "even_fit_minimal_counterexample": (
                    "At L=6 the fit predicts 55, but an exact seven-dimensional common A=B=0 "
                    "space plus the disjoint eight-dimensional geometric algebra proves dimension at least 57."
                ),
                "odd_4L_status": "[CONJECTURE] from two data points only; no induction or local recurrence",
                "all_L_status": "[UNRESOLVED]",
            },
            "monomial_symmetry_classification": {
                "claim_tag": "[LEMMA]",
                "statement": (
                    "A computational-basis monomial permutation induced by hypercube coordinate permutations "
                    "and bit translations commutes with A and B iff its site permutation is a ladder graph "
                    "automorphism and, because the ladder is connected, its bit translation is either identity "
                    "or global spin flip."
                ),
                "consequence": (
                    "For L>=3 all monomial spatial symmetries lie in the eight-dimensional C2^3 image. "
                    "The remaining exact elements are non-monomial module projectors/intertwiners."
                ),
            },
            "bipartite_antisymmetry_mechanism": {
                "claim_tag": "[LEMMA]",
                "checkerboard_X": "U_X=product of X on one ladder bipartition: U_X A U_X=A and U_X B U_X=-B",
                "global_Z": "U_Z=product_v Z_v: U_Z A U_Z=-A and U_Z B U_Z=B",
                "reflection_action": (
                    "row reflection sends U_X to P U_X; column reflection sends U_X to P U_X "
                    "for even L and fixes U_X for odd L"
                ),
                "role": (
                    "These exact all-L spectral dualities preserve the common A=B=0 space and organize "
                    "some character/eigenvalue pairings, but neither is a joint conserved operator on the full space."
                ),
            },
            "integrability_scope": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "These are finite-layer, coupling-independent module symmetries.  No compatible local-density "
                    "recurrence or induction in L is certified, so the finite commutants are not promoted to an "
                    "extensive family of local charges and do not establish integrability."
                ),
            },
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    failed = [entry["name"] for entry in checks if not entry["passed"]]
    if failed:
        print("FAIL: " + ", ".join(failed))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
