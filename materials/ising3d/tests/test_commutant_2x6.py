#!/usr/bin/env python3
"""Clean-room exact verifier for the open-2x6 joint-commutant theorem.

This file imports neither the producer nor its helpers.  It deliberately uses
reverse orbit traversal, maximum orbit representatives, reverse Hom-variable
order, and rightmost anchor columns before rebuilding every Hom certificate.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterator, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "commutant_2x6.json"
CHARS = tuple(itertools.product((0, 1), repeat=3))
LABELS = tuple("".join(str(bit) for bit in character) for character in CHARS)
P1 = 1_000_000_007
P2 = 2_147_483_647


class IndependentSector:
    def __init__(self, columns: int, character: tuple[int, int, int]) -> None:
        sites = 2 * columns
        seen: set[int] = set()
        vectors: list[dict[int, int]] = []
        representatives: list[int] = []
        for state in reversed(range(1 << sites)):
            if state in seen:
                continue
            orbit = {transform(state, element, columns) for element in CHARS}
            seen.update(orbit)
            projected: dict[int, int] = {}
            for element in CHARS:
                phase = -1 if sum(x * y for x, y in zip(character, element)) & 1 else 1
                image = transform(state, element, columns)
                projected[image] = projected.get(image, 0) + phase
            projected = {key: value for key, value in projected.items() if value}
            if not projected:
                continue
            representative = max(projected)
            normalizer = projected[representative]
            vector = {key: value // normalizer for key, value in projected.items()}
            assert all(value in (-1, 1) for value in vector.values())
            representatives.append(representative)
            vectors.append(vector)

        self.character = character
        self.representatives = representatives
        self.dimension = len(representatives)
        self.a_columns: list[dict[int, int]] = [dict() for _ in vectors]
        self.a_rows: list[dict[int, int]] = [dict() for _ in vectors]
        self.b_diagonal: list[int] = []
        self.norms = [sum(value * value for value in vector.values()) for vector in vectors]
        graph_edges = edges(columns)
        for source, vector in enumerate(vectors):
            image: dict[int, int] = {}
            for state, coefficient in vector.items():
                for site in reversed(range(sites)):
                    flipped = state ^ (1 << site)
                    image[flipped] = image.get(flipped, 0) + coefficient
            for target, representative in enumerate(representatives):
                coefficient = image.get(representative, 0)
                if coefficient:
                    self.a_columns[source][target] = coefficient
                    self.a_rows[target][source] = coefficient
            state = representatives[source]
            cut = sum(((state >> left) ^ (state >> right)) & 1 for left, right in graph_edges)
            self.b_diagonal.append(len(graph_edges) - 2 * cut)
        assert all(
            self.norms[row] * coefficient
            == self.norms[column] * self.a_columns[row].get(column, 0)
            for column, entries in enumerate(self.a_columns)
            for row, coefficient in entries.items()
        )


def edges(columns: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(2)
        for column in range(columns - 1)
    ) + tuple((column, columns + column) for column in range(columns))


def transform(state: int, element: tuple[int, int, int], columns: int) -> int:
    row_flip, column_flip, bit_flip = element
    output = 0
    for old in range(2 * columns):
        row, column = divmod(old, columns)
        new_row = 1 - row if row_flip else row
        new_column = columns - 1 - column if column_flip else column
        output |= ((state >> old) & 1) << (new_row * columns + new_column)
    if bit_flip:
        output ^= (1 << (2 * columns)) - 1
    return output


def block_system(
    target: IndependentSector, source: IndependentSector
) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int], Iterator[dict[int, int]]]:
    variables = [
        (row, column)
        for row in reversed(range(target.dimension))
        for column in reversed(range(source.dimension))
        if target.b_diagonal[row] == source.b_diagonal[column]
    ]
    index = {pair: position for position, pair in enumerate(variables)}

    def rows() -> Iterator[dict[int, int]]:
        for row in reversed(range(target.dimension)):
            for column in reversed(range(source.dimension)):
                equation: dict[int, int] = {}
                for middle, coefficient in source.a_columns[column].items():
                    position = index.get((row, middle))
                    if position is not None:
                        equation[position] = equation.get(position, 0) + coefficient
                for middle, coefficient in target.a_rows[row].items():
                    position = index.get((middle, column))
                    if position is not None:
                        equation[position] = equation.get(position, 0) - coefficient
                equation = {position: value for position, value in equation.items() if value}
                if equation:
                    yield equation

    return variables, index, rows()


def modular_nullspace(
    rows: Iterator[dict[int, int]], columns: int, prime: int
) -> tuple[int, list[list[int]], int]:
    domain = sp.GF(prime)
    entries: dict[int, dict[int, object]] = {}
    row_count = 0
    for source in rows:
        row = {
            column: domain.convert(value % prime)
            for column, value in source.items()
            if value % prime
        }
        if row:
            entries[row_count] = row
            row_count += 1
    matrix = sp.polys.matrices.DomainMatrix(entries, (row_count, columns), domain)
    rref, pivots = matrix.rref()
    nullspace = rref.nullspace_from_rref(pivots).to_Matrix()
    basis = [
        [int(nullspace[row, column]) % prime for column in range(columns)]
        for row in range(nullspace.rows)
    ]
    return len(pivots), basis, row_count


def inverse_mod(matrix: Sequence[Sequence[int]], prime: int) -> list[list[int]]:
    size = len(matrix)
    work = [
        [entry % prime for entry in row]
        + [1 if row_index == column else 0 for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for pivot in range(size):
        source = next((row for row in range(pivot, size) if work[row][pivot]), None)
        assert source is not None
        if source != pivot:
            work[pivot], work[source] = work[source], work[pivot]
        inverse = pow(work[pivot][pivot], prime - 2, prime)
        work[pivot] = [entry * inverse % prime for entry in work[pivot]]
        for row in range(size):
            if row == pivot:
                continue
            factor = work[row][pivot]
            if factor:
                work[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(work[row], work[pivot])
                ]
    return [row[size:] for row in work]


def choose_rightmost_anchors(basis: Sequence[Sequence[int]], prime: int) -> list[int]:
    if not basis:
        return []
    dimension = len(basis)
    pivots: dict[int, list[int]] = {}
    anchors: list[int] = []
    for column in reversed(range(len(basis[0]))):
        vector = [basis[row][column] % prime for row in range(dimension)]
        while any(vector):
            lead = max(index for index, value in enumerate(vector) if value)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = [value * inverse % prime for value in vector]
                anchors.append(column)
                break
            factor = vector[lead]
            vector = [(value - factor * old_value) % prime for value, old_value in zip(vector, old)]
        if len(anchors) == dimension:
            return anchors
    raise AssertionError("modular nullspace basis is dependent")


def normalized_basis(
    basis: Sequence[Sequence[int]], anchors: Sequence[int], prime: int
) -> list[list[int]]:
    if not basis:
        return []
    inverse = inverse_mod([[row[column] for column in anchors] for row in basis], prime)
    return [
        [
            sum(inverse[row][middle] * basis[middle][column] for middle in range(len(basis)))
            % prime
            for column in range(len(basis[0]))
        ]
        for row in range(len(basis))
    ]


def crt(first: int, second: int) -> int:
    return first + P1 * ((second - first) * pow(P1 % P2, P2 - 2, P2) % P2)


def reconstruct(residue: int) -> Fraction:
    modulus = P1 * P2
    residue %= modulus
    if not residue:
        return Fraction(0)
    bound = math.isqrt(modulus // 2)
    old_remainder, remainder = modulus, residue
    old_denominator, denominator = 0, 1
    while abs(remainder) > bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_denominator, denominator = denominator, old_denominator - quotient * denominator
    assert denominator and abs(denominator) <= bound
    if denominator < 0:
        remainder, denominator = -remainder, -denominator
    value = Fraction(remainder, denominator)
    assert (value.denominator * residue - value.numerator) % modulus == 0
    return value


def exact_basis(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]]
) -> tuple[list[list[Fraction]], list[int]]:
    assert len(first) == len(second)
    if not first:
        return [], []
    anchors = choose_rightmost_anchors(first, P1)
    first_normalized = normalized_basis(first, anchors, P1)
    second_normalized = normalized_basis(second, anchors, P2)
    result = [
        [reconstruct(crt(first_normalized[row][column], second_normalized[row][column])) for column in range(len(first[0]))]
        for row in range(len(first))
    ]
    assert all(
        result[row][anchor] == (Fraction(1) if row == anchor_index else Fraction(0))
        for row in range(len(result))
        for anchor_index, anchor in enumerate(anchors)
    )
    return result, anchors


def verify_exact(
    target: IndependentSector,
    source: IndependentSector,
    index: dict[tuple[int, int], int],
    basis: Sequence[Sequence[Fraction]],
) -> None:
    for equation in block_equations(target, source, index):
        for vector in basis:
            assert sum(coefficient * vector[column] for column, coefficient in equation.items()) == 0


def block_equations(
    target: IndependentSector, source: IndependentSector, index: dict[tuple[int, int], int]
) -> Iterator[dict[int, int]]:
    for row in reversed(range(target.dimension)):
        for column in reversed(range(source.dimension)):
            equation: dict[int, int] = {}
            for middle, coefficient in source.a_columns[column].items():
                position = index.get((row, middle))
                if position is not None:
                    equation[position] = equation.get(position, 0) + coefficient
            for middle, coefficient in target.a_rows[row].items():
                position = index.get((middle, column))
                if position is not None:
                    equation[position] = equation.get(position, 0) - coefficient
            equation = {position: value for position, value in equation.items() if value}
            if equation:
                yield equation


def fraction_echelon(rows: Sequence[dict[int, int]]) -> dict[int, dict[int, Fraction]]:
    pivots: dict[int, dict[int, Fraction]] = {}
    for source in reversed(rows):
        row = {column: Fraction(value) for column, value in source.items() if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                value = row[lead]
                pivots[lead] = {column: entry / value for column, entry in row.items() if entry}
                break
            factor = row[lead]
            for column, entry in old.items():
                updated = row.get(column, Fraction(0)) - factor * entry
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
    return pivots


def common_zero(sector: IndependentSector) -> int:
    selected = [column for column, energy in enumerate(sector.b_diagonal) if energy == 0]
    variable = {column: position for position, column in enumerate(reversed(selected))}
    equations = [
        {
            variable[column]: coefficient
            for column, coefficient in row.items()
            if column in variable
        }
        for row in sector.a_rows
    ]
    equations = [row for row in equations if row]
    pivots = fraction_echelon(equations)
    free = [column for column in range(len(selected)) if column not in pivots]
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
    assert all(
        sum(coefficient * vector.get(column, Fraction(0)) for column, coefficient in row.items()) == 0
        for row in equations
        for vector in basis
    )
    return len(basis)


def independently_rebuild(columns: int) -> tuple[list[int], list[int], list[list[int]], dict[tuple[str, str], tuple[int, int, int]]]:
    sectors = [IndependentSector(columns, character) for character in CHARS]
    table = [[0] * 8 for _ in range(8)]
    certificates: dict[tuple[str, str], tuple[int, int, int]] = {}
    for target_index, target in enumerate(sectors):
        for source_index in range(target_index, len(sectors)):
            source = sectors[source_index]
            variables, index, first_rows = block_system(target, source)
            rank_first, first_basis, first_count = modular_nullspace(first_rows, len(variables), P1)
            _, _, second_rows = block_system(target, source)
            rank_second, second_basis, second_count = modular_nullspace(second_rows, len(variables), P2)
            assert first_count == second_count
            assert rank_first == rank_second
            assert rank_first + len(first_basis) == len(variables)
            rational_basis, anchors = exact_basis(first_basis, second_basis)
            assert len(rational_basis) == len(first_basis)
            assert len(anchors) == len(rational_basis)
            verify_exact(target, source, index, rational_basis)
            table[target_index][source_index] = len(rational_basis)
            table[source_index][target_index] = len(rational_basis)
            certificates[(LABELS[target_index], LABELS[source_index])] = (
                len(variables),
                rank_first,
                rank_second,
            )
    return [sector.dimension for sector in sectors], [common_zero(sector) for sector in sectors], table, certificates


def check_type_2x6(table: list[list[int]], zero: list[int]) -> None:
    assert zero == [1, 0, 2, 0, 2, 0, 2, 0]
    residual = [
        [table[row][column] - zero[row] * zero[column] for column in range(8)]
        for row in range(8)
    ]
    assert residual == [
        [2, 0, 1, 0, 0, 0, 0, 0],
        [0, 3, 0, 0, 0, 0, 0, 0],
        [1, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, 0, 3, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 3],
    ]
    assert 7**2 + 2**2 + 2**2 + 14 == 71


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["provenance"]["script"] == "experiments/e89_commutant_2x6.py"
    assert all(entry["passed"] for entry in payload["checks"])
    stored_control = payload["data"]["control_2x5"]
    stored_theorem = payload["data"]["exact_2x6"]

    sizes_5, zero_5, table_5, certificates_5 = independently_rebuild(5)
    assert sizes_5 == [152, 136, 120, 120, 120, 136, 120, 120]
    assert zero_5 == [0] * 8
    assert table_5 == stored_control["hom_dimension_matrix_Q_rows_target_columns_source"]
    assert sum(map(sum, table_5)) == 20
    assert stored_control["wedderburn_type_over_C"] == "M_2(C)^2 direct-sum C^12"

    sizes_6, zero_6, table_6, certificates_6 = independently_rebuild(6)
    assert sizes_6 == [560, 512, 496, 512, 496, 512, 496, 512]
    assert table_6 == stored_theorem["hom_dimension_matrix_Q_rows_target_columns_source"]
    assert sum(map(sum, table_6)) == 71
    check_type_2x6(table_6, zero_6)
    assert stored_theorem["claim_tag"] == "[THEOREM]"
    assert stored_theorem["commutant_dimension_Q_equals_dimension_C"] == 71
    assert stored_theorem["wedderburn_type_over_C"] == "M_7(C) direct-sum M_2(C)^2 direct-sum C^14"
    assert stored_theorem["centre_dimension_over_C"] == 17
    assert stored_theorem["prior_bounds_cross_check"] == {
        "prior_exact_lower_bound": 57,
        "prior_modular_upper_bound": 71,
        "closed_exact_value": 71,
    }

    stored_records = {
        (record["target_character"], record["source_character"]): record
        for record in stored_theorem["upper_triangle_Hom_certificates"]
    }
    for key, (variables, rank_first, rank_second) in certificates_6.items():
        record = stored_records[key]
        assert record["variables_after_B"] == variables
        assert record["rank_mod_primary"] == rank_first
        assert record["rank_mod_secondary"] == rank_second
        assert record["nullity_Q"] == variables - rank_first
        assert record["exact_substitution_of_reconstructed_kernel"]
    print("clean-room 2x5=20; 2x6=71; C-type=M_7(C)+M_2(C)^2+C^14")
    print("PASS")


if __name__ == "__main__":
    main()
