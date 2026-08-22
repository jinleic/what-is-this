#!/usr/bin/env python3
"""Independent exact regression for the 2 x L commutant-family certificate.

This test does not import experiments/e60_commutant_family.py. It rebuilds
character orbit bases from the largest representative (the experiment uses the
smallest), orders variables and equations backwards, and uses greatest-column
exact finite-field pivots. It reconstructs every Hom dimension at L=3 and at
the new largest theorem case L=5; stored exact-Q bases/projectors are checked
by the experiment itself.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "commutant_family.json"
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))
PRIME = 2_147_483_647


def edges(columns: int) -> tuple[tuple[int, int], ...]:
    horizontal = tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(2)
        for column in range(columns - 1)
    )
    rungs = tuple((column, columns + column) for column in range(columns))
    return horizontal + rungs


def transform(state: int, element: tuple[int, int, int], columns: int) -> int:
    row_flip, column_flip, bit_flip = element
    image = 0
    for old in range(2 * columns):
        row, column = divmod(old, columns)
        new_row = 1 - row if row_flip else row
        new_column = columns - 1 - column if column_flip else column
        new = new_row * columns + new_column
        image |= ((state >> old) & 1) << new
    if bit_flip:
        image ^= (1 << (2 * columns)) - 1
    return image


class IndependentSector:
    def __init__(self, columns: int, character: tuple[int, int, int]) -> None:
        n = 2 * columns
        visited: set[int] = set()
        vectors: list[dict[int, int]] = []
        representatives: list[int] = []
        # Reverse traversal and maximum representatives deliberately differ
        # from the generator's canonical orbit basis.
        for state in reversed(range(1 << n)):
            if state in visited:
                continue
            orbit = {transform(state, element, columns) for element in CHARACTERS}
            visited.update(orbit)
            projected: dict[int, int] = {}
            for element in CHARACTERS:
                sign = (-1) ** (sum(x * y for x, y in zip(character, element)) & 1)
                image = transform(state, element, columns)
                projected[image] = projected.get(image, 0) + sign
            projected = {configuration: value for configuration, value in projected.items() if value}
            if not projected:
                continue
            representative = max(projected)
            scale = projected[representative]
            vector = {configuration: value // scale for configuration, value in projected.items()}
            assert all(value in (-1, 1) for value in vector.values())
            representatives.append(representative)
            vectors.append(vector)

        self.dimension = len(vectors)
        self.a_columns: list[dict[int, int]] = [dict() for _ in vectors]
        self.a_rows: list[dict[int, int]] = [dict() for _ in vectors]
        self.b_diagonal: list[int] = []
        graph_edges = edges(columns)
        for source, vector in enumerate(vectors):
            a_image: dict[int, int] = {}
            for state, coefficient in vector.items():
                for site in range(n):
                    target = state ^ (1 << site)
                    a_image[target] = a_image.get(target, 0) + coefficient
            for target, representative in enumerate(representatives):
                coefficient = a_image.get(representative, 0)
                if coefficient:
                    self.a_columns[source][target] = coefficient
                    self.a_rows[target][source] = coefficient
            representative = representatives[source]
            cut = sum(
                ((representative >> left) ^ (representative >> right)) & 1
                for left, right in graph_edges
            )
            self.b_diagonal.append(len(graph_edges) - 2 * cut)


def independent_system(
    target: IndependentSector, source: IndependentSector
) -> tuple[int, list[dict[int, int]]]:
    pairs = [
        (row, column)
        for row in reversed(range(target.dimension))
        for column in reversed(range(source.dimension))
        if target.b_diagonal[row] == source.b_diagonal[column]
    ]
    variable = {pair: index for index, pair in enumerate(pairs)}
    equations: list[dict[int, int]] = []
    for row in reversed(range(target.dimension)):
        for column in reversed(range(source.dimension)):
            equation: dict[int, int] = {}
            for middle, coefficient in source.a_columns[column].items():
                index = variable.get((row, middle))
                if index is not None:
                    equation[index] = equation.get(index, 0) + coefficient
            for middle, coefficient in target.a_rows[row].items():
                index = variable.get((middle, column))
                if index is not None:
                    equation[index] = equation.get(index, 0) - coefficient
            equation = {index: value for index, value in equation.items() if value}
            if equation:
                equations.append(equation)
    return len(pairs), equations


def reverse_rank_mod(rows: list[dict[int, int]], prime: int) -> int:
    """Independent greatest-column sparse elimination over a prime field."""
    pivots: dict[int, dict[int, int]] = {}
    for source in reversed(rows):
        row = {
            column: value % prime
            for column, value in source.items()
            if value % prime
        }
        while row:
            lead = max(row)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(row[lead], prime - 2, prime)
                pivots[lead] = {
                    column: value * inverse % prime
                    for column, value in row.items()
                }
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = (row.get(column, 0) - coefficient * value) % prime
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return len(pivots)

def rebuild_hom_table(columns: int) -> tuple[list[int], list[list[int]]]:
    sectors = [IndependentSector(columns, character) for character in CHARACTERS]
    table = [[0] * 8 for _ in range(8)]
    for target_index in range(8):
        for source_index in range(target_index, 8):
            variable_count, equations = independent_system(
                sectors[target_index], sectors[source_index]
            )
            nullity = variable_count - reverse_rank_mod(equations, PRIME)
            table[target_index][source_index] = nullity
            table[source_index][target_index] = nullity
    return [sector.dimension for sector in sectors], table


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["provenance"]["script"] == "experiments/e60_commutant_family.py"
    assert all(entry["passed"] for entry in payload["checks"])
    exact = payload["data"]["exact_finite_cases"]

    sizes_3, hom_3 = rebuild_hom_table(3)
    assert sizes_3 == exact["2x3"]["sector_dimensions"] == [14, 10, 6, 6, 6, 10, 6, 6]
    assert hom_3 == exact["2x3"]["hom_dimension_matrix_Q_rows_target_columns_source"]
    assert sum(map(sum, hom_3)) == 12
    assert hom_3 == [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
    ]

    sizes_5, hom_5 = rebuild_hom_table(5)
    assert sizes_5 == exact["2x5"]["sector_dimensions"] == [152, 136, 120, 120, 120, 136, 120, 120]
    assert hom_5 == exact["2x5"]["hom_dimension_matrix_Q_rows_target_columns_source"]
    assert sum(map(sum, hom_5)) == 20
    assert hom_5 == [
        [2, 0, 1, 0, 0, 0, 0, 0],
        [0, 3, 0, 0, 0, 0, 0, 0],
        [1, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 2],
    ]

    assert exact["2x2"]["wedderburn_type_over_C"] == "M_4(C) direct-sum M_2(C)^2 direct-sum C^3"
    assert exact["2x3"]["wedderburn_type_over_C"] == "M_2(C) direct-sum C^8"
    assert exact["2x4"]["wedderburn_type_over_C"] == "M_5(C) direct-sum M_2(C) direct-sum C^12"
    assert exact["2x5"]["wedderburn_type_over_C"] == "M_2(C)^2 direct-sum C^12"
    assert exact["2x5"]["centre_dimension_over_Q_R_C"] == 14
    for case in exact.values():
        assert case["structure_checks"]["all_minimal_projector_corners_have_dimension_one_over_Q"]
        assert case["structure_checks"]["all_primitive_projectors_exact_complete_orthogonal"]
        assert sum(
            component["multiplicity"] ** 2
            for component in case["wedderburn_components"]
        ) == case["commutant_dimension_Q_equals_dimension_C"]

    scout = payload["data"]["modular_scout_2x6"]
    assert scout["status_over_Q"] == "[UNRESOLVED]"
    assert scout["exact_Q_lower_bound"]["common_A0_B0_total_dimension"] == 7
    assert scout["exact_Q_lower_bound"]["commutant_dimension_at_least"] == 57
    assert scout["full_commutant_nullity_over_Fp"] == 71
    assert "<=" in scout["rigorous_Q_statement"]
    print("independent exact Hom dimensions:", {3: sum(map(sum, hom_3)), 5: sum(map(sum, hom_5))})
    print("PASS")


if __name__ == "__main__":
    main()
