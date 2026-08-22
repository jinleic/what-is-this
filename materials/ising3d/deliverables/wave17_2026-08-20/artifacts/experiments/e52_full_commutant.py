#!/usr/bin/env python3
"""Exact full-matrix commutants of finite Ising layer generator pairs.

For A=sum_i X_i and B=sum_{ij in E} Z_i Z_j, this script computes

    Comm(A,B) = {T in End(C^(2^n)) : [T,A]=[T,B]=0}.

The 2x3 calculation directly row-reduces the 8192 by 4096 sparse integer
commutator map in the ordered Pauli basis Q_(a,b)=X^a Z^b over Q.  Larger
layers are exactly block-diagonalized by two reflections and global spin flip,
then row-reduced over Q in every intertwiner block.  Two large prime fields are
independent upper-bound checks; rational conclusions use exact-Q elimination
and explicit exact kernel vectors, never prime agreement alone.
"""
from __future__ import annotations

import itertools
import json
import platform
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "full_commutant.json"
SCRIPT = "experiments/e52_full_commutant.py"
PRIMES = (2_147_483_647, 2_147_483_629)
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))

SparseRow = dict[int, int]
RationalRow = dict[int, Fraction]


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    vertical = tuple(
        (row * columns + column, (row + 1) * columns + column)
        for row in range(rows - 1)
        for column in range(columns)
    )
    horizontal = tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(rows)
        for column in range(columns - 1)
    )
    return vertical + horizontal


def chain_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, site + 1) for site in range(n - 1))


def character_label(character: tuple[int, int, int]) -> str:
    return "".join(str(bit) for bit in character)


def pauli_commutator_rows(
    n: int, edges: tuple[tuple[int, int], ...]
) -> tuple[list[SparseRow], dict[str, int]]:
    """Rows of (ad_A/2, ad_B/2) in Q_(a,b)=X^a Z^b coordinates."""
    dimension = 1 << (2 * n)
    mask = (1 << n) - 1
    a_rows: list[SparseRow] = [{} for _ in range(dimension)]
    b_rows: list[SparseRow] = [{} for _ in range(dimension)]
    for column in range(dimension):
        a_mask = column & mask
        b_mask = column >> n
        for site in range(n):
            if (b_mask >> site) & 1:
                target = (a_mask ^ (1 << site)) | (b_mask << n)
                a_rows[target][column] = a_rows[target].get(column, 0) + 1
        for left, right in edges:
            if ((a_mask >> left) ^ (a_mask >> right)) & 1:
                target_b = b_mask ^ (1 << left) ^ (1 << right)
                target = a_mask | (target_b << n)
                b_rows[target][column] = b_rows[target].get(column, 0) - 1
    rows = [row for row in a_rows + b_rows if row]
    return rows, {
        "rows": 2 * dimension,
        "columns": dimension,
        "nonzero_rows": len(rows),
        "nonzero_entries": sum(len(row) for row in rows),
    }


def echelon_q(rows: Iterable[dict[int, int | Fraction]]) -> tuple[dict[int, RationalRow], dict[str, int]]:
    """Sparse exact-Q row echelon form, always using the least column pivot."""
    pivots: dict[int, RationalRow] = {}
    maximum_row_support = 0
    maximum_coefficient_bits = 0
    for source in rows:
        row = {column: Fraction(value) for column, value in source.items() if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                coefficient = row[lead]
                if coefficient != 1:
                    row = {column: value / coefficient for column, value in row.items()}
                pivots[lead] = row
                maximum_row_support = max(maximum_row_support, len(row))
                for value in row.values():
                    maximum_coefficient_bits = max(
                        maximum_coefficient_bits,
                        value.numerator.bit_length(),
                        value.denominator.bit_length(),
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
        "rank": len(pivots),
        "maximum_pivot_row_support": maximum_row_support,
        "maximum_coefficient_bits": maximum_coefficient_bits,
    }


def rank_mod(
    rows: Iterable[dict[int, int | Fraction]], prime: int
) -> tuple[int, dict[str, int]]:
    """Sparse exact rank over F_p, with the same deterministic pivot order."""
    pivots: dict[int, dict[int, int]] = {}
    maximum_row_support = 0
    for source in rows:
        row: dict[int, int] = {}
        for column, value in source.items():
            if isinstance(value, Fraction):
                residue = value.numerator * pow(value.denominator, prime - 2, prime) % prime
            else:
                residue = value % prime
            if residue:
                row[column] = residue
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(row[lead], prime - 2, prime)
                if inverse != 1:
                    row = {column: value * inverse % prime for column, value in row.items()}
                pivots[lead] = row
                maximum_row_support = max(maximum_row_support, len(row))
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = (row.get(column, 0) - coefficient * value) % prime
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return len(pivots), {
        "rank": len(pivots),
        "maximum_pivot_row_support": maximum_row_support,
    }


def nullspace_q(pivots: dict[int, RationalRow], n_columns: int) -> list[RationalRow]:
    """Back-substitute a canonical basis, one vector per nonpivot column."""
    free_columns = [column for column in range(n_columns) if column not in pivots]
    descending_pivots = sorted(pivots, reverse=True)
    basis: list[RationalRow] = []
    for free in free_columns:
        vector: RationalRow = {free: Fraction(1)}
        for pivot in descending_pivots:
            row = pivots[pivot]
            value = sum(
                coefficient * vector.get(column, Fraction(0))
                for column, coefficient in row.items()
                if column != pivot
            )
            if value:
                vector[pivot] = -value
        basis.append(vector)
    return basis


def kernel_failures(
    rows: Iterable[dict[int, int | Fraction]], basis: list[RationalRow]
) -> list[dict[str, int]]:
    failures: list[dict[str, int]] = []
    row_list = list(rows)
    for vector_index, vector in enumerate(basis):
        nonzero_equations = 0
        for row in row_list:
            value = sum(Fraction(coefficient) * vector.get(column, Fraction(0)) for column, coefficient in row.items())
            nonzero_equations += int(value != 0)
        if nonzero_equations:
            failures.append({"vector": vector_index, "nonzero_equations": nonzero_equations})
    return failures


def basis_rank_q(basis: list[RationalRow]) -> int:
    pivots, _ = echelon_q(basis)
    return len(pivots)


def encode_pauli_basis(basis: list[RationalRow], n: int) -> list[dict[str, object]]:
    mask = (1 << n) - 1
    encoded: list[dict[str, object]] = []
    for index, vector in enumerate(basis):
        if index == 0:
            name = "identity"
        elif index == 1:
            name = "global_spin_flip"
        else:
            name = f"kernel_{index:02d}"
        terms = []
        for column, value in sorted(vector.items()):
            terms.append(
                [
                    column & mask,
                    column >> n,
                    value.numerator,
                    value.denominator,
                ]
            )
        encoded.append({"name": name, "terms_a_b_numerator_denominator": terms})
    return encoded


def analyze_pauli_case(
    name: str,
    n: int,
    edges: tuple[tuple[int, int], ...],
) -> tuple[dict[str, object], list[RationalRow]]:
    started = time.monotonic()
    rows, matrix = pauli_commutator_rows(n, edges)
    pivots, q_stats = echelon_q(rows)
    basis = nullspace_q(pivots, matrix["columns"])
    failures = kernel_failures(rows, basis)
    modular = []
    for prime in PRIMES:
        rank, stats = rank_mod(rows, prime)
        modular.append({"prime": prime, "rank": rank, "nullity": matrix["columns"] - rank, **stats})
    identity_ok = basis and basis[0] == {0: Fraction(1)}
    parity_column = (1 << n) - 1
    parity_ok = len(basis) >= 2 and basis[1] == {parity_column: Fraction(1)}
    record: dict[str, object] = {
        "claim_tag": "[THEOREM]",
        "name": name,
        "n": n,
        "edges": [list(edge) for edge in edges],
        "map": matrix,
        "arithmetic": "direct sparse Gaussian elimination over Q in ordered Pauli coordinates",
        "rank_Q": q_stats["rank"],
        "kernel_dimension_Q": len(basis),
        "exact_elimination_stats": q_stats,
        "modular_cross_checks": modular,
        "basis_support_sizes": [len(vector) for vector in basis],
        "basis_coefficient_set": sorted(
            {str(value) for vector in basis for value in vector.values()}
        ),
        "basis_exact_rank_Q": basis_rank_q(basis),
        "basis_vectors_commute_exactly": not failures,
        "basis_failures": failures,
        "identity_and_parity_are_first": bool(identity_ok and parity_ok),
        "exact_pauli_basis": encode_pauli_basis(basis, n),
        "elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
    }
    return record, basis


def act_configuration(
    configuration: int,
    symmetry: tuple[int, int, int],
    rows: int,
    columns: int,
) -> int:
    flip_rows, flip_columns, spin_flip = symmetry
    output = 0
    n = rows * columns
    for row in range(rows):
        for column in range(columns):
            old = row * columns + column
            new_row = rows - 1 - row if flip_rows else row
            new_column = columns - 1 - column if flip_columns else column
            new = new_row * columns + new_column
            if (configuration >> old) & 1:
                output |= 1 << new
    if spin_flip:
        output ^= (1 << n) - 1
    return output


@dataclass
class Sector:
    character: tuple[int, int, int]
    representatives: list[int]
    orbit_vectors: list[dict[int, int]]
    a: list[list[int]]
    b_diagonal: list[int]

    @property
    def dimension(self) -> int:
        return len(self.representatives)

    @property
    def orbit_norms(self) -> list[int]:
        return [sum(value * value for value in vector.values()) for vector in self.orbit_vectors]


def symmetry_sector(shape: tuple[int, int], character: tuple[int, int, int]) -> Sector:
    rows, columns = shape
    n = rows * columns
    seen: set[int] = set()
    representatives: list[int] = []
    orbit_vectors: list[dict[int, int]] = []
    for configuration in range(1 << n):
        if configuration in seen:
            continue
        seen.update(
            act_configuration(configuration, symmetry, rows, columns)
            for symmetry in CHARACTERS
        )
        values: dict[int, int] = {}
        for symmetry in CHARACTERS:
            sign = -1 if sum(x * y for x, y in zip(character, symmetry)) & 1 else 1
            image = act_configuration(configuration, symmetry, rows, columns)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(state for state, value in values.items() if value)
        scale = values[representative]
        vector = {
            state: value // scale
            for state, value in values.items()
            if value
        }
        assert all(value in (-1, 1) for value in vector.values())
        representatives.append(representative)
        orbit_vectors.append(vector)

    dimension = len(representatives)
    a = [[0] * dimension for _ in range(dimension)]
    b_diagonal: list[int] = []
    edges = grid_edges(rows, columns)
    for source, vector in enumerate(orbit_vectors):
        image: dict[int, int] = {}
        for configuration, coefficient in vector.items():
            for site in range(n):
                target = configuration ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for target, representative in enumerate(representatives):
            a[target][source] = image.get(representative, 0)
        representative = representatives[source]
        b_diagonal.append(
            sum(
                1 if ((representative >> left) & 1) == ((representative >> right) & 1) else -1
                for left, right in edges
            )
        )
    return Sector(character, representatives, orbit_vectors, a, b_diagonal)


def block_commutator_system(
    target: Sector, source: Sector
) -> tuple[list[tuple[int, int]], list[SparseRow]]:
    """Equations X A_source=A_target X after imposing X B_source=B_target X."""
    variables = [
        (row, column)
        for row in range(target.dimension)
        for column in range(source.dimension)
        if target.b_diagonal[row] == source.b_diagonal[column]
    ]
    variable_index = {pair: index for index, pair in enumerate(variables)}
    equations: list[SparseRow] = []
    for row in range(target.dimension):
        for column in range(source.dimension):
            equation: SparseRow = {}
            for middle in range(source.dimension):
                variable = variable_index.get((row, middle))
                coefficient = source.a[middle][column]
                if variable is not None and coefficient:
                    equation[variable] = equation.get(variable, 0) + coefficient
            for middle in range(target.dimension):
                variable = variable_index.get((middle, column))
                coefficient = target.a[row][middle]
                if variable is not None and coefficient:
                    equation[variable] = equation.get(variable, 0) - coefficient
            equation = {index: value for index, value in equation.items() if value}
            if equation:
                equations.append(equation)
    return variables, equations


def encode_block_basis(
    variables: list[tuple[int, int]], basis: list[RationalRow]
) -> list[list[list[int]]]:
    return [
        [
            [
                variables[column][0],
                variables[column][1],
                value.numerator,
                value.denominator,
            ]
            for column, value in sorted(vector.items())
        ]
        for vector in basis
    ]


def matrix_vector(matrix: list[list[int]], vector: list[Fraction]) -> list[Fraction]:
    return [
        sum(Fraction(value) * coefficient for value, coefficient in zip(row, vector))
        for row in matrix
    ]


def common_line_certificate(
    sector: Sector,
    coordinates: dict[int, int],
    a_eigenvalue: int,
    b_eigenvalue: int,
) -> dict[str, object]:
    vector = [Fraction(coordinates.get(index, 0)) for index in range(sector.dimension)]
    a_image = matrix_vector(sector.a, vector)
    b_image = [Fraction(value) * coefficient for value, coefficient in zip(sector.b_diagonal, vector)]
    passed = (
        a_image == [a_eigenvalue * value for value in vector]
        and b_image == [b_eigenvalue * value for value in vector]
    )
    return {
        "character": character_label(sector.character),
        "coordinates_index_coefficient": [[index, value] for index, value in sorted(coordinates.items())],
        "A_eigenvalue": a_eigenvalue,
        "B_eigenvalue": b_eigenvalue,
        "passed": passed,
    }


def analyze_rectangle(shape: tuple[int, int]) -> dict[str, object]:
    started = time.monotonic()
    sectors = {character: symmetry_sector(shape, character) for character in CHARACTERS}
    hom_matrix = [[0] * len(CHARACTERS) for _ in CHARACTERS]
    block_certificates: list[dict[str, object]] = []
    total_b_compatible = 0
    total_restricted_rank_q = 0
    total_modular_ranks = {prime: 0 for prime in PRIMES}
    all_basis_verified = True
    maximum_q_bits = 0
    maximum_q_row_support = 0

    for target_index, target_character in enumerate(CHARACTERS):
        target = sectors[target_character]
        for source_index, source_character in enumerate(CHARACTERS):
            source = sectors[source_character]
            variables, equations = block_commutator_system(target, source)
            total_b_compatible += len(variables)
            pivots, q_stats = echelon_q(equations)
            basis = nullspace_q(pivots, len(variables))
            failures = kernel_failures(equations, basis)
            independent = basis_rank_q(basis) == len(basis)
            all_basis_verified &= not failures and independent
            hom_matrix[target_index][source_index] = len(basis)
            total_restricted_rank_q += q_stats["rank"]
            maximum_q_bits = max(maximum_q_bits, q_stats["maximum_coefficient_bits"])
            maximum_q_row_support = max(maximum_q_row_support, q_stats["maximum_pivot_row_support"])
            modular_rows = []
            for prime in PRIMES:
                modular_rank, modular_stats = rank_mod(equations, prime)
                total_modular_ranks[prime] += modular_rank
                modular_rows.append(
                    {
                        "prime": prime,
                        "rank": modular_rank,
                        "nullity": len(variables) - modular_rank,
                        "maximum_pivot_row_support": modular_stats["maximum_pivot_row_support"],
                    }
                )
            if basis:
                block_certificates.append(
                    {
                        "target_character": character_label(target_character),
                        "source_character": character_label(source_character),
                        "variable_count_after_B": len(variables),
                        "equation_count_for_A": len(equations),
                        "rank_Q": q_stats["rank"],
                        "hom_dimension_Q": len(basis),
                        "basis_verified_exactly": not failures,
                        "basis_independent_over_Q": independent,
                        "modular": modular_rows,
                        "basis_entries_target_source_numerator_denominator": encode_block_basis(variables, basis),
                    }
                )

    dimension = sum(sum(row) for row in hom_matrix)
    n = shape[0] * shape[1]
    end_dimension = 1 << (2 * n)
    sector_dimensions = [sectors[character].dimension for character in CHARACTERS]
    modular_nullities = {
        str(prime): total_b_compatible - total_modular_ranks[prime]
        for prime in PRIMES
    }
    record: dict[str, object] = {
        "claim_tag": "[THEOREM]",
        "shape": list(shape),
        "n": n,
        "edges": [list(edge) for edge in grid_edges(*shape)],
        "original_joint_map_shape": [2 * end_dimension, end_dimension],
        "full_joint_map_rank_Q": end_dimension - dimension,
        "kernel_dimension_Q": dimension,
        "method": (
            "exact rational C2^3 character decomposition, exact B-eigenspace restriction, "
            "and sparse Fraction row reduction of every A-intertwiner block"
        ),
        "character_order": [character_label(character) for character in CHARACTERS],
        "sector_dimensions": sector_dimensions,
        "sector_basis": {
            character_label(character): {
                "representatives": sectors[character].representatives,
                "orbit_norms": sectors[character].orbit_norms,
            }
            for character in CHARACTERS
        },
        "sector_basis_formula": (
            "For representative r and character chi, sum_g (-1)^(chi dot g)|g r>, "
            "divided by its coefficient at r; g=(row reflection,column reflection,global spin flip)."
        ),
        "B_commutant_dimension_before_A": total_b_compatible,
        "restricted_A_map_rank_Q": total_restricted_rank_q,
        "hom_dimension_matrix_Q_rows_target_columns_source": hom_matrix,
        "exact_block_basis": block_certificates,
        "exact_basis_count": sum(
            certificate["hom_dimension_Q"] for certificate in block_certificates
        ),
        "every_basis_element_verified_exactly": all_basis_verified,
        "modular_cross_checks": [
            {
                "prime": prime,
                "restricted_rank": total_modular_ranks[prime],
                "kernel_dimension_upper_bound_over_Q": modular_nullities[str(prime)],
                "interpretation": (
                    "rank over F_p is a lower bound for rank over Q; this prime alone therefore "
                    "gives the displayed upper bound on rational nullity"
                ),
            }
            for prime in PRIMES
        ],
        "exact_elimination_stats": {
            "maximum_coefficient_bits": maximum_q_bits,
            "maximum_pivot_row_support": maximum_q_row_support,
        },
        "elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
    }

    if shape == (2, 3):
        line_001 = common_line_certificate(sectors[(0, 0, 1)], {7: 1, 4: -1}, 0, -1)
        line_101 = common_line_certificate(sectors[(1, 0, 1)], {6: 1, 1: -1}, 0, 1)
        expected_nonzero = {
            (0, 0): 1,
            (1, 1): 2,
            (2, 2): 1,
            (3, 3): 1,
            (4, 4): 1,
            (4, 6): 1,
            (5, 5): 2,
            (6, 4): 1,
            (6, 6): 1,
            (7, 7): 1,
        }
        actual_nonzero = {
            (row, column): value
            for row, values in enumerate(hom_matrix)
            for column, value in enumerate(values)
            if value
        }
        structure_ok = actual_nonzero == expected_nonzero and line_001["passed"] and line_101["passed"]
        record["commutant_algebra_structure"] = {
            "claim_tag": "[THEOREM]",
            "over_C": "M_2(C) direct-sum C^8",
            "dimension": 12,
            "centre_dimension": 9,
            "geometric_group_algebra": "C[C2(row reflection) x C2(column reflection) x C2(global spin flip)]",
            "geometric_group_algebra_dimension": 8,
            "dimensions_beyond_geometric_group_algebra": 4,
            "equivalent_sector_pair_forming_M2": ["100", "110"],
            "split_sectors": [line_001, line_101],
            "justification": (
                "Diagonal Hom dimensions one certify irreducibility; diagonal dimension two in "
                "001 and 101 is split by the displayed common eigenline; the sole off-diagonal "
                "Hom pair 100<->110 gives two equivalent irreducibles and hence an M2 block."
            ),
            "structure_check_passed": structure_ok,
        }
    return record


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    total_started = time.monotonic()

    chain_records: list[dict[str, object]] = []
    chain_expected = {4: 5, 5: 6, 6: 7}
    for n in (4, 5, 6):
        record, _ = analyze_pauli_case(f"open_chain_{n}", n, chain_edges(n))
        record["control_interpretation"] = {
            "claim_tag": "[COMPUTATION]",
            "free_fermion_number_sector_count": n + 1,
            "matches_exact_commutant_dimension": record["kernel_dimension_Q"] == n + 1,
            "scope": "finite n only; the dimensions match the n+1 eigenspaces of the conserved complex-structure number operator",
        }
        chain_records.append(record)

    cycle_record, _ = analyze_pauli_case("2x2_open_grid_is_C4", 4, grid_edges(2, 2))
    direct_23, basis_23 = analyze_pauli_case("2x3_open_grid", 6, grid_edges(2, 3))

    rectangles = {
        "2x2": analyze_rectangle((2, 2)),
        "2x3": analyze_rectangle((2, 3)),
        "2x4": analyze_rectangle((2, 4)),
        "3x3": analyze_rectangle((3, 3)),
    }

    dimensions = {
        "chain_4": chain_records[0]["kernel_dimension_Q"],
        "chain_5": chain_records[1]["kernel_dimension_Q"],
        "chain_6": chain_records[2]["kernel_dimension_Q"],
        "2x2": rectangles["2x2"]["kernel_dimension_Q"],
        "2x3": direct_23["kernel_dimension_Q"],
        "2x4": rectangles["2x4"]["kernel_dimension_Q"],
        "3x3": rectangles["3x3"]["kernel_dimension_Q"],
    }

    direct_mod_ok = all(row["nullity"] == 12 for row in direct_23["modular_cross_checks"])
    rectangle_mod_ok = all(
        all(
            row["kernel_dimension_upper_bound_over_Q"] == record["kernel_dimension_Q"]
            for row in record["modular_cross_checks"]
        )
        for record in rectangles.values()
    )
    checks = [
        check(
            "2x3_direct_exact_Q_rank",
            direct_23["map"] == {
                "rows": 8192,
                "columns": 4096,
                "nonzero_rows": 8000,
                "nonzero_entries": 26624,
            }
            and direct_23["rank_Q"] == 4084
            and direct_23["kernel_dimension_Q"] == 12,
            "direct 8192x4096 ordered-Pauli sparse map has exact-Q rank 4084 and nullity 12",
        ),
        check(
            "2x3_explicit_basis",
            len(basis_23) == 12
            and direct_23["basis_exact_rank_Q"] == 12
            and direct_23["basis_vectors_commute_exactly"]
            and direct_23["identity_and_parity_are_first"],
            "all 12 stored ordered-Pauli vectors are exactly independent and commute term-by-term",
        ),
        check(
            "2x3_independent_block_agreement",
            rectangles["2x3"]["kernel_dimension_Q"] == 12
            and rectangles["2x3"]["full_joint_map_rank_Q"] == 4084,
            "independent computational-basis C2^3 block reduction also gives exact dimension 12",
        ),
        check(
            "2x3_dimension_two_conjecture_falsified",
            direct_23["kernel_dimension_Q"] == 12
            and rectangles["2x3"]["commutant_algebra_structure"]["geometric_group_algebra_dimension"] == 8,
            "row reflection, column reflection, and parity already span 8 dimensions; full dimension is 12",
        ),
        check(
            "2x3_commutant_algebra_structure",
            rectangles["2x3"]["commutant_algebra_structure"]["structure_check_passed"],
            "exact Hom table and common eigenlines certify M2(C) direct-sum C^8",
        ),
        check(
            "chain_free_fermion_controls",
            all(
                record["kernel_dimension_Q"] == chain_expected[record["n"]]
                and record["basis_vectors_commute_exactly"]
                for record in chain_records
            ),
            "open-chain n=4,5,6 dimensions are exactly 5,6,7=n+1",
        ),
        check(
            "2x2_large_commutant_verified",
            cycle_record["kernel_dimension_Q"] == 27
            and cycle_record["basis_vectors_commute_exactly"]
            and rectangles["2x2"]["kernel_dimension_Q"] == 27
            and rectangles["2x2"]["every_basis_element_verified_exactly"],
            "direct Pauli and independent symmetry-block exact-Q calculations both give 27",
        ),
        check(
            "2x4_exact_full_commutant",
            rectangles["2x4"]["kernel_dimension_Q"] == 41
            and rectangles["2x4"]["exact_basis_count"] == 41
            and rectangles["2x4"]["every_basis_element_verified_exactly"],
            "41 exact-Q block-basis elements plus exact rank certify the 65536-variable kernel",
        ),
        check(
            "3x3_exact_full_commutant",
            rectangles["3x3"]["kernel_dimension_Q"] == 98
            and rectangles["3x3"]["exact_basis_count"] == 98
            and rectangles["3x3"]["every_basis_element_verified_exactly"],
            "98 exact-Q block-basis elements plus exact rank certify the 262144-variable kernel",
        ),
        check(
            "modular_upper_bound_cross_checks",
            direct_mod_ok and rectangle_mod_ok,
            "both listed primes reproduce every exact nullity; modular ranks are used only as rational upper bounds on nullity",
        ),
    ]

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python_version": platform.python_version(),
            "method": (
                "exact sparse integer/rational commutator kernels; direct ordered-Pauli Q elimination "
                "for 2x3; exact C2^3 symmetry and B-energy blocks for larger layers; two prime-field rank cross-checks"
            ),
            "certifying_arithmetic": "integers and fractions.Fraction only; elapsed times are noncertifying",
            "total_elapsed_seconds_noncertifying": round(time.monotonic() - total_started, 6),
        },
        "data": {
            "definition": {
                "claim_tag": "[LEMMA]",
                "commutant": "Comm(A,B)={T in End(C^(2^n)) : [T,A]=[T,B]=0}",
                "A": "sum_v X_v",
                "B": "sum_{(u,v) in E} Z_u Z_v",
                "ordered_pauli_basis": "Q_(a,b)=X^a Z^b, column index a+(b<<n)",
                "associative_algebra_identity": "Comm(A,B)=Comm(<A,B>_assoc)",
                "all_couplings_identity": "Comm(A,B)=intersection_lambda Comm(A+lambda B)",
            },
            "finite_exact_dimensions": {
                "claim_tag": "[THEOREM]",
                **dimensions,
                "scope": "exactly the listed finite open graphs; no all-size grid theorem is asserted",
            },
            "two_by_three_direct_pauli": direct_23,
            "chain_controls": chain_records,
            "two_by_two_direct_pauli": cycle_record,
            "rectangular_symmetry_block_certificates": rectangles,
            "burnside_corollary": {
                "claim_tag": "[LEMMA]",
                "conditional_statement": (
                    "If a finite Hermitian pair commuting with a central involution P had Comm(A,B)=span{I,P}, "
                    "then <A,B>_assoc would equal End(V_+) direct-sum End(V_-)."
                ),
                "reason": (
                    "The unital *-algebra equals its double commutant; the commutant of span{I,P} is exactly "
                    "the full block algebra preserving the two P eigenspaces."
                ),
                "application_here": (
                    "The antecedent is false for 2x3: its commutant has dimension 12 and structure M2(C) direct-sum C^8."
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
