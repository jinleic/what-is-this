#!/usr/bin/env python3
"""Independent regression for the exact full-commutant certificate.

This test deliberately does not import experiments/e52_full_commutant.py.  It
uses computational matrix units after diagonalizing B, processes equations in
reverse order with greatest-column pivots, and checks two prime fields.  The
stored 2x3 Pauli basis is then verified term-by-term with exact Fractions.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "full_commutant.json"
PRIMES = (2_147_483_647, 2_147_483_629)


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (row * columns + column, (row + 1) * columns + column)
        for row in range(rows - 1)
        for column in range(columns)
    ) + tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(rows)
        for column in range(columns - 1)
    )


def chain_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, site + 1) for site in range(n - 1))


def computational_commutant_system(
    n: int, edges: tuple[tuple[int, int], ...]
) -> tuple[list[tuple[int, int]], list[dict[int, int]]]:
    """Parameterize [T,B]=0, then return the exact equations [T,A]=0."""
    dimension = 1 << n
    energies = []
    for state in range(dimension):
        cut = sum(((state >> left) ^ (state >> right)) & 1 for left, right in edges)
        energies.append(len(edges) - 2 * cut)
    variables = [
        (row, column)
        for row in range(dimension)
        for column in range(dimension)
        if energies[row] == energies[column]
    ]
    index = {pair: position for position, pair in enumerate(variables)}
    equations: list[dict[int, int]] = []
    for row in range(dimension):
        for column in range(dimension):
            equation: dict[int, int] = {}
            for site in range(n):
                # (T A)_(row,column)
                variable = index.get((row, column ^ (1 << site)))
                if variable is not None:
                    equation[variable] = equation.get(variable, 0) + 1
                # -(A T)_(row,column)
                variable = index.get((row ^ (1 << site), column))
                if variable is not None:
                    equation[variable] = equation.get(variable, 0) - 1
            equation = {position: value for position, value in equation.items() if value}
            if equation:
                equations.append(equation)
    return variables, equations


def reverse_echelon_mod(
    rows: list[dict[int, int]], prime: int
) -> dict[int, dict[int, int]]:
    """Independent modular elimination: reverse rows and use greatest pivots."""
    pivots: dict[int, dict[int, int]] = {}
    for source in reversed(rows):
        row = {column: value % prime for column, value in source.items() if value % prime}
        while row:
            lead = max(row)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(row[lead], prime - 2, prime)
                row = {column: value * inverse % prime for column, value in row.items()}
                pivots[lead] = row
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = (row.get(column, 0) - coefficient * value) % prime
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return pivots


def reverse_nullspace_mod(
    pivots: dict[int, dict[int, int]], n_columns: int, prime: int
) -> list[dict[int, int]]:
    free_columns = [column for column in range(n_columns) if column not in pivots]
    basis = []
    for free in free_columns:
        vector = {free: 1}
        for pivot in sorted(pivots):
            value = sum(
                coefficient * vector.get(column, 0)
                for column, coefficient in pivots[pivot].items()
                if column != pivot
            ) % prime
            if value:
                vector[pivot] = -value % prime
        basis.append(vector)
    return basis


def verify_modular_kernel(
    equations: list[dict[int, int]], basis: list[dict[int, int]], prime: int
) -> bool:
    return all(
        sum(coefficient * vector.get(column, 0) for column, coefficient in equation.items()) % prime == 0
        for vector in basis
        for equation in equations
    )


def decode_pauli_basis(record: dict[str, object]) -> list[dict[int, Fraction]]:
    n = int(record["n"])
    basis = []
    for encoded in record["exact_pauli_basis"]:
        vector = {}
        for a_mask, b_mask, numerator, denominator in encoded["terms_a_b_numerator_denominator"]:
            vector[int(a_mask) | (int(b_mask) << n)] = Fraction(int(numerator), int(denominator))
        basis.append(vector)
    return basis


def exact_pauli_commutators(
    vector: dict[int, Fraction], n: int, edges: tuple[tuple[int, int], ...]
) -> tuple[dict[int, Fraction], dict[int, Fraction]]:
    mask = (1 << n) - 1
    image_a: dict[int, Fraction] = {}
    image_b: dict[int, Fraction] = {}
    for column, coefficient in vector.items():
        a_mask = column & mask
        b_mask = column >> n
        for site in range(n):
            if (b_mask >> site) & 1:
                target = (a_mask ^ (1 << site)) | (b_mask << n)
                image_a[target] = image_a.get(target, Fraction(0)) + coefficient
        for left, right in edges:
            if ((a_mask >> left) ^ (a_mask >> right)) & 1:
                target_b = b_mask ^ (1 << left) ^ (1 << right)
                target = a_mask | (target_b << n)
                image_b[target] = image_b.get(target, Fraction(0)) - coefficient
    return (
        {column: value for column, value in image_a.items() if value},
        {column: value for column, value in image_b.items() if value},
    )


def sparse_vector_rank_mod(vectors: list[dict[int, Fraction]], prime: int) -> int:
    rows: list[dict[int, int]] = []
    for vector in vectors:
        row = {
            column: value.numerator * pow(value.denominator, prime - 2, prime) % prime
            for column, value in vector.items()
        }
        rows.append(row)
    return len(reverse_echelon_mod(rows, prime))


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert all(check["passed"] for check in payload["checks"])
    assert payload["provenance"]["script"] == "experiments/e52_full_commutant.py"

    # Independent basis/order: B-diagonal computational matrix units, not Paulis.
    variables_23, equations_23 = computational_commutant_system(6, grid_edges(2, 3))
    assert len(variables_23) == 944
    assert len(equations_23) == 3776
    modular_bases = []
    for prime in PRIMES:
        pivots = reverse_echelon_mod(equations_23, prime)
        basis = reverse_nullspace_mod(pivots, len(variables_23), prime)
        assert len(pivots) == 932
        assert len(basis) == 12
        assert verify_modular_kernel(equations_23, basis, prime)
        modular_bases.append(basis)

    # Exact verification of every stored ordered-Pauli basis element.
    direct = payload["data"]["two_by_three_direct_pauli"]
    stored_basis = decode_pauli_basis(direct)
    assert len(stored_basis) == 12
    assert stored_basis[0] == {0: Fraction(1)}
    assert stored_basis[1] == {(1 << 6) - 1: Fraction(1)}
    for vector in stored_basis:
        image_a, image_b = exact_pauli_commutators(vector, 6, grid_edges(2, 3))
        assert image_a == {}
        assert image_b == {}
    assert sparse_vector_rank_mod(stored_basis, PRIMES[0]) == 12
    assert {value for vector in stored_basis for value in vector.values()} == {Fraction(-1), Fraction(1)}

    # The one-dimensional free-fermion controls are independently rebuilt in
    # computational matrix units and compared with the stored exact values.
    stored_chains = {
        int(record["n"]): int(record["kernel_dimension_Q"])
        for record in payload["data"]["chain_controls"]
    }
    observed_chains = {}
    for n in (4, 5, 6):
        variables, equations = computational_commutant_system(n, chain_edges(n))
        nullities = []
        for prime in PRIMES:
            pivots = reverse_echelon_mod(equations, prime)
            nullities.append(len(variables) - len(pivots))
        assert nullities == [n + 1, n + 1]
        observed_chains[n] = nullities[0]
    assert observed_chains == stored_chains == {4: 5, 5: 6, 6: 7}

    dimensions = payload["data"]["finite_exact_dimensions"]
    assert {key: dimensions[key] for key in ("2x2", "2x3", "2x4", "3x3")} == {
        "2x2": 27,
        "2x3": 12,
        "2x4": 41,
        "3x3": 98,
    }
    assert direct["kernel_dimension_Q"] == 12
    structure = payload["data"]["rectangular_symmetry_block_certificates"]["2x3"]["commutant_algebra_structure"]
    assert structure["over_C"] == "M_2(C) direct-sum C^8"
    assert structure["centre_dimension"] == 9

    print("independent 2x3 modular nullspaces:", [len(basis) for basis in modular_bases])
    print("open-chain controls:", Counter(observed_chains.values()))
    print("PASS")


if __name__ == "__main__":
    main()
