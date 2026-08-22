"""Standalone exact verifier for the one-line-cabled 16x16 RLLL calculation.

This verifier intentionally imports no producer module.  It rebuilds the seven-bit
operators and every stored fixed-q maximal minor from the JSON certificate.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "tetra16.json"
DIMENSION = 128
AUXILIARY_MASK = (1 << 4) - 1


def _bits(state: int, width: int) -> tuple[int, ...]:
    return tuple((state >> position) & 1 for position in range(width))


def _state(bits: tuple[int, ...]) -> int:
    return sum(int(bit) << position for position, bit in enumerate(bits))


def _replace_local_bits(global_state: int, positions: tuple[int, ...], local_state: int) -> int:
    result = global_state
    for local_position, global_position in enumerate(positions):
        if (local_state >> local_position) & 1:
            result |= 1 << global_position
        else:
            result &= ~(1 << global_position)
    return result


def _integer_local_matrix(q_numerator: int, q_denominator: int) -> sp.SparseMatrix:
    """Return denominator^3 times the exact binary Ising L(q)."""

    entries: dict[tuple[int, int], int] = {}
    for input_state in range(8):
        for output_state in range(8):
            exponent = input_state.bit_count() + output_state.bit_count()
            if exponent % 2 == 0:
                half_exponent = exponent // 2
                entries[(output_state, input_state)] = (
                    q_numerator**half_exponent
                    * q_denominator ** (3 - half_exponent)
                )
    return sp.SparseMatrix(8, 8, entries)


def _embedded_three_space(local: sp.SparseMatrix, positions: tuple[int, ...]) -> sp.SparseMatrix:
    entries: dict[tuple[int, int], int] = {}
    for global_input in range(DIMENSION):
        local_input = _state(tuple((global_input >> position) & 1 for position in positions))
        for local_output in range(8):
            value = local[local_output, local_input]
            if value:
                entries[(_replace_local_bits(global_input, positions, local_output), global_input)] = int(value)
    return sp.SparseMatrix(DIMENSION, DIMENSION, entries)


def _products(q_numerator: int, q_denominator: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Build the scaled cabled forward/reverse products independently."""

    local = _integer_local_matrix(q_numerator, q_denominator)
    first = _embedded_three_space(local, (0, 4, 5))
    second = _embedded_three_space(local, (1, 4, 5))
    cable = first * second
    reverse_cable = second * first
    middle = _embedded_three_space(local, (2, 4, 6))
    last = _embedded_three_space(local, (3, 5, 6))
    return cable * middle * last, last * middle * cable, cable, reverse_cable


def _coordinates() -> tuple[tuple[int, int], ...]:
    preserving = []
    reversing = []
    for input_state in range(16):
        for output_state in range(16):
            target = preserving if (output_state.bit_count() - input_state.bit_count()) % 2 == 0 else reversing
            target.append((output_state, input_state))
    coordinates = tuple(preserving + reversing)
    assert len(coordinates) == 256 and len(set(coordinates)) == 256
    return coordinates


def _component_row(
    output_state: int,
    input_state: int,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    coordinates: tuple[tuple[int, int], ...],
) -> list[int]:
    output_auxiliary = output_state & AUXILIARY_MASK
    input_auxiliary = input_state & AUXILIARY_MASK
    row: list[int] = []
    for coordinate_output, coordinate_input in coordinates:
        value = 0
        if output_auxiliary == coordinate_output:
            intermediate = (output_state & ~AUXILIARY_MASK) | coordinate_input
            value += int(forward[intermediate, input_state])
        if input_auxiliary == coordinate_input:
            intermediate = (input_state & ~AUXILIARY_MASK) | coordinate_output
            value -= int(reverse[output_state, intermediate])
        row.append(value)
    return row


def _factor_value(record: dict) -> int:
    value = int(record["sign"])
    for prime, exponent in record["prime_factors"]:
        prime = int(prime)
        exponent = int(exponent)
        assert sp.isprime(prime), f"non-prime purported factor {prime}"
        assert exponent > 0
        value *= prime**exponent
    return value


def _minor(
    rows: list[int],
    columns: range,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    coordinates: tuple[tuple[int, int], ...],
) -> sp.Matrix:
    return sp.Matrix(
        [
            _component_row(row // DIMENSION, row % DIMENSION, forward, reverse, coordinates)[columns.start : columns.stop]
            for row in rows
        ]
    )


def _modular_rank(rows: list[list[int]], prime: int) -> int:
    """Independent dense-in-the-selected-minor elimination over F_p."""

    matrix = [[int(value) % prime for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    rank = 0
    for column in range(column_count):
        pivot = next((index for index in range(rank, row_count) if matrix[index][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(row_count):
            if index == rank:
                continue
            coefficient = matrix[index][column]
            if coefficient:
                matrix[index] = [
                    (value - coefficient * pivot_value) % prime
                    for value, pivot_value in zip(matrix[index], matrix[rank], strict=True)
                ]
        rank += 1
        if rank == row_count:
            break
    return rank


def _run_check(name: str, function) -> None:
    try:
        detail = function()
    except Exception as exc:
        print(f"{name}: FAIL ({exc})")
        raise
    print(f"{name}: PASS ({detail})")


def main() -> None:
    artifact = json.loads(ARTIFACT.read_text())
    data = artifact["data"]
    coordinates = _coordinates()
    cache: dict[str, object] = {}

    def check_envelope_and_geometry() -> str:
        assert artifact["meta"]["provenance"] == "experiments/e83_tetra16.py"
        ansatz = data["ansatz"]
        assert ansatz["true_system_shape"] == [16384, 256]
        assert ansatz["global_operator_shape"] == [128, 128]
        assert ansatz["auxiliary_R_shape"] == [16, 16]
        assert ansatz["independent_unknowns"] == 256
        assert ansatz["grading_assumption_on_R"] == "none"
        cabling = data["cabling"]
        assert cabling["bit_order"] == ["i1", "i2", "2", "3", "4", "5", "6"]
        assert cabling["factor_order"] == "L_(i1,4,5) * L_(i2,4,5)"
        assert cabling["uncabled_C4_leg_status"] == "out_of_scope_different_bilinear_problem"
        assert data["sector_decomposition"]["column_counts"] == {
            "parity_preserving": 128,
            "parity_reversing": 128,
        }
        assert len(data["fixed_q_certificates"]) >= 3
        assert all(check["passed"] for check in artifact["checks"])
        return "128-dimensional seven-bit cabled geometry and >=3 exact fixed-q certificates"

    def check_cable_order_and_sector_support() -> str:
        first = data["fixed_q_certificates"][0]
        numerator, denominator = map(int, first["q"].split("/")) if "/" in first["q"] else (int(first["q"]), 1)
        forward, reverse, cable, reverse_cable = _products(numerator, denominator)
        assert cable != reverse_cable, "the two shared-physical-leg factors must make cabling order material"
        off_sector_nonzeros = {"parity_preserving_rows_to_reversing_columns": 0, "parity_reversing_rows_to_preserving_columns": 0}
        nonzero_rows = 0
        for output_state in range(DIMENSION):
            for input_state in range(DIMENSION):
                row = _component_row(output_state, input_state, forward, reverse, coordinates)
                if any(row):
                    nonzero_rows += 1
                row_grade = (output_state.bit_count() - input_state.bit_count()) % 2
                for column, value in enumerate(row):
                    if not value:
                        continue
                    column_grade = (coordinates[column][0].bit_count() - coordinates[column][1].bit_count()) % 2
                    assert row_grade == column_grade
                    if row_grade == 0 and column >= 128:
                        off_sector_nonzeros["parity_preserving_rows_to_reversing_columns"] += 1
                    if row_grade == 1 and column < 128:
                        off_sector_nonzeros["parity_reversing_rows_to_preserving_columns"] += 1
        assert nonzero_rows == 16384
        assert off_sector_nonzeros == {
            "parity_preserving_rows_to_reversing_columns": 0,
            "parity_reversing_rows_to_preserving_columns": 0,
        }
        cache["first_products"] = (forward, reverse)
        return "fixed-order cable is noncommuting; all 16384 rows respect the exact two-sector split"

    def check_fixed_q_maximal_minors() -> str:
        details = []
        for fixed in data["fixed_q_certificates"]:
            numerator, denominator = map(int, fixed["q"].split("/")) if "/" in fixed["q"] else (int(fixed["q"]), 1)
            forward, reverse, _, _ = _products(numerator, denominator)
            sector_determinants = []
            for name, columns in (("parity_preserving", range(0, 128)), ("parity_reversing", range(128, 256))):
                certificate = fixed["sectors"][name]
                rows = [int(row) for row in certificate["rows"]]
                assert len(rows) == 128 and len(set(rows)) == 128
                minor = _minor(rows, columns, forward, reverse, coordinates)
                determinant = int(minor.det(method="domain-ge"))
                assert determinant != 0
                factor = certificate["integer_factorization"]
                assert determinant == int(factor["value"])
                assert _factor_value(factor) == determinant
                assert certificate["rank"] == 128
                sector_determinants.append(determinant)
            assert int(fixed["full_block_diagonal_determinant"]) == sector_determinants[0] * sector_determinants[1]
            assert fixed["exact_full_rank"] == 256
            details.append(f"q={fixed['q']}")
        return "nonzero 128x128 exact sector minors at " + ", ".join(details)

    def check_modular_generic_evidence() -> str:
        evidence = data["modular_generic_evidence"]
        assert len(evidence) >= 2
        details = []
        for row in evidence:
            prime = int(row["prime"])
            residue = int(row["q_residue"])
            assert sp.isprime(prime)
            local = _integer_local_matrix(residue, 1)
            first = _embedded_three_space(local, (0, 4, 5))
            second = _embedded_three_space(local, (1, 4, 5))
            cable = first * second
            middle = _embedded_three_space(local, (2, 4, 6))
            last = _embedded_three_space(local, (3, 5, 6))
            forward = cable * middle * last
            reverse = last * middle * cable
            for name, columns in (("parity_preserving", range(0, 128)), ("parity_reversing", range(128, 256))):
                rows = [int(value) for value in row["sector_rows"][name]]
                minor_rows = [
                    _component_row(index // DIMENSION, index % DIMENSION, forward, reverse, coordinates)[columns.start : columns.stop]
                    for index in rows
                ]
                assert _modular_rank(minor_rows, prime) == 128
            details.append(f"F_{prime}@{residue}")
        return "generic full-rank modular witnesses " + ", ".join(details)

    def check_generic_characteristic_zero_rank() -> str:
        certificate = data["generic_rank_over_Qq"]
        assert certificate["rank_over_Qq"] == 256
        assert certificate["nullity_over_Qq"] == 0
        assert certificate["proof_status"] == "[THEOREM]"
        assert certificate["exceptional_set_status"] == "UNRESOLVED_NOT_FACTORED"
        witnesses = certificate["modular_witnesses"]
        assert len(witnesses) >= 2
        for witness in witnesses:
            assert witness in data["modular_generic_evidence"]
        assert "[THEOREM]" in data["status"]
        return "nonzero fixed-row minors modulo p prove rank 256 over Q(q), with roots still unclassified"

    def check_symbolic_wall_record() -> str:
        wall = data["symbolic_all_q_status"]
        assert wall["status"] == "UNRESOLVED"
        assert wall["matrix_shape"] == [128, 128]
        assert wall["peak_rss_bytes"] > 0
        assert wall["elapsed_seconds"]
        sectors = wall["sectors"]
        assert set(sectors) == {"parity_preserving", "parity_reversing"}
        for row in sectors.values():
            assert row["status"] == "timed_out"
            assert row["timeout_seconds"] == 120
            assert float(row["elapsed_seconds"]) >= 120
            assert row["peak_rss_bytes"] > 0
        return "two 128x128 exact symbolic determinant walls carry elapsed time and RSS"

    _run_check("artifact envelope and corrected geometry", check_envelope_and_geometry)
    _run_check("cabling order and exact parity sector support", check_cable_order_and_sector_support)
    _run_check("fixed-q exact maximal minors", check_fixed_q_maximal_minors)
    _run_check("modular generic evidence", check_modular_generic_evidence)
    _run_check("generic characteristic-zero rank theorem", check_generic_characteristic_zero_rank)
    _run_check("symbolic all-q wall record", check_symbolic_wall_record)
    print("PASS")


if __name__ == "__main__":
    main()
