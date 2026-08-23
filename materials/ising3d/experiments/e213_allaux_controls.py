#!/usr/bin/env python3
"""Independent exact algebra controls for the all-auxiliary Lax audit.

This producer uses only Fraction arithmetic and packed Pauli labels.  It checks
Schur/intertwiner dimensions, tensor-factor commutants, local L-entry envelopes,
and two small local-term Lie families.  It imports no other experiment.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import time
from collections import deque
from fractions import Fraction
from typing import Iterable, Sequence

CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 2_000_000_000
DenseMatrix = tuple[tuple[Fraction, ...], ...]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    used = time.process_time() - started
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(f"process-time budget exceeded at {stage}: {used}")
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {rss}")


def dense(rows: Iterable[Iterable[int | Fraction]]) -> DenseMatrix:
    answer = tuple(tuple(Fraction(value) for value in row) for row in rows)
    if answer and any(len(row) != len(answer[0]) for row in answer):
        raise ValueError("ragged matrix")
    return answer


def zero_matrix(rows: int, columns: int) -> DenseMatrix:
    return tuple(tuple(Fraction(0) for _ in range(columns)) for _ in range(rows))


def identity_matrix(dimension: int) -> DenseMatrix:
    return tuple(
        tuple(Fraction(row == column) for column in range(dimension))
        for row in range(dimension)
    )


def matrix_unit(dimension: int, row: int, column: int) -> DenseMatrix:
    return tuple(
        tuple(
            Fraction(row_index == row and column_index == column)
            for column_index in range(dimension)
        )
        for row_index in range(dimension)
    )


def matrix_add(left: DenseMatrix, right: DenseMatrix) -> DenseMatrix:
    return tuple(
        tuple(a + b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def matrix_scale(matrix: DenseMatrix, scalar: int | Fraction) -> DenseMatrix:
    factor = Fraction(scalar)
    return tuple(tuple(factor * value for value in row) for row in matrix)


def matrix_multiply(left: DenseMatrix, right: DenseMatrix) -> DenseMatrix:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix shape mismatch")
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def kronecker(left: DenseMatrix, right: DenseMatrix) -> DenseMatrix:
    return tuple(
        tuple(
            left[left_row][left_column] * right[right_row][right_column]
            for left_column in range(len(left[0]))
            for right_column in range(len(right[0]))
        )
        for left_row in range(len(left))
        for right_row in range(len(right))
    )


def flatten(matrix: DenseMatrix) -> tuple[Fraction, ...]:
    return tuple(value for row in matrix for value in row)


def matrix_digest(matrices: Sequence[DenseMatrix]) -> str:
    payload = [
        [[value.numerator, value.denominator] for value in flatten(matrix)]
        for matrix in matrices
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode()
    ).hexdigest()


def rational_rank(rows: Iterable[Sequence[int | Fraction]], columns: int) -> int:
    work = [
        [Fraction(value) for value in row]
        for row in rows
        if any(Fraction(value) != 0 for value in row)
    ]
    rank = 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row_index in range(len(work)):
            if row_index == rank or not work[row_index][column]:
                continue
            factor = work[row_index][column]
            work[row_index] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(work[row_index], work[rank])
            ]
        rank += 1
        if rank == len(work) or rank == columns:
            break
    return rank


def intertwiner_nullity(
    source_representation: Sequence[DenseMatrix],
    target_representation: Sequence[DenseMatrix],
) -> tuple[int, int, int]:
    if len(source_representation) != len(target_representation):
        raise ValueError("representation generator mismatch")
    source_dimension = len(source_representation[0])
    target_dimension = len(target_representation[0])
    unknowns = source_dimension * target_dimension
    equations: list[list[Fraction]] = []
    for source, target in zip(source_representation, target_representation):
        for target_row in range(target_dimension):
            for source_column in range(source_dimension):
                equation = [Fraction(0)] * unknowns
                for source_row in range(source_dimension):
                    equation[target_row * source_dimension + source_row] += source[
                        source_row
                    ][source_column]
                for target_column in range(target_dimension):
                    equation[target_column * source_dimension + source_column] -= target[
                        target_row
                    ][target_column]
                equations.append(equation)
    rank = rational_rank(equations, unknowns)
    return rank, unknowns - rank, len(equations)


def commutator_zero(left: DenseMatrix, right: DenseMatrix) -> bool:
    return matrix_multiply(left, right) == matrix_multiply(right, left)


def full_matrix_units(dimension: int) -> tuple[DenseMatrix, ...]:
    return tuple(
        matrix_unit(dimension, row, column)
        for row in range(dimension)
        for column in range(dimension)
    )


def schur_controls() -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for dimension in (2, 3):
        representation = full_matrix_units(dimension)
        rank, nullity, equation_count = intertwiner_nullity(
            representation, representation
        )
        rows.append(
            {
                "tag": "[COMPUTATION]",
                "case": f"full_M{dimension}_commutant",
                "module_dimension": dimension,
                "equation_count": equation_count,
                "coefficient_rank_Q": rank,
                "intertwiner_nullity_Q": nullity,
                "expected_nullity": 1,
                "passed": nullity == 1,
            }
        )

    rho = full_matrix_units(2)
    similarity = dense(((1, 1), (0, 1)))
    similarity_inverse = dense(((1, -1), (0, 1)))
    sigma = tuple(
        matrix_multiply(matrix_multiply(similarity, generator), similarity_inverse)
        for generator in rho
    )
    rank, nullity, equation_count = intertwiner_nullity(rho, sigma)
    witness_passed = all(
        matrix_multiply(similarity, source)
        == matrix_multiply(target, similarity)
        for source, target in zip(rho, sigma)
    )
    rows.append(
        {
            "tag": "[COMPUTATION]",
            "case": "equivalent_absolutely_irreducible_M2_modules",
            "equation_count": equation_count,
            "coefficient_rank_Q": rank,
            "intertwiner_nullity_Q": nullity,
            "expected_nullity": 1,
            "invertible_similarity_witness": [[1, 1], [0, 1]],
            "witness_determinant": 1,
            "witness_passed": witness_passed,
            "passed": nullity == 1 and witness_passed,
        }
    )

    source = (dense(((1,),)), dense(((0,),)))
    target = (dense(((1,),)), dense(((1,),)))
    rank, nullity, equation_count = intertwiner_nullity(source, target)
    rows.append(
        {
            "tag": "[COMPUTATION]",
            "case": "nonisomorphic_one_dimensional_k_plus_k_modules",
            "equation_count": equation_count,
            "coefficient_rank_Q": rank,
            "intertwiner_nullity_Q": nullity,
            "expected_nullity": 0,
            "passed": nullity == 0,
        }
    )

    identity = identity_matrix(2)
    single_rank, single_nullity, single_equations = intertwiner_nullity(
        (identity,), (identity,)
    )
    full_rank, full_nullity, full_equations = intertwiner_nullity(rho, rho)
    rows.append(
        {
            "tag": "[COMPUTATION]",
            "case": "one_operator_equation_is_not_full_algebra_intertwining",
            "single_operator_equations": single_equations,
            "single_operator_rank_Q": single_rank,
            "single_operator_nullity_Q": single_nullity,
            "full_M2_equations": full_equations,
            "full_M2_rank_Q": full_rank,
            "full_M2_nullity_Q": full_nullity,
            "passed": single_nullity == 4 and full_nullity == 1,
        }
    )
    return {"rows": rows, "all_passed": all(bool(row["passed"]) for row in rows)}


def tensor_commutant_controls() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for auxiliary_dimension, physical_dimension in ((2, 2), (3, 2), (2, 3)):
        physical_units = full_matrix_units(physical_dimension)
        identity_auxiliary = identity_matrix(auxiliary_dimension)
        representation = tuple(
            kronecker(identity_auxiliary, unit) for unit in physical_units
        )
        rank, nullity, equation_count = intertwiner_nullity(
            representation, representation
        )
        expected_witnesses = tuple(
            kronecker(unit, identity_matrix(physical_dimension))
            for unit in full_matrix_units(auxiliary_dimension)
        )
        witness_rank = rational_rank(
            (flatten(witness) for witness in expected_witnesses),
            (auxiliary_dimension * physical_dimension) ** 2,
        )
        witnesses_commute = all(
            commutator_zero(witness, generator)
            for witness in expected_witnesses
            for generator in representation
        )
        rows.append(
            {
                "tag": "[COMPUTATION]",
                "auxiliary_dimension": auxiliary_dimension,
                "physical_dimension": physical_dimension,
                "total_module_dimension": auxiliary_dimension * physical_dimension,
                "equation_count": equation_count,
                "coefficient_rank_Q": rank,
                "commutant_nullity_Q": nullity,
                "expected_auxiliary_square": auxiliary_dimension**2,
                "explicit_witness_rank_Q": witness_rank,
                "all_auxiliary_matrix_unit_witnesses_commute": witnesses_commute,
                "passed": nullity == auxiliary_dimension**2
                and witness_rank == auxiliary_dimension**2
                and witnesses_commute,
            }
        )
    return rows


def swap_two_factor_matrix(dimension: int) -> DenseMatrix:
    size = dimension**2
    rows = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for auxiliary_input in range(dimension):
        for physical_input in range(dimension):
            column = auxiliary_input * dimension + physical_input
            row = physical_input * dimension + auxiliary_input
            rows[row][column] = Fraction(1)
    return dense(rows)


def extract_auxiliary_block(
    matrix: DenseMatrix, dimension: int, auxiliary_row: int, auxiliary_column: int
) -> DenseMatrix:
    return tuple(
        tuple(
            matrix[auxiliary_row * dimension + physical_row][
                auxiliary_column * dimension + physical_column
            ]
            for physical_column in range(dimension)
        )
        for physical_row in range(dimension)
    )


def local_entry_control(dimension: int, spectral_value: Fraction) -> dict[str, object]:
    size = dimension**2
    local_l = matrix_add(
        matrix_scale(identity_matrix(size), spectral_value),
        swap_two_factor_matrix(dimension),
    )
    blocks = {
        (row, column): extract_auxiliary_block(local_l, dimension, row, column)
        for row in range(dimension)
        for column in range(dimension)
    }
    formulas_pass = all(
        block
        == matrix_add(
            matrix_scale(
                identity_matrix(dimension),
                spectral_value if row == column else Fraction(0),
            ),
            matrix_unit(dimension, column, row),
        )
        for (row, column), block in blocks.items()
    )
    generated: set[tuple[Fraction, ...]] = set()
    off_diagonal = [
        blocks[(row, column)]
        for row in range(dimension)
        for column in range(dimension)
        if row != column
    ]
    generated.update(flatten(matrix) for matrix in off_diagonal)
    for row in range(dimension):
        other = 0 if row != 0 else 1
        left = blocks[(other, row)]
        right = blocks[(row, other)]
        generated.add(flatten(matrix_multiply(left, right)))
    canonical = {flatten(unit) for unit in full_matrix_units(dimension)}
    generated_matrices = tuple(
        tuple(
            tuple(vector[row * dimension + column] for column in range(dimension))
            for row in range(dimension)
        )
        for vector in sorted(generated)
    )
    span_rank = rational_rank(
        (flatten(matrix) for matrix in generated_matrices), dimension**2
    )
    return {
        "tag": "[COMPUTATION]",
        "dimension_d": dimension,
        "spectral_value": str(spectral_value),
        "local_L_matrix_size": size,
        "block_formula_passed": formulas_pass,
        "generated_matrix_unit_count": len(generated),
        "associative_span_rank_Q": span_rank,
        "expected_full_matrix_dimension": dimension**2,
        "generated_set_is_all_matrix_units": generated == canonical,
        "local_L_sha256": matrix_digest((local_l,)),
        "passed": formulas_pass
        and generated == canonical
        and span_rank == dimension**2,
    }


def tensor_unit_control(dimension: int, sites: int) -> dict[str, object]:
    local_units = full_matrix_units(dimension)
    products: list[DenseMatrix] = []
    for choices in itertools.product(range(dimension**2), repeat=sites):
        product = local_units[choices[0]]
        for choice in choices[1:]:
            product = kronecker(product, local_units[choice])
        products.append(product)
    hilbert_dimension = dimension**sites
    rank = rational_rank(
        (flatten(matrix) for matrix in products), hilbert_dimension**2
    )
    canonical = set(flatten(unit) for unit in full_matrix_units(hilbert_dimension))
    actual = set(flatten(matrix) for matrix in products)
    return {
        "tag": "[COMPUTATION]",
        "dimension_d": dimension,
        "sites_n": sites,
        "hilbert_dimension": hilbert_dimension,
        "tensor_product_count": len(products),
        "span_rank_Q": rank,
        "expected_endomorphism_dimension": dimension ** (2 * sites),
        "products_are_all_global_matrix_units": actual == canonical,
        "tensor_units_sha256": matrix_digest(products),
        "passed": rank == dimension ** (2 * sites) and actual == canonical,
    }


def symplectic_product(left: int, right: int, qubits: int) -> int:
    mask = (1 << qubits) - 1
    x_left, z_left = left & mask, left >> qubits
    x_right, z_right = right & mask, right >> qubits
    return (
        (x_left & z_right).bit_count() + (z_left & x_right).bit_count()
    ) & 1


def pauli_lie_closure(generators: Sequence[int], qubits: int) -> frozenset[int]:
    raw = tuple(dict.fromkeys(generators))
    seen = set(raw)
    queue = deque(raw)
    while queue:
        parent = queue.popleft()
        for generator in raw:
            if symplectic_product(parent, generator, qubits):
                child = parent ^ generator
                if child not in seen:
                    seen.add(child)
                    queue.append(child)
    return frozenset(seen)


def pauli_digest(labels: Iterable[int]) -> str:
    payload = ",".join(str(label) for label in sorted(labels)).encode()
    return hashlib.sha256(payload).hexdigest()


def full_control_generators(qubits: int) -> tuple[int, ...]:
    generators: list[int] = []
    for site in range(qubits):
        generators.extend((1 << site, 1 << (qubits + site)))
    for site in range(qubits - 1):
        generators.append(
            (1 << (qubits + site)) | (1 << (qubits + site + 1))
        )
    return tuple(generators)


def full_control_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for qubits in range(1, 5):
        generators = full_control_generators(qubits)
        closure = pauli_lie_closure(generators, qubits)
        expected = frozenset(range(1, 1 << (2 * qubits)))
        rows.append(
            {
                "tag": "[COMPUTATION]",
                "qubits_n": qubits,
                "raw_local_generator_count": len(generators),
                "lie_dimension_Q": len(closure),
                "expected_sl_dimension": 4**qubits - 1,
                "closure_is_every_nonidentity_pauli_label": closure == expected,
                "closure_sha256": pauli_digest(closure),
                "passed": closure == expected,
            }
        )
    return rows


def xxx_split_generators(qubits: int) -> tuple[int, ...]:
    generators: list[int] = []
    for site in range(qubits - 1):
        xx = (1 << site) | (1 << (site + 1))
        zz = xx << qubits
        yy = xx | zz
        generators.extend((xx, yy, zz))
    return tuple(generators)


def xxx_expected_noncentral_even_even(qubits: int) -> frozenset[int]:
    mask = (1 << qubits) - 1
    radical = {
        x | (z << qubits)
        for x in ((0, mask) if qubits % 2 == 0 else (0,))
        for z in ((0, mask) if qubits % 2 == 0 else (0,))
    }
    return frozenset(
        x | (z << qubits)
        for x in range(1 << qubits)
        if x.bit_count() % 2 == 0
        for z in range(1 << qubits)
        if z.bit_count() % 2 == 0
        if (x | (z << qubits)) not in radical
    )


def xxx_split_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for qubits in range(3, 7):
        generators = xxx_split_generators(qubits)
        closure = pauli_lie_closure(generators, qubits)
        expected = xxx_expected_noncentral_even_even(qubits)
        expected_dimension = 4 ** (qubits - 1) - (4 if qubits % 2 == 0 else 1)
        rows.append(
            {
                "tag": "[COMPUTATION]",
                "qubits_n": qubits,
                "raw_Pauli_summand_count": len(generators),
                "lie_dimension_Q": len(closure),
                "expected_noncentral_even_even_dimension": expected_dimension,
                "closure_equals_noncentral_even_even_labels": closure == expected,
                "closure_sha256": pauli_digest(closure),
                "finite_scope_only": True,
                "passed": closure == expected and len(closure) == expected_dimension,
            }
        )
    return rows


def run_control_audit() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    schur = schur_controls()
    tensor_commutants = tensor_commutant_controls()
    budget_tick(started, "Schur and tensor commutants")

    local_entries = [
        local_entry_control(dimension, Fraction(2)) for dimension in (2, 3, 4)
    ]
    tensor_units = [tensor_unit_control(2, 3), tensor_unit_control(3, 2)]
    budget_tick(started, "local L-entry envelopes")

    control_rows = full_control_rows()
    split_rows = xxx_split_rows()
    budget_tick(started, "packed Pauli controls")

    checks = [
        {
            "name": "Schur controls distinguish equivalent, inequivalent, and one-operator cases",
            "passed": bool(schur["all_passed"]),
            "detail": ", ".join(
                f"{row['case']}:{row.get('intertwiner_nullity_Q', row.get('single_operator_nullity_Q'))}"
                for row in schur["rows"]  # type: ignore[index]
            ),
        },
        {
            "name": "tensor-factor commutants have exact auxiliary-square dimension",
            "passed": all(bool(row["passed"]) for row in tensor_commutants),
            "detail": ", ".join(
                f"a={row['auxiliary_dimension']},p={row['physical_dimension']}:dim={row['commutant_nullity_Q']}"
                for row in tensor_commutants
            ),
        },
        {
            "name": "rational Yang local entries generate full one-site matrix algebras",
            "passed": all(bool(row["passed"]) for row in local_entries),
            "detail": ", ".join(
                f"d={row['dimension_d']}:rank={row['associative_span_rank_Q']}"
                for row in local_entries
            ),
        },
        {
            "name": "local matrix units tensor to full exact chain endomorphism bases",
            "passed": all(bool(row["passed"]) for row in tensor_units),
            "detail": ", ".join(
                f"d={row['dimension_d']},n={row['sites_n']}:rank={row['span_rank_Q']}"
                for row in tensor_units
            ),
        },
        {
            "name": "entry-derived onsite controls plus nearest-neighbour coupling close to full sl in small qubit controls",
            "passed": all(bool(row["passed"]) for row in control_rows),
            "detail": ", ".join(
                f"n={row['qubits_n']}:dim={row['lie_dimension_Q']}"
                for row in control_rows
            ),
        },
        {
            "name": "XXX Pauli-summand decomposition has exact exponential finite closures",
            "passed": all(bool(row["passed"]) for row in split_rows),
            "detail": ", ".join(
                f"n={row['qubits_n']}:dim={row['lie_dimension_Q']}"
                for row in split_rows
            ),
        },
    ]
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)
    return {
        "schur_controls": schur,
        "tensor_factor_commutants": tensor_commutants,
        "local_l_entry_envelopes": local_entries,
        "tensor_product_entry_envelopes": tensor_units,
        "full_control_lie_closures": control_rows,
        "xxx_pauli_summand_closures": split_rows,
        "process_time_seconds": str(time.process_time() - started),
        "peak_rss_bytes": max_rss_bytes(),
    }, checks


def main() -> int:
    data, checks = run_control_audit()
    print(
        f"PASS e213 ({len(checks)} checks, "
        f"peak_rss={data['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
