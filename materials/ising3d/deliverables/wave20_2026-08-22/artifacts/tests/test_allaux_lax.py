#!/usr/bin/env python3
"""Clean-room exact verifier for proofs/allaux_lax.md.

No experiment module is imported.  The verifier independently rebuilds the
universal S3 identity, finite RLL matrices, exact intertwiner ranks, tensor
commutants, local-entry envelopes, and packed-Pauli controls.
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
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "allaux_lax.json"
CPU_BUDGET_SECONDS = 60.0
RSS_CAP_BYTES = 2_000_000_000
STARTED = time.process_time()
PASSED: list[str] = []
FAILED: list[str] = []

Dense = tuple[tuple[Fraction, ...], ...]
Permutation = tuple[int, int, int]
PermutationN = tuple[int, ...]
Polynomial = dict[tuple[int, int], int]
GroupElement = dict[Permutation, Polynomial]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(stage: str) -> None:
    used = time.process_time() - STARTED
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(f"process-time budget exceeded at {stage}: {used}")
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {rss}")


def check(name: str, condition: bool, detail: str = "") -> None:
    print(
        f"[{'PASS' if condition else 'FAIL'}] {name}"
        + (f": {detail}" if detail else ""),
        flush=True,
    )
    (PASSED if condition else FAILED).append(name)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Independent polynomial group-algebra implementation.
def permutation_product(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[position]] for position in range(3))  # type: ignore[return-value]


def polynomial_sum(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for power, coefficient in right.items():
        result[power] = result.get(power, 0) + coefficient
        if result[power] == 0:
            del result[power]
    return result


def polynomial_product(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            power = (a + c, b + d)
            result[power] = result.get(power, 0) + x * y
    return {power: coefficient for power, coefficient in result.items() if coefficient}


def group_sum(left: GroupElement, right: GroupElement) -> GroupElement:
    result = {permutation: dict(polynomial) for permutation, polynomial in left.items()}
    for permutation, polynomial in right.items():
        combined = polynomial_sum(result.get(permutation, {}), polynomial)
        if combined:
            result[permutation] = combined
        else:
            result.pop(permutation, None)
    return result


def group_product(left: GroupElement, right: GroupElement) -> GroupElement:
    result: GroupElement = {}
    for p, f in left.items():
        for q, g in right.items():
            permutation = permutation_product(p, q)
            result[permutation] = polynomial_sum(
                result.get(permutation, {}), polynomial_product(f, g)
            )
            if not result[permutation]:
                del result[permutation]
    return result


def atom(permutation: Permutation, polynomial: Polynomial) -> GroupElement:
    return {permutation: polynomial}


def serialized_group(element: GroupElement) -> list[dict[str, object]]:
    return [
        {
            "permutation": list(permutation),
            "polynomial_low_to_high": [
                {
                    "u_power": a,
                    "v_power": b,
                    "coefficient": coefficient,
                }
                for (a, b), coefficient in sorted(element[permutation].items())
            ],
        }
        for permutation in sorted(element)
    ]


def rebuild_universal_identity() -> tuple[GroupElement, GroupElement]:
    identity = (0, 1, 2)
    ab = (1, 0, 2)
    aq = (2, 1, 0)
    bq = (0, 2, 1)
    one = {(0, 0): 1}
    r = group_sum(atom(identity, {(1, 0): 1, (0, 1): -1}), atom(ab, one))
    la = group_sum(atom(identity, {(1, 0): 1}), atom(aq, one))
    lb = group_sum(atom(identity, {(0, 1): 1}), atom(bq, one))
    return group_product(group_product(r, la), lb), group_product(
        group_product(lb, la), r
    )


def generic_product(left: PermutationN, right: PermutationN) -> PermutationN:
    return tuple(left[right[position]] for position in range(len(left)))


def generic_swap(degree: int, first: int, second: int) -> PermutationN:
    permutation = list(range(degree))
    permutation[first], permutation[second] = (
        permutation[second],
        permutation[first],
    )
    return tuple(permutation)


def rebuild_strict_rlll_permutations() -> tuple[PermutationN, PermutationN]:
    r = generic_swap(6, 0, 1)
    l145 = generic_swap(6, 0, 3)
    l246 = generic_swap(6, 1, 3)
    l356 = generic_swap(6, 2, 4)
    lhs = generic_product(
        generic_product(generic_product(r, l145), l246),
        l356,
    )
    rhs = generic_product(
        generic_product(generic_product(l356, l246), l145),
        r,
    )
    return lhs, rhs


def generic_digits(index: int, dimension: int, factors: int) -> list[int]:
    answer = [0] * factors
    for position in range(factors - 1, -1, -1):
        answer[position] = index % dimension
        index //= dimension
    return answer


def generic_index(digits: Sequence[int], dimension: int) -> int:
    answer = 0
    for digit in digits:
        answer = dimension * answer + digit
    return answer


def permutation_operator_digest(
    dimension: int, permutation: PermutationN
) -> str:
    payload: list[list[int]] = []
    for column in range(dimension ** len(permutation)):
        source = generic_digits(column, dimension, len(permutation))
        target = [0] * len(permutation)
        for source_position, target_position in enumerate(permutation):
            target[target_position] = source[source_position]
        payload.append([generic_index(target, dimension), column, 1, 1])
    payload.sort()
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode()
    ).hexdigest()


def strict_rlll_finite_result(dimension: int) -> dict[str, object]:
    lhs, rhs = rebuild_strict_rlll_permutations()
    determinant = (-1) ** (dimension * dimension * (dimension - 1) // 2)
    return {
        "dimension_d": dimension,
        "six_factor_matrix_size": dimension**6,
        "auxiliary_R_matrix_size": dimension**3,
        "local_L_matrix_size": dimension**3,
        "residual_nonzero_entries": 0 if lhs == rhs else dimension**6,
        "lhs_sha256": permutation_operator_digest(dimension, lhs),
        "rhs_sha256": permutation_operator_digest(dimension, rhs),
        "auxiliary_R_determinant": determinant,
        "R_invertible": abs(determinant) == 1,
        "R_nonscalar": dimension >= 2,
        "local_L_active_entry_span_dimension": dimension**2,
        "passed": lhs == rhs and abs(determinant) == 1,
    }


# Sparse integer matrices, represented differently from the producer.
def base_d_digits(index: int, dimension: int) -> tuple[int, int, int]:
    return (
        index // (dimension * dimension),
        (index // dimension) % dimension,
        index % dimension,
    )


def base_d_index(digits: Sequence[int], dimension: int) -> int:
    return digits[0] * dimension * dimension + digits[1] * dimension + digits[2]


def swap_map(dimension: int, first: int, second: int) -> dict[tuple[int, int], int]:
    result: dict[tuple[int, int], int] = {}
    for column in range(dimension**3):
        output = list(base_d_digits(column, dimension))
        output[first], output[second] = output[second], output[first]
        result[(base_d_index(output, dimension), column)] = 1
    return result


def identity_map(size: int) -> dict[tuple[int, int], int]:
    return {(index, index): 1 for index in range(size)}


def add_maps(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
    left_scale: int = 1,
    right_scale: int = 1,
) -> dict[tuple[int, int], int]:
    result: dict[tuple[int, int], int] = {}
    for key, value in left.items():
        if left_scale * value:
            result[key] = left_scale * value
    for key, value in right.items():
        result[key] = result.get(key, 0) + right_scale * value
        if result[key] == 0:
            del result[key]
    return result


def multiply_maps(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
) -> dict[tuple[int, int], int]:
    right_by_row: dict[int, list[tuple[int, int]]] = {}
    for (row, column), value in right.items():
        right_by_row.setdefault(row, []).append((column, value))
    result: dict[tuple[int, int], int] = {}
    for (row, middle), left_value in left.items():
        for column, right_value in right_by_row.get(middle, ()):
            key = (row, column)
            result[key] = result.get(key, 0) + left_value * right_value
            if result[key] == 0:
                del result[key]
    return result


def sparse_map_digest(matrix: dict[tuple[int, int], int]) -> str:
    payload = [
        [row, column, value, 1]
        for (row, column), value in sorted(matrix.items())
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode()
    ).hexdigest()


def finite_rll(dimension: int) -> dict[str, object]:
    size = dimension**3
    identity = identity_map(size)
    r = add_maps(identity, swap_map(dimension, 0, 1), 3, 1)
    la = add_maps(identity, swap_map(dimension, 0, 2), 2, 1)
    lb = add_maps(identity, swap_map(dimension, 1, 2), -1, 1)
    lhs = multiply_maps(multiply_maps(r, la), lb)
    rhs = multiply_maps(multiply_maps(lb, la), r)
    symmetric = dimension * (dimension + 1) // 2
    antisymmetric = dimension * (dimension - 1) // 2
    return {
        "dimension_d": dimension,
        "three_factor_matrix_size": size,
        "residual_nonzero_entries": len(add_maps(lhs, rhs, 1, -1)),
        "lhs_nonzero_entries": len(lhs),
        "lhs_sha256": sparse_map_digest(lhs),
        "rhs_sha256": sparse_map_digest(rhs),
        "auxiliary_R_determinant": str(4**symmetric * 2**antisymmetric),
        "R_invertible": True,
        "R_nonscalar": any(row != column for row, column in r),
        "passed": lhs == rhs,
    }


# Exact rational linear algebra and representations.
def qmatrix(rows: Iterable[Iterable[int | Fraction]]) -> Dense:
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


def identity(dimension: int) -> Dense:
    return tuple(
        tuple(Fraction(row == column) for column in range(dimension))
        for row in range(dimension)
    )


def unit(dimension: int, selected_row: int, selected_column: int) -> Dense:
    return tuple(
        tuple(
            Fraction(row == selected_row and column == selected_column)
            for column in range(dimension)
        )
        for row in range(dimension)
    )


def matrix_units(dimension: int) -> tuple[Dense, ...]:
    return tuple(
        unit(dimension, row, column)
        for row in range(dimension)
        for column in range(dimension)
    )


def add(left: Dense, right: Dense) -> Dense:
    return tuple(
        tuple(x + y for x, y in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def scale(matrix: Dense, scalar: int | Fraction) -> Dense:
    return tuple(
        tuple(Fraction(scalar) * value for value in row) for row in matrix
    )


def multiply(left: Dense, right: Dense) -> Dense:
    transposed = tuple(zip(*right))
    return tuple(
        tuple(sum(x * y for x, y in zip(row, column)) for column in transposed)
        for row in left
    )


def kron(left: Dense, right: Dense) -> Dense:
    return tuple(
        tuple(
            left[i][j] * right[k][ell]
            for j in range(len(left[0]))
            for ell in range(len(right[0]))
        )
        for i in range(len(left))
        for k in range(len(right))
    )


def vectorize(matrix: Dense) -> tuple[Fraction, ...]:
    return tuple(value for row in matrix for value in row)


def rank_q(rows: Iterable[Sequence[int | Fraction]], column_count: int) -> int:
    matrix = [
        [Fraction(value) for value in row]
        for row in rows
        if any(Fraction(value) for value in row)
    ]
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot_value for value in matrix[pivot_row]]
        for row in range(pivot_row + 1, len(matrix)):
            if matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    value - factor * pivot_entry
                    for value, pivot_entry in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(matrix) or pivot_row == column_count:
            break
    return pivot_row


def hom_nullity(source: Sequence[Dense], target: Sequence[Dense]) -> int:
    source_dimension = len(source[0])
    target_dimension = len(target[0])
    unknown_count = source_dimension * target_dimension
    equations: list[list[Fraction]] = []
    for a, b in zip(source, target):
        for i in range(target_dimension):
            for j in range(source_dimension):
                row = [Fraction(0)] * unknown_count
                for k in range(source_dimension):
                    row[i * source_dimension + k] += a[k][j]
                for ell in range(target_dimension):
                    row[ell * source_dimension + j] -= b[i][ell]
                equations.append(row)
    return unknown_count - rank_q(equations, unknown_count)


def tensor_commutant_nullity(auxiliary: int, physical: int) -> tuple[int, int]:
    representation = tuple(
        kron(identity(auxiliary), generator) for generator in matrix_units(physical)
    )
    nullity = hom_nullity(representation, representation)
    witnesses = tuple(
        kron(generator, identity(physical)) for generator in matrix_units(auxiliary)
    )
    witness_rank = rank_q(
        (vectorize(witness) for witness in witnesses), (auxiliary * physical) ** 2
    )
    return nullity, witness_rank


def dense_digest(matrices: Sequence[Dense]) -> str:
    payload = [
        [[value.numerator, value.denominator] for value in vectorize(matrix)]
        for matrix in matrices
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode()
    ).hexdigest()


def local_entry_result(dimension: int) -> dict[str, object]:
    size = dimension**2
    permutation = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for a in range(dimension):
        for q in range(dimension):
            permutation[q * dimension + a][a * dimension + q] = Fraction(1)
    local_l = add(scale(identity(size), 2), qmatrix(permutation))
    blocks: dict[tuple[int, int], Dense] = {}
    for a_out in range(dimension):
        for a_in in range(dimension):
            blocks[(a_out, a_in)] = tuple(
                tuple(
                    local_l[a_out * dimension + q_out][a_in * dimension + q_in]
                    for q_in in range(dimension)
                )
                for q_out in range(dimension)
            )
    formula = all(
        block
        == add(
            scale(identity(dimension), 2 if a_out == a_in else 0),
            unit(dimension, a_in, a_out),
        )
        for (a_out, a_in), block in blocks.items()
    )
    generated = {
        vectorize(blocks[(a_out, a_in)])
        for a_out in range(dimension)
        for a_in in range(dimension)
        if a_out != a_in
    }
    for row in range(dimension):
        other = 0 if row else 1
        generated.add(
            vectorize(multiply(blocks[(other, row)], blocks[(row, other)]))
        )
    expected = {vectorize(matrix) for matrix in matrix_units(dimension)}
    return {
        "dimension_d": dimension,
        "local_L_matrix_size": size,
        "block_formula_passed": formula,
        "generated_matrix_unit_count": len(generated),
        "associative_span_rank_Q": rank_q(generated, dimension**2),
        "expected_full_matrix_dimension": dimension**2,
        "generated_set_is_all_matrix_units": generated == expected,
        "local_L_sha256": dense_digest((local_l,)),
        "passed": formula and generated == expected,
    }


def tensor_unit_result(dimension: int, sites: int) -> dict[str, object]:
    locals_ = matrix_units(dimension)
    products: list[Dense] = []
    for choices in itertools.product(range(dimension**2), repeat=sites):
        result = locals_[choices[0]]
        for choice in choices[1:]:
            result = kron(result, locals_[choice])
        products.append(result)
    hilbert = dimension**sites
    vectors = [vectorize(product) for product in products]
    return {
        "dimension_d": dimension,
        "sites_n": sites,
        "hilbert_dimension": hilbert,
        "tensor_product_count": len(products),
        "span_rank_Q": rank_q(vectors, hilbert**2),
        "expected_endomorphism_dimension": dimension ** (2 * sites),
        "products_are_all_global_matrix_units": set(vectors)
        == {vectorize(matrix) for matrix in matrix_units(hilbert)},
        "tensor_units_sha256": dense_digest(products),
    }


# Packed-Pauli controls.
def anticommutes(left: int, right: int, qubits: int) -> bool:
    mask = (1 << qubits) - 1
    return (
        ((left & mask) & (right >> qubits)).bit_count()
        + ((left >> qubits) & (right & mask)).bit_count()
    ) % 2 == 1


def pauli_closure(raw_generators: Sequence[int], qubits: int) -> frozenset[int]:
    generators = tuple(dict.fromkeys(raw_generators))
    reached = set(generators)
    queue = deque(generators)
    while queue:
        parent = queue.popleft()
        for generator in generators:
            if anticommutes(parent, generator, qubits):
                child = parent ^ generator
                if child not in reached:
                    reached.add(child)
                    queue.append(child)
    return frozenset(reached)


def label_digest(labels: Iterable[int]) -> str:
    return hashlib.sha256(
        ",".join(str(label) for label in sorted(labels)).encode()
    ).hexdigest()


def full_control_closure(qubits: int) -> frozenset[int]:
    raw: list[int] = []
    for site in range(qubits):
        raw.extend((1 << site, 1 << (qubits + site)))
    raw.extend(
        (1 << (qubits + site)) | (1 << (qubits + site + 1))
        for site in range(qubits - 1)
    )
    return pauli_closure(raw, qubits)


def xxx_split_closure(qubits: int) -> frozenset[int]:
    raw: list[int] = []
    for site in range(qubits - 1):
        xx = (1 << site) | (1 << (site + 1))
        zz = xx << qubits
        raw.extend((xx, xx | zz, zz))
    return pauli_closure(raw, qubits)


def expected_xxx_labels(qubits: int) -> frozenset[int]:
    mask = (1 << qubits) - 1
    central_options = (0, mask) if qubits % 2 == 0 else (0,)
    radical = {
        x | (z << qubits) for x in central_options for z in central_options
    }
    return frozenset(
        x | (z << qubits)
        for x in range(1 << qubits)
        if x.bit_count() % 2 == 0
        for z in range(1 << qubits)
        if z.bit_count() % 2 == 0
        if (x | (z << qubits)) not in radical
    )


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check(
        "artifact envelope",
        set(artifact) == {"meta", "data", "checks"}
        and len(artifact["checks"]) == 16
        and all(row["passed"] is True for row in artifact["checks"]),
    )

    lhs, rhs = rebuild_universal_identity()
    stored_universal = artifact["data"]["exact_intertwiner_controls"][
        "universal_yang_rll"
    ]
    check(
        "universal Z[u,v][S3] Yang identity",
        lhs == rhs
        and serialized_group(lhs) == stored_universal["lhs_group_algebra_normal_form"]
        and serialized_group(rhs) == stored_universal["rhs_group_algebra_normal_form"]
        and stored_universal["residual_term_count"] == 0,
        f"normal-form terms={len(lhs)}",
    )

    strict_lhs, strict_rhs = rebuild_strict_rlll_permutations()
    stored_strict = artifact["data"]["exact_intertwiner_controls"][
        "universal_permutation_rlll"
    ]
    check(
        "universal strict permutation RLLL identity",
        strict_lhs == strict_rhs
        and list(strict_lhs) == stored_strict["lhs_factor_permutation"]
        and list(strict_rhs) == stored_strict["rhs_factor_permutation"]
        and stored_strict["identity_passed"] is True
        and strict_lhs != tuple(range(6)),
        f"factor permutation={list(strict_lhs)}",
    )
    budget_tick("universal identity")

    stored_rll = {
        row["dimension_d"]: row
        for row in artifact["data"]["exact_intertwiner_controls"][
            "finite_sparse_controls"
        ]
    }
    for dimension in (2, 3, 4):
        rebuilt = finite_rll(dimension)
        stored = stored_rll[dimension]
        keys = (
            "three_factor_matrix_size",
            "residual_nonzero_entries",
            "lhs_nonzero_entries",
            "lhs_sha256",
            "rhs_sha256",
            "auxiliary_R_determinant",
            "R_invertible",
            "R_nonscalar",
        )
        check(
            f"d={dimension} exact nonscalar invertible RLL control",
            rebuilt["passed"] is True
            and all(rebuilt[key] == stored[key] for key in keys),
            f"size={rebuilt['three_factor_matrix_size']}, det={rebuilt['auxiliary_R_determinant']}",
        )

    stored_strict_controls = {
        row["dimension_d"]: row
        for row in artifact["data"]["exact_intertwiner_controls"][
            "finite_permutation_rlll_controls"
        ]
    }
    for dimension in (2, 3, 4):
        rebuilt = strict_rlll_finite_result(dimension)
        stored = stored_strict_controls[dimension]
        keys = (
            "six_factor_matrix_size",
            "auxiliary_R_matrix_size",
            "local_L_matrix_size",
            "residual_nonzero_entries",
            "lhs_sha256",
            "rhs_sha256",
            "auxiliary_R_determinant",
            "R_invertible",
            "R_nonscalar",
            "local_L_active_entry_span_dimension",
        )
        check(
            f"d={dimension} exact nonscalar invertible strict RLLL control",
            rebuilt["passed"] is True
            and all(rebuilt[key] == stored[key] for key in keys),
            f"size={rebuilt['six_factor_matrix_size']}, det={rebuilt['auxiliary_R_determinant']}",
        )
    budget_tick("finite RLL matrices")

    algebra = artifact["data"]["exact_algebra_controls"]
    stored_entries = {row["dimension_d"]: row for row in algebra["local_l_entry_envelopes"]}
    for dimension in (2, 3, 4):
        rebuilt = local_entry_result(dimension)
        stored = stored_entries[dimension]
        keys = (
            "local_L_matrix_size",
            "block_formula_passed",
            "generated_matrix_unit_count",
            "associative_span_rank_Q",
            "expected_full_matrix_dimension",
            "generated_set_is_all_matrix_units",
            "local_L_sha256",
        )
        check(
            f"d={dimension} local L entries generate full M_d",
            rebuilt["passed"] is True
            and all(rebuilt[key] == stored[key] for key in keys),
            f"rank={rebuilt['associative_span_rank_Q']}",
        )

    stored_tensors = {
        (row["dimension_d"], row["sites_n"]): row
        for row in algebra["tensor_product_entry_envelopes"]
    }
    for dimension, sites in ((2, 3), (3, 2)):
        rebuilt = tensor_unit_result(dimension, sites)
        stored = stored_tensors[(dimension, sites)]
        keys = (
            "hilbert_dimension",
            "tensor_product_count",
            "span_rank_Q",
            "expected_endomorphism_dimension",
            "products_are_all_global_matrix_units",
            "tensor_units_sha256",
        )
        check(
            f"d={dimension},n={sites} tensor entries span full chain End",
            all(rebuilt[key] == stored[key] for key in keys),
            f"rank={rebuilt['span_rank_Q']}",
        )
    budget_tick("entry envelopes")

    units2 = matrix_units(2)
    units3 = matrix_units(3)
    similarity = qmatrix(((1, 1), (0, 1)))
    inverse = qmatrix(((1, -1), (0, 1)))
    conjugate = tuple(multiply(multiply(similarity, generator), inverse) for generator in units2)
    recomputed_nullities = {
        "full_M2_commutant": hom_nullity(units2, units2),
        "full_M3_commutant": hom_nullity(units3, units3),
        "equivalent_absolutely_irreducible_M2_modules": hom_nullity(units2, conjugate),
        "nonisomorphic_one_dimensional_k_plus_k_modules": hom_nullity(
            (qmatrix(((1,),)), qmatrix(((0,),))),
            (qmatrix(((1,),)), qmatrix(((1,),))),
        ),
    }
    schur_rows = {row["case"]: row for row in algebra["schur_controls"]["rows"]}
    schur_pass = all(
        schur_rows[name]["intertwiner_nullity_Q"] == value
        for name, value in recomputed_nullities.items()
    )
    one_operator = schur_rows[
        "one_operator_equation_is_not_full_algebra_intertwining"
    ]
    schur_pass = (
        schur_pass
        and hom_nullity((identity(2),), (identity(2),))
        == one_operator["single_operator_nullity_Q"]
        == 4
        and recomputed_nullities["full_M2_commutant"]
        == one_operator["full_M2_nullity_Q"]
        == 1
    )
    check(
        "Schur and one-operator intertwiner nullities",
        schur_pass,
        str(recomputed_nullities),
    )

    tensor_rows = {
        (row["auxiliary_dimension"], row["physical_dimension"]): row
        for row in algebra["tensor_factor_commutants"]
    }
    tensor_pass = True
    tensor_details: list[str] = []
    for auxiliary, physical in ((2, 2), (3, 2), (2, 3)):
        nullity, witness_rank = tensor_commutant_nullity(auxiliary, physical)
        stored = tensor_rows[(auxiliary, physical)]
        tensor_pass &= (
            nullity == stored["commutant_nullity_Q"] == auxiliary**2
            and witness_rank == stored["explicit_witness_rank_Q"] == auxiliary**2
        )
        tensor_details.append(f"a={auxiliary},p={physical}:dim={nullity}")
    check(
        "tensor-factor commutants are End(auxiliary) tensor I",
        tensor_pass,
        ", ".join(tensor_details),
    )
    budget_tick("intertwiner ranks")

    stored_full = {row["qubits_n"]: row for row in algebra["full_control_lie_closures"]}
    full_pass = True
    for qubits in range(1, 5):
        closure = full_control_closure(qubits)
        stored = stored_full[qubits]
        full_pass &= (
            closure == frozenset(range(1, 1 << (2 * qubits)))
            and len(closure) == stored["lie_dimension_Q"] == 4**qubits - 1
            and label_digest(closure) == stored["closure_sha256"]
        )
    check("small full-control Lie closures", full_pass, "dimensions 3,15,63,255")

    stored_xxx = {row["qubits_n"]: row for row in algebra["xxx_pauli_summand_closures"]}
    xxx_pass = True
    xxx_details: list[str] = []
    for qubits in range(3, 7):
        closure = xxx_split_closure(qubits)
        expected = expected_xxx_labels(qubits)
        stored = stored_xxx[qubits]
        xxx_pass &= (
            closure == expected
            and len(closure) == stored["lie_dimension_Q"]
            and label_digest(closure) == stored["closure_sha256"]
            and stored["finite_scope_only"] is True
        )
        xxx_details.append(f"n={qubits}:dim={len(closure)}")
    check("finite XXX Pauli-summand closures", xxx_pass, ", ".join(xxx_details))
    budget_tick("Pauli closures")

    source_hashes = artifact["meta"]["source_sha256"]
    sources_pass = all(
        sha256_file(ROOT / relative) == digest
        for relative, digest in source_hashes.items()
    )
    scope = artifact["data"]["existing_theorem_scope_guard"]
    audited_paths = {
        "tetra_ungraded": ROOT / "proofs" / "tetra_ungraded.md",
        "tetra16_locus": ROOT / "proofs" / "tetra16_locus.md",
        "local_term_trichotomy": ROOT / "proofs" / "clifford_grade_classification.md",
    }
    scope_hashes_pass = all(
        sha256_file(path) == scope["source_sha256"][name]
        for name, path in audited_paths.items()
    )
    ungraded = audited_paths["tetra_ungraded"].read_text().lower()
    cabled = audited_paths["tetra16_locus"].read_text().lower()
    trichotomy = audited_paths["local_term_trichotomy"].read_text().lower()
    phrase_pass = (
        "does not classify higher auxiliary dimensions" in ungraded
        and "fixed-order one-line-cabled" in cabled
        and "does not establish an rlll solution for an unrestricted" in cabled
        and "g_gamma = lie_r" in trichotomy
        and "act on the full `2^n`-dimensional spin hilbert space" in trichotomy
    )
    check(
        "source hashes and prior theorem scope guards",
        sources_pass and scope_hashes_pass and phrase_pass,
    )

    verdict = artifact["data"]["verdict"]
    formal = artifact["data"]["formal_rlll_intertwiner"]
    check(
        "scope-safe method-limitation verdict",
        verdict["tag"] == "[THEOREM]"
        and verdict["ising_conclusion"].startswith("[UNRESOLVED]")
        and "intersection" in formal["one_operator_intertwiner_identity"]
        and formal["missing_implication"].startswith("[UNRESOLVED]")
        and "neither reproves nor enlarges"
        in scope["preserved_8x8_statement"],
    )

    budget_tick("final")
    if FAILED:
        print(f"FAIL ({len(FAILED)} failures): {', '.join(FAILED)}")
        return 1
    print(
        f"PASS test_allaux_lax ({len(PASSED)} checks, "
        f"peak_rss={max_rss_bytes()})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
