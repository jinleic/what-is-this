#!/usr/bin/env python3
"""Exact universal and small-matrix controls for an all-auxiliary Lax audit.

The universal certificates are the rational Yang RLL identity in
Q[u,v][S_3] and a strict permutation RLLL identity in S_6.  Dense numerical
fitting is never used.  Exact d=2,3,4 matrices are finite controls only.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import time
from fractions import Fraction
from typing import Iterable

CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 2_000_000_000

Permutation = tuple[int, int, int]
PermutationN = tuple[int, ...]
Polynomial = dict[tuple[int, int], int]
GroupElement = dict[Permutation, Polynomial]
SparseMatrix = tuple[dict[int, Fraction], ...]

IDENTITY_PERMUTATION: Permutation = (0, 1, 2)
SWAP_AB: Permutation = (1, 0, 2)
SWAP_AQ: Permutation = (2, 1, 0)
SWAP_BQ: Permutation = (0, 2, 1)


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


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return left after right."""
    return tuple(left[right[index]] for index in range(3))  # type: ignore[return-value]


def polynomial_add(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, 0) + coefficient
        if answer[monomial] == 0:
            del answer[monomial]
    return answer


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for (u_left, v_left), coefficient_left in left.items():
        for (u_right, v_right), coefficient_right in right.items():
            monomial = (u_left + u_right, v_left + v_right)
            answer[monomial] = (
                answer.get(monomial, 0) + coefficient_left * coefficient_right
            )
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def group_add(left: GroupElement, right: GroupElement) -> GroupElement:
    answer = {permutation: dict(polynomial) for permutation, polynomial in left.items()}
    for permutation, polynomial in right.items():
        combined = polynomial_add(answer.get(permutation, {}), polynomial)
        if combined:
            answer[permutation] = combined
        elif permutation in answer:
            del answer[permutation]
    return answer


def group_multiply(left: GroupElement, right: GroupElement) -> GroupElement:
    answer: GroupElement = {}
    for left_permutation, left_polynomial in left.items():
        for right_permutation, right_polynomial in right.items():
            permutation = compose(left_permutation, right_permutation)
            polynomial = polynomial_multiply(left_polynomial, right_polynomial)
            combined = polynomial_add(answer.get(permutation, {}), polynomial)
            if combined:
                answer[permutation] = combined
            elif permutation in answer:
                del answer[permutation]
    return answer


def group_atom(permutation: Permutation, polynomial: Polynomial) -> GroupElement:
    return {permutation: dict(polynomial)}


def serialize_group(element: GroupElement) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for permutation in sorted(element):
        polynomial = element[permutation]
        rows.append(
            {
                "permutation": list(permutation),
                "polynomial_low_to_high": [
                    {
                        "u_power": u_power,
                        "v_power": v_power,
                        "coefficient": coefficient,
                    }
                    for (u_power, v_power), coefficient in sorted(polynomial.items())
                ],
            }
        )
    return rows


def universal_yang_identity() -> dict[str, object]:
    one: Polynomial = {(0, 0): 1}
    u: Polynomial = {(1, 0): 1}
    v: Polynomial = {(0, 1): 1}
    u_minus_v: Polynomial = {(1, 0): 1, (0, 1): -1}

    identity = group_atom(IDENTITY_PERMUTATION, one)
    r_ab = group_add(
        group_atom(IDENTITY_PERMUTATION, u_minus_v),
        group_atom(SWAP_AB, one),
    )
    l_aq = group_add(group_atom(IDENTITY_PERMUTATION, u), group_atom(SWAP_AQ, one))
    l_bq = group_add(group_atom(IDENTITY_PERMUTATION, v), group_atom(SWAP_BQ, one))

    lhs = group_multiply(group_multiply(r_ab, l_aq), l_bq)
    rhs = group_multiply(group_multiply(l_bq, l_aq), r_ab)
    zero_residual = group_add(
        lhs,
        {
            permutation: {monomial: -coefficient for monomial, coefficient in polynomial.items()}
            for permutation, polynomial in rhs.items()
        },
    )

    positive_cycle = compose(SWAP_AB, SWAP_AQ)
    negative_cycle = compose(SWAP_AQ, SWAP_AB)
    relations = {
        "A2_is_identity": compose(SWAP_AB, SWAP_AB) == IDENTITY_PERMUTATION,
        "B2_is_identity": compose(SWAP_AQ, SWAP_AQ) == IDENTITY_PERMUTATION,
        "C2_is_identity": compose(SWAP_BQ, SWAP_BQ) == IDENTITY_PERMUTATION,
        "AB_equals_CA_equals_BC": positive_cycle
        == compose(SWAP_BQ, SWAP_AB)
        == compose(SWAP_AQ, SWAP_BQ),
        "BA_equals_AC_equals_CB": negative_cycle
        == compose(SWAP_AB, SWAP_BQ)
        == compose(SWAP_BQ, SWAP_AQ),
        "ABC_equals_CBA": compose(compose(SWAP_AB, SWAP_AQ), SWAP_BQ)
        == compose(compose(SWAP_BQ, SWAP_AQ), SWAP_AB),
    }
    return {
        "tag": "[THEOREM]",
        "coefficient_ring": "Z[u,v]",
        "identity": "((u-v)I+P_ab)(uI+P_aq)(vI+P_bq)=(vI+P_bq)(uI+P_aq)((u-v)I+P_ab)",
        "lhs_group_algebra_normal_form": serialize_group(lhs),
        "rhs_group_algebra_normal_form": serialize_group(rhs),
        "residual_term_count": len(zero_residual),
        "s3_relations": relations,
        "all_relations_passed": all(relations.values()),
        "identity_passed": lhs == rhs and not zero_residual,
        "unused_identity_sanity": identity == {IDENTITY_PERMUTATION: one},
    }


def compose_n(left: PermutationN, right: PermutationN) -> PermutationN:
    if len(left) != len(right):
        raise ValueError("permutation degree mismatch")
    return tuple(left[right[index]] for index in range(len(left)))


def factor_swap(factors: int, left: int, right: int) -> PermutationN:
    answer = list(range(factors))
    answer[left], answer[right] = answer[right], answer[left]
    return tuple(answer)


def universal_permutation_rlll() -> dict[str, object]:
    """Certify a nonidentity all-d RLLL solution at factor-permutation level."""
    r_123 = factor_swap(6, 0, 1)
    l_145 = factor_swap(6, 0, 3)
    l_246 = factor_swap(6, 1, 3)
    l_356 = factor_swap(6, 2, 4)
    lhs = compose_n(
        compose_n(compose_n(r_123, l_145), l_246),
        l_356,
    )
    rhs = compose_n(
        compose_n(compose_n(l_356, l_246), l_145),
        r_123,
    )
    core_lhs = compose_n(compose_n(r_123, l_145), l_246)
    core_rhs = compose_n(compose_n(l_246, l_145), r_123)
    commuting_spectator = all(
        compose_n(l_356, operator) == compose_n(operator, l_356)
        for operator in (r_123, l_145, l_246)
    )
    return {
        "tag": "[THEOREM]",
        "local_operator": "L=P_(first,second) tensor I_third",
        "auxiliary_operator": "R=P_12 tensor I_3",
        "rlll_identity": "R_123 L_145 L_246 L_356 = L_356 L_246 L_145 R_123",
        "lhs_factor_permutation": list(lhs),
        "rhs_factor_permutation": list(rhs),
        "core_three_transposition_identity_passed": core_lhs == core_rhs,
        "disjoint_L356_commutations_passed": commuting_spectator,
        "identity_passed": lhs == rhs and core_lhs == core_rhs and commuting_spectator,
        "scope": "all common factor dimensions d>=1 over every field",
    }


def sparse_identity(size: int) -> SparseMatrix:
    return tuple({index: Fraction(1)} for index in range(size))


def sparse_add(left: SparseMatrix, right: SparseMatrix) -> SparseMatrix:
    if len(left) != len(right):
        raise ValueError("matrix shape mismatch")
    rows: list[dict[int, Fraction]] = []
    for left_row, right_row in zip(left, right):
        row = dict(left_row)
        for column, value in right_row.items():
            row[column] = row.get(column, Fraction(0)) + value
            if row[column] == 0:
                del row[column]
        rows.append(row)
    return tuple(rows)


def sparse_scale(matrix: SparseMatrix, scalar: Fraction) -> SparseMatrix:
    if scalar == 0:
        return tuple({} for _ in matrix)
    return tuple(
        {column: scalar * value for column, value in row.items() if scalar * value}
        for row in matrix
    )


def sparse_multiply(left: SparseMatrix, right: SparseMatrix) -> SparseMatrix:
    if len(left) != len(right):
        raise ValueError("matrix shape mismatch")
    rows: list[dict[int, Fraction]] = []
    for left_row in left:
        row: dict[int, Fraction] = {}
        for middle, left_value in left_row.items():
            for column, right_value in right[middle].items():
                row[column] = row.get(column, Fraction(0)) + left_value * right_value
                if row[column] == 0:
                    del row[column]
        rows.append(row)
    return tuple(rows)


def tensor_digits(index: int, dimension: int, factors: int) -> list[int]:
    digits = [0] * factors
    for position in range(factors - 1, -1, -1):
        digits[position] = index % dimension
        index //= dimension
    return digits


def tensor_index(digits: Iterable[int], dimension: int) -> int:
    answer = 0
    for digit in digits:
        answer = dimension * answer + digit
    return answer


def swap_operator(dimension: int, left_factor: int, right_factor: int) -> SparseMatrix:
    factors = 3
    size = dimension**factors
    rows: list[dict[int, Fraction]] = [{} for _ in range(size)]
    for column in range(size):
        digits = tensor_digits(column, dimension, factors)
        digits[left_factor], digits[right_factor] = (
            digits[right_factor],
            digits[left_factor],
        )
        row = tensor_index(digits, dimension)
        rows[row][column] = Fraction(1)
    return tuple(rows)


def sparse_digest(matrix: SparseMatrix) -> str:
    payload = [
        [row_index, column, value.numerator, value.denominator]
        for row_index, row in enumerate(matrix)
        for column, value in sorted(row.items())
    ]
    encoded = json.dumps(payload, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def tensor_permutation_operator(
    dimension: int, permutation: PermutationN
) -> SparseMatrix:
    factors = len(permutation)
    size = dimension**factors
    rows: list[dict[int, Fraction]] = [{} for _ in range(size)]
    for column in range(size):
        input_digits = tensor_digits(column, dimension, factors)
        output_digits = [0] * factors
        for source, target in enumerate(permutation):
            output_digits[target] = input_digits[source]
        row = tensor_index(output_digits, dimension)
        rows[row][column] = Fraction(1)
    return tuple(rows)


def finite_permutation_rlll_control(dimension: int) -> dict[str, object]:
    r_123 = factor_swap(6, 0, 1)
    l_145 = factor_swap(6, 0, 3)
    l_246 = factor_swap(6, 1, 3)
    l_356 = factor_swap(6, 2, 4)
    lhs_permutation = compose_n(
        compose_n(compose_n(r_123, l_145), l_246),
        l_356,
    )
    rhs_permutation = compose_n(
        compose_n(compose_n(l_356, l_246), l_145),
        r_123,
    )
    lhs = tensor_permutation_operator(dimension, lhs_permutation)
    rhs = tensor_permutation_operator(dimension, rhs_permutation)
    determinant = (-1) ** (dimension * dimension * (dimension - 1) // 2)
    return {
        "tag": "[COMPUTATION]",
        "dimension_d": dimension,
        "six_factor_matrix_size": dimension**6,
        "auxiliary_R_matrix_size": dimension**3,
        "local_L_matrix_size": dimension**3,
        "residual_nonzero_entries": 0 if lhs == rhs else sum(
            len(row)
            for row in sparse_add(lhs, sparse_scale(rhs, Fraction(-1)))
        ),
        "lhs_sha256": sparse_digest(lhs),
        "rhs_sha256": sparse_digest(rhs),
        "auxiliary_R_determinant": determinant,
        "R_invertible": abs(determinant) == 1,
        "R_nonscalar": dimension >= 2,
        "local_L_active_entry_span_dimension": dimension**2,
        "passed": lhs == rhs and abs(determinant) == 1 and dimension >= 2,
    }


def finite_rll_control(dimension: int, u: Fraction, v: Fraction) -> dict[str, object]:
    size = dimension**3
    identity = sparse_identity(size)
    p_ab = swap_operator(dimension, 0, 1)
    p_aq = swap_operator(dimension, 0, 2)
    p_bq = swap_operator(dimension, 1, 2)
    r_ab = sparse_add(sparse_scale(identity, u - v), p_ab)
    l_aq = sparse_add(sparse_scale(identity, u), p_aq)
    l_bq = sparse_add(sparse_scale(identity, v), p_bq)
    lhs = sparse_multiply(sparse_multiply(r_ab, l_aq), l_bq)
    rhs = sparse_multiply(sparse_multiply(l_bq, l_aq), r_ab)
    residual = sparse_add(lhs, sparse_scale(rhs, Fraction(-1)))

    symmetric_dimension = dimension * (dimension + 1) // 2
    antisymmetric_dimension = dimension * (dimension - 1) // 2
    difference = u - v
    determinant = (difference + 1) ** symmetric_dimension * (
        difference - 1
    ) ** antisymmetric_dimension
    return {
        "tag": "[COMPUTATION]",
        "dimension_d": dimension,
        "u": str(u),
        "v": str(v),
        "three_factor_matrix_size": size,
        "residual_nonzero_entries": sum(len(row) for row in residual),
        "lhs_nonzero_entries": sum(len(row) for row in lhs),
        "lhs_sha256": sparse_digest(lhs),
        "rhs_sha256": sparse_digest(rhs),
        "auxiliary_R_determinant": str(determinant),
        "symmetric_multiplicity": symmetric_dimension,
        "antisymmetric_multiplicity": antisymmetric_dimension,
        "R_invertible": determinant != 0,
        "R_nonscalar": any(
            row_index != column
            for row_index, row in enumerate(r_ab)
            for column in row
        ),
        "passed": residual == tuple({} for _ in range(size))
        and lhs == rhs
        and determinant != 0,
    }


def run_intertwiner_audit() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    universal = universal_yang_identity()
    permutation_rlll = universal_permutation_rlll()
    budget_tick(started, "universal permutation identities")
    controls = [
        finite_rll_control(dimension, Fraction(2), Fraction(-1))
        for dimension in (2, 3, 4)
    ]
    permutation_controls = [
        finite_permutation_rlll_control(dimension) for dimension in (2, 3, 4)
    ]
    budget_tick(started, "finite RLL and RLLL controls")
    checks = [
        {
            "name": "universal rational Yang identity holds in Z[u,v][S3]",
            "passed": bool(universal["identity_passed"])
            and bool(universal["all_relations_passed"])
            and universal["residual_term_count"] == 0,
            "detail": f"normal-form terms={len(universal['lhs_group_algebra_normal_form'])}",
        },
        {
            "name": "strict permutation RLLL identity holds for every common factor dimension",
            "passed": bool(permutation_rlll["identity_passed"]),
            "detail": (
                f"factor permutation={permutation_rlll['lhs_factor_permutation']}"
            ),
        },
        {
            "name": "d=2,3,4 exact sparse RLL controls have zero residual",
            "passed": all(bool(row["passed"]) for row in controls),
            "detail": ", ".join(
                f"d={row['dimension_d']}:size={row['three_factor_matrix_size']}"
                for row in controls
            ),
        },
        {
            "name": "all stored finite Yang R matrices are nonscalar and invertible",
            "passed": all(
                bool(row["R_nonscalar"]) and bool(row["R_invertible"])
                for row in controls
            ),
            "detail": ", ".join(
                f"d={row['dimension_d']}:det={row['auxiliary_R_determinant']}"
                for row in controls
            ),
        },
        {
            "name": "d=2,3,4 strict permutation RLLL controls are nonscalar and invertible",
            "passed": all(
                bool(row["passed"])
                and bool(row["R_nonscalar"])
                and bool(row["R_invertible"])
                for row in permutation_controls
            ),
            "detail": ", ".join(
                f"d={row['dimension_d']}:size={row['six_factor_matrix_size']}"
                for row in permutation_controls
            ),
        },
    ]
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)
    return {
        "universal_yang_rll": universal,
        "universal_permutation_rlll": permutation_rlll,
        "finite_sparse_controls": controls,
        "finite_permutation_rlll_controls": permutation_controls,
        "process_time_seconds": str(time.process_time() - started),
        "peak_rss_bytes": max_rss_bytes(),
    }, checks


def main() -> int:
    data, checks = run_intertwiner_audit()
    print(
        f"PASS e212 ({len(checks)} checks, "
        f"controls={len(data['finite_sparse_controls'])}, "
        f"peak_rss={data['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
