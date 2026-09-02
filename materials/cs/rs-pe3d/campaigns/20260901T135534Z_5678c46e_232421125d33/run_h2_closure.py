#!/usr/bin/env python3
"""Exact controls for gate H-MIX-H2-CLOSURE.

Standalone by design: imports nothing from a frozen campaign directory.
Mathematical decisions use Python integer arithmetic modulo p only.
"""
import sys

sys.dont_write_bytecode = True

import itertools
import json
import os
from pathlib import Path
import time
import traceback

RUN_ID = "20260901T135534Z_5678c46e_232421125d33"
GATE = "H-MIX-H2-CLOSURE"
CPU_CAP_SECONDS = 600.0
WALL_CAP_SECONDS = 900.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT_DIR = Path(__file__).resolve().parent
RESULTS = {
    "schema": "h-mix-h2-closure-controls/1",
    "run_id": RUN_ID,
    "gate": GATE,
    "arithmetic": "exact Python integers modulo p; no numerical algebra package",
    "prereg_commit": "b3e61f76fdb72c06a65cb90a61069273ef4644f3",
    "source_frozen_run": "20260901T132343Z_4f169cbb_2f81dd8f2d8d",
    "controls": {},
}
ASSERTION_COUNT = 0


class BudgetExceeded(RuntimeError):
    pass


def cpu_seconds():
    return time.process_time() - START_CPU


def wall_seconds():
    return time.monotonic() - START_WALL


def check_budget(stage):
    cpu = cpu_seconds()
    wall = wall_seconds()
    if cpu >= CPU_CAP_SECONDS:
        raise BudgetExceeded(
            f"hard CPU cap {CPU_CAP_SECONDS}s reached at {stage}: {cpu:.6f}s"
        )
    if wall >= WALL_CAP_SECONDS:
        raise BudgetExceeded(
            f"hard wall cap {WALL_CAP_SECONDS}s reached at {stage}: {wall:.6f}s"
        )


def check(condition, tag, detail=None):
    global ASSERTION_COUNT
    ASSERTION_COUNT += 1
    if not condition:
        suffix = "" if detail is None else f": {detail}"
        raise AssertionError(f"{tag}{suffix}")
    RESULTS["controls"][tag] = {"pass": True, "detail": detail}


def pin(got, want, tag):
    check(got == want, tag, {"got": got, "want": want})


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return [safe(v) for v in sorted(value)]
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    return value


def rref_with_pivots(rows, p):
    matrix = [[x % p for x in row] for row in rows]
    if not matrix:
        return matrix, []
    nrows = len(matrix)
    ncols = len(matrix[0])
    pivot_row = 0
    pivots = []
    for col in range(ncols):
        chosen = None
        for row in range(pivot_row, nrows):
            if matrix[row][col] % p:
                chosen = row
                break
        if chosen is None:
            continue
        matrix[pivot_row], matrix[chosen] = matrix[chosen], matrix[pivot_row]
        inv = pow(matrix[pivot_row][col], -1, p)
        matrix[pivot_row] = [(x * inv) % p for x in matrix[pivot_row]]
        for row in range(nrows):
            if row == pivot_row or matrix[row][col] % p == 0:
                continue
            factor = matrix[row][col]
            matrix[row] = [
                (matrix[row][j] - factor * matrix[pivot_row][j]) % p
                for j in range(ncols)
            ]
        pivots.append(col)
        pivot_row += 1
        if pivot_row == nrows:
            break
    return matrix, pivots


def rank_cols(columns, p):
    if not columns:
        return 0
    dim = len(columns[0])
    rows = [[column[i] % p for column in columns] for i in range(dim)]
    _, pivots = rref_with_pivots(rows, p)
    return len(pivots)


def dependent(columns, p):
    return rank_cols(columns, p) < len(columns)


def spark_cols(columns, p):
    for size in range(1, len(columns) + 1):
        for indices in itertools.combinations(range(len(columns)), size):
            if dependent([columns[i] for i in indices], p):
                return size
    return None


class Factor:
    def __init__(self, name, rows, p):
        self.name = name
        self.p = p
        self.rows = [[x % p for x in row] for row in rows]
        check(bool(self.rows), f"{name}.nonempty_rows")
        self.r = len(self.rows)
        self.s = len(self.rows[0])
        check(
            all(len(row) == self.s for row in self.rows),
            f"{name}.rectangular",
        )
        self.cols = [
            [self.rows[i][j] for i in range(self.r)] for j in range(self.s)
        ]
        check(
            all(any(x % p for x in column) for column in self.cols),
            f"{name}.no_zero_columns",
        )
        self.spark = spark_cols(self.cols, p)

    def circuit_supports(self, size):
        supports = []
        for support in itertools.combinations(range(self.s), size):
            columns = [self.cols[i] for i in support]
            if not dependent(columns, self.p):
                continue
            if all(
                not dependent([self.cols[i] for i in proper], self.p)
                for proper in itertools.combinations(support, size - 1)
            ):
                supports.append(tuple(support))
        return supports


def vandermonde_2(name, points, p):
    return Factor(name, [[1 for _ in points], list(points)], p)


def kernel_vector_on_circuit(factor, support):
    support = tuple(sorted(support))
    columns = [factor.cols[i] for i in support]
    matrix = [
        [columns[j][row] for j in range(len(support))]
        for row in range(factor.r)
    ]
    reduced, pivots = rref_with_pivots(matrix, factor.p)
    free = [j for j in range(len(support)) if j not in pivots]
    if len(pivots) != len(support) - 1 or len(free) != 1:
        raise AssertionError(
            f"{factor.name}: support {support} does not have one-dimensional relation space"
        )
    free_col = free[0]
    local = [0] * len(support)
    local[free_col] = 1
    for row, pivot_col in enumerate(pivots):
        local[pivot_col] = (-reduced[row][free_col]) % factor.p
    if not all(local):
        raise AssertionError(
            f"{factor.name}: circuit relation on {support} lacks full support: {local}"
        )
    full = [0] * factor.s
    for index, coefficient in zip(support, local):
        full[index] = coefficient
    for row in range(factor.r):
        total = sum(
            factor.rows[row][index] * full[index] for index in support
        ) % factor.p
        if total:
            raise AssertionError(
                f"{factor.name}: solved relation fails at row {row}: {total}"
            )
    return full


def product_columns(factor_a, factor_b):
    if factor_a.p != factor_b.p:
        raise AssertionError("factor fields differ")
    p = factor_a.p
    return {
        (u, v): [
            factor_a.cols[u][i] * factor_b.cols[v][j] % p
            for i in range(factor_a.r)
            for j in range(factor_b.r)
        ]
        for u in range(factor_a.s)
        for v in range(factor_b.s)
    }


def profile(support):
    return (
        len({u for u, _ in support}),
        len({v for _, v in support}),
    )


def zero_matrix(rows, cols):
    return [[0 for _ in range(cols)] for _ in range(rows)]


def add_matrices(left, right, p):
    return [
        [(left[u][v] + right[u][v]) % p for v in range(len(left[0]))]
        for u in range(len(left))
    ]


def matrix_relation(factor_a, factor_b, gamma):
    p = factor_a.p
    return [
        [
            sum(
                factor_a.rows[i][u]
                * gamma[u][v]
                * factor_b.rows[j][v]
                for u in range(factor_a.s)
                for v in range(factor_b.s)
            )
            % p
            for j in range(factor_b.r)
        ]
        for i in range(factor_a.r)
    ]


def coefficient_sum(factor_a, factor_b, gamma):
    p = factor_a.p
    columns = product_columns(factor_a, factor_b)
    total = [0] * (factor_a.r * factor_b.r)
    for (u, v), column in columns.items():
        coefficient = gamma[u][v] % p
        for index, value in enumerate(column):
            total[index] = (total[index] + coefficient * value) % p
    return total


def support_of(gamma, p):
    return frozenset(
        (u, v)
        for u, row in enumerate(gamma)
        for v, value in enumerate(row)
        if value % p
    )


def is_circuit_support(factor_a, factor_b, support):
    keys = tuple(sorted(support))
    columns = product_columns(factor_a, factor_b)
    if rank_cols([columns[key] for key in keys], factor_a.p) != len(keys) - 1:
        return False
    return all(
        rank_cols([columns[keys[i]] for i in positions], factor_a.p)
        == len(keys) - 1
        for positions in itertools.combinations(range(len(keys)), len(keys) - 1)
    )


def crossing_family(factor_a, factor_b):
    family = set()
    circuits_a = factor_a.circuit_supports(3)
    circuits_b = factor_b.circuit_supports(3)
    for circuit_a in circuits_a:
        for center_a in circuit_a:
            for circuit_b in circuits_b:
                for center_b in circuit_b:
                    support = frozenset(
                        [(u, center_b) for u in circuit_a if u != center_a]
                        + [(center_a, v) for v in circuit_b if v != center_b]
                    )
                    family.add(support)
    return family


def owner_accept_and_rejects():
    p = 13
    factor_a = vandermonde_2("T1.owner.A", (1, 2, 3, 4), p)
    factor_b = vandermonde_2("T1.owner.B", (1, 2, 3, 4), p)
    pin(factor_a.spark, 3, "T1.owner.spark_A")
    pin(factor_b.spark, 3, "T1.owner.spark_B")

    c = [1, 11, 1, 0]
    delta = [1, 11, 1, 0]
    for factor, relation, tag in (
        (factor_a, c, "A"),
        (factor_b, delta, "B"),
    ):
        row_sums = [
            sum(factor.rows[row][i] * relation[i] for i in range(factor.s)) % p
            for row in range(factor.r)
        ]
        pin(row_sums, [0, 0], f"T1.owner.factor_relation_{tag}")

    u1, u2, w = 0, 1, 2
    v1, v2, x = 0, 1, 2
    alpha = (delta[v1], delta[v2])
    beta = ((-c[u1]) % p, (-c[u2]) % p)
    c1 = zero_matrix(factor_a.s, factor_b.s)
    c2 = zero_matrix(factor_a.s, factor_b.s)
    for column, scalar in zip((v1, v2), alpha):
        for row, coefficient in enumerate(c):
            c1[row][column] = scalar * coefficient % p
    for row, scalar in zip((u1, u2), beta):
        for column, coefficient in enumerate(delta):
            c2[row][column] = scalar * coefficient % p
    gamma = add_matrices(c1, c2, p)
    observed = support_of(gamma, p)
    crossing = frozenset({(w, v1), (w, v2), (u1, x), (u2, x)})
    expected = frozenset({(0, 2), (1, 2), (2, 0), (2, 1)})
    pin(observed, expected, "T1.owner.support_exact")
    pin(observed, crossing, "T1.owner.support_equals_crossing")
    pin(matrix_relation(factor_a, factor_b, gamma), [[0, 0], [0, 0]],
        "T1.owner.AGammaBT_zero")
    pin(coefficient_sum(factor_a, factor_b, gamma), [0, 0, 0, 0],
        "T1.owner.coefficient_sum_zero")
    check(
        is_circuit_support(factor_a, factor_b, observed),
        "T1.owner.genuine_circuit",
    )
    columns = product_columns(factor_a, factor_b)
    pin(
        rank_cols([columns[key] for key in sorted(observed)], p),
        3,
        "T1.owner.rank_three",
    )

    perturbed = [list(row) for row in gamma]
    perturbed[w][v1] = (perturbed[w][v1] + 1) % p
    perturbed_matrix = matrix_relation(factor_a, factor_b, perturbed)
    perturbed_sum = coefficient_sum(factor_a, factor_b, perturbed)
    check(
        any(value for row in perturbed_matrix for value in row),
        "T4.coefficient_perturbation_reject_matrix",
        {"perturbed_cell": [w, v1], "AGammaBT": perturbed_matrix},
    )
    check(
        any(perturbed_sum),
        "T4.coefficient_perturbation_reject_sum",
        {"coefficient_sum": perturbed_sum},
    )

    duplicate_rows = [row + [row[0]] for row in factor_b.rows]
    duplicate_factor = Factor("T4.duplicate.B", duplicate_rows, p)
    pin(duplicate_factor.spark, 2, "T4.duplicate.spark_drop")
    dependent_pairs = [
        pair
        for pair in itertools.combinations(range(duplicate_factor.s), 2)
        if dependent([duplicate_factor.cols[i] for i in pair], p)
    ]
    check(bool(dependent_pairs), "T4.duplicate.dependent_pair_detected")

    RESULTS["T1_owner_GF13"] = {
        "factor_points": [1, 2, 3, 4],
        "T": [0, 1, 2],
        "Z": [0, 1, 2],
        "c": c[:3],
        "delta": delta[:3],
        "alpha": list(alpha),
        "beta": list(beta),
        "support": observed,
        "crossing_support": crossing,
        "AGammaBT": matrix_relation(factor_a, factor_b, gamma),
        "rank": 3,
        "all_three_subsets_independent": True,
        "circuit": True,
    }
    RESULTS["T4_rejects"] = {
        "coefficient_perturbation": {
            "cell": [w, v1],
            "delta": 1,
            "AGammaBT": perturbed_matrix,
            "coefficient_sum": perturbed_sum,
            "reject_fired": True,
        },
        "duplicated_column": {
            "original_spark": factor_b.spark,
            "planted_spark": duplicate_factor.spark,
            "dependent_pairs": dependent_pairs,
            "reject_fired": True,
        },
    }
    return factor_a, factor_b


def sweep_22(factor_a, factor_b, label, expected_crossings, shape_constant):
    p = factor_a.p
    pin(factor_a.spark, 3, f"T2.{label}.spark_A")
    pin(factor_b.spark, 3, f"T2.{label}.spark_B")
    circuit_vectors_a = {
        support: kernel_vector_on_circuit(factor_a, support)
        for support in factor_a.circuit_supports(3)
    }
    circuit_vectors_b = {
        support: kernel_vector_on_circuit(factor_b, support)
        for support in factor_b.circuit_supports(3)
    }
    expected_all_a = len(list(itertools.combinations(range(factor_a.s), 3)))
    expected_all_b = len(list(itertools.combinations(range(factor_b.s), 3)))
    pin(len(circuit_vectors_a), expected_all_a, f"T2.{label}.all_A_triples_circuits")
    pin(len(circuit_vectors_b), expected_all_b, f"T2.{label}.all_B_triples_circuits")

    constructed = crossing_family(factor_a, factor_b)
    pin(len(constructed), expected_crossings, f"T2.{label}.crossing_family_size")
    expected_shapes = (
        len(list(itertools.combinations(range(factor_a.s), 2)))
        * len(list(itertools.combinations(range(factor_b.s), 2)))
        * (factor_a.s - 2) ** 2
        * (factor_b.s - 2) ** 2
    )
    pin(expected_shapes, shape_constant, f"T2.{label}.shape_formula")
    expected_normalized_tuples = expected_shapes * (p - 1) ** 3

    shape_count = 0
    normalized_tuples = 0
    cancellation_solutions = 0
    size_four_solutions = 0
    profile_counts = {}
    profile33_solution_count = 0
    profile33_supports = set()
    new_supports = set()
    product = product_columns(factor_a, factor_b)

    for active_rows in itertools.combinations(range(factor_a.s), 2):
        u1, u2 = active_rows
        remaining_rows = tuple(i for i in range(factor_a.s) if i not in active_rows)
        for active_columns in itertools.combinations(range(factor_b.s), 2):
            v1, v2 = active_columns
            remaining_columns = tuple(
                j for j in range(factor_b.s) if j not in active_columns
            )
            for w1, w2 in itertools.product(remaining_rows, repeat=2):
                support_c1 = tuple(sorted((u1, u2, w1)))
                support_c2 = tuple(sorted((u1, u2, w2)))
                c1_direction = circuit_vectors_a[support_c1]
                c2_direction = circuit_vectors_a[support_c2]
                for x1, x2 in itertools.product(remaining_columns, repeat=2):
                    check_budget(f"T2.{label}.shape{shape_count}")
                    support_d1 = tuple(sorted((v1, v2, x1)))
                    support_d2 = tuple(sorted((v1, v2, x2)))
                    d1_direction = circuit_vectors_b[support_d1]
                    d2_direction = circuit_vectors_b[support_d2]
                    shape_count += 1
                    for alpha2 in range(1, p):
                        for beta1 in range(1, p):
                            for beta2 in range(1, p):
                                normalized_tuples += 1
                                if (
                                    c1_direction[u1]
                                    + beta1 * d1_direction[v1]
                                ) % p:
                                    continue
                                if (
                                    c1_direction[u2]
                                    + beta2 * d2_direction[v1]
                                ) % p:
                                    continue
                                if (
                                    alpha2 * c2_direction[u1]
                                    + beta1 * d1_direction[v2]
                                ) % p:
                                    continue
                                if (
                                    alpha2 * c2_direction[u2]
                                    + beta2 * d2_direction[v2]
                                ) % p:
                                    continue
                                cancellation_solutions += 1
                                c1 = zero_matrix(factor_a.s, factor_b.s)
                                c2 = zero_matrix(factor_a.s, factor_b.s)
                                for row in support_c1:
                                    c1[row][v1] = c1_direction[row]
                                for row in support_c2:
                                    c1[row][v2] = (
                                        alpha2 * c2_direction[row]
                                    ) % p
                                for column in support_d1:
                                    c2[u1][column] = (
                                        beta1 * d1_direction[column]
                                    ) % p
                                for column in support_d2:
                                    c2[u2][column] = (
                                        beta2 * d2_direction[column]
                                    ) % p
                                gamma = add_matrices(c1, c2, p)
                                support = support_of(gamma, p)
                                if any(gamma[u][v] % p for u in active_rows for v in active_columns):
                                    raise AssertionError(
                                        f"T2.{label}: retained tuple did not cancel all of X"
                                    )
                                if len(support) != 4:
                                    raise AssertionError(
                                        f"T2.{label}: retained tuple has support size {len(support)}: {sorted(support)}"
                                    )
                                if any(
                                    value
                                    for row in matrix_relation(factor_a, factor_b, gamma)
                                    for value in row
                                ):
                                    raise AssertionError(
                                        f"T2.{label}: retained decomposition is not in the product kernel"
                                    )
                                size_four_solutions += 1
                                prof = profile(support)
                                profile_counts[prof] = profile_counts.get(prof, 0) + 1
                                if prof != (3, 3):
                                    continue
                                profile33_solution_count += 1
                                profile33_supports.add(support)
                                if support not in constructed:
                                    new_supports.add(support)
                                if rank_cols([product[key] for key in sorted(support)], p) != 3:
                                    raise AssertionError(
                                        f"T2.{label}: profile-(3,3) candidate rank is not three"
                                    )
                                if not all(
                                    rank_cols([product[key] for key in proper], p) == 3
                                    for proper in itertools.combinations(sorted(support), 3)
                                ):
                                    raise AssertionError(
                                        f"T2.{label}: profile-(3,3) candidate is not a circuit"
                                    )

    pin(shape_count, expected_shapes, f"T2.{label}.all_shapes_swept")
    pin(
        normalized_tuples,
        expected_normalized_tuples,
        f"T2.{label}.all_normalized_scalar_tuples_swept",
    )
    pin(cancellation_solutions, size_four_solutions,
        f"T2.{label}.every_X_cancellation_has_size_four")
    pin(len(new_supports), 0, f"T2.{label}.no_new_profile33_support")
    pin(profile33_solution_count, expected_crossings,
        f"T2.{label}.profile33_solution_count")
    pin(len(profile33_supports), expected_crossings,
        f"T2.{label}.profile33_distinct_support_count")
    pin(profile33_supports, constructed, f"T2.{label}.set_equality_crossing")

    return {
        "prime": p,
        "factor_A": {"name": factor_a.name, "r": factor_a.r, "s": factor_a.s},
        "factor_B": {"name": factor_b.name, "r": factor_b.r, "s": factor_b.s},
        "shape_choices_swept": shape_count,
        "normalized_scalar_tuples_swept": normalized_tuples,
        "raw_scalar_tuples_represented": normalized_tuples * (p - 1),
        "cancellation_solutions": cancellation_solutions,
        "size_four_solutions": size_four_solutions,
        "retained_profile_counts": profile_counts,
        "profile33_solution_count": profile33_solution_count,
        "profile33_distinct_supports": len(profile33_supports),
        "crossing_family_size": len(constructed),
        "new_profile33_supports": new_supports,
        "set_level_containment": True,
        "set_equality": profile33_supports == constructed,
        "all_profile33_candidates_circuits": True,
    }


def exhaustive_parameter_sweeps():
    rows = []
    for p in (7, 11, 13):
        check_budget(f"T2.prime{p}.start")
        a4 = vandermonde_2(f"T2.p{p}.V4.A", (1, 2, 3, 4), p)
        b4 = vandermonde_2(f"T2.p{p}.V4.B", (1, 2, 3, 4), p)
        rows.append(sweep_22(a4, b4, f"p{p}.V4xV4", 144, 576))

        a3 = vandermonde_2(f"T2.p{p}.V3.A", (0, 1, 3), p)
        b5 = vandermonde_2(f"T2.p{p}.V5.B", (0, 1, 2, 4, 5), p)
        rows.append(sweep_22(a3, b5, f"p{p}.V3xV5", 90, 270))
    RESULTS["T2_exhaustive_22_sweeps"] = rows


def previous_campaign_census():
    p = 13
    factor_a = Factor(
        "T3.previous.Awit",
        [[1, 1, 1, 1], [1, 2, 3, 4]],
        p,
    )
    factor_b = Factor(
        "T3.previous.Bwit",
        [[1, -2, 1, 0], [0, 1, -2, 1]],
        p,
    )
    pin(factor_a.spark, 3, "T3.previous.spark_A")
    pin(factor_b.spark, 3, "T3.previous.spark_B")
    columns = product_columns(factor_a, factor_b)
    keys = tuple(sorted(columns))
    circuits = []
    by_profile = {}
    subsets = 0
    for support in itertools.combinations(keys, 4):
        subsets += 1
        if subsets % 128 == 0:
            check_budget(f"T3.previous.subset{subsets}")
        if rank_cols([columns[key] for key in support], p) >= 4:
            continue
        if not all(
            rank_cols([columns[key] for key in proper], p) == 3
            for proper in itertools.combinations(support, 3)
        ):
            continue
        frozen_support = frozenset(support)
        circuits.append(frozen_support)
        prof = profile(frozen_support)
        by_profile.setdefault(prof, set()).add(frozen_support)
    pin(subsets, 1820, "T3.previous.all_four_subsets_swept")
    pin(len(circuits), 156, "T3.previous.total_circuits")
    measured33 = by_profile.get((3, 3), set())
    measured44 = by_profile.get((4, 4), set())
    pin(len(measured33), 144, "T3.previous.profile33_count")
    pin(len(measured44), 12, "T3.previous.profile44_count")
    constructed = crossing_family(factor_a, factor_b)
    pin(len(constructed), 144, "T3.previous.constructed_crossing_count")
    pin(measured33, constructed, "T3.previous.profile33_set_equality")
    pin(len(measured33 - constructed), 0, "T3.previous.zero_measured_residual")
    pin(len(constructed - measured33), 0, "T3.previous.zero_constructed_residual")
    RESULTS["T3_previous_campaign_census"] = {
        "prime": p,
        "subsets_swept": subsets,
        "total_circuits": len(circuits),
        "profile33_measured": len(measured33),
        "profile33_constructed": len(constructed),
        "profile44": len(measured44),
        "measured_minus_crossing": measured33 - constructed,
        "crossing_minus_measured": constructed - measured33,
        "set_equality": True,
        "zero_residual": True,
    }


def tensor_layout_guard():
    p = 13
    factor_a = vandermonde_2("T5.guard.A2", (1, 2, 4, 7), p)
    points_b = (1, 2, 3, 5)
    factor_b = Factor(
        "T5.guard.B3",
        [
            [1 for _ in points_b],
            list(points_b),
            [x * x % p for x in points_b],
        ],
        p,
    )
    actual = product_columns(factor_a, factor_b)
    expected = {
        (u, v): [
            factor_a.cols[u][i] * factor_b.cols[v][j] % p
            for i in range(factor_a.r)
            for j in range(factor_b.r)
        ]
        for u in range(factor_a.s)
        for v in range(factor_b.s)
    }
    pin(actual, expected, "T5.row_major_Kronecker_coordinates")
    actual_dimensions = {len(column) for column in actual.values()}
    pin(actual_dimensions, {6}, "T5.Kronecker_ambient_dimension")
    block_concat = {
        (u, v): factor_a.cols[u] + factor_b.cols[v]
        for u in range(factor_a.s)
        for v in range(factor_b.s)
    }
    concat_dimensions = {len(column) for column in block_concat.values()}
    pin(concat_dimensions, {5}, "T5.planted_concat_dimension")
    check(actual != block_concat, "T5.concat_reject_coordinates_differ")
    RESULTS["T5_tensor_layout_guard"] = {
        "factor_row_dimensions": [factor_a.r, factor_b.r],
        "Kronecker_dimension": 6,
        "block_concat_dimension": 5,
        "sample_key": [1, 2],
        "Kronecker_sample": actual[(1, 2)],
        "block_concat_sample": block_concat[(1, 2)],
        "guard_fired": True,
    }


def runtime_controls():
    pin(os.environ.get("PYTHONDONTWRITEBYTECODE"), "1", "runtime.no_bytecode_env")
    check(sys.dont_write_bytecode, "runtime.no_bytecode_interpreter")
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        pin(os.environ.get(variable), "1", f"runtime.{variable}")
    if hasattr(os, "getpriority") and hasattr(os, "PRIO_PROCESS"):
        priority = os.getpriority(os.PRIO_PROCESS, 0)
        check(priority >= 15, "runtime.nice_at_least_15", {"nice": priority})
    else:
        priority = "platform-unavailable"
    RESULTS["runtime"] = {
        "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
        "sys_dont_write_bytecode": sys.dont_write_bytecode,
        "thread_env": {
            variable: os.environ.get(variable)
            for variable in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
        "nice": priority,
        "CPU_cap_seconds": CPU_CAP_SECONDS,
        "wall_cap_seconds": WALL_CAP_SECONDS,
    }


def main():
    runtime_controls()
    print(f"== {GATE} run {RUN_ID}")
    print(
        f"== exact GF(p), CPU cap {CPU_CAP_SECONDS:.0f}s, wall cap {WALL_CAP_SECONDS:.0f}s"
    )
    owner_accept_and_rejects()
    print("[T1/T4] owner GF(13) ACCEPT and both planted REJECTs passed")
    exhaustive_parameter_sweeps()
    print("[T2] six exhaustive (2,2) parameter sweeps passed set-level containment")
    previous_campaign_census()
    print("[T3] previous 1820-subset census reproduced; profile-(3,3) residual zero")
    tensor_layout_guard()
    print("[T5] Kronecker-vs-block-concatenation fail-loud guard passed")
    check_budget("finalize")
    RESULTS["assertion_count"] = ASSERTION_COUNT
    RESULTS["failures"] = 0
    RESULTS["budget"] = {
        "CPU_seconds": round(cpu_seconds(), 6),
        "wall_seconds": round(wall_seconds(), 6),
        "CPU_cap_seconds": CPU_CAP_SECONDS,
        "wall_cap_seconds": WALL_CAP_SECONDS,
        "within_cap": True,
    }
    output = OUT_DIR / "controls_results.json"
    output.write_text(json.dumps(safe(RESULTS), indent=2, sort_keys=True) + "\n")
    print(
        f"== ALL CONTROLS PASS: {ASSERTION_COUNT} aggregate asserts; "
        f"CPU {cpu_seconds():.3f}s, wall {wall_seconds():.3f}s; {output.name}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        RESULTS["failure"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
            "CPU_seconds": cpu_seconds(),
            "wall_seconds": wall_seconds(),
        }
        broken = OUT_DIR / "broken_results.json"
        broken.write_text(json.dumps(safe(RESULTS), indent=2, sort_keys=True) + "\n")
        print(f"FAIL: {type(exc).__name__}: {exc}; preserved {broken}", file=sys.stderr)
        raise
