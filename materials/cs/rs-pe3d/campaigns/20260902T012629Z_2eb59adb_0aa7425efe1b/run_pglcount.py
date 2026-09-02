#!/usr/bin/env python3
"""Exact controls for H-ALLDISTINCT-PGLCOUNT.

Standalone: imports no frozen campaign code. Mathematical decisions use only
Python integers modulo p.
"""
import sys

sys.dont_write_bytecode = True

import itertools
import json
import math
import os
from pathlib import Path
import time
import traceback

RUN_ID = "20260902T012629Z_2eb59adb_0aa7425efe1b"
GATE = "H-ALLDISTINCT-PGLCOUNT"
PREREG_COMMIT = "bddaa259bc2e7f175453c42e978d91159fc8ac73"
CPU_CAP_SECONDS = 600.0
WALL_CAP_SECONDS = 900.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT_DIR = Path(__file__).resolve().parent
ASSERTION_COUNT = 0
RESULTS = {
    "schema": "h-alldistinct-pglcount-controls/1",
    "run_id": RUN_ID,
    "gate": GATE,
    "prereg_commit": PREREG_COMMIT,
    "arithmetic": "exact Python integers modulo p; canonical PGL matrices",
    "controls": {},
}


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
        raise BudgetExceeded(f"CPU cap reached at {stage}: {cpu:.6f}s")
    if wall >= WALL_CAP_SECONDS:
        raise BudgetExceeded(f"wall cap reached at {stage}: {wall:.6f}s")


def check(condition, tag, detail=None):
    global ASSERTION_COUNT
    ASSERTION_COUNT += 1
    if not condition:
        raise AssertionError(f"{tag}: {detail}")
    RESULTS["controls"][tag] = {"pass": True, "detail": detail}


def pin(got, want, tag):
    check(got == want, tag, {"got": got, "want": want})


def safe(value):
    if isinstance(value, dict):
        return {str(key): safe(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [safe(item) for item in sorted(value)]
    if isinstance(value, (list, tuple)):
        return [safe(item) for item in value]
    return value


def choose(n, k):
    return math.comb(n, k) if n >= k else 0


def rref_with_pivots(rows, p):
    matrix = [[value % p for value in row] for row in rows]
    if not matrix:
        return matrix, []
    nrows = len(matrix)
    ncols = len(matrix[0])
    pivot_row = 0
    pivots = []
    for col in range(ncols):
        chosen = next(
            (row for row in range(pivot_row, nrows) if matrix[row][col] % p),
            None,
        )
        if chosen is None:
            continue
        matrix[pivot_row], matrix[chosen] = matrix[chosen], matrix[pivot_row]
        inv = pow(matrix[pivot_row][col], -1, p)
        matrix[pivot_row] = [(value * inv) % p for value in matrix[pivot_row]]
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
    rows = [
        [column[row] % p for column in columns]
        for row in range(len(columns[0]))
    ]
    _, pivots = rref_with_pivots(rows, p)
    return len(pivots)


def determinant(matrix, p):
    work = [[value % p for value in row] for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise AssertionError("determinant matrix is not square")
    det = 1
    for col in range(size):
        chosen = next((row for row in range(col, size) if work[row][col]), None)
        if chosen is None:
            return 0
        if chosen != col:
            work[col], work[chosen] = work[chosen], work[col]
            det = (-det) % p
        pivot = work[col][col]
        det = det * pivot % p
        inv = pow(pivot, -1, p)
        for row in range(col + 1, size):
            if work[row][col] == 0:
                continue
            factor = work[row][col] * inv % p
            for j in range(col, size):
                work[row][j] = (work[row][j] - factor * work[col][j]) % p
    return det


def dependent(columns, p):
    return rank_cols(columns, p) < len(columns)


def spark_cols(columns, p):
    for size in range(1, len(columns) + 1):
        for indices in itertools.combinations(range(len(columns)), size):
            if dependent([columns[index] for index in indices], p):
                return size
    return None


class Factor:
    def __init__(self, name, rows, p):
        self.name = name
        self.p = p
        self.rows = [[value % p for value in row] for row in rows]
        check(bool(self.rows), f"{name}.nonempty")
        self.r = len(self.rows)
        self.s = len(self.rows[0])
        check(all(len(row) == self.s for row in self.rows), f"{name}.rectangular")
        self.cols = [
            [self.rows[row][col] for row in range(self.r)]
            for col in range(self.s)
        ]
        check(
            all(any(value for value in column) for column in self.cols),
            f"{name}.no_zero_columns",
        )
        self.spark = spark_cols(self.cols, p)


def grs2(name, points, p):
    return Factor(name, [[1] * len(points), list(points)], p)


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
    return len({u for u, _ in support}), len({v for _, v in support})


def support_matrix(columns, support):
    return [
        [columns[cell][row] for cell in support]
        for row in range(len(columns[support[0]]))
    ]


def is_circuit(columns, support, p):
    if rank_cols([columns[cell] for cell in support], p) != len(support) - 1:
        return False
    return all(
        rank_cols([columns[cell] for cell in proper], p) == len(support) - 1
        for proper in itertools.combinations(support, len(support) - 1)
    )


def cross_ratio(values, p):
    a, b, c, d = (value % p for value in values)
    numerator = (a - c) * (b - d) % p
    denominator = (a - d) * (b - c) % p
    if denominator == 0:
        raise AssertionError(f"cross-ratio denominator zero: {values} mod {p}")
    return numerator * pow(denominator, -1, p) % p


def all_distinct_supports(n_a, n_b):
    for rows in itertools.combinations(range(n_a), 4):
        for columns in itertools.combinations(range(n_b), 4):
            for permutation in itertools.permutations(columns):
                yield tuple(zip(rows, permutation))


def cross_ratio_support_set(x_points, y_points, p, stage):
    accepted = set()
    scanned = 0
    for support in all_distinct_supports(len(x_points), len(y_points)):
        scanned += 1
        if scanned % 65536 == 0:
            check_budget(f"{stage}.cross_ratio.{scanned}")
        xs = [x_points[u] for u, _ in support]
        ys = [y_points[v] for _, v in support]
        if cross_ratio(xs, p) == cross_ratio(ys, p):
            accepted.add(support)
    return accepted, scanned


def all_distinct_product_set(x_points, y_points, p, stage):
    factor_a = grs2(f"{stage}.A", x_points, p)
    factor_b = grs2(f"{stage}.B", y_points, p)
    pin(factor_a.spark, 3, f"{stage}.spark_A")
    pin(factor_b.spark, 3, f"{stage}.spark_B")
    columns = product_columns(factor_a, factor_b)
    expected_columns = {
        (u, v): [1, y_points[v] % p, x_points[u] % p,
                 x_points[u] * y_points[v] % p]
        for u in range(len(x_points))
        for v in range(len(y_points))
    }
    pin(columns, expected_columns, f"{stage}.pure_tensor_coordinates")
    accepted = set()
    scanned = 0
    for support in all_distinct_supports(len(x_points), len(y_points)):
        scanned += 1
        if scanned % 65536 == 0:
            check_budget(f"{stage}.product.{scanned}")
        if determinant(support_matrix(columns, support), p) != 0:
            continue
        if not is_circuit(columns, support, p):
            raise AssertionError(f"{stage}: det-zero support is not a circuit")
        accepted.add(support)
    return accepted, scanned, factor_a, factor_b, columns


def full_product_census(x_points, y_points, p, stage):
    factor_a = grs2(f"{stage}.A", x_points, p)
    factor_b = grs2(f"{stage}.B", y_points, p)
    pin(factor_a.spark, 3, f"{stage}.spark_A")
    pin(factor_b.spark, 3, f"{stage}.spark_B")
    columns = product_columns(factor_a, factor_b)
    keys = tuple(sorted(columns))
    accepted44 = set()
    profile_counts = {}
    scanned = 0
    for support in itertools.combinations(keys, 4):
        scanned += 1
        if scanned % 32768 == 0:
            check_budget(f"{stage}.full_product.{scanned}")
        if determinant(support_matrix(columns, support), p) != 0:
            continue
        if not is_circuit(columns, support, p):
            continue
        prof = profile(support)
        profile_counts[prof] = profile_counts.get(prof, 0) + 1
        if prof == (4, 4):
            accepted44.add(tuple(sorted(support)))
    return {
        "accepted44": accepted44,
        "scanned": scanned,
        "profile_counts": profile_counts,
        "factor_a": factor_a,
        "factor_b": factor_b,
        "columns": columns,
    }


def canonical_matrix(matrix, p):
    values = tuple(value % p for value in matrix)
    first = next((value for value in values if value), None)
    if first is None:
        raise AssertionError("zero matrix has no PGL class")
    inverse = pow(first, -1, p)
    return tuple(value * inverse % p for value in values)


def pgl_matrices(p):
    classes = set()
    for a, b, c, d in itertools.product(range(p), repeat=4):
        if (a * d - b * c) % p == 0:
            continue
        classes.add(canonical_matrix((a, b, c, d), p))
    return tuple(sorted(classes))


def mobius_apply(matrix, point, p):
    """Projective points use 0..p-1 for finite values and p for infinity."""
    a, b, c, d = matrix
    infinity = p
    if point == infinity:
        if c == 0:
            return infinity
        return a * pow(c, -1, p) % p
    numerator = (a * point + b) % p
    denominator = (c * point + d) % p
    if denominator == 0:
        return infinity
    return numerator * pow(denominator, -1, p) % p


def pole_zero(matrix, p):
    points = range(p + 1)
    poles = [point for point in points if mobius_apply(matrix, point, p) == p]
    zeros = [point for point in points if mobius_apply(matrix, point, p) == 0]
    if len(poles) != 1 or len(zeros) != 1:
        raise AssertionError(f"PGL map lacks unique pole/zero: {matrix}")
    return poles[0], zeros[0]


def pgl_structure(p):
    matrices = pgl_matrices(p)
    pin(len(matrices), p * (p * p - 1), f"T1.p{p}.PGL_size")
    pair_counts = {}
    classes = {0: 0, 1: 0, 2: 0}
    multiplicative = set(range(1, p))
    for index, matrix in enumerate(matrices):
        if index % 512 == 0:
            check_budget(f"T1.p{p}.PGL.{index}")
        pole, zero = pole_zero(matrix, p)
        if pole == zero:
            raise AssertionError(f"T1.p{p}: pole equals zero")
        pair_counts[(pole, zero)] = pair_counts.get((pole, zero), 0) + 1
        k = int(pole in multiplicative) + int(zero in multiplicative)
        classes[k] += 1
    expected_pairs = p * (p + 1)
    pin(len(pair_counts), expected_pairs, f"T1.p{p}.ordered_pole_zero_pairs")
    pin(set(pair_counts.values()), {p - 1}, f"T1.p{p}.maps_per_pole_zero_pair")
    q = p - 1
    expected_classes = {0: 2 * q, 1: 4 * q * q, 2: q * q * (q - 1)}
    pin(classes, expected_classes, f"T1.p{p}.multiplicative_class_sizes")
    pin(sum(classes.values()), len(matrices), f"T1.p{p}.class_size_identity")
    return matrices, {
        "prime": p,
        "PGL_size": len(matrices),
        "ordered_pole_zero_pairs": len(pair_counts),
        "maps_per_pair": p - 1,
        "class_sizes": classes,
        "class_sum": sum(classes.values()),
    }


def pgl_graph_support_set(matrices, x_points, y_points, p, stage):
    y_index = {value % p: index for index, value in enumerate(y_points)}
    if len(y_index) != len(y_points):
        raise AssertionError(f"{stage}: Y has duplicate field points")
    supports = set()
    emitted_pairs = 0
    valid_histogram = {}
    for index, matrix in enumerate(matrices):
        if index % 256 == 0:
            check_budget(f"{stage}.PGL_graph.{index}")
        valid = []
        for u, x in enumerate(x_points):
            image = mobius_apply(matrix, x % p, p)
            if image == p or image not in y_index:
                continue
            valid.append((u, y_index[image]))
        if len({v for _, v in valid}) != len(valid):
            raise AssertionError(f"{stage}: Möbius map is not injective")
        valid_histogram[len(valid)] = valid_histogram.get(len(valid), 0) + 1
        for chosen in itertools.combinations(valid, 4):
            support = tuple(chosen)
            emitted_pairs += 1
            if support in supports:
                raise AssertionError(
                    f"{stage}: two (M,U) pairs emitted the same support {support}"
                )
            supports.add(support)
    pin(len(supports), emitted_pairs, f"{stage}.PGL_support_uniqueness")
    return supports, valid_histogram, emitted_pairs


def assert_set_equal(left, right, tag):
    check(
        left == right,
        tag,
        {
            "left": len(left),
            "right": len(right),
            "left_only": list(sorted(left - right))[:3],
            "right_only": list(sorted(right - left))[:3],
        },
    )


def full_field_formula(p):
    return p * (p - 1) * choose(p, 4) + p * p * (p - 1) * choose(p - 1, 4)


def multiplicative_formula(p):
    q = p - 1
    return (
        2 * q * choose(q, 4)
        + 4 * q * q * choose(q - 1, 4)
        + q * q * (q - 1) * choose(q - 2, 4)
    )


def full_field_controls(p, matrices, expected):
    stage = f"T2.p{p}.full_field"
    points = tuple(range(p))
    census = full_product_census(points, points, p, stage)
    product_set = census["accepted44"]
    pin(census["scanned"], choose(p * p, 4), f"{stage}.all_product_subsets_scanned")
    cr_set, cr_scanned = cross_ratio_support_set(points, points, p, stage)
    pin(cr_scanned, choose(p, 4) ** 2 * 24, f"{stage}.all_pairings_CR_scanned")
    pgl_set, histogram, emitted = pgl_graph_support_set(
        matrices, points, points, p, stage
    )
    assert_set_equal(product_set, cr_set, f"{stage}.product_equals_CR")
    assert_set_equal(product_set, pgl_set, f"{stage}.product_equals_PGL")
    formula = full_field_formula(p)
    pin(formula, expected, f"{stage}.formula_anchor")
    pin(len(product_set), formula, f"{stage}.three_route_count")
    pin(histogram.get(p, 0), p * (p - 1), f"{stage}.affine_map_count")
    pin(histogram.get(p - 1, 0), p * p * (p - 1),
        f"{stage}.finite_pole_map_count")
    return {
        "prime": p,
        "product_census_subsets": census["scanned"],
        "product_profile_counts": census["profile_counts"],
        "all_distinct_pairings": cr_scanned,
        "product_count": len(product_set),
        "cross_ratio_count": len(cr_set),
        "PGL_count": len(pgl_set),
        "formula": formula,
        "PGL_valid_domain_histogram": histogram,
        "set_equality": True,
    }


def multiplicative_controls(p, matrices, expected):
    stage = f"T3.p{p}.multiplicative"
    points = tuple(range(1, p))
    product_set, product_scanned, factor_a, factor_b, columns = (
        all_distinct_product_set(points, points, p, stage)
    )
    universe = choose(p - 1, 4) ** 2 * 24
    pin(product_scanned, universe, f"{stage}.all_product_pairings_scanned")
    cr_set, cr_scanned = cross_ratio_support_set(points, points, p, stage)
    pin(cr_scanned, universe, f"{stage}.all_CR_pairings_scanned")
    pgl_set, histogram, emitted = pgl_graph_support_set(
        matrices, points, points, p, stage
    )
    assert_set_equal(product_set, cr_set, f"{stage}.product_equals_CR")
    assert_set_equal(product_set, pgl_set, f"{stage}.product_equals_PGL")
    formula = multiplicative_formula(p)
    pin(formula, expected, f"{stage}.formula_anchor")
    pin(len(product_set), formula, f"{stage}.three_route_count")
    q = p - 1
    pin(histogram.get(q, 0), 2 * q, f"{stage}.k0_valid_domain_class")
    pin(histogram.get(q - 1, 0), 4 * q * q, f"{stage}.k1_valid_domain_class")
    pin(histogram.get(q - 2, 0), q * q * (q - 1),
        f"{stage}.k2_valid_domain_class")
    return {
        "prime": p,
        "all_distinct_pairings": universe,
        "product_count": len(product_set),
        "cross_ratio_count": len(cr_set),
        "PGL_count": len(pgl_set),
        "formula": formula,
        "PGL_valid_domain_histogram": histogram,
        "set_equality": True,
    }


def planted_controls(matrices_by_p):
    p = 7
    q = 6
    true_size = p * (p * p - 1)
    wrong_k2 = q * q * (q - 2)
    wrong_class_sum = 2 * q + 4 * q * q + wrong_k2
    pin(wrong_class_sum, 300, "T4.wrong_class.planted_sum")
    check(wrong_class_sum != true_size, "T4.wrong_class.identity_reject",
          {"planted": wrong_class_sum, "PGL_size": true_size})
    wrong_formula = (
        2 * q * choose(q, 4)
        + 4 * q * q * choose(q - 1, 4)
        + wrong_k2 * choose(q - 2, 4)
    )
    pin(wrong_formula, 1044, "T4.wrong_class.planted_formula")
    pin(multiplicative_formula(p), 1080, "T4.wrong_class.true_formula")
    check(wrong_formula != multiplicative_formula(p), "T4.wrong_class.count_reject")

    reciprocal = canonical_matrix((0, 1, 1, 0), p)
    check(reciprocal in matrices_by_p[p], "T4.selected_pole.map_in_PGL")
    pole, zero = pole_zero(reciprocal, p)
    pin(pole, 0, "T4.selected_pole.at_zero")
    pin(mobius_apply(reciprocal, 0, p), p, "T4.selected_pole.maps_to_infinity")
    planted_u = {0, 1, 2, 3}
    valid_x = {
        x for x in range(p)
        if mobius_apply(reciprocal, x, p) != p
    }
    check(not planted_u.issubset(valid_x), "T4.selected_pole.U_rejected",
          {"U": sorted(planted_u), "valid_domain": sorted(valid_x)})

    points = tuple(range(p))
    factor_a = grs2("T4.tensor.A", points, p)
    factor_b = grs2("T4.tensor.B", points, p)
    columns = product_columns(factor_a, factor_b)
    identity = tuple((index, index) for index in (0, 1, 2, 3))
    pin(rank_cols([columns[cell] for cell in identity], p), 3,
        "T4.tensor.identity_accept")
    planted = {cell: list(column) for cell, column in columns.items()}
    planted[(0, 0)] = [
        (left + right) % p
        for left, right in zip(columns[(0, 0)], columns[(0, 1)])
    ]
    pin(rank_cols([planted[cell] for cell in identity], p), 4,
        "T4.tensor.non_tensor_reject")
    xs = [points[u] for u, _ in identity]
    ys = [points[v] for _, v in identity]
    pin(cross_ratio(xs, p), cross_ratio(ys, p), "T4.tensor.labels_still_accept")

    duplicate_rows = [row + [row[0]] for row in factor_b.rows]
    duplicate = Factor("T4.duplicate.B", duplicate_rows, p)
    pin(duplicate.spark, 2, "T4.duplicate.spark_drop")
    dependent_pairs = [
        pair
        for pair in itertools.combinations(range(duplicate.s), 2)
        if dependent([duplicate.cols[index] for index in pair], p)
    ]
    check(bool(dependent_pairs), "T4.duplicate.dependent_pair_detected")

    guard_a = grs2("T4.guard.A2", (1, 2, 4, 6), p)
    guard_points = (1, 2, 3, 5)
    guard_b = Factor(
        "T4.guard.B3",
        [[1] * 4, list(guard_points), [x * x % p for x in guard_points]],
        p,
    )
    kron = product_columns(guard_a, guard_b)
    independent = {
        (u, v): [
            guard_a.cols[u][i] * guard_b.cols[v][j] % p
            for i in range(guard_a.r)
            for j in range(guard_b.r)
        ]
        for u in range(guard_a.s)
        for v in range(guard_b.s)
    }
    pin(kron, independent, "T4.concat_guard.coordinates")
    pin({len(column) for column in kron.values()}, {6}, "T4.concat_guard.kron_dim")
    concat = {
        (u, v): guard_a.cols[u] + guard_b.cols[v]
        for u in range(guard_a.s)
        for v in range(guard_b.s)
    }
    pin({len(column) for column in concat.values()}, {5}, "T4.concat_guard.concat_dim")
    check(kron != concat, "T4.concat_guard.rejects_block_concat")

    RESULTS["T4_plants"] = {
        "wrong_class": {
            "q": q,
            "wrong_k2": wrong_k2,
            "wrong_class_sum": wrong_class_sum,
            "true_PGL_size": true_size,
            "wrong_formula": wrong_formula,
            "true_formula": multiplicative_formula(p),
            "rejected": True,
        },
        "selected_pole": {
            "map": reciprocal,
            "pole": pole,
            "zero": zero,
            "planted_U": sorted(planted_u),
            "valid_domain": sorted(valid_x),
            "rejected": True,
        },
        "non_tensor": {"pristine_rank": 3, "planted_rank": 4, "rejected": True},
        "duplicate": {"spark": duplicate.spark, "dependent_pairs": dependent_pairs},
        "concat": {"Kronecker_dim": 6, "concat_dim": 5, "rejected": True},
    }


def runtime_controls():
    pin(os.environ.get("PYTHONDONTWRITEBYTECODE"), "1", "runtime.no_bytecode_env")
    check(sys.dont_write_bytecode, "runtime.no_bytecode_interpreter")
    thread_vars = (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    )
    for variable in thread_vars:
        pin(os.environ.get(variable), "1", f"runtime.{variable}")
    if hasattr(os, "getpriority") and hasattr(os, "PRIO_PROCESS"):
        priority = os.getpriority(os.PRIO_PROCESS, 0)
        check(priority >= 15, "runtime.nice_at_least_15", {"nice": priority})
    else:
        priority = "platform-unavailable"
    RESULTS["runtime"] = {
        "bytecode_env": os.environ.get("PYTHONDONTWRITEBYTECODE"),
        "sys_dont_write_bytecode": sys.dont_write_bytecode,
        "threads": {variable: os.environ.get(variable) for variable in thread_vars},
        "nice": priority,
        "CPU_cap_seconds": CPU_CAP_SECONDS,
        "wall_cap_seconds": WALL_CAP_SECONDS,
    }


def main():
    runtime_controls()
    print(f"== {GATE} run {RUN_ID}")
    matrices_by_p = {}
    structure_rows = []
    for p in (5, 7, 11, 13):
        matrices, record = pgl_structure(p)
        matrices_by_p[p] = matrices
        structure_rows.append(record)
    RESULTS["T1_PGL_structure"] = structure_rows
    print("[T1] PGL sizes, pole/zero fibers, and k-class identities passed")

    full_rows = [
        full_field_controls(5, matrices_by_p[5], 200),
        full_field_controls(7, matrices_by_p[7], 5880),
    ]
    RESULTS["T2_full_field"] = full_rows
    print("[T2] full-field three-route set equality: p5=200, p7=5880")

    multiplicative_rows = [
        multiplicative_controls(5, matrices_by_p[5], 8),
        multiplicative_controls(7, matrices_by_p[7], 1080),
        multiplicative_controls(11, matrices_by_p[11], 117600),
    ]
    RESULTS["T3_multiplicative"] = multiplicative_rows
    print("[T3] multiplicative three-route set equality: 8/1080/117600")

    planted_controls(matrices_by_p)
    print("[T4] wrong coefficient, pole, non-tensor, duplicate, concat REJECTs passed")
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
