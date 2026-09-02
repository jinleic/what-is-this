#!/usr/bin/env python3
"""Exact controls for H-ALLDISTINCT-CROSSRATIO.

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

RUN_ID = "20260902T011204Z_7ece99fe_a108b246f87d"
GATE = "H-ALLDISTINCT-CROSSRATIO"
PREREG_COMMIT = "44a491993b229ecf1789a9f7d64b4b8130e2dcb9"
CPU_CAP_SECONDS = 600.0
WALL_CAP_SECONDS = 900.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT_DIR = Path(__file__).resolve().parent
ASSERTION_COUNT = 0
RESULTS = {
    "schema": "h-alldistinct-crossratio-controls/1",
    "run_id": RUN_ID,
    "gate": GATE,
    "prereg_commit": PREREG_COMMIT,
    "arithmetic": "exact Python integers modulo p; no numerical algebra package",
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
        raise AssertionError("determinant matrix not square")
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


def null_vector(matrix, p):
    reduced, pivots = rref_with_pivots(matrix, p)
    ncols = len(matrix[0])
    free = [col for col in range(ncols) if col not in pivots]
    if not free:
        return None
    free_col = free[-1]
    vector = [0] * ncols
    vector[free_col] = 1
    for row, pivot_col in enumerate(pivots):
        vector[pivot_col] = (-reduced[row][free_col]) % p
    if not any(vector):
        raise AssertionError("null solver returned zero vector")
    if any(sum(row[j] * vector[j] for j in range(ncols)) % p for row in matrix):
        raise AssertionError("null solver produced non-null vector")
    return vector


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
    ordered = tuple(sorted(support))
    return [
        [columns[key][row] for key in ordered]
        for row in range(len(columns[ordered[0]]))
    ]


def is_circuit(columns, support, p):
    ordered = tuple(sorted(support))
    if rank_cols([columns[key] for key in ordered], p) != len(ordered) - 1:
        return False
    return all(
        rank_cols([columns[key] for key in proper], p) == len(ordered) - 1
        for proper in itertools.combinations(ordered, len(ordered) - 1)
    )


def cross_ratio(values, p):
    a, b, c, d = (value % p for value in values)
    numerator = (a - c) * (b - d) % p
    denominator = (a - d) * (b - c) % p
    if denominator == 0:
        raise AssertionError(f"cross-ratio denominator zero for {values} mod {p}")
    return numerator * pow(denominator, -1, p) % p


def canonical_pairing(support):
    ordered = sorted(support)
    if len({u for u, _ in ordered}) != len(ordered):
        raise AssertionError("support is not a matching on A indices")
    return ordered


def cr_values(factor_a, factor_b, support):
    ordered = canonical_pairing(support)
    xs = [factor_a.cols[u][1] for u, _ in ordered]
    ys = [factor_b.cols[v][1] for _, v in ordered]
    return cross_ratio(xs, factor_a.p), cross_ratio(ys, factor_a.p)


def all_distinct_supports(factor_a, factor_b):
    for rows in itertools.combinations(range(factor_a.s), 4):
        for columns in itertools.combinations(range(factor_b.s), 4):
            for permutation in itertools.permutations(columns):
                yield frozenset(zip(rows, permutation))


def bilinear_witness(factor_a, factor_b, support):
    p = factor_a.p
    points = [
        (factor_a.cols[u][1], factor_b.cols[v][1])
        for u, v in canonical_pairing(support)
    ]
    evaluation = [[1, y, x, x * y % p] for x, y in points]
    coefficients = null_vector(evaluation, p)
    if coefficients is None:
        raise AssertionError("det-zero support has no bilinear left-null vector")
    c0, c1, c2, c3 = coefficients
    values = [
        (c0 + c1 * y + c2 * x + c3 * x * y) % p for x, y in points
    ]
    if any(values):
        raise AssertionError(f"bilinear witness does not vanish: {values}")
    delta = (c0 * c3 - c1 * c2) % p
    if delta == 0:
        raise AssertionError(
            f"accepted all-distinct support has degenerate Möbius determinant: {coefficients}"
        )
    denominators = [(c3 * x + c1) % p for x, _ in points]
    if any(value == 0 for value in denominators):
        raise AssertionError(
            f"accepted all-distinct support hits a Möbius pole: {denominators}"
        )
    predicted = [
        (-(c2 * x + c0) * pow(c3 * x + c1, -1, p)) % p
        for x, _ in points
    ]
    if predicted != [y % p for _, y in points]:
        raise AssertionError("Möbius witness does not reproduce y values")
    if c3 == 0:
        if c1 == 0 or c2 == 0:
            raise AssertionError("c3=0 accepted curve is axis-parallel/degenerate")
        curve = {
            "kind": "line",
            "slope": (-c2 * pow(c1, -1, p)) % p,
            "intercept": (-c0 * pow(c1, -1, p)) % p,
        }
    else:
        inverse_c3 = pow(c3, -1, p)
        shift_x = c1 * inverse_c3 % p
        shift_y = c2 * inverse_c3 % p
        r = (c1 * c2 - c0 * c3) * pow(c3 * c3 % p, -1, p) % p
        if r == 0:
            raise AssertionError("accepted hyperbola has r=0")
        if any((x + shift_x) % p == 0 for x, _ in points):
            raise AssertionError("accepted hyperbola contains selected pole")
        if any(((x + shift_x) * (y + shift_y) - r) % p for x, y in points):
            raise AssertionError("hyperbola factorization check failed")
        curve = {
            "kind": "hyperbola",
            "shift_x": shift_x,
            "shift_y": shift_y,
            "r": r,
        }
    return {
        "coefficients": coefficients,
        "mobius_determinant": delta,
        "denominators": denominators,
        "curve": curve,
    }


def census_size_four(factor_a, factor_b):
    p = factor_a.p
    columns = product_columns(factor_a, factor_b)
    keys = tuple(sorted(columns))
    circuits = set()
    profiles = {}
    scanned = 0
    for support in itertools.combinations(keys, 4):
        scanned += 1
        if scanned % 256 == 0:
            check_budget(f"census.{factor_a.name}.{factor_b.name}.{scanned}")
        if rank_cols([columns[key] for key in support], p) >= 4:
            continue
        if not all(
            rank_cols([columns[key] for key in proper], p) == 3
            for proper in itertools.combinations(support, 3)
        ):
            continue
        frozen = frozenset(support)
        circuits.add(frozen)
        prof = profile(frozen)
        profiles.setdefault(prof, set()).add(frozen)
    return {
        "scanned": scanned,
        "circuits": circuits,
        "profiles": profiles,
    }


def evaluate_configuration(name, x_points, y_points, p, expected_all=None,
                           expected_total=None, expected_crossing=None):
    check_budget(f"{name}.start")
    factor_a = grs2(f"{name}.A", x_points, p)
    factor_b = grs2(f"{name}.B", y_points, p)
    pin(factor_a.spark, 3, f"{name}.spark_A")
    pin(factor_b.spark, 3, f"{name}.spark_B")
    columns = product_columns(factor_a, factor_b)
    expected_columns = {
        (u, v): [1, y_points[v] % p, x_points[u] % p,
                 x_points[u] * y_points[v] % p]
        for u in range(len(x_points))
        for v in range(len(y_points))
    }
    pin(columns, expected_columns, f"{name}.pure_tensor_coordinates")

    universe = list(all_distinct_supports(factor_a, factor_b))
    universe_expected = (
        math.comb(factor_a.s, 4) * math.comb(factor_b.s, 4) * math.factorial(4)
    )
    pin(len(universe), universe_expected, f"{name}.all_pairings_enumerated")
    pin(len(set(universe)), universe_expected, f"{name}.pairings_unique")

    det_zero = set()
    cr_equal = set()
    witness_rows = {}
    curve_counts = {"line": 0, "hyperbola": 0}
    for index, support in enumerate(universe):
        if index % 128 == 0:
            check_budget(f"{name}.pairing{index}")
        matrix = support_matrix(columns, support)
        det = determinant(matrix, p)
        cr_x, cr_y = cr_values(factor_a, factor_b, support)
        if det == 0:
            det_zero.add(support)
            if not is_circuit(columns, support, p):
                raise AssertionError(f"{name}: det-zero all-distinct support not circuit")
            witness = bilinear_witness(factor_a, factor_b, support)
            curve_counts[witness["curve"]["kind"]] += 1
            witness_rows[tuple(sorted(support))] = witness
        if cr_x == cr_y:
            cr_equal.add(support)

    pin(det_zero, cr_equal, f"{name}.det_zero_equals_cross_ratio")
    census = census_size_four(factor_a, factor_b)
    measured44 = census["profiles"].get((4, 4), set())
    pin(measured44, det_zero, f"{name}.census44_equals_determinant")
    pin(measured44, cr_equal, f"{name}.census44_equals_cross_ratio")
    if expected_all is not None:
        pin(len(det_zero), expected_all, f"{name}.all_distinct_count")
    if expected_total is not None:
        pin(len(census["circuits"]), expected_total, f"{name}.total_circuit_count")
    measured33 = census["profiles"].get((3, 3), set())
    if expected_crossing is not None:
        pin(len(measured33), expected_crossing, f"{name}.profile33_count")
    other_profiles = {
        prof: len(supports)
        for prof, supports in census["profiles"].items()
        if prof not in ((3, 3), (4, 4))
    }
    pin(other_profiles, {}, f"{name}.no_other_circuit_profiles")

    return {
        "factor_a": factor_a,
        "factor_b": factor_b,
        "columns": columns,
        "det_zero": det_zero,
        "cr_equal": cr_equal,
        "witness_rows": witness_rows,
        "record": {
            "name": name,
            "prime": p,
            "x_points": list(x_points),
            "y_points": list(y_points),
            "all_distinct_pairings": len(universe),
            "det_zero": len(det_zero),
            "cross_ratio_equal": len(cr_equal),
            "census_profile44": len(measured44),
            "census_profile33": len(measured33),
            "census_total_circuits": len(census["circuits"]),
            "census_subsets_scanned": census["scanned"],
            "curve_counts": curve_counts,
            "set_equality": True,
            "zero_disagreements": True,
        },
    }


def symmetric_sweep():
    counts = {7: 8, 11: 4, 13: 12, 17: 4, 31: 4}
    totals = {7: 152, 11: 148, 13: 156, 17: 148, 31: 148}
    rows = []
    p13_data = None
    for p in (7, 11, 13, 17, 31):
        data = evaluate_configuration(
            f"T1.p{p}.symmetric4",
            (1, 2, 3, 4),
            (1, 2, 3, 4),
            p,
            expected_all=counts[p],
            expected_total=totals[p],
            expected_crossing=144,
        )
        identity = frozenset((i, i) for i in range(4))
        reversal = frozenset((i, 3 - i) for i in range(4))
        check(identity in data["det_zero"], f"T1.p{p}.identity_accept")
        check(reversal in data["det_zero"], f"T1.p{p}.reversal_accept")
        identity_witness = data["witness_rows"][tuple(sorted(identity))]
        reversal_witness = data["witness_rows"][tuple(sorted(reversal))]
        pin(identity_witness["curve"]["kind"], "line", f"T1.p{p}.identity_line")
        pin(reversal_witness["curve"]["kind"], "line", f"T1.p{p}.reversal_line")
        rows.append(data["record"])
        if p == 13:
            p13_data = data
    hyperbola = frozenset({(0, 0), (1, 2), (2, 3), (3, 1)})
    check(hyperbola in p13_data["det_zero"], "T1.p13.hyperbola_accept")
    hyperbola_witness = p13_data["witness_rows"][tuple(sorted(hyperbola))]
    pin(hyperbola_witness["curve"]["kind"], "hyperbola",
        "T1.p13.hyperbola_not_line")
    check(hyperbola_witness["curve"]["r"] != 0, "T1.p13.hyperbola_r_nonzero")
    RESULTS["T1_symmetric_five_prime"] = {
        "rows": rows,
        "known_accepts": {
            "identity": sorted(frozenset((i, i) for i in range(4))),
            "reversal": sorted(frozenset((i, 3 - i) for i in range(4))),
            "GF13_hyperbola": sorted(hyperbola),
            "GF13_hyperbola_witness": hyperbola_witness,
        },
        "det_CR_agreements": 24 * 5,
        "all_set_equalities": True,
    }
    return p13_data


def asymmetric_sweep():
    rows = []
    for p in (13, 17):
        data = evaluate_configuration(
            f"T2.p{p}.asymmetric5",
            (1, 2, 3, 4, 5),
            (2, 3, 5, 7, 11),
            p,
            expected_crossing=900,
        )
        pin(data["record"]["all_distinct_pairings"], 600,
            f"T2.p{p}.exact_600_pairings")
        rows.append(data["record"])
    RESULTS["T2_asymmetric_two_prime"] = {
        "rows": rows,
        "pairings_checked_total": 1200,
        "all_set_equalities": True,
        "zero_disagreements": True,
    }


def planted_rejects_and_layout(p13_data):
    p = 13
    factor_a = p13_data["factor_a"]
    factor_b = p13_data["factor_b"]
    columns = p13_data["columns"]

    mismatch = frozenset({(0, 0), (1, 1), (2, 3), (3, 2)})
    cr_x, cr_y = cr_values(factor_a, factor_b, mismatch)
    mismatch_det = determinant(support_matrix(columns, mismatch), p)
    check(cr_x != cr_y, "T3.cross_ratio_mismatch_reject",
          {"CR_x": cr_x, "CR_y": cr_y})
    check(mismatch_det != 0, "T3.cross_ratio_mismatch_det_nonzero",
          {"determinant": mismatch_det})
    check(not is_circuit(columns, mismatch, p), "T3.cross_ratio_mismatch_not_circuit")

    identity = frozenset((i, i) for i in range(4))
    pin(rank_cols([columns[key] for key in sorted(identity)], p), 3,
        "T3.nontensor_pristine_identity_rank")
    planted = {key: list(column) for key, column in columns.items()}
    planted[(0, 0)] = [
        (a + b) % p for a, b in zip(columns[(0, 0)], columns[(0, 1)])
    ]
    planted_rank = rank_cols([planted[key] for key in sorted(identity)], p)
    pin(planted_rank, 4, "T3.nontensor_plant_rejected")
    identity_cr = cr_values(factor_a, factor_b, identity)
    pin(identity_cr[0], identity_cr[1], "T3.nontensor_label_CR_still_accepts")
    check(
        determinant(support_matrix(planted, identity), p) != 0,
        "T3.nontensor_actual_determinant_disagrees",
    )

    duplicate_rows = [row + [row[0]] for row in factor_b.rows]
    duplicate = Factor("T3.duplicate.B", duplicate_rows, p)
    pin(duplicate.spark, 2, "T3.duplicate.spark_drop")
    dependent_pairs = [
        pair
        for pair in itertools.combinations(range(duplicate.s), 2)
        if dependent([duplicate.cols[index] for index in pair], p)
    ]
    check(bool(dependent_pairs), "T3.duplicate.dependent_pair_detected")

    guard_a = grs2("T3.guard.A2", (1, 2, 4, 7), p)
    guard_points = (1, 2, 3, 5)
    guard_b = Factor(
        "T3.guard.B3",
        [[1] * 4, list(guard_points), [x * x % p for x in guard_points]],
        p,
    )
    kron = product_columns(guard_a, guard_b)
    independently_expected = {
        (u, v): [
            guard_a.cols[u][i] * guard_b.cols[v][j] % p
            for i in range(guard_a.r)
            for j in range(guard_b.r)
        ]
        for u in range(guard_a.s)
        for v in range(guard_b.s)
    }
    pin(kron, independently_expected, "T3.concat_guard_exact_coordinates")
    pin({len(column) for column in kron.values()}, {6},
        "T3.concat_guard_kron_dimension")
    concat = {
        (u, v): guard_a.cols[u] + guard_b.cols[v]
        for u in range(guard_a.s)
        for v in range(guard_b.s)
    }
    pin({len(column) for column in concat.values()}, {5},
        "T3.concat_guard_buggy_dimension")
    check(kron != concat, "T3.concat_guard_rejects_buggy_coordinates")

    RESULTS["T3_rejects_and_layout"] = {
        "cross_ratio_mismatch": {
            "support": mismatch,
            "CR_x": cr_x,
            "CR_y": cr_y,
            "determinant": mismatch_det,
            "rejected": True,
        },
        "non_tensor_plant": {
            "modified_column": [0, 0],
            "added_column": [0, 1],
            "pristine_rank": 3,
            "planted_rank": planted_rank,
            "label_CR_equal": True,
            "rejected": True,
        },
        "duplicated_factor_column": {
            "original_spark": factor_b.spark,
            "planted_spark": duplicate.spark,
            "dependent_pairs": dependent_pairs,
            "rejected": True,
        },
        "concat_guard": {
            "Kronecker_dimension": 6,
            "concat_dimension": 5,
            "sample_Kronecker": kron[(1, 2)],
            "sample_concat": concat[(1, 2)],
            "rejected": True,
        },
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
        "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
        "sys_dont_write_bytecode": sys.dont_write_bytecode,
        "threads": {variable: os.environ.get(variable) for variable in thread_vars},
        "nice": priority,
        "CPU_cap_seconds": CPU_CAP_SECONDS,
        "wall_cap_seconds": WALL_CAP_SECONDS,
    }


def main():
    runtime_controls()
    print(f"== {GATE} run {RUN_ID}")
    p13_data = symmetric_sweep()
    print("[T1] symmetric 4x4: five-prime det/CR/census set equality and 8/4/12/4/4 passed")
    asymmetric_sweep()
    print("[T2] asymmetric 5x5: 600 pairings per prime, det/CR/census set equality passed")
    planted_rejects_and_layout(p13_data)
    print("[T3] mismatch, non-tensor, duplicate, and concat REJECT controls passed")
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
