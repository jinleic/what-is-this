#!/usr/bin/env python3
"""Exact instrument for gate H-DP1-THREEFACTOR.

Standalone by design: never imports from a frozen campaign directory.
"""

import sys

sys.dont_write_bytecode = True

import itertools
import json
import math
import os
import resource
import time
import traceback
from collections import Counter
from pathlib import Path

P = 7
CPU_CAP = 600.0
WALL_CAP = 900.0
RUN_ID = "20260902T014115Z_25757c71_9c6065a4173c"
HERE = Path(__file__).resolve().parent
START_WALL = time.monotonic()
START_CPU = time.process_time()
CONTROLS = {}
ASSERTIONS = 0
RESULTS = {
    "run_id": RUN_ID,
    "gate": "H-DP1-THREEFACTOR",
    "arithmetic": "exact Python integers modulo 7; true nested row-major Kronecker columns",
}


def check(name, condition, detail=None):
    global ASSERTIONS
    ASSERTIONS += 1
    passed = bool(condition)
    CONTROLS[name] = {"pass": passed, "detail": detail}
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def rank_cols(cols, p=P):
    if not cols:
        return 0
    n = len(cols)
    dim = len(cols[0])
    if any(len(c) != dim for c in cols):
        raise ValueError("unequal column dimensions")
    a = [[cols[c][r] % p for c in range(n)] for r in range(dim)]
    row = 0
    for col in range(n):
        pivot = row
        while pivot < dim and a[pivot][col] == 0:
            pivot += 1
        if pivot == dim:
            continue
        if pivot != row:
            a[row], a[pivot] = a[pivot], a[row]
        inv = pow(a[row][col], -1, p)
        a[row] = [(x * inv) % p for x in a[row]]
        for rr in range(dim):
            if rr != row and a[rr][col]:
                scale = a[rr][col]
                a[rr] = [(a[rr][cc] - scale * a[row][cc]) % p for cc in range(n)]
        row += 1
        if row == dim:
            break
    return row


def is_circuit(cols, p=P):
    m = len(cols)
    if rank_cols(cols, p) != m - 1:
        return False
    return all(rank_cols(cols[:i] + cols[i + 1 :], p) == m - 1 for i in range(m))


def kron(a, b, p=P):
    return [(x * y) % p for x in a for y in b]


def vandermonde(name, rows, points, p=P):
    return {
        "name": name,
        "cols": [[pow(x, e, p) for e in range(rows)] for x in points],
        "rows": rows,
        "points": list(points),
    }


def measure_spark(factor, p=P):
    cols = factor["cols"]
    for k in range(1, len(cols) + 1):
        for inds in itertools.combinations(range(len(cols)), k):
            chosen = [cols[i] for i in inds]
            if rank_cols(chosen, p) < k:
                return k
    return len(cols) + 1


def factor_circuits(factor, m, p=P):
    out = set()
    for inds in itertools.combinations(range(len(factor["cols"])), m):
        if is_circuit([factor["cols"][i] for i in inds], p):
            out.add(tuple(inds))
    return out


def product_columns(factors, p=P):
    sizes = [len(f["cols"]) for f in factors]
    coords = []
    cols = []
    for cell in itertools.product(*(range(s) for s in sizes)):
        vec = [1]
        for axis, idx in enumerate(cell):
            vec = kron(vec, factors[axis]["cols"][idx], p)
        coords.append(tuple(cell))
        cols.append(vec)
    return coords, cols


def profile(support, n_axes):
    return tuple(len({cell[a] for cell in support}) for a in range(n_axes))


def circuit_census(factors, m, p=P):
    coords, cols = product_columns(factors, p)
    circuits = set()
    profiles = Counter()
    for inds in itertools.combinations(range(len(cols)), m):
        chosen = [cols[i] for i in inds]
        if not is_circuit(chosen, p):
            continue
        support = tuple(sorted(coords[i] for i in inds))
        circuits.add(support)
        profiles[profile(support, len(factors))] += 1
    return {
        "circuits": circuits,
        "profiles": profiles,
        "universe": math.comb(len(cols), m),
        "n_columns": len(cols),
        "ambient_dim": len(cols[0]),
    }


def pair_reference(left, right, m, p=P):
    census = circuit_census([left, right], m, p)
    nonfiber = {
        support
        for support in census["circuits"]
        if profile(support, 2)[0] > 1 and profile(support, 2)[1] > 1
    }
    fiber = census["circuits"] - nonfiber
    return {
        "all": census["circuits"],
        "nonfiber": nonfiber,
        "fiber": fiber,
        "profiles": census["profiles"],
        "universe": census["universe"],
    }


def construct_family(factors, m, include_pair_fibers=False, p=P):
    n = len(factors)
    multiplicity = Counter()
    components = Counter()
    factor_data = {}
    pair_data = {}

    for axis, factor in enumerate(factors):
        circuits = factor_circuits(factor, m, p)
        factor_data[str(axis)] = len(circuits)
        other_axes = [a for a in range(n) if a != axis]
        fixed_ranges = [range(len(factors[a]["cols"])) for a in other_axes]
        for circ in circuits:
            for fixed_values in itertools.product(*fixed_ranges):
                fixed = dict(zip(other_axes, fixed_values))
                support = []
                for idx in circ:
                    cell = tuple(idx if a == axis else fixed[a] for a in range(n))
                    support.append(cell)
                support = tuple(sorted(support))
                multiplicity[support] += 1
                components[f"single:{axis}"] += 1

    for i in range(n):
        for j in range(i + 1, n):
            ref = pair_reference(factors[i], factors[j], m, p)
            selected = ref["all"] if include_pair_fibers else ref["nonfiber"]
            pair_data[f"{i},{j}"] = {
                "all": len(ref["all"]),
                "nonfiber": len(ref["nonfiber"]),
                "fiber": len(ref["fiber"]),
                "profiles": {str(k): v for k, v in sorted(ref["profiles"].items())},
            }
            other_axes = [a for a in range(n) if a not in (i, j)]
            fixed_ranges = [range(len(factors[a]["cols"])) for a in other_axes]
            for pair_support in selected:
                for fixed_values in itertools.product(*fixed_ranges):
                    fixed = dict(zip(other_axes, fixed_values))
                    support = []
                    for left_idx, right_idx in pair_support:
                        cell = tuple(
                            left_idx if a == i else right_idx if a == j else fixed[a]
                            for a in range(n)
                        )
                        support.append(cell)
                    support = tuple(sorted(support))
                    multiplicity[support] += 1
                    components[f"pair:{i},{j}"] += 1

    return {
        "supports": set(multiplicity),
        "multiplicity": multiplicity,
        "components": components,
        "factor_circuits": factor_data,
        "pair_data": pair_data,
    }


def summarize_config(name, factors, measured, constructed):
    measured_set = measured["circuits"]
    constructed_set = constructed["supports"]
    check(f"{name}.support_set_equality", measured_set == constructed_set, {
        "measured_only": len(measured_set - constructed_set),
        "constructed_only": len(constructed_set - measured_set),
    })
    check(f"{name}.constructed_unique", all(v == 1 for v in constructed["multiplicity"].values()), {
        "max_multiplicity": max(constructed["multiplicity"].values(), default=0),
    })
    all_vary = sum(all(x > 1 for x in profile(s, len(factors))) for s in measured_set)
    return {
        "factor_names": [f["name"] for f in factors],
        "factor_sparks": [measure_spark(f) for f in factors],
        "factor_sizes": [len(f["cols"]) for f in factors],
        "n_columns": measured["n_columns"],
        "ambient_dim": measured["ambient_dim"],
        "universe": measured["universe"],
        "circuit_count": len(measured_set),
        "profiles": {str(k): v for k, v in sorted(measured["profiles"].items())},
        "all_coordinates_vary_count": all_vary,
        "constructed_components": dict(sorted(constructed["components"].items())),
        "factor_circuits_at_m": constructed["factor_circuits"],
        "pair_references": constructed["pair_data"],
        "set_equality": True,
    }


def write_failure(exc):
    cpu = time.process_time() - START_CPU
    wall = time.monotonic() - START_WALL
    payload = {
        "run_id": RUN_ID,
        "gate": RESULTS["gate"],
        "error_type": type(exc).__name__,
        "error": str(exc),
        "traceback": traceback.format_exc(),
        "partial_results": RESULTS,
        "controls": CONTROLS,
        "assertion_count_before_failure": ASSERTIONS,
        "budget": {"CPU_seconds": cpu, "wall_seconds": wall},
    }
    (HERE / "broken_results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run():
    print(f"== H-DP1-THREEFACTOR run {RUN_ID}", flush=True)

    check("runtime.dont_write_bytecode", sys.dont_write_bytecode is True, sys.dont_write_bytecode)
    check("runtime.env.PYTHONDONTWRITEBYTECODE", os.environ.get("PYTHONDONTWRITEBYTECODE") == "1", os.environ.get("PYTHONDONTWRITEBYTECODE"))
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        check(f"runtime.env.{var}", os.environ.get(var) == "1", os.environ.get(var))
    nice_value = os.getpriority(os.PRIO_PROCESS, 0)
    check("runtime.nice", nice_value == 15, nice_value)

    v3a = vandermonde("V3a", 2, [1, 2, 3])
    v3b = vandermonde("V3b", 2, [1, 2, 3])
    v3c = vandermonde("V3c", 2, [1, 2, 3])
    v4a = vandermonde("V4a", 2, [1, 2, 3, 4])
    v4b = vandermonde("V4b", 2, [1, 2, 3, 4])
    v4c = vandermonde("V4c", 2, [1, 2, 3, 4])
    m34 = vandermonde("M34", 3, [1, 2, 3, 4])

    for factor, want in ((v3a, 3), (v4a, 3), (m34, 4)):
        got = measure_spark(factor)
        check(f"factors.{factor['name']}.spark", got == want, {"got": got, "want": want})
    check("factors.M34.four_circuit_count", len(factor_circuits(m34, 4)) == 1, len(factor_circuits(m34, 4)))

    # T1: 3x V3, both spark-level and d+1 censuses.
    t1_factors = [v3a, v3b, v3c]
    t1_m3 = circuit_census(t1_factors, 3)
    check("T1.m3.universe", t1_m3["universe"] == 2925, t1_m3["universe"])
    check("T1.m3.count", len(t1_m3["circuits"]) == 27, len(t1_m3["circuits"]))
    expected_m3_profiles = {(1, 1, 3): 9, (1, 3, 1): 9, (3, 1, 1): 9}
    check("T1.m3.profiles", dict(t1_m3["profiles"]) == expected_m3_profiles, {str(k): v for k, v in t1_m3["profiles"].items()})
    check("T1.m3.all_fibers", all(sum(x > 1 for x in profile(s, 3)) == 1 for s in t1_m3["circuits"]), True)

    t1_m4 = circuit_census(t1_factors, 4)
    t1_constructed = construct_family(t1_factors, 4)
    t1_summary = summarize_config("T1.m4", t1_factors, t1_m4, t1_constructed)
    check("T1.m4.universe", t1_m4["universe"] == 17550, t1_m4["universe"])
    check("T1.m4.count", len(t1_m4["circuits"]) == 81, len(t1_m4["circuits"]))
    expected_t1_profiles = {(1, 3, 3): 27, (3, 1, 3): 27, (3, 3, 1): 27}
    check("T1.m4.profiles", dict(t1_m4["profiles"]) == expected_t1_profiles, {str(k): v for k, v in t1_m4["profiles"].items()})
    check("T1.m4.no_genuine_3D", t1_summary["all_coordinates_vary_count"] == 0, t1_summary["all_coordinates_vary_count"])
    check("T1.m4.pair_formula", sum(t1_constructed["components"].values()) == 3 * 3 * 9 == 81, dict(t1_constructed["components"]))
    for pair_key, pdata in t1_constructed["pair_data"].items():
        check(f"T1.m4.pair_{pair_key}_nonfiber", pdata["nonfiber"] == 9, pdata)
    RESULTS["T1_V3_cubed"] = {
        "m3": {"universe": t1_m3["universe"], "count": 27, "profiles": {str(k): v for k, v in sorted(t1_m3["profiles"].items())}},
        "m4": t1_summary,
    }
    print("[T1] V3^3 exact m3/m4 censuses: 27 fibers; 81 pairwise circuits", flush=True)

    # T2: 3x V4 anchor.
    t2_factors = [v4a, v4b, v4c]
    t2_m4 = circuit_census(t2_factors, 4)
    t2_constructed = construct_family(t2_factors, 4)
    t2_summary = summarize_config("T2.m4", t2_factors, t2_m4, t2_constructed)
    check("T2.m4.universe", t2_m4["universe"] == 635376, t2_m4["universe"])
    check("T2.m4.count", len(t2_m4["circuits"]) == 1824, len(t2_m4["circuits"]))
    expected_t2_profiles = {
        (1, 3, 3): 576,
        (1, 4, 4): 32,
        (3, 1, 3): 576,
        (3, 3, 1): 576,
        (4, 1, 4): 32,
        (4, 4, 1): 32,
    }
    check("T2.m4.profiles", dict(t2_m4["profiles"]) == expected_t2_profiles, {str(k): v for k, v in t2_m4["profiles"].items()})
    check("T2.m4.no_genuine_3D", t2_summary["all_coordinates_vary_count"] == 0, t2_summary["all_coordinates_vary_count"])
    check("T2.m4.formula", len(t2_m4["circuits"]) == 3 * 4 * 152, len(t2_m4["circuits"]))
    for pair_key, pdata in t2_constructed["pair_data"].items():
        check(f"T2.m4.pair_{pair_key}_total", pdata["nonfiber"] == 152, pdata)
        check(f"T2.m4.pair_{pair_key}_split", pdata["profiles"] == {"(3, 3)": 144, "(4, 4)": 8}, pdata)
    RESULTS["T2_V4_cubed"] = t2_summary
    print("[T2] V4^3 exact census: 1824 = 3 x 4 x (144+8); zero genuine 3D", flush=True)

    # T3: one size-4 factor circuit and unequal sparks.
    t3_factors = [m34, v4b, v4c]
    t3_m4 = circuit_census(t3_factors, 4)
    t3_constructed = construct_family(t3_factors, 4)
    t3_summary = summarize_config("T3.m4", t3_factors, t3_m4, t3_constructed)
    check("T3.m4.universe", t3_m4["universe"] == 635376, t3_m4["universe"])
    check("T3.m4.count", len(t3_m4["circuits"]) == 624, len(t3_m4["circuits"]))
    expected_t3_profiles = {(1, 3, 3): 576, (1, 4, 4): 32, (4, 1, 1): 16}
    check("T3.m4.profiles", dict(t3_m4["profiles"]) == expected_t3_profiles, {str(k): v for k, v in t3_m4["profiles"].items()})
    check("T3.m4.no_genuine_3D", t3_summary["all_coordinates_vary_count"] == 0, t3_summary["all_coordinates_vary_count"])
    check("T3.m4.single_factor_component", t3_constructed["components"]["single:0"] == 16, dict(t3_constructed["components"]))
    check("T3.m4.unequal_pair_01_no_nonfiber", t3_constructed["pair_data"]["0,1"]["nonfiber"] == 0, t3_constructed["pair_data"]["0,1"])
    check("T3.m4.unequal_pair_02_no_nonfiber", t3_constructed["pair_data"]["0,2"]["nonfiber"] == 0, t3_constructed["pair_data"]["0,2"])
    check("T3.m4.tied_pair_12", t3_constructed["pair_data"]["1,2"]["nonfiber"] == 152, t3_constructed["pair_data"]["1,2"])
    RESULTS["T3_M34_V4_V4"] = t3_summary
    print("[T3] M34 x V4 x V4 exact census: 16 fibers + 608 pair circuits = 624", flush=True)

    # T4: unequal sizes and sparks.
    t4_factors = [m34, v3b, v4c]
    t4_m4 = circuit_census(t4_factors, 4)
    t4_constructed = construct_family(t4_factors, 4)
    t4_summary = summarize_config("T4.m4", t4_factors, t4_m4, t4_constructed)
    check("T4.m4.universe", t4_m4["universe"] == 194580, t4_m4["universe"])
    check("T4.m4.count", len(t4_m4["circuits"]) == 156, len(t4_m4["circuits"]))
    expected_t4_profiles = {(1, 3, 3): 144, (4, 1, 1): 12}
    check("T4.m4.profiles", dict(t4_m4["profiles"]) == expected_t4_profiles, {str(k): v for k, v in t4_m4["profiles"].items()})
    check("T4.m4.no_genuine_3D", t4_summary["all_coordinates_vary_count"] == 0, t4_summary["all_coordinates_vary_count"])
    check("T4.m4.pair_12", t4_constructed["pair_data"]["1,2"]["nonfiber"] == 36, t4_constructed["pair_data"]["1,2"])
    check("T4.m4.formula", len(t4_m4["circuits"]) == 3 * 4 + 4 * 36, dict(t4_constructed["components"]))
    RESULTS["T4_M34_V3_V4"] = t4_summary
    print("[T4] M34 x V3 x V4 exact census: 12 fibers + 144 pair circuits = 156", flush=True)

    # T5 controls in both directions.
    known_support = tuple(sorted(((0, 1, 0), (0, 2, 0), (1, 0, 0), (2, 0, 0))))
    coords_t1, cols_t1 = product_columns(t1_factors)
    col_by_coord_t1 = dict(zip(coords_t1, cols_t1))
    known_cols = [col_by_coord_t1[c] for c in known_support]
    check("T5.known_crossing_accept.rank", rank_cols(known_cols) == 3, rank_cols(known_cols))
    check("T5.known_crossing_accept.circuit", is_circuit(known_cols), known_support)
    check("T5.known_crossing_accept.constructed", known_support in t1_constructed["supports"], known_support)

    corrupted = [list(c) for c in known_cols]
    extra = col_by_coord_t1[(0, 1, 1)]
    corrupted[0] = [(x + y) % P for x, y in zip(corrupted[0], extra)]
    check("T5.non_tensor_reject.rank4", rank_cols(corrupted) == 4, rank_cols(corrupted))
    check("T5.non_tensor_reject.labels_still_predicted", known_support in t1_constructed["supports"], known_support)
    check("T5.non_tensor_reject.not_circuit", not is_circuit(corrupted), rank_cols(corrupted))

    diagonal_support = tuple((i, i, i) for i in range(4))
    coords_t2, cols_t2 = product_columns(t2_factors)
    col_by_coord_t2 = dict(zip(coords_t2, cols_t2))
    diagonal_cols = [col_by_coord_t2[c] for c in diagonal_support]
    check("T5.genuine_3D_reject.profile", profile(diagonal_support, 3) == (4, 4, 4), profile(diagonal_support, 3))
    check("T5.genuine_3D_reject.rank4", rank_cols(diagonal_cols) == 4, rank_cols(diagonal_cols))
    check("T5.genuine_3D_reject.not_measured", diagonal_support not in t2_m4["circuits"], diagonal_support)
    check("T5.genuine_3D_reject.not_constructed", diagonal_support not in t2_constructed["supports"], diagonal_support)

    dup = {"name": "V4_duplicate", "cols": [list(c) for c in v4a["cols"]] + [list(v4a["cols"][0])], "rows": 2, "points": [1, 2, 3, 4, 1]}
    duplicate_spark = measure_spark(dup)
    check("T5.duplicate_reject.spark2", duplicate_spark == 2, duplicate_spark)
    dependent_pairs = [inds for inds in itertools.combinations(range(5), 2) if rank_cols([dup["cols"][i] for i in inds]) < 2]
    check("T5.duplicate_reject.pair", dependent_pairs == [(0, 4)], dependent_pairs)

    a, b, c = [1, 2], [3, 4], [5, 6]
    true_three = kron(kron(a, b), c)
    direct_three = [(x * y * z) % P for x in a for y in b for z in c]
    block_concat = a + b + c
    dropped_third = kron(a, b)
    check("T5.kronecker_guard.true_dim", len(true_three) == 8, len(true_three))
    check("T5.kronecker_guard.coordinates", true_three == direct_three, {"nested": true_three, "direct": direct_three})
    check("T5.kronecker_guard.concat_reject", len(block_concat) == 6 and block_concat != true_three, {"concat": block_concat, "true": true_three})
    check("T5.kronecker_guard.dropped_reject", len(dropped_third) == 4 and dropped_third != true_three, {"dropped": dropped_third, "true": true_three})

    wrong_t3 = construct_family(t3_factors, 4, include_pair_fibers=True)
    raw_wrong = sum(wrong_t3["multiplicity"].values())
    duplicate_excess = raw_wrong - len(wrong_t3["supports"])
    triple_count = sum(v == 3 for v in wrong_t3["multiplicity"].values())
    check("T5.formula_semantics_reject.raw_count", raw_wrong == 656, raw_wrong)
    check("T5.formula_semantics_reject.unique_count", len(wrong_t3["supports"]) == 624, len(wrong_t3["supports"]))
    check("T5.formula_semantics_reject.duplicate_excess", duplicate_excess == 32, duplicate_excess)
    check("T5.formula_semantics_reject.triple_count", triple_count == 16, triple_count)

    RESULTS["T5_plants"] = {
        "known_crossing": {"support": known_support, "rank": 3, "accepted": True},
        "non_tensor": {"labeled_support": known_support, "actual_rank": 4, "rejected": True},
        "genuine_3D": {"support": diagonal_support, "profile": (4, 4, 4), "rank": 4, "rejected": True},
        "duplicate": {"spark": duplicate_spark, "dependent_pairs": dependent_pairs},
        "kronecker_guard": {"true_dim": 8, "concat_dim": 6, "dropped_dim": 4},
        "formula_semantics": {"correct_count": 624, "wrong_raw_count": raw_wrong, "wrong_unique_count": len(wrong_t3["supports"]), "duplicate_excess": duplicate_excess, "supports_with_multiplicity_3": triple_count},
    }
    print("[T5] ACCEPT and all planted non-tensor/3D/duplicate/concat/count REJECTs passed", flush=True)

    cpu = time.process_time() - START_CPU
    wall = time.monotonic() - START_WALL
    check("budget.CPU", cpu <= CPU_CAP, {"used": cpu, "cap": CPU_CAP})
    check("budget.wall", wall <= WALL_CAP, {"used": wall, "cap": WALL_CAP})
    RESULTS["controls"] = CONTROLS
    RESULTS["assertion_count"] = ASSERTIONS
    RESULTS["budget"] = {
        "CPU_cap_seconds": CPU_CAP,
        "CPU_seconds": cpu,
        "wall_cap_seconds": WALL_CAP,
        "wall_seconds": wall,
        "within_cap": True,
    }
    RESULTS["runtime"] = {
        "nice": nice_value,
        "thread_caps": {var: os.environ.get(var) for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
        "dont_write_bytecode": sys.dont_write_bytecode,
        "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
        "max_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (HERE / "controls_results.json").write_text(json.dumps(RESULTS, indent=2, sort_keys=True) + "\n")
    print(f"== ALL CONTROLS PASS: {ASSERTIONS} aggregate asserts; CPU {cpu:.3f}s, wall {wall:.3f}s; controls_results.json", flush=True)


if __name__ == "__main__":
    try:
        run()
    except Exception as exc:
        write_failure(exc)
        print(f"FAILED: {type(exc).__name__}: {exc}; preserved broken_results.json", file=sys.stderr, flush=True)
        raise
