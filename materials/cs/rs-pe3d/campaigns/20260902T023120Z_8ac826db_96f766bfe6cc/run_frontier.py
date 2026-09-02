#!/usr/bin/env python3
"""Exact controls for H-2ROW-COMPLETE-3ROW-OPEN; no frozen imports."""

import sys

sys.dont_write_bytecode = True

import functools
import itertools
import json
import math
import os
import resource
import time
import traceback
from collections import Counter, defaultdict, deque
from pathlib import Path

CPU_CAP = 600.0
WALL_CAP = 900.0
RUN_ID = "20260902T023120Z_8ac826db_96f766bfe6cc"
HERE = Path(__file__).resolve().parent
START_CPU = time.process_time()
START_WALL = time.monotonic()
CONTROLS = {}
ASSERTIONS = 0
RESULTS = {
    "run_id": RUN_ID,
    "gate": "H-2ROW-COMPLETE-3ROW-OPEN",
    "arithmetic": "exact Python integers modulo p; true row-major Kronecker columns",
}


def check(name, condition, detail=None):
    global ASSERTIONS
    ASSERTIONS += 1
    passed = bool(condition)
    CONTROLS[name] = {"pass": passed, "detail": detail}
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def rank_cols(cols, p):
    if not cols:
        return 0
    n = len(cols)
    dim = len(cols[0])
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


def kron(a, b, p):
    return [(x * y) % p for x in a for y in b]


def mask_from_indices(indices):
    mask = 0
    for i in indices:
        mask |= 1 << i
    return mask


def mask_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def mask_support(mask, cells):
    return tuple(cells[i] for i in mask_indices(mask))


def support_mask(support, cell_index):
    return mask_from_indices(cell_index[cell] for cell in support)


def make_rank_oracle(cols, p):
    @functools.lru_cache(maxsize=None)
    def rank_mask(mask):
        return rank_cols([cols[i] for i in mask_indices(mask)], p)

    return rank_mask


def is_circuit_mask(mask, m, rank_mask):
    if rank_mask(mask) != m - 1:
        return False
    return all(rank_mask(mask ^ (1 << i)) == m - 1 for i in mask_indices(mask))


def profile(support):
    return (len({i for i, _ in support}), len({j for _, j in support}))


def degrees(support):
    return Counter(i for i, _ in support), Counter(j for _, j in support)


def component_signature(support):
    adjacency = defaultdict(set)
    for i, j in support:
        r, c = ("A", i), ("B", j)
        adjacency[r].add(c)
        adjacency[c].add(r)
    seen = set()
    comps = []
    for start in adjacency:
        if start in seen:
            continue
        queue = deque([start])
        nodes = set()
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            nodes.add(node)
            queue.extend(adjacency[node] - seen)
        rows = {v for side, v in nodes if side == "A"}
        cols = {v for side, v in nodes if side == "B"}
        edges = [(i, j) for i, j in support if i in rows and j in cols]
        row_degs = sorted((sum(1 for i2, _ in edges if i2 == i) for i in rows), reverse=True)
        col_degs = sorted((sum(1 for _, j2 in edges if j2 == j) for j in cols), reverse=True)
        cycle_rank = len(edges) - len(rows) - len(cols) + 1
        comps.append((len(edges), len(rows), len(cols), cycle_rank, tuple(row_degs), tuple(col_degs)))
    comps.sort(reverse=True)
    dr, dc = degrees(support)
    return str((profile(support), tuple(sorted(dr.values(), reverse=True)), tuple(sorted(dc.values(), reverse=True)), tuple(comps)))


def unequal_columns(n, p):
    points = list(range(1, n + 1))
    cells = [(i, j) for i in range(n) for j in range(n)]
    cols = []
    for i, j in cells:
        x, y = points[i] % p, points[j] % p
        cols.append([1, y, y * y % p, x, x * y % p, x * y * y % p])
    return points, cells, cols


def two_row_columns(n, p):
    points = list(range(1, n + 1))
    cells = [(i, j) for i in range(n) for j in range(n)]
    cols = []
    for i, j in cells:
        x, y = points[i] % p, points[j] % p
        cols.append([1, y, x, x * y % p])
    return points, cells, cols


def construct_size3_fibers(n, cell_index):
    out = set()
    for rows in itertools.combinations(range(n), 3):
        for fixed_b in range(n):
            out.add(support_mask(tuple((i, fixed_b) for i in rows), cell_index))
    return out


def construct_size4_fibers(n, cell_index):
    out = set()
    for fixed_a in range(n):
        for cols in itertools.combinations(range(n), 4):
            out.add(support_mask(tuple((fixed_a, j) for j in cols), cell_index))
    return out


def construct_unequal_crossings(n, cell_index):
    out = set()
    multiplicity = Counter()
    for R in itertools.combinations(range(n), 3):
        for Z in itertools.combinations(range(n), 4):
            for i0 in R:
                for j0 in Z:
                    support = tuple(sorted(
                        {(i0, j) for j in Z if j != j0}
                        | {(i, j0) for i in R if i != i0}
                    ))
                    mask = support_mask(support, cell_index)
                    multiplicity[mask] += 1
                    out.add(mask)
    return out, multiplicity


def construct_fieldfree_size6(n, cell_index):
    out34 = set()
    out44 = set()
    out54 = set()
    for rows in itertools.combinations(range(n), 3):
        for cols in itertools.combinations(range(n), 4):
            for c4_rows in itertools.combinations(rows, 2):
                remaining_row = next(r for r in rows if r not in c4_rows)
                for c4_cols in itertools.combinations(cols, 2):
                    remaining_cols = [c for c in cols if c not in c4_cols]
                    support = {(r, c) for r in c4_rows for c in c4_cols}
                    support |= {(remaining_row, c) for c in remaining_cols}
                    out34.add(support_mask(tuple(sorted(support)), cell_index))
    for rows in itertools.combinations(range(n), 4):
        for cols in itertools.combinations(range(n), 4):
            for p3_row in rows:
                p5_rows = [r for r in rows if r != p3_row]
                for p3_cols in itertools.combinations(cols, 2):
                    p5_cols = [c for c in cols if c not in p3_cols]
                    for central_row in p5_rows:
                        leaf_rows = [r for r in p5_rows if r != central_row]
                        for perm in itertools.permutations(p5_cols):
                            support = {(p3_row, c) for c in p3_cols}
                            support |= {(central_row, c) for c in p5_cols}
                            support |= set(zip(leaf_rows, perm))
                            out44.add(support_mask(tuple(sorted(support)), cell_index))
    for rows in itertools.combinations(range(n), 5):
        for cols in itertools.combinations(range(n), 4):
            for row_center in rows:
                remaining_rows = [r for r in rows if r != row_center]
                for row_center_cols in itertools.combinations(cols, 2):
                    col_centers = sorted(c for c in cols if c not in row_center_cols)
                    for first_pair in itertools.combinations(remaining_rows, 2):
                        second_pair = [r for r in remaining_rows if r not in first_pair]
                        support = {(row_center, c) for c in row_center_cols}
                        support |= {(r, col_centers[0]) for r in first_pair}
                        support |= {(r, col_centers[1]) for r in second_pair}
                        out54.add(support_mask(tuple(sorted(support)), cell_index))
    return {(3, 4): out34, (4, 4): out44, (5, 4): out54}


def construct_fieldfree_size7_34(n, cell_index):
    out = set()
    for rows in itertools.combinations(range(n), 3):
        for cols in itertools.combinations(range(n), 4):
            for isolated_row in rows:
                block_rows = [r for r in rows if r != isolated_row]
                for isolated_col in cols:
                    block_cols = [c for c in cols if c != isolated_col]
                    support = {(r, c) for r in block_rows for c in block_cols}
                    support.add((isolated_row, isolated_col))
                    out.add(support_mask(tuple(sorted(support)), cell_index))
    return out


def unequal_census(n, p):
    points, cells, cols = unequal_columns(n, p)
    cell_index = {cell: i for i, cell in enumerate(cells)}
    rank_mask = make_rank_oracle(cols, p)
    circuits = {}
    profiles = {}
    signatures = {}
    degree_violations = {}
    for m in range(3, 8):
        cset = set()
        pc = Counter()
        sc = Counter()
        violations = 0
        for inds in itertools.combinations(range(n * n), m):
            mask = mask_from_indices(inds)
            if not is_circuit_mask(mask, m, rank_mask):
                continue
            support = tuple(cells[i] for i in inds)
            cset.add(mask)
            pc[profile(support)] += 1
            sc[component_signature(support)] += 1
            if m >= 5:
                dr, dc = degrees(support)
                if max(dr.values()) > 3 or max(dc.values()) > 2:
                    violations += 1
        circuits[m] = cset
        profiles[m] = pc
        signatures[m] = sc
        degree_violations[m] = violations
    return {
        "n": n,
        "p": p,
        "points": points,
        "cells": cells,
        "cols": cols,
        "cell_index": cell_index,
        "rank_mask": rank_mask,
        "circuits": circuits,
        "profiles": profiles,
        "signatures": signatures,
        "degree_violations": degree_violations,
    }


def two_row_size6_check(n, p):
    points, cells, cols = two_row_columns(n, p)
    rank_mask = make_rank_oracle(cols, p)
    total = 0
    circuits = 0
    all_have_dep5 = True
    for inds in itertools.combinations(range(n * n), 6):
        total += 1
        mask = mask_from_indices(inds)
        if is_circuit_mask(mask, 6, rank_mask):
            circuits += 1
        if not any(rank_mask(mask ^ (1 << i)) < 5 for i in inds):
            all_have_dep5 = False
    return {"n": n, "p": p, "universe": total, "circuits": circuits, "all_have_dependent_5subset": all_have_dep5}


def measure_spark(cols, p):
    for k in range(1, len(cols) + 1):
        for inds in itertools.combinations(range(len(cols)), k):
            if rank_cols([cols[i] for i in inds], p) < k:
                return k
    return len(cols) + 1


def summarize_unequal(data):
    n, p = data["n"], data["p"]
    label = f"unequal.n{n}.p{p}"
    cell_index = data["cell_index"]
    size3 = construct_size3_fibers(n, cell_index)
    size4 = construct_size4_fibers(n, cell_index)
    crossing5, crossing_mult = construct_unequal_crossings(n, cell_index)
    check(f"{label}.size3_fiber_set_equality", data["circuits"][3] == size3, {"measured": len(data["circuits"][3]), "constructed": len(size3)})
    check(f"{label}.size3_formula", len(size3) == n * math.comb(n, 3), len(size3))
    check(f"{label}.size4_fiber_set_equality", data["circuits"][4] == size4, {"measured": len(data["circuits"][4]), "constructed": len(size4)})
    check(f"{label}.size4_formula", len(size4) == n * math.comb(n, 4), len(size4))
    check(f"{label}.crossing_parameter_injective", all(v == 1 for v in crossing_mult.values()), max(crossing_mult.values(), default=0))
    measured34 = {mask for mask in data["circuits"][5] if profile(mask_support(mask, data["cells"])) == (3, 4)}
    check(f"{label}.size5_crossing_set_equality", measured34 == crossing5, {"measured": len(measured34), "constructed": len(crossing5)})
    want_crossing = 12 * math.comb(n, 3) * math.comb(n, 4)
    check(f"{label}.size5_crossing_formula", len(crossing5) == want_crossing, {"got": len(crossing5), "want": want_crossing})
    fieldfree6 = construct_fieldfree_size6(n, cell_index)
    formulas6 = {
        (3, 4): 18 * math.comb(n, 3) * math.comb(n, 4),
        (4, 4): 144 * math.comb(n, 4) ** 2,
        (5, 4): 180 * math.comb(n, 5) * math.comb(n, 4),
    }
    for pr, constructed in fieldfree6.items():
        measured = {mask for mask in data["circuits"][6] if profile(mask_support(mask, data["cells"])) == pr}
        check(f"{label}.size6_fieldfree_{pr}_set_equality", measured == constructed, {"measured": len(measured), "constructed": len(constructed)})
        check(f"{label}.size6_fieldfree_{pr}_formula", len(constructed) == formulas6[pr], {"got": len(constructed), "want": formulas6[pr]})
    fieldfree7_34 = construct_fieldfree_size7_34(n, cell_index)
    measured7_34 = {mask for mask in data["circuits"][7] if profile(mask_support(mask, data["cells"])) == (3, 4)}
    formula7_34 = 12 * math.comb(n, 3) * math.comb(n, 4)
    check(f"{label}.size7_fieldfree_34_set_equality", measured7_34 == fieldfree7_34, {"measured": len(measured7_34), "constructed": len(fieldfree7_34)})
    check(f"{label}.size7_fieldfree_34_formula", len(fieldfree7_34) == formula7_34, {"got": len(fieldfree7_34), "want": formula7_34})
    for m in (5, 6, 7):
        check(f"{label}.m{m}.degree_bounds", data["degree_violations"][m] == 0, data["degree_violations"][m])
        check(f"{label}.m{m}.projection_floor", all(
            profile(mask_support(mask, data["cells"]))[0] >= 3 and profile(mask_support(mask, data["cells"]))[1] >= 4
            for mask in data["circuits"][m]
        ), len(data["circuits"][m]))
    return {
        "n": n,
        "prime": p,
        "ambient_dim": 6,
        "sizes": {
            str(m): {
                "universe": math.comb(n * n, m),
                "circuit_count": len(data["circuits"][m]),
                "profiles": {str(k): v for k, v in sorted(data["profiles"][m].items())},
                "template_signature_count": len(data["signatures"][m]),
                "template_signatures": dict(sorted(data["signatures"][m].items())),
                "degree_violations": data["degree_violations"][m],
            }
            for m in range(3, 8)
        },
        "size5_crossing": {"count": len(crossing5), "formula": want_crossing, "set_equality": True},
        "fieldfree_size6": {
            str(pr): {"count": len(fieldfree6[pr]), "formula": formulas6[pr], "set_equality": True}
            for pr in sorted(fieldfree6)
        },
        "fieldfree_size7_34": {"count": len(fieldfree7_34), "formula": formula7_34, "set_equality": True},
        "rank_cache": data["rank_mask"].cache_info()._asdict(),
    }


def compare_primes(left, right):
    n = left["n"]
    out = {"n": n, "primes": [left["p"], right["p"]], "sizes": {}}
    for m in range(3, 8):
        A, B = left["circuits"][m], right["circuits"][m]
        by_profile = {}
        profiles = set(left["profiles"][m]) | set(right["profiles"][m])
        for pr in sorted(profiles):
            Ap = {mask for mask in A if profile(mask_support(mask, left["cells"])) == pr}
            Bp = {mask for mask in B if profile(mask_support(mask, right["cells"])) == pr}
            by_profile[str(pr)] = {
                "p_left": len(Ap),
                "p_right": len(Bp),
                "intersection": len(Ap & Bp),
                "left_only": len(Ap - Bp),
                "right_only": len(Bp - Ap),
                "sets_equal": Ap == Bp,
            }
        out["sizes"][str(m)] = {
            "left": len(A),
            "right": len(B),
            "intersection": len(A & B),
            "left_only": len(A - B),
            "right_only": len(B - A),
            "sets_equal": A == B,
            "by_profile": by_profile,
        }
    return out


def write_failure(exc):
    payload = {
        "run_id": RUN_ID,
        "gate": RESULTS["gate"],
        "error_type": type(exc).__name__,
        "error": str(exc),
        "traceback": traceback.format_exc(),
        "partial_results": RESULTS,
        "controls": CONTROLS,
        "assertion_count_before_failure": ASSERTIONS,
        "budget": {"CPU_seconds": time.process_time() - START_CPU, "wall_seconds": time.monotonic() - START_WALL},
    }
    (HERE / "broken_results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run():
    print(f"== H-2ROW-COMPLETE-3ROW-OPEN run {RUN_ID}", flush=True)
    check("runtime.dont_write_bytecode", sys.dont_write_bytecode is True, sys.dont_write_bytecode)
    check("runtime.env.PYTHONDONTWRITEBYTECODE", os.environ.get("PYTHONDONTWRITEBYTECODE") == "1", os.environ.get("PYTHONDONTWRITEBYTECODE"))
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        check(f"runtime.env.{var}", os.environ.get(var) == "1", os.environ.get(var))
    nice_value = os.getpriority(os.PRIO_PROCESS, 0)
    check("runtime.nice", nice_value == 15, nice_value)

    two_row = []
    for n in (4, 5):
        for p in (7, 13):
            result = two_row_size6_check(n, p)
            check(f"two_row.n{n}.p{p}.universe", result["universe"] == math.comb(n * n, 6), result["universe"])
            check(f"two_row.n{n}.p{p}.zero_size6_circuits", result["circuits"] == 0, result["circuits"])
            check(f"two_row.n{n}.p{p}.every6_has_dependent5", result["all_have_dependent_5subset"], result)
            two_row.append(result)
            print(f"[P1] two-row n={n}, p={p}: {result['universe']} six-sets, zero circuits", flush=True)

    unequal_data = {}
    unequal_summaries = []
    for n in (4, 5):
        for p in (7, 13):
            data = unequal_census(n, p)
            unequal_data[(n, p)] = data
            summary = summarize_unequal(data)
            unequal_summaries.append(summary)
            counts = {m: len(data["circuits"][m]) for m in range(3, 8)}
            print(f"[P2] unequal n={n}, p={p}: counts {counts}", flush=True)

    prime_comparisons = [compare_primes(unequal_data[(n, 7)], unequal_data[(n, 13)]) for n in (4, 5)]

    # Both-direction controls.
    # Two-row R+1 accept: C4 rectangle plus opposite isolated edge.
    points4, cells22, cols22 = two_row_columns(4, 13)
    index22 = {cell: i for i, cell in enumerate(cells22)}
    rank22 = make_rank_oracle(cols22, 13)
    accept5_support = tuple(sorted(((0, 0), (0, 1), (1, 0), (1, 1), (2, 2))))
    accept5_mask = support_mask(accept5_support, index22)
    check("plants.two_row_size5_accept", is_circuit_mask(accept5_mask, 5, rank22), accept5_support)
    six_support = tuple(cells22[:6])
    six_mask = support_mask(six_support, index22)
    check("plants.two_row_six_reject", not is_circuit_mask(six_mask, 6, rank22), {"rank": rank22(six_mask), "support": six_support})
    check("plants.two_row_six_has_dependent5", any(rank22(six_mask ^ (1 << i)) < 5 for i in mask_indices(six_mask)), six_support)

    data5 = unequal_data[(5, 13)]
    crossing5, _ = construct_unequal_crossings(5, data5["cell_index"])
    known_cross = next(iter(crossing5))
    check("plants.unequal_crossing_accept", known_cross in data5["circuits"][5], mask_support(known_cross, data5["cells"]))

    fiber3 = next(iter(construct_size3_fibers(5, data5["cell_index"])))
    fiber_cells = set(mask_support(fiber3, data5["cells"]))
    extras = [cell for cell in data5["cells"] if cell not in fiber_cells][:2]
    fiber_plant_support = tuple(sorted(fiber_cells | set(extras)))
    fiber_plant_mask = support_mask(fiber_plant_support, data5["cell_index"])
    check("plants.proper_fiber_reject", not is_circuit_mask(fiber_plant_mask, 5, data5["rank_mask"]), fiber_plant_support)

    seven_support = tuple(data5["cells"][:7])
    seven_mask = support_mask(seven_support, data5["cell_index"])
    seven_rank = data5["rank_mask"](seven_mask)
    seven_is_circuit = is_circuit_mask(seven_mask, 7, data5["rank_mask"])
    check("plants.ambient7_automatic_dependence", seven_rank <= 6, {"rank": seven_rank, "support": seven_support})
    check("plants.ambient7_minimality_separate", not seven_is_circuit, {"rank": seven_rank, "is_circuit": seven_is_circuit})
    check("plants.ambient7_positive_population", len(data5["circuits"][7]) > 0, len(data5["circuits"][7]))
    seven_accept = next(iter(data5["circuits"][7]))
    check("plants.ambient7_accept_all_deletions", all(data5["rank_mask"](seven_accept ^ (1 << i)) == 6 for i in mask_indices(seven_accept)), mask_support(seven_accept, data5["cells"]))

    known_cross_support = mask_support(known_cross, data5["cells"])
    known_cols = [data5["cols"][data5["cell_index"][cell]] for cell in known_cross_support]
    corrupt_cols = [list(c) for c in known_cols]
    corrupt_cols[0] = list(corrupt_cols[1])
    corrupt_circuit = rank_cols(corrupt_cols, 13) == 4 and all(rank_cols(corrupt_cols[:i] + corrupt_cols[i + 1 :], 13) == 4 for i in range(5))
    check("plants.non_tensor_label_accept", known_cross in data5["circuits"][5], known_cross_support)
    check("plants.non_tensor_actual_reject", not corrupt_circuit, rank_cols(corrupt_cols, 13))

    factorA = [[1, x] for x in range(1, 6)]
    dupA = factorA + [list(factorA[0])]
    check("plants.duplicate.spark2", measure_spark(dupA, 13) == 2, measure_spark(dupA, 13))
    dep_pairs = [pair for pair in itertools.combinations(range(6), 2) if rank_cols([dupA[i] for i in pair], 13) < 2]
    check("plants.duplicate.pair", dep_pairs == [(0, 5)], dep_pairs)

    a, b = [1, 2], [3, 4, 5]
    true_kron = kron(a, b, 13)
    concat = a + b
    check("plants.kron_guard.true_dim", len(true_kron) == 6, true_kron)
    check("plants.kron_guard.concat_reject", len(concat) == 5 and concat != true_kron, {"true": true_kron, "concat": concat})
    check("plants.kron_guard.coordinates", true_kron == [3, 4, 5, 6, 8, 10], true_kron)

    RESULTS["two_row_size6"] = two_row
    RESULTS["unequal_runs"] = unequal_summaries
    RESULTS["prime_comparisons"] = prime_comparisons
    RESULTS["plants"] = {
        "two_row_size5_accept": accept5_support,
        "two_row_six_reject": six_support,
        "unequal_crossing_accept": known_cross_support,
        "proper_fiber_reject": fiber_plant_support,
        "ambient7": {"dependent_plant_rank": seven_rank, "dependent_plant_is_circuit": seven_is_circuit, "positive_population": len(data5["circuits"][7]), "accepted_example": mask_support(seven_accept, data5["cells"])},
        "non_tensor": {"labels_accept": True, "actual_rejected": True},
        "duplicate": {"spark": 2, "dependent_pairs": dep_pairs},
        "kronecker_guard": {"true_dim": 6, "concat_dim": 5, "true_coordinates": true_kron, "concat": concat},
    }

    cpu = time.process_time() - START_CPU
    wall = time.monotonic() - START_WALL
    check("budget.CPU", cpu <= CPU_CAP, {"used": cpu, "cap": CPU_CAP})
    check("budget.wall", wall <= WALL_CAP, {"used": wall, "cap": WALL_CAP})
    RESULTS["controls"] = CONTROLS
    RESULTS["assertion_count"] = ASSERTIONS
    RESULTS["budget"] = {"CPU_cap_seconds": CPU_CAP, "CPU_seconds": cpu, "wall_cap_seconds": WALL_CAP, "wall_seconds": wall, "within_cap": True}
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
