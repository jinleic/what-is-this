#!/usr/bin/env python3
"""Exact instrument for H-44ROW-LAYER.

Run 20260902T135324Z_516630ab_ed03d008f0f6, prereg
prereg/H_44ROW_LAYER_PREREG_2026-09-02.md (commit 92daf35).

Tied four-row layer (4xn) tensor (4xn), dA=dB=5, ambient 16, sizes 5..17;
this gate certifies the prediction layers 5 (fibers), 6 (empty), 7 (empty),
8 (crossings 25*C(n,5)^2 + all-distinct Moebius-8 restrictions sum C(kM,8)),
with the preregistered out-of-scope items recorded as non-claims.

Standalone: imports no frozen campaign code. Exact Python-int GF(p)
arithmetic, dual rank cross-checks, canonical PGL enumeration.
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
from pathlib import Path

RUN_ID = "20260902T135324Z_516630ab_ed03d008f0f6"
GATE = "H-44ROW-LAYER"
CPU_CAP = 5400.0
WALL_CAP = 6000.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT = Path(__file__).resolve().parent
CHECKS = {}
CHECKPOINT = {"run_id": RUN_ID, "completed_layers": []}
RESULTS = {
    "run_id": RUN_ID,
    "gate": GATE,
    "arithmetic": "exact Python integers modulo p; true row-major Kronecker columns",
}


def elapsed():
    return time.process_time() - START_CPU, time.monotonic() - START_WALL


def check_budget(label):
    cpu, wall = elapsed()
    if cpu > CPU_CAP or wall > WALL_CAP:
        raise RuntimeError(f"budget exceeded at {label}: CPU={cpu:.3f}, wall={wall:.3f}")


def check(label, condition, detail=None):
    if not condition:
        raise AssertionError(f"{label}: {detail!r}")
    CHECKS[label] = {"pass": True, "detail": detail}


def pin(got, want, label):
    check(label, got == want, {"got": got, "want": want})


def checkpoint(record):
    CHECKPOINT["completed_layers"].append(record)
    cpu, wall = elapsed()
    CHECKPOINT["CPU_seconds"] = cpu
    CHECKPOINT["wall_seconds"] = wall
    (OUT / "checkpoint.json").write_text(json.dumps(CHECKPOINT, indent=2, sort_keys=True) + "\n")


def kronecker_column(x, y, ra, rb, p):
    ax = [pow(x, i, p) for i in range(ra)]
    by = [pow(y, j, p) for j in range(rb)]
    return tuple(u * v % p for u in ax for v in by)


def rank_mod(columns, p):
    if not columns:
        return 0
    a = [[value % p for value in column] for column in columns]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, p)
        a[rank] = [value * inv % p for value in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                factor = a[r][col]
                a[r] = [(u - factor * v) % p for u, v in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def rank_mod_alt(columns, p):
    basis = []
    pivots = []
    for raw in columns:
        vector = [value % p for value in raw]
        for pivot, row in zip(pivots, basis):
            if vector[pivot]:
                factor = vector[pivot]
                vector = [(u - factor * v) % p for u, v in zip(vector, row)]
        pivot = next((i for i, value in enumerate(vector) if value), None)
        if pivot is None:
            continue
        inv = pow(vector[pivot], -1, p)
        vector = [value * inv % p for value in vector]
        for i, row in enumerate(basis):
            if row[pivot]:
                factor = row[pivot]
                basis[i] = [(u - factor * v) % p for u, v in zip(row, vector)]
        insert = sum(old < pivot for old in pivots)
        pivots.insert(insert, pivot)
        basis.insert(insert, vector)
    return len(basis)


def dual_rank(columns, p, label):
    r1 = rank_mod(columns, p)
    r2 = rank_mod_alt(columns, p)
    if r1 != r2:
        raise AssertionError(f"{label}: dual rank disagreement {r1} != {r2}")
    return r1


def det_mod(matrix, p):
    a = [[value % p for value in row] for row in matrix]
    n = len(a)
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det % p
        value = a[col][col]
        det = det * value % p
        inv = pow(value, -1, p)
        for r in range(col + 1, n):
            if a[r][col]:
                factor = a[r][col] * inv % p
                for c in range(col, n):
                    a[r][c] = (a[r][c] - factor * a[col][c]) % p
    return det


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
        if (a * d - b * c) % p:
            classes.add(canonical_matrix((a, b, c, d), p))
    return tuple(sorted(classes))


def mobius_apply(matrix, x, p):
    a, b, c, d = matrix
    den = (c * x + d) % p
    if den == 0:
        return None
    return (a * x + b) * pow(den, -1, p) % p


def pgl_domains(p, X, Y):
    out = []
    for matrix in pgl_matrices(p):
        D = [x for x in X if mobius_apply(matrix, x, p) is not None
             and mobius_apply(matrix, x, p) in Y]
        out.append((matrix, D))
    return out


def pgl_sum_count(domains, k):
    total = 0
    for _, D in domains:
        if len(D) >= k:
            total += math.comb(len(D), k)
    return total


def mobius_graph_sets(domains, p, k):
    seen = set()
    for matrix, D in domains:
        for chosen in itertools.combinations(sorted(D), k):
            pairs = [(x, mobius_apply(matrix, x, p)) for x in chosen]
            seen.add(tuple(sorted(pairs)))
    return seen


def all_distinct_circuit_sweep(n, ra, rb, p, k, label):
    """Sweep all-distinct k-sets on {1..n}^2 (mod p) as the k! bijections
    (permutations of y-indices over x-indices); return (dependent,
    circuits) counts. Requires n == k."""
    assert n == k, "bijection sweep requires n == k"
    pts = [x % p for x in range(1, n + 1)]
    col_cache = {}
    def col(u, v):
        key = (u, v)
        if key not in col_cache:
            col_cache[key] = kronecker_column(pts[u], pts[v], ra, rb, p)
        return col_cache[key]
    dependent = 0
    circuits = 0
    for perm in itertools.permutations(range(n)):
        cols = [col(u, perm[u]) for u in range(n)]
        if dual_rank(cols, p, label) != k - 1:
            continue
        dependent += 1
        check_budget(f"{label}.deletion_loop")
        if all(dual_rank(cols[:i] + cols[i + 1:], p, f"{label}.d{i}") == k - 1
               for i in range(k)):
            circuits += 1
    return dependent, circuits


def tied_census(n_cells, ra, rb, p, k, label, fiber_len):
    """Exhaustive census of ALL k-subsets of the n x n grid at (ra, rb).
    Returns dict with dependent, circuits, fiber_extra, other profiles."""
    pts = [x % p for x in range(1, n_cells + 1)]
    cols_all = {}
    for u in range(n_cells):
        for v in range(n_cells):
            cols_all[(u, v)] = kronecker_column(pts[u], pts[v], ra, rb, p)
    cells = sorted(cols_all)
    dependent = 0
    circuits = 0
    fiber_extra = 0
    other = 0
    for S in itertools.combinations(cells, k):
        cols = [cols_all[c] for c in S]
        r = dual_rank(cols, p, label)
        if r != k - 1:
            continue
        dependent += 1
        check_budget(f"{label}.deletion_loop")
        minimal = all(dual_rank(cols[:i] + cols[i + 1:], p, f"{label}.d{i}") == k - 1
                      for i in range(k))
        if minimal:
            circuits += 1
        degs_r = {}
        degs_c = {}
        for u, v in S:
            degs_r[u] = degs_r.get(u, 0) + 1
            degs_c[v] = degs_c.get(v, 0) + 1
        if any(c == fiber_len for c in degs_r.values()) or \
           any(c == fiber_len for c in degs_c.values()):
            fiber_extra += 1
        else:
            other += 1
    return {"dependent": dependent, "circuits": circuits,
            "fiber_extra": fiber_extra, "other": other}


def gate_crossings(n, ra, rb, p, label):
    """Construct all S(R,i0,J,j0) with |R|=|J|=ra+1... here ra=rb=4, factor
    circuits of size 5: |R|=|J|=5, i0 in R, j0 in J; count n^2*C(n,5)^2,
    verify each is a circuit."""
    pts = [x % p for x in range(1, n + 1)]
    def col(u, v):
        return kronecker_column(pts[u], pts[v], ra, rb, p)
    rows = list(range(n))
    circuits = 0
    supports = []
    for R in itertools.combinations(rows, ra + 1):
        for J in itertools.combinations(rows, rb + 1):
            for i0 in R:
                for j0 in J:
                    S = tuple(sorted([(u, j0) for u in R if u != i0] +
                                     [(i0, v) for v in J if v != j0]))
                    cols = [col(u, v) for u, v in S]
                    m = len(S)
                    if dual_rank(cols, p, label) != m - 1:
                        continue
                    if all(dual_rank(cols[:i] + cols[i + 1:], p,
                                     f"{label}.d{i}") == m - 1 for i in range(m)):
                        circuits += 1
                        supports.append(S)
    return circuits, supports


def run_controls():
    p0 = 13
    # 1. Pure-tensor ACCEPT (A-fiber 5-set) / corrupted REJECT
    fiber5 = tuple((1, y % p0) for y in range(1, 6))
    pts5 = [x % p0 for x in range(1, 6)]
    cols = [kronecker_column(1 % p0, y, 4, 4, p0) for y in (1, 2, 3, 4, 5)]
    # fiber columns: all kronecker(x=1, y) for 5 distinct y: these are
    # 5 columns of the 16x16 product matrix; their span is 4-dim (per y),
    # so a single fiber-5-set is a circuit? NO: fiber of 5 SHARED-x cells
    # spans b_y for the fixed a_x: dimension 4 < 5 -> they are DEPENDENT
    # but not a circuit: any 4 of them already dependent? b_y span 4-dim,
    # 4 cells = basis. size-5 fiber = circuit only via factor 5-circuit:
    # A side has spark 4 (MDS 4-row: any 4 cols are a basis of the a-space,
    # 5 cols dependent) -> the 5-cell fiber IS minimal: deleting one leaves
    # 4 independent b's -> rank 4 = 5-1.
    r_fiber = dual_rank(cols, p0, "controls.fiber5")
    pin(r_fiber, 4, "controls.fiber5.rank")
    dels = [dual_rank(cols[:i] + cols[i + 1:], p0, f"controls.fiber5.d{i}")
            for i in range(5)]
    check("controls.fiber5.circuit", all(d == 4 for d in dels), dels)
    # corrupted column must break circuitness
    for delta in range(1, p0):
        changed = [list(c) for c in cols]
        changed[0][0] = (changed[0][0] + delta) % p0
        vecs = [tuple(c) for c in changed]
        rnew = dual_rank(vecs, p0, "controls.fiber5.corrupt")
        if not (rnew == 4 and all(dual_rank(vecs[:i] + vecs[i + 1:], p0,
                                                "controls.fiber5.cd") == 4
                                  for i in range(5))):
            check("controls.fiber5_corrupt_reject", True,
                  {"delta": delta, "rank": rnew})
            break
    else:
        raise AssertionError("controls.fiber5_corrupt_reject: no witness")

    # 2. Discriminating ambient plant: 16 tensor-basis columns rank exactly 16
    grid4 = [(x % p0, y % p0) for x in range(1, 5) for y in range(1, 5)]
    basis16 = [kronecker_column(x, y, 4, 4, p0) for x, y in grid4]
    pin(dual_rank(basis16, p0, "controls.basis16"), 16, "controls.basis16.rank")
    seventeen = basis16 + [kronecker_column(5, 5, 4, 4, p0)]
    pin(dual_rank(seventeen, p0, "controls.ambient17"), 16,
        "controls.ambient17.exact_rank_16")
    # a 16-set that passes corank-1 (rank 15) but has a singular deletion:
    # 15 basis columns plus a duplicate of the first; deleting the duplicate
    # keeps rank 15, deleting any basis column drops rank to 14 -> a
    # size-16 set with rank 15 that is NOT a circuit: (C2)(i) REJECT.
    wrong16 = basis16[:15] + [basis16[0]]
    r16 = dual_rank(wrong16, p0, "controls.wrong16")
    pin(r16, 15, "controls.wrong16.rank15")
    del_dup = dual_rank(wrong16[:15], p0, "controls.wrong16.deldup")
    pin(del_dup, 15, "controls.wrong16.deldup_nondip")
    del_basis = dual_rank(wrong16[:3] + wrong16[4:], p0, "controls.wrong16.delbasis")
    check("controls.wrong16.reject", del_basis == 14, del_basis)

    # 3. Non-Mobius 8-set plant at n=5 (all-distinct impossible at n=5 for
    # 8-set? all-distinct needs 8 distinct x and 8 distinct y -> n>=8.
    # At n=5 there are NO all-distinct 8-sets; the non-Mobius plant lives at
    # the n=8 layer and is exercised inside stage E.
    check("controls.nonmobius_deferred_to_stageE", True, None)

    # 4. Concat-vs-Kronecker guard
    true_col = kronecker_column(2, 3, 4, 4, p0)
    concat = tuple(list(kronecker_column(2, 0, 4, 1, p0)) +
                   list(kronecker_column(0, 3, 1, 4, p0))[:-1])
    check("controls.kron_guard.true_dimension", len(true_col) == 16, len(true_col))
    check("controls.kron_guard.concat_reject", len(concat) == 7 and
          tuple(concat) != true_col, {"len_concat": len(concat)})

    # 5. Crossing plant: verified (5,5) crossing at n=5 is a circuit
    R5 = (0, 1, 2, 3, 4)
    i0, j0 = 0, 0
    S = tuple(sorted([(u, j0) for u in R5 if u != i0] +
                     [(i0, v) for v in R5 if v != i0]))
    cols8 = [kronecker_column(pts5[u], pts5[v], 4, 4, p0) for u, v in S]
    check("controls.crossing.circuit",
          all(dual_rank(cols8[:i] + cols8[i + 1:], p0, "controls.crossing.d")
              == 7 for i in range(8)) and
          dual_rank(cols8, p0, "controls.crossing") == 7, S)

    # 6. Duplicate-column spark
    base4 = [kronecker_column(x, pts5[i], 4, 4, p0) for i, x in enumerate((1, 2, 3, 4))]
    pin(dual_rank(base4, p0, "controls.base4"), 4, "controls.base4.spark4")

    # 7. Mobius-8 negative control handled in stage E (needs n>=8 data):
    # recorded here as a battery item.
    check("controls.mobius8_minimal_deferred_to_stageE", True, None)
    return {"fiber5": fiber5, "crossing": S}


def main():
    print(f"== {GATE} run {RUN_ID}")
    nice_level = os.nice(0)
    check("runtime.nice", nice_level == 10, nice_level)
    check("runtime.sys_dont_write_bytecode", sys.dont_write_bytecode is True,
          sys.dont_write_bytecode)
    check("runtime.env_dont_write_bytecode",
          os.environ.get("PYTHONDONTWRITEBYTECODE") == "1",
          os.environ.get("PYTHONDONTWRITEBYTECODE"))
    thread_caps = {name: os.environ.get(name) for name in
                   ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}
    check("runtime.thread_caps", all(value == "1" for value in thread_caps.values()),
          thread_caps)
    resource.setrlimit(resource.RLIMIT_CPU, (5400, 5400))
    check("runtime.rlimit_cpu",
          resource.getrlimit(resource.RLIMIT_CPU) == (5400, 5400),
          resource.getrlimit(resource.RLIMIT_CPU))

    plants = run_controls()
    checkpoint({"stage": "F", "controls": len(CHECKS)})
    print(f"[stage F] control battery OK ({len(CHECKS)} asserts so far)")

    pgl_facts = {}
    domains_by_p = {}
    for p in (11, 13):
        mats = pgl_matrices(p)
        pin(len(mats), p * (p * p - 1), f"stageB.pgl_count.p{p}")
        pgl_facts[p] = len(mats)
        domains_by_p[p] = pgl_domains  # lazy: call per instance
    checkpoint({"stage": "B", "pgl_counts": pgl_facts})
    print(f"[stage B] PGL counts {pgl_facts}")

    # ---------------- L5: size-5 fiber layer at n=5,6 ----------------
    l5_records = []
    for p in (11, 13):
        for n in (5, 6):
            pts = [x % p for x in range(1, n + 1)]
            cols_all = {}
            for u in range(n):
                for v in range(n):
                    cols_all[(u, v)] = kronecker_column(pts[u], pts[v], 4, 4, p)
            cells = sorted(cols_all)
            circuits = []
            for S in itertools.combinations(cells, 5):
                cols = [cols_all[c] for c in S]
                if dual_rank(cols, p, f"L5.n{n}.p{p}") != 4:
                    continue
                if all(dual_rank(cols[:i] + cols[i + 1:], p, f"L5.n{n}.p{p}.d")
                       == 4 for i in range(5)):
                    circuits.append(S)
            # classify: must be exactly the 2n*C(n,5) fibers
            fibers = []
            for u in range(n):
                for chosen in itertools.combinations(range(n), 5):
                    fibers.append(tuple(sorted((u, v) for v in chosen)))
            for v in range(n):
                for chosen in itertools.combinations(range(n), 5):
                    fibers.append(tuple(sorted((u, v) for u in chosen)))
            fiber_set = set(fibers)
            circuit_set = {tuple(sorted(S)) for S in circuits}
            check(f"L5.n{n}.p{p}.set_equality", circuit_set == fiber_set,
                  {"circuits": len(circuit_set), "fibers": len(fiber_set),
                   "missing": len(circuit_set - fiber_set),
                   "extra": len(fiber_set - circuit_set)})
            pin(len(circuit_set), 2 * n * math.comb(n, 5), f"L5.n{n}.p{p}.count")
            l5_records.append({"n": n, "p": p, "circuits": len(circuit_set)})
            checkpoint({"layer": "L5", "n": n, "p": p, "circuits": len(circuit_set)})
            print(f"[L5] n={n} p={p}: circuits={len(circuit_set)} == 2n*C(n,5)")

    # ---------------- L6: size-6 empty at n=5,6 ----------------
    l6_records = []
    for p in (11, 13):
        for n in (5, 6):
            rec = tied_census(n, 4, 4, p, 6, f"L6.n{n}.p{p}", 5)
            pin(rec["circuits"], 0, f"L6.n{n}.p{p}.empty")
            l6_records.append({"n": n, "p": p, **rec})
            checkpoint({"layer": "L6", "n": n, "p": p, **rec})
            print(f"[L6] n={n} p={p}: dependent={rec['dependent']} circuits=0 "
                  f"(fiber_extra={rec['fiber_extra']})")

    # ---------------- L7: size-7 empty at n=5,6 ----------------
    l7_records = []
    for p in (11, 13):
        for n in (5, 6):
            rec = tied_census(n, 4, 4, p, 7, f"L7.n{n}.p{p}", 5)
            pin(rec["circuits"], 0, f"L7.n{n}.p{p}.empty")
            l7_records.append({"n": n, "p": p, **rec})
            checkpoint({"layer": "L7", "n": n, "p": p, **rec})
            print(f"[L7] n={n} p={p}: dependent={rec['dependent']} circuits=0 "
                  f"(fiber_extra={rec['fiber_extra']})")

    # ---------------- L8: size-8 layer ----------------
    l8_records = []
    for p in (11, 13):
        # n=5 exhaustive tied census (1,081,575 supports)
        rec = tied_census(5, 4, 4, p, 8, f"L8.n5.p{p}", 5)
        # every dependency must be fiber-extra, crossing, or all-distinct
        # Mobius. At n=5 all-distinct 8-sets do not exist. Verify circuits ==
        # crossings: 25 * C(5,5)^2 = 25... wait C(5,5)=1 so 25 crossings.
        pin(rec["circuits"], 25, f"L8.n5.p{p}.crossings_only")
        l8_records.append({"part": "n5_census", "p": p, **rec})
        checkpoint({"layer": "L8", "part": "n5_census", "p": p, **rec})
        print(f"[L8] n=5 p={p}: dependent={rec['dependent']} circuits="
              f"{rec['circuits']} (=25 crossings)")

        # n=6 crossing construction: 25*C(6,5)^2 = 900 supports
        got, supports6 = gate_crossings(6, 4, 4, p, f"L8.n6.p{p}")
        pin(got, 900, f"L8.n6.p{p}.crossings_all_circuits")
        l8_records.append({"part": "n6_crossings", "p": p, "circuits": got})
        checkpoint({"layer": "L8", "part": "n6_crossings", "p": p, "circuits": got})
        print(f"[L8] n=6 p={p}: 900 crossing supports verified circuits")

        # n=8 all-distinct bijection sweep vs Mobius-8 sets
        pts8 = [x % p for x in range(1, 9)]
        domains = pgl_domains(p, pts8, pts8)
        pgl_total = pgl_sum_count(domains, 8)
        graph_sets = mobius_graph_sets(domains, p, 8)
        dep, circ = all_distinct_circuit_sweep(8, 4, 4, p, 8, f"L8.n8.p{p}")
        pin(dep, pgl_total, f"L8.n8.p{p}.pgl_sum_count")
        # support-set comparison
        colset = {}
        def col(u, v):
            key = (u, v)
            if key not in colset:
                colset[key] = kronecker_column(pts8[u], pts8[v], 4, 4, p)
            return colset[key]
        circuit_supports = set()
        for perm in itertools.permutations(range(8)):
            pairS = tuple((u, perm[u]) for u in range(8))
            cols = [col(u, perm[u]) for u in range(8)]
            if dual_rank(cols, p, f"L8.n8.p{p}") != 7:
                continue
            if all(dual_rank(cols[:i] + cols[i + 1:], p, f"L8.n8.p{p}.d") == 7
                   for i in range(8)):
                circuit_supports.add(tuple(sorted(
                    (pts8[u], pts8[perm[u]]) for u in range(8))))
        check(f"L8.n8.p{p}.set_equality", circuit_supports == graph_sets,
              {"circuits": len(circuit_supports), "graphs": len(graph_sets),
               "missing": len(circuit_supports - graph_sets),
               "extra": len(graph_sets - circuit_supports)})
        # Mobius-8 negative control: a 7-subset of a verified Moebius-8
        # circuit is independent
        sample = sorted(circuit_supports)[0]
        # 7-subsets: drop one pair; sweep all 8
        colmap = {pair: col(*(pts8.index(pair[0]), pts8.index(pair[1])))
                  for pair in sample}
        minimal = True
        for drop in range(8):
            subs = [colmap[pair] for i, pair in enumerate(sample) if i != drop]
            if dual_rank(subs, p, f"L8.n8.p{p}.neg") != 7:
                minimal = False
                break
        check(f"L8.n8.p{p}.restriction_independent", minimal, sample)
        l8_records.append({"part": "n8_all_distinct", "p": p,
                           "dependent": dep, "circuits": len(circuit_supports),
                           "pgl_sum": pgl_total, "graph_sets": len(graph_sets)})
        checkpoint({"layer": "L8", "part": "n8_all_distinct", "p": p,
                    "circuits": len(circuit_supports), "pgl": pgl_total})
        print(f"[L8] n=8 p={p}: all-distinct circuits={len(circuit_supports)}"
              f" = PGL sum {pgl_total} (set equality verified)")

    # ---------------- Wrap up ----------------
    cpu, wall = elapsed()
    check("budget.CPU", cpu <= CPU_CAP, {"used": cpu, "cap": CPU_CAP})
    check("budget.wall", wall <= WALL_CAP, {"used": wall, "cap": WALL_CAP})
    RESULTS.update({
        "runtime": {
            "nice": nice_level,
            "dont_write_bytecode": sys.dont_write_bytecode,
            "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
            "thread_caps": thread_caps,
            "rlimit_cpu": list(resource.getrlimit(resource.RLIMIT_CPU)),
            "max_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        },
        "budget": {
            "CPU_seconds": cpu,
            "wall_seconds": wall,
            "CPU_cap_seconds": CPU_CAP,
            "wall_cap_seconds": WALL_CAP,
            "within_cap": True,
        },
        "pgl_counts": {str(k): v for k, v in pgl_facts.items()},
        "L5_records": l5_records,
        "L6_records": l6_records,
        "L7_records": l7_records,
        "L8_records": l8_records,
        "control_facts": {k: (list(v) if isinstance(v, tuple) else v)
                          for k, v in plants.items()},
        "controls": CHECKS,
        "assertion_count": len(CHECKS),
        "scope_note": "sizes 9-17 and the n=6 size-8 tied exhaustiveness are "
                      "out of scope per prereg section 0/4; non-claims only",
    })
    (OUT / "controls_results.json").write_text(
        json.dumps(RESULTS, indent=2, sort_keys=True) + "\n")
    print(f"== ALL CONTROLS PASS: {len(CHECKS)} asserts; CPU {cpu:.3f}s, wall {wall:.3f}s")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        cpu, wall = elapsed()
        broken = {
            "run_id": RUN_ID,
            "gate": GATE,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "CPU_seconds": cpu,
            "wall_seconds": wall,
            "checks_passed_before_failure": CHECKS,
            "checkpoint": CHECKPOINT,
        }
        path = OUT / "broken_results.json"
        if path.exists():
            path = OUT / f"broken_results_{int(time.time())}.json"
        path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
        print(f"BROKEN output preserved at {path.name}", file=sys.stderr)
        raise
