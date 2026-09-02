#!/usr/bin/env python3
"""Exact instrument for H-DP2-SIZE5; standalone, no frozen imports."""

import sys

sys.dont_write_bytecode = True

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
RUN_ID = "20260902T021721Z_7629001e_215016bf40c3"
HERE = Path(__file__).resolve().parent
START_CPU = time.process_time()
START_WALL = time.monotonic()
CONTROLS = {}
ASSERTIONS = 0
RESULTS = {
    "run_id": RUN_ID,
    "gate": "H-DP2-SIZE5",
    "arithmetic": "exact Python integers modulo p; product columns (1,y,x,xy)",
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


def det_cols(cols, p):
    if len(cols) != len(cols[0]):
        raise ValueError("determinant requires square column matrix")
    n = len(cols)
    a = [[cols[c][r] % p for c in range(n)] for r in range(n)]
    det = 1
    for col in range(n):
        pivot = col
        while pivot < n and a[pivot][col] == 0:
            pivot += 1
        if pivot == n:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = (-det) % p
        pv = a[col][col]
        det = (det * pv) % p
        inv = pow(pv, -1, p)
        for rr in range(col + 1, n):
            if a[rr][col]:
                scale = (a[rr][col] * inv) % p
                for cc in range(col, n):
                    a[rr][cc] = (a[rr][cc] - scale * a[col][cc]) % p
    return det


def product_col(cell, points, p):
    i, j = cell
    x, y = points[i] % p, points[j] % p
    return [1, y, x, (x * y) % p]


def support_cols(support, points, p):
    return [product_col(cell, points, p) for cell in support]


def profile(support):
    return (len({i for i, _ in support}), len({j for _, j in support}))


def is_circuit5(support, points, p):
    cols = support_cols(support, points, p)
    if rank_cols(cols, p) != 4:
        return False
    return all(rank_cols(cols[:i] + cols[i + 1 :], p) == 4 for i in range(5))


def det_support4(support, points, p):
    return det_cols(support_cols(tuple(support), points, p), p)


def degrees(support):
    dr = Counter(i for i, _ in support)
    dc = Counter(j for _, j in support)
    return dr, dc


def component_data(support):
    adjacency = defaultdict(set)
    edge_set = set(support)
    for i, j in support:
        r, c = ("r", i), ("c", j)
        adjacency[r].add(c)
        adjacency[c].add(r)
    seen = set()
    comps = []
    for start in list(adjacency):
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
        rows = {v for side, v in nodes if side == "r"}
        cols = {v for side, v in nodes if side == "c"}
        edges = {(i, j) for i, j in edge_set if i in rows and j in cols}
        comps.append({
            "rows": rows,
            "cols": cols,
            "edges": edges,
            "n_edges": len(edges),
            "n_vertices": len(nodes),
            "cycle": len(edges) == len(nodes),
        })
    return comps


def is_crossing4(support):
    if len(support) != 4 or profile(support) != (3, 3):
        return False
    comps = component_data(support)
    return sorted(c["n_edges"] for c in comps) == [2, 2]


def structural_accept5(support, points, p):
    dr, dc = degrees(support)
    if max(dr.values()) > 2 or max(dc.values()) > 2:
        return False
    for sub in itertools.combinations(support, 4):
        pr = profile(sub)
        if pr == (3, 3) and is_crossing4(sub):
            return False
        if pr == (4, 4) and det_support4(sub, points, p) == 0:
            return False
    return True


def accepted_template(support):
    pr = profile(support)
    comps = component_data(support)
    edge_parts = tuple(sorted((c["n_edges"] for c in comps), reverse=True))
    cycles = tuple(sorted((c["n_edges"] for c in comps if c["cycle"]), reverse=True))
    if pr == (3, 3) and edge_parts == (4, 1) and cycles == (4,):
        return "33:C4+K2"
    if pr == (3, 4) and edge_parts == (4, 1):
        return "34:P5+K2"
    if pr == (4, 3) and edge_parts == (4, 1):
        return "43:P5+K2"
    if pr == (3, 5) and edge_parts == (2, 2, 1):
        return "35:2P3(row)+K2"
    if pr == (5, 3) and edge_parts == (2, 2, 1):
        return "53:2P3(col)+K2"
    if pr == (4, 4) and edge_parts == (3, 1, 1):
        return "44:P4+2K2"
    if pr == (4, 5) and edge_parts == (2, 1, 1, 1):
        return "45:P3(row)+3K2"
    if pr == (5, 4) and edge_parts == (2, 1, 1, 1):
        return "54:P3(col)+3K2"
    if pr == (5, 5) and edge_parts == (1, 1, 1, 1, 1):
        return "55:5K2"
    return "UNEXPECTED"


def full_census(n, p):
    points = list(range(1, n + 1))
    cells = [(i, j) for i in range(n) for j in range(n)]
    direct = set()
    structural = set()
    profiles = Counter()
    templates = Counter()
    for comb in itertools.combinations(cells, 5):
        support = tuple(comb)
        if is_circuit5(support, points, p):
            direct.add(support)
            profiles[profile(support)] += 1
            templates[accepted_template(support)] += 1
        if structural_accept5(support, points, p):
            structural.add(support)
    return {
        "points": points,
        "direct": direct,
        "structural": structural,
        "profiles": profiles,
        "templates": templates,
        "universe": math.comb(n * n, 5),
    }


def triple_relation(indices, points, p):
    a, b, c = tuple(sorted(indices))
    xa, xb, xc = points[a] % p, points[b] % p, points[c] % p
    rel = {
        a: (xb - xc) % p,
        b: (xc - xa) % p,
        c: (xa - xb) % p,
    }
    if any(v == 0 for v in rel.values()):
        raise AssertionError("triple relation has zero coefficient")
    if sum(rel.values()) % p or sum(rel[i] * (points[i] % p) for i in rel) % p:
        raise AssertionError("invalid triple relation")
    return rel


def add_entry(matrix, key, value, p):
    matrix[key] = (matrix.get(key, 0) + value) % p


def relation_is_product_zero(gamma, points, p):
    total = [0, 0, 0, 0]
    for cell, coeff in gamma.items():
        col = product_col(cell, points, p)
        total = [(x + coeff * y) % p for x, y in zip(total, col)]
    return total == [0, 0, 0, 0]


def verify_12_decomposition(R, w, j0, Z1, Z2, support, points, p):
    active_rows = sorted(set(R) - {w})
    u1, u2 = active_rows
    c_rel = triple_relation(R, points, p)
    gamma = {}
    C1 = {}
    C2 = {}
    for u, coeff in c_rel.items():
        C1[(u, j0)] = coeff
        add_entry(gamma, (u, j0), coeff, p)
    for u, Z in ((u1, Z1), (u2, Z2)):
        z_rel = triple_relation(Z, points, p)
        beta = (-c_rel[u] * pow(z_rel[j0], -1, p)) % p
        for v, coeff in z_rel.items():
            value = (beta * coeff) % p
            C2[(u, v)] = value
            add_entry(gamma, (u, v), value, p)
    nonzero_support = tuple(sorted(k for k, v in gamma.items() if v % p))
    return (
        nonzero_support == support
        and gamma[(u1, j0)] == 0
        and gamma[(u2, j0)] == 0
        and relation_is_product_zero({k: v for k, v in gamma.items() if v}, points, p)
        and len(C1) == 3
        and len(C2) == 6
    )


def construct_12(n, p):
    points = list(range(1, n + 1))
    multiplicity = Counter()
    by_profile = defaultdict(set)
    relation_ok = 0
    parameter_count = 0
    triples = list(itertools.combinations(range(n), 3))
    for R in triples:
        for w in R:
            active_rows = sorted(set(R) - {w})
            u1, u2 = active_rows
            for j0 in range(n):
                z_options = [Z for Z in triples if j0 in Z]
                for Z1 in z_options:
                    for Z2 in z_options:
                        inter = len(set(Z1) & set(Z2))
                        if inter not in (1, 2, 3):
                            continue
                        support = tuple(sorted(
                            {(w, j0)}
                            | {(u1, v) for v in Z1 if v != j0}
                            | {(u2, v) for v in Z2 if v != j0}
                        ))
                        if len(support) != 5:
                            raise AssertionError("constructed (1,2) support is not size five")
                        parameter_count += 1
                        multiplicity[support] += 1
                        by_profile[profile(support)].add(support)
                        if verify_12_decomposition(R, w, j0, Z1, Z2, support, points, p):
                            relation_ok += 1
    return {
        "supports": set(multiplicity),
        "multiplicity": multiplicity,
        "by_profile": by_profile,
        "parameter_count": parameter_count,
        "relation_ok": relation_ok,
    }


def transpose_support(support):
    return tuple(sorted((j, i) for i, j in support))


def all_distinct_matchings(n, k):
    for rows in itertools.combinations(range(n), k):
        for cols in itertools.combinations(range(n), k):
            for perm in itertools.permutations(cols):
                yield tuple(sorted(zip(rows, perm)))


def all_distinct4_data(n, points, p):
    dependent = set()
    independent = set()
    for support in all_distinct_matchings(n, 4):
        if det_support4(support, points, p) == 0:
            dependent.add(support)
        else:
            independent.add(support)
    return dependent, independent


def construct_44(n, points, p):
    dependent4, independent4 = all_distinct4_data(n, points, p)
    multiplicity = Counter()
    for Q in independent4:
        for first in Q:
            for second in Q:
                if first == second:
                    continue
                cross = (second[0], first[1])
                if cross in Q:
                    raise AssertionError("ordered distinct matching edges yielded existing cross")
                support = tuple(sorted(set(Q) | {cross}))
                multiplicity[support] += 1
    return {
        "supports": set(multiplicity),
        "multiplicity": multiplicity,
        "dependent4": dependent4,
        "independent4": independent4,
    }


def construct_45(n, points, p):
    out = set()
    candidate_count = 0
    for rows in itertools.combinations(range(n), 4):
        for cols in itertools.combinations(range(n), 5):
            for repeated_row in rows:
                other_rows = tuple(r for r in rows if r != repeated_row)
                for repeated_cols in itertools.combinations(cols, 2):
                    other_cols = tuple(c for c in cols if c not in repeated_cols)
                    for perm in itertools.permutations(other_cols):
                        candidate_count += 1
                        base = set(zip(other_rows, perm))
                        support = tuple(sorted(base | {(repeated_row, repeated_cols[0]), (repeated_row, repeated_cols[1])}))
                        q1 = tuple(sorted(base | {(repeated_row, repeated_cols[0])}))
                        q2 = tuple(sorted(base | {(repeated_row, repeated_cols[1])}))
                        if det_support4(q1, points, p) and det_support4(q2, points, p):
                            out.add(support)
    return out, candidate_count


def construct_55(n, points, p):
    out = set()
    candidates = set(all_distinct_matchings(n, 5)) if n >= 5 else set()
    for support in candidates:
        if all(det_support4(sub, points, p) for sub in itertools.combinations(support, 4)):
            out.add(support)
    return out, len(candidates)


def build_decomposition(C1, C2, p):
    gamma = {}
    for key, val in C1.items():
        add_entry(gamma, key, val, p)
    for key, val in C2.items():
        add_entry(gamma, key, val, p)
    return gamma


def verify_axis_kernels(C1, C2, points, p):
    by_col = defaultdict(dict)
    by_row = defaultdict(dict)
    for (u, v), coeff in C1.items():
        if coeff % p:
            by_col[v][u] = coeff % p
    for (u, v), coeff in C2.items():
        if coeff % p:
            by_row[u][v] = coeff % p
    for rel in by_col.values():
        if sum(rel.values()) % p or sum(rel[i] * points[i] for i in rel) % p:
            return False
    for rel in by_row.values():
        if sum(rel.values()) % p or sum(rel[j] * points[j] for j in rel) % p:
            return False
    return True


def decompose_44(support, points, p):
    dr, dc = degrees(support)
    u1 = next((u for u, d in dr.items() if d == 2), None)
    v1 = next((v for v, d in dc.items() if d == 2), None)
    if u1 is None or v1 is None or (u1, v1) not in support:
        return None
    w1 = next(u for u, v in support if v == v1 and u != u1)
    x1 = next(v for u, v in support if u == u1 and v != v1)
    isolates = [next(iter(comp["edges"])) for comp in component_data(support) if comp["n_edges"] == 1]
    if len(isolates) != 2:
        return None
    for first, second in (isolates, list(reversed(isolates))):
        w2, v2 = first
        u2, x2 = second
        R1 = triple_relation((u1, u2, w1), points, p)
        R2 = triple_relation((u1, u2, w2), points, p)
        Z1 = triple_relation((v1, v2, x1), points, p)
        Z2 = triple_relation((v1, v2, x2), points, p)
        alpha2 = 1
        beta1 = (-alpha2 * R2[u1] * pow(Z1[v2], -1, p)) % p
        beta2 = (-alpha2 * R2[u2] * pow(Z2[v2], -1, p)) % p
        alpha1 = (-beta2 * Z2[v1] * pow(R1[u2], -1, p)) % p
        C1 = {}
        C2 = {}
        for u, val in R1.items():
            C1[(u, v1)] = (alpha1 * val) % p
        for u, val in R2.items():
            C1[(u, v2)] = (alpha2 * val) % p
        for v, val in Z1.items():
            C2[(u1, v)] = (beta1 * val) % p
        for v, val in Z2.items():
            C2[(u2, v)] = (beta2 * val) % p
        gamma = build_decomposition(C1, C2, p)
        nz = tuple(sorted(k for k, v in gamma.items() if v))
        cancellations = sum(1 for k in set(C1) & set(C2) if gamma.get(k, 0) == 0)
        if nz == support and cancellations == 3 and verify_axis_kernels(C1, C2, points, p):
            return {"C1": C1, "C2": C2, "T1": 6, "T2": 6, "r": 4, "c": 3}
    return None


def solve_pair(col1, col2, target, points, p):
    x1, x2 = points[col1] % p, points[col2] % p
    # y1 + y2 = target[0], x1*y1 + x2*y2 = target[1]
    y2 = ((target[1] - x1 * target[0]) * pow((x2 - x1) % p, -1, p)) % p
    y1 = (target[0] - y2) % p
    return y1, y2


def decompose_45(support, points, p):
    dr, _ = degrees(support)
    u1 = next((u for u, d in dr.items() if d == 2), None)
    if u1 is None:
        return None
    repeat_cols = sorted(v for u, v in support if u == u1)
    isolates = [next(iter(comp["edges"])) for comp in component_data(support) if comp["n_edges"] == 1]
    if len(isolates) != 3:
        return None
    for chosen in range(3):
        u2, x3 = isolates[chosen]
        remaining = [isolates[i] for i in range(3) if i != chosen]
        for first, second in (remaining, list(reversed(remaining))):
            w1, v1 = first
            w2, v2 = second
            R1 = triple_relation((u1, u2, w1), points, p)
            R2 = triple_relation((u1, u2, w2), points, p)
            Z2 = triple_relation((v1, v2, x3), points, p)
            beta2 = 1
            alpha1 = (-beta2 * Z2[v1] * pow(R1[u2], -1, p)) % p
            alpha2 = (-beta2 * Z2[v2] * pow(R2[u2], -1, p)) % p
            yv1 = (-alpha1 * R1[u1]) % p
            yv2 = (-alpha2 * R2[u1]) % p
            x1, x2 = repeat_cols
            target = (
                (-(yv1 + yv2)) % p,
                (-(yv1 * points[v1] + yv2 * points[v2])) % p,
            )
            yx1, yx2 = solve_pair(x1, x2, target, points, p)
            C1 = {}
            C2 = {}
            for u, val in R1.items():
                C1[(u, v1)] = (alpha1 * val) % p
            for u, val in R2.items():
                C1[(u, v2)] = (alpha2 * val) % p
            row1 = {v1: yv1, v2: yv2, x1: yx1, x2: yx2}
            for v, val in row1.items():
                C2[(u1, v)] = val % p
            for v, val in Z2.items():
                C2[(u2, v)] = (beta2 * val) % p
            gamma = build_decomposition(C1, C2, p)
            nz = tuple(sorted(k for k, v in gamma.items() if v))
            cancellations = sum(1 for k in set(C1) & set(C2) if gamma.get(k, 0) == 0)
            if nz == support and cancellations == 4 and len([v for v in row1.values() if v]) == 4 and verify_axis_kernels(C1, C2, points, p):
                return {"C1": C1, "C2": C2, "T1": 6, "T2": 7, "r": 4, "c": 4}
    return None


def cofactor_relation(support, points, p):
    cols = support_cols(support, points, p)
    gamma = []
    for i in range(5):
        minor = det_cols(cols[:i] + cols[i + 1 :], p)
        gamma.append(minor if i % 2 == 0 else (-minor) % p)
    total = [sum(gamma[i] * cols[i][r] for i in range(5)) % p for r in range(4)]
    if total != [0, 0, 0, 0]:
        gamma = [(-x) % p for x in gamma]
        total = [sum(gamma[i] * cols[i][r] for i in range(5)) % p for r in range(4)]
    if total != [0, 0, 0, 0]:
        raise AssertionError("cofactor relation failed")
    return gamma


def decompose_55(support, points, p):
    if profile(support) != (5, 5):
        return None
    ordered = list(support)
    gamma_vals = cofactor_relation(ordered, points, p)
    if any(v == 0 for v in gamma_vals):
        return None
    gamma = {cell: gamma_vals[i] for i, cell in enumerate(ordered)}
    (u1, v1), (u2, v2) = ordered[0], ordered[1]
    C2 = {}
    C1 = {}
    for idx, (u, v) in enumerate(ordered):
        l1, l2 = solve_pair(u1, u2, (1, points[u] % p), points, p)
        C2[(u1, v)] = (l1 * gamma_vals[idx]) % p
        C2[(u2, v)] = (l2 * gamma_vals[idx]) % p
    for key in set(gamma) | set(C2):
        value = (gamma.get(key, 0) - C2.get(key, 0)) % p
        if value:
            C1[key] = value
    combined = build_decomposition(C1, C2, p)
    nz = tuple(sorted(k for k, v in combined.items() if v))
    by_col = Counter(v for u, v in C1 if C1[(u, v)])
    by_row = Counter(u for u, v in C2 if C2[(u, v)])
    r = len(set(C1) & set(C2))
    c = sum(1 for key in set(C1) & set(C2) if combined.get(key, 0) == 0)
    if (
        nz == support
        and sorted(by_col.values()) == [3, 3, 3]
        and sorted(by_row.values()) == [4, 4]
        and len(C1) == 9
        and len([v for v in C2.values() if v]) == 8
        and r == c == 6
        and verify_axis_kernels(C1, C2, points, p)
    ):
        return {"C1": C1, "C2": C2, "T1": 9, "T2": 8, "r": 6, "c": 6}
    return None


def measure_spark(cols, p):
    for k in range(1, len(cols) + 1):
        for inds in itertools.combinations(range(len(cols)), k):
            if rank_cols([cols[i] for i in inds], p) < k:
                return k
    return len(cols) + 1


def summarize_run(n, p, census, family12):
    points = census["points"]
    direct = census["direct"]
    structural = census["structural"]
    label = f"n{n}.p{p}"
    check(f"{label}.direct_structural_set_equality", direct == structural, {
        "direct_only": len(direct - structural),
        "structural_only": len(structural - direct),
    })
    check(f"{label}.templates_expected", "UNEXPECTED" not in census["templates"], dict(census["templates"]))

    expected_rows3 = set().union(*(family12["by_profile"].get(pr, set()) for pr in ((3, 3), (3, 4), (3, 5))))
    measured_rows3 = {s for s in direct if profile(s)[0] == 3}
    check(f"{label}.fieldfree_rows3_set_equality", measured_rows3 == expected_rows3, {
        "measured": len(measured_rows3), "constructed": len(expected_rows3),
    })
    mirrored = {transpose_support(s) for s in family12["supports"]}
    measured_cols3 = {s for s in direct if profile(s)[1] == 3}
    check(f"{label}.fieldfree_cols3_set_equality", measured_cols3 == mirrored, {
        "measured": len(measured_cols3), "constructed": len(mirrored),
    })
    check(f"{label}.profile33_orientation_coincides", family12["by_profile"][(3, 3)] == {transpose_support(s) for s in family12["by_profile"][(3, 3)]}, len(family12["by_profile"][(3, 3)]))
    check(f"{label}.family12_parameter_injective", all(v == 1 for v in family12["multiplicity"].values()), max(family12["multiplicity"].values(), default=0))
    check(f"{label}.family12_all_relations", family12["relation_ok"] == family12["parameter_count"], {"ok": family12["relation_ok"], "parameters": family12["parameter_count"]})

    want33 = 9 * math.comb(n, 3) ** 2
    want34 = 72 * math.comb(n, 3) * math.comb(n, 4)
    want35 = 18 * n * math.comb(n, 3) * math.comb(n - 1, 4)
    for pr, want in (((3, 3), want33), ((3, 4), want34), ((4, 3), want34), ((3, 5), want35), ((5, 3), want35)):
        check(f"{label}.formula.{pr}", census["profiles"].get(pr, 0) == want, {"got": census["profiles"].get(pr, 0), "want": want})

    c44 = construct_44(n, points, p)
    measured44 = {s for s in direct if profile(s) == (4, 4)}
    check(f"{label}.P44_set_equality", measured44 == c44["supports"], {"measured": len(measured44), "constructed": len(c44["supports"])})
    check(f"{label}.P44_parameter_injective", all(v == 1 for v in c44["multiplicity"].values()), max(c44["multiplicity"].values(), default=0))
    M4 = 24 * math.comb(n, 4) ** 2
    want44 = 12 * (M4 - len(c44["dependent4"]))
    check(f"{label}.P44_closed_reduction", len(measured44) == want44, {"got": len(measured44), "M4": M4, "N_all4": len(c44["dependent4"]), "want": want44})
    decomp44_ok = sum(decompose_44(s, points, p) is not None for s in measured44)
    check(f"{label}.P44_canonical_decompositions", decomp44_ok == len(measured44), {"ok": decomp44_ok, "total": len(measured44)})

    c45, candidates45 = construct_45(n, points, p)
    measured45 = {s for s in direct if profile(s) == (4, 5)}
    check(f"{label}.P45_set_equality", measured45 == c45, {"measured": len(measured45), "constructed": len(c45), "candidates": candidates45})
    decomp45_ok = sum(decompose_45(s, points, p) is not None for s in measured45)
    check(f"{label}.P45_canonical_decompositions", decomp45_ok == len(measured45), {"ok": decomp45_ok, "total": len(measured45)})

    c54 = {transpose_support(s) for s in c45}
    measured54 = {s for s in direct if profile(s) == (5, 4)}
    check(f"{label}.P54_set_equality", measured54 == c54, {"measured": len(measured54), "constructed": len(c54)})
    decomp54_ok = sum(decompose_45(transpose_support(s), points, p) is not None for s in measured54)
    check(f"{label}.P54_canonical_decompositions", decomp54_ok == len(measured54), {"ok": decomp54_ok, "total": len(measured54)})

    c55, candidates55 = construct_55(n, points, p)
    measured55 = {s for s in direct if profile(s) == (5, 5)}
    check(f"{label}.P55_set_equality", measured55 == c55, {"measured": len(measured55), "constructed": len(c55), "candidates": candidates55})
    decomp55_ok = sum(decompose_55(s, points, p) is not None for s in measured55)
    check(f"{label}.P55_canonical_decompositions", decomp55_ok == len(measured55), {"ok": decomp55_ok, "total": len(measured55)})

    return {
        "n": n,
        "prime": p,
        "universe": census["universe"],
        "total": len(direct),
        "profiles": {str(k): v for k, v in sorted(census["profiles"].items())},
        "templates": dict(sorted(census["templates"].items())),
        "fieldfree": {"N33": want33, "N34=N43": want34, "N35=N53": want35, "constructed_supports_12": len(family12["supports"])},
        "P44": {"count": len(measured44), "M4": M4, "N_all4": len(c44["dependent4"]), "formula": want44, "decompositions": decomp44_ok},
        "P45": {"count": len(measured45), "candidate_graphs": candidates45, "decompositions": decomp45_ok},
        "P54": {"count": len(measured54), "candidate_graphs": candidates45, "decompositions": decomp54_ok},
        "P55": {"count": len(measured55), "candidate_matchings": candidates55, "decompositions": decomp55_ok},
        "full_set_equality": True,
    }


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
    print(f"== H-DP2-SIZE5 run {RUN_ID}", flush=True)
    check("runtime.dont_write_bytecode", sys.dont_write_bytecode is True, sys.dont_write_bytecode)
    check("runtime.env.PYTHONDONTWRITEBYTECODE", os.environ.get("PYTHONDONTWRITEBYTECODE") == "1", os.environ.get("PYTHONDONTWRITEBYTECODE"))
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        check(f"runtime.env.{var}", os.environ.get(var) == "1", os.environ.get(var))
    nice_value = os.getpriority(os.PRIO_PROCESS, 0)
    check("runtime.nice", nice_value == 15, nice_value)

    run_summaries = []
    cached = {}
    for n in (4, 5):
        for p in (7, 11, 13):
            census = full_census(n, p)
            family12 = construct_12(n, p)
            summary = summarize_run(n, p, census, family12)
            cached[(n, p)] = (census, family12, summary)
            run_summaries.append(summary)
            print(f"[T] n={n}, p={p}: {summary['total']} circuits; profiles {summary['profiles']}", flush=True)

    # Supplied anchors.
    expected_n4_p13 = {(3, 3): 144, (3, 4): 288, (4, 3): 288, (4, 4): 144}
    expected_n4_p7 = {(3, 3): 144, (3, 4): 288, (4, 3): 288, (4, 4): 192}
    for p, want_total, want_profiles in ((13, 864, expected_n4_p13), (7, 912, expected_n4_p7)):
        census = cached[(4, p)][0]
        check(f"anchors.n4.p{p}.total", len(census["direct"]) == want_total, {"got": len(census["direct"]), "want": want_total})
        check(f"anchors.n4.p{p}.profiles", dict(census["profiles"]) == want_profiles, {"got": {str(k): v for k, v in census["profiles"].items()}, "want": {str(k): v for k, v in want_profiles.items()}})

    expected_n5_p13 = {
        (3, 3): 900,
        (3, 4): 3600,
        (3, 5): 900,
        (4, 3): 3600,
        (4, 4): 6192,
        (4, 5): 864,
        (5, 3): 900,
        (5, 4): 864,
        (5, 5): 60,
    }
    census_5_13 = cached[(5, 13)][0]
    check("anchors.n5.p13.total", len(census_5_13["direct"]) == 17880, {"got": len(census_5_13["direct"]), "want": 17880})
    check("anchors.n5.p13.profiles", dict(census_5_13["profiles"]) == expected_n5_p13, {"got": {str(k): v for k, v in census_5_13["profiles"].items()}, "want": {str(k): v for k, v in expected_n5_p13.items()}})

    # Plants and both directions.
    points5 = list(range(1, 6))
    census5 = cached[(5, 13)][0]
    family5 = cached[(5, 13)][1]
    accepts = {}
    for pr in ((3, 3), (3, 4), (3, 5)):
        support = next(iter(family5["by_profile"][pr]))
        accepts[str(pr)] = support
        check(f"plants.accept.{pr}.direct", support in census5["direct"], support)
        check(f"plants.accept.{pr}.structural", support in census5["structural"], support)
        check(f"plants.accept.{pr}.all_minors", all(det_support4(sub, points5, 13) for sub in itertools.combinations(support, 4)), support)

    crossing4 = tuple(sorted(((0, 1), (0, 2), (1, 0), (2, 0))))
    crossing_plant = tuple(sorted(set(crossing4) | {(3, 3)}))
    check("plants.inserted_crossing.proper_minor", is_crossing4(crossing4) and det_support4(crossing4, points5, 13) == 0, crossing4)
    check("plants.inserted_crossing.direct_reject", not is_circuit5(crossing_plant, points5, 13), crossing_plant)
    check("plants.inserted_crossing.structural_reject", not structural_accept5(crossing_plant, points5, 13), crossing_plant)

    all55 = set(all_distinct_matchings(5, 5))
    circuits55 = {s for s in census5["direct"] if profile(s) == (5, 5)}
    check("plants.rank_sum_only.all_candidates", len(all55) == 120, len(all55))
    check("plants.rank_sum_only.factor_ranks", all(rank_cols([[1, points5[i]] for i, _ in s], 13) + rank_cols([[1, points5[j]] for _, j in s], 13) == 4 <= 5 for s in all55), True)
    check("plants.rank_sum_only.insufficient", len(circuits55) == 60 and len(all55 - circuits55) == 60, {"accepted_by_rank_sum": 120, "circuits": len(circuits55)})

    label_support = accepts["(3, 4)"]
    true_cols = support_cols(label_support, points5, 13)
    corrupt_cols = [list(c) for c in true_cols]
    corrupt_cols[0] = list(corrupt_cols[1])
    corrupt_is_circuit = rank_cols(corrupt_cols, 13) == 4 and all(rank_cols(corrupt_cols[:i] + corrupt_cols[i + 1 :], 13) == 4 for i in range(5))
    check("plants.non_tensor.labels_accept", structural_accept5(label_support, points5, 13), label_support)
    check("plants.non_tensor.actual_reject", not corrupt_is_circuit, rank_cols(corrupt_cols, 13))

    factor_cols = [[1, x] for x in points5]
    dup_cols = factor_cols + [list(factor_cols[0])]
    check("plants.duplicate.spark2", measure_spark(dup_cols, 13) == 2, measure_spark(dup_cols, 13))
    dep_pairs = [inds for inds in itertools.combinations(range(6), 2) if rank_cols([dup_cols[i] for i in inds], 13) < 2]
    check("plants.duplicate.pair", dep_pairs == [(0, 5)], dep_pairs)

    a, b = [1, 2], [3, 4]
    true_kron = [(x * y) % 13 for x in a for y in b]
    concat = a + b
    check("plants.kron_guard.dimension_collision", len(true_kron) == len(concat) == 4, {"true": true_kron, "concat": concat})
    check("plants.kron_guard.coordinate_reject", true_kron != concat and true_kron == [3, 4, 6, 8], {"true": true_kron, "concat": concat})

    example44 = next(iter({s for s in census_5_13["direct"] if profile(s) == (4, 4)}))
    d44 = decompose_44(example44, points5, 13)
    correct_law = d44["T1"] + d44["T2"] - d44["r"] - d44["c"]
    wrong_law = d44["T1"] + d44["T2"] - 2 * d44["c"]
    check("plants.general_law.correct", correct_law == 5, {"stats": {k: d44[k] for k in ("T1", "T2", "r", "c")}, "value": correct_law})
    check("plants.general_law.tight_misuse_reject", wrong_law == 6 and wrong_law != 5, wrong_law)

    RESULTS["runs"] = run_summaries
    RESULTS["plants"] = {
        "known_accepts": {k: v for k, v in accepts.items()},
        "inserted_crossing": {"support": crossing_plant, "rejected": True},
        "rank_sum_only": {"matching_candidates": 120, "circuits": 60, "rejected_as_sufficient": True},
        "non_tensor": {"support": label_support, "labels_accept": True, "actual_rejected": True},
        "duplicate": {"spark": 2, "dependent_pairs": dep_pairs},
        "kronecker_guard": {"true": true_kron, "concat": concat, "same_dimension": True, "coordinates_differ": True},
        "general_law": {"T1": d44["T1"], "T2": d44["T2"], "r": d44["r"], "c": d44["c"], "correct": correct_law, "wrong_minus_2c": wrong_law},
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
