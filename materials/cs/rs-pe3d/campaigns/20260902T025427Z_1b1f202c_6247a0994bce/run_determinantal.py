#!/usr/bin/env python3
"""Exact in-run instrument for H-3ROW-DETERMINANTAL.

No code is imported from any frozen campaign.  Direct circuit sets are computed
from exact ranks of the actual six-coordinate Kronecker columns.  Independent
predictions use graph structure, reduced Vandermonde determinants, and explicit
PGL(2,p) constructions.
"""
import sys
sys.dont_write_bytecode = True

import gc
import itertools
import json
import math
import os
import resource
import time
import traceback
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

RUN_ID = "20260902T025427Z_1b1f202c_6247a0994bce"
GATE = "H-3ROW-DETERMINANTAL"
CPU_CAP = 900.0
WALL_CAP = 1200.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT = Path(__file__).resolve().parent
CHECKS = {}
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


def encode(obj):
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(type(obj).__name__)


def check(name, condition, detail=None):
    if not condition:
        raise AssertionError(f"{name}: {detail!r}")
    CHECKS[name] = {"pass": True, "detail": detail}


def rank_mod(columns, p):
    if not columns:
        return 0
    a = [list(row) for row in zip(*columns)]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col] % p), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col] % p, -1, p)
        a[rank] = [(v * inv) % p for v in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col] % p:
                factor = a[r][col] % p
                a[r] = [(a[r][c] - factor * a[rank][c]) % p for c in range(cols)]
        rank += 1
        if rank == rows:
            break
    return rank


def det_mod(matrix, p):
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("determinant requires a square matrix")
    a = [[v % p for v in row] for row in matrix]
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        pv = a[col][col]
        det = det * pv % p
        inv = pow(pv, -1, p)
        for r in range(col + 1, n):
            if a[r][col]:
                factor = a[r][col] * inv % p
                for c in range(col, n):
                    a[r][c] = (a[r][c] - factor * a[col][c]) % p
    return det % p


def column_matrix(vectors):
    return [list(row) for row in zip(*vectors)]


def is_circuit_vectors(vectors, p):
    m = len(vectors)
    if rank_mod(vectors, p) != m - 1:
        return False
    return all(rank_mod(vectors[:k] + vectors[k + 1 :], p) == m - 1 for k in range(m))


def product_column(x, y, p):
    return (1, y % p, y * y % p, x % p, x * y % p, x * y * y % p)


def bilinear_column(x, y, p):
    return (1, y % p, x % p, x * y % p)


def d4(points, p):
    points = tuple(sorted(points))
    return det_mod(column_matrix([bilinear_column(x, y, p) for x, y in points]), p)


def vandermonde3(ys, p):
    y0, y1, y2 = ys
    return (y1 - y0) * (y2 - y0) * (y2 - y1) % p


def d6_laplace(points, p):
    points = tuple(sorted(points))
    if len(points) != 6:
        raise ValueError("D6 requires six points")
    total = 0
    indices = range(6)
    for I in itertools.combinations(indices, 3):
        Iset = set(I)
        J = tuple(j for j in indices if j not in Iset)
        vi = vandermonde3(tuple(points[i][1] for i in I), p)
        vj = vandermonde3(tuple(points[j][1] for j in J), p)
        xprod = math.prod(points[j][0] for j in J) % p
        sign = -1 if sum(i + 1 for i in I) % 2 else 1
        total = (total + sign * xprod * vi * vj) % p
    return total


def support_mask(support, cell_index):
    mask = 0
    for cell in support:
        mask |= 1 << cell_index[cell]
    return mask


def mask_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def mask_support(mask, cells):
    return tuple(cells[i] for i in mask_indices(mask))


def degree_data(support):
    da = Counter(i for i, _ in support)
    db = Counter(j for _, j in support)
    return da, db


def profile(support):
    da, db = degree_data(support)
    return len(da), len(db)


def component_signature(support):
    da, db = degree_data(support)
    adjacency = defaultdict(set)
    vertices = set()
    for i, j in support:
        a, b = (0, i), (1, j)
        adjacency[a].add(b)
        adjacency[b].add(a)
        vertices.add(a)
        vertices.add(b)
    components = []
    unseen = set(vertices)
    while unseen:
        root = min(unseen)
        stack = [root]
        comp = set()
        while stack:
            v = stack.pop()
            if v in comp:
                continue
            comp.add(v)
            unseen.discard(v)
            stack.extend(adjacency[v] - comp)
        av = sorted(v[1] for v in comp if v[0] == 0)
        bv = sorted(v[1] for v in comp if v[0] == 1)
        edges = sum(1 for i, j in support if i in av and j in bv)
        components.append((edges, len(av), len(bv), edges - len(av) - len(bv) + 1,
                           tuple(sorted((da[i] for i in av), reverse=True)),
                           tuple(sorted((db[j] for j in bv), reverse=True))))
    return (profile(support), tuple(sorted((da.values()), reverse=True)),
            tuple(sorted((db.values()), reverse=True)), tuple(sorted(components, reverse=True)))


def comb0(n, k):
    return math.comb(n, k) if n >= k >= 0 else 0


def inclusion_pair_cover(q, power):
    return sum((-1) ** t * math.comb(q, t) * comb0(q - t, 2) ** power
               for t in range(q + 1))


def canonical_pgl_matrices(p):
    maps = []
    for raw in itertools.product(range(p), repeat=4):
        a, b, c, d = raw
        if (a * d - b * c) % p == 0:
            continue
        first = next(v for v in raw if v)
        inv = pow(first, -1, p)
        normalized = tuple(v * inv % p for v in raw)
        if normalized == raw:
            maps.append(raw)
    return maps


def pgl_value(M, y, p):
    a, b, c, d = M
    den = (c + d * y) % p
    if den == 0:
        return None
    return (a + b * y) * pow(den, -1, p) % p


def build_pgl_families(points, p, cell_index):
    value_to_index = {x: i for i, x in enumerate(points)}
    matching5_mult = Counter()
    profile45_mult = Counter()
    profile55_mult = Counter()
    khist = Counter()
    maps = canonical_pgl_matrices(p)
    for M in maps:
        mapping = {}
        for j, y in enumerate(points):
            value = pgl_value(M, y, p)
            if value in value_to_index:
                mapping[j] = value_to_index[value]
        K = tuple(sorted(mapping))
        k = len(K)
        khist[k] += 1
        for Q in itertools.combinations(K, 5):
            support = tuple(sorted((mapping[j], j) for j in Q))
            matching5_mult[support_mask(support, cell_index)] += 1
        for Q in itertools.combinations(K, 4):
            Qset = set(Q)
            base = {(mapping[j], j) for j in Q}
            outside_y = [j for j in range(len(points)) if j not in Qset]
            for y0 in outside_y:
                for chosen in itertools.combinations(Q, 2):
                    support = set(base)
                    support |= {(mapping[j], y0) for j in chosen}
                    profile45_mult[support_mask(tuple(sorted(support)), cell_index)] += 1
                image_Q = {mapping[j] for j in Q}
                for distinguished in Q:
                    for xnew in range(len(points)):
                        if xnew in image_Q:
                            continue
                        if y0 in mapping and xnew == mapping[y0]:
                            continue
                        support = set(base)
                        support.add((mapping[distinguished], y0))
                        support.add((xnew, y0))
                        profile55_mult[support_mask(tuple(sorted(support)), cell_index)] += 1
    return {
        "maps": maps,
        "k_histogram": khist,
        "matching5_mult": matching5_mult,
        "profile45_mult": profile45_mult,
        "profile55_mult": profile55_mult,
    }


def construct_crossings(n, cell_index):
    out = set()
    for rows in itertools.combinations(range(n), 3):
        for center_row in rows:
            for cols in itertools.combinations(range(n), 4):
                for center_col in cols:
                    support = {(i, center_col) for i in rows if i != center_row}
                    support |= {(center_row, j) for j in cols if j != center_col}
                    out.add(support_mask(tuple(sorted(support)), cell_index))
    return out


def construct_q4_size6(n, cell_index):
    by_profile = defaultdict(set)
    mult = Counter()
    for psize in range(3, min(5, n) + 1):
        for rows in itertools.combinations(range(n), psize):
            for cols in itertools.combinations(range(n), 4):
                for singleton_cols in itertools.combinations(cols, 2):
                    doubled_cols = [j for j in cols if j not in singleton_cols]
                    for singleton_row in rows:
                        other_rows = [i for i in rows if i != singleton_row]
                        for pair0 in itertools.combinations(other_rows, 2):
                            for pair1 in itertools.combinations(other_rows, 2):
                                if set(pair0) | set(pair1) != set(other_rows):
                                    continue
                                support = {(singleton_row, j) for j in singleton_cols}
                                support |= {(i, doubled_cols[0]) for i in pair0}
                                support |= {(i, doubled_cols[1]) for i in pair1}
                                mask = support_mask(tuple(sorted(support)), cell_index)
                                by_profile[(psize, 4)].add(mask)
                                mult[mask] += 1
    return by_profile, mult


def construct_q4_size7(n, cell_index):
    by_profile = defaultdict(set)
    mult = Counter()
    for psize in range(3, min(7, n) + 1):
        for rows in itertools.combinations(range(n), psize):
            for cols in itertools.combinations(range(n), 4):
                for singleton_col in cols:
                    doubled_cols = [j for j in cols if j != singleton_col]
                    for singleton_row in rows:
                        other_rows = [i for i in rows if i != singleton_row]
                        pair_choices = list(itertools.combinations(other_rows, 2))
                        for pairs in itertools.product(pair_choices, repeat=3):
                            if set().union(*map(set, pairs)) != set(other_rows):
                                continue
                            support = {(singleton_row, singleton_col)}
                            for j, pair in zip(doubled_cols, pairs):
                                support |= {(i, j) for i in pair}
                            mask = support_mask(tuple(sorted(support)), cell_index)
                            by_profile[(psize, 4)].add(mask)
                            mult[mask] += 1
    return by_profile, mult


def q5_singletons(support):
    _, db = degree_data(support)
    doubled = [j for j, d in db.items() if d == 2]
    singletons = [(i, j) for i, j in support if db[j] == 1]
    return doubled, singletons


def dependent6_reduced(support, p):
    da, db = degree_data(support)
    q = len(db)
    if q == 4 and sorted(db.values()) == [1, 1, 2, 2]:
        singleton_cells = [(i, j) for i, j in support if db[j] == 1]
        return singleton_cells[0][0] == singleton_cells[1][0]
    if q == 5 and sorted(db.values()) == [1, 1, 1, 1, 2]:
        _, singletons = q5_singletons(support)
        return d4(singletons, p) == 0
    if q == 6:
        return d6_laplace(support, p) == 0
    return False


def predict_size7(support, p):
    da, db = degree_data(support)
    if len(da) < 3 or len(db) < 4 or max(da.values()) > 3 or max(db.values()) > 2:
        return False
    q = len(db)
    if q == 4:
        if sorted(db.values()) != [1, 2, 2, 2]:
            return False
        singleton = next((i, j) for i, j in support if db[j] == 1)
        r0 = singleton[0]
        return all(i != r0 for i, j in support if db[j] == 2)
    if q == 5:
        if sorted(db.values()) != [1, 1, 1, 2, 2]:
            return False
        doubled, singletons = q5_singletons(support)
        if len({i for i, _ in singletons}) != 3:
            return False
        for j in doubled:
            for cell in (cell for cell in support if cell[1] == j):
                if d4(singletons + [cell], p) == 0:
                    return False
        return True
    if q == 6:
        if sorted(db.values()) != [1, 1, 1, 1, 1, 2]:
            return False
        doubled, singletons = q5_singletons(support)
        for deleted in singletons:
            remaining_singletons = [cell for cell in singletons if cell != deleted]
            if d4(remaining_singletons, p) == 0:
                return False
        doubled_cells = [cell for cell in support if cell[1] == doubled[0]]
        for deleted in doubled_cells:
            remaining = tuple(cell for cell in support if cell != deleted)
            if d6_laplace(remaining, p) == 0:
                return False
        return True
    if q == 7:
        return all(d6_laplace(tuple(cell for cell in support if cell != deleted), p) != 0
                   for deleted in support)
    return False


def census_configuration(n, p):
    label = f"n{n}.p{p}"
    points = tuple(range(1, n + 1))
    cells = tuple((i, j) for i in range(n) for j in range(n))
    cell_index = {cell: k for k, cell in enumerate(cells)}
    columns = tuple(product_column(points[i], points[j], p) for i, j in cells)

    @lru_cache(maxsize=None)
    def rank_mask(mask):
        return rank_mod([columns[i] for i in mask_indices(mask)], p)

    direct = {5: set(), 6: set(), 7: set()}
    factor_signs = defaultdict(set)
    factor_checked = 0
    laplace_checked = 0

    for m in (5, 6):
        for combo in itertools.combinations(range(n * n), m):
            mask = sum(1 << i for i in combo)
            rank = rank_mask(mask)
            if rank == m - 1 and all(rank_mask(mask ^ (1 << i)) == m - 1 for i in combo):
                direct[m].add(mask)
            if m == 6:
                support = tuple(cells[i] for i in combo)
                direct_det = det_mod(column_matrix([columns[i] for i in combo]), p)
                laplace_det = d6_laplace(tuple((points[i], points[j]) for i, j in support), p)
                if direct_det != laplace_det:
                    raise AssertionError(
                        f"{label}.laplace.{mask}: direct={direct_det}, laplace={laplace_det}"
                    )
                laplace_checked += 1
                _, db = degree_data(support)
                if len(db) == 5 and sorted(db.values()) == [1, 1, 1, 1, 2]:
                    doubled_j = next(j for j, d in db.items() if d == 2)
                    duplicated = [cell for cell in support if cell[1] == doubled_j]
                    singleton_cells = [cell for cell in support if cell[1] != doubled_j]
                    xa, xb = (points[duplicated[0][0]], points[duplicated[1][0]])
                    y0 = points[doubled_j]
                    singleton_points = [(points[i], points[j]) for i, j in singleton_cells]
                    rhs = (xb - xa) % p
                    for _, y in singleton_points:
                        rhs = rhs * (y - y0) % p
                    rhs = rhs * d4(singleton_points, p) % p
                    if rhs == 0:
                        if laplace_det != 0:
                            raise AssertionError(
                                f"{label}.factor.zero.{mask}: D6={laplace_det}, rhs={rhs}"
                            )
                    else:
                        ratio = laplace_det * pow(rhs, -1, p) % p
                        if ratio not in (1, p - 1):
                            raise AssertionError(f"{label}.factor.unit.{mask}: ratio={ratio}")
                        positions = tuple(k for k, cell in enumerate(support) if cell[1] == doubled_j)
                        factor_signs[positions].add(ratio)
                    factor_checked += 1
        check_budget(f"{label}.m{m}")
    check(f"{label}.laplace.exhaustive", laplace_checked == math.comb(n * n, 6),
          {"checked": laplace_checked, "universe": math.comb(n * n, 6)})
    factor_ok = factor_checked == 0 if n < 5 else factor_checked > 0 and all(len(v) == 1 for v in factor_signs.values())
    check(f"{label}.factor.exhaustive", factor_ok,
          {"checked": factor_checked, "sign_classes": len(factor_signs), "applicable": n >= 5})

    for combo in itertools.combinations(range(n * n), 7):
        mask = sum(1 << i for i in combo)
        if all(rank_mask(mask ^ (1 << i)) == 6 for i in combo):
            direct[7].add(mask)
    check_budget(f"{label}.m7")

    crossings = construct_crossings(n, cell_index)
    pgl = build_pgl_families(points, p, cell_index)
    matching5 = set(pgl["matching5_mult"])
    predicted5 = crossings | matching5
    check(f"{label}.pgl.matching.injective", all(v == 1 for v in pgl["matching5_mult"].values()),
          max(pgl["matching5_mult"].values(), default=0))
    check(f"{label}.size5.set_equality", direct[5] == predicted5,
          {"direct": len(direct[5]), "predicted": len(predicted5),
           "direct_only": len(direct[5] - predicted5), "predicted_only": len(predicted5 - direct[5])})
    formula_cross = 12 * math.comb(n, 3) * math.comb(n, 4)
    formula_match = sum(count * comb0(k, 5) for k, count in pgl["k_histogram"].items())
    check(f"{label}.size5.crossing.formula", len(crossings) == formula_cross,
          {"got": len(crossings), "want": formula_cross})
    check(f"{label}.size5.matching.formula", len(matching5) == formula_match,
          {"got": len(matching5), "want": formula_match})

    # Barycentric kernel identity on every five-point B set available here.
    barycentric_checked = 0
    for ys in itertools.combinations(points, 5):
        weights = []
        for y in ys:
            denom = math.prod((y - z) % p for z in ys if z != y) % p
            weights.append(pow(denom, -1, p))
        for exponent in range(4):
            check(f"{label}.barycentric.{ys}.{exponent}",
                  sum(w * pow(y, exponent, p) for w, y in zip(weights, ys)) % p == 0,
                  {"ys": ys, "exponent": exponent})
        barycentric_checked += 1

    predicted6 = set()
    wrong55 = set()
    for combo in itertools.combinations(range(n * n), 6):
        mask = sum(1 << i for i in combo)
        support = tuple(cells[i] for i in combo)
        da, db = degree_data(support)
        if len(da) < 3 or len(db) < 4 or max(da.values()) > 3 or max(db.values()) > 2:
            continue
        dependency = dependent6_reduced(tuple((points[i], points[j]) for i, j in support), p)
        if not dependency:
            continue
        if profile(support) == (5, 5):
            degree2_a = [i for i, d in da.items() if d == 2]
            degree2_b = [j for j, d in db.items() if d == 2]
            if len(degree2_a) == len(degree2_b) == 1 and (degree2_a[0], degree2_b[0]) in support:
                wrong55.add(mask)
        if all((mask ^ (1 << i)) not in predicted5 for i in combo):
            predicted6.add(mask)
    check(f"{label}.size6.set_equality", direct[6] == predicted6,
          {"direct": len(direct[6]), "predicted": len(predicted6),
           "direct_only": len(direct[6] - predicted6), "predicted_only": len(predicted6 - direct[6])})

    constructed45 = set(pgl["profile45_mult"])
    constructed55 = set(pgl["profile55_mult"])
    direct45 = {mask for mask in direct[6] if profile(mask_support(mask, cells)) == (4, 5)}
    direct55 = {mask for mask in direct[6] if profile(mask_support(mask, cells)) == (5, 5)}
    check(f"{label}.pgl.45.injective", all(v == 1 for v in pgl["profile45_mult"].values()),
          max(pgl["profile45_mult"].values(), default=0))
    check(f"{label}.pgl.55.injective", all(v == 1 for v in pgl["profile55_mult"].values()),
          max(pgl["profile55_mult"].values(), default=0))
    check(f"{label}.size6.45.set_equality", direct45 == constructed45,
          {"direct": len(direct45), "constructed": len(constructed45)})
    check(f"{label}.size6.55.set_equality", direct55 == constructed55,
          {"direct": len(direct55), "constructed": len(constructed55)})
    formula45 = 6 * (n - 4) * sum(count * comb0(k, 4) for k, count in pgl["k_histogram"].items())
    formula55 = 4 * sum(count * comb0(k, 4) * ((n - 4) ** 2 - (k - 4))
                       for k, count in pgl["k_histogram"].items())
    check(f"{label}.size6.45.formula", len(constructed45) == formula45,
          {"got": len(constructed45), "want": formula45})
    check(f"{label}.size6.55.formula", len(constructed55) == formula55,
          {"got": len(constructed55), "want": formula55})

    q4_6, q4_6_mult = construct_q4_size6(n, cell_index)
    constructed_q4_6 = set().union(*q4_6.values()) if q4_6 else set()
    direct_q4_6 = {mask for mask in direct[6] if profile(mask_support(mask, cells))[1] == 4}
    check(f"{label}.q4.size6.injective", all(v == 1 for v in q4_6_mult.values()),
          max(q4_6_mult.values(), default=0))
    check(f"{label}.q4.size6.set_equality", direct_q4_6 == constructed_q4_6,
          {"direct": len(direct_q4_6), "constructed": len(constructed_q4_6)})
    q4_6_formula = {}
    for psize in range(3, min(5, n) + 1):
        h = inclusion_pair_cover(psize - 1, 2)
        want = 6 * psize * math.comb(n, psize) * math.comb(n, 4) * h
        got = len(q4_6[(psize, 4)])
        check(f"{label}.q4.size6.profile{psize}.formula", got == want, {"got": got, "want": want, "h": h})
        q4_6_formula[str((psize, 4))] = {"count": got, "formula": want, "h": h}

    predicted7 = set()
    for combo in itertools.combinations(range(n * n), 7):
        support = tuple(cells[i] for i in combo)
        point_support = tuple((points[i], points[j]) for i, j in support)
        if predict_size7(point_support, p):
            predicted7.add(sum(1 << i for i in combo))
    check(f"{label}.size7.set_equality", direct[7] == predicted7,
          {"direct": len(direct[7]), "predicted": len(predicted7),
           "direct_only": len(direct[7] - predicted7), "predicted_only": len(predicted7 - direct[7])})

    direct_signatures7 = defaultdict(set)
    predicted_signatures7 = defaultdict(set)
    for mask in direct[7]:
        direct_signatures7[str(component_signature(mask_support(mask, cells)))].add(mask)
    for mask in predicted7:
        predicted_signatures7[str(component_signature(mask_support(mask, cells)))].add(mask)
    check(f"{label}.size7.signature.keys", set(direct_signatures7) == set(predicted_signatures7),
          {"direct": len(direct_signatures7), "predicted": len(predicted_signatures7)})
    for signature in direct_signatures7:
        check(f"{label}.size7.signature.{signature}",
              direct_signatures7[signature] == predicted_signatures7[signature],
              len(direct_signatures7[signature]))

    q4_7, q4_7_mult = construct_q4_size7(n, cell_index)
    constructed_q4_7 = set().union(*q4_7.values()) if q4_7 else set()
    direct_q4_7 = {mask for mask in direct[7] if profile(mask_support(mask, cells))[1] == 4}
    check(f"{label}.q4.size7.injective", all(v == 1 for v in q4_7_mult.values()),
          max(q4_7_mult.values(), default=0))
    check(f"{label}.q4.size7.set_equality", direct_q4_7 == constructed_q4_7,
          {"direct": len(direct_q4_7), "constructed": len(constructed_q4_7)})
    q4_7_formula = {}
    for psize in range(3, min(7, n) + 1):
        g = inclusion_pair_cover(psize - 1, 3)
        want = 4 * psize * math.comb(n, psize) * math.comb(n, 4) * g
        got = len(q4_7[(psize, 4)])
        check(f"{label}.q4.size7.profile{psize}.formula", got == want, {"got": got, "want": want, "g": g})
        q4_7_formula[str((psize, 4))] = {"count": got, "formula": want, "g": g}

    # The wrong predicate omits the central matching deletion at profile (5,5).
    check(f"{label}.wrong55.contains_right", direct55 <= wrong55,
          {"right": len(direct55), "wrong": len(wrong55)})
    if n == 5:
        identity_R = {(i, i) for i in range(4)}
        wrong_plant = identity_R | {(0, 4), (4, 4)}
        wrong_mask = support_mask(tuple(sorted(wrong_plant)), cell_index)
        check(f"{label}.wrong55.plant_wrong_accept", wrong_mask in wrong55, sorted(wrong_plant))
        check(f"{label}.wrong55.plant_actual_reject", wrong_mask not in direct[6], sorted(wrong_plant))
        check(f"{label}.wrong55.set_inequality", wrong55 != direct55,
              {"overaccept": len(wrong55 - direct55), "underaccept": len(direct55 - wrong55)})

    # Direct template assertions for the two named residual six-channels.
    expected45 = ((4, 5), (2, 2, 1, 1), (2, 1, 1, 1, 1),
                  ((4, 2, 3, 0, (2, 2), (2, 1, 1)),
                   (1, 1, 1, 0, (1,), (1,)), (1, 1, 1, 0, (1,), (1,))))
    expected55 = ((5, 5), (2, 1, 1, 1, 1), (2, 1, 1, 1, 1),
                  ((3, 2, 2, 0, (2, 1), (2, 1)),
                   (1, 1, 1, 0, (1,), (1,)), (1, 1, 1, 0, (1,), (1,)),
                   (1, 1, 1, 0, (1,), (1,))))
    check(f"{label}.template45", all(component_signature(mask_support(mask, cells)) == expected45 for mask in direct45), len(direct45))
    check(f"{label}.template55", all(component_signature(mask_support(mask, cells)) == expected55 for mask in direct55), len(direct55))

    summary_sizes = {}
    for m in (5, 6, 7):
        profile_counts = Counter(profile(mask_support(mask, cells)) for mask in direct[m])
        signature_counts = Counter(str(component_signature(mask_support(mask, cells))) for mask in direct[m])
        summary_sizes[str(m)] = {
            "universe": math.comb(n * n, m),
            "circuit_count": len(direct[m]),
            "profiles": {str(k): v for k, v in sorted(profile_counts.items())},
            "signature_counts": dict(sorted(signature_counts.items())),
            "direct_predicted_set_equality": True,
        }

    record = {
        "n": n,
        "prime": p,
        "sizes": summary_sizes,
        "pgl": {
            "group_size": len(pgl["maps"]),
            "k_histogram": {str(k): v for k, v in sorted(pgl["k_histogram"].items())},
            "size5_matching": {"count": len(matching5), "formula": formula_match, "set_equality": True},
            "size6_profile45": {"count": len(constructed45), "formula": formula45, "set_equality": True},
            "size6_profile55": {"count": len(constructed55), "formula": formula55, "set_equality": True},
        },
        "q4_size6": q4_6_formula,
        "q4_size7": q4_7_formula,
        "factorization": {
            "laplace_six_sets_checked": laplace_checked,
            "one_repeated_B_sets_checked": factor_checked,
            "sign_by_duplicate_positions": {str(k): sorted(v) for k, v in sorted(factor_signs.items())},
        },
        "wrong_profile55": {"wrong_count": len(wrong55), "circuit_count": len(direct55), "overaccept": len(wrong55 - direct55)},
        "barycentric_five_sets_checked": barycentric_checked,
        "rank_cache": rank_mask.cache_info()._asdict(),
    }
    retained = {m: set(direct[m]) for m in direct}
    del rank_mask
    gc.collect()
    return record, retained, cells, columns


def cross_prime_comparison(n, entries, cells):
    out = {"n": n, "pairs": []}
    for (p1, sets1), (p2, sets2) in itertools.combinations(entries, 2):
        pair = {"primes": [p1, p2], "sizes": {}}
        for m in (5, 6, 7):
            profiles = sorted({profile(mask_support(mask, cells)) for mask in sets1[m] | sets2[m]})
            by_profile = {}
            for pr in profiles:
                A = {mask for mask in sets1[m] if profile(mask_support(mask, cells)) == pr}
                B = {mask for mask in sets2[m] if profile(mask_support(mask, cells)) == pr}
                by_profile[str(pr)] = {
                    "left": len(A), "right": len(B), "intersection": len(A & B),
                    "left_only": len(A - B), "right_only": len(B - A), "sets_equal": A == B,
                }
            pair["sizes"][str(m)] = by_profile
        out["pairs"].append(pair)
    return out


def run_plants(census_by_key):
    p = 13
    n = 5
    record, direct, cells, columns = census_by_key[(n, p)]
    cell_index = {cell: i for i, cell in enumerate(cells)}

    crossings = construct_crossings(n, cell_index)
    crossing_mask = min(crossings)
    crossing_support = mask_support(crossing_mask, cells)
    crossing_vectors = [columns[i] for i in mask_indices(crossing_mask)]
    check("plants.crossing_accept", is_circuit_vectors(crossing_vectors, p), crossing_support)

    identity = tuple((i, i) for i in range(5))
    identity_mask = support_mask(identity, cell_index)
    identity_vectors = [columns[i] for i in mask_indices(identity_mask)]
    check("plants.mobius_matching_accept", identity_mask in direct[5] and is_circuit_vectors(identity_vectors, p), identity)

    added = next(cell for cell in cells if cell not in identity)
    inserted = tuple(sorted(identity + (added,)))
    inserted_vectors = [columns[cell_index[cell]] for cell in inserted]
    check("plants.inserted_dependent_reject", not is_circuit_vectors(inserted_vectors, p), inserted)

    nonminimal7 = tuple((0, j) for j in range(4)) + ((1, 0), (2, 1), (3, 2))
    nonminimal7_vectors = [columns[cell_index[cell]] for cell in nonminimal7]
    check("plants.ambient7_automatic_dependence", rank_mod(nonminimal7_vectors, p) <= 6,
          rank_mod(nonminimal7_vectors, p))
    check("plants.ambient7_nonminimal_reject", not is_circuit_vectors(nonminimal7_vectors, p), nonminimal7)

    accept7_mask = min(direct[7])
    accept7_support = mask_support(accept7_mask, cells)
    accept7_vectors = [columns[i] for i in mask_indices(accept7_mask)]
    check("plants.ambient7_actual_accept", is_circuit_vectors(accept7_vectors, p), accept7_support)

    true = product_column(2, 3, p)
    concat = (1, 2, 1, 3, 9)
    check("plants.kron_guard.true_dimension", len(true) == 6, true)
    check("plants.kron_guard.concat_reject", len(concat) == 5 and concat != true, {"true": true, "concat": concat})

    profile45 = [mask for mask in direct[6] if profile(mask_support(mask, cells)) == (4, 5)]
    check("plants.profile45_population", bool(profile45), len(profile45))
    chosen_mask = profile45[0]
    support = mask_support(chosen_mask, cells)
    original = [list(columns[i]) for i in mask_indices(chosen_mask)]
    corrupted = None
    for coordinate in range(6):
        for delta in range(1, p):
            candidate = [row[:] for row in original]
            candidate[0][coordinate] = (candidate[0][coordinate] + delta) % p
            vectors = [tuple(row) for row in candidate]
            if not is_circuit_vectors(vectors, p):
                corrupted = {"coordinate": coordinate, "delta": delta, "rank": rank_mod(vectors, p)}
                break
        if corrupted:
            break
    check("plants.corrupted_actual_reject", corrupted is not None, {"support": support, "corruption": corrupted})

    wrong_support = tuple(sorted({(i, i) for i in range(4)} | {(0, 4), (4, 4)}))
    wrong_points = tuple((i + 1, j + 1) for i, j in wrong_support)
    _, wrong_singletons = q5_singletons(wrong_points)
    check("plants.wrong_minor_accepts", d4(wrong_singletons, p) == 0, wrong_support)
    wrong_vectors = [columns[cell_index[cell]] for cell in wrong_support]
    check("plants.wrong_minor_actual_reject", not is_circuit_vectors(wrong_vectors, p), wrong_support)

    return {
        "crossing_accept": list(crossing_support),
        "mobius_matching_accept": list(identity),
        "inserted_dependent_reject": list(inserted),
        "ambient7_accept": list(accept7_support),
        "ambient7_nonminimal_reject": list(nonminimal7),
        "wrong_minor_reject": list(wrong_support),
        "corrupted_actual_reject": corrupted,
        "kron_guard": {"true": true, "concat": concat},
    }


def main():
    print(f"== {GATE} run {RUN_ID}")
    nice_level = os.nice(0)
    check("runtime.nice", nice_level == 15, nice_level)
    check("runtime.sys_dont_write_bytecode", sys.dont_write_bytecode is True, sys.dont_write_bytecode)
    check("runtime.env_dont_write_bytecode", os.environ.get("PYTHONDONTWRITEBYTECODE") == "1",
          os.environ.get("PYTHONDONTWRITEBYTECODE"))
    thread_caps = {name: os.environ.get(name) for name in
                   ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    check("runtime.thread_caps", all(value == "1" for value in thread_caps.values()), thread_caps)

    summaries = []
    census_by_key = {}
    retained_by_n = defaultdict(list)
    for n in (4, 5):
        for p in (7, 11, 13):
            record, direct, cells, columns = census_configuration(n, p)
            summaries.append(record)
            census_by_key[(n, p)] = (record, direct, cells, columns)
            retained_by_n[n].append((p, direct))
            counts = {m: len(direct[m]) for m in (5, 6, 7)}
            print(f"[census] n={n}, p={p}: {counts}")
            check_budget(f"after n{n}p{p}")

    comparisons = []
    for n in (4, 5):
        cells = census_by_key[(n, 7)][2]
        comparisons.append(cross_prime_comparison(n, retained_by_n[n], cells))

    plants = run_plants(census_by_key)
    cpu, wall = elapsed()
    check("budget.CPU", cpu <= CPU_CAP, {"used": cpu, "cap": CPU_CAP})
    check("budget.wall", wall <= WALL_CAP, {"used": wall, "cap": WALL_CAP})

    RESULTS.update({
        "runtime": {
            "nice": nice_level,
            "dont_write_bytecode": sys.dont_write_bytecode,
            "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
            "thread_caps": thread_caps,
            "max_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        },
        "budget": {"CPU_seconds": cpu, "wall_seconds": wall, "CPU_cap_seconds": CPU_CAP,
                   "wall_cap_seconds": WALL_CAP, "within_cap": True},
        "configurations": summaries,
        "prime_comparisons": comparisons,
        "plants": plants,
        "controls": CHECKS,
        "assertion_count": len(CHECKS),
    })
    (OUT / "controls_results.json").write_text(json.dumps(RESULTS, indent=2, sort_keys=True, default=encode) + "\n")
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
        }
        path = OUT / "broken_results.json"
        if path.exists():
            path = OUT / f"broken_results_{int(time.time())}.json"
        path.write_text(json.dumps(broken, indent=2, sort_keys=True, default=encode) + "\n")
        print(f"BROKEN output preserved at {path.name}", file=sys.stderr)
        raise
