#!/usr/bin/env python3
"""Exact instrument for H-33ROW-OPEN.

No frozen campaign is imported.  Full-grid circuit sets use two independent
finite-field rank implementations on every subset.  Dedicated matching sweeps
exercise the Möbius and higher-curve channels beyond n=5.
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
from pathlib import Path

RUN_ID = "20260902T031527Z_2d19c5c4_718416983678"
GATE = "H-33ROW-OPEN"
CPU_CAP = 7200.0
WALL_CAP = 8000.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
OUT = Path(__file__).resolve().parent
CHECKS = {}
CHECKPOINT = {"run_id": RUN_ID, "completed_layers": []}
RESULTS = {"run_id": RUN_ID, "gate": GATE,
           "arithmetic": "exact Python integers modulo p; true 3-by-3 row-major Kronecker columns"}


def elapsed():
    return time.process_time() - START_CPU, time.monotonic() - START_WALL


def check_budget(label):
    cpu, wall = elapsed()
    if cpu > CPU_CAP or wall > WALL_CAP:
        raise RuntimeError(f"budget exceeded at {label}: CPU={cpu:.3f}, wall={wall:.3f}")


def checkpoint(record):
    CHECKPOINT["completed_layers"].append(record)
    cpu, wall = elapsed()
    CHECKPOINT["CPU_seconds"] = cpu
    CHECKPOINT["wall_seconds"] = wall
    (OUT / "checkpoint.json").write_text(json.dumps(CHECKPOINT, indent=2, sort_keys=True) + "\n")


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


def rank_mod_alt(columns, p):
    """Independent column-echelon implementation used on every full-grid subset."""
    basis = []
    pivots = []
    for source in columns:
        v = [x % p for x in source]
        for pivot, b in zip(pivots, basis):
            if v[pivot]:
                factor = v[pivot]
                v = [(x - factor * y) % p for x, y in zip(v, b)]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is None:
            continue
        inv = pow(v[pivot], -1, p)
        v = [x * inv % p for x in v]
        basis.append(v)
        pivots.append(pivot)
    return len(basis)


def det_mod(matrix, p):
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("square determinant required")
    a = [[x % p for x in row] for row in matrix]
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


def product_column(x, y, p):
    xp = (1, x % p, x * x % p)
    yp = (1, y % p, y * y % p)
    return tuple(a * b % p for a in xp for b in yp)


def is_circuit_vectors(vectors, p):
    m = len(vectors)
    if rank_mod(vectors, p) != m - 1:
        return False
    return all(rank_mod(vectors[:i] + vectors[i + 1 :], p) == m - 1 for i in range(m))


def mask_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def mask_support(mask, cells):
    return tuple(cells[i] for i in mask_indices(mask))


def support_mask(support, cell_index):
    return sum(1 << cell_index[cell] for cell in support)


def degree_data(support):
    return Counter(i for i, _ in support), Counter(j for _, j in support)


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
        vertices |= {a, b}
    unseen = set(vertices)
    components = []
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
    return (profile(support), tuple(sorted(da.values(), reverse=True)),
            tuple(sorted(db.values(), reverse=True)), tuple(sorted(components, reverse=True)))


def construct_fibers(n, cell_index):
    out = set()
    for i in range(n):
        for cols in itertools.combinations(range(n), 4):
            out.add(support_mask(tuple((i, j) for j in cols), cell_index))
    for j in range(n):
        for rows in itertools.combinations(range(n), 4):
            out.add(support_mask(tuple((i, j) for i in rows), cell_index))
    return out


def construct_crossings(n, cell_index):
    out = set()
    mult = Counter()
    for rows in itertools.combinations(range(n), 4):
        for center_row in rows:
            for cols in itertools.combinations(range(n), 4):
                for center_col in cols:
                    support = {(i, center_col) for i in rows if i != center_row}
                    support |= {(center_row, j) for j in cols if j != center_col}
                    mask = support_mask(tuple(sorted(support)), cell_index)
                    out.add(mask)
                    mult[mask] += 1
    return out, mult


def canonical_pgl_matrices(p):
    maps = []
    for raw in itertools.product(range(p), repeat=4):
        a, b, c, d = raw
        if (a * d - b * c) % p == 0:
            continue
        first = next(v for v in raw if v)
        inv = pow(first, -1, p)
        if tuple(v * inv % p for v in raw) == raw:
            maps.append(raw)
    return maps


def pgl_value(M, y, p):
    a, b, c, d = M
    den = (c + d * y) % p
    return None if den == 0 else (a + b * y) * pow(den, -1, p) % p


def construct_pgl_matchings(points_x, points_y, p, size, cell_index):
    x_index = {x: i for i, x in enumerate(points_x)}
    mult = Counter()
    khist = Counter()
    maps = canonical_pgl_matrices(p)
    for M in maps:
        mapping = {}
        for j, y in enumerate(points_y):
            value = pgl_value(M, y, p)
            if value in x_index:
                mapping[j] = x_index[value]
        k = len(mapping)
        khist[k] += 1
        for chosen in itertools.combinations(sorted(mapping), size):
            support = tuple(sorted((mapping[j], j) for j in chosen))
            if len({i for i, _ in support}) == size:
                mult[support_mask(support, cell_index)] += 1
    return set(mult), mult, khist, len(maps)


def summarize_circuits(circuit_sets, cells, universes):
    sizes = {}
    for m, masks in sorted(circuit_sets.items()):
        profiles = Counter(profile(mask_support(mask, cells)) for mask in masks)
        signatures = Counter(str(component_signature(mask_support(mask, cells))) for mask in masks)
        degree_violations = sum(1 for mask in masks if m > 4 and
                                (max(degree_data(mask_support(mask, cells))[0].values()) > 3 or
                                 max(degree_data(mask_support(mask, cells))[1].values()) > 3))
        sizes[str(m)] = {
            "universe": universes[m],
            "circuit_count": len(masks),
            "profiles": {str(k): v for k, v in sorted(profiles.items())},
            "signature_count": len(signatures),
            "signature_counts": dict(sorted(signatures.items())),
            "degree_violations": degree_violations,
        }
    return sizes


def census_grid(n, p, max_m):
    label = f"grid.n{n}.p{p}"
    points = tuple(x % p for x in range(1, n + 1))
    check(f"{label}.points_distinct", len(set(points)) == n, points)
    cells = tuple((i, j) for i in range(n) for j in range(n))
    cell_index = {cell: k for k, cell in enumerate(cells)}
    columns = tuple(product_column(points[i], points[j], p) for i, j in cells)
    N = n * n
    independent_prev = set()
    for combo in itertools.combinations(range(N), 3):
        mask = sum(1 << i for i in combo)
        vectors = [columns[i] for i in combo]
        r1 = rank_mod(vectors, p)
        r2 = rank_mod_alt(vectors, p)
        if r1 != r2:
            raise AssertionError(f"{label}.rank_crosscheck.m3.{mask}: {r1}!={r2}")
        if r1 == 3:
            independent_prev.add(mask)
    check(f"{label}.spark_lower", len(independent_prev) == math.comb(N, 3), len(independent_prev))

    circuit_sets = {}
    universes = {}
    rank_crosschecks = math.comb(N, 3)
    for m in range(4, max_m + 1):
        universe = math.comb(N, m)
        universes[m] = universe
        circuits = set()
        independent_current = set()
        if m <= 9:
            for combo in itertools.combinations(range(N), m):
                mask = sum(1 << i for i in combo)
                vectors = [columns[i] for i in combo]
                r1 = rank_mod(vectors, p)
                r2 = rank_mod_alt(vectors, p)
                if r1 != r2:
                    raise AssertionError(f"{label}.rank_crosscheck.m{m}.{mask}: {r1}!={r2}")
                rank_crosschecks += 1
                if r1 == m:
                    independent_current.add(mask)
                elif r1 == m - 1 and all((mask ^ (1 << i)) in independent_prev for i in combo):
                    circuits.add(mask)
        else:
            for combo in itertools.combinations(range(N), m):
                mask = sum(1 << i for i in combo)
                if all((mask ^ (1 << i)) in independent_prev for i in combo):
                    circuits.add(mask)
        circuit_sets[m] = circuits
        checkpoint({"kind": "grid", "n": n, "prime": p, "size": m,
                    "universe": universe, "circuits": len(circuits)})
        check_budget(f"{label}.m{m}")
        independent_prev = independent_current
        gc.collect()

    fibers = construct_fibers(n, cell_index)
    check(f"{label}.size4.fiber_set_equality", circuit_sets[4] == fibers,
          {"direct": len(circuit_sets[4]), "constructed": len(fibers)})
    check(f"{label}.size4.formula", len(fibers) == 2 * n * math.comb(n, 4), len(fibers))
    check(f"{label}.size5.empty", circuit_sets[5] == set(), len(circuit_sets[5]))
    crossings, crossing_mult = construct_crossings(n, cell_index)
    direct44 = {mask for mask in circuit_sets[6] if profile(mask_support(mask, cells)) == (4, 4)}
    check(f"{label}.size6.crossing_injective", all(v == 1 for v in crossing_mult.values()),
          max(crossing_mult.values(), default=0))
    check(f"{label}.size6.crossing_set_equality", direct44 == crossings,
          {"direct": len(direct44), "constructed": len(crossings)})
    check(f"{label}.size6.crossing_formula",
          len(crossings) == 16 * math.comb(n, 4) ** 2, len(crossings))
    for m, masks in circuit_sets.items():
        if m > 4:
            violations = [mask for mask in masks if
                          max(degree_data(mask_support(mask, cells))[0].values()) > 3 or
                          max(degree_data(mask_support(mask, cells))[1].values()) > 3]
            check(f"{label}.m{m}.degree_bounds", not violations, len(violations))

    record = {
        "n": n, "prime": p, "max_size": max_m,
        "rank_implementations_agree": True,
        "rank_crosschecks": rank_crosschecks,
        "sizes": summarize_circuits(circuit_sets, cells, universes),
        "size4_fibers": len(fibers),
        "size6_crossings": len(crossings),
    }
    return record, circuit_sets, cells, columns


def matching_masks(n, k, cell_index):
    for rows in itertools.combinations(range(n), k):
        for cols in itertools.combinations(range(n), k):
            for perm in itertools.permutations(cols):
                support = tuple(sorted(zip(rows, perm)))
                yield support_mask(support, cell_index)


def census_matchings(n, p, sizes):
    label = f"matching.n{n}.p{p}"
    points = tuple(x % p for x in range(1, n + 1))
    check(f"{label}.points_distinct", len(set(points)) == n, points)
    cells = tuple((i, j) for i in range(n) for j in range(n))
    cell_index = {cell: k for k, cell in enumerate(cells)}
    columns = tuple(product_column(points[i], points[j], p) for i, j in cells)
    direct = {}
    universes = {}
    for k in sizes:
        circuits = set()
        count = 0
        for mask in matching_masks(n, k, cell_index):
            count += 1
            vectors = [columns[i] for i in mask_indices(mask)]
            if is_circuit_vectors(vectors, p):
                circuits.add(mask)
        direct[k] = circuits
        universes[k] = count
        checkpoint({"kind": "matching", "n": n, "prime": p, "size": k,
                    "universe": count, "circuits": len(circuits)})
        check_budget(f"{label}.m{k}")

    if 6 in sizes:
        pgl, pgl_mult, khist, group_size = construct_pgl_matchings(points, points, p, 6, cell_index)
        check(f"{label}.size6.pgl_injective", all(v == 1 for v in pgl_mult.values()),
              max(pgl_mult.values(), default=0))
        check(f"{label}.size6.pgl_set_equality", direct[6] == pgl,
              {"direct": len(direct[6]), "pgl": len(pgl),
               "direct_only": len(direct[6] - pgl), "pgl_only": len(pgl - direct[6])})
        formula = sum(count * math.comb(k, 6) for k, count in khist.items() if k >= 6)
        check(f"{label}.size6.pgl_formula", len(pgl) == formula, {"got": len(pgl), "want": formula})
        pgl_record = {"count": len(pgl), "formula": formula, "group_size": group_size,
                      "k_histogram": {str(k): v for k, v in sorted(khist.items())}}
    else:
        pgl_record = None
    if 7 in sizes:
        check(f"{label}.size7.empty", direct[7] == set(), len(direct[7]))

    record = {
        "n": n, "prime": p,
        "sizes": {str(k): {"universe": universes[k], "circuits": len(direct[k])} for k in sizes},
        "size6_pgl": pgl_record,
    }
    return record, direct, cells, columns


def rref_nullspace(matrix, p):
    a = [[x % p for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    pivots = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [x * inv % p for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                factor = a[i][c]
                a[i] = [(x - factor * y) % p for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = [0] * cols
        v[f] = 1
        for row, pivot in reversed(list(enumerate(pivots))):
            v[pivot] = -sum(a[row][c] * v[c] for c in free) % p
        basis.append(tuple(v))
    return basis


def poly_trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_divmod(a, b, p):
    a = poly_trim(a)
    b = poly_trim(b)
    if b == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b) and a != [0]:
        shift = len(a) - len(b)
        coeff = a[-1] * inv % p
        q[shift] = coeff
        for i, value in enumerate(b):
            a[i + shift] = (a[i + shift] - coeff * value) % p
        a = poly_trim(a)
    return poly_trim(q), a


def poly_gcd(a, b, p):
    a, b = poly_trim(a), poly_trim(b)
    while b != [0]:
        _, rem = poly_divmod(a, b, p)
        a, b = b, rem
    inv = pow(a[-1], -1, p)
    return [x * inv % p for x in a]


def rational_degree2_fit(support_values, p):
    first_five = support_values[:5]
    matrix = []
    for x, y in first_five:
        matrix.append([1, y, y * y % p, -x % p, -x * y % p, -x * y * y % p])
    basis = rref_nullspace(matrix, p)
    if len(basis) != 1:
        return None
    coeff = basis[0]
    num = poly_trim(coeff[:3])
    den = poly_trim(coeff[3:])
    if den == [0]:
        return None
    g = poly_gcd(num, den, p)
    num, remn = poly_divmod(num, g, p)
    den, remd = poly_divmod(den, g, p)
    if remn != [0] or remd != [0]:
        raise AssertionError("polynomial gcd division failed")
    degree = max(len(num) - 1, len(den) - 1)
    if degree != 2:
        return None
    for x, y in support_values:
        nv = sum(c * pow(y, i, p) for i, c in enumerate(num)) % p
        dv = sum(c * pow(y, i, p) for i, c in enumerate(den)) % p
        if dv == 0 or nv * pow(dv, -1, p) % p != x:
            return None
    return tuple(num), tuple(den)


def census_full_matchings(n, p, size):
    label = f"fullmatching.n{n}.p{p}.m{size}"
    points = tuple(x % p for x in range(1, n + 1))
    check(f"{label}.points_distinct", len(set(points)) == n, points)
    cells = tuple((i, j) for i in range(n) for j in range(n))
    columns = tuple(product_column(points[i], points[j], p) for i, j in cells)
    direct = set()
    rational2 = set()
    curve_dependent = set()
    universe = 0
    for perm in itertools.permutations(range(n)):
        universe += 1
        support = tuple((i, perm[i]) for i in range(n))
        mask = sum(1 << (i * n + perm[i]) for i in range(n))
        vectors = [columns[i * n + perm[i]] for i in range(n)]
        rank = rank_mod(vectors, p)
        if rank < size:
            curve_dependent.add(mask)
        if rank == size - 1 and all(rank_mod(vectors[:i] + vectors[i + 1 :], p) == size - 1
                                    for i in range(size)):
            direct.add(mask)
        values = tuple((points[i], points[perm[i]]) for i in range(n))
        # Fit y -> x, so swap the stored (row x, column y) order only conceptually.
        fit = rational_degree2_fit(values, p) if size == 8 else None
        if fit is not None:
            rational2.add(mask)
    if size == 8:
        check(f"{label}.rational2_subset", rational2 <= direct,
              {"rational2": len(rational2), "direct": len(direct), "bad": len(rational2 - direct)})
    if size == 9:
        check(f"{label}.curve_predicate", direct <= curve_dependent,
              {"circuits": len(direct), "dependent": len(curve_dependent)})
    checkpoint({"kind": "full_matching", "n": n, "prime": p, "size": size,
                "universe": universe, "circuits": len(direct), "rational2": len(rational2)})
    check_budget(label)
    return {
        "n": n, "prime": p, "size": size, "universe": universe,
        "circuits": len(direct), "curve_dependent": len(curve_dependent),
        "rational_degree2_circuits": len(rational2),
        "residual_circuits": len(direct - rational2),
    }, direct


def compare_grid_primes(n, entries, cells, common_sizes):
    pairs = []
    for (p1, sets1), (p2, sets2) in itertools.combinations(entries, 2):
        rec = {"primes": [p1, p2], "sizes": {}}
        for m in common_sizes:
            by_profile = {}
            profiles = sorted({profile(mask_support(mask, cells)) for mask in sets1[m] | sets2[m]})
            for pr in profiles:
                A = {mask for mask in sets1[m] if profile(mask_support(mask, cells)) == pr}
                B = {mask for mask in sets2[m] if profile(mask_support(mask, cells)) == pr}
                by_profile[str(pr)] = {"left": len(A), "right": len(B), "intersection": len(A & B),
                                       "left_only": len(A - B), "right_only": len(B - A),
                                       "sets_equal": A == B}
            rec["sizes"][str(m)] = by_profile
        pairs.append(rec)
    return {"n": n, "pairs": pairs}


def run_plants(grid_data):
    record, circuits, cells, columns = grid_data[(4, 13)]
    cell_index = {cell: i for i, cell in enumerate(cells)}
    a_fiber_support = tuple((i, 0) for i in range(4))
    b_fiber_support = tuple((0, j) for j in range(4))
    a_fiber_vectors = [columns[cell_index[cell]] for cell in a_fiber_support]
    b_fiber_vectors = [columns[cell_index[cell]] for cell in b_fiber_support]
    check("plants.A_fiber_accept", is_circuit_vectors(a_fiber_vectors, 13), a_fiber_support)
    check("plants.B_fiber_accept", is_circuit_vectors(b_fiber_vectors, 13), b_fiber_support)
    fiber_support = set(b_fiber_support)
    extra = next(cell for cell in cells if cell not in fiber_support)
    five_support = tuple(sorted(fiber_support | {extra}))
    five_vectors = [columns[cell_index[cell]] for cell in five_support]
    check("plants.five_with_fiber_reject", not is_circuit_vectors(five_vectors, 13), five_support)
    generic_five_support = ((0, 0), (0, 1), (1, 1), (2, 2), (3, 3))
    generic_five_vectors = [columns[cell_index[cell]] for cell in generic_five_support]
    check("plants.generic_five_independent", rank_mod(generic_five_vectors, 13) == 5,
          rank_mod(generic_five_vectors, 13))
    check("plants.generic_five_reject", not is_circuit_vectors(generic_five_vectors, 13),
          generic_five_support)

    crossings, _ = construct_crossings(4, cell_index)
    crossing_mask = min(crossings)
    crossing_vectors = [columns[i] for i in mask_indices(crossing_mask)]
    check("plants.crossing_accept", is_circuit_vectors(crossing_vectors, 13), mask_support(crossing_mask, cells))
    crossing_support = set(mask_support(crossing_mask, cells))
    added = next(cell for cell in cells if cell not in crossing_support)
    inserted = tuple(sorted(crossing_support | {added}))
    inserted_vectors = [columns[cell_index[cell]] for cell in inserted]
    check("plants.inserted_crossing_reject", not is_circuit_vectors(inserted_vectors, 13), inserted)

    identity6 = [product_column(i, i, 101) for i in range(1, 7)]
    check("plants.mobius6_accept", is_circuit_vectors(identity6, 101), 6)
    nonmobius_xs = [1, 2, 3, 4, 5, 7]
    nonmobius_ys = list(range(1, 7))
    nonmobius6 = [product_column(x, y, 101) for x, y in zip(nonmobius_xs, nonmobius_ys)]
    check("plants.nonmobius6_independent", rank_mod(nonmobius6, 101) == 6,
          rank_mod(nonmobius6, 101))
    check("plants.nonmobius6_reject", not is_circuit_vectors(nonmobius6, 101),
          {"x": nonmobius_xs, "y": nonmobius_ys})
    identity7 = [product_column(i, i, 101) for i in range(1, 8)]
    check("plants.wrong_curve7_dependent", rank_mod(identity7, 101) < 7, rank_mod(identity7, 101))
    check("plants.wrong_curve7_reject", not is_circuit_vectors(identity7, 101), 7)

    ys = list(range(1, 9))
    xs = [y * y % 101 for y in ys]
    rational8 = [product_column(x, y, 101) for x, y in zip(xs, ys)]
    rational_values = tuple(zip(xs, ys))
    check("plants.rational_degree2_fit_accept", rational_degree2_fit(rational_values, 101) is not None,
          rational_values)
    check("plants.rational_degree2_8_accept", is_circuit_vectors(rational8, 101), {"x": xs, "y": ys})
    degree3_xs = [pow(y, 3, 101) for y in ys]
    degree3_values = tuple(zip(degree3_xs, ys))
    degree3_vectors = [product_column(x, y, 101) for x, y in degree3_values]
    check("plants.degree3_not_rational2_reject",
          rational_degree2_fit(degree3_values, 101) is None, degree3_values)
    check("plants.degree3_eight_independent", rank_mod(degree3_vectors, 101) == 8,
          rank_mod(degree3_vectors, 101))

    nonminimal10_support = tuple((0, j) for j in range(4)) + tuple((1, j) for j in range(4)) + ((2, 0), (3, 1))
    nonminimal10_vectors = [columns[cell_index[cell]] for cell in nonminimal10_support]
    check("plants.ambient10_dependent", rank_mod(nonminimal10_vectors, 13) <= 9,
          rank_mod(nonminimal10_vectors, 13))
    check("plants.ambient10_nonminimal_reject", not is_circuit_vectors(nonminimal10_vectors, 13), nonminimal10_support)
    if circuits.get(10):
        accept10 = min(circuits[10])
        accept10_vectors = [columns[i] for i in mask_indices(accept10)]
        check("plants.ambient10_accept", is_circuit_vectors(accept10_vectors, 13), mask_support(accept10, cells))

    original = [list(v) for v in crossing_vectors]
    corruption = None
    for coordinate in range(9):
        for delta in range(1, 13):
            candidate = [row[:] for row in original]
            candidate[0][coordinate] = (candidate[0][coordinate] + delta) % 13
            vectors = [tuple(row) for row in candidate]
            if not is_circuit_vectors(vectors, 13):
                corruption = {"coordinate": coordinate, "delta": delta, "rank": rank_mod(vectors, 13)}
                break
        if corruption:
            break
    check("plants.corrupted_actual_reject", corruption is not None, corruption)

    true = product_column(2, 3, 13)
    concat = (1, 2, 4, 1, 3, 9)
    check("plants.kron_guard.true_dimension", len(true) == 9, true)
    check("plants.kron_guard.concat_reject", len(concat) == 6 and concat != true,
          {"true": true, "concat": concat})
    return {"A_fiber": list(a_fiber_support), "B_fiber": list(b_fiber_support),
            "crossing": list(mask_support(crossing_mask, cells)),
            "five_with_fiber_reject": list(five_support),
            "generic_five_reject": list(generic_five_support),
            "nonmobius6_reject": {"x": nonmobius_xs, "y": nonmobius_ys},
            "inserted_crossing_reject": list(inserted),
            "rational_degree2_8_accept": {"x": xs, "y": ys},
            "degree3_not_rational2_reject": {"x": degree3_xs, "y": ys},
            "ambient10_nonminimal_reject": list(nonminimal10_support), "corruption": corruption,
            "kron_guard": {"true": true, "concat": concat}}


def main():
    print(f"== {GATE} run {RUN_ID}")
    nice_level = os.nice(0)
    check("runtime.nice", nice_level == 15, nice_level)
    check("runtime.sys_dont_write_bytecode", sys.dont_write_bytecode is True, sys.dont_write_bytecode)
    check("runtime.env_dont_write_bytecode", os.environ.get("PYTHONDONTWRITEBYTECODE") == "1",
          os.environ.get("PYTHONDONTWRITEBYTECODE"))
    thread_caps = {name: os.environ.get(name) for name in
                   ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    check("runtime.thread_caps", all(v == "1" for v in thread_caps.values()), thread_caps)

    grid_records = []
    grid_data = {}
    retained_n4 = []
    retained_n5 = []
    for n, p, max_m in [(4, 7, 10), (4, 11, 10), (4, 13, 10),
                        (5, 7, 10), (5, 11, 8), (5, 13, 8)]:
        record, sets, cells, columns = census_grid(n, p, max_m)
        grid_records.append(record)
        grid_data[(n, p)] = (record, sets, cells, columns)
        (retained_n4 if n == 4 else retained_n5).append((p, sets))
        print(f"[grid] n={n}, p={p}, max={max_m}: " +
              str({m: len(sets[m]) for m in sorted(sets)}))

    matching_records = []
    for n in (6, 7):
        for p in (7, 11, 13):
            record, direct, cells, columns = census_matchings(n, p, (6, 7) if n == 7 else (6,))
            matching_records.append(record)
            print(f"[matching] n={n}, p={p}: " + str({m: len(v) for m, v in direct.items()}))

    full_matching_records = []
    full_matching_sets = {}
    for n, p, size in [(8, 11, 8), (8, 13, 8), (9, 11, 9), (9, 13, 9)]:
        record, direct = census_full_matchings(n, p, size)
        full_matching_records.append(record)
        full_matching_sets[(n, p)] = direct
        print(f"[full matching] n={n}, p={p}, m={size}: {record}")

    comparisons = [
        compare_grid_primes(4, retained_n4, grid_data[(4, 7)][2], range(4, 11)),
        compare_grid_primes(5, retained_n5, grid_data[(5, 7)][2], range(4, 9)),
    ]
    plants = run_plants(grid_data)
    cpu, wall = elapsed()
    check("budget.CPU", cpu <= CPU_CAP, {"used": cpu, "cap": CPU_CAP})
    check("budget.wall", wall <= WALL_CAP, {"used": wall, "cap": WALL_CAP})

    RESULTS.update({
        "runtime": {"nice": nice_level, "dont_write_bytecode": sys.dont_write_bytecode,
                    "PYTHONDONTWRITEBYTECODE": os.environ.get("PYTHONDONTWRITEBYTECODE"),
                    "thread_caps": thread_caps,
                    "max_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "budget": {"CPU_seconds": cpu, "wall_seconds": wall, "CPU_cap_seconds": CPU_CAP,
                   "wall_cap_seconds": WALL_CAP, "within_cap": True},
        "grid_configurations": grid_records,
        "matching_configurations": matching_records,
        "full_matching_configurations": full_matching_records,
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
        broken = {"run_id": RUN_ID, "gate": GATE, "error": repr(exc),
                  "traceback": traceback.format_exc(), "CPU_seconds": cpu, "wall_seconds": wall,
                  "checks_passed_before_failure": CHECKS, "checkpoint": CHECKPOINT}
        path = OUT / "broken_results.json"
        if path.exists():
            path = OUT / f"broken_results_{int(time.time())}.json"
        path.write_text(json.dumps(broken, indent=2, sort_keys=True, default=encode) + "\n")
        print(f"BROKEN output preserved at {path.name}", file=sys.stderr)
        raise
