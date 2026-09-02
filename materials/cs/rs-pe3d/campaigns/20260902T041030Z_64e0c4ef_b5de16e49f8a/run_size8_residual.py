#!/usr/bin/env python3
"""Exact instrument for H-33ROW-SIZE8-RESIDUAL.

Standalone: imports no frozen campaign code. Mathematical decisions use only
Python integers modulo p, two independent ranks, exact nullspaces, exact
Sylvester resultants, and exact deletion predicates.
"""
import sys
sys.dont_write_bytecode = True

import hashlib
import itertools
import json
import math
import os
import resource
import time
import traceback
from collections import Counter
from pathlib import Path

RUN_ID = "20260902T041030Z_64e0c4ef_b5de16e49f8a"
GATE = "H-33ROW-SIZE8-RESIDUAL"
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
    "arithmetic": "exact Python integers modulo p; true 3-by-3 row-major Kronecker columns",
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


def product_column(x, y, p):
    return tuple(pow(x, i, p) * pow(y, j, p) % p for i in range(3) for j in range(3))


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
    """Incremental basis with explicit pivot normalization; independent of rank_mod."""
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


def det_mod(matrix, p):
    a = [[value % p for value in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("determinant requires square matrix")
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


def rref_nullspace(matrix, p):
    a = [[value % p for value in row] for row in matrix]
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
        a[r] = [value * inv % p for value in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                factor = a[i][c]
                a[i] = [(u - factor * v) % p for u, v in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        vector = [0] * cols
        vector[f] = 1
        for row, pivot in reversed(list(enumerate(pivots))):
            vector[pivot] = -sum(a[row][c] * vector[c] for c in free) % p
        basis.append(tuple(vector))
    return basis


def is_circuit(columns, p):
    m = len(columns)
    rank = rank_mod(columns, p)
    return rank == m - 1 and all(rank_mod(columns[:i] + columns[i + 1 :], p) == m - 1
                                 for i in range(m))


def poly_trim(poly, p=None):
    out = list(poly)
    if p is not None:
        out = [value % p for value in out]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [0]


def poly_is_zero(poly):
    return poly_trim(poly) == [0]


def poly_degree(poly):
    poly = poly_trim(poly)
    return -1 if poly == [0] else len(poly) - 1


def poly_add(a, b, p):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return poly_trim(out)


def poly_scale(a, scalar, p):
    return poly_trim([scalar * value % p for value in a])


def poly_mul(a, b, p):
    if poly_is_zero(a) or poly_is_zero(b):
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i + j] = (out[i + j] + u * v) % p
    return poly_trim(out)


def poly_divmod(a, b, p):
    a = poly_trim(a, p)
    b = poly_trim(b, p)
    if b == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, p)
    while a != [0] and len(a) >= len(b):
        shift = len(a) - len(b)
        coeff = a[-1] * inv % p
        q[shift] = coeff
        for i, value in enumerate(b):
            a[i + shift] = (a[i + shift] - coeff * value) % p
        a = poly_trim(a)
    return poly_trim(q), a


def poly_gcd(a, b, p):
    a, b = poly_trim(a, p), poly_trim(b, p)
    while b != [0]:
        _, remainder = poly_divmod(a, b, p)
        a, b = b, remainder
    if a == [0]:
        return [0]
    return poly_scale(a, pow(a[-1], -1, p), p)


def poly_monic(poly, p):
    poly = poly_trim(poly, p)
    if poly == [0]:
        return [0]
    return poly_scale(poly, pow(poly[-1], -1, p), p)


def poly_det(matrix, p):
    n = len(matrix)
    total = [0]
    for perm in itertools.permutations(range(n)):
        inversions = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        term = [1]
        for i, j in enumerate(perm):
            term = poly_mul(term, matrix[i][j], p)
        if inversions % 2:
            term = poly_scale(term, -1, p)
        total = poly_add(total, term, p)
    return poly_trim(total)


def form_coefficients_in_x(form):
    return [poly_trim(form[3 * i : 3 * i + 3]) for i in range(3)]


def transpose_form(form):
    return tuple(form[3 * i + j] for j in range(3) for i in range(3))


def resultant_x(form_f, form_g, p):
    a0, a1, a2 = form_coefficients_in_x(form_f)
    b0, b1, b2 = form_coefficients_in_x(form_g)
    zero = [0]
    matrix = [
        [a2, a1, a0, zero],
        [zero, a2, a1, a0],
        [b2, b1, b0, zero],
        [zero, b2, b1, b0],
    ]
    return poly_det(matrix, p)


def homogeneous_content_positive(forms, p):
    coefficients = []
    for form in forms:
        coefficients.extend(poly for poly in form_coefficients_in_x(form) if not poly_is_zero(poly))
    if not coefficients:
        raise AssertionError("zero pencil basis")
    gcd = coefficients[0]
    for poly in coefficients[1:]:
        gcd = poly_gcd(gcd, poly, p)
    finite = poly_degree(gcd) > 0
    infinity_multiplicity = min(2 - poly_degree(poly) for poly in coefficients)
    return finite or infinity_multiplicity > 0, {
        "finite_content_degree": max(0, poly_degree(gcd)),
        "infinity_content_multiplicity": max(0, infinity_multiplicity),
    }


def common_component_info(form_f, form_g, p):
    res_x = resultant_x(form_f, form_g, p)
    y_content, y_detail = homogeneous_content_positive((form_f, form_g), p)
    trans_f, trans_g = transpose_form(form_f), transpose_form(form_g)
    res_y = resultant_x(trans_f, trans_g, p)
    x_content, x_detail = homogeneous_content_positive((trans_f, trans_g), p)
    route_x = poly_is_zero(res_x) or y_content
    route_y = poly_is_zero(res_y) or x_content
    if route_x != route_y:
        raise AssertionError({"common_component_routes_disagree": [route_x, route_y],
                              "res_x": res_x, "res_y": res_y,
                              "x_content": x_detail, "y_content": y_detail})
    return {
        "has_common_component": route_x,
        "resultant_x": res_x,
        "resultant_y": res_y,
        "x_route_resultant_zero": poly_is_zero(res_x),
        "y_route_resultant_zero": poly_is_zero(res_y),
        "x_content": x_detail,
        "y_content": y_detail,
    }


def node_polynomial(values, p):
    result = [1]
    for value in values:
        result = poly_mul(result, [-value % p, 1], p)
    return poly_monic(result, p)


def form_eval(form, x, y, p):
    return sum(form[3 * i + j] * pow(x, i, p) * pow(y, j, p)
               for i in range(3) for j in range(3)) % p


def form_mul(left, right, p):
    out = [0] * 9
    for i in range(3):
        for j in range(3):
            u = left[3 * i + j]
            if not u:
                continue
            for k in range(3 - i):
                for ell in range(3 - j):
                    v = right[3 * k + ell]
                    if v:
                        out[3 * (i + k) + j + ell] = (out[3 * (i + k) + j + ell] + u * v) % p
    return tuple(out)


def form_from_terms(terms, p):
    out = [0] * 9
    for (i, j), value in terms.items():
        out[3 * i + j] = value % p
    return tuple(out)


def pencil_basis(columns, p):
    basis = rref_nullspace([list(column) for column in columns], p)
    return basis


def maximal_row_minor_zero_mask(columns, p):
    if len(columns) != 8 or any(len(column) != 9 for column in columns):
        raise ValueError("expected 9-by-8 evaluation columns")
    mask = 0
    for omitted in range(9):
        rows = [r for r in range(9) if r != omitted]
        matrix = [[columns[c][r] for c in range(8)] for r in rows]
        if det_mod(matrix, p) == 0:
            mask |= 1 << omitted
    return mask


def rational_degree2_fit(values, p):
    matrix = []
    for x, y in values[:5]:
        matrix.append([1, y, y * y % p, -x % p, -x * y % p, -x * y * y % p])
    basis = rref_nullspace(matrix, p)
    if len(basis) != 1:
        return None
    coeff = basis[0]
    num = poly_trim(coeff[:3])
    den = poly_trim(coeff[3:])
    if den == [0]:
        return None
    gcd = poly_gcd(num, den, p)
    num, remainder_num = poly_divmod(num, gcd, p)
    den, remainder_den = poly_divmod(den, gcd, p)
    if remainder_num != [0] or remainder_den != [0]:
        raise AssertionError("polynomial gcd division failed")
    if max(poly_degree(num), poly_degree(den)) != 2:
        return None
    for x, y in values:
        numerator = sum(c * pow(y, i, p) for i, c in enumerate(num)) % p
        denominator = sum(c * pow(y, i, p) for i, c in enumerate(den)) % p
        if denominator == 0 or numerator * pow(denominator, -1, p) % p != x:
            return None
    return tuple(num), tuple(den)


def has_mobius_six(columns, p):
    return any(rank_mod([columns[i] for i in chosen], p) == 5
               for chosen in itertools.combinations(range(8), 6))


def set_digest(permutations):
    payload = json.dumps(sorted(permutations), separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def fingerprint(columns, common_info, graph12, graph21, p):
    rank = rank_mod(columns, p)
    deletion_ranks = tuple(rank_mod(columns[:i] + columns[i + 1 :], p) for i in range(8))
    minor_mask = maximal_row_minor_zero_mask(columns, p)
    return str((
        (8, 8),
        "8K2_isolated_perfect_matching",
        (1,) * 8,
        (1,) * 8,
        8 - rank,
        9 - rank,
        deletion_ranks,
        minor_mask,
        common_info["has_common_component"],
        graph12,
        graph21,
        poly_degree(common_info["resultant_x"]),
        poly_degree(common_info["resultant_y"]),
    ))


def census_prime(p):
    label = f"n8.p{p}"
    points = tuple(range(1, 9))
    check(f"{label}.points_distinct", len(set(value % p for value in points)) == 8, points)
    rank7_set = set()
    direct_set = set()
    lower6_set = set()
    graph12_set = set()
    graph21_set = set()
    ci_set = set()
    common_set = set()
    fingerprints = Counter()
    direct_minor_masks = Counter()
    rank_crosschecks = 0
    resultant_certificates = 0
    common_route_checks = 0

    want_x_nodes = node_polynomial(points, p)
    want_y_nodes = node_polynomial(points, p)
    for index, perm in enumerate(itertools.permutations(range(8))):
        values = tuple((points[i] % p, points[perm[i]] % p) for i in range(8))
        columns = [product_column(x, y, p) for x, y in values]
        r1 = rank_mod(columns, p)
        r2 = rank_mod_alt(columns, p)
        if r1 != r2:
            raise AssertionError(f"{label}.rank disagreement at {perm}: {r1}!={r2}")
        rank_crosschecks += 1
        fit12 = rational_degree2_fit(values, p)
        fit21 = rational_degree2_fit(tuple((y, x) for x, y in values), p)
        if fit12 is not None:
            graph12_set.add(perm)
        if fit21 is not None:
            graph21_set.add(perm)
        if r1 == 7:
            rank7_set.add(perm)
            deletion_ranks = tuple(rank_mod(columns[:i] + columns[i + 1 :], p) for i in range(8))
            direct = all(value == 7 for value in deletion_ranks)
            lower6 = has_mobius_six(columns, p)
            basis = pencil_basis(columns, p)
            if len(basis) != 2:
                raise AssertionError(f"{label}.rank7 nullity is not 2: {perm}, {len(basis)}")
            common = common_component_info(basis[0], basis[1], p)
            common_route_checks += 1
            if common["has_common_component"]:
                common_set.add(perm)
            else:
                ci_set.add(perm)
                got_y_nodes = poly_monic(common["resultant_x"], p)
                got_x_nodes = poly_monic(common["resultant_y"], p)
                if got_y_nodes != want_y_nodes or got_x_nodes != want_x_nodes:
                    raise AssertionError({"stage": label, "perm": perm,
                                          "got_y": got_y_nodes, "want_y": want_y_nodes,
                                          "got_x": got_x_nodes, "want_x": want_x_nodes})
                resultant_certificates += 1
            if direct:
                direct_set.add(perm)
                minor_mask = maximal_row_minor_zero_mask(columns, p)
                direct_minor_masks[minor_mask] += 1
                if minor_mask != (1 << 9) - 1 or deletion_ranks != (7,) * 8:
                    raise AssertionError({"stage": label, "perm": perm,
                                          "minor_mask": minor_mask,
                                          "deletion_ranks": deletion_ranks})
                fingerprints[fingerprint(columns, common, fit12 is not None, fit21 is not None, p)] += 1
            if lower6:
                lower6_set.add(perm)
        if index % 2048 == 0:
            check_budget(f"{label}.perm{index}")

    graph_union = graph12_set | graph21_set
    rank7_noncircuit = rank7_set - direct_set
    check(f"{label}.graph12_subset_direct", graph12_set <= direct_set,
          {"graph12": len(graph12_set), "bad": len(graph12_set - direct_set)})
    check(f"{label}.graph21_subset_direct", graph21_set <= direct_set,
          {"graph21": len(graph21_set), "bad": len(graph21_set - direct_set)})
    check(f"{label}.graph_orientations_disjoint", not (graph12_set & graph21_set),
          len(graph12_set & graph21_set))
    check(f"{label}.classification_set_equality", direct_set == graph_union | ci_set,
          {"direct": len(direct_set), "graphs": len(graph_union), "ci": len(ci_set),
           "missing": len(direct_set - (graph_union | ci_set)),
           "extra": len((graph_union | ci_set) - direct_set)})
    check(f"{label}.classification_disjoint", not (graph_union & ci_set), len(graph_union & ci_set))
    check(f"{label}.rank7_noncircuit_lower6", rank7_noncircuit == lower6_set,
          {"noncircuit": len(rank7_noncircuit), "lower6": len(lower6_set),
           "symmetric_difference": len(rank7_noncircuit ^ lower6_set)})
    check(f"{label}.rank7_noncircuit_common", rank7_noncircuit == common_set - graph_union,
          {"noncircuit": len(rank7_noncircuit), "common_non_graph": len(common_set - graph_union),
           "symmetric_difference": len(rank7_noncircuit ^ (common_set - graph_union))})
    check(f"{label}.ci_resultant_all", resultant_certificates == len(ci_set),
          {"certificates": resultant_certificates, "ci": len(ci_set)})
    check(f"{label}.fingerprint_population", sum(fingerprints.values()) == len(direct_set),
          {"fingerprinted": sum(fingerprints.values()), "direct": len(direct_set)})
    pin(len(direct_set), 560 if p == 11 else 416, f"{label}.direct_anchor")
    pin(len(graph12_set), 0, f"{label}.graph12_zero")
    pin(len(graph21_set), 0, f"{label}.graph21_zero")
    pin(len(ci_set), 560 if p == 11 else 416, f"{label}.ci_anchor")
    check(f"{label}.all_direct_minor_masks",
          direct_minor_masks == Counter({(1 << 9) - 1: len(direct_set)}),
          dict(direct_minor_masks))

    record = {
        "n": 8,
        "prime": p,
        "universe": math.factorial(8),
        "rank_crosschecks": rank_crosschecks,
        "rank7_pencils": len(rank7_set),
        "rank7_wrong_predicate_overacceptance": len(rank7_noncircuit),
        "direct_circuits": len(direct_set),
        "six_mobius_nonminimal": len(lower6_set),
        "graph_12": len(graph12_set),
        "graph_21": len(graph21_set),
        "complete_intersections": len(ci_set),
        "common_component_pencils": len(common_set),
        "resultant_certificates": resultant_certificates,
        "direct_8x8_row_minor_zero_masks": {str(key): value for key, value in sorted(direct_minor_masks.items())},
        "common_component_route_checks": common_route_checks,
        "set_digests": {
            "rank7": set_digest(rank7_set),
            "direct": set_digest(direct_set),
            "lower6": set_digest(lower6_set),
            "graph12": set_digest(graph12_set),
            "graph21": set_digest(graph21_set),
            "complete_intersection": set_digest(ci_set),
            "common_component": set_digest(common_set),
        },
        "fingerprint_definition": [
            "profile", "bipartite_template", "A_degrees", "B_degrees",
            "rank_drop_from_8", "left_kernel_dimension", "seven_deletion_ranks",
            "8x8_row_minor_zero_mask", "common_component", "graph12", "graph21",
            "resultant_x_degree", "resultant_y_degree"
        ],
        "fingerprint_histogram": dict(sorted(fingerprints.items())),
        "first_ci_permutation": list(min(ci_set)) if ci_set else None,
        "first_rank7_noncircuit": list(min(rank7_noncircuit)) if rank7_noncircuit else None,
    }
    checkpoint({"kind": "full_matching_pencil", "n": 8, "prime": p,
                "universe": math.factorial(8), "rank7": len(rank7_set),
                "direct": len(direct_set), "ci": len(ci_set),
                "nonminimal": len(rank7_noncircuit)})
    return record


def find_ci_plant(p):
    one = form_from_terms({(0, 0): 1}, p)
    del one
    roots = {}
    c_values = range(1, 21)
    s_values = range(0, 21)
    for c in c_values:
        for s in s_values:
            roots[c, s] = tuple((x, c * pow(x, -1, p) % p)
                                for x in range(1, p)
                                if (c * pow(x, -1, p) - x - s) % p == 0)
    for c1, c2 in itertools.combinations(c_values, 2):
        for s1, s2 in itertools.combinations(s_values, 2):
            points = set()
            valid = True
            for c in (c1, c2):
                for s in (s1, s2):
                    pair = roots[c, s]
                    if len(pair) != 2:
                        valid = False
                        break
                    points.update(pair)
                if not valid:
                    break
            if not valid or len(points) != 8:
                continue
            points = tuple(sorted(points))
            if len({x for x, _ in points}) != 8 or len({y for _, y in points}) != 8:
                continue
            xy_c1 = form_from_terms({(1, 1): 1, (0, 0): -c1}, p)
            xy_c2 = form_from_terms({(1, 1): 1, (0, 0): -c2}, p)
            diag_s1 = form_from_terms({(0, 1): 1, (1, 0): -1, (0, 0): -s1}, p)
            diag_s2 = form_from_terms({(0, 1): 1, (1, 0): -1, (0, 0): -s2}, p)
            form_f = form_mul(xy_c1, xy_c2, p)
            form_g = form_mul(diag_s1, diag_s2, p)
            scanned = tuple((x, y) for x in range(p) for y in range(p)
                            if form_eval(form_f, x, y, p) == 0 and form_eval(form_g, x, y, p) == 0)
            if tuple(sorted(scanned)) != points:
                continue
            columns = [product_column(x, y, p) for x, y in points]
            if not is_circuit(columns, p):
                continue
            common = common_component_info(form_f, form_g, p)
            if common["has_common_component"]:
                continue
            return {
                "parameters": {"c": [c1, c2], "s": [s1, s2]},
                "points": points,
                "form_f": form_f,
                "form_g": form_g,
                "columns": columns,
                "common": common,
            }
    raise AssertionError("deterministic CI plant search found no witness")


def run_plants():
    p = 101
    ci = find_ci_plant(p)
    ci_columns = ci["columns"]
    ci_basis = pencil_basis(ci_columns, p)
    pin(len(ci_basis), 2, "plants.ci.left_kernel_dimension")
    ci_common = common_component_info(ci_basis[0], ci_basis[1], p)
    check("plants.ci.no_common_component", not ci_common["has_common_component"], ci["parameters"])
    check("plants.ci.accept", is_circuit(ci_columns, p), ci["points"])
    pin(maximal_row_minor_zero_mask(ci_columns, p), (1 << 9) - 1, "plants.ci.all_8x8_minors_zero")
    pin(poly_monic(ci_common["resultant_x"], p),
        node_polynomial([y for _, y in ci["points"]], p), "plants.ci.resultant_y_nodes")
    pin(poly_monic(ci_common["resultant_y"], p),
        node_polynomial([x for x, _ in ci["points"]], p), "plants.ci.resultant_x_nodes")

    ys = list(range(1, 9))
    graph12_points = tuple((y * y % p, y) for y in ys)
    graph21_points = tuple((y, y * y % p) for y in ys)
    for name, points, fit_values in (
        ("graph12", graph12_points, graph12_points),
        ("graph21", graph21_points, tuple((y, x) for x, y in graph21_points)),
    ):
        columns = [product_column(x, y, p) for x, y in points]
        check(f"plants.{name}.coordinates_distinct",
              len({x for x, _ in points}) == 8 and len({y for _, y in points}) == 8, points)
        check(f"plants.{name}.degree2_fit", rational_degree2_fit(fit_values, p) is not None, points)
        check(f"plants.{name}.accept", is_circuit(columns, p), points)
        basis = pencil_basis(columns, p)
        pin(len(basis), 2, f"plants.{name}.left_kernel_dimension")
        check(f"plants.{name}.common_component",
              common_component_info(basis[0], basis[1], p)["has_common_component"], points)

    degree3_points = tuple((pow(y, 3, p), y) for y in ys)
    degree3_columns = [product_column(x, y, p) for x, y in degree3_points]
    pin(rank_mod(degree3_columns, p), 8, "plants.wrong_degree.rank8")
    check("plants.wrong_degree.no_degree2_fit", rational_degree2_fit(degree3_points, p) is None,
          degree3_points)
    check("plants.wrong_degree.reject", not is_circuit(degree3_columns, p), degree3_points)

    nonminimal_points = tuple((i, i) for i in range(1, 7)) + ((7, 8), (8, 7))
    nonminimal_columns = [product_column(x, y, p) for x, y in nonminimal_points]
    pin(rank_mod(nonminimal_columns, p), 7, "plants.rank7_wrong.rank")
    pin(maximal_row_minor_zero_mask(nonminimal_columns, p), (1 << 9) - 1,
        "plants.rank7_wrong.all_8x8_minors_zero")
    check("plants.rank7_wrong.lower6", has_mobius_six(nonminimal_columns, p), nonminimal_points)
    check("plants.rank7_wrong.reject", not is_circuit(nonminimal_columns, p), nonminimal_points)
    nonminimal_basis = pencil_basis(nonminimal_columns, p)
    check("plants.rank7_wrong.common_component",
          common_component_info(nonminimal_basis[0], nonminimal_basis[1], p)["has_common_component"],
          nonminimal_points)

    # Common-component detector unit controls.
    h11 = form_from_terms({(1, 0): 1, (0, 1): -1}, p)
    lx = form_from_terms({(0, 0): 1, (1, 0): 1}, p)
    ly = form_from_terms({(0, 0): 1, (0, 1): 1}, p)
    check("plants.detector.common11", common_component_info(form_mul(h11, lx, p),
                                                             form_mul(h11, ly, p), p)["has_common_component"], None)
    h12 = form_from_terms({(1, 0): 1, (0, 2): -1}, p)
    check("plants.detector.common12", common_component_info(h12, form_mul(h12, lx, p), p)["has_common_component"], None)
    hx = form_from_terms({(1, 0): 1, (0, 0): -3}, p)
    check("plants.detector.x_only", common_component_info(form_mul(hx, ly, p),
                                                           form_mul(hx, h11, p), p)["has_common_component"], None)
    hy = form_from_terms({(0, 1): 1, (0, 0): -4}, p)
    check("plants.detector.y_only", common_component_info(form_mul(hy, lx, p),
                                                           form_mul(hy, h11, p), p)["has_common_component"], None)
    check("plants.detector.coprime", not common_component_info(ci["form_f"], ci["form_g"], p)["has_common_component"],
          ci["parameters"])

    # Ambient ten-set tests over the 4-by-4 GF(13) grid.
    q = 13
    grid_points = tuple((x, y) for x in range(1, 5) for y in range(1, 5))
    grid_columns = tuple(product_column(x, y, q) for x, y in grid_points)
    accept10 = None
    for chosen in itertools.combinations(range(16), 10):
        columns = [grid_columns[i] for i in chosen]
        if is_circuit(columns, q):
            accept10 = chosen
            break
    check("plants.ambient10.accept_exists", accept10 is not None, accept10)
    check("plants.ambient10.accept", is_circuit([grid_columns[i] for i in accept10], q), accept10)
    nonminimal10_points = tuple((1, y) for y in range(1, 5)) + tuple((2, y) for y in range(1, 5)) + ((3, 1), (4, 2))
    nonminimal10_columns = [product_column(x, y, q) for x, y in nonminimal10_points]
    check("plants.ambient10.dependent", rank_mod(nonminimal10_columns, q) <= 9,
          rank_mod(nonminimal10_columns, q))
    check("plants.ambient10.reject", not is_circuit(nonminimal10_columns, q), nonminimal10_points)

    corruption = None
    for coordinate in range(9):
        for delta in range(1, p):
            changed = [list(column) for column in ci_columns]
            changed[0][coordinate] = (changed[0][coordinate] + delta) % p
            vectors = [tuple(column) for column in changed]
            if not is_circuit(vectors, p):
                corruption = {"coordinate": coordinate, "delta": delta, "rank": rank_mod(vectors, p)}
                break
        if corruption:
            break
    check("plants.corrupted_ci_reject", corruption is not None, corruption)

    true = product_column(2, 3, 13)
    concat = (1, 2, 4, 1, 3, 9)
    check("plants.kron_guard.true_dimension", len(true) == 9, true)
    check("plants.kron_guard.concat_reject", len(concat) == 6 and concat != true,
          {"true": true, "concat": concat})

    return {
        "ci_accept": {"parameters": ci["parameters"], "points": ci["points"]},
        "graph12_accept": graph12_points,
        "graph21_accept": graph21_points,
        "wrong_degree_reject": degree3_points,
        "rank7_wrong_reject": nonminimal_points,
        "ambient10_accept_indices": accept10,
        "ambient10_reject": nonminimal10_points,
        "corruption": corruption,
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

    plants = run_plants()
    records = []
    for p in (11, 13):
        record = census_prime(p)
        records.append(record)
        print(f"[n=8,p={p}] rank7={record['rank7_pencils']} direct={record['direct_circuits']} "
              f"CI={record['complete_intersections']} nonminimal={record['rank7_wrong_predicate_overacceptance']}")

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
        "budget": {
            "CPU_seconds": cpu,
            "wall_seconds": wall,
            "CPU_cap_seconds": CPU_CAP,
            "wall_cap_seconds": WALL_CAP,
            "within_cap": True,
        },
        "configurations": records,
        "plants": plants,
        "controls": CHECKS,
        "assertion_count": len(CHECKS),
    })
    (OUT / "controls_results.json").write_text(json.dumps(RESULTS, indent=2, sort_keys=True) + "\n")
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
