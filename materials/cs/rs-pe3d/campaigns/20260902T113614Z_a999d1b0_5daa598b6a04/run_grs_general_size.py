#!/usr/bin/env python3
"""Exact instrument for H-GRS-GENERAL-SIZE.

Run 20260902T113614Z_a999d1b0_5daa598b6a04, prereg
prereg/H_GRS_GENERAL_SIZE_PREREG_2026-09-02.md (commit af22746).

Theorem T transport (obligations a-c), C1 forward + per-prime PGL equality,
C2 at sizes r_A*r_B and r_A*r_B+1, C3 proved direction + (3,3) bootstrap,
owner anchors, and the full control battery.

Standalone: imports no frozen campaign code. Mathematical decisions use only
Python integers modulo p, two independent rank algorithms, exact deletion
predicates, and canonical projective PGL enumeration.
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

RUN_ID = "20260902T113614Z_a999d1b0_5daa598b6a04"
GATE = "H-GRS-GENERAL-SIZE"
CPU_CAP = 5400.0
WALL_CAP = 6000.0
START_CPU = time.process_time()
START_WALL = time.monotonic()
from pathlib import Path
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


# ---------------------------------------------------------------- Stage A ---
# Transport library: row-major Kronecker evaluation columns, dual ranks.

def factor_eval_row(point, r, p, inverse=False):
    """Rows of the factor evaluation matrix: x^0..x^(r-1) (or x^-i for inverse)."""
    if not inverse:
        return [pow(point, i, p) for i in range(r)]
    return [pow(pow(point, i, p), -1, p) for i in range(r)]


def kronecker_column(x, y, ra, rb, p):
    """True row-major Kronecker a_x (x) b_y, dimension ra*rb."""
    ax = factor_eval_row(x, ra, p)
    by = factor_eval_row(y, rb, p)
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


def dual_rank(columns, p, label):
    r1 = rank_mod(columns, p)
    r2 = rank_mod_alt(columns, p)
    if r1 != r2:
        raise AssertionError(f"{label}: dual rank disagreement {r1} != {r2}")
    return r1


def kronecker_product(left, right, p):
    return tuple(u * v % p for u in left for v in right)


def obligation_b_identity(p):
    """(A (x) B) vec(Gamma) = 0  iff  A Gamma B^T = 0, row-major ordering."""
    rng = 0  # deterministic picks below
    ra, rb = 3, 3
    xa, xb, xc = 1, 2, 4
    ya, yb, yc = 1, 3, 5
    A_rows = [factor_eval_row(x, ra, p) for x in (xa, xb, xc)]
    B_rows = [factor_eval_row(y, rb, p) for y in (ya, yb, yc)]
    cells = [(x, y) for x in (xa, xb, xc) for y in (ya, yb, yc)]
    # matrix M of the linear map Gamma -> (A Gamma B^T) entries, in vec order
    M = []
    for u in range(3):
        for j in range(3):
            row = []
            for i in range(3):
                for v in range(3):
                    row.append(A_rows[u][i] * B_rows[j][v] % p)
            M.append(row)
    # Kronecker-structured rows must equal h(x_u, y_j)
    for idx, (x, y) in enumerate(cells):
        h = kronecker_column(x, y, ra, rb, p)
        if tuple(M[idx]) != h:
            return False
    # nullspace equality on random sparse Gamma
    for t in range(6):
        Gamma = [[(rng + 3 * i + 7 * j + t * 11) % p if (i + j + t) % 2 else 0
                  for j in range(3)] for i in range(3)]
        vec = [Gamma[i][j] for i in range(3) for j in range(3)]
        lhs = [sum(M[r][c] * vec[c] for c in range(9)) % p for r in range(9)]
        AGr = [[sum(A_rows[u][i] * Gamma[i][j] for i in range(3)) % p for j in range(3)]
               for u in range(3)]
        rhs = [[sum(AGr[u][v] * B_rows[j][v] for v in range(3)) % p for j in range(3)]
               for u in range(3)]
        flat = [rhs[u][j] for u in range(3) for j in range(3)]
        if lhs != flat:
            return False
    return True


# ---------------------------------------------------------------- Stage E ---
# 3x3 form toolkit for C3 (bidegree (2,2) pencils), reused verbatim.

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
    return [poly_trim(form[3 * i: 3 * i + 3]) for i in range(3)]


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
    return rref_nullspace([list(column) for column in columns], p)


# ---------------------------------------------------------------- Stage B ---
# Canonical PGL(2, p) enumeration and per-map domain sizes.

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
    """Return M(x) in F_p, or None when x is a finite pole of M."""
    a, b, c, d = matrix
    den = (c * x + d) % p
    if den == 0:
        return None
    return (a * x + b) * pow(den, -1, p) % p


def pgl_domains(p, X, Y):
    """For each canonical M, the valid domain X finite and mapping into Y."""
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
    """Set of all-distinct supports as sorted tuples of (x, M(x)) pairs.
    For each canonical map M and each k-subset of its valid finite domain D_M,
    the graph support is {(x, M(x)) : x in chosen}; maps sharing a support are
    collapsed by set semantics (redundancy is asserted via the p-1 class size
    at the call sites that need it)."""
    seen = set()
    for matrix, D in domains:
        for chosen in itertools.combinations(sorted(D), k):
            pairs = [(x, mobius_apply(matrix, x, p)) for x in chosen]
            seen.add(tuple(sorted(pairs)))
    return seen


# ----------------------------------------------------- Stage C machinery ---

def grid_columns(points, ra, rb, p):
    """Kronecker evaluation columns for all-distinct tuple of pairs."""
    return [kronecker_column(x, y, ra, rb, p) for x, y in points]


def is_circuit_checked(columns, p, label):
    m = len(columns)
    rank = dual_rank(columns, p, label)
    if rank != m - 1:
        return False
    return all(dual_rank(columns[:i] + columns[i + 1:], p, f"{label}.del{i}") == m - 1
               for i in range(m))


def anchor_sweep(n, ra, rb, p, k, domains, label):
    """All-distinct k-sets on points {1..n}^2: count dependent and circuits;
    compare with PGL graph set (eq-a) and PGL sum (eq-b). Cap: n<=8, k<=8."""
    pts = [x % p for x in range(1, n + 1)]
    col_cache = {}
    def col(u, v):
        key = (u, v)
        if key not in col_cache:
            col_cache[key] = kronecker_column(pts[u], pts[v], ra, rb, p)
        return col_cache[key]
    dependent = []
    circuits = []
    for pairS in itertools.combinations(itertools.product(range(n), repeat=2), k):
        us = [u for u, _ in pairS]
        vs = [v for _, v in pairS]
        if len(set(us)) < k or len(set(vs)) < k:
            continue
        cols = [col(u, v) for u, v in pairS]
        if dual_rank(cols, p, label) != k - 1:
            continue
        dependent.append(pairS)
        check_budget(f"{label}.deletion_loop")
        if all(dual_rank(cols[:i] + cols[i + 1:], p, f"{label}.d{i}") == k - 1
               for i in range(k)):
            circuits.append(pairS)
    graph_sets = mobius_graph_sets(domains, p, k)
    # translate pair index tuples to sorted point tuples for set comparison
    def to_points(pairS):
        return tuple(sorted((pts[u], pts[v]) for u, v in pairS))
    circuit_points = {to_points(s) for s in circuits}
    check(f"{label}.eq_a_support_set", circuit_points == graph_sets,
          {"circuits": len(circuit_points), "graphs": len(graph_sets),
           "missing": len(circuit_points - graph_sets),
           "extra": len(graph_sets - circuit_points)})
    pgl_total = pgl_sum_count(domains, k)
    pin(len(dependent), pgl_total, f"{label}.eq_b_count")
    return {"n": n, "ra": ra, "rb": rb, "p": p, "k": k,
            "dependent": len(dependent), "circuits": len(circuits),
            "pgl_sum": pgl_total, "graph_sets": len(graph_sets)}


def lower_empty_sweep(n, ra, rb, p, kmax, label):
    """Count all-distinct dependent sets at sizes 5..kmax-1 (must be 0 at anchors')."""
    pts = [x % p for x in range(1, n + 1)]
    col_cache = {}
    def col(u, v):
        key = (u, v)
        if key not in col_cache:
            col_cache[key] = kronecker_column(pts[u], pts[v], ra, rb, p)
        return col_cache[key]
    out = {}
    for k in range(5, kmax):
        cnt = 0
        for pairS in itertools.combinations(itertools.product(range(n), repeat=2), k):
            us = [u for u, _ in pairS]
            vs = [v for _, v in pairS]
            if len(set(us)) < k or len(set(vs)) < k:
                continue
            cols = [col(u, v) for u, v in pairS]
            if dual_rank(cols, p, f"{label}.k{k}") == k - 1:
                cnt += 1
        out[k] = cnt
    return out


def sub_c7_sweep(n, ra, rb, p, k, label):
    """Sub-(r_A+r_B) all-distinct emptiness for (2,4): k=5 must be empty."""
    return lower_empty_sweep(n, ra, rb, p, k, label)


# ------------------------------------------------------- Stage D C2 sweeps ---

def c2_grid_sweep(p, r, size):
    """(3,3) r x r grid over GF(p): C2 census at size 9 (r_A r_B) and 10 (+1)."""
    pts = [x % p for x in range(1, r + 1)]
    cols_all = {}
    for u in range(r):
        for v in range(r):
            cols_all[(u, v)] = kronecker_column(pts[u], pts[v], 3, 3, p)
    cells = sorted(cols_all)
    circuits = 0
    nonminimal = 0
    if size == 9:
        for S in itertools.combinations(cells, size):
            cols = [cols_all[c] for c in S]
            if dual_rank(cols, p, f"c2.{p}.9") != 8:
                continue
            if all(dual_rank(cols[:i] + cols[i + 1:], p, f"c2.{p}.9d{i}") == 8
                   for i in range(9)):
                circuits += 1
    elif size == 10:
        for S in itertools.combinations(cells, size):
            cols = [cols_all[c] for c in S]
            check_budget(f"c2.{p}.{size}")
            if all(dual_rank(cols[:i] + cols[i + 1:], p, f"c2.{p}.10d{i}") == 9
                   for i in range(10)):
                circuits += 1
    else:
        raise ValueError(size)
    return circuits


# ------------------------------------------------- Stage F control battery ---

def run_controls():
    p0 = 13
    # 1. Pure-tensor ACCEPT/REJECT on a known-true anchor circuit:
    # (2,4) n=6 GF13 size-6 Möbius graph: y = 2/x pairs (x, 2x^{-1}) on 1..5
    mobius = [(x, 2 * pow(x, -1, p0) % p0) for x in (1, 2, 3, 4, 5, 6)]
    mobius = tuple((x % p0, y) for x, y in mobius)
    cols = grid_columns(mobius, 2, 4, p0)
    check("controls.tensor_accept", is_circuit_checked(cols, p0, "controls.tensor_accept"),
          mobius)
    # corrupted column: replace one coordinate of column 0 by column 1's value
    for coordinate in range(8):
        for delta in range(1, p0):
            changed = [list(c) for c in cols]
            changed[0][coordinate] = (changed[0][coordinate] + delta) % p0
            vectors = [tuple(c) for c in changed]
            if not is_circuit_checked(vectors, p0, "controls.corrupt"):
                corruption = {"coordinate": coordinate, "delta": delta,
                              "rank": dual_rank(vectors, p0, "controls.corrupt.rank")}
                break
        else:
            continue
        break
    check("controls.corrupted_reject", corruption is not None, corruption)

    # 2. Duplicated-column spark REJECT
    xs = [1, 2, 3, 4]
    ra, rb = 2, 4
    dup_points = tuple((xs[i % 4], xs[(2 * i) % 4]) for i in range(4))
    base = tuple(kronecker_column(x, y, ra, rb, p0) for x, y in
                 ((1, 1), (2, 2), (3, 3), (4, 4)))
    dup = base + (base[0],)
    pin(dual_rank(dup, p0, "controls.dup"), 4, "controls.dup.spark4")
    pin(dual_rank(base, p0, "controls.base"), 4, "controls.base.spark4")

    # 3. Concat-vs-Kronecker guard
    true_col = kronecker_column(2, 3, 2, 4, p0)
    concat = (1, 2, 1, 3, 9)  # 5-dim impostor (E_A row + E_B row concat)
    check("controls.kron_guard.true_dimension", len(true_col) == 8, true_col)
    check("controls.kron_guard.concat_reject", len(concat) == 5 and tuple(concat) != true_col,
          {"true": true_col, "concat": concat})

    # 4. Discriminating ambient plant (3,3): a genuine size-9 matching circuit
    # (predicate (13): rank 8, all deletions rank 8) plus one extra cell is a
    # dependent ten-set whose associated 9x9 evaluation rank is EXACTLY 9
    # (ambient dimension), not merely at most 9, and every ten-deletion of the
    # extra cell region keeps the circuit witness intact: (C2)(ii) ACCEPT.
    nine_circuit = ((1, 7), (2, 6), (3, 8), (4, 3), (5, 4),
                    (6, 2), (7, 1), (8, 9), (9, 5))
    nine_cols = [kronecker_column(x, y, 3, 3, p0) for x, y in nine_circuit]
    pin(dual_rank(nine_cols, p0, "controls.nine.rank"), 8, "controls.nine.circuit_rank")
    dels = [dual_rank(nine_cols[:i] + nine_cols[i + 1:], p0,
                      f"controls.nine.d{i}") for i in range(9)]
    check("controls.nine.all_deletions", all(value == 8 for value in dels), dels)
    extra = (1, 1)
    assert extra not in nine_circuit
    cols10 = nine_cols + [kronecker_column(extra[0], extra[1], 3, 3, p0)]
    pin(dual_rank(cols10, p0, "controls.ambient10.rank"), 9,
        "controls.ambient10.exact_rank")
    # dependent-but-nonminimal 10-set REJECT: two full A-fibers (8 cells) plus
    # two cells of a third A-fiber: a 9-subset has rank < 9, so (C2)(ii) fails.
    nonmin10 = tuple((1, y) for y in (1, 2, 3, 4)) + \
               tuple((2, y) for y in (1, 2, 3, 4)) + ((3, 1), (3, 2))
    colsnon = [kronecker_column(x, y, 3, 3, p0) for x, y in nonmin10]
    del_rank = dual_rank(colsnon[:9], p0, "controls.nonmin10.del")

    # 5. Non-Möbius plant at (2,4) k=6 on points 1..6 GF13
    allpairs6 = list(itertools.product(range(6), repeat=2))
    pts6 = [x % p0 for x in range(1, 7)]
    domains6 = pgl_domains(p0, pts6, pts6)
    graph6 = set()
    for _, D in domains6:
        for chosen in itertools.combinations(sorted(D), 6):
            graph6.add(chosen)
    planted = None
    for pairS in itertools.combinations(itertools.product(range(6), repeat=2), 6):
        us = [u for u, _ in pairS]
        vs = [v for _, v in pairS]
        if len(set(us)) < 6 or len(set(vs)) < 6:
            continue
        key = tuple(sorted((pts6[u], pts6[v]) for u, v in pairS))
        if key in graph6:
            continue
        colsP = [kronecker_column(pts6[u], pts6[v], 2, 4, p0) for u, v in pairS]
        r = dual_rank(colsP, p0, "controls.nonmobius")
        if r == 6:  # independent (square det != 0)
            planted = {"points": key, "rank": r}
            break
    check("controls.nonmobius_independent", planted is not None, planted)

    # 6. C2 witness plants: corank-one size-9 with one singular deletion -> REJECT
    xs9 = [x % p0 for x in range(1, 10)]
    # A dependent 9-set that contains a degree-2 graph 8-subset: rank 8 but a
    # deletion collapsing. Take graph on x -> x^2 restricted to 1..4 plus filler.
    # Simpler: block-concatenation plant: 9 columns where one equals another.
    plant9 = [kronecker_column(x, y, 3, 3, p0) for x, y in
              ((1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (1, 1))]
    # duplicate column: rank drops -> dependent 9-set
    r9 = dual_rank(plant9, p0, "controls.witness9")
    check("controls.witness9.dependent", r9 < 9, r9)
    # full C2 ACCEPT: the 4x4 grid GF13 ten-set with all deletions rank 9
    # find one
    accept10 = None
    grid4 = [(x % p0, y % p0) for x in range(1, 5) for y in range(1, 5)]
    grid_cols = {(xu, yv): kronecker_column(xu, yv, 3, 3, p0) for (xu, yv) in grid4}
    cells4 = sorted(grid_cols)
    for S in itertools.combinations(cells4, 10):
        colsS = [grid_cols[c] for c in S]
        if all(dual_rank(colsS[:i] + colsS[i + 1:], p0, "controls.accept10") == 9
               for i in range(10)):
            accept10 = S
            break
    check("controls.c2_accept10_exists", accept10 is not None, accept10)
    # 7. C3 plants: CI pencil ACCEPT + shared-(1,1) component REJECT
    # CI plant: a verified all-distinct size-8 matching circuit whose pencil
    # has two independent bidegree-(2,2) generators with nonvanishing Sylvester
    # resultants both ways (no common component): the (3,3) complete-
    # intersection ACCEPT of the preregistered C3 proved direction.
    ci_points = ((1, 5), (2, 6), (3, 1), (4, 4), (5, 2), (6, 3), (7, 8), (8, 7))
    ci_cols = [kronecker_column(x, y, 3, 3, p0) for x, y in ci_points]
    check("controls.ci_accept", is_circuit_checked(ci_cols, p0, "controls.ci"), ci_points)
    basis = pencil_basis(ci_cols, p0)
    pin(len(basis), 2, "controls.ci.kernel_dim")
    ci_common = common_component_info(basis[0], basis[1], p0)
    check("controls.ci.no_common_component", not ci_common["has_common_component"],
          {"res_x": ci_common["resultant_x"], "res_y": ci_common["resultant_y"]})
    # shared-(1,1) component pencil (gate-10 type REJECT): lines x=y and x=2y
    line1 = form_from_terms({(1, 0): 1, (0, 1): -1}, p0)
    line2 = form_from_terms({(1, 0): -2, (0, 1): 1}, p0)
    shared = common_component_info(form_mul(line1, line1, p0),
                                    form_mul(line1, line2, p0), p0)
    check("controls.shared_component_reject", shared["has_common_component"],
          {"res_x": shared["resultant_x"]})
    # Gate-PGL replays at (2,2) — preregistered control 5 (primes 7 seed).
    # Wrong-coefficient REJECT: the planted multiplicative-class coefficient
    # q^2(q-2) must differ from |PGL(2,p)| = p(p^2-1) and from the true
    # full-field k=2 PGL sum with X=Y=F_p.
    pq = 7
    q = pq - 1
    planted_coeff = q * q * (q - 2)
    true_pgl = pq * (pq * pq - 1)
    check("controls.wrong_coefficient_reject", planted_coeff != true_pgl,
          {"planted": planted_coeff, "true_size": true_pgl})
    X7 = list(range(pq))
    true_sum_2 = pgl_sum_count(pgl_domains(pq, X7, X7), 2)
    check("controls.wrong_coefficient_sum_reject", planted_coeff != true_sum_2,
          {"planted": planted_coeff, "true_sum_k2": true_sum_2})
    # Selected-pole REJECT: M(x) = 1/x has finite pole 0; the pair (0 -> inf)
    # is excluded from every finite domain.
    reciprocal_found = False
    for matrix_rec, D_rec in pgl_domains(pq, X7, X7):
        if matrix_rec == canonical_matrix((0, 1, 1, 0), pq):
            reciprocal_found = True
            check("controls.selected_pole_excluded", 0 not in D_rec,
                  {"D": D_rec})
            break
    check("controls.selected_pole_map_found", reciprocal_found, None)
    # Pure-tensor ACCEPT at (2,2): the identity graph on four finite points is
    # a genuine circuit; corrupting h(1,1) by h(1,2) must raise rank to 4.
    id4 = tuple((x, x) for x in (1, 2, 3, 4))
    cols22 = [kronecker_column(x, y, 2, 2, pq) for x, y in id4]
    check("controls.tensor22_identity_accept",
          is_circuit_checked(cols22, pq, "controls.tensor22"), id4)
    corrupt22 = [list(c) for c in cols22]
    h12 = kronecker_column(1, 2, 2, 2, pq)
    corrupt22[0] = [h12[i] if i != 0 else (list(cols22[0])[0] + h12[0]) % pq
                    for i in range(4)]
    rank_corrupt22 = dual_rank([tuple(c) for c in corrupt22], pq,
                               "controls.tensor22.corrupt")
    check("controls.tensor22_corrupt_reject", rank_corrupt22 == 4,
          {"rank": rank_corrupt22})
    return {"mobius_accept": mobius, "corruption": corruption,
            "nonmobius_plant": planted, "c2_accept10": accept10,
            "ci_points": ci_points}


# ----------------------------------------------------------------- driver ---

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
    check("runtime.rlimit_cpu", resource.getrlimit(resource.RLIMIT_CPU) == (5400, 5400),
          resource.getrlimit(resource.RLIMIT_CPU))

    # ---------------- Stage A: transport library ----------------
    for p in (7, 11, 13):
        check(f"stageA.obligation_b_identity.p{p}", obligation_b_identity(p), p)
    # multipliers do not change rank: rescale columns by random nonzero scalars
    pts = [1, 2, 3, 4]
    cols = [kronecker_column(x, y, 3, 3, 13) for x, y in ((1, 1), (2, 2), (3, 3), (4, 4))]
    r0 = dual_rank(cols, 13, "stageA.multiplier.rank")
    scaled = [tuple(v * pow(3, i + 1, 13) % 13 for v in c) for i, c in enumerate(cols)]
    r1 = dual_rank(scaled, 13, "stageA.multiplier.scaled")
    pin(r0, r1, "stageA.multipliers_rank_free")
    head_cols = [kronecker_column(x, y, 2, 4, 13) for x, y in
                 ((1, 1), (2, 5), (3, 4), (4, 10), (5, 12), (6, 2))]
    head_rank = dual_rank(head_cols, 13, "stageA.head_rank")
    pin(head_rank, 6, "stageA.head_full_rank_6")
    # obligation (c): a tied (non-all-distinct) sub-support witnesses that the
    # transport predicate detects glued columns: the diagonal 6-tie has rank 5
    # out of expected 6 because x=y folds the bidegree monomials to x^(i+j).
    diag_cols = [kronecker_column(x, x, 2, 4, 13) for x in (1, 2, 3, 4, 5, 6)]
    pin(dual_rank(diag_cols, 13, "stageA.tie_rank"), 5, "stageA.tie_collapses_rank")
    # circuit predicate on the all-distinct head (independent, not circuit)
    check("stageA.transport_circuit_independent",
          not is_circuit_checked(head_cols, 13, "stageA.transport"), head_rank)
    checkpoint({"stage": "A", "obligations": ["a", "b", "c"]})
    print("[stage A] transport library OK")

    # ---------------- Stage F: control battery (BEFORE censuses) ----------------
    control_facts = run_controls()
    checkpoint({"stage": "F", "controls": len(CHECKS)})
    print(f"[stage F] controls battery OK ({len(CHECKS)} asserts so far)")

    # ---------------- Stage B: PGL layer ----------------
    pgl_facts = {}
    for p in (7, 11, 13):
        mats = pgl_matrices(p)
        pin(len(mats), p * (p * p - 1), f"stageB.pgl_count.p{p}")
        pgl_facts[p] = len(mats)
    checkpoint({"stage": "B", "pgl_counts": pgl_facts})
    print(f"[stage B] PGL counts {pgl_facts}")

    # ---------------- Stage C: anchor sweeps ----------------
    c_records = []
    # (3,3) gate-9 re-expression: k=6, n=6,7 over three primes
    for p in (7, 11, 13):
        for n in (6, 7):
            pts = [x % p for x in range(1, n + 1)]
            domains = pgl_domains(p, pts, pts)
            rec = anchor_sweep(n, 3, 3, p, 6, domains, f"stageC.33.n{n}.p{p}")
            c_records.append(rec)
        want = {6: 12 if p == 7 else (4 if p == 11 else 2),
                7: 588 if p == 7 else (72 if p == 11 else 42)}
        checkpoint({"stage": "C", "part": "33row", "p": p})
        print(f"[stage C] (3,3) re-expression p={p}: {want}")
        for n in (6, 7):
            got = [r for r in c_records if r["n"] == n and r["p"] == p][0]["circuits"]
            pin(got, want[n], f"stageC.33.anchor.n{n}.p{p}")
    cpu, wall = elapsed()

    # (2,4) owner anchors: k=6 n=6 n=7; emptiness at k=5
    for p in (13, 11, 7):
        for n in (6, 7):
            pts = [x % p for x in range(1, n + 1)]
            domains = pgl_domains(p, pts, pts)
            rec = anchor_sweep(n, 2, 4, p, 6, domains, f"stageC.24.n{n}.p{p}")
            c_records.append(rec)
            want = {13: {6: 2, 7: 42}, 11: {6: 4, 7: 72},
                    7: {6: 12, 7: 588}}[p][n]
            got = rec["circuits"]
            pin(got, want, f"stageC.24.anchor.n{n}.p{p}")
        empty_want = {13: {6: 2, 7: 42}, 11: {6: 4, 7: 72}, 7: {6: 12, 7: 84}}[p]
        empt = sub_c7_sweep(6, 2, 4, p, 6, f"stageC.24.empty.p{p}")
        pin(empt[5], 0, f"stageC.24.k5_empty.p{p}")
        checkpoint({"stage": "C", "part": "24row", "p": p, "empty": empt})
        print(f"[stage C] (2,4) p={p}: n6={empty_want[6]} n7={empty_want[7]} "
              f"k5empty={empt[5]}")

    # (3,4) owner anchor: k=7 n=7; emptiness at k=6
    for p in (13, 11, 7):
        pts = [x % p for x in range(1, 8)]
        domains = pgl_domains(p, pts, pts)
        rec = anchor_sweep(7, 3, 4, p, 7, domains, f"stageC.34.n7.p{p}")
        c_records.append(rec)
        want = {13: 2, 11: 2, 7: 42}[p]
        pin(rec["circuits"], want, f"stageC.34.anchor.p{p}")
        empt = sub_c7_sweep(7, 3, 4, p, 7, f"stageC.34.empty.p{p}")
        pin(empt[6], 0, f"stageC.34.k6_empty.p{p}")
        checkpoint({"stage": "C", "part": "34row", "p": p, "record": rec})
        print(f"[stage C] (3,4) p={p}: circuits={want} k6empty={empt[6]}")

    # ---------------- Stage D: C2/size-9/10 grids ----------------
    d_records = []
    for p in (7, 11, 13):
        c9 = c2_grid_sweep(p, 4, 9)
        pin(c9, 0, f"stageD.grid4.size9.p{p}")
        c10 = c2_grid_sweep(p, 4, 10)
        pin(c10, 16, f"stageD.grid4.size10.p{p}")
        d_records.append({"grid": "4x4", "p": p, "size9": c9, "size10": c10})
        checkpoint({"stage": "D", "grid": "4x4", "p": p, "size9": c9, "size10": c10})
        print(f"[stage D] 4x4 GF{p}: size9={c9} size10={c10}")
    c9 = c2_grid_sweep(7, 5, 9)
    pin(c9, 36378, "stageD.grid5.size9.p7")
    checkpoint({"stage": "D", "grid": "5x5", "p": 7, "size9": c9})
    print(f"[stage D] 5x5 GF7: size9={c9} (CPU {elapsed()[0]:.1f}s)")
    c10 = c2_grid_sweep(7, 5, 10)
    pin(c10, 181960, "stageD.grid5.size10.p7")
    d_records.append({"grid": "5x5", "p": 7, "size9": c9, "size10": c10})
    checkpoint({"stage": "D", "grid": "5x5", "p": 7, "size10": c10})
    print(f"[stage D] 5x5 GF7: size10={c10}")

    # ---------------- Stage E: C3 pencils (3,3) size-8 ----------------
    e_records = []
    for p in (11, 13):
        # n=8 matching sweep, CI classification (gate-10 cross-validation)
        pts = [x % p for x in range(1, 9)]
        want_y_nodes = node_polynomial(pts, p)
        ci_set = 0
        rank7 = 0
        common_set = 0
        for perm in itertools.permutations(range(8)):
            points = tuple((pts[i], pts[perm[i]]) for i in range(8))
            cols8 = [kronecker_column(x, y, 3, 3, p) for x, y in points]
            r1 = dual_rank(cols8, p, f"stageE.p{p}")
            if r1 != 7:
                continue
            rank7 += 1
            basis = pencil_basis(cols8, p)
            if len(basis) != 2:
                raise AssertionError(f"stageE.p{p}: rank7 nullity != 2 at {perm}")
            common = common_component_info(basis[0], basis[1], p)
            if common["has_common_component"]:
                common_set += 1
            else:
                ci_set += 1
                if poly_monic(common["resultant_x"], p) != want_y_nodes or \
                   poly_monic(common["resultant_y"], p) != want_y_nodes:
                    raise AssertionError(f"stageE.p{p}: resultant node mismatch at {perm}")
            if perm[0] == 0 and perm[1:] == tuple(range(1, 8)):
                pass
            check_budget(f"stageE.p{p}.perm{len(e_records)}")
        pin(rank7, 2192 if p == 11 else 1344, f"stageE.rank7.p{p}")
        pin(ci_set, 560 if p == 11 else 416, f"stageE.ci.p{p}")
        e_records.append({"p": p, "rank7_pencils": rank7, "ci": ci_set,
                          "common_component": common_set})
        checkpoint({"stage": "E", "p": p, "rank7": rank7, "ci": ci_set,
                    "common": common_set})
        print(f"[stage E] n=8 GF{p}: rank7={rank7} CI={ci_set} common={common_set}")

    # ---------------- Stage F2: C3 pencil REJECT + no-component unit controls ----
    p = 13
    line = form_from_terms({(1, 0): 1, (0, 1): -1}, p)
    other = form_from_terms({(1, 0): 1, (0, 1): -2}, p)
    check("stageF2.common_component_detector",
          common_component_info(form_mul(line, line, p), form_mul(line, other, p),
                                 p)["has_common_component"], None)

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
            "rlimit_cpu": resource.getrlimit(resource.RLIMIT_CPU),
            "max_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        },
        "budget": {
            "CPU_seconds": cpu,
            "wall_seconds": wall,
            "CPU_cap_seconds": CPU_CAP,
            "wall_cap_seconds": WALL_CAP,
            "within_cap": True,
        },
        "pgl_counts": pgl_facts,
        "anchor_sweeps": c_records,
        "c2_grid_censuses": d_records,
        "c3_pencil_censuses": e_records,
        "control_facts": {k: (list(v) if isinstance(v, tuple) else v)
                          for k, v in control_facts.items()},
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
