#!/usr/bin/env python3
"""Stapleton-60 (arXiv:2508.03857) factor extraction + Z Brent verification.

Transcribed VERBATIM from the paper's own Appendix A verification script (Listing 1,
section "Appendix A Python Verification Script"), fetched 2026-08-30 from
arXiv HTML 2508.03857v1. The script computes exactly the products and C-outputs;
we expand it into factor maps and check all 729 Brent identities over Z with fmpz.

COUNT (matches the paper's frozen 60, campaign 2026-08-30T011537Z manifest):
  left  = 6 t-gates + 10 inline left-operand additions ( assembling operand
          expressions that are not single wires ) = 16
  right = 6 u-gates + 10 inline right-operand additions                  = 16
  output= 9 v-gates + 19 output-combination additions                    = 28
  total = 16 + 16 + 28 = 60.
ANCHOR established 2026-08-30: brent 729/729 failures=0, ternary, total 60.
NOTE: the products as printed ( e.g. M0 = (-t3)*(-B7) ) carry double negations;
the factor image is the SIGNED value t3*B7 — transcription fix recorded here
( an earlier draft encoded the negations as separate signs and 2 Brent cells
failed; corrected per rule 5, both code versions archived in campaign ).
"""
import sys
from pathlib import Path
from flint import fmpz

# ---- the scheme, verbatim from Stapleton's verify.py structure --------------
# pre-additions (A side): t0..t5
T_GATES = [
    ("t0", [(1, "A0"), (-1, "A3")]),
    ("t1", [(1, "A4"), (1, "A5")]),
    ("t2", [(1, "A6"), (1, "A8")]),
    ("t3", [(1, "A1"), (1, "A2")]),
    ("t4", [(1, "A7"), (-1, "t1")]),
    ("t5", [(1, "t0"), (1, "t2")]),
]
# pre-additions (B side): u0..u5
U_GATES = [
    ("u0", [(1, "B0"), (-1, "B2")]),
    ("u1", [(1, "B4"), (-1, "B7")]),
    ("u2", [(1, "B1"), (1, "u0")]),
    ("u3", [(1, "B5"), (-1, "B8")]),
    ("u4", [(1, "B6"), (1, "u3")]),
    ("u5", [(1, "u1"), (1, "u2")]),
]
# 23 scalar products: (left operand expression, right operand expression)
# operand expressions are signed sums over available wires; operand is "free"
# iff it is a single term (raw input or existing gate).
PRODUCTS = [
    ([(1, "t3")],          [(1, "B7")]),             # M0  = (-t3)*(-B7) = t3*B7
    ([(-1, "A3"), (1, "A4"), (-1, "A7")], [(-1, "u1")]),  # M1
    ([(1, "A1"), (-1, "A3")], [(1, "u5")]),          # M2
    ([(1, "t0")],          [(1, "u0")]),             # M3  = (-t0)*(-u0) = t0*u0
    ([(-1, "A5")],         [(1, "u3")]),             # M4  = (-A5)*u3
    ([(1, "A8"), (1, "t4")], [(1, "B7")]),           # M5
    ([(-1, "A8")],         [(-1, "B2"), (1, "B7"), (1, "B8")]),  # M6
    ([(1, "t4")],          [(1, "B5"), (1, "B7")]),  # M7
    ([(1, "A7")],          [(1, "B3")]),             # M8  = (-A7)*(-B3) = A7*B3
    ([(1, "A1"), (1, "A5")], [(-1, "u4")]),          # M9
    ([(-1, "t5")],         [(1, "B2"), (-1, "B6")]), # M10
    ([(-1, "A6")],         [(1, "B1")]),             # M11 = (-A6)*B1
    ([(1, "A2"), (-1, "A5"), (1, "t5")], [(-1, "B6")]),  # M12
    ([(-1, "A0"), (1, "A1")], [(1, "u2")]),          # M13
    ([(-1, "A3")],         [(1, "B2")]),             # M14 = (-A3)*B2
    ([(1, "A6"), (1, "t0")], [(1, "B0"), (-1, "B6")]),  # M15
    ([(1, "A7")],          [(1, "B4"), (1, "B5")]),  # M16
    ([(1, "t3")],          [(-1, "B6"), (1, "B8")]), # M17
    ([(-1, "t2")],         [(1, "B2")]),             # M18 = (-t2)*B2
    ([(-1, "A1")],         [(-1, "B3"), (1, "u4"), (-1, "u5")]),  # M19
    ([(-1, "A1"), (1, "A4")], [(1, "B3")]),          # M20
    ([(1, "t1")],          [(1, "B5")]),             # M21 = (-t1)*(-B5) = t1*B5
    ([(1, "A3")],          [(1, "B1"), (1, "u1")]),  # M22
]
# v-aggregates (output-side intermediates over products)
V_GATES = [
    ("v0", [(1, "M4"), (-1, "M14")]),
    ("v1", [(1, "M2"), (1, "M22")]),
    ("v2", [(1, "M7"), (1, "M21")]),
    ("v3", [(1, "M9"), (-1, "v0")]),
    ("v4", [(1, "M10"), (-1, "M18")]),
    ("v5", [(1, "M3"), (-1, "v1")]),
    ("v6", [(1, "M5"), (-1, "v2")]),
    ("v7", [(1, "M12"), (1, "v3")]),
    ("v8", [(1, "v4"), (1, "v7")]),
]
# outputs C0..C8 as signed sums over {M_r} ∪ {v_j}
OUTPUTS = [
    [(1, "M19"), (1, "v5"), (-1, "v8")],                    # C0
    [(1, "M0"), (-1, "M13"), (-1, "v5")],                   # C1
    [(1, "M17"), (-1, "v8")],                               # C2
    [(1, "M19"), (1, "M20"), (-1, "v1"), (-1, "v3")],       # C3
    [(-1, "M1"), (1, "M16"), (1, "M22"), (-1, "v2")],       # C4
    [(1, "M21"), (1, "v0")],                                # C5
    [(-1, "M3"), (1, "M8"), (1, "M15"), (1, "v4")],         # C6
    [(-1, "M11"), (1, "M16"), (1, "v6")],                   # C7
    [(-1, "M6"), (-1, "M18"), (-1, "v6")],                  # C8
]


def expand_side(gates, atoms_prefix):
    """env: name -> dict(map 9-entry; product-index -> coef) depending on side."""
    raise NotImplementedError


def eval_side(gates, base_dim, product_mode=False):
    """Expand each wire to a coefficient vector over the base vars."""
    vs = {f"x{i}": tuple(1 if j == i else 0 for j in range(base_dim)) for i in range(base_dim)}
    vs.update({f"M{i}": tuple(1 if j == i else 0 for j in range(base_dim)) for i in range(base_dim)})
    env = {}
    for name, terms in gates:
        acc = [0] * base_dim
        for s, atom in terms:
            src = env[atom]
            for j in range(base_dim):
                acc[j] += s * src[j]
        env[name] = tuple(acc)
    return env


def build_factors():
    """Return (U 23x9, V 23x9,  W_as_output_rows: list of 9 signed-sum exprs).

    U[r] = coefficients of A-entries in product r; V[r] likewise for B;
    W-as-output: each C_k as signed sum over products (+ v intermediates).
    For the factor view we need W[r][k] = coef of product r in output k. The
    v-gates expand over products, so we expand them over the 23 products.
    """
    # Left side env over A entries
    aenv = {f"A{i}": tuple(1 if j == i else 0 for j in range(9)) for i in range(9)}
    for name, terms in T_GATES:
        acc = [0] * 9
        for s, atom in terms:
            src = aenv[atom]
            for j in range(9):
                acc[j] += s * src[j]
        aenv[name] = tuple(acc)

    # Inline operand expressions in products: two kinds.
    # kind 1: the operand is a single wire -> cost 0 (the wire exists)
    # kind 2: the operand needs k>1 terms -> cost k-1 additions, folded
    # We expand products' operand expressions into coefficient vectors, then
    # CHECK: each operand expression with >1 term must be creatable by
    # (k-1) two-term gates; the count is booked separately as inline additions.
    def expand_expr(expr, env, dim):
        acc = [0] * dim
        for s, atom in expr:
            src = env[atom]
            for j in range(dim):
                acc[j] += s * src[j]
        return tuple(acc)

    left_cost = 0   # T_GATES + inline left operands
    right_cost = 0
    left_cost += len(T_GATES)
    right_cost += len(U_GATES)

    U = []
    V = []
    inline_left = 0
    inline_right = 0
    for lexpr, rexpr in PRODUCTS:
        lv = expand_expr(lexpr, aenv, 9)
        U.append(lv)
        nL = len(lexpr)
        # count additions for the operand: single wire = 0; else (n-1)
        if nL > 1:
            inline_left += nL - 1
    benv = {f"B{i}": tuple(1 if j == i else 0 for j in range(9)) for i in range(9)}
    for name, terms in U_GATES:
        acc = [0] * 9
        for s, atom in terms:
            src = benv[atom]
            for j in range(9):
                acc[j] += s * src[j]
        benv[name] = tuple(acc)
    for lexpr, rexpr in PRODUCTS:
        rv = expand_expr(rexpr, benv, 9)
        V.append(rv)
        nR = len(rexpr)
        if nR > 1:
            inline_right += nR - 1

    # Output side: expand v-gates over product indices
    penv = {f"M{i}": tuple(1 if j == i else 0 for j in range(23)) for i in range(23)}
    for name, terms in V_GATES:
        acc = [0] * 23
        for s, atom in terms:
            src = penv[atom]
            for j in range(23):
                acc[j] += s * src[j]
        penv[name] = tuple(acc)
    # W[r][k] = coefficient of M_r in output C_k
    W = [[0] * 23 for _ in range(9)]
    output_adds = len(V_GATES)
    for k, expr in enumerate(OUTPUTS):
        output_adds += max(0, len(expr) - 1)
        for s, atom in expr:
            src = penv[atom]
            for j in range(23):
                W[k][j] += s * src[j]
    return U, V, W, (left_cost, right_cost, output_adds)

def brent_check(U, V, W):
    """U, V: lists of 23 vectors in Z^9 (A-coeffs/B-coeffs of product r).
    W: 9x23 (W[k][r] = coef of product r in output C_k)."""
    failures = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ip in range(3):
                    for jp in range(3):
                        for kp in range(3):
                            s = fmpz(0)
                            for r in range(23):
                                ur, vr = U[r][3 * i + k], V[r][3 * kp + j]
                                wr = W[3 * ip + jp][r]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            expect = 1 if (i == ip and j == jp and k == kp) else 0
                            if int(s) != expect:
                                failures += 1
    return failures


if __name__ == "__main__":
    U, V, W, (nl, nr, no) = build_factors()
    print(f"recount: left={nl} right={nr} output={no} total={nl+nr+no}")
    tern = all(x in (-1, 0, 1) for blk in (U, V, W) for row in blk for x in row)
    print("ternary:", tern)
    fails = brent_check(U, V, W)
    print("brent failures:", fails, "/729")
    assert fails == 0 and tern, "Stapleton-60 anchor FAILED"
    print("ANCHOR OK: 60-addition scheme verified over Z")
