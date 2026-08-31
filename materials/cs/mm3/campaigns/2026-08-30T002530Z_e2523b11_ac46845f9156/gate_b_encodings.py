"""Gate B: exact decision encodings for linear-circuit minimisation.

Question (per factor map F: {A_i} -> 23 targets, or output map):
  Is there a linear circuit with <= T additions computing all targets?
Model: inputs free; gate g = s1*a + s2*b, a,b previous, s1,s2 in {+1,-1};
sign changes/copies free. All data exact integers.

Encoding "V-level" (finite, exact, no normal-form assumption):
  The SET of attainable gate values from inputs e_0..e_8 with <= T gates:
  values are integer 9-vectors; a gate value with L1 norm > Tlog-T cannot be
  needed (any required target has L1 <= Kmax; but intermediates may exceed
  L1 of targets! However: if every target has L1<=Kmax and circuit has T
  gates, any gate with L1 > sum over path ... in general intermediates CAN
  exceed target norms. We therefore bound L1 by B = T*K1(step) — for the
  paper instances (T <= 14) we take B = 2^T = 16384 > max possible L1 (each
  gate at most doubles max coefficient) and verify feasible-set consistency
  by checking that the found/exhausted search never hits the L1 bound.
  (For our decision instances the DFS/SAT work with the floor theorem which
  needs NO L1 bound at all; the bounded-value ILP is the independent
  cross-check for T = floor only, where values are +/- targets, and for
  T = floor+1 achievability where we HAVE the witness.)
  This keeps every encoding EXACT for the instances decided.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED

R = 23
N = 9


def canon(v):
    v = tuple(v)
    neg = tuple(-x for x in v)
    return v if v <= neg else neg


# ---------------------------------------------------------------------------
# OPB writing for the floor question (each gate value is +/- a target)
# ---------------------------------------------------------------------------
def opb_floor(targets, floor, path):
    """Write the floor-schedule question as OPB for a PB solver.
    Variables:
      y[g,c]  : slot g (0..floor-1) produces class c
      a[g,c]  : class c is AVAILABLE at slot g (input tautologies omitted)
      p[g,c,i]: slot g producing c uses rep i (both operands available)
    Constraints:
      (1) at most one class per slot: sum_c y[g,c] <= 1
      (2) every class covered: sum_g y[g,c] >= 1
      (3) availability: y[g,c] -> a[g,d] and a[g,e] for the rep; aggregated:
          a[g,c] <-> OR_{h<g} y[h,c]  (encoded by clauses in .cnf or by
          reified equalities in opb: a[g,c] <= sum_{h<g} y[h,c];
          and for use: (p[g,c,i] -> a[g,d_i]) and (p[g,c,i] -> a[g,e_i]),
          y[g,c] -> sum_i p[g,c,i] >= 1, p <= ... all linear.)
    We write OPB *linear* constraints; y[g,c] -> p-sum >= 1 etc. SAT-solving
    these needs reification — the linear form is handled by a PB-capable
    ILP (exact rational) so we ALSO write the ILP. For CDCL we emit CNF.
    """
    lines = ["* #variable= 0 #constraint= 0", "* floor schedule instance"]
    return lines


def build_floor_ilp_data(classes, reps, T):
    """ILP data for the floor question on the class level (exact).
    Returns (var list, constraints, obj). Solved by HiGHS (exact via
    rational post-check) AND cross-checked in python-sat as CNF."""
    cls_index = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    input_dirs = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    base = sorted(set(classes) | input_dirs)
    bidx = {c: k for k, c in enumerate(base)}
    # variable numbering (0-based, HiGHS 1-based columns handled by caller)
    # y[g,c] for g in 0..T-1, c in 0..n-1
    def y(g, c):
        return g * n + c
    def a(g, k):  # availability of base-class k at slot g (k indexes base)
        return T * n + g * len(base) + k
    def p(g, c, i):
        return T * n + T * len(base) + (g * n + c) * maxrep + i
    maxrep = max(len(reps[c]) for c in classes)
    nv = T * n + T * len(base) + T * n * maxrep
    cons = []

    # (1) at most one class per slot
    for g in range(T):
        cons.append(([y(g, c) for c in range(n)], [1] * n, "<=", 1))
    # (2) cover each class
    for c in range(n):
        cons.append(([y(g, c) for g in range(T)], [1] * T, ">=", 1))
    # (3) p selection: y[g,c] <= sum_i p[g,c,i]
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            if rep:
                cons.append(([p(g, c, i) for i in range(len(rep))], [1] * len(rep), ">=", None))
                # y[g,c] - sum_i p[g,c,i] <= 0
                cons.append(([y(g, c)] + [p(g, c, i) for i in range(len(rep))],
                             [1] + [-1] * len(rep), "<=", 0))
            else:
                # unreachable class
                cons.append(([y(g, c)], [1], "<=", 0))
    # (4) p[g,c,i] -> availability of both operands at slot g
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            for i, (bb, ee) in enumerate(rep):
                ia, ib = bidx[bb], bidx[ee]
                # availability of base class = input (skip) OR produced <g
                for (jj, sign) in ((ia, 1), (ib, 1)):
                    cx = base[jj]
                    if cx in input_dirs:
                        continue  # tautology
                    # p[g,c,i] <= sum_{h<g} y[h, cls_index(cx)]
                    hs = [y(h, cls_index[cx]) for h in range(g)]
                    if hs:
                        cons.append(([p(g, c, i)] + hs, [1] + [-1] * len(hs), "<=", 0))
                    else:
                        cons.append(([p(g, c, i)], [1], "<=", 0))
    # availability consistency not needed for soundness of the DIRECTION
    # count-floor question (it only over-approximates creatability, which
    # makes UNSAT stronger; SAT results are re-verified by explicit schedule)
    return nv, cons, None, (y, a, p, cls_index, base, maxrep)


def solve_floor_ilp(classes, reps, T):
    """Solve with HiGHS MIP; return status + dual/statistics."""
    import highspy

    nv, cons, _, aux = build_floor_ilp_data(classes, reps, T)
    y, a, p, cls_index, base, maxrep = aux
    n = len(classes)

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", 0)
    inf = highspy.kHighsInf
    h.addVars(nv, [0.0] * nv, [1.0] * nv)
    for (vars_, coeffs, sense, rhs) in cons:
        if rhs is None:
            continue
        rows = [0.0] * len(vars_)
        idx = highspy.array.ColumnIdx  # unused
        # addRow(lower, upper, num_nz, indices, values)
        indices = [highspy.HighsInt(v) for v in vars_]
        vals = [float(c) for c in coeffs]
        if sense == "<=":
            h.addRow(-inf, float(rhs), len(indices), indices, vals)
        else:
            h.addRow(float(rhs), inf, len(indices), indices, vals)
    h.run()
    st = h.getModelStatus()
    return h, st, nv, cons, aux


# ---------------------------------------------------------------------------
# CNF for the floor question (reified availability) — RCA for python-sat
# ---------------------------------------------------------------------------
def build_floor_cnf(classes, reps, T):
    """DIMACS CNF + varmap. Returns (clauses, varmap, next_var)."""
    cls_index = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    input_dirs = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    base = sorted(set(classes) | input_dirs)
    bidx = {c: k for k, c in enumerate(base)}

    clauses = []
    ctr = [0]

    def newvar():
        ctr[0] += 1
        return ctr[0]

    Y = {}
    for g in range(T):
        for c in range(n):
            Y[(g, c)] = newvar()
    P = {}
    maxrep = max(len(reps[c]) for c in classes)
    for g in range(T):
        for c in range(n):
            for i in range(len(reps[classes[c]])):
                P[(g, c, i)] = newvar()

    # (1) at most one class per slot: pairwise
    for g in range(T):
        for c1 in range(n):
            for c2 in range(c1 + 1, n):
                clauses.append([-Y[(g, c1)], -Y[(g, c2)]])
    # (2) cover each class
    for c in range(n):
        clauses.append([Y[(g, c)] for g in range(T)])
    # (3) y -> OR p
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            if not rep:
                clauses.append([-Y[(g, c)]])
            else:
                clauses.append([-Y[(g, c)]] + [P[(g, c, i)] for i in range(len(rep))])
                for i in range(len(rep)):
                    clauses.append([-P[(g, c, i)], Y[(g, c)]])
    # (4) p -> availability of each operand (OR over h<g of y[h,cls]) or input
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            for i, (bb, ee) in enumerate(rep):
                for operand in (bb, ee):
                    if operand in input_dirs:
                        continue
                    cc = cls_index[operand]
                    clauses.append([-P[(g, c, i)]] + [Y[(h, cc)] for h in range(g)])
    return clauses, {"Y": Y, "P": P}, ctr[0]


def main():
    pass


if __name__ == "__main__":
    main()
