"""task1_enc_only.py — Task 1 cross-encodings (ILP HiGHS + SAT python-sat)
on the three 55-paper floor questions. Independent of the DFS module logic
except for the shared class/REP construction (prep), which is a plain
integer computation."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED
from gate_b_floor import prep


def ilp_floor(classes, reps, T):
    """Exact ILP of the floor question solved by HiGHS.
    same encoding semantics as build_floor_ilp_data (=rebuilt here to fix
    an old stub). Mirrors the DFS model exactly:
      y[g,c] slot g produces class c; at most one per slot (the floor means
      every used gate makes a new needed class — unused slots allowed);
      every class covered; slot g may only produce c if a representation
      (d,e) exists with d,e available OR produced at slots < g (reified)."""
    import highspy

    cls_index = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    input_dirs = {canon9(i) for i in range(9)}
    base = sorted(set(classes) | input_dirs)
    bidx = {c: k for k, c in enumerate(base)}

    def y(g, c):
        return g * n + c

    maxrep = max(len(reps[c]) for c in classes)
    # p variable per (slot, class, rep-id)
    P = {}
    nextid = T * n
    for g in range(T):
        for c in range(n):
            for i in range(len(reps[classes[c]])):
                P[(g, c, i)] = nextid
                nextid += 1
    nv = nextid

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", 0)
    inf = highspy.kHighsInf
    h.addVars(nv, [0.0] * nv, [1.0] * nv)
    for i in range(nv):
        h.changeColIntegrality(i, highspy.HighsVarType.kInteger)

    def row(vars_, coeffs, lo, hi):
        h.addRow(lo, hi, len(vars_), [int(v) for v in vars_],
                 [float(c) for c in coeffs])

    # at most one class per slot
    for g in range(T):
        row([y(g, c) for c in range(n)], [1] * n, -inf, 1.0)
    # coverage
    for c in range(n):
        row([y(g, c) for g in range(T)], [1] * T, 1.0, inf)
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            if not rep:
                row([y(g, c)], [1], -inf, 0.0)
                continue
            # y <= sum p
            row([y(g, c)] + [P[(g, c, i)] for i in range(len(rep))],
                [1] + [-1] * len(rep), -inf, 0.0)
            for i, (bb, ee) in enumerate(rep):
                for opnd in (bb, ee):
                    opc = canon9t(opnd)
                    if any(all(abs(a) == b for a, b in zip(opc, e)) for e in input_dirs):
                        continue
                    cc = cls_index[opc]
                    if g == 0:
                        row([P[(g, c, i)]], [1], -inf, 0.0)
                    else:
                        row([P[(g, c, i)]] + [y(h2, cc) for h2 in range(g)],
                            [1] + [-1] * g, -inf, 0.0)
    h.run()
    status = h.getModelStatus()
    name = h.modelStatusToString(status)
    return name


def canon9(i):
    v = [0] * 9
    v[i] = 1
    return tuple(v)


def canon9t(v):
    v = tuple(v)
    neg = tuple(-x for x in v)
    return v if v <= neg else neg


def cnf_floor(classes, reps, T):
    """DIMACS CNF of the same floor question; solved by python-sat CaDiCaL153."""
    from pysat.solvers import Cadical153

    cls_index = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    input_dirs = {canon9(i) for i in range(9)}
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
    for g in range(T):
        for c in range(n):
            for i in range(len(reps[classes[c]])):
                P[(g, c, i)] = newvar()

    for g in range(T):
        for c1 in range(n):
            for c2 in range(c1 + 1, n):
                clauses.append([-Y[(g, c1)], -Y[(g, c2)]])
    for c in range(n):
        clauses.append([Y[(g, c)] for g in range(T)])
    for g in range(T):
        for c in range(n):
            rep = reps[classes[c]]
            if not rep:
                clauses.append([-Y[(g, c)]])
                continue
            clauses.append([-Y[(g, c)]] + [P[(g, c, i)] for i in range(len(rep))])
            for i in range(len(rep)):
                clauses.append([-P[(g, c, i)], Y[(g, c)]])
                for opnd in rep[i]:
                    opc = canon9t(opnd)
                    if any(all(abs(a) == b for a, b in zip(opc, e)) for e in input_dirs):
                        continue
                    cc = cls_index[opc]
                    if g == 0:
                        clauses.append([-P[(g, c, i)]])
                    else:
                        clauses.append([-P[(g, c, i)]] + [Y[(h2, cc)] for h2 in range(g)])
    s = Cadical153()
    for cl in clauses:
        s.add_clause(cl)
    sat = s.solve()
    s.delete()
    return bool(sat), ctr[0], len(clauses)


def main():
    out = {}
    for name, blk in (("U", U_BLOCK_PRINTED), ("V", V_BLOCK_PRINTED), ("Wf", W_BLOCK_PRINTED)):
        targets = [tuple(blk[i][r] for i in range(9)) for r in range(23)]
        classes, reps = prep(targets)
        d = len(classes)
        ilp_st = ilp_floor(classes, reps, d)
        sat, nv, nc = cnf_floor(classes, reps, d)
        out[name] = {
            "d(F)": d,
            "ILP_HiGHS_at_T=d": {"status": ilp_st,
                                  "feasible": "Optimal" in ilp_st and None,
                                  "interpretation": "INFEASIBLE if status=Infeasible"},
            "SAT_CaDiCaL153_at_T=d": {"SAT": sat, "vars": nv, "clauses": nc,
                                       "interpretation": "UNSAT iff no schedule"},
        }
    print(json.dumps(out, indent=2, default=str))
    Path(__file__).parent.joinpath("task1_enc_result.json").write_text(
        json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
