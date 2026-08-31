"""task1_crosscheck.py — Task 1: independent decision procedures for gate B.

(a) Expand Sun's printed 56-addition SLP; run 729 Brent identities; identify
    which factor block his "transformed left" circuit realizes; cross-check
    against the V-floor conclusion.
(b) ILP route via HiGHS on the floor questions (U<=12, V<=13, W-factor<=13,
    and Sun's claimed transformed-left<=12).
(c) SAT route via python-sat CaDiCaL153 on the same.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from flint import fmpz

sys.path.insert(0, str(Path(__file__).parent))

from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED
from tensor_data_sun56 import (SUN_LEFT_SLP, SUN_RIGHT_SLP, SUN_PRODUCTS,
                               SUN_OUT_SHARED, SUN_C, SUN_BLOCKS_TEXT,
                               parse_blocks)
from gate_b_floor import canon, N, R, prep, subset_dfs
import gate_b_encodings as enc


def expand_env(slp, atom_env, input_prefix):
    env = dict(atom_env)
    for name, terms in slp:
        acc = None
        for sign, atom in terms:
            val = env[atom]
            val = val if sign == 1 else -val
            acc = val if acc is None else acc + val
        env[f"{input_prefix}{name}"] = acc
    return env


def sun_factors():
    R_ = R
    U = [[0] * R_ for _ in range(9)]
    V = [[0] * R_ for _ in range(9)]

    # expand left/right gates into A/B terms (symbolic over integers)
    def expand(slp, base_var):
        env = {}
        for i in range(9):
            env[f"{base_var}{i}"] = {f"{base_var}{i}": 1}
        for name, terms in slp:
            acc: dict[str, int] = {}
            for sign, atom in terms:
                src = env[atom]
                for k, v in src.items():
                    acc[k] = acc.get(k, 0) + sign * v
            env[name] = acc
        return env

    L = expand(SUN_LEFT_SLP, "A")
    Rr = expand(SUN_RIGHT_SLP, "B")
    for r, ((ls, la), (rs, rb)) in enumerate(SUN_PRODUCTS):
        lterms = L[la] if la.startswith("u") else {la: 1}
        rterms = Rr[rb] if rb.startswith("v") else {rb: 1}
        if ls == -1:
            lterms = {k: -v for k, v in lterms.items()}
        if rs == -1:
            rterms = {k: -v for k, v in rterms.items()}
        for a, c in lterms.items():
            U[int(a[1:])][r] += c
        for b, c in rterms.items():
            V[int(b[1:])][r] += c

    # output side
    out_env: dict[str, dict[int, int]] = {}
    for name, terms in SUN_OUT_SHARED:
        acc: dict[int, int] = {}
        for sign, atom in terms:
            idx = int(atom[1:])
            acc[idx] = acc.get(idx, 0) + sign
        out_env[name] = acc
    W = [[0] * R_ for _ in range(9)]
    for k, (cname, terms) in enumerate(SUN_C):
        # SUN_C order = C0..C8 (row-major)
        acc: dict[int, int] = {}
        for sign, atom in terms:
            if atom.startswith("M"):
                idx = int(atom[1:])
                acc[idx] = acc.get(idx, 0) + sign
            else:
                for j, v in out_env[atom].items():
                    acc[j] = acc.get(j, 0) + sign * v
        for midx, coef in acc.items():
            W[k][midx] += coef
    return U, V, W


def brent(U, V, W):
    fails = 0
    Uf = [[fmpz(x) for x in row] for row in U]
    Vf = [[fmpz(x) for x in row] for row in V]
    Wf = [[fmpz(x) for x in row] for row in W]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Ur = Uf[3 * i + k]
                for ip in range(3):
                    for jp in range(3):
                        Wr = Wf[3 * ip + jp]
                        for kp in range(3):
                            Vr = Vf[3 * kp + j]
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = Ur[r], Vr[r], Wr[r]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            expect = 1 if (i == ip and j == jp and k == kp) else 0
                            if int(s) != expect:
                                fails += 1
    return fails


def main():
    out = {}

    # ---- (a) Sun SLP -> factors; Brent identity test --------------------
    Usun, Vsun, Wsun = sun_factors()
    out["sun_brent_failures"] = brent(Usun, Vsun, Wsun)
    out["sun_ternary"] = all(x in (-1, 0, 1) for blk in (Usun, Vsun, Wsun)
                             for row in blk for x in row)

    blocks = parse_blocks(SUN_BLOCKS_TEXT)
    B_U, B_V, B_W = blocks
    out["sun_slp_vs_blocks"] = {
        "U": all(Usun[i][r] == B_U[i][r] for i in range(9) for r in range(R)),
        "V": all(Vsun[i][r] == B_V[i][r] for i in range(9) for r in range(R)),
        "W": all(Wsun[i][r] == B_W[i][r] for i in range(9) for r in range(R)),
    }

    # is Sun's U-block the 55-paper's V-block (cyclic sigma)? test all
    # mappings of the three 55-paper blocks onto Sun's U block (with
    # permutations of columns? no: columns are aligned by product pairing,
    # but Sun's product ORDER may differ. Compare as SETS of column vectors.)
    paper_cols = {
        "U": {canon(tuple(U_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)},
        "V": {canon(tuple(V_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)},
        "W": {canon(tuple(W_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)},
    }
    for name, blk in (("U", B_U), ("V", B_V), ("W", B_W)):
        cols = {canon(tuple(blk[i][r] for i in range(9))) for r in range(R)}
        out[f"sun_{name}_block_as_set_equals"] = {
            pn: (cols == pc) for pn, pc in paper_cols.items()
        }

    # which block does Sun's LEFT SLP realize?
    L = {}
    for i in range(9):
        L[f"A{i}"] = {f"A{i}": 1}
    def expand(slp, base):
        env = {}
        for i in range(9):
            env[f"{base}{i}"] = {f"{base}{i}": 1}
        for name, terms in slp:
            acc = {}
            for sign, atom in terms:
                for k, v in env[atom].items():
                    acc[k] = acc.get(k, 0) + sign * v
            env[name] = acc
        return env
    Lexp = expand(SUN_LEFT_SLP, "A")
    realized = []
    for r, ((ls, la), _) in enumerate(SUN_PRODUCTS):
        t = Lexp[la] if la.startswith("u") else {la: 1}
        if ls == -1:
            t = {k: -v for k, v in t.items()}
        vec = tuple(t.get(f"A{i}", 0) for i in range(9))
        realized.append(vec)
    out["sun_left_target_set"] = sorted({canon(v) for v in realized if any(v)})

    # compare against paper V block columns (as sets)
    v_cols = {canon(tuple(V_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)}
    u_cols = {canon(tuple(U_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)}
    w_cols = {canon(tuple(W_BLOCK_PRINTED[i][r] for i in range(9))) for r in range(R)}
    s = set(out["sun_left_target_set"])
    out["sun_left_meets"] = {"U_cols": s == u_cols, "V_cols": s == v_cols, "W_cols": s == w_cols}
    out["sun_left_subset_of"] = {"U": s <= u_cols, "V": s <= v_cols, "W": s <= w_cols}

    # ---- (b) HiGHS ILP on the floor questions ---------------------------
    res_ilp = {}
    for name, blk in (("U", U_BLOCK_PRINTED), ("V", V_BLOCK_PRINTED), ("Wf", W_BLOCK_PRINTED)):
        targets = [tuple(blk[i][r] for i in range(9)) for r in range(R)]
        classes, reps = prep(targets)
        d = len(classes)
        try:
            h, st, nv, cons, aux = enc.solve_floor_ilp(classes, reps, d)
            stname = str(st)
            res_ilp[name] = {"d": d, "T": d, "highs_status": stname}
        except Exception as e:  # noqa: BLE001
            res_ilp[name] = {"d": d, "error": repr(e)}
    out["ilp_floor"] = res_ilp

    # ---- (c) python-sat CNF on the floor questions ----------------------
    res_sat = {}
    from pysat.solvers import Cadical153
    for name, blk in (("U", U_BLOCK_PRINTED), ("V", V_BLOCK_PRINTED), ("Wf", W_BLOCK_PRINTED)):
        targets = [tuple(blk[i][r] for i in range(9)) for r in range(R)]
        classes, reps = prep(targets)
        d = len(classes)
        clauses, varmap, nv = enc.build_floor_cnf(classes, reps, d)
        with Cadical153(bootstrap_clauses=clauses) as s:
            sat = s.solve()
            res_sat[name] = {"d": d, "T": d, "SAT": bool(sat),
                             "n_vars": nv, "n_clauses": len(clauses)}
    out["sat_floor"] = res_sat

    print(json.dumps(out, indent=2, default=str))
    Path(__file__).parent.joinpath("task1_result.json").write_text(
        json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
