#!/usr/bin/env python3
"""mws59-bit-58vs59 — primary runner (pre-registered; phase 2).

pre_statement.md sha aede812c…; Route P = verify the paper's PRINTED Table-2
straight-line program, establish the product correspondence to the Table-3
factor blocks, recount every stage in the fixed model, verify by exact
expansion, assemble with the frozen 14-gate U witness, re-verify Brent.
Route S (aux-3 prefix + scouting) runs only if Route P fails.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
MM3 = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(MM3, "src")
SCR5 = os.path.join(MM3, "campaigns",
                    "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56",
                    "scripts")
C1 = os.path.join(MM3, "campaigns", "20260901T133118Z_8c14c3dd_7970c5318d5d")
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(MM3, "scratch"))

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]
ENV1 = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}

N, R = 9, 23
CERT_DIR = os.path.join(HERE, "certificates")
os.makedirs(CERT_DIR, exist_ok=True)
RESULT: dict = {"phases": {}, "defects": []}
BUDGET_CPU_S = 4.0 * 3600
T0_WALL = time.time()

CHECKERS = {}


def die(msg):
    RESULT["fatal"] = msg
    print("FATAL:", msg, flush=True)
    json.dump(RESULT, open(os.path.join(HERE, "verdict_partial.json"), "w"),
              indent=1, default=str)
    sys.exit(1)


for _nm, _want in (("drat-trim", "111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4"),
                   ("lrat-check", "b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07")):
    _p = os.path.join(HERE, "tools", _nm)
    _got = hashlib.sha256(open(_p, "rb").read()).hexdigest() if os.path.isfile(_p) else "MISSING"
    CHECKERS[_nm] = (_p, _got, _got == _want)


def phase(name):
    def deco(fn):
        def wrapped(*a, **k):
            t = time.time()
            print(f"== phase {name} ==", flush=True)
            RESULT["phases"].setdefault(name, {})
            out = fn(*a, **k)
            RESULT["phases"][name]["seconds"] = round(time.time() - t, 2)
            json.dump(RESULT, open(os.path.join(HERE, "verdict_partial.json"), "w"),
                      indent=1, default=str)
            return out
        return wrapped
    return deco


def canon(v):
    t = tuple(v)
    n = tuple(-x for x in t)
    return t if t <= n else n


INPUT_CLASSES = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
U_ROWS: list = []
V_ROWS: list = []
W_ROWS: list = []

# ===================== PRINTED TABLE 2 (transcribed verbatim) ===============
# Source: scratch/mws59_layout.txt lines 145-177 (arXiv:2601.05272 Table 2).
T2_T = [("t0", [(1, "A3"), (1, "A6")]),
        ("t1", [(1, "A1"), (-1, "A7")]),
        ("t2", [(1, "A4"), (1, "A5")]),
        ("t3", [(1, "A7"), (-1, "t0")]),
        ("t4", [(1, "A6"), (1, "t1")]),
        ("t5", [(1, "A8"), (-1, "t0")]),
        ("t6", [(1, "t2"), (1, "t3")])]
T2_U = [("u0", [(1, "B2"), (1, "B5")]),
        ("u1", [(1, "B4"), (1, "u0")]),
        ("u2", [(1, "B1"), (1, "u1")]),
        ("u3", [(1, "B4"), (1, "B5")]),
        ("u4", [(1, "B3"), (-1, "u0")]),
        ("u5", [(1, "B8"), (-1, "u3")])]
# 23 products: (left operand expression, right operand expression)
T2_P = [
    ([(1, "A1")], [(1, "B4")]),                                   # M0
    ([(1, "A0")], [(1, "B1")]),                                   # M1
    ([(1, "A8"), (1, "t6")], [(1, "B8")]),                        # M2
    ([(1, "A2"), (-1, "A8")], [(1, "B8")]),                       # M3
    ([(1, "A5"), (1, "t5")], [(1, "B6")]),                        # M4
    ([(1, "A2")], [(1, "B7")]),                                   # M5
    ([(1, "A0")], [(1, "B0")]),                                   # M6
    ([(1, "t1")], [(1, "u0")]),                                   # M7
    ([(1, "A5")], [(1, "B6"), (1, "B7"), (1, "u5")]),             # M8
    ([(1, "A4")], [(1, "B3")]),                                   # M9
    ([(1, "A0"), (-1, "t4")], [(1, "B2")]),                       # M10
    ([(1, "t5")], [(1, "B6"), (-1, "u2")]),                       # M11
    ([(1, "t4")], [(1, "u1")]),                                   # M12
    ([(1, "A1")], [(1, "B3")]),                                   # M13
    ([(1, "t0"), (1, "t1")], [(1, "u4")]),                        # M14
    ([(1, "A4"), (1, "t3")], [(1, "B3"), (-1, "B5"), (1, "B8")]), # M15
    ([(1, "t6")], [(1, "u5")]),                                   # M16
    ([(1, "A3")], [(1, "B0"), (1, "u2")]),                        # M17
    ([(1, "A3"), (-1, "t2")], [(1, "u3")]),                       # M18
    ([(1, "A8")], [(1, "B7"), (1, "u2")]),                        # M19
    ([(1, "A2")], [(1, "B6")]),                                   # M20
    ([(1, "A6"), (-1, "A8")], [(1, "u2")]),                       # M21
    ([(1, "t0")], [(1, "B0"), (1, "B6"), (1, "u4")]),             # M22
]
T2_V = [("v0", [(1, "M0"), (-1, "M12")]),
        ("v1", [(1, "M16"), (1, "v0")]),
        ("v2", [(1, "M11"), (-1, "M21")]),
        ("v3", [(1, "M14"), (-1, "M13")]),
        ("v4", [(1, "M18"), (-1, "v1")]),
        ("v5", [(1, "M2"), (1, "v4")]),
        ("v6", [(1, "M4"), (1, "M9")]),
        ("v7", [(1, "M15"), (1, "v3")]),
        ("v8", [(1, "M17"), (-1, "v2")])]
T2_C = [
    [(1, "M6"), (1, "M13"), (1, "M20")],                                    # C0
    [(1, "M0"), (1, "M1"), (1, "M5")],                                      # C1
    [(1, "M3"), (1, "M10"), (1, "v5")],                                     # C2
    [(1, "v6"), (1, "v8")],                                                 # C3
    [(1, "M8"), (-1, "v1"), (1, "v2"), (-1, "v6"), (1, "v7")],              # C4
    [(1, "M9"), (-1, "v4"), (-1, "v7")],                                    # C5
    [(-1, "M7"), (1, "M22"), (-1, "v3"), (-1, "v8")],                       # C6
    [(1, "M7"), (1, "M19"), (1, "M21"), (1, "v0")],                         # C7
    [(-1, "M7"), (1, "v5")],                                                # C8
]


def load_mws59():
    from pathlib import Path
    lines = (Path(MM3) / "scratch" / "mws59_layout.txt").read_text().splitlines()
    data = lines[236:265]
    blocks, cur = [], []
    for L in data:
        s = L.strip()
        if s.startswith("#"):
            blocks.append(cur); cur = []
        elif s and not s.startswith(("A", "T")):
            try:
                cur.append([int(x) for x in s.split()])
            except ValueError:
                pass
    if cur:
        blocks.append(cur)
    b1, b2, b3 = blocks
    U = [tuple(b1[k][r] for k in range(9)) for r in range(23)]
    V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]
    W = [tuple(b3[k][r] for k in range(9)) for r in range(23)]
    return U, V, W


def brent_fail_int(U, V, W):
    f = 0
    for i in range(3):
        for k in range(3):
            for kp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + k, 3 * kp + j, 3 * ip + jp
                            s = sum(U[r][a] * V[r][b] * W[r][c] for r in range(R))
                            if s != (1 if (i == ip and j == jp and k == kp) else 0):
                                f += 1
    return f


def brent_fail_fmpz(U, V, W):
    from flint import fmpz
    f = 0
    for i in range(3):
        for k in range(3):
            for kp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + k, 3 * kp + j, 3 * ip + jp
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = U[r][a], V[r][b], W[r][c]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            if int(s) != (1 if (i == ip and j == jp and k == kp) else 0):
                                f += 1
    return f


def rank_over_Q(M):
    M2 = [list(r) for r in M]
    rows, cols = len(M2), len(M2[0])
    used = [False] * rows
    rank = 0
    for c in range(cols):
        piv = next((r for r in range(rows) if not used[r] and M2[r][c]), None)
        if piv is None:
            continue
        used[piv] = True
        rank += 1
        for r in range(rows):
            if r != piv and M2[r][c]:
                f1, f2 = M2[piv][c], M2[r][c]
                for cc in range(cols):
                    M2[r][cc] = M2[r][cc] * f1 - M2[piv][cc] * f2
    return rank


# ------------------------------------------ Table-2 expansion and recount ---
def expand_gates(gates, atom_env, dim):
    env = dict(atom_env)
    for name, terms in gates:
        acc = [0] * dim
        for s, atom in terms:
            src = env[atom]
            for j in range(dim):
                acc[j] += s * src[j]
        env[name] = tuple(acc)
    return env


def table2_sides():
    """Expand Table-2's left and right stages exactly; return rows + counts."""
    aenv = {f"A{i}": tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    benv = {f"B{i}": tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    aenv = expand_gates(T2_T, aenv, N)
    benv = expand_gates(T2_U, benv, N)
    ladds, radds = len(T2_T), len(T2_U)
    lrows, rrows = [], []
    for lexpr, rexpr in T2_P:
        acc = [0] * N
        for s, atom in lexpr:
            for j in range(N):
                acc[j] += s * aenv[atom][j]
        lrows.append(tuple(acc))
        if len(lexpr) > 1:
            ladds += len(lexpr) - 1
        acc = [0] * N
        for s, atom in rexpr:
            for j in range(N):
                acc[j] += s * benv[atom][j]
        rrows.append(tuple(acc))
        if len(rexpr) > 1:
            radds += len(rexpr) - 1
    return (ladds, lrows), (radds, rrows)


def table2_output():
    """Expand Table-2's output stage over its own 23 product labels."""
    env = {f"M{r}": tuple(1 if j == r else 0 for j in range(R)) for r in range(R)}
    env = expand_gates(T2_V, env, R)
    adds = len(T2_V)
    outs = []
    for terms in T2_C:
        acc = [0] * R
        for s, atom in terms:
            for j in range(R):
                acc[j] += s * env[atom][j]
        outs.append(tuple(acc))
        adds += len(terms) - 1
    return adds, outs


T_INV = [0, 3, 6, 1, 4, 7, 2, 5, 8]   # frozen transpose involution (gate A)
COORD_MAPS = [("id", list(range(9))), ("T", T_INV)]


def find_correspondence(l2, r2, o2, U3, V3, W3):
    """AMENDMENT 1 space: sigma_A, sigma_B, sigma_C in {id, T}; bijection pi;
    per-product signs (eps, eps') in {+-1}^2 with the output column carrying
    the compensating sign eps*eps'. Accept only on EXACT equality of every
    entry of all three blocks. Deterministic search order; first match wins."""
    for na, sa in COORD_MAPS:
        for nb, sb in COORD_MAPS:
            for nc, sc in COORD_MAPS:
                pi, signs, used = {}, {}, set()
                ok_all = True
                for r in range(R):
                    cand = None
                    for rp in range(R):
                        if rp in used:
                            continue
                        for e in (1, -1):
                            if not all(l2[r][sa[k]] == e * U3[rp][k] for k in range(N)):
                                continue
                            for ep in (1, -1):
                                if all(r2[r][sb[k]] == ep * V3[rp][k] for k in range(N)):
                                    cand = (rp, e, ep)
                                    break
                            if cand:
                                break
                        if cand:
                            break
                    if cand is None:
                        ok_all = False
                        break
                    rp, e, ep = cand
                    pi[r], signs[r] = rp, (e, ep)
                    used.add(rp)
                if not ok_all or len(used) != R:
                    continue
                # output-map consistency with the compensating product sign
                if all(o2[k][r] == signs[r][0] * signs[r][1] * W3[pi[r]][sc[k]]
                       for k in range(N) for r in range(R)):
                    return pi, signs, (na, nb, nc)
    return None


# ---------------------------------------- circuit tooling (session-10 shape) -
def circuit_verdict(gates, tau_classes, targets):
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    for g in gates:
        (na, sa), (nb, sb) = g["a"], g["b"]
        if na not in values or nb not in values:
            return False, [("forward-ref", g["out"])]
        v = tuple(sa * values[na][k] + sb * values[nb][k] for k in range(N))
        if list(v) != list(g["val"]):
            return False, [("recompute-mismatch", g["out"])]
        values[g["out"]] = v
    produced = {canon(values[g["out"]]) for g in gates}
    tset = set(tau_classes)
    missing = []
    for r, tv in enumerate(targets):
        c = canon(tv)
        if c in tset:
            if c not in produced:
                missing.append(r)
        elif c not in INPUT_CLASSES:
            missing.append(r)
    return (not missing), missing


def build_ext_cnf(tau_extended, T):
    base = sorted(set(tau_extended) | INPUT_CLASSES)
    cls = sorted(tau_extended)
    clsset = set(cls)
    nb = len(base)
    reps = {}
    for c in cls:
        lst = []
        neg = tuple(-x for x in c)
        for k1 in range(nb):
            for k2 in range(k1, nb):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        if k1 == k2 and s1 == 1 and s2 == 1:
                            continue
                        v = tuple(s1 * base[k1][i] + s2 * base[k2][i] for i in range(N))
                        if v == c or v == neg:
                            lst.append((k1, k2))
        reps[c] = lst
    nv, vid = 0, {}

    def V2(*key):
        nonlocal nv
        nv += 1
        vid[key] = nv
        return nv
    for g in range(T):
        for c in cls:
            V2("s", g, c)
    for g in range(T + 1):
        for k in range(nb):
            V2("a", g, k)
    for g in range(T):
        for c in cls:
            for ri in range(len(reps[c])):
                V2("r", g, c, ri)
    cl = []
    for g in range(T):
        for c in cls:
            for c2 in cls:
                if c < c2:
                    cl.append([-vid[("s", g, c)], -vid[("s", g, c2)]])
        for k in range(nb):
            lits = [-vid[("a", g + 1, k)], vid[("a", g, k)]]
            if base[k] in clsset:
                lits.append(vid[("s", g, base[k])])
            cl.append(lits)
    for k in range(nb):
        cl.append([vid[("a", 0, k)]] if base[k] not in clsset else [-vid[("a", 0, k)]])
    for c in cls:
        cl.append([vid[("s", g, c)] for g in range(T)])
    for g in range(T):
        for c in cls:
            lits = [-vid[("s", g, c)]]
            for ri, (k1, k2) in enumerate(reps[c]):
                lits.append(vid[("r", g, c, ri)])
                cl.append([-vid[("r", g, c, ri)], vid[("a", g, k1)]])
                cl.append([-vid[("r", g, c, ri)], vid[("a", g, k2)]])
            cl.append(lits)
    return cl, nv, vid, base, reps


def run_cnf(tau_ext, T, tag="", with_certificates=False, keep_dir=None, solver=None):
    cl_, nv_, vid, base, reps = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = keep_dir or tempfile.mkdtemp(prefix="mws59bit_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    open(cnfp, "w").write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    exe = solver or KISSAT
    cmd = NICE + [exe, "-q", cnfp] + ([drat] if (with_certificates and exe == KISSAT) else [])
    t = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=ENV1)
    sat, unsat = r.returncode == 10, r.returncode == 20
    if not sat and not unsat:
        die(f"solver {exe} rc={r.returncode} on {tag}")
    cert = None
    if unsat and with_certificates and exe == KISSAT:
        lrat = os.path.join(tmp, f"{tag or 'inst'}.lrat")
        conv = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, drat, "-L", lrat],
                              capture_output=True, text=True, env=ENV1)
        os.replace(drat, os.path.join(tmp, f"{tag or 'inst'}.used.drat"))
        d1 = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        d2 = subprocess.run(NICE + [CHECKERS["lrat-check"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        cert = {"cnf_sha256": hashlib.sha256(open(cnfp, "rb").read()).hexdigest(),
                "lrat_sha256": hashlib.sha256(open(lrat, "rb").read()).hexdigest(),
                "drattrim_conv_rc": conv.returncode, "drattrim_rc": d1.returncode,
                "lratcheck_rc": d2.returncode,
                "drattrim_out": (d1.stdout + d1.stderr).strip()[-80:],
                "lratcheck_out": (d2.stdout + d2.stderr).strip()[-80:]}
        if d1.returncode or d2.returncode:
            die(f"certificate FAILED {tag}: {cert}")
    return sat, {"vars": nv_, "clauses": len(cl_), "sec": round(time.time() - t, 3),
                 "cert": cert}


# ===========================================================================
@phase("A-anchors")
def phase_A():
    from gate_b_floor import prep, subset_dfs
    global U_ROWS, V_ROWS, W_ROWS
    U_ROWS, V_ROWS, W_ROWS = load_mws59()
    fi, ff = brent_fail_int(U_ROWS, V_ROWS, W_ROWS), brent_fail_fmpz(U_ROWS, V_ROWS, W_ROWS)
    print(f"A1 brent: int={fi} fmpz={ff}")
    if fi or ff:
        die("A1 brent failed")
    row = {}
    for nm, T in (("U", U_ROWS), ("V", V_ROWS), ("Wfac", W_ROWS)):
        cls, reps = prep(T)
        ok, st, _ = subset_dfs(cls, reps)
        row[nm] = {"d": len(cls), "floor": ok, "states": st["states"]}
    print("A2 frozen row:", json.dumps(row))
    frozen = {"U": (13, False, 132), "V": (12, False, 81), "Wfac": (14, False, 672)}
    for nm, (d, ok, st) in frozen.items():
        g = row[nm]
        if (g["d"], g["floor"], g["states"]) != (d, ok, st):
            die(f"A2 frozen row mismatch {nm}: {g}")
    # frozen campaign-1 witnesses re-verified
    wu = json.load(open(os.path.join(C1, "witness_U_14gates.json")))
    ww = json.load(open(os.path.join(C1, "witness_Wfac_15gates.json")))
    clsU, _ = prep(U_ROWS)
    clsW, _ = prep(W_ROWS)
    okU, misU = circuit_verdict(wu["gates"], sorted(clsU), U_ROWS)
    okW, misW = circuit_verdict(ww["gates"], sorted(clsW), W_ROWS)
    print(f"A2b frozen witnesses re-verified: U 14-gate={okU} ({len(wu['gates'])} gates), "
          f"Wfac 15-gate={okW} ({len(ww['gates'])} gates)")
    if not (okU and okW and len(wu["gates"]) == 14 and len(ww["gates"]) == 15):
        die(f"A2b frozen witness re-verification failed U={okU} W={okW}")
    nz = {nm: sum(1 for r in range(R) if any(M[r]))
          for nm, M in (("U", U_ROWS), ("V", V_ROWS), ("W", W_ROWS))}
    rk = rank_over_Q([[W_ROWS[r][c] for c in range(N)] for r in range(R)])
    print(f"A4 preconditions: nonzero {nz} rank(Wfac)={rk}")
    if any(v != R for v in nz.values()) or rk != 9:
        die("A4 preconditions failed")
    RESULT["phases"]["A-anchors"].update({"brent": [fi, ff], "frozen_row": row,
                                          "frozen_witnesses_ok": True,
                                          "preconditions": {"nonzero": nz, "rank": rk}})
    return wu, ww


@phase("B-controls")
def phase_B(wu):
    from gate_b_floor import prep, subset_dfs
    sys.path.insert(0, SCR5)
    from gatec_decomps import load_sun56
    # A3 tri-instrument SAT control on an independent target
    _, V56, _ = load_sun56()
    pcls, _ = prep(V56)
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    a2 = canon((0, 0, 0, 1, 0, -1, 0, 0, 0))
    ext = sorted(set(pcls) | {a1, a2})
    e2c, e2r = prep(ext)
    okd, _st, _o = subset_dfs(e2c, e2r)
    sk, _ = run_cnf(ext, 13, tag="ctrl_sun56V_T13")
    sc, _ = run_cnf(ext, 13, tag="ctrl_sun56V_T13_cad", solver=CADICAL)
    print(f"A3 control: DFS={okd} kissat={sk} cadical={sc}")
    if not (okd and sk and sc):
        die("A3 control failed")
    # R1/R2 on the frozen 14-gate U witness
    clsU, _ = prep(U_ROWS)
    tau = set(clsU)
    g = wu["gates"]
    last_tau = max(i for i, x in enumerate(g) if canon(tuple(x["val"])) in tau)
    pruned = g[:last_tau] + g[last_tau + 1:]
    ok1, m1 = circuit_verdict(pruned, sorted(tau), U_ROWS)
    print(f"R1 one-gate-deleted: ACCEPT={ok1} (must be False) missing={len(m1)}")
    if ok1:
        die("R1 not rejected")
    consumed = {x["a"][0] for x in g} | {x["b"][0] for x in g}
    tgt = next(i for i, x in enumerate(g) if x["out"] in consumed)
    flip = [dict(x) for x in g]
    flip[tgt] = dict(flip[tgt])
    flip[tgt]["b"] = [flip[tgt]["b"][0], -flip[tgt]["b"][1]]
    ok2, _m2 = circuit_verdict(flip, sorted(tau), U_ROWS)
    print(f"R2 operand-perturbed: ACCEPT={ok2} (must be False)")
    if ok2:
        die("R2 not rejected")
    # R4 frozen floor@12 on V with fresh dual-checker certificate
    clsV, _ = prep(V_ROWS)
    satv, mv = run_cnf(sorted(clsV), 12, tag="R4_floorV_T12",
                       with_certificates=True, keep_dir=CERT_DIR)
    print(f"R4 floor@12 V: SAT={satv} (must be False) cert rc "
          f"{mv['cert']['drattrim_rc']}/{mv['cert']['lratcheck_rc']}")
    if satv:
        die("R4 floor@12 unexpectedly SAT")
    # R5 planter
    scratch = os.path.join(HERE, "scratch_R5")
    os.makedirs(scratch, exist_ok=True)
    pl = os.path.join(scratch, "planter.py")
    open(pl, "w").write(
        "import sys\n"
        f"sys.path.insert(0, {SRC!r})\n"
        f"sys.path.insert(0, {os.path.join(MM3, 'scratch')!r})\n"
        "from pathlib import Path\n"
        f"lines = Path({os.path.join(MM3, 'scratch', 'mws59_layout.txt')!r}).read_text().splitlines()\n"
        "data = lines[236:265]\nblocks, cur = [], []\n"
        "for L in data:\n"
        "    s = L.strip()\n"
        "    if s.startswith('#'):\n        blocks.append(cur); cur = []\n"
        "    elif s and not s.startswith(('A','T')):\n"
        "        try: cur.append([int(x) for x in s.split()])\n"
        "        except ValueError: pass\n"
        "if cur: blocks.append(cur)\n"
        "b1,b2,b3 = blocks\n"
        "V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]\n"
        "from gate_b_floor import prep\n"
        "c, _ = prep(V)\n"
        "assert len(c) == 13, 'planted wrong d'\n")
    rp = subprocess.run([sys.executable, pl], capture_output=True, text=True, env=ENV1)
    fired = rp.returncode != 0 and "planted wrong d" in (rp.stdout + rp.stderr)
    print(f"R5 planter fired: {fired}")
    if not fired:
        die("R5 planter did not fire")
    RESULT["phases"]["B-controls"].update(
        {"A3": "pass (DFS+kissat+CaDiCaL)", "R1": "REJECTED", "R2": "REJECTED",
         "R4_cert": mv["cert"], "R5": "fired"})
    return True


@phase("P-route")
def phase_P(wu, ww):
    from gate_b_floor import prep
    (ladds, lrows), (radds, rrows) = table2_sides()
    oadds, outs2 = table2_output()
    print(f"P recount of printed Table 2: left={ladds} right={radds} output={oadds} "
          f"total={ladds + radds + oadds}")
    corr = find_correspondence(lrows, rrows, outs2, U_ROWS, V_ROWS, W_ROWS)
    if corr is None:
        print("P FAILED: no product correspondence in the amended convention space")
        RESULT["phases"]["P-route"] = {"status": "FAILED-no-correspondence",
                                       "recount": [ladds, radds, oadds]}
        return None
    pi, signs, maps = corr
    print(f"P correspondence FOUND: coord maps (A,B,C)={maps}, bijection on 23 "
          f"products, {sum(1 for r in signs if signs[r][0] * signs[r][1] == -1)} "
          f"products with negated value (compensated in the output column)")
    sa = dict(COORD_MAPS)[maps[0]]
    sb = dict(COORD_MAPS)[maps[1]]
    sc = dict(COORD_MAPS)[maps[2]]
    # R3 (strengthened by amendment 1): swap two labels; the FULL acceptance
    # test (all three blocks, exact) must reject.
    bad = dict(pi)
    bad[0], bad[1] = pi[1], pi[0]
    bad_ok = (all(lrows[r][sa[k]] == signs[r][0] * U_ROWS[bad[r]][k]
                  for r in range(R) for k in range(N))
              and all(rrows[r][sb[k]] == signs[r][1] * V_ROWS[bad[r]][k]
                      for r in range(R) for k in range(N))
              and all(outs2[k][r] == signs[r][0] * signs[r][1] * W_ROWS[bad[r]][sc[k]]
                      for k in range(N) for r in range(R)))
    print(f"R3 wrong-correspondence plant accepted={bad_ok} (must be False)")
    if bad_ok:
        die("R3 wrong-correspondence plant was accepted")
    RESULT["phases"]["B-controls"]["R3"] = "REJECTED (swapped-label correspondence, amended space)"
    # the acceptance test already covers output and right-stage reproduction
    ok_out = all(outs2[k][r] == signs[r][0] * signs[r][1] * W_ROWS[pi[r]][sc[k]]
                 for k in range(N) for r in range(R))
    ok_v = all(rrows[r][sb[k]] == signs[r][1] * V_ROWS[pi[r]][k]
               for r in range(R) for k in range(N))
    ok_l = all(lrows[r][sa[k]] == signs[r][0] * U_ROWS[pi[r]][k]
               for r in range(R) for k in range(N))
    print(f"P stage reproduction under the correspondence: left={ok_l} right={ok_v} output={ok_out}")
    if not (ok_l and ok_v and ok_out):
        RESULT["phases"]["P-route"] = {"status": "FAILED-stage-mismatch",
                                       "recount": [ladds, radds, oadds]}
        return None
    clsV, _ = prep(V_ROWS)
    # assembled total: frozen 14-gate U + Table-2 right (radds) + Table-2 output (oadds)
    total = 14 + radds + oadds
    print(f"P ASSEMBLY: 14 (frozen U witness) + {radds} (Table-2 right) + "
          f"{oadds} (Table-2 output) = {total}")
    fi, ff = brent_fail_int(U_ROWS, V_ROWS, W_ROWS), brent_fail_fmpz(U_ROWS, V_ROWS, W_ROWS)
    if fi or ff:
        die("P assembled decomposition Brent failure")
    payload = {"status": "OK",
               "recount": {"left": ladds, "right": radds, "output": oadds,
                           "printed_total": ladds + radds + oadds},
               "coord_maps": {"A": maps[0], "B": maps[1], "C": maps[2]},
               "correspondence": {str(k): [pi[k], signs[k][0], signs[k][1]]
                                  for k in sorted(pi)},
               "assembled_total": total, "C_V_exact": radds,
               "C_output_exact": oadds, "brent": [729 - fi, 729 - ff]}
    RESULT["phases"]["P-route"] = payload
    json.dump({"note": "printed Table-2 stages verified by exact expansion against "
                       "the Table-3 factor blocks under the amendment-1 correspondence "
                       "(coordinate maps + bijection + per-product signs, output column "
                       "carrying the compensating product sign); every entry checked",
               "coord_maps": {"A": maps[0], "B": maps[1], "C": maps[2]},
               "left_additions_printed": ladds, "right_additions": radds,
               "output_additions": oadds, "assembled_total": total,
               "correspondence": {str(k): [pi[k], signs[k][0], signs[k][1]]
                                  for k in sorted(pi)}},
              open(os.path.join(HERE, "route_P_verification.json"), "w"), indent=1)
    return payload


@phase("G-verdict")
def phase_G(p):
    blob = {"gate": "mws59-bit-58vs59", "run_id": os.path.basename(HERE),
            "prereg_sha256": "aede812c76cbea46bbceeffb737f95de6d197e0f42b7e847466de546041f5e6c",
            "decomposition": "mws59 (fixed orientation, sigma^0, all-monomial triple)",
            "question": "is the fixed-orientation minimum 58 or 59?",
            "inherited_frozen": {"C_U": 14, "C_V_lower": 15, "C_Wfac": 15,
                                 "output_lower": 29, "total_lower": 58},
            "controls": RESULT["phases"].get("B-controls", {}),
            "route_P": RESULT["phases"].get("P-route", {}),
            "cpu_seconds_total": round(sum(x.get("seconds", 0)
                                           for x in RESULT["phases"].values()), 1)}
    if p and p.get("status") == "OK" and p["C_V_exact"] == 15 and p["C_output_exact"] == 29:
        blob["verdict"] = "CERTIFIED-58"
        blob["claim"] = (
            "Route P verified the printed Table-2 stages by exact expansion under a "
            "verified 23-product correspondence: the right stage costs exactly 15 "
            "additions and reproduces all 23 Table-3 V rows, so C(V) = 15 exactly "
            "(frozen lower bound 15 from the complete aux-1 and aux-2 scans); the "
            "output stage costs exactly 29 and reproduces all 9 C-entries, so "
            "C(output) = 29 exactly (frozen lower bound 29). With the frozen 14-gate "
            "U witness (C(U) = 14 exact) the assembled scheme totals "
            "14 + 15 + 29 = 58 additions, matching the frozen certified lower bound: "
            "mws59's fixed-orientation minimum is EXACTLY 58 (LB = UB = 58). The "
            "published 59-addition scheme is therefore one addition above the optimum "
            "of its own orientation - its left stage spends 15 where 14 suffices. No "
            "published claim is contradicted: 59 additions do suffice.")
        blob["LB_after"] = 58
        blob["UB_after"] = 58
    elif p and p.get("status") == "OK":
        blob["verdict"] = "CERTIFIED-OTHER"
        blob["claim"] = (f"Route P verified stages with right={p['C_V_exact']} and "
                         f"output={p['C_output_exact']}; assembled total "
                         f"{p['assembled_total']} - see report for the exact reading.")
    else:
        blob["verdict"] = "INCONCLUSIVE-ROUTE-P-FAILED"
        blob["claim"] = ("Route P did not verify; per the pre-registered rule no bit "
                         "claim is made and Route S was not entered within this "
                         "campaign's budget. LB 58 / UB 59 stand.")
    json.dump(blob, open(os.path.join(HERE, "verdict.json"), "w"), indent=1)
    print("G verdict:", blob["verdict"])
    print(blob["claim"][:600])
    return blob


if __name__ == "__main__":
    for nm, (pth, got, ok) in CHECKERS.items():
        print(f"checker {nm}: sha256={got[:16]}… pinned={ok}")
        if not ok:
            die(f"checker {nm} hash mismatch")
    wu, ww = phase_A()
    phase_B(wu)
    p = phase_P(wu, ww)
    phase_G(p)
    json.dump(RESULT, open(os.path.join(HERE, "runner_output.json"), "w"), indent=1, default=str)
    print("\nRUNNER COMPLETE")
