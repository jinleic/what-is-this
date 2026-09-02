#!/usr/bin/env python3
"""stapleton60-gap-total60 — primary runner (pre-registered; phase 2).

Implements pre_statement.md (sha c0203bb0…) exactly:
  A anchors + printed-SLP recount -> B controls (A1-A4 / R1-R4) ->
  C ScanU (T=15) -> D ScanV (T=15) -> E ScanW (T=14) -> F certificates ->
  G assembly + 729 Brent -> H verdict.
Dual instrument on every instance (frozen subset-DFS + fresh slot-availability
CNF via kissat 4.0.4); disagreement aborts. Written AFTER campaign init.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
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
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(MM3, "scratch"))

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]
ENV1 = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}

N = 9
R = 23
CERT_DIR = os.path.join(HERE, "certificates")
os.makedirs(CERT_DIR, exist_ok=True)

RESULT: dict = {"phases": {}, "defects": []}

PINS = {
    "U":    {"d": 14, "tau": "3b1ec60de0117cc24f2b5e7100e2deac5f72320afd39743508e3b37974528c72",
             "au_n": 428, "au": "9cb18bc76cc897e5c1de0995b37008418df062deebd034f7f4bfeeded681eb68",
             "T": 15},
    "V":    {"d": 14, "tau": "3d8092f303264f77cf4e83a81cf97a63fc0a5675c9c4b4c78d10821baa3c0dc4",
             "au_n": 426, "au": "b01988694358e18db31feb9cd3c782ca5614d1f4bd09bd6720e8fc00cf93256d",
             "T": 15},
    "Wfac": {"d": 13, "tau": "c7c79f0aa65f0e3242d243023cd1f6ae54e42a7fcd11040c91c887a9dfa6a7c3",
             "au_n": 398, "au": "29f85ab1084e9799c0f864ca903c2c83181bc7d0abd4dca25ef092d3558f5693",
             "T": 14},
}

CHECKERS = {}


def die(msg):
    RESULT["fatal"] = msg
    print("FATAL:", msg, flush=True)
    with open(os.path.join(HERE, "verdict_partial.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    sys.exit(1)


for _nm, _want in (("drat-trim", "111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4"),
                   ("lrat-check", "b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07")):
    _p = os.path.join(HERE, "tools", _nm)
    _got = hashlib.sha256(open(_p, "rb").read()).hexdigest() if os.path.isfile(_p) else "MISSING"
    CHECKERS[_nm] = (_p, _got, _got == _want)


def phase(name):
    def deco(fn):
        def wrapped(*a, **k):
            t0 = time.time()
            print(f"== phase {name} ==", flush=True)
            RESULT["phases"].setdefault(name, {})
            out = fn(*a, **k)
            RESULT["phases"][name]["seconds"] = round(time.time() - t0, 2)
            with open(os.path.join(HERE, "verdict_partial.json"), "w") as f:
                json.dump(RESULT, f, indent=1, default=str)
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


# ------------------------------------------------------------------ anchors
def brent_failures_int(U, V, W):
    fails = 0
    for i in range(3):
        for k in range(3):
            for kp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + k, 3 * kp + j, 3 * ip + jp
                            s = sum(U[r][a] * V[r][b] * W[r][c] for r in range(R))
                            if s != (1 if (i == ip and j == jp and k == kp) else 0):
                                fails += 1
    return fails


def brent_failures_fmpz(U, V, W):
    from flint import fmpz
    fails = 0
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
                                fails += 1
    return fails


def rank_over_Q(M):
    M2 = [list(row) for row in M]
    rows, cols = len(M2), len(M2[0])
    used = [False] * rows
    rank = 0
    for c in range(cols):
        piv = None
        for r in range(rows):
            if not used[r] and M2[r][c]:
                piv = r
                break
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


# ------------------------------------- printed-SLP expansion and recount ----
def expand_printed_side(gates, prefix, dim):
    """Expand t/u gates over the 9 input coordinates. Returns env name->vec."""
    env = {f"{prefix}{i}": tuple(1 if j == i else 0 for j in range(dim))
           for i in range(dim)}
    for name, terms in gates:
        acc = [0] * dim
        for s, atom in terms:
            src = env[atom]
            for j in range(dim):
                acc[j] += s * src[j]
        env[name] = tuple(acc)
    return env


def recount_printed_side(gates, prefix, operand_exprs):
    """Exact recount of one input-side stage: explicit gates + inline operand
    additions. Returns (adds, produced rows, distinct value classes)."""
    env = expand_printed_side(gates, prefix, N)
    adds = len(gates)
    rows = []
    for expr in operand_exprs:
        acc = [0] * N
        for s, atom in expr:
            src = env[atom]
            for j in range(N):
                acc[j] += s * src[j]
        rows.append(tuple(acc))
        if len(expr) > 1:
            adds += len(expr) - 1
    return adds, rows, env


def recount_printed_output(v_gates, outputs):
    """Exact recount of the output stage over the 23 products."""
    env = {f"M{r}": tuple(1 if j == r else 0 for j in range(R)) for r in range(R)}
    adds = len(v_gates)
    for name, terms in v_gates:
        acc = [0] * R
        for s, atom in terms:
            src = env[atom]
            for j in range(R):
                acc[j] += s * src[j]
        env[name] = tuple(acc)
    outs = []
    for terms in outputs:
        acc = [0] * R
        for s, atom in terms:
            src = env[atom]
            for j in range(R):
                acc[j] += s * src[j]
        outs.append(tuple(acc))
        adds += len(terms) - 1
    return adds, outs


# --------------------------------------------- fresh CNF (instrument 2) -----
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
    nv = 0
    vid = {}

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
    clauses = []
    for g in range(T):
        for c in cls:
            for c2 in cls:
                if c < c2:
                    clauses.append([-vid[("s", g, c)], -vid[("s", g, c2)]])
        for k in range(nb):
            lits = [-vid[("a", g + 1, k)], vid[("a", g, k)]]
            if base[k] in clsset:
                lits.append(vid[("s", g, base[k])])
            clauses.append(lits)
    for k in range(nb):
        if base[k] not in clsset:
            clauses.append([vid[("a", 0, k)]])
        else:
            clauses.append([-vid[("a", 0, k)]])
    for c in cls:
        clauses.append([vid[("s", g, c)] for g in range(T)])
    for g in range(T):
        for c in cls:
            lits = [-vid[("s", g, c)]]
            for ri, (k1, k2) in enumerate(reps[c]):
                lits.append(vid[("r", g, c, ri)])
                clauses.append([-vid[("r", g, c, ri)], vid[("a", g, k1)]])
                clauses.append([-vid[("r", g, c, ri)], vid[("a", g, k2)]])
            clauses.append(lits)
    return clauses, nv, vid, base, reps


def run_cnf(tau_ext, T, tag="", with_certificates=False, keep_dir=None,
            solver=None):
    cl_, nv_, vid, base, reps = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = keep_dir or tempfile.mkdtemp(prefix="stap60gap_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    with open(cnfp, "w") as f:
        f.write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    exe = solver or KISSAT
    cmd = NICE + [exe, "-q", cnfp] + ([drat] if (with_certificates and exe == KISSAT) else [])
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=ENV1)
    sat, unsat = (r.returncode == 10), (r.returncode == 20)
    if not sat and not unsat:
        die(f"solver {exe} returned {r.returncode} on {tag}")
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
                "drattrim_conv_rc": conv.returncode,
                "drattrim_rc": d1.returncode,
                "drattrim_out": (d1.stdout + d1.stderr).strip()[-90:],
                "lratcheck_rc": d2.returncode,
                "lratcheck_out": (d2.stdout + d2.stderr).strip()[-90:]}
        if d1.returncode != 0 or d2.returncode != 0:
            die(f"certificate verification FAILED for {tag}: {cert}")
    return sat, {"vars": nv_, "clauses": len(cl_), "sec": round(time.time() - t0, 3),
                 "cert": cert, "path": cnfp}


# ------------------------------------------------------- universes + verify --
def aux_universe(tau_classes):
    forms = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)] + sorted(tau_classes)
    ts = set(tau_classes)
    u = set()
    for i in range(len(forms)):
        for j in range(i + 1):
            for si in (1, -1):
                for sj in (1, -1):
                    if i == j and si == 1 and sj == 1:
                        continue
                    c = canon(tuple(si * forms[i][k] + sj * forms[j][k] for k in range(N)))
                    if c != (0,) * N and c not in INPUT_CLASSES and c not in ts:
                        u.add(c)
    return sorted(u)


def csv_hash(rows):
    return hashlib.sha256(("\n".join(",".join(str(x) for x in r) for r in rows)).encode()).hexdigest()


def replay_schedule(order):
    """Replay a class-creation order as an explicit gate list over exact
    vectors. Returns the gate list or None."""
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    gates = []
    for idx, c in enumerate(order):
        names = sorted(values.keys())
        found = None
        for i1 in range(len(names)):
            for i2 in range(i1, len(names)):
                va, vb = values[names[i1]], values[names[i2]]
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        if names[i1] == names[i2] and s1 == 1 and s2 == 1:
                            continue
                        v = tuple(s1 * va[k] + s2 * vb[k] for k in range(N))
                        if canon(v) == c:
                            found = ((names[i1], s1), (names[i2], s2), v)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if not found:
            return None
        (n1, s1), (n2, s2), v = found
        out = 100 + idx
        values[out] = v
        gates.append({"out": out, "a": [n1, s1], "b": [n2, s2],
                      "cls": list(c), "val": list(v)})
    return gates


def circuit_verdict(gates, tau_classes, targets):
    """Exact-expansion verifier, independent of the producing search: recompute
    every gate from its operand structure and require every target row to be
    ± a produced value (or ± a free input). Returns (ok, missing)."""
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


# ===========================================================================
@phase("A-anchors")
def phase_A():
    from gate_b_floor import prep, subset_dfs
    sys.path.insert(0, SCR5)
    from gatec_decomps import load_stapleton60, META
    import stapleton60_data as sd

    global U_ROWS, V_ROWS, W_ROWS
    U, V, W = load_stapleton60()
    U_ROWS, V_ROWS, W_ROWS = U, V, W

    f_int, f_fmpz = brent_failures_int(U, V, W), brent_failures_fmpz(U, V, W)
    print(f"A1 brent: int={f_int} fmpz={f_fmpz}")
    if f_int or f_fmpz:
        die(f"A1 brent failed: {f_int}/{f_fmpz}")
    if not all(x in (-1, 0, 1) for blk in (U, V, W) for row in blk for x in row):
        die("A1 factors not ternary")

    # A2 printed-SLP recount by exact expansion
    left_exprs = [p[0] for p in sd.PRODUCTS]
    right_exprs = [p[1] for p in sd.PRODUCTS]
    ladds, lrows, _ = recount_printed_side(sd.T_GATES, "A", left_exprs)
    radds, rrows, _ = recount_printed_side(sd.U_GATES, "B", right_exprs)
    oadds, outs = recount_printed_output(sd.V_GATES, sd.OUTPUTS)
    print(f"A2 printed recount: left={ladds} right={radds} output={oadds} total={ladds+radds+oadds}")
    if (ladds, radds, oadds) != (16, 16, 28):
        die(f"A2 recount mismatch: {(ladds, radds, oadds)}")
    # the printed rows must equal the loader's factor rows up to sign
    lok = all(canon(lrows[r]) == canon(U[r]) for r in range(R))
    rok = all(canon(rrows[r]) == canon(V[r]) for r in range(R))
    ook = all(all(outs[k][r] == W[r][k] or outs[k][r] == -W[r][k] for r in range(R))
              for k in range(N))
    ook_exact = all(all(outs[k][r] == W[r][k] for r in range(R)) for k in range(N))
    print(f"A2 rows reproduced: left={lok} right={rok} output(up to sign)={ook} exact={ook_exact}")
    if not (lok and rok and ook):
        die("A2 printed stages do not reproduce the factor rows")
    RESULT["phases"]["A-anchors"].update(
        {"brent": [f_int, f_fmpz], "printed_split": [ladds, radds, oadds],
         "printed_total": ladds + radds + oadds,
         "rows_reproduced": {"left": lok, "right": rok, "output_up_to_sign": ook,
                             "output_exact": ook_exact}})

    # A1b frozen landscape row
    row = {}
    for nm, T in (("U", U), ("V", V), ("Wfac", W)):
        cls, reps = prep(T)
        ok, stats, _ = subset_dfs(cls, reps)
        row[nm] = {"d": len(cls), "floor": ok, "states": stats["states"]}
    print("A1b frozen row:", json.dumps(row))
    frozen = {"U": (14, False, 1152), "V": (14, False, 1344), "Wfac": (13, False, 460)}
    for nm, (d, ok, st) in frozen.items():
        got = row[nm]
        if (got["d"], got["floor"], got["states"]) != (d, ok, st):
            die(f"A1b frozen row mismatch {nm}: {got} vs {frozen[nm]}")
    if 15 + 15 + 14 + 14 != 58:
        die("A1b arithmetic")
    RESULT["phases"]["A-anchors"]["frozen_row"] = row

    # A4 transposition preconditions
    nz = {nm: sum(1 for r in range(R) if any(M[r])) for nm, M in (("U", U), ("V", V), ("W", W))}
    rk = rank_over_Q([[W[r][c] for c in range(N)] for r in range(R)])
    print(f"A4 preconditions: nonzero rows {nz}, rank(Wfac)={rk}")
    if any(v != R for v in nz.values()) or rk != 9:
        die(f"A4 transposition precondition failed: {nz} rank={rk}")
    RESULT["phases"]["A-anchors"]["preconditions"] = {"nonzero_rows": nz, "rank_Wfac": rk}
    return True


@phase("B-controls")
def phase_B():
    from gate_b_floor import prep, subset_dfs
    from gatec_decomps import load_sun56
    import stapleton60_data as sd

    # A3 extension-positive control on an INDEPENDENT target (sun56-V)
    _, V56, _ = load_sun56()
    pcls, _ = prep(V56)
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    a2 = canon((0, 0, 0, 1, 0, -1, 0, 0, 0))
    ext = sorted(set(pcls) | {a1, a2})
    e2c, e2r = prep(ext)
    ok13, st13, _ = subset_dfs(e2c, e2r)
    sat13, _m = run_cnf(ext, 13, tag="ctrl_sun56V_aux2_T13")
    sat13b, _m2 = run_cnf(ext, 13, tag="ctrl_sun56V_aux2_T13_cadical", solver=CADICAL)
    print(f"A3 control: DFS={ok13} kissat={sat13} cadical={sat13b} (all must be True)")
    if not (ok13 and sat13 and sat13b):
        die("A3 extension-positive control failed / instruments disagree")
    RESULT["phases"]["B-controls"]["A3"] = "pass (DFS + kissat + CaDiCaL)"

    # R1/R2 plants on a synthesized Wfac-side witness shape: use the FIRST
    # feasible aux found on the Wfac side later; here use the printed-left
    # stage as the plant carrier (exactly 16 additions, verified in A2).
    left_exprs = [p[0] for p in sd.PRODUCTS]
    ladds, lrows, env = recount_printed_side(sd.T_GATES, "A", left_exprs)
    tauU = set(prep(U_ROWS)[0])
    # explicit gate list of the printed left stage (t-gates + inline sums)
    gates = []
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    name_to_id = {f"A{i}": i for i in range(N)}
    nxt = 100
    for name, terms in sd.T_GATES:
        acc_id = None
        cur = None
        for s, atom in terms:
            if cur is None:
                cur = tuple(s * x for x in values[name_to_id[atom]])
                acc_id = (name_to_id[atom], s)
                continue
            v = tuple(cur[k] + s * values[name_to_id[atom]][k] for k in range(N))
            gates.append({"out": nxt, "a": [acc_id[0], acc_id[1]],
                          "b": [name_to_id[atom], s], "cls": list(canon(v)),
                          "val": list(v)})
            values[nxt] = v
            cur = v
            acc_id = (nxt, 1)
            nxt += 1
        name_to_id[name] = acc_id[0]
    for expr in left_exprs:
        if len(expr) == 1:
            continue
        cur = None
        acc_id = None
        for s, atom in expr:
            if cur is None:
                cur = tuple(s * x for x in values[name_to_id[atom]])
                acc_id = (name_to_id[atom], s)
                continue
            v = tuple(cur[k] + s * values[name_to_id[atom]][k] for k in range(N))
            gates.append({"out": nxt, "a": [acc_id[0], acc_id[1]],
                          "b": [name_to_id[atom], s], "cls": list(canon(v)),
                          "val": list(v)})
            values[nxt] = v
            cur = v
            acc_id = (nxt, 1)
            nxt += 1
    print(f"B printed-left explicit gate list: {len(gates)} gates (expect 16)")
    if len(gates) != 16:
        die(f"B plant carrier has {len(gates)} gates, expected 16")
    okc, missc = circuit_verdict(gates, sorted(tauU), U_ROWS)
    print(f"B carrier verifies: {okc} missing={missc[:3]}")
    if not okc:
        die("B plant carrier fails exact expansion — cannot host the plants")
    RESULT["phases"]["B-controls"]["carrier_gates"] = len(gates)

    # R1 delete the last tau-creating gate
    tau_pos = [i for i, g in enumerate(gates) if canon(tuple(g["val"])) in tauU]
    pruned = gates[:tau_pos[-1]] + gates[tau_pos[-1] + 1:]
    ok_r1, miss_r1 = circuit_verdict(pruned, sorted(tauU), U_ROWS)
    print(f"R1 one-gate-deleted: verifier ACCEPT={ok_r1} (must be False), missing={len(miss_r1)}")
    if ok_r1:
        die("R1 deletion plant not rejected")

    # R2 operand-sign perturbation on a consumed gate
    consumed = {g["a"][0] for g in gates} | {g["b"][0] for g in gates}
    tgt = next(i for i, g in enumerate(gates) if g["out"] in consumed)
    flipped = [dict(g) for g in gates]
    flipped[tgt] = dict(flipped[tgt])
    flipped[tgt]["b"] = [flipped[tgt]["b"][0], -flipped[tgt]["b"][1]]
    ok_r2, miss_r2 = circuit_verdict(flipped, sorted(tauU), U_ROWS)
    print(f"R2 operand-perturbed: verifier ACCEPT={ok_r2} (must be False)")
    if ok_r2:
        die("R2 perturbation plant not rejected")
    RESULT["phases"]["B-controls"]["R1"] = "REJECTED by exact expansion"
    RESULT["phases"]["B-controls"]["R2"] = "REJECTED by exact expansion"

    # R3 weaker-T negative with fresh dual-checker certificate
    wcls, _ = prep(W_ROWS)
    sat12, m12 = run_cnf(sorted(wcls), 12, tag="R3_tauW_T12",
                         with_certificates=True, keep_dir=CERT_DIR)
    print(f"R3 tau(Wfac) at T=12: SAT={sat12} (must be False), cert rc "
          f"{m12['cert']['drattrim_rc'] if m12['cert'] else 'NA'}/"
          f"{m12['cert']['lratcheck_rc'] if m12['cert'] else 'NA'}")
    if sat12:
        die("R3 T=12 plant unexpectedly SAT")
    RESULT["phases"]["B-controls"]["R3_cert"] = m12["cert"]

    # R4 assertion planter
    scratch = os.path.join(HERE, "scratch_R4")
    os.makedirs(scratch, exist_ok=True)
    planter = os.path.join(scratch, "planter.py")
    with open(planter, "w") as f:
        f.write("import sys\n"
                f"sys.path.insert(0, {SRC!r})\n"
                f"sys.path.insert(0, {SCR5!r})\n"
                f"sys.path.insert(0, {os.path.join(MM3, 'scratch')!r})\n"
                "from gatec_decomps import load_stapleton60\n"
                "from gate_b_floor import prep\n"
                "U, V, W = load_stapleton60()\n"
                "c, _ = prep(U)\n"
                "assert len(c) == 15, 'planted wrong d'\n")
    rp = subprocess.run([sys.executable, planter], capture_output=True, text=True, env=ENV1)
    fired = rp.returncode != 0 and "planted wrong d" in (rp.stdout + rp.stderr)
    print(f"R4 planter fired: {fired} (rc={rp.returncode})")
    if not fired:
        die("R4 planter did not fire")
    RESULT["phases"]["B-controls"]["R4"] = "fired"
    return gates


def run_scan(side, rows, phase_name):
    """Dual-instrument aux-1 scan at T = d+1 over the pinned universe."""
    from gate_b_floor import prep, subset_dfs
    pin = PINS[side]
    cls, _ = prep(rows)
    if len(cls) != pin["d"] or csv_hash(sorted(cls)) != pin["tau"]:
        die(f"{side}: tau mismatch (d={len(cls)}, sha={csv_hash(sorted(cls))})")
    au = aux_universe(cls)
    if len(au) != pin["au_n"] or csv_hash(au) != pin["au"]:
        die(f"{side}: AU mismatch (n={len(au)}, sha={csv_hash(au)})")
    T = pin["T"]
    print(f"{phase_name} universe {side}: d={len(cls)} |AU|={len(au)} T={T} pins OK")
    ck = open(os.path.join(HERE, f"scan_{side}_checkpoint.jsonl"), "a")
    hits = []
    t0 = time.time()
    for idx, a in enumerate(au):
        ext = sorted(set(cls) | {a})
        e2c, e2r = prep(ext)
        ok, stats, order = subset_dfs(e2c, e2r)
        sat, meta = run_cnf(ext, T, tag=f"{side}{idx:03d}")
        if ok != sat:
            die(f"{phase_name} instrument disagreement at {side}[{idx}]={a}: DFS={ok} kissat={sat}")
        ck.write(json.dumps({"i": idx, "aux": list(a), "dfs": ok,
                             "states": stats["states"],
                             "t": round(time.time() - t0, 1)}) + "\n")
        ck.flush()
        os.fsync(ck.fileno())
        if ok:
            hits.append((a, order))
        if (idx + 1) % 100 == 0:
            print(f"{phase_name} {side} {idx+1}/{len(au)} hits={len(hits)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    ck.close()
    print(f"{phase_name} {side} COMPLETE: {len(au)}/{len(au)} agreed, feasible={len(hits)}")
    RESULT["phases"][phase_name].update(
        {"side": side, "instances": len(au), "agreed": len(au),
         "hits": [list(a) for a, _ in hits], "T": T})
    return cls, au, hits


@phase("C-ScanU")
def phase_C():
    return run_scan("U", U_ROWS, "C-ScanU")


@phase("D-ScanV")
def phase_D():
    return run_scan("V", V_ROWS, "D-ScanV")


@phase("E-ScanW")
def phase_E():
    return run_scan("Wfac", W_ROWS, "E-ScanW")


@phase("F-certificates")
def phase_F(scans):
    from gate_b_floor import prep
    out = {}
    for side, rows, T_floor in (("U", U_ROWS, 14), ("V", V_ROWS, 14), ("Wfac", W_ROWS, 13)):
        cls, _ = prep(rows)
        sat, meta = run_cnf(sorted(cls), T_floor, tag=f"floor_{side}_T{T_floor}",
                            with_certificates=True, keep_dir=CERT_DIR)
        print(f"F floor_{side}_T{T_floor}: SAT={sat} (must be False)")
        if sat:
            die(f"F floor_{side} unexpectedly SAT — contradicts the frozen landscape")
        out[f"floor_{side}_T{T_floor}"] = meta["cert"]
    # per-scan: first infeasible + seeded sample of 10 certified
    random.seed(20260902)
    for side, (cls, au, hits) in scans.items():
        hitset = {tuple(a) for a, _ in hits}
        infeas = [a for a in au if tuple(a) not in hitset]
        if not infeas:
            continue
        picks = [infeas[0]] + random.sample(infeas, min(10, len(infeas)))
        seen = set()
        for j, a in enumerate(picks):
            if tuple(a) in seen:
                continue
            seen.add(tuple(a))
            ext = sorted(set(cls) | {a})
            sat, meta = run_cnf(ext, PINS[side]["T"], tag=f"{side}_infeas_{j:02d}",
                                with_certificates=True, keep_dir=CERT_DIR)
            if sat:
                die(f"F {side} sampled instance unexpectedly SAT")
            out[f"{side}_infeas_{j:02d}"] = {"aux": list(a), **(meta["cert"] or {})}
        print(f"F {side}: {len(seen)} infeasible instances certified (first + seeded sample)")
    RESULT["phases"]["F-certificates"] = out
    return out


@phase("G-assembly")
def phase_G(scans, carrier_gates):
    """Assemble the best circuit the scans support and verify it end to end."""
    from gate_b_floor import prep
    import stapleton60_data as sd
    sides = {}
    for side, (cls, au, hits) in scans.items():
        if hits:
            a, order = hits[0]
            gates = replay_schedule(order)
            if gates is None:
                die(f"G replay failed for {side} aux {a}")
            rows = {"U": U_ROWS, "V": V_ROWS, "Wfac": W_ROWS}[side]
            ok, missing = circuit_verdict(gates, sorted(cls), rows)
            if not ok:
                die(f"G {side} witness fails exact expansion: missing={missing[:4]}")
            sides[side] = {"gates": len(gates), "aux": list(a), "verified": True,
                           "circuit": gates}
            print(f"G {side}: {len(gates)}-gate witness verified by exact expansion")
        else:
            sides[side] = {"gates": None, "verified": False,
                           "note": "no admission; value pinned by printed witness + floor"}
    # exact stage values
    cU = 15 if scans["U"][2] else 16
    cV = 15 if scans["V"][2] else 16
    cW = 14 if scans["Wfac"][2] else None
    out_stage = 28
    total = cU + cV + out_stage
    print(f"G stage values: C(U)={cU} C(V)={cV} C(Wfac)={cW} output={out_stage} total={total}")
    # full-scheme Brent (the decomposition is unchanged by circuit choice, so
    # this re-verifies the object the assembled circuit computes)
    f_int = brent_failures_int(U_ROWS, V_ROWS, W_ROWS)
    f_fmpz = brent_failures_fmpz(U_ROWS, V_ROWS, W_ROWS)
    if f_int or f_fmpz:
        die("G assembled scheme Brent failure")
    payload = {"C_U": cU, "C_V": cV, "C_Wfac": cW, "output": out_stage,
               "total": total, "brent": [729 - f_int, 729 - f_fmpz],
               "witnesses": {k: {kk: vv for kk, vv in v.items() if kk != "circuit"}
                             for k, v in sides.items()}}
    for k, v in sides.items():
        if v.get("circuit"):
            with open(os.path.join(HERE, f"witness_{k}_{v['gates']}gates.json"), "w") as f:
                json.dump({"aux": v["aux"], "gates": v["circuit"],
                           "note": f"{v['gates']}-gate {k} circuit replayed from the admitted "
                                   f"aux schedule and verified by exact expansion "
                                   f"independent of the search"}, f, indent=1)
    RESULT["phases"]["G-assembly"] = payload
    return payload


@phase("H-verdict")
def phase_H(scans, asm):
    cU, cV, cW = asm["C_U"], asm["C_V"], asm["C_Wfac"]
    total = asm["total"]
    complete = all(len(scans[s][1]) == PINS[s]["au_n"] for s in scans)
    blob = {
        "gate": "stapleton60-gap-total60",
        "run_id": os.path.basename(HERE),
        "prereg_sha256": "c0203bb05df042570d2f3da40a0d59109f4603618de4f27c86300fd2e494f707",
        "decomposition": "stapleton60 (Stapleton's fixed orientation, sigma^0, all-monomial triple)",
        "question": "decide the fixed-orientation minimum in the frozen gap [58, 60]",
        "reduction": ("total = C(U) + C(V) + C(output); C(output) = 28 exact "
                      "(frozen floor@13 on Wfac => >= 28; printed output stage "
                      "recounted in-run by exact expansion = 28); the gap is "
                      "exactly the two single-aux questions C(U) <= 15?, C(V) <= 15?"),
        "stage_values": {"C_U": cU, "C_V": cV, "C_Wfac": cW, "output": 28},
        "search_space": {s: {"instances": PINS[s]["au_n"], "T": PINS[s]["T"],
                             "tau_sha256": PINS[s]["tau"], "au_sha256": PINS[s]["au"],
                             "agreed": RESULT["phases"][p]["agreed"],
                             "hits": RESULT["phases"][p]["hits"]}
                         for s, p in (("U", "C-ScanU"), ("V", "D-ScanV"), ("Wfac", "E-ScanW"))},
        "instruments": ["subset_dfs (complete census)",
                        "fresh slot-availability CNF + kissat 4.0.4",
                        "CaDiCaL 3.0.1 (control cross-solve)"],
        "certificates": {k: v for k, v in RESULT["phases"]["F-certificates"].items()},
        "controls": RESULT["phases"]["B-controls"],
        "brent": asm["brent"],
        "cpu_seconds_total": round(sum(p.get("seconds", 0) for p in RESULT["phases"].values()), 1),
    }
    if not complete:
        blob["verdict"] = "INCONCLUSIVE-PARTIAL"
        blob["claim"] = "a scan did not complete; only the completed prefix is reported"
    elif cU == 16 and cV == 16:
        blob["verdict"] = "NEGATIVE-PUBLISHED-OPTIMAL"
        blob["claim"] = (
            "C(U) = C(V) = 16 exactly (floors@14 impossible, re-certified; "
            "0 of 428 and 0 of 426 single-aux extensions admitted at T = 15, "
            "both instruments agreeing on every instance; printed 16-gate "
            "stages verified by exact expansion) and output = 28 exactly, so "
            "total = 16 + 16 + 28 = 60 EXACTLY: Stapleton's published "
            "60-addition scheme is optimal for its own fixed orientation, and "
            "the frozen certified LB 58 was not tight. LB = UB = 60.")
    else:
        blob["verdict"] = "CERTIFIED-IMPROVEMENT"
        blob["claim"] = (
            f"C(U) = {cU}, C(V) = {cV} (each exact: floor + admitted-aux "
            f"witness replayed and verified by exact expansion), output = 28 "
            f"exact => total = {total} EXACTLY, improving the published 60 by "
            f"{60 - total}. LB = UB = {total}.")
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump(blob, f, indent=1)
    print("H verdict:", blob["verdict"])
    print(blob["claim"])
    return blob


if __name__ == "__main__":
    for nm, (p, got, ok) in CHECKERS.items():
        print(f"checker {nm}: sha256={got[:16]}… pinned={ok}")
        if not ok:
            die(f"checker {nm} hash mismatch: {got}")
    phase_A()
    carrier = phase_B()
    scans = {}
    scans["U"] = phase_C()
    scans["V"] = phase_D()
    scans["Wfac"] = phase_E()
    phase_F(scans)
    asm = phase_G(scans, carrier)
    phase_H(scans, asm)
    with open(os.path.join(HERE, "runner_output.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    print("\nRUNNER COMPLETE")
