#!/usr/bin/env python3
"""paper55-sigma12-total55 — primary runner (pre-registered; phase 2).

prereg sha 10b3f32f…. Decides the paper55 sigma^1 and sigma^2 orientations,
each frozen at [55, 57], by (i) inheriting the frozen gate-B exact side costs
through the proved fixed-coordinate-relabeling symmetry, (ii) synthesizing the
missing 14-gate Wfac circuit with the aux-1 dual instrument, (iii) building each
orientation's output stage by reverse-mode transposition as an explicit gate
list, and (iv) verifying every stage by exact expansion plus an end-to-end
symbolic expansion of the assembled scheme over all 81 monomials.
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
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(MM3, "scratch"))

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]
ENV1 = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}
N, R = 9, 23
T_INV = [0, 3, 6, 1, 4, 7, 2, 5, 8]
CERT_DIR = os.path.join(HERE, "certificates")
os.makedirs(CERT_DIR, exist_ok=True)
RESULT: dict = {"phases": {}, "defects": []}
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
    _g = hashlib.sha256(open(_p, "rb").read()).hexdigest() if os.path.isfile(_p) else "MISSING"
    CHECKERS[_nm] = (_p, _g, _g == _want)


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


def Tblocks(rows):
    """Per-product 3x3 transpose of 9-vectors (the frozen T of sigma_orbit)."""
    out = []
    for r in rows:
        b = [r[0:3], r[3:6], r[6:9]]
        bt = [list(x) for x in zip(*b)]
        out.append(tuple(bt[0] + bt[1] + bt[2]))
    return out


def brent_fail(U, V, W, exact=False):
    if exact:
        from flint import fmpz
    f = 0
    for i in range(3):
        for k in range(3):
            for kp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + k, 3 * kp + j, 3 * ip + jp
                            if exact:
                                s = fmpz(0)
                                for r in range(R):
                                    ur, vr, wr = U[r][a], V[r][b], W[r][c]
                                    if ur and vr and wr:
                                        s += ur * vr * wr
                                s = int(s)
                            else:
                                s = sum(U[r][a] * V[r][b] * W[r][c] for r in range(R))
                            if s != (1 if (i == ip and j == jp and k == kp) else 0):
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


# ------------------------------------------------- printed-SLP expansion -----
def expand_printed(slp, prefix):
    """Expand a printed input-side SLP; returns (env, gate_count, gate_list)
    where gate_list has explicit two-operand structure over 9 inputs."""
    env = {f"{prefix}{i}": tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    ids = {f"{prefix}{i}": i for i in range(N)}
    gates = []
    nxt = 100
    for name, terms in slp:
        assert len(terms) == 2, f"{name}: expected two terms"
        (s1, a1), (s2, a2) = terms
        v = tuple(s1 * env[a1][k] + s2 * env[a2][k] for k in range(N))
        env[name] = v
        gates.append({"out": nxt, "a": [ids[a1], s1], "b": [ids[a2], s2],
                      "val": list(v), "name": name})
        ids[name] = nxt
        nxt += 1
    return env, len(gates), gates, ids


def eval_circuit(gates, input_perm=None):
    """Recompute every gate from its operand structure over 9 free inputs.
    `input_perm` feeds input wire i with coordinate perm[i] (a FREE relabeling
    of input wires); with perm=None the stored `val` fields are additionally
    asserted for integrity. Returns the value dict or None on inconsistency."""
    if input_perm is None:
        values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    else:
        values = {i: tuple(1 if j == input_perm[i] else 0 for j in range(N))
                  for i in range(N)}
    for g in gates:
        (na, sa), (nb, sb) = g["a"], g["b"]
        if na not in values or nb not in values:
            return None
        v = tuple(sa * values[na][k] + sb * values[nb][k] for k in range(N))
        if input_perm is None and list(v) != list(g["val"]):
            return None
        values[g["out"]] = v
    return values


def circuit_verdict(gates, targets, input_perm=None):
    """Exact-expansion verifier: recompute each gate; every target row must be
    +/- a produced value or +/- a free input."""
    values = eval_circuit(gates, input_perm)
    if values is None:
        return False, [("recompute-inconsistent", None)]
    avail = {canon(v) for v in values.values()}
    missing = [r for r, tv in enumerate(targets) if canon(tv) not in avail]
    return (not missing), missing


# --------------------------------------------- aux-1 instrument (dual) -------
def aux_universe(tau):
    forms = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)] + sorted(tau)
    ts = set(tau)
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


def aux_universe_alt(tau):
    """Independent second enumeration loop (different shape) for the pin."""
    ts = set(tau)
    elems = sorted(INPUT_CLASSES) + sorted(tau)
    u = set()
    for x in elems:
        for y in elems:
            for a_ in (x, tuple(-t for t in x)):
                for b_ in (y, tuple(-t for t in y)):
                    c = canon(tuple(a_[k] + b_[k] for k in range(N)))
                    if c != (0,) * N and c not in INPUT_CLASSES and c not in ts:
                        u.add(c)
    return sorted(u)


def csv_hash(rows):
    return hashlib.sha256(("\n".join(",".join(str(x) for x in r) for r in rows)).encode()).hexdigest()


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
    return cl, nv


def run_cnf(tau_ext, T, tag="", with_certificates=False, keep_dir=None, solver=None):
    cl_, nv_ = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = keep_dir or tempfile.mkdtemp(prefix="p55sig_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    open(cnfp, "w").write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    exe = solver or KISSAT
    cmd = NICE + [exe, "-q", cnfp] + ([drat] if (with_certificates and exe == KISSAT) else [])
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=ENV1)
    sat, unsat = r.returncode == 10, r.returncode == 20
    if not sat and not unsat:
        die(f"solver rc {r.returncode} on {tag}")
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
                "drattrim_rc": d1.returncode, "lratcheck_rc": d2.returncode,
                "drattrim_conv_rc": conv.returncode}
        if d1.returncode or d2.returncode:
            die(f"certificate FAILED {tag}: {cert}")
    return sat, {"vars": nv_, "clauses": len(cl_), "cert": cert}


def replay_schedule(order):
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
        gates.append({"out": out, "a": [n1, s1], "b": [n2, s2], "val": list(v)})
    return gates


# ------------------------------- reverse-mode transposition (construction) ---
def transpose_circuit(gates, targets, input_perm=None):
    """Given an explicit 9->23 circuit (gate list over 9 free inputs, optionally
    with a free input relabeling) whose produced values cover every target row
    up to sign, build the transposed 23->9 output-stage circuit as an EXPLICIT
    gate list and return (adds, outputs_over_products, meta).

    Construction: each target row r is read off some wire with a sign; seed that
    wire's adjoint with the product variable r (signed), then process gates in
    reverse, pushing each gate's adjoint into its two operands with the operand
    signs; the adjoints of the 9 inputs are the 9 outputs. Additions counted
    exactly: every accumulation into an already-nonempty adjoint costs 1."""
    values = eval_circuit(gates, input_perm)
    if values is None:
        return None
    # read-off map: target r -> (wire, sign)
    reads = []
    for r, tv in enumerate(targets):
        pick = None
        for wid, v in values.items():
            if list(v) == list(tv):
                pick = (wid, 1)
                break
            if all(v[k] == -tv[k] for k in range(N)):
                pick = (wid, -1)
                break
        if pick is None:
            return None
        reads.append(pick)
    # adjoints as vectors over the 23 product variables
    adj = {wid: [0] * R for wid in values}
    nonempty = set()
    adds = 0
    for r, (wid, sg) in enumerate(reads):
        if wid in nonempty:
            adds += 1                      # accumulate another product into the wire
        else:
            nonempty.add(wid)
        adj[wid][r] += sg
    for g in reversed(gates):
        out = g["out"]
        if out not in nonempty:
            continue                        # dead wire: contributes nothing
        (na, sa), (nb, sb) = g["a"], g["b"]
        for (nn, ss) in ((na, sa), (nb, sb)):
            if nn in nonempty:
                adds += 1                  # accumulate into an existing adjoint
            else:
                nonempty.add(nn)
            for j in range(R):
                if adj[out][j]:
                    adj[nn][j] += ss * adj[out][j]
    # The 9 transposed outputs are indexed by INPUT WIRES; when the forward
    # circuit was fed under a coordinate relabeling, wire k carries coordinate
    # input_perm[k], so the outputs must be reported in coordinate order.
    # Re-indexing output wires is a free relabeling (defect D2, caught by the
    # pre-registered output-map equality check).
    if input_perm is None:
        outs = [tuple(adj[i]) for i in range(N)]
    else:
        outs = [None] * N
        for k in range(N):
            outs[input_perm[k]] = tuple(adj[k])
    return adds, outs, {"reads": reads}


# ===========================================================================
@phase("A-anchors")
def phase_A():
    from gate_b_floor import prep, subset_dfs
    import tensor_data as td
    sys.path.insert(0, SCR5)
    from gatec_decomps import load_paper55
    U, V, W = load_paper55()
    fi, ff = brent_fail(U, V, W), brent_fail(U, V, W, exact=True)
    print(f"A1 sigma^0 brent: int={fi} fmpz={ff}")
    if fi or ff:
        die("A1 brent failed")
    # printed circuits re-verified + recounted
    envU, nU, gU, _ = expand_printed(td.LEFT_SLP, "A")
    envV, nV, gV, _ = expand_printed(td.RIGHT_SLP, "B")
    okU, misU = circuit_verdict(gU, U)
    okV, misV = circuit_verdict(gV, V)
    print(f"A2 printed left: {nU} gates, covers 23 U rows={okU}; "
          f"printed right: {nV} gates, covers 23 V rows={okV}")
    if not (okU and okV and nU == 13 and nV == 14):
        die(f"A2 printed circuits failed: {nU},{okU},{nV},{okV}")
    # printed output stage: expand over the 23 products, count, check the 9 outputs
    envW = {f"M{r}": tuple(1 if j == r else 0 for j in range(R)) for r in range(R)}
    for name, terms in td.OUTPUT_SLP:
        acc = [0] * R
        for s, atom in terms:
            for j in range(R):
                acc[j] += s * envW[atom][j]
        envW[name] = tuple(acc)
    outs0 = [envW[a] for a in td.C_ALIASES]
    nW = len(td.OUTPUT_SLP)
    ok_out = all(outs0[k][r] == W[r][k] for k in range(N) for r in range(R))
    print(f"A2 printed output: {nW} gates, reproduces W columns exactly={ok_out}")
    if not (ok_out and nW == 28):
        die(f"A2 printed output failed: {nW}, {ok_out}")
    print(f"A2 sigma^0 total recount = {nU} + {nV} + {nW} = {nU + nV + nW}")
    if nU + nV + nW != 55:
        die("A2 sigma^0 total is not 55")
    # A3 frozen rows for the sigma classes + A4 relabeling symmetry
    orbit = {"sigma^0": (U, V, W), "sigma^1": (V, Tblocks(W), Tblocks(U)),
             "sigma^2": (Tblocks(W), U, Tblocks(V))}
    frozen = {"sigma^0": [(12, 33), (13, 116), (13, 66)],
              "sigma^1": [(13, 116), (13, 66), (12, 33)],
              "sigma^2": [(13, 66), (12, 33), (13, 116)]}
    rows = {}
    for nm, (X, Y, Z) in orbit.items():
        got = []
        for M in (X, Y, Z):
            cls, reps = prep(M)
            ok, st, _ = subset_dfs(cls, reps)
            if ok:
                die(f"A3 floor unexpectedly achievable on {nm}")
            got.append((len(cls), st["states"]))
        rows[nm] = got
        if got != frozen[nm]:
            die(f"A3 frozen row mismatch {nm}: {got} vs {frozen[nm]}")
        bi, bf = brent_fail(X, Y, Z), brent_fail(X, Y, Z, exact=True)
        if bi or bf:
            die(f"A1 brent failed for {nm}")
        print(f"A3 {nm}: d/states {got} floors all impossible, brent 729/729 (int+fmpz)")
    # A4 relabeling symmetry numerically: d and floor verdict invariant under T
    for nm, M in (("U", U), ("V", V), ("W", W)):
        c1, r1 = prep(M)
        c2, r2 = prep(Tblocks(M))
        o1, s1, _ = subset_dfs(c1, r1)
        o2, s2, _ = subset_dfs(c2, r2)
        if (len(c1), o1) != (len(c2), o2):
            die(f"A4 relabeling symmetry violated on {nm}")
    print("A4 relabeling symmetry: d and floor verdicts invariant under T on U,V,W")
    # A5 transposition preconditions per orientation
    pre = {}
    for nm, (X, Y, Z) in orbit.items():
        nz = all(any(M[r]) for M in (X, Y, Z) for r in range(R))
        rk = rank_over_Q([[Z[r][c] for c in range(N)] for r in range(R)])
        pre[nm] = {"all_rows_nonzero": nz, "rank_outfac": rk}
        if not nz or rk != 9:
            die(f"A5 preconditions failed for {nm}: {pre[nm]}")
    print("A5 preconditions ok for all three orientations (rank 9, 23/23 rows)")
    RESULT["phases"]["A-anchors"].update(
        {"sigma0_recount": [nU, nV, nW], "frozen_rows": rows, "preconditions": pre})
    return (U, V, W), (gU, gV), orbit


@phase("B-controls")
def phase_B(blocks, printed):
    from gate_b_floor import prep, subset_dfs
    from gatec_decomps import load_sun56
    U, V, W = blocks
    gU, gV = printed
    # A3-style tri-instrument SAT control on an independent target
    _, V56, _ = load_sun56()
    pcls, _ = prep(V56)
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    a2 = canon((0, 0, 0, 1, 0, -1, 0, 0, 0))
    ext = sorted(set(pcls) | {a1, a2})
    e2c, e2r = prep(ext)
    okd, _, _ = subset_dfs(e2c, e2r)
    sk, _ = run_cnf(ext, 13, tag="ctrl_sun56V_T13")
    sc, _ = run_cnf(ext, 13, tag="ctrl_sun56V_T13_cad", solver=CADICAL)
    print(f"CTRL tri-instrument SAT control: DFS={okd} kissat={sk} cadical={sc}")
    if not (okd and sk and sc):
        die("SAT control failed")
    # R1/R2 plants on the printed left circuit
    okp, _ = circuit_verdict(gU, U)
    if not okp:
        die("R-plants: carrier does not verify")
    pruned = gU[:-1]
    ok1, m1 = circuit_verdict(pruned, U)
    print(f"R1 one-gate-deleted: ACCEPT={ok1} (must be False), missing={len(m1)}")
    if ok1:
        die("R1 not rejected")
    flip = [dict(g) for g in gU]
    flip[3] = dict(flip[3])
    flip[3]["b"] = [flip[3]["b"][0], -flip[3]["b"][1]]
    ok2, _ = circuit_verdict(flip, U)
    print(f"R2 operand-perturbed: ACCEPT={ok2} (must be False)")
    if ok2:
        die("R2 not rejected")
    # R4 frozen Wfac floor@13 with a fresh dual-checker certificate
    wcls, _ = prep(W)
    satw, mw = run_cnf(sorted(wcls), len(wcls), tag="R4_floorW_T13",
                       with_certificates=True, keep_dir=CERT_DIR)
    print(f"R4 floor Wfac at T={len(wcls)}: SAT={satw} (must be False) cert rc "
          f"{mw['cert']['drattrim_rc']}/{mw['cert']['lratcheck_rc']}")
    if satw:
        die("R4 floor unexpectedly SAT")
    # R5 planter
    scratch = os.path.join(HERE, "scratch_R5")
    os.makedirs(scratch, exist_ok=True)
    pl = os.path.join(scratch, "planter.py")
    open(pl, "w").write("import sys\n"
                        f"sys.path.insert(0, {SRC!r})\n"
                        f"sys.path.insert(0, {SCR5!r})\n"
                        f"sys.path.insert(0, {os.path.join(MM3, 'scratch')!r})\n"
                        "from gatec_decomps import load_paper55\n"
                        "from gate_b_floor import prep\n"
                        "U, V, W = load_paper55()\n"
                        "c, _ = prep(W)\n"
                        "assert len(c) == 14, 'planted wrong d'\n")
    rp = subprocess.run([sys.executable, pl], capture_output=True, text=True, env=ENV1)
    fired = rp.returncode != 0 and "planted wrong d" in (rp.stdout + rp.stderr)
    print(f"R5 planter fired: {fired}")
    if not fired:
        die("R5 planter did not fire")
    RESULT["phases"]["B-controls"].update({"SAT_control": "pass (DFS+kissat+CaDiCaL)",
                                           "R1": "REJECTED", "R2": "REJECTED",
                                           "R4_cert": mw["cert"], "R5": "fired"})
    return True


@phase("C-Wfac-synthesis")
def phase_C(blocks):
    from gate_b_floor import prep, subset_dfs
    U, V, W = blocks
    cls, _ = prep(W)
    au1, au2 = aux_universe(cls), aux_universe_alt(cls)
    if au1 != au2:
        die("C universe: the two independent enumeration loops disagree")
    h = csv_hash(au1)
    print(f"C Wfac universe: d={len(cls)} |AU|={len(au1)} sha={h[:16]}… "
          f"(two loops agree)")
    ck = open(os.path.join(HERE, "scan_Wfac_checkpoint.jsonl"), "a")
    T = len(cls) + 1
    hit = None
    t0 = time.time()
    for idx, a in enumerate(au1):
        ext = sorted(set(cls) | {a})
        e2c, e2r = prep(ext)
        ok, st, order = subset_dfs(e2c, e2r)
        sat, _m = run_cnf(ext, T, tag=f"W{idx:03d}")
        if ok != sat:
            die(f"C instrument disagreement at W[{idx}]={a}: DFS={ok} kissat={sat}")
        ck.write(json.dumps({"i": idx, "aux": list(a), "dfs": ok,
                             "states": st["states"], "t": round(time.time() - t0, 1)}) + "\n")
        ck.flush()
        os.fsync(ck.fileno())
        if ok and hit is None:
            hit = (a, order, idx)
        if (idx + 1) % 100 == 0:
            print(f"C progress {idx+1}/{len(au1)} first_hit={hit[2] if hit else None} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    ck.close()
    if hit is None:
        print(f"C COMPLETE: {len(au1)}/{len(au1)} decided, NO 14-gate Wfac circuit")
        RESULT["phases"]["C-Wfac-synthesis"] = {"universe": len(au1), "sha256": h,
                                                "T": T, "admitted": 0}
        return None
    a, order, idx = hit
    gates = replay_schedule(order)
    if gates is None:
        die("C replay failed for the admitted Wfac aux")
    ok, missing = circuit_verdict(gates, W)
    print(f"C admitted aux[{idx}]={list(a)}: {len(gates)}-gate Wfac circuit, "
          f"exact-expansion verified={ok}")
    if not (ok and len(gates) == T):
        die(f"C Wfac witness failed verification ok={ok} gates={len(gates)}")
    json.dump({"aux": list(a), "gates": gates, "T": T,
               "note": "14-gate Wfac circuit synthesized by the aux-1 dual instrument "
                       "and verified by exact expansion; NOT transposition-derived"},
              open(os.path.join(HERE, "witness_Wfac_14gates.json"), "w"), indent=1)
    RESULT["phases"]["C-Wfac-synthesis"] = {"universe": len(au1), "sha256": h, "T": T,
                                            "first_hit_index": idx, "aux": list(a),
                                            "gates": len(gates), "verified": True}
    return gates


def assemble_and_verify(name, left, left_target, right, right_target,
                        ofac, ofac_target, triple):
    """Assemble an orientation's scheme and verify it end to end.
    `left`, `right`, `ofac` are (gate_list, input_perm) pairs; input_perm=None
    means the identity assignment. Input relabeling is free (proved symmetry)."""
    (left_gates, lp), (right_gates, rp), (ofac_gates, op) = left, right, ofac
    okL, misL = circuit_verdict(left_gates, left_target, lp)
    okR, misR = circuit_verdict(right_gates, right_target, rp)
    if not (okL and okR):
        die(f"{name}: stage verification failed L={okL} ({misL[:3]}) "
            f"R={okR} ({misR[:3]})")
    tr = transpose_circuit(ofac_gates, ofac_target, op)
    if tr is None:
        die(f"{name}: output-stage construction failed (target row not readable)")
    oadds, outs, meta = tr
    X, Y, Z = triple
    ok_out = all(outs[k][r] == Z[r][k] for k in range(N) for r in range(R))
    nL, nR = len(left_gates), len(right_gates)
    total = nL + nR + oadds
    print(f"{name}: left={nL} right={nR} output={oadds} total={total} "
          f"output-map exact={ok_out}")
    if not ok_out:
        return {"status": "OUTPUT-MAP-MISMATCH", "left": nL, "right": nR,
                "output": oadds, "total": total}
    # END-TO-END: expand the assembled scheme over the 81 monomials A_i*B_j
    lvals = eval_circuit(left_gates, lp)
    rvals = eval_circuit(right_gates, rp)
    if lvals is None or rvals is None:
        die(f"{name}: stage re-evaluation inconsistent")

    def pick(vals, tv):
        for wid, v in vals.items():
            if list(v) == list(tv):
                return v
            if all(v[k] == -tv[k] for k in range(N)):
                return tuple(-x for x in v)
        return None
    prods = []
    for r in range(R):
        lv, rv = pick(lvals, X[r]), pick(rvals, Y[r])
        if lv is None or rv is None:
            die(f"{name}: product {r} operand not realisable from the stage circuits")
        M = [[0] * N for _ in range(N)]
        for i in range(N):
            if lv[i]:
                for j in range(N):
                    if rv[j]:
                        M[i][j] += lv[i] * rv[j]
        prods.append(M)
    # outputs from the constructed output stage's coefficient vectors
    couts = []
    for k in range(N):
        acc = [[0] * N for _ in range(N)]
        for r in range(R):
            c = outs[k][r]
            if c:
                for i in range(N):
                    for j in range(N):
                        if prods[r][i][j]:
                            acc[i][j] += c * prods[r][i][j]
        couts.append(acc)
    # target bilinear map defined by the orientation's triple
    tgt = []
    for k in range(N):
        acc = [[0] * N for _ in range(N)]
        for r in range(R):
            if Z[r][k]:
                for i in range(N):
                    if X[r][i]:
                        for j in range(N):
                            if Y[r][j]:
                                acc[i][j] += Z[r][k] * X[r][i] * Y[r][j]
        tgt.append(acc)
    bad = [k for k in range(N) if couts[k] != tgt[k]]
    print(f"{name}: END-TO-END outputs matching the orientation's bilinear map: "
          f"{N - len(bad)}/9 (mismatches {bad})")
    if bad:
        return {"status": "END-TO-END-MISMATCH", "left": nL, "right": nR,
                "output": oadds, "total": total, "bad": bad}
    return {"status": "OK", "left": nL, "right": nR, "output": oadds, "total": total,
            "outputs_exact": 9}


@phase("D-sigma1")
def phase_D(blocks, printed, wfac_gates, orbit):
    U, V, W = blocks
    gU, gV = printed
    X, Y, Z = orbit["sigma^1"]          # (V, T(W), T(U))
    # left map = V  -> printed RIGHT_SLP (14 gates), inputs renamed B->A (free)
    # right map = T(W) -> synthesized 14-gate Wfac circuit with T-permuted inputs
    # ofac map  = T(U) -> printed LEFT_SLP (13 gates) with T-permuted inputs
    res = assemble_and_verify("sigma^1",
                              (gV, None), X,               # left  = V map
                              (wfac_gates, T_INV), Y,      # right = T(W) map
                              (gU, T_INV), Z,              # ofac  = T(U) map
                              (X, Y, Z))
    RESULT["phases"]["D-sigma1"] = res
    if res["status"] == "OK":
        json.dump({"orientation": "sigma^1 = (V, T(W), T(U))",
                   "left_additions": res["left"], "right_additions": res["right"],
                   "output_additions": res["output"], "total": res["total"],
                   "note": "left = printed RIGHT_SLP under input renaming; right = "
                           "in-run synthesized Wfac circuit with T-permuted inputs; "
                           "output = reverse-mode transposition of the printed "
                           "LEFT_SLP with T-permuted inputs, verified by exact "
                           "expansion and by end-to-end 81-monomial expansion"},
                  open(os.path.join(HERE, "assembled_sigma1.json"), "w"), indent=1)
    return res


@phase("E-sigma2")
def phase_E(blocks, printed, wfac_gates, orbit):
    U, V, W = blocks
    gU, gV = printed
    X, Y, Z = orbit["sigma^2"]          # (T(W), U, T(V))
    res = assemble_and_verify("sigma^2",
                              (wfac_gates, T_INV), X,      # left  = T(W) map
                              (gU, None), Y,               # right = U map
                              (gV, T_INV), Z,              # ofac  = T(V) map
                              (X, Y, Z))
    RESULT["phases"]["E-sigma2"] = res
    if res["status"] == "OK":
        json.dump({"orientation": "sigma^2 = (T(W), U, T(V))",
                   "left_additions": res["left"], "right_additions": res["right"],
                   "output_additions": res["output"], "total": res["total"],
                   "note": "left = in-run synthesized Wfac circuit with T-permuted "
                           "inputs; right = printed LEFT_SLP under input renaming; "
                           "output = reverse-mode transposition of the printed "
                           "RIGHT_SLP with T-permuted inputs, verified by exact "
                           "expansion and by end-to-end 81-monomial expansion"},
                  open(os.path.join(HERE, "assembled_sigma2.json"), "w"), indent=1)
    return res



@phase("G-verdict")
def phase_G(cwfac, d1, d2):
    blob = {"gate": "paper55-sigma12-total55", "run_id": os.path.basename(HERE),
            "prereg_sha256": "10b3f32f1f4c3c61fa06a94307242d33c03951e55b2520513e0fd581c42f4b64",
            "decomposition": "paper55 (arXiv:2607.28676) sigma^1 and sigma^2 orientations",
            "question": "decide each orientation's minimum inside the frozen [55, 57]",
            "reduction": ("all three side costs exact per orientation by the proved "
                          "fixed-coordinate-relabeling symmetry applied to the frozen "
                          "gate-B values C(U)=13, C(V)=14, C(Wfac)=14; the open "
                          "quantity was whether each output stage attains its bound"),
            "wfac_synthesis": RESULT["phases"].get("C-Wfac-synthesis", {}),
            "sigma1": d1, "sigma2": d2,
            "controls": RESULT["phases"].get("B-controls", {}),
            "cpu_seconds_total": round(sum(x.get("seconds", 0)
                                           for x in RESULT["phases"].values()), 1)}
    ok1 = d1 and d1.get("status") == "OK"
    ok2 = d2 and d2.get("status") == "OK"
    if ok1 and ok2 and d1["total"] == 55 and d2["total"] == 55:
        blob["verdict"] = "CERTIFIED-55-BOTH"
        blob["claim"] = (
            "Both sigma^1 and sigma^2 attain 55: with every side cost exact "
            "(14/14/13 and 14/13/14 respectively, inherited from the frozen gate-B "
            "exact values through the proved free relabeling C(T(X)) = C(X)) and an "
            "output stage constructed and VERIFIED at 27 and 28 additions, the "
            "assembled schemes total 14+14+27 = 55 and 14+13+28 = 55, each verified "
            "end to end by exact symbolic expansion over all 81 monomials against "
            "that orientation's bilinear map. Both orientations therefore have "
            "minimum EXACTLY 55 (LB 55 frozen; UB 55 attained), superseding their "
            "'best known 57' rows. This is the record VALUE at two further "
            "orientations of the same tensor, NOT a new record: paper55 sigma^0's 55 "
            "remains the record.")
        blob["LB_UB_after"] = {"sigma^1": [55, 55], "sigma^2": [55, 55]}
    elif ok1 or ok2:
        blob["verdict"] = "CERTIFIED-PARTIAL"
        blob["claim"] = (f"sigma^1 status {d1 and d1.get('status')} total "
                         f"{d1 and d1.get('total')}; sigma^2 status "
                         f"{d2 and d2.get('status')} total {d2 and d2.get('total')} — "
                         "see report for the exact per-orientation reading.")
    else:
        blob["verdict"] = "INCONCLUSIVE"
        blob["claim"] = ("neither orientation was assembled and verified; per the "
                         "pre-registered rule no total claim is made and both "
                         "intervals stay [55, 57].")
    json.dump(blob, open(os.path.join(HERE, "verdict.json"), "w"), indent=1)
    print("G verdict:", blob["verdict"])
    print(blob["claim"][:700])
    return blob


if __name__ == "__main__":
    for nm, (pth, got, ok) in CHECKERS.items():
        print(f"checker {nm}: sha256={got[:16]}… pinned={ok}")
        if not ok:
            die(f"checker {nm} hash mismatch")
    blocks, printed, orbit = phase_A()
    phase_B(blocks, printed)
    wfac = phase_C(blocks)
    if wfac is None:
        phase_G(None, None, None)
    else:
        d1 = phase_D(blocks, printed, wfac, orbit)
        d2 = phase_E(blocks, printed, wfac, orbit)
        phase_G(wfac, d1, d2)
    json.dump(RESULT, open(os.path.join(HERE, "runner_output.json"), "w"),
              indent=1, default=str)
    print("\nRUNNER COMPLETE")
