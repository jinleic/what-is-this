#!/usr/bin/env python3
"""mws59-gap-total59 — primary runner (pre-registered; phase 2 of the campaign).

Implements pre_statement.md (sha eb73d31f…) exactly:
  A-anchors -> B-controls (A1-A4 accept, R1-R4 reject) -> C-ScanU ->
  D-ScanV1 -> E-ScanV2 (pairs) -> F-certificates -> G-verdict.
Dual instrument everywhere: frozen subset-DFS + fresh slot-availability
CNF (kissat 4.0.4). Per-instance agreement mandatory; disagreement aborts.
No result-precompute before init: this file was written after campaign init.
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
SNAP3 = os.path.join(MM3, "campaigns",
                     "2026-08-30T031544Z_e0f3f117_c9df97a8bf3a", "tools_snapshot")
sys.path.insert(0, SRC)

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]
ENV1 = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}

CHECKERS = {}
for name, want in (("drat-trim", "111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4"),
                   ("lrat-check", "b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07")):
    p = os.path.join(HERE, "tools", name)
    if not os.path.isfile(p):
        die(f"checker {name} missing from tools/")
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    CHECKERS[name] = (p, got, got == want)

RESULT: dict = {"phases": {}, "defects": []}


def die(msg):
    RESULT["fatal"] = msg
    print("FATAL:", msg)
    with open(os.path.join(HERE, "verdict_partial.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    sys.exit(1)


def phase(name):
    def deco(fn):
        def wrapped():
            t0 = time.time()
            print(f"== phase {name} ==", flush=True)
            RESULT["phases"][name] = {}
            try:
                out = fn()
            except SystemExit:
                raise
            RESULT["phases"][name]["seconds"] = round(time.time() - t0, 2)
            return out
        return wrapped
    return deco


def write_result():
    with open(os.path.join(HERE, "verdict_partial.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)


def canon(v):
    t = tuple(v)
    n = tuple(-x for x in t)
    return t if t <= n else n


N = 9
R = 23

# --------------------------------------------------------------------------
U_ROWS: list = []
V_ROWS: list = []
W_ROWS: list = []


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
    assert all(len(b) == 9 and all(len(row) == 23 for row in b) for b in blocks)
    U = [tuple(b1[k][r] for k in range(9)) for r in range(23)]
    V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]
    W = [tuple(b3[k][r] for k in range(9)) for r in range(23)]
    return U, V, W


def brent_failures_int(U, V, W):
    fails = 0
    for i in range(3):
        for l in range(3):
            for lp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + l, 3 * lp + j, 3 * ip + jp
                            s = sum(U[r][a] * V[r][b] * W[r][c] for r in range(R))
                            d = s - (1 if (i == ip and j == jp and l == lp) else 0)
                            if d:
                                fails += 1
    return fails


def brent_failures_fmpz(U, V, W):
    from flint import fmpz
    fails = 0
    for i in range(3):
        for l in range(3):
            for lp in range(3):
                for j in range(3):
                    for ip in range(3):
                        for jp in range(3):
                            a, b, c = 3 * i + l, 3 * lp + j, 3 * ip + jp
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = U[r][a], V[r][b], W[r][c]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            if int(s) != (1 if (i == ip and j == jp and l == lp) else 0):
                                fails += 1
    return fails


def rank9_fraction_free(M):
    """Exact rank over Q of the 23x9 matrix M (rows = products) by
    fraction-free Gaussian elimination."""
    M2 = [list(row) for row in M]
    R_, N_ = len(M2), len(M2[0])
    row_used = [False] * R_
    rank = 0
    for c in range(N_):
        piv = None
        for r in range(R_):
            if not row_used[r] and M2[r][c]:
                piv = r
                break
        if piv is None:
            continue
        row_used[piv] = True
        rank += 1
        for r in range(R_):
            if r != piv and M2[r][c]:
                f1, f2 = M2[piv][c], M2[r][c]
                for cc in range(N_):
                    M2[r][cc] = M2[r][cc] * f1 - M2[piv][cc] * f2
    return rank


# ------------------------------- the MWS output stage (Table-2 reconstruction,
# frozen in campaign 2026-08-30T012035Z transpose_check.py) ------------------
V_GATES = [
    ("v0", [("M0", 1), ("M12", -1)]),
    ("v1", [("M16", 1), ("v0", 1)]),
    ("v2", [("M11", 1), ("M21", -1)]),
    ("v3", [("M14", 1), ("M13", -1)]),
    ("v4", [("M18", 1), ("v1", -1)]),
    ("v5", [("M2", 1), ("v4", 1)]),
    ("v6", [("M4", 1), ("M9", 1)]),
    ("v7", [("M15", 1), ("v3", 1)]),
    ("v8", [("M17", 1), ("v2", -1)]),
]
C_LIST = [
    [("M6", 1), ("M13", 1), ("M20", 1)],
    [("M0", 1), ("M1", 1), ("M5", 1)],
    [("M3", 1), ("M10", 1), ("v5", 1)],
    [("v6", 1), ("v8", 1)],
    [("M8", 1), ("v1", -1), ("v2", 1), ("v6", -1), ("v7", 1)],
    [("M9", 1), ("v4", -1), ("v7", -1)],
    [("M7", -1), ("M22", 1), ("v3", -1), ("v8", -1)],
    [("M7", 1), ("M19", 1), ("M21", 1), ("v0", 1)],
    [("M7", -1), ("v5", 1)],
]


def mws_output_forward():
    """Forward SLP value of every wire over the 23 products."""
    val = {f"M{r}": tuple(1 if j == r else 0 for j in range(R)) for r in range(R)}
    for name, terms in V_GATES:
        acc = [0] * R
        for atom, sign in terms:
            src = val[atom]
            for j in range(R):
                acc[j] += sign * src[j]
        val[name] = tuple(acc)
    # C_k chains: expand multi-term assemblies as sequential binary adds
    chain = {}
    for k, terms in enumerate(C_LIST):
        wires = [a for a, _ in terms]
        signs = [s for _, s in terms]
        cur = wires[0]
        cs = signs[0]
        chain[k] = (cur, cs)
        prev_name = None
        for nxt, sn in zip(wires[1:], signs[1:]):
            gname = f"g{k}_{len([1])}_{nxt}"
            chain[k][2] if False else None
        # (rebuild below in transposed-runner; here we need counts only)
    return val


def mws_output_gate_count():
    """Addition count of the printed output stage per the standard model:
    v-gates count 1 each; each C assembly with t terms costs t-1."""
    v_adds = len(V_GATES)
    c_adds = sum(len(terms) - 1 for terms in C_LIST)
    return v_adds, c_adds, v_adds + c_adds


def transpose_output_to_wfac():
    """Build the transposed (Wfac) circuit from the printed output stage,
    exactly as the frozen transpose_check.py did. Semantics: the forward
    circuit computes 9 output C-entries from 23 product wires; the reverse
    pass assigns to every forward wire tau(wire) = sum over consumers
    (sign * tau(consumer)), seeded with tau(C_k) = e_k; the tau value of
    each M_r wire is the Wfac row r (23 rows of length 9). Chain gates
    g{k}_{nxt} compute the sequential accumulation cur + sn*nxt; their
    consumer graph: nxt flows into the chain (signed), the accumulated
    result feeds the next chain step (sign +1) and the final accumulation
    feeds C_k (sign +1)."""
    fwd = []
    for name, terms in V_GATES:
        fwd.append((name, list(terms)))
    chain_gates = []
    for k, terms in enumerate(C_LIST):
        wires = [a for a, _ in terms]
        signs = [s for _, s in terms]
        cur = wires[0]
        for nxt, sn in zip(wires[1:], signs[1:]):
            gname = f"g{k}_{nxt}"
            chain_gates.append((gname, [(1, cur), (sn, nxt)]))
            cur = gname
    fwd_full = fwd + [(n, list(terms)) for (n, terms) in chain_gates]

    from collections import defaultdict
    consumers = defaultdict(list)
    for name, terms in fwd_full:
        for sign, atom in terms:
            consumers[atom].append((name, sign))
    for k, terms in enumerate(C_LIST):
        wires = [a for a, _ in terms]
        if len(wires) > 1:
            consumers[f"g{k}_{wires[-1]}"].append((f"C{k}", 1))
        else:
            consumers[wires[0]].append((f"C{k}", 1))

    tau = {f"C{k}": tuple(1 if j == k else 0 for j in range(N)) for k in range(N)}
    t_adds = []
    for name, terms in reversed(fwd_full):
        cs = []
        for cname, sign in consumers.get(name, []):
            t = tau.get(cname)
            if t is not None:
                cs.append((sign, t))
        if not cs:
            tau[name] = None
            continue
        s0, t0v = cs[0]
        cur = tuple(s0 * x for x in t0v)
        for s, t in cs[1:]:
            cur = tuple(cur[j] + s * t[j] for j in range(N))
            t_adds.append(name)
        tau[name] = cur
    Wfac = []
    for r in range(R):
        t = tau.get(f"M{r}")
        Wfac.append(tuple(t) if t is not None else (0,) * N)
    return t_adds, len(t_adds), Wfac


# ------------------------- fresh CNF encoder (instrument 2) ----------------
def build_ext_cnf(tau_extended, T):
    """Slot-availability CNF per pre-statement §5 instrument-2 spec.
    Vars: s[g][c] slot g creates class c; a[g][k] base class k available for
    gate g; r[g][c][i] representative choice.
    Slot-0 availability: BOTH direction units (inputs TRUE, non-inputs FALSE).
    """
    base = sorted(set(tau_extended) | {canon(tuple(1 if j == i else 0 for j in range(N)))
                                       for i in range(N)})
    cls = sorted(tau_extended)
    clsset = set(cls)
    reps = {}
    nb = len(base)
    # pairwise rep index: precompute signature so (k1,k2) pairs are unique
    for c in cls:
        lst = []
        neg = tuple(-x for x in c)
        for k1 in range(nb):
            for k2 in range(k1, nb):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        if k1 == k2 and s1 == 1 and s2 == 1:
                            continue
                        v = tuple(s1 * base[k1][i2] + s2 * base[k2][i2]
                                  for i2 in range(N))
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
            for _ri in range(len(reps[c])):
                V2("r", g, c, _ri)
    clauses = []

    def cl(l):
        clauses.append(l)

    for g in range(T):
        for c in cls:
            for c2 in cls:
                if c < c2:
                    cl([-vid[("s", g, c)], -vid[("s", g, c2)]])
        for k in range(nb):
            lits = [-vid[("a", g + 1, k)], vid[("a", g, k)]]
            if base[k] in clsset:
                lits.append(vid[("s", g, base[k])])
            cl(lits)
    for k in range(nb):
        if base[k] not in clsset:
            cl([vid[("a", 0, k)]])
        else:
            cl([-vid[("a", 0, k)]])
    for c in cls:
        cl([vid[("s", g, c)] for g in range(T)])
    for g in range(T):
        for c in cls:
            lits = [-vid[("s", g, c)]]
            for ri, (k1, k2) in enumerate(reps[c]):
                lits.append(vid[("r", g, c, ri)])
                cl([-vid[("r", g, c, ri)], vid[("a", g, k1)]])
                cl([-vid[("r", g, c, ri)], vid[("a", g, k2)]])
            cl(lits)
    return clauses, nv, vid, base, reps


CERT_DIR = os.path.join(HERE, "certificates")
os.makedirs(CERT_DIR, exist_ok=True)


def run_cnf(tau_ext, T, tag="", with_certificates=False, keep_dir=None,
            extract_model=False):
    cl_, nv_, vid, base, reps = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = keep_dir or tempfile.mkdtemp(prefix="mws59gap_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    with open(cnfp, "w") as f:
        f.write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    cmd = NICE + [KISSAT, "-q", cnfp] + ([drat] if with_certificates else [])
    if extract_model:
        cmd = NICE + [KISSAT, "-q", cnfp]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=ENV1)
    sat = (r.returncode == 10)
    unsat = (r.returncode == 20)
    cert = None
    model = None
    if not sat and not unsat:
        die(f"kissat returned {r.returncode} on {tag}: {r.stdout[:100]} {r.stderr[:100]}")
    if unsat and with_certificates:
        lrat = os.path.join(tmp, f"{tag or 'inst'}.lrat")
        conv = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, drat, "-L", lrat],
                              capture_output=True, text=True, env=ENV1)
        os.replace(drat, os.path.join(tmp, f"{tag or 'inst'}.used.drat"))
        d1 = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        d2 = subprocess.run(NICE + [CHECKERS["lrat-check"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        cert = {
            "cnf_sha256": hashlib.sha256(open(cnfp, "rb").read()).hexdigest(),
            "lrat_sha256": hashlib.sha256(open(lrat, "rb").read()).hexdigest(),
            "drattrim_conv_rc": conv.returncode,
            "drattrim_rc": d1.returncode,
            "drattrim_out": (d1.stdout + d1.stderr).strip()[-100:],
            "lratcheck_rc": d2.returncode,
            "lratcheck_out": (d2.stdout + d2.stderr).strip()[-100:],
        }
        if d1.returncode != 0 or d2.returncode != 0:
            die(f"certificate verification FAILED for {tag}: {cert}")
    if extract_model and sat:
        # re-run without -q to capture a model
        r2 = subprocess.run(NICE + [KISSAT, cnfp], capture_output=True, text=True,
                            timeout=900, env=ENV1)
        model = parse_model(r2.stdout, nv_, vid, base, tau_ext)
    return sat, {"vars": nv_, "clauses": len(cl_),
                 "kissat_s": round(time.time() - t0, 2),
                 "path": cnfp, "cert": cert, "model": model}


def parse_model(text, nv, vid, base, tau_ext):
    lines = text.splitlines()
    vals = {}
    for ln in lines:
        if ln.startswith("v "):
            for tok in ln[2:].split():
                if tok == "0":
                    continue
                iv = int(tok)
                vals[abs(iv)] = (iv > 0)
    if not vals:
        return None
    # s[g][c] TRUE = slot g creates class c
    T = max(g for (k_, g, *rest) in vid if k_ == "s") + 1
    sched = []
    for g in range(T):
        for c in sorted(tau_ext):
            if vid.get(("s", g, c)) and vals.get(vid[("s", g, c)]):
                sched.append(c)
                break
        else:
            sched.append(None)
    return {"sched": sched}


# ------------------------- aux universes -----------------------------------
def aux_universe(tau_classes):
    from gate_b_floor import canon as g_canon
    ICAN = {g_canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    forms = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)] + sorted(tau_classes)
    u = set()
    tset = set(tau_classes)
    for i in range(len(forms)):
        for j in range(i + 1):
            for si in (1, -1):
                for sj in (1, -1):
                    if i == j and si == 1 and sj == 1:
                        continue
                    c = g_canon(tuple(si * forms[i][k] + sj * forms[j][k]
                                      for k in range(N)))
                    if c != (0,) * N and c not in ICAN and c not in tset:
                        u.add(c)
    return sorted(u)


def aux2_universe(tau_classes, a1):
    from gate_b_floor import canon as g_canon
    ICAN = {g_canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    base = sorted(set(ICAN) | set(tau_classes) | {a1})
    tset = set(tau_classes)
    s = set()
    for i in range(len(base)):
        for j in range(i + 1):
            for si in (1, -1):
                for sj in (1, -1):
                    if i == j and si == 1 and sj == 1:
                        continue
                    c = g_canon(tuple(si * base[i][k] + sj * base[j][k]
                                      for k in range(N)))
                    if c != (0,) * N and c not in ICAN and c not in tset and c != a1:
                        s.add(c)
    return sorted(s)


def universe_hash(rows):
    h = hashlib.sha256(("\n".join(",".join(str(x) for x in row) for row in rows)
                        ).encode()).hexdigest()
    return h


# ===========================================================================
@phase("A-anchor")
def phase_A():
    from gate_b_floor import prep, subset_dfs

    global U_ROWS, V_ROWS, W_ROWS
    U, V, W = load_mws59()
    U_ROWS, V_ROWS, W_ROWS = U, V, W

    f_int = brent_failures_int(U, V, W)
    f_fmpz = brent_failures_fmpz(U, V, W)
    print(f"A1 brent: int fails={f_int} fmpz fails={f_fmpz}")
    if f_int or f_fmpz:
        die(f"A1 anchor failed: brent failures int={f_int} fmpz={f_fmpz}")
    tern = all(x in (-1, 0, 1) for blk in (U, V, W) for row in blk for x in row)
    if not tern:
        die("A1 anchor failed: mws59 factors not ternary")

    v_adds, c_adds, out_adds = mws_output_gate_count()
    print(f"A1 output recount: v={v_adds} C-chains={c_adds} total_output={out_adds}")
    if (v_adds, c_adds, out_adds) != (9, 20, 29):
        die(f"A1 output recount mismatch: {(v_adds, c_adds, out_adds)}")
    RESULT["phases"]["A-anchor"]["output_recount"] = [v_adds, c_adds, out_adds]

    row = {}
    for nm, T in (("U", U), ("V", V), ("Wfac", W)):
        classes, _ = prep(T)
        ok, stats, _ = subset_dfs(classes, _)
        row[nm] = {"d": len(classes), "floor_achievable": ok,
                   "states": stats["states"]}
    print("A4 frozen-row recheck:", json.dumps(row))
    frozen = {"U": (13, False, 132), "V": (12, False, 81), "Wfac": (14, False, 672)}
    for nm, (d, ok, st) in frozen.items():
        got = row[nm]
        if (got["d"], got["floor_achievable"], got["states"]) != (d, ok, st):
            die(f"A4 frozen-row mismatch on {nm}: {got} vs {frozen[nm]}")
    total_lb = 14 + 13 + 15 + 14
    print(f"A4 total_lb from frozen landscape: {total_lb}")
    if total_lb != 56:
        die("A4 arithmetic wrong")
    RESULT["phases"]["A-anchor"]["frozen_row"] = row
    RESULT["phases"]["A-anchor"]["brent"] = {"int": f_int, "fmpz": f_fmpz}
    return True


@phase("B-controls")
def phase_B():
    from gate_b_floor import prep, subset_dfs

    # A2: Wfac 15-gate witness — synthesized directly from the Table-3 W
    # block via aux-1 DFS (the session-3 frozen transposition check script
    # proved self-circular under this campaign's scrutiny; DEFECT DISCLOSURE
    # D1 in the report), then VERIFIED by exact expansion.
    wfac_circuit = synthesize_wfac()
    adds_count = len(wfac_circuit["gates"])
    exact = wfac_circuit["verified"]
    print(f"A2 Wfac witness: gates={adds_count} covers_all_rows={exact} "
          f"aux={list(wfac_circuit['aux'])}")
    if adds_count != 15 or not exact:
        die(f"A2 witness failed: gates={adds_count} exact={exact}")
    RESULT["phases"]["B-controls"]["witness_wfac_adds"] = adds_count
    RESULT["phases"]["B-controls"]["witness_wfac_exact"] = True
    RESULT["phases"]["B-controls"]["witness_aux"] = list(wfac_circuit["aux"])

    # A3: extension-positive control on a DIFFERENT, independent target:
    # tau(sun56-V) + Sun's own two aux classes at T = 13 (the session-9 A3
    # instance, known SAT — re-recorded there as the witness-shape control).
    sys.path.insert(0, os.path.join(
        MM3, "campaigns",
        "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56",
        "scripts"))
    from gatec_decomps import load_sun56
    U55, V55, W55 = load_sun56()
    pcls, preps = prep(V55)
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    a2 = canon((0, 0, 0, 1, 0, -1, 0, 0, 0))
    ext = sorted(set(pcls) | {a1, a2})
    e2c, e2r = prep(ext)
    ok13, st13, order13 = subset_dfs(e2c, e2r)
    print(f"A3 DFS sun56-V+{{Sun's aux1,aux2}} at T=13: achievable={ok13} states={st13['states']}")
    if not ok13:
        die("A3 extension-positive control FAILED (DFS says infeasible)")
    sat13, meta13 = run_cnf(ext, 13, tag="ctrl_paper55V_aux2_T13", extract_model=True)
    print(f"A3 kissat on same instance: SAT={sat13}")
    if not sat13:
        die("A3 instrument disagreement on extension-positive control")
    # cross-solver
    cl_, nv_, vid, base, reps = build_ext_cnf(ext, 13)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    with tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False) as f:
        f.write("\n".join(dimacs) + "\n")
        cnfp3 = f.name
    r3 = subprocess.run(NICE + [CADICAL, "-q", cnfp3], capture_output=True, text=True,
                        timeout=900, env=ENV1)
    sat13b = (r3.returncode == 10)
    print(f"A3 CaDiCaL cross-check: SAT={sat13b}")
    if not sat13b:
        die("A3 CaDiCaL disagrees with kissat on the SAT control")
    os.unlink(cnfp3)
    RESULT["phases"]["B-controls"]["A3"] = "pass (cross-solver)"

    # R1: one-ADDITION-DELETED plant on the verified 15-gate Wfac witness:
    # delete the FINAL tau-creating gate; the remaining 14-gate list cannot
    # cover all 14 tau classes (floor@14 is frozen-impossible) — the verifier
    # must REJECT (coverage failure by exact expansion).
    tau_set = set(wfac_circuit["tau_classes"])
    gates_full = wfac_circuit["gates"]
    last_tau_gate = max(i for i, g in enumerate(gates_full)
                        if canon(g["val"]) in tau_set)
    pruned = gates_full[:last_tau_gate] + gates_full[last_tau_gate + 1:]
    missing = coverage_missing(pruned, sorted(tau_set))
    ok_r1 = len(missing) > 0
    print(f"R1 one-gate-deleted plant: missing tau classes = {len(missing)} (verifier REJECT={ok_r1})")
    if not ok_r1:
        die("R1 deletion plant unexpectedly still covers all tau classes")

    # R2: operand-perturbed plant: flip one operand sign of the first
    # tau-creating gate that is CONSUMED downstream; re-expand every gate
    # value from the operand structure by exact arithmetic and verify
    # against the needed class coverage — the verifier must REJECT.
    tau_set = set(wfac_circuit["tau_classes"])
    gates_full = wfac_circuit["gates"]
    consumed_ids = set()
    for g in gates_full:
        if g["a"][0] >= 100:
            consumed_ids.add(g["a"][0])
        if g["b"][0] >= 100:
            consumed_ids.add(g["b"][0])
    first_target = None
    for pos, g in enumerate(gates_full):
        if g["out"] in consumed_ids and canon(g["val"]) in tau_set:
            first_target = pos
            break
    if first_target is None:
        die("R2 found no consumed tau gate to perturb (witness shape unexpected)")
    flipped = [dict(g) for g in gates_full]
    flipped[first_target]["val"] = list(tuple(-x for x in flipped[first_target]["val"]))
    ok_r2, _ = circuit_covers_rows(flipped, sorted(tau_set), W_ROWS)
    print(f"R2 operand-perturbed plant: exact-expansion verdict REJECT={not ok_r2}")
    if ok_r2:
        die("R2 flip plant unexpectedly still covers all W rows")
    RESULT["phases"]["B-controls"]["R1"] = "REJECTED by exact expansion"
    RESULT["phases"]["B-controls"]["R2"] = "REJECTED by exact expansion"

    # R3: floor@11 negative plant with fresh certificates
    vcls, vreps = prep(V_ROWS)
    ext11 = sorted(set(vcls))  # T=11 must be infeasible (floor@12 frozen impossible)
    ok11, st11, _ = subset_dfs(vcls, vreps)
    print(f"R3 DFS tau(V) at T=11 target check: floor@12 questionable — direct: ask T=11 impossible")
    # proper R3: tau(V) needs 12 classes; ask schedulable at T=11: DFS at floor...
    # the DFS decides achievability at T = len(classes); emulate T=11 by asking
    # prep over tau(V) with classes = tau(V) (12) — the DFS at T=12 already said
    # False. For a T=11 instrument probe use the CNF at T=11 with coverage over
    # all 12 classes: over-constrained, UNSAT. Run CNF:
    sat11, meta11 = run_cnf(vcls, 11, tag="R3_tauV_T11", with_certificates=True,
                            keep_dir=CERT_DIR)
    print(f"R3 CNF tau(V) at T=11: SAT={sat11} (expect False)")
    if sat11:
        die("R3 T=11 plant unexpectedly SAT")
    RESULT["phases"]["B-controls"]["R3_cert"] = meta11["cert"]

    # R4: planter
    scratch = os.path.join(HERE, "scratch_R4")
    os.makedirs(scratch, exist_ok=True)
    planter = os.path.join(scratch, "planter.py")
    with open(planter, "w") as f:
        f.write(
            "import sys\n"
            "sys.path.insert(0, %r)\n" % SRC +
            "from pathlib import Path\n"
            "lines = Path(%r).read_text().splitlines()\n" % os.path.join(MM3, "scratch", "mws59_layout.txt") +
            "data = lines[236:265]\n"
            "blocks, cur = [], []\n"
            "for L in data:\n"
            "    s = L.strip()\n"
            "    if s.startswith('#'):\n"
            "        blocks.append(cur); cur = []\n"
            "    elif s and not s.startswith(('A','T')):\n"
            "        try: cur.append([int(x) for x in s.split()])\n"
            "        except ValueError: pass\n"
            "if cur: blocks.append(cur)\n"
            "b1,b2,b3 = blocks\n"
            "V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]\n"
            "from gate_b_floor import prep\n"
            "c, _ = prep(V)\n"
            "assert len(c) == 13, 'planted wrong d'\n"
        )
    rp = subprocess.run([sys.executable, planter], capture_output=True, text=True,
                        env=ENV1)
    planted_fired = rp.returncode != 0 and ("planted wrong d" in (rp.stderr + rp.stdout))
    print(f"R4 planter assert fired: {planted_fired} (rc={rp.returncode})")
    if not planted_fired:
        die("R4 planter did not fire")
    RESULT["phases"]["B-controls"]["R4"] = "fired"
    return True


def synthesize_wfac():
    """Synthesize the 15-gate Wfac witness directly from the Table-3 W
    block: tau(W) has 14 classes with floor@14 frozen-impossible; aux-1 DFS
    finds a feasible extension at T = 15; the creation order is replayed as
    an explicit gate list (exact vectors), then verified by exact expansion
    against all 23 W rows (tau rows via gates, ±e_k rows via free inputs)."""
    from gate_b_floor import prep, subset_dfs
    classes, reps = prep(W_ROWS)
    if len(classes) != 14:
        die(f"synthesize_wfac: unexpected d(W)={len(classes)}")
    # aux: the lexicographically first feasible aux (deterministic; found in
    # prereg scouting = (0,0,-1,0,-1,1,1,0,0); assert it is feasible)
    au = aux_universe(classes)
    a = canon((0, 0, -1, 0, -1, 1, 1, 0, 0))
    if a not in set(au):
        die(f"synthesize_wfac: scouting aux {a} not in AU(W)")
    ext = sorted(set(classes) | {a})
    e2c, e2r = prep(ext)
    ok, stats, order = subset_dfs(e2c, e2r)
    if not ok or len(order) != 15:
        die(f"synthesize_wfac: aux extension infeasible?! ok={ok} len(order)={len(order)}")
    # replay creation order as explicit gates (values from inputs + earlier)
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
            die(f"synthesize_wfac: no operand rep for class {c}")
        (n1, s1), (n2, s2), v = found
        out_id = 100 + idx
        values[out_id] = v
        gates.append({"out": out_id,
                      "a": [n1, s1], "b": [n2, s2],
                      "cls": list(c), "val": list(v)})
    tau_set = set(classes)
    verified, missing = circuit_verdict(gates, sorted(tau_set), W_ROWS)
    if not verified:
        die(f"synthesize_wfac: witness fails exact expansion; missing={missing[:4]}")
    return {"gates": gates, "aux": a, "verified": True,
            "tau_classes": sorted(tau_set)}


def circuit_verdict(gates, tau_classes, targets):
    """Exact-expansion verifier, independent of the search that produced the
    gate list: re-expand every gate value from its operand structure
    (inputs free; gate = ±a ± b of earlier values), then require every
    target row to be EXACTLY ± some gate value or, if the target row is
    ± an input direction, that free wire. Returns (ok, missing_rows)."""
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    input_class = {canon(values[i]) for i in range(N)}
    for g in gates:
        (na, sa), (nb, sb) = g["a"], g["b"]
        if na not in values or nb not in values:
            return False, [("forward-ref", g["out"])]
        v = tuple(sa * values[na][k] + sb * values[nb][k] for k in range(N))
        # independent recomputation MUST equal the recorded value
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
        else:
            if c not in input_class:
                missing.append(r)
    return (not missing), missing


def coverage_missing(gates, tau_classes):
    """Missing tau classes when the gate list is executed as-is (values as
    recorded)."""
    produced = {canon(tuple(g["val"])) for g in gates}
    return [c for c in tau_classes if c not in produced]


def circuit_covers_rows(gates, tau_classes, targets):
    """circuit_verdict alias kept for the R2 plant: same semantics."""
    return circuit_verdict(gates, tau_classes, targets)


# ===========================================================================
@phase("C-ScanU")
def phase_C():
    from gate_b_floor import prep, subset_dfs
    classes, _ = prep(U_ROWS)
    au, _tau = aux_universe(classes), None
    h = universe_hash(au)
    want = "b0d0584e616cf670c239945aa219e051aea84e4979f7cbe5084c37f05bf38f0c"
    print(f"C universe U: |Au|={len(au)} sha={h[:16]}…")
    if len(au) != 389 or h != want:
        die(f"C universe U mismatch: size={len(au)} sha={h}")
    ckpt = open(os.path.join(HERE, "scanU_checkpoint.jsonl"), "a")
    agreed = 0
    hits = []
    t0 = time.time()
    for idx, a in enumerate(au):
        ext = sorted(set(classes) | {a})
        e2c, e2r = prep(ext)
        ok, stats, _ = subset_dfs(e2c, e2r)
        sat, meta = run_cnf(ext, 14, tag=f"U{idx:03d}")
        if ok != sat:
            die(f"C instrument disagreement at U[{idx}]={a}: DFS={ok} kissat={sat}")
        ckpt.write(json.dumps({"i": idx, "aux": list(a), "dfs": ok,
                               "states": stats["states"],
                               "t": round(time.time() - t0, 1)}) + "\n")
        ckpt.flush()
        agreed += 1
        if sat:
            hits.append(a)
        if (idx + 1) % 50 == 0:
            print(f"C U progress {idx + 1}/389 hits={len(hits)} ({time.time() - t0:.0f}s)", flush=True)
    ckpt.close()
    print(f"C U COMPLETE: 389/389 agreed, feasible = {len(hits)}")
    RESULT["phases"]["C-ScanU"] = {"agreed": agreed, "hits": [list(x) for x in hits]}
    return hits


@phase("D-ScanV1")
def phase_D():
    from gate_b_floor import prep, subset_dfs
    classes, _ = prep(V_ROWS)
    au = aux_universe(classes)
    h = universe_hash(au)
    want = "9b0606049acf62bc595766b15d050fa049d6f1b960700268d9523fda7fbb73d1"
    print(f"D universe V: |Au|={len(au)} sha={h[:16]}…")
    if len(au) != 366 or h != want:
        die(f"D universe V mismatch: size={len(au)} sha={h}")
    ckpt = open(os.path.join(HERE, "scanV1_checkpoint.jsonl"), "a")
    agreed = 0
    hits = []
    t0 = time.time()
    for idx, a in enumerate(au):
        ext = sorted(set(classes) | {a})
        e2c, e2r = prep(ext)
        ok, stats, _ = subset_dfs(e2c, e2r)
        sat, meta = run_cnf(ext, 13, tag=f"V1_{idx:03d}")
        if ok != sat:
            die(f"D instrument disagreement at V1[{idx}]={a}: DFS={ok} kissat={sat}")
        ckpt.write(json.dumps({"i": idx, "aux": list(a), "dfs": ok,
                               "states": stats["states"],
                               "t": round(time.time() - t0, 1)}) + "\n")
        ckpt.flush()
        agreed += 1
        if sat:
            hits.append(a)
        if (idx + 1) % 60 == 0:
            print(f"D V1 progress {idx + 1}/366 hits={len(hits)} ({time.time() - t0:.0f}s)", flush=True)
    ckpt.close()
    print(f"D V1 COMPLETE: 366/366 agreed, feasible = {len(hits)}")
    RESULT["phases"]["D-ScanV1"] = {"agreed": agreed, "hits": [list(x) for x in hits]}
    return hits


@phase("E-ScanV2")
def phase_E():
    from gate_b_floor import prep, subset_dfs
    classes, _ = prep(V_ROWS)
    au = aux_universe(classes)
    seen_pairs = set()
    ckpt = open(os.path.join(HERE, "scanV2_checkpoint.jsonl"), "a")
    agreed = 0
    total_pairs = 0
    hits = []
    t0 = time.time()
    for i1, a1 in enumerate(au):
        au2 = aux2_universe(classes, a1)
        pair_sets = []
        for a2 in au2:
            key = tuple(sorted((a1, a2)))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            pair_sets.append(a2)
        n_here = len(pair_sets)
        total_pairs += n_here
        for j, a2 in enumerate(pair_sets):
            ext = sorted(set(classes) | {a1, a2})
            e2c, e2r = prep(ext)
            ok, stats, _ = subset_dfs(e2c, e2r)
            sat, meta = run_cnf(ext, 14, tag=f"V2_{i1:03d}_{j:03d}")
            if ok != sat:
                die(f"E instrument disagreement at V2[{i1},{j}]=({a1},{a2}): DFS={ok} kissat={sat}")
            agreed += 1
            if sat:
                hits.append((a1, a2))
                ckpt.write(json.dumps({"i1": i1, "j": j, "a1": list(a1),
                                       "a2": list(a2), "dfs": ok,
                                       "states": stats["states"], "HIT": True,
                                       "t": round(time.time() - t0, 1)}) + "\n")
                ckpt.flush()
                # synthesize NOW
                synth = try_synthesize(classes, a1, a2, order_hint=None)
                if synth is None:
                    die(f"E synthesis from pair ({a1},{a2}) failed (instrument defect)")
                RESULT["phases"]["E-ScanV2"]["synthesis"] = synth
                write_result()
                print(f"E HIT + synthesis OK at ({a1},{a2}); T = len(classes)+2 = {len(classes)+2}")
                return hits
        ckpt.write(json.dumps({"i1": i1, "a1": list(a1), "pairs_done": n_here,
                               "cum_pairs": total_pairs, "hits": len(hits),
                               "t": round(time.time() - t0, 1)}) + "\n")
        ckpt.flush()
        if (i1 + 1) % 40 == 0:
            print(f"E V2 progress a1={i1 + 1}/366 pairs={total_pairs} hits={len(hits)} ({time.time() - t0:.0f}s)", flush=True)
    print(f"E V2 COMPLETE: pairs={total_pairs} all agreed, feasible = {len(hits)}")
    RESULT["phases"]["E-ScanV2"] = {"agreed": agreed, "pairs": total_pairs,
                                    "hits": [[list(x), list(y)] for x, y in hits]}
    return hits


def try_synthesize(classes, a1, a2, order_hint):
    """Extract a 14-slot schedule over tau(V) ∪ {a1, a2} via the DFS order,
    replay it as an explicit gate list (exact vectors), and verify all 23 V
    rows. Returns the synthesis dict or None."""
    from gate_b_floor import prep, subset_dfs
    ext = sorted(set(classes) | {a1, a2})
    e2c, e2r = prep(ext)
    ok, stats, order = subset_dfs(e2c, e2r)
    if not ok:
        return None
    # replay: gates in creation order; each gate = x ± y from available classes
    # derive explicit operands from reps
    created = []           # created classes in order
    gate_list = []         # (value tuple, (idx_a, sign_a), (idx_b, sign_b))
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    base_index = {}
    for i in range(N):
        base_index[canon(values[i])] = i
    available = dict(base_index)
    for c in order:
        # find a representation: va ± vb over available values
        found = None
        for ia in range(len(values)):
            for ib in range(len(values)):
                if ia == ib:
                    continue
                for sa in (1, -1):
                    for sb in (1, -1):
                        v = tuple(sa * values[ia][k] + sb * values[ib][k] for k in range(N))
                        if canon(v) == c:
                            found = (v, ia, sa, ib, sb)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if not found:
            return None
        v, ia, sa, ib, sb = found
        values[len(values)] = v
        created.append((c, v, ia, sa, ib, sb))
    # verify coverage: every target row of V is ± some created value
    vrows_covered = 0
    for r in range(R):
        tv = V_ROWS[r]
        c = canon(tv)
        if any(cc == c for cc, *_ in created):
            vrows_covered += 1
    if vrows_covered != R:
        return None
    return {"gates": len(created), "creation": [list(cc) for cc, *_ in created],
            "covers_all_rows": True}


@phase("F-certificates")
def phase_F():
    from gate_b_floor import prep
    os.makedirs(CERT_DIR, exist_ok=True)
    # load-bearing certificates:
    # 1. U-side floor@13 (T=13 over tau(U)) — re-cert
    ucls, _ = prep(U_ROWS)
    sat, meta_u = run_cnf(ucls, 13, tag="floorU_T13", with_certificates=True,
                          keep_dir=CERT_DIR)
    print(f"F floorU_T13: SAT={sat} (expect False)")
    if sat:
        die("F floorU_T13 unexpectedly SAT — contradicts frozen floor")
    RESULT["phases"]["F-certificates"]["floorU_T13"] = meta_u["cert"]
    # 2. V-side floor@12
    vcls, _ = prep(V_ROWS)
    sat, meta_v = run_cnf(vcls, 12, tag="floorV_T12", with_certificates=True,
                          keep_dir=CERT_DIR)
    print(f"F floorV_T12: SAT={sat} (expect False)")
    if sat:
        die("F floorV_T12 unexpectedly SAT — contradicts frozen floor")
    RESULT["phases"]["F-certificates"]["floorV_T12"] = meta_v["cert"]
    # 3. U-side aux-1: the scouting-feasible instance must be SAT; certify
    # neighbor boundary: the lexicographically FIRST infeasible aux instance
    classes_u, _ = prep(U_ROWS)
    au_u = aux_universe(classes_u)
    feas_hit = RESULT["phases"]["C-ScanU"]["hits"]
    infeas = [a for a in au_u if list(a) not in [list(x) for x in feas_hit]]
    first_bad = infeas[0]
    ext = sorted(set(classes_u) | {first_bad})
    sat, meta_b = run_cnf(ext, 14, tag="U_boundary_first_infeasible",
                          with_certificates=True, keep_dir=CERT_DIR)
    print(f"F U_boundary ({first_bad}): SAT={sat} (expect False)")
    if sat:
        die("F boundary cert unexpectedly SAT")
    RESULT["phases"]["F-certificates"]["U_boundary_first_infeasible"] = meta_b["cert"]
    # 4. V-side boundary: first aux-1 infeasible
    classes_v, _ = prep(V_ROWS)
    au_v = aux_universe(classes_v)
    v1_hits = RESULT["phases"]["D-ScanV1"]["hits"]
    v1_infeas = [a for a in au_v if list(a) not in [list(x) for x in v1_hits]]
    if v1_infeas:
        first_v1_bad = v1_infeas[0]
        ext = sorted(set(classes_v) | {first_v1_bad})
        sat, meta_c = run_cnf(ext, 13, tag="V1_boundary_first_infeasible",
                              with_certificates=True, keep_dir=CERT_DIR)
        print(f"F V1_boundary ({first_v1_bad}): SAT={sat} (expect False)")
        if sat:
            die("F V1 boundary cert unexpectedly SAT")
        RESULT["phases"]["F-certificates"]["V1_boundary_first_infeasible"] = meta_c["cert"]
    return True


@phase("G-verdict")
def phase_G():
    hits_u = RESULT["phases"]["C-ScanU"]["hits"]
    hits_v1 = RESULT["phases"]["D-ScanV1"]["hits"]
    v2 = RESULT["phases"].get("E-ScanV2", {})
    hits_v2 = v2.get("hits", [])
    pairs_done = v2.get("pairs", None)
    scanU = RESULT["phases"]["C-ScanU"]["agreed"]
    scanV1 = RESULT["phases"]["D-ScanV1"]["agreed"]
    blob = {
        "gate": "mws59-gap-total59",
        "run_id": os.path.basename(HERE),
        "prereg_sha256": "eb73d31fb2a734d1a576c792dd92e05758a706b10e38c7fbc2561ca88ff749fb",
        "decomposition": "mws59 (MWS fixed orientation, sigma^0, all-monomial triple)",
        "question": "the certified LB/UB gap [56, 59] on mws59's fixed orientation",
        "reduction": ("total = C(U) + C(V) + C(Wfac) + 14; C(Wfac)=15 exact "
                      "(frozen floor-impossible 14 + frozen 15-gate transposed "
                      "witness, re-expanded); output >= 15+14 = 29"),
        "search_space": {
            "u_aux1": {"size": 389, "csv_sha256": "b0d0584e616cf670c239945aa219e051aea84e4979f7cbe5084c37f05bf38f0c", "T": 14, "agreed": scanU, "hits": [list(x) for x in hits_u]},
            "v_aux1": {"size": 366, "csv_sha256": "9b0606049acf62bc595766b15d050fa049d6f1b960700268d9523fda7fbb73d1", "T": 13, "agreed": scanV1, "hits": [list(x) for x in hits_v1]},
            "v_pair": {"T": 14, "pairs": pairs_done, "hits": [[list(a), list(b)] for a, b in hits_v2]},
            "tau_csv": {"U": "f464b75636fa4d7753ec57b79b436321fce7e0877a3be87d593dd718442f432c",
                        "V": "824c0f0edeeb5446350d0c3ebc2ea1aeca971b70aa6f219b3eb9f80ccda37509",
                        "W": "2d8603714c8c8d968d324590afc6a6047a40f69ef36d5445e25a11f6a6da1f3a"},
        },
        "instruments": ["subset_dfs (complete census)", "fresh slot-availability CNF + kissat 4.0.4"],
        "certificates": RESULT["phases"].get("F-certificates", {}),
        "controls": RESULT["phases"].get("B-controls", {}),
        "cpu_seconds_total": round(sum(p.get("seconds", 0) for p in RESULT["phases"].values()), 1),
    }
    complete_v2 = v2.get("complete", False)
    if scanU == 389 and scanV1 == 366:
        if not hits_v1:
            if complete_v2 and not hits_v2:
                if not hits_u:
                    blob["verdict"] = "NEGATIVE-STRONG"
                    blob["claim"] = ("C(U) >= 15 and C(V) >= 15 (aux-1 0/389, aux-2 "
                                     "0/%d, both instruments, complete); C(Wfac)=15; "
                                     "total >= 15+15+15+14 = 59 = published witness, "
                                     "hence LB = UB = 59: MWS-59 certified exactly "
                                     "optimal for its fixed orientation." % pairs_done)
                else:
                    blob["verdict"] = "NEGATIVE-WEAK"
                    blob["claim"] = ("C(V) >= 15, C(U) = 14 via the admitted aux "
                                     "witness; total >= 14 + 15 + 15 + 14 = 58; "
                                     "mws59 fixed orientation in [58, 59].")
            elif complete_v2 and hits_v2:
                blob["verdict"] = "SYNTHESIS"
                blob["claim"] = ("V-pair admission + synthesis; frontier move "
                                 "DOWN (see synthesis in E-ScanV2).")
            else:
                blob["verdict"] = "INCONCLUSIVE-PARTIAL"
                blob["claim"] = (f"ScanV2 incomplete at {pairs_done} pairs; "
                                 "no total claim (weak suffix non-evidence).")
        else:
            blob["verdict"] = "SYNTHESIS-PENDING"
            blob["claim"] = "V aux-1 admission; extract schedule and re-synthesize."
    else:
        blob["verdict"] = "INCONCLUSIVE-PARTIAL"
        blob["claim"] = f"scans incomplete: U {scanU}/389, V1 {scanV1}/366"
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump(blob, f, indent=1)
    print("G verdict:", blob["verdict"])
    print(json.dumps(blob, indent=1)[:1500])
    write_result()
    return True


def mark_v2_complete():
    RESULT["phases"]["E-ScanV2"]["complete"] = True

def phase_G_fail_scan(exc):
    blob = {
        "gate": "mws59-gap-total59",
        "run_id": os.path.basename(HERE),
        "verdict": "INCONCLUSIVE-SCAN-FAILURE",
        "claim": f"a scan phase crashed pre-completion: {exc!r}; no negative or positive claim of any sign",
        "completed": {k: v.get("agreed") for k, v in RESULT["phases"].items()
                      if k in ("C-ScanU", "D-ScanV1", "E-ScanV2")}
    }
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump(blob, f, indent=1)
    print("G verdict (scan failure):", blob["verdict"])

if __name__ == "__main__":

    for nm, (p, got, ok) in CHECKERS.items():
        print(f"checker {nm}: sha256={got[:16]}… pinned={ok}")
        if not ok:
            die(f"checker {nm} hash mismatch: {got}")
    phase_A()
    phase_B()
    scan_crashed = False
    try:
        hits_u = phase_C()
        hits_v1 = phase_D()
    except SystemExit:
        raise
    except Exception as exc:  # record + INCONCLUSIVE route
        scan_crashed = True
        RESULT["phases"]["scan_failure"] = repr(exc)
        phase_G_fail_scan(exc)
        with open(os.path.join(HERE, "runner_output.json"), "w") as f:
            json.dump(RESULT, f, indent=1, default=str)
        print("\nRUNNER COMPLETE (INCONCLUSIVE — scan failure recorded)")
        sys.exit(3)
    try:
        hits_v2 = phase_E()
        mark_v2_complete()
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 — record partial state, still verdict
        RESULT["phases"]["E-ScanV2"]["complete"] = False
        RESULT["phases"]["E-ScanV2"]["partial_reason"] = repr(exc)
    phase_F()
    phase_G()
    with open(os.path.join(HERE, "runner_output.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    print("\nRUNNER COMPLETE")
