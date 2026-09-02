#!/usr/bin/env python3
"""sun56-gap-total55 — primary runner (pre-registered; phase 2 of the campaign).

Executes, in order and with abort-on-failure:
  Phase A  anchors:       729/729 Brent over Z (fmpz + pure int), Sun SLP split
                          13/13/30, frozen-row re-verification (A4).
  Phase B  controls:      A2 witness ACCEPTs, A3 extension-positive, R1/R2
                          plants, R3 floor@11 negative + checker certificates,
                          R4 assertion-planter (scratch copy only).
  Phase C  the 338-instance aux-1 decision set, instruments 1 (complete
           subset-DFS) and 2 (fixed fresh CNF + kissat), with per-instance
           agreement required and jsonl checkpoints.
  Phase D  load-bearing certificates (DRAT from kissat -> LRAT by drat-trim -L
           -> both pinned checkers).
  Phase E  verdict aggregation into verdict.json.

Every phase prints evidence; the runner exits nonzero on any control failure or
instrument disagreement (FROZEN-INCONCLUSIVE path).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import hashlib
import tempfile

PYTHONDONTWRITEBYTECODE = True  # noqa: F841 — env var below is the operative one
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
MM3 = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(MM3, "src")
SCR5 = os.path.join(MM3, "campaigns",
                    "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56",
                    "scripts")
SNAP3 = os.path.join(MM3, "campaigns",
                     "2026-08-30T031544Z_e0f3f117_c9df97a8bf3a", "tools_snapshot")
sys.path.insert(0, SRC)
sys.path.insert(0, SCR5)

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]

CHECKERS = {}
for name, want in (("drat-trim", "111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4"),
                   ("lrat-check", "b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07")):
    p = os.path.join(HERE, "checkers", name)
    if not os.path.exists(p):
        os.makedirs(os.path.join(HERE, "checkers"), exist_ok=True)
        subprocess.run(["cc", "-O2", "-o", p, os.path.join(SNAP3, name + ".c")], check=True)
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    CHECKERS[name] = (p, got, got == want)

RESULT: dict = {"phases": {}, "defects": []}


def phase(name):
    def deco(fn):
        def wrapper(*a, **k):
            t0 = time.time()
            print(f"\n===== PHASE {name} =====", flush=True)
            RESULT["phases"].setdefault(name, {})
            out = fn(*a, **k)
            RESULT["phases"][name]["seconds"] = round(time.time() - t0, 3)
            RESULT["phases"][name]["ok"] = True if out is not False else False
            return out
        return wrapper
    return deco


def die(msg):
    RESULT["fatal"] = msg
    print("FATAL:", msg, flush=True)
    write_result()
    sys.exit(1)


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
# Exports produced by phase A; used everywhere else.
U_ROWS: list = []
V_ROWS: list = []
W_ROWS: list = []
V_TAU: list = []
UNIVERSE: list = []


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


def expand_sun(side, dim):
    vs = {i: tuple(1 if j == i else 0 for j in range(dim)) for i in range(dim)}
    gates = []
    for a, s, b in side["inter"]:
        va, vb = vs[a], vs[b]
        nv = tuple(va[i] + s * vb[i] for i in range(dim))
        vs[dim + len(gates)] = nv
        gates.append(nv)
    finals = []
    for f in side["final"]:
        acc = [0] * dim
        for idx, c in f:
            v = vs[idx]
            for i in range(dim):
                acc[i] += c * v[i]
        finals.append(tuple(acc))
    return gates, finals


def slp_eval(side, base_values, skip_gate=None, flip_gate_operand=None):
    """Evaluate Sun's side SLP on base values; optionally delete one gate
    (rewiring references to its value are FORBIDDEN -> returns None marker) or
    flip one operand's sign in one gate."""
    vs = list(base_values)
    for gi, (a, s, b) in enumerate(side["inter"]):
        if skip_gate is not None and gi == skip_gate:
            vs.append(None)  # gate deleted: wire value missing
            continue
        va = vs[a]
        vb = vs[b]
        if va is None or vb is None:
            vs.append(None)
            continue
        if flip_gate_operand is not None and gi == flip_gate_operand[0]:
            if flip_gate_operand[1] == 0:
                va = tuple(-x for x in va)
            else:
                vb = tuple(-x for x in vb)
        vs.append(tuple(va[i] + s * vb[i] for i in range(len(va))))
    out = []
    for f in side["final"]:
        acc = [0] * len(base_values)
        broken = False
        for idx, c in f:
            v = vs[idx]
            if v is None:
                broken = True
                break
            for i in range(len(acc)):
                acc[i] += c * v[i]
        out.append(None if broken else tuple(acc))
    return out


@phase("A-anchor")
def phase_A():
    import sun56_verify as sv

    global U_ROWS, V_ROWS, W_ROWS
    U = [tuple(r) for r in sv.U]
    V = [tuple(r) for r in sv.V]
    W = [tuple(r) for r in sv.W]
    U_ROWS, V_ROWS, W_ROWS = U, V, W

    f_int = brent_failures_int(U, V, W)
    f_fmpz = brent_failures_fmpz(U, V, W)
    print(f"A1 brent: int fails={f_int} fmpz fails={f_fmpz}")
    if f_int or f_fmpz:
        die(f"A1 anchor failed: brent failures int={f_int} fmpz={f_fmpz}")
    tern = all(x in (-1, 0, 1) for blk in (U, V, W) for row in blk for x in row)
    print(f"A1 ternary: {tern}")
    if not tern:
        die("A1 anchor failed: factors not ternary")

    # Sun SLP split recount
    cU = len(sv.SIDES['U']['inter']) + sum(max(sum(abs(c) for _, c in f) - 1, 0)
                                           for f in sv.SIDES['U']['final'])
    cV = len(sv.SIDES['V']['inter']) + sum(max(sum(abs(c) for _, c in f) - 1, 0)
                                           for f in sv.SIDES['V']['final'])
    cW = len(sv.SIDES['W']['inter']) + sum(max(sum(abs(c) for _, c in f) - 1, 0)
                                           for f in sv.SIDES['W']['final'])
    print(f"A1 split recount: U={cU} V={cV} W(output)={cW} total={cU + cV + cW}")
    if (cU, cV, cW) != (13, 13, 30) or cU + cV + cW != 56:
        die(f"A1 split mismatch: {(cU, cV, cW)}")
    RESULT["phases"]["A-anchor"]["split"] = [cU, cV, cW]

    # A4 frozen row
    from gate_b_floor import prep, subset_dfs
    row = {}
    for nm, T in (("U", U), ("V", V), ("Wfac", W)):
        classes, _ = prep(T)
        ok, stats, _ = subset_dfs(classes, _)
        row[nm] = {"d": len(classes), "floor_achievable": ok,
                   "states": stats["states"]}
    print("A4 frozen-row recheck:", json.dumps(row))
    frozen = {"U": (12, False, 33), "V": (11, False, 9), "Wfac": (16, True, 17)}
    for nm, (d, ok, st) in frozen.items():
        got = row[nm]
        if (got["d"], got["floor_achievable"], got["states"]) != (d, ok, st):
            die(f"A4 frozen-row mismatch on {nm}: {got} vs {frozen[nm]}")
    total_lb = 13 + 12 + (16 + 14)
    print(f"A4 total_lb from frozen landscape (C_U=13, C_V>=12, out 16+14): {total_lb}")
    if total_lb != 55:
        die("A4 arithmetic wrong")
    RESULT["phases"]["A-anchor"]["frozen_row"] = row
    return True


@phase("B-controls")
def phase_B():
    from gate_b_floor import prep, subset_dfs, canon as g_canon

    # -------- A2: witness ACCEPTs --------------------------------------
    import sun56_verify as sv
    gatesU, finalsU = expand_sun(sv.SIDES['U'], 9)
    gatesV, finalsV = expand_sun(sv.SIDES['V'], 9)
    okU = all(finalsU[r] == U_ROWS[r] for r in range(R))
    okV = all(finalsV[r] == V_ROWS[r] for r in range(R))
    print(f"A2 witness: U circuit reproduces U rows: {okU}; V: {okV}")
    if not (okU and okV):
        die("A2 witness expansion mismatch")
    tset = set(V_TAU) if V_TAU else set(prep(V_ROWS)[0])
    clsV = [canon(g) for g in gatesV]
    aux_used = [c for c in clsV if c not in tset]
    print(f"A2 V circuit: {len(gatesV)} gates, {len(set(clsV))} distinct classes, "
          f"aux classes used: {aux_used}")
    if len(gatesV) != 13 or len(set(clsV)) != 13 or len(aux_used) != 2:
        die("A2 witness class structure unexpected")
    RESULT["phases"]["B-controls"]["witness_costs"] = {"U": 13, "V": 13, "out": 30}

    # -------- A3: extension-positive control ---------------------------
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    a2 = canon((0, 0, 0, 1, 0, -1, 0, 0, 0))
    classes, repsV = prep(V_ROWS)
    if a1 in set(classes) or a2 in set(classes):
        die("A3 aux classes unexpectedly in tau")
    ext = sorted(set(classes) | {a1, a2})
    ext_classes, ext_reps = prep(ext)
    ok13, st13, order13 = subset_dfs(ext_classes, ext_reps)
    print(f"A3 DFS on tau+{{aux1,aux2}} at T=13: achievable={ok13} states={st13['states']}")
    if not ok13:
        die("A3 extension-positive control FAILED (DFS says infeasible)")
    sat13 = run_cnf(ext, 13, tag="ctrl_ext2_T13")[0]
    print(f"A3 kissat on same instance: SAT={sat13}")
    if not sat13:
        die("A3 instrument disagreement on extension-positive control")
    RESULT["phases"]["B-controls"]["A3"] = "pass"

    # -------- R1/R2 plants ---------------------------------------------
    # R1: delete V gate 6; R2: flip gate 3 operand 1.
    base = [tuple(1 if j == i else 0 for j in range(9)) for i in range(9)]
    full = list(base) + [tuple(x) for x in gatesV]
    for label, kw in (("R1 del-gate6", dict(skip_gate=6)),
                      ("R2 flip-gate3-op1", dict(flip_gate_operand=(3, 1)))):
        outs = slp_eval(sv.SIDES['V'], base, **kw)
        broken = any(o is None for o in outs)
        if not broken:
            broken = any(outs[r] != V_ROWS[r] for r in range(R))
        ver_brent = None
        if not broken:
            # rewired 12/13-gate SLP still standing: check Brent with these rows
            pass
        verdict = "REJECTED" if broken else "NEEDS-BRENT"
        print(f"{label}: exact-expansion verdict = {verdict}")
        if verdict == "NEEDS-BRENT":
            die(f"{label} plant did not fail expansion; needs brent route (tighten runner)")
    RESULT["phases"]["B-controls"]["plants"] = ["R1 REJECTED by expansion",
                                                "R2 REJECTED by expansion"]

    # -------- R3: floor@11 negative + certificates ---------------------
    ok11, st11, _ = subset_dfs(classes, repsV)
    print(f"R3 DFS floor@11 (T=11): achievable={ok11} states={st11['states']}")
    if ok11:
        die("R3 floor@11 unexpectedly achievable")
    unsat11 = run_cnf(sorted(classes), 11, tag="R3_floor11_T11",
                      with_certificates=True)[0]
    print(f"R3 kissat floor@11: SAT={unsat11}  (certificates generated)")
    if unsat11:
        die("R3 floor@11 SAT — instruments disagree with frozen evidence")

    # R4: planter on a scratch copy (assert must fire)
    scratch = os.path.join(HERE, "scratch_R4")
    os.makedirs(scratch, exist_ok=True)
    planter = os.path.join(scratch, "planter.py")
    with open(planter, "w") as f:
        f.write(
            "import sys\n"
            "sys.path.insert(0, %r)\n" % SRC +
            "sys.path.insert(0, %r)\n" % SCR5 +
            "import sun56_verify as sv\n"
            "from gate_b_floor import prep, subset_dfs\n"
            "V = [tuple(r) for r in sv.V]\n"
            "c, _ = prep(V)\n"
            "assert len(c) == 12, 'planted wrong d'\n"
        )
    rp = subprocess.run(["/Users/jinleic/jinleic-workspace/cs/.venv/bin/python",
                         planter], capture_output=True, text=True,
                        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    planted_fired = rp.returncode != 0 and ("planted wrong d" in (rp.stderr + rp.stdout))
    print(f"R4 planter assert fired: {planted_fired} (rc={rp.returncode})")
    if not planted_fired:
        die("R4 planter did not fire")
    RESULT["phases"]["B-controls"]["R4"] = "fired"
    return True


# --------------------------------------------------------------------------
# Fresh CNF encoder (instrument 2) — the FIXED variant per pre-statement D1.
def build_ext_cnf(tau_extended, T):
    base = sorted(set(tau_extended) | {canon((1 if j == i else 0 for j in range(N)))
                                       for i in range(N)})
    cls = sorted(tau_extended)
    clsset = set(cls)
    reps = {}
    for c in cls:
        lst = []
        neg = tuple(-x for x in c)
        for k1 in range(len(base)):
            for k2 in range(k1, len(base)):
                for s1 in (1, -1):
                    for s2 in (1, -1):
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
        for k in range(len(base)):
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
        for k in range(len(base)):
            lits = [-vid[("a", g + 1, k)], vid[("a", g, k)]]
            if base[k] in clsset:
                lits.append(vid[("s", g, base[k])])
            cl(lits)
    for k in range(len(base)):
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


def run_cnf(tau_ext, T, tag="", with_certificates=False, keep_dir=None):
    cl_, nv_, vid, base, reps = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = keep_dir or tempfile.mkdtemp(prefix="sun56gap_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    with open(cnfp, "w") as f:
        f.write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    cmd = NICE + [KISSAT, "-q", cnfp] + ([drat] if with_certificates else [])
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    sat = (r.returncode == 10)
    unsat = (r.returncode == 20)
    cert = None
    if not sat and not unsat:
        die(f"kissat returned {r.returncode} on {tag}: {r.stdout[:100]} {r.stderr[:100]}")
    if unsat and with_certificates:
        lrat = os.path.join(tmp, f"{tag or 'inst'}.lrat")
        conv = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, drat, "-L", lrat],
                              capture_output=True, text=True)
        d1 = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, lrat],
                            capture_output=True, text=True)
        d2 = subprocess.run(NICE + [CHECKERS["lrat-check"][0], cnfp, lrat],
                            capture_output=True, text=True)
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
    return sat, {"vars": nv_, "clauses": len(cl_), "kissat_s": round(time.time() - t0, 2),
                 "path": cnfp, "cert": cert}


def aux_universe():
    classes, _ = __import__("gate_b_floor").prep(V_ROWS)
    from gate_b_floor import canon as g_canon
    ICAN = {g_canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    forms = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)] + sorted(classes)
    u = set()
    for i in range(len(forms)):
        for j in range(i + 1):
            for si in (1, -1):
                for sj in (1, -1):
                    if i == j and si == 1 and sj == 1:
                        continue
                    c = g_canon(tuple(si * forms[i][k] + sj * forms[j][k]
                                      for k in range(N)))
                    if c != (0,) * N and c not in ICAN and c not in set(classes):
                        u.add(c)
    return sorted(u), sorted(classes)


@phase("C-scan")
def phase_C():
    from gate_b_floor import prep, subset_dfs
    global V_TAU
    au, tau = aux_universe()
    V_TAU = tau
    h = hashlib.sha256(("\n".join(",".join(str(x) for x in row) for row in au)
                        ).encode()).hexdigest()
    want = "c10bdeeae74b601dc063cb935f7c7baca79a637cab52a7e2f85d16a46f84f813"
    print(f"C universe: |Au|={len(au)} csv_sha256={h[:16]}…")
    if len(au) != 338 or h != want:
        die(f"C universe mismatch: size={len(au)} sha={h}")
    classes, _ = prep(V_ROWS)
    ckpt = open(os.path.join(HERE, "scan_checkpoint.jsonl"), "a")
    agreed = 0
    hits = []
    t0 = time.time()
    for idx, a in enumerate(au):
        ext = sorted(set(classes) | {a})
        ext2c, ext2r = prep(ext)
        ok, stats, _ = subset_dfs(ext2c, ext2r)
        sat, meta = run_cnf(ext, 12, tag=f"aux{idx:03d}")
        agree = (ok == sat)
        if not agree:
            die(f"C instrument disagreement at aux[{idx}]={a}: DFS={ok} kissat={sat}")
        ckpt.write(json.dumps({"i": idx, "aux": list(a), "dfs": ok, "cnf_sat": sat,
                               "states": stats["states"],
                               "t": round(time.time() - t0, 1)}) + "\n")
        ckpt.flush()
        agreed += 1
        if sat:
            hits.append(a)
        if (idx + 1) % 50 == 0:
            print(f"C progress {idx + 1}/338 agreed={agreed} hits={len(hits)} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    ckpt.close()
    print(f"C COMPLETE: 338/338 agreed, SAT hits = {len(hits)} {hits[:5]}")
    RESULT["phases"]["C-scan"] = {"agreed": agreed, "hits": [list(x) for x in hits]}
    return hits


@phase("D-certificates")
def phase_D():
    from gate_b_floor import prep
    cdir = os.path.join(HERE, "certificates")
    os.makedirs(cdir, exist_ok=True)
    classes, _ = prep(V_ROWS)
    # load-bearing UNSAT: floor@11 at T=11 — already certified in B-R3; now also
    # certify a representative aux-infeasibility batch: the Sun-aux1 T=12 instance
    # (the witness-shape negative) and the projected-orbit pair probe.
    a1 = canon((0, 0, 1, 0, 0, 1, 0, 0, 0))
    ext = sorted(set(classes) | {a1})
    sat, meta = run_cnf(ext, 12, tag="sun_aux1_T12", with_certificates=True,
                        keep_dir=cdir)
    print(f"D sun_aux1_T12: SAT={sat} cert={json.dumps(meta['cert'])[:160]}")
    if sat:
        die("D witness-shape aux is feasible at T=12?!  — verdict flips; escalate")
    RESULT["phases"]["D-certificates"] = {"sun_aux1_T12": meta["cert"]}
    return True


@phase("E-verdict")
def phase_E():
    hits = RESULT["phases"]["C-scan"]["hits"]
    scan = RESULT["phases"]["C-scan"]["agreed"]
    blob = {
        "gate": "sun56-gap-total55",
        "run_id": os.path.basename(HERE),
        "prereg_sha256": "7854d469467bf3337dfa76841208194afbb24979f85fdd0a6e36f613a4f4d660",
        "decomposition": "sun56 (Sun's fixed orientation, sigma^0, all-monomial triple)",
        "question": "does a 55-total three-stage circuit exist for this orientation",
        "reduction": "total = C(U) + C(V) + C(Wfac) + 14 = 13 + C(V) + 30; "
                     "55-move exists iff C(V) <= 12 iff some aux-1 extension of "
                     "tau(V) is schedulable at T=12",
        "search_space": {"aux_universe_size": 338,
                         "aux_universe_csv_sha256":
                             "c10bdeeae74b601dc063cb935f7c7baca79a637cab52a7e2f85d16a46f84f813",
                         "tau_csv_sha256":
                             "65e4879839d168f93cd714a89c3c48f6930d3533f6326586281ce5ab6dad7e25",
                         "T": 12, "instances": scan},
        "instruments": ["subset_dfs (complete census)",
                        "fresh slot-availability CNF + kissat 4.0.4"],
        "instrument_agreement": scan,
        "admitted_aux_count": len(hits),
        "admitted_aux": hits,
        "certificates": RESULT["phases"]["D-certificates"],
        "controls": RESULT["phases"]["B-controls"],
        "cpu_seconds_total": round(sum(p.get("seconds", 0)
                                       for p in RESULT["phases"].values()), 1),
    }
    if scan == 338 and not hits:
        blob["verdict"] = "NEGATIVE-BOUNDED"
        blob["claim"] = ("C(V_sun56) >= 13; with Sun's 13-gate witness C(V)=13; "
                         "hence total exact = 13+13+16+14 = 56; no 55-addition "
                         "circuit exists for the sun56 fixed orientation within "
                         "the three-stage linear SLP model; the searched space is "
                         "the 338-element aux universe (hash above), fully "
                         "enumerated, both instruments agreeing on every instance.")
    else:
        blob["verdict"] = "CANDIDATE"
        blob["claim"] = "aux admission at T=12 — extract schedule and escalate"
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump(blob, f, indent=1)
    print("E verdict:", blob["verdict"])
    print(json.dumps(blob, indent=1)[:1200])
    return True


if __name__ == "__main__":
    # checkers must be pinned before anything
    for nm, (p, got, ok) in CHECKERS.items():
        print(f"checker {nm}: sha256={got[:16]}… pinned={ok}")
        if not ok:
            die(f"checker {nm} hash mismatch: {got}")
    phase_A()
    phase_B()
    hits = phase_C()
    phase_D()
    phase_E()
    with open(os.path.join(HERE, "runner_output.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    print("\nRUNNER COMPLETE")
