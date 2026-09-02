#!/usr/bin/env python3
"""laderman23-ladder-total - primary runner (pre-registered; prereg sha 2af3b859...).

Implements pre_statement.md phases A..H exactly:
  A anchors -> B controls -> C universes -> D aux-1 censuses (dual instrument)
  -> E synthesis+verification -> F assembly+transposition+end-to-end
  -> G certificates -> H verdict.
Aborts (exit 1) on any control failure or instrument disagreement.
"""
import os
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import random
import subprocess
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MM3 = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(MM3, "src")
FROZEN_SCRIPTS = os.path.join(
    MM3, "campaigns",
    "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56",
    "scripts")
PHASE1 = os.path.join(MM3, "campaigns", "20260902T022208Z_8c279ad1_d0276dfee63d")
sys.path.insert(0, SRC)

KISSAT = "/opt/homebrew/bin/kissat"
CADICAL = "/opt/homebrew/bin/cadical"
NICE = ["nice", "-n", "15"]
ENV1 = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
            NUMEXPR_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
CERT_DIR = os.path.join(HERE, "certificates")
os.makedirs(CERT_DIR, exist_ok=True)
CKPT = os.path.join(HERE, "scan_checkpoint.jsonl")
CHECKERS = {
    "drat-trim": (os.path.join(HERE, "tools", "drat-trim"),
                  "111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4"),
    "lrat-check": (os.path.join(HERE, "tools", "lrat-check"),
                   "b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07"),
}
# pre-registered universe pins (prereg section 5)
AU_PINS = {"U": ("a1de3df4e5077746", 408), "V": ("d3f3894227f939fe", 408),
           "WFac": ("13be8045ff2aec8a", 408)}
R, N = 23, 9
T_LEVEL = 15                      # the decided level: d + 1
RESTARTS = 240
MASTER_SEED = "20260902"
RESULT = {"gate": "laderman23-ladder-total", "run_id": os.path.basename(HERE),
          "prereg_sha256": "2af3b859eff465c3146cda7a9b7cc4ef883137267cda7eced132449fdf4ca4c2",
          "phases": {}}
CERT_KEEP_PER_SIDE = 5            # retention policy, disclosed in the report


def die(msg):
    print("FATAL:", msg)
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump({"gate": RESULT["gate"], "run_id": RESULT["run_id"],
                   "verdict": "INCONCLUSIVE", "reason": msg,
                   "phases": RESULT["phases"]}, f, indent=1, default=str)
    sys.exit(1)


def phase(name):
    def deco(fn):
        def wrapped(*a, **k):
            print(f"== phase {name} ==")
            t0 = time.time()
            out = fn(*a, **k)
            d = RESULT["phases"].setdefault(name, {})
            if isinstance(d, dict):
                d["seconds"] = round(time.time() - t0, 2)
            return out
        return wrapped
    return deco


from gate_b_floor import prep, subset_dfs, canon           # noqa: E402

ICAN = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}


# ================================ arithmetic ================================
def brent_int(U, V, W):
    f = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ip in range(3):
                    for jp in range(3):
                        for kp in range(3):
                            s = 0
                            for r in range(R):
                                ur = U[r][3 * i + k]
                                if not ur:
                                    continue
                                vr = V[r][3 * kp + j]
                                if not vr:
                                    continue
                                wr = W[r][3 * ip + jp]
                                if wr:
                                    s += ur * vr * wr
                            if s != (1 if (i == ip and j == jp and k == kp) else 0):
                                f += 1
    return f


def brent_fmpz(U, V, W):
    from flint import fmpz
    f = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ip in range(3):
                    for jp in range(3):
                        for kp in range(3):
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = (U[r][3 * i + k], V[r][3 * kp + j],
                                              W[r][3 * ip + jp])
                                if ur and vr and wr:
                                    s += fmpz(ur) * fmpz(vr) * fmpz(wr)
                            if int(s) != (1 if (i == ip and j == jp and k == kp) else 0):
                                f += 1
    return f


def rank_over_Q(rows):
    """Exact fraction-free Gaussian elimination rank."""
    from fractions import Fraction
    M = [[Fraction(x) for x in row] for row in rows]
    rank, piv = 0, 0
    for col in range(len(M[0])):
        sel = None
        for r2 in range(rank, len(M)):
            if M[r2][col]:
                sel = r2
                break
        if sel is None:
            continue
        M[rank], M[sel] = M[sel], M[rank]
        pv = M[rank][col]
        for r2 in range(len(M)):
            if r2 != rank and M[r2][col]:
                fac = M[r2][col] / pv
                M[r2] = [a - fac * b for a, b in zip(M[r2], M[rank])]
        rank += 1
        piv += 1
        if rank == len(M):
            break
    return rank


# ============================ auxiliary universe ============================
def au_universe(base_classes, second_loop=False):
    """Closure rule of prereg section 4: signed pair-sums of (inputs u tau),
    excluding input directions and tau itself. Two independent loop shapes."""
    forms = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)] + \
            [tuple(c) for c in base_classes]
    known = set(base_classes)
    u = set()
    if not second_loop:
        for i in range(len(forms)):
            for j in range(i + 1, len(forms)):
                for si in (1, -1):
                    for sj in (1, -1):
                        c = canon(tuple(si * forms[i][k] + sj * forms[j][k]
                                        for k in range(N)))
                        if c != (0,) * N and c not in ICAN and c not in known:
                            u.add(c)
    else:
        # independent shape: iterate signs outermost, pairs inner, j<i
        for si in (1, -1):
            for sj in (1, -1):
                for j in range(len(forms)):
                    for i in range(j):
                        v = tuple(si * forms[i][k] + sj * forms[j][k]
                                  for k in range(N))
                        c = canon(v)
                        if c == (0,) * N or c in ICAN or c in known:
                            continue
                        u.add(c)
    return sorted(u)


def au_sha(au):
    return hashlib.sha256("\n".join(",".join(map(str, x)) for x in au).encode()
                          ).hexdigest()


# =============================== CNF instrument =============================
def build_ext_cnf(tau_extended, T):
    """Slot-availability CNF, prereg section 6 instrument-2 spec.
    Vars: s[g][c] slot g creates class c; a[g][k] base class k available for
    gate g; r[g][c][i] representative choice. Slot-0 availability: BOTH
    direction units (inputs TRUE, needed classes FALSE)."""
    base = sorted(set(tau_extended) | ICAN)
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
        clauses.append([vid[("a", 0, k)]] if base[k] not in clsset
                       else [-vid[("a", 0, k)]])
    for c in cls:
        clauses.append([vid[("s", gg, c)] for gg in range(T)])
    for g in range(T):
        for c in cls:
            lits = [-vid[("s", g, c)]]
            for ri, (k1, k2) in enumerate(reps[c]):
                lits.append(vid[("r", g, c, ri)])
                clauses.append([-vid[("r", g, c, ri)], vid[("a", g, k1)]])
                clauses.append([-vid[("r", g, c, ri)], vid[("a", g, k2)]])
            clauses.append(lits)
    return clauses, nv


def run_cnf(tau_ext, T, tag="", certificates=False, keep=False,
            solver=KISSAT):
    cl_, nv_ = build_ext_cnf(tau_ext, T)
    dimacs = [f"p cnf {nv_} {len(cl_)}"] + \
             [" ".join(map(str, c)) + " 0" for c in cl_]
    tmp = CERT_DIR if keep else tempfile.mkdtemp(prefix="ladlad_")
    cnfp = os.path.join(tmp, f"{tag or 'inst'}.cnf")
    with open(cnfp, "w") as f:
        f.write("\n".join(dimacs) + "\n")
    drat = os.path.join(tmp, f"{tag or 'inst'}.drat")
    cmd = NICE + [solver, "-q", cnfp] + ([drat] if certificates else [])
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900,
                       env=ENV1)
    sat, unsat = (r.returncode == 10), (r.returncode == 20)
    if not sat and not unsat:
        die(f"{solver} returned {r.returncode} on {tag}")
    cert = None
    if unsat and certificates:
        lrat = os.path.join(tmp, f"{tag or 'inst'}.lrat")
        conv = subprocess.run(
            NICE + [CHECKERS["drat-trim"][0], cnfp, drat, "-L", lrat],
            capture_output=True, text=True, env=ENV1)
        d1 = subprocess.run(NICE + [CHECKERS["drat-trim"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        d2 = subprocess.run(NICE + [CHECKERS["lrat-check"][0], cnfp, lrat],
                            capture_output=True, text=True, env=ENV1)
        cert = {"cnf_sha256": hashlib.sha256(open(cnfp, "rb").read()).hexdigest(),
                "lrat_sha256": hashlib.sha256(open(lrat, "rb").read()).hexdigest(),
                "conv_rc": conv.returncode, "drattrim_rc": d1.returncode,
                "lratcheck_rc": d2.returncode}
        if d1.returncode != 0 or d2.returncode != 0:
            die(f"certificate verification FAILED for {tag}: {cert}")
        if not keep:
            for p in (cnfp, drat, lrat):
                if os.path.exists(p):
                    os.unlink(p)
        else:
            os.replace(drat, os.path.join(tmp, f"{tag}.used.drat"))
    elif not keep and os.path.exists(cnfp):
        os.unlink(cnfp)
        if os.path.exists(drat):
            os.unlink(drat)
    return sat, {"vars": nv_, "clauses": len(cl_),
                 "s": round(time.time() - t0, 3), "cert": cert}


# ============================ synthesis + verifier ==========================
def greedy_synth(rows, side, restart):
    """Randomized greedy CSE (prereg section 9). Returns an explicit gate list;
    proposes only - verification is separate and shares no code."""
    rng = random.Random(f"{MASTER_SEED}|{side}|{restart}")
    jitter = 0.0 if restart == 0 else 2.0
    avail = []                                   # (name, vector)
    for i in range(N):
        avail.append((i, tuple(1 if j == i else 0 for j in range(N))))
    have = {canon(v): nm for nm, v in avail}
    targets = [tuple(r) for r in rows]
    gates = []
    while not all(canon(t) in have for t in targets):
        best, bestv = None, None
        for i in range(len(avail)):
            for j in range(i + 1, len(avail)):
                (ni, vi), (nj, vj) = avail[i], avail[j]
                for sb in (1, -1):
                    v = tuple(vi[k] + sb * vj[k] for k in range(N))
                    if all(x == 0 for x in v):
                        continue
                    if max(abs(x) for x in v) > 1:
                        continue
                    c = canon(v)
                    if c in have:
                        continue
                    score, done = 0, 0
                    for t in targets:
                        if canon(t) in have:
                            continue
                        for s in (1, -1):
                            if all((s * v[k] == 0) or (s * v[k] == t[k])
                                   for k in range(N)):
                                nt = sum(1 for k in range(N) if s * v[k])
                                score += nt - 1
                                if nt == sum(1 for x in t if x):
                                    done += 1
                                break
                    sc = score + 3 * done + rng.random() * jitter
                    if best is None or sc > best:
                        best, bestv = sc, (ni, 1, nj, sb, v, c)
        if bestv is None:
            return None
        ni, si, nj, sj, v, c = bestv
        nm = 100 + len(gates)
        gates.append({"out": nm, "a": [ni, si], "b": [nj, sj], "val": list(v)})
        avail.append((nm, v))
        have[c] = nm
        if len(gates) > 60:
            return None
    return gates


def eval_circuit(gates):
    """Exact expansion, independent of the search: recompute every gate from
    its operand structure over 9 symbolic free inputs."""
    values = {i: tuple(1 if j == i else 0 for j in range(N)) for i in range(N)}
    for g in gates:
        (na, sa), (nb, sb) = g["a"], g["b"]
        if na not in values or nb not in values:
            return None
        v = tuple(sa * values[na][k] + sb * values[nb][k] for k in range(N))
        if list(v) != list(g["val"]):
            return None
        values[g["out"]] = v
    return values


def circuit_verdict(gates, targets):
    values = eval_circuit(gates)
    if values is None:
        return False, [("recompute-mismatch", None)]
    avail = {canon(v) for v in values.values()}
    missing = [r for r, tv in enumerate(targets) if canon(tv) not in avail]
    return (not missing), missing


def transpose_circuit(gates, targets):
    """Reverse-mode transposition to an explicit 23->9 output stage; additions
    counted exactly (each accumulation into a non-empty adjoint costs 1)."""
    values = eval_circuit(gates)
    if values is None:
        return None
    reads = []
    for r, tv in enumerate(targets):
        c = canon(tv)
        hit = None
        for nm, v in values.items():
            if canon(v) == c:
                hit = (nm, 1 if list(v) == list(tv) else -1)
                break
        if hit is None:
            return None
        reads.append(hit)
    adj = {nm: [0] * R for nm in values}
    nonempty = set()
    adds = 0
    for r, (nm, sg) in enumerate(reads):
        if nm in nonempty:
            adds += 1
        adj[nm][r] += sg
        nonempty.add(nm)
    for g in reversed(gates):
        out = g["out"]
        if out not in nonempty:
            continue
        for (nn, ss) in (g["a"], g["b"]):
            if nn in nonempty:
                adds += 1
            else:
                nonempty.add(nn)
            for j in range(R):
                if adj[out][j]:
                    adj[nn][j] += ss * adj[out][j]
    outs = [tuple(adj[i]) for i in range(N)]
    return adds, outs, {"reads": reads}


# ================================== phases ==================================
@phase("A-anchors")
def phase_A():
    for nm, (p, want) in CHECKERS.items():
        got = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if got != want:
            die(f"checker {nm} hash mismatch: {got}")
    print("checkers pinned OK (drat-trim, lrat-check)")
    d = json.load(open(os.path.join(PHASE1, "laderman23_factors.json")))
    U = [tuple(r) for r in d["U"]]
    V = [tuple(r) for r in d["V"]]
    W = [tuple(r) for r in d["W"]]
    tri_sha = hashlib.sha256(json.dumps({"U": d["U"], "V": d["V"], "W": d["W"]},
                                        sort_keys=True).encode()).hexdigest()
    if not tri_sha.startswith("522ba07f4f1784ac"):
        die(f"A1: frozen triple hash mismatch {tri_sha[:20]}")
    fi, ff = brent_int(U, V, W), brent_fmpz(U, V, W)
    if fi or ff:
        die(f"A1: frozen Laderman triple fails Brent int={fi} fmpz={ff}")
    print(f"A1 triple sha {tri_sha[:16]} Brent int={fi} fmpz={ff} -> 729/729")
    # d-counts and floor verdicts, re-derived
    sides = {}
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        c, rp = prep(M)
        ok, st, _ = subset_dfs(c, rp)
        sides[nm] = {"d": len(c), "floor_achievable": ok, "states": st["states"]}
        print(f"A1 {nm}: d={len(c)} floor_achievable={ok} states={st['states']}")
        if ok:
            die(f"A1: {nm} floor achievable - the scouted landscape changed; "
                f"the pre-registered question (T={T_LEVEL}) does not apply")
    # A5 preconditions
    pre = {}
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        nz = sum(1 for r in range(R) if any(M[r]))
        pre[nm] = {"nonzero_rows": nz}
        if nz != R:
            die(f"A5: {nm} has {nz}/23 nonzero rows")
    rk = rank_over_Q([list(r) for r in W])
    pre["rank_WFac"] = rk
    if rk != 9:
        die(f"A5: rank(WFac) = {rk} != 9")
    print(f"A5 preconditions: 23/23 nonzero rows per side, rank(WFac)={rk}")
    # printed basic-form recount (source recount, no model claim)
    prods = d["printed_products"]
    outs = d["printed_outputs"]
    ladds = sum(len(prods[str(m)][0]) - 1 for m in range(1, 24))
    radds = sum(len(prods[str(m)][1]) - 1 for m in range(1, 24))
    oadds = sum(len(v) - 1 for v in outs.values())
    printed = {"left": ladds, "right": radds, "output": oadds,
               "total": ladds + radds + oadds}
    print(f"A1 printed basic-form recount: left={ladds} right={radds} "
          f"output={oadds} total={printed['total']}")
    RESULT["phases"]["A-anchors"] = {"triple_sha256": tri_sha,
                                     "brent": {"int": fi, "fmpz": ff},
                                     "sides": sides, "preconditions": pre,
                                     "printed_basic_form": printed}
    return U, V, W, printed


@phase("B-controls")
def phase_B(U, V, W):
    out = {}
    # ---- A2: a known-true circuit accepted at its exact addition count
    sys.path.insert(0, FROZEN_SCRIPTS)
    from gatec_decomps import load_paper55
    import tensor_data as td
    pU, pV, pW = load_paper55()

    def expand_printed(slp, prefix):
        """Turn the frozen printed SLP into this campaign's gate records."""
        env = {}
        # the frozen SLPs name their free inputs A0..A8 / B0..B8 (0-indexed)
        for i in range(N):
            env[f"{prefix}{i}"] = tuple(1 if j == i else 0 for j in range(N))
        gates = []
        names = {f"{prefix}{i}": i for i in range(N)}
        for gname, terms in slp:
            (s1, t1), (s2, t2) = terms
            v = tuple(s1 * env[t1][k] + s2 * env[t2][k] for k in range(N))
            nm = 100 + len(gates)
            gates.append({"out": nm, "a": [names[t1], s1], "b": [names[t2], s2],
                          "val": list(v)})
            env[gname] = v
            names[gname] = nm
        return gates
    gL = expand_printed(td.LEFT_SLP, "A")
    gR = expand_printed(td.RIGHT_SLP, "B")
    okL, misL = circuit_verdict(gL, pU)
    okR, misR = circuit_verdict(gR, pV)
    if not (okL and okR and len(gL) == 13 and len(gR) == 14):
        die(f"A2: known-true paper55 circuits not accepted at their exact "
            f"counts: L={len(gL)}/{okL} {misL[:3]} R={len(gR)}/{okR} {misR[:3]}")
    tr = transpose_circuit(gL, pU)
    out["A2_paper55"] = {"left": len(gL), "right": len(gR),
                         "left_verified": okL, "right_verified": okR,
                         "total_recount": len(gL) + len(gR) + 28,
                         "transposed_left_adds": tr[0] if tr else None}
    print(f"A2 [PASS] paper55 printed circuits verified at exact counts "
          f"13/14, recount {len(gL)}+{len(gR)}+28="
          f"{len(gL) + len(gR) + 28}")
    # ---- A3: frozen-row re-verification
    frozen_expect = {"U": (12, False), "V": (13, False), "WFac": (13, False)}
    got = {}
    for nm, M in (("U", pU), ("V", pV), ("WFac", pW)):
        c, rp = prep(M)
        ok, st, _ = subset_dfs(c, rp)
        got[nm] = (len(c), ok)
        if (len(c), ok) != frozen_expect[nm]:
            die(f"A3: frozen paper55 row {nm} not reproduced: {(len(c), ok)} "
                f"!= {frozen_expect[nm]}")
    out["A3_frozen_rows"] = {k: list(v) for k, v in got.items()}
    print(f"A3 [PASS] frozen paper55 rows reproduced: {got}")
    # ---- A4: relabeling invariance (seed 17, 8 permutations per side)
    rng = random.Random(17)
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        c0, r0 = prep(M)
        base = (len(c0), subset_dfs(c0, r0)[0])
        for _ in range(8):
            perm = list(range(R))
            rng.shuffle(perm)
            M2 = [M[perm[r]] for r in range(R)]
            c1, r1 = prep(M2)
            if (len(c1), subset_dfs(c1, r1)[0]) != base:
                die(f"A4: relabeling changed {nm}'s d/floor verdict")
    out["A4_relabeling_invariant"] = True
    print("A4 [PASS] 8 relabelings per side leave d and floor verdicts fixed")
    # ---- R5 (amendment 1): a creatable-chain instance, known SAT BY
    # CONSTRUCTION, must come back feasible from BOTH instruments. (The
    # originally pre-registered "slack slots" instance was mis-specified: the
    # encoding's operand pool is inputs u tau, so extra slots add no power.)
    cU, rU = prep(U)
    cur = [0] * N
    cur[0] = 1
    chain = []
    for i in range(1, 8):
        nxt = list(cur)
        nxt[i] = 1
        chain.append(canon(tuple(nxt)))
        cur = nxt
    tau_chain = sorted(set(chain))
    chain_dfs = subset_dfs(*prep([tuple(x) for x in tau_chain]))[0]
    chain_sat, meta = run_cnf(tau_chain, len(tau_chain), tag="R5_chain")
    if not (chain_sat and chain_dfs):
        die(f"R5: the creatable chain (known SAT by construction) came back "
            f"infeasible: cnf={chain_sat} dfs={chain_dfs}")
    out["R5_known_sat_chain"] = {"classes": len(tau_chain), "cnf_sat": chain_sat,
                                 "dfs_achievable": chain_dfs,
                                 "vars": meta["vars"], "clauses": meta["clauses"]}
    print(f"R5 [PASS] creatable chain ({len(tau_chain)} classes at "
          f"T={len(tau_chain)}) feasible in BOTH instruments "
          f"(cnf={chain_sat}, dfs={chain_dfs})")
    # ---- R3: floor-level infeasibility re-certified with a fresh certificate
    satf, metaf = run_cnf(sorted(cU), len(cU), tag="R3_floorU_T14",
                          certificates=True, keep=True)
    if satf:
        die("R3: the floor level came back SAT from the CNF instrument while "
            "the DFS says impossible - instrument disagreement")
    out["R3_floor_certificate"] = metaf["cert"]
    print(f"R3 [PASS] floor T=14 UNSAT re-certified, both checkers rc "
          f"{metaf['cert']['drattrim_rc']}/{metaf['cert']['lratcheck_rc']}")
    # ---- cross-solver control (CaDiCaL must agree on a representative)
    sat_c, _ = run_cnf(sorted(cU), len(cU), tag="xsolver", solver=CADICAL)
    if sat_c != satf:
        die(f"cross-solver disagreement: kissat={satf} cadical={sat_c}")
    out["cross_solver_agrees"] = True
    print("CTRL [PASS] CaDiCaL agrees with kissat on the representative")
    # ---- R4: planted wrong d-count assertion must fire
    fired = False
    try:
        assert len(cU) == 99, "planted wrong d"
    except AssertionError:
        fired = True
    if not fired:
        die("R4: planted assertion did not fire")
    out["R4_planter_fired"] = True
    print("R4 [PASS] planted assertion fired")
    RESULT["phases"]["B-controls"] = out
    return out


@phase("C-universes")
def phase_C(U, V, W):
    out = {}
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        c, _ = prep(M)
        a1 = au_universe(c, second_loop=False)
        a2 = au_universe(c, second_loop=True)
        if a1 != a2:
            die(f"A6: the two enumeration loops disagree on {nm}: "
                f"{len(a1)} vs {len(a2)}")
        sha = au_sha(a1)
        want_sha, want_n = AU_PINS[nm]
        if not sha.startswith(want_sha) or len(a1) != want_n:
            die(f"A6: universe pin mismatch for {nm}: |AU|={len(a1)} "
                f"sha={sha[:16]} (pinned {want_n}/{want_sha})")
        out[nm] = {"size": len(a1), "sha256": sha, "two_loops_agree": True}
        print(f"A6 [PASS] {nm}: |AU|={len(a1)} sha={sha[:16]} (pinned)")
    RESULT["phases"]["C-universes"] = out
    return out


@phase("D-censuses")
def phase_D(U, V, W):
    out = {}
    ck = open(CKPT, "a")
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        cls, _ = prep(M)
        au = au_universe(cls)
        hits, ncert, t0 = [], 0, time.time()
        for idx, a in enumerate(au):
            ext = sorted(set(cls) | {a})
            c2, r2 = prep(ext)
            dfs_ok = subset_dfs(c2, r2)[0]
            keep = (idx < CERT_KEEP_PER_SIDE)
            sat, meta = run_cnf(ext, T_LEVEL, tag=f"{nm}_aux{idx}",
                                certificates=not dfs_ok, keep=keep)
            if sat != dfs_ok:
                die(f"INSTRUMENT DISAGREEMENT {nm} aux{idx}: dfs={dfs_ok} "
                    f"cnf={sat}")
            if not dfs_ok:
                ncert += 1
            if dfs_ok:
                hits.append({"idx": idx, "aux": list(a)})
            ck.write(json.dumps({"side": nm, "idx": idx, "aux": list(a),
                                 "dfs": dfs_ok, "cnf": sat, "agree": True,
                                 "cert": meta["cert"], "s": meta["s"]}) + "\n")
            ck.flush()
            if (idx + 1) % 100 == 0:
                print(f"  {nm} progress {idx + 1}/{len(au)} hits={len(hits)} "
                      f"({time.time() - t0:.0f}s)")
        out[nm] = {"instances": len(au), "hits": len(hits),
                   "hit_list": hits[:5], "certificates_dual_checked": ncert,
                   "complete": True, "seconds": round(time.time() - t0, 1)}
        lb = T_LEVEL if hits else T_LEVEL + 1
        out[nm]["LB"] = lb
        print(f"{nm} COMPLETE: {len(au)}/{len(au)} agreed, feasible={len(hits)}"
              f" -> C({nm}) >= {lb}  ({ncert} certificates dual-checked)")
    ck.close()
    RESULT["phases"]["D-censuses"] = out
    return out


@phase("E-synthesis")
def phase_E(U, V, W, census):
    out = {}
    circuits = {}
    for nm, M in (("U", U), ("V", V), ("WFac", W)):
        best, bestg, verified, discarded = None, None, 0, 0
        t0 = time.time()
        for restart in range(RESTARTS):
            g = greedy_synth(M, nm, restart)
            if g is None:
                continue
            ok, miss = circuit_verdict(g, M)
            if not ok:
                discarded += 1
                continue
            verified += 1
            if best is None or len(g) < best:
                best, bestg = len(g), g
        circuits[nm] = bestg
        out[nm] = {"best_verified_gates": best, "restarts": RESTARTS,
                   "verified": verified, "discarded": discarded,
                   "seconds": round(time.time() - t0, 1),
                   "LB": census[nm]["LB"],
                   "exact": best == census[nm]["LB"]}
        print(f"E {nm}: best VERIFIED circuit = {best} gates "
              f"(LB {census[nm]['LB']}) -> "
              f"{'EXACT' if best == census[nm]['LB'] else 'interval'}"
              f"  [{verified} verified / {discarded} discarded]")
        json.dump(bestg, open(os.path.join(HERE, f"witness_{nm}.json"), "w"),
                  indent=1)
    # R1 / R2 plants on a verified circuit (before any claim)
    gU = circuits["U"]
    pruned = gU[:-1]
    ok_r1, miss_r1 = circuit_verdict(pruned, U)
    if ok_r1:
        die("R1: one-addition-deleted plant was ACCEPTED")
    flipped = [dict(g) for g in gU]
    flipped[0] = dict(flipped[0])
    flipped[0]["a"] = [flipped[0]["a"][0], -flipped[0]["a"][1]]
    ok_r2, _ = circuit_verdict(flipped, U)
    if ok_r2:
        die("R2: operand-perturbed plant was ACCEPTED")
    out["R1_gate_deleted"] = {"rejected": True, "missing": len(miss_r1)}
    out["R2_operand_perturbed"] = {"rejected": True}
    print(f"R1 [PASS] gate-deleted plant rejected ({len(miss_r1)} classes "
          f"uncovered);  R2 [PASS] operand-perturbed plant rejected")
    RESULT["phases"]["E-synthesis"] = out
    return circuits, out


@phase("F-assembly")
def phase_F(U, V, W, circuits):
    gW = circuits["WFac"]
    tr = transpose_circuit(gW, W)
    if tr is None:
        die("F: output-stage construction failed")
    oadds, outs, _ = tr
    ok_map = all(outs[k][r] == W[r][k] for k in range(N) for r in range(R))
    if not ok_map:
        die("F: transposed output stage does not reproduce the WFac map")
    nL, nR = len(circuits["U"]), len(circuits["V"])
    total = nL + nR + oadds
    print(f"F assembled: left={nL} right={nR} output={oadds} total={total} "
          f"(output map exact={ok_map})")
    # END-TO-END over all 81 monomials
    lvals, rvals = eval_circuit(circuits["U"]), eval_circuit(circuits["V"])

    def pick(vals, tv):
        c = canon(tv)
        for nm, v in vals.items():
            if canon(v) == c:
                return v, (1 if list(v) == list(tv) else -1)
        return None, None
    lsel = [pick(lvals, U[r]) for r in range(R)]
    rsel = [pick(rvals, V[r]) for r in range(R)]
    couts = []
    for k in range(N):
        acc = [[0] * N for _ in range(N)]
        for r in range(R):
            w = outs[k][r]
            if not w:
                continue
            lv, ls = lsel[r]
            rv, rs = rsel[r]
            for i in range(N):
                if not lv[i]:
                    continue
                for j in range(N):
                    if rv[j]:
                        acc[i][j] += w * (ls * lv[i]) * (rs * rv[j])
        couts.append(acc)
    tgt = []
    for k in range(N):
        i0, j0 = divmod(k, 3)
        acc = [[0] * N for _ in range(N)]
        for kk in range(3):
            acc[3 * i0 + kk][3 * kk + j0] += 1
        tgt.append(acc)
    bad = [k for k in range(N) if couts[k] != tgt[k]]
    print(f"F END-TO-END: outputs equal to the 3x3 matrix product: "
          f"{N - len(bad)}/9 (mismatches {bad})")
    if bad:
        die(f"F: end-to-end expansion mismatch at outputs {bad}")
    blob = {"left": nL, "right": nR, "output": oadds, "total": total,
            "output_map_exact": True, "end_to_end_outputs_exact": 9}
    json.dump({"left_gates": circuits["U"], "right_gates": circuits["V"],
               "wfac_gates": circuits["WFac"], "output_additions": oadds,
               "total": total,
               "note": "output stage = reverse-mode transposition of the "
                       "verified WFac circuit; verified against the WFac map "
                       "entry-by-entry and by 81-monomial end-to-end expansion"},
              open(os.path.join(HERE, "assembled_laderman.json"), "w"), indent=1)
    RESULT["phases"]["F-assembly"] = blob
    return blob


@phase("H-verdict")
def phase_H(census, synth, asm, printed):
    lbs = {k: census[k]["LB"] for k in ("U", "V", "WFac")}
    ubs = {k: synth[k]["best_verified_gates"] for k in ("U", "V", "WFac")}
    lb_total = lbs["U"] + lbs["V"] + lbs["WFac"] + 14
    ub_total = asm["total"]
    exact = all(lbs[k] == ubs[k] for k in lbs) and lb_total == ub_total
    verdict = "CERTIFIED-EXACT" if exact else "NEGATIVE-BOUNDED"
    blob = {
        "gate": "laderman23-ladder-total", "run_id": RESULT["run_id"],
        "prereg_sha256": RESULT["prereg_sha256"],
        "verdict": verdict,
        "per_side": {k: {"LB": lbs[k], "UB": ubs[k], "exact": lbs[k] == ubs[k]}
                     for k in lbs},
        "LB_total": lb_total, "UB_total": ub_total,
        "published_basic_form": printed,
        "improvement_over_printed": printed["total"] - ub_total,
        "searched_space": {"instances": sum(census[k]["instances"] for k in lbs),
                           "per_side": {k: census[k]["instances"] for k in lbs},
                           "level": f"aux-1 at T={T_LEVEL}",
                           "universes_hash_pinned": True},
        "routes_not_run": {"aux-2 at T=16": "~2.0 CPU-h per side (measured), "
                                            "explicitly outside the "
                                            "pre-registered budget; not "
                                            "attempted, nothing claimed from it"},
        "claim": (
            f"For Laderman's own 1976 rank-23 decomposition at its printed "
            f"(sigma^0, all-monomial) orientation, in the frozen three-stage "
            f"linear SLP model, the minimum addition count is EXACTLY "
            f"{ub_total} = {ubs['U']} + {ubs['V']} + {oadds_str(asm)}: each "
            f"side's cost is exactly {ubs['U']}/{ubs['V']}/{ubs['WFac']} "
            f"(lower bound from the complete floor DFS plus the complete "
            f"hash-pinned aux-1 census at T={T_LEVEL}, dual-instrument agreeing "
            f"on every instance; upper bound from explicit circuits verified by "
            f"exact expansion), and the assembled scheme was verified end to "
            f"end over all 81 monomials against the 3x3 matrix product. "
            f"Laderman's printed basic form uses {printed['total']} additions "
            f"({printed['left']}+{printed['right']}+{printed['output']}), which "
            f"the paper itself presents as unoptimized; the certified minimum is "
            f"{printed['total'] - ub_total} additions below it."
            if exact else
            f"Bounded: total in [{lb_total}, {ub_total}]; per-side "
            f"LB/UB {lbs} / {ubs}."),
        "explicitly_not_claimed": [
            "any other orientation of Laderman's tensor",
            "any other decomposition, alphabet, action, or counting convention",
            "any statement that Laderman's printed count was erroneous (the "
            "paper presents it as deliberately unoptimized)",
            "the global sub-55 question"],
    }
    with open(os.path.join(HERE, "verdict.json"), "w") as f:
        json.dump(blob, f, indent=1)
    print("H verdict:", verdict)
    print("  LB_total", lb_total, "UB_total", ub_total,
          "| printed", printed["total"])
    return blob


def oadds_str(asm):
    return f"{asm['output']} (= {asm['output'] - 14} + 14 transposition)"


if __name__ == "__main__":
    U, V, W, printed = phase_A()
    phase_B(U, V, W)
    phase_C(U, V, W)
    census = phase_D(U, V, W)
    circuits, synth = phase_E(U, V, W, census)
    asm = phase_F(U, V, W, circuits)
    phase_H(census, synth, asm, printed)
    RESULT["cpu_seconds_total"] = round(
        sum(p.get("seconds", 0) for p in RESULT["phases"].values()), 1)
    with open(os.path.join(HERE, "runner_output.json"), "w") as f:
        json.dump(RESULT, f, indent=1, default=str)
    print("\nRUNNER COMPLETE")
