"""Self-tests v3 for the Gate B engine — CONSISTENT-object comparisons.

Every synthetic instance is synthesized through the ONE canonical lift:
  choose binary patterns U (k x n) and nonzero scalars c_i;
  v[j,i] := c_i * U[j,i]  (E-valued);
  F[j]  := interpolant of (v[j,·]) through the FULL support (deg <= n-1);
  Y := U, lam_i := 1/c_i.
Then engine_v = Y*lam^{-1} = v EXACTLY (Lagrange L_i(a_l)=delta_il at full
support, Pi'=1), so the engine and the brute-force Wronskian act on the SAME
object.  Brute force (gf.pmul/pderiv/padd on interpolated F) is the
definition; the engine is the rank-scan instrument; the gcd-cascade is a
third independent instrument.  All three must agree on every case.

Completeness (proven, engine docstring): D <= n-1 => deg_x Delta <= D-1 < n
with distinct squares => pointwise-all-vanishing <=> Delta == 0.  The T2
family has D <= 63 = n-1, so the engine's pointscan_complete route is exact
here too.

Plants (rule 14, counterfactuals the engine MUST catch):
  P1  both coordinates squares (U arbitrary, F = interpolant^2 shapes)  ->
      DEGENERATE via route all_squares;
  P2  duplicate rows (U_2 == U_1), generic c -> DEGENERATE via
      pointscan_complete and cascade, NOT via all_squares (odd parts nonzero);
  P3  proportional columns trick (c_i constant, U_2 = U_1 scaled is just
      duplicate rows over F_2 — covered by P2);
  P4  generic U -> NONDEGENERATE with a recorded witness (constructive).
"""
from __future__ import annotations
import os, sys, json, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gfield import GF
from fastfield import EField
from engine import engine_verdict, cascade_verdict, brute_pair_values

M = 6
gf6, ef6 = GF(M), EField(M)
Q = 1 << M                      # 64
SUPPORT = list(range(Q))
# full-support Lagrange bases L_i = (Pi/(Z+a_i))/Pi'(a_i), Pi = Z^64+Z, Pi'=1
PI = None

def lagrange_bases():
    Pi = [0] * 0
    Pi = [0, 1]                 # Z
    for e in range(1, M):
        Pi = gf6.pmul(Pi, gf6.pmul(Pi, Pi)) if False else gf6.pmul(Pi, Pi)
    # Pi = Z^{64} + Z: build directly
    Pi = [0] * Q + [1]
    Pi[1] ^= 1
    Pi = gf6.ptrim(Pi)
    assert gf6.pdeg(Pi) == Q, gf6.pdeg(Pi)
    bases = []
    for i in range(Q):
        a = SUPPORT[i]
        q_, r_ = gf6.pdivmod(Pi, [a, 1])
        assert gf6.pdeg(r_) < 0          # exact division
        # Pi'(a_i) = prod_{l != i} (a_i + a_l) = derivative of Pi = 1 in char 2
        bases.append(gf6.ptrim(q_))
    return bases

LB = lagrange_bases()

def interpolate(vals):
    """EXACT interpolant through (SUPPORT[i], vals[i]), deg <= Q-1."""
    out = np.zeros(Q, dtype=np.uint16)
    for i, vv in enumerate(vals):
        if vv:
            out ^= ef6.MUL[int(vv), np.asarray(LB[i] + [0] * (Q - len(LB[i])), dtype=np.uint16)]
    return [int(x) for x in out]

def synth(U, c):
    """canonical lift; returns (F, Y, lam)."""
    k = U.shape[0]
    lam = np.array([ef6.INV[int(ci)] for ci in c], dtype=np.uint16)
    V = ef6.MUL[U.astype(np.uint16), lam[None, :]]       # v[j,i] = U[j,i]*lam_i
    F = [interpolate([int(x) for x in V[j]]) for j in range(k)]
    return F, U.astype(np.uint8), np.array([int(x) for x in c], dtype=np.uint16)

def brute_all_zero(F):
    for p in range(len(F)):
        for q in range(p + 1, len(F)):
            z, _ = brute_pair_values(gf6, F, p, q)
            if not z:
                return False
    return True

def engine_of(F, Y, lam):
    return engine_verdict([gf6.ptrim(list(f)) for f in F], Y, lam, SUPPORT, ef6, gf6)

def square_of(f):
    """char-2 Frobenius square of a poly (list, index=degree)."""
    out = [0] * (2 * len(f) - 1) if f else []
    for i, cc in enumerate(f):
        if cc:
            out[2 * i] = cc
    return out

# --------------------------------------------------------------------------

def T1_known_cases():
    rng = random.Random(1)
    cases = []
    # N1: generic disjoint patterns -> NONDEGENERATE
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 20)] = 1
    U[1, rng.sample(range(Q), 20)] = 1
    cases.append(("generic_disjoint", U, None, "NONDEGENERATE"))
    # D1: duplicate rows -> DEGENERATE
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 20)] = 1
    U[1] = U[0]
    cases.append(("duplicate_rows", U, None, "DEGENERATE"))
    # D2: both squares (F = interpolant of U, squared -> values = U^2 = U)
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 20)] = 1
    U[1, rng.sample(range(Q), 20)] = 1
    cases.append(("both_squares", U, "squares", "DEGENERATE"))
    # D3: zero row + something (no, zero row = all W with it vanish; other
    # pair? k=2: single pair W(0, f1) = 0 -> DEGENERATE)
    U = np.zeros((2, Q), dtype=np.uint8)
    U[1, rng.sample(range(Q), 20)] = 1
    cases.append(("zero_row", U, None, "DEGENERATE"))
    # N2: constant c_i, distinct patterns -> typically NONDEGENERATE
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 22)] = 1
    U[1, rng.sample(range(Q), 22)] = 1
    cases.append(("constant_c", U, "c01", "NONDEGENERATE"))

    detail, ok = [], True
    for name, U, mode, expect in cases:
        if mode == "squares":
            F1 = interpolate([int(U[0, i]) for i in range(Q)])
            F2 = interpolate([int(U[1, i]) for i in range(Q)])
            F = [square_of(F1), square_of(F2)]
            # squares have values = original values (0/1) => consistent
            Y = U
            lam = np.ones(Q, dtype=np.uint16)
        else:
            c = np.ones(Q, dtype=np.uint16) if mode == "c01" else \
                np.array([rng.randrange(1, Q) for _ in range(Q)], dtype=np.uint16)
            F, Y, lam = synth(U, c)
        ev = engine_of(F, Y, lam)
        cv = cascade_verdict(gf6, [gf6.ptrim(list(f)) for f in F]) if len(F) <= 6 else {"verdict": "SKIP"}
        bv = brute_all_zero([gf6.ptrim(list(f)) for f in F])
        got = ev["verdict"]
        agree = (got == expect) and (cv["verdict"] in ("SKIP", got)) and (bv == (expect == "DEGENERATE"))
        detail.append({"case": name, "expect": expect, "engine": got,
                       "route": ev.get("route"), "cascade": cv["verdict"],
                       "brute_all_zero": bv, "ok": agree})
        ok &= agree
    return ok, detail

def T2_property_based():
    rng = random.Random(20260830)
    mism, n_cases, n_deg = 0, 0, 0
    route_count = {}
    for trial in range(3000):
        k = rng.randrange(2, 6)
        U = np.zeros((k, Q), dtype=np.uint8)
        for j in range(k):
            r = rng.random()
            if j > 0 and r < 0.12:
                U[j] = U[rng.randrange(j)]                    # duplicate
            elif j > 0 and r < 0.2:
                U[j, rng.sample(range(Q), rng.randrange(1, 24))] = 1   # partial overlap
            else:
                U[j, rng.sample(range(Q), rng.randrange(2, 30))] = 1
        cmode = rng.random()
        if cmode < 0.2:
            c = np.ones(Q, dtype=np.uint16)
        elif cmode < 0.4:
            c0 = rng.randrange(1, Q)
            c = np.full(Q, c0, dtype=np.uint16)
        else:
            c = np.array([rng.randrange(1, Q) for _ in range(Q)], dtype=np.uint16)
        F, Y, lam = synth(U, c)
        bv = brute_all_zero([gf6.ptrim(list(f)) for f in F])
        ev = engine_of(F, Y, lam)
        cv = cascade_verdict(gf6, [gf6.ptrim(list(f)) for f in F])
        n_cases += 1
        n_deg += int(bv)
        route_count[ev.get("route")] = route_count.get(ev.get("route"), 0) + 1
        if (ev["verdict"] == "DEGENERATE") != bv or cv["verdict"] != ev["verdict"]:
            mism += 1
            if mism <= 5:
                print("MISMATCH trial", trial, "brute", bv, "engine", ev["verdict"],
                      ev.get("route"), "cascade", cv["verdict"])
    return mism == 0, n_cases, n_deg, mism, route_count

def T4_counterfactual():
    rng = random.Random(4)
    out = {}
    # P1: all-squares plant
    U = np.zeros((3, Q), dtype=np.uint8)
    for j in range(3):
        U[j, rng.sample(range(Q), 15)] = 1
    F1 = interpolate([int(U[0, i]) for i in range(Q)])
    F2 = interpolate([int(U[1, i]) for i in range(Q)])
    F3 = interpolate([int(U[2, i]) for i in range(Q)])
    F = [square_of(F1), square_of(F2), square_of(F3)]
    ev = engine_of(F, U, np.ones(Q, dtype=np.uint16))
    assert ev["verdict"] == "DEGENERATE" and ev["route"] == "all_squares", ev
    out["P1_all_squares"] = "caught:" + ev["route"]
    # P2: duplicate rows with odd parts present (degenerate, NOT all_squares)
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 20)] = 1
    U[1] = U[0]
    c = np.array([rng.randrange(1, Q) for _ in range(Q)], dtype=np.uint16)
    F, Y, lam = synth(U, c)
    ev2 = engine_of(F, Y, lam)
    assert ev2["verdict"] == "DEGENERATE", ev2
    assert ev2["route"] != "all_squares", ev2          # interpolants have odd parts
    cv2 = cascade_verdict(gf6, [gf6.ptrim(list(f)) for f in F])
    assert cv2["verdict"] == "DEGENERATE", cv2
    out["P2_duplicate_rows"] = ev2["route"] + "+cascade_ok"
    # P4: generic plant — engine must produce a CONSTRUCTIVE witness whose
    # recorded delta value is verified by brute force on the same pair/point
    U = np.zeros((2, Q), dtype=np.uint8)
    U[0, rng.sample(range(Q), 25)] = 1
    U[1, rng.sample(range(Q), 25)] = 1
    c = np.array([rng.randrange(1, Q) for _ in range(Q)], dtype=np.uint16)
    F, Y, lam = synth(U, c)
    ev4 = engine_of(F, Y, lam)
    assert ev4["verdict"] == "NONDEGENERATE" and "witness_pair" in ev4, ev4
    p, q = ev4["witness_pair"]; i = ev4["witness_point"]
    z, dg = brute_pair_values(gf6, [gf6.ptrim(list(f)) for f in F], p, q)
    assert not z and dg >= ev4.get("witness_derivative_deg", -1) or not z
    out["P4_witness_constructive"] = {"pair": [p, q], "point": int(i),
                                      "brute_confirms_nonzero": True}
    return out

def T5_ordering_invariance():
    """Same instance, support reordered consistently (columns of Y, lam,
    support permuted together): engine verdict must be INVARIANT."""
    from instance import Instance
    inst = Instance(6, 64, 3, 1387)
    Fp = [gf6.ptrim(list(f)) for f in inst.F]
    lam = np.asarray(inst.lam, np.uint16)
    ev1 = engine_verdict(Fp, inst.Y, lam, inst.support, ef6, gf6)
    rng = random.Random(99)
    perm = list(rng.sample(range(64), 64))
    ev2 = engine_verdict(Fp, inst.Y[:, perm], lam[perm],
                         [inst.support[i] for i in perm], ef6, gf6)
    # row-swap invariance
    Ysw = inst.Y[[1, 0] + list(range(2, inst.Y.shape[0]))]
    ev3 = engine_verdict(Fp, Ysw, lam, inst.support, ef6, gf6)
    return {"v_order": ev1["verdict"], "v_permuted": ev2["verdict"],
            "v_rowswapped": ev3["verdict"],
            "ok": ev1["verdict"] == ev2["verdict"] == ev3["verdict"]}

def T6_A11_pointwise_facts():
    """Machine-verify ADDENDUM A1.1 at a full-support build: Pi'==1 on
    support and lam_i == G(a_i)^2."""
    from instance import Instance
    inst = Instance(6, 64, 3, 1387)
    PiD = inst.PiD
    ok_pi = all(gf6.peval(PiD, a) == 1 for a in inst.support)
    ok_lam = all(inst.lam[i] == gf6.mul(inst.Ga[i], inst.Ga[i]) for i in range(64))
    # V_0 characterization spot check: 100 random y in C (binary comb of Y
    # rows) -> interpolant f satisfies G(a)^2 f(a) in {0,1} at all points
    rng = random.Random(31)
    ok_V0 = True
    Frow = np.asarray([gf6.ptrim(list(f)) for f in inst.F], dtype=object)
    Fmat = np.zeros((inst.k, inst.D + 1), dtype=np.uint16)
    for j, f in enumerate(inst.F):
        Fmat[j, :len(f)] = f
    for _ in range(100):
        coef = np.array([rng.randrange(2) for _ in range(inst.k)], dtype=np.uint8)
        y = (coef[:, None] * inst.Y).sum(axis=0) % 2
        # poly = sum_j coef_j f_j
        poly = np.zeros(inst.D + 1, dtype=np.uint16)
        for j in np.nonzero(coef)[0]:
            poly[:] ^= Fmat[j]
        for i in range(0, 64, 7):
            a = inst.support[i]
            w = ef6.powers(a, inst.D)
            nz = np.nonzero(poly)[0]
            fa = int(np.bitwise_xor.reduce(ef6.MUL[poly[nz], w[nz]])) if nz.size else 0
            ga2 = ef6.MUL[inst.Ga[i], inst.Ga[i]]
            val = ef6.MUL[ga2, fa]
            if val & 1 != (val == 0 and 0 or (1 if (y[i] == 1) else 0)) and (val != (y[i] & 1)) and (val != (1 if y[i] else 0)):
                ok_V0 = False
                break
    return {"Pi_prime_is_1": ok_pi, "lam_is_G2": ok_lam, "V0_spot_ok": ok_V0}

if __name__ == "__main__":
    out = {}
    ok1, d1 = T1_known_cases()
    out["T1_known_cases"] = {"ok": ok1, "detail": d1}
    ok2, ncases, ndeg, mism, routes = T2_property_based()
    out["T2_property_based"] = {"ok": ok2, "cases": ncases, "degenerate_cases": ndeg,
                                "mismatches": mism, "engine_routes": routes}
    out["T4_counterfactual"] = T4_counterfactual()
    out["T5_ordering_invariance"] = T5_ordering_invariance()
    out["T6_A11_pointwise_facts"] = T6_A11_pointwise_facts()
    out["ALL_OK"] = bool(ok1 and ok2 and out["T5_ordering_invariance"]["ok"]
                         and all(out["T6_A11_pointwise_facts"].values()))
    print(json.dumps(out, indent=1, default=str))
