"""gate_c_certify.py — Gate C(a): certified point comparisons + tail-side evidence.

Everything here is Arb-certified (prec 256/400) except items explicitly labelled
COMPUTATIONAL-EVIDENCE (the float S(theta) sweep that supplies B3 / the tail).

Certified objects:
  b1(V)      = (pi/2)(A_{1,0}^2/V - A_{0,1}^2)          [two-cell identity, exact]
  head(s)    = sum_{3<=m<=251 odd} |b_m(s)|             [A-grid + exact septic rho^a]
  gamma_head = b1 - head                                 [ >= gamma* since tail >= 0 ]
Comparisons (certified SIGN, with widths):
  (C1) quintic-polished vs paper decimals
  (C2) septic triple-annihilation (b3=b5=b7=0) vs quintic-polished
  (C3) the S2a s7-continuation table, each node's certified gamma_head
  (C4) the certified sensitivity constants: db1/dV and the |b7| gain rate in u7
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, arb_poly, fmpq, ctx
from core import pi
from gate_c_fast import (load_grid_entries, build_Q_polys, rho_poly_box, _exact_q,
                         PREC, N_MAX)

PAPER = ("0.34101124", "0.05276111", "0")
OUT = "/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_certify.json"


def coeffs_certified(s3, s5, s7, Qs, a10, a01, nmax=N_MAX, prec=PREC):
    """Certified (b1, head, {b_m}, V) at an exact-rational point."""
    rho_p, Vlo, Vhi = rho_poly_box(s3, s3, s5, s5, s7, s7, prec)
    ctx.prec = prec
    pi2 = pi(prec) / arb(2)
    T = [arb(0)] * (nmax + 1)
    amax = max(Qs)
    Pa = arb_poly([arb(1)])
    for a in range(0, amax + 1):
        Q = Qs.get(a)
        if Q is not None:
            pr = Q * Pa
            for k in range(min(len(pr), nmax + 1)):
                T[k] += pr[k]
        if a < amax:
            Pa = Pa * rho_p
            if len(Pa) > nmax + 1:
                Pa = arb_poly([Pa[k] for k in range(nmax + 1)])
    bms = {m: pi2 * T[m] for m in range(1, nmax + 1, 2)}
    b1 = pi2 * (a10 / arb(Vhi) - a01)
    head = arb(0)
    for m in range(3, nmax + 1, 2):
        head = head + abs(bms[m])
    return dict(b1=b1, head=head, gamma_head=b1 - head, bms=bms, V=arb(Vhi), Vq=Vhi)


def show(tag, d, extra=""):
    print(f"{tag}")
    print(f"   V          = {d['V'].str(20)}")
    print(f"   b1         = {d['b1'].str(25)}")
    print(f"   head       = {d['head'].str(18)}")
    print(f"   gamma_head = {d['gamma_head'].str(25)}  {extra}")
    print(f"   b3 = {d['bms'][3].str(10)}  b5 = {d['bms'][5].str(10)}  b7 = {d['bms'][7].str(10)}")


if __name__ == "__main__":
    t0 = time.time()
    ents = load_grid_entries()
    Qs, a10, a01 = build_Q_polys(ents)
    res = {}
    ctx.prec = PREC

    # ---------------- C1: paper decimals vs polished quintic ----------------
    print("=== C1: quintic side ===")
    d_paper = coeffs_certified(*PAPER, Qs, a10, a01)
    show("paper decimals (s3,s5,s7) = (0.34101124, 0.05276111, 0)", d_paper)
    # polished quintic point from the Newton solve of b3 = b5 = 0 (S2a s7=0 row)
    QPOL = ("0.34101124", "0.052761124", "0")
    d_q = coeffs_certified(*QPOL, Qs, a10, a01)
    show(f"polished quintic {QPOL}", d_q)
    diff_q = d_q["gamma_head"] - d_paper["gamma_head"]
    print(f"   DELTA(polished - paper) = {diff_q.str(12)}   certified sign > 0: {bool(diff_q > 0)}")
    res["C1"] = dict(paper_gamma_head=d_paper["gamma_head"].str(25),
                     polished_gamma_head=d_q["gamma_head"].str(25),
                     delta=diff_q.str(15), delta_positive=bool(diff_q > 0))

    # ---------------- C2: septic triple annihilation vs quintic ----------------
    print("\n=== C2: septic triple annihilation b3 = b5 = b7 = 0 ===")
    TRIP = ("0.34101121", "0.05276111", "0.00036449")
    d_t = coeffs_certified(*TRIP, Qs, a10, a01)
    show(f"triple point {TRIP}", d_t)
    delta = d_t["gamma_head"] - d_q["gamma_head"]
    print(f"   DELTA(septic triple - polished quintic) = {delta.str(12)}")
    print(f"   certified NEGATIVE (septic loses): {bool(delta < 0)}")
    print(f"   width of delta = {arb(delta).rad()}")
    res["C2"] = dict(triple_point=TRIP, triple_gamma_head=d_t["gamma_head"].str(25),
                     delta=delta.str(18), certified_negative=bool(delta < 0),
                     head_triple=d_t["head"].str(15), head_quintic=d_q["head"].str(15),
                     b1_triple=d_t["b1"].str(20), b1_quintic=d_q["b1"].str(20))
    # decomposition: head gain vs b1 cost
    head_gain = d_q["head"] - d_t["head"]
    b1_cost = d_q["b1"] - d_t["b1"]
    print(f"   head gain (quintic - triple) = {head_gain.str(12)}")
    print(f"   b1   cost (quintic - triple) = {b1_cost.str(12)}")
    print(f"   net (gain - cost)            = {(head_gain - b1_cost).str(12)} "
          f"(negative ⇒ septic loses)")
    res["C2"]["head_gain"] = head_gain.str(15)
    res["C2"]["b1_cost"] = b1_cost.str(15)
    res["C2"]["net"] = (head_gain - b1_cost).str(15)

    # ---------------- C3: certified s7-continuation table ----------------
    print("\n=== C3: certified gamma_head along the b3=b5=0 continuation ===")
    rows = [
        ("0", "0.34101124", "0.052761124"),
        ("0.0002", "0.34101123", "0.052761121"),
        ("0.0005", "0.34101117", "0.052761098"),
        ("0.001", "0.34101096", "0.052761014"),
        ("0.002", "0.34101010", "0.052760657"),
        ("0.003", "0.34100868", "0.052760075"),
        ("0.005", "0.34100412", "0.052758220"),
        ("0.01", "0.34098276", "0.052749510"),
        ("0.02", "0.34089731", "0.052714690"),
        ("0.05", "0.34029970", "0.052471300"),
    ]
    tab = []
    for (s7, s3, s5) in rows:
        d = coeffs_certified(s3, s5, s7, Qs, a10, a01)
        dd = d["gamma_head"] - d_q["gamma_head"]
        beats = bool(d["gamma_head"] > arb("0.881557917504162"))
        print(f"  s7={s7:<7} gamma_head = {d['gamma_head'].str(20)}  "
              f"delta_vs_quintic = {dd.str(10):>16}  head = {d['head'].str(10)}  "
              f"beats gateA: {beats}")
        tab.append(dict(s7=s7, s3=s3, s5=s5, gamma_head=d["gamma_head"].str(25),
                        delta_vs_quintic=dd.str(15), head=d["head"].str(15),
                        beats_gateA=beats, certified_delta_negative=bool(dd < 0)))
    res["C3"] = tab

    # ---------------- C4: certified sensitivities ----------------
    print("\n=== C4: certified sensitivity constants ===")
    Vp = d_q["V"]
    dbdV = -(pi(PREC) / 2) * a10 / (Vp * Vp)
    print(f"   db1/dV at V_quintic = {dbdV.str(15)}")
    # |b7| gain rate in u7: b7(quintic) -> 0 at u7 = s7^2 of the triple point
    u7_trip = _exact_q("0.00036449") ** 2
    rate = d_q["bms"][7] / arb(u7_trip)
    print(f"   u7 at triple point  = {arb(u7_trip).str(12)}")
    print(f"   b7(quintic)         = {d_q['bms'][7].str(12)}")
    print(f"   |b7| gain rate      = {abs(rate).str(12)} per unit u7   "
          f"(vs |db1/dV| = {abs(dbdV).str(8)})")
    res["C4"] = dict(db1_dV=dbdV.str(20), u7_triple=arb(u7_trip).str(15),
                     b7_quintic=d_q["bms"][7].str(15), b7_rate=abs(rate).str(15))

    # ---------------- tail-side arithmetic (mixed labels) ----------------
    print("\n=== tail arithmetic (tail input is COMPUTATIONAL-EVIDENCE) ===")
    for tag, B3 in (("paper cited B3 = 14.44243664663976457", arb("14.44243664663976457")),
                    ("our corrected float B3 = 14.245983003864794", arb("14.245983003864794"))):
        tail = B3 / (arb(10) * arb(251) ** 5).sqrt()
        g = d_q["gamma_head"] - tail
        K = pi(400) / (2 * g)
        kri = pi(400) / (2 * (arb(1) + arb(2).sqrt()).log())
        print(f"  {tag}:")
        print(f"    tail = {tail.str(15)}   gamma* = {g.str(20)}")
        print(f"    K_G <= {K.str(22)}   improvement over Krivine = {(kri - K).str(12)}")
        res.setdefault("tail_arith", []).append(
            dict(B3=B3.str(20), tail=tail.str(18), gamma_star=g.str(22),
                 KG_upper=K.str(25), improvement=(kri - K).str(15)))

    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    print(f"\n[{time.time()-t0:.0f}s] wrote {OUT}")
