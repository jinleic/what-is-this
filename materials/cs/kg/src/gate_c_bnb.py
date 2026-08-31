"""gate_c_bnb.py — Gate C(a): certified branch-and-bound over the pre-registered box.

Pre-registered (pre_statement.md Addendum 1 + Addendum 2): box (s3,s5,s7) in [0,0.6]^3,
depth 7 per axis, <= 200,000 boxes, 8-way bisection, M-truncated head_lo with M = 23.

CERTIFIED PRUNING BOUND (needs no tail input, since tail >= 0 and every |b_m| >= 0):
    for all s in B:   gamma*(s) = b1(V(s)) - head(s) - tail(s)
                              <= b1(V_min(B)) - head_lo^{(M)}(B) =: UB(B)
  b1(V) = (pi/2)(A_{1,0}^2/V - A_{0,1}^2) is strictly decreasing in V, so b1(V_min(B))
  bounds b1 over B; head_lo^{(M)}(B) = sum_{3<=m<=M odd} inf_B |b_m| <= head(s) because
  the discarded terms are non-negative. A box with UB(B) <= gamma_ref provably contains
  no scheme beating gamma_ref.

References (fixed in the pre-statement, no post-hoc change):
  gamma_paper  = 0.881545409            paper's own gamma (the "beat 3.47e-4" threshold)
  gamma_gateA  = 0.881557917504162      repo Gate A gamma* (paper's cited B3 tail)
  gamma_qhead  = 0.881562493211119      polished quintic head-objective optimum (s7 = 0)
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, arb_poly, ctx
from core import pi
from gate_c_fast import (load_grid_entries, build_Q_polys, rho_poly_box, PREC)

GAMMA = {
    "paper": arb("0.881545409"),
    "gateA": arb("0.881557917504162"),
    "qhead": arb("0.881562493211119"),
}
BOX0 = ((0.0, 0.6), (0.0, 0.6), (0.0, 0.6))
MAX_DEPTH = 7
BUDGET = 200000
M_TRUNC = 23


def box_bound(box, Qs, a10, a01, mtrunc=M_TRUNC, prec=PREC):
    """(UB, diagnostics) — certified upper bound on gamma* over the box."""
    (s3lo, s3hi), (s5lo, s5hi), (s7lo, s7hi) = box
    rho_p, Vlo, Vhi = rho_poly_box(s3lo, s3hi, s5lo, s5hi, s7lo, s7hi, prec)
    ctx.prec = prec
    pi2 = pi(prec) / arb(2)
    # b1 upper over box = b1 at V_min = Vlo
    b1_hi = pi2 * (a10 / arb(Vlo) - a01)
    # truncated head lower bound
    T = [arb(0)] * (mtrunc + 1)
    Pa = arb_poly([arb(1)])
    for a in range(0, mtrunc + 1):
        Q = Qs.get(a)
        if Q is not None:
            pr = Q * Pa
            for k in range(min(len(pr), mtrunc + 1)):
                T[k] += pr[k]
        if a < mtrunc:
            Pa = Pa * rho_p
            if len(Pa) > mtrunc + 1:
                Pa = arb_poly([Pa[k] for k in range(mtrunc + 1)])
    head_lo = arb(0)
    for m in range(3, mtrunc + 1, 2):
        lo = abs(pi2 * T[m]).lower()
        if lo > 0:
            head_lo += lo
    ub = b1_hi.upper() - head_lo
    return ub, dict(b1_hi=b1_hi, head_lo=head_lo, Vlo=Vlo, Vhi=Vhi)


def bisect(box):
    (a0, a1), (b0, b1), (c0, c1) = box
    am, bm, cm = (a0 + a1) / 2, (b0 + b1) / 2, (c0 + c1) / 2
    out = []
    for xa in ((a0, am), (am, a1)):
        for xb in ((b0, bm), (bm, b1)):
            for xc in ((c0, cm), (cm, c1)):
                out.append((xa, xb, xc))
    return out


def run(ref_key="gateA", max_depth=MAX_DEPTH, budget=BUDGET, logpath=None, verbose=True):
    gamma_ref = GAMMA[ref_key]
    ents = load_grid_entries()
    Qs, a10, a01 = build_Q_polys(ents)
    t0 = time.time()
    stack = [(BOX0, 0)]
    n_eval = n_pruned = 0
    pruned_vol = 0.0
    survivors = []
    logf = open(logpath, "w") if logpath else None
    while stack:
        if n_eval >= budget:
            survivors.extend([(b, d, None) for (b, d) in stack])
            if verbose:
                print(f"BUDGET EXHAUSTED at {n_eval}; {len(stack)} boxes unexplored")
            break
        box, depth = stack.pop()
        ub, info = box_bound(box, Qs, a10, a01)
        n_eval += 1
        vol = 1.0
        for (lo, hi) in box:
            vol *= (hi - lo)
        if ub <= gamma_ref:
            n_pruned += 1
            pruned_vol += vol
            if logf:
                logf.write(json.dumps(dict(box=[[l, h] for l, h in box], depth=depth,
                                           ub=ub.str(18), verdict="PRUNED")) + "\n")
            continue
        if depth >= max_depth:
            survivors.append((box, depth, ub))
            if logf:
                logf.write(json.dumps(dict(box=[[l, h] for l, h in box], depth=depth,
                                           ub=ub.str(18), verdict="SURVIVOR")) + "\n")
            continue
        stack.extend((c, depth + 1) for c in bisect(box))
        if verbose and n_eval % 5000 == 0:
            print(f"  [{n_eval}] depth {depth} stack {len(stack)} pruned {n_pruned} "
                  f"vol {pruned_vol/0.216*100:.4f}%  surv {len(survivors)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    if logf:
        logf.close()
    return dict(ref=ref_key, gamma_ref=gamma_ref.str(18), n_eval=n_eval, n_pruned=n_pruned,
                pruned_frac=pruned_vol / 0.216, survivors=survivors,
                seconds=time.time() - t0)


def summarize(out):
    surv = out["survivors"]
    print(f"\nref {out['ref']} = {out['gamma_ref']}")
    print(f"evals {out['n_eval']}  pruned {out['n_pruned']}  "
          f"certified-pruned volume {out['pruned_frac']*100:.5f}%  [{out['seconds']:.0f}s]")
    print(f"survivors: {len(surv)}")
    if surv:
        s7max = max(b[2][1] for (b, d, u) in surv)
        s7min = min(b[2][0] for (b, d, u) in surv)
        s3rng = (min(b[0][0] for (b, d, u) in surv), max(b[0][1] for (b, d, u) in surv))
        s5rng = (min(b[1][0] for (b, d, u) in surv), max(b[1][1] for (b, d, u) in surv))
        print(f"  survivor hull: s3 in [{s3rng[0]:.6f},{s3rng[1]:.6f}]  "
              f"s5 in [{s5rng[0]:.6f},{s5rng[1]:.6f}]  s7 in [{s7min:.6f},{s7max:.6f}]")
        vol = sum((b[0][1]-b[0][0])*(b[1][1]-b[1][0])*(b[2][1]-b[2][0]) for (b,d,u) in surv)
        print(f"  survivor volume {vol:.3e} ({vol/0.216*100:.5f}% of box)")
    return surv


if __name__ == "__main__":
    key = sys.argv[1] if len(sys.argv) > 1 else "gateA"
    log = f"/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_bnb_{key}.jsonl"
    out = run(ref_key=key, logpath=log)
    surv = summarize(out)
    with open(f"/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_bnb_{key}_summary.json", "w") as f:
        json.dump(dict(ref=out["ref"], gamma_ref=out["gamma_ref"], n_eval=out["n_eval"],
                       n_pruned=out["n_pruned"], pruned_frac=out["pruned_frac"],
                       seconds=out["seconds"], n_survivors=len(surv),
                       survivors=[dict(box=[[l, h] for l, h in b], depth=d,
                                       ub=(u.str(20) if u is not None else None))
                                  for (b, d, u) in surv]), f, indent=1)
    print(f"wrote gc_bnb_{key}_summary.json")
