"""gate_c_valley.py — Gate C(a) stages S2/S3: valley continuation + free local search.

Pre-registered in pre_statement.md ADDENDUM 2 (committed before this ran).

Objective (certified per point, prec 256):
    gamma_head(s3,s5,s7) = b1(V) - head,   head = sum_{3<=m<=251 odd} |b_m|,
    b1(V) = (pi/2)(A_{1,0}^2/V - A_{0,1}^2),  V = 1+s3^2+s5^2+s7^2.
gamma*(s) = gamma_head(s) - tail(s) with tail(s) > 0, so gamma_head is a certified UPPER
bound for the achievable gamma at that point: if gamma_head(s) <= gamma_ref then s cannot
beat gamma_ref, no tail needed.

S2: on the paper's own construction manifold b_3 = b_5 = 0, continued in s7 (Newton in
    (s3,s5) at fixed s7), plus the septic triple-annihilation b_3 = b_5 = b_7 = 0.
S3: Nelder-Mead on gamma_head from the paper point and from the S2 winner.
"""
import os, sys, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, arb_poly, ctx
from core import pi
from gate_c_fast import (load_grid_entries, build_Q_polys, rho_poly_box, PREC, N_MAX)

GAMMA_PAPER = 0.881545409
GAMMA_GATEA = 0.881557917504162
S3_P, S5_P = 0.34101124, 0.05276111

_ENTS = None
_QS = None
_A10 = None
_A01 = None


def setup():
    global _ENTS, _QS, _A10, _A01
    if _QS is None:
        _ENTS = load_grid_entries()
        _QS, _A10, _A01 = build_Q_polys(_ENTS)
    return _QS, _A10, _A01


def coeffs_at(s3: float, s5: float, s7: float, nmax=N_MAX):
    """Certified (b1, [b_m for odd m<=nmax], head) at a point. Returns arbs."""
    Qs, a10, a01 = setup()
    s3s, s5s, s7s = repr(float(s3)), repr(float(s5)), repr(float(s7))
    rho_p, Vlo, Vhi = rho_poly_box(s3s, s3s, s5s, s5s, s7s, s7s)
    ctx.prec = PREC
    pi2 = pi(PREC) / arb(2)
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
        head += abs(bms[m]).upper()
    return b1, bms, head, Vhi


def gamma_head(s3, s5, s7):
    b1, bms, head, V = coeffs_at(s3, s5, s7)
    g = b1 - head
    return float(g.mid()), dict(b1=float(b1.mid()), head=float(head.mid()),
                                V=float(arb(V).mid()),
                                b3=float(bms[3].mid()), b5=float(bms[5].mid()),
                                b7=float(bms[7].mid()), b9=float(bms[9].mid()),
                                b11=float(bms[11].mid()),
                                g_arb=g)


def bvec(s3, s5, s7, which=(3, 5)):
    _, bms, _, _ = coeffs_at(s3, s5, s7, nmax=31)
    return [float(bms[m].mid()) for m in which]


def newton_solve(s7, which=(3, 5), start=None, iters=40, tol=1e-16, verbose=False):
    """Solve b_m = 0 for m in `which` by Newton in the first len(which) of (s3,s5,s7).
    For which=(3,5): unknowns (s3,s5) at fixed s7. For which=(3,5,7): unknowns all three
    (s7 start given, all three move)."""
    nvar = len(which)
    x = list(start) if start else ([S3_P, S5_P] if nvar == 2 else [S3_P, S5_P, s7])
    for it in range(iters):
        if nvar == 2:
            F = bvec(x[0], x[1], s7, which)
        else:
            F = bvec(x[0], x[1], x[2], which)
        nrm = max(abs(f) for f in F)
        if nrm < tol:
            break
        # numeric Jacobian
        J = [[0.0] * nvar for _ in range(nvar)]
        h = 1e-7
        for j in range(nvar):
            xp = list(x)
            xp[j] += h
            if nvar == 2:
                Fp = bvec(xp[0], xp[1], s7, which)
            else:
                Fp = bvec(xp[0], xp[1], xp[2], which)
            for i in range(nvar):
                J[i][j] = (Fp[i] - F[i]) / h
        # solve J dx = -F (Gauss)
        A = [row[:] + [-F[i]] for i, row in enumerate(J)]
        for c in range(nvar):
            p = max(range(c, nvar), key=lambda r: abs(A[r][c]))
            A[c], A[p] = A[p], A[c]
            if abs(A[c][c]) < 1e-300:
                return None, None
            for r in range(nvar):
                if r != c:
                    f = A[r][c] / A[c][c]
                    for k in range(c, nvar + 1):
                        A[r][k] -= f * A[c][k]
        dx = [A[i][nvar] / A[i][i] for i in range(nvar)]
        damp = 1.0
        for j in range(nvar):
            x[j] += damp * dx[j]
            if x[j] < 0:
                x[j] = 0.0
        if verbose:
            print(f"    newton it{it}: x = {x}, |F| = {nrm:.3e}")
    if nvar == 2:
        return (x[0], x[1], s7), max(abs(f) for f in bvec(x[0], x[1], s7, which))
    return (x[0], x[1], x[2]), max(abs(f) for f in bvec(x[0], x[1], x[2], which))


def nelder_mead(f, x0, step=1e-3, maxit=400, tol=1e-15):
    """Maximize f (float) on 3 vars, coordinates clipped to >= 0."""
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] = max(0.0, p[i] + step * (1 + abs(p[i])))
        pts.append(p)
    vals = [f(p) for p in pts]
    hist = []
    for it in range(maxit):
        order = sorted(range(n + 1), key=lambda i: -vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        hist.append((it, vals[0], list(pts[0])))
        if abs(vals[0] - vals[-1]) < tol:
            break
        centroid = [sum(p[i] for p in pts[:-1]) / n for i in range(n)]
        worst = pts[-1]
        refl = [max(0.0, centroid[i] + 1.0 * (centroid[i] - worst[i])) for i in range(n)]
        fr = f(refl)
        if fr > vals[0]:
            exp = [max(0.0, centroid[i] + 2.0 * (centroid[i] - worst[i])) for i in range(n)]
            fe = f(exp)
            pts[-1], vals[-1] = (exp, fe) if fe > fr else (refl, fr)
        elif fr > vals[-2]:
            pts[-1], vals[-1] = refl, fr
        else:
            con = [centroid[i] + 0.5 * (worst[i] - centroid[i]) for i in range(n)]
            fc = f(con)
            if fc > vals[-1]:
                pts[-1], vals[-1] = con, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [(pts[i][j] + pts[0][j]) / 2 for j in range(n)]
                    vals[i] = f(pts[i])
    order = sorted(range(n + 1), key=lambda i: -vals[i])
    return pts[order[0]], vals[order[0]], hist


if __name__ == "__main__":
    setup()
    out = {"stage": "S2+S3", "gamma_paper": GAMMA_PAPER, "gamma_gateA": GAMMA_GATEA}
    t0 = time.time()

    # ---- reference: the paper point ----
    g0, i0 = gamma_head(S3_P, S5_P, 0.0)
    print(f"reference (paper, s7=0): gamma_head = {g0:.15f}")
    print(f"   b1 = {i0['b1']:.15f}  head = {i0['head']:.6e}  V = {i0['V']:.10f}")
    print(f"   b3 = {i0['b3']:+.3e}  b5 = {i0['b5']:+.3e}  b7 = {i0['b7']:+.3e}  b11 = {i0['b11']:+.3e}")
    out["reference"] = dict(s=(S3_P, S5_P, 0.0), gamma_head=g0, **{k: v for k, v in i0.items() if k != 'g_arb'})

    # ---- S2a: b3 = b5 = 0 manifold continued in s7 (EVERY node reported) ----
    S7_GRID = [0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01,
               0.015, 0.02, 0.03, 0.05]
    print("\nS2a: b3 = b5 = 0 manifold, s7 continuation (all nodes reported)")
    rows = []
    xstart = [S3_P, S5_P]
    for s7 in S7_GRID:
        sol, resid = newton_solve(s7, which=(3, 5), start=xstart)
        if sol is None:
            print(f"  s7 = {s7}: newton FAILED")
            rows.append(dict(s7=s7, status="newton-failed"))
            continue
        xstart = [sol[0], sol[1]]
        g, info = gamma_head(*sol)
        flag = "BEATS-gateA" if g > GAMMA_GATEA else ("beats-paper" if g > GAMMA_PAPER else "-")
        print(f"  s7 = {s7:<7}: s3 = {sol[0]:.8f} s5 = {sol[1]:.8f}  V = {info['V']:.9f}  "
              f"head = {info['head']:.6e}  gamma_head = {g:.15f}  {flag}")
        rows.append(dict(s7=s7, s3=sol[0], s5=sol[1], resid=resid, gamma_head=g,
                         head=info['head'], b1=info['b1'], V=info['V'],
                         b7=info['b7'], b9=info['b9'], b11=info['b11'], flag=flag))
    out["S2a"] = rows

    # ---- S2b: triple annihilation b3 = b5 = b7 = 0 ----
    print("\nS2b: triple annihilation b3 = b5 = b7 = 0")
    best_triple = None
    for s7_0 in (0.001, 0.005, 0.02, 0.05):
        sol, resid = newton_solve(s7_0, which=(3, 5, 7), start=[S3_P, S5_P, s7_0])
        if sol is None:
            print(f"  start s7={s7_0}: FAILED")
            continue
        g, info = gamma_head(*sol)
        print(f"  start s7={s7_0}: s = ({sol[0]:.8f}, {sol[1]:.8f}, {sol[2]:.8f})  "
              f"resid = {resid:.2e}  V = {info['V']:.9f}  head = {info['head']:.6e}  "
              f"gamma_head = {g:.15f}")
        out.setdefault("S2b", []).append(dict(start_s7=s7_0, s=sol, resid=resid,
                                             gamma_head=g, head=info['head'], V=info['V']))
        if best_triple is None or g > best_triple[1]:
            best_triple = (sol, g)

    # ---- S3: free Nelder-Mead ----
    print("\nS3: free Nelder-Mead on gamma_head")
    def obj(p):
        s3, s5, s7 = [max(0.0, v) for v in p]
        return gamma_head(s3, s5, s7)[0]
    starts = [(S3_P, S5_P, 0.0)]
    if best_triple:
        starts.append(tuple(best_triple[0]))
    best_nm = None
    for st in starts:
        xb, vb, hist = nelder_mead(obj, list(st), step=2e-4, maxit=250)
        print(f"  from {tuple(round(v,8) for v in st)}: best gamma_head = {vb:.15f} at "
              f"({xb[0]:.9f}, {xb[1]:.9f}, {xb[2]:.9f})")
        out.setdefault("S3", []).append(dict(start=list(st), best=list(xb), gamma_head=vb,
                                             n_steps=len(hist)))
        if best_nm is None or vb > best_nm[1]:
            best_nm = (xb, vb)

    print(f"\n[{time.time()-t0:.0f}s]")
    print(f"paper-point gamma_head = {g0:.15f}")
    if best_nm:
        print(f"best found gamma_head   = {best_nm[1]:.15f}  at s7 = {best_nm[0][2]:.9f}")
        print(f"delta vs paper point    = {best_nm[1] - g0:+.6e}")
        print(f"beats gamma_gateA ({GAMMA_GATEA})? {best_nm[1] > GAMMA_GATEA}")
    with open("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_valley.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("wrote scratch/logs/gc_valley.json")
