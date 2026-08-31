"""septic_legs.py — Gate C(a): the septic (D={3,5,7}) S(theta) legs evaluator.

One family, one evaluator, parameters as function arguments:
  rho(t) = (t - s3^2 t^3 + s5^2 t^5 - s7^2 t^7)/V,  V = 1+s3^2+s5^2+s7^2   [eq (6)/(8), line 672/687]
  rho_j(t) = (t - 3^j s3^2 t^3 + 5^j s5^2 t^5 - 7^j s7^2 t^7)/V            [Euler D_t^j, lines 857/1826]
c-table rows exactly as printed (paper 1836-1852); eta fixed 0.136419125 (line 682) so
vartheta = eta/sqrt(V). s7=0 with the paper decimals reproduces the corrected quintic
assembly (anchor: avg(S)=123.377609, avg(leg1)=1.939156, avg(leg2)=30.359613,
B3=14.2459830).

Machinery: same symbolic entries as d3h_sweep49_asm (d3h_phi3.G_pa_terms,
d3h_kernels.mult_terms, d3h_dxdy.weight_eval) but parameterized — no module-level s3/s5.

All values COMPUTATIONAL-EVIDENCE (float64) unless run through gate_c_head (certified).
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

ETA = 0.136419125
# vartheta convention (pre_statement Addendum 3): FIXED at the paper's value so that the
# certified A-grid (computed at this vartheta) stays valid across the septic sweep.
# vartheta_paper = eta/sqrt(V_paper), V_paper = 1 + 0.34101124^2 + 0.05276111^2.
VT_PAPER = 0.136419125 / math.sqrt(1.0 + 0.34101124 ** 2 + 0.05276111 ** 2)

# GL grid — identical to d3h_sweep49_asm.py (window 6.0, 18 panels x 24 nodes)
XW = 6.0
P = 18
GM = 24
_edges = np.linspace(-XW, XW, P + 1)
_xg, _wg = leggauss(GM)
nl, wl = [], []
for _i in range(P):
    _a, _b = _edges[_i], _edges[_i + 1]
    _xm, _xw = 0.5 * (_a + _b), 0.5 * (_b - _a)
    _nds = _xm + _xw * _xg
    nl.append(_nds)
    wl.append(_xw * _wg * np.exp(-0.5 * _nds ** 2) / math.sqrt(2 * math.pi))
XS = np.concatenate(nl)
WS = np.concatenate(wl)
WIJ = np.outer(WS, WS)

IO = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def scheme(s3, s5, s7=0.0, vt=None):
    """Derived constants. vartheta FIXED at VT_PAPER by default (Addendum 3 convention);
    pass vt to override (e.g. vt='eta' reproduces the eta-fixed convention)."""
    V = 1.0 + s3 * s3 + s5 * s5 + s7 * s7
    if vt is None:
        VT = VT_PAPER
    elif vt == "eta":
        VT = ETA / math.sqrt(V)
    else:
        VT = float(vt)
    return V, VT


def rho0(t, s3, s5, s7=0.0, vt=None):
    V = scheme(s3, s5, s7, vt)[0]
    return (t - s3 * s3 * t ** 3 + s5 * s5 * t ** 5 - s7 * s7 * t ** 7) / V


def rhoj(j, t, s3, s5, s7=0.0, vt=None):
    """Euler D_t^j rho (paper lines 857/1826): exponent-m^j weights, sign per eq (6) sigma_m."""
    V = scheme(s3, s5, s7, vt)[0]
    return (t - (3 ** j) * s3 * s3 * t ** 3 + (5 ** j) * s5 * s5 * t ** 5
            - (7 ** j) * s7 * s7 * t ** 7) / V


def ctab(t, s3, s5, s7=0.0, vt=None):
    """Paper c-table verbatim (lines 1836-1852), rho_j = D_t^j rho."""
    r1 = rhoj(1, t, s3, s5, s7, vt)
    r2 = rhoj(2, t, s3, s5, s7, vt)
    r3 = rhoj(3, t, s3, s5, s7, vt)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3,
            (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2,
            (1, 2): 3 * t * t * r1}


def legs_at(t, s3, s5, s7=0.0, vt=None):
    """S-legs at unit-circle point t: (||Phi3||^2, ||dxdy Phi3||^2)."""
    from d3h_phi3 import G_pa_terms
    from d3h_kernels import mult_terms
    from d3h_dxdy import weight_eval
    from flint import acb

    V, VT = scheme(s3, s5, s7, vt)
    r = rho0(t, s3, s5, s7, vt)
    s = 1 - r * r
    sq_s = np.sqrt(complex(s))
    ux = VT * (XS ** 3 - 3 * XS) / math.sqrt(6)
    Km = np.exp(-(ux[:, None] ** 2 - 2 * r * ux[:, None] * ux[None, :] + ux[None, :] ** 2) / (2 * s)) \
        / (2 * math.pi * sq_s)
    ct = ctab(t, s3, s5, s7, vt)
    _wv = {}

    def wvec(wt):
        if wt not in _wv:
            _wv[wt] = np.array([float(weight_eval(wt, acb(x)).mid().real) for x in XS])
        return _wv[wt]

    F1 = np.zeros((len(XS), len(XS)), dtype=complex)
    F2 = np.zeros((len(XS), len(XS)), dtype=complex)
    for (p, a) in IO:
        cval = ct[(p, a)]
        if cval == 0:
            continue
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            Mmat = np.zeros((len(XS), len(XS)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (s ** (-pw)) if pw else 1.0
                up = (ux ** au) if au else np.ones_like(ux)
                vp = (ux ** av) if av else np.ones_like(ux)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F1 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
        for (wx, wy, kind, tp) in G_pa_terms(p, a + 1):
            Mmat = np.zeros((len(XS), len(XS)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (s ** (-pw)) if pw else 1.0
                up = (ux ** au) if au else np.ones_like(ux)
                vp = (ux ** av) if av else np.ones_like(ux)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F2 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
    l1 = float(np.real(np.sum(WIJ * np.abs(F1) ** 2)))
    l2 = float(np.real(np.sum(WIJ * np.abs(F2) ** 2)))
    return l1, l2


def assert_symmetries(S, th, tol=1e-10):
    """Tripwires T1 (pi-periodicity via S(0)==S(pi) endpoint check on [0,pi]), T2 (reflection
    about pi/2 — 24 independent equalities on the 49-pt grid), T3 (min at theta=pi/2)."""
    n = len(S)
    ok = True
    # T1 endpoints
    d0 = abs(S[0] - S[-1]) / max(1.0, abs(S[-1]))
    if d0 > tol:
        print(f"TRIPWIRE T1 (pi-periodicity) FAIL: |S(0)-S(pi)|/|S| = {d0:.3e}")
        ok = False
    # T2 reflection: S[k] == S[n-1-k]
    worst = 0.0
    for k in range(1, n // 2):
        d = abs(S[k] - S[n - 1 - k]) / max(1.0, abs(S[n - 1 - k]))
        worst = max(worst, d)
    if worst > tol:
        print(f"TRIPWIRE T2 (reflection) FAIL: worst rel = {worst:.3e}")
        ok = False
    # T3 min at midpoint
    mid = n // 2
    is_min = S[mid] <= S.min() + 1e-12 * max(1.0, abs(S.min()))
    if not is_min:
        print(f"TRIPWIRE T3 (min at pi/2) FAIL: S[mid] = {S[mid]:.9f} vs min = {S.min():.9f}")
        ok = False
    if ok:
        print(f"TRIPWIRES OK: T1 rel = {d0:.2e}, T2 worst rel = {worst:.2e}, "
              f"T3: min at pi/2 = {S[mid]:.6f}")
    return ok


def sweep(s3, s5, s7=0.0, vt=None, n=48, quiet=False, save_prefix=None):
    """49-pt uniform trapezoid stern sweep on [0, pi]; returns diagnostics dict."""
    import time
    th = np.linspace(0, math.pi, n + 1)
    L1 = np.zeros(n + 1)
    L2 = np.zeros(n + 1)
    t0 = time.time()
    for i, x in enumerate(th):
        L1[i], L2[i] = legs_at(complex(math.cos(x), math.sin(x)), s3, s5, s7, vt)
        if (not quiet) and (i % 8 == 0 or i == n):
            print(f"  th[{i:2d}/{n}] = {x:.4f}: leg1 = {L1[i]:12.5f}  leg2 = {L2[i]:12.5f}",
                  flush=True)
    h = th[1] - th[0]
    a1 = h * (L1[0] / 2 + L1[-1] / 2 + L1[1:-1].sum()) / math.pi
    a2 = h * (L2[0] / 2 + L2[-1] / 2 + L2[1:-1].sum()) / math.pi
    Sv = a1 + 4 * a2
    out = dict(s3=s3, s5=s5, s7=s7, avg_leg1=a1, avg_leg2=a2, avg_S=Sv,
               B3sq=(math.pi ** 2 / 6) * Sv, B3=math.sqrt((math.pi ** 2 / 6) * Sv),
               seconds=time.time() - t0)
    out["sym_ok"] = assert_symmetries(L1 + 4 * L2, th)
    out["min_S_theta"] = float(th[np.argmin(L1 + 4 * L2)])
    out["max_S"] = float((L1 + 4 * L2).max())
    out["S_at_0"] = float(L1[0] + 4 * L2[0])
    out["S_at_pi"] = float(L1[-1] + 4 * L2[-1])
    if save_prefix:
        np.save(save_prefix + "_th.npy", th)
        np.save(save_prefix + "_leg1.npy", L1)
        np.save(save_prefix + "_leg2.npy", L2)
    if not quiet:
        print(f"s(s3={s3}, s5={s5}, s7={s7}): avg(leg1) = {a1:.6f}  avg(leg2) = {a2:.6f}  "
              f"avg(S) = {Sv:.6f}  B3 = {out['B3']:.7f}  [{out['seconds']:.0f}s]")
    return out


if __name__ == "__main__":
    import sys
    # ---- ANCHOR: s7=0 at paper decimals ----
    s3, s5 = 0.34101124, 0.05276111
    V, VT = scheme(s3, s5)
    print(f"anchor params: V = {V!r}  vartheta = {VT!r}")
    print(f"rho(1) = {rho0(1.0, s3, s5)!r}  rho_2(1) = {rhoj(2, 1.0, s3, s5)!r}")
    out = sweep(s3, s5, 0.0, n=48)
    # checks (6-digit targets from the frozen corrected run)
    tgt = dict(avg_S=123.377609, avg_leg1=1.939156, avg_leg2=30.359613, B3=14.2459830)
    ok = True
    for k, v in tgt.items():
        got = out[k]
        rel = abs(got - v) / max(1.0, abs(v))
        stat = "OK " if rel < 5e-7 else "MISS"
        if rel >= 5e-7:
            ok = False
        print(f"  anchor {k}: got {got:.6f}  target {v}  rel {rel:.2e}  {stat}")
    print("ANCHOR", "PASS" if (ok and out["sym_ok"]) else "FAIL")
    sys.exit(0 if (ok and out["sym_ok"]) else 1)
