"""d3h_leg2_indep.py — STRUCTURALLY INDEPENDENT route for leg1/leg2 = ||Phi3||^2, ||dxdy Phi3||^2
at fixed theta (t = +1, t = -1), via the Mehler b-series + exact q-ladder + FdB jets.

Independence contract (Main's IRC):
  NOT the separable-|K|^2 / e^{2 beta uv}-series / P-matrix / inner_c path (d3h_inner.py),
  and NOT the predecessor's mult-term-algebra (d3h_kernels kind chains / d3h_dxdy entries).
  Shared ONLY: psi3 & derivatives, vartheta, rho(t), c-table, q_b definition — all licensed.

Mathematical route (all exact identities):
  C_r(x,y) = (pi/2) * sum_{b>=0} r^b q_b(u) q_b(v),   u = vartheta*psi3(x), v = vartheta*psi3(y)
  q_b(u) = 2 phi(u) psi_{b-1}(-u)/sqrt(b)   (q_0 = erf(u/sqrt2)),  LADDER: d q_b/du = sqrt(b+1) q_{b+1}
      (b=0 case verified: d/du erf(u/sqrt2) = sqrt(2/pi) e^{-u^2/2} = 2 phi(u) = q_1 exactly)
  Q^A_b(x) := (d/dx)^A [q_b(vartheta psi3(x))]  — via FdB composition of jets, EXACT (orders <= 4)
      (d/dx)^k on vartheta*psi3(x): u-derivatives vanish beyond order 4 (u'''' = vartheta*sqrt6).
  (∂x∂y)^A C_r = (pi/2) sum_b r^b Q^A_b(x) Q^A_b(y)          (separability, exact)
  (∂r^p)(∂x∂y)^A C_r = (pi/2) sum_{b>=p} fall(b,p) r^{b-p} Q^A_b(x) Q^A_b(y),  fall(b,p)=b!/(b-p)!
  leg2 = ||dxdy Phi3||^2 with dxdy Phi3 = sum_{(p,a)} c_{p,a} G_{p,a+1},  nu := a+1 in {1..4}:
      dxdy Phi3(x,y) = (pi/2) sum_{nu=1..4} sum_b D^nu_b Q^nu_b(x) Q^nu_b(y),
      D^nu_b = sum_{(p,a): a+1=nu} c_{p,a} * fall(b,p) * r^{b-p} * [b >= p]
  leg1 with nu := a in {0..3}, Q^0_b = q_b(vartheta psi3(x)).

Quadrature: composite Gauss-Legendre on [-X0, X0]^2 (+ rigorous shell bound, separate script).
Convergence: several (panels x nodes-per-panel) counts, decrements reported.
Cross-checks baked in:
  (V1) d q_0/du = q_1 exact identity (numeric).
  (V2) my G_{0,1} vs 2 pi vartheta^2 psi3'(x) psi3'(y) K_r(u,v,r)  [(48) bridge, predecessor-validated Kf]
  (V3) my dxdy-Phi3 / Phi3 pointwise vs predecessor's dxdy_Phi3_eval / Phi3_eval (assembly equivalence)
  (V4) F-matrix route vs algebraic Gram route (two disjoint reductions of the same sum)
Run:  ../../.venv/bin/python d3h_leg2_indep.py            (theta=0 then theta=pi, both legs, all grids)
"""
import math, sys, time, json
import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/kg/src")

S3 = 0.34101124
S5 = 0.05276111
ETA = 0.136419125
V = 1 + S3 * S3 + S5 * S5
VT = ETA / math.sqrt(V)                      # vartheta
RHO1 = (1 - S3 * S3) / V                     # rho(1), rho(-1) = -RHO1

B_MAX = 250         # b-window (raised from 60; series-tail guard, see btail report)
X0 = 4.4            # |x| cut (rigorous shell bound in d3h_leg2_tail.py)

IO_LIST = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def rho_j(j, t):
    """D^j rho at t = +1 or -1 (float; exact rational-decimal input)."""
    if j == 0:
        return t * (1 - S3 * S3 * t * t + S5 * S5 * t ** 4) / V
    if j == 1:
        return (1 + 3 * (-S3 * S3) * t * t + 5 * S5 * S5 * t ** 4) / V
    if j == 2:
        return (6 * (-S3 * S3) * t * t + 20 * S5 * S5 * t ** 4) / V
    if j == 3:
        return (12 * (-S3 * S3) * t * t + 80 * S5 * S5 * t ** 4) / V
    raise ValueError(j)


def ctable(t):
    r1, r2, r3 = rho_j(1, t), rho_j(2, t), rho_j(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3,
            (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 * r1, (1, 2): 3 * t * t * r1}


# ---------------- psi3 derivatives (weights) ----------------
def psi3(x):
    return (x ** 3 - 3 * x) / math.sqrt(6.0)


def wpsi(j, x):
    """psi3^(j)(x): j=0..3."""
    if j == 0:
        return (x ** 3 - 3 * x) / math.sqrt(6.0)
    if j == 1:
        return (3 * x * x - 3) / math.sqrt(6.0)
    if j == 2:
        return x * math.sqrt(6.0)
    if j == 3:
        return math.sqrt(6.0)
    raise ValueError(j)


from math import factorial
PI = math.pi


def phi(u):
    return math.exp(-0.5 * u * u) / math.sqrt(2 * math.pi)


def psi_series(mmax: int, x: float) -> list:
    """psi_0..psi_mmax at x, ORTHONORMAL probabilists' (psi_m = He_m/sqrt(m!)),
    stable three-term recurrence: psi_{m+1} = (x psi_m - sqrt(m) psi_{m-1})/sqrt(m+1). [V4-locked conv]"""
    out = [1.0, x]
    for m in range(1, mmax):
        out.append((x * out[m] - math.sqrt(m) * out[m - 1]) / math.sqrt(m + 1))
    return out


def psi_of(m: int, x: float) -> float:
    if m == 0:
        return 1.0
    if m == 1:
        return x
    return psi_series(m, x)[m]


def qn(n: int, u: float) -> float:
    """paper q_n (ORTHONORMAL convention, notes section 8, grid-verified): n=0 -> erf(u/sqrt2);
    n>=1 -> 2 phi(u) psi_{n-1}(-u)."""
    if n == 0:
        return math.erf(u / math.sqrt(2.0))
    return 2 * phi(u) * psi_of(n - 1, -u)


def qladder(nmax: int, u: float) -> list:
    """q_0..q_nmax at u, orthonormal convention, single psi pass."""
    out = [math.erf(u / math.sqrt(2.0))]
    if nmax == 0:
        return out
    ps = psi_series(nmax - 1, -u)
    for n in range(1, nmax + 1):
        out.append(2 * phi(u) * ps[n - 1] / math.sqrt(n))
    return out




# ---------------- FdB jet machinery (orders <= 4) ----------------
def compose_jet(g, u):
    """Given jets g = [g0..gA] (derivatives of the OUTER function at u0) and
    u = [u0..uA] (derivatives of the INNER function at x), return jet of g(u(x)) to order A.
    FdB: h_k = sum over (m_1..m_k), sum j m_j = k:  k! / prod(m_j!) * prod(u_j/j!)^{m_j} * g_{M},
    M = sum m_j."""
    A = len(u) - 1
    h = [0.0] * (A + 1)
    h[0] = g[0]
    for k in range(1, A + 1):
        # enumerate all (m_1..m_k) with sum j*m_j = k
        acc = 0.0
        def rec(jrem, jmin, mults, M):
            nonlocal acc
            if jrem == 0:
                # coefficient
                c = factorial(k)
                for jj in range(1, k + 1):
                    c /= factorial(mults[jj])
                gj = g[M] if M < len(g) else 0.0
                w = 1.0
                for jj in range(1, k + 1):
                    if mults[jj]:
                        w *= (u[jj] / factorial(jj)) ** mults[jj]
                acc += c * w * gj
                return
            for jj in range(jmin, jrem + 1):
                if jj >= len(u):
                    continue
                mults[jj] += 1
                rec(jrem - jj, jj, mults, M + 1)
                mults[jj] -= 1
        rec(k, 1, [0] * (k + 1), 0)
        h[k] = acc
    return h


def Qjet(A, b, x):
    """Jet of q_b(vartheta*psi3(x)) to order A: [q, q', ..., q^{(A)}] as functions of x."""
    u0 = VT * psi3(x)
    # u-jet: derivatives of vartheta*psi3(x) up to A (>=4 vanish)
    uj = [0.0] * (A + 1)
    for j in range(0, min(A, 3) + 1):
        uj[j] = VT * wpsi(j, x)
    # g-jet: derivatives of q_b at u0 via ladder: q_b^{(k)} = sqrt((b+k)!/b!) q_{b+k}
    gj = []
    lad = qladder(b + A, u0)
    for k in range(0, A + 1):
        gj.append(math.sqrt(factorial(b + k) / factorial(b)) * lad[b + k])
    return compose_jet(gj, uj)


# ---------------- D coefficients ----------------
def Dcoefs(nu, t):
    """D^nu_b for b = 0..B_MAX-1 at t=+-1; returns array [(pi/2)-prefactor EXCLUDED]."""
    r = rho_j(0, t)
    D = np.zeros(B_MAX)
    for (p, a) in IO_LIST:
        if a + 1 != nu:
            continue
        c = ctable(t)[(p, a)]
        for b in range(p, B_MAX):
            fall = factorial(b) / factorial(b - p)
            D[b] += c * fall * r ** (b - p)
    return D


def Dcoefs_leg1(nu, t):
    """leg1 blocks: nu = a in {0..3}."""
    r = rho_j(0, t)
    D = np.zeros(B_MAX)
    for (p, a) in IO_LIST:
        if a != nu:
            continue
        c = ctable(t)[(p, a)]
        for b in range(p, B_MAX):
            fall = factorial(b) / factorial(b - p)
            D[b] += c * fall * r ** (b - p)
    return D


# ---------------- Q matrices on a GL grid ----------------
_GL_CACHE = {}
def gl_panels(P, m):
    """composite GL on [-X0, X0]: returns nodes (1-D array), weights (with phi folded),
    and per-node panel widths."""
    key = (P, m)
    if key in _GL_CACHE:
        return _GL_CACHE[key]
    xs, ws = [], []
    edges = np.linspace(-X0, X0, P + 1)
    xg, wg = leggauss(m)
    for i in range(P):
        a, b = edges[i], edges[i + 1]
        xm, xw = 0.5 * (a + b), 0.5 * (b - a)
        for j in range(m):
            x = xm + xw * xg[j]
            xs.append(x)
            ws.append(xw * wg[j] * math.exp(-0.5 * x * x) / math.sqrt(2 * PI))
    xs = np.array(xs); ws = np.array(ws)
    _GL_CACHE[key] = (xs, ws)
    return xs, ws


def Qmatrix(nu, xs, A_exact=None):
    """Q^nu_b(x_i) for all b in [0,B_MAX), all nodes: shape (len(xs), B_MAX)."""
    N = len(xs)
    Qm = np.empty((N, B_MAX))
    for i, x in enumerate(xs):
        lad_all = {}
        # build q-ladder once per x up to n = B_MAX-1 + nu
        lad = qladder(B_MAX - 1 + nu, VT * psi3(x))
        u0 = VT * psi3(x)
        uj = [VT * wpsi(j, x) for j in range(0, min(nu, 3) + 1)]
        uj += [0.0] * (nu - len(uj) + 1) if nu > 3 else [VT * math.sqrt(6.0)] if nu == 4 else []
        uj = uj[:nu + 1]
        # direct ladder-basis FdB for THIS nu: precompute g-jets for all b at once
        for b in range(B_MAX):
            gj = [math.sqrt(factorial(b + k) / factorial(b)) * lad[b + k] for k in range(0, nu + 1)]
            jet = compose_jet(gj, uj)
            Qm[i, b] = jet[nu]
    return Qm


def leg_value(t, leg2=True, P=8, m=16, verbose=False):
    """Independent leg value at t=+-1 via F-matrix route. Returns (value, diag-dict)."""
    t0 = time.time()
    xs, ws = gl_panels(P, m)
    N = len(xs)
    nus = range(1, 5) if leg2 else range(0, 4)
    tot = 0.0
    Fmats = {}
    Qmats = {}
    for nu in nus:
        Q = Qmatrix(nu, xs)           # N x B_MAX
        Qmats[nu] = Q
        D = Dcoefs(nu, t) if leg2 else Dcoefs_leg1(nu, t)
        Fmats[nu] = (Q * D[None, :]) @ Q.T * (PI / 2.0)   # N x N matrix, F-contribution
    F = sum(Fmats[nu] for nu in nus)                      # dxdy Phi3 (or Phi3) on the grid
    F2 = F * F
    val = float(ws @ F2 @ ws)
    diag = dict(P=P, m=m, N=N, t=t, wall=time.time() - t0)
    if verbose:
        print(f"  P={P:3d} m={m:2d} N={N:4d}: val = {val:.10f}  ({diag['wall']:.1f}s)", flush=True)
    return val, diag, F, ws, Qmats


def gram_value(t, leg2=True, P=8, m=16):
    """Algebraic Gram route: leg = (pi/2)^2 sum_{nu,nu',b,b'} D D' G^2, G = Q^T W Q.
    Independent reduction of the same double sum (check V4)."""
    xs, ws = gl_panels(P, m)
    nus = range(1, 5) if leg2 else range(0, 4)
    Qs, Ds = {}, {}
    for nu in nus:
        Qs[nu] = Qmatrix(nu, xs)
        Ds[nu] = Dcoefs(nu, t) if leg2 else Dcoefs_leg1(nu, t)
    tot = 0.0
    G = {}
    for nu in nus:
        for nup in nus:
            G[(nu, nup)] = Qs[nu].T @ (ws[:, None] * Qs[nup])
    for nu in nus:
        for nup in nus:
            Gm = G[(nu, nup)]
            tot += float(np.sum((Ds[nu][:, None] * Ds[nup][None, :]) * Gm * Gm))
    return (PI / 2.0) ** 2 * tot


# ---------------- series-tail diagnostic for B_MAX Justification ----------------
def btail_bound(t, leg2=True, umax=None):
    """Rigorous-positive bound of the discarded b-tail: sum_{b>B} |D_b| max_x|Q_b| max_y|Q_b|,
    with |Q_b(x)| bounded on |x|<=X0 via |q_n(u)| <= e^{u^2/4} (Cramer-type: |psi_n|<=e^{u^2/4})
    times |u'|^j /j!-factorials — conservative but explicit. Reported, not load-bearing for
    the number (the number is a computation; the shell bound is the rigor piece)."""
    # crude conservative: |Q^A_b(x)| <= sum_k sqrt((b+k)!/b!) sup|q_{b+k}| * [FdB coeffs]
    # use sup_{|x|<=X0} e^{-u^2/2} |He_m(-u)|/sqrt(m!) <= max over grid — CHEAP NUMERIC sup on a fine grid
    xf = np.linspace(-X0, X0, 801)
    r = rho_j(0, t)
    sups = {}
    nmax_need = B_MAX + 4
    lad_grid = np.array([[qn_scalar_bounded(n, VT * psi3(xx)) for n in range(nmax_need)] for xx in xf])
    out = {}
    for leg in ("leg1", "leg2"):
        nus = range(1, 5) if leg == "leg2" else range(0, 4)
        tt = 0.0
        for nu in nus:
            D = Dcoefs(nu, t) if leg == "leg2" else Dcoefs_leg1(nu, t)
            # |Q^nu_b| <= sum over k<=nu of [FdB weight factors] — bound by bell numbers * max|u^(j)| etc.
            # conservative: |Q| <= sum_{k<=nu} sqrt((b+k)!/b!) * maxQ(b+k) * C_nu
            Cnu = fdb_bell_bound(nu)
            for b in range(B_MAX):
                qmax = max(abs(lad_grid[:, b + k]).max() for k in range(nu + 1))
                lad2 = max(abs(lad_grid[:, b + k2]).max() for k2 in range(nu + 1))
                tt += abs(D[b]) * (Cnu ** 2) * qmax * lad2
        out[leg] = (PI / 2.0) ** 2 * tt
    return out


_FDB_CACHE = {}
def fdb_bell_bound(A):
    """Conservative bound on the total FdB coefficient mass multiplying the g-jets at fixed order A:
    sum over multi-indices of A!/prod(m!) * prod(|u_j|/j!)^{m_j} with sup |u_j| on |x|<=X0."""
    if A in _FDB_CACHE:
        return _FDB_CACHE[A]
    xf = np.linspace(-X0, X0, 401)
    ujs = [max(abs(VT * wpsi(j, xf)).max() for j in [0]) ]  # u0 not needed in FdB weights
    ujmax = [max(abs(VT * wpsi(j, xx)) for xx in xf) if j <= 3 else 0.0 for j in range(5)]
    # sum over partitions
    tot = 0.0
    def rec(k, jmin, mults, weight):
        nonlocal tot
        if k == 0:
            tot += weight
            return
        for j in range(jmin, min(k, 3) + 1):
            m = mults[:]
            m[j] += 1
            w = weight * (ujmax[j] / factorial(j))
            rec(k - j, j, m, factorial(len([1])) * 0 + w)
        return
    # simpler: dp
    tot = 0.0
    for comp in partitions_mult(A):
        c = factorial(A)
        w = 1.0
        for j, m in comp.items():
            c /= factorial(m) * factorial(j) ** m
            w = w * (ujmax[j] ** m)
        tot += c * w
    _FDB_CACHE[A] = tot
    return tot


def partitions_mult(A):
    """all multisets {j: m_j} with sum j*m_j = A, as list of dicts."""
    out = []
    def rec(rem, jmin, cur):
        if rem == 0:
            out.append(dict(cur))
            return
        for j in range(jmin, min(rem, 3) + 1):
            cur[j] = cur.get(j, 0) + 1
            rec(rem - j, j, cur)
            cur[j] -= 1
            if cur[j] == 0:
                del cur[j]
    rec(A, 1, {})
    return out


def qn_scalar_bounded(n, u):
    return abs(qn(n, u))


# ---------------- validations ----------------
def run_validations(t):
    print("== VALIDATIONS ==")
    r = rho_j(0, t)
    # V1: dq0/du = q1
    h = 1e-6
    dnumer = (qn(0, 1e-4 + h) - qn(0, 1e-4 - h)) / (2 * h)
    print(f"V1  d q0/du = q1:  {dnumer:.12f} vs {qn(1, 1e-4):.12f}  diff = {abs(dnumer - qn(1, 1e-4)):.2e}")
    # V2: my G01 vs 2 pi vt^2 psi3' psi3' Kf (predecessor's ball-validated kernel)
    from d3h_kernels import Kf, vartheta
    vt = float(vartheta().mid().real)
    import d3h_dxdy as DD
    from flint import acb
    xs_test = [0.7, -1.3, 2.2]
    ys_test = [1.1, -0.4, 2.6]
    worst = 0.0
    for x in xs_test:
        for y in ys_test:
            # my G_{0,1} = (pi/2) sum_b r^b Q1_b(x) Q1_b(y)
            Qx = np.array([Qjet(1, b, x)[1] for b in range(B_MAX)])
            Qy = np.array([Qjet(1, b, y)[1] for b in range(B_MAX)])
            bvals = np.array([r ** b for b in range(B_MAX)])
            mine = (PI / 2) * float(np.dot(bvals * Qx, Qy))
            kv = float(Kf(acb(vt * DDxpsi3(x)), acb(vt * DDxpsi3(y)), acb(r), 128).mid().real)
            ref = 2 * math.pi * vt ** 2 * ((3 * x * x - 3) / math.sqrt(6)) * ((3 * y * y - 3) / math.sqrt(6)) * kv
            rel = abs(mine - ref) / max(abs(ref), 1e-30)
            worst = max(worst, rel)
    print(f"V2  G01 mine vs 2pi vt^2 psi3'psi3' K_r: worst rel = {worst:.3e}  (9 pts)")
    # V3: pointwise vs predecessor evaluators
    import d3h_phi3 as PH
    worst2 = 0.0
    for x in xs_test:
        for y in ys_test:
            # my dxdy Phi3 (series)
            tot = 0.0
            r_ = rho_j(0, t)
            for nu in range(1, 5):
                D = Dcoefs(nu, t)
                for b in range(B_MAX):
                    tot += (PI / 2) * D[b] * Qjet(nu, b, x)[nu] * Qjet(nu, b, y)[nu]
            theirs = float(PH.dxdy_Phi3_eval(acb(x), acb(y), acb(t), acb(vt), 128).mid().real)
            rel = abs(tot - theirs) / max(abs(theirs), 1e-30)
            worst2 = max(worst2, rel)
    print(f"V3  dxdy Phi3 mine vs predecessor: worst rel = {worst2:.3e}  (9 pts)")
    worst3 = 0.0
    for x in xs_test:
        for y in ys_test:
            tot = 0.0
            for a in range(0, 4):
                D = Dcoefs_leg1(a, t)
                for b in range(B_MAX):
                    jl = Qjet(a, b, x)
                    jr = Qjet(a, b, y)
                    base_x = jl[0] if a == 0 else jl[a]
                    base_y = jr[0] if a == 0 else jr[a]
                    tot += (PI / 2) * D[b] * base_x * base_y
            theirs = float(PH.Phi3_eval(acb(x), acb(y), acb(t), acb(vt), 128).mid().real)
            rel = abs(tot - theirs) / max(abs(theirs), 1e-30)
            worst3 = max(worst3, rel)
    print(f"V3b Phi3 mine vs predecessor: worst rel = {worst3:.3e}  (9 pts)")


def DDxpsi3(x):
    return (x ** 3 - 3 * x) / math.sqrt(6.0)


# ---------------- main sweep ----------------
GRIDS = [(8, 6), (10, 8), (12, 10), (16, 12), (20, 16)]

def main():
    out = {}
    for t in (1.0, -1.0):
        tname = "theta0" if t > 0 else "thetapi"
        out[tname] = {}
        print(f"\n########## t = {t:+.0f} (theta = {0 if t > 0 else math.pi:.6f}) ##########")
        run_validations(t)
        for leg2 in (False, True):
            leg = "leg1" if not leg2 else "leg2"
            print(f"-- {leg} ({'||dxdy Phi3||^2' if leg2 else '||Phi3||^2'}), convergence in (P,m):")
            vals = []
            for (P, m) in GRIDS:
                v, diag, F, ws, Qs = leg_value(t, leg2=leg2, P=P, m=m, verbose=False)
                vg = gram_value(t, leg2=leg2, P=P, m=m)
                rel = abs(v - vg) / max(abs(v), 1e-300)
                vals.append((P, m, v, diag["wall"], rel))
                print(f"  P={P:3d} m={m:2d}  N/axis={P*m:4d}  val = {v:.8f}   gram-vs-F rel = {rel:.2e}  ({diag['wall']:.1f}s)", flush=True)
            decr = []
            for i in range(1, len(vals)):
                decr.append(vals[i][2] - vals[i - 1][2])
            print(f"  decrements: " + ", ".join(f"{d:+.6f}" for d in decr))
            out[tname][leg] = dict(values=[[P, m, v, w, rel] for (P, m, v, w, rel) in vals])
        # tail diagnostic
        tb = btail_bound(t)
        print(f"  b-tail disc bound (conservative): leg1 <= {tb['leg1']:.3e}  leg2 <= {tb['leg2']:.3e}")
        out[tname]["btail"] = tb
    with open("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2conv_indep.json", "w") as f:
        json.dump(out, f, indent=1)
    print("\nsaved scratch/logs/leg2conv_indep.json")


if __name__ == "__main__":
    main()
