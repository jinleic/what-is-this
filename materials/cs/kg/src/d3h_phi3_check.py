"""phi3_check — AUTHORISED: independent Phi3 derivation from paper eq (52) recursion.
Derives Phi0 = C_{rho(t)} and Phi_{k+1} = D_t Phi_k - t dxdy Phi_k step by step with an
INDEPENDENT c-table read from the paper text (line numbers cited), comparing pointwise
against the current d3h_phi3 assembly at theta = 0, 3pi/4, pi.

Derivation sketch (all symbolic, sympy-checked where feasible, exact rational coefficients):
  Phi0(t;x,y) = C_{rho(t)}(x,y)
  D_t = t d/dt.
  Phi1 = D_t Phi0 - t dxdy Phi0
       = rho1 * dC/dr|_{r=rho} - t dxdy C_{rho}
       = rho1 * G_{1,0} - t * G_{0,1}
  Phi2 = D_t Phi1 - t dxdy Phi1
  Phi3 = D_t Phi2 - t dxdy Phi2
Each G_{p,a} = [dr^p dxdy^a C_r]_{r=rho(t)} with C_r = (pi/2) sum r^b q_b q_b (b-series) — the
RULE here: products D_t[G] and dxdy[G] follow the product rule through BOTH the coefficient
functions (rho_j) and the kernel parameter r = rho(t).

c-table from paper lines 1835-1852 (verbatim):
  (1,0): rho3              (0,1): -t
  (2,0): 3 rho1 rho2       (0,2): 3 t^2
  (3,0): rho1^3            (0,3): -t^3
  (1,1): -3t(rho1+rho2)    (2,1): -3t rho1^2    (1,2): 3t^2 rho1
  with rho_j := D_t^j rho   (paper line 1824: "Write rho_j = D_t^j rho")
"""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import mpmath
from mpmath import mp, mpf, mpc
mp.dps = 50

# rho(t) and its D_t derivatives, exactly as decimals (all inputs exact terminating rationals):
S3V = mpf('0.34101124')
S5V = mpf('0.05276111')
V = 1 + S3V * S3V + S5V * S5V
C1 = 1 / V
C3 = -(S3V * S3V) / V
C5 = (S5V * S5V) / V
VT = mpf('0.136419125') / mpmath.sqrt(V)


def rho(t):
    return C1 * t + C3 * t ** 3 + C5 * t ** 5


def rho1(t):
    return C1 + 3 * C3 * t ** 2 + 5 * C5 * t ** 4


def rho2(t):
    return 6 * C3 * t ** 2 + 20 * C5 * t ** 4


def rho3(t):
    return 12 * C3 * t ** 2 + 80 * C5 * t ** 4


# ---------------- C_r derivatives via b-series (converges for |r| < 1) ----------------
def qq(x, y, r, bmax=900):
    """sum_{b>=0} r^b q_b(u(x)) q_b(v(y)); u = VT*psi3."""
    lx = qlad(x, bmax)
    ly = qlad(y, bmax)
    tot = mpf(0)
    for b in range(bmax):
        tot += (r ** b) * lx[b] * ly[b]
    return tot


def qlad(x, bmax):
    out = [mpmath.erf(VT * (x ** 3 - 3 * x) / mpmath.sqrt(6) / mpmath.sqrt(2))]
    u = VT * (x ** 3 - 3 * x) / mpmath.sqrt(6)
    ph = mpmath.e ** (-u * u / 2) / mpmath.sqrt(2 * mpmath.pi)
    # psi_{m}(-u), m = 0..bmax-2
    p0, p1 = mpf(1), -u
    out.append(2 * ph * p1 / mpmath.sqrt(1))
    for n in range(2, bmax):
        # psi_n(-u) = ((-u) psi_{n-1} - sqrt(n-1) psi_{n-2})/sqrt(n)
        p0, p1 = p1, ((-u) * p1 - mpmath.sqrt(n - 1) * p0) / mpmath.sqrt(n)
        out.append(2 * ph * p1 / mpmath.sqrt(n))
    return out


def G(p, a, x, y, t):
    """[dr^p dxdy^a C_r](x,y) at r=rho(t), b-series form, mpmath."""
    r = rho(t)
    # dxdy^a: apply the derivative chain symbolically via the Q-family jets? To stay
    # INDEPENDENT of prior code, build the b-series derivative directly:
    # dxdy^a [r^b q_b(u(x)) q_b(u(y))] = r^b Qb_a(x) Qb_a(y)
    # Qb_a via FdB jets of q_b(u(x)) — same formulas as validated Qjet but local copy.
    tot = mpf(0)
    hr = mpf('0.02')
    if p == 0:
        return dxdy_series(a, x, y, r)
    # r-derivative: use falling factorial identity — d^p/dr^p r^b = b!/(b-p)! r^{b-p}
    lx = qlad(x, 1100)
    ly = qlad(y, 1100)
    for b in range(p, 1100):
        fall = mpf(1)
        for j in range(p):
            fall *= (b - j)
        tot += fall * (r ** (b - p)) * dxdy_part(a, b, x, y, lx, ly)
    return tot


def dxdy_series(a, x, y, r):
    lx = qlad(x, 1100)
    ly = qlad(y, 1100)
    tot = mpf(0)
    for b in range(1100):
        tot += (r ** b) * dxdy_part(a, b, x, y, lx, ly)
    return tot


def dxdy_part(a, b, x, y, lx, ly):
    """d^a/dx^a d^a/dy^a of q_b(u(x)) q_b(u(y)) at (x,y): product of 1-D jets.
    Q^A_b built via FdB with local ladder (independent copy of validated machinery)."""
    u = VT * (x ** 3 - 3 * x) / mpmath.sqrt(6)
    v = VT * (y ** 3 - 3 * y) / mpmath.sqrt(6)
    # u-jet up to order a
    uj = [VT * mpmath.sqrt(6) * (m * 0) for m in range(0)]
    uj = [VT * ((3 * x * x - 3) / mpmath.sqrt(6)), VT * (x * mpmath.sqrt(6)),
          VT * mpmath.sqrt(6), VT * 0]
    # 0th, 1st, 2nd, 3rd derivative of u(x) = VT*psi3(x):
    uj = [VT * (x ** 3 - 3 * x) / mpmath.sqrt(6),
          VT * (3 * x * x - 3) / mpmath.sqrt(6),
          VT * x * mpmath.sqrt(6),
          VT * mpmath.sqrt(6),
          mpf(0)]
    vj = [VT * (y ** 3 - 3 * y) / mpmath.sqrt(6),
          VT * (3 * y * y - 3) / mpmath.sqrt(6),
          VT * y * mpmath.sqrt(6),
          VT * mpmath.sqrt(6),
          mpf(0)]
    # g-jet: q_b derivatives at u: q_b^{(k)}(u) = sqrt((b+k)!/b!) q_{b+k}(u)
    def qn_at(n, uu, lad_at):
        # evaluate q_n AT u uses its own ladder at uu — must build separate ladder per argument.
        return lad_at(uu, n)
    return jet_pair(a, x, y, uj, vj, b)


def lad_at(u, n):
    """q_n(u)."""
    if n == 0:
        return mpmath.erf(u / mpmath.sqrt(2))
    ph = mpmath.e ** (-u * u / 2) / mpmath.sqrt(2 * mpmath.pi)
    p0, p1 = mpf(1), -u
    if n == 1:
        return 2 * ph * p1
    for m in range(2, n + 1):
        p0, p1 = p1, ((-u) * p1 - mpmath.sqrt(m - 1) * p0) / mpmath.sqrt(m)
    return 2 * ph * p1 / mpmath.sqrt(n)


def jet_pair(A, x, y, uj, vj, b):
    """d^A/dx^A q_b(u(x)) * d^A/dy^A q_b(u(y)) via FdB (order A on each side)."""
    def side(z, zj):
        # h_A = sum over partitions of A
        u0 = zj[0]
        tot = mpf(0)
        # g_j = q_b^{(j)}(u0) = sqrt((b+j)!/b!) q_{b+j}(u0)
        gs = [mpmath.sqrt(mpmath.factorial(b + j) / mpmath.factorial(b)) * lad_at(u0, b + j)
              for j in range(A + 1)]
        # FdB enumerate
        h = [mpf(0)] * (A + 1)
        h[0] = gs[0]
        for k in range(1, A + 1):
            # multi-index enumeration
            def rec(rem, jmin, mults, M):
                if rem == 0:
                    c = mpmath.factorial(k)
                    w = mpf(1)
                    for jj in range(1, k + 1):
                        c /= mpmath.factorial(mults[jj])
                        w *= (zj[jj] / mpmath.factorial(jj)) ** mults[jj] if mults[jj] else mpf(1)
                    hj = gs[M] if M < len(gs) else mpf(0)
                    return c * w * hj
                s = mpf(0)
                for jj in range(jmin, min(rem, 3) + 1):
                    mults[jj] += 1
                    s += rec(rem - jj, jj, mults, M + 1)
                    mults[jj] -= 1
                return s
            h[k] = rec(k, 1, [0] * (k + 1), 0)
        return h[A]
    return side(x, uj) * side(y, vj)


# ---------------- Phi recursion (from (52), independent assembly) ----------------
def Phi_stage(k, x, y, t):
    """Phi_k at t, as {c-coefficient dict} accumulated over G_{p,a} rows: returns dict[(p,a)] = coeff."""
    T = {(0, 0): mpf(1)}  # Phi0 = G_{0,0} = C_r itself
    stages = [T]
    for _ in range(k):
        # D_t stage: every row (p,a) with coefficient c: D_t[c*G_{p,a}] =
        #   (D_t c)*G_{p,a} + c*rho1*(p+1)*G_{p+1,a}    [D_t G_{p,a} = rho1*(p+1) G_{p+1,a}]
        new = {}
        for (p, a), c in T.items():
            dc = dt_c(p, a, t)
            if dc != 0:
                new[(p, a)] = new.get((p, a), mpf(0)) + dc
            if p + 1 <= 3:
                new[(p + 1, a)] = new.get((p + 1, a), mpf(0)) + c * rho1(t) * (p + 1)
        # -t dxdy stage: (p,a) -> -t*c*G_{p,a+1}
        for (p, a), c in T.items():
            if a + 1 <= 3:
                new[(p, a + 1)] = new.get((p, a + 1), mpf(0)) - t * c
        T = new
        stages.append(T)
    return T


def dt_c(p, a, t):
    """D_t of the c_{p,a} polynomial: symbolic; returns value at t."""
    r1, r2, r3 = rho1(t), rho2(t), rho3(t)
    r4 = 24 * C3 * t ** 2 + 240 * C5 * t ** 4
    r5 = 48 * C3 * t ** 2 + 1200 * C5 * t ** 4
    r6 = 96 * C3 * t ** 2 + 3360 * C5 * t ** 4
    # D_t[c] by direct formula per row:
    if (p, a) == (0, 0):
        return mpf(0)
    if (p, a) == (1, 0):
        return r3
    if (p, a) == (2, 0):
        return 3 * (r2 * r2 + r1 * r4)
    if (p, a) == (3, 0):
        return 3 * r1 * r1 * rho2(t)
    if (p, a) == (0, 1):
        return -t
    if (p, a) == (0, 2):
        return 6 * t
    if (p, a) == (0, 3):
        return -3 * t ** 2
    if (p, a) == (1, 1):
        return -3 * (r1 + r2) - 3 * t * (r2 + r4)
    if (p, a) == (2, 1):
        return -3 * r1 ** 2 - 6 * t * r1 * r2
    if (p, a) == (1, 2):
        return 6 * t * r1 + 3 * t * t * r2
    raise KeyError((p, a))


def Dt_rho2(t):
    return 24 * C3 * t ** 2 + 240 * C5 * t ** 4


if __name__ == "__main__":
    for th in (0, 3 * mpmath.pi / 4, mpmath.pi):
        t = mpmath.cos(th) + mpmath.pi * 0  # real-t only works for th=0/pi; complex for others
        t = mpmath.exp(1j * th)
        x, y = mpf('0.7'), mpf('1.1')
        T = Phi_stage(3, x, y, t)
        tot = mpc(0)
        for (p, a), c in T.items():
            if c != 0:
                tot += c * G(p, a, x, y, t)
        # compare against d3h_phi3.Phi3_eval:
        import d3h_phi3 as PH
        import d3h_kernels as K
        from flint import acb
        vta = K.vartheta(160)
        thf = float(th)
        tc = complex(math.cos(thf), math.sin(thf))
        ref = complex(PH.Phi3_eval(acb(x), acb(y), acb(tc), vta, 192).mid().real)
        r = float(rho(t).real)
        mv = float(tot.real) + 0j + float(tot.imag)*1j
        rel = abs(mv - ref)/max(abs(ref), 1e-12)
        print(f"theta = {thf:.4f}: my-indep Phi3 = {mv:+.8f}   predecessor = {ref:+.8f}   ratio = {abs(mv/ref) if abs(ref) > 1e-12 else float('nan'):.6f}   rel = {rel:.3e}")
