"""fd_oracle — independent high-precision evaluation of the nine G blocks and Phi3 at a point
by finite differences of the RAW q-series (49): C_r(x,y) = (pi/2) sum_b r^b q_b(u(x)) q_b(v(y)).
No Q-jets, no D-coefficient grouping, no mult algebra, no predecessor evaluators.
q_b from the locked orthonormal convention, ladder-validated.
"""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import mpmath
from mpmath import mp, mpf, mpc

mp.dps = 50
sq = mpmath.sqrt
mexp = mpmath.exp

ETA = mpf('0.136419125')
S3 = mpf('0.34101124')
S5 = mpf('0.05276111')
V = 1 + S3 * S3 + S5 * S5
VT = ETA / sq(V)
BNEW = 1400


def u_of(x):
    return VT * (x ** 3 - 3 * x) / sq(6)


_psi_cache = {}
def psi_at(m, x):
    key = (m, x)
    if key in _psi_cache:
        return _psi_cache[key]
    p0, p1 = mpf(1), x
    if m == 0:
        r = mpf(1)
    elif m == 1:
        r = x
    else:
        for k in range(1, m):
            p0, p1 = p1, (x * p1 - sq(k) * p0) / sq(k + 1)
        r = p1
    _psi_cache[key] = r
    return r


_erf_factor = {}
def q_ladder_x(x, bmax):
    """q_0..q_bmax at u = VT psi3(x)."""
    u = u_of(x)
    out = [mpmath.erf(u / sq(2))]
    ph = mexp(-u * u / 2) / sq(2 * mp.pi)
    for n in range(1, bmax + 1):
        out.append(2 * ph * psi_at(n - 1, -u) / sq(n))
    return out


def C_series(x, y, r, bmax=BNEW):
    lx = q_ladder_x(x, bmax)
    ly = q_ladder_x(y, bmax)
    tot = mpf(0)
    for b in range(bmax + 1):
        tot += (r ** b) * lx[b] * ly[b]
    return (mp.pi / 2) * tot


def fd(f, pts, h):
    """finite-difference combination: sum pts (offset, weight) / h^order."""
    return sum(w * f(off) for (off, w) in pts) / (h ** 0)


# central-difference stencils
ST1 = [(-1, -mpf(1) / 2), (1, mpf(1) / 2)]
ST2 = [(-1, mpf(1)), (0, mpf(-2)), (1, mpf(1))]
ST3 = [(-2, mpf(1) / 2), (-1, mpf(-1)), (1, mpf(1)), (2, mpf(-1) / 2)]

HR = mpf('0.02')
HX = mpf('0.05')


def G_block(p, a, x, y, r):
    """[dr^p dxdy^a C_r](x,y) at r via FD stencils on the raw series."""
    if a == 0 and p == 0:
        return C_series(x, y, r)
    # d/dx d/dy first (a>=1), then r-differentiate that
    if a == 0:
        g = lambda rr: C_series(x, y, rr)
        if p == 1:
            return sum(w * g(r + off * HR) for (off, w) in ST1) / HR
        if p == 2:
            return sum(w * g(r + off * HR) for (off, w) in ST2) / HR ** 2
        if p == 3:
            return sum(w * g(r + off * HR) for (off, w) in ST3) / HR ** 3
    if a == 1:
        gxy = lambda rr: (C_series(x + HX, y + HX, rr) - C_series(x + HX, y - HX, rr)
                          - C_series(x - HX, y + HX, rr) + C_series(x - HX, y - HX, rr)) / (4 * HX * HX)
        if p == 0:
            return gxy(r)
        if p == 1:
            return sum(w * gxy(r + off * HR) for (off, w) in ST1) / HR
        if p == 2:
            return sum(w * gxy(r + off * HR) for (off, w) in ST2) / HR ** 2
        if p == 3:
            return sum(w * gxy(r + off * HR) for (off, w) in ST3) / HR ** 3
    if a == 2:
        gxy = lambda xx, yy, rr: (
            C_series(xx + HX, yy + HX, rr) - C_series(xx + HX, yy - HX, rr)
            - C_series(xx - HX, yy + HX, rr) + C_series(xx - HX, yy - HX, rr)) / (4 * HX * HX)
        g2 = lambda xx, yy, rr: (
            gxy(xx + HX, yy + HX, rr) - gxy(xx + HX, yy - HX, rr)
            - gxy(xx - HX, yy + HX, rr) + gxy(xx - HX, yy - HX, rr)) / (4 * HX * HX)
        if p == 0:
            return g2(x, y, r)
        if p == 1:
            return sum(w * g2(x, y, r + off * HR) for (off, w) in ST1) / HR
        if p == 2:
            return sum(w * g2(x, y, r + off * HR) for (off, w) in ST2) / HR ** 2
        if p == 3:
            return sum(w * g2(x, y, r + off * HR) for (off, w) in ST3) / HR ** 3
    if a == 3:
        gxy = lambda xx, yy, rr: (
            C_series(xx + HX, yy + HX, rr) - C_series(xx + HX, yy - HX, rr)
            - C_series(xx - HX, yy + HX, rr) + C_series(xx - HX, yy - HX, rr)) / (4 * HX * HX)
        g2 = lambda xx, yy, rr: (
            gxy(xx + HX, yy + HX, rr) - gxy(xx + HX, yy - HX, rr)
            - gxy(xx - HX, yy + HX, rr) + gxy(xx - HX, yy - HX, rr)) / (4 * HX * HX)
        g3 = lambda xx, yy, rr: (
            g2(xx + HX, yy + HX, rr) - g2(xx + HX, yy - HX, rr)
            - g2(xx - HX, yy + HX, rr) + g2(xx - HX, yy - HX, rr)) / (4 * HX * HX)
        if p == 0:
            return g3(x, y, r)
        if p == 1:
            return sum(w * g3(x, y, r + off * HR) for (off, w) in ST1) / HR
        if p == 2:
            return sum(w * g3(x, y, r + off * HR) for (off, w) in ST2) / HR ** 2
        if p == 3:
            return sum(w * g3(x, y, r + off * HR) for (off, w) in ST3) / HR ** 3
    raise ValueError((p, a))


def rho_j(j, t):
    if j == 0:
        return (t - S3 ** 2 * t ** 3 + S5 ** 2 * t ** 5) / V
    if j == 1:
        return (1 + 3 * (-S3 ** 2) * t ** 2 + 5 * S5 ** 2 * t ** 4) / V
    if j == 2:
        return (6 * (-S3 ** 2) * t ** 2 + 20 * S5 ** 2 * t ** 4) / V
    if j == 3:
        return (12 * (-S3 ** 2) * t ** 2 + 80 * S5 ** 2 * t ** 4) / V


def ctable(t):
    r1, r2, r3 = rho_j(1, t), rho_j(2, t), rho_j(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3,
            (0, 1): -t, (0, 2): 3 * t ** 2, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2, (1, 2): 3 * t ** 2 * r1}


IO = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def Phi3_fd(x, y, t):
    c = ctable(t)
    r = rho_j(0, t)
    tot = mpf(0)
    for (p, a) in IO:
        tot += c[(p, a)] * G_block(p, a, x, y, r)
    return tot


def dxdyPhi3_fd(x, y, t):
    c = ctable(t)
    r = rho_j(0, t)
    tot = mpf(0)
    for (p, a) in IO:
        tot += c[(p, a)] * G_block(p, a + 1, x, y, r)
    return tot


if __name__ == "__main__":
    t = mpf(1)
    x = mpf('0') if len(sys.argv) < 2 else mpf(sys.argv[1])
    y = mpf('0') if len(sys.argv) < 3 else mpf(sys.argv[2])
    print(f"point (x,y) = ({float(x)}, {float(y)}), t = 1, r = rho(1) = {float(rho_j(0, t)):.6f}")
    print("block   G_(p,a)   [FD oracle on raw series]")
    phi = mpf(0)
    for (p, a) in IO:
        g = G_block(p, a, x, y, rho_j(0, t))
        c = ctable(t)[(p, a)]
        phi += c * g
        print(f"  ({p},{a}): G = {float(g):+.8f}   c = {float(c):+.6f}   c*G = {float(c * g):+.6f}")
    print(f"Phi3(0,0) [FD oracle]  = {float(phi):+.8f}")
    d2 = dxdyPhi3_fd(x, y, t)
    print(f"dxdyPhi3 [FD oracle]   = {float(d2):+.8f}")
