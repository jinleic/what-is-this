"""Validation of the Cayley-map Hilbert transform.

Section 1: rapidly decaying data (Poisson kernels). The naive `hilbert_line` is exact.
Section 2: the closed-form pair H[M_c]=M_s, H[M_s]=-M_c, against mpmath PV quadrature.
Section 3: algebraically decaying data -- the case the CCF profile actually needs, and
           the one that exposed the constant-offset defect. Naive vs tail-corrected.
"""

import os
import sys

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hilbert import (  # noqa: E402
    const_pow_tail, hilbert_line, hilbert_line_tail, make_grid,
    odd_tail, odd_tail_amplitude, odd_tail_hilbert, pow_tail,
)

mp.mp.dps = 40
LAM = 1.1807776628998          # arXiv:2509.14185 Fig.2, CCF stable
P = 1.0 / (1.0 + LAM)

ok = True


def check(name, err, tol):
    global ok
    good = bool(err < tol)
    ok = ok and good
    print(f"   [{'PASS' if good else 'FAIL'}] {name:48s} {err:.3e}  (tol {tol:.0e})")


def pv_hilbert(f, x):
    """H f at one point by symmetrised principal-value quadrature."""
    return float(mp.quad(lambda s: (f(x - s) - f(x + s)) / s,
                         [0, 1, 10, 100, mp.inf]) / mp.pi)


print("1. exact transform pairs, rapidly decaying data")
for n in (256, 4096, 65536):
    theta, y, half = make_grid(n)
    e1 = np.max(np.abs(hilbert_line(1.0 / (1.0 + y**2), half) - y / (1.0 + y**2)))
    e2 = np.max(np.abs(hilbert_line(y / (1.0 + y**2), half) + 1.0 / (1.0 + y**2)))
    a, b = 2.3, -0.7
    e3 = np.max(np.abs(hilbert_line(a / (a**2 + (y - b) ** 2), half)
                       - (y - b) / (a**2 + (y - b) ** 2)))
    check(f"N={n:<6d} Poisson / conjugate / shifted", max(e1, e2, e3), 1e-14)

theta, y, half = make_grid(8192)
f = 1.0 / (1.0 + y**2)
check("H^2 + I on decaying data",
      np.max(np.abs(hilbert_line(hilbert_line(f, half), half) + f)), 1e-14)

print("\n2. closed forms")
c_ref = float(-mp.quad(lambda x: x**2 * (1 + x**2) ** (-(P + 3) / 2),
                       [-mp.inf, 0, mp.inf]) / mp.pi)
check(f"C[pow_tail] = -B(3/2,p/2)/pi = {const_pow_tail(P):.12f}",
      abs(const_pow_tail(P) / c_ref - 1), 1e-13)

# The exact pair, checked pointwise against principal-value quadrature.
m_s = lambda s: (1 + s**2) ** (-P / 2) * mp.sin(P * mp.atan(s))
errs = [abs(pv_hilbert(m_s, t) - float(odd_tail_hilbert(np.array([t]), P)[0]))
        for t in (0.3, 1.0, 4.0, 30.0)]
check("H[M_s] = -M_c at y=0.3,1,4,30", max(errs), 1e-12)

theta, y, half = make_grid(4096)
amp = odd_tail_amplitude(P)
far = np.argmin(np.abs(y - 3000.0))
# The approach to the asymptote is O(1/y), not O(1/y^2): expanding
# sin(p(pi/2 - 1/y)) gives a relative defect cot(p*pi/2)*p/y ~ 1.6e-4 at y=3000.
check("odd_tail ~ sin(p*pi/2)|y|^{-p} far field",
      abs(odd_tail(y, P)[far] * y[far] ** P / amp - 1.0), 5e-4)
print("\n3. algebraically decaying data: naive vs tail-corrected")
print(f"   f = y(1+y^2)^-(p+1)/2,  p = {P:.10f}   (NOT the model: remainder is nonzero)")
print(f"   {'N':>7}  {'naive':>12}  {'tail-corrected':>16}")
f_pow = lambda s: s * (1 + s**2) ** (-(P + 1) / 2)
c_tail = 1.0 / amp                      # match the |y|^{-p} amplitude of f
rows = []
for n in (1024, 4096, 16384, 65536):
    theta, y, half = make_grid(n)
    g = pow_tail(y, P)
    idx = [np.argmin(np.abs(y - t)) for t in (0.5, 2.0, 20.0)]
    ref = [pv_hilbert(f_pow, y[i]) for i in idx]
    naive = hilbert_line(g, half)
    fixed = hilbert_line_tail(g, half, y, P, c_tail)
    en = max(abs(naive[i] - r) for i, r in zip(idx, ref))
    ef = max(abs(fixed[i] - r) for i, r in zip(idx, ref))
    rows.append((n, en, ef))
    print(f"   {n:>7}  {en:12.3e}  {ef:16.3e}")

# pow_tail - c*odd_tail still carries an O(|y|^{-p-1}) corner, so the FFT converges at
# O(N^{-(p+1)}) = O(N^{-1.46}).  Measured: 4.8e-6 -> 1.0e-8 over 64x.  Machine precision
# would need a multi-term tail; see the (1-iz)^{-q} family noted in hilbert.py.
check("tail-corrected error at N=65536", rows[-1][2], 5e-8)
check("tail-corrected beats naive by >1e5", rows[-1][2] / rows[-1][1], 1e-5)
# Self-consistency: on data that IS the model, the transform must be exact.
theta, y, half = make_grid(2048)
check("exact on f = odd_tail (remainder identically 0)",
      np.max(np.abs(hilbert_line_tail(odd_tail(y, P), half, y, P, 1.0)
                    - odd_tail_hilbert(y, P))), 1e-15)

print("\n" + ("ALL PASS" if ok else "FAILURES PRESENT"))
sys.exit(0 if ok else 1)
