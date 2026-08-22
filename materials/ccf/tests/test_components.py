"""Component tests for the CCF profile machinery. Each piece against an exact answer.

NOTE: the Hilbert-transform sections here are superseded by test_hilbert.py, which
covers the tail-corrected operator. What remains uniquely useful in this file is the
quadrature and closed-form antiderivative validation for ccf_profile.py.
"""

import os
import sys

import numpy as np
import mpmath as mp
from scipy.special import hyp2f1

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hilbert import hilbert_line, make_grid
from ccf_profile import cumulative_theta, int_even_model, int_odd_model
mp.mp.dps = 40
LAM = 1.1807776628998
P = 1.0 / (1.0 + LAM)

print(f"lambda = {LAM},  p = 1/(1+lambda) = {P:.13f}\n")

# ---------------------------------------------------------------- 1. cumulative_theta
# g = 1/(1+y^2)^2  ->  int_0^y g = y/(2(1+y^2)) + arctan(y)/2   (fast decay, clean test)
print("1. cumulative_theta on a fast-decaying integrand")
for n in (256, 1024, 4096, 16384):
    theta, y, half_t = make_grid(n)
    g = 1.0 / (1.0 + y**2) ** 2
    exact = y / (2.0 * (1.0 + y**2)) + np.arctan(y) / 2.0
    got = cumulative_theta(g, y, n)
    print(f"   N={n:6d}  max err = {np.max(np.abs(got - exact)):.3e}")

# ---------------------------------------------------------------- 2. closed forms
print("\n2. closed-form antiderivatives vs mpmath")
for yv in (0.5, 2.0, 50.0, 1000.0, 1e5):
    ie = int_even_model(np.array([yv]), P)[0]
    ie_ref = float(mp.quad(lambda s: (1 + s**2) ** (-P / 2), [0, yv]))
    io = int_odd_model(np.array([yv]), P)[0]
    io_ref = float(mp.quad(lambda s: s * (1 + s**2) ** (-(P + 1) / 2), [0, yv]))
    print(f"   y={yv:9.1f}  even rel err {abs(ie/ie_ref-1):.3e}   "
          f"odd rel err {abs(io/io_ref-1):.3e}")

# ---------------------------------------------------------------- 3. full U pipeline
# Pick Omega = y(1+y^2)^{-(p+1)/2} exactly (i.e. V == 1). Compute H Omega numerically,
# then U = int_0^y H Omega two ways: our split routine vs direct mpmath double integral.
print("\n3. U = int_0^y H(Omega) pipeline, Omega = y(1+y^2)^{-(p+1)/2}")


def u_pipeline(n):
    theta, y, half_t = make_grid(n)
    om = y * (1.0 + y**2) ** (-(P + 1.0) / 2.0)
    hom = hilbert_line(om, half_t)
    c_u = hom[-1] * (1.0 + y[-1] ** 2) ** (P / 2.0)
    rem = hom - c_u * (1.0 + y**2) ** (-P / 2.0)
    return y, hom, c_u, c_u * int_even_model(y, P) + cumulative_theta(rem, y, n)


def hom_exact(x):
    """H Omega by direct principal-value quadrature at a single point."""
    f = lambda s: s * (1 + s**2) ** (-(P + 1) / 2)
    g = lambda s: (f(x - s) - f(x + s)) / s  # PV-symmetrised
    return float(mp.quad(g, [0, 1, 10, 100, mp.inf]) / mp.pi)


for n in (1024, 4096, 16384):
    y, hom, c_u, big_u = u_pipeline(n)
    idx = [np.argmin(np.abs(y - t)) for t in (0.5, 2.0, 20.0)]
    errs = [abs(hom[i] - hom_exact(y[i])) for i in idx]
    print(f"   N={n:6d}  c_u={c_u: .6f}  H-Omega errs at y~0.5,2,20: "
          + ", ".join(f"{e:.2e}" for e in errs))

# convergence of U itself between resolutions, sampled at fixed y
print("\n   U self-convergence at fixed y (Richardson between N):")
tab = {}
for n in (1024, 4096, 16384, 65536):
    y, hom, c_u, big_u = u_pipeline(n)
    tab[n] = [big_u[np.argmin(np.abs(y - t))] for t in (0.5, 2.0, 20.0)]
ns = sorted(tab)
for a, b in zip(ns, ns[1:]):
    d = [abs(x - z) for x, z in zip(tab[a], tab[b])]
    print(f"   N {a:6d}->{b:6d}   dU at y~0.5,2,20: " + ", ".join(f"{e:.2e}" for e in d))

# ---------------------------------------------------------------- 4. decay of H Omega
print("\n4. is H(Omega) ~ c|y|^{-p} in the far field? (checks the subtraction model)")
n = 16384
theta, y, half_t = make_grid(n)
om = y * (1.0 + y**2) ** (-(P + 1.0) / 2.0)
hom = hilbert_line(om, half_t)
for t in (10.0, 100.0, 1000.0, 5000.0):
    i = np.argmin(np.abs(y - t))
    print(f"   y={y[i]:10.2f}   H_Omega*|y|^p = {hom[i]*abs(y[i])**P: .8f}")
