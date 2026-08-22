"""Verify the (1-iz)^{-q} algebra against mpmath: H, d/dy, and the antiderivatives.

Part of the PAUSED ccf/ subproject (technique validation only; see ../../README.md).
Kept because the algebra is exact and reusable, and because it is the piece that
removed the divergent quadrature for U.
"""

import os
import sys

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tail_basis import Basis, C, S  # noqa: E402

mp.mp.dps = 40
ok = True


def check(name, err, tol):
    global ok
    good = bool(err < tol)
    ok = ok and good
    print(f"   [{'PASS' if good else 'FAIL'}] {name:50s} {err:.3e}  (tol {tol:.0e})")


def mp_S(q, x):
    return (1 + x**2) ** (-mp.mpf(q) / 2) * mp.sin(mp.mpf(q) * mp.atan(x))


def mp_C(q, x):
    return (1 + x**2) ** (-mp.mpf(q) / 2) * mp.cos(mp.mpf(q) * mp.atan(x))


def pv_hilbert(f, x):
    return float(mp.quad(lambda t: (f(x - t) - f(x + t)) / t,
                         [0, 1, 10, 100, mp.inf]) / mp.pi)


QS = [0.4585520188566, 0.9171040377132, 1.4585520188566, 2.4585520188566]
PTS = [0.25, 1.0, 3.7, 40.0]
h = 1e-6

print("1. H[C_q] = S_q  and  H[S_q] = -C_q, vs principal-value quadrature")
for q in QS:
    e1 = max(abs(pv_hilbert(lambda t: mp_C(q, t), x) - S(q, np.array([x]))[0]) for x in PTS)
    e2 = max(abs(pv_hilbert(lambda t: mp_S(q, t), x) + C(q, np.array([x]))[0]) for x in PTS)
    check(f"q={q:.6f}  H[C_q]=S_q, H[S_q]=-C_q", max(e1, e2), 1e-10)

print("\n2. antiderivatives, vs mpmath quadrature")
for q in QS:
    e1 = max(abs(float(mp.quad(lambda t: mp_C(q, t), [0, x]))
                 - S(q - 1, np.array([x]))[0] / (q - 1)) for x in PTS)
    e2 = max(abs(float(mp.quad(lambda t: mp_S(q, t), [0, x]))
                 - (1 - C(q - 1, np.array([x]))[0]) / (q - 1)) for x in PTS)
    check(f"q={q:.6f}  int C_q, int S_q", max(e1, e2), 1e-11)

print("\n3. derivatives, vs central differences")
y = np.array(PTS)
for q in QS:
    d_s = (S(q, y + h) - S(q, y - h)) / (2 * h)
    d_c = (C(q, y + h) - C(q, y - h)) / (2 * h)
    check(f"q={q:.6f}  dS=qC_(q+1), dC=-qS_(q+1)",
          max(np.max(np.abs(d_s - q * C(q + 1, y))),
              np.max(np.abs(d_c + q * S(q + 1, y)))), 1e-8)

print("\n4. Basis wiring, including the dilation scale")
a = np.array([0.7, -0.4, 0.25, 0.1])
for s in (1.0, 2.5):
    b = Basis(QS, s=s)
    y = np.array([0.3, 2.0, 15.0])
    om = a @ b.omega(y)
    e_h = np.max(np.abs(a @ b.h_omega(y)
                        + sum(ak * C(qk, y / s) for ak, qk in zip(a, QS))))
    e_u = np.max(np.abs(a @ b.big_u(y)
                        + sum(ak * s * S(qk - 1, y / s) / (qk - 1)
                              for ak, qk in zip(a, QS))))
    check(f"s={s}  h_omega / big_u consistency", max(e_h, e_u), 1e-15)
    du = (a @ b.big_u(y + h) - a @ b.big_u(y - h)) / (2 * h)
    check(f"s={s}  d/dy U = H(Omega)", np.max(np.abs(du - a @ b.h_omega(y))), 1e-7)
    df = (a @ b.big_f(y + h) - a @ b.big_f(y - h)) / (2 * h)
    check(f"s={s}  d/dy F = Omega", np.max(np.abs(df - om)), 1e-7)
    check(f"s={s}  F(0) = U(0) = 0",
          max(abs((a @ b.big_f(np.array([0.0])))[0]),
              abs((a @ b.big_u(np.array([0.0])))[0])), 1e-15)
    check(f"s={s}  Omega'(0) vector matches pointwise",
          abs(a @ b.omega_prime_at_zero() - (a @ b.omega_prime(np.array([0.0])))[0]), 1e-14)

print("\n5. parity: S_q odd, C_q even")
y = np.array(PTS)
for q in QS:
    check(f"q={q:.6f}  S odd, C even",
          max(np.max(np.abs(S(q, -y) + S(q, y))), np.max(np.abs(C(q, -y) - C(q, y)))), 1e-15)

print("\n" + ("ALL PASS" if ok else "FAILURES PRESENT"))
sys.exit(0 if ok else 1)
