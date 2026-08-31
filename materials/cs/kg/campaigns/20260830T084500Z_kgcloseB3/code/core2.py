"""core2.py — shared certified Arb machinery for the Gate-B close-out campaign
(agent KgGateBClose, 2026-08-30).

 provenance: `gauss_legendre_32`, `_legendre_pd`, `Prec` are byte-equivalent
 re-implementations of the predecessor campaign's core.py (campaign
 20260830T013720Z_e24aa159_afa514.../code/core.py), audited line-by-line:
  - GL nodes bracketed by sign-change scan of P_32 (20,000 subintervals),
    Newton-polished to ball plateau, degree-63 exactness self-check with
    |err| < 1e-6 (and odd moments ~ 0);
  - all arithmetic inside `Prec(bits)` runs with python-flint arb ball
    arithmetic = Arb outward rounding, PREC = 256 bits.
"""
from __future__ import annotations

import math
from typing import Callable, List, Optional, Sequence, Tuple

from flint import fmpz, fmpq, arb, acb, ctx

PREC: int = 256


class Prec:
    def __init__(self, bits_: int = PREC):
        self.bits_ = bits_

    def __enter__(self):
        self.old = ctx.prec
        ctx.prec = self.bits_
        return self

    def __exit__(self, *a):
        ctx.prec = self.old
        return False


# ---------------------------------------------------------------------------
# Certified 32-point Gauss-Legendre on [-1, 1] (audited port of predecessor)
# ---------------------------------------------------------------------------

_gl32: Optional[Tuple[Sequence[arb], Sequence[arb]]] = None


def _legendre_pd(m: int, x: arb) -> Tuple[arb, arb]:
    if m == 0:
        return arb(1), arb(0)
    if m == 1:
        return x, arb(1)
    p0, p1 = arb(1), x
    dp0, dp1 = arb(0), arb(1)
    for k in range(2, m + 1):
        kk = arb(k)
        p2 = ((2 * kk - 1) * x * p1 - (kk - 1) * p0) / kk
        dp2 = ((2 * kk - 1) * (p1 + x * dp1) - (kk - 1) * dp0) / kk
        p0, p1 = p1, p2
        dp0, dp1 = dp1, dp2
    return p1, dp1


def gauss_legendre_32() -> Tuple[Sequence[arb], Sequence[arb]]:
    global _gl32
    if _gl32 is not None:
        return _gl32
    n = 32
    nodes: List[arb] = []
    weights: List[arb] = []
    with Prec(500):
        roots: List[Tuple[arb, arb]] = []
        prev_pt = arb(-1)
        p0, _ = _legendre_pd(n, prev_pt)
        prev_sign = 1 if p0 > 0 else -1
        n_scan = 20000
        for i in range(1, n_scan + 1):
            b = arb(-1) + fmpq(i, n_scan)
            p, _ = _legendre_pd(n, b)
            s = 1 if p > 0 else -1
            if s != prev_sign:
                roots.append((prev_pt, b))
            prev_pt, prev_sign = b, s
        for (u, v) in roots:
            xa = (u + v) / arb(2)
            prev_rad = None
            for _it in range(120):
                p, dp = _legendre_pd(n, xa)
                if prev_rad is not None and p.rad() >= prev_rad:
                    break
                prev_rad = p.rad()
                xa = xa - p / dp
            _, dp_end = _legendre_pd(n, xa)
            nodes.append(xa)
            weights.append(arb(2) / ((arb(1) - xa * xa) * dp_end * dp_end))
        full_nodes: List[arb] = []
        full_weights: List[arb] = []
        for k in range(n // 2):
            full_nodes.append(-nodes[n // 2 - 1 - k])
            full_weights.append(weights[n // 2 - 1 - k])
        for k in range(n // 2):
            full_nodes.append(nodes[k])
            full_weights.append(weights[k])
        for m in range(2 * n):
            lhs = arb(0)
            for xi, wi in zip(full_nodes, full_weights):
                lhs += wi * xi ** m
            if m % 2 == 0:
                rhs = (arb(1) - arb(-1) ** (m + 1)) / arb(m + 1)
                err = lhs - rhs
                if not abs(err) < arb(10) ** -6:
                    raise RuntimeError(f"GL32 exactness failed at m={m}: {err.str(5)}")
            else:
                if not abs(lhs) < arb(10) ** -120:
                    raise RuntimeError(f"GL32 odd moment failed at m={m}: {lhs.str(5)}")
        _gl32 = (full_nodes, full_weights)
    return _gl32


# ---------------------------------------------------------------------------
# psi_j(s), K_e, K_o  (orthonormal probabilist Hermite, paper Section 8)
# ---------------------------------------------------------------------------

def psi(n: int, s: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        if n == 0:
            return arb(1)
        if n == 1:
            return s
        if n == 2:
            return (s * s - arb(1)) / arb(2).sqrt()
        if n == 3:
            return (s ** 3 - arb(3) * s) / arb(6).sqrt()
        raise ValueError(n)


def K_e(s: arb, bits: int = PREC) -> arb:
    """psi0^2 + psi2^2 = 3/2 - s^2 + s^4/2  (paper line 1980)."""
    with Prec(bits):
        return arb(3)/2 - s * s + (s ** 4) / 2


def K_o(s: arb, bits: int = PREC) -> arb:
    """psi1^2 + psi3^2 = (5/2)s^2 - s^4 + s^6/6  (paper line 1976)."""
    with Prec(bits):
        return (arb(5)/2) * s * s - s ** 4 + (s ** 6) / 6


def phi_density(s: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return (-s * s / 2).exp() / (arb(2) * arb.pi()).sqrt()


def nu(bits: int = PREC) -> arb:
    with Prec(bits):
        return (arb(2) / arb.pi()).sqrt()


def gauss_cdf(z: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return (arb(1) + (z / arb(2).sqrt()).erf()) / arb(2)


def tail_account(S_MAX, bits: int = PREC) -> arb:
    """Certified enclosure of int_{S_MAX}^inf (1 + (s^2+1)/sqrt(2) + s^3) phi(s) ds.

    Majorant per paper lines 1982-1984.  With S_MAX = 8 the paper certifies
    the value <= 3.7e-13 (true value ~ 3.64e-13).  Computed here in exact-ish
    ball arithmetic: split into termwise closed forms with certified integrals:
      int_a^inf phi = erfc(a/sqrt2)/2
      int_a^inf s^2 phi = a phi(a) + erfc(a/sqrt2)/2
      int_a^inf s^3 phi = (a^2+2) phi(a)
    each evaluated as balls (erfc carries its own rigorous radius).
    """
    with Prec(bits):
        a = arb(S_MAX)
        p = phi_density(a)
        e2 = (a / arb(2).sqrt()).erfc() / 2  # int_a^inf phi
        s2 = a * p + e2
        s3 = (a * a + 2) * p
        sq2 = arb(2).sqrt()
        return e2 + (s2 * 1 + e2) / sq2 + s3  # 1*e2 + (s2+1*e2)/sqrt2 + s3


def d_of_c(c: arb, bits: int = PREC) -> arb:
    """d(c) = (3 nu / 2) (sqrt(1+c^2) - c)   (paper line 1923)."""
    with Prec(bits):
        return (arb(3) * nu() / 2) * ((arb(1) + c * c).sqrt() - c)


sqrt_25 = (arb(5) / 2).sqrt()   # skip threshold used by gate_B2 (K decreasing for s >= sqrt(5/2))


def K(s: arb, bits: int = PREC) -> arb:
    """K(s) = sum_{j<=3} psi_j(s)^2 = 3/2 + (3/2)s^2 - (1/2)s^4 + (1/6)s^6
    (paper line 2051)."""
    with Prec(bits):
        return arb(3)/2 + (arb(3)/2) * s * s - (s ** 4) / 2 + (s ** 6) / 6
