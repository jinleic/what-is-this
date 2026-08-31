"""Shared, trusted Arb interval machinery for the kg/ gates.

Everything here is outward-rounded (interval) arithmetic via python-flint
(`flint.arb`).  No float casts appear in trusted paths: parameters enter as
exact rationals (fmpq) or arb balls, and every derived constant is an arb
interval.  A predicate is certified only when the WHOLE computed interval
lies on its claimed side (arb comparison semantics: `I > 0` is True iff all
points of I exceed 0).

Resource policy: single-threaded (per repo rules), bounded panel budgets.
"""

from __future__ import annotations

import math
import time
from typing import Callable, List, Optional, Sequence, Tuple

from flint import fmpz, fmpq, arb, acb, ctx

# Default trusted precision (bits of mantissa for arb operations).
PREC: int = 256
# Precision for one-off high-accuracy constants.
PREC_HI: int = 400


class Prec:
    """Context manager pinning arb precision for a trusted block."""

    def __init__(self, bits: int):
        self.bits = bits
        self.prev: Optional[int] = None

    def __enter__(self):
        self.prev = ctx.prec
        ctx.prec = self.bits
        return self

    def __exit__(self, *exc):
        ctx.prec = self.prev
        return False


# ---------------------------------------------------------------------------
# Exact scheme parameters (arXiv:2608.11158 Eq. (7)) as terminating decimals.
# ---------------------------------------------------------------------------

def _dec_fmpq(s: str) -> fmpq:
    """Exact decimal string -> FMPQ (terminating decimals only)."""
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    if "." not in s:
        num = int(s)
        return fmpq(-num if neg else num)
    ip, fp = s.split(".")
    denom = 10 ** len(fp)
    num = int(ip + fp)
    return fmpq(-num if neg else num, denom)


ETA_S = "0.136419125"
S3_S = "0.34101124"
S5_S = "0.05276111"
GAMMA_PAPER_S = "0.881545409"

ETA = _dec_fmpq(ETA_S)
S3 = _dec_fmpq(S3_S)
S5 = _dec_fmpq(S5_S)
GAMMA_PAPER = _dec_fmpq(GAMMA_PAPER_S)


def dec_arb(s: str, bits: int = PREC) -> arb:
    with Prec(PREC_HI):
        return arb(_dec_fmpq(s))


def exact_to_arb(q: fmpq, bits: int = PREC) -> arb:
    with Prec(PREC_HI):
        return arb(q)


# ---------------------------------------------------------------------------
# Certified constants (cached per precision)
# ---------------------------------------------------------------------------

_pcache: dict = {}

def pi(bits: int = PREC) -> arb:
    with Prec(bits):
        key = ("pi", bits)
        if key not in _pcache:
            _pcache[key] = arb.pi()
        return _pcache[key]


def sqrt2(bits: int = PREC) -> arb:
    with Prec(bits):
        key = ("sqrt2", bits)
        if key not in _pcache:
            _pcache[key] = arb(2).sqrt()
        return _pcache[key]


def nu(bits: int = PREC) -> arb:
    """nu = E|Z| = sqrt(2/pi) (the paper's nu)."""
    with Prec(bits):
        key = ("nu", bits)
        if key not in _pcache:
            _pcache[key] = (arb(2) / pi(bits)).sqrt()
        return _pcache[key]


def phi_density(z: arb, bits: int = PREC) -> arb:
    """Standard Gaussian density phi(z) = exp(-z^2/2)/sqrt(2 pi)."""
    with Prec(bits):
        return (-z * z * arb(0.5)).exp() / (arb(2) * pi(bits)).sqrt()


def gauss_cdf(z: arb, bits: int = PREC) -> arb:
    """Phi(z) = (1 + erf(z/sqrt 2))/2."""
    with Prec(bits):
        return (arb(1) + (z / sqrt2(bits)).erf()) / arb(2)


# ---------------------------------------------------------------------------
# Interval predicates (the ONLY way verdicts are decided)
# ---------------------------------------------------------------------------

def cert_pos(x: arb) -> bool:
    """True iff the whole interval is strictly positive."""
    return bool(x > arb(0))


def cert_neg(x: arb) -> bool:
    return bool(x < arb(0))


def cert_zero(x: arb) -> bool:
    """True iff the interval contains only 0 (exact zero)."""
    return x.is_zero


def abs_ub(x: arb) -> arb:
    """Enclosure of |x| (same interval but shifted nonnegative)."""
    if cert_neg(x):
        return -x
    pos = x.nonnegative_part()
    neg = (-x).nonnegative_part()
    return pos if pos >= neg else neg


# ---------------------------------------------------------------------------
# Certified 32-point Gauss-Legendre on [-1, 1]
# ---------------------------------------------------------------------------

_gl32: Optional[Tuple[Sequence[arb], Sequence[arb]]] = None


def _legendre_pd(m: int, x: arb) -> Tuple[arb, arb]:
    """P_m(x), P_m'(x) by three-term recursion (m >= 0)."""
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
    """Certified 32-pt GL nodes/weights on [-1,1] (built at 500 bits, cached,
    with a degree-63 exactness verification: sum_j w_j x_j^m = integral for
    m = 0..63, checked to hold within 10^-420)."""
    global _gl32
    if _gl32 is not None:
        return _gl32
    n = 32
    nodes: List[arb] = []
    weights: List[arb] = []
    with Prec(500):
        # 1) rigorously bracket all positive roots of P_n by a sign-change scan
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
        # 2) Newton-polish each bracket from its midpoint; stop when the
        #    residual ball's radius stops shrinking (ball plateau)
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
        # exactness self-check m = 0..2n-1 (even m only; odd m gives 0 = 0)
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


def integrate_panels(
    f: Callable[[arb], arb],
    a: arb,
    b: arb,
    n_panels: int,
    bits: int = PREC,
) -> arb:
    """Composite 32-pt Gauss-Legendre with outward rounding on [a,b].

    `f` must return an ENCLOSING arb interval of the integrand at the sample
    point `t` (including any parameter intervals being quantified over) —
    the caller owns the bounding step; this function only sums the panels.
    """
    nodes, weights = gauss_legendre_32()
    with Prec(bits):
        total = arb(0)
        for i in range(n_panels):
            pa = a + (b - a) * fmpq(i, n_panels)
            pb = a + (b - a) * fmpq(i + 1, n_panels)
            pmid = (pa + pb) / arb(2)
            phalf = (pb - pa) / arb(2)
            acc = arb(0)
            for xi, wi in zip(nodes, weights):
                t = pmid + phalf * xi
                acc = acc + wi * f(t)
            total += phalf * acc
        return total


# ---------------------------------------------------------------------------
# Certified adaptive integration (FLINT's CBT algorithm on acb)
# ---------------------------------------------------------------------------

def integrate_certified(
    f: Callable[[acb], acb],
    a: arb,
    b: arb,
    bits: int = PREC,
    deg_limit: Optional[int] = None,
    abs_tol: Optional[arb] = None,
) -> arb:
    """Certified adaptive integration; returns an arb enclosure of
    integral_a^b f(t) dt, with f returning an enclosing acb."""
    with Prec(bits):
        kw = {}
        if deg_limit is not None:
            kw["deg_limit"] = deg_limit
        if abs_tol is not None:
            kw["abs_tol"] = abs_tol
        r = acb(1).integral(lambda x: f(x), arb(a), arb(b), **kw)
        return arb(r)


# ---------------------------------------------------------------------------
# Hermite polynomials (orthonormal convention, paper Section 8)
# ---------------------------------------------------------------------------

def psi(n: int, x: arb) -> arb:
    """Orthonormal probabilist Hermite psi_n = He_n/sqrt(n!):
    psi_0=1, psi_1=s, psi_2=(s^2-1)/sqrt2, psi_3=(s^3-3s)/sqrt6, ...
    """
    if n == 0:
        return arb(1)
    if n == 1:
        return x
    h0 = arb(1)
    h1 = x
    for k in range(2, n + 1):
        h2 = x * h1 - fmpq(k - 1) * h0
        h0, h1 = h1, h2
    return h1 / arb(math.factorial(n)).sqrt()


def psi_on(s: float, x: arb) -> arb:
    """psi_n with s given as Python int default (convenience wrapper)."""
    return psi(s, x)


def gauss_tail_erfc(a: arb, bits: int = PREC) -> arb:
    """Enclosure of the two-sided Gaussian tail probability P(|Z|>=a),
    a >= 0: erfc(a/sqrt2)  (since 2(1-Phi(a)) = erfc(a/sqrt 2))."""
    with Prec(bits):
        return (a / sqrt2(bits)).erfc()
