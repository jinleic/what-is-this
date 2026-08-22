"""Shared rigorous Arb enclosures for the union-closed entropy certifiers.

Importing this module only defines constants and functions.  In particular, the
moderately expensive global ``rh`` cover is computed lazily by
:func:`get_rh_gmax`, and validation is confined to the ``__main__`` block.
"""

from flint import arb, ctx

ctx.prec = 160

ALPHA = arb("0.0356069")
ONE = arb(1)
TWO = arb(2)
HALF = ONE / TWO
LOG2 = TWO.log()


def hull(*balls):
    """Return the outward-rounded interval hull of one or more Arb balls."""
    if not balls:
        raise ValueError("hull needs at least one ball")
    out = arb(balls[0])
    for ball in balls[1:]:
        out = out.union(arb(ball))
    return out


def h_pt(x):
    """Binary entropy at an Arb point, with exact endpoint values."""
    x = arb(x)
    if x <= 0 or x >= 1:
        return arb(0)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def h_encl(u, v):
    """Exact range enclosure of binary entropy on ``[u, v]``."""
    u, v = arb(u), arb(v)
    e = hull(h_pt(u), h_pt(v))
    if v <= HALF or u >= HALF:
        return e
    return hull(e, ONE)


def amin(a, b):
    """Certified interval extension of ``min`` for Arb balls."""
    a, b = arb(a), arb(b)
    if a <= b:
        return a
    if b <= a:
        return b
    return hull(a, b)


def amax(a, b):
    """Certified interval extension of ``max`` for Arb balls."""
    a, b = arb(a), arb(b)
    if a >= b:
        return a
    if b >= a:
        return b
    return hull(a, b)


def sstar_pt(p, q):
    """The clipped median ``s*(p,q)`` at Arb points."""
    p, q = arb(p), arb(q)
    return amin(amax(HALF, amax(p, q)), amin(p + q, ONE))


def sstar_encl(pl, ph, ql, qh):
    """Range enclosure of the coordinatewise nondecreasing ``s*``."""
    return hull(sstar_pt(pl, ql), sstar_pt(ph, qh))


def rh_encl(u, v, gmax=None):
    """Enclose ``max(0, 2(1-p)h(p)-h((1-p)^2))`` on ``[u,v]``."""
    u, v = arb(u), arb(v)
    p = hull(u, v)
    hp = h_encl(u, v)
    x_lo = (ONE - v) * (ONE - v)
    x_hi = (ONE - u) * (ONE - u)
    h2 = h_encl(x_lo, x_hi)
    raw = TWO * (ONE - p) * hp - h2
    lo = amax(raw.lower(), arb(0))
    hi = amax(raw.upper(), arb(0))
    if gmax is not None:
        hi = amin(hi, arb(gmax))
    return hull(lo, hi)


def sqrt_rh_enc(u, v, gmax=None):
    """NaN-free two-sided enclosure of ``sqrt(rh)`` on ``[u,v]``.

    A square root is taken only from a certified positive lower endpoint and a
    certified nonnegative upper endpoint.  Thus narrow positive boxes retain a
    useful lower bound, while boxes reaching zero conservatively start at zero.
    """
    e = rh_encl(u, v, gmax)
    lo = e.lower()
    sqrt_lo = lo.sqrt() if lo > 0 else arb(0)
    hi = amax(e.upper(), arb(0)).upper()
    return hull(sqrt_lo, hi.sqrt())


def _rh_global(n=20000):
    """Prove an upper bound for global ``rh`` by covering ``[0,1]``."""
    if n <= 0:
        raise ValueError("n must be positive")
    best = arb(0)
    for i in range(n):
        e = rh_encl(arb(i) / n, arb(i + 1) / n)
        upper = e.upper()
        if not (upper <= best):
            best = amax(best, upper)
    return best.upper()


RH_GMAX = None
_RH_GMAX_N = None


def get_rh_gmax(n=20000):
    """Return the lazily computed, process-local global ``rh`` upper bound."""
    global RH_GMAX, _RH_GMAX_N
    if RH_GMAX is None or _RH_GMAX_N != n:
        RH_GMAX = _rh_global(n)
        _RH_GMAX_N = n
    return RH_GMAX


if __name__ == "__main__":
    assert not (arb(0.0, 1.0) >= 0)
    assert not (arb(0.5, 0.5) >= 0)
    assert arb(0.5, 0.25) >= 0
    bound = get_rh_gmax()
    assert bound >= arb("0.2342294") and bound <= arb("0.2350")
    endpoint = sqrt_rh_enc(arb("0.95"), arb(1), bound)
    assert endpoint.is_finite()
    assert endpoint.lower() <= 0 and endpoint.upper() <= arb("0.17")
    print("arbcore validation PASS: certified comparisons and finite endpoint enclosure")
    print("certified 20,000-cell global rh upper bound: %.12f" % float(bound.upper()))
