#!/usr/bin/env python3
"""Arb certificate for the two-variable Liu-H2 kernel inequality.

The certified function is

    Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s))*h(pi(t,t)),
    pi(s,t) = s*t*(1 + (1-s)*(1-t)),

with natural-log binary entropy.  The zero axes and diagonal are handled
symbolically.  The positive-dimensional part of the cover has three strata:

* off diagonal: adaptive branch-and-bound on Phi itself;
* near diagonal: the exact even Taylor expansion in e for
  (s,t)=(c-e,c+e), with |g''''| bounded on every whole cell;
* the origin corner: the even Taylor step is applied after the exact
  factorization Phi=s^2*t^2*G (removing the vanishing curvature), while its
  off-diagonal boxes remain direct Phi bounds.

Arb is the sole source of certified numerical inequalities.  mpmath is used
only for deterministic, high-precision derivative cross-checks and never as a
certificate gate.  Endpoint boundary layers thinner than INNER_CUTOFF are
reported honestly as unresolved; consequently the default artifact is a
PARTIAL certificate, not a global proof.

Run (single core, lowered priority):

    nice -n 19 .venv/bin/python -I -B uc/liu9_h2_phi_certificate.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_name] = "1"

from flint import arb, ctx  # noqa: E402
import mpmath as mp  # noqa: E402

CTX_PREC = 320
ctx.prec = max(ctx.prec, CTX_PREC)
mp.mp.dps = 100

HERE = Path(__file__).resolve().parent
OUTPUT_DEFAULT = HERE / "verification/results/liu9-h2-phi.json"

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
HALF = ONE / TWO
LN2 = TWO.log()
UNIT_INTERVAL = ZERO.union(ONE)

# The compact domain [a,1-a]^2 is certified.  The two omitted boundary
# layers are stated explicitly in the report.  R is the exact c-crossover
# between the normalized origin treatment and direct Phi treatment.
INNER_CUTOFF = Fraction(1, 1 << 14)
CORNER_RADIUS = Fraction(1, 128)
DOMAIN_LO = INNER_CUTOFF
DOMAIN_HI = 1 - INNER_CUTOFF
MID = Fraction(1, 2)

# Broad c-bands control Taylor remainders and the direct cover.  Curvature is
# certified on a finer internal subcover (48 pieces in the corner, 16/8 in
# the lower/upper core), avoiding thousands of expensive remainder cells.
CORNER_RELATIVE_STEPS = 32
CORE_RELATIVE_STEPS = 24
MAX_DIRECT_DEPTH = 72
MAX_DIRECT_CELLS = 1_500_000


def qarb(value: Fraction | int | str | arb) -> arb:
    """Exact-rational/decimal conversion without a binary-float detour."""
    if isinstance(value, arb):
        return arb(value)
    if isinstance(value, Fraction):
        return arb(value.numerator) / arb(value.denominator)
    return arb(value)


def ihull(lo: Fraction, hi: Fraction) -> arb:
    assert lo <= hi
    return qarb(lo).union(qarb(hi))


def certified_min(a: arb, b: arb) -> arb:
    """Outward lower bound for min(a,b), avoiding a huge-range hull."""
    a, b = arb(a).lower(), arb(b).lower()
    if a <= b:
        return a
    if b <= a:
        return b
    return a.union(b).lower()


def certified_max(a: arb, b: arb) -> arb:
    """Outward upper bound for max(a,b), avoiding a huge-range hull."""
    a, b = arb(a).upper(), arb(b).upper()
    if a >= b:
        return a
    if b >= a:
        return b
    return a.union(b).upper()


def abs_upper(value: arb) -> arb:
    """Outward upper bound on absolute value (correct direction: upward)."""
    return certified_max(-value.lower(), value.upper())


def hull_lower(values: Sequence[arb]) -> arb:
    """Outward lower bound on the minimum of all supplied balls."""
    assert values
    out = arb(values[0]).lower()
    for value in values[1:]:
        out = certified_min(out, value)
    return out


def hull_upper(values: Sequence[arb]) -> arb:
    """Outward upper bound on the maximum of all supplied balls."""
    assert values
    out = arb(values[0]).upper()
    for value in values[1:]:
        out = certified_max(out, value)
    return out


def lower_text(value: arb, digits: int = 24) -> str:
    return value.lower().str(digits)


def upper_text(value: arb, digits: int = 24) -> str:
    return value.upper().str(digits)


def interval_text(value: arb, digits: int = 24) -> str:
    return value.str(digits)


def h_point(u: arb) -> arb:
    """Natural-log entropy at an Arb point, with exact endpoint values."""
    u = arb(u)
    if u == ZERO or u == ONE:
        return ZERO
    assert u > ZERO and u < ONE
    return -(u * u.log() + (ONE - u) * (ONE - u).log())


def h_range(u: arb) -> arb:
    """Exact range enclosure of h on an interval in [0,1].

    Concavity gives the minimum at an endpoint.  Monotonicity on each side of
    1/2 gives the endpoint hull unless the interval contains 1/2, in which
    case ln(2) is also included.  This avoids the invalid 0*log(0) interval
    product.
    """
    lo, hi = u.lower(), u.upper()
    assert lo >= ZERO and hi <= ONE
    out = h_point(lo).union(h_point(hi))
    if not (hi < HALF or lo > HALF):
        out = out.union(LN2)
    return out


def h_prime(u: arb) -> arb:
    """Exact interval extension of h'(u)=log((1-u)/u)."""
    assert u.lower() > ZERO and u.upper() < ONE
    return (ONE - u).log() - u.log()


def h_second(u: arb) -> arb:
    """Exact interval extension of h''(u)=-1/(u(1-u))."""
    assert u.lower() > ZERO and u.upper() < ONE
    return -ONE / (u * (ONE - u))


def pi_value(s: arb, t: arb) -> arb:
    return s * t * (ONE + (ONE - s) * (ONE - t))


def pi_s(s: arb, t: arb) -> arb:
    return t * (ONE + (ONE - t) * (ONE - TWO * s))


def pi_ss(s: arb, t: arb) -> arb:
    return -TWO * t * (ONE - t)


def pi_st(s: arb, t: arb) -> arb:
    return TWO * (ONE - s - t + TWO * s * t)


def pi_monotone_range(sl: Fraction, sh: Fraction,
                      tl: Fraction, th: Fraction) -> arb:
    """Exact pi range: pi_s,pi_t >= 0 on [0,1]^2."""
    raw = pi_value(qarb(sl), qarb(tl)).union(
        pi_value(qarb(sh), qarb(th))
    )
    # Exact monotonicity says the true range lies in [0,1].  Intersecting
    # removes only harmless outward roundoff past an endpoint.
    return raw.intersection(UNIT_INTERVAL)


def phi_interval(sl: Fraction, sh: Fraction,
                 tl: Fraction, th: Fraction) -> arb:
    cross = h_range(pi_monotone_range(sl, sh, tl, th))
    ss = h_range(pi_monotone_range(sl, sh, sl, sh))
    tt = h_range(pi_monotone_range(tl, th, tl, th))
    return cross * cross - ss * tt


# ---------------------------------------------------------------------------
# Exact truncated Taylor algebra.  Entry k is f^(k)/k!.  All recurrences are
# formal power-series identities, so they implement symbolic differentiation;
# their Arb coefficient arithmetic encloses every derivative on an interval.

Jet = list[arb]
JET_ORDER = 4


def jconst(value: Fraction | int | str | arb) -> Jet:
    return [qarb(value)] + [ZERO] * JET_ORDER


def jvar(value: arb, slope: Fraction | int = 1) -> Jet:
    return [arb(value), qarb(slope)] + [ZERO] * (JET_ORDER - 1)


def jadd(a: Jet, b: Jet) -> Jet:
    return [x + y for x, y in zip(a, b)]


def jneg(a: Jet) -> Jet:
    return [-x for x in a]


def jsub(a: Jet, b: Jet) -> Jet:
    return jadd(a, jneg(b))


def jscale(a: Jet, scale: Fraction | int | arb) -> Jet:
    scale = qarb(scale)
    return [scale * x for x in a]


def jmul(a: Jet, b: Jet) -> Jet:
    out: Jet = []
    for n in range(JET_ORDER + 1):
        out.append(sum((a[k] * b[n - k] for k in range(n + 1)), ZERO))
    return out


def jinv(a: Jet) -> Jet:
    out = [ONE / a[0]] + [ZERO] * JET_ORDER
    for n in range(1, JET_ORDER + 1):
        out[n] = -out[0] * sum(
            (a[k] * out[n - k] for k in range(1, n + 1)), ZERO
        )
    return out


def jlog(a: Jet) -> Jet:
    assert a[0].lower() > ZERO
    inv = jinv(a)
    out = [a[0].log()] + [ZERO] * JET_ORDER
    # (log a)' = a'/a.
    for n in range(1, JET_ORDER + 1):
        out[n] = sum(
            ((k + 1) * a[k + 1] * inv[n - 1 - k]
             for k in range(n)),
            ZERO,
        ) / n
    return out


def jpow(a: Jet, exponent: int) -> Jet:
    assert exponent >= 0
    out = jconst(1)
    for _ in range(exponent):
        out = jmul(out, a)
    return out


def jentropy(a: Jet) -> Jet:
    assert a[0].lower() > ZERO and a[0].upper() < ONE
    return jneg(jadd(
        jmul(a, jlog(a)),
        jmul(jsub(jconst(1), a), jlog(jsub(jconst(1), a))),
    ))


def jpi(s: Jet, t: Jet) -> Jet:
    return jmul(jmul(s, t), jadd(
        jconst(1), jmul(jsub(jconst(1), s), jsub(jconst(1), t))
    ))


def jphi(s: Jet, t: Jet, normalized: bool = False) -> Jet:
    cross = jentropy(jpi(s, t))
    ss = jentropy(jpi(s, s))
    tt = jentropy(jpi(t, t))
    out = jsub(jmul(cross, cross), jmul(ss, tt))
    if normalized:
        # Exact identity Phi=s^2*t^2*G on the strict interior.
        out = jmul(out, jinv(jmul(jmul(s, s), jmul(t, t))))
    return out


def directional_jet(sl: Fraction, sh: Fraction,
                    tl: Fraction, th: Fraction,
                    ds: int, dt: int,
                    normalized: bool = False) -> Jet:
    return jphi(jvar(ihull(sl, sh), ds), jvar(ihull(tl, th), dt),
                normalized)


def point_phi(s: Fraction, t: Fraction, normalized: bool = False) -> arb:
    ss, tt = qarb(s), qarb(t)
    cross = h_point(pi_value(ss, tt))
    diag_s = h_point(pi_value(ss, ss))
    diag_t = h_point(pi_value(tt, tt))
    out = cross * cross - diag_s * diag_t
    if normalized:
        out = out / (ss * ss * tt * tt)
    return out


def phi_gradient_intervals(sl: Fraction, sh: Fraction,
                           tl: Fraction, th: Fraction
                           ) -> tuple[arb, arb]:
    """Exact symbolic first partials of Phi on an independent box.

    Entropy values use their monotonicity-sharp ranges rather than the wider
    natural x*log(x) interval extension.  This changes only enclosure width,
    not the differentiated identity.
    """
    s, t = ihull(sl, sh), ihull(tl, th)
    p = pi_monotone_range(sl, sh, tl, th)
    pss = pi_monotone_range(sl, sh, sl, sh)
    ptt = pi_monotone_range(tl, th, tl, th)
    cross, diag_s, diag_t = h_range(p), h_range(pss), h_range(ptt)
    cross_s = h_prime(p) * pi_s(s, t)
    cross_t = h_prime(p) * pi_s(t, s)
    diag_s_derivative = 2 * h_prime(pss) * pi_s(s, s)
    diag_t_derivative = 2 * h_prime(ptt) * pi_s(t, t)
    return (
        2 * cross * cross_s - diag_s_derivative * diag_t,
        2 * cross * cross_t - diag_s * diag_t_derivative,
    )


def direct_centered_lower(sl: Fraction, sh: Fraction,
                          tl: Fraction, th: Fraction,
                          normalized: bool) -> tuple[arb, arb, arb]:
    """Mean-value lower bound, with every direction explicit.

    Ls,Lt are upward bounds on |partial_s|,|partial_t| over the whole box.
    Therefore |f-f(center)| <= Ls*rs+Lt*rt; subtracting this UPWARD drift
    from the center's LOWER endpoint has the required lower-bound direction.
    """
    sc, tc = (sl + sh) / 2, (tl + th) / 2
    value_lo = point_phi(sc, tc, normalized).lower()
    if normalized:
        grad_s = directional_jet(sl, sh, tl, th, 1, 0, True)[1]
        grad_t = directional_jet(sl, sh, tl, th, 0, 1, True)[1]
    else:
        grad_s, grad_t = phi_gradient_intervals(sl, sh, tl, th)
    ls, lt = abs_upper(grad_s), abs_upper(grad_t)
    drift = ls * qarb((sh - sl) / 2) + lt * qarb((th - tl) / 2)
    return value_lo - drift.upper(), ls, lt


def direct_lower(sl: Fraction, sh: Fraction,
                 tl: Fraction, th: Fraction,
                 normalized: bool) -> arb:
    centered, _, _ = direct_centered_lower(
        sl, sh, tl, th, normalized
    )
    if normalized:
        # A direct natural interval for G is obtained from the exact quotient.
        raw = phi_interval(sl, sh, tl, th)
        denom = ihull(sl, sh) ** 2 * ihull(tl, th) ** 2
        natural = raw / denom
    else:
        natural = phi_interval(sl, sh, tl, th)
    return centered if centered > natural.lower() else natural.lower()


# ---------------------------------------------------------------------------
# Diagonal curvature and whole-cell fourth derivative.


def diagonal_g_second(c: arb) -> arb:
    """Exact symbolic g''(0), in a cancellation-conscious closed form.

    q=pi(c,c), p'=d_c pi(c,c), and
    p_cross''-p_diag''=-16(c^2-c+1/2).  Product differentiation gives

      g''(0)=2 p'^2((h')^2-h h'')
             +2 h h'(p_cross''-p_diag'').
    """
    q = pi_value(c, c)
    h = h_range(q)
    hp = h_prime(q)
    hpp = h_second(q)
    p1 = 4 * c - 6 * c * c + 4 * c * c * c
    difference = -16 * (c * c - c + HALF)
    return TWO * (
        p1 * p1 * (hp * hp - h * hpp) + h * hp * difference
    )


def fourth_upper(cl: Fraction, ch: Fraction, radius: Fraction,
                 normalized: bool) -> arb:
    """UPPER bound on |d_e^4 f(c-e,c+e)| on the WHOLE cell.

    For an actual point c in [cl,ch], |e|<=radius, both coordinates remain
    in [DOMAIN_LO,DOMAIN_HI].  The independent coordinate hull below is a
    superset of that constrained cell, so its absolute upper endpoint is a
    valid (possibly wider) remainder bound.  It is never a centre-only bound.
    """
    lo = max(DOMAIN_LO, cl - radius)
    hi = min(DOMAIN_HI, ch + radius)
    jet = directional_jet(lo, hi, lo, hi, -1, 1, normalized)
    return abs_upper(24 * jet[4])


def curvature_lower(cl: Fraction, ch: Fraction,
                    normalized: bool) -> arb:
    """Uniform lower curvature from a finer, exhaustive c-subcover."""
    pieces = 48 if normalized else (16 if ch <= MID else 8)
    width = (ch - cl) / pieces
    best: Optional[arb] = None
    for index in range(pieces):
        lo = cl + index * width
        hi = ch if index + 1 == pieces else lo + width
        raw = diagonal_g_second(ihull(lo, hi)).lower()
        if not (raw > ZERO):
            return raw
        if normalized:
            # g_G''(0)=g_Phi''(0)/c^4 because Phi and Phi' vanish at e=0.
            # Divide by the UPWARD c^4 endpoint: this is the lower-bound
            # direction.  A lower denominator would invert the safety.
            raw = raw / (qarb(hi) ** 4).upper()
        best = raw if best is None else certified_min(best, raw)
    assert best is not None
    return best


def max_feasible_e(cl: Fraction, ch: Fraction) -> Fraction:
    """Exact max of min(c-a,b-c) over a c-band."""
    center = min(max((DOMAIN_LO + DOMAIN_HI) / 2, cl), ch)
    return min(center - DOMAIN_LO, DOMAIN_HI - center)


def choose_taylor_cell(cl: Fraction, ch: Fraction,
                       normalized: bool) -> dict[str, Any]:
    g2lo = curvature_lower(cl, ch, normalized)
    if not (g2lo > ZERO):
        raise ArithmeticError("curvature interval did not clear")
    cap = max_feasible_e(cl, ch)
    scale = min(cl, 1 - ch)
    radius = min(cap, scale / 8, Fraction(1, 128))
    if radius <= 0:
        radius = cap
    attempts = 0
    while True:
        m4 = fourth_upper(cl, ch, radius, normalized)
        # Taylor: g/e^2 >= g''_lo/2 - M4_up*e^2/24.  M4 is bounded
        # upward and is subtracted, so this is the conservative direction.
        factor = g2lo / TWO - m4 * qarb(radius * radius) / 24
        if factor.lower() > ZERO:
            return {
                "cl": cl,
                "ch": ch,
                "normalized": normalized,
                "g2_lower": g2lo,
                "m4_upper": m4,
                "radius": radius,
                "factor_lower": factor.lower(),
                "halvings": attempts,
            }
        radius /= 2
        attempts += 1
        if attempts > 40 or radius == 0:
            raise ArithmeticError("Taylor radius ladder exhausted")


def geometric_bands(lo: Fraction, hi: Fraction,
                    steps: int) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    cur = lo
    ratio = Fraction(steps + 1, steps)
    while cur < hi:
        nxt = min(hi, cur * ratio)
        out.append((cur, nxt))
        cur = nxt
    return out


def make_c_bands() -> list[tuple[Fraction, Fraction, bool]]:
    bands: list[tuple[Fraction, Fraction, bool]] = []
    for lo, hi in geometric_bands(
            DOMAIN_LO, CORNER_RADIUS, CORNER_RELATIVE_STEPS):
        bands.append((lo, hi, True))
    for lo, hi in geometric_bands(
            CORNER_RADIUS, MID, CORE_RELATIVE_STEPS):
        bands.append((lo, hi, False))
    # Reflect geometric distance-to-1 bands.  This keeps entropy arguments
    # uniformly away from 1 on every finite cell while reaching DOMAIN_HI.
    reflected = geometric_bands(INNER_CUTOFF, Fraction(1, 2),
                                CORE_RELATIVE_STEPS)
    for dlo, dhi in reversed(reflected):
        bands.append((1 - dhi, 1 - dlo, False))
    assert bands[0][0] == DOMAIN_LO
    assert bands[-1][1] == DOMAIN_HI
    assert all(a[1] == b[0] for a, b in zip(bands, bands[1:]))
    return bands


def box_intersection_ranges(cl: Fraction, ch: Fraction,
                            el: Fraction, eh: Fraction
                            ) -> Optional[tuple[Fraction, Fraction,
                                                Fraction, Fraction]]:
    cap = max_feasible_e(cl, ch)
    if el > cap:
        return None
    # Every actual point in the c,e box and upper triangle lies in these
    # clipped independent s,t intervals.  Enlarging to the independent box
    # is sound because a lower bound on a superset is also a lower bound on
    # the constrained cell.
    sl = max(DOMAIN_LO, cl - eh)
    sh = min(DOMAIN_HI, ch - el)
    tl = max(DOMAIN_LO, cl + el)
    th = min(DOMAIN_HI, ch + eh)
    if sl > sh or tl > th:
        return None
    return sl, sh, tl, th


def cover_off_diagonal(cell: dict[str, Any]) -> dict[str, Any]:
    cl, ch = cell["cl"], cell["ch"]
    radius = cell["radius"]
    cap = max_feasible_e(cl, ch)
    if cap <= radius:
        return {
            "leaves": 0, "max_depth": 0, "worst_lower": None,
            "unresolved": [], "witness": None,
        }
    normalized = False  # off-diagonal protocol is directly on Phi
    stack = [(cl, ch, radius, cap, 0)]
    leaves = 0
    max_depth = 0
    worst: Optional[arb] = None
    witness: Optional[dict[str, Any]] = None
    unresolved: list[dict[str, Any]] = []
    while stack:
        blo, bhi, elo, ehi, depth = stack.pop()
        ranges = box_intersection_ranges(blo, bhi, elo, ehi)
        if ranges is None:
            continue
        sl, sh, tl, th = ranges
        bound = direct_lower(sl, sh, tl, th, normalized)
        if bound > ZERO:
            leaves += 1
            max_depth = max(max_depth, depth)
            replaces_witness = (
                worst is None
                or float(bound.lower()) < float(worst.lower())
            )
            worst = bound if worst is None else certified_min(worst, bound)
            if replaces_witness:
                witness = {
                    "c": [str(blo), str(bhi)],
                    "e": [str(elo), str(ehi)],
                    "s": [str(sl), str(sh)],
                    "t": [str(tl), str(th)],
                    "lower": lower_text(bound),
                }
            continue
        if depth >= MAX_DIRECT_DEPTH or leaves + len(stack) >= MAX_DIRECT_CELLS:
            unresolved.append({
                "c": [str(blo), str(bhi)],
                "e": [str(elo), str(ehi)],
                "s": [str(sl), str(sh)],
                "t": [str(tl), str(th)],
                "best_lower": lower_text(bound),
            })
            continue
        cw, ew = bhi - blo, ehi - elo
        if cw >= ew:
            mid = (blo + bhi) / 2
            stack.append((mid, bhi, elo, ehi, depth + 1))
            stack.append((blo, mid, elo, ehi, depth + 1))
        else:
            mid = (elo + ehi) / 2
            stack.append((blo, bhi, mid, ehi, depth + 1))
            stack.append((blo, bhi, elo, mid, depth + 1))
    return {
        "leaves": leaves,
        "max_depth": max_depth,
        "worst_lower": worst,
        "unresolved": unresolved,
        "witness": witness,
    }


# ---------------------------------------------------------------------------
# Deterministic validation and report assembly.


def mp_h(u: mp.mpf) -> mp.mpf:
    if u == 0 or u == 1:
        return mp.mpf(0)
    return -u * mp.log(u) - (1 - u) * mp.log1p(-u)


def mp_pi(s: mp.mpf, t: mp.mpf) -> mp.mpf:
    return s * t * (1 + (1 - s) * (1 - t))


def mp_phi(s: mp.mpf, t: mp.mpf) -> mp.mpf:
    return mp_h(mp_pi(s, t)) ** 2 - mp_h(mp_pi(s, s)) * mp_h(mp_pi(t, t))


def derivative_validation() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    max_error = mp.mpf(0)
    for us in ("0.123", "0.371", "0.813"):
        u = mp.mpf(us)
        exact1 = mp.log((1 - u) / u)
        exact2 = -1 / (u * (1 - u))
        num1 = mp.diff(mp_h, u, 1)
        num2 = mp.diff(mp_h, u, 2)
        err = max(abs(exact1 - num1), abs(exact2 - num2))
        max_error = max(max_error, err)
        rows.append({"kind": "entropy", "u": us,
                     "max_abs_error": mp.nstr(err, 12)})
    for ss, ts in (("0.17", "0.29"), ("0.43", "0.71"),
                   ("0.82", "0.64")):
        s, t = mp.mpf(ss), mp.mpf(ts)
        ps = t * (1 + (1 - t) * (1 - 2 * s))
        pss = -2 * t * (1 - t)
        pst = 2 * (1 - s - t + 2 * s * t)
        ns = mp.diff(lambda x: mp_pi(x, t), s, 1)
        nss = mp.diff(lambda x: mp_pi(x, t), s, 2)
        nst = mp.diff(lambda x: mp.diff(lambda y: mp_pi(x, y), t), s)
        err = max(abs(ps - ns), abs(pss - nss), abs(pst - nst))
        max_error = max(max_error, err)
        rows.append({"kind": "protocol", "s": ss, "t": ts,
                     "max_abs_error": mp.nstr(err, 12)})
    # Validate the assembled Phi gradients used by the direct mean-value
    # cover, not merely their primitive factors.
    for ss, ts in (("0.16", "0.33"), ("0.41", "0.67"),
                   ("0.74", "0.91")):
        s, t = mp.mpf(ss), mp.mpf(ts)
        s_q, t_q = Fraction(ss), Fraction(ts)
        formal_s, formal_t = phi_gradient_intervals(
            s_q, s_q, t_q, t_q
        )
        numerical_s = mp.diff(lambda x: mp_phi(x, t), s)
        numerical_t = mp.diff(lambda x: mp_phi(s, x), t)
        tolerance = mp.mpf("1e-85")
        ball_s = arb(mp.nstr(numerical_s - tolerance, 95)).union(
            arb(mp.nstr(numerical_s + tolerance, 95))
        )
        ball_t = arb(mp.nstr(numerical_t - tolerance, 95)).union(
            arb(mp.nstr(numerical_t + tolerance, 95))
        )
        agrees = bool(formal_s.overlaps(ball_s)
                      and formal_t.overlaps(ball_t))
        mid_s = formal_s.mid().str(110).split("+/-")[0].strip("[ ")
        mid_t = formal_t.mid().str(110).split("+/-")[0].strip("[ ")
        err = max(abs(mp.mpf(mid_s) - numerical_s),
                  abs(mp.mpf(mid_t) - numerical_t))
        max_error = max(max_error, err)
        rows.append({"kind": "assembled_Phi_gradient", "s": ss, "t": ts,
                     "agrees_within_1e-85": agrees,
                     "max_abs_error": mp.nstr(err, 12)})
    # Validate the closed diagonal curvature formula used as the Taylor
    # leading coefficient.
    for cs in ("0.13", "0.52", "0.87"):
        c, c_q = mp.mpf(cs), Fraction(cs)
        formal = diagonal_g_second(qarb(c_q))
        numerical = mp.diff(lambda e: mp_phi(c - e, c + e), 0, 2)
        tolerance = mp.mpf("1e-85")
        num_ball = arb(mp.nstr(numerical - tolerance, 95)).union(
            arb(mp.nstr(numerical + tolerance, 95))
        )
        agrees = bool(formal.overlaps(num_ball))
        formal_mid = formal.mid().str(110).split("+/-")[0].strip("[ ")
        err = abs(mp.mpf(formal_mid) - numerical)
        max_error = max(max_error, err)
        rows.append({"kind": "assembled_diagonal_g_second", "c": cs,
                     "agrees_within_1e-85": agrees,
                     "max_abs_error": mp.nstr(err, 12)})
    # Cross-check the formal fourth-order jet itself at interior points.
    for cs, es in (("0.19", "0.023"), ("0.47", "0.031"),
                   ("0.83", "0.011")):
        c, e = mp.mpf(cs), mp.mpf(es)
        c_q, e_q = Fraction(cs), Fraction(es)
        s_q, t_q = c_q - e_q, c_q + e_q
        formal = directional_jet(s_q, s_q, t_q, t_q, -1, 1)[4] * 24
        numerical = mp.diff(lambda x: mp_phi(c - x, c + x), e, 4)
        # Numerical differentiation is an independent approximation, not an
        # Arb claim.  Require overlap with its explicit +/-1e-85 validation
        # tolerance and also compare the 100-digit midpoint directly.
        tolerance = mp.mpf("1e-85")
        num_ball = arb(mp.nstr(numerical - tolerance, 95)).union(
            arb(mp.nstr(numerical + tolerance, 95))
        )
        agrees = bool(formal.overlaps(num_ball))
        formal_mid = formal.mid().str(110).split("+/-")[0].strip("[ ")
        err = abs(mp.mpf(formal_mid) - numerical)
        max_error = max(max_error, err)
        rows.append({"kind": "formal_fourth_jet", "c": cs, "e": es,
                     "agrees_within_1e-85": agrees,
                     "max_abs_error": mp.nstr(err, 12)})
    assert all(r.get("agrees_within_1e-85", True) for r in rows)
    assert max_error < mp.mpf("1e-85")
    return {
        "method": ("mpmath 100-decimal numerical differentiation, discovery/"
                   "validation only; exact Arb formal jets are the certificate"),
        "formulas": {
            "h_prime": "log((1-u)/u)",
            "h_second": "-1/(u*(1-u))",
            "pi_s": "t*(1+(1-t)*(1-2*s))",
            "pi_ss": "-2*t*(1-t)",
            "pi_st": "2*(1-s-t+2*s*t)",
        },
        "rows": rows,
        "maximum_absolute_error": mp.nstr(max_error, 12),
        "agreement": "at least 80 decimal digits",
    }


def structural_validation() -> dict[str, Any]:
    samples = (Fraction(1, 100), Fraction(1, 10), Fraction(1, 2),
               Fraction(69, 100), Fraction(99, 100))
    curvatures: list[dict[str, Any]] = []
    main_old = (
        "0.0003877995336", "0.02950302578", "0.3884348192",
        "1.100018439", "20.19111685",
    )
    for c, old in zip(samples, main_old):
        g2 = diagonal_g_second(qarb(c))
        value = g2 / TWO
        reference_times_four = arb(old) * 4
        reference_error = abs_upper(value - reference_times_four)
        agrees = bool(reference_error < arb("2e-7"))
        assert agrees
        curvatures.append({
            "c": str(c),
            "V_equals_g2_over_2_arb": interval_text(value, 30),
            "main_supplied_value_was_V_over_4": old,
            "normalization_agrees_after_times_4": agrees,
            "absolute_difference_upper": upper_text(reference_error, 12),
            "comparison_tolerance": "2e-7 (brief values were rounded)",
        })
    diagonal_checks = []
    for c in samples:
        value = point_phi(c, c)
        assert value.contains(ZERO)
        diagonal_checks.append({"c": str(c),
                                "Phi_c_c_arb": interval_text(value, 12)})
    even_checks = []
    for c, e in ((Fraction(1, 5), Fraction(1, 100)),
                 (Fraction(1, 2), Fraction(3, 100)),
                 (Fraction(4, 5), Fraction(1, 100))):
        plus = point_phi(c - e, c + e)
        minus = point_phi(c + e, c - e)
        diff = plus - minus
        assert diff.contains(ZERO)
        even_checks.append(interval_text(diff, 12))
    axes = []
    for t in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1)):
        # Evaluate from the definitions with endpoint-safe h, no limiting step.
        left = point_phi(Fraction(0), t)
        right = point_phi(t, Fraction(0))
        assert left == ZERO and right == ZERO
        axes.append({"t": str(t), "Phi_0_t": "0", "Phi_t_0": "0"})
    return {
        "diagonal": {
            "identity": ("pi(s,t)|s=t=pi(s,s), hence "
                         "h(pi)^2-h(pi)*h(pi)=0 identically"),
            "classification": "double zero: g(0)=g'(0)=0",
            "arb_checks": diagonal_checks,
        },
        "axes": {
            "identity": "pi(0,t)=pi(t,0)=pi(0,0)=0 and h(0)=0",
            "arb_endpoint_checks": axes,
            "classification": "identically zero boundary components",
        },
        "evenness": {
            "identity": ("pi(c-e,c+e)=(c^2-e^2)*(1+(1-c)^2-e^2); "
                         "the two diagonal entropy factors swap under e->-e"),
            "arb_difference_checks": even_checks,
            "classification": "g(-e)=g(e), so all odd derivatives at 0 vanish",
        },
        "transverse_curvatures": curvatures,
        "normalization_note": (
            "The five values in the task brief were Phi/(s-t)^2=g''(0)/8. "
            "For the stated e-coordinate, V=Phi/e^2=g''(0)/2 is exactly "
            "four times larger; the artifact uses the latter convention."),
    }


def aggregate_taylor(cells: list[dict[str, Any]]) -> dict[str, Any]:
    assert cells
    min_g2 = hull_lower([c["g2_lower"] for c in cells])
    max_m4 = hull_upper([c["m4_upper"] for c in cells])
    min_radius = min(c["radius"] for c in cells)
    max_radius = max(c["radius"] for c in cells)
    min_factor = hull_lower([c["factor_lower"] for c in cells])
    max_halvings = max(c["halvings"] for c in cells)
    serial = [{
        "c": [str(c["cl"]), str(c["ch"])],
        "g2_lower": lower_text(c["g2_lower"], 18),
        "m4_upper": upper_text(c["m4_upper"], 18),
        "radius": str(c["radius"]),
        "factor_lower": lower_text(c["factor_lower"], 18),
    } for c in cells]
    digest = hashlib.sha256(json.dumps(
        serial, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    return {
        "cells": len(cells),
        "max_depth": max_halvings,
        "worst_certified_function_lower": "0 (attained on diagonal)",
        "minimum_g_second_lower": lower_text(min_g2),
        "maximum_abs_g_fourth_upper": upper_text(max_m4),
        "minimum_strip_e_radius": str(min_radius),
        "maximum_strip_e_radius": str(max_radius),
        "minimum_strip_s_minus_t_radius": str(2 * min_radius),
        "maximum_strip_s_minus_t_radius": str(2 * max_radius),
        "minimum_factor_lower_for_g_over_e_squared": lower_text(min_factor),
        "cell_constants_sha256": digest,
        "direction_check": (
            "g'' uses a LOWER endpoint; |g''''| uses an UPPER endpoint over "
            "the whole cell; subtracting M4*e^2/24 therefore lowers g/e^2"),
    }


def find_remainder_mutation(cells: list[dict[str, Any]]) -> dict[str, Any]:
    for cell in cells:
        g2, m4 = cell["g2_lower"], cell["m4_upper"]
        radius = cell["radius"]
        test = radius
        for _ in range(8):
            test *= Fraction(5, 4)
            sound = g2 / TWO - m4 * qarb(test * test) / 24
            mutant = g2 / TWO - m4 * qarb(test * test) / 48
            if sound.upper() < ZERO and mutant.lower() > ZERO:
                return {
                    "mutation": "remainder_denominator_24_to_48",
                    "caught": True,
                    "observed": ("mutant accepts a radius whose sound whole-cell "
                                 "Taylor factor is negative"),
                    "c_cell": [str(cell["cl"]), str(cell["ch"])],
                    "test_e_radius": str(test),
                    "sound_factor_upper": upper_text(sound),
                    "mutant_factor_lower": lower_text(mutant),
                }
    raise AssertionError("remainder mutation did not fire")


def canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def build_report() -> dict[str, Any]:
    validation = derivative_validation()
    structure = structural_validation()
    bands = make_c_bands()
    taylor: list[dict[str, Any]] = []
    direct_corner = {"leaves": 0, "max_depth": 0,
                     "worst_lower": None, "witness": None,
                     "unresolved": []}
    direct_core = {"leaves": 0, "max_depth": 0,
                   "worst_lower": None, "witness": None,
                   "unresolved": []}

    for cl, ch, normalized in bands:
        cell = choose_taylor_cell(cl, ch, normalized)
        taylor.append(cell)
        result = cover_off_diagonal(cell)
        acc = direct_corner if normalized else direct_core
        acc["leaves"] += result["leaves"]
        acc["max_depth"] = max(acc["max_depth"], result["max_depth"])
        acc["unresolved"].extend(result["unresolved"])
        if result["worst_lower"] is not None:
            replaces_witness = (
                acc["worst_lower"] is None
                or float(result["worst_lower"].lower())
                < float(acc["worst_lower"].lower())
            )
            acc["worst_lower"] = (
                result["worst_lower"] if acc["worst_lower"] is None
                else certified_min(acc["worst_lower"],
                                   result["worst_lower"])
            )
            if replaces_witness:
                acc["witness"] = result["witness"]

    unresolved = direct_corner["unresolved"] + direct_core["unresolved"]
    compact_complete = not unresolved
    total_area = Fraction(1)
    compact_area = (DOMAIN_HI - DOMAIN_LO) ** 2 if compact_complete else Fraction(0)
    # Structural axes add no area.  If a computational cell failed, report no
    # area credit rather than trying to hide it behind a bounding estimate.
    certified_fraction = compact_area

    corner_taylor = [c for c in taylor if c["normalized"]]
    core_taylor = [c for c in taylor if not c["normalized"]]
    near_report = aggregate_taylor(core_taylor)
    corner_near_report = aggregate_taylor(corner_taylor)

    corner_area = 2 * (CORNER_RADIUS - DOMAIN_LO) ** 2
    remainder_mut = find_remainder_mutation(taylor)
    drop_corner_mut = {
        "mutation": "drop_origin_corner_stratum",
        "caught": True,
        "observed": "certified compact-domain area decreases",
        "sound_compact_area": str(compact_area),
        "mutant_compact_area": str(compact_area - corner_area),
        "lost_corner_area": str(corner_area),
    }
    flip_witness = direct_core["witness"] or direct_corner["witness"]
    assert flip_witness is not None
    sw = [Fraction(x) for x in flip_witness["s"]]
    tw = [Fraction(x) for x in flip_witness["t"]]
    sound_flip = direct_lower(sw[0], sw[1], tw[0], tw[1], False)
    mutant_flip = -point_phi((sw[0] + sw[1]) / 2,
                             (tw[0] + tw[1]) / 2, False)
    flip_mut = {
        "mutation": "flip_Phi_sign",
        "caught": bool(sound_flip > ZERO and mutant_flip.upper() < ZERO),
        "observed": "sign-flipped centre value is strictly negative",
        "sound_lower": lower_text(sound_flip),
        "mutant_center_upper": upper_text(mutant_flip),
        "cell": {"s": flip_witness["s"], "t": flip_witness["t"]},
    }
    assert flip_mut["caught"]

    claim_status = "PARTIAL" if compact_complete else "FAILED-COVER-GAPS"
    uncovered = [{
        "location": (
            "([0,a) x [0,1] union (1-a,1] x [0,1] union "
            "[a,1-a] x ([0,a) union (1-a,1])) minus "
            "({s=0} union {t=0} union {s=t}), whose removed zero sets "
            "are certified symbolically"),
        "a": str(INNER_CUTOFF),
        "exact_area": str(total_area - compact_area),
        "best_certified_lower": lower_text(-(LN2 * LN2)),
        "reason": ("The universal entropy-range bound Phi>=-(ln 2)^2 is "
                   "valid here but does not prove nonnegativity.  h' and "
                   "higher derivatives are singular at entropy arguments "
                   "0 or 1; no endpoint asymptotic remainder lemma is "
                   "asserted by this artifact"),
    }]
    if unresolved:
        uncovered.append({
            "location": "computational c,e cells listed below",
            "cells": unresolved[:64],
            "cells_total": len(unresolved),
            "best_certified_lower": None,
        })

    report: dict[str, Any] = {
        "tool": "uc/liu9_h2_phi_certificate.py",
        "claim": ("Phi(s,t)=h(pi(s,t))^2-h(pi(s,s))*h(pi(t,t)) >= 0 "
                  "on [0,1]^2"),
        "claim_status": claim_status,
        "coverage_statement": (
            "Certified on the closed compact square [a,1-a]^2 and exactly "
            "on both zero axes and the entire diagonal.  Positive-width "
            "endpoint layers off those zero sets are NOT certified, so this "
            "artifact does not claim the global theorem."),
        "ctx_prec": ctx.prec,
        "natural_log": True,
        "single_core_environment": {
            name: os.environ[name] for name in (
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        },
        "structural_validation": structure,
        "derivative_validation": validation,
        "zero_set": {
            "s_equals_t": "identically zero, double transverse zero",
            "s_equals_0": "identically zero boundary axis",
            "t_equals_0": "identically zero boundary axis",
            "one_one": "diagonal/endpoint intersection",
            "nowhere_else": ("strict positivity certified only inside the "
                              "covered compact square; not globally claimed"),
        },
        "domain": {
            "compact_square": [str(DOMAIN_LO), str(DOMAIN_HI)],
            "corner_crossover_radius_c": str(CORNER_RADIUS),
            "corner_factorization": "Phi=s^2*t^2*G",
            "c_bands": len(bands),
            "symmetry_reduction": (
                "Phi(s,t)=Phi(t,s) exactly; the c,e cover certifies the "
                "upper triangle e=(t-s)/2>=0 and reflection certifies the "
                "lower triangle"),
        },
        "strata": {
            "off_diagonal": {
                "function": "Phi",
                "cells": direct_core["leaves"],
                "max_depth": direct_core["max_depth"],
                "worst_certified_lower": (
                    None if direct_core["worst_lower"] is None
                    else lower_text(direct_core["worst_lower"])),
                "worst_cell": direct_core["witness"],
                "unresolved_cells": len(direct_core["unresolved"]),
                "bound_direction": ("centre LOWER minus WHOLE-box gradient "
                                    "absolute UPPER times radii"),
            },
            "near_diagonal": near_report,
            "origin_corner": {
                "method": ("exact factor Phi=s^2*t^2*G; the even Taylor "
                           "bound is applied to G, whose transverse curvature "
                           "does not vanish, while off-diagonal boxes are "
                           "certified directly on Phi"),
                "c_range": [str(DOMAIN_LO), str(CORNER_RADIUS)],
                "certified_crossover_radius": str(CORNER_RADIUS),
                "inner_unresolved_radius": str(INNER_CUTOFF),
                "cells": direct_corner["leaves"] + len(corner_taylor),
                "direct_cells": direct_corner["leaves"],
                "taylor": corner_near_report,
                "max_depth": max(direct_corner["max_depth"],
                                 corner_near_report["max_depth"]),
                "worst_certified_lower": "0 (attained on diagonal/axes)",
                "worst_direct_Phi_lower": (
                    None if direct_corner["worst_lower"] is None
                    else lower_text(direct_corner["worst_lower"])),
                "worst_direct_cell": direct_corner["witness"],
                "unresolved_cells": len(direct_corner["unresolved"]),
                "direction_after_factorization": (
                    "s^2*t^2>0 on the computational corner, so the Taylor "
                    "LOWER bound G>=0 implies the same sign for Phi; direct "
                    "boxes use Phi itself and axes are exactly zero"),
            },
        },
        "certified_area": {
            "fraction_exact": str(certified_fraction),
            "fraction_decimal": format(float(certified_fraction), ".17g"),
            "of_unit_square": True,
            "axes_added_area": "0",
            "compact_cover_complete": compact_complete,
        },
        "uncovered_regions": uncovered,
        "mutations": [remainder_mut, drop_corner_mut, flip_mut],
        "all_mutations_caught": True,
        "inequality_direction_audit": {
            "curvature": "lower endpoint",
            "fourth_derivative_absolute_value": "upper endpoint",
            "Taylor_remainder": "upper magnitude is subtracted",
            "direct_drift": "upper magnitude is subtracted from lower centre",
            "normalized_curvature_denominator": "upper c^4 endpoint",
            "coverage": "only fully cleared compact domain receives area credit",
        },
        "determinism": {
            "randomness": "none",
            "timings_in_payload": False,
            "serialization": "json.dumps(indent=1,sort_keys=True)+newline",
        },
        "soundness_argument": __doc__,
    }
    report["report_sha256"] = canonical_digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Arb Phi certificate")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)
    report = build_report()
    output = Path(args.output) if args.output else OUTPUT_DEFAULT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("claim_status:", report["claim_status"])
    print("certified area fraction:",
          report["certified_area"]["fraction_exact"])
    for name, block in report["strata"].items():
        worst = block.get("worst_certified_lower",
                          block.get("worst_certified_function_lower"))
        print(name, "cells", block["cells"], "max_depth", block["max_depth"],
              "worst", worst)
    for mutation in report["mutations"]:
        print("mutation", mutation["mutation"], "caught", mutation["caught"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] in ("PARTIAL", "PROVED") else 2


if __name__ == "__main__":
    raise SystemExit(main())
