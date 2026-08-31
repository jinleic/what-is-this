from core import Prec, PREC
from d3h_kernels import (psi3f, psi3pf, psi3ppf, psi3pppf, Kf, eval_terms, mult_terms,
                         _KIND_CHAINS)
from flint import acb, arb


def bump_kind(kind: str, ch: str) -> str:
    p = kind.count('r'); i = kind.count('u'); j = kind.count('v')
    assert p + i + j == len(kind) - 1 and kind[0] == 'K'
    if ch == 'u':
        i += 1
    else:
        j += 1
    return "K" + "r" * p + "u" * i + "v" * j


def d_weight(w: tuple) -> tuple:
    """d/dx of a weight-list (levels 1..3): sum-over-factors; returns list of new weight tuples.
    Level-3 factor (psi3''' = sqrt6 const) drops under d/dx."""
    out = []
    for idx in range(len(w)):
        lvl = w[idx]
        if lvl >= 3:
            continue
        out.append(w[:idx] + (lvl + 1,) + w[idx + 1:])
    return out


def weight_eval(w: tuple, x: acb, bits: int = PREC) -> acb:
    tot = acb(1)
    for lvl in w:
        if lvl == 1:
            tot = tot * psi3pf(x)
        elif lvl == 2:
            tot = tot * psi3ppf(x)
        elif lvl == 3:
            tot = tot * psi3pppf(x)
    return tot


def dxdy_terms(a: int):
    if a == 0:
        return [((), (), "K", 0)]
    """Symbolic entries for (d/dx d/dy)^a C_r(x,y), a >= 1, as list of (wx, wy, kind, theta_power).
    Value form:  2pi * theta^(theta_power) * Wx(x) Wy(y) * mult(kind)(u,v) * K_r(u,v).
    Base (48) at a=1: 2pi theta^2 psi3'(x) psi3'(y) K  -> entry weights ((1,), (1,)) tp=2."""
    entries = [((1,), (1,), "K", 2)]
    for _ in range(a - 1):
        # apply ONE d/dx to each entry (2 branches), then ONE d/dy to each result (2 branches)
        after_x = []
        for (wx, wy, kind, tp) in entries:
            for nw in d_weight(wx):
                after_x.append((nw, wy, kind, tp))
            after_x.append((tuple(wx) + (1,), wy, bump_kind(kind, 'u'), tp + 1))
        entries2 = []
        for (wx, wy, kind, tp) in after_x:
            for nw in d_weight(wy):
                entries2.append((wx, nw, kind, tp))
            entries2.append((wx, tuple(wy) + (1,), bump_kind(kind, 'v'), tp + 1))
        entries = entries2
    return entries

def dxdy_eval(a: int, x: acb, y: acb, r: acb, vt: arb, bits: int = PREC) -> acb:
    """Numerically evaluate (d/dx d/dy)^a C_r(x,y) from symbolic entries."""
    with Prec(bits):
        u = acb(vt) * psi3f(x)
        v = acb(vt) * psi3f(y)
        tot = acb(0)
        for (wx, wy, kind, tp) in dxdy_terms(a):
            terms = mult_terms(kind)
            mval = eval_terms(terms, u, v, r, bits)
            wpair = weight_eval(wx, x, bits) * weight_eval(wy, y, bits)
            tot = tot + acb(2) * acb.pi() * (acb(vt) ** tp) * wpair * mval * Kf(u, v, r, bits)
        return tot
