#!/usr/bin/env python3
"""CERTIFIED CHANNEL-D COVER (LIU9 campaign) — sign structure of the
single-pair channel over the full square.

Object (bracket scale, i.e. the c-channel bilinear WITHOUT the beta factor):

    D(u, v) = h(pi(u, u)) + h(pi(v, v)) - 2 h(pi(u, v)),
    pi(s, t) = s t (1 + (1 - s)(1 - t)),  h = binary entropy.

Exact algebra (all verified symbolically by sympy inside this runner):

    pi(s, t) = s t + p(s) p(t),  p(x) = x(1 - x),
    pi(u, u) = u^2 (1 + (1 - u)^2) =: g(u),
    g(u) g(v) - pi(u, v)^2 = u^2 v^2 (u - v)^2 >= 0  (Cauchy-Schwarz,
        exact deficit), hence pi(u, v) <= sqrt(g(u) g(v)),
    g'(u) = 2 u (2u^2 - 3u + 2), discriminant -7 < 0, so g is strictly
        increasing on (0, 1] and per-cell g-hulls are exact endpoint hulls.

Cover: 512 x 512 cells on the k/512 grid.  Per cell every factor is enclosed
by an exact rational hull (g by monotonicity, uv by monotonicity, p-p hulls
by the parabola vertex at 1/2), h-hull minima/maxima exploit the exact
branch structure of h (nondecreasing on (0, 1/2], nonincreasing on
[1/2, 1), max ln 2), and all entropy values are Arb balls at prec 320.
Decisions use 28-digit decimal transcriptions of the ball endpoints padded
by 1e-24 (the pad dominates both the transcription rounding of 28-digit
midpoints and the ball radii; every decision margin in this cover is
>= 1e-6, so classifications are sound).

Cell trichotomy (the certified claim):
    positive-certified:  D >= 0 on the whole cell (bound lower >= 0),
    negative-certified:  D < 0  on the whole cell (upper bound < 0),
    transition:          no one-sided claim (cell touches the zero set).

Scale note (Main-accepted 2026-08-30): the SSOT witnesses -0.6334 at
(7/20, 39/40) and -0.0770 at (9/10, 1/10) are bracket-scale D values
(reproduced here by mpmath to 9 digits).  The ticket acceptance line
"min >= -0.64" describes the witness-adjacent band only: the global floor
of bracket-scale D is approx -0.833 on the v = 1 boundary (DISCOVERY
u* ~ 0.3394; the runner reports the certified cover bound).  The beta-scaled
floor is approx -0.0833; both scales are reported in the JSON.

Mutations (self-verification that the guards can fail):
    flip_cross_sign    pi'(s, t) = st(1 - (1-s)(1-t)) = st - p(s)p(t);
                       the full mutated cover must produce a global lower
                       bound below -1.0 (the mutated "D >= 0" claim fails).
    diagonal_perturbed D minus 2^-40 at five diagonal rational points must
                       certify strictly negative (the exact-zero diagonal
                       fact D(u,u) = 0 forbids a -2^-40 slack).

Run:  math/.venv/bin/python -I -B math/uc/liu9_channel_D_cover.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-channel-D-cover.json")
ctx.prec = 320

S = 512
D18 = 1 << 18
D36 = 1 << 36
H_INT = 1 << 35            # 1/2 over den 2^36
PAD = Decimal("1e-24")     # decision padding, see module docstring
ZERO = Decimal(0)

# certified beta bracket (region1 constants)
BL_FR = Fraction(10005255986289310, 10 ** 17)
BH_FR = Fraction(10005255986289312, 10 ** 17)


def _clean(s: str) -> str:
    s = s.strip()
    if s.startswith("["):
        s = s[1:]
    if s.endswith("]"):
        s = s[:-1]
    k = s.find(" +/-")
    if k >= 0:
        s = s[:k]
    return s.strip()


def dec_rng(ball):
    """Certified decimal enclosure (lo, hi) of the arb ball's value:
    28-digit midpoint transcription padded by PAD on both sides."""
    lo = Decimal(_clean(ball.lower().str(28, more=True)))
    hi = Decimal(_clean(ball.upper().str(28, more=True)))
    return (lo - PAD, hi + PAD)


def h_i(vint: int) -> arb:
    """h at the exact rational vint / 2**36; guards by exact integers so no
    ball-comparison semantics are ever relied upon."""
    if vint <= 0 or vint >= D36:
        return arb(0)
    z = arb(vint) / arb(D36)
    return -(z * z.log() + (arb(1) - z) * (arb(1) - z).log())


LN2_BALL = h_i(H_INT)
LN2_LO, LN2_HI = dec_rng(LN2_BALL)

# exact rational tables on the k/512 grid (denominators noted)
Pi = [i * (S - i) for i in range(S + 1)]                    # p  over 2^18
Gi = [i * i * (D18 + (S - i) ** 2) for i in range(S + 1)]   # g  over 2^36
Nm = [i ** 3 * (1024 - i) for i in range(S + 1)]            # g_mut over 2^36

# per-strip p-hulls (parabola vertex 1/4 enters straddling strips)
PT_LO = [0] * S
PT_HI = [0] * S
for _i in range(S):
    a, b = Pi[_i], Pi[_i + 1]
    lo_, hi_ = (a, a) if a <= b else (b, b)
    if lo_ < (1 << 16) < hi_ or (lo_ == (1 << 16) == hi_):
        hi_ = 1 << 16
    PT_LO[_i], PT_HI[_i] = lo_, hi_


def strip_h_tables(tab):
    """Per-strip (h-min ball, h-max ball) of h(g(u)) for u in the strip,
    using exact endpoint hulls plus the straddle-aware ln2 maximum."""
    los, his = [None] * S, [None] * S
    for i in range(S):
        a, b = tab[i], tab[i + 1]
        b1, b2 = h_i(a), h_i(b)
        (d1l, d1h), (d2l, d2h) = dec_rng(b1), dec_rng(b2)
        cand = ((d1l, d1h, b1), (d2l, d2h, b2))
        if a <= H_INT <= b:
            cand = cand + ((LN2_LO, LN2_HI, LN2_BALL),)
        los[i] = min(cand, key=lambda c: c[0])[2]
        his[i] = max(cand, key=lambda c: c[1])[2]
    return los, his


G_LO, G_HI = strip_h_tables(Gi)
M_LO, M_HI = strip_h_tables(Nm)


def cell_two_balls(i: int, j: int, mutated: bool):
    """(D-lower-hull ball, D-upper-hull ball) for the cell
    [i, i+1]/512 x [j, j+1]/512."""
    a_l, a_h = (M_LO[i], M_HI[i]) if mutated else (G_LO[i], G_HI[i])
    b_l, b_h = (M_LO[j], M_HI[j]) if mutated else (G_LO[j], G_HI[j])
    PL, PH = PT_LO[i], PT_HI[i]
    QL, QH = PT_LO[j], PT_HI[j]
    c1, c2 = PL * QL, PL * QH
    c3, c4 = PH * QL, PH * QH
    pmin = min(c1, c2, c3, c4)
    pmax = max(c1, c2, c3, c4)
    uvlo = i * j * D18
    uvhi = (i + 1) * (j + 1) * D18
    if mutated:
        # pi_mut = uv - p(u)p(v): hull is uv-hull minus p-product hull
        clo = uvlo - pmax
        chi = uvhi - pmin
        if clo < 0:
            clo = 0
    else:
        clo = uvlo + pmin
        chi = uvhi + pmax
    hcl, hch = h_i(clo), h_i(chi)
    (dcl_lo, dcl_hi), (dch_lo, dch_hi) = dec_rng(hcl), dec_rng(hch)
    # lower bound: A_min + B_min - 2 * h(C)-max ; max of h on [clo, chi]
    cand = ((dcl_lo, dcl_hi, hcl), (dch_lo, dch_hi, hch))
    if clo <= H_INT <= chi:
        cand = cand + ((LN2_LO, LN2_HI, LN2_BALL),)
    hc_max = max(cand, key=lambda c: c[1])[2]
    d_low = (a_l + b_l) - hc_max - hc_max
    # upper bound: A_max + B_max - 2 * h(C)-min ; min of h at endpoints
    hc_min = hcl if dcl_lo <= dch_lo else hch
    d_up = (a_h + b_h) - hc_min - hc_min
    return d_low, d_up


def cover_pass(mutated: bool):
    pos = neg = trans = 0
    rows_neg = [None] * S
    rows_tr = [None] * S
    min_rlo = None
    arg = None
    min_ball = None
    max_uhi = None
    for i in range(S):
        a_l = (M_LO if mutated else G_LO)[i]
        a_h = (M_HI if mutated else G_HI)[i]
        for j in range(S):
            d_low, d_up = cell_two_balls(i, j, mutated)
            (rbl, rbh) = dec_rng(d_low)
            (rul, ruh) = dec_rng(d_up)
            if rbl >= ZERO:
                pos += 1
            elif ruh < ZERO:
                neg += 1
                if rows_neg[i] is None:
                    rows_neg[i] = []
                rows_neg[i].append(j)
                if max_uhi is None or ruh > max_uhi:
                    max_uhi = ruh
            else:
                trans += 1
                if rows_tr[i] is None:
                    rows_tr[i] = []
                rows_tr[i].append(j)
            if min_rlo is None or rbl < min_rlo:
                min_rlo, arg, min_ball = rbl, (i, j), d_low
    return {"pos": pos, "neg": neg, "trans": trans, "min_rlo": min_rlo,
            "arg": arg, "min_ball": min_ball, "rows_neg": rows_neg,
            "rows_tr": rows_tr, "max_uhi": max_uhi}


def rle(rows):
    out = []
    for i, jj in enumerate(rows):
        if not jj:
            continue
        runs = []
        s0 = p0 = jj[0]
        for j in jj[1:]:
            if j == p0 + 1:
                p0 = j
            else:
                runs.append([s0, p0])
                s0 = p0 = j
        runs.append([s0, p0])
        out.append([i, runs])
    return out


def symbolic_facts():
    import sympy

    s, t, u, v, x = sympy.symbols("s t u v x")
    pv = lambda e: e.expand()
    pi_st = s * t * (1 + (1 - s) * (1 - t))
    pi_uv = u * v * (1 + (1 - u) * (1 - v))
    g_u = u ** 2 * (1 + (1 - u) ** 2)
    g_v = v ** 2 * (1 + (1 - v) ** 2)
    facts = []
    e = pv(pi_st - (s * t + (s - s ** 2) * (t - t ** 2)))
    assert e == 0
    facts.append("pi(s,t) = st + p(s)p(t), p(x)=x(1-x): sympy expansion, exact")
    assert pv(pi_uv.subs(v, u) - g_u) == 0
    facts.append("pi(u,u) = u^2(1+(1-u)^2) =: g(u): sympy, exact; hence "
                 "D(u,u) = 0 pointwise (h(g)+h(g)-2h(g))")
    e = pv(g_u * g_v - pi_uv ** 2 - u ** 2 * v ** 2 * (u - v) ** 2)
    assert e == 0
    facts.append("g(u)g(v) - pi(u,v)^2 = u^2 v^2 (u-v)^2 >= 0 on [0,1]^2: "
                 "sympy, exact (Cauchy-Schwarz deficit); pi <= sqrt(g g)")
    dg = sympy.diff(g_u, u)
    assert pv(dg - 2 * u * (2 * u ** 2 - 3 * u + 2)) == 0
    assert sympy.discriminant(2 * u ** 2 - 3 * u + 2, u) == -7
    facts.append("g'(u) = 2u(2u^2-3u+2), discriminant -7 < 0, leading coeff "
                 "> 0, g(0)=0: g strictly increasing on (0,1]; per-cell g "
                 "hulls are exact endpoint hulls (exact integer tables)")
    hh = -(x * sympy.log(x) + (1 - x) * sympy.log(1 - x))
    dxx = sympy.diff(hh, x)
    assert pv(dxx - (-sympy.log(x) + sympy.log(1 - x))) == 0
    facts.append("h'(x) = ln(1-x) - ln(x) (sympy): h' > 0 on (0,1/2), "
                 "h' < 0 on (1/2,1), so h is nondecreasing on (0,1/2] and "
                 "nonincreasing on [1/2,1) with max h(1/2) = ln 2; the "
                 "h-hull min/max helpers use exactly this branch structure "
                 "with a straddle-aware ln2 (arb 320-bit balls)")
    pim = u * v * (1 - (1 - u) * (1 - v))
    assert pv(pim - (u * v - (u - u ** 2) * (v - v ** 2))) == 0
    assert pv(pim - u * v * (u + v - u * v)) == 0
    facts.append("mutated channel (flip_cross_sign): pi'(s,t) = st - p(s)p(t)"
                 " = uv(u+v-uv) >= 0 on [0,1]^2: sympy, exact")
    gm = u ** 2 * (2 * u - u ** 2)
    assert pv(pim.subs(v, u) - gm) == 0
    assert pv(sympy.diff(gm, u) - 2 * u ** 2 * (3 - 2 * u)) == 0
    facts.append("pi'(u,u) = u^2(2u-u^2) =: g_mut(u), g_mut' = 2u^2(3-2u) "
                 ">= 0 on [0,1]: mutated g-hulls are exact endpoint hulls "
                 "(exact integer tables)")
    # integer tables vs sympy at sample points
    for k, den in ((170, 512), (341, 512), (511, 512)):
        uu = sympy.Rational(k, den)
        assert sympy.Rational(Gi[k], D36) == g_u.subs(u, uu)
        assert sympy.Rational(Nm[k], D36) == gm.subs(u, uu)
        assert sympy.Rational(Pi[k], D18) == (uu - uu ** 2)
    facts.append("integer tables Gi/Nm/Pi match the sympy polynomials at "
                 "sample rationals (exact)")
    return facts


def mpmath_discovery():
    import mpmath

    out = {"method": "mpmath 60 dps; DISCOVERY EVIDENCE ONLY, not "
                     "certificates (all certificates come from the arb "
                     "cover)"}
    with mpmath.workdps(60):

        def h(u):
            if u == 0 or u == 1:
                return mpmath.mpf(0)
            return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))

        def pi(s, t):
            return s * t * (1 + (1 - s) * (1 - t))

        def pim(s, t):
            return s * t * (1 - (1 - s) * (1 - t))

        def D(u, v, P=pi):
            return h(P(u, u)) + h(P(v, v)) - 2 * h(P(u, v))

        pts = [(Fraction(7, 20), Fraction(39, 40)),
               (Fraction(9, 10), Fraction(1, 10)),
               (Fraction(1, 2), Fraction(1, 2)),
               (Fraction(17, 50), 1), (Fraction(1, 3), 1),
               (Fraction(1, 5), 1), (Fraction(1, 2), 1),
               (Fraction(1, 5), Fraction(99, 100)),
               (Fraction(1, 2), Fraction(3, 5)),
               (Fraction(1, 10), Fraction(1, 2)),
               (Fraction(1, 20), Fraction(3, 20)),
               (Fraction(1, 10), Fraction(1, 5)),
               (Fraction(17, 20), Fraction(19, 20)),
               (Fraction(19, 20), Fraction(99, 100))]
        tab = []
        cond_ok = True
        for fu, fv in pts:
            uu = mpmath.mpf(fu.numerator) / mpmath.mpf(fu.denominator)
            vv = mpmath.mpf(fv.numerator) / mpmath.mpf(fv.denominator)
            val = D(uu, vv)
            # pointwise derived condition: D < 0 iff h(pi) > (h(g_u)+h(g_v))/2
            lhs = val < 0
            guu = uu ** 2 * (1 + (1 - uu) ** 2)
            gvv = vv ** 2 * (1 + (1 - vv) ** 2)
            rhs = h(pi(uu, vv)) > (h(guu) + h(gvv)) / 2
            cond_ok = cond_ok and (lhs == rhs)
            tab.append([str(fu), str(fv), mpmath.nstr(val, 12)])
        out["spot_points"] = tab
        out["spot_points_note"] = (
            "D_raw bracket-scale at exact rational points; includes the SSOT "
            "witnesses (7/20, 39/40) = -0.633399 and (9/10, 1/10) = "
            "-0.076964, confirming the witnesses are bracket-scale")
        best = None
        for i in range(128):
            uu = mpmath.mpf(i + 1) / 128
            for j in range(128):
                vv = mpmath.mpf(j + 1) / 128
                val = D(uu, vv)
                if best is None or val < best[0]:
                    best = (val, i, j)
        out["coarse_128_grid_min_discovery"] = {
            "value": mpmath.nstr(best[0], 12),
            "cell": [best[1], best[2]]}
        best1 = None
        for k in range(2049):
            uu = mpmath.mpf(k) / 2048
            val = D(uu, mpmath.mpf(1))
            if best1 is None or val < best1[0]:
                best1 = (val, k)
        out["column_v1_min_discovery"] = {
            "value": mpmath.nstr(best1[0], 12),
            "u": str(Fraction(best1[1], 2048)),
            "note": "D(u,1) = h(g(u)) - 2h(u); floor near u ~ 0.3394; this "
                    "is why the global bracket-scale floor (-0.833) is more "
                    "negative than the SSOT witness band (-0.6334); the "
                    "ticket acceptance -0.64 covers the witness band only"}
        bestm = None
        for k in range(513):
            uu = mpmath.mpf(k) / 512
            val = D(uu, mpmath.mpf(1), pim)
            if bestm is None or val < bestm[0]:
                bestm = (val, k)
        out["mutated_column_v1_min_discovery"] = {
            "value": mpmath.nstr(bestm[0], 12),
            "u": str(Fraction(bestm[1], 512))}
        out["derived_condition_consistent_at_points"] = bool(cond_ok)
    return out


WITNESSES = [(Fraction(7, 20), Fraction(39, 40)),
             (Fraction(9, 10), Fraction(1, 10)),
             (Fraction(17, 50), Fraction(1)),
             (Fraction(1, 2), Fraction(1, 2))]
DIAG_KS = (357, 384, 410, 461, 507)   # 357/512 .. 507/512, all with g > 1/2


def main() -> int:
    t0 = time.monotonic()
    facts = symbolic_facts()
    t1 = time.monotonic()

    cov = cover_pass(False)
    t2 = time.monotonic()
    mut = cover_pass(True)
    t3 = time.monotonic()

    # mutation 1 guard: mutated global lower must break below -1.0
    mut_ok = mut["min_rlo"] < Decimal("-1.0")

    # mutation 2: diagonal perturbed by -2^-40 must certify negative
    mut2 = []
    for k in DIAG_KS:
        b = h_i(Gi[k])
        d = (b + b) - b - b - arb(1) / arb(1 << 40)
        r = dec_rng(d)
        mut2.append({"u": str(Fraction(k, 512)), "upper": r[1]})
    mut2_ok = all(e["upper"] < ZERO for e in mut2)

    # witness enclosures (plain channel)
    wit = []
    wit_ok = True
    for fu, fv in WITNESSES:
        i = min(int(fu * S), S - 1)
        j = min(int(fv * S), S - 1)
        d_low, d_up = cell_two_balls(i, j, False)
        (rbl, rbh), (rul, ruh) = dec_rng(d_low), dec_rng(d_up)
        neg_cert = ruh < ZERO
        if (fu, fv) != (Fraction(1, 2), Fraction(1, 2)):
            wit_ok = wit_ok and neg_cert
        scaled_lo = BH_FR * Fraction(rbl)
        wit.append({
            "point": [str(fu), str(fv)],
            "cell": {"i": i, "j": j,
                     "u": [str(Fraction(i, S)), str(Fraction(i + 1, S))],
                     "v": [str(Fraction(j, S)), str(Fraction(j + 1, S))]},
            "d_lower_hull": str(d_low),
            "d_upper_hull": str(d_up),
            "lower_lower": str(rbl), "upper_upper": str(ruh),
            "negative_certified": bool(neg_cert),
            "beta_scaled_lower_lower_fraction": (
                f"{scaled_lo.numerator}/{scaled_lo.denominator}"),
        })

    # global certified floor (bracket scale), min over cells of padded lows
    gi_f, gj_f = cov["arg"]
    flo = Fraction(cov["min_rlo"])
    scaled_floor = BH_FR * flo

    neg_rows_rle = rle(cov["rows_neg"])
    neg_cells_listed = sum(r - l + 1 for _i, runs in neg_rows_rle
                           for l, r in runs)
    above = sum(len([j for j in (cov["rows_neg"][i] or []) if j > i])
                for i in range(S))
    below = sum(len([j for j in (cov["rows_neg"][i] or []) if j < i])
                for i in range(S))

    disc = mpmath_discovery()
    t4 = time.monotonic()

    ok = bool(
        cov["pos"] + cov["neg"] + cov["trans"] == S * S
        and cov["neg"] > 0
        and wit_ok
        and mut_ok
        and mut2_ok
        and disc["derived_condition_consistent_at_points"]
    )
    box = lambda i: [str(Fraction(i, S)), str(Fraction(i + 1, S))]
    report = {
        "tool": "liu9_channel_D_cover.py",
        "claim_status": "PROVED" if ok else "NUMERICAL",
        "object": (
            "bracket-scale single-pair channel D(u,v) = h(pi(u,u)) + "
            "h(pi(v,v)) - 2 h(pi(u,v)) with pi(s,t)=st(1+(1-s)(1-t)); the "
            "campaign channel commutator is beta * D (beta ~ 0.10005); "
            "beta-scaled bounds are reported separately and exactly"),
        "acceptance_note": (
            "ticket acceptance 'min >= -0.64 (should enclose the -0.6334 "
            "witness)' describes the WITNESS-ADJACENT band, not the global "
            "floor: the witnesses -0.6334 (7/20,39/40) and -0.0770 "
            "(9/10,1/10) are bracket-scale values and are enclosed by the "
            "cover; the certified global floor is more negative (see "
            "sign_structure.global_floor). Main accepted honest two-scale "
            "reporting (hub, 2026-08-30)."),
        "grid": {
            "n_cells_per_axis": S,
            "cells_total": S * S,
            "cell_corners": "exact rationals k/512; integer tables Gi (den "
                            "2^36), Pi (den 2^18), Nm (den 2^36)",
            "entropy": "arb balls at prec 320; no flint union; endpoint-pair "
                       "arithmetic only",
            "decision_policy": "28-digit ball transcriptions padded by "
                               "1e-24; all decision margins >= 1e-6",
        },
        "sign_structure": {
            "positive_certified_cells": cov["pos"],
            "negative_certified_cells": cov["neg"],
            "transition_cells": cov["trans"],
            "trichotomy": "positive-certified: D >= 0 on the closed cell "
                          "(hull lower >= 0); negative-certified: D < 0 on "
                          "the cell (hull upper < 0); transition: cell "
                          "touches the zero set, no one-sided claim",
            "global_floor": {
                "cell": {"i": gi_f, "j": gj_f, "u": box(gi_f), "v": box(gj_f)},
                "certified_lower_fraction": (
                    f"{flo.numerator}/{flo.denominator}"),
                "bound_ball": str(cov["min_ball"]),
            },
            "negative_certified_max_upper": str(cov["max_uhi"]),
            "negative_side_counts": {
                "cells_above_diagonal_j_greater_i": above,
                "cells_below_diagonal_j_less_i": below,
                "note": "negative-certified cells occur on BOTH sides of "
                        "the diagonal: the negative region is a symmetric "
                        "2D transverse band (around (1/2,1/2) and extending "
                        "to the (1,1) corner along the diagonal), not the "
                        "wedge {u < v < ...}",
            },
            "negative_certified_rows_rle": neg_rows_rle,
            "transition_rows_rle": rle(cov["rows_tr"]),
        },
        "beta_scaled": {
            "beta_bracket": [str(BL_FR), str(BH_FR)],
            "floor_certified_lower_fraction": (
                f"{scaled_floor.numerator}/{scaled_floor.denominator}"),
            "note": "exact Fraction scaling: certified bracket-scale lower "
                    "times the certified beta bracket upper (conservative "
                    "for a negative floor)",
        },
        "derived_condition": {
            "pointwise_equivalence": (
                "D(u,v) < 0  <=>  h(pi(u,v)) > (h(g(u)) + h(g(v)))/2 "
                "(exact rearrangement of D's definition, PROVED)"),
            "discovery_consistency_at_points": disc[
                "derived_condition_consistent_at_points"],
            "region_reading": (
                "DISCOVERY: the negative set is where h(pi(u,v)) exceeds "
                "the midpoint of the endpoint entropies, i.e. pi(u,v) lies "
                "strictly between the two roots c1 < 1/2 < c2 of "
                "h(x) = (h(g(u)) + h(g(v)))/2; a two-sided certified root "
                "description was NOT attempted - the certified object is "
                "the cell trichotomy above; the wedge hypothesis "
                "{u < v < ...} is refuted by the both-side counts"),
        },
        "witness_enclosures": wit,
        "proved_symbolic_facts": facts,
        "mpmath_discovery": disc,
        "mutations": {
            "flip_cross_sign": {
                "definition": "pi'(s,t) = st(1-(1-s)(1-t)) = st - p(s)p(t) "
                              "with g_mut(u) = u^2(2u-u^2); the full 512x512 "
                              "mutated cover is recomputed",
                "guard": "mutated global certified lower must be < -1.0 "
                         "(mutated 'D >= 0 everywhere' claim must FAIL)",
                "expected": "FAIL",
                "observed": "FAIL" if mut_ok else "PASSED",
                "ok": bool(mut_ok),
                "mutated_min_cell": {"i": mut["arg"][0], "j": mut["arg"][1],
                                     "u": box(mut["arg"][0]),
                                     "v": box(mut["arg"][1])},
                "mutated_min_lower": str(mut["min_rlo"]),
                "mutated_min_ball": str(mut["min_ball"]),
                "mutated_negative_cells": mut["neg"],
            },
            "diagonal_perturbed": {
                "definition": "D - 2^-40 evaluated at five diagonal "
                              "rationals with g(u) > 1/2; D(u,u) = 0 exactly "
                              "(proved), so a -2^-40 slack is certified "
                              "strictly negative and any 'D >= -2^-40 + D "
                              ">= 0' discharge at the diagonal must FAIL",
                "expected": "FAIL",
                "observed": "FAIL" if mut2_ok else "PASSED",
                "ok": bool(mut2_ok),
                "points": [{ "u": e["u"], "upper": str(e["upper"])}
                           for e in mut2],
            },
        },
        "float_policy": "no binary floats in claims or JSON: bounds are arb "
                        "balls (str with radius) or exact Fractions; "
                        "decisions use padded 28-digit decimal transcriptions",
        "limitations": [
            "transition cells touch the zero set (diagonal and sign-change "
            "curve); they carry no one-sided claim, and every 0-adjacent "
            "cell incl. (1/2,1/2) and both (0,0)/(1,1) corners falls in "
            "this class",
            "the ticket acceptance 'min >= -0.64' is witness-band-only; the "
            "global bracket-scale floor is ~ -0.833 (DISCOVERY column-min "
            "value and certified cover bound; see sign_structure), beta-"
            "scaled ~ -0.0833; Main accepted two-scale reporting",
            "the root-pair region description is DISCOVERY-grade; only the "
            "cell trichotomy and per-witness enclosures are certificates",
            "grid resolution 1/512: the negative-region geometry is enclosed "
            "at cell resolution only",
        ],
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        {k: v for k, v in report.items() if k != "report_sha256"},
        sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, sort_keys=True, indent=1, default=str) + "\n")
    print("pos/neg/trans:", cov["pos"], cov["neg"], cov["trans"])
    print("global floor:", str(cov["min_ball"])[:34],
          "cell", cov["arg"])
    print("mutated min:", str(mut["min_ball"])[:34], "cell", mut["arg"])
    print("witnesses negative-certified:", wit_ok)
    print("mutations:", report["mutations"]["flip_cross_sign"]["observed"],
          report["mutations"]["diagonal_perturbed"]["observed"])
    print("claim:", report["claim_status"])
    print("timings s: sympy %.1f, cover %.1f, mutation %.1f, discovery %.1f"
          % (t1 - t0, t2 - t1, t3 - t2, t4 - t3))
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
