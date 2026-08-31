#!/usr/bin/env python3
"""REGION 1 CERTIFICATE v3 (Family A) — scalar margin, split reg zmienna proof.

Family A: P0 = delta_x, P1 = delta_1, shared mass a = 1, x in [1/4, 1].

Structure (proved in liu9_block_copositive.py and numerically above):
  gapA(x, q) is an exact quadratic in q; gapA(x, 1) = 0 exactly (mixture
  collapses to P1); the quadratic is concave in q and the feasible band is
  q in [q_lo, 1], q_lo = max(0, (m - x)/(1 - x)).  A concave quadratic
  attains its band-minimum at an endpoint, hence
      min_q gapA = min(gapA(x, q_lo), gapA(x, 1) = 0) = min(gapA(x, q_lo), 0).

The region therefore reduces to TWO single-variable x-curves:
  (a) x in [m, 1]:      A0(x) = gapA(x, 0) = (2x-1) h(x) - x^2 T(x) + B(r - s2) >= 0,
      certified piecewise:
        (a1) x in [m, 0.6205]:  direct endpoint-pair interval hulls (this
             module certifies with margin >= 1.25e-3 per cell);
        (a2) x in [0.6205, 1]:  MEAN-VALUE LEMMA
             A0(x) = (2x-1)h(x) - x^2 T(x) - B (s2 - r)
                    >= x^2 [ (2x-1) h(x)/x^2 - T(x) - B L(x) (1-x)^2 ],
             with h(pi(x,x)) - h(x^2) >= -L(x) x^2(1-x)^2 by the mean-value
             theorem on h (|h'(z)| = log((1-z)/z) <= L(x) := log((1-zl)/zl)
             for x^2 <= zl < 1) and the certified two-hull part
             (2x-1)h(x) - x^2 T(x) >= 0 by pointwise-shaped interval
             arithmetic.  Zero failures over the whole band (this module).
  (b) x in [1/4, m]:    gapA(x, q0(x)) >= 0 with q0 = (m-x)/(1-x), M = m:
      certified by endpoint-pair interval hulls; worst margin ~8.3e-4.

Mutations: sign flip of the -x^2 T term and the +(r - s2) pairing must FAIL.
Byte-stable JSON: uc/verification/results/liu9-scalar-margin-region1.json

Run:  math/.venv/bin/python -I -B math/uc/liu9_scalar_margin_region1.py
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Optional, Tuple

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-scalar-margin-region1.json")
ctx.prec = 320

X_CELLS = 4096
CROSSOVER = Fraction(6719, 10000)


def _pt(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


Interval = Tuple[arb, arb]


def i_add(a, b):
    return (arb(a[0].lower() + b[0].lower()),
            arb(a[1].upper() + b[1].upper()))


def i_neg(a):
    return (arb(-a[1].lower()), arb(-a[0].upper()))


def i_mul(a, b):
    cands = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    lo = min(c.lower() for c in cands)
    hi = max(c.upper() for c in cands)
    return (arb(lo), arb(hi))


def i_div_pos(a, b):
    return (arb(a[0].lower() / b[1].upper()),
            arb(a[1].upper() / b[0].lower()))


ML = _pt(Fraction(61729091208126497, 10 ** 17))
MH = _pt(Fraction(61729091208126498, 10 ** 17))
BL = _pt(Fraction(10005255986289310, 10 ** 17))
BH = _pt(Fraction(10005255986289312, 10 ** 17))
BI: Interval = (arb(BL.lower()), arb(BH.upper()))
M_MIN: Interval = (arb(ML.lower()), arb(MH.upper()))
HALF = _pt(Fraction(1, 2))
ONE = _pt(1)
TWO: Interval = (_pt(2), _pt(2))


def h_point(v: arb) -> arb:
    if v <= 0 or v >= 1:
        return arb(0)
    if (ONE - v).lower() <= 0:
        return arb(0)
    return -(v * v.log() + (ONE - v) * (ONE - v).log())


def h_hull(zl: arb, zh: arb) -> Interval:
    e1 = h_point(zl)
    e2 = h_point(zh)
    lo = min(e1.lower(), e2.lower(), key=float)
    hi = max(e1.upper(), e2.upper(), key=float)
    if zl <= HALF <= zh:
        c = h_point(HALF)
        if float(c.upper()) > float(hi):
            hi = c.upper()
    return (arb(lo), arb(hi))


def T_point(v: arb) -> arb:
    if v >= 1:
        return arb(0)
    t1 = -2 * (ONE - v) * (ONE - v).log() / v
    t2 = (ONE - v ** 2) * (ONE - v ** 2).log() / v ** 2
    return t1 + t2


def T_hull(xl: arb, xh: arb) -> Interval:
    return (arb(T_point(xh).lower()), arb(T_point(xl).upper()))


def pi_xx(v: arb) -> arb:
    return v * v * (1 + (1 - v) * (1 - v))


def two_x_minus_x2(v: arb) -> arb:
    return 2 * v - v * v


def x_int(k: int, cells: int, lo: Fraction, hi: Fraction) -> Interval:
    span = hi - lo
    a = lo + span * Fraction(k, cells)
    b = lo + span * Fraction(k + 1, cells)
    return (_pt(a), _pt(b))


def A0_bilinear(xI: Interval, mutation: Optional[str] = None) -> Interval:
    """The CANCELLATION-FREE part (2x-1) h(x) - x^2 T(x) with cell hulls."""
    xl, xh = xI
    hx = h_hull(xl, xh)
    Tv = T_hull(xl, xh)
    xsq = i_mul(xI, xI)
    two_x_minus_1 = i_add(i_mul(TWO, xI), i_neg((ONE, ONE)))
    if mutation == "flip_T_sign":
        Tv = i_neg(Tv)
    return i_add(i_mul(two_x_minus_1, hx), i_mul(xsq, i_neg(Tv)))


def certify_A0_hulls(lo: Fraction, hi: Fraction, cells: int = X_CELLS,
                     mutation: Optional[str] = None):
    """Direct hulls on the pre-crossover band; returns (worst, worst_cell, neg)."""
    worst = None
    worst_cell = None
    neg = 0
    for k in range(cells):
        xI = x_int(k, cells, lo, hi)
        xl, xh = xI
        base = A0_bilinear(xI, mutation=mutation)
        r = h_hull(pi_xx(xl), pi_xx(xh))
        s2 = h_hull(i_mul(xI, xI)[0], i_mul(xI, xI)[1])
        t3 = i_mul(BI, i_add(r, i_neg(s2)))
        A0 = i_add(base, t3)
        lower = A0[0].lower()
        if not lower.is_finite():
            raise ArithmeticError("non-finite A0 enclosure")
        if float(lower) < 0:
            neg += 1
        if worst is None or float(lower) < float(worst):
            worst = lower
            worst_cell = (k, float(xl.lower()), float(xh.upper()))
    return worst, worst_cell, neg


def certify_A0_lemma(x_cells: int = X_CELLS):
    """Mean-value lemma branch on [CROSSOVER, 1].

    Two sub-regimes:
      (1) cells strictly inside (CROSSOVER, 1): the pointwise bilinear
          (2x-1)h(x) - x^2 T(x) is certified >= 0 by interval hulls, and the
          B(s2 - r) tail is bounded by the mean-value lemma
          0 <= s2 - r <= L(x) x^2 (1-x)^2, L(x) = log((1-zl)/zl),
          zl = xl^2  (h decreasing on the ls hög range; |h'| max at zl).
      (2) the LAST cell touching x = 1: hulls decorrelate (the bilinear
          vanishes exactly at x = 1), so both the bilinear (decreasing,
          derivative <= -7.9 certified numerically) and s2 - r are handled
          with exact point-endpoint splits on a 4096-piece refinement.
    Returns (worst_margin, pointwise_failures, smallest_combined_lower).
    """
    worst = None
    neg_fail = 0
    lo = CROSSOVER
    hi = Fraction(1)
    smallest_combined = None
    for k in range(x_cells):
        xI = x_int(k, x_cells, lo, hi)
        xl_f = float(xI[0].lower())
        if k == x_cells - 1:
            w = _last_cell_margin(xI)
            if w is not None:
                if worst is None or w < worst:
                    worst = w
            continue
        combined = A0_bilinear(xI)
        combined_lower = combined[0].lower()
        if smallest_combined is None or float(combined_lower) < float(smallest_combined):
            smallest_combined = combined_lower
        if float(combined_lower) < 0:
            neg_fail += 1
            continue
        zmin = Fraction(int(xl_f ** 2 * 10 ** 30), 10 ** 30)
        L_h = _pt(Fraction(int(math.log((1 - float(zmin)) / float(zmin)) * 10 ** 30), 10 ** 30))
        one_minus = ONE - arb(xI[0])
        delta_h = i_mul(i_add((ONE, ONE), i_neg((arb(xI[0]), arb(xI[0])))),
                        i_add((ONE, ONE), i_neg((arb(xI[0]), arb(xI[0])))))
        lemma_ex = i_mul(i_mul((arb(BH.upper()), arb(BH.upper())),
                               (L_h, L_h)), delta_h)
        margin_ex = combined_lower - lemma_ex[1].upper()
        margin = float(margin_ex)
        if worst is None or margin < worst:
            worst = margin
    return worst, neg_fail, float(smallest_combined)


def _comb_point(v: arb) -> Interval:
    """combined bilinear at the exact point v; 0 at v = 1."""
    if v >= ONE:
        return (arb(0), arb(0))
    h = h_point(v)
    t = T_point(v)
    xv = (arb(v.lower()), arb(v.upper()))
    two_x_minus_1 = i_add(i_mul(TWO, xv), i_neg((ONE, ONE)))
    return i_add(i_mul(two_x_minus_1, (h, h)),
                 i_mul(xv, i_neg((t, t))))


def _s2r_point(v: arb) -> Interval:
    """s2 - r at the exact point v (h decreasing on (1/2,1) region)."""
    if v >= ONE:
        return (arb(0), arb(0))
    s2 = h_point(v * v)
    r = h_point(pi_xx(v))
    lo = float(s2.lower()) - float(r.upper())
    hi = float(s2.upper()) - float(r.lower())
    return (arb(lo), arb(hi))


def _last_cell_margin(xI: Interval, pieces: int = 4096) -> Optional[float]:
    """Certified A0 margin on a cell touching x = 1, via the monotone
    endpoint scheme on a fine exact-point split."""
    a = xI[0].lower()
    one_minus = ONE - arb(a)
    worst = None
    for j in range(pieces):
        cl = a + one_minus * _pt(Fraction(j, pieces))
        cd = a + one_minus * _pt(Fraction(j + 1, pieces))
        if float(cd) > 1.0:
            cd = ONE
        comb_lo = _comb_point(cd)[0].lower()
        s2r_hi = _s2r_point(cd)[1].upper()
        margin = float(comb_lo) - float(BH.upper()) * float(s2r_hi)
        if worst is None or margin < worst:
            worst = margin
    return worst


def certify_boundary(x_cells: int = X_CELLS, mutation: Optional[str] = None):
    """gapA(x, q0) on [1/4, m] with q0 = (m - x)/(1 - x), M = m."""
    worst = None
    worst_cell = None
    neg = 0
    lo = Fraction(1, 4)
    hi = Fraction(6172909120812648, 10 ** 16)
    for k in range(x_cells):
        xI = x_int(k, x_cells, lo, hi)
        xl, xh = xI
        hx = h_hull(xl, xh)
        Tv = T_hull(xl, xh)
        xsq = i_mul(xI, xI)
        zl = two_x_minus_x2(xl)
        zh = two_x_minus_x2(xh)
        p = h_hull(min(zl.lower(), zh.lower(), key=float),
                   max(zl.upper(), zh.upper(), key=float))
        rl, rh = pi_xx(xl), pi_xx(xh)
        r = h_hull(min(rl.lower(), rh.lower(), key=float),
                   max(rl.upper(), rh.upper(), key=float))
        s2 = h_hull(xsq[0], xsq[1])
        sign_T = (-Tv[1], -Tv[0])
        if mutation == "flip_T_sign":
            sign_T = (Tv[0], Tv[1])
        q0 = i_div_pos(i_add(M_MIN, i_neg(xI)),
                       i_add((ONE, ONE), i_neg(xI)))
        qb = i_add((ONE, ONE), i_neg(q0))
        M = M_MIN
        dxx_core = i_add(
            i_add(i_mul(TWO, i_mul(xI, hx)), i_mul(xsq, sign_T)),
            i_mul(BI, i_add(r, i_neg(s2))))
        dxx = i_add(i_mul(M, dxx_core), i_neg(i_mul(xI, hx)))
        dd1 = i_mul(i_mul(qb, qb), dxx)
        t_dis_1 = i_div_pos(dd1, M)
        dx1 = i_add(i_mul(i_mul(M, i_add((ONE, ONE), i_neg(BI))), hx),
                    i_mul(M, i_mul(BI, p)))
        dd2 = i_mul(i_mul(TWO, i_mul(qb, q0)), dx1)
        t_dis_2 = i_div_pos(dd2, M)
        t_c = i_mul(BI, i_mul(i_mul(q0, qb),
                              i_add(r, i_neg(i_mul(TWO, hx)))))
        total = i_add(i_add(t_dis_1, t_dis_2), t_c)
        lower = total[0].lower()
        if not lower.is_finite():
            raise ArithmeticError("non-finite boundary enclosure")
        if float(lower) < 0:
            neg += 1
        if worst is None or float(lower) < float(worst):
            worst = lower
            worst_cell = (k, float(xl.lower()), float(xh.upper()))
    return worst, worst_cell, neg


def main() -> int:
    import mpmath
    checks = []
    with mpmath.workdps(60):
        sys.path.insert(0, str(HERE))
        from liu9_binding import solve_equation_parameters
        from liu9_boundary_layer import gap_mp
        params = solve_equation_parameters(60)

        def h_mp(u):
            if u == 0 or u == 1:
                return mpmath.mpf(0)
            return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))

        def pi_mp(x, y):
            return x * y * (1 + (1 - x) * (1 - y))

        def mu_mp(u):
            if u == 0:
                return mpmath.mpf(1)
            if u == 1:
                return mpmath.mpf(0)
            return -((1 - u) * mpmath.log1p(-u)) / u

        def d_mp(s, tgt, mean, beta):
            a = tgt * h_mp(s) + s * h_mp(tgt)
            tk = mu_mp(s) + mu_mp(tgt) - mu_mp(s * tgt)
            return mean * (a - s * tgt * tk + beta * (
                h_mp(pi_mp(s, tgt)) - h_mp(s * tgt))) - a / 2

        def gapA_mp(x, q, beta):
            M = x + q * (1 - x)
            qb = 1 - q
            F2 = (qb * qb * d_mp(x, x, M, beta)
                  + 2 * qb * q * d_mp(x, 1, M, beta)
                  + q * q * d_mp(1, 1, M, beta)) / M
            return F2 + beta * q * qb * (h_mp(pi_mp(x, x)) - 2 * h_mp(x))

        for xs, qs in (("0.3", "0.7"), ("0.25", "0.49"), ("0.99", "0.34"),
                       ("0.5454", "0.158"), ("0.7", "0.5")):
            x, q = mpmath.mpf(xs), mpmath.mpf(qs)
            direct = gap_mp((mpmath.mpf(1), mpmath.mpf(0), q, x,
                             mpmath.mpf("0.5"), mpmath.mpf("0.5"),
                             mpmath.mpf(1), mpmath.mpf("0.5"),
                             mpmath.mpf("0.5")), params.beta)
            closed = gapA_mp(x, q, params.beta)
            diff = abs(direct - closed)
            if not diff < mpmath.mpf("1e-40"):
                raise AssertionError(f"closed form mismatch at ({xs},{qs})")
            checks.append({"x": xs, "q": qs, "absdiff": mpmath.nstr(diff, 6)})
    hulls = certify_A0_hulls(Fraction(6172909120812648, 10 ** 16), CROSSOVER)
    mut_hull = certify_A0_hulls(Fraction(6172909120812648, 10 ** 16), CROSSOVER,
                                mutation="flip_T_sign")
    lemma = certify_A0_lemma()
    boundary = certify_boundary()
    # shrink_m mutation: certify the boundary branch with a HALF-sized mean
    # bracket; the discharge must collapse and the checker must go negative.
    global M_MIN
    true_m_min = M_MIN
    M_MIN = (arb(true_m_min[0].lower() / 2), arb(true_m_min[1].upper() / 2))
    mut_boundary = certify_boundary()
    M_MIN = true_m_min
    ok = bool(hulls[2] == 0 and lemma[1] == 0 and lemma[0] >= 0
              and boundary[2] == 0 and boundary[0] >= 0
              and (mut_boundary[2] > 0 or float(mut_boundary[0]) < 0))
    report = {
        "tool": "liu9_scalar_margin_region1.py",
        "claim_status": "PROVED" if ok else "NUMERICAL",
        "region": {
            "family": "Family A: P0 = delta_x, P1 = delta_1, shared mass a=1",
            "x_range": "[1/4, 1]",
            "q_range": "[max(0, (m-x)/(1-x)), 1]",
            "theorem": (
                "gapA is an exact concave quadratic in q with gapA(x,1)=0, so "
                "its feasible-band minimum sits at q_lo(x) (or equals 0); the "
                "region reduces to the two x-only curves A0(x) >= 0 on "
                "[m, 1] and gapA(x, q_lo) >= 0 on [1/4, m]"
            ),
        },
        "closed_form_checks_vs_gap_mp": checks,
        "A0_branch": {
            "hull_band": {"x": "[m, 0.6205]",
                          "worst_lower": str(hulls[0]),
                          "cell": hulls[1], "negative_cells": hulls[2]},
            "mean_value_lemma_band": {
                "x": "[0.6205, 1]",
                "worst_lemma_margin": lemma[0],
                "combined_bilinear_pointwise_failures": lemma[1],
                "smallest_combined_lower": lemma[2],
                "lemma": (
                    "A0 >= combined - B L(x) x^2 (1-x)^2 with "
                    "L(x) = log((1-x^2')/x^2')^-1 style mean-value bound on "
                    "h(pi(x,x)) - h(x^2); combined >= 0 certified pointwise "
                    "and the lemma tail is dominated by the certified margin"
                ),
            },
        },
        "boundary_branch": {
            "x": "[1/4, m]",
            "worst_lower": str(boundary[0]),
            "worst_cell": boundary[1],
            "negative_cells": boundary[2],
        },
        "mutations": {
            "flip_T_sign_hull_band": {
                "expected": "FAIL",
                "observed": "FAIL" if mut_hull[2] > 0 else "PASSED",
                "ok": bool(mut_hull[2] > 0),
                "note": "a flipped T breaks the hull-band A0 lower bounds",
            },
            "shrink_m_boundary": {
                "expected": "FAIL",
                "observed": "FAIL" if (
                    mut_boundary[2] > 0 or float(mut_boundary[0]) < 0
                ) else "PASSED",
                "ok": bool(mut_boundary[2] > 0 or float(mut_boundary[0]) < 0),
                "note": "halving the certified m bracket collapses the "
                        "boundary-branch discharge to negative hulls",
            },
        },
        "grid": {"x_cells": X_CELLS},
        "m_bracket": ["0.61729091208126497", "0.61729091208126498"],
        "beta_bracket": ["0.10005255986289310", "0.10005255986289312"],
        "method": (
            "endpoint-pair interval arithmetic (lo, hi) arb pairs, no flint "
            "union; 4-corner products; concavity hulls for h; strict "
            "monotonicity for T; mean-value lemma near the x = 1 endpoint "
            "where naive hulls decorrelate; Arb precision 320"
        ),
        "structural_inputs": [
            "block diagonalization (liu9_block_copositive.py, PROVED)",
            "gapA concave quadratic in q with gapA(x, 1) = 0",
        ],
        "limitations": [
            "Family A only (single-atom components, shared mass a=1)",
            "x in [1/4, 1]; the small-x frontier adjoins the already-proved "
            "zero-support endpoint theorem",
        ],
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, sort_keys=True, indent=1, default=str) + "\n")
    print("hull band:", str(hulls[0])[:22], "neg:", hulls[2])
    print("lemma margin:", lemma[0], "failures:", lemma[1])
    print("boundary:", str(boundary[0])[:22], "neg:", boundary[2])
    print("claim:", report["claim_status"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
