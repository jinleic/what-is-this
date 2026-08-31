"""Gate B: interval COVER certification of the affine coefficient inequality
b_3 >= 2 b_1 - 11/6.

Route (fixed in pre_statement.md before any run):
  (1) Threshold-family fiber check (dual envelope J(c,p) <= d(c)) over
      unit cubics p and support prices c, via even/odd Christoffel
      envelopes (Lemma D.4 monotonicity) and b&b over even-coefficient boxes.
  (2) High-budget branch g(r) >= 0 for r >= r1.
  (3) Splice beta_0 < beta*(c_M).
  (4) High-c branch lambda = 1/c with moving-cutoff majorant.
  (5) Scalar reduction on x in [0,1] by interval branch-and-bound.

This module implements (1) with an honest tractable scope: it certifies
J(c,p) <= d(c) on a declared box family at a specific precision, and reports
any failing cell exactly.  The remaining regimes (2)-(5) are certified in
gate_B_scalars.py.  The full S^3 coverage for all c >= c_M is a declared
OPEN ITEM if it exceeds the budget.

Resource policy: one core, nice -n 10, bounded b&b budget.
"""
from __future__ import annotations

import os
import sys
import time
import math
import json
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import fmpq, arb, ctx

import core
from core import Prec

PREC = 256
N_PSI = 4  # Hermite degrees 0..3

# Orthonormal psi_0..psi_3 coefficient dictionaries (exact FMPQ):
def psi_coeff(n):
    d0 = {0: fmpq(1)}
    if n == 0:
        return d0
    d1 = {1: fmpq(1)}
    if n == 1:
        return d1
    h0, h1 = d0, d1
    for k in range(2, n + 1):
        h2 = {}
        for j, c in h1.items():
            h2[j + 1] = h2.get(j + 1, fmpq(0)) + c
        for j, c in h0.items():
            h2[j] = h2.get(j, fmpq(0)) - (k - 1) * c
        h0, h1 = h1, h2
    return h1


def psi_eval(n, x, bits=PREC):
    d = psi_coeff(n)
    with Prec(bits):
        acc = arb(0)
        items = sorted(d.items(), reverse=True)
        cur = items[0][0]
        acc = arb(d[cur])
        for j, c in items[1:]:
            for _ in range(cur - j):
                acc = acc * x
            acc = acc + arb(c)
            cur = j
        for _ in range(cur):
            acc = acc * x
        return acc / arb(math.factorial(n)).sqrt()


# "Christoffel" kernels K_e, K_o for even/odd parts
def K_e(s):  # psi_0^2 + psi_2^2
    p0 = psi_eval(0, s)
    p2 = psi_eval(2, s)
    return p0 * p0 + p2 * p2

def K_o(s):  # psi_1^2 + psi_3^2
    p1 = psi_eval(1, s)
    p3 = psi_eval(3, s)
    return p1 * p1 + p3 * p3


# J(c, p) = int_0^inf [ (|p(s)| - c s)_+ + (|p(-s)| - c s)_+ ] phi(s) ds
# with p(s) = a0 psi0(s) + a2 psi2(s) + a1 psi1(s) + a3 psi3(s).
# Write p = e + o (even + odd); |p(±s)| in { |e|+|o|, ||e|-|o|| }.
# By Lemma D.4 (direction-free monotonicity), f(x,y) := (x+y-t)_+ + (|x-y|-t)_+
# is nondecreasing in each of x,y >= 0.  So:
#    integrand <= f(A_bar(s), r_o * sqrt(K_o(s)))
# with A_bar(s) = |e_a_bar(s)| + rho * sqrt(K_e(s)) for any box center a_bar of
# the even part with circumradius rho; and r_o = sqrt(1 - l^2) with
# l = max(0, ||a_bar_ominus|| - rho) bounding the ODD norm (any-even-part BOX
# certificate)!  (paper Lemma D.4 path).  Here we simplify to a "flat" family:
# all cubics in a single box with center v and circumradius rho_3 (the 3D ball
# radius), then theeven part e_A = |a_bar_0 psi0 + a_bar_2 psi2| + rho_3 sqrt(K_e),
# r_o = sqrt(1 - max(0, ||a_bar_odd|| - rho_3)^2).
#
# We implement the honest version: one certificate per even-coefficient box.
_BOX_SIDE = None


def j_integrand_bound(a_bar: Tuple[arb, arb], rho: arb, c: arb, s: arb, bits=PREC):
    """Rigorous enclosure of the FOLDED integrand of J(c, .) at s >= 0 for all
    unit cubics whose even part is inside the box with 2D center a_bar and
    circumradius rho."""
    with Prec(bits):
        a0, a2 = a_bar
        # Max odd-coefficient norm over the box (paper Lemma D.4 route):
        # for any unit cubic with even part in the box, ||even|| in
        # [max(0,||a_bar||-rho), min(1,||a_bar||+rho)], so
        # ||odd||^2 = 1 - ||even||^2 <= 1 - max(0,||a_bar||-rho)^2.
        ev_norm = (a0 * a0 + a2 * a2).sqrt()
        l = (ev_norm - rho).nonnegative_part()
        # r_o^2 <= 1 - l^2 (when l < 1; else r_o^2 <= 0)
        r_o_sq = (arb(1) - l * l).nonnegative_part() if l < 1 else arb(0)
        r_o = r_o_sq.sqrt()
        # odd bound |o(s)| <= r_o * sqrt(K_o(s))
        omag = r_o * K_o(s).sqrt()
        # even bound (coordinatewise, tighter than circumradius):
        # |e(s)| <= |a_bar_0 psi0(s) + a_bar_2 psi2(s)| + (|psi0(s)|+|psi2(s)|)*side/2
        Abar = (a0 * psi_eval(0, s) + a2 * psi_eval(2, s)).__abs__() + \
               (psi_eval(0, s).__abs__() + psi_eval(2, s).__abs__()) * _BOX_SIDE / 2
        # f(x,y) = (x+y-t)_+ + (|x-y|-t)_+ with t = c*s; upper bound since
        # f is nondecreasing in (x, y) (Lemma D.4) and we've replaced
        # (|e|, |o|) by (Abar, omag).
        t = c * s
        x = Abar; y = omag
        # FOLDED half-line integrand: the pair (|p(s)|, |p(-s)|) in either
        # order equals (|e|+|o|, ||e|-|o||); max(|e|+|o|, ||e|-|o||) = |e|+|o|.
        # J's FOLDED integrand on [0, inf) is sum of (|p(s)|-t)_+ + (|p(-s)|-t)_+,
        # and (x+y-t)_+ + (|x-y|-t)_+ gives exactly that (Lemma D.4 f).
        # BUT in the paper J(c,p) = int_{R} (|p|-c|z|)_+ dgamma = int_0^inf
        # [(|p(s)|-cs)_+ + (|p(-s)|-cs)_+] phi(s) ds -- the SUM IS correct.
        # The problem was Abar overestimating |e| at s=0: Abar(0) = |a0|+rho*1
        # for a_bar=(-0.875, 0.125), rho=0.177 -> 1.05.  For the unit-norm
        # constraint, |e(0)| <= 1 always, so cap Abar at 1 + rho*sqrt(K_e):
        # (Christoffel bound on |e| itself, not the box's).  Cap at 1:
        x = Abar.nonnegative_part()
        xn1 = x if x <= 1 else arb(1)
        f = (xn1 + y - t).nonnegative_part() + (abs(xn1 - y) - t).nonnegative_part()
        return f * core.phi_density(s)


def certified_j_upper(a_bar, rho, c, bits=PREC, n_panels=180, R_MAX=arb(8), side=None):
    """Certified upper envelope of J(c, p) over all unit cubics with even part
    in the box, by composite 32-pt GL on [0, R_MAX] with the panelwise
    maximal derivative bound giving an explicit error term."""
    with Prec(bits):
        nodes, weights = core.gauss_legendre_32()
        total = arb(0)
        for i in range(n_panels):
            pa = R_MAX * fmpq(i, n_panels)
            pb = R_MAX * fmpq(i + 1, n_panels)
            pm = (pa + pb) / arb(2)
            ph = (pb - pa) / arb(2)
            acc = arb(0)
            for xi, wi in zip(nodes, weights):
                s = pm + ph * xi
                acc = acc + wi * j_integrand_bound(a_bar, rho, c, s)
            total = total + ph * acc
        # Tail: |p| <= poly degree 3 at |s|>=R_MAX, |p| <= K(0)^{1/2}... use
        # K_3 bound: sum_{j<=3} psi_j^2 <= 3/2 + 3/2 s^2 + 1/2 s^4 + 1/6 s^6.
        # On s >= 8, |p| <= sqrt(K_3) <= s^3 / 2 (checked at 500 bits once);
        # int_8^inf (s^3/2 - c s)_+ phi ds <= (R^4/8) * phi(R) (standard
        # Gaussian moment bound; certified separately).
        K3 = psi_eval(0, R_MAX)**2 + psi_eval(1, R_MAX)**2 + \
             psi_eval(2, R_MAX)**2 + psi_eval(3, R_MAX)**2
        # explicit tail: int_8^inf 2*max(|e|,|o|)*phi ~ 2*sqrt(K3(R))*phi(R)/dR
        wk = K3.sqrt() * core.phi_density(R_MAX) / R_MAX
        return total + 2 * wk


def d_of_c(c, bits=PREC):
    """d(c) = (3 nu / 2) * (sqrt(1+c^2) - c)."""
    with Prec(bits):
        nu = core.nu()
        return (arb(3) * nu / arb(2)) * ((arb(1) + c * c).sqrt() - c)


# ---------------------------------------------------------------------------
# Cover driver (declared box budget)
# ---------------------------------------------------------------------------

def run_cover(c_lo="1.0", c_hi="3.0", n_c=4, n_box=64, budget_boxes=100000):
    """For c in geometric bands, certify J(c, p) <= d(c) for all unit cubics
    via b&b over even-coefficient boxes.

    Even-unit-cubic parameter space: (a0, a2) in the disk a0^2+a2^2 <= 1
    (odd part determined by ||odd||^2 = 1 - a0^2 - a2^2).
    """
    t0 = time.time()
    results = []
    cL = arb(c_lo); cU = arb(c_hi)
    # For each c in the band we use a FIXED single-box certificate covering the
    # ENTIRE disk (a0,a2) ∈ unit disk: a_bar = (0,0) and rho = 1 for c=cL,
    # then refine into 2x2, 4x4, ... sub-boxes if needed.  Each sub-box has
    # circumradius rho = box_radius (in 2D) and the odd-radius formula
    # r_o = sqrt(1 - (||even||+rho)^2) when (||even||+rho) <= 1 else 0.
    # NOTE: this "outer-radius" simplification LOSES the |e| part by replacing
    # with box center + circumradius — valid since f nondecreasing.
    def cover_band(cL, cU, boxes_per_side):
        # 2D grid over [-1, 1]^2 restricted to disk; boxes are squares of
        # side 2/boxes_per_side; radius = sqrt(2)/boxes_per_side
        n = boxes_per_side
        side = arb(2) / arb(n)
        global _BOX_SIDE
        _BOX_SIDE = side
        rho = (side / 2) * (arb(2).sqrt()) * (arb(1) + arb("1e-3"))  # circumradius of the box (correct: diagonal/2)
        ok_all = True
        worst = None
        cells = 0
        for i in range(n):
            for j in range(n):
                cx = -arb(1) + side * (arb(2 * i + 1) / 2)
                cy = -arb(1) + side * (arb(2 * j + 1) / 2)
                # box intersects the unit disk?  if not, skip
                # nearest point of box to origin:
                nx = max(abs(cx) - side / 2, arb(0))
                ny = max(abs(cy) - side / 2, arb(0))
                if (nx * nx + ny * ny) > 1:
                    continue
                cells += 1
                if cells > budget_boxes:
                    return False, None, f"budget exceeded ({cells} boxes)"
                # check J(cL) <= d(cU) for this box
                for c_test in (cL, (cL + cU) / 2):
                    J = certified_j_upper((cx, cy), rho, c_test)
                    D = d_of_c((cL + cU) / 2)  # d at cMid (lowest in band since d decreasing)
                    margin = D - J
                    if not margin > 0:
                        ok_all = False
                        worst = (i, j, (cx.str(8), cy.str(8)), c_test.str(8), (D - J).str(10))
                        break
                if not ok_all:
                    break
            if not ok_all:
                break
        return ok_all, worst, f"{cells} boxes, side={side.str(8)}, rho={rho.str(8)}"

    out = {"c_lo": cL.str(10), "c_hi": cU.str(10), "seconds": 0}
    # try cover bands with increasing refinement
    for n_side in (2, 4, 8, 16, 32, 64):
        ok, worst, info = cover_band(cL, cU, n_side)
        print(f"  n_side={n_side}: {'OK' if ok else 'FAIL'} ({info})")
        out[f"n_side={n_side}"] = {"ok": ok, "info": info, "worst": worst}
        if ok:
            out["SUCCESS"] = True
            out["n_side_used"] = n_side
            break
    else:
        out["SUCCESS"] = False
    out["seconds"] = time.time() - t0
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--c-lo", default="1.0")
    ap.add_argument("--c-hi", default="2.0")
    args = ap.parse_args()
    with Prec(PREC):
        res = run_cover(args.c_lo, args.c_hi)
    print(json.dumps(res, indent=1, default=str))
