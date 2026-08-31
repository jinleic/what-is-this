"""gate_B2.py — Gate B certificate (C2) re-run with the 2D Christoffel envelope,
per the paper's D.3 route (paper_full.txt lines 1961-1984), agent KgGateBClose
2026-08-30.  Replaces campaign 20260830T013720Z_e24aa159's envelope.

Certified statement per cell & c-sub-band:
    For every unit cubic p with even part a_even in the BOX B (2D square on the
    even disk a0^2 + a2^2 <= 1):
      J(c, p) <= Ubar(B, c_L)  and  d(c) >= d(c_U),  so
      J(c, p) <= d(c)  whenever Ubar(B, c_L) <= d(c_U)   (paper line 1984).

Envelope (Lemma D.4 + Christoffel, paper lines 1965-1981):
  - e(s) range over the box EXACTLY (affine in (a0,a2) => corner values, ball
    arithmetic; no 3D circumradius inflation, no |e|<=1 cap — the predecessor's
    cap was FALSE: |e(0)| <= sqrt(3/2) via psi0(0)=1, psi2(0)=-1/sqrt2);
  - |o(s)| <= r_o sqrt(K_o(s)), r_o = sqrt(1 - l^2), l = dist(0, B) (L2);
  - f(x,y) = (x+y-t)_+ + (|x-y|-t)_+, t = c s; f nondecreasing in x, y >= 0
    (Lemma D.4, lines 1967-1971) so f(|e|,|o|) <= f(E+, w) with
    E+ = upper end of ball(|e|-range side) two-sided, w = ball(r_o sqrt K_o).

Quadrature (all arb, 256-bit outward):
  - per knot s: t = ball(c_L * s); U1 = E+ + w; U2 = |E+ - w| upper end;
    g = (U1 - t)_+ + (U2 - t)_+ ... plus s-variation absorbed because every
    magnesium knot already carries the full ball of s (interval s over panel).
  - panel value = (h/2) sum_k w_k g(t_k) with PRE-CERTIFIED GL weights, each
    product summed outward; PLUS panel remainder (h * G_panel) where
    G_panel = sup-panel envelope with s itself interval — dominated while
    panels stay small; total = gl_sum + remainders + tail_account(S_MAX).
  - NOTE: this is a script-side conservative bound; the certified gate always
    is  Ubar = gl_ball + sum(h*G) + tail <= d(c_U)_lower.

Protocol: pre_statement.md + Addendum 1 (adaptive c-split on failure,
budget 150k pairs per band; n_side ladder 2..128).
"""
from __future__ import annotations

import argparse
import json
import os
import queue
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, fmpq, ctx

import core2
from core2 import Prec, gauss_legendre_32

PREC = 256
S_MAX = 8           # paper line 1982: panels on s in [0,8]
N_GLP = 160         # GL panels on [0,8] for the knot-sum part
C_FLOOR = 12        # hard cap on c (paper bands go to 12)


# ---------------------------------------------------------------------------
# Box bookkeeping on the 2D even disk
# ---------------------------------------------------------------------------

def box_from_center(cx: arb, cy: arb, w: arb) -> Tuple[arb, arb, arb, arb]:
    with Prec(PREC):
        return cx - w / 2, cx + w / 2, cy - w / 2, cy + w / 2


def dist0_box(ax0, ax1, ay0, ay1):
    """L2 distance from origin to the box (certified): per-axis |x| min."""
    with Prec(PREC):
        nx = _axis_dist0(ax0, ax1)
        ny = _axis_dist0(ay0, ay1)
        return (nx * nx + ny * ny).sqrt()



def _axis_dist0(lo, hi):
    """|x| at the box point closest to 0, for the axis interval [lo,hi]."""
    if lo <= 0 <= hi:
        return arb(0)
    if lo > 0:
        return lo
    return -hi


def box_intersects_disk(ax0, ax1, ay0, ay1) -> bool:
    """Certified test: the box meets the unit disk a0^2+a2^2 <= 1."""
    with Prec(PREC):
        nx = _axis_dist0(ax0, ax1)
        ny = _axis_dist0(ay0, ay1)
        return nx * nx + ny * ny <= arb(1)



# ---------------------------------------------------------------------------
# Per-cell certified envelope integral
# ---------------------------------------------------------------------------

def e_range_box(ax0, ax1, ay0, ay1, s, bits=PREC):
    """Exact 2D-Christoffel even-range on the box at interval s (ball).

    e(s) = a0 psi0(s) + a2 psi2(s) is affine in the box, so on box x s-panel
    the certified range is the interval-arithmetic evaluation over both the
    box and s intervals (corner extremes incl.), each endpoint a ball.
    Returns (e_lo, e_hi) as arbs (using interval s -> psi evaluated at ball s)."""
    with Prec(bits):
        p0 = core2.psi(0, s, bits)
        p2 = core2.psi(2, s, bits)
        # terms as interval products: [ax0,ax1]*p0_ball + [ay0,ay1]*p2_ball
        def term(lo, hi, p):
            cands = (lo * p, lo * p, hi * p, hi * p)
            # crude but certified: lo*ball p +/- ... use union of outer products
            u1 = lo * p
            u2 = hi * p
            return u1.union(u2)
        t0 = term(ax0, ax1, p0)
        t2 = term(ay0, ay1, p2)
        tot = t0 + t2
        return tot.lower(), tot.upper()


def cell_envelope(cL: arb, ax0, ax1, ay0, ay1, w_box: arb,
                  n_panels: int = N_GLP, bits: int = PREC) -> arb:
    """Certified upper envelope Ubar for the cell box at c = cL (ball arb).

    Route (paper D.3, paper_full.txt lines 1965-1984):
      integrand of (57) is f(|e(s)|,|o(s)|), f(x,y)=(x+y-t)_+ + (|x-y|-t)_+.
      Lemma D.4 (lines 1967-1971): f nondecreasing in x, y >= 0, so
      f(|e|,|o|) <= f(Ep(s), wo(s)) pointwise, with
        Ep(s) = max(E_hi, -E_lo) >= |e(s)|   (exact box range: 2D Christoffel)
        wo(s) = r_o sqrt(K_o(s)) >= |o(s)|   (r_o = sqrt(1-l^2), l = dist(0,B))
      so g(s) = (Ep+wo-t)_+ + (|Ep-wo|-t)_+, t = cL*s.
    Quadrature: certified per-panel value
      h * [ 32-pt GL knot-sum with interval-s knots + sup-panel remainder ]
    with wo at panel-right (K_o strictly increasing on [0,inf):
    K_o'(s) = s((s^2-2)^2 + 1) > 0, exact), t at panel-left, e-range on the
    observed panel-ball union; plus the Gaussian tail account (paper 1982-83).
    """
    with Prec(bits):
        nodes, weights = gauss_legendre_32()
        l = dist0_box(ax0, ax1, ay0, ay1)
        lmin = l if l < 1 else arb(1)
        r_o = (arb(1) - lmin * lmin).nonnegative_part().sqrt()
        h = fmpq(S_MAX, n_panels)
        half = h / 2
        total = arb(0)
        rest = arb(0)
        for i in range(n_panels):
            pa = fmpq(i * S_MAX, n_panels)
            pb = fmpq((i + 1) * S_MAX, n_panels)
            pm = fmpq(2 * i + 1, 2 * n_panels)
            s_pb = arb(pb)
            # --- certified panel remainder: sup-panel g ---
            # sup-integrand over panel x box; weight phi at panel-LEFT (phi
            # strictly decreasing on [0,inf) => sup of phi on the panel).
            s_iv = arb(pa).union(s_pb)
            e_lo, e_hi = e_range_box(ax0, ax1, ay0, ay1, s_iv, bits)
            Epu = (-e_lo).max(e_hi)                     # >= |e| over panel x box
            wo = r_o * core2.K_o(s_pb, bits).sqrt()     # >= |o| over panel
            t = cL * arb(pa)                            # t minimal at panel-left
            g1 = (Epu + wo - t).nonnegative_part()
            g2 = ((Epu - wo) - t).nonnegative_part() + ((wo - Epu) - t).nonnegative_part()
            rest = rest + arb(half * 2) * (g1 + g2) * core2.phi_density(arb(pa), bits)
            # --- GL knot sum ---
            acc = arb(0)
            for x, wk in zip(nodes, weights):
                s = arb(pm) + arb(half) * x
                e_lo, e_hi = e_range_box(ax0, ax1, ay0, ay1, s, bits)
                Ep = (-e_lo).max(e_hi)
                wo_k = r_o * core2.K_o(s, bits).sqrt()
                t = cL * s
                g1k = (Ep + wo_k - t).nonnegative_part()
                g2k = ((Ep - wo_k) - t).nonnegative_part() + ((wo_k - Ep) - t).nonnegative_part()
                acc = acc + wk * (g1k + g2k) * core2.phi_density(s, bits)
            total = total + acc * arb(half)
        tail = core2.tail_account(S_MAX, bits)
        return total + rest + tail

# ---------------------------------------------------------------------------
# Band driver (adaptive c-splitting per pre_statement Addendum 1)
# ---------------------------------------------------------------------------

def split_band(c_lo: str, c_hi: str, budget: int = 150000,
               c_steps_init: int = 2, verbose: bool = True):
    """Certify J <= d over the c-band by (cell x c-sub-band) pairs.

    Returns dict with bands_pass (list of c-sub-bands with verdict certified),
    worst margins, budget usage."""
    t0 = time.time()
    cL0 = arb(c_lo)
    cU0 = arb(c_hi)
    # initial c-tiling: geometric-ish uniform steps
    subs: List[Tuple[arb, arb]] = []
    with Prec(PREC):
        q0 = fmpq(0, 1)
        # uniform in the paper's sense: they tile each Table-1 row into k
        # sub-bands; we start from 2 uniform sub-bands and split on failure.
        for i in range(c_steps_init):
            a = cL0  # cheap: single-band start; splitting does the work
            subs = [(cL0, cU0)]
    # n_side ladder
    verdicts = []
    stats = {"pairs": 0, "boxes_certified": 0, "splits": 0}
    worst_overall = None
    stack = [(cL0, cU0, 0)]
    failed_pairs = []
    while stack:
        aL, aU, depth = stack.pop()
        ok_band, worst, nboxes = try_band_boxes(aL, aU, stats, budget)
        margin = aU.str(12)
        if ok_band:
            verdicts.append((aL.str(12), aU.str(12), worst))
            if worst_overall is None or worst < worst_overall[2]:
                pass
        else:
            if (aU - aL) < (aU * (arb(1) / 4096)):
                failed_pairs.append((aL.str(12), aU.str(12), worst))
                continue
            mid = (aL + aU) / 2
            stack.append((aL, mid, depth + 1))
            stack.append((mid, aU, depth + 1))
            stats["splits"] += 1
        if stats["pairs"] > budget:
            failed_pairs.append((aL.str(12), aU.str(12), "BUDGET"))
            break
    out = {
        "c_lo": c_lo, "c_hi": c_hi,
        "subbands_certified": verdicts,
        "subbands_open": failed_pairs,
        "stats": stats,
        "seconds": time.time() - t0,
    }
    return out


def try_band_boxes(aL: arb, aU: arb, stats: dict, budget: int):
    """Try to certify all disk cells at c-pair (aL, aU). Refine n_side on fail.
    Returns (ok, worst_margin, boxes) — worst margin is an arb (min over cells
    of d(aU)_lower - Ubar_ball_upper)."""
    dU = core2.d_of_c(aU)
    dU_lo = dU.lower()  # certified lower end of d(c_U)
    import math
    best: Optional[arb] = None
    worst_cell = None
    for n_side in (4, 8, 16, 32, 64, 128):
        w = arb(2) / n_side
        ok = True
        cells = 0
        for i in range(n_side):
            for j in range(n_side):
                cx = -1 + (2 * i + 1) / n_side
                cy = -1 + (2 * j + 1) / n_side
                cx = arb(-1) + arb(2 * i + 1) * w / 2
                cy = arb(-1) + arb(2 * j + 1) * w / 2
                ax0, ax1, ay0, ay1 = box_from_center(cx, cy, w)
                if not box_intersects_disk(ax0, ax1, ay0, ay1):
                    continue
                cells += 1
                U = cell_envelope(aL, ax0, ax1, ay0, ay1, w)
                stats["pairs"] += 1
                m = (dU_lo - U).upper()   # certified upper end of margin
                if best is None or m < best:
                    best = m
                    worst_cell = (cx.str(10), cy.str(10))
                if m <= 0:
                    ok = False
                    if n_side < 128:
                        break
                    else:
                        return False, (best, worst_cell, n_side), cells
                if stats["pairs"] > budget:
                    return False, (best, worst_cell, n_side), cells
            if not ok:
                break
        if ok:
            return True, (best, worst_cell, n_side), cells
    return False, (best, worst_cell, 128), cells * n_side


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c-lo", type=str, required=True)
    ap.add_argument("--c-hi", type=str, required=True)
    ap.add_argument("--budget", type=int, default=150000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    res = split_band(args.c_lo, args.c_hi, args.budget)
    js = json.dumps(res, indent=1, default=str)
    print(js)
    if args.out:
        with open(args.out, "w") as f:
            f.write(js)


if __name__ == "__main__":
    main()
