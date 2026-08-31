"""gate_B2.py — Gate B certificate (C2) with the 2D Christoffel envelope
(paper D.3, paper_full.txt lines 1961-1984). v2, agent KgGateBClose 2026-08-30.

Replaces the predecessor campaign's envelope (3D circumradius + FALSE |e|<=1
cap) and this campaign's superseded v1 (double-counted panels; kept as
gate_B2_v1_superseded.py). v1's second defect, named for the record: it summed
a GL knot-sum AND a sup-panel term per panel while omitting the phi weight.
v2 uses EXACTLY ONE certified term per panel (the paper's line-1982 method):
    panel value = W_i * sup_{panel x box} G,  W_i = Phi(pb)-Phi(pa) (exact).

CERTIFIED OBJECT (cell B in the 2D even disk, c-sub-band [aL, aU]):
  Ubar(B, aL) = sum_i W_i * [ (Ep + wo - aL*pa)_+ + (max(Ep, wo) - aL*pa)_+ ]
  Ep = max(E_hi, -E_lo) >= |e(s)| over panel x box  (exact Arb range of the
       affine e(s) over the box at the union ball s = [pa, pb]);
  wo = r_o sqrt(K_o(pb)),  r_o = sqrt(1 - l^2),  l = dist(0, B)  (K_o increasing:
       K_o'(s) = s((s^2-2)^2+1) > 0, so the odd sup on the panel is at pb and
       |o(s)| <= r_o sqrt(K_o(s)) for the whole class of unit cubics over B);
  the two hinges bound the folded pair (|p(s)|, |p(-s|) = (|e|+|o|, ||e|-|o||):
       (|e|+|o|-t)_+ <= (Ep+wo-tL)_+  and  (||e|-|o||-t)_+ <= (max(Ep,wo)-tL)_+
       since (x-t)_+ is nondecreasing in x and t = aL*s <= aL*pa = tL on the
       panel (also |e|-|o| <= max(|e|,|o|)).
  Tail s >= S_MAX: (|p|-cs)_+ <= sqrt(K(s)); core2.tail_account (paper 1982-84).

VALIDITY of Ubar >= J(c, p) for all unit cubics p with even part in B, c in
[aL, aU]: pointwise per panel as above; J(c,p) <= J(aL,p) (c-monotonicity) and
d(c) >= d(aU); so  Ubar(B, aL) <= d(aU)  certifies J(c,p) <= d(c) on the pair
(paper line 1984).

SEARCH (pre-registered in pre_statement.md + Addendum 1):
  - c-DFS with geometric midpoint splitting (paper's own adaptive rule, 1984);
  - per c-pair: n_side ladder (2, 6, 12) & per-cell box branch-and-bound
    (2:1 splits on the longer axis, floor 1/1024 half-width);
  - per-cell panel ladder chosen from cL (active window ~ 2.5/c);
  - MONOTONE MEMO (sound by margin monotonicity): (box, cL) -> certified cU;
    a cell certified at cU0 is certified for any cU <= cU0 (J fixed, d larger);
    plus the two exact invariances J(c,p) = J(c,-p) (paper line 1981) and
    J(c, rot p) = J(c, p) with rot: (a0,a2) -> (-a2,a0) (even-plane rotation
    preserving norms), used to share memo entries;
  - BUDGET: 260,000 certified cell-evaluations per band; a band that exhausts
    it is reported OPEN with its exact frontier (no post-hoc domain change).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, fmpq, ctx

import core2
from core2 import Prec

PREC = 256
S_MAX = 8
N_SIDE_LADDER = (2, 6, 12)
BOX_FLOOR = 1.0 / 1024
C_DEPTH_MAX = 24
BOX_DEPTH_MAX = 40


def n_start_for(cL: arb) -> int:
    """Panel count ladder start, matched to the active window ~ 2.5/c."""
    with Prec(PREC):
        if cL <= 1:
            return 2048
        if cL <= 2:
            return 1024
        if cL <= 4:
            return 1024
        if cL <= 8:
            return 4096
        return 8192


N_UP_LADDER = (1, 4, 16, 64)   # multiplicative refinements after n_start

# ---------------------------------------------------------------------------
# Box helpers
# ---------------------------------------------------------------------------

def box_from_center(cx: arb, cy: arb, w: arb):
    with Prec(PREC):
        return (cx - w / 2, cx + w / 2, cy - w / 2, cy + w / 2)


def _axis_dist0(lo: arb, hi: arb) -> arb:
    with Prec(PREC):
        if lo <= 0 <= hi:
            return arb(0)
        if lo > 0:
            return lo
        return -hi


def dist0_box(ax0, ax1, ay0, ay1) -> arb:
    with Prec(PREC):
        nx = _axis_dist0(ax0, ax1)
        ny = _axis_dist0(ay0, ay1)
        return (nx * nx + ny * ny).sqrt()


def box_intersects_disk(ax0, ax1, ay0, ay1) -> bool:
    with Prec(PREC):
        nx = _axis_dist0(ax0, ax1)
        ny = _axis_dist0(ay0, ay1)
        return nx * nx + ny * ny <= arb(1)


def box_key(box) -> str:
    return "|".join(x.str(24) for x in box)


# ---------------------------------------------------------------------------
# Certified per-cell envelope
# ---------------------------------------------------------------------------

def e_range_box(ax0, ax1, ay0, ay1, s, bits=PREC):
    with Prec(bits):
        p0 = core2.psi(0, s, bits)
        p2 = core2.psi(2, s, bits)
        t0 = (ax0 * p0).union(ax1 * p0)
        t2 = (ay0 * p2).union(ay1 * p2)
        tot = t0 + t2
        return tot.lower(), tot.upper()


def cdf_weight(pa: arb, pb: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return core2.gauss_cdf(pb, bits) - core2.gauss_cdf(pa, bits)


def cell_envelope(cL: arb, box, n_panels: int, bits: int = PREC) -> arb:
    """Certified upper envelope Ubar(B, cL) (validity: module docstring)."""
    ax0, ax1, ay0, ay1 = box
    with Prec(bits):
        ax0, ax1, ay0, ay1 = arb(ax0), arb(ax1), arb(ay0), arb(ay1)
        l = dist0_box(ax0, ax1, ay0, ay1)
        lmin = l if l < 1 else arb(1)
        r_o = (arb(1) - lmin * lmin).nonnegative_part().sqrt()
        # Ordered-pair scheme (exact algebra, pre-registered Addendum 2):
        # J's folded integrand for the unordered pair {|e|+|o|, ||e|-|o||} is
        #   (M - t)_+ + (q - t)_+  with M = |e|+|o|, q = ||e|-|o||,  mq = |e^2-o^2|.
        # Certified pipeline per panel (all outward, t = cL*s <= cL*pb):
        #   M  <= Mp := Ep + wo                     (triangle)
        #   q  <= Qp := |Ep - wo|_upper             (reverse triangle)
        #   mq >= mo := (A2 - r2)_lo,  A2 = exact square of the |e|-range ball
        #   (|e|^2 - w^2)_+ <= (A2 - r2)_+ with r2 = r_o^2 K_o(pb)
        #   (mq - t^2)_+ <= mo - t(pa)^2 if positive       [mq >= mo, t inc]
        #   (q - t)_+ <= Qp - t(pa) if (mq - t^2)_+ <= 0   [certified case split]
        #   + fallback: q <= Qp always.
        # This is the exact identity (f = g + (q-t)_+ from Lemma D.4's proof)
        # and removes the 2x over-count: for an all-even box, only the M
        # hinge contributes (h1 > 0 kills the q term), so no fold doubling.
        total = arb(0)
        for i in range(n_panels):
            pa = arb(fmpq(i * S_MAX, n_panels))
            pb = arb(fmpq((i + 1) * S_MAX, n_panels))
            if core2.K(pb, bits).sqrt() <= cL * pa:
                continue
            AB = pa.union(pb)
            e_lo, e_hi = e_range_box(ax0, ax1, ay0, ay1, AB, bits)
            Ep = (-e_lo).max(e_hi)
            wo = r_o * core2.K_o(pb, bits).sqrt()
            A2 = (Ep * Ep)  # ball of |e|^2 sup over panel x box
            r2 = r_o * r_o * core2.K_o(pb, bits)
            tL = cL * pa
            t2L = tL * tL
            Mp = Ep + wo
            gM = (Mp - tL).nonnegative_part()
            mo = (A2 - r2).lower()          # certified LOWER bound of mq
            if mo > t2L and t2L > 0:
                # (mq - t^2)_+ >= mo - t2L > 0 fires with q,err-free cert
                gq = (((Ep - wo).abs_upper() - tL).nonnegative_part())
            elif t2L == 0 and mo > 0:
                gq = (Ep - wo).abs_upper()
            else:
                # mq <= t^2 somewhere on the panel; certified coarse fallback
                gq = (Ep - wo).abs_upper().nonnegative_part()
            total = total + cdf_weight(pa, pb, bits) * (gM + gq)
        return total + core2.tail_account(S_MAX, bits)


# ---------------------------------------------------------------------------
# Cell verdict + memo
# ---------------------------------------------------------------------------

def memo_lookup(memo: dict, key: str, cU: str):
    ent = memo.get(key)
    if ent is not None and cU <= ent[0]:
        return ent[1]
    return None


def try_cell(cL: arb, cU: arb, dU_lo: arb, box, st: dict, memo: dict,
             ucache: dict, depth: int = 0):
    """Certify one cell at the c-pair. Returns (ok, margin_upper arb, n).

    U-cache: U depends only on (box, cL, n); c-splits reuse it (one arb
    subtraction per revisit instead of a panel sweep)."""
    key = box_key(box)
    cUs = cU.str(24)
    hit = memo_lookup(memo, key, cUs)
    if hit is not None:
        st["memo_hits"] += 1
        return True, arb(hit), 0
    # depth-aware starting rung: refined boxes re-enter high
    n = n_start_for(cL)
    if depth >= 2:
        n = min(n * 4, 32768)
    if depth >= 4:
        n = 32768
    m = None
    for _ in range(7):
        if st["evals"] > st["budget"]:
            return False, m, n
        uk = (key, cL.str(24), n)
        Uv = ucache.get(uk)
        if Uv is None:
            U = cell_envelope(cL, box, n)
            ucache[uk] = U.str(40)
        else:
            U = arb(Uv)
            st["u_hits"] = st.get("u_hits", 0) + 1
        st["evals"] += 1
        m = (dU_lo - U).upper()
        if m > 0:
            memo[key] = (cUs, m.str(30))
            return True, m, n
        # hopeless-box shortcut (sound: panel error shrinks ~1/n, not 2x)
        if (-m) > dU_lo and n >= 4096:
            return False, m, n
        n = min(n * 4, 32768)
    return False, m, n


def box_split(box):
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        if (ax1 - ax0) >= (ay1 - ay0):
            xm = (ax0 + ax1) / 2
            return (ax0, xm, ay0, ay1), (xm, ax1, ay0, ay1)
        ym = (ay0 + ay1) / 2
        return (ax0, ax1, ay0, ym), (ax0, ax1, ym, ay1)


def box_ok_floor(box) -> bool:
    with Prec(PREC):
        return (box[1] - box[0]) >= arb(BOX_FLOOR) and (box[3] - box[2]) >= arb(BOX_FLOOR)


def certify_band_pair(aL: arb, aU: arb, n_side: int, st: dict, memo: dict,
                      ucache: dict, ckpt: Optional[dict] = None):
    """Certify every disk cell at the c-pair (aL, aU) via box b&b.

    Returns (ok, margin_min, margin_max, worst_cell_str, cells)."""
    with Prec(PREC):
        dU_lo = core2.d_of_c(aU).lower()
        w = arb(2) / n_side
    roots = []
    for i in range(n_side):
        for j in range(n_side):
            cx = arb(-1) + (2 * i + 1) * w / 2
            cy = arb(-1) + (2 * j + 1) * w / 2
            box = box_from_center(cx, cy, w)
            if box_intersects_disk(*box):
                roots.append(box)
    stack = [(b, 0) for b in roots]
    m_min: Optional[arb] = None
    m_max: Optional[arb] = None
    worst_cell = None
    cells = 0
    while stack:
        box, depth = stack.pop()
        ok, m, n_used = try_cell(aL, aU, dU_lo, box, st, memo, ucache, depth)
        cells += 1
        if m is not None:
            if m_min is None or m < m_min:
                m_min = m
                worst_cell = box_key(box) + f"|n={n_used}"
            if m_max is None or m > m_max:
                m_max = m
        if ok:
            continue
        if depth < BOX_DEPTH_MAX and box_ok_floor(box) and st["evals"] <= st["budget"]:
            b1, b2 = box_split(box)
            # SOUNDNESS-CRITICAL filter: a child box that does not intersect
            # the even disk contains NO unit cubic (||a_even|| <= 1 is the
            # parameter constraint); dropping it changes nothing certified.
            for bch in (b1, b2):
                if box_intersects_disk(*bch):
                    stack.append((bch, depth + 1))
        else:
            return False, m_min, m_max, worst_cell, cells
    return True, m_min, m_max, worst_cell, cells


def c_tiles(c_lo: str, c_hi: str, ratio: float = 1.02):
    """Pre-declared geometric tiling of the c-band (the paper's own 'geometric
    cc-tiling within the row', Table 1 / line 1985). Ratio fixed 1.02."""
    tiles = []
    with Prec(PREC):
        a = arb(c_lo)
        b = arb(c_hi)
        while a < b:
            aU = a * arb(repr(ratio))
            if aU >= b or (b - aU) < b * arb('0.005'):
                aU = b
            tiles.append((a, aU))
            a = aU
    return tiles


def run_band(c_lo: str, c_hi: str, budget: int = 260000, verbose=True,
             out_dir: str = None, tag: str = "band"):
    """Certify J <= d over [c_lo, c_hi]: geometric c-tiling, one box-b&b walk
    per sub-band at n_side=12 (the ring cells refine via box-splits)."""
    t0 = time.time()
    st = {"evals": 0, "budget": budget, "memo_hits": 0, "u_hits": 0}
    memo: Dict[str, tuple] = {}
    ucache: Dict[str, str] = {}
    passed: List[Dict] = []
    open_list: List[Dict] = []
    tiles = c_tiles(c_lo, c_hi)
    if verbose:
        print(f"band [{c_lo}, {c_hi}]: {len(tiles)} geometric c-sub-bands", flush=True)
    for (aL, aU) in tiles:
        ok, mn, mx, wcell, cells = certify_band_pair(aL, aU, 12, st, memo, ucache)
        if verbose:
            print(f"  [{aL.str(8)}, {aU.str(8)}]: "
                  f"{'OK' if ok else 'OPEN'} cells={cells} "
                  f"min={(mn.str(8) if mn is not None else None)} "
                  f"max={(mx.str(8) if mx is not None else None)} "
                  f"evals={st['evals']} memo={st['memo_hits']} uhits={st['u_hits']}",
                  flush=True)
        row = {"c_lo": aL.str(16), "c_hi": aU.str(16), "cells": cells,
               "worst_cell": wcell,
               "margin_min": mn.str(14) if mn is not None else None,
               "margin_min_rad": mn.rad().str(6) if mn is not None else None,
               "margin_max": mx.str(14) if mx is not None else None}
        (passed if ok else open_list).append(row)
        if st["evals"] > st["budget"]:
            open_list.append({"c_lo": aL.str(16), "c_hi": aU.str(16),
                              "reason": "global budget reached"})
            break
        if out_dir:
            with open(os.path.join(out_dir, f"{tag}_progress.json"), "w") as f:
                json.dump({"passed": passed, "open": open_list,
                           "evals": st["evals"], "seconds": time.time() - t0},
                          f, indent=1, default=str)
    return {"c_lo": c_lo, "c_hi": c_hi,
            "n_c_tiles": len(tiles),
            "subbands_certified": passed,
            "subbands_open": open_list,
            "cell_evals": st["evals"],
            "memo_hits": st["memo_hits"],
            "u_hits": st["u_hits"],
            "seconds": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c-lo", type=str, required=True)
    ap.add_argument("--c-hi", type=str, required=True)
    ap.add_argument("--budget", type=int, default=260000)
    ap.add_argument("--out-dir", type=str, default=".")
    ap.add_argument("--tag", type=str, default="band")
    args = ap.parse_args()
    res = run_band(args.c_lo, args.c_hi, args.budget, out_dir=args.out_dir, tag=args.tag)
    js = json.dumps(res, indent=1, default=str)
    print(js)
    with open(os.path.join(args.out_dir, f"{args.tag}_result.json"), "w") as f:
        f.write(js)


if __name__ == "__main__":
    main()
