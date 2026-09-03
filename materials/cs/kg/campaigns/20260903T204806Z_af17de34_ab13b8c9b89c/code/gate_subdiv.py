#!/usr/bin/env python3
"""Campaign B — Phase-1 instrument: the pre-registered mq = t^2 panel subdivision.

Envelope V(B, c_L; n): exactly Campaign A's `cell_envelope` (byte-shared
`gate_corr_env`, unmodified, imported as `gce`) EXCEPT the per-panel fold
bound. On every panel (or subpanel after subdivision):

    odd2   = (1 - l^2)_+ outward, l = dist(0, B) 2-D (A's `dist0_box`)
    M-bound:  M_up = max(E_P_upper, W_b),  W_b = sqrt(odd2) * sqrt(K_o(b))
              E_P_upper from A's corner-exact `e_panel_range` (no cap, no 3-D)
    q-bound:  q_up = g2_outward + t_L + c_qlin (see case_exact; >= |||e|-|o||| -
              t on the subpanel; g2_outward = A's outward |E - W| - t_L bound)
    t^2-crossing kill (per subpanel, outward all 四):
        kill iff  max(E_P_upper^2, W_b^2) <= (c_L * a_sub * E_min)_lower
        E_min = outward min over subpanel of (E_lo if E_lo>0 else -E_hi else 0)
        E_min = 0 or the inequality failing => A's exact bound G_P.
    Contribution: GaussWeight(sub-panel) * (M_up - t_L)_+ + kill*0 +
        (1 - kill) * q_full where q_full = V's certified q-hinge ball.

Crossing refinement: the set {s in [0, S_MAX] : odd2*K_o(s) = (c_L*s)^2}
(above the s^2 factor: roots of odd2*(5/2 - s^2 + s^4/6) = c_L^2). On
[sqrt(3), inf) the left side is strictly increasing, so there is exactly one
root > sqrt(3) whenever it exists in-range; it is bracketed by outward
bisection from a doubling bracket (abort if the bracket cannot double past
the root) and splits every panel containing it (no gap, no overlap: left
piece ends at bracket lower, right begins at bracket upper, `gauss_cdf`
evaluated outward per piece). h(s) DECREASES on [0, sqrt(3)], so there is a
second potential root there, bracketed the same way (abort if the h-monotone
domain brackering fails). Every root found is a certified enclosure.

Resume: `--resume CELL_JSON` re-runs the exact leaf recorded in a
Phase-1 tile-3 run (open_cell/open_depth/open_panels/margin fields) instead
of the frozen A frontier; the run log is identical otherwise. Default cell =
A's frozen frontier cell (§2 of the pre-statement), no resume.

Everything soundness-relevant is asserted: no float ever reaches a trusted
path (endpoints exact fmpq via rational literals only), every outward ball is
built with midpoint/radius semantics (never `arb(lo, hi)`), and the V vs U
nonnegativity check (pre-statement §1 item 3) runs on planted cells before
the tile compute.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import importlib.util as _ilu

_HERE = os.path.dirname(os.path.abspath(__file__))
_A_CODE = os.path.join(_HERE, '..', '..', '20260831T231500Z_KgEvenBoxRerun_campA_corr_env', 'code')
_CORE_DIR = os.path.join(_HERE, '..', '..', '20260831T082425Z_kg_direct_d4_restart', 'code')

_spec = _ilu.spec_from_file_location('d4core', os.path.join(_CORE_DIR, 'core.py'))
_dc = _ilu.module_from_spec(_spec)
sys.modules['d4core'] = _dc
_spec.loader.exec_module(_dc)
import d4core
from d4core import Prec

_spec2 = _ilu.spec_from_file_location('gce', os.path.join(_A_CODE, 'gate_corr_env.py'))
gce = _ilu.module_from_spec(_spec2)
_spec2.loader.exec_module(gce)

from flint import arb, fmpq  # noqa: E402

PREC = 256
S_MAX = 8
DEFICIT_ANCHOR = '-1.76291787028287693244667476546e-5'


# ---------------------------------------------------------------------------
# frontier cell (rule 15 exact; A's frozen open_cell display, exact fmpq encode)
# ---------------------------------------------------------------------------

CELL_BOX_FMPQ = (
    (fmpq(-638, 768), fmpq(-637, 768)),   # a0
    (fmpq(429, 768), fmpq(430, 768)),     # a2
)


def frontier_box():
    (lo0, hi0), (lo2, hi2) = CELL_BOX_FMPQ
    with Prec(PREC):
        return (arb(lo0), arb(hi0), arb(lo2), arb(hi2))


def resume_box_from_json(path: str):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    if 'box' in data:
        b = data['box']
        (lo0, hi0), (lo2, hi2) = b[0], b[1]
        with Prec(PREC):
            return (arb(fmpq(lo0[0], lo0[1])), arb(fmpq(hi0[0], hi0[1])),
                    arb(fmpq(lo2[0], lo2[1])), arb(fmpq(hi2[0], hi2[1])))
    raise KeyError('resume json lacks "box" field with fmpq pairs')


# ---------------------------------------------------------------------------
# outward helpers (all midpoint/radius or explicit endpooint arithmetic; NO arb(lo,hi))
# ---------------------------------------------------------------------------

def outward_max2(x: arb, y: arb, bits: int) -> arb:
    """Outward max of two balls: use upper ends for comparison, return the
    ball with larger upper end (an upper bound of max(x, y) is needed)."""
    with Prec(bits):
        return x if bool(x.upper() >= y.upper()) else y


def outward_abs_diff_up(x: arb, y: arb, bits: int) -> arb:
    """Outward upper ball of |x - y| (A's explicit case logic, verbatim)."""
    with Prec(bits):
        diff = x - y
        du = diff.upper() if bool(diff.upper() >= 0) else None
        dl = (-diff).upper() if bool(diff.upper() <= 0 or diff.lower() <= 0) else None
        # replicate A's semantics exactly: they used diff_up = x - y as ball,
        # du = upper if upper >= 0, dl = -lower if lower <= 0
        dxu = diff.upper() if bool(diff.upper() >= 0) else None
        dxl = (-diff).upper() if bool(diff.lower() <= 0) else None
        if dxu is None:
            return dxl
        if dxl is None:
            return dxu
        return dxu if bool(dxu >= dxl) else dxl


def emin_outward(rng: arb, bits: int) -> arb:
    """Outward lower bound of min|e| over the subpanel given the even-range ball."""
    with Prec(bits):
        lo = rng.lower()
        hi = rng.upper()
        if bool(lo > 0):
            return lo.lower() if hasattr(lo, 'lower') else lo
        if bool(hi < 0):
            return (-hi).lower() if hasattr(-hi, 'lower') else (-hi)
        return arb(0)


# ---------------------------------------------------------------------------
# crossing solve (certified brackets)
# ---------------------------------------------------------------------------

def h_mono(s: arb, odd2: arb, bits: int) -> arb:
    with Prec(bits):
        return odd2 * (arb(5) / 2 - s * s + (s ** 4) / 6)


def root_above_sqrt3(c: arb, odd2: arb, bits: int) -> Optional[arb]:
    """Certified outward root of odd2*(5/2 - s^2 + s^4/6) = c^2 on [sqrt(3), inf)."""
    with Prec(bits):
        target = c * c
        lo = arb(3).sqrt()
        hi = lo * 2
        it = 0
        while bool(h_mono(hi, odd2, bits) < target):
            lo_used = hi
            hi = hi * 2
            it += 1
            if it > 200:
                raise AssertionError('crossing bracket failed to double past root (above sqrt3)')
        flo = lo; fhi = hi
        for _ in range(300):
            mid = (flo + fhi) / 2
            if bool(h_mono(mid, odd2, bits) < target):
                flo = mid
            else:
                fhi = mid
            if bool((fhi - flo) < arb(2) ** -212):
                break
        return (flo + fhi) / 2


def root_below_sqrt3(c: arb, odd2: arb, bits: int) -> Optional[arb]:
    """Root of odd2*(5/2 - s^2 + s^4/6) = c^2 in [0, sqrt(3)) if any. h decreases
    strictly there from 5/2*odd2; in the band c^2 >= 16 > 5/2*odd2, so NO root
    exists in-range — the function returns None in that case, by outward test."""
    with Prec(bits):
        target = c * c
        if bool(h_mono(arb(0), odd2, bits) < target):
            return None
        lo = arb(0)
        hi = arb(3).sqrt()
        flo, fhi = lo, hi
        for _ in range(300):
            mid = (flo + fhi) / 2
            if bool(h_mono(mid, odd2, bits) >= target):
                flo = mid
            else:
                fhi = mid
            if bool((fhi - flo) < arb(2) ** -212):
                break
        return (flo + fhi) / 2


def crossings(c: arb, box, bits: int) -> List[arb]:
    """All certified roots of odd2*K_o(s) = (c*s)^2 inside [0, S_MAX] for the box."""
    with Prec(bits):
        ell = gce.dist0_box(box)
        ell = ell.min(arb(1))
        odd2 = (arb(1) - ell * ell).nonnegative_part()
        roots = []
        r0 = root_below_sqrt3(c, odd2, bits)
        if r0 is not None:
            roots.append(r0)
        r1 = root_above_sqrt3(c, odd2, bits)
        if r1 is not None and not bool(r1 > S_MAX):
            roots.append(r1)
        return roots


# ---------------------------------------------------------------------------
# the V envelope (subdivided)
# ---------------------------------------------------------------------------

def _subpanel_pieces(a: arb, b: arb, roots: List[arb], bits: int) -> List[Tuple[arb, arb]]:
    with Prec(bits):
        pieces = []
        cur = a
        for s_star in roots:
            if bool(s_star <= cur) or bool(s_star >= b):
                continue
            slo = s_star.lower()
            shi = s_star.upper()
            if bool(slo <= cur):
                slo = cur
            if bool(shi >= b):
                shi = b
            pieces.append((cur, slo))
            cur = shi
        pieces.append((cur, b))
        return pieces


def cell_envelope_V(c_left: arb, box, n: int, roots: List[arb], bits: int = PREC) -> arb:
    """V(B, c_L; n): A's envelope with the mq = t^2 case/kill per subpanel."""
    with Prec(bits):
        ell = gce.dist0_box(box)
        ell = ell.min(arb(1))
        odd2 = (arb(1) - ell * ell).nonnegative_part()
        odd_r = odd2.sqrt()
        total = arb(0)
        for i in range(n):
            a = arb(fmpq(i * S_MAX, n))
            b = arb(fmpq((i + 1) * S_MAX, n))
            if bool(d4core.K(b, bits).sqrt() <= c_left * a):
                continue
            for (pa, pb) in _subpanel_pieces(a, b, roots, bits):
                panel = pa.union(pb)
                rng, elo, ehi, Eup = gce.e_panel_range(box, panel, bits)
                Wo = odd_r * d4core.K_o(pb, bits).sqrt()
                tL = c_left * pa
                wt = d4core.gauss_cdf(pb, bits) - d4core.gauss_cdf(pa, bits)
                M_up = outward_max2(Eup, Wo, bits)
                # Lemma D.4: f nondecreasing in x AND y — the first hinge keeps the
                # additive odd term (A's form verbatim); the kill zeroes ONLY the
                # second hinge (certified q <= t case), never the first.
                first = (Eup + Wo - tL).nonnegative_part()
                # kill test (outward, midpoint/radius discipline)
                Emin_lo = emin_outward(rng, bits)
                killed = False
                if bool(Emin_lo > 0):
                    thr = (c_left * pa * Emin_lo).lower()
                    if bool(Eup.upper() ** 2 <= thr) and bool(Wo.upper() ** 2 <= thr):
                        killed = True
                if killed:
                    total += wt * first
                else:
                    adu = outward_abs_diff_up(Eup, Wo, bits)
                    second = (adu - tL).nonnegative_part()
                    total += wt * (first + second)
        return total + d4core.tail_account(S_MAX, bits)


# ---------------------------------------------------------------------------
# soundness insurance: V <= U on planted cells (pre-statement §1 item 3)
# ---------------------------------------------------------------------------

PLANTED_CELLS = [
    # (a0 lo/hi, a2 lo/hi) as fmpq pairs — kernel point box, mid-disk box, ring cell
    ((fmpq(2, 3) if False else fmpq(799, 1000), fmpq(801, 1000)),
     (fmpq(-601, 1000), fmpq(-599, 1000))),
    ((fmpq(-638, 768), fmpq(-637, 768)), (fmpq(429, 768), fmpq(430, 768))),
    ((fmpq(-1, 8), fmpq(1, 8)), (fmpq(-1, 8), fmpq(1, 8))),
    ((fmpq(4, 5), fmpq(81, 100)), (fmpq(-61, 100), fmpq(-59, 100))),
]


def check_insurance(c_left: arb) -> None:
    with Prec(PREC):
        roots = crossings(c_left, frontier_box(), PREC)
        for (x0, x1), (y0, y1) in PLANTED_CELLS:
            box = (arb(x0), arb(x1), arb(y0), arb(y1))
            for n in (64, 4096, 65536):
                U = gce.cell_envelope(c_left, box, n, PREC)
                V = cell_envelope_V(c_left, box, n, roots, PREC)
                tol = arb(2) ** -200
                if not bool((V - U).upper() <= arb(0) + tol):
                    raise AssertionError(
                        f'INSURANCE FAILED: V > U on box {x0},{x1},{y0},{y1} n={n}: '
                        f'V-U = {(V - U).str(12)}')


# ---------------------------------------------------------------------------
# tile runner (byte-shared logic with A's certify_tile; V instead of U)
# ---------------------------------------------------------------------------

def evaluate_cell_V(c_left: arb, c_right: arb, box, roots: List[arb], stats: Dict) -> Dict:
    last_upper: Optional[arb] = None
    last_margin: Optional[arb] = None
    last_n = 0
    for n in gce.panel_ladder(c_left):
        if stats['envelope_evals'] >= stats['budget']:
            return {'passed': False, 'reason': 'band budget reached',
                    'margin': last_margin, 'upper': last_upper, 'n_panels': last_n}
        upper = cell_envelope_V(c_left, box, n, roots, PREC)
        margin = gce.strict_margin_lower(c_right, upper)
        stats['envelope_evals'] += 1
        stats['panel_histogram'][str(n)] = stats['panel_histogram'].get(str(n), 0) + 1
        last_upper, last_margin, last_n = upper, margin, n
        if bool(margin > 0):
            return {'passed': True, 'reason': 'strict endpoint separation',
                    'margin': margin, 'upper': upper, 'n_panels': n}
    return {'passed': False, 'reason': 'panel cap reached',
            'margin': last_margin, 'upper': last_upper, 'n_panels': last_n}

def certify_tile_V(c_left: arb, c_right: arb, cell_roots: List[arb], stats: Dict) -> Dict:
    cover = gce.root_boxes(gce.N_SIDE)
    assert len(cover) == 132
    stack = [(box, 0) for box in reversed(cover)]
    tile_stats = {
        'root_boxes': len(cover),
        'boxes_evaluated': 0,
        'certified_leaves': 0,
        'margin_min_ball': None,
        'margin_max_ball': None,
        'worst_certified_cell': None,
        'worst_certified_panels': None,
    }
    while stack:
        box, depth = stack.pop()
        roots = crossings(c_left, box, PREC)
        result = evaluate_cell_V(c_left, c_right, box, roots, stats)
        tile_stats['boxes_evaluated'] += 1
        if result['passed']:
            tile_stats['certified_leaves'] += 1
            margin = result['margin']
            if (tile_stats['margin_min_ball'] is None
                    or bool(margin < tile_stats['margin_min_ball'])):
                tile_stats['margin_min_ball'] = margin
                tile_stats['worst_certified_cell'] = gce.box_key(box)
                tile_stats['worst_certified_panels'] = result['n_panels']
            if (tile_stats['margin_max_ball'] is None
                    or bool(margin > tile_stats['margin_max_ball'])):
                tile_stats['margin_max_ball'] = margin
            continue
        if result['reason'] != 'band budget reached' and gce.can_split(box, depth):
            children = gce.split_box(box)
            for child in reversed(children):
                if gce.box_intersects_disk(child):
                    stack.append((child, depth + 1))
            continue
        margin = result['margin']
        return {
            **tile_stats,
            'passed': False,
            'open_reason': result['reason'],
            'open_cell': gce.box_key(box),
            'open_depth': depth,
            'open_panels': result['n_panels'],
            'open_margin': margin.str(30) if margin is not None else None,
            'open_margin_radius': margin.rad().str(12) if margin is not None else None,
            'remaining_stack_boxes': len(stack),
        }
    return {
        **tile_stats,
        'passed': True,
        'open_reason': None,
        'open_cell': None,
        'open_depth': None,
        'open_panels': None,
        'open_margin': None,
        'open_margin_radius': None,
        'remaining_stack_boxes': 0,
    }


def serialise_tile(c_left: arb, c_right: arb, tile: Dict) -> Dict:
    result = {k: v for k, v in tile.items()
              if k not in ('margin_min_ball', 'margin_max_ball')}
    result['c_lo'] = c_left.str(30)
    result['c_hi'] = c_right.str(30)
    minimum = tile['margin_min_ball']
    maximum = tile['margin_max_ball']
    result['certified_margin_min'] = minimum.str(30) if minimum is not None else None
    result['certified_margin_min_radius'] = (minimum.rad().str(12)
                                             if minimum is not None else None)
    result['certified_margin_max'] = maximum.str(30) if maximum is not None else None
    return result


def run_phase1_frontier(out_dir: str, budget: int, resume: Optional[str],
                        verbose: bool = True) -> Dict:
    """Mechanism validation on the known-failing cell only (Phase 1, step A)."""
    started = time.process_time()
    with Prec(PREC):
        tiles = gce.c_tiles('4.083', '6.0')
        gce.assert_tile_cover('4.083', '6.0', tiles)
        c_left, c_right = tiles[2]
        box = resume_box_from_json(resume) if resume else frontier_box()
        roots = crossings(c_left, box, PREC)
        check_insurance(c_left)
        stats = {'budget': budget, 'envelope_evals': 0, 'panel_histogram': {}}
        if verbose:
            print(f'tile-3 pair ({c_left.str(30)}, {c_right.str(30)}) '
                  f'crossings={[r.str(24) for r in roots]}', flush=True)
        cell_result = evaluate_cell_V(c_left, c_right, box, roots, stats)
        margin = cell_result['margin']
        A_margin = arb(DEFICIT_ANCHOR)
        recovery = None if margin is None else (margin - A_margin)
        out = {
            'phase': '1A-frontier-cell',
            'c_pair': [c_left.str(30), c_right.str(30)],
            'box_key': gce.box_key(box),
            'box_fmpq': [[[str(CELL_BOX_FMPQ[0][0]), str(CELL_BOX_FMPQ[0][1])],
                          [str(CELL_BOX_FMPQ[1][0]), str(CELL_BOX_FMPQ[1][1])]]]
            if resume is None else {'resume': resume},
            'crossings_s': [r.str(30) for r in roots],
            'insurance': 'V <= U on planted cells at n in {64,4096,65536}: PASS',
            'panel_histogram': stats['panel_histogram'],
            'envelope_evals': stats['envelope_evals'],
            'verdict_cell': {
                'passed': cell_result['passed'],
                'reason': cell_result['reason'],
                'n_panels': cell_result['n_panels'],
                'margin': margin.str(30) if margin is not None else None,
                'margin_radius': margin.rad().str(12) if margin is not None else None,
                'recovery_vs_A': recovery.str(30) if recovery is not None else None,
            },
            'a_frozen_open_margin': DEFICIT_ANCHOR,
            'process_seconds': time.process_time() - started,
        }
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'phase1_frontier_cell.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    return out


def run_phase1_tile(out_dir: str, budget: int, verbose: bool = True) -> Dict:
    """Full tile-3 verdict under V (Phase 1, step B)."""
    started = time.process_time()
    with Prec(PREC):
        tiles = gce.c_tiles('4.083', '6.0')
        gce.assert_tile_cover('4.083', '6.0', tiles)
        c_left, c_right = tiles[2]
        check_insurance(c_left)
        stats = {'budget': budget, 'envelope_evals': 0, 'panel_histogram': {}}
        tile = certify_tile_V(c_left, c_right, [], stats)
        row = serialise_tile(c_left, c_right, tile)
        out = {
            'phase': '1B-tile3',
            'row': row,
            'envelope_evals': stats['envelope_evals'],
            'budget': budget,
            'panel_histogram': stats['panel_histogram'],
            'process_seconds': time.process_time() - started,
        }
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'phase1_tile3.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-dir', default='../logs')
    parser.add_argument('--budget', type=int, default=260000)
    parser.add_argument('--mode', choices=['frontier', 'tile'], default='frontier')
    parser.add_argument('--resume', default=None)
    args = parser.parse_args()
    if args.mode == 'frontier':
        result = run_phase1_frontier(args.out_dir, args.budget, args.resume)
    else:
        result = run_phase1_tile(args.out_dir, args.budget)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

