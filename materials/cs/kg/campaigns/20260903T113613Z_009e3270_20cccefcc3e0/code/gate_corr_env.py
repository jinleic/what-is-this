#!/usr/bin/env python3
"""Campaign A: corrected-envelope recertification of band [4.083, 6.0].

Derived from the paper text (pre-statement quotes, lines 1972-1981), NOT from
predecessor code: corner-exact 2-D even-box envelope E_P (affine vertex range
on B x P), paper odd radius r_o = sqrt(1-l^2) with l = dist(0,B) 2-D, Lemma D.4
direct combined-hinge integrand G_P, exact CDF panel weights, closed tail.

Skill discipline: no |e|-cap anywhere; no 3-D circumradius anywhere.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

from flint import arb, fmpq

PREC = 256
ANALYTIC_BITS = 320
S_MAX = 8
N_SIDE = 12
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    'd4core', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', '..', '20260831T082425Z_kg_direct_d4_restart', 'code', 'core.py'))
_d4core_mod = _ilu.module_from_spec(_spec)
sys.modules['d4core'] = _d4core_mod
_spec.loader.exec_module(_d4core_mod)
import d4core
from d4core import Prec

BOX_FLOOR_Q = fmpq(1, 1024)
PANEL_CAP = 131072
BOX_DEPTH_MAX = 40
KERNEL_EPS = arb(2) ** -260 if False else None  # set inside Prec context below


# --------------------------------------------------------------------------
# outward Arb helpers (all inside Prec)
# --------------------------------------------------------------------------

def axis_dist0(lo: arb, hi: arb) -> arb:
    """Distance from 0 to [lo, hi] as an outward arb (exact if endpoints exact)."""
    with Prec(PREC):
        if lo <= 0 <= hi:
            return arb(0)
        if lo > 0:
            return lo
        return -hi


def dist0_box(box) -> arb:
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        dx = axis_dist0(ax0, ax1)
        dy = axis_dist0(ay0, ay1)
        return (dx * dx + dy * dy).sqrt()


def box_intersects_disk(box) -> bool:
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        dx = axis_dist0(ax0, ax1)
        dy = axis_dist0(ay0, ay1)
        # nearest point of box to origin inside closed unit disk  <=>  intersects
        if dx * dx + dy * dy <= 1:
            return True
        # else check farthest is irrelevant: box misses disk iff nearest point
        # outside OR farthest in the all-negative "both axes straddle 0" cases
        # the nearest-point test is exactly correct for L2 distance to a box
        return False


def root_boxes(n_side: int = N_SIDE) -> List:
    boxes = []
    step = fmpq(2, n_side)
    for i in range(n_side):
        for j in range(n_side):
            ax0 = fmpq(-1) + step * i
            ax1 = fmpq(-1) + step * (i + 1)
            ay0 = fmpq(-1) + step * j
            ay1 = fmpq(-1) + step * (j + 1)
            with Prec(PREC):
                box = (arb(ax0), arb(ax1), arb(ay0), arb(ay1))
            if box_intersects_disk(box):
                boxes.append(box)
    return boxes


def box_key(box) -> str:
    return "|".join(x.str(24) for x in box)


def split_box(box):
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        if (ax1 - ax0) >= (ay1 - ay0):
            mid = (ax0 + ax1) / 2
            return (ax0, mid, ay0, ay1), (mid, ax1, ay0, ay1)
        mid = (ay0 + ay1) / 2
        return (ax0, ax1, ay0, mid), (ax0, ax1, mid, ay1)


def can_split(box, depth: int) -> bool:
    if depth >= BOX_DEPTH_MAX:
        return False
    with Prec(PREC):
        ax0, ax1, ay0, ay1 = box
        wx = ax1 - ax0
        wy = ay1 - ay0
        floor = arb(BOX_FLOOR_Q)
        if wx >= wy:
            return wx / 2 >= floor
        return wy / 2 >= floor


def e_panel_range(box, panel: arb, bits: int = PREC):
    """Corner-exact outward enclosure [E_lo, E_hi] of e(s) over B x panel.

    e(s) = a0*psi0(s) + a2*psi2(s), psi0 = 1, psi2(s) = (s^2-1)/sqrt(2).
    Affine in (a0, a2, psi2(s)); psi2 is monotone increasing on [0, S_MAX]
    (psi2'(s) = sqrt(2)*s >= 0), so psi2(s) ranges over
    [psi2(a), psi2(b)] on panel P = [a, b]; extrema of the affine function on
    the 3-D vertex set {a0 endpoints} x {a2 endpoints} x {psi2 endpoints}.
    Returns (lo, hi, upper_of_abs) as outward Arb balls.
    """
    ax0, ax1, ay0, ay1 = box
    with Prec(bits):
        a = panel.lower() if hasattr(panel, 'lower') else panel
        b = panel.upper() if hasattr(panel, 'lower') else panel
        p2a = d4core.psi(2, a, bits)
        p2b = d4core.psi(2, b, bits)
        verts = []
        for a0 in (ax0, ax1):
            for a2 in (ay0, ay1):
                for p2 in (p2a, p2b):
                    verts.append(a0 + a2 * p2, )
        lo = verts[0]
        hi = verts[0]
        for v in verts[1:]:
            lo = lo.union(v) if False else lo  # placeholder, replaced below
        # proper outward enclosure: union of vertex balls
        rng = None
        for v in verts:
            rng = v if rng is None else rng.union(v)
        # vertices are exact-ish balls; the true range on the box is
        # [min vertex value, max vertex value]. Outward enclosure: use the
        # union's lower and upper as the interval ends (contains all vertices;
        # by affinity contains the whole range up to rounding).
        elo = arb(min(v.lower() for v in verts))
        ehi = arb(max(v.upper() for v in verts))
        # outward |e| upper: max(|elo|, |ehi|) computed outward
        up = ehi if ehi >= 0 and elo >= 0 else (
             -elo if ehi <= 0 and elo <= 0 else (ehi if ehi >= -elo else -elo))
        # safer outward: ehi.abs handled by explicit comparisons above; ensure
        # up is an upper bound of |x| for all x in [elo, ehi]:
        cand1 = ehi.upper() if hasattr(ehi, 'upper') else ehi
        cand2 = (-elo).upper() if hasattr(elo, 'upper') else -elo
        up = cand1 if cand1 >= cand2 else cand2
        return rng, elo, ehi, up


def cdf_weight(a: arb, b: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return d4core.gauss_cdf(b, bits) - d4core.gauss_cdf(a, bits)


def cell_envelope(c_left: arb, box, n_panels: int, bits: int = PREC) -> arb:
    """Certified U(B, c_left) >= J(c_left, p) for every unit cubic with even
    part in B (corrected envelope — no |e| cap, no 3-D circumradius)."""
    with Prec(bits):
        ell = dist0_box(box)
        ell = ell.min(arb(1))
        odd_radius_sq = (arb(1) - ell * ell).nonnegative_part()
        odd_radius = odd_radius_sq.sqrt()
        total = arb(0)
        for i in range(n_panels):
            a = arb(fmpq(i * S_MAX, n_panels))
            b = arb(fmpq((i + 1) * S_MAX, n_panels))
            # panel skip: no active fold unless sqrt(K(b)) > c_left * a
            if d4core.K(b, bits).sqrt() <= c_left * a:
                continue
            panel = a.union(b)
            _, _, _, even_upper = e_panel_range(box, panel, bits)
            odd_upper = odd_radius * d4core.K_o(b, bits).sqrt()
            threshold_lower = c_left * a
            first = (even_upper + odd_upper - threshold_lower).nonnegative_part()
            diff_up = (even_upper - odd_upper)
            # outward upper of |E - W|: two cases
            du = diff_up.upper() if bool(diff_up.upper() >= 0) else None
            dl = (-diff_up.lower()) if bool(diff_up.lower() <= 0) else None
            if du is None:
                abs_diff_up = dl
            elif dl is None:
                abs_diff_up = du
            else:
                abs_diff_up = du if du >= dl else dl
            second = (abs_diff_up - threshold_lower).nonnegative_part()
            total += cdf_weight(a, b, bits) * (first + second)
        return total + d4core.tail_account(S_MAX, bits)


def strict_margin_lower(c_right: arb, upper: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return d4core.d_of_c(c_right, bits).lower() - upper.upper()


def panel_ladder(c_left: arb) -> List[int]:
    with Prec(PREC):
        if c_left <= 1:
            start = 2048
        elif c_left <= 4:
            start = 1024
        elif c_left <= 8:
            start = 4096
        else:
            start = 8192
    values: List[int] = []
    n = start
    while True:
        n = min(n, PANEL_CAP)
        if not values or values[-1] != n:
            values.append(n)
        if n == PANEL_CAP:
            return values
        n *= 4


TILE_RATIO = '1.008'


def c_tiles(c_lo: str, c_hi: str) -> List[Tuple[arb, arb]]:
    tiles: List[Tuple[arb, arb]] = []
    with Prec(PREC):
        left = arb(c_lo)
        top = arb(c_hi)
        ratio = arb(TILE_RATIO)
        while left < top:
            right = left * ratio
            if right >= top:
                right = top
            tiles.append((left, right))
            left = right
    return tiles


def assert_tile_cover(c_lo: str, c_hi: str, tiles) -> None:
    assert tiles
    with Prec(PREC):
        assert (tiles[0][0] - arb(c_lo)).contains(0)
        assert (tiles[-1][1] - arb(c_hi)).contains(0)
        for previous, current in zip(tiles, tiles[1:]):
            assert (previous[1] - current[0]).contains(0)
            assert previous[0] < previous[1]
        assert tiles[-1][0] < tiles[-1][1]


def evaluate_cell(c_left: arb, c_right: arb, box, stats: Dict) -> Dict:
    last_upper: Optional[arb] = None
    last_margin: Optional[arb] = None
    last_n = 0
    for n in panel_ladder(c_left):
        if stats['envelope_evals'] >= stats['budget']:
            return {
                'passed': False,
                'reason': 'band budget reached',
                'margin': last_margin,
                'upper': last_upper,
                'n_panels': last_n,
            }
        upper = cell_envelope(c_left, box, n)
        margin = strict_margin_lower(c_right, upper)
        stats['envelope_evals'] += 1
        stats['panel_histogram'][str(n)] = stats['panel_histogram'].get(str(n), 0) + 1
        last_upper, last_margin, last_n = upper, margin, n
        if margin > 0:
            return {
                'passed': True,
                'reason': 'strict endpoint separation',
                'margin': margin,
                'upper': upper,
                'n_panels': n,
            }
    return {
        'passed': False,
        'reason': 'panel cap reached',
        'margin': last_margin,
        'upper': last_upper,
        'n_panels': last_n,
    }


def _update_leaf_extrema(tile_stats: Dict, box, result: Dict) -> None:
    margin = result['margin']
    if margin is None:
        return
    if tile_stats['margin_min_ball'] is None or margin < tile_stats['margin_min_ball']:
        tile_stats['margin_min_ball'] = margin
        tile_stats['worst_certified_cell'] = box_key(box)
        tile_stats['worst_certified_panels'] = result['n_panels']
    if tile_stats['margin_max_ball'] is None or margin > tile_stats['margin_max_ball']:
        tile_stats['margin_max_ball'] = margin


def certify_tile(c_left: arb, c_right: arb, stats: Dict) -> Dict:
    roots = root_boxes(N_SIDE)
    assert len(roots) == 132
    stack = [(box, 0) for box in reversed(roots)]
    tile_stats = {
        'root_boxes': len(roots),
        'boxes_evaluated': 0,
        'certified_leaves': 0,
        'margin_min_ball': None,
        'margin_max_ball': None,
        'worst_certified_cell': None,
        'worst_certified_panels': None,
    }
    while stack:
        box, depth = stack.pop()
        result = evaluate_cell(c_left, c_right, box, stats)
        tile_stats['boxes_evaluated'] += 1
        if result['passed']:
            tile_stats['certified_leaves'] += 1
            _update_leaf_extrema(tile_stats, box, result)
            continue
        if result['reason'] != 'band budget reached' and can_split(box, depth):
            children = split_box(box)
            for child in reversed(children):
                if box_intersects_disk(child):
                    stack.append((child, depth + 1))
            continue
        margin = result['margin']
        return {
            **tile_stats,
            'passed': False,
            'open_reason': result['reason'],
            'open_cell': box_key(box),
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
    result['certified_margin_min_radius'] = minimum.rad().str(12) if minimum is not None else None
    result['certified_margin_max'] = maximum.str(30) if maximum is not None else None
    return result


def run_band(c_lo: str, c_hi: str, budget: int, out_dir: str,
             tag: str, verbose: bool = True) -> Dict:
    started = time.time()
    tiles = c_tiles(c_lo, c_hi)
    assert_tile_cover(c_lo, c_hi, tiles)
    stats: Dict = {
        'budget': budget,
        'envelope_evals': 0,
        'panel_histogram': {},
    }
    certified: List[Dict] = []
    open_tiles: List[Dict] = []
    if verbose:
        print(f'band [{c_lo},{c_hi}] tiles={len(tiles)}', flush=True)
    for index, (c_left, c_right) in enumerate(tiles):
        if stats['envelope_evals'] >= budget:
            open_tiles.append({
                'c_lo': c_left.str(30),
                'c_hi': c_hi,
                'open_reason': 'band budget reached before tile',
                'remaining_named_tiles': len(tiles) - index,
            })
            break
        tile = certify_tile(c_left, c_right, stats)
        row = serialise_tile(c_left, c_right, tile)
        (certified if tile['passed'] else open_tiles).append(row)
        if verbose:
            margin = row['certified_margin_min'] if tile['passed'] else row['open_margin']
            print(f"  {index + 1}/{len(tiles)} {'PASS' if tile['passed'] else 'OPEN'} "
                  f"boxes={row['boxes_evaluated']} margin={margin} "
                  f"evals={stats['envelope_evals']}", flush=True)
        progress = {
            'c_lo': c_lo,
            'c_hi': c_hi,
            'tiles_total': len(tiles),
            'tiles_processed': index + 1,
            'tiles_certified': certified,
            'tiles_open': open_tiles,
            'envelope_evals': stats['envelope_evals'],
            'panel_histogram': stats['panel_histogram'],
            'seconds': time.time() - started,
        }
        with open(os.path.join(out_dir, f'{tag}_progress.json'), 'w', encoding='utf-8') as f:
            json.dump(progress, f, indent=2)
            f.write('\n')

    result = {
        'c_lo': c_lo,
        'c_hi': c_hi,
        'precision_bits': PREC,
        'root_boxes_known_and_asserted': 132,
        'tiles_total': len(tiles),
        'tiles_certified_count': len(certified),
        'tiles_open_count': len(open_tiles),
        'band_passed': len(certified) == len(tiles) and not open_tiles,
        'tiles_certified': certified,
        'tiles_open': open_tiles,
        'envelope_evals': stats['envelope_evals'],
        'budget': budget,
        'panel_histogram': stats['panel_histogram'],
        'seconds': time.time() - started,
    }
    with open(os.path.join(out_dir, f'{tag}_result.json'), 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--c-lo', required=True)
    parser.add_argument('--c-hi', required=True)
    parser.add_argument('--budget', type=int, default=260000)
    parser.add_argument('--out-dir', default='.')
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    result = run_band(args.c_lo, args.c_hi, args.budget, args.out_dir, args.tag)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
