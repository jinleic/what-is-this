#!/usr/bin/env python3
"""Direct Lemma-D.4 interval envelope for paper inequality (56).

For every panel [a,b], coefficient box B, and c >= cL,

  f(|e|, |o|; c s)
    <= (E + W - cL*a)_+ + (|E-W| - cL*a)_+,

where E and W are certified panel upper bounds.  The Gaussian panel mass is
Phi(b)-Phi(a), exactly enclosed in Arb.  A c-tile [cL,cU] passes only when
U(B,cL).upper() < d(cU).lower() for every even-disk leaf.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

from flint import arb, fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import core
from core import Prec

PREC = 256
S_MAX = 8
N_SIDE = 12
BOX_FLOOR_Q = fmpq(1, 1024)
BOX_DEPTH_MAX = 40
PANEL_CAP = 32768


def axis_dist0(lo: arb, hi: arb) -> arb:
    with Prec(PREC):
        if lo <= 0 <= hi:
            return arb(0)
        if lo > 0:
            return lo
        return -hi


def dist0_box(box: Tuple[arb, arb, arb, arb]) -> arb:
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        dx = axis_dist0(ax0, ax1)
        dy = axis_dist0(ay0, ay1)
        return (dx * dx + dy * dy).sqrt()


def box_intersects_disk(box: Tuple[arb, arb, arb, arb]) -> bool:
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        dx = axis_dist0(ax0, ax1)
        dy = axis_dist0(ay0, ay1)
        return dx * dx + dy * dy <= 1


def root_boxes(n_side: int = N_SIDE) -> List[Tuple[arb, arb, arb, arb]]:
    boxes: List[Tuple[arb, arb, arb, arb]] = []
    with Prec(PREC):
        for i in range(n_side):
            ax0 = arb(-1) + arb(fmpq(2 * i, n_side))
            ax1 = arb(-1) + arb(fmpq(2 * (i + 1), n_side))
            for j in range(n_side):
                ay0 = arb(-1) + arb(fmpq(2 * j, n_side))
                ay1 = arb(-1) + arb(fmpq(2 * (j + 1), n_side))
                box = (ax0, ax1, ay0, ay1)
                if box_intersects_disk(box):
                    boxes.append(box)
    return boxes


def box_key(box: Tuple[arb, arb, arb, arb]) -> str:
    return "|".join(x.str(24) for x in box)


def split_box(box: Tuple[arb, arb, arb, arb]):
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        if ax1 - ax0 >= ay1 - ay0:
            mid = (ax0 + ax1) / 2
            return (ax0, mid, ay0, ay1), (mid, ax1, ay0, ay1)
        mid = (ay0 + ay1) / 2
        return (ax0, ax1, ay0, mid), (ax0, ax1, mid, ay1)


def can_split(box: Tuple[arb, arb, arb, arb], depth: int) -> bool:
    if depth >= BOX_DEPTH_MAX:
        return False
    ax0, ax1, ay0, ay1 = box
    with Prec(PREC):
        wx = ax1 - ax0
        wy = ay1 - ay0
        floor = arb(BOX_FLOOR_Q)
        if wx >= wy:
            return wx / 2 >= floor and wy >= floor
        return wy / 2 >= floor and wx >= floor


def e_panel_upper(box: Tuple[arb, arb, arb, arb], panel: arb,
                  bits: int = PREC) -> arb:
    """Upper bound for |a0 psi0(s)+a2 psi2(s)| over B x panel."""
    ax0, ax1, ay0, ay1 = box
    with Prec(bits):
        a0 = arb(ax0).union(arb(ax1))
        a2 = arb(ay0).union(arb(ay1))
        e_ball = a0 * core.psi(0, panel, bits) + a2 * core.psi(2, panel, bits)
        return e_ball.abs_upper()


def cdf_weight(a: arb, b: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return core.gauss_cdf(b, bits) - core.gauss_cdf(a, bits)


def cell_envelope(c_left: arb, box: Tuple[arb, arb, arb, arb],
                  n_panels: int, bits: int = PREC) -> arb:
    """Certified U(B,c_left) for every compatible unit cubic."""
    with Prec(bits):
        ell = dist0_box(box)
        ell = ell.min(arb(1))
        odd_radius = (arb(1) - ell * ell).nonnegative_part().sqrt()
        total = arb(0)
        for i in range(n_panels):
            a = arb(fmpq(i * S_MAX, n_panels))
            b = arb(fmpq((i + 1) * S_MAX, n_panels))
            # K=sum_{j<=3} psi_j^2 is increasing for s>=0 because
            # K'(s)=s((s^2-1)^2+2)>0.  This safely skips inactive panels.
            if core.K(b, bits).sqrt() <= c_left * a:
                continue
            panel = a.union(b)
            even_upper = e_panel_upper(box, panel, bits)
            odd_upper = odd_radius * core.K_o(b, bits).sqrt()
            threshold_lower = c_left * a
            first = (even_upper + odd_upper - threshold_lower).nonnegative_part()
            difference_upper = (even_upper - odd_upper).abs_upper()
            second = (difference_upper - threshold_lower).nonnegative_part()
            total += cdf_weight(a, b, bits) * (first + second)
        return total + core.tail_account(S_MAX, bits)


def strict_margin_lower(c_right: arb, upper: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return core.d_of_c(c_right, bits).lower() - upper.upper()


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


def c_tiles(c_lo: str, c_hi: str) -> List[Tuple[arb, arb]]:
    tiles: List[Tuple[arb, arb]] = []
    with Prec(PREC):
        left = arb(c_lo)
        top = arb(c_hi)
        ratio = arb('1.02')
        while left < top:
            right = left * ratio
            if right >= top:
                right = top
            tiles.append((left, right))
            left = right
    return tiles


def assert_tile_cover(c_lo: str, c_hi: str,
                      tiles: List[Tuple[arb, arb]]) -> None:
    assert tiles
    with Prec(PREC):
        assert (tiles[0][0] - arb(c_lo)).contains(0)
        assert (tiles[-1][1] - arb(c_hi)).contains(0)
        for previous, current in zip(tiles, tiles[1:]):
            assert (previous[1] - current[0]).contains(0)
            assert previous[0] < previous[1]
        assert tiles[-1][0] < tiles[-1][1]


def evaluate_cell(c_left: arb, c_right: arb,
                  box: Tuple[arb, arb, arb, arb], stats: Dict) -> Dict:
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
