#!/usr/bin/env python3
"""Direct Lemma-D.4 interval envelope for paper inequality (56).

For every panel [a,b], coefficient box B, and c >= cL,

  f(|e|, |o|; c s)
    <= (E + W - cL*a)_+ + (|E-W| - cL*a)_+,

where E and W are certified panel upper bounds.  The Gaussian panel mass is
Phi(b)-Phi(a), exactly enclosed in Arb.  A c-tile [cL,cU] passes only when
U(B,cL).upper() < d(cU).lower() for every even-disk leaf.

Provenance (campaign 20260831T111226Z_kg_direct_d4_highband_panelcache):
surgical cache-only refinement of the frozen as-run source
cs/kg/campaigns/20260831T082425Z_kg_direct_d4_restart/asrun/gate_direct.py.asrun
(sha256 8b364078e83223f5ac9ff8c3001d370ffea588376fce3c0f2abea2e1e26a4cb5;
instrumented by agent KgHighRefine, 2026-08-31).  Same fixed domain
c in [4.083, 6.0], same full even unit disk, same fixed geometric ratio-1.02
c tiling, same S_MAX=8 / Arb-256 / monotone combined-hinge envelope
f(E,W;t).  Every changed line is caching of panel-only Arb quantities,
a finer subdivision floor, provenance, or reporting; the inequality, the
c tiling, the panel cap, the tail, the root cover, and the acceptance rule
are byte-for-byte the frozen instrument's.  This is same-domain instrument
refinement after FAILURE TO CERTIFY, not a new result.

Panel-cache invariants (panel-only Arb quantities, built once per
key (bits, n_panels) and reused for every box and tile at that n):
  cache entry (field under key (bits, n))  | stored content (outward Arb
  enclosure, never a midpoint float)       | consumer
  -----------------------------------------+-----------------------------+
  'panels'          | panel union [a_i, b_{i+1}], outward        | e_panel_upper
  'psi0'/'psi2'     | psi_0/psi_2 on that union, outward Arb     | e_panel_upper
  'kk'              | sqrt(K(b)) at every panel-right b          | cell_envelope skip
  'kko'             | sqrt(K_o(b)) at every panel-right b        | odd_upper
  'left'/'right'    | fmpq-exact endpoints as arb balls          | thresholds/skip
  'mass'            | [Phi(b)-Phi(a)] outward Arb                | total
Cache keying: (precision bits, panel count).  Cache population is
Arb-outward; cached intervals are only ever intersected/shifted/compared,
never widened by float inputs.
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
BOX_FLOOR_Q = fmpq(1, 16384)
BOX_DEPTH_MAX = 40
PANEL_CAP = 32768

# --- static fixed-domain / cache invariants (no numerical evaluation) -----
# Same fixed domain as the frozen d4-restart instrument: band [4.083, 6.0],
# full even unit disk, fixed geometric ratio-1.02 c tiling, S_MAX = 8,
# Arb 256 bits, monotone combined-hinge envelope f(E,W;t).  Only the
# subdivision floor is finer; the cache below never changes any of these.
assert PREC == 256 and S_MAX == 8 and N_SIDE == 12
assert PANEL_CAP == 32768 and BOX_DEPTH_MAX == 40
assert BOX_FLOOR_Q == fmpq(1, 16384)
assert BOX_FLOOR_Q * 16 == fmpq(1, 1024), 'floor must refine, not exceed, 1/1024'
DOMAIN_C_LO = '4.083'
DOMAIN_C_HI = '6.0'
assert arb('4.083') < arb('6.0') and arb('4.083') > arb(4) and arb('6.0') == 6
assert arb(fmpq(4083, 1000)).overlaps(arb('4.083')) and fmpq(4083, 1000) < 4100
TILE_RATIO_STR = '1.02'
assert 4 < arb('4.083') and arb('6.0') <= 8, 'panel ladder start stays 4096 band'
# Panel cache: per-resolution outward Arb enclosures only.  Actual key is
# (bits, n_panels) mapping to a dict of per-kind arb-ball lists; populated
# lazily exactly once per (bits, n).
PANEL_CACHE: Dict[Tuple[int, int], Dict[str, List[arb]]] = {}
ENVELOPE_FUNCTION_ID = 'f(E,W;t) = (E+W-t)_+ + (|E-W|-t)_+; monotone combined hinge, frozen d4-restart'
RELEASE_ENV = 'KG_HIGHBAND_RELEASED'
RELEASE_BUDGET = 260000
RELEASE_TAG = 'band_4p083_6_panelcache'
EXPECTED_OUT_DIR = os.path.realpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs'))


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


def panel_cache(n_panels: int, bits: int = PREC) -> Dict[str, List[arb]]:
    """Build panel-only Arb enclosures exactly once per (bits, n_panels).

    Every entry is an outward Arb enclosure computed with the same
    expressions and fmpq-exact endpoints as the frozen instrument; caching
    changes WHEN the value is computed, never HOW it is computed.  No
    midpoint float is ever used: endpoints come from fmpq, unions and
    square roots via core.psi/core.K/core.K_o/core.gauss_cdf under Prec.
    """
    key = (bits, n_panels)
    cached = PANEL_CACHE.get(key)
    if cached is not None:
        return cached
    assert len(PANEL_CACHE) < 16, 'panel-cache table exceeded pre-set bounds'
    assert 0 < n_panels <= PANEL_CAP
    with Prec(bits):
        left = [arb(fmpq(i * S_MAX, n_panels)) for i in range(n_panels)]
        right = [arb(fmpq((i + 1) * S_MAX, n_panels)) for i in range(n_panels)]
        panels = [left[i].union(right[i]) for i in range(n_panels)]
        psi0 = [core.psi(0, panels[i], bits) for i in range(n_panels)]
        psi2 = [core.psi(2, panels[i], bits) for i in range(n_panels)]
        kk_list = [core.K(right[i], bits).sqrt() for i in range(n_panels)]
        kko_list = [core.K_o(right[i], bits).sqrt() for i in range(n_panels)]
        mass = [core.gauss_cdf(right[i], bits) - core.gauss_cdf(left[i], bits)
                for i in range(n_panels)]
    entry = {
        'left': left,
        'right': right,
        'panels': panels,
        'psi0': psi0,
        'psi2': psi2,
        'kk': kk_list,
        'kko': kko_list,
        'mass': mass,
    }
    PANEL_CACHE[key] = entry
    return entry


def e_panel_upper(box: Tuple[arb, arb, arb, arb],
                  cache_entry: Dict[str, List[arb]], i: int,
                  bits: int = PREC) -> arb:
    """Upper bound for |a0 psi0(s)+a2 psi2(s)| over B x panel (cached).

    Unchanged frozen semantics: the affine even polynomial is evaluated
    outward over the SAME per-panel union enclosure the frozen instrument
    recomputed per call (`panel = a.union(b)`, psi0 = 1, psi2 from core);
    only the per-resolution panel ball and psi2 ball are now built once in
    panel_cache and reused for every box and tile at that resolution.
    """
    ax0, ax1, ay0, ay1 = box
    with Prec(bits):
        a0 = arb(ax0).union(arb(ax1))
        a2 = arb(ay0).union(arb(ay1))
        e_ball = (a0 * cache_entry['psi0'][i]
                  + a2 * cache_entry['psi2'][i])
        return e_ball.abs_upper()


def cdf_weight(a: arb, b: arb, bits: int = PREC) -> arb:
    with Prec(bits):
        return core.gauss_cdf(b, bits) - core.gauss_cdf(a, bits)


def cell_envelope(c_left: arb, box: Tuple[arb, arb, arb, arb],
                  n_panels: int, bits: int = PREC) -> arb:
    """Certified U(B,c_left) for every compatible unit cubic (cached).

    Identical arithmetic to the frozen instrument: same monotone
    combined-hinge terms, same skip rule, same weights and tail.  All
    panel-only quantities (endpoints a,b; panel union; sqrt(K(b));
    sqrt(K_o(b)); [Phi(b)-Phi(a)]) come from panel_cache, which stores the
    outward Arb results of the very same expressions, built once per
    (bits, n_panels) and reused for every box and c-tile.
    """
    with Prec(bits):
        ell = dist0_box(box)
        ell = ell.min(arb(1))
        odd_radius = (arb(1) - ell * ell).nonnegative_part().sqrt()
        cache_entry = panel_cache(n_panels, bits)
        total = arb(0)
        for i in range(n_panels):
            a = cache_entry['left'][i]
            b = cache_entry['right'][i]
            # K=sum_{j<=3} psi_j^2 is increasing for s>=0 because
            # K'(s)=s((s^2-1)^2+2)>0.  This safely skips inactive panels.
            if cache_entry['kk'][i] <= c_left * a:
                continue
            even_upper = e_panel_upper(box, cache_entry, i, bits)
            odd_upper = odd_radius * cache_entry['kko'][i]
            threshold_lower = c_left * a
            first = (even_upper + odd_upper - threshold_lower).nonnegative_part()
            difference_upper = (even_upper - odd_upper).abs_upper()
            second = (difference_upper - threshold_lower).nonnegative_part()
            total += cache_entry['mass'][i] * (first + second)
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
        ratio = arb(TILE_RATIO_STR)
        assert ratio.overlaps(arb('1.02')) and not ratio.overlaps(arb('1.0199'))
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


def release_paths(out_dir: str, tag: str) -> Tuple[str, str]:
    resolved = os.path.realpath(out_dir)
    if resolved != EXPECTED_OUT_DIR:
        raise RuntimeError(
            f'out-dir must be the campaign logs directory {EXPECTED_OUT_DIR!r}')
    if tag != RELEASE_TAG:
        raise RuntimeError(f'tag must be exactly {RELEASE_TAG!r}')
    if not os.path.isdir(resolved):
        raise RuntimeError(f'campaign logs directory does not exist: {resolved!r}')
    return (
        os.path.join(resolved, f'{tag}_progress.json'),
        os.path.join(resolved, f'{tag}_result.json'),
    )


def preflight_release_outputs(out_dir: str, tag: str) -> None:
    progress_path, result_path = release_paths(out_dir, tag)
    for path in (progress_path, result_path):
        for reserved in (path, f'{path}.tmp'):
            if os.path.lexists(reserved):
                raise RuntimeError(
                    f'refusing to overwrite preserved campaign artifact: {reserved}')


def atomic_write_json(path: str, payload: Dict) -> None:
    temporary = f'{path}.tmp'
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(temporary, flags, 0o600)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
            json.dump(payload, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_descriptor = os.open(os.path.dirname(path), os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        # Preserve a partial exclusive temp as evidence and block accidental rerun.
        raise


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
        atomic_write_json(
            os.path.join(out_dir, f'{tag}_progress.json'), progress)

    result = {
        'c_lo': c_lo,
        'c_hi': c_hi,
        'precision_bits': PREC,
        'root_boxes_known_and_asserted': 132,
        'instrument': 'gate_direct_cached panel-cache rerun; same-domain refinement of 20260831T082425Z_kg_direct_d4_restart',
        'envelope_function': ENVELOPE_FUNCTION_ID,
        'panel_cache': {
            'entries': len(PANEL_CACHE),
            'precisions': sorted({bits for (bits, _n) in PANEL_CACHE}),
            'resolutions': sorted({n for (_bits, n) in PANEL_CACHE}),
        },
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
    atomic_write_json(os.path.join(out_dir, f'{tag}_result.json'), result)
    return result


def assert_frozen_resource_limits() -> None:
    """Hard release and single-thread gate before any envelope computation."""
    if not __debug__:
        raise RuntimeError('optimized Python is forbidden for the certified release')
    if os.environ.get(RELEASE_ENV) != '1':
        raise RuntimeError(f'{RELEASE_ENV}=1 is required for the certified release')
    if os.getpriority(os.PRIO_PROCESS, 0) < 10:
        raise RuntimeError('process niceness must be at least 10')
    required = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
    for name in required:
        if os.environ.get(name) != '1':
            raise RuntimeError(
                f'{name} must be exactly 1 (frozen single-thread release command)')


def main() -> None:
    parser = argparse.ArgumentParser(
        description=('Same-domain cached rerun instrument for the fixed '
                     'KG band [4.083,6.0]; computes no result unless '
                     'explicitly invoked with a band command.'))
    parser.add_argument('--c-lo', required=True)
    parser.add_argument('--c-hi', required=True)
    parser.add_argument('--budget', type=int, default=260000)
    parser.add_argument('--out-dir', default='.')
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    assert_frozen_resource_limits()
    if args.c_lo != DOMAIN_C_LO or args.c_hi != DOMAIN_C_HI:
        raise RuntimeError(
            'this cached instrument is domain-frozen to [4.083, 6.0]')
    if args.budget != RELEASE_BUDGET:
        raise RuntimeError(f'budget must be exactly {RELEASE_BUDGET}')
    preflight_release_outputs(args.out_dir, args.tag)
    result = run_band(
        args.c_lo, args.c_hi, args.budget, os.path.realpath(args.out_dir), args.tag)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()


# ---------------------------------------------------------------------------
# Frozen release protocol (COMPUTE PENDING).
#
#   cd cs/kg/campaigns/20260831T111226Z_kg_direct_d4_highband_panelcache/code
#   set -o noclobber
#   nice -n10 env KG_HIGHBAND_RELEASED=1 PYTHONDONTWRITEBYTECODE=1 \
#     OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 \
#     ../../../../.venv/bin/python gate_direct_cached.py \
#     --c-lo 4.083 --c-hi 6.0 --budget 260000 --out-dir ../logs \
#     --tag band_4p083_6_panelcache \
#     > ../logs/band_4p083_6_panelcache_stdout.log 2>&1
#
# The release gate, niceness >= 10, normal Python assertions, fixed band,
# budget, tag, output directory, and five single-thread variables are checked
# before run_band.  Pre-existing progress/result/temp files are refused; JSON
# checkpoints use an exclusive temp, file fsync, atomic rename, and directory
# fsync.  Shell noclobber preserves the separately captured stdout/stderr log.
# The python-flint arb context is single-threaded and the instrument uses
# no BLAS.  Startup controls (startup_controls.py) are NOT re-run by this
# instrument; their PASS is inherited from the parent campaign
# 20260831T082425Z_kg_direct_d4_restart, whose domain/tiling/panel/tail/
# cover/acceptance semantics this file reproduces bit-for-bit.
# ---------------------------------------------------------------------------
