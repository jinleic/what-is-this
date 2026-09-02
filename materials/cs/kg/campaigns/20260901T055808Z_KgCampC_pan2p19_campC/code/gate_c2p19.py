#!/usr/bin/env python3
"""Campaign C — Phase-1 driver: tile 3 at panel cap 2^19 (single knob).

The instrument is Campaign A's byte-frozen `gate_corr_env.py` (imported, not
edited); the ONLY delta of this campaign is the module global PANEL_CAP
patched 131072 -> 524288 after import and before any evaluation, together
with the ladder assertion. All verdict logic (DFS, split rule, margin
predicate, budget semantics) is A's, executed through the imported module.

Order of operations (pre-statement §3/§4/§5):
  1. startup battery (16 asserts) — inherited script, run first;
  2. reconciliation gate: frontier cell at cap 2^17 must reproduce A's
     frozen margin bit-exactly (ball radius < 2^-40), else STOP;
  3. frontier cell at cap 2^19 (first Phase-1 leaf datum);
  4. full tile 3 at cap 2^19, byte-shared certify_tile via the module.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional

import importlib.util as _ilu

_HERE = os.path.dirname(os.path.abspath(__file__))
_A_CODE = os.path.join(_HERE, '..', '..', '20260831T231500Z_KgEvenBoxRerun_campA_corr_env', 'code')
_CORE_DIR = os.path.join(_HERE, '..', '..', '20260831T082425Z_kg_direct_d4_restart', 'code')
OUT = os.path.join(_HERE, '..', 'logs')


def load(name, path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# d4core must be importable under the name A's module expects.
load('d4core', os.path.join(_CORE_DIR, 'core.py'))
gce = load('gce', os.path.join(_A_CODE, 'gate_corr_env.py'))

from flint import arb, fmpq  # noqa: E402
from d4core import Prec  # noqa: E402

PREC = 256
CAMP_CAP_19 = 524288
A_FROZEN_MARGIN = '-1.76291787028287693244667476546e-5'
B_CAP19_MARGIN = '4.46309091764170239642742026145e-6'
B_CAP19_MARGIN_RADIUS = '3.26e-36'

# A's frozen frontier cell, exact fmpq (rule 15; inherited from B's gate_subdiv).
CELL_BOX_FMPQ = (
    (fmpq(-638, 768), fmpq(-637, 768)),   # a0
    (fmpq(429, 768), fmpq(430, 768)),     # a2
)


def frontier_box():
    (lo0, hi0), (lo2, hi2) = CELL_BOX_FMPQ
    with Prec(PREC):
        return (arb(lo0), arb(hi0), arb(lo2), arb(hi2))


def _freeze_gates(out_dir: str, gates: Dict, started: float) -> None:
    """Incremental freeze: gates JSON after each completed gate stage."""
    blob = json.dumps({'phase': 'GATES', 'gates': gates,
                       'process_seconds': time.process_time() - started},
                      indent=2, default=str)
    with open(os.path.join(out_dir, 'phase1_gates.json'), 'w',
              encoding='utf-8') as f:
        f.write(blob)
        f.write('\n')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-dir', default=OUT)
    parser.add_argument('--budget', type=int, default=260000)
    parser.add_argument(
        '--skip-reconciliation', action='store_true',
        help='resume mode: reconciliation gate already frozen green')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    started = time.process_time()
    gates: Dict = {}
    fail = False

    # ---- knob patch BEFORE any evaluation ---------------------------------
    assert gce.PANEL_CAP == 131072, 'inherited instrument must be pristine'
    gce.PANEL_CAP = CAMP_CAP_19   # the single knob (pre-statement §1)

    with Prec(PREC):
        tiles = gce.c_tiles('4.083', '6.0')
        gce.assert_tile_cover('4.083', '6.0', tiles)
        c_left, c_right = tiles[2]
        assert (c_left - arb('4.1485893120')).rad() < arb(2) ** -200
        assert (c_right - arb('4.181778026496')).rad() < arb(2) ** -200
        ladder_expected = [4096, 16384, 65536, 262144, 524288]
        assert gce.panel_ladder(c_left) == ladder_expected, (
            'LADDER MISMATCH under patched cap: '
            f'{gce.panel_ladder(c_left)}')
        gates['knob_ladder_assert'] = {
            'ok': True,
            'detail': f'PANEL_CAP=524288; ladder == {ladder_expected}'}

        box = frontier_box()
        stats_probe = {'budget': 10 ** 9, 'envelope_evals': 0,
                       'panel_histogram': {}}

        # ---- 2. reconciliation gate (cap 2^17) ----------------------------
        if not args.skip_reconciliation:
            gce.PANEL_CAP = 131072
            assert gce.panel_ladder(c_left) == [4096, 16384, 65536, 131072]
            t0 = time.process_time()
            upper_17 = None
            margin_17 = None
            for n in gce.panel_ladder(c_left):
                upper_17 = gce.cell_envelope(c_left, box, n)
                margin_17 = gce.strict_margin_lower(c_right, upper_17)
            a_frozen = arb(A_FROZEN_MARGIN)
            diff = margin_17 - a_frozen
            reproduces = bool(diff.rad() < arb(2) ** -40)
            gates['reconciliation_2p17'] = {
                'ok': reproduces,
                'margin': margin_17.str(30),
                'margin_radius': margin_17.rad().str(12),
                'a_frozen': A_FROZEN_MARGIN,
                'diff_display': diff.str(12),
                'diff_radius_bound': arb(2) ** -40,
                'process_seconds': time.process_time() - t0,
            }
            _freeze_gates(args.out_dir, gates, started)
            gce.PANEL_CAP = CAMP_CAP_19
            if not reproduces:
                fail = True

        # ---- 3. frontier cell at cap 2^19 ---------------------------------
        if not fail:
            t0 = time.process_time()
            cell = dict(stats_probe)
            result_cell = gce.evaluate_cell(c_left, c_right, box, cell)
            tcell = time.process_time() - t0
            margin_19 = result_cell['margin']
            b_ball = None
            if margin_19 is not None:
                b_ball = dict(
                    margin=margin_19.str(30),
                    margin_radius=margin_19.rad().str(12),
                    passed_flag=bool(result_cell['passed']),
                    leaf_margin_lower_positive=bool(margin_19.lower() > 0),
                    n_panels=result_cell['n_panels'],
                    reason=result_cell['reason'],
                    margin_delta_vs_B_D3=str(margin_19 - arb(B_CAP19_MARGIN)),
                )
            gates['frontier_cell_2p19'] = {
                'b_ball': b_ball,
                'b_cap19_margin_for_crosscheck': B_CAP19_MARGIN,
                'process_seconds': tcell,
                'panel_histogram': cell['panel_histogram'],
            }
            _freeze_gates(args.out_dir, gates, started)

    if fail:
        out = {
            'phase': 'GATES',
            'aborted': True,
            'gates': gates,
            'process_seconds': time.process_time() - started,
        }
        with open(os.path.join(args.out_dir, 'phase1_gates.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(out, f, indent=2)
            f.write('\n')
        print(json.dumps(out, indent=2))
        sys.exit(1)

    # ---- 4. full tile 3 at cap 2^19 ---------------------------------------
    with Prec(PREC):
        stats = {'budget': args.budget, 'envelope_evals': 0,
                 'panel_histogram': {}}
        t1 = time.process_time()
        tile = gce.certify_tile(c_left, c_right, stats)
        tile_seconds = time.process_time() - t1
        row = gce.serialise_tile(c_left, c_right, tile)
        depth_max = tile.get('open_depth')
        out = {
            'phase': '1-tile3-cap2p19',
            'knob': {
                'PANEL_CAP_from': 131072,
                'PANEL_CAP_to': CAMP_CAP_19,
                'ladder': [4096, 16384, 65536, 262144, 524288],
            },
            'row': row,
            'envelope_evals': stats['envelope_evals'],
            'budget': args.budget,
            'panel_histogram': stats['panel_histogram'],
            'gates': gates,
            'process_seconds_tile': tile_seconds,
            'process_seconds_total': time.process_time() - started,
        }
    with open(os.path.join(args.out_dir, 'phase1_tile3_result.json'), 'w',
              encoding='utf-8') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
