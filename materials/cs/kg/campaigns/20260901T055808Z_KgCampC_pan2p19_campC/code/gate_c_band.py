#!/usr/bin/env python3
"""Campaign C — Phase-2 driver: resume tiles 4-49 at cap 2^19.

Same single knob (gce.PANEL_CAP = 524288). Tile-level incremental freezes:
logs/inproc_tile_NNN.json after each tile; running summary
logs/band_resume_progress.json. Start tile configurable; previously frozen
tiles are skipped (resumable at tile granularity). Stop conditions
(pre-statement §6): first failing tile; tile 49 certified; global CPU budget.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List

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


load('d4core', os.path.join(_CORE_DIR, 'core.py'))
gce = load('gce', os.path.join(_A_CODE, 'gate_corr_env.py'))

from d4core import Prec  # noqa: E402

PREC = 256
CAMP_CAP_19 = 524288
GLOBAL_CPU_BUDGET_S = 40 * 3600  # pre-statement §6: 40 CPU-hours


def load_frozen_tile(out_dir: str, idx: int) -> Dict:
    path = os.path.join(out_dir, f'inproc_tile_{idx:03d}.json')
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-dir', default=OUT)
    parser.add_argument('--budget', type=int, default=260000)
    parser.add_argument('--start-tile', type=int, default=4)
    parser.add_argument('--cpu-budget', type=float, default=GLOBAL_CPU_BUDGET_S)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    assert gce.PANEL_CAP == 131072, 'inherited instrument must be pristine'
    gce.PANEL_CAP = CAMP_CAP_19

    started = time.process_time()
    tiles = None
    with Prec(PREC):
        tiles = gce.c_tiles('4.083', '6.0')
        gce.assert_tile_cover('4.083', '6.0', tiles)
        assert len(tiles) == 49
        assert gce.panel_ladder(tiles[0][0]) == [4096, 16384, 65536, 262144, 524288]

    certified: List[Dict] = []
    open_tiles: List[Dict] = []
    cpu_used = 0.0
    stop_reason = None

    # resume: adopt previously frozen tiles
    for idx in range(1, args.start_tile):
        path = os.path.join(args.out_dir, f'inproc_tile_{idx:03d}.json')
        if os.path.exists(path):
            row = load_frozen_tile(args.out_dir, idx)
            (certified if row.get('passed') else open_tiles).append(row)
            cpu_used += float(row.get('process_seconds', 0.0))

    for index in range(args.start_tile - 1, len(tiles)):
        c_left, c_right = tiles[index]
        if cpu_used >= args.cpu_budget:
            stop_reason = 'global CPU budget exhausted before tile'
            open_tiles.append({
                'index': index + 1,
                'c_lo': c_left.str(30), 'c_hi': c_right.str(30),
                'open_reason': stop_reason,
                'remaining_named_tiles': len(tiles) - index,
            })
            break
        with Prec(PREC):
            stats = {'budget': args.budget, 'envelope_evals': 0,
                     'panel_histogram': {}}
            t0 = time.process_time()
            tile = gce.certify_tile(c_left, c_right, stats)
            tile_seconds = time.process_time() - t0
            row = gce.serialise_tile(c_left, c_right, tile)
            row['tile_index_1based'] = index + 1
            row['process_seconds'] = tile_seconds
            row['envelope_evals'] = stats['envelope_evals']
            row['panel_histogram'] = stats['panel_histogram']
            row['c_label'] = f'{c_left.str(24)}..{c_right.str(24)}'
        cpu_used += tile_seconds
        (certified if tile['passed'] else open_tiles).append(row)
        tile_path = os.path.join(args.out_dir, f'inproc_tile_{index + 1:03d}.json')
        with open(tile_path, 'w', encoding='utf-8') as f:
            json.dump(row, f, indent=2)
            f.write('\n')
        progress = {
            'band': '[4.083, 6.0] resume tiles %d-49 at cap 2^19' % args.start_tile,
            'cap': CAMP_CAP_19,
            'tiles_total': 49,
            'tiles_completed': index + 1,
            'tiles_certified': certified,
            'tiles_open': open_tiles,
            'envelope_evals_total': sum(r.get('envelope_evals', 0)
                                        for r in certified + open_tiles
                                        if isinstance(r, dict)),
            'cpu_seconds_used': cpu_used,
            'cpu_budget_seconds': args.cpu_budget,
            'stop_reason': stop_reason,
        }
        with open(os.path.join(args.out_dir, 'band_resume_progress.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(progress, f, indent=2)
            f.write('\n')
        print(f"  tile {index + 1}/49 {'PASS' if tile['passed'] else 'OPEN'} "
              f"margin={row.get('certified_margin_min') or row.get('open_margin')} "
              f"cpu={tile_seconds:.1f}s cum={cpu_used:.1f}s", flush=True)
        if not tile['passed']:
            stop_reason = 'first failing tile (pre-registered stop rule)'
            break
    else:
        stop_reason = 'tile 49 certified' if not open_tiles else stop_reason

    result = {
        'band': '[4.083, 6.0]',
        'tiles_total': 49,
        'cap': CAMP_CAP_19,
        'tiles_certified_count': len([r for r in certified if r.get('passed')]),
        'band_passed': bool(len(certified) == 49 and not
                            [r for r in certified + open_tiles
                             if isinstance(r, dict) and not r.get('passed')]),
        'stop_reason': stop_reason,
        'tiles_certified': certified,
        'tiles_open': open_tiles,
        'cpu_seconds_used': cpu_used,
        'cpu_budget_seconds': args.cpu_budget,
        'process_seconds_this_run': time.process_time() - started,
    }
    with open(os.path.join(args.out_dir, 'band_resume_result.json'), 'w',
              encoding='utf-8') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
