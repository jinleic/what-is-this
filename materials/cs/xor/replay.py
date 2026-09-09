#!/usr/bin/env python3
"""Replay archived XOR scripts into fresh output; the parent owns failure receipts."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS')


def save_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def supervise(command, receipt, wall_seconds, *, environment=None):
    """A killed child cannot own its terminal receipt; retain incomplete output."""
    receipt = Path(receipt)
    started = time.monotonic()
    reason = None
    code = None
    with receipt.with_suffix('.log').open('w') as log:
        try:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       env=environment)
        except OSError as error:
            reason = f'LAUNCH-ERROR: {error}'
        else:
            try:
                code = process.wait(timeout=wall_seconds)
            except subprocess.TimeoutExpired:
                reason = 'WALL-TIMEOUT'
                process.terminate()
                try:
                    code = process.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    process.kill()
                    code = process.wait()
            if reason is None and code != 0:
                reason = ('SIGNAL: ' + signal.Signals(-code).name
                          if code < 0 else f'EXIT: {code}')
    payload = None
    if receipt.exists():
        try:
            payload = json.loads(receipt.read_text())
            if not isinstance(payload, dict):
                raise ValueError('receipt is not an object')
        except (ValueError, UnicodeError):
            payload = None
            reason = reason or 'INVALID-RECEIPT'
    elif reason is None:
        reason = 'MISSING-RECEIPT'
    # Preserve native inconclusive evidence; never preserve PASS after a failed process.
    if reason is not None and (payload is None or
            payload.get('search_verdict', payload.get('status')) != 'INCONCLUSIVE'):
        if receipt.exists():
            receipt.with_name(receipt.name + '.incomplete').write_bytes(receipt.read_bytes())
        payload = {'schema': 'xor-interrupted-stage/1', 'status': 'INCONCLUSIVE',
                   'search_verdict': 'INCONCLUSIVE', 'reason': reason,
                   'returncode': code, 'primary_completed': False}
        save_json(receipt, payload)
    return {'receipt': receipt.name, 'guard_status': 'INCONCLUSIVE' if reason else 'COMPLETED',
            'reason': reason, 'returncode': code, 'wall_seconds': time.monotonic() - started,
            'primary_status': payload.get('search_verdict', payload.get('status'))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path, help='read-only archived campaign')
    parser.add_argument('output', type=Path, help='new output directory; existing paths are refused')
    parser.add_argument('--verify-python', default=str(Path(__file__).resolve().parents[2] /
                                                     'physics/.venv/bin/python'))
    parser.add_argument('--wall-seconds', type=float, default=95,
                        help='per-stage outer deadline, at most 95s; raw scripts retain their 60 CPU/90 wall limits')
    args = parser.parse_args()
    if not 0 < args.wall_seconds <= 95:
        parser.error('--wall-seconds must be positive and at most 95')
    root, output = args.run.resolve(), args.output.resolve()
    for name in ('input.json', 'preregistration.json', 'headroom.py', 'verify.py'):
        if not (root / name).is_file():
            parser.error(f'missing archived input: {name}')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('input.json', 'preregistration.json'):
        shutil.copyfile(root / name, output / name)
    environment = dict(os.environ, **{key: '1' for key in THREAD_KEYS})
    result = {'schema': 'xor-replay/1', 'source': str(root), 'status': 'INCONCLUSIVE',
              'stages': [], 'frozen_source_modified': False}
    for script, interpreter, receipt in (
            ('headroom.py', sys.executable, 'result.json'),
            ('verify.py', args.verify_python, 'verification.json')):
        stage = supervise([interpreter, '-I', '-B', str(root / script), str(output)],
                          output / receipt, args.wall_seconds, environment=environment)
        result['stages'].append(stage)
        save_json(output / 'replay.json', result)
        if stage['guard_status'] != 'COMPLETED':
            break
    if (len(result['stages']) == 2
            and result['stages'][-1]['guard_status'] == 'COMPLETED'
            and result['stages'][-1]['primary_status'] == 'PASS'):
        result['status'] = 'COMPLETED'
    save_json(output / 'replay.json', result)
    print(json.dumps(result))
    return 0 if result['status'] == 'COMPLETED' else 2


if __name__ == '__main__':
    raise SystemExit(main())
