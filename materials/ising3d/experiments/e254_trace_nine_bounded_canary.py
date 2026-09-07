#!/usr/bin/env python3
"""H681: two exact QQ(q) columns with released lift data and streamed hashing.

A successful producer is not resource acceptance: the external guard and
independent verifier must both exit successfully before a wider sweep.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import lzma
import math
import os
from pathlib import Path
import sys
import time

from sympy.core.cache import clear_cache

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
SCRIPT = Path(__file__).resolve()
VERIFIER = ROOT / 'tests' / 'test_trace_nine_bounded_canary.py'
GUARD = ROOT / 'tools' / 'resource_guard.py'
OUTPUT = ROOT / 'results' / 'spectral' / 'trace_nine_bounded_canary.json'
RULE_DIGEST = 'f4cfe038d7625c6a411f2e01b9af12d9ae9ca06d980e64ad38d8eff030fe22cc'
COLUMN_DIGESTS = {
    (0, 0, 0, 0, 0): '69065ff900271f44f3f34856b1bad85fd148182b8a045aefd1cc66cb1b8bdfb9',
    (0, 0, 0, 0, 8): '63b8a477b5d9dd25a354d1bd2f7d18666cbcf1bf21cb6060d0c2b1b51441f47b',
}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


base = load_module('h681_wave27_producer', ROOT / 'experiments' / 'e252_trace_nine_norm_envelope.py')
guard = load_module('h681_resource_guard', GUARD)


def hash_integer_array(digest, values) -> None:
    digest.update(b'[')
    for index, value in enumerate(values):
        if index:
            digest.update(b',')
        digest.update(str(value).encode('ascii'))
    digest.update(b']')


def stream_column_sha256(column) -> str:
    """Hash the existing canonical nested JSON without copying its coefficient vector."""
    digest = hashlib.sha256()
    digest.update(b'[')
    for index, monomial in enumerate(sorted(column)):
        if index:
            digest.update(b',')
        value = column[monomial]
        digest.update(b'[')
        hash_integer_array(digest, monomial)
        digest.update(b',[')
        hash_integer_array(digest, value.numerator)
        digest.update(b',')
        hash_integer_array(digest, value.denominator)
        digest.update(b']]')
    digest.update(b']')
    return digest.hexdigest()


def smoke() -> None:
    # Distinct canonical byte boundaries: empty, zero polynomial, negative and
    # large integers, nonconstant denominators, and differing insertion order.
    rational = base.Rat.make(base.fmpz_poly([-(1 << 1024), 0, 7]), base.fmpz_poly([2, 0, 3]))
    cases = [
        {},
        {(0, 0, 0, 0, 0): base.Rat.scalar(0)},
        {(0, 0, 0, 0, 1): rational, (0, 0, 0, 0, 0): base.Rat.scalar(-17)},
    ]
    for column in cases:
        canonical = [[list(term), [list(map(int, value.numerator)),
                                   list(map(int, value.denominator))]]
                     for term, value in sorted(column.items())]
        assert stream_column_sha256(column) == base.canonical_sha256(canonical)
        assert stream_column_sha256(dict(reversed(list(column.items())))) == stream_column_sha256(column)
    print('H681_STREAMING_HASH_SMOKE PASS cases=3', flush=True)


def report_memory(stage: str, started: float) -> None:
    current = guard.sample_process(os.getpid())
    print('H681_MEMORY ' + json.dumps({
        'stage': stage,
        'cpu_seconds': round(time.process_time() - started, 6),
        'current_resident_bytes': current['resident_bytes'],
        'current_footprint_bytes': current['footprint_bytes'],
        'native_peak_rss_bytes': base.peak_rss_bytes(),
    }, sort_keys=True, separators=(',', ':')), flush=True)


def source_hashes() -> dict[str, str]:
    result = base.source_hashes()
    for path in (SCRIPT, VERIFIER, GUARD, base.TEMPLATE, base.OUTPUT,
                 ROOT / 'experiments' / 'e248_replica_trace_nine.py',
                 ROOT / 'tests' / 'test_replica_trace_nine.py'):
        result[str(path.relative_to(WORKSPACE))] = base.file_sha256(path)
    return result


def load_pinned_template(certificate: dict) -> dict:
    """Reuse certified bytes, not a fresh preset-9 recompression of identical data."""
    compressed = base.TEMPLATE.read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == certificate['compressed_sha256']
    raw = lzma.decompress(compressed, format=lzma.FORMAT_XZ)
    assert hashlib.sha256(raw).hexdigest() == certificate['payload_sha256']
    payload = json.loads(raw)
    assert payload['schema'] == 'ising3d.trace-nine-lift-template/v1'
    return payload


def build_inputs(started: float, template_certificate: dict):
    """Lift-only objects die with this frame; no retained runpy namespace."""
    template = load_pinned_template(template_certificate)
    equations, _ = base.mod.abstract_mode_equations()
    trace = json.loads(base.TRACE_ARTIFACT.read_text())
    trace_rows = [[int(value) for value in row]
                  for row in trace['data']['trace_coefficients_ascending']]
    targets = tuple(base.Rat.from_sympy(value)
                    for value in base.mod.normalized_targets(trace_rows).values())
    physical = [{monomial: base.evaluate_target_polynomial(terms, targets)
                 for monomial, terms in row}
                for row in base.mod.compile_equations(equations)]
    base.guard_resources(started, 'H681 exact targets')
    remainders = [
        {monomial: value for monomial, value in row.items()
         if sum(monomial) < equation.total_degree()}
        for equation, row in zip(equations[:5], physical[:5])
    ]
    lifted = []
    zero, one = base.Rat.scalar(0), base.Rat.scalar(1)
    for record in template['records']:
        polynomial = {monomial: base.Rat.scalar(value)
                      for monomial, value in base.parse_q_polynomial(record['basis']).items()}
        for source, entries in zip(remainders, record['transform']):
            for left, scalar in base.parse_q_polynomial(entries).items():
                for right, coefficient in source.items():
                    monomial = base.add_monomials(left, right)
                    value = polynomial.get(monomial, zero) + coefficient * scalar
                    if value:
                        polynomial[monomial] = value
                    else:
                        polynomial.pop(monomial, None)
        leading = tuple(record['leading_monomial'])
        assert polynomial[leading] == one
        lifted.append((leading, {term: value for term, value in polynomial.items()
                                 if term != leading}))
        base.guard_resources(started, 'H681 exact lift')
    processed = []
    by_index = {}
    for index in sorted(range(len(lifted)), key=lambda index: (sum(lifted[index][0]), lifted[index][0])):
        leading, tail = lifted[index]
        reduced, _ = base.reduce_exact(tail, processed)
        by_index[index] = (leading, reduced)
        processed.append((leading, reduced))
        base.guard_resources(started, 'H681 exact interreduction')
    rules = [by_index[index] for index in range(len(lifted))]
    leaders = [leading for leading, _ in rules]
    _, standard = base.standard_monomials(leaders)
    assert len(standard) == 96
    assert base.canonical_sha256([stream_column_sha256(tail) for _, tail in rules]) == RULE_DIGEST
    for row in physical[:5]:
        assert not base.reduce_exact(row, rules)[0]
    report_memory('lift_before_frame_release', started)
    return rules, physical[5], standard


def compute_column(monomial, physical_f9, rules, started: float) -> dict:
    report_memory('column_start:' + ','.join(map(str, monomial)), started)
    phase_started = time.process_time()
    initial = {base.add_monomials(term, monomial): value
               for term, value in physical_f9.items()}
    column, steps = base.reduce_exact(initial, rules)
    assert all(not any(base.divides(leader, term) for leader, _ in rules) for term in column)
    scalar_lcm = math.lcm(*(abs(int(value.denominator.content())) for value in column.values()))
    record = {
        'basis_monomial': list(monomial),
        'entries': len(column),
        'ordered_reduction_steps': steps,
        'maximum_numerator_degree': max(value.numerator.degree() for value in column.values()),
        'maximum_denominator_degree': max(value.denominator.degree() for value in column.values()),
        'maximum_numerator_l1_log2': max(base.ceil_log2(sum(abs(int(c)) for c in value.numerator))
                                       for value in column.values()),
        'column_scalar_lcm_bits': scalar_lcm.bit_length(),
        'column_sha256': stream_column_sha256(column),
    }
    assert record['column_sha256'] == COLUMN_DIGESTS[monomial]
    base.guard_resources(started, 'H681 complete streamed column')
    print('H681_COLUMN ' + json.dumps(record, sort_keys=True, separators=(',', ':')), flush=True)
    print(f'H681_COLUMN_CPU seconds={time.process_time() - phase_started:.6f}', flush=True)
    report_memory('column_before_frame_release', started)
    return record


def produce_data(started: float, template_certificate: dict) -> dict:
    rules, physical_f9, standard = build_inputs(started, template_certificate)
    clear_cache()
    gc.collect()
    report_memory('lift_after_frame_release', started)
    largest = max(standard, key=lambda term: (sum(term), term))
    assert standard[0] == (0, 0, 0, 0, 0)
    assert largest == (0, 0, 0, 0, 8)
    columns = []
    for monomial in (standard[0], largest):
        columns.append(compute_column(monomial, physical_f9, rules, started))
        gc.collect()
        report_memory('column_after_frame_release', started)
    return {
        'scope': 'BOUNDED_TWO_COLUMN_CANARY_ONLY',
        'quotient_rank': len(standard),
        'final_rational_rule_sha256': RULE_DIGEST,
        'columns': columns,
        'full_matrix_status': 'NOT_COMPUTED',
        'physical_root_status': 'UNRESOLVED',
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    smoke()
    if args.smoke:
        report_memory('streaming_smoke', time.process_time())
        return 0
    if args.output.exists():
        raise FileExistsError(f'refusing to overwrite existing evidence: {args.output}')
    started = time.process_time()
    base.runtime_manifest()
    inherited = json.loads(base.OUTPUT.read_text())
    assert inherited['meta']['source_sha256'] == base.source_hashes()
    assert inherited['meta']['data_sha256'] == base.canonical_sha256(inherited['data'])
    assert all(check['passed'] for check in inherited['checks'])
    sources = source_hashes()
    data = produce_data(started, inherited['data']['template'])
    gc.collect()
    report_memory('all_canary_data_released', started)
    assert source_hashes() == sources, 'source changed during exact reconstruction'
    base.guard_resources(started, 'H681 before artifact write')
    artifact = {'meta': {'schema_version': 1, 'source_sha256': sources,
                         'data_sha256': base.canonical_sha256(data)}, 'data': data}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as stream:
        stream.write(base.canonical_bytes(artifact))
    print('H681_BOUNDED_CANARY PASS artifact=' + str(args.output), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
