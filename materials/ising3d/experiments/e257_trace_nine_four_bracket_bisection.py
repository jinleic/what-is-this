#!/usr/bin/env python3
"""H688: exact bisection of all four H684 real-root brackets of the trace-nine norm.

The confirmed census (e255, results/spectral/trace_nine_exact_sign_census.json,
data digest 009d8d7e...a55e9) records exact determinant signs -1,+1,-1,+1,-1,-1
at q=5/4,3/2,5/3,2,3,5, so the characteristic-zero norm numerator N has at least
one real root in each of the FOUR open brackets (5/4,3/2), (3/2,5/3), (5/3,2),
(2,3).  The rejected e256 hard-coded only three of them; this producer derives
the bracket list from the census signs and asserts it against the artifact's
recorded sign_change_brackets, so a miscount fails before any heavy work.

It replays the pinned lift once through e254.build_inputs, reuses
e255.census_point unchanged (scalar specialization of the 35 monic rules and
F9, 96 columns over QQ, integer column clearing, one FLINT integer
determinant), and bisects every bracket DEPTH=3 times at exact rational
midpoints (twelve points).  A midpoint sign opposite to the lower endpoint keeps
the lower half, otherwise the upper half; an exact zero stops that bracket with
an exact rational root.  Depth 3 rather than 4 keeps the PREFLIGHT CPU estimate
(lift plus twelve points) well inside the 1800 CPU s self-limit; deeper
refinement is a separate later task on the refined brackets.

What a refined bracket proves: at least one real root of N lies in the final
open interval.  It does not prove uniqueness, multiplicity, emptiness outside
the sampled intervals, the shifted H0/H1 physical branch, any polynomial
coefficient, aggregate height, or thermodynamic behaviour.  Live memory and
timings go to stdout only; they never enter the artifact data.
"""
from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import os
from fractions import Fraction
from pathlib import Path
import sys
import time

from sympy.core.cache import clear_cache

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
SCRIPT = Path(__file__).resolve()
E255 = ROOT / 'experiments' / 'e255_trace_nine_exact_sign_census.py'
CENSUS_ARTIFACT = ROOT / 'results' / 'spectral' / 'trace_nine_exact_sign_census.json'
OUTPUT = ROOT / 'results' / 'spectral' / 'trace_nine_four_bracket_bisection.json'
SCOPE = 'EXACT_RATIONAL_BRACKET_BISECTION_ONLY'
# Confirmed H684 census data digest (supervisor evaluation 0002_H684).
CENSUS_DATA_SHA256 = '009d8d7e647184a0cc13052bcbc63eb477eacc0fe86d2c121645ec675b6a55e9'
CENSUS_LABELS = ('5/4', '3/2', '5/3', '2', '3', '5')
CENSUS_SIGNS = (-1, 1, -1, 1, -1, -1)
# All consecutive sign changes of CENSUS_SIGNS: four brackets, including (2,3).
EXPECTED_BRACKETS = (('5/4', '3/2'), ('3/2', '5/3'), ('5/3', '2'), ('2', '3'))
DEPTH = 3


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


e255 = load_module('h688_sign_census', E255)
e254 = e255.e254
base = e255.base
guard = e255.guard


def report_memory(stage: str, started: float) -> None:
    current = guard.sample_process(os.getpid())
    print('H688_MEMORY ' + json.dumps({
        'stage': stage,
        'cpu_seconds': round(time.process_time() - started, 6),
        'current_resident_bytes': current['resident_bytes'],
        'current_footprint_bytes': current['footprint_bytes'],
        'native_peak_rss_bytes': base.peak_rss_bytes(),
    }, sort_keys=True, separators=(',', ':')), flush=True)


def source_hashes() -> dict[str, str]:
    result = e255.source_hashes()
    for path in (SCRIPT, CENSUS_ARTIFACT):
        result[str(path.relative_to(WORKSPACE))] = base.file_sha256(path)
    return result


def brackets_from_signs(labels, signs) -> tuple[tuple[str, str], ...]:
    """Every consecutive pair of sampled points with strictly opposite signs."""
    return tuple((left, right) for (left, lsign), (right, rsign) in zip(zip(labels, signs), zip(labels[1:], signs[1:]))
                 if lsign * rsign < 0)


def bisect_bracket(lower: Fraction, upper: Fraction, lower_sign: int, upper_sign: int,
                   depth: int, evaluate) -> dict:
    """Deterministic exact bisection; evaluate(Fraction) -> record with 'determinant_sign'."""
    assert lower < upper and lower_sign * upper_sign < 0
    initial_lower, initial_upper = lower, upper
    steps = []
    exact_root = None
    for _ in range(depth):
        midpoint = (lower + upper) / 2
        assert lower < midpoint < upper
        record = evaluate(midpoint)
        sign = record['determinant_sign']
        if sign == 0:
            exact_root = midpoint
            steps.append({'q': str(midpoint), 'determinant_sign': 0, 'replaced': 'exact_root'})
            break
        if sign * lower_sign < 0:
            upper, upper_sign = midpoint, sign
            steps.append({'q': str(midpoint), 'determinant_sign': sign, 'replaced': 'upper'})
        else:
            lower, lower_sign = midpoint, sign
            steps.append({'q': str(midpoint), 'determinant_sign': sign, 'replaced': 'lower'})
    assert exact_root is not None or lower_sign * upper_sign < 0
    return {
        'initial_lower': str(initial_lower),
        'initial_upper': str(initial_upper),
        'initial_width': str(initial_upper - initial_lower),
        'steps': steps,
        'final_lower': str(lower),
        'final_upper': str(upper),
        'final_lower_sign': lower_sign,
        'final_upper_sign': upper_sign,
        'final_width': str(upper - lower),
        'width_ratio': str((initial_upper - initial_lower) / (upper - lower)),
        'exact_root': None if exact_root is None else str(exact_root),
        'root_status': ('EXACT_RATIONAL_ROOT_FOUND' if exact_root is not None
                        else 'AT_LEAST_ONE_REAL_ROOT_IN_FINAL_OPEN_BRACKET'),
    }


def smoke() -> None:
    e255.smoke()
    assert brackets_from_signs(CENSUS_LABELS, CENSUS_SIGNS) == EXPECTED_BRACKETS
    assert brackets_from_signs(('1', '2', '3'), (1, 1, -1)) == (('2', '3'),)
    assert brackets_from_signs(('1', '2'), (1, 1)) == ()
    target = Fraction(137, 100)

    def synthetic(point: Fraction) -> dict:
        value = point - target
        return {'determinant_sign': (value > 0) - (value < 0)}

    result = bisect_bracket(Fraction(5, 4), Fraction(3, 2), -1, 1, 3, synthetic)
    assert Fraction(result['final_lower']) < target < Fraction(result['final_upper'])
    assert Fraction(result['final_width']) == Fraction(1, 32) and result['width_ratio'] == '8'
    # 11/8 > 1.37 replaces the upper end; 21/16, 43/32 < 1.37 replace the lower end.
    assert [step['q'] for step in result['steps']] == ['11/8', '21/16', '43/32']
    assert [step['replaced'] for step in result['steps']] == ['upper', 'lower', 'lower']
    assert (result['final_lower'], result['final_upper']) == ('43/32', '11/8')
    assert result['exact_root'] is None and len(result['steps']) == 3
    flipped = bisect_bracket(Fraction(5, 4), Fraction(3, 2), 1, -1, 3, lambda p: {
        'determinant_sign': -synthetic(p)['determinant_sign']})
    assert flipped['final_lower'] == result['final_lower'] and flipped['final_upper'] == result['final_upper']
    wide = bisect_bracket(Fraction(2), Fraction(3), 1, -1, 3, lambda p: {
        'determinant_sign': (p < Fraction(23, 10)) - (p > Fraction(23, 10))})
    # 5/2 > 2.3 replaces the upper end; 9/4 < 2.3 replaces the lower end; 19/8 > 2.3 replaces the upper end.
    assert [step['q'] for step in wide['steps']] == ['5/2', '9/4', '19/8']
    assert (wide['final_lower'], wide['final_upper'], wide['final_width']) == ('9/4', '19/8', '1/8')
    exact = bisect_bracket(Fraction(1), Fraction(2), -1, 1, 3, lambda p: {
        'determinant_sign': (p > Fraction(3, 2)) - (p < Fraction(3, 2))})
    assert exact['exact_root'] == '3/2' and len(exact['steps']) == 1
    assert exact['root_status'] == 'EXACT_RATIONAL_ROOT_FOUND'
    print('H688_BISECTION_SMOKE PASS cases=7', flush=True)


def load_census() -> dict:
    artifact = json.loads(CENSUS_ARTIFACT.read_text())
    data = artifact['data']
    assert artifact['meta']['schema_version'] == 1
    assert artifact['meta']['data_sha256'] == base.canonical_sha256(data) == CENSUS_DATA_SHA256
    assert artifact['meta']['source_sha256'] == e255.source_hashes(), 'census provenance drift'
    assert data['scope'] == e255.SCOPE and data['quotient_rank'] == 96
    assert data['final_rational_rule_sha256'] == e254.RULE_DIGEST
    assert data['point_count'] == 6 and data['summary']['zero_count'] == 0
    labels = tuple(record['q'] for record in data['points'])
    signs = tuple(record['determinant_sign'] for record in data['points'])
    assert labels == CENSUS_LABELS and signs == CENSUS_SIGNS
    assert tuple(str(Fraction(label)) for label in labels) == labels
    assert data['summary']['signs'] == list(CENSUS_SIGNS)
    recorded = tuple((entry['lower'], entry['upper']) for entry in data['summary']['sign_change_brackets'])
    derived = brackets_from_signs(labels, signs)
    assert recorded == derived == EXPECTED_BRACKETS, (recorded, derived)
    sign_of = dict(zip(labels, signs))
    return {'labels': labels, 'signs': signs, 'sign_of': sign_of, 'brackets': derived}


def produce_data(started: float, template_certificate: dict, modular_norm, census: dict) -> dict:
    rules, physical_f9, standard = e254.build_inputs(started, template_certificate)
    clear_cache()
    gc.collect()
    report_memory('lift_after_frame_release', started)
    assert len(standard) == 96
    points: list[dict] = []

    def evaluate(point: Fraction) -> dict:
        assert point > 1
        record = e255.census_point(point, rules, physical_f9, standard, modular_norm, started)
        points.append(record)
        gc.collect()
        return record

    input_brackets = census['brackets']
    brackets = []
    for index, (lower_label, upper_label) in enumerate(input_brackets):
        lower, upper = Fraction(lower_label), Fraction(upper_label)
        result = bisect_bracket(lower, upper, census['sign_of'][lower_label],
                                census['sign_of'][upper_label], DEPTH, evaluate)
        result = {'index': index, 'initial_lower_sign': census['sign_of'][lower_label],
                  'initial_upper_sign': census['sign_of'][upper_label], **result}
        brackets.append(result)
        print('H688_BRACKET ' + json.dumps(result, sort_keys=True, separators=(',', ':')), flush=True)
    zero_count = sum(1 for record in points if record['determinant_sign'] == 0)
    agreement = sum(1 for record in points if record['matches_reduced_modular_norm'] is True)
    return {
        'scope': SCOPE,
        'quotient_rank': len(standard),
        'final_rational_rule_sha256': e254.RULE_DIGEST,
        'source_census_artifact': str(CENSUS_ARTIFACT.relative_to(WORKSPACE)),
        'source_census_data_sha256': CENSUS_DATA_SHA256,
        'input_labels': list(census['labels']),
        'input_signs': list(census['signs']),
        'input_bracket_count': len(input_brackets),
        'input_brackets': [{'lower': lower, 'upper': upper, 'lower_sign': census['sign_of'][lower],
                            'upper_sign': census['sign_of'][upper]} for lower, upper in input_brackets],
        'bisection_depth': DEPTH,
        'planned_point_count': DEPTH * len(input_brackets),
        'computed_point_count': len(points),
        'points': points,
        'brackets': brackets,
        'summary': {
            'zero_count': zero_count,
            'reduced_modular_norm_agreement_count': agreement,
            'all_points_match_reduced_modular_norm': agreement == len(points),
            'final_brackets': [{'lower': entry['final_lower'], 'upper': entry['final_upper'],
                                'width': entry['final_width']} for entry in brackets],
            'exact_roots': [entry['exact_root'] for entry in brackets if entry['exact_root'] is not None],
        },
        'physical_root_status': (
            'EXACT_RATIONAL_ROOT_AT_BISECTION_POINT' if zero_count else
            'FOUR_REFINED_REAL_ROOT_BRACKETS_NO_UNIQUENESS_OR_BRANCH_CLAIM'),
        'full_polynomial_status': 'NOT_RECONSTRUCTED',
        'interpretation': (
            'Each point record is the exact rational determinant of multiplication by F9 on the fixed '
            'rank-96 quotient at one rational q>1, computed by the unchanged e255 specialization; every '
            'recorded e252 denominator atom is positive for q>1, so the sign is the sign of the norm '
            'numerator N there. Each final open bracket contains at least one real root of N. Nothing '
            'here proves uniqueness, multiplicity, emptiness elsewhere, the shifted H0/H1 physical '
            'branch, any characteristic-zero coefficient, an aggregate height, or thermodynamic behaviour.'),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    smoke()
    if args.smoke:
        report_memory('api_smoke', time.process_time())
        return 0
    if args.output.exists():
        raise FileExistsError(f'refusing to overwrite existing evidence: {args.output}')
    started = time.process_time()
    base.runtime_manifest()
    inherited = json.loads(base.OUTPUT.read_text())
    assert inherited['meta']['source_sha256'] == base.source_hashes()
    assert inherited['meta']['data_sha256'] == base.canonical_sha256(inherited['data'])
    assert all(check['passed'] for check in inherited['checks'])
    modular_norm = e255.load_reduced_modular_norm()
    census = load_census()
    sources = source_hashes()
    data = produce_data(started, inherited['data']['template'], modular_norm, census)
    gc.collect()
    report_memory('all_bisection_data_released', started)
    assert source_hashes() == sources, 'source changed during exact computation'
    base.guard_resources(started, 'H688 before artifact write')
    artifact = {'meta': {'schema_version': 1, 'source_sha256': sources,
                         'data_sha256': base.canonical_sha256(data)}, 'data': data}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as stream:
        stream.write(base.canonical_bytes(artifact))
    print('H688_FOUR_BRACKET_BISECTION PASS artifact=' + str(args.output), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
