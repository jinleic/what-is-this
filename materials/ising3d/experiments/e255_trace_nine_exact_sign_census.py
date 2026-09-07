#!/usr/bin/env python3
"""H682: exact characteristic-zero sign census of the trace-nine norm at rational q>1.

Replays the pinned Wave-27 QQ(q) lift through e254.build_inputs (lift-only data
dies with that frame), specializes the 35 monic rules and F9 at fixed rational
points q0>1, computes all 96 multiplication columns over QQ by the same ordered
reduction, clears every column by its integer denominator LCM, and evaluates the
exact integer determinant with FLINT.  The rational determinant at q0 is the
characteristic-zero norm N(q0)/D(q0) of the fixed rank-96 quotient; every
recorded denominator atom is positive for q>1, so a sign change between two
consecutive points proves a real root of N in that open bracket.  Equal signs
prove nothing.  No polynomial reconstruction, aggregate height, root list,
shifted H0/H1 branch disposition, or thermodynamic claim follows.

Live memory and timings go to stdout only; they never enter the artifact data.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import heapq
import importlib.util
import json
import math
import os
from fractions import Fraction
from pathlib import Path
import sys
import time

from flint import fmpq, fmpz_mat
from sympy.core.cache import clear_cache

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
SCRIPT = Path(__file__).resolve()
E254 = ROOT / 'experiments' / 'e254_trace_nine_bounded_canary.py'
E253_ARTIFACT = ROOT / 'results' / 'spectral' / 'trace_nine_structured_determinant.json'
E251_ARTIFACT = ROOT / 'results' / 'spectral' / 'trace_nine_projection.json'
OUTPUT = ROOT / 'results' / 'spectral' / 'trace_nine_exact_sign_census.json'
SCOPE = 'EXACT_RATIONAL_POINT_SIGN_CENSUS_ONLY'
PRIMARY_PRIME = 2147483647
SECONDARY_PRIME = 2147483629
POINTS = (Fraction(5, 4), Fraction(3, 2), Fraction(5, 3), Fraction(2), Fraction(3), Fraction(5))
# Landed modular witnesses: e251 (q=2 mod 2147483647; q=5/3 mod 2147483629) and the
# e253 direct point checks (q=2, 3, 5/3 mod 2147483647).  The exact determinant reduced
# modulo the prime must reproduce each of them.
PINNED_RESIDUES = {
    ('2', PRIMARY_PRIME): 1660951362,
    ('3', PRIMARY_PRIME): 2143576240,
    ('5/3', PRIMARY_PRIME): 843620109,
    ('5/3', SECONDARY_PRIME): 1100728375,
}
# e251 recorded 26131 ordered reduction steps at q=2 and q=5/3; e253 recorded the same
# generic count over F_p(q).  Specialization can only drop terms, so both exact counts
# are forced to equal 26131.
PINNED_STEPS = {'2': 26131, '5/3': 26131}
REDUCED_MODULAR_NUMERATOR_SHA256 = '45834f272f6e171ae6072008bfd55556d2b3f39f0621c2f57684872999e61686'


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


e254 = load_module('h682_bounded_canary', E254)
base = e254.base
guard = e254.guard
ZERO = fmpq(0)


def nonzero(value: fmpq) -> bool:
    return int(value.p) != 0


def as_pair(value: fmpq) -> tuple[int, int]:
    return int(value.p), int(value.q)


def evaluate_poly(polynomial, point: fmpq) -> fmpq:
    """Horner evaluation of an fmpz_poly (ascending coefficients) at a rational point."""
    value = fmpq(0)
    for coefficient in reversed([int(c) for c in polynomial]):
        value = value * point + coefficient
    return value


def rat_at(value, point: fmpq) -> fmpq:
    denominator = evaluate_poly(value.denominator, point)
    if not nonzero(denominator):
        raise ZeroDivisionError('a lifted denominator vanishes at the census point')
    return evaluate_poly(value.numerator, point) / denominator


def bareiss_det(rows) -> int:
    matrix = [[int(entry) for entry in row] for row in rows]
    size = len(matrix)
    sign, previous = 1, 1
    for k in range(size - 1):
        if matrix[k][k] == 0:
            swap = next((i for i in range(k + 1, size) if matrix[i][k]), None)
            if swap is None:
                return 0
            matrix[k], matrix[swap] = matrix[swap], matrix[k]
            sign = -sign
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                matrix[i][j] = (matrix[i][j] * matrix[k][k] - matrix[i][k] * matrix[k][j]) // previous
        previous = matrix[k][k]
    return sign * matrix[size - 1][size - 1]


def smoke() -> None:
    """Fail before the expensive lift if the pinned python-flint scalar/matrix API drifts."""
    left, right = fmpq(-7, 12), fmpq(5, 18)
    assert as_pair(left * right) == (Fraction(-7, 12) * Fraction(5, 18)).as_integer_ratio()
    assert as_pair(left - right) == (Fraction(-7, 12) - Fraction(5, 18)).as_integer_ratio()
    assert as_pair(left / right) == (Fraction(-7, 12) / Fraction(5, 18)).as_integer_ratio()
    assert as_pair(left * right + 3) == (Fraction(-7, 12) * Fraction(5, 18) + 3).as_integer_ratio()
    assert nonzero(left) and not nonzero(left - left) and as_pair(ZERO) == (0, 1)
    rows = [[(-1) ** ((i * j) % 3) * (3 ** (40 * i + 7 * j) + i - j) for j in range(5)] for i in range(5)]
    assert int(fmpz_mat(rows).det()) == bareiss_det(rows)
    assert int(fmpz_mat([[2, 1], [1, 3]]).det()) == 5
    assert int(fmpz_mat([[1, 2], [2, 4]]).det()) == 0
    polynomial = base.fmpz_poly([3, -2, 0, 5])
    expected = Fraction(3) - 2 * Fraction(5, 3) + 5 * Fraction(5, 3) ** 3
    assert as_pair(evaluate_poly(polynomial, fmpq(5, 3))) == expected.as_integer_ratio()
    assert as_pair(evaluate_poly(base.fmpz_poly([]), fmpq(5, 3))) == (0, 1)
    rational = base.Rat.make(base.fmpz_poly([-(1 << 70), 0, 7]), base.fmpz_poly([2, 0, 3]))
    expected = Fraction(-(1 << 70) + 7 * Fraction(4, 9)) / Fraction(2 + 3 * Fraction(4, 9))
    assert as_pair(rat_at(rational, fmpq(2, 3))) == expected.as_integer_ratio()
    assert residue_mod(Fraction(-7, 12), 101) == (-7 * pow(12, -1, 101)) % 101
    assert residue_mod(Fraction(5, 101), 101) is None
    print('H682_API_SMOKE PASS cases=14', flush=True)


def residue_mod(value: Fraction, prime: int) -> int | None:
    if value.denominator % prime == 0:
        return None
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def report_memory(stage: str, started: float) -> None:
    current = guard.sample_process(os.getpid())
    print('H682_MEMORY ' + json.dumps({
        'stage': stage,
        'cpu_seconds': round(time.process_time() - started, 6),
        'current_resident_bytes': current['resident_bytes'],
        'current_footprint_bytes': current['footprint_bytes'],
        'native_peak_rss_bytes': base.peak_rss_bytes(),
    }, sort_keys=True, separators=(',', ':')), flush=True)


def source_hashes() -> dict[str, str]:
    result = e254.source_hashes()
    for path in (SCRIPT, E253_ARTIFACT, E251_ARTIFACT):
        result[str(path.relative_to(WORKSPACE))] = base.file_sha256(path)
    return result


def reduce_point(polynomial: dict, rules: list) -> tuple[dict, int]:
    """Ordered reduction over QQ with the priority and rule order of the pinned lift."""
    priority = base.mod.reduction_priority
    work = {monomial: coefficient for monomial, coefficient in polynomial.items() if nonzero(coefficient)}
    heap = [priority(monomial) + (monomial,) for monomial in work]
    heapq.heapify(heap)
    queued = set(work)
    steps = 0
    while heap:
        *_, monomial = heapq.heappop(heap)
        queued.discard(monomial)
        coefficient = work.pop(monomial, None)
        if coefficient is None:
            continue
        selected = next(((leading, tail) for leading, tail in rules if base.divides(leading, monomial)), None)
        if selected is None:
            work[monomial] = coefficient
            continue
        leading, tail = selected
        shift = tuple(a - b for a, b in zip(monomial, leading))
        for tail_monomial, tail_coefficient in tail.items():
            output = base.add_monomials(shift, tail_monomial)
            updated = work.get(output, ZERO) - coefficient * tail_coefficient
            if nonzero(updated):
                work[output] = updated
                if output not in queued:
                    heapq.heappush(heap, priority(output) + (output,))
                    queued.add(output)
            else:
                work.pop(output, None)
        steps += 1
    return work, steps


def load_reduced_modular_norm() -> tuple[list[int], list[int]]:
    artifact = json.loads(E253_ARTIFACT.read_text())
    assert artifact['meta']['data_sha256'] == base.canonical_sha256(artifact['data'])
    record = artifact['data']['reduced_modular_norm']
    numerator = [int(value) for value in record['numerator_coefficients_ascending']]
    denominator = [int(value) for value in record['denominator_coefficients_ascending']]
    assert base.canonical_sha256(numerator) == record['numerator_coefficients_sha256'] == REDUCED_MODULAR_NUMERATOR_SHA256
    assert base.canonical_sha256(denominator) == record['denominator_coefficients_sha256']
    assert len(numerator) - 1 == 19846 and len(denominator) - 1 == 16048
    assert artifact['data']['field'] == {'parameter': 'q', 'prime': PRIMARY_PRIME}
    checks = {(check['q'], PRIMARY_PRIME): check['determinant_residue']
              for check in artifact['data']['direct_point_checks']}
    for key, residue in PINNED_RESIDUES.items():
        if key[1] == PRIMARY_PRIME:
            assert checks[key] == residue, f'e253 point check drift at q={key[0]}'
    witnesses = json.loads(E251_ARTIFACT.read_text())['data']['modular_certificates']['physical_witnesses']
    recorded = {(witness['q'], witness['prime']): witness['determinant_residue'] for witness in witnesses}
    assert recorded[('2', PRIMARY_PRIME)] == PINNED_RESIDUES[('2', PRIMARY_PRIME)]
    assert recorded[('5/3', SECONDARY_PRIME)] == PINNED_RESIDUES[('5/3', SECONDARY_PRIME)]
    return numerator, denominator


def horner_mod(coefficients: list[int], point: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * point + coefficient) % prime
    return value


def census_point(q_value: Fraction, rules, physical_f9, standard, modular_norm, started: float) -> dict:
    label = str(q_value)
    report_memory('point_start:' + label, started)
    phase_started = time.process_time()
    point = fmpq(q_value.numerator, q_value.denominator)
    rules_point = []
    for leading, tail in rules:
        specialized = {}
        for monomial, value in tail.items():
            scalar = rat_at(value, point)
            if nonzero(scalar):
                specialized[monomial] = scalar
        rules_point.append((leading, specialized))
    f9_point = {monomial: scalar for monomial, scalar in
                ((monomial, rat_at(value, point)) for monomial, value in physical_f9.items()) if nonzero(scalar)}
    row_index = {monomial: index for index, monomial in enumerate(standard)}
    size = len(standard)
    columns: list[dict[int, fmpq]] = []
    total_steps = 0
    maximum_steps = 0
    for basis_monomial in standard:
        initial = {base.add_monomials(term, basis_monomial): value for term, value in f9_point.items()}
        reduced, steps = reduce_point(initial, rules_point)
        for monomial in reduced:
            if monomial not in row_index:
                raise AssertionError(f'nonstandard matrix monomial {monomial} at q={label}')
        columns.append({row_index[monomial]: value for monomial, value in reduced.items()})
        total_steps += steps
        maximum_steps = max(maximum_steps, steps)
    del rules_point, f9_point
    base.guard_resources(started, f'H682 exact columns at q={label}')
    column_lcms = [math.lcm(*(as_pair(value)[1] for value in column.values())) if column else 1
                   for column in columns]
    cleared = [[0] * size for _ in range(size)]
    for column_index, column in enumerate(columns):
        scale = column_lcms[column_index]
        for row, value in column.items():
            numerator, denominator = as_pair(value)
            cleared[row][column_index] = numerator * (scale // denominator)
    nonzero_entries = sum(1 for row in cleared for entry in row if entry)
    del columns
    cleared_digest = base.canonical_sha256([[str(entry) for entry in row] for row in cleared])
    determinant_started = time.process_time()
    cleared_determinant = int(fmpz_mat(cleared).det())
    determinant_cpu = time.process_time() - determinant_started
    del cleared
    scale_product = math.prod(column_lcms)
    determinant = Fraction(cleared_determinant, scale_product)
    sign = (determinant > 0) - (determinant < 0)
    q_residue = residue_mod(q_value, PRIMARY_PRIME)
    modular_denominator = horner_mod(modular_norm[1], q_residue, PRIMARY_PRIME)
    modular_norm_residue = None
    if modular_denominator:
        modular_norm_residue = (horner_mod(modular_norm[0], q_residue, PRIMARY_PRIME)
                                * pow(modular_denominator, -1, PRIMARY_PRIME)) % PRIMARY_PRIME
    primary_residue = residue_mod(determinant, PRIMARY_PRIME)
    secondary_residue = residue_mod(determinant, SECONDARY_PRIME)
    for prime, residue in ((PRIMARY_PRIME, primary_residue), (SECONDARY_PRIME, secondary_residue)):
        pinned = PINNED_RESIDUES.get((label, prime))
        if pinned is not None and residue != pinned:
            raise AssertionError(f'exact determinant at q={label} misses the landed residue mod {prime}')
    if label in PINNED_STEPS and total_steps != PINNED_STEPS[label]:
        raise AssertionError(f'ordered reduction count at q={label} differs from the landed witness')
    record = {
        'q': label,
        'q_numerator': q_value.numerator,
        'q_denominator': q_value.denominator,
        'ordered_reduction_steps_total': total_steps,
        'ordered_reduction_steps_maximum': maximum_steps,
        'matrix_nonzero_entries': nonzero_entries,
        'column_lcm_bits_maximum': max(value.bit_length() for value in column_lcms),
        'column_lcm_product_bits': scale_product.bit_length(),
        'cleared_integer_matrix_sha256': cleared_digest,
        'cleared_integer_determinant_sign': (cleared_determinant > 0) - (cleared_determinant < 0),
        'cleared_integer_determinant_bits': abs(cleared_determinant).bit_length(),
        'determinant_sign': sign,
        'determinant_numerator_bits': abs(determinant.numerator).bit_length(),
        'determinant_denominator_bits': determinant.denominator.bit_length(),
        'determinant_numerator_decimal': str(determinant.numerator),
        'determinant_denominator_decimal': str(determinant.denominator),
        'determinant_sha256': hashlib.sha256(
            f'{determinant.numerator}/{determinant.denominator}'.encode('ascii')).hexdigest(),
        'residue_mod_2147483647': primary_residue,
        'residue_mod_2147483629': secondary_residue,
        'reduced_modular_norm_residue_mod_2147483647': modular_norm_residue,
        'matches_reduced_modular_norm': (None if modular_norm_residue is None or primary_residue is None
                                         else modular_norm_residue == primary_residue),
    }
    base.guard_resources(started, f'H682 exact determinant at q={label}')
    print('H682_POINT ' + json.dumps({key: value for key, value in record.items()
                                      if not key.endswith('_decimal')},
                                     sort_keys=True, separators=(',', ':')), flush=True)
    print(f'H682_POINT_CPU q={label} total_seconds={time.process_time() - phase_started:.6f} '
          f'determinant_seconds={determinant_cpu:.6f}', flush=True)
    report_memory('point_done:' + label, started)
    return record


def produce_data(started: float, template_certificate: dict, modular_norm) -> dict:
    rules, physical_f9, standard = e254.build_inputs(started, template_certificate)
    clear_cache()
    gc.collect()
    report_memory('lift_after_frame_release', started)
    assert len(standard) == 96
    points = []
    for q_value in POINTS:
        assert q_value > 1
        points.append(census_point(q_value, rules, physical_f9, standard, modular_norm, started))
        gc.collect()
    signs = [record['determinant_sign'] for record in points]
    brackets = [{'lower': left['q'], 'upper': right['q']}
                for left, right in zip(points, points[1:])
                if left['determinant_sign'] * right['determinant_sign'] < 0]
    zero_count = sum(1 for sign in signs if sign == 0)
    return {
        'scope': SCOPE,
        'quotient_rank': len(standard),
        'final_rational_rule_sha256': e254.RULE_DIGEST,
        'point_count': len(points),
        'points': points,
        'summary': {
            'signs': signs,
            'zero_count': zero_count,
            'sign_change_brackets': brackets,
            'reduced_modular_norm_agreement_count': sum(
                1 for record in points if record['matches_reduced_modular_norm'] is True),
            'landed_residue_matches': len(PINNED_RESIDUES),
        },
        'physical_root_status': (
            'EXACT_ROOT_AT_SAMPLED_POINT' if zero_count else
            'REAL_ROOT_BRACKET_PROVED_BY_SIGN_CHANGE' if brackets else
            'NO_SIGN_CHANGE_AMONG_SAMPLED_POINTS_NOT_EVIDENCE_OF_EMPTINESS'),
        'full_polynomial_status': 'NOT_RECONSTRUCTED',
        'interpretation': (
            'Each record is the exact rational determinant of multiplication by F9 on the fixed '
            'rank-96 quotient in the standard-monomial basis at one rational q>1, i.e. the '
            'characteristic-zero norm value N(q)/D(q). Recorded denominator atoms are positive for '
            'q>1, so a sign change between consecutive points proves a real root of the numerator in '
            'the open bracket. Equal signs, and the sampled points themselves, prove nothing about '
            'emptiness, multiplicity, or the shifted H0/H1 physical branch.'),
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
    modular_norm = load_reduced_modular_norm()
    sources = source_hashes()
    data = produce_data(started, inherited['data']['template'], modular_norm)
    gc.collect()
    report_memory('all_census_data_released', started)
    assert source_hashes() == sources, 'source changed during exact computation'
    base.guard_resources(started, 'H682 before artifact write')
    artifact = {'meta': {'schema_version': 1, 'source_sha256': sources,
                         'data_sha256': base.canonical_sha256(data)}, 'data': data}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as stream:
        stream.write(base.canonical_bytes(artifact))
    print('H682_SIGN_CENSUS PASS artifact=' + str(args.output), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
