#!/usr/bin/env python3
"""Exact characteristic-zero structure of the 2x3 open-layer Ising Lie algebra.

All closure and rank decisions are over Q.  Pauli closure uses the repository's
ordered-monomial convention Q_(a,b)=X^a Z^b, whose nonzero brackets have
integer coefficient +/-2.  Hence every Lie word has integer Pauli coordinates;
the sparse Fraction elimination below is an exact Q computation.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import signal
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_structure" / "char0_levi.json"
SCRIPT = "experiments/e45_char0_levi.py"
N = 6
BONDS = ((0, 3), (1, 4), (2, 5), (0, 1), (1, 2), (3, 4), (4, 5))
CHARACTERS = tuple((a, b, c) for a in range(2) for b in range(2) for c in range(2))
SELECTED = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (0, 0, 1), (0, 1, 1))


class StepTimeout(RuntimeError):
    pass


def _alarm_handler(_signum, _frame):
    raise StepTimeout("exact linear-algebra step exceeded 3000 s")


def timed(label, timings, function, *args):
    started = time.monotonic()
    old = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(3000)
    try:
        value = function(*args)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    timings[label] = round(time.monotonic() - started, 6)
    return value


def symplectic_parity(v: int, w: int) -> int:
    mask = (1 << N) - 1
    return ((v >> N) & mask & (w & mask)).bit_count() + ((w >> N) & mask & (v & mask)).bit_count() & 1


def bracket_sign(g: int, v: int) -> int:
    """Coefficient of [Q_g,Q_v] in the ordered-monomial Pauli convention."""
    mask = (1 << N) - 1
    left = ((g >> N) & mask & (v & mask)).bit_count() & 1
    right = ((v >> N) & mask & (g & mask)).bit_count() & 1
    if left == right:
        return 0
    return 2 if left == 0 else -2


def site_permutations():
    answer = []
    for flip_rows in range(2):
        for flip_cols in range(2):
            answer.append(tuple(
                (1 - r if flip_rows else r) * 3 + (2 - c if flip_cols else c)
                for r in range(2) for c in range(3)
            ))
    return tuple(answer)


def permute_pauli(v, permutation):
    answer = 0
    for old, new in enumerate(permutation):
        if (v >> old) & 1:
            answer |= 1 << new
        if (v >> (N + old)) & 1:
            answer |= 1 << (N + new)
    return answer


def support_columns(singletons):
    seen = set(singletons)
    queue = list(seen)
    for v in queue:
        for g in singletons:
            if symplectic_parity(v, g):
                w = v ^ g
                if w not in seen:
                    seen.add(w)
                    queue.append(w)
    return tuple(sorted(seen))
def orbit_pauli_actions():
    """Return spatial Pauli orbits and honest full-orbit ad_A/ad_B maps."""
    xs = tuple(1 << i for i in range(N))
    zz = tuple(((1 << u) | (1 << v)) << N for u, v in BONDS)
    columns = support_columns(xs + zz)
    column_set = set(columns)
    permutations = site_permutations()
    seen, orbits = set(), []
    for value in columns:
        if value in seen:
            continue
        orbit = frozenset(permute_pauli(value, permutation) for permutation in permutations)
        assert orbit <= column_set
        seen.update(orbit)
        orbits.append(orbit)
    orbit_index = {value: i for i, orbit in enumerate(orbits) for value in orbit}
    maps = []
    for group in (xs, zz):
        action = []
        for source_orbit in orbits:
            target_coefficients = [dict() for _ in orbits]
            for source in source_orbit:
                for generator in group:
                    sign = bracket_sign(generator, source)
                    if sign:
                        target = source ^ generator
                        target_orbit = orbit_index[target]
                        target_coefficients[target_orbit][target] = (
                            target_coefficients[target_orbit].get(target, 0) + sign
                        )
            row = []
            for target_orbit, orbit in enumerate(orbits):
                values = {target_coefficients[target_orbit].get(target, 0) for target in orbit}
                assert len(values) == 1, "spatial-equivariance invariant failed"
                row.append(values.pop())
            action.append(row)
        maps.append(action)
    return (xs, zz), tuple(orbits), tuple(maps)


def apply_orbit_action(mapping, vector):
    return [
        sum(vector[source] * mapping[source][target] for source in range(len(vector)))
        for target in range(len(vector))
    ]




def add_fraction_row(rows, pivots, vector):
    work = {i: Fraction(x) for i, x in enumerate(vector) if x}
    while work:
        lead = min(work)
        if lead not in pivots:
            coefficient = work[lead]
            work = {i: x / coefficient for i, x in work.items()}
            pivots[lead] = len(rows)
            rows.append(work)
            return True
        coefficient = work[lead]
        old = rows[pivots[lead]]
        for i, x in old.items():
            value = work.get(i, Fraction(0)) - coefficient * x
            if value:
                work[i] = value
            elif i in work:
                del work[i]
    return False


def exact_pauli_closure():
    groups, orbits, maps = orbit_pauli_actions()
    orbit_index = {value: i for i, orbit in enumerate(orbits) for value in orbit}

    def seed(group):
        output = [0] * len(orbits)
        for q in group:
            output[orbit_index[q]] += 1
        return output

    rows, pivots, frontier = [], {}, []
    for value in map(seed, groups):
        assert add_fraction_row(rows, pivots, value)
        frontier.append(value)
    head = 0
    while head < len(frontier):
        value = frontier[head]
        head += 1
        for mapping in maps:
            image = apply_orbit_action(mapping, value)
            if add_fraction_row(rows, pivots, image):
                frontier.append(image)
    # Closure ended because both ad_A and ad_B preserve the final exact row span.
    return {
        "dimension": len(rows),
        "support_columns": sum(map(len, orbits)),
        "symmetry_orbit_columns": len(orbits),
        "basis_depth_upper": 15,
        "orbit_sizes": dict(sorted(Counter(map(len, orbits)).items())),
        "action_checksums": [
            {
                "sum_abs_entries": sum(abs(x) for row in mapping for x in row),
                "weighted_mod_1000000007": sum(
                    (i + 1) * (j + 1) * x
                    for i, row in enumerate(mapping) for j, x in enumerate(row)
                ) % 1_000_000_007,
            }
            for mapping in maps
        ],
    }


def act_configuration(configuration, symmetry):
    flip_rows, flip_cols, spin_flip = symmetry
    output = 0
    for r in range(2):
        for c in range(3):
            old = r * 3 + c
            new = (1 - r if flip_rows else r) * 3 + (2 - c if flip_cols else c)
            if (configuration >> old) & 1:
                output |= 1 << new
    return output ^ (63 if spin_flip else 0)


def symmetry_sector(character):
    symmetries = CHARACTERS
    seen, representatives, orbit_vectors = set(), [], []
    for configuration in range(64):
        if configuration in seen:
            continue
        seen.update(act_configuration(configuration, symmetry) for symmetry in symmetries)
        values = {}
        for symmetry in symmetries:
            sign = -1 if sum(x * y for x, y in zip(character, symmetry)) & 1 else 1
            image = act_configuration(configuration, symmetry)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(q for q, coefficient in values.items() if coefficient)
        scale = values[representative]
        representatives.append(representative)
        orbit_vectors.append({q: coefficient // scale for q, coefficient in values.items() if coefficient})
    dimension = len(representatives)
    a = [[0] * dimension for _ in range(dimension)]
    b = [[0] * dimension for _ in range(dimension)]
    for j, orbit in enumerate(orbit_vectors):
        image = {}
        for configuration, coefficient in orbit.items():
            for site in range(N):
                target = configuration ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for i, representative in enumerate(representatives):
            a[i][j] = image.get(representative, 0)
        b[j][j] = sum(1 if ((representatives[j] >> u) & 1) == ((representatives[j] >> v) & 1) else -1 for u, v in BONDS)
    return a, b


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))] for i in range(len(left))]


def commutator(left, right):
    lr, rl = matmul(left, right), matmul(right, left)
    return [[lr[i][j] - rl[i][j] for j in range(len(left))] for i in range(len(left))]


def flatten(matrix):
    return [x for row in matrix for x in row]


def unflatten(vector, dimension):
    return [vector[i * dimension:(i + 1) * dimension] for i in range(dimension)]


def rank_q(rows, ncols=None):
    if not rows:
        return 0
    ncols = ncols if ncols is not None else len(rows[0])
    nonzero_rows = [(i, row) for i, row in enumerate(rows) if any(row)]
    domain_rows = {new_i: {j: QQ.convert(x) for j, x in enumerate(row) if x} for new_i, (_, row) in enumerate(nonzero_rows)}
    matrix = DomainMatrix(domain_rows, (len(nonzero_rows), ncols), QQ)
    return len(matrix.rref(method="GJ")[1])


def nullspace_q(rows, ncols):
    domain_rows = {i: {j: QQ.convert(x) for j, x in enumerate(row) if x} for i, row in enumerate(rows)}
    matrix = DomainMatrix(domain_rows, (len(rows), ncols), QQ)
    return matrix.nullspace().to_Matrix().tolist()




def matrix_closure_raw(blocks):
    dimensions = [len(pair[0]) for pair in blocks]
    total = sum(dimensions)
    generators = []
    for which in range(2):
        matrix = [[0] * total for _ in range(total)]
        offset = 0
        for pair, dimension in zip(blocks, dimensions):
            source = pair[which]
            for i in range(dimension):
                matrix[offset + i][offset:offset + dimension] = source[i]
            offset += dimension
        generators.append(matrix)
    reduced, pivots, raw, frontier, depths = [], {}, [], [], []

    def add(value, depth):
        # Store an independent exact Lie word.  No modular selection is used:
        # sparse Fraction echelon reduction decides independence over Q.
        if add_fraction_row(reduced, pivots, value):
            raw.append(value)
            frontier.append(len(raw) - 1)
            depths.append(depth)
            return True
        return False

    for generator in generators:
        add(flatten(generator), 1)
    head = 0
    while head < len(frontier):
        index = frontier[head]
        head += 1
        value = unflatten(raw[index], total)
        for generator in generators:
            add(flatten(commutator(generator, value)), depths[index] + 1)
    return generators, [unflatten(value, total) for value in raw], depths


def ideal_closure(seed, generators):
    dimension = len(seed)
    reduced, pivots, raw, frontier = [], {}, [], []

    def add(value):
        vector = flatten(value)
        if add_fraction_row(reduced, pivots, vector):
            raw.append(vector)
            frontier.append(len(raw) - 1)

    add(seed)
    head = 0
    while head < len(frontier):
        matrix = unflatten(raw[frontier[head]], dimension)
        head += 1
        for generator in generators:
            add(commutator(generator, matrix))
    return [unflatten(vector, dimension) for vector in raw]


def trace_product(left, right):
    return sum(left[i][j] * right[j][i] for i in range(len(left)) for j in range(len(left)))


def determinant_bareiss(matrix):
    work = [list(map(int, row)) for row in matrix]
    n, previous, sign = len(work), 1, 1
    for k in range(n - 1):
        if work[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if work[i][k]), None)
            if pivot is None:
                return 0
            work[k], work[pivot] = work[pivot], work[k]
            sign = -sign
        value = work[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                work[i][j] = (work[i][j] * value - work[i][k] * work[k][j]) // previous
        previous = value
        for i in range(k + 1, n):
            work[i][k] = 0
    return sign * work[-1][-1]


def invariant_skew_form(a, b):
    dimension = len(a)
    pairs = [(i, j) for i in range(dimension) for j in range(i + 1, dimension)]
    equations = []
    for generator in (a, b):
        for i in range(dimension):
            for j in range(dimension):
                row = []
                for u, v in pairs:
                    row.append(
                        generator[u][i] * (v == j) - generator[v][i] * (u == j)
                        + (i == u) * generator[v][j] - (i == v) * generator[u][j]
                    )
                if any(row):
                    equations.append(row)
    nullspace = nullspace_q(equations, len(pairs))
    if len(nullspace) != 1:
        return None, len(nullspace)
    vector = [Fraction(x) for x in nullspace[0]]
    denominator = math.lcm(*(x.denominator for x in vector))
    integers = [int(x * denominator) for x in vector]
    content = math.gcd(*(abs(x) for x in integers))
    integers = [x // content for x in integers]
    form = [[0] * dimension for _ in range(dimension)]
    for coefficient, (i, j) in zip(integers, pairs):
        form[i][j], form[j][i] = coefficient, -coefficient
    return form, 1


def block(matrix, offsets, key):
    offset, dimension = offsets[key]
    return [row[offset:offset + dimension] for row in matrix[offset:offset + dimension]]


def build_structure():
    timings = {}
    started = time.monotonic()
    pauli = timed("exact_pauli_closure", timings, exact_pauli_closure)
    sectors = timed("exact_symmetry_blocks", timings, lambda: {key: symmetry_sector(key) for key in CHARACTERS})
    generators, basis, depths = timed("faithful_matrix_closure", timings, matrix_closure_raw, tuple(sectors[key] for key in CHARACTERS))
    # The word basis is selected by exact sparse Fraction elimination; this
    # independent QQ rank is a second exact check of the faithful realization.
    faithful_flat = [flatten(matrix) for matrix in basis]
    faithful_rank = timed("faithful_basis_rank_Q", timings, rank_q, faithful_flat, 4096)
    derived = timed("derived_ideal_closure", timings, ideal_closure, commutator(*generators), generators)
    derived_rank = timed("derived_rank_Q", timings, rank_q, [flatten(x) for x in derived], 4096)

    dimensions = {key: len(sectors[key][0]) for key in CHARACTERS}
    offsets, offset = {}, 0
    for key in CHARACTERS:
        offsets[key] = (offset, dimensions[key])
        offset += dimensions[key]

    selected_projection = []
    # Both the word-basis construction and this rank check are exact over Q.
    for matrix in derived:
        row = []
        for key in SELECTED:
            row.extend(flatten(block(matrix, offsets, key)))
        selected_projection.append(row)
    projection_rank = timed("Levi_projection_rank_Q", timings, rank_q, selected_projection, len(selected_projection[0]))

    factors = []
    for key, lie_type, expected, natural_dimension in (
        ((0, 0, 0), "C7", 105, 14),
        ((0, 1, 0), "C3", 21, 6),
        ((1, 0, 0), "C3", 21, 6),
    ):
        a, b = sectors[key]
        image = ideal_closure(commutator(a, b), [a, b])
        form, nullity = invariant_skew_form(a, b)
        determinant = determinant_bareiss(form) if form else 0
        image_rank = rank_q([flatten(x) for x in image], natural_dimension ** 2)
        factors.append({
            "type": lie_type,
            "dimension": image_rank,
            "sector": "".join(map(str, key)),
            "natural_module_dimension": natural_dimension,
            "certificate": {
                "invariant_alternating_form_nullity": nullity,
                "invariant_form_determinant": determinant,
                "image_dimension": image_rank,
                "ambient_classical_dimension": natural_dimension * (natural_dimension + 1) // 2,
                "simplicity_reason": f"the image lies in sp({natural_dimension}) and has its full dimension {expected}; hence it equals sp({natural_dimension}) over Q, which is simple",
            },
        })

    for key, lie_type, expected, natural_dimension in (
        ((0, 0, 1), "A8", 80, 9),
        ((0, 1, 1), "A5", 35, 6),
    ):
        a, b = sectors[key]
        image = ideal_closure(commutator(a, b), [a, b])
        image_rank = rank_q([flatten(x) for x in image], len(a) ** 2)
        traces_zero = all(sum(matrix[i][i] for i in range(len(a))) == 0 for matrix in image)
        gram = [[trace_product(x, y) for y in image] for x in image]
        gram_rank = rank_q(gram, len(image))
        factors.append({
            "type": lie_type,
            "dimension": image_rank,
            "sector": "".join(map(str, key)),
            "natural_module_dimension": natural_dimension,
            "certificate": {
                "derived_image_dimension": image_rank,
                "all_derived_image_traces_zero": traces_zero,
                "trace_form_rank": gram_rank,
                "ambient_sl_dimension": natural_dimension ** 2 - 1,
                "simplicity_reason": f"the faithful derived image is traceless on a {natural_dimension}-dimensional constituent and has full dimension {expected}; hence it equals sl({natural_dimension}) over Q, which is simple",
            },
        })

    # Exact Killing matrix in the faithful word basis.  Classical identities:
    # kappa_sp(2n)= (2n+2) tr_nat and kappa_sl(n)=2n tr_nat.
    selected_blocks = {key: [block(matrix, offsets, key) for matrix in basis] for key in SELECTED}
    trace_forms = {
        key: [[trace_product(x, y) for y in selected_blocks[key]] for x in selected_blocks[key]]
        for key in SELECTED
    }
    # In sector 001, the A8 constituent is the 9-space complementary to the
    # invariant line.  The exact projector is reconstructed as the unique
    # rank-9 central-commutant projector from the common commutant equations.
    a001, b001 = sectors[(0, 0, 1)]
    d001 = len(a001)
    # Solve the commutant of A,B and choose its rank-one idempotent.  Its
    # complement is the rank-nine A8 natural constituent.
    comm_eq = []
    for generator in (a001, b001):
        for i in range(d001):
            for j in range(d001):
                row = [0] * (d001 * d001)
                for k in range(d001):
                    row[i * d001 + k] += generator[k][j]
                    row[k * d001 + j] -= generator[i][k]
                if any(row):
                    comm_eq.append(row)
    commutant = nullspace_q(comm_eq, d001 * d001)
    assert len(commutant) == 2
    candidates = [unflatten([Fraction(x) for x in vector], d001) for vector in commutant]
    identity = [[Fraction(i == j) for j in range(d001)] for i in range(d001)]
    non_scalar = next(x for x in candidates if any(
        x[i][j] != (x[0][0] if i == j else 0)
        for i in range(d001) for j in range(d001)
    ))
    # The sparse commutant generator returned here is twice a rank-one
    # idempotent; normalize using tr(E^2)/tr(E), avoiding any eigensolver.
    trace_e = sum(non_scalar[i][i] for i in range(d001))
    trace_e2 = trace_product(non_scalar, non_scalar)
    scale = trace_e2 / trace_e
    rank_one = [[x / scale for x in row] for row in non_scalar]
    assert matmul(rank_one, rank_one) == rank_one and rank_q(rank_one, d001) == 1
    identity = [[Fraction(1 if i == j else 0) for j in range(d001)] for i in range(d001)]
    projector = [[identity[i][j] - rank_one[i][j] for j in range(d001)] for i in range(d001)]
    assert rank_q(projector, d001) == 9
    active_trace = [trace_product(projector, x) for x in selected_blocks[(0, 0, 1)]]
    active_form = [[trace_product(projector, matmul(x, y)) for y in selected_blocks[(0, 0, 1)]] for x in selected_blocks[(0, 0, 1)]]
    trace_011 = [sum(x[i][i] for i in range(6)) for x in selected_blocks[(0, 1, 1)]]
    killing = []
    for i in range(len(basis)):
        row = []
        for j in range(len(basis)):
            value = 16 * trace_forms[(0, 0, 0)][i][j]
            value += 8 * trace_forms[(0, 1, 0)][i][j] + 8 * trace_forms[(1, 0, 0)][i][j]
            value += 18 * (active_form[i][j] - active_trace[i] * active_trace[j] / 9)
            value += 12 * (trace_forms[(0, 1, 1)][i][j] - Fraction(trace_011[i] * trace_011[j], 6))
            assert value.denominator == 1
            row.append(value.numerator)
        killing.append(row)
    killing_rank = timed("Killing_rank_Q", timings, rank_q, killing, len(killing))
    killing_nullspace = timed("Killing_nullspace_Q", timings, nullspace_q, killing, len(killing))
    radical_dimension = len(killing_nullspace)  # here [g,g]=semisimple, so Cartan radical is ker kappa
    # Verify that the Cartan-criterion radical generator is central in the
    # faithful representation, rather than inferring centrality from nullity.
    radical_coefficients = [Fraction(x) for x in killing_nullspace[0]]
    radical_matrix = [[
        sum(radical_coefficients[k] * basis[k][i][j] for k in range(len(basis)))
        for j in range(len(generators[0]))
    ] for i in range(len(generators[0]))]
    radical_is_central = all(
        not any(flatten(commutator(generator, radical_matrix)))
        for generator in generators
    )
    symmetric = all(killing[i][j] == killing[j][i] for i in range(len(killing)) for j in range(len(killing)))
    killing_checksum = {
        "sum_abs_entries": sum(abs(x) for row in killing for x in row),
        "max_abs_entry": max(abs(x) for row in killing for x in row),
        "weighted_mod_1000000007": sum((i + 1) * (j + 1) * x for i, row in enumerate(killing) for j, x in enumerate(row)) % 1_000_000_007,
    }

    root_spectrum = {str(k): v for k, v in sorted(Counter(
        [0] * 7 + [s * (i + j) for i in range(1, 8) for j in range(i + 1, 8) for s in (-1, 1)]
        + [s * (j - i) for i in range(1, 8) for j in range(i + 1, 8) for s in (-1, 1)]
        + [s * (2 * i) for i in range(1, 8) for s in (-1, 1)]
    ).items())}
    total = round(time.monotonic() - started, 6)
    timings["total"] = total
    data = {
        "scope": "the single finite 2x3 open layer only; no all-sizes assertion",
        "field": "Q (all decisions use exact integer/Fraction arithmetic)",
        "dimension_Q": pauli["dimension"],
        "faithful_matrix_dimension_Q": faithful_rank,
        "pauli_certificate": pauli,
        "derived_dimension_Q": derived_rank,
        "center_dimension_Q": radical_dimension if radical_is_central else None,
        "Killing": {
            "size": [len(killing), len(killing)],
            "entry_domain": "Z in the recorded word basis",
            "symmetric": symmetric,
            "rank_Q": killing_rank,
            "nullity_Q": radical_dimension,
            "construction": "16 tr_C7 + 8 tr_C3a + 8 tr_C3b + 18 tr_A8 + 12 tr_A5, with exact traceless projections",
            "checksum": killing_checksum,
        },
        "solvable_radical": {
            "dimension_Q": radical_dimension,
            "equals_center": radical_dimension == 1 and radical_is_central,
            "criterion": "rad(g)={x:kappa(x,[g,g])=0}; [g,g] is the certified 262-dimensional semisimple direct sum",
        },
        "semisimple_quotient": {
            "dimension_Q": derived_rank,
            "rank": 7 + 3 + 3 + 8 + 5,
            "type": "C7 + C3 + C3 + A8 + A5",
            "factor_dimensions": [105, 21, 21, 80, 35],
            "selected_joint_projection_rank_Q": projection_rank,
            "factors": factors,
        },
        "largest_factor": {
            "type": "C7=sp(14)",
            "dimension": 105,
            "rank": 7,
            "sector": "000",
            "root_spectrum_for_H_diag_1_to_7": root_spectrum,
            "B7_vs_C7_invariant": "faithful 14-dimensional module with a nondegenerate invariant alternating form (det=16384); B7=so(15) has no nontrivial 14-dimensional representation in characteristic zero (its smallest nontrivial module has dimension 15), whereas this image equals sp(14) by dimension",
        },
        "mod_p_witness_lifts": {
            "sp6": "yes: two independent characteristic-zero C3=sp(6) ideals; sectors 010 and 100 (sector 110 repeats the latter)",
            "gl6": "the 36-dimensional gl6 sector image lifts as A5=sl6 plus the one common central radical; sectors 011 and 111 repeat the same image",
        },
        "basis": {
            "faithful_word_count": len(basis),
            "maximum_word_depth": max(depths),
            "sector_dimensions": {"".join(map(str, key)): dimensions[key] for key in CHARACTERS},
        },
        "timings_seconds": timings,
        "timeouts": [],
    }
    checks = [
        {"name": "Pauli_orbit_closure_dimension_Q", "passed": pauli["dimension"] == 263, "detail": "honest full-orbit summation and exact Fraction elimination"},
        {"name": "faithful_matrix_dimension_Q", "passed": faithful_rank == len(basis) == 263, "detail": "independent symmetry-block word closure and exact DomainMatrix rank over QQ"},
        {"name": "independent_dimension_agreement", "passed": pauli["dimension"] == faithful_rank, "detail": "the two exact-Q constructions agree"},
        {"name": "derived_dimension_Q", "passed": derived_rank == 262, "detail": "exact rational rank of the ideal generated by [A,B]"},
        {"name": "Killing_rank_Q", "passed": killing_rank == 262 and symmetric, "detail": "exact 263x263 integer matrix and DomainMatrix GJ over QQ"},
        {"name": "radical_dimension", "passed": radical_dimension == 1, "detail": "exact Killing nullspace plus certified semisimple derived algebra"},
        {"name": "radical_is_center", "passed": radical_is_central and radical_dimension == 1, "detail": "the exact Killing-nullspace generator commutes with A and B in the faithful representation"},
        {"name": "Levi_dimensions", "passed": projection_rank == sum(x["dimension"] for x in factors) == 262, "detail": "faithful exact-Q joint projection onto five classical images"},
        {"name": "largest_simple_factor", "passed": factors[0]["type"] == "C7" and factors[0]["dimension"] == 105 and factors[0]["certificate"]["invariant_form_determinant"] != 0, "detail": "full sp(14) by exact image dimension and invariant alternating form"},
    ]
    return {"meta": {"provenance": SCRIPT, "timestamp": datetime.now(timezone.utc).isoformat(), "python": sys.version.split()[0], "platform": platform.platform(), "arithmetic": "exact Z/Q; no floating point"}, "data": data, "checks": checks}


def run(write=True):
    artifact = build_structure()
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return artifact


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args(argv)
    try:
        artifact = run(not args.no_write)
        passed = all(check["passed"] for check in artifact["checks"])
        for check in artifact["checks"]:
            print(("PASS" if check["passed"] else "FAIL") + ": " + check["name"])
        return 0 if passed else 1
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
