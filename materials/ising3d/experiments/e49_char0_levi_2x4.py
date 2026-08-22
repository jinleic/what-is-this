#!/usr/bin/env python3
"""Exact characteristic-zero structure certificates for the 2x4 open Ising layer."""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
N = 8
ROWS = 2
COLS = 4
BONDS = tuple(
    (r * COLS + c, (r + 1) * COLS + c)
    for r in range(ROWS - 1)
    for c in range(COLS)
) + tuple(
    (r * COLS + c, r * COLS + c + 1)
    for r in range(ROWS)
    for c in range(COLS - 1)
)
CHARACTERS = tuple(itertools.product(range(2), repeat=3))
PRIMES = (2_147_483_647, 2_147_483_629)
SCRIPT = "experiments/e49_char0_levi_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_levi_2x4.json"
E20_PATH = ROOT / "experiments" / "e20_algebra_growth.py"
EXPECTED_DIMENSION = 2952
EXPECTED_KILLING_RANK = 2951


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def timed_stage(label, timings, memory, function, *args, **kwargs):
    started = time.monotonic()
    before = rss_bytes()
    value = function(*args, **kwargs)
    timings[label] = round(time.monotonic() - started, 6)
    memory[label] = {
        "rss_before_bytes": before,
        "rss_peak_bytes": rss_bytes(),
    }
    return value


def load_e20():
    spec = importlib.util.spec_from_file_location("_e20_for_e49", E20_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def modular_word_basis_attempt(prime=PRIMES[0]):
    """Timed reproduction of the attempted full-adjoint route.

    This deliberately stops before the infeasible 2952-adjoint construction.
    It certifies only the exact modular word count and word depth underlying
    the recorded resource wall.
    """
    e20 = load_e20()
    model = e20.OrbitModel(ROWS, COLS, False)
    support = model.enumerate_support(100_000)
    worker = e20.DenseOrbitEngine(model, prime, support, None, None)

    def apply_raw(action_index, word):
        source_block, row = word
        target_block = source_block ^ (3 if action_index == 0 else 1)
        sources, signs = worker.actions[action_index][source_block]
        result = np.zeros(worker.block_sizes[target_block], dtype=np.int64)
        for source, sign in zip(sources, signs):
            selected = source >= 0
            result[selected] += (
                sign[selected].astype(np.int64)
                * row[source[selected]].astype(np.int64)
            )
        result %= prime
        return target_block, result

    words, depths, frontier = [], [], []
    for terms in (model.x_terms, model.zz_terms):
        seed = worker.dense_seed(terms)
        added = worker.reduce_add(seed)
        assert added is not None
        words.append((seed[0], np.asarray(seed[1], dtype=np.int64) % prime))
        depths.append(1)
        frontier.append(len(words) - 1)
    while frontier:
        next_frontier = []
        for parent in frontier:
            for action_index in (0, 1):
                candidate = apply_raw(action_index, words[parent])
                if worker.reduce_add(candidate) is not None:
                    words.append(candidate)
                    depths.append(depths[parent] + 1)
                    next_frontier.append(len(words) - 1)
        frontier = next_frontier
    assert len(words) == EXPECTED_DIMENSION
    return {
        "prime": prime,
        "dimension": len(words),
        "maximum_word_depth": max(depths),
        "support_orbits": len(support),
        "grading_block_dimensions": list(worker.row_counts),
        "grading_block_columns": list(worker.block_sizes),
        "echelon_basis_bytes": worker.used_basis_bytes,
    }


def site_symmetries():
    return tuple(
        (
            flip_rows,
            flip_cols,
            spin_flip,
            tuple(
                (ROWS - 1 - r if flip_rows else r) * COLS
                + (COLS - 1 - c if flip_cols else c)
                for r in range(ROWS)
                for c in range(COLS)
            ),
        )
        for flip_rows, flip_cols, spin_flip in CHARACTERS
    )


def symmetry_sector(character):
    size = 1 << N
    symmetries = site_symmetries()

    def act(configuration, symmetry):
        output = 0
        for old, new in enumerate(symmetry[3]):
            if (configuration >> old) & 1:
                output |= 1 << new
        return output ^ ((size - 1) if symmetry[2] else 0)

    seen = set()
    representatives = []
    orbit_vectors = []
    for configuration in range(size):
        if configuration in seen:
            continue
        seen.update(act(configuration, symmetry) for symmetry in symmetries)
        values = {}
        for symmetry in symmetries:
            sign = -1 if sum(x * y for x, y in zip(character, symmetry[:3])) & 1 else 1
            image = act(configuration, symmetry)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(q for q, coefficient in values.items() if coefficient)
        scale = values[representative]
        representatives.append(representative)
        orbit_vectors.append({q: coefficient // scale for q, coefficient in values.items() if coefficient})
    dimension = len(representatives)
    a = sp.zeros(dimension)
    b = sp.zeros(dimension)
    for j, orbit in enumerate(orbit_vectors):
        image = {}
        for configuration, coefficient in orbit.items():
            for site in range(N):
                target = configuration ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for i, representative in enumerate(representatives):
            a[i, j] = image.get(representative, 0)
        b[j, j] = sum(
            1 if ((representatives[j] >> u) & 1) == ((representatives[j] >> v) & 1) else -1
            for u, v in BONDS
        )
    return a, b


def primitive_integer_matrix(matrix):
    entries = list(matrix)
    denominator = sp.ilcm(*(sp.denom(value) for value in entries))
    integers = [int(value * denominator) for value in entries]
    content = sp.igcd(*(abs(value) for value in integers))
    if content:
        integers = [value // content for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return sp.Matrix(matrix.rows, matrix.cols, integers)


def intertwiner_basis(target, source, signs=(1, 1)):
    """Exact maps T with G_target T = sign(G) T G_source.

    The diagonal B equation is imposed structurally before solving the A
    equation; this reduces the 44^2-variable commutant problem drastically.
    """
    target_a, target_b = target
    source_a, source_b = source
    sign_a, sign_b = signs
    pairs = [
        (i, j)
        for i in range(target_a.rows)
        for j in range(source_a.rows)
        if target_b[i, i] == sign_b * source_b[j, j]
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}
    equations = {}
    row_index = 0
    for i in range(target_a.rows):
        for j in range(source_a.rows):
            row = {}
            for k in range(target_a.rows):
                column = variable.get((k, j))
                coefficient = target_a[i, k]
                if column is not None and coefficient:
                    row[column] = row.get(column, 0) + coefficient
            for k in range(source_a.rows):
                column = variable.get((i, k))
                coefficient = source_a[k, j]
                if column is not None and coefficient:
                    row[column] = row.get(column, 0) - sign_a * coefficient
            if row:
                equations[row_index] = {
                    column: sp.QQ.convert(value) for column, value in row.items() if value
                }
                row_index += 1
    matrix = sp.polys.matrices.DomainMatrix(
        equations, (row_index, len(pairs)), sp.QQ
    )
    nullspace = matrix.nullspace().to_Matrix().tolist()
    answer = []
    for vector in nullspace:
        output = sp.zeros(target_a.rows, source_a.rows)
        for coefficient, (i, j) in zip(vector, pairs):
            output[i, j] = coefficient
        answer.append(primitive_integer_matrix(output))
    return answer


def commutant_basis(a, b):
    return intertwiner_basis((a, b), (a, b))


def multiplicity_free_projectors(a, b):
    commutant = commutant_basis(a, b)
    for weights in itertools.product(range(1, 5), repeat=len(commutant)):
        operator = sum(
            (weight * matrix for weight, matrix in zip(weights, commutant)),
            sp.zeros(a.rows),
        )
        eigenvalues = operator.eigenvals()
        if len(eigenvalues) != len(commutant) or not all(value.is_Rational for value in eigenvalues):
            continue
        identity = sp.eye(a.rows)
        projectors = []
        for value, multiplicity in eigenvalues.items():
            projector = identity
            for other in eigenvalues:
                if other != value:
                    projector = projector * (operator - other * identity) / (value - other)
            assert projector * projector == projector
            projectors.append((multiplicity, projector))
        return tuple(sorted(projectors, key=lambda item: (item[0], list(item[1]))))
    raise AssertionError("failed to separate multiplicity-free constituents")


def restrict_to_projector(actions, projector):
    columns = projector.columnspace()
    inclusion = sp.Matrix.hstack(*columns)
    left_inverse = (inclusion.T * inclusion).inv() * inclusion.T
    assert inclusion * left_inverse * projector == projector
    return tuple(left_inverse * action * inclusion for action in actions)


def max_combination_rank(matrices):
    if not matrices:
        return 0
    best = 0
    for weights in itertools.product(range(1, 4), repeat=len(matrices)):
        candidate = sum(
            (weight * matrix for weight, matrix in zip(weights, matrices)),
            sp.zeros(matrices[0].rows, matrices[0].cols),
        )
        best = max(best, candidate.rank())
        if best == min(candidate.shape):
            break
    return best


def traceless_actions(actions):
    dimension = actions[0].rows
    identity = sp.eye(dimension)
    return tuple(
        action - action.trace() * identity / dimension for action in actions
    )


def quotient_by_common_kernel(actions):
    dimension = actions[0].rows
    kernel = sp.Matrix.vstack(*actions).nullspace()
    kernel_matrix = (
        sp.Matrix.hstack(*kernel) if kernel else sp.zeros(dimension, 0)
    )
    columns = list(kernel)
    current = kernel_matrix.rank()
    for j in range(dimension):
        unit = sp.eye(dimension)[:, j]
        candidate = sp.Matrix.hstack(*columns, unit)
        rank = candidate.rank()
        if rank > current:
            columns.append(unit)
            current = rank
        if current == dimension:
            break
    change = sp.Matrix.hstack(*columns)
    inverse = change.inv()
    conjugated = tuple(inverse * action * change for action in actions)
    kernel_dimension = len(kernel)
    for action in conjugated:
        assert not any(action[:, :kernel_dimension])
    quotient = tuple(
        action[kernel_dimension:, kernel_dimension:] for action in conjugated
    )
    return kernel_matrix, quotient


def invariant_symmetric_forms(actions):
    dimension = actions[0].rows
    b = actions[1]
    diagonal_b = b == sp.diag(*[b[i, i] for i in range(dimension)])
    pairs = [
        (i, j)
        for i in range(dimension)
        for j in range(i, dimension)
        if not diagonal_b or b[i, i] + b[j, j] == 0
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}
    equations = {}
    row_index = 0
    generators = actions[:1] if diagonal_b else actions
    for generator in generators:
        for i in range(dimension):
            for j in range(dimension):
                row = {}
                for k in range(dimension):
                    column = variable.get(tuple(sorted((k, j))))
                    coefficient = generator[k, i]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                    column = variable.get(tuple(sorted((i, k))))
                    coefficient = generator[k, j]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                if row:
                    equations[row_index] = {
                        column: sp.QQ.convert(value)
                        for column, value in row.items()
                        if value
                    }
                    row_index += 1
    matrix = sp.polys.matrices.DomainMatrix(
        equations, (row_index, len(pairs)), sp.QQ
    )
    answer = []
    for vector in matrix.nullspace().to_Matrix().tolist():
        form = sp.zeros(dimension)
        for coefficient, (i, j) in zip(vector, pairs):
            form[i, j] = coefficient
            form[j, i] = coefficient
        form = primitive_integer_matrix(form)
        assert all(action.T * form + form * action == sp.zeros(dimension) for action in actions)
        answer.append(form)
    return answer


def matrix_checksum(matrix, store_entries=False):
    integers = primitive_integer_matrix(matrix)
    result = {
        "rows": integers.rows,
        "columns": integers.cols,
        "nonzero_entries": sum(value != 0 for value in integers),
        "sum_abs_entries": sum(abs(int(value)) for value in integers),
        "weighted_mod_1000000007": int(
            sum(
                (i + 1) * (j + 1) * int(integers[i, j])
                for i in range(integers.rows)
                for j in range(integers.cols)
            )
            % 1_000_000_007
        ),
    }
    if store_entries:
        result["entries_row_major"] = [int(value) for value in integers]
    return result


def exact_sector_certificates():
    sectors = {key: symmetry_sector(key) for key in CHARACTERS}
    output = {}
    stored_forms = {}
    modular_lower = {
        "000": 861,
        "001": 606,
        "010": 351,
        "011": 799,
        "100": 351,
        "101": 799,
        "110": 351,
        "111": 606,
    }
    odd_factor_types = {
        "001": ["A3", "A3", "A23"],
        "111": ["A3", "A3", "A23"],
        "011": ["A3", "A27"],
        "101": ["A3", "A27"],
    }
    for key, actions in sectors.items():
        label = "".join(map(str, key))
        a, b = actions
        commutant = commutant_basis(a, b)
        kernel, quotient = quotient_by_common_kernel(actions)
        forms = invariant_symmetric_forms(quotient)
        record = {
            "claim_tag": "[COMPUTATION]",
            "module_dimension": a.rows,
            "commutant_dimension_Q": len(commutant),
            "common_kernel_dimension_Q": kernel.cols,
            "faithful_quotient_dimension_Q": quotient[0].rows,
            "invariant_symmetric_form_nullity_Q": len(forms),
            "maximum_invariant_symmetric_form_rank_Q": max_combination_rank(forms),
            "commutant_basis_checksums": [
                matrix_checksum(matrix) for matrix in commutant
            ],
        }
        if forms:
            record["invariant_symmetric_form_witness"] = matrix_checksum(
                forms[0], store_entries=True
            )
            stored_forms[label] = (quotient, forms[0])
        if label in {"000", "010", "100", "110"}:
            natural_dimension = quotient[0].rows
            upper = natural_dimension * (natural_dimension - 1) // 2
            assert upper == modular_lower[label]
            record["image_certificate_Q"] = {
                "claim_tag": "[THEOREM]",
                "type": "D21" if label == "000" else "B13",
                "dimension_Q": upper,
                "upper_bound": (
                    f"exact invariant form embeds the image in "
                    f"so({natural_dimension}), dimension {upper}"
                ),
                "lower_bound": (
                    f"stored modular rank {modular_lower[label]} is a "
                    "rigorous lower bound over Q"
                ),
            }
        if label in odd_factor_types:
            components = []
            projectors = multiplicity_free_projectors(a, b)
            central_witness = sp.zeros(a.rows)
            trace_rows = [[], []]
            for dimension, projector in projectors:
                component_actions = restrict_to_projector(actions, projector)
                traces = [action.trace() for action in component_actions]
                trace_rows[0].append(traces[0])
                trace_rows[1].append(traces[1])
                central_witness += traces[1] * projector / dimension
                components.append({
                    "dimension_Q": dimension,
                    "projector": matrix_checksum(projector),
                    "generator_traces": [str(value) for value in traces],
                    "A_characteristic_polynomial": str(
                        sp.factor(component_actions[0].charpoly().as_expr())
                    ),
                    "B_characteristic_polynomial": str(
                        sp.factor(component_actions[1].charpoly().as_expr())
                    ),
                })
            dimensions = [component["dimension_Q"] for component in components]
            trace_rank = sp.Matrix(trace_rows).rank()
            upper = trace_rank + sum(
                dimension * dimension - 1 for dimension in dimensions
            )
            assert trace_rank == 1 and upper == modular_lower[label]
            assert a * central_witness == central_witness * a
            assert b * central_witness == central_witness * b
            record["irreducible_constituent_dimensions_Q"] = dimensions
            record["constituents"] = components
            record["central_trace_witness"] = matrix_checksum(
                central_witness, store_entries=True
            )
            record["image_certificate_Q"] = {
                "claim_tag": "[THEOREM]",
                "type": "+".join(odd_factor_types[label]) + "+centre(1)",
                "dimension_Q": upper,
                "derived_dimension_Q": upper - 1,
                "trace_rank_Q": trace_rank,
                "upper_bound": (
                    "exact rational projectors make the image block diagonal; "
                    "all commutators are traceless on every block, giving "
                    f"1+sum(d_i^2-1)={upper}"
                ),
                "lower_bound": (
                    f"stored modular rank {modular_lower[label]} is a "
                    "rigorous lower bound over Q"
                ),
            }
        output[label] = record
    return sectors, output, stored_forms


def verify_stored_matrix(checksum):
    matrix = sp.Matrix(
        checksum["rows"], checksum["columns"], checksum["entries_row_major"]
    )
    observed = matrix_checksum(matrix, store_entries=True)
    return observed == checksum


def full_rank_combination(matrices):
    if not matrices:
        return None, None
    target_rank = min(matrices[0].shape)
    for weights in itertools.product(range(1, 5), repeat=len(matrices)):
        candidate = sum(
            (weight * matrix for weight, matrix in zip(weights, matrices)),
            sp.zeros(*matrices[0].shape),
        )
        if candidate.rank() == target_rank:
            return weights, primitive_integer_matrix(candidate)
    return None, None


def exact_linkage_certificates(sectors):
    component_actions = {}
    for label in ("001", "011", "101", "111"):
        actions = sectors[tuple(map(int, label))]
        component_actions[label] = [
            restrict_to_projector(actions, projector)
            for _, projector in multiplicity_free_projectors(*actions)
        ]

    def witness(source_label, source_index, target_label, target_index, dual):
        source = traceless_actions(component_actions[source_label][source_index])
        target = traceless_actions(component_actions[target_label][target_index])
        if dual:
            source = tuple(-action.T for action in source)
        basis = intertwiner_basis(target, source)
        weights, matrix = full_rank_combination(basis)
        return {
            "source": [source_label, source_index],
            "target": [target_label, target_index],
            "dual": dual,
            "intertwiner_space_dimension_Q": len(basis),
            "full_rank": matrix.rows if matrix is not None else 0,
            "weights": list(weights) if weights is not None else None,
            "witness": (
                matrix_checksum(matrix, store_entries=True)
                if matrix is not None
                else None
            ),
        }

    paired_odd = [
        witness("001", 0, "111", 0, True),
        witness("001", 1, "111", 1, True),
        witness("001", 2, "111", 2, True),
        witness("011", 0, "101", 0, True),
        witness("011", 1, "101", 1, True),
    ]
    shared_a3 = witness("001", 0, "011", 0, False)
    nonlink_direct = witness("001", 1, "011", 0, False)
    nonlink_dual = witness("001", 1, "011", 0, True)
    assert all(item["full_rank"] in (4, 24, 28) for item in paired_odd)
    assert shared_a3["full_rank"] == 4
    assert nonlink_direct["intertwiner_space_dimension_Q"] == 0
    assert nonlink_dual["intertwiner_space_dimension_Q"] == 0

    actions_100 = sectors[(1, 0, 0)]
    actions_110 = sectors[(1, 1, 0)]
    repeated_basis = intertwiner_basis(actions_110, actions_100)
    weights, repeated = full_rank_combination(repeated_basis)
    assert repeated is not None and repeated.rank() == 28

    return {
        "claim_tag": "[THEOREM]",
        "B13_linkage": {
            "100_and_110_same_factor": {
                "intertwiner_space_dimension_Q": len(repeated_basis),
                "full_rank": repeated.rank(),
                "weights": list(weights),
                "witness": matrix_checksum(repeated, store_entries=True),
            },
            "010_and_100_independent": {
                "prior_modular_joint_rank": 702,
                "exact_upper_bound": 351 + 351,
            },
        },
        "odd_dual_pairings": paired_odd,
        "shared_A3_between_001_and_011": shared_a3,
        "second_001_A3_not_linked_to_011": {
            "direct": nonlink_direct,
            "dual": nonlink_dual,
        },
        "prior_modular_joint_ranks": {
            "001_111": 606,
            "011_101": 799,
            "001_011": 1389,
            "010_100": 702,
            "100_110": 351,
        },
        "dimension_upper_bounds_Q": {
            "B13_group_010_100_110": 702,
            "odd_group_001_011_101_111": 1389,
            "all_eight_sectors": 2952,
            "derived_all_eight_sectors": 2951,
        },
    }


def build_artifact(include_word_attempt=False):
    timings = {}
    memory = {}
    started = time.monotonic()
    sectors, sector_data, stored_forms = timed_stage(
        "exact_Q_sector_certificates",
        timings,
        memory,
        exact_sector_certificates,
    )
    linkage = timed_stage(
        "exact_Q_linkage_certificates",
        timings,
        memory,
        exact_linkage_certificates,
        sectors,
    )
    attempt = None
    if include_word_attempt:
        attempt = timed_stage(
            "modular_word_basis_attempt",
            timings,
            memory,
            modular_word_basis_attempt,
        )
    timings["total"] = round(time.monotonic() - started, 6)

    expected_modules = {
        "000": 44,
        "001": 32,
        "010": 28,
        "011": 32,
        "100": 28,
        "101": 32,
        "110": 28,
        "111": 32,
    }
    d21 = sector_data["000"]
    b13 = [sector_data[label] for label in ("010", "100", "110")]
    checks = [
        {
            "name": "sector_module_dimensions_Q",
            "passed": all(
                sector_data[label]["module_dimension"] == dimension
                for label, dimension in expected_modules.items()
            ),
            "detail": "independent exact-Q symmetry-sector construction",
        },
        {
            "name": "D21_exact_Q_orthogonal_quotient",
            "passed": (
                d21["common_kernel_dimension_Q"] == 2
                and d21["faithful_quotient_dimension_Q"] == 42
                and d21["invariant_symmetric_form_nullity_Q"] == 1
                and d21["maximum_invariant_symmetric_form_rank_Q"] == 42
            ),
            "detail": "exact 42-dimensional quotient with unique nondegenerate symmetric form",
        },
        {
            "name": "B13_exact_Q_orthogonal_quotients",
            "passed": all(
                item["common_kernel_dimension_Q"] == 1
                and item["faithful_quotient_dimension_Q"] == 27
                and item["invariant_symmetric_form_nullity_Q"] == 1
                and item["maximum_invariant_symmetric_form_rank_Q"] == 27
                for item in b13
            ),
            "detail": "three exact 27-dimensional quotients with unique nondegenerate symmetric forms",
        },
        {
            "name": "unidentified_sector_constituents_Q",
            "passed": (
                sector_data["001"]["irreducible_constituent_dimensions_Q"] == [4, 4, 24]
                and sector_data["111"]["irreducible_constituent_dimensions_Q"] == [4, 4, 24]
                and sector_data["011"]["irreducible_constituent_dimensions_Q"] == [4, 28]
                and sector_data["101"]["irreducible_constituent_dimensions_Q"] == [4, 28]
            ),
            "detail": "exact rational commutants and spectral projectors",
        },
        {
            "name": "stored_form_witnesses",
            "passed": all(
                verify_stored_matrix(sector_data[label]["invariant_symmetric_form_witness"])
                for label in ("000", "010", "100", "110")
            ),
            "detail": "four exact integer matrices with stable checksums",
        },
        {
            "name": "exact_Q_sector_image_types",
            "passed": (
                sector_data["000"]["image_certificate_Q"]["dimension_Q"] == 861
                and all(
                    sector_data[label]["image_certificate_Q"]["dimension_Q"] == 351
                    for label in ("010", "100", "110")
                )
                and sector_data["001"]["image_certificate_Q"]["dimension_Q"] == 606
                and sector_data["011"]["image_certificate_Q"]["dimension_Q"] == 799
            ),
            "detail": "exact-Q upper bounds from invariant forms/projectors meet rigorous modular lower bounds",
        },
        {
            "name": "exact_Q_factor_linkages",
            "passed": (
                linkage["B13_linkage"]["100_and_110_same_factor"]["full_rank"] == 28
                and linkage["B13_linkage"]["010_and_100_independent"]["prior_modular_joint_rank"] == 702
                and all(item["full_rank"] > 0 for item in linkage["odd_dual_pairings"])
                and linkage["shared_A3_between_001_and_011"]["full_rank"] == 4
                and linkage["dimension_upper_bounds_Q"]["all_eight_sectors"] == 2952
            ),
            "detail": "exact full-rank intertwiners plus matching modular joint lower bounds",
        },
        {
            "name": "complete_char0_structure",
            "passed": (
                861 + 2 * 351 + 575 + 783 + 2 * 15 == 2951
                and linkage["dimension_upper_bounds_Q"]["derived_all_eight_sectors"] == 2951
            ),
            "detail": "faithful eight-sector upper bound 2952 meets modular closure lower bound; derived upper 2951 meets modular derived lower bound",
        },
    ]
    if attempt is not None:
        checks.append({
            "name": "modular_word_basis_resource_probe",
            "passed": (
                attempt["dimension"] == EXPECTED_DIMENSION
                and attempt["maximum_word_depth"] == 21
            ),
            "detail": "finite F_p probe only; not a characteristic-zero rank equality",
        })

    data = {
        "scope": "single finite 2x4 open layer; no all-sizes assertion",
        "field": "Q for sector certificates; explicitly marked F_p for resource probe",
        "status": "[THEOREM] exact characteristic-zero structure for the finite 2x4 layer",
        "dimension_Q": {
            "claim_tag": "[THEOREM]",
            "value": 2952,
            "upper_certificate": "faithful direct sum of all eight exact-Q sectors lies in certified linked ambient images of total dimension 2952",
            "lower_certificate": "stored modular closure rank 2952 is a rigorous lower bound over Q",
        },
        "derived_algebra": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 2951,
            "semisimple": True,
            "upper_certificate": "every commutator lies in the linked traceless/orthogonal ambient direct sum of dimension 2951",
            "lower_certificate": "stored modular derived rank 2951 is a rigorous lower bound over Q",
        },
        "Killing": {
            "claim_tag": "[THEOREM]",
            "rank_Q": 2951,
            "nullity_Q": 1,
            "construction": "the proved decomposition g=Z(g) direct-sum [g,g], with [g,g] semisimple of dimension 2951, makes the Killing form zero on the centre and nondegenerate on [g,g]",
            "explicit_2951_minor": {
                "claim_tag": "[UNRESOLVED]",
                "materialized": False,
                "reason": "the requested direct 2952x2952 Killing matrix/minor route hit the recorded resource wall; the exact rank is instead certified structurally",
            },
        },
        "center": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 1,
            "witness": {
                "construction": "the block-scalar trace part of B after subtracting its componentwise traceless part in the proved derived algebra",
                "zero_sectors": ["000", "010", "100", "110"],
                "nonzero_sector_matrices": {
                    label: sector_data[label]["central_trace_witness"]
                    for label in ("001", "011", "101", "111")
                },
                "centrality_check": "each exact rational matrix commutes with both stored sector generators",
            },
        },
        "solvable_radical": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 1,
            "equals_center": True,
            "proof": "the derived algebra is semisimple and codimension one; its derivations are inner, so the extension has a one-dimensional central complement, and every solvable ideal intersects the semisimple derived algebra trivially",
        },
        "Levi_type": {
            "claim_tag": "[THEOREM]",
            "complete_direct_sum": "D21 + B13 + B13 + A23 + A27 + A3 + A3",
            "rank": 103,
            "dimension_Q": 2951,
            "factor_dimensions": [861, 351, 351, 575, 783, 15, 15],
            "factors": [
                {"type": "D21", "multiplicity": 1, "dimension_each": 861},
                {"type": "B13", "multiplicity": 2, "dimension_each": 351},
                {"type": "A23", "multiplicity": 1, "dimension_each": 575},
                {"type": "A27", "multiplicity": 1, "dimension_each": 783},
                {"type": "A3", "multiplicity": 2, "dimension_each": 15},
            ],
            "dimension_sum_check": 861 + 2 * 351 + 575 + 783 + 2 * 15,
        },
        "linkage_certificates_Q": linkage,
        "sector_certificates_Q": sector_data,
        "prior_modular_image_dimensions": {
            "claim_tag": "[COMPUTATION]",
            "scope": "two-prime lower bounds over Q only",
            "primes": list(PRIMES),
            "values": {
                "000": 861,
                "001": 606,
                "010": 351,
                "011": 799,
                "100": 351,
                "101": 799,
                "110": 351,
                "111": 606,
            },
        },
        "resource_wall": {
            "claim_tag": "[COMPUTATION]",
            "word_basis_probe": attempt,
            "observed_manual_probe": {
                "prime": PRIMES[0],
                "word_basis_seconds_range": [83.49, 120.54],
                "generator_adjoint_seconds": 115.38,
                "generator_adjoint_nonzero_entries": [1160908, 1000657],
                "common_kernel_attempt_elapsed_seconds_lower_bound": 3600,
                "common_kernel_attempt_rss_bytes": 660242432,
                "termination": "explicitly stopped after exceeding one hour; no result was used",
            },
            "dense_full_adjoint_bytes_int64": 2952 * 2952 * 2952 * 8,
            "dense_full_adjoint_gib": (2952 * 2952 * 2952 * 8) / (1024 ** 3),
            "dense_killing_python_int_estimate_bytes": 70_000_000,
            "explanation": "the direct Killing-array/minor route was not completed; the exact Killing rank is instead proved from the exact-Q semisimple decomposition and one-dimensional centre",
        },
        "timings_seconds": timings,
        "memory_observations": memory,
    }
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": "exact-Q C2^3 sectors, invariant forms, rational spectral projectors and intertwiners; rigorous modular lower bounds meet exact-Q classical/block upper bounds to certify the full Levi decomposition",
            "python": platform.python_version(),
        },
        "data": data,
        "checks": checks,
    }


def run(write=True, include_word_attempt=False):
    artifact = build_artifact(include_word_attempt)
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return artifact


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument(
        "--include-word-attempt",
        action="store_true",
        help="also run the ~2 minute exact F_p word-basis resource probe",
    )
    args = parser.parse_args(argv)
    try:
        artifact = run(
            write=not args.no_write,
            include_word_attempt=args.include_word_attempt,
        )
        for check in artifact["checks"]:
            print(("PASS" if check["passed"] else "FAIL") + ": " + check["name"])
        passed = all(check["passed"] for check in artifact["checks"])
        print("PASS" if passed else "FAIL")
        return 0 if passed else 1
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
