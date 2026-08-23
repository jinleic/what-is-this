#!/usr/bin/env python3
"""Clean-room verifier for the Kac--Ward section-orbit and unit-census result.

This test imports none of e169/e170/e171.  It independently rebuilds the raw
42-equation system, derives the group action with a custom Laurent-monomial
implementation, enumerates all 4^13 unit points using generator 3 (the producer
uses generator 2), recomputes every section status, and replays every exact
holdout or lift certificate.  A universal classification row is accepted only
when this independent census covers it.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    ALLOWED_DIRECTION_PAIRS,
    DIRECTION_GAUGE_TREE,
    DIRECTIONS_3D,
    _transition_indices,
    directed_edges,
    formal_log_coefficients,
    full_weight_finite_system,
)
from ising.lattices import cubic

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/orbit_quotient.json"
E169 = ROOT / "experiments/e169_kw_section_group.py"
E170 = ROOT / "experiments/e170_kw_unit_census.py"
E171 = ROOT / "experiments/e171_kw_orbit_quotient.py"
E139 = ROOT / "experiments/e139_kw_components.py"
UNIT_VALUES = (1, 2, 3, 4)
GENERATOR3_VALUES = (1, 3, 4, 2)
TOTAL_POINTS = 4**13
TOTAL_SECTIONS = 4**7
CAP_SECONDS = 180.0
CHUNK_SIZE = 1 << 19
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(set(itertools.permutations((2, 2, 1))) | set(itertools.permutations((3, 2, 1))))
)
FULL_SHAPES = ((3, 3, 2), (2, 2, 3), (2, 2, 2), (3, 2, 2), (2, 3, 2))
CONSTRUCTION_SHAPES = PRIMARY_SHAPES + FULL_SHAPES
SELECTED_INDICES = (28, 29, 32, 35, 38, 41)
DIAGONAL_NAMES = ("u_px_px", "u_py_py", "u_pz_pz")
NONDIAGONAL_NAMES = (
    "u_py_px", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx", "u_pz_px", "u_pz_mx",
    "u_pz_py", "u_pz_my", "u_mz_px", "u_mz_mx", "u_mz_py", "u_mz_my",
)
FREE_NAMES = (
    "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx",
    "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx",
    "u_mz_py", "u_mz_my",
)
DIRECTION_NAMES = ("px", "mx", "py", "my", "pz", "mz")


def require_time(deadline: float) -> None:
    if time.process_time() >= deadline:
        raise AssertionError("independent verifier exceeded its declared process-time cap")


def pair_name(pair: tuple[int, int]) -> str:
    return f"u_{DIRECTION_NAMES[pair[0]]}_{DIRECTION_NAMES[pair[1]]}"


def branch_equations(system):
    fixed = set(DIRECTION_GAUGE_TREE)
    substitutions = {symbol: 1 for pair, symbol in system.all_symbols if pair in fixed}
    variables = [symbol for pair, symbol in system.all_symbols if pair not in fixed]
    equations = [sp.expand(expression.subs(substitutions)) for expression in system.equations]
    return variables, equations


def thin_substitutions(symbol):
    return {
        symbol["u_mx_my"]: -1 / symbol["u_my_px"],
        symbol["u_mx_py"]: -1 / (symbol["u_my_mx"] * symbol["u_py_px"]),
        symbol["u_mx_mx"]: 1 / symbol["u_px_px"],
        symbol["u_mx_mz"]: -1 / (symbol["u_mz_px"] * symbol["u_pz_mx"]),
        symbol["u_mx_pz"]: -1 / (symbol["u_mz_mx"] * symbol["u_pz_px"]),
        symbol["u_mz_mz"]: 1 / symbol["u_pz_pz"],
        symbol["u_my_mz"]: -1 / (symbol["u_mz_py"] * symbol["u_py_pz"] * symbol["u_pz_my"]),
        symbol["u_my_pz"]: -1 / (symbol["u_mz_my"] * symbol["u_py_mz"] * symbol["u_pz_py"]),
        symbol["u_my_my"]: 1 / symbol["u_py_py"],
    }


def build_raw_system():
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False)
    variables, equations = branch_equations(system)
    symbol = {str(value): value for value in variables}
    free = [symbol[name] for name in FREE_NAMES]
    non_diagonal = [symbol[name] for name in NONDIAGONAL_NAMES]
    thin = thin_substitutions(symbol)
    diagonal = {symbol[name]: 1 for name in DIAGONAL_NAMES}
    primitives = []
    for index in SELECTED_INDICES:
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        polynomial = sp.Poly(sp.expand(numerator).subs(diagonal), *non_diagonal, domain=sp.ZZ)
        _content, primitive = polynomial.primitive()
        primitives.append(primitive)
    assert [len(polynomial.terms()) for polynomial in primitives] == [33, 277, 155, 7, 155, 155]
    return system, variables, equations, symbol, free, non_diagonal, thin, primitives


# A Laurent monomial is (integer sign, 13 integer exponents).  All group-action
# expressions have coefficient +/-1, so this tiny implementation is independent
# of the producer's SymPy expression manipulations.
Monomial = tuple[int, tuple[int, ...]]


def mono_mul(left: Monomial, right: Monomial) -> Monomial:
    return left[0] * right[0], tuple(a + b for a, b in zip(left[1], right[1], strict=True))


def mono_inv(value: Monomial) -> Monomial:
    return value[0], tuple(-entry for entry in value[1])


def mono_div(left: Monomial, right: Monomial) -> Monomial:
    return mono_mul(left, mono_inv(right))


def mono_variable(index: int) -> Monomial:
    exponents = [0] * 13
    exponents[index] = 1
    return 1, tuple(exponents)


def raw_chart_monomials() -> dict[tuple[int, int], Monomial]:
    one: Monomial = (1, (0,) * 13)
    value = {name: mono_variable(index) for index, name in enumerate(NONDIAGONAL_NAMES)}
    raw: dict[tuple[int, int], Monomial] = {}
    for pair in ALLOWED_DIRECTION_PAIRS:
        name = pair_name(pair)
        if pair in set(DIRECTION_GAUGE_TREE) or name in DIAGONAL_NAMES:
            raw[pair] = one
        elif name in value:
            raw[pair] = value[name]
    dependencies = {
        "u_mx_my": ("u_my_px",),
        "u_mx_py": ("u_my_mx", "u_py_px"),
        "u_mx_mx": (),
        "u_mx_mz": ("u_mz_px", "u_pz_mx"),
        "u_mx_pz": ("u_mz_mx", "u_pz_px"),
        "u_mz_mz": (),
        "u_my_mz": ("u_mz_py", "u_py_pz", "u_pz_my"),
        "u_my_pz": ("u_mz_my", "u_py_mz", "u_pz_py"),
        "u_my_my": (),
    }
    for pair in ALLOWED_DIRECTION_PAIRS:
        name = pair_name(pair)
        if name not in dependencies:
            continue
        if not dependencies[name]:
            raw[pair] = one
            continue
        denominator = one
        for dependency in dependencies[name]:
            denominator = mono_mul(denominator, value[dependency])
        raw[pair] = mono_mul((-1, (0,) * 13), mono_inv(denominator))
    assert set(raw) == set(ALLOWED_DIRECTION_PAIRS)
    return raw


def direction_permutation(axis_permutation, signs):
    index = {tuple(vector): position for position, vector in enumerate(DIRECTIONS_3D)}
    return tuple(
        index[tuple(signs[axis] * vector[axis_permutation[axis]] for axis in range(3))]
        for vector in DIRECTIONS_3D
    )


def normalized_point_action(raw, permutation):
    inverse = [0] * 6
    for old, new in enumerate(permutation):
        inverse[new] = old
    transformed = {
        pair: raw[(inverse[pair[0]], inverse[pair[1]])] for pair in ALLOWED_DIRECTION_PAIRS
    }
    one: Monomial = (1, (0,) * 13)
    potential = {0: one}
    pending = list(DIRECTION_GAUGE_TREE)
    while pending:
        next_pending = []
        progress = False
        for source, target in pending:
            if source in potential and target not in potential:
                potential[target] = mono_div(potential[source], transformed[(source, target)])
                progress = True
            elif target in potential and source not in potential:
                potential[source] = mono_mul(potential[target], transformed[(source, target)])
                progress = True
            else:
                next_pending.append((source, target))
        assert progress or not next_pending
        pending = next_pending
    action = {
        pair: mono_mul(
            transformed[pair], mono_div(potential[pair[1]], potential[pair[0]])
        )
        for pair in ALLOWED_DIRECTION_PAIRS
    }
    assert all(action[pair] == one for pair in DIRECTION_GAUGE_TREE)
    return action


def verify_group(equations, thin, free):
    shapes = set(CONSTRUCTION_SHAPES)
    raw = raw_chart_monomials()
    full_order = 0
    catalog_order = 0
    section_order = 0
    action_records = []
    for axis_permutation in itertools.permutations(range(3)):
        mapped = {tuple(shape[axis_permutation[axis]] for axis in range(3)) for shape in shapes}
        preserves = mapped == shapes
        for signs in itertools.product((1, -1), repeat=3):
            full_order += 1
            if not preserves:
                continue
            catalog_order += 1
            permutation = direction_permutation(axis_permutation, signs)
            action = normalized_point_action(raw, permutation)
            coordinate_action = {
                pair_name(pair): monomial for pair, monomial in action.items()
                if pair_name(pair) in NONDIAGONAL_NAMES
            }
            induces = all(
                not any(coordinate_action[name][1][coordinate] for coordinate in range(6))
                for name in NONDIAGONAL_NAMES[6:]
            )
            section_order += int(induces)
            action_records.append({
                "axis_permutation": tuple(axis_permutation),
                "axis_signs": tuple(signs),
                "coordinate_action": coordinate_action,
                "induces_section_map": induces,
            })
    assert (full_order, catalog_order, section_order) == (48, 16, 1)

    num16 = []
    for index in SELECTED_INDICES:
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        num16.append(sp.Poly(sp.expand(numerator), *free, domain=sp.ZZ))
    rows = []
    for polynomial in num16:
        terms = polynomial.terms()
        origin = terms[0][0]
        rows.extend(
            [left - right for left, right in zip(monomial, origin, strict=True)]
            for monomial, _ in terms[1:]
        )
    matrix = sp.Matrix(rows)
    assert matrix.rank() == 13
    kernel = matrix.nullspace()
    assert len(kernel) == 3
    diagonal_positions = [0, 2, 11]
    diagonal_minor = sp.Matrix(
        [[kernel[row][column] for column in diagonal_positions] for row in range(3)]
    )
    assert abs(int(diagonal_minor.det())) == 1
    # Unimodularity makes the diagonal-normalized torus stabilizer trivial; the
    # five-edge tree on six vertices likewise leaves only the ineffective common gauge.
    assert len(DIRECTION_GAUGE_TREE) == 5
    orbit_representatives = set(itertools.product(UNIT_VALUES, repeat=7))
    assert len(orbit_representatives) == TOTAL_SECTIONS
    return full_order, catalog_order, len(orbit_representatives), action_records


def evaluate_custom_action(point, coordinate_action):
    image = []
    for name in NONDIAGONAL_NAMES:
        coefficient, exponents = coordinate_action[name]
        value = coefficient % 5
        for coordinate, exponent in enumerate(exponents):
            if exponent:
                value = value * pow(point[coordinate], exponent % 4, 5) % 5
        image.append(value)
    return tuple(image)


def verify_variety_fiber_counterexamples(payload, points, action_records):
    stored = payload["data"]["symmetry"]["section_action"][
        "construction_variety_counterexamples_mod5"
    ]
    assert len(stored) == 15
    stored_by_element = {
        (tuple(row["axis_permutation"]), tuple(row["axis_signs"])): row for row in stored
    }
    nonidentity_count = 0
    for action_record in action_records:
        action = action_record["coordinate_action"]
        images = {point: evaluate_custom_action(point, action) for point in points}
        assert set(images.values()) == points
        if action_record["induces_section_map"]:
            continue
        nonidentity_count += 1
        key = (action_record["axis_permutation"], action_record["axis_signs"])
        witness = stored_by_element[key]
        first = tuple(witness["first_point13_mod5"])
        second = tuple(witness["second_point13_mod5"])
        assert first in points and second in points and first[6:] == second[6:]
        assert images[first] == tuple(witness["first_image13_mod5"])
        assert images[second] == tuple(witness["second_image13_mod5"])
        assert images[first][6:] != images[second][6:]
        assert witness["input_section"] == list(first[6:])
        assert witness["distinct_target_sections"] == [
            list(images[first][6:]), list(images[second][6:])
        ]
    assert nonidentity_count == 15


def enumerate_unit_points(primitives, deadline: float) -> set[tuple[int, ...]]:
    exponent_tables = []
    coefficient_tables = []
    for polynomial in primitives:
        terms = polynomial.terms()
        exponent_tables.append(
            np.array([[int(entry) % 4 for entry in monomial] for monomial, _ in terms], dtype=np.uint8)
        )
        coefficient_tables.append(np.array([int(value) % 5 for _, value in terms], dtype=np.int16))
    # Deliberately differ from the producer: generator 3, 2^19 chunks, and a
    # different tie order among the three 155-term equations.
    order = (3, 0, 4, 2, 5, 1)
    powers = np.array(GENERATOR3_VALUES, dtype=np.int16)
    result: set[tuple[int, ...]] = set()
    for start in range(0, TOTAL_POINTS, CHUNK_SIZE):
        require_time(deadline)
        count = min(CHUNK_SIZE, TOTAL_POINTS - start)
        indices = np.arange(start, start + count, dtype=np.uint32)
        logs = np.empty((count, 13), dtype=np.uint8)
        for coordinate in range(13):
            logs[:, coordinate] = (indices >> (2 * coordinate)) & 3
        live = np.arange(count, dtype=np.int32)
        for polynomial_index in order:
            phases = (logs[live] @ exponent_tables[polynomial_index].T) & 3
            residues = (
                powers[phases] * coefficient_tables[polynomial_index]
            ).sum(axis=1, dtype=np.int32) % 5
            live = live[residues == 0]
            if not len(live):
                break
        for index in indices[live]:
            result.add(
                tuple(GENERATOR3_VALUES[(int(index) >> (2 * coordinate)) & 3] for coordinate in range(13))
            )
    assert len(result) == 2_960
    return result


def poly_terms(polynomial):
    return [(tuple(map(int, monomial)), int(coefficient)) for monomial, coefficient in polynomial.terms()]


def derivative_terms(terms, coordinate):
    result = []
    for exponents, coefficient in terms:
        if exponents[coordinate]:
            derivative = list(exponents)
            derivative[coordinate] -= 1
            result.append((tuple(derivative), coefficient * exponents[coordinate]))
    return result


def eval_terms(terms, values, modulus):
    total = 0
    for exponents, coefficient in terms:
        term = coefficient % modulus
        for value, exponent in zip(values, exponents, strict=True):
            if exponent:
                term = term * pow(int(value) % modulus, exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def numeric_raw_weights(point, modulus):
    value = dict(zip(NONDIAGONAL_NAMES, point, strict=True))
    result = {}
    for pair in ALLOWED_DIRECTION_PAIRS:
        name = pair_name(pair)
        if pair in set(DIRECTION_GAUGE_TREE) or name in DIAGONAL_NAMES or name in {"u_mx_mx", "u_my_my", "u_mz_mz"}:
            result[pair] = 1
        elif name in value:
            result[pair] = int(value[name]) % modulus
    inverse = {
        "u_mx_my": ("u_my_px",),
        "u_mx_py": ("u_my_mx", "u_py_px"),
        "u_mx_mz": ("u_mz_px", "u_pz_mx"),
        "u_mx_pz": ("u_mz_mx", "u_pz_px"),
        "u_my_mz": ("u_mz_py", "u_py_pz", "u_pz_my"),
        "u_my_pz": ("u_mz_my", "u_py_mz", "u_pz_py"),
    }
    for pair in ALLOWED_DIRECTION_PAIRS:
        name = pair_name(pair)
        if name in inverse:
            denominator = 1
            for dependency in inverse[name]:
                denominator = denominator * int(value[dependency]) % modulus
            result[pair] = -pow(denominator, -1, modulus) % modulus
    assert set(result) == set(ALLOWED_DIRECTION_PAIRS)
    return result


def transfer_structure():
    lattice = cubic(3, 3, 3, periodic=False)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    direction_index = {tuple(vector): index for index, vector in enumerate(DIRECTIONS_3D)}
    directions = [direction_index[tuple(edge.direction)] for edge in edges]
    return edges, transitions, directions


def holdout_residue(point, order, modulus, logs, structure):
    _edges, transitions, directions = structure
    weights = numeric_raw_weights(point, modulus)
    size = len(transitions)
    matrix = np.zeros((size, size), dtype=np.int64)
    for edge, next_edges in enumerate(transitions):
        for next_edge in next_edges:
            matrix[edge, next_edge] = weights[(directions[edge], directions[next_edge])]
    result = np.eye(size, dtype=np.int64)
    base = matrix
    exponent = order
    while exponent:
        if exponent & 1:
            result = (result @ base) % modulus
        exponent //= 2
        if exponent:
            base = (base @ base) % modulus
    return (int(np.trace(result) % modulus) + 2 * order * logs[order]) % modulus


def rref_data(matrix, rhs):
    rows = len(matrix)
    columns = len(matrix[0])
    augmented = [
        [int(matrix[row][column]) % 5 for column in range(columns)] + [int(rhs[row]) % 5]
        for row in range(rows)
    ]
    pivots = []
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows) if augmented[row][column]), None)
        if pivot is None:
            continue
        augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
        inverse = pow(augmented[pivot_row][column], -1, 5)
        augmented[pivot_row] = [value * inverse % 5 for value in augmented[pivot_row]]
        for row in range(rows):
            if row != pivot_row and augmented[row][column]:
                factor = augmented[row][column]
                augmented[row] = [
                    (left - factor * right) % 5
                    for left, right in zip(augmented[row], augmented[pivot_row], strict=True)
                ]
        pivots.append(column)
        pivot_row += 1
    inconsistent = any(
        all(augmented[row][column] == 0 for column in range(columns)) and augmented[row][-1]
        for row in range(pivot_row, rows)
    )
    return augmented, pivots, inconsistent


def affine_solutions(matrix, rhs):
    reduced, pivots, inconsistent = rref_data(matrix, rhs)
    if inconsistent:
        return []
    columns = len(matrix[0])
    free = [column for column in range(columns) if column not in pivots]
    particular = [0] * columns
    for row, column in enumerate(pivots):
        particular[column] = reduced[row][-1]
    basis = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, column in enumerate(pivots):
            vector[column] = -reduced[row][free_column] % 5
        basis.append(vector)
    result = []
    for coefficients in itertools.product(range(5), repeat=len(basis)):
        vector = particular.copy()
        for coefficient, basis_vector in zip(coefficients, basis, strict=True):
            vector = [
                (left + coefficient * right) % 5
                for left, right in zip(vector, basis_vector, strict=True)
            ]
        result.append(vector)
    return result


def determinant(matrix):
    reduced, pivots, inconsistent = rref_data(matrix, [0] * len(matrix))
    del reduced
    if inconsistent or len(pivots) < len(matrix):
        return 0
    # A separate small elimination retains determinant factors.
    work = [[int(value) % 5 for value in row] for row in matrix]
    result = 1
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result = result * pivot_value % 5
        inverse = pow(pivot_value, -1, 5)
        for row in range(column + 1, len(work)):
            factor = work[row][column] * inverse % 5
            work[row] = [
                (left - factor * right) % 5
                for left, right in zip(work[row], work[column], strict=True)
            ]
    return result % 5


def reconstruct25(system, variables, point, modulus):
    weights = numeric_raw_weights(point, modulus)
    by_name = {pair_name(pair): value for pair, value in weights.items()}
    return [by_name[str(variable)] for variable in variables]


def full42_zero(equation_terms, values25, modulus):
    return all(eval_terms(terms, values25, modulus) == 0 for terms in equation_terms)


def nonsingular_lift(primitive_terms, point, jacobian):
    solve = list(point[:6])
    held = tuple(point[6:])
    modulus = 5
    for _ in range(3):
        next_modulus = 5 * modulus
        residues = [eval_terms(terms, tuple(solve) + held, next_modulus) for terms in primitive_terms]
        rhs = [-(residue // modulus) % 5 for residue in residues]
        corrections = affine_solutions(jacobian, rhs)
        assert len(corrections) == 1
        solve = [
            (value + modulus * digit) % next_modulus
            for value, digit in zip(solve, corrections[0], strict=True)
        ]
        modulus = next_modulus
    lifted = tuple(solve) + held
    assert all(eval_terms(terms, lifted, 625) == 0 for terms in primitive_terms)
    return lifted


def singular_tree(
    point,
    jacobian,
    primitive_terms,
    system,
    variables,
    equation_terms,
    logs,
    structure,
    deadline,
):
    current = [tuple(point[:6])]
    held = tuple(point[6:])
    modulus = 5
    levels = []
    for _ in range(3):
        require_time(deadline)
        next_modulus = 5 * modulus
        next_points = []
        branches = []
        for solve in current:
            residues = [eval_terms(terms, solve + held, next_modulus) for terms in primitive_terms]
            rhs = [-(residue // modulus) % 5 for residue in residues]
            corrections = affine_solutions(jacobian, rhs)
            candidates = []
            for correction in corrections:
                lifted_solve = tuple(
                    (value + modulus * digit) % next_modulus
                    for value, digit in zip(solve, correction, strict=True)
                )
                lifted = lifted_solve + held
                assert all(eval_terms(terms, lifted, next_modulus) == 0 for terms in primitive_terms)
                values25 = reconstruct25(system, variables, lifted, next_modulus)
                assert full42_zero(equation_terms, values25, next_modulus)
                residue8 = holdout_residue(lifted, 8, next_modulus, logs, structure)
                candidates.append({
                    "correction_mod5": correction,
                    "solve6": list(lifted_solve),
                    "construction_residues_all_zero": True,
                    "holdout_order8_residue": residue8,
                })
                if residue8 == 0:
                    next_points.append(lifted_solve)
            branches.append({
                "input_solve6": list(solve),
                "linear_rhs_mod5": rhs,
                "correction_count": len(corrections),
                "candidates": candidates,
            })
        levels.append({
            "modulus": next_modulus,
            "input_count": len(current),
            "construction_lift_count": sum(branch["correction_count"] for branch in branches),
            "holdout_zero_lift_count": len(next_points),
            "branches": branches,
        })
        current = next_points
        modulus = next_modulus
        if not current:
            break
    return {
        "complete": True,
        "levels": levels,
        "survivors": [list(value) for value in current],
        "last_modulus": modulus,
    }


def verify_every_certificate(
    payload,
    points,
    system,
    variables,
    equations,
    primitives,
    deadline,
):
    polynomial333 = even_subgraph_polynomial(cubic(3, 3, 3, periodic=False))
    log_coefficients = formal_log_coefficients(polynomial333, 12)
    logs = {order: int(log_coefficients[order]) for order in (4, 6, 8, 10, 12)}
    assert logs == {4: 36, 6: 164, 8: 663, 10: 2280, 12: 972}
    assert payload["data"]["census"]["holdout_log_P_coefficients"] == {
        str(order): value for order, value in logs.items()
    }
    structure = transfer_structure()
    primitive_terms = [poly_terms(polynomial) for polynomial in primitives]
    derivative_table = [
        [derivative_terms(primitive_terms[row], column) for column in range(6)]
        for row in range(6)
    ]
    equation_terms = [poly_terms(sp.Poly(equation, *variables, domain=sp.ZZ)) for equation in equations]
    stored_records = payload["data"]["census"]["holdout_classification"]["point_records"]
    stored_by_point = {tuple(record["point13_mod5"]): record for record in stored_records}
    assert set(stored_by_point) == points
    recomputed_counts = Counter()

    for number, point in enumerate(sorted(points)):
        require_time(deadline)
        stored = stored_by_point[point]
        assert stored["point_id"] == f"P{number:04d}"
        residue8 = holdout_residue(point, 8, 5, logs, structure)
        if residue8:
            assert stored["disposition"] == "KILLED_H8_MOD5"
            assert stored["holdout_residues"] == {"8_mod_5": residue8}
            recomputed_counts[stored["disposition"]] += 1
            continue
        residue12 = holdout_residue(point, 12, 5, logs, structure)
        if residue12:
            assert stored["disposition"] == "KILLED_H12_MOD5"
            assert stored["holdout_residues"] == {"8_mod_5": 0, "12_mod_5": residue12}
            recomputed_counts[stored["disposition"]] += 1
            continue
        residue10 = holdout_residue(point, 10, 5, logs, structure)
        assert residue10 == 0
        jacobian = [
            [eval_terms(derivative_table[row][column], point, 5) for column in range(6)]
            for row in range(6)
        ]
        determinant_mod5 = determinant(jacobian)
        assert stored["jacobian_mod5"] == jacobian
        assert stored["jacobian_determinant_mod5"] == determinant_mod5
        if determinant_mod5:
            lifted = nonsingular_lift(primitive_terms, point, jacobian)
            values25 = reconstruct25(system, variables, lifted, 625)
            assert full42_zero(equation_terms, values25, 625)
            lifted_residue8 = holdout_residue(lifted, 8, 625, logs, structure)
            assert lifted_residue8 != 0
            assert stored["disposition"] == "KILLED_H8_MOD625_NONSINGULAR"
            assert stored["lift_solve6_mod625"] == list(lifted[:6])
            assert stored["holdout_residues"]["8_mod_625"] == lifted_residue8
        else:
            tree = singular_tree(
                point,
                jacobian,
                primitive_terms,
                system,
                variables,
                equation_terms,
                logs,
                structure,
                deadline,
            )
            assert not tree["survivors"]
            assert stored["disposition"] == "KILLED_BY_SINGULAR_LIFT_OBSTRUCTION"
            assert stored["singular_lift_tree"] == tree
        recomputed_counts[stored["disposition"]] += 1

    assert recomputed_counts == Counter({
        "KILLED_H8_MOD5": 2232,
        "KILLED_H12_MOD5": 664,
        "KILLED_H8_MOD625_NONSINGULAR": 45,
        "KILLED_BY_SINGULAR_LIFT_OBSTRUCTION": 19,
    })
    assert payload["data"]["classification_summary"]["point_disposition_counts"] == dict(
        sorted(recomputed_counts.items())
    )
    return stored_by_point, recomputed_counts


def bitmap(statuses):
    packed = bytearray((len(statuses) + 7) // 8)
    for index, status in enumerate(statuses):
        if status:
            packed[index // 8] |= 1 << (index % 8)
    return packed.hex()


def verify_all_section_rows(payload, points, stored_by_point):
    by_section = defaultdict(list)
    for point in points:
        by_section[point[6:]].append(point)
    rows = payload["data"]["classification_table"]
    assert len(rows) == TOTAL_SECTIONS
    statuses = []
    undecided_indices = []
    for index, (section, row) in enumerate(
        zip(itertools.product(UNIT_VALUES, repeat=7), rows, strict=True)
    ):
        assert row["lex_index"] == index
        assert row["section"] == list(section)
        assert row["orbit_size"] == 1
        section_points = sorted(by_section.get(section, []))
        assert row["solution_count_mod5"] == len(section_points)
        assert row["point_ids"] == [stored_by_point[point]["point_id"] for point in section_points]
        expected = "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT" if section_points else "UNDECIDED"
        assert row["status"] == expected
        statuses.append(bool(section_points))
        if not section_points:
            undecided_indices.append(index)
    assert sum(statuses) == 2_437
    assert len(undecided_indices) == 13_947
    summary = payload["data"]["classification_summary"]
    assert summary["status_counts"] == {
        "EMPTY_OVER_Q": 0,
        "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT": 2_437,
        "UNDECIDED": 13_947,
    }
    assert summary["remaining_undecided_lex_indices"] == undecided_indices
    assert summary["empty_over_Q_certificates"] == []
    stored_bitmap = payload["data"]["census"]["section_classification"]["status_bitmap"]
    assert stored_bitmap["hex"] == bitmap(statuses)
    return by_section, undecided_indices


def verify_provenance(payload):
    hashes = payload["meta"]["source_sha256"]
    for path in (E169, E170, E171, E139):
        relative = str(path.relative_to(ROOT))
        assert hashes[relative] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert payload["meta"]["floating_point_decides_no_statement"] is True
    assert payload["data"]["status"] == "COMPLETE_F5_UNIT_CENSUS_BRANCH_UNRESOLVED"
    assert all(check["passed"] for check in payload["data"]["checks"])


def main() -> None:
    started = time.process_time()
    deadline = started + CAP_SECONDS
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    verify_provenance(payload)
    system, variables, equations, _symbol, free, _non_diagonal, thin, primitives = build_raw_system()
    full_order, catalog_order, orbit_count, action_records = verify_group(equations, thin, free)
    assert payload["data"]["symmetry"]["cubic_point_group"]["full_lattice_group_order"] == full_order == 48
    assert payload["data"]["symmetry"]["cubic_point_group"]["construction_catalog_group_order"] == catalog_order == 16
    assert payload["data"]["symmetry"]["section_action"]["group_order"] == 1
    assert payload["data"]["symmetry"]["section_action"]["orbit_count"] == orbit_count == 16_384

    points = enumerate_unit_points(primitives, deadline)
    stored_points = {
        tuple(record["point13_mod5"])
        for record in payload["data"]["census"]["holdout_classification"]["point_records"]
    }
    assert points == stored_points
    verify_variety_fiber_counterexamples(payload, points, action_records)
    stored_by_point, mechanisms = verify_every_certificate(
        payload, points, system, variables, equations, primitives, deadline
    )
    by_section, undecided_indices = verify_all_section_rows(payload, points, stored_by_point)

    # Explicitly expose five independently replayed certificate types, exceeding
    # the requested minimum of three different types.
    examples = {}
    for disposition in mechanisms:
        examples[disposition] = next(
            record["point_id"]
            for record in stored_by_point.values()
            if record["disposition"] == disposition
        )
    undecided_section = tuple(payload["data"]["classification_table"][undecided_indices[0]]["section"])
    assert undecided_section not in by_section
    assert len(examples) == 4
    print("independent orbit count", orbit_count, "catalog group", catalog_order, "section group 1")
    print("independent complete census", len(points), "points in", len(by_section), "sections")
    print("certificate mechanisms", examples, "plus UNDECIDED no-point section", undecided_section)
    print("process_seconds", round(time.process_time() - started, 6), "cap", CAP_SECONDS)
    print("PASS")


if __name__ == "__main__":
    main()
