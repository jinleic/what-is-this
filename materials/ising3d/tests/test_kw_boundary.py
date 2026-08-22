"""Clean-room exact verifier for the wave-15 Kac--Ward boundary-lift census.

Recomputes, with exact integer arithmetic only (no producer import):

* the complete 5^6 = 15625 anchor-section census of the six primitive
  numerators (count 1553, unit stratum 3) and of the full 42-numerator system
  (set-equal, 1553);
* the full 42x13 Jacobian rank histogram mod 5 over all 1553 solutions and
  per-point ranks for representative rows (the three unit points, the boundary
  point q = (1,3,0,3,2,3), and the all-zero solve point);
* one exhaustive Hensel lift: q along the solve-6 chart to mod 625 by Newton
  iteration, all ten nonzero numerator residues at the lift, and the exact
  kept-zero coordinate;
* one holdout evaluation: the anchor unit lift, reconstructed into the
  25-weight chart, with order-8/10/12 residues by BOTH the walk-monomial
  symbolic route and the 108x108 matrix-power route (cross-validated, and
  equal to the e139 controls 250/55/266);
* the unit-theorem certificate: each thin relation has constant unit
  right-hand side and monomial denominator, and the nine denominator monomials
  jointly involve all sixteen free coordinates; every one of the 1550 boundary
  census points has a vanishing thin denominator, no unit point has.

Run: PYTHONPATH=src .venv/bin/python tests/test_kw_boundary.py
"""

from __future__ import annotations

import itertools
import json
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import (
    ALLOWED_DIRECTION_PAIRS,
    DIRECTION_GAUGE_TREE,
    _transition_indices,
    directed_edges,
    formal_log_coefficients,
    full_weight_finite_system,
)
from ising.exact_enumeration import even_subgraph_polynomial
from ising.lattices import cubic


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/boundary_lifts.json"
PRIOR_COMPONENTS = ROOT / "results/kac_ward/components.json"
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(set(itertools.permutations((2, 2, 1))) | set(itertools.permutations((3, 2, 1))))
)
FULL_SHAPES = ((3, 3, 2), (2, 2, 3), (2, 2, 2), (3, 2, 2), (2, 3, 2))
CONSTRUCTION_SHAPES = PRIMARY_SHAPES + FULL_SHAPES
HOLDOUT_SHAPE = (3, 3, 3)
SELECTED_NUMERATOR_INDICES = [28, 29, 32, 35, 38, 41]
DIAGONAL_NAMES = ["u_px_px", "u_py_py", "u_pz_pz"]
NONDIAGONAL_NAMES = [
    "u_py_px",
    "u_py_pz",
    "u_py_mz",
    "u_my_px",
    "u_my_mx",
    "u_pz_px",
    "u_pz_mx",
    "u_pz_py",
    "u_pz_my",
    "u_mz_px",
    "u_mz_mx",
    "u_mz_py",
    "u_mz_my",
]
FREE_NAMES = [
    "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx",
    "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx",
    "u_mz_py", "u_mz_my",
]
ANCHOR_HELD = [1, 1, 1, 1, 1, 4, 4]
PRIME = 5
MODULUS = 625
HOLDOUT_ORDERS = (4, 6, 8, 10, 12)
DIRECTIONS_3D = (
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
)


def branch_11111_equations(system):
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


def sparse_poly(expression, variables):
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return {
        tuple(map(int, exponents)): int(coefficient)
        for exponents, coefficient in polynomial.terms()
        if coefficient
    }


def eval_sparse(polynomial, values, modulus):
    total = 0
    for exponents, coefficient in polynomial.items():
        term = coefficient % modulus
        for position, exponent in enumerate(exponents):
            if exponent:
                term = term * pow(values[position], exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def substitute_held(sparse, held_values):
    out = {}
    for exponents, coefficient in sparse.items():
        factor = coefficient
        for position, exponent in enumerate(exponents[6:]):
            if exponent:
                factor *= pow(held_values[position], exponent)
        key = exponents[:6]
        out[key] = out.get(key, 0) + factor
    return {key: value for key, value in out.items() if value}


def eval_grid_mod5(sparse, nvars, radix=PRIME, modulus=PRIME):
    if nvars == 0:
        return [sum(sparse.values()) % modulus]
    buckets = {}
    for exponents, coefficient in sparse.items():
        bucket = buckets.setdefault(exponents[-1], {})
        reduced = exponents[:-1]
        bucket[reduced] = bucket.get(reduced, 0) + coefficient
    base = {e: eval_grid_mod5(sub, nvars - 1, radix, modulus) for e, sub in buckets.items()}
    out = [0] * (radix**nvars)
    stride = radix ** (nvars - 1)
    for value in range(radix):
        offset = value * stride
        for exponent, grid in base.items():
            factor = pow(value, exponent, modulus)
            if not factor:
                continue
            for i in range(stride):
                out[offset + i] = (out[offset + i] + factor * grid[i]) % modulus
    return out


def derive_sparse(sparse, position):
    out = {}
    for exponents, coefficient in sparse.items():
        if exponents[position]:
            key = list(exponents)
            exponent = key[position]
            key[position] -= 1
            key = tuple(key)
            out[key] = out.get(key, 0) + coefficient * exponent
    return {key: value for key, value in out.items() if value}


def rank_mod5(rows, columns):
    work = [list(row) for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, len(work)) if work[r][column] % PRIME), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, PRIME)
        work[rank] = [value * inverse % PRIME for value in work[rank]]
        for r in range(len(work)):
            if r != rank and work[r][column]:
                factor = work[r][column]
                work[r] = [(a - factor * b) % PRIME for a, b in zip(work[r], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def solve_linear_system(matrix, rhs, prime):
    size = len(matrix)
    augmented = [
        [entry % prime for entry in row] + [value % prime]
        for row, value in zip(matrix, rhs, strict=True)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column] % prime), None
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, prime)
        augmented[column] = [value * inverse % prime for value in augmented[column]]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (left - factor * right) % prime
                for left, right in zip(augmented[row], augmented[column], strict=True)
            ]
    return [row[-1] for row in augmented]


def hensel_lift_chart(slice_sparse6, initial, prime=PRIME, exponent=4):
    jacobian = [
        [eval_sparse(derive_sparse(row, column), initial, prime) for column in range(6)]
        for row in slice_sparse6
    ]
    values = list(initial)
    modulus = prime
    for _ in range(exponent - 1):
        residuals = [eval_sparse(row, values, modulus * prime) for row in slice_sparse6]
        assert all(residual % modulus == 0 for residual in residuals)
        correction = solve_linear_system(
            jacobian, [-(residual // modulus) % prime for residual in residuals], prime
        )
        assert correction is not None
        values = [v + modulus * d for v, d in zip(values, correction, strict=True)]
        modulus *= prime
    return values


def symbolic_trace_monomials(shape, max_length):
    lattice = cubic(*shape, periodic=False)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    n = len(edges)
    directions = [DIRECTIONS_3D.index(tuple(edge.direction)) for edge in edges]
    gauge = set(DIRECTION_GAUGE_TREE)
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge]
    weight_index = {pair: index for index, pair in enumerate(nongauge)}
    width = len(nongauge)
    rows = [{start: {(0,) * width: 1}} for start in range(n)]
    traces = {}
    for length in range(1, max_length + 1):
        next_rows = [dict() for _ in range(n)]
        for start in range(n):
            bucket = next_rows[start]
            for end, monomials in rows[start].items():
                for next_edge in transitions[end]:
                    key = (directions[end], directions[next_edge])
                    target = bucket.get(next_edge)
                    if key in gauge:
                        if target is None:
                            bucket[next_edge] = dict(monomials)
                        else:
                            for monomial, count in monomials.items():
                                merged = target.get(monomial, 0) + count
                                if merged:
                                    target[monomial] = merged
                        continue
                    position = weight_index[key]
                    if target is None:
                        expanded = {}
                        for monomial, count in monomials.items():
                            exponents = list(monomial)
                            exponents[position] += 1
                            expanded[tuple(exponents)] = count
                        bucket[next_edge] = expanded
                    else:
                        for monomial, count in monomials.items():
                            exponents = list(monomial)
                            exponents[position] += 1
                            monomial = tuple(exponents)
                            merged = target.get(monomial, 0) + count
                            if merged:
                                target[monomial] = merged
        rows = next_rows
        if length in HOLDOUT_ORDERS:
            trace = {}
            for start in range(n):
                for monomial, count in rows[start].get(start, {}).items():
                    trace[monomial] = trace.get(monomial, 0) + count
            traces[length] = trace
    return traces


def matrix_mul(left, right, modulus):
    size = len(left)
    output = [[0] * size for _ in range(size)]
    for row in range(size):
        for middle in range(size):
            if not left[row][middle]:
                continue
            factor = left[row][middle]
            for column in range(size):
                output[row][column] = (output[row][column] + factor * right[middle][column]) % modulus
    return output


def matrix_trace_power(nongauge_values, structure, length, modulus):
    n, transitions, directions = structure
    gauge = set(DIRECTION_GAUGE_TREE)
    matrix = [[0] * n for _ in range(n)]
    for edge in range(n):
        for next_edge in transitions[edge]:
            key = (directions[edge], directions[next_edge])
            weight = 1 if key in gauge else nongauge_values[
                [p for p in ALLOWED_DIRECTION_PAIRS if p not in gauge].index(key)
            ] % modulus
            matrix[edge][next_edge] = weight
    result = [[int(i == j) for j in range(n)] for i in range(n)]
    base = matrix
    factor = length
    while factor:
        if factor & 1:
            result = matrix_mul(result, base, modulus)
        base = matrix_mul(base, base, modulus)
        factor >>= 1
    return sum(result[i][i] for i in range(n)) % modulus


def main() -> None:
    started = time.process_time()

    # ---------------- raw build (independent of the producer) ----------------
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False)
    variables, equations = branch_11111_equations(system)
    assert len(equations) == 42
    symbol = {str(value): value for value in variables}
    pair_sym = dict(system.all_symbols)
    thin = thin_substitutions(symbol)
    free = [symbol[name] for name in FREE_NAMES]
    non_diagonal = [symbol[name] for name in NONDIAGONAL_NAMES]
    diagonal = {symbol[name]: 1 for name in DIAGONAL_NAMES}

    numerators13 = []
    for index in range(42):
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        numerator = sp.expand(sp.expand(numerator).subs(diagonal))
        _content, primitive = sp.Poly(numerator, *non_diagonal, domain=sp.ZZ).primitive()
        numerators13.append(sparse_poly(primitive.as_expr(), non_diagonal))
    nonzero_indices = [i for i, s in enumerate(numerators13) if s]
    assert nonzero_indices == [28, 29, 31, 32, 34, 35, 37, 38, 40, 41]
    primitives = [numerators13[i] for i in SELECTED_NUMERATOR_INDICES]
    assert [len(s) for s in primitives] == [33, 277, 155, 7, 155, 155]

    # ---------------- full 15625-point census ----------------
    prim6 = [substitute_held(s, ANCHOR_HELD) for s in primitives]
    grids_prim = [eval_grid_mod5(s, 6) for s in prim6]
    grids_all = [
        eval_grid_mod5(substitute_held(numerators13[i], ANCHOR_HELD), 6) for i in nonzero_indices
    ]
    census = [i for i in range(PRIME**6) if all(g[i] == 0 for g in grids_prim)]
    census42 = [i for i in range(PRIME**6) if all(g[i] == 0 for g in grids_all)]
    assert len(census) == 1553, len(census)
    assert set(census42) == set(census)
    solutions = [
        [(i // PRIME**k) % PRIME for k in range(6)] + list(ANCHOR_HELD) for i in census
    ]
    unit = [p for p in solutions if all(v for v in p[:6])]
    assert len(unit) == 3
    assert sorted(p[:6] for p in unit) == [
        [1, 2, 3, 3, 1, 1], [2, 4, 1, 3, 3, 1], [3, 2, 3, 1, 1, 1],
    ]
    print("census 1553 reproduced (grid 15625, unit 3, 42-numerator set-equal)")
    # ---------------- rank histogram (all 1553) and representative rows ------
    deriv_grids = {}
    for row, index in enumerate(nonzero_indices):
        for column in range(13):
            deriv_grids[(row, column)] = eval_grid_mod5(
                substitute_held(derive_sparse(numerators13[index], column), ANCHOR_HELD), 6
            )
    primitive_rows = sorted(nonzero_indices.index(i) for i in SELECTED_NUMERATOR_INDICES)
    histogram = {}
    ranks = {}
    for flat, point in zip(census, solutions, strict=True):
        rows = [
            [deriv_grids[(row, column)][flat] for column in range(13)]
            for row in range(len(nonzero_indices))
        ]
        tau = rank_mod5(rows, 13)
        rho = rank_mod5([rows[row] for row in primitive_rows], 13)
        assert tau == rho
        histogram[tau] = histogram.get(tau, 0) + 1
        ranks[tuple(point[:6])] = tau
    assert {str(k): v for k, v in sorted(histogram.items())} == {
        "0": 49, "1": 349, "2": 364, "3": 630, "4": 157, "6": 4,
    }
    assert ranks[(1, 2, 3, 3, 1, 1)] == 6
    assert ranks[(2, 4, 1, 3, 3, 1)] == 6
    assert ranks[(3, 2, 3, 1, 1, 1)] == 6
    assert ranks[(1, 3, 0, 3, 2, 3)] == 6
    assert ranks[(0, 0, 0, 0, 0, 0)] == 0  # the all-zero solve point (rank-0 stratum)
    print(f"rank histogram reproduced: {dict(sorted(histogram.items()))}")

    # ---------------- one exhaustive Hensel lift: the boundary point q --------
    q6 = [1, 3, 0, 3, 2, 3]
    lift_q = hensel_lift_chart(prim6, q6, PRIME, 4)
    assert all(v % PRIME == x % PRIME for v, x in zip(lift_q, q6, strict=True))
    lift_q13 = [int(v) for v in lift_q] + list(ANCHOR_HELD)
    assert lift_q13[:6] == [276, 278, 0, 288, 282, 83], lift_q13[:6]
    assert lift_q13[2] == 0, "the lifted u_py_mz must stay exactly zero"
    residues10 = {
        index: eval_sparse(numerators13[index], lift_q13, MODULUS)
        for index in nonzero_indices
    }
    assert all(value == 0 for value in residues10.values()), residues10
    denominators_q = {
        name: any(lift_q13[pos] % PRIME == 0 for pos in positions)
        for name, positions in {
            "u_mx_my": [3], "u_mx_py": [4, 0], "u_mx_pz": [10, 5],
            "u_my_mz": [11, 1, 8], "u_my_pz": [12, 2, 7],
        }.items()
    }
    assert denominators_q == {
        "u_mx_my": False, "u_mx_py": False, "u_mx_pz": False,
        "u_my_mz": False, "u_my_pz": True,
    }
    print(f"boundary lift q -> {lift_q13[:6]} mod 625, ten residues all zero, u_py_mz kept 0")

    # ---------------- one holdout evaluation: the anchor unit lift -----------
    p0 = [1, 2, 3, 3, 1, 1]
    lift_p0 = hensel_lift_chart(prim6, p0, PRIME, 4)
    assert [int(v) for v in lift_p0] == [226, 197, 98, 63, 161, 216]
    values = {
        name: int(v) for name, v in zip(NONDIAGONAL_NAMES, lift_p0 + ANCHOR_HELD, strict=True)
    }
    values.update({name: 1 for name in DIAGONAL_NAMES})
    free_values = [values[name] for name in FREE_NAMES]
    for name, expression in thin.items():
        numerator, denominator = sp.fraction(sp.cancel(expression))
        num = eval_sparse(sparse_poly(sp.Poly(numerator, *free, domain=sp.ZZ), free), free_values, MODULUS)
        den = eval_sparse(sparse_poly(sp.Poly(denominator, *free, domain=sp.ZZ), free), free_values, MODULUS)
        assert den % PRIME != 0
        values[str(name)] = num * pow(den, -1, MODULUS) % MODULUS
    values25 = [values[str(value)] for value in variables]
    sparse42 = [sparse_poly(expression, variables) for expression in equations]
    construction = [eval_sparse(sparse, values25, MODULUS) for sparse in sparse42]
    assert construction == [0] * 42

    traces = symbolic_trace_monomials(HOLDOUT_SHAPE, 12)
    log_p = formal_log_coefficients(
        even_subgraph_polynomial(cubic(*HOLDOUT_SHAPE, periodic=False)), 12
    )
    lattice = cubic(*HOLDOUT_SHAPE, periodic=False)
    edges = directed_edges(lattice)
    structure = (
        len(edges),
        _transition_indices(lattice),
        [DIRECTIONS_3D.index(tuple(edge.direction)) for edge in edges],
    )
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in set(DIRECTION_GAUGE_TREE)]
    walk_values = [values[str(pair_sym[pair])] for pair in nongauge]
    walk_residues = []
    matrix_residues = []
    for order in HOLDOUT_ORDERS:
        total = 0
        for monomial, count in traces[order].items():
            term = count % MODULUS
            for position, exponent in enumerate(monomial):
                if exponent:
                    term = term * pow(int(walk_values[position]), exponent, MODULUS) % MODULUS
            total = (total + term) % MODULUS
        required = -2 * order * log_p[order]
        walk_residues.append(
            (total - int(required.numerator) * pow(int(required.denominator) % MODULUS, -1, MODULUS))
            % MODULUS
        )
        trace = matrix_trace_power(walk_values, structure, order, MODULUS)
        matrix_residues.append(
            (trace - int(required.numerator) * pow(int(required.denominator) % MODULUS, -1, MODULUS))
            % MODULUS
        )
    assert walk_residues == matrix_residues
    assert walk_residues == [0, 0, 250, 55, 266]
    assert walk_residues[2:5] == [250, 55, 266]
    print(f"anchor holdout residues (walk == matrix): {walk_residues}")

    # ---------------- unit theorem certificate ----------------
    denominator_variables = {}
    for variable, expression in thin.items():
        _numerator, denominator = sp.fraction(sp.cancel(expression))
        monomial = sp.Poly(denominator, *free, domain=sp.ZZ)
        terms = monomial.terms()
        assert len(terms) == 1 and terms[0][1] in (1, -1)
        numerator = sp.fraction(sp.cancel(expression))[0]
        assert numerator.is_number and int(numerator) in (1, -1)
        denominator_variables[str(variable)] = [str(free[j]) for j, e in enumerate(terms[0][0]) if e]
    union = sorted({n for names in denominator_variables.values() for n in names})
    assert union == sorted(FREE_NAMES)
    solve_position = {name: j for j, name in enumerate(NONDIAGONAL_NAMES)}
    boundary_with_vanishing = 0
    unit_with_vanishing = 0
    for point in solutions:
        vanished = any(
            any(point[solve_position[name]] == 0 for name in names)
            for names in denominator_variables.values()
            if any(n in solve_position for n in names)
        )
        if all(v for v in point[:6]):
            unit_with_vanishing += int(vanished)
        else:
            boundary_with_vanishing += int(vanished)
    assert boundary_with_vanishing == 1550 and unit_with_vanishing == 0
    print("unit theorem: denominators cover all 16 free coordinates; 1550/1550 boundary, 0/3 unit")

    # ---------------- artifact agreement ----------------
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["census"]["primitive_solution_count"] == 1553
    assert payload["census"]["full42_numerator_solution_count"] == 1553
    assert payload["census"]["full42_set_equals_primitive_set"] is True
    assert payload["census"]["unit_stratum_count"] == 3
    assert payload["census"]["boundary_count"] == 1550
    assert sum(payload["census"]["zero_coordinate_histogram"].values()) == 1553
    assert payload["rank_analysis"]["tau_histogram_mod5"] == {
        "0": 49, "1": 349, "2": 364, "3": 630, "4": 157, "6": 4,
    }
    artifact_ranks = {tuple(row[0][:6]): row[1] for row in payload["census"]["points_with_tau"]}
    assert len(artifact_ranks) == 1553
    assert all(artifact_ranks[key] == ranks[key] for key in ranks)
    chart_by_point = {tuple(r["point13_mod5"][:6]): r for r in payload["chart_analysis"]["records"]}
    assert chart_by_point[(1, 3, 0, 3, 2, 3)]["nonsingular_charts_of_1716_primitive_rows"] == 588
    assert chart_by_point[(1, 2, 3, 3, 1, 1)]["nonsingular_charts_of_1716_primitive_rows"] == 1225
    assert chart_by_point[(2, 4, 1, 3, 3, 1)]["nonsingular_charts_of_1716_primitive_rows"] == 1023
    assert chart_by_point[(3, 2, 3, 1, 1, 1)]["nonsingular_charts_of_1716_primitive_rows"] == 1153
    assert chart_by_point[(1, 3, 0, 3, 2, 3)]["solve6_minor_determinant_mod5"] == 3
    lift_by_point = {tuple(r["point13_mod5"][:6]): r for r in payload["hensel_lifts"]["records"]}
    q_record = lift_by_point[(1, 3, 0, 3, 2, 3)]
    assert q_record["lift_solve6_mod625"] == [276, 278, 0, 288, 282, 83]
    assert q_record["chart_degenerate"] is True
    assert q_record["vanishing_thin_denominators"] == ["u_my_pz"]
    assert q_record["holdout"]["status"] == "NOT_EVALUABLE_CHART_DEGENERATE"
    assert q_record["holdout"]["lift_kept_zero_coordinates"] == ["u_py_mz"]
    assert all(v == 0 for v in q_record["all_ten_numerator_residues_mod625"].values())
    anchor_record = lift_by_point[(1, 2, 3, 3, 1, 1)]
    assert anchor_record["lift_solve6_mod625"] == [226, 197, 98, 63, 161, 216]
    assert anchor_record["holdout"]["walk_residues_mod625"] == [0, 0, 250, 55, 266]
    assert anchor_record["holdout"]["first_failing_order"] == 8
    assert payload["controls"]["anchor_lift_matches_e139"] is True
    assert payload["controls"]["anchor_holdout_matches_e139"] is True
    assert payload["controls"]["boundary_lift_all_ten_numerator_residues_zero_mod625"] is True
    assert payload["unit_theorem"]["union_covers_all_16_free_coordinates"] is True
    assert payload["unit_theorem"]["empirical_check"] == {
        "boundary_points_with_vanishing_thin_denominator": 1550,
        "boundary_points_total": 1550,
        "unit_points_with_vanishing_thin_denominator": 0,
    }
    assert all(check["passed"] for check in payload["checks"])
    prior = json.loads(PRIOR_COMPONENTS.read_text(encoding="utf-8"))
    assert prior["thinvariants"]["column"]["solutions_total_mod5"] == 3
    assert payload["controls"]["unit_lift_residues_orders_4_to_12"] == [
        [0, 0, 430, 480, 286],
        [0, 0, 250, 55, 266],
        [0, 0, 456, 350, 132],
    ]
    print(f"artifact agreement complete; event seconds {round(time.process_time() - started, 1)}")
    print("PASS")


if __name__ == "__main__":
    main()
