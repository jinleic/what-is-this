"""Boundary stratification of the Kac--Ward branch ``11111`` anchor-section census.

Wave-15 continuation of ``experiments/e139_kw_components.py`` /
``proofs/kw_components.md``.  e139 left the 1549 non-unit F_5 solutions of the
anchor thin-torus section ``(1,1,1,1,1,4,4)`` unanalysed.  This producer
certifies, with exact integer arithmetic only (mod 5 and mod 5^4 = 625):

* Full-grid census.  Over the complete ``5^6 = 15625`` anchor-section grid the
  six primitive quotient numerators have EXACTLY 1553 solutions (unit-chart
  subset: 3, reproducing e139); the full 42-equation thin-reduced numerator
  system (32 of the 42 numerators are identically zero on the thin torus; the
  ten nonzero ones are the five full shapes at orders 6 and 8) has the SAME
  1553-point zero set, set-equal to the primitive zero set.
* Rank stratification.  For each of the 1553 solutions the full 42x13 Jacobian
  rank mod 5 (13 variables = 6 solve + 7 held) is computed exactly; the rank
  histogram is recorded, and equals at every point the rank of the 6x13
  primitive Jacobian.  Exactly four points reach the maximal rank 6 (the three
  unit points and the boundary point q = (1,3,0,3,2,3)); the other 1549 are
  singular (rank <= 4), so ALL their C(13,6) = 1716 coordinate charts have
  every 6x6 minor singular.
* Chart census.  For the maximal-rank stratum every C(13,6) = 1716 coordinate
  charts are tested for a nonsingular 6x6 minor (primitive rows and full 42
  rows); exact counts recorded per point.
* Hensel lifts.  Every maximal-rank point is Hensel-lifted to Z/625 along the
  solve-6 chart by exact Newton iteration (mod 25, 125, 625); the lifted
  13-tuples satisfy the six primitives identically (Newton) and the four
  remaining nonzero numerators are evaluated mod 625 and recorded.
* Chart exclusion (the boundary verdict).  Each of the nine thin-torus
  relations is the equation ``w * (monomial in the 16 free coordinates) = +-1``
  with a constant unit right-hand side, and the union of the denominator
  monomials covers ALL SIXTEEN free coordinates; hence ANY point of ANY
  commutative ring satisfying the nine relations has all sixteen free
  coordinates invertible (unit theorem).  Consequently no lift of any of the
  1550 boundary solutions can be reconstructed into the 25-weight chart: each
  keeps a coordinate = 0 mod 5 (verified per point), the corresponding thin
  denominator is a zero divisor mod 625, and the 3x3x3 holdout is not
  evaluable there.  Empirically: every boundary point has >= 1 vanishing thin
  denominator, every unit point none.
* Controls.  The three unit lifts reproduce e139 exactly: L0 =
  (226,197,98,63,161,216) with construction residues [0]*42 and holdout
  residues (0,0,250,55,266) at orders 4,6,8,10,12, cross-validated between
  the walk-monomial symbolic route and the 108x108 matrix-power route.

No floating-point calculation decides any statement in this file; every residue
is an exact integer in 0..624.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from datetime import datetime, timezone
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
OUTPUT = ROOT / "results/kac_ward/boundary_lifts.json"
PRIOR_COMPONENTS = ROOT / "results/kac_ward/components.json"
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(set(itertools.permutations((2, 2, 1))) | set(itertools.permutations((3, 2, 1))))
)
FULL_SHAPES = ((3, 3, 2), (2, 2, 3), (2, 2, 2), (3, 2, 2), (2, 3, 2))
CONSTRUCTION_SHAPES = PRIMARY_SHAPES + FULL_SHAPES
HOLDOUT_SHAPE = (3, 3, 3)
SELECTED_NUMERATOR_INDICES = [28, 29, 32, 35, 38, 41]
EXPECTED_NONZERO_NUMERATOR_INDICES = [28, 29, 31, 32, 34, 35, 37, 38, 40, 41]
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
SLICE6_NAMES = NONDIAGONAL_NAMES[:6]
ANCHOR_HELD = [1, 1, 1, 1, 1, 4, 4]
PRIME = 5
MODULUS = PRIME**4  # 625
HOLDOUT_ORDERS = (4, 6, 8, 10, 12)
DIRECTIONS_3D = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
CHARTS_13_6 = list(itertools.combinations(range(13), 6))


def branch_11111_equations(system) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    fixed = set(DIRECTION_GAUGE_TREE)
    substitutions = {symbol: 1 for pair, symbol in system.all_symbols if pair in fixed}
    variables = [symbol for pair, symbol in system.all_symbols if pair not in fixed]
    equations = [sp.expand(expression.subs(substitutions)) for expression in system.equations]
    return variables, equations


def thin_torus_substitutions(symbol: dict[str, sp.Symbol]) -> dict[sp.Symbol, sp.Expr]:
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


def sparse_integer_poly(expression: sp.Expr, variables: list[sp.Symbol]) -> dict:
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return {
        tuple(map(int, exponents)): int(coefficient)
        for exponents, coefficient in polynomial.terms()
        if coefficient
    }


def eval_sparse_mod(polynomial: dict, values: list[int], modulus: int) -> int:
    total = 0
    for exponents, coefficient in polynomial.items():
        term = coefficient % modulus
        for position, exponent in enumerate(exponents):
            if exponent:
                term = term * pow(values[position], exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


# ---------------------------------------------------------------------------
# exact grid machinery over the F_5^6 anchor-section grid
# ---------------------------------------------------------------------------


def substitute_held(sparse: dict, held_values: list[int]) -> dict:
    """Fix the seven held coordinates; keep a sparse poly over the solve six."""
    out: dict = {}
    for exponents, coefficient in sparse.items():
        factor = coefficient
        for position, exponent in enumerate(exponents[6:]):
            if exponent:
                factor *= pow(held_values[position], exponent)
        key = exponents[:6]
        out[key] = out.get(key, 0) + factor
    return {key: value for key, value in out.items() if value}


def eval_grid_mod5(sparse: dict, nvars: int, radix: int = PRIME, modulus: int = PRIME) -> list[int]:
    """Evaluate a sparse integer poly at ALL radix^nvars points of the grid.

    Index convention: flat[i] with i = v0 + 5*v1 + ... + 5^(nvars-1)*v_(last),
    v_j in 0..4 the value of solve variable j (NONDIAGONAL_NAMES order).
    """

    if nvars == 0:
        return [sum(sparse.values()) % modulus]
    buckets: dict[int, dict] = {}
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


def derive_sparse(sparse: dict, position: int) -> dict:
    out: dict = {}
    for exponents, coefficient in sparse.items():
        if exponents[position]:
            key = list(exponents)
            exponent = key[position]
            key[position] -= 1
            key = tuple(key)
            out[key] = out.get(key, 0) + coefficient * exponent
    return {key: value for key, value in out.items() if value}


def rank_mod5(rows: list[list[int]], columns: int) -> int:
    """Exact Gaussian-elimination rank over F_5."""
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


def det_mod5(matrix: list[list[int]]) -> int:
    size = len(matrix)
    work = [row[:] for row in matrix]
    det = 1
    for column in range(size):
        pivot = next((r for r in range(column, size) if work[r][column] % PRIME), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            det = -det
        det = det * work[column][column] % PRIME
        inverse = pow(work[column][column], -1, PRIME)
        work[column] = [value * inverse % PRIME for value in work[column]]
        for r in range(column + 1, size):
            if work[r][column]:
                factor = work[r][column]
                work[r] = [(a - factor * b) % PRIME for a, b in zip(work[r], work[column])]
    return det % PRIME


def decode_index(index: int, nvars: int = 6) -> list[int]:
    return [(index // PRIME**k) % PRIME for k in range(nvars)]


# ---------------------------------------------------------------------------
# exact Hensel lifting (Newton iteration mod 25, 125, 625 in integers)
# ---------------------------------------------------------------------------


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


def hensel_lift_chart(
    slice_sparse6: list[dict],
    initial: list[int],
    prime: int = PRIME,
    exponent: int = 4,
):
    """Newton lift of a nonsingular F_p zero of six polynomials in six unknowns.

    ``slice_sparse6[r]`` is the sparse dict of equation r in the six chart
    variables.  Iterates exactly in integers through moduli p, p^2, p^3, p^4
    (Newton steps to mod 25, 125, 625).  Returns (values mod p**exponent,
    jacobian mod p) or (None, jacobian).
    """

    jacobian = [
        [eval_sparse_mod(derive_sparse(row, column), initial, prime) for column in range(6)]
        for row in slice_sparse6
    ]
    values = list(initial)
    modulus = prime
    for _ in range(exponent - 1):
        residuals = [eval_sparse_mod(row, values, modulus * prime) for row in slice_sparse6]
        assert all(residual % modulus == 0 for residual in residuals), "Hensel invariant failed"
        correction = solve_linear_system(
            jacobian,
            [-(residual // modulus) % prime for residual in residuals],
            prime,
        )
        if correction is None:
            return None, jacobian
        values = [value + modulus * delta for value, delta in zip(values, correction, strict=True)]
        modulus *= prime
    return values, jacobian


# ---------------------------------------------------------------------------
# transfer-matrix routes for the 3x3x3 open box (e139 conventions)
# ---------------------------------------------------------------------------


def transfer_matrix_structure(shape: tuple[int, int, int]):
    lattice = cubic(*shape, periodic=False)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    directions = [DIRECTIONS_3D.index(tuple(edge.direction)) for edge in edges]
    gauge = set(DIRECTION_GAUGE_TREE)
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge]
    return {
        "n": len(edges),
        "transitions": transitions,
        "directions": directions,
        "nongauge": nongauge,
    }


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


def trace_power_mod(weights25: list[int], structure: dict, length: int, modulus: int, nongauge):
    """Matrix-power route: Tr(T^length) mod modulus."""
    n = structure["n"]
    transitions = structure["transitions"]
    directions = structure["directions"]
    matrix = [[0] * n for _ in range(n)]
    for edge in range(n):
        for next_edge in transitions[edge]:
            key = (directions[edge], directions[next_edge])
            weight = 1 if key in set(DIRECTION_GAUGE_TREE) else weights25[nongauge.index(key)] % modulus
            matrix[edge][next_edge] = weight
    power_matrix = [[int(edge == column) for column in range(n)] for edge in range(n)]
    base = matrix
    factor = length
    while factor:
        if factor & 1:
            power_matrix = matrix_mul(power_matrix, base, modulus)
        base = matrix_mul(base, base, modulus)
        factor >>= 1
    return sum(power_matrix[edge][edge] for edge in range(n)) % modulus


def symbolic_trace_monomials(shape: tuple[int, int, int], max_length: int) -> dict[int, dict]:
    """Exact Tr(T^length) as a monomial multiset over the nongauge weights."""
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
    traces: dict[int, dict] = {}
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
            trace: dict = {}
            for start in range(n):
                for monomial, count in rows[start].get(start, {}).items():
                    trace[monomial] = trace.get(monomial, 0) + count
            traces[length] = trace
    return traces


def holdout_residue_walk(
    trace_monomials: dict,
    nongauge_values25: list[int],
    log_p: dict[int, Fraction],
    order: int,
    modulus: int,
) -> int:
    """Exact equation residue: Tr(T^order) - (-2*order*logP[order]) mod modulus."""
    total = 0
    for monomial, count in trace_monomials.items():
        term = count % modulus
        for position, exponent in enumerate(monomial):
            if exponent:
                term = term * pow(int(nongauge_values25[position]), exponent, modulus) % modulus
        total = (total + term) % modulus
    required = -2 * order * log_p[order]
    numerator = int(required.numerator)
    denominator = int(required.denominator) % modulus
    assert denominator, f"logP denominator not a 5-adic unit at order {order}"
    return (total - numerator * pow(denominator, -1, modulus)) % modulus


def main() -> None:
    started = time.process_time()

    def note(msg):
        print(msg, flush=True)
        return time.process_time()

    t_stage = started
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False)
    variables, equations = branch_11111_equations(system)
    labels = [list(label) for label in system.labels]
    assert len(equations) == 42 and len(labels) == 42
    symbol = {str(value): value for value in variables}
    pair_sym = dict(system.all_symbols)
    thin = thin_torus_substitutions(symbol)
    free = [symbol[name] for name in FREE_NAMES]
    non_diagonal = [symbol[name] for name in NONDIAGONAL_NAMES]
    assert [str(value) for value in non_diagonal] == NONDIAGONAL_NAMES
    diagonal = {symbol[name]: 1 for name in DIAGONAL_NAMES}
    sparse42 = [sparse_integer_poly(expression, variables) for expression in equations]

    # ---- all 42 thin-reduced primitive numerators in the 13 non-diagonal vars ----
    numerators13_sparse: list[dict] = []
    numerator_poly_objects: dict[int, sp.Poly] = {}
    numerator_catalogue = []
    for index in range(42):
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        numerator = sp.expand(sp.expand(numerator).subs(diagonal))
        _content, primitive = sp.Poly(numerator, *non_diagonal, domain=sp.ZZ).primitive()
        numerators13_sparse.append(sparse_integer_poly(primitive.as_expr(), non_diagonal))
        if numerators13_sparse[-1]:
            numerator_poly_objects[index] = primitive
        numerator_catalogue.append(
            {
                "index": index,
                "label": [labels[index][0], labels[index][1]],
                "terms": len(numerators13_sparse[-1]),
                "identically_zero_on_thin_torus": not numerators13_sparse[-1],
            }
        )
    nonzero_indices = [i for i, s in enumerate(numerators13_sparse) if s]
    assert nonzero_indices == EXPECTED_NONZERO_NUMERATOR_INDICES, nonzero_indices
    primitive_sparse13 = [numerators13_sparse[i] for i in SELECTED_NUMERATOR_INDICES]
    assert [len(s) for s in primitive_sparse13] == [33, 277, 155, 7, 155, 155]
    note(f"stage numerators {round(time.process_time() - t_stage, 1)}s")
    t_stage = time.process_time()

    # ---- full-grid census over F_5^6 ----
    prim6 = [substitute_held(s, ANCHOR_HELD) for s in primitive_sparse13]
    grids_prim = [eval_grid_mod5(s, 6) for s in prim6]
    grids_all = [
        eval_grid_mod5(substitute_held(numerators13_sparse[i], ANCHOR_HELD), 6)
        for i in nonzero_indices
    ]
    census_primitive = [i for i in range(PRIME**6) if all(g[i] == 0 for g in grids_prim)]
    census_all42 = [i for i in range(PRIME**6) if all(g[i] == 0 for g in grids_all)]
    assert len(census_primitive) == 1553, len(census_primitive)
    assert set(census_all42) == set(census_primitive)
    solutions = [decode_index(i) + list(ANCHOR_HELD) for i in census_primitive]
    unit_solutions = [p for p in solutions if all(v for v in p[:6])]
    assert len(unit_solutions) == 3
    unit_sorted = sorted(unit_solutions)
    assert unit_sorted == [
        [1, 2, 3, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4],
        [2, 4, 1, 3, 3, 1, 1, 1, 1, 1, 1, 4, 4],
        [3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4],
    ]
    zero_histogram: dict[str, int] = {}
    for point in solutions:
        key = ",".join(str(j) for j in range(6) if not point[j]) or "none"
        zero_histogram[key] = zero_histogram.get(key, 0) + 1
    note(
        f"stage census {round(time.process_time() - t_stage, 1)}s: "
        f"{len(census_primitive)} solutions, 3 unit"
    )
    t_stage = time.process_time()

    # ---- 42x13 Jacobian rank mod 5 at every solution ----
    deriv_grids = {}
    for row, index in enumerate(nonzero_indices):
        for column in range(13):
            deriv_grids[(row, column)] = eval_grid_mod5(
                substitute_held(derive_sparse(numerators13_sparse[index], column), ANCHOR_HELD),
                6,
            )
    primitive_rows = sorted(nonzero_indices.index(i) for i in SELECTED_NUMERATOR_INDICES)
    n_nonzero = len(nonzero_indices)
    tau_histogram: dict[int, int] = {}
    per_point: list[dict] = []
    for flat, point in zip(census_primitive, solutions, strict=True):
        full_rows = [
            [deriv_grids[(row, column)][flat] for column in range(13)]
            for row in range(n_nonzero)
        ]
        tau = rank_mod5(full_rows, 13)
        prim_rows = [full_rows[row] for row in primitive_rows]
        rho = rank_mod5(prim_rows, 13)
        assert tau == rho, (point, tau, rho)
        tau_histogram[tau] = tau_histogram.get(tau, 0) + 1
        per_point.append({"flat": flat, "point13": point, "tau": tau, "rho": rho})
    top_rank = max(tau_histogram)
    top_entries = [pp for pp in per_point if pp["tau"] == top_rank]
    assert top_rank == 6 and len(top_entries) == 4
    note(
        f"stage ranks {round(time.process_time() - t_stage, 1)}s: histogram "
        f"{dict(sorted(tau_histogram.items()))}"
    )
    t_stage = time.process_time()

    # ---- exhaustive C(13,6) chart test on the maximal-rank stratum ----
    chart_records = []
    for pp in top_entries:
        full_rows = [
            [deriv_grids[(row, column)][pp["flat"]] for column in range(13)]
            for row in range(n_nonzero)
        ]
        prim_rows = [full_rows[row] for row in primitive_rows]
        solve_det = det_mod5([row[:6] for row in prim_rows])
        nonsingular_primitive_charts = 0
        nonsingular_full_charts = 0
        first_chart = None
        for chart in CHARTS_13_6:
            det = det_mod5([[row[c] for c in chart] for row in prim_rows])
            if det:
                nonsingular_primitive_charts += 1
                if first_chart is None:
                    first_chart = list(chart)
            if rank_mod5([[row[c] for c in chart] for row in full_rows], 6) == 6:
                nonsingular_full_charts += 1
        assert solve_det != 0
        assert nonsingular_primitive_charts >= 1 and nonsingular_full_charts >= 1
        chart_records.append(
            {
                "point13_mod5": pp["point13"],
                "tau_42x13_mod5": top_rank,
                "rho_6x13_mod5": pp["rho"],
                "solve6_minor_determinant_mod5": solve_det,
                "nonsingular_charts_of_1716_primitive_rows": nonsingular_primitive_charts,
                "nonsingular_charts_of_1716_full42_rows": nonsingular_full_charts,
                "first_nonsingular_chart_positions": first_chart,
            }
        )
    note(f"stage charts {round(time.process_time() - t_stage, 1)}s")
    t_stage = time.process_time()

    # ---- unit theorem: denominator monomials of the nine thin relations ----
    denominator_variables: dict[str, list[str]] = {}
    for variable, expression in thin.items():
        _numerator, denominator = sp.fraction(sp.cancel(expression))
        monomial = sp.Poly(denominator, *free, domain=sp.ZZ)
        terms = monomial.terms()
        assert len(terms) == 1 and terms[0][1] in (1, -1), "thin denominator is not a monomial"
        reduced_relation = sp.cancel(expression)
        relation_numerator = sp.fraction(reduced_relation)[0]
        assert relation_numerator.is_number and int(relation_numerator) in (1, -1), (
            "thin relation right-hand side is not the constant unit +-1"
        )
        denominator_variables[str(variable)] = [
            str(free[j]) for j, e in enumerate(terms[0][0]) if e
        ]
    union_denominator_variables = sorted(
        {name for names in denominator_variables.values() for name in names}
    )
    assert union_denominator_variables == sorted(FREE_NAMES)

    denominator_positions13: dict[str, list[int]] = {
        name: [
            j for j, coord in enumerate(NONDIAGONAL_NAMES) if coord in denominator_variables[name]
        ]
        for name in denominator_variables
    }

    def vanishing_denominators(point13: list[int]) -> list[str]:
        return [
            name
            for name, positions in denominator_positions13.items()
            if any(point13[j] % PRIME == 0 for j in positions)
        ]

    boundary_with_vanishing = 0
    unit_with_vanishing = 0
    for point in solutions:
        vanished = vanishing_denominators(point)
        if all(v for v in point[:6]):
            unit_with_vanishing += int(bool(vanished))
        else:
            boundary_with_vanishing += int(bool(vanished))
    boundary_count = len(solutions) - len(unit_solutions)
    assert boundary_with_vanishing == boundary_count == 1550
    assert unit_with_vanishing == 0
    note(f"stage unit theorem {round(time.process_time() - t_stage, 1)}s")
    t_stage = time.process_time()

    # ---- Hensel lifts of the maximal-rank stratum along the solve-6 chart ----
    lift_records = []
    for record in chart_records:
        point = record["point13_mod5"]
        lift, _jacobian = hensel_lift_chart(prim6, point[:6], PRIME, 4)
        assert lift is not None
        assert all(
            value % PRIME == initial % PRIME
            for value, initial in zip(lift, point[:6], strict=True)
        )
        lift13 = [int(v) for v in lift] + list(ANCHOR_HELD)
        residues10 = {
            index: eval_sparse_mod(numerators13_sparse[index], lift13, MODULUS)
            for index in nonzero_indices
        }
        nonzero_residues = {index: value for index, value in residues10.items() if value}
        vanished = vanishing_denominators(lift13)
        unit_lift = not any(v % PRIME == 0 for v in lift13)
        primitive_residues = {
            index: residues10[index] for index in SELECTED_NUMERATOR_INDICES
        }
        assert all(value == 0 for value in primitive_residues.values())
        lift_records.append(
            {
                "point13_mod5": point,
                "chart_positions": [0, 1, 2, 3, 4, 5],
                "chart_variable_names": SLICE6_NAMES,
                "pinned_values": list(ANCHOR_HELD),
                "solve6_minor_determinant_mod5": record["solve6_minor_determinant_mod5"],
                "lift_solve6_mod625": [int(v) for v in lift],
                "lift13_mod625": lift13,
                "all_ten_numerator_residues_mod625": {str(k): v for k, v in residues10.items()},
                "nonprimitive_numerator_residues_mod625": {
                    str(k): v
                    for k, v in residues10.items()
                    if k not in SELECTED_NUMERATOR_INDICES
                },
                "unit_lift": unit_lift,
                "chart_degenerate": bool(vanished),
                "vanishing_thin_denominators": vanished,
            }
        )
    note(f"stage hensel {round(time.process_time() - t_stage, 1)}s")
    t_stage = time.process_time()

    # ---- holdout machinery (only evaluable on unit lifts) ----
    structural = transfer_matrix_structure(HOLDOUT_SHAPE)
    gauge = set(DIRECTION_GAUGE_TREE)
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge]

    thin_spec = {}
    for variable, expression in thin.items():
        numerator, denominator = sp.fraction(sp.cancel(expression))
        thin_spec[str(variable)] = (
            sparse_integer_poly(sp.Poly(numerator, *free, domain=sp.ZZ), free),
            sparse_integer_poly(sp.Poly(denominator, *free, domain=sp.ZZ), free),
        )

    def reconstruct_values25(point13_mod625: list[int]) -> dict[str, int]:
        values = {
            name: int(value) % MODULUS
            for name, value in zip(NONDIAGONAL_NAMES, point13_mod625, strict=True)
        }
        values.update({name: 1 for name in DIAGONAL_NAMES})
        free_values = [values[name] for name in FREE_NAMES]
        for name, (numerator, denominator) in thin_spec.items():
            num = eval_sparse_mod(numerator, free_values, MODULUS)
            den = eval_sparse_mod(denominator, free_values, MODULUS)
            assert den % PRIME != 0, f"thin denominator for {name} is not a 5-adic unit"
            values[name] = num * pow(den, -1, MODULUS) % MODULUS
        return values

    symbolic_traces = symbolic_trace_monomials(HOLDOUT_SHAPE, 12)
    polynomial = even_subgraph_polynomial(cubic(*HOLDOUT_SHAPE, periodic=False))
    log_p = formal_log_coefficients(polynomial, 12)
    log_p_component = {
        k: [int(log_p[k].numerator), int(log_p[k].denominator)] for k in HOLDOUT_ORDERS
    }

    def holdout_both_routes(values25: dict[str, int]) -> list[int]:
        nongauge_values = [values25[str(pair_sym[pair])] for pair in nongauge]
        walk = [
            holdout_residue_walk(symbolic_traces[order], nongauge_values, log_p, order, MODULUS)
            for order in HOLDOUT_ORDERS
        ]
        matrix = []
        for order in HOLDOUT_ORDERS:
            trace = trace_power_mod(nongauge_values, structural, order, MODULUS, nongauge)
            required = -2 * order * log_p[order]
            matrix.append(
                (
                    trace
                    - int(required.numerator) * pow(int(required.denominator) % MODULUS, -1, MODULUS)
                )
                % MODULUS
            )
        assert walk == matrix, "walk/matrix holdout mismatch"
        return walk

    for entry in lift_records:
        if not entry["unit_lift"]:
            entry["holdout"] = {
                "status": "NOT_EVALUABLE_CHART_DEGENERATE",
                "detail": (
                    "the lift keeps a solve coordinate = 0 mod 5, so a thin-torus "
                    "denominator is a zero divisor mod 625 and no 25-weight chart "
                    "point exists over this lift (unit theorem)"
                ),
                "orders": list(HOLDOUT_ORDERS),
                "first_failing_order": None,
                "walk_residues_mod625": None,
                "matrix_residues_mod625": None,
            }
            continue
        values25 = reconstruct_values25(entry["lift13_mod625"])
        value_list25 = [values25[str(value)] for value in variables]
        construction = [eval_sparse_mod(sparse, value_list25, MODULUS) for sparse in sparse42]
        assert construction == [0] * 42
        residues = holdout_both_routes(values25)
        failing = [
            order for order, residue in zip(HOLDOUT_ORDERS, residues, strict=True) if residue
        ]
        entry["lifted_values25_mod625"] = {
            name: values25[name] for name in NONDIAGONAL_NAMES + sorted(thin_spec)
        }
        entry["construction_residues_mod625"] = construction
        entry["holdout"] = {
            "status": "EVALUATED_UNIT_LIFT",
            "orders": list(HOLDOUT_ORDERS),
            "walk_residues_mod625": residues,
            "matrix_residues_mod625": residues,
            "first_failing_order": next(iter(failing), None),
        }
    note(f"stage holdout {round(time.process_time() - t_stage, 1)}s")

    # ---- controls against the wave-14 artifact ----
    prior = json.loads(PRIOR_COMPONENTS.read_text(encoding="utf-8"))
    prior_lifts = prior["hensel_components"]["lift_records"]
    prior_map = {tuple(record["point13_mod5"]): record for record in prior_lifts}
    for entry in lift_records:
        key = tuple(entry["point13_mod5"])
        if key in prior_map:
            assert entry["lift_solve6_mod625"] == prior_map[key]["lift_solve6_mod625"]
            assert entry["holdout"]["walk_residues_mod625"] == prior_map[key][
                "holdout_residues_symbolic_mod625"
            ]
    unit_lift_entries = [e for e in lift_records if e["unit_lift"]]
    assert len(unit_lift_entries) == 3
    anchor_entry = next(e for e in unit_lift_entries if e["point13_mod5"][:6] == [1, 2, 3, 3, 1, 1])
    assert anchor_entry["lift_solve6_mod625"] == [226, 197, 98, 63, 161, 216]
    assert anchor_entry["holdout"]["walk_residues_mod625"] == [0, 0, 250, 55, 266]
    boundary_lift_entries = [e for e in lift_records if not e["unit_lift"]]
    assert len(boundary_lift_entries) == 1
    boundary_lift_entries[0]["holdout"]["lift_kept_zero_coordinates"] = [
        NONDIAGONAL_NAMES[j] for j in range(6) if boundary_lift_entries[0]["lift13_mod625"][j] % PRIME == 0
    ]

    census_digest = hashlib.sha256(
        json.dumps([[pp["point13"], pp["tau"]] for pp in per_point], separators=(",", ":")).encode()
    ).hexdigest()
    singular_count = sum(count for rank, count in tau_histogram.items() if rank < 6)

    payload = {
        "schema_version": 1,
        "headline": (
            "[COMPUTATION] Over the FULL 5^6 = 15625 anchor-section grid of Kac--Ward branch "
            "11111 the six primitive numerators have exactly 1553 F_5 solutions (unit stratum "
            "exactly 3, reproducing wave-14); the full 42-numerator system (32 numerators "
            "identically zero on the thin torus, ten nonzero) has the same 1553-point zero "
            f"set. The 42x13 Jacobian rank histogram mod 5 is {dict(sorted(tau_histogram.items()))} "
            "as a map rank->count: exactly four points reach the maximal rank 6 (the three unit "
            "points and the boundary point (1,3,0,3,2,3)); the other 1549 are singular (rank "
            "<= 4, so every one of their C(13,6)=1716 coordinate charts has all 6x6 minors "
            "singular). [THEOREM] Unit theorem: every thin-torus relation is the equation "
            "w*(monomial in the 16 free coordinates) = +-1 with constant unit right-hand side, "
            "and the denominator monomials jointly involve all 16 free coordinates, so ANY "
            "point of ANY commutative ring satisfying the nine relations has all sixteen free "
            "coordinates invertible; hence no boundary (non-unit) census point, smooth or "
            "singular, is the reduction of ANY thin-torus chart point over Z/625 or Q_5 -- the "
            "3x3x3 holdout is not evaluable there and no new components arise. The three unit "
            "lifts remain the complete chart-point set over the anchor section and all fail "
            "the order-8/10/12 3x3x3 holdout mod 625 (e139 controls reproduced exactly)."
        ),
        "status": "UNRESOLVED",
        "claim_tag": "[THEOREM]",
        "section": {
            "held_variable_names": NONDIAGONAL_NAMES[6:],
            "held_values_mod5": ANCHOR_HELD,
            "solve_variable_names": SLICE6_NAMES,
        },
        "numerator_catalogue": {
            "nondiagonal_variable_order": NONDIAGONAL_NAMES,
            "rows": numerator_catalogue,
            "identically_zero_indices": [
                row["index"] for row in numerator_catalogue if row["identically_zero_on_thin_torus"]
            ],
            "nonzero_indices": nonzero_indices,
            "selected_primitive_indices": SELECTED_NUMERATOR_INDICES,
            "primitive_term_counts": [len(s) for s in primitive_sparse13],
        },
        "census": {
            "grid_size": PRIME**6,
            "primitive_solution_count": len(census_primitive),
            "full42_numerator_solution_count": len(census_all42),
            "full42_set_equals_primitive_set": set(census_all42) == set(census_primitive),
            "unit_stratum_count": len(unit_solutions),
            "unit_solutions": unit_sorted,
            "boundary_count": boundary_count,
            "zero_coordinate_histogram": zero_histogram,
            "census_digest_sha256": census_digest,
            "points_with_tau": [[pp["point13"], pp["tau"]] for pp in per_point],
        },
        "rank_analysis": {
            "jacobian": (
                "42 x 13 (rows = thin-reduced primitive numerators, 32 of them identically "
                "zero; columns = the 13 non-diagonal variables = 6 solve + 7 held)"
            ),
            "tau_histogram_mod5": {str(rank): count for rank, count in sorted(tau_histogram.items())},
            "per_point_tau_equals_rho": True,
            "max_rank": top_rank,
            "max_rank_points": len(top_entries),
            "singular_count": singular_count,
            "singular_stratum_certificate": (
                "rank < 6 at each of these points forces every 6x6 minor of every one of the "
                "C(13,6) = 1716 coordinate charts to be singular mod 5"
            ),
            "charts_tested_per_max_rank_point": len(CHARTS_13_6),
        },
        "chart_analysis": {
            "charts_enumerated": len(CHARTS_13_6),
            "records": chart_records,
        },
        "hensel_lifts": {
            "prime": PRIME,
            "modulus": MODULUS,
            "newton_steps": "mod 25, 125, 625 (exact integers)",
            "records": lift_records,
        },
        "unit_theorem": {
            "statement": (
                "Each of the nine thin-torus substitutions u = c/(monomial in the 16 free "
                "coordinates), c = +-1, is the defining equation u * monomial = c of the "
                "chart. In any commutative ring, a product equal to a unit forces every "
                "factor to be a unit; the nine denominator monomials jointly involve all "
                "sixteen free coordinates, so every chart point has all sixteen free "
                "coordinates invertible. In particular no point with any non-diagonal "
                "coordinate = 0 mod 5 -- every one of the 1550 boundary census points -- is "
                "the reduction of a chart point over Z/625 or Q_5, whatever its smoothness."
            ),
            "relations": [
                {
                    "name": name,
                    "denominator_variables": denominator_variables[name],
                    "denominator_positions_in_13": denominator_positions13[name],
                }
                for name in sorted(denominator_variables)
            ],
            "union_covers_all_16_free_coordinates": union_denominator_variables == sorted(FREE_NAMES),
            "empirical_check": {
                "boundary_points_with_vanishing_thin_denominator": boundary_with_vanishing,
                "boundary_points_total": boundary_count,
                "unit_points_with_vanishing_thin_denominator": unit_with_vanishing,
            },
        },
        "controls": {
            "prior_artifact": "results/kac_ward/components.json",
            "anchor_lift_matches_e139": anchor_entry["lift_solve6_mod625"] == [226, 197, 98, 63, 161, 216],
            "anchor_holdout_matches_e139": anchor_entry["holdout"]["walk_residues_mod625"]
            == [0, 0, 250, 55, 266],
            "unit_lift_residues_orders_4_to_12": [
                entry["holdout"]["walk_residues_mod625"] for entry in unit_lift_entries
            ],
            "boundary_lift_all_ten_numerator_residues_zero_mod625": all(
                all(value == 0 for value in entry["all_ten_numerator_residues_mod625"].values())
                and all(
                    value == 0
                    for value in entry["nonprimitive_numerator_residues_mod625"].values()
                )
                for entry in boundary_lift_entries
            ),
        },
        "checks": [
            {
                "name": "full_grid_census_1553",
                "passed": len(census_primitive) == 1553,
                "detail": {"count": len(census_primitive), "grid": 15625},
            },
            {
                "name": "full42_numerator_census_set_equal",
                "passed": set(census_all42) == set(census_primitive) and len(census_all42) == 1553,
                "detail": {"count42": len(census_all42)},
            },
            {
                "name": "unit_stratum_exactly_three",
                "passed": len(unit_solutions) == 3,
                "detail": {"unit_solutions": unit_sorted},
            },
            {
                "name": "rank_histogram_recorded",
                "passed": sum(tau_histogram.values()) == 1553 and top_rank == 6,
                "detail": {"histogram": {str(r): c for r, c in sorted(tau_histogram.items())}},
            },
            {
                "name": "maximal_rank_stratum_size_four",
                "passed": len(top_entries) == 4,
                "detail": {"points": [pp["point13"] for pp in top_entries]},
            },
            {
                "name": "every_max_rank_point_lifted",
                "passed": sorted(
                    (r["point13_mod5"] for r in lift_records), key=lambda p: tuple(p)
                )
                == sorted((pp["point13"] for pp in top_entries), key=lambda p: tuple(p)),
                "detail": {
                    "lifted": [r["point13_mod5"] for r in lift_records],
                    "all_ten_numerator_residues_zero_mod625": all(
                        all(value == 0 for value in r["all_ten_numerator_residues_mod625"].values())
                        for r in lift_records
                    ),
                },
            },
            {
                "name": "boundary_chart_exclusion",
                "passed": boundary_with_vanishing == 1550 and unit_with_vanishing == 0,
                "detail": {
                    "boundary_with_vanishing_denominator": boundary_with_vanishing,
                    "unit_with_vanishing_denominator": unit_with_vanishing,
                },
            },
            {
                "name": "nonunit_lifts_chart_degenerate",
                "passed": all(
                    r["chart_degenerate"]
                    and r["holdout"]["status"] == "NOT_EVALUABLE_CHART_DEGENERATE"
                    for r in lift_records
                    if not r["unit_lift"]
                ),
                "detail": {
                    "nonunit_lifts": [r["point13_mod5"] for r in lift_records if not r["unit_lift"]],
                    "kept_zero_coordinates": boundary_lift_entries[0]["holdout"][
                        "lift_kept_zero_coordinates"
                    ],
                },
            },
            {
                "name": "unit_lifts_reproduce_e139",
                "passed": anchor_entry["lift_solve6_mod625"] == [226, 197, 98, 63, 161, 216]
                and anchor_entry["holdout"]["walk_residues_mod625"] == [0, 0, 250, 55, 266]
                and all(r["construction_residues_mod625"] == [0] * 42 for r in unit_lift_entries),
                "detail": {
                    "anchor_lift": anchor_entry["lift_solve6_mod625"],
                    "anchor_holdout": anchor_entry["holdout"]["walk_residues_mod625"],
                },
            },
        ],
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact integer modular arithmetic over Z/5Z and Z/625Z: full 42-numerator "
                "thin-torus reduction, complete F_5^6 grid census by recursive grid "
                "substitution, exact 42x13 Jacobian ranks by Gaussian elimination over F_5, "
                "exhaustive C(13,6)=1716 chart minor census on the maximal-rank stratum, "
                "exact Newton--Hensel lifting mod 25/125/625, denominator-monomial unit "
                "theorem, walk-monomial and matrix-power holdout residues"
            ),
            "no_float_decisions": True,
            "script": "experiments/e146_kw_boundary.py",
            "prior_artifact": "results/kac_ward/components.json",
            "log_p_coefficients_3x3x3": log_p_component,
            "process_seconds": round(time.process_time() - started, 3),
        },
    }
    if not all(check["passed"] for check in payload["checks"]):
        raise AssertionError("one boundary-lift check failed")
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS census 1553 (unit 3), tau histogram {dict(sorted(tau_histogram.items()))}")
    print(
        "PASS max-rank stratum: 4 points; 1716-chart census per point; "
        f"nonsingular charts {[r['nonsingular_charts_of_1716_primitive_rows'] for r in chart_records]}"
    )
    for record in lift_records:
        if record["unit_lift"]:
            print(
                f"PASS unit lift {record['point13_mod5'][:6]} -> {record['lift_solve6_mod625']} "
                f"holdout {record['holdout']['walk_residues_mod625']}"
            )
        else:
            print(
                f"PASS boundary lift {record['point13_mod5'][:6]} -> {record['lift_solve6_mod625']} "
                f"chart-degenerate (denominators {record['vanishing_thin_denominators']}), "
                f"holdout {record['holdout']['status']}"
            )
    print(f"PASS boundary chart exclusion: {boundary_with_vanishing}/1550 vanishing denominators")
    print(f"PASS event seconds: {round(time.process_time() - started, 1)}")


if __name__ == "__main__":
    main()
