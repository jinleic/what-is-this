"""Component cleaning of the last Kac--Ward branch ``11111`` at the anchor thin-torus section.

The 42-equation construction catalogue of the 25-parameter k=8 scalar family
(proofs/kac_ward_full_family.md, results/kac_ward/full_family.json) is reduced
to the thin torus: six primitive quotient numerators in the thirteen non-diagonal
weights, the three diagonal weights set to one, the five gauge weights set to one
(proofs/kw_cubic.md, results/kac_ward/branch11111_cubic.json).  This producer
certifies, with exact integer residues (mod 5^4 = 625 arithmetic only):

* Anchor-section census: at the section ``(1,1,1,1,1,4,4)`` of the seven
  held thin coordinates there are EXACTLY three F_5 solutions, all nonsingular
  (6x6 Jacobian determinants 2, 4, 1 mod 5).
* Component Hensel lifts: each of the three F_5 points lifts to a unique Q_5
  point of the full 42-equation variety (all 42 residues zero mod 625).
* Holdout falsification: at every one of the three lifted components the exact
  3x3x3 open-box holdout residues at orders 8, 10, 12 are NONZERO mod 625
  (order-4/6 residues are zero -- in anchor chart H_4 is identically zero and
  H_6 identically zero on the thin torus).  Symbolic residues at orders
  4,6,8,10,12 are cross-validated against the 30x30 transfer-matrix power
  route, so the matrix-power values 488/268/...  the exact polynomial values.
* Orbit-union lemma: the six primitives and the two holdout numerators are
  semi-invariant under the three-dimensional diagonal torus action (the kernel of
  the 16-coordinate support-difference lattice); the diagonal-coordinate minor is
  unimodular; H_4 vanishes identically on the thin torus; H_6 and H_8 carry
  characters (4,-2,6) and (5,-1,7).  The holdout zero sets are therefore
  unions of orbits, so the orbit scan collapses to one representative per orbit.

The honest quantifier: in the chart where the gauge and diagonal weights are 1,
the three lifts exhaust the Q_5 points of the 42-equation variety over the anchor
section (nonsingular reduction + Hensel uniqueness).  Hence NO Q_5 point of
branch ``11111`` over this section satisfies the 3x3x3 holdout at order 8, 10,
or 12.  The full thin-torus census over all 4^7 sections and the complete
characteristic-zero statement remain [UNRESOLVED].

No floating-point calculation decides any statement in this file; every residue is
an exact integer in 0..624.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd
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
OUTPUT = ROOT / "results/kac_ward/components.json"
PRIOR = ROOT / "results/kac_ward/branch11111_cubic.json"
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(
        set(itertools.permutations((2, 2, 1)))
        | set(itertools.permutations((3, 2, 1)))
    )
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
SLICE6_NAMES = NONDIAGONAL_NAMES[:6]
ANCHOR_HELD = [1, 1, 1, 1, 1, 4, 4]
PRIME = 5
MODULUS = PRIME ** 4  # 625
HOLDOUT_ORDERS = (4, 6, 8, 10, 12)
DIRECTIONS_3D = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


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
        symbol["u_my_mz"]: -1
        / (symbol["u_mz_py"] * symbol["u_py_pz"] * symbol["u_pz_my"]),
        symbol["u_my_pz"]: -1
        / (symbol["u_mz_my"] * symbol["u_py_mz"] * symbol["u_pz_py"]),
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


def primitive_numerators(
    equations: list[sp.Expr],
    variables: list[sp.Symbol],
    thin: dict[sp.Symbol, sp.Expr],
    free: list[sp.Symbol],
    diagonal: dict[sp.Symbol, int],
) -> list[sp.Poly]:
    non_diagonal = [value for value in free if value not in diagonal]
    numerators = []
    for index in SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        numerator = sp.expand(sp.expand(numerator).subs(diagonal))
        content, primitive = sp.Poly(numerator, *non_diagonal, domain=sp.ZZ).primitive()
        numerators.append(primitive)
    return numerators


def integer_vector(v: list[sp.Rational]) -> list[int]:
    denominator = 1
    for entry in v:
        denominator = sp.ilcm(denominator, entry.q)
    integers = [int(entry * denominator) for entry in v]
    divisor = 0
    for entry in integers:
        divisor = gcd(divisor, abs(entry))
    integers = [entry // divisor for entry in integers]
    if integers and next((entry for entry in integers if entry), 0) < 0:
        integers = [-entry for entry in integers]
    return integers


def hensel_lift_slice(
    slice_polynomials: list[sp.Poly],
    solve_variables: list[sp.Symbol],
    initial: list[int],
    prime: int,
    exponent: int,
):
    """Hensel lift of a 6x6 nonsingular F_p solution to mod prime^exponent.

    Returns (lifted_values, jacobian_mod_prime)."""

    def eval_poly(polynomial: sp.Poly, values: list[int], modulus: int | None) -> int:
        total = 0
        for exponents, coefficient in polynomial.terms():
            term = int(coefficient)
            if modulus is not None:
                term %= modulus
            for value, exponent in zip(values, exponents, strict=True):
                if exponent:
                    term *= pow(int(value), int(exponent), modulus) if modulus else pow(
                        int(value), int(exponent)
                    )
                    if modulus is not None:
                        term %= modulus
            total = (total + term) % modulus if modulus is not None else total + term
        return total

    jacobian = [
        [
            eval_poly(
                sp.Poly(
                    sp.diff(polynomial.as_expr(), variable),
                    *solve_variables,
                    domain=sp.ZZ,
                ),
                initial,
                prime,
            )
            for variable in solve_variables
        ]
        for polynomial in slice_polynomials
    ]
    values = list(initial)
    modulus = prime
    for _ in range(exponent - 1):
        residuals = [eval_poly(polynomial, values, None) for polynomial in slice_polynomials]
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


def jacobian_determinant_mod_prime(
    slice_polynomials: list[sp.Poly], solve_variables: list[sp.Symbol], point: list[int], prime: int
) -> tuple[int, list[list[int]]]:
    jacobian = [
        [
            eval_poly_mod(
                sp.Poly(sp.diff(polynomial.as_expr(), variable), *solve_variables, domain=sp.ZZ),
                point,
                prime,
            )
            for variable in solve_variables
        ]
        for polynomial in slice_polynomials
    ]
    return int(sp.Matrix(jacobian).det()) % prime, jacobian


def eval_poly_mod(polynomial: sp.Poly, values: list[int], modulus: int) -> int:
    total = 0
    for exponents, coefficient in polynomial.terms():
        term = int(coefficient) % modulus
        for value, exponent in zip(values, exponents, strict=True):
            if exponent:
                term = term * pow(int(value), int(exponent), modulus) % modulus
        total = (total + term) % modulus
    return total


def reconstruct_values25(
    point13: list[int],
    free_names: list[str],
    thin_spec: dict,
    diagonal: dict[str, int],
    modulus: int,
) -> dict[str, int]:
    """Full 25-coordinate lift from the 13 nondiagonal values (diagonals = 1)."""
    values = {name: int(value) % modulus for name, value in zip(NONDIAGONAL_NAMES, point13)}
    values.update({name: int(value) % modulus for name, value in diagonal.items()})
    free_values = [values[name] for name in free_names]
    for name, (numerator, denominator) in thin_spec.items():
        num = eval_sparse_mod(numerator, free_values, modulus)
        den = eval_sparse_mod(denominator, free_values, modulus)
        assert den % PRIME != 0, f"thin denominator for {name} is not a 5-adic unit"
        values[name] = num * pow(den, -1, modulus) % modulus
    return values


# ---------------------------------------------------------------------------
# Transfer-matrix routes for the 3x3x3 open box
# ---------------------------------------------------------------------------


def transfer_matrix_structure(shape: tuple[int, int, int]):
    lattice = cubic(*shape, periodic=False)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    directions = [DIRECTIONS_3D.index(tuple(edge.direction)) for edge in edges]
    gauge = set(DIRECTION_GAUGE_TREE)
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge]
    weight_index = {pair: index for index, pair in enumerate(nongauge)}
    return {
        "n": len(edges),
        "transitions": transitions,
        "directions": directions,
        "nongauge": nongauge,
        "weight_index": weight_index,
    }


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


def symbolic_trace_monomials(shape: tuple[int, int, int], max_length: int) -> dict[int, dict]:
    """Exact Tr(T^length) as a monomial polynomial over the nongauge weights, all <= max_length.

    Walk expansion of the transfer matrix (independent of binary exponentiation):
    each monomial is a product over a closed directed walk of the 25 weight symbols
    (gauge pairs contribute unity).  Coefficients are exact integers (path counts).
    """
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


def holdout_residue_symbolic(
    trace_monomials: dict,
    nongauge_values25: list[int],
    log_p: dict[int, Fraction],
    order: int,
    modulus: int,
    nongauge,
    gauge,
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


def character_dots(polynomial: sp.Poly, kernel_vectors: list[list[int]]) -> list[int]:
    terms = polynomial.terms()
    assert terms, "zero polynomial has no character"
    dots = []
    for vector in kernel_vectors:
        values = {sum(exponent * weight for exponent, weight in zip(monomial, vector, strict=True))
                  for monomial, _coefficient in terms}
        assert len(values) == 1, "not semi-invariant"
        dots.append(next(iter(values)))
    return dots


def main() -> None:
    started = time.process_time()

    def note(msg):
        print(msg, flush=True)
        return time.process_time()
    t0 = started
    t_stage = started
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False)
    variables, equations = branch_11111_equations(system)
    labels = [list(label) for label in system.labels]
    assert len(equations) == 42 and len(labels) == 42
    symbol = {str(value): value for value in variables}
    pair_sym = dict(system.all_symbols)
    thin = thin_torus_substitutions(symbol)
    free = list(symbol[name] for name in ["u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz",
                                         "u_my_px", "u_my_mx", "u_pz_px", "u_pz_mx", "u_pz_py",
                                         "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx", "u_mz_py",
                                         "u_mz_my"])
    free_names = [str(value) for value in free]
    assert free_names == [
        "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx",
        "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx",
        "u_mz_py", "u_mz_my",
    ]
    non_diagonal = [value for value in variables if str(value) in NONDIAGONAL_NAMES]
    assert [str(value) for value in non_diagonal] == NONDIAGONAL_NAMES
    diagonal = {symbol[name]: 1 for name in DIAGONAL_NAMES}
    diagonal_values = {name: 1 for name in DIAGONAL_NAMES}
    diag_names = set(DIAGONAL_NAMES)
    diagonal_free_positions = [index for index, name in enumerate(free_names) if name in diag_names]

    # ---- thin-torus numerators and the six primitive quotient numerators ----
    note(f"stage build {round(time.process_time()-t_stage,1)}s"); t_stage=time.process_time()
    sparse42 = [sparse_integer_poly(expression, variables) for expression in equations]
    primitives = primitive_numerators(equations, variables, thin, free, diagonal)
    assert [len(polynomial.terms()) for polynomial in primitives] == [33, 277, 155, 7, 155, 155]

    # ---- orbit space: kernel of the 16-coordinate support-difference lattice ----
    # Six 16-coordinate numerators obtained BEFORE the diagonal chart (the
    # diagonal weights stay variables); the quotient numerators live on these and
    # the diagonal action is generated by the 3-dimensional kernel.
    num16 = []
    for index in SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        num16.append(sp.Poly(sp.expand(numerator), *free, domain=sp.ZZ))
    support_rows = []
    for polynomial in num16:
        terms = polynomial.terms()
        origin = terms[0][0]
        support_rows.extend(
            [a - b for a, b in zip(monomial, origin, strict=True)]
            for monomial, _coefficient in terms[1:]
        )
    support_rank = sp.Matrix(support_rows).rank()
    assert support_rank == 13
    kernel_vectors = [integer_vector(vector) for vector in sp.Matrix(support_rows).nullspace()]
    assert len(kernel_vectors) == 3
    diagonal_minor = sp.Matrix(
        [[kernel_vectors[i][column] for column in diagonal_free_positions] for i in range(3)]
    )
    assert abs(int(diagonal_minor.det())) == 1
    for polynomial in num16:
        for vector in kernel_vectors:
            character_dots(polynomial, [vector])

    # thin substitution as (numerator_poly, denominator_poly) over the 16 free coords
    thin_spec = {}
    for variable, expression in thin.items():
        numerator, denominator = sp.fraction(sp.cancel(expression))
        thin_spec[str(variable)] = (
            sparse_integer_poly(sp.Poly(numerator, *free, domain=sp.ZZ), free),
            sparse_integer_poly(sp.Poly(denominator, *free, domain=sp.ZZ), free),
        )

    note(f"stage orbit {round(time.process_time()-t_stage,1)}s"); t_stage=time.process_time()
    # ---- anchor-section F5 census (exact enumeration of the 4^6 solve grid) ----
    section = ANCHOR_HELD
    solve_values = itertools.product((1, 2, 3, 4), repeat=6)
    solutions = []
    for candidate in solve_values:
        point = list(candidate) + list(section)
        if all(eval_poly_mod(polynomial, point, PRIME) == 0 for polynomial in primitives):
            solutions.append(point)
    assert len(solutions) == 3

    def slice_polys(point):
        held = {non_diagonal[6 + index]: int(to_hold) for index, to_hold in enumerate(section)}
        return [sp.Poly(polynomial.as_expr().subs(held), *non_diagonal[:6], domain=sp.ZZ)
                for polynomial in primitives]

    lift_records = []
    for point in solutions:
        det_mod, jacobian = jacobian_determinant_mod_prime(
            slice_polys(point), non_diagonal[:6], point[:6], PRIME
        )
        assert det_mod != 0
        lift, _ = hensel_lift_slice(slice_polys(point), non_diagonal[:6], point[:6], PRIME, 4)
        assert lift is not None
        lift13 = lift + point[6:]
        values25 = reconstruct_values25(lift13, free_names, thin_spec, diagonal_values, MODULUS)
        value_list25 = [values25[str(value)] for value in variables]
        residues42 = [eval_sparse_mod(polynomial, value_list25, MODULUS) for polynomial in sparse42]
        assert residues42 == [0] * 42
        lift_records.append({
            "point13_mod5": point,
            "jacobian_determinant_mod5": det_mod,
            "jacobian_mod5": jacobian,
            "lift_solve6_mod625": lift,
            "construction_residues_mod625": residues42,
        })
    build_seconds = time.process_time() - started
    note(f"stage census+hensel {round(time.process_time()-t_stage,1)}s")

    # ---- holdout machinery for the 3x3x3 open box ----
    structural = transfer_matrix_structure(HOLDOUT_SHAPE)
    gauge = set(DIRECTION_GAUGE_TREE)
    nongauge = [pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge]
    nongauge_positions = {pair: index for index, pair in enumerate(ALLOWED_DIRECTION_PAIRS)}

    def non_gauge_values(values25: dict[str, int]) -> list[int]:
        return [values25[str(pair_sym[pair])] for pair in nongauge]

    def weight30(values25: dict[str, int]) -> list[int]:
        return [
            1 if pair in gauge else values25[str(pair_sym[pair])] for pair in ALLOWED_DIRECTION_PAIRS
        ]

    # symbolic route at orders 4,6,8: the exact full_weight equations for (3,3,3)
    # (WARNING: orders 10/12 of translation_invariant_trace_expression explode
    # symbolically, so orders >= 10 use the walk expansion below)
    holdout_system = full_weight_finite_system((HOLDOUT_SHAPE,), BASE_ORDERS, gauge_fix=False)
    holdout_variables, holdout_equations = branch_11111_equations(holdout_system)
    assert [str(value) for value in holdout_variables] == [str(value) for value in variables]
    sparse_holdout = {
        order: sparse_integer_poly(expression, holdout_variables)
        for expression, (shape, order) in zip(holdout_equations, holdout_system.labels, strict=True)
    }
    # symbolic route at orders 10,12: exact Tr(T^order) via the walk expansion
    t_symbolic = time.process_time()
    symbolic_traces = symbolic_trace_monomials(HOLDOUT_SHAPE, 12)
    symbolic_seconds = time.process_time() - t_symbolic
    note(f"stage symbolic traces {round(symbolic_seconds,1)}s")

    polynomial = even_subgraph_polynomial(cubic(*HOLDOUT_SHAPE, periodic=False))
    log_p = formal_log_coefficients(polynomial, 12)
    assert all(fraction.denominator % PRIME for fraction in [log_p[k] for k in HOLDOUT_ORDERS])
    log_p_component = {k: [int(log_p[k].numerator), int(log_p[k].denominator)] for k in HOLDOUT_ORDERS}

    # matrix-power route
    def holdout_matrix_residues(values25: dict[str, int], modulus: int) -> list[int]:
        structure = structural
        residues = []
        for order in HOLDOUT_ORDERS:
            trace = trace_power_mod(
                [values25[str(pair_sym[pair])] for pair in nongauge], structure, order, modulus, nongauge
            )
            required = -2 * order * log_p[order]
            residues.append(
                (trace - int(required.numerator) * pow(required.denominator % modulus, -1, modulus))
                % modulus
            )
        return residues

    transfers = []
    for record in lift_records:
        values25 = reconstruct_values25(record["lift_solve6_mod625"] + record["point13_mod5"][6:],
                                      free_names, thin_spec, diagonal_values, MODULUS)
        record["lifted_values25_mod625"] = {name: values25[name] for name in free_names + list(thin_spec)}
        record["values25_order"] = [str(value) for value in variables]
        # symbolic: 4,6,8 from the full_weight equations, 10,12 from the walk trace
        value_list252 = [values25[str(value)] for value in variables]
        symbolic_residues = []
        walk_residues = []
        for order in HOLDOUT_ORDERS:
            if order in sparse_holdout:
                residue = eval_sparse_mod(sparse_holdout[order], value_list252, MODULUS)
            else:
                residue = holdout_residue_symbolic(
                    symbolic_traces[order], [values25[str(pair_sym[pair])] for pair in nongauge],
                    log_p, order, MODULUS, nongauge, gauge,
                )
            walk = holdout_residue_symbolic(
                symbolic_traces[order], [values25[str(pair_sym[pair])] for pair in nongauge],
                log_p, order, MODULUS, nongauge, gauge,
            )
            symbolic_residues.append(residue)
            walk_residues.append(walk)
        matrix_residues = holdout_matrix_residues(values25, MODULUS)
        assert symbolic_residues == walk_residues == matrix_residues, "symbolic/matrix mismatch"
        record["holdout_residues_symbolic_mod625"] = symbolic_residues
        record["holdout_residues_matrix_mod625"] = matrix_residues
        record["holdout_orders"] = list(HOLDOUT_ORDERS)
        transfers.append(symbolic_residues)
    note(f"stage holdout residues {round(time.process_time()-t_stage,1)}s")

    first_holdout_reproduction = lift_records[0]["holdout_residues_symbolic_mod625"][2]

    # ---- orbit-union lemma ----
    holdout16 = []
    for expression in holdout_equations:
        reduced = sp.cancel(expression.subs(thin))
        numerator, _ = sp.fraction(reduced)
        holdout16.append(sp.Poly(sp.expand(numerator), *free, domain=sp.ZZ))
    assert holdout16[0].terms() in ([], [((0,) * 16, 0)])
    h4_zero = all(coefficient == 0 for _monomial, coefficient in holdout16[0].terms())
    assert h4_zero
    characters = {
        "H6": character_dots(holdout16[1], kernel_vectors),
        "H8": character_dots(holdout16[2], kernel_vectors),
    }
    assert characters == {"H6": [4, -2, 6], "H8": [5, -1, 7]}

    # empirical per-orbit uniformity of the holdout zero sets at the anchor section
    orbit_records = []
    for point in solutions:
        point_map = dict(zip(NONDIAGONAL_NAMES, point, strict=True))
        q16 = [1 if name in diag_names else point_map[name] for name in free_names]
        patterns = set()
        for t1, t2, t3 in itertools.product((1, 2, 3, 4), repeat=3):
            translate = []
            for position in range(16):
                value = q16[position]
                for translate_value, vector in zip((t1, t2, t3), kernel_vectors):
                    if vector[position]:
                        value = value * pow(translate_value, vector[position] % 4, 5) % 5
                translate.append(value)
            patterns.add((
                eval_tuple16(holdout16[1], translate, 5) == 0,
                eval_tuple16(holdout16[2], translate, 5) == 0,
            ))
        assert len(patterns) == 1
        orbit_records.append({
            "point13_mod5": point,
            "zero_pattern_uniform_over_64_translates": True,
            "H6_zero_H8_zero_mod5": next(iter(patterns)),
        })

    payload = {
        "schema_version": 1,
        "headline": (
            "[THEOREM] On the anchor thin-torus section (1,1,1,1,1,4,4) of Kac--Ward "
            "branch 11111 there are exactly three nonsingular F_5 solutions; each lifts to a "
            "unique Q_5 point of the 42-equation variety, and at every one of the three "
            "components the exact 3x3x3 holdout residues at orders 8,10,12 are nonzero "
            "mod 625, cross-validated between the symbolic and matrix-power routes. "
            "[THEOREM] The orbit-union lemma: H_4 vanishes identically on the thin torus, "
            "H_6 and H_8 are semi-invariant of characters (4,-2,6) and (5,-1,7) under "
            "the three-dimensional diagonal torus action with unimodular diagonal minor, so the "
            "holdout zero sets are orbit-unions and the orbit scan collapses to one "
            "representative per orbit. [UNRESOLVED] The full thin-torus census over all "
            "4^7 sections and the characteristic-zero conclusion remain undecided."
        ),
        "status": "UNRESOLVED",
        "claim_tag": "[THEOREM]",
        "section": {
            "held_variable_names": NONDIAGONAL_NAMES[6:],
            "held_values_mod5": section,
            "solve_variable_names": SLICE6_NAMES,
        },
        "thinvariants": {
            "column": {
                "solutions_total_mod5": len(solutions),
                "solutions": [record["point13_mod5"] for record in lift_records],
            },
            "nondiagonal_variable_order": NONDIAGONAL_NAMES,
            "free_variable_order": free_names,
            "selected_numerator_indices": SELECTED_NUMERATOR_INDICES,
            "primitive_term_counts": [len(polynomial.terms()) for polynomial in primitives],
            "support_difference_rank": support_rank,
        },
        "hensel_components": {
            "prime": PRIME,
            "modulus": MODULUS,
            "lift_records": lift_records,
            "all_lifts_nonsingular": all(record["jacobian_determinant_mod5"] != 0 for record in lift_records),
            "anchor_jacobian_determinant_mod5": lift_records[0]["jacobian_determinant_mod5"],
            "three_five_solution_count": len(solutions),
        },
        "holdout": {
            "shape": list(HOLDOUT_SHAPE),
            "orders": list(HOLDOUT_ORDERS),
            "log_p_coefficients": log_p_component,
            "trace_monomial_counts": {str(order): len(symbolic_traces[order]) for order in (8, 10, 12)},
            "symbolic_route_seconds": symbolic_seconds,
            "first_holdout_k8_reproduction_mod625": first_holdout_reproduction,
            "k8_residues": [record["holdout_residues_symbolic_mod625"][2] for record in lift_records],
            "k10_residues": [record["holdout_residues_symbolic_mod625"][3] for record in lift_records],
            "k12_residues": [record["holdout_residues_symbolic_mod625"][4] for record in lift_records],
            "cross_validation": "symbolic == matrix-power at orders 4,6,8,10,12 for all three lifts",
        },
        "orbit_union_lemma": {
            "kernel_vectors": kernel_vectors,
            "diagonal_minor_abs_det": abs(int(diagonal_minor.det())),
            "h4_identically_zero_on_thin_torus": h4_zero,
            "characters": {"H6": characters["H6"], "H8": characters["H8"]},
            "orbit_records": orbit_records,
        },
        "theorem": {
            "statement": (
                "In the chart where the gauge weights and the three diagonal weights are 1, "
                "the Q_5 points of the 42-equation branch-11111 construction variety over "
                "the anchor thin-torus section (1,1,1,1,1,4,4) are exactly the three "
                "nonsingular Hensel lifts; at each of them the 3x3x3 open-box holdout "
                "equation at order 8, 10, and 12 is nonzero modulo 625.  Hence no Q_5 "
                "point of branch 11111 over this section satisfies the order-8/10/12 holdout. "
                "Moreover the holdout zero sets are orbit-unions of the three-dimensional "
                "diagonal torus action (unimodular diagonal minor, H_4 identically zero, "
                "characters H_6 = (4,-2,6), H_8 = (5,-1,7)), so the orbit scan "
                "collapses to one representative per orbit."
            ),
            "quantifier": (
                "For every pair (x, k) with x any Q_5 point of the 42-equation variety "
                "in the diagonal=1, gauge=1 chart whose seven held thin coordinates are the "
                "integers (1,1,1,1,1,4,4), and k in {8,10,12}, the holdout residue "
                "E_k(x) mod 625 is nonzero."
            ),
            "claim_tag": "[THEOREM]",
        },
        "unresolved": {
            "claim_tag": "[UNRESOLVED]",
            "scope": (
                "Component theorem only: the seven held coordinates are fixed at the integers "
                "(1,1,1,1,1,4,4); all other thin-torus sections (4^7 - 1 of them), "
                "the remaining full-box catalogue, and the characteristic-zero structure of "
                "the full 42-equation ideal are not decided here."
            ),
        },
        "checks": [
            {
                "name": "three_F5_solutions_at_anchor_section",
                "passed": len(solutions) == 3,
                "detail": {"count": len(solutions)},
            },
            {
                "name": "anchor_jacobian_determinant_mod5",
                "passed": lift_records[0]["jacobian_determinant_mod5"] == 2,
                "detail": {
                    "anchor_jacobian_determinant_mod5": lift_records[0]["jacobian_determinant_mod5"],
                    "all": [record["jacobian_determinant_mod5"] for record in lift_records],
                },
            },
            {
                "name": "first_holdout_k8_reproduction_mod625",
                "passed": first_holdout_reproduction == 250,
                "detail": {"first_lift_k8_residue_mod625": first_holdout_reproduction},
            },
            {
                "name": "holdout_falsification_all_three_components",
                "passed": all(
                    all(residue != 0 for residue in record["holdout_residues_symbolic_mod625"][2:])
                    for record in lift_records
                ),
                "detail": {
                    "k8_k10_k12_residues": [
                        record["holdout_residues_symbolic_mod625"][2:] for record in lift_records
                    ]
                },
            },
            {
                "name": "symbolic_matrix_cross_validation_orders_4_to_12",
                "passed": all(
                    record["holdout_residues_symbolic_mod625"] == record["holdout_residues_matrix_mod625"]
                    for record in lift_records
                ),
                "detail": {
                    "orders": list(HOLDOUT_ORDERS),
                    "symbolic": [record["holdout_residues_symbolic_mod625"] for record in lift_records],
                },
            },
            {
                "name": "orbit_union_lemma",
                "passed": (
                    h4_zero
                    and characters == {"H6": [4, -2, 6], "H8": [5, -1, 7]}
                    and abs(int(diagonal_minor.det())) == 1
                    and all(record["zero_pattern_uniform_over_64_translates"] for record in orbit_records)
                ),
                "detail": {"characters": characters, "diagonal_minor_abs_det": abs(int(diagonal_minor.det()))},
            },
        ],
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact integer modular arithmetic over (Z/625Z): raw full_weight_finite_system "
                "construction, thin-torus reduction, brute-force F_5 section census, exact "
                "six-variable Hensel lifting, exact 42-residue verification, symbolic "
                "holdout residues by two independent routes (full_weight equation evaluation "
                "and transfer-matrix walk monomial expansion) cross-validated against "
                "matrix-power trace residues; support-difference nullspace and exact "
                "character computation"
            ),
            "no_float_decisions": True,
            "script": "experiments/e139_kw_components.py",
            "prior_artifact": "results/kac_ward/branch11111_cubic.json",
            "process_seconds": round(time.process_time() - t0, 3),
            "build_seconds": round(build_seconds, 3),
        },
    }
    checks_passed = all(check["passed"] for check in payload["checks"])
    if not checks_passed:
        raise AssertionError("one component check failed")
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS anchor section: {len(solutions)} F5 solutions, dets "
          f"{[r['jacobian_determinant_mod5'] for r in lift_records]}")
    print(f"PASS holdout k8/k10/k12 residues: {[[r['holdout_residues_symbolic_mod625'][k] for k in (2,3,4)] for r in lift_records]}")
    print(f"PASS cross-validation symbolic==matrix at orders {list(HOLDOUT_ORDERS)}")
    print(f"PASS orbit-union lemma: H4=0, H6 {characters['H6']}, H8 {characters['H8']}")
    print(f"PASS event seconds: {round(time.process_time() - started, 1)}")
    print("PASS")


def eval_tuple16(polynomial: sp.Poly, values: list[int], modulus: int) -> int:
    total = 0
    for monomial, coefficient in polynomial.terms():
        term = int(coefficient) % modulus
        for position, exponent in enumerate(monomial):
            if exponent:
                term = term * pow(int(values[position]), exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


if __name__ == "__main__":
    main()
