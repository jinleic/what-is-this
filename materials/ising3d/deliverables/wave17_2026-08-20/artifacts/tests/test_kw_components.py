"""Clean-room exact verifier for the wave-14 Kac--Ward branch-``11111`` component theorem.

This file imports no producer experiment.  It reconstructs the finite
Kac--Ward construction system from stable ``ising`` APIs, recomputes the three
mandatory surfaces from raw inputs:

* the F_5 solution count of the anchor thin-torus section ``(1,1,1,1,1,4,4)``
  (exactly three, all nonsingular, anchor Jacobian determinant 2 mod 5),
* the ``250 mod 625`` order-eight holdout reproduction at the anchor lift
  (the wave-10 value stored in ``results/kac_ward/branch11111_cubic.json``),
* the order-10/12 holdout residues at all three lifts by an independent transfer-
  matrix walk expansion, cross-validated against a matrix-power trace route,

and the orbit-union lemma (kernel dimension, unimodular diagonal minor,
H_4 = 0 identically, characters ``(4,-2,6)`` and ``(5,-1,7)``, and the
per-orbit zero-pattern uniformity at the anchor).  All arithmetic is exact integer
modular arithmetic over ``(Z/625Z)``; no floating point decides anything.
"""

from __future__ import annotations

import itertools
import json
from fractions import Fraction
from math import gcd
from pathlib import Path

import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    ALLOWED_DIRECTION_PAIRS,
    DIRECTION_GAUGE_TREE,
    _transition_indices,
    directed_edges,
    formal_log_coefficients,
    full_weight_finite_system,
)
from ising.lattices import cubic


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/components.json"
PRIOR = ROOT / "results/kac_ward/branch11111_cubic.json"
BASE_ORDERS = (4, 6, 8, 10, 12)
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
    "u_py_px", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx", "u_pz_px",
    "u_pz_mx", "u_pz_py", "u_pz_my", "u_mz_px", "u_mz_mx", "u_mz_py",
    "u_mz_my",
]
FREE_NAMES = [
    "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px",
    "u_my_mx", "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz",
    "u_mz_px", "u_mz_mx", "u_mz_py", "u_mz_my",
]
ANCHOR_HELD = [1, 1, 1, 1, 1, 4, 4]
PRIME = 5
MODULUS = 625
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
        symbol["u_my_mz"]: -1
        / (symbol["u_mz_py"] * symbol["u_py_pz"] * symbol["u_pz_my"]),
        symbol["u_my_pz"]: -1
        / (symbol["u_mz_my"] * symbol["u_py_mz"] * symbol["u_pz_py"]),
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


def integer_vector(entries):
    denominator = 1
    for entry in entries:
        denominator = sp.ilcm(denominator, entry.q)
    integers = [int(entry * denominator) for entry in entries]
    divisor = 0
    for entry in integers:
        divisor = gcd(divisor, abs(entry))
    integers = [entry // divisor for entry in integers]
    if integers and next((entry for entry in integers if entry), 0) < 0:
        integers = [-entry for entry in integers]
    return integers


def eval_poly_values(polynomial, values, modulus):
    total = 0
    for exponents, coefficient in polynomial.terms():
        term = int(coefficient) % modulus
        for value, exponent in zip(values, exponents, strict=True):
            if exponent:
                term = term * pow(int(value), int(exponent), modulus) % modulus
        total = (total + term) % modulus
    return total


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


def hensel_lift_slice(slice_polynomials, solve_variables, initial, prime, exponent):
    """Hensel lift of a nonsingular F_prime zero to mod prime**exponent (exact integers)."""
    jacobian = [
        [
            eval_poly_values(
                sp.Poly(sp.diff(polynomial.as_expr(), variable), *solve_variables, domain=sp.ZZ),
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
        residuals = []
        for polynomial in slice_polynomials:
            evaluated = 0
            for exponents, coefficient in polynomial.terms():
                term = int(coefficient)
                for value, entry in zip(values, exponents, strict=True):
                    if entry:
                        term *= pow(value, entry)
                evaluated += term
            residuals.append(evaluated)
        assert all(residual % modulus == 0 for residual in residuals), "Hensel invariant failed"
        correction = solve_linear_system(
            jacobian, [-(residual // modulus) % prime for residual in residuals], prime
        )
        if correction is None:
            return None, jacobian
        values = [value + modulus * delta for value, delta in zip(values, correction, strict=True)]
        modulus *= prime
    return values, jacobian


def symbolic_trace_monomials(shape, max_length):
    """Exact Tr(T^length) monomial polys by transfer-walk expansion (no matrix squaring)."""
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
        traces[length] = {}
        for start in range(n):
            for monomial, count in rows[start].get(start, {}).items():
                traces[length][monomial] = traces[length].get(monomial, 0) + count
    return traces


def matrix_trace_power(nongauge_values, structure, length, modulus):
    """Matrix-power route: Tr(T^length) mod modulus by binary exponentiation."""
    gauge = set(DIRECTION_GAUGE_TREE)
    n = structure["n"]
    transitions = structure["transitions"]
    directions = structure["directions"]
    index = {pair: position for position, pair in enumerate([
        pair for pair in ALLOWED_DIRECTION_PAIRS if pair not in gauge
    ])}
    matrix = [[0] * n for _ in range(n)]
    for edge in range(n):
        for next_edge in transitions[edge]:
            key = (directions[edge], directions[next_edge])
            value = 1 if key in gauge else nongauge_values[index[key]]
            matrix[edge][next_edge] = value
    result = [[int(row == column) for column in range(n)] for row in range(n)]
    base = matrix
    while length:
        if length & 1:
            result = mul_mod(result, base, n, modulus)
        base = mul_mod(base, base, n, modulus)
        length >>= 1
    return sum(result[edge][edge] for edge in range(n)) % modulus


def mul_mod(left, right, n, modulus):
    output = [[0] * n for _ in range(n)]
    for row in range(n):
        for middle in range(n):
            if not left[row][middle]:
                continue
            factor = left[row][middle]
            for column in range(n):
                output[row][column] = (output[row][column] + factor * right[middle][column]) % modulus
    return output


def holdout_residues_at_lift(
    traces, walk_nongauge_values, structure, nongauge, log_p, modulus
):
    """Symbolic residues at orders 4,6,8,10,12 by the walk expansion."""
    residues = []
    for order in BASE_ORDERS:
        total = 0
        for monomial, count in traces[order].items():
            term = count % modulus
            for position, exponent in enumerate(monomial):
                if exponent:
                    term = term * pow(int(walk_nongauge_values[position]), exponent, modulus) % modulus
            total = (total + term) % modulus
        required = -2 * order * log_p[order]
        residues.append(
            (total - int(required.numerator) * pow(int(required.denominator) % modulus, -1, modulus))
            % modulus
        )
    return residues


def main() -> None:
    # ---------------- raw build ----------------
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, (4, 6, 8), gauge_fix=False)
    variables, equations = branch_11111_equations(system)
    assert len(equations) == 42
    symbol = {str(value): value for value in variables}
    pair_sym = dict(system.all_symbols)
    free = [symbol[name] for name in FREE_NAMES]
    non_diagonal = [symbol[name] for name in NONDIAGONAL_NAMES]
    diagonal = {symbol[name]: 1 for name in DIAGONAL_NAMES}
    thin = thin_substitutions(symbol)
    sparse42 = [sparse_poly(expression, variables) for expression in equations]

    primitives = []
    for index in SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(equations[index].subs(thin))
        numerator, _ = sp.fraction(reduced)
        _content, primitive = sp.Poly(
            sp.expand(numerator).subs(diagonal), *non_diagonal, domain=sp.ZZ
        ).primitive()
        primitives.append(primitive)
    assert [len(polynomial.terms()) for polynomial in primitives] == [33, 277, 155, 7, 155, 155]

    # ---------------- control 1: the three F_5 solutions at the anchor section ----
    solutions = []
    for candidate in itertools.product((1, 2, 3, 4), repeat=6):
        point = list(candidate) + list(ANCHOR_HELD)
        if all(eval_poly_values(polynomial, point, PRIME) == 0 for polynomial in primitives):
            solutions.append(point)
    assert len(solutions) == 3, f"expected 3 anchor F5 solutions, got {len(solutions)}"
    assert solutions[0] == [1, 2, 3, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4]

    # ---------------- control 2: Jacobian determinant 2 mod 5 at the anchor ----
    slice_polynomials = []
    for polynomial in primitives:
        held = {non_diagonal[6 + index]: value for index, value in enumerate(ANCHOR_HELD)}
        slice_polynomials.append(
            sp.Poly(polynomial.as_expr().subs(held), *non_diagonal[:6], domain=sp.ZZ)
        )
    jacobian = [
        [
            eval_poly_values(
                sp.Poly(sp.diff(polynomial.as_expr(), variable), *non_diagonal[:6], domain=sp.ZZ),
                solutions[0][:6],
                PRIME,
            )
            for variable in non_diagonal[:6]
        ]
        for polynomial in slice_polynomials
    ]
    anchor_det = int(sp.Matrix(jacobian).det()) % PRIME
    assert anchor_det == 2, f"anchor Jacobian determinant {anchor_det} != 2 mod 5"

    # ---------------- Hensel lifts + construction + holdout residues ---------------
    dets = []
    for point in solutions:
        held = {non_diagonal[6 + index]: value for index, value in enumerate(ANCHOR_HELD)}
        slice_polynomials = [
            sp.Poly(polynomial.as_expr().subs(held), *non_diagonal[:6], domain=sp.ZZ)
            for polynomial in primitives
        ]
        jacobian = [
            [
                eval_poly_values(
                    sp.Poly(sp.diff(polynomial.as_expr(), variable), *non_diagonal[:6], domain=sp.ZZ),
                    point[:6],
                    PRIME,
                )
                for variable in non_diagonal[:6]
            ]
            for polynomial in slice_polynomials
        ]
        dets.append(int(sp.Matrix(jacobian).det()) % PRIME)
        assert dets[-1] != 0
        lift, _ = hensel_lift_slice(slice_polynomials, non_diagonal[:6], point[:6], PRIME, 4)
        assert lift is not None
        lift13 = lift + point[6:]

        # 25-coordinate reconstruction via the thin substitution chart (diagonals 1)
        values = {name: value for name, value in zip(NONDIAGONAL_NAMES, lift13)}
        values.update({name: 1 for name in DIAGONAL_NAMES})
        free_values = [values[name] for name in FREE_NAMES]
        for name, expression in thin.items():
            numerator, denominator = sp.fraction(sp.cancel(expression))
            poly_num = sparse_poly(sp.Poly(numerator, *free, domain=sp.ZZ), free)
            poly_den = sparse_poly(sp.Poly(denominator, *free, domain=sp.ZZ), free)
            value_num = eval_sparse(poly_num, free_values, MODULUS)
            value_den = eval_sparse(poly_den, free_values, MODULUS)
            assert value_den % PRIME != 0
            values[str(name)] = value_num * pow(value_den, -1, MODULUS) % MODULUS
        values25 = [values[str(value)] for value in variables]
        construction = [eval_sparse(polynomial, values25, MODULUS) for polynomial in sparse42]
        assert construction == [0] * 42

        # control 3 / cross-validation: order-8 residue 250 at the anchor via the
        # full_weight 3x3x3 equations, then the walk + matrix residues at 4..12
        walk_values = [values[str(pair_sym[pair])] for pair in ALLOWED_DIRECTION_PAIRS if pair not in set(DIRECTION_GAUGE_TREE)]
        traces = symbolic_trace_monomials(HOLDOUT_SHAPE, 12)
        polynomial = even_subgraph_polynomial(cubic(*HOLDOUT_SHAPE, periodic=False))
        log_p = formal_log_coefficients(polynomial, 12)
        lattice = cubic(*HOLDOUT_SHAPE, periodic=False)
        edges = directed_edges(lattice)
        structure = {
            "n": len(edges),
            "transitions": _transition_indices(lattice),
            "directions": [
                DIRECTIONS_3D.index(tuple(edge.direction)) for edge in edges
            ],
        }
        walk_residues = []
        for order in BASE_ORDERS:
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
        matrix_residues = []
        for order in BASE_ORDERS:
            required = -2 * order * log_p[order]
            trace = matrix_trace_power(walk_values, structure, order, MODULUS)
            matrix_residues.append(
                (trace - int(required.numerator) * pow(int(required.denominator) % MODULUS, -1, MODULUS))
                % MODULUS
            )
        assert walk_residues == matrix_residues
        recorded = {
            "point13_mod5": point,
            "lift_solve6_mod625": lift,
            "residues": walk_residues,
        }
        if point == solutions[0]:
            assert anchor_det == 2
            assert recorded["residues"] == [0, 0, 250, 55, 266]
            assert walk_residues[2] == 250
        print(
            f"recomputed lift {point[:6]} det {dets[-1]} residues 4,6,8,10,12 "
            f"{walk_residues} (walk == matrix)"
        )

    # ---------------- orbit-union lemma ----------------
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
    kernel = [integer_vector(vector) for vector in sp.Matrix(support_rows).nullspace()]
    assert len(kernel) == 3
    diagonal_positions = [index for index, name in enumerate(FREE_NAMES) if name in DIAGONAL_NAMES]
    assert abs(int(sp.Matrix(kernel)[:, diagonal_positions].det())) == 1

    holdout_system = full_weight_finite_system((HOLDOUT_SHAPE,), (4, 6, 8), gauge_fix=False)
    holdout_variables, holdout_equations = branch_11111_equations(holdout_system)
    holdout16 = []
    for expression in holdout_equations:
        reduced = sp.cancel(expression.subs(thin))
        numerator, _ = sp.fraction(reduced)
        holdout16.append(sp.Poly(sp.expand(numerator), *free, domain=sp.ZZ))
    for polynomial in num16:
        for vector in kernel:
            dots = {
                sum(exponent * weight for exponent, weight in zip(monomial, vector, strict=True))
                for monomial, _coefficient in polynomial.terms()
            }
            assert len(dots) == 1
    assert holdout16[0].terms() in ([], [((0,) * 16, 0)])
    characters = {}
    for label, polynomial in (("H6", holdout16[1]), ("H8", holdout16[2])):
        dots = []
        for vector in kernel:
            values = {
                sum(exponent * weight for exponent, weight in zip(monomial, vector, strict=True))
                for monomial, _coefficient in polynomial.terms()
            }
            assert len(values) == 1
            dots.append(next(iter(values)))
        characters[label] = dots
    assert characters == {"H6": [4, -2, 6], "H8": [5, -1, 7]}

    # empirical zero-pattern uniformity at the anchor section
    for point in solutions:
        point_map = dict(zip(NONDIAGONAL_NAMES, point, strict=True))
        patterns = set()
        for t1, t2, t3 in itertools.product((1, 2, 3, 4), repeat=3):
            translate = [0] * 16
            for position in range(16):
                value = point_map.get(FREE_NAMES[position], 1)
                for translate_value, vector in zip((t1, t2, t3), kernel):
                    if vector[position]:
                        value = value * pow(translate_value, vector[position] % 4, 5) % 5
                translate[position] = value
            residues = []
            for polynomial in (holdout16[1], holdout16[2]):
                total = 0
                for monomial, coefficient in polynomial.terms():
                    term = coefficient % 5
                    for position, exponent in enumerate(monomial):
                        if exponent:
                            term = term * pow(translate[position], exponent, 5) % 5
                    total = (total + term) % 5
                residues.append(total == 0)
            patterns.add(tuple(residues))
        assert len(patterns) == 1

    # ---------------- artifact agreement ----------------
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["status"] == "UNRESOLVED"
    assert payload["hensel_components"]["three_five_solution_count"] == 3
    assert payload["hensel_components"]["anchor_jacobian_determinant_mod5"] == 2
    assert payload["holdout"]["first_holdout_k8_reproduction_mod625"] == 250
    assert payload["holdout"]["k8_residues"] == [250, 456, 430]
    assert payload["holdout"]["k10_residues"] == [55, 350, 480]
    assert payload["holdout"]["k12_residues"] == [266, 132, 286]
    assert payload["orbit_union_lemma"]["characters"] == {
        "H6": [4, -2, 6],
        "H8": [5, -1, 7],
    }
    assert payload["orbit_union_lemma"]["diagonal_minor_abs_det"] == 1
    assert payload["orbit_union_lemma"]["h4_identically_zero_on_thin_torus"]
    assert all(check["passed"] for check in payload["checks"])
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    assert (
        prior["data"]["hensel_five_adic_section"]["holdout_residues_mod_625"] == [0, 0, 250]
    )
    assert (
        prior["data"]["hensel_five_adic_section"]["jacobian_determinant_mod_5"] == 2
    )
    assert (
        prior["data"]["hensel_five_adic_section"]["lifted_solve_values_mod_625"]
        == [226, 197, 98, 63, 161, 216]
    )
    print("three F5 solutions, dets", dets, "cross-validated residues and orbit-union lemma agree")
    print("PASS")


if __name__ == "__main__":
    main()
