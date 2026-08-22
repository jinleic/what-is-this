"""Exact thin-box closure of the remaining Kac--Ward gauge-tree branches.

This experiment imports the 32-branch map from e50, searches a larger exact
thin-box constraint span through order ten, and records only independently
checkable characteristic-zero certificates.  Modular points and selected
p-adic paths are diagnostics and never decide Q-emptiness.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from e42_kw_hensel import (
    derivative_table,
    eval_system_mod,
    jacobian_mod,
    rank_deficient_first_lift_attempt,
    rref_mod,
)
from e50_kw_branches import (
    branch_equations,
    constant_certificate_search,
    instantiate_linear_solution,
    solve_rectangular_mod,
    sparse_polynomials,
)
from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    formal_log_coefficients,
    full_weight_finite_system,
)
from ising.lattices import cubic

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "results/kac_ward/branches.json"
OUTPUT = ROOT / "results/kac_ward/close.json"

BASE_ORDERS = (4, 6, 8)
SCAN_ORDERS = (4, 6, 8, 10)
SCAN_CANONICAL_THIN_SHAPES = (
    (2, 2, 1),
    (3, 2, 1),
    (4, 2, 1),
    (3, 3, 1),
    (4, 3, 1),
    (4, 4, 1),
)
SCAN_THIN_SHAPES = tuple(
    sorted(
        set().union(
            *(set(itertools.permutations(shape)) for shape in SCAN_CANONICAL_THIN_SHAPES)
        )
    )
)
PRIMARY_THIN_SHAPES = tuple(
    sorted(
        set(itertools.permutations((2, 2, 1)))
        | set(itertools.permutations((3, 2, 1)))
    )
)
MACAULAY_FULL_SHAPES = (
    (3, 3, 2),
    (2, 2, 3),
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
)
MACAULAY_SHAPES = PRIMARY_THIN_SHAPES + MACAULAY_FULL_SHAPES
MACAULAY_PRIMES = (101, 1009, 10007)
HOLDOUT_SHAPE = (3, 3, 3)

MODULAR_POINTS = {
    3: [1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1],
    5: [1, 4, 2, 2, 1, 1, 2, 2, 1, 1, 2, 3, 2, 4, 1, 1, 3, 2, 3, 3, 3, 4, 2, 1, 3],
    7: [1, 3, 5, 4, 6, 1, 5, 5, 2, 4, 1, 6, 5, 1, 3, 1, 3, 5, 3, 1, 4, 6, 1, 3, 3],
    11: [4, 8, 9, 1, 8, 6, 9, 3, 4, 4, 8, 4, 2, 7, 6, 3, 5, 9, 5, 7, 4, 2, 6, 7, 6],
}
MODULAR_SEARCH_TRIALS = {3: 43, 5: 7498, 7: 50803, 11: 535055}
MODULAR_SEARCH_SEED_RULE = "62000 + prime"


def pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def q_record(value: sp.Rational | Fraction | int) -> dict[str, int]:
    rational = sp.Rational(value)
    return {"numerator": int(rational.p), "denominator": int(rational.q)}


def q_text(value: sp.Rational | Fraction | int) -> str:
    rational = sp.Rational(value)
    return str(int(rational.p)) if rational.q == 1 else f"{int(rational.p)}/{int(rational.q)}"


def constraint_id(shape: Iterable[int], order: int) -> str:
    return "x".join(str(int(side)) for side in shape) + f":k{int(order)}"


def raw_lattice_constant(shape: tuple[int, int, int], order: int) -> dict:
    polynomial = even_subgraph_polynomial(cubic(*shape, periodic=False))
    log_coefficient = formal_log_coefficients(polynomial, order)[order]
    required_trace = -2 * order * log_coefficient
    return {
        "constraint_id": constraint_id(shape, order),
        "shape": list(shape),
        "order": order,
        "even_subgraph_coefficients_through_order": [
            int(value) for value in polynomial[: order + 1]
        ],
        "log_even_subgraph_coefficient": q_record(log_coefficient),
        "required_trace": q_record(required_trace),
        "equation_convention": "Tr_shape(U^k) - required_trace = 0",
    }


def canonical_hash(labels: list[dict], equations: Iterable[sp.Expr]) -> str:
    text = "\n".join(
        f"{tuple(label['shape'])}:{label['order']}:{sp.srepr(sp.expand(equation))}"
        for label, equation in zip(labels, equations, strict=True)
    )
    return hashlib.sha256(text.encode()).hexdigest()


def exact_span_record(equations: list[sp.Expr], variables: list[sp.Symbol]) -> dict:
    polynomials = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
    zero = (0,) * len(variables)
    nonconstant = sorted(
        set().union(*(set(polynomial.monoms()) - {zero} for polynomial in polynomials))
    )
    all_monomials = sorted(set(nonconstant) | {zero})
    nonconstant_matrix = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
            for monomial in nonconstant
        ]
    )
    full_matrix = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
            for monomial in all_monomials
        ]
    )
    nonconstant_rank = int(nonconstant_matrix.rank())
    full_rank = int(full_matrix.rank())
    return {
        "field": "Q",
        "equation_count": len(equations),
        "variable_count": len(variables),
        "nonconstant_monomial_count": len(nonconstant),
        "nonconstant_coefficient_rank": nonconstant_rank,
        "full_coefficient_rank": full_rank,
        "nullspace_dimension_for_nonconstant_rows": len(equations) - nonconstant_rank,
        "constant_in_equation_span": full_rank > nonconstant_rank,
    }


def exact_linear_span_coefficients(
    target: sp.Expr, equations: list[sp.Expr], variables: list[sp.Symbol]
) -> list[sp.Rational]:
    polynomials = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
    target_polynomial = sp.Poly(target, *variables, domain=sp.QQ)
    monomials = sorted(
        set(target_polynomial.monoms()).union(*(set(polynomial.monoms()) for polynomial in polynomials))
    )
    matrix = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
            for monomial in monomials
        ]
    )
    right = sp.Matrix([target_polynomial.coeff_monomial(monomial) for monomial in monomials])
    solutions = sp.linsolve((matrix, right))
    if solutions is sp.EmptySet:
        raise AssertionError(f"target is not in the exact equation span: {target}")
    solution = next(iter(solutions))
    parameters = sorted(set().union(*(value.free_symbols for value in solution)), key=str)
    substitutions = {parameter: sp.Integer(0) for parameter in parameters}
    coefficients = [sp.Rational(value.subs(substitutions)) for value in solution]
    if sp.expand(sum(c * e for c, e in zip(coefficients, equations, strict=True)) - target) != 0:
        raise AssertionError("linear-span reconstruction failed")
    return coefficients


def sparse_linear_terms(coefficients: list[sp.Rational], labels: list[dict]) -> list[dict]:
    return [
        {
            "coefficient": q_text(coefficient),
            "constraint_id": constraint_id(label["shape"], label["order"]),
            "shape": label["shape"],
            "order": label["order"],
        }
        for coefficient, label in zip(coefficients, labels, strict=True)
        if coefficient
    ]


def primary_radical_record(
    system, variables: list[sp.Symbol], equations: list[sp.Expr]
) -> tuple[dict, bool]:
    labels = [
        {"shape": list(shape), "order": order} for shape, order in system.labels
    ]
    symbols = {str(variable): variable for variable in variables}
    planes = (
        (
            "xy",
            symbols["u_mx_my"] * symbols["u_my_px"],
            symbols["u_mx_py"] * symbols["u_my_mx"] * symbols["u_py_px"],
            symbols["u_mx_mx"] * symbols["u_px_px"],
        ),
        (
            "xz",
            symbols["u_mx_mz"] * symbols["u_mz_px"] * symbols["u_pz_mx"],
            symbols["u_mx_pz"] * symbols["u_mz_mx"] * symbols["u_pz_px"],
            symbols["u_mz_mz"] * symbols["u_pz_pz"],
        ),
        (
            "yz",
            symbols["u_my_mz"]
            * symbols["u_mz_py"]
            * symbols["u_py_pz"]
            * symbols["u_pz_my"],
            symbols["u_my_pz"]
            * symbols["u_mz_my"]
            * symbols["u_py_mz"]
            * symbols["u_pz_py"],
            symbols["u_my_my"] * symbols["u_py_py"],
        ),
    )
    plane_records = []
    radical_binomials: list[sp.Expr] = []
    all_verified = True
    for name, left, right, straight in planes:
        targets = {
            "sum_relation": sp.expand(left + right + 2),
            "product_relation": sp.expand(left * right - 1),
            "straight_relation": sp.expand(straight * (left + right) + 2),
        }
        target_records = {}
        for target_name, target in targets.items():
            coefficients = exact_linear_span_coefficients(target, equations, variables)
            reconstructed = sp.expand(
                sum(
                    coefficient * equation
                    for coefficient, equation in zip(coefficients, equations, strict=True)
                )
            )
            verified = reconstructed == target
            all_verified &= verified
            target_records[target_name] = {
                "polynomial": str(target),
                "linear_combination_of_raw_constraints": sparse_linear_terms(
                    coefficients, labels
                ),
                "identity_verified": verified,
            }
        square_left = sp.expand(left * targets["sum_relation"] - targets["product_relation"])
        square_right = sp.expand(right * targets["sum_relation"] - targets["product_relation"])
        straight_identity = sp.expand(
            (straight * targets["sum_relation"] - targets["straight_relation"]) / 2
        )
        consequences_verified = (
            square_left == sp.expand((left + 1) ** 2)
            and square_right == sp.expand((right + 1) ** 2)
            and straight_identity == sp.expand(straight - 1)
        )
        all_verified &= consequences_verified
        radical_binomials.extend((sp.expand(left + 1), sp.expand(right + 1), sp.expand(straight - 1)))
        plane_records.append(
            {
                "plane": name,
                "left_turn_monomial": str(left),
                "right_turn_monomial": str(right),
                "straight_product": str(straight),
                **target_records,
                "radical_consequences": {
                    "left": f"({left} + 1)^2 = left*sum_relation - product_relation",
                    "right": f"({right} + 1)^2 = right*sum_relation - product_relation",
                    "straight": "straight_product - 1 = (straight_product*sum_relation - straight_relation)/2",
                    "identities_verified": consequences_verified,
                },
            }
        )

    exponent_rows = []
    for binomial in radical_binomials:
        polynomial = sp.Poly(binomial, *variables, domain=sp.ZZ)
        nonconstant_terms = [monomial for monomial, _ in polynomial.terms() if any(monomial)]
        if len(nonconstant_terms) != 1:
            raise AssertionError("radical relation is not monomial-minus-constant")
        exponent_rows.append(nonconstant_terms[0])
    smith = smith_normal_form(sp.Matrix(exponent_rows), domain=ZZ)
    diagonal = [abs(int(value)) for value in smith.diagonal() if value]
    saturated = diagonal == [1] * len(exponent_rows)
    all_verified &= saturated
    radical_basis = sp.groebner(
        radical_binomials, *variables, order="grevlex", domain=sp.QQ
    )
    raw_constraints_in_radical = all(
        sp.expand(radical_basis.reduce(equation)[1]) == 0 for equation in equations
    )
    all_verified &= raw_constraints_in_radical
    forced_nonzero = len(
        {
            index
            for row in exponent_rows
            for index, exponent in enumerate(row)
            if exponent
        }
    )
    return (
        {
            "claim_tag": "[LEMMA]",
            "claim": (
                "[LEMMA] On branch 11111 the nine-plane thin-box ideal has the "
                "displayed prime radical: six oriented planar turn products equal -1 "
                "and three opposite straight products equal 1."
            ),
            "raw_constraint_count": len(equations),
            "raw_constraint_labels": labels,
            "planes": plane_records,
            "radical_binomials": [str(binomial) for binomial in radical_binomials],
            "all_raw_constraints_reduce_to_zero_mod_radical": raw_constraints_in_radical,
            "exponent_lattice_smith_diagonal": diagonal,
            "exponent_lattice_rank": len(diagonal),
            "saturated_exponent_lattice": saturated,
            "prime_radical_dimension": len(variables) - len(diagonal),
            "variables_forced_nonzero_by_monomial_relations": forced_nonzero,
            "scope": (
                "This is an exact radical reduction with one identified prime component. "
                "It is not an emptiness certificate after the genuinely three-dimensional "
                "constraints are added."
            ),
        },
        all_verified,
    )


def independent_polynomial_indices(
    expressions: list[sp.Expr], variables: list[sp.Symbol]
) -> list[int]:
    polynomials = [sp.Poly(expression, *variables, domain=sp.QQ) for expression in expressions]
    monomials = sorted(set().union(*(set(polynomial.monoms()) for polynomial in polynomials)))
    matrix = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for monomial in monomials]
            for polynomial in polynomials
        ]
    )
    return list(matrix.T.rref()[1])


def macaulay_columns(
    polynomials: list[sp.Poly], variable_count: int
) -> list[dict[tuple[int, ...], int]]:
    zero = (0,) * variable_count
    multipliers = [zero] + [
        tuple(1 if index == coordinate else 0 for index in range(variable_count))
        for coordinate in range(variable_count)
    ]
    columns = []
    for polynomial in polynomials:
        terms = {monomial: int(coefficient) for monomial, coefficient in polynomial.terms()}
        for multiplier in multipliers:
            columns.append(
                {
                    tuple(a + b for a, b in zip(monomial, multiplier, strict=True)): coefficient
                    for monomial, coefficient in terms.items()
                }
            )
    return columns


def sparse_matrix_hash(columns: list[dict[tuple[int, ...], int]]) -> str:
    digest = hashlib.sha256()
    for column_index, column in enumerate(columns):
        digest.update(f"C{column_index}:".encode())
        for monomial, coefficient in sorted(column.items()):
            digest.update(
                (f"{coefficient}@{','.join(map(str, monomial))};").encode()
            )
    return digest.hexdigest()


def modular_column_rank_and_constant(
    columns: list[dict[tuple[int, ...], int]], variable_count: int, prime: int
) -> dict:
    basis: dict[tuple[int, ...], dict[tuple[int, ...], int]] = {}
    pivot_product = 1
    pivot_digest = hashlib.sha256()
    for source in columns:
        vector = {
            monomial: coefficient % prime
            for monomial, coefficient in source.items()
            if coefficient % prime
        }
        while vector:
            pivot = max(vector, key=lambda monomial: (sum(monomial), monomial))
            value = vector[pivot]
            if pivot not in basis:
                pivot_product = pivot_product * value % prime
                pivot_digest.update((",".join(map(str, pivot)) + ";").encode())
                inverse = pow(value, -1, prime)
                vector = {
                    monomial: coefficient * inverse % prime
                    for monomial, coefficient in vector.items()
                    if coefficient % prime
                }
                basis[pivot] = vector
                break
            pivot_column = basis[pivot]
            for monomial, coefficient in pivot_column.items():
                reduced = (vector.get(monomial, 0) - value * coefficient) % prime
                if reduced:
                    vector[monomial] = reduced
                elif monomial in vector:
                    del vector[monomial]

    target = {(0,) * variable_count: 1}
    while target:
        pivot = max(target, key=lambda monomial: (sum(monomial), monomial))
        value = target[pivot]
        if pivot not in basis:
            return {
                "prime": prime,
                "column_rank": len(basis),
                "augmented_rank": len(basis) + 1,
                "constant_target_independent": True,
                "constant_residual_pivot": list(pivot),
                "constant_residual_coefficient": value,
                "pivot_product_mod_p": pivot_product,
                "pivot_monomial_sha256": pivot_digest.hexdigest(),
            }
        pivot_column = basis[pivot]
        for monomial, coefficient in pivot_column.items():
            reduced = (target.get(monomial, 0) - value * coefficient) % prime
            if reduced:
                target[monomial] = reduced
            elif monomial in target:
                del target[monomial]
    return {
        "prime": prime,
        "column_rank": len(basis),
        "augmented_rank": len(basis),
        "constant_target_independent": False,
        "pivot_product_mod_p": pivot_product,
        "pivot_monomial_sha256": pivot_digest.hexdigest(),
    }


def selected_zero_free_lift(
    prime: int,
    initial: list[int],
    polynomials,
    derivatives,
) -> dict:
    point = [value % prime for value in initial]
    jacobian = jacobian_mod(derivatives, point, prime)
    rank, pivots, _ = rref_mod(jacobian, prime)
    steps = []
    precision = 1
    while precision < 3:
        old_modulus = prime**precision
        new_modulus = old_modulus * prime
        residuals = eval_system_mod(polynomials, point, new_modulus)
        if any(residual % old_modulus for residual in residuals):
            raise AssertionError("selected p-adic path lost its current precision")
        right_hand_side = [-(residual // old_modulus) % prime for residual in residuals]
        solve = solve_rectangular_mod(jacobian, right_hand_side, prime)
        if not solve["solvable"]:
            left_kernel = solve["left_kernel_vector"]
            kernel_times_jacobian = [
                sum(left_kernel[row] * jacobian[row][column] for row in range(len(jacobian)))
                % prime
                for column in range(len(point))
            ]
            pairing = (
                sum(
                    left_kernel[row] * right_hand_side[row]
                    for row in range(len(jacobian))
                )
                % prime
            )
            if any(kernel_times_jacobian) or pairing == 0:
                raise AssertionError("invalid selected-path nonlift certificate")
            return {
                "status": "selected_residue_has_no_next_digit",
                "jacobian_rank_mod_p": rank,
                "pivot_columns": pivots,
                "last_solved_exponent": precision,
                "attempted_exponent": precision + 1,
                "modulus": old_modulus,
                "residues": point,
                "steps": steps,
                "certificate": {
                    "left_kernel_vector": left_kernel,
                    "left_kernel_times_jacobian_mod_p": kernel_times_jacobian,
                    "pairing_with_rhs_mod_p": pairing,
                    "linearized_rhs_mod_p": right_hand_side,
                },
                "scope": (
                    "[THEOREM] This displayed residue has no next p-adic digit. "
                    "[UNRESOLVED] Other points and other singular sections are not excluded."
                ),
            }
        correction = instantiate_linear_solution(
            solve, [0] * len(solve["free_columns"]), prime, len(point)
        )
        point = [
            (value + old_modulus * digit) % new_modulus
            for value, digit in zip(point, correction, strict=True)
        ]
        if any(eval_system_mod(polynomials, point, new_modulus)):
            raise AssertionError("selected zero-free digit did not lift exactly")
        steps.append(
            {
                "from_exponent": precision,
                "to_exponent": precision + 1,
                "free_coordinates_set_to_zero": len(solve["free_columns"]),
                "residuals_zero": True,
            }
        )
        precision += 1
    return {
        "status": "selected_path_reached_budget",
        "jacobian_rank_mod_p": rank,
        "pivot_columns": pivots,
        "last_solved_exponent": precision,
        "modulus": prime**precision,
        "residues": point,
        "steps": steps,
        "scope": "[UNRESOLVED] A finite selected p-adic section is not a rational point.",
    }


def main() -> None:
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    prior_rows = prior["data"]["branches"]
    prior_by_branch = {row["branch"]: row for row in prior_rows}
    labels = [format(index, "05b") for index in range(32)]
    if [row["branch"] for row in prior_rows] != labels:
        raise AssertionError("e50 branch ordering changed")

    scan_system = full_weight_finite_system(
        SCAN_THIN_SHAPES, SCAN_ORDERS, gauge_fix=False
    )
    scan_labels = [
        {"shape": list(shape), "order": order} for shape, order in scan_system.labels
    ]
    branch_records = []
    new_certificate_count = 0
    new_identities_verified = True
    used_constraints: set[tuple[tuple[int, int, int], int]] = set()
    generic_scan = None
    equivalence_groups: dict[tuple[tuple[tuple[int, ...], int, int], ...], list[str]] = {}

    for branch in labels:
        prior_row = prior_by_branch[branch]
        if prior_row["status"] == "EMPTY_OVER_Q":
            branch_records.append(
                {
                    "branch": branch,
                    "status": "EMPTY_OVER_Q",
                    "claim_tag": "[THEOREM]",
                    "claim": prior_row["claim"],
                    "certificate_origin": "e50_imported_unchanged",
                    "empty_certificate": prior_row["empty_certificate"],
                }
            )
            continue

        variables, equations = branch_equations(scan_system, branch)
        search = constant_certificate_search(equations, variables, scan_labels)
        if search["certificate_found"]:
            certificate = search["certificate"]
            coefficients = certificate["coefficient_vector"]
            combination = sp.expand(
                sum(
                    coefficient * equation
                    for coefficient, equation in zip(coefficients, equations, strict=True)
                )
            )
            verified = (
                not combination.free_symbols
                and combination == certificate["constant"]
                and certificate["constant"] != 0
            )
            new_identities_verified &= verified
            sparse_terms = []
            signature = []
            for index, coefficient in enumerate(coefficients):
                if not coefficient:
                    continue
                shape = tuple(scan_labels[index]["shape"])
                order = int(scan_labels[index]["order"])
                used_constraints.add((shape, order))
                sparse_terms.append(
                    {
                        "coefficient": int(coefficient),
                        "constraint_id": constraint_id(shape, order),
                        "shape": list(shape),
                        "order": order,
                    }
                )
                signature.append((shape, order, int(coefficient)))
            equivalence_groups.setdefault(tuple(signature), []).append(branch)
            branch_records.append(
                {
                    "branch": branch,
                    "status": "EMPTY_OVER_Q",
                    "claim_tag": "[THEOREM]",
                    "claim": (
                        f"[THEOREM] Branch {branch} is empty over Q because the displayed "
                        "integer combination of raw thin-box equations is the nonzero "
                        f"constant {certificate['constant']}."
                    ),
                    "certificate_origin": "e62_new_thin_box_identity",
                    "empty_certificate": {
                        "claim_tag": "[THEOREM]",
                        "type": "integer_constant_identity",
                        "terms": sparse_terms,
                        "constant": int(certificate["constant"]),
                        "identity_verified": verified,
                        "consequence": (
                            "At a common zero the left side is zero, contradicting the "
                            "displayed nonzero integer constant."
                        ),
                    },
                }
            )
            new_certificate_count += 1
        else:
            generic_scan = exact_span_record(equations, variables)
            generic_scan.update(
                {
                    "claim_tag": "[THEOREM]",
                    "claim": (
                        "[THEOREM] Within the recorded order<=10 thin-box equation span, "
                        "branch 11111 has no degree-zero rational Nullstellensatz certificate."
                    ),
                    "system_sha256": canonical_hash(scan_labels, equations),
                }
            )
            branch_records.append(
                {
                    "branch": branch,
                    "status": "UNRESOLVED",
                    "claim_tag": "[UNRESOLVED]",
                    "claim": (
                        "[UNRESOLVED] Branch 11111 has neither an exact emptiness certificate "
                        "nor an exactly verified rational solution of all construction and "
                        "holdout equations."
                    ),
                }
            )

    if generic_scan is None:
        generic_scan = {"claim_tag": "[COMPUTATION]", "claim": "No survivor remained to scan."}

    raw_constants = [
        raw_lattice_constant(shape, order)
        for shape, order in sorted(used_constraints)
    ]

    primary_system = full_weight_finite_system(
        PRIMARY_THIN_SHAPES, BASE_ORDERS, gauge_fix=False
    )
    primary_variables, primary_equations = branch_equations(primary_system, "11111")
    primary_record, primary_verified = primary_radical_record(
        primary_system, primary_variables, primary_equations
    )

    macaulay_system = full_weight_finite_system(
        MACAULAY_SHAPES, BASE_ORDERS, gauge_fix=False
    )
    macaulay_variables, macaulay_equations = branch_equations(macaulay_system, "11111")
    macaulay_labels = [
        {"shape": list(shape), "order": order} for shape, order in macaulay_system.labels
    ]
    basis_indices = independent_polynomial_indices(macaulay_equations, macaulay_variables)
    basis_polynomials = [
        sp.Poly(macaulay_equations[index], *macaulay_variables, domain=sp.ZZ)
        for index in basis_indices
    ]
    columns = macaulay_columns(basis_polynomials, len(macaulay_variables))
    matrix_hash = sparse_matrix_hash(columns)
    modular_ranks = [
        modular_column_rank_and_constant(columns, len(macaulay_variables), prime)
        for prime in MACAULAY_PRIMES
    ]
    column_count = len(columns)
    degree_one_certified_absent = all(
        record["column_rank"] == column_count
        and record["augmented_rank"] == column_count + 1
        and record["constant_target_independent"]
        and record["pivot_product_mod_p"]
        and record["constant_residual_coefficient"]
        for record in modular_ranks
    )
    macaulay_record = {
        "claim_tag": "[THEOREM]",
        "claim": (
            "[THEOREM] For the exact recorded 42-constraint system on branch 11111, "
            "no rational Nullstellensatz certificate with multiplier degree at most one exists."
        ),
        "constraint_count": len(macaulay_equations),
        "constraint_labels": macaulay_labels,
        "exact_Q_constraint_span_rank": len(basis_indices),
        "basis_constraint_indices": basis_indices,
        "basis_constraint_labels": [macaulay_labels[index] for index in basis_indices],
        "variable_order": [str(variable) for variable in macaulay_variables],
        "multiplier_basis": ["1"] + [str(variable) for variable in macaulay_variables],
        "multiplier_degree_bound": 1,
        "column_count": column_count,
        "integer_macaulay_matrix_sha256": matrix_hash,
        "modular_same_matrix_rank_certificates": modular_ranks,
        "characteristic_zero_rank_argument": (
            "Each reduction has full column rank 546, while the constant-augmented "
            "matrix has rank 547. Column counts are exact upper bounds over Q for "
            "the same integer matrices, so their Q-ranks are exactly 546 and 547."
        ),
        "certificate_absent_at_budget": degree_one_certified_absent,
        "scope": (
            "This excludes only multiplier degree <=1 for this finite exact constraint "
            "catalogue. It is not a Q-emptiness or nonemptiness theorem."
        ),
    }

    macaulay_sparse = sparse_polynomials(macaulay_equations, macaulay_variables)
    derivatives = derivative_table(macaulay_sparse, len(macaulay_variables))
    holdout_system = full_weight_finite_system((HOLDOUT_SHAPE,), BASE_ORDERS, gauge_fix=False)
    holdout_variables, holdout_equations = branch_equations(holdout_system, "11111")
    if [str(variable) for variable in holdout_variables] != [
        str(variable) for variable in macaulay_variables
    ]:
        raise AssertionError("holdout variable order differs from construction order")
    holdout_sparse = sparse_polynomials(holdout_equations, holdout_variables)
    modular_records = []
    modular_checks_ok = True
    for prime, initial in MODULAR_POINTS.items():
        construction_residues = eval_system_mod(macaulay_sparse, initial, prime)
        exact_construction = [
            int(sp.Poly(equation, *macaulay_variables).eval(dict(zip(macaulay_variables, initial, strict=True))))
            for equation in macaulay_equations
        ]
        exact_holdout = [
            int(sp.Poly(equation, *holdout_variables).eval(dict(zip(holdout_variables, initial, strict=True))))
            for equation in holdout_equations
        ]
        lift = selected_zero_free_lift(
            prime, initial, macaulay_sparse, derivatives
        )
        modulus = int(lift["modulus"])
        holdout_residues = eval_system_mod(
            holdout_sparse, [int(value) for value in lift["residues"]], modulus
        )
        initial_verified = not any(construction_residues)
        certificate = lift.get("certificate", {})
        lift_verified = (
            lift["status"] == "selected_residue_has_no_next_digit"
            and not any(certificate.get("left_kernel_times_jacobian_mod_p", [1]))
            and certificate.get("pairing_with_rhs_mod_p", 0) != 0
        )
        modular_checks_ok &= initial_verified and lift_verified
        modular_records.append(
            {
                "prime": prime,
                "search": {
                    "method": "deterministic nonzero radical-parameter sampling",
                    "seed_rule": MODULAR_SEARCH_SEED_RULE,
                    "accepted_trial": MODULAR_SEARCH_TRIALS[prime],
                },
                "values_in_variable_order": initial,
                "construction_residues_mod_p": construction_residues,
                "construction_solution_verified": initial_verified,
                "jacobian_rank_mod_p": lift["jacobian_rank_mod_p"],
                "selected_zero_free_padic_path": lift,
                "holdout_3x3x3_residues_at_last_precision": holdout_residues,
                "least_residue_integer_candidate": {
                    "construction_residuals": exact_construction,
                    "holdout_residuals": exact_holdout,
                    "passed_all_exact_equations": not any(exact_construction + exact_holdout),
                },
                "scope": (
                    "The finite-field point proves only nonemptiness of this reduction. "
                    "The selected p-adic nonlift certificate excludes only the displayed "
                    "residue section, not the branch."
                ),
            }
        )

    generic_branch = next(record for record in branch_records if record["branch"] == "11111")
    generic_branch["degree_zero_thin_span"] = generic_scan
    generic_branch["thin_primary_radical"] = primary_record
    generic_branch["bounded_nullstellensatz"] = macaulay_record
    generic_branch["modular_and_padic_diagnostics"] = modular_records
    generic_branch["exact_rational_candidate_count"] = 0
    generic_branch["candidate_gate"] = (
        "Any reconstructed rational candidate must vanish on all 42 construction "
        "equations and all three exact 3x3x3 holdout equations; no candidate reached "
        "rational reconstruction precision."
    )

    status_counts = {
        status: sum(record["status"] == status for record in branch_records)
        for status in ("EMPTY_OVER_Q", "RATIONAL_POINT_FOUND", "UNRESOLVED")
    }
    all_empty = status_counts["EMPTY_OVER_Q"] == 32
    raw_ids = {record["constraint_id"] for record in raw_constants}
    raw_coverage = all(
        term["constraint_id"] in raw_ids
        for record in branch_records
        if record.get("certificate_origin") == "e62_new_thin_box_identity"
        for term in record["empty_certificate"]["terms"]
    )
    imported_count = sum(
        record.get("certificate_origin") == "e50_imported_unchanged"
        for record in branch_records
    )
    checks = [
        {
            "name": "all_32_e50_branches_imported_in_order",
            "passed": [record["branch"] for record in branch_records] == labels,
            "detail": {"first": labels[0], "last": labels[-1], "count": len(labels)},
        },
        {
            "name": "seventeen_prior_empty_certificates_imported_unchanged",
            "passed": imported_count == 17,
            "detail": imported_count,
        },
        {
            "name": "fourteen_new_integer_constant_identities_verified_exactly",
            "passed": new_certificate_count == 14 and new_identities_verified,
            "detail": {
                "count": new_certificate_count,
                "branches": [
                    record["branch"]
                    for record in branch_records
                    if record.get("certificate_origin") == "e62_new_thin_box_identity"
                ],
            },
        },
        {
            "name": "every_new_certificate_has_raw_lattice_constant_data",
            "passed": raw_coverage,
            "detail": sorted(raw_ids),
        },
        {
            "name": "generic_thin_span_has_no_degree_zero_constant",
            "passed": generic_scan.get("constant_in_equation_span") is False
            and generic_scan.get("full_coefficient_rank")
            == generic_scan.get("nonconstant_coefficient_rank"),
            "detail": generic_scan,
        },
        {
            "name": "thin_primary_radical_relations_and_smith_form_verified",
            "passed": primary_verified,
            "detail": {
                "smith_diagonal": primary_record["exponent_lattice_smith_diagonal"],
                "dimension": primary_record["prime_radical_dimension"],
            },
        },
        {
            "name": "degree_one_nullstellensatz_space_certified_empty_over_Q",
            "passed": degree_one_certified_absent,
            "detail": {
                "matrix_sha256": matrix_hash,
                "ranks": modular_ranks,
            },
        },
        {
            "name": "four_modular_points_and_selected_padic_nonlifts_verified",
            "passed": modular_checks_ok,
            "detail": {
                str(record["prime"]): {
                    "rank": record["jacobian_rank_mod_p"],
                    "last_solved_exponent": record["selected_zero_free_padic_path"][
                        "last_solved_exponent"
                    ],
                }
                for record in modular_records
            },
        },
        {
            "name": "status_map_is_31_empty_and_one_honest_survivor",
            "passed": status_counts
            == {"EMPTY_OVER_Q": 31, "RATIONAL_POINT_FOUND": 0, "UNRESOLVED": 1},
            "detail": status_counts,
        },
        {
            "name": "global_no_go_withheld_while_branch_11111_survives",
            "passed": all_empty or generic_branch["status"] == "UNRESOLVED",
            "detail": {"all_empty": all_empty, "survivor": generic_branch["branch"]},
        },
    ]

    output = {
        "provenance": {
            "script": "experiments/e62_kw_close.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "prior_artifact": "results/kac_ward/branches.json",
            "prior_artifact_sha256": hashlib.sha256(PRIOR.read_bytes()).hexdigest(),
            "method": (
                "exact ZZ/QQ thin-box coefficient spans through k=10; primitive integer "
                "constant identities; exact radical/SNF reduction; same-matrix modular "
                "Macaulay rank with characteristic-zero column-count upper bounds; exact "
                "finite-field substitution and selected singular p-adic nonlift certificates"
            ),
            "no_float_decisions": True,
        },
        "data": {
            "headline": (
                "[COMPUTATION] The exact gauge-tree map now has 31 EMPTY_OVER_Q branches; "
                "[UNRESOLVED] branch 11111 remains."
            ),
            "total_branch_count": 32,
            "status_counts": status_counts,
            "gauge_tree_order": [pair_name(pair) for pair in DIRECTION_GAUGE_TREE],
            "bit_convention": "0=fixed tree weight vanishes; 1=normalized to one",
            "systematic_thin_scan": {
                "shapes": [list(shape) for shape in SCAN_THIN_SHAPES],
                "orders": list(SCAN_ORDERS),
                "equation_count": len(scan_system.equations),
                "generic_branch_degree_zero_result": generic_scan,
            },
            "raw_lattice_constants_for_new_certificates": raw_constants,
            "new_certificate_equivalence_classes": [
                {
                    "branches": branches,
                    "terms": [
                        {
                            "shape": list(shape),
                            "order": order,
                            "coefficient": coefficient,
                        }
                        for shape, order, coefficient in signature
                    ],
                    "interpretation": (
                        "These branches have exactly the same branch-reduced constant "
                        "difference, so one raw identity is reused without a new inference."
                    ),
                }
                for signature, branches in sorted(equivalence_groups.items())
            ],
            "branches": branch_records,
            "global_claim": (
                "[THEOREM] Every enumerated gauge-tree branch is empty over Q."
                if all_empty
                else "[UNRESOLVED] No global no-go is claimed because branch 11111 remains unresolved."
            ),
            "scope": (
                "[THEOREM] Each EMPTY_OVER_Q status is only for its named gauge-tree chart. "
                "[UNRESOLVED] The 32 charts are not all support patterns of the 30-weight "
                "family, and the full three-dimensional scalar Kac--Ward problem remains open."
            ),
        },
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"failed checks: {failed}")
    print("PASS")


if __name__ == "__main__":
    main()
