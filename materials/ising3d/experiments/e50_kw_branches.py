"""Exact branch ledger for the gauge-tree charts of the scalar Kac--Ward family.

The experiment separates decisive characteristic-zero certificates from bounded
computations.  A timeout and a failed rational reconstruction remain explicitly
non-decisive.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import multiprocessing as mp
import queue
import random
import time
from datetime import datetime, timezone
from fractions import Fraction
from functools import reduce
from pathlib import Path

import sympy as sp
from sympy.polys.orderings import ProductOrder, grevlex, lex

from e42_kw_hensel import (
    derivative_table,
    eval_system_mod,
    hensel_lift_full_row_rank,
    jacobian_mod,
    load_stored_system,
    rank_deficient_first_lift_attempt,
    rational_reconstruct,
    reconstruct_lift,
    rref_mod,
    target_exponent,
)
from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    full_weight_finite_system,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/kac_ward/full_family.json"
OUTPUT = ROOT / "results/kac_ward/branches.json"
CONSTRUCTION_SHAPES = ((3, 3, 2), (2, 2, 3))
# The first two shapes preserve the wave-6 certificate as the first pair tried.
# The ten boxes separate every axis permutation used in the wave-4 ten-box fit.
CERTIFICATE_SHAPES = (
    (3, 3, 2),
    (2, 2, 3),
    (3, 2, 3),
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
    (4, 2, 2),
    (2, 4, 2),
    (2, 2, 4),
    (2, 3, 3),
)
HOLDOUT_SHAPES = ((3, 3, 3), (4, 2, 3))
ORDERS = (4, 6, 8)
MINIMUM_PADIC_EXPONENT = 60
MINIMUM_DECIMAL_DIGITS = 60
GROEBNER_TIMEOUT_SECONDS = 60
SINGULAR_TRIALS_PER_DIGIT = 256

SparsePolynomial = tuple[tuple[int, tuple[int, ...]], ...]


def pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def all_branch_labels() -> list[str]:
    return [
        "".join(map(str, bits))
        for bits in itertools.product((0, 1), repeat=len(DIRECTION_GAUGE_TREE))
    ]


def branch_equations(ungauged, branch: str) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    if len(branch) != len(DIRECTION_GAUGE_TREE) or set(branch) - {"0", "1"}:
        raise ValueError(f"invalid branch label {branch!r}")
    symbol_map = dict(ungauged.all_symbols)
    substitutions = {
        symbol_map[pair]: int(bit)
        for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
    }
    equations = [sp.expand(equation.subs(substitutions)) for equation in ungauged.equations]
    active = sorted(set().union(*(equation.free_symbols for equation in equations)), key=str)
    return active, equations


def sparse_polynomials(
    expressions: list[sp.Expr] | tuple[sp.Expr, ...], variables: list[sp.Symbol]
) -> list[SparsePolynomial]:
    records: list[SparsePolynomial] = []
    for expression in expressions:
        polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
        records.append(
            tuple(
                (int(coefficient), tuple(int(exponent) for exponent in monomial))
                for monomial, coefficient in polynomial.terms()
            )
        )
    return records


def primitive_integer_certificate(
    coefficients: list[sp.Rational], constant: sp.Rational
) -> tuple[list[int], int]:
    denominators = [int(value.q) for value in coefficients] + [int(constant.q)]
    scale = int(sp.ilcm(*denominators))
    integers = [int(value * scale) for value in coefficients]
    integer_constant = int(constant * scale)
    divisor = reduce(
        math.gcd,
        (abs(value) for value in integers + [integer_constant] if value),
        0,
    )
    if not divisor:
        raise AssertionError("a certificate cannot be identically zero")
    integers = [value // divisor for value in integers]
    integer_constant //= divisor
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
        integer_constant = -integer_constant
    return integers, integer_constant


def constant_certificate_search(
    equations: list[sp.Expr],
    variables: list[sp.Symbol],
    labels: list[dict],
) -> dict:
    polynomials = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
    zero_monomial = (0,) * len(variables)
    coefficient_maps = []
    constants = []
    all_nonconstant_monomials: set[tuple[int, ...]] = set()
    for polynomial in polynomials:
        mapping = {
            monomial: sp.Rational(coefficient)
            for monomial, coefficient in polynomial.terms()
            if monomial != zero_monomial
        }
        coefficient_maps.append(mapping)
        all_nonconstant_monomials.update(mapping)
        constants.append(sp.Rational(polynomial.coeff_monomial(zero_monomial)))

    monomials = sorted(all_nonconstant_monomials)
    matrix = sp.Matrix(
        [
            [mapping.get(monomial, sp.Rational(0)) for mapping in coefficient_maps]
            for monomial in monomials
        ]
    )
    nullspace = matrix.nullspace()

    candidate: tuple[list[sp.Rational], sp.Rational, str] | None = None
    # Pair certificates are preferred because they are shortest to audit.  This
    # deterministic scan recovers 2*F_(3x3x2,4)-3*F_(2x2x3,4)=56 first where valid.
    for left in range(len(polynomials)):
        if candidate is not None:
            break
        left_map = coefficient_maps[left]
        if not left_map:
            continue
        for right in range(left + 1, len(polynomials)):
            right_map = coefficient_maps[right]
            if left_map.keys() != right_map.keys() or not right_map:
                continue
            first_monomial = next(iter(left_map))
            ratio = left_map[first_monomial] / right_map[first_monomial]
            if any(left_map[key] != ratio * right_map[key] for key in left_map):
                continue
            constant = constants[left] - ratio * constants[right]
            if not constant:
                continue
            coefficients = [sp.Rational(0)] * len(polynomials)
            coefficients[left] = sp.Rational(1)
            coefficients[right] = -ratio
            candidate = (coefficients, constant, "two_equation_constant_identity")
            break

    if candidate is None:
        for index, mapping in enumerate(coefficient_maps):
            if not mapping and constants[index]:
                coefficients = [sp.Rational(0)] * len(polynomials)
                coefficients[index] = sp.Rational(1)
                candidate = (
                    coefficients,
                    constants[index],
                    "single_equation_nonzero_constant",
                )
                break

    if candidate is None:
        constant_vector = sp.Matrix(constants)
        for vector in nullspace:
            constant = (constant_vector.T * vector)[0]
            if constant:
                candidate = (
                    [sp.Rational(value) for value in vector],
                    sp.Rational(constant),
                    "full_Q_linear_span_constant_identity",
                )
                break

    record = {
        "claim_tag": "[COMPUTATION]",
        "method": "exact Q-nullspace of all nonconstant monomial coefficients",
        "equation_count": len(equations),
        "nonconstant_monomial_count": len(monomials),
        "nonconstant_coefficient_rank": len(equations) - len(nullspace),
        "nullspace_dimension": len(nullspace),
        "certificate_found": candidate is not None,
    }
    if candidate is None:
        record["conclusion"] = (
            "[UNRESOLVED] No nonzero constant lies in the Q-linear span of these "
            "selected constraints; polynomial-multiplier certificates are not excluded."
        )
        return record

    rational_coefficients, rational_constant, certificate_type = candidate
    coefficients, constant = primitive_integer_certificate(
        rational_coefficients, rational_constant
    )
    combination = sp.expand(
        sum(
            coefficient * equation
            for coefficient, equation in zip(coefficients, equations, strict=True)
        )
    )
    if combination.free_symbols or combination != constant or constant == 0:
        raise AssertionError("invalid constant certificate produced by exact elimination")
    record["certificate"] = {
        "claim_tag": "[THEOREM]",
        "type": certificate_type,
        "coefficient_vector": coefficients,
        "terms": [
            {
                "coefficient": coefficient,
                "shape": labels[index]["shape"],
                "order": labels[index]["order"],
            }
            for index, coefficient in enumerate(coefficients)
            if coefficient
        ],
        "constant": constant,
        "identity_verified": True,
        "consequence": (
            "[THEOREM] At a common zero the left side would be zero, contradicting "
            "the displayed nonzero integer constant; this branch is EMPTY_OVER_Q."
        ),
    }
    return record


def independent_equation_indices(
    equations: list[sp.Expr], variables: list[sp.Symbol]
) -> list[int]:
    polynomials = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
    monomials = sorted(set().union(*(set(polynomial.monoms()) for polynomial in polynomials)))
    rows = [
        [polynomial.coeff_monomial(monomial) for monomial in monomials]
        for polynomial in polynomials
    ]
    return list(sp.Matrix(rows).T.rref()[1])


def groebner_worker(branch: str, result_queue) -> None:
    started = time.monotonic()
    try:
        ungauged = full_weight_finite_system(
            CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
        )
        variables, equations = branch_equations(ungauged, branch)
        indices = independent_equation_indices(equations, variables)
        independent = [equations[index] for index in indices]
        linear = [
            variable
            for variable in variables
            if max(sp.degree(equation, variable) for equation in independent) <= 1
            and any(sp.degree(equation, variable) == 1 for equation in independent)
        ]
        nonlinear = [variable for variable in variables if variable not in linear]
        ordered = linear + nonlinear
        split = len(linear)
        order = ProductOrder(
            (lex, lambda monomial: monomial[:split]),
            (grevlex, lambda monomial: monomial[split:]),
        )
        basis = sp.groebner(
            independent,
            *ordered,
            domain=sp.QQ,
            order=order,
            method="f5b",
        )
        basis_strings = [str(polynomial.as_expr()) for polynomial in basis.polys]
        result_queue.put(
            {
                "status": "completed",
                "elapsed_seconds": time.monotonic() - started,
                "field": "Q",
                "input_equation_indices": indices,
                "input_equation_count_after_exact_row_reduction": len(indices),
                "active_variables": [str(variable) for variable in variables],
                "linear_elimination_block": [str(variable) for variable in linear],
                "remaining_block": [str(variable) for variable in nonlinear],
                "order": "ProductOrder(lex(linear block), grevlex(remaining block))",
                "basis": basis_strings,
                "basis_size": len(basis_strings),
                "basis_term_counts": [len(polynomial.terms()) for polynomial in basis.polys],
                "basis_sha256": hashlib.sha256("\n".join(basis_strings).encode()).hexdigest(),
                "contains_one": basis_strings == ["1"],
            }
        )
    except BaseException as error:
        result_queue.put(
            {
                "status": "error",
                "elapsed_seconds": time.monotonic() - started,
                "error_type": type(error).__name__,
                "error": str(error),
            }
        )


def bounded_groebner_attempt(branch: str) -> dict:
    context = mp.get_context("spawn")
    result_queue = context.Queue()
    process = context.Process(target=groebner_worker, args=(branch, result_queue))
    started = time.monotonic()
    process.start()
    try:
        record = result_queue.get(timeout=GROEBNER_TIMEOUT_SECONDS)
    except queue.Empty:
        if process.is_alive():
            process.terminate()
            process.join(10)
            return {
                "claim_tag": "[UNRESOLVED]",
                "branch": branch,
                "status": "timeout",
                "field": "Q",
                "limit_seconds": GROEBNER_TIMEOUT_SECONDS,
                "elapsed_seconds": time.monotonic() - started,
                "scope": "[UNRESOLVED] A hard timeout is not evidence of nonemptiness.",
            }
        record = {
            "status": "error",
            "error_type": "MissingWorkerResult",
            "error": f"worker exited with code {process.exitcode} without a result",
        }
    process.join(10)
    if process.is_alive():
        process.terminate()
        process.join(10)
    record.update(
        {
            "claim_tag": "[COMPUTATION]",
            "branch": branch,
            "limit_seconds": GROEBNER_TIMEOUT_SECONDS,
        }
    )
    if record.get("status") == "completed":
        record["scope"] = (
            "[THEOREM] The exact Q-basis [1] proves this six-equation branch empty."
            if record.get("contains_one")
            else "[UNRESOLVED] A proper finite basis describes only the selected ideal."
        )
    return record


def solve_rectangular_mod(
    matrix: list[list[int]], right_hand_side: list[int], prime: int
) -> dict:
    row_count = len(matrix)
    column_count = len(matrix[0])
    augmented = [
        [value % prime for value in row]
        + [right_hand_side[index] % prime]
        + [1 if index == other else 0 for other in range(row_count)]
        for index, row in enumerate(matrix)
    ]
    rank = 0
    pivots: list[int] = []
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if augmented[row][column]),
            None,
        )
        if pivot is None:
            continue
        augmented[rank], augmented[pivot] = augmented[pivot], augmented[rank]
        inverse = pow(augmented[rank][column], -1, prime)
        augmented[rank] = [(value * inverse) % prime for value in augmented[rank]]
        for row in range(row_count):
            if row == rank or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (left - factor * right) % prime
                for left, right in zip(augmented[row], augmented[rank], strict=True)
            ]
        pivots.append(column)
        rank += 1
        if rank == row_count:
            break
    for row in augmented:
        if all(row[column] == 0 for column in range(column_count)) and row[column_count]:
            return {
                "solvable": False,
                "rank": rank,
                "left_kernel_vector": row[column_count + 1 :],
                "pairing_with_rhs_mod_p": row[column_count],
            }
    free_columns = [column for column in range(column_count) if column not in pivots]
    reduced_rows = [
        {
            "pivot": pivot,
            "coefficients": augmented[row][:column_count],
            "right_hand_side": augmented[row][column_count],
        }
        for row, pivot in enumerate(pivots)
    ]
    return {
        "solvable": True,
        "rank": rank,
        "pivots": pivots,
        "free_columns": free_columns,
        "reduced_rows": reduced_rows,
    }


def instantiate_linear_solution(
    solve: dict, free_values: list[int], prime: int, column_count: int
) -> list[int]:
    if len(free_values) != len(solve["free_columns"]):
        raise ValueError("wrong free-coordinate count")
    solution = [0] * column_count
    for column, value in zip(solve["free_columns"], free_values, strict=True):
        solution[column] = value % prime
    for row in solve["reduced_rows"]:
        solution[row["pivot"]] = (
            row["right_hand_side"]
            - sum(
                row["coefficients"][column] * solution[column]
                for column in solve["free_columns"]
            )
        ) % prime
    return solution


def singular_hensel_section(
    prime: int,
    initial: list[int],
    polynomials: list[SparsePolynomial],
    derivatives: list[list[SparsePolynomial]],
    exponent: int,
    *,
    trials_per_digit: int = SINGULAR_TRIALS_PER_DIGIT,
) -> dict:
    """Lift one deterministic singular path with exact one-digit lookahead.

    A finite lookahead search is not a classification of all singular lifts.
    Each accepted digit is nevertheless an exact solution modulo its recorded
    power, and the final residue vector is independently rechecked.
    """

    point = [value % prime for value in initial]
    jacobian = jacobian_mod(derivatives, point, prime)
    steps = []
    for precision in range(1, exponent):
        old_modulus = prime**precision
        new_modulus = old_modulus * prime
        residuals = eval_system_mod(polynomials, point, new_modulus)
        if any(value % old_modulus for value in residuals):
            raise AssertionError("singular path lost its current precision")
        right_hand_side = [-(value // old_modulus) % prime for value in residuals]
        solve = solve_rectangular_mod(jacobian, right_hand_side, prime)
        if not solve["solvable"]:
            return {
                "status": "selected_path_blocked",
                "last_solved_exponent": precision,
                "attempted_exponent": precision + 1,
                "certificate": solve,
                "residues": point,
                "modulus": old_modulus,
                "scope": (
                    "[THEOREM] This selected residue has no next digit; "
                    "[UNRESOLVED] other singular paths remain possible."
                ),
            }
        generator = random.Random(50000 + 1000 * prime + precision)
        chosen_point = None
        attempts_used = 0
        free_count = len(solve["free_columns"])
        for attempt in range(trials_per_digit):
            free_values = (
                [0] * free_count
                if attempt == 0
                else [generator.randrange(prime) for _ in range(free_count)]
            )
            correction = instantiate_linear_solution(
                solve, free_values, prime, len(point)
            )
            candidate = [
                (value + old_modulus * digit) % new_modulus
                for value, digit in zip(point, correction, strict=True)
            ]
            if any(eval_system_mod(polynomials, candidate, new_modulus)):
                raise AssertionError("singular digit correction failed")
            attempts_used = attempt + 1
            if precision + 1 < exponent:
                lookahead_modulus = new_modulus * prime
                lookahead_residuals = eval_system_mod(
                    polynomials, candidate, lookahead_modulus
                )
                lookahead_rhs = [
                    -(value // new_modulus) % prime for value in lookahead_residuals
                ]
                if not solve_rectangular_mod(
                    jacobian, lookahead_rhs, prime
                )["solvable"]:
                    continue
            chosen_point = candidate
            break
        if chosen_point is None:
            return {
                "status": "lookahead_search_exhausted",
                "last_solved_exponent": precision,
                "attempted_exponent": precision + 1,
                "trials": trials_per_digit,
                "residues": point,
                "modulus": old_modulus,
                "scope": (
                    "[UNRESOLVED] Finite exact search failure on this selected "
                    "path is not a nonlift theorem."
                ),
            }
        point = chosen_point
        steps.append(
            {
                "from_exponent": precision,
                "to_exponent": precision + 1,
                "lookahead_trials": attempts_used,
                "free_dimension": free_count,
            }
        )
    homogeneous = solve_rectangular_mod(jacobian, [0] * len(jacobian), prime)
    return {
        "status": "lifted",
        "jacobian_rank_mod_p": homogeneous["rank"],
        "pivot_columns": homogeneous["pivots"],
        "precision_exponent": exponent,
        "modulus": prime**exponent,
        "residues": point,
        "steps": steps,
    }


def serializable_lift(record: dict) -> dict:
    converted = dict(record)
    if "modulus" in converted:
        converted["modulus"] = str(converted["modulus"])
    if "residues" in converted:
        converted["residues"] = [str(value) for value in converted["residues"]]
    return converted


def padic_valuation_mod(residue: int, prime: int, exponent: int) -> int:
    residue %= prime**exponent
    if residue == 0:
        return exponent
    valuation = 0
    while residue % prime == 0:
        residue //= prime
        valuation += 1
    return valuation


def holdout_modular_record(
    prime: int,
    lift: dict,
    holdout_polynomials: list[SparsePolynomial],
) -> dict:
    modulus = int(lift["modulus"])
    exponent = int(lift["precision_exponent"])
    residues = eval_system_mod(
        holdout_polynomials,
        [int(value) for value in lift["residues"]],
        modulus,
    )
    valuations = [padic_valuation_mod(value, prime, exponent) for value in residues]
    passed = not any(residues)
    return {
        "claim_tag": "[COMPUTATION]",
        "shape": [3, 3, 3],
        "orders": list(ORDERS),
        "modulus": str(modulus),
        "residues": [str(value) for value in residues],
        "p_adic_valuations_capped_at_precision": valuations,
        "passed_to_full_precision": passed,
        "conclusion": (
            "[UNRESOLVED] The selected construction lift passes this holdout only to finite precision."
            if passed
            else "[COMPUTATION] A nonzero exact modular residual falsifies this selected p-adic construction lift on the 3x3x3 identity."
        ),
    }


def exact_candidate_record(
    fractions: list[Fraction],
    variables: list[sp.Symbol],
    construction_expressions: list[sp.Expr],
) -> dict:
    substitutions = {
        variable: sp.Rational(value.numerator, value.denominator)
        for variable, value in zip(variables, fractions, strict=True)
    }
    construction_residuals = [
        sp.cancel(equation.subs(substitutions)) for equation in construction_expressions
    ]
    record = {
        "claim_tag": "[COMPUTATION]",
        "coordinates": [str(value) for value in fractions],
        "construction_residuals": [str(value) for value in construction_residuals],
        "construction_passed": not any(construction_residuals),
        "holdouts": [],
        "rational_point_found": False,
    }
    if any(construction_residuals):
        record["conclusion"] = (
            "[COMPUTATION] Exact substitution falsifies the reconstructed vector "
            "already on the construction equations."
        )
        return record
    for shape in HOLDOUT_SHAPES:
        system = full_weight_finite_system((shape,), ORDERS, gauge_fix=True)
        if [str(value) for value in system.variables] != [str(value) for value in variables]:
            raise AssertionError("holdout variable order changed")
        residuals = [sp.cancel(equation.subs(substitutions)) for equation in system.equations]
        record["holdouts"].append(
            {
                "shape": list(shape),
                "orders": list(ORDERS),
                "residuals": [str(value) for value in residuals],
                "passed": not any(residuals),
            }
        )
        if any(residuals):
            record["conclusion"] = (
                "[COMPUTATION] Exact rational evaluation falsifies the candidate "
                f"on the {shape} lattice identity."
            )
            return record
    record["rational_point_found"] = len(record["holdouts"]) == 2
    record["conclusion"] = (
        "[COMPUTATION] The rational vector passes the construction equations and two exact holdouts."
    )
    return record


def analyze_padic_witnesses(
    source: dict,
    variables: list[sp.Symbol],
    construction_expressions: list[sp.Expr],
    polynomials: list[SparsePolynomial],
) -> tuple[list[dict], list[dict], dict, list[SparsePolynomial]]:
    derivatives = derivative_table(polynomials, len(variables))
    holdout = full_weight_finite_system(((3, 3, 3),), ORDERS, gauge_fix=True)
    if [str(value) for value in holdout.variables] != [str(value) for value in variables]:
        raise AssertionError("3x3x3 variable order differs from construction system")
    holdout_polynomials = sparse_polynomials(list(holdout.equations), variables)
    witnesses = {
        int(row["prime"]): [int(value) for value in row["values_in_variable_order"]]
        for row in source["data"]["exact_modular_solution_witnesses"]
    }
    records = []
    exact_candidates = []
    internal_lifts: dict[int, dict] = {}
    for prime in sorted(witnesses):
        initial = witnesses[prime]
        construction_residues = eval_system_mod(polynomials, initial, prime)
        if any(construction_residues):
            raise AssertionError(f"stored F_{prime} point is no longer a construction solution")
        rank, pivots, _ = rref_mod(jacobian_mod(derivatives, initial, prime), prime)
        record = {
            "claim_tag": "[COMPUTATION]",
            "branch": "11111",
            "prime": prime,
            "stored_solution_verified": True,
            "stored_construction_residues": construction_residues,
            "jacobian_rank_mod_p": rank,
            "initial_pivot_columns": pivots,
        }
        requested_exponent = max(
            MINIMUM_PADIC_EXPONENT,
            target_exponent(prime, MINIMUM_DECIMAL_DIGITS),
        )
        if prime == 3:
            lift = singular_hensel_section(
                prime,
                initial,
                polynomials,
                derivatives,
                requested_exponent,
            )
            record["local_section"] = {
                "method": "singular one-digit lift with deterministic exact one-digit lookahead",
                "random_seed_rule": "50000 + 1000*prime + current_exponent",
                "trials_per_digit": SINGULAR_TRIALS_PER_DIGIT,
                "scope": "[UNRESOLVED] One selected singular path is not the full F_3 residue class.",
            }
        elif rank == len(polynomials):
            lift = hensel_lift_full_row_rank(
                prime,
                initial,
                polynomials,
                derivatives,
                requested_exponent,
            )
            record["local_section"] = {
                "method": "full-row-rank Newton/Hensel section",
                "newton_variables": [str(variables[index]) for index in lift["pivot_columns"]],
                "held_coordinates": [
                    str(variables[index])
                    for index in range(len(variables))
                    if index not in lift["pivot_columns"]
                ],
                "held_coordinate_rule": "least nonnegative F_p representative held exactly",
                "scope": "[UNRESOLVED] This selected section is not the full residue class.",
            }
        else:
            lift = rank_deficient_first_lift_attempt(
                prime, initial, polynomials, derivatives
            )
        record.update(serializable_lift(lift))
        if lift["status"] == "lifted":
            internal_lifts[prime] = lift
            exponent = int(lift["precision_exponent"])
            modulus = int(lift["modulus"])
            record["decimal_digits_floor"] = len(str(modulus)) - 1
            reconstruction = reconstruct_lift(lift, [str(variable) for variable in variables])
            fractions = reconstruction.pop("fractions")
            reconstruction["failed_variables"] = [
                row["variable"]
                for row in reconstruction["coordinates"]
                if row["status"] == "failure"
            ]
            reconstruction["method"] = (
                "exact extended-Euclidean/continued-fraction rational reconstruction "
                "with symmetric numerator and denominator bound"
            )
            reconstruction["claim"] = (
                "[COMPUTATION] Every coordinate was attempted; each failure excludes a "
                "bounded fraction for this residue, not other p-adic sections."
            )
            record["rational_reconstruction"] = reconstruction
            record["holdout_3x3x3"] = holdout_modular_record(
                prime, lift, holdout_polynomials
            )
            if reconstruction["all_coordinates_reconstructed"]:
                if any(value is None for value in fractions):
                    raise AssertionError("all-coordinate reconstruction flag is inconsistent")
                candidate = exact_candidate_record(
                    [value for value in fractions if value is not None],
                    variables,
                    construction_expressions,
                )
                candidate["prime"] = prime
                exact_candidates.append(candidate)
            if exponent < MINIMUM_PADIC_EXPONENT:
                raise AssertionError("lift did not reach the required p-adic exponent")
        records.append(record)

    p3_initial = witnesses[3]
    combined_expressions = construction_expressions + list(holdout.equations)
    combined_polynomials = sparse_polynomials(combined_expressions, variables)
    combined_derivatives = derivative_table(combined_polynomials, len(variables))
    combined_attempt = rank_deficient_first_lift_attempt(
        3, p3_initial, combined_polynomials, combined_derivatives
    )
    combined_attempt.update(
        {
            "claim_tag": "[THEOREM]",
            "branch": "11111",
            "prime": 3,
            "system": "six construction equations plus 3x3x3 k=4,6,8",
            "equation_count": len(combined_polynomials),
            "conclusion": (
                "[THEOREM] This stored F_3 point has no simultaneous lift modulo 9 "
                "after the 3x3x3 identities are imposed; other F_3 points remain unresolved."
            ),
        }
    )
    return records, exact_candidates, combined_attempt, holdout_polynomials


def main() -> None:
    source, variables, construction_expressions, construction_polynomials = load_stored_system()
    source_rows = source["data"]["zero_pattern_branches"]
    source_by_branch = {row["branch"]: row for row in source_rows}
    labels = all_branch_labels()
    gauge_names = [pair_name(pair) for pair in DIRECTION_GAUGE_TREE]

    construction_ungauged = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
    )
    certificate_ungauged = full_weight_finite_system(
        CERTIFICATE_SHAPES, ORDERS, gauge_fix=False
    )
    certificate_labels = [
        {"shape": list(shape), "order": order}
        for shape, order in certificate_ungauged.labels
    ]
    remaining_symbols = [
        symbol
        for pair, symbol in construction_ungauged.all_symbols
        if pair not in DIRECTION_GAUGE_TREE
    ]

    branch_records = []
    constant_identities_verified = True
    metadata_matches = True
    for branch in labels:
        active, construction_equations = branch_equations(construction_ungauged, branch)
        globally_linear = [
            variable
            for variable in remaining_symbols
            if max(sp.degree(equation, variable) for equation in construction_equations) <= 1
            and any(sp.degree(equation, variable) == 1 for equation in construction_equations)
        ]
        certificate_variables, certificate_equations = branch_equations(
            certificate_ungauged, branch
        )
        search = constant_certificate_search(
            certificate_equations,
            certificate_variables,
            certificate_labels,
        )
        source_row = source_by_branch[branch]
        zero_pattern = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "0"
        ]
        normalized_nonzero = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "1"
        ]
        row_matches = (
            source_row["zero_pattern"] == zero_pattern
            and source_row["normalized_nonzero"] == normalized_nonzero
            and source_row["active_variable_count_in_six_equations"] == len(active)
            and source_row["globally_linear_variable_count"] == len(globally_linear)
        )
        metadata_matches &= row_matches
        record = {
            "branch": branch,
            "bits_in_gauge_tree_order": [int(bit) for bit in branch],
            "zero_pattern": zero_pattern,
            "normalized_nonzero": normalized_nonzero,
            "active_variable_count_in_six_equations": len(active),
            "globally_linear_variable_count": len(globally_linear),
            "wave4_metadata_recomputed_exactly": row_matches,
            "constant_certificate_search": search,
            "finite_field_primes_with_stored_points": (
                [3, 5, 7, 11] if branch == "11111" else []
            ),
        }
        if search["certificate_found"]:
            certificate = search["certificate"]
            record.update(
                {
                    "status": "EMPTY_OVER_Q",
                    "claim_tag": "[THEOREM]",
                    "claim": (
                        f"[THEOREM] Branch {branch} is empty over Q by the stored "
                        "nonzero integer constant identity."
                    ),
                    "empty_certificate": certificate,
                    "groebner": {
                        "status": "not_needed",
                        "reason": "exact constant identity already decides the branch",
                    },
                }
            )
            coefficients = certificate["coefficient_vector"]
            combination = sp.expand(
                sum(
                    coefficient * equation
                    for coefficient, equation in zip(
                        coefficients, certificate_equations, strict=True
                    )
                )
            )
            constant_identities_verified &= (
                not combination.free_symbols
                and combination == certificate["constant"]
                and certificate["constant"] != 0
            )
        else:
            record["status"] = "PENDING_GROEBNER"
        branch_records.append(record)

    for record in branch_records:
        if record["status"] != "PENDING_GROEBNER":
            continue
        groebner = bounded_groebner_attempt(record["branch"])
        if groebner["status"] == "error":
            raise RuntimeError(
                f"Groebner worker failed on {record['branch']}: {groebner.get('error')}"
            )
        record["groebner"] = groebner
        if groebner.get("status") == "completed" and groebner.get("contains_one"):
            record.update(
                {
                    "status": "EMPTY_OVER_Q",
                    "claim_tag": "[THEOREM]",
                    "claim": (
                        f"[THEOREM] Branch {record['branch']} is empty over Q because "
                        "its independently recorded exact Groebner basis is [1]."
                    ),
                    "empty_certificate": {
                        "claim_tag": "[THEOREM]",
                        "type": "exact_Q_groebner_basis_one",
                        "basis_sha256": groebner["basis_sha256"],
                        "basis": groebner["basis"],
                    },
                }
            )
        else:
            record.update(
                {
                    "status": "UNRESOLVED",
                    "claim_tag": "[UNRESOLVED]",
                    "claim": (
                        f"[UNRESOLVED] Branch {record['branch']} has no exact emptiness "
                        "certificate and no verified rational point."
                    ),
                }
            )

    padic_records, exact_candidates, p3_combined, holdout_polynomials = (
        analyze_padic_witnesses(
            source,
            variables,
            construction_expressions,
            construction_polynomials,
        )
    )
    breakthrough = any(row["rational_point_found"] for row in exact_candidates)
    generic_branch = next(row for row in branch_records if row["branch"] == "11111")
    generic_branch["p_adic_evidence"] = {
        str(row["prime"]): {
            "status": row["status"],
            "precision_exponent": row.get("precision_exponent"),
            "height_bound": row.get("rational_reconstruction", {}).get("height_bound"),
            "failed_variables": row.get("rational_reconstruction", {}).get(
                "failed_variables", []
            ),
            "holdout_3x3x3_passed": row.get("holdout_3x3x3", {}).get(
                "passed_to_full_precision"
            ),
        }
        for row in padic_records
    }
    generic_branch["p3_combined_holdout_nonlift"] = {
        "status": p3_combined["status"],
        "attempted_modulus": 9,
    }
    if breakthrough:
        generic_branch.update(
            {
                "status": "RATIONAL_POINT_FOUND",
                "claim_tag": "[COMPUTATION]",
                "claim": (
                    "[COMPUTATION] A rational point passed the construction system and "
                    "two exact holdout lattice identities."
                ),
            }
        )

    status_counts = {
        status: sum(row["status"] == status for row in branch_records)
        for status in ("EMPTY_OVER_Q", "RATIONAL_POINT_FOUND", "UNRESOLVED")
    }
    all_empty = status_counts["EMPTY_OVER_Q"] == len(branch_records)
    if all_empty:
        headline = (
            "[THEOREM] Every one of the 32 gauge-tree branches is empty over Q; "
            "the selected 25-coordinate k<=8 family has no rational realization."
        )
    else:
        headline = (
            f"[COMPUTATION] Exact partial branch map: {status_counts['EMPTY_OVER_Q']} "
            f"of {len(branch_records)} branches are EMPTY_OVER_Q; "
            f"[UNRESOLVED] {status_counts['UNRESOLVED']} branches remain."
        )

    lifted = [row for row in padic_records if row["status"] == "lifted"]
    padic_residuals_ok = all(
        not any(
            eval_system_mod(
                construction_polynomials,
                [int(value) for value in row["residues"]],
                int(row["modulus"]),
            )
        )
        for row in lifted
    )
    reconstruction_ok = True
    for row in lifted:
        modulus = int(row["modulus"])
        bound = int(row["rational_reconstruction"]["height_bound"])
        reconstruction_ok &= 2 * bound * bound < modulus
        reconstruction_ok &= 2 * (bound + 1) * (bound + 1) >= modulus
        for coordinate, residue in zip(
            row["rational_reconstruction"]["coordinates"],
            row["residues"],
            strict=True,
        ):
            observed = rational_reconstruct(int(residue), modulus, bound)
            if coordinate["status"] == "failure":
                reconstruction_ok &= observed is None
            else:
                reconstruction_ok &= observed == Fraction(
                    int(coordinate["numerator"]), int(coordinate["denominator"])
                )

    p7 = next(row for row in padic_records if row["prime"] == 7)
    p7_certificate = p7.get("certificate", {})
    p3_certificate = p3_combined.get("certificate", {})
    checks = [
        {
            "name": "all_32_binary_branches_reconstructed_in_wave4_order",
            "passed": (
                len(labels) == 32
                and labels == [row["branch"] for row in source_rows]
                and labels == [row["branch"] for row in branch_records]
                and metadata_matches
            ),
            "detail": {
                "total": len(branch_records),
                "first": labels[0],
                "last": labels[-1],
                "gauge_tree": gauge_names,
            },
        },
        {
            "name": "stored_constant_certificates_recomputed_exactly",
            "passed": constant_identities_verified,
            "detail": {
                "certificate_count": sum(
                    row["constant_certificate_search"]["certificate_found"]
                    for row in branch_records
                ),
                "constraint_count": len(certificate_ungauged.equations),
            },
        },
        {
            "name": "seventeen_branches_have_integer_constant_certificates",
            "passed": sum(
                row["constant_certificate_search"]["certificate_found"]
                for row in branch_records
            )
            == 17,
            "detail": [
                row["branch"]
                for row in branch_records
                if row["constant_certificate_search"]["certificate_found"]
            ],
        },
        {
            "name": "every_nonconstant_branch_received_bounded_exact_Q_groebner_attempt",
            "passed": all(
                row["groebner"]["status"] in {"completed", "timeout"}
                and row["groebner"]["limit_seconds"] <= 900
                and (
                    row["groebner"]["status"] != "completed"
                    or bool(row["groebner"].get("basis"))
                )
                for row in branch_records
                if not row["constant_certificate_search"]["certificate_found"]
            ),
            "detail": {
                row["branch"]: row["groebner"]["status"]
                for row in branch_records
                if not row["constant_certificate_search"]["certificate_found"]
            },
        },
        {
            "name": "per_branch_status_ledger_is_exhaustive",
            "passed": (
                sum(status_counts.values()) == 32
                and all(
                    row["status"]
                    in {"EMPTY_OVER_Q", "RATIONAL_POINT_FOUND", "UNRESOLVED"}
                    for row in branch_records
                )
            ),
            "detail": status_counts,
        },
        {
            "name": "stored_F3_F5_F7_F11_points_reverified",
            "passed": all(row["stored_solution_verified"] for row in padic_records),
            "detail": {
                str(row["prime"]): row["stored_construction_residues"]
                for row in padic_records
            },
        },
        {
            "name": "all_successful_lifts_reach_at_least_p_power_60",
            "passed": (
                {row["prime"] for row in lifted} == {3, 5, 11}
                and all(row["precision_exponent"] >= MINIMUM_PADIC_EXPONENT for row in lifted)
                and padic_residuals_ok
            ),
            "detail": {
                str(row["prime"]): {
                    "exponent": row["precision_exponent"],
                    "decimal_digits_floor": row["decimal_digits_floor"],
                }
                for row in lifted
            },
        },
        {
            "name": "every_lifted_coordinate_received_exact_rational_reconstruction",
            "passed": reconstruction_ok
            and all(
                len(row["rational_reconstruction"]["coordinates"]) == len(variables)
                for row in lifted
            ),
            "detail": {
                str(row["prime"]): {
                    "successes": row["rational_reconstruction"][
                        "successful_coordinate_count"
                    ],
                    "failures": row["rational_reconstruction"]["failed_variables"],
                    "height_bound": row["rational_reconstruction"]["height_bound"],
                }
                for row in lifted
            },
        },
        {
            "name": "stored_F7_point_nonlift_mod_49_certificate",
            "passed": (
                p7["status"] == "no_lift_mod_p_squared"
                and not any(p7_certificate.get("left_kernel_times_jacobian_mod_p", [1]))
                and p7_certificate.get("pairing_with_rhs_mod_p", 0) != 0
            ),
            "detail": p7_certificate,
        },
        {
            "name": "stored_F3_point_combined_3x3x3_nonlift_mod_9_certificate",
            "passed": (
                p3_combined["status"] == "no_lift_mod_p_squared"
                and not any(p3_certificate.get("left_kernel_times_jacobian_mod_p", [1]))
                and p3_certificate.get("pairing_with_rhs_mod_p", 0) != 0
            ),
            "detail": p3_certificate,
        },
        {
            "name": "selected_padic_holdout_falsifications_recorded_exactly",
            "passed": all(
                not row["holdout_3x3x3"]["passed_to_full_precision"] for row in lifted
            ),
            "detail": {
                str(row["prime"]): row["holdout_3x3x3"][
                    "p_adic_valuations_capped_at_precision"
                ]
                for row in lifted
            },
        },
        {
            "name": "rational_point_claim_requires_two_exact_holdouts",
            "passed": (
                not breakthrough
                or any(
                    row["construction_passed"]
                    and len(row["holdouts"]) == 2
                    and all(holdout_row["passed"] for holdout_row in row["holdouts"])
                    for row in exact_candidates
                )
            ),
            "detail": {
                "candidate_count": len(exact_candidates),
                "breakthrough": breakthrough,
            },
        },
        {
            "name": "global_theorem_withheld_unless_all_branches_empty",
            "passed": all_empty or headline.startswith("[COMPUTATION]"),
            "detail": {"all_empty": all_empty, "headline": headline},
        },
    ]

    output = {
        "provenance": {
            "script": "experiments/e50_kw_branches.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact SymPy ZZ/QQ branch substitution and coefficient nullspaces; "
                "hard-timeout exact QQ Groebner; exact integer p-adic lifting and "
                "continued-fraction rational reconstruction; no float64 decisions"
            ),
            "source_artifact": "results/kac_ward/full_family.json",
            "source_system_sha256": source["data"]["minimal_unresolved_system"][
                "sha256_sympy_srepr"
            ],
            "groebner_timeout_seconds_per_branch": GROEBNER_TIMEOUT_SECONDS,
            "minimum_padic_exponent": MINIMUM_PADIC_EXPONENT,
            "minimum_decimal_digits_when_liftable": MINIMUM_DECIMAL_DIGITS,
        },
        "data": {
            "headline": headline,
            "total_branch_count": len(branch_records),
            "status_counts": status_counts,
            "gauge_tree_order": gauge_names,
            "bit_convention": "0=fixed tree weight vanishes; 1=normalized to one",
            "constraint_families": {
                "construction_shapes": [list(shape) for shape in CONSTRUCTION_SHAPES],
                "constant_certificate_shapes": [
                    list(shape) for shape in CERTIFICATE_SHAPES
                ],
                "orders": list(ORDERS),
                "constant_certificate_equation_count": len(
                    certificate_ungauged.equations
                ),
            },
            "branches": branch_records,
            "p_adic_records_for_branch_11111": padic_records,
            "p3_combined_3x3x3_nonlift": p3_combined,
            "exact_rational_reconstruction_candidates": exact_candidates,
            "global_claim": (
                "[THEOREM] All enumerated branches are empty over Q."
                if all_empty
                else "[UNRESOLVED] No global no-go is claimed because at least one branch remains unresolved."
            ),
            "scope": (
                "[UNRESOLVED] EMPTY_OVER_Q is decisive for the named gauge-tree branch. "
                "Finite-field points, selected p-adic lifts, reconstruction bounds, and "
                "Groebner timeouts do not decide any other rational point or support chart."
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
