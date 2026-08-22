"""p-adic lifting diagnostics for the stored 25-variable Kac--Ward system.

The construction variety is underdetermined (six equations in 25 variables).
At a full-row-rank finite-field point we therefore choose a reproducible local
section: the first six pivot coordinates are Newton variables and the other
19 coordinates are held at their least nonnegative representatives.  Failure
of rational reconstruction concerns this selected p-adic section only.
"""

from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import queue
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy.polys.orderings import ProductOrder, grevlex, lex

from ising.fermions.kac_ward import DIRECTION_GAUGE_TREE, full_weight_finite_system

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/kac_ward/full_family.json"
OUTPUT = ROOT / "results/kac_ward/hensel_lift.json"
TARGET_PRIMES = (5, 7, 11)
GROEBNER_PRIMES = (101, 32003)
GROEBNER_TIMEOUT_SECONDS = 2400
MIN_DECIMAL_DIGITS = 60
CONSTRUCTION_SHAPES = ((3, 3, 2), (2, 2, 3))
ORDERS = (4, 6, 8)

SparsePolynomial = tuple[tuple[int, tuple[int, ...]], ...]


def _canonical_system_hash(
    labels: list[dict], equations: Iterable[sp.Expr]
) -> str:
    canonical = "\n".join(
        f"{tuple(label['shape'])}:{label['order']}:{sp.srepr(equation)}"
        for label, equation in zip(labels, equations, strict=True)
    )
    return hashlib.sha256(canonical.encode()).hexdigest()


def load_stored_system() -> tuple[dict, list[sp.Symbol], list[sp.Expr], list[SparsePolynomial]]:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    stored = payload["data"]["minimal_unresolved_system"]
    variables = list(sp.symbols(" ".join(stored["variables"])))
    namespace = dict(zip(stored["variables"], variables, strict=True))
    expressions = [sp.expand(sp.sympify(text, locals=namespace)) for text in stored["equations"]]
    observed_hash = _canonical_system_hash(stored["labels"], expressions)
    if observed_hash != stored["sha256_sympy_srepr"]:
        raise AssertionError(
            f"stored system hash mismatch: {observed_hash} != {stored['sha256_sympy_srepr']}"
        )
    sparse = []
    for expression in expressions:
        polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
        sparse.append(
            tuple(
                (int(coefficient), tuple(int(exponent) for exponent in monomial))
                for monomial, coefficient in polynomial.terms()
            )
        )
    return payload, variables, expressions, sparse


def eval_sparse_mod(poly: SparsePolynomial, point: list[int], modulus: int) -> int:
    total = 0
    for coefficient, exponents in poly:
        term = coefficient % modulus
        for value, exponent in zip(point, exponents, strict=True):
            if exponent:
                term = (term * pow(value, exponent, modulus)) % modulus
        total = (total + term) % modulus
    return total


def eval_system_mod(
    polynomials: list[SparsePolynomial], point: list[int], modulus: int
) -> list[int]:
    return [eval_sparse_mod(poly, point, modulus) for poly in polynomials]


def derivative_sparse(poly: SparsePolynomial, coordinate: int) -> SparsePolynomial:
    derivative = []
    for coefficient, exponents in poly:
        exponent = exponents[coordinate]
        if not exponent:
            continue
        reduced = list(exponents)
        reduced[coordinate] -= 1
        derivative.append((coefficient * exponent, tuple(reduced)))
    return tuple(derivative)


def derivative_table(
    polynomials: list[SparsePolynomial], variable_count: int
) -> list[list[SparsePolynomial]]:
    return [
        [derivative_sparse(poly, coordinate) for coordinate in range(variable_count)]
        for poly in polynomials
    ]


def jacobian_mod(
    derivatives: list[list[SparsePolynomial]], point: list[int], modulus: int
) -> list[list[int]]:
    return [
        [eval_sparse_mod(derivative, point, modulus) for derivative in row]
        for row in derivatives
    ]


def rref_mod(matrix: list[list[int]], prime: int) -> tuple[int, list[int], list[list[int]]]:
    reduced = [[value % prime for value in row] for row in matrix]
    rank = 0
    pivots: list[int] = []
    column_count = len(reduced[0]) if reduced else 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, len(reduced)) if reduced[row][column]), None
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        inverse = pow(reduced[rank][column], -1, prime)
        reduced[rank] = [(value * inverse) % prime for value in reduced[rank]]
        for row in range(len(reduced)):
            if row == rank or not reduced[row][column]:
                continue
            factor = reduced[row][column]
            reduced[row] = [
                (left - factor * right) % prime
                for left, right in zip(reduced[row], reduced[rank], strict=True)
            ]
        pivots.append(column)
        rank += 1
        if rank == len(reduced):
            break
    return rank, pivots, reduced


def solve_square_mod(
    matrix: list[list[int]], right_hand_side: list[int], modulus: int
) -> list[int]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("a nonempty square matrix is required")
    augmented = [
        [value % modulus for value in row] + [right_hand_side[index] % modulus]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if math.gcd(augmented[row][column], modulus) == 1
            ),
            None,
        )
        if pivot is None:
            raise ArithmeticError("matrix is not invertible over Z/modulus")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, modulus)
        augmented[column] = [
            (value * inverse) % modulus for value in augmented[column]
        ]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (left - factor * right) % modulus
                for left, right in zip(augmented[row], augmented[column], strict=True)
            ]
    return [augmented[index][-1] for index in range(size)]


def target_exponent(prime: int, minimum_decimal_digits: int = MIN_DECIMAL_DIGITS) -> int:
    exponent = 1
    while len(str(prime**exponent)) - 1 < minimum_decimal_digits:
        exponent += 1
    return exponent


def reconstruction_height_bound(modulus: int) -> int:
    """Largest symmetric H with the rational-reconstruction uniqueness 2 H^2 < M."""
    if modulus <= 2:
        raise ValueError("modulus must exceed two")
    return math.isqrt((modulus - 1) // 2)


def rational_reconstruct(residue: int, modulus: int, bound: int) -> Fraction | None:
    """Return the unique bounded rational representative, if it exists.

    The caller must impose ``2*bound**2 < modulus``.  This is the classical
    extended-Euclidean rational reconstruction with symmetric numerator and
    positive denominator bounds.
    """
    if bound < 1 or 2 * bound * bound >= modulus:
        raise ValueError("bound must be positive and satisfy 2*bound^2 < modulus")
    residue %= modulus
    if residue == 0:
        return Fraction(0, 1)
    old_r, r = modulus, residue
    old_t, t = 0, 1
    while r > bound:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_t, t = t, old_t - quotient * t
    numerator, denominator = r, t
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        abs(numerator) > bound
        or denominator < 1
        or denominator > bound
        or math.gcd(numerator, denominator) != 1
        or math.gcd(denominator, modulus) != 1
        or (numerator - residue * denominator) % modulus
    ):
        return None
    return Fraction(numerator, denominator)


def hensel_lift_full_row_rank(
    prime: int,
    initial: list[int],
    polynomials: list[SparsePolynomial],
    derivatives: list[list[SparsePolynomial]],
    exponent: int,
) -> dict:
    jacobian_at_base = jacobian_mod(derivatives, initial, prime)
    rank, pivot_columns, _ = rref_mod(jacobian_at_base, prime)
    equation_count = len(polynomials)
    if rank != equation_count:
        raise ValueError(f"Jacobian row rank {rank}, expected {equation_count}")
    pivot_columns = pivot_columns[:equation_count]
    point = [value % prime for value in initial]
    precision = 1
    steps = [{"from_exponent": 0, "to_exponent": 1, "residuals_zero": True}]
    while precision < exponent:
        new_precision = min(2 * precision, exponent)
        old_modulus = prime**precision
        correction_modulus = prime ** (new_precision - precision)
        new_modulus = prime**new_precision
        residuals = eval_system_mod(polynomials, point, new_modulus)
        if any(residual % old_modulus for residual in residuals):
            raise AssertionError("input ceased to solve at the current Hensel precision")
        right_hand_side = [
            (-(residual // old_modulus)) % correction_modulus for residual in residuals
        ]
        full_jacobian = jacobian_mod(derivatives, point, correction_modulus)
        square_jacobian = [
            [row[column] for column in pivot_columns] for row in full_jacobian
        ]
        correction = solve_square_mod(
            square_jacobian, right_hand_side, correction_modulus
        )
        for column, value in zip(pivot_columns, correction, strict=True):
            point[column] = (point[column] + old_modulus * value) % new_modulus
        point = [value % new_modulus for value in point]
        lifted_residuals = eval_system_mod(polynomials, point, new_modulus)
        if any(lifted_residuals):
            raise AssertionError("Newton/Hensel correction did not annihilate residuals")
        steps.append(
            {
                "from_exponent": precision,
                "to_exponent": new_precision,
                "correction_modulus": str(correction_modulus),
                "residuals_zero": True,
            }
        )
        precision = new_precision
    modulus = prime**exponent
    return {
        "status": "lifted",
        "jacobian_rank_mod_p": rank,
        "pivot_columns": pivot_columns,
        "precision_exponent": exponent,
        "modulus": modulus,
        "residues": point,
        "steps": steps,
    }


def rank_deficient_first_lift_attempt(
    prime: int,
    initial: list[int],
    polynomials: list[SparsePolynomial],
    derivatives: list[list[SparsePolynomial]],
) -> dict:
    """Attempt x+p*t modulo p^2 and return an exact inconsistency certificate."""
    jacobian = jacobian_mod(derivatives, initial, prime)
    rank, pivot_columns, _ = rref_mod(jacobian, prime)
    residuals_mod_p2 = eval_system_mod(polynomials, initial, prime * prime)
    if any(residual % prime for residual in residuals_mod_p2):
        raise AssertionError("stored point is not a solution modulo p")
    right_hand_side = [
        (-(residual // prime)) % prime for residual in residuals_mod_p2
    ]

    equation_count = len(jacobian)
    variable_count = len(jacobian[0])
    augmented = []
    for row_index, (row, rhs) in enumerate(zip(jacobian, right_hand_side, strict=True)):
        augmented.append(
            [value % prime for value in row]
            + [rhs % prime]
            + [1 if row_index == index else 0 for index in range(equation_count)]
        )
    reduced_row = 0
    for column in range(variable_count):
        pivot = next(
            (
                row
                for row in range(reduced_row, equation_count)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        augmented[reduced_row], augmented[pivot] = (
            augmented[pivot],
            augmented[reduced_row],
        )
        inverse = pow(augmented[reduced_row][column], -1, prime)
        augmented[reduced_row] = [
            (value * inverse) % prime for value in augmented[reduced_row]
        ]
        for row in range(equation_count):
            if row == reduced_row or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (left - factor * right) % prime
                for left, right in zip(
                    augmented[row], augmented[reduced_row], strict=True
                )
            ]
        reduced_row += 1
        if reduced_row == equation_count:
            break

    certificate = None
    for row in augmented:
        if (
            all(row[column] % prime == 0 for column in range(variable_count))
            and row[variable_count] % prime
        ):
            certificate = {
                "left_kernel_vector": row[variable_count + 1 :],
                "pairing_with_rhs_mod_p": row[variable_count] % prime,
            }
            break
    if certificate is None:
        return {
            "status": "first_digit_lift_solvable_but_not_continued",
            "jacobian_rank_mod_p": rank,
            "pivot_columns": pivot_columns,
            "residuals_mod_p2": residuals_mod_p2,
            "linearized_rhs_mod_p": right_hand_side,
            "scope": "Only one deterministic linearized correction was examined.",
        }

    left_kernel = certificate["left_kernel_vector"]
    kernel_products = [
        sum(left_kernel[row] * jacobian[row][column] for row in range(equation_count))
        % prime
        for column in range(variable_count)
    ]
    rhs_pairing = (
        sum(left_kernel[row] * right_hand_side[row] for row in range(equation_count))
        % prime
    )
    if any(kernel_products) or rhs_pairing == 0:
        raise AssertionError("invalid first-order nonlift certificate")
    certificate["left_kernel_times_jacobian_mod_p"] = kernel_products
    return {
        "status": "no_lift_mod_p_squared",
        "jacobian_rank_mod_p": rank,
        "pivot_columns": pivot_columns,
        "residuals_mod_p2": residuals_mod_p2,
        "linearized_rhs_mod_p": right_hand_side,
        "certificate": certificate,
        "scope": (
            "This proves that this particular stored F_p point has no lift to "
            "a solution modulo p^2; it says nothing about other F_p points."
        ),
    }


def reconstruct_lift(record: dict, variable_names: list[str]) -> dict:
    modulus = int(record["modulus"])
    bound = reconstruction_height_bound(modulus)
    coordinates = []
    fractions: list[Fraction | None] = []
    for name, residue in zip(variable_names, record["residues"], strict=True):
        candidate = rational_reconstruct(int(residue), modulus, bound)
        fractions.append(candidate)
        if candidate is None:
            coordinates.append({"variable": name, "status": "failure"})
        else:
            coordinates.append(
                {
                    "variable": name,
                    "status": "success",
                    "numerator": str(candidate.numerator),
                    "denominator": str(candidate.denominator),
                }
            )
    return {
        "height_bound": str(bound),
        "uniqueness_inequality": f"2*H^2 < {modulus}",
        "successful_coordinate_count": sum(value is not None for value in fractions),
        "all_coordinates_reconstructed": all(value is not None for value in fractions),
        "coordinates": coordinates,
        "fractions": fractions,
    }


def _exact_candidate_checks(
    fractions: list[Fraction],
    expressions: list[sp.Expr],
    variables: list[sp.Symbol],
) -> dict:
    substitutions = {
        variable: sp.Rational(value.numerator, value.denominator)
        for variable, value in zip(variables, fractions, strict=True)
    }
    construction_residuals = [sp.cancel(expression.subs(substitutions)) for expression in expressions]
    construction_passed = not any(construction_residuals)
    record = {
        "construction_residuals": [str(value) for value in construction_residuals],
        "construction_passed": construction_passed,
        "holdouts": [],
        "breakthrough_candidate": False,
    }
    if not construction_passed:
        return record

    for shape in ((3, 3, 3), (4, 2, 3)):
        holdout = full_weight_finite_system((shape,), ORDERS, gauge_fix=True)
        if [str(variable) for variable in holdout.variables] != [str(v) for v in variables]:
            raise AssertionError("holdout variable order differs from stored system")
        residuals = [sp.cancel(equation.subs(substitutions)) for equation in holdout.equations]
        record["holdouts"].append(
            {
                "shape": list(shape),
                "orders": list(ORDERS),
                "residuals": [str(value) for value in residuals],
                "passed": not any(residuals),
            }
        )
        if any(residuals):
            break
    record["breakthrough_candidate"] = (
        len(record["holdouts"]) == 2
        and all(holdout["passed"] for holdout in record["holdouts"])
    )
    return record


def _smallest_branch_data() -> tuple[str, list[str]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    branches = source["data"]["zero_pattern_branches"]
    minimum = min(row["active_variable_count_in_six_equations"] for row in branches)
    ties = sorted(
        row["branch"]
        for row in branches
        if row["active_variable_count_in_six_equations"] == minimum
    )
    return ties[0], ties


def exact_branch_constant_certificate(branch: str) -> dict:
    """Verify the short integer Nullstellensatz certificate found by elimination."""
    ungauged = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
    )
    symbol_map = dict(ungauged.all_symbols)
    substitutions = {
        symbol_map[pair]: int(bit)
        for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
    }
    equations = [
        sp.expand(equation.subs(substitutions)) for equation in ungauged.equations
    ]
    coefficients = [2, 0, 0, -3, 0, 0]
    combination = sp.expand(
        sum(
            coefficient * equation
            for coefficient, equation in zip(coefficients, equations, strict=True)
        )
    )
    constant = int(combination) if not combination.free_symbols else None
    return {
        "branch": branch,
        "equation_labels": [
            {"shape": list(shape), "order": order}
            for shape, order in ungauged.labels
        ],
        "combination_coefficients": coefficients,
        "constant": constant,
        "identity_verified": constant == 56,
        "identity": "2*F_(3x3x2,4) - 3*F_(2x2x3,4) = 56",
        "consequence": (
            "The branch has no common zero over Q or any field of characteristic "
            "not dividing 56. This exact integer certificate, unlike modular [1] "
            "alone, proves characteristic-zero emptiness of this one branch."
        ),
        "scope": (
            "Only tree-zero pattern 00000 is excluded. The other 31 recorded tree "
            "patterns and support charts outside this tree stratification are not decided."
        ),
    }


def _groebner_worker(prime: int, branch: str, result_queue) -> None:
    started = time.monotonic()
    try:
        ungauged = full_weight_finite_system(
            CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
        )
        symbol_map = dict(ungauged.all_symbols)
        bits = tuple(int(bit) for bit in branch)
        substitutions = {
            symbol_map[pair]: value
            for pair, value in zip(DIRECTION_GAUGE_TREE, bits, strict=True)
        }
        equations = [sp.expand(equation.subs(substitutions)) for equation in ungauged.equations]
        equations = [equation for equation in equations if equation != 0]
        active = sorted(
            set().union(*(equation.free_symbols for equation in equations)), key=str
        )
        linear = [
            variable
            for variable in active
            if max(sp.degree(equation, variable) for equation in equations) <= 1
            and any(sp.degree(equation, variable) == 1 for equation in equations)
        ]
        nonlinear = [variable for variable in active if variable not in linear]
        ordered_variables = linear + nonlinear
        split = len(linear)
        order = ProductOrder(
            (lex, lambda monomial: monomial[:split]),
            (grevlex, lambda monomial: monomial[split:]),
        )
        basis = sp.groebner(
            equations,
            *ordered_variables,
            modulus=prime,
            order=order,
            method="f5b",
        )
        basis_strings = [str(poly.as_expr()) for poly in basis.polys]
        canonical = "\n".join(basis_strings)
        result_queue.put(
            {
                "status": "completed",
                "elapsed_seconds": time.monotonic() - started,
                "active_variables": [str(variable) for variable in active],
                "linear_elimination_block": [str(variable) for variable in linear],
                "remaining_block": [str(variable) for variable in nonlinear],
                "order": "ProductOrder(lex(linear block), grevlex(remaining block))",
                "basis_size": len(basis.polys),
                "basis_term_counts": [len(poly.terms()) for poly in basis.polys],
                "basis_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
                "contains_one": basis_strings == ["1"],
            }
        )
    except BaseException as error:  # returned to the parent as an honest failed attempt
        result_queue.put(
            {
                "status": "error",
                "elapsed_seconds": time.monotonic() - started,
                "error_type": type(error).__name__,
                "error": str(error),
            }
        )


def bounded_groebner_attempt(prime: int, branch: str) -> dict:
    context = mp.get_context("spawn")
    result_queue = context.Queue()
    process = context.Process(target=_groebner_worker, args=(prime, branch, result_queue))
    started = time.monotonic()
    process.start()
    process.join(GROEBNER_TIMEOUT_SECONDS)
    if process.is_alive():
        process.terminate()
        process.join(10)
        return {
            "prime": prime,
            "branch": branch,
            "status": "timeout",
            "limit_seconds": GROEBNER_TIMEOUT_SECONDS,
            "elapsed_seconds": time.monotonic() - started,
            "scope": "A timeout is non-decisive and is not evidence of nonemptiness.",
        }
    try:
        record = result_queue.get(timeout=5)
    except queue.Empty:
        record = {
            "status": "error",
            "error_type": "MissingWorkerResult",
            "error": f"worker exited with code {process.exitcode} without a result",
        }
    record.update(
        {
            "prime": prime,
            "branch": branch,
            "limit_seconds": GROEBNER_TIMEOUT_SECONDS,
        }
    )
    if record.get("status") == "completed":
        record["scope"] = (
            "A basis [1] decides only this affine branch after reduction modulo this "
            "prime and excludes p-integral lifts on this chart. It does not imply "
            "emptiness over Q; rational points may have denominators divisible by p."
            if record.get("contains_one")
            else "A completed proper basis is structural finite-field information only."
        )
    return record


def _serializable_lift(record: dict) -> dict:
    converted = dict(record)
    converted["modulus"] = str(record["modulus"])
    converted["residues"] = [str(value) for value in record["residues"]]
    return converted


def main() -> None:
    source, variables, expressions, polynomials = load_stored_system()
    variable_names = [str(variable) for variable in variables]
    derivatives = derivative_table(polynomials, len(variables))
    stored_witnesses = {
        int(row["prime"]): [int(value) for value in row["values_in_variable_order"]]
        for row in source["data"]["exact_modular_solution_witnesses"]
    }

    lift_records = []
    reconstruction_internal: dict[int, list[Fraction | None]] = {}
    exact_candidates = []
    for prime in TARGET_PRIMES:
        initial = stored_witnesses[prime]
        construction_residues = eval_system_mod(polynomials, initial, prime)
        if any(construction_residues):
            raise AssertionError(f"stored F_{prime} witness no longer solves the system")
        rank, pivots, _ = rref_mod(jacobian_mod(derivatives, initial, prime), prime)
        base = {
            "prime": prime,
            "stored_solution_verified": True,
            "stored_construction_residues": construction_residues,
            "jacobian_rank_mod_p": rank,
            "jacobian_row_rank_is_full": rank == len(polynomials),
            "initial_pivot_columns": pivots,
        }
        if rank == len(polynomials):
            exponent = target_exponent(prime)
            lifted = hensel_lift_full_row_rank(
                prime, initial, polynomials, derivatives, exponent
            )
            reconstruction = reconstruct_lift(lifted, variable_names)
            reconstruction_internal[prime] = reconstruction.pop("fractions")
            base.update(_serializable_lift(lifted))
            base["decimal_digits_floor"] = len(str(prime**exponent)) - 1
            base["local_section"] = {
                "newton_variables": [variable_names[index] for index in lifted["pivot_columns"]],
                "held_coordinates": [
                    variable_names[index]
                    for index in range(len(variables))
                    if index not in lifted["pivot_columns"]
                ],
                "held_coordinate_rule": "least nonnegative F_p representative, held exactly",
            }
            base["rational_reconstruction"] = reconstruction
            if reconstruction["all_coordinates_reconstructed"]:
                fractions = reconstruction_internal[prime]
                if any(value is None for value in fractions):
                    raise AssertionError("reconstruction flag is inconsistent")
                exact = _exact_candidate_checks(
                    [value for value in fractions if value is not None],
                    expressions,
                    variables,
                )
                exact["prime"] = prime
                exact_candidates.append(exact)
        else:
            base.update(
                rank_deficient_first_lift_attempt(
                    prime, initial, polynomials, derivatives
                )
            )
        lift_records.append(base)

    common_primes = sorted(reconstruction_internal)
    cross_prime = {
        "primes": common_primes,
        "consistently_reconstructed_coordinates": [],
        "inconsistently_reconstructed_coordinates": [],
        "unavailable_coordinates": [],
        "scope": (
            "Equality is compared only where both independent selected sections reconstruct. "
            "Agreement of a coordinate is partial arithmetic structure, not a common rational point."
        ),
    }
    if len(common_primes) >= 2:
        left_prime, right_prime = common_primes[:2]
        left = reconstruction_internal[left_prime]
        right = reconstruction_internal[right_prime]
        for name, left_value, right_value in zip(variable_names, left, right, strict=True):
            if left_value is None or right_value is None:
                cross_prime["unavailable_coordinates"].append(name)
            elif left_value == right_value:
                cross_prime["consistently_reconstructed_coordinates"].append(
                    {
                        "variable": name,
                        "value": str(left_value),
                        "note": "may be forced by the coordinate-holding section",
                    }
                )
            else:
                cross_prime["inconsistently_reconstructed_coordinates"].append(name)

    smallest_branch, smallest_ties = _smallest_branch_data()
    branch_certificate = exact_branch_constant_certificate(smallest_branch)
    groebner_attempts = [
        bounded_groebner_attempt(prime, smallest_branch)
        for prime in GROEBNER_PRIMES
    ]

    breakthrough = any(row.get("breakthrough_candidate") for row in exact_candidates)
    pattern_evidence = {
        "stored_solution_counts_by_pattern_available": False,
        "reason": (
            "The predecessor artifact stores one generic-chart witness per prime, not "
            "solution counts for all 32 patterns. Random-point Jacobian ranks are not solutions."
        ),
        "proven_nonempty_patterns_from_stored_witnesses": {
            str(prime): ["11111"] for prime in TARGET_PRIMES
        },
        "groebner_tested_pattern": smallest_branch,
        "exact_characteristic_zero_certificate": branch_certificate,
        "smallest_active_variable_count_ties": smallest_ties,
        "fresh_prime_search": [
            {
                "prime": prime,
                "status": "not_run_not_cheap",
                "reason": (
                    f"Uniform search for six equations has heuristic hit cost p^6={prime**6}; "
                    "no triangular finite-field solver is stored, so a tiny random sample would "
                    "supply no meaningful negative evidence."
                ),
            }
            for prime in (13, 17)
        ],
    }

    checks = []
    checks.append(
        {
            "name": "stored_system_hash_recomputed",
            "passed": _canonical_system_hash(
                source["data"]["minimal_unresolved_system"]["labels"], expressions
            )
            == source["data"]["minimal_unresolved_system"]["sha256_sympy_srepr"],
            "detail": source["data"]["minimal_unresolved_system"]["sha256_sympy_srepr"],
        }
    )
    checks.append(
        {
            "name": "stored_target_witnesses_reverified",
            "passed": all(row["stored_solution_verified"] for row in lift_records),
            "detail": {str(row["prime"]): row["stored_construction_residues"] for row in lift_records},
        }
    )
    observed_ranks = {str(row["prime"]): row["jacobian_rank_mod_p"] for row in lift_records}
    checks.append(
        {
            "name": "exact_jacobian_ranks_computed",
            "passed": observed_ranks == {"5": 6, "7": 5, "11": 6},
            "detail": observed_ranks,
        }
    )
    lifted_records = [row for row in lift_records if row["status"] == "lifted"]
    checks.append(
        {
            "name": "full_rank_points_lifted_beyond_60_decimal_digits",
            "passed": (
                {row["prime"] for row in lifted_records} == {5, 11}
                and all(row["decimal_digits_floor"] >= MIN_DECIMAL_DIGITS for row in lifted_records)
                and all(
                    not any(
                        eval_system_mod(
                            polynomials,
                            [int(value) for value in row["residues"]],
                            int(row["modulus"]),
                        )
                    )
                    for row in lifted_records
                )
            ),
            "detail": {
                str(row["prime"]): {
                    "exponent": row["precision_exponent"],
                    "decimal_digits_floor": row["decimal_digits_floor"],
                }
                for row in lifted_records
            },
        }
    )
    seven = next(row for row in lift_records if row["prime"] == 7)
    certificate = seven.get("certificate", {})
    checks.append(
        {
            "name": "F7_stored_point_nonlift_mod_49_certificate",
            "passed": (
                seven["status"] == "no_lift_mod_p_squared"
                and not any(certificate.get("left_kernel_times_jacobian_mod_p", [1]))
                and certificate.get("pairing_with_rhs_mod_p", 0) != 0
            ),
            "detail": certificate,
        }
    )
    reconstruction_bookkeeping_ok = True
    for row in lifted_records:
        modulus = int(row["modulus"])
        bound = int(row["rational_reconstruction"]["height_bound"])
        reconstruction_bookkeeping_ok &= 2 * bound * bound < modulus
        for coordinate, residue in zip(
            row["rational_reconstruction"]["coordinates"], row["residues"], strict=True
        ):
            if coordinate["status"] != "success":
                continue
            numerator = int(coordinate["numerator"])
            denominator = int(coordinate["denominator"])
            reconstruction_bookkeeping_ok &= (
                abs(numerator) <= bound
                and 1 <= denominator <= bound
                and (numerator - int(residue) * denominator) % modulus == 0
            )
    checks.append(
        {
            "name": "rational_reconstruction_bounds_and_congruences",
            "passed": bool(reconstruction_bookkeeping_ok),
            "detail": {
                str(row["prime"]): {
                    "height_bound": row["rational_reconstruction"]["height_bound"],
                    "successes": row["rational_reconstruction"]["successful_coordinate_count"],
                    "all_coordinates": row["rational_reconstruction"]["all_coordinates_reconstructed"],
                }
                for row in lifted_records
            },
        }
    )
    checks.append(
        {
            "name": "smallest_branch_exact_constant_certificate",
            "passed": (
                branch_certificate["identity_verified"]
                and branch_certificate["constant"] != 0
            ),
            "detail": {
                "branch": branch_certificate["branch"],
                "coefficients": branch_certificate["combination_coefficients"],
                "constant": branch_certificate["constant"],
            },
        }
    )
    checks.append(
        {
            "name": "bounded_groebner_retries_accounted",
            "passed": all(row["status"] in {"completed", "timeout"} for row in groebner_attempts),
            "detail": [
                {
                    "prime": row["prime"],
                    "status": row["status"],
                    "contains_one": row.get("contains_one"),
                }
                for row in groebner_attempts
            ],
        }
    )
    checks.append(
        {
            "name": "breakthrough_requires_three_exact_system_checks",
            "passed": (
                not breakthrough
                or any(
                    row["construction_passed"]
                    and len(row["holdouts"]) == 2
                    and all(holdout["passed"] for holdout in row["holdouts"])
                    for row in exact_candidates
                )
            ),
            "detail": "construction system plus exact 3x3x3 and 4x2x3 holdouts",
        }
    )

    output = {
        "meta": {
            "provenance": "experiments/e42_kw_hensel.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "source_artifact": "results/kac_ward/full_family.json",
            "source_system_sha256": source["data"]["minimal_unresolved_system"][
                "sha256_sympy_srepr"
            ],
            "arithmetic": "exact Python integers, Fraction, SymPy ZZ/GF(p); no floating point",
            "minimum_decimal_digits": MIN_DECIMAL_DIGITS,
            "groebner_timeout_seconds_per_prime": GROEBNER_TIMEOUT_SECONDS,
        },
        "data": {
            "headline": (
                "BREAKTHROUGH_CANDIDATE pending parent audit"
                if breakthrough
                else "UNRESOLVED: selected p-adic sections do not decide the characteristic-zero family"
            ),
            "method_scope": (
                "Six equations in 25 variables define a positive-dimensional construction variety. "
                "At full-row-rank points, six pivot coordinates are lifted while 19 coordinates "
                "are held. Reconstruction bounds apply only to these selected residue classes, "
                "not to every rational point of the variety."
            ),
            "lifts": lift_records,
            "cross_prime_reconstruction": cross_prime,
            "exact_reconstruction_candidates": exact_candidates,
            "pattern_evidence": pattern_evidence,
            "exact_smallest_branch_certificate": branch_certificate,
            "groebner_retries": groebner_attempts,
            "epistemic_status": (
                "A failed bounded reconstruction is exact evidence about one selected p-adic "
                "residue class, not a full-family no-go theorem. The exact constant certificate "
                "does prove branch 00000 empty over Q, while the other branches remain open. "
                "Finite-field data alone do not decide Q. The full scalar Kac--Ward family "
                "remains UNRESOLVED."
                if not breakthrough
                else "An exact candidate passed construction and two independent holdouts; label "
                "THEOREM-CANDIDATE only until an independent parent audit reproduces it."
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
