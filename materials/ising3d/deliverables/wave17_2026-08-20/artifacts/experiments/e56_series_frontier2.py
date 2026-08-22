"""Exact second-order finite-prefix structure frontiers for 3D Ising series.

The searches are deliberately finite.  They certify ranks of coefficient
matrices over Q; they do not assert non-D-finiteness or differential
transcendence of either infinite free-energy series.

Run from the repository root with

    .venv/bin/python experiments/e56_series_frontier2.py
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

SCRIPT = "experiments/e56_series_frontier2.py"
INTERPRETER = ".venv/bin/python"
ROOT = Path(__file__).resolve().parents[1]
HT_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
LT_PATH = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
RESULT_PATH = ROOT / "results" / "series" / "frontier2.json"

RationalSeries = tuple[Fraction, ...]
MatrixRows = list[list[Fraction]]

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _fraction_strings(values: Iterable[Fraction | int]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _load_series(path: Path) -> tuple[RationalSeries, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    coefficients = tuple(Fraction(value) for value in payload["data"]["coefficients"])
    achieved_order = int(payload["data"]["achieved_order"])
    if len(coefficients) != achieved_order + 1:
        raise AssertionError(f"{path.name}: coefficient count does not match achieved order")
    if any(coefficients[order] for order in range(1, len(coefficients), 2)):
        raise AssertionError(f"{path.name}: expected the stored free-energy series to be even")
    return coefficients, {
        "artifact": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "variable": payload["data"]["variable"],
        "achieved_order": achieved_order,
        "coefficient_count": len(coefficients),
        "coefficients": _fraction_strings(coefficients),
    }


def _valuation(series: Sequence[Fraction]) -> int:
    for order, coefficient in enumerate(series):
        if coefficient:
            return order
    raise ValueError("finite series is identically zero")


def _differentiate(series: RationalSeries) -> RationalSeries:
    return tuple(Fraction(order) * series[order] for order in range(1, len(series)))


def _convolve(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> RationalSeries:
    result = [Fraction(0) for _ in range(order + 1)]
    for i, first in enumerate(left[: order + 1]):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                result[i + j] += first * second
    return tuple(result)


def _powers(series: RationalSeries, maximum: int, order: int) -> tuple[RationalSeries, ...]:
    values = [(Fraction(1),) + (Fraction(0),) * order]
    truncated = tuple(series[: order + 1])
    if len(truncated) < order + 1:
        truncated += (Fraction(0),) * (order + 1 - len(truncated))
    for _ in range(maximum):
        values.append(_convolve(values[-1], truncated, order))
    return tuple(values)


def _primitive_vector(vector: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = math.lcm(*(value.denominator for value in vector))
    integers = [value.numerator * (denominator // value.denominator) for value in vector]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor:
        integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def _dot(row: Sequence[Fraction], vector: Sequence[int]) -> Fraction:
    return sum(
        (value * coefficient for value, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _integer_rows(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int]
) -> list[list[int]]:
    integer_rows: list[list[int]] = []
    for row in rows:
        selected = [Fraction(row[column]) for column in columns]
        denominator = math.lcm(*(value.denominator for value in selected)) if selected else 1
        integer_rows.append(
            [value.numerator * (denominator // value.denominator) for value in selected]
        )
    return integer_rows


def _fraction_free_rank_and_rows(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int]
) -> tuple[int, tuple[int, ...]]:
    """Exact Q-rank by rowwise denominator clearing and Bareiss elimination."""

    matrix = _integer_rows(rows, columns)
    if not matrix or not columns:
        return 0, ()
    row_ids = list(range(len(matrix)))
    row_count = len(matrix)
    column_count = len(columns)
    pivot_row = 0
    previous_pivot = 1
    selected_rows: list[int] = []
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        row_ids[pivot_row], row_ids[pivot] = row_ids[pivot], row_ids[pivot_row]
        pivot_value = matrix[pivot_row][column]
        for row in range(pivot_row + 1, row_count):
            multiplier = matrix[row][column]
            for later_column in range(column + 1, column_count):
                numerator = (
                    pivot_value * matrix[row][later_column]
                    - multiplier * matrix[pivot_row][later_column]
                )
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise ArithmeticError("non-exact division in Bareiss rank elimination")
                matrix[row][later_column] = quotient
            matrix[row][column] = 0
        selected_rows.append(row_ids[pivot_row])
        previous_pivot = pivot_value
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row, tuple(selected_rows)


def _determinant(rows: Sequence[Sequence[Fraction]]) -> Fraction:
    """Exact determinant by row scaling followed by integer Bareiss elimination."""

    size = len(rows)
    if size == 0:
        return Fraction(1)
    if any(len(row) != size for row in rows):
        raise ValueError("determinant requires a square matrix")
    row_scales = []
    matrix = []
    for row in rows:
        values = [Fraction(value) for value in row]
        scale = math.lcm(*(value.denominator for value in values))
        row_scales.append(scale)
        matrix.append([value.numerator * (scale // value.denominator) for value in values])
    sign = 1
    previous_pivot = 1
    for diagonal in range(size - 1):
        pivot = next(
            (row for row in range(diagonal, size) if matrix[row][diagonal]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != diagonal:
            matrix[diagonal], matrix[pivot] = matrix[pivot], matrix[diagonal]
            sign = -sign
        pivot_value = matrix[diagonal][diagonal]
        for row in range(diagonal + 1, size):
            multiplier = matrix[row][diagonal]
            for column in range(diagonal + 1, size):
                numerator = (
                    pivot_value * matrix[row][column]
                    - multiplier * matrix[diagonal][column]
                )
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise ArithmeticError("non-exact division in Bareiss determinant")
                matrix[row][column] = quotient
            matrix[row][diagonal] = 0
        previous_pivot = pivot_value
    return Fraction(sign * matrix[-1][-1], math.prod(row_scales))


def _nullspace(rows: Sequence[Sequence[Fraction]], column_count: int) -> tuple[tuple[int, ...], ...]:
    """Primitive integer basis of an exact rational kernel."""

    matrix = [[Fraction(value) for value in row] for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                value - multiplier * pivot_value
                for value, pivot_value in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    free_columns = [column for column in range(column_count) if column not in pivot_columns]
    basis = []
    for free_column in free_columns:
        vector = [Fraction(0) for _ in range(column_count)]
        vector[free_column] = Fraction(1)
        for row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = -matrix[row][free_column]
        primitive = _primitive_vector(vector)
        if any(_dot(source_row, primitive) for source_row in rows):
            raise AssertionError("constructed nullspace vector has a nonzero residual")
        basis.append(primitive)
    return tuple(basis)


def _minor_certificate(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int], phase: str
) -> dict:
    rank, selected_rows = _fraction_free_rank_and_rows(rows, columns)
    if rank != len(columns):
        raise AssertionError("requested maximal-minor certificate for a rank-deficient matrix")
    minor = [[rows[row][column] for column in columns] for row in selected_rows]
    determinant = _determinant(minor)
    if not determinant:
        raise AssertionError("fraction-free pivots produced a zero certifying determinant")
    return {
        "type": "NONZERO_MAXIMAL_MINOR_OVER_Q",
        "phase": phase,
        "matrix_shape": [len(rows), len(columns)],
        "row_orders": list(selected_rows),
        "column_indices": list(columns),
        "determinant": str(determinant),
        "method": "rowwise denominator clearing and exact integer Bareiss elimination",
    }


def _store_certificate(certificates: dict[str, dict], certificate: dict, context: dict) -> str:
    certificate_id = f"C{len(certificates) + 1:04d}"
    certificates[certificate_id] = {**context, **certificate}
    return f"/data/certificates/{certificate_id}"


def _independent_row_indices(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int]
) -> tuple[int, ...]:
    rank, selected_rows = _fraction_free_rank_and_rows(rows, columns)
    if rank != len(columns):
        raise AssertionError("columns are not independent on the supplied rows")
    return selected_rows


def _support_progression(
    rows: Sequence[Sequence[Fraction]], formal_valuations: Sequence[int]
) -> int | None:
    """Return a common arithmetic support modulus, if one is exact.

    Only a progression actually respected by every nonzero matrix entry is
    returned.  Formal valuations fix the expected residue of each column.
    """

    equation_count = len(rows)
    for modulus in range(2, equation_count + 1):
        if all(
            all(
                (not row[column]) or order % modulus == formal_valuations[column] % modulus
                for column in range(len(formal_valuations))
            )
            for order, row in enumerate(rows)
        ):
            return modulus
    return None


def _unobservable_subspace_certificate(
    rows: MatrixRows,
    monomials: Sequence[tuple[int, ...]],
    formal_valuations: Sequence[int],
    full_basis: Sequence[tuple[int, ...]],
    certificates: dict[str, dict],
    context: dict,
) -> tuple[dict, str] | None:
    """Certify that every surviving direction starts after the known prefix.

    Besides literal zero columns, parity-supported linear systems can have
    observable combinations whose first *unfixed* coefficient is beyond the
    truncation: there are fewer populated coefficient orders than parameters.
    We prove this by selecting exact independent known rows, verifying the
    kernel, and recording the first unavailable order in each exact support
    progression.
    """

    equation_count = len(rows)
    unknown_count = len(monomials)
    modulus = _support_progression(rows, formal_valuations)
    if modulus is None:
        return None
    groups: dict[int, list[int]] = {}
    for column, valuation in enumerate(formal_valuations):
        groups.setdefault(valuation % modulus, []).append(column)
    components = []
    kernel_dimension = 0
    for residue, group_columns_list in sorted(groups.items()):
        group_columns = tuple(group_columns_list)
        populated_orders = tuple(
            order for order in range(equation_count) if order % modulus == residue
        )
        group_rows = [rows[order] for order in populated_orders]
        group_rank, _ = _fraction_free_rank_and_rows(group_rows, group_columns)
        group_nullity = len(group_columns) - group_rank
        kernel_dimension += group_nullity
        if group_nullity == 0:
            selected_local = _independent_row_indices(group_rows, group_columns)
            selected_orders = [populated_orders[index] for index in selected_local]
            certificate = _minor_certificate(
                [rows[order] for order in selected_orders],
                group_columns,
                f"support_residue_{residue}",
            )
            determinant = certificate["determinant"]
        else:
            selected_orders = []
            determinant = None
        next_unfixed_order = (
            residue
            if not populated_orders
            else populated_orders[-1] + modulus
        )
        components.append(
            {
                "support_modulus": modulus,
                "support_residue": residue,
                "column_indices": list(group_columns),
                "monomials": [list(monomials[column]) for column in group_columns],
                "known_populated_orders": list(populated_orders),
                "known_rank": group_rank,
                "known_nullity": group_nullity,
                "certifying_row_orders_when_full_rank": selected_orders,
                "certifying_minor_determinant_when_full_rank": determinant,
                "first_unfixed_supported_order": next_unfixed_order,
                "last_available_order": equation_count - 1,
            }
        )
    if kernel_dimension != len(full_basis):
        return None
    if any(
        component["known_nullity"] > 0
        and component["first_unfixed_supported_order"] < equation_count
        for component in components
    ):
        return None
    if any(any(_dot(row, vector) for row in rows) for vector in full_basis):
        raise AssertionError("purported unobservable kernel has a known nonzero residual")
    certificate = {
        "type": "TRUNCATION_UNOBSERVABLE_SUPPORT_KERNEL",
        "phase": "all_available_equations",
        "matrix_shape": [equation_count, unknown_count],
        "full_kernel_basis": [list(vector) for vector in full_basis],
        "support_components": components,
        "proof": (
            "every matrix column is supported in its displayed residue class; exact "
            "component ranks account for the full kernel, and each component's next "
            "unfixed supported order exceeds the last available coefficient order"
        ),
    }
    pointer = _store_certificate(certificates, certificate, context)
    return certificate, pointer


def _analyze_budget(
    rows: MatrixRows,
    monomials: Sequence[tuple[int, ...]],
    formal_valuations: Sequence[int],
    certificates: dict[str, dict],
    context: dict,
) -> dict:
    """Fit U+2 prefix equations, then apply every remaining exact holdout."""

    equation_count = len(rows)
    unknown_count = len(monomials)
    training_count = unknown_count + 2
    if training_count > equation_count:
        raise AssertionError("frontier admitted a dimensionally unsupported budget")
    columns = tuple(range(unknown_count))
    training_rows = rows[:training_count]
    training_rank, _ = _fraction_free_rank_and_rows(training_rows, columns)
    base = {
        "unknown_count": unknown_count,
        "available_equation_count": equation_count,
        "training_equation_count": training_count,
        "holdout_equation_count": equation_count - training_count,
        "training_orders": [0, training_count - 1],
        "holdout_orders": (
            [training_count, equation_count - 1] if training_count < equation_count else []
        ),
        "monomials": [list(monomial) for monomial in monomials],
        "formal_monomial_valuations": list(formal_valuations),
        "training_rank": training_rank,
        "training_nullity": unknown_count - training_rank,
    }
    if training_rank == unknown_count:
        certificate = _minor_certificate(training_rows, columns, "training_prefix")
        return {
            **base,
            "verdict": "NO_RELATION_AT_BUDGET",
            "certificate_pointer": _store_certificate(certificates, certificate, context),
            "detail": (
                "the U+2 prefix coefficient matrix has exact full column rank over Q"
            ),
        }

    training_basis = _nullspace(training_rows, unknown_count)
    if len(training_basis) != unknown_count - training_rank:
        raise AssertionError("fraction-free rank and Fraction nullspace disagree")
    residual_constraints: MatrixRows = []
    holdout_checks = []
    first_refuting_order = None
    for order in range(training_count, equation_count):
        residuals = [_dot(rows[order], vector) for vector in training_basis]
        residual_constraints.append(residuals)
        residual_rank, _ = _fraction_free_rank_and_rows(
            residual_constraints, tuple(range(len(training_basis)))
        )
        surviving_dimension = len(training_basis) - residual_rank
        holdout_checks.append(
            {
                "order": order,
                "residual_on_training_kernel_basis": _fraction_strings(residuals),
                "surviving_kernel_dimension": surviving_dimension,
            }
        )
        if surviving_dimension == 0 and first_refuting_order is None:
            first_refuting_order = order

    full_rank, _ = _fraction_free_rank_and_rows(rows, columns)
    full_nullity = unknown_count - full_rank
    candidate_fields = {
        "training_kernel_basis": [list(vector) for vector in training_basis],
        "holdout_checks": holdout_checks,
        "first_refuting_holdout_order": first_refuting_order,
        "full_rank": full_rank,
        "full_nullity": full_nullity,
    }
    if full_rank == unknown_count:
        if first_refuting_order is None:
            raise AssertionError("full rank after a training kernel lacks a refuting holdout")
        certificate = _minor_certificate(rows, columns, "all_available_equations")
        return {
            **base,
            **candidate_fields,
            "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
            "certificate_pointer": _store_certificate(certificates, certificate, context),
            "detail": (
                "the training kernel is eliminated by withheld exact coefficient equations"
            ),
        }

    zero_columns = [
        column for column in columns if all(row[column] == 0 for row in rows)
    ]
    if any(formal_valuations[column] < equation_count for column in zero_columns):
        raise AssertionError("a zero monomial column has an observable formal valuation")

    full_basis = _nullspace(rows, unknown_count)
    if len(full_basis) != full_nullity:
        raise AssertionError("candidate kernel dimension disagrees with fraction-free rank")
    visible_columns = tuple(column for column in columns if column not in zero_columns)
    visible_rank, _ = _fraction_free_rank_and_rows(rows, visible_columns)
    if full_nullity == len(zero_columns) and visible_rank == len(visible_columns):
        certificate = _minor_certificate(rows, visible_columns, "observable_columns")
        unobservable = [
            {
                "column_index": column,
                "monomial": list(monomials[column]),
                "lowest_possible_relation_order": formal_valuations[column],
                "last_available_relation_order": equation_count - 1,
            }
            for column in zero_columns
        ]
        return {
            **base,
            **candidate_fields,
            "verdict": "TRUNCATION_UNOBSERVABLE",
            "zero_column_indices_at_known_order": zero_columns,
            "unobservable_monomials": unobservable,
            "all_survivors_are_exactly_unobservable_zero_columns": True,
            "certificate_pointer": _store_certificate(certificates, certificate, context),
            "detail": (
                "the full kernel is exactly the coordinate span of monomials whose "
                "first possible coefficient lies after every derivative-safe equation"
            ),
        }

    support_certificate = _unobservable_subspace_certificate(
        rows,
        monomials,
        formal_valuations,
        full_basis,
        certificates,
        context,
    )
    if support_certificate is not None:
        certificate, pointer = support_certificate
        return {
            **base,
            **candidate_fields,
            "verdict": "TRUNCATION_UNOBSERVABLE",
            "zero_column_indices_at_known_order": zero_columns,
            "unobservable_monomials": [],
            "support_components": certificate["support_components"],
            "all_survivors_are_exactly_unobservable_zero_columns": False,
            "all_survivors_are_truncation_unobservable": True,
            "full_kernel_basis": [list(vector) for vector in full_basis],
            "certificate_pointer": pointer,
            "detail": (
                "the exact surviving kernel is forced by sparse support: its first "
                "unfixed supported relation order lies beyond the available prefix"
            ),
        }

    certificate = {
        "type": "EXACT_FULL_KERNEL",
        "phase": "all_available_equations",
        "matrix_shape": [equation_count, unknown_count],
        "kernel_basis": [list(vector) for vector in full_basis],
        "residuals": [
            _fraction_strings(_dot(row, vector) for row in rows) for vector in full_basis
        ],
    }
    return {
        **base,
        **candidate_fields,
        "verdict": "CANDIDATE_SURVIVES_PARENT_REVIEW",
        "zero_column_indices_at_known_order": zero_columns,
        "full_kernel_basis": [list(vector) for vector in full_basis],
        "certificate_pointer": _store_certificate(certificates, certificate, context),
        "audit_warning": (
            "POTENTIAL FINITE-PREFIX DISCOVERY: exact candidate passed every known "
            "coefficient; do not promote without parent review and new coefficients"
        ),
        "detail": "an observable exact kernel remains after all available coefficients",
    }


def _jet_matrix(
    series: RationalSeries, monomials: Sequence[tuple[int, int, int, int]]
) -> tuple[MatrixRows, list[int]]:
    equation_count = len(series) - 2
    derivatives = (
        tuple(series[:equation_count]),
        tuple(_differentiate(series)[:equation_count]),
        tuple(_differentiate(_differentiate(series))[:equation_count]),
    )
    maxima = [max(monomial[index] for monomial in monomials) for index in (1, 2, 3)]
    powers = tuple(
        _powers(derivatives[index], maxima[index], equation_count - 1)
        for index in range(3)
    )
    columns = []
    for degree_t, degree_f, degree_f_prime, degree_f_second in monomials:
        product = _convolve(
            _convolve(
                powers[0][degree_f], powers[1][degree_f_prime], equation_count - 1
            ),
            powers[2][degree_f_second],
            equation_count - 1,
        )
        columns.append(
            [
                product[order - degree_t] if order >= degree_t else Fraction(0)
                for order in range(equation_count)
            ]
        )
    rows = [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]
    valuations = tuple(_valuation(derivative) for derivative in derivatives)
    formal_valuations = [
        degree_t
        + degree_f * valuations[0]
        + degree_f_prime * valuations[1]
        + degree_f_second * valuations[2]
        for degree_t, degree_f, degree_f_prime, degree_f_second in monomials
    ]
    return rows, formal_valuations


def _total_degree_monomials(degree: int) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (degree_t, degree_f, degree_f_prime, degree_f_second)
        for degree_f_second in range(degree + 1)
        for degree_f_prime in range(degree + 1 - degree_f_second)
        for degree_f in range(degree + 1 - degree_f_second - degree_f_prime)
        for degree_t in range(
            degree + 1 - degree_f_second - degree_f_prime - degree_f
        )
    )


def _jet_bidegree_monomials(
    degree_t: int, jet_total_degree: int
) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (power_t, degree_f, degree_f_prime, degree_f_second)
        for degree_f_second in range(jet_total_degree + 1)
        for degree_f_prime in range(jet_total_degree + 1 - degree_f_second)
        for degree_f in range(
            jet_total_degree + 1 - degree_f_second - degree_f_prime
        )
        for power_t in range(degree_t + 1)
    )


def _second_order_frontier(
    series_name: str,
    series: RationalSeries,
    certificates: dict[str, dict],
) -> list[dict]:
    equation_count = len(series) - 2
    maximum_unknowns = equation_count - 2
    specifications = []
    total_degree = 1
    while True:
        monomials = _total_degree_monomials(total_degree)
        if len(monomials) > maximum_unknowns:
            break
        specifications.append(
            ("second_order_total_degree", {"total_degree": total_degree}, monomials)
        )
        total_degree += 1
    jet_degree = 1
    while True:
        base_monomials = _jet_bidegree_monomials(0, jet_degree)
        if len(base_monomials) > maximum_unknowns:
            break
        degree_t = 0
        while True:
            monomials = _jet_bidegree_monomials(degree_t, jet_degree)
            if len(monomials) > maximum_unknowns:
                break
            specifications.append(
                (
                    "second_order_jet_bidegree",
                    {"degree_t": degree_t, "jet_total_degree": jet_degree},
                    monomials,
                )
            )
            degree_t += 1
        jet_degree += 1
    specifications.sort(
        key=lambda item: (len(item[2]), item[0], tuple(item[1].values()))
    )
    table = []
    for family, parameters, monomials in specifications:
        if not any(monomial[3] for monomial in monomials):
            raise AssertionError("second-order search space omitted F''")
        rows, valuations = _jet_matrix(series, monomials)
        context = {
            "series": series_name,
            "ansatz_family": family,
            "parameters": parameters,
        }
        table.append(
            {
                "ansatz_family": family,
                "parameters": parameters,
                "highest_derivative_order": 2,
                **_analyze_budget(rows, monomials, valuations, certificates, context),
            }
        )
    return table


def _theta_columns(series: RationalSeries, order: int) -> tuple[RationalSeries, ...]:
    return tuple(
        tuple(Fraction(index**power) * coefficient for index, coefficient in enumerate(series))
        for power in range(order + 1)
    )


def _linear_ode_matrix(
    series: RationalSeries, order: int, degree: int
) -> tuple[MatrixRows, tuple[tuple[int, int], ...], list[int]]:
    equation_count = len(series)
    theta = _theta_columns(series, order)
    monomials = tuple(
        (degree_t, theta_power)
        for theta_power in range(order + 1)
        for degree_t in range(degree + 1)
    )
    columns = [
        [
            theta[theta_power][coefficient_order - degree_t]
            if coefficient_order >= degree_t
            else Fraction(0)
            for coefficient_order in range(equation_count)
        ]
        for degree_t, theta_power in monomials
    ]
    rows = [
        [columns[column][coefficient_order] for column in range(len(columns))]
        for coefficient_order in range(equation_count)
    ]
    valuation = _valuation(series)
    valuations = [degree_t + valuation for degree_t, _ in monomials]
    return rows, monomials, valuations


def _linear_ode_frontier(
    series_name: str,
    series: RationalSeries,
    certificates: dict[str, dict],
) -> list[dict]:
    equation_count = len(series)
    table = []
    for order in range(3):
        degree = 0
        while (order + 1) * (degree + 1) <= equation_count - 2:
            rows, monomials, valuations = _linear_ode_matrix(series, order, degree)
            parameters = {"order": order, "polynomial_degree": degree}
            context = {
                "series": series_name,
                "ansatz_family": "linear_euler_ode",
                "parameters": parameters,
            }
            table.append(
                {
                    "ansatz_family": "linear_euler_ode",
                    "parameters": parameters,
                    **_analyze_budget(
                        rows, monomials, valuations, certificates, context
                    ),
                }
            )
            degree += 1
    return table


def _mahler_matrix(
    series: RationalSeries,
    degree_t: int,
    dependent_total_degree: int,
) -> tuple[MatrixRows, tuple[tuple[int, int, int], ...], list[int]]:
    equation_count = len(series)
    composed = tuple(
        series[order // 2] if order % 2 == 0 else Fraction(0)
        for order in range(equation_count)
    )
    monomials = tuple(
        (power_t, degree_f, degree_composed)
        for degree_composed in range(dependent_total_degree + 1)
        for degree_f in range(dependent_total_degree + 1 - degree_composed)
        for power_t in range(degree_t + 1)
    )
    f_powers = _powers(series, dependent_total_degree, equation_count - 1)
    composed_powers = _powers(composed, dependent_total_degree, equation_count - 1)
    columns = []
    for power_t, degree_f, degree_composed in monomials:
        product = _convolve(
            f_powers[degree_f], composed_powers[degree_composed], equation_count - 1
        )
        columns.append(
            [
                product[order - power_t] if order >= power_t else Fraction(0)
                for order in range(equation_count)
            ]
        )
    rows = [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]
    valuation = _valuation(series)
    valuations = [
        power_t + degree_f * valuation + degree_composed * (2 * valuation)
        for power_t, degree_f, degree_composed in monomials
    ]
    return rows, monomials, valuations


def _mahler_frontier(
    series_name: str,
    series: RationalSeries,
    certificates: dict[str, dict],
) -> list[dict]:
    equation_count = len(series)
    maximum_unknowns = equation_count - 2
    table = []
    for dependent_degree in (1, 2):
        degree_t = 0
        while True:
            rows, monomials, valuations = _mahler_matrix(
                series, degree_t, dependent_degree
            )
            if len(monomials) > maximum_unknowns:
                break
            parameters = {
                "degree_t": degree_t,
                "dependent_total_degree": dependent_degree,
            }
            context = {
                "series": series_name,
                "ansatz_family": "mahler_bidegree",
                "parameters": parameters,
            }
            table.append(
                {
                    "ansatz_family": "mahler_bidegree",
                    "parameters": parameters,
                    **_analyze_budget(
                        rows, monomials, valuations, certificates, context
                    ),
                }
            )
            degree_t += 1
    table.sort(key=lambda row: (row["unknown_count"], tuple(row["parameters"].values())))
    return table


def _verdict_counts(table: Sequence[dict]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in table:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    return counts


def _series_table(table: list[dict]) -> dict:
    return {
        "equation_count": table[0]["available_equation_count"] if table else 0,
        "frontier_row_count": len(table),
        "largest_unknown_count": max(row["unknown_count"] for row in table),
        "verdict_counts": _verdict_counts(table),
        "frontier": table,
    }


def _summary_rows(frontiers: dict) -> list[dict]:
    rows = []
    for frontier_name, frontier in frontiers.items():
        for series_name in ("HT", "LT"):
            table = frontier[series_name]
            rows.append(
                {
                    "ansatz_family": frontier_name,
                    "series": series_name,
                    "available_equation_count": table["equation_count"],
                    "frontier_row_count": table["frontier_row_count"],
                    "largest_unknown_count": table["largest_unknown_count"],
                    "verdict_counts": table["verdict_counts"],
                }
            )
    return rows


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    ht, ht_metadata = _load_series(HT_PATH)
    lt, lt_metadata = _load_series(LT_PATH)
    certificates: dict[str, dict] = {}

    second_ht = _second_order_frontier("HT", ht, certificates)
    second_lt = _second_order_frontier("LT", lt, certificates)
    linear_ht = _linear_ode_frontier("HT", ht, certificates)
    linear_lt = _linear_ode_frontier("LT", lt, certificates)
    mahler_ht = _mahler_frontier("HT", ht, certificates)
    mahler_lt = _mahler_frontier("LT", lt, certificates)

    frontiers = {
        "second_order_differential_algebraic": {
            "ansatz": "P(t,F,F',F'')=0",
            "families": {
                "second_order_total_degree": "a+b+c+d <= D",
                "second_order_jet_bidegree": "0<=a<=d_t and b+c+d<=D_jet",
            },
            "derivative_safe_rule": (
                "N coefficients through t^(N-1) give N-2 safe equations through "
                "t^(N-3) when F'' occurs"
            ),
            "HT": _series_table(second_ht),
            "LT": _series_table(second_lt),
            "scope": (
                "finite total-degree and jet-bidegree polynomial spaces only; this is "
                "not a proof of differential-algebraic independence"
            ),
        },
        "linear_euler_ode_order_at_most_2": {
            "ansatz": "sum_{j=0}^r Q_j(t) theta^j F(t)=0, theta=t*d/dt, r<=2",
            "degree_rule": "all common polynomial degrees d with (r+1)(d+1)<=N-2",
            "HT": _series_table(linear_ht),
            "LT": _series_table(linear_lt),
            "scope": (
                "finite polynomial-coefficient Euler-ODE spaces only; this is not a "
                "proof of non-D-finiteness"
            ),
        },
        "mahler_F_t2_vs_F_t": {
            "ansatz": (
                "M(t,F(t),F(t^2))=sum m_abc t^a F(t)^b F(t^2)^c with "
                "b+c<=2"
            ),
            "degree_rule": (
                "dependent total degree 1 or 2 and every t-degree allowed by U+2<=N"
            ),
            "HT": _series_table(mahler_ht),
            "LT": _series_table(mahler_lt),
            "scope": (
                "low dependent-degree finite Mahler-type polynomial spaces only; no "
                "general functional-transcendence conclusion"
            ),
        },
    }
    all_tables = (second_ht, second_lt, linear_ht, linear_lt, mahler_ht, mahler_lt)
    all_rows = [row for table in all_tables for row in table]
    unexplained = [
        row for row in all_rows if row["verdict"] == "CANDIDATE_SURVIVES_PARENT_REVIEW"
    ]
    unobservable = [
        row for row in all_rows if row["verdict"] == "TRUNCATION_UNOBSERVABLE"
    ]
    refuted = [
        row for row in all_rows if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ]

    checks: list[dict] = []
    _record(
        checks,
        "input orders",
        len(ht) == 23 and len(lt) == 33,
        f"loaded exact HT v^22 ({len(ht)} coefficients) and LT x^32 ({len(lt)} coefficients)",
    )
    _record(
        checks,
        "full budget enumeration",
        all(
            row["training_equation_count"] == row["unknown_count"] + 2
            and row["training_equation_count"] <= row["available_equation_count"]
            for row in all_rows
        ),
        f"all {len(all_rows)} frontier rows use the predetermined U+2 training rule",
    )
    _record(
        checks,
        "second derivative included",
        all(
            row["highest_derivative_order"] == 2
            and any(monomial[3] for monomial in row["monomials"])
            for row in second_ht + second_lt
        ),
        f"all {len(second_ht) + len(second_lt)} second-order spaces contain F''",
    )
    max_linear_degrees = {
        "HT": {
            order: max(
                row["parameters"]["polynomial_degree"]
                for row in linear_ht
                if row["parameters"]["order"] == order
            )
            for order in range(3)
        },
        "LT": {
            order: max(
                row["parameters"]["polynomial_degree"]
                for row in linear_lt
                if row["parameters"]["order"] == order
            )
            for order in range(3)
        },
    }
    _record(
        checks,
        "linear ODE true-budget endpoints",
        max_linear_degrees == {"HT": {0: 20, 1: 9, 2: 6}, "LT": {0: 30, 1: 14, 2: 9}},
        f"maximal common polynomial degrees are {max_linear_degrees}",
    )
    _record(
        checks,
        "holdout refutations present",
        bool(refuted)
        and all(
            row["holdout_equation_count"] > 0
            and row["first_refuting_holdout_order"] is not None
            for row in refuted
        ),
        f"{len(refuted)} prefix candidates are eliminated by exact withheld coefficients",
    )
    _record(
        checks,
        "truncation-unobservable classification",
        all(
            (
                row["all_survivors_are_exactly_unobservable_zero_columns"]
                and all(
                    item["lowest_possible_relation_order"]
                    > item["last_available_relation_order"]
                    for item in row["unobservable_monomials"]
                )
            )
            or (
                row.get("all_survivors_are_truncation_unobservable", False)
                and all(
                    component["first_unfixed_supported_order"]
                    > component["last_available_order"]
                    for component in row["support_components"]
                    if component["known_nullity"] > 0
                )
            )
            for row in unobservable
        ),
        f"{len(unobservable)} full-system kernels have exact beyond-prefix observability certificates",
    )
    _record(
        checks,
        "zero unexplained survivors",
        not unexplained,
        (
            "no observable candidate passes every known coefficient"
            if not unexplained
            else f"{len(unexplained)} observable candidates require parent review"
        ),
    )
    _record(
        checks,
        "certificate coverage",
        len(certificates) == len(all_rows)
        and all(row["certificate_pointer"].startswith("/data/certificates/C") for row in all_rows)
        and all(Fraction(certificate.get("determinant", "1")) != 0 for certificate in certificates.values()),
        f"all {len(all_rows)} rows point to one exact certificate; {len(certificates)} certificates stored",
    )

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": INTERPRETER,
            "method": (
                "exact Fraction coefficient matrices; rowwise denominator clearing; "
                "fraction-free integer Bareiss ranks and maximal minors; exact Fraction "
                "kernel and holdout residual checks"
            ),
            "python_version": platform.python_version(),
            "benchmark_used_for_fit_selection_or_validation": False,
        },
        "data": {
            "claim_tag": "[COMPUTATION]",
            "input_series": {"HT": ht_metadata, "LT": lt_metadata},
            "search_protocol": {
                "budget_rule": (
                    "for U unknown homogeneous coefficients use equations 0 through U+1 "
                    "as training and every later safe equation as strict holdout"
                ),
                "selection_independence": (
                    "degree grids and the U+2 split depend only on series length and ansatz "
                    "dimension, never on coefficient values, residuals, or K_c"
                ),
                "verdict_meanings": {
                    "NO_RELATION_AT_BUDGET": (
                        "a stored nonzero exact maximal minor proves full column rank already "
                        "on the training prefix"
                    ),
                    "CANDIDATE_REFUTED_BY_HOLDOUT": (
                        "the training kernel is killed by later exact coefficients; a full "
                        "system nonzero maximal minor is stored"
                    ),
                    "TRUNCATION_UNOBSERVABLE": (
                        "the full kernel is exactly accounted for by zero monomial "
                        "columns beyond the safe order or by exact sparse-support "
                        "components whose first unfixed supported order is beyond the prefix"
                    ),
                    "CANDIDATE_SURVIVES_PARENT_REVIEW": (
                        "an observable exact kernel passes all known data and must be audited"
                    ),
                },
            },
            "frontier_summary": _summary_rows(frontiers),
            "frontiers": frontiers,
            "survivor_audit": {
                "unexplained_survivor_count": len(unexplained),
                "truncation_unobservable_row_count": len(unobservable),
                "holdout_refuted_row_count": len(refuted),
                "parent_review_rows": unexplained,
            },
            "certificates": certificates,
            "scope_warning": (
                "These exact negative frontiers bound finite structure searches.  They do "
                "not prove non-D-finiteness, non-differential-algebraicity, differential "
                "transcendence, or the absence of relations outside the displayed budgets."
            ),
        },
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more exact frontier checks failed")
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e56_series_frontier2: {error}")
        raise
