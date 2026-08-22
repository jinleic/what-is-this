"""Independent exhaustive exact replay of the expanded series-structure grid.

Run from the repository root:

    PYTHONPATH=src .venv/bin/python tests/test_series_grid.py

This test does not import the producer experiment.  It rebuilds every matrix
from raw coefficient artifacts, uses separate exact Fraction/Bareiss routines,
and verifies every stored certificate, holdout ladder, survivor support proof,
and legacy LT replay.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "series" / "grid.json"
HT_INPUT = ROOT / "results" / "series" / "ht_v28.json"
LT_INPUT = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"

ALLOWED_VERDICTS = {
    "NO_RELATION_AT_BUDGET",
    "CANDIDATE_REFUTED_BY_HOLDOUT",
    "CANDIDATE_SURVIVES",
}
FAMILIES = (
    "algebraic",
    "first_order_differential_algebraic",
    "linear_euler_ode",
    "mahler",
)
LEGACY_LT_ROWS = (
    ("algebraic", {"degree_x": 0, "degree_f": 6}),
    ("algebraic", {"degree_x": 0, "degree_f": 7}),
    ("algebraic", {"degree_x": 0, "degree_f": 8}),
    ("algebraic", {"degree_x": 0, "degree_f": 9}),
    ("algebraic", {"degree_x": 0, "degree_f": 10}),
    ("algebraic", {"degree_x": 0, "degree_f": 11}),
    ("algebraic", {"degree_x": 0, "degree_f": 12}),
    ("algebraic", {"degree_x": 0, "degree_f": 13}),
    ("algebraic", {"degree_x": 1, "degree_f": 6}),
    ("algebraic", {"degree_x": 0, "degree_f": 14}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 8}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 9}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 10}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 11}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 5, "degree_f_prime": 1}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 12}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 0, "degree_f_prime": 13}),
    ("differential_algebraic_order_1", {"degree_x": 0, "degree_f": 6, "degree_f_prime": 1}),
)

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


# ---------------------------------------------------------------------------
# Independent exact linear algebra over Q.  This is intentionally separate
# from e56: it clears each row's denominators and performs integer Bareiss.
# ---------------------------------------------------------------------------


def _integer_rows(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int]
) -> list[list[int]]:
    result = []
    for row in rows:
        selected = [Fraction(row[column]) for column in columns]
        denominator = math.lcm(*(value.denominator for value in selected)) if selected else 1
        result.append(
            [value.numerator * (denominator // value.denominator) for value in selected]
        )
    return result


def _bareiss_rank_and_rows(
    rows: Sequence[Sequence[Fraction]], columns: Sequence[int]
) -> tuple[int, tuple[int, ...]]:
    matrix = _integer_rows(rows, columns)
    if not matrix or not columns:
        return 0, ()
    row_ids = list(range(len(matrix)))
    pivot_row = 0
    previous_pivot = 1
    for column in range(len(columns)):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        row_ids[pivot_row], row_ids[pivot] = row_ids[pivot], row_ids[pivot_row]
        pivot_value = matrix[pivot_row][column]
        for row in range(pivot_row + 1, len(matrix)):
            multiple = matrix[row][column]
            for later_column in range(column + 1, len(columns)):
                numerator = (
                    pivot_value * matrix[row][later_column]
                    - multiple * matrix[pivot_row][later_column]
                )
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise AssertionError("independent Bareiss division was not exact")
                matrix[row][later_column] = quotient
            matrix[row][column] = 0
        previous_pivot = pivot_value
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row, tuple(row_ids[:pivot_row])


def _rank(rows: Sequence[Sequence[Fraction]], column_count: int) -> int:
    return _bareiss_rank_and_rows(rows, tuple(range(column_count)))[0]


def _rref(rows: Sequence[Sequence[Fraction]], column_count: int):
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivots: list[int] = []
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
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiple = matrix[row][column]
            matrix[row] = [
                entry - multiple * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix, pivots


def _primitive(vector: Sequence[Fraction]) -> list[int]:
    denominator = math.lcm(*(Fraction(value).denominator for value in vector))
    integers = [
        Fraction(value).numerator * (denominator // Fraction(value).denominator)
        for value in vector
    ]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor:
        integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    return [-value for value in integers] if first < 0 else integers


def _nullspace(rows: Sequence[Sequence[Fraction]], column_count: int) -> list[list[int]]:
    matrix, pivots = _rref(rows, column_count)
    free = [column for column in range(column_count) if column not in pivots]
    basis: list[list[int]] = []
    for free_column in free:
        vector = [Fraction(0) for _ in range(column_count)]
        vector[free_column] = Fraction(1)
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -matrix[row][free_column]
        primitive = _primitive(vector)
        if any(_dot(row, primitive) for row in rows):
            raise AssertionError("independent nullspace construction has nonzero residual")
        basis.append(primitive)
    return basis


def _determinant(rows: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(rows)
    if size == 0:
        return Fraction(1)
    if any(len(row) != size for row in rows):
        raise AssertionError("independent determinant requires a square matrix")
    scales = []
    matrix = []
    for row in rows:
        values = [Fraction(value) for value in row]
        scale = math.lcm(*(value.denominator for value in values))
        scales.append(scale)
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
            multiple = matrix[row][diagonal]
            for column in range(diagonal + 1, size):
                numerator = (
                    pivot_value * matrix[row][column]
                    - multiple * matrix[diagonal][column]
                )
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise AssertionError("independent Bareiss determinant division was not exact")
                matrix[row][column] = quotient
            matrix[row][diagonal] = 0
        previous_pivot = pivot_value
    return Fraction(sign * matrix[-1][-1], math.prod(scales))


def _dot(row: Sequence[Fraction], vector: Sequence[int]) -> Fraction:
    return sum(
        (Fraction(entry) * coefficient for entry, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _selected_rows(rows: Sequence[Sequence[Fraction]], columns: Sequence[int]):
    return [[row[column] for column in columns] for row in rows]


# ---------------------------------------------------------------------------
# Independent matrix builders from the published ansatz definitions.
# ---------------------------------------------------------------------------


def _convolve(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    result = [Fraction(0) for _ in range(order + 1)]
    for left_order, left_value in enumerate(left[: order + 1]):
        if not left_value:
            continue
        for right_order, right_value in enumerate(right[: order + 1 - left_order]):
            if right_value:
                result[left_order + right_order] += left_value * right_value
    return tuple(result)


def _powers(series: Sequence[Fraction], maximum: int, order: int):
    values = [(Fraction(1),) + (Fraction(0),) * order]
    truncated = tuple(series[: order + 1])
    for _ in range(maximum):
        values.append(_convolve(values[-1], truncated, order))
    return values


def _valuation(series: Sequence[Fraction]) -> int:
    for order, coefficient in enumerate(series):
        if coefficient:
            return order
    raise AssertionError("series is identically zero")


def _algebraic_matrix(series: Sequence[Fraction], degree_x: int, degree_f: int):
    count = len(series)
    monomials = [
        (power_x, power_f)
        for power_f in range(degree_f + 1)
        for power_x in range(degree_x + 1)
    ]
    powers = _powers(series, degree_f, count - 1)
    rows = [
        [
            powers[power_f][order - power_x] if order >= power_x else Fraction(0)
            for power_x, power_f in monomials
        ]
        for order in range(count)
    ]
    leading = _valuation(series)
    valuations = [power_x + power_f * leading for power_x, power_f in monomials]
    return rows, monomials, valuations


def _first_derivative_matrix(
    series: Sequence[Fraction], degree_x: int, degree_f: int, degree_f_prime: int
):
    count = len(series) - 1
    truncated = tuple(series[:count])
    derivative = tuple(Fraction(order + 1) * series[order + 1] for order in range(count))
    f_powers = _powers(truncated, degree_f, count - 1)
    derivative_powers = _powers(derivative, degree_f_prime, count - 1)
    monomials = [
        (power_x, power_f, power_f_prime)
        for power_f_prime in range(degree_f_prime + 1)
        for power_f in range(degree_f + 1)
        for power_x in range(degree_x + 1)
    ]
    columns = []
    for power_x, power_f, power_f_prime in monomials:
        product = _convolve(f_powers[power_f], derivative_powers[power_f_prime], count - 1)
        columns.append(
            [product[order - power_x] if order >= power_x else Fraction(0) for order in range(count)]
        )
    rows = [[column[order] for column in columns] for order in range(count)]
    f_valuation = _valuation(series)
    valuations = [
        power_x + power_f * f_valuation + power_f_prime * (f_valuation - 1)
        for power_x, power_f, power_f_prime in monomials
    ]
    return rows, monomials, valuations


def _euler_matrix(series: Sequence[Fraction], order: int, degree: int):
    count = len(series)
    monomials = [
        (power_t, theta_power)
        for theta_power in range(order + 1)
        for power_t in range(degree + 1)
    ]
    rows = [
        [
            Fraction(coefficient_order - power_t) ** theta_power
            * series[coefficient_order - power_t]
            if coefficient_order >= power_t
            else Fraction(0)
            for power_t, theta_power in monomials
        ]
        for coefficient_order in range(count)
    ]
    leading = _valuation(series)
    valuations = [power_t + leading for power_t, _ in monomials]
    return rows, monomials, valuations


def _mahler_matrix(series: Sequence[Fraction], degree_t: int, dependent_total_degree: int):
    count = len(series)
    composed = tuple(
        series[order // 2] if order % 2 == 0 else Fraction(0) for order in range(count)
    )
    monomials = [
        (power_t, power_f, power_composed)
        for power_composed in range(dependent_total_degree + 1)
        for power_f in range(dependent_total_degree + 1 - power_composed)
        for power_t in range(degree_t + 1)
    ]
    f_powers = _powers(series, dependent_total_degree, count - 1)
    composed_powers = _powers(composed, dependent_total_degree, count - 1)
    columns = []
    for power_t, power_f, power_composed in monomials:
        product = _convolve(f_powers[power_f], composed_powers[power_composed], count - 1)
        columns.append(
            [product[order - power_t] if order >= power_t else Fraction(0) for order in range(count)]
        )
    rows = [[column[order] for column in columns] for order in range(count)]
    leading = _valuation(series)
    valuations = [
        power_t + power_f * leading + power_composed * (2 * leading)
        for power_t, power_f, power_composed in monomials
    ]
    return rows, monomials, valuations


def _build_matrix(family: str, series: Sequence[Fraction], parameters: dict):
    if family == "algebraic":
        return _algebraic_matrix(series, parameters["degree_x"], parameters["degree_f"])
    if family in {"first_order_differential_algebraic", "differential_algebraic_order_1"}:
        return _first_derivative_matrix(
            series,
            parameters["degree_x"],
            parameters["degree_f"],
            parameters["degree_f_prime"],
        )
    if family == "linear_euler_ode":
        return _euler_matrix(series, parameters["order"], parameters["polynomial_degree"])
    if family == "mahler":
        return _mahler_matrix(
            series,
            parameters["degree_t"],
            parameters["dependent_total_degree"],
        )
    raise AssertionError(f"unknown family {family}")


def _load_raw_series() -> dict[str, tuple[Fraction, ...]]:
    ht_payload = json.loads(HT_INPUT.read_text(encoding="utf-8"))
    lt_payload = json.loads(LT_INPUT.read_text(encoding="utf-8"))
    ht = tuple(Fraction(value) for value in ht_payload["data"]["series"]["coefficients"])
    lt = tuple(Fraction(value) for value in lt_payload["data"]["coefficients"])
    _check(
        "raw input sizes and parity",
        len(ht) == 29
        and len(lt) == 33
        and not any(ht[order] for order in range(1, len(ht), 2))
        and not any(lt[order] for order in range(1, len(lt), 2)),
        "HT is v^0..v^28 and LT is x^0..x^32 with exact odd zeros",
    )
    return {"HT": ht, "LT": lt}


# ---------------------------------------------------------------------------
# Certificate and survivor replay.
# ---------------------------------------------------------------------------


def _certificate(payload: dict, pointer: str) -> tuple[str, dict]:
    prefix = "/data/certificates/"
    if not pointer.startswith(prefix):
        raise AssertionError(f"bad certificate pointer {pointer}")
    certificate_id = pointer.removeprefix(prefix)
    certificate = payload["data"]["certificates"].get(certificate_id)
    if certificate is None:
        raise AssertionError(f"missing certificate {certificate_id}")
    return certificate_id, certificate


def _replay_minor(
    certificate: dict,
    rows: Sequence[Sequence[Fraction]],
    *,
    expected_columns: Sequence[int] | None = None,
    maximum_row_exclusive: int | None = None,
) -> None:
    columns = certificate["column_indices"]
    row_orders = certificate["row_orders"]
    if expected_columns is not None and columns != list(expected_columns):
        raise AssertionError("certificate columns differ from the independently expected columns")
    if len(columns) != len(row_orders) or len(set(columns)) != len(columns):
        raise AssertionError("certificate is not square with distinct columns")
    if any(order < 0 or order >= len(rows) for order in row_orders):
        raise AssertionError("certificate uses an unavailable row")
    if maximum_row_exclusive is not None and any(order >= maximum_row_exclusive for order in row_orders):
        raise AssertionError("training-prefix certificate used a holdout row")
    minor = [
        [rows[row][column] for column in columns]
        for row in row_orders
    ]
    determinant = _determinant(minor)
    if determinant == 0 or str(determinant) != certificate["determinant"]:
        raise AssertionError("independent exact maximal-minor replay disagrees")


def _verify_basis(
    rows: Sequence[Sequence[Fraction]], basis: Sequence[Sequence[int]], column_count: int
) -> None:
    if any(len(vector) != column_count for vector in basis):
        raise AssertionError("kernel vector has the wrong dimension")
    if any(_dot(row, vector) for vector in basis for row in rows):
        raise AssertionError("stored kernel vector has a nonzero exact residual")
    if _rank(basis, column_count) != len(basis):
        raise AssertionError("stored kernel basis is not independent")


def _support_modulus(rows: Sequence[Sequence[Fraction]], valuations: Sequence[int]) -> int | None:
    for modulus in range(2, len(rows) + 1):
        if all(
            all(
                (not entry) or order % modulus == valuations[column] % modulus
                for column, entry in enumerate(row)
            )
            for order, row in enumerate(rows)
        ):
            return modulus
    return None


def _replay_survivor_support(
    payload: dict,
    series_name: str,
    family: str,
    cell: dict,
    rows: Sequence[Sequence[Fraction]],
    monomials: Sequence[tuple[int, ...]],
    valuations: Sequence[int],
) -> int:
    support = cell["support_analysis"]
    if support["classification"] != "TRUNCATION_UNOBSERVABLE":
        raise AssertionError("an observable finite-prefix survivor is unexplained")
    pointer = support["support_certificate_pointer"]
    _, certificate = _certificate(payload, pointer)
    if (
        certificate["series"] != series_name
        or certificate["ansatz_family"] != family
        or certificate["parameters"] != cell["parameters"]
    ):
        raise AssertionError("support certificate context does not match its cell")

    column_count = len(monomials)
    zero_columns = [
        column for column in range(column_count) if all(row[column] == 0 for row in rows)
    ]
    full_basis = cell["full_kernel_basis"]
    _verify_basis(rows, full_basis, column_count)
    expected_full_nullity = column_count - _rank(rows, column_count)
    if expected_full_nullity != len(full_basis):
        raise AssertionError("stored full survivor kernel has the wrong dimension")
    residuals = [
        [str(_dot(row, vector)) for row in rows]
        for vector in full_basis
    ]
    if residuals != cell["full_kernel_residuals"]:
        raise AssertionError("stored full-kernel residual table disagrees with raw coefficients")
    if zero_columns != cell["zero_column_indices_at_known_order"]:
        raise AssertionError("stored zero-column support set disagrees")

    components = support["support_components"]
    if components and all(component.get("type") == "ZERO_COLUMN" for component in components):
        visible_columns = tuple(
            column for column in range(column_count) if column not in zero_columns
        )
        if certificate["phase"] != "observable_columns":
            raise AssertionError("zero-column survivor has the wrong certificate phase")
        if certificate.get("zero_column_indices") != zero_columns:
            raise AssertionError("zero-column certificate omits or changes a zero column")
        _replay_minor(
            certificate,
            rows,
            expected_columns=visible_columns,
        )
        if _rank(_selected_rows(rows, visible_columns), len(visible_columns)) != len(visible_columns):
            raise AssertionError("visible columns were not independently full rank")
        expected_components = [
            {
                "type": "ZERO_COLUMN",
                "column_index": column,
                "monomial": list(monomials[column]),
                "first_unfixed_supported_order": valuations[column],
                "last_available_order": len(rows) - 1,
            }
            for column in zero_columns
        ]
        if components != expected_components:
            raise AssertionError("zero-column support components disagree with raw matrix")
        expected_unobservable = [
            {
                "column_index": column,
                "monomial": list(monomials[column]),
                "lowest_possible_relation_order": valuations[column],
                "last_available_relation_order": len(rows) - 1,
            }
            for column in zero_columns
        ]
        if certificate.get("unobservable_monomials") != expected_unobservable:
            raise AssertionError("zero-column certificate monomials disagree")
        first_testable = min(valuations[column] for column in zero_columns)
    else:
        if certificate["type"] != "TRUNCATION_UNOBSERVABLE_SUPPORT_KERNEL":
            raise AssertionError("support-progression survivor has the wrong certificate type")
        if certificate["phase"] != "all_available_equations":
            raise AssertionError("support-progression survivor has the wrong phase")
        _verify_basis(rows, certificate["full_kernel_basis"], column_count)
        modulus = _support_modulus(rows, valuations)
        if modulus is None:
            raise AssertionError("stored support progression lacks an independently valid modulus")
        expected_components = []
        total_nullity = 0
        for residue in range(modulus):
            columns = tuple(
                column for column, valuation in enumerate(valuations)
                if valuation % modulus == residue
            )
            if not columns:
                continue
            populated_orders = tuple(
                order for order in range(len(rows)) if order % modulus == residue
            )
            component_rows = [rows[order] for order in populated_orders]
            rank, selected_local_rows = _bareiss_rank_and_rows(component_rows, columns)
            nullity = len(columns) - rank
            total_nullity += nullity
            if nullity == 0:
                selected_orders = [populated_orders[index] for index in selected_local_rows]
                determinant = _determinant(
                    [[rows[order][column] for column in columns] for order in selected_orders]
                )
                if determinant == 0:
                    raise AssertionError("support component full-rank determinant vanished")
                determinant_string = str(determinant)
            else:
                selected_orders = []
                determinant_string = None
            first_unfixed = residue if not populated_orders else populated_orders[-1] + modulus
            expected_components.append(
                {
                    "support_modulus": modulus,
                    "support_residue": residue,
                    "column_indices": list(columns),
                    "monomials": [list(monomials[column]) for column in columns],
                    "known_populated_orders": list(populated_orders),
                    "known_rank": rank,
                    "known_nullity": nullity,
                    "certifying_row_orders_when_full_rank": selected_orders,
                    "certifying_minor_determinant_when_full_rank": determinant_string,
                    "first_unfixed_supported_order": first_unfixed,
                    "last_available_order": len(rows) - 1,
                }
            )
        if components != expected_components or certificate["support_components"] != expected_components:
            raise AssertionError("support progression components disagree with raw matrix")
        if total_nullity != len(full_basis):
            raise AssertionError("support components do not account for the full kernel")
        first_testable = min(
            component["first_unfixed_supported_order"]
            for component in expected_components
            if component["known_nullity"] > 0
        )

    if first_testable != support["first_testable_order"] or first_testable != cell["first_testable_order"]:
        raise AssertionError("survivor first testable order disagrees with support replay")
    if first_testable < len(rows):
        raise AssertionError("purported truncation-unobservable survivor is already testable")
    return first_testable


def _replay_cell(
    payload: dict,
    series_name: str,
    family: str,
    cell: dict,
    series: Sequence[Fraction],
) -> tuple[dict, bool]:
    rows, monomials, valuations = _build_matrix(family, series, cell["parameters"])
    column_count = len(monomials)
    training_count = column_count + 2
    if cell["monomials"] != [list(monomial) for monomial in monomials]:
        raise AssertionError("stored monomial ordering differs from the independent builder")
    if cell["formal_monomial_valuations"] != valuations:
        raise AssertionError("stored formal valuations differ from the independent builder")
    if (
        cell["unknown_count"] != column_count
        or cell["available_equation_count"] != len(rows)
        or cell["training_equation_count"] != training_count
        or cell["training_orders"] != [0, training_count - 1]
        or cell["holdout_orders"]
        != ([training_count, len(rows) - 1] if training_count < len(rows) else [])
    ):
        raise AssertionError("strict U+2 budget metadata disagrees with the raw matrix")

    strict_training_rows = rows[:training_count]
    strict_rank = _rank(strict_training_rows, column_count)
    strict_basis = _nullspace(strict_training_rows, column_count)
    if strict_rank != cell["training_rank"] or len(strict_basis) != cell["training_nullity"]:
        raise AssertionError("independent strict-prefix rank or nullity disagrees")
    short_rank = _rank(rows[: column_count + 1], column_count)
    split_sensitive = short_rank != strict_rank

    if cell["verdict"] == "NO_RELATION_AT_BUDGET":
        if strict_rank != column_count or _rank(rows, column_count) != column_count:
            raise AssertionError("stored no-relation cell is not independently full rank")
        _, certificate = _certificate(payload, cell["certificate_pointer"])
        if (
            certificate["series"] != series_name
            or certificate["ansatz_family"] != family
            or certificate["parameters"] != cell["parameters"]
            or certificate["phase"] != "training_prefix"
        ):
            raise AssertionError("no-relation certificate context or phase disagrees")
        _replay_minor(
            certificate,
            rows,
            expected_columns=tuple(range(column_count)),
            maximum_row_exclusive=training_count,
        )
        return {"verdict": cell["verdict"], "split_sensitive": split_sensitive}, split_sensitive

    stored_training_basis = cell["training_kernel_basis"]
    if len(stored_training_basis) != len(strict_basis):
        raise AssertionError("stored training kernel dimension disagrees")
    _verify_basis(strict_training_rows, stored_training_basis, column_count)
    residual_constraints: list[list[Fraction]] = []
    expected_holdout_checks = []
    first_refuting_order = None
    first_refuting_residual = None
    for order in range(training_count, len(rows)):
        residuals = [_dot(rows[order], vector) for vector in stored_training_basis]
        residual_constraints.append(residuals)
        residual_rank = _rank(residual_constraints, len(stored_training_basis))
        surviving_dimension = len(stored_training_basis) - residual_rank
        expected_holdout_checks.append(
            {
                "order": order,
                "residual_on_training_kernel_basis": [str(value) for value in residuals],
                "surviving_kernel_dimension": surviving_dimension,
            }
        )
        if surviving_dimension == 0 and first_refuting_order is None:
            first_refuting_order = order
            first_refuting_residual = [str(value) for value in residuals]
    if expected_holdout_checks != cell["holdout_checks"]:
        raise AssertionError("full exact holdout ladder disagrees with raw coefficients")
    full_rank = _rank(rows, column_count)
    full_nullity = column_count - full_rank
    if full_rank != cell["full_rank"] or full_nullity != cell["full_nullity"]:
        raise AssertionError("full raw matrix rank or nullity disagrees")
    if first_refuting_order != cell["first_refuting_holdout_order"]:
        raise AssertionError("first exact refuting holdout order disagrees")

    if cell["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT":
        if full_rank != column_count or first_refuting_residual is None:
            raise AssertionError("stored holdout refutation is not independently complete")
        if first_refuting_residual != cell["first_refuting_residual_on_training_kernel_basis"]:
            raise AssertionError("stored first refuting residual disagrees exactly")
        _, certificate = _certificate(payload, cell["certificate_pointer"])
        if (
            certificate["series"] != series_name
            or certificate["ansatz_family"] != family
            or certificate["parameters"] != cell["parameters"]
            or certificate["phase"] != "all_available_equations"
        ):
            raise AssertionError("holdout-refutation certificate context or phase disagrees")
        _replay_minor(
            certificate,
            rows,
            expected_columns=tuple(range(column_count)),
        )
        return {"verdict": cell["verdict"], "split_sensitive": split_sensitive}, split_sensitive

    if cell["verdict"] != "CANDIDATE_SURVIVES" or full_nullity == 0:
        raise AssertionError("forbidden or inconsistent verdict")
    first_testable_order = _replay_survivor_support(
        payload, series_name, family, cell, rows, monomials, valuations
    )
    return {
        "verdict": cell["verdict"],
        "split_sensitive": split_sensitive,
        "first_testable_order": first_testable_order,
        "support_classification": cell["support_analysis"]["classification"],
    }, split_sensitive


# ---------------------------------------------------------------------------
# Legacy e43 LT survivor replay in its original reduced u=x^2 convention.
# ---------------------------------------------------------------------------


def _replay_legacy_lt(payload: dict, lt: Sequence[Fraction]) -> None:
    legacy = payload["data"]["legacy_lt_survivor_reaudit"]
    rows = legacy["rows"]
    if legacy["stored_row_count"] != 18 or legacy["replayed_row_count"] != 18 or len(rows) != 18:
        raise AssertionError("legacy LT audit does not contain exactly 18 rows")
    reduced = tuple(lt[order] for order in range(0, len(lt), 2))
    for index, ((relation_class, parameters), stored) in enumerate(zip(LEGACY_LT_ROWS, rows, strict=True)):
        if (
            stored["legacy_row_index"] != index
            or stored["relation_class"] != relation_class
            or stored["degrees"] != parameters
        ):
            raise AssertionError("legacy LT row descriptor differs from the fixed e43 survivor list")
        matrix_rows, monomials, valuations = _build_matrix(relation_class, reduced, parameters)
        column_count = len(monomials)
        rank = _rank(matrix_rows, column_count)
        basis = _nullspace(matrix_rows, column_count)
        zero_columns = [
            column for column in range(column_count) if all(row[column] == 0 for row in matrix_rows)
        ]
        if (
            stored["full_rank"] != rank
            or stored["full_nullity"] != len(basis)
            or stored["zero_column_indices_at_known_order"] != zero_columns
            or stored["support_classification"] != "TRUNCATION_UNOBSERVABLE"
            or not stored["unchanged"]
        ):
            raise AssertionError("legacy LT raw reduced-series replay disagrees")
        _verify_basis(matrix_rows, basis, column_count)
        visible_columns = tuple(
            column for column in range(column_count) if column not in zero_columns
        )
        if len(basis) != len(zero_columns):
            raise AssertionError("legacy full kernel is not exactly the zero-column span")
        if _rank(_selected_rows(matrix_rows, visible_columns), len(visible_columns)) != len(visible_columns):
            raise AssertionError("legacy visible columns are not full rank")
        first_testable = min(valuations[column] for column in zero_columns)
        if first_testable != stored["first_testable_order_in_u"] or first_testable < len(matrix_rows):
            raise AssertionError("legacy first testable order disagrees")
        _, certificate = _certificate(payload, stored["support_certificate_pointer"])
        if (
            certificate["series"] != "LT"
            or certificate["ansatz_family"] != f"legacy_e43_{relation_class}"
            or certificate["parameters"] != parameters
            or certificate["legacy_row_index"] != index
            or certificate["phase"] != "observable_columns"
            or certificate.get("zero_column_indices") != zero_columns
        ):
            raise AssertionError("legacy zero-column certificate context disagrees")
        _replay_minor(certificate, matrix_rows, expected_columns=visible_columns)
    if legacy["changed_rows"] != []:
        raise AssertionError("legacy audit claims changed rows despite independent replay")
    if legacy["newly_available_equation_orders"] != []:
        raise AssertionError("identical x^32 input unexpectedly reports new legacy equations")


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    raw_series = _load_raw_series()
    observed_counts = {verdict: 0 for verdict in ALLOWED_VERDICTS}
    replayed_survivors = []
    split_sensitive_witnesses: dict[str, dict] = {}
    evaluated = 0

    for series_name in ("HT", "LT"):
        families = payload["data"]["frontiers"][series_name]["families"]
        if tuple(families) != FAMILIES:
            raise AssertionError("frontier family ordering or vocabulary changed")
        for family in FAMILIES:
            family_data = families[family]
            family_counts = {verdict: 0 for verdict in ALLOWED_VERDICTS}
            for cell in family_data["cells"]:
                if cell["verdict"] not in ALLOWED_VERDICTS:
                    raise AssertionError("evaluated cell uses a forbidden verdict")
                replay, split_sensitive = _replay_cell(
                    payload, series_name, family, cell, raw_series[series_name]
                )
                observed_counts[replay["verdict"]] += 1
                family_counts[replay["verdict"]] += 1
                evaluated += 1
                if split_sensitive and family not in split_sensitive_witnesses:
                    split_sensitive_witnesses[family] = {
                        "series": series_name,
                        "parameters": cell["parameters"],
                        "strict_training_rank": cell["training_rank"],
                    }
                if replay["verdict"] == "CANDIDATE_SURVIVES":
                    replayed_survivors.append(
                        {
                            "claim_tag": "[COMPUTATION]",
                            "series": series_name,
                            "ansatz_family": family,
                            "parameters": cell["parameters"],
                            "full_nullity": cell["full_nullity"],
                            "support_classification": replay["support_classification"],
                            "first_testable_order": replay["first_testable_order"],
                        }
                    )
            if family_counts != family_data["verdict_counts"]:
                raise AssertionError(f"{series_name} {family} verdict counts disagree")

    summary = payload["data"]["summary"]
    _check(
        "exhaustive raw-series replay",
        evaluated == summary["evaluated_cell_count"] == 633
        and observed_counts == summary["verdict_counts"]
        and set(observed_counts) == ALLOWED_VERDICTS,
        f"replayed all {evaluated} evaluated HT/LT cells with exact counts {observed_counts}",
    )
    _check(
        "exact verdict vocabulary",
        summary["verdict_vocabulary"] == sorted(ALLOWED_VERDICTS),
        "only NO_RELATION_AT_BUDGET, CANDIDATE_REFUTED_BY_HOLDOUT, and CANDIDATE_SURVIVES occur",
    )
    _check(
        "split-sensitive U+2 witnesses",
        set(split_sensitive_witnesses) == set(FAMILIES),
        f"independent rank changes from 0..U to 0..U+1 for {split_sensitive_witnesses}",
    )
    _check(
        "survivor support replay",
        replayed_survivors == summary["survivors"]
        and summary["truncation_unobservable_survivor_count"] == len(replayed_survivors)
        and summary["unexplained_survivor_count"] == 0
        and summary["unexplained_survivors"] == [],
        f"all {len(replayed_survivors)} full kernels have independently replayed support certificates",
    )
    _replay_legacy_lt(payload, raw_series["LT"])
    _check(
        "all 18 legacy LT survivors",
        True,
        "rebuilt reduced-u matrices, zero columns, visible minors, and first testable orders",
    )
    print("PASS test_series_grid")


if __name__ == "__main__":
    main()
