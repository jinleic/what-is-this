"""Independent exact checks for the second-order series frontier artifact."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "series" / "frontier2.json"
HT = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
LT = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _load(path: Path) -> tuple[Fraction, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return tuple(Fraction(value) for value in payload["data"]["coefficients"])


def _differentiate(series: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(Fraction(order) * series[order] for order in range(1, len(series)))


def _convolve(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    result = [Fraction(0) for _ in range(order + 1)]
    for i, first in enumerate(left[: order + 1]):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                result[i + j] += first * second
    return tuple(result)


def _powers(
    series: Sequence[Fraction], maximum: int, order: int
) -> tuple[tuple[Fraction, ...], ...]:
    padded = tuple(series[: order + 1])
    if len(padded) < order + 1:
        padded += (Fraction(0),) * (order + 1 - len(padded))
    result = [(Fraction(1),) + (Fraction(0),) * order]
    for _ in range(maximum):
        result.append(_convolve(result[-1], padded, order))
    return tuple(result)


def _determinant(rows: Sequence[Sequence[Fraction]]) -> Fraction:
    """Independent exact Fraction Gaussian-elimination determinant."""

    matrix = [[Fraction(value) for value in row] for row in rows]
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    result = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            result = -result
        pivot_value = matrix[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            if not matrix[row][column]:
                continue
            multiplier = matrix[row][column] / pivot_value
            for later_column in range(column + 1, size):
                matrix[row][later_column] -= multiplier * matrix[column][later_column]
            matrix[row][column] = Fraction(0)
    return result


def _rank(rows: Sequence[Sequence[Fraction]], column_count: int) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        for row in range(pivot_row + 1, len(matrix)):
            if not matrix[row][column]:
                continue
            multiplier = matrix[row][column] / pivot_value
            for later_column in range(column + 1, column_count):
                matrix[row][later_column] -= multiplier * matrix[pivot_row][later_column]
            matrix[row][column] = Fraction(0)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def _jet_rows(
    series: tuple[Fraction, ...], monomials: Sequence[Sequence[int]]
) -> list[list[Fraction]]:
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
    return [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]


def _linear_rows(
    series: tuple[Fraction, ...], monomials: Sequence[Sequence[int]]
) -> list[list[Fraction]]:
    maximum_theta = max(theta_power for _, theta_power in monomials)
    theta = tuple(
        tuple(Fraction(index**power) * coefficient for index, coefficient in enumerate(series))
        for power in range(maximum_theta + 1)
    )
    columns = [
        [
            theta[theta_power][order - degree_t]
            if order >= degree_t
            else Fraction(0)
            for order in range(len(series))
        ]
        for degree_t, theta_power in monomials
    ]
    return [
        [columns[column][order] for column in range(len(columns))]
        for order in range(len(series))
    ]


def _mahler_rows(
    series: tuple[Fraction, ...], monomials: Sequence[Sequence[int]]
) -> list[list[Fraction]]:
    equation_count = len(series)
    maximum = max(max(degree_f, degree_composed) for _, degree_f, degree_composed in monomials)
    composed = tuple(
        series[order // 2] if order % 2 == 0 else Fraction(0)
        for order in range(equation_count)
    )
    f_powers = _powers(series, maximum, equation_count - 1)
    composed_powers = _powers(composed, maximum, equation_count - 1)
    columns = []
    for degree_t, degree_f, degree_composed in monomials:
        product = _convolve(
            f_powers[degree_f], composed_powers[degree_composed], equation_count - 1
        )
        columns.append(
            [
                product[order - degree_t] if order >= degree_t else Fraction(0)
                for order in range(equation_count)
            ]
        )
    return [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]


def _table(result: dict, family: str, series: str) -> list[dict]:
    return result["data"]["frontiers"][family][series]["frontier"]


def _row(result: dict, family: str, series: str, parameters: dict) -> dict:
    return next(
        row
        for row in _table(result, family, series)
        if row["parameters"] == parameters
    )


def _matrix_for_row(
    series: tuple[Fraction, ...], row: dict
) -> list[list[Fraction]]:
    family = row["ansatz_family"]
    if family.startswith("second_order_"):
        return _jet_rows(series, row["monomials"])
    if family == "linear_euler_ode":
        return _linear_rows(series, row["monomials"])
    if family == "mahler_bidegree":
        return _mahler_rows(series, row["monomials"])
    raise ValueError(f"unsupported family {family}")


def _certificate(result: dict, row: dict) -> dict:
    certificate_id = row["certificate_pointer"].rsplit("/", 1)[-1]
    return result["data"]["certificates"][certificate_id]


def _recompute_minor(
    result: dict, series: tuple[Fraction, ...], row: dict
) -> tuple[Fraction, Fraction]:
    matrix = _matrix_for_row(series, row)
    certificate = _certificate(result, row)
    minor = [
        [matrix[order][column] for column in certificate["column_indices"]]
        for order in certificate["row_orders"]
    ]
    return _determinant(minor), Fraction(certificate["determinant"])


def main() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    ht = _load(HT)
    lt = _load(LT)

    sample_specs = (
        (
            "second_order_differential_algebraic",
            "HT",
            {"degree_t": 1, "jet_total_degree": 1},
            ht,
        ),
        (
            "second_order_differential_algebraic",
            "LT",
            {"degree_t": 6, "jet_total_degree": 1},
            lt,
        ),
        (
            "linear_euler_ode_order_at_most_2",
            "HT",
            {"order": 2, "polynomial_degree": 5},
            ht,
        ),
        (
            "mahler_F_t2_vs_F_t",
            "LT",
            {"degree_t": 8, "dependent_total_degree": 1},
            lt,
        ),
    )
    minor_checks = []
    for family, series_name, parameters, series in sample_specs:
        row = _row(result, family, series_name, parameters)
        if row["verdict"] != "NO_RELATION_AT_BUDGET":
            raise AssertionError("minor sample is not a stored NO_RELATION row")
        recomputed, stored = _recompute_minor(result, series, row)
        minor_checks.append(recomputed == stored and recomputed != 0)
    _check(
        "sample stored NO_RELATION minors",
        all(minor_checks),
        "four independently reconstructed rational maximal minors equal the stored nonzero determinants",
    )

    holdout = _row(
        result,
        "second_order_differential_algebraic",
        "HT",
        {"degree_t": 0, "jet_total_degree": 2},
    )
    holdout_rows = _matrix_for_row(ht, holdout)
    training_count = holdout["training_equation_count"]
    unknown_count = holdout["unknown_count"]
    training_rank = _rank(holdout_rows[:training_count], unknown_count)
    full_rank = _rank(holdout_rows, unknown_count)
    first_refuting = None
    for stop in range(training_count + 1, len(holdout_rows) + 1):
        if _rank(holdout_rows[:stop], unknown_count) == unknown_count:
            first_refuting = stop - 1
            break
    _check(
        "fresh second-order holdout refutation",
        holdout["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
        and training_rank < unknown_count
        and full_rank == unknown_count
        and first_refuting == holdout["first_refuting_holdout_order"] == 12,
        (
            f"fresh Fraction ranks are training={training_rank}, full={full_rank}; "
            f"the first refuting withheld order is {first_refuting}"
        ),
    )

    support_row = _row(
        result,
        "linear_euler_ode_order_at_most_2",
        "HT",
        {"order": 2, "polynomial_degree": 6},
    )
    support_matrix = _matrix_for_row(ht, support_row)
    support_certificate = _certificate(result, support_row)
    vector = tuple(support_certificate["full_kernel_basis"][0])
    residuals = [
        sum(
            (value * coefficient for value, coefficient in zip(row, vector, strict=True)),
            Fraction(0),
        )
        for row in support_matrix
    ]
    supported_component = next(
        component
        for component in support_certificate["support_components"]
        if component["known_nullity"] > 0
    )
    _check(
        "truncation-unobservable support certificate",
        support_row["verdict"] == "TRUNCATION_UNOBSERVABLE"
        and all(residual == 0 for residual in residuals)
        and supported_component["first_unfixed_supported_order"] == 24
        and supported_component["last_available_order"] == 22,
        "the exact HT order-2 degree-6 kernel vanishes through v^22 and its next parity-supported order is v^24",
    )

    provenance = result["provenance"]
    inputs = result["data"]["input_series"]
    summaries = result["data"]["frontier_summary"]
    audit = result["data"]["survivor_audit"]
    _check(
        "artifact envelope and provenance",
        set(result) == {"provenance", "data", "checks"}
        and set(provenance) >= {"script", "generated_utc", "interpreter", "method"}
        and provenance["script"] == "experiments/e56_series_frontier2.py"
        and provenance["interpreter"] == ".venv/bin/python"
        and inputs["HT"]["sha256"] == hashlib.sha256(HT.read_bytes()).hexdigest()
        and inputs["LT"]["sha256"] == hashlib.sha256(LT.read_bytes()).hexdigest()
        and all(check["passed"] for check in result["checks"]),
        "required JSON envelope, exact input hashes, provenance, and embedded checks are valid",
    )
    _check(
        "complete frontiers and zero unexplained survivors",
        {(row["ansatz_family"], row["series"]) for row in summaries}
        == {
            (family, series)
            for family in (
                "second_order_differential_algebraic",
                "linear_euler_ode_order_at_most_2",
                "mahler_F_t2_vs_F_t",
            )
            for series in ("HT", "LT")
        }
        and sum(row["frontier_row_count"] for row in summaries) == 138
        and audit["unexplained_survivor_count"] == 0
        and not audit["parent_review_rows"],
        "all six HT/LT family tables are present (138 rows) and every kernel is refuted or certified unobservable",
    )
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_series_frontier2: {error}")
        raise
