"""Exact finite-prefix structure certificates for the extended 3D Ising series.

This experiment performs two deliberately finite searches:

* algebraic relations P(t, f(t)) = 0 in rectangular bidegree budgets;
* first-order differential-algebraic relations Q(t, f(t), f'(t)) = 0
  in rectangular tridegree budgets.

It also extends the independent interlayer c2 comparison from v^6 to v^8.
Run from the repository root with

    .venv/bin/python experiments/e43_series_structure.py
"""

from __future__ import annotations

import hashlib
import json
import math
import signal
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Callable, Iterable, Sequence

from sympy import Matrix, Rational

from ising.interlayer import (
    anisotropic_box_even_subgraph,
    anisotropic_flm_c2_series,
)
from ising.lattices import square

SCRIPT = "experiments/e43_series_structure.py"
ROOT = Path(__file__).resolve().parents[1]
HT_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
LT_PATH = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
RESULT_PATH = ROOT / "results" / "series" / "structure_certificates.json"
SOLVE_TIMEOUT_SECONDS = 1800
INTERLAYER_V_ORDER = 8

RationalSeries = tuple[Fraction, ...]


class ExactSolveTimeout(TimeoutError):
    """Raised when one exact linear-algebra operation exceeds its hard limit."""


@contextmanager
def _solve_deadline(seconds: int = SOLVE_TIMEOUT_SECONDS):
    """Apply a real-time SIGALRM deadline to one exact solve."""

    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_timer = signal.getitimer(signal.ITIMER_REAL)

    def _raise_timeout(_signum, _frame):
        raise ExactSolveTimeout(f"exact linear solve exceeded {seconds} seconds")

    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_timer[0] > 0:
            signal.setitimer(signal.ITIMER_REAL, *previous_timer)


def _sympy(value: Fraction | int) -> Rational:
    exact = Fraction(value)
    return Rational(exact.numerator, exact.denominator)


def _fraction(value) -> Fraction:
    exact = Rational(value)
    return Fraction(int(exact.p), int(exact.q))


def _fraction_strings(values: Iterable[Fraction | int]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _primitive_vector(vector: Sequence) -> tuple[int, ...]:
    fractions = [_fraction(value) for value in vector]
    denominator = math.lcm(*(value.denominator for value in fractions))
    integers = [
        value.numerator * (denominator // value.denominator) for value in fractions
    ]
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
        (Fraction(entry) * coefficient for entry, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _exact_rank(rows: Sequence[Sequence[Fraction]], column_count: int) -> int:
    if not rows:
        return 0
    with _solve_deadline():
        matrix = Matrix([[_sympy(value) for value in row] for row in rows])
        return int(matrix.rank())


def _matrix_analysis(
    rows: Sequence[Sequence[Fraction]], column_count: int
) -> dict:
    """Return an exact nullspace and, at full rank, a nonzero-minor witness."""

    with _solve_deadline():
        matrix = Matrix([[_sympy(value) for value in row] for row in rows])
        nullspace = matrix.nullspace()
        rank = column_count - len(nullspace)
        basis = [_primitive_vector(tuple(vector)) for vector in nullspace]
        certificate = None
        if rank == column_count:
            independent_rows = tuple(int(index) for index in matrix.T.rref()[1])
            if len(independent_rows) != column_count:
                raise AssertionError("full column rank did not produce enough pivot rows")
            minor = matrix.extract(independent_rows, tuple(range(column_count)))
            determinant = _fraction(minor.det())
            if determinant == 0:
                raise AssertionError("claimed full-rank minor has zero determinant")
            certificate = {
                "type": "NONZERO_MAXIMAL_MINOR",
                "row_orders": list(independent_rows),
                "column_indices": list(range(column_count)),
                "determinant": str(determinant),
            }
    return {
        "rank": rank,
        "nullity": len(basis),
        "kernel_basis": basis,
        "full_rank_certificate": certificate,
    }


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
    one = (Fraction(1),) + (Fraction(0),) * order
    values = [one]
    for _ in range(maximum):
        values.append(_convolve(values[-1], series, order))
    return tuple(values)


def _valuation(series: Sequence[Fraction]) -> int:
    for degree, coefficient in enumerate(series):
        if coefficient:
            return degree
    raise ValueError("the supplied finite series is identically zero")


def _load_reduced_even_series(path: Path) -> tuple[RationalSeries, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    coefficients = tuple(Fraction(value) for value in payload["data"]["coefficients"])
    if any(coefficients[degree] for degree in range(1, len(coefficients), 2)):
        raise AssertionError(f"{path.name} has a nonzero odd-power coefficient")
    reduced = tuple(coefficients[degree] for degree in range(0, len(coefficients), 2))
    metadata = {
        "artifact": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "original_variable": payload["data"]["variable"],
        "original_achieved_order": payload["data"]["achieved_order"],
        "reduced_coefficient_count": len(reduced),
        "reduced_coefficients": _fraction_strings(reduced),
    }
    return reduced, metadata


def _algebraic_matrix(
    series: RationalSeries, degree_x: int, degree_f: int
) -> tuple[list[list[Fraction]], list[tuple[int, int]], list[int]]:
    equation_count = len(series)
    power = _powers(series, degree_f, equation_count - 1)
    monomials = [
        (a, b)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    columns = []
    for a, b in monomials:
        columns.append(
            [
                power[b][order - a] if order >= a else Fraction(0)
                for order in range(equation_count)
            ]
        )
    rows = [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]
    leading = _valuation(series)
    formal_valuations = [a + b * leading for a, b in monomials]
    return rows, monomials, formal_valuations


def _differential_algebraic_matrix(
    series: RationalSeries,
    degree_x: int,
    degree_f: int,
    degree_derivative: int,
) -> tuple[list[list[Fraction]], list[tuple[int, int, int]], list[int]]:
    # f_N is needed to know [t^(N-1)] f', so N known coefficients supply
    # exactly N-1 derivative-safe equations, at orders 0,...,N-2.
    equation_count = len(series) - 1
    truncated_f = tuple(series[:equation_count])
    derivative = tuple(
        Fraction(order + 1) * series[order + 1] for order in range(equation_count)
    )
    f_powers = _powers(truncated_f, degree_f, equation_count - 1)
    derivative_powers = _powers(derivative, degree_derivative, equation_count - 1)
    monomials = [
        (a, b, c)
        for c in range(degree_derivative + 1)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    columns = []
    for a, b, c in monomials:
        product = _convolve(f_powers[b], derivative_powers[c], equation_count - 1)
        columns.append(
            [
                product[order - a] if order >= a else Fraction(0)
                for order in range(equation_count)
            ]
        )
    rows = [
        [columns[column][order] for column in range(len(columns))]
        for order in range(equation_count)
    ]
    f_valuation = _valuation(series)
    derivative_valuation = f_valuation - 1
    formal_valuations = [
        a + b * f_valuation + c * derivative_valuation for a, b, c in monomials
    ]
    return rows, monomials, formal_valuations


def _analyze_budget(
    rows: list[list[Fraction]],
    monomials: Sequence[tuple[int, ...]],
    formal_valuations: Sequence[int],
) -> dict:
    """Use U+2 prefix equations for fitting and every later equation as holdout."""

    equation_count = len(rows)
    unknown_count = len(monomials)
    training_count = unknown_count + 2
    base = {
        "unknown_count": unknown_count,
        "available_equation_count": equation_count,
        "training_equation_count": training_count,
        "holdout_equation_count": equation_count - training_count,
        "training_orders": [0, training_count - 1],
        "holdout_orders": (
            [training_count, equation_count - 1]
            if training_count < equation_count
            else []
        ),
        "monomials": [list(monomial) for monomial in monomials],
        "formal_monomial_valuations": list(formal_valuations),
    }
    try:
        training = _matrix_analysis(rows[:training_count], unknown_count)
    except ExactSolveTimeout as error:
        return {
            **base,
            "verdict": "SKIPPED_TIMEOUT",
            "timeout_phase": "training",
            "timeout_seconds": SOLVE_TIMEOUT_SECONDS,
            "detail": str(error),
        }

    base.update(
        {
            "training_rank": training["rank"],
            "training_nullity": training["nullity"],
        }
    )
    if training["nullity"] == 0:
        base.update(
            {
                "verdict": "NO_RELATION_AT_BUDGET",
                "certificate": training["full_rank_certificate"],
                "detail": (
                    "the minimally overdetermined prefix matrix has full column rank; "
                    "therefore no polynomial in this ansatz can annihilate all known coefficients"
                ),
            }
        )
        return base

    training_basis = training["kernel_basis"]
    holdout_rows = []
    residual_constraints: list[list[Fraction]] = []
    first_refuting_order = None
    try:
        for order in range(training_count, equation_count):
            residuals = [_dot(rows[order], vector) for vector in training_basis]
            residual_constraints.append(residuals)
            constraint_rank = _exact_rank(residual_constraints, len(training_basis))
            surviving_dimension = len(training_basis) - constraint_rank
            holdout_rows.append(
                {
                    "order": order,
                    "residual_on_training_kernel_basis": _fraction_strings(residuals),
                    "surviving_kernel_dimension": surviving_dimension,
                }
            )
            if surviving_dimension == 0 and first_refuting_order is None:
                first_refuting_order = order
        full = _matrix_analysis(rows, unknown_count)
    except ExactSolveTimeout as error:
        return {
            **base,
            "verdict": "SKIPPED_TIMEOUT",
            "timeout_phase": "holdout_or_full_system",
            "timeout_seconds": SOLVE_TIMEOUT_SECONDS,
            "training_kernel_basis": [list(vector) for vector in training_basis],
            "holdout_checks_completed": holdout_rows,
            "detail": str(error),
        }

    zero_columns = [
        column
        for column in range(unknown_count)
        if all(row[column] == 0 for row in rows)
    ]
    base.update(
        {
            "training_kernel_basis": [list(vector) for vector in training_basis],
            "holdout_checks": holdout_rows,
            "first_refuting_holdout_order": first_refuting_order,
            "full_rank": full["rank"],
            "full_nullity": full["nullity"],
            "zero_column_indices_at_known_order": zero_columns,
        }
    )
    if full["nullity"] == 0:
        base.update(
            {
                "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
                "certificate": full["full_rank_certificate"],
                "detail": (
                    "the overdetermined training prefix has a kernel, but later exact "
                    "coefficients eliminate every vector in that kernel"
                ),
            }
        )
        return base

    full_basis = full["kernel_basis"]
    survivor_residuals = [
        _fraction_strings(_dot(row, vector) for row in rows) for vector in full_basis
    ]
    survivors_are_unobservable = (
        full["rank"] == unknown_count - len(zero_columns)
        and full["nullity"] == len(zero_columns)
    )
    base.update(
        {
            "verdict": "CANDIDATE_SURVIVES(!)",
            "full_kernel_basis": [list(vector) for vector in full_basis],
            "full_kernel_residuals": survivor_residuals,
            "all_survivors_are_truncation_unobservable_monomials": survivors_are_unobservable,
            "audit_warning": (
                "SURVIVOR REQUIRES PARENT AUDIT; this is not an algebraic or "
                "differential-algebraic discovery"
            ),
            "detail": (
                "a nonzero kernel remains after every available coefficient; "
                "the observability fields identify whether it consists only of "
                "monomials whose first possible term lies beyond the truncation"
            ),
        }
    )
    return base


def _algebraic_frontier(series: RationalSeries) -> list[dict]:
    equation_count = len(series)
    pairs = sorted(
        (
            (degree_x, degree_f)
            for degree_x in range(equation_count)
            for degree_f in range(equation_count)
            if (degree_x + 1) * (degree_f + 1) <= equation_count - 2
        ),
        key=lambda pair: ((pair[0] + 1) * (pair[1] + 1), pair),
    )
    table = []
    for degree_x, degree_f in pairs:
        rows, monomials, valuations = _algebraic_matrix(series, degree_x, degree_f)
        table.append(
            {
                "degree_x": degree_x,
                "degree_f": degree_f,
                **_analyze_budget(rows, monomials, valuations),
            }
        )
    return table


def _differential_algebraic_frontier(series: RationalSeries) -> list[dict]:
    equation_count = len(series) - 1
    triples = sorted(
        (
            (degree_x, degree_f, degree_derivative)
            for degree_x in range(equation_count)
            for degree_f in range(equation_count)
            for degree_derivative in range(1, equation_count)
            if (degree_x + 1) * (degree_f + 1) * (degree_derivative + 1)
            <= equation_count - 2
        ),
        key=lambda triple: (
            (triple[0] + 1) * (triple[1] + 1) * (triple[2] + 1),
            triple,
        ),
    )
    table = []
    for degree_x, degree_f, degree_derivative in triples:
        rows, monomials, valuations = _differential_algebraic_matrix(
            series, degree_x, degree_f, degree_derivative
        )
        table.append(
            {
                "degree_x": degree_x,
                "degree_f": degree_f,
                "degree_f_prime": degree_derivative,
                **_analyze_budget(rows, monomials, valuations),
            }
        )
    return table


def _parity_polynomials(width: int, height: int, order: int):
    """Count open-square-lattice edge subsets by boundary parity and size."""

    lattice = square(width, height, periodic=False)
    coefficients = [[0] * (order + 1) for _ in range(1 << lattice.n_sites)]
    coefficients[0][0] = 1
    for left, right in lattice.bonds:
        toggle = (1 << left) | (1 << right)
        previous = [row[:] for row in coefficients]
        for mask, polynomial in enumerate(previous):
            target = coefficients[mask ^ toggle]
            for degree in range(order):
                target[degree + 1] += polynomial[degree]
    return lattice, coefficients


def _series_divide(
    numerator: Sequence[int | Fraction],
    denominator: Sequence[int | Fraction],
    order: int,
) -> RationalSeries:
    if not denominator[0]:
        raise ValueError("formal denominator has zero constant coefficient")
    quotient = [Fraction(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        convolution = sum(
            (
                Fraction(denominator[shift]) * quotient[degree - shift]
                for shift in range(1, degree + 1)
            ),
            Fraction(0),
        )
        quotient[degree] = (Fraction(numerator[degree]) - convolution) / denominator[0]
    return tuple(quotient)


def _finite_box_overlap_series(width: int, height: int, order: int) -> RationalSeries:
    """Return sum_{u,v in box} G_box(u,v)^2 as an exact formal series."""

    lattice, parity = _parity_polynomials(width, height, order)
    partition = parity[0]
    total = [Fraction(0) for _ in range(order + 1)]
    for source in range(lattice.n_sites):
        for target in range(lattice.n_sites):
            if source == target:
                correlation = (Fraction(1),) + (Fraction(0),) * order
            else:
                boundary = (1 << source) | (1 << target)
                correlation = _series_divide(parity[boundary], partition, order)
            square_series = _convolve(correlation, correlation, order)
            for degree, coefficient in enumerate(square_series):
                total[degree] += coefficient
    return tuple(total)


def _overlap_flm_series(order: int, bound_slack: int = 0) -> dict:
    """Independent 2D finite-lattice expansion of sum_r G(r)^2."""

    span_budget = order // 2 + bound_slack
    boxes = tuple(
        sorted(
            (
                (width, height)
                for width in range(1, span_budget + 2)
                for height in range(1, span_budget + 2)
                if (width - 1) + (height - 1) <= span_budget
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    weights: dict[tuple[int, int], RationalSeries] = {}
    bulk = [Fraction(0) for _ in range(order + 1)]
    for width, height in boxes:
        weight = list(_finite_box_overlap_series(width, height, order))
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for degree in range(order + 1):
                    weight[degree] -= placements * subweight[degree]
        exact_weight = tuple(weight)
        weights[(width, height)] = exact_weight
        for degree in range(order + 1):
            bulk[degree] += exact_weight[degree]
    return {
        "series": tuple(bulk),
        "boxes": boxes,
        "weights": weights,
        "bound_slack": bound_slack,
    }


def _interlayer_extension() -> dict:
    overlap = _overlap_flm_series(INTERLAYER_V_ORDER)
    enlarged_overlap = _overlap_flm_series(INTERLAYER_V_ORDER, bound_slack=1)
    anisotropic = anisotropic_flm_c2_series(INTERLAYER_V_ORDER)
    enlarged_anisotropic = anisotropic_flm_c2_series(
        INTERLAYER_V_ORDER, bound_slack=1
    )

    total = tuple(value / 2 for value in overlap["series"])
    residual = list(total)
    residual[0] -= Fraction(1, 2)
    residual = tuple(residual)

    c3_box_checks = []
    for shape in anisotropic.boxes:
        polynomial = anisotropic_box_even_subgraph(shape, INTERLAYER_V_ORDER, 3)
        odd_columns_zero = all(
            row[1] == 0 and row[3] == 0 for row in polynomial
        )
        c3_box_checks.append(
            {
                "shape": list(shape),
                "w1_and_w3_columns_zero": odd_columns_zero,
            }
        )

    return {
        "v_order": INTERLAYER_V_ORDER,
        "two_dimensional_overlap_flm": {
            "quantity": "sum_r G(r)^2",
            "series": _fraction_strings(overlap["series"]),
            "minimal_box_count": len(overlap["boxes"]),
            "minimal_boxes": [list(shape) for shape in overlap["boxes"]],
            "enlarged_box_count": len(enlarged_overlap["boxes"]),
            "bound_stable": overlap["series"] == enlarged_overlap["series"],
            "method": (
                "exact edge-subset parity polynomials on open 2D rectangles, exact "
                "correlation ratios, ordered source-pair sum, and rectangular "
                "finite-lattice Moebius inversion"
            ),
            "completeness_bound": (
                "after the two partition denominators are removed by finite-lattice "
                "inversion, a contributing doubled-correlation cluster is a connected "
                "even multigraph. Spanning an a-by-b rectangle costs at least "
                "2*((a-1)+(b-1)) edges, so span <= 4 is complete through v^8"
            ),
        },
        "c2_total": _fraction_strings(total),
        "c2_residual": _fraction_strings(residual),
        "anisotropic_flm_c2_total": _fraction_strings(anisotropic.total),
        "anisotropic_flm_c2_residual": _fraction_strings(anisotropic.residual),
        "new_v8_coefficient": {
            "overlap_sum": str(overlap["series"][8]),
            "c2_total": str(total[8]),
            "c2_residual": str(residual[8]),
        },
        "exact_crosscheck_through_v8": residual == anisotropic.residual,
        "anisotropic_bound_stable": anisotropic.residual == enlarged_anisotropic.residual,
        "c3": {
            "thermodynamic_coefficient": "0 identically as a formal function of v",
            "reason": (
                "for open or even-layer stacks, flipping every other layer sends "
                "w to -w; equivalently every even subgraph crosses each interlayer "
                "cut an even number of times. Fixed odd periodic layer counts have "
                "wrapping terms and are excluded before the thermodynamic limit"
            ),
            "first_two_even_v_coefficients": {"v^0": "0", "v^2": "0"},
            "exact_coefficients_through_v8": ["0"] * (INTERLAYER_V_ORDER + 1),
            "anisotropic_box_w1_w3_checks": c3_box_checks,
            "all_anisotropic_checks_pass": all(
                row["w1_and_w3_columns_zero"] for row in c3_box_checks
            ),
        },
        "scope": (
            "c2 is extended by one exact order, from v^6 to v^8. The calculation "
            "is a local expansion around w=tanh(K_z)=0 and gives no continuation "
            "to the isotropic critical point"
        ),
    }


def _verdict_counts(table: Sequence[dict]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in table:
        verdict = row["verdict"]
        counts[verdict] = counts.get(verdict, 0) + 1
    return counts


def _survivor_audit(series_name: str, relation_class: str, table: Sequence[dict]):
    rows = []
    for row in table:
        if row["verdict"] != "CANDIDATE_SURVIVES(!)":
            continue
        degrees = {
            key: row[key]
            for key in ("degree_x", "degree_f", "degree_f_prime")
            if key in row
        }
        rows.append(
            {
                "series": series_name,
                "relation_class": relation_class,
                "degrees": degrees,
                "full_nullity": row["full_nullity"],
                "zero_column_indices_at_known_order": row[
                    "zero_column_indices_at_known_order"
                ],
                "all_survivors_are_truncation_unobservable_monomials": row[
                    "all_survivors_are_truncation_unobservable_monomials"
                ],
                "warning": row["audit_warning"],
            }
        )
    return rows


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    ht, ht_metadata = _load_reduced_even_series(HT_PATH)
    lt, lt_metadata = _load_reduced_even_series(LT_PATH)

    algebraic_ht = _algebraic_frontier(ht)
    algebraic_lt = _algebraic_frontier(lt)
    da_ht = _differential_algebraic_frontier(ht)
    da_lt = _differential_algebraic_frontier(lt)
    interlayer = _interlayer_extension()

    survivor_audit = [
        *_survivor_audit("HT", "algebraic", algebraic_ht),
        *_survivor_audit("LT", "algebraic", algebraic_lt),
        *_survivor_audit("HT", "differential_algebraic_order_1", da_ht),
        *_survivor_audit("LT", "differential_algebraic_order_1", da_lt),
    ]
    all_tables = [algebraic_ht, algebraic_lt, da_ht, da_lt]
    all_rows = [row for table in all_tables for row in table]
    allowed_verdicts = {
        "NO_RELATION_AT_BUDGET",
        "CANDIDATE_REFUTED_BY_HOLDOUT",
        "CANDIDATE_SURVIVES(!)",
        "SKIPPED_TIMEOUT",
    }

    checks: list[dict] = []
    _record(
        checks,
        "input coefficient counts",
        len(ht) == 12 and len(lt) == 17,
        f"HT has {len(ht)} z=v^2 coefficients and LT has {len(lt)} u=x^2 coefficients",
    )
    _record(
        checks,
        "frontier design rule",
        all(
            row["unknown_count"] + 2 <= row["training_equation_count"]
            <= row["available_equation_count"]
            for row in all_rows
        ),
        "every fitted matrix has at least two more equations than unknowns",
    )
    _record(
        checks,
        "verdict vocabulary and completed solves",
        all(row["verdict"] in allowed_verdicts for row in all_rows)
        and not any(row["verdict"] == "SKIPPED_TIMEOUT" for row in all_rows),
        f"{len(all_rows)} exact budget rows completed; timeout limit {SOLVE_TIMEOUT_SECONDS} s per solve",
    )
    _record(
        checks,
        "stored exact linear certificates self-check",
        all(
            (
                Fraction(row["certificate"]["determinant"]) != 0
                if row["verdict"]
                in {"NO_RELATION_AT_BUDGET", "CANDIDATE_REFUTED_BY_HOLDOUT"}
                else all(
                    Fraction(value) == 0
                    for residual in row.get("full_kernel_residuals", [])
                    for value in residual
                )
            )
            for row in all_rows
        ),
        "every negative row stores a nonzero exact maximal minor; every survivor basis annihilates every available row",
    )
    _record(
        checks,
        "survivors are loudly scoped",
        all(
            row["all_survivors_are_truncation_unobservable_monomials"]
            and "PARENT AUDIT" in row["warning"]
            for row in survivor_audit
        ),
        f"{len(survivor_audit)} budget rows survive, all solely through monomials beyond the known truncation",
    )
    _record(
        checks,
        "interlayer c2 independent v8 identity",
        interlayer["exact_crosscheck_through_v8"]
        and interlayer["new_v8_coefficient"]["c2_residual"] == "778",
        "2D overlap FLM and anisotropic 3D FLM both give residual c2 coefficient [v^8]=778",
    )
    _record(
        checks,
        "interlayer box-bound stability",
        interlayer["two_dimensional_overlap_flm"]["bound_stable"]
        and interlayer["anisotropic_bound_stable"],
        "adding one span unit changes neither exact v^0,...,v^8 array",
    )
    _record(
        checks,
        "interlayer c3 odd-parity control",
        interlayer["c3"]["all_anisotropic_checks_pass"],
        "every required open anisotropic box has zero w^1 and w^3 columns through v^8",
    )

    payload = {
        "meta": {
            "provenance": SCRIPT,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "arithmetic": "exact Python Fraction and SymPy Rational only",
            "linear_solve_timeout_seconds": SOLVE_TIMEOUT_SECONDS,
            "benchmark_used_for_fit_selection_or_validation": False,
        },
        "data": {
            "input_series": {
                "HT": {
                    **ht_metadata,
                    "reduced_variable": "z=v^2",
                    "series_definition": "f_HT(z)=sum_n [v^(2n)] phi_HT * z^n",
                    "valuation": _valuation(ht),
                },
                "LT": {
                    **lt_metadata,
                    "reduced_variable": "u=x^2",
                    "series_definition": "f_LT(u)=sum_n [x^(2n)] (phi_LT-3K) * u^n",
                    "valuation": _valuation(lt),
                },
            },
            "design_rule": {
                "rule": (
                    "for U monomial coefficients use the first U+2 exact coefficient "
                    "equations as training; test every later known coefficient as a "
                    "strict holdout. Sweep every rectangular degree tuple with U <= E-2"
                ),
                "reason": (
                    "U+2 makes the training system overdetermined by two equations. "
                    "Budgets with fewer than U equations have a guaranteed kernel and "
                    "are not evidence of structure"
                ),
                "holdout_selection": (
                    "prefix length depends only on U; no coefficient value, holdout "
                    "residual, critical benchmark, or surviving candidate affects it"
                ),
            },
            "algebraicity": {
                "ansatz": "P(t,f)=sum_{a=0}^dx sum_{b=0}^df p_ab t^a f^b",
                "HT": {
                    "equation_count": len(ht),
                    "frontier_row_count": len(algebraic_ht),
                    "verdict_counts": _verdict_counts(algebraic_ht),
                    "frontier": algebraic_ht,
                },
                "LT": {
                    "equation_count": len(lt),
                    "frontier_row_count": len(algebraic_lt),
                    "verdict_counts": _verdict_counts(algebraic_lt),
                    "frontier": algebraic_lt,
                },
                "scope": (
                    "a NO_RELATION verdict is an exact statement only for the displayed "
                    "bidegree and known prefix, not a proof of transcendence"
                ),
            },
            "differential_algebraicity_order_1": {
                "ansatz": (
                    "Q(t,f,f')=sum q_abc t^a f^b (f')^c with rectangular "
                    "multidegrees; c>=1 is included in every searched space"
                ),
                "derivative_safe_equations": (
                    "N coefficients of f provide N-1 exact equations because the next "
                    "coefficient would be needed for the next derivative coefficient"
                ),
                "HT": {
                    "equation_count": len(ht) - 1,
                    "frontier_row_count": len(da_ht),
                    "verdict_counts": _verdict_counts(da_ht),
                    "frontier": da_ht,
                },
                "LT": {
                    "equation_count": len(lt) - 1,
                    "frontier_row_count": len(da_lt),
                    "verdict_counts": _verdict_counts(da_lt),
                    "frontier": da_lt,
                },
                "scope": (
                    "this probes first derivative order in finite rectangular polynomial "
                    "budgets only; it neither decides differential algebraicity nor tests "
                    "higher derivatives"
                ),
            },
            "survivor_parent_audit": {
                "warning": (
                    "CANDIDATE_SURVIVES(!) IS NOT A DISCOVERY. Every survivor below "
                    "must be audited against new coefficients before interpretation"
                ),
                "row_count": len(survivor_audit),
                "rows": survivor_audit,
                "observed_explanation": (
                    "all current survivors are exactly the span of zero matrix columns: "
                    "high powers of the positive-valuation LT series (or derivative) "
                    "cannot begin before the available truncation ends"
                ),
            },
            "interlayer_extension": interlayer,
        },
        "checks": checks,
    }
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more exact checks failed")
    print("PASS e43_series_structure")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e43_series_structure: {error}")
        raise
