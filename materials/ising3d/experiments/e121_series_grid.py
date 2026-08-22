"""Expanded exact finite-prefix structure grid for the 3D Ising free-energy series.

The four ansatz matrix conventions are imported from e43/e56 without changing
monomial order: algebraic and first-order differential-algebraic matrices from
e43, Euler-ODE and Mahler matrices from e56.  This experiment only adds the
strict current-prefix classification and certificates.

Run from the repository root:

    PYTHONPATH=src .venv/bin/python experiments/e121_series_grid.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import signal
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Callable, Iterable, Sequence

SCRIPT = "experiments/e121_series_grid.py"
ROOT = Path(__file__).resolve().parents[1]
HT_PATH = ROOT / "results" / "series" / "ht_v28.json"
HT_PREDECESSOR_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
LT_PATH = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
E43_PATH = ROOT / "experiments" / "e43_series_structure.py"
E56_PATH = ROOT / "experiments" / "e56_series_frontier2.py"
E43_RESULT_PATH = ROOT / "results" / "series" / "structure_certificates.json"
RESULT_PATH = ROOT / "results" / "series" / "grid.json"
PROOF_PATH = ROOT / "proofs" / "series_grid.md"

CELL_WALL_SECONDS = 60
ALLOWED_VERDICTS = (
    "NO_RELATION_AT_BUDGET",
    "CANDIDATE_REFUTED_BY_HOLDOUT",
    "CANDIDATE_SURVIVES",
)

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


class CellWallExceeded(TimeoutError):
    """Raised when one requested exact grid cell exceeds its declared wall."""


@contextmanager
def _cell_deadline(seconds: int = CELL_WALL_SECONDS):
    """Bound a single exact cell without silently changing its mathematical result."""

    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_timer = signal.getitimer(signal.ITIMER_REAL)

    def _raise_timeout(_signum, _frame):
        raise CellWallExceeded(f"exact grid cell exceeded its {seconds}-second wall")

    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_timer[0] > 0:
            signal.setitimer(signal.ITIMER_REAL, *previous_timer)


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ImportError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_strings(values: Iterable[Fraction | int]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _load_ht_series() -> tuple[tuple[Fraction, ...], dict]:
    payload = json.loads(HT_PATH.read_text(encoding="utf-8"))
    data = payload["data"]
    series_data = data["series"]
    coefficients = tuple(Fraction(value) for value in series_data["coefficients"])
    achieved_order = int(series_data["achieved_order"])
    if series_data["variable"] != "v" or achieved_order != 28:
        raise AssertionError("ht_v28 must provide the canonical v^0,...,v^28 series")
    if len(coefficients) != achieved_order + 1:
        raise AssertionError("ht_v28 coefficient length does not match achieved order")
    if any(coefficients[order] for order in range(1, len(coefficients), 2)):
        raise AssertionError("canonical HT series must be even through v^28")

    predecessor = json.loads(HT_PREDECESSOR_PATH.read_text(encoding="utf-8"))
    predecessor_coefficients = tuple(
        Fraction(value) for value in predecessor["data"]["coefficients"]
    )
    if coefficients[: len(predecessor_coefficients)] != predecessor_coefficients:
        raise AssertionError("ht_v28 does not extend the frozen v^22 prefix exactly")

    return coefficients, {
        "claim_tag": "[COMPUTATION]",
        "artifact": str(HT_PATH.relative_to(ROOT)),
        "sha256": _sha256(HT_PATH),
        "variable": "v",
        "achieved_order": achieved_order,
        "coefficient_count": len(coefficients),
        "coefficients": _fraction_strings(coefficients),
        "series_definition": "F_HT(v)=phi_HT(v)-log(2)",
        "predecessor_prefix": {
            "artifact": str(HT_PREDECESSOR_PATH.relative_to(ROOT)),
            "sha256": _sha256(HT_PREDECESSOR_PATH),
            "last_order": len(predecessor_coefficients) - 1,
            "matches_exactly": True,
        },
    }


def _load_lt_series() -> tuple[tuple[Fraction, ...], dict]:
    payload = json.loads(LT_PATH.read_text(encoding="utf-8"))
    data = payload["data"]
    coefficients = tuple(Fraction(value) for value in data["coefficients"])
    achieved_order = int(data["achieved_order"])
    if data["variable"] != "x" or achieved_order != 32:
        raise AssertionError("LT input must provide x^0,...,x^32")
    if len(coefficients) != achieved_order + 1:
        raise AssertionError("LT coefficient length does not match achieved order")
    if any(coefficients[order] for order in range(1, len(coefficients), 2)):
        raise AssertionError("canonical LT series must be even through x^32")
    return coefficients, {
        "claim_tag": "[COMPUTATION]",
        "artifact": str(LT_PATH.relative_to(ROOT)),
        "sha256": _sha256(LT_PATH),
        "variable": "x",
        "achieved_order": achieved_order,
        "coefficient_count": len(coefficients),
        "coefficients": _fraction_strings(coefficients),
        "series_definition": "F_LT(x)=phi_LT(x)-3K",
    }


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append(
        {
            "claim_tag": "[COMPUTATION]",
            "name": name,
            "passed": bool(passed),
            "detail": detail,
        }
    )
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def _algebraic_cells(equation_count: int) -> list[tuple[int, int]]:
    return sorted(
        (
            (degree_x, degree_f)
            for degree_x in range(equation_count)
            for degree_f in range(equation_count)
            if (degree_x + 1) * (degree_f + 1) + 2 <= equation_count
        ),
        key=lambda degrees: ((degrees[0] + 1) * (degrees[1] + 1), degrees),
    )


def _first_order_da_cells(equation_count: int) -> list[tuple[int, int, int]]:
    return sorted(
        (
            (degree_x, degree_f, degree_f_prime)
            for degree_x in range(equation_count)
            for degree_f in range(equation_count)
            for degree_f_prime in range(1, equation_count)
            if (degree_x + 1) * (degree_f + 1) * (degree_f_prime + 1) + 2
            <= equation_count
        ),
        key=lambda degrees: (
            (degrees[0] + 1) * (degrees[1] + 1) * (degrees[2] + 1),
            degrees,
        ),
    )


def _euler_cells(equation_count: int) -> list[tuple[int, int]]:
    return [
        (order, degree)
        for order in range(3)
        for degree in range(equation_count)
        if (order + 1) * (degree + 1) + 2 <= equation_count
    ]


def _mahler_cells(equation_count: int) -> list[tuple[int, int]]:
    cells = []
    for dependent_total_degree in (1, 2):
        for degree_t in range(equation_count):
            unknown_count = (
                (degree_t + 1)
                * (dependent_total_degree + 1)
                * (dependent_total_degree + 2)
                // 2
            )
            if unknown_count + 2 <= equation_count:
                cells.append((degree_t, dependent_total_degree))
    return sorted(
        cells,
        key=lambda degrees: (
            (degrees[0] + 1) * (degrees[1] + 1) * (degrees[1] + 2) // 2,
            degrees,
        ),
    )


def _budget_skip_summary(family: str, equation_count: int, admissible_count: int) -> dict:
    """Record compact monotone frontiers for unsupported cells, not fake verdicts.

    e43 bounds its rectangular degree scans by range(E), and e56 scans r=0,1,2
    or D=1,2 with degree_t below E.  The records below preserve those finite
    candidate boxes while avoiding tens of thousands of redundant no-solve rows.
    """

    records: list[dict] = []
    if family == "algebraic":
        candidate_count = equation_count * equation_count
        for degree_x in range(equation_count):
            first_degree_f = (equation_count - 2) // (degree_x + 1)
            if first_degree_f < equation_count:
                unknown_count = (degree_x + 1) * (first_degree_f + 1)
                records.append(
                    {
                        "record_type": "SKIPPED_INSUFFICIENT_EQUATIONS",
                        "parameters": {
                            "degree_x": degree_x,
                            "degree_f": first_degree_f,
                        },
                        "unknown_count": unknown_count,
                        "required_training_equation_count": unknown_count + 2,
                        "available_equation_count": equation_count,
                        "excluded_monotone_region": (
                            f"degree_x={degree_x}, degree_f>={first_degree_f} "
                            f"within 0..{equation_count - 1}"
                        ),
                    }
                )
    elif family == "first_order_differential_algebraic":
        candidate_count = equation_count * equation_count * (equation_count - 1)
        for degree_x in range(equation_count):
            for degree_f in range(equation_count):
                first_degree_f_prime = max(
                    1, (equation_count - 2) // ((degree_x + 1) * (degree_f + 1))
                )
                if first_degree_f_prime < equation_count:
                    unknown_count = (
                        (degree_x + 1)
                        * (degree_f + 1)
                        * (first_degree_f_prime + 1)
                    )
                    records.append(
                        {
                            "record_type": "SKIPPED_INSUFFICIENT_EQUATIONS",
                            "parameters": {
                                "degree_x": degree_x,
                                "degree_f": degree_f,
                                "degree_f_prime": first_degree_f_prime,
                            },
                            "unknown_count": unknown_count,
                            "required_training_equation_count": unknown_count + 2,
                            "available_equation_count": equation_count,
                            "excluded_monotone_region": (
                                f"degree_x={degree_x}, degree_f={degree_f}, "
                                f"degree_f_prime>={first_degree_f_prime} "
                                f"within 1..{equation_count - 1}"
                            ),
                        }
                    )
    elif family == "linear_euler_ode":
        candidate_count = 3 * equation_count
        for order in range(3):
            first_degree = (equation_count - 2) // (order + 1)
            if first_degree < equation_count:
                unknown_count = (order + 1) * (first_degree + 1)
                records.append(
                    {
                        "record_type": "SKIPPED_INSUFFICIENT_EQUATIONS",
                        "parameters": {"order": order, "polynomial_degree": first_degree},
                        "unknown_count": unknown_count,
                        "required_training_equation_count": unknown_count + 2,
                        "available_equation_count": equation_count,
                        "excluded_monotone_region": (
                            f"order={order}, polynomial_degree>={first_degree} "
                            f"within 0..{equation_count - 1}"
                        ),
                    }
                )
    elif family == "mahler":
        candidate_count = 2 * equation_count
        for dependent_total_degree in (1, 2):
            multiplier = (dependent_total_degree + 1) * (dependent_total_degree + 2) // 2
            first_degree_t = (equation_count - 2) // multiplier
            if first_degree_t < equation_count:
                unknown_count = (first_degree_t + 1) * multiplier
                records.append(
                    {
                        "record_type": "SKIPPED_INSUFFICIENT_EQUATIONS",
                        "parameters": {
                            "degree_t": first_degree_t,
                            "dependent_total_degree": dependent_total_degree,
                        },
                        "unknown_count": unknown_count,
                        "required_training_equation_count": unknown_count + 2,
                        "available_equation_count": equation_count,
                        "excluded_monotone_region": (
                            f"dependent_total_degree={dependent_total_degree}, "
                            f"degree_t>={first_degree_t} within 0..{equation_count - 1}"
                        ),
                    }
                )
    else:
        raise ValueError(f"unsupported family {family}")

    skipped_count = candidate_count - admissible_count
    if skipped_count < 0:
        raise AssertionError("budget summary has negative skipped count")
    if not all(
        record["required_training_equation_count"]
        > record["available_equation_count"]
        for record in records
    ):
        raise AssertionError("a skip frontier does not actually exceed the available equations")
    return {
        "claim_tag": "[COMPUTATION]",
        "candidate_parameter_box_count": candidate_count,
        "admissible_cell_count": admissible_count,
        "skipped_cell_count": skipped_count,
        "frontier_records": records,
        "interpretation": (
            "Each record is the first dimensionally unsupported point of a monotone "
            "degree region. Unsupported points are recorded separately and never "
            "receive a mathematical verdict."
        ),
    }


def _parameters(family: str, degrees: tuple[int, ...]) -> dict:
    if family == "algebraic":
        return {"degree_x": degrees[0], "degree_f": degrees[1]}
    if family == "first_order_differential_algebraic":
        return {
            "degree_x": degrees[0],
            "degree_f": degrees[1],
            "degree_f_prime": degrees[2],
        }
    if family == "linear_euler_ode":
        return {"order": degrees[0], "polynomial_degree": degrees[1]}
    if family == "mahler":
        return {"degree_t": degrees[0], "dependent_total_degree": degrees[1]}
    raise ValueError(f"unsupported family {family}")


def _family_definition(family: str) -> tuple[str, str]:
    definitions = {
        "algebraic": (
            "P(t,F)=sum_(a,b) p_ab t^a F(t)^b; rectangular bidegrees",
            "experiments/e43_series_structure.py::_algebraic_matrix",
        ),
        "first_order_differential_algebraic": (
            "Q(t,F,F')=sum_(a,b,c) q_abc t^a F(t)^b (F'(t))^c; "
            "rectangular tridegrees with degree_f_prime>=1",
            "experiments/e43_series_structure.py::_differential_algebraic_matrix",
        ),
        "linear_euler_ode": (
            "sum_(j=0)^r Q_j(t) theta^j F=0, theta=t*d/dt, r<=2, deg Q_j<=d",
            "experiments/e56_series_frontier2.py::_linear_ode_matrix",
        ),
        "mahler": (
            "M(t,F(t),F(t^2))=sum_(a,b,c) m_abc t^a F(t)^b F(t^2)^c, b+c<=D, D in {1,2}",
            "experiments/e56_series_frontier2.py::_mahler_matrix",
        ),
    }
    return definitions[family]


def _survivor_support_analysis(
    *,
    rows: Sequence[Sequence[Fraction]],
    monomials: Sequence[tuple[int, ...]],
    formal_valuations: Sequence[int],
    full_basis: Sequence[tuple[int, ...]],
    e56,
    certificates: dict[str, dict],
    context: dict,
) -> tuple[list[int], dict]:
    """Classify a full kernel as literal-zero-column or sparse-support induced.

    e56 first handles literal zero columns before its progression certificate.
    That branch matters for the old e43 reduced LT rows: their zero columns
    need not share one global congruence support progression.
    """

    columns = tuple(range(len(monomials)))
    zero_columns = [
        column for column in columns if all(row[column] == 0 for row in rows)
    ]
    visible_columns = tuple(column for column in columns if column not in zero_columns)
    full_nullity = len(full_basis)
    if zero_columns and full_nullity == len(zero_columns):
        visible_rank, _ = e56._fraction_free_rank_and_rows(rows, visible_columns)
        if visible_rank == len(visible_columns):
            certificate = e56._minor_certificate(
                rows, visible_columns, "observable_columns"
            )
            certificate["zero_column_indices"] = zero_columns
            certificate["unobservable_monomials"] = [
                {
                    "column_index": column,
                    "monomial": list(monomials[column]),
                    "lowest_possible_relation_order": formal_valuations[column],
                    "last_available_relation_order": len(rows) - 1,
                }
                for column in zero_columns
            ]
            pointer = e56._store_certificate(certificates, certificate, context)
            first_testable_order = min(
                formal_valuations[column] for column in zero_columns
            )
            if first_testable_order < len(rows):
                raise AssertionError("zero-column survivor has an observable formal valuation")
            return zero_columns, {
                "classification": "TRUNCATION_UNOBSERVABLE",
                "support_certificate_pointer": pointer,
                "support_components": [
                    {
                        "type": "ZERO_COLUMN",
                        "column_index": column,
                        "monomial": list(monomials[column]),
                        "first_unfixed_supported_order": formal_valuations[column],
                        "last_available_order": len(rows) - 1,
                    }
                    for column in zero_columns
                ],
                "first_testable_order": first_testable_order,
                "detail": (
                    "The exact full kernel is the coordinate span of literal zero "
                    "columns whose formal first terms lie beyond the available prefix."
                ),
            }

    support = e56._unobservable_subspace_certificate(
        rows,
        monomials,
        formal_valuations,
        full_basis,
        certificates,
        context,
    )
    if support is None:
        return zero_columns, {
            "classification": "OBSERVABLE_FINITE_PREFIX_KERNEL",
            "support_certificate_pointer": None,
            "support_components": [],
            "first_testable_order": len(rows),
            "detail": (
                "The exact full kernel is not accounted for by a literal-zero-column "
                "or inherited support-progression certificate."
            ),
        }
    support_certificate, support_pointer = support
    unresolved_components = [
        component
        for component in support_certificate["support_components"]
        if component["known_nullity"] > 0
    ]
    if not unresolved_components:
        raise AssertionError("support certificate has a kernel but no unresolved component")
    first_testable_order = min(
        component["first_unfixed_supported_order"]
        for component in unresolved_components
    )
    return zero_columns, {
        "classification": "TRUNCATION_UNOBSERVABLE",
        "support_certificate_pointer": support_pointer,
        "support_components": support_certificate["support_components"],
        "first_testable_order": first_testable_order,
        "detail": (
            "Every full-kernel direction is accounted for by exact support "
            "components whose next supported equation lies beyond the available prefix."
        ),
    }


def _analyze_cell(
    *,
    series_name: str,
    family: str,
    parameters: dict,
    rows: Sequence[Sequence[Fraction]],
    monomials: Sequence[tuple[int, ...]],
    formal_valuations: Sequence[int],
    e56,
    certificates: dict[str, dict],
) -> dict:
    """Classify one strictly predeclared dimensionally admissible cell."""

    started = time.perf_counter()
    unknown_count = len(monomials)
    equation_count = len(rows)
    training_count = unknown_count + 2
    if training_count > equation_count:
        raise AssertionError("unsupported cell reached exact analysis")
    columns = tuple(range(unknown_count))
    context = {
        "claim_tag": "[COMPUTATION]",
        "series": series_name,
        "ansatz_family": family,
        "parameters": parameters,
    }
    base = {
        "claim_tag": "[COMPUTATION]",
        "ansatz_family": family,
        "parameters": parameters,
        "unknown_count": unknown_count,
        "available_equation_count": equation_count,
        "training_equation_count": training_count,
        "holdout_equation_count": equation_count - training_count,
        "training_orders": [0, training_count - 1],
        "holdout_orders": [training_count, equation_count - 1]
        if training_count < equation_count
        else [],
        "monomials": [list(monomial) for monomial in monomials],
        "formal_monomial_valuations": list(formal_valuations),
    }

    try:
        with _cell_deadline():
            training_rows = rows[:training_count]
            training_rank, _ = e56._fraction_free_rank_and_rows(training_rows, columns)
            base["training_rank"] = training_rank
            base["training_nullity"] = unknown_count - training_rank
            if training_rank == unknown_count:
                certificate = e56._minor_certificate(
                    training_rows, columns, "training_prefix"
                )
                pointer = e56._store_certificate(certificates, certificate, context)
                return {
                    **base,
                    "verdict": "NO_RELATION_AT_BUDGET",
                    "certificate_pointer": pointer,
                    "detail": (
                        "The U+2 training prefix has full exact Q-column rank; "
                        "the stored nonzero maximal minor certifies this finite ansatz exclusion."
                    ),
                    "elapsed_seconds": round(time.perf_counter() - started, 9),
                }

            training_basis = e56._nullspace(training_rows, unknown_count)
            if len(training_basis) != unknown_count - training_rank:
                raise AssertionError("Fraction nullspace and exact Bareiss rank disagree")
            residual_constraints: list[list[Fraction]] = []
            holdout_checks = []
            first_refuting_order = None
            first_refuting_residual = None
            for order in range(training_count, equation_count):
                residuals = [e56._dot(rows[order], vector) for vector in training_basis]
                residual_constraints.append(residuals)
                residual_rank, _ = e56._fraction_free_rank_and_rows(
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
                    first_refuting_residual = _fraction_strings(residuals)

            full_rank, _ = e56._fraction_free_rank_and_rows(rows, columns)
            full_nullity = unknown_count - full_rank
            candidate_fields = {
                "training_kernel_basis": [list(vector) for vector in training_basis],
                "holdout_checks": holdout_checks,
                "first_refuting_holdout_order": first_refuting_order,
                "full_rank": full_rank,
                "full_nullity": full_nullity,
            }
            if full_rank == unknown_count:
                if first_refuting_order is None or first_refuting_residual is None:
                    raise AssertionError("full rank after a kernel lacks a first holdout refutation")
                certificate = e56._minor_certificate(
                    rows, columns, "all_available_equations"
                )
                pointer = e56._store_certificate(certificates, certificate, context)
                return {
                    **base,
                    **candidate_fields,
                    "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
                    "first_refuting_residual_on_training_kernel_basis": first_refuting_residual,
                    "certificate_pointer": pointer,
                    "detail": (
                        "An exact kernel of the prescribed training prefix is eliminated "
                        "by withheld coefficient equations."
                    ),
                    "elapsed_seconds": round(time.perf_counter() - started, 9),
                }

            full_basis = e56._nullspace(rows, unknown_count)
            if len(full_basis) != full_nullity:
                raise AssertionError("full Fraction kernel and exact Bareiss rank disagree")
            zero_columns, support_analysis = _survivor_support_analysis(
                rows=rows,
                monomials=monomials,
                formal_valuations=formal_valuations,
                full_basis=full_basis,
                e56=e56,
                certificates=certificates,
                context=context,
            )
            first_testable_order = support_analysis["first_testable_order"]
            return {
                **base,
                **candidate_fields,
                "verdict": "CANDIDATE_SURVIVES",
                "full_kernel_basis": [list(vector) for vector in full_basis],
                "full_kernel_residuals": [
                    _fraction_strings(e56._dot(row, vector) for row in rows)
                    for vector in full_basis
                ],
                "zero_column_indices_at_known_order": zero_columns,
                "support_analysis": support_analysis,
                "first_testable_order": first_testable_order,
                "detail": (
                    "This is not a discovery: the displayed support analysis classifies "
                    "whether the surviving finite kernel is truncation-unobservable."
                ),
                "elapsed_seconds": round(time.perf_counter() - started, 9),
            }
    except CellWallExceeded as error:
        raise RuntimeError(
            f"[UNRESOLVED] resource wall in {series_name} {family} {parameters}: {error}"
        ) from error


def _build_family(
    *,
    series_name: str,
    series: Sequence[Fraction],
    family: str,
    e43,
    e56,
    certificates: dict[str, dict],
) -> dict:
    definition, matrix_builder = _family_definition(family)
    if family == "algebraic":
        equation_count = len(series)
        degree_cells = _algebraic_cells(equation_count)
        build: Callable[[tuple[int, ...]], tuple] = (
            lambda degrees: e43._algebraic_matrix(series, *degrees)
        )
    elif family == "first_order_differential_algebraic":
        equation_count = len(series) - 1
        degree_cells = _first_order_da_cells(equation_count)
        build = lambda degrees: e43._differential_algebraic_matrix(series, *degrees)
    elif family == "linear_euler_ode":
        equation_count = len(series)
        degree_cells = _euler_cells(equation_count)
        build = lambda degrees: e56._linear_ode_matrix(series, *degrees)
    elif family == "mahler":
        equation_count = len(series)
        degree_cells = _mahler_cells(equation_count)
        build = lambda degrees: e56._mahler_matrix(series, *degrees)
    else:
        raise ValueError(f"unsupported family {family}")

    cells = []
    for degrees in degree_cells:
        parameters = _parameters(family, degrees)
        rows, monomials, valuations = build(degrees)
        if len(rows) != equation_count:
            raise AssertionError("shared matrix builder returned an unexpected equation count")
        cells.append(
            _analyze_cell(
                series_name=series_name,
                family=family,
                parameters=parameters,
                rows=rows,
                monomials=monomials,
                formal_valuations=valuations,
                e56=e56,
                certificates=certificates,
            )
        )

    verdict_counts: dict[str, int] = {verdict: 0 for verdict in ALLOWED_VERDICTS}
    for cell in cells:
        if cell["verdict"] not in ALLOWED_VERDICTS:
            raise AssertionError("grid emitted a forbidden verdict")
        verdict_counts[cell["verdict"]] += 1
    if sum(verdict_counts.values()) != len(cells):
        raise AssertionError("verdict counts do not cover every admissible cell")
    skip_summary = _budget_skip_summary(family, equation_count, len(cells))
    return {
        "claim_tag": "[COMPUTATION]",
        "ansatz": definition,
        "matrix_builder": matrix_builder,
        "available_equation_count": equation_count,
        "cells": cells,
        "verdict_counts": verdict_counts,
        "budget_skip_summary": skip_summary,
    }


def _legacy_lt_reaudit(
    lt: Sequence[Fraction], e43, e56, certificates: dict[str, dict]
) -> dict:
    """Replay the 18 e43 LT zero-column survivors in their original u=x^2 convention."""

    stored = json.loads(E43_RESULT_PATH.read_text(encoding="utf-8"))
    audit = stored["data"]["survivor_parent_audit"]
    stored_rows = audit["rows"]
    if len(stored_rows) != 18 or any(row["series"] != "LT" for row in stored_rows):
        raise AssertionError("e43 stored survivor audit no longer has the expected 18 LT rows")
    old_input = stored["data"]["input_series"]["LT"]
    reduced = tuple(lt[order] for order in range(0, len(lt), 2))
    current_hash = _sha256(LT_PATH)
    if int(old_input["reduced_coefficient_count"]) != len(reduced):
        raise AssertionError("e43 LT reduction length changed unexpectedly")

    replayed = []
    changed_rows = []
    for index, stored_row in enumerate(stored_rows):
        degrees = stored_row["degrees"]
        relation_class = stored_row["relation_class"]
        if relation_class == "algebraic":
            rows, monomials, valuations = e43._algebraic_matrix(
                reduced, degrees["degree_x"], degrees["degree_f"]
            )
        elif relation_class == "differential_algebraic_order_1":
            rows, monomials, valuations = e43._differential_algebraic_matrix(
                reduced,
                degrees["degree_x"],
                degrees["degree_f"],
                degrees["degree_f_prime"],
            )
        else:
            raise AssertionError(f"unexpected e43 relation class {relation_class}")
        columns = tuple(range(len(monomials)))
        full_rank, _ = e56._fraction_free_rank_and_rows(rows, columns)
        full_basis = e56._nullspace(rows, len(monomials))
        context = {
            "claim_tag": "[COMPUTATION]",
            "series": "LT",
            "ansatz_family": f"legacy_e43_{relation_class}",
            "parameters": degrees,
            "legacy_row_index": index,
        }
        zero_columns, support_analysis = _survivor_support_analysis(
            rows=rows,
            monomials=monomials,
            formal_valuations=valuations,
            full_basis=full_basis,
            e56=e56,
            certificates=certificates,
            context=context,
        )
        support_pointer = support_analysis["support_certificate_pointer"]
        first_testable_order = support_analysis["first_testable_order"]
        support_classification = support_analysis["classification"]
        unchanged = (
            len(full_basis) == stored_row["full_nullity"]
            and zero_columns == stored_row["zero_column_indices_at_known_order"]
            and support_classification == "TRUNCATION_UNOBSERVABLE"
        )
        replayed_row = {
            "claim_tag": "[COMPUTATION]",
            "legacy_row_index": index,
            "relation_class": relation_class,
            "degrees": degrees,
            "full_rank": full_rank,
            "full_nullity": len(full_basis),
            "zero_column_indices_at_known_order": zero_columns,
            "support_classification": support_classification,
            "support_certificate_pointer": support_pointer,
            "first_testable_order_in_u": first_testable_order,
            "unchanged": unchanged,
        }
        replayed.append(replayed_row)
        if not unchanged:
            changed_rows.append(replayed_row)

    prior_hash = old_input["sha256"]
    same_input = prior_hash == current_hash
    return {
        "claim_tag": "[COMPUTATION]",
        "source_artifact": str(E43_RESULT_PATH.relative_to(ROOT)),
        "source_sha256": _sha256(E43_RESULT_PATH),
        "stored_row_count": len(stored_rows),
        "replayed_row_count": len(replayed),
        "prior_lt_input_sha256": prior_hash,
        "current_lt_input_sha256": current_hash,
        "same_lt_input_artifact": same_input,
        "newly_available_equation_orders": [] if same_input else None,
        "rows": replayed,
        "changed_rows": changed_rows,
        "scope": (
            "The original e43 variable u=x^2 and its original algebraic/first-order "
            "differential-algebraic builders are replayed exactly; this is separate "
            "from the raw-x grid."
        ),
    }


def _verdict_counts(frontiers: dict) -> dict[str, int]:
    counts = {verdict: 0 for verdict in ALLOWED_VERDICTS}
    for series in frontiers.values():
        for family in series["families"].values():
            for verdict, count in family["verdict_counts"].items():
                counts[verdict] += count
    return counts


def _survivor_summary(frontiers: dict) -> tuple[list[dict], list[dict]]:
    all_survivors = []
    unexplained = []
    for series_name, series in frontiers.items():
        for family_name, family in series["families"].items():
            for cell in family["cells"]:
                if cell["verdict"] != "CANDIDATE_SURVIVES":
                    continue
                support = cell["support_analysis"]
                entry = {
                    "claim_tag": "[COMPUTATION]",
                    "series": series_name,
                    "ansatz_family": family_name,
                    "parameters": cell["parameters"],
                    "full_nullity": cell["full_nullity"],
                    "support_classification": support["classification"],
                    "first_testable_order": cell["first_testable_order"],
                }
                all_survivors.append(entry)
                if support["classification"] != "TRUNCATION_UNOBSERVABLE":
                    unexplained.append(entry)
    return all_survivors, unexplained


def _proof_text(payload: dict) -> str:
    data = payload["data"]
    summary = data["summary"]
    ht = data["frontiers"]["HT"]
    lt = data["frontiers"]["LT"]
    return rf"""# Expanded exact series-structure grid

## Fixed ansatz spaces and protocol

[LEMMA] For every displayed ansatz cell with $U$ homogeneous unknown
coefficients, this experiment uses precisely coefficient equations
$0,\ldots,U+1$ as the training prefix and reserves every later safe equation
as a holdout.  Algebraic, Euler, and Mahler matrices use all available stored
orders.  A first-order differential-algebraic row at order $n$ is safe only
through $n=N-2$ when $N$ coefficients of $F$ are known, because the next
coefficient is needed for $[t^n]F'$.

[THEOREM] Let $A$ be one of these finite exact training matrices over
$\mathbb{{Q}}$ for a fixed displayed ansatz space.  If a stored $U\times U$
minor of $A$ has nonzero determinant, then no nonzero coefficient vector in
that ansatz space annihilates the prescribed training prefix.  Indeed, the
minor gives $\operatorname{{rank}}_\mathbb{{Q}} A=U$, so its nullspace is zero.
This theorem is only about the stated finite matrix and ansatz space.

[COMPUTATION] The canonical HT input is `ht_v28.json`, with 29 exact raw
coefficients through $v^{{28}}$, and the LT input has 33 exact raw coefficients
through $x^{{32}}$.  The HT prefix through $v^{{22}}$ agrees exactly with the
frozen predecessor.  The matrix conventions are imported directly from e43
(algebraic and first-order differential-algebraic) and e56 (Euler and Mahler);
no numerical benchmark enters selection or validation.

[COMPUTATION] The evaluated raw-variable grid has
{summary['evaluated_cell_count']} admissible cells: HT has
{ht['admissible_cell_count']} and LT has {lt['admissible_cell_count']}.  Its
verdict counts are `{json.dumps(summary['verdict_counts'], sort_keys=True)}`.
The separately recorded bounded-degree envelopes contain
{summary['skipped_by_budget_count']} dimensionally unsupported cells; their
monotone boundary records identify why no training solve was attempted.  Such
skips are not verdicts.

[COMPUTATION] Every `NO_RELATION_AT_BUDGET` row stores a nonzero exact maximal
minor computed after rowwise denominator clearing by fraction-free integer
Bareiss elimination.  Every `CANDIDATE_REFUTED_BY_HOLDOUT` row stores its first
withheld refuting order and exact residual vector on the primitive training
kernel basis.  Exact modular rank was not used to decide any row.  In general
a modular rank can only be a lower bound on rational rank; it would require an
exact confirmation before supporting a no-relation claim.

[COMPUTATION] There are {len(summary['survivors'])} rows labelled
`CANDIDATE_SURVIVES`, but all are certified
`TRUNCATION_UNOBSERVABLE` support kernels, with their complete kernel bases,
residue components, and first testable supported orders in `grid.json`.  Thus
there are {summary['unexplained_survivor_count']} unexplained finite-prefix
survivors.  None of the truncation-unobservable rows is called a discovery.

## Re-audit of the e43 LT survivors

[COMPUTATION] The 18 stored e43 LT truncation-unobservable survivors were
replayed in their original $u=x^2$ convention.  The current LT artifact hash
is identical to the one recorded by e43, so no additional equation order has
become available.  All {data['legacy_lt_survivor_reaudit']['replayed_row_count']}
full nullities and zero-column sets agree; the changed-row list is
`{data['legacy_lt_survivor_reaudit']['changed_rows']}`.

## Scope

[UNRESOLVED] These are finite exact frontier facts only.  Finite-budget
emptiness is not a non-existence proof for the true infinite free energy, and
a full-rank cell does not prove transcendence, non-algebraicity,
non-differential-algebraicity, or non-D-finiteness.  The calculation neither
claims an exact solution of the simple-cubic 3D Ising model nor excludes
relations outside the displayed ansatz spaces or degree budgets.
"""


def main() -> None:
    e43 = _load_module("series_grid_e43", E43_PATH)
    e56 = _load_module("series_grid_e56", E56_PATH)
    ht, ht_metadata = _load_ht_series()
    lt, lt_metadata = _load_lt_series()
    certificates: dict[str, dict] = {}
    started = time.perf_counter()

    frontiers = {}
    for series_name, series, metadata in (("HT", ht, ht_metadata), ("LT", lt, lt_metadata)):
        families = {}
        for family in (
            "algebraic",
            "first_order_differential_algebraic",
            "linear_euler_ode",
            "mahler",
        ):
            families[family] = _build_family(
                series_name=series_name,
                series=series,
                family=family,
                e43=e43,
                e56=e56,
                certificates=certificates,
            )
        frontiers[series_name] = {
            "claim_tag": "[COMPUTATION]",
            "input": metadata,
            "admissible_cell_count": sum(len(family["cells"]) for family in families.values()),
            "families": families,
        }

    legacy_reaudit = _legacy_lt_reaudit(lt, e43, e56, certificates)
    counts = _verdict_counts(frontiers)
    survivors, unexplained = _survivor_summary(frontiers)
    evaluated_count = sum(counts.values())
    skipped_count = sum(
        family["budget_skip_summary"]["skipped_cell_count"]
        for series in frontiers.values()
        for family in series["families"].values()
    )
    candidate_envelope_count = evaluated_count + skipped_count
    checks: list[dict] = []
    all_cells = [
        cell
        for series in frontiers.values()
        for family in series["families"].values()
        for cell in family["cells"]
    ]
    _record(
        checks,
        "strict training and holdout rule",
        all(
            cell["training_equation_count"] == cell["unknown_count"] + 2
            and cell["training_orders"] == [0, cell["unknown_count"] + 1]
            and cell["training_equation_count"] <= cell["available_equation_count"]
            for cell in all_cells
        ),
        f"all {len(all_cells)} admissible cells train on exactly 0..U+1",
    )
    _record(
        checks,
        "exact three-label verdict vocabulary",
        all(cell["verdict"] in ALLOWED_VERDICTS for cell in all_cells),
        "unsupported budget regions are separate records, never extra verdict labels",
    )
    _record(
        checks,
        "nonzero minor certificates",
        all(
            Fraction(certificates[cell["certificate_pointer"].rsplit("/", 1)[-1]]["determinant"])
            != 0
            for cell in all_cells
            if cell["verdict"] in {
                "NO_RELATION_AT_BUDGET",
                "CANDIDATE_REFUTED_BY_HOLDOUT",
            }
        ),
        "every finite negative verdict has an exact nonzero maximal minor",
    )
    _record(
        checks,
        "holdout refutation residuals",
        all(
            cell["first_refuting_holdout_order"] is not None
            and any(
                Fraction(value) != 0
                for value in cell["first_refuting_residual_on_training_kernel_basis"]
            )
            for cell in all_cells
            if cell["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
        ),
        "every refuted training kernel has a recorded nonzero first holdout residual",
    )
    _record(
        checks,
        "survivor support completeness",
        not unexplained
        and all(
            survivor["first_testable_order"] is not None for survivor in survivors
        ),
        f"{len(survivors)} surviving kernels are all truncation-unobservable",
    )
    _record(
        checks,
        "legacy LT survivor re-audit",
        legacy_reaudit["stored_row_count"] == 18
        and legacy_reaudit["replayed_row_count"] == 18
        and legacy_reaudit["same_lt_input_artifact"]
        and not legacy_reaudit["changed_rows"],
        "all stored e43 LT zero-column survivors reproduce exactly at x^32",
    )

    payload = {
        "meta": {
            "provenance": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "arithmetic": (
                "exact Python Fraction; rowwise denominator clearing and fraction-free "
                "integer Bareiss elimination; Fraction RREF only for exact kernel bases"
            ),
            "cell_wall_seconds": CELL_WALL_SECONDS,
            "benchmark_used_for_fit_selection_or_validation": False,
            "matrix_builder_sources": {
                "e43": {
                    "path": str(E43_PATH.relative_to(ROOT)),
                    "sha256": _sha256(E43_PATH),
                },
                "e56": {
                    "path": str(E56_PATH.relative_to(ROOT)),
                    "sha256": _sha256(E56_PATH),
                },
            },
        },
        "data": {
            "claim_tag": "[COMPUTATION]",
            "input_series": {"HT": ht_metadata, "LT": lt_metadata},
            "protocol": {
                "claim_tag": "[LEMMA]",
                "training_rule": (
                    "For U unknowns, use exactly orders 0..U+1 as the training prefix; "
                    "every later safe equation is a strict holdout."
                ),
                "first_order_derivative_safety": (
                    "With N coefficients of F, Q(t,F,F') has derivative-safe equation "
                    "orders 0..N-2 only."
                ),
                "verdict_vocabulary": list(ALLOWED_VERDICTS),
                "certificate_theorem": {
                    "claim_tag": "[THEOREM]",
                    "statement": (
                        "A nonzero U-by-U minor over Q of a displayed finite training "
                        "matrix proves that its U-column nullspace is zero."
                    ),
                },
                "modular_rank_scope": (
                    "No modular rank decides a verdict here. A modular rank is only a "
                    "lower bound on rational rank and would need exact confirmation "
                    "before a no-relation classification."
                ),
            },
            "frontiers": frontiers,
            "certificates": certificates,
            "legacy_lt_survivor_reaudit": legacy_reaudit,
            "summary": {
                "claim_tag": "[COMPUTATION]",
                "verdict_vocabulary": sorted(ALLOWED_VERDICTS),
                "candidate_envelope_cell_count": candidate_envelope_count,
                "evaluated_cell_count": evaluated_count,
                "skipped_by_budget_count": skipped_count,
                "verdict_counts": counts,
                "survivors": survivors,
                "truncation_unobservable_survivor_count": sum(
                    survivor["support_classification"] == "TRUNCATION_UNOBSERVABLE"
                    for survivor in survivors
                ),
                "unexplained_survivor_count": len(unexplained),
                "unexplained_survivors": unexplained,
                "total_elapsed_seconds": round(time.perf_counter() - started, 9),
            },
            "scope": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "This finite grid neither proves non-D-finiteness nor proves the "
                    "absence of relations for the true infinite free energy."
                ),
            },
        },
        "checks": checks,
    }
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    PROOF_PATH.write_text(_proof_text(payload), encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    print(f"WROTE {PROOF_PATH.relative_to(ROOT)}")
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more exact grid checks failed")
    print("PASS e121_series_grid")


if __name__ == "__main__":
    main()
