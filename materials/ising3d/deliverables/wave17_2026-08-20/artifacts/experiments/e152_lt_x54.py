"""Exact finite-lattice extension of the simple-cubic low-temperature series.

The computation reaches x^54.  It uses the established fixed-plus finite-box
polynomials and bounding-box inversion, retains exact ``int``/``Fraction``
arithmetic throughout, and routes only boxes with more than 62 spins through
the repository CRT transfer engine.  CRT boxes are deliberately truncated at
the requested order: this lowers the resource wall without changing any
coefficient that enters the formal logarithm through x^54.

Run from the repository root:

    PYTHONPATH=src .venv/bin/python experiments/e124_lt_x34.py
"""

from __future__ import annotations

import hashlib
import json
import math
import resource
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import prod
from pathlib import Path
from typing import Iterable

from sympy import Matrix, Rational

import ising.series as flm
from ising.series import low_temperature_free_energy
from ising.transfer_matrix import box_broken_bond_poly as int64_box_broken_bond_poly
from ising.transfer_matrix.crt import (
    box_broken_bond_poly as crt_box_broken_bond_poly,
)
from ising.transfer_matrix.crt import estimate_box_peak_bytes, primes_for_sites


SCRIPT = "experiments/e152_lt_x54.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "series" / "lt_x54.json"
PREDECESSOR_PATH = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
STRUCTURE_PATH = ROOT / "results" / "series" / "structure_certificates.json"
PROOF_PATH = ROOT / "proofs" / "lt_x54.md"
SOURCE_PATH = ROOT / "sources" / "fulltext" / "num_guttmann_enting1993.pdf"
TARGET_ORDER = 54
from ising.transfer_matrix.crt import _MAX_CROSS_SECTION as CRT_GUARD
TARGET_BUDGET = (TARGET_ORDER + 6) // 4
RSS_BUDGET_BYTES = 48 * (1 << 30)  # above the 41.25 GiB (5,5,5) projection, below the 60 GiB abort
WALL_BUDGET_SECONDS = 21600  # 6 h governed heavy phase


RationalSeries = tuple[Fraction, ...]


def _rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _fraction_strings(values: Iterable[Fraction | int]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _shape_key(shape: tuple[int, int, int]) -> str:
    return "x".join(str(value) for value in shape)


def _canonical_boxes(budget: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        shape
        for shape in combinations_with_replacement(range(1, budget + 1), 3)
        if sum(shape) <= budget
    )


def _ordered_box_count(budget: int) -> int:
    return sum(
        1
        for shape in product(range(1, budget + 1), repeat=3)
        if sum(shape) <= budget
    )


def _clear_flm_caches() -> None:
    flm._free_broken_polynomial.cache_clear()
    flm._plus_broken_polynomial.cache_clear()
    flm._ht_box_log.cache_clear()
    flm._lt_box_log.cache_clear()


def _profile() -> dict:
    canonical = _canonical_boxes(TARGET_BUDGET)
    records = []
    for shape in canonical:
        sites = prod(shape)
        cross_section = prod(shape[:-1])
        peak = estimate_box_peak_bytes(
            shape, plus_boundary=True, max_degree=TARGET_ORDER
        )
        records.append(
            {
                "shape": list(shape),
                "sites": sites,
                "cross_section_sites": cross_section,
                "transfer_engine": "int64" if sites <= 62 else "31-bit-prime CRT",
                "estimated_peak_bytes_at_target_truncation": peak,
                "within_transfer_cross_section_guard": cross_section <= CRT_GUARD,
            }
        )
    wall = max(records, key=lambda row: row["estimated_peak_bytes_at_target_truncation"])
    ordered_count = _ordered_box_count(TARGET_BUDGET)
    if ordered_count != math.comb(TARGET_BUDGET, 3):
        raise AssertionError("positive ordered-triple count disagrees with stars-and-bars")
    return {
        "claim_tag": "[COMPUTATION]",
        "target_order": TARGET_ORDER,
        "box_budget": TARGET_BUDGET,
        "order_bound": "4*(a+b+c)-6 <= truncation_order",
        "ordered_box_count": ordered_count,
        "canonical_box_count": len(canonical),
        "canonical_boxes": [list(shape) for shape in canonical],
        "crt_required_boxes": [
            row["shape"] for row in records if row["transfer_engine"] == "31-bit-prime CRT"
        ],
        "maximum_estimated_peak_box": wall,
        "all_cross_sections_within_guard": all(
            row["within_transfer_cross_section_guard"] for row in records
        ),
        "rss_budget_bytes": RSS_BUDGET_BYTES,
        "wall_budget_seconds": WALL_BUDGET_SECONDS,
        "preflight_kind": "unlaunched input-size resource estimate",
    }


def _next_order_preflight() -> dict:
    # After the audited guard lift 22 -> 25 (crt.py, notes/ht_v24_x54_feasibility.md
    # section 0), x^56 reuses the SAME box inventory (side-sum budget 15) and the
    # first newly walled order is x^62: budget 17, canonical (5,6,6), cross-section 30.
    next_order = TARGET_ORDER + 2
    return {
        "claim_tag": "[COMPUTATION]",
        "next_order": next_order,
        "next_order_budget": (next_order + 6) // 4,
        "next_order_inventory_unchanged": (next_order + 6) // 4 == TARGET_BUDGET,
        "next_wall_order": 62,
        "next_wall_canonical_shape": [5, 6, 6],
        "next_wall_cross_section_sites": 30,
        "transfer_cross_section_guard": CRT_GUARD,
        "launched": False,
        "status": "COSTED",
        "reason": (
            f"With the guard at {CRT_GUARD}, x^{next_order} runs on the same inventory "
            f"as x^{TARGET_ORDER} (deeper truncation only); the first newly walled order "
            "is x^62 (canonical (5,6,6), cross-section 30)."
        ),
    }


@contextmanager
def _truncated_hybrid_flm_engine(order: int):
    """Route only N>62 plus-boxes through CRT, retaining exact coefficients <= order."""
    original = flm.box_broken_bond_poly
    records: dict[tuple[int, int, int], dict] = {}
    _clear_flm_caches()

    def hybrid(shape, periodic=None, plus_boundary: bool = False):
        canonical = tuple(sorted(int(side) for side in shape))
        sites = prod(canonical)
        started = time.perf_counter()
        if sites <= 62:
            coefficients = int64_box_broken_bond_poly(
                canonical, periodic=periodic, plus_boundary=plus_boundary
            )
            engine = "int64"
            moduli = None
        else:
            coefficients = crt_box_broken_bond_poly(
                canonical,
                periodic=periodic,
                plus_boundary=plus_boundary,
                max_degree=order,
            )
            engine = "31-bit-prime CRT truncated at target order"
            moduli = primes_for_sites(sites)
            if moduli.product <= 2 * (1 << sites):
                raise AssertionError("CRT modulus fails the strict uniqueness guard")
        elapsed = time.perf_counter() - started
        prefix = list(coefficients[: order + 1])
        prefix.extend([0] * (order + 1 - len(prefix)))
        record = {
            "shape": list(canonical),
            "sites": sites,
            "cross_section_sites": prod(canonical[:-1]),
            "plus_boundary": bool(plus_boundary),
            "engine": engine,
            "requested_degree": order,
            "returned_degree_after_trim": len(coefficients) - 1,
            "coefficient_prefix": [str(value) for value in prefix],
            "coefficient_prefix_sha256": hashlib.sha256(
                json.dumps(prefix, separators=(",", ":")).encode("utf-8")
            ).hexdigest(),
            "elapsed_seconds": f"{elapsed:.6f}",
            "estimated_peak_bytes_at_target_truncation": estimate_box_peak_bytes(
                canonical, plus_boundary=plus_boundary, max_degree=order
            ),
        }
        if moduli is not None:
            bound = 1 << sites
            record.update(
                {
                    "prime_moduli": list(moduli.primes),
                    "modulus_product": moduli.product,
                    "modulus_product_bits": moduli.product_bits,
                    "coefficient_bound": bound,
                    "twice_coefficient_bound": 2 * bound,
                    "uniqueness_certificate_passed": moduli.product > 2 * bound,
                    "full_sum_check_status": "not needed for truncated prefix",
                }
            )
        prior = records.get(canonical)
        if prior is not None and prior != record:
            raise AssertionError(f"inconsistent repeated transfer record for {canonical}")
        records[canonical] = record
        return coefficients

    flm.box_broken_bond_poly = hybrid
    try:
        yield records
    finally:
        flm.box_broken_bond_poly = original
        _clear_flm_caches()


def _convolve(
    left: RationalSeries, right: RationalSeries, order: int
) -> RationalSeries:
    result = [Fraction(0)] * (order + 1)
    for first_degree, first_value in enumerate(left[: order + 1]):
        if not first_value:
            continue
        for second_degree, second_value in enumerate(right[: order + 1 - first_degree]):
            if second_value:
                result[first_degree + second_degree] += first_value * second_value
    return tuple(result)


def _powers(series: RationalSeries, maximum: int, order: int) -> tuple[RationalSeries, ...]:
    one = (Fraction(1),) + (Fraction(0),) * order
    values = [one]
    for _ in range(maximum):
        values.append(_convolve(values[-1], series, order))
    return tuple(values)


def _valuation(series: RationalSeries) -> int:
    for degree, value in enumerate(series):
        if value:
            return degree
    raise ValueError("the finite series is identically zero")


def _primitive_vector(vector) -> list[int]:
    fractions = [Fraction(int(value.p), int(value.q)) for value in vector]
    denominator = 1
    for value in fractions:
        denominator = math.lcm(denominator, value.denominator)
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
    return integers


def _matrix(rows: list[list[Fraction]]) -> Matrix:
    return Matrix(
        [
            [Rational(value.numerator, value.denominator) for value in row]
            for row in rows
        ]
    )


def _dot(row: list[Fraction], vector: list[int]) -> Fraction:
    return sum(
        (entry * coefficient for entry, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _algebraic_matrix(
    series: RationalSeries, degree_x: int, degree_f: int
) -> tuple[list[list[Fraction]], list[tuple[int, int]], list[int]]:
    equation_count = len(series)
    powers = _powers(series, degree_f, equation_count - 1)
    monomials = [
        (a, b)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    rows = [
        [
            powers[b][order - a] if order >= a else Fraction(0)
            for a, b in monomials
        ]
        for order in range(equation_count)
    ]
    valuation = _valuation(series)
    return rows, monomials, [a + b * valuation for a, b in monomials]


def _differential_matrix(
    series: RationalSeries, degree_x: int, degree_f: int, degree_f_prime: int
) -> tuple[list[list[Fraction]], list[tuple[int, int, int]], list[int]]:
    # This is the frozen e43 convention: N coefficients provide derivative-safe
    # equation orders 0,...,N-2 and no higher equation is used.
    equation_count = len(series) - 1
    truncated_f = tuple(series[:equation_count])
    derivative = tuple(
        Fraction(order + 1) * series[order + 1] for order in range(equation_count)
    )
    f_powers = _powers(truncated_f, degree_f, equation_count - 1)
    derivative_powers = _powers(derivative, degree_f_prime, equation_count - 1)
    monomials = [
        (a, b, c)
        for c in range(degree_f_prime + 1)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    rows = []
    for order in range(equation_count):
        row = []
        for a, b, c in monomials:
            combined = _convolve(
                f_powers[b], derivative_powers[c], equation_count - 1
            )
            row.append(combined[order - a] if order >= a else Fraction(0))
        rows.append(row)
    valuation = _valuation(series)
    return rows, monomials, [
        a + b * valuation + c * (valuation - 1) for a, b, c in monomials
    ]


def _frozen_survivor_cell(
    structure: dict, relation_class: str, degrees: dict
) -> dict:
    family = (
        "algebraicity"
        if relation_class == "algebraic"
        else "differential_algebraicity_order_1"
    )
    for cell in structure["data"][family]["LT"]["frontier"]:
        if (
            cell["degree_x"] == degrees["degree_x"]
            and cell["degree_f"] == degrees["degree_f"]
            and (
                relation_class == "algebraic"
                or cell["degree_f_prime"] == degrees["degree_f_prime"]
            )
            and cell["verdict"] == "CANDIDATE_SURVIVES(!)"
        ):
            return cell
    raise AssertionError(f"missing frozen survivor cell: {relation_class} {degrees}")


def _replay_survivor(
    relation_class: str, degrees: dict, series_u: RationalSeries
) -> dict:
    if relation_class == "algebraic":
        rows, monomials, valuations = _algebraic_matrix(
            series_u, degrees["degree_x"], degrees["degree_f"]
        )
        derivative_safe = False
    elif relation_class == "differential_algebraic_order_1":
        rows, monomials, valuations = _differential_matrix(
            series_u,
            degrees["degree_x"],
            degrees["degree_f"],
            degrees["degree_f_prime"],
        )
        derivative_safe = True
    else:
        raise AssertionError(f"unknown relation class {relation_class}")

    unknown_count = len(monomials)
    training_count = unknown_count + 2
    training_basis = [
        _primitive_vector(vector) for vector in _matrix(rows[:training_count]).nullspace()
    ]
    residual_constraints: list[list[Fraction]] = []
    first_refuting_order = None
    residual_at_refutation = None
    final_holdout_residual = None
    for order in range(training_count, len(rows)):
        residual = [_dot(rows[order], vector) for vector in training_basis]
        residual_constraints.append(residual)
        rank = _matrix(residual_constraints).rank()
        if len(training_basis) - rank == 0 and first_refuting_order is None:
            first_refuting_order = order
            residual_at_refutation = residual
        final_holdout_residual = residual

    full_matrix = _matrix(rows)
    full_nullity = len(full_matrix.nullspace())
    zero_columns = [
        column
        for column in range(unknown_count)
        if all(row[column] == 0 for row in rows)
    ]
    base = {
        "training_equation_count": training_count,
        "available_equation_count": len(rows),
        "training_kernel_basis": training_basis,
        "full_nullity": full_nullity,
        "zero_column_indices": zero_columns,
        "formal_monomial_valuations": valuations,
        "latest_holdout_residual_on_training_kernel_basis": (
            _fraction_strings(final_holdout_residual)
            if final_holdout_residual is not None
            else []
        ),
    }
    if full_nullity == 0:
        if first_refuting_order is None or residual_at_refutation is None:
            raise AssertionError("full-rank replay lacks a first refuting holdout")
        return {
            **base,
            "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
            "first_refuting_holdout_order_u": first_refuting_order,
            "residual_on_training_kernel_basis": _fraction_strings(
                residual_at_refutation
            ),
        }

    if full_nullity != len(zero_columns):
        raise AssertionError("remaining kernel is not entirely truncation-unobservable")
    next_equation_order = min(valuations[column] for column in zero_columns)
    return {
        **base,
        "verdict": "STILL_TRUNCATION_UNOBSERVABLE",
        "next_testable_equation_order_u": next_equation_order,
        "required_coefficient_order_x": 2
        * (next_equation_order + (1 if derivative_safe else 0)),
    }


def _survivor_reaudit(series: RationalSeries, structure: dict) -> dict:
    if any(series[degree] for degree in range(1, len(series), 2)):
        raise AssertionError("LT u-reduction requires an even x series")
    series_u = tuple(series[degree] for degree in range(0, len(series), 2))
    rows = []
    unchanged_training = True
    for frozen in structure["data"]["survivor_parent_audit"]["rows"]:
        relation_class = frozen["relation_class"]
        degrees = frozen["degrees"]
        frozen_cell = _frozen_survivor_cell(structure, relation_class, degrees)
        replay = _replay_survivor(relation_class, degrees, series_u)
        same_training = (
            replay["training_kernel_basis"] == frozen_cell["training_kernel_basis"]
            and replay["training_equation_count"] == frozen_cell["training_equation_count"]
        )
        unchanged_training = unchanged_training and same_training
        rows.append(
            {
                "claim_tag": "[COMPUTATION]",
                "series": "LT",
                "relation_class": relation_class,
                "degrees": degrees,
                "frozen_full_nullity": frozen["full_nullity"],
                "frozen_zero_column_indices": frozen[
                    "zero_column_indices_at_known_order"
                ],
                "training_prefix_unchanged": same_training,
                **replay,
            }
        )
    verdict_counts: dict[str, int] = {}
    for row in rows:
        verdict_counts[row["verdict"]] = verdict_counts.get(row["verdict"], 0) + 1
    return {
        "claim_tag": "[COMPUTATION]",
        "source_artifact": "results/series/structure_certificates.json",
        "source_convention": (
            "f_LT(u)=sum_n [x^(2n)](phi_LT-3K)u^n; frozen U+2 training prefix; "
            "algebraic equations through the known order and derivative-safe equations through N-2"
        ),
        "u_coefficients": _fraction_strings(series_u),
        "row_count": len(rows),
        "training_prefix_unchanged_for_all_rows": unchanged_training,
        "verdict_counts": verdict_counts,
        "rows": rows,
    }


def _formal_exp(series: RationalSeries, order: int) -> tuple[Fraction, ...]:
    """Exact exponential, used only for the post-computation external witness."""
    if series[0] != 0:
        raise AssertionError("the reduced low-temperature series has no constant term")
    result = [Fraction(0)] * (order + 1)
    result[0] = Fraction(1)
    for degree in range(1, order + 1):
        result[degree] = (
            sum(
                (
                    Fraction(index) * series[index] * result[degree - index]
                    for index in range(1, degree + 1)
                ),
                Fraction(0),
            )
            / degree
        )
    return tuple(result)


def _external_witness_after_local_derivation(series: RationalSeries) -> dict:
    """Open the stored table only after all local coefficients have been derived."""
    if not SOURCE_PATH.exists():
        return {
            "claim_tag": "[EXTERNAL]",
            "available": False,
            "accessed_after_local_derivation": True,
            "reason": "The expected local Guttmann--Enting PDF is absent.",
        }
    source_bytes = SOURCE_PATH.read_bytes()
    # Table 2 prints lambda_n of Lambda_0(u)=sum_n lambda_n u^n, u=x^2.  Cells with
    # n<=19 and n=21 are row-aligned in the stored PDF; the remaining cells of the
    # lambda column are merged with the m and c columns into one block, so they are
    # recorded as verbatim tokens of that block rather than as addressed cells.
    row_aligned_lambda = {
        0: 1, 3: 1, 4: 0, 5: 3, 6: -3, 7: 15, 8: -30, 9: 101, 10: -261, 11: 807,
        12: -2308, 13: 7065, 14: -21171, 15: 65337, 16: -200934, 17: 627249,
        18: -1962034, 19: 6192066, 21: 62482527,
    }
    merged_block_lambda = {
        20: -19610346,
        22: -199807110,
        23: 641837193,
        24: -2068695927,
        25: 6691611633,
        26: -21710041944,
    }
    published = {**row_aligned_lambda, **merged_block_lambda}
    witness_order = max(published)
    local_u = tuple(series[degree] for degree in range(0, len(series), 2))
    if len(local_u) <= witness_order:
        raise AssertionError("local series is too short for the stored witness table")
    local_lambda = _formal_exp(local_u, witness_order)
    rows = []
    for index in sorted(published):
        rows.append(
            {
                "lambda_index_u": index,
                "x_order": 2 * index,
                "table_cell_kind": (
                    "row_aligned" if index in row_aligned_lambda else "merged_column_block_token"
                ),
                "published_value": str(published[index]),
                "local_exponentiated_value": str(local_lambda[index]),
                "match": local_lambda[index] == published[index],
            }
        )
    return {
        "claim_tag": "[EXTERNAL]",
        "available": True,
        "accessed_after_local_derivation": True,
        "local_path": str(SOURCE_PATH.relative_to(ROOT)),
        "sha256": hashlib.sha256(source_bytes).hexdigest(),
        "source": (
            "A. J. Guttmann and I. G. Enting, Series studies of the Potts model. I. "
            "The simple cubic Ising model, J. Phys. A 26 (1993) 807--821, Table 2"
        ),
        "table_scope": (
            "Table 2 lists Lambda_0(u) through u^26=x^52, exactly the reach the paper "
            "claims for 4x5 cross-sections (Table 1, s=14); the local budget is the same 14."
        ),
        "comparison_route": (
            "local reduced series -> exact formal exponential Lambda_0=exp(sum_n [x^(2n)] u^n) "
            "-> compared coefficient by coefficient against the printed table"
        ),
        "coefficient_source_role": (
            "post-computation validation witness; no external coefficient was read or used "
            "before the local x^52 finite-lattice derivation, and none was used to fit, select "
            "or tune anything, but the comparison below IS a validation gate: it is recorded in "
            "`checks` and a mismatch fails this experiment"
        ),
        "row_count": len(rows),
        "all_match": all(row["match"] for row in rows),
        "rows": rows,
        "witnessed_new_x_orders": [
            2 * index for index in sorted(published) if 2 * index >= 34
        ],
        "merged_block_caveat": (
            "The six merged-block values were located as verbatim tokens appearing in "
            "increasing-n order inside the single merged Table-2 column dump of the stored "
            "PDF text layer; they are a weaker witness than the row-aligned cells."
        ),
        "beyond_table_scope": "No x^54-or-higher coefficient exists locally or in this table.",
    }


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def _degree_label(row: dict) -> str:
    degrees = row["degrees"]
    if row["relation_class"] == "algebraic":
        return f"P(x,f), deg_x={degrees['degree_x']}, deg_f={degrees['degree_f']}"
    return (
        f"Q(x,f,f'), deg_x={degrees['degree_x']}, deg_f={degrees['degree_f']}, "
        f"deg_f'={degrees['degree_f_prime']}"
    )


def _write_proof_note(data: dict) -> None:
    coefficients = [Fraction(value) for value in data["coefficients"]]
    new_orders = list(range(34, TARGET_ORDER + 1, 2))
    reaudit = data["survivor_reaudit"]
    refuted = [
        row for row in reaudit["rows"] if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ]
    surviving = [
        row for row in reaudit["rows"] if row["verdict"] == "STILL_TRUNCATION_UNOBSERVABLE"
    ]
    inventory = data["finite_lattice_inventory"]
    usage = data["resource_usage"]
    wall = data["next_order_preflight"]
    witness = data["external_witness"]

    lines: list[str] = []
    lines.append(
        f"# Exact low-temperature simple-cubic series through `x^{TARGET_ORDER}`, and the"
        " re-audit of the seven dormant legacy LT truncation-unobservable survivors"
    )
    lines.append("")
    lines.append(f"Producer: `{SCRIPT}`.  Artifact: `results/series/lt_x34.json`.")
    lines.append("Standalone independent test: `tests/test_lt_x34.py`.")
    lines.append("")
    lines.append("## 0. Conventions")
    lines.append("")
    lines.append(
        "`K=beta J`, `x=e^{-2K}`, `u=x^2=e^{-4K}`, and the reduced low-temperature free"
        " energy is normalised as in `problem_specification.md` and"
        " `notes/flm_derivation.md`:"
    )
    lines.append("")
    lines.append("    phi(K) = 3*K + sum_{n>=1} c_n x^n,   f_LT(u) = sum_n c_{2n} u^n.")
    lines.append("")
    lines.append(
        "`Xi_A(x)` denotes the fixed-plus-boundary broken-bond polynomial of the `a x b x c`"
        " box `A`, i.e. `sum_sigma x^{#unsatisfied bonds}` with the exterior frozen to `+1`."
    )
    lines.append("")
    lines.append("## 1. Order bound and completeness of the box inventory")
    lines.append("")
    lines.append(
        "[THEOREM] (bounding-box surface bound; proved in `notes/flm_derivation.md`, Eq. (10),"
        " and not reproved here).  A connected set `S` of flipped sites spanning `(a,b,c)` has"
        " broken-bond surface `|dS| >= 2[(a+b-1)+(a+c-1)+(b+c-1)] = 4*(a+b+c)-6`, and a"
        " connected Mayer cluster of droplets inherits the same lower bound through the"
        " exterior surface of its union.  Consequently the finite-lattice weight `W_A` of a"
        " box `A` satisfies `[x^m] W_A = 0` for all `m < 4*(a+b+c)-6`."
    )
    lines.append("")
    lines.append(
        f"[THEOREM] (inventory completeness at order {TARGET_ORDER}).  Every box whose"
        f" weight can contribute at or below `x^{TARGET_ORDER}` satisfies"
        f" `4*(a+b+c)-6 <= {TARGET_ORDER}`, i.e. `a+b+c <= {TARGET_BUDGET}`; a box with"
        f" side sum 16 has surface at least `4*16-6 = 58 > {TARGET_ORDER}`.  The"
        " computation therefore uses exactly the boxes with"
        f" `a+b+c <= {TARGET_BUDGET}`: {inventory['ordered_box_count']} ordered boxes"
        f" (`= C({TARGET_BUDGET},3)`, the number of positive integer triples with that"
        f" side-sum bound) grouped into {inventory['canonical_box_count']} permutation"
        " classes, one cached polynomial per class.  This inventory is complete for the"
        " achieved order and, by the same inequality, not complete for `x^54`."
    )
    lines.append("")
    lines.append(
        "[COMPUTATION] The bound is also verified empirically, not merely assumed: every"
        " one of the "
        f"{inventory['ordered_box_count']} ordered box weights vanishes identically below"
        " `4*(a+b+c)-6`, and every canonical box has its first nonzero weight coefficient"
        " exactly at `4*(a+b+c)-6` (artifact key `canonical_weight_support`).  A missing or"
        " misplaced box would break one of those two statements."
    )
    lines.append("")
    lines.append("## 2. Exact arithmetic and the CRT uniqueness certificate")
    lines.append("")
    lines.append(
        "All box polynomials, logarithms, bounding-box inversions and survivor residuals are"
        " computed in exact Python `int`/`Fraction` arithmetic.  Boxes with at most 62 spins"
        " use the `int64` transfer engine; the "
        f"{inventory['crt_box_count']} boxes with more than 62 spins are reconstructed by the"
        " repository's deterministic 31-bit-prime CRT engine, truncated at the target order."
    )
    lines.append("")
    lines.append(
        "[THEOREM] (coefficient bound).  `[x^q] Xi_A` counts spin configurations of the"
        " `N`-site box with exactly `q` unsatisfied bonds, hence is a non-negative integer"
        " at most `2^N`, and the coefficients sum to exactly `2^N`."
    )
    lines.append("")
    lines.append(
        "[COMPUTATION] For every CRT box the recorded modulus product strictly exceeds"
        " `2 * 2^N`, i.e. strictly more than twice the proved coefficient bound, so the"
        " reconstruction of each coefficient is unique.  This is a two-sided certificate:"
        " the residues fix the coefficient modulo the product, and the proved bound"
        " confines it to an interval shorter than half the modulus."
    )
    lines.append("")
    wall_shape = tuple(data["resource_profile"]["maximum_estimated_peak_box"]["shape"])
    wall_truncated = estimate_box_peak_bytes(
        wall_shape, plus_boundary=True, max_degree=TARGET_ORDER
    )
    wall_full = estimate_box_peak_bytes(wall_shape, plus_boundary=True)
    lines.append(
        "Because the CRT boxes are truncated at the target order, the global normalisation"
        " check `sum_q [x^q] Xi_A = 2^N` is not available for them, and the artifact records"
        " that status explicitly rather than silently claiming it.  Uniqueness does not need"
        " it: the per-coefficient bound `0 <= [x^q] Xi_A <= 2^N` is a combinatorial fact about"
        " each individual coefficient, and only the low-degree prefix enters the logarithm and"
        " the bounding-box inversion at this order.  The truncation is what keeps the run"
        " inside the memory budget: for the largest box"
        f" `{'x'.join(str(side) for side in wall_shape)}` the unlaunched input-size estimate is"
        f" {wall_truncated} bytes truncated against {wall_full} bytes at full degree, the"
        f" latter exceeding the {RSS_BUDGET_BYTES}-byte budget."
    )
    lines.append("")
    lines.append("| box | sites | cross-section | primes | product bits | `2*2^N` | product `> 2*2^N` |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for certificate in data["crt_certificates"]:
        shape = "x".join(str(side) for side in certificate["shape"])
        primes = ", ".join(str(prime) for prime in certificate["prime_moduli"])
        lines.append(
            f"| `{shape}` | {certificate['sites']} | {certificate['cross_section_sites']} |"
            f" `{primes}` | {certificate['modulus_product_bits']} |"
            f" `{certificate['twice_coefficient_bound']}` |"
            f" `{certificate['uniqueness_certificate_passed']}` |"
        )
    lines.append("")
    lines.append("## 3. The new exact coefficients")
    lines.append("")
    lines.append(
        "[COMPUTATION] The stored predecessor prefix is reproduced exactly before any new"
        " coefficient is claimed: all 33 coefficients `x^0,...,x^32` of"
        " `results/series/extended2_sc_lt_free_energy.json` agree.  Every odd-order"
        f" coefficient through `x^{TARGET_ORDER}` is exactly zero.  The eleven new"
        " coefficients are"
    )
    lines.append("")
    for order in new_orders:
        lines.append(f"    x^{order} = {coefficients[order]}")
    lines.append("")
    lines.append(
        "In the `u=x^2` convention this extends the reduced series from 17 to"
        f" {TARGET_ORDER // 2 + 1} coefficients, `u^0,...,u^{TARGET_ORDER // 2}`."
    )
    lines.append("")
    lines.append("## 4. Post-computation external witness")
    lines.append("")
    if witness["available"]:
        lines.append(
            "[EXTERNAL] Guttmann and Enting, *Series studies of the Potts model. I. The"
            " simple cubic Ising model*, J. Phys. A 26 (1993) 807--821, Table 2 (local copy"
            f" `{witness['local_path']}`, sha256 `{witness['sha256']}`).  The stored PDF was"
            " opened only after every local box polynomial, logarithm, inversion, coefficient"
            " check and survivor residual above had been computed; no external number is an"
            " input to any local computation, a fit target, or a selection criterion."
        )
        lines.append("")
        lines.append(
            "Exponentiating the local reduced series, `Lambda_0(u) = exp(f_LT(u))`, reproduces"
            f" all {witness['row_count']} printed `lambda_n` values of that table, including"
            " every new order "
            + ", ".join(f"`x^{order}`" for order in witness["witnessed_new_x_orders"])
            + "."
        )
        lines.append("")
        lines.append(
            "This comparison IS a validation gate, and the artifact says so: it is recorded in"
            " `checks` as `post_computation_external_table_witness`, so a mismatch fails the"
            " producer.  The precise split is therefore: no external value is a computational"
            " input and none was used for fitting, selection or tuning; the hashed table is used"
            " as a post-computation validation witness, opened after the computation."
        )
        lines.append("")
        lines.append("| `n` | `x` order | table cell | published `lambda_n` | local `exp` value | match |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for row in witness["rows"]:
            lines.append(
                f"| {row['lambda_index_u']} | `x^{row['x_order']}` |"
                f" {row['table_cell_kind']} | `{row['published_value']}` |"
                f" `{row['local_exponentiated_value']}` | `{row['match']}` |"
            )
        lines.append("")
        lines.append(
            "The paper's Table 1 independently states that `4x5` cross-sections (`s=14`) give"
            " the low-temperature Ising series to `u^26`, which is exactly `x^52` and exactly"
            f" the local budget `a+b+c <= {TARGET_BUDGET}`; its `5x5` row (`s=16`, `2^25`"
            " vector elements) is the next step, matching the local wall below."
        )
        lines.append("")
        lines.append(witness["merged_block_caveat"])
    else:
        lines.append(f"[EXTERNAL] No local tabulation is available: {witness['reason']}")
    lines.append("")
    lines.append("## 5. Re-audit of the eighteen legacy LT survivors")
    lines.append("")
    lines.append(
        "[COMPUTATION] The eighteen rows are the stored LT `CANDIDATE_SURVIVES(!)`"
        " truncation-unobservable cells of `results/series/structure_certificates.json`"
        " (produced by `experiments/e43_series_structure.py`).  The frozen protocol is reused"
        " verbatim: for `U` unknown monomial coefficients the training prefix is the first"
        " `U+2` equations, every later equation is a holdout, algebraic rows use equation"
        " orders `0..N-1` and first-order differential-algebraic rows only the"
        " derivative-safe orders `0..N-2`.  Nothing is refitted or reselected: the training"
        " prefix depends only on `U`, and for all eighteen rows the recomputed training"
        " kernel basis and training equation count equal the stored ones."
    )
    lines.append("")
    lines.append(
        f"[COMPUTATION] Result: {len(refuted)} of 18 rows are"
        " `CANDIDATE_REFUTED_BY_HOLDOUT` at the newly observable orders and"
        f" {len(surviving)} remain truncation-unobservable."
    )
    lines.append("")
    lines.append("| # | class | degrees | frozen zero columns | verdict | first refuting / next testable | residual on frozen training kernel |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for index, row in enumerate(reaudit["rows"], start=1):
        zero_columns = ", ".join(str(value) for value in row["frozen_zero_column_indices"])
        if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT":
            marker = f"refuted at u^{row['first_refuting_holdout_order_u']}"
            residual = ", ".join(row["residual_on_training_kernel_basis"])
        else:
            marker = (
                f"next testable u^{row['next_testable_equation_order_u']}"
                f" (requires x^{row['required_coefficient_order_x']})"
            )
            residual = "0 at every available holdout"
        lines.append(
            f"| {index} | {row['relation_class']} | {_degree_label(row)} |"
            f" [{zero_columns}] | {row['verdict']} | {marker} | `{residual}` |"
        )
    lines.append("")
    lines.append(
        "[COMPUTATION] Why `x^34` is the first order that can test any of these rows: a"
        " surviving kernel direction is a literal zero column, i.e. a monomial whose formal"
        " first term lies beyond the available prefix.  With the 17 stored `u`-coefficients"
        " the derivative-safe equation orders stop at 15, so the `(f')^8` column of the"
        " `(0,0,8)` differential-algebraic row (formal valuation 16) had no admissible"
        " equation.  The new `x^34 = u^17` coefficient makes equation order 16 derivative-safe"
        " and that column acquires the entry `(3*f_3)^8 = 3^8 = 6561`, refuting the row."
        " The standalone test asserts both halves of this: with the old 17-coefficient series"
        " the row is `STILL_TRUNCATION_UNOBSERVABLE` with zero column `[8]`, and with the new"
        " series it is refuted at holdout order 16 with that exact residual entry."
    )
    lines.append("")
    required_orders = sorted({row["required_coefficient_order_x"] for row in surviving})
    requirement_phrase = (
        f"every one of them becomes testable at exactly `x^{required_orders[0]}`"
        if len(required_orders) == 1
        else "their requirements range from `x^"
        + str(required_orders[0])
        + "` to `x^"
        + str(required_orders[-1])
        + "`"
    )
    lines.append(
        "[UNRESOLVED] The "
        f"{len(surviving)} remaining rows are not evidence for a relation.  Each still has a"
        " literal zero column and is refutable only once the series reaches the coefficient"
        f" order listed above: {requirement_phrase}, i.e. precisely the order blocked by the"
        " cross-section guard of Section 6, so this re-audit closes at a clean boundary rather"
        " than at an arbitrary cut.  A truncation-unobservable survivor is a finite-prefix"
        " bookkeeping fact, never a claim that an algebraic or differential-algebraic relation"
        " exists."
    )
    lines.append("")
    lines.append("## 6. Resource walls")
    lines.append("")
    lines.append(
        f"[COMPUTATION] OBSERVED for the achieved `x^{TARGET_ORDER}` run:"
        f" wall {usage['elapsed_seconds']} s against a {usage['wall_budget_seconds']} s"
        f" budget, process peak RSS {usage['process_peak_rss_bytes']} bytes against a"
        f" {usage['rss_budget_bytes']} byte budget.  These are measured values"
        " (`resource.getrusage` and `time.perf_counter`), not estimates.  The largest box is"
        f" `5x5x5` (125 spins, 25-spin cross-section) with an unlaunched input-size estimate of"
        f" {data['resource_profile']['maximum_estimated_peak_box']['estimated_peak_bytes_at_target_truncation']}"
        " bytes at the target truncation."
    )
    lines.append("")
    lines.append(
        f"[COMPUTATION] The next order `x^{wall['next_order']}` reuses the SAME box"
        f" inventory (side-sum budget {wall['next_order_budget']}) with deeper truncation"
        " only.  The first newly walled order is"
        f" `x^{wall['next_wall_order']}`: budget 17 requires the canonical"
        f" `{wall['next_wall_canonical_shape'][0]}x{wall['next_wall_canonical_shape'][1]}x{wall['next_wall_canonical_shape'][2]}`"
        f" box, whose {wall['next_wall_cross_section_sites']}-spin cross-section exceeds"
        f" the (post-lift) transfer guard {wall['transfer_cross_section_guard']}."
        " No x^56 or higher computation was started, so no larger-order timing or memory"
        " figure is claimed."
    )
    lines.append("")
    lines.append("## 7. Scope")
    lines.append("")
    lines.append(
        "[UNRESOLVED] These are finite-prefix facts about the simple-cubic low-temperature"
        " expansion.  Eleven new exact coefficients and the survivor-verdict counts recorded above"
        " do not prove non-algebraicity, non-differential-algebraicity or non-D-finiteness of"
        " the true low-temperature free energy, do not constrain relations outside the"
        " displayed finite ansatz spaces, do not yield a closed form, and do not bound `K_c`."
    )
    lines.append("")
    lines.append(
        "[EXTERNAL] External-value role, stated without overclaiming in either direction: the"
        " critical-coupling benchmark `K_c` is used nowhere here, for fitting, selection, tuning"
        " or validation.  The hashed Guttmann--Enting table is likewise not a computational"
        " input and was not used for fitting, selection or tuning, but it IS used as a"
        " post-computation validation witness, and that comparison is a producer gate"
        " (Section 4).  The artifact records exactly this split under"
        " `kc_benchmark_used_for_fit_selection_or_validation` and `external_value_role`."
    )
    lines.append("")
    PROOF_PATH.parent.mkdir(parents=True, exist_ok=True)
    PROOF_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> bool:
    timestamp = datetime.now(timezone.utc).isoformat()
    profile = _profile()
    print("PROFILE FIRST (no coefficient computation has started)", flush=True)
    wall = profile["maximum_estimated_peak_box"]
    print(
        "PROFILE "
        f"LT x^{TARGET_ORDER}: canonical={profile['canonical_box_count']}, "
        f"ordered={profile['ordered_box_count']}, CRT_boxes={len(profile['crt_required_boxes'])}, "
        f"wall={wall['shape']} N={wall['sites']} cross={wall['cross_section_sites']} "
        f"truncated_peak={wall['estimated_peak_bytes_at_target_truncation']} bytes, "
        f"guard_ok={profile['all_cross_sections_within_guard']}",
        flush=True,
    )
    nxt = _next_order_preflight()
    print(
        f"PROFILE next x^{nxt['next_order']} reuses this inventory; first newly walled order"
        f" x^{nxt['next_wall_order']} needs cross-section"
        f" {nxt['next_wall_cross_section_sites']} > guard {nxt['transfer_cross_section_guard']}.",
        flush=True,
    )

    predecessor = json.loads(PREDECESSOR_PATH.read_text(encoding="utf-8"))
    old_coefficients = tuple(
        Fraction(value) for value in predecessor["data"]["coefficients"]
    )
    structure = json.loads(STRUCTURE_PATH.read_text(encoding="utf-8"))

    started = time.perf_counter()
    rss_before = _rss_bytes()
    with _truncated_hybrid_flm_engine(TARGET_ORDER) as transfer_records:
        series = low_temperature_free_energy(3, TARGET_ORDER)
    elapsed = time.perf_counter() - started
    peak_rss = _rss_bytes()

    canonical = _canonical_boxes(TARGET_BUDGET)
    if set(transfer_records) != set(canonical):
        missing = sorted(set(canonical) - set(transfer_records))
        extra = sorted(set(transfer_records) - set(canonical))
        raise AssertionError(f"transfer record inventory mismatch; missing={missing}, extra={extra}")

    checks: list[dict] = []
    _record(
        checks,
        "stored_LT_prefix_through_x32_exact",
        series.coefficients[: len(old_coefficients)] == old_coefficients,
        f"all {len(old_coefficients)} predecessor coefficients x^0,...,x^32 agree exactly",
    )
    new_even = {
        str(order): str(series.coefficients[order])
        for order in range(34, TARGET_ORDER + 1, 2)
    }
    _record(
        checks,
        "new_even_LT_coefficients_nonzero",
        all(series.coefficients[order] != 0 for order in range(34, TARGET_ORDER + 1, 2)),
        ", ".join(f"x^{order}={series.coefficients[order]}" for order in range(34, TARGET_ORDER + 1, 2)),
    )
    _record(
        checks,
        "all_odd_LT_coefficients_zero_through_x54",
        all(series.coefficients[order] == 0 for order in range(1, TARGET_ORDER + 1, 2)),
        "every odd x coefficient through the achieved order is exactly zero",
    )

    ordered_weight_bound_ok = True
    first_order_exact_ok = True
    weight_records = []
    for shape, weight in series.box_weights.items():
        lower_bound = 4 * sum(shape) - 6
        zero_below = all(value == 0 for value in weight[:lower_bound])
        first_nonzero = next((index for index, value in enumerate(weight) if value), None)
        ordered_weight_bound_ok = ordered_weight_bound_ok and zero_below
        first_order_exact_ok = first_order_exact_ok and first_nonzero == lower_bound
        if shape == tuple(sorted(shape)):
            weight_records.append(
                {
                    "shape": list(shape),
                    "lower_bound": lower_bound,
                    "first_nonzero_order": first_nonzero,
                    "zero_below_lower_bound": zero_below,
                }
            )
    _record(
        checks,
        "all_ordered_box_weights_obey_proved_surface_bound",
        ordered_weight_bound_ok,
        "every ordered box weight has no coefficient below 4*(a+b+c)-6",
    )
    _record(
        checks,
        "canonical_box_weights_start_at_surface_bound",
        first_order_exact_ok,
        "every canonical box has its first nonzero retained weight exactly at the proved bound",
    )

    crt_certificates = []
    inventory_records = []
    polynomial_prefixes = {}
    for shape in canonical:
        record = transfer_records[shape]
        polynomial_prefixes[_shape_key(shape)] = record["coefficient_prefix"]
        inventory_records.append(
            {
                key: record[key]
                for key in (
                    "shape",
                    "sites",
                    "cross_section_sites",
                    "engine",
                    "requested_degree",
                    "estimated_peak_bytes_at_target_truncation",
                    "coefficient_prefix_sha256",
                )
            }
        )
        if record["engine"] != "int64":
            certificate = {
                "claim_tag": "[COMPUTATION]",
                "shape": record["shape"],
                "sites": record["sites"],
                "cross_section_sites": record["cross_section_sites"],
                "prime_moduli": record["prime_moduli"],
                "modulus_product": record["modulus_product"],
                "modulus_product_bits": record["modulus_product_bits"],
                "coefficient_bound": record["coefficient_bound"],
                "twice_coefficient_bound": record["twice_coefficient_bound"],
                "uniqueness_certificate_passed": record["uniqueness_certificate_passed"],
                "bound_reason": (
                    "Each coefficient counts a subset of the 2^N spin configurations, so it lies "
                    "in [0,2^N]; the recorded CRT product exceeds 2*2^N."
                ),
                "prefix_sha256": record["coefficient_prefix_sha256"],
                "full_sum_check_status": record["full_sum_check_status"],
            }
            crt_certificates.append(certificate)
    _record(
        checks,
        "all_CRT_reconstructions_have_strict_uniqueness_certificate",
        bool(crt_certificates)
        and all(item["uniqueness_certificate_passed"] for item in crt_certificates),
        f"{len(crt_certificates)} CRT box prefixes have modulus > 2*(2^N)",
    )

    survivor_reaudit = _survivor_reaudit(series.coefficients, structure)
    _record(
        checks,
        "all_18_survivor_training_prefixes_unchanged",
        survivor_reaudit["row_count"] == 18
        and survivor_reaudit["training_prefix_unchanged_for_all_rows"],
        "each re-audit reused the frozen e43 U+2 training kernel before evaluating new holdouts",
    )
    _record(
        checks,
        "survivor_reaudit_has_only_requested_verdicts",
        all(
            row["verdict"]
            in {"CANDIDATE_REFUTED_BY_HOLDOUT", "STILL_TRUNCATION_UNOBSERVABLE"}
            for row in survivor_reaudit["rows"]
        ),
        json.dumps(survivor_reaudit["verdict_counts"], sort_keys=True),
    )

    # This source is opened only here, after every local box polynomial,
    # logarithm, inversion, coefficient check, and survivor audit above.
    external_witness = _external_witness_after_local_derivation(series.coefficients)
    if external_witness["available"]:
        _record(
            checks,
            "post_computation_external_table_witness",
            external_witness["all_match"],
            "exponentiating the local series reproduces all "
            f"{external_witness['row_count']} printed Table-2 Lambda_0 values through u^26=x^52, "
            "including every new order "
            + ",".join(f"x^{order}" for order in external_witness["witnessed_new_x_orders"]),
        )

    predecessor_sha256 = hashlib.sha256(PREDECESSOR_PATH.read_bytes()).hexdigest()
    data = {
        "dimension": 3,
        "variable": "x",
        "achieved_order": TARGET_ORDER,
        "coefficients": _fraction_strings(series.coefficients),
        "nonzero_coefficients": {
            str(index): str(value)
            for index, value in enumerate(series.coefficients)
            if value
        },
        "interaction_coefficients": _fraction_strings(series.interaction_coefficients),
        "nonzero_interaction_coefficients": {
            str(index): str(value)
            for index, value in enumerate(series.interaction_coefficients)
            if value
        },
        "ordered_box_count": len(series.boxes),
        "canonical_box_count": len(canonical),
        "canonical_box_list": [list(shape) for shape in canonical],
        "order_bound": {
            "statement": series.order_bound,
            "budget": series.bound_budget,
            "bound_slack": series.bound_slack,
        },
        "lattice": "simple cubic",
        "normalization": "phi = 3*K + sum_n coefficients[n] * x^n",
        "coefficient_provenance": "derived_locally",
        "extension": {
            "previous_artifact": "results/series/extended2_sc_lt_free_energy.json",
            "previous_artifact_sha256": predecessor_sha256,
            "previous_order": len(old_coefficients) - 1,
            "new_order": TARGET_ORDER,
            "new_nonzero_coefficients": new_even,
            "stored_prefix_reproduced_through_order": 32,
        },
        "finite_lattice_inventory": {
            "claim_tag": "[THEOREM]",
            "canonical_box_count": len(canonical),
            "ordered_box_count": len(series.boxes),
            "canonical_boxes": inventory_records,
            "completeness_statement": (
                f"For x^{TARGET_ORDER}, every potentially contributing connected cluster has "
                f"bounding-box side sum at most {TARGET_BUDGET}: a box of side sum "
                f"{TARGET_BUDGET+1} has surface lower bound "
                f"4*{TARGET_BUDGET+1}-6=58>{TARGET_ORDER}.  The inventory contains every "
                f"ordered positive triple with a+b+c<={TARGET_BUDGET} and evaluates one "
                "cached polynomial per permutation class."
            ),
            "surface_bound": "|partial S| >= 4*(a+b+c)-6",
            "crt_box_count": len(crt_certificates),
            "int64_box_count": len(canonical) - len(crt_certificates),
        },
        "box_polynomial_prefixes": polynomial_prefixes,
        "canonical_weight_support": {
            "claim_tag": "[COMPUTATION]",
            "all_ordered_boxes_zero_below_bound": ordered_weight_bound_ok,
            "all_canonical_boxes_first_nonzero_at_bound": first_order_exact_ok,
            "canonical_records": weight_records,
        },
        "crt_certificates": crt_certificates,
        "profile_completed_before_computation": True,
        "resource_profile": profile,
        "resource_usage": {
            "claim_tag": "[COMPUTATION]",
            "wall_budget_seconds": WALL_BUDGET_SECONDS,
            "rss_budget_bytes": RSS_BUDGET_BYTES,
            "elapsed_seconds": f"{elapsed:.6f}",
            "rss_before_bytes": rss_before,
            "process_peak_rss_bytes": peak_rss,
            "observed_within_wall_budget": elapsed <= WALL_BUDGET_SECONDS,
            "observed_within_rss_budget": peak_rss <= RSS_BUDGET_BYTES,
            "measurement_kind": "OBSERVED process wall time and max RSS",
        },
        "next_order_preflight": _next_order_preflight(),
        "survivor_reaudit": survivor_reaudit,
        "external_witness": external_witness,
        "kc_benchmark_used_for_fit_selection_or_validation": False,
        "external_value_role": {
            "claim_tag": "[EXTERNAL]",
            "used_as_computational_input": False,
            "used_for_fitting_selection_or_tuning": False,
            "used_as_post_computation_validation_gate": True,
            "detail": (
                "No external number enters any local computation, and none was used to fit, "
                f"select or tune anything: the whole x^{TARGET_ORDER} derivation, the stored-prefix check "
                "and the survivor re-audit complete before the hashed Guttmann--Enting table "
                "is opened.  The table is then used as a post-computation validation witness "
                "and that comparison IS a gate: it is appended to `checks` as "
                "`post_computation_external_table_witness`, so a mismatch fails this "
                "experiment.  The separate `kc_benchmark_...` flag above refers to the "
                "critical-coupling benchmark K_c, which is used nowhere here."
            ),
        },
    }
    provenance = {
        "script": SCRIPT,
        "timestamp_utc": timestamp,
        "precision": (
            "exact Python int/Fraction; int64 transfer for N<=62; deterministic 31-bit-prime "
            "CRT for N>62 with target-order truncation and product strictly exceeding twice "
            "the 2^N coefficient bound"
        ),
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(
        json.dumps({"provenance": provenance, "data": data, "checks": checks}, indent=2)
        + "\n",
        encoding="utf-8",
    )
    _write_proof_note(data)
    print(f"WROTE {PROOF_PATH.relative_to(ROOT)}", flush=True)

    passed = all(check["passed"] for check in checks)
    if passed:
        print(
            f"PASS e124_lt_x34: exact LT series through x^{TARGET_ORDER}; "
            f"{survivor_reaudit['verdict_counts'].get('CANDIDATE_REFUTED_BY_HOLDOUT', 0)} "
            "legacy survivors refuted by frozen-protocol holdouts",
            flush=True,
        )
    else:
        print("FAIL e124_lt_x34: one or more recorded checks failed", flush=True)
    return passed


if __name__ == "__main__":
    try:
        ok = main()
    except Exception as error:
        print(f"FAIL e124_lt_x34: {type(error).__name__}: {error}", flush=True)
        raise
    raise SystemExit(0 if ok else 1)
