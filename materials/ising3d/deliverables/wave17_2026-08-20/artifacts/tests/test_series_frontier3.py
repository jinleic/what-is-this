"""Independent exact replay of the wave-11 third-frontier holdout refutation.

Everything decisive is rebuilt from raw artifact data:

* the canonical exact HT reduced-free-energy coefficients `f_0..f_28` are read
  from `results/series/ht_v28.json` (cross-checked against `ht_v24.json`,
  `ht_v26.json`, and the frozen e56 training prefix);
* the frozen e56 candidate (sole HT truncation-unobservable Euler-ODE row) is
  read from `results/series/frontier2.json` as *input*, never as a conclusion;
* no code is imported from `experiments/e113_series_frontier3.py` or from any
  producer; all linear algebra is an independent exact `Fraction`
  implementation inside this file.

Run: `PYTHONPATH=src .venv/bin/python tests/test_series_frontier3.py`.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "series" / "frontier3.json"
FRONTIER2 = ROOT / "results" / "series" / "frontier2.json"
HT24 = ROOT / "results" / "series" / "ht_v24.json"
HT26 = ROOT / "results" / "series" / "ht_v26.json"
HT28 = ROOT / "results" / "series" / "ht_v28.json"
EXT2 = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
E43 = ROOT / "results" / "series" / "structure_certificates.json"

ALLOWED_VERDICTS = {
    "NO_RELATION_AT_BUDGET",
    "CANDIDATE_REFUTED_BY_HOLDOUT",
    "CANDIDATE_SURVIVES",
}

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# Independent exact linear algebra over Q.
# --------------------------------------------------------------------------


def _rref(rows: Sequence[Sequence[Fraction]], column_count: int):
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivot_row = 0
    pivot_columns: list[int] = []
    for column in range(column_count):
        selected = None
        for row in range(pivot_row, len(matrix)):
            if matrix[row][column] != 0:
                selected = row
                break
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row != pivot_row and matrix[row][column] != 0:
                factor = matrix[row][column]
                matrix[row] = [
                    a - factor * b
                    for a, b in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix, pivot_columns


def _rank(rows: Sequence[Sequence[Fraction]], column_count: int) -> int:
    return len(_rref(rows, column_count)[1])


def _nullspace(
    rows: Sequence[Sequence[Fraction]], column_count: int
) -> list[list[Fraction]]:
    matrix, pivot_columns = _rref(rows, column_count)
    free_columns = [
        column for column in range(column_count) if column not in pivot_columns
    ]
    basis = []
    for free in free_columns:
        vector = [Fraction(0)] * column_count
        vector[free] = Fraction(1)
        for row, column in enumerate(pivot_columns):
            vector[column] = -matrix[row][free]
        basis.append(vector)
    return basis


def _determinant(rows: Sequence[Sequence[Fraction]]) -> Fraction:
    matrix = [[Fraction(value) for value in row] for row in rows]
    size = len(matrix)
    sign = 1
    for index in range(size):
        selected = None
        for row in range(index, size):
            if matrix[row][index] != 0:
                selected = row
                break
        if selected is None:
            return Fraction(0)
        if selected != index:
            matrix[index], matrix[selected] = matrix[selected], matrix[index]
            sign = -sign
        pivot = matrix[index][index]
        for row in range(index + 1, size):
            factor = matrix[row][index] / pivot
            if factor:
                matrix[row] = [
                    a - factor * b for a, b in zip(matrix[row], matrix[index])
                ]
    result = Fraction(sign)
    for index in range(size):
        result *= matrix[index][index]
    return result


def _primitive(vector: Sequence[Fraction]) -> list[int]:
    denominator = 1
    for value in vector:
        denominator = (
            denominator * value.denominator
            // math.gcd(denominator, value.denominator)
        )
    integers = [int(value * denominator) for value in vector]
    divisor = math.gcd(*[abs(value) for value in integers])
    integers = [value // divisor for value in integers]
    for value in integers:
        if value:
            if value < 0:
                integers = [-item for item in integers]
            break
    return integers


# --------------------------------------------------------------------------
# Independent rebuilders of the two ansatz matrix conventions.
# --------------------------------------------------------------------------


def _euler_entry(
    order: int, degree_t: int, theta_power: int, coefficients: Sequence[Fraction]
) -> Fraction:
    """[v^order] (v^degree_t * theta^theta_power F) with theta = v d/dv."""
    shifted = order - degree_t
    if shifted < 0:
        return Fraction(0)
    return Fraction(shifted) ** theta_power * coefficients[shifted]

def _euler_row(
    order: int, monomials: Sequence[Sequence[int]], coefficients
) -> list[Fraction]:
    return [
        _euler_entry(order, degree_t, theta_power, coefficients)
        for degree_t, theta_power in monomials
    ]


def _euler_residual(
    order: int,
    monomials: Sequence[Sequence[int]],
    kernel: Sequence[int],
    coefficients: Sequence[Fraction],
) -> Fraction:
    return sum(
        Fraction(weight) * entry
        for weight, entry in zip(kernel, _euler_row(order, monomials, coefficients))
    )


def _convolve(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(order + 1)]
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if i + j <= order:
                    result[i + j] += a * b
    return result


def _powers(series: Sequence[Fraction], maximum: int, order: int):
    values = [(Fraction(1),) + (Fraction(0),) * order]
    for _ in range(maximum):
        previous = values[-1]
        values.append(tuple(_convolve(previous, series, order)))
    return values


def _da_matrix(
    series: Sequence[Fraction], degree_x: int, degree_f: int, degree_derivative: int
):
    """e43 convention: Q(t,f,f') = sum q_abc t^a f^b (f')^c, c >= 1 included."""
    equation_count = len(series) - 1
    truncated = tuple(series[:equation_count])
    derivative = tuple(
        Fraction(order + 1) * series[order + 1] for order in range(equation_count)
    )
    f_powers = _powers(truncated, degree_f, equation_count - 1)
    d_powers = _powers(derivative, degree_derivative, equation_count - 1)
    monomials = [
        (a, b, c)
        for c in range(degree_derivative + 1)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    columns = []
    for a, b, c in monomials:
        product = _convolve(f_powers[b], d_powers[c], equation_count - 1)
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
    return rows, monomials


def _da_holdout_ladder(rows, monomials):
    """e43/e56 budget rule: train on orders 0..U+1, hold out the rest."""
    unknown_count = len(monomials)
    training_count = unknown_count + 2
    training = rows[:training_count]
    training_rank = _rank(training, unknown_count)
    basis = _nullspace(training, unknown_count)
    first_refuting = None
    ladder = []
    for order in range(training_count, len(rows)):
        residuals = [
            sum(a * b for a, b in zip(rows[order], vector)) for vector in basis
        ]
        constraint_rank = _rank(
            [list(item) for item in ladder + [residuals]], len(basis)
        ) if basis else 0
        surviving = len(basis) - constraint_rank
        ladder.append(residuals)
        if surviving == 0 and first_refuting is None:
            first_refuting = order
    full_rank = _rank(rows, unknown_count)
    if training_rank == unknown_count:
        verdict = "NO_RELATION_AT_BUDGET"
    elif full_rank == unknown_count:
        verdict = "CANDIDATE_REFUTED_BY_HOLDOUT"
    else:
        verdict = "CANDIDATE_SURVIVES"
    return {
        "training_rank": training_rank,
        "training_nullity": unknown_count - training_rank,
        "first_refuting_holdout_order": first_refuting,
        "full_rank": full_rank,
        "verdict": verdict,
    }


# --------------------------------------------------------------------------
# Replay.
# --------------------------------------------------------------------------


def main() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    data = result["data"]
    frontier2 = json.loads(FRONTIER2.read_text(encoding="utf-8"))["data"]
    ht28 = json.loads(HT28.read_text(encoding="utf-8"))["data"]
    ht24 = json.loads(HT24.read_text(encoding="utf-8"))["data"]
    ht26 = json.loads(HT26.read_text(encoding="utf-8"))["data"]
    ext2 = json.loads(EXT2.read_text(encoding="utf-8"))["data"]
    e43 = json.loads(E43.read_text(encoding="utf-8"))["data"]

    # --- 1. provenance: recorded hashes match the actual source artifacts ----
    recorded = {
        entry["path"]: entry["sha256"] for entry in result["provenance"]["sources"]
    }
    actual = {
        "results/series/frontier2.json": _sha256(FRONTIER2),
        "results/series/ht_v28.json": _sha256(HT28),
        "results/series/ht_v24.json": _sha256(HT24),
        "results/series/ht_v26.json": _sha256(HT26),
        "results/series/extended2_sc_ht_free_energy.json": _sha256(EXT2),
        "results/series/structure_certificates.json": _sha256(E43),
    }
    _check(
        "source_hashes",
        recorded == actual,
        f"all {len(actual)} recorded canonical source SHA-256 digests match the files",
    )
    _check(
        "no_benchmark_use",
        result["provenance"]["benchmark_used_for_fit_selection_or_validation"] is False,
        "artifact records that no benchmark value entered selection or validation",
    )

    # --- 2. canonical coefficient chain through v^28 ------------------------
    coefficients = [
        Fraction(value) for value in ht28["series"]["coefficients"]
    ]
    _check(
        "canonical_length",
        len(coefficients) == 29 and ht28["series"]["achieved_order"] == 28,
        "ht_v28.json supplies 29 exact coefficients f_0..f_28",
    )
    _check(
        "e56_training_prefix_unchanged",
        [Fraction(value) for value in ext2["coefficients"]] == coefficients[:23]
        and [Fraction(value) for value in frontier2["input_series"]["HT"]["coefficients"]]
        == coefficients[:23],
        "the frozen e56 training prefix v^0..v^22 is byte-identical in all three artifacts",
    )
    _check(
        "holdout_coefficient_crosscheck",
        Fraction(ht24["series"]["v24"]) == coefficients[24]
        and Fraction(ht26["series"]["v26"]) == coefficients[26],
        "v^24 and v^26 agree between their own artifacts and the canonical v^28 list",
    )

    # --- 3. the frozen candidate, taken from frontier2.json as input ---------
    candidates = [
        row
        for row in frontier2["frontiers"]["linear_euler_ode_order_at_most_2"]["HT"][
            "frontier"
        ]
        if row["verdict"] == "TRUNCATION_UNOBSERVABLE"
    ]
    _check("sole_frozen_survivor", len(candidates) == 1, "exactly one frozen e56 HT Euler-ODE survivor row")
    frozen = candidates[0]
    monomials = [tuple(item) for item in frozen["monomials"]]
    kernel = list(frozen["full_kernel_basis"][0])
    _check(
        "candidate_identity",
        frozen["parameters"] == {"order": 2, "polynomial_degree": 6}
        and len(monomials) == 21
        and data["candidate_provenance"]["parameters"]
        == frozen["parameters"],
        "candidate is the order-2 degree-6 Euler row with 21 monomial columns",
    )
    _check(
        "first_unobservable_order_is_24",
        frozen["support_components"][0]["first_unfixed_supported_order"] == 24,
        "the frozen support certificate names v^24 as the first unfixed supported order",
    )

    # --- 4. independent kernel reconstruction on the training prefix --------
    training_rows = [
        _euler_row(order, monomials, coefficients) for order in range(23)
    ]
    training_rank = _rank(training_rows, 21)
    basis = _nullspace(training_rows, 21)
    _check(
        "training_prefix_rank",
        training_rank == 20 and len(basis) == 1,
        f"independent elimination gives training rank {training_rank}, nullity {len(basis)}",
    )
    rebuilt = _primitive(basis[0])
    _check(
        "kernel_reconstruction",
        rebuilt == kernel,
        "independently reconstructed primitive kernel equals the frozen vector exactly",
    )
    _check(
        "frozen_kernel_satisfies_training",
        all(
            _euler_residual(order, monomials, kernel, coefficients) == 0
            for order in range(23)
        ),
        "frozen kernel annihilates every training order 0..22 (rebuilt)",
    )

    # --- 5. decisive holdout evaluation --------------------------------------
    holdout = data["holdout_evaluation"]
    residuals = {
        entry["order"]: entry["residual"]
        for entry in holdout["residuals"]
    }
    replayed = {
        order: _euler_residual(order, monomials, kernel, coefficients)
        for order in range(23, 29)
    }
    _check(
        "holdout_residuals_replay",
        all(Fraction(residuals[order]) == replayed[order] for order in range(23, 29)),
        "all six holdout residuals replay exactly from the raw coefficient list",
    )
    _check(
        "odd_orders_parity_zero",
        all(replayed[order] == 0 for order in (23, 25, 27)),
        "odd holdout orders are exactly zero (even effective support, even series support)",
    )
    _check(
        "decisive_residual_nonzero",
        replayed[24] != 0
        and holdout["decisive_order"] == 24
        and Fraction(holdout["first_nonzero_residual"]) == replayed[24],
        f"[v^24] residual = {replayed[24]} != 0",
    )
    _check(
        "confirmation_only",
        holdout["confirmation_only_orders"] == [26, 28]
        and replayed[26] != 0
        and replayed[28] != 0,
        "v^26 and v^28 residuals are nonzero and flagged confirmation-only",
    )

    # surviving-dimension ladder: dimension 1 through order 23, killed at 24.
    surviving = [1]
    for order in range(23, 29):
        residual = replayed[order]
        surviving.append(0 if (residual != 0 and surviving[-1] == 1) else surviving[-1])
    _check(
        "surviving_dimension_ladder",
        surviving == [1, 1, 0, 0, 0, 0, 0],
        "kernel line survives order 23 and is first killed at order 24",
    )

    # --- 6. forced-value decomposition ---------------------------------------
    forced = holdout["forced_value_decomposition"]
    # S collects exactly the columns that touch f_24, i.e. degree_t == 0.
    big_s = sum(
        Fraction(weight) * Fraction(24) ** theta_power
        for weight, (degree_t, theta_power) in zip(kernel, monomials)
        if degree_t == 0
    )
    little_t = sum(
        Fraction(weight)
        * _euler_entry(24, degree_t, theta_power, coefficients)
        for weight, (degree_t, theta_power) in zip(kernel, monomials)
        if degree_t > 0
    )
    _check(
        "forced_value_decomposition",
        Fraction(forced["S"]) == big_s
        and Fraction(forced["T"]) == little_t
        and big_s * coefficients[24] + little_t == replayed[24]
        and Fraction(forced["forced_v24"]) == -little_t / big_s
        and Fraction(forced["actual_v24"]) == coefficients[24],
        "residual = S*f_24 + T with the stored exact S, T and forced value",
    )

    # --- 7. nonzero maximal-minor refutation certificate ----------------------
    certificate = holdout["certificate"]
    full_rows = [
        _euler_row(order, monomials, coefficients) for order in range(29)
    ]
    minor = [
        [full_rows[order][column] for column in certificate["column_indices"]]
        for order in certificate["row_orders"]
    ]
    determinant = _determinant(minor)
    _check(
        "minor_certificate_replay",
        determinant == Fraction(certificate["determinant"]) and determinant != 0,
        f"stored {len(certificate['row_orders'])}x{len(certificate['column_indices'])} minor replays to a nonzero exact determinant",
    )
    _check(
        "full_system_rank",
        _rank(full_rows, 21) == 21 == holdout["full_system"]["rank"],
        "all 29 exact equations have full column rank 21 over Q",
    )

    # --- 8. classification semantics ------------------------------------------
    _check(
        "semantics_not_no_relation",
        21 - training_rank == 1,
        "training nullity is 1, so the row is not NO_RELATION_AT_BUDGET",
    )
    _check(
        "semantics_not_survives",
        _rank(full_rows, 21) == 21 and replayed[24] != 0,
        "full-rank system plus nonzero holdout residual excludes CANDIDATE_SURVIVES",
    )
    _check(
        "stored_verdict",
        holdout["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
        and holdout["first_refuting_holdout_order"] == 24,
        "stored decisive verdict is CANDIDATE_REFUTED_BY_HOLDOUT at order 24",
    )

    def _verdict_values(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "verdict":
                    yield value
                yield from _verdict_values(value)
        elif isinstance(node, list):
            for item in node:
                yield from _verdict_values(item)

    verdict_values = sorted(set(_verdict_values(data)))
    _check(
        "verdict_vocabulary",
        set(verdict_values) <= ALLOWED_VERDICTS
        and "CANDIDATE_SURVIVES" not in verdict_values,
        f"verdict labels {verdict_values} stay inside the three-way vocabulary",
    )
    _check(
        "branch_taken",
        data["holdout_protocol"]["branch_taken"] == "REFUTED_AT_FIRST_UNOBSERVABLE_ORDER"
        and data["holdout_protocol"]["predeclared_expanded_prefix_rank_grid"]["status"]
        == "NOT_EXECUTED",
        "predeclared decision tree took the refutation branch; expanded grid not executed",
    )

    # --- 9. e43 strict holdout controls ---------------------------------------
    controls = {
        (control["parameters"]["degree_f"], control["parameters"]["degree_f_prime"]): control
        for control in data["controls"]["e43_strict_holdout_rows"]
    }
    lt_reduced = [
        Fraction(value)
        for value in e43["input_series"]["LT"]["reduced_coefficients"]
    ]
    strict_rows = [
        row
        for row in e43["differential_algebraicity_order_1"]["LT"]["frontier"]
        if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ]
    matched = 0
    for row in strict_rows:
        rows, mons = _da_matrix(
            lt_reduced,
            row["degree_x"],
            row["degree_f"],
            row["degree_f_prime"],
        )
        _check(
            f"e43_monomials_{row['degree_f']}_{row['degree_f_prime']}",
            [list(item) for item in mons] == row["monomials"],
            f"monomial enumeration matches for (d_f={row['degree_f']}, d_f'={row['degree_f_prime']})",
        )
        ladder = _da_holdout_ladder(rows, mons)
        stored = controls[(row["degree_f"], row["degree_f_prime"])]
        subminor = [
            [rows[order][column] for column in row["certificate"]["column_indices"]]
            for order in row["certificate"]["row_orders"]
        ]
        determinant = _determinant(subminor)
        ok = (
            ladder["first_refuting_holdout_order"] == row["first_refuting_holdout_order"]
            and ladder["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
            and determinant == Fraction(row["certificate"]["determinant"])
            and stored["recomputed_first_refuting_holdout_order"]
            == ladder["first_refuting_holdout_order"]
            and Fraction(stored["recomputed_determinant"]) == determinant
        )
        matched += int(ok)
        _check(
            f"e43_control_{row['degree_f']}_{row['degree_f_prime']}",
            ok,
            f"strict LT control reproduces order {row['first_refuting_holdout_order']} and its nonzero minor",
        )
    _check(
        "e43_controls_complete",
        matched == len(strict_rows) == 8,
        f"all {matched} frozen e43 strict holdout refutations replay exactly",
    )

    ht_reduced = [
        Fraction(value)
        for value in e43["input_series"]["HT"]["reduced_coefficients"]
    ]
    no_relation = e43["differential_algebraicity_order_1"]["HT"]["frontier"][0]
    rows, mons = _da_matrix(
        ht_reduced,
        no_relation["degree_x"],
        no_relation["degree_f"],
        no_relation["degree_f_prime"],
    )
    ladder = _da_holdout_ladder(rows, mons)
    subminor = [
        [
            rows[order][column]
            for column in no_relation["certificate"]["column_indices"]
        ]
        for order in no_relation["certificate"]["row_orders"]
    ]
    _check(
        "no_relation_control",
        ladder["verdict"] == "NO_RELATION_AT_BUDGET"
        and ladder["training_rank"] == len(mons)
        and _determinant(subminor)
        == Fraction(no_relation["certificate"]["determinant"]),
        "a frozen e43 full-rank row replays as NO_RELATION_AT_BUDGET with its minor",
    )

    # --- 10. envelope ----------------------------------------------------------
    _check(
        "all_producer_checks_passed",
        all(item["passed"] for item in result["checks"]),
        f"all {len(result['checks'])} producer-recorded checks passed",
    )
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        raise
