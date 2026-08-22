"""Third series-structure frontier: holdout refutation of the frozen e56 survivor.

The frozen wave-6/7 frontier `results/series/frontier2.json` left exactly one
HT candidate undecided: the order-two Euler-ODE row (polynomial degree six,
21 monomial columns) whose exact kernel was certified truncation-unobservable
because its first unfixed supported order, v^24, lay beyond the then-available
prefix. Waves 7-9 later computed the exact coefficients v^24 (e68), v^26
(e86), and v^28 (e94) without ever using them for structure fitting.

This experiment does the one thing those coefficients are reserved for:

1. load the frozen candidate and the canonical exact coefficient list f_0..f_28
   with full hash-chain provenance;
2. reconstruct the candidate's exact row/polynomial relation from the raw
   v^0..v^22 training prefix alone (no refit, no reselection);
3. evaluate the previously unobservable holdout orders 23..28 strictly in
   increasing order -- order 24 is decisive, orders 26/28 are recorded only as
   independent confirmation;
4. classify under the predeclared three-way semantics and store exact
   residuals, a nonzero maximal-minor certificate, forced-value decomposition,
   the full replay matrix, and re-audited e43 strict holdout controls.

All arithmetic is exact over Q. No critical benchmark enters any selection or
validation step, and no survival discovery or D-finiteness claim is made.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Sequence

SCRIPT = "experiments/e113_series_frontier3.py"
INTERPRETER = ".venv/bin/python"
ROOT = Path(__file__).resolve().parents[1]
FRONTIER2_PATH = ROOT / "results" / "series" / "frontier2.json"
EXT2_HT_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
HT24_PATH = ROOT / "results" / "series" / "ht_v24.json"
HT26_PATH = ROOT / "results" / "series" / "ht_v26.json"
HT28_PATH = ROOT / "results" / "series" / "ht_v28.json"
E43_PATH = ROOT / "results" / "series" / "structure_certificates.json"
RESULT_PATH = ROOT / "results" / "series" / "frontier3.json"
PROOF_PATH = ROOT / "proofs" / "series_frontier3.md"

E56_FREEZE_ORDER = 22          # last training order available to e56
CANONICAL_ORDER = 28           # last exact HT coefficient known
DECISIVE_ORDER = 24            # first unfixed supported order of the candidate
ALLOWED_VERDICTS = (
    "NO_RELATION_AT_BUDGET",
    "CANDIDATE_REFUTED_BY_HOLDOUT",
    "CANDIDATE_SURVIVES",
)

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

# Independent crosscheck values reported by the coordinating parent agent
# (evaluated from the artifacts before this script ran).  They are recorded
# for auditability only and are never inputs to any computation below.
PARENT_REPORTED_RESIDUALS = {
    24: "499752451221372610239349461368647432702945848808/11",
    26: "500071378398815869213964078020051854907899278078520/143",
    28: "22667554452529238570074088272879080635241164323868544/143",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_strings(values) -> list[str]:
    return [str(Fraction(value)) for value in values]


# --------------------------------------------------------------------------
# Exact linear algebra over Q (Fraction Gaussian elimination; the test file
# carries its own independent copy).
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


def _independent_row_orders(
    rows: Sequence[Sequence[Fraction]], column_count: int
) -> list[int]:
    """Deterministic sorted greedy selection of independent rows."""
    chosen: list[int] = []
    kept: list[list[Fraction]] = []
    current_rank = 0
    for order, row in enumerate(rows):
        candidate_rank = _rank(kept + [list(row)], column_count)
        if candidate_rank > current_rank:
            chosen.append(order)
            kept.append(list(row))
            current_rank = candidate_rank
        if current_rank == column_count:
            break
    return chosen


# --------------------------------------------------------------------------
# Euler-ODE ansatz matrix (e56 convention).
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


# --------------------------------------------------------------------------
# e43 first-order differential-algebraic ansatz (control rebuild).
# --------------------------------------------------------------------------


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
    """e43 convention: Q(t,f,f') = sum q_abc t^a f^b (f')^c."""
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


def _budget_ladder(rows, unknown_count: int) -> dict:
    """e43/e56 U+2 training rule with strict later-coefficient holdouts."""
    training_count = unknown_count + 2
    training = rows[:training_count]
    training_rank = _rank(training, unknown_count)
    basis = _nullspace(training, unknown_count)
    first_refuting = None
    constraints: list[list[Fraction]] = []
    for order in range(training_count, len(rows)):
        residuals = [
            sum(a * b for a, b in zip(rows[order], vector)) for vector in basis
        ]
        constraints.append(residuals)
        if basis:
            constraint_rank = _rank(constraints, len(basis))
        else:
            constraint_rank = 0
        if basis and len(basis) - constraint_rank == 0 and first_refuting is None:
            first_refuting = order
    full_rank = _rank(rows, unknown_count)
    if training_rank == unknown_count:
        verdict = "NO_RELATION_AT_BUDGET"
    elif full_rank == unknown_count:
        verdict = "CANDIDATE_REFUTED_BY_HOLDOUT"
    else:
        verdict = "CANDIDATE_SURVIVES"
    return {
        "training_count": training_count,
        "training_rank": training_rank,
        "training_nullity": unknown_count - training_rank,
        "first_refuting_holdout_order": first_refuting,
        "full_rank": full_rank,
        "full_nullity": unknown_count - full_rank,
        "verdict": verdict,
    }


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> bool:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)
    return bool(passed)


def _proof_text(payload: dict) -> str:
    data = payload["data"]
    holdout = data["holdout_evaluation"]
    candidate = data["candidate_provenance"]
    reconstruction = candidate["relation_reconstruction"]
    forced = holdout["forced_value_decomposition"]
    controls = data["controls"]
    q0, q1, q2 = (
        candidate["polynomial_coefficients"][key]
        for key in ("Q0_v_powers_0_to_6", "Q1_v_powers_0_to_6", "Q2_v_powers_0_to_6")
    )
    return f"""# Third-frontier holdout evaluation of the frozen e56 Euler-ODE survivor

## Exact inputs and frozen candidate

[COMPUTATION] The canonical input is the exact simple-cubic high-temperature
reduced free energy `phi = log 2 + sum_n f_n v^n` through `v^28`
(`results/series/ht_v28.json`, 29 exact rational coefficients). The stored
prefix agrees coefficient-for-coefficient with the frozen e56 training prefix
`results/series/extended2_sc_ht_free_energy.json` through `v^22`, and the new
holdout coefficients agree with their own production artifacts:
`f_24 = {data['input_series']['HT']['v24']}`
(`results/series/ht_v24.json`) and
`f_26 = {data['input_series']['HT']['v26']}`
(`results/series/ht_v26.json`). SHA-256 hashes of all six source artifacts
are recorded in `results/series/frontier3.json`.

[COMPUTATION] The evaluated candidate is the sole HT row of the frozen e56
frontier not closed at the `v^22` budget: the homogeneous Euler-form linear
ODE `sum_j Q_j(v) theta^j F = 0` with `theta = v d/dv`, order `r = 2` and
`deg Q_j <= 6` (21 monomial columns). Its e56 verdict was
TRUNCATION_UNOBSERVABLE with first unfixed supported order `v^24`. The
producer verified it is the unique HT truncation-unobservable row across all
three e56 HT families ({candidate['sole_survivor_audit']['ht_truncation_unobservable_rows_total']} row in total).

## Predeclared holdout protocol

[LEMMA] The orders `23..28` were unavailable at the e56 freeze, so under the
inherited U+2 discipline they are strict holdouts: the training prefix of the
candidate remains exactly the equations at orders `0..22`, and no coefficient
beyond `v^22` may enter any refit, reselection, or kernel reconstruction.

[COMPUTATION] The decision tree was fixed before any holdout evaluation:
(i) a nonzero residual at order 24 classifies
`CANDIDATE_REFUTED_BY_HOLDOUT` with `v^26`/`v^28` recorded as confirmation
only; (ii) a zero residual at 24 triggers sequential evaluation at 26 then
28; (iii) if all three vanish, a predeclared expanded prefix-rank grid
(cells `(r,d)` with `(r+1)(d+1)+2 <= 29`, training orders `0..U+1`, holdouts
`U+2..28` preserved) is recomputed with no survival discovery and no
D-finiteness claims. The branch taken was
`{data['holdout_protocol']['branch_taken']}`; the expanded grid was not
executed.

## Reconstruction of the exact relation

[COMPUTATION] Rebuilding the 23x21 training matrix from the raw `v^0..v^22`
coefficients with the e56 entry rule
`[v^n](v^p theta^q F) = (n-p)^q f_(n-p)` gives exact rank
{reconstruction['recomputed_training_rank']} and a one-dimensional kernel
spanned by the primitive integer vector stored in the frozen artifact
(span agreement certified; the vector annihilates all 23 training equations
exactly). In polynomial form the candidate relation is
`Q0(v) F + Q1(v) theta F + Q2(v) theta^2 F = 0` with

* `Q0 = [{q0[0]}, {q0[2]}, {q0[4]}, {q0[6]}]` on `v^0, v^2, v^4, v^6`,
* `Q1 = [{q1[0]}, {q1[2]}, {q1[4]}, {q1[6]}]`,
* `Q2 = [{q2[0]}, {q2[2]}, {q2[4]}, {q2[6]}]`.

## Decisive holdout evaluation

[LEMMA] (Parity.) Every column with a nonzero kernel weight has even `p`,
and every odd-index coefficient of the HT series is exactly zero, so the
residual of the candidate at any odd order is identically zero before any
arithmetic: orders 23, 25, 27 carry no information.

[COMPUTATION] The first informative holdout order is 24. The exact residual
of the frozen kernel there is

    [v^24] (Q0 + Q1 theta + Q2 theta^2) F = {holdout['first_nonzero_residual']}

which is a nonzero rational number: the candidate is refuted by the first
previously unobservable coefficient. Writing the residual as `S f_24 + T`,
where `S = sum_q c_(0,q) 24^q = {forced['S']}` is the exact weight the
relation puts on `f_24` and `T = {forced['T']}` is fixed by the training
prefix alone, the relation would force
`f_24 = {forced['forced_v24']}`, whereas the exact coefficient is
`f_24 = {forced['actual_v24']}`.

[COMPUTATION] The surviving-kernel dimension is 1 through order 23 and 0
from order 24 onward, so `first_refuting_holdout_order = 24`. The full
29-equation system has exact rank 21 (full column rank); the stored
nonzero 21x21 maximal minor on row orders
{holdout['certificate']['row_orders']} has exact determinant
`{holdout['certificate']['determinant']}`. Orders 26 and 28 give further
nonzero exact residuals ({holdout['confirmation_residuals']['26']} and
{holdout['confirmation_residuals']['28']}) but are flagged
confirmation-only: they played no role in the classification.

[COMPUTATION] An independent evaluator (the coordinating parent agent)
recomputed the three residuals from the artifacts before this script ran;
all three match the values above exactly (recorded in
`results/series/frontier3.json` under `independent_crosscheck`).

## Controls on the strict holdout machinery

[COMPUTATION] All {controls['e43_strict_holdout_row_count']} strict
`CANDIDATE_REFUTED_BY_HOLDOUT` rows of the frozen e43 first-order
differential-algebraic LT frontier were re-audited from the stored reduced
coefficients: the e43 U+2 budget ladder reproduces every stored first
refuting holdout order and every stored nonzero maximal-minor determinant
exactly. One frozen e43 HT `NO_RELATION_AT_BUDGET` row was likewise
reproduced with full training-prefix rank and its stored minor. These
controls demonstrate the exact three-way classification semantics on known
data.

## Verdict semantics and scope

* `NO_RELATION_AT_BUDGET`: a nonzero exact maximal minor of the training
  prefix proves no nonzero relation exists in the ansatz space at this
  budget.
* `CANDIDATE_REFUTED_BY_HOLDOUT`: a training kernel exists but later exact
  coefficients kill it; here first killed at withheld order `v^24`.
* `CANDIDATE_SURVIVES`: an observable exact kernel remains after every
  known coefficient; no row of this experiment carries this verdict.

[UNRESOLVED] This is a finite-prefix statement about one frozen ansatz row.
It does not prove non-D-finiteness or non-differential-algebraicity of the
true 3D Ising free energy, does not exclude relations outside the searched
spaces, and is not an exact solution of the model. No survival discovery
was performed and none is claimed.
"""


def main() -> None:
    checks: list[dict] = []

    frontier2_payload = json.loads(FRONTIER2_PATH.read_text(encoding="utf-8"))
    frontier2 = frontier2_payload["data"]
    ht28_payload = json.loads(HT28_PATH.read_text(encoding="utf-8"))
    ht24_payload = json.loads(HT24_PATH.read_text(encoding="utf-8"))
    ht26_payload = json.loads(HT26_PATH.read_text(encoding="utf-8"))
    ext2_payload = json.loads(EXT2_HT_PATH.read_text(encoding="utf-8"))
    e43_payload = json.loads(E43_PATH.read_text(encoding="utf-8"))
    e43 = e43_payload["data"]

    sources = [
        {
            "role": "frozen_candidate_artifact",
            "path": "results/series/frontier2.json",
            "sha256": _sha256(FRONTIER2_PATH),
        },
        {
            "role": "canonical_coefficients_v0_to_v28",
            "path": "results/series/ht_v28.json",
            "sha256": _sha256(HT28_PATH),
        },
        {
            "role": "v24_production_artifact",
            "path": "results/series/ht_v24.json",
            "sha256": _sha256(HT24_PATH),
        },
        {
            "role": "v26_production_artifact",
            "path": "results/series/ht_v26.json",
            "sha256": _sha256(HT26_PATH),
        },
        {
            "role": "e56_frozen_training_prefix",
            "path": "results/series/extended2_sc_ht_free_energy.json",
            "sha256": _sha256(EXT2_HT_PATH),
        },
        {
            "role": "e43_control_artifact",
            "path": "results/series/structure_certificates.json",
            "sha256": _sha256(E43_PATH),
        },
    ]

    # --- canonical coefficient chain ----------------------------------------
    coefficients = [
        Fraction(value) for value in ht28_payload["data"]["series"]["coefficients"]
    ]
    ext2_coefficients = [
        Fraction(value) for value in ext2_payload["data"]["coefficients"]
    ]
    frozen_input = [
        Fraction(value) for value in frontier2["input_series"]["HT"]["coefficients"]
    ]
    v24 = Fraction(ht24_payload["data"]["series"]["v24"])
    v26 = Fraction(ht26_payload["data"]["series"]["v26"])
    chain_ok = (
        len(coefficients) == CANONICAL_ORDER + 1
        and ht28_payload["data"]["series"]["achieved_order"] == CANONICAL_ORDER
        and ext2_coefficients == coefficients[: E56_FREEZE_ORDER + 1]
        and frozen_input == coefficients[: E56_FREEZE_ORDER + 1]
        and v24 == coefficients[DECISIVE_ORDER]
        and v26 == coefficients[26]
        and ht28_payload["data"]["series"]["v24_artifact_sha256"]
        == _sha256(HT24_PATH)
        and ht28_payload["data"]["series"]["v26_artifact_sha256"]
        == _sha256(HT26_PATH)
        and frontier2["input_series"]["HT"]["sha256"] == _sha256(EXT2_HT_PATH)
    )
    _record(
        checks,
        "canonical_source_chain",
        chain_ok,
        "ht_v28 prefix, extended2 prefix, frontier2 input, v24/v26 artifacts, and recorded hashes all agree exactly",
    )
    tu_rows = []
    for family, family_data in frontier2["frontiers"].items():
        for series_name in ("HT", "LT"):
            series_data = family_data.get(series_name)
            if not isinstance(series_data, dict):
                continue
            for row in series_data["frontier"]:
                if series_name == "HT" and row["verdict"] == "TRUNCATION_UNOBSERVABLE":
                    tu_rows.append((family, series_name, row))
    euler_tu = [
        (family, row)
        for family, series_name, row in tu_rows
        if row["ansatz_family"] == "linear_euler_ode"
    ]
    _record(
        checks,
        "sole_survivor_loaded",
        len(euler_tu) == 1
        and euler_tu[0][1]["parameters"] == {"order": 2, "polynomial_degree": 6},
        f"frontier2 holds {len(tu_rows)} HT truncation-unobservable row(s); the Euler one is the order-2 degree-6 candidate",
    )
    frozen = euler_tu[0][1]
    monomials = [tuple(item) for item in frozen["monomials"]]
    frozen_kernel = list(frozen["full_kernel_basis"][0])
    unknown_count = len(monomials)

    # --- predeclared protocol (fixed before any holdout evaluation) ------------
    predeclared_protocol = {
        "training_orders": list(range(E56_FREEZE_ORDER + 1)),
        "holdout_orders": list(range(E56_FREEZE_ORDER + 1, CANONICAL_ORDER + 1)),
        "refit_or_reselection_using_new_orders": False,
        "decision_tree": [
            {
                "branch": "REFUTED_AT_FIRST_UNOBSERVABLE_ORDER",
                "condition": "residual of the frozen kernel at order 24 is nonzero",
                "action": "classify CANDIDATE_REFUTED_BY_HOLDOUT with exact residual and nonzero maximal minor; record v26/v28 residuals as confirmation only",
            },
            {
                "branch": "ZERO_THROUGH_V24",
                "condition": "residual at order 24 is zero",
                "action": "evaluate order 26, then order 28, sequentially",
            },
            {
                "branch": "ALL_ZERO_THROUGH_V28",
                "condition": "residuals at orders 24, 26, 28 are all zero",
                "action": "recompute the predeclared expanded prefix-rank grid preserving holdouts; report only rank/unobservability accounting; no survival discovery; no D-finiteness claims",
            },
        ],
        "predeclared_expanded_prefix_rank_grid": {
            "cells": "(r, d) with 0 <= r <= 2 and (r+1)(d+1)+2 <= 29",
            "budget_rule": "U = (r+1)(d+1); training orders 0..U+1; holdouts U+2..28 (nonempty by construction)",
            "verdict_rule": "training-prefix rank, strict later-order holdouts, three-way semantics only",
            "status": "NOT_EXECUTED",
            "reason": None,
        },
    }

    # --- reconstruction from the training prefix alone -------------------------
    training_rows = [
        _euler_row(order, monomials, coefficients)
        for order in range(E56_FREEZE_ORDER + 1)
    ]
    training_rank = _rank(training_rows, unknown_count)
    basis = _nullspace(training_rows, unknown_count)
    rebuilt_primitive = _primitive(basis[0]) if basis else None
    frozen_residuals_zero = all(
        _euler_residual(order, monomials, frozen_kernel, coefficients) == 0
        for order in range(E56_FREEZE_ORDER + 1)
    )
    reconstruction_ok = (
        training_rank == unknown_count - 1
        and len(basis) == 1
        and rebuilt_primitive == frozen_kernel
        and frozen_residuals_zero
    )
    _record(
        checks,
        "kernel_reconstruction",
        reconstruction_ok,
        f"training rank {training_rank}, one-dimensional kernel; reconstructed primitive vector equals the frozen kernel and annihilates all 23 training orders",
    )

    polynomial_coefficients = {
        "Q0_v_powers_0_to_6": [
            frozen_kernel[0 + p] for p in range(7)
        ],
        "Q1_v_powers_0_to_6": [
            frozen_kernel[7 + p] for p in range(7)
        ],
        "Q2_v_powers_0_to_6": [
            frozen_kernel[14 + p] for p in range(7)
        ],
    }

    # --- parity support lemma ---------------------------------------------------
    effective_even = all(
        degree_t % 2 == 0
        for weight, (degree_t, _) in zip(frozen_kernel, monomials)
        if weight != 0
    )
    odd_support_zero = all(coefficients[odd] == 0 for odd in range(1, 29, 2))

    # --- holdout ladder ----------------------------------------------------------
    holdout_orders = list(range(E56_FREEZE_ORDER + 1, CANONICAL_ORDER + 1))
    holdout_rows = {
        str(order): _fraction_strings(
            _euler_row(order, monomials, coefficients)
        )
        for order in holdout_orders
    }
    residual_values = {
        order: _euler_residual(order, monomials, frozen_kernel, coefficients)
        for order in holdout_orders
    }
    surviving_ladder = []
    dimension = 1
    for order in holdout_orders:
        if dimension == 1 and residual_values[order] != 0:
            dimension = 0
        surviving_ladder.append({"order": order, "surviving_kernel_dimension": dimension})
    decisive_order = next(
        (
            entry["order"]
            for entry in surviving_ladder
            if entry["surviving_kernel_dimension"] == 0
        ),
        None,
    )

    roles = {}
    for order in holdout_orders:
        if order % 2 == 1:
            roles[order] = "parity_forced_zero"
        elif decisive_order is not None and order == decisive_order:
            roles[order] = "decisive"
        elif order in (26, 28) and decisive_order is not None and order > decisive_order:
            roles[order] = "confirmation_only"
        else:
            roles[order] = "evaluated"

    residual_records = [
        {
            "order": order,
            "residual": str(residual_values[order]),
            "nonzero": residual_values[order] != 0,
            "parity_forced_zero": order % 2 == 1 and effective_even and odd_support_zero,
            "role": roles[order],
        }
        for order in holdout_orders
    ]

    # forced-value decomposition at the decisive order
    big_s = sum(
        Fraction(weight) * Fraction(decisive_order) ** theta_power
        for weight, (degree_t, theta_power) in zip(frozen_kernel, monomials)
        if degree_t == 0
    )
    little_t = sum(
        Fraction(weight)
        * _euler_entry(decisive_order, degree_t, theta_power, coefficients)
        for weight, (degree_t, theta_power) in zip(frozen_kernel, monomials)
        if degree_t > 0
    )
    forced_value = -little_t / big_s
    forced_identity = (
        big_s * coefficients[decisive_order] + little_t
        == residual_values[decisive_order]
    )

    # full system rank and nonzero maximal minor
    full_rows = [
        _euler_row(order, monomials, coefficients)
        for order in range(CANONICAL_ORDER + 1)
    ]
    full_rank = _rank(full_rows, unknown_count)
    minor_orders = _independent_row_orders(full_rows, unknown_count)
    minor = [
        [full_rows[order][column] for column in range(unknown_count)]
        for order in minor_orders
    ]
    minor_determinant = _determinant(minor)

    branch_taken = (
        "REFUTED_AT_FIRST_UNOBSERVABLE_ORDER"
        if decisive_order == DECISIVE_ORDER
        else "ZERO_THROUGH_V24"
    )
    _record(
        checks,
        "decisive_residual_nonzero",
        decisive_order == DECISIVE_ORDER
        and residual_values[DECISIVE_ORDER] != 0
        and full_rank == unknown_count
        and minor_determinant != 0
        and forced_identity,
        f"[v^24] residual {residual_values[DECISIVE_ORDER]} != 0; full 29-equation rank {full_rank}; nonzero 21x21 minor stored",
    )
    _record(
        checks,
        "parity_support_lemma",
        effective_even
        and odd_support_zero
        and all(residual_values[odd] == 0 for odd in (23, 25, 27)),
        "even effective monomial support plus even series support forces odd-order residuals to vanish identically",
    )
    parent_match = all(
        Fraction(PARENT_REPORTED_RESIDUALS[order]) == residual_values[order]
        for order in (24, 26, 28)
    )
    _record(
        checks,
        "independent_crosscheck_agreement",
        parent_match,
        "parent-agent independent recomputation of the v^24/v^26/v^28 residuals matches exactly",
    )

    # --- e43 strict holdout controls ---------------------------------------------
    lt_reduced = [
        Fraction(value) for value in e43["input_series"]["LT"]["reduced_coefficients"]
    ]
    strict_rows = [
        row
        for row in e43["differential_algebraicity_order_1"]["LT"]["frontier"]
        if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ]
    control_records = []
    controls_ok = True
    for row in strict_rows:
        rows, mons = _da_matrix(
            lt_reduced, row["degree_x"], row["degree_f"], row["degree_f_prime"]
        )
        monomials_match = [list(item) for item in mons] == row["monomials"]
        ladder = _budget_ladder(rows, len(mons))
        subminor = [
            [
                rows[order][column]
                for column in row["certificate"]["column_indices"]
            ]
            for order in row["certificate"]["row_orders"]
        ]
        determinant = _determinant(subminor)
        ok = (
            monomials_match
            and ladder["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
            and ladder["first_refuting_holdout_order"]
            == row["first_refuting_holdout_order"]
            and determinant == Fraction(row["certificate"]["determinant"])
        )
        controls_ok = controls_ok and ok
        control_records.append(
            {
                "series": "LT",
                "ansatz_family": "first_order_differential_algebraic",
                "parameters": {
                    "degree_x": row["degree_x"],
                    "degree_f": row["degree_f"],
                    "degree_f_prime": row["degree_f_prime"],
                },
                "unknown_count": row["unknown_count"],
                "stored_first_refuting_holdout_order": row["first_refuting_holdout_order"],
                "recomputed_first_refuting_holdout_order": ladder[
                    "first_refuting_holdout_order"
                ],
                "stored_determinant": row["certificate"]["determinant"],
                "recomputed_determinant": str(determinant),
                "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
                "monomial_enumeration_matches": monomials_match,
                "match": ok,
            }
        )

    ht_reduced = [
        Fraction(value) for value in e43["input_series"]["HT"]["reduced_coefficients"]
    ]
    no_relation_row = e43["differential_algebraicity_order_1"]["HT"]["frontier"][0]
    nr_rows, nr_mons = _da_matrix(
        ht_reduced,
        no_relation_row["degree_x"],
        no_relation_row["degree_f"],
        no_relation_row["degree_f_prime"],
    )
    nr_ladder = _budget_ladder(nr_rows, len(nr_mons))
    nr_subminor = [
        [
            nr_rows[order][column]
            for column in no_relation_row["certificate"]["column_indices"]
        ]
        for order in no_relation_row["certificate"]["row_orders"]
    ]
    nr_determinant = _determinant(nr_subminor)
    no_relation_ok = (
        nr_ladder["verdict"] == "NO_RELATION_AT_BUDGET"
        and nr_ladder["training_rank"] == len(nr_mons)
        and nr_determinant == Fraction(no_relation_row["certificate"]["determinant"])
    )
    controls_ok = controls_ok and no_relation_ok
    _record(
        checks,
        "e43_controls_reproduced",
        controls_ok,
        f"{len(strict_rows)} strict LT holdout refutations and 1 HT NO_RELATION row replay exactly (first refuting orders and minor determinants)",
    )

    # --- assemble payload -----------------------------------------------------------
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "interpreter": INTERPRETER,
            "python_version": sys.version.split()[0],
            "method": "exact holdout evaluation of one frozen e56 Euler-ODE candidate against the reserved v^24/v^26/v^28 coefficients",
            "benchmark_used_for_fit_selection_or_validation": False,
            "sources": sources,
        },
        "data": {
            "claim_tag": "[COMPUTATION]",
            "input_series": {
                "HT": {
                    "artifact": "results/series/ht_v28.json",
                    "sha256": _sha256(HT28_PATH),
                    "variable": "v",
                    "achieved_order": CANONICAL_ORDER,
                    "coefficient_count": len(coefficients),
                    "normalization": ht28_payload["data"]["series"]["normalization"],
                    "coefficients": _fraction_strings(coefficients),
                    "v24": str(coefficients[24]),
                    "v26": str(coefficients[26]),
                    "v24_artifact": "results/series/ht_v24.json",
                    "v24_artifact_sha256": _sha256(HT24_PATH),
                    "v26_artifact": "results/series/ht_v26.json",
                    "v26_artifact_sha256": _sha256(HT26_PATH),
                    "e56_training_prefix_artifact": "results/series/extended2_sc_ht_free_energy.json",
                    "e56_training_prefix_sha256": _sha256(EXT2_HT_PATH),
                }
            },
            "candidate_provenance": {
                "frozen_artifact": "results/series/frontier2.json",
                "frozen_artifact_sha256": _sha256(FRONTIER2_PATH),
                "frozen_certificate_pointer": frozen["certificate_pointer"],
                "ansatz_family": frozen["ansatz_family"],
                "ansatz_definition": "sum_{q<=2} Q_q(v) theta^q F = 0 with theta = v d/dv, deg Q_q <= 6, 21 monomial columns v^p theta^q F",
                "parameters": frozen["parameters"],
                "e56_frozen_verdict": frozen["verdict"],
                "sole_survivor_audit": {
                    "ht_truncation_unobservable_rows_total": len(tu_rows),
                    "euler_ht_truncation_unobservable_rows": len(euler_tu),
                    "families_scanned": sorted(frontier2["frontiers"].keys()),
                },
                "monomials": [list(item) for item in monomials],
                "column_labels": [
                    f"v^{degree_t}*theta^{theta_power}*F"
                    for degree_t, theta_power in monomials
                ],
                "frozen_kernel_basis": [frozen_kernel],
                "polynomial_coefficients": polynomial_coefficients,
                "frozen_support": {
                    "support_components": frozen["support_components"],
                    "first_unfixed_supported_order": frozen["support_components"][0][
                        "first_unfixed_supported_order"
                    ],
                },
                "relation_reconstruction": {
                    "matrix_entry_rule": "[v^n](v^p theta^q F) = (n-p)^q f_(n-p) if n >= p else 0",
                    "training_row_orders": list(range(E56_FREEZE_ORDER + 1)),
                    "recomputed_training_rank": training_rank,
                    "recomputed_training_nullity": unknown_count - training_rank,
                    "recomputed_kernel_primitive": rebuilt_primitive,
                    "span_agrees_with_frozen_kernel": rebuilt_primitive == frozen_kernel,
                    "frozen_kernel_annihilates_training": frozen_residuals_zero,
                    "new_orders_used_in_reconstruction": False,
                },
            },
            "holdout_protocol": {
                **predeclared_protocol,
                "branch_taken": branch_taken,
            },
            "holdout_evaluation": {
                "matrix": {
                    "row_orders": list(range(CANONICAL_ORDER + 1)),
                    "column_labels": [
                        f"v^{degree_t}*theta^{theta_power}*F"
                        for degree_t, theta_power in monomials
                    ],
                    "rows": {
                        str(order): _fraction_strings(
                            _euler_row(order, monomials, coefficients)
                        )
                        for order in range(CANONICAL_ORDER + 1)
                    },
                },
                "residuals": residual_records,
                "surviving_dimension_ladder": surviving_ladder,
                "decisive_order": decisive_order,
                "first_nonzero_residual": str(residual_values[decisive_order]),
                "forced_value_decomposition": {
                    "S": str(big_s),
                    "T": str(little_t),
                    "forced_v24": str(forced_value),
                    "actual_v24": str(coefficients[decisive_order]),
                    "identity_S_times_actual_plus_T_equals_residual": forced_identity,
                },
                "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
                "first_refuting_holdout_order": decisive_order,
                "full_system": {
                    "row_orders": list(range(CANONICAL_ORDER + 1)),
                    "rank": full_rank,
                    "nullity": unknown_count - full_rank,
                },
                "certificate": {
                    "type": "NONZERO_MAXIMAL_MINOR",
                    "phase": "all_available_equations_orders_0_to_28",
                    "row_orders": minor_orders,
                    "column_indices": list(range(unknown_count)),
                    "determinant": str(minor_determinant),
                },
                "confirmation_only_orders": [26, 28],
                "confirmation_residuals": {
                    "26": str(residual_values[26]),
                    "28": str(residual_values[28]),
                },
                "independent_crosscheck": {
                    "reported_by": "coordinating parent agent, computed independently from the artifacts before this script ran",
                    "residuals": {key: value for key, value in PARENT_REPORTED_RESIDUALS.items()},
                    "all_match": parent_match,
                    "used_in_computation": False,
                },
            },
            "controls": {
                "description": "re-audit of the frozen e43 strict holdout refutations (and one NO_RELATION row) as classification-semantics controls",
                "e43_strict_holdout_row_count": len(strict_rows),
                "e43_strict_holdout_rows": control_records,
                "no_relation_control": {
                    "series": "HT",
                    "ansatz_family": "first_order_differential_algebraic",
                    "parameters": {
                        "degree_x": no_relation_row["degree_x"],
                        "degree_f": no_relation_row["degree_f"],
                        "degree_f_prime": no_relation_row["degree_f_prime"],
                    },
                    "recomputed_training_rank": nr_ladder["training_rank"],
                    "unknown_count": len(nr_mons),
                    "recomputed_determinant": str(nr_determinant),
                    "verdict": "NO_RELATION_AT_BUDGET",
                    "match": no_relation_ok,
                },
                "all_controls_match": controls_ok,
            },
            "verdict_semantics": {
                "NO_RELATION_AT_BUDGET": "a nonzero exact maximal minor of the training prefix proves no nonzero relation in the ansatz space at this budget",
                "CANDIDATE_REFUTED_BY_HOLDOUT": "a training kernel exists but later exact coefficients kill it (here first at withheld order v^24)",
                "CANDIDATE_SURVIVES": "an observable exact kernel remains after every known coefficient; no row of this experiment carries this verdict",
                "finite_budget_caveat": "all statements are finite-prefix; none proves non-D-finiteness or non-differential-algebraicity",
            },
            "scope_warning": "finite-budget exact computation on one frozen ansatz row; no survival discovery performed; no D-finiteness or exact-solution claim; benchmarks never used",
        },
        "checks": checks,
    }

    RESULT_PATH.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    PROOF_PATH.write_text(_proof_text(payload), encoding="utf-8")

    vocabulary_ok = True

    def _walk(node):
        nonlocal vocabulary_ok
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "verdict" and value not in ALLOWED_VERDICTS:
                    vocabulary_ok = False
                _walk(value)
        elif isinstance(node, list):
            for item in node:
                _walk(item)

    _walk(payload["data"])
    _record(
        checks,
        "verdict_vocabulary",
        vocabulary_ok,
        "every stored verdict label stays inside the three-way semantics and no CANDIDATE_SURVIVES is emitted",
    )
    RESULT_PATH.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    if not all(item["passed"] for item in checks):
        raise AssertionError("a recorded check failed")
    print("PASS e113_series_frontier3")


if __name__ == "__main__":
    main()
