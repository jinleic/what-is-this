"""Exact spectral non-Gaussianity certificates for the open 2x5 Ising layer.

The transfer operator and every symmetry decision are exact over the integers.
Float64 is used only to propose an integer congruence basis and rational
spectral windows.  For every stored shift, an exact Frobenius perturbation
inequality certifies the inertia of every symmetry-sector pencil.
"""

from __future__ import annotations

import gc
import json
import math
import os
import platform
import resource
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.e54_spectral_3x3 import (  # noqa: E402
    Projector,
    ScaledOperator,
    Sector,
    build_scaled_operator,
    build_sectors,
    chain_actions,
    chain_projectors,
    interval_json,
    projector_columns,
    state_permutation,
    verify_projectors,
)
from ising.transfer_matrix import layer_bonds  # noqa: E402


DYADIC_BITS = 36
ENCLOSURE_RADIUS = Fraction(1, 1000)


@dataclass
class RectangleAction:
    name: str
    state_map: list[int]
    row_reflected: bool
    column_reflected: bool
    flipped: bool


@dataclass
class PreparedSector:
    label: str
    dimension: int
    operator_congruence: np.ndarray
    gram_congruence: np.ndarray
    guide_eigenvalues: list[float]
    preparation_seconds: float
    operator_congruence_symmetric: bool
    gram_congruence_symmetric: bool


@dataclass
class ExactCount:
    total: int
    by_sector: dict[str, int]
    sector_witnesses: dict[str, dict[str, object]]
    seconds: float


def ceiling_sqrt(value: int) -> int:
    root = math.isqrt(value)
    return root if root * root == value else root + 1


def rectangle_actions(rows: int = 2, columns: int = 5) -> list[RectangleAction]:
    """Return the faithful C2(row) x C2(column) x C2(spin-flip) action."""
    n = rows * columns
    full_mask = (1 << n) - 1
    actions: list[RectangleAction] = []
    for row_reflected in (False, True):
        for column_reflected in (False, True):
            sites: list[int] = []
            for row in range(rows):
                for column in range(columns):
                    target_row = rows - 1 - row if row_reflected else row
                    target_column = columns - 1 - column if column_reflected else column
                    sites.append(target_row * columns + target_column)
            spatial = state_permutation(n, sites)
            for flipped in (False, True):
                mapping = [value ^ full_mask if flipped else value for value in spatial]
                actions.append(
                    RectangleAction(
                        (
                            f"row{'-' if row_reflected else '+'}_"
                            f"column{'-' if column_reflected else '+'}_"
                            f"flip{'-' if flipped else '+'}"
                        ),
                        mapping,
                        row_reflected,
                        column_reflected,
                        flipped,
                    )
                )
    if len(actions) != 8 or len({tuple(action.state_map) for action in actions}) != 8:
        raise AssertionError("the rectangle symmetry action is not faithful")
    return actions


def rectangle_projectors(
    dimension: int, actions: list[RectangleAction]
) -> list[Projector]:
    """Build all eight exact one-dimensional character projectors."""
    projectors: list[Projector] = []
    for row_sign in (1, -1):
        for column_sign in (1, -1):
            for flip_sign in (1, -1):
                terms: list[tuple[int, list[int]]] = []
                characters: list[int] = []
                for action in actions:
                    character = row_sign if action.row_reflected else 1
                    if action.column_reflected:
                        character *= column_sign
                    if action.flipped:
                        character *= flip_sign
                    characters.append(character)
                    terms.append((character, action.state_map))
                label = (
                    f"row{'+' if row_sign == 1 else '-'}_"
                    f"column{'+' if column_sign == 1 else '-'}_"
                    f"flip{'+' if flip_sign == 1 else '-'}"
                )
                projectors.append(
                    Projector(
                        label,
                        8,
                        terms,
                        projector_columns(dimension, terms),
                        1,
                        flip_sign,
                        characters,
                    )
                )
    return projectors


def exact_projector_verification(
    matrix: list[list[int]], actions: list[object], projectors: list[Projector]
) -> dict[str, object]:
    verification = verify_projectors(matrix, actions, projectors)
    required = (
        "common_denominator_verified",
        "projector_self_adjoint_exact",
        "projector_idempotence_exact",
        "projector_orthogonality_exact",
        "projector_completeness_exact",
        "projector_traces_integral",
        "matrix_invariant_under_every_group_action_exact",
        "projectors_commute_with_matrix_exact",
    )
    if not all(bool(verification[key]) for key in required):
        raise AssertionError("an exact projector identity or commutation check failed")
    return verification


def prepare_sector_witness(sector: Sector, dyadic_bits: int) -> PreparedSector:
    """Prepare an integer congruence basis; float64 proposes, integers certify."""
    started = time.perf_counter()
    size = sector.dimension
    if size == 0:
        raise AssertionError("empty symmetry sector")
    gram_diagonal = [sector.gram[index][index] for index in range(size)]
    if any(value <= 0 for value in gram_diagonal):
        raise AssertionError("sector Gram matrix is not positive diagonal")
    if any(
        sector.gram[row][column] != 0
        for row in range(size)
        for column in range(size)
        if row != column
    ):
        raise AssertionError("the abelian orbit basis should have diagonal Gram matrix")

    largest = max(abs(value) for row in sector.matrix for value in row)
    if largest == 0:
        raise AssertionError("zero restricted operator")
    inverse_roots = np.array([1.0 / math.sqrt(value) for value in gram_diagonal])
    proposal = np.empty((size, size), dtype=float)
    for row in range(size):
        for column in range(size):
            proposal[row, column] = (
                (sector.matrix[row][column] / largest)
                * inverse_roots[row]
                * inverse_roots[column]
            )
    proposal = 0.5 * (proposal + proposal.T)
    eigenvalues_scaled, eigenvectors = np.linalg.eigh(proposal)
    guides = [float(value * largest) for value in eigenvalues_scaled]

    # Columns approximate generalized eigenvectors in the original sector basis.
    generalized_vectors = inverse_roots[:, None] * eigenvectors
    dyadic_denominator = 1 << dyadic_bits
    q_integer_i64 = np.rint(generalized_vectors * dyadic_denominator).astype(np.int64)
    q_integer = q_integer_i64.astype(object)

    operator = np.array(sector.matrix, dtype=object)
    operator_congruence = q_integer.T @ (operator @ q_integer)
    weighted_q = q_integer.copy()
    for row, gram_value in enumerate(gram_diagonal):
        weighted_q[row, :] *= gram_value
    gram_congruence = q_integer.T @ weighted_q

    operator_symmetric = all(
        int(operator_congruence[row, column])
        == int(operator_congruence[column, row])
        for row in range(size)
        for column in range(row)
    )
    gram_symmetric = all(
        int(gram_congruence[row, column]) == int(gram_congruence[column, row])
        for row in range(size)
        for column in range(row)
    )
    if not operator_symmetric or not gram_symmetric:
        raise AssertionError("an exact congruence product lost symmetry")
    return PreparedSector(
        sector.label,
        size,
        operator_congruence,
        gram_congruence,
        guides,
        time.perf_counter() - started,
        operator_symmetric,
        gram_symmetric,
    )


class ExactSectorCounter:
    """Exact sector-summed inertia from strict integer congruence inequalities."""

    def __init__(self, sectors: list[Sector], dyadic_bits: int = DYADIC_BITS):
        started = time.perf_counter()
        self.dyadic_bits = dyadic_bits
        self.prepared = [prepare_sector_witness(sector, dyadic_bits) for sector in sectors]
        self.preparation_seconds = time.perf_counter() - started
        self.cache: dict[Fraction, ExactCount] = {}
        self.count_seconds = 0.0
        self.max_count_seconds = 0.0

    def guides(self) -> list[float]:
        return sorted(
            value
            for prepared in self.prepared
            for value in prepared.guide_eigenvalues
        )[:3]

    def count(self, shift: Fraction) -> ExactCount:
        shift = Fraction(shift)
        if shift in self.cache:
            return self.cache[shift]
        started = time.perf_counter()
        by_sector: dict[str, int] = {}
        witnesses: dict[str, dict[str, object]] = {}
        total = 0
        for prepared in self.prepared:
            size = prepared.dimension
            diagonal: list[int] = []
            off_diagonal_entries: list[tuple[int, int, int]] = []
            off_diagonal_square_sum = 0
            for row in range(size):
                for column in range(row + 1):
                    entry = (
                        shift.denominator
                        * int(prepared.operator_congruence[row, column])
                        - shift.numerator
                        * int(prepared.gram_congruence[row, column])
                    )
                    if row == column:
                        diagonal.append(entry)
                    else:
                        off_diagonal_entries.append((row, column, entry))
                        off_diagonal_square_sum += 2 * entry * entry
            if any(value == 0 for value in diagonal):
                raise ArithmeticError(
                    f"zero proposed diagonal for {prepared.label} at shift {shift}"
                )

            # Scale each coordinate by sqrt(abs(D_ii)).  If the Frobenius norm
            # of the resulting relative off-diagonal matrix is < 1, then
            # J+tF is nonsingular for 0<=t<=1, so C and D have equal inertia.
            # Dyadic lower bounds for abs(D_ii) keep the comparison integral.
            diagonal_exponents = [abs(value).bit_length() - 1 for value in diagonal]
            relative_common_exponent = max(
                (
                    diagonal_exponents[row] + diagonal_exponents[column]
                    for row, column, _ in off_diagonal_entries
                ),
                default=0,
            )
            relative_square_numerator = sum(
                2
                * entry
                * entry
                * (
                    1
                    << (
                        relative_common_exponent
                        - diagonal_exponents[row]
                        - diagonal_exponents[column]
                    )
                )
                for row, column, entry in off_diagonal_entries
            )
            relative_square_denominator = 1 << relative_common_exponent
            relative_safety_margin = (
                relative_square_denominator - relative_square_numerator
            )
            certified = relative_safety_margin > 0
            off_diagonal_bound = ceiling_sqrt(off_diagonal_square_sum)
            minimum_diagonal = min(abs(value) for value in diagonal)
            global_safety_margin = minimum_diagonal - off_diagonal_bound
            if not certified:
                raise ArithmeticError(
                    f"integer relative-congruence witness failed for {prepared.label} "
                    f"at shift {shift}: relative margin {relative_safety_margin}, "
                    f"global margin {global_safety_margin}"
                )
            inertia = sum(value < 0 for value in diagonal)
            by_sector[prepared.label] = inertia
            total += inertia
            witnesses[prepared.label] = {
                "claim_tag": "[COMPUTATION]",
                "dimension": size,
                "inertia_count_below_shift": inertia,
                "minimum_abs_diagonal_numerator": str(minimum_diagonal),
                "offdiagonal_frobenius_upper_numerator": str(off_diagonal_bound),
                "global_frobenius_safety_margin_numerator": str(global_safety_margin),
                "relative_frobenius_squared_upper_numerator": str(
                    relative_square_numerator
                ),
                "relative_frobenius_squared_upper_denominator": str(
                    relative_square_denominator
                ),
                "strict_safety_margin_numerator": str(relative_safety_margin),
                "strict_relative_frobenius_inequality": certified,
                "congruence_basis_invertible_exact": certified,
                "inertia_certified_exactly": certified,
            }
        elapsed = time.perf_counter() - started
        record = ExactCount(total, by_sector, witnesses, elapsed)
        self.cache[shift] = record
        self.count_seconds += elapsed
        self.max_count_seconds = max(self.max_count_seconds, elapsed)
        return record

    def resource_record(self) -> dict[str, object]:
        return {
            "claim_tag": "[COMPUTATION]",
            "method": (
                "For each exact sector pencil B-sigma G, float64 proposes an integer "
                "matrix Q. Exact integer products form C=Q^T(dB-nG)Q=D+E. After "
                "scaling coordinate i by sqrt(abs(D_ii)), exact dyadic lower bounds "
                "on abs(D_ii) prove the relative off-diagonal Frobenius norm is <1. "
                "Thus diag(sign(D_ii))+tF is nonsingular for 0<=t<=1, proving Q "
                "invertible and fixing inertia(C)=inertia(D) by Sylvester congruence. "
                "Counts are summed over every sector in the complete projector family."
            ),
            "dyadic_proposal_bits": self.dyadic_bits,
            "dyadic_proposal_denominator": str(1 << self.dyadic_bits),
            "sector_preparation": {
                prepared.label: {
                    "dimension": prepared.dimension,
                    "seconds": round(prepared.preparation_seconds, 6),
                    "operator_congruence_symmetric_exact": prepared.operator_congruence_symmetric,
                    "gram_congruence_symmetric_exact": prepared.gram_congruence_symmetric,
                }
                for prepared in self.prepared
            },
            "preparation_seconds": round(self.preparation_seconds, 6),
            "count_seconds": round(self.count_seconds, 6),
            "max_single_count_seconds": round(self.max_count_seconds, 6),
            "exact_shifts_certified": len(self.cache),
            "shift_witnesses": {
                str(shift): {
                    "shift_A_units": str(shift),
                    "shift_numerator": str(shift.numerator),
                    "shift_denominator": str(shift.denominator),
                    "total_inertia_count_below_shift": record.total,
                    "by_sector": record.by_sector,
                    "sectors": record.sector_witnesses,
                }
                for shift, record in self.cache.items()
            },
        }


def endpoint_record(record: ExactCount) -> dict[str, object]:
    return {"total": record.total, "by_sector": record.by_sector}


def certify_enclosure(
    counter: ExactSectorCounter, index: int, guide: float
) -> tuple[Fraction, Fraction, ExactCount, ExactCount]:
    if not math.isfinite(guide) or guide <= 0:
        raise ValueError("a positive finite guide is required")
    center = Fraction(format(guide, ".16e"))
    radius = ENCLOSURE_RADIUS
    lower = center * (1 - radius)
    upper = center * (1 + radius)
    lower_record = counter.count(lower)
    upper_record = counter.count(upper)
    if lower_record.total != index or upper_record.total != index + 1:
        raise AssertionError(
            f"the fixed rational bracket did not isolate eigenvalue {index}: "
            f"counts {lower_record.total}, {upper_record.total}"
        )
    return lower, upper, lower_record, upper_record


def case_actions_and_projectors(
    kind: str, dimension: int
) -> tuple[list[object], list[Projector], str]:
    if kind == "3D-layer":
        actions = rectangle_actions()
        projectors = rectangle_projectors(dimension, actions)
        symmetry = "C2(row reflection) x C2(column reflection) x C2(global spin flip)"
    elif kind == "1D-control":
        actions = chain_actions(10)
        projectors = chain_projectors(10, actions)
        symmetry = "C2(chain reflection) x C2(global spin flip)"
    else:
        raise ValueError(kind)
    return list(actions), projectors, symmetry


def certify_case(
    name: str, kind: str, bonds: list[tuple[int, int]], t: Fraction
) -> dict[str, object]:
    started = time.perf_counter()
    timings: dict[str, float] = {}

    began = time.perf_counter()
    scaled: ScaledOperator = build_scaled_operator(10, bonds, t)
    timings["matrix_build_seconds"] = time.perf_counter() - began

    actions, projectors, symmetry_name = case_actions_and_projectors(kind, 1024)
    began = time.perf_counter()
    verification = exact_projector_verification(scaled.matrix, actions, projectors)
    timings["projector_verification_seconds"] = time.perf_counter() - began

    began = time.perf_counter()
    sectors = build_sectors(scaled.matrix, actions, projectors)
    timings["sector_build_seconds"] = time.perf_counter() - began
    if sum(sector.dimension for sector in sectors) != 1024:
        raise AssertionError("sector dimensions are not complete")

    began = time.perf_counter()
    counter = ExactSectorCounter(sectors)
    timings["congruence_preparation_seconds"] = time.perf_counter() - began

    positive = counter.count(Fraction(0))
    if positive.total != 0:
        raise AssertionError("the exact transfer operator is not positive definite")

    enclosures: list[tuple[Fraction, Fraction]] = []
    endpoint_records: list[tuple[ExactCount, ExactCount]] = []
    for index, guide in enumerate(counter.guides()):
        lower, upper, lower_record, upper_record = certify_enclosure(counter, index, guide)
        enclosures.append((lower, upper))
        endpoint_records.append((lower_record, upper_record))

    (l0, h0), (l1, h1), (l2, h2) = enclosures
    distinct = h0 < l1 and h1 < l2
    if not distinct:
        raise AssertionError("the three lowest certified enclosures overlap")
    predicted_lower = l1 * l2 / h0
    predicted_upper = h1 * h2 / l0
    forced_lower_record = counter.count(predicted_lower)
    forced_upper_record = counter.count(predicted_upper)
    absent = forced_lower_record.total == forced_upper_record.total

    lambda_fields: dict[str, object] = {}
    for index, ((lower, upper), (lower_record, upper_record)) in enumerate(
        zip(enclosures, endpoint_records)
    ):
        lower_by_sector = lower_record.by_sector
        upper_by_sector = upper_record.by_sector
        lambda_fields[f"lambda{index}"] = {
            **interval_json(lower / scaled.scale, upper / scaled.scale),
            "endpoint_inertia": {
                "lower": endpoint_record(lower_record),
                "upper": endpoint_record(upper_record),
            },
            "sector_multiplicity_in_enclosure": {
                label: upper_by_sector[label] - lower_by_sector[label]
                for label in lower_by_sector
                if upper_by_sector[label] != lower_by_sector[label]
            },
        }

    timings["exact_count_seconds"] = counter.count_seconds
    timings["total_seconds"] = time.perf_counter() - started
    is_layer = kind == "3D-layer"
    statement = (
        f"The exact open 2x5 layer transfer operator at t={t} is not a "
        "positive-spectrum ten-mode Gaussian subset-product operator: its three "
        "lowest eigenvalues are distinct and the Gaussian-forced interval is empty."
        if is_layer
        else f"For the exact open n=10 chain at t={t}, the certified Gaussian-forced "
        "window is occupied. This is a consistency control, not an exact equality proof."
    )
    row = {
        "name": name,
        "kind": kind,
        "claim_tag": "[THEOREM]" if is_layer else "[COMPUTATION]",
        "statement": statement,
        "status": "CERTIFIED",
        "t_tanh_half_Kstar": str(t),
        "exp_2K": str(scaled.q),
        "n_sites": 10,
        "n_bonds": len(bonds),
        "dimension": 1024,
        "integer_scale_for_R": str(scaled.scale),
        "bond_parity": scaled.bond_parity,
        "energy_exponent_range": list(scaled.energy_exponents),
        "positive_definite_count_below_zero": positive.total,
        "lambda_enclosures": lambda_fields,
        "three_lowest_distinct": distinct,
        "relative_gap_lower_bounds": {
            "gap01": str((l1 - h0) / l1),
            "gap12": str((l2 - h1) / l2),
        },
        "forced_value_window": {
            "predicted": interval_json(
                predicted_lower / scaled.scale, predicted_upper / scaled.scale
            ),
            "actually_tested": interval_json(
                predicted_lower / scaled.scale, predicted_upper / scaled.scale
            ),
            "tested_window_encloses_prediction": True,
            "endpoint_inertia": {
                "lower": endpoint_record(forced_lower_record),
                "upper": endpoint_record(forced_upper_record),
            },
        },
        "inertia_counts_at_forced_window_ends": [
            forced_lower_record.total,
            forced_upper_record.total,
        ],
        "predicted_eigenvalue_absent": absent,
        "symmetry": {
            "claim_tag": "[COMPUTATION]",
            "group": symmetry_name,
            "verification": verification,
            "sector_dimensions": {
                sector.label: sector.dimension for sector in sectors
            },
            "sector_dimensions_sum": sum(sector.dimension for sector in sectors),
            "sector_irrep_dimensions": {
                sector.label: sector.irrep_dimension for sector in sectors
            },
            "basis_convention": (
                "Exact independent integer columns of each character-projector numerator; "
                "B=C^T A C and G=C^T C. Every stored inertia is the exact sum over all sectors."
            ),
        },
        "exact_sector_inertia_certificate": counter.resource_record(),
        "exact_inertia_shifts": len(counter.cache),
        "elapsed": {key: round(value, 6) for key, value in timings.items()},
    }
    return row


def check_record(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    grid_bonds = list(layer_bonds((2, 5), (False, False)))
    chain_bonds = list(layer_bonds((10,), (False,)))
    definitions = [
        ("layer_2x5_t_1_3", "3D-layer", grid_bonds, Fraction(1, 3)),
        ("layer_2x5_t_1_2", "3D-layer", grid_bonds, Fraction(1, 2)),
        ("control_chain_n10_t_1_3", "1D-control", chain_bonds, Fraction(1, 3)),
        ("control_chain_n10_t_1_2", "1D-control", chain_bonds, Fraction(1, 2)),
    ]
    rows: list[dict[str, object]] = []
    all_started = time.perf_counter()
    for name, kind, bonds, t in definitions:
        print(f"certifying {name}", flush=True)
        row = certify_case(name, kind, bonds, t)
        rows.append(row)
        print(
            f"  counts={row['inertia_counts_at_forced_window_ends']} "
            f"absent={row['predicted_eigenvalue_absent']} "
            f"seconds={row['elapsed']['total_seconds']}",
            flush=True,
        )
        gc.collect()

    layer_rows = [row for row in rows if row["kind"] == "3D-layer"]
    control_rows = [row for row in rows if row["kind"] == "1D-control"]
    checks = [
        check_record(
            "two_exact_open_2x5_absence_certificates",
            len(layer_rows) == 2
            and all(
                row["three_lowest_distinct"]
                and row["predicted_eigenvalue_absent"]
                and len(set(row["inertia_counts_at_forced_window_ends"])) == 1
                for row in layer_rows
            ),
            "; ".join(
                f"{row['name']} counts={row['inertia_counts_at_forced_window_ends']}"
                for row in layer_rows
            ),
        ),
        check_record(
            "n10_chain_forced_windows_occupied",
            len(control_rows) == 2
            and all(
                row["three_lowest_distinct"]
                and not row["predicted_eigenvalue_absent"]
                and row["inertia_counts_at_forced_window_ends"][1]
                > row["inertia_counts_at_forced_window_ends"][0]
                for row in control_rows
            ),
            "; ".join(
                f"{row['name']} counts={row['inertia_counts_at_forced_window_ends']}"
                for row in control_rows
            ),
        ),
        check_record(
            "all_projector_identities_completeness_and_commutation_exact",
            all(
                row["symmetry"]["verification"][key]
                for row in rows
                for key in (
                    "projector_self_adjoint_exact",
                    "projector_idempotence_exact",
                    "projector_orthogonality_exact",
                    "projector_completeness_exact",
                    "projectors_commute_with_matrix_exact",
                )
            ),
            "every projector flag passed in all four independently constructed operators",
        ),
        check_record(
            "all_sector_sums_are_complete_1024",
            all(row["symmetry"]["sector_dimensions_sum"] == 1024 for row in rows),
            "; ".join(
                f"{row['name']} sum={row['symmetry']['sector_dimensions_sum']}"
                for row in rows
            ),
        ),
        check_record(
            "every_stored_sector_inertia_has_strict_exact_witness",
            all(
                sector["inertia_certified_exactly"]
                and int(sector["strict_safety_margin_numerator"]) > 0
                for row in rows
                for shift in row["exact_sector_inertia_certificate"]["shift_witnesses"].values()
                for sector in shift["sectors"].values()
            ),
            "all integer congruence Frobenius margins are strictly positive",
        ),
        check_record(
            "open_2x5_graph_and_untuned_couplings",
            len(grid_bonds) == 13
            and max(sum(vertex in bond for bond in grid_bonds) for vertex in range(10)) == 3
            and {row["t_tanh_half_Kstar"] for row in layer_rows} == {"1/3", "1/2"},
            "open 2x5 has 13 bonds and degree-3 vertices; t=1/3 and t=1/2 were predeclared simple rationals",
        ),
    ]

    total_seconds = time.perf_counter() - all_started
    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {
        "provenance": {
            "script": "experiments/e67_spectral_2x5.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "method": (
                "Exact integer transfer operators using e39/e54 conventions; exact character "
                "projectors; exact sector restrictions; exact integer congruence/Frobenius "
                "inertia witnesses. Float64 only proposes bases and brackets and decides no claim."
            ),
            "coupling_selection": (
                "t=1/3 and t=1/2 were predeclared simple rational carry-over/stress points. "
                "No critical-coupling benchmark was used to tune, select, or fit them."
            ),
        },
        "data": {
            "claim_tag": "[THEOREM]",
            "scope": (
                "Exactly the finite open 2x5 layer at t in {1/3,1/2}, dimension 1024; "
                "a whole-spectrum positive-spectrum ten-mode Gaussian subset-product obstruction only."
            ),
            "parameterization": {
                "claim_tag": "[LEMMA]",
                "t": "tanh(K*/2)",
                "exp_2K": "(1+t^2)/(2t)",
                "operator": "R=P_t diag(q^((b-eps)/2)) P_t, P_t[k,l]=t^hamming(k,l)",
                "forced_value": "lambda1*lambda2/lambda0 must be spectral when lambda0<lambda1<lambda2 in a positive Gaussian subset-product spectrum",
            },
            "rows": rows,
            "resource_record": {
                "claim_tag": "[COMPUTATION]",
                "total_wall_seconds": round(total_seconds, 6),
                "process_max_rss_raw": max_rss,
                "process_max_rss_unit": "bytes on Darwin; platform-defined getrusage units otherwise",
                "cpu": platform.processor() or platform.machine(),
                "dyadic_proposal_bits": DYADIC_BITS,
                "case_total_seconds": {
                    row["name"]: row["elapsed"]["total_seconds"] for row in rows
                },
                "dense_1024_check_performed": False,
                "dense_1024_check_status": "[UNRESOLVED] optional; complete exact sector sums decide every stored inertia",
            },
            "limitations": {
                "claim_tag": "[UNRESOLVED]",
                "finite_scope": (
                    "No coupling interval, all-size theorem, embedding theorem, thermodynamic-limit "
                    "claim, or exact solution of the three-dimensional Ising model follows."
                ),
                "escape_routes": [
                    "parity-projected or paired Gaussian descriptions not equal to the full 1024-value subset-product spectrum",
                    "restriction of a larger Gaussian operator",
                    "different isolated couplings outside t in {1/3,1/2}",
                    "non-positive or complex-spectrum Gaussian notions outside the lemma's hypotheses",
                ],
                "control_caveat": (
                    "An occupied chain window proves only that some eigenvalue lies in the finite "
                    "window; it does not prove exact equality to lambda1*lambda2/lambda0."
                ),
            },
        },
        "checks": checks,
    }

    output = ROOT / "results" / "spectral" / "spectral_2x5.json"
    os.makedirs(output.parent, exist_ok=True)
    with output.open("w") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    for check in checks:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    if not all(check["passed"] for check in checks):
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
