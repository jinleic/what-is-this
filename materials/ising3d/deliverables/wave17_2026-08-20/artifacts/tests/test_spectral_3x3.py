"""Standalone exact re-derivation of the 3x3, t=1/2 spectral certificate.

This test intentionally does not import ``experiments.e54_spectral_3x3``.  It
independently rebuilds the 512-dimensional integer transfer matrix, all exact
D4 x spin-flip projectors (including the two E rows), their generalized
sector pencils, and the forced-value inertia certificate.  It then validates
the stored two-coupling artifact.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


@dataclass
class Projector:
    label: str
    denominator: int
    columns: list[dict[int, int]]


@dataclass
class Vector:
    coefficients: dict[int, int]
    orbit: int


@dataclass
class Sector:
    label: str
    matrix: list[list[int]]
    gram: list[list[int]]


# ---------------------------------------------------------- independent exact construction


def build_integer_operator(t: Fraction) -> tuple[list[list[int]], int, Fraction]:
    n = 9
    bonds = list(layer_bonds((3, 3), (False, False)))
    dimension = 1 << n
    q = (1 + t * t) / (2 * t)
    energies: list[int] = []
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    assert all((energy - parity) % 2 == 0 for energy in energies)
    exponents = [(energy - parity) // 2 for energy in energies]
    weight_scale = (
        q.denominator ** max(max(exponents), 0)
        * q.numerator ** max(-min(exponents), 0)
    )
    weights = []
    for exponent in exponents:
        weight = q**exponent * weight_scale
        assert weight.denominator == 1
        weights.append(weight.numerator)

    numerator, denominator = t.numerator, t.denominator
    numerator_powers = [numerator**power for power in range(n + 1)]
    denominator_powers = [denominator**power for power in range(n + 1)]
    p_numerator = [
        [
            numerator_powers[(row ^ col).bit_count()]
            * denominator_powers[n - (row ^ col).bit_count()]
            for col in range(dimension)
        ]
        for row in range(dimension)
    ]
    matrix = [[weights[row] * entry for entry in p_numerator[row]] for row in range(dimension)]
    for bit in range(n):
        mask = 1 << bit
        for base in range(dimension):
            if base & mask:
                continue
            mate = base | mask
            first, second = matrix[base], matrix[mate]
            matrix[base] = [
                denominator * left + numerator * right
                for left, right in zip(first, second)
            ]
            matrix[mate] = [
                numerator * left + denominator * right
                for left, right in zip(first, second)
            ]
    assert all(
        matrix[row][col] == matrix[col][row]
        for row in range(dimension)
        for col in range(row)
    )
    return matrix, denominator ** (2 * n) * weight_scale, q


def state_permutation(site_permutation: list[int]) -> list[int]:
    result: list[int] = []
    for state in range(512):
        image = 0
        for source, target in enumerate(site_permutation):
            if (state >> (8 - source)) & 1:
                image |= 1 << (8 - target)
        result.append(image)
    return result


def d4_actions() -> list[tuple[int, bool, bool, list[int]]]:
    def rotate(row: int, col: int, power: int) -> tuple[int, int]:
        for _ in range(power):
            row, col = col, 2 - row
        return row, col

    result: list[tuple[int, bool, bool, list[int]]] = []
    for reflected in (False, True):
        for power in range(4):
            sites = []
            for row in range(3):
                for col in range(3):
                    target_row, target_col = rotate(row, col, power)
                    if reflected:
                        target_col = 2 - target_col
                    sites.append(3 * target_row + target_col)
            spatial = state_permutation(sites)
            for flipped in (False, True):
                mapping = [value ^ 511 if flipped else value for value in spatial]
                result.append((power, reflected, flipped, mapping))
    assert len({tuple(action[3]) for action in result}) == 16
    return result


def character(irrep: str, power: int, reflected: bool) -> int:
    if irrep == "A1":
        return 1
    if irrep == "A2":
        return -1 if reflected else 1
    if irrep == "B1":
        return -1 if power % 2 else 1
    if irrep == "B2":
        rotation_sign = -1 if power % 2 else 1
        return -rotation_sign if reflected else rotation_sign
    if irrep == "E":
        return 0 if reflected or power % 2 else (2 if power == 0 else -2)
    raise ValueError(irrep)


def make_columns(terms: list[tuple[int, list[int]]]) -> list[dict[int, int]]:
    columns = []
    for state in range(512):
        column: dict[int, int] = {}
        for coefficient, mapping in terms:
            image = mapping[state]
            column[image] = column.get(image, 0) + coefficient
        columns.append({image: value for image, value in column.items() if value})
    return columns


def exact_projectors(
    actions: list[tuple[int, bool, bool, list[int]]]
) -> list[Projector]:
    projectors: list[Projector] = []
    axial = next(
        mapping
        for power, reflected, flipped, mapping in actions
        if power == 0 and reflected and not flipped
    )
    for irrep, irrep_dimension in (
        ("A1", 1),
        ("A2", 1),
        ("B1", 1),
        ("B2", 1),
        ("E", 2),
    ):
        for flip_sign in (1, -1):
            central_terms = []
            for power, reflected, flipped, mapping in actions:
                coefficient = irrep_dimension * character(irrep, power, reflected)
                if flipped:
                    coefficient *= flip_sign
                if coefficient:
                    central_terms.append((coefficient, mapping))
            if irrep != "E":
                terms = [(2 * coefficient, mapping) for coefficient, mapping in central_terms]
                label = f"{irrep}_flip{'+' if flip_sign == 1 else '-'}"
                projectors.append(Projector(label, 32, make_columns(terms)))
                continue
            for axial_sign in (1, -1):
                terms = central_terms + [
                    (
                        axial_sign * coefficient,
                        [axial[value] for value in mapping],
                    )
                    for coefficient, mapping in central_terms
                ]
                label = (
                    f"E_axis{'+' if axial_sign == 1 else '-'}_"
                    f"flip{'+' if flip_sign == 1 else '-'}"
                )
                projectors.append(Projector(label, 32, make_columns(terms)))
    return projectors


def compose(
    left_columns: list[dict[int, int]], right: dict[int, int]
) -> dict[int, int]:
    result: dict[int, int] = {}
    for middle, right_value in right.items():
        for image, left_value in left_columns[middle].items():
            result[image] = result.get(image, 0) + left_value * right_value
    return {image: value for image, value in result.items() if value}


def verify_projector_family(
    matrix: list[list[int]],
    actions: list[tuple[int, bool, bool, list[int]]],
    projectors: list[Projector],
) -> dict[str, bool]:
    self_adjoint = all(
        projector.columns[target].get(source, 0) == coefficient
        for projector in projectors
        for source, column in enumerate(projector.columns)
        for target, coefficient in column.items()
    )
    idempotent = all(
        compose(projector.columns, column)
        == {image: 32 * value for image, value in column.items()}
        for projector in projectors
        for column in projector.columns
    )
    orthogonal = all(
        not compose(left.columns, column)
        for left_index, left in enumerate(projectors)
        for right in projectors[left_index + 1 :]
        for column in right.columns
    )
    complete = True
    for state in range(512):
        total: dict[int, int] = {}
        for projector in projectors:
            for image, coefficient in projector.columns[state].items():
                total[image] = total.get(image, 0) + coefficient
        if {image: value for image, value in total.items() if value} != {state: 32}:
            complete = False
            break
    invariant = all(
        matrix[mapping[row]][mapping[col]] == matrix[row][col]
        for _, _, _, mapping in actions
        for row in range(512)
        for col in range(512)
    )
    return {
        "self_adjoint": self_adjoint,
        "idempotent": idempotent,
        "orthogonal": orthogonal,
        "complete": complete,
        "commutes": invariant,
    }


# ---------------------------------------------------------- independent sector restrictions


def normalize(vector: dict[int, int]) -> dict[int, int]:
    if not vector:
        return {}
    common = 0
    for value in vector.values():
        common = math.gcd(common, abs(value))
    result = {state: value // common for state, value in vector.items()}
    if result[min(result)] < 0:
        result = {state: -value for state, value in result.items()}
    return result


def independent_columns(
    orbit: list[int], columns: list[dict[int, int]]
) -> list[dict[int, int]]:
    position = {state: index for index, state in enumerate(orbit)}
    echelon: dict[int, list[Fraction]] = {}
    selected: list[dict[int, int]] = []
    for state in orbit:
        candidate = normalize(columns[state])
        row = [Fraction(candidate.get(site, 0)) for site in orbit]
        for pivot in sorted(echelon):
            if row[pivot]:
                factor = row[pivot]
                row = [left - factor * right for left, right in zip(row, echelon[pivot])]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        pivot_value = row[pivot]
        echelon[pivot] = [value / pivot_value for value in row]
        selected.append(candidate)
    return selected


def projected_basis(
    actions: list[tuple[int, bool, bool, list[int]]], projector: Projector
) -> list[Vector]:
    seen: set[int] = set()
    result: list[Vector] = []
    orbit_index = 0
    for representative in range(512):
        if representative in seen:
            continue
        orbit = sorted({mapping[representative] for _, _, _, mapping in actions})
        seen.update(orbit)
        result.extend(
            Vector(coefficients, orbit_index)
            for coefficients in independent_columns(orbit, projector.columns)
        )
        orbit_index += 1
    return result


def build_sectors(
    matrix: list[list[int]],
    actions: list[tuple[int, bool, bool, list[int]]],
    projectors: list[Projector],
) -> list[Sector]:
    sectors: list[Sector] = []
    for projector in projectors:
        basis = projected_basis(actions, projector)
        size = len(basis)
        transformed = [
            [
                sum(matrix[row][site] * value for site, value in vector.coefficients.items())
                for row in range(512)
            ]
            for vector in basis
        ]
        restricted = [[0] * size for _ in range(size)]
        gram = [[0] * size for _ in range(size)]
        for row in range(size):
            for col in range(row + 1):
                left, right = basis[row], basis[col]
                restricted_value = sum(
                    value * transformed[col][site]
                    for site, value in left.coefficients.items()
                )
                gram_value = (
                    0
                    if left.orbit != right.orbit
                    else sum(
                        value * right.coefficients.get(site, 0)
                        for site, value in left.coefficients.items()
                    )
                )
                restricted[row][col] = restricted[col][row] = restricted_value
                gram[row][col] = gram[col][row] = gram_value
        trace = sum(
            projector.columns[state].get(state, 0) for state in range(512)
        ) // 32
        assert size == trace
        sectors.append(Sector(projector.label, restricted, gram))
    assert sum(len(sector.gram) for sector in sectors) == 512
    return sectors


# --------------------------------------------------------------- independent exact inertia


def bareiss_inertia(lower: list[list[int]]) -> int | None:
    work = [row[:] for row in lower]
    previous = 1
    negatives = 0
    for pivot_index in range(len(work)):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        if pivot_index + 1 == len(work):
            break
        for row_index in range(pivot_index + 1, len(work)):
            row = work[row_index]
            left = row[pivot_index]
            for col_index in range(pivot_index + 1, row_index + 1):
                numerator = pivot * row[col_index] - left * work[col_index][pivot_index]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                row[col_index] = quotient
        previous = pivot
    return negatives


def sector_count(sectors: list[Sector], shift: Fraction) -> tuple[int, dict[str, int]]:
    total = 0
    details: dict[str, int] = {}
    for sector in sectors:
        size = len(sector.gram)
        shifted = [[0] * size for _ in range(size)]
        common = 0
        for row in range(size):
            for col in range(row + 1):
                value = (
                    shift.denominator * sector.matrix[row][col]
                    - shift.numerator * sector.gram[row][col]
                )
                shifted[row][col] = value
                common = math.gcd(common, abs(value))
        if common > 1:
            for row in range(size):
                for col in range(row + 1):
                    shifted[row][col] //= common
        count = bareiss_inertia(shifted)
        assert count is not None
        details[sector.label] = count
        total += count
    return total, details


def rederive_t_one_half() -> dict[str, object]:
    matrix, scale, q = build_integer_operator(Fraction(1, 2))
    actions = d4_actions()
    projectors = exact_projectors(actions)
    projector_checks = verify_projector_family(matrix, actions, projectors)
    assert all(projector_checks.values())
    sectors = build_sectors(matrix, actions, projectors)

    largest = max(abs(value) for row in matrix for value in row)
    proposal = np.array([[value / largest for value in row] for row in matrix], dtype=float)
    guides = [float(value * largest) for value in np.linalg.eigvalsh(proposal)[:3]]
    enclosures: list[tuple[Fraction, Fraction]] = []
    endpoint_counts: list[tuple[tuple[int, dict[str, int]], tuple[int, dict[str, int]]]] = []
    for index, guide in enumerate(guides):
        center = Fraction(format(guide, ".16e"))
        lower, upper = center * Fraction(999, 1000), center * Fraction(1001, 1000)
        lower_count, upper_count = sector_count(sectors, lower), sector_count(sectors, upper)
        assert lower_count[0] <= index < upper_count[0]
        enclosures.append((lower, upper))
        endpoint_counts.append((lower_count, upper_count))

    (l0, h0), (l1, h1), (l2, h2) = enclosures
    predicted_lower, predicted_upper = l1 * l2 / h0, h1 * h2 / l0
    count_lower = sector_count(sectors, predicted_lower)
    count_upper = sector_count(sectors, predicted_upper)
    return {
        "q": q,
        "scale": scale,
        "projector_checks": projector_checks,
        "sector_dimensions": {sector.label: len(sector.gram) for sector in sectors},
        "positive_count": sector_count(sectors, Fraction(0))[0],
        "enclosures": enclosures,
        "endpoint_counts": endpoint_counts,
        "distinct": h0 < l1 and h1 < l2,
        "predicted": (predicted_lower, predicted_upper),
        "window_counts": (count_lower[0], count_upper[0]),
        "window_sector_counts": (count_lower[1], count_upper[1]),
    }


print("1. independent exact 3x3, t=1/2 re-derivation")
independent = rederive_t_one_half()
check("independent exp(2K)=5/4", independent["q"] == Fraction(5, 4))
check("independent positive definiteness", independent["positive_count"] == 0)
check("independent exact projector identities", all(independent["projector_checks"].values()))
check(
    "independent refined sector dimensions",
    independent["sector_dimensions"]
    == {
        "A1_flip+": 51,
        "A1_flip-": 51,
        "A2_flip+": 19,
        "A2_flip-": 19,
        "B1_flip+": 33,
        "B1_flip-": 33,
        "B2_flip+": 33,
        "B2_flip-": 33,
        "E_axis+_flip+": 60,
        "E_axis-_flip+": 60,
        "E_axis+_flip-": 60,
        "E_axis-_flip-": 60,
    },
    str(independent["sector_dimensions"]),
)
check("independent three lowest distinct", independent["distinct"])
check(
    "independent forced window absent",
    independent["window_counts"] == (11, 11),
    f"exact counts={independent['window_counts']}",
)
check(
    "independent per-sector absence",
    independent["window_sector_counts"][0] == independent["window_sector_counts"][1],
)


print("2. stored artifact validation")
artifact_path = ROOT / "results" / "spectral" / "spectral_3x3.json"
with artifact_path.open() as handle:
    artifact = json.load(handle)
check(
    "artifact envelope",
    set(artifact) == {"provenance", "data", "checks"}
    and artifact["provenance"]["script"] == "experiments/e54_spectral_3x3.py"
    and artifact["provenance"]["interpreter"] == ".venv/bin/python",
)
rows = {row["name"]: row for row in artifact["data"]["rows"]}
check(
    "artifact contains exactly four named cases",
    set(rows)
    == {
        "layer_3x3_t_1_3",
        "layer_3x3_t_1_2",
        "control_chain_n9_t_1_3",
        "control_chain_n9_t_1_2",
    },
)
for name, expected_counts, expected_absent, expected_tag in (
    ("layer_3x3_t_1_3", [8, 8], True, "[THEOREM]"),
    ("layer_3x3_t_1_2", [11, 11], True, "[THEOREM]"),
    ("control_chain_n9_t_1_3", [10, 11], False, "[COMPUTATION]"),
    ("control_chain_n9_t_1_2", [10, 11], False, "[COMPUTATION]"),
):
    row = rows[name]
    check(
        f"stored verdict {name}",
        row["claim_tag"] == expected_tag
        and row["three_lowest_distinct"]
        and row["inertia_counts_at_forced_window_ends"] == expected_counts
        and row["predicted_eigenvalue_absent"] is expected_absent,
    )
    check(
        f"stored exact intervals {name}",
        all(
            Fraction(interval["lower"]) < Fraction(interval["upper"])
            for interval in row["lambda_enclosures"].values()
        )
        and Fraction(row["forced_value_window"]["predicted"]["lower"])
        <= Fraction(row["forced_value_window"]["predicted"]["upper"])
        and row["forced_value_window"]["tested_window_encloses_prediction"],
    )
    check(
        f"stored projector and dimension checks {name}",
        row["symmetry"]["sector_dimensions_sum"] == 512
        and all(
            row["symmetry"]["verification"][key]
            for key in (
                "projector_self_adjoint_exact",
                "projector_idempotence_exact",
                "projector_orthogonality_exact",
                "projector_completeness_exact",
                "projectors_commute_with_matrix_exact",
            )
        ),
    )

stored_half = rows["layer_3x3_t_1_2"]
check(
    "stored t=1/2 agrees with independent sector dimensions",
    stored_half["symmetry"]["sector_dimensions"] == independent["sector_dimensions"],
)
check(
    "stored t=1/2 exact scale and q",
    int(stored_half["integer_scale_for_R"]) == independent["scale"]
    and Fraction(stored_half["exp_2K"]) == independent["q"],
)

dense = artifact["data"]["full_dense_verification"]
check(
    "stored full-dense exact residual inequalities",
    dense["dense_counts_match_sector_sums"]
    and dense["coarse_window_contains_predicted_window"]
    and all(
        witness["inertia_certified_exactly"]
        and int(witness["strict_safety_margin_numerator"]) > 0
        and Fraction(witness["residual_frobenius_upper_bound"])
        < Fraction(witness["sigma_min_model_lower_bound"])
        for witness in dense["witnesses"]
    ),
)
check(
    "all stored top-level checks pass",
    all(item["passed"] for item in artifact["checks"]),
)
check(
    "honest finite scope and resource record",
    artifact["data"]["claim_tag"] == "[THEOREM]"
    and artifact["data"]["limitations"]["claim_tag"] == "[UNRESOLVED]"
    and dense["bareiss_resource_record"]["claim_tag"] == "[UNRESOLVED]"
    and dense["bareiss_resource_record"]["observed_wall_seconds_lower_bound"] == 1200,
)

if FAILURES:
    print(f"FAIL: {FAILURES}")
    raise SystemExit(1)
print("PASS")
