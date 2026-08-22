"""Standalone exact re-derivation of the open 2x5, t=1/3 certificate.

This test deliberately does not import either spectral experiment.  It rebuilds
the 1024-dimensional integer operator, the eight rectangle character
projectors, their exact sector pencils, and an independent integer-congruence
inertia witness at every stored endpoint for one coupling.
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
DYADIC_BITS = 36


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


@dataclass
class IndependentProjector:
    label: str
    columns: list[dict[int, int]]


@dataclass
class IndependentVector:
    coefficients: dict[int, int]
    orbit_index: int


@dataclass
class IndependentSector:
    label: str
    matrix: list[list[int]]
    gram_diagonal: list[int]


@dataclass
class PreparedSector:
    label: str
    operator_congruence: np.ndarray
    gram_congruence: np.ndarray


# -------------------------------------------------------- independent exact operator


def build_integer_operator(t: Fraction) -> tuple[list[list[int]], int, Fraction]:
    n = 10
    dimension = 1 << n
    bonds = list(layer_bonds((2, 5), (False, False)))
    q = (1 + t * t) / (2 * t)
    energies: list[int] = []
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    assert parity == 1
    assert all((energy - parity) % 2 == 0 for energy in energies)
    exponents = [(energy - parity) // 2 for energy in energies]
    weight_scale = (
        q.denominator ** max(max(exponents), 0)
        * q.numerator ** max(-min(exponents), 0)
    )
    weights: list[int] = []
    for exponent in exponents:
        weight = q**exponent * weight_scale
        assert weight.denominator == 1
        weights.append(weight.numerator)

    numerator, denominator = t.numerator, t.denominator
    numerator_powers = [numerator**power for power in range(n + 1)]
    denominator_powers = [denominator**power for power in range(n + 1)]
    p_numerator = [
        [
            numerator_powers[(row ^ column).bit_count()]
            * denominator_powers[n - (row ^ column).bit_count()]
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    matrix = [
        [weights[row] * entry for entry in p_numerator[row]]
        for row in range(dimension)
    ]
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
        matrix[row][column] == matrix[column][row]
        for row in range(dimension)
        for column in range(row)
    )
    return matrix, denominator ** (2 * n) * weight_scale, q


# ------------------------------------------------ independent C2 x C2 x C2 projectors


def state_permutation(site_permutation: list[int]) -> list[int]:
    assert sorted(site_permutation) == list(range(10))
    result: list[int] = []
    for state in range(1024):
        image = 0
        for source, target in enumerate(site_permutation):
            if (state >> (9 - source)) & 1:
                image |= 1 << (9 - target)
        result.append(image)
    return result


def rectangle_actions() -> list[tuple[bool, bool, bool, list[int]]]:
    result: list[tuple[bool, bool, bool, list[int]]] = []
    for row_reflected in (False, True):
        for column_reflected in (False, True):
            sites: list[int] = []
            for row in range(2):
                for column in range(5):
                    target_row = 1 - row if row_reflected else row
                    target_column = 4 - column if column_reflected else column
                    sites.append(5 * target_row + target_column)
            spatial = state_permutation(sites)
            for flipped in (False, True):
                mapping = [value ^ 1023 if flipped else value for value in spatial]
                result.append((row_reflected, column_reflected, flipped, mapping))
    assert len({tuple(action[3]) for action in result}) == 8
    return result


def make_columns(
    terms: list[tuple[int, list[int]]]
) -> list[dict[int, int]]:
    columns: list[dict[int, int]] = []
    for state in range(1024):
        column: dict[int, int] = {}
        for coefficient, mapping in terms:
            image = mapping[state]
            column[image] = column.get(image, 0) + coefficient
        columns.append({image: value for image, value in column.items() if value})
    return columns


def character_projectors(
    actions: list[tuple[bool, bool, bool, list[int]]]
) -> list[IndependentProjector]:
    result: list[IndependentProjector] = []
    for row_sign in (1, -1):
        for column_sign in (1, -1):
            for flip_sign in (1, -1):
                terms: list[tuple[int, list[int]]] = []
                for row_reflected, column_reflected, flipped, mapping in actions:
                    value = row_sign if row_reflected else 1
                    if column_reflected:
                        value *= column_sign
                    if flipped:
                        value *= flip_sign
                    terms.append((value, mapping))
                label = (
                    f"row{'+' if row_sign == 1 else '-'}_"
                    f"column{'+' if column_sign == 1 else '-'}_"
                    f"flip{'+' if flip_sign == 1 else '-'}"
                )
                result.append(IndependentProjector(label, make_columns(terms)))
    return result


def compose(
    left_columns: list[dict[int, int]], right_column: dict[int, int]
) -> dict[int, int]:
    result: dict[int, int] = {}
    for middle, right_value in right_column.items():
        for image, left_value in left_columns[middle].items():
            result[image] = result.get(image, 0) + left_value * right_value
    return {image: value for image, value in result.items() if value}


def verify_projectors(
    matrix: list[list[int]],
    actions: list[tuple[bool, bool, bool, list[int]]],
    projectors: list[IndependentProjector],
) -> dict[str, bool]:
    self_adjoint = all(
        projector.columns[target].get(source, 0) == coefficient
        for projector in projectors
        for source, column in enumerate(projector.columns)
        for target, coefficient in column.items()
    )
    idempotent = all(
        compose(projector.columns, column)
        == {image: 8 * value for image, value in column.items()}
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
    for state in range(1024):
        total: dict[int, int] = {}
        for projector in projectors:
            for image, value in projector.columns[state].items():
                total[image] = total.get(image, 0) + value
        if {image: value for image, value in total.items() if value} != {state: 8}:
            complete = False
            break
    commutes = all(
        matrix[mapping[row]][mapping[column]] == matrix[row][column]
        for _, _, _, mapping in actions
        for row in range(1024)
        for column in range(1024)
    )
    return {
        "self_adjoint": self_adjoint,
        "idempotent": idempotent,
        "orthogonal": orthogonal,
        "complete": complete,
        "commutes": commutes,
    }


# ---------------------------------------------------- independent sector restrictions


def normalize(vector: dict[int, int]) -> dict[int, int]:
    common = 0
    for value in vector.values():
        common = math.gcd(common, abs(value))
    result = {state: value // common for state, value in vector.items()}
    if result and result[min(result)] < 0:
        result = {state: -value for state, value in result.items()}
    return result


def projected_basis(
    actions: list[tuple[bool, bool, bool, list[int]]],
    projector: IndependentProjector,
) -> list[IndependentVector]:
    """One nonzero character vector spans each orbit/character intersection."""
    seen: set[int] = set()
    result: list[IndependentVector] = []
    orbit_index = 0
    for representative in range(1024):
        if representative in seen:
            continue
        orbit = sorted({mapping[representative] for _, _, _, mapping in actions})
        seen.update(orbit)
        vector = normalize(projector.columns[representative])
        if vector:
            assert set(vector).issubset(orbit)
            result.append(IndependentVector(vector, orbit_index))
        orbit_index += 1
    return result


def sparse_dot(left: dict[int, int], right: dict[int, int]) -> int:
    return sum(value * right.get(state, 0) for state, value in left.items())


def build_sectors(
    matrix: list[list[int]],
    actions: list[tuple[bool, bool, bool, list[int]]],
    projectors: list[IndependentProjector],
) -> list[IndependentSector]:
    sectors: list[IndependentSector] = []
    for projector in projectors:
        basis = projected_basis(actions, projector)
        size = len(basis)
        transformed = [
            [
                sum(matrix[row][state] * value for state, value in vector.coefficients.items())
                for row in range(1024)
            ]
            for vector in basis
        ]
        restricted = [[0] * size for _ in range(size)]
        gram_diagonal: list[int] = []
        for row, left in enumerate(basis):
            gram_diagonal.append(sparse_dot(left.coefficients, left.coefficients))
            for column in range(row + 1):
                value = sum(
                    coefficient * transformed[column][state]
                    for state, coefficient in left.coefficients.items()
                )
                restricted[row][column] = restricted[column][row] = value
        trace = sum(
            projector.columns[state].get(state, 0) for state in range(1024)
        ) // 8
        assert trace == size
        sectors.append(IndependentSector(projector.label, restricted, gram_diagonal))
    assert sum(len(sector.gram_diagonal) for sector in sectors) == 1024
    return sectors


# ------------------------------------------ independent exact congruence inertia witness


def prepare_sector(sector: IndependentSector) -> PreparedSector:
    size = len(sector.gram_diagonal)
    largest = max(abs(value) for row in sector.matrix for value in row)
    inverse_roots = np.array(
        [1.0 / math.sqrt(value) for value in sector.gram_diagonal]
    )
    proposal = np.empty((size, size), dtype=float)
    for row in range(size):
        for column in range(size):
            proposal[row, column] = (
                sector.matrix[row][column]
                / largest
                * inverse_roots[row]
                * inverse_roots[column]
            )
    _, eigenvectors = np.linalg.eigh(0.5 * (proposal + proposal.T))
    q_integer = np.rint(
        inverse_roots[:, None] * eigenvectors * (1 << DYADIC_BITS)
    ).astype(np.int64).astype(object)
    operator = np.array(sector.matrix, dtype=object)
    operator_congruence = q_integer.T @ (operator @ q_integer)
    weighted_q = q_integer.copy()
    for row, gram_value in enumerate(sector.gram_diagonal):
        weighted_q[row, :] *= gram_value
    gram_congruence = q_integer.T @ weighted_q
    assert all(
        int(operator_congruence[row, column])
        == int(operator_congruence[column, row])
        and int(gram_congruence[row, column])
        == int(gram_congruence[column, row])
        for row in range(size)
        for column in range(row)
    )
    return PreparedSector(sector.label, operator_congruence, gram_congruence)


def exact_sector_count(
    prepared: list[PreparedSector], shift: Fraction
) -> tuple[int, dict[str, int], dict[str, int]]:
    total = 0
    details: dict[str, int] = {}
    margins: dict[str, int] = {}
    for sector in prepared:
        size = len(sector.operator_congruence)
        diagonal: list[int] = []
        off_diagonal: list[tuple[int, int, int]] = []
        for row in range(size):
            for column in range(row + 1):
                entry = (
                    shift.denominator * int(sector.operator_congruence[row, column])
                    - shift.numerator * int(sector.gram_congruence[row, column])
                )
                if row == column:
                    diagonal.append(entry)
                else:
                    off_diagonal.append((row, column, entry))
        assert all(diagonal)
        exponents = [abs(value).bit_length() - 1 for value in diagonal]
        common_exponent = max(
            (exponents[row] + exponents[column] for row, column, _ in off_diagonal),
            default=0,
        )
        numerator = sum(
            2
            * entry
            * entry
            * (1 << (common_exponent - exponents[row] - exponents[column]))
            for row, column, entry in off_diagonal
        )
        denominator = 1 << common_exponent
        margin = denominator - numerator
        assert margin > 0
        inertia = sum(value < 0 for value in diagonal)
        details[sector.label] = inertia
        margins[sector.label] = margin
        total += inertia
    return total, details, margins


def artifact_raw_shifts(row: dict[str, object]) -> list[Fraction]:
    scale = int(row["integer_scale_for_R"])
    shifts = [Fraction(0)]
    for interval in row["lambda_enclosures"].values():
        shifts.extend(
            [
                Fraction(interval["lower"]) * scale,
                Fraction(interval["upper"]) * scale,
            ]
        )
    forced = row["forced_value_window"]["actually_tested"]
    shifts.extend(
        [Fraction(forced["lower"]) * scale, Fraction(forced["upper"]) * scale]
    )
    return shifts


print("1. load stored endpoint definitions")
artifact_path = ROOT / "results" / "spectral" / "spectral_2x5.json"
with artifact_path.open() as handle:
    artifact = json.load(handle)
rows = {row["name"]: row for row in artifact["data"]["rows"]}
stored = rows["layer_2x5_t_1_3"]

print("2. independent exact 2x5, t=1/3 re-derivation")
matrix, scale, q = build_integer_operator(Fraction(1, 3))
actions = rectangle_actions()
projectors = character_projectors(actions)
projector_checks = verify_projectors(matrix, actions, projectors)
sectors = build_sectors(matrix, actions, projectors)
prepared = [prepare_sector(sector) for sector in sectors]
sector_dimensions = {
    sector.label: len(sector.gram_diagonal) for sector in sectors
}
independent_counts = {
    shift: exact_sector_count(prepared, shift) for shift in artifact_raw_shifts(stored)
}

check("independent exp(2K)=5/3", q == Fraction(5, 3))
check(
    "independent integer scale",
    scale == 198583267838203125 == int(stored["integer_scale_for_R"]),
)
check(
    "independent projector self-adjointness/idempotence/orthogonality/completeness/commutation",
    all(projector_checks.values()),
    str(projector_checks),
)
check(
    "independent sector dimensions",
    sector_dimensions
    == {
        "row+_column+_flip+": 152,
        "row+_column+_flip-": 136,
        "row+_column-_flip+": 120,
        "row+_column-_flip-": 120,
        "row-_column+_flip+": 120,
        "row-_column+_flip-": 136,
        "row-_column-_flip+": 120,
        "row-_column-_flip-": 120,
    }
    == stored["symmetry"]["sector_dimensions"],
    str(sector_dimensions),
)
check(
    "independent positivity",
    independent_counts[Fraction(0)][0] == 0,
)

raw_shifts = artifact_raw_shifts(stored)
for index in range(3):
    lower_shift, upper_shift = raw_shifts[1 + 2 * index : 3 + 2 * index]
    lower_count = independent_counts[lower_shift]
    upper_count = independent_counts[upper_shift]
    stored_endpoint = stored["lambda_enclosures"][f"lambda{index}"]["endpoint_inertia"]
    check(
        f"independent lambda{index} exact endpoint counts",
        lower_count[0] == index
        and upper_count[0] == index + 1
        and lower_count[1] == stored_endpoint["lower"]["by_sector"]
        and upper_count[1] == stored_endpoint["upper"]["by_sector"],
        f"counts={lower_count[0]},{upper_count[0]}",
    )

forced_lower, forced_upper = raw_shifts[-2:]
forced_counts = independent_counts[forced_lower], independent_counts[forced_upper]
check(
    "independent forced window absent",
    (forced_counts[0][0], forced_counts[1][0]) == (8, 8)
    and forced_counts[0][1] == forced_counts[1][1]
    and stored["inertia_counts_at_forced_window_ends"] == [8, 8],
    f"exact counts={forced_counts[0][0], forced_counts[1][0]}",
)
check(
    "independent strict integer congruence margins",
    all(
        margin > 0
        for _, _, margins in independent_counts.values()
        for margin in margins.values()
    ),
)
check(
    "independent rational interval logic",
    all(
        Fraction(interval["lower"]) < Fraction(interval["upper"])
        for interval in stored["lambda_enclosures"].values()
    )
    and Fraction(stored["lambda_enclosures"]["lambda0"]["upper"])
    < Fraction(stored["lambda_enclosures"]["lambda1"]["lower"])
    and Fraction(stored["lambda_enclosures"]["lambda1"]["upper"])
    < Fraction(stored["lambda_enclosures"]["lambda2"]["lower"])
    and Fraction(stored["forced_value_window"]["predicted"]["lower"])
    <= Fraction(stored["forced_value_window"]["actually_tested"]["lower"])
    <= Fraction(stored["forced_value_window"]["actually_tested"]["upper"])
    <= Fraction(stored["forced_value_window"]["predicted"]["upper"]),
)

print("3. stored four-case artifact validation")
check(
    "artifact envelope and provenance",
    set(artifact) == {"provenance", "data", "checks"}
    and artifact["provenance"]["script"] == "experiments/e67_spectral_2x5.py"
    and artifact["provenance"]["interpreter"] == ".venv/bin/python",
)
check(
    "artifact contains exactly four requested cases",
    set(rows)
    == {
        "layer_2x5_t_1_3",
        "layer_2x5_t_1_2",
        "control_chain_n10_t_1_3",
        "control_chain_n10_t_1_2",
    },
)
for name, expected_counts, expected_absent, expected_tag in (
    ("layer_2x5_t_1_3", [8, 8], True, "[THEOREM]"),
    ("layer_2x5_t_1_2", [12, 12], True, "[THEOREM]"),
    ("control_chain_n10_t_1_3", [11, 12], False, "[COMPUTATION]"),
    ("control_chain_n10_t_1_2", [11, 12], False, "[COMPUTATION]"),
):
    row = rows[name]
    check(
        f"stored exact verdict {name}",
        row["claim_tag"] == expected_tag
        and row["dimension"] == 1024
        and row["three_lowest_distinct"]
        and row["inertia_counts_at_forced_window_ends"] == expected_counts
        and row["predicted_eigenvalue_absent"] is expected_absent,
    )
    check(
        f"stored complete sector certificate {name}",
        row["symmetry"]["sector_dimensions_sum"] == 1024
        and sum(row["symmetry"]["sector_dimensions"].values()) == 1024
        and all(
            row["symmetry"]["verification"][key]
            for key in (
                "projector_self_adjoint_exact",
                "projector_idempotence_exact",
                "projector_orthogonality_exact",
                "projector_completeness_exact",
                "projectors_commute_with_matrix_exact",
            )
        )
        and all(
            witness["inertia_certified_exactly"]
            and int(witness["relative_frobenius_squared_upper_numerator"])
            < int(witness["relative_frobenius_squared_upper_denominator"])
            and int(witness["strict_safety_margin_numerator"]) > 0
            for shift in row["exact_sector_inertia_certificate"]["shift_witnesses"].values()
            for witness in shift["sectors"].values()
        ),
    )

check(
    "all top-level checks pass",
    all(item["passed"] for item in artifact["checks"]),
)
check(
    "honest finite scope resource record and escape routes",
    artifact["data"]["claim_tag"] == "[THEOREM]"
    and artifact["data"]["limitations"]["claim_tag"] == "[UNRESOLVED]"
    and len(artifact["data"]["limitations"]["escape_routes"]) >= 4
    and artifact["data"]["resource_record"]["total_wall_seconds"] > 0
    and artifact["data"]["resource_record"]["process_max_rss_raw"] > 0
    and artifact["data"]["resource_record"]["dense_1024_check_performed"] is False,
)

if FAILURES:
    print(f"FAIL: {FAILURES}")
    raise SystemExit(1)
print("PASS")
