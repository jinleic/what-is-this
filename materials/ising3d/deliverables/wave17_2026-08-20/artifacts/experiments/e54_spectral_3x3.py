"""Exact spectral-Gaussianity obstruction for the open 3x3 Ising layer.

For rational ``t = tanh(K*/2)``, this script constructs the integer multiple
``A = scale * P_t diag(q**e) P_t`` of the symmetrized layer transfer operator,
where ``q = (1+t*t)/(2*t)``.  Every spectral decision is an exact integer or
rational certificate.

The 512-dimensional layer is split by all ten central projectors of
``D4 x C2(spin flip)``.  Projector idempotence, mutual orthogonality,
completeness, and invariance of the matrix are checked entry by entry over the
integers.  Exact fraction-free LDL inertia in the resulting generalized
sector pencils encloses the three lowest eigenvalues and proves that the
Gaussian-forced value ``lambda1*lambda2/lambda0`` is absent.

NumPy float64 eigenvectors/eigenvalues only propose rational brackets and a
rational change of basis.  Every accepted bracket has exact inertia endpoint
counts.  The independent full-512 bookkeeping check is an exact dyadic
residual certificate: after rounding the proposed eigenbasis to integers, an
integer Frobenius bound proves the inertia of the full dense shifted matrix.
Thus no float64 result decides a stored claim.
"""

from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402


@dataclass
class ScaledOperator:
    matrix: list[list[int]]
    scale: int
    q: Fraction
    bond_parity: int
    energy_exponents: tuple[int, int]


@dataclass
class Action:
    name: str
    state_map: list[int]
    rotation_power: int = 0
    reflected: bool = False
    flipped: bool = False


@dataclass
class Projector:
    label: str
    denominator: int
    terms: list[tuple[int, list[int]]]
    columns: list[dict[int, int]]
    irrep_dimension: int
    flip_sign: int
    character_values: list[int]


@dataclass
class BasisVector:
    coefficients: dict[int, int]
    orbit_index: int


@dataclass
class Sector:
    label: str
    matrix: list[list[int]]
    gram: list[list[int]]
    dimension: int
    projector_trace: int
    irrep_dimension: int
    flip_sign: int


@dataclass
class CountRecord:
    total: int | None
    by_sector: dict[str, int] | None
    seconds: float


def exact_decimal(value: Fraction, digits: int = 18) -> str:
    """Display only; exact rational strings are always stored beside it."""
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".12E")


def interval_json(lower: Fraction, upper: Fraction) -> dict[str, object]:
    return {
        "lower": str(lower),
        "upper": str(upper),
        "lower_decimal": exact_decimal(lower),
        "upper_decimal": exact_decimal(upper),
        "relative_width": str((upper - lower) / lower),
    }


# ---------------------------------------------------------------- exact transfer operator


def build_scaled_operator(
    n: int, bonds: list[tuple[int, int]], t: Fraction
) -> ScaledOperator:
    """Return the exact integer matrix ``A = scale * R`` from e38/e39."""
    if not 0 < t < 1:
        raise ValueError("the certificate requires 0 < t < 1")
    dim = 1 << n
    q = (1 + t * t) / (2 * t)

    energies: list[int] = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    if not all((energy - parity) % 2 == 0 for energy in energies):
        raise AssertionError("bond-count parity assumption failed")
    exponents = [(energy - parity) // 2 for energy in energies]
    minimum, maximum = min(exponents), max(exponents)

    weight_scale = (
        q.denominator ** max(maximum, 0) * q.numerator ** max(-minimum, 0)
    )
    weights: list[int] = []
    for exponent in exponents:
        weight = q**exponent * weight_scale
        if weight.denominator != 1:
            raise AssertionError("weight denominator was not cleared")
        weights.append(weight.numerator)

    numerator, denominator = t.numerator, t.denominator
    numerator_powers = [numerator**power for power in range(n + 1)]
    denominator_powers = [denominator**power for power in range(n + 1)]
    p_numerator = [
        [
            numerator_powers[(row ^ col).bit_count()]
            * denominator_powers[n - (row ^ col).bit_count()]
            for col in range(dim)
        ]
        for row in range(dim)
    ]

    # Apply P on the left as n local two-by-two transforms to diag(weights) P.
    matrix = [[weights[row] * value for value in p_numerator[row]] for row in range(dim)]
    for bit in range(n):
        mask = 1 << bit
        for base in range(dim):
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

    if any(matrix[row][col] != matrix[col][row] for row in range(dim) for col in range(row)):
        raise AssertionError("the exact transfer matrix is not symmetric")
    scale = denominator ** (2 * n) * weight_scale
    return ScaledOperator(matrix, scale, q, parity, (minimum, maximum))


# ---------------------------------------------------------- exact finite group projectors


def state_permutation(n: int, site_permutation: list[int]) -> list[int]:
    if sorted(site_permutation) != list(range(n)):
        raise ValueError("invalid site permutation")
    result: list[int] = []
    for state in range(1 << n):
        image = 0
        for source, target in enumerate(site_permutation):
            if (state >> (n - 1 - source)) & 1:
                image |= 1 << (n - 1 - target)
        result.append(image)
    return result


def d4_actions() -> list[Action]:
    """All 16 actions of D4 x spin flip on the row-major 3x3 basis."""

    def rotate(row: int, col: int, power: int) -> tuple[int, int]:
        for _ in range(power):
            row, col = col, 2 - row
        return row, col

    actions: list[Action] = []
    full_mask = (1 << 9) - 1
    for reflected in (False, True):
        for power in range(4):
            sites: list[int] = []
            for row in range(3):
                for col in range(3):
                    target_row, target_col = rotate(row, col, power)
                    if reflected:
                        target_col = 2 - target_col
                    sites.append(3 * target_row + target_col)
            spatial = state_permutation(9, sites)
            for flipped in (False, True):
                mapping = [value ^ full_mask if flipped else value for value in spatial]
                prefix = "s" if reflected else "r"
                actions.append(
                    Action(
                        f"{prefix}{power}_flip{'-' if flipped else '+'}",
                        mapping,
                        power,
                        reflected,
                        flipped,
                    )
                )
    if len({tuple(action.state_map) for action in actions}) != 16:
        raise AssertionError("D4 x spin-flip action is not faithful")
    return actions


def d4_character(irrep: str, power: int, reflected: bool) -> int:
    if irrep == "A1":
        return 1
    if irrep == "A2":
        return -1 if reflected else 1
    if irrep == "B1":
        value = -1 if power % 2 else 1
        return value
    if irrep == "B2":
        value = -1 if power % 2 else 1
        return -value if reflected else value
    if irrep == "E":
        if reflected or power % 2:
            return 0
        return 2 if power == 0 else -2
    raise ValueError(f"unknown D4 irrep {irrep}")


def projector_columns(
    dimension: int, terms: list[tuple[int, list[int]]]
) -> list[dict[int, int]]:
    columns: list[dict[int, int]] = []
    for state in range(dimension):
        column: dict[int, int] = {}
        for coefficient, mapping in terms:
            image = mapping[state]
            column[image] = column.get(image, 0) + coefficient
        columns.append({image: value for image, value in column.items() if value})
    return columns


def d4_projectors(actions: list[Action]) -> list[Projector]:
    """A complete exact refinement of the ten central D4 x C2 projectors.

    The four one-dimensional D4 irreps give eight projectors.  Each E-isotypic
    projector is further split by the two eigenvalues of an axial reflection.
    This replaces two 120-dimensional pencils by four 60-dimensional pencils.
    The transfer operator acts identically on the two E rows, but both rows are
    retained and counted independently.  A common denominator 32 makes the
    complete family directly checkable.
    """
    dimensions = {"A1": 1, "A2": 1, "B1": 1, "B2": 1, "E": 2}
    projectors: list[Projector] = []
    axial_reflection = next(
        action.state_map
        for action in actions
        if action.reflected
        and action.rotation_power == 0
        and not action.flipped
    )
    for irrep, irrep_dimension in dimensions.items():
        characters = [
            d4_character(irrep, action.rotation_power, action.reflected)
            for action in actions
            if not action.flipped
        ]
        for flip_sign in (1, -1):
            central_terms: list[tuple[int, list[int]]] = []
            for action in actions:
                character = d4_character(irrep, action.rotation_power, action.reflected)
                coefficient = irrep_dimension * character
                if action.flipped:
                    coefficient *= flip_sign
                if coefficient:
                    central_terms.append((coefficient, action.state_map))

            if irrep != "E":
                terms = [(2 * coefficient, mapping) for coefficient, mapping in central_terms]
                label = f"{irrep}_flip{'+' if flip_sign == 1 else '-'}"
                projectors.append(
                    Projector(
                        label,
                        32,
                        terms,
                        projector_columns(512, terms),
                        irrep_dimension,
                        flip_sign,
                        characters,
                    )
                )
                continue

            for axial_sign in (1, -1):
                reflected_terms = [
                    (
                        axial_sign * coefficient,
                        [axial_reflection[value] for value in mapping],
                    )
                    for coefficient, mapping in central_terms
                ]
                terms = central_terms + reflected_terms
                label = (
                    f"E_axis{'+' if axial_sign == 1 else '-'}_"
                    f"flip{'+' if flip_sign == 1 else '-'}"
                )
                projectors.append(
                    Projector(
                        label,
                        32,
                        terms,
                        projector_columns(512, terms),
                        irrep_dimension,
                        flip_sign,
                        characters,
                    )
                )
    return projectors


def chain_actions(n: int) -> list[Action]:
    identity = list(range(n))
    reflection = list(reversed(identity))
    full_mask = (1 << n) - 1
    actions: list[Action] = []
    for reflected, sites in ((False, identity), (True, reflection)):
        spatial = state_permutation(n, sites)
        for flipped in (False, True):
            mapping = [value ^ full_mask if flipped else value for value in spatial]
            actions.append(
                Action(
                    f"reflection{'-' if reflected else '+'}_flip{'-' if flipped else '+'}",
                    mapping,
                    reflected=reflected,
                    flipped=flipped,
                )
            )
    return actions


def chain_projectors(n: int, actions: list[Action]) -> list[Projector]:
    projectors: list[Projector] = []
    for reflection_sign in (1, -1):
        for flip_sign in (1, -1):
            terms: list[tuple[int, list[int]]] = []
            for action in actions:
                coefficient = reflection_sign if action.reflected else 1
                if action.flipped:
                    coefficient *= flip_sign
                terms.append((coefficient, action.state_map))
            label = (
                f"reflection{'+' if reflection_sign == 1 else '-'}_"
                f"flip{'+' if flip_sign == 1 else '-'}"
            )
            projectors.append(
                Projector(
                    label,
                    4,
                    terms,
                    projector_columns(1 << n, terms),
                    1,
                    flip_sign,
                    [1, reflection_sign],
                )
            )
    return projectors


def compose_columns(
    left: list[dict[int, int]], right_column: dict[int, int]
) -> dict[int, int]:
    result: dict[int, int] = {}
    for middle, right_coefficient in right_column.items():
        for image, left_coefficient in left[middle].items():
            result[image] = result.get(image, 0) + left_coefficient * right_coefficient
    return {image: value for image, value in result.items() if value}


def matrix_invariant(matrix: list[list[int]], actions: list[Action]) -> bool:
    dimension = len(matrix)
    for action in actions:
        mapping = action.state_map
        for row in range(dimension):
            source_row = matrix[row]
            target_row = matrix[mapping[row]]
            for col in range(dimension):
                if target_row[mapping[col]] != source_row[col]:
                    return False
    return True


def verify_projectors(
    matrix: list[list[int]], actions: list[Action], projectors: list[Projector]
) -> dict[str, object]:
    """Check every projector identity on every computational-basis vector."""
    dimension = len(matrix)
    denominator = projectors[0].denominator
    same_denominator = all(projector.denominator == denominator for projector in projectors)

    self_adjoint = True
    for projector in projectors:
        for source, column in enumerate(projector.columns):
            for target, coefficient in column.items():
                if projector.columns[target].get(source, 0) != coefficient:
                    self_adjoint = False
                    break
            if not self_adjoint:
                break
        if not self_adjoint:
            break

    idempotent = True
    for projector in projectors:
        for state, column in enumerate(projector.columns):
            composed = compose_columns(projector.columns, column)
            expected = {image: denominator * value for image, value in column.items()}
            if composed != expected:
                idempotent = False
                break
        if not idempotent:
            break

    orthogonal = True
    pairs_checked = 0
    for left_index, left in enumerate(projectors):
        for right in projectors[left_index + 1 :]:
            pairs_checked += 1
            for column in right.columns:
                if compose_columns(left.columns, column):
                    orthogonal = False
                    break
            if not orthogonal:
                break
        if not orthogonal:
            break

    complete = True
    for state in range(dimension):
        total: dict[int, int] = {}
        for projector in projectors:
            for image, coefficient in projector.columns[state].items():
                total[image] = total.get(image, 0) + coefficient
        total = {image: value for image, value in total.items() if value}
        if total != {state: denominator}:
            complete = False
            break

    traces: dict[str, int] = {}
    integral_traces = True
    for projector in projectors:
        numerator_trace = sum(
            projector.columns[state].get(state, 0) for state in range(dimension)
        )
        quotient, remainder = divmod(numerator_trace, projector.denominator)
        integral_traces &= remainder == 0
        traces[projector.label] = quotient

    invariant = matrix_invariant(matrix, actions)
    return {
        "claim_tag": "[COMPUTATION]",
        "group_action_count": len(actions),
        "projector_count": len(projectors),
        "common_projector_denominator": denominator,
        "common_denominator_verified": same_denominator,
        "projector_self_adjoint_exact": self_adjoint,
        "projector_idempotence_exact": idempotent,
        "projector_orthogonality_exact": orthogonal,
        "orthogonal_pairs_checked": pairs_checked,
        "projector_completeness_exact": complete,
        "projector_traces_integral": integral_traces,
        "projector_traces": traces,
        "projector_trace_sum": sum(traces.values()),
        "matrix_invariant_under_every_group_action_exact": invariant,
        "projectors_commute_with_matrix_exact": invariant,
    }


# ------------------------------------------------------ sector bases and exact restrictions


def normalise_vector(vector: dict[int, int]) -> dict[int, int]:
    if not vector:
        return {}
    common = 0
    for value in vector.values():
        common = math.gcd(common, abs(value))
    if common > 1:
        vector = {state: value // common for state, value in vector.items()}
    first = vector[min(vector)]
    if first < 0:
        vector = {state: -value for state, value in vector.items()}
    return vector


def independent_projected_columns(
    orbit: list[int], columns: list[dict[int, int]]
) -> list[dict[int, int]]:
    """Exact rational echelon selection inside one group orbit (size at most 16)."""
    position = {state: index for index, state in enumerate(orbit)}
    echelon: dict[int, list[Fraction]] = {}
    selected: list[dict[int, int]] = []
    for state in orbit:
        candidate = normalise_vector(dict(columns[state]))
        row = [Fraction(0)] * len(orbit)
        for image, coefficient in candidate.items():
            row[position[image]] = Fraction(coefficient)
        for pivot in sorted(echelon):
            if row[pivot]:
                factor = row[pivot]
                basis_row = echelon[pivot]
                row = [left - factor * right for left, right in zip(row, basis_row)]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        pivot_value = row[pivot]
        echelon[pivot] = [value / pivot_value for value in row]
        selected.append(candidate)
    return selected


def projector_basis(
    actions: list[Action], projector: Projector
) -> list[BasisVector]:
    dimension = len(projector.columns)
    seen: set[int] = set()
    result: list[BasisVector] = []
    orbit_index = 0
    for representative in range(dimension):
        if representative in seen:
            continue
        orbit = sorted({action.state_map[representative] for action in actions})
        seen.update(orbit)
        for vector in independent_projected_columns(orbit, projector.columns):
            result.append(BasisVector(vector, orbit_index))
        orbit_index += 1
    return result


def dot_sparse(left: dict[int, int], right: dict[int, int]) -> int:
    if len(left) > len(right):
        left, right = right, left
    return sum(coefficient * right.get(state, 0) for state, coefficient in left.items())


def build_sectors(
    matrix: list[list[int]], actions: list[Action], projectors: list[Projector]
) -> list[Sector]:
    dimension = len(matrix)
    sectors: list[Sector] = []
    for projector in projectors:
        basis = projector_basis(actions, projector)
        size = len(basis)
        transformed: list[list[int]] = []
        for vector in basis:
            coefficients = vector.coefficients
            transformed.append(
                [
                    sum(matrix[row][state] * coefficient for state, coefficient in coefficients.items())
                    for row in range(dimension)
                ]
            )

        restricted = [[0] * size for _ in range(size)]
        gram = [[0] * size for _ in range(size)]
        for row in range(size):
            left = basis[row]
            for col in range(row + 1):
                right = basis[col]
                matrix_value = sum(
                    coefficient * transformed[col][state]
                    for state, coefficient in left.coefficients.items()
                )
                gram_value = (
                    0
                    if left.orbit_index != right.orbit_index
                    else dot_sparse(left.coefficients, right.coefficients)
                )
                restricted[row][col] = restricted[col][row] = matrix_value
                gram[row][col] = gram[col][row] = gram_value

        trace = sum(
            projector.columns[state].get(state, 0) for state in range(dimension)
        ) // projector.denominator
        if trace != size:
            raise AssertionError(f"basis dimension does not equal projector trace for {projector.label}")
        sectors.append(
            Sector(
                projector.label,
                restricted,
                gram,
                size,
                trace,
                projector.irrep_dimension,
                projector.flip_sign,
            )
        )
    if sum(sector.dimension for sector in sectors) != dimension:
        raise AssertionError("sector dimensions do not sum to the full dimension")
    return sectors


# --------------------------------------------------------------- exact inertia machinery


def fraction_free_inertia(lower_matrix: list[list[int]]) -> int | None:
    """Exact symmetric Bareiss/LDL inertia, with zero pivots reported."""
    size = len(lower_matrix)
    work = [row[:] for row in lower_matrix]
    previous = 1
    negatives = 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        if pivot_index + 1 == size:
            break
        for row_index in range(pivot_index + 1, size):
            row = work[row_index]
            left = row[pivot_index]
            for col_index in range(pivot_index + 1, row_index + 1):
                numerator = pivot * row[col_index] - left * work[col_index][pivot_index]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Bareiss division was not exact")
                row[col_index] = quotient
        previous = pivot
    return negatives


class SectorCounter:
    def __init__(self, sectors: list[Sector]):
        self.sectors = sectors
        self.cache: dict[Fraction, CountRecord] = {}
        self.count_seconds = 0.0
        self.max_count_seconds = 0.0

    def count(self, shift: Fraction) -> CountRecord:
        shift = Fraction(shift)
        if shift in self.cache:
            return self.cache[shift]
        started = time.perf_counter()
        by_sector: dict[str, int] = {}
        total = 0
        for sector in self.sectors:
            size = sector.dimension
            shifted = [[0] * size for _ in range(size)]
            common = 0
            for row in range(size):
                for col in range(row + 1):
                    value = shift.denominator * sector.matrix[row][col]
                    value -= shift.numerator * sector.gram[row][col]
                    shifted[row][col] = value
                    common = math.gcd(common, abs(value))
            if common > 1:
                for row in range(size):
                    for col in range(row + 1):
                        shifted[row][col] //= common
            count = fraction_free_inertia(shifted)
            if count is None:
                record = CountRecord(None, None, time.perf_counter() - started)
                break
            by_sector[sector.label] = count
            total += count
        else:
            record = CountRecord(total, by_sector, time.perf_counter() - started)
        self.count_seconds += record.seconds
        self.max_count_seconds = max(self.max_count_seconds, record.seconds)
        self.cache[shift] = record
        return record


def outward_count(
    counter: SectorCounter, shift: Fraction, direction: int
) -> tuple[Fraction, CountRecord]:
    record = counter.count(shift)
    if record.total is not None:
        return shift, record
    step = abs(shift) if shift else Fraction(1)
    for attempt in range(1, 40):
        moved = shift + direction * step * Fraction(1, 10 ** (18 - attempt // 2))
        record = counter.count(moved)
        if record.total is not None:
            return moved, record
    raise RuntimeError("no non-degenerate outward shift found")


def guide_eigenvalues(matrix: list[list[int]]) -> list[float]:
    largest = max(abs(value) for row in matrix for value in row)
    scaled = np.array([[value / largest for value in row] for row in matrix], dtype=float)
    values = np.linalg.eigvalsh(0.5 * (scaled + scaled.T))[:3]
    return [float(value * largest) for value in values]


def certified_guide_enclosure(
    counter: SectorCounter,
    index: int,
    guide: float,
    initial_radius: Fraction = Fraction(1, 1000),
) -> tuple[Fraction, Fraction, CountRecord, CountRecord]:
    """Use float only as a proposal; exact counts certify both endpoints."""
    if not math.isfinite(guide) or guide <= 0:
        raise ValueError("positive finite guide required")
    center = Fraction(format(guide, ".16e"))
    radius = initial_radius
    for _ in range(30):
        lower = max(Fraction(0), center * (1 - radius))
        upper = center * (1 + radius)
        lower, lower_record = outward_count(counter, lower, -1)
        upper, upper_record = outward_count(counter, upper, +1)
        assert lower_record.total is not None and upper_record.total is not None
        if lower_record.total <= index and upper_record.total > index:
            return lower, upper, lower_record, upper_record
        radius *= 2
    raise RuntimeError(f"could not certify an enclosure for eigenvalue {index}")


# ------------------------------------------------------- exact full-dense residual witness


def ceiling_sqrt(value: int) -> int:
    root = math.isqrt(value)
    return root if root * root == value else root + 1


def dense_dyadic_inertia_witness(
    matrix: list[list[int]],
    shifts: list[int],
    dyadic_bits: int = 36,
    include_bareiss_resource_record: bool = False,
) -> dict[str, object]:
    """Certify dense 512x512 inertias by an exact dyadic residual inequality.

    Float64 supplies only candidate eigenvectors/eigenvalues.  With
    ``Q = Q_integer / 2**dyadic_bits`` and integral diagonal ``D``, all Gram
    matrices, residuals, Frobenius bounds, and final strict inequalities below
    are recomputed over Python integers.
    """
    started = time.perf_counter()
    dimension = len(matrix)
    if dimension != 512:
        raise ValueError("the dense witness is configured for the 512-dimensional layer")
    if max(abs(value) for row in matrix for value in row) >= 1 << 53:
        raise AssertionError("float proposal matrix entries are not exactly representable")
    proposed = np.array(matrix, dtype=float)
    if any(int(proposed[row, col]) != matrix[row][col] for row in range(dimension) for col in range(dimension)):
        raise AssertionError("integer-to-float proposal conversion was not exact")

    eigenvalues, eigenvectors = np.linalg.eigh(proposed)
    dyadic_denominator = 1 << dyadic_bits
    q_integer = (
        np.rint(eigenvectors * dyadic_denominator)
        .astype(np.int64)
        .astype(object)
    )
    diagonal = np.array([int(round(float(value))) for value in eigenvalues], dtype=object)

    gram_started = time.perf_counter()
    gram = q_integer.T @ q_integer
    gram_seconds = time.perf_counter() - gram_started
    product_started = time.perf_counter()
    model_numerator = (q_integer * diagonal) @ q_integer.T
    model_seconds = time.perf_counter() - product_started

    square_denominator = dyadic_denominator * dyadic_denominator
    gram_square_sum = 0
    for row in range(dimension):
        for col in range(dimension):
            difference = int(gram[row, col]) - (square_denominator if row == col else 0)
            gram_square_sum += difference * difference
    gram_frobenius_numerator = ceiling_sqrt(gram_square_sum)
    invertibility_margin_numerator = square_denominator - gram_frobenius_numerator
    if invertibility_margin_numerator <= 0:
        raise AssertionError("rounded dyadic eigenbasis was not certified invertible")

    witnesses: list[dict[str, object]] = []
    for shift in shifts:
        residual_square_sum = 0
        for row in range(dimension):
            for col in range(dimension):
                exact_entry = matrix[row][col] * square_denominator
                if row == col:
                    exact_entry -= shift * square_denominator
                model_entry = int(model_numerator[row, col]) - shift * int(gram[row, col])
                residual = exact_entry - model_entry
                residual_square_sum += residual * residual
        residual_numerator = ceiling_sqrt(residual_square_sum)
        minimum_diagonal_distance = min(abs(int(value) - shift) for value in diagonal)
        safety_margin = (
            invertibility_margin_numerator * minimum_diagonal_distance
            - residual_numerator
        )
        certified = safety_margin > 0
        witnesses.append(
            {
                "claim_tag": "[COMPUTATION]",
                "shift_integer_A_units": str(shift),
                "inertia_count_below_shift": sum(int(value) < shift for value in diagonal),
                "minimum_abs_diagonal_D_minus_shift": str(minimum_diagonal_distance),
                "residual_frobenius_upper_bound": str(
                    Fraction(residual_numerator, square_denominator)
                ),
                "sigma_min_model_lower_bound": str(
                    Fraction(
                        invertibility_margin_numerator * minimum_diagonal_distance,
                        square_denominator,
                    )
                ),
                "strict_safety_margin_numerator": str(safety_margin),
                "inertia_certified_exactly": certified,
            }
        )

    result: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "method": (
            "Exact integer/dyadic a-posteriori congruence certificate on the full dense "
            "512x512 matrix; float64 proposes Q,D only."
        ),
        "dimension": dimension,
        "dyadic_basis_denominator": str(dyadic_denominator),
        "dyadic_square_denominator": str(square_denominator),
        "matrix_entries_exactly_representable_in_float_proposal": True,
        "gram_frobenius_upper_bound": str(
            Fraction(gram_frobenius_numerator, square_denominator)
        ),
        "rounded_basis_invertibility_margin": str(
            Fraction(invertibility_margin_numerator, square_denominator)
        ),
        "witnesses": witnesses,
        "elapsed": {
            "exact_gram_product_seconds": round(gram_seconds, 6),
            "exact_model_product_seconds": round(model_seconds, 6),
            "total_seconds": round(time.perf_counter() - started, 6),
        },
    }
    if include_bareiss_resource_record:
        result["bareiss_resource_record"] = {
            "claim_tag": "[UNRESOLVED]",
            "algorithm": "full 512x512 fraction-free Bareiss LDL at integer shift 3000000000",
            "observed_wall_seconds_lower_bound": 1200,
            "outcome": "timed out before completion",
            "replacement": "the exact dyadic residual witness above certifies the same full-dense inertia",
        }
    return result


# ----------------------------------------------------------------------- case driver


def spectrum_endpoint_record(record: CountRecord) -> dict[str, object]:
    if record.total is None or record.by_sector is None:
        raise AssertionError("a degenerate count cannot be stored as a certificate")
    return {"total": record.total, "by_sector": record.by_sector}


def certify_case(
    name: str,
    kind: str,
    bonds: list[tuple[int, int]],
    t: Fraction,
) -> tuple[dict[str, object], ScaledOperator, list[Sector], SectorCounter]:
    started = time.perf_counter()
    timings: dict[str, float] = {}

    began = time.perf_counter()
    scaled = build_scaled_operator(9, bonds, t)
    timings["matrix_build_seconds"] = time.perf_counter() - began

    if kind == "3D-layer":
        actions = d4_actions()
        projectors = d4_projectors(actions)
        symmetry_name = "D4 x C2(global spin flip)"
    else:
        actions = chain_actions(9)
        projectors = chain_projectors(9, actions)
        symmetry_name = "C2(chain reflection) x C2(global spin flip)"

    began = time.perf_counter()
    verification = verify_projectors(scaled.matrix, actions, projectors)
    timings["projector_verification_seconds"] = time.perf_counter() - began
    required_projector_flags = (
        "common_denominator_verified",
        "projector_self_adjoint_exact",
        "projector_idempotence_exact",
        "projector_orthogonality_exact",
        "projector_completeness_exact",
        "projector_traces_integral",
        "matrix_invariant_under_every_group_action_exact",
        "projectors_commute_with_matrix_exact",
    )
    if not all(verification[key] for key in required_projector_flags):
        raise AssertionError(f"an exact symmetry check failed for {name}")

    began = time.perf_counter()
    sectors = build_sectors(scaled.matrix, actions, projectors)
    timings["sector_build_seconds"] = time.perf_counter() - began
    counter = SectorCounter(sectors)
    positive = counter.count(Fraction(0))
    if positive.total != 0:
        raise AssertionError("the exact transfer operator is not positive definite")

    guides = guide_eigenvalues(scaled.matrix)
    raw_enclosures: list[tuple[Fraction, Fraction]] = []
    endpoint_records: list[tuple[CountRecord, CountRecord]] = []
    for index, guide in enumerate(guides):
        lower, upper, lower_record, upper_record = certified_guide_enclosure(
            counter, index, guide
        )
        raw_enclosures.append((lower, upper))
        endpoint_records.append((lower_record, upper_record))

    (l0, h0), (l1, h1), (l2, h2) = raw_enclosures
    distinct = h0 < l1 and h1 < l2
    if not distinct:
        raise AssertionError("the three lowest exact enclosures are not distinct")
    predicted_lower = l1 * l2 / h0
    predicted_upper = h1 * h2 / l0
    tested_lower, lower_window_record = outward_count(counter, predicted_lower, -1)
    tested_upper, upper_window_record = outward_count(counter, predicted_upper, 1)
    assert lower_window_record.total is not None and upper_window_record.total is not None
    absent = lower_window_record.total == upper_window_record.total

    lambda_fields: dict[str, object] = {}
    for index, ((lower, upper), (lower_record, upper_record)) in enumerate(
        zip(raw_enclosures, endpoint_records)
    ):
        normalised_lower, normalised_upper = lower / scaled.scale, upper / scaled.scale
        lower_by_sector = lower_record.by_sector or {}
        upper_by_sector = upper_record.by_sector or {}
        lambda_fields[f"lambda{index}"] = {
            **interval_json(normalised_lower, normalised_upper),
            "endpoint_inertia": {
                "lower": spectrum_endpoint_record(lower_record),
                "upper": spectrum_endpoint_record(upper_record),
            },
            "sector_multiplicity_in_enclosure": {
                label: upper_by_sector[label] - lower_by_sector[label]
                for label in lower_by_sector
                if upper_by_sector[label] != lower_by_sector[label]
            },
        }

    normalized_prediction = (
        predicted_lower / scaled.scale,
        predicted_upper / scaled.scale,
    )
    normalized_tested = (tested_lower / scaled.scale, tested_upper / scaled.scale)
    timings["inertia_seconds"] = counter.count_seconds
    timings["max_single_inertia_seconds"] = counter.max_count_seconds
    timings["total_seconds"] = time.perf_counter() - started

    theorem_statement = (
        f"The exact open 3x3 layer transfer operator at t={t} is not a positive-spectrum "
        "nine-mode Gaussian subset-product operator: the three lowest eigenvalues are "
        "distinct and the forced-value interval is empty."
        if kind == "3D-layer"
        else f"For the exact open n=9 chain control at t={t}, the certified forced-value "
        "window is occupied; this is a consistency control, not an equality proof."
    )
    row = {
        "name": name,
        "kind": kind,
        "claim_tag": "[THEOREM]" if kind == "3D-layer" else "[COMPUTATION]",
        "statement": theorem_statement,
        "status": "CERTIFIED",
        "t_tanh_half_Kstar": str(t),
        "exp_2K": str(scaled.q),
        "n_sites": 9,
        "n_bonds": len(bonds),
        "dimension": 512,
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
            "predicted": interval_json(*normalized_prediction),
            "actually_tested": interval_json(*normalized_tested),
            "tested_window_encloses_prediction": (
                tested_lower <= predicted_lower <= predicted_upper <= tested_upper
            ),
            "endpoint_inertia": {
                "lower": spectrum_endpoint_record(lower_window_record),
                "upper": spectrum_endpoint_record(upper_window_record),
            },
        },
        "inertia_counts_at_forced_window_ends": [
            lower_window_record.total,
            upper_window_record.total,
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
                "exact independent integer columns of the refined projector numerator; "
                "B=C^T A C and G=C^T C"
            ),
        },
        "exact_inertia_shifts": len(counter.cache),
        "elapsed": {key: round(value, 6) for key, value in timings.items()},
    }
    return row, scaled, sectors, counter


def fraction_floor(value: Fraction) -> int:
    return value.numerator // value.denominator


def fraction_ceiling(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def certify_control_case(
    name: str,
    bonds: list[tuple[int, int]],
    t: Fraction,
) -> tuple[dict[str, object], ScaledOperator, list[Sector], dict[str, object]]:
    """Certify the n=9 chain control with one exact full-dense dyadic witness.

    A single rounded rational basis certifies all nine shifts at once.  This is
    substantially faster than repeated 136-dimensional Bareiss pencils while
    retaining an integer-only decision certificate.
    """
    started = time.perf_counter()
    timings: dict[str, float] = {}
    began = time.perf_counter()
    scaled = build_scaled_operator(9, bonds, t)
    timings["matrix_build_seconds"] = time.perf_counter() - began

    actions = chain_actions(9)
    projectors = chain_projectors(9, actions)
    began = time.perf_counter()
    verification = verify_projectors(scaled.matrix, actions, projectors)
    timings["projector_verification_seconds"] = time.perf_counter() - began
    required_flags = (
        "projector_self_adjoint_exact",
        "projector_idempotence_exact",
        "projector_orthogonality_exact",
        "projector_completeness_exact",
        "projectors_commute_with_matrix_exact",
    )
    if not all(verification[key] for key in required_flags):
        raise AssertionError(f"an exact chain-projector check failed for {name}")

    began = time.perf_counter()
    sectors = build_sectors(scaled.matrix, actions, projectors)
    timings["sector_build_seconds"] = time.perf_counter() - began

    guides = guide_eigenvalues(scaled.matrix)
    radius = Fraction(1, 100)
    raw_enclosures: list[tuple[Fraction, Fraction]] = []
    for guide in guides:
        center = Fraction(format(guide, ".16e"))
        lower = Fraction(fraction_floor(center * (1 - radius)))
        upper = Fraction(fraction_ceiling(center * (1 + radius)))
        raw_enclosures.append((lower, upper))

    (l0, h0), (l1, h1), (l2, h2) = raw_enclosures
    if not (h0 < l1 and h1 < l2):
        raise AssertionError("the proposed chain enclosures are not pairwise disjoint")
    predicted_lower = l1 * l2 / h0
    predicted_upper = h1 * h2 / l0
    tested_lower = Fraction(fraction_floor(predicted_lower))
    tested_upper = Fraction(fraction_ceiling(predicted_upper))
    requested_shifts = [
        0,
        *[
            int(endpoint)
            for enclosure in raw_enclosures
            for endpoint in enclosure
        ],
        int(tested_lower),
        int(tested_upper),
    ]
    shifts = list(dict.fromkeys(requested_shifts))
    dense = dense_dyadic_inertia_witness(
        scaled.matrix,
        shifts,
        dyadic_bits=40,
    )
    if not all(witness["inertia_certified_exactly"] for witness in dense["witnesses"]):
        raise AssertionError("an exact dense chain inertia witness failed its residual bound")
    witness_by_shift = {
        int(witness["shift_integer_A_units"]): witness
        for witness in dense["witnesses"]
    }

    def record_at(shift: Fraction) -> CountRecord:
        count = int(witness_by_shift[int(shift)]["inertia_count_below_shift"])
        return CountRecord(count, {"full_dense": count}, 0.0)

    positive = record_at(Fraction(0))
    if positive.total != 0:
        raise AssertionError("the chain transfer operator is not positive definite")

    endpoint_records = [
        (record_at(lower), record_at(upper))
        for lower, upper in raw_enclosures
    ]
    for index, (lower_record, upper_record) in enumerate(endpoint_records):
        assert lower_record.total is not None and upper_record.total is not None
        if not (lower_record.total <= index and upper_record.total > index):
            raise AssertionError(f"dense witness does not bracket chain eigenvalue {index}")

    lower_window_record = record_at(tested_lower)
    upper_window_record = record_at(tested_upper)
    assert lower_window_record.total is not None and upper_window_record.total is not None
    absent = lower_window_record.total == upper_window_record.total

    lambda_fields: dict[str, object] = {}
    for index, ((lower, upper), (lower_record, upper_record)) in enumerate(
        zip(raw_enclosures, endpoint_records)
    ):
        lambda_fields[f"lambda{index}"] = {
            **interval_json(lower / scaled.scale, upper / scaled.scale),
            "endpoint_inertia": {
                "lower": spectrum_endpoint_record(lower_record),
                "upper": spectrum_endpoint_record(upper_record),
            },
            "sector_multiplicity_in_enclosure": {
                "full_dense": upper_record.total - lower_record.total
            },
        }

    timings["exact_dense_witness_seconds"] = dense["elapsed"]["total_seconds"]
    timings["total_seconds"] = time.perf_counter() - started
    row = {
        "name": name,
        "kind": "1D-control",
        "claim_tag": "[COMPUTATION]",
        "statement": (
            f"For the exact open n=9 chain control at t={t}, the certified "
            "forced-value window is occupied; this is a consistency control, "
            "not an equality proof."
        ),
        "status": "CERTIFIED",
        "t_tanh_half_Kstar": str(t),
        "exp_2K": str(scaled.q),
        "n_sites": 9,
        "n_bonds": len(bonds),
        "dimension": 512,
        "integer_scale_for_R": str(scaled.scale),
        "bond_parity": scaled.bond_parity,
        "energy_exponent_range": list(scaled.energy_exponents),
        "positive_definite_count_below_zero": positive.total,
        "lambda_enclosures": lambda_fields,
        "three_lowest_distinct": True,
        "relative_gap_lower_bounds": {
            "gap01": str((l1 - h0) / l1),
            "gap12": str((l2 - h1) / l2),
        },
        "forced_value_window": {
            "predicted": interval_json(
                predicted_lower / scaled.scale,
                predicted_upper / scaled.scale,
            ),
            "actually_tested": interval_json(
                tested_lower / scaled.scale,
                tested_upper / scaled.scale,
            ),
            "tested_window_encloses_prediction": (
                tested_lower <= predicted_lower <= predicted_upper <= tested_upper
            ),
            "endpoint_inertia": {
                "lower": spectrum_endpoint_record(lower_window_record),
                "upper": spectrum_endpoint_record(upper_window_record),
            },
        },
        "inertia_counts_at_forced_window_ends": [
            lower_window_record.total,
            upper_window_record.total,
        ],
        "predicted_eigenvalue_absent": absent,
        "symmetry": {
            "claim_tag": "[COMPUTATION]",
            "group": "C2(chain reflection) x C2(global spin flip)",
            "verification": verification,
            "sector_dimensions": {
                sector.label: sector.dimension for sector in sectors
            },
            "sector_dimensions_sum": sum(sector.dimension for sector in sectors),
            "sector_irrep_dimensions": {
                sector.label: sector.irrep_dimension for sector in sectors
            },
            "basis_convention": (
                "exact integer projector bases were constructed; chain spectral "
                "counts use the independently certified full-dense dyadic witness"
            ),
        },
        "full_dense_exact_witness": dense,
        "exact_inertia_shifts": len(shifts),
        "elapsed": {key: round(value, 6) for key, value in timings.items()},
    }
    return row, scaled, sectors, dense


def check_record(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    grid_bonds = list(layer_bonds((3, 3), (False, False)))
    chain_bonds = list(layer_bonds((9,), (False,)))
    rows: list[dict[str, object]] = []
    retained: dict[str, tuple[ScaledOperator, list[Sector], object]] = {}

    definitions = [
        ("layer_3x3_t_1_3", "3D-layer", grid_bonds, Fraction(1, 3)),
        ("layer_3x3_t_1_2", "3D-layer", grid_bonds, Fraction(1, 2)),
        ("control_chain_n9_t_1_3", "1D-control", chain_bonds, Fraction(1, 3)),
        ("control_chain_n9_t_1_2", "1D-control", chain_bonds, Fraction(1, 2)),
    ]
    for name, kind, bonds, t in definitions:
        print(f"certifying {name}", flush=True)
        if kind == "3D-layer":
            row, scaled, sectors, certificate = certify_case(name, kind, bonds, t)
        else:
            row, scaled, sectors, certificate = certify_control_case(name, bonds, t)
        rows.append(row)
        retained[name] = (scaled, sectors, certificate)
        counts = row["inertia_counts_at_forced_window_ends"]
        print(
            f"  exact counts={counts}; absent={row['predicted_eigenvalue_absent']}; "
            f"seconds={row['elapsed']['total_seconds']}",
            flush=True,
        )

    # A deliberately coarse, simple-integer window contains the t=1/2 prediction and
    # lies in its large spectral gap.  It makes the independent full-dense check compact.
    dense_name = "layer_3x3_t_1_2"
    dense_row = next(row for row in rows if row["name"] == dense_name)
    dense_scaled, _, dense_counter = retained[dense_name]
    if not isinstance(dense_counter, SectorCounter):
        raise AssertionError("the retained layer certificate is not a sector counter")
    dense_shifts = [3_000_000_000, 3_100_000_000]
    predicted = dense_row["forced_value_window"]["predicted"]
    predicted_raw_lower = Fraction(predicted["lower"]) * dense_scaled.scale
    predicted_raw_upper = Fraction(predicted["upper"]) * dense_scaled.scale
    if not dense_shifts[0] <= predicted_raw_lower <= predicted_raw_upper <= dense_shifts[1]:
        raise AssertionError("the configured dense verification window misses the prediction")
    dense_sector_records = [dense_counter.count(Fraction(shift)) for shift in dense_shifts]
    if any(record.total is None for record in dense_sector_records):
        raise AssertionError("a dense verification endpoint is a degenerate sector pivot")
    print("building exact full-dense dyadic residual witness", flush=True)
    dense_verification = dense_dyadic_inertia_witness(
        dense_scaled.matrix,
        dense_shifts,
        include_bareiss_resource_record=True,
    )
    dense_verification["coarse_window_contains_predicted_window"] = True
    dense_verification["sector_counts"] = [
        spectrum_endpoint_record(record) for record in dense_sector_records
    ]
    dense_verification["dense_counts_match_sector_sums"] = all(
        witness["inertia_count_below_shift"] == record.total
        for witness, record in zip(dense_verification["witnesses"], dense_sector_records)
    )

    checks: list[dict[str, object]] = []
    layer_rows = [row for row in rows if row["kind"] == "3D-layer"]
    control_rows = [row for row in rows if row["kind"] == "1D-control"]
    checks.append(
        check_record(
            "two_exact_3x3_absence_certificates",
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
        )
    )
    checks.append(
        check_record(
            "n9_chain_controls_occupied_at_both_couplings",
            len(control_rows) == 2
            and all(
                row["three_lowest_distinct"]
                and not row["predicted_eigenvalue_absent"]
                and len(set(row["inertia_counts_at_forced_window_ends"])) == 2
                for row in control_rows
            ),
            "; ".join(
                f"{row['name']} counts={row['inertia_counts_at_forced_window_ends']}"
                for row in control_rows
            ),
        )
    )
    checks.append(
        check_record(
            "all_exact_projector_checks",
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
            "all four cases passed self-adjointness, idempotence, orthogonality, completeness, and commutation",
        )
    )
    checks.append(
        check_record(
            "sector_dimensions_sum_to_512",
            all(row["symmetry"]["sector_dimensions_sum"] == 512 for row in rows),
            "; ".join(
                f"{row['name']} sum={row['symmetry']['sector_dimensions_sum']}"
                for row in rows
            ),
        )
    )
    checks.append(
        check_record(
            "degree_four_vertex_present",
            max(sum(vertex in bond for bond in grid_bonds) for vertex in range(9)) == 4,
            "the centre site of the open 3x3 grid has degree 4",
        )
    )
    checks.append(
        check_record(
            "full_dense_exact_inertia_matches_sector_bookkeeping",
            dense_verification["dense_counts_match_sector_sums"]
            and all(
                witness["inertia_certified_exactly"]
                for witness in dense_verification["witnesses"]
            ),
            (
                f"dense counts={[w['inertia_count_below_shift'] for w in dense_verification['witnesses']]}; "
                f"sector sums={[record.total for record in dense_sector_records]}"
            ),
        )
    )

    result = {
        "provenance": {
            "script": "experiments/e54_spectral_3x3.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python_version": platform.python_version(),
            "method": (
                "Exact integer transfer matrices; exact D4 x spin-flip central projectors; "
                "fraction-free Bareiss LDL in generalized sector pencils; exact rational "
                "forced-value interval arithmetic; exact full-dense dyadic residual inertia "
                "witness. Float64 proposes brackets and dyadic bases only and never decides a claim."
            ),
        },
        "data": {
            "claim_tag": "[THEOREM]",
            "scope": (
                "Exactly the open 3x3 layer at t in {1/3,1/2}; finite 512-dimensional "
                "positive-spectrum Gaussian subset-product obstruction only."
            ),
            "parameterization": {
                "claim_tag": "[LEMMA]",
                "t": "tanh(K*/2)",
                "exp_2K": "(1+t^2)/(2t)",
                "operator": "R=P_t diag(q^((b-eps)/2)) P_t, P_t[k,l]=t^hamming(k,l)",
            },
            "rows": rows,
            "full_dense_verification": dense_verification,
            "limitations": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "No interval in coupling, all-size theorem, thermodynamic-limit statement, "
                    "larger-system embedding obstruction, or parity-projected escape-route "
                    "classification is proved here. Occupied chain windows do not prove exact "
                    "algebraic equality with the forced value."
                ),
            },
        },
        "checks": checks,
    }

    output = ROOT / "results" / "spectral" / "spectral_3x3.json"
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
