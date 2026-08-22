"""Standalone independent audit for the e152 low-temperature extension (x^54).

This file deliberately does not import ``experiments/e124_lt_x34.py`` or
``ising.series.low_temperature_free_energy``.  It rebuilds all finite-box
inputs through x^54: an exact Python state-DP handles every canonical box
whose transfer cross-section has at most twelve spins, while a fresh-prime
modular transfer implementation reconstructs the remaining CRT boxes.  It
then uses formal inversion/integration and mixed second differences, rather
than the producer's recursive box-weight subtraction, to obtain the bulk
series.
"""

from __future__ import annotations

import gc
import hashlib
import json
import math
from fractions import Fraction
from itertools import combinations_with_replacement, product
from pathlib import Path
from typing import Iterable

import numpy as np
from sympy import Matrix, Rational

from ising.transfer_matrix import box_broken_bond_poly


ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "series" / "lt_x54.json"
PREDECESSOR_PATH = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
STRUCTURE_PATH = ROOT / "results" / "series" / "structure_certificates.json"
PROOF_PATH = ROOT / "proofs" / "lt_x54.md"
TARGET_ORDER = 54
TARGET_BUDGET = (TARGET_ORDER + 6) // 4
PREFIX_LENGTH = TARGET_ORDER + 1
DEFAULT_PRODUCER_PRIMES = {
    2147483647,
    2147483629,
    2147483587,
    2147483579,
}


def _shape_key(shape: tuple[int, int, int]) -> str:
    return "x".join(str(value) for value in shape)


def _canonical_boxes(budget: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        shape
        for shape in combinations_with_replacement(range(1, budget + 1), 3)
        if sum(shape) <= budget
    )


def _ordered_boxes(budget: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (a, b, c)
        for a in range(1, budget + 1)
        for b in range(1, budget + 1)
        for c in range(1, budget + 1)
        if a + b + c <= budget
    )


def _layer_bonds(shape: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    rows, columns = shape
    bonds: list[tuple[int, int]] = []
    for row in range(rows):
        for column in range(columns):
            index = row * columns + column
            if column + 1 < columns:
                bonds.append((index, index + 1))
            if row + 1 < rows:
                bonds.append((index, index + columns))
    return tuple(bonds)


def _layer_data(shape: tuple[int, int]) -> tuple[list[int], list[int]]:
    """Return physical-plus-inplane-ghost and transfer-ghost counts per state."""
    cross_sites = math.prod(shape)
    bonds = _layer_bonds(shape)
    degree = [0] * cross_sites
    for left, right in bonds:
        degree[left] += 1
        degree[right] += 1
    inplane_ghosts = [4 - value for value in degree]
    state_count = 1 << cross_sites
    inplane = [0] * state_count
    transfer_ghost = [0] * state_count
    for state in range(state_count):
        broken = sum(
            ((state >> left) ^ (state >> right)) & 1 for left, right in bonds
        )
        down = 0
        for site, ghost_count in enumerate(inplane_ghosts):
            if not ((state >> site) & 1):
                broken += ghost_count
                down += 1
        inplane[state] = broken
        transfer_ghost[state] = down
    return inplane, transfer_ghost


def _independent_plus_polynomial(
    shape: tuple[int, int, int], order: int
) -> list[int]:
    """Exact Python-int layer DP for a frozen-plus box, truncated at ``order``.

    The implementation is intentionally independent of the producer's NumPy
    transfer engine.  It stores one Python polynomial for each layer state and
    applies each vertical bond factor by explicit state-pair butterflies.
    """
    first, second, layers = shape
    cross_sites = first * second
    state_count = 1 << cross_sites
    width = order + 1
    inplane, transfer_ghost = _layer_data((first, second))

    vectors = [[0] * width for _ in range(state_count)]
    for state in range(state_count):
        degree = inplane[state] + transfer_ghost[state]
        if layers == 1:
            degree += transfer_ghost[state]
        if degree <= order:
            vectors[state][degree] = 1

    for layer in range(1, layers):
        for bit_index in range(cross_sites):
            block = 1 << bit_index
            period = block << 1
            for base in range(0, state_count, period):
                for offset in range(block):
                    left = base + offset
                    right = left + block
                    left_poly = vectors[left]
                    right_poly = vectors[right]
                    vectors[left] = [
                        left_poly[degree]
                        + (right_poly[degree - 1] if degree else 0)
                        for degree in range(width)
                    ]
                    vectors[right] = [
                        right_poly[degree]
                        + (left_poly[degree - 1] if degree else 0)
                        for degree in range(width)
                    ]

        on_last_layer = layer == layers - 1
        for state, polynomial in enumerate(vectors):
            shift = inplane[state] + (
                transfer_ghost[state] if on_last_layer else 0
            )
            vectors[state] = (
                [0] * min(shift, width) + polynomial[: max(0, width - shift)]
            )

    result = [0] * width
    for polynomial in vectors:
        for degree, value in enumerate(polynomial):
            result[degree] += value
    return result


def _brute_plus_polynomial(
    shape: tuple[int, int, int], order: int
) -> list[int]:
    """Independent spin-configuration count including all frozen exterior bonds."""
    first, second, layers = shape
    cross_sites = first * second
    site_count = cross_sites * layers
    bonds = _layer_bonds((first, second))
    degree = [0] * cross_sites
    for left, right in bonds:
        degree[left] += 1
        degree[right] += 1
    inplane_ghosts = [4 - value for value in degree]
    result = [0] * (order + 1)

    for configuration in range(1 << site_count):
        broken = 0
        for layer in range(layers):
            base = layer * cross_sites
            for left, right in bonds:
                broken += ((configuration >> (base + left)) ^ (configuration >> (base + right))) & 1
        for layer in range(layers - 1):
            base = layer * cross_sites
            next_base = base + cross_sites
            for site in range(cross_sites):
                broken += ((configuration >> (base + site)) ^ (configuration >> (next_base + site))) & 1
        for layer in range(layers):
            base = layer * cross_sites
            for site, ghost_count in enumerate(inplane_ghosts):
                if not ((configuration >> (base + site)) & 1):
                    broken += ghost_count
                    if layer == 0:
                        broken += 1
                    if layer == layers - 1:
                        broken += 1
        if broken <= order:
            result[broken] += 1
    return result


def _is_prime_32(value: int) -> bool:
    if value < 2:
        return False
    for divisor in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if value == divisor:
            return True
        if value % divisor == 0:
            return False
    odd_part = value - 1
    exponent = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        exponent += 1
    for base in (2, 3, 5, 7, 11, 13, 17):
        if base >= value:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(exponent - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def _fresh_primes_for_sites(site_count: int) -> tuple[int, ...]:
    """Pick deterministic 31-bit primes disjoint from the producer's moduli."""
    target = 1 << (site_count + 1)  # strictly larger than twice 2**site_count
    primes: list[int] = []
    product_value = 1
    candidate = 2147483501
    if candidate % 2 == 0:
        candidate -= 1
    while product_value <= target:
        if candidate not in DEFAULT_PRODUCER_PRIMES and _is_prime_32(candidate):
            primes.append(candidate)
            product_value *= candidate
        candidate -= 2
    if set(primes) & DEFAULT_PRODUCER_PRIMES:
        raise AssertionError("fresh CRT moduli overlap the producer's moduli")
    return tuple(primes)


def _independent_modular_residue(
    shape: tuple[int, int, int], order: int, modulus: int
) -> list[int]:
    """Fresh NumPy tensor-layout propagation modulo ``modulus``.

    This does not import the producer's CRT helper or any private transfer
    primitive.  State-dependent shifts are grouped by their explicitly
    constructed broken-bond values; every vertical-bond butterfly is reduced
    after the bit update.
    """
    first, second, layers = shape
    cross_sites = first * second
    state_count = 1 << cross_sites
    width = order + 1
    bonds = _layer_bonds((first, second))
    degree = [0] * cross_sites
    for left, right in bonds:
        degree[left] += 1
        degree[right] += 1
    inplane_ghosts = [4 - value for value in degree]

    states = np.arange(state_count, dtype=np.uint32)
    inplane = np.zeros(state_count, dtype=np.int64)
    transfer_ghost = np.zeros(state_count, dtype=np.int64)
    for left, right in bonds:
        inplane += (
            ((states >> np.uint32(left)) ^ (states >> np.uint32(right)))
            & np.uint32(1)
        ).astype(np.int64)
    for site, ghost_count in enumerate(inplane_ghosts):
        down = (((states >> np.uint32(site)) & np.uint32(1)) == 0).astype(np.int64)
        inplane += ghost_count * down
        transfer_ghost += down

    def shift_rows(vector: np.ndarray, shifts: np.ndarray) -> np.ndarray:
        shifted = np.zeros_like(vector)
        for shift in np.unique(shifts):
            offset = int(shift)
            if offset > order:
                continue
            mask = shifts == offset
            shifted[mask, offset:] = vector[mask, : width - offset]
        return shifted

    vector = np.zeros((state_count, width), dtype=np.int64)
    first_shift = inplane + transfer_ghost
    if layers == 1:
        first_shift = first_shift + transfer_ghost
    valid = first_shift <= order
    vector[np.nonzero(valid)[0], first_shift[valid]] = 1

    for layer in range(1, layers):
        for bit_index in range(cross_sites):
            block = 1 << bit_index
            groups = state_count // (2 * block)
            tensor = vector.reshape(groups, 2, block, width)
            updated = tensor.copy()
            updated[:, 0, :, 1:] += tensor[:, 1, :, :-1]
            updated[:, 1, :, 1:] += tensor[:, 0, :, :-1]
            np.remainder(updated, modulus, out=updated)
            vector = updated.reshape(state_count, width)
        layer_shift = inplane + (transfer_ghost if layer == layers - 1 else 0)
        vector = shift_rows(vector, layer_shift)
        np.remainder(vector, modulus, out=vector)

    result = vector.sum(axis=0, dtype=np.int64)
    np.remainder(result, modulus, out=result)
    output = [int(value) for value in result]
    del states, inplane, transfer_ghost, vector, result
    gc.collect()
    return output


def _independent_crt_prefix(
    shape: tuple[int, int, int], order: int
) -> tuple[list[int], dict]:
    site_count = math.prod(shape)
    primes = _fresh_primes_for_sites(site_count)
    reconstructed: list[int] | None = None
    product_value = 1
    for prime in primes:
        residues = _independent_modular_residue(shape, order, prime)
        if reconstructed is None:
            reconstructed = residues
        else:
            inverse = pow(product_value % prime, -1, prime)
            for index, residue in enumerate(residues):
                correction = ((residue - reconstructed[index]) % prime) * inverse % prime
                reconstructed[index] += product_value * correction
        product_value *= prime
    if reconstructed is None:
        raise AssertionError("CRT selected no moduli")
    configuration_bound = 1 << site_count
    if product_value <= 2 * configuration_bound:
        raise AssertionError("fresh CRT modulus does not exceed twice the coefficient bound")
    if any(value < 0 or value > configuration_bound for value in reconstructed):
        raise AssertionError("reconstructed coefficient violates the configuration-count bound")
    return reconstructed, {
        "primes": list(primes),
        "modulus_product": product_value,
        "coefficient_bound": configuration_bound,
    }


def _formal_log_via_inverse(polynomial: list[int]) -> tuple[Fraction, ...]:
    """Compute log(A) from A'/A, independently of the producer's recurrence."""
    if not polynomial or polynomial[0] != 1:
        raise AssertionError("finite-box polynomial must have unit constant term")
    order = len(polynomial) - 1
    values = [Fraction(value) for value in polynomial]
    inverse = [Fraction(0)] * (order + 1)
    inverse[0] = Fraction(1)
    for degree in range(1, order + 1):
        inverse[degree] = -sum(
            values[offset] * inverse[degree - offset]
            for offset in range(1, degree + 1)
        )
    derivative = [
        Fraction(degree + 1) * values[degree + 1]
        if degree < order
        else Fraction(0)
        for degree in range(order + 1)
    ]
    logarithmic_derivative = [Fraction(0)] * order
    for degree in range(order):
        logarithmic_derivative[degree] = sum(
            derivative[offset] * inverse[degree - offset]
            for offset in range(degree + 1)
        )
    result = [Fraction(0)] * (order + 1)
    for degree in range(1, order + 1):
        result[degree] = logarithmic_derivative[degree - 1] / degree
    return tuple(result)


def _mixed_second_difference(
    shape: tuple[int, int, int],
    logs: dict[tuple[int, int, int], tuple[Fraction, ...]],
) -> tuple[Fraction, ...]:
    """Compute W_A = Delta_1^2 Delta_2^2 Delta_3^2 log Xi_A exactly."""
    order = len(next(iter(logs.values()))) - 1
    result = [Fraction(0)] * (order + 1)
    coefficients = (1, -2, 1)
    for deltas in product(range(3), repeat=3):
        reduced = tuple(shape[index] - deltas[index] for index in range(3))
        if min(reduced) <= 0:
            continue
        factor = math.prod(coefficients[delta] for delta in deltas)
        logarithm = logs[tuple(sorted(reduced))]
        for degree, value in enumerate(logarithm):
            result[degree] += factor * value
    return tuple(result)


def _primitive_vector(vector: Iterable) -> list[int]:
    fractions = [Fraction(int(value.p), int(value.q)) for value in vector]
    denominator = 1
    for value in fractions:
        denominator = math.lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in fractions]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor:
        integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return integers


def _to_sympy_rows(rows: list[list[Fraction]]) -> Matrix:
    return Matrix(
        [[Rational(value.numerator, value.denominator) for value in row] for row in rows]
    )


def _convolve(
    left: list[Fraction] | tuple[Fraction, ...],
    right: list[Fraction] | tuple[Fraction, ...],
    order: int,
) -> list[Fraction]:
    result = [Fraction(0)] * (order + 1)
    for first_degree, first_value in enumerate(left[: order + 1]):
        if not first_value:
            continue
        for second_degree, second_value in enumerate(right[: order + 1 - first_degree]):
            if second_value:
                result[first_degree + second_degree] += first_value * second_value
    return result


def _powers(
    series: tuple[Fraction, ...], maximum: int, order: int
) -> list[list[Fraction]]:
    one = [Fraction(1)] + [Fraction(0)] * order
    values = [one]
    for _ in range(maximum):
        values.append(_convolve(values[-1], series, order))
    return values


def _formal_exp(series: tuple[Fraction, ...], order: int) -> list[Fraction]:
    """exp(series) by summing f^k/k!, not by the producer's ODE recurrence."""
    if series[0] != 0:
        raise AssertionError("exponential needs a series without constant term")
    truncated = list(series[: order + 1])
    truncated.extend([Fraction(0)] * (order + 1 - len(truncated)))
    result = [Fraction(0)] * (order + 1)
    result[0] = Fraction(1)
    term = [Fraction(1)] + [Fraction(0)] * order
    for power in range(1, order + 1):
        term = _convolve(term, truncated, order)
        factorial = Fraction(math.factorial(power))
        for degree, value in enumerate(term):
            if value:
                result[degree] += value / factorial
    return result


def _valuation(series: tuple[Fraction, ...]) -> int:
    for degree, value in enumerate(series):
        if value:
            return degree
    raise AssertionError("series unexpectedly vanishes")


def _algebraic_matrix(
    series: tuple[Fraction, ...], degree_x: int, degree_f: int
) -> tuple[list[list[Fraction]], list[tuple[int, int]], list[int]]:
    equation_count = len(series)
    powers = _powers(series, degree_f, equation_count - 1)
    monomials = [
        (a, b)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    rows = []
    for order in range(equation_count):
        rows.append(
            [
                powers[b][order - a] if order >= a else Fraction(0)
                for a, b in monomials
            ]
        )
    valuation = _valuation(series)
    return rows, monomials, [a + b * valuation for a, b in monomials]


def _differential_matrix(
    series: tuple[Fraction, ...], degree_x: int, degree_f: int, degree_f_prime: int
) -> tuple[list[list[Fraction]], list[tuple[int, int, int]], list[int]]:
    # Exactly e43's derivative-safety convention: N known f coefficients give
    # orders 0,...,N-2 for f'.
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
            combined = _convolve(f_powers[b], derivative_powers[c], equation_count - 1)
            row.append(combined[order - a] if order >= a else Fraction(0))
        rows.append(row)
    valuation = _valuation(series)
    return rows, monomials, [
        a + b * valuation + c * (valuation - 1) for a, b, c in monomials
    ]


def _dot(row: list[Fraction], vector: list[int]) -> Fraction:
    return sum(
        (entry * coefficient for entry, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _replay_survivor(
    relation_class: str,
    degrees: dict,
    series_u: tuple[Fraction, ...],
) -> dict:
    if relation_class == "algebraic":
        rows, monomials, formal_valuations = _algebraic_matrix(
            series_u, degrees["degree_x"], degrees["degree_f"]
        )
        derivative_safe = False
    elif relation_class == "differential_algebraic_order_1":
        rows, monomials, formal_valuations = _differential_matrix(
            series_u,
            degrees["degree_x"],
            degrees["degree_f"],
            degrees["degree_f_prime"],
        )
        derivative_safe = True
    else:
        raise AssertionError(f"unexpected survivor family {relation_class}")

    unknown_count = len(monomials)
    training_count = unknown_count + 2
    training_matrix = _to_sympy_rows(rows[:training_count])
    training_basis = [_primitive_vector(vector) for vector in training_matrix.nullspace()]
    residual_constraints: list[list[Fraction]] = []
    first_refuting_order = None
    residual_at_refutation = None
    for order in range(training_count, len(rows)):
        residual = [_dot(rows[order], vector) for vector in training_basis]
        residual_constraints.append(residual)
        rank = _to_sympy_rows(residual_constraints).rank()
        if len(training_basis) - rank == 0 and first_refuting_order is None:
            first_refuting_order = order
            residual_at_refutation = residual

    full_matrix = _to_sympy_rows(rows)
    full_nullity = len(full_matrix.nullspace())
    zero_columns = [
        column
        for column in range(unknown_count)
        if all(row[column] == 0 for row in rows)
    ]
    result = {
        "training_kernel_basis": training_basis,
        "training_equation_count": training_count,
        "available_equation_count": len(rows),
        "full_nullity": full_nullity,
        "zero_column_indices": zero_columns,
        "formal_monomial_valuations": formal_valuations,
    }
    if full_nullity == 0:
        if first_refuting_order is None or residual_at_refutation is None:
            raise AssertionError("full-rank extension has no recorded first refuting order")
        result.update(
            {
                "verdict": "CANDIDATE_REFUTED_BY_HOLDOUT",
                "first_refuting_holdout_order_u": first_refuting_order,
                "residual_on_training_kernel_basis": [
                    str(value) for value in residual_at_refutation
                ],
            }
        )
    else:
        if full_nullity != len(zero_columns):
            raise AssertionError("survivor is not a literal-zero-column truncation kernel")
        next_equation_order = min(formal_valuations[column] for column in zero_columns)
        result.update(
            {
                "verdict": "STILL_TRUNCATION_UNOBSERVABLE",
                "next_testable_equation_order_u": next_equation_order,
                "required_coefficient_order_x": 2
                * (next_equation_order + (1 if derivative_safe else 0)),
            }
        )
    return result


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
    raise AssertionError(f"frozen survivor cell missing: {relation_class} {degrees}")


def main() -> None:
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    predecessor = json.loads(PREDECESSOR_PATH.read_text(encoding="utf-8"))
    structure = json.loads(STRUCTURE_PATH.read_text(encoding="utf-8"))
    data = payload["data"]

    assert data["achieved_order"] == TARGET_ORDER
    assert data["variable"] == "x"
    assert data["order_bound"]["budget"] == TARGET_BUDGET
    assert data["ordered_box_count"] == math.comb(TARGET_BUDGET, 3)
    assert all(check["passed"] for check in payload["checks"])

    canonical = _canonical_boxes(TARGET_BUDGET)
    ordered = _ordered_boxes(TARGET_BUDGET)
    assert len(ordered) == math.comb(TARGET_BUDGET, 3)
    assert data["canonical_box_list"] == [list(shape) for shape in canonical]
    assert data["finite_lattice_inventory"]["canonical_box_count"] == len(canonical)
    assert data["finite_lattice_inventory"]["ordered_box_count"] == len(ordered)
    assert (5, 5, 5) in canonical  # post-lift: the former wall box is now required
    wall = data["next_order_preflight"]
    assert wall["next_order"] == 56
    assert wall["next_order_inventory_unchanged"] is True
    assert wall["next_wall_order"] == 62
    assert wall["next_wall_canonical_shape"] == [5, 6, 6]
    assert wall["next_wall_cross_section_sites"] == 30
    assert wall["transfer_cross_section_guard"] == 25
    assert wall["launched"] is False

    producer_prefixes = {
        key: [int(value) for value in values]
        for key, values in data["box_polynomial_prefixes"].items()
    }
    independent_prefixes: dict[tuple[int, int, int], list[int]] = {}
    fresh_crt_records: dict[str, dict] = {}

    for shape in canonical:
        key = _shape_key(shape)
        cross_section = shape[0] * shape[1]
        if cross_section <= 12:
            independent = _independent_plus_polynomial(shape, TARGET_ORDER)
            if math.prod(shape) <= 16:
                assert independent == _brute_plus_polynomial(shape, TARGET_ORDER), shape
            if math.prod(shape) <= 62:
                engine = list(box_broken_bond_poly(shape, plus_boundary=True)[:PREFIX_LENGTH])
                engine.extend([0] * (PREFIX_LENGTH - len(engine)))
                assert independent == engine, shape
        else:
            assert math.prod(shape) > 62, shape
            independent, certificate = _independent_crt_prefix(shape, TARGET_ORDER)
            fresh_crt_records[key] = certificate
        assert independent == producer_prefixes[key], shape
        independent_prefixes[shape] = independent

    # The independently reconstructed high-cross-section boxes use fresh CRT
    # primes, while every lower-cross-section box was obtained exactly in
    # Python integers.  This is the raw-input audit before formal logarithms.
    expected_big = {
        _shape_key(shape) for shape in canonical if shape[0] * shape[1] > 12
    }
    assert set(fresh_crt_records) == expected_big
    for key, certificate in fresh_crt_records.items():
        assert set(certificate["primes"]).isdisjoint(DEFAULT_PRODUCER_PRIMES), key
        assert certificate["modulus_product"] > 2 * certificate["coefficient_bound"], key

    logs = {
        shape: _formal_log_via_inverse(prefix)
        for shape, prefix in independent_prefixes.items()
    }
    weights = {
        shape: _mixed_second_difference(shape, logs)
        for shape in ordered
    }
    bulk = [Fraction(0)] * PREFIX_LENGTH
    for shape, weight in weights.items():
        bound = 4 * sum(shape) - 6
        assert all(value == 0 for value in weight[:bound]), shape
        for degree, value in enumerate(weight):
            bulk[degree] += value

    artifact_coefficients = tuple(Fraction(value) for value in data["coefficients"])
    assert tuple(bulk) == artifact_coefficients
    predecessor_coefficients = tuple(
        Fraction(value) for value in predecessor["data"]["coefficients"]
    )
    assert tuple(bulk[: len(predecessor_coefficients)]) == predecessor_coefficients
    assert all(bulk[degree] == 0 for degree in range(1, PREFIX_LENGTH, 2))

    # Replay every legacy survivor in the frozen u=x^2 convention without
    # changing the U+2 training prefix.  The extended coefficients are only
    # holdouts.
    series_u = tuple(bulk[degree] for degree in range(0, PREFIX_LENGTH, 2))
    survivor_rows = structure["data"]["survivor_parent_audit"]["rows"]
    recorded_rows = data["survivor_reaudit"]["rows"]
    assert len(survivor_rows) == 18
    assert len(recorded_rows) == 18
    recorded_by_key = {
        (
            row["relation_class"],
            tuple(sorted(row["degrees"].items())),
        ): row
        for row in recorded_rows
    }
    replayed_rows: list[dict] = []
    for frozen in survivor_rows:
        relation_class = frozen["relation_class"]
        degrees = frozen["degrees"]
        frozen_cell = _frozen_survivor_cell(structure, relation_class, degrees)
        replayed = _replay_survivor(relation_class, degrees, series_u)
        assert replayed["training_kernel_basis"] == frozen_cell["training_kernel_basis"]
        assert replayed["training_equation_count"] == frozen_cell["training_equation_count"]
        key = (relation_class, tuple(sorted(degrees.items())))
        recorded = recorded_by_key[key]
        assert recorded["verdict"] == replayed["verdict"]
        assert recorded["full_nullity"] == replayed["full_nullity"]
        assert recorded["zero_column_indices"] == replayed["zero_column_indices"]
        if replayed["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT":
            assert (
                recorded["first_refuting_holdout_order_u"]
                == replayed["first_refuting_holdout_order_u"]
            )
            assert (
                recorded["residual_on_training_kernel_basis"]
                == replayed["residual_on_training_kernel_basis"]
            )
        else:
            assert (
                recorded["next_testable_equation_order_u"]
                == replayed["next_testable_equation_order_u"]
            )
            assert (
                recorded["required_coefficient_order_x"]
                == replayed["required_coefficient_order_x"]
            )
        replayed_rows.append(replayed)

    # Plausible scope/off-by-one defects must fail these gates.  With the
    # predecessor's 17 u-coefficients, the derivative-safe order 16 is not
    # available, so (f')^8 remains a literal zero column.  It becomes
    # refutable only after the new x^34 coefficient makes that equation safe.
    old_series_u = tuple(
        predecessor_coefficients[degree]
        for degree in range(0, len(predecessor_coefficients), 2)
    )
    old_da8 = _replay_survivor(
        "differential_algebraic_order_1",
        {"degree_x": 0, "degree_f": 0, "degree_f_prime": 8},
        old_series_u,
    )
    new_da8 = _replay_survivor(
        "differential_algebraic_order_1",
        {"degree_x": 0, "degree_f": 0, "degree_f_prime": 8},
        series_u,
    )
    assert old_da8["verdict"] == "STILL_TRUNCATION_UNOBSERVABLE"
    assert old_da8["zero_column_indices"] == [8]
    assert new_da8["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    assert new_da8["first_refuting_holdout_order_u"] == 16
    assert new_da8["residual_on_training_kernel_basis"][-1] == str(3**8)

    # Re-derive every CRT certificate from its own raw fields instead of
    # trusting the stored scalar inequality or Boolean: the covered shape set,
    # the site count, the coefficient bound, the modulus product, its bit
    # length, distinctness and primality of the moduli, and the strict
    # uniqueness inequality itself.  A swapped prime list, an inconsistent
    # product, or a dropped certificate row must fail here.
    certificates = data["crt_certificates"]
    expected_crt_shapes = {shape for shape in canonical if math.prod(shape) > 62}
    assert {tuple(item["shape"]) for item in certificates} == expected_crt_shapes
    assert len(certificates) == len(expected_crt_shapes)
    assert data["finite_lattice_inventory"]["crt_box_count"] == len(expected_crt_shapes)
    assert data["finite_lattice_inventory"]["int64_box_count"] == len(canonical) - len(
        expected_crt_shapes
    )
    for certificate in certificates:
        shape = tuple(certificate["shape"])
        site_count = math.prod(shape)
        bound = 1 << site_count
        moduli = certificate["prime_moduli"]
        product = math.prod(moduli)
        assert certificate["claim_tag"] == "[COMPUTATION]"
        assert certificate["sites"] == site_count, shape
        assert certificate["cross_section_sites"] == shape[0] * shape[1], shape
        assert certificate["coefficient_bound"] == bound, shape
        assert certificate["twice_coefficient_bound"] == 2 * bound, shape
        assert len(set(moduli)) == len(moduli), shape
        for modulus in moduli:
            assert 1 << 30 < modulus < 1 << 31, (shape, modulus)
            assert _is_prime_32(modulus), (shape, modulus)
        assert certificate["modulus_product"] == product, shape
        assert certificate["modulus_product_bits"] == product.bit_length(), shape
        assert product > 2 * bound, shape
        assert certificate["uniqueness_certificate_passed"] == (product > 2 * bound), shape
        assert certificate["full_sum_check_status"] == "not needed for truncated prefix"
        # This box was independently recomputed either by the exact pure-Python
        # DP (cross-section <= 12) or by a disjoint fresh prime set, so the
        # certificate is never checked against itself.
        if shape[0] * shape[1] > 12:
            fresh = fresh_crt_records[_shape_key(shape)]
            assert set(fresh["primes"]).isdisjoint(moduli), shape
            assert fresh["coefficient_bound"] == bound, shape
        else:
            assert _shape_key(shape) not in fresh_crt_records, shape
            assert independent_prefixes[shape] == producer_prefixes[_shape_key(shape)]

    # Independently replay the post-computation external witness: exponentiate
    # the series this test derived itself and compare against every published
    # value the artifact records, plus the recorded ordering of the two cell
    # kinds.  A dropped or mislabelled row must fail here.
    witness = data["external_witness"]
    assert witness["available"] is True
    assert witness["accessed_after_local_derivation"] is True
    lambda_local = _formal_exp(series_u, max(row["lambda_index_u"] for row in witness["rows"]))
    assert len(witness["rows"]) == 25
    for row in witness["rows"]:
        index = row["lambda_index_u"]
        assert row["x_order"] == 2 * index
        assert lambda_local[index].denominator == 1
        assert str(lambda_local[index]) == row["local_exponentiated_value"], index
        assert lambda_local[index] == Fraction(row["published_value"]), index
        assert row["match"] is True
        assert row["table_cell_kind"] in {"row_aligned", "merged_column_block_token"}
    # Table 2 prints Lambda_0 through u^26 = x^52 only; for the x^54 run the
    # witnessed new orders end at 52 and x^54 is explicitly beyond the printed
    # literature.  Assert exactly that, not full coverage.
    assert witness["witnessed_new_x_orders"] == list(range(34, 52 + 1, 2))
    assert 54 not in witness["witnessed_new_x_orders"]

    # The witness comparison is a producer gate, so the metadata must not deny
    # a validation role while still denying input/fit/selection use.  The
    # K_c benchmark claim is kept separate from the external-table claim.
    gate = next(
        check
        for check in payload["checks"]
        if check["name"] == "post_computation_external_table_witness"
    )
    assert gate["passed"] is True
    role = data["external_value_role"]
    assert data["kc_benchmark_used_for_fit_selection_or_validation"] is False
    assert "benchmark_used_for_fit_selection_or_validation" not in data
    assert role["used_as_computational_input"] is False
    assert role["used_for_fitting_selection_or_tuning"] is False
    assert role["used_as_post_computation_validation_gate"] is True
    assert "post_computation_external_table_witness" in role["detail"]
    assert "validation witness" in witness["coefficient_source_role"]
    assert "validation gate" in witness["coefficient_source_role"]

    # The proof note must carry the independently recomputed decisive numbers,
    # so narrative and computation cannot drift apart silently.
    proof_text = PROOF_PATH.read_text(encoding="utf-8")
    for tag in ("[THEOREM]", "[COMPUTATION]", "[UNRESOLVED]"):
        assert tag in proof_text, tag
    assert "[EXTERNAL]" in proof_text
    assert "post-computation validation witness" in proof_text
    assert "is a producer gate" in proof_text
    assert "No benchmark value was used to fit, select, tune or validate" not in proof_text
    assert "4*(a+b+c)-6" in proof_text
    assert f"a+b+c <= {TARGET_BUDGET}" in proof_text
    assert "5x5x5" in proof_text
    for degree in range(34, PREFIX_LENGTH, 2):
        assert f"x^{degree} = {bulk[degree]}" in proof_text, degree
    refuted = [
        row
        for row in recorded_rows
        if row["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ]
    assert f"{len(refuted)} of 18" in proof_text

    # Parse the note's re-audit table row by row and align it with the
    # independently replayed verdicts: a misaligned, dropped or edited row
    # must fail, which a bare substring search would not catch.
    verdict_labels = {"CANDIDATE_REFUTED_BY_HOLDOUT", "STILL_TRUNCATION_UNOBSERVABLE"}
    all_note_rows = [
        [cell.strip() for cell in line.strip().strip("|").split("|")]
        for line in proof_text.splitlines()
        if line.startswith("| ")
    ]
    table_rows = [
        cells
        for cells in all_note_rows
        if len(cells) == 7 and cells[0].isdigit() and cells[4] in verdict_labels
    ]
    assert len(table_rows) == 18, len(table_rows)
    for position, (cells, replayed) in enumerate(
        zip(table_rows, replayed_rows, strict=True), start=1
    ):
        index, relation_class, _degrees, zero_columns, verdict, marker, residual = cells
        assert int(index) == position
        assert relation_class == recorded_rows[position - 1]["relation_class"]
        assert zero_columns == "[" + ", ".join(
            str(value)
            for value in recorded_rows[position - 1]["frozen_zero_column_indices"]
        ) + "]"
        assert verdict == replayed["verdict"]
        if replayed["verdict"] == "CANDIDATE_REFUTED_BY_HOLDOUT":
            assert marker == f"refuted at u^{replayed['first_refuting_holdout_order_u']}"
            assert residual == "`" + ", ".join(
                replayed["residual_on_training_kernel_basis"]
            ) + "`"
        else:
            assert marker == (
                f"next testable u^{replayed['next_testable_equation_order_u']}"
                f" (requires x^{replayed['required_coefficient_order_x']})"
            )
            assert residual == "`0 at every available holdout`"
    assert sum(
        1 for cells in table_rows if cells[4] == "CANDIDATE_REFUTED_BY_HOLDOUT"
    ) == len(refuted)

    digest = hashlib.sha256(RESULT_PATH.read_bytes()).hexdigest()
    print(
        "OK test_lt_x54: independently rebuilt x^0..x^54; "
        f"replayed 18 survivors; artifact_sha256={digest}"
    )


if __name__ == "__main__":
    main()
