#!/usr/bin/env python3
"""Exact local Callen expansion and finite-graph identity systems.

The producer uses the high-temperature variable ``v=tanh(K)`` throughout.
It proves the permutation-symmetric six-neighbour identity from the three
positive local-field orbits 2, 4, and 6, then builds every Callen equation on
three small open graphs in the complete monomial-correlator basis.  Exact
spin enumeration checks those equations as identities in Q(v).  Modular ranks
at several rational specializations are lower bounds only; an explicit
symbolic Gibbs null vector supplies the matching rational ceiling.
"""

from __future__ import annotations

import itertools
import json
import math
import platform
import resource
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e240_callen_identity_system.py"
OUTPUT = ROOT / "results" / "correlations" / "callen_identity_system.json"
RSS_CAP_BYTES = 2 * 1024**3
CPU_CAP_SECONDS = 120.0
RANK_SAMPLES = (
    (Fraction(1, 2), 1_000_003),
    (Fraction(1, 3), 1_000_033),
    (Fraction(2, 5), 1_000_037),
)

Poly = tuple[int, ...]
QPoly = list[Fraction]


def trim(values: Iterable[int]) -> Poly:
    out = [int(value) for value in values]
    if not out:
        return (0,)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_add(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def poly_sub(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return trim(
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def poly_scale(poly: Poly, scalar: int) -> Poly:
    return trim(scalar * coefficient for coefficient in poly)


def poly_multiply(left: Poly, right: Poly) -> Poly:
    if left == (0,) or right == (0,):
        return (0,)
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def poly_eval_mod(poly: Poly, value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(poly):
        result = (result * value + coefficient) % prime
    return result


def qtrim(values: Iterable[Fraction]) -> QPoly:
    out = list(values)
    if not out:
        return [Fraction(0)]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def q_divmod(dividend: QPoly, divisor: QPoly) -> tuple[QPoly, QPoly]:
    dividend = qtrim(dividend)
    divisor = qtrim(divisor)
    if divisor == [0]:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(dividend) < len(divisor):
        return [Fraction(0)], dividend
    quotient = [Fraction(0)] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while remainder != [0] and len(remainder) >= len(divisor):
        degree = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[degree] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + degree] -= coefficient * value
        remainder = qtrim(remainder)
    return qtrim(quotient), remainder


def q_gcd(left: QPoly, right: QPoly) -> QPoly:
    left = qtrim(left)
    right = qtrim(right)
    while right != [0]:
        _, remainder = q_divmod(left, right)
        left, right = right, remainder
    if left == [0]:
        return [Fraction(1)]
    leading = left[-1]
    return qtrim(value / leading for value in left)


def normalize_rational(numerator: Poly, denominator: Poly) -> tuple[Poly, Poly]:
    numerator = trim(numerator)
    denominator = trim(denominator)
    if denominator == (0,):
        raise ZeroDivisionError("zero rational-function denominator")
    if numerator == (0,):
        return (0,), (1,)

    qnumerator = [Fraction(value) for value in numerator]
    qdenominator = [Fraction(value) for value in denominator]
    common = q_gcd(qnumerator, qdenominator)
    qnumerator, numerator_remainder = q_divmod(qnumerator, common)
    qdenominator, denominator_remainder = q_divmod(qdenominator, common)
    if numerator_remainder != [0] or denominator_remainder != [0]:
        raise AssertionError("polynomial gcd did not divide exactly")

    multiplier = 1
    for value in qnumerator + qdenominator:
        multiplier = math.lcm(multiplier, value.denominator)
    integer_numerator = [int(value * multiplier) for value in qnumerator]
    integer_denominator = [int(value * multiplier) for value in qdenominator]
    content = 0
    for value in integer_numerator + integer_denominator:
        content = math.gcd(content, abs(value))
    integer_numerator = [value // content for value in integer_numerator]
    integer_denominator = [value // content for value in integer_denominator]
    if integer_denominator[-1] < 0:
        integer_numerator = [-value for value in integer_numerator]
        integer_denominator = [-value for value in integer_denominator]
    return trim(integer_numerator), trim(integer_denominator)


@dataclass(frozen=True)
class RationalFunction:
    numerator: Poly
    denominator: Poly = (1,)

    def __post_init__(self) -> None:
        numerator, denominator = normalize_rational(
            self.numerator, self.denominator
        )
        object.__setattr__(self, "numerator", numerator)
        object.__setattr__(self, "denominator", denominator)

    @staticmethod
    def zero() -> "RationalFunction":
        return RationalFunction((0,))

    @staticmethod
    def one() -> "RationalFunction":
        return RationalFunction((1,))

    def __add__(self, other: "RationalFunction") -> "RationalFunction":
        return RationalFunction(
            poly_add(
                poly_multiply(self.numerator, other.denominator),
                poly_multiply(other.numerator, self.denominator),
            ),
            poly_multiply(self.denominator, other.denominator),
        )

    def __neg__(self) -> "RationalFunction":
        return RationalFunction(poly_scale(self.numerator, -1), self.denominator)

    def __sub__(self, other: "RationalFunction") -> "RationalFunction":
        return self + (-other)

    def __mul__(self, other: "RationalFunction") -> "RationalFunction":
        return RationalFunction(
            poly_multiply(self.numerator, other.numerator),
            poly_multiply(self.denominator, other.denominator),
        )

    def scaled(self, scalar: Fraction | int) -> "RationalFunction":
        scalar = Fraction(scalar)
        return RationalFunction(
            poly_scale(self.numerator, scalar.numerator),
            poly_scale(self.denominator, scalar.denominator),
        )

    def eval_mod(self, value: Fraction, prime: int) -> int:
        denominator_mod = value.denominator % prime
        if denominator_mod == 0:
            raise ZeroDivisionError("specialization denominator vanished modulo prime")
        v_mod = value.numerator % prime * pow(denominator_mod, -1, prime) % prime
        function_denominator = poly_eval_mod(self.denominator, v_mod, prime)
        if function_denominator == 0:
            raise ZeroDivisionError("local coefficient pole modulo prime")
        return (
            poly_eval_mod(self.numerator, v_mod, prime)
            * pow(function_denominator, -1, prime)
            % prime
        )

    def as_json(self) -> dict[str, list[int]]:
        return {
            "numerator_coefficients_ascending": list(self.numerator),
            "denominator_coefficients_ascending": list(self.denominator),
        }


def tanh_multiple(field: int) -> RationalFunction:
    """Return tanh(field*K) in Q(v), where v=tanh(K)."""

    if field == 0:
        return RationalFunction.zero()
    sign = 1 if field > 0 else -1
    degree = abs(field)
    numerator = [0] * (degree + 1)
    denominator = [0] * (degree + 1)
    for power in range(degree + 1):
        if power % 2:
            numerator[power] = sign * math.comb(degree, power)
        else:
            denominator[power] = math.comb(degree, power)
    return RationalFunction(trim(numerator), trim(denominator))


def spin(state: int, site: int) -> int:
    return 1 if (state >> site) & 1 else -1


def elementary_symmetric(spins: tuple[int, ...], degree: int) -> int:
    return sum(
        math.prod(spins[index] for index in subset)
        for subset in itertools.combinations(range(len(spins)), degree)
    )


def walsh_field_counts(degree: int, subset_size: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for state in range(1 << degree):
        spins = tuple(spin(state, site) for site in range(degree))
        character = math.prod(spins[:subset_size])
        field = sum(spins)
        counts[field] = counts.get(field, 0) + character
    return dict(sorted(counts.items()))


def local_coefficient(degree: int, subset_size: int) -> RationalFunction:
    """Walsh-interpolate one symmetric odd coefficient from raw spin states."""

    total = RationalFunction.zero()
    for field, multiplicity in walsh_field_counts(degree, subset_size).items():
        total = total + tanh_multiple(field).scaled(multiplicity)
    return total.scaled(Fraction(1, 1 << degree))


def local_coefficients(degree: int) -> dict[int, RationalFunction]:
    return {
        subset_size: local_coefficient(degree, subset_size)
        for subset_size in range(1, degree + 1, 2)
    }


def determinant_three(matrix: list[list[int]]) -> int:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def inverse_fraction(matrix: list[list[int]]) -> list[list[Fraction]]:
    dimension = len(matrix)
    augmented = [
        [Fraction(value) for value in row]
        + [Fraction(int(row_index == column)) for column in range(dimension)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(dimension):
        pivot = next(
            row for row in range(column, dimension) if augmented[row][column]
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [value / pivot_value for value in augmented[column]]
        for row in range(dimension):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    left - factor * right
                    for left, right in zip(augmented[row], augmented[column])
                ]
    return [row[dimension:] for row in augmented]


def combine_rational_functions(
    coefficients: list[Fraction], functions: list[RationalFunction]
) -> RationalFunction:
    total = RationalFunction.zero()
    for coefficient, function in zip(coefficients, functions):
        total = total + function.scaled(coefficient)
    return total


@dataclass(frozen=True)
class Graph:
    name: str
    coordinates: tuple[tuple[int, ...], ...]
    edges: tuple[tuple[int, int], ...]

    @property
    def order(self) -> int:
        return len(self.coordinates)

    @property
    def neighbours(self) -> tuple[tuple[int, ...], ...]:
        adjacency = [[] for _ in range(self.order)]
        for left, right in self.edges:
            adjacency[left].append(right)
            adjacency[right].append(left)
        return tuple(tuple(sorted(row)) for row in adjacency)


def grid_graph(name: str, shape: tuple[int, ...]) -> Graph:
    coordinates = tuple(itertools.product(*(range(length) for length in shape)))
    index = {coordinate: site for site, coordinate in enumerate(coordinates)}
    edges: list[tuple[int, int]] = []
    for coordinate in coordinates:
        for axis in range(len(shape)):
            neighbour = list(coordinate)
            neighbour[axis] += 1
            candidate = tuple(neighbour)
            if candidate in index:
                edges.append((index[coordinate], index[candidate]))
    return Graph(name, coordinates, tuple(edges))


def state_weight_polynomial(state: int, graph: Graph) -> Poly:
    weight = (1,)
    for left, right in graph.edges:
        edge_spin = spin(state, left) * spin(state, right)
        weight = poly_multiply(weight, (1, edge_spin))
    return weight


def correlation_numerators(graph: Graph) -> list[Poly]:
    """All unnormalised spin moments via an exact Walsh-Hadamard transform."""

    values = [state_weight_polynomial(state, graph) for state in range(1 << graph.order)]
    stride = 1
    while stride < len(values):
        for start in range(0, len(values), 2 * stride):
            for offset in range(stride):
                left = start + offset
                right = left + stride
                a, b = values[left], values[right]
                values[left] = poly_add(a, b)
                values[right] = poly_sub(a, b)
        stride *= 2
    for mask in range(len(values)):
        if mask.bit_count() % 2:
            values[mask] = poly_scale(values[mask], -1)
    return values


def neighbour_subset_masks(
    neighbours: tuple[int, ...], subset_size: int
) -> tuple[int, ...]:
    return tuple(
        sum(1 << site for site in subset)
        for subset in itertools.combinations(neighbours, subset_size)
    )


def exact_identity_audit(
    graph: Graph, numerators: list[Poly], coefficients: dict[int, dict[int, RationalFunction]]
) -> dict[str, object]:
    zero_rows = 0
    nonzero_rows: list[tuple[int, int, RationalFunction]] = []
    negative_nonzero = 0
    first_negative: dict[str, object] | None = None
    for vertex, neighbours in enumerate(graph.neighbours):
        local = coefficients[len(neighbours)]
        masks_by_size = {
            size: neighbour_subset_masks(neighbours, size) for size in local
        }
        for observable in range(1 << graph.order):
            if (observable >> vertex) & 1:
                continue
            residual = RationalFunction(numerators[observable ^ (1 << vertex)])
            for size, coefficient in local.items():
                orbit_sum = (0,)
                for subset_mask in masks_by_size[size]:
                    orbit_sum = poly_add(
                        orbit_sum, numerators[observable ^ subset_mask]
                    )
                residual = residual - coefficient * RationalFunction(orbit_sum)
            if residual.numerator == (0,):
                zero_rows += 1
            else:
                nonzero_rows.append((vertex, observable, residual))

            # Negative control: perturb only the one-neighbour coefficient by +1.
            singleton_sum = (0,)
            for subset_mask in masks_by_size[1]:
                singleton_sum = poly_add(
                    singleton_sum, numerators[observable ^ subset_mask]
                )
            if singleton_sum != (0,):
                negative_nonzero += 1
                if first_negative is None:
                    first_negative = {
                        "vertex": vertex,
                        "observable_mask": observable,
                        "residual_numerator_coefficients_ascending": list(
                            poly_scale(singleton_sum, -1)
                        ),
                    }

    row_count = graph.order * (1 << (graph.order - 1))
    if nonzero_rows:
        vertex, observable, residual = nonzero_rows[0]
        first_failure = {
            "vertex": vertex,
            "observable_mask": observable,
            **residual.as_json(),
        }
    else:
        first_failure = None
    return {
        "identity_rows_checked": row_count,
        "symbolic_zero_rows": zero_rows,
        "symbolic_nonzero_rows": len(nonzero_rows),
        "first_failure": first_failure,
        "negative_control": {
            "mutation": "a_(degree,1)(v) -> a_(degree,1)(v) + 1 at every vertex",
            "nonzero_rows": negative_nonzero,
            "rejected": negative_nonzero > 0,
            "first_nonzero_residual": first_negative,
        },
    }


def add_mod_entry(row: dict[int, int], column: int, value: int, prime: int) -> None:
    updated = (row.get(column, 0) + value) % prime
    if updated:
        row[column] = updated
    else:
        row.pop(column, None)


def modular_rows(
    graph: Graph,
    coefficients: dict[int, dict[int, RationalFunction]],
    value: Fraction,
    prime: int,
    correlator_parity: int,
) -> list[dict[int, int]]:
    rows: list[dict[int, int]] = []
    for vertex, neighbours in enumerate(graph.neighbours):
        local = coefficients[len(neighbours)]
        evaluated = {
            size: coefficient.eval_mod(value, prime)
            for size, coefficient in local.items()
        }
        masks_by_size = {
            size: neighbour_subset_masks(neighbours, size) for size in local
        }
        for observable in range(1 << graph.order):
            if (observable >> vertex) & 1:
                continue
            if (observable.bit_count() + 1) % 2 != correlator_parity:
                continue
            row: dict[int, int] = {observable ^ (1 << vertex): 1}
            for size, coefficient in evaluated.items():
                for subset_mask in masks_by_size[size]:
                    add_mod_entry(
                        row, observable ^ subset_mask, -coefficient, prime
                    )
            rows.append(row)
    return rows


def sparse_modular_rank(
    rows: list[dict[int, int]], prime: int
) -> tuple[int, list[int], list[int]]:
    pivots: dict[int, dict[int, int]] = {}
    pivot_sources: dict[int, int] = {}
    for source_index, source in enumerate(rows):
        row = {column: value % prime for column, value in source.items() if value % prime}
        while row:
            pivot = min(row)
            basis = pivots.get(pivot)
            if basis is None:
                inverse = pow(row[pivot], -1, prime)
                row = {
                    column: value * inverse % prime
                    for column, value in row.items()
                    if value * inverse % prime
                }
                pivots[pivot] = row
                pivot_sources[pivot] = source_index
                break
            factor = row[pivot]
            for column, value in basis.items():
                add_mod_entry(row, column, -factor * value, prime)
    ordered = sorted(pivots)
    return (
        len(ordered),
        ordered,
        [pivot_sources[pivot] for pivot in ordered],
    )


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def rank_audit(
    graph: Graph, coefficients: dict[int, dict[int, RationalFunction]]
) -> dict[str, object]:
    columns = 1 << graph.order
    parity_columns = columns // 2
    specialization_rows: list[dict[str, object]] = []
    for value, prime in RANK_SAMPLES:
        if not is_prime(prime):
            raise AssertionError(f"rank modulus {prime} is not prime")
        even_rows = modular_rows(graph, coefficients, value, prime, 0)
        odd_rows = modular_rows(graph, coefficients, value, prime, 1)
        even_rank, even_pivots, even_sources = sparse_modular_rank(even_rows, prime)
        odd_rank, odd_pivots, odd_sources = sparse_modular_rank(odd_rows, prime)
        full_rank = even_rank + odd_rank
        specialization_rows.append(
            {
                "v": {
                    "numerator": value.numerator,
                    "denominator": value.denominator,
                },
                "prime": prime,
                "good_reduction": True,
                "modular_lower_bounds": {
                    "even_correlator_block": even_rank,
                    "odd_correlator_block": odd_rank,
                    "full_system": full_rank,
                },
                "matching_rational_ceilings": {
                    "even_correlator_block": parity_columns - 1,
                    "odd_correlator_block": parity_columns,
                    "full_system": columns - 1,
                },
                "exact_rank_at_this_rational_v": full_rank
                if full_rank == columns - 1
                else None,
                "pivot_certificate": {
                    "even_columns": even_pivots,
                    "even_source_row_indices": even_sources,
                    "odd_columns": odd_pivots,
                    "odd_source_row_indices": odd_sources,
                },
            }
        )
    lower_bound = max(
        int(row["modular_lower_bounds"]["full_system"])
        for row in specialization_rows
    )
    exact = lower_bound == columns - 1
    return {
        "claim_tag": "[COMPUTATION]",
        "field": "Q(v)",
        "basis": "all 2^|V| monomial correlators C_A=<prod_(i in A) sigma_i>, including C_empty",
        "equations": "one row for every x and every monomial O not containing x",
        "matrix_shape": [graph.order * (1 << (graph.order - 1)), columns],
        "parity_block_columns": {
            "even": parity_columns,
            "odd": parity_columns,
        },
        "rational_ceiling_certificate": {
            "full_rank_at_most": columns - 1,
            "even_block_rank_at_most": parity_columns - 1,
            "odd_block_rank_at_most": parity_columns,
            "reason": "exact symbolic spin enumeration gives a nonzero Q(v) Gibbs moment vector in the kernel; its empty-set component is the nonzero partition polynomial",
        },
        "specialization_lower_bounds": specialization_rows,
        "generic_rank": columns - 1 if exact else None,
        "generic_nullity": 1 if exact else None,
        "even_block_generic_rank": parity_columns - 1 if exact else None,
        "even_block_generic_nullity": 1 if exact else None,
        "odd_block_generic_rank": parity_columns if exact else None,
        "odd_block_generic_nullity": 0 if exact else None,
        "rank_closed_exactly": exact,
        "classification": "generically_saturated_up_to_the_single_normalization_scale"
        if exact
        else "unresolved",
    }


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    rss = peak_rss_bytes()
    if elapsed > CPU_CAP_SECONDS:
        raise RuntimeError(f"CPU guard exceeded at {stage}: {elapsed:.3f}s")
    if rss > RSS_CAP_BYTES:
        raise MemoryError(f"RSS guard exceeded at {stage}: {rss} bytes")


def record_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    six_counts = {
        size: walsh_field_counts(6, size) for size in (1, 3, 5)
    }
    six_coefficients = local_coefficients(6)
    expected_coefficients = {
        1: RationalFunction(
            (0, 1, 0, 16, 0, 46, 0, 16, 0, 1),
            (1, 0, 21, 0, 106, 0, 106, 0, 21, 0, 1),
        ),
        3: RationalFunction(
            (0, 0, 0, -2), (1, 0, 15, 0, 15, 0, 1)
        ),
        5: RationalFunction(
            (0, 0, 0, 0, 0, 16),
            (1, 0, 21, 0, 106, 0, 106, 0, 21, 0, 1),
        ),
    }
    record_check(
        checks,
        "six_neighbour_walsh_counts_use_only_requested_orbits",
        set().union(*(set(row) for row in six_counts.values()))
        == {-6, -4, -2, 0, 2, 4, 6}
        and all(row[0] == 0 for row in six_counts.values()),
        str(six_counts),
    )
    record_check(
        checks,
        "six_neighbour_coefficients_reduced_exactly",
        six_coefficients == expected_coefficients,
        "c1,c3,c5 match the reduced Q(v) rational functions",
    )

    orbit_fields = (2, 4, 6)
    orbit_spins = {
        field: tuple([1] * ((6 + field) // 2) + [-1] * ((6 - field) // 2))
        for field in orbit_fields
    }
    orbit_matrix = [
        [elementary_symmetric(orbit_spins[field], size) for size in (1, 3, 5)]
        for field in orbit_fields
    ]
    orbit_determinant = determinant_three(orbit_matrix)
    orbit_inverse = inverse_fraction(orbit_matrix)
    interpolated_from_orbits = [
        combine_rational_functions(
            row, [tanh_multiple(field) for field in orbit_fields]
        )
        for row in orbit_inverse
    ]
    orbit_residuals = {}
    for field, row in zip(orbit_fields, orbit_matrix):
        value = combine_rational_functions(
            [Fraction(entry) for entry in row],
            [six_coefficients[size] for size in (1, 3, 5)],
        )
        orbit_residuals[str(field)] = value - tanh_multiple(field)
    record_check(
        checks,
        "three_positive_orbits_uniquely_interpolate_master_identity",
        orbit_determinant == 512
        and interpolated_from_orbits
        == [six_coefficients[size] for size in (1, 3, 5)]
        and all(value.numerator == (0,) for value in orbit_residuals.values()),
        f"orbit determinant={orbit_determinant}",
    )

    truth_residuals: list[RationalFunction] = []
    for state in range(64):
        spins = tuple(spin(state, site) for site in range(6))
        value = RationalFunction.zero()
        for size in (1, 3, 5):
            value = value + six_coefficients[size].scaled(
                elementary_symmetric(spins, size)
            )
        truth_residuals.append(value - tanh_multiple(sum(spins)))
    record_check(
        checks,
        "master_identity_all_64_spin_assignments",
        all(value.numerator == (0,) for value in truth_residuals),
        "64/64 exact Q(v) residuals vanish",
    )

    perturbation = RationalFunction(
        (0, 0, 0, 0, 0, 1), six_coefficients[5].denominator
    )
    negative_orbit_residuals: dict[str, RationalFunction] = {}
    for field, row in zip(orbit_fields, orbit_matrix):
        mutated = dict(six_coefficients)
        mutated[5] = mutated[5] + perturbation
        value = combine_rational_functions(
            [Fraction(entry) for entry in row],
            [mutated[size] for size in (1, 3, 5)],
        )
        negative_orbit_residuals[str(field)] = value - tanh_multiple(field)
    record_check(
        checks,
        "master_identity_negative_control_rejected",
        all(value.numerator != (0,) for value in negative_orbit_residuals.values()),
        "adding v^5/D to c5 leaves nonzero residuals on fields 2,4,6",
    )
    guard(started, "local master identity")

    graphs = (
        grid_graph("chain_3", (3,)),
        grid_graph("square_open_2x2", (2, 2)),
        grid_graph("cube_open_2x2x2", (2, 2, 2)),
    )
    required_degrees = sorted(
        {len(neighbours) for graph in graphs for neighbours in graph.neighbours}
    )
    coefficients = {degree: local_coefficients(degree) for degree in required_degrees}
    graph_rows: list[dict[str, object]] = []
    for graph in graphs:
        numerators = correlation_numerators(graph)
        partition = numerators[0]
        identity = exact_identity_audit(graph, numerators, coefficients)
        odd_numerators_zero = all(
            numerator == (0,)
            for mask, numerator in enumerate(numerators)
            if mask.bit_count() % 2
        )
        rank = rank_audit(graph, coefficients)
        row = {
            "name": graph.name,
            "claim_tag": "[COMPUTATION]",
            "vertex_count": graph.order,
            "edge_count": len(graph.edges),
            "coordinates": [list(coordinate) for coordinate in graph.coordinates],
            "edges": [list(edge) for edge in graph.edges],
            "degree_sequence": sorted(len(row) for row in graph.neighbours),
            "complete_correlator_basis": {
                "column_count": 1 << graph.order,
                "column_encoding": "integer subset mask 0..2^|V|-1",
                "includes_empty_set_normalization": True,
            },
            "exact_enumeration": {
                "arithmetic": "integer polynomials in v from raw products prod_edges(1+v*sigma_u*sigma_v)",
                "state_count": 1 << graph.order,
                "partition_polynomial_coefficients_ascending": list(partition),
                "partition_constant_term": partition[0],
                "all_odd_moment_numerators_zero": odd_numerators_zero,
                **identity,
            },
            "rank": rank,
        }
        graph_rows.append(row)
        row_count = graph.order * (1 << (graph.order - 1))
        record_check(
            checks,
            f"{graph.name}_every_callen_row_exact_enumeration",
            identity["symbolic_zero_rows"] == row_count
            and identity["symbolic_nonzero_rows"] == 0
            and odd_numerators_zero,
            f"{identity['symbolic_zero_rows']}/{row_count} Q(v) rows vanish",
        )
        record_check(
            checks,
            f"{graph.name}_negative_control_rejected",
            bool(identity["negative_control"]["rejected"]),
            f"nonzero mutated rows={identity['negative_control']['nonzero_rows']}",
        )
        record_check(
            checks,
            f"{graph.name}_generic_rank_closed",
            bool(rank["rank_closed_exactly"])
            and rank["generic_rank"] == (1 << graph.order) - 1
            and all(
                sample["exact_rank_at_this_rational_v"]
                == (1 << graph.order) - 1
                for sample in rank["specialization_lower_bounds"]
            ),
            f"rank/nullity={rank['generic_rank']}/{rank['generic_nullity']}",
        )
        guard(started, graph.name)

    if not all(row["passed"] for row in checks):
        failures = [row["name"] for row in checks if not row["passed"]]
        raise AssertionError(f"producer checks failed: {failures}")

    theorem_coefficients = {
        "c1": {
            **six_coefficients[1].as_json(),
            "factored": "v*(v^8+16*v^6+46*v^4+16*v^2+1)/((1+v^2)*(1+6*v^2+v^4)*(1+14*v^2+v^4))",
        },
        "c3": {
            **six_coefficients[3].as_json(),
            "factored": "-2*v^3/((1+v^2)*(1+14*v^2+v^4))",
        },
        "c5": {
            **six_coefficients[5].as_json(),
            "factored": "16*v^5/((1+v^2)*(1+6*v^2+v^4)*(1+14*v^2+v^4))",
        },
    }
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": SCRIPT,
            "interpreter": sys.executable,
            "arithmetic": "exact integer polynomials, rational functions over Q(v), and finite-field rank lower bounds closed by rational ceilings",
            "single_process": True,
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": peak_rss_bytes(),
            "rss_cap_bytes": RSS_CAP_BYTES,
            "cpu_cap_seconds": CPU_CAP_SECONDS,
            "benchmark_used": False,
        },
        "data": {
            "local_master_identity": {
                "claim_tag": "[THEOREM]",
                "variable": "v=tanh(K)",
                "statement": "For every six spins sigma_i in {+1,-1}, tanh(K*sum_i sigma_i)=c1(v)*e1+c3(v)*e3+c5(v)*e5, where e_j is the elementary symmetric spin polynomial of degree j.",
                "coefficients": theorem_coefficients,
                "interpolation_proof": {
                    "positive_local_fields_only": list(orbit_fields),
                    "negative_fields": "follow by oddness",
                    "zero_field": "both sides vanish by oddness",
                    "basis_order": ["e1", "e3", "e5"],
                    "orbit_matrix_rows_fields_2_4_6": orbit_matrix,
                    "orbit_matrix_determinant": orbit_determinant,
                    "orbit_matrix_inverse": [
                        [
                            {
                                "numerator": value.numerator,
                                "denominator": value.denominator,
                            }
                            for value in row
                        ]
                        for row in orbit_inverse
                    ],
                    "walsh_field_counts": {
                        f"e{size}": {str(field): value for field, value in counts.items()}
                        for size, counts in six_counts.items()
                    },
                    "orbit_residual_rational_functions": {
                        field: residual.as_json()
                        for field, residual in orbit_residuals.items()
                    },
                    "all_spin_assignments_checked": len(truth_residuals),
                    "all_spin_residuals_zero": all(
                        value.numerator == (0,) for value in truth_residuals
                    ),
                },
                "negative_control": {
                    "mutation": "c5 -> c5 + v^5/D, with D the stored c5 denominator",
                    "orbit_residuals": {
                        field: residual.as_json()
                        for field, residual in negative_orbit_residuals.items()
                    },
                    "rejected_on_every_positive_orbit": all(
                        value.numerator != (0,)
                        for value in negative_orbit_residuals.values()
                    ),
                },
                "scope": "local six-spin identity in Q(v), equivalently for real finite K after v=tanh(K)",
            },
            "callen_system_definition": {
                "claim_tag": "[EXACT IDENTITY]",
                "statement": "On a finite loopless undirected zero-field graph with uniform internal coupling K, every vertex spin dynamical, and free boundary, for each vertex x and each spin monomial O independent of sigma_x, <sigma_x O>=<tanh(K*sum_(y~x) sigma_y) O>.",
                "derivation": "sum sigma_x conditionally at fixed spins off x; the two weights are proportional to exp(+/- K*sum_(y~x) sigma_y)",
                "linear_row": "C_(O symmetric_difference {x}) - sum_(B subset N(x), |B| odd) a_(deg(x),|B|)(v) C_(O symmetric_difference B)=0",
                "basis_convention": "all subset monomials; products reduce by symmetric difference because sigma_i^2=1",
                "local_coefficients_used": {
                    str(degree): {
                        str(size): function.as_json()
                        for size, function in local.items()
                    }
                    for degree, local in coefficients.items()
                },
            },
            "finite_all_graph_reconstruction": {
                "claim_tag": "[THEOREM]",
                "statement": "Let G be any finite loopless undirected graph, possibly disconnected, with n dynamical spins, uniform coupling K on every internal edge, free boundary, and no external field. Over Q(v), v=tanh(K), the complete Callen system with all monomials O independent of x has rank 2^n-1 and its nullspace is the unnormalised Gibbs moment vector up to scale.",
                "proof": [
                    "Monomials O on the n-1 spins off x form a basis of every function of those spins, so the Callen rows are equivalent to the pointwise formal relation q(+1,eta)*(1-tanh(K*h_x(eta)))=q(-1,eta)*(1+tanh(K*h_x(eta)).",
                    "The nonzero formal ratio is q(+1,eta)/q(-1,eta)=exp(2K*h_x(eta)); no positivity or probability interpretation is needed for the Q(v) nullspace.",
                    "These rational ratios determine all 2^n configuration weights from one reference weight because the spin-configuration hypercube is connected, even if G is disconnected; consistency is witnessed by the unnormalised finite Gibbs weights.",
                    "The Walsh transform from configuration weights to all monomial correlators is invertible in characteristic zero, so the correlator-system nullity is one and the rank is 2^n-1.",
                ],
                "field_scope": "The rank statement is over Q(v). At a specialization v0, require all local coefficients regular and all flip factors 1+/-tanh(K*h) nonzero; in the uniform model this includes v0!=+/-1. Require Z_G(v0)!=0 before normalizing C_empty=1 or calling the vector a Gibbs probability. Every real finite K has |v0|<1 and satisfies all conditions.",
                "model_scope": "Fixed boundary spins, external/site-dependent fields, self-loops, or nonuniform couplings change h_x or its coefficient functions and are not included.",
                "complexity_scope": "the theorem is an exact finite-volume reconstruction but uses 2^n correlators and does not compress the thermodynamic problem",
            },
            "finite_graph_computations": graph_rows,
            "finite_rank_conclusion": {
                "claim_tag": "[COMPUTATION]",
                "graphs": [graph.name for graph in graphs],
                "result": "Each declared complete finite correlator system has generic nullity one, exactly the overall normalization scale; after C_empty=1 it is generically saturated rather than a smaller correlator compression.",
                "scope": "only the three explicitly enumerated finite open graphs",
            },
            "scope": {
                "proved": [
                    "the local six-neighbour master identity over Q(v)",
                    "the finite-graph Callen conditional-expectation identity",
                    "rank 2^n-1 and nullity one for the complete Callen correlator system on every finite zero-field graph",
                ],
                "computed_exactly": [
                    "all Callen rows and independent finite rank controls for the open three-site chain, open 2x2 square, and open 2x2x2 cube",
                ],
                "not_claimed": [
                    "a subexponential correlator compression",
                    "an infinite-volume closure or solution",
                    "any endpoint, transition point, or criticality statement",
                ],
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in payload["checks"]:
        status = "PASS" if row["passed"] else "FAIL"
        print(f"[{status}] {row['name']}: {row['detail']}")
    print(
        f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks; "
        f"cpu={payload['meta']['process_cpu_seconds']:.3f}s; "
        f"rss={payload['meta']['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
