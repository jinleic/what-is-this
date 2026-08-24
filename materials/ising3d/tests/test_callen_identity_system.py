#!/usr/bin/env python3
"""Independent exact verifier for the Callen identity-system artifact.

This verifier does not import the producer.  It reconstructs the six-neighbour
coefficients from raw local-spin interpolation, checks the complete local truth
table, rebuilds every finite-graph identity from raw spin sums, and independently
recomputes the modular rank lower bounds recorded in the artifact.
"""

from __future__ import annotations

import itertools
import json
import math
import platform
import resource
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "correlations" / "callen_identity_system.json"
RSS_CAP_BYTES = 2 * 1024**3
CPU_CAP_SECONDS = 120.0
FAILURES: list[str] = []

Poly = tuple[int, ...]


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def trim(values: Iterable[int]) -> Poly:
    out = [int(value) for value in values]
    if not out:
        return (0,)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def add(left: Poly, right: Poly) -> Poly:
    return trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    )


def multiply(left: Poly, right: Poly) -> Poly:
    if left == (0,) or right == (0,):
        return (0,)
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)
def poly_power(poly: Poly, exponent: int) -> Poly:
    result: Poly = (1,)
    base = poly
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        exponent >>= 1
        if exponent:
            base = multiply(base, base)
    return result




def scale(poly: Poly, scalar: int) -> Poly:
    return trim(scalar * value for value in poly)


def evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def evaluate_mod(poly: Poly, value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(poly):
        result = (result * value + coefficient) % prime
    return result


@dataclass(frozen=True)
class RawRationalFunction:
    """Uncancelled rational function; equality is exact cross multiplication."""

    numerator: Poly
    denominator: Poly = (1,)

    def __post_init__(self) -> None:
        numerator = trim(self.numerator)
        denominator = trim(self.denominator)
        if denominator == (0,):
            raise ZeroDivisionError("zero polynomial denominator")
        if denominator[-1] < 0:
            numerator = scale(numerator, -1)
            denominator = scale(denominator, -1)
        object.__setattr__(self, "numerator", numerator)
        object.__setattr__(self, "denominator", denominator)

    @staticmethod
    def zero() -> "RawRationalFunction":
        return RawRationalFunction((0,))

    def plus(self, other: "RawRationalFunction") -> "RawRationalFunction":
        return RawRationalFunction(
            add(
                multiply(self.numerator, other.denominator),
                multiply(other.numerator, self.denominator),
            ),
            multiply(self.denominator, other.denominator),
        )

    def times_fraction(self, scalar: Fraction | int) -> "RawRationalFunction":
        scalar = Fraction(scalar)
        return RawRationalFunction(
            scale(self.numerator, scalar.numerator),
            scale(self.denominator, scalar.denominator),
        )

    def times(self, other: "RawRationalFunction") -> "RawRationalFunction":
        return RawRationalFunction(
            multiply(self.numerator, other.numerator),
            multiply(self.denominator, other.denominator),
        )

    def equals(self, other: "RawRationalFunction") -> bool:
        return multiply(self.numerator, other.denominator) == multiply(
            other.numerator, self.denominator
        )

    def is_zero(self) -> bool:
        return self.numerator == (0,)


def function_from_record(record: dict[str, object]) -> RawRationalFunction:
    return RawRationalFunction(
        tuple(record["numerator_coefficients_ascending"]),
        tuple(record["denominator_coefficients_ascending"]),
    )


def tanh_symbolic(field: int) -> RawRationalFunction:
    if field == 0:
        return RawRationalFunction.zero()
    sign = 1 if field > 0 else -1
    degree = abs(field)
    numerator = [0] * (degree + 1)
    denominator = [0] * (degree + 1)
    for power in range(degree + 1):
        if power % 2:
            numerator[power] = sign * math.comb(degree, power)
        else:
            denominator[power] = math.comb(degree, power)
    return RawRationalFunction(trim(numerator), trim(denominator))


def tanh_fraction(field: int, value: Fraction) -> Fraction:
    function = tanh_symbolic(field)
    return evaluate(function.numerator, value) / evaluate(function.denominator, value)


def spin(state: int, site: int) -> int:
    return 1 if (state >> site) & 1 else -1


def elementary(spins: tuple[int, ...], degree: int) -> int:
    return sum(
        math.prod(spins[index] for index in subset)
        for subset in itertools.combinations(range(len(spins)), degree)
    )


def invert_integer_matrix(matrix: list[list[int]]) -> list[list[Fraction]]:
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


def determinant_three(matrix: list[list[int]]) -> int:
    first, second, third = matrix
    return (
        first[0] * (second[1] * third[2] - second[2] * third[1])
        - first[1] * (second[0] * third[2] - second[2] * third[0])
        + first[2] * (second[0] * third[1] - second[1] * third[0])
    )


def interpolate_six_coefficients() -> tuple[
    list[list[int]], list[list[Fraction]], dict[int, RawRationalFunction]
]:
    fields = (2, 4, 6)
    basis_degrees = (1, 3, 5)
    representatives: list[tuple[int, ...]] = []
    for field in fields:
        representative = next(
            tuple(spin(state, site) for site in range(6))
            for state in range(64)
            if sum(spin(state, site) for site in range(6)) == field
        )
        representatives.append(representative)
    matrix = [
        [elementary(spins, degree) for degree in basis_degrees]
        for spins in representatives
    ]
    inverse = invert_integer_matrix(matrix)
    field_functions = [tanh_symbolic(field) for field in fields]
    coefficients: dict[int, RawRationalFunction] = {}
    for degree, inverse_row in zip(basis_degrees, inverse):
        function = RawRationalFunction.zero()
        for scalar, field_function in zip(inverse_row, field_functions):
            function = function.plus(field_function.times_fraction(scalar))
        coefficients[degree] = function
    return matrix, inverse, coefficients


def raw_walsh_counts(degree: int, subset_size: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for state in range(1 << degree):
        spins = tuple(spin(state, site) for site in range(degree))
        field = sum(spins)
        character = math.prod(spins[:subset_size])
        counts[field] = counts.get(field, 0) + character
    return dict(sorted(counts.items()))


def numeric_local_coefficients(degree: int, value: Fraction) -> dict[int, Fraction]:
    coefficients: dict[int, Fraction] = {}
    for subset_size in range(1, degree + 1, 2):
        total = Fraction(0)
        for state in range(1 << degree):
            spins = tuple(spin(state, site) for site in range(degree))
            total += tanh_fraction(sum(spins), value) * math.prod(
                spins[:subset_size]
            )
        coefficients[subset_size] = total / (1 << degree)
    return coefficients


def raw_moment_numerators(
    vertex_count: int, edges: tuple[tuple[int, int], ...], value: Fraction
) -> tuple[list[Fraction], Fraction]:
    weights: list[Fraction] = []
    for state in range(1 << vertex_count):
        weight = Fraction(1)
        for left, right in edges:
            weight *= 1 + value * spin(state, left) * spin(state, right)
        weights.append(weight)
    moments: list[Fraction] = []
    for mask in range(1 << vertex_count):
        total = Fraction(0)
        sites = tuple(site for site in range(vertex_count) if (mask >> site) & 1)
        for state, weight in enumerate(weights):
            character = math.prod(spin(state, site) for site in sites)
            total += weight * character
        moments.append(total)
    return moments, sum(weights)


def adjacency(vertex_count: int, edges: tuple[tuple[int, int], ...]) -> tuple[tuple[int, ...], ...]:
    rows = [[] for _ in range(vertex_count)]
    for left, right in edges:
        rows[left].append(right)
        rows[right].append(left)
    return tuple(tuple(sorted(row)) for row in rows)


def subset_masks(neighbours: tuple[int, ...], size: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << site for site in subset)
        for subset in itertools.combinations(neighbours, size)
    )


def exact_graph_identity_audit(
    graph_record: dict[str, object], value: Fraction
) -> tuple[int, int, int, Fraction]:
    vertex_count = int(graph_record["vertex_count"])
    edges = tuple(tuple(edge) for edge in graph_record["edges"])
    neighbours = adjacency(vertex_count, edges)
    moments, partition = raw_moment_numerators(vertex_count, edges, value)
    degree_coefficients = {
        degree: numeric_local_coefficients(degree, value)
        for degree in {len(row) for row in neighbours}
    }
    zero_rows = 0
    negative_nonzero_rows = 0
    row_count = 0
    for vertex, local_neighbours in enumerate(neighbours):
        coefficients = degree_coefficients[len(local_neighbours)]
        masks_by_size = {
            size: subset_masks(local_neighbours, size) for size in coefficients
        }
        for observable in range(1 << vertex_count):
            if (observable >> vertex) & 1:
                continue
            row_count += 1
            residual = moments[observable ^ (1 << vertex)]
            for size, coefficient in coefficients.items():
                residual -= coefficient * sum(
                    moments[observable ^ mask] for mask in masks_by_size[size]
                )
            if residual == 0:
                zero_rows += 1
            mutated_residual = residual - sum(
                moments[observable ^ mask] for mask in masks_by_size[1]
            )
            if mutated_residual != 0:
                negative_nonzero_rows += 1
    return row_count, zero_rows, negative_nonzero_rows, partition


def tanh_mod(field: int, v_mod: int, prime: int) -> int:
    function = tanh_symbolic(field)
    denominator = evaluate_mod(function.denominator, v_mod, prime)
    if denominator == 0:
        raise ZeroDivisionError("bad local-field reduction")
    return evaluate_mod(function.numerator, v_mod, prime) * pow(
        denominator, -1, prime
    ) % prime


def modular_local_coefficients(
    degree: int, v_mod: int, prime: int
) -> dict[int, int]:
    inverse_states = pow(1 << degree, -1, prime)
    coefficients: dict[int, int] = {}
    for subset_size in range(1, degree + 1, 2):
        total = 0
        for state in range(1 << degree):
            spins = tuple(spin(state, site) for site in range(degree))
            character = math.prod(spins[:subset_size])
            total += character * tanh_mod(sum(spins), v_mod, prime)
        coefficients[subset_size] = total % prime * inverse_states % prime
    return coefficients


def update_mod(row: dict[int, int], column: int, value: int, prime: int) -> None:
    result = (row.get(column, 0) + value) % prime
    if result:
        row[column] = result
    else:
        row.pop(column, None)


def independent_modular_rows(
    vertex_count: int,
    neighbours: tuple[tuple[int, ...], ...],
    value: Fraction,
    prime: int,
    parity: int,
) -> list[dict[int, int]]:
    v_mod = value.numerator % prime * pow(value.denominator, -1, prime) % prime
    by_degree = {
        degree: modular_local_coefficients(degree, v_mod, prime)
        for degree in {len(row) for row in neighbours}
    }
    rows: list[dict[int, int]] = []
    for vertex, local_neighbours in enumerate(neighbours):
        coefficients = by_degree[len(local_neighbours)]
        masks_by_size = {
            size: subset_masks(local_neighbours, size) for size in coefficients
        }
        for observable in range(1 << vertex_count):
            if (observable >> vertex) & 1:
                continue
            if (observable.bit_count() + 1) % 2 != parity:
                continue
            row = {observable ^ (1 << vertex): 1}
            for size, coefficient in coefficients.items():
                for mask in masks_by_size[size]:
                    update_mod(row, observable ^ mask, -coefficient, prime)
            rows.append(row)
    return rows


def sparse_rank(rows: list[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for source in rows:
        row = {column: value % prime for column, value in source.items() if value % prime}
        while row:
            pivot = min(row)
            basis = pivots.get(pivot)
            if basis is None:
                inverse = pow(row[pivot], -1, prime)
                pivots[pivot] = {
                    column: value * inverse % prime
                    for column, value in row.items()
                    if value * inverse % prime
                }
                break
            factor = row[pivot]
            for column, value in basis.items():
                update_mod(row, column, -factor * value, prime)
    return len(pivots)


def rank_record_audit(graph_record: dict[str, object]) -> bool:
    vertex_count = int(graph_record["vertex_count"])
    edges = tuple(tuple(edge) for edge in graph_record["edges"])
    neighbours = adjacency(vertex_count, edges)
    rank_record = graph_record["rank"]
    columns = 1 << vertex_count
    parity_columns = columns // 2
    if (
        rank_record["claim_tag"] != "[COMPUTATION]"
        or rank_record["matrix_shape"]
        != [vertex_count * (1 << (vertex_count - 1)), columns]
        or rank_record["generic_rank"] != columns - 1
        or rank_record["generic_nullity"] != 1
        or rank_record["even_block_generic_rank"] != parity_columns - 1
        or rank_record["odd_block_generic_rank"] != parity_columns
        or not rank_record["rank_closed_exactly"]
    ):
        return False
    for sample in rank_record["specialization_lower_bounds"]:
        value = Fraction(sample["v"]["numerator"], sample["v"]["denominator"])
        prime = int(sample["prime"])
        even_rank = sparse_rank(
            independent_modular_rows(vertex_count, neighbours, value, prime, 0),
            prime,
        )
        odd_rank = sparse_rank(
            independent_modular_rows(vertex_count, neighbours, value, prime, 1),
            prime,
        )
        stored = sample["modular_lower_bounds"]
        ceilings = sample["matching_rational_ceilings"]
        pivots = sample["pivot_certificate"]
        if not (
            even_rank == stored["even_correlator_block"] == parity_columns - 1
            and odd_rank == stored["odd_correlator_block"] == parity_columns
            and even_rank + odd_rank == stored["full_system"] == columns - 1
            and ceilings
            == {
                "even_correlator_block": parity_columns - 1,
                "odd_correlator_block": parity_columns,
                "full_system": columns - 1,
            }
            and sample["exact_rank_at_this_rational_v"] == columns - 1
            and len(pivots["even_columns"]) == even_rank
            and len(pivots["odd_columns"]) == odd_rank
            and len(pivots["even_source_row_indices"]) == even_rank
            and len(pivots["odd_source_row_indices"]) == odd_rank
        ):
            return False
    return True


def collect_claim_tags(value: object) -> list[str]:
    tags: list[str] = []
    if isinstance(value, dict):
        if "claim_tag" in value:
            tags.append(str(value["claim_tag"]))
        for child in value.values():
            tags.extend(collect_claim_tags(child))
    elif isinstance(value, list):
        for child in value:
            tags.extend(collect_claim_tags(child))
    return tags


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard(started: float, stage: str) -> None:
    cpu = time.process_time() - started
    rss = peak_rss_bytes()
    if cpu > CPU_CAP_SECONDS:
        raise RuntimeError(f"CPU guard exceeded at {stage}: {cpu:.3f}s")
    if rss > RSS_CAP_BYTES:
        raise MemoryError(f"RSS guard exceeded at {stage}: {rss}")


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact has meta data checks", {"meta", "data", "checks"} <= set(artifact))
    names = [row["name"] for row in artifact["checks"]]
    check(
        "producer checks unique and true",
        len(names) == len(set(names)) and all(row["passed"] for row in artifact["checks"]),
    )
    check(
        "producer obeyed exact resource contract",
        artifact["meta"]["single_process"]
        and not artifact["meta"]["benchmark_used"]
        and artifact["meta"]["peak_rss_bytes"] < artifact["meta"]["rss_cap_bytes"]
        and artifact["meta"]["process_cpu_seconds"] < artifact["meta"]["cpu_cap_seconds"],
    )

    master = artifact["data"]["local_master_identity"]
    matrix, inverse, derived = interpolate_six_coefficients()
    stored_coefficients = {
        degree: function_from_record(master["coefficients"][f"c{degree}"])
        for degree in (1, 3, 5)
    }
    check(
        "raw orbit interpolation matrix rebuilt",
        matrix == [[2, -4, 2], [4, 0, -4], [6, 20, 6]]
        and determinant_three(matrix) == 512
        and master["interpolation_proof"]["orbit_matrix_rows_fields_2_4_6"] == matrix
        and master["interpolation_proof"]["orbit_matrix_determinant"] == 512,
    )
    stored_inverse = [
        [Fraction(value["numerator"], value["denominator"]) for value in row]
        for row in master["interpolation_proof"]["orbit_matrix_inverse"]
    ]
    check("orbit inverse independently solved", inverse == stored_inverse)
    check(
        "c1 c3 c5 independently rederived in Q(v)",
        all(derived[degree].equals(stored_coefficients[degree]) for degree in (1, 3, 5)),
    )

    rebuilt_counts = {
        f"e{degree}": {
            str(field): count
            for field, count in raw_walsh_counts(6, degree).items()
        }
        for degree in (1, 3, 5)
    }
    check(
        "six-spin interpolation uses only fields plus-minus 2 4 6",
        rebuilt_counts == master["interpolation_proof"]["walsh_field_counts"]
        and master["interpolation_proof"]["positive_local_fields_only"] == [2, 4, 6],
    )

    local_truth_ok = True
    for state in range(64):
        spins = tuple(spin(state, site) for site in range(6))
        rhs = RawRationalFunction.zero()
        for degree in (1, 3, 5):
            rhs = rhs.plus(derived[degree].times_fraction(elementary(spins, degree)))
        if not rhs.equals(tanh_symbolic(sum(spins))):
            local_truth_ok = False
            break
    check(
        "local master identity rebuilt on all raw spin states",
        local_truth_ok
        and master["interpolation_proof"]["all_spin_assignments_checked"] == 64
        and master["interpolation_proof"]["all_spin_residuals_zero"],
    )
    check(
        "local negative control is genuinely rejected",
        master["negative_control"]["rejected_on_every_positive_orbit"]
        and all(
            tuple(record["numerator_coefficients_ascending"]) != (0,)
            for record in master["negative_control"]["orbit_residuals"].values()
        ),
    )
    guard(started, "master identity")

    graph_records = artifact["data"]["finite_graph_computations"]
    expected_graphs = {
        "chain_3": (3, 2, [1, 1, 2]),
        "square_open_2x2": (4, 4, [2, 2, 2, 2]),
        "cube_open_2x2x2": (8, 12, [3] * 8),
    }
    check(
        "finite graph controls exactly declared",
        {
            row["name"]: (
                row["vertex_count"],
                row["edge_count"],
                row["degree_sequence"],
            )
            for row in graph_records
        }
        == expected_graphs,
    )
    for graph_record in graph_records:
        name = graph_record["name"]
        first_sample = graph_record["rank"]["specialization_lower_bounds"][0]
        value = Fraction(
            first_sample["v"]["numerator"], first_sample["v"]["denominator"]
        )
        row_count, zero_rows, negative_nonzero, raw_partition = exact_graph_identity_audit(
            graph_record, value
        )
        stored_enumeration = graph_record["exact_enumeration"]
        stored_partition_value = evaluate(
            tuple(stored_enumeration["partition_polynomial_coefficients_ascending"]),
            value,
        )
        check(
            f"{name} every Callen row rebuilt from raw spin sums",
            zero_rows == row_count == stored_enumeration["identity_rows_checked"]
            and stored_enumeration["symbolic_zero_rows"] == row_count
            and stored_enumeration["symbolic_nonzero_rows"] == 0
            and raw_partition == stored_partition_value,
            f"rows={zero_rows}/{row_count}",
        )
        check(
            f"{name} mutated coefficient negative control rebuilt",
            negative_nonzero
            == stored_enumeration["negative_control"]["nonzero_rows"]
            and negative_nonzero > 0
            and stored_enumeration["negative_control"]["rejected"],
            f"nonzero rows={negative_nonzero}",
        )
        check(
            f"{name} exact modular ranks independently rebuilt",
            rank_record_audit(graph_record),
        )
        guard(started, name)
    ratio_audit = True
    for field in range(7):
        function = tanh_symbolic(field)
        conditional_ratio = RawRationalFunction(
            add(function.denominator, function.numerator),
            add(function.denominator, scale(function.numerator, -1)),
        )
        expected_ratio = RawRationalFunction(
            poly_power((1, 1), field),
            poly_power((1, -1), field),
        )
        ratio_audit &= conditional_ratio.equals(expected_ratio)
    check(
        "Callen conditional relation gives exact Gibbs edge ratios",
        ratio_audit,
    )
    walsh_size = 16
    walsh_invertible = all(
        sum(
            (-1) ** ((left & state).bit_count() + (right & state).bit_count())
            for state in range(walsh_size)
        )
        == (walsh_size if left == right else 0)
        for left in range(walsh_size)
        for right in range(walsh_size)
    )
    reconstruction = artifact["data"]["finite_all_graph_reconstruction"]
    check(
        "all-finite-graph reconstruction theorem has load-bearing steps",
        walsh_invertible
        and reconstruction["claim_tag"] == "[THEOREM]"
        and len(reconstruction["proof"]) == 4
        and "2^n-1" in reconstruction["statement"]
        and "does not compress" in reconstruction["complexity_scope"],
    )
    k2_neighbours = ((1,), (0,))
    prime = 1_000_003
    k2_ranks = {
        str(value): sum(
            sparse_rank(
                independent_modular_rows(
                    2, k2_neighbours, value, prime, parity
                ),
                prime,
            )
            for parity in (0, 1)
        )
        for value in (Fraction(-1), Fraction(1, 2), Fraction(1))
    }
    check(
        "K2 endpoint control enforces specialization scope",
        k2_ranks
        == {"-1": 2, "1/2": 3, "1": 2}
        and "v0!=+/-1" in reconstruction["field_scope"]
        and "Z_G(v0)!=0" in reconstruction["field_scope"],
        str(k2_ranks),
    )
    check(
        "all-graph theorem declares interaction and boundary hypotheses",
        all(
            phrase in reconstruction["statement"]
            for phrase in (
                "loopless undirected",
                "uniform coupling",
                "free boundary",
                "no external field",
            )
        )
        and "Fixed boundary spins" in reconstruction["model_scope"],
    )

    tags = collect_claim_tags(artifact["data"])
    check(
        "theorem tags and finite-computation scope are separated",
        tags.count("[THEOREM]") == 2
        and master["claim_tag"] == "[THEOREM]"
        and reconstruction["claim_tag"] == "[THEOREM]"
        and artifact["data"]["callen_system_definition"]["claim_tag"]
        == "[EXACT IDENTITY]"
        and artifact["data"]["finite_rank_conclusion"]["claim_tag"]
        == "[COMPUTATION]",
        str(tags),
    )
    limitations = artifact["data"]["scope"]["not_claimed"]
    check(
        "scope excludes compression infinite volume and criticality",
        any("subexponential" in item for item in limitations)
        and any("infinite-volume" in item for item in limitations)
        and any("criticality" in item for item in limitations),
    )
    check(
        "verifier resource guard respected",
        time.process_time() - started < CPU_CAP_SECONDS
        and peak_rss_bytes() < RSS_CAP_BYTES,
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks")
        return 1
    print(f"PASS: independent Callen verifier ({time.process_time() - started:.3f}s CPU)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
