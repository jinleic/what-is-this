#!/usr/bin/env python3
"""Clean-room verifier for the selected-row Callen Walsh certificate.

This file deliberately imports neither e246 nor e240.  It rebuilds the Callen
coefficients as rational functions from the defining Walsh sum, specializes at
v=1/3, reconstructs every raw and folded matrix, and checks all source, scope,
and resource gates in the stored artifact.
"""

from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "correlations" / "callen_walsh_compression.json"
PRODUCER = ROOT / "experiments" / "e246_callen_walsh_compression.py"
E240_PRODUCER = ROOT / "experiments" / "e240_callen_identity_system.py"
E243_PRODUCER = ROOT / "experiments" / "e243_callen_termwise_closure.py"
RSS_LIMIT_BYTES = 2 * 1024**3
V = Fraction(1, 3)
DIRECTIONS = ("+x", "-x", "+y", "-y", "+z", "-z")
INVERSION = (1, 0, 3, 2, 5, 4)
EXPECTED_SCOPE = (
    "six exact inversion-invariant higher-template relations and the all-even-degree "
    "obstruction to pair-row leakage cancellation in the fixed selected-pivot/far-mark "
    "Callen block; no susceptibility closure and no contact-inclusive or radius-expanded no-go"
)
EXPECTED_VALUES = {
    2: {1: Fraction(3, 10)},
    4: {1: Fraction(177, 680), 3: Fraction(-27, 680)},
    6: {
        1: Fraction(4143, 17680),
        3: Fraction(-27, 1040),
        5: Fraction(243, 17680),
    },
}
FAILURES: list[str] = []

Poly = tuple[Fraction, ...]
Matrix = list[list[Fraction]]


def check(name: str, passed: bool, detail: str = "") -> None:
    suffix = f": {detail}" if detail else ""
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}{suffix}")
    if not passed:
        FAILURES.append(name)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        class TimeValue(ctypes.Structure):
            _fields_ = [
                ("seconds", ctypes.c_int32),
                ("microseconds", ctypes.c_int32),
            ]

        class MachTaskBasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_time", TimeValue),
                ("system_time", TimeValue),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        system = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        system.mach_task_self.restype = ctypes.c_uint32
        info = MachTaskBasicInfo()
        count = ctypes.c_uint32(
            ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_int32)
        )
        result = system.task_info(
            system.mach_task_self(),
            20,
            ctypes.byref(info),
            ctypes.byref(count),
        )
        if result != 0:
            raise OSError(f"mach task_info failed with status {result}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def peak_rss_measurement() -> str:
    if platform.system() == "Darwin":
        return "mach_task_basic_info.resident_size_max (task_info flavor 20)"
    return "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"


def fraction_record(value: Fraction | int) -> dict[str, int]:
    value = Fraction(value)
    return {"numerator": value.numerator, "denominator": value.denominator}


def read_fraction(record: dict[str, object]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def trim(poly: Poly) -> Poly:
    end = len(poly)
    while end > 1 and poly[end - 1] == 0:
        end -= 1
    return poly[:end] if end else (Fraction(0),)


def poly_add(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return trim(
        tuple(
            (left[index] if index < len(left) else Fraction(0))
            + (right[index] if index < len(right) else Fraction(0))
            for index in range(size)
        )
    )


def poly_scale(poly: Poly, scalar: Fraction) -> Poly:
    return trim(tuple(scalar * coefficient for coefficient in poly))


def poly_multiply(left: Poly, right: Poly) -> Poly:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_power, left_value in enumerate(left):
        for right_power, right_value in enumerate(right):
            result[left_power + right_power] += left_value * right_value
    return trim(tuple(result))


def poly_evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def tanh_polynomials(field: int) -> tuple[Poly, Poly]:
    if field == 0:
        return (Fraction(0),), (Fraction(1),)
    degree = abs(field)
    sign = 1 if field > 0 else -1
    numerator = tuple(
        Fraction(sign * math.comb(degree, power) if power % 2 else 0)
        for power in range(degree + 1)
    )
    denominator = tuple(
        Fraction(math.comb(degree, power) if power % 2 == 0 else 0)
        for power in range(degree + 1)
    )
    return trim(numerator), trim(denominator)


def field_multiplicities(degree: int, subset_size: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for state in itertools.product((-1, 1), repeat=degree):
        field = sum(state)
        sign = math.prod(state[:subset_size])
        counts[field] = counts.get(field, 0) + sign
    return dict(sorted(counts.items()))


def reconstructed_coefficient(
    degree: int, subset_size: int
) -> tuple[Poly, Poly, dict[int, int]]:
    numerator: Poly = (Fraction(0),)
    denominator: Poly = (Fraction(1),)
    counts = field_multiplicities(degree, subset_size)
    for field, multiplicity in counts.items():
        term_numerator, term_denominator = tanh_polynomials(field)
        term_numerator = poly_scale(
            term_numerator, Fraction(multiplicity, 1 << degree)
        )
        numerator = poly_add(
            poly_multiply(numerator, term_denominator),
            poly_multiply(term_numerator, denominator),
        )
        denominator = poly_multiply(denominator, term_denominator)
    return numerator, denominator, counts


def formal_coefficient_matches(
    stored: dict[str, object], numerator: Poly, denominator: Poly
) -> bool:
    rational = stored["rational_function"]
    stored_numerator = tuple(
        Fraction(int(value))
        for value in rational["numerator_coefficients_ascending"]
    )
    stored_denominator = tuple(
        Fraction(int(value))
        for value in rational["denominator_coefficients_ascending"]
    )
    return trim(poly_multiply(numerator, stored_denominator)) == trim(
        poly_multiply(stored_numerator, denominator)
    )


def tanh_multiple(field: int, value: Fraction) -> Fraction:
    numerator, denominator = tanh_polynomials(field)
    return poly_evaluate(numerator, value) / poly_evaluate(denominator, value)


def matrix_rank(matrix: Matrix) -> int:
    if not matrix:
        return 0
    reduced = [row[:] for row in matrix]
    rows = len(reduced)
    columns = len(reduced[0])
    rank = 0
    for column in range(columns):
        pivot = None
        for candidate in range(rank, rows):
            if reduced[candidate][column] != 0:
                pivot = candidate
                break
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        divisor = reduced[rank][column]
        for index in range(column, columns):
            reduced[rank][index] /= divisor
        for candidate in range(rank + 1, rows):
            factor = reduced[candidate][column]
            if factor == 0:
                continue
            for index in range(column, columns):
                reduced[candidate][index] -= factor * reduced[rank][index]
        rank += 1
        if rank == rows:
            break
    return rank


def matrix_record(matrix: Matrix, column_count: int) -> dict[str, object]:
    rank = matrix_rank(matrix)
    shape = [len(matrix), column_count]
    entries = [[fraction_record(value) for value in row] for row in matrix]
    return {
        "shape": shape,
        "rank": rank,
        "left_nullity": len(matrix) - rank,
        "sha256": canonical_sha256({"shape": shape, "entries": entries}),
    }


def character(mask: int, query: int) -> int:
    return 1 if (mask & query).bit_count() % 2 == 0 else -1


def masks_with_parity(degree: int, parity: int) -> list[int]:
    return [
        mask
        for mask in range(1 << degree)
        if mask.bit_count() % 2 == parity
    ]


def rebuild_matrix(
    rows: list[int], columns: list[int], coefficients: dict[int, Fraction]
) -> Matrix:
    result: Matrix = []
    for row_mask in rows:
        result.append(
            [
                coefficients[(row_mask ^ column_mask).bit_count()]
                for column_mask in columns
            ]
        )
    return result


def rebuild_control(
    degree: int, coefficients: dict[int, Fraction]
) -> tuple[dict[str, object], bool]:
    full_mask = (1 << degree) - 1
    even = masks_with_parity(degree, 0)
    odd = masks_with_parity(degree, 1)
    selected = [mask for mask in even if mask != full_mask]
    no_empty = [mask for mask in selected if mask != 0]
    full_matrix = rebuild_matrix(even, odd, coefficients)
    selected_matrix = rebuild_matrix(selected, odd, coefficients)
    no_empty_matrix = rebuild_matrix(no_empty, odd, coefficients)

    m = degree // 2
    middle_pairs = math.comb(degree, m) // 2
    rank_formula = (1 << (degree - 1)) - middle_pairs
    selected_nullity = middle_pairs - 1
    spectrum_holds = True
    layers = []
    for query_size in range(degree + 1):
        expected = tanh_multiple(degree - 2 * query_size, V)
        observed = []
        for query in range(1 << degree):
            if query.bit_count() == query_size:
                observed.append(
                    sum(
                        coefficients[mask.bit_count()] * character(mask, query)
                        for mask in odd
                    )
                )
        spectrum_holds &= len(observed) == math.comb(degree, query_size)
        spectrum_holds &= all(value == expected for value in observed)
        layers.append(
            {
                "query_size": query_size,
                "multiplicity": len(observed),
                "eigenvalue": fraction_record(expected),
            }
        )
    matrices = {
        "full_even_convolution": matrix_record(full_matrix, len(odd)),
        "selected_proper_even": matrix_record(selected_matrix, len(odd)),
        "selected_without_empty": matrix_record(no_empty_matrix, len(odd)),
    }
    ranks_hold = (
        matrices["full_even_convolution"]["rank"] == rank_formula
        and matrices["selected_proper_even"]["rank"] == rank_formula
        and matrices["selected_proper_even"]["left_nullity"]
        == selected_nullity
        and matrices["selected_without_empty"]["rank"] == rank_formula - 1
        and matrices["selected_without_empty"]["left_nullity"]
        == selected_nullity
    )
    return (
        {
            "degree": degree,
            "m": m,
            "row_set": "all even U except the full neighbour set N",
            "column_set": "all odd W",
            "predicted_rank": rank_formula,
            "predicted_selected_left_nullity": selected_nullity,
            "middle_layer_size": math.comb(degree, m),
            "walsh_layers": layers,
            "matrices": matrices,
        },
        spectrum_holds and ranks_hold,
    )


def invert(mask: int) -> int:
    image = 0
    for source in range(6):
        if mask & (1 << source):
            image |= 1 << INVERSION[source]
    return image


def partition_under_inversion(masks: list[int]) -> list[tuple[int, ...]]:
    unused = set(masks)
    classes = []
    while unused:
        representative = min(unused)
        members = tuple(sorted(set((representative, invert(representative)))))
        if not set(members).issubset(unused):
            raise AssertionError("invalid inversion class")
        unused -= set(members)
        classes.append(members)
    return classes


def names(mask: int) -> list[str]:
    return [DIRECTIONS[index] for index in range(6) if mask & (1 << index)]


def orbit_record(orbit: tuple[int, ...]) -> dict[str, object]:
    return {"masks": list(orbit), "subsets": [names(mask) for mask in orbit]}


def pair_orbits() -> list[tuple[tuple[int, int], ...]]:
    all_mask = 63
    pair_classes = {
        tuple(sorted((mask, all_mask ^ mask)))
        for mask in range(64)
        if mask.bit_count() == 3
    }
    unused = set(pair_classes)
    classes = []
    while unused:
        representative = min(unused)
        image_mask = invert(representative[0])
        image = tuple(sorted((image_mask, all_mask ^ image_mask)))
        members = tuple(sorted(set((representative, image))))
        if not set(members).issubset(unused):
            raise AssertionError("invalid complement-pair class")
        unused -= set(members)
        classes.append(members)
    classes.sort(key=lambda item: (len(item), item))
    return classes


def pair_orbit_record(orbit: tuple[tuple[int, int], ...]) -> dict[str, object]:
    return {
        "weight": len(orbit),
        "complement_pairs": [
            {
                "masks": list(pair),
                "subsets": [names(mask) for mask in pair],
            }
            for pair in orbit
        ],
    }


def kernel_value(
    alpha: list[int], mask: int, central_orbits: list[tuple[tuple[int, int], ...]]
) -> int:
    value = 0
    for coefficient, central_orbit in zip(alpha, central_orbits):
        value += coefficient * sum(
            character(pair[0], mask) for pair in central_orbit
        )
    return value


def rebuild_folded(coefficients: dict[int, Fraction]) -> tuple[dict[str, object], bool]:
    rows = [mask for mask in masks_with_parity(6, 0) if mask != 63]
    columns = masks_with_parity(6, 1)
    row_classes = partition_under_inversion(rows)
    column_classes = partition_under_inversion(columns)

    folded: Matrix = []
    invariant = True
    for row_class in row_classes:
        entries = []
        for column_class in column_classes:
            candidate_values = []
            for column in column_class:
                candidate_values.append(
                    sum(
                        coefficients[(row ^ column).bit_count()]
                        for row in row_class
                    )
                )
            invariant &= all(value == candidate_values[0] for value in candidate_values)
            entries.append(candidate_values[0])
        folded.append(entries)
    without_empty = [
        row for row_class, row in zip(row_classes, folded) if row_class != (0,)
    ]
    selected_record = matrix_record(folded, len(column_classes))
    without_empty_record = matrix_record(without_empty, len(column_classes))

    central_orbits = pair_orbits()
    weights = [len(orbit) for orbit in central_orbits]
    vectors = []
    vector_checks = True
    for relation_index in range(1, 7):
        alpha = [0] * 7
        alpha[0] = -weights[relation_index]
        alpha[relation_index] = weights[0]
        row_values = [
            kernel_value(alpha, row_class[0], central_orbits)
            for row_class in row_classes
        ]
        class_invariance = all(
            kernel_value(alpha, member, central_orbits) == row_value
            for row_class, row_value in zip(row_classes, row_values)
            for member in row_class
        )
        residual = [
            sum(
                Fraction(row_values[row_index]) * folded[row_index][column_index]
                for row_index in range(len(row_classes))
            )
            for column_index in range(len(column_classes))
        ]
        weighted_sum = sum(
            weight * coefficient for weight, coefficient in zip(weights, alpha)
        )
        lambda_empty = kernel_value(alpha, 0, central_orbits)
        lambda_full = kernel_value(alpha, 63, central_orbits)
        vector_checks &= (
            class_invariance
            and weighted_sum == 0
            and lambda_empty == 0
            and lambda_full == 0
            and not any(residual)
        )
        vectors.append(
            {
                "name": f"higher_relation_{relation_index}",
                "alpha": alpha,
                "weighted_alpha_sum": weighted_sum,
                "folded_row_values": row_values,
                "lambda_empty": lambda_empty,
                "lambda_full": lambda_full,
                "leakage_residual": [fraction_record(value) for value in residual],
            }
        )
    vector_rank = matrix_rank(
        [
            [Fraction(value) for value in vector["folded_row_values"]]
            for vector in vectors
        ]
    )
    result = {
        "inversion": {
            "direction_order": list(DIRECTIONS),
            "permutation": list(INVERSION),
            "definition": "the single cubic inversion simultaneously swaps +e_i and -e_i for i=x,y,z",
        },
        "row_classes": [orbit_record(orbit) for orbit in row_classes],
        "column_classes": [orbit_record(orbit) for orbit in column_classes],
        "folding_convention": "one row variable per inversion orbit; the folded row is the sum over raw rows and one representative of each invariant column equation",
        "matrices": {
            "selected_proper_even": selected_record,
            "selected_without_empty": without_empty_record,
        },
        "central_complement_pair_orbits": [
            pair_orbit_record(orbit) for orbit in central_orbits
        ],
        "weighted_alpha_constraint": "alpha1+alpha2+alpha3+alpha4+2*alpha5+2*alpha6+2*alpha7=0",
        "weights": weights,
        "kernel_vectors": vectors,
        "kernel_vector_rank": vector_rank,
    }
    passed = (
        invariant
        and len(row_classes) == 19
        and len(column_classes) == 16
        and selected_record["shape"] == [19, 16]
        and selected_record["rank"] == 13
        and selected_record["left_nullity"] == 6
        and without_empty_record["shape"] == [18, 16]
        and without_empty_record["rank"] == 12
        and without_empty_record["left_nullity"] == 6
        and weights == [1, 1, 1, 1, 2, 2, 2]
        and vector_checks
        and vector_rank == 6
        and all(vector["folded_row_values"][0] == 0 for vector in vectors)
    )
    return result, passed


def rebuild_far_marks() -> dict[str, object]:
    shape = (3, 3, 3)
    pivot = (0, 0, 0)
    displacement = (
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    )
    sites = {
        (x, y, z)
        for x in range(shape[0])
        for y in range(shape[1])
        for z in range(shape[2])
    }
    neighbours = tuple(
        tuple((pivot[axis] + delta[axis]) % shape[axis] for axis in range(3))
        for delta in displacement
    )
    forbidden = set(neighbours)
    forbidden.add(pivot)
    far = sorted(sites.difference(forbidden))
    pivot_contacts = len([site for site in far if site == pivot])
    neighbour_contacts = len([site for site in far if site in neighbours])
    return {
        "shape": list(shape),
        "pivot": list(pivot),
        "direction_sites": {
            direction: list(site)
            for direction, site in zip(DIRECTIONS, neighbours)
        },
        "site_count": len(sites),
        "distinct_neighbour_count": len(set(neighbours)),
        "closed_neighbourhood_count": len(forbidden),
        "far_mark_count": len(far),
        "far_marks": [list(site) for site in far],
        "far_mark_pivot_contacts": pivot_contacts,
        "far_mark_neighbour_contacts": neighbour_contacts,
        "no_contact_terms": pivot_contacts == 0 and neighbour_contacts == 0,
        "contact_definition": "a contact occurs only when the marked site equals the selected pivot or one of its six neighbour variables",
    }


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    meta = artifact["meta"]
    data = artifact["data"]
    certificate = data["certificate"]

    expected_source_hashes = {
        "experiments/e246_callen_walsh_compression.py": file_sha256(PRODUCER),
        "tests/test_callen_walsh_compression.py": file_sha256(Path(__file__).resolve()),
        "experiments/e240_callen_identity_system.py": file_sha256(E240_PRODUCER),
        "experiments/e243_callen_termwise_closure.py": file_sha256(E243_PRODUCER),
    }
    check(
        "current source hash gate",
        meta["source_sha256"] == expected_source_hashes,
        "producer, verifier, e240, and e243 must all match the generated artifact",
    )
    check(
        "inherited producer hash bindings",
        data["inherited_sources"]
        == {
            "e240_coefficients": {
                "path": "experiments/e240_callen_identity_system.py",
                "sha256": expected_source_hashes[
                    "experiments/e240_callen_identity_system.py"
                ],
            },
            "e243_scope_predecessor": {
                "path": "experiments/e243_callen_termwise_closure.py",
                "sha256": expected_source_hashes[
                    "experiments/e243_callen_termwise_closure.py"
                ],
            },
        },
    )
    check(
        "certificate content hash",
        meta["certificate_sha256"] == canonical_sha256(certificate),
    )

    check_names = [str(row["name"]) for row in artifact["checks"]]
    required_checks = {
        "source_hashes_bound",
        "exact_degree_2_4_6_coefficients",
        "finite_walsh_spectra_and_rank_formula",
        "corrected_raw_d6_rank_criterion",
        "cubic_inversion_folded_ranks",
        "six_weighted_folded_kernel_vectors",
        "periodic_3x3x3_far_marks",
        "scope_gate",
        "benchmark_not_used",
        "declared_resource_limits",
    }
    check(
        "producer checks are hard gates",
        len(check_names) == len(set(check_names))
        and set(check_names) == required_checks
        and all(row["passed"] is True for row in artifact["checks"]),
    )
    check(
        "producer resource and benchmark gates",
        meta["benchmark_used"] is False
        and certificate["parameter"]["benchmark_used"] is False
        and meta["process_cpu_seconds"] < meta["process_cpu_budget_seconds"]
        and meta["peak_rss_bytes"] < meta["rss_limit_bytes"]
        and meta["rss_limit_bytes"] == RSS_LIMIT_BYTES
        and meta["peak_rss_measurement"] == peak_rss_measurement()
        and meta["arithmetic"]
        == "exact integers, fractions, and rational functions over Z[v]",
    )

    rebuilt_values: dict[int, dict[int, Fraction]] = {}
    coefficients_hold = True
    coefficient_payload = certificate["coefficients"]
    coefficients_hold &= set(coefficient_payload) == {"2", "4", "6"}
    for degree in (2, 4, 6):
        stored_degree = coefficient_payload[str(degree)]
        expected_sizes = set(range(1, degree + 1, 2))
        coefficients_hold &= {int(size) for size in stored_degree} == expected_sizes
        rebuilt_values[degree] = {}
        for subset_size in sorted(expected_sizes):
            numerator, denominator, counts = reconstructed_coefficient(
                degree, subset_size
            )
            value = poly_evaluate(numerator, V) / poly_evaluate(denominator, V)
            stored = stored_degree[str(subset_size)]
            rebuilt_values[degree][subset_size] = value
            coefficients_hold &= formal_coefficient_matches(
                stored, numerator, denominator
            )
            coefficients_hold &= stored["walsh_field_multiplicities"] == {
                str(field): multiplicity for field, multiplicity in counts.items()
            }
            coefficients_hold &= read_fraction(stored["value_at_v"]) == value
        coefficients_hold &= rebuilt_values[degree] == EXPECTED_VALUES[degree]
    check(
        "independent exact coefficient reconstruction",
        coefficients_hold,
        "formal Walsh sums over Q(v), then the exact v=1/3 specialization",
    )

    controls_hold = True
    rebuilt_controls = {}
    for degree in (2, 4, 6):
        control, internal_holds = rebuild_control(degree, rebuilt_values[degree])
        rebuilt_controls[str(degree)] = control
        controls_hold &= internal_holds
    controls_hold &= rebuilt_controls == certificate["finite_controls"]
    check(
        "independent d=2,4,6 Walsh and matrix controls",
        controls_hold,
    )

    d6_raw = rebuilt_controls["6"]["matrices"]
    check(
        "corrected d6 raw ranks and nullities",
        d6_raw["selected_proper_even"]["shape"] == [31, 32]
        and d6_raw["selected_proper_even"]["rank"] == 22
        and d6_raw["selected_proper_even"]["left_nullity"] == 9
        and d6_raw["selected_without_empty"]["shape"] == [30, 32]
        and d6_raw["selected_without_empty"]["rank"] == 21
        and d6_raw["selected_without_empty"]["left_nullity"] == 9,
    )

    folded, folded_holds = rebuild_folded(rebuilt_values[6])
    check(
        "independent folded matrices and six kernels",
        folded_holds and folded == certificate["cubic_inversion"],
        "19x16 rank 13; 18x16 rank 12; six-dimensional explicit left kernel",
    )
    check(
        "weighted alpha and pair-row obstruction",
        folded["weighted_alpha_constraint"]
        == "alpha1+alpha2+alpha3+alpha4+2*alpha5+2*alpha6+2*alpha7=0"
        and folded["weights"] == [1, 1, 1, 1, 2, 2, 2]
        and all(
            vector["weighted_alpha_sum"] == 0
            and vector["lambda_empty"] == 0
            and vector["lambda_full"] == 0
            and not any(read_fraction(value) for value in vector["leakage_residual"])
            for vector in folded["kernel_vectors"]
        ),
    )

    far_marks = rebuild_far_marks()
    check(
        "independent 3x3x3 far-mark audit",
        far_marks == certificate["periodic_3x3x3"]
        and far_marks["far_mark_count"] == 20
        and far_marks["no_contact_terms"] is True,
    )

    expected_theorem = {
        "tag": "[THEOREM]",
        "assumptions": "d=2m with integer m>=1 and 0<v=tanh(K)<1",
        "row_set": "E^o={even U subset N: U!=N}",
        "column_set": "O={odd W subset N}",
        "matrix": "L_d(U,W)=c_|U symmetric_difference W|(v)",
        "walsh_eigenvalue": "tanh(K*(d-2*|Q|))",
        "zero_criterion": "the eigenvalue is zero exactly when |Q|=m",
        "rank": "2^(d-1)-binom(d,m)/2",
        "selected_left_nullity": "binom(d,m)/2-1",
        "fourier_symmetry": "lambda_hat(Q)=lambda_hat(Q^c) on even coordinates",
        "endpoint_identity": "lambda(N)=(-1)^m*lambda(empty)",
        "pair_row_obstruction": "extending a selected-row left-kernel vector by lambda(N)=0 forces lambda(empty)=0",
    }
    check("general even-degree theorem schema", data["theorem"] == expected_theorem)

    expected_scope = {
        "tag": "[SCOPE]",
        "claim": EXPECTED_SCOPE,
        "proved": [
            "six exact inversion-invariant higher-template relations in the d=6 folded block",
            "the all-even-degree pair-row obstruction for the fixed selected-pivot/far-mark block",
            "the contact-free count of 20 far marks for one pivot on the periodic 3x3x3 cube",
        ],
        "not_proved": [
            "pair or susceptibility closure",
            "a no-go for contact-inclusive selected-row systems",
            "a no-go after expanding the radius, pivot set, or observable family",
            "an exact critical coupling or thermodynamic solution",
        ],
    }
    check(
        "scope is a hard gate",
        data["scope"] == expected_scope
        and data["scope"]["claim"] == EXPECTED_SCOPE,
    )

    current_rss = peak_rss_bytes()
    check(
        "verifier RSS wall",
        current_rss < RSS_LIMIT_BYTES,
        f"{current_rss}/{RSS_LIMIT_BYTES} bytes via {peak_rss_measurement()}",
    )
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: selected-row Callen Walsh compression")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
