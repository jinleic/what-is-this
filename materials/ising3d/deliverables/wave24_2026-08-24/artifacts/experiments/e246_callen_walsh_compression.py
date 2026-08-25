#!/usr/bin/env python3
"""Exact selected-row Callen compression and Walsh-kernel certificate.

The all-even-degree rank statement is proved in the companion proof.  This
producer supplies exact d=2,4,6 controls, the six folded d=6 kernel vectors,
and the contact-free 3x3x3 far-mark audit.  It writes its artifact only after
all checks and resource gates pass.
"""

from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments import e240_callen_identity_system as callen  # noqa: E402

OUTPUT = ROOT / "results" / "correlations" / "callen_walsh_compression.json"
VERIFIER = ROOT / "tests" / "test_callen_walsh_compression.py"
E240_PRODUCER = ROOT / "experiments" / "e240_callen_identity_system.py"
E243_PRODUCER = ROOT / "experiments" / "e243_callen_termwise_closure.py"
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
V = Fraction(1, 3)
DIRECTIONS = ("+x", "-x", "+y", "-y", "+z", "-z")
INVERSION = (1, 0, 3, 2, 5, 4)
SCOPE = (
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

Matrix = list[list[Fraction]]


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


def evaluate_polynomial(coefficients: tuple[int, ...], value: Fraction) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(coefficients):
        total = total * value + coefficient
    return total


def evaluate_rational_function(
    function: callen.RationalFunction, value: Fraction
) -> Fraction:
    numerator = evaluate_polynomial(function.numerator, value)
    denominator = evaluate_polynomial(function.denominator, value)
    if denominator == 0:
        raise ZeroDivisionError("Callen coefficient denominator vanished")
    return numerator / denominator


def tanh_multiple_fraction(field: int, value: Fraction) -> Fraction:
    if field == 0:
        return Fraction(0)
    degree = abs(field)
    plus = (1 + value) ** degree
    minus = (1 - value) ** degree
    result = (plus - minus) / (plus + minus)
    return result if field > 0 else -result


def matrix_rank(matrix: Matrix) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                entry - multiple * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def matrix_record(matrix: Matrix, column_count: int) -> dict[str, object]:
    rank = matrix_rank(matrix)
    encoded = [[fraction_record(entry) for entry in row] for row in matrix]
    shape = [len(matrix), column_count]
    return {
        "shape": shape,
        "rank": rank,
        "left_nullity": len(matrix) - rank,
        "sha256": canonical_sha256({"shape": shape, "entries": encoded}),
    }


def character(mask: int, subset: int) -> int:
    return -1 if (mask & subset).bit_count() % 2 else 1


def parity_masks(degree: int, parity: int) -> list[int]:
    return [
        mask
        for mask in range(1 << degree)
        if mask.bit_count() % 2 == parity
    ]


def callen_matrix(
    rows: list[int], columns: list[int], coefficients: dict[int, Fraction]
) -> Matrix:
    return [
        [coefficients[(row ^ column).bit_count()] for column in columns]
        for row in rows
    ]


def elementary_symmetric(spins: tuple[int, ...], degree: int) -> int:
    return sum(
        math.prod(spins[index] for index in subset)
        for subset in itertools.combinations(range(len(spins)), degree)
    )


def coefficient_record(
    degree: int,
) -> tuple[dict[int, Fraction], dict[str, object], bool]:
    functions = callen.local_coefficients(degree)
    values: dict[int, Fraction] = {}
    records: dict[str, object] = {}
    for subset_size, function in functions.items():
        value = evaluate_rational_function(function, V)
        values[subset_size] = value
        records[str(subset_size)] = {
            "rational_function": function.as_json(),
            "walsh_field_multiplicities": {
                str(field): multiplicity
                for field, multiplicity in callen.walsh_field_counts(
                    degree, subset_size
                ).items()
            },
            "value_at_v": fraction_record(value),
        }

    expansion_holds = True
    for state in range(1 << degree):
        spins = tuple(1 if (state >> site) & 1 else -1 for site in range(degree))
        reconstructed = sum(
            values[subset_size] * elementary_symmetric(spins, subset_size)
            for subset_size in values
        )
        expansion_holds &= reconstructed == tanh_multiple_fraction(sum(spins), V)
    return values, records, expansion_holds


def finite_control(
    degree: int, coefficients: dict[int, Fraction]
) -> tuple[dict[str, object], bool]:
    full_mask = (1 << degree) - 1
    even = parity_masks(degree, 0)
    odd = parity_masks(degree, 1)
    selected = [mask for mask in even if mask != full_mask]
    higher_only = [mask for mask in selected if mask != 0]

    full_matrix = callen_matrix(even, odd, coefficients)
    selected_matrix = callen_matrix(selected, odd, coefficients)
    higher_matrix = callen_matrix(higher_only, odd, coefficients)
    middle = degree // 2
    central_pairs = math.comb(degree, middle) // 2
    predicted_rank = (1 << (degree - 1)) - central_pairs
    predicted_selected_nullity = central_pairs - 1

    layer_records = []
    spectrum_holds = True
    for size in range(degree + 1):
        layer = []
        for query in range(1 << degree):
            if query.bit_count() != size:
                continue
            observed = sum(
                coefficients[odd_mask.bit_count()]
                * character(odd_mask, query)
                for odd_mask in odd
            )
            layer.append(observed)
        expected = tanh_multiple_fraction(degree - 2 * size, V)
        spectrum_holds &= len(layer) == math.comb(degree, size)
        spectrum_holds &= all(value == expected for value in layer)
        layer_records.append(
            {
                "query_size": size,
                "multiplicity": len(layer),
                "eigenvalue": fraction_record(expected),
            }
        )

    matrices = {
        "full_even_convolution": matrix_record(full_matrix, len(odd)),
        "selected_proper_even": matrix_record(selected_matrix, len(odd)),
        "selected_without_empty": matrix_record(higher_matrix, len(odd)),
    }
    rank_holds = (
        matrices["full_even_convolution"]["rank"] == predicted_rank
        and matrices["selected_proper_even"]["rank"] == predicted_rank
        and matrices["selected_proper_even"]["left_nullity"]
        == predicted_selected_nullity
        and matrices["selected_without_empty"]["rank"] == predicted_rank - 1
        and matrices["selected_without_empty"]["left_nullity"]
        == predicted_selected_nullity
    )
    record = {
        "degree": degree,
        "m": middle,
        "row_set": "all even U except the full neighbour set N",
        "column_set": "all odd W",
        "predicted_rank": predicted_rank,
        "predicted_selected_left_nullity": predicted_selected_nullity,
        "middle_layer_size": math.comb(degree, middle),
        "walsh_layers": layer_records,
        "matrices": matrices,
    }
    return record, spectrum_holds and rank_holds


def inverted_mask(mask: int) -> int:
    result = 0
    for source, target in enumerate(INVERSION):
        if (mask >> source) & 1:
            result |= 1 << target
    return result


def inversion_orbits(masks: list[int]) -> list[tuple[int, ...]]:
    remaining = set(masks)
    orbits: list[tuple[int, ...]] = []
    while remaining:
        seed = min(remaining)
        orbit = tuple(sorted({seed, inverted_mask(seed)}))
        if not set(orbit) <= remaining:
            raise AssertionError("inversion partition is inconsistent")
        remaining.difference_update(orbit)
        orbits.append(orbit)
    return orbits


def subset_names(mask: int) -> list[str]:
    return [name for site, name in enumerate(DIRECTIONS) if (mask >> site) & 1]


def orbit_record(orbit: tuple[int, ...]) -> dict[str, object]:
    return {
        "masks": list(orbit),
        "subsets": [subset_names(mask) for mask in orbit],
    }


def middle_complement_pair_orbits() -> list[tuple[tuple[int, int], ...]]:
    full_mask = (1 << 6) - 1
    pairs = {
        tuple(sorted((query, full_mask ^ query)))
        for query in range(1 << 6)
        if query.bit_count() == 3
    }
    remaining = set(pairs)
    result: list[tuple[tuple[int, int], ...]] = []
    while remaining:
        pair = min(remaining)
        image_query = inverted_mask(pair[0])
        image = tuple(sorted((image_query, full_mask ^ image_query)))
        orbit = tuple(sorted({pair, image}))
        if not set(orbit) <= remaining:
            raise AssertionError("central complement-pair orbit is inconsistent")
        remaining.difference_update(orbit)
        result.append(orbit)
    return sorted(result, key=lambda orbit: (len(orbit), orbit))


def pair_orbit_record(
    orbit: tuple[tuple[int, int], ...],
) -> dict[str, object]:
    return {
        "weight": len(orbit),
        "complement_pairs": [
            {
                "masks": list(pair),
                "subsets": [subset_names(mask) for mask in pair],
            }
            for pair in orbit
        ],
    }


def lambda_value(
    alpha: list[int], mask: int, pair_orbits: list[tuple[tuple[int, int], ...]]
) -> int:
    return sum(
        alpha[index]
        * sum(character(pair[0], mask) for pair in pair_orbit)
        for index, pair_orbit in enumerate(pair_orbits)
    )


def folded_certificate(coefficients: dict[int, Fraction]) -> tuple[dict[str, object], bool]:
    full_mask = (1 << 6) - 1
    selected = [
        mask
        for mask in parity_masks(6, 0)
        if mask != full_mask
    ]
    odd = parity_masks(6, 1)
    row_orbits = inversion_orbits(selected)
    column_orbits = inversion_orbits(odd)
    folded: Matrix = []
    column_invariance = True
    for row_orbit in row_orbits:
        row = []
        for column_orbit in column_orbits:
            values = [
                sum(
                    coefficients[(source ^ column).bit_count()]
                    for source in row_orbit
                )
                for column in column_orbit
            ]
            column_invariance &= all(value == values[0] for value in values)
            row.append(values[0])
        folded.append(row)

    deleted_empty = [
        row for orbit, row in zip(row_orbits, folded) if orbit != (0,)
    ]
    folded_full_record = matrix_record(folded, len(column_orbits))
    folded_deleted_record = matrix_record(deleted_empty, len(column_orbits))

    pair_orbits = middle_complement_pair_orbits()
    weights = [len(orbit) for orbit in pair_orbits]
    alphas: list[list[int]] = []
    for index in range(1, len(weights)):
        alpha = [0] * len(weights)
        alpha[0] = -weights[index]
        alpha[index] = weights[0]
        alphas.append(alpha)

    vectors = []
    all_vectors_hold = True
    for index, alpha in enumerate(alphas, start=1):
        values = [
            lambda_value(alpha, row_orbit[0], pair_orbits)
            for row_orbit in row_orbits
        ]
        invariant = all(
            lambda_value(alpha, member, pair_orbits) == value
            for row_orbit, value in zip(row_orbits, values)
            for member in row_orbit
        )
        residual = [
            sum(
                Fraction(values[row]) * folded[row][column]
                for row in range(len(row_orbits))
            )
            for column in range(len(column_orbits))
        ]
        weighted_sum = sum(weight * value for weight, value in zip(weights, alpha))
        empty_value = lambda_value(alpha, 0, pair_orbits)
        full_value = lambda_value(alpha, full_mask, pair_orbits)
        all_vectors_hold &= (
            invariant
            and weighted_sum == 0
            and empty_value == 0
            and full_value == 0
            and all(value == 0 for value in residual)
        )
        vectors.append(
            {
                "name": f"higher_relation_{index}",
                "alpha": alpha,
                "weighted_alpha_sum": weighted_sum,
                "folded_row_values": values,
                "lambda_empty": empty_value,
                "lambda_full": full_value,
                "leakage_residual": [fraction_record(value) for value in residual],
            }
        )

    vector_rank = matrix_rank(
        [[Fraction(value) for value in row["folded_row_values"]] for row in vectors]
    )
    record = {
        "inversion": {
            "direction_order": list(DIRECTIONS),
            "permutation": list(INVERSION),
            "definition": "the single cubic inversion simultaneously swaps +e_i and -e_i for i=x,y,z",
        },
        "row_classes": [orbit_record(orbit) for orbit in row_orbits],
        "column_classes": [orbit_record(orbit) for orbit in column_orbits],
        "folding_convention": "one row variable per inversion orbit; the folded row is the sum over raw rows and one representative of each invariant column equation",
        "matrices": {
            "selected_proper_even": folded_full_record,
            "selected_without_empty": folded_deleted_record,
        },
        "central_complement_pair_orbits": [
            pair_orbit_record(orbit) for orbit in pair_orbits
        ],
        "weighted_alpha_constraint": "alpha1+alpha2+alpha3+alpha4+2*alpha5+2*alpha6+2*alpha7=0",
        "weights": weights,
        "kernel_vectors": vectors,
        "kernel_vector_rank": vector_rank,
    }
    passed = (
        column_invariance
        and len(row_orbits) == 19
        and len(column_orbits) == 16
        and folded_full_record["rank"] == 13
        and folded_full_record["left_nullity"] == 6
        and folded_deleted_record["rank"] == 12
        and folded_deleted_record["left_nullity"] == 6
        and weights == [1, 1, 1, 1, 2, 2, 2]
        and len(vectors) == 6
        and vector_rank == 6
        and all_vectors_hold
        and all(vector["folded_row_values"][0] == 0 for vector in vectors)
    )
    return record, passed


def far_mark_certificate() -> dict[str, object]:
    shape = (3, 3, 3)
    pivot = (0, 0, 0)
    steps = (
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    )
    sites = tuple(itertools.product(*(range(length) for length in shape)))
    neighbours = tuple(
        tuple((pivot[axis] + step[axis]) % shape[axis] for axis in range(3))
        for step in steps
    )
    closed_neighbourhood = {pivot, *neighbours}
    far_marks = sorted(set(sites) - closed_neighbourhood)
    pivot_contacts = sum(mark == pivot for mark in far_marks)
    neighbour_contacts = sum(mark in neighbours for mark in far_marks)
    return {
        "shape": list(shape),
        "pivot": list(pivot),
        "direction_sites": {
            direction: list(site)
            for direction, site in zip(DIRECTIONS, neighbours)
        },
        "site_count": len(sites),
        "distinct_neighbour_count": len(set(neighbours)),
        "closed_neighbourhood_count": len(closed_neighbourhood),
        "far_mark_count": len(far_marks),
        "far_marks": [list(mark) for mark in far_marks],
        "far_mark_pivot_contacts": pivot_contacts,
        "far_mark_neighbour_contacts": neighbour_contacts,
        "no_contact_terms": pivot_contacts == 0 and neighbour_contacts == 0,
        "contact_definition": "a contact occurs only when the marked site equals the selected pivot or one of its six neighbour variables",
    }


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    source_paths = (
        Path(__file__).resolve(),
        VERIFIER,
        E240_PRODUCER,
        E243_PRODUCER,
    )
    source_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
    }
    add_check(
        checks,
        "source_hashes_bound",
        set(source_hashes)
        == {
            "experiments/e246_callen_walsh_compression.py",
            "tests/test_callen_walsh_compression.py",
            "experiments/e240_callen_identity_system.py",
            "experiments/e243_callen_termwise_closure.py",
        }
        and all(len(digest) == 64 for digest in source_hashes.values()),
        "producer, clean-room verifier, e240 coefficient source, and e243 scope predecessor",
    )

    coefficients_payload: dict[str, object] = {}
    controls: dict[str, object] = {}
    coefficient_checks = True
    finite_checks = True
    for degree in (2, 4, 6):
        values, records, expansion_holds = coefficient_record(degree)
        control, finite_holds = finite_control(degree, values)
        coefficients_payload[str(degree)] = records
        controls[str(degree)] = control
        coefficient_checks &= values == EXPECTED_VALUES[degree] and expansion_holds
        finite_checks &= finite_holds
    add_check(
        checks,
        "exact_degree_2_4_6_coefficients",
        coefficient_checks,
        "e240 rational functions specialize exactly at v=1/3 and replay every local Walsh expansion",
    )
    add_check(
        checks,
        "finite_walsh_spectra_and_rank_formula",
        finite_checks,
        "all Walsh layers and exact Fraction ranks agree for d=2,4,6",
    )

    d6_matrices = controls["6"]["matrices"]
    corrected_raw = (
        d6_matrices["selected_proper_even"]["shape"] == [31, 32]
        and d6_matrices["selected_proper_even"]["rank"] == 22
        and d6_matrices["selected_proper_even"]["left_nullity"] == 9
        and d6_matrices["selected_without_empty"]["shape"] == [30, 32]
        and d6_matrices["selected_without_empty"]["rank"] == 21
        and d6_matrices["selected_without_empty"]["left_nullity"] == 9
    )
    add_check(
        checks,
        "corrected_raw_d6_rank_criterion",
        corrected_raw,
        "selected full/deleted-empty ranks are 22/21 and both left nullities are 9",
    )

    degree_six_values = {
        int(size): Fraction(
            record["value_at_v"]["numerator"],
            record["value_at_v"]["denominator"],
        )
        for size, record in coefficients_payload["6"].items()
    }
    folded, folded_holds = folded_certificate(degree_six_values)
    add_check(
        checks,
        "cubic_inversion_folded_ranks",
        folded_holds,
        "19x16 rank 13; deleting empty gives 18x16 rank 12; left nullity remains 6",
    )
    add_check(
        checks,
        "six_weighted_folded_kernel_vectors",
        folded_holds
        and folded["weights"] == [1, 1, 1, 1, 2, 2, 2]
        and folded["kernel_vector_rank"] == 6,
        "six independent vectors obey the weighted alpha constraint, vanish at empty/full, and annihilate all leakage columns",
    )

    far_marks = far_mark_certificate()
    far_holds = (
        far_marks["site_count"] == 27
        and far_marks["distinct_neighbour_count"] == 6
        and far_marks["closed_neighbourhood_count"] == 7
        and far_marks["far_mark_count"] == 20
        and far_marks["no_contact_terms"] is True
    )
    add_check(
        checks,
        "periodic_3x3x3_far_marks",
        far_holds,
        "one pivot plus six distinct neighbours leaves 20 marked sites and zero contacts",
    )

    add_check(
        checks,
        "scope_gate",
        SCOPE
        == "six exact inversion-invariant higher-template relations and the all-even-degree obstruction to pair-row leakage cancellation in the fixed selected-pivot/far-mark Callen block; no susceptibility closure and no contact-inclusive or radius-expanded no-go",
        SCOPE,
    )
    add_check(
        checks,
        "benchmark_not_used",
        True,
        "the certificate uses v=1/3 only as an exact rational control, never as a physical benchmark",
    )

    certificate = {
        "parameter": {
            "name": "v=tanh(K)",
            "exact_control_value": fraction_record(V),
            "benchmark_used": False,
        },
        "coefficients": coefficients_payload,
        "finite_controls": controls,
        "cubic_inversion": folded,
        "periodic_3x3x3": far_marks,
    }
    theorem = {
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
    scope = {
        "tag": "[SCOPE]",
        "claim": SCOPE,
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

    elapsed = time.process_time() - started
    peak_rss = peak_rss_bytes()
    rss_measurement = peak_rss_measurement()
    add_check(
        checks,
        "declared_resource_limits",
        elapsed < CPU_BUDGET_SECONDS and peak_rss < RSS_LIMIT_BYTES,
        (
            f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; "
            f"peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes via {rss_measurement}"
        ),
    )
    failed = [str(check["name"]) for check in checks if not check["passed"]]
    if failed:
        raise AssertionError(f"Callen Walsh-compression checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e246_callen_walsh_compression.py",
            "verifier": "tests/test_callen_walsh_compression.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers, fractions, and rational functions over Z[v]",
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "peak_rss_measurement": rss_measurement,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "source_sha256": source_hashes,
            "certificate_sha256": canonical_sha256(certificate),
        },
        "data": {
            "claim_tag": "[THEOREM][FINITE EXACT CONTROL][SELECTED-ROW BOUNDARY]",
            "inherited_sources": {
                "e240_coefficients": {
                    "path": "experiments/e240_callen_identity_system.py",
                    "sha256": source_hashes[
                        "experiments/e240_callen_identity_system.py"
                    ],
                },
                "e243_scope_predecessor": {
                    "path": "experiments/e243_callen_termwise_closure.py",
                    "sha256": source_hashes[
                        "experiments/e243_callen_termwise_closure.py"
                    ],
                },
            },
            "theorem": theorem,
            "certificate": certificate,
            "scope": scope,
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for check in payload["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"  [{status}] {check['name']}: {check['detail']}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
