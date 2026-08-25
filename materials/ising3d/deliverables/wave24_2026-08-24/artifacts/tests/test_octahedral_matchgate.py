#!/usr/bin/env python3
"""Clean-room exact verifier for the octahedral-star matchgate theorem."""

from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "octahedral_matchgate.json"
PRODUCER = ROOT / "experiments" / "e244_octahedral_matchgate.py"
E234_PRODUCER = ROOT / "experiments" / "e234_star_decimation.py"
DLMF_SOURCE = ROOT / "sources" / "dlmf_4_36.html"
DIMENSION = 6
AUDIT_POINTS = (Fraction(1, 3), Fraction(1, 2))
RSS_LIMIT_BYTES = 2 * 1024**3
FAILURES: list[str] = []

CONTROL_GRAPHS: dict[str, tuple[int, int, tuple[tuple[int, int], ...]]] = {
    "C4": (2, 2, ((0, 0), (1, 0), (0, 1), (1, 1))),
    "K2_3": (
        2,
        3,
        tuple((left, right) for right in range(3) for left in range(2)),
    ),
}


def check(name: str, passed: bool, detail: str = "") -> None:
    print(
        f"  [{'PASS' if passed else 'FAIL'}] {name}"
        + (f": {detail}" if detail else "")
    )
    if not passed:
        FAILURES.append(name)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def expected_signature(mask: int, v: Fraction) -> Fraction:
    if mask.bit_count() % 2:
        return Fraction(0)
    return v ** mask.bit_count()


def center_sum_value(spin_mask: int, v: Fraction) -> Fraction:
    total = Fraction(0)
    for center in (-1, 1):
        product = Fraction(1)
        for site in range(DIMENSION):
            boundary_spin = -1 if (spin_mask >> site) & 1 else 1
            product *= 1 + v * center * boundary_spin
        total += product
    return total / 2


def independent_walsh_signature(v: Fraction) -> tuple[Fraction, ...]:
    values = [center_sum_value(state, v) for state in range(1 << DIMENSION)]
    coefficients: list[Fraction] = []
    for support in range(1 << DIMENSION):
        numerator = Fraction(0)
        for state, value in enumerate(values):
            sign = -1 if (state & support).bit_count() % 2 else 1
            numerator += sign * value
        coefficients.append(numerator / (1 << DIMENSION))
    return tuple(coefficients)


def combinatorial_log_counts(degree: int, subset_size: int) -> dict[int, int]:
    """Orbit-count reconstruction, independent of the producer's spin loop."""

    counts: Counter[int] = Counter()
    for inside_minus in range(subset_size + 1):
        sign = -1 if inside_minus % 2 else 1
        inside_count = math.comb(subset_size, inside_minus)
        for outside_minus in range(degree - subset_size + 1):
            multiplicity = inside_count * math.comb(
                degree - subset_size, outside_minus
            )
            magnetization = abs(
                degree - 2 * (inside_minus + outside_minus)
            )
            counts[magnetization] += sign * multiplicity
    return dict(sorted(counts.items()))


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            out[left_degree + right_degree] += left_value * right_value
    return out


def polynomial_power(base: list[int], exponent: int) -> list[int]:
    out = [1]
    for _ in range(exponent):
        out = polynomial_multiply(out, base)
    return out


def polynomial_add(left: list[int], right: list[int]) -> list[int]:
    width = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(width)
    ]


def verify_log_sign_controls(stored: dict[str, object]) -> bool:
    expected_top = {
        2: {0: -2, 2: 2},
        4: {0: 6, 2: -8, 4: 2},
        6: {0: -20, 2: 30, 4: -12, 6: 2},
    }
    passed = set(stored["degrees"]) == {"2", "4", "6"}
    for degree in (2, 4, 6):
        degree_row = stored["degrees"][str(degree)]
        passed &= degree_row["spin_states"] == 1 << degree
        sectors = degree_row["nonempty_even_sectors"]
        passed &= set(sectors) == {
            str(size) for size in range(2, degree + 1, 2)
        }
        for subset_size in range(2, degree + 1, 2):
            counts = combinatorial_log_counts(degree, subset_size)
            row = sectors[str(subset_size)]
            stored_counts = {
                int(magnetization): int(coefficient)
                for magnetization, coefficient in row[
                    "scaled_integer_multiplicities_by_abs_magnetization"
                ].items()
            }
            expected_sign = (
                "positive"
                if (subset_size // 2 + 1) % 2 == 0
                else "negative"
            )
            passed &= (
                stored_counts == counts
                and row["Walsh_normalization_denominator"] == 1 << degree
                and row["theorem_sign"] == expected_sign
            )
            if subset_size == degree:
                passed &= counts == expected_top[degree]

    certificates = stored["top_sector_ratio_certificates"]
    d4_gap = certificates["d4"][
        "denominator_minus_numerator_coefficients_ascending"
    ]
    passed &= d4_gap == polynomial_power([-1, 1], 2) == [1, -2, 1]

    left = polynomial_multiply([0] * 8 + [1], [-3, 4])
    right = polynomial_power([-1, 2], 6)
    width = max(len(left), len(right))
    difference = [
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(width)
    ]
    shifted_q = certificates["d6"][
        "Q6_at_1_plus_y_coefficients_ascending"
    ]
    reconstructed = [0]
    for degree, coefficient in enumerate(shifted_q):
        term = [
            coefficient * value
            for value in polynomial_power([-1, 1], degree + 3)
        ]
        reconstructed = polynomial_add(reconstructed, term)
    passed &= (
        shifted_q == [8, 54, 144, 188, 120, 33, 4]
        and all(coefficient > 0 for coefficient in shifted_q)
        and reconstructed == difference
    )
    return bool(passed)


def ordered_pairings(
    vertices: tuple[int, ...],
) -> Iterator[tuple[tuple[int, int], ...]]:
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for partner_index in range(1, len(vertices)):
        partner = vertices[partner_index]
        remainder = vertices[1:partner_index] + vertices[partner_index + 1 :]
        for rest in ordered_pairings(remainder):
            yield ((first, partner), *rest)


def crossing_sign(pairs: tuple[tuple[int, int], ...]) -> int:
    crossings = 0
    for first_index, (left, right) in enumerate(pairs):
        for other_left, other_right in pairs[first_index + 1 :]:
            crossings += int(
                left < other_left < right < other_right
                or other_left < left < other_right < right
            )
    return -1 if crossings % 2 else 1


def pairing_pfaffian(matrix: list[list[Fraction]]) -> Fraction:
    if len(matrix) % 2:
        return Fraction(0)
    total = Fraction(0)
    for pairs in ordered_pairings(tuple(range(len(matrix)))):
        term = Fraction(crossing_sign(pairs))
        for left, right in pairs:
            term *= matrix[left][right]
        total += term
    return total


def principal_matrix(
    matrix: list[list[Fraction]], mask: int
) -> list[list[Fraction]]:
    selected = [index for index in range(len(matrix)) if (mask >> index) & 1]
    return [[matrix[row][column] for column in selected] for row in selected]


def make_generator(v: Fraction) -> list[list[Fraction]]:
    return [
        [
            Fraction(0)
            if row == column
            else (v * v if row < column else -(v * v))
            for column in range(DIMENSION)
        ]
        for row in range(DIMENSION)
    ]


def matchgate_inventory(
    values: tuple[Fraction, ...]
) -> tuple[int, dict[str, int], bool]:
    total = 0
    active: Counter[int] = Counter()
    all_zero = True
    for alpha, beta in itertools.product(range(1 << DIMENSION), repeat=2):
        positions = tuple(
            site
            for site in range(DIMENSION)
            if ((alpha ^ beta) >> site) & 1
        )
        residual = Fraction(0)
        nonzero_term = False
        for one_based, site in enumerate(positions, start=1):
            term = (
                (-1 if one_based % 2 else 1)
                * values[alpha ^ (1 << site)]
                * values[beta ^ (1 << site)]
            )
            residual += term
            nonzero_term |= term != 0
        total += 1
        if nonzero_term:
            active[len(positions)] += 1
        all_zero &= residual == 0
    return (
        total,
        {str(key): value for key, value in sorted(active.items())},
        all_zero,
    )


def embed(local_state: int, sites: tuple[int, ...]) -> int:
    mask = 0
    for local_index, site in enumerate(sites):
        if (local_state >> local_index) & 1:
            mask |= 1 << site
    return mask


def flattening(left_mask: int, v: Fraction) -> list[list[Fraction]]:
    left = tuple(
        site for site in range(DIMENSION) if (left_mask >> site) & 1
    )
    right = tuple(
        site for site in range(DIMENSION) if not ((left_mask >> site) & 1)
    )
    return [
        [
            expected_signature(embed(row, left) | embed(column, right), v)
            for column in range(1 << len(right))
        ]
        for row in range(1 << len(left))
    ]


def row_content(row: list[int]) -> int:
    content = 0
    for value in row:
        content = math.gcd(content, abs(value))
    return content


def integer_fraction_free_rank(matrix: list[list[Fraction]]) -> int:
    """Clear row denominators, then perform exact integer cross-elimination."""

    if not matrix:
        return 0
    width = len(matrix[0])
    work: list[list[int]] = []
    for source_row in matrix:
        denominator = 1
        for value in source_row:
            denominator = math.lcm(
                denominator, Fraction(value).denominator
            )
        row = [int(Fraction(value) * denominator) for value in source_row]
        content = row_content(row)
        work.append([value // content for value in row] if content else row)

    rank = 0
    for column in range(width):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row_index in range(rank + 1, len(work)):
            target_value = work[row_index][column]
            if target_value == 0:
                continue
            common = math.gcd(abs(pivot_value), abs(target_value))
            pivot_multiplier = pivot_value // common
            target_multiplier = target_value // common
            work[row_index] = [
                pivot_multiplier * work[row_index][index]
                - target_multiplier * work[rank][index]
                for index in range(width)
            ]
            content = row_content(work[row_index])
            if content > 1:
                work[row_index] = [
                    value // content for value in work[row_index]
                ]
        rank += 1
        if rank == len(work):
            break
    return rank


def parity_gram(
    v: Fraction,
) -> tuple[list[list[Fraction]], Fraction, Fraction]:
    even = [
        v ** state.bit_count()
        if state.bit_count() % 2 == 0
        else Fraction(0)
        for state in range(8)
    ]
    odd = [
        v ** state.bit_count()
        if state.bit_count() % 2
        else Fraction(0)
        for state in range(8)
    ]
    matrix = [
        [
            even[row] * even[column] + odd[row] * odd[column]
            for column in range(8)
        ]
        for row in range(8)
    ]
    return (
        matrix,
        sum((entry * entry for entry in even), Fraction(0)),
        sum((entry * entry for entry in odd), Fraction(0)),
    )


def ternary_rectangle_inventory() -> dict[str, object]:
    free_histogram: Counter[int] = Counter()
    even_singletons = 0
    odd_singletons = 0
    mixed_non_singletons = 0
    even_only = 0
    largest_even_only = 0
    for code in range(3**DIMENSION):
        residual = code
        choices: list[tuple[int, ...]] = []
        free = 0
        for _site in range(DIMENSION):
            digit = residual % 3
            residual //= 3
            if digit == 0:
                choices.append((0,))
            elif digit == 1:
                choices.append((1,))
            else:
                choices.append((0, 1))
                free += 1
        states = [
            sum(bit << site for site, bit in enumerate(bits))
            for bits in itertools.product(*choices)
        ]
        parities = {state.bit_count() % 2 for state in states}
        free_histogram[free] += 1
        if len(states) == 1:
            if parities == {0}:
                even_singletons += 1
            else:
                odd_singletons += 1
        elif parities == {0, 1}:
            mixed_non_singletons += 1
        if parities == {0}:
            even_only += 1
            largest_even_only = max(largest_even_only, len(states))
    return {
        "dimension": DIMENSION,
        "nonempty_support_choices_per_leg": 3,
        "rectangles_checked": 3**DIMENSION,
        "by_free_coordinates": {
            str(key): value
            for key, value in sorted(free_histogram.items())
        },
        "singleton_even": even_singletons,
        "singleton_odd": odd_singletons,
        "non_singleton": 3**DIMENSION - even_singletons - odd_singletons,
        "non_singleton_containing_both_parities": mixed_non_singletons,
        "rectangles_contained_in_even_support": even_only,
        "maximum_even_only_rectangle_size": largest_even_only,
    }


def boundary_dp_polynomial(
    left_order: int,
    right_order: int,
    edges: tuple[tuple[int, int], ...],
) -> dict[int, int]:
    states: dict[int, Counter[int]] = {0: Counter({0: 1})}
    for left, right in edges:
        boundary_bit = (1 << left) | (1 << (left_order + right))
        enlarged: dict[int, Counter[int]] = {
            boundary: Counter(counts) for boundary, counts in states.items()
        }
        for boundary, counts in states.items():
            target = enlarged.setdefault(boundary ^ boundary_bit, Counter())
            for degree, count in counts.items():
                target[degree + 1] += count
        states = enlarged
    return dict(sorted(states.get(0, Counter()).items()))


def spin_polynomial_value(
    left_order: int,
    right_order: int,
    edges: tuple[tuple[int, int], ...],
    v: Fraction,
) -> Fraction:
    order = left_order + right_order
    total = Fraction(0)
    for state in range(1 << order):
        term = Fraction(1)
        for left, right in edges:
            left_spin = -1 if (state >> left) & 1 else 1
            right_spin = (
                -1 if (state >> (left_order + right)) & 1 else 1
            )
            term *= 1 + v * left_spin * right_spin
        total += term
    return total / (1 << order)


def spin_basis_matrix(v: Fraction) -> tuple[list[list[Fraction]], bool]:
    values: list[Fraction] = []
    positive = True
    for state in range(1 << DIMENSION):
        atoms: list[Fraction] = []
        for center in (-1, 1):
            atom = Fraction(1)
            for site in range(DIMENSION):
                spin = -1 if (state >> site) & 1 else 1
                factor = 1 + v * center * spin
                positive &= factor > 0
                atom *= factor
            atoms.append(atom)
        values.append(sum(atoms, Fraction(0)))
    matrix = [
        [
            values[row | (column << 1)]
            for column in range(1 << (DIMENSION - 1))
        ]
        for row in range(2)
    ]
    return matrix, positive


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check(
        "artifact top-level schema",
        set(artifact) == {"meta", "data", "checks"},
    )
    meta = artifact["meta"]
    data = artifact["data"]
    hashes = meta["source_sha256"]
    check(
        "producer/verifier/artifact paths",
        meta["producer"] == "experiments/e244_octahedral_matchgate.py"
        and meta["verifier"] == "tests/test_octahedral_matchgate.py"
        and meta["artifact"]
        == "results/integrability/octahedral_matchgate.json",
    )
    check(
        "current producer hash",
        hashes["experiments/e244_octahedral_matchgate.py"]
        == file_sha256(PRODUCER),
    )
    check(
        "current verifier hash",
        hashes["tests/test_octahedral_matchgate.py"]
        == file_sha256(Path(__file__).resolve()),
    )
    e234_hash = file_sha256(E234_PRODUCER)
    check(
        "current e234 producer hash",
        hashes["experiments/e234_star_decimation.py"] == e234_hash
        and data["inherited_source"]["sha256"] == e234_hash
        and data["inherited_source"]["producer"]
        == "experiments/e234_star_decimation.py",
    )
    check(
        "archived DLMF product source",
        hashes["sources/dlmf_4_36.html"] == file_sha256(DLMF_SOURCE)
        and "4.36.2" in DLMF_SOURCE.read_text()
        and "cosh" in DLMF_SOURCE.read_text(),
    )

    expected_j = [
        [
            0 if row == column else (1 if row < column else -1)
            for column in range(DIMENSION)
        ]
        for row in range(DIMENSION)
    ]
    check(
        "explicit J6 generator",
        data["pfaffian_matchgate"]["J_6"] == expected_j
        and data["pfaffian_matchgate"]["even_principal_subset_count"]
        == 32,
    )

    log_theorem = data["all_even_star_log_walsh_sign"]
    check(
        "independent all-even-star log-Walsh controls",
        log_theorem["claim_tag"] == "[THEOREM]"
        and "sign(c_A)=(-1)^(|A|/2+1)" in log_theorem["statement"]
        and log_theorem["analytic_certificate"]["product_source"]
        == "NIST DLMF 4.36.E2, https://dlmf.nist.gov/4.36.E2"
        and "cos(t)^(d-2r)" in log_theorem["analytic_certificate"][
            "Gaussian_Walsh_transform"
        ]
        and verify_log_sign_controls(log_theorem["finite_exact_controls"]),
    )

    all_walsh = True
    all_pfaffians = True
    all_identities = True
    all_ranks = True
    all_psd = True
    all_nonnegative_upper = True
    all_spin_positive = True
    for v in AUDIT_POINTS:
        key = str(v)
        values = independent_walsh_signature(v)
        all_walsh &= values == tuple(
            expected_signature(mask, v) for mask in range(1 << DIMENSION)
        )
        stored_walsh = data["star_walsh"]["audits"][key]
        all_walsh &= (
            stored_walsh["spin_states"] == 64
            and stored_walsh["walsh_entries"] == 64
            and stored_walsh["walsh_summands"] == 4096
            and stored_walsh["all_entries_match_g"] is True
        )
        for size in range(DIMENSION + 1):
            stored = stored_walsh["coefficients_by_subset_size"][str(size)]
            expected_value = v**size if size % 2 == 0 else Fraction(0)
            all_walsh &= (
                stored["entry_count"] == math.comb(DIMENSION, size)
                and Fraction(stored["value"]) == expected_value
                and stored["constant_on_weight_class"] is True
            )

        generator = make_generator(v)
        by_size: Counter[int] = Counter()
        for mask in range(1 << DIMENSION):
            if mask.bit_count() % 2:
                continue
            observed = pairing_pfaffian(principal_matrix(generator, mask))
            all_pfaffians &= observed == values[mask]
            by_size[mask.bit_count()] += 1
        stored_pfaffian = data["pfaffian_matchgate"][
            "principal_pfaffian_audits"
        ][key]
        all_pfaffians &= (
            Fraction(stored_pfaffian["matrix_scale_v_squared"]) == v * v
            and stored_pfaffian["even_principal_pfaffians_checked"] == 32
            and stored_pfaffian["all_equal_signature"] is True
        )
        for size, count in sorted(by_size.items()):
            stored = stored_pfaffian["by_subset_size"][str(size)]
            all_pfaffians &= (
                stored["count"] == count
                and Fraction(stored["value"]) == v**size
                and stored["constant"] is True
            )

        pair_count, active, identities_zero = matchgate_inventory(values)
        stored_identities = data["pfaffian_matchgate"][
            "matchgate_identity_audits"
        ][key]
        all_identities &= (
            pair_count == 4096
            and active == {"2": 480, "4": 480, "6": 32}
            and identities_zero
            and stored_identities["ordered_bitstring_pairs"] == pair_count
            and stored_identities["algebraically_active_by_term_count"]
            == active
            and stored_identities["algebraically_active_total"] == 992
            and stored_identities["all_residuals_zero"] is True
        )

        ranks_by_size: dict[str, dict[str, object]] = {}
        all_flattening_ranks: list[int] = []
        for left_size in range(1, DIMENSION):
            ranks = [
                integer_fraction_free_rank(flattening(left_mask, v))
                for left_mask in range(1, (1 << DIMENSION) - 1)
                if left_mask.bit_count() == left_size
            ]
            all_flattening_ranks.extend(ranks)
            ranks_by_size[str(left_size)] = {
                "ordered_left_subsets": len(ranks),
                "rank_histogram": {"2": len(ranks)},
            }
        signed_formula = all(
            values[mask]
            == (v ** mask.bit_count() + (-v) ** mask.bit_count()) / 2
            for mask in range(1 << DIMENSION)
        )
        stored_cp = data["ordinary_real_CP_rank"]["audits"][key]
        all_ranks &= (
            signed_formula
            and len(all_flattening_ranks) == 62
            and all(rank == 2 for rank in all_flattening_ranks)
            and stored_cp["signature_entries_checked"] == 64
            and stored_cp["signed_two_atom_decomposition_exact"] is True
            and stored_cp["ordered_nontrivial_flattenings"] == 62
            and stored_cp["all_flattening_ranks_two"] is True
            and stored_cp["rank_histogram"] == {"2": 62}
            and stored_cp["by_left_size"] == ranks_by_size
        )

        gram, even_eigenvalue, odd_eigenvalue = parity_gram(v)
        actual_three = flattening(0b000111, v)
        stored_three = stored_cp["three_by_three_operator"]
        all_psd &= (
            actual_three == gram
            and integer_fraction_free_rank(actual_three) == 2
            and even_eigenvalue == 1 + 3 * v**4 > 0
            and odd_eigenvalue == 3 * v**2 + v**6 > 0
            and stored_three["rank"] == 2
            and stored_three["symmetric"] is True
            and stored_three["equals_two_parity_outer_products"] is True
            and Fraction(
                stored_three["nonzero_eigenvalues"]["even_1_plus_3v4"]
            )
            == even_eigenvalue
            and Fraction(
                stored_three["nonzero_eigenvalues"]["odd_3v2_plus_v6"]
            )
            == odd_eigenvalue
            and stored_three["zero_eigenvalue_multiplicity"] == 6
            and stored_three["positive_semidefinite_certificate"] is True
        )

        even_states = {
            state
            for state in range(1 << DIMENSION)
            if state.bit_count() % 2 == 0
        }
        singleton_reconstruction = tuple(
            values[state] if state in even_states else Fraction(0)
            for state in range(1 << DIMENSION)
        )
        stored_nonnegative = data["nonnegative_all_leg_CP_rank"][
            "audits"
        ][key]
        all_nonnegative_upper &= (
            singleton_reconstruction == values
            and all(values[state] > 0 for state in even_states)
            and stored_nonnegative["singleton_atoms"] == 32
            and stored_nonnegative["all_atom_weights_positive"] is True
            and stored_nonnegative["reconstructs_signature"] is True
        )

        spin_matrix, positive = spin_basis_matrix(v)
        stored_spin = data["spin_basis_positive_rank"]["audits"][key]
        all_spin_positive &= (
            positive
            and integer_fraction_free_rank(spin_matrix) == 2
            and stored_spin["spin_states_checked"] == 64
            and stored_spin["two_positive_atoms_exact"] is True
            and stored_spin["one_by_five_flattening_rank"] == 2
        )

    check("independent star Walsh transforms", all_walsh)
    check("pairing-enumeration principal Pfaffians", all_pfaffians)
    check("exhaustive rational matchgate identities", all_identities)
    check("fraction-free ranks for all flattenings", all_ranks)
    check("exact 3|3 PSD parity decomposition", all_psd)
    check("nonnegative singleton upper decomposition", all_nonnegative_upper)
    check("spin-basis positive rank two", all_spin_positive)

    cycle_data = data["cycle_space_factorization"]
    stored_controls = {
        row["name"]: row for row in cycle_data["controls"]
    }
    controls_pass = (
        set(stored_controls) == set(CONTROL_GRAPHS)
        and cycle_data["retained_spin_average_normalization"] == "2^(-|V|)"
        and cycle_data["normalized_identity"]
        == "2^(-|V|)*sum_(sigma_V) product_(u in U) g_u(sigma_delta(u)) = sum_(F subset E: boundary(F)=empty) v^|F|"
        and cycle_data["raw_sum_identity"]
        == "sum_(sigma_V) product_(u in U) g_u(sigma_delta(u)) = 2^|V|*sum_(F subset E: boundary(F)=empty) v^|F|"
    )
    expected_polynomials = {
        "C4": {0: 1, 4: 1},
        "K2_3": {0: 1, 4: 3},
    }
    for name, (left_order, right_order, edges) in CONTROL_GRAPHS.items():
        polynomial = boundary_dp_polynomial(left_order, right_order, edges)
        row = stored_controls[name]
        controls_pass &= (
            polynomial == expected_polynomials[name]
            and row["left_vertices"] == left_order
            and row["right_vertices"] == right_order
            and {tuple(edge) for edge in row["edges"]} == set(edges)
            and row["edge_count"] == len(edges)
            and row["retained_spin_sum_factor"] == 1 << right_order
            and row["coefficient_by_edge_count"]
            == {str(degree): count for degree, count in polynomial.items()}
        )
        for v in AUDIT_POINTS:
            polynomial_value = sum(
                count * v**degree for degree, count in polynomial.items()
            )
            controls_pass &= (
                spin_polynomial_value(left_order, right_order, edges, v)
                == polynomial_value
                and Fraction(row["values"][str(v)]) == polynomial_value
            )
    check("independent cycle-space dynamic programs", controls_pass)

    rectangles = ternary_rectangle_inventory()
    stored_rectangles = data["nonnegative_all_leg_CP_rank"][
        "rectangle_enumeration"
    ]
    expected_even_states = [
        format(state, "06b")
        for state in range(1 << DIMENSION)
        if state.bit_count() % 2 == 0
    ]
    check(
        "independent 3^6 rectangle barrier",
        rectangles == stored_rectangles
        and rectangles["rectangles_checked"] == 729
        and rectangles["non_singleton_containing_both_parities"] == 665
        and rectangles["rectangles_contained_in_even_support"] == 32
        and rectangles["maximum_even_only_rectangle_size"] == 1
        and data["nonnegative_all_leg_CP_rank"]["rank"] == 32
        and data["nonnegative_all_leg_CP_rank"]["even_singleton_states"]
        == expected_even_states,
    )

    bounds = data["finite_enumeration_bounds"]
    check(
        "finite enumeration bounds",
        bounds
        == {
            "Walsh_summands_per_audit_point": 4096,
            "principal_Pfaffians_per_audit_point": 32,
            "largest_Pfaffian_pairing_expansion": 15,
            "ordered_matchgate_identities_per_audit_point": 4096,
            "ordered_nontrivial_flattenings_per_audit_point": 62,
            "largest_flattening_core": "8x8",
            "Cartesian_rectangles": 729,
            "maximum_points_in_one_rectangle": 64,
            "largest_cycle_control_edge_subset_space": 64,
            "largest_cycle_control_local_star_space": 16,
            "largest_cycle_control_spin_space": 32,
            "log_weight_spin_states_d2_d4_d6_total": 84,
            "log_weight_nonempty_even_sector_tables": 6,
        },
    )

    scope = data["scope"]
    check(
        "honest hidden-auxiliary and global scope",
        data["parameter"]["domain"] == "0<v<1"
        and data["parameter"]["benchmark_or_Kc_input"] is False
        and data["ordinary_real_CP_rank"]["rank"] == 2
        and data["spin_basis_positive_rank"]["rank"] == 2
        and any(
            "every even star degree d" in statement
            for statement in scope["proved"]
        )
        and "a no-go for arbitrary hidden auxiliary variables or unrestricted tensor networks"
        in scope["not_proved"]
        and "that a nonplanar or three-dimensional contraction of the local signatures is one global Pfaffian"
        in scope["not_proved"]
        and "global Pfaffian solvability or nonsolvability of the three-dimensional Ising partition function"
        in scope["not_proved"],
    )

    expected_check_names = {
        "exact_star_Walsh_signature",
        "all_even_star_log_Walsh_sign_controls",
        "bipartite_cycle_space_factorization_controls",
        "explicit_principal_Pfaffian_certificate",
        "exhaustive_degree_six_matchgate_identities",
        "signed_CP_rank_exactly_two",
        "three_by_three_PSD_parity_blocks",
        "all_nonempty_Cartesian_support_rectangles",
        "nonnegative_all_leg_CP_rank_32",
        "spin_basis_positive_rank_two_undoes_decimation",
        "scope_keeps_hidden_auxiliary_and_global_boundaries",
        "benchmark_not_used",
        "declared_resource_limits",
    }
    check(
        "stored producer checks",
        {row["name"] for row in artifact["checks"]}
        == expected_check_names
        and all(row["passed"] is True for row in artifact["checks"]),
    )
    check(
        "producer metadata and resource gate",
        meta["benchmark_used"] is False
        and meta["rss_limit_bytes"] == RSS_LIMIT_BYTES
        and meta["peak_rss_bytes"] < RSS_LIMIT_BYTES
        and meta["peak_rss_measurement"] == peak_rss_measurement()
        and meta["process_cpu_seconds"]
        < meta["process_cpu_budget_seconds"],
    )
    verifier_rss = peak_rss_bytes()
    check(
        "verifier RSS wall",
        verifier_rss < RSS_LIMIT_BYTES,
        f"{verifier_rss}/{RSS_LIMIT_BYTES} bytes via {peak_rss_measurement()}",
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: octahedral star matchgate and positivity barrier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
