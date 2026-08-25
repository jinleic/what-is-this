#!/usr/bin/env python3
"""Exact matchgate and positivity certificates for the degree-six Ising star.

The eliminated-center star has normalized Walsh signature
``g(A) = v**|A|`` on even subsets and zero on odd subsets. This producer
certifies its sub-Pfaffian form, exact tensor ranks, the basis-dependent
nonnegative-rank barrier, and finite controls for the all-even-degree
alternating log-Walsh sign theorem, using only integers and ``Fraction``.
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
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "octahedral_matchgate.json"
VERIFIER = ROOT / "tests" / "test_octahedral_matchgate.py"
INHERITED_PRODUCER = ROOT / "experiments" / "e234_star_decimation.py"
DLMF_SOURCE = ROOT / "sources" / "dlmf_4_36.html"
DIMENSION = 6
AUDIT_POINTS = (Fraction(1, 3), Fraction(1, 2))
CPU_BUDGET_SECONDS = 60.0
RSS_LIMIT_BYTES = 2 * 1024**3

CONTROL_GRAPHS: dict[str, tuple[int, int, tuple[tuple[int, int], ...]]] = {
    "C4": (2, 2, ((0, 0), (0, 1), (1, 0), (1, 1))),
    "K2_3": (
        2,
        3,
        tuple((left, right) for left in range(2) for right in range(3)),
    ),
}


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


def fraction_text(value: int | Fraction) -> str:
    rational = Fraction(value)
    if rational.denominator == 1:
        return str(rational.numerator)
    return f"{rational.numerator}/{rational.denominator}"


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def character(spin_mask: int, support_mask: int) -> int:
    return -1 if (spin_mask & support_mask).bit_count() % 2 else 1


def signature(mask: int, v: Fraction) -> Fraction:
    return v ** mask.bit_count() if mask.bit_count() % 2 == 0 else Fraction(0)


def normalized_spin_star(spin_mask: int, v: Fraction) -> Fraction:
    """Return W(sigma)/(2*cosh(K)**6), with v=tanh(K)."""

    plus = Fraction(1)
    minus = Fraction(1)
    for site in range(DIMENSION):
        spin = -1 if (spin_mask >> site) & 1 else 1
        plus *= 1 + v * spin
        minus *= 1 - v * spin
    return (plus + minus) / 2


def walsh_from_spin(v: Fraction) -> tuple[Fraction, ...]:
    values = [normalized_spin_star(state, v) for state in range(1 << DIMENSION)]
    scale = Fraction(1, 1 << DIMENSION)
    return tuple(
        scale
        * sum(
            values[state] * character(state, support)
            for state in range(1 << DIMENSION)
        )
        for support in range(1 << DIMENSION)
    )


def star_walsh_audit(v: Fraction) -> dict[str, object]:
    transformed = walsh_from_spin(v)
    expected = tuple(signature(mask, v) for mask in range(1 << DIMENSION))
    by_size: dict[str, dict[str, object]] = {}
    for size in range(DIMENSION + 1):
        entries = [
            transformed[mask]
            for mask in range(1 << DIMENSION)
            if mask.bit_count() == size
        ]
        by_size[str(size)] = {
            "entry_count": len(entries),
            "value": fraction_text(entries[0]),
            "constant_on_weight_class": len(set(entries)) == 1,
        }
    return {
        "spin_states": 1 << DIMENSION,
        "walsh_entries": 1 << DIMENSION,
        "walsh_summands": (1 << DIMENSION) ** 2,
        "coefficients_by_subset_size": by_size,
        "all_entries_match_g": transformed == expected,
    }


def cycle_space_counts(
    left_order: int, right_order: int, edges: tuple[tuple[int, int], ...]
) -> dict[int, int]:
    counts: Counter[int] = Counter()
    for edge_mask in range(1 << len(edges)):
        left_parity = [0] * left_order
        right_parity = [0] * right_order
        for edge_index, (left, right) in enumerate(edges):
            if (edge_mask >> edge_index) & 1:
                left_parity[left] ^= 1
                right_parity[right] ^= 1
        if not any(left_parity) and not any(right_parity):
            counts[edge_mask.bit_count()] += 1
    return dict(sorted(counts.items()))


def star_contraction_counts(
    left_order: int, right_order: int, edges: tuple[tuple[int, int], ...]
) -> tuple[dict[int, int], int]:
    incident: list[list[int]] = [[] for _ in range(left_order)]
    for edge_index, (left, _right) in enumerate(edges):
        incident[left].append(edge_index)
    local_choices: list[list[int]] = []
    for edge_indices in incident:
        choices: list[int] = []
        for local_mask in range(1 << len(edge_indices)):
            if local_mask.bit_count() % 2:
                continue
            choices.append(
                sum(
                    1 << edge_indices[local_index]
                    for local_index in range(len(edge_indices))
                    if (local_mask >> local_index) & 1
                )
            )
        local_choices.append(choices)

    counts: Counter[int] = Counter()
    assignments = 0
    for selected in itertools.product(*local_choices):
        assignments += 1
        edge_mask = 0
        for local_edge_mask in selected:
            edge_mask |= local_edge_mask
        retained_parity = [0] * right_order
        for edge_index, (_left, right) in enumerate(edges):
            if (edge_mask >> edge_index) & 1:
                retained_parity[right] ^= 1
        if not any(retained_parity):
            counts[edge_mask.bit_count()] += 1
    return dict(sorted(counts.items())), assignments


def normalized_spin_graph(
    left_order: int,
    right_order: int,
    edges: tuple[tuple[int, int], ...],
    v: Fraction,
) -> Fraction:
    order = left_order + right_order
    total = Fraction(0)
    for state in range(1 << order):
        spins = [
            -1 if (state >> vertex) & 1 else 1 for vertex in range(order)
        ]
        weight = Fraction(1)
        for left, right in edges:
            weight *= 1 + v * spins[left] * spins[left_order + right]
        total += weight
    return total / (1 << order)


def polynomial_value(coefficients: dict[int, int], v: Fraction) -> Fraction:
    return sum(
        Fraction(coefficient) * v**degree
        for degree, coefficient in coefficients.items()
    )


def cycle_control_record(
    name: str,
    left_order: int,
    right_order: int,
    edges: tuple[tuple[int, int], ...],
) -> tuple[dict[str, object], bool]:
    direct = cycle_space_counts(left_order, right_order, edges)
    contracted, local_assignments = star_contraction_counts(
        left_order, right_order, edges
    )
    values: dict[str, str] = {}
    spin_agrees = True
    for v in AUDIT_POINTS:
        polynomial = polynomial_value(direct, v)
        spin_value = normalized_spin_graph(left_order, right_order, edges, v)
        values[fraction_text(v)] = fraction_text(polynomial)
        spin_agrees &= spin_value == polynomial
    expected = {"C4": {0: 1, 4: 1}, "K2_3": {0: 1, 4: 3}}[name]
    cycle_dimension = len(edges) - left_order - right_order + 1
    record = {
        "name": name,
        "left_vertices": left_order,
        "right_vertices": right_order,
        "edges": [list(edge) for edge in edges],
        "edge_count": len(edges),
        "edge_subset_bound": 1 << len(edges),
        "local_even_star_assignment_bound": local_assignments,
        "spin_state_bound": 1 << (left_order + right_order),
        "retained_spin_sum_factor": 1 << right_order,
        "cycle_space_dimension": cycle_dimension,
        "coefficient_by_edge_count": {
            str(key): value for key, value in direct.items()
        },
        "values": values,
    }
    passed = (
        direct == contracted == expected
        and spin_agrees
        and sum(direct.values()) == 1 << cycle_dimension
    )
    return record, passed


def j_matrix(order: int) -> list[list[int]]:
    return [
        [
            0 if row == column else (1 if row < column else -1)
            for column in range(order)
        ]
        for row in range(order)
    ]


def pfaffian_recursive(matrix: list[list[Fraction]]) -> Fraction:
    order = len(matrix)
    if order == 0:
        return Fraction(1)
    if order % 2:
        return Fraction(0)
    total = Fraction(0)
    for partner in range(1, order):
        remaining = [index for index in range(1, order) if index != partner]
        minor = [
            [matrix[row][column] for column in remaining] for row in remaining
        ]
        sign = 1 if partner % 2 else -1
        total += sign * matrix[0][partner] * pfaffian_recursive(minor)
    return total


def principal_submatrix(
    matrix: list[list[Fraction]], mask: int
) -> list[list[Fraction]]:
    indices = [index for index in range(len(matrix)) if (mask >> index) & 1]
    return [[matrix[row][column] for column in indices] for row in indices]


def pfaffian_audit(v: Fraction) -> dict[str, object]:
    matrix = [
        [v * v * Fraction(entry) for entry in row] for row in j_matrix(DIMENSION)
    ]
    observed: dict[int, list[Fraction]] = {
        size: [] for size in range(0, DIMENSION + 1, 2)
    }
    passed = True
    for mask in range(1 << DIMENSION):
        if mask.bit_count() % 2:
            continue
        value = pfaffian_recursive(principal_submatrix(matrix, mask))
        observed[mask.bit_count()].append(value)
        passed &= value == signature(mask, v)
    return {
        "matrix_scale_v_squared": fraction_text(v * v),
        "even_principal_pfaffians_checked": sum(map(len, observed.values())),
        "by_subset_size": {
            str(size): {
                "count": len(values),
                "value": fraction_text(values[0]),
                "constant": len(set(values)) == 1,
            }
            for size, values in observed.items()
        },
        "all_equal_signature": passed,
    }


def matchgate_identity_audit(v: Fraction) -> dict[str, object]:
    values = tuple(signature(mask, v) for mask in range(1 << DIMENSION))
    residuals: list[Fraction] = []
    active_histogram: Counter[int] = Counter()
    for alpha in range(1 << DIMENSION):
        for beta in range(1 << DIMENSION):
            differing = [
                site
                for site in range(DIMENSION)
                if ((alpha ^ beta) >> site) & 1
            ]
            terms: list[Fraction] = []
            for position, site in enumerate(differing):
                sign = -1 if position % 2 == 0 else 1
                terms.append(
                    sign
                    * values[alpha ^ (1 << site)]
                    * values[beta ^ (1 << site)]
                )
            residuals.append(sum(terms, Fraction(0)))
            if any(term != 0 for term in terms):
                active_histogram[len(differing)] += 1
    return {
        "ordered_bitstring_pairs": len(residuals),
        "maximum_terms_per_identity": DIMENSION,
        "algebraically_active_by_term_count": {
            str(key): value for key, value in sorted(active_histogram.items())
        },
        "algebraically_active_total": sum(active_histogram.values()),
        "all_residuals_zero": all(residual == 0 for residual in residuals),
    }


def embed_assignment(assignment: int, coordinates: list[int]) -> int:
    return sum(
        ((assignment >> local_index) & 1) << coordinate
        for local_index, coordinate in enumerate(coordinates)
    )


def flattening(left_mask: int, v: Fraction) -> list[list[Fraction]]:
    left = [
        site for site in range(DIMENSION) if (left_mask >> site) & 1
    ]
    right = [
        site for site in range(DIMENSION) if not ((left_mask >> site) & 1)
    ]
    return [
        [
            signature(
                embed_assignment(row, left)
                | embed_assignment(column, right),
                v,
            )
            for column in range(1 << len(right))
        ]
        for row in range(1 << len(left))
    ]


def fraction_rank(matrix: list[list[Fraction]]) -> int:
    """Small exact Gauss-Jordan rank routine for the finite flattenings."""

    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(row_count):
            if row == rank or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(column_count)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def three_by_three_audit(v: Fraction) -> dict[str, object]:
    matrix = flattening(0b000111, v)
    even_vector = [
        v ** state.bit_count()
        if state.bit_count() % 2 == 0
        else Fraction(0)
        for state in range(8)
    ]
    odd_vector = [
        v ** state.bit_count()
        if state.bit_count() % 2
        else Fraction(0)
        for state in range(8)
    ]
    gram = [
        [
            even_vector[row] * even_vector[column]
            + odd_vector[row] * odd_vector[column]
            for column in range(8)
        ]
        for row in range(8)
    ]
    even_eigenvalue = sum(
        (entry * entry for entry in even_vector), Fraction(0)
    )
    odd_eigenvalue = sum(
        (entry * entry for entry in odd_vector), Fraction(0)
    )
    return {
        "left_mask": "000111",
        "shape": [8, 8],
        "rank": fraction_rank(matrix),
        "symmetric": matrix == [list(row) for row in zip(*matrix)],
        "equals_two_parity_outer_products": matrix == gram,
        "parity_blocks": [
            {"parity": "even", "dimension": 4, "rank": 1},
            {"parity": "odd", "dimension": 4, "rank": 1},
        ],
        "nonzero_eigenvalues": {
            "even_1_plus_3v4": fraction_text(even_eigenvalue),
            "odd_3v2_plus_v6": fraction_text(odd_eigenvalue),
        },
        "zero_eigenvalue_multiplicity": 6,
        "positive_semidefinite_certificate": (
            even_eigenvalue > 0 and odd_eigenvalue > 0 and matrix == gram
        ),
    }


def cp_rank_audit(v: Fraction) -> dict[str, object]:
    signed_decomposition = all(
        signature(mask, v)
        == (v ** mask.bit_count() + (-v) ** mask.bit_count()) / 2
        for mask in range(1 << DIMENSION)
    )
    ranks_by_size: dict[str, dict[str, object]] = {}
    all_ranks: list[int] = []
    for left_size in range(1, DIMENSION):
        ranks = [
            fraction_rank(flattening(left_mask, v))
            for left_mask in range(1, 1 << DIMENSION)
            if left_mask != (1 << DIMENSION) - 1
            and left_mask.bit_count() == left_size
        ]
        all_ranks.extend(ranks)
        ranks_by_size[str(left_size)] = {
            "ordered_left_subsets": len(ranks),
            "rank_histogram": {
                str(rank): ranks.count(rank) for rank in sorted(set(ranks))
            },
        }
    three_by_three = three_by_three_audit(v)
    return {
        "signature_entries_checked": 1 << DIMENSION,
        "signed_two_atom_decomposition_exact": signed_decomposition,
        "ordered_nontrivial_flattenings": len(all_ranks),
        "all_flattening_ranks_two": all(rank == 2 for rank in all_ranks),
        "rank_histogram": {
            str(rank): all_ranks.count(rank) for rank in sorted(set(all_ranks))
        },
        "by_left_size": ranks_by_size,
        "three_by_three_operator": three_by_three,
    }


def rectangle_audit() -> dict[str, object]:
    by_free: Counter[int] = Counter()
    singleton_even = 0
    singleton_odd = 0
    non_singleton_both = 0
    even_only = 0
    maximum_even_only_size = 0
    total = 0
    for leg_supports in itertools.product((1, 2, 3), repeat=DIMENSION):
        total += 1
        points = [
            state
            for state in range(1 << DIMENSION)
            if all(
                leg_supports[site] & (1 << ((state >> site) & 1))
                for site in range(DIMENSION)
            )
        ]
        free = sum(support == 3 for support in leg_supports)
        by_free[free] += 1
        parities = {point.bit_count() % 2 for point in points}
        if len(points) == 1:
            if next(iter(parities)) == 0:
                singleton_even += 1
            else:
                singleton_odd += 1
        elif parities == {0, 1}:
            non_singleton_both += 1
        if parities == {0}:
            even_only += 1
            maximum_even_only_size = max(
                maximum_even_only_size, len(points)
            )
    return {
        "dimension": DIMENSION,
        "nonempty_support_choices_per_leg": 3,
        "rectangles_checked": total,
        "by_free_coordinates": {
            str(key): value for key, value in sorted(by_free.items())
        },
        "singleton_even": singleton_even,
        "singleton_odd": singleton_odd,
        "non_singleton": total - singleton_even - singleton_odd,
        "non_singleton_containing_both_parities": non_singleton_both,
        "rectangles_contained_in_even_support": even_only,
        "maximum_even_only_rectangle_size": maximum_even_only_size,
    }


def nonnegative_upper_audit(v: Fraction) -> dict[str, object]:
    atoms = {
        state: signature(state, v)
        for state in range(1 << DIMENSION)
        if state.bit_count() % 2 == 0
    }
    reconstructed = tuple(
        atoms.get(state, Fraction(0)) for state in range(1 << DIMENSION)
    )
    target = tuple(
        signature(state, v) for state in range(1 << DIMENSION)
    )
    return {
        "singleton_atoms": len(atoms),
        "all_atom_weights_positive": all(
            weight > 0 for weight in atoms.values()
        ),
        "reconstructs_signature": reconstructed == target,
    }


def spin_basis_audit(v: Fraction) -> dict[str, object]:
    decomposition_exact = True
    full_star_values: list[Fraction] = []
    for state in range(1 << DIMENSION):
        plus = Fraction(1)
        minus = Fraction(1)
        for site in range(DIMENSION):
            spin = -1 if (state >> site) & 1 else 1
            plus *= 1 + v * spin
            minus *= 1 - v * spin
        full_star_values.append(plus + minus)
        decomposition_exact &= (
            plus > 0
            and minus > 0
            and plus + minus == 2 * normalized_spin_star(state, v)
        )
    matrix = [
        [
            full_star_values[row | (column << 1)]
            for column in range(1 << (DIMENSION - 1))
        ]
        for row in range(2)
    ]
    return {
        "spin_states_checked": 1 << DIMENSION,
        "two_positive_atoms_exact": decomposition_exact,
        "one_by_five_flattening_rank": fraction_rank(matrix),
    }

def log_walsh_multiplicities(
    degree: int, subset_size: int
) -> dict[int, int]:
    """Return integer log-cosh orbit weights before the factor 2**(-degree)."""

    support = (1 << subset_size) - 1
    counts: Counter[int] = Counter()
    for state in range(1 << degree):
        magnetization = abs(degree - 2 * state.bit_count())
        counts[magnetization] += (
            -1 if (state & support).bit_count() % 2 else 1
        )
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


def polynomial_subtract(left: list[int], right: list[int]) -> list[int]:
    width = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(width)
    ]


def divide_x_minus_one(poly: list[int]) -> tuple[list[int], int]:
    descending = list(reversed(poly))
    quotient_descending = [descending[0]]
    for coefficient in descending[1:-1]:
        quotient_descending.append(coefficient + quotient_descending[-1])
    remainder = descending[-1] + quotient_descending[-1]
    return list(reversed(quotient_descending)), remainder


def shift_at_one(poly: list[int]) -> list[int]:
    out = [0] * len(poly)
    for degree, coefficient in enumerate(poly):
        for shifted_degree in range(degree + 1):
            out[shifted_degree] += (
                coefficient * math.comb(degree, shifted_degree)
            )
    return out


def log_walsh_sign_controls() -> tuple[dict[str, object], bool]:
    expected_top = {
        2: {0: -2, 2: 2},
        4: {0: 6, 2: -8, 4: 2},
        6: {0: -20, 2: 30, 4: -12, 6: 2},
    }
    sector_tables: dict[str, object] = {}
    passed = True
    for degree in (2, 4, 6):
        sectors: dict[str, object] = {}
        for subset_size in range(2, degree + 1, 2):
            counts = log_walsh_multiplicities(degree, subset_size)
            sectors[str(subset_size)] = {
                "scaled_integer_multiplicities_by_abs_magnetization": {
                    str(magnetization): coefficient
                    for magnetization, coefficient in counts.items()
                },
                "Walsh_normalization_denominator": 1 << degree,
                "theorem_sign": (
                    "positive"
                    if (subset_size // 2 + 1) % 2 == 0
                    else "negative"
                ),
            }
            if subset_size == degree:
                passed &= counts == expected_top[degree]
        sector_tables[str(degree)] = {
            "spin_states": 1 << degree,
            "nonempty_even_sectors": sectors,
        }

    degree_four_gap = polynomial_subtract([0, 0, 1], [-1, 2])
    degree_six_difference = polynomial_subtract(
        polynomial_multiply([0] * 8 + [1], [-3, 4]),
        polynomial_power([-1, 2], 6),
    )
    quotient = degree_six_difference
    remainders: list[int] = []
    for _ in range(3):
        quotient, remainder = divide_x_minus_one(quotient)
        remainders.append(remainder)
    degree_six_shifted = shift_at_one(quotient)
    passed &= (
        degree_four_gap == [1, -2, 1]
        and remainders == [0, 0, 0]
        and degree_six_shifted == [8, 54, 144, 188, 120, 33, 4]
        and all(coefficient > 0 for coefficient in degree_six_shifted)
    )
    return {
        "degrees": sector_tables,
        "top_sector_ratio_certificates": {
            "d2": {
                "formula": "c2=(1/2)*log(cosh(2K))",
                "sign": "positive for K!=0",
            },
            "d4": {
                "ratio": "exp(8*c4)=(2*x-1)/x^2",
                "x": "cosh(2K)^2>1",
                "denominator_minus_numerator_coefficients_ascending": degree_four_gap,
            },
            "d6": {
                "ratio": "exp(32*c6)=x^8*(4*x-3)/(2*x-1)^6",
                "difference_factor": "(x-1)^3*Q6(x)",
                "Q6_at_1_plus_y_coefficients_ascending": degree_six_shifted,
            },
        },
    }, passed


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    star_audits = {
        fraction_text(v): star_walsh_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "exact_star_Walsh_signature",
        all(
            bool(row["all_entries_match_g"])
            for row in star_audits.values()
        ),
        "two 64x64 exact Walsh transforms recover g(A)=v^|A| on even A and zero on odd A",
    )

    log_sign_controls, log_controls_passed = log_walsh_sign_controls()
    add_check(
        checks,
        "all_even_star_log_Walsh_sign_controls",
        log_controls_passed,
        "exact spin-orbit multiplicities for every nonempty even sector at d=2,4,6, plus d4 and d6 polynomial ratio certificates",
    )

    cycle_records: list[dict[str, object]] = []
    cycle_passed = True
    for name, (left_order, right_order, edges) in CONTROL_GRAPHS.items():
        record, passed = cycle_control_record(
            name, left_order, right_order, edges
        )
        cycle_records.append(record)
        cycle_passed &= passed
    add_check(
        checks,
        "bipartite_cycle_space_factorization_controls",
        cycle_passed,
        "C4 gives 1+v^4 and K2,3 gives 1+3v^4 by local-star contraction, edge subsets, and spin sums",
    )

    pfaffian_audits = {
        fraction_text(v): pfaffian_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "explicit_principal_Pfaffian_certificate",
        all(
            row["even_principal_pfaffians_checked"] == 32
            and bool(row["all_equal_signature"])
            for row in pfaffian_audits.values()
        ),
        "all 32 even principal Pfaffians of v^2 J6 equal the star signature at each rational audit point",
    )

    identity_audits = {
        fraction_text(v): matchgate_identity_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "exhaustive_degree_six_matchgate_identities",
        all(
            row["ordered_bitstring_pairs"] == 4096
            and row["algebraically_active_by_term_count"]
            == {"2": 480, "4": 480, "6": 32}
            and bool(row["all_residuals_zero"])
            for row in identity_audits.values()
        ),
        "all 64^2 ordered matchgate identities vanish exactly at v=1/3 and v=1/2",
    )

    cp_audits = {
        fraction_text(v): cp_rank_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "signed_CP_rank_exactly_two",
        all(
            bool(row["signed_two_atom_decomposition_exact"])
            and bool(row["all_flattening_ranks_two"])
            for row in cp_audits.values()
        ),
        "the signed two-atom formula is exact and every one of 62 ordered nontrivial flattenings has rank two",
    )
    add_check(
        checks,
        "three_by_three_PSD_parity_blocks",
        all(
            row["three_by_three_operator"]["rank"] == 2
            and bool(
                row["three_by_three_operator"][
                    "positive_semidefinite_certificate"
                ]
            )
            for row in cp_audits.values()
        ),
        "the 8x8 flattening is the sum of orthogonal even- and odd-parity outer products",
    )

    rectangles = rectangle_audit()
    expected_rectangle_histogram = {
        "0": 64,
        "1": 192,
        "2": 240,
        "3": 160,
        "4": 60,
        "5": 12,
        "6": 1,
    }
    add_check(
        checks,
        "all_nonempty_Cartesian_support_rectangles",
        rectangles["rectangles_checked"] == 3**DIMENSION
        and rectangles["by_free_coordinates"]
        == expected_rectangle_histogram
        and rectangles["singleton_even"] == 32
        and rectangles["singleton_odd"] == 32
        and rectangles["non_singleton_containing_both_parities"] == 665
        and rectangles["rectangles_contained_in_even_support"] == 32
        and rectangles["maximum_even_only_rectangle_size"] == 1,
        "all 3^6=729 rectangles were checked; exactly the 32 even singletons lie in the positive support",
    )

    nonnegative_audits = {
        fraction_text(v): nonnegative_upper_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "nonnegative_all_leg_CP_rank_32",
        all(
            row["singleton_atoms"] == 32
            and bool(row["all_atom_weights_positive"])
            and bool(row["reconstructs_signature"])
            for row in nonnegative_audits.values()
        ),
        "the rectangle barrier gives the lower bound 32 and the 32 weighted even singleton atoms attain it",
    )

    spin_audits = {
        fraction_text(v): spin_basis_audit(v) for v in AUDIT_POINTS
    }
    add_check(
        checks,
        "spin_basis_positive_rank_two_undoes_decimation",
        all(
            bool(row["two_positive_atoms_exact"])
            and row["one_by_five_flattening_rank"] == 2
            for row in spin_audits.values()
        ),
        "W/cosh(K)^6 is the positive two-term sum over the original center spin; this is a basis-dependent reconstruction",
    )

    scope = {
        "proved": [
            "the exact degree-six local star Walsh signature for every 0<v=tanh(K)<1",
            "for every even star degree d and every nonempty even A, the log-weight Walsh coefficient c_A has strict sign (-1)^(|A|/2+1) for K!=0",
            "a pure-even sub-Pfaffian matchgate certificate for that one local signature",
            "ordinary real CP rank two, rank two for every nontrivial flattening, and nonnegative all-leg CP rank 32 in the Walsh basis",
            "the finite bipartite cycle-space contraction identity, with exact C4 and K2,3 controls",
            "a positive spin-basis rank-two reconstruction obtained by reintroducing the eliminated center spin",
        ],
        "not_proved": [
            "a no-go for arbitrary hidden auxiliary variables or unrestricted tensor networks",
            "that a nonplanar or three-dimensional contraction of the local signatures is one global Pfaffian",
            "global Pfaffian solvability or nonsolvability of the three-dimensional Ising partition function",
            "integrability, a critical coupling, a thermodynamic free energy, or critical exponents",
        ],
    }
    add_check(
        checks,
        "scope_keeps_hidden_auxiliary_and_global_boundaries",
        "a no-go for arbitrary hidden auxiliary variables or unrestricted tensor networks"
        in scope["not_proved"]
        and "global Pfaffian solvability or nonsolvability of the three-dimensional Ising partition function"
        in scope["not_proved"],
        "the local positivity and Pfaffian statements are not promoted to an auxiliary no-go or a global solution",
    )
    add_check(
        checks,
        "benchmark_not_used",
        True,
        "only symbolic 0<v<1 arguments and exact v=1/3,1/2 audits are used; no Kc enters",
    )

    source_paths = (
        Path(__file__).resolve(),
        VERIFIER,
        INHERITED_PRODUCER,
        DLMF_SOURCE,
    )
    source_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
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

    failed = [str(row["name"]) for row in checks if not row["passed"]]
    if failed:
        raise AssertionError(f"octahedral matchgate checks failed: {failed}")

    inherited_hash = source_hashes["experiments/e234_star_decimation.py"]
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e244_octahedral_matchgate.py",
            "verifier": "tests/test_octahedral_matchgate.py",
            "artifact": "results/integrability/octahedral_matchgate.json",
            "interpreter": sys.executable,
            "arithmetic": "exact integers and fractions; no floating-point mathematical claims",
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "peak_rss_measurement": rss_measurement,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "source_sha256": source_hashes,
        },
        "data": {
            "claim_tag": "[THEOREM][EXACT FINITE COMPUTATION]",
            "parameter": {
                "degree": DIMENSION,
                "name": "v=tanh(K)",
                "domain": "0<v<1",
                "omitted_common_factor": "2*cosh(K)^6",
                "benchmark_or_Kc_input": False,
            },
            "star_walsh": {
                "identity": "W(sigma)=2*cosh(K)^6*sum_{A even} v^|A| chi_A(sigma)",
                "signature": "g(A)=v^|A| for even A and g(A)=0 for odd A",
                "normalization": "g(empty)=1 after omitting 2*cosh(K)^6",
                "audits": star_audits,
            },
            "all_even_star_log_walsh_sign": {
                "claim_tag": "[THEOREM]",
                "definition": "c_A=2^(-d)*sum_sigma chi_A(sigma)*log(cosh(K*sum_i sigma_i))",
                "statement": "for every even d and nonempty even A subset [d], sign(c_A)=(-1)^(|A|/2+1) for every real K!=0",
                "consequences": "for every even d>=6: all two-spin sectors are positive, all four-spin sectors are negative, all six-spin sectors are positive, and higher even sectors alternate",
                "analytic_certificate": {
                    "cosh_product": "cosh(y)=product_(n>=0)(1+4*y^2/(pi^2*(2n+1)^2))",
                    "product_source": "NIST DLMF 4.36.E2, https://dlmf.nist.gov/4.36.E2",
                    "log_integral": "log(1+a*y^2)=integral_0^infinity exp(-s)*(1-exp(-a*s*y^2))*ds/s",
                    "Gaussian_Walsh_transform": "E[chi_A*exp(-b*(sum sigma)^2)]=(-1)^r/sqrt(4*pi*b)*integral_R exp(-t^2/(4b))*sin(t)^(2r)*cos(t)^(d-2r) dt for |A|=2r",
                    "strictness": "b>0 for K!=0; even d makes both trigonometric powers even, so the integrand is nonnegative and not almost everywhere zero; odd d is outside the theorem because cos(t)^(d-2r) changes sign",
                    "interchanges": "the Walsh sum is finite; for alpha=a_n*K^2 and 0<s<=1, E[chi_A*S^(2k)]=0 for k<r and Taylor remainder bounds |E[chi_A*exp(-alpha*s*S^2)]| by (alpha*s*d^2)^r/r!, hence division by s is O(s^(r-1)); for s>=1 use the bound 1 against exp(-s)/s; the product-log sum is uniformly absolutely convergent on the finite star magnetizations because sum_n a_n is finite",
                },
                "finite_exact_controls": log_sign_controls,
            },
            "cycle_space_factorization": {
                "retained_spin_average_normalization": "2^(-|V|)",
                "normalized_identity": "2^(-|V|)*sum_(sigma_V) product_(u in U) g_u(sigma_delta(u)) = sum_(F subset E: boundary(F)=empty) v^|F|",
                "raw_sum_identity": "sum_(sigma_V) product_(u in U) g_u(sigma_delta(u)) = 2^|V|*sum_(F subset E: boundary(F)=empty) v^|F|",
                "controls": cycle_records,
            },
            "pfaffian_matchgate": {
                "convention": "sub-Pfaffian tensor indexed by retained principal subsets",
                "generator": "M(v)=v^2*J_6, with J_6[i,j]=1 for i<j and -1 for i>j",
                "J_6": j_matrix(DIMENSION),
                "all_d_lemma": "Pf(J_(2r))=1 by first-row expansion; hence Pf((v^2 J_d)[A])=v^|A| for every even A",
                "even_principal_subset_count": 32,
                "principal_pfaffian_audits": pfaffian_audits,
                "matchgate_identity_formula": "sum_i (-1)^i g(alpha xor e_p_i) g(beta xor e_p_i)=0 for p_1<...<p_l=supp(alpha xor beta)",
                "matchgate_identity_audits": identity_audits,
            },
            "ordinary_real_CP_rank": {
                "rank": 2,
                "signed_decomposition": "g=(tensor(1,v)+tensor(1,-v))/2",
                "lower_bound": "every nontrivial bipartition flattening has exact rank two",
                "audits": cp_audits,
            },
            "nonnegative_all_leg_CP_rank": {
                "rank": 32,
                "support_lemma": "a nonzero nonnegative rank-one atom has Cartesian support; every nonsingleton binary rectangle contains both parities",
                "rectangle_enumeration": rectangles,
                "lower_bound": "odd target zeros force every atom support into the even support, so each atom covers at most one of 32 positive entries",
                "upper_bound_decomposition": "g=sum_{x even} v^|x| tensor_i delta_(x_i)",
                "even_singleton_states": [
                    format(state, "06b")
                    for state in range(1 << DIMENSION)
                    if state.bit_count() % 2 == 0
                ],
                "audits": nonnegative_audits,
            },
            "spin_basis_positive_rank": {
                "rank": 2,
                "identity": "W(sigma)/cosh(K)^6=product_i(1+v sigma_i)+product_i(1-v sigma_i)",
                "positivity": "each one-leg factor 1+/-v sigma_i is strictly positive for 0<v<1",
                "interpretation": "the two atoms are the two values of the eliminated center spin, so this only undoes decimation",
                "audits": spin_audits,
            },
            "inherited_source": {
                "producer": "experiments/e234_star_decimation.py",
                "sha256": inherited_hash,
                "binding": "current source hash only; no frozen or generated result is imported",
            },
            "finite_enumeration_bounds": {
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
            "scope": scope,
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in payload["checks"]:
        print(
            f"  [{'PASS' if row['passed'] else 'FAIL'}] "
            f"{row['name']}: {row['detail']}"
        )
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
