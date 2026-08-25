#!/usr/bin/env python3
"""Exact support-defect and Smith-block certificate for the e237 W-law cores.

The producer rebuilds the literal residual minors at (L,m)=(8,4),(9,4),
(9,5) from the stored exact annihilator bases and the deterministic e237 leaf
rule.  It then decodes every surviving orbit coordinate, performs direct and
Burnside orbit counts, and computes exact Smith forms for the resulting
support-defect blocks.  It writes only after provenance, semantics, scope, CPU,
and RSS gates pass.
"""

from __future__ import annotations

import ctypes
import hashlib
import heapq
import itertools
import json
import math
import platform
import resource
import sys
import time
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "algebra_growth" / "wlaw_support_defect.json"
VERIFIER = ROOT / "tests" / "test_wlaw_support_defect.py"
PROOF = ROOT / "proofs" / "wlaw_support_defect.md"
E230_PRODUCER = ROOT / "experiments" / "e230_wlaw_rank_mechanism.py"
E237_PRODUCER = ROOT / "experiments" / "e237_wlaw_exact_sequence.py"
E230_ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw_rank_mechanism.json"
E237_ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw_exact_sequence.json"
L8_SOURCE = ROOT / "results" / "ladder" / "l8_saturation.json"
W9_METADATA = ROOT / "results" / "ladder" / "w9_saturation.json"
W9_BASIS = ROOT / "results" / "ladder" / "w9_basis.json"
PYPROJECT = ROOT / "pyproject.toml"
UV_LOCK = ROOT / "uv.lock"

TARGETS = ((8, 4), (9, 4), (9, 5))
GROUP_ELEMENTS = ("id", "tau", "rho", "tau_rho")
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
EXPECTED_UPSTREAM_SHA256 = {
    "experiments/e230_wlaw_rank_mechanism.py": "c2c8a3644777796f5f5b15e1643d87918fdfe45a0f911c2bbb91ecc3c2f20ce4",
    "experiments/e237_wlaw_exact_sequence.py": "4eb19f853d3e6acafad5f52bc2771d1cb87ddda8bb9044b14920e3ca60a52bbd",
    "results/algebra_growth/wlaw_rank_mechanism.json": "f277440fa05be9685f3b2e02b23b0f08cbdd1786ac8918ba4c8f407ca971480e",
    "results/algebra_growth/wlaw_exact_sequence.json": "50a9685c22679369b413fd9386319875fb02b231c77f8dc18ab5f346dc19647a",
    "results/ladder/l8_saturation.json": "d9109f4c11102f948dfbaa3b0972a98f02023fc463a26916d81a094aa7b74ccc",
    "results/ladder/w9_saturation.json": "06d680faf20f1c865f9736f145ad10cef2d1c4f760a744ec196180f481fbbe7e",
    "results/ladder/w9_basis.json": "f2b1b202b1e05ceee91856ab3ed832fbe58dabb97fea42f96da88e739e285413",
    "pyproject.toml": "3b779706037618fec9d5d9da29987e45292b1fe48b3856f37caa80a3499c36f6",
    "uv.lock": "4726afabb3d04a24194d5196ebc5b5331ae1bd2ad653987df10d5f80aaee3020",
}
EXPECTED_DETERMINANTS = {
    (8, 4): -(2**32) * 3,
    (9, 4): -(2**288) * 3**5,
    (9, 5): -(2**288) * 3**5,
}
EXPECTED_FAMILY_COUNTS = {
    (8, 4): {"(6,2)": 16, "(8,0)": 1},
    (9, 4): {"(6,2)": 128, "(8,0)": 5},
    (9, 5): {"(6,2)": 128, "(8,0)": 5},
}
EXPECTED_FULL_SMITH = {
    (8, 4): {2: 4, 4: 11, 8: 1, 24: 1},
    (9, 4): {2: 4, 4: 103, 8: 21, 24: 5},
    (9, 5): {2: 4, 4: 103, 8: 21, 24: 5},
}
SCOPE_CLAIM = (
    "the elementary rung-defect inequality for every even-particle coordinate, "
    "and the exact support, Burnside, particle-hole, determinant, and Smith-block "
    "classification of only the stored e237 cores (8,4), (9,4), and (9,5); "
    "no recurrence or rank statement for L>=10"
)
EXPECTED_PYTHON_VERSION = "3.14.3"
EXPECTED_SYMPY_VERSION = "1.14.0"

MatrixZ = list[list[int]]


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


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def tau_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 1) << (2 * rung + 1)
        image |= ((config >> (2 * rung + 1)) & 1) << (2 * rung)
    return image


def rho_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return image


def group_action(config: int, L: int, element: str) -> int:
    if element == "id":
        return config
    if element == "tau":
        return tau_config(config, L)
    if element == "rho":
        return rho_config(config, L)
    if element == "tau_rho":
        return tau_config(rho_config(config, L), L)
    raise ValueError(f"unknown group element {element}")


def orbit(config: int, L: int) -> set[int]:
    return {group_action(config, L, element) for element in GROUP_ELEMENTS}


def canonical_representative(config: int, L: int) -> int:
    return min(orbit(config, L))


def complement_config(config: int, L: int) -> int:
    return config ^ ((1 << (2 * L)) - 1)


def decode_sets(config: int, L: int) -> dict[str, list[int]]:
    decoded = {"T": [], "B": [], "E": [], "D": []}
    labels = ("E", "T", "B", "D")
    for rung in range(L):
        decoded[labels[(config >> (2 * rung)) & 3]].append(rung)
    return decoded


def decode_record(config: int, L: int, m: int) -> dict[str, object]:
    decoded = decode_sets(config, L)
    counts = {label: len(rungs) for label, rungs in decoded.items()}
    q = counts["E"] + counts["D"]
    signed_defect = counts["E"] - counts["D"]
    single_counts = sorted((counts["T"], counts["B"]), reverse=True)
    defect_rungs = decoded["E"] + decoded["D"]
    defect_position = (
        min(defect_rungs[0], L - 1 - defect_rungs[0])
        if len(defect_rungs) == 1
        else None
    )
    return {
        "representative": config,
        "T": decoded["T"],
        "B": decoded["B"],
        "E": decoded["E"],
        "D": decoded["D"],
        "counts": counts,
        "particle_count": counts["T"] + counts["B"] + 2 * counts["D"],
        "absolute_leg_imbalance": abs(counts["T"] - counts["B"]),
        "signed_defect_E_minus_D": signed_defect,
        "q_E_plus_D": q,
        "packing_boundary": q == abs(L - 2 * m),
        "no_cancelling_empty_double_pair": min(counts["E"], counts["D"]) == 0,
        "family": f"({single_counts[0]},{single_counts[1]})",
        "defect_position_mod_reflection": defect_position,
    }


def configuration_orbits(L: int, sector: int) -> list[dict[str, int]]:
    seen: set[int] = set()
    records: list[dict[str, int]] = []
    for config in range(1 << (2 * L)):
        if config.bit_count() != sector or config in seen:
            continue
        values = orbit(config, L)
        seen.update(values)
        representative = min(values)
        decoded = decode_sets(representative, L)
        records.append(
            {
                "representative": representative,
                "orbit_size": len(values),
                "particle_sector": sector,
                "absolute_leg_imbalance": abs(
                    len(decoded["T"]) - len(decoded["B"])
                ),
            }
        )
    return records


def canonical_basis_sha256(basis: list[list[dict[str, int]]]) -> str:
    return canonical_sha256(basis)


def source_exact(record: dict[str, object]) -> bool:
    validation = record["validation"]
    return (
        len(record["basis"]) == int(record["W_dim"])
        and sum(int(value) for value in record["W_sectors"].values())
        == int(record["W_dim"])
        and int(record["cyclic_Q_dim"]) + int(record["W_dim"])
        == int(record["K_dim"])
        and bool(validation["A_invariant_over_Q"])
        and bool(validation["basis_independent_over_Q"])
        and bool(validation["four_B_invariant_over_Q"])
        and bool(validation["sector_homogeneous"])
        and bool(validation["vacuum_orthogonal"])
        and not validation["outside_image_indices"]
    )


def weighted_restriction_rows(
    L: int,
    m: int,
    orbit_records: list[dict[str, int]],
    annihilator_record: dict[str, object],
) -> tuple[list[dict[str, int]], list[dict[str, object]]]:
    sector = 2 * m
    high_orbits = [
        row
        for row in orbit_records
        if row["absolute_leg_imbalance"] >= 4
    ]
    high_index = {
        row["representative"]: index for index, row in enumerate(high_orbits)
    }
    high_sizes = {
        row["representative"]: row["orbit_size"] for row in high_orbits
    }
    rows: list[dict[str, object]] = []
    for basis_index, vector in enumerate(annihilator_record["basis"]):
        if not vector or int(vector[0]["representative"]).bit_count() != sector:
            continue
        entries: dict[int, int] = {}
        for term in vector:
            representative = int(term["representative"])
            column = high_index.get(representative)
            if column is not None:
                entries[column] = (
                    int(term["coefficient"]) * high_sizes[representative]
                )
        rows.append({"basis_index": basis_index, "entries": entries})
    return high_orbits, rows


def leaf_certificate(
    rows: list[dict[str, object]], high_orbits: list[dict[str, int]]
) -> tuple[list[dict[str, int]], set[int], set[int]]:
    supports = [set(row["entries"]) for row in rows]
    by_column: list[list[int]] = [[] for _ in high_orbits]
    for row_index, support in enumerate(supports):
        for column in support:
            by_column[column].append(row_index)
    remaining = set(range(len(high_orbits)))
    used_rows: set[int] = set()
    queue = [index for index, support in enumerate(supports) if len(support) == 1]
    heapq.heapify(queue)
    selected: list[dict[str, int]] = []
    while queue:
        row_index = heapq.heappop(queue)
        if row_index in used_rows:
            continue
        active = supports[row_index] & remaining
        if len(active) != 1:
            continue
        pivot = next(iter(active))
        entries = rows[row_index]["entries"]
        selected.append(
            {
                "source_basis_index": int(rows[row_index]["basis_index"]),
                "pivot_representative": int(
                    high_orbits[pivot]["representative"]
                ),
                "diagonal_entry": int(entries[pivot]),
            }
        )
        used_rows.add(row_index)
        remaining.remove(pivot)
        for neighbour in by_column[pivot]:
            if neighbour not in used_rows and len(
                supports[neighbour] & remaining
            ) == 1:
                heapq.heappush(queue, neighbour)
    return selected, remaining, used_rows


def dense_submatrix(
    rows: list[dict[str, object]], row_indices: Iterable[int], columns: list[int]
) -> MatrixZ:
    position = {column: index for index, column in enumerate(columns)}
    result: MatrixZ = []
    for row_index in row_indices:
        dense = [0] * len(columns)
        for column, coefficient in rows[row_index]["entries"].items():
            target = position.get(column)
            if target is not None:
                dense[target] = int(coefficient)
        result.append(dense)
    return result


def bareiss_determinant(matrix: MatrixZ) -> int:
    dimension = len(matrix)
    if dimension == 0:
        return 1
    if any(len(row) != dimension for row in matrix):
        raise ValueError("Bareiss determinant requires a square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous_pivot = 1
    for column in range(dimension - 1):
        pivot_row = next(
            (row for row in range(column, dimension) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        pivot_data = work[column]
        for row_index in range(column + 1, dimension):
            row = work[row_index]
            entry = row[column]
            for target in range(column + 1, dimension):
                numerator = row[target] * pivot - entry * pivot_data[target]
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise AssertionError("non-exact Bareiss division")
                row[target] = quotient
            row[column] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


def integer_factorization(value: int) -> dict[str, int]:
    sign = -1 if value < 0 else 1
    residual = abs(value)
    exponent_2 = 0
    exponent_3 = 0
    while residual and residual % 2 == 0:
        residual //= 2
        exponent_2 += 1
    while residual and residual % 3 == 0:
        residual //= 3
        exponent_3 += 1
    return {
        "sign": sign,
        "exponent_2": exponent_2,
        "exponent_3": exponent_3,
        "residual_after_2_and_3": residual,
    }


def smith_record(matrix: MatrixZ) -> dict[str, object]:
    form = smith_normal_form(Matrix(matrix), domain=ZZ)
    diagonal = [
        abs(int(form[index, index]))
        for index in range(min(form.rows, form.cols))
        if form[index, index]
    ]
    divisibility_chain = all(
        right % left == 0 for left, right in zip(diagonal, diagonal[1:])
    )
    multiplicities = Counter(diagonal)
    determinant_abs = math.prod(diagonal) if len(matrix) == len(matrix[0]) else None
    return {
        "rank_Z_equals_rank_Q": len(diagonal),
        "diagonal": diagonal,
        "multiplicities": {
            str(value): multiplicities[value] for value in sorted(multiplicities)
        },
        "divisibility_chain": divisibility_chain,
        "absolute_determinant_from_diagonal": determinant_abs,
    }


def sparse_rows(matrix: MatrixZ) -> list[list[list[int]]]:
    return [
        [[column, value] for column, value in enumerate(row) if value]
        for row in matrix
    ]


def matrix_sha256(matrix: MatrixZ) -> str:
    column_count = len(matrix[0]) if matrix else 0
    return canonical_sha256(
        {"shape": [len(matrix), column_count], "sparse_rows": sparse_rows(matrix)}
    )


def matrix_record(matrix: MatrixZ, include_sparse_rows: bool = False) -> dict[str, object]:
    column_count = len(matrix[0]) if matrix else 0
    determinant = bareiss_determinant(matrix) if len(matrix) == column_count else None
    smith = smith_record(matrix)
    record: dict[str, object] = {
        "shape": [len(matrix), column_count],
        "nnz": sum(value != 0 for row in matrix for value in row),
        "sha256": matrix_sha256(matrix),
        "rank_exact_Q": smith["rank_Z_equals_rank_Q"],
        "smith_normal_form": smith,
    }
    if determinant is not None:
        record["determinant_Z"] = determinant
        record["determinant_factorization"] = integer_factorization(determinant)
        record["determinant_matches_smith"] = (
            abs(determinant) == smith["absolute_determinant_from_diagonal"]
        )
    if include_sparse_rows:
        record["sparse_rows"] = sparse_rows(matrix)
    return record


def rectangular_rank(matrix: MatrixZ) -> int:
    if not matrix:
        return 0
    return int(Matrix(matrix).rank())


def permutation_sign(indices: list[int]) -> int:
    inversions = sum(
        indices[left] > indices[right]
        for left in range(len(indices))
        for right in range(left + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def signed_row_equivalence(
    source: MatrixZ, target: MatrixZ
) -> tuple[list[int], list[int]]:
    by_row: dict[tuple[int, ...], deque[int]] = defaultdict(deque)
    for index, row in enumerate(target):
        by_row[tuple(row)].append(index)
    target_indices: list[int] = []
    signs: list[int] = []
    used: set[int] = set()
    for row in source:
        key = tuple(row)
        sign = 1
        candidates = by_row.get(key)
        while candidates and candidates[0] in used:
            candidates.popleft()
        if not candidates:
            key = tuple(-value for value in row)
            sign = -1
            candidates = by_row.get(key)
            while candidates and candidates[0] in used:
                candidates.popleft()
        if not candidates:
            raise AssertionError("matrices are not signed-row equivalent")
        target_index = candidates.popleft()
        used.add(target_index)
        target_indices.append(target_index)
        signs.append(sign)
    if len(used) != len(target):
        raise AssertionError("signed-row match is not bijective")
    return target_indices, signs


def orbit_census(representatives: list[int], L: int) -> dict[str, object]:
    closure = set().union(*(orbit(value, L) for value in representatives))
    fixed = {
        element: sum(
            group_action(config, L, element) == config for config in closure
        )
        for element in GROUP_ELEMENTS
    }
    numerator = sum(fixed.values())
    if numerator % len(GROUP_ELEMENTS):
        raise AssertionError("Burnside numerator is not divisible by group order")
    direct_representatives = sorted(
        {canonical_representative(config, L) for config in closure}
    )
    return {
        "configuration_count": len(closure),
        "fixed_points": fixed,
        "burnside_numerator": numerator,
        "burnside_orbit_count": numerator // len(GROUP_ELEMENTS),
        "direct_orbit_count": len(direct_representatives),
        "direct_representatives": direct_representatives,
        "closure_sha256": canonical_sha256(sorted(closure)),
        "invariant_under_all_actions": all(
            group_action(config, L, element) in closure
            for config in closure
            for element in GROUP_ELEMENTS
        ),
    }


def exhaustive_family_configs(L: int, m: int) -> set[int]:
    expected: set[int] = set()
    for config in range(1 << (2 * L)):
        if config.bit_count() != 2 * m:
            continue
        record = decode_record(config, L, m)
        if record["packing_boundary"] and record["family"] in {
            "(6,2)",
            "(8,0)",
        }:
            expected.add(config)
    return expected


def exhaustive_defect_control(max_L: int) -> dict[str, object]:
    checked = 0
    passed = True
    for L in range(1, max_L + 1):
        for config in range(1 << (2 * L)):
            if config.bit_count() % 2:
                continue
            m = config.bit_count() // 2
            decoded = decode_sets(config, L)
            E = len(decoded["E"])
            D = len(decoded["D"])
            q = E + D
            passed &= E - D == L - 2 * m
            passed &= q >= abs(L - 2 * m)
            passed &= (q == abs(L - 2 * m)) == (min(E, D) == 0)
            checked += 1
    return {
        "checked_even_particle_coordinates_L1_through_L9": checked,
        "all_identities_hold": passed,
    }


def first_nonzero_core_choice(
    rows: list[dict[str, object]], remaining: set[int], used_rows: set[int]
) -> tuple[list[int], list[int], int]:
    columns = sorted(remaining)
    active = [
        index
        for index, row in enumerate(rows)
        if index not in used_rows and set(row["entries"]) & remaining
    ]
    candidate_count = math.comb(len(active), len(columns))
    if candidate_count > 64:
        raise AssertionError(
            f"residual core search would require {candidate_count} candidates"
        )
    for choice in itertools.combinations(active, len(columns)):
        determinant = bareiss_determinant(dense_submatrix(rows, choice, columns))
        if determinant:
            return list(active), list(choice), determinant
    raise AssertionError("no nonzero residual core minor")


def position_of_record(record: dict[str, object]) -> int:
    position = record["defect_position_mod_reflection"]
    if position is None:
        raise AssertionError("position block requires one support defect")
    return int(position)


def build_case(
    L: int,
    m: int,
    annihilator_record: dict[str, object],
    e237_record: dict[str, object],
) -> tuple[dict[str, object], dict[str, bool], MatrixZ, list[int]]:
    orbit_records = configuration_orbits(L, 2 * m)
    high_orbits, rows = weighted_restriction_rows(
        L, m, orbit_records, annihilator_record
    )
    leaf, remaining, used_rows = leaf_certificate(rows, high_orbits)
    active_local, selected_local, first_determinant = first_nonzero_core_choice(
        rows, remaining, used_rows
    )
    columns = sorted(remaining)
    core_representatives = [
        int(high_orbits[column]["representative"]) for column in columns
    ]
    active_source_indices = [
        int(rows[index]["basis_index"]) for index in active_local
    ]
    selected_source_indices = [
        int(rows[index]["basis_index"]) for index in selected_local
    ]
    matrix = dense_submatrix(rows, selected_local, columns)
    determinant = bareiss_determinant(matrix)
    source_core = e237_record["residual_core"]

    coordinate_records = [
        decode_record(representative, L, m)
        for representative in core_representatives
    ]
    family_counts = Counter(
        str(record["family"]) for record in coordinate_records
    )
    family_censuses: dict[str, object] = {}
    complete_closure: set[int] = set()
    for family in ("(6,2)", "(8,0)"):
        representatives = [
            int(record["representative"])
            for record in coordinate_records
            if record["family"] == family
        ]
        census = orbit_census(representatives, L)
        family_censuses[family] = census
        complete_closure.update(
            config for representative in representatives for config in orbit(representative, L)
        )
    expected_family_configs = exhaustive_family_configs(L, m)

    position_censuses: list[dict[str, object]] = []
    if L == 9:
        for family in ("(6,2)", "(8,0)"):
            for position in range(5):
                representatives = [
                    int(record["representative"])
                    for record in coordinate_records
                    if record["family"] == family
                    and position_of_record(record) == position
                ]
                position_censuses.append(
                    {
                        "family": family,
                        "defect_position_mod_reflection": position,
                        **orbit_census(representatives, L),
                    }
                )

    packing_columns = [
        index
        for index, record in enumerate(coordinate_records)
        if record["family"] == "(6,2)"
    ]
    polar_columns = [
        index
        for index, record in enumerate(coordinate_records)
        if record["family"] == "(8,0)"
    ]
    packing_rows = [
        index
        for index, row in enumerate(matrix)
        if all(row[column] == 0 for column in polar_columns)
    ]
    polar_rows = [
        index for index in range(len(matrix)) if index not in set(packing_rows)
    ]
    packing_matrix = [
        [matrix[row][column] for column in packing_columns]
        for row in packing_rows
    ]
    upper_right = [
        [matrix[row][column] for column in polar_columns]
        for row in packing_rows
    ]
    lower_left = [
        [matrix[row][column] for column in packing_columns]
        for row in polar_rows
    ]
    polar_matrix = [
        [matrix[row][column] for column in polar_columns]
        for row in polar_rows
    ]
    row_reordering = packing_rows + polar_rows
    column_reordering = packing_columns + polar_columns
    reordering_sign = permutation_sign(row_reordering) * permutation_sign(
        column_reordering
    )
    packing_record = matrix_record(packing_matrix)
    polar_record = matrix_record(polar_matrix)
    literal_record = matrix_record(matrix, include_sparse_rows=True)

    position_blocks: list[dict[str, object]] = []
    position_factorization_holds = True
    if L == 9:
        column_positions = [position_of_record(record) for record in coordinate_records]
        row_positions: list[int] = []
        for row in matrix:
            positions = {
                column_positions[column]
                for column, value in enumerate(row)
                if value
            }
            if len(positions) != 1:
                position_factorization_holds = False
                row_positions.append(-1)
            else:
                row_positions.append(next(iter(positions)))
        position_row_order: list[int] = []
        position_column_order: list[int] = []
        block_determinants: list[int] = []
        for position in range(5):
            block_rows = [
                index for index, value in enumerate(row_positions) if value == position
            ]
            block_columns = [
                index
                for index, value in enumerate(column_positions)
                if value == position
            ]
            block = [
                [matrix[row][column] for column in block_columns]
                for row in block_rows
            ]
            block_record = matrix_record(block)
            position_blocks.append(
                {
                    "defect_position_mod_reflection": position,
                    "row_indices": block_rows,
                    "column_indices": block_columns,
                    "matrix": block_record,
                }
            )
            position_row_order.extend(block_rows)
            position_column_order.extend(block_columns)
            block_determinants.append(int(block_record["determinant_Z"]))
        position_sign = permutation_sign(position_row_order) * permutation_sign(
            position_column_order
        )
        position_factorization_holds &= (
            sorted(position_row_order) == list(range(len(matrix)))
            and sorted(position_column_order) == list(range(len(matrix)))
            and determinant == position_sign * math.prod(block_determinants)
        )
    else:
        position_sign = None

    determinant_identity = (
        determinant
        == reordering_sign
        * int(packing_record["determinant_Z"])
        * int(polar_record["determinant_Z"])
    )
    expected_family_count = EXPECTED_FAMILY_COUNTS[(L, m)]
    defect_identities_hold = all(
        int(record["particle_count"]) == 2 * m
        and int(record["signed_defect_E_minus_D"]) == L - 2 * m
        and int(record["q_E_plus_D"]) >= abs(L - 2 * m)
        and bool(record["packing_boundary"])
        and bool(record["no_cancelling_empty_double_pair"])
        for record in coordinate_records
    )
    burnside_holds = all(
        census["invariant_under_all_actions"]
        and census["burnside_orbit_count"] == census["direct_orbit_count"]
        and census["direct_representatives"]
        == sorted(
            int(record["representative"])
            for record in coordinate_records
            if record["family"] == family
        )
        for family, census in family_censuses.items()
    )
    if position_censuses:
        burnside_holds &= all(
            census["burnside_orbit_count"] == census["direct_orbit_count"]
            for census in position_censuses
        )

    validations = {
        "leaf_matches_e237": leaf == e237_record["leaf_certificate"]["steps"],
        "remaining_matches_e237": core_representatives
        == source_core["remaining_coordinate_representatives"],
        "active_rows_match_e237": active_source_indices
        == source_core["active_source_basis_indices"],
        "selected_rows_match_e237": selected_source_indices
        == source_core["selected_minor_source_basis_indices"],
        "determinant_matches_e237": determinant
        == int(source_core["determinant_Z"])
        == first_determinant,
        "determinant_matches_named_value": determinant
        == EXPECTED_DETERMINANTS[(L, m)],
        "full_smith_matches_named_value": {
            int(value): int(count)
            for value, count in literal_record["smith_normal_form"][
                "multiplicities"
            ].items()
        }
        == EXPECTED_FULL_SMITH[(L, m)],
        "coordinate_defect_identities": defect_identities_hold,
        "family_counts": dict(family_counts) == expected_family_count,
        "family_coverage_no_extras": complete_closure == expected_family_configs,
        "burnside_equals_direct": burnside_holds,
        "packing_polar_square_blocks": len(packing_rows) == len(packing_columns)
        and len(polar_rows) == len(polar_columns),
        "upper_right_zero": not any(value for row in upper_right for value in row),
        "packing_polar_determinant_identity": determinant_identity,
        "position_factorization": position_factorization_holds,
    }

    certificate = {
        "L": L,
        "m": m,
        "particle_sector": 2 * m,
        "e237_reconstruction": {
            "high_coordinate_count": len(high_orbits),
            "annihilator_row_count": len(rows),
            "leaf_dimension": len(leaf),
            "leaf_steps_sha256": canonical_sha256(leaf),
            "active_source_basis_indices": active_source_indices,
            "selected_minor_source_basis_indices": selected_source_indices,
            "remaining_coordinate_representatives": core_representatives,
        },
        "decoded_coordinates": coordinate_records,
        "support_defect_classification": {
            "family_counts": dict(sorted(family_counts.items())),
            "all_coordinates_on_packing_boundary": defect_identities_hold,
            "family_closure_equals_exhaustive_coordinate_set": complete_closure
            == expected_family_configs,
            "family_orbit_censuses": family_censuses,
            "defect_position_orbit_censuses": position_censuses,
        },
        "literal_core_matrix": literal_record,
        "packing_polar_partition": {
            "packing_family": "(6,2)",
            "polar_family": "(8,0)",
            "packing_row_indices": packing_rows,
            "polar_row_indices": polar_rows,
            "packing_column_indices": packing_columns,
            "polar_column_indices": polar_columns,
            "row_permutation_sign": permutation_sign(row_reordering),
            "column_permutation_sign": permutation_sign(column_reordering),
            "combined_reordering_sign": reordering_sign,
            "upper_right_zero": not any(
                value for row in upper_right for value in row
            ),
            "packing_matrix": packing_record,
            "lower_left_coupling": {
                "shape": [
                    len(lower_left),
                    len(lower_left[0]) if lower_left else 0,
                ],
                "nnz": sum(value != 0 for row in lower_left for value in row),
                "rank_exact_Q": rectangular_rank(lower_left),
                "sha256": matrix_sha256(lower_left),
            },
            "polar_matrix": polar_record,
            "literal_determinant_from_blocks": reordering_sign
            * int(packing_record["determinant_Z"])
            * int(polar_record["determinant_Z"]),
        },
        "defect_position_blocks": {
            "applies": L == 9,
            "combined_reordering_sign": position_sign,
            "blocks": position_blocks,
            "exact_direct_sum_after_row_column_permutations": position_factorization_holds,
        },
    }
    return certificate, validations, matrix, core_representatives


def delete_rung(config: int, rung: int) -> int:
    lower_mask = (1 << (2 * rung)) - 1
    lower = config & lower_mask
    upper = config >> (2 * (rung + 1))
    return lower | (upper << (2 * rung))


def central_inheritance(
    case8: dict[str, object],
    case9: dict[str, object],
    matrix8: MatrixZ,
    matrix9: MatrixZ,
    m: int,
) -> dict[str, object]:
    block = next(
        record
        for record in case9["defect_position_blocks"]["blocks"]
        if record["defect_position_mod_reflection"] == 4
    )
    row_indices = [int(value) for value in block["row_indices"]]
    column_indices = [int(value) for value in block["column_indices"]]
    representatives8 = [
        int(value)
        for value in case8["e237_reconstruction"][
            "remaining_coordinate_representatives"
        ]
    ]
    representatives9 = [
        int(value)
        for value in case9["e237_reconstruction"][
            "remaining_coordinate_representatives"
        ]
    ]
    index8 = {representative: index for index, representative in enumerate(representatives8)}
    target_columns: list[int] = []
    for column in column_indices:
        reduced = delete_rung(representatives9[column], 4)
        if m == 5:
            reduced = complement_config(reduced, 8)
        target_columns.append(index8[canonical_representative(reduced, 8)])
    column_order = sorted(range(len(target_columns)), key=target_columns.__getitem__)
    central = [
        [matrix9[row][column_indices[column]] for column in column_order]
        for row in row_indices
    ]
    target_rows, row_signs = signed_row_equivalence(central, matrix8)
    identity_holds = all(
        central[row]
        == [row_signs[row] * value for value in matrix8[target_rows[row]]]
        for row in range(len(central))
    )
    return {
        "source_case": [9, m],
        "central_defect_position": 4,
        "operation_on_coordinates": (
            "delete the central empty rung"
            if m == 4
            else "delete the central double rung, then particle-hole complement"
        ),
        "central_row_indices_in_L9_core": row_indices,
        "central_column_indices_in_L9_core": column_indices,
        "column_target_indices_in_L8_core": target_columns,
        "row_target_indices_in_L8_core": target_rows,
        "row_signs": row_signs,
        "exact_signed_permutation_identity": identity_holds,
        "central_matrix_sha256": matrix_sha256(central),
        "L8_matrix_sha256": matrix_sha256(matrix8),
    }


def particle_hole_certificate(
    case4: dict[str, object],
    case5: dict[str, object],
    matrix4: MatrixZ,
    matrix5: MatrixZ,
) -> dict[str, object]:
    representatives4 = [
        int(value)
        for value in case4["e237_reconstruction"][
            "remaining_coordinate_representatives"
        ]
    ]
    representatives5 = [
        int(value)
        for value in case5["e237_reconstruction"][
            "remaining_coordinate_representatives"
        ]
    ]
    index5 = {representative: index for index, representative in enumerate(representatives5)}
    column_permutation = [
        index5[
            canonical_representative(complement_config(representative, 9), 9)
        ]
        for representative in representatives4
    ]
    coordinate_bijection = sorted(column_permutation) == list(range(len(matrix4)))
    commutes = all(
        complement_config(group_action(config, 9, element), 9)
        == group_action(complement_config(config, 9), 9, element)
        for config in representatives4
        for element in GROUP_ELEMENTS
    )
    matrix_identity = all(
        matrix5[row][column_permutation[column]] == matrix4[row][column]
        for row in range(len(matrix4))
        for column in range(len(matrix4))
    )
    decoded_swap = True
    for index, representative4 in enumerate(representatives4):
        record4 = decode_record(representative4, 9, 4)
        record5 = decode_record(
            representatives5[column_permutation[index]], 9, 5
        )
        counts4 = record4["counts"]
        counts5 = record5["counts"]
        decoded_swap &= (
            counts5["E"] == counts4["D"]
            and counts5["D"] == counts4["E"]
            and sorted((counts5["T"], counts5["B"]))
            == sorted((counts4["T"], counts4["B"]))
            and record5["defect_position_mod_reflection"]
            == record4["defect_position_mod_reflection"]
        )
    return {
        "map": "bitwise complement on all 18 sites, followed by the canonical <tau,rho> representative",
        "m4_to_m5_column_permutation": column_permutation,
        "coordinate_bijection": coordinate_bijection,
        "commutes_with_tau_and_rho": commutes,
        "empty_double_counts_swap_and_defect_position_is_preserved": decoded_swap,
        "selected_core_rows_remain_in_the_same_order": True,
        "exact_matrix_identity_after_column_permutation": matrix_identity,
    }


def normalized_shell_record(block: MatrixZ) -> dict[str, object]:
    if any(value % 4 for row in block for value in row):
        raise AssertionError("off-centre shell is not entrywise divisible by four")
    normalized = [[value // 4 for value in row] for row in block]
    return matrix_record(normalized)


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    source_paths = (
        Path(__file__).resolve(),
        VERIFIER,
        PROOF,
        E230_PRODUCER,
        E237_PRODUCER,
        E230_ARTIFACT,
        E237_ARTIFACT,
        L8_SOURCE,
        W9_METADATA,
        W9_BASIS,
        PYPROJECT,
        UV_LOCK,
    )
    source_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
    }
    expected_source_paths = {
        "experiments/e249_wlaw_support_defect.py",
        "tests/test_wlaw_support_defect.py",
        "proofs/wlaw_support_defect.md",
        *EXPECTED_UPSTREAM_SHA256,
    }
    add_check(
        checks,
        "current_source_hashes_bound",
        set(source_hashes) == expected_source_paths
        and all(len(digest) == 64 for digest in source_hashes.values()),
        "producer, clean-room verifier, proof, both predecessors, all exact source artifacts, and locked dependency manifests",
    )
    add_check(
        checks,
        "immutable_upstream_sha256_gate",
        all(source_hashes[path] == digest for path, digest in EXPECTED_UPSTREAM_SHA256.items()),
        "every e230/e237 producer and exact artifact equals its Wave-25 frozen SHA-256",
    )
    add_check(
        checks,
        "locked_runtime_versions",
        platform.python_version() == EXPECTED_PYTHON_VERSION
        and sp.__version__ == EXPECTED_SYMPY_VERSION,
        (
            f"Python {platform.python_version()} and SymPy {sp.__version__}; "
            f"expected {EXPECTED_PYTHON_VERSION}/{EXPECTED_SYMPY_VERSION}"
        ),
    )

    e230 = json.loads(E230_ARTIFACT.read_text())
    e237 = json.loads(E237_ARTIFACT.read_text())
    l8 = json.loads(L8_SOURCE.read_text())
    w9_metadata = json.loads(W9_METADATA.read_text())
    w9_basis_envelope = json.loads(W9_BASIS.read_text())
    w9_basis = w9_basis_envelope["basis"]
    l8_record = l8["data"]["L8_record"]
    w9_source = w9_metadata["data"]["W9"]
    w9_record = {
        "L": 9,
        "K_dim": w9_metadata["data"]["K9_dimension_record"]["K_dim_formula"],
        "W_dim": w9_source["value"],
        "W_sectors": w9_source["W_sectors"],
        "cyclic_Q_dim": w9_source["cyclic_Q_dim"],
        "basis": w9_basis,
        "validation": w9_source["validation"],
    }
    records_by_L = {8: l8_record, 9: w9_record}

    e230_row = next(
        row
        for row in e230["data"]["two_slice_observability_exact_L3_8"]
        if (int(row["L"]), int(row["m"])) == (8, 4)
    )
    e237_records = {
        (int(row["L"]), int(row["m"])): row
        for row in e237["data"]["two_slice_observability_exact_L3_9"]
    }
    e237_nonempty = {
        key
        for key, row in e237_records.items()
        if int(row["residual_core"]["dimension"]) > 0
    }
    predecessor_fields = (
        "K_sector_dimension",
        "retained_abs_leg_imbalance_0_or_2_dimension",
        "omitted_abs_leg_imbalance_at_least_4_dimension",
        "annihilator_W_sector_dimension_exact_Q",
        "cyclic_U_sector_dimension_exact_Q",
        "retained_restriction_rank_exact_Q",
    )
    predecessor_semantics = (
        all(check["passed"] is True for check in e230["checks"])
        and all(check["passed"] is True for check in e237["checks"])
        and e230["meta"]["source_artifacts"]
        == {"results/ladder/l8_saturation.json": source_hashes["results/ladder/l8_saturation.json"]}
        and e237["meta"]["source_artifacts"]
        == {
            "results/ladder/l8_saturation.json": source_hashes["results/ladder/l8_saturation.json"],
            "results/ladder/w9_basis.json": source_hashes["results/ladder/w9_basis.json"],
            "results/ladder/w9_saturation.json": source_hashes["results/ladder/w9_saturation.json"],
        }
        and e237_nonempty == set(TARGETS)
        and all(
            e230_row[field] == e237_records[(8, 4)][field]
            for field in predecessor_fields
        )
    )
    add_check(
        checks,
        "predecessor_semantic_and_provenance_gate",
        predecessor_semantics,
        "e230 and e237 pass internally, bind the current basis artifacts, agree at (8,4), and name exactly three nonempty cores",
    )
    exact_sources = (
        source_exact(l8_record)
        and source_exact(w9_record)
        and canonical_basis_sha256(w9_basis)
        == w9_source["basis_sha256"]
        == w9_basis_envelope["sha256"]
        == e237["meta"]["canonical_w9_basis_sha256"]
    )
    add_check(
        checks,
        "stored_annihilator_bases_exact_Q",
        exact_sources,
        "the literal L8 and L9 bases satisfy the inherited exact-Q annihilator predicates and the W9 canonical digest",
    )

    cases: dict[str, dict[str, object]] = {}
    validations: dict[tuple[int, int], dict[str, bool]] = {}
    matrices: dict[tuple[int, int], MatrixZ] = {}
    representatives: dict[tuple[int, int], list[int]] = {}
    for L, m in TARGETS:
        certificate, validation, matrix, reps = build_case(
            L, m, records_by_L[L], e237_records[(L, m)]
        )
        cases[f"L{L}_m{m}"] = certificate
        validations[(L, m)] = validation
        matrices[(L, m)] = matrix
        representatives[(L, m)] = reps

    reconstruction_keys = {
        "leaf_matches_e237",
        "remaining_matches_e237",
        "active_rows_match_e237",
        "selected_rows_match_e237",
        "determinant_matches_e237",
        "determinant_matches_named_value",
        "full_smith_matches_named_value",
    }
    add_check(
        checks,
        "literal_e237_cores_and_determinants_rebuilt",
        all(
            all(validation[name] for name in reconstruction_keys)
            for validation in validations.values()
        ),
        "deterministic leaf replay gives the exact e237 coordinates, active/selected rows, literal matrices, determinants, and full Smith forms",
    )

    defect_control = exhaustive_defect_control(9)
    add_check(
        checks,
        "elementary_defect_identity_exhaustive_control_L1_9",
        defect_control["all_identities_hold"] is True,
        f"checked {defect_control['checked_even_particle_coordinates_L1_through_L9']} even-particle coordinates; the proof is algebraic for every L",
    )
    add_check(
        checks,
        "finite_support_defect_classification_no_extras",
        all(
            validation["coordinate_defect_identities"]
            and validation["family_counts"]
            and validation["family_coverage_no_extras"]
            for validation in validations.values()
        ),
        "the 17 and 133 residual coordinate sets are exactly the named packing-boundary (6,2) and (8,0) orbit families",
    )
    add_check(
        checks,
        "explicit_group_actions_burnside_equals_direct",
        all(validation["burnside_equals_direct"] for validation in validations.values()),
        "fixed points are enumerated under id, tau, rho, and tau*rho on each explicit orbit closure; no Burnside count is inserted as a constant",
    )

    particle_hole = particle_hole_certificate(
        cases["L9_m4"], cases["L9_m5"], matrices[(9, 4)], matrices[(9, 5)]
    )
    add_check(
        checks,
        "particle_hole_bijection_and_literal_matrix_identity",
        all(
            particle_hole[name]
            for name in (
                "coordinate_bijection",
                "commutes_with_tau_and_rho",
                "empty_double_counts_swap_and_defect_position_is_preserved",
                "exact_matrix_identity_after_column_permutation",
            )
        ),
        "bitwise complement maps the complete m=4 core to m=5 and the selected literal matrices agree after its column permutation",
    )

    add_check(
        checks,
        "packing_polar_block_triangularity_and_exact_ranks",
        all(
            validation["packing_polar_square_blocks"]
            and validation["upper_right_zero"]
            and validation["packing_polar_determinant_identity"]
            for validation in validations.values()
        ),
        "each literal core is [A 0; C D] after recorded permutations, with full-rank exact Smith blocks",
    )
    add_check(
        checks,
        "L9_defect_position_direct_sum",
        validations[(9, 4)]["position_factorization"]
        and validations[(9, 5)]["position_factorization"],
        "each L9 row and column has one reflection-class defect position, giving four 29x29 shells and one central 17x17 block",
    )

    central4 = central_inheritance(
        cases["L8_m4"], cases["L9_m4"], matrices[(8, 4)], matrices[(9, 4)], 4
    )
    central5 = central_inheritance(
        cases["L8_m4"], cases["L9_m5"], matrices[(8, 4)], matrices[(9, 5)], 5
    )
    add_check(
        checks,
        "central_L8_core_inherited_by_signed_permutations",
        central4["exact_signed_permutation_identity"]
        and central5["exact_signed_permutation_identity"],
        "deleting the central defect identifies each L9 central 17x17 block with the L8 core up to one recorded signed row permutation",
    )

    shell_records: dict[str, list[dict[str, object]]] = {}
    shell_gate = True
    for m in (4, 5):
        case = cases[f"L9_m{m}"]
        records: list[dict[str, object]] = []
        for block in case["defect_position_blocks"]["blocks"][:4]:
            rows = [int(value) for value in block["row_indices"]]
            columns = [int(value) for value in block["column_indices"]]
            matrix = [
                [matrices[(9, m)][row][column] for column in columns]
                for row in rows
            ]
            normalized = normalized_shell_record(matrix)
            records.append(
                {
                    "defect_position_mod_reflection": block[
                        "defect_position_mod_reflection"
                    ],
                    "divisor": 4,
                    "normalized_matrix": normalized,
                }
            )
            shell_gate &= normalized["smith_normal_form"]["multiplicities"] == {
                "1": 23,
                "2": 5,
                "6": 1,
            }
        shell_records[f"L9_m{m}"] = records

    packing8 = cases["L8_m4"]["packing_polar_partition"]["packing_matrix"]
    packing9 = cases["L9_m4"]["packing_polar_partition"]["packing_matrix"]
    eight_copy_determinant = int(packing8["determinant_Z"]) ** 8
    eight_copy_obstruction = {
        "proposed_recurrence": "A_9 is unimodularly row-column equivalent to the direct sum of eight copies of A_8 because 128=8*16",
        "A8_dimension": 16,
        "A9_dimension": 128,
        "eight_copy_dimension_matches": 8 * 16 == 128,
        "absolute_determinant_of_eight_A8_copies": abs(eight_copy_determinant),
        "eight_copy_determinant_factorization": integer_factorization(
            eight_copy_determinant
        ),
        "actual_A9_absolute_determinant": abs(int(packing9["determinant_Z"])),
        "actual_A9_determinant_factorization": packing9[
            "determinant_factorization"
        ],
        "determinants_disagree": abs(eight_copy_determinant)
        != abs(int(packing9["determinant_Z"])),
        "first_shell_failure": {
            "defect_position_mod_reflection": 0,
            "normalized_smith_multiplicities": shell_records["L9_m4"][0][
                "normalized_matrix"
            ]["smith_normal_form"]["multiplicities"],
            "reason": "the normalized shell already has invariant 6 and determinant divisible by 3, so it is not a dyadic/unimodular shell",
        },
        "scope": "this refutes only the named count-driven eight-copy/dyadic-shell recurrence, not every possible block recurrence",
    }
    add_check(
        checks,
        "first_exact_eight_copy_recurrence_obstruction",
        shell_gate
        and eight_copy_obstruction["eight_copy_dimension_matches"]
        and eight_copy_obstruction["determinants_disagree"]
        and eight_copy_obstruction["first_shell_failure"][
            "normalized_smith_multiplicities"
        ]
        == {"1": 23, "2": 5, "6": 1},
        "the tempting 128=8*16 direct-sum recurrence fails by determinant, and its first normalized shell has Smith invariant 6",
    )

    determinant_mechanism = {
        "L8": {
            "formula": "|det M_8|=|det A_8|*|det D_8|=(3*2^31)*2=3*2^32",
            "full_smith_multiplicities": cases["L8_m4"]["literal_core_matrix"][
                "smith_normal_form"
            ]["multiplicities"],
        },
        "L9": {
            "formula": "|det M_9|=(3*2^64)^4*(3*2^32)=3^5*2^288",
            "four_off_centre_shell_smith_multiplicities": cases["L9_m4"][
                "defect_position_blocks"
            ]["blocks"][0]["matrix"]["smith_normal_form"]["multiplicities"],
            "inherited_central_smith_multiplicities": cases["L9_m4"][
                "defect_position_blocks"
            ]["blocks"][4]["matrix"]["smith_normal_form"]["multiplicities"],
            "full_smith_multiplicities": cases["L9_m4"]["literal_core_matrix"][
                "smith_normal_form"
            ]["multiplicities"],
        },
        "all_claimed_matrices_have_full_exact_rank": all(
            case["literal_core_matrix"]["rank_exact_Q"]
            == case["literal_core_matrix"]["shape"][1]
            and case["packing_polar_partition"]["packing_matrix"][
                "rank_exact_Q"
            ]
            == case["packing_polar_partition"]["packing_matrix"]["shape"][1]
            and case["packing_polar_partition"]["polar_matrix"]["rank_exact_Q"]
            == case["packing_polar_partition"]["polar_matrix"]["shape"][1]
            for case in cases.values()
        ),
    }
    add_check(
        checks,
        "smith_blocks_explain_all_powers_of_two_and_three",
        determinant_mechanism["all_claimed_matrices_have_full_exact_rank"]
        and cases["L8_m4"]["literal_core_matrix"]["determinant_factorization"]
        == {"sign": -1, "exponent_2": 32, "exponent_3": 1, "residual_after_2_and_3": 1}
        and cases["L9_m4"]["literal_core_matrix"]["determinant_factorization"]
        == {"sign": -1, "exponent_2": 288, "exponent_3": 5, "residual_after_2_and_3": 1}
        and cases["L9_m5"]["literal_core_matrix"]["determinant_factorization"]
        == {"sign": -1, "exponent_2": 288, "exponent_3": 5, "residual_after_2_and_3": 1},
        "exact Smith diagonals and permutation signs reproduce both literal determinants without modular promotion",
    )

    elementary_theorem = {
        "tag": "[LEMMA]",
        "assumptions": "L>=1, 0<=m<=L, and a rung word with exactly 2m occupied sites",
        "partition": "{0,...,L-1}=T disjoint_union B disjoint_union E disjoint_union D",
        "particle_identity": "2m=|T|+|B|+2|D|",
        "rung_identity": "L=|T|+|B|+|E|+|D|",
        "defect_identity": "|E|-|D|=L-2m",
        "inequality": "q=|E|+|D|>=|L-2m|",
        "equality_criterion": "q=|L-2m| iff min(|E|,|D|)=0; equivalently there is no cancelling empty-double pair and every rung beyond the forced defect set is singly occupied",
    }
    finite_theorem = {
        "tag": "[THEOREM][COMPUTATION]",
        "scope": "the deterministic e237 raw-leaf residual minors for exactly (8,4),(9,4),(9,5)",
        "classification": {
            "L8_m4": "16 (6,2) no-empty/no-double orbits plus one (8,0) orbit",
            "L9_m4": "128 (6,2) one-empty/no-double orbits plus five (8,0) one-empty orbits",
            "L9_m5": "the particle-hole mirror: 128 (6,2) one-double/no-empty orbits plus five (8,0) one-double orbits",
        },
        "matrix_statement": "the literal cores have the recorded exact Smith forms; L9 is a defect-position direct sum of four 29x29 shells and one signed-permutation copy of the L8 17x17 core",
        "determinants": {
            "L8_m4": "-3*2^32",
            "L9_m4": "-3^5*2^288",
            "L9_m5": "-3^5*2^288",
        },
    }
    scope = {
        "tag": "[SCOPE]",
        "claim": SCOPE_CLAIM,
        "proved": [
            "the elementary support-defect identity and inequality for every even-particle rung coordinate",
            "the complete finite residual classifications 17=16+1 and 133=128+5 for the three stored e237 cores",
            "the explicit Burnside/direct-orbit agreement and the L9 particle-hole bijection",
            "the exact packing/polar and defect-position Smith-block determinant mechanism",
            "the determinant obstruction to the named eight-copy/dyadic-shell recurrence",
        ],
        "not_proved": [
            "that the e237 residual core equals a support-defect boundary for any L>=10",
            "a uniform recurrence, contracting homotopy, or all-L two-slice rank theorem",
            "a no-go in an arbitrary basis or a Pfaffian-number lower bound",
            "an exact solution or thermodynamic statement for the three-dimensional Ising model",
        ],
        "L_greater_equal_10_status": "[UNRESOLVED] no exact source is loaded and every residual or rank statement remains open",
        "sequence_fit_used": False,
    }
    add_check(
        checks,
        "scope_and_no_recurrence_fit_gate",
        scope["claim"] == SCOPE_CLAIM
        and scope["sequence_fit_used"] is False
        and "[UNRESOLVED]" in scope["L_greater_equal_10_status"],
        SCOPE_CLAIM,
    )
    add_check(
        checks,
        "no_L10_closure_no_benchmark",
        True,
        "only stored L8/L9 bases are read; no L>=10 closure, physical coupling, floating point, or determinant-sequence fit is used",
    )

    certificate = {
        "group_action": {
            "group": "K=<tau,rho> is the Klein four group",
            "elements": list(GROUP_ELEMENTS),
            "tau": "swap top and bottom on every rung: T<->B, E->E, D->D",
            "rho": "reverse rung order r->L-1-r",
            "tau_rho": "apply rho then tau (the generators commute)",
            "burnside_formula": "number of K-orbits=(Fix(id)+Fix(tau)+Fix(rho)+Fix(tau_rho))/4",
        },
        "finite_decoder_control": defect_control,
        "cases": cases,
        "particle_hole_L9": particle_hole,
        "central_L8_inheritance": [central4, central5],
        "normalized_off_centre_shells": shell_records,
        "determinant_mechanism": determinant_mechanism,
        "eight_copy_recurrence_obstruction": eight_copy_obstruction,
    }
    data_payload = {
        "claim_tag": "[LEMMA][THEOREM][FINITE EXACT COMPUTATION][UNRESOLVED L>=10]",
        "elementary_defect_theorem": elementary_theorem,
        "finite_residual_theorem": finite_theorem,
        "certificate": certificate,
        "scope": scope,
    }

    elapsed = time.process_time() - started
    peak_rss = peak_rss_bytes()
    rss_measurement = peak_rss_measurement()
    add_check(
        checks,
        "declared_resource_limits",
        elapsed < CPU_BUDGET_SECONDS and peak_rss < RSS_LIMIT_BYTES,
        f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes via {rss_measurement}",
    )
    names = [str(check["name"]) for check in checks]
    if len(names) != len(set(names)):
        raise AssertionError("producer check names are not unique")
    failed = [str(check["name"]) for check in checks if not check["passed"]]
    if failed:
        raise AssertionError(f"W-law support-defect checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e249_wlaw_support_defect.py",
            "verifier": "tests/test_wlaw_support_defect.py",
            "proof": "proofs/wlaw_support_defect.md",
            "interpreter": sys.executable,
            "arithmetic": "Python integers and exact Smith normal forms over ZZ; no floating point or modular rank promotion",
            "runtime_versions": {
                "python": platform.python_version(),
                "sympy": sp.__version__,
            },
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "peak_rss_measurement": rss_measurement,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "L10_closure_run": False,
            "source_sha256": source_hashes,
            "frozen_upstream_sha256": EXPECTED_UPSTREAM_SHA256,
            "certificate_sha256": canonical_sha256(certificate),
            "data_sha256": canonical_sha256(data_payload),
        },
        "data": data_payload,
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT)
    for check in payload["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"  [{status}] {check['name']}: {check['detail']}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
