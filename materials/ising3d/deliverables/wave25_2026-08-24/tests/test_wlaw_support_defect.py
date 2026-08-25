#!/usr/bin/env python3
"""Clean-room verifier for the e249 W-law support-defect certificate.

This file does not import e249 or e237.  It independently enumerates the
K=<tau,rho> coordinate orbits, replays the raw singleton leaf rule against the
stored exact W bases, rebuilds all three literal integer cores, and checks the
decisive determinant and Smith-block identities.
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
from collections import Counter
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw_support_defect.json"
PRODUCER = ROOT / "experiments" / "e249_wlaw_support_defect.py"
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
EXPECTED_FAMILIES = {
    (8, 4): {"(6,2)": 16, "(8,0)": 1},
    (9, 4): {"(6,2)": 128, "(8,0)": 5},
    (9, 5): {"(6,2)": 128, "(8,0)": 5},
}
EXPECTED_FULL_SMITH = {
    (8, 4): {"2": 4, "4": 11, "8": 1, "24": 1},
    (9, 4): {"2": 4, "4": 103, "8": 21, "24": 5},
    (9, 5): {"2": 4, "4": 103, "8": 21, "24": 5},
}
EXPECTED_SCOPE = (
    "the elementary rung-defect inequality for every even-particle coordinate, "
    "and the exact support, Burnside, particle-hole, determinant, and Smith-block "
    "classification of only the stored e237 cores (8,4), (9,4), and (9,5); "
    "no recurrence or rank statement for L>=10"
)
EXPECTED_PYTHON_VERSION = "3.14.3"
EXPECTED_SYMPY_VERSION = "1.14.0"
EXPECTED_CLAIM_TAG = "[LEMMA][THEOREM][FINITE EXACT COMPUTATION][UNRESOLVED L>=10]"
FAILURES: list[str] = []

MatrixZ = list[list[int]]


def check(name: str, passed: bool, detail: str = "") -> None:
    suffix = f": {detail}" if detail else ""
    print(f"{'PASS' if passed else 'FAIL'} {name}{suffix}")
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


def tau(config: int, L: int) -> int:
    result = 0
    for rung in range(L):
        result |= ((config >> (2 * rung)) & 1) << (2 * rung + 1)
        result |= ((config >> (2 * rung + 1)) & 1) << (2 * rung)
    return result


def rho(config: int, L: int) -> int:
    result = 0
    for rung in range(L):
        result |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return result


def act(config: int, L: int, element: str) -> int:
    if element == "id":
        return config
    if element == "tau":
        return tau(config, L)
    if element == "rho":
        return rho(config, L)
    if element == "tau_rho":
        return tau(rho(config, L), L)
    raise ValueError(element)


def orbit(config: int, L: int) -> set[int]:
    return {act(config, L, element) for element in GROUP_ELEMENTS}


def canonical(config: int, L: int) -> int:
    return min(orbit(config, L))


def complement(config: int, L: int) -> int:
    return config ^ ((1 << (2 * L)) - 1)


def decode(config: int, L: int) -> dict[str, list[int]]:
    result = {"T": [], "B": [], "E": [], "D": []}
    labels = ("E", "T", "B", "D")
    for rung in range(L):
        result[labels[(config >> (2 * rung)) & 3]].append(rung)
    return result


def coordinate_record(config: int, L: int, m: int) -> dict[str, object]:
    sets = decode(config, L)
    counts = {name: len(values) for name, values in sets.items()}
    q = counts["E"] + counts["D"]
    singles = sorted((counts["T"], counts["B"]), reverse=True)
    defects = sets["E"] + sets["D"]
    position = (
        min(defects[0], L - 1 - defects[0]) if len(defects) == 1 else None
    )
    return {
        "representative": config,
        "T": sets["T"],
        "B": sets["B"],
        "E": sets["E"],
        "D": sets["D"],
        "counts": counts,
        "particle_count": counts["T"] + counts["B"] + 2 * counts["D"],
        "absolute_leg_imbalance": abs(counts["T"] - counts["B"]),
        "signed_defect_E_minus_D": counts["E"] - counts["D"],
        "q_E_plus_D": q,
        "packing_boundary": q == abs(L - 2 * m),
        "no_cancelling_empty_double_pair": min(counts["E"], counts["D"]) == 0,
        "family": f"({singles[0]},{singles[1]})",
        "defect_position_mod_reflection": position,
    }


def config_orbits(L: int, sector: int) -> list[dict[str, int]]:
    seen: set[int] = set()
    records: list[dict[str, int]] = []
    for config in range(1 << (2 * L)):
        if config.bit_count() != sector or config in seen:
            continue
        values = orbit(config, L)
        seen.update(values)
        representative = min(values)
        sets = decode(representative, L)
        records.append(
            {
                "representative": representative,
                "orbit_size": len(values),
                "absolute_leg_imbalance": abs(len(sets["T"]) - len(sets["B"])),
            }
        )
    return records


def weighted_rows(
    L: int,
    m: int,
    orbit_records: list[dict[str, int]],
    basis: list[list[dict[str, int]]],
) -> tuple[list[dict[str, int]], list[dict[str, object]]]:
    high = [
        record
        for record in orbit_records
        if record["absolute_leg_imbalance"] >= 4
    ]
    index = {record["representative"]: j for j, record in enumerate(high)}
    sizes = {record["representative"]: record["orbit_size"] for record in high}
    rows: list[dict[str, object]] = []
    for basis_index, vector in enumerate(basis):
        if not vector or int(vector[0]["representative"]).bit_count() != 2 * m:
            continue
        entries: dict[int, int] = {}
        for term in vector:
            representative = int(term["representative"])
            column = index.get(representative)
            if column is not None:
                entries[column] = int(term["coefficient"]) * sizes[representative]
        rows.append({"basis_index": basis_index, "entries": entries})
    return high, rows


def peel(
    rows: list[dict[str, object]], high: list[dict[str, int]]
) -> tuple[list[dict[str, int]], set[int], set[int]]:
    supports = [set(row["entries"]) for row in rows]
    neighbours: list[list[int]] = [[] for _ in high]
    for row_index, support in enumerate(supports):
        for column in support:
            neighbours[column].append(row_index)
    remaining = set(range(len(high)))
    used: set[int] = set()
    queue = [index for index, support in enumerate(supports) if len(support) == 1]
    heapq.heapify(queue)
    steps: list[dict[str, int]] = []
    while queue:
        row_index = heapq.heappop(queue)
        if row_index in used:
            continue
        active = supports[row_index] & remaining
        if len(active) != 1:
            continue
        pivot = next(iter(active))
        entries = rows[row_index]["entries"]
        steps.append(
            {
                "source_basis_index": int(rows[row_index]["basis_index"]),
                "pivot_representative": int(high[pivot]["representative"]),
                "diagonal_entry": int(entries[pivot]),
            }
        )
        used.add(row_index)
        remaining.remove(pivot)
        for adjacent in neighbours[pivot]:
            if adjacent not in used and len(supports[adjacent] & remaining) == 1:
                heapq.heappush(queue, adjacent)
    return steps, remaining, used


def dense(
    rows: list[dict[str, object]], row_indices: Iterable[int], columns: list[int]
) -> MatrixZ:
    positions = {column: index for index, column in enumerate(columns)}
    matrix: MatrixZ = []
    for row_index in row_indices:
        target = [0] * len(columns)
        for column, value in rows[row_index]["entries"].items():
            if column in positions:
                target[positions[column]] = int(value)
        matrix.append(target)
    return matrix


def determinant(matrix: MatrixZ) -> int:
    n = len(matrix)
    if n == 0:
        return 1
    if any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(n - 1):
        pivot_row = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, n):
            entry = work[row][column]
            for target in range(column + 1, n):
                numerator = work[row][target] * pivot - entry * work[column][target]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise AssertionError("Bareiss division failed")
                work[row][target] = quotient
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def smith(matrix: MatrixZ) -> dict[str, object]:
    form = smith_normal_form(Matrix(matrix), domain=ZZ)
    diagonal = [
        abs(int(form[index, index]))
        for index in range(min(form.rows, form.cols))
        if form[index, index]
    ]
    counts = Counter(diagonal)
    return {
        "rank": len(diagonal),
        "diagonal": diagonal,
        "multiplicities": {str(value): counts[value] for value in sorted(counts)},
        "determinant_abs": math.prod(diagonal),
    }


def sparse_rows(matrix: MatrixZ) -> list[list[list[int]]]:
    return [
        [[column, value] for column, value in enumerate(row) if value]
        for row in matrix
    ]


def matrix_hash(matrix: MatrixZ) -> str:
    return canonical_sha256(
        {
            "shape": [len(matrix), len(matrix[0]) if matrix else 0],
            "sparse_rows": sparse_rows(matrix),
        }
    )


def first_core(
    rows: list[dict[str, object]], remaining: set[int], used: set[int]
) -> tuple[list[int], list[int], MatrixZ, int]:
    columns = sorted(remaining)
    active = [
        index
        for index, row in enumerate(rows)
        if index not in used and set(row["entries"]) & remaining
    ]
    if math.comb(len(active), len(columns)) > 64:
        raise AssertionError("clean-room core search exceeded frozen bound")
    for choice in itertools.combinations(active, len(columns)):
        matrix = dense(rows, choice, columns)
        value = determinant(matrix)
        if value:
            return active, list(choice), matrix, value
    raise AssertionError("clean-room core is singular")


def census(representatives: list[int], L: int) -> dict[str, object]:
    closure = set().union(*(orbit(value, L) for value in representatives))
    fixed = {
        element: sum(act(config, L, element) == config for config in closure)
        for element in GROUP_ELEMENTS
    }
    direct = sorted({canonical(config, L) for config in closure})
    numerator = sum(fixed.values())
    return {
        "configuration_count": len(closure),
        "fixed_points": fixed,
        "burnside_numerator": numerator,
        "burnside_orbit_count": numerator // 4,
        "direct_orbit_count": len(direct),
        "direct_representatives": direct,
        "closure_sha256": canonical_sha256(sorted(closure)),
        "invariant_under_all_actions": all(
            act(config, L, element) in closure
            for config in closure
            for element in GROUP_ELEMENTS
        ),
    }


def expected_family_closure(L: int, m: int) -> set[int]:
    result: set[int] = set()
    for config in range(1 << (2 * L)):
        if config.bit_count() != 2 * m:
            continue
        record = coordinate_record(config, L, m)
        if record["packing_boundary"] and record["family"] in {"(6,2)", "(8,0)"}:
            result.add(config)
    return result


def permutation_sign(indices: list[int]) -> int:
    inversions = sum(
        indices[left] > indices[right]
        for left in range(len(indices))
        for right in range(left + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def delete_rung(config: int, rung: int) -> int:
    lower = config & ((1 << (2 * rung)) - 1)
    upper = config >> (2 * (rung + 1))
    return lower | (upper << (2 * rung))


def source_records() -> tuple[
    dict[int, list[list[dict[str, int]]]], dict[tuple[int, int], dict[str, object]]
]:
    l8 = json.loads(L8_SOURCE.read_text())
    w9_metadata = json.loads(W9_METADATA.read_text())
    w9_basis = json.loads(W9_BASIS.read_text())
    e237 = json.loads(E237_ARTIFACT.read_text())
    expected_w9_hash = w9_metadata["data"]["W9"]["basis_sha256"]
    check(
        "independent W9 canonical basis digest",
        canonical_sha256(w9_basis["basis"])
        == w9_basis["sha256"]
        == expected_w9_hash
        == e237["meta"]["canonical_w9_basis_sha256"],
    )
    bases = {
        8: l8["data"]["L8_record"]["basis"],
        9: w9_basis["basis"],
    }
    records = {
        (int(row["L"]), int(row["m"])): row
        for row in e237["data"]["two_slice_observability_exact_L3_9"]
    }
    return bases, records


def rebuild_case(
    L: int,
    m: int,
    basis: list[list[dict[str, int]]],
    e237_record: dict[str, object],
    stored: dict[str, object],
) -> tuple[MatrixZ, list[int], list[dict[str, object]], bool]:
    orbits = config_orbits(L, 2 * m)
    high, rows = weighted_rows(L, m, orbits, basis)
    steps, remaining, used = peel(rows, high)
    active, selected, matrix, value = first_core(rows, remaining, used)
    columns = sorted(remaining)
    representatives = [int(high[column]["representative"]) for column in columns]
    active_sources = [int(rows[index]["basis_index"]) for index in active]
    selected_sources = [int(rows[index]["basis_index"]) for index in selected]
    source_core = e237_record["residual_core"]
    reconstruction_holds = (
        steps == e237_record["leaf_certificate"]["steps"]
        and canonical_sha256(steps)
        == stored["e237_reconstruction"]["leaf_steps_sha256"]
        and representatives == source_core["remaining_coordinate_representatives"]
        == stored["e237_reconstruction"]["remaining_coordinate_representatives"]
        and active_sources == source_core["active_source_basis_indices"]
        == stored["e237_reconstruction"]["active_source_basis_indices"]
        and selected_sources == source_core["selected_minor_source_basis_indices"]
        == stored["e237_reconstruction"]["selected_minor_source_basis_indices"]
        and value == int(source_core["determinant_Z"])
        == int(stored["literal_core_matrix"]["determinant_Z"])
        == EXPECTED_DETERMINANTS[(L, m)]
        and sparse_rows(matrix) == stored["literal_core_matrix"]["sparse_rows"]
        and matrix_hash(matrix) == stored["literal_core_matrix"]["sha256"]
    )

    coordinates = [coordinate_record(value, L, m) for value in representatives]
    coordinate_holds = coordinates == stored["decoded_coordinates"]
    coordinate_holds &= all(
        int(record["signed_defect_E_minus_D"]) == L - 2 * m
        and int(record["q_E_plus_D"]) >= abs(L - 2 * m)
        and record["packing_boundary"] is True
        and record["no_cancelling_empty_double_pair"] is True
        for record in coordinates
    )
    family_counts = Counter(str(record["family"]) for record in coordinates)
    coordinate_holds &= dict(family_counts) == EXPECTED_FAMILIES[(L, m)]

    closure: set[int] = set()
    family_censuses = {}
    for family in ("(6,2)", "(8,0)"):
        family_reps = [
            int(record["representative"])
            for record in coordinates
            if record["family"] == family
        ]
        rebuilt = census(family_reps, L)
        family_censuses[family] = rebuilt
        coordinate_holds &= rebuilt == stored["support_defect_classification"][
            "family_orbit_censuses"
        ][family]
        coordinate_holds &= rebuilt["burnside_orbit_count"] == rebuilt["direct_orbit_count"]
        closure.update(
            config for representative in family_reps for config in orbit(representative, L)
        )
    coordinate_holds &= closure == expected_family_closure(L, m)
    coordinate_holds &= stored["support_defect_classification"][
        "family_closure_equals_exhaustive_coordinate_set"
    ] is True

    if L == 9:
        rebuilt_position_censuses = []
        for family in ("(6,2)", "(8,0)"):
            for position in range(5):
                position_reps = [
                    int(record["representative"])
                    for record in coordinates
                    if record["family"] == family
                    and record["defect_position_mod_reflection"] == position
                ]
                rebuilt_position_censuses.append(
                    {
                        "family": family,
                        "defect_position_mod_reflection": position,
                        **census(position_reps, L),
                    }
                )
        coordinate_holds &= rebuilt_position_censuses == stored[
            "support_defect_classification"
        ]["defect_position_orbit_censuses"]

    full_smith = smith(matrix)
    matrix_holds = (
        full_smith["rank"] == len(matrix)
        and full_smith["determinant_abs"] == abs(value)
        and full_smith["multiplicities"] == EXPECTED_FULL_SMITH[(L, m)]
        == stored["literal_core_matrix"]["smith_normal_form"]["multiplicities"]
        and full_smith["diagonal"]
        == stored["literal_core_matrix"]["smith_normal_form"]["diagonal"]
    )

    packing_columns = [
        index for index, record in enumerate(coordinates) if record["family"] == "(6,2)"
    ]
    polar_columns = [
        index for index, record in enumerate(coordinates) if record["family"] == "(8,0)"
    ]
    packing_rows = [
        row
        for row in range(len(matrix))
        if all(matrix[row][column] == 0 for column in polar_columns)
    ]
    polar_rows = [row for row in range(len(matrix)) if row not in set(packing_rows)]
    A = [[matrix[row][column] for column in packing_columns] for row in packing_rows]
    B = [[matrix[row][column] for column in polar_columns] for row in packing_rows]
    C = [[matrix[row][column] for column in packing_columns] for row in polar_rows]
    D = [[matrix[row][column] for column in polar_columns] for row in polar_rows]
    partition = stored["packing_polar_partition"]
    block_holds = (
        packing_rows == partition["packing_row_indices"]
        and polar_rows == partition["polar_row_indices"]
        and packing_columns == partition["packing_column_indices"]
        and polar_columns == partition["polar_column_indices"]
        and not any(value for row in B for value in row)
        and determinant(A) == int(partition["packing_matrix"]["determinant_Z"])
        and determinant(D) == int(partition["polar_matrix"]["determinant_Z"])
        and smith(A)["diagonal"]
        == partition["packing_matrix"]["smith_normal_form"]["diagonal"]
        and smith(D)["diagonal"]
        == partition["polar_matrix"]["smith_normal_form"]["diagonal"]
        and matrix_hash(C) == partition["lower_left_coupling"]["sha256"]
        and int(Matrix(C).rank())
        == partition["lower_left_coupling"]["rank_exact_Q"]
        == (1 if L == 8 else 5)
        and value
        == partition["combined_reordering_sign"] * determinant(A) * determinant(D)
    )
    expected_D = [[2]] if L == 8 else [
        [4 if row == column and row < 4 else 2 if row == column else 0 for column in range(5)]
        for row in range(5)
    ]
    block_holds &= D == expected_D

    if L == 9:
        positions = [int(record["defect_position_mod_reflection"]) for record in coordinates]
        row_positions = []
        for row in matrix:
            support_positions = {
                positions[column] for column, value in enumerate(row) if value
            }
            block_holds &= len(support_positions) == 1
            row_positions.append(next(iter(support_positions)))
        row_order: list[int] = []
        column_order: list[int] = []
        block_determinants: list[int] = []
        stored_blocks = stored["defect_position_blocks"]["blocks"]
        for position in range(5):
            block_rows = [index for index, value in enumerate(row_positions) if value == position]
            block_columns = [index for index, value in enumerate(positions) if value == position]
            submatrix = [
                [matrix[row][column] for column in block_columns]
                for row in block_rows
            ]
            stored_block = stored_blocks[position]
            sub_smith = smith(submatrix)
            block_holds &= (
                block_rows == stored_block["row_indices"]
                and block_columns == stored_block["column_indices"]
                and matrix_hash(submatrix) == stored_block["matrix"]["sha256"]
                and determinant(submatrix) == stored_block["matrix"]["determinant_Z"]
                and sub_smith["diagonal"]
                == stored_block["matrix"]["smith_normal_form"]["diagonal"]
                and sub_smith["rank"] == len(submatrix)
            )
            if position < 4:
                normalized = [[entry // 4 for entry in row] for row in submatrix]
                block_holds &= all(entry % 4 == 0 for row in submatrix for entry in row)
                block_holds &= smith(normalized)["multiplicities"] == {
                    "1": 23,
                    "2": 5,
                    "6": 1,
                }
            row_order.extend(block_rows)
            column_order.extend(block_columns)
            block_determinants.append(determinant(submatrix))
        block_holds &= value == (
            permutation_sign(row_order)
            * permutation_sign(column_order)
            * math.prod(block_determinants)
        )

    return matrix, representatives, coordinates, (
        reconstruction_holds and coordinate_holds and matrix_holds and block_holds
    )


def elementary_control() -> tuple[int, bool]:
    checked = 0
    passed = True
    for L in range(1, 10):
        for config in range(1 << (2 * L)):
            if config.bit_count() % 2:
                continue
            m = config.bit_count() // 2
            sets = decode(config, L)
            E = len(sets["E"])
            D = len(sets["D"])
            q = E + D
            passed &= E - D == L - 2 * m
            passed &= q >= abs(L - 2 * m)
            passed &= (q == abs(L - 2 * m)) == (min(E, D) == 0)
            checked += 1
    return checked, passed


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    meta = artifact["meta"]
    data = artifact["data"]
    certificate = data["certificate"]

    source_paths = (
        PRODUCER,
        Path(__file__).resolve(),
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
    current_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
    }
    check(
        "current source hash gate",
        meta["source_sha256"] == current_hashes
        and meta["frozen_upstream_sha256"] == EXPECTED_UPSTREAM_SHA256
        and all(current_hashes[path] == digest for path, digest in EXPECTED_UPSTREAM_SHA256.items()),
        "producer, verifier, proof, predecessors, exact bases, artifacts, pyproject, and uv.lock",
    )
    check(
        "certificate content hash",
        meta["certificate_sha256"] == canonical_sha256(certificate),
    )
    check(
        "complete data content hash",
        meta["data_sha256"] == canonical_sha256(data),
    )
    check(
        "locked runtime versions",
        meta["runtime_versions"]
        == {"python": EXPECTED_PYTHON_VERSION, "sympy": EXPECTED_SYMPY_VERSION}
        and platform.python_version() == EXPECTED_PYTHON_VERSION
        and sp.__version__ == EXPECTED_SYMPY_VERSION,
    )
    required_checks = {
        "current_source_hashes_bound",
        "immutable_upstream_sha256_gate",
        "locked_runtime_versions",
        "predecessor_semantic_and_provenance_gate",
        "stored_annihilator_bases_exact_Q",
        "literal_e237_cores_and_determinants_rebuilt",
        "elementary_defect_identity_exhaustive_control_L1_9",
        "finite_support_defect_classification_no_extras",
        "explicit_group_actions_burnside_equals_direct",
        "particle_hole_bijection_and_literal_matrix_identity",
        "packing_polar_block_triangularity_and_exact_ranks",
        "L9_defect_position_direct_sum",
        "central_L8_core_inherited_by_signed_permutations",
        "first_exact_eight_copy_recurrence_obstruction",
        "smith_blocks_explain_all_powers_of_two_and_three",
        "scope_and_no_recurrence_fit_gate",
        "no_L10_closure_no_benchmark",
        "declared_resource_limits",
    }
    names = [str(record["name"]) for record in artifact["checks"]]
    check(
        "producer checks are hard gates",
        len(names) == len(set(names))
        and set(names) == required_checks
        and all(record["passed"] is True for record in artifact["checks"]),
    )
    check(
        "producer exact arithmetic and resource wall",
        meta["arithmetic"]
        == "Python integers and exact Smith normal forms over ZZ; no floating point or modular rank promotion"
        and meta["runtime_versions"]
        == {"python": EXPECTED_PYTHON_VERSION, "sympy": EXPECTED_SYMPY_VERSION}
        and meta["process_cpu_seconds"] < meta["process_cpu_budget_seconds"]
        and meta["process_cpu_budget_seconds"] == CPU_BUDGET_SECONDS
        and meta["peak_rss_bytes"] < meta["rss_limit_bytes"]
        and meta["rss_limit_bytes"] == RSS_LIMIT_BYTES
        and meta["peak_rss_measurement"] == peak_rss_measurement()
        and meta["benchmark_used"] is False
        and meta["L10_closure_run"] is False,
    )

    bases, e237_records = source_records()
    matrices: dict[tuple[int, int], MatrixZ] = {}
    representatives: dict[tuple[int, int], list[int]] = {}
    all_cases_hold = True
    for L, m in TARGETS:
        matrix, reps, _coordinates, passed = rebuild_case(
            L,
            m,
            bases[L],
            e237_records[(L, m)],
            certificate["cases"][f"L{L}_m{m}"],
        )
        matrices[(L, m)] = matrix
        representatives[(L, m)] = reps
        all_cases_hold &= passed
    check(
        "independent literal cores, coordinates, Burnside counts, and Smith blocks",
        all_cases_hold,
        "all 283 residual orbit coordinates and all three selected integer matrices rebuilt without importing a producer",
    )

    reps4 = representatives[(9, 4)]
    reps5 = representatives[(9, 5)]
    index5 = {representative: index for index, representative in enumerate(reps5)}
    column_permutation = [
        index5[canonical(complement(representative, 9), 9)]
        for representative in reps4
    ]
    mirror = certificate["particle_hole_L9"]
    mirror_holds = (
        sorted(column_permutation) == list(range(133))
        and column_permutation == mirror["m4_to_m5_column_permutation"]
        and all(
            complement(act(config, 9, element), 9)
            == act(complement(config, 9), 9, element)
            for config in reps4
            for element in GROUP_ELEMENTS
        )
        and all(
            matrices[(9, 5)][row][column_permutation[column]]
            == matrices[(9, 4)][row][column]
            for row in range(133)
            for column in range(133)
        )
    )
    for column, representative in enumerate(reps4):
        record4 = coordinate_record(representative, 9, 4)
        record5 = coordinate_record(reps5[column_permutation[column]], 9, 5)
        mirror_holds &= (
            record5["counts"]["E"] == record4["counts"]["D"]
            and record5["counts"]["D"] == record4["counts"]["E"]
            and sorted((record5["counts"]["T"], record5["counts"]["B"]))
            == sorted((record4["counts"]["T"], record4["counts"]["B"]))
            and record5["defect_position_mod_reflection"]
            == record4["defect_position_mod_reflection"]
        )
    check("independent particle-hole bijection and matrix identity", mirror_holds)

    central_holds = True
    index8 = {representative: index for index, representative in enumerate(representatives[(8, 4)])}
    for stored in certificate["central_L8_inheritance"]:
        m = int(stored["source_case"][1])
        rows = [int(value) for value in stored["central_row_indices_in_L9_core"]]
        columns = [int(value) for value in stored["central_column_indices_in_L9_core"]]
        expected_rows = [
            row
            for row, matrix_row in enumerate(matrices[(9, m)])
            if {
                coordinate_record(representatives[(9, m)][column], 9, m)[
                    "defect_position_mod_reflection"
                ]
                for column, value in enumerate(matrix_row)
                if value
            }
            == {4}
        ]
        expected_columns = [
            column
            for column, representative in enumerate(representatives[(9, m)])
            if coordinate_record(representative, 9, m)[
                "defect_position_mod_reflection"
            ]
            == 4
        ]
        central_holds &= rows == expected_rows and columns == expected_columns
        target_columns = []
        for column in columns:
            reduced = delete_rung(representatives[(9, m)][column], 4)
            if m == 5:
                reduced = complement(reduced, 8)
            target_columns.append(index8[canonical(reduced, 8)])
        central_holds &= target_columns == stored["column_target_indices_in_L8_core"]
        order = sorted(range(len(target_columns)), key=target_columns.__getitem__)
        central = [
            [matrices[(9, m)][row][columns[column]] for column in order]
            for row in rows
        ]
        targets = [int(value) for value in stored["row_target_indices_in_L8_core"]]
        signs = [int(value) for value in stored["row_signs"]]
        central_holds &= sorted(targets) == list(range(17))
        central_holds &= all(sign in (-1, 1) for sign in signs)
        central_holds &= all(
            central[row]
            == [signs[row] * value for value in matrices[(8, 4)][targets[row]]]
            for row in range(17)
        )
    check("independent central L8 signed-permutation inheritance", central_holds)

    A8_record = certificate["cases"]["L8_m4"]["packing_polar_partition"]["packing_matrix"]
    A9_record = certificate["cases"]["L9_m4"]["packing_polar_partition"]["packing_matrix"]
    obstruction = certificate["eight_copy_recurrence_obstruction"]
    recurrence_holds = (
        int(A8_record["shape"][0]) == 16
        and int(A9_record["shape"][0]) == 128
        and abs(int(A8_record["determinant_Z"])) ** 8
        != abs(int(A9_record["determinant_Z"]))
        and obstruction["determinants_disagree"] is True
        and obstruction["first_shell_failure"]["normalized_smith_multiplicities"]
        == {"1": 23, "2": 5, "6": 1}
        and obstruction["scope"]
        == "this refutes only the named count-driven eight-copy/dyadic-shell recurrence, not every possible block recurrence"
    )
    check("independent first recurrence obstruction", recurrence_holds)

    checked_coordinates, defect_holds = elementary_control()
    check(
        "independent elementary all-L algebra control through L9",
        defect_holds
        and checked_coordinates
        == certificate["finite_decoder_control"][
            "checked_even_particle_coordinates_L1_through_L9"
        ],
        str(checked_coordinates),
    )
    expected_elementary = {
        "tag": "[LEMMA]",
        "assumptions": "L>=1, 0<=m<=L, and a rung word with exactly 2m occupied sites",
        "partition": "{0,...,L-1}=T disjoint_union B disjoint_union E disjoint_union D",
        "particle_identity": "2m=|T|+|B|+2|D|",
        "rung_identity": "L=|T|+|B|+|E|+|D|",
        "defect_identity": "|E|-|D|=L-2m",
        "inequality": "q=|E|+|D|>=|L-2m|",
        "equality_criterion": "q=|L-2m| iff min(|E|,|D|)=0; equivalently there is no cancelling empty-double pair and every rung beyond the forced defect set is singly occupied",
    }
    expected_finite = {
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
    expected_scope_object = {
        "tag": "[SCOPE]",
        "claim": EXPECTED_SCOPE,
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
    check(
        "theorem and finite/open scope are pinned",
        data["claim_tag"] == EXPECTED_CLAIM_TAG
        and data["elementary_defect_theorem"] == expected_elementary
        and data["finite_residual_theorem"] == expected_finite
        and data["scope"] == expected_scope_object,
    )

    elapsed = time.process_time() - started
    current_rss = peak_rss_bytes()
    check(
        "clean-room verifier resource wall",
        elapsed < CPU_BUDGET_SECONDS and current_rss < RSS_LIMIT_BYTES,
        f"CPU={elapsed:.6f}/{CPU_BUDGET_SECONDS}s RSS={current_rss}/{RSS_LIMIT_BYTES} via {peak_rss_measurement()}",
    )
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: W-law support-defect certificate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
