#!/usr/bin/env python3
"""Exact contact-inclusive cubic-orbit obstruction for Callen rows on C3^3."""
from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import platform
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments import e240_callen_identity_system as callen  # noqa: E402

SCRIPT = "experiments/e250_callen_contact_orbit.py"
OUTPUT = ROOT / "results" / "correlations" / "callen_contact_orbit.json"
VERIFIER = ROOT / "tests" / "test_callen_contact_orbit.py"
PROOF = ROOT / "proofs" / "callen_contact_orbit.md"
E240_PRODUCER = ROOT / "experiments" / "e240_callen_identity_system.py"
E240_ARTIFACT = ROOT / "results" / "correlations" / "callen_identity_system.json"
E243_PRODUCER = ROOT / "experiments" / "e243_callen_termwise_closure.py"
E246_PRODUCER = ROOT / "experiments" / "e246_callen_walsh_compression.py"
CPU_LIMIT_SECONDS = 60.0
RSS_LIMIT_BYTES = 2 * 1024**3
V = sp.symbols("v")
C1, C3, C5 = sp.symbols("c1 c3 c5")
FORMAL_INDEX = {1: 1, 3: 2, 5: 3}
EXPECTED_INHERITED_SHA256 = {
    "experiments/e240_callen_identity_system.py": "f810172970a797d3717f66825d39a60fcc6755d8ec13decdde1c951840ad7ed9",
    "results/correlations/callen_identity_system.json": "c3346266fce6a03e6bb2c8466c50cbcef95f1265700bfdcf72c76e0aced06861",
    "experiments/e243_callen_termwise_closure.py": "71c5f884c8597b4890e3699ba6d19cbf7cfdf987aebc8315a0b334eefd252a28",
    "experiments/e246_callen_walsh_compression.py": "899744f1a9c6898e65c18c52c0b8d7c7b1fab34747c40c0dba030c04b6b58d52",
}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        system = ctypes.CDLL("/usr/lib/libSystem.B.dylib", use_errno=True)

        class BasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_seconds", ctypes.c_int32),
                ("user_microseconds", ctypes.c_int32),
                ("system_seconds", ctypes.c_int32),
                ("system_microseconds", ctypes.c_int32),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        system.mach_task_self.restype = ctypes.c_uint32
        system.task_info.argtypes = [
            ctypes.c_uint32,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        system.task_info.restype = ctypes.c_int
        info = BasicInfo()
        count = ctypes.c_uint32(ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_uint32))
        result = system.task_info(
            system.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count)
        )
        if result:
            raise RuntimeError(f"task_info failed with kern_return={result}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def coordinate(index: int) -> tuple[int, int, int]:
    return index // 9, (index // 3) % 3, index % 3


def index_of(value: tuple[int, int, int]) -> int:
    return (value[0] % 3) * 9 + (value[1] % 3) * 3 + value[2] % 3


def space_group_maps() -> tuple[tuple[int, ...], ...]:
    coordinates = tuple(coordinate(index) for index in range(27))
    linear_maps: list[tuple[int, ...]] = []
    for axes in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            linear_maps.append(
                tuple(
                    index_of(tuple(signs[axis] * point[axes[axis]] for axis in range(3)))
                    for point in coordinates
                )
            )
    linear_maps = list(dict.fromkeys(linear_maps))
    if len(linear_maps) != 48:
        raise AssertionError(f"origin stabilizer has {len(linear_maps)} maps")
    maps = []
    for linear in linear_maps:
        for translation in coordinates:
            maps.append(
                tuple(
                    index_of(
                        tuple(
                            coordinate(linear[site])[axis] + translation[axis]
                            for axis in range(3)
                        )
                    )
                    for site in range(27)
                )
            )
    maps = tuple(dict.fromkeys(maps))
    if len(maps) != 1296:
        raise AssertionError(f"space group has {len(maps)} maps")
    return maps


def apply_map(mask: int, mapping: tuple[int, ...]) -> int:
    image = 0
    while mask:
        bit = mask & -mask
        site = bit.bit_length() - 1
        image |= 1 << mapping[site]
        mask -= bit
    return image


def canonical_mask(mask: int, maps: tuple[tuple[int, ...], ...]) -> int:
    return min(apply_map(mask, mapping) for mapping in maps)


def neighbours_of_origin() -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                index_of((1, 0, 0)),
                index_of((-1, 0, 0)),
                index_of((0, 1, 0)),
                index_of((0, -1, 0)),
                index_of((0, 0, 1)),
                index_of((0, 0, -1)),
            }
        )
    )


def odd_neighbour_masks(neighbours: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(1 << neighbours[index] for index in range(6) if subset >> index & 1)
        for subset in range(1 << 6)
        if subset.bit_count() % 2
    )


def formal_rows(
    all_maps: tuple[tuple[int, ...], ...],
    stabilizer: tuple[tuple[int, ...], ...],
) -> tuple[list[int], list[dict[int, tuple[int, int, int, int]]], list[int], list[int]]:
    neighbours = neighbours_of_origin()
    odd_masks = odd_neighbour_masks(neighbours)
    row_representatives = sorted(
        {min(apply_map(mask, mapping) for mapping in stabilizer) for mask in odd_masks}
    )
    rows: list[dict[int, tuple[int, int, int, int]]] = []
    columns: set[int] = set()
    for factor_mask in row_representatives:
        mutable: dict[int, list[int]] = {}
        target = canonical_mask(1 | factor_mask, all_maps)
        mutable[target] = [1, 0, 0, 0]
        for local_mask in odd_masks:
            column = canonical_mask(factor_mask ^ local_mask, all_maps)
            record = mutable.setdefault(column, [0, 0, 0, 0])
            record[FORMAL_INDEX[local_mask.bit_count()]] -= 1
        row = {
            column: tuple(coefficients)
            for column, coefficients in mutable.items()
            if any(coefficients)
        }
        rows.append(row)
        columns.update(row)
    row_orbit_sizes = []
    for representative in row_representatives:
        row_orbit_sizes.append(
            len({apply_map(representative, mapping) for mapping in stabilizer}) * 27
        )
    return row_representatives, rows, sorted(columns), row_orbit_sizes


def formal_expression(record: tuple[int, int, int, int]) -> sp.Expr:
    return record[0] + record[1] * C1 + record[2] * C3 + record[3] * C5


def evaluate_record(
    record: tuple[int, int, int, int], coefficients: dict[int, Fraction]
) -> Fraction:
    return (
        Fraction(record[0])
        + record[1] * coefficients[1]
        + record[2] * coefficients[3]
        + record[3] * coefficients[5]
    )


def fraction_rank(matrix: list[list[Fraction]]) -> tuple[int, list[int]]:
    work = [row[:] for row in matrix]
    rank = 0
    pivots: list[int] = []
    if not work:
        return 0, pivots
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [left - value * right for left, right in zip(work[row], work[rank])]
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return rank, pivots


def rational_function_expression(function: callen.RationalFunction) -> sp.Expr:
    numerator = sum(value * V**degree for degree, value in enumerate(function.numerator))
    denominator = sum(value * V**degree for degree, value in enumerate(function.denominator))
    return sp.cancel(numerator / denominator)


def fraction_at(function: callen.RationalFunction, value: Fraction) -> Fraction:
    def evaluate(coefficients: tuple[int, ...]) -> Fraction:
        total = Fraction(0)
        for coefficient in reversed(coefficients):
            total = total * value + coefficient
        return total

    return evaluate(function.numerator) / evaluate(function.denominator)


def support_record(mask: int) -> dict[str, object]:
    return {
        "mask": mask,
        "size": mask.bit_count(),
        "coordinates": [list(coordinate(site)) for site in range(27) if mask >> site & 1],
    }


def add_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    source_paths = (
        ROOT / SCRIPT,
        VERIFIER,
        PROOF,
        E240_PRODUCER,
        E240_ARTIFACT,
        E243_PRODUCER,
        E246_PRODUCER,
    )
    source_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
    }
    add_check(
        checks,
        "C0_provenance_gate",
        set(source_hashes)
        == {
            SCRIPT,
            "tests/test_callen_contact_orbit.py",
            "proofs/callen_contact_orbit.md",
            *EXPECTED_INHERITED_SHA256,
        }
        and all(
            source_hashes[path] == digest
            for path, digest in EXPECTED_INHERITED_SHA256.items()
        ),
        "current producer/verifier/proof plus frozen e240/e243/e246 dependencies",
    )
    all_maps = space_group_maps()
    stabilizer = tuple(mapping for mapping in all_maps if mapping[0] == 0)
    row_representatives, rows, columns, row_orbit_sizes = formal_rows(all_maps, stabilizer)
    add_check(
        checks,
        "C1_complete_row_orbits",
        len(row_representatives) == 4 and sum(row_orbit_sizes) == 864,
        f"row orbits {row_orbit_sizes}",
    )
    add_check(
        checks,
        "C2_complete_column_orbits",
        len(columns) == 9,
        f"column support sizes {[column.bit_count() for column in columns]}",
    )

    coefficients_rf = callen.local_coefficients(6)
    coefficients_at_control = {
        size: fraction_at(function, Fraction(1, 3))
        for size, function in coefficients_rf.items()
    }
    expected_control = {
        1: Fraction(4143, 17680),
        3: Fraction(-27, 1040),
        5: Fraction(243, 17680),
    }
    add_check(
        checks,
        "C3_exact_coefficient_control",
        coefficients_at_control == expected_control,
        str(coefficients_at_control),
    )

    empty_column = canonical_mask(0, all_maps)
    nearest_pair_column = canonical_mask(1 | (1 << neighbours_of_origin()[0]), all_maps)
    nonallowed_columns = [
        column for column in columns if column not in (empty_column, nearest_pair_column)
    ]
    numerical_matrix = [
        [evaluate_record(row.get(column, (0, 0, 0, 0)), coefficients_at_control)
         for column in nonallowed_columns]
        for row in rows
    ]
    rank, pivot_indices = fraction_rank(numerical_matrix)
    pivot_columns = [nonallowed_columns[index] for index in pivot_indices]
    add_check(
        checks,
        "C4_control_nonallowed_rank",
        rank == 4 and len(pivot_columns) == 4,
        f"rank {rank}; pivots {pivot_columns}",
    )

    symbolic_minor = sp.Matrix(
        [
            [formal_expression(row.get(column, (0, 0, 0, 0))) for column in pivot_columns]
            for row in rows
        ]
    )
    control_substitution = {
        C1: sp.Rational(coefficients_at_control[1].numerator, coefficients_at_control[1].denominator),
        C3: sp.Rational(coefficients_at_control[3].numerator, coefficients_at_control[3].denominator),
        C5: sp.Rational(coefficients_at_control[5].numerator, coefficients_at_control[5].denominator),
    }
    if symbolic_minor.det().subs(control_substitution) < 0:
        pivot_columns[-2], pivot_columns[-1] = pivot_columns[-1], pivot_columns[-2]
        symbolic_minor = sp.Matrix(
            [
                [
                    formal_expression(row.get(column, (0, 0, 0, 0)))
                    for column in pivot_columns
                ]
                for row in rows
            ]
        )
    formal_determinant = sp.factor(symbolic_minor.det())
    expected_formal = 4 * (C1 - C5) * (C1 + 4 * C3 + C5)
    add_check(
        checks,
        "C5_formal_minor_factor",
        sp.expand(formal_determinant - expected_formal) == 0,
        str(formal_determinant),
    )

    substitutions = {
        C1: rational_function_expression(coefficients_rf[1]),
        C3: rational_function_expression(coefficients_rf[3]),
        C5: rational_function_expression(coefficients_rf[5]),
    }
    determinant_v = sp.factor(sp.cancel(formal_determinant.subs(substitutions)))
    numerator_v, denominator_v = determinant_v.as_numer_denom()
    numerator_v = sp.factor(numerator_v)
    denominator_v = sp.factor(denominator_v)
    positive_factors = (
        V**8 + 16 * V**6 + 30 * V**4 + 16 * V**2 + 1,
        V**8 + 8 * V**6 + 14 * V**4 + 8 * V**2 + 1,
    )
    common_denominator = (
        (1 + V**2) * (1 + 6 * V**2 + V**4) * (1 + 14 * V**2 + V**4)
    )
    expected_v = sp.factor(4 * V**2 * positive_factors[0] * positive_factors[1] / common_denominator**2)
    add_check(
        checks,
        "C6_all_positive_coupling_factor",
        sp.cancel(determinant_v - expected_v) == 0
        and all(
            all(coefficient > 0 for _monomial, coefficient in sp.Poly(factor, V).terms())
            for factor in positive_factors
        ),
        str(determinant_v),
    )
    add_check(
        checks,
        "C7_full_neighbour_row_included",
        sorted(mask.bit_count() for mask in row_representatives) == [1, 3, 3, 5],
        "the size-five factor is the neighbour-contact reduction of e246's U=N row",
    )
    data_payload = {
        "theorem": {
            "statement": (
                "For every 0<v=tanh(K)<1, the translation-and-full-cubic-symmetry quotient "
                "of all 864 contact Callen rows on C3^3 has no nonzero combination supported "
                "only on normalization and the nearest-neighbour pair orbit."
            ),
            "scope": (
                "linear contact rows with factor B an odd subset of one pivot neighbourhood; "
                "all pivots, the size-five neighbour-contact reduction of e246's U=N row, and "
                "cubic orbit aggregation are included; radius-expanded, nonlinear, auxiliary, "
                "and nonsymmetric zero-average relations are not"
            ),
        },
        "geometry": {
            "shape": [3, 3, 3],
            "vertices": 27,
            "degree": 6,
            "space_group_order": len(all_maps),
            "origin_stabilizer_order": len(stabilizer),
            "raw_row_count": 27 * 32,
            "row_orbit_sizes": row_orbit_sizes,
        },
        "row_representatives": [support_record(mask) for mask in row_representatives],
        "column_representatives": [support_record(mask) for mask in columns],
        "allowed_columns": {
            "empty": support_record(empty_column),
            "nearest_pair": support_record(nearest_pair_column),
        },
        "nonallowed_pivot_columns": [support_record(mask) for mask in pivot_columns],
        "formal_minor_matrix": [
            [str(symbolic_minor[row, column]) for column in range(4)] for row in range(4)
        ],
        "formal_determinant": str(formal_determinant),
        "coupling_determinant": {
            "expression": str(determinant_v),
            "numerator": str(numerator_v),
            "denominator": str(denominator_v),
            "positive_factors": [str(factor) for factor in positive_factors],
        },
        "control_v": "1/3",
        "control_rank": rank,
    }

    cpu = time.process_time() - started
    rss = peak_rss_bytes()
    add_check(checks, "C8_cpu_budget", cpu < CPU_LIMIT_SECONDS, f"{cpu:.6f}s")
    add_check(checks, "C9_rss_budget", rss < RSS_LIMIT_BYTES, f"{rss} bytes")
    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"refusing to write failed artifact: {failed}")

    payload = {
        "meta": {
            "experiment": "e250",
            "script": SCRIPT,
            "source_sha256": source_hashes,
            "frozen_inherited_sha256": EXPECTED_INHERITED_SHA256,
            "data_sha256": canonical_sha256(data_payload),
            "process_cpu_seconds": cpu,
            "peak_rss_bytes": rss,
            "peak_rss_measurement": (
                "mach_task_basic_info.resident_size_max (task_info flavor 20)"
                if platform.system() == "Darwin"
                else "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"
            ),
            "cpu_limit_seconds": CPU_LIMIT_SECONDS,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
        },
        "checks": checks,
        "data": data_payload,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> int:
    payload = run()
    print(json.dumps(payload["meta"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
