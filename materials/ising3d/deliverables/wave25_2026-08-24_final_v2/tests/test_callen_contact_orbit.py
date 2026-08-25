#!/usr/bin/env python3
"""Clean-room verifier for the e250 contact-inclusive Callen orbit theorem."""
from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import platform
import resource
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "correlations" / "callen_contact_orbit.json"
PRODUCER = ROOT / "experiments" / "e250_callen_contact_orbit.py"
SELF = ROOT / "tests" / "test_callen_contact_orbit.py"
PROOF = ROOT / "proofs" / "callen_contact_orbit.md"
E240_PRODUCER = ROOT / "experiments" / "e240_callen_identity_system.py"
E240_ARTIFACT = ROOT / "results" / "correlations" / "callen_identity_system.json"
E243_PRODUCER = ROOT / "experiments" / "e243_callen_termwise_closure.py"
E246_PRODUCER = ROOT / "experiments" / "e246_callen_walsh_compression.py"
CPU_LIMIT_SECONDS = 60.0
RSS_LIMIT_BYTES = 2 * 1024**3
EXPECTED_INHERITED_SHA256 = {
    "experiments/e240_callen_identity_system.py": "f810172970a797d3717f66825d39a60fcc6755d8ec13decdde1c951840ad7ed9",
    "results/correlations/callen_identity_system.json": "c3346266fce6a03e6bb2c8466c50cbcef95f1265700bfdcf72c76e0aced06861",
    "experiments/e243_callen_termwise_closure.py": "71c5f884c8597b4890e3699ba6d19cbf7cfdf987aebc8315a0b334eefd252a28",
    "experiments/e246_callen_walsh_compression.py": "899744f1a9c6898e65c18c52c0b8d7c7b1fab34747c40c0dba030c04b6b58d52",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coord(index: int) -> tuple[int, int, int]:
    return index // 9, index // 3 % 3, index % 3


def site(value: tuple[int, int, int]) -> int:
    return value[0] * 9 + value[1] * 3 + value[2]


def independent_group() -> tuple[tuple[int, ...], ...]:
    """Build S3 wr S3 directly, rather than signed affine maps."""
    points = tuple(coord(index) for index in range(27))
    value_permutations = tuple(itertools.permutations(range(3)))
    maps = []
    for axes in itertools.permutations(range(3)):
        for coordinate_maps in itertools.product(value_permutations, repeat=3):
            maps.append(
                tuple(
                    site(
                        tuple(
                            coordinate_maps[output_axis][point[axes[output_axis]]]
                            for output_axis in range(3)
                        )
                    )
                    for point in points
                )
            )
    result = tuple(dict.fromkeys(maps))
    if len(result) != 1296:
        raise AssertionError(f"independent group order is {len(result)}")
    return result


def image(mask: int, mapping: tuple[int, ...]) -> int:
    output = 0
    for index in range(27):
        if mask >> index & 1:
            output |= 1 << mapping[index]
    return output


def canonical(mask: int, maps: tuple[tuple[int, ...], ...]) -> int:
    return min(image(mask, mapping) for mapping in maps)


def neighbourhood() -> tuple[int, ...]:
    origin = coord(0)
    return tuple(
        index
        for index in range(1, 27)
        if sum(left != right for left, right in zip(origin, coord(index))) == 1
    )


def odd_masks(neighbours: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(1 << neighbours[index] for index in range(6) if subset >> index & 1)
        for subset in range(64)
        if subset.bit_count() & 1
    )


def rebuild_rows(
    maps: tuple[tuple[int, ...], ...]
) -> tuple[list[int], list[dict[int, tuple[int, int, int, int]]], list[int], list[int]]:
    stabilizer = tuple(mapping for mapping in maps if mapping[0] == 0)
    local = odd_masks(neighbourhood())
    representatives = sorted(
        {min(image(mask, mapping) for mapping in stabilizer) for mask in local}
    )
    rows = []
    columns: set[int] = set()
    orbit_sizes = []
    coefficient_slot = {1: 1, 3: 2, 5: 3}
    for factor in representatives:
        mutable: dict[int, list[int]] = {canonical(1 | factor, maps): [1, 0, 0, 0]}
        for expansion in local:
            column = canonical(factor ^ expansion, maps)
            mutable.setdefault(column, [0, 0, 0, 0])[coefficient_slot[expansion.bit_count()]] -= 1
        row = {key: tuple(value) for key, value in mutable.items() if any(value)}
        rows.append(row)
        columns.update(row)
        orbit_sizes.append(
            27 * len({image(factor, mapping) for mapping in stabilizer})
        )
    return representatives, rows, sorted(columns), orbit_sizes


def support(mask: int) -> dict[str, object]:
    return {
        "mask": mask,
        "size": mask.bit_count(),
        "coordinates": [list(coord(index)) for index in range(27) if mask >> index & 1],
    }


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        library = ctypes.CDLL("/usr/lib/libSystem.B.dylib")

        class Info(ctypes.Structure):
            _fields_ = [
                ("virtual", ctypes.c_uint64),
                ("resident", ctypes.c_uint64),
                ("resident_max", ctypes.c_uint64),
                ("u_s", ctypes.c_int32),
                ("u_us", ctypes.c_int32),
                ("s_s", ctypes.c_int32),
                ("s_us", ctypes.c_int32),
                ("policy", ctypes.c_int32),
                ("suspend", ctypes.c_int32),
            ]

        library.mach_task_self.restype = ctypes.c_uint32
        library.task_info.argtypes = [
            ctypes.c_uint32,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        library.task_info.restype = ctypes.c_int
        info = Info()
        count = ctypes.c_uint32(ctypes.sizeof(info) // 4)
        result = library.task_info(library.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count))
        if result:
            raise RuntimeError(f"task_info failed: {result}")
        return int(info.resident_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    failures: list[str] = []
    expected_hashes = {
        "experiments/e250_callen_contact_orbit.py": digest(PRODUCER),
        "tests/test_callen_contact_orbit.py": digest(SELF),
        "proofs/callen_contact_orbit.md": digest(PROOF),
        "experiments/e240_callen_identity_system.py": digest(E240_PRODUCER),
        "results/correlations/callen_identity_system.json": digest(E240_ARTIFACT),
        "experiments/e243_callen_termwise_closure.py": digest(E243_PRODUCER),
        "experiments/e246_callen_walsh_compression.py": digest(E246_PRODUCER),
    }
    if artifact["meta"]["source_sha256"] != expected_hashes:
        failures.append("source hash map mismatch")
    if artifact["meta"]["frozen_inherited_sha256"] != EXPECTED_INHERITED_SHA256:
        failures.append("frozen inherited hash map mismatch")
    if any(
        expected_hashes[path] != expected
        for path, expected in EXPECTED_INHERITED_SHA256.items()
    ):
        failures.append("current inherited dependency hash mismatch")
    if artifact["meta"]["data_sha256"] != hashlib.sha256(
        json.dumps(artifact["data"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest():
        failures.append("complete data hash mismatch")

    maps = independent_group()
    representatives, rows, columns, orbit_sizes = rebuild_rows(maps)
    if [support(mask) for mask in representatives] != artifact["data"]["row_representatives"]:
        failures.append("row representatives mismatch")
    if [support(mask) for mask in columns] != artifact["data"]["column_representatives"]:
        failures.append("column representatives mismatch")
    if orbit_sizes != [162, 324, 216, 162] or sum(orbit_sizes) != 864:
        failures.append(f"row orbit census mismatch: {orbit_sizes}")
    expected_geometry = {
        "shape": [3, 3, 3],
        "vertices": 27,
        "degree": 6,
        "space_group_order": 1296,
        "origin_stabilizer_order": 48,
        "raw_row_count": 864,
        "row_orbit_sizes": [162, 324, 216, 162],
    }
    if artifact["data"]["geometry"] != expected_geometry:
        failures.append("geometry semantic record mismatch")

    empty = canonical(0, maps)
    nearest_pair = canonical(1 | (1 << neighbourhood()[0]), maps)
    expected_allowed = {
        "empty": support(empty),
        "nearest_pair": support(nearest_pair),
    }
    if artifact["data"]["allowed_columns"] != expected_allowed:
        failures.append("allowed-column semantic record mismatch")
    nonallowed = set(columns) - {empty, nearest_pair}

    pivots = [entry["mask"] for entry in artifact["data"]["nonallowed_pivot_columns"]]
    if len(pivots) != 4 or len(set(pivots)) != 4 or not set(pivots) <= nonallowed:
        failures.append("pivot columns are not four distinct independently nonallowed columns")
    c1, c3, c5, v = sp.symbols("c1 c3 c5 v")
    symbols = (1, c1, c3, c5)
    matrix = sp.Matrix(
        [
            [sum(value * symbol for value, symbol in zip(row.get(column, (0, 0, 0, 0)), symbols))
             for column in pivots]
            for row in rows
        ]
    )
    expected_matrix_record = [
        [str(matrix[row, column]) for column in range(4)] for row in range(4)
    ]
    if artifact["data"]["formal_minor_matrix"] != expected_matrix_record:
        failures.append("stored literal formal minor matrix mismatch")
    determinant = sp.factor(matrix.det())
    expected_formal = 4 * (c1 - c5) * (c1 + 4 * c3 + c5)
    if sp.expand(determinant - expected_formal):
        failures.append(f"formal determinant mismatch: {determinant}")
    if str(determinant) != artifact["data"]["formal_determinant"]:
        failures.append("stored formal determinant mismatch")

    coefficient_data = json.loads(E240_ARTIFACT.read_text())["data"]["local_master_identity"]["coefficients"]

    def expression(name: str) -> sp.Expr:
        record = coefficient_data[name]
        numerator = sum(value * v**degree for degree, value in enumerate(record["numerator_coefficients_ascending"]))
        denominator = sum(value * v**degree for degree, value in enumerate(record["denominator_coefficients_ascending"]))
        return sp.cancel(numerator / denominator)

    determinant_v = sp.factor(
        sp.cancel(determinant.subs({c1: expression("c1"), c3: expression("c3"), c5: expression("c5")}))
    )
    positive_one = v**8 + 16 * v**6 + 30 * v**4 + 16 * v**2 + 1
    positive_two = v**8 + 8 * v**6 + 14 * v**4 + 8 * v**2 + 1
    denominator = (1 + v**2) * (1 + 6 * v**2 + v**4) * (1 + 14 * v**2 + v**4)
    expected_v = sp.factor(4 * v**2 * positive_one * positive_two / denominator**2)
    if sp.cancel(determinant_v - expected_v):
        failures.append("all-v determinant factor mismatch")
    if str(determinant_v) != artifact["data"]["coupling_determinant"]["expression"]:
        failures.append("stored all-v determinant mismatch")
    expected_numerator, expected_denominator = expected_v.as_numer_denom()
    coupling_record = artifact["data"]["coupling_determinant"]
    if coupling_record != {
        "expression": str(determinant_v),
        "numerator": str(sp.factor(expected_numerator)),
        "denominator": str(sp.factor(expected_denominator)),
        "positive_factors": [str(positive_one), str(positive_two)],
    }:
        failures.append("complete coupling-determinant record mismatch")
    expected_theorem = {
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
    }
    if artifact["data"]["theorem"] != expected_theorem:
        failures.append("theorem/scope semantic record mismatch")
    if artifact["data"]["control_v"] != "1/3" or artifact["data"]["control_rank"] != 4:
        failures.append("control value/rank mismatch")
    if any(
        coefficient <= 0
        for factor in (positive_one, positive_two)
        for _monomial, coefficient in sp.Poly(factor, v).terms()
    ):
        failures.append("claimed positive factor has a nonpositive nonzero coefficient")
    expected_check_names = {
        "C0_provenance_gate",
        "C1_complete_row_orbits",
        "C2_complete_column_orbits",
        "C3_exact_coefficient_control",
        "C4_control_nonallowed_rank",
        "C5_formal_minor_factor",
        "C6_all_positive_coupling_factor",
        "C7_full_neighbour_row_included",
        "C8_cpu_budget",
        "C9_rss_budget",
    }
    if {check["name"] for check in artifact["checks"]} != expected_check_names:
        failures.append("producer semantic check-name set mismatch")

    if not all(check["passed"] for check in artifact["checks"]):
        failures.append("producer artifact contains a failed check")
    elapsed = time.process_time() - started
    rss = peak_rss_bytes()
    if elapsed >= CPU_LIMIT_SECONDS:
        failures.append(f"CPU {elapsed:.6f}s reached {CPU_LIMIT_SECONDS}s")
    if rss >= RSS_LIMIT_BYTES:
        failures.append(f"RSS {rss} reached {RSS_LIMIT_BYTES}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print(f"PASS: e250 contact Callen orbit theorem; CPU={elapsed:.6f}s RSS={rss} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
