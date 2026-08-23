#!/usr/bin/env python3
"""Independent verifier for the exact Pluecker/Hirota certificate.

[COMPUTATION] This file imports no experiment/producer module.  It rebuilds all
five fixed-boundary tensors from ``ising.lattices`` using a binomial expansion
that differs from the producer's convolution engine, independently repeats the
terminal Walsh transform, and redoes the bounded multivariate coefficient
search by recursive edge inclusion.

[UNRESOLVED] Passing this verifier establishes only the finite tensors and the
explicit 125-candidate ansatz recorded in the artifact.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import sys
import time
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from ising.lattices import Lattice, hyperrect  # noqa: E402

ARTIFACT = ROOT / "results" / "integrability" / "plucker_hirota.json"
COEFFICIENTS = (1, -1, 1, -1)
BOUND = 2
CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 2_000_000_000
STARTED = time.process_time()
PASSED: list[str] = []
FAILED: list[str] = []

CASES = (
    (
        "planar_training",
        (2, 3),
        ((0, 0), (0, 2), (1, 2), (1, 0)),
    ),
    (
        "planar_holdout",
        (3, 3),
        ((0, 0), (0, 2), (2, 2), (2, 0)),
    ),
    (
        "cofacial_cube_scope_control",
        (2, 2, 2),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
    ),
    (
        "skew_cube_certificate",
        (2, 2, 2),
        ((0, 0, 0), (0, 0, 1), (1, 1, 1), (1, 1, 0)),
    ),
    (
        "surface_slab_holdout",
        (2, 2, 3),
        ((0, 0, 0), (0, 0, 2), (0, 1, 2), (0, 1, 0)),
    ),
)


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget(stage: str) -> None:
    cpu = time.process_time() - STARTED
    rss = max_rss_bytes()
    if cpu >= CPU_BUDGET_SECONDS:
        raise RuntimeError(f"{stage}: process CPU {cpu:.3f}>={CPU_BUDGET_SECONDS}")
    if rss >= RSS_CAP_BYTES:
        raise RuntimeError(f"{stage}: RSS {rss}>={RSS_CAP_BYTES}")


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(
        f"[{'PASS' if condition else 'FAIL'}] {name}" + (f": {detail}" if detail else ""),
        flush=True,
    )
    (PASSED if condition else FAILED).append(name)


def trim(poly: Sequence[int]) -> list[int]:
    last = len(poly)
    while last > 1 and int(poly[last - 1]) == 0:
        last -= 1
    return [int(value) for value in poly[:last]]


def multiply(a: Sequence[int], b: Sequence[int]) -> list[int]:
    answer = [0] * (len(a) + len(b) - 1)
    for total_degree in range(len(answer)):
        lower = max(0, total_degree - len(b) + 1)
        upper = min(len(a) - 1, total_degree)
        answer[total_degree] = sum(
            int(a[left_degree]) * int(b[total_degree - left_degree])
            for left_degree in range(lower, upper + 1)
        )
    return trim(answer)


def power(base: Sequence[int], exponent: int) -> list[int]:
    answer = [1]
    for _ in range(exponent):
        answer = multiply(answer, base)
    return answer


def combine(channels: Sequence[Sequence[int]], coefficients: Sequence[int]) -> list[int]:
    size = max(len(channel) for channel in channels)
    answer = [0] * size
    for channel, coefficient in zip(channels, coefficients, strict=True):
        for degree, value in enumerate(channel):
            answer[degree] += int(coefficient) * int(value)
    return trim(answer)


def configuration_kernel(aligned: int, disagreed: int) -> list[int]:
    """Coefficient formula for ``(1+v)^aligned (1-v)^disagreed``."""

    degree = aligned + disagreed
    return [
        sum(
            comb(aligned, plus_degree)
            * comb(disagreed, total_degree - plus_degree)
            * (-1) ** (total_degree - plus_degree)
            for plus_degree in range(
                max(0, total_degree - disagreed), min(aligned, total_degree) + 1
            )
        )
        for total_degree in range(degree + 1)
    ]


def fixed_boundary_tensor(
    lat: Lattice, terminal_coordinates: Sequence[Sequence[int]]
) -> tuple[list[list[int]], tuple[int, ...]]:
    terminals = tuple(lat.index(tuple(coord)) for coord in terminal_coordinates)
    terminal_set = set(terminals)
    interior = tuple(site for site in range(lat.n_sites) if site not in terminal_set)
    kernels = [
        configuration_kernel(aligned, lat.n_bonds - aligned)
        for aligned in range(lat.n_bonds + 1)
    ]
    tensor: list[list[int]] = []

    for boundary_mask in range(16):
        row = [0] * (lat.n_bonds + 1)
        fixed_plus = {
            terminal: bool((boundary_mask >> position) & 1)
            for position, terminal in enumerate(terminals)
        }
        for interior_mask in range(1 << len(interior)):
            plus_sites = set(
                site
                for position, site in enumerate(interior)
                if (interior_mask >> position) & 1
            )
            plus_sites.update(site for site, is_plus in fixed_plus.items() if is_plus)
            aligned = sum(
                (left in plus_sites) == (right in plus_sites)
                for left, right in lat.bonds
            )
            for degree, value in enumerate(kernels[aligned]):
                row[degree] += value
        tensor.append(trim(row))
        budget(f"fixed boundary tensor {lat.shape}: mask={boundary_mask}")
    return tensor, terminals


def normalized_walsh(tensor: Sequence[Sequence[int]], site_count: int) -> list[list[int]]:
    divisor = 1 << site_count
    width = max(len(row) for row in tensor)
    answer: list[list[int]] = []
    for subset in range(16):
        row = [0] * width
        for boundary_mask, source in enumerate(tensor):
            minus_selected = subset.bit_count() - (subset & boundary_mask).bit_count()
            sign = -1 if minus_selected & 1 else 1
            for degree, value in enumerate(source):
                row[degree] += sign * int(value)
        if any(value % divisor for value in row):
            raise AssertionError("nonintegral normalized Walsh coefficient")
        answer.append(trim([value // divisor for value in row]))
    return answer


def four_channels(currents: Sequence[Sequence[int]]) -> list[list[int]]:
    return [
        multiply(currents[0], currents[15]),
        multiply(currents[3], currents[12]),
        multiply(currents[5], currents[10]),
        multiply(currents[9], currents[6]),
    ]


def residual(currents: Sequence[Sequence[int]]) -> list[int]:
    return combine(four_channels(currents), COEFFICIENTS)


def evaluate(poly: Sequence[int], value: Fraction) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(poly):
        answer = answer * value + int(coefficient)
    return answer


def monomial(degree: int, coefficient: int) -> list[int]:
    return [0] * degree + [coefficient]


def independent_multivariate_currents(
    lat: Lattice, terminals: Sequence[int]
) -> list[dict[int, int]]:
    """Recursive include/exclude construction with independent base-three edge variables."""

    if lat.n_bonds > 12:
        raise ValueError("independent multivariate verifier is intentionally capped at 12 edges")
    terminal_position = {site: position for position, site in enumerate(terminals)}
    terminal_set = set(terminals)
    interior_mask = sum(1 << site for site in range(lat.n_sites) if site not in terminal_set)
    answer: list[dict[int, int]] = [{} for _ in range(16)]
    edge_codes = [3**index for index in range(lat.n_bonds)]

    def visit(edge_index: int, parity: int, code: int) -> None:
        if edge_index == lat.n_bonds:
            if parity & interior_mask:
                return
            boundary_mask = sum(
                1 << position
                for site, position in terminal_position.items()
                if (parity >> site) & 1
            )
            answer[boundary_mask][code] = answer[boundary_mask].get(code, 0) + 1
            return
        visit(edge_index + 1, parity, code)
        left, right = lat.bonds[edge_index]
        visit(
            edge_index + 1,
            parity ^ (1 << left) ^ (1 << right),
            code + edge_codes[edge_index],
        )

    visit(0, 0, 0)
    budget(f"multivariate recursion {lat.shape}")
    return answer


def sparse_product(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    answer: dict[int, int] = {}
    for left_code, left_value in a.items():
        for right_code, right_value in b.items():
            code = left_code + right_code
            answer[code] = answer.get(code, 0) + left_value * right_value
    return {code: value for code, value in answer.items() if value}


def sparse_channels(currents: Sequence[dict[int, int]]) -> list[dict[int, int]]:
    return [
        sparse_product(currents[0], currents[15]),
        sparse_product(currents[3], currents[12]),
        sparse_product(currents[5], currents[10]),
        sparse_product(currents[9], currents[6]),
    ]


def sparse_combine(
    channels: Sequence[dict[int, int]], coefficients: Sequence[int]
) -> dict[int, int]:
    keys = set().union(*(channel.keys() for channel in channels))
    return {
        code: value
        for code in keys
        if (
            value := sum(
                int(coefficient) * channel.get(code, 0)
                for coefficient, channel in zip(coefficients, channels, strict=True)
            )
        )
    }


def ternary_digits(code: int, count: int) -> list[int]:
    digits = []
    for _ in range(count):
        digits.append(code % 3)
        code //= 3
    if code:
        raise AssertionError("monomial exceeds edge-variable count")
    return digits


def specialize_uniform(poly: dict[int, int], edge_count: int) -> list[int]:
    answer = [0] * (2 * edge_count + 1)
    for code, coefficient in poly.items():
        answer[sum(ternary_digits(code, edge_count))] += coefficient
    return trim(answer)


def sparse_digest(poly: dict[int, int]) -> str:
    payload = "\n".join(f"{code}:{poly[code]}" for code in sorted(poly)).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    ok(
        "[COMPUTATION] artifact shape",
        set(artifact) == {"meta", "data", "checks"},
    )
    ok(
        "[COMPUTATION] stored producer checks",
        bool(artifact["checks"]) and all(bool(row["passed"]) for row in artifact["checks"]),
    )

    rebuilt: dict[str, dict[str, object]] = {}
    for key, shape, terminal_coordinates in CASES:
        lat = hyperrect(shape, periodic=False)
        tensor, terminals = fixed_boundary_tensor(lat, terminal_coordinates)
        currents = normalized_walsh(tensor, lat.n_sites)
        case_residual = residual(currents)
        stored = artifact["data"][key]
        rebuilt[key] = {
            "lat": lat,
            "terminals": terminals,
            "tensor": tensor,
            "currents": currents,
            "residual": case_residual,
        }
        ok(
            f"[COMPUTATION] {key}: exact tensor",
            tensor == stored["boundary_partition_tensor_coefficients_ascending"],
        )
        ok(
            f"[COMPUTATION] {key}: exact Walsh currents",
            currents == stored["normalized_walsh_current_coefficients_ascending"],
        )
        ok(
            f"[COMPUTATION] {key}: exact residual",
            case_residual == stored["plucker_residual_coefficients_ascending"],
        )

    cube_factor = multiply(monomial(6, 4), power([1, 0, -1], 6))
    slab_factor = multiply(
        multiply(
            multiply(monomial(8, 4), power([1, 0, -1], 8)),
            [1, 0, 1],
        ),
        [3, 0, 6, 0, -1],
    )
    ok(
        "[COMPUTATION] two planar controls vanish",
        rebuilt["planar_training"]["residual"] == [0]
        and rebuilt["planar_holdout"]["residual"] == [0],
    )
    ok(
        "[COMPUTATION] cofacial cube scope control vanishes",
        rebuilt["cofacial_cube_scope_control"]["residual"] == [0],
    )
    ok(
        "[THEOREM] skew cube factor certificate",
        rebuilt["skew_cube_certificate"]["residual"] == cube_factor
        and evaluate(cube_factor, Fraction(1, 2)) == Fraction(729, 65536),
    )
    ok(
        "[THEOREM] surface slab factor certificate",
        rebuilt["surface_slab_holdout"]["residual"] == slab_factor
        and evaluate(slab_factor, Fraction(1, 2))
        == Fraction(2329155, 268435456),
    )

    training_lat = rebuilt["planar_training"]["lat"]
    training_terminals = rebuilt["planar_training"]["terminals"]
    training_sparse = independent_multivariate_currents(training_lat, training_terminals)
    training_channels = sparse_channels(training_sparse)
    candidates = [(1, *tail) for tail in itertools.product(range(-BOUND, BOUND + 1), repeat=3)]
    survivors = [
        coefficients
        for coefficients in candidates
        if not sparse_combine(training_channels, coefficients)
    ]
    ok(
        "[COMPUTATION] independent 125-candidate planar search",
        len(candidates) == 125 and survivors == [COEFFICIENTS],
    )

    holdout_lat = rebuilt["planar_holdout"]["lat"]
    holdout_terminals = rebuilt["planar_holdout"]["terminals"]
    holdout_sparse = independent_multivariate_currents(holdout_lat, holdout_terminals)
    ok(
        "[COMPUTATION] independent-edge planar holdout",
        not sparse_combine(sparse_channels(holdout_sparse), COEFFICIENTS),
    )

    cube_lat = rebuilt["skew_cube_certificate"]["lat"]
    cube_terminals = rebuilt["skew_cube_certificate"]["terminals"]
    cube_sparse = independent_multivariate_currents(cube_lat, cube_terminals)
    cube_sparse_channels = sparse_channels(cube_sparse)
    cube_sparse_residual = sparse_combine(cube_sparse_channels, COEFFICIENTS)
    joint_survivors = [
        coefficients
        for coefficients in candidates
        if not sparse_combine(training_channels, coefficients)
        and not sparse_combine(cube_sparse_channels, coefficients)
    ]
    stored_multi = artifact["data"]["skew_cube_multivariate_certificate"]
    ok(
        "[COMPUTATION] independent multivariate cube digest",
        len(cube_sparse_residual) == stored_multi["nonzero_term_count"]
        and sparse_digest(cube_sparse_residual) == stored_multi["canonical_sha256"]
        and specialize_uniform(cube_sparse_residual, cube_lat.n_bonds) == cube_factor,
    )
    stored_witness = stored_multi["first_minimum_degree_witness"]
    witness_code = int(stored_witness["monomial_code_base3"])
    ok(
        "[COMPUTATION] stored nonzero cube monomial witness",
        cube_sparse_residual.get(witness_code) == int(stored_witness["coefficient"])
        and ternary_digits(witness_code, cube_lat.n_bonds)
        == stored_witness["edge_exponents"],
    )
    ok(
        "[THEOREM] independently bounded planar/cube no-go",
        not joint_survivors
        and artifact["data"]["bounded_identity_search"][
            "joint_planar_training_and_cube_survivors"
        ]
        == [],
    )

    scope = artifact["data"]["conclusion"]["scope"]
    ok(
        "[UNRESOLVED] scope guard",
        "arbitrary determinant/Pfaffian" in scope
        and "all-size identity" in scope
        and "other terminal orders" in scope,
    )
    ok(
        "[COMPUTATION] verifier resource budget",
        time.process_time() - STARTED < CPU_BUDGET_SECONDS
        and max_rss_bytes() < RSS_CAP_BYTES,
        (
            f"CPU={time.process_time() - STARTED:.3f}s; "
            f"RSS={max_rss_bytes()}"
        ),
    )

    if FAILED:
        raise SystemExit(f"FAIL: {', '.join(FAILED)}")
    print(f"OK: {len(PASSED)} independent Pluecker/Hirota checks", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
