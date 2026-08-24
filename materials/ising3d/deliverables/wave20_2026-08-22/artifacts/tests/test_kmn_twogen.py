#!/usr/bin/env python3
"""[COMPUTATION] Clean-room verifier for the K_{m,n} collective-spin certificate."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import resource
import time
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "kmn_twogen.json"
VERIFY_PRIME = 65_537
CPU_CAP_SECONDS = 300.0
RSS_CAP_BYTES = 1_900_000_000
SOURCE_PATHS = (
    "experiments/e194_kmn_reduction.py",
    "experiments/e195_kmn_closure.py",
    "experiments/e196_kmn_theorem.py",
    "tests/test_kmn_twogen.py",
    "proofs/kmn_twogen.md",
)


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    if elapsed > CPU_CAP_SECONDS:
        raise TimeoutError(f"verifier CPU cap at {stage}: {elapsed:.3f}s")
    if peak_rss_bytes() >= RSS_CAP_BYTES:
        raise MemoryError(f"verifier RSS cap at {stage}: {peak_rss_bytes()}")


def check(name: str, condition: object, detail: object = "") -> None:
    if not bool(condition):
        raise AssertionError(f"FAIL {name}: {detail}")
    print(f"PASS {name}: {detail}")


def labels_for_qubits(count: int) -> tuple[int, ...]:
    return tuple(range(count & 1, count + 1, 2))


def dual_spin_pair(label: int) -> tuple[np.ndarray, np.ndarray]:
    """Dual polynomial-basis matrices, deliberately transposed from the producer convention."""
    dimension = label + 1
    x = np.zeros((dimension, dimension), dtype=np.int64)
    for row in range(dimension):
        if row + 1 < dimension:
            x[row, row + 1] = label - row
        if row:
            x[row, row - 1] = row
    z = np.diag([label - 2 * row for row in range(dimension)]).astype(np.int64)
    return x, z


def clean_blocks(left_size: int, right_size: int) -> tuple[tuple[int, int], ...]:
    left = labels_for_qubits(left_size)
    right = labels_for_qubits(right_size)
    if left_size != right_size:
        return tuple((r, s) for r in left for s in right if r or s)
    return tuple(
        (r, s)
        for position, r in enumerate(left)
        for s in right[position:]
        if r or s
    )


def clean_generators(
    left_size: int, right_size: int
) -> tuple[tuple[tuple[int, int], ...], tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    blocks = clean_blocks(left_size, right_size)
    field: list[np.ndarray] = []
    bond: list[np.ndarray] = []
    for left_label, right_label in blocks:
        left_x, left_z = dual_spin_pair(left_label)
        right_x, right_z = dual_spin_pair(right_label)
        left_identity = np.eye(left_label + 1, dtype=np.int64)
        right_identity = np.eye(right_label + 1, dtype=np.int64)
        field.append(
            np.kron(left_x, right_identity) + np.kron(left_identity, right_x)
        )
        bond.append(np.kron(left_z, right_z))
    return blocks, tuple(field), tuple(bond)


def flatten(matrices: tuple[np.ndarray, ...]) -> np.ndarray:
    return np.concatenate([matrix.reshape(-1) for matrix in matrices]).astype(np.int64)


def clean_modular_rank(left_size: int, right_size: int, started: float) -> dict[str, object]:
    blocks, field_raw, bond_raw = clean_generators(left_size, right_size)
    field = tuple(matrix % VERIFY_PRIME for matrix in field_raw)
    bond = tuple(matrix % VERIFY_PRIME for matrix in bond_raw)
    actions = (field, bond)
    pivots: dict[int, np.ndarray] = {}
    words: list[tuple[np.ndarray, ...]] = []
    depths: list[int] = []
    frontier: list[int] = []

    def insert(candidate: tuple[np.ndarray, ...], depth: int) -> bool:
        row = flatten(candidate) % VERIFY_PRIME
        while True:
            support = np.flatnonzero(row)
            if not support.size:
                return False
            pivot = int(support[-1])
            old = pivots.get(pivot)
            if old is None:
                row = row * pow(int(row[pivot]), VERIFY_PRIME - 2, VERIFY_PRIME) % VERIFY_PRIME
                pivots[pivot] = row
                words.append(tuple(matrix.copy() % VERIFY_PRIME for matrix in candidate))
                depths.append(depth)
                return True
            row = (row - int(row[pivot]) * old) % VERIFY_PRIME

    for seed in actions:
        if not insert(seed, 1):
            raise AssertionError("dependent clean-room seeds")
        frontier.append(len(words) - 1)
    while frontier:
        following: list[int] = []
        for parent in frontier:
            for generators in actions:
                candidate = tuple(
                    (generator @ value - value @ generator) % VERIFY_PRIME
                    for generator, value in zip(generators, words[parent], strict=True)
                )
                if insert(candidate, depths[parent] + 1):
                    following.append(len(words) - 1)
        frontier = following
        if len(words) % 64 < len(following):
            guard(started, f"clean K_{{{left_size},{right_size}}} rank {len(words)}")
    return {
        "rank": len(words),
        "maximum_depth": max(depths),
        "coordinates": sum(matrix.size for matrix in field),
        "blocks": blocks,
        "saturated": not frontier,
    }


def so_size(dimension: int) -> int:
    return dimension * (dimension - 1) // 2


def symplectic_size(vector_dimension: int) -> int:
    if vector_dimension & 1:
        raise AssertionError(vector_dimension)
    return vector_dimension * (vector_dimension + 1) // 2


def clean_sector_bound(left_label: int, right_label: int) -> tuple[int, int]:
    total_dimension = (left_label + 1) * (right_label + 1)
    if (left_label + right_label) & 1:
        half = total_dimension // 2
        return half * half, 1
    if left_label != right_label:
        if left_label & 1:
            half = total_dimension // 2
            return 2 * symplectic_size(half), 0
        plus = (total_dimension + 1) // 2
        minus = (total_dimension - 1) // 2
        return so_size(plus) + so_size(minus), 0
    one = left_label + 1
    if left_label & 1:
        a = one * (one + 2) // 4
        b = one * (one - 2) // 4
        c = one * one // 4
        return symplectic_size(a) + symplectic_size(b) + c * c, 1
    a = ((one + 1) // 2) ** 2
    b = ((one - 1) // 2) ** 2
    c = (one * one - 1) // 4
    return so_size(a) + so_size(b) + c * c, 1


def clean_upper(left_size: int, right_size: int) -> int:
    if (left_size - right_size) & 1:
        semisimple = 0
        for r in labels_for_qubits(left_size):
            for s in labels_for_qubits(right_size):
                if r and s:
                    bound, center = clean_sector_bound(r, s)
                    semisimple += bound - center
        return semisimple + 1
    if left_size != right_size:
        raise ValueError((left_size, right_size))
    positive = tuple(label for label in labels_for_qubits(left_size) if label)
    semisimple = 0
    for index, r in enumerate(positive):
        for s in positive[index:]:
            bound, center = clean_sector_bound(r, s)
            semisimple += bound - center
    return semisimple + (2 if left_size % 2 == 0 else 1)


def clean_local_dimension(left_size: int, right_size: int) -> int:
    vertices = left_size + right_size
    if vertices & 1:
        return 2 ** (2 * vertices - 2) - 1
    return 2 ** (2 * vertices - 2) - (-1) ** (vertices // 2) * 2 ** (vertices - 1)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact top-level shape", set(artifact) == {"meta", "data", "checks"})
    check(
        "stored checks pass",
        artifact["checks"] and all(row["passed"] for row in artifact["checks"]),
        len(artifact["checks"]),
    )
    source_hashes = artifact["meta"]["source_sha256"]
    rebuilt_hashes = {path: file_sha256(ROOT / path) for path in SOURCE_PATHS}
    check("artifact source hashes", source_hashes == rebuilt_hashes)
    check(
        "Schur square-sum identity",
        all(
            sum((label + 1) ** 2 for label in labels_for_qubits(size))
            == comb(size + 3, 3)
            for size in range(1, 9)
        ),
        "m=1..8",
    )
    check(
        "verifier modulus prime",
        all(VERIFY_PRIME % divisor for divisor in range(2, math.isqrt(VERIFY_PRIME) + 1)),
        VERIFY_PRIME,
    )

    decisive: dict[str, dict[str, object]] = {}
    for left_size, right_size in ((4, 4), (4, 5)):
        graph = f"K_{{{left_size},{right_size}}}"
        closure = clean_modular_rank(left_size, right_size, started)
        upper = clean_upper(left_size, right_size)
        decisive[graph] = {"closure": closure, "upper": upper}
        check(
            f"{graph} clean modular lower meets clean upper",
            closure["saturated"] and closure["rank"] == upper,
            {"rank": closure["rank"], "upper": upper, "prime": VERIFY_PRIME},
        )

    expected = {"K_{4,4}": 137, "K_{4,5}": 471}
    check(
        "decisive exact dimensions",
        {graph: row["upper"] for graph, row in decisive.items()} == expected,
        expected,
    )
    artifact_cases = {
        row["graph"]: row for row in artifact["data"]["complete_bipartite_cases"]
    }
    check(
        "artifact decisive rows match clean rebuild",
        all(
            artifact_cases[graph]["dimension_Q"] == row["upper"]
            and artifact_cases[graph]["structural_upper"]["upper_dimension"] == row["upper"]
            and all(rank == row["upper"] for rank in artifact_cases[graph]["modular_ranks"])
            for graph, row in decisive.items()
        ),
    )
    check(
        "local-term contrast rebuilt",
        all(
            artifact_cases[graph]["local_term_dimension"]
            == clean_local_dimension(*artifact_cases[graph]["part_sizes"])
            > row["upper"]
            for graph, row in decisive.items()
        ),
        {
            graph: [
                row["upper"],
                clean_local_dimension(*artifact_cases[graph]["part_sizes"]),
            ]
            for graph, row in decisive.items()
        },
    )
    check(
        "new anchors clear quadratic ceiling",
        all(
            row["upper"]
            > sum(artifact_cases[graph]["part_sizes"])
            * (2 * sum(artifact_cases[graph]["part_sizes"]) - 1)
            for graph, row in decisive.items()
        ),
    )
    check(
        "opposite-parity closed-form controls",
        clean_upper(2, 3) == 44
        and clean_upper(3, 4) == 167
        and clean_upper(4, 5) == 471,
        "44,167,471",
    )
    check(
        "symmetry-container threshold limitation",
        all(
            ((m + 1) * (m + 2) // 2) ** 2
            > (2 * m + 1) * (4 * m + 1)
            for m in range(3, 65)
        ),
        "largest block already above the ceiling for audited symbolic range m=3..64",
    )
    guard(started, "verifier completion")
    check("verifier RSS cap", peak_rss_bytes() < RSS_CAP_BYTES, peak_rss_bytes())
    print(
        f"PASS test_kmn_twogen: cpu={time.process_time() - started:.6f}s, "
        f"rss={peak_rss_bytes()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
