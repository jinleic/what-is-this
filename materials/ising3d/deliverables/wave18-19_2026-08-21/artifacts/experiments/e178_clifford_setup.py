"""Exact Jordan--Wigner conventions and finite Pauli-closure primitives.

This module is the convention layer for the wave-18 Clifford-grade audit.  It
uses only integer/finite combinatorics.  In particular, Pauli phases are stored
as powers of ``i`` and are never inferred numerically.
"""

from __future__ import annotations

import json
import platform
import resource
import time
from collections import deque
from dataclasses import dataclass
from math import comb
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
CPU_BUDGET_SECONDS = 240.0
RSS_CAP_BYTES = 4_000_000_000

# A word is i**phase times a tensor product of the Hermitian I/X/Y/Z matrices.
PauliWord = tuple[int, tuple[str, ...]]
_LOCAL_PRODUCT: dict[tuple[str, str], tuple[int, str]] = {
    ("I", "I"): (0, "I"),
    ("I", "X"): (0, "X"),
    ("I", "Y"): (0, "Y"),
    ("I", "Z"): (0, "Z"),
    ("X", "I"): (0, "X"),
    ("Y", "I"): (0, "Y"),
    ("Z", "I"): (0, "Z"),
    ("X", "X"): (0, "I"),
    ("Y", "Y"): (0, "I"),
    ("Z", "Z"): (0, "I"),
    ("X", "Y"): (1, "Z"),
    ("Y", "X"): (3, "Z"),
    ("Y", "Z"): (1, "X"),
    ("Z", "Y"): (3, "X"),
    ("Z", "X"): (1, "Y"),
    ("X", "Z"): (3, "Y"),
}


@dataclass(frozen=True)
class ClosureResult:
    """Exact monomial Lie closure reached by generator-adjoint BFS."""

    dimension: int
    elements: frozenset[int]
    process_time_seconds: float
    peak_rss_bytes: int


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str, budget: float = CPU_BUDGET_SECONDS) -> None:
    used = time.process_time() - started
    if used > budget:
        raise RuntimeError(
            f"NON-DECISIVE process-time expiry at {stage}: {used:.3f}s > {budget:.3f}s"
        )
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise RuntimeError(
            f"NON-DECISIVE RSS expiry at {stage}: {rss} >= {RSS_CAP_BYTES} bytes"
        )


def multiply_words(left: PauliWord, right: PauliWord) -> PauliWord:
    if len(left[1]) != len(right[1]):
        raise ValueError("Pauli words have different site counts")
    phase = (left[0] + right[0]) % 4
    letters: list[str] = []
    for a, b in zip(left[1], right[1], strict=True):
        local_phase, letter = _LOCAL_PRODUCT[(a, b)]
        phase = (phase + local_phase) % 4
        letters.append(letter)
    return phase, tuple(letters)


def scale_word(word: PauliWord, i_power: int) -> PauliWord:
    return (word[0] + i_power) % 4, word[1]


def local_word(n: int, labels: dict[int, str], phase: int = 0) -> PauliWord:
    letters = ["I"] * n
    for site, label in labels.items():
        if not 0 <= site < n or label not in {"X", "Y", "Z"}:
            raise ValueError((site, label))
        letters[site] = label
    return phase % 4, tuple(letters)


def jw_majoranas(n: int) -> tuple[PauliWord, ...]:
    """Return gamma_(2r)=X_<r Z_r and gamma_(2r+1)=X_<r Y_r."""
    if n < 1:
        raise ValueError("n must be positive")
    values: list[PauliWord] = []
    for site in range(n):
        prefix = {prior: "X" for prior in range(site)}
        values.append(local_word(n, prefix | {site: "Z"}))
        values.append(local_word(n, prefix | {site: "Y"}))
    return tuple(values)


def majorana_product(majoranas: Sequence[PauliWord], support: Iterable[int]) -> PauliWord:
    answer: PauliWord = (0, tuple("I" for _ in range(len(majoranas) // 2)))
    for index in sorted(support):
        answer = multiply_words(answer, majoranas[index])
    return answer


def jw_edge_support(p: int, q: int) -> tuple[int, ...]:
    """Majorana support of Z_p Z_q for zero-based path positions p<q."""
    if not 0 <= p < q:
        raise ValueError((p, q))
    return tuple(range(2 * p + 1, 2 * q + 1))


def verify_jw_conventions(max_n: int = 9) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []
    for n in range(1, max_n + 1):
        gammas = jw_majoranas(n)
        identity = local_word(n, {})
        squares = all(multiply_words(gamma, gamma) == identity for gamma in gammas)
        anticommute = True
        for a in range(2 * n):
            for b in range(a + 1, 2 * n):
                ab = multiply_words(gammas[a], gammas[b])
                ba = multiply_words(gammas[b], gammas[a])
                if ab != scale_word(ba, 2):
                    anticommute = False
        checks.append(
            {
                "name": f"JW Clifford relations n={n}",
                "passed": squares and anticommute,
                "detail": "gamma_a^2=1 and gamma_a gamma_b=-gamma_b gamma_a",
            }
        )
        for site in range(n):
            image = scale_word(
                multiply_words(gammas[2 * site], gammas[2 * site + 1]), 1
            )
            checks.append(
                {
                    "name": f"JW field phase n={n} site={site}",
                    "passed": image == local_word(n, {site: "X"}),
                    "detail": "X_r=i gamma_(2r) gamma_(2r+1)",
                }
            )
        for p in range(n):
            for q in range(p + 1, n):
                support = jw_edge_support(p, q)
                image = scale_word(majorana_product(gammas, support), q - p)
                target = local_word(n, {p: "Z", q: "Z"})
                checks.append(
                    {
                        "name": f"JW edge phase n={n} edge={p}-{q}",
                        "passed": image == target and len(support) == 2 * (q - p),
                        "detail": (
                            f"Z_p Z_q=i^{q-p} gamma_[2p+1..2q], grade={len(support)}"
                        ),
                    }
                )
    return checks


def normalize_edges(edges: Iterable[tuple[int, int]], n: int) -> tuple[tuple[int, int], ...]:
    answer: set[tuple[int, int]] = set()
    for u, v in edges:
        if not 0 <= u < n or not 0 <= v < n or u == v:
            raise ValueError((u, v))
        answer.add((u, v) if u < v else (v, u))
    return tuple(sorted(answer))


def grid_edges(rows: int, cols: int) -> tuple[tuple[int, int], ...]:
    if rows < 1 or cols < 1:
        raise ValueError((rows, cols))
    edges: list[tuple[int, int]] = []
    index = lambda row, col: row * cols + col  # noqa: E731
    for row in range(rows):
        for col in range(cols):
            if row + 1 < rows:
                edges.append((index(row, col), index(row + 1, col)))
            if col + 1 < cols:
                edges.append((index(row, col), index(row, col + 1)))
    return normalize_edges(edges, rows * cols)


def snake_path(rows: int, cols: int) -> tuple[int, ...]:
    path: list[int] = []
    for row in range(rows):
        columns = range(cols) if row % 2 == 0 else range(cols - 1, -1, -1)
        path.extend(row * cols + col for col in columns)
    return tuple(path)


def cycle_edges(n: int) -> tuple[tuple[int, int], ...]:
    if n < 3:
        raise ValueError("a simple cycle needs at least three vertices")
    return normalize_edges(
        [(site, site + 1) for site in range(n - 1)] + [(n - 1, 0)], n
    )


def validate_hamiltonian_path(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> None:
    if len(path) != n or set(path) != set(range(n)):
        raise ValueError("path is not a permutation of the vertices")
    edge_set = set(normalize_edges(edges, n))
    for u, v in zip(path, path[1:]):
        edge = (u, v) if u < v else (v, u)
        if edge not in edge_set:
            raise ValueError(f"missing Hamiltonian-path edge {edge}")


def bipartite_coloring(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[bool, tuple[int, ...]]:
    adjacency = [[] for _ in range(n)]
    for u, v in normalize_edges(edges, n):
        adjacency[u].append(v)
        adjacency[v].append(u)
    colors = [-1] * n
    for root in range(n):
        if colors[root] >= 0:
            continue
        colors[root] = 0
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if colors[v] < 0:
                    colors[v] = colors[u] ^ 1
                    queue.append(v)
                elif colors[v] == colors[u]:
                    return False, tuple(colors)
    return True, tuple(colors)


def nonpath_edge_rows(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> tuple[dict[str, object], ...]:
    validate_hamiltonian_path(n, edges, path)
    position = {vertex: place for place, vertex in enumerate(path)}
    path_edges = {
        tuple(sorted((u, v))) for u, v in zip(path, path[1:])
    }
    rows: list[dict[str, object]] = []
    for edge in normalize_edges(edges, n):
        if edge in path_edges:
            continue
        p, q = sorted((position[edge[0]], position[edge[1]]))
        distance = q - p
        rows.append(
            {
                "edge": list(edge),
                "path_positions_zero_based": [p, q],
                "path_distance": distance,
                "majorana_grade": 2 * distance,
                "grade_mod_4": (2 * distance) % 4,
                "support_zero_based": list(jw_edge_support(p, q)),
                "phase_i_power": distance % 4,
            }
        )
    return tuple(rows)


def raw_local_generators(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    fields = [1 << site for site in range(n)]
    bonds = [(1 << (n + u)) | (1 << (n + v)) for u, v in normalize_edges(edges, n)]
    return tuple(fields + bonds)


def symplectic_form(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    x_left, z_left = left & mask, left >> n
    x_right, z_right = right & mask, right >> n
    return ((x_left & z_right).bit_count() + (z_left & x_right).bit_count()) & 1


def exact_pauli_closure(
    n: int,
    edges: Sequence[tuple[int, int]],
    *,
    budget: float = CPU_BUDGET_SECONDS,
) -> ClosureResult:
    """Compute the exact monomial DLA by generator-adjoint (left-normed) BFS."""
    started = time.process_time()
    generators = raw_local_generators(n, edges)
    seen = set(generators)
    queue = deque(generators)
    steps = 0
    while queue:
        current = queue.popleft()
        for generator in generators:
            if symplectic_form(current, generator, n):
                candidate = current ^ generator
                if candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
        steps += 1
        if steps % 16_384 == 0:
            budget_tick(started, f"Pauli closure n={n}, seen={len(seen)}", budget)
    budget_tick(started, f"completed Pauli closure n={n}", budget)
    return ClosureResult(
        dimension=len(seen),
        elements=frozenset(seen),
        process_time_seconds=time.process_time() - started,
        peak_rss_bytes=max_rss_bytes(),
    )


def grade_two_mod_four_dimension(n: int) -> int:
    """Target noncentral grade-class dimension on 2n Majoranas."""
    total = sum(comb(2 * n, grade) for grade in range(2, 2 * n + 1, 4))
    return total - (1 if n % 2 else 0)


def grade_two_mod_four_closed_form(n: int) -> int:
    if n % 2:
        return (1 << (2 * n - 2)) - 1
    return (1 << (2 * n - 2)) - ((-1) ** (n // 2)) * (1 << (n - 1))


def full_even_derived_dimension(n: int) -> int:
    """Dimension of [Cl^even,Cl^even]: omit identity and the volume element."""
    return (1 << (2 * n - 1)) - 2


def run_setup() -> dict[str, object]:
    started = time.process_time()
    checks = verify_jw_conventions()
    for n in range(2, 25):
        direct = grade_two_mod_four_dimension(n)
        closed = grade_two_mod_four_closed_form(n)
        checks.append(
            {
                "name": f"root-of-unity dimension formula n={n}",
                "passed": direct == closed,
                "detail": f"direct={direct}, closed_form={closed}",
            }
        )
    budget_tick(started, "JW setup checks")
    if not all(bool(row["passed"]) for row in checks):
        raise AssertionError("Jordan--Wigner setup check failed")
    return {
        "tag": "[COMPUTATION]",
        "convention": {
            "majoranas": "gamma_(2r)=X_0...X_(r-1) Z_r; gamma_(2r+1)=X_0...X_(r-1) Y_r",
            "field": "X_r=i gamma_(2r) gamma_(2r+1)",
            "edge": "Z_p Z_q=i^(q-p) gamma_(2p+1)...gamma_(2q), p<q",
            "antihermitian_generators": "iX_r and iZ_pZ_q",
        },
        "checks": checks,
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_setup()
    print(json.dumps({"checks": len(result["checks"]), "process_time_seconds": result["process_time_seconds"]}))
    print("PASS e178 Clifford setup")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
