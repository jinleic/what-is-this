#!/usr/bin/env python3
"""Exact finite-memory SAW automata and Collatz--Wielandt certificates.

A memory-k walk on Z^d may revisit a vertex only after more than k steps.
Every self-avoiding walk is therefore a memory-k walk.  Following Poenitz--
Tittmann, a state retains the longest suffix that can still be completed to a
loop of length at most k.  Lattice symmetries identify equivalent suffixes.

The transition matrix is never materialised densely.  A positive integer
vector w satisfying q*A*w < p*w componentwise is an exact certificate that
rho(A) < p/q and hence that the SAW connective constant is below p/q.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import resource
import sys
import time
from collections import defaultdict, deque
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import mpmath as mp
import numpy as np
from scipy.sparse import csr_matrix, identity

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "bounds" / "finite_memory_saw.json"
DIMENSION = 3
DIRECTIONS = tuple(
    tuple(sign if coordinate == axis else 0 for coordinate in range(DIMENSION))
    for axis in range(DIMENSION)
    for sign in (1, -1)
)
SYMMETRIES = tuple(
    (permutation, signs)
    for permutation in itertools.permutations(range(DIMENSION))
    for signs in itertools.product((-1, 1), repeat=DIMENSION)
)
CPU_BUDGET_SECONDS = 1_800.0
RSS_CAP_BYTES = 2_000_000_000
CONTROL_TABLE = {
    4: (3, 4.864536512317585),
    6: (20, 4.807410857993084),
    8: (205, 4.7779698931188275),
    10: (2722, 4.759837615343451),
}

Point = tuple[int, ...]
State = tuple[Point, ...]
Row = tuple[tuple[int, int], ...]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    cpu = time.process_time() - started
    rss = max_rss_bytes()
    if cpu > CPU_BUDGET_SECONDS:
        raise TimeoutError(f"PREFLIGHT/RESOURCE: {stage}: process CPU {cpu:.1f}s")
    if rss > RSS_CAP_BYTES:
        raise MemoryError(f"PREFLIGHT/RESOURCE: {stage}: RSS {rss} bytes")


def add(left: Point, right: Point) -> Point:
    return tuple(a + b for a, b in zip(left, right, strict=True))


@lru_cache(maxsize=None)
def canonical(path: State) -> State:
    base = path[0]
    relative = tuple(
        tuple(point[axis] - base[axis] for axis in range(DIMENSION)) for point in path
    )
    return min(
        tuple(
            tuple(signs[out] * point[permutation[out]] for out in range(DIMENSION))
            for point in relative
        )
        for permutation, signs in SYMMETRIES
    )


@lru_cache(maxsize=None)
def closable(path: State, memory: int) -> bool:
    """Whether `path` can be completed to a simple loop of length <= memory.

    The completion may touch either endpoint but no interior vertex of the
    retained path.  A bounded BFS is exact because only paths of length at
    most `memory - (len(path)-1)` matter.
    """

    length = len(path) - 1
    budget = memory - length
    if budget < 1:
        return False
    start = path[0]
    endpoint = path[-1]
    blocked = frozenset(path[1:-1])
    queue: deque[tuple[Point, int]] = deque([(endpoint, 0)])
    seen = {endpoint}
    while queue:
        point, distance = queue.popleft()
        if distance >= budget:
            continue
        for direction in DIRECTIONS:
            candidate = add(point, direction)
            next_distance = distance + 1
            if candidate == start:
                return True
            if candidate in blocked or candidate in seen:
                continue
            remaining_l1 = sum(abs(candidate[axis] - start[axis]) for axis in range(DIMENSION))
            if next_distance + remaining_l1 > budget:
                continue
            seen.add(candidate)
            queue.append((candidate, next_distance))
    return False


def reduce_state(path: State, memory: int) -> State:
    """Drop exactly the prefix that can no longer participate in a short loop."""

    current = path
    while not closable(current, memory):
        current = current[1:]
    return canonical(current)


def build_automaton(memory: int, *, started: float | None = None) -> tuple[list[State], list[Row]]:
    if memory < 2 or memory % 2:
        raise ValueError("memory must be a positive even integer")
    began = time.process_time() if started is None else started
    origin = (0,) * DIMENSION
    seed = canonical((origin, add(origin, DIRECTIONS[0])))
    states = [seed]
    index = {seed: 0}
    rows: list[Row] = []
    queue: deque[State] = deque([seed])
    while queue:
        state = queue.popleft()
        targets: dict[int, int] = defaultdict(int)
        endpoint = state[-1]
        for direction in DIRECTIONS:
            candidate = add(endpoint, direction)
            if candidate in state:
                # This closes a loop whose length is <= memory.
                continue
            target = reduce_state(state + (candidate,), memory)
            target_index = index.get(target)
            if target_index is None:
                target_index = len(states)
                index[target] = target_index
                states.append(target)
                queue.append(target)
            targets[target_index] += 1
        rows.append(tuple(sorted(targets.items())))
        if len(rows) % 4096 == 0:
            budget_tick(began, f"memory {memory}, {len(rows)}/{len(states)} states")
    if len(rows) != len(states):
        raise AssertionError("transition rows and states disagree")
    return states, rows


def sparse_matrix(rows: list[Row]) -> csr_matrix:
    indptr = [0]
    indices: list[int] = []
    data: list[int] = []
    for row in rows:
        for target, multiplicity in row:
            indices.append(target)
            data.append(multiplicity)
        indptr.append(len(indices))
    return csr_matrix(
        (
            np.asarray(data, dtype=np.float64),
            np.asarray(indices, dtype=np.int32),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=(len(rows), len(rows)),
    )


def apply(rows: list[Row], vector: np.ndarray) -> np.ndarray:
    output = np.empty(len(rows), dtype=np.float64)
    for source, row in enumerate(rows):
        output[source] = math.fsum(multiplicity * vector[target] for target, multiplicity in row)
    return output


def perron_approximation(rows: list[Row]) -> tuple[float, float, np.ndarray, int]:
    """Collatz bracket for rho(A) via power iteration on B = A + I.

    B has the same Perron vector as A with root rho+1, and B has no zero rows
    (dead-end suffix states give zero rows in A), so iterates stay strictly
    positive.  Nothing here is proof; the integer replay below is.
    """

    matrix = sparse_matrix(rows) + identity(len(rows), format="csr")
    live = np.array([bool(row) for row in rows])
    if not live.any():
        raise AssertionError("automaton has no live states")
    vector = np.ones(len(rows), dtype=np.float64)
    lower, upper = 0.0, math.inf
    for iteration in range(1, 200_001):
        image = matrix @ vector
        image /= float(np.max(image))
        vector = image
        if iteration % 25 == 0:
            masked = live & (vector > 0)
            ratios = (matrix @ vector)[masked] / vector[masked]
            lower = float(np.min(ratios)) - 1.0
            upper = float(np.max(ratios)) - 1.0
            if upper - lower < 1e-12:
                return lower, upper, vector, iteration
    return lower, upper, vector, 200_000


def exact_certificate(
    rows: list[Row], numerator: int, denominator: int, approximation: np.ndarray
) -> tuple[list[int], dict[str, int]]:
    """Convert a strict floating Collatz vector into a checked integer one.

    Dead-end states carry Perron weight zero, so their float components decay
    to 0; floor them at half the smallest positive component.  Their rows have
    empty right-hand sides, so the floor only adds negligible mass to their
    predecessors, and the final inequality is replayed in exact integers.
    """

    positive = approximation[approximation > 0]
    if positive.size == 0:
        raise ArithmeticError("Collatz proposal vanished")
    approximation = np.maximum(approximation, float(np.min(positive)) / 2.0)
    image = apply(rows, approximation)
    slack = numerator * approximation - denominator * image
    minimum_slack = float(np.min(slack))
    if not minimum_slack > 0:
        raise ArithmeticError(
            f"approximation is not below {numerator}/{denominator}; min slack={minimum_slack}"
        )
    maximum_row_sum = max(sum(multiplicity for _, multiplicity in row) for row in rows)
    rounding_error = (numerator + denominator * maximum_row_sum) / 2
    scale = max(2, math.ceil(4 * rounding_error / minimum_slack))
    while True:
        weights = [max(1, int(round(scale * float(value)))) for value in approximation]
        residuals = [
            numerator * weights[source]
            - denominator
            * sum(multiplicity * weights[target] for target, multiplicity in row)
            for source, row in enumerate(rows)
        ]
        minimum = min(residuals)
        if minimum > 0:
            return weights, {
                "numerator": numerator,
                "denominator": denominator,
                "scale": scale,
                "minimum_residual": minimum,
                "maximum_residual": max(residuals),
                "maximum_weight": max(weights),
                "minimum_weight": min(weights),
            }
        scale *= 2
        if scale > 10**30:
            raise ArithmeticError("failed to integerize the Collatz certificate")


def state_digest(states: Iterable[State]) -> str:
    payload = ";".join(
        "/".join(",".join(str(coordinate) for coordinate in point) for point in state)
        for state in states
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def transition_digest(rows: Iterable[Row]) -> str:
    payload = ";".join(
        ",".join(f"{target}:{multiplicity}" for target, multiplicity in row) for row in rows
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def choose_decimal_upper(rho: float, digits: int = 6) -> tuple[int, int]:
    denominator = 10**digits
    numerator = math.ceil((rho + 5e-12) * denominator)
    return numerator, denominator


def run(memories: tuple[int, ...]) -> dict[str, object]:
    started = time.process_time()
    rows_out: list[dict[str, object]] = []
    target_payload: dict[str, object] | None = None
    for memory in memories:
        states, rows = build_automaton(memory, started=started)
        bracket_low, bracket_high, vector, iterations = perron_approximation(rows)
        numerator, denominator = choose_decimal_upper(bracket_high)
        weights, certificate = exact_certificate(rows, numerator, denominator, vector)
        exact_ok = all(
            numerator * weights[source]
            > denominator
            * sum(multiplicity * weights[target] for target, multiplicity in row)
            for source, row in enumerate(rows)
        )
        if not exact_ok:
            raise AssertionError("stored Collatz certificate failed exact replay")
        row = {
            "memory": memory,
            "state_count": len(states),
            "transition_count": sum(len(row) for row in rows),
            "weighted_transition_count": sum(
                sum(multiplicity for _, multiplicity in row) for row in rows
            ),
            "state_sha256": state_digest(states),
            "transition_sha256": transition_digest(rows),
            "perron_float_approximation": bracket_high,
            "perron_collatz_float_bracket": [bracket_low, bracket_high],
            "power_iterations": iterations,
            "exact_upper": certificate,
            "upper_decimal": numerator / denominator,
        }
        rows_out.append(row)
        target_payload = {
            "memory": memory,
            "state_sha256": row["state_sha256"],
            "transition_sha256": row["transition_sha256"],
            "weights": weights,
            "exact_upper": certificate,
        }
        budget_tick(started, f"completed memory {memory}")

    if target_payload is None:
        raise ValueError("no memories requested")
    numerator = int(target_payload["exact_upper"]["numerator"])
    denominator = int(target_payload["exact_upper"]["denominator"])
    mp.mp.dps = 100
    kc_value = mp.atanh(mp.mpf(denominator) / numerator)
    kc_floor_digits = 60
    kc_floor = mp.floor(kc_value * mp.mpf(10) ** kc_floor_digits) / mp.mpf(10) ** kc_floor_digits
    incumbent = mp.mpf("0.2122159753270231627267517174278577806020")

    checks: list[dict[str, object]] = []
    for row in rows_out:
        memory = int(row["memory"])
        if memory in CONTROL_TABLE:
            expected_states, expected_rho = CONTROL_TABLE[memory]
            checks.append(
                {
                    "name": f"memory_{memory}_published_control",
                    "passed": row["state_count"] == expected_states
                    and abs(row["perron_float_approximation"] - expected_rho) < 2e-12,
                    "detail": f"states={row['state_count']}, rho~{row['perron_float_approximation']:.15f}",
                }
            )
    checks.extend(
        [
            {
                "name": "target_exact_collatz_strict",
                "passed": int(target_payload["exact_upper"]["minimum_residual"]) > 0,
                "detail": f"minimum integer residual {target_payload['exact_upper']['minimum_residual']}",
            },
            {
                "name": "ising_lower_endpoint_strictly_improves_incumbent",
                "passed": kc_floor > incumbent,
                "detail": f"new floor {mp.nstr(kc_floor, 70)} > {mp.nstr(incumbent, 50)}",
            },
        ]
    )
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError("one or more producer checks failed")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e227_finite_memory_saw.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integer transitions and Collatz certificate; float64 only proposes the positive vector; 100-dps mpmath only reports the directed decimal floor",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
        },
        "data": {
            "dimension": DIMENSION,
            "memories": rows_out,
            "target": target_payload,
            "theorem": {
                "tag": "[THEOREM][COMPUTATION][EXTERNAL]",
                "statement": f"Every simple-cubic SAW is accepted by the memory-{target_payload['memory']} automaton, whose exact integer Collatz certificate proves mu < {numerator}/{denominator}.",
                "ising_consequence": f"K_c > atanh({denominator}/{numerator}) > {mp.nstr(kc_floor, kc_floor_digits + 2)}.",
                "directed_floor": mp.nstr(kc_floor, kc_floor_digits + 2),
                "mpmath_dps": mp.mp.dps,
                "benchmark_used_for_selection": False,
                "source": {
                    "citation": "A. Poenitz and P. Tittmann, Improved Upper Bounds for Self-Avoiding Walks in Z^d, EJC 7 (2000) R21",
                    "doi": "10.37236/1499",
                    "exact_locations": ["Section 2 automaton", "Section 3 upper bounds", "Table 2"],
                },
            },
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-memory", type=int, default=12, choices=(4, 6, 8, 10, 12, 14))
    args = parser.parse_args()
    memories = tuple(memory for memory in (4, 6, 8, 10, 12, 14) if memory <= args.max_memory)
    payload = run(memories)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in payload["data"]["memories"]:
        exact = row["exact_upper"]
        print(
            f"[PASS] k={row['memory']:2d}: states={row['state_count']:7d}, "
            f"rho < {exact['numerator']}/{exact['denominator']} = {row['upper_decimal']:.7f}, "
            f"min residual={exact['minimum_residual']}"
        )
    theorem = payload["data"]["theorem"]
    print("[PASS]", theorem["statement"])
    print("[PASS]", theorem["ising_consequence"])
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
