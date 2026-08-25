"""EXP-065: exact H-orbit quotient projection for the hard N=105 k=8 class.

For a pole stabilizer H, partition both physical blocks into H-orbits and map a
word to the parity on each orbit.  Hamming weight cannot increase under this
projection.  Projecting the Z-stabilizer rowspace gives a short binary code;
its complete coset-leader table is computed by a solver-free syndrome DP in
O(L 2^h), where L is the projected length and h its codimension.  Projecting
all 2^k logical pole classes then yields a rigorous physical distance lower
bound (possibly weak, never overstated).

Run:
  python experiments/exp065_quotient_projection.py run
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp065", "exp056_odd_distance.py")

from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402

SCHEMA = "exp065-h-orbit-projection-v1"
OUT = ROOT / "results" / "processed" / "exp065_n105_k8_projection.json"
TARGET: dict[str, Any] = {
    "ell": 15,
    "m": 7,
    "expected_k": 8,
    "A": [[0, 0], [1, 0], [4, 1]],
    "B": [[0, 0], [2, 3], [8, 1]],
}


def coordinate_orbits(problem: dict[str, Any], shifts: list[list[int]]) -> list[list[int]]:
    unseen = set(range(problem["n"]))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = sorted({
            int(E56._permutation_map(problem, a, b, False)[seed])
            for a, b in shifts
        })
        if len(orbit) != len(shifts):
            raise RuntimeError("H action is not free on physical coordinates")
        orbits.append(orbit)
        unseen.difference_update(orbit)
    if sorted(index for orbit in orbits for index in orbit) != list(range(problem["n"])):
        raise RuntimeError("H orbits do not partition physical coordinates")
    return orbits


def projection_matrix(problem: dict[str, Any], orbits: list[list[int]]) -> np.ndarray:
    matrix = np.zeros((len(orbits), problem["n"]), dtype=np.uint8)
    for row, orbit in enumerate(orbits):
        matrix[row, orbit] = 1
    return matrix


def syndrome_leader_table(checks: np.ndarray) -> list[int]:
    h, length = checks.shape
    columns = []
    for column in range(length):
        value = 0
        for row in range(h):
            value |= int(checks[row, column]) << row
        columns.append(value)
    infinity = length + 1
    dp = [infinity] * (1 << h)
    dp[0] = 0
    for column in columns:
        prior = list(dp)
        for syndrome, weight in enumerate(prior):
            if weight == infinity:
                continue
            target = syndrome ^ column
            if weight + 1 < dp[target]:
                dp[target] = weight + 1
    if any(weight == infinity for weight in dp):
        raise RuntimeError("projected syndrome table is incomplete")
    return dp


def syndrome_bits(checks: np.ndarray, vector: np.ndarray) -> int:
    values = (checks.astype(np.int64) @ vector.astype(np.int64) % 2).astype(np.uint8)
    return sum(int(bit) << index for index, bit in enumerate(values))


def run() -> int:
    problem = E56.build_problem(TARGET)
    cover = E56.build_orbit_cover(problem)
    shifts = cover["sectors"][0]["stabilizer_shifts"]
    orbits = coordinate_orbits(problem, shifts)
    projection = projection_matrix(problem, orbits)
    projected_stabilizers = (
        problem["HZ"].astype(np.int64) @ projection.T.astype(np.int64) % 2
    ).astype(np.uint8)
    stabilizer_rank = rank_np(projected_stabilizers)
    checks = nullspace_np(projected_stabilizers)
    if rank_np(checks) != len(orbits) - stabilizer_rank:
        raise RuntimeError("projected parity-check rank is inconsistent")
    leaders = syndrome_leader_table(checks)

    logical_histogram: Counter[int] = Counter()
    logical_syndromes = set()
    zero_projected_classes = 0
    for mask in range(1, 1 << problem["k"]):
        coefficients = np.asarray(
            [(mask >> bit) & 1 for bit in range(problem["k"])], dtype=np.uint8
        )
        representative = (coefficients @ cover["z_pole"]) % 2
        projected = (projection.astype(np.int64) @ representative.astype(np.int64) % 2).astype(np.uint8)
        syndrome = syndrome_bits(checks, projected)
        logical_syndromes.add(syndrome)
        leader = leaders[syndrome]
        logical_histogram[leader] += 1
        zero_projected_classes += int(syndrome == 0)
    lower_bound = min(logical_histogram)
    payload = {
        "schema": SCHEMA,
        "target": TARGET,
        "n": problem["n"],
        "k": problem["k"],
        "cover_complete": cover["complete"],
        "stabilizer_size": len(shifts),
        "projected_length": len(orbits),
        "orbit_sizes": sorted({len(orbit) for orbit in orbits}),
        "projected_stabilizer_rank": stabilizer_rank,
        "projected_codimension": checks.shape[0],
        "logical_classes": (1 << problem["k"]) - 1,
        "distinct_projected_logical_syndromes": len(logical_syndromes),
        "zero_projected_classes": zero_projected_classes,
        "projected_coset_leader_histogram": {
            str(weight): count for weight, count in sorted(logical_histogram.items())
        },
        "rigorous_physical_distance_lower_bound": lower_bound,
        "argument": (
            "For every physical coset vector v, each odd H-orbit contains at least "
            "one support bit, so wt(v) >= wt(projection(v)). The syndrome DP is "
            "complete for every projected stabilizer coset."
        ),
    }
    E56.atomic_write_json(OUT, payload)
    print(json.dumps({k: payload[k] for k in (
        "stabilizer_size", "projected_length", "projected_stabilizer_rank",
        "projected_codimension", "zero_projected_classes",
        "projected_coset_leader_histogram",
        "rigorous_physical_distance_lower_bound",
    )}, indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run",))
    parser.parse_args()
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
