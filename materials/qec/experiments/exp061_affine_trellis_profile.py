"""EXP-061: exact affine-trellis width preflight for N=105 residuals.

For an affine sector Mv=b and a fixed column order, the exact conventional
binary-trellis width at cut i is

  w_i = rank(M[:, prefix_i]) + rank(M[:, suffix_i]) - rank(M).

The number of completable syndromes is exactly 2^w_i.  This preflight computes
that profile before implementing any dynamic program.  If every deterministic
CRT/grid order is too wide, the method fails closed without a distance claim.

Run:
  python experiments/exp061_affine_trellis_profile.py run
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp061", "exp056_odd_distance.py")
E59 = _load("exp059_for_exp061", "exp059_coprime_transport.py")

from qec_research.gf2.linalg import rank_bitset, rank_np, rows_to_bitsets  # noqa: E402

SCHEMA = "exp061-affine-trellis-profile-v1"
SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "15x7.json"
OUT = ROOT / "results" / "processed" / "exp061_affine_trellis_profile.json"
STATE_LIMIT = 5_000_000


def problem_of(record: dict[str, Any]) -> dict[str, Any]:
    return E56.build_problem(
        {
            "ell": record["ell"],
            "m": record["m"],
            "expected_k": record["k_parent"],
            "A": record["A"],
            "B": record["B"],
        }
    )


def affine_sector_matrices(problem: dict[str, Any]) -> list[np.ndarray]:
    cover = E56.build_orbit_cover(problem)
    basis = E56._independent_rows(problem["HX"])
    matrices = []
    for sector in cover["sectors"]:
        matrix = E56._independent_rows(
            np.vstack([basis, sector["detector"]])
        )
        if matrix.shape[0] != basis.shape[0] + 1:
            raise RuntimeError("detector did not define an affine functional")
        matrices.append(matrix)
    return matrices


def _source_order_for_grid(target: tuple[int, int]) -> np.ndarray:
    mapping = E59.coordinate_transport(15, 7, *target)
    return np.argsort(mapping)


def order_family() -> dict[str, np.ndarray]:
    N = 105
    pi = np.asarray([(e % 15) * 7 + (e % 7) for e in range(N)], dtype=np.int64)
    orders: dict[str, np.ndarray] = {}
    bases = {
        "grid15x7": np.arange(N, dtype=np.int64),
        "grid21x5": _source_order_for_grid((21, 5)),
        "grid35x3": _source_order_for_grid((35, 3)),
        "pi": pi,
    }
    for name, base in bases.items():
        orders[f"{name}_blocks"] = np.concatenate([base, N + base])
        orders[f"{name}_interleaved"] = np.column_stack(
            [base, N + base]
        ).reshape(-1)
        reverse = base[::-1]
        orders[f"{name}_reverse_interleaved"] = np.column_stack(
            [reverse, N + reverse]
        ).reshape(-1)
    for name, order in orders.items():
        if not np.array_equal(np.sort(order), np.arange(2 * N)):
            raise RuntimeError(f"{name} is not a coordinate permutation")
    return orders


def _column_rank_numpy(matrix: np.ndarray, columns: np.ndarray) -> int:
    return rank_np(matrix[:, columns]) if len(columns) else 0


def _column_rank_bitset(matrix: np.ndarray, columns: np.ndarray) -> int:
    if not len(columns):
        return 0
    vectors = matrix[:, columns].T
    return rank_bitset(rows_to_bitsets(vectors), matrix.shape[0])


def rank_profile(
    matrix: np.ndarray, order: np.ndarray, *, bitset_crosscheck: bool
) -> dict[str, Any]:
    q, n = matrix.shape
    total_rank = rank_np(matrix)
    if total_rank != q:
        raise RuntimeError("affine matrix is not row independent")
    widths = []
    for cut in range(n + 1):
        prefix, suffix = order[:cut], order[cut:]
        left = _column_rank_numpy(matrix, prefix)
        right = _column_rank_numpy(matrix, suffix)
        if bitset_crosscheck:
            if left != _column_rank_bitset(matrix, prefix):
                raise RuntimeError("prefix rank paths disagree")
            if right != _column_rank_bitset(matrix, suffix):
                raise RuntimeError("suffix rank paths disagree")
        widths.append(left + right - q)
    maximum = max(widths)
    return {
        "q": q,
        "n": n,
        "max_width": maximum,
        "peak_unpruned_states": 1 << maximum,
        "within_state_limit": (1 << maximum) <= STATE_LIMIT,
        "widths": widths,
        "bitset_crosschecked": bitset_crosscheck,
    }


def selected_targets(shard: dict[str, Any]) -> list[tuple[str, int, dict[str, Any]]]:
    # One representative for each hard k plus the no-reference k=24 class.
    return [
        ("undecided", 0, shard["undecided"][0]),  # k=14
        ("undecided", 1, shard["undecided"][1]),  # k=10
        ("undecided", 4, shard["undecided"][4]),  # k=8
        ("no_reference", 0, shard["no_reference"][0]),  # k=24
    ]


def run() -> int:
    shard = json.loads(SHARD.read_text(encoding="utf-8"))
    orders = order_family()
    outputs = []
    for group, index, record in selected_targets(shard):
        problem = problem_of(record)
        sectors = affine_sector_matrices(problem)
        per_order = []
        for name, order in orders.items():
            profiles = [rank_profile(matrix, order, bitset_crosscheck=False)
                        for matrix in sectors]
            per_order.append(
                {
                    "name": name,
                    "max_width": max(p["max_width"] for p in profiles),
                    "sectors": profiles,
                }
            )
        per_order.sort(key=lambda item: (item["max_width"], item["name"]))
        best = per_order[0]
        # The selected performance claim is independently rank-checked.
        best_profiles = [
            rank_profile(matrix, orders[best["name"]], bitset_crosscheck=True)
            for matrix in sectors
        ]
        best = {
            "name": best["name"],
            "max_width": max(p["max_width"] for p in best_profiles),
            "sectors": best_profiles,
        }
        out = {
            "group": group,
            "index": index,
            "A": record["A"],
            "B": record["B"],
            "n": problem["n"],
            "k": problem["k"],
            "threshold": record["threshold"],
            "best": best,
            "orders": per_order,
            "viable": all(
                profile["within_state_limit"] for profile in best_profiles
            ),
        }
        outputs.append(out)
        print(json.dumps({
            "group": group, "index": index, "k": problem["k"],
            "best_order": best["name"], "max_width": best["max_width"],
            "viable": out["viable"],
        }), flush=True)
    payload = {
        "schema": SCHEMA,
        "state_limit": STATE_LIMIT,
        "targets": outputs,
        "viable_targets": sum(item["viable"] for item in outputs),
        "falsified_targets": sum(not item["viable"] for item in outputs),
        "interpretation": (
            "Only a complete future trellis recurrence can certify distance. "
            "This artifact certifies exact order widths and fails the method "
            "closed when peak unpruned states exceed the state limit."
        ),
    }
    E56.atomic_write_json(OUT, payload)
    print(json.dumps({
        "viable_targets": payload["viable_targets"],
        "falsified_targets": payload["falsified_targets"],
    }, indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run",))
    parser.parse_args()
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
