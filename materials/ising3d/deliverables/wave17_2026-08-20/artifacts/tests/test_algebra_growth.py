"""Standalone acceptance checks for the depth-graded Onsager-algebra experiment."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.e20_algebra_growth import PRIME_31, compute_profile


RESULT_DIR = ROOT / "results" / "algebra_growth"


def _assert_profile(record: dict) -> None:
    dimensions = record["dimensions"]
    assert dimensions, record
    assert dimensions[0] == 2, record
    assert all(a <= b for a, b in zip(dimensions, dimensions[1:])), record
    equal_steps = [i for i, (a, b) in enumerate(zip(dimensions, dimensions[1:]), start=2) if a == b]
    assert not equal_steps or (record["saturated"] and equal_steps == [len(dimensions)]), record
    if record["incomplete_depth"] is None:
        assert record["lower_bound"] == dimensions[-1], record
    else:
        assert record["lower_bound"] >= dimensions[-1], record
    assert record["prime"] == PRIME_31, record
    if record["saturated"]:
        assert record["exact_dimension"] == dimensions[-1], record


def main() -> None:
    # Recompute the finite Onsager controls rather than trusting stored data.
    for n in range(2, 9):
        record = compute_profile(1, n, max_depth=None)
        _assert_profile(record)
        assert record["saturated"] and record["exact_dimension"] == n * n, record
    for n in range(3, 9):
        record = compute_profile(1, n, periodic_cols=True, max_depth=None)
        _assert_profile(record)
        assert record["saturated"] and record["exact_dimension"] == 3 * n - 1, record

    # Reproduce every previously known exact two-generator dimension with the new engine.
    expected = {(2, 2): 11, (2, 3): 263, (2, 4): 2952}
    for shape, dimension in expected.items():
        record = compute_profile(*shape, max_depth=None)
        _assert_profile(record)
        assert record["saturated"] and record["exact_dimension"] == dimension, record

    dense_prefix = compute_profile(2, 3, max_depth=10, engine="dense")
    sparse_prefix = compute_profile(2, 3, max_depth=10, engine="sparse")
    assert dense_prefix["dimensions"] == sparse_prefix["dimensions"]
    assert sum(dense_prefix["grading_block_columns"]) == dense_prefix["support_orbits"]

    with (RESULT_DIR / "profiles.json").open() as handle:
        profiles_artifact = json.load(handle)
    with (RESULT_DIR / "exact_dimensions.json").open() as handle:
        exact_artifact = json.load(handle)

    assert profiles_artifact["provenance"]["script"] == "experiments/e20_algebra_growth.py"
    requested = {
        *(f"chain_open_{n}" for n in range(2, 33)),
        *(f"ring_{n}" for n in range(3, 33)),
        *(f"grid_{r}x{c}" for r, c in ((2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 3), (3, 4))),
    }
    stored = {row["case"] for row in profiles_artifact["data"]["profiles"]}
    assert requested <= stored, requested - stored
    for row in profiles_artifact["data"]["profiles"]:
        _assert_profile(row)

    dimensions = {
        row["case"]: row["exact_dimension"]
        for row in exact_artifact["data"]["dimensions"]
        if row["exact_dimension"] is not None
    }
    assert dimensions["grid_2x2"] == 11
    assert dimensions["grid_2x3"] == 263
    assert dimensions["grid_2x4"] == 2952
    assert dimensions["grid_3x3"] == 8034
    support = {
        row["case"]: row["support_orbits"]
        for row in exact_artifact["data"]["support_statistics"]
        if row["complete"]
    }
    assert support["grid_2x5"] == 66_304
    assert support["grid_2x6"] == 1_049_088
    assert support["grid_3x3"] == 8_739
    assert support["grid_3x4"] == 1_052_672
    certificate = profiles_artifact["data"]["finite_ladder_certificate"]
    assert [row["length"] for row in certificate] == [2, 3, 4, 5, 6]
    assert all(row["passed"] for row in certificate)
    assert all(check["passed"] for check in profiles_artifact["checks"])
    assert all(check["passed"] for check in exact_artifact["checks"])
    print("PASS")


if __name__ == "__main__":
    main()
