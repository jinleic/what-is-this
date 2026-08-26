"""Regression checks for EXP-070's n=270 connected-cluster certificates."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp070", ROOT / "experiments" / "exp070_n270_connected_cluster.py"
)
E70 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E70
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E70)


def test_n270_targets_have_verified_weight20_upper_bounds() -> None:
    for key in E70.TARGETS:
        problem = E70.problem_of(key)
        upper = E70.upper_certificate(key)

        assert problem["n"] == 270
        assert problem["k"] == 8
        assert upper["weight"] == 20
        assert upper["z_valid_numpy"] is True
        assert upper["z_valid_bitset"] is True
        assert upper["x_valid_numpy"] is True
        assert upper["x_valid_bitset"] is True


def test_prefix_partition_and_pure_block_enumeration_cover_all_logicals() -> None:
    expected = {
        "liang270a": [55, 62, 135, 150, 159],
        "liang270b": [55, 56, 135, 175, 269],
    }
    for key in E70.TARGETS:
        reduction = E70.root_reduction(key)
        partition = E70.prefix_partition(key)
        pure_second = E70.pure_second_logical_enumeration(key)

        assert reduction["translation_group_size"] == 135
        assert reduction["presentations"] == ["original"]
        assert reduction["root_coordinate"] == 0
        assert reduction["valid"] is True
        assert partition["prefix_columns"] == expected[key]
        assert partition["complete"] is True
        assert pure_second["complete"] is True
        assert pure_second["minimum_weight"] >= 20


def test_prefix_solver_build_is_hash_bound() -> None:
    binding = E70.solver_build_binding()
    assert binding["binary_sha256"] == E70.PREFIX_BINARY_SHA256
    assert binding["source_commit"] == E70.E67.DIST_M4RI_COMMIT


def test_published_transport_is_part_of_target_identity() -> None:
    for key in E70.TARGETS:
        identity = E70.certificate_identity(key)
        assert identity["published_transport"] == E70.E69.published_transport(key)
        assert identity["target"]["expected_d"] == 20
