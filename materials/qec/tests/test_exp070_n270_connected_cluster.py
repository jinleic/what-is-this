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


def test_two_root_reduction_covers_both_qubit_blocks() -> None:
    for key in E70.TARGETS:
        reduction = E70.root_reduction(key)

        assert reduction["translation_group_size"] == 135
        assert reduction["presentations"] == ["original", "block_swapped"]
        assert reduction["root_coordinate"] == 0
        assert reduction["valid"] is True


def test_published_transport_is_part_of_target_identity() -> None:
    for key in E70.TARGETS:
        identity = E70.certificate_identity(key)
        assert identity["published_transport"] == E70.E69.published_transport(key)
        assert identity["target"]["expected_d"] == 20
