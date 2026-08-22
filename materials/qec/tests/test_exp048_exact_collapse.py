"""Regression checks for the exact EXP-048 collapse certificate."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp048_exact_collapse_sweep", ROOT / "experiments" / "exp048_exact_collapse_sweep.py"
)
assert SPEC and SPEC.loader
EXP048 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP048)


def _clone_record(record: dict, **changes):
    import copy

    cloned = copy.deepcopy(record)
    cloned.update(changes)
    return cloned


def _item(cohort, index=37):
    return cohort[index]


def test_flagship_sibling_is_exact_distance_six() -> None:
    _, cohort = EXP048.load_cohort()
    record = EXP048.certify_row(_item(cohort))  # 12_6_0193
    assert record["label"] == "12_6_0193"
    assert record["k_parent"] == record["k_q"] == 12
    assert record["x0_witness"]["verified_logical"]
    assert record["exact_full_distance"] == 6
    assert record["exact_x_distance"] == 6
    assert record["strict_parent_drop_certified"]
    assert record["calls"][-1]["status"] == "UNSAT"


def test_identity_gate_rejects_tampered_certificate() -> None:
    _, cohort = EXP048.load_cohort()
    record = EXP048.certify_row(_item(cohort))
    assert EXP048.expected_identity(_item(cohort))["parent_fingerprint"] == record[
        "parent_fingerprint"
    ]
    EXP048.assert_record_identity(_item(cohort), record)
    forged = _clone_record(record, parent_fingerprint="0" * 64)
    import pytest

    with pytest.raises(RuntimeError, match="row identity mismatch"):
        EXP048.assert_record_identity(_item(cohort), forged)


def test_fixed_cap_protocol_is_enforced() -> None:
    import pytest

    _, cohort = EXP048.load_cohort()
    with pytest.raises(ValueError, match="fixed cap=5"):
        EXP048.certify_row(_item(cohort), cap=4)
