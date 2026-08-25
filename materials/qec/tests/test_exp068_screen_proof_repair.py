"""Regression checks for EXP-068's record-level domination proofs."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp068", ROOT / "experiments" / "exp068_screen_proof_repair.py"
)
E68 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E68
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E68)


def test_all_fallback_evidence_is_physical_and_identity_bound() -> None:
    screen = json.loads(E68.E55.SCREEN_OUT.read_text(encoding="utf-8"))
    records = [
        record
        for shard in screen["lattices"]
        for record in shard["records"]
        if record["verdict"] == "dominated"
    ]
    assert len(records) == 57
    independent = transported = 0
    for record in records:
        path = E68.evidence_path(record)
        payload = json.loads(path.read_text(encoding="utf-8"))
        E68.validate_evidence(payload, record)
        assert record["fallback_witness"]["sha256"] == E68.file_sha256(path)
        assert E68.E55._fallback_witness_evidence_valid(record)
        if payload["status"] == "TRANSPORTED":
            transported += 1
        else:
            independent += 1
    assert (independent, transported) == (29, 28)


def test_independent_recovery_set_has_exactly_twenty_nine_records() -> None:
    items = E68.source_records()
    assert len(items) == 29
    assert all(E68.evidence_path(item["record"]).is_file() for item in items)
    assert {tuple((item["record"]["ell"], item["record"]["m"])) for item in items} == (
        E68.SOURCE_LATTICES
    )


def test_fallback_evidence_tampering_fails_closed() -> None:
    item = E68.source_records()[0]
    record = item["record"]
    payload = json.loads(E68.evidence_path(record).read_text(encoding="utf-8"))
    forged = json.loads(json.dumps(payload))
    forged["support"] = [0]
    forged["weight"] = 1
    with pytest.raises(RuntimeError, match="failed validation"):
        E68.validate_evidence(forged, record)


def test_repair_worker_budget_is_capped() -> None:
    with pytest.raises(ValueError, match="one or two"):
        E68.run(3, 1.0)
