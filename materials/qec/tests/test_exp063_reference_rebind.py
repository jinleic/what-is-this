"""Regression checks for EXP-063's proof-preserving reference ratchet."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp063", ROOT / "experiments" / "exp063_reference_rebind.py"
)
E63 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E63
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E63)

SCREEN_DIR = ROOT / "results" / "partial_runs" / "exp055_screen"


def test_n210_reference_ratchet_closes_every_class() -> None:
    shard = json.loads((SCREEN_DIR / "15x7.json").read_text(encoding="utf-8"))
    assert shard["protocol"]["reference_validation_version"] == (
        E63.E55.REFERENCE_VALIDATION_VERSION
    )
    assert shard["candidates_after_symmetry"] == len(shard["records"]) == 419
    assert shard["survivors"] == shard["undecided"] == shard["no_reference"] == []
    assert shard["verdicts"] == {
        "dominated": 14,
        "dominated_by_cdcl_witness": 68,
        "dominated_by_witness": 337,
    }
    rebind = shard["reference_rebind"]
    assert rebind["schema"] == "exp063-monotone-reference-rebind-v1"
    assert rebind["all_thresholds_monotone"] is True
    assert rebind["all_carried_witnesses_rechecked"] is True
    assert (ROOT / rebind["archive"]).exists()
    E63._validate_bound_archive(shard)
    forged = json.loads(json.dumps(shard))
    forged["reference_rebind"]["old_shard_sha256"] = "f" * 64
    with pytest.raises(RuntimeError, match="archive hash"):
        E63._validate_bound_archive(forged)

    ratcheted = [record for record in shard["records"] if record.get("reference_rebind_ratchet")]
    assert len(ratcheted) == 6
    assert {(record["k_parent"], record["witness_bound"]) for record in ratcheted} == {
        (24, 4), (14, 12), (10, 16), (8, 16)
    }
    assert all(
        record["witness_bound"] <= record["threshold"]
        for record in ratcheted
    )


def test_rebind_missing_witness_support_fails_closed() -> None:
    shard = json.loads((SCREEN_DIR / "15x7.json").read_text(encoding="utf-8"))
    physical = next(
        record for record in shard["records"]
        if record["verdict"] in {
            "dominated_by_witness", "dominated_by_cdcl_witness"
        }
    )
    physical = json.loads(json.dumps(physical))
    physical.pop("witness_support", None)
    assert E63._physical_witness_valid(physical) is False

    ceiling = next(
        record for record in shard["records"]
        if record.get("ceiling_witness_support")
    )
    ceiling = json.loads(json.dumps(ceiling))
    ceiling.pop("ceiling_witness_support", None)
    assert E63._ceiling_witness_valid(ceiling) is False


def test_all_final_shards_share_the_current_reference_hash() -> None:
    current_hash = E63.E55._reference_fingerprint()
    names = [
        "3x3", "5x3", "7x3", "9x3", "9x5", "15x3", "7x7", "9x7",
        "21x3", "25x3", "17x5", "9x9", "15x5", "27x3", "31x3",
        "11x9", "33x3", "15x7", "21x5", "35x3",
    ]
    for name in names:
        shard = json.loads((SCREEN_DIR / f"{name}.json").read_text(encoding="utf-8"))
        assert shard["protocol"]["reference_sha256"] == current_hash
        assert shard["protocol"]["reference_validation_version"] == (
            E63.E55.REFERENCE_VALIDATION_VERSION
        )
