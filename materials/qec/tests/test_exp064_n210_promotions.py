"""Regression checks for EXP-064's replayed n=210 exact references."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp064", ROOT / "experiments" / "exp064_n210_promotions.py"
)
E64 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E64
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E64)


@pytest.mark.parametrize(
    "key,k,d",
    [("k24d4", 24, 4), ("k14d12", 14, 12), ("k10d16", 10, 16)],
)
def test_exact_promotion_rebuilds_current_protocol(key: str, k: int, d: int) -> None:
    problem = E64.problem_of(key)
    assert (problem["n"], problem["k"]) == (210, k)
    certificate = json.loads(
        E64.certificate_path(key).read_text(encoding="utf-8")
    )
    E64.validate_exact_certificate_payload(certificate)
    assert certificate["verdict"] == {
        "classification": "CERTIFIED_EXACT",
        "d": d,
        "d_X": d,
        "d_Z": d,
        "exact": True,
    }
    lower = certificate["lower_bound"]
    assert lower["initial"]["status"] == lower["replay"]["status"] == "UNSAT"
    assert lower["initial"]["cnf_sha256"] == lower["replay"]["cnf_sha256"]
    assert certificate["upper"]["weight"] == d
    assert certificate["parity"]["valid"] is True
    assert certificate["duality"]["valid"] is True


def test_promotion_tampering_fails_closed() -> None:
    certificate = json.loads(
        E64.certificate_path("k14d12").read_text(encoding="utf-8")
    )
    certificate["lower_bound"]["initial"]["cnf_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="does not validate"):
        E64.validate_exact_certificate_payload(certificate)
