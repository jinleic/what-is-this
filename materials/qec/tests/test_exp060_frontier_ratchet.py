"""Regression checks for EXP-060's adaptive n=210 frontier ratchet."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp060", ROOT / "experiments" / "exp060_frontier_ratchet.py"
)
E60 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E60
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E60)


def test_promoted_target_rebuilds_and_reduces_to_one_even_sector() -> None:
    problem = E60.promoted_problem()
    assert (problem["ell"], problem["m"], problem["n"], problem["k"]) == (
        15,
        7,
        210,
        18,
    )
    parity = E60.E56.kernel_weight_parity_proof(
        problem, E60.PROMOTED_TARGET["lower_cap"]
    )
    assert parity["valid"] is True
    assert parity["effective_even_cap"] == 6
    assert E60.E56.bb_duality_proof(problem)["valid"] is True


def test_all_thirteen_screen_survivors_are_exact_210_18_8() -> None:
    summary = json.loads(E60.SUMMARY.read_text(encoding="utf-8"))
    assert summary["groups"] == ["survivors"]
    assert summary["all_exact"] is True
    assert summary["exact_parameter_histogram"] == {"[[210,18,8]]": 13}
    assert len(summary["records"]) == 13
    for record in summary["records"]:
        assert record["k"] == 18
        assert record["threshold"] == 6
        assert record["lower_bound"] == record["exact_distance"] == 8
        assert record["upper"]["weight"] == 8
        assert record["decisions"][0]["status"] == "UNSAT"
        assert record["duality"]["valid"] is True
        assert record["parity"]["valid"] is True


def test_promoted_exact_certificate_rebuilds_current_protocol() -> None:
    certificate = json.loads(E60.CERTIFICATE.read_text(encoding="utf-8"))
    E60.validate_exact_certificate_payload(certificate)
    assert certificate["verdict"] == {
        "classification": "CERTIFIED_EXACT",
        "d": 8,
        "d_X": 8,
        "d_Z": 8,
        "exact": True,
    }
    assert certificate["ratchet_family"] == {
        "exact_classes": 13,
        "parameters": "[[210,18,8]]",
        "all_exact": True,
    }
    lower = certificate["lower_bound"]
    assert lower["initial"]["status"] == lower["replay"]["status"] == "UNSAT"
    assert lower["initial"]["cnf_sha256"] == lower["replay"]["cnf_sha256"]

    tampered = json.loads(json.dumps(certificate))
    tampered["upper"]["weight"] = 6
    with pytest.raises(RuntimeError, match="does not validate"):
        E60.validate_exact_certificate_payload(tampered)
