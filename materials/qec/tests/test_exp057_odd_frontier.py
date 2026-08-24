"""Regression checks for EXP-057's exact [[170,16,10]] frontier certificate."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp057", ROOT / "experiments" / "exp057_odd_frontier.py"
)
E57 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E57
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E57)


def test_target_rebuild_and_structural_reductions() -> None:
    problem = E57.build_problem()
    assert (problem["ell"], problem["m"], problem["n"], problem["k"]) == (
        17,
        5,
        170,
        16,
    )
    assert not (problem["HX"] @ problem["HZ"].T % 2).any()

    parity = E57.E56.kernel_weight_parity_proof(
        problem, E57.TARGET["lower_cap"]
    )
    assert parity["valid"] is True
    assert parity["effective_even_cap"] == 8
    assert E57.E56.bb_duality_proof(problem)["valid"] is True

    upper = E57.E56.verify_upper_witness(
        problem, E57.E56.default_witness(problem)
    )
    assert upper["weight"] == 10
    assert upper["z_valid_numpy"] and upper["z_valid_bitset"]
    assert upper["x_valid_numpy"] and upper["x_valid_bitset"]


def test_persisted_exact_certificate_rebuilds_current_protocol() -> None:
    certificate = json.loads(E57.CERTIFICATE.read_text(encoding="utf-8"))
    E57.validate_exact_certificate_payload(certificate)

    assert certificate["verdict"] == {
        "classification": "CERTIFIED_EXACT",
        "d": 10,
        "d_X": 10,
        "d_Z": 10,
        "exact": True,
    }
    lower = certificate["lower_bound"]
    assert lower["requested_cap"] == 9
    assert lower["effective_even_cap"] == 8
    assert lower["initial"]["status"] == "UNSAT"
    assert lower["replay"]["status"] == "UNSAT"
    assert lower["initial"]["cnf_sha256"] == lower["replay"]["cnf_sha256"]

    tampered = json.loads(json.dumps(certificate))
    tampered["lower_bound"]["replay"]["cnf_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="does not validate"):
        E57.validate_exact_certificate_payload(tampered)


def test_certificate_witness_is_physical_independently() -> None:
    problem = E57.build_problem()
    certificate = json.loads(E57.CERTIFICATE.read_text(encoding="utf-8"))
    support = certificate["upper_bound"]["support"]
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[support] = 1
    assert int(vector.sum()) == 10
    assert not (problem["HX"] @ vector % 2).any()
    assert E57.rank_np(np.vstack([problem["HZ"], vector])) == E57.rank_np(
        problem["HZ"]
    ) + 1
