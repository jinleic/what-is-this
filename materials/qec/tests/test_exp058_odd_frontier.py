"""Regression checks for EXP-058's exact [[186,10,14]] frontier certificate."""
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
    "exp058", ROOT / "experiments" / "exp058_odd_frontier.py"
)
E58 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E58
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E58)


def test_target_rebuild_and_structural_reductions() -> None:
    problem = E58.build_problem()
    assert (problem["ell"], problem["m"], problem["n"], problem["k"]) == (
        31,
        3,
        186,
        10,
    )
    assert not (problem["HX"] @ problem["HZ"].T % 2).any()

    parity = E58.E56.kernel_weight_parity_proof(
        problem, E58.TARGET["lower_cap"]
    )
    assert parity["valid"] is True
    assert parity["effective_even_cap"] == 12
    assert E58.E56.bb_duality_proof(problem)["valid"] is True

    upper = E58.E56.verify_upper_witness(
        problem, E58.E56.default_witness(problem)
    )
    assert upper["weight"] == 14
    assert upper["z_valid_numpy"] and upper["z_valid_bitset"]
    assert upper["x_valid_numpy"] and upper["x_valid_bitset"]


def test_persisted_exact_certificate_rebuilds_current_protocol() -> None:
    certificate = json.loads(E58.CERTIFICATE.read_text(encoding="utf-8"))
    E58.validate_exact_certificate_payload(certificate)

    assert certificate["verdict"] == {
        "classification": "CERTIFIED_EXACT",
        "d": 14,
        "d_X": 14,
        "d_Z": 14,
        "exact": True,
    }
    lower = certificate["lower_bound"]
    assert lower["requested_cap"] == 13
    assert lower["effective_even_cap"] == 12
    assert lower["initial"]["status"] == "UNSAT"
    assert lower["replay"]["status"] == "UNSAT"
    assert lower["initial"]["cnf_sha256"] == lower["replay"]["cnf_sha256"]
    assert lower["initial"]["solver"]["name"] == "PySAT kissat404"

    sibling = certificate["sibling_classes"]
    assert sibling["count"] == 6
    assert sibling["all_both_sector_unsat_cap12"] is True
    assert sibling["all_weight14_witness_found"] is True
    for relative in sibling["evidence_files"]:
        record = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        if record["cap"] == 12:
            assert record["status"] == "UNSAT"
        else:
            assert record["cap"] == 14 and record["weight"] == 14

    tampered = json.loads(json.dumps(certificate))
    tampered["lower_bound"]["replay"]["cnf_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="does not validate"):
        E58.validate_exact_certificate_payload(tampered)


def test_certificate_witness_is_physical_independently() -> None:
    problem = E58.build_problem()
    certificate = json.loads(E58.CERTIFICATE.read_text(encoding="utf-8"))
    support = certificate["upper_bound"]["support"]
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[support] = 1
    assert int(vector.sum()) == 14
    assert not (problem["HX"] @ vector % 2).any()
    assert E58.rank_np(np.vstack([problem["HZ"], vector])) == E58.rank_np(
        problem["HZ"]
    ) + 1
