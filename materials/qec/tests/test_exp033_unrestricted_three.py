"""Machine-check the persisted EXP-033 unrestricted schedule certificates."""
from __future__ import annotations

import importlib.util
import hashlib
import json
import sys
from pathlib import Path

from qec_research.circuits.bicycle_schedule import pbb_supports_and_orbits
from qec_research.circuits.scheduling import generator_supports, verify_schedule

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "processed" / "exp033_unrestricted_three.json"
NOISELESS = ROOT / "results" / "processed" / "exp033_noiseless_check.json"

_SPEC = importlib.util.spec_from_file_location(
    "exp031_depth_criterion_for_exp033_test",
    ROOT / "experiments" / "exp031_depth_criterion.py",
)
_E31 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _E31
_SPEC.loader.exec_module(_E31)

_NOISELESS_SPEC = importlib.util.spec_from_file_location(
    "exp033_noiseless_check_for_test",
    ROOT / "experiments" / "exp033_noiseless_check.py",
)
_E33N = importlib.util.module_from_spec(_NOISELESS_SPEC)
sys.modules[_NOISELESS_SPEC.name] = _E33N
_NOISELESS_SPEC.loader.exec_module(_E33N)


def test_all_three_persist_exact_unrestricted_depth_nine_witnesses() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["experiment"] == "EXP-033"
    assert payload["verdict"] == "POSITIVE"
    assert payload["n_survive_two_value"] == 3
    assert payload["n_refuted_unrestricted"] == 0
    assert payload["n_undecided"] == 0

    instances = _E31._pbb_instances()
    assert len(payload["results"]) == 3
    for result in payload["results"]:
        index = int(result["exp031_instance"].removeprefix("pbb144-"))
        code, _, _ = pbb_supports_and_orbits(instances[index]["spec"])
        supports = generator_supports(code.H)
        w_max = max(support.weight for support in supports)
        slot = {(check, qubit): depth for check, qubit, depth in result["slot"]}
        verification = verify_schedule(supports, slot)

        # w_max is an unconditional lower bound in the one-ancilla model; this
        # valid depth-w_max witness therefore proves the exact minimum without
        # trusting CP-SAT's status string.
        assert w_max == result["w_max"] == 9
        assert max(slot.values()) + 1 == result["certified_depth"] == w_max
        assert verification["valid"]
        assert verification == result["verification"]
        assert result["trace"][-1]["status"] == "OPTIMAL"
        assert result["trace"][-1]["valid"] is True
        assert result["two_value_survives"] is True
        assert result["two_value_refuted"] is False


def test_noiseless_stim_check_covers_all_three_witnesses() -> None:
    """Prop C5 evidence standard: 0 firings / 0 flips noiselessly (12 rounds)."""
    payload = json.loads(NOISELESS.read_text())
    assert payload["verdict"] == "POSITIVE"
    assert payload["canonical_input"] == "results/processed/exp033_unrestricted_three.json"
    assert payload["canonical_input_sha256"] == hashlib.sha256(
        ARTIFACT.read_bytes()
    ).hexdigest()
    assert payload["artifact_route"] == "canonical"
    assert payload["artifact_path"] == "results/processed/exp033_noiseless_check.json"
    canonical = json.loads(ARTIFACT.read_text())
    checked = {row["exp031_instance"] for row in payload["results"]}
    assert checked == {r["exp031_instance"] for r in canonical["results"]}
    canonical_by_id = {
        row["exp031_instance"]: row for row in canonical["results"]
    }
    for row in payload["results"]:
        assert row["depth"] == 9
        source = canonical_by_id[row["exp031_instance"]]
        assert row["certified_depth"] == source["certified_depth"]
        assert row["w_max"] == source["w_max"]
        assert row["slot_sha256"] == _E33N.canonical_json_sha256(source["slot"])
        assert row["noiseless_shots"] >= 4000
        assert row["rounds"] == 12
        assert row["noiseless_detector_firings"] == 0
        assert row["noiseless_observable_flips"] == 0
        assert row["passes"] is True


def test_failing_noiseless_check_cannot_overwrite_canonical(
    tmp_path: Path, monkeypatch
) -> None:
    canonical = tmp_path / "results" / "processed" / "exp033_noiseless_check.json"
    canonical.parent.mkdir(parents=True)
    canonical.write_text("trusted canonical evidence")
    monkeypatch.setattr(_E33N, "ROOT", tmp_path)
    monkeypatch.setattr(_E33N, "OUT", canonical)

    destination = _E33N.write_result(
        {"verdict": "NEGATIVE"},
        clean=False,
        full_coverage=True,
        stamp="regression",
    )

    assert canonical.read_text() == "trusted canonical evidence"
    assert destination == (
        tmp_path / "results" / "quarantine"
        / "exp033_noiseless_check-regression.json"
    )
    quarantined = json.loads(destination.read_text())
    assert quarantined["artifact_route"] == "quarantine"
    assert quarantined["artifact_path"] == (
        "results/quarantine/exp033_noiseless_check-regression.json"
    )
