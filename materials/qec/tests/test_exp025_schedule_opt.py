from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp025", ROOT / "experiments" / "exp025_schedule_opt.py"
)
exp = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = exp
SPEC.loader.exec_module(exp)


def _candidate(schedule_hash: str) -> dict:
    return {"schedule_hash": schedule_hash}


def _catalog() -> dict:
    return {
        "protocol_sha256": exp.protocol_hash(),
        "schedule_class": "translation-invariant",
        "codes": {
            "gross": {
                "schedule_class": "translation-invariant",
                "population_count": 8496,
                "population_exhaustive_scope": "translation-invariant only",
                "depth_claim": "exact(OPTIMAL) within translation-invariant class",
                "default_schedule_hash": "gd",
                "candidates": [_candidate("g1"), _candidate("g2"), _candidate("gd")],
            },
            "pbb": {
                "schedule_class": "translation-invariant",
                "population_count": 9968,
                "population_exhaustive_scope": "translation-invariant only",
                "depth_claim": "class-free exact: max-check-weight=8 plus witness",
                "default_schedule_hash": None,
                "candidates": [_candidate("p1"), _candidate("p2"), _candidate("p3")],
            },
        },
    }


def _row(schedule_hash: str, failures: int, shots: int) -> dict:
    return {
        "schedule_hash": schedule_hash,
        "operational_failures": failures,
        "shots": shots,
        "protocol_sha256": exp.protocol_hash(),
        "fatal_errors": 0,
        "circuit_sha256": f"circuit-{schedule_hash}",
        "paired": False,
    }


def test_protocol_v2_is_ti_scoped_and_fresh():
    protocol = exp.protocol()
    assert protocol["version"] == 2
    assert protocol["schedule_class"] == "translation-invariant"
    assert protocol["revision"]["v1_data_reused"] is False
    assert protocol["sampling"]["K_per_code"] == 16
    assert protocol["stage1"]["shots_per_schedule"] == 1000
    assert protocol["stage2"]["selected_per_code"] == 3
    assert protocol["stage2"]["shots_per_selected_schedule"] == 10000
    assert protocol["statistics"]["cross_circuit_pairing"] is False
    assert "within translation-invariant class" in (
        protocol["schedule_populations"]["gross"]["depth_claim"]
    )
    assert "class-free exact" in protocol["schedule_populations"]["pbb"]["depth_claim"]


def test_frozen_sampler_is_deterministic_without_replacement():
    first = exp.frozen_sample_indices()
    second = exp.frozen_sample_indices()
    assert first == second
    assert len(first["gross"]) == len(set(first["gross"])) == 16
    assert len(first["pbb"]) == len(set(first["pbb"])) == 16
    assert all(0 <= index < 8496 for index in first["gross"])
    assert all(0 <= index < 9968 for index in first["pbb"])


def test_scope_regression_rejects_unqualified_ti_exhaustiveness():
    catalog = _catalog()
    exp.validate_catalog_scope(catalog)
    catalog["codes"]["gross"]["population_exhaustive_scope"] = "all schedules"
    with pytest.raises(RuntimeError, match="scope widened"):
        exp.validate_catalog_scope(catalog)


def test_selection_uses_failure_count_then_hash_and_requires_complete_codes():
    catalog = _catalog()
    stage1 = {
        "gross": {
            "g1": _row("g1", 3, 1000),
            "g2": _row("g2", 1, 1000),
            "gd": _row("gd", 1, 1000),
        },
        "pbb": {
            "p1": _row("p1", 4, 1000),
            "p2": _row("p2", 2, 1000),
            "p3": _row("p3", 3, 1000),
        },
    }
    selection = exp.select_stage2(stage1, catalog)
    assert selection["gross"] == ["g2", "gd", "g1"]
    assert selection["pbb"] == ["p2", "p3", "p1"]
    del stage1["pbb"]["p1"]
    with pytest.raises(RuntimeError, match="incomplete/asymmetric"):
        exp.select_stage2(stage1, catalog)


def test_cross_circuit_comparison_is_independent_and_rejects_pairing():
    first = _row("g1", 4, 1000)
    second = _row("p1", 20, 1000)
    result = exp.independent_exact_test(first, second, alternative="first_less")
    assert result["test"] == "one-sided Fisher exact on independent circuits"
    assert 0 <= result["p_value"] <= 1
    first["paired"] = True
    with pytest.raises(ValueError, match="McNemar/paired"):
        exp.independent_exact_test(first, second, alternative="first_less")
    first["paired"] = False
    second["circuit_sha256"] = first["circuit_sha256"]
    with pytest.raises(ValueError, match="distinct circuit"):
        exp.independent_exact_test(first, second, alternative="first_less")


def test_holm_stepdown_and_directional_verdicts():
    correction = exp.holm({"a": 0.001, "b": 0.02, "c": 0.9})
    assert correction["a"]["rejected"]
    assert correction["b"]["rejected"]
    assert not correction["c"]["rejected"]
    pbb = _row("p1", 1, 1000)
    gross = _row("g1", 100, 1000)
    assert exp.determine_verdict(pbb, gross)["verdict"] == "BREAKTHROUGH_CANDIDATE"
    assert exp.determine_verdict(gross, pbb)["verdict"] == "NEGATIVE"
    pbb["operational_failures"] = gross["operational_failures"]
    assert exp.determine_verdict(pbb, gross)["verdict"] == "INCONCLUSIVE"


def test_canonicalization_rejects_asymmetric_incomplete_coverage(tmp_path, monkeypatch):
    catalog = _catalog()
    selection = {"gross": ["g1", "g2", "gd"], "pbb": ["p1", "p2", "p3"]}
    state = {
        "protocol_sha256": exp.protocol_hash(),
        "selection": selection,
        "stage1": {
            "gross": {key: _row(key, 1, 1000) for key in ("g1", "g2", "gd")},
            "pbb": {key: _row(key, 1, 1000) for key in ("p1", "p2", "p3")},
        },
        "stage2": {
            "gross": {key: _row(key, 10, 10000) for key in ("g1", "g2", "gd")},
            "pbb": {key: _row(key, 10, 10000) for key in ("p1", "p2")},
        },
    }
    monkeypatch.setattr(exp, "PROCESSED_RESULT_PATH", tmp_path / "canonical.json")
    with pytest.raises(RuntimeError, match="incomplete/asymmetric stage2 coverage for pbb"):
        exp.canonicalize(state, catalog)
    assert not (tmp_path / "canonical.json").exists()


def test_partial_route_cannot_target_canonical_paths():
    assert exp.canonical_route(True, False) == "partial_runs"
    assert exp.canonical_route(False, True) == "quarantine"
    assert exp.canonical_route(True, True) == "canonical"
    gate = exp.harness_gate_status()
    assert gate["passed"] is False
    with pytest.raises(RuntimeError, match="production is blocked"):
        exp.require_production_gate()
