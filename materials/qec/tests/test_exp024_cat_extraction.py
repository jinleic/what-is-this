from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp024", ROOT / "experiments" / "exp024_cat_extraction.py"
)
exp = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = exp
SPEC.loader.exec_module(exp)


def _structural() -> dict:
    return json.loads(exp.STRUCTURAL_OUTPUT.read_text())


def _schedule(payload: dict) -> dict:
    data_slot = {(int(c), int(q)): int(t) for c, q, _half, t in payload["data_slot"]}
    half = {(int(c), int(q)): int(h) for c, q, h, _t in payload["data_slot"]}
    return {
        "depth": int(payload["depth"]),
        "data_slot": data_slot,
        "half": half,
        "prep_slot": {int(c): int(t) for c, t in payload["prep_slot"]},
        "merge_slot": {int(c): int(t) for c, t in payload["merge_slot"]},
    }


def test_exact_target_and_independent_resource_accounting() -> None:
    row, spec = exp._load_target()
    assert row["code_id"] == "12_6_0193"
    assert (row["ell"], row["m"], row["n"], row["k"], row["d"]) == (
        12, 6, 144, 12, 12
    )
    assert row["A_terms"] == [[1, 2], [10, 3], [10, 4]]
    assert row["B_terms"] == [[0, 0], [1, 5], [11, 4]]
    assert row["C_terms"] == [[1, 3], [10, 3]]
    assert row["D_terms"] == [[1, 5], [10, 5]]

    code, supports, _ = exp.pbb_supports_and_orbits(spec)
    weights = [support.weight for support in supports]
    assert code.n == 144
    assert weights.count(8) == 72
    assert weights.count(6) == 72
    assert sum(weights) == 1008
    assert len(supports) + weights.count(8) == 216
    assert sum(weights) + 2 * weights.count(8) == 1152

    rows = {row["label"]: row for row in _structural()["resource_table"]}
    cat = rows["PBB two-ancilla unverified cat"]
    assert (
        cat["data_qubits"], cat["ancilla_qubits"],
        cat["total_physical_qubits"], cat["total_two_qubit_gates_per_round"],
        cat["total_two_qubit_gate_layer_depth_per_round"],
    ) == (144, 216, 360, 1152, 9)


def test_independent_schedule_and_scope_verification() -> None:
    result = _structural()
    row, spec = exp._load_target()
    code, supports, _ = exp.pbb_supports_and_orbits(spec)
    schedule = _schedule(result["schedule_search"]["actual_schedule"])
    verification = exp._verify_cat_schedule(supports, code.n, schedule)
    assert verification["valid"]
    assert verification["physical_qubit_conflict_free"]
    assert verification["semantic_parity_defects"] == 0
    assert verification["all_weight8_splits_are_4_plus_4"]
    assert schedule["depth"] == 9

    search = result["schedule_search"]
    assert search["integrated_depth7"]["status"] == "INFEASIBLE"
    assert search["depth7_claim"] == "exact(INFEASIBLE)"
    assert search["integrated_depth8_budgeted"]["status"] == "UNKNOWN"
    assert search["integrated_optimum_bracket"] == [8, 9]
    assert search["data_only"]["depth"] == 7
    assert search["data_only"]["status"] == "OPTIMAL"

    pbb_slot, _ = exp._exp016_slot("nonCSS-PBB [[144,12,12]] 12_6_0193")
    assert exp.verify_schedule(supports, pbb_slot)["valid"]
    resource_rows = {row["label"]: row for row in result["resource_table"]}
    assert "class-free" in resource_rows["PBB one ancilla"]["total_two_qubit_depth_status"]
    assert "translation-invariant" in resource_rows["CSS Gross one ancilla"]["total_two_qubit_depth_status"]


def test_independent_4096_shot_semantics_and_observables() -> None:
    result = _structural()
    _, spec = exp._load_target()
    code, supports, _ = exp.pbb_supports_and_orbits(spec)
    schedule = _schedule(result["schedule_search"]["actual_schedule"])
    observables = exp.pure_z_logical_basis(code)
    cat, _ = exp.build_cat_memory_circuit(
        code.H, observables, exp.ROUNDS, 0.0, schedule
    )
    det, obs = cat.compile_detector_sampler(seed=exp.NOISELESS_SEED).sample(
        4096, separate_observables=True
    )
    assert int(det.sum()) == 0
    assert int(obs.sum()) == 0
    assert obs.shape[1] == 12

    pbb_slot, source = exp._exp016_slot(
        "nonCSS-PBB [[144,12,12]] 12_6_0193"
    )
    pbb, _ = exp.build_memory_circuit(
        exp.CircuitSpec(
            H=code.H,
            observables=observables,
            rounds=exp.ROUNDS,
            p=0.0,
            basis="Z",
            layers=exp.slots_to_layers(pbb_slot, source["depth"]),
        )
    )
    assert exp._observable_definitions(cat) == exp._observable_definitions(pbb)
    digest = hashlib.sha256(
        "\n".join(exp._observable_definitions(cat)).encode()
    ).hexdigest()
    assert digest == result["semantics"]["observable_definition_sha256"]


def test_routing_checkpoint_and_frozen_v2_guard(tmp_path, monkeypatch) -> None:
    assert exp.canonical_route(True, True) == "canonical"
    assert exp.canonical_route(True, False) == "partial_runs"
    assert exp.canonical_route(False, False) == "quarantine"
    assert exp.OUTPUT.exists() and exp.OUTPUT.stat().st_size == 0
    sentinel_hash = hashlib.sha256(exp.OUTPUT.read_bytes()).hexdigest()

    monkeypatch.setattr(exp, "V2_GATE", tmp_path / "missing-gate.json")
    try:
        exp.require_benchmark_v2_execution_gate()
    except RuntimeError as error:
        assert "frozen but blocked" in str(error)
    else:
        raise AssertionError("missing shared-harness gate must block v2 execution")

    legacy_called = False

    def legacy_benchmark(*_args, **_kwargs):
        nonlocal legacy_called
        legacy_called = True
        raise AssertionError("legacy benchmark must be unreachable")

    monkeypatch.setattr(exp, "_benchmark", legacy_benchmark)
    try:
        exp.run_benchmark_v2_production()
    except RuntimeError as error:
        assert "frozen but blocked" in str(error)
    else:
        raise AssertionError("normal production entry must fail closed")
    assert not legacy_called

    target_load_called = False

    def load_target():
        nonlocal target_load_called
        target_load_called = True
        raise AssertionError("production gate must precede structural work")

    monkeypatch.setattr(exp, "_load_target", load_target)
    try:
        exp.run(workers=1, validate_only=False)
    except RuntimeError as error:
        assert "frozen but blocked" in str(error)
    else:
        raise AssertionError("normal run must fail closed before structural work")
    assert not target_load_called

    monkeypatch.setattr(exp, "V2_CHECKPOINT_DIR", tmp_path / "checkpoints")
    binding = exp.benchmark_v2_binding(
        "development", "cat_pbb", "circuit", 7, "sample", 3
    )
    records = [{"shot_id": 0}, {"shot_id": 1}]
    path = exp.save_benchmark_v2_arm_checkpoint(binding, records)
    assert path.exists()
    assert not path.with_suffix(path.suffix + ".tmp").exists()
    assert exp.load_benchmark_v2_arm_checkpoint(binding) == records
    changed = {**binding, "seed": 8}
    try:
        exp.load_benchmark_v2_arm_checkpoint(changed)
    except RuntimeError as error:
        assert "binding mismatch" in str(error)
    else:
        raise AssertionError("checkpoint binding drift must be rejected")

    assert hashlib.sha256(exp.OUTPUT.read_bytes()).hexdigest() == sentinel_hash
    amendment = json.loads(exp.FEASIBILITY_AMENDMENT.read_text())
    assert amendment["status"] == "V1_RETIRED_WITHOUT_RESULT_V2_FROZEN_NOT_EXECUTED"
    assert amendment["v2_frozen_protocol"]["execution_gate"]["status"] == "BLOCKED"
    assert amendment["v2_frozen_protocol"]["result"] is None
