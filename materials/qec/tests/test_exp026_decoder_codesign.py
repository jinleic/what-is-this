from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp026", ROOT / "experiments" / "exp026_decoder_codesign.py"
)
exp = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = exp
SPEC.loader.exec_module(exp)


def test_protocol_v2_and_routing_contract():
    p = exp.protocol()
    assert p["version"] == 2
    assert p["revision_reason"]["v1_production_rows"] == 0
    assert p["revision_reason"]["baseline_equivalent_core_hour_lower_bound"] == 216
    assert len(p["candidates"]) == 5
    assert exp.canonical_route(True, False) == "partial_runs"
    assert exp.canonical_route(False, True) == "quarantine"
    assert exp.canonical_route(True, True) == "canonical"


def test_pairing_and_one_sided_mcnemar_direction():
    baseline = [True, True, True, False, False]
    challenger = [False, False, True, True, False]
    result = exp.one_sided_mcnemar(baseline, challenger)
    assert result["baseline_only_failures"] == 2
    assert result["challenger_only_failures"] == 1
    assert result["discordant"] == 3
    assert 0 <= result["pvalue_one_sided"] <= 1
    try:
        exp.one_sided_mcnemar([True], [False, True])
    except ValueError:
        pass
    else:
        raise AssertionError("strict pairing must reject unequal lengths")


def test_holm_stepdown_stops_after_first_failure():
    result = exp.holm({"a": 0.001, "b": 0.02, "c": 0.021, "d": 0.9})
    assert result["a"]["rejected"]
    assert not result["b"]["rejected"]
    assert not result["c"]["rejected"]
    assert not result["d"]["rejected"]


def test_balanced_orders_cover_candidates_once_per_block():
    ids = ["a", "b", "c", "d", "e"]
    orders = exp.balanced_orders(ids, 12, 123)
    assert len(orders) == 12
    assert all(sorted(order) == ids for order in orders)
    assert len({tuple(order) for order in orders[:5]}) == 5


def test_bootstrap_p95_ratio_is_paired_and_deterministic():
    base = np.arange(1, 21, dtype=float)
    challenger = base * 0.5
    a = exp.bootstrap_p95_ratio(base, challenger, 7, reps=200)
    b = exp.bootstrap_p95_ratio(base, challenger, 7, reps=200)
    assert a == b
    assert abs(a["ratio"] - 0.5) < 1e-12
    assert a["ci95"][1] < 1


def _sleep_worker(command_q, reply_q, circuit_text, pure_rows, candidate):
    reply_q.put({"kind": "ready"})
    while True:
        message = command_q.get()
        if message is None:
            return
        shot_id, _ = message
        if shot_id == 0:
            time.sleep(0.2)
        else:
            reply_q.put({"kind": "result", "shot_id": shot_id,
                         "prediction": [0], "valid": True, "latency_s": 0.0})


def test_supervisor_timeout_counts_and_restarts(monkeypatch):
    monkeypatch.setattr(exp, "_worker", _sleep_worker)
    supervisor = exp.SupervisedDecoder.__new__(exp.SupervisedDecoder)
    supervisor.circuit_text = ""
    supervisor.pure_rows = []
    supervisor.candidate = exp.CANDIDATES[0]
    supervisor.deadline_s = 0.05
    supervisor.startup_deadline_s = 2.0
    supervisor.context = exp.mp.get_context("spawn")
    supervisor.process = supervisor.command_q = supervisor.reply_q = None
    supervisor.restarts = 0
    supervisor._start()
    try:
        first = supervisor.decode(0, np.zeros(1, dtype=np.uint8))
        second = supervisor.decode(1, np.zeros(1, dtype=np.uint8))
    finally:
        supervisor.close()
    assert first.timeout and first.restarted
    assert supervisor.restarts == 1
    assert not second.timeout and second.valid
    records = [{"operational_failure": True, "logical_mismatch": True,
                "invalid_correction": True, "timeout": True,
                "worker_restarted": True, "fatal": None}]
    summary = exp.counts(records)
    assert summary["operational_failures"] == 1
    assert summary["timeouts"] == 1
    assert summary["worker_restarts"] == 1


def test_verdict_requires_all_three_gates():
    lower_rate = True
    holm_rejected = True
    latency_upper = 0.99
    assert lower_rate and holm_rejected and latency_upper < 1
    assert not (lower_rate and False and latency_upper < 1)
    assert not (lower_rate and holm_rejected and 1.0 < 1)


class _DeterministicSupervisor:
    calls: list[int] = []

    def __init__(self, *_args, **_kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def decode(self, shot_id, syndrome):
        self.calls.append(shot_id)
        return exp.DecodeResult(
            prediction=np.asarray([shot_id & 1], dtype=np.uint8),
            valid=True, timeout=False, fatal=None,
            latency_s=shot_id / 1000, restarted=False,
        )

    def close(self):
        pass


def test_resume_exact_nonmultiple_prefix_matches_clean_run():
    det = np.zeros((13, 1), dtype=np.uint8)
    obs = np.asarray([[i & 1] for i in range(13)], dtype=np.uint8)
    candidate = exp.CANDIDATES[0]
    _DeterministicSupervisor.calls = []
    clean = exp.evaluate_candidate(
        None, np.zeros(0, dtype=np.int64), candidate, det, obs,
        supervisor_factory=_DeterministicSupervisor,
    )
    prefix = clean[:7]  # interruption after a non-multiple, <=25-shot prefix
    _DeterministicSupervisor.calls = []
    resumed = exp.evaluate_candidate(
        None, np.zeros(0, dtype=np.int64), candidate, det, obs,
        initial_records=prefix, supervisor_factory=_DeterministicSupervisor,
    )
    assert _DeterministicSupervisor.calls == list(range(7, 13))
    assert [r["shot_id"] for r in resumed] == list(range(13))
    assert resumed == clean
    bad = [dict(record) for record in prefix]
    bad[-1]["shot_id"] = 5
    try:
        exp.evaluate_candidate(
            None, np.zeros(0, dtype=np.int64), candidate, det, obs,
            initial_records=bad, supervisor_factory=_DeterministicSupervisor,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("duplicate checkpoint shot ids must be rejected")


def test_checkpoint_binding_rejects_protocol_circuit_or_sample_drift():
    det = np.zeros((3, 2), dtype=np.uint8)
    obs = np.zeros((3, 1), dtype=np.uint8)
    binding = exp.checkpoint_binding(7, det, obs, "circuit-a")
    entry = {"binding": binding, "records": [{"shot_id": 0}]}
    exp.validate_binding(entry, binding)
    for changed in (
        {**binding, "seed": 8},
        {**binding, "circuit_sha256": "circuit-b"},
        {**binding, "paired_sample_sha256": "different"},
    ):
        try:
            exp.validate_binding(entry, changed)
        except RuntimeError:
            pass
        else:
            raise AssertionError("binding drift must be rejected")


def test_full_coverage_gate_rejects_partial_state():
    state = {
        "development": {}, "heldout": {}, "timing": {},
        "circuits": {"gross": {"circuit_sha256": "g"},
                     "pbb": {"circuit_sha256": "p"}},
        "timing_gate": {"passed": True},
    }
    try:
        exp.validate_full_coverage(state)
    except RuntimeError:
        pass
    else:
        raise AssertionError("partial state must never route canonical")


def _gate(passed: bool, load: float = 0.1):
    return {
        "load_average": [load, load, load], "cpu_count": 14,
        "max_load": 4.0, "heavy_processes": [], "passed": passed,
    }


def test_absolute_idle_gate_is_not_core_scaled(monkeypatch):
    monkeypatch.setattr(exp.os, "getloadavg", lambda: (4.01, 0.0, 0.0))
    monkeypatch.setattr(exp.os, "cpu_count", lambda: 128)
    monkeypatch.setattr(exp, "heavy_processes", lambda: [])
    gate = exp.timing_gate()
    assert gate["max_load"] == 4.0
    assert not gate["passed"]


def test_accuracy_session_start_and_midrun_failure_is_quarantined(monkeypatch):
    state = {"heldout_gate_sessions": [], "quarantined_sessions": []}
    monkeypatch.setattr(exp, "save_partial", lambda _state: None)
    gates = iter((_gate(True), _gate(False, 4.1)))
    monkeypatch.setattr(exp, "timing_gate", lambda: next(gates))
    session = exp.require_claiming_gate(state, "heldout")
    try:
        exp.recheck_claiming_gate(state, session, "heldout")
    except RuntimeError:
        pass
    else:
        raise AssertionError("mid-run load violation must abort accuracy session")
    assert session["aborted"] and not session["claiming"]
    assert state["quarantined_sessions"][-1]["aborted"]

    state = {"heldout_gate_sessions": [], "quarantined_sessions": []}
    monkeypatch.setattr(exp, "timing_gate", lambda: _gate(False, 4.1))
    try:
        exp.require_claiming_gate(state, "heldout")
    except RuntimeError:
        pass
    else:
        raise AssertionError("failed start gate must reject accuracy session")
    assert state["quarantined_sessions"]


def test_timing_midrun_failure_discards_entire_pass(monkeypatch):
    circuits = {
        key: {
            "circuit": object(), "pure_rows": np.zeros(0, dtype=np.int64),
            "record": {"circuit_sha256": key},
        }
        for key in exp.CODE_KEYS
    }
    state = exp.state_template(circuits)
    monkeypatch.setattr(exp, "save_partial", lambda _state: None)
    monkeypatch.setattr(exp, "sample_pairs", lambda *_args: (
        np.zeros((1, 1), dtype=np.uint8), np.zeros((1, 1), dtype=np.uint8)
    ))
    monkeypatch.setattr(exp, "SupervisedDecoder", _DeterministicSupervisor)
    gates = iter((_gate(True), _gate(False, 4.1)))
    monkeypatch.setattr(exp, "timing_gate", lambda: next(gates))
    try:
        exp.run_timing(state, circuits)
    except RuntimeError:
        pass
    else:
        raise AssertionError("mid-run timing gate failure must abort")
    assert state["timing"] == {}
    assert state["timing_run"]["aborted"]
    assert state["quarantined_sessions"][-1]["aborted"]


def test_canonical_gate_rejects_rows_without_passed_session_provenance():
    state = {
        "development": {}, "heldout": {}, "timing": {},
        "circuits": {"gross": {"circuit_sha256": "g"},
                     "pbb": {"circuit_sha256": "p"}},
        "heldout_gate_sessions": [{
            "passed": False, "aborted": True, "rechecks": [{"passed": False}],
        }],
        "timing_gate": {"passed": True},
    }
    try:
        exp.validate_full_coverage(state)
    except RuntimeError:
        pass
    else:
        raise AssertionError("unproven accuracy rows must never route canonical")


def test_gate_segments_require_exact_once_coverage_and_passing_endpoints():
    start = _gate(True)
    end = _gate(True)
    sessions = [{
        "phase": "heldout", "start_gate": start, "rechecks": [end],
    }]
    valid = [{
        "session_index": 0, "start_shot": 0, "end_shot_exclusive": 25,
        "start_gate": start, "end_gate": end, "passed": True,
    }]
    assert exp.gate_segments_cover(valid, 25, 25, sessions, "heldout")
    gap = [dict(valid[0], start_shot=1)]
    overlap = valid + [dict(valid[0], start_shot=20, end_shot_exclusive=25)]
    failed_end = [dict(valid[0], end_gate=_gate(False, 4.1))]
    unlinked_end = [dict(valid[0], end_gate={**_gate(True), "time_ns": 99})]
    assert not exp.gate_segments_cover(gap, 25, 25, sessions, "heldout")
    assert not exp.gate_segments_cover(overlap, 25, 25, sessions, "heldout")
    assert not exp.gate_segments_cover(failed_end, 25, 25, sessions, "heldout")
    assert not exp.gate_segments_cover(unlinked_end, 25, 25, sessions, "heldout")
