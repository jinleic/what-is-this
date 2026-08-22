"""Machine checks for EXP-036's SAT comparison protocol and identity gates."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp036_delta_closure", ROOT / "experiments" / "exp036_delta_closure.py"
)
assert SPEC and SPEC.loader
EXP036 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP036)

from qec_research.codes.bicycle import BRAVYI_BB, build_bb  # noqa: E402
from qec_research.distance.sat_decide import (  # noqa: E402
    css_logical_bases,
    css_side_instance,
    decide_weight_bounded,
    symplectic_instance,
    verify_witness_two_paths,
)


def test_css_logical_bases_shapes_and_pairing() -> None:
    HX, HZ = build_bb(BRAVYI_BB["[[72,12,6]]"])
    LX, LZ = css_logical_bases(HX, HZ)
    assert LX.shape == (12, 72)
    assert LZ.shape == (12, 72)


def test_decision_semantics_reproduce_known_exact_distance() -> None:
    # [[72,12,6]]: published exact d=6, previously double-verified in EXP-027.
    HX, HZ = build_bb(BRAVYI_BB["[[72,12,6]]"])
    for side in ("x", "z"):
        instance = css_side_instance(HX, HZ, side)
        below = decide_weight_bounded(instance, 5)
        at = decide_weight_bounded(instance, 6)
        assert below["status"] == "UNSAT"          # d > 5
        assert at["status"] == "SAT"               # d <= 6, verified witness
        assert at["weight"] == 6
        assert at["verification"]["valid"]


def test_witness_verification_rejects_corrupted_vector() -> None:
    HX, HZ = build_bb(BRAVYI_BB["[[72,12,6]]"])
    instance = css_side_instance(HX, HZ, "x")
    record = decide_weight_bounded(instance, 6)
    vector = np.asarray(record["vector"], dtype=np.uint8)
    corrupted = vector.copy()
    corrupted[0] ^= 1
    verification = verify_witness_two_paths(instance, corrupted)
    assert not verification["valid"]


def _row_by_code_id(code_id: str) -> tuple[int, dict]:
    rows = EXP036.E27.load_catalogue()
    for index, row in enumerate(rows):
        if str(row.get("code_id", "")) == code_id:
            return index, row
    raise AssertionError(f"code_id {code_id} not in catalogue")


def test_compare_certifies_known_reversal(tmp_path, monkeypatch) -> None:
    # phase2_58: parent [[72,8,4]] (exact d=4), PBB [[72,4,6]] (exact d=6) -
    # EXP-027's double-verified certified reversal.  The frozen protocol must
    # reproduce CERTIFIED_REVERSAL from scratch.
    monkeypatch.setattr(EXP036, "STATE_DIR", tmp_path / "state")
    index, row = _row_by_code_id("phase2_58")
    solver = EXP036.RowSolver(index, row, conflict_budget=0)
    state = solver.compare()
    assert state["verdict"] == "CERTIFIED_REVERSAL"
    assert state["U_p"] == 4
    assert state["U_b"] == 6
    assert state["bounds"]["pbb_lower_gt"] == 4  # d_pbb > 4 >= d_parent
    # tie-safety: a reversal claim must always show pbb strictly above parent
    assert state["U_b"] > state["U_p"]


def test_compare_certifies_known_domination(tmp_path, monkeypatch) -> None:
    # Smallest domination case: an n=36 delta>0 row EXP-027 already proved.
    artifact = json.loads(EXP036.EXP027_ARTIFACT.read_text())
    chosen = None
    for item in artifact["per_item_results"]:
        if item["n"] == 36 and item["verdict"] == "DOMINATION_PROVED":
            chosen = item
            break
    assert chosen is not None
    rows = EXP036.E27.load_catalogue()
    index = chosen["catalogue_index"]
    row = rows[index]
    monkeypatch.setattr(EXP036, "STATE_DIR", tmp_path / "state")
    solver = EXP036.RowSolver(index, row, conflict_budget=0)
    state = solver.compare()
    assert state["verdict"] == "DOMINATION_PROVED"
    assert state["U_p"] >= state["U_b"]


def test_identity_gate_fail_stops_on_tampered_state(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(EXP036, "STATE_DIR", tmp_path / "state")
    index, row = _row_by_code_id("phase2_58")
    identity = EXP036.row_identity(row, index)
    path = EXP036.row_state_path(index)
    tampered = dict(identity)
    tampered["pbb_fingerprint"] = "0" * 64
    tampered["verdict"] = "DOMINATION_PROVED"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(RuntimeError, match="identity mismatch"):
        EXP036.verified_row_state(index, row)


def test_forged_unsat_cannot_survive_replay(tmp_path, monkeypatch) -> None:
    # Threat model: field edits manufacture a verdict.  Two forgery classes
    # on the known reversal row (truth: parent d = 4, PBB d = 6):
    #   crude  - forged parent UNSAT@5 coexisting with the genuine weight-4
    #            parent witnesses: caught arithmetically, no solving;
    #   subtle - forged UNSAT@5 with the parent witnesses scrubbed: derives
    #            DOMINATION from the call log, is denied replay-verified
    #            status, and fail-stops under proof-grade replay.
    monkeypatch.setattr(EXP036, "STATE_DIR", tmp_path / "state")
    index, row = _row_by_code_id("phase2_58")
    solver = EXP036.RowSolver(index, row, conflict_budget=0)
    state = solver.compare()
    assert state["verdict"] == "CERTIFIED_REVERSAL"
    forged_calls = [
        {
            "target": target,
            "kind": f"css_{target[-1]}",
            "weight_cap": 5,
            "status": "UNSAT",
            "cnf_sha256": "0" * 64,
            "encoding_version": "sat-decide-v1",
            "solver": {"name": "PySAT cadical195", "wall_time_s": 0.0},
            "utc": "forged",
        }
        for target in ("parent_x", "parent_z")
    ]
    instances = EXP036.build_instances(row)

    crude = json.loads(json.dumps(state))
    crude["calls"].extend(forged_calls)
    with pytest.raises(RuntimeError, match="internal contradiction"):
        EXP036.validated_bounds(instances, crude)

    subtle = json.loads(json.dumps(state))
    subtle["calls"] = [
        entry
        for entry in subtle["calls"]
        if not (entry["target"].startswith("parent") and entry["status"] == "SAT")
    ] + forged_calls
    subtle["witnesses"] = {
        key: value
        for key, value in subtle["witnesses"].items()
        if not key.startswith("parent")
    }
    subtle["verdict"] = "DOMINATION_PROVED"
    EXP036.atomic_write_json(EXP036.row_state_path(index), subtle)
    bounds = EXP036.validated_bounds(instances, subtle)
    derived, _ = EXP036.derive_verdict(bounds)
    # The scrubbed forgery derives DOMINATION - which is exactly why replay
    # stamps are mandatory before anything counts as certified.
    assert derived == "DOMINATION_PROVED"
    assert not EXP036.replay_stamps_valid(instances, derived, bounds, subtle)
    with pytest.raises(RuntimeError, match="replay refuted stored UNSAT"):
        EXP036.reverify_row(index, row)


def test_reverify_stamps_genuine_row(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(EXP036, "STATE_DIR", tmp_path / "state")
    index, row = _row_by_code_id("phase2_58")
    solver = EXP036.RowSolver(index, row, conflict_budget=0)
    state = solver.compare()
    assert state["verdict"] == "CERTIFIED_REVERSAL"
    stamped = EXP036.reverify_row(index, row)
    assert stamped["replay"]["verdict"] == "CERTIFIED_REVERSAL"
    assert len(stamped["replay"]["stamps"]) == 1
    stamp = stamped["replay"]["stamps"][0]
    assert stamp["target"] == "pbb"
    assert stamp["weight_cap"] == 4
    instances = EXP036.build_instances(row)
    bounds = EXP036.validated_bounds(instances, stamped)
    assert EXP036.replay_stamps_valid(
        instances, "CERTIFIED_REVERSAL", bounds, stamped
    )


def test_envelope_check_artifact_is_machine_rederivable() -> None:
    # The "CSS envelope stands" claim for the new reversals must re-derive
    # from replay-stamped evidence, never from prose.
    artifact_path = ROOT / "results" / "processed" / "exp036_envelope_check.json"
    if not artifact_path.exists():
        return
    spec = importlib.util.spec_from_file_location(
        "exp036_envelope_check", ROOT / "experiments" / "exp036_envelope_check.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    payload, all_dominated = module.build_payload()
    assert all_dominated
    stored = json.loads(artifact_path.read_text())
    assert stored["all_reversals_css_dominated_at_equal_n"]
    assert stored["num_reversals_checked"] == payload["num_reversals_checked"]
    for fresh, persisted in zip(payload["checks"], stored["checks"]):
        assert fresh["dominated"] and persisted["dominated"]
        assert fresh["catalogue_index"] == persisted["catalogue_index"]
        assert fresh["dominator"]["n"] == fresh["n"]
        assert fresh["dominator"]["k_css"] >= fresh["k_pbb"]
        assert fresh["dominator"]["d_css_exact"] >= fresh["d_pbb_upper"]


def test_envelope_check_covers_both_reversal_sources() -> None:
    # Regression: the checker once iterated only EXP-027's *undecided* rows, so
    # EXP-027's own two reversals (331/339) were structurally invisible and the
    # claim silently covered 5 of 7.  Scope must never shrink back.
    spec = importlib.util.spec_from_file_location(
        "exp036_envelope_check", ROOT / "experiments" / "exp036_envelope_check.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    payload, all_dominated = module.build_payload()
    assert all_dominated

    by_index = {c["catalogue_index"]: c for c in payload["checks"]}
    assert {331, 339} <= set(by_index), "EXP-027 reversals must be checked"
    assert payload["num_reversals_checked"] >= 7
    sources = {c["source"] for c in payload["checks"]}
    assert sources == {"exp027_artifact", "exp036_tie_branch"}

    # The n=72 pair needs a genuine same-length dominator, not a d_Z read as a
    # code distance: k_css >= 4 and certified-exact d_css >= 6.
    for index in (331, 339):
        check = by_index[index]
        assert check["n"] == 72 and check["d_pbb_upper"] == 6
        dom = check["dominator"]
        assert dom["n"] == 72 and dom["k_css"] >= check["k_pbb"]
        assert dom["d_css_exact"] >= check["d_pbb_upper"]

    # Dominator choice is canonical (strongest by k*d^2), so the reported
    # comparison cannot drift with iteration order.
    for check in payload["checks"]:
        dom = check["dominator"]
        best = max(
            c["k_css"] * c["d_css_exact"] ** 2
            for c in payload["exact_css_candidates"]
            if c["n"] == check["n"]
            and c["k_css"] >= check["k_pbb"]
            and c["d_css_exact"] >= check["d_pbb_upper"]
        )
        assert dom["k_css"] * dom["d_css_exact"] ** 2 == best
        assert check["kd2_over_n_css"] >= check["kd2_over_n_pbb"]
