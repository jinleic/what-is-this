"""Focused machine checks for EXP-035's PBB distance certificate protocol."""

from __future__ import annotations

import importlib.util
import os
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp035_pbb_distance", ROOT / "experiments" / "exp035_pbb_exact_distance.py"
)
assert SPEC and SPEC.loader
EXP035 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP035)


def target() -> tuple[dict, str, object, dict, np.ndarray, dict]:
    row, raw = EXP035.load_target_catalogue_row()
    code, metadata = EXP035.rebuild_metadata(row, raw)
    logicals, logical_metadata = EXP035.logical_basis_metadata(code)
    return row, raw, code, metadata, logicals, logical_metadata


def test_exact_source_terms_rebuild_general_non_css_code() -> None:
    row, raw, code, metadata, logicals, logical_metadata = target()
    assert {key: row[key] for key in EXP035.EXPECTED_TERMS} == EXP035.EXPECTED_TERMS
    assert metadata["source"]["line"] == 38
    assert metadata["source"]["raw_line_sha256"] == EXP035.sha256_bytes(
        (raw + "\n").encode()
    )
    assert metadata["convention"]["symplectic_columns"] == "(x|z)"
    assert metadata["convention"]["check_matrix"] == "H=(A B | C D ; 0 0 | B^T A^T)"
    assert code.n == 144
    assert metadata["rebuild"]["rank_numpy"] == 132
    assert metadata["rebuild"]["rank_bitset"] == 132
    assert metadata["rebuild"]["k"] == 12
    assert metadata["rebuild"]["commutes_numpy"]
    assert metadata["rebuild"]["commutes_bitset"]
    assert metadata["rebuild"]["non_css_as_stored"]
    assert metadata["rebuild"]["independent_builder_matrix_matches_existing_builder"]
    assert logicals.shape == (24, 288)
    assert logical_metadata["basis_is_complete_quotient_basis"]


def test_translation_orbits_preserve_code_and_span_full_dual() -> None:
    _, _, code, _, logicals, _ = target()
    proof = EXP035.translation_orbit_proof(code, logicals)
    assert proof["translation_group"] == "Z_12 x Z_6"
    assert proof["translation_count"] == 72
    assert proof["every_translation_is_bijection"]
    assert proof["every_translation_preserves_stabilizer_rowspace_bitset"]
    assert proof["every_translation_preserves_logical_pairing_numpy"]
    assert proof["logical_pairing_rank_numpy"] == 24
    assert proof["logical_pairing_rank_bitset"] == 24
    assert proof["selected_orbit_union_rank_numpy"] == 24
    assert proof["selected_orbit_union_rank_bitset"] == 24
    assert proof["selected_orbits_span_full_dual"]
    assert [
        representative["orbit_quotient_rank_numpy"]
        for representative in proof["representatives"]
    ] == [8, 6, 8, 8]
    assert [
        representative["orbit_quotient_rank_bitset"]
        for representative in proof["representatives"]
    ] == [8, 6, 8, 8]
    for shift in proof["stabilizer_translation_checks"]:
        assert shift["permutation_is_bijection"]
        assert shift["maps_stabilizer_rowspace_to_itself_bitset"]
        assert shift["logical_pairing_matrix_preserved_numpy"]
    for representative in proof["representatives"]:
        assert len(representative["translations"]) == 72
        for translated in representative["translations"]:
            assert translated["numpy_and_bitset_coordinates_agree"]
            assert translated["residual_is_in_stabilizer_rowspace_bitset"]
            assert translated["translated_logical_is_centralizer_bitset"]
            assert translated["symplectic_pairing_transport_identity_numpy"]
            assert translated["symplectic_weight_preserved"]


def test_weight_twelve_witness_is_independently_nontrivial() -> None:
    _, _, code, _, _, _ = target()
    witness = EXP035.witness_with_crosschecks(code, EXP035.default_witness(code))
    assert witness["support"] == EXP035.DEFAULT_WITNESS_X_SUPPORT
    assert witness["symplectic_weight"] == 12
    assert witness["symplectic_weight_numpy"] == 12
    assert witness["commutes_with_all_stabilizers_bitset"]
    assert witness["commutes_with_all_stabilizers_numpy"]
    assert not witness["in_stabilizer_row_space_bitset"]
    assert not witness["in_stabilizer_row_space_numpy"]
    assert witness["nontrivial_logical_bitset"]
    assert witness["nontrivial_logical_numpy"]
    assert witness["independent_paths_agree"]


def test_negative_controls_reject_noncentralizer_and_stabilizer() -> None:
    _, _, code, _, _, _ = target()
    stabilizer = code.H[0].copy()
    stabilizer_check = EXP035.verify_witness_bitset(code, stabilizer)
    assert stabilizer_check["commutes_with_all_stabilizers_bitset"]
    assert stabilizer_check["in_stabilizer_row_space_bitset"]
    assert not stabilizer_check["nontrivial_logical_bitset"]

    single_qubit_x = np.zeros(2 * code.n, dtype=np.uint8)
    single_qubit_x[0] = 1
    noncentralizer = EXP035.verify_witness_bitset(code, single_qubit_x)
    assert not noncentralizer["commutes_with_all_stabilizers_bitset"]
    assert not noncentralizer["nontrivial_logical_bitset"]


def test_small_weight_exact_enumeration_certifies_lower_bound_six() -> None:
    _, _, code, _, _, _ = target()
    proof = EXP035.small_weight_centralizer_exclusion(code)
    assert proof["atom_count"] == 432
    assert proof["distinct_qubit_pair_count"] == 92_664
    assert proof["distinct_qubit_triple_count"] == 13_158_288
    assert proof["zero_syndrome_witnesses"] == {
        "weight_1": None,
        "weight_2": None,
        "weight_3": None,
        "weight_4": None,
        "weight_5": None,
    }
    assert proof["no_centralizer_of_weight_at_most_5"]
    assert proof["certified_distance_lower_bound"] == 6


def test_timeout_or_missing_sector_cannot_route_canonical(tmp_path: Path, monkeypatch) -> None:
    _, _, code, metadata, logicals, logical_metadata = target()
    monkeypatch.setattr(EXP035, "STATE_DIR", tmp_path / "state")
    monkeypatch.setattr(EXP035, "CANONICAL", tmp_path / "canonical.json")
    monkeypatch.setattr(EXP035, "PARTIAL_SUMMARY", tmp_path / "partial.json")
    canonical_report = tmp_path / "canonical_report.json"
    partial_report = tmp_path / "partial_report.json"
    monkeypatch.setattr(EXP035, "REPORT", canonical_report)
    monkeypatch.setattr(EXP035, "PARTIAL_REPORT", partial_report)
    result = EXP035.assemble_result(
        code,
        metadata,
        logicals,
        logical_metadata,
        records=[],
        missing=list(range(len(EXP035.ORBIT_REPRESENTATIVE_MASKS))),
        unresolved=[],
        base_seed=EXP035.DEFAULT_SEED,
        time_limit_s=10.0,
        workers=1,
        memory_mb=512,
    )
    route, path = EXP035.route_result(result)
    assert route == "partial"
    assert path != EXP035.CANONICAL
    assert result["bounds"]["distance"] is None
    assert result["bounds"]["certified_lower_bound"] == 7
    assert result["bounds"]["certified_upper_bound"] == 12
    assert not result["bounds"]["CERTIFIED_EXACT"]
    assert result["resume_command"] is not None
    result["artifact_route"] = route
    result["artifact_path"] = str(path)
    EXP035.write_concise_report(result, [])
    assert partial_report.exists()
    assert not canonical_report.exists()


def test_weight_six_classification_target_and_counterexample_control() -> None:
    _, _, code, _, _, _ = target()
    weight_six = EXP035.weight_six_zero_syndrome_classification(code)
    assert weight_six["complete"]
    assert weight_six["no_weight6_logical"]
    assert weight_six["distinct_weight_six_zero_syndrome_vectors"] == 72
    assert weight_six["stabilizer_members"] == 72
    # Ten 3+3 splits per weight-6 element is an exact combinatorial identity.
    assert weight_six["qubit_disjoint_equal_syndrome_pairs"] == 10 * 72
    assert weight_six["certified_lower_bound_with_weight5_exclusion"] == 7

    # Control with a certified weight-6 logical: phase2_58 is the EXP-027
    # [[72,4,6]] reversal with exact d=6, so the same classification must
    # exhibit at least one verified nonmember of weight exactly 6.
    control_row = None
    with EXP035.CATALOGUE.open(encoding="utf-8") as stream:
        for line in stream:
            candidate = json.loads(line)
            if candidate.get("code_id") == "phase2_58":
                control_row = candidate
                break
    assert control_row is not None
    control = EXP035.rebuild_target_independently(
        control_row
    )
    classified = EXP035.weight_six_zero_syndrome_classification(control)
    assert classified["complete"]
    assert not classified["no_weight6_logical"]
    assert classified["certified_lower_bound_with_weight5_exclusion"] == 6
    for logical in classified["nonmember_logicals"]:
        assert logical["symplectic_weight"] == 6
        assert logical["nontrivial_logical_bitset"]
        assert logical["nontrivial_logical_numpy"]
        assert logical["independent_paths_agree"]


def test_verified_weight6_counterexample_routes_to_certificate_not_canonical(
    monkeypatch,
) -> None:
    artifact = json.loads(EXP035.PARTIAL_SUMMARY.read_text())
    shaped = json.loads(json.dumps(artifact))
    fabricated = dict(shaped["witness"])
    fabricated["symplectic_weight"] = 6
    fabricated["symplectic_weight_numpy"] = 6
    shaped["unexpected_low_weight_logical"] = fabricated
    shaped["bounds"] = {
        "certified_lower_bound": 6,
        "certified_upper_bound": 6,
        "distance": 6,
        "CERTIFIED_EXACT": True,
    }
    # Scalar mutation of the genuine weight-12 witness must FAIL router
    # reverification (fresh algebra recomputes weight 12) and quarantine.
    route, path = EXP035.route_result(shaped)
    assert route == "quarantine"
    assert path == EXP035.QUARANTINE

    # Dispatch mechanics for a candidate that passes reverification: no
    # genuine target counterexample can exist (weight-6 classification proves
    # none), so acceptance is exercised by stubbing only the reverifier while
    # the rejection direction above uses the real algebra.
    monkeypatch.setattr(
        EXP035, "reverify_low_weight_counterexample", lambda candidate: True
    )
    route, path = EXP035.route_result(shaped)
    assert route == "counterexample"
    assert path == EXP035.COUNTEREXAMPLE
    assert path != EXP035.CANONICAL

    corrupted = json.loads(json.dumps(shaped))
    corrupted["rebuild"]["k"] = 11
    route, path = EXP035.route_result(corrupted)
    assert route == "quarantine"
    assert path == EXP035.QUARANTINE


def test_pysat_cnf_backend_matches_algebra_on_decisive_control() -> None:
    """Both backends must agree on the phase2_58 control where d=6 exactly."""

    control_row = None
    with EXP035.CATALOGUE.open(encoding="utf-8") as stream:
        for line in stream:
            candidate = json.loads(line)
            if candidate.get("code_id") == "phase2_58":
                control_row = candidate
                break
    assert control_row is not None
    control = EXP035.rebuild_target_independently(control_row)
    logicals = control.logical_basis()
    assert logicals.shape[0] == 8  # k = 4

    from pysat.solvers import Solver

    sat_sector = None
    for sector in range(logicals.shape[0]):
        detector = EXP035.lambda_swap(logicals[sector : sector + 1])[0]
        # Weight cap 5 must be UNSAT for every sector: d = 6 exactly, so no
        # logical of weight <= 5 exists at all.
        cnf5 = EXP035.build_sector_cnf(control, detector, 5)
        with Solver(name="cadical195", bootstrap_with=cnf5.clauses) as engine:
            assert engine.solve() is False
        if sat_sector is None:
            cnf6 = EXP035.build_sector_cnf(control, detector, 6)
            with Solver(name="cadical195", bootstrap_with=cnf6.clauses) as engine:
                if engine.solve():
                    model = engine.get_model()
                    assignment = {abs(l): l > 0 for l in model}
                    n = control.n
                    vector = np.asarray(
                        [1 if assignment.get(i + 1, False) else 0 for i in range(2 * n)],
                        dtype=np.uint8,
                    )
                    checked = EXP035.witness_with_crosschecks(control, vector)
                    assert checked["symplectic_weight"] <= 6
                    assert checked["commutes_with_all_stabilizers_bitset"]
                    assert checked["commutes_with_all_stabilizers_numpy"]
                    # pairing 1 with this sector => outside S
                    assert checked["nontrivial_logical_bitset"]
                    sat_sector = sector
    # d = 6 means some class has a weight-6 member; its pairing is nonzero
    # against at least one dual functional, so at least one sector is SAT.
    assert sat_sector is not None


def test_only_complete_matching_infeasible_orbit_cover_routes_canonical(
    tmp_path: Path, monkeypatch
) -> None:
    _, _, code, metadata, logicals, logical_metadata = target()
    state = tmp_path / "state"
    monkeypatch.setattr(EXP035, "STATE_DIR", state)
    monkeypatch.setattr(EXP035, "CANONICAL", tmp_path / "canonical.json")
    monkeypatch.setattr(EXP035, "PARTIAL_SUMMARY", tmp_path / "partial.json")
    records = []
    for sector in range(len(EXP035.ORBIT_REPRESENTATIVE_MASKS)):
        trace = state / "traces" / f"sector_{sector:02d}.txt"
        trace.parent.mkdir(parents=True, exist_ok=True)
        trace.write_text("status: INFEASIBLE\n")
        record = {
            **EXP035.sector_identity(
                code,
                logicals,
                sector,
                EXP035.DEFAULT_SEED,
                EXP035.LOWER_EXCLUSION_WEIGHT,
            ),
            "solver": {
                "name": "OR-Tools CP-SAT",
                "status": "INFEASIBLE",
                "proof_complete_for_sector": True,
                "wall_time_s": 0.1,
                "trace_path": str(trace.relative_to(ROOT))
                if trace.is_relative_to(ROOT)
                else str(trace),
                "trace_sha256": EXP035.sha256_bytes(trace.read_bytes()),
            },
            "solution": None,
        }
        records.append(record)
    result = EXP035.assemble_result(
        code,
        metadata,
        logicals,
        logical_metadata,
        records=records,
        missing=[],
        unresolved=[],
        base_seed=EXP035.DEFAULT_SEED,
        time_limit_s=10.0,
        workers=1,
        memory_mb=512,
    )
    route, path = EXP035.route_result(result)
    assert route == "canonical"
    assert path == EXP035.CANONICAL
    assert result["bounds"] == {
        "certified_lower_bound": 12,
        "certified_upper_bound": 12,
        "distance": 12,
        "CERTIFIED_EXACT": True,
    }


def test_committed_partial_artifact_claims_are_machine_checked() -> None:
    assert EXP035.PARTIAL_SUMMARY.exists()
    artifact = json.loads(EXP035.PARTIAL_SUMMARY.read_text())
    _, _, code, metadata, logicals, logical_metadata = target()
    assert artifact["schema"] == "pbb-exact-distance-certificate-v3"
    assert artifact["artifact_route"] == "partial"
    if EXP035.CANONICAL.exists():
        # The canonical certificate supersedes the frozen partial snapshot;
        # the partial file is retained as historical evidence only.  Its
        # sector-progress fields describe an earlier campaign state, so only
        # the state-independent claims are re-checked here.
        canonical = json.loads(EXP035.CANONICAL.read_text())
        assert canonical["schema"] == artifact["schema"]
        assert canonical["source"] == artifact["source"]
        assert canonical["bounds"]["CERTIFIED_EXACT"]
        assert artifact["bounds"]["certified_lower_bound"] <= (
            canonical["bounds"]["certified_lower_bound"]
        )
        assert artifact["bounds"]["certified_upper_bound"] == (
            canonical["bounds"]["certified_upper_bound"]
        )
        return
    assert artifact["source"] == metadata["source"]
    assert artifact["rebuild"] == metadata["rebuild"]
    assert artifact["logical_quotient"] == logical_metadata
    assert artifact["translation_orbit_reduction"] == (
        EXP035.translation_orbit_proof(code, logicals)
    )
    assert artifact["small_weight_exclusion"] == (
        EXP035.small_weight_centralizer_exclusion(code)
    )
    assert artifact["weight_six_classification"] == (
        EXP035.weight_six_zero_syndrome_classification(code)
    )
    assert artifact["weight_six_classification"]["no_weight6_logical"]
    assert artifact["weight_six_classification"][
        "matches_stored_pure_z_weight6_rows"
    ]
    assert artifact["unexpected_low_weight_logical"] is None
    witness_vector = np.asarray(
        artifact["witness"]["x"] + artifact["witness"]["z"], dtype=np.uint8
    )
    checked_witness = EXP035.witness_with_crosschecks(code, witness_vector)
    assert checked_witness["support"] == artifact["witness"]["support"]
    assert checked_witness["nontrivial_logical_bitset"]
    assert checked_witness["nontrivial_logical_numpy"]
    problem = artifact["distance_problem"]
    covered = problem["covered_infeasible_representatives"]
    covered_coordinates = [
        int(translation["quotient_coefficient_mask_hex"], 16)
        for representative in covered
        for translation in artifact["translation_orbit_reduction"][
            "representatives"
        ][representative]["translations"]
    ]
    recomputed_rank = EXP035.rank_bitset(covered_coordinates, logicals.shape[0])
    assert problem["covered_dual_rank_numpy"] == recomputed_rank
    assert problem["covered_dual_rank_bitset"] == recomputed_rank
    pending = sorted(
        problem["missing_representatives"]
        + problem["unresolved_representatives"]
    )
    assert sorted(covered + pending) == [0, 1, 2, 3]
    statuses = problem["representative_statuses"]
    for sector in covered:
        assert statuses[str(sector)] in {"INFEASIBLE", "UNSAT"}
    assert artifact["bounds"]["certified_lower_bound"] == 7
    assert artifact["bounds"]["certified_upper_bound"] == 12
    assert artifact["bounds"]["distance"] is None
    assert not artifact["bounds"]["CERTIFIED_EXACT"]
    resume = artifact["resume_command"]
    assert "exp035_pbb_exact_distance.py run" in resume
    for sector in pending:
        assert f" {sector}" in resume.split("--sectors", 1)[1].split("--seed", 1)[0]


def test_committed_certificate_claims_are_machine_checked_if_present() -> None:
    if not EXP035.CANONICAL.exists():
        return
    certificate = json.loads(EXP035.CANONICAL.read_text())
    _, _, code, metadata, logicals, logical_metadata = target()
    assert certificate["schema"] == "pbb-exact-distance-certificate-v3"
    assert certificate["source"]["terms"] == EXP035.EXPECTED_TERMS
    assert certificate["source"]["catalogue_file_sha256"] == metadata["source"][
        "catalogue_file_sha256"
    ]
    assert certificate["rebuild"]["matrix_sha256"] == EXP035.matrix_sha256(code.H)
    assert certificate["logical_quotient"]["basis_sha256"] == logical_metadata[
        "basis_sha256"
    ]
    assert certificate["translation_orbit_reduction"] == (
        EXP035.translation_orbit_proof(code, logicals)
    )
    assert certificate["small_weight_exclusion"] == (
        EXP035.small_weight_centralizer_exclusion(code)
    )
    vector = np.array(certificate["witness"]["x"] + certificate["witness"]["z"], dtype=np.uint8)
    witness = EXP035.witness_with_crosschecks(code, vector)
    assert witness["support"] == certificate["witness"]["support"]
    assert witness["symplectic_weight"] == certificate["witness"]["symplectic_weight"] == 12
    assert witness["nontrivial_logical_bitset"]
    assert witness["nontrivial_logical_numpy"]
    problem = certificate["distance_problem"]
    assert problem["covered_infeasible_representatives"] == list(
        range(len(EXP035.ORBIT_REPRESENTATIVE_MASKS))
    )
    assert problem["covered_dual_rank_numpy"] == 24
    assert problem["covered_dual_rank_bitset"] == 24
    assert problem["missing_representatives"] == []
    assert problem["unresolved_representatives"] == []
    statuses = set(problem["representative_statuses"].values())
    assert statuses and statuses <= {"INFEASIBLE", "UNSAT"}
    weight_six = certificate["weight_six_classification"]
    assert weight_six["complete"] and weight_six["no_weight6_logical"]
    assert certificate["unexpected_low_weight_logical"] is None
    assert certificate["bounds"] == {
        "certified_lower_bound": 12,
        "certified_upper_bound": 12,
        "distance": 12,
        "CERTIFIED_EXACT": True,
    }


def _write_sector_record(
    state_dir: Path,
    identity: dict,
    *,
    backend: str,
    status: str,
    proof_complete: bool,
    solution: dict | None = None,
) -> dict:
    """Persist a synthetic record satisfying the identity/trace binding."""

    sector = identity["sector"]
    suffix = "_pysat" if backend.startswith("PySAT") else ""
    trace = state_dir / f"trace_{sector:02d}{suffix}.txt"
    trace.write_text(f"synthetic trace {backend} {status}\n", encoding="utf-8")
    record = dict(identity)
    record["solver"] = {
        "name": backend,
        "status": status,
        "proof_complete_for_sector": proof_complete,
        "wall_time_s": 0.0,
        "trace_path": os.path.relpath(trace, EXP035.ROOT),
        "trace_sha256": EXP035.sha256_bytes(trace.read_bytes()),
    }
    if solution is not None:
        record["solution"] = solution
    path = (
        EXP035.pysat_sector_path(sector)
        if suffix
        else EXP035.sector_path(sector)
    )
    path.write_text(json.dumps(record), encoding="utf-8")
    return record


def _isolated_state_dir(tmp_path: Path, monkeypatch) -> Path:
    state = tmp_path / "state"
    state.mkdir()
    monkeypatch.setattr(EXP035, "STATE_DIR", state)
    return state


def test_record_proof_complete_accepts_only_frozen_backend_status_pairs() -> None:
    ok_cpsat = {"solver": {"name": "OR-Tools CP-SAT", "status": "INFEASIBLE", "proof_complete_for_sector": True}}
    ok_pysat = {"solver": {"name": "PySAT cadical195", "status": "UNSAT", "proof_complete_for_sector": True}}
    unclaimed = {"solver": {"name": "OR-Tools CP-SAT", "status": "UNKNOWN", "proof_complete_for_sector": False}}
    assert EXP035.record_proof_complete(0, ok_cpsat)
    assert EXP035.record_proof_complete(1, ok_pysat)
    assert not EXP035.record_proof_complete(2, unclaimed)
    import pytest

    with pytest.raises(RuntimeError, match="frozen backend/"):
        EXP035.record_proof_complete(
            3,
            {"solver": {"name": "PySAT cadical195", "status": "SAT", "proof_complete_for_sector": True}},
        )
    with pytest.raises(RuntimeError, match="frozen backend/"):
        EXP035.record_proof_complete(
            0,
            {"solver": {"name": "OR-Tools CP-SAT", "status": "UNKNOWN", "proof_complete_for_sector": True}},
        )


def test_collect_covers_sector_from_independent_pysat_unsat(tmp_path, monkeypatch) -> None:
    _, _, code, _, logicals, _ = target()
    state = _isolated_state_dir(tmp_path, monkeypatch)
    base_seed, cap = 777, EXP035.LOWER_EXCLUSION_WEIGHT
    identity = EXP035.sector_identity(code, logicals, 1, base_seed, cap)
    _write_sector_record(
        state, identity, backend="OR-Tools CP-SAT", status="UNKNOWN", proof_complete=False
    )
    _write_sector_record(
        state, identity, backend="PySAT cadical195", status="UNSAT", proof_complete=True
    )
    records, missing, unresolved = EXP035.collect_sector_records(
        code, logicals, base_seed=base_seed, weight_cap=cap
    )
    assert [record["sector"] for record in records] == [1]
    assert records[0]["solver"]["name"].startswith("PySAT")
    assert missing == [0, 2, 3]
    assert unresolved == []


def test_collect_fail_stops_on_unverifiable_claimed_solution(tmp_path, monkeypatch) -> None:
    _, _, code, _, logicals, _ = target()
    state = _isolated_state_dir(tmp_path, monkeypatch)
    base_seed, cap = 778, EXP035.LOWER_EXCLUSION_WEIGHT
    identity = EXP035.sector_identity(code, logicals, 0, base_seed, cap)
    witness = EXP035.witness_with_crosschecks(
        code, EXP035.default_witness(code)
    )
    # Genuine weight-12 logical: exceeds the weight-11 cap, so a SAT record
    # claiming it as a sector solution must fail-stop, not rank as evidence.
    _write_sector_record(
        state,
        identity,
        backend="PySAT cadical195",
        status="SAT",
        proof_complete=False,
        solution={"vector": witness["x"] + witness["z"]},
    )
    import pytest

    with pytest.raises(RuntimeError, match="failed reverification"):
        EXP035.collect_sector_records(
            code, logicals, base_seed=base_seed, weight_cap=cap
        )


def test_collect_fail_stops_on_cross_backend_contradiction(tmp_path, monkeypatch) -> None:
    _, _, code, _, logicals, _ = target()
    state = _isolated_state_dir(tmp_path, monkeypatch)
    base_seed, cap = 779, EXP035.LOWER_EXCLUSION_WEIGHT
    identity = EXP035.sector_identity(code, logicals, 2, base_seed, cap)
    witness = EXP035.witness_with_crosschecks(
        code, EXP035.default_witness(code)
    )
    _write_sector_record(
        state, identity, backend="OR-Tools CP-SAT", status="INFEASIBLE", proof_complete=True
    )
    _write_sector_record(
        state,
        identity,
        backend="PySAT cadical195",
        status="SAT",
        proof_complete=False,
        solution={"vector": witness["x"] + witness["z"]},
    )
    import pytest

    with pytest.raises(RuntimeError, match="backend contradiction"):
        EXP035.collect_sector_records(
            code, logicals, base_seed=base_seed, weight_cap=cap
        )
