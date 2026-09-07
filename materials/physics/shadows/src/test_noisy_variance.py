# -*- coding: utf-8 -*-
"""Exact regression contract for the noisy global-shadow variance campaign.

Run: python3 src/test_noisy_variance.py  (or pytest src/test_noisy_variance.py)
"""

import contextlib
import io
import sys
import json
import tempfile
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import shadows_noisy_variance as N  # noqa: E402


def test_paper_threshold_and_ratio_anchors():
    assert [N.rho_s(d) for d in N.DIMENSIONS] == [F(10, 9), F(27, 25), F(85, 81)]
    assert [N.beta_star(d) for d in N.DIMENSIONS] == [F(13, 4), F(128, 23), F(752, 77)]
    assert [N.depolarizing_p_star(d) for d in N.DIMENSIONS] == [F(3, 4), F(15, 23), F(45, 77)]
    assert N.rho_l(4, N.depolarizing_beta(4, F(1, 2))) == F(50, 27)
    for d in N.DIMENSIONS:
        assert N.rho_l(d, N.beta_star(d)) == 2


def test_anchor_invariants_are_exact():
    assert N.anchor_invariants(4, "maximally_mixed_pauli_z") == (F(4), F(1), F(0))
    assert N.anchor_invariants(4, "ghz_projector") == (F(3, 4), F(9, 16), F(3, 4))
    assert N.anchor_invariants(4, "ghz_pauli_x") == (F(4), F(1), F(1))


def test_grid_is_deduplicated_and_stays_in_invertible_scope():
    rows = N.build_rows()
    assert len(rows) == 132
    assert sum(row["noise_model"] == "depolarizing" for row in rows) == 60
    assert sum(row["noise_model"] == "amplitude_damping" for row in rows) == 72
    keys = {(row["noise_model"], row["d"], row["p"], row["anchor"]) for row in rows}
    assert len(keys) == len(rows)
    assert all(F(1) < row["beta"] <= row["d"] for row in rows)


def test_grid_matches_every_preregistered_key_and_crossing_role():
    rows = N.build_rows()
    dep_base = {F(1), F(3, 4), F(1, 2), F(1, 4), F(1, 10), F(1, 100)}
    dep = {
        4: dep_base,
        8: dep_base | {F(15, 23)},
        16: dep_base | {F(45, 77)},
    }
    amplitude = {
        F(1), F(9, 10), F(4, 5), F(3, 4),
        F(1, 2), F(1, 4), F(1, 10), F(1, 100),
    }
    expected_keys = set()
    for d in (4, 8, 16):
        for p in dep[d]:
            for anchor in (
                "maximally_mixed_pauli_z", "ghz_projector", "ghz_pauli_x"
            ):
                expected_keys.add(("depolarizing", d, p, anchor))
        for p in amplitude:
            for anchor in (
                "maximally_mixed_pauli_z", "ghz_projector", "ghz_pauli_x"
            ):
                expected_keys.add(("amplitude_damping", d, p, anchor))
    actual_keys = {
        (row["noise_model"], row["d"], row["p"], row["anchor"]) for row in rows
    }
    assert actual_keys == expected_keys
    crossings = {4: F(3, 4), 8: F(15, 23), 16: F(45, 77)}
    for row in rows:
        expected_role = (
            row["noise_model"] == "depolarizing"
            and row["p"] == crossings[row["d"]]
        )
        assert ("crossing_point" in row["p_roles"]) == expected_role
    grid = N.grid_summary(rows)
    assert grid["exact_grid_ok"]
    assert grid["crossing_roles_ok"]


def test_runner_binds_loaded_source_and_requires_revision_marker():
    assert N.MODULE_BYTES_AT_IMPORT == N.MODULE_PATH.read_bytes()
    assert N.prereg_revision_from_bytes(N.PREREG_PATH.read_bytes()) == 7
    assert N.prereg_revision_from_bytes(b"# no revision marker\n") == -1


def test_direct_moments_match_independent_ratio_paths_and_never_flip():
    for row in N.build_rows():
        assert row["second_moment_ratio"] == row["criterion_ratio"]
        assert row["variance_ratio"] == (
            row["variance_unitary"] / row["variance_orthogonal"]
        )
        assert row["variance_ratio"] == row["eq44_variance_ratio"]
        assert row["variance_ratio"] >= row["second_moment_ratio"] > 1
        assert row["variance_orthogonal"] > 0
        assert row["variance_unitary"] > 0


def test_terminal_mapping_uses_control_plane_verdicts():
    assert N.classify_verdict(False, True, True, False) == (
        "INCONCLUSIVE-IDENTITIES-INVALID", "FROZEN-INCONCLUSIVE"
    )
    assert N.classify_verdict(True, True, True, True) == (
        "FLIP-FOUND", "FROZEN-NEGATIVE"
    )
    assert N.classify_verdict(True, True, True, False) == (
        "NO-FLIP-CONFIRMED-IN-SCOPE", "FROZEN-CERTIFIED"
    )


def test_runner_accepts_control_plane_manifest_and_binds_prereg_hash():
    run_dir = Path(tempfile.mkdtemp(prefix="shadows-manifest-test-"))
    try:
        manifest = {
            "run_id": run_dir.name,
            "target": "physics/shadows",
            "gate": "noisy-variance-r7",
            "agent": "test",
            "created_utc": "2026-09-04T00:00:00Z",
            "prereg_sha256": N.sha256_file(N.PREREG_PATH),
            "status": "RUNNING",
        }
        (run_dir / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        assert N.load_manifest(run_dir) == manifest
        manifest["prereg_sha256"] = "0" * 64
        (run_dir / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        try:
            N.load_manifest(run_dir)
        except N.RunnerRefused:
            pass
        else:
            raise AssertionError("runner accepted a preregistration hash mismatch")
        (run_dir / "manifest.json").write_text("[]", encoding="utf-8")
        try:
            N.load_manifest(run_dir)
        except N.RunnerRefused:
            pass
        else:
            raise AssertionError("runner accepted a non-object manifest")
    finally:
        # Workspace policy forbids deleting even test-temporary files. The OS
        # owns eventual cleanup of this private system-temporary directory.
        pass


def test_end_to_end_artifacts_are_bound_deterministic_and_no_overwrite():
    run_dir = Path(tempfile.mkdtemp(prefix="shadows-artifact-test-"))
    manifest = {
        "run_id": run_dir.name,
        "target": "physics/shadows",
        "gate": "noisy-variance-r7",
        "agent": "test",
        "created_utc": "2026-09-04T00:00:00Z",
        "prereg_sha256": N.sha256_file(N.PREREG_PATH),
        "status": "RUNNING",
    }
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    first = N.run_payloads(run_dir)[3]
    second = N.run_payloads(run_dir)[3]
    assert first == second
    with contextlib.redirect_stdout(io.StringIO()):
        assert N.write_artifacts(run_dir) == 0

    assert (
        run_dir / N.MODULE_SNAPSHOT_NAME
    ).read_bytes() == N.MODULE_BYTES_AT_IMPORT
    assert (
        run_dir / N.PREREG_SNAPSHOT_NAME
    ).read_bytes() == N.PREREG_PATH.read_bytes()
    results = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    assert results["snapshots"]["module"]["path"] == (
        "src/shadows_noisy_variance.py"
    )
    assert results["preregistration"]["path"] == "pre_statement.md"
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    inventory = json.loads(
        (run_dir / "inventory.json").read_text(encoding="utf-8")
    )
    assert results["campaign_manifest"] == manifest
    assert results["module_sha256"] == N.sha256_file(
        run_dir / N.MODULE_SNAPSHOT_NAME
    )
    assert summary["results_sha256"] == N.sha256_file(run_dir / "results.json")
    assert inventory["summary_sha256"] == N.sha256_file(run_dir / "summary.json")
    for artifact in inventory["artifacts"]:
        path = run_dir / artifact["name"]
        assert artifact["sha256"] == N.sha256_file(path)
        assert artifact["bytes"] == path.stat().st_size

    before = {
        path.name: (N.sha256_file(path), path.stat().st_size)
        for path in run_dir.iterdir()
    }
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            N.write_artifacts(run_dir)
    except N.RunnerRefused:
        pass
    else:
        raise AssertionError("runner overwrote an existing artifact set")
    after = {
        path.name: (N.sha256_file(path), path.stat().st_size)
        for path in run_dir.iterdir()
    }
    assert after == before

    refused_dir = Path(tempfile.mkdtemp(prefix="shadows-refusal-test-"))
    refused_manifest = dict(manifest)
    refused_manifest["run_id"] = refused_dir.name
    refused_manifest["prereg_sha256"] = "0" * 64
    (refused_dir / "manifest.json").write_text(
        json.dumps(refused_manifest), encoding="utf-8"
    )
    initial_names = {path.name for path in refused_dir.iterdir()}
    try:
        N.write_artifacts(refused_dir)
    except N.RunnerRefused:
        pass

    else:
        raise AssertionError("runner accepted invalid provenance")
    assert {path.name for path in refused_dir.iterdir()} == initial_names

def test_zero_denominator_is_a_documented_refusal_exit():
    original = N.write_artifacts
    N.write_artifacts = lambda _: 1 / 0
    stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(stderr):
            assert N.main(["--run-dir", "/not-used-by-test"]) == 2
    finally:
        N.write_artifacts = original
    assert stderr.getvalue().startswith(
        "REFUSED: invalid input or arithmetic domain:"
    )

if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
    print(f"PASS: {len(tests)} exact noisy-shadow regression tests")
