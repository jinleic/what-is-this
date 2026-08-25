"""Regression checks for EXP-067's rooted connected-cluster certificates."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp067", ROOT / "experiments" / "exp067_n234_connected_cluster.py"
)
E67 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E67
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E67)


def _synthetic_complete_run(target: str, presentation: str) -> dict:
    timeout = 3500.0
    paths = E67.input_paths(target, presentation)
    stdout = "-16\n"
    stderr = "# exhaustive connected-cluster control\n"
    return {
        "schema": E67.RUN_SCHEMA,
        "utc": "2026-08-24T00:00:00+00:00",
        "target": target,
        "presentation": presentation,
        "status": "COMPLETE",
        "returncode": 0,
        "wall_time_s": 1.0,
        "input": E67.input_identity(target, presentation),
        "input_files": E67.input_file_bindings(target, presentation),
        "invocation": [
            "/tmp/dist_m4ri",
            "method=2",
            f"finH={paths['HX'].resolve()}",
            f"finG={paths['HZ'].resolve()}",
            "start=0",
            "dmin=1",
            f"wmax={E67.CAP}",
            "smax=0",
            "debug=3",
            f"timeout={timeout:g}",
        ],
        "solver": {
            "name": "QEC-pages/dist-m4ri legacy STANDALONE method 2",
            "source": "https://github.com/QEC-pages/dist-m4ri",
            "source_commit": E67.DIST_M4RI_COMMIT,
            "binary_sha256": E67.DIST_M4RI_BINARY_SHA256,
        },
        "protocol": {
            "engine_variant": "legacy_single_thread_do_CC_dist",
            "method": 2,
            "root": 0,
            "dmin": 1,
            "wmax": E67.CAP,
            "smax": 0,
            "timeout_s": timeout,
        },
        "stdout": stdout,
        "stderr": stderr,
        "stdout_sha256": hashlib.sha256(stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(stderr.encode()).hexdigest(),
        "result": {
            "status": "NO_LOGICAL_THROUGH_CAP",
            "exhaustive_through": 16,
            "found_distance": None,
        },
    }


def test_liang_table_iii_rows_match_the_two_open_bundles_uniquely() -> None:
    row2 = E67.published_transport("bundle19")
    row1 = E67.published_transport("bundle22")

    assert row2["published_parameters"] == [234, 8, 18]
    assert row2["unique_maps"] == 1
    assert row2["map"]["x_image"] == [23, 2]
    assert row2["map"]["y_image"] == [34, 0]
    assert row1["unique_maps"] == 1
    assert row1["map"]["x_image"] == [2, 2]
    assert row1["map"]["y_image"] == [5, 0]
    for row in (row1, row2):
        assert row["matrix_transport"]["HX_rowspace_equal"] is True
        assert row["matrix_transport"]["HZ_rowspace_equal"] is True


def test_solver_build_is_source_and_binary_bound() -> None:
    build = E67.solver_build_binding()
    assert build["source_commit"] == E67.DIST_M4RI_COMMIT
    assert build["binary_sha256"] == E67.DIST_M4RI_BINARY_SHA256
    assert len(build["source_archive_sha256"]) == 64
    assert len(build["m4ri_source_archive_sha256"]) == 64


def test_two_root_presentations_are_exact_block_permutations() -> None:
    for key in E67.TARGETS:
        problem = E67.problem_of(key)
        hx, hz, lx = E67._presentation_matrices(problem, "original")
        sx, sz, sl = E67._presentation_matrices(problem, "block_swapped")
        block = problem["block"]
        permutation = np.r_[np.arange(block, 2 * block), np.arange(block)]

        assert np.array_equal(sx, hx[:, permutation])
        assert np.array_equal(sz, hz[:, permutation])
        assert np.array_equal(sl, lx[:, permutation])
        assert E67.rank_np(hx) == E67.rank_np(sx) == 113
        assert E67.rank_np(hz) == E67.rank_np(sz) == 113
        assert E67.rank_np(lx) == E67.rank_np(sl) == 8


def test_root_reduction_uses_zero_not_the_unsound_second_block_index() -> None:
    for key in E67.TARGETS:
        reduction = E67.root_reduction(key)
        assert reduction["valid"] is True
        assert reduction["translation_group_size"] == 117
        assert reduction["presentations"] == ["original", "block_swapped"]
        assert reduction["root_coordinate"] == 0
        assert "start=117 is not used" in reduction["dist_m4ri_ordering_guard"]
        assert "same logical class" in reduction["irreducible_subset_lemma"]
        assert "recursive branch" in reduction["enumeration_argument"]


def test_run_record_contract_rejects_partial_or_unbound_searches() -> None:
    complete = _synthetic_complete_run("bundle19", "original")
    E67.validate_run_record(complete, "bundle19", "original")

    partial = json.loads(json.dumps(complete))
    partial["result"]["exhaustive_through"] = 15
    partial["stdout"] = "-15\n"
    partial["stdout_sha256"] = hashlib.sha256(partial["stdout"].encode()).hexdigest()
    with pytest.raises(RuntimeError, match="invalid EXP-067 run record"):
        E67.validate_run_record(partial, "bundle19", "original")

    race_affected_format = json.loads(json.dumps(complete))
    race_affected_format["stdout"] = "17 0 0\n"
    race_affected_format["stdout_sha256"] = hashlib.sha256(
        race_affected_format["stdout"].encode()
    ).hexdigest()
    race_affected_format["result"] = {
        "dmin": 17,
        "dmax": 0,
        "rw_steps": 0,
    }
    with pytest.raises(RuntimeError, match="invalid EXP-067 run record"):
        E67.validate_run_record(race_affected_format, "bundle19", "original")

    wrong_commit = json.loads(json.dumps(complete))
    wrong_commit["solver"]["source_commit"] = "f" * 40
    with pytest.raises(RuntimeError, match="invalid EXP-067 run record"):
        E67.validate_run_record(wrong_commit, "bundle19", "original")

    wrong_binary = json.loads(json.dumps(complete))
    wrong_binary["solver"]["binary_sha256"] = "f" * 64
    with pytest.raises(RuntimeError, match="invalid EXP-067 run record"):
        E67.validate_run_record(wrong_binary, "bundle19", "original")

    wrong_input = json.loads(json.dumps(complete))
    wrong_input["input_files"]["HX"]["sha256"] = "f" * 64
    with pytest.raises(RuntimeError, match="invalid EXP-067 run record"):
        E67.validate_run_record(wrong_input, "bundle19", "original")


def test_weight_eighteen_witnesses_and_even_kernel_close_seventeen() -> None:
    for key in E67.TARGETS:
        problem = E67.problem_of(key)
        upper = E67.upper_certificate(key)
        parity = E67.E56.kernel_weight_parity_proof(problem, 17)
        duality = E67.E56.bb_duality_proof(problem)

        assert upper["weight"] == 18
        assert upper["z_valid_numpy"] is True
        assert upper["z_valid_bitset"] is True
        assert upper["x_valid_numpy"] is True
        assert upper["x_valid_bitset"] is True
        assert parity["valid"] is True
        assert parity["effective_even_cap"] == 16
        assert duality["valid"] is True


def test_persisted_exp067_certificates_rebuild_exactly() -> None:
    for key in E67.TARGETS:
        payload = json.loads(E67.certificate_path(key).read_text(encoding="utf-8"))
        E67.validate_exact_certificate_payload(payload)
        assert payload["verdict"] == {
            "classification": "CERTIFIED_EXACT",
            "d": 18,
            "d_X": 18,
            "d_Z": 18,
            "exact": True,
        }
