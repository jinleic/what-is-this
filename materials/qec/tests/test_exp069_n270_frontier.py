"""Regression checks for EXP-069's canonical n=270 frontier."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp069", ROOT / "experiments" / "exp069_n270_frontier.py"
)
E69 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E69
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E69)


def test_noncyclic_presentation_transport_is_exact_and_deduplicated() -> None:
    audit = E69.presentation_transport_audit()

    assert audit["source_group"] == "Z_15 x Z_9"
    assert audit["target_group"] == "Z_45 x Z_3"
    assert audit["group_order"] == 135
    assert audit["all_coordinates_bijective"] is True
    assert audit["inverse_exact"] is True
    assert audit["homomorphism_exhaustive"] is True
    assert audit["candidate_classes"] == {"15x9": 5024, "45x3": 5024}
    assert audit["represented_pairs"] == {"15x9": 192834, "45x3": 192834}
    assert audit["screened_presentation"] == [15, 9]
    assert audit["excluded_duplicate_presentation"] == [45, 3]


def test_liang_n270_rows_are_bound_to_local_rectangular_codes() -> None:
    for key in E69.LIANG_TARGETS:
        binding = E69.published_transport(key)

        assert binding["quotient_order"] == 135
        assert binding["lattice_relations_vanish"] is True
        assert binding["quotient_map_surjective"] is True
        assert binding["target_45x3"]["k"] == 8
        assert binding["local_15x9"]["k"] == 8
        assert binding["matrix_transport"]["HX_exact"] is True
        assert binding["matrix_transport"]["HZ_exact"] is True
        assert binding["matrix_transport"]["qubit_permutation_bijective"] is True


def test_deep_reduction_recovers_both_published_weight20_witnesses() -> None:
    for target in E69.LIANG_TARGETS.values():
        local = target["local_15x9"]
        record = {
            "ell": 15,
            "m": 9,
            "n": 270,
            "A": local["A"],
            "B": local["B"],
            "threshold": 20,
        }
        result = E69.deepen_record(record, tries=2_000)

        assert result["bound"] == 20
        assert len(result["witness_support"]) == 20
        assert result["verification"]["weight"] == 20
        assert result["verification"]["z_valid_numpy"] is True
        assert result["verification"]["z_valid_bitset"] is True
        stored = E69.verify_witness(record, target["witness_support"])
        assert stored["weight"] == 20


def test_deep_resolution_is_fail_closed() -> None:
    record = {
        "verdict": "solver_required",
        "threshold": 18,
        "solver_calls": 0,
    }
    dominated = E69.finalize_deep_record(
        record,
        {"bound": 18, "witness_support": [1, 2], "wall_time_s": 0.1},
    )
    assert dominated["verdict"] == "dominated_by_deep_reduction"
    assert dominated["witness_bound"] == 18
    assert dominated["witness_support"] == [1, 2]

    open_record = E69.finalize_deep_record(
        record,
        {"bound": 20, "witness_support": [3, 4], "wall_time_s": 0.1},
    )
    assert open_record["verdict"] == "undecided"
    assert open_record["witness_bound"] == 20


def test_exact_fallback_persists_a_physical_z_witness(monkeypatch) -> None:
    record = {
        "ell": 3,
        "m": 3,
        "n": 18,
        "A": [[0, 0], [0, 1], [0, 2]],
        "B": [[0, 0], [1, 0], [0, 1]],
        "k_parent": 4,
        "orbit": 1,
        "threshold": 2,
        "threshold_source": "small-control",
        "verdict": "solver_required",
        "solver_calls": 0,
    }
    monkeypatch.setattr(
        E69.E55,
        "_cdcl_witness_bound",
        lambda *_args: {
            "status": "UNSAT",
            "cnf_sha256": "control",
            "encoding_version": "control",
            "solver": {"name": "control"},
        },
    )

    resolved = E69.solve_residual_record(record, time_limit_s=30.0)

    assert resolved["verdict"] == "dominated_by_exact_witness"
    assert resolved["screen_decided"] is True
    assert resolved["witness_bound"] == 2
    assert len(resolved["witness_support"]) == 2
    assert E69.verify_witness(
        record, resolved["witness_support"]
    )["weight"] == 2



def test_k8_resolution_waits_for_a_bound_n270_reference(monkeypatch) -> None:
    record = {
        "k_parent": 8,
        "threshold": 18,
        "threshold_source": "EXP-067 [[234,8,18]]",
        "verdict": "solver_required",
        "solver_calls": 0,
    }
    monkeypatch.setattr(
        E69,
        "deepen_record",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("premature deep search")
        ),
    )
    monkeypatch.setattr(
        E69,
        "solve_residual_record",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("premature exact search")
        ),
    )

    deferred = E69.resolve_candidate_record(record, tries=2_000)

    assert deferred["verdict"] == "undecided"
    assert deferred["reference_pending"]["required_threshold"] == 20
    assert deferred["solver_calls"] == 0

    hard = E69.resolve_candidate_record(
        {
            **record,
            "k_parent": 12,
            "threshold": 12,
            "threshold_source": "EXP-037 [[180,12,12]]",
        },
        tries=2_000,
    )
    assert hard["verdict"] == "undecided"
    assert hard["exact_pending"]["resource_reason"] == (
        "serial per-sector fallback exceeded the bounded screen budget"
    )

def test_initial_screen_is_serial_checkpointed_and_resumable(
    tmp_path: Path, monkeypatch
) -> None:
    candidates = [
        {"A": [(0, 0)], "B": [(0, 0)], "k_parent": 8, "orbit": 1},
        {"A": [(0, 0)], "B": [(1, 0)], "k_parent": 10, "orbit": 2},
    ]
    calls: list[int] = []

    monkeypatch.setattr(E69, "PARTIAL_DIR", tmp_path)
    monkeypatch.setattr(
        E69.E55, "enumerate_candidates", lambda *_args: candidates
    )

    def fake_screen(task):
        calls.append(task[2]["k_parent"])
        return {**task[2], "verdict": "solver_required", "solver_calls": 0}

    monkeypatch.setattr(E69.E55, "_screen_one", fake_screen)
    payload = E69.screen_initial((15, 9), checkpoint_every=1)

    assert calls == [8, 10]
    assert payload["status"] == "COMPLETE"
    assert payload["next_index"] == 2
    assert payload["resource_policy"] == {
        "processes": 1,
        "threads": 1,
        "external_solvers": False,
    }
    assert len(payload["records"]) == 2

    cached = E69.screen_initial((15, 9), checkpoint_every=1)
    assert cached == payload
    assert calls == [8, 10]
