"""Focused contract tests for EXP-027's certified-bound bookkeeping."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp027", ROOT / "experiments" / "exp027_delta_audit.py"
)
exp027 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(exp027)


def test_catalogue_scope_and_rank_relation() -> None:
    rows = exp027.load_catalogue()
    items = [
        item
        for index, row in enumerate(rows)
        if (item := exp027.structural_item(row, index)) is not None
    ]
    assert len(rows) == 368
    assert len(items) == 155
    assert all(item["delta"] > 0 for item in items)
    assert all(item["parent_k"] == item["k"] + item["delta"] for item in items)
    assert len({item["label"] for item in items}) == len(items)


def test_catalogue_distance_is_never_promoted_to_exact() -> None:
    catalogue_exact_row = next(row for row in exp027.load_catalogue() if row["d_is_exact"])
    bound = exp027.catalogue_pbb_bound(catalogue_exact_row)
    assert bound["value"] == catalogue_exact_row["d"]
    assert bound["upper_bound"] == catalogue_exact_row["d"]
    assert bound["lower_bound"] == 1
    assert bound["exact"] is False
    assert bound["bound_direction"] == "UPPER_BOUND"
    assert bound["catalogue_d_is_exact_ignored"] is True


def test_certificate_and_bravyi_matching_use_matrices_without_promoting_labels() -> None:
    parent_certificates, pbb_certificates = exp027.load_certificate_indexes()
    assert len(parent_certificates) == 4
    assert pbb_certificates == {}
    references = exp027.bravyi_index()
    assert len(references) == len(exp027.BRAVYI_BB)
    for name in ("[[72,12,6]]", "[[108,8,10]]", "[[144,12,12]]"):
        spec = exp027.BRAVYI_BB[name]
        HX, HZ = exp027.build_bb(spec)
        fingerprint = exp027.matrix_fingerprint(HX, HZ)
        assert references[fingerprint]["exact"] is False
        assert references[fingerprint]["lower_bound"] == 1
        assert references[fingerprint]["upper_bound"] is None
        assert references[fingerprint]["lower_bound_certified"] is True
        assert references[fingerprint]["upper_bound_certified"] is False
        assert references[fingerprint]["literature_reported_value"] > 1
        assert references[fingerprint]["provenance"].startswith(
            "literature_reported_uncertified"
        )


def test_partial_parent_solver_result_is_unsettled() -> None:
    spec = exp027.BRAVYI_BB["[[72,12,6]]"]
    raw = {
        "d_X": 6,
        "d_X_lower_bound": 6,
        "d_X_witness": [3, 6, 12, 15, 18, 24],
        "d_Z": None,
        "d_Z_lower_bound": 4,
        "d_Z_witness": None,
        "d_exact": False,
        "wall_time_s": 1.0,
    }
    normalized = exp027.normalize_parent_solver_result(
        spec, raw, seed=17, search_cap=6
    )
    assert normalized["lower_bound"] == 1
    assert normalized["upper_bound"] is None
    assert normalized["exact"] is False
    assert normalized["status"] == "UNSETTLED_SOLVER_BOUND"
    assert normalized["upper_bound_certified"] is False


def test_partial_pbb_solver_evidence_is_diagnostic_only() -> None:
    row = next(row for row in exp027.load_catalogue() if row["n"] == 36)
    _, code = exp027.pbb_code(row)
    logical = code.logical_basis()[0]
    raw = {
        "value": exp027.symplectic_weight(logical),
        "lower_bound": 3,
        "exact": False,
        "status": "PARTIAL",
        "witness_vector": np.asarray(logical, dtype=np.uint8).tolist(),
        "wall_time_s": 1.0,
    }
    normalized = exp027.normalize_pbb_solver_result(
        row, raw, seed=23, search_cap=int(row["d"])
    )
    assert normalized["upper_bound"] == int(row["d"])
    assert normalized["lower_bound"] == 1
    assert normalized["exact"] is False
    assert normalized["status"] == "UNSETTLED_PARTIAL_CATALOGUE_UPPER"
    assert normalized["witness_verification_diagnostic_only"]["valid"] is True


def test_partial_route_cannot_be_canonical() -> None:
    assert exp027.canonical_route(clean=True, full_coverage=False) == "partial_runs"
    assert exp027.CHECKPOINT_PATH.parent.name == "partial_runs"
    assert exp027.PROCESSED_PATH.parent.name == "processed"


def test_json_safe_rejects_cycles_and_repr_unknowns() -> None:
    cycle: list[object] = []
    cycle.append(cycle)
    try:
        exp027.json_safe(cycle)
    except TypeError as exc:
        message = str(exc)
    else:
        raise AssertionError("self-referential list was silently serialized")
    assert "cycle detected" in message
    assert "builtins.list" in message

    class Unknown:
        def __repr__(self) -> str:
            return "<unknown-safe>"

    assert exp027.json_safe(Unknown()) == "<unknown-safe>"
