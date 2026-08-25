"""Regression checks for EXP-066's fail-closed n=234 frontier screen."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp066", ROOT / "experiments" / "exp066_n234_frontier.py"
)
E66 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E66
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E66)


def test_full_z39_z3_automorphism_audit_is_complete() -> None:
    audit = E66.automorphism_audit()
    assert audit["decomposition"] == "Z_39 x Z_3 ~= Z_13 x F_3^2"
    assert audit["expected_automorphisms"] == 12 * 48
    assert audit["automorphisms"] == 576
    assert audit["unique_coordinate_permutations"] == 576
    assert audit["all_bijective"] is True
    assert audit["all_homomorphisms"] is True


def test_full_automorphism_transport_preserves_both_css_rowspaces() -> None:
    source = {
        "ell": 39,
        "m": 3,
        "n": 234,
        "A": [[0, 0], [1, 0], [5, 0]],
        "B": [[0, 0], [1, 1], [23, 2]],
    }
    target = {
        "ell": 39,
        "m": 3,
        "n": 234,
        "A": [[0, 0], [1, 0], [17, 0]],
        "B": [[0, 0], [1, 1], [29, 2]],
    }
    support = [
        18, 31, 33, 34, 36, 39, 53, 65, 68,
        77, 136, 138, 145, 151, 156, 191, 201, 205,
    ]
    assert E66.verify_witness(source, support)["weight"] == 18
    transported = E66.transport_witness(source, target, support)
    assert transported["verification"]["weight"] == 18
    assert transported["transport"]["HX_rowspace_preserved"] is True
    assert transported["transport"]["HZ_rowspace_preserved"] is True


def test_persisted_n234_screen_closes_after_exp067_promotion() -> None:
    summary = json.loads(E66.SUMMARY.read_text(encoding="utf-8"))
    E66.validate_summary(summary)

    assert summary["scope"] == {
        "lattices": [[13, 9], [39, 3]],
        "n": 234,
        "k_range": [8, 24],
        "classes": 842,
        "represented_pairs": 27748,
    }
    assert summary["initial"] == {
        "dominated": 660,
        "solver_required": 182,
    }
    assert summary["automorphism_bundles"]["residual_classes"] == 182
    assert summary["automorphism_bundles"]["bundles"] == 30
    assert summary["routes"] == {
        "dominated_by_automorphism_transport": 158,
        "dominated_by_witness": 684,
    }
    assert summary["verdict"]["all_referenced_dominated"] is True
    assert summary["verdict"]["dominated"] == 842
    assert summary["verdict"]["survivors"] == 0
    assert summary["verdict"]["undecided"] == 0
    assert summary["verdict"]["exact_promotions"] == []

    promoted = json.loads(json.dumps(summary))
    promoted["verdict"]["exact_promotions"] = ["[[234,8,18]]"]
    with pytest.raises(RuntimeError, match="no replayable exact promotion"):
        E66.validate_summary(promoted)

    falsely_incomplete = json.loads(json.dumps(summary))
    falsely_incomplete["verdict"]["complete"] = False
    with pytest.raises(RuntimeError, match="verdict summary"):
        E66.validate_summary(falsely_incomplete)



def test_global_screen_records_the_n234_exact_closure() -> None:
    screen = json.loads(E66.E55.SCREEN_OUT.read_text(encoding="utf-8"))

    assert screen["scope"]["n_max"] == 234
    assert screen["scope"]["lattices_completed"] == 22
    assert screen["scope"]["missing"] == []
    assert screen["verdict"]["complete"] is True
    assert screen["verdict"]["candidates_after_symmetry"] == 4862
    assert screen["verdict"]["orbits_represented"] == 150581
    assert screen["verdict"]["with_reference"] == 4658
    assert screen["verdict"]["dominated"] == 4658
    assert screen["verdict"]["no_reference"] == 204
    assert screen["verdict"]["survivors"] == 0
    assert screen["verdict"]["undecided"] == 0
    assert screen["verdict"]["all_referenced_dominated"] is True

    assert E66._validate_global_screen(screen) == screen["verdict"]
    stale = json.loads(json.dumps(screen))
    stale["lattices"][0]["verdicts"]["no_reference"] = 0
    with pytest.raises(RuntimeError, match="shard aggregate"):
        E66._validate_global_screen(stale)

    forged = json.loads(json.dumps(screen))
    forged["lattices"][-1]["records"][-1]["verdict"] = "dominated_by_fake"
    with pytest.raises(RuntimeError, match="validated EXP-066 shard"):
        E66._validate_global_screen(forged)
