"""Regression checks for the N=105 trellis falsifier and pole audit."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TRELLIS = ROOT / "results" / "processed" / "exp061_affine_trellis_profile.json"
POLES = ROOT / "results" / "processed" / "exp062_n105_pole_stabilizers.json"
PROJECTION = (
    ROOT / "results" / "processed" / "exp065_n105_k8_projection.json"
)


def test_affine_trellis_preflight_fails_closed_on_all_hard_types() -> None:
    payload = json.loads(TRELLIS.read_text(encoding="utf-8"))
    assert payload["schema"] == "exp061-affine-trellis-profile-v1"
    assert payload["state_limit"] == 5_000_000
    assert payload["viable_targets"] == 0
    assert payload["falsified_targets"] == 4
    widths = {target["k"]: target["best"]["max_width"] for target in payload["targets"]}
    assert widths == {14: 53, 10: 80, 8: 99, 24: 48}
    for target in payload["targets"]:
        assert target["viable"] is False
        assert all(
            sector["bitset_crosschecked"] is True
            for sector in target["best"]["sectors"]
        )
        assert all(
            sector["peak_unpruned_states"] > payload["state_limit"]
            for sector in target["best"]["sectors"]
        )


def test_pole_stabilizer_prediction_is_recorded_not_retrofitted() -> None:
    payload = json.loads(POLES.read_text(encoding="utf-8"))
    assert payload["schema"] == "exp062-pole-stabilizer-audit-v1"
    assert payload["all_covers_complete"] is True
    assert len(payload["rows"]) == 20
    assert payload["histogram"] == {
        "survivors:H=5": 13,
        "undecided:H=1": 1,
        "undecided:H=5": 3,
        "undecided:H=7": 2,
        "no_reference:H=1": 1,
    }
    assert payload["prediction"] == {
        "easy_survivors_large_H": True,
        "hard_residuals_H1": False,
    }


def test_h_orbit_projection_is_exact_but_too_weak_for_k8() -> None:
    payload = json.loads(PROJECTION.read_text(encoding="utf-8"))
    assert payload["schema"] == "exp065-h-orbit-projection-v1"
    assert payload["cover_complete"] is True
    assert payload["stabilizer_size"] == 7
    assert payload["projected_length"] == 30
    assert payload["projected_stabilizer_rank"] == 11
    assert payload["projected_codimension"] == 19
    assert payload["zero_projected_classes"] == 0
    assert payload["projected_coset_leader_histogram"] == {
        "4": 45,
        "6": 195,
        "8": 15,
    }
    assert payload["rigorous_physical_distance_lower_bound"] == 4
