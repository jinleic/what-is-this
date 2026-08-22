"""Regression guards for EXP-016 claim semantics (FR-013/FR-016)."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "raw" / "exp016_schedule_controlled.json"


def test_point_rank_and_familywise_separation_are_distinct_claims() -> None:
    payload = json.loads(ARTIFACT.read_text())
    css, pbb = payload["codes"]
    verdict = payload["verdict"]

    point_separation = (
        min(row["ler_per_shot"] for row in pbb["runs"])
        > max(row["ler_per_shot"] for row in css["runs"])
    )
    familywise_separation = (
        pbb["min_simultaneous_lo"] > css["max_simultaneous_hi"]
    )

    assert point_separation is True
    assert familywise_separation is False
    assert verdict["every_sampled_pbb_point_estimate_exceeds_every_css"] is True
    assert verdict["every_pair_separated_at_familywise_95"] is False
    assert "every_sampled_pbb_worse_than_every_sampled_css" not in verdict


def test_claim_names_translation_invariant_sample_scope() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert "TRANSLATION-INVARIANT" in payload["verdict"]["scope"]
    assert "non-TI" in payload["verdict"]["scope"]
    for code in payload["codes"]:
        assert code["enumeration_complete"] is True
        assert code["total_valid_schedules"] > 0
        assert len(code["runs"]) == 5
