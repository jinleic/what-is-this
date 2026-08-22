"""Regression checks for EXP-051 channel queries and the assembled artifact."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp051_xmonotone_lightscan", ROOT / "experiments" / "exp051_xmonotone_lightscan.py"
)
assert SPEC and SPEC.loader
EXP051 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP051)

ARTIFACT = ROOT / "results" / "processed" / "exp051_xmonotone_lightscan.json"


def test_channel_query_finds_multirow_light() -> None:
    """synthetic SAT case: r3 = r1+r2 must be found at its weight."""
    HX = np.array(
        [[1, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0]], np.uint8
    )
    result = EXP051.channel_query("synthetic-sat", HX, cap=2)
    assert result["status"] in ("OPTIMAL", "FEASIBLE")
    # any feasible witness: weight 2 == r1+r2 (only multi-row light at cap 2)
    assert result["witness_weight"] == 2
    assert result["witness"] is not None
    word = np.array(result["witness"], dtype=np.uint8)
    # not an original row
    assert not any(np.array_equal(word, HX[i]) for i in range(3))


def test_channel_query_infeasible_when_no_multirow_light() -> None:
    """synthetic INFEASIBLE case: catches the XOR-parity bug class.

    Two disjoint rows: rowspace = {0, r1, r2, r1+r2}; r1+r2 has weight 4 > cap,
    so no non-original codeword with weight <= 1, and the query must be
    INFEASIBLE rather than wrongly returning a parity-flipped codeword.
    """
    HX = np.array(
        [[1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0]], np.uint8
    )
    result = EXP051.channel_query("synthetic-infeasible", HX, cap=1)
    assert result["status"] == "INFEASIBLE"
    assert result["complete"]


def test_verdict_vocab_and_headline() -> None:
    assert EXP051.SCAN_TIME_LIMIT_S == 900.0
    d = json.loads(ARTIFACT.read_text())
    records = {r["label"]: r for r in d["records"]}
    # exactly the declared target set
    assert set(records) == set(EXP051.TARGETS) | set(EXP051.OPEN_TARGETS)
    certified = [r for r in d["records"] if r["verdict"] == "CERTIFIED_NO_MULTIROW_LIGHT"]
    channels = [r for r in d["records"] if r["verdict"] == "CHANNEL_PRESENT"]
    assert len(certified) == 8, sorted(records)
    assert len(channels) == 2
    assert all(c["independent_verified"] for c in channels)
    assert all(c["channel_witness_weight"] == 8 for c in channels)
    assert d["all_xm_certified"] is True
