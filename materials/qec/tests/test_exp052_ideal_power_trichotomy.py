"""Regression checks for EXP-052 ideal-power trichotomy.

Layers: (1) artifact schema contract; (2) independent recomputation of the six
enumeration-infeasible (k_P>20) parents and one 3-step cascade parent through
the production code path --- these are the only genuinely new, non-enumerated
conclusions; (3) synthetic immune/demote chains hand-derived from the quotient
action; (4) the case-mapping unit guard for the never-observed mixed branch.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "exp052_ideal_power_trichotomy", ROOT / "experiments" / "exp052_ideal_power_trichotomy.py"
)
assert _SPEC and _SPEC.loader
E52 = importlib.util.module_from_spec(_SPEC)
import sys as _sys

_sys.modules[_SPEC.name] = E52
_SPEC.loader.exec_module(E52)

ARTIFACT = ROOT / "results" / "processed" / "exp052_ideal_power_trichotomy.json"

HEAVY = {"phase2_90", "30_6_0007", "30_6_0008", "30_6_0012", "30_6_0013", "30_6_0039"}


@pytest.fixture(scope="module")
def artifact() -> dict:
    assert ARTIFACT.exists(), f"required artifact missing: {ARTIFACT}"
    d = json.loads(ARTIFACT.read_text())
    assert d.get("schema") == "exp052-ideal-power-trichotomy-v1-assembled", d.get("schema")
    return d


@pytest.fixture(scope="module")
def parents_and_rows():
    rows = E52.E27.load_catalogue()
    return E52.E39.distinct_parents(rows), rows


def test_summary_counts(artifact: dict) -> None:
    assert artifact["parents_total"] == 202
    assert artifact["cases"] == {"demote_full": 192, "immune": 10}
    assert artifact["mixed_labels"] == []


def test_cross_validations_are_clean(artifact: dict) -> None:
    assert artifact["immune_mismatch_vs_exp047"] == []
    assert artifact["fraction_mismatches_vs_exp050"] == []
    assert artifact["crosschecked_fractions"] == 196


def test_immune_set_is_the_exp047_ten(artifact: dict) -> None:
    assert artifact["immune_labels"] == sorted(
        ["9_6_0175", "15_6_0256", "30_6_0289", "phase2_109",
         "phase2_75", "phase2_76", "phase2_77", "phase2_83",
         "phase2_84", "phase2_87"]
    )


def test_chains_are_monotone_and_short(artifact: dict) -> None:
    for rec in artifact["records"].values():
        chain = rec["chain"]
        assert chain[0] == rec["k_parent"]
        assert chain[-1] == rec["dim_S"]
        if rec["case"] != "immune":
            assert chain[-1] == 0
            zpos = chain.index(0)
            assert all(chain[i] > chain[i + 1] for i in range(zpos))
    assert artifact["max_chain_steps"] == 3


def _fp_of(parents, label):
    return next(f for f, e in parents.items() if e["members"][0]["label"] == label)


def test_heavy_parents_recomputed_independently(artifact: dict, parents_and_rows) -> None:
    parents, rows = parents_and_rows
    for label in sorted(HEAVY):
        fp = _fp_of(parents, label)
        rec = E52.parent_chain(fp, parents[fp], rows[parents[fp]["members"][0]["catalogue_index"]])
        stored = artifact["records"][fp]
        assert rec["chain"] == stored["chain"], label
        assert rec["case"] == stored["case"] == "demote_full", label
        assert rec["dim_S"] == 0, label


def test_three_step_cascade_recomputed(artifact: dict, parents_and_rows) -> None:
    parents, rows = parents_and_rows
    fp = _fp_of(parents, "phase2_49")
    rec = E52.parent_chain(fp, parents[fp], rows[parents[fp]["members"][0]["catalogue_index"]])
    stored = artifact["records"][fp]
    assert rec["chain"] == stored["chain"] == [12, 4, 0, 0]
    assert rec["case"] == "demote_full"
    assert rec["fraction_pred"] == 1.0


def test_flagship_chain_recomputed(artifact: dict, parents_and_rows) -> None:
    parents, rows = parents_and_rows
    fp = _fp_of(parents, "12_6_0193")
    rec = E52.parent_chain(fp, parents[fp], rows[parents[fp]["members"][0]["catalogue_index"]])
    stored = artifact["records"][fp]
    assert rec["chain"] == stored["chain"] == [12, 0, 0]
    assert rec["case"] == "demote_full"
    assert rec["fraction_pred"] == 1.0


def _chain_for(ell: int, m: int, a_terms, b_terms, k_p: int):
    entry = {"ell": ell, "m": m, "k_parent": k_p,
             "members": [{"label": "synthetic", "catalogue_index": 0}]}
    row = {"ell": ell, "m": m, "A_terms": a_terms, "B_terms": b_terms}
    return E52.parent_chain("synthetic", entry, row)


def test_synthetic_immune_parent() -> None:
    # R = GF(2)[x]/(x^3-1), A = B = 1+x+x^2: the pair vanishes on the GF(4)
    # factor, so M is carried entirely by I-full factors: I M = M, k_P = 4.
    rec = _chain_for(3, 1, [[0, 0], [1, 0], [2, 0]], [[0, 0], [1, 0], [2, 0]], 4)
    assert rec["chain"] == [4, 4]
    assert rec["case"] == "immune"
    assert rec["fraction_pred"] == 0.0


def test_synthetic_demote_parent_one_descent() -> None:
    # R = GF(2)[x]/(x^6-1), A = B = x^2+x+1, k_P = 4: M is carried by the
    # (x^2+x+1)^2 block where the pair is nilpotent; every image t*y of the
    # ideal (x+1)^2 (x^2+x+1) already lies in S_Z, so I M = 0 in one step.
    rec = _chain_for(6, 1, [[2, 0], [1, 0], [0, 0]], [[2, 0], [1, 0], [0, 0]], 4)
    assert rec["chain"] == [4, 0, 0]
    assert rec["case"] == "demote_full"


def test_synthetic_demote_parent_nilpotent_block() -> None:
    # R = GF(2)[x,y]/((x+1)^2,(y+1)^2), A = 1+x, B = 1+y, k_P = 2:
    # I = (1+x)(1+y) R annihilates M immediately.
    rec = _chain_for(2, 2, [[0, 0], [1, 0]], [[0, 0], [0, 1]], 2)
    assert rec["chain"] == [2, 0, 0]
    assert rec["case"] == "demote_full"


def test_case_mapping_covers_mixed_branch() -> None:
    assert E52.classify_case(4, 4) == ("immune", 0.0)
    assert E52.classify_case(0, 4) == ("demote_full", 1.0)
    case, frac = E52.classify_case(2, 4)
    assert case == "mixed"
    assert abs(frac - (1.0 - 3 / 15)) < 1e-12
