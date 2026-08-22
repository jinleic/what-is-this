"""Machine checks for the PBB structural survival criterion (EXP-038).

The criterion claims ``d_PBB <= d_Z(parent)`` from linear algebra alone.  These
tests pin the three algebraic facts it rests on and, most importantly, the
falsification gate: the criterion must be *silent* on every certified reversal.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.pbb_survival import (  # noqa: E402
    dressing_space,
    quotient_rank,
    survival_test,
    translation_orbit,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402
from qec_research.symplectic.core import symplectic_product_matrix  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E27
_SPEC.loader.exec_module(E27)

# Two reversals and two dominations, spanning n = 108 and n = 144.
REVERSAL_INDEXES = [332, 366]
SAMPLE_INDEXES = [332, 366, 74, 95]


@pytest.fixture(scope="module")
def catalogue():
    return E27.load_catalogue()


def _parts(row):
    _, HX, HZ = E27.parent_matrices(row)
    _, code = E27.pbb_code(row)
    n = HX.shape[1]
    return HX, HZ, code, n, code.H[: HX.shape[0], n:]


@pytest.mark.parametrize("index", SAMPLE_INDEXES)
def test_pure_z_centralizer_is_unchanged(catalogue, index):
    """Fact 1: the PBB's pure-Z centralizer is exactly ker[A B]."""

    HX, HZ, code, n, _ = _parts(catalogue[index])
    rng = np.random.default_rng(20260817 + index)
    from qec_research.gf2.linalg import nullspace_np

    kernel = nullspace_np(HX.copy())
    assert kernel.size, "ker[A B] must be nontrivial"
    for _ in range(8):
        coeffs = rng.integers(0, 2, size=kernel.shape[0], dtype=np.uint8)
        z = (coeffs @ kernel) % 2
        vector = np.concatenate([np.zeros(n, dtype=np.uint8), z])
        # A pure-Z vector built from ker[A B] must commute with every PBB row.
        assert not (symplectic_product_matrix(code.H, vector.reshape(1, -1)) % 2).any()


@pytest.mark.parametrize("index", SAMPLE_INDEXES)
def test_dimension_identity(catalogue, index):
    """Fact 2: ``k_parent - k_pbb == dim Delta`` exactly."""

    HX, HZ, code, n, CD = _parts(catalogue[index])
    k_parent = n - rank_np(HX) - rank_np(HZ)
    delta = dressing_space(HX, CD)
    dim_delta = quotient_rank(delta, HZ)
    assert k_parent - code.k == dim_delta


@pytest.mark.parametrize("index", SAMPLE_INDEXES)
def test_translation_orbit_preserves_weight_and_logicality(catalogue, index):
    """Translations are automorphisms: same weight, still nontrivial logicals."""

    row = catalogue[index]
    HX, HZ, _, n, _ = _parts(row)
    from qec_research.gf2.linalg import nullspace_np

    kernel = nullspace_np(HX.copy())
    rng = np.random.default_rng(7 + index)
    logical = None
    for _ in range(64):
        coeffs = rng.integers(0, 2, size=kernel.shape[0], dtype=np.uint8)
        candidate = (coeffs @ kernel) % 2
        if quotient_rank(candidate.reshape(1, -1), HZ) > 0:
            logical = candidate
            break
    assert logical is not None

    orbit = translation_orbit(logical, int(row["ell"]), int(row["m"]))
    assert orbit.shape[0] == int(row["ell"]) * int(row["m"])
    weights = orbit.sum(axis=1)
    assert (weights == int(logical.sum())).all(), "translation changed weight"
    for member in orbit:
        # still in the centralizer ...
        assert not ((HX.astype(np.int64) @ member.astype(np.int64)) % 2).any()
        # ... and still a nontrivial logical of the parent
        assert quotient_rank(member.reshape(1, -1), HZ) > 0


@pytest.mark.parametrize("index", REVERSAL_INDEXES)
def test_criterion_is_silent_on_certified_reversals(catalogue, index):
    """The falsification gate.

    These rows are independently certified reversals: ``d_PBB > d_parent``.  A
    criterion that produced a surviving minimum-weight parent logical here would
    contradict a machine-checked UNSAT proof, so it must find none.
    """

    row = catalogue[index]
    HX, HZ, code, n, CD = _parts(row)
    state = ROOT / "results" / "partial_runs" / "exp036" / f"row_{index:04d}.json"
    if state.exists():
        stored = json.loads(state.read_text(encoding="utf-8"))
        assert stored["verdict"] == "CERTIFIED_REVERSAL", "fixture drifted"

    from qec_research.distance.sat_decide import (
        css_side_instance,
        decide_weight_bounded,
    )

    block = int(row["ell"]) * int(row["m"])
    instance = css_side_instance(HX, HZ, "z", block_length=block)
    weight, witness = None, None
    for cap in range(2, 13):
        record = decide_weight_bounded(instance, cap, conflict_budget=0)
        if record["status"] == "SAT":
            weight = int(record["weight"])
            witness = np.asarray(record["vector"], dtype=np.uint8)
            break
    assert witness is not None, "no parent Z-logical found below weight 13"

    verdict = survival_test(
        HX, HZ, CD, witness, int(row["ell"]), int(row["m"]),
        k_parent=n - rank_np(HX) - rank_np(HZ), k_pbb=code.k,
    )
    assert verdict.dimension_identity_holds
    assert not verdict.dominated, (
        f"criterion claimed domination on certified reversal {index} "
        f"(weight {verdict.survivor_weight}) - refuted"
    )
    # The sharp structural fact: reversals sit exactly at margin zero, i.e. the
    # dressing space is precisely large enough to absorb the whole orbit span.
    assert verdict.counting_margin == 0


def test_survivor_is_a_genuine_pbb_logical(catalogue):
    """When the criterion fires, its survivor really is a PBB logical."""

    index = 95
    row = catalogue[index]
    HX, HZ, code, n, CD = _parts(row)
    from qec_research.distance.sat_decide import (
        css_side_instance,
        decide_weight_bounded,
    )

    block = int(row["ell"]) * int(row["m"])
    instance = css_side_instance(HX, HZ, "z", block_length=block)
    witness = None
    for cap in range(2, 15):
        record = decide_weight_bounded(instance, cap, conflict_budget=0)
        if record["status"] == "SAT":
            witness = np.asarray(record["vector"], dtype=np.uint8)
            break
    assert witness is not None

    verdict = survival_test(
        HX, HZ, CD, witness, int(row["ell"]), int(row["m"]),
        k_parent=n - rank_np(HX) - rank_np(HZ), k_pbb=code.k,
    )
    assert verdict.dominated, "expected the criterion to fire on this row"

    survivor = np.concatenate(
        [np.zeros(n, dtype=np.uint8), verdict.survivor.astype(np.uint8)]
    )
    # Commutes with every PBB stabilizer generator ...
    assert not (symplectic_product_matrix(code.H, survivor.reshape(1, -1)) % 2).any()
    # ... and is not itself in the PBB stabilizer group.
    assert rank_np(np.vstack([code.H, survivor.reshape(1, -1)])) > rank_np(code.H)
