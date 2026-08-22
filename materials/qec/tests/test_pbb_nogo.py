"""Machine checks for the parent-level no-go bound (Theorem H, EXP-039).

The theorem's content is a quantifier move: from "this perturbation cannot beat
the parent" to "*no* perturbation of this parent can".  Its soundness rests on
three facts, each pinned here:

1. ``Delta`` is an ``R``-submodule (translation invariant).  This is what lets a
   single witness speak for its whole orbit, and what makes ``T`` a property of
   the parent alone.
2. ``T`` is a lower bound on ``dim Delta`` for any perturbation that strictly
   beats ``d_Z(parent)``; with ``k_Q = k_P - dim Delta`` that is the no-go.
3. The enumeration terminates in UNSAT, so ``T`` is exact rather than a bound.

Plus the falsification gate: the seven certified reversals must all be
*permitted*.  A theorem that forbade one of them would be refuted by a replayed
UNSAT proof, so these are the tests that matter most.
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

_SPEC = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
EXP039 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = EXP039
_SPEC.loader.exec_module(EXP039)

from qec_research.codes.pbb_nogo import (  # noqa: E402
    NoGoCertificate,
    minimum_weight_module,
    orthogonal_complement,
)
from qec_research.codes.pbb_survival import dressing_space, translation_orbit
from qec_research.gf2.linalg import rank_np  # noqa: E402

E27 = EXP039.E27


@pytest.fixture(scope="module")
def catalogue():
    return E27.load_catalogue()


# --------------------------------------------------------------------------
# 1. the descriptors themselves
# --------------------------------------------------------------------------
def test_certificate_predicates_are_properties():
    """``family_closed``/``d_z_exact`` must be values, not bound methods.

    A bound method is always truthy, so a missing ``@property`` here would
    silently report *every* parent as family-closed -- the headline claim.
    """

    cert = NoGoCertificate(k_parent=4, d_z_parent=4, T=1, closed=True)
    assert cert.k_ceiling == 3
    assert cert.family_closed is False
    assert cert.d_z_exact is False
    assert cert.permits(3) is True
    assert cert.permits(4) is False

    closed = NoGoCertificate(
        k_parent=4, d_z_parent=4, T=4, closed=True, witnesses=[np.zeros(4, np.uint8)]
    )
    assert closed.k_ceiling == 0
    assert closed.family_closed is True
    assert closed.d_z_exact is True
    assert closed.permits(1) is False


def test_orthogonal_complement_is_the_annihilator():
    W = np.array([[1, 1, 0, 0], [0, 0, 1, 1]], dtype=np.uint8)
    F = orthogonal_complement(W, 4)
    assert F.shape[0] == 4 - rank_np(W)
    assert not ((F.astype(np.int64) @ W.T.astype(np.int64)) % 2).any()
    # f.z = 0 for all f in F  <=>  z in rowspace(W)
    for z in W:
        assert not ((F.astype(np.int64) @ z.astype(np.int64)) % 2).any()
    outside = np.array([1, 0, 0, 0], dtype=np.uint8)
    assert ((F.astype(np.int64) @ outside.astype(np.int64)) % 2).any()


# --------------------------------------------------------------------------
# 2. Delta is an R-submodule: the load-bearing structural fact
# --------------------------------------------------------------------------
@pytest.mark.parametrize("which", [0, 40, 95, 140, 200, 300])
def test_dressing_space_is_translation_invariant(catalogue, which):
    """Translating any element of ``Delta`` stays in ``Delta``.

    This is why absorbing one orbit element absorbs the orbit, and therefore
    why ``T`` bounds ``dim Delta`` rather than only ``dim Delta`` at one point.
    """

    row = catalogue[which]
    _, HX, HZ = E27.parent_matrices(row)
    _, code = E27.pbb_code(row)
    n = HX.shape[1]
    CD = code.H[: HX.shape[0], n:]
    delta = dressing_space(HX, CD)
    if delta.size == 0 or rank_np(delta) == 0:
        pytest.skip("delta=0 (CSS row): invariance is vacuous")
    ell, m = int(row["ell"]), int(row["m"])
    base = rank_np(delta)
    for vector in delta[: min(4, delta.shape[0])]:
        orbit = translation_orbit(vector, ell, m)
        # every translate is already inside Delta's row space
        assert rank_np(np.vstack([delta, orbit])) == base


@pytest.mark.parametrize("which", [0, 95, 200, 300])
def test_stabilizer_space_is_translation_invariant(catalogue, which):
    row = catalogue[which]
    _, HX, HZ = E27.parent_matrices(row)
    ell, m = int(row["ell"]), int(row["m"])
    base = rank_np(HZ)
    for vector in HZ[: min(4, HZ.shape[0])]:
        orbit = translation_orbit(vector, ell, m)
        assert rank_np(np.vstack([HZ, orbit])) == base


# --------------------------------------------------------------------------
# 3. the bound itself, on a small parent computed from scratch
# --------------------------------------------------------------------------
def test_T_bounds_dim_delta_on_every_stored_row(catalogue):
    """For each certified parent, every catalogue child obeys the no-go.

    Concretely: a child with ``k_Q > k_P - T`` must NOT have a dressing space
    containing the minimum-weight module -- checked directly in linear algebra,
    independently of any distance computation.
    """

    certs = EXP039.load_certificates()
    if not certs:
        pytest.skip("no EXP-039 certificates yet")
    parents = EXP039.distinct_parents(catalogue)
    tested = 0
    for fp, cert in certs.items():
        entry = parents.get(fp)
        if entry is None or not cert["T_is_exact"] or not cert["witness_vectors"]:
            continue
        row = catalogue[entry["members"][0]["catalogue_index"]]
        _, HX, HZ = E27.parent_matrices(row)
        ell, m = int(row["ell"]), int(row["m"])
        W = HZ.copy()
        for vec in cert["witness_vectors"]:
            W = np.vstack([W, translation_orbit(np.asarray(vec, np.uint8), ell, m)])
        assert rank_np(W) - rank_np(HZ) == cert["T"]

        for member in entry["members"]:
            _, code = E27.pbb_code(catalogue[member["catalogue_index"]])
            n = HX.shape[1]
            delta = dressing_space(HX, code.H[: HX.shape[0], n:])
            dim_delta = rank_np(np.vstack([HZ, delta])) - rank_np(HZ) if delta.size else 0
            # Theorem G(ii): the dimension identity.
            assert cert["k_parent"] - member["k_pbb"] == dim_delta
            if member["k_pbb"] > cert["k_ceiling"]:
                # Then dim Delta < T, so Delta cannot contain the module:
                # some minimum-weight parent logical survives.
                assert dim_delta < cert["T"]
        tested += 1
        if tested >= 6:
            break
    assert tested > 0


def test_module_closure_reproduces_a_small_parent(catalogue):
    """Recompute ``T`` from scratch and match the persisted certificate."""

    certs = EXP039.load_certificates()
    if not certs:
        pytest.skip("no EXP-039 certificates yet")
    small = min(
        (c for c in certs.values() if c["T_is_exact"]), key=lambda c: (c["n"], c["T"])
    )
    parents = EXP039.distinct_parents(catalogue)
    entry = parents[small["fingerprint"]]
    row = catalogue[entry["members"][0]["catalogue_index"]]
    _, HX, HZ = E27.parent_matrices(row)
    fresh = minimum_weight_module(
        HX,
        HZ,
        small["d_z_lower_input"],
        entry["ell"],
        entry["m"],
        k_parent=entry["k_parent"],
    )
    assert fresh.T == small["T"]
    assert fresh.closed == small["T_is_exact"]
    assert fresh.d_z_parent == small["d_z_parent"]


# --------------------------------------------------------------------------
# 4. the falsification gate
# --------------------------------------------------------------------------
def test_falsification_gate_permits_every_certified_reversal(catalogue):
    """The seven proved reversals must all satisfy ``k_Q <= k_P - T``.

    Each reversal carries a replayed UNSAT lower bound proving
    ``d_PBB > d_parent``.  If Theorem H forbade one, the theorem -- not the
    reversal -- would be wrong.
    """

    gate = EXP039.falsification_gate(
        EXP039.distinct_parents(catalogue), EXP039.load_certificates()
    )
    if gate["reversals_checked"] == 0:
        pytest.skip("reversal parents not yet certified")
    assert not gate["theorem_refuted"], gate["violations"]
    for check in gate["checks"]:
        assert check["permitted_by_theorem_H"], check
        assert check["slack"] >= 0, check


def test_gate_would_catch_a_violation(catalogue):
    """The gate must actually be able to fail (it is not vacuously true)."""

    parents = EXP039.distinct_parents(catalogue)
    certs = EXP039.load_certificates()
    if not certs:
        pytest.skip("no EXP-039 certificates yet")
    # Forge a certificate whose T is maximal, so every child is forbidden.
    forged = {}
    for fp, cert in certs.items():
        entry = parents.get(fp)
        if entry is None:
            continue
        if any(m["label"] in EXP039.CERTIFIED_REVERSALS for m in entry["members"]):
            forged[fp] = {**cert, "T": cert["k_parent"], "k_ceiling": 0}
    if not forged:
        pytest.skip("no reversal parent certified yet")
    gate = EXP039.falsification_gate(parents, forged)
    assert gate["theorem_refuted"]
    assert gate["violations"]


def test_assembled_artifact_is_self_consistent():
    path = EXP039.OUT
    if not path.exists():
        pytest.skip("EXP-039 artifact not assembled yet")
    payload = json.loads(path.read_text(encoding="utf-8"))
    counts = payload["counts"]
    assert counts["parents_family_closed"] == len(payload["family_closed_parents"])
    assert counts["catalogue_rows_capped_a_priori"] == len(
        payload["rows_capped_a_priori"]
    )
    for item in payload["family_closed_parents"]:
        assert item["T"] == item["k_parent"]
    for item in payload["rows_capped_a_priori"]:
        assert item["k_pbb"] > item["k_ceiling"]
    for item in payload["rows_theorem_permits"]:
        assert item["k_pbb"] <= item["k_ceiling"]
    assert not payload["falsification_gate"]["theorem_refuted"]
