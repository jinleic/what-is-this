"""Machine checks for Theorems H and I and their lemmas (EXP-039/040).

Theorem I (forced saturation) is proved on paper in proofs/pbb_structure.md;
these tests verify its two computational legs and its scope claim:

  Lemma 2  k_P = 2 dim L  for every CSS BB parent, where L is the left
           kernel of [A B]  (checked on ALL distinct parents, certified or
           not — the proof is elementary and the check is cheap);
  Lemma 3  dim(Delta_bar) <= k_P/2 on every catalogue row, via the Theorem-G
           identity dim(Delta_bar) = k_P - k_Q;
  scope    T >= k_P/2 on every certified parent except the single known
           exception 9a7638586033* (n=144, k_P=12, T=4);
  probe    the EXP-040 small-lattice aggregate: zero upper-law violations,
           all strict increases at equality.
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

from qec_research.gf2.linalg import nullspace_np  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
assert _SPEC and _SPEC.loader
E39 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(E39)
E27 = E39.E27

SATURATION = ROOT / "results" / "processed" / "exp040_saturation_probe.json"
SCOPE_EXCEPTION_PREFIX = "9a7638586033"  # n=144, k_P=12, T=4


@pytest.fixture(scope="module")
def catalogue():
    return E27.load_catalogue()


def test_lemma2_left_kernel_dimension_all_parents(catalogue) -> None:
    # k_P = 2 dim L is elementary (proofs/pbb_structure.md, Lemma 2); this
    # verifies it on every distinct parent of the whole catalogue, not just
    # the certified ones.
    parents = E39.distinct_parents(catalogue)
    by_index = {i: row for i, row in enumerate(catalogue)}
    assert len(parents) == 202
    for entry in parents.values():
        idx = entry["members"][0]["catalogue_index"]
        _, HX, _ = E27.parent_matrices(by_index[idx])
        # Left kernel of [A B]: nullspace of its transpose.
        dim_L = nullspace_np(np.asarray(HX, dtype=np.uint8).T).shape[0]
        assert 2 * dim_L == int(entry["k_parent"]), (
            f"parent {entry['fingerprint'][:12]}: dim L={dim_L}, "
            f"k_P={entry['k_parent']} -- Lemma 2 violated"
        )


def test_lemma3_dressing_ceiling_every_row(catalogue) -> None:
    # dim(Delta_bar) = k_P - k_Q (Theorem G(ii)); Lemma 3 caps it at k_P/2.
    parents = E39.distinct_parents(catalogue)
    for entry in parents.values():
        k_p = int(entry["k_parent"])
        for member in entry["members"]:
            dim_bar = k_p - int(member["k_pbb"])
            assert dim_bar <= k_p / 2, (
                f"{member['label']}: dim(Delta_bar)={dim_bar} > k_P/2={k_p / 2} "
                "-- Lemma 3 violated"
            )


def test_theorem_I_scope_single_exception() -> None:
    # T >= k_P/2 must hold on every certified parent except the one known
    # exception.  A second exception appearing as the sweep progresses is a
    # genuine discovery (the theorem's hypothesis class would widen) and must
    # fail here rather than pass silently.
    certs = E39.load_certificates()
    exceptions = [
        fp
        for fp, c in certs.items()
        if 2 * int(c["T"]) < int(c["k_parent"])
    ]
    assert all(fp.startswith(SCOPE_EXCEPTION_PREFIX) for fp in exceptions), (
        f"new T < k_P/2 parents appeared: "
        f"{[fp[:12] for fp in exceptions if not fp.startswith(SCOPE_EXCEPTION_PREFIX)]}"
    )


@pytest.mark.skipif(not SATURATION.exists(), reason="EXP-040 probe artifact absent")
def test_saturation_probe_aggregate() -> None:
    payload = json.loads(SATURATION.read_text(encoding="utf-8"))
    small = payload["small_lattice"]
    agg = small["aggregate"]
    assert agg["counterexamples"] == 0
    assert agg["d_Q_gt_d_Z"] == 496
    assert agg["delta_positive_perturbations"] == 26898

    # Per-parent relation tallies: strictly-below, at, never above.
    below = at = above = 0
    for lattice in small["lattices"]:
        for parent in lattice["parent_records"]:
            rel = parent["perturbations"]["delta_bar_relation_to_T"]
            below += rel.get("less_than_T", 0)
            at += rel.get("equal_T", 0)
            above += rel.get("greater_than_T", 0)
    assert above == 0
    assert below + at == agg["delta_positive_perturbations"]

    # Every probed parent satisfies the Theorem-I hypothesis class: T >= k_P/2
    # (so its strict examples at equality were theorem-forced, not lucky).
    for lattice in small["lattices"]:
        for parent in lattice["parent_records"]:
            assert 2 * int(parent["T"]) >= int(parent["k_P"]), (
                f"{parent['id']}: T={parent['T']} < k_P/2={parent['k_P'] / 2} "
                "-- probe left the hypothesis class of Theorem I"
            )

    # The 7 catalogue reversals: not just dim(Delta_bar) == T but actual
    # submodule containment, both directions.
    for row in payload["certified_reversals"]:
        assert int(row["dim_Delta_bar"]) == int(row["T"]) == int(row["k_drop"])
        assert row["equality_dim_Delta_bar_eq_T"] is True
        assert row["M_bar_contained_in_Delta_bar"] is True
        assert row["Delta_bar_equals_M_bar"] is True
        assert int(row["catalogue_k_Q"]) == int(row["k_Q"])
