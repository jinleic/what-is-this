"""Regression checks for EXP-055: odd-lattice census, Theorem K, literature battery.

Layers:
  (1) artifact schema contracts for the census and the literature battery;
  (2) INDEPENDENT in-test recomputation of Theorem K's *mechanism* on published
      instances -- not just the recorded number: the witness (u,0) is rebuilt and
      verified to lie in ker H_X and outside rowspace(H_Z);
  (3) an independent GF(2) gcd implementation for the published rate law;
  (4) a guard that the two unreproducible literature rows stay excluded from
      every pass/fail count, so the battery cannot silently "improve".
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
    "exp055", ROOT / "experiments" / "exp055_odd_lattice_sweep.py")
E55 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E55
_SPEC.loader.exec_module(E55)

from qec_research.gf2.linalg import rank_np  # noqa: E402

CENSUS = ROOT / "results" / "processed" / "exp055_odd_lattice_sweep.json"
LIT = ROOT / "results" / "processed" / "exp055_literature_validation.json"


@pytest.fixture(scope="module")
def census() -> dict:
    assert CENSUS.exists(), f"required artifact missing: {CENSUS}"
    d = json.loads(CENSUS.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-odd-lattice-sweep-v1"
    return d


@pytest.fixture(scope="module")
def lit() -> dict:
    assert LIT.exists(), f"required artifact missing: {LIT}"
    d = json.loads(LIT.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-literature-v2"
    return d


# --------------------------------------------------------------------------- #
# (1) artifact contracts
# --------------------------------------------------------------------------- #
def test_census_is_exhaustive_and_clean(census: dict) -> None:
    v = census["verdict"]
    assert v["lattices_swept"] == 65
    assert v["pairs_total"] > 4e9
    # k by the ideal route vs k = n - rank H_X - rank H_Z must never disagree
    assert v["k_mismatches"] == 0
    # semisimplicity: every annihilator is idempotent on an odd lattice (J-G1)
    assert v["idempotence_violations"] == 0
    assert v["idempotence_tested"] >= 200
    # the sweep must actually cover lattices outside our own catalogue's m in {3,6}
    lats = {(r["ell"], r["m"]) for r in census["lattices"]}
    assert (9, 9) in lats and (7, 7) in lats and (13, 13) in lats
    assert all(r["ell"] % 2 and r["m"] % 2 for r in census["lattices"])


def test_literature_battery_contract(lit: dict) -> None:
    assert lit["instances"] == 27
    assert lit["reproduced_instances"] == 25
    assert lit["k_all_agree_on_reproduced"] is True
    assert lit["our_routes_always_agree"] is True
    assert lit["ceiling_never_violated"] is True
    assert lit["ceiling_violations"] == []
    assert lit["ceiling_certified_count"] == 25
    assert lit["no_certificate"] == []
    # the ceiling is sound but loose: keep the measured looseness honest
    assert lit["slack_min"] >= 0
    assert lit["slack_median"] >= 2


def test_unreproducible_rows_stay_flagged(lit: dict) -> None:
    """Two App.C rows give k=0 by both our routes; they must stay excluded."""
    bad = lit["unreproducible_rows"]
    assert len(bad) == 2
    for r in bad:
        assert (r["ell"], r["m"]) in {(5, 9), (7, 11)}
        assert r["k_ideal"] == r["k_matrices"] == 0
        assert r["k_published"] in (4, 6)
        assert r["our_routes_agree"] is True     # our two routes agree with each other


# --------------------------------------------------------------------------- #
# (2) Theorem K's mechanism, recomputed in-test
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "ell,m,A,B,k_pub,d_pub",
    [
        # Bravyi et al. arXiv:2308.07915 Table 3
        (15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)], 8, 10),
        # Eberhardt-Steffan arXiv:2407.03973v1 Table 1
        (9, 9, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 8, 12),
        (9, 9, [(0, 0), (1, 0), (0, 1)], [(3, 0), (0, 1), (0, 2)], 4, 16),
        (9, 15, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 8, 18),
    ],
)
def test_theorem_k_mechanism_on_published_codes(ell, m, A, B, k_pub, d_pub) -> None:
    """The ceiling must hold AND its witness must really be a nontrivial logical."""
    n = 2 * ell * m
    HX, HZ = E55.E53.bb_from_terms(ell, m, A, B)
    assert int(n - rank_np(HX) - rank_np(HZ)) == k_pub          # rate law, matrices
    I = E55.intersect(E55.ann_basis(A, ell, m), E55.ann_basis(B, ell, m))
    assert 2 * I.shape[0] == k_pub                              # rate law, ideal
    bar = E55.bar_permutation(ell, m)
    cel = E55.certified_ceiling(I, ell, m, bar)
    assert cel["ceiling"] is not None and cel["dim_I0"] > 0
    assert cel["ceiling"] >= d_pub, (cel["ceiling"], d_pub)     # THE falsification test

    # rebuild the minimising witness and verify the theorem's two claims directly
    I0 = I if cel["bar_invariant"] else E55.intersect(I, I[:, bar])
    vecs = E55.span_vectors(I0)
    w = vecs.sum(axis=1)
    u = vecs[int(np.argmin(w))]
    assert int(u.sum()) == cel["ceiling"]
    full = np.zeros(n, np.uint8)
    full[: ell * m] = u
    assert not (HX @ full % 2).any()                            # in ker H_X
    rz = rank_np(HZ)
    assert rank_np(np.vstack([HZ, full[None, :]])) == rz + 1    # outside S_Z


def test_bar_involution_is_an_involution() -> None:
    for (ell, m) in [(9, 9), (15, 3), (7, 7), (15, 9)]:
        bar = E55.bar_permutation(ell, m)
        assert np.array_equal(bar[bar], np.arange(ell * m))


def test_member_fallback_when_pole_ideal_vanishes() -> None:
    """(7,7) [[98,6,12]] has Z cap Z^-1 empty: Theorem K is vacuous, Cor K1 is not."""
    ell, m = 7, 7
    A, B = [(1, 0), (0, 3), (0, 4)], [(0, 1), (3, 0), (4, 0)]
    I = E55.intersect(E55.ann_basis(A, ell, m), E55.ann_basis(B, ell, m))
    bar = E55.bar_permutation(ell, m)
    cel = E55.certified_ceiling(I, ell, m, bar)
    assert cel["dim_I0"] == 0 and cel["ceiling"] is None        # theorem vacuous here
    mc = E55.member_ceiling(A, B, ell, m)
    assert mc["nontrivial_found"] is True
    assert mc["ceiling"] >= 12                                  # published exact d


# --------------------------------------------------------------------------- #
# (3) the published rate law, independent implementation
# --------------------------------------------------------------------------- #
def test_gcd_rate_law_matches_annihilator_route() -> None:
    """k = 2 deg gcd(a(pi), b(pi), pi^N-1) on coprime lattices (PK Prop. 1),
    computed with bitmask polynomials that share no code with the ring route."""
    def deg(p: int) -> int:
        return p.bit_length() - 1

    def pmod(a: int, b: int) -> int:
        db = deg(b)
        while a and deg(a) >= db:
            a ^= b << (deg(a) - db)
        return a

    def pgcd(a: int, b: int) -> int:
        while b:
            a, b = b, pmod(a, b)
        return a

    cases = [  # (ell, m, pi_A, pi_B, k) from arXiv:2408.10001v4 Table 2
        (3, 5, [0, 1, 2], [1, 3, 8], 4),
        (3, 7, [0, 2, 3], [1, 3, 11], 6),
        (5, 7, [0, 1, 5], [0, 1, 12], 6),
        (7, 9, [0, 1, 58], [3, 16, 44], 12),
    ]
    for ell, m, ea, eb, k_pub in cases:
        N = ell * m
        a = 0
        b = 0
        for e in ea:
            a ^= 1 << (e % N)
        for e in eb:
            b ^= 1 << (e % N)
        k_gcd = 2 * deg(pgcd(pgcd(a, b), (1 << N) | 1))
        assert k_gcd == k_pub, (ell, m, k_gcd, k_pub)
        # and the module's own coprime route must agree
        terms_a = [(e % ell, e % m) for e in ea]
        terms_b = [(e % ell, e % m) for e in eb]
        assert E55.cyclic_k(terms_a, terms_b, ell, m) == k_pub


def test_symmetry_reduction_preserves_k_and_ceiling() -> None:
    """Every generator of the reduction group is a qubit permutation, so k and the
    certified ceiling are invariant along an orbit."""
    ell, m = 9, 9
    A, B = [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]
    bar = E55.bar_permutation(ell, m)

    def invariants(P, Q):
        I = E55.intersect(E55.ann_basis(P, ell, m), E55.ann_basis(Q, ell, m))
        cel = E55.certified_ceiling(I, ell, m, bar)
        HX, HZ = E55.E53.bb_from_terms(ell, m, P, Q)
        return (2 * I.shape[0], cel["ceiling"],
                int(2 * ell * m - rank_np(HX) - rank_np(HZ)))

    base = invariants(A, B)
    variants = [
        ([((a + 4) % ell, b) for a, b in A], B),                  # translate A
        (A, [(a, (b + 5) % m) for a, b in B]),                    # translate B
        ([((a * 2) % ell, (b * 5) % m) for a, b in A],
         [((a * 2) % ell, (b * 5) % m) for a, b in B]),           # unit map
        ([(b, a) for a, b in A], [(b, a) for a, b in B]),         # x <-> y
        (B, A),                                                    # block swap
    ]
    for P, Q in variants:
        assert invariants(P, Q) == base, (P, Q)
