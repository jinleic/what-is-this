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
DISCOVERED = ROOT / "results" / "certificates" / "exp055_discovered_references.json"
SCREEN = ROOT / "results" / "processed" / "exp055_odd_lattice_screen.json"
SURVIVORS = ROOT / "results" / "certificates" / "exp055_odd_lattice_survivors.json"


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
    assert d.get("schema") == "exp055-literature-v3"
    return d


@pytest.fixture(scope="module")
def discovered() -> dict:
    assert DISCOVERED.exists(), f"required artifact missing: {DISCOVERED}"
    d = json.loads(DISCOVERED.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-discovered-references-v1"
    return d


@pytest.fixture(scope="module")
def screen() -> dict:
    assert SCREEN.exists(), f"required artifact missing: {SCREEN}"
    d = json.loads(SCREEN.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-odd-lattice-screen-v2"
    return d


@pytest.fixture(scope="module")
def survivor_cert() -> dict:
    assert SURVIVORS.exists(), f"required artifact missing: {SURVIVORS}"
    d = json.loads(SURVIVORS.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-odd-lattice-survivors-v2"
    return d

# --------------------------------------------------------------------------- #
# (1) artifact contracts
# --------------------------------------------------------------------------- #
def test_census_is_exhaustive_and_clean(census: dict) -> None:
    v = census["verdict"]
    assert v["lattices_swept"] == 65
    assert census["scope"]["complete"] is True
    assert census["scope"]["missing"] == []
    assert census["verdict"]["complete"] is True
    assert all(r["schema"] == "exp055-census-lattice-v2"
               and r["protocol"] == census["protocol"] for r in census["lattices"])
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
    assert lit["pole_isomorphism_all"] is True
    assert all(r["pole_isomorphism"] for r in lit["records"])
    assert lit["reported_ceiling_sanity_holds"] is True
    assert lit["reported_ceiling_violations"] == []
    assert lit["reported_distance_instances"] == 25
    assert lit["exact_ceiling_never_violated"] is True
    assert lit["exact_ceiling_violations"] == []
    assert lit["exact_distance_instances"] == 5
    assert lit["no_certificate"] == []
    # Source provenance is explicit: every W-M distance is an estimate and is
    # forbidden from the exact-reference set.
    wm = [r for r in lit["records"] if r["source"].startswith("2408.10001")]
    assert wm and all(r["source_distance_estimate"] for r in wm)
    assert sum(bool(r["distance_exact_certified_here"]) for r in wm) == 1


def test_discovered_reference_certificates(discovered: dict) -> None:
    assert discovered["all_exact"] is True
    assert discovered["all_indecomposable"] is True
    params = {(r["n"], r["k"], r["d"]) for r in discovered["records"]}
    assert params == {(30, 8, 4), (54, 8, 6), (126, 12, 10)}
    for r in discovered["records"]:
        assert r["d_X_exact"] and r["d_Z_exact"]
        assert r["d_X"] == r["d_Z"] == r["d"]
        assert r["direct_sum_component_sizes"] == [r["n"]]


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
    assert cel["ceiling"] is not None and cel["dim_pole"] == I.shape[0]
    assert cel["ceiling"] >= d_pub, (cel["ceiling"], d_pub)     # THE falsification test

    # Raw I is the LEFT annihilator; the physical right-kernel pole is J=bar(I).
    J = I[:, bar]
    vecs = E55.span_vectors(J)
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


def test_bar_convention_on_noninvariant_code() -> None:
    """Raw I is not physical; J=bar(I) is.  (7,7) catches the hidden bar bug."""
    ell, m = 7, 7
    A, B = [(1, 0), (0, 3), (0, 4)], [(0, 1), (3, 0), (4, 0)]
    HX, HZ = E55.E53.bb_from_terms(ell, m, A, B)
    I = E55.intersect(E55.ann_basis(A, ell, m), E55.ann_basis(B, ell, m))
    bar = E55.bar_permutation(ell, m)
    assert any((HX @ np.r_[u, np.zeros(ell * m, np.uint8)] % 2).any() for u in I)
    J = I[:, bar]
    assert all(not (HX @ np.r_[u, np.zeros(ell * m, np.uint8)] % 2).any() for u in J)
    cel = E55.certified_ceiling(I, ell, m, bar)
    assert cel["ceiling"] is not None and cel["ceiling"] >= 12
    mc = E55.member_ceiling(A, B, ell, m)
    assert mc["nontrivial_found"] is True and mc["ceiling"] == cel["ceiling"]


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


@pytest.mark.parametrize(
    "ell,m,A,B,d_pub",
    [
        (15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)], 10),
        (9, 9, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
        (7, 7, [(1, 0), (0, 3), (0, 4)], [(0, 1), (3, 0), (4, 0)], 12),
    ],
)
def test_reduced_witness_is_a_real_logical(ell, m, A, B, d_pub) -> None:
    """The fast screen returns a physical witness, not a decoder estimate."""
    HX, HZ = E55.E53.bb_from_terms(ell, m, A, B)
    I = E55.intersect(E55.ann_basis(A, ell, m), E55.ann_basis(B, ell, m))
    J = I[:, E55.bar_permutation(ell, m)]
    rec = E55.reduced_witness_bound(J, HZ, ell, m)
    assert rec["bound"] is not None and rec["bound"] >= d_pub
    support = rec["witness_support"]
    assert len(support) == rec["bound"]
    word = np.zeros(2 * ell * m, np.uint8)
    word[support] = 1
    assert not (HX @ word % 2).any()                # physical representative in ker H_X
    assert rank_np(np.vstack([HZ, word[None, :]])) == rank_np(HZ) + 1


def test_exp055_survivor_is_exact_and_indecomposable() -> None:
    """Independent recomputation of the first Pareto survivor: [[30,8,4]]."""
    from qec_research.distance.exact import exact_distance_css
    from qec_research.equivalence.css import stabilizer_direct_sum_certificate

    ell, m = 5, 3
    A = [(0, 0), (1, 0), (3, 1)]
    B = [(0, 0), (1, 1), (4, 1)]
    HX, HZ = E55.E53.bb_from_terms(ell, m, A, B)
    n = 2 * ell * m
    assert n - rank_np(HX) - rank_np(HZ) == 8
    dist = exact_distance_css(HX, HZ, time_limit_s=60, workers=4)
    assert dist["d_exact"] and dist["d_X"] == dist["d_Z"] == dist["d"] == 4
    zeros = np.zeros_like(HX)
    H = np.vstack([np.hstack([HX, zeros]), np.hstack([zeros, HZ])])
    cert = stabilizer_direct_sum_certificate(H)
    assert cert["component_sizes"] == [30]
    assert cert["is_direct_sum"] is False


def test_orbit_accounting_and_representative_witnesses() -> None:
    """5x3 has 11 valid unordered pairs modulo translations before unit/swap
    quotienting; every stored representative must own its stored pole witness."""
    cands = E55.enumerate_candidates(5, 3, 8, 24)
    assert len(cands) == 4
    assert sum(r["orbit"] for r in cands) == 11
    for r in cands:
        if not r.get("ceiling_witness_support"):
            continue
        HX, HZ = E55.E53.bb_from_terms(5, 3, r["A"], r["B"])
        word = np.zeros(30, np.uint8)
        word[r["ceiling_witness_support"]] = 1
        assert not (HX @ word % 2).any(), r
        assert rank_np(np.vstack([HZ, word[None, :]])) == rank_np(HZ) + 1, r


def test_screen_protocol_gate_rejects_inadequate_census() -> None:
    base = {"scope": {"complete": True},
            "protocol": {"max_weight": 3, "k_min": 8}}
    E55._validate_census_for_screen(base, 8, 24)
    with pytest.raises(RuntimeError, match="weight-3"):
        E55._validate_census_for_screen(
            {"scope": {"complete": True},
             "protocol": {"max_weight": 2, "k_min": 8}}, 8, 24)
    with pytest.raises(RuntimeError, match="above"):
        E55._validate_census_for_screen(
            {"scope": {"complete": True},
             "protocol": {"max_weight": 3, "k_min": 10}}, 8, 24)
    with pytest.raises(RuntimeError, match="exceeds"):
        E55._validate_census_for_screen(base, 8, 26)


def test_estimated_wang_mueller_distance_never_sets_threshold() -> None:
    # The n=126 value is admissible only because EXP-055 exact-certified it.
    threshold126, source126 = E55.domination_threshold(126, 12)
    assert threshold126 == 10
    assert "EXP-055" not in source126 and "2408.10001" in source126
    # The estimated [[162,8,14]] is still excluded; exact Gross gives cap 12.
    threshold162, source162 = E55.domination_threshold(162, 8)
    assert threshold162 == 12
    assert source162 == "[[144,12,12]]"
    threshold54, source54 = E55.domination_threshold(54, 8)
    assert threshold54 == 6
    assert source54 == "EXP-055 [[54,8,6]]"


def test_screen_hash_and_completeness_gates() -> None:
    protocol = E55._screen_protocol("census-a", "refs-a", 8, 24, 90.0)
    good = {
        "schema": "exp055-odd-lattice-screen-v2",
        "protocol": protocol,
        "scope": {"missing": []},
        "verdict": {"complete": True, "all_referenced_decided": True},
    }
    E55._validate_screen_for_certification(good, "census-a", "refs-a")

    for field, value, message in (
        ("census_sha256", "census-b", "current census"),
        ("reference_sha256", "refs-b", "reference set"),
    ):
        bad = json.loads(json.dumps(good))
        bad["protocol"][field] = value
        with pytest.raises(RuntimeError, match=message):
            E55._validate_screen_for_certification(bad, "census-a", "refs-a")

    incomplete = json.loads(json.dumps(good))
    incomplete["scope"]["missing"] = ["9x9"]
    incomplete["verdict"]["complete"] = False
    with pytest.raises(RuntimeError, match="incomplete"):
        E55._validate_screen_for_certification(
            incomplete, "census-a", "refs-a")

    # A shard from any different search scope must compare unequal and cannot be
    # silently reused by `_screen_lattice` / `assemble_screen`.
    assert protocol != E55._screen_protocol("census-a", "refs-a", 10, 24, 90.0)
    assert protocol != E55._screen_protocol("census-a", "refs-a", 8, 24, 120.0)


def test_fixed_point_screen_contract(screen: dict, survivor_cert: dict) -> None:
    scope, verdict = screen["scope"], screen["verdict"]
    assert scope["lattices_expected"] == scope["lattices_completed"] == 9
    assert scope["missing"] == [] and scope["n_max"] == 126
    assert scope["k_range"] == [8, 24]
    assert verdict["complete"] and verdict["all_referenced_decided"]
    assert verdict["survivors"] == verdict["undecided"] == 0
    assert verdict["all_referenced_dominated"] is True
    assert verdict["candidates_after_symmetry"] == 715
    assert verdict["orbits_represented"] == 10_247
    assert verdict["with_reference"] == verdict["dominated"] == 598
    assert verdict["no_reference"] == 117
    assert verdict["solver_calls"] == 10
    assert verdict["verdicts"] == {
        "dominated": 10, "dominated_by_witness": 588, "no_reference": 117}
    assert screen["protocol"]["census_sha256"] == E55._file_sha256(CENSUS)
    assert screen["protocol"]["reference_sha256"] == E55._reference_fingerprint()
    assert all(r["schema"] == "exp055-screen-v3"
               and r["protocol"] == screen["protocol"] for r in screen["lattices"])

    assert survivor_cert["screen_sha256"] == E55._file_sha256(SCREEN)
    assert survivor_cert["survivors"] == 0 and survivor_cert["records"] == []
    assert survivor_cert["all_exact"] and survivor_cert["all_beat_reference"]
    assert survivor_cert["all_k_match_screen"]


def test_every_persisted_screen_witness_is_physical(screen: dict) -> None:
    """The 588 solver-free domination decisions each carry a real logical word."""
    checked = 0
    for lattice in screen["lattices"]:
        for r in lattice["records"]:
            if r["verdict"] != "dominated_by_witness":
                continue
            HX, HZ = E55.E53.bb_from_terms(r["ell"], r["m"], r["A"], r["B"])
            word = np.zeros(r["n"], np.uint8)
            word[r["witness_support"]] = 1
            assert int(word.sum()) == r["witness_bound"] <= r["threshold"], r
            assert not (HX @ word % 2).any(), r
            assert rank_np(np.vstack([HZ, word[None, :]])) == rank_np(HZ) + 1, r
            checked += 1
    assert checked == 588
