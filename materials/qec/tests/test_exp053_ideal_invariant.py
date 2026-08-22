"""Regression checks for Theorems J-G / J-H / J-I / J-J (EXP-053, EXP-054).

Layers:
 (1) artifact schema contracts for both experiments;
 (2) INDEPENDENT in-test recomputation of the genuinely new conclusions --
     the mixed witness (all three routes), the published-baseline table
     (including the odd-lattice [[90,8,10]] immunity), and the sharpened
     catalogue facts (I^2 = 0 on demoting parents, I^2 = I on immune ones);
 (3) direct unit checks of the two structural lemmas behind Theorem J-I
     (exact-vanishing coset law, scalar law) plus the weight-4 tightness case;
 (4) prose-independent guards on the census verdict.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


E53 = _load("exp053_ideal_classification", "exp053_ideal_classification.py")

ART53 = ROOT / "results" / "processed" / "exp053_ideal_classification.json"
ART54 = ROOT / "results" / "processed" / "exp054_mixed_census.json"


@pytest.fixture(scope="module")
def art53() -> dict:
    assert ART53.exists(), f"required artifact missing: {ART53}"
    d = json.loads(ART53.read_text())
    assert d["schema"] == "exp053-ideal-classification-v1"
    return d


@pytest.fixture(scope="module")
def art54() -> dict:
    assert ART54.exists(), f"required artifact missing: {ART54}"
    d = json.loads(ART54.read_text())
    assert d["schema"] == "exp054-mixed-census-v1"
    return d


# ----------------------------------------------------------- artifact contract

def test_catalogue_cases_and_agreement(art53: dict) -> None:
    v = art53["verdict"]
    assert v["parents"] == 202
    assert v["cases"] == {"demote_full": 192, "immune": 10, "mixed": 0, "degenerate": 0}
    assert v["ideal_route_matches_exp052"] is True
    assert art53["catalogue"]["mismatches"] == []


def test_sharpened_catalogue_extremes(art53: dict) -> None:
    """Every demoting parent has I^2 = 0; every immune parent has I^2 = I."""
    for r in art53["catalogue"]["records"]:
        if r["case"] == "demote_full":
            assert r["power_dims"][1] == 0, r["label"]
            assert r["nilpotency_index"] == 2, r["label"]
        else:
            assert r["case"] == "immune"
            assert r["power_dims"] == [r["dim_I"], r["dim_I"]], r["label"]
            assert 2 * r["dim_I"] == r["k_parent"], r["label"]


def test_dim_S_equals_twice_stable_ideal_power(art53: dict) -> None:
    for r in art53["catalogue"]["records"]:
        assert r["dim_S_pred"] == 2 * r["dim_I_infty"] == r["dim_S_exp052"], r["label"]
        assert r["k_pred"] == r["k_parent"] == 2 * r["dim_I"], r["label"]


# --------------------------------------------- independent recomputation layer

def test_mixed_witness_recomputed_all_three_routes() -> None:
    """Theorem J-J from scratch: ideal chain, module F-chain, brute-force census."""
    ell, m = 2, 3
    A = [(0, 1), (0, 2), (1, 1), (1, 2)]      # (1+x)(y+y^2)
    B = [(0, 2), (0, 0), (1, 2), (1, 0)]      # A*y
    HX, HZ = E53.bb_from_terms(ell, m, A, B)
    from qec_research.gf2.linalg import rank_np

    k_P = HX.shape[1] - rank_np(HX) - rank_np(HZ)
    assert (HX.shape[1], k_P) == (12, 8)

    cls = E53.classify_ideal(E53.ideal_of_parent(HX, ell, m), ell, m)
    assert cls["case"] == "mixed"
    assert cls["power_dims"] == [4, 2, 2]
    assert (cls["dim_I"], cls["dim_I_infty"], cls["dim_S_pred"]) == (4, 2, 4)

    chain = E53.module_chain(HX, HZ, ell, m, k_P)
    assert chain == [8, 4, 4]                      # route 2 agrees with route 1

    census = E53.demote_census(HX, HZ, ell, m, k_P)
    assert census == {"classes_total": 255, "classes_demote": 240,
                      "fraction": 240 / 255}       # route 3 agrees: 2^8 - 2^4
    assert E53.coset_support(A, ell, m) == {(0, 0), (1, 0)}   # two G_2-cosets


def test_published_baselines_recomputed() -> None:
    """Only the odd x odd [[90,8,10]] lattice is structurally immune."""
    cases = {b["label"]: b for b in E53.semisimple_baseline()["baselines"]}
    assert cases["[[90,8,10]]"]["case"] == "immune"
    assert cases["[[90,8,10]]"]["odd_lattice"] is True
    assert cases["[[90,8,10]]"]["power_dims"] == [4, 4]
    for label in ("[[72,12,6]]", "[[108,8,10]]", "[[144,12,12]]", "[[288,12,18]]",
                  "[[360,12,<=24]]", "[[756,16,<=34]]"):
        b = cases[label]
        assert b["case"] == "demote_full", label
        assert b["power_dims"][1] == 0, label      # I^2 = 0
        assert b["odd_lattice"] is False, label
    assert cases["[[144,12,12]]"]["k_parent"] == 12


def test_odd_lattice_corollary_on_fresh_parents() -> None:
    """Corollary J-G1: semisimple R => no demotion; test lattices not in the scan."""
    for ell, m in ((5, 7), (11, 3), (9, 7)):
        HX, _ = E53.bb_from_terms(ell, m, [(1, 0), (0, 1), (0, 2)], [(0, 1), (1, 0), (2, 0)])
        cls = E53.classify_ideal(E53.ideal_of_parent(HX, ell, m), ell, m)
        assert cls["case"] in ("immune", "degenerate"), (ell, m, cls)
        assert cls["dim_I_infty"] == cls["dim_I"], (ell, m, cls)


# ------------------------------------------------------------- lemma unit layer

def test_coset_law_blocks_exact_vanishing_at_weight_three() -> None:
    """J-H(3)/(4): a multi-coset weight-3 polynomial cannot vanish exactly."""
    ell, m = 4, 3
    multi = [(0, 0), (1, 0), (0, 1)]              # 2-parts {0,1}: two cosets
    assert len(E53.coset_support(multi, ell, m)) == 2
    HX, _ = E53.bb_from_terms(ell, m, multi, multi)
    assert E53.classify_ideal(E53.ideal_of_parent(HX, ell, m), ell, m)["dim_I_infty"] == 0


def test_scalar_law_gives_idempotent_annihilator() -> None:
    """J-H(4): a single-coset polynomial is zero-or-unit at every factor."""
    ell, m = 4, 3
    single = [(1, 0), (1, 1), (1, 2)]             # x*(1+y+y^2): one x-exponent
    assert len(E53.coset_support(single, ell, m)) == 1
    HX, _ = E53.bb_from_terms(ell, m, single, single)
    cls = E53.classify_ideal(E53.ideal_of_parent(HX, ell, m), ell, m)
    assert cls["dim_I_infty"] == cls["dim_I"] > 0   # Ann is idempotent


def test_coset_criterion_separates_the_catalogue(art53: dict) -> None:
    c = art53["coset_criterion"]
    assert c["criterion_exact"] is True
    assert c["exceptions"] == []
    assert c["weight_shapes"] == {"3,3": 202}
    assert c["table"] == {"demote_full|single_coset=False": 192,
                          "immune|single_coset=True": 10}


def test_lemma_falsification_battery_clean(art53: dict) -> None:
    lc = art53["lemma_checks"]
    assert lc["L1_holds"] and lc["L2_holds"]
    assert lc["L1_failures"] == [] and lc["L2_failures"] == []
    assert lc["tested"] >= 200


def test_dual_ideal_audit(art53: dict) -> None:
    """Ann_R(M) = (bar a, bar b) by two independent computations, plus the
    Frobenius identity dim I + dim Ann_R(M) = dim R, on every parent."""
    da = art53["duality_audit"]
    assert da["audit_holds"] and da["bad"] == []
    assert da["parents"] == 202
    bar_invariant = 0
    for r in da["records"]:
        assert r["ann_M_equals_reciprocal_ideal"], r["label"]
        assert r["frobenius_identity"], r["label"]
        assert r["dim_I"] + r["dim_ann_M"] == r["dim_R"], r["label"]
        assert r["zero_part_dim"] + r["nilpotent_part_dim"] == r["dim_I"], r["label"]
        if r["case"] == "immune":
            assert r["nilpotent_part_dim"] == 0, r["label"]
        else:
            assert r["zero_part_dim"] == 0, r["label"]
        bar_invariant += bool(r["ideal_AB_is_bar_invariant"])
    # the bar genuinely matters: comparing against (a,b) would fail on 46 parents
    assert bar_invariant == 156, bar_invariant


def test_the_two_annihilators_are_not_interchangeable(art53: dict) -> None:
    """Substituting Ann_R(M) into dim S = 2 dim I^infty must give a WRONG answer.
    This pins the mislabel that would otherwise pass unnoticed."""
    da = art53["duality_audit"]
    assert da["annihilators_interchangeable"] == []
    k_of = {r["label"]: r["k_parent"] for r in art53["catalogue"]["records"]}
    for r in da["records"]:
        assert r["wrong_route_dim_S"] != r["dim_S_correct_route"], r["label"]
        if r["case"] == "immune":
            assert r["dim_S_correct_route"] == k_of[r["label"]], r["label"]
            assert r["wrong_route_dim_S"] == 2 * (r["dim_R"] - r["dim_I"]), r["label"]
        else:
            assert r["dim_S_correct_route"] == 0, r["label"]


def test_dual_ideal_recomputed_on_the_witness() -> None:
    """In-test recomputation on the (2,3) witness: the wrong route would report
    demote-full instead of mixed, and Ann_R(M) is the reciprocal ideal."""
    from qec_research.gf2.linalg import rank_np

    ell, m = 2, 3
    A = [(0, 1), (0, 2), (1, 1), (1, 2)]
    B = [(0, 2), (0, 0), (1, 2), (1, 0)]
    HX, HZ = E53.bb_from_terms(ell, m, A, B)
    k_P = HX.shape[1] - rank_np(HX) - rank_np(HZ)
    aud = E53.duality_audit(HX, HZ, ell, m, k_P)
    assert aud["dim_S_correct_route"] == 4          # = 2 dim I^infty, matches J-J
    assert aud["wrong_route_dim_S"] == 0            # Ann_R(M) is nilpotent here
    assert aud["ann_M_equals_reciprocal_ideal"] and aud["frobenius_identity"]
    assert (aud["zero_part_dim"], aud["nilpotent_part_dim"]) == (2, 2)  # mixed
    ann_M = E53.module_annihilator(HX, HZ, ell, m, k_P)
    recip = E53.reciprocal_ideal(HZ, ell, m)
    assert rank_np(ann_M) == rank_np(recip) == 2
    assert rank_np(np.vstack([ann_M, recip])) == 2   # equal subspaces


def test_weight4_witness_is_tight(art53: dict) -> None:
    w4 = art53["lemma_checks"]["L3_weight4_witness"]
    assert w4["multi_coset"] and w4["vanishes_exactly"]   # tightness at weight 4


def test_ring_multiplication_matches_matrix_action() -> None:
    """poly_mult must agree with the circulant matrix action (convention lock)."""
    from qec_research.codes.bicycle import poly_matrix
    ell, m = 4, 3
    rng = np.random.default_rng(7)
    for _ in range(6):
        u = rng.integers(0, 2, ell * m, dtype=np.uint8)
        terms = [(1, 0), (0, 2), (3, 1)]
        want = (u @ poly_matrix(ell, m, terms)) % 2
        v = np.zeros(ell * m, np.uint8)
        for a, b in terms:
            v[a * m + b] ^= 1
        assert np.array_equal(E53.poly_mult(u, v, ell, m), want)


# ------------------------------------------------------------- census contract

def test_census_has_no_mixed_pair_anywhere(art54: dict) -> None:
    v = art54["verdict"]
    assert art54["max_weight"] == 3
    assert v["lattices_with_mixed"] == []
    assert v["crosscheck_mismatches"] == 0
    assert v["pairs_total"] > 6.5e8
    assert v["lattices_scanned"] == 18
    for name in ("9x6", "12x6", "15x6", "30x6", "12x12", "15x12", "6x6", "3x6", "6x3"):
        assert name in v["lattices_without_mixed"], name


def test_census_covers_every_catalogue_lattice(art54: dict) -> None:
    scanned = {(r["ell"], r["m"]) for r in art54["lattices"]}
    catalogue_lattices = {(3, 6), (6, 3), (6, 6), (9, 6), (12, 6), (15, 6), (30, 6)}
    assert catalogue_lattices <= scanned
    assert {(12, 12), (15, 12)} <= scanned          # the underexplored pair
    for r in art54["lattices"]:
        assert r["counts"]["mixed"] == 0
        assert r["crosscheck_mismatches"] == []
