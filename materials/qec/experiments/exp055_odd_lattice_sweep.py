"""EXP-055: exhaustive weight-<=3 sweep of CSS BB codes on ODD x ODD lattices.

WHY THIS REGION.  Theorem J-G1 (EXP-053/054) proves that when ell and m are both
odd, R is semisimple and the single-row X-collapse channel is structurally empty.
Our 202-parent PBB catalogue never searched the region (m in {3,6}), but the
broader BB literature DID: a 2026-08-21 novelty audit found at least 27 sourced
odd x odd instances.  The earlier claim that [[90,8,10]] was the only one is
retracted in notes/next_breakthroughs.md.

WHAT IS DECIDED HERE.

(1) Published rate law, independently reproduced.  Semisimplicity gives

        k_P = 2 dim I = 2 sum_{common-root orbits chi} [F_chi:F_2],
        I   = Ann_left(a,b).

    Route 1 computes dim I by GF(2) rank; route 2 computes
    k = n-rank(H_X)-rank(H_Z); on coprime lattices route 3 uses the published
    gcd formula.  The result is already in Panteleev--Kalachev Prop. 1,
    Lin--Pryadko Eq. 47, Wang--Mueller Eq. 11, Postema--Kokkelmans Thm. 2.6,
    and Eberhardt--Steffan Cor. 2.11--2.12.  It is not claimed as new.

(2) Exact reciprocal-pole logical transversal and certified bounds.  The stored
    I is a LEFT annihilator.  The physical right-kernel pole is J = bar(I):

        P = J x J <= ker H_X,    P cap S_Z = 0,    dim P = k,
        therefore P ~= ker H_X / S_Z.

    This is the Eberhardt--Steffan principal-code isomorphism in our matrix/bar
    convention.  It yields the solver-free ceiling

        d <= min{wt(u) : 0 != u in J},

    plus much tighter self-certifying witnesses by reducing selected pole
    vectors modulo S_Z.  Random information-set orders affect tightness only;
    every returned vector is verified in ker H_X outside S_Z.  The raw ceiling
    is a short corollary of published structure, not oversold as a deep theorem.

(3) Exhaustive census and Pareto screen.  Grouping polynomials by canon(Ann)
    collapses O(P^2) to O(G^2), exactly for k and the pole ceiling.  Translation
    normalisation is exhaustive up to independent qubit permutations.  The
    screen further canonicalises unit maps, block swap and x/y interchange,
    rejects with explicit witnesses first, and uses capped CP-SAT only on the
    residual.  Every screen survivor is exact-distance certified and checked for
    direct-sum decomposition.

Artifacts:
  results/processed/exp055_odd_lattice_sweep.json
  results/processed/exp055_literature_validation.json
  results/processed/exp055_odd_lattice_screen.json
  results/certificates/exp055_odd_lattice_survivors.json
  results/partial_runs/exp055/*.json
  results/partial_runs/exp055_screen/*.json

Run:
  python experiments/exp055_odd_lattice_sweep.py run --max-dim 180
  python experiments/exp055_odd_lattice_sweep.py literature
  python experiments/exp055_odd_lattice_sweep.py screen --min-n 18
  python experiments/exp055_odd_lattice_sweep.py screen-certify
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from functools import lru_cache
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


E53 = _load("exp053", "exp053_ideal_classification.py")

from qec_research.codes.bicycle import BRAVYI_BB, poly_matrix  # noqa: E402
from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.distance.sat_decide import (  # noqa: E402
    ENCODING_VERSION,
    css_side_instance,
    decide_weight_bounded,
    decision_cnf_digest,
    verify_witness_two_paths,
)
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

SCHEMA = "exp055-odd-lattice-sweep-v1"
OUT = ROOT / "results" / "processed" / "exp055_odd_lattice_sweep.json"
SHARDS = ROOT / "results" / "partial_runs" / "exp055"
SCREEN_DIR = ROOT / "results" / "partial_runs" / "exp055_screen"
LIT_OUT = ROOT / "results" / "processed" / "exp055_literature_validation.json"
SCREEN_OUT = ROOT / "results" / "processed" / "exp055_odd_lattice_screen.json"
SURVIVOR_CERT = ROOT / "results" / "certificates" / "exp055_odd_lattice_survivors.json"

MAX_DIM = 180            # ell*m <= MAX_DIM, i.e. n = 2*ell*m <= 360
K_MIN = 8                # the interesting regime; [[90,8,10]] has k = 8
K_CEIL_MAX = 24          # above this k the codes are degenerate; skip the ceiling
WITNESS_TRIES = 8
WITNESS_KEEP = 64
SCREEN_RESIDUAL_WORKERS = 8
DIM_I_ENUM_CAP = 16      # 2^16 vectors per enumeration keeps memory bounded
SCREEN_CDCL_CONFLICT_BUDGET = 1_000_000
SCREEN_CDCL_SOLVER = "cadical195"
REFERENCE_VALIDATION_VERSION = "exp055-exact-reference-validator-v11"
CERT_TIME_LIMIT_S = 900.0
CERT_WORKERS = 8


# --------------------------------------------------------------------------- #
# lattices
# --------------------------------------------------------------------------- #
def odd_lattices(max_dim: int = MAX_DIM) -> list[tuple[int, int]]:
    """All (ell, m) with ell >= m >= 3 both odd and ell*m <= max_dim.

    ell >= m is exhaustive: x <-> y interchange is a qubit permutation, so
    (ell, m) and (m, ell) carry the same set of codes.
    """
    out = []
    for m in range(3, max_dim + 1, 2):
        for ell in range(m, max_dim + 1, 2):
            if ell * m <= max_dim:
                out.append((ell, m))
    return sorted(out, key=lambda t: (t[0] * t[1], t[0]))


# --------------------------------------------------------------------------- #
# ring-level helpers
# --------------------------------------------------------------------------- #
def bar_permutation(ell: int, m: int) -> np.ndarray:
    """Index permutation of the involution x -> x^-1, y -> y^-1 on R."""
    idx = np.arange(ell * m).reshape(ell, m)
    return idx[(-np.arange(ell)) % ell][:, (-np.arange(m)) % m].reshape(-1)


def ann_basis(terms, ell: int, m: int) -> np.ndarray:
    """Ann(a) = {lambda : lambda a = 0} as a GF(2) row basis.

    Same convention as EXP-053/054: the left null space of poly_matrix(a).
    """
    M = poly_matrix(ell, m, [tuple(t) for t in terms])
    a = nullspace_np(M.T)
    if a.size == 0:
        return np.zeros((0, ell * m), np.uint8)
    r, _ = rref_np(a)
    return r[: rank_np(r)]


def canon(basis: np.ndarray) -> bytes:
    if basis.size == 0:
        return b""
    rr, _ = rref_np(basis)
    return rr[: rank_np(rr)].tobytes()


def intersect(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Basis of the intersection of two GF(2) row spaces (Zassenhaus)."""
    dx = rank_np(X) if X.size else 0
    dy = rank_np(Y) if Y.size else 0
    if dx == 0 or dy == 0:
        return np.zeros((0, X.shape[1] if X.size else Y.shape[1]), np.uint8)
    n = X.shape[1]
    top = np.hstack([X, X])
    bot = np.hstack([Y, np.zeros_like(Y)])
    rr, _ = rref_np(np.vstack([top, bot]))
    rows = [r for r in rr if not r[:n].any() and r[n:].any()]
    if not rows:
        return np.zeros((0, n), np.uint8)
    out = np.array([r[n:] for r in rows], np.uint8)
    rr2, _ = rref_np(out)
    return rr2[: rank_np(rr2)]


def span_vectors(basis: np.ndarray) -> np.ndarray:
    """All 2^d - 1 nonzero vectors of the span (d = rows of a reduced basis)."""
    d = basis.shape[0]
    coeffs = np.array(
        [[(c >> i) & 1 for i in range(d)] for c in range(1, 1 << d)], np.uint8
    )
    return (coeffs @ basis) % 2


def normalised_supports(ell: int, m: int, max_weight: int):
    """Supports containing (0,0) of weight 1..max_weight."""
    pts = [(a, b) for a in range(ell) for b in range(m) if (a, b) != (0, 0)]
    yield [(0, 0)]
    for w in range(1, max_weight):
        for extra in itertools.combinations(pts, w):
            yield [(0, 0), *extra]


# --------------------------------------------------------------------------- #
# the certified ceiling
# --------------------------------------------------------------------------- #
def certified_ceiling(I: np.ndarray, ell: int, m: int, bar: np.ndarray) -> dict:
    """Rigorous solver-free pole bound from the LEFT annihilator I.

    The coefficient representation used by ``nullspace(poly_matrix(a).T)`` is
    the left-annihilator convention.  A physical column vector lies in the
    right kernel of H_X only after the reciprocal involution:

        J := bar(I),       J x J <= ker H_X.

    On an odd lattice R is semisimple.  Character by character the Z-stabilizer
    block vanishes exactly on the support of J, so

        (J x J) cap S_Z = 0,     dim(J x J) = k,

    and J x J -> ker(H_X)/S_Z is an isomorphism (the principal-code logical
    structure of Eberhardt--Steffan, Cor. 2.11--2.12, in our bar convention).
    Thus every nonzero u in J makes (u,0) a nontrivial logical and

        d(P) <= min{wt(u) : 0 != u in J}.

    Bar preserves Hamming weight, so the numerical minimum equals min wt(I),
    but the physical witness MUST use J.  This distinction is regression-locked:
    using raw I fails on the published (7,7) [[98,6,12]] code.
    """
    d = I.shape[0]
    out: dict = {"dim_I": d, "k_parent": 2 * d}
    if d == 0:
        out.update({"dim_pole": 0, "ceiling": None, "bar_invariant": None})
        return out
    J = I[:, bar]
    out["dim_pole"] = d
    out["bar_invariant"] = bool(rank_np(np.vstack([I, J])) == d)
    if d <= DIM_I_ENUM_CAP:
        vecs = span_vectors(J)
        weights = vecs.sum(axis=1)
        idx = int(np.argmin(weights))
        out["ceiling"] = int(weights[idx])
        out["witness_support"] = [int(j) for j in np.flatnonzero(vecs[idx])]
        out["ceiling_is_exact_min"] = True
    else:
        rng = np.random.default_rng(0xC0FFEE ^ (d * 7919))
        combos = rng.integers(0, 2, size=(4096, d), dtype=np.uint8)
        cand = np.vstack([J, (combos @ J) % 2])
        w = cand.sum(axis=1)
        nz = np.flatnonzero(w > 0)
        idx = int(nz[int(np.argmin(w[nz]))])
        out["ceiling"] = int(w[idx])
        out["witness_support"] = [int(j) for j in np.flatnonzero(cand[idx])]
        out["ceiling_is_exact_min"] = False
    return out


def member_ceiling(A_terms, B_terms, ell: int, m: int) -> dict:
    """Directly verify the smallest reciprocal-pole witness for one BB code."""
    I = intersect(ann_basis(A_terms, ell, m), ann_basis(B_terms, ell, m))
    d = I.shape[0]
    if d == 0 or d > DIM_I_ENUM_CAP:
        return {"ceiling": None, "reason": "dim out of range"}
    bar = bar_permutation(ell, m)
    J = I[:, bar]                                  # physical right-kernel pole
    HX, HZ = E53.bb_from_terms(ell, m, A_terms, B_terms)
    Rz, _ = rref_np(HZ)
    Rz = Rz[: rank_np(Rz)]
    r0 = Rz.shape[0]
    n = 2 * ell * m
    vecs = span_vectors(J)
    for idx in np.argsort(vecs.sum(axis=1), kind="stable"):
        u = vecs[idx]
        full = np.zeros(n, np.uint8)
        full[: ell * m] = u
        if (HX @ full % 2).any():
            raise AssertionError("reciprocal-pole vector is not in ker H_X")
        if rank_np(np.vstack([Rz, full[None, :]])) > r0:
            return {
                "ceiling": int(u.sum()),
                "witness_support": [int(j) for j in np.flatnonzero(full)],
                "nontrivial_found": True,
            }
    return {"ceiling": None, "nontrivial_found": False}

def reduced_witness_bound(J: np.ndarray, HZ: np.ndarray, ell: int, m: int,
                          tries: int = WITNESS_TRIES, keep: int = WITNESS_KEEP,
                          seed: int = 0x55AA,
                          target: int | None = None) -> dict:
    """Tighten the pole ceiling by reducing explicit logicals modulo S_Z.

    Here J = bar(Ann(a,b)) is the physical right-kernel pole.  Every nonzero
    (u,v) in J x J is a nontrivial logical on an odd lattice.  Reducing it against
    rref(H_Z) produces a concrete representative of the same class; ANY nonzero
    vector returned is therefore a self-certifying logical witness and its weight
    is a rigorous upper bound on d.

    Extra random information sets only lighten representatives of already
    nontrivial classes, so soundness never depends on the randomisation.
    """
    d = J.shape[0]
    n = 2 * ell * m
    if d == 0 or d > DIM_I_ENUM_CAP:
        return {"bound": None, "classes_tested": 0}
    vecs = span_vectors(J)
    order = np.argsort(vecs.sum(axis=1), kind="stable")[:keep]
    selected = vecs[order]
    count = len(selected)
    # Both one-block slices and the Cartesian grid of light pole pairs.
    targets = np.zeros((2 * count + count * count, n), np.uint8)
    targets[:count, : ell * m] = selected
    targets[count:2 * count, ell * m:] = selected
    targets[2 * count:, : ell * m] = np.repeat(selected, count, axis=0)
    targets[2 * count:, ell * m:] = np.tile(selected, (count, 1))

    HZ = np.asarray(HZ, np.uint8) & 1

    def reduce_against(perm: np.ndarray | None) -> np.ndarray:
        M = HZ if perm is None else HZ[:, perm]
        R, _ = rref_np(M)
        R = R[: rank_np(R)]
        T = targets if perm is None else targets[:, perm]
        T = T.copy()
        for row in R:
            j = int(np.flatnonzero(row)[0])
            hit = T[:, j] == 1
            if hit.any():
                T[hit] ^= row
        return T

    base = reduce_against(None)
    nontrivial = base.any(axis=1)                       # class is nontrivial
    if not nontrivial.any():
        return {"bound": None, "classes_tested": int(targets.shape[0]),
                "all_trivial": True}
    idxs = np.flatnonzero(nontrivial)
    base_weights = base[idxs].sum(axis=1)
    j0 = int(idxs[int(np.argmin(base_weights))])
    best_vec = base[j0].copy()
    best = int(best_vec.sum())
    rng = np.random.default_rng(seed ^ (ell * 131 + m))
    if target is not None and best <= target:
        return {
            "bound": best,
            "witness_support": [int(j) for j in np.flatnonzero(best_vec)],
            "classes_tested": int(nontrivial.sum()),
            "all_trivial": False,
        }
    for _ in range(tries):
        perm = rng.permutation(n)
        T = reduce_against(perm)
        candidate_rows = T[idxs]
        weights = candidate_rows.sum(axis=1)
        j = int(np.argmin(weights))
        if int(weights[j]) < best:
            best = int(weights[j])
            best_vec = np.zeros(n, np.uint8)
            best_vec[perm] = candidate_rows[j]
        if target is not None and best <= target:
            break
    return {
        "bound": best,
        "witness_support": [int(j) for j in np.flatnonzero(best_vec)],
        "classes_tested": int(nontrivial.sum()),
        "all_trivial": False,
    }


# --------------------------------------------------------------------------- #
# literature validation battery
# --------------------------------------------------------------------------- #
# Sourced odd x odd BB instances with their REPORTED distances.  Provenance is
# classified per row below: Wang--Mueller uses BP-OSD upper bounds and Postema
# uses Monte-Carlo estimates; only ``distance_exact_certified_here`` rows may
# enter the domination reference set.  "pi" means the coprime form pi=xy.
# Recorded from the primary texts and checked 2026-08-22.
LITERATURE_ODD: list[dict] = [
    {"ell": 15, "m": 3, "k": 8, "d": 10, "A": [(9, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (2, 0), (7, 0)], "source": "2308.07915 Table 3",
     "distance_exact_certified_here": True, "certificate_kind": "legacy_css_exact",
     "distance_certificate": "results/certificates/bb_distance_90_8_10.json"},
    {"ell": 3, "m": 3, "k": 4, "d": 4, "A": [(0, 0), (1, 0), (0, 1)],
     "B": [(0, 0), (2, 0), (0, 2)], "source": "2408.10001v4 App.B Table 3"},
    {"ell": 3, "m": 3, "k": 4, "d": 2, "A": [(0, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (1, 0), (0, 1)], "source": "2502.17052v4 Table 2",
     "source_distance_estimate": True, "distance_exact_certified_here": True,
     "certificate_kind": "legacy_css_exact",
     "distance_certificate": "results/certificates/exp055_postema_rows.json"},
    {"ell": 3, "m": 3, "k": 8, "d": 2, "A": [(0, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (1, 0), (2, 0)], "source": "2502.17052v4 Table 2",
     "source_distance_estimate": True, "distance_exact_certified_here": True,
     "certificate_kind": "legacy_css_exact",
     "distance_certificate": "results/certificates/exp055_postema_rows.json"},
    {"ell": 3, "m": 5, "k": 4, "d": 6, "pi_A": [0, 1, 2], "pi_B": [1, 3, 8],
     "source": "2408.10001v4 Table 2"},
    {"ell": 3, "m": 7, "k": 6, "d": 6, "pi_A": [0, 2, 3], "pi_B": [1, 3, 11],
     "source": "2408.10001v4 Table 2"},
    {"ell": 3, "m": 7, "k": 10, "d": 4, "pi_A": [0, 1, 5], "pi_B": [0, 2, 10],
     "source": "2408.10001v4 App.C Table 4"},
    {"ell": 3, "m": 9, "k": 8, "d": 6, "A": [(0, 0), (0, 2), (0, 4)],
     "B": [(0, 3), (1, 0), (2, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 3, "m": 9, "k": 4, "d": 8, "A": [(1, 0), (0, 1), (0, 3)],
     "B": [(0, 0), (0, 2), (2, 0)], "source": "2408.10001v4 App.B Table 3"},
    {"ell": 3, "m": 11, "k": 4, "d": 10, "pi_A": [0, 1, 5], "pi_B": [0, 1, 23],
     "source": "2408.10001v4 App.C Table 4"},
    {"ell": 5, "m": 7, "k": 6, "d": 8, "pi_A": [0, 1, 5], "pi_B": [0, 1, 12],
     "source": "2408.10001v4 Table 2"},
    {"ell": 5, "m": 9, "k": 4, "d": 12, "pi_A": [0, 1, 4], "pi_B": [0, 8, 34],
     "source": "2408.10001v4 App.C Table 4", "unreproduced": True},
    {"ell": 5, "m": 9, "k": 8, "d": 8, "pi_A": [0, 1, 12], "pi_B": [0, 2, 9],
     "source": "2408.10001v4 App.C Table 4"},
    {"ell": 7, "m": 7, "k": 6, "d": 12, "A": [(1, 0), (0, 3), (0, 4)],
     "B": [(0, 1), (3, 0), (4, 0)], "source": "2407.03973v1 Table 1"},
    {"ell": 7, "m": 7, "k": 6, "d": 12, "A": [(3, 0), (0, 5), (0, 6)],
     "B": [(0, 2), (3, 0), (5, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 7, "m": 7, "k": 6, "d": 8, "A": [(4, 0), (0, 1), (0, 3)],
     "B": [(0, 4), (1, 0), (3, 0)], "source": "2502.17052v4 Table 2",
     "source_distance_estimate": True, "distance_exact_certified_here": True,
     "certificate_kind": "legacy_css_exact",
     "distance_certificate": "results/certificates/exp055_postema_rows.json"},
    {"ell": 7, "m": 9, "k": 12, "d": 10, "pi_A": [0, 1, 58], "pi_B": [3, 16, 44],
     "source": "2408.10001v4 Table 2", "source_distance_estimate": True,
     "distance_exact_certified_here": True, "certificate_kind": "legacy_css_exact",
     "distance_certificate": "results/certificates/exp055_discovered_references.json"},
    {"ell": 7, "m": 9, "k": 6, "d": 14, "pi_A": [0, 4, 19], "pi_B": [0, 6, 16],
     "source": "2408.10001v4 App.C Table 4"},
    {"ell": 3, "m": 21, "k": 8, "d": 10, "A": [(0, 0), (0, 2), (0, 10)],
     "B": [(0, 3), (1, 0), (2, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 5, "m": 15, "k": 16, "d": 8, "A": [(0, 0), (0, 6), (0, 8)],
     "B": [(0, 5), (1, 0), (4, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 3, "m": 27, "k": 8, "d": 14, "A": [(0, 0), (0, 10), (0, 14)],
     "B": [(0, 12), (1, 0), (2, 0)], "source": "2408.10001v4 Table 1",
     "source_distance_estimate": True, "distance_exact_certified_here": True,
     "certificate_kind": "exp056_odd_exact",
     "distance_certificate":
         "results/certificates/exp056_wm_162_8_14_distance.json"},
    {"ell": 7, "m": 11, "k": 6, "d": 16, "pi_A": [0, 4, 31], "pi_B": [0, 19, 53],
     "source": "2408.10001v4 App.C Table 4", "unreproduced": True},
    {"ell": 9, "m": 9, "k": 4, "d": 16, "A": [(0, 0), (1, 0), (0, 1)],
     "B": [(3, 0), (0, 1), (0, 2)], "source": "2407.03973v1 Table 1"},
    {"ell": 9, "m": 9, "k": 12, "d": 8, "A": [(0, 0), (1, 0), (0, 6)],
     "B": [(0, 3), (2, 0), (3, 0)], "source": "2407.03973v1 Table 1"},
    {"ell": 9, "m": 9, "k": 24, "d": 6, "A": [(0, 0), (0, 1), (0, 2)],
     "B": [(0, 3), (3, 0), (6, 0)], "source": "2407.03973v1 Table 1"},
    {"ell": 9, "m": 9, "k": 8, "d": 12, "A": [(3, 0), (0, 1), (0, 2)],
     "B": [(0, 3), (1, 0), (2, 0)], "source": "2407.03973v1 Table 1"},
    {"ell": 9, "m": 15, "k": 8, "d": 18, "A": [(3, 0), (0, 1), (0, 2)],
     "B": [(0, 3), (1, 0), (2, 0)], "source": "2407.03973v1 Table 1"},
    {"ell": 39, "m": 3, "k": 8, "d": 18,
     "A": [(0, 0), (1, 0), (5, 0)],
     "B": [(0, 0), (1, 1), (23, 2)],
     "source": "2503.03827v3 Table III row 2",
     "source_distance_exact_claim": True,
     "distance_exact_certified_here": True,
     "certificate_kind": "exp067_odd_exact",
     "distance_certificate":
         "results/certificates/exp067_234_8_18_bundle19_distance.json"},
    {"ell": 39, "m": 3, "k": 8, "d": 18,
     "A": [(0, 0), (1, 0), (5, 0)],
     "B": [(0, 0), (2, 2), (22, 1)],
     "source": "2503.03827v3 Table III row 1",
     "source_distance_exact_claim": True,
     "distance_exact_certified_here": True,
     "certificate_kind": "exp067_odd_exact",
     "distance_certificate":
         "results/certificates/exp067_234_8_18_bundle22_distance.json"},
]


def _certificate_stat(path: Path) -> tuple[int, int, int, int, int]:
    stat = path.stat()
    return (
        int(stat.st_dev), int(stat.st_ino), int(stat.st_size),
        int(stat.st_mtime_ns), int(stat.st_ctime_ns),
    )


@lru_cache(maxsize=None)
def _read_distance_certificate(
    relative_path: str, snapshot: tuple[int, int, int, int, int]
) -> tuple[dict, str]:
    path = ROOT / relative_path
    payload = path.read_bytes()
    if _certificate_stat(path) != snapshot:
        raise RuntimeError(f"distance certificate changed while reading: {relative_path}")
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


def _distance_certificate_snapshot(relative_path: str) -> tuple[dict, str]:
    path = ROOT / relative_path
    if not path.exists():
        raise RuntimeError(f"required distance certificate missing: {relative_path}")
    return _read_distance_certificate(relative_path, _certificate_stat(path))


@lru_cache(maxsize=1)
def _exp056_validator_module():
    return _load("exp056_certificate_validator", "exp056_odd_distance.py")


@lru_cache(maxsize=None)
def _validated_exp056_snapshot(
    relative_path: str, snapshot: tuple[int, int, int, int, int]
) -> tuple[dict, str]:
    certificate, certificate_sha256 = _read_distance_certificate(
        relative_path, snapshot
    )
    _exp056_validator_module().validate_exact_certificate_payload(certificate)
    return certificate, certificate_sha256


def _validated_distance_reference(rec: dict) -> tuple[bool, str | None]:
    if not rec.get("distance_exact_certified_here", False):
        return False, None
    relative = rec.get("distance_certificate")
    if relative is None:
        raise RuntimeError("exact literature reference has no distance certificate")
    path = ROOT / str(relative)
    if not path.exists():
        raise RuntimeError(f"required distance certificate missing: {relative}")
    kind = rec.get("certificate_kind")
    A, B = _terms_of(rec)
    if kind == "legacy_css_exact":
        normalized = {
            **rec,
            "n": 2 * int(rec["ell"]) * int(rec["m"]),
            "A": [list(term) for term in A],
            "B": [list(term) for term in B],
        }
        return _validated_exact_reference(normalized)
    if kind == "exp067_odd_exact":
        normalized = {
            **rec,
            "n": 2 * int(rec["ell"]) * int(rec["m"]),
            "A": [list(term) for term in A],
            "B": [list(term) for term in B],
        }
        return _validated_exact_reference(normalized)
    if kind != "exp056_odd_exact":
        raise RuntimeError(f"unsupported literature certificate kind: {kind}")
    certificate, certificate_sha256 = _validated_exp056_snapshot(
        str(relative), _certificate_stat(path)
    )
    verdict = certificate.get("verdict", {})
    target = certificate.get("target", {})
    identity = certificate.get("identity", {})
    lower = certificate.get("lower_bound", {}).get("class_route", {})
    upper = certificate.get("upper_bound", {})
    HX, HZ = E53.bb_from_terms(rec["ell"], rec["m"], A, B)
    matrix_sha = lambda matrix: hashlib.sha256(
        np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1).tobytes()
    ).hexdigest()
    valid = bool(
        certificate.get("schema") == "exp056-odd-orbit-distance-v1"
        and verdict.get("exact") is True
        and int(verdict.get("d", -1)) == int(rec["d"])
        and int(verdict.get("d_X", -1)) == int(rec["d"])
        and int(verdict.get("d_Z", -1)) == int(rec["d"])
        and int(target.get("ell", -1)) == int(rec["ell"])
        and int(target.get("m", -1)) == int(rec["m"])
        and int(target.get("expected_k", -1)) == int(rec["k"])
        and target.get("A") == [list(term) for term in A]
        and target.get("B") == [list(term) for term in B]
        and identity.get("HX_sha256") == matrix_sha(HX)
        and identity.get("HZ_sha256") == matrix_sha(HZ)
        and lower.get("all_classes_unsat") is True
        and lower.get("all_classes_replayed") is True
        and int(upper.get("weight", -1)) == int(rec["d"])
        and upper.get("z_valid_numpy") is True
        and upper.get("z_valid_bitset") is True
        and upper.get("x_valid_numpy") is True
        and upper.get("x_valid_bitset") is True
    )
    if not valid:
        raise RuntimeError(
            f"distance certificate does not bind [[{2 * rec['ell'] * rec['m']},"
            f"{rec['k']},{rec['d']}]]: {relative}"
        )
    return True, certificate_sha256


def _distance_exact_here(rec: dict) -> bool:
    return _validated_distance_reference(rec)[0]

def _exp027_matrix_fingerprint(*matrices: np.ndarray) -> str:
    """Hash matrices with the exact EXP-027/037 shape-aware convention."""
    digest = hashlib.sha256()
    for matrix in matrices:
        array = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1)
        digest.update(len(array.shape).to_bytes(1, "big"))
        for dimension in array.shape:
            digest.update(int(dimension).to_bytes(8, "big"))
        digest.update(array.tobytes())
    return digest.hexdigest()

def _legacy_witness_valid(
    HX: np.ndarray,
    HZ: np.ndarray,
    support: list[int] | None,
    side: str,
    distance: int,
) -> bool:
    if support is None:
        return False
    vector = np.zeros(HX.shape[1], dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    if int(vector.sum()) != int(distance):
        return False
    if side == "x":
        constraints, stabilizers = HZ, HX
    else:
        constraints, stabilizers = HX, HZ
    return bool(
        not np.any(constraints @ vector % 2)
        and rank_np(np.vstack([stabilizers, vector])) == rank_np(stabilizers) + 1
    )


def validate_legacy_reference_payload(rec: dict, certificate: dict) -> None:
    """Bind aggregate/legacy CSS exact records to rebuilt matrices and witnesses."""
    ell, m = int(rec["ell"]), int(rec["m"])
    n, k, distance = int(rec["n"]), int(rec["k"]), int(rec["d"])
    A = [list(map(int, term)) for term in rec["A"]]
    B = [list(map(int, term)) for term in rec["B"]]
    candidates = certificate.get("records", [certificate])
    matching = [
        row
        for row in candidates
        if int(row.get("ell", -1)) == ell
        and int(row.get("m", -1)) == m
        and int(row.get("n", -1)) == n
        and int(row.get("k", -1)) == k
        and int(row.get("d", -1)) == distance
        and row.get("A") == A
        and row.get("B") == B
    ]
    if len(matching) != 1:
        raise RuntimeError(
            f"legacy distance certificate does not uniquely bind [[{n},{k},{distance}]]"
        )
    row = matching[0]
    HX, HZ = E53.bb_from_terms(ell, m, A, B)
    rx, rz = rank_np(HX), rank_np(HZ)
    schema = certificate.get("schema")
    if schema == "exp055-discovered-references-v1":
        exact_flags = bool(
            certificate.get("all_exact") is True
            and row.get("d_exact") is True
            and row.get("d_X_exact") is True
            and row.get("d_Z_exact") is True
        )
        x_support, z_support = row.get("d_X_witness"), row.get("d_Z_witness")
        witness_gate = bool(
            _legacy_witness_valid(HX, HZ, x_support, "x", distance)
            and _legacy_witness_valid(HX, HZ, z_support, "z", distance)
        )
    elif schema == "exp055-postema-rows-v1":
        exact_flags = bool(
            certificate.get("all_exact") is True and row.get("exact") is True
        )
        # These historical CP-SAT records predate persisted witnesses.
        witness_gate = True
    else:
        exact_flags = bool(
            row.get("CERTIFIED_EXACT") is True
            and row.get("solver_reports_exact") is True
            and row.get("witnesses_independently_verified") is True
            and row.get("all_sectors_decided") == {"X": True, "Z": True}
            and int(row.get("d_lower_bound", -1)) == distance
        )
        witness_gate = bool(
            _legacy_witness_valid(
                HX, HZ, row.get("witness_X", {}).get("support"), "x", distance
            )
            and _legacy_witness_valid(
                HX, HZ, row.get("witness_Z", {}).get("support"), "z", distance
            )
        )
    valid = bool(
        n == 2 * ell * m
        and n - rx - rz == k
        and int(row.get("d_X", -1)) == distance
        and int(row.get("d_Z", -1)) == distance
        and exact_flags
        and witness_gate
    )
    if not valid:
        raise RuntimeError(
            f"legacy distance certificate does not validate [[{n},{k},{distance}]]"
        )


@lru_cache(maxsize=None)
def _validated_legacy_reference_snapshot(
    relative_path: str,
    snapshot: tuple[int, int, int, int, int],
    reference_json: str,
) -> tuple[dict, str]:
    certificate, certificate_sha256 = _read_distance_certificate(
        relative_path, snapshot
    )
    validate_legacy_reference_payload(json.loads(reference_json), certificate)
    return certificate, certificate_sha256


def validate_exp037_reference_payload(rec: dict, certificate: dict) -> None:
    """Rebuild and bind one legacy EXP-037 exact CSS distance certificate."""
    try:
        ell, m = int(rec["ell"]), int(rec["m"])
        n, k, distance = int(rec["n"]), int(rec["k"]), int(rec["d"])
        A = [tuple(map(int, term)) for term in rec["A"]]
        B = [tuple(map(int, term)) for term in rec["B"]]
        source_index = int(rec["source_catalogue_index"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("EXP-037 reference identity is incomplete") from exc

    HX, HZ = E53.bb_from_terms(ell, m, A, B)
    instances = {
        side: css_side_instance(HX, HZ, side, block_length=ell * m)
        for side in ("x", "z")
    }
    lower_hashes = {
        side: decision_cnf_digest(instances[side], distance - 1)
        for side in ("x", "z")
    }
    calls = certificate.get("calls", [])
    lower_bound_matches = all(
        any(
            call.get("side") == side
            and int(call.get("cap", -1)) == distance - 1
            and call.get("status") == "UNSAT"
            and call.get("cnf_sha256") == lower_hashes[side]
            for call in calls
        )
        for side in ("x", "z")
    )

    witness_side = certificate.get("witness_side")
    witness = np.asarray(certificate.get("witness_vector", []), dtype=np.uint8)
    witness_check = (
        verify_witness_two_paths(instances[witness_side], witness)
        if witness_side in instances and witness.shape == (n,)
        else {"valid": False}
    )
    rx, rz = rank_np(HX), rank_np(HZ)
    valid = bool(
        certificate.get("schema") == "exp037-envelope-classification-v1"
        and certificate.get("encoding_version") == ENCODING_VERSION
        and certificate.get("exact") is True
        and int(certificate.get("n", -1)) == n == 2 * ell * m
        and int(certificate.get("k", -1)) == k == n - rx - rz
        and int(certificate.get("distance", -1)) == distance
        and int(certificate.get("certified_upper_bound", -1)) == distance
        and int(certificate.get("rank_HX", -1)) == rx
        and int(certificate.get("rank_HZ", -1)) == rz
        and int(certificate.get("source_catalogue_index", -1)) == source_index
        and certificate.get("css_fingerprint")
        == _exp027_matrix_fingerprint(HX, HZ)
        and lower_bound_matches
        and witness_check.get("valid") is True
        and int(witness.sum()) == distance
    )
    if not valid:
        raise RuntimeError(
            f"EXP-037 distance certificate does not bind [[{n},{k},{distance}]]"
        )


@lru_cache(maxsize=None)
def _validated_exp037_reference_snapshot(
    relative_path: str,
    snapshot: tuple[int, int, int, int, int],
    reference_json: str,
) -> tuple[dict, str]:
    certificate, certificate_sha256 = _read_distance_certificate(
        relative_path, snapshot
    )
    validate_exp037_reference_payload(json.loads(reference_json), certificate)
    return certificate, certificate_sha256

@lru_cache(maxsize=None)
def _frontier_validator_module(module_name: str, filename: str):
    return _load(module_name, filename)


@lru_cache(maxsize=None)
def _validated_frontier_reference_snapshot(
    relative_path: str,
    snapshot: tuple[int, int, int, int, int],
    reference_json: str,
    label: str,
    module_file: str,
) -> tuple[dict, str]:
    certificate, certificate_sha256 = _read_distance_certificate(
        relative_path, snapshot
    )
    _frontier_validator_module(
        f"{label.lower()}_certificate_validator", module_file
    ).validate_exact_certificate_payload(certificate)
    rec = json.loads(reference_json)
    identity = certificate.get("identity", {})
    target = identity.get("target", {})
    valid = bool(
        int(identity.get("n", -1)) == int(rec["n"])
        and int(identity.get("k", -1)) == int(rec["k"])
        and int(target.get("expected_d", -1)) == int(rec["d"])
        and int(target.get("ell", -1)) == int(rec["ell"])
        and int(target.get("m", -1)) == int(rec["m"])
        and target.get("A") == rec["A"]
        and target.get("B") == rec["B"]
    )
    if not valid:
        raise RuntimeError(
            f"{label} distance certificate does not bind "
            f"[[{rec['n']},{rec['k']},{rec['d']}]]"
        )
    return certificate, certificate_sha256


FRONTIER_REFERENCE_KINDS = {
    "exp057_odd_exact": ("EXP-057", "exp057_odd_frontier.py"),
    "exp058_odd_exact": ("EXP-058", "exp058_odd_frontier.py"),
    "exp060_odd_exact": ("EXP-060", "exp060_frontier_ratchet.py"),
    "exp064_odd_exact": ("EXP-064", "exp064_n210_promotions.py"),
    "exp067_odd_exact": ("EXP-067", "exp067_n234_connected_cluster.py"),
}


def _validated_exact_reference(rec: dict) -> tuple[bool, str | None]:
    relative = rec.get("distance_certificate")
    if relative is None:
        raise RuntimeError("exact reference has no distance certificate")
    certificate_kind = rec.get("certificate_kind")
    path = ROOT / str(relative)
    if not path.exists():
        raise RuntimeError(f"required distance certificate missing: {relative}")
    reference_json = json.dumps(rec, sort_keys=True, separators=(",", ":"))
    if certificate_kind == "legacy_css_exact":
        _certificate, certificate_sha256 = _validated_legacy_reference_snapshot(
            str(relative), _certificate_stat(path), reference_json
        )
        return True, certificate_sha256
    if certificate_kind == "exp037_css_exact":
        _certificate, certificate_sha256 = _validated_exp037_reference_snapshot(
            str(relative), _certificate_stat(path), reference_json
        )
        return True, certificate_sha256
    if certificate_kind not in FRONTIER_REFERENCE_KINDS:
        raise RuntimeError(f"unsupported exact-reference certificate: {relative}")
    label, module_file = FRONTIER_REFERENCE_KINDS[certificate_kind]
    _certificate, certificate_sha256 = _validated_frontier_reference_snapshot(
        str(relative), _certificate_stat(path), reference_json,
        label, module_file,
    )
    return True, certificate_sha256


def _terms_of(rec: dict) -> tuple[list, list]:
    """Polynomial supports, converting the pi = xy coprime form when present."""
    ell, m = rec["ell"], rec["m"]
    if "pi_A" in rec:
        return ([(e % ell, e % m) for e in rec["pi_A"]],
                [(e % ell, e % m) for e in rec["pi_B"]])
    return rec["A"], rec["B"]


def _gf2_deg(p: int) -> int:
    return p.bit_length() - 1


def _gf2_mod(a: int, b: int) -> int:
    db = _gf2_deg(b)
    while a and _gf2_deg(a) >= db:
        a ^= b << (_gf2_deg(a) - db)
    return a


def _gf2_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, _gf2_mod(a, b)
    return a


def cyclic_k(terms, other, ell: int, m: int) -> int | None:
    """k by the PUBLISHED gcd formula, valid when gcd(ell, m) = 1.

    Panteleev-Kalachev Prop. 1 / Wang-Mueller Eq. (11): a coprime BB code is a
    generalised bicycle code over GF(2)[pi]/(pi^N - 1) with N = ell*m, and
    k = 2 deg gcd(a(pi), b(pi), pi^N - 1).  Shares no code with the annihilator
    route, so it is a genuinely independent third check.
    """
    from math import gcd as _g
    if _g(ell, m) != 1:
        return None
    N = ell * m
    inv_m, inv_l = pow(m, -1, ell), pow(ell, -1, m)

    def to_pi(ts) -> int:
        v = 0
        for (a, b) in ts:                      # CRT: pi^e <-> x^(e mod l) y^(e mod m)
            e = (a * m * inv_m + b * ell * inv_l) % N
            v ^= 1 << e
        return v

    g = _gf2_gcd(_gf2_gcd(to_pi(terms), to_pi(other)), (1 << N) | 1)
    return 2 * _gf2_deg(g)


def validate_literature(args: argparse.Namespace) -> int:
    """Independent checks against every sourced published odd x odd instance.

    (1) k by three routes where available: annihilator dimension, matrix ranks,
        and (coprime lattices) the published gcd formula.
    (2) exact pole isomorphism: J=bar(I), J x J <= ker H_X, intersection with
        S_Z zero, and dim(J x J)=k.
    (3) the pole ceiling must be >= the published exact d.  One violation
        falsifies the claimed upper bound.
    """
    recs = []
    for rec in LITERATURE_ODD:
        ell, m = rec["ell"], rec["m"]
        A, B = _terms_of(rec)
        n = 2 * ell * m
        bar = bar_permutation(ell, m)
        HX, HZ = E53.bb_from_terms(ell, m, A, B)
        k_mat = int(n - rank_np(HX) - rank_np(HZ))
        I = intersect(ann_basis(A, ell, m), ann_basis(B, ell, m))
        J = I[:, bar]
        k_ideal = 2 * I.shape[0]
        k_gcd = cyclic_k(A, B, ell, m)
        cel = certified_ceiling(I, ell, m, bar)
        ceiling = cel.get("ceiling")
        pole = np.vstack([
            np.hstack([J, np.zeros_like(J)]),
            np.hstack([np.zeros_like(J), J]),
        ])
        pole_rank = rank_np(pole)
        rank_hz = rank_np(HZ)
        pole_in_kernel = not (HX @ pole.T % 2).any()
        pole_intersection_zero = (
            rank_np(np.vstack([HZ, pole])) == rank_hz + pole_rank
        )
        pole_isomorphism = bool(
            pole_in_kernel and pole_intersection_zero
            and pole_rank == k_mat == k_ideal
        )
        routes_agree = k_ideal == k_mat and (k_gcd is None or k_gcd == k_ideal)
        exact_here, distance_certificate_sha256 = _validated_distance_reference(rec)
        source_estimate = bool(
            rec.get("source_distance_estimate", False)
            or rec["source"].startswith("2408.10001")
            or rec["source"].startswith("2502.17052")
        )
        out = {
            "source_distance_estimate": source_estimate,
            "distance_exact_certified_here": exact_here,
            "distance_certificate_sha256": distance_certificate_sha256,
            "source": rec["source"], "ell": ell, "m": m, "n": n,
            "k_published": rec["k"], "k_ideal": k_ideal, "k_matrices": k_mat,
            "k_gcd_formula": k_gcd, "our_routes_agree": bool(routes_agree),
            "k_agrees_with_paper": bool(routes_agree and k_ideal == rec["k"]),
            "pole_rank": pole_rank, "pole_in_kernel": pole_in_kernel,
            "pole_intersection_zero": pole_intersection_zero,
            "pole_isomorphism": pole_isomorphism,
            "d_published": rec["d"], "ceiling": ceiling, "ceiling_route": "pole",
            "dim_pole": cel.get("dim_pole"), "bar_invariant": cel.get("bar_invariant"),
            "ceiling_respects_d": (None if ceiling is None else
                                   bool(ceiling >= rec["d"])),
            "slack": (None if ceiling is None else int(ceiling - rec["d"])),
        }
        recs.append(out)
        print(json.dumps({k: out[k] for k in
                          ("source", "ell", "m", "k_published", "k_ideal",
                           "k_gcd_formula", "k_agrees_with_paper", "d_published",
                           "ceiling", "ceiling_route", "ceiling_respects_d")}),
              flush=True)
    # A row whose transcribed exponents give k = 0 by BOTH of our independent
    # routes is a transcription/printing problem, not a convention problem: it is
    # reported separately and excluded from the pass/fail count.
    unreproducible = [r for r in recs
                      if r["our_routes_agree"] and not r["k_agrees_with_paper"]]
    core = [r for r in recs if r not in unreproducible]
    reported_bad = [r for r in core if r["ceiling_respects_d"] is False]
    reported_slacks = [r["slack"] for r in core if r["slack"] is not None]
    exact_core = [r for r in core if r["distance_exact_certified_here"]]
    exact_bad = [r for r in exact_core if r["ceiling_respects_d"] is False]
    exact_slacks = [r["slack"] for r in exact_core if r["slack"] is not None]
    payload = {
        "schema": "exp055-literature-v3", "utc": E53.E52.utc_now(),
        "instances": len(recs), "records": recs,
        "reproduced_instances": len(core),
        "k_all_agree_on_reproduced": all(r["k_agrees_with_paper"] for r in core),
        "our_routes_always_agree": all(r["our_routes_agree"] for r in recs),
        "pole_isomorphism_all": all(r["pole_isomorphism"] for r in recs),
        "unreproducible_rows": unreproducible,
        "reported_distance_instances": len(core),
        "reported_ceiling_violations": reported_bad,
        "reported_ceiling_sanity_holds": not reported_bad,
        "reported_slack_min": min(reported_slacks) if reported_slacks else None,
        "reported_slack_median": (
            sorted(reported_slacks)[len(reported_slacks) // 2]
            if reported_slacks else None
        ),
        "reported_slack_max": max(reported_slacks) if reported_slacks else None,
        "exact_distance_instances": len(exact_core),
        "exact_ceiling_violations": exact_bad,
        "exact_ceiling_never_violated": not exact_bad,
        "exact_slack_min": min(exact_slacks) if exact_slacks else None,
        "exact_slack_median": (
            sorted(exact_slacks)[len(exact_slacks) // 2] if exact_slacks else None
        ),
        "exact_slack_max": max(exact_slacks) if exact_slacks else None,
        "no_certificate": [r["source"] for r in core if r["ceiling"] is None],
    }
    E53.E52.atomic_write_json(LIT_OUT, payload)
    print(json.dumps({k: payload[k] for k in (
        "instances", "reproduced_instances", "k_all_agree_on_reproduced",
        "our_routes_always_agree", "pole_isomorphism_all",
        "reported_ceiling_sanity_holds", "reported_distance_instances",
        "exact_ceiling_never_violated", "exact_distance_instances",
        "exact_slack_min", "exact_slack_median", "exact_slack_max",
        "no_certificate",
    )}, indent=1))
    return 0

# --------------------------------------------------------------------------- #
# candidate enumeration, symmetry reduction, domination screen
# --------------------------------------------------------------------------- #
# Local two-sided certificates only. EXP-055 contributes exact
# $[[30,8,4]]$ and $[[54,8,6]]$ references; the four standard baselines and
# every promoted reference must have a local certificate. EXP-057/058/060
# supply the new $[[170,16,10]]$, $[[186,10,14]]$, and $[[210,18,8]]$ fixed
# points; EXP-037's $n=180$ frontier is rebound to its constructors, current
# CNFs, and physical witnesses.
EXACT_REFERENCES = [
    {
        "name": "EXP-055 [[30,8,4]]", "n": 30, "k": 8, "d": 4,
        "ell": 5, "m": 3,
        "A": [[0, 0], [1, 0], [3, 1]],
        "B": [[0, 0], [1, 1], [4, 1]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/exp055_discovered_references.json",
    },
    {
        "name": "EXP-055 [[54,8,6]]", "n": 54, "k": 8, "d": 6,
        "ell": 9, "m": 3,
        "A": [[0, 0], [0, 1], [3, 2]],
        "B": [[0, 0], [1, 0], [2, 0]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/exp055_discovered_references.json",
    },
    {
        "name": "[[72,12,6]]", "n": 72, "k": 12, "d": 6,
        "ell": 6, "m": 6,
        "A": [[3, 0], [0, 1], [0, 2]],
        "B": [[0, 3], [1, 0], [2, 0]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/bb_distance_72_12_6.json",
    },
    {
        "name": "[[90,8,10]]", "n": 90, "k": 8, "d": 10,
        "ell": 15, "m": 3,
        "A": [[9, 0], [0, 1], [0, 2]],
        "B": [[0, 0], [2, 0], [7, 0]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/bb_distance_90_8_10.json",
    },
    {
        "name": "[[108,8,10]]", "n": 108, "k": 8, "d": 10,
        "ell": 9, "m": 6,
        "A": [[3, 0], [0, 1], [0, 2]],
        "B": [[0, 3], [1, 0], [2, 0]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/bb_distance_108_8_10.json",
    },
    {
        "name": "[[144,12,12]]", "n": 144, "k": 12, "d": 12,
        "ell": 12, "m": 6,
        "A": [[3, 0], [0, 1], [0, 2]],
        "B": [[0, 3], [1, 0], [2, 0]],
        "certificate_kind": "legacy_css_exact",
        "distance_certificate":
            "results/certificates/bb_distance_144_12_12.json",
    },
    {
        "name": "EXP-057 [[170,16,10]]",
        "n": 170,
        "k": 16,
        "d": 10,
        "ell": 17,
        "m": 5,
        "A": [[0, 0], [1, 0], [9, 1]],
        "B": [[0, 0], [2, 1], [15, 1]],
        "certificate_kind": "exp057_odd_exact",
        "distance_certificate":
            "results/certificates/exp057_170_16_10_distance.json",
    },
    {
        "name": "EXP-058 [[186,10,14]]",
        "n": 186,
        "k": 10,
        "d": 14,
        "ell": 31,
        "m": 3,
        "A": [[0, 0], [1, 0], [12, 0]],
        "B": [[0, 0], [3, 1], [8, 1]],
        "certificate_kind": "exp058_odd_exact",
        "distance_certificate":
            "results/certificates/exp058_186_10_14_distance.json",
    },
    {
        "name": "EXP-060 [[210,18,8]]",
        "n": 210,
        "k": 18,
        "d": 8,
        "ell": 15,
        "m": 7,
        "A": [[0, 0], [0, 1], [3, 3]],
        "B": [[0, 0], [0, 1], [12, 3]],
        "certificate_kind": "exp060_odd_exact",
        "distance_certificate":
            "results/certificates/exp060_210_18_8_distance.json",
    },
    {
        "name": "EXP-064 [[210,24,4]]",
        "n": 210,
        "k": 24,
        "d": 4,
        "ell": 15,
        "m": 7,
        "A": [[0, 0], [0, 1], [0, 3]],
        "B": [[0, 0], [1, 0], [4, 0]],
        "certificate_kind": "exp064_odd_exact",
        "distance_certificate":
            "results/certificates/exp064_210_24_4_distance.json",
    },
    {
        "name": "EXP-064 [[210,14,12]]",
        "n": 210,
        "k": 14,
        "d": 12,
        "ell": 15,
        "m": 7,
        "A": [[0, 0], [1, 1], [4, 3]],
        "B": [[0, 0], [1, 2], [4, 6]],
        "certificate_kind": "exp064_odd_exact",
        "distance_certificate":
            "results/certificates/exp064_210_14_12_distance.json",
    },
    {
        "name": "EXP-064 [[210,10,16]]",
        "n": 210,
        "k": 10,
        "d": 16,
        "ell": 15,
        "m": 7,
        "A": [[0, 0], [1, 1], [2, 3]],
        "B": [[0, 0], [1, 6], [11, 2]],
        "certificate_kind": "exp064_odd_exact",
        "distance_certificate":
            "results/certificates/exp064_210_10_16_distance.json",
    },
    {
        "name": "EXP-067 [[234,8,18]] row 2",
        "n": 234,
        "k": 8,
        "d": 18,
        "ell": 39,
        "m": 3,
        "A": [[0, 0], [1, 0], [5, 0]],
        "B": [[0, 0], [1, 1], [23, 2]],
        "certificate_kind": "exp067_odd_exact",
        "distance_certificate":
            "results/certificates/exp067_234_8_18_bundle19_distance.json",
    },
    {
        "name": "EXP-067 [[234,8,18]] row 1",
        "n": 234,
        "k": 8,
        "d": 18,
        "ell": 39,
        "m": 3,
        "A": [[0, 0], [1, 0], [5, 0]],
        "B": [[0, 0], [2, 2], [22, 1]],
        "certificate_kind": "exp067_odd_exact",
        "distance_certificate":
            "results/certificates/exp067_234_8_18_bundle22_distance.json",
    },
    {
        "name": "EXP-037 [[180,8,16]]",
        "n": 180,
        "k": 8,
        "d": 16,
        "ell": 15,
        "m": 6,
        "A": [[2, 5], [3, 5], [4, 2]],
        "B": [[0, 0], [5, 5], [4, 1]],
        "source_catalogue_index": 28,
        "certificate_kind": "exp037_css_exact",
        "distance_certificate":
            "results/partial_runs/exp037/css_855e8bf135ba4e05.json",
    },
    {
        "name": "EXP-037 [[180,20,6]]",
        "n": 180,
        "k": 20,
        "d": 6,
        "ell": 15,
        "m": 6,
        "A": [[10, 2], [5, 3], [5, 4]],
        "B": [[0, 0], [10, 5], [0, 4]],
        "source_catalogue_index": 211,
        "certificate_kind": "exp037_css_exact",
        "distance_certificate":
            "results/partial_runs/exp037/css_15654a1cebe3e639.json",
    },
]


def domination_threshold(n: int, k: int) -> tuple[int, str]:
    """Largest locally exact distance that could dominate candidate (n,k,.).

    A candidate [[n,k,d]] is Pareto-dominated by [[n',k',d']] when
    n' <= n, k' >= k and d' >= d.  Only independent two-sided certificates are
    admitted.  Wang--Mueller's BP-OSD ``distance_upperbound`` outputs, Postema
    Monte-Carlo estimates, and un-replayed source values are all excluded.
    """
    best, who = 0, "none"
    for r in EXACT_REFERENCES:
        if r["n"] > n or r["k"] < k or r["d"] <= best:
            continue
        exact, _certificate_sha256 = _validated_exact_reference(r)
        if exact:
            best, who = r["d"], r["name"]
    for r in LITERATURE_ODD:
        if not _distance_exact_here(r):
            continue
        nn = 2 * r["ell"] * r["m"]
        if nn <= n and r["k"] >= k and r["d"] > best:
            best, who = r["d"], f"[[{nn},{r['k']},{r['d']}]] {r['source']}"
    return best, who


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_fingerprint() -> str:
    rows = []
    for rec in EXACT_REFERENCES:
        exact, certificate_sha256 = _validated_exact_reference(rec)
        if exact:
            rows.append({**rec, "certificate_sha256": certificate_sha256})
    for rec in LITERATURE_ODD:
        exact, certificate_sha256 = _validated_distance_reference(rec)
        if not exact:
            continue
        rows.append(
            {
                "ell": rec["ell"], "m": rec["m"],
                "k": rec["k"], "d": rec["d"],
                "source": rec["source"],
                "certificate_sha256": certificate_sha256,
            }
        )
    return hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

def _coprime(nmod: int) -> list[int]:
    from math import gcd
    return [u for u in range(1, nmod) if gcd(u, nmod) == 1]


def _canon_translate(terms, ell: int, m: int) -> tuple:
    """Lexicographic-least translate of a support (translations are qubit maps)."""
    return min(tuple(sorted(((a - p) % ell, (b - q) % m) for a, b in terms))
               for (p, q) in terms)


def _orbit_tables(ell: int, m: int, weight: int = 3):
    """Support ids under translation, plus the induced action of the rest of the
    BB symmetry group.

    Generators: independent monomial multiplication of A and of B (folded into
    the id by taking the canonical translate), the ring automorphisms
    x -> x^u, y -> y^v for units u, v, the transpose/swap (A,B) -> (B,A), and
    x <-> y when ell == m.  All are qubit permutations, so [[n,k,d]] is exactly
    preserved.  Returns (ids, reps, maps) with maps[t][i] the image of id i.
    """
    ids: dict[tuple, int] = {}
    reps: list[tuple] = []
    for s in normalised_supports(ell, m, weight):
        if len(s) != weight:
            continue
        c = _canon_translate(s, ell, m)
        if c not in ids:
            ids[c] = len(reps)
            reps.append(c)
    transforms = []
    for u in _coprime(ell):
        for v in _coprime(m):
            transforms.append(lambda t, u=u, v=v: [((a * u) % ell, (b * v) % m)
                                                   for a, b in t])
    if ell == m:
        transforms += [lambda t, f=f: [(b, a) for a, b in f(t)] for f in transforms]
    maps = []
    for f in transforms:
        maps.append([ids[_canon_translate(f(r), ell, m)] for r in reps])
    return ids, reps, maps


def enumerate_candidates(ell: int, m: int, k_min: int, k_max: int) -> list[dict]:
    """Symmetry-reduced weight-3 x weight-3 pairs with k in [k_min, k_max].

    k and the certified ceiling depend only on (Ann(a), Ann(b)), so both are
    computed once per GROUP pair; the symmetry reduction then runs on integer
    support ids, which is what makes an exhaustive pass affordable.
    """
    bar = bar_permutation(ell, m)
    ids, reps, maps = _orbit_tables(ell, m, 3)
    groups: dict[bytes, dict] = {}
    for c, i in ids.items():
        A = ann_basis(list(c), ell, m)
        if 2 * A.shape[0] < k_min:
            continue
        key = canon(A)
        g = groups.setdefault(key, {"ann": A, "members": []})
        g["members"].append(i)
    keys = list(groups)
    seen: dict[tuple, dict] = {}
    visited_pairs: set[tuple[int, int]] = set()
    for i, ka in enumerate(keys):
        ga = groups[ka]
        for kb in keys[i:]:
            gb = groups[kb]
            I = intersect(ga["ann"], gb["ann"])
            k = 2 * I.shape[0]
            if not k_min <= k <= k_max:
                continue
            cel = certified_ceiling(I, ell, m, bar)
            for ia in ga["members"]:
                for ib in gb["members"]:
                    raw = (ia, ib) if ia <= ib else (ib, ia)
                    if raw in visited_pairs:
                        continue
                    orbit = set()
                    for mp in maps:
                        x, y = mp[ia], mp[ib]
                        orbit.add((x, y) if x <= y else (y, x))
                    visited_pairs.update(orbit)
                    best = min(orbit)  # class key only; witness is for original pair
                    seen[best] = {
                        "A": list(reps[ia]), "B": list(reps[ib]),
                        "k_parent": k, "orbit": len(orbit),
                        "ceiling": cel.get("ceiling"),
                        "ceiling_witness_support": cel.get("witness_support"),
                        "dim_pole": cel.get("dim_pole"),
                        "bar_invariant": cel.get("bar_invariant"),
                    }
    return sorted(seen.values(), key=lambda r: (-r["k_parent"], r["A"], r["B"]))


def _screen_lattice(task: tuple) -> str:
    """One lattice = one config-bound coarse task; writes its own shard."""
    ell, m, protocol, force = task
    k_min, k_max = protocol["k_range"]
    time_limit = protocol["time_limit_s"]
    p = SCREEN_DIR / f"{ell}x{m}.json"
    if p.exists() and not force:
        old = json.loads(p.read_text())
        if old.get("schema") == "exp055-screen-v3" and old.get("protocol") == protocol:
            return f"({ell},{m}) cached"
    t0 = time.time()
    cands = enumerate_candidates(ell, m, k_min, k_max)
    recs = [_screen_one((ell, m, c, time_limit, False)) for c in cands]
    residual = [i for i, r in enumerate(recs) if r["verdict"] == "solver_required"]
    if residual:
        workers = min(protocol["residual_workers"], len(residual))
        tasks = [(ell, m, cands[i], time_limit, True) for i in residual]
        with ThreadPoolExecutor(max_workers=workers) as ex:
            solved = list(ex.map(_screen_one, tasks))
        for i, rec in zip(residual, solved):
            recs[i] = rec
    verdicts: dict[str, int] = {}
    for r in recs:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    payload = {
        "schema": "exp055-screen-v3", "utc": E53.E52.utc_now(),
        "protocol": protocol,
        "ell": ell, "m": m, "n": 2 * ell * m,
        "k_range": [k_min, k_max], "time_limit_s": time_limit,
        "candidates_after_symmetry": len(cands),
        "orbit_total": sum(c["orbit"] for c in cands),
        "verdicts": verdicts,
        "survivors": [r for r in recs if r["verdict"] == "survivor"],
        "no_reference": [r for r in recs if r["verdict"] == "no_reference"],
        "undecided": [r for r in recs if r["verdict"] == "undecided"],
        "solver_calls": sum(r.get("solver_calls", 0) for r in recs),
        "records": recs, "wall_s": round(time.time() - t0, 1),
    }
    E53.E52.atomic_write_json(p, payload)
    return (f"({ell},{m}) n={2*ell*m} cands={len(cands)} "
            f"orbits={payload['orbit_total']} solver={payload['solver_calls']} "
            f"{verdicts} {payload['wall_s']}s")
def _screen_protocol(census_sha256: str, reference_sha256: str,
                     k_min: int, k_max: int, time_limit: float) -> dict:
    return {
        "census_sha256": census_sha256,
        "reference_sha256": reference_sha256,
        "reference_validation_version": REFERENCE_VALIDATION_VERSION,
        "k_range": [k_min, k_max],
        "time_limit_s": float(time_limit),
        "witness_tries": WITNESS_TRIES,
        "witness_keep": WITNESS_KEEP,
        "cdcl_solver": SCREEN_CDCL_SOLVER,
        "cdcl_conflict_budget": SCREEN_CDCL_CONFLICT_BUDGET,
        "residual_workers": SCREEN_RESIDUAL_WORKERS,
    }


def _validate_screen_for_certification(screen_payload: dict,
                                       census_sha256: str,
                                       reference_sha256: str) -> None:
    if screen_payload.get("schema") != "exp055-odd-lattice-screen-v2":
        raise RuntimeError("unexpected or stale screen schema")
    verdict = screen_payload.get("verdict", {})
    if (not verdict.get("complete") or not verdict.get("all_referenced_decided")
            or screen_payload.get("scope", {}).get("missing")):
        raise RuntimeError("screen is incomplete/undecided; refusing survivor claims")
    protocol = screen_payload.get("protocol", {})
    if protocol.get("census_sha256") != census_sha256:
        raise RuntimeError("screen is not bound to the current census")
    if protocol.get("reference_sha256") != reference_sha256:
        raise RuntimeError("screen reference set has changed")




def _validate_census_for_screen(census: dict, k_min: int, k_max: int) -> None:
    if not census.get("scope", {}).get("complete"):
        raise RuntimeError("EXP-055 census is incomplete; refusing the screen")
    protocol = census.get("protocol", {})
    if int(protocol.get("max_weight", 0)) < 3:
        raise RuntimeError("census did not include weight-3 polynomials")
    if int(protocol.get("k_min", 10**9)) > k_min:
        raise RuntimeError("census k_min is above the requested screen k_min")
    if k_max > K_CEIL_MAX:
        raise RuntimeError("screen k_max exceeds the census ceiling range")


def screen(args: argparse.Namespace) -> int:
    """Rigorous Pareto screen: is any odd-lattice weight-3 BB code undominated?"""
    census = json.loads(OUT.read_text())
    _validate_census_for_screen(census, args.k_min, args.k_max)
    lats = ([tuple(int(v) for v in tok.split("x")) for tok in args.lattices.split(",")]
            if args.lattices else
            [(r["ell"], r["m"]) for r in census["lattices"]
             if r["frontier_weight3"] and r["n"] >= args.min_n])
    protocol = _screen_protocol(
        _file_sha256(OUT), _reference_fingerprint(),
        args.k_min, args.k_max, args.time_limit,
    )
    SCREEN_DIR.mkdir(parents=True, exist_ok=True)
    tasks = [(ell, m, protocol, args.force) for ell, m in lats]
    print(f"{len(tasks)} lattices, {args.workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for line in ex.map(_screen_lattice, tasks):
            print(line, flush=True)
    return assemble_screen(args)


def _validate_screen_shard_aggregates(shard: dict) -> None:
    """Reject cached shard summaries that disagree with their record payload."""
    records = shard.get("records", [])
    verdicts: dict[str, int] = {}
    for record in records:
        verdict = record["verdict"]
        threshold = int(record.get("threshold", -1))
        if verdict == "no_reference":
            if threshold != 0:
                raise RuntimeError(
                    "no_reference screen record has an admissible reference"
                )
        elif (
            verdict in {"survivor", "undecided"}
            or verdict.startswith("dominated")
        ):
            if threshold <= 0:
                raise RuntimeError(
                    "referenced screen verdict lacks an admissible reference"
                )
        else:
            raise RuntimeError(f"nonterminal screen verdict: {verdict}")
        verdicts[verdict] = verdicts.get(verdict, 0) + 1
    expected = {
        "candidates_after_symmetry": len(records),
        "orbit_total": sum(int(record["orbit"]) for record in records),
        "verdicts": verdicts,
        "survivors": [
            record for record in records if record["verdict"] == "survivor"
        ],
        "no_reference": [
            record for record in records if record["verdict"] == "no_reference"
        ],
        "undecided": [
            record for record in records if record["verdict"] == "undecided"
        ],
        "solver_calls": sum(
            int(record.get("solver_calls", 0)) for record in records
        ),
    }
    if any(shard.get(key) != value for key, value in expected.items()):
        raise RuntimeError("screen shard aggregate does not match its records")
def _screen_record_identity(record: dict) -> dict:
    return {
        "ell": int(record["ell"]),
        "m": int(record["m"]),
        "A": record["A"],
        "B": record["B"],
        "k_parent": int(record["k_parent"]),
        "threshold": int(record["threshold"]),
    }


def _screen_witness_valid(
    record: dict, support: list[int] | None, claimed_weight: int | None
) -> bool:
    if support is None or claimed_weight is None:
        return False
    n = int(record["n"])
    indexes = [int(index) for index in support]
    if len(indexes) != len(set(indexes)) or any(
        index < 0 or index >= n for index in indexes
    ):
        return False
    HX, HZ = E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(n, dtype=np.uint8)
    vector[np.asarray(indexes, dtype=int)] = 1
    return bool(
        int(vector.sum()) == int(claimed_weight)
        and not np.any(HX @ vector % 2)
        and rank_np(np.vstack([HZ, vector])) == rank_np(HZ) + 1
    )


def _fallback_witness_evidence_valid(record: dict) -> bool:
    binding = record.get("fallback_witness")
    if not isinstance(binding, dict):
        return False
    relative = binding.get("path")
    if not isinstance(relative, str):
        return False
    path = ROOT / relative
    if not path.is_file() or _file_sha256(path) != binding.get("sha256"):
        return False
    evidence = json.loads(path.read_text(encoding="utf-8"))
    verification = evidence.get("verification", {})
    return bool(
        evidence.get("schema") == "exp068-screen-fallback-witness-v1"
        and evidence.get("identity") == _screen_record_identity(record)
        and evidence.get("support") == record.get("witness_support")
        and int(evidence.get("weight", -1))
        == int(record.get("witness_bound", -2))
        and verification.get("commutes_with_HX") is True
        and verification.get("outside_Z_stabilizer") is True
    )


def _validate_screen_shard_records(shard: dict) -> None:
    """Rebuild every record identity, threshold, and claimed domination proof."""
    if shard.get("schema") != "exp055-screen-v3":
        raise RuntimeError("unexpected screen shard schema")
    ell, m = int(shard["ell"]), int(shard["m"])
    k_min, k_max = map(int, shard["protocol"]["k_range"])
    records = shard.get("records", [])
    def identity_key(value: dict) -> str:
        return json.dumps(
            {
                "A": [list(term) for term in value["A"]],
                "B": [list(term) for term in value["B"]],
                "k_parent": int(value["k_parent"]),
                "orbit": int(value["orbit"]),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    transport = shard.get("transport")
    if isinstance(transport, dict):
        E59 = _load("exp059_screen_shard_validator", "exp059_coprime_transport.py")
        source = tuple(map(int, transport.get("source_lattice", [])))
        target = (ell, m)
        if len(source) != 2 or transport.get("target_lattice") != [ell, m]:
            raise RuntimeError("transported shard lattice binding is stale")
        source_path = SCREEN_DIR / f"{source[0]}x{source[1]}.json"
        source_shard = json.loads(source_path.read_text(encoding="utf-8"))
        _validate_screen_shard_aggregates(source_shard)
        _validate_screen_shard_records(source_shard)
        mapping = E59.coordinate_transport(*source, *target)
        expected_records = [
            E59._transport_record(record, source, target, mapping)
            for record in source_shard["records"]
        ]
        theorem_path = ROOT / str(transport.get("theorem_certificate"))
        transport_valid = bool(
            transport.get("schema") == E59.SCHEMA
            and transport.get("valid") is True
            and transport.get("mapping_sha256")
            == E59.matrix_sha256(mapping[None, :])
            and transport.get("source_shard_sha256")
            == _file_sha256(source_path)
            and theorem_path.is_file()
            and transport.get("theorem_certificate_sha256")
            == _file_sha256(theorem_path)
        )
        if not transport_valid:
            raise RuntimeError("transported shard proof binding is stale")
    else:
        expected_records = enumerate_candidates(ell, m, k_min, k_max)
    expected_keys = {identity_key(candidate) for candidate in expected_records}
    record_keys = [identity_key(record) for record in records]
    if (
        len(records) != len(expected_records)
        or len(set(record_keys)) != len(record_keys)
        or set(record_keys) != expected_keys
    ):
        raise RuntimeError("screen shard comparison-class identity mismatch")
    for record in records:
        if (
            int(record.get("ell", -1)),
            int(record.get("m", -1)),
            int(record.get("n", -1)),
        ) != (ell, m, 2 * ell * m):
            raise RuntimeError("screen shard record lattice identity mismatch")
        threshold, source = domination_threshold(
            int(record["n"]), int(record["k_parent"])
        )
        if (record.get("threshold"), record.get("threshold_source")) != (
            threshold,
            source,
        ):
            raise RuntimeError("screen shard record threshold is stale")
        verdict = record["verdict"]
        if verdict == "no_reference":
            if threshold != 0:
                raise RuntimeError("no_reference record has an exact reference")
            continue
        if verdict in {"survivor", "undecided"}:
            if threshold <= 0:
                raise RuntimeError("referenced open record has no exact reference")
            if verdict == "survivor" and record.get("screen_decided") is not True:
                raise RuntimeError("survivor lacks a complete lower decision")
            continue
        if threshold <= 0:
            raise RuntimeError("domination record has no exact reference")
        if verdict == "dominated_by_ceiling":
            valid = bool(
                record.get("ceiling") is not None
                and int(record["ceiling"]) <= threshold
                and _screen_witness_valid(
                    record,
                    record.get("ceiling_witness_support"),
                    record.get("ceiling"),
                )
            )
        elif verdict in {
            "dominated_by_witness",
            "dominated_by_cdcl_witness",
            "dominated_by_automorphism_transport",
            "dominated_by_deep_reduction",
            "dominated_by_exact_witness",
        }:
            valid = bool(
                record.get("witness_bound") is not None
                and int(record["witness_bound"]) <= threshold
                and _screen_witness_valid(
                    record,
                    record.get("witness_support"),
                    record.get("witness_bound"),
                )
            )
            if verdict == "dominated_by_cdcl_witness":
                valid = valid and bool(
                    record.get("cdcl", {}).get("status") == "SAT"
                    or isinstance(record.get("reference_rebind_ratchet"), dict)
                )
        elif verdict == "dominated":
            valid = bool(
                record.get("witness_bound") is not None
                and int(record["witness_bound"]) <= threshold
                and _screen_witness_valid(
                    record,
                    record.get("witness_support"),
                    record.get("witness_bound"),
                )
                and _fallback_witness_evidence_valid(record)
            )
        else:
            raise RuntimeError(f"unknown screen verdict: {verdict}")
        if not valid:
            raise RuntimeError(
                "screen domination lacks a physical proof: "
                f"{_screen_record_identity(record)} verdict={verdict}"
            )


def assemble_screen(args: argparse.Namespace) -> int:
    """Assemble only shards hash-bound to this exact screen protocol."""
    census = json.loads(OUT.read_text())
    k_min = int(getattr(args, "k_min", K_MIN))
    k_max = int(getattr(args, "k_max", K_CEIL_MAX))
    _validate_census_for_screen(census, k_min, k_max)
    protocol = _screen_protocol(
        _file_sha256(OUT), _reference_fingerprint(),
        k_min, k_max, float(getattr(args, "time_limit", 120.0)),
    )
    if getattr(args, "lattices", ""):
        expected = [tuple(int(v) for v in tok.split("x"))
                    for tok in args.lattices.split(",")]
    else:
        min_n = int(getattr(args, "min_n", 18))
        expected = [(r["ell"], r["m"]) for r in census["lattices"]
                    if r["frontier_weight3"] and r["n"] >= min_n]
    shards: dict[tuple[int, int], dict] = {}
    for p in SCREEN_DIR.glob("*x*.json"):
        d = json.loads(p.read_text())
        if d.get("schema") == "exp055-screen-v3" and d.get("protocol") == protocol:
            _validate_screen_shard_aggregates(d)
            _validate_screen_shard_records(d)
            shards[(d["ell"], d["m"])] = d
    missing = [f"{ell}x{m}" for ell, m in expected if (ell, m) not in shards]
    recs = [shards[pair] for pair in expected if pair in shards]
    verdicts: dict[str, int] = {}
    survivors, no_reference, undecided = [], [], []
    for r in recs:
        for key, value in r["verdicts"].items():
            verdicts[key] = verdicts.get(key, 0) + value
        survivors.extend(r["survivors"])
        no_reference.extend(r["no_reference"])
        undecided.extend(r["undecided"])
    with_reference = sum(value for key, value in verdicts.items()
                         if key != "no_reference")
    dominated = sum(value for key, value in verdicts.items()
                    if key.startswith("dominated"))
    payload = {
        "schema": "exp055-odd-lattice-screen-v2",
        "utc": E53.E52.utc_now(),
        "protocol": protocol,
        "scope": {
            "lattices_expected": len(expected),
            "lattices_completed": len(recs),
            "missing": missing,
            "weight_A": 3,
            "weight_B": 3,
            "k_range": protocol["k_range"],
            "n_max": max((2 * ell * m for ell, m in expected), default=0),
            "reference_rule": (
                "max independently exact-certified d with n_ref <= n_candidate "
                "and k_ref >= k_candidate; estimates excluded"
            ),
        },
        "verdict": {
            "complete": not missing,
            "candidates_after_symmetry": sum(r["candidates_after_symmetry"] for r in recs),
            "orbits_represented": sum(r["orbit_total"] for r in recs),
            "solver_calls": sum(r["solver_calls"] for r in recs),
            "verdicts": verdicts,
            "with_reference": with_reference,
            "dominated": dominated,
            "survivors": len(survivors),
            "no_reference": len(no_reference),
            "undecided": len(undecided),
            "all_referenced_decided": (
                dominated + len(survivors) == with_reference
            ),
            "all_referenced_dominated": dominated == with_reference,
        },
        "survivors": survivors,
        "no_reference": no_reference,
        "undecided": undecided,
        "lattices": recs,
    }
    E53.E52.atomic_write_json(SCREEN_OUT, payload)
    print(json.dumps(payload["verdict"], indent=1))
    if missing:
        print("missing:", ", ".join(missing), file=sys.stderr)
    return 0


def certify_screen_survivors(args: argparse.Namespace) -> int:
    """Exact distance and decomposition certificates for a complete screen."""
    screen_payload = json.loads(SCREEN_OUT.read_text())
    _validate_screen_for_certification(
        screen_payload, _file_sha256(OUT), _reference_fingerprint()
    )
    from qec_research.equivalence.css import (
        stabilizer_direct_sum_certificate,
        stabilizer_incidence_components,
    )
    out = []
    for row in screen_payload["survivors"]:
        ell, m, n = row["ell"], row["m"], row["n"]
        HX, HZ = E53.bb_from_terms(ell, m, row["A"], row["B"])
        k_direct = int(n - rank_np(HX) - rank_np(HZ))
        t0 = time.time()
        dist = exact_distance_css(
            HX, HZ, time_limit_s=args.time_limit, workers=CERT_WORKERS
        )
        zeros = np.zeros_like(HX)
        H = np.vstack([np.hstack([HX, zeros]), np.hstack([zeros, HZ])])
        direct = stabilizer_direct_sum_certificate(H)
        incidence = stabilizer_incidence_components(H)
        rec = {
            **row,
            "k_from_matrices": k_direct,
            "k_matches_screen": k_direct == row["k_parent"],
            "d": dist["d"],
            "d_X": dist["d_X"],
            "d_Z": dist["d_Z"],
            "d_exact": bool(dist["d_exact"]),
            "d_X_exact": bool(dist["d_X_exact"]),
            "d_Z_exact": bool(dist["d_Z_exact"]),
            "d_X_witness": dist["d_X_witness"],
            "d_Z_witness": dist["d_Z_witness"],
            "beats_reference": bool(
                dist["d_exact"] and dist["d"] is not None
                and dist["d"] > row["threshold"]
            ),
            "incidence_component_sizes": [len(c) for c in incidence],
            "direct_sum_component_sizes": direct["component_sizes"],
            "is_direct_sum": bool(direct["is_direct_sum"]),
            "max_check_weight": int(max(HX.sum(axis=1).max(),
                                        HZ.sum(axis=1).max())),
            "max_qubit_degree": int((HX.sum(axis=0) + HZ.sum(axis=0)).max()),
            "wall_s": round(time.time() - t0, 2),
        }
        out.append(rec)
        print(json.dumps({k: rec[k] for k in
                          ("ell", "m", "n", "k_from_matrices", "d",
                           "d_exact", "beats_reference", "is_direct_sum",
                           "wall_s")}), flush=True)
    payload = {
        "schema": "exp055-odd-lattice-survivors-v2",
        "utc": E53.E52.utc_now(),
        "screen_schema": screen_payload["schema"],
        "screen_sha256": _file_sha256(SCREEN_OUT),
        "survivors": len(out),
        "all_k_match_screen": all(r["k_matches_screen"] for r in out),
        "all_exact": all(r["d_exact"] for r in out),
        "all_beat_reference": all(r["beats_reference"] for r in out),
        "all_indecomposable": all(not r["is_direct_sum"] for r in out),
        "records": out,
    }
    E53.E52.atomic_write_json(SURVIVOR_CERT, payload)
    if (not payload["all_exact"] or not payload["all_beat_reference"]
            or not payload["all_k_match_screen"]):
        raise RuntimeError("survivor certification contradicted the screen")
    return 0


def _cdcl_witness_bound(
    HX: np.ndarray, HZ: np.ndarray, block: int, threshold: int
) -> dict:
    """Find one physical Z-logical at or below threshold; UNSAT is diagnostic."""
    instance = css_side_instance(HX, HZ, "z", block_length=block)
    record = decide_weight_bounded(
        instance,
        int(threshold),
        solver_name=SCREEN_CDCL_SOLVER,
        conflict_budget=SCREEN_CDCL_CONFLICT_BUDGET,
    )
    output = {
        "status": record["status"],
        "cnf_sha256": record["cnf_sha256"],
        "encoding_version": record["encoding_version"],
        "solver": record["solver"],
    }
    if record["status"] == "SAT":
        vector = np.asarray(record["vector"], dtype=np.uint8)
        weight = int(vector.sum())
        outside = rank_np(np.vstack([HZ, vector])) == rank_np(HZ) + 1
        if (
            np.any(HX @ vector % 2)
            or not outside
            or weight != int(record["weight"])
            or weight > int(threshold)
        ):
            raise RuntimeError("CDCL screen witness failed physical verification")
        output.update(
            {
                "weight": weight,
                "witness_support": [
                    int(index) for index in np.flatnonzero(vector)
                ],
                "outside_z_stabilizer": True,
            }
        )
    return output


def _screen_one(task: tuple) -> dict:
    """Classify one candidate against the Pareto threshold.

    Order matters: the certified ceiling is free, so it is consulted first and
    rejects a candidate outright when ceiling <= threshold (then
    d <= ceiling <= threshold, so the candidate is dominated with no solver
    call at all). Survivors first reach a bounded CDCL witness query; only its
    non-SAT outcomes reach the older CP-SAT exact fallback.
    """
    ell, m, cand, time_limit, allow_solver = task
    n = 2 * ell * m
    thr, who = domination_threshold(n, cand["k_parent"])
    out = {**cand, "ell": ell, "m": m, "n": n,
           "threshold": thr, "threshold_source": who}
    if thr == 0:
        # No published code with n' <= n and k' >= k: nothing can dominate this
        # (n,k), so no solver call is warranted.  The ceiling is still recorded.
        out.update({"verdict": "no_reference", "solver_calls": 0, "wall_s": 0.0})
        return out
    ceil = cand.get("ceiling")
    if ceil is not None and ceil <= thr:
        # d <= ceiling <= threshold: dominated, proved without any solver call.
        out.update({"verdict": "dominated_by_ceiling", "solver_calls": 0,
                    "wall_s": 0.0})
        return out
    HX, HZ = E53.bb_from_terms(ell, m, cand["A"], cand["B"])
    out["k_from_matrices"] = int(n - rank_np(HX) - rank_np(HZ))
    # Certified reciprocal-pole witnesses first: microseconds, no solver.
    I = intersect(ann_basis(cand["A"], ell, m), ann_basis(cand["B"], ell, m))
    J = I[:, bar_permutation(ell, m)]
    t0 = time.time()
    rw = reduced_witness_bound(J, HZ, ell, m, target=thr)
    out["witness_bound"] = rw["bound"]
    out["witness_classes"] = rw.get("classes_tested")
    out["witness_support"] = rw.get("witness_support")
    out["witness_wall_s"] = round(time.time() - t0, 2)
    if rw["bound"] is not None and rw["bound"] <= thr:
        out.update({"verdict": "dominated_by_witness", "solver_calls": 0})
        return out
    if not allow_solver:
        out.update({"verdict": "solver_required", "solver_calls": 0})
        return out
    t0 = time.time()
    cdcl = _cdcl_witness_bound(HX, HZ, ell * m, thr)
    out["cdcl"] = cdcl
    if cdcl["status"] == "SAT":
        out.update({
            "verdict": "dominated_by_cdcl_witness",
            "solver_calls": 1,
            "witness_bound": cdcl["weight"],
            "witness_support": cdcl["witness_support"],
            "wall_s": round(time.time() - t0, 1),
        })
        return out
    res = exact_distance_css(
        HX, HZ, time_limit_s=time_limit, workers=2, upper_bound=thr
    )
    decided = bool(
        res["d_X_all_sectors_decided"] and res["d_Z_all_sectors_decided"]
    )
    out.update({
        "solver_calls": 2, "screen_decided": decided,
        "d_found": res["d"], "d_exact": bool(res["d_exact"]),
        "wall_s": round(time.time() - t0, 1),
    })
    if res["d"] is not None:
        out["verdict"] = "dominated"         # certified witness of weight <= thr
    elif decided:
        out["verdict"] = "survivor"          # proven d > threshold
    else:
        out["verdict"] = "undecided"
    return out


# --------------------------------------------------------------------------- #
# per-lattice sweep
# --------------------------------------------------------------------------- #
def sweep_lattice(ell: int, m: int, max_weight: int, k_min: int, sample: int,
                  seed: int) -> dict:
    t0 = time.time()
    rng = random.Random(hash((seed, ell, m)) & 0xFFFFFFFF)
    bar = bar_permutation(ell, m)
    supports = list(normalised_supports(ell, m, max_weight))

    groups: dict[bytes, dict] = {}
    for terms in supports:
        A = ann_basis(terms, ell, m)
        key = canon(A)
        g = groups.get(key)
        if g is None:
            groups[key] = {"ann": A, "count": 1, "rep": terms, "min_weight": len(terms)}
        else:
            g["count"] += 1
            if len(terms) < g["min_weight"]:
                g["min_weight"], g["rep"] = len(terms), terms

    keys = list(groups)
    k_hist: dict[int, int] = {}
    qualifying: list[dict] = []
    pairs_total = 0
    for i, ka in enumerate(keys):
        ga = groups[ka]
        for kb in keys[i:]:
            gb = groups[kb]
            mult = ga["count"] * gb["count"] * (1 if ka == kb else 2)
            pairs_total += mult
            I = intersect(ga["ann"], gb["ann"])
            k_par = 2 * I.shape[0]
            k_hist[k_par] = k_hist.get(k_par, 0) + mult
            if k_par < k_min or k_par > K_CEIL_MAX:
                continue
            cel = certified_ceiling(I, ell, m, bar)
            qualifying.append({
                "A": ga["rep"], "B": gb["rep"], "k_parent": k_par,
                "pairs_represented": mult, "wt_A": ga["min_weight"],
                "wt_B": gb["min_weight"], **cel,
            })


    # frontier: per k, the largest certified ceiling.  Two versions: over the
    # whole weight-<=max_weight family, and restricted to the catalogue's own
    # shape (both polynomials of weight exactly 3).
    def _frontier(pool) -> dict[int, dict]:
        best: dict[int, dict] = {}
        for q in pool:
            c = q.get("ceiling")
            if c is None:
                continue
            cur = best.get(q["k_parent"])
            if cur is None or c > cur["ceiling"]:
                best[q["k_parent"]] = {
                    "ceiling": c, "A": q["A"], "B": q["B"],
                    "bar_invariant": q["bar_invariant"],
                    "dim_pole": q.get("dim_pole"),
                    "exact_min": bool(q.get("ceiling_is_exact_min", False)),
                    "route": "reciprocal_pole",
                }
        return best

    frontier = _frontier(qualifying)
    frontier_w3 = _frontier([q for q in qualifying
                             if q["wt_A"] == 3 and q["wt_B"] == 3])

    # cross-check 1: k from the matrices, independent of the ideal theory
    checks, mismatches = 0, []
    pool = qualifying if qualifying else []
    for q in pool[: max(1, sample)]:
        HX, HZ = E53.bb_from_terms(ell, m, q["A"], q["B"])
        k_direct = 2 * ell * m - rank_np(HX) - rank_np(HZ)
        checks += 1
        if k_direct != q["k_parent"]:
            mismatches.append({"A": q["A"], "B": q["B"], "ideal": q["k_parent"],
                               "matrices": int(k_direct)})
    for _ in range(sample):
        a, b = rng.choice(supports), rng.choice(supports)
        I = intersect(ann_basis(a, ell, m), ann_basis(b, ell, m))
        HX, HZ = E53.bb_from_terms(ell, m, a, b)
        k_direct = 2 * ell * m - rank_np(HX) - rank_np(HZ)
        checks += 1
        if k_direct != 2 * I.shape[0]:
            mismatches.append({"A": a, "B": b, "ideal": 2 * I.shape[0],
                               "matrices": int(k_direct)})

    # cross-check 2: semisimplicity -- every Ann is idempotent (the J-G1 mechanism)
    idem_tested, idem_bad = 0, []
    for _ in range(min(sample, 12)):
        t = rng.choice(supports)
        A = ann_basis(t, ell, m)
        if A.shape[0] == 0:
            continue
        st = E53.stable_ideal_power(A, ell, m)
        idem_tested += 1
        if rank_np(st) != rank_np(A):
            idem_bad.append({"terms": t, "dim_ann": int(rank_np(A)),
                             "dim_stable": int(rank_np(st))})

    return {
        "schema": "exp055-census-lattice-v2",
        "protocol": {
            "max_weight": max_weight, "k_min": k_min,
            "sample": sample, "seed": seed,
        },
        "ell": ell, "m": m, "n": 2 * ell * m, "max_weight": max_weight,
        "supports": len(supports), "groups": len(keys), "pairs_total": pairs_total,
        "k_histogram": {str(k): v for k, v in sorted(k_hist.items())},
        "k_max": max(k_hist) if k_hist else 0,
        "qualifying_group_pairs": len(qualifying),
        "qualifying_pairs_represented": sum(q["pairs_represented"] for q in qualifying),
        "no_group_certificate_pairs": 0,
        "no_group_certificate_represented": 0,
        "frontier": {str(k): v for k, v in sorted(frontier.items())},
        "frontier_weight3": {str(k): v for k, v in sorted(frontier_w3.items())},
        "ceiling_max": max((v["ceiling"] for v in frontier.values()), default=None),
        "ceiling_max_weight3": max((v["ceiling"] for v in frontier_w3.values()),
                                   default=None),
        "k_crosschecks": checks, "k_mismatches": mismatches,
        "idempotence_tested": idem_tested, "idempotence_violations": idem_bad,
        "wall_s": round(time.time() - t0, 2),
    }


def _shard_path(ell: int, m: int) -> Path:
    return SHARDS / f"{ell}x{m}.json"


def _worker(task: tuple) -> str:
    ell, m, max_weight, k_min, sample, seed, force = task
    p = _shard_path(ell, m)
    identity = {"max_weight": max_weight, "k_min": k_min,
                "sample": sample, "seed": seed}
    if p.exists() and not force:
        old = json.loads(p.read_text())
        if old.get("schema") == "exp055-census-lattice-v2" and old.get("protocol") == identity:
            return f"({ell},{m}) cached"
    rec = sweep_lattice(ell, m, max_weight, k_min, sample, seed)
    E53.E52.atomic_write_json(p, rec)
    return (f"({ell},{m}) n={rec['n']} groups={rec['groups']} "
            f"kmax={rec['k_max']} qual={rec['qualifying_group_pairs']} "
            f"ceil={rec['ceiling_max']} bad={len(rec['k_mismatches'])} "
            f"{rec['wall_s']}s")


def run(args: argparse.Namespace) -> int:
    SHARDS.mkdir(parents=True, exist_ok=True)
    lats = odd_lattices(args.max_dim)
    if args.lattices:
        lats = [tuple(int(v) for v in tok.split("x")) for tok in args.lattices.split(",")]
    tasks = [(ell, m, args.max_weight, args.k_min, args.sample, args.seed, args.force)
             for ell, m in lats]
    workers = args.workers or min(len(tasks), max(1, (os.cpu_count() or 4) - 2))
    print(f"{len(tasks)} odd lattices, {workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for line in ex.map(_worker, tasks):
            print(line, flush=True)
    return assemble(args)


def assemble(args: argparse.Namespace) -> int:
    identity = {
        "max_weight": int(getattr(args, "max_weight", 3)),
        "k_min": int(getattr(args, "k_min", K_MIN)),
        "sample": int(getattr(args, "sample", 12)),
        "seed": int(getattr(args, "seed", 55)),
    }
    if getattr(args, "lattices", ""):
        expected = [tuple(int(v) for v in tok.split("x"))
                    for tok in args.lattices.split(",")]
    else:
        expected = odd_lattices(int(getattr(args, "max_dim", MAX_DIM)))
    by_lattice: dict[tuple[int, int], dict] = {}
    for p in SHARDS.glob("*x*.json"):
        r = json.loads(p.read_text())
        if (r.get("schema") == "exp055-census-lattice-v2"
                and r.get("protocol") == identity):
            by_lattice[(r["ell"], r["m"])] = r
    missing = [f"{ell}x{m}" for ell, m in expected
               if (ell, m) not in by_lattice]
    recs = [by_lattice[pair] for pair in expected if pair in by_lattice]
    recs.sort(key=lambda r: (r["n"], r["ell"]))
    k_hist: dict[int, int] = {}
    for r in recs:
        for k, v in r["k_histogram"].items():
            k_hist[int(k)] = k_hist.get(int(k), 0) + v

    # global frontier: best certified ceiling per (n, k), both families
    def _global(field: str) -> dict[str, dict]:
        out: dict[str, dict] = {}
        for r in recs:
            for k, v in r[field].items():
                key = f"{r['n']}_{k}"
                cur = out.get(key)
                if cur is None or v["ceiling"] > cur["ceiling"]:
                    out[key] = {"n": r["n"], "k": int(k), "ell": r["ell"],
                                "m": r["m"], **v}
        return out

    frontier = _global("frontier")
    frontier_w3 = _global("frontier_weight3")
    published = {}
    for name, spec in BRAVYI_BB.items():
        if spec.ell % 2 and spec.m % 2:
            published[name] = {"ell": spec.ell, "m": spec.m, "n": 2 * spec.ell * spec.m}
    payload = {
        "schema": SCHEMA, "utc": E53.E52.utc_now(),
        "protocol": identity,
        "scope": {
            "lattices_expected": len(expected),
            "lattices_completed": len(recs),
            "missing": missing,
            "complete": not missing,
        },
        "max_weight": identity["max_weight"], "k_min": identity["k_min"],
        "lattices": recs,
        "verdict": {
            "lattices_swept": len(recs),
            "complete": not missing,
            "pairs_total": sum(r["pairs_total"] for r in recs),
            "k_histogram_global": {str(k): v for k, v in sorted(k_hist.items())},
            "k_max_global": max(k_hist) if k_hist else 0,
            "k_mismatches": sum(len(r["k_mismatches"]) for r in recs),
            "idempotence_violations": sum(len(r["idempotence_violations"]) for r in recs),
            "idempotence_tested": sum(r["idempotence_tested"] for r in recs),
            "lattices_with_k_ge_kmin": [f"{r['ell']}x{r['m']}" for r in recs
                                        if r["qualifying_group_pairs"] > 0],
            "no_group_certificate_pairs": sum(r["no_group_certificate_pairs"]
                                              for r in recs),
            "frontier": frontier,
            "frontier_weight3": frontier_w3,
            "ceiling_max_global": max((v["ceiling"] for v in frontier.values()),
                                      default=None),
            "ceiling_max_weight3": max((v["ceiling"] for v in frontier_w3.values()),
                                       default=None),
            "published_odd_instances": published,
        },
    }
    E53.E52.atomic_write_json(OUT, payload)
    v = payload["verdict"]
    print(json.dumps({k: v[k] for k in
                      ("lattices_swept", "pairs_total", "k_max_global", "k_mismatches",
                       "idempotence_violations", "idempotence_tested")}, indent=1))
    print("weight-3 frontier:", json.dumps(v["frontier_weight3"], indent=1)[:1500])
    return 0


def certify(args: argparse.Namespace) -> int:
    """Exact CP-SAT distance for the strongest frontier candidates."""
    payload = json.loads(OUT.read_text())
    rows = sorted(payload["verdict"]["frontier_weight3"].values(),
                  key=lambda r: (-r["ceiling"], r["n"]))
    rows = [r for r in rows if r["n"] <= args.max_n][: args.top]
    out = []
    for r in rows:
        ell, m = r["ell"], r["m"]
        HX, HZ = E53.bb_from_terms(ell, m, r["A"], r["B"])
        k_direct = 2 * ell * m - rank_np(HX) - rank_np(HZ)
        t0 = time.time()
        res = exact_distance_css(HX, HZ, time_limit_s=args.time_limit,
                                 workers=CERT_WORKERS)
        rec = {
            "ell": ell, "m": m, "n": r["n"], "k_frontier": r["k"],
            "k_from_matrices": int(k_direct), "A": r["A"], "B": r["B"],
            "certified_ceiling": r["ceiling"], "bar_invariant": r["bar_invariant"],
            "d": res["d"], "d_X": res["d_X"], "d_Z": res["d_Z"],
            "d_exact": bool(res["d_exact"]),
            "d_X_exact": bool(res["d_X_exact"]), "d_Z_exact": bool(res["d_Z_exact"]),
            "ceiling_respected": (res["d"] is None or res["d"] <= r["ceiling"]),
            "wall_s": round(time.time() - t0, 1),
        }
        out.append(rec)
        print(json.dumps({k: rec[k] for k in
                          ("ell", "m", "n", "k_frontier", "certified_ceiling", "d",
                           "d_exact", "ceiling_respected", "wall_s")}), flush=True)
    payload["certification"] = {
        "utc": E53.E52.utc_now(), "time_limit_s": args.time_limit,
        "workers": CERT_WORKERS, "records": out,
        "all_ceilings_respected": all(r["ceiling_respected"] for r in out),
        "best_certified": max(
            ((r["d"], r["n"], r["k_frontier"]) for r in out
             if r["d"] is not None and r["d_exact"]), default=None),
    }
    E53.E52.atomic_write_json(OUT, payload)
    print(json.dumps(payload["certification"]["best_certified"]))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(required=True)
    r = sub.add_parser("run")
    r.add_argument("--max-dim", type=int, default=MAX_DIM)
    r.add_argument("--lattices", default="")
    r.add_argument("--max-weight", type=int, default=3)
    r.add_argument("--k-min", type=int, default=K_MIN)
    r.add_argument("--sample", type=int, default=12)
    r.add_argument("--seed", type=int, default=55)
    r.add_argument("--workers", type=int, default=0)
    r.add_argument("--force", action="store_true")
    r.set_defaults(fn=run)
    a = sub.add_parser("assemble")
    a.add_argument("--max-dim", type=int, default=MAX_DIM)
    a.add_argument("--lattices", default="")
    a.add_argument("--max-weight", type=int, default=3)
    a.add_argument("--k-min", type=int, default=K_MIN)
    a.add_argument("--sample", type=int, default=12)
    a.add_argument("--seed", type=int, default=55)
    a.set_defaults(fn=assemble)
    lit = sub.add_parser("literature")
    lit.set_defaults(fn=validate_literature)
    s_ = sub.add_parser("screen")
    s_.add_argument("--lattices", default="")
    s_.add_argument("--min-n", type=int, default=90)
    s_.add_argument("--k-min", type=int, default=K_MIN)
    s_.add_argument("--k-max", type=int, default=K_CEIL_MAX)
    s_.add_argument("--time-limit", type=float, default=120.0)
    s_.add_argument("--workers", type=int, default=12)
    s_.add_argument("--force", action="store_true")
    sa = sub.add_parser("screen-assemble")
    sa.add_argument("--lattices", default="")
    sa.add_argument("--min-n", type=int, default=18)
    sa.add_argument("--k-min", type=int, default=K_MIN)
    sa.add_argument("--k-max", type=int, default=K_CEIL_MAX)
    sa.add_argument("--time-limit", type=float, default=120.0)
    sa.set_defaults(fn=assemble_screen)
    sc = sub.add_parser("screen-certify")
    sc.add_argument("--time-limit", type=float, default=CERT_TIME_LIMIT_S)
    sc.set_defaults(fn=certify_screen_survivors)
    s_.set_defaults(fn=screen)
    c = sub.add_parser("certify")
    c.add_argument("--top", type=int, default=6)
    c.add_argument("--max-n", type=int, default=200)
    c.add_argument("--time-limit", type=float, default=CERT_TIME_LIMIT_S)
    c.set_defaults(fn=certify)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
