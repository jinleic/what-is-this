"""EXP-055: exhaustive weight-<=3 sweep of CSS BB codes on ODD x ODD lattices.

WHY THIS REGION.  Theorem J-G1 (EXP-053/054) proves that when ell and m are both
odd the group algebra R = GF(2)[x,y]/(x^ell-1, y^m-1) is semisimple, so every
ideal is idempotent and the single-row X-collapse channel is structurally EMPTY:
a PBB perturbation there cannot silently demote a parent stabilizer to a logical.
The published catalogue never searched this region (all 202 parents have
m in {3,6}), and the only odd x odd BB instance in the literature is [[90,8,10]]
on (15,3).  So theory hands us an unexplored, provably collapse-free region.

WHAT IS DECIDED HERE.

(1) Rate law (two independent routes, all lattices).  R semisimple means
    R = prod_chi F_chi over the Frobenius orbits of characters, multiplication by
    a is diagonal, so

        Ann(a) = sum over {chi : a(chi) = 0} F_chi,
        k_P    = 2 dim (Ann(a) cap Ann(b)) = 2 #{chi : a(chi) = b(chi) = 0},

    i.e. k is exactly twice the number of COMMON ROOTS.  Route 1 computes
    dim(Ann(a) cap Ann(b)) by GF(2) rank; route 2 computes
    k = n - rank(H_X) - rank(H_Z) from the matrices themselves.  They must agree.

(2) A solver-free CERTIFIED distance ceiling.  For u, v in I := Ann(a) cap Ann(b)
    we have H_X (u,v)^T = a u + b v = 0, so I x I <= ker H_X, and dim(I x I) =
    2 dim I = k_P = dim(ker H_X / S_Z).  The trivial part is
    (I x I) cap S_Z, whose dimension is |Z \\ Z^{-1}| in character terms; hence

        if I is bar-invariant (I = conj(I))  ==>  (I x I) cap S_Z = 0
        ==>  EVERY nonzero element of I x I is a nontrivial logical
        ==>  d(P) <= min{ wt(u) : 0 != u in I }  =: b0,

    a rigorous upper bound obtained by enumerating 2^dim I - 1 vectors -- no
    solver, microseconds.  When I is not bar-invariant the trivial subspace is
    computed explicitly and the smallest NONTRIVIAL weight is used instead, so
    the ceiling is rigorous either way.  b0 depends only on I, hence only on the
    pair of annihilators: it is a GROUP invariant, which is what makes an
    exhaustive sweep affordable.

(3) The frontier.  For every (n, k) reached, the maximum of the certified
    ceiling over ALL weight-<=3 pairs is an upper bound on what this family can
    achieve; the best pairs are then certified exactly with CP-SAT
    (exact_distance_css) and compared against the published BB instances.

Grouping polynomials by canon(Ann) collapses the O(P^2) pair sweep to O(G^2)
exactly as in EXP-054, and is exact for k and for b0 because both depend on the
pair (Ann(a), Ann(b)) alone.  Translation normalisation (supports containing the
monomial 1) is exhaustive up to independent translations of A and B, which are
qubit permutations.

Artifacts: results/processed/exp055_odd_lattice_sweep.json
           results/partial_runs/exp055/<ell>x<m>.json  (per-lattice shards)
Run: python experiments/exp055_odd_lattice_sweep.py run [--max-dim 180]
     python experiments/exp055_odd_lattice_sweep.py assemble
     python experiments/exp055_odd_lattice_sweep.py certify [--top 6]
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor
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
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

SCHEMA = "exp055-odd-lattice-sweep-v1"
OUT = ROOT / "results" / "processed" / "exp055_odd_lattice_sweep.json"
SHARDS = ROOT / "results" / "partial_runs" / "exp055"
SCREEN_DIR = ROOT / "results" / "partial_runs" / "exp055_screen"
LIT_OUT = ROOT / "results" / "processed" / "exp055_literature_validation.json"

MAX_DIM = 180            # ell*m <= MAX_DIM, i.e. n = 2*ell*m <= 360
K_MIN = 8                # the interesting regime; [[90,8,10]] has k = 8
K_CEIL_MAX = 24          # above this k the codes are degenerate; skip the ceiling
DIM_I_ENUM_CAP = 16      # 2^16 vectors per enumeration keeps memory bounded
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
    """Rigorous solver-free upper bound on d(P) from I = Ann(a) cap Ann(b).

    I x I <= ker H_X always (H_X (u,v)^T = a u + b v = 0).  The trivial part is
    (I x I) cap S_Z, which at a character chi is nonzero only when chi is in
    Z := Z(a) cap Z(b) while chi^{-1} is NOT -- because (S_Z)_chi is spanned by
    (conj(b)(chi), conj(a)(chi)) = (b(chi^{-1}), a(chi^{-1})).  Restricting to

        I_0 := I cap conj(I)     (characters Z cap Z^{-1}, inversion-closed)

    therefore kills the trivial part outright: (I_0 x I_0) cap S_Z = 0, so EVERY
    nonzero u in I_0 makes (u, 0) a nontrivial logical operator and

        d(P) <= min{ wt(u) : 0 != u in I_0 }.

    No bar-invariance hypothesis is needed, and the bound depends only on the
    pair of annihilators, so it is a group invariant.  For dim I_0 <= the
    enumeration cap the minimum is exact; above it we minimise over the basis
    rows and random combinations, which is still a valid upper bound on the
    minimum weight and hence on d.
    """
    d = I.shape[0]
    out: dict = {"dim_I": d, "k_parent": 2 * d}
    if d == 0:
        out.update({"dim_I0": 0, "ceiling": None, "bar_invariant": None})
        return out
    Ib = I[:, bar]
    out["bar_invariant"] = bool(rank_np(np.vstack([I, Ib])) == d)
    I0 = I if out["bar_invariant"] else intersect(I, Ib)
    d0 = I0.shape[0]
    out["dim_I0"] = d0
    if d0 == 0:
        # Z cap Z^{-1} empty: no group-level certificate exists, flagged for
        # per-member treatment.  Counted in the verdict, never silently dropped.
        out["ceiling"] = None
        return out
    if d0 <= DIM_I_ENUM_CAP:
        out["ceiling"] = int(span_vectors(I0).sum(axis=1).min())
        out["ceiling_is_exact_min"] = True
    else:
        rng = np.random.default_rng(0xC0FFEE ^ (d0 * 7919))
        combos = rng.integers(0, 2, size=(4096, d0), dtype=np.uint8)
        cand = np.vstack([I0, (combos @ I0) % 2])
        w = cand.sum(axis=1)
        w = w[w > 0]
        out["ceiling"] = int(w.min())
        out["ceiling_is_exact_min"] = False
    return out


def member_ceiling(A_terms, B_terms, ell: int, m: int) -> dict:
    """Ceiling for one concrete pair when I cap conj(I) = 0.

    Falls back to the smallest weight of a NONTRIVIAL (u,0), u in I, tested
    directly against rowspace(H_Z).
    """
    I = intersect(ann_basis(A_terms, ell, m), ann_basis(B_terms, ell, m))
    d = I.shape[0]
    if d == 0 or d > DIM_I_ENUM_CAP:
        return {"ceiling": None, "reason": "dim out of range"}
    HX, HZ = E53.bb_from_terms(ell, m, A_terms, B_terms)
    Rz, _ = rref_np(HZ)
    Rz = Rz[: rank_np(Rz)]
    r0 = Rz.shape[0]
    n = 2 * ell * m
    vecs = span_vectors(I)
    for idx in np.argsort(vecs.sum(axis=1), kind="stable"):
        u = vecs[idx]
        full = np.zeros(n, np.uint8)
        full[: ell * m] = u
        if rank_np(np.vstack([Rz, full[None, :]])) > r0:
            return {"ceiling": int(u.sum()), "nontrivial_found": True}
    return {"ceiling": None, "nontrivial_found": False}


# --------------------------------------------------------------------------- #
# literature validation battery
# --------------------------------------------------------------------------- #
# Every published odd x odd BB instance we could source with an EXACT distance.
# "pi" rows are the coprime construction 1 + x^a y^a + x^b y^b (pi = xy).
# Recorded 2026-08-21 from the primary texts named in `source`.
LITERATURE_ODD: list[dict] = [
    {"ell": 15, "m": 3, "k": 8, "d": 10, "A": [(9, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (2, 0), (7, 0)], "source": "2308.07915 Table 3"},
    {"ell": 3, "m": 3, "k": 4, "d": 4, "A": [(0, 0), (1, 0), (0, 1)],
     "B": [(0, 0), (2, 0), (0, 2)], "source": "2408.10001v4 App.B Table 3"},
    {"ell": 3, "m": 3, "k": 4, "d": 2, "A": [(0, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (1, 0), (0, 1)], "source": "2502.17052v4 Table 2"},
    {"ell": 3, "m": 3, "k": 8, "d": 2, "A": [(0, 0), (0, 1), (0, 2)],
     "B": [(0, 0), (1, 0), (2, 0)], "source": "2502.17052v4 Table 2"},
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
     "B": [(0, 4), (1, 0), (3, 0)], "source": "2502.17052v4 Table 2"},
    {"ell": 7, "m": 9, "k": 12, "d": 10, "pi_A": [0, 1, 58], "pi_B": [3, 16, 44],
     "source": "2408.10001v4 Table 2"},
    {"ell": 7, "m": 9, "k": 6, "d": 14, "pi_A": [0, 4, 19], "pi_B": [0, 6, 16],
     "source": "2408.10001v4 App.C Table 4"},
    {"ell": 3, "m": 21, "k": 8, "d": 10, "A": [(0, 0), (0, 2), (0, 10)],
     "B": [(0, 3), (1, 0), (2, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 5, "m": 15, "k": 16, "d": 8, "A": [(0, 0), (0, 6), (0, 8)],
     "B": [(0, 5), (1, 0), (4, 0)], "source": "2408.10001v4 Table 1"},
    {"ell": 3, "m": 27, "k": 8, "d": 14, "A": [(0, 0), (0, 10), (0, 14)],
     "B": [(0, 12), (1, 0), (2, 0)], "source": "2408.10001v4 Table 1"},
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
]


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

    (1) k by three routes where available: our k = 2 dim(Ann(a) cap Ann(b)), the
        matrix route k = n - rank H_X - rank H_Z, and (coprime lattices only) the
        published gcd formula.  All must equal the published k.
    (2) ceiling: the certified ceiling must be >= the published exact d.  When
        I cap conj(I) = 0 the group-level certificate does not exist and the
        per-member route is used instead.  One violation falsifies the theorem.
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
        k_ideal = 2 * I.shape[0]
        k_gcd = cyclic_k(A, B, ell, m)
        cel = certified_ceiling(I, ell, m, bar)
        route = "ideal"
        if cel.get("ceiling") is None and cel.get("dim_I0") == 0:
            cel.update(member_ceiling(A, B, ell, m))
            route = "member"
        ceiling = cel.get("ceiling")
        routes_agree = k_ideal == k_mat and (k_gcd is None or k_gcd == k_ideal)
        out = {
            "source": rec["source"], "ell": ell, "m": m, "n": n,
            "k_published": rec["k"], "k_ideal": k_ideal, "k_matrices": k_mat,
            "k_gcd_formula": k_gcd, "our_routes_agree": bool(routes_agree),
            "k_agrees_with_paper": bool(routes_agree and k_ideal == rec["k"]),
            "d_published": rec["d"], "ceiling": ceiling, "ceiling_route": route,
            "dim_I0": cel.get("dim_I0"), "bar_invariant": cel.get("bar_invariant"),
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
    c_bad = [r for r in core if r["ceiling_respects_d"] is False]
    slacks = [r["slack"] for r in core if r["slack"] is not None]
    payload = {
        "schema": "exp055-literature-v2", "utc": E53.E52.utc_now(),
        "instances": len(recs), "records": recs,
        "reproduced_instances": len(core),
        "k_all_agree_on_reproduced": all(r["k_agrees_with_paper"] for r in core),
        "our_routes_always_agree": all(r["our_routes_agree"] for r in recs),
        "unreproducible_rows": unreproducible,
        "ceiling_violations": c_bad,
        "ceiling_never_violated": not c_bad,
        "ceiling_certified_count": len(slacks),
        "slack_min": min(slacks) if slacks else None,
        "slack_median": sorted(slacks)[len(slacks) // 2] if slacks else None,
        "slack_max": max(slacks) if slacks else None,
        "no_certificate": [r["source"] for r in core if r["ceiling"] is None],
    }
    E53.E52.atomic_write_json(LIT_OUT, payload)
    print(json.dumps({k: payload[k] for k in
                      ("instances", "reproduced_instances",
                       "k_all_agree_on_reproduced", "our_routes_always_agree",
                       "ceiling_never_violated", "ceiling_certified_count",
                       "slack_min", "slack_median", "slack_max",
                       "no_certificate")}, indent=1))
    return 0

# --------------------------------------------------------------------------- #
# candidate enumeration, symmetry reduction, domination screen
# --------------------------------------------------------------------------- #
# Published BB instances with an exact, primary-source distance.  [[360,12,<=24]]
# is deliberately absent: its 24 is an upper bound, and an upper bound cannot
# dominate anything.  A larger reference set can only strengthen a negative
# result, never weaken it.
PUBLISHED_EXACT = [
    {"name": "[[72,12,6]]", "n": 72, "k": 12, "d": 6},
    {"name": "[[90,8,10]]", "n": 90, "k": 8, "d": 10},
    {"name": "[[108,8,10]]", "n": 108, "k": 8, "d": 10},
    {"name": "[[144,12,12]]", "n": 144, "k": 12, "d": 12},
    {"name": "[[288,12,18]]", "n": 288, "k": 12, "d": 18},
    {"name": "[[784,24,24]]", "n": 784, "k": 24, "d": 24},
]


def domination_threshold(n: int, k: int) -> tuple[int, str]:
    """Largest exact distance among known codes that would dominate (n,k,.).

    A candidate [[n,k,d]] is Pareto-dominated by a known [[n',k',d']] when
    n' <= n, k' >= k and d' >= d, so the candidate is interesting only if
    d > threshold.  The reference set is the published BB table PLUS every
    odd x odd instance whose (k, d) we independently reproduced in
    `validate_literature` -- the codes a new odd-lattice code actually has to
    beat.  Monte-Carlo distance estimates are excluded (an estimate cannot
    dominate), and so are the two rows our routes could not reproduce.
    """
    best, who = 0, "none"
    for r in PUBLISHED_EXACT:
        if r["n"] <= n and r["k"] >= k and r["d"] > best:
            best, who = r["d"], r["name"]
    for r in LITERATURE_ODD:
        nn = 2 * r["ell"] * r["m"]
        if r.get("unreproduced"):
            continue
        if nn <= n and r["k"] >= k and r["d"] > best:
            best, who = r["d"], f"[[{nn},{r['k']},{r['d']}]] {r['source']}"
    return best, who

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
                    best = None
                    for mp in maps:
                        x, y = mp[ia], mp[ib]
                        cand = (x, y) if (x, y) <= (y, x) else (y, x)
                        if best is None or cand < best:
                            best = cand
                    row = seen.get(best)
                    if row is not None:
                        row["orbit"] += 1
                        continue
                    seen[best] = {
                        "A": list(reps[ia]), "B": list(reps[ib]),
                        "k_parent": k, "orbit": 1, "ceiling": cel.get("ceiling"),
                        "dim_I0": cel.get("dim_I0"),
                        "bar_invariant": cel.get("bar_invariant"),
                    }
    return sorted(seen.values(), key=lambda r: (-r["k_parent"], r["A"], r["B"]))


def screen(args: argparse.Namespace) -> int:
    """Rigorous Pareto screen: is any odd-lattice weight-3 BB code undominated?"""
    lats = ([tuple(int(v) for v in tok.split("x")) for tok in args.lattices.split(",")]
            if args.lattices else
            [(r["ell"], r["m"]) for r in json.loads(OUT.read_text())["lattices"]
             if r["frontier_weight3"] and 2 * r["ell"] * r["m"] >= args.min_n])
    SCREEN_DIR.mkdir(parents=True, exist_ok=True)
    for (ell, m) in lats:
        p = SCREEN_DIR / f"{ell}x{m}.json"
        if p.exists() and not args.force:
            print(f"({ell},{m}) cached", flush=True)
            continue
        t0 = time.time()
        cands = enumerate_candidates(ell, m, args.k_min, args.k_max)
        tasks = [(ell, m, c, args.time_limit) for c in cands]
        recs = []
        if tasks:
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                recs = list(ex.map(_screen_one, tasks))
        verdicts: dict[str, int] = {}
        for r in recs:
            verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
        payload = {
            "schema": "exp055-screen-v1", "utc": E53.E52.utc_now(),
            "ell": ell, "m": m, "n": 2 * ell * m,
            "k_range": [args.k_min, args.k_max], "time_limit_s": args.time_limit,
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
        print(f"({ell},{m}) n={2*ell*m} cands={len(cands)} "
              f"orbits={payload['orbit_total']} solver_calls={payload['solver_calls']} "
              f"{verdicts} {payload['wall_s']}s", flush=True)
    return 0


def _screen_one(task: tuple) -> dict:
    """Classify one candidate against the Pareto threshold.

    Order matters: the certified ceiling is free, so it is consulted first and
    rejects a candidate outright when ceiling <= threshold (then
    d <= ceiling <= threshold, so the candidate is dominated with no solver
    call at all).  Only the survivors of that test reach CP-SAT.
    """
    ell, m, cand, time_limit = task
    n = 2 * ell * m
    thr, who = domination_threshold(n, cand["k_parent"])
    out = {**cand, "ell": ell, "m": m, "n": n,
           "threshold": thr, "threshold_source": who}
    ceil = cand.get("ceiling")
    if thr > 0 and ceil is not None and ceil <= thr:
        out.update({"verdict": "dominated_by_ceiling", "solver_calls": 0,
                    "wall_s": 0.0})
        return out
    HX, HZ = E53.bb_from_terms(ell, m, cand["A"], cand["B"])
    out["k_from_matrices"] = int(n - rank_np(HX) - rank_np(HZ))
    t0 = time.time()
    res = exact_distance_css(HX, HZ, time_limit_s=time_limit, workers=2,
                             upper_bound=(thr if thr > 0 else None))
    decided = bool(res["d_X_all_sectors_decided"] and res["d_Z_all_sectors_decided"])
    out.update({
        "solver_calls": 1, "screen_decided": decided,
        "d_found": res["d"], "d_exact": bool(res["d_exact"]),
        "wall_s": round(time.time() - t0, 1),
    })
    if thr == 0:
        out["verdict"] = "no_reference"      # nothing published dominates this (n,k)
    elif res["d"] is not None:
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

    # Pairs with I cap conj(I) = 0 carry no group-level certificate; resolve the
    # representative of each directly so nothing is silently dropped.
    no_cert = [q for q in qualifying if q.get("ceiling") is None
               and q.get("dim_I0") == 0]
    for q in no_cert[:64]:
        q.update(member_ceiling(q["A"], q["B"], ell, m))
        q["ceiling_route"] = "member"

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
                    "bar_invariant": q["bar_invariant"], "dim_I0": q.get("dim_I0"),
                    "exact_min": bool(q.get("ceiling_is_exact_min", False)),
                    "route": q.get("ceiling_route", "ideal"),
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
        "ell": ell, "m": m, "n": 2 * ell * m, "max_weight": max_weight,
        "supports": len(supports), "groups": len(keys), "pairs_total": pairs_total,
        "k_histogram": {str(k): v for k, v in sorted(k_hist.items())},
        "k_max": max(k_hist) if k_hist else 0,
        "qualifying_group_pairs": len(qualifying),
        "qualifying_pairs_represented": sum(q["pairs_represented"] for q in qualifying),
        "no_group_certificate_pairs": len(no_cert),
        "no_group_certificate_represented": sum(q["pairs_represented"] for q in no_cert),
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
    if p.exists() and not force:
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
    recs = [json.loads(p.read_text()) for p in sorted(SHARDS.glob("*x*.json"))]
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
        "max_weight": args.max_weight, "k_min": args.k_min,
        "lattices": recs,
        "verdict": {
            "lattices_swept": len(recs),
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
    a.add_argument("--max-weight", type=int, default=3)
    a.add_argument("--k-min", type=int, default=K_MIN)
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
