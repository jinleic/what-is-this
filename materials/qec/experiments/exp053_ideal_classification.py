"""EXP-053: the demote law is an invariant of the ideal I = Ann_R(A,B).

Theorem J-G (proved in notes/theorem_jg_ideal_invariant.md, machine-checked here).
For a CSS bivariate-bicycle parent P with H_X = [A B] over
R = GF(2)[x,y]/(x^ell-1, y^m-1), put

    I  := L_pre = {lambda in R : lambda A = lambda B = 0} = Ann_R(A,B)   (an ideal)
    M  := ker H_X / S_Z                                                  (dim = k_P)
    S  := {y in M : y in I y}    (the demote-immune fixed set, Theorem J-E'')

Then, writing I^infty for the stable ideal power (I^r = I^{r+1}):

    (1)  k_P    = 2 dim_GF(2) I            [Theorem I, already in the ledger]
    (2)  dim S  = 2 dim_GF(2) I^infty      [NEW: S is an ideal invariant]
    (3)  demote-full  <=>  I nilpotent  (I^infty = 0)
         immune       <=>  I idempotent (I^2 = I, i.e. I = e R)
         mixed        <=>  0 != I^infty != I

The module-level F-chain of EXP-052 is therefore never needed: the trichotomy
is decided by a chain of ideal powers inside R (dim <= ell*m = n/2), not by a
chain of submodules of M.  This is a strictly cheaper and structurally sharper
route, and it makes two things immediate:

  * COROLLARY (odd lattices).  If ell and m are both odd then R is semisimple
    (x^ell-1 and y^m-1 are squarefree over GF(2)), so every ideal is generated
    by an idempotent: I^infty = I always.  Hence NO parent on an odd x odd
    lattice admits a single-row syzygy demotion: every such parent with k_P > 0
    is demote-immune.  The X-distance collapse mechanism of EXP-046/047 needs
    an even lattice dimension.  (The catalogue is entirely even in m.)

  * EXISTENCE (open problem J-A resolved).  Mixed parents exist.  The smallest
    one built here lives on (ell,m) = (2,3): A = (1+x)(y+y^2), B = A*y, giving
    I = R_1 x (u) with R_1 = GF(2)[u]/(u^2) the y=1 factor and (u) the radical of
    the GF(4)[u]/(u^2) factor: dim I = 4, dim I^2 = 2 = dim I^infty, so
    k_P = 8, dim S = 4, demote fraction = 1 - 15/255.  Verified three ways:
    ideal chain, EXP-052's module F-chain, and brute-force enumeration of all
    255 nonzero quotient classes through EXP-047's exact rank test.

Artifacts: results/processed/exp053_ideal_classification.json
Run:  python experiments/exp053_ideal_classification.py run
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
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


E27 = _load("exp027", "exp027_delta_audit.py")
E39 = _load("exp039", "exp039_nogo_module.py")
E47 = _load("exp047", "exp047_exact_demotion_decision.py")
E52 = _load("exp052", "exp052_ideal_power_trichotomy.py")

from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

SCHEMA = "exp053-ideal-classification-v1"
OUT = ROOT / "results" / "processed" / "exp053_ideal_classification.json"
EXP52 = ROOT / "results" / "processed" / "exp052_ideal_power_trichotomy.json"
MAX_POWERS = 200


# ---------------------------------------------------------------- ring algebra

def poly_mult(u: np.ndarray, v: np.ndarray, ell: int, m: int) -> np.ndarray:
    """Product in R = GF(2)[x,y]/(x^ell-1, y^m-1); vectors are flat (ell*m,)
    with index x*m + y, matching EXP-047's shift convention."""
    U = u.reshape(ell, m)
    V = v.reshape(ell, m)
    out = np.zeros((ell, m), np.uint8)
    for a, b in np.argwhere(U):
        out ^= np.roll(np.roll(V, int(a), axis=0), int(b), axis=1)
    return out.reshape(-1)


def ideal_product(P: np.ndarray, Q: np.ndarray, ell: int, m: int) -> np.ndarray:
    """GF(2) row basis of the ideal product (span P)*(span Q).  Bilinearity of
    the ring multiplication makes pairwise products of GF(2) bases a spanning
    set; no R-module generators are needed."""
    if P.size == 0 or Q.size == 0:
        return np.zeros((0, ell * m), np.uint8)
    rows = [poly_mult(p, q, ell, m) for p in P for q in Q]
    basis, _ = rref_np(np.array(rows, np.uint8))
    return basis[: rank_np(basis)]


def ideal_power_chain(I: np.ndarray, ell: int, m: int) -> list[int]:
    """dims of I, I^2, I^3, ... until the power stabilizes (last two equal)."""
    cur = I
    dims = [rank_np(cur)]
    if dims[0] == 0:
        return dims + [0]
    while True:
        nxt = ideal_product(I, cur, ell, m)
        d = rank_np(nxt)
        if d > dims[-1]:
            raise RuntimeError(f"ideal powers not monotone: {dims} -> {d}")
        dims.append(d)
        if d == dims[-2]:
            return dims
        cur = nxt
        if len(dims) > MAX_POWERS:
            raise RuntimeError("ideal-power chain did not stabilize")


def is_shift_closed(I: np.ndarray, ell: int, m: int) -> bool:
    """Tripwire: the GF(2) row space must be an R-ideal (shift-closed)."""
    if I.size == 0:
        return True
    r = rank_np(I)
    for a, b in ((1, 0), (0, 1)):
        shifted = np.array([np.roll(np.roll(v.reshape(ell, m), a, axis=0), b, axis=1).reshape(-1)
                            for v in I], np.uint8)
        if rank_np(np.vstack([I, shifted])) != r:
            return False
    return True


def classify_ideal(I: np.ndarray, ell: int, m: int) -> dict:
    """Trichotomy from the ideal alone (Theorem J-G)."""
    if not is_shift_closed(I, ell, m):
        raise RuntimeError("I is not shift-closed: convention error")
    dims = ideal_power_chain(I, ell, m)
    dim_I, dim_inf = dims[0], dims[-1]
    if dim_I == 0:
        case = "degenerate"          # M = 0: no quotient classes at all
    elif dim_inf == 0:
        case = "demote_full"
    elif dim_inf == dim_I:
        case = "immune"
    else:
        case = "mixed"
    return {
        "dim_I": int(dim_I),
        "power_dims": [int(d) for d in dims],
        "dim_I_infty": int(dim_inf),
        "k_pred": 2 * int(dim_I),
        "dim_S_pred": 2 * int(dim_inf),
        "case": case,
        "nilpotency_index": len(dims) - 1 if dim_inf == 0 else None,
    }


def ideal_of_parent(HX: np.ndarray, ell: int, m: int) -> np.ndarray:
    """I = L_pre = left null space of [A B], a GF(2) basis of the ideal."""
    L = nullspace_np(HX.T)
    if L.size == 0:
        return np.zeros((0, ell * m), np.uint8)
    basis, _ = rref_np(L)
    return basis[: rank_np(basis)]


# ------------------------------------------------------- exact demote semantics

def demote_census(HX: np.ndarray, HZ: np.ndarray, ell: int, m: int, k_P: int) -> dict:
    """Brute-force EXP-047 semantics: for every nonzero class y in M decide
    demote(y) <=> y not in I*y.  Returns exact counts (k_P must be small)."""
    Kb = nullspace_np(HX)
    Q, R, kq = E47.quotient_and_projector(HZ, Kb)
    if kq != k_P:
        raise RuntimeError(f"quotient dim {kq} != k_P {k_P}")
    Rk = R[:, :k_P]
    Lpre = nullspace_np(HX.T)
    demote = 0
    total = (1 << k_P) - 1
    for bits in range(1, 1 << k_P):
        y = np.array([(bits >> i) & 1 for i in range(k_P)], np.uint8)
        z = (y @ Q) % 2
        proj = (Lpre @ E47.shift_matrix_of(z, ell, m) % 2) @ Rk % 2
        r0 = rank_np(proj)
        if rank_np(np.vstack([proj, y[None, :]])) > r0:
            demote += 1
    return {"classes_total": total, "classes_demote": demote,
            "fraction": demote / total}


def module_chain(HX: np.ndarray, HZ: np.ndarray, ell: int, m: int, k_P: int) -> list[int]:
    """EXP-052's module F-chain F_{r+1} = I F_r, F_0 = M (independent route)."""
    Kb = nullspace_np(HX)
    Q, R, kq = E47.quotient_and_projector(HZ, Kb)
    if kq != k_P:
        raise RuntimeError(f"quotient dim {kq} != k_P {k_P}")
    Rk = R[:, :k_P]
    Lpre = nullspace_np(HX.T)
    basis = np.eye(k_P, dtype=np.uint8)
    chain = [k_P]
    while True:
        rows = []
        for b in basis:
            z = (b @ Q) % 2
            rows.append((Lpre @ E47.shift_matrix_of(z, ell, m) % 2) @ Rk % 2)
        nxt, _ = rref_np(np.vstack(rows))
        r_next, r_cur = rank_np(nxt), rank_np(basis)
        if r_next > r_cur:
            raise RuntimeError(f"module chain not monotone: {r_cur} -> {r_next}")
        chain.append(int(r_next))
        if r_next == r_cur:
            return chain
        basis = nxt[:r_next]


def stable_ideal_power(I: np.ndarray, ell: int, m: int) -> np.ndarray:
    """I^infty as a GF(2) row basis (stable power of the ideal)."""
    if I.size == 0 or rank_np(I) == 0:
        return np.zeros((0, ell * m), np.uint8)
    cur, _ = rref_np(I)
    cur = cur[: rank_np(cur)]
    for _ in range(MAX_POWERS):
        rows = [poly_mult(g, b, ell, m) for g in I for b in cur]
        if not rows:
            return np.zeros((0, ell * m), np.uint8)
        nxt, _ = rref_np(np.array(rows, np.uint8))
        nxt = nxt[: rank_np(nxt)]
        if nxt.shape[0] in (0, cur.shape[0]):
            return nxt
        cur = nxt
    raise RuntimeError("ideal power did not stabilize")


def module_annihilator(HX: np.ndarray, HZ: np.ndarray, ell: int, m: int, k_P: int) -> np.ndarray:
    """Ann_R(M) = {lambda in R : lambda acts as 0 on M = ker H_X / S_Z}.

    Structurally Ann_R(M) = (A,B), the ideal GENERATED by a and b, because
    M_beta = (R_beta/(A_beta,B_beta))^2 on every local factor and Ann(R/J) = J.
    It is therefore the Frobenius DUAL of Theorem J-G's ideal
    I = L_pre = Ann_R((A,B)), so dim I + dim Ann_R(M) = dim R.

    The two are NOT interchangeable.  I is full exactly on the zero factors
    (A_beta = B_beta = 0) and zero where either is a unit; Ann_R(M) is full
    exactly on the unit factors and zero on the zero factors -- and on a
    deficient factor it is a proper NONZERO ideal (A_beta,B_beta), which is why
    dim Ann_R(M) is not a "unit-factor dimension".  Theorem J-E uses both ideals
    at once (immunity iff 1 in L_pre + Ann_R(M)); Theorem J-G's
    dim S = 2 dim I^infty holds for I = L_pre ONLY -- substituting Ann_R(M)
    lands the stable image on the unit factors and returns a wrong dim S
    (measured per parent in duality_audit's wrong_route_dim_S).
    """
    Kb = nullspace_np(HX)
    Q, R, kq = E47.quotient_and_projector(HZ, Kb)
    if kq != k_P:
        raise RuntimeError(f"quotient dim {kq} != k_P {k_P}")
    Rk = R[:, :k_P]
    blocks = [E47.shift_matrix_of((b @ Q) % 2, ell, m) @ Rk % 2
              for b in np.eye(k_P, dtype=np.uint8)]        # each lm x k_P
    stack = np.hstack(blocks) % 2
    ann = nullspace_np(stack.T)
    if ann.size == 0:
        return np.zeros((0, ell * m), np.uint8)
    basis, _ = rref_np(ann)
    return basis[: rank_np(basis)]


def reciprocal_ideal(HZ: np.ndarray, ell: int, m: int) -> np.ndarray:
    """The ideal (bar a, bar b) <= R as a GF(2) row basis.

    H_Z = [B^T A^T], so the shifts of the reciprocal polynomials bar b and bar a
    are exactly the rows of its two blocks.  This is the ideal that annihilates
    M = ker H_X / S_Z, because H_X acts through the transposed (reciprocal)
    circulants: M_beta = (R_beta/(bar a_beta, bar b_beta))^2.
    """
    dim = ell * m
    basis, _ = rref_np(np.vstack([HZ[:, :dim], HZ[:, dim:]]))
    return basis[: rank_np(basis)]


def ideal_generated_by_AB(HX: np.ndarray, ell: int, m: int) -> np.ndarray:
    """The ideal (a,b) <= R as a GF(2) row basis (all shifts of a and of b)."""
    dim = ell * m
    basis, _ = rref_np(np.vstack([HX[:, :dim], HX[:, dim:]]))
    return basis[: rank_np(basis)]

def duality_audit(HX: np.ndarray, HZ: np.ndarray, ell: int, m: int, k_P: int) -> dict:
    """Three independent facts about the pair of ideals, plus the mislabel probe.

    (i) Ann_R(M) = (bar a, bar b), the RECIPROCAL ideal -- computed two
        completely different ways: the left side from the quotient ACTION on
        M = ker H_X/S_Z, the right side from the H_Z blocks.  The bar is
        essential: H_X acts through transposed circulants, so the naive ideal
        (a,b) has the right DIMENSION but is a different subspace whenever
        (a,b) is not bar-invariant (measured: 46 of the 202 catalogue parents).
    (ii) Frobenius duality: dim I + dim Ann_R(M) = dim R = ell*m, since
         I = Ann_R((a,b)) and F_2[G] is a Frobenius algebra (bar preserves
         dimensions).
    (iii) Mislabel probe: substituting Ann_R(M) for I in dim S = 2 dim I^infty
          gives a different (wrong) number -- recorded per parent so the
          substitution can never pass unnoticed.

    What is NOT claimed: dim Ann_R(M) is not a "unit-factor dimension" -- on a
    deficient factor Ann(M_beta) is the proper nonzero ideal
    (bar a_beta, bar b_beta).  The zero/deficient/unit factor partition is a
    statement about which factors occur; its machine proxies are dim I^infty
    (zero part) and dim I - dim I^infty (nilpotent part on deficient factors).
    """
    I = ideal_of_parent(HX, ell, m)
    ann_M = module_annihilator(HX, HZ, ell, m, k_P)
    recip = reciprocal_ideal(HZ, ell, m)
    ab = ideal_generated_by_AB(HX, ell, m)
    dim_I, dim_ann = int(rank_np(I)), int(rank_np(ann_M))
    dim_recip, dim_ab = int(rank_np(recip)), int(rank_np(ab))

    def _same(P, Q, dP, dQ):
        return dP == dQ and (dP == 0 or int(rank_np(np.vstack([P, Q]))) == dP)

    zero_part = int(rank_np(stable_ideal_power(I, ell, m)))
    return {
        "dim_R": ell * m,
        "dim_I": dim_I,
        "dim_ann_M": dim_ann,
        "dim_reciprocal_ideal": dim_recip,
        "dim_ideal_AB": dim_ab,
        "ann_M_equals_reciprocal_ideal": bool(_same(ann_M, recip, dim_ann, dim_recip)),
        "ideal_AB_is_bar_invariant": bool(_same(ab, recip, dim_ab, dim_recip)),
        "frobenius_identity": dim_I + dim_ann == ell * m,
        "zero_part_dim": zero_part,
        "nilpotent_part_dim": dim_I - zero_part,
        "dim_S_correct_route": 2 * zero_part,
        "wrong_route_dim_S": 2 * int(rank_np(stable_ideal_power(ann_M, ell, m))),
    }


def duality_check() -> dict:
    """Run the duality audit over every catalogue parent and the mixed witness.

    Verified per parent: Ann_R(M) = (bar a, bar b) computed two independent ways, the
    Frobenius identity dim I + dim Ann_R(M) = dim R, the case signature in terms
    of the zero/nilpotent parts of I, and that swapping in Ann_R(M) changes the
    answer (so the mislabel cannot pass silently)."""
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    recs, bad, interchangeable = [], [], []
    for fp, entry in sorted(parents.items(), key=lambda kv: (kv[1]["n"], kv[0])):
        ell, m = int(entry["ell"]), int(entry["m"])
        k_P = int(entry["k_parent"])
        _, HX, HZ = E27.parent_matrices(rows[entry["members"][0]["catalogue_index"]])
        aud = duality_audit(HX, HZ, ell, m, k_P)
        case = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)["case"]
        label = entry["members"][0]["label"]
        ok = (aud["ann_M_equals_reciprocal_ideal"] and aud["frobenius_identity"]
              and 2 * aud["dim_I"] == k_P
              and (case != "immune" or aud["nilpotent_part_dim"] == 0)
              and (case != "demote_full" or aud["zero_part_dim"] == 0))
        if not ok:
            bad.append({"label": label, "case": case, **aud})
        if aud["wrong_route_dim_S"] == aud["dim_S_correct_route"]:
            interchangeable.append({"label": label, "case": case, **aud})
        recs.append({"label": label, "case": case, **aud})
    ell, m = MIXED_WITNESS["ell"], MIXED_WITNESS["m"]
    HX, HZ = bb_from_terms(ell, m, MIXED_WITNESS["A_terms"], MIXED_WITNESS["B_terms"])
    k_P = HX.shape[1] - rank_np(HX) - rank_np(HZ)
    witness = duality_audit(HX, HZ, ell, m, k_P)
    return {
        "parents": len(recs), "records": recs, "bad": bad,
        "annihilators_interchangeable": interchangeable,
        "witness": witness,
        "witness_has_zero_and_nilpotent_parts": witness["zero_part_dim"] > 0
        and witness["nilpotent_part_dim"] > 0,
        "audit_holds": not bad,
    }


# ------------------------------------------------------------------ experiments

def bb_from_terms(ell: int, m: int, a_terms, b_terms):
    spec = E27.BBSpec(ell=ell, m=m, A=[tuple(t) for t in a_terms],
                      B=[tuple(t) for t in b_terms])
    HX, HZ = E27.build_bb(spec)
    return HX, HZ


def catalogue_check() -> dict:
    """Predict dim S for all 202 catalogue parents from the ideal alone and
    compare against EXP-052's module chain."""
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    ref = json.loads(EXP52.read_text())["records"]
    out, mismatches = [], []
    for fp, entry in sorted(parents.items(), key=lambda kv: (kv[1]["n"], kv[0])):
        ell, m = int(entry["ell"]), int(entry["m"])
        _, HX, _ = E27.parent_matrices(rows[entry["members"][0]["catalogue_index"]])
        cls = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)
        k_P = int(entry["k_parent"])
        r = ref[fp]
        rec = {"fingerprint": fp, "label": entry["members"][0]["label"], "ell": ell, "m": m,
               "n": int(entry["n"]), "k_parent": k_P, **cls,
               "dim_S_exp052": int(r["dim_S"]), "case_exp052": r["case"]}
        if cls["k_pred"] != k_P or cls["dim_S_pred"] != int(r["dim_S"]) or cls["case"] != r["case"]:
            mismatches.append(rec)
        out.append(rec)
    return {"parents": len(out), "records": out, "mismatches": mismatches,
            "cases": {c: sum(1 for r in out if r["case"] == c)
                      for c in ("demote_full", "immune", "mixed", "degenerate")}}


MIXED_WITNESS = {
    "ell": 2, "m": 3,
    # A = (1+x)(y+y^2);  B = A*y  (both weight 4)
    "A_terms": [[0, 1], [0, 2], [1, 1], [1, 2]],
    "B_terms": [[0, 2], [0, 0], [1, 2], [1, 0]],
}


def mixed_witness() -> dict:
    """Explicit mixed parent, verified by three independent routes."""
    ell, m = MIXED_WITNESS["ell"], MIXED_WITNESS["m"]
    HX, HZ = bb_from_terms(ell, m, MIXED_WITNESS["A_terms"], MIXED_WITNESS["B_terms"])
    n = HX.shape[1]
    k_P = n - rank_np(HX) - rank_np(HZ)
    cls = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)
    chain = module_chain(HX, HZ, ell, m, k_P)
    census = demote_census(HX, HZ, ell, m, k_P)
    s = cls["dim_S_pred"]
    predicted_demote = (1 << k_P) - (1 << s)
    return {
        **MIXED_WITNESS, "n": int(n), "k_parent": int(k_P), **cls,
        "module_chain": chain, "dim_S_module": chain[-1],
        "census": census, "predicted_classes_demote": int(predicted_demote),
        "agree_ideal_vs_module": chain[-1] == s,
        "agree_ideal_vs_census": census["classes_demote"] == predicted_demote,
        "fraction_pred": 1.0 - ((1 << s) - 1) / ((1 << k_P) - 1),
    }


ODD_LATTICES = [(3, 3), (5, 3), (7, 3), (5, 5), (9, 3), (9, 5), (15, 3), (7, 7), (15, 5)]


def _weight3_terms(ell: int, m: int):
    """Translation-normalised weight-3 supports containing the monomial 1."""
    pts = [(a, b) for a in range(ell) for b in range(m)][1:]
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            yield [(0, 0), pts[i], pts[j]]


def odd_lattice_scan(limit_per_lattice: int = 400) -> dict:
    """Corollary check: on odd x odd lattices every parent with k_P > 0 must be
    demote-immune (I idempotent), because R is semisimple."""
    recs, violations = [], []
    for ell, m in ODD_LATTICES:
        tested = 0
        cases: dict[str, int] = {}
        terms = list(_weight3_terms(ell, m))
        step = max(1, len(terms) // max(1, int(limit_per_lattice ** 0.5)))
        picks = terms[::step]
        for a in picks:
            for b in picks:
                if tested >= limit_per_lattice:
                    break
                HX, _ = bb_from_terms(ell, m, a, b)
                cls = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)
                cases[cls["case"]] = cases.get(cls["case"], 0) + 1
                if cls["case"] not in ("immune", "degenerate"):
                    violations.append({"ell": ell, "m": m, "A": a, "B": b, **cls})
                tested += 1
            if tested >= limit_per_lattice:
                break
        recs.append({"ell": ell, "m": m, "tested": tested, "cases": cases})
    return {"lattices": recs, "violations": violations,
            "corollary_holds": not violations}


def semisimple_baseline() -> dict:
    """The published BB baseline instances, classified by the ideal route.
    [[90,8,10]] has (ell,m) = (15,3): odd x odd, so the corollary predicts it is
    demote-immune -- i.e. no PBB sibling of it can lose X-distance through the
    single-row syzygy channel."""
    specs = {
        "[[72,12,6]]": (6, 6, [[3, 0], [0, 1], [0, 2]], [[0, 3], [1, 0], [2, 0]]),
        "[[90,8,10]]": (15, 3, [[9, 0], [0, 1], [0, 2]], [[0, 0], [2, 0], [7, 0]]),
        "[[108,8,10]]": (9, 6, [[3, 0], [0, 1], [0, 2]], [[0, 3], [1, 0], [2, 0]]),
        "[[144,12,12]]": (12, 6, [[3, 0], [0, 1], [0, 2]], [[0, 3], [1, 0], [2, 0]]),
        "[[288,12,18]]": (12, 12, [[3, 0], [0, 2], [0, 7]], [[0, 3], [1, 0], [2, 0]]),
        "[[360,12,<=24]]": (30, 6, [[9, 0], [0, 1], [0, 2]], [[0, 3], [25, 0], [26, 0]]),
        "[[756,16,<=34]]": (21, 18, [[3, 0], [0, 10], [0, 17]], [[0, 5], [3, 0], [19, 0]]),
    }
    out = []
    for label, (ell, m, a, b) in specs.items():
        HX, HZ = bb_from_terms(ell, m, a, b)
        k_P = HX.shape[1] - rank_np(HX) - rank_np(HZ)
        cls = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)
        out.append({"label": label, "ell": ell, "m": m, "n": int(HX.shape[1]),
                    "k_parent": int(k_P), "odd_lattice": ell % 2 == 1 and m % 2 == 1, **cls})
    return {"baselines": out}


def two_part(v: int, n: int) -> int:
    """Coordinate of v in the 2-primary factor of Z_n (CRT: v mod 2^s)."""
    s = 0
    while n % 2 == 0:
        n //= 2
        s += 1
    return v % (1 << s) if s else 0


def coset_support(terms, ell: int, m: int) -> set:
    """The set of G_2-cosets occupied by the support of a polynomial."""
    return {(two_part(int(a), ell), two_part(int(b), m)) for a, b in terms}


def coset_criterion() -> dict:
    """Theorem J-I corollary, tested on the catalogue: for a trinomial pair the
    parent is demote-immune iff both A and B occupy a SINGLE G_2-coset (and then
    a common character zero exists).  Perfect separation is the prediction."""
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    table: dict[str, int] = {}
    exceptions = []
    weights: dict[str, int] = {}
    for fp, entry in parents.items():
        row = rows[entry["members"][0]["catalogue_index"]]
        ell, m = int(row["ell"]), int(row["m"])
        A = [tuple(t) for t in E27.terms(row, "A_terms")]
        B = [tuple(t) for t in E27.terms(row, "B_terms")]
        weights[f"{len(set(A))},{len(set(B))}"] = weights.get(f"{len(set(A))},{len(set(B))}", 0) + 1
        single = len(coset_support(A, ell, m)) == 1 and len(coset_support(B, ell, m)) == 1
        _, HX, _ = E27.parent_matrices(row)
        case = classify_ideal(ideal_of_parent(HX, ell, m), ell, m)["case"]
        key = f"{case}|single_coset={single}"
        table[key] = table.get(key, 0) + 1
        if (case == "immune") != single:
            exceptions.append({"label": entry["members"][0]["label"], "case": case,
                               "single_coset": single, "ell": ell, "m": m, "A": A, "B": B})
    return {"table": table, "exceptions": exceptions, "weight_shapes": weights,
            "criterion_exact": not exceptions}


LEMMA_LATTICES = [(6, 6), (9, 6), (12, 6), (4, 3), (2, 3), (8, 6), (4, 9), (6, 5)]


def lemma_checks(trials: int = 260, seed: int = 53) -> dict:
    """Falsification attempts against the two structural lemmas.

    L1 (exact-vanishing coset law): for weight <= 3, if A vanishes exactly on
       some local factor (Ann(A)^infty != 0) then supp(A) lies in ONE G_2-coset.
    L2 (scalar law): if supp(A) lies in one G_2-coset then A is zero-or-unit at
       every local factor, i.e. Ann(A) is an idempotent ideal.
    L3 (tightness): the weight-4 witness A = (1+x)(y+y^2) on (2,3) vanishes
       exactly somewhere while occupying TWO cosets -- so 4 is the true
       threshold and L1 fails at weight 4.
    """
    import random
    rng = random.Random(seed)
    l1_fail, l2_fail, tested = [], [], 0
    for ell, m in LEMMA_LATTICES:
        pts = [(a, b) for a in range(ell) for b in range(m)]
        for _ in range(trials // len(LEMMA_LATTICES)):
            w = rng.choice([1, 2, 3])
            terms = rng.sample(pts, w)
            HX, _ = bb_from_terms(ell, m, terms, terms)   # B = A: probes A alone
            ann = ideal_of_parent(HX, ell, m)
            dims = ideal_power_chain(ann, ell, m) if ann.size else [0, 0]
            single = len(coset_support(terms, ell, m)) == 1
            if dims[-1] > 0 and not single:
                l1_fail.append({"ell": ell, "m": m, "terms": terms, "dims": dims})
            if single and dims[-1] != dims[0]:
                l2_fail.append({"ell": ell, "m": m, "terms": terms, "dims": dims})
            tested += 1
    wt4 = MIXED_WITNESS["A_terms"]
    HX, _ = bb_from_terms(MIXED_WITNESS["ell"], MIXED_WITNESS["m"], wt4, wt4)
    dims4 = ideal_power_chain(ideal_of_parent(HX, MIXED_WITNESS["ell"], MIXED_WITNESS["m"]),
                              MIXED_WITNESS["ell"], MIXED_WITNESS["m"])
    cos4 = coset_support(wt4, MIXED_WITNESS["ell"], MIXED_WITNESS["m"])
    return {
        "tested": tested, "L1_failures": l1_fail, "L2_failures": l2_fail,
        "L1_holds": not l1_fail, "L2_holds": not l2_fail,
        "L3_weight4_witness": {"terms": wt4, "cosets": sorted(cos4),
                               "power_dims": dims4,
                               "vanishes_exactly": dims4[-1] > 0,
                               "multi_coset": len(cos4) > 1},
    }


def run(args: argparse.Namespace) -> int:
    t0 = time.time()
    payload = {
        "schema": SCHEMA,
        "utc": E52.utc_now(),
        "catalogue": catalogue_check(),
        "mixed_witness": mixed_witness(),
        "odd_lattice_scan": odd_lattice_scan(args.odd_budget),
        "baselines": semisimple_baseline(),
        "coset_criterion": coset_criterion(),
        "lemma_checks": lemma_checks(),
        "duality_audit": duality_check(),
    }
    cat = payload["catalogue"]
    payload["verdict"] = {
        "ideal_route_matches_exp052": not cat["mismatches"],
        "parents": cat["parents"],
        "cases": cat["cases"],
        "nilpotency_index_2_on_all_demoting": all(
            r["nilpotency_index"] == 2 for r in cat["records"] if r["case"] == "demote_full"),
        "mixed_exists": payload["mixed_witness"]["case"] == "mixed"
        and payload["mixed_witness"]["agree_ideal_vs_module"]
        and payload["mixed_witness"]["agree_ideal_vs_census"],
        "odd_corollary_holds": payload["odd_lattice_scan"]["corollary_holds"],
        "coset_criterion_exact": payload["coset_criterion"]["criterion_exact"],
        "lemmas_hold": payload["lemma_checks"]["L1_holds"] and payload["lemma_checks"]["L2_holds"],
        "duality_audit_holds": payload["duality_audit"]["audit_holds"],
        "annihilators_never_interchangeable":
            payload["duality_audit"]["annihilators_interchangeable"] == [],
        "wall_s": round(time.time() - t0, 2),
    }
    E52.atomic_write_json(OUT, payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(required=True)
    r = sub.add_parser("run")
    r.add_argument("--odd-budget", type=int, default=400)
    r.set_defaults(fn=run)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
