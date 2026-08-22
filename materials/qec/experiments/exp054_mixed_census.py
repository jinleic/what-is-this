"""EXP-054: exhaustive weight census of the demote trichotomy over BB lattices.

Theorem J-G (EXP-053) reduces the demote law to the ideal I = Ann_R(A,B):

    dim S = 2 dim I^infty,   I^infty = V(A) cap V(B),   V(A) := Ann_R(A)^infty,

so the trichotomy of a pair (A,B) is decided by two per-polynomial invariants:

    Ann(A) = {lambda : lambda A = 0}        (a GF(2) subspace of R)
    V(A)   = Ann(A)^infty = e_{V(A)} R     (the idempotent ideal of the local
                                            factors on which A vanishes exactly)

  dim I       = dim(Ann(A) cap Ann(B))          -> k_P = 2 dim I
  dim I^infty = dim(V(A) cap V(B))              -> dim S = 2 dim I^infty
  case: degenerate (dim I = 0) / demote_full (dim I^infty = 0)
        / immune (dim I^infty = dim I) / mixed (otherwise)

Both invariants are translation-invariant (Ann(x^a y^b A) = Ann(A)), so a census
over supports normalised to contain the monomial 1 is exhaustive over all pairs
up to independent translations of A and B.  Grouping polynomials by the pair of
canonical bases collapses the O(P^2) pair sweep to O(G^2) with G = #groups.

Questions decided here, per lattice:
  Q1  does any weight-<=3 pair (the catalogue's own shape) give a MIXED parent?
  Q2  what is the exact case census over all weight-<=3 pairs?
  Q3  if mixed occurs, what is the minimal-weight witness (and its (k_P, dim S))?

Cross-validation: a random sample of pairs per lattice is re-decided by EXP-053's
direct route (build [A B], take L_pre, iterate ideal powers) and must agree.

Artifacts: results/processed/exp054_mixed_census.json
Run: python experiments/exp054_mixed_census.py run [--lattices 9x6,12x12] [--max-weight 3]
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import random
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


E53 = _load("exp053", "exp053_ideal_classification.py")

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

SCHEMA = "exp054-mixed-census-v1"
OUT = ROOT / "results" / "processed" / "exp054_mixed_census.json"

# catalogue lattices + the two "underexplored" ones + a small parity grid
DEFAULT_LATTICES = [
    (3, 6), (6, 3), (6, 6), (9, 6), (12, 6), (15, 6), (30, 6), (12, 12), (15, 12),
    (2, 3), (4, 3), (2, 5), (4, 5), (6, 5), (2, 9), (4, 9), (8, 3), (6, 9),
]


def canon(basis: np.ndarray) -> bytes:
    if basis.size == 0:
        return b""
    rr, _ = rref_np(basis)
    r = rank_np(rr)
    return rr[:r].tobytes()


def ideal_generators(I: np.ndarray, ell: int, m: int) -> np.ndarray:
    """Small R-module generating set of the ideal spanned by rows of I."""
    if I.size == 0:
        return I
    target = rank_np(I)
    gens: list[np.ndarray] = []
    cur = np.zeros((0, ell * m), np.uint8)
    for v in I:
        if rank_np(np.vstack([cur, v[None, :]])) == rank_np(cur):
            continue
        gens.append(v)
        shifts = np.array([np.roll(np.roll(v.reshape(ell, m), a, axis=0), b, axis=1).reshape(-1)
                           for a in range(ell) for b in range(m)], np.uint8)
        cur, _ = rref_np(np.vstack([cur, shifts]))
        cur = cur[: rank_np(cur)]
        if rank_np(cur) == target:
            break
    return np.array(gens, np.uint8)


def stable_power(I: np.ndarray, ell: int, m: int) -> np.ndarray:
    """I^infty as a GF(2) row basis (the idempotent ideal e_full R)."""
    if I.size == 0 or rank_np(I) == 0:
        return np.zeros((0, ell * m), np.uint8)
    gens = ideal_generators(I, ell, m)
    cur, _ = rref_np(I)
    cur = cur[: rank_np(cur)]
    while True:
        rows = [E53.poly_mult(g, b, ell, m) for g in gens for b in cur]
        if not rows:
            return np.zeros((0, ell * m), np.uint8)
        nxt, _ = rref_np(np.array(rows, np.uint8))
        nxt = nxt[: rank_np(nxt)]
        if nxt.shape[0] == 0 or nxt.shape[0] == cur.shape[0]:
            return nxt
        cur = nxt


def poly_invariants(terms, ell: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """(Ann(A), V(A)) as GF(2) row bases."""
    M = poly_matrix(ell, m, [tuple(t) for t in terms])
    ann = nullspace_np(M.T)
    if ann.size:
        ann, _ = rref_np(ann)
        ann = ann[: rank_np(ann)]
    else:
        ann = np.zeros((0, ell * m), np.uint8)
    return ann, stable_power(ann, ell, m)


def dim_intersection(X: np.ndarray, Y: np.ndarray) -> int:
    dx, dy = rank_np(X) if X.size else 0, rank_np(Y) if Y.size else 0
    if dx == 0 or dy == 0:
        return 0
    return dx + dy - rank_np(np.vstack([X, Y]))


def case_of(dim_I: int, dim_inf: int) -> str:
    if dim_I == 0:
        return "degenerate"
    if dim_inf == 0:
        return "demote_full"
    if dim_inf == dim_I:
        return "immune"
    return "mixed"


def normalised_supports(ell: int, m: int, max_weight: int):
    """Supports containing (0,0), of weight 1..max_weight (translation-normalised)."""
    pts = [(a, b) for a in range(ell) for b in range(m) if (a, b) != (0, 0)]
    yield [(0, 0)]
    for w in range(1, max_weight):
        for extra in itertools.combinations(pts, w):
            yield [(0, 0), *extra]


def census_lattice(ell: int, m: int, max_weight: int, sample: int, rng: random.Random) -> dict:
    t0 = time.time()
    supports = list(normalised_supports(ell, m, max_weight))
    groups: dict[tuple[bytes, bytes], dict] = {}
    for terms in supports:
        ann, V = poly_invariants(terms, ell, m)
        key = (canon(ann), canon(V))
        g = groups.get(key)
        if g is None:
            groups[key] = {"ann": ann, "V": V, "count": 1, "rep": terms,
                           "min_weight": len(terms)}
        else:
            g["count"] += 1
            if len(terms) < g["min_weight"]:
                g["min_weight"], g["rep"] = len(terms), terms
    keys = list(groups)
    counts = {"degenerate": 0, "demote_full": 0, "immune": 0, "mixed": 0}
    best_mixed = None
    for i, ka in enumerate(keys):
        ga = groups[ka]
        for kb in keys[i:]:
            gb = groups[kb]
            dim_I = dim_intersection(ga["ann"], gb["ann"])
            dim_inf = dim_intersection(ga["V"], gb["V"])
            c = case_of(dim_I, dim_inf)
            mult = ga["count"] * gb["count"] * (1 if ka == kb else 2)
            counts[c] += mult
            if c == "mixed":
                w = ga["min_weight"] + gb["min_weight"]
                cand = {"A": ga["rep"], "B": gb["rep"], "weight_total": w,
                        "k_parent": 2 * dim_I, "dim_S": 2 * dim_inf,
                        "wt_A": ga["min_weight"], "wt_B": gb["min_weight"]}
                if best_mixed is None or w < best_mixed["weight_total"]:
                    best_mixed = cand
    # independent cross-check: re-decide sampled pairs through EXP-053's route
    checked, mismatches = 0, []
    for _ in range(sample):
        a, b = rng.choice(supports), rng.choice(supports)
        ann_a, Va = poly_invariants(a, ell, m)
        ann_b, Vb = poly_invariants(b, ell, m)
        fast = case_of(dim_intersection(ann_a, ann_b), dim_intersection(Va, Vb))
        HX, _ = E53.bb_from_terms(ell, m, a, b)
        direct = E53.classify_ideal(E53.ideal_of_parent(HX, ell, m), ell, m)
        if direct["case"] != fast:
            mismatches.append({"A": a, "B": b, "fast": fast, "direct": direct["case"]})
        checked += 1
    return {
        "ell": ell, "m": m, "n": 2 * ell * m, "max_weight": max_weight,
        "supports": len(supports), "groups": len(keys),
        "pairs_total": sum(counts.values()), "counts": counts,
        "mixed_exists": counts["mixed"] > 0, "minimal_mixed": best_mixed,
        "crosscheck_pairs": checked, "crosscheck_mismatches": mismatches,
        "wall_s": round(time.time() - t0, 2),
    }


def run(args: argparse.Namespace) -> int:
    rng = random.Random(args.seed)
    lats = DEFAULT_LATTICES
    if args.lattices:
        lats = [tuple(int(v) for v in tok.split("x")) for tok in args.lattices.split(",")]
    recs = []
    for ell, m in lats:
        rec = census_lattice(ell, m, args.max_weight, args.sample, rng)
        recs.append(rec)
        print(f"({ell},{m}) supports={rec['supports']} groups={rec['groups']} "
              f"pairs={rec['pairs_total']} mixed={rec['counts']['mixed']} "
              f"minimal={rec['minimal_mixed']} xchk={rec['crosscheck_pairs']}/"
              f"{len(rec['crosscheck_mismatches'])} bad {rec['wall_s']}s", flush=True)
    payload = {
        "schema": SCHEMA, "utc": E53.E52.utc_now(), "max_weight": args.max_weight,
        "seed": args.seed, "lattices": recs,
        "verdict": {
            "lattices_scanned": len(recs),
            "lattices_with_mixed": [f"{r['ell']}x{r['m']}" for r in recs if r["mixed_exists"]],
            "lattices_without_mixed": [f"{r['ell']}x{r['m']}" for r in recs if not r["mixed_exists"]],
            "odd_lattices_with_mixed": [f"{r['ell']}x{r['m']}" for r in recs
                                        if r["mixed_exists"] and r["ell"] % 2 and r["m"] % 2],
            "crosscheck_mismatches": sum(len(r["crosscheck_mismatches"]) for r in recs),
            "pairs_total": sum(r["pairs_total"] for r in recs),
        },
    }
    E53.E52.atomic_write_json(OUT, payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(required=True)
    r = sub.add_parser("run")
    r.add_argument("--lattices", default="")
    r.add_argument("--max-weight", type=int, default=3)
    r.add_argument("--sample", type=int, default=25)
    r.add_argument("--seed", type=int, default=54)
    r.set_defaults(fn=run)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
