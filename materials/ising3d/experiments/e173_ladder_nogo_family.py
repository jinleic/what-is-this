#!/usr/bin/env python3
"""Partition-family lower bound for the two-generator ladder algebra: exact
arithmetic, finite-instance verification of the triangularity lemmas, and the
all-L >= 32 crossing of the free-fermion threshold.

THE THEOREM (proved in proofs/ladder_w8.md, verified in instances here).
For the open 2xL ladder let

    A_L = sum_v X_v,  B_L = sum_{(u,v) in E} Z_u Z_v,  psi = |+>^(2L),
    U_L = alg(A_L, B_L) psi   (vacuum cyclic module),
    g_L = <A_L, B_L>_Lie.

Hard-core pieces D (pair creation), F (hop), D+ (pair annihilation) of B_L
lie in the associative algebra alg(A_L, B_L) (Theorem 1 of
proofs/ladder_alll_proof.md and the sector-projector Lagrange argument of
proofs/ladder_l9.md sec. 2).  For a partition
lambda = (a_1 >= ... >= a_m >= 0), m >= 1, define the walk word

    w(lambda) = F^{a_m} D F^{a_m-1} ... D F^{a_1} D psi  in U_L.

Witness configuration C(lambda): m particle pairs on the TOP leg, pair j
with rung distance d_j = a_j + 1, packed left-to-right with particle gap
exactly one rung:

    x_1 = 0,  pair j occupies rungs {x_j, x_j + d_j},  x_{j+1} = x_j + d_j + 1.

C(lambda) exists iff span(lambda) := |lambda| + 2m <= L.

Lemmas (proved in the note, machine-verified here in instances):
  (P)  All matrix entries of F and D in the configuration basis are >= 0,
       so every coefficient of w(lambda') is a count of histories.
  (C)  cost lemma: a history reaching a config whose 2m particles sit on one
       leg at rungs p_1 < ... < p_2m pays at least min over pairings of
       sum(dist) - #(leg-created dimers); rung-created or wrong-leg-created
       dimers pay strictly more; the minimum over pairings is the
       consecutive pairing (p_1p_2)(p_3p_4)..., every other pairing pays
       at least +2.
  (S)  scheduling lemma: hops of the j-th created dimer occur in phases
       >= j; with the exact budget sum(a_i') this forces (Hall/greedy) the
       suffix-domination condition, i.e. C(lambda) is reachable from
       w(lambda') with |lambda'| = |lambda| only if lambda dominates
       lambda'.
  (T)  Therefore, ordering words by (sector = 2m, total = |lambda|, any
       linear extension of dominance), the witness-coefficient matrix is
       block-triangular with strictly positive diagonal, and

           dim_Q g_L >= dim_Q U_L >= Q(L),
           Q(L) := #{lambda : |lambda| + 2 len(lambda) <= L}.

  (A)  Arithmetic: Q(L) > 8L^2 - 2L + 1 for every 32 <= L <= 300 by exact
       DP, and for L > 300 by the explicit bound
       Q(L) >= C(L-8+4, 4)/24 >= (L-8)^4/576 > 8L^2 (the first inequality
       counts partitions into at most 4 parts of total <= L-8, each a valid
       index with m = 4).

This producer verifies (P), (C), (S), (T) exhaustively in finite instances
(L = 8, 9, 10 over F_p) and (A) exactly, and emits the no-go table.

Usage:
  PYTHONPATH=src .venv/bin/python experiments/e173_ladder_nogo_family.py
"""
from __future__ import annotations

import json
import resource
import sys
import time
from functools import lru_cache
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e173_ladder_nogo_family.py"
OUT_PATH = ROOT / "results" / "algebra_growth" / "ladder_w8_family.json"
P1 = 2_147_483_647
VERIFY_L = (8, 9, 10)
DP_CHECK_MAX = 300


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def log(msg: str) -> None:
    print(f"[e173 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def edges(L: int) -> list[tuple[int, int]]:
    E = [(2 * r, 2 * r + 1) for r in range(L)]
    for r in range(L - 1):
        E += [(2 * r, 2 * r + 2), (2 * r + 1, 2 * r + 3)]
    return E


# ---------------------------------------------------------------------------
# Index set and counting.
# ---------------------------------------------------------------------------
def family_index(L: int) -> list[tuple[int, ...]]:
    """All partitions lambda = (a_1 >= ... >= a_m >= 0), m >= 1, |lambda| + 2m <= L."""
    out = []

    def rec(prefix, budget, maxpart):
        for a in range(min(maxpart, budget - 2), -1, -1):
            cur = prefix + (a,)
            out.append(cur)
            rec(cur, budget - (a + 2), a)

    rec(tuple(), L, L)
    return out


def Q_count(L: int) -> int:
    @lru_cache(maxsize=None)
    def cnt(budget, maxpart):
        total = 1
        for a in range(0, min(maxpart, budget - 2) + 1):
            total += cnt(budget - (a + 2), a)
        return total

    return cnt(L, L) - 1


def dominates(lam: tuple, mu: tuple) -> bool:
    """lam >= mu in dominance order (equal length, equal total assumed)."""
    s1 = s2 = 0
    for a, b in zip(lam, mu):
        s1 += a
        s2 += b
        if s1 < s2:
            return False
    return True


# ---------------------------------------------------------------------------
# Finite-instance verification of (P), (C), (S), (T) over F_p.
# ---------------------------------------------------------------------------
def apply_D(vec: dict, E, p: int) -> dict:
    out: dict = {}
    for cfg, c in vec.items():
        s = set(cfg)
        for (u, v) in E:
            if u in s or v in s:
                continue
            key = tuple(sorted(s | {u, v}))
            out[key] = (out.get(key, 0) + c) % p
    return {k: v for k, v in out.items() if v}


def apply_F(vec: dict, E, p: int) -> dict:
    out: dict = {}
    for cfg, c in vec.items():
        s = set(cfg)
        for (u, v) in E:
            iu, iv = u in s, v in s
            if iu == iv:
                continue
            ns = set(s)
            if iu:
                ns.remove(u)
                ns.add(v)
            else:
                ns.remove(v)
                ns.add(u)
            key = tuple(sorted(ns))
            out[key] = (out.get(key, 0) + c) % p
    return {k: v for k, v in out.items() if v}


def walk_vector(lam: tuple, L: int, p: int) -> dict:
    """w(lam) psi with the j-th created dimer followed by phase F^{a_j}:
    the FIRST-created dimer receives the LARGEST phase a_1.  As an operator
    word (leftmost acts last) this is  F^{a_m} D ... F^{a_2} D F^{a_1} D."""
    E = edges(L)
    vec = {tuple(): 1}
    for a in lam:  # a_1 (largest) first: D then F^{a_1}, then D F^{a_2}, ...
        vec = apply_D(vec, E, p)
        for _ in range(a):
            vec = apply_F(vec, E, p)
    return vec


def witness_config(lam: tuple) -> tuple:
    sites = []
    x = 0
    for a in lam:
        d = a + 1
        sites.append(2 * x)          # top-leg site of rung x
        sites.append(2 * (x + d))    # top-leg site of rung x + d
        x = x + d + 1
    return tuple(sorted(sites))


def verify_instance(L: int, p: int) -> dict:
    started = time.process_time()
    lambdas = family_index(L)
    lambdas.sort(key=lambda lam: (len(lam), sum(lam), lam))
    vectors = {lam: walk_vector(lam, L, p) for lam in lambdas}
    diag_ok = True
    tri_ok = True
    dom_ok = True
    matrix_entries = 0
    for lam in lambdas:
        C = witness_config(lam)
        for mu in lambdas:
            coeff = vectors[mu].get(C, 0)
            matrix_entries += 1
            if mu == lam and coeff == 0:
                diag_ok = False
            if coeff:
                if len(mu) != len(lam) or sum(mu) < sum(lam):
                    tri_ok = False  # sector mismatch impossible; budget < cost impossible
                elif sum(mu) == sum(lam) and not dominates(lam, mu):
                    dom_ok = False  # dominance necessity (S)+(C): coeff != 0 => lam >= mu
    # rank of the family over F_p (echelon on full config coordinates)
    rows: dict = {}

    def add(vec: dict):
        v = dict(vec)
        while v:
            piv = max(v)
            got = rows.get(piv)
            if got is None:
                inv = pow(v[piv], p - 2, p)
                rows[piv] = {k: val * inv % p for k, val in v.items()}
                return piv
            c = v[piv]
            for k, val in got.items():
                nv = (v.get(k, 0) - c * val) % p
                if nv:
                    v[k] = nv
                elif k in v:
                    del v[k]
        return None

    for lam in lambdas:
        add(vectors[lam])
    rank = len(rows)
    return {
        "L": L,
        "prime": p,
        "family_size_Q": len(lambdas),
        "Q_count_dp": Q_count(L),
        "diagonal_positive": diag_ok,
        "budget_sector_block_triangular": tri_ok,
        "dominance_necessity": dom_ok,
        "family_rank_mod_p": rank,
        "rank_equals_Q": rank == len(lambdas),
        "matrix_entries_checked": matrix_entries,
        "cpu_seconds": round(time.process_time() - started, 2),
    }


# ---------------------------------------------------------------------------
# Arithmetic (A): crossing and the tail bound.
# ---------------------------------------------------------------------------
def tail_bound_ok(L: int) -> bool:
    r"""Q(L) >= C(L-4, 4)/24 > 8L^2 - 2L + 1 via <=4-part partitions, L > 300.

    #{a_1>=..>=a_4>=0, sum <= N} >= C(N+4,4)/24 with N = L-8."""
    N = L - 8
    return comb(N + 4, 4) >= 24 * (8 * L * L - 2 * L + 2)


def main() -> int:
    started = time.process_time()
    checks = []

    instances = []
    for L in VERIFY_L:
        rec = verify_instance(L, P1)
        instances.append(rec)
        log(f"instance L={L}: Q={rec['family_size_Q']} rank={rec['family_rank_mod_p']} "
            f"diag={rec['diagonal_positive']} tri={rec['budget_sector_block_triangular']} "
            f"dom={rec['dominance_necessity']}")
        checks.append({
            "name": f"instance_L{L}_triangularity",
            "passed": rec["diagonal_positive"] and rec["budget_sector_block_triangular"]
            and rec["dominance_necessity"] and rec["rank_equals_Q"]
            and rec["family_size_Q"] == rec["Q_count_dp"],
            "detail": f"witness matrix block-triangular, diagonal > 0, rank {rec['family_rank_mod_p']} = Q({L})",
        })

    crossing_rows = []
    all_dp = True
    for L in range(3, DP_CHECK_MAX + 1):
        q = Q_count(L)
        thr = 8 * L * L - 2 * L + 1
        ok = q > thr
        if L >= 32 and not ok:
            all_dp = False
        if L <= 45 or L % 25 == 0:
            crossing_rows.append({"L": L, "Q": q, "threshold_8LL_2L_1": thr, "exceeds": ok})
    checks.append({
        "name": "dp_crossing_32_to_300",
        "passed": all_dp and Q_count(31) <= 8 * 31 * 31 - 2 * 31 + 1,
        "detail": "Q(L) > 8L^2-2L+1 for every 32 <= L <= 300; Q(31) does not exceed (family alone starts at 32)",
    })

    tail_ok = all(tail_bound_ok(L) for L in range(301, 2001))
    tail_math = all(
        comb(L - 4, 4) >= 24 * (8 * L * L - 2 * L + 2) for L in range(301, 5001)
    )
    checks.append({
        "name": "tail_bound_L_gt_300",
        "passed": tail_ok and tail_math,
        "detail": "C(L-4,4)/24 (a lower bound for Q via <=4-part partitions of total <= L-8) "
                  "exceeds 8L^2-2L+1 for 301 <= L <= 5000; quartic vs quadratic beyond",
    })

    envelope = {
        "meta": {
            "script": SCRIPT,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version.split()[0],
            "cpu_seconds_total": round(time.process_time() - started, 2),
            "peak_rss_bytes": peak_rss_bytes(),
            "theorem": (
                "dim_Q g_L >= dim_Q U_L >= Q(L) = #{partitions lambda, |lambda| + 2 len(lambda) <= L} "
                "for every L >= 2; Q(L) > 8L^2 - 2L + 1 for every L >= 32."
            ),
            "proof_location": "proofs/ladder_w8.md sections 4-5",
        },
        "data": {
            "instances": instances,
            "crossing_table": crossing_rows,
            "Q_values_32_to_60": {str(L): Q_count(L) for L in range(32, 61)},
            "checks": checks,
        },
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT_PATH.with_suffix(".tmp")
    tmp.write_text(json.dumps(envelope, indent=1, sort_keys=True))
    tmp.replace(OUT_PATH)
    log(f"wrote {OUT_PATH.relative_to(ROOT)}")
    ok = all(c["passed"] for c in checks)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
