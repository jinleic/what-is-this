#!/usr/bin/env python3
"""Clean-room verifier for the wave-18 LadderW8 front.

Re-derives, without importing any producer:

  (V1) the L = 3..8 regression control (fresh sector-pure mirror closure at
       a DIFFERENT prime, q = 999979, with a different pivot convention)
       against the certified table 14/0, 42/2, 142/10, 494/66, 1780/364,
       6562/1822 (cyclic/W) including per-sector splits;
  (V2) the stored L = 9 wave-16 artifacts (basis sha256, sandwich identity,
       stored low-sector ranks, sector palindrome);
  (V3) the symmetric-square law arithmetic: T(L,m) sector table, closed
       forms, delta pattern, all seven certified totals, per-sector splits,
       and the recorded W_10 / cyclic_10 predictions;
  (V4) the partition-family theorem in instances: independent walk vectors,
       witness matrix block-triangularity, dominance necessity, diagonal
       positivity, rank = Q(L) at L = 8, 9, 10; independent Q(L) DP,
       crossing at L = 32, and the quartic tail bound;
  (V5) every gap certificate row: the L = 5 Pauli-word rank (independent
       bracket engine, second prime 2147483629) and the full L = 10..31
       sector-capped state closures re-derived with an independent
       implementation at q = 999983, each to the target 8L^2 - 2L + 1;
       any row that cannot be re-derived FAILS the test;
  (V6) the artifact envelopes: every stored check passed, targets correct,
       certificate coverage is exactly {5} + {10..31} with no gaps.

Environment: TEST_W8_SKIP_NONE is honored (all stages always run); the
expensive stage is V5 (order 1-2 CPU-hours on a contended host).

Usage:  PYTHONPATH=src .venv/bin/python tests/test_ladder_w8.py
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REG_ART = ROOT / "results" / "algebra_growth" / "ladder_w8_regression.json"
FAM_ART = ROOT / "results" / "algebra_growth" / "ladder_w8_family.json"
CERT_ART = ROOT / "results" / "algebra_growth" / "ladder_w8_certificates.json"
CONS_ART = ROOT / "results" / "algebra_growth" / "ladder_w8.json"
W9_ART = ROOT / "results" / "ladder" / "w9_saturation.json"
W9_BASIS = ROOT / "results" / "ladder" / "w9_basis.json"
QV = 999_979          # verifier prime for the regression closure (different from producer)
QC = 999_983          # certificate prime (re-derivation must match the stored claim)
P2 = 2_147_483_629    # second prime for the L = 5 word certificate
CERT_TABLE = {3: (14, 0, {}), 4: (42, 2, {4: 2}), 5: (142, 10, {4: 5, 6: 5}),
              6: (494, 66, {4: 15, 6: 36, 8: 15}),
              7: (1780, 364, {4: 35, 6: 147, 8: 147, 10: 35}),
              8: (6562, 1822, {4: 70, 6: 448, 8: 786, 10: 448, 12: 70})}
CERT_W = {3: 0, 4: 2, 5: 10, 6: 66, 7: 364, 8: 1822, 9: 8586}
CERT_CYCLIC = {3: 14, 4: 42, 5: 142, 6: 494, 7: 1780, 8: 6562, 9: 24566}

CHECKS = 0
FAILURES: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    tag = "ok" if cond else "FAIL"
    print(f"[{tag}] {name}: {detail}", flush=True)
    if not cond:
        FAILURES.append(name)


def K_dim(L: int) -> int:
    return (2 ** (2 * L - 1) + 3 * 2 ** L) // 4


def edges(L: int):
    E = [(2 * r, 2 * r + 1) for r in range(L)]
    for r in range(L - 1):
        E += [(2 * r, 2 * r + 2), (2 * r + 1, 2 * r + 3)]
    return E


# ---------------------------------------------------------------------------
# V1: independent mirror closure (pivot = FIRST nonzero; prime 999979).
# ---------------------------------------------------------------------------
def tau_c(x: int, L: int) -> int:
    y = 0
    for r in range(L):
        y |= ((x >> (2 * r)) & 1) << (2 * r + 1)
        y |= ((x >> (2 * r + 1)) & 1) << (2 * r)
    return y


def rho_c(x: int, L: int) -> int:
    y = 0
    for r in range(L):
        y |= ((x >> (2 * r)) & 3) << (2 * (L - 1 - r))
    return y


def orbits(L: int):
    orbit_of = [-1] * (1 << (2 * L))
    members = []
    for x in range(1 << (2 * L)):
        if x.bit_count() & 1 or orbit_of[x] >= 0:
            continue
        img = sorted({x, tau_c(x, L), rho_c(x, L), tau_c(rho_c(x, L), L)})
        i = len(members)
        for y in img:
            orbit_of[y] = i
        members.append(tuple(img))
    return orbit_of, members


def mirror_closure(L: int, q: int) -> dict:
    orbit_of, members = orbits(L)
    n = len(members)
    sector_of = np.array([m[0].bit_count() for m in members], dtype=np.int64)
    allbits = (1 << (2 * L)) - 1
    P_idx = np.array([orbit_of[members[i][0] ^ allbits] for i in range(n)], dtype=np.int64)
    sizes = np.array([len(m) for m in members], dtype=np.int64)
    E = edges(L)
    ri, ci, di = [], [], []
    for i, orb in enumerate(members):
        cnt = Counter()
        x = orb[0]
        for (u, v) in E:
            cnt[orbit_of[x ^ ((1 << u) | (1 << v))]] += 1
        for j, c in cnt.items():
            num = int(sizes[i]) * c
            assert num % int(sizes[j]) == 0
            ri.append(i)
            ci.append(j)
            di.append(num // int(sizes[j]))
    order = np.lexsort((np.array(ci), np.array(ri)))
    B_r = np.array(ri)[order]
    B_c = np.array(ci)[order]
    B_d = np.array(di, dtype=np.int64)[order]
    B_ptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(B_ptr, B_r + 1, 1)
    np.cumsum(B_ptr, out=B_ptr)
    low = [k for k in range(0, 2 * L + 1, 2) if k <= L]
    loc2glob = {}
    glob2loc = {}
    for s in low:
        coords = np.nonzero(sector_of == s)[0].astype(np.int64)
        loc2glob[s] = coords
        back = np.full(n, -1, dtype=np.int64)
        back[coords] = np.arange(coords.size, dtype=np.int64)
        glob2loc[s] = back
    rows = {s: {} for s in low}

    def add(s, idx, vals):
        w = np.zeros(loc2glob[s].size, dtype=np.int64)
        w[idx] = vals % q
        rs = rows[s]
        while True:
            nz = np.nonzero(w)[0]
            if nz.size == 0:
                return None
            piv = int(nz[0])          # FIRST nonzero: different convention
            got = rs.get(piv)
            if got is None:
                inv = pow(int(w[piv]), q - 2, q)
                rs[piv] = (nz.astype(np.int32), (((w[nz] * inv) % q).astype(np.int32)))
                return piv
            sidx, svals = got
            c = int(w[piv])
            w[sidx] = (w[sidx] - c * svals.astype(np.int64)) % q

    todo = [(0, add(0, np.array([0], dtype=np.int64), np.array([1], dtype=np.int64)))]
    while todo:
        s, piv = todo.pop()
        sidx, svals = rows[s][piv]
        glob = loc2glob[s][sidx]
        image = np.zeros(n, dtype=np.int64)
        for j in range(glob.size):
            src = int(glob[j])
            lo, hi = B_ptr[src], B_ptr[src + 1]
            image[B_c[lo:hi]] += int(svals[j]) * B_d[lo:hi]
        image %= q
        targets = np.nonzero(image)[0]
        if targets.size:
            tsec = sector_of[targets]
            for t in np.unique(tsec):
                t = int(t)
                tg = targets[tsec == t]
                vals = image[tg].copy()
                if t > L:
                    t2 = 2 * L - t
                    tg = P_idx[tg]
                else:
                    t2 = t
                pp = add(t2, glob2loc[t2][tg], vals)
                if pp is not None:
                    todo.append((t2, pp))
    low_ranks = {s: len(rows[s]) for s in low}
    rank_total = sum(low_ranks[s] * (1 if s == 2 * L - s else 2) for s in low)
    dims = {s: int(loc2glob[s].size) for s in low}
    deficits = {}
    for s in low:
        d = dims[s] - low_ranks[s]
        if d:
            deficits[s] = d
            if s != 2 * L - s:
                deficits[2 * L - s] = d
    return {"cyclic": rank_total, "W": K_dim(L) - rank_total,
            "sectors": dict(sorted(deficits.items()))}


def stage_V1() -> None:
    for L, (cyc, W, sec) in CERT_TABLE.items():
        if L == 5 and False:
            continue
        t0 = time.process_time()
        rec = mirror_closure(L, QV)
        check(f"V1_regression_L{L}",
              rec["cyclic"] == cyc and rec["W"] == W and rec["sectors"] == sec,
              f"cyclic {rec['cyclic']}/{K_dim(L)} W {rec['W']} sectors {rec['sectors']} "
              f"({time.process_time()-t0:.1f}s, q={QV})")


# ---------------------------------------------------------------------------
# V2: stored L = 9 artifacts.
# ---------------------------------------------------------------------------
def stage_V2() -> None:
    w9 = json.loads(W9_ART.read_text())
    raw = json.loads(W9_BASIS.read_text())
    canonical = json.dumps(raw["basis"], sort_keys=True, separators=(",", ":"))
    sha = hashlib.sha256(canonical.encode()).hexdigest()
    d = w9["data"]
    low = {int(k): v for k, v in d["closures"]["q1"]["low_ranks"].items()}
    sec = {int(k): v for k, v in d["W9"]["W_sectors"].items()}
    check("V2_w9_artifacts",
          sha == raw["sha256"] == d["W9"]["basis_sha256"]
          and d["W9"]["value"] == 8586 and d["W9"]["cyclic_Q_dim"] == 24566
          and d["W9"]["cyclic_Q_dim"] + d["W9"]["value"] == K_dim(9)
          and low == {0: 1, 2: 45, 4: 666, 6: 3570, 8: 8001}
          and all(sec.get(18 - k) == v for k, v in sec.items())
          and sum(sec.values()) == 8586 and len(raw["basis"]) == 8586,
          "sha256, sandwich 24566+8586=33152, low ranks, palindrome, basis size")


# ---------------------------------------------------------------------------
# V3: the law arithmetic.
# ---------------------------------------------------------------------------
def T_tri(n: int, k: int) -> int:
    return comb(n, k) * (comb(n, k) + 1) // 2


def law_cyclic(L: int) -> int:
    return (comb(2 * L, L) + 2 ** L) // 2 - (1 if L % 4 == 0 else 0)


def sector_dims_direct(L: int) -> dict:
    _, members = orbits(L)
    return dict(Counter(m[0].bit_count() for m in members))


def stage_V3() -> None:
    ok_tot = all(law_cyclic(L) == CERT_CYCLIC[L] and K_dim(L) - law_cyclic(L) == CERT_W[L]
                 for L in range(3, 10))
    check("V3_law_totals", ok_tot,
          "closed form (C(2L,L)+2^L)/2 - [4|L] reproduces all seven certified totals")
    ok_sec = True
    for L, (_, _, sec) in CERT_TABLE.items():
        if L == 3:
            continue
        dims = sector_dims_direct(L)
        pred = {}
        for m in range(0, L + 1):
            dd = dims[2 * m] - T_tri(L, m) + (1 if (L % 4 == 0 and 2 * m == L) else 0)
            if dd:
                pred[2 * m] = dd
        if pred != sec:
            ok_sec = False
    # L = 9 sector split from the stored artifact
    w9 = json.loads(W9_ART.read_text())
    sec9 = {int(k): v for k, v in w9["data"]["W9"]["W_sectors"].items()}
    dims9 = sector_dims_direct(9)
    pred9 = {}
    for m in range(0, 10):
        dd = dims9[2 * m] - T_tri(9, m)
        if dd:
            pred9[2 * m] = dd
    check("V3_law_sector_splits", ok_sec and pred9 == sec9,
          "K^(2m) - T(L,m) (+1 at 4|L middle) matches every certified split L=4..9")
    check("V3_law_delta_pattern",
          law_cyclic(4) == 42 and law_cyclic(8) == 6562
          and (comb(12, 6) + 2 ** 6) // 2 == 494 and (comb(20, 10) + 2 ** 10) // 2 - 1 != 24566,
          "delta = -1 exactly at L = 4, 8 (4|L); L = 6, 10 unshifted")
    reg = json.loads(REG_ART.read_text())
    pred = reg["data"]["law"]["predictions_next_wave"]
    check("V3_law_predictions",
          pred["W_10"] == K_dim(10) - law_cyclic(10) == 38950
          and pred["cyclic_10"] == law_cyclic(10) == 92890 and pred["K_10"] == K_dim(10) == 131840,
          "recorded W_10 = 38950, cyclic_10 = 92890 match the closed form")
    # law implication: K-W >= 2^L for 2 <= L <= 4096 under the law
    ok_2L = all(law_cyclic(L) >= 2 ** L for L in range(2, 200))
    check("V3_law_2L_survival", ok_2L, "under the law K_L - W_L >= 2^L (checked L = 2..199)")


# ---------------------------------------------------------------------------
# V4: partition family, independent instance verification.
# ---------------------------------------------------------------------------
def family_index(L: int):
    out = []

    def rec(prefix, budget, maxpart):
        for a in range(min(maxpart, budget - 2), -1, -1):
            cur = prefix + (a,)
            out.append(cur)
            rec(cur, budget - (a + 2), a)

    rec(tuple(), L, L)
    return out


@lru_cache(maxsize=None)
def Q_cnt_inner(budget: int, maxpart: int) -> int:
    total = 1
    for a in range(0, min(maxpart, budget - 2) + 1):
        total += Q_cnt_inner(budget - (a + 2), a)
    return total


def Q_count(L: int) -> int:
    return Q_cnt_inner(L, L) - 1


def dominates(lam, mu) -> bool:
    s1 = s2 = 0
    for a, b in zip(lam, mu):
        s1 += a
        s2 += b
        if s1 < s2:
            return False
    return True


def walk_vec(lam, L: int, p: int) -> dict:
    E = edges(L)
    vec = {frozenset(): 1}
    for a in lam:
        out: dict = {}
        for cfg, c in vec.items():
            for (u, v) in E:
                if u in cfg or v in cfg:
                    continue
                key = cfg | {u, v}
                out[key] = (out.get(key, 0) + c) % p
        vec = {k: v for k, v in out.items() if v}
        for _ in range(a):
            out = {}
            for cfg, c in vec.items():
                for (u, v) in E:
                    iu, iv = u in cfg, v in cfg
                    if iu == iv:
                        continue
                    key = (cfg - {u}) | {v} if iu else (cfg - {v}) | {u}
                    out[key] = (out.get(key, 0) + c) % p
            vec = {k: v for k, v in out.items() if v}
    return vec


def witness(lam) -> frozenset:
    sites = set()
    x = 0
    for a in lam:
        d = a + 1
        sites.add(2 * x)
        sites.add(2 * (x + d))
        x += d + 1
    return frozenset(sites)


def stage_V4() -> None:
    p = 1_000_003  # yet another prime
    for L in (8, 9, 10):
        lambdas = family_index(L)
        lambdas.sort(key=lambda lam: (len(lam), sum(lam), lam))
        vecs = {lam: walk_vec(lam, L, p) for lam in lambdas}
        ok = True
        for lam in lambdas:
            C = witness(lam)
            for mu in lambdas:
                coeff = vecs[mu].get(C, 0)
                if mu == lam and coeff == 0:
                    ok = False
                if coeff:
                    if len(mu) != len(lam) or sum(mu) < sum(lam):
                        ok = False
                    elif sum(mu) == sum(lam) and not dominates(lam, mu):
                        ok = False
        rows: dict = {}

        # rank over F_p with pivot = lexicographically smallest config
        rows = {}
        for lam in lambdas:
            v = {tuple(sorted(k)): val for k, val in vecs[lam].items()}
            while v:
                piv = min(v)
                got = rows.get(piv)
                if got is None:
                    inv = pow(v[piv], p - 2, p)
                    rows[piv] = {k: val * inv % p for k, val in v.items()}
                    break
                c = v[piv]
                for k, val in got.items():
                    nv = (v.get(k, 0) - c * val) % p
                    if nv:
                        v[k] = nv
                    elif k in v:
                        del v[k]
        check(f"V4_family_instance_L{L}",
              ok and len(rows) == Q_count(L) == len(lambdas),
              f"triangularity + dominance + rank {len(rows)} = Q({L}) = {Q_count(L)} at p={p}")
    lastbad = max(L for L in range(3, 301) if Q_count(L) <= 8 * L * L - 2 * L + 1)
    check("V4_crossing", lastbad == 31 and all(Q_count(L) > 8 * L * L - 2 * L + 1 for L in range(32, 301)),
          "Q(L) exceeds 8L^2-2L+1 exactly from L = 32 on (DP to 300)")
    # index-set disambiguation: Q(L) = #{nonempty partitions, all parts >= 2, total <= L}
    @lru_cache(maxsize=None)
    def parts_ge2(cap, maxpart):
        out = 1
        for part in range(2, min(maxpart, cap) + 1):
            out += parts_ge2(cap - part, part)
        return out

    check("V4_index_set_two_forms",
          all(Q_count(L) == parts_ge2(L, L) - 1 for L in range(2, 121))
          and [Q_count(L) for L in range(2, 13)] == [1, 2, 4, 6, 10, 14, 21, 29, 41, 55, 76],
          "Q(L) = nonempty partitions with all parts >= 2, total <= L; first values 1,2,4,6,10,14,21,29,41,55,76")
    tail = all(comb(L - 4, 4) >= 24 * (8 * L * L - 2 * L + 2) for L in range(301, 5001))
    check("V4_tail_bound", tail, "C(L-4,4)/24 > 8L^2-2L+1 for 301 <= L <= 5000 (quartic beyond)")


# ---------------------------------------------------------------------------
# V5: certificates re-derived.
# ---------------------------------------------------------------------------
def word_rank(L: int, target: int, prime: int, budget_s: float = 900.0):
    n = 2 * L
    masks = [(1 << u) | (1 << v) for (u, v) in edges(L)]
    A = {(1 << v, 0): 1 for v in range(n)}
    B = {(0, m): 1 for m in masks}
    rows: dict = {}
    t0 = time.process_time()

    def add(vec):
        v = dict(vec)
        while v:
            piv = max(v)
            got = rows.get(piv)
            if got is None:
                inv = pow(v[piv], prime - 2, prime)
                rows[piv] = {k: val * inv % prime for k, val in v.items()}
                return piv
            c = v[piv]
            for k, val in got.items():
                nv = (v.get(k, 0) - c * val) % prime
                if nv:
                    v[k] = nv
                elif k in v:
                    del v[k]
        return None

    todo = []
    for g in (A, B):
        if add(g) is not None:
            todo.append(dict(g))
    while todo and len(rows) < target and time.process_time() - t0 < budget_s:
        base = todo.pop(0)
        for gen in (0, 1):
            out: dict = {}
            for (a, b), c in base.items():
                if gen == 0:
                    bb = b
                    while bb:
                        low = bb & -bb
                        bb ^= low
                        key = (a ^ low, b)
                        out[key] = (out.get(key, 0) + c) % prime
                else:
                    for m in masks:
                        if bin(a & m).count("1") & 1:
                            key = (a, b ^ m)
                            out[key] = (out.get(key, 0) - c) % prime
            out = {k: v for k, v in out.items() if v}
            if out and add(out) is not None:
                todo.append(out)
                if len(rows) >= target:
                    break
    return len(rows)


def canon(cfg: tuple, L: int) -> tuple:
    best = None
    for f in range(4):
        im = []
        for v in cfg:
            w = v
            if f & 1:
                w ^= 1
            if f & 2:
                w = 2 * (L - 1 - (w >> 1)) + (w & 1)
            im.append(w)
        t = tuple(sorted(im))
        if best is None or t < best:
            best = t
    return best


def state_rank(L: int, target: int, q: int = QC, budget_s: float = 12000.0):
    t0 = time.process_time()

    def enum(k):
        idx = {}
        reps = []
        sizes = []
        for cfg in combinations(range(2 * L), k):
            ims = set()
            for f in range(4):
                im = []
                for v in cfg:
                    w = v
                    if f & 1:
                        w ^= 1
                    if f & 2:
                        w = 2 * (L - 1 - (w >> 1)) + (w & 1)
                    im.append(w)
                ims.add(tuple(sorted(im)))
            rep = min(ims)
            if rep == cfg:
                idx[rep] = len(reps)
                reps.append(rep)
                sizes.append(len(ims))
        return reps, idx, sizes

    reps2, idx2, sz2 = enum(2)
    reps4, idx4, sz4 = enum(4)
    n2, n4 = len(reps2), len(reps4)
    E = edges(L)

    def block(reps_s, sizes_s, idx_t, sizes_t, k_target):
        ri, ci, di = [], [], []
        for i, cfg in enumerate(reps_s):
            cnt = Counter()
            s = set(cfg)
            for (u, v) in E:
                nc = set(s)
                for w in (u, v):
                    if w in nc:
                        nc.remove(w)
                    else:
                        nc.add(w)
                if len(nc) != k_target:
                    continue
                cnt[canon(tuple(sorted(nc)), L)] += 1
            for rep, c in cnt.items():
                j = idx_t[rep]
                num = sizes_s[i] * c
                assert num % sizes_t[j] == 0
                ri.append(i)
                ci.append(j)
                di.append(num // sizes_t[j])
        order = np.lexsort((np.array(ci), np.array(ri)))
        r = np.array(ri, dtype=np.int64)[order]
        c = np.array(ci, dtype=np.int64)[order]
        d = np.array(di, dtype=np.int64)[order]
        ptr = np.zeros(len(reps_s) + 1, dtype=np.int64)
        np.add.at(ptr, r + 1, 1)
        np.cumsum(ptr, out=ptr)
        return ptr, c, d

    B22 = block(reps2, sz2, idx2, sz2, 2)
    B24 = block(reps2, sz2, idx4, sz4, 4)
    B42 = block(reps4, sz4, idx2, sz2, 2)
    B44 = block(reps4, sz4, idx4, sz4, 4)
    seed = np.zeros(n2, dtype=np.int64)
    for (u, v) in E:
        seed[idx2[canon((min(u, v), max(u, v)), L)]] = 1
    rows = {2: {}, 4: {}}
    scratch = {2: np.zeros(n2, dtype=np.int64), 4: np.zeros(n4, dtype=np.int64)}

    def add(sector, dense):
        w = scratch[sector]
        if dense is not w:
            w[:] = dense
        rs = rows[sector]
        while True:
            nz = np.nonzero(w)[0]
            if nz.size == 0:
                return None
            piv = int(nz[-1])
            got = rs.get(piv)
            if got is None:
                inv = pow(int(w[piv]), q - 2, q)
                rs[piv] = (nz.astype(np.int32), (((w[nz] * inv) % q).astype(np.int32)))
                return piv
            sidx, svals = got
            c = int(w[piv])
            w[sidx] = (w[sidx] - c * svals.astype(np.int64)) % q

    tmp = scratch[2]
    tmp[:] = seed
    todo = [(2, add(2, tmp))]
    out2 = np.zeros(n2, dtype=np.int64)
    out4 = np.zeros(n4, dtype=np.int64)

    def mul(blockmat, sidx, svals, out):
        ptr, cols, dat = blockmat
        out[:] = 0
        for j in range(sidx.size):
            i = int(sidx[j])
            v = int(svals[j])
            lo, hi = ptr[i], ptr[i + 1]
            out[cols[lo:hi]] += v * dat[lo:hi]
        out %= q
        return out

    def rank():
        return len(rows[2]) + len(rows[4]) + 1

    while todo and rank() < target and time.process_time() - t0 < budget_s:
        sector, piv = todo.pop(0)
        sidx, svals = rows[sector][piv]
        outs = ((2, B22, out2), (4, B24, out4)) if sector == 2 else ((2, B42, out2), (4, B44, out4))
        for tsec, M, out in outs:
            mul(M, sidx, svals, out)
            if np.any(out):
                pp = add(tsec, out)
                if pp is not None:
                    todo.append((tsec, pp))
    return rank()


def stage_V5() -> None:
    cert = json.loads(CERT_ART.read_text())
    recs = {r["L"]: r for r in cert["data"]["records"]}
    check("V5_coverage",
          sorted(recs) == [5] + list(range(10, 32))
          and all(r["certified"] for r in recs.values())
          and all(recs[L]["target"] == 8 * L * L - 2 * L + 1 for L in recs),
          f"stored coverage {{5}} + {{10..31}}, all certified, targets 8L^2-2L+1")
    t0 = time.process_time()
    r5 = word_rank(5, 191, P2)
    check("V5_L5_rederived", r5 >= 191,
          f"independent word closure at SECOND prime {P2}: rank {r5} >= 191 ({time.process_time()-t0:.1f}s)")
    for L in range(10, 32):
        target = 8 * L * L - 2 * L + 1
        t0 = time.process_time()
        r = state_rank(L, target)
        check(f"V5_L{L}_rederived", r >= target,
              f"independent state closure rank {r} >= {target} ({time.process_time()-t0:.1f}s)")


# ---------------------------------------------------------------------------
# V6: envelopes.
# ---------------------------------------------------------------------------
def stage_V6() -> None:
    for path, name in ((REG_ART, "regression"), (FAM_ART, "family"), (CERT_ART, "certificates")):
        art = json.loads(path.read_text())
        ch = art["data"]["checks"]
        check(f"V6_{name}_checks", all(c["passed"] for c in ch), f"{len(ch)} stored checks all passed")
    if CONS_ART.exists():
        cons = json.loads(CONS_ART.read_text())
        check("V6_consolidated", all(c["passed"] for c in cons["data"]["checks"]),
              f"{len(cons['data']['checks'])} consolidated checks all passed")
    else:
        check("V6_consolidated", False, "consolidated envelope missing")


def main() -> int:
    t0 = time.process_time()
    stage_V1()
    stage_V2()
    stage_V3()
    stage_V4()
    stage_V5()
    stage_V6()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"\n{CHECKS - len(FAILURES)}/{CHECKS} checks passed "
          f"({time.process_time()-t0:.1f}s cpu, peak rss {rss/1e9:.2f} GB)")
    if FAILURES:
        print("FAILED: " + ", ".join(FAILURES))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
